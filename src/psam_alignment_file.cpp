#include <iostream>

#include "psam_alignment_file.h"
#include "assert.h"
#include "string_utils.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;

bool pSamAlignementFileDebug=false;

// =====================================================

ThreadArray::ThreadArray(int num){
  threadNum = num;
  threads = new pthread_t[threadNum];
  threadDataPtr = new ThreadData[threadNum];
  assert(threads != NULL);
  assert(threadDataPtr != NULL);
}

ThreadArray::~ThreadArray(){
  delete[] threads;
  delete[] threadDataPtr;
  threads = NULL;
  threadDataPtr = NULL;
}

void ThreadArray::init_boundaries(long units){
  begins.clear();
  ends.clear();

  long units_per_thread = units / threadNum;
  long cur_begin = 0;
  for(int i=0; i<threadNum; i++){
    if(i!=threadNum-1){
      begins.push_back(cur_begin);
      ends.push_back(cur_begin + units_per_thread);
      cur_begin += units_per_thread;
    } else {
      begins.push_back(cur_begin);
      ends.push_back(units);
    }
  }
}

int ThreadArray::join(){
  int res = 0;
  for(int i=0; i<threadNum; i++){
    int cur_join_res = pthread_join((*(threads+i)), NULL);
    if(cur_join_res != 0) res = cur_join_res;
  }
  return res;
}

void* ParseSamStrings(void* ptr){
  ThreadData* threadData = (ThreadData*) ptr;
  int threadInd = threadData->ind;
  int& counter = *(threadData->helper1);

  vector<string>* dummyPtr1 = (vector<string>*) threadData->inData;
  vector<string>& samStrings = (*dummyPtr1);

  vector<SamAlignment>* dummyPtr2 = (vector<SamAlignment>*) threadData->outData;
  vector<SamAlignment>& samAlignments = (*dummyPtr2);

  //cout<<endl;
  //cout<<endl<<"curThread: "<<threadInd<<" [ "<<threadData->begin<<","<<threadData->end<<" ]"<<endl;
  //cout<<"sam strings: "<<samStrings.size()<<endl;
  
  if(samStrings.size()==0) return 0;
  
  assert(threadData->begin < (int) samStrings.size());
  assert(threadData->end <= (int) samStrings.size());
  assert(samStrings.size() == samAlignments.size());

  if(pSamAlignementFileDebug){
    if(threadInd == 0) cout<<endl; }
  for(int i=threadData->begin; i<threadData->end; i++){ 
    if(pSamAlignementFileDebug){
      if(counter%10000 == 0 && threadInd == 0){ 
	cout<<"\rParsed: "<<commify(counter)<<" sam strings            "; cout.flush();
      }
    }
    string sam_error_msg;
    SamAlignment tmpSamAlignment(samStrings[i], sam_error_msg);
    samAlignments[i] = tmpSamAlignment;
    //samAlignments[i].read(samStrings[i].c_str());
    counter++;
  }
  if(pSamAlignementFileDebug){ if(threadInd == 0) cout<<endl;}
  return 0;
}

void* _find_new_lines(void *ptr){
  ThreadData* threadData = (ThreadData*) (ptr);
  int threadInd = threadData->ind;
  char* buf = (char*) (threadData->inData);

  vector< vector<long> >* dummy_ptr1 = (vector< vector<long> >*) (threadData->outData);  
  vector< vector<long> >& dummy_ref1 = (*dummy_ptr1);

  assert(threadData->ind < (int) dummy_ref1.size());
  vector<long>& new_lines = dummy_ref1[threadInd];

  
  for(long i=threadData->begin; i<threadData->end; i++){
    if(buf[i] == '\n'){
      new_lines.push_back(i);
    }
  }
  return 0;
}

void ThreadArray::create(void *(*start_routine)(void*)){
  for(int i=0; i<threadNum; i++){
    threadDataPtr[i].ind = i;
    threadDataPtr[i].inData = inData;
    threadDataPtr[i].outData = outData;
    threadDataPtr[i].begin = begins[i];
    threadDataPtr[i].end = ends[i];
    threadDataPtr[i].helper1 = &helper1;

    int iret = pthread_create( threads+i, NULL, start_routine, (void*) (threadDataPtr+i));
    assert(iret == 0);
  }
  
}

// =====================================================

PBufferedFile::PBufferedFile(){
  buf = NULL;
  bSize = -1;

  page1 = NULL;
  page2 = NULL;

  filepos = -1;
  new_line_ind = -1;
  threadNum = 8;
}

bool PBufferedFile::open(const char* _fname, long _bSize){
  assert(_bSize>0);
  if(page1 != NULL){ delete[] page1; page1 = NULL;}
  if(page2 != NULL){ delete[] page2; page2 = NULL;}

  bSize = _bSize;
  bufOffset = -1;
  filepos = 0;
  input = fopen(_fname, "r");

  if(input == NULL){
    perror("Error openeing file");
    return false;
  }

  fseek(input, 0, SEEK_END);
  filesize = ftell(input);
  rewind(input);

  page1 = new char[bSize];
  page2 = new char[bSize];

  page1_reader_params.input = input;
  page1_reader_params.bSize = bSize;
  page1_reader_params.page = page1;
  page1_reader_params.find_new_lines_threadNum = threadNum;
  page1_reader_params.ptrNewLines = &(new_lines1);
  
  page2_reader_params.input = input;
  page2_reader_params.bSize = bSize;
  page2_reader_params.page = page2;
  page2_reader_params.find_new_lines_threadNum = threadNum;
  page2_reader_params.ptrNewLines = &(new_lines2);

  cur_page = 1; //1 for page1, 2 for page2  
  
  if(cur_page == 1) buf = page1;
  if(cur_page == 2) buf = page2;
  
  return true;
}

void PBufferedFile::SetThreadNum(int _threadNum){
  assert(_threadNum > 0);
  threadNum = _threadNum;
}

bool PBufferedFile::close(){
  new_lines1.clear();
  new_lines2.clear();
  
  if(page1 != NULL){ delete[] page1; page1 = NULL;}
  if(page2 != NULL){ delete[] page2; page2 = NULL;}

  fclose(input);

  input = NULL;  
  return true;
}

void* read_page(void* ptr){
  page_reader_params* pars = (page_reader_params*)(ptr);  
  //cout<<"Reading page from file..."; cout.flush();
  size_t fread_result = fread(pars->page, 1, pars->bSize, pars->input);
  //cout<<"done!"<<endl;

  //find new lines
  int threadNum = pars->find_new_lines_threadNum;
  vector<long>* ptrNewLines = pars->ptrNewLines;
  ptrNewLines->clear();

  vector< vector<long> > threadResults(threadNum); //this will hold new lines per thread

  ThreadArray threads(threadNum);
  threads.inData = (void*) (pars->page);
  threads.outData = (void*)(&threadResults);

  threads.init_boundaries((long)fread_result);
  threads.create(_find_new_lines);
  threads.join();

  for(unsigned int i=0; i<threadResults.size(); i++){
    ptrNewLines->insert(ptrNewLines->end(), threadResults[i].begin(), threadResults[i].end());
  }

  return 0;
}

void PBufferedFile::read_page1(){
  int iret = pthread_create(&page1_reader, NULL, read_page, (void*)(&page1_reader_params));
  assert(iret == 0);
}
void PBufferedFile::read_page2(){
  int iret = pthread_create(&page2_reader, NULL, read_page, (void*)(&page2_reader_params));
  assert(iret == 0);
}
void PBufferedFile::wait_read_page1(){
  pthread_join(page1_reader, NULL);
}

void PBufferedFile::wait_read_page2(){
  pthread_join(page2_reader, NULL);
}

void PBufferedFile::getline(string& str){
  str = "";
  if(buf == NULL){
    printf("Insfufficient memory to read the file.\n");
    return;
  }

  if(bufOffset == -1){
    //first time to read the file
    if(pSamAlignementFileDebug)
      cout<<endl<<"Buffering first time from the file..."<<endl;
    read_page1();
    wait_read_page1();
    buf = page1;
    cur_page = 1;
    ptrNewLines = &new_lines1;

    read_page2();

    bufOffset = 0;
    filepos = 0;
    new_line_ind = -1;
  }

  bool _continue = true;
  while(_continue){
    if(new_line_ind < (int)ptrNewLines->size()-1){ //more lines are present
      long line_length = (*ptrNewLines)[new_line_ind+1] - bufOffset;
      string tmp_string(buf + bufOffset, line_length);
      if(str.length() == 0) str = tmp_string;
      else str = str + tmp_string;

      bufOffset += line_length+1; filepos += line_length+1;

      //cout<<"read: "<<endl;
      //cout<<str<<"--->"<<endl;
      //cout<<"filepos:  "<<filepos<<endl;
      //cout<<"filesize: "<<filesize<<endl;
      new_line_ind++;
      _continue = false;
    } else {
      //no more new lines. 
      //check if eof
      if(filepos == filesize){ filepos++; bufOffset++; _continue=false; return; }
      
      //copy the remainder of the buffer and read the new buffer
      long remainder_length;
      if(filepos + bSize - bufOffset < filesize){
	remainder_length = bSize - bufOffset;
      } else { //EOF
	remainder_length = filesize - filepos;
	_continue = false;
      }

      //part of the line is in one buffer, part in another
      string tmp_string(buf + bufOffset, remainder_length);
      str = str + tmp_string;
      bufOffset=0; filepos += remainder_length;
      if(filepos < filesize){

	if(cur_page == 1){
	  wait_read_page2();
	  buf = page2;
	  cur_page = 2;
	  ptrNewLines = &new_lines2;
	  read_page1();
	}
	else if(cur_page == 2){
	  wait_read_page1();
	  buf = page1;
	  cur_page = 1;
	  ptrNewLines = &new_lines1;
	  read_page2();
	}
	
	bufOffset = 0;
	new_line_ind = -1;
      }
    }
  }

  return;
}

/*
void PBufferedFile::find_new_lines(){
  new_lines.clear();

  vector< vector<long> > threadResults(threadNum); //this will hold new lines per thread
  
  ThreadArray threads(threadNum);
  threads.inData = (void*) buf;
  threads.outData = (void*)(&threadResults);

  threads.init_boundaries((long)fread_result);
  threads.create(_find_new_lines);
  threads.join();

  for(unsigned int i=0; i<threadResults.size(); i++){
    new_lines.insert(new_lines.end(), threadResults[i].begin(), threadResults[i].end());
  }
  //cout<<endl;
  cout<<"New lines: "<<new_lines.size()<<endl;
}
*/
bool PBufferedFile::good(){
  if(input == NULL) return false;
  if(filepos <= filesize) return true;
  else return false;
}

void PBufferedFile::SetBufferSize(long _bSize){
  assert(buf == NULL); //cannot change the size of the buffer mid-way through reading the file
  bSize = _bSize;
}

PSamAlignmentFile::PSamAlignmentFile(string _fname){
  fname = _fname;
  threadNum = 8;
  samBufferSize = 100000;
  bSize = samBufferSize*244;
  samAlignmentIndex = -1;

  pBufferedFile.open(fname.c_str(), bSize);
  eof = !pBufferedFile.good();
  stringBuffer.reserve(samBufferSize);
  samAlignmentBuffer.reserve(samBufferSize);
}
PSamAlignmentFile::PSamAlignmentFile(string _fname, int _samBufferSize){
  fname = _fname;
  threadNum = 8;
  samBufferSize = _samBufferSize;
  bSize = samBufferSize*244;
  samAlignmentIndex = -1;

  pBufferedFile.open(fname.c_str(), bSize);
  eof = !pBufferedFile.good();
  stringBuffer.reserve(samBufferSize);
  samAlignmentBuffer.reserve(samBufferSize);
}

void PSamAlignmentFile::SetSamBufferSize(int _samBufferSize){
  samBufferSize = _samBufferSize;
}

void PSamAlignmentFile::SetFileBufferSize(long _bSize){
  bSize = _bSize;
  pBufferedFile.SetBufferSize(bSize);
}

void PSamAlignmentFile::SetThreadNum(int _threadNum){
  threadNum = _threadNum;
  pBufferedFile.SetThreadNum(threadNum);
}

SamAlignment& PSamAlignmentFile::read_next_alignment(){
  if(eof) return dummyNullSamAlignment;

  if(samAlignmentIndex >= 0 && samAlignmentIndex < (int) samAlignmentBuffer.size()-1){
    samAlignmentIndex++;
    return samAlignmentBuffer[samAlignmentIndex];
  }
  else if(samAlignmentIndex != -1 && samAlignmentIndex == (int) samAlignmentBuffer.size()-1){
    //End of buffer, check for EOF and buffer more alignments

    if(pBufferedFile.good()){
      BufferAlignments();
      samAlignmentIndex = 0;
      if(samAlignmentBuffer.size() !=0 ) return samAlignmentBuffer[samAlignmentIndex];
      else  return dummyNullSamAlignment;
    } else {
      eof = true; return dummyNullSamAlignment;
    }
  }
  else if(samAlignmentIndex == -1){
    //first time reading the alignment    
    if(pBufferedFile.good()){
      BufferAlignments();
    }
    else eof = true;
    samAlignmentIndex = 0;
    if(samAlignmentBuffer.size() !=0 ){
      return samAlignmentBuffer[samAlignmentIndex];
    }
    else return dummyNullSamAlignment;
  }
  
  cout<<"samAlignmentIndex: "<<samAlignmentIndex<<endl;
  cout<<"samAlignmentBuffer.size(): "<<samAlignmentBuffer.size()<<endl;
  assert(false); //you should not be here
  return dummyNullSamAlignment;
}

void PSamAlignmentFile::BufferAlignments(){
  if(pSamAlignementFileDebug){
    cout<<endl<<"resizing..."; cout.flush();}

  //SamAlignment dummy_al;
  stringBuffer.resize(samBufferSize); 
  samAlignmentBuffer.resize(samBufferSize); //, dummy_al); 
  if(pSamAlignementFileDebug)cout<<" done!"<<endl;

  if(!pBufferedFile.good()) return;

  if(pSamAlignementFileDebug){
    cout<<endl;
    cout<<"Buffering sam strings..."<<endl;
  }
  for(unsigned int i=0; i<samBufferSize && pBufferedFile.good(); i++){
    if(pSamAlignementFileDebug){
      if(i%100000 == 0)
	cout<<"\rBuffered "<<commify(i)<<" sam strings"; cout.flush();
    }
    string curString;
    pBufferedFile.getline(curString);

    if(pBufferedFile.good()){
      //stringBuffer.push_back(curString);
      stringBuffer[i] = curString;
    } else {
      stringBuffer.resize(i);
      samAlignmentBuffer.resize(i);
    }
  }
  if(pSamAlignementFileDebug){ cout<<endl;
    cout<<"Buffering alignments..."<<endl;
  }

  samAlignmentBuffer.resize(stringBuffer.size());

  //parse sam strings
  ThreadArray threads(threadNum);
  threads.inData = (void*) (&stringBuffer);
  threads.outData = (void*)(&samAlignmentBuffer);
  threads.helper1 = 0;

  threads.init_boundaries((long)stringBuffer.size());
  //cout<<endl;
  threads.create(ParseSamStrings);
  int join_res = threads.join();
  //cout<<"join_res: "<<join_res<<endl;
  assert(join_res == 0);
  //cout<<endl;
}

bool PSamAlignmentFile::good(){
  if(eof) return false;
  if(samAlignmentIndex == -1) return true;
  if(samAlignmentIndex >= 0 && samAlignmentIndex < (int) samAlignmentBuffer.size()) return true;
  if(samAlignmentBuffer.size() == 0) return false; //alignments have been read, but the buffer is empty -> EOF
  assert(false); //you should not be here
  return false;
}

void PSamAlignmentFile::close(){
  stringBuffer.clear();
  samAlignmentBuffer.clear();
  pBufferedFile.close();
}

PSamAlignmentFile::~PSamAlignmentFile(){
  //file should be closed explicitly using close, no need for destructor
}
