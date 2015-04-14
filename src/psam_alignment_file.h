#ifndef PSAM_ALIGNMENT_FILE_H
#define PSAM_ALIGNMENT_FILE_H

#include <stdio.h>
#include <string>
#include <vector>

#include "sam_alignment.h"
//#include "sam_alignment_file.h"

struct ThreadData{
  int ind;
  int *helper1;
  long begin;
  long end;

  void* inData;
  void* outData;
};

class ThreadArray{
 public:
  int threadNum;
  pthread_t* threads;
  void* inData; //data to be analyzed
  void* outData; //data structure to be filled in
  std::vector<long> begins; //block begin for each thread
  std::vector<long> ends; //block end for each thread
  
  //std::vector<ThreadData> threadData;
  ThreadData* threadDataPtr;

  int helper1;

  ThreadArray(int num);
  ~ThreadArray();
  void init_boundaries(long units);
  void create( void * (*start_routine)(void*) );
  int join();

};

struct page_reader_params{
  FILE* input;
  char* page;
  long bSize;
  std::vector<long>* ptrNewLines; //pointer to new lines vector
  int find_new_lines_threadNum;
};

class PBufferedFile{
  //parallel implementation of the file allowing to 
  //split lines in parallel and extract lines

 private:
  FILE* input;
  char* buf; //alias for the current page

  char* page1; //buffers to which the data is read from the file
  char* page2;
  
  pthread_t page1_reader;
  pthread_t page2_reader;

  page_reader_params page1_reader_params;
  page_reader_params page2_reader_params;

  int cur_page; //keeps track of a current page

  long bSize; //size of the buffer
  long bufOffset;
  long filesize;
  long filepos; //file position corresponding to bufOffset

  size_t fread_result;

  int new_line_ind; //keeps track of position in the new_lines
  std::vector<long>* ptrNewLines; //pointer to positions of new line characters in a current page
  std::vector<long> new_lines1;
  std::vector<long> new_lines2; 

 public:
  int threadNum; //number of threads used for parallelization

  PBufferedFile();
  
  void getline(std::string& str);
  bool open(const char* _fname, long _bSize);
  bool close();
  bool good();
  void SetThreadNum(int _threadNume);
  void SetBufferSize(long _bSize);

 private:
  void find_new_lines();
  void read_page1();
  void read_page2();
  void wait_read_page1();
  void wait_read_page2();
}; 

class PSamAlignmentFile{
 private:
  PBufferedFile pBufferedFile;
  std::string fname;
  std::string sam_error_msg;
  int threadNum; //number of threads to use
  long bSize; //buffer size for pBufferedFile
  long samBufferSize; //buffer size for sam entries
  int samAlignmentIndex; // internal index of a current alignment in the buffer
  bool eof; //EOF marker

  std::vector<std::string> stringBuffer; //strings buffered from the sam file
  std::vector<SamAlignment> samAlignmentBuffer;

  void BufferAlignments(); //buffer the alignments

  SamAlignment dummyNullSamAlignment; //dummy to return in case of error
 public:
  PSamAlignmentFile(std::string _fname);
  PSamAlignmentFile(std::string _fname, int _samBufferSize);
  ~PSamAlignmentFile();

  void SetSamBufferSize(int _samBufferSize);
  void SetFileBufferSize(long _bSize);

  bool good();

  void SetThreadNum(int _threadNum);
  SamAlignment& read_next_alignment(); //the reference is only guaranteed until the next alignment is read
  void close();

};

#endif
