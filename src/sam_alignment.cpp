#include <string>
#include <iostream>
#include <vector>
#include <fstream>

#include <stdlib.h>
#include "assert.h"

//#include "cigar.h"
//#include "alignment.h"
#include "sam_alignment.h"
#include "string_utils.h"

using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::cout;
using std::ofstream;

bool number(char c){
  if(c=='0' || c=='1' || c=='2' || c=='3' || c=='4' || 
     c=='5' || c=='6' || c=='7' || c=='8' || c=='9'){
    return true;
  } 
  return false;
}

bool nucl(char c){
  if(c=='A' || c=='C' || c=='G' || c=='T' || c=='N'){
    return true;
  }
  return false;
}

/*
bool __number(char c){
  if(c=='1' || c=='2' || c=='3' || c=='4' || c=='5' ||
     c=='6' || c=='7' || c=='8' || c=='9' || c=='0'){
    return true;
  }
  return false;
  }
*/

int StringViewCigar::reference_length(){
  //returns the length of the alignments spanned by this cigar string                                          
  int res = 0;
  for(unsigned int i=0; i<numbers.size(); i++){
    if(chars[i] != 'S' && chars[i] != 'H' && chars[i] != 'I'){
      res += get_cigar_element_length(i);
    }
  }
  return res;
}

int StringViewCigar::get_cigar_element_length(int ind){
  //returns the reference length of the element with index ind                                                 
  switch(chars[ind]){
  case 'M': return numbers[ind];
  case 'I': return 0;
  case 'D': return numbers[ind];
  case 'N': return numbers[ind];
  case 'S': return 0;
  case 'H': return 0;
  case 'P': return numbers[ind];
  case '=': return numbers[ind];
  case 'X': return numbers[ind];
  default: cerr<<"Unexpected CIGAR char"<<endl; exit(1);
  }
}

StringViewCigar::StringViewCigar(const StringView& inp_string){//cigar string parser             
  init(inp_string);
}
int StringViewCigar::init(const StringView& inp_string){
  int curi = 0;

  while(curi < (int)inp_string.length()){
    if(!number(inp_string[curi])){
      cout<<endl;
      cout<<"cigar: "<<inp_string<<endl;
      cout<<"curi: "<<curi<<" char: "<<inp_string[curi]<<endl;
      cout<<"sam_string: "<<endl;
      cout<<(*(inp_string.ptrStr))<<endl;
      return -1;
      //assert(false);
    }
    //assert(number(inp_string[curi]));

    //string cur_num = "";
    StringView cur_num;
    cur_num.ptrStr = inp_string.ptrStr;
    cur_num.begin = inp_string.begin + curi;
    cur_num._length = 0;

    bool continue_reading_num = true;
    while(continue_reading_num){
      if(number(inp_string[curi])){
        //cur_num = cur_num + inp_string[curi];
	cur_num._length++;
        curi++;
        if(curi>=(int)inp_string.length()) continue_reading_num = false;
      } else {
        continue_reading_num = false;
      }
    }
    char cur_op = inp_string[curi];
    if(cur_op != 'M' && cur_op != 'S' && cur_op != 'D' && cur_op != 'I' && cur_op != 'H'){
      cerr<<"Could not understand CIGAR char "<<cur_op<<" in "<<inp_string<<endl;
      return -1;
      //exit(1);
    }
    curi++;

    numbers.push_back(atoi(cur_num.str().c_str()));
    chars.push_back(cur_op);
  }
  return 0;
}

/*
void StringViewCigar::init(const StringView& inp_string){
  int curi = 0;

  while(curi < (int)inp_string.length()){
    assert(number(inp_string[curi]));

    string cur_num = "";
    bool continue_reading_num = true;
    while(continue_reading_num){
      if(number(inp_string[curi])){
        cur_num = cur_num + inp_string[curi];
        curi++;
        if(curi>=(int)inp_string.length()) continue_reading_num = false;
      } else {
        continue_reading_num = false;
      }
    }
    char cur_op = inp_string[curi];
    if(cur_op != 'M' && cur_op != 'S' && cur_op != 'D' && cur_op != 'I' && cur_op != 'H'){
      cerr<<"Could not understand CIGAR char "<<cur_op<<" in "<<inp_string<<endl;
      exit(1);
    }
    curi++;

    numbers.push_back(atoi(cur_num.c_str()));
    chars.push_back(cur_op);
  }
}
*/
/*
cigar::cigar(const string& inp_string){//cigar string parser                                                   
  int curi = 0;

  while(curi < (int)inp_string.length()){
    assert(__number(inp_string[curi]));

    string cur_num = "";
    bool continue_reading_num = true;
    while(continue_reading_num){
      if(__number(inp_string[curi])){
        cur_num = cur_num + inp_string[curi];
        curi++;
        if(curi>=(int)inp_string.length()) continue_reading_num = false;
      } else {
        continue_reading_num = false;
      }
    }
    char cur_op = inp_string[curi];
    if(cur_op != 'M' && cur_op != 'S' && cur_op != 'D' && cur_op != 'I' && cur_op != 'H'){
      cerr<<"Could not understand CIGAR char "<<cur_op<<" in "<<inp_string<<endl;
      exit(1);
    }
    curi++;

    //cout<<inp_string<<" parsed: "<<cur_op<<" "<<cur_num<<endl;                                               
    numbers.push_back(atoi(cur_num.c_str()));
    chars.push_back(cur_op);
  }
}
*/
/*
mdfield::mdfield(const string& inp_string){
  if(inp_string.length() == 0)  return;

  assert(inp_string.length() >= 3);
  assert(inp_string.substr(0,5) == "MD:Z:"); //the string starts with MD:Z:

  //parse mdfield string
  //unsigned int last_i = 5;
  unsigned int cur_i = 5;
  while(cur_i < inp_string.length()){
    bool operation_parsed = false;
    if(number(inp_string[cur_i])){
      operation_parsed = true;
      string cur_num;
      unsigned int ni = cur_i;
      while(ni<inp_string.length() && number(inp_string[ni])){
	cur_num = cur_num + inp_string[ni];
	ni++;
      }
      cur_i = ni;
      numbers.push_back(atoi(cur_num.c_str()));
      ops.push_back('M');
      bases.push_back("");
    }
    if(nucl(inp_string[cur_i])){
      operation_parsed = true;
      numbers.push_back(1);
      ops.push_back('X');
      bases.push_back(inp_string.substr(cur_i,1));
      //cout<<"detected: X: "<<inp_string.substr(cur_i,1)<<endl;
      cur_i++;
    }
    if(inp_string[cur_i] == '^'){
      operation_parsed = true;
      string inserted_seq;
      while(inp_string[cur_i] == '^'){
	unsigned int ni = cur_i+1;
	while(ni < inp_string.length() && !number(inp_string[ni])){
	  inserted_seq = inserted_seq + inp_string[ni];
	  ni++;
	}
	cur_i = ni;
	if(cur_i + 1 < inp_string.length()){
	  if(inp_string[cur_i] == '0' && inp_string[cur_i+1] == '^'){
	    cur_i++; //shift to the next insertion base
	    //exit(1);
	  }
	}
      }
      numbers.push_back(0);
      ops.push_back('I');
      bases.push_back(inserted_seq);
      //cout<<"detected:"<<inserted_seq<<endl;
    }
    
    if(!operation_parsed){
      cout<<"could not understand mdfield operation: "<<inp_string[cur_i]<<" in "<<inp_string<<endl;
      exit(1);
      //cur_i++;
    }
  }
}
*/
mdfield::mdfield(const StringView& inp_string){
  init(inp_string);
}
void mdfield::init(const StringView& inp_string){
  assert(inp_string[inp_string.length()-1] != char(10));

  if(inp_string.length() == 0)  return;

  //assert(inp_string.length() >= 3);
  //assert(inp_string.substr(0,5) == "MD:Z:"); //the string starts with MD:Z:

  //parse mdfield string
  //unsigned int last_i = 5;
  unsigned int cur_i = 5;
  while(cur_i < inp_string.length()){
    bool operation_parsed = false;
    if(number(inp_string[cur_i])){
      operation_parsed = true;
      string cur_num;
      unsigned int ni = cur_i;
      while(ni<inp_string.length() && number(inp_string[ni])){
	cur_num = cur_num + inp_string[ni];
	ni++;
      }
      cur_i = ni;
      numbers.push_back(atoi(cur_num.c_str()));
      ops.push_back('M');
      //StringView empty(inp_string, 0, 0); //emptry StringView
      //bases.push_back(empty);
      bases.push_back("");
      //cout<<"detected: M:"<<atoi(cur_num.c_str())<<endl;
    }
    if(nucl(inp_string[cur_i])){
      operation_parsed = true;
      numbers.push_back(1);
      ops.push_back('X');
      //StringView strView(inp_string, cur_i, 1);
      //bases.push_back(strView);
      bases.push_back(inp_string.substr(cur_i,1));
      cur_i++;
    }
    if(inp_string[cur_i] == '^'){
      operation_parsed = true;
      string inserted_seq;
      while(inp_string[cur_i] == '^'){
	unsigned int ni = cur_i+1;
	while(ni < inp_string.length() && !number(inp_string[ni])){
	  inserted_seq = inserted_seq + inp_string[ni];
	  ni++;
	}
	cur_i = ni;
	if(cur_i + 1 < inp_string.length()){
	  if(inp_string[cur_i] == '0' && inp_string[cur_i+1] == '^'){
	    cur_i++; //shift to the next insertion base
	    //exit(1);
	  }
	}
      }
      numbers.push_back(inserted_seq.length());
      ops.push_back('I');
      bases.push_back(inserted_seq);
      //cout<<"detected:"<<inserted_seq<<endl;
    }
    
    if(!operation_parsed){
      
      cout<<"could not understand mdfield operation: ["<<inp_string[cur_i]<<"] in "<<inp_string<<endl;
      cout<<"int(inp_string[cur_i]): "<<int(inp_string[cur_i])<<endl;
      exit(1);
      //cur_i++;
    }
  }
}

int SamAlignment::init_cigar(){
  //cigar tmp_cigar(cigar_string);
  //_cigar = tmp_cigar;
  return _cigar.init(cigar_string);
}

void SamAlignment::init_mdfield(){
  //mdfield tmp_mdfield(mdfield_string);
  //_mdfield = tmp_mdfield;
  if(ptr_mdfield != NULL){
    //for(unsigned int i=0; i<ptr_mdfield->ops.size(); i++){
    //  cout<<" "<<ptr_mdfield->ops[i]<<"-"<<ptr_mdfield->numbers[i];
    //}
    //cout<<endl;
    //cout<<sam_string<<endl;
    //assert(false);
    delete ptr_mdfield; ptr_mdfield = NULL;
  }
  ptr_mdfield = new mdfield;
  ptr_mdfield->init(mdfield_string);
}

bool SamAlignment::mapped(){
  return _mapped;
}

SamAlignment& SamAlignment::operator=(const SamAlignment& rhs){
  sam_string = rhs.sam_string;
  flag = rhs.flag;
  _cigar = rhs._cigar;
  if(rhs.ptr_mdfield != NULL){
    if(ptr_mdfield != NULL){
      assert(false);
    }
    ptr_mdfield = new mdfield;
    ptr_mdfield->numbers = rhs.ptr_mdfield->numbers;
    ptr_mdfield->ops = rhs.ptr_mdfield->ops;
    ptr_mdfield->bases = rhs.ptr_mdfield->bases;
  }

  begin = rhs.begin;
  end = rhs.end;
  mapq = rhs.mapq;

  chrom = rhs.chrom;
  strand = rhs.strand;

  qname = rhs.qname;
  seq = rhs.seq;
  qual = rhs.qual;
  cigar_string = rhs.cigar_string;
  mdfield_string = rhs.mdfield_string;

  _valid = rhs._valid;
  _mapped = rhs._mapped;

  qname.ptrStr = &sam_string;
  chrom.ptrStr = &sam_string;
  seq.ptrStr = &sam_string;
  qual.ptrStr = &sam_string;
  cigar_string.ptrStr = &sam_string;
  mdfield_string.ptrStr = &sam_string;
  
  return *this;
}

bool SamAlignment::valid(){
  return _valid;
}

SamAlignment::~SamAlignment(){
  if(ptr_mdfield != NULL){
    delete ptr_mdfield;
    ptr_mdfield = NULL;
    //cout<<"Memory freed"<<endl;
    return;
  }
  
  //cout<<"deleting"<<endl;
}

SamAlignment::SamAlignment(){
  _valid = false;
  strand = '!'; 
  flag = 4;
  ptr_mdfield = NULL;
}

SamAlignment::SamAlignment(const SamAlignment& al){
  //cout<<"copy constructor"<<endl;
  sam_string = al.sam_string;
  flag = al.flag;
  _cigar = al._cigar;

  if(al.ptr_mdfield != NULL){
    ptr_mdfield = NULL;
    /*
    assert(ptr_mdfield == NULL);
    if(ptr_mdfield != NULL){
      delete ptr_mdfield;
    }
    ptr_mdfield = new mdfield;
    ptr_mdfield->numbers = al.ptr_mdfield->numbers;
    ptr_mdfield->ops = al.ptr_mdfield->ops;
    ptr_mdfield->bases = al.ptr_mdfield->bases;
    */
  } else {
    ptr_mdfield = NULL;
  }
  begin = al.begin;
  end = al.end;
  mapq = al.mapq;
 
  chrom = al.chrom;
  strand = al.strand;

  qname = al.qname;
  seq = al.seq;
  qual = al.qual;
  cigar_string = al.cigar_string;
  mdfield_string = al.mdfield_string;

  _valid = al._valid;
  _mapped = al._mapped;

  qname.ptrStr = &sam_string;
  chrom.ptrStr = &sam_string;
  seq.ptrStr = &sam_string;
  qual.ptrStr = &sam_string;
  cigar_string.ptrStr = &sam_string;
  mdfield_string.ptrStr = &sam_string;
}

SamAlignment::SamAlignment(std::string& _sam_string, std::string& error_msg){
  //error is written into error_msg
  sam_string = _sam_string;
  //cout<<endl;
  //cout<<"sam_string: "<<sam_string<<endl;
  strand = '!'; //default it's not defined
  flag = 4;
  ptr_mdfield = NULL;
  if(sam_string.size() == 0){
    _valid = false; 
    return;
  }
  
  if(sam_string[0] == '@' || sam_string[0] == '#'){
    _valid = false;
    return;
  }
  
  char tab_char = char(9);
  //vector<string> sam_string_fields = split(sam_string, tab_char);

  vector<StringView> fields;
  split(sam_string, tab_char, fields);
  
  //cout<<"Fields: "<<fields.size()<<endl;
  //for(unsigned int i=0; i<fields.size(); i++){
  //  cout<<i+1<<": "<<fields[i]<<endl;
  //}
  
  if(fields.size() < 11){
    // fewer fiedls than expected
    //cerr<<"Error in SamAlignment::SamAlignment(std::string& sam_string)"<<endl;
    //cerr<<"too few fields"<<endl;
    //cerr<<"sam_string: "<<sam_string<<endl;
    //cerr<<"expected at least 12 fields, but found "<<fields.size()<<endl;
    char buffer[50];
    sprintf(buffer,"%d", (int) fields.size()); string buffer_string(buffer);
    error_msg = "expected at least 12 fields, but found " + buffer_string; //itoa(fields.size());
    //exit(0);
    _valid = false; 
    _mapped = false;
    return;
  } else { 
    qname = fields[0];
    //cout<<"mapq: "<<fields[4]<<endl;
    mapq = atoi(fields[4].str().c_str());
    //cout<<"interpreted: "<<mapq<<endl;
    seq = fields[9];
    qual = fields[10];

    //assert(qual.size() == seq.size());
    /*
    if(qual.size() != seq.size()){
      cout<<endl<<"Error: qual.size() != seq.size()"<<endl;
      cout<<sam_string<<endl;
      _valid = false;
      _mapped = false;
      return;
    }
    */
    cigar_string = fields[5];

    if(fields.size() >= 16){
      mdfield_string = fields[15];

      if(qual.size() != seq.size()){
	cout<<endl<<"Error: qual.size() != seq.size()"<<endl;
	cout<<sam_string<<endl;
	_valid = false;
	_mapped = false;
	return;
      }
    }
  }
  

  flag = (unsigned int) (atoi(fields[1].str().c_str()));
  if( (flag & 0x4) == 0x4){
    _mapped = false;
    _valid = true;
    return;
  } else {
    _mapped = true;
  }

  if( !(this->mapped()) ){
    _valid = true;
    return;
  } else {
    if(fields.size() >= 11){
      //found the right number of fields
      _valid = true;

      // initialize alignment
      if( (flag & 16) == 16){
	strand = '-';      
      } else {
	strand = '+';
      }
      //cout<<"strand: "<<strand<<endl;
      
      chrom = fields[2];
      string begin_string = fields[3].str();
      begin = atoi(begin_string.c_str());//atoi(fields[3].c_str());
    } else {
      cerr<<sam_string<<endl;
      cerr<<"Flag "<<(flag&0x4)<<" said mapped, but "<<fields.size()<<" fields were found instead of 14, which is a nonsense"<<endl;
      _valid = false;      
    }
  }

  if(_valid && _mapped){
    int res = init_cigar();
    if(res != 0){ 
      _mapped =  false; _valid = false; 
      cout<<"Error in cigar: "<<cigar_string<<endl;
      cout<<"Sam: "<<sam_string<<endl;

      return;
      //cout<<"flag: ["<<fields[1].c_str()<<"] (uint)flag: "<<flag<<" flag&0x4: "<<(flag & 0x4)<<endl;}
    }
    else{
      end = begin + _cigar.reference_length(); 
    }
  }
  
  return;
}

/*
bool SamAlignment::read(const char* _sam_str){
  //sam_string = string(_sam_str);
  sam_string = _sam_str;
  strand = '!'; //default it's not defined
  flag = 4;
  if(sam_string.size() == 0){
    _valid = false; 
    return false;
  }
  
  if(sam_string[0] == '@' || sam_string[0] == '#'){
    _valid = false;
    return false;
  }
  
  //char tab_char = char(9);
  char* ssf = strtok((char*)_sam_str, "\t");

  vector<string> sam_string_fields;
  while(ssf != NULL){
    sam_string_fields.push_back(string(ssf));
    ssf = strtok(NULL, "\t");
  }

  //cout<<"sam string: "<<sam_string<<endl;
  //cout<<"fileds: "<<endl;
  //for(unsigned int i=0; i<sam_string_fields.size(); i++){
  //  cout<<sam_string_fields[i]<<endl;
  //}
  //exit(0);

  if(sam_string_fields.size() < 12){
    // fewer fiedls than expected
    cerr<<sam_string<<endl;
    cerr<<"Error in SamAlignment::SamAlignment(std::string& sam_string)"<<endl;
    cerr<<"too few fields"<<endl;
    cerr<<"sam_string: "<<sam_string<<endl;
    cerr<<"expected at least 12 fields, but found "<<sam_string_fields.size()<<endl;
    _valid = false; 
    return false;
  } else { 
    qname = sam_string_fields[0];
    mapq = atoi(sam_string_fields[4].c_str());
    seq = sam_string_fields[9];
    qual = sam_string_fields[10];

    //assert(qual.size() == seq.size());
    if(qual.size() != seq.size()){
      _valid = false;
    }
    cigar_string = sam_string_fields[5];

    if(sam_string_fields.size() >= 16){
      mdfield_string = sam_string_fields[15];
    }
  }
  

  flag = (unsigned int) (atoi(sam_string_fields[1].c_str()));
  if( (flag & 4) == 4){
    _mapped = false;
  } else {
    _mapped = true;
  }
  
	      
  if( !(this->mapped()) ){
    _valid = true;
    
  } else {
    if(sam_string_fields.size() >= 11){
      //found the right number of fields
      _valid = true;

      //qname = sam_string_fields[0];
      //seq = sam_string_fields[9];
      //qual = sam_string_fields[10];

      // initialize alignment
      if( (flag & 16) == 16){
	strand = '-';      
      } else {
	strand = '+';
      }
      //cout<<"strand: "<<strand<<endl;
      
      chrom = sam_string_fields[2];
      begin = atoi(sam_string_fields[3].c_str());
    } else {
      cerr<<sam_string<<endl;
      cerr<<"Flag "<<(flag&4)<<" said mapped, but "<<sam_string_fields.size()<<" fields were found instead of 14, which is a nonsense"<<endl;
      _valid = false;      
    }
  }
  
  if(_valid && this->mapped()){
    cigar _cigar(cigar_string);
    end = begin + _cigar.reference_length(); //seq.length();
    //cout<<"begin: "<<begin<<" end: "<<end<<endl;
    //mapped_length = _cigar.mapped_length();
  }
  
  return true;
}
*/

/*
bool SamAlignment::output(ofstream& ofstr){
  if(!ofstr.good()) return false;
  char sep = char(9);

  ofstr<<qname<<sep<<flag<<sep<<chrom<<sep<<begin<<sep<<mapq<<sep<<cigar_string;
  ofstr<<sep<<"*"<<sep<<"0"<<sep<<"0"<<sep<<seq<<sep<<qual<<sep<<"*"<<sep<<"*"<<sep<<"*";
  ofstr<<endl;
  return true;
}
*/

/*
bool SamMPAlignment::valid(){
	return _valid;
}

SamMPAlignment::SamMPAlignment(){
  //mapped = false;
  _valid = false;
  strand1 = '!'; //not defined
	strand2 = '!';
  flag = 4;
}

SamMPAlignment::SamMPAlignment(std::string& _sam_string){
	//cout<<"cur_sam_string:"<<endl;
	//cout<<_sam_string<<endl;
	
	if(_sam_string.length() == 0){
		flag=4;
		_valid=false;
		return;
	}
	if(_sam_string[0] == '@'){
		flag=4;
		_valid=false;
		return;
	}

	
	sam_string = _sam_string;
  strand1 = '!'; //default it's not defined
	strand2 = '!';
	
  flag = 4;
  if(sam_string.size() == 0){
    _valid = false; 
    return;
  }
  
  if(sam_string[0] == '@' || sam_string[0] == '#'){
    _valid = false;
    return;
  }
  
  char tab_char = char(9);
  vector<string> sam_string_fields = split(sam_string, tab_char);

  if(sam_string_fields.size() < 12){
    // fewer fiedls than expected
    cerr<<sam_string<<endl;
    cerr<<"Error in SamAlignment::SamAlignment(std::string& sam_string)"<<endl;
    cerr<<"too few fields"<<endl;
    _valid = false; 
    return;
  } else { 
    qname = sam_string_fields[0];
    mapq = atoi(sam_string_fields[4].c_str());
    seq = sam_string_fields[9];
    qual = sam_string_fields[10];
    cigar_string = sam_string_fields[5];
		_valid = true;	
  }


  flag = (unsigned int) (atoi(sam_string_fields[1].c_str()));
  if( (flag & 4) == 4){
    mapped1 = false;
  } else {
    mapped1 = true;
  }

  if( (flag & 8) == 8){
    mapped2 = false;
  } else {
    mapped2 = true;
  }
  
  
  if( (flag & 0x10) == 0x10){ // + strand
    strand1 = '-';      
  } else {
    strand1 = '+';
  }
  
  if( (flag & 0x0020) == 0x0020){
    strand2 = '-';
  } else {
    strand2 = '+';
  }
  
  chrom1 = sam_string_fields[2];
  begin1 = atoi(sam_string_fields[3].c_str());
  
  chrom2 = sam_string_fields[6];
  if(chrom2 == "="){
    chrom2 = chrom1;
  }
  begin2 = atoi(sam_string_fields[7].c_str());
  
  if( (flag & 0x40) == 0x40){
    string tmp_chrom = chrom1;
    int tmp_begin = begin1;	
    int tmp_strand = strand1;
    
    chrom1 = chrom2;
    begin1 = begin2;
    strand1 = strand1;
    
    chrom2 = tmp_chrom;
    begin2 = tmp_begin;
    strand2 = tmp_strand;
    
    flag &= ~(1 << 6); // now 0x40 is set to 0
    flag |= (1 << 7); // now 0x80 is set to 1
    assert(flag & 0x80 == 0x80);
  }
  
  
  
  //cigar _cigar(cigar_string);
  
  //cout<<chrom1<<" "<<strand1<<" "<<begin1<<endl;
  //cout<<chrom2<<" "<<strand2<<" "<<begin2<<endl;
  //exit(1);
  //end = begin + _cigar.reference_length();
  return;

}
//bool SamMPAlignment::mapped(){
//	return _mapped;
//}
*/
