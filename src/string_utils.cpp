#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <cstring>

#include <assert.h>
#include "string_utils.h"

using std::vector;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;

const char* StringView::c_str() const{
  return (ptrStr->substr(begin, _length)).c_str();
}

string StringView::str() const {
  assert(ptrStr != NULL);
  return ptrStr->substr(begin, _length);
}

const char StringView::at(int i) const{
  assert(i<_length); return (*ptrStr)[begin+i];
}

bool operator==(const StringView& lhs, const char* rhs){
  for(int i=0; i<(int)strlen(rhs); i++){
    if(lhs[i] != rhs[i]) return false;
  }
  return true;
}

bool operator==(const StringView& lhs, const string& rhs){
  if(lhs.length() != rhs.length()) return false;

  for(int i=0; i<(int)lhs.length(); i++){
    if(lhs[i] != rhs[i]) return false;
  }
  return true;
}

bool operator==(const StringView& lhs, const StringView& rhs){
  if(lhs._length != rhs._length) return false;
  
  if(lhs.ptrStr == rhs.ptrStr && lhs.begin == rhs.begin && lhs._length == rhs._length) return true;

  //lengths are the same
  for(int i=0; i<lhs._length; i++){
    if(lhs[i] != rhs[i]) return false;
  }
  return true;
}

bool operator<(const StringView& lhs, const StringView& rhs) {
  string lhs_string = lhs.str(); 
  string rhs_string = rhs.str();

  return lhs_string < rhs_string;
}

//#include "utils.h"
/*
std::vector<string> split(const string& inp_string, const string& regexp){
  std::vector<string> res;
  if(regexp[0] == '[' && regexp[regexp.length()-1] == ']'){
    //split by any character within []
    
    
    if(inp_string.length() == 0)  return res;
    
    unsigned int last_i = 0;
    for(unsigned int i=0; i<inp_string.length(); i++){
      //if(inp_string.substr(i,sep.length()) == sep){
      if(inp_string[i] == sep){
	if(i>last_i){ 
	  res.push_back(inp_string.substr(last_i, i-last_i));
	}
	last_i = i+1;
	//}
      }
    }
    if(last_i < inp_string.length()) 
      res.push_back(inp_string.substr(last_i, inp_string.length()-last_i));
    
    if(res.size() == 0 && inp_string[0] != sep) res.push_back(inp_string);
  }

  return res;
}
*/

std::vector<string> split(const string& inp_string, char sep){
  std::vector<string> res;
  if(inp_string.length() == 0)  return res;

  unsigned int last_i = 0;
  for(unsigned int i=0; i<inp_string.length(); i++){
    if(inp_string[i] == sep){
      if(i>last_i){ 
	res.push_back(inp_string.substr(last_i, i-last_i));
      }
      last_i = i+1;
    }
  }
  if(last_i < inp_string.length()) 
    res.push_back(inp_string.substr(last_i, inp_string.length()-last_i));

  if(res.size() == 0 && inp_string[0] != sep) res.push_back(inp_string);
  return res;
}

void split(const string& inp_string, char sep, vector<string>& res){
  //std::vector<string> res;
  if(inp_string.length() == 0)  return;

  unsigned int last_i = 0;
  for(unsigned int i=0; i<inp_string.length()/*-1*/; i++){
    if(inp_string[i] == sep){
      if(i>last_i){ 
	res.push_back(inp_string.substr(last_i, i-last_i));
      }
      last_i = i+1;
      //}
    }
  }
  if(last_i < inp_string.length()) 
    res.push_back(inp_string.substr(last_i, inp_string.length()-last_i));

  if(res.size() == 0 && inp_string[0] != sep) res.push_back(inp_string);
  return;
}

void split(string& inp_string, char sep, vector<StringView>& res){
  if(inp_string.length() == 0)  return;

  unsigned int last_i = 0;
  for(unsigned int i=0; i<inp_string.length()/*-1*/; i++){
    if(inp_string[i] == sep){
      if(i>last_i){ 
	StringView stringView(inp_string, last_i, (int) i-last_i);
	res.push_back(stringView);
      }
      last_i = i+1;
    }
  }
  if(last_i < inp_string.length()) {
    StringView stringView(inp_string, (int)last_i, (int)(inp_string.length()-last_i));
    res.push_back(stringView);
  }

  if(res.size() == 0 && inp_string[0] != sep){
    StringView stringView(inp_string, 0, (int) inp_string.length());
  }
  return;
}

void split(const StringView& inpStringView, char sep, std::vector<StringView>& res){
  if(inpStringView.length() == 0)  return;

  unsigned int last_i = 0;
  for(unsigned int i=0; i<inpStringView.length(); i++){
    if(inpStringView[i] == sep){
      if(i>last_i){
        StringView stringView(inpStringView, last_i, (int) i-last_i);
        res.push_back(stringView);
      }
      last_i = i+1;
    }
  }
  if(last_i < inpStringView.length()) {
    StringView stringView(inpStringView, (int)last_i, (int)(inpStringView.length()-last_i));
    res.push_back(stringView);
  }

  if(res.size() == 0 && inpStringView[0] != sep){
    StringView stringView(inpStringView, 0, (int) inpStringView.length());
  }
  return;
}


unsigned int get_string_count(std::ifstream& ifs){
  using std::ios;

  assert(ifs.good());
  std::ios::pos_type cur_pos = ifs.tellg();
  ifs.seekg(0, ios::beg);

  string dummy_string;
  unsigned int string_counter = 0;
  while(ifs.good()){
    getline(ifs,dummy_string);
    if(ifs.good()) string_counter++;    
  }

  ifs.clear();
  ifs.seekg(cur_pos);
  return string_counter;
}

void chomp(string& str){
  //erases bad characters
  bool continue_chomp = true;

  while(continue_chomp){
    switch(str[str.length()-1]){
    case char(13):{ str.erase(str.length()-1,1); break;}
    case char(10):{ str.erase(str.length()-1,1); break;}
    default:{ continue_chomp = false; break; }
    }
  }  
}


std::ios::pos_type perc_file_length(std::ifstream& ifs){
  //assert(ifs.good());
  if(!ifs.good()){ cerr<<"bad file in perc_file_length"<<endl; std::terminate();}
  std::ios::pos_type cur_pos = ifs.tellg();
  ifs.seekg(0, std::ios::end);
  std::ios::pos_type _length = ifs.tellg();
  ifs.seekg(cur_pos);
  return (_length - _length%100) / 100;
}


int k_pow(int k){
  //returns power of 1000
  assert(k >= 0);
  if(k==0) return 1;
  if(k==1) return 1000;
  if(k==2) return 1000000;
  if(k==3) return 1000000000;
  assert(false); //you are asking for too much!
  return -1;
}

long long_k_pow(int k){
  assert(k >= 0);
  if(k==0) return 1;
  if(k==1) return 1000;
  if(k==2) return 1000000;
  if(k==3) return 1000000000;
  //if(k==4) return 1000000000000;
  //if(k==5) return 1000000000000000;
  cerr<<"k = "<<k<<" is too big to handle"<<endl;
  assert(false); //you are asking for too much!
  return -1;
}

std::string commify(int i){
  std::string res;
  
  bool __continue = true;
  int cur_pow = 1; // power of 1000

  int prev_div = 0;
  while(__continue){
    
    int cur_div = long_k_pow(cur_pow);
    if(i < cur_div){
      __continue = false;
    }
    int cur_num = (i%cur_div - prev_div ) / long_k_pow(cur_pow-1);
    std::stringstream tmp1;
    tmp1<<cur_num;    
    std::string cur_num_string =tmp1.str();
    

    if(cur_pow > 1){      
      res = "," + res;
    }
    if(i>cur_div){ //there will be numbers in front, so ok to add zeroes
      if(cur_num_string.length() == 2){
	cur_num_string = "0" + cur_num_string;
      }
      if(cur_num_string.length() == 1){
	cur_num_string = "00" + cur_num_string;
      }
      if(cur_num_string.length() == 0){
	cur_num_string = "000";// + cur_num_string;
      }
    }

    res = cur_num_string + res;
    cur_pow++;
    prev_div = cur_num;
  }
  return res;
}

std::string commify_long(long i){
  std::string res;
  
  bool __continue = true;
  int cur_pow = 1; // power of 1000

  long prev_div = 0;
  while(__continue){
    
    long cur_div = long_k_pow(cur_pow);
    if(i < cur_div){
      __continue = false;
    }
    long cur_num = (i%cur_div - prev_div ) / long_k_pow(cur_pow-1);
    std::stringstream tmp1;
    tmp1<<cur_num;    
    std::string cur_num_string =tmp1.str();
    

    if(cur_pow > 1){      
      res = "," + res;
    }
    if(i>cur_div){ //there will be numbers in front, so ok to add zeroes
      if(cur_num_string.length() == 2){
	cur_num_string = "0" + cur_num_string;
      }
      if(cur_num_string.length() == 1){
	cur_num_string = "00" + cur_num_string;
      }
      if(cur_num_string.length() == 0){
	cur_num_string = "000";// + cur_num_string;
      }
    }

    res = cur_num_string + res;
    cur_pow++;
    prev_div = cur_num;
  }
  return res;
}

std::string int2string(int i){
  std::string s;
  std::stringstream out;
  out<<i;
  s = out.str();
  return s;  
}
