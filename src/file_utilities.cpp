#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "assert.h"
//#include "string_utils.h"
#include "file_utilities.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;

std::vector<std::string> split(const std::string& inp_string, char sep);


bool mkdir(const string& path){
  int status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  return status; //true;                                                                                                            
}

bool directory_exists(const string& path){
  struct stat sb;
  if(stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) return true;
  return false;
}

bool create_file_path(const string& file_path){
  char delim = '/';
  vector<string> path_fields = split(file_path, delim);
  if(path_fields.size() == 1){
    cout<<"Warning, only one field in the path, no action taken"<<endl;
    return true;
  }
  string cur_path = "";
  if(file_path.length() > 0){ if(file_path[0] == '/') cur_path="/";}
  for(unsigned int i=0; i<path_fields.size()-1; i++){
    cur_path = cur_path + path_fields[i] + "/";
    if(!directory_exists(cur_path)){
      cout<<"Creating directory "<<cur_path<<endl;
      mkdir(cur_path);
      if(!directory_exists(cur_path)){
        cout<<"Failed to create directory "<<cur_path<<endl;
	exit(1);
      }
    } else {
      cout<<cur_path<<" already exists"<<endl;
    }
  }
  return true;
}

bool file_exists(const std::string& filename){
  struct stat buf;
  if (stat(filename.c_str(), &buf) != -1)
    {
      return true;
    }
  return false;
}

string get_file_path(const string& file_path){
  //returns a path to directory containing a provided file path
  char delim = '/';
  vector<string> path_fields = split(file_path, delim);
  if(path_fields.size() == 1){
    cout<<"Warning, only one field in the path, no action taken"<<endl;
    cout<<"path: "<<file_path<<endl;
    exit(1);
  }
  string cur_path = "";
  for(unsigned int i=0; i<path_fields.size()-1; i++){
    if(i==0){
      if(path_fields[i] != ".") cur_path = "/";
    }
    cur_path = cur_path + path_fields[i] + "/";
    if(!directory_exists(cur_path)){
      cout<<"No such directory "<<cur_path<<endl;
      exit(1);
    }
  }
  return cur_path;
}
