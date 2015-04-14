#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <exception>
#include <stdio.h>
#include <assert.h>
#include <algorithm>

#include "params.h"
#include "string_utils.h"
#include "genome_table.h"
//#include "math_functions.h"
//#include "paths.h"

using std::ifstream;
using std::ofstream;

using std::istringstream;
using std::cerr;
using std::endl;
using std::cout;

struct genome_coord{
  string chrom;
  int pos;
  double value;
  double helper1;
  //bool value_loaded;
  string entry;
};

bool coord_less(const genome_coord& l, const genome_coord& r){
  if(l.chrom == r.chrom) return l.pos < r.pos;
  else return l.chrom < r.chrom;
}

bool value_less(const genome_coord& l, const genome_coord& r){
  return l.value < r.value;
}

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("dyad_calls_file", "dyad_calls_file", STRING_TYPE);
  pars.require("mock_dyad_calls_file", "mock_dyad_calls_file", STRING_TYPE);
  pars.require("output_file", "output_file", STRING_TYPE);

  //pars.optional("max_dist", "max_dist", "3000", INT_TYPE);
  //pars.optional("peak_value", "peak_value", "0.3", DOUBLE_TYPE);
  //pars.optional("trough_drop", "trough_drop", "0.1", DOUBLE_TYPE);
  //pars.optional("mock", "yes/no", "no", STRING_TYPE);
  
  if(!pars.enforce()){
    exit(1);
  }

  string dyad_calls_fname = pars.get_string_value("dyad_calls_file");
  string mock_dyad_calls_fname = pars.get_string_value("mock_dyad_calls_file");
  string output_fname = pars.get_string_value("output_file");

  ofstream ofstr(output_fname.c_str());
  assert(ofstr.good());

  //read the positions
  vector<genome_coord> dyad_calls;
  ifstream dyad_calls_ifstr(dyad_calls_fname.c_str());
  assert(dyad_calls_ifstr.good());
  string cur_line;
  char delim = '\t';
  while(dyad_calls_ifstr.good()){
    getline(dyad_calls_ifstr, cur_line);
    if(dyad_calls_ifstr.good()){
      vector<string> cur_line_fields = split(cur_line, delim);
      assert(cur_line_fields.size() >= 2);
      if(cur_line_fields.size() >= 2){
	genome_coord cur_dyad_call;
	cur_dyad_call.chrom = cur_line_fields[0];
	cur_dyad_call.pos = atoi(cur_line_fields[1].c_str());
	cur_dyad_call.value = atof(cur_line_fields[2].c_str());
	dyad_calls.push_back(cur_dyad_call);
      }
    }
  }
  dyad_calls_ifstr.close();
  cout<<"Read "<<dyad_calls.size()<<" dyad calls"<<endl;
  sort(dyad_calls.begin(), dyad_calls.end(), value_less);

  vector<genome_coord> mock_dyad_calls;
  ifstream mock_dyad_calls_ifstr(mock_dyad_calls_fname.c_str());
  assert(mock_dyad_calls_ifstr.good());

  while(mock_dyad_calls_ifstr.good()){
    getline(mock_dyad_calls_ifstr, cur_line);
    if(mock_dyad_calls_ifstr.good()){
      vector<string> cur_line_fields = split(cur_line, delim);
      assert(cur_line_fields.size() >= 2);
      if(cur_line_fields.size() >= 2){
	genome_coord cur_dyad_call;
	cur_dyad_call.chrom = cur_line_fields[0];
	cur_dyad_call.pos = atoi(cur_line_fields[1].c_str());
	cur_dyad_call.value = atof(cur_line_fields[2].c_str());
	mock_dyad_calls.push_back(cur_dyad_call);
      }
    }
  }
  mock_dyad_calls_ifstr.close();

  cout<<"Read "<<mock_dyad_calls.size()<<" mock_dyad calls"<<endl;
  sort(mock_dyad_calls.begin(), mock_dyad_calls.end(), value_less);
  
  unsigned int j=0; //to track position within the mock dyad call list
  unsigned int ii=0; //to track position in dyad_calls
  int last_ii = -1;
  for(unsigned int i=0; i<dyad_calls.size(); i++){
    //cout<<i<<"\t"<<ii<<"\t"<<j<<endl;
    //if(i%10000 == 0){
    //cout<<"\r"<<((double i))/
    // }
    double cur_value = dyad_calls[i].value;
    bool _continue = true;
    while(_continue){
      if(ii<dyad_calls.size()){
        if(dyad_calls[ii].value >= cur_value) _continue = false;
        else ii++;
      } else {
        _continue = false;
      }
    }
    //move j
    
    //bool _continue = true;
    _continue = true;
    while(_continue){
      if(j<mock_dyad_calls.size()){
	if(mock_dyad_calls[j].value >= cur_value) _continue = false;
        else j++;
      } else {
	_continue = false;
	//no need to move, already at the end
      }
    }
    double cur_qvalue;
    if( j!= mock_dyad_calls.size())
      //cur_qvalue = ((double)((int)(mock_dyad_calls.size())-(int)(j)))/((double)dyad_calls.size());
      cur_qvalue = ((double)((int)(mock_dyad_calls.size())-(int)(j)))/((double)dyad_calls.size()-(int)ii);
    else{
      //cur_qvalue = 1/((double)dyad_calls.size()+1);
      if(last_ii == -1) last_ii = ii;
      cur_qvalue = 1.0/((double)dyad_calls.size()-(int)last_ii);
    }
    dyad_calls[i].helper1 = cur_qvalue;
  }

  j=0;
  ii=0;
  last_ii = 0;

  string curve_fname = output_fname + ".curve";
  ofstream curve_ofstr(curve_fname.c_str());



  for(int i=0; i<=1000; i++){
    double cur_value = ((double)i)/1000;

    bool _continue = true;

    while(_continue){
      if(ii<dyad_calls.size()){
        if(dyad_calls[ii].value >= cur_value) _continue = false;
        else ii++;
      } else {
        _continue = false;
      }
    }

    _continue = true;
    while(_continue){
      if(j<mock_dyad_calls.size()){
	if(mock_dyad_calls[j].value >= cur_value) _continue = false;
        else j++;
      } else {
	_continue = false;
	//no need to move, already at the end
      }
    }
    double cur_qvalue;
    if( j!= mock_dyad_calls.size())
      cur_qvalue = ((double)((int)(mock_dyad_calls.size())-(int)(j)))/((double)dyad_calls.size()-(int)ii);
    else{
      if(last_ii == -1)last_ii= ii;
      cur_qvalue = 1.0/((double)dyad_calls.size()-(int)last_ii);
    }
    curve_ofstr<<cur_value<<"\t"<<cur_qvalue<<endl;
  }
  curve_ofstr.close();

  sort(dyad_calls.begin(), dyad_calls.end(), coord_less);
  

  for(unsigned int i=0; i<dyad_calls.size(); i++){
    ofstr<<dyad_calls[i].chrom<<"\t"<<dyad_calls[i].pos<<"\t"<<dyad_calls[i].value<<" "<<dyad_calls[i].helper1<<endl;
  }
  ofstr.close();

  return 0;
}
