#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <exception>
#include "assert.h"

#include <stdio.h>
#include "params.h"
#include "string_utils.h"
#include "math_functions.h"

double kernel(double x, double bw){
  if(x<=-bw || x>=bw) return 0;
  //return  exp( - x*x/(2.0*bw*bw) ) * 0.3989423 / bw;
  //return  exp( - x*x/(2.0*bw*bw) );// * 0.3989423 / bw;
  
  double u=1.0-x*x/(bw*bw); 
  return u*u*u*35/(32*bw);
}

using std::ifstream;
using std::ofstream;

using std::istringstream;
using std::cerr;
using std::endl;
using std::cout;

int main(int argc, char* argv[]){

  params pars(argc, argv);
  pars.require("dist_file", "dist_file", STRING_TYPE);
  pars.require("output_file", "output_file", STRING_TYPE);

  pars.optional("bw","bw","5",DOUBLE_TYPE);
  pars.optional("field", "field", "1", INT_TYPE);
  pars.optional("normalize", "yes/no", "no", STRING_TYPE);


  if(!pars.enforce()){
    exit(1);
  }
   
  double bw = pars.get_double_value("bw");
  string dist_fname = pars.get_string_value("dist_file");
  int field = pars.get_int_value("field");

  string output_fname = pars.get_string_value("output_file");
  string normalize_string = pars.get_string_value("normalize");

  assert(normalize_string == "yes" || normalize_string == "no");

  bool normalize;
  if(normalize_string == "yes") normalize = true;
  else normalize = false;

  cout<<"bw: "<<bw<<endl;
  cout<<"dist_file: "<<dist_fname<<endl;
  cout<<"field to smooth: "<<field<<endl;

  ifstream dist_ifstr(dist_fname.c_str());
  if(!dist_ifstr.good()){
    cerr<<"Bad file name: "<<dist_fname<<endl;
    exit(1);
  }

  vector<string> entries;
  vector<int> positions;
  vector<int> raw_values;
  vector<int> pos; //positions
  vector<double> values;

  char delim = '\t';

  while(dist_ifstr.good()){
    string cur_line;
    getline(dist_ifstr, cur_line);
    if(dist_ifstr.good()){
      vector<string> cur_line_fields = split(cur_line, delim);
      
      int cur_pos;
      double cur_value;
      
      cur_pos = atoi(cur_line_fields.at(0).c_str());
      positions.push_back(cur_pos);
      cur_value = atof(cur_line_fields.at(field).c_str());
      raw_values.push_back(cur_value);
      
      entries.push_back(cur_line);
      pos.push_back(cur_pos);
      values.push_back(cur_value);
      
      //cout<<"pos: "<<cur_pos<<" value: "<<cur_value<<endl;
    }
  }

  dist_ifstr.close();

  //now smooth

  vector<double> smoothed_values(values.size(), 0);
  vector<double> baseline(values.size(), 1.0);

  vector<double> smoothed_baseline(values.size(), 0);
  for(int i=0; i<(int)values.size(); i++){
    for(int j=i-bw; j<=i+bw; j++){
      if(j>=0 && j<(int)smoothed_baseline.size()){
	smoothed_baseline[i] += baseline[j]*kernel((double)(i-j), bw);
      }
    }
  }
  for(int i=0; i<(int)values.size(); i++){
    smoothed_values.at(i) = 0;
    int half_win = (int)(bw);
    for(int j=(int)i-half_win; j<= (int)i+half_win; j++){
      if(j>=0 && j< (int)smoothed_values.size()){
	smoothed_values.at(i) += values.at(j)*kernel((double)(i-j), bw);
      }
    }
  }
  for(int i=0; i<(int) smoothed_values.size(); i++){
    smoothed_values[i] = smoothed_values[i]/smoothed_baseline[i];
    //smoothed_values[i] = smoothed_baseline[i];
  }

  //finding peaks / valleys

  string peak_fname = output_fname + ".peak_pos.txt";
  ofstream peak_ofstr(peak_fname.c_str());
  assert(peak_ofstr.good());

  cout<<"peaks: "<<endl<<endl;
  int peak_counter = 0;

  for(int i=1; i<(int) smoothed_values.size()-1; i++){
    if(smoothed_values[i] > smoothed_values[i-1] && smoothed_values[i] > smoothed_values[i+1]){
      //check if the peak is the local maximum in a nucleosome-sized window
      double cur_max = smoothed_values[i];
      bool local_maximum = true;
      for(int j=max(0,i-75); j<min(smoothed_values.size(), i+75); j++){
	//cout<<i<<" "<<j; //<<endl;
	if(smoothed_values.at(j) > cur_max){
	  local_maximum = false;
	  //cout<<" *";
	}
      }
      //cout<<"Max: "<<i<<endl;
      
      if(local_maximum){
	//check that the value is sufficiently above the local trough
	double left_trough_value = cur_max;
	double right_trough_value = cur_max;

	for(int j=max(i-120,0); j<i; j++){
	  if(smoothed_values[j] < left_trough_value) left_trough_value = smoothed_values[j];
	}
	for(int j=i; j<min(i+120, (int)smoothed_values.size()); j++){
	  if(smoothed_values[j] < right_trough_value) right_trough_value = smoothed_values[j];
	}
	//cout<<"left_trough: "<<left_trough_value<<endl;
	//cout<<"right_trough: "<<right_trough_value<<endl;
	//cout<<"cur_max: "<<cur_max<<endl;
	if(cur_max / left_trough_value >= 1.01 &&
	   cur_max / right_trough_value >= 1.01){
	  cout<<pos[i]<<" "<<smoothed_values[i]<<endl;
	  //cout<<pos[i]<<" "<<smoothed_values[i]<<endl;
	  peak_ofstr<<peak_counter<<"\t"<<pos[i]<<endl;
	  peak_counter++;
	}
      }
    }
  }
  peak_ofstr.close();
  
  string trough_fname = output_fname + ".trough_pos.txt";
  ofstream trough_ofstr(trough_fname.c_str());
  assert(trough_ofstr.good());
  
  cout<<"trouphs: "<<endl<<endl;
  int trough_counter=0;
  
  
  for(int i=50; i<(int)smoothed_values.size()-1; i++){
    if(smoothed_values[i] < smoothed_values[i-1] && smoothed_values[i] < smoothed_values[i+1]){
      //cout<<pos[i]<<" "<<smoothed_values[i]<<endl;
      //check if this is the trough in the nucleosome-sized window
      bool cur_min = true;
      double min_value = smoothed_values[i];
      for(int j=max(0, i-75); j<min(i+75, smoothed_values.size()); j++){
	if(smoothed_values[j] < min_value) cur_min = false;
      }
      if(cur_min){
	cout<<trough_counter<<" "<<i<<endl;
	trough_ofstr<<trough_counter<<"\t"<<pos[i]<<endl;
	trough_counter++;
      }
    }
  }
  trough_ofstr.close();

  double norm = 1;
  
  if(normalize){
    double sum = 0;
    for(unsigned int i=0; i<smoothed_values.size(); i++){
      sum += smoothed_values[i];
    }
    norm = sum / ((double)smoothed_values.size());
  }

  ofstream out_str(output_fname.c_str());
  for(unsigned int i=0; i<smoothed_values.size(); i++){
    //out_str<<entries.at(i)<<" "<<smoothed_values.at(i)<<endl;
    out_str<<positions.at(i)<<"\t"<<raw_values.at(i)<<"\t"<<smoothed_values.at(i)/norm<<endl;
  }

  out_str.close();

  return 0;
}
