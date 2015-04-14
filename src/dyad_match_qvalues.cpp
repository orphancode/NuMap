#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <exception>
#include <stdio.h>
#include <assert.h>

#include "params.h"
#include "string_utils.h"
#include "genome_table.h"
#include "math_functions.h"
#include "paths.h"
#include "file_utilities.h"
#include <algorithm>    // std::sort

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
};

struct genome_coord_pair{
  string chrom;
  int pos1;
  int pos2;
  double value1;
  double value2;
  int dist;
  double qvalue;
  string entry;
};

bool dist_less(const genome_coord_pair& l, const genome_coord_pair& r){
  return l.dist < r.dist;
}

bool coord_less(const genome_coord_pair& l, const genome_coord_pair& r){
  if(l.chrom != r.chrom) return l.chrom < r.chrom;
  else return l.pos1 < r.pos1;
}

struct region{
  string chrom;
  int begin;
  int end;
};

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("dyad_match_file", "dyad_match_file", STRING_TYPE);
  pars.require("mock_dyad_match_file", "mock_dyad_match_file", STRING_TYPE);
  pars.require("output_file", "output_file", STRING_TYPE);

  if(!pars.enforce()){
    exit(1);
  }

  string dyad_match_fname = pars.get_string_value("dyad_match_file");
  string mock_dyad_match_fname = pars.get_string_value("mock_dyad_match_file");
  string output_fname = pars.get_string_value("output_file");

  string line;
  char delim = '\t';

  ifstream ifstr1(dyad_match_fname.c_str());

  vector<genome_coord_pair> dyad_matches;

  while(ifstr1.good()){
    getline(ifstr1, line);
    vector<string> line_fields = split(line, delim);
    if(line_fields.size() == 5){
      genome_coord_pair cur_match;
      cur_match.chrom = line_fields[0];
      cur_match.pos1 = atoi(line_fields[1].c_str());
      cur_match.pos2 = atoi(line_fields[2].c_str());
      cur_match.value1 = atof(line_fields[3].c_str());
      cur_match.value2 = atof(line_fields[4].c_str());

      cur_match.entry = line;
      cur_match.dist = abs(cur_match.pos1 - cur_match.pos2);

      //cout<<cur_match.dist<<endl;

      dyad_matches.push_back(cur_match);
    }
  }

  ifstr1.close();
  sort(dyad_matches.begin(), dyad_matches.end(), dist_less);

  cout<<"Read "<<dyad_matches.size()<<" dyad matches"<<endl;

  ifstream ifstr2(mock_dyad_match_fname.c_str());

  vector<genome_coord_pair> mock_dyad_matches;

  while(ifstr2.good()){
    getline(ifstr2, line);
    vector<string> line_fields = split(line, delim);
    if(line_fields.size() == 5){
      genome_coord_pair cur_match;
      cur_match.chrom = line_fields[0];
      cur_match.pos1 = atoi(line_fields[1].c_str());
      cur_match.pos2 = atoi(line_fields[2].c_str());

      cur_match.entry = line;
      cur_match.dist = abs(cur_match.pos1 - cur_match.pos2);
      //cout<<cur_match.dist<<endl;

      mock_dyad_matches.push_back(cur_match);
    }
  }

  ifstr2.close();

  cout<<"Read "<<mock_dyad_matches.size()<<" mock dyad matches"<<endl;

  sort(mock_dyad_matches.begin(), mock_dyad_matches.end(), dist_less);

  string curve_fname = output_fname + ".curve";
  ofstream curve_ofstr(curve_fname.c_str());

  //  int j=0;
  int j_mock = 0; //mock_dyad_matches.size()-1;
  int ii = 0; //dyad_matches.size()-1;

  for(int i=0; i<=100; i++){
    int cur_dist = i;

    bool _continue = true;
    while(_continue){
      if(ii>=(int)dyad_matches.size()){ _continue = false; }
      else{
	//if(cur_dist <= dyad_matches[ii].dist){
	if(dyad_matches[ii].dist <= cur_dist){
	  ii++;
	}
	else _continue = false;
      }
    }

    _continue = true;
    while(_continue){

      if(j_mock >= (int)mock_dyad_matches.size()){ _continue = false; }
      else{
	//if(cur_dist <= mock_dyad_matches[j_mock].dist){
	if(mock_dyad_matches[j_mock].dist <= cur_dist){
	  j_mock++;
	} else _continue = false;
      }
      /*
      if(j_mock<(int) mock_dyad_matches.size()){
	if(cur_dist < mock_dyad_matches[j_mock].dist) _continue = false;

	else j_mock++;
      } else{ _continue = false; }
      */
    }

    /*
    _continue = true;
    while(_continue){
      if(j<(int) dyad_matches.size()){
	if(cur_dist < dyad_matches[j].dist) _continue = false;
	else j++;
      } else { _continue = false; }
      }*/

    double cur_qvalue;
    if(j_mock>=0)
      //cur_qvalue = ((double)j_mock)/((double)j);
      cur_qvalue = ((double)j_mock)/((double)ii);
    else
      assert(false); //need to add more code here
      //cur_qvalue = 1;
      //else cur_qvalue = 1/((double)dyad_matches.size()+1);

    curve_ofstr<<cur_dist<<"\t"<<cur_qvalue<<endl;
  }

  curve_ofstr.close();

  exit(1);

  ofstream ofstr(output_fname.c_str());
  assert(ofstr.good());

  j_mock = 0;
  ii=0;

  for(unsigned int i=0; i<dyad_matches.size(); i++){
    int cur_dist = dyad_matches[i].dist;

    bool _continue = true;
    while(_continue){
      if(ii>=(int)dyad_matches.size()) _continue = false;
      if(cur_dist >= dyad_matches[ii].dist){
	ii++;
      }
      else _continue = false;
    }

    _continue = true;

    while(_continue){
      if(j_mock<(int) mock_dyad_matches.size()){
	if(cur_dist < mock_dyad_matches[j_mock].dist) _continue = false;
	else j_mock++;
      } else{ _continue = false; }
    }

    /*
    _continue = true;
    while(_continue){
      if(j<(int) dyad_matches.size()){
	if(cur_dist < dyad_matches[j].dist) _continue = false;
	else j++;
      } else { _continue = false; }
    }
    */

    double cur_qvalue;
    if(j_mock!=0) cur_qvalue = ((double)mock_dyad_matches.size() - (double)j_mock)/((double)dyad_matches.size()-(double)ii);
    else cur_qvalue = 1.0/((double)(dyad_matches.size()+1));
    //if(j_mock<(int) mock_dyad_matches.size())
      //cur_qvalue = (double)(((int)mock_dyad_matches.size()-j))/((double)dyad_matches.size());
      //cur_qvalue = ((double)j)/((double)dyad_matches.size());
    //if(j!=0)
    //  cur_qvalue = ((double)j_mock)/((double)j);
    //else
    //  cur_qvalue = 1;
      //else cur_qvalue = 1/((double)dyad_matches.size()+1);
    dyad_matches[i].qvalue = cur_qvalue;
    //curve_ofstr<<cur_dist<<"\t"<<cur_qvalue<<endl;

  }

  sort(dyad_matches.begin(), dyad_matches.end(), coord_less);
  for(unsigned int i=0; i<dyad_matches.size(); i++){
    ofstr<<dyad_matches[i].chrom<<"\t"<<dyad_matches[i].pos1<<"\t"<<dyad_matches[i].pos2;
    ofstr<<"\t"<<dyad_matches[i].value1<<"\t"<<dyad_matches[i].value2<<"\t"<<dyad_matches[i].qvalue<<endl;
  }

  return 0;
}
