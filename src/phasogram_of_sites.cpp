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
#include "math_functions.h"
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
  bool value_loaded;
};

bool coord_less(const genome_coord& l, const genome_coord& r){
  if(l.chrom == r.chrom) return l.pos < r.pos;
  else return l.chrom < r.chrom;
}

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("positions_file", "positions_file", STRING_TYPE);
  pars.require("output_file", "outptu_file", STRING_TYPE);
  pars.optional("max_dist", "max_dist", "3000", INT_TYPE);
  
  if(!pars.enforce()){
    exit(1);
  }

  string genomic_sites_fname = pars.get_string_value("positions_file");
  string output_fname = pars.get_string_value("output_file");
  int max_dist = pars.get_int_value("max_dist");

  //string gt_fname = analysis_path + genome_table_suffix;
  //genome_table gt(gt_fname);

  ofstream ofstr(output_fname.c_str());
  assert(ofstr.good());

  //read the positions
  vector<genome_coord> genomic_sites;
  ifstream genomic_sites_ifstr(genomic_sites_fname.c_str());
  assert(genomic_sites_ifstr.good());
  string cur_line;
  char delim = '\t';
  while(genomic_sites_ifstr.good()){
    getline(genomic_sites_ifstr, cur_line);
    if(genomic_sites_ifstr.good()){
      vector<string> cur_line_fields = split(cur_line, delim);
      //assert(cur_line_fields.size() >= 2);
      if(cur_line_fields.size() >= 2){
	genome_coord new_site;
	new_site.chrom = cur_line_fields[0];
	new_site.pos = atoi(cur_line_fields[1].c_str());
	if(cur_line_fields.size() >= 3){ 
	  //new_site.value = atof(cur_line_fields[2].c_str());
	  //new_site.value_loaded = true;
	}
	genomic_sites.push_back(new_site);
      }
    }
  }
  genomic_sites_ifstr.close();

  cout<<"Read "<<genomic_sites.size()<<" genomic sites"<<endl;
  sort(genomic_sites.begin(), genomic_sites.end(), coord_less);

  vector<int> phase(max_dist+1,0);

  for(unsigned int i=0; i<genomic_sites.size(); i++){
    string cur_chrom = genomic_sites[i].chrom;
    int cur_pos = genomic_sites[i].pos;
    
    bool _continue = true;
    int cur_i = i+1;
    while(_continue){
      if(cur_i < (int)genomic_sites.size()){
	if(cur_chrom == genomic_sites[cur_i].chrom){
	  int dist = genomic_sites[cur_i].pos - cur_pos;
	  if(dist <= max_dist && dist>=0 ) phase[dist]++;
	  else _continue = false;
	} 
	else _continue = false;
      } 
      else _continue = false;
      cur_i++;
    }
  }

  for(unsigned int i=0; i<phase.size(); i++){
    ofstr<<i<<"\t"<<phase[i]<<endl;
  }
  

  ofstr.close();

  return 0;
}
