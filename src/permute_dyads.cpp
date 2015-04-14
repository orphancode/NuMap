#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <exception>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <algorithm>

#include "params.h"
#include "string_utils.h"
#include "genome_table.h"
#include "math_functions.h"
#include "paths.h"

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

bool pos_less(const genome_coord& l, const genome_coord& r){
  return l.pos < r.pos;
}

int assign_matches(int dist, vector<bool>& dyad_mask1, 
		   vector<bool>& dyad_mask2, vector<int>& match_mask){ 
  int res = 0;
  int cur_chrom_size = (int) dyad_mask1.size();
  for(int i=0; i<(int)dyad_mask1.size(); i++){
    //if(dyad_mask1[i] == true && dyad_mask2[i] == true)
    //  cout<<"Bingo!"<<endl;
    if(dyad_mask1[i] == true){
      //cout<<"here"<<endl;
      int closest_peak = -1;
      bool _continue = true;
      int cur_offset = 0;
      while(_continue){
	if(i-cur_offset >= 0 && i-cur_offset < cur_chrom_size){
	  if(dyad_mask2[i - cur_offset] == true){
	    closest_peak = i-cur_offset; _continue = false;
	  }
	}
	if(i+cur_offset >= 0 && i+cur_offset < cur_chrom_size){
	  if(dyad_mask2[i + cur_offset] == true){
	    closest_peak = i+cur_offset; _continue = false;
	  }
	}
	cur_offset++;
	if(cur_offset > dist) _continue = false;
      }
      if(closest_peak != -1){ match_mask[i] = closest_peak; res++; }
    }
  }
  return res;
}

int main(int argc, char* argv[]){
  srand(time(NULL));

  params pars(argc, argv);

  pars.require("called_dyads_file", "called_dyads_file", STRING_TYPE);
  pars.require("genome_table", "genome_table", STRING_TYPE);
  pars.require("output_file", "outptu_file", STRING_TYPE);
  pars.optional("window", "window", "1000", STRING_TYPE); //permute window
  
  if(!pars.enforce()){
    exit(1);
  }


  string dyads_fname1 = pars.get_string_value("called_dyads_file");

  string output_fname = pars.get_string_value("output_file");
  int w = pars.get_int_value("window");

  string gt_fname = pars.get_string_value("genome_table");//analysis_path + genome_table_suffix;
  genome_table gt(gt_fname);

  ofstream ofstr(output_fname.c_str());
  assert(ofstr.good());

  //read the positions
  vector<genome_coord> dyads1;
  ifstream dyads_ifstr1(dyads_fname1.c_str());
  assert(dyads_ifstr1.good());

  string cur_line;
  char delim = '\t';
  while(dyads_ifstr1.good()){
    getline(dyads_ifstr1, cur_line);
    if(dyads_ifstr1.good()){
      vector<string> cur_line_fields = split(cur_line, delim);
      assert(cur_line_fields.size() >= 2);
      genome_coord new_site;
      new_site.chrom = cur_line_fields[0];
      new_site.pos = atoi(cur_line_fields[1].c_str());
      new_site.value = atof(cur_line_fields[2].c_str());
      dyads1.push_back(new_site);
    }
  }
  dyads_ifstr1.close();
  cout<<"Read "<<dyads1.size()<<" dyads from file 1"<<endl;

  
  for(unsigned int k=0; k<gt.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];
    cout<<"Processing "<<cur_chrom<<endl;
    
    vector<bool> dyad_mask(cur_chrom_size, false);
    vector<int> dyad_inds(cur_chrom_size, -1);
    int peaks = 0;
    for(unsigned int i=0; i<dyads1.size(); i++){
      if(dyads1[i].chrom == cur_chrom){
	dyad_mask.at(dyads1[i].pos) = true;
	dyad_inds.at(dyads1[i].pos) = i;
	peaks++;
      }
    }
    cout<<"found "<<commify(peaks)<<" peak in file"<<endl;

    //permuting

    vector<bool> permuted_dyad_mask(cur_chrom_size, false);
    vector<genome_coord> permuted_dyads;
    for(int i=0; i<cur_chrom_size-w; i+=w){
      //int loc_dyads = 0;
      for(int j=i; j<i+w; j++){

	if(dyad_mask[j] == true){
	  bool _continue = true;
	  while(_continue){
	    int new_pos = i+rand()%w;
	    if(new_pos >= 0 && new_pos < cur_chrom_size){
	      if(permuted_dyad_mask[new_pos] == false){
		permuted_dyad_mask[new_pos] = true;
		_continue = false;
		genome_coord cur_dyad;
		cur_dyad.chrom = cur_chrom;
		cur_dyad.pos = new_pos;
		cur_dyad.value = dyads1.at(dyad_inds[j]).value;
		permuted_dyads.push_back(cur_dyad);
	      }
	    }
	  }
	}
      }
    }
    sort(permuted_dyads.begin(), permuted_dyads.end(), pos_less);
    for(unsigned int i=0; i<permuted_dyads.size(); i++){
      ofstr<<permuted_dyads[i].chrom<<"\t"<<permuted_dyads[i].pos<<"\t"<<permuted_dyads[i].value<<endl;
    }
  }


  ofstr.close();
  return 0;
}
