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
//#include "genome.h"
#include "string_utils.h"
#include "genome_table.h"
#include "math_functions.h"
//#include "sequence_utility_functions.h"
#include "kernels.h"
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
  //bool using_value;
  //int helper;
  string record;
};

bool pos_less(const genome_coord& l, const genome_coord& r){
  return l.pos<r.pos;
}

//bool value_less(const genome_coord& l, const genome_coord& r){
//  assert(l.using_value && r.using_value);
//  if(l.value < r.value) return true;
//  else return false;
//}

//compare positions file 1 to positions file 2
//report pos1 overlap pos2
//report pos1 - pos2 

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("dyad_match_file", "dyad_match_file", STRING_TYPE);
  pars.require("unmatched_file1", "unmatched_file1", STRING_TYPE);
  pars.require("unmatched_file2", "unmatched_file2", STRING_TYPE);

  pars.require("positions_file" , "positions_file2", STRING_TYPE);
  pars.require("genome_table", "genome_table", STRING_TYPE);
  pars.require("output_file", "output_file", STRING_TYPE);
  //pars.optional("match_dist", "match_dist", "50", INT_TYPE);
  pars.optional("dist", "dist", "1500", INT_TYPE);


  if(!pars.enforce()){
    exit(1);
  }

  string positions_fname1 = pars.get_string_value("dyad_match_file");
  string positions_fname2 = pars.get_string_value("positions_file");

  string unm_fname1 = pars.get_string_value("unmatched_file1");
  string unm_fname2 = pars.get_string_value("unmatched_file2");


  int dist = pars.get_int_value("dist");

  string gt_fname = pars.get_string_value("genome_table");
  genome_table gt(gt_fname);
  
  string output_fname = pars.get_string_value("output_file");

  ofstream ofstr(output_fname.c_str());
  assert(ofstr.good());

  string detailed_output_fname = output_fname + ".detailed";
  ofstream d_ofstr(detailed_output_fname.c_str());
  assert(d_ofstr.good());
  
  
  ifstream positions_ifstr1(positions_fname1.c_str());
  assert(positions_ifstr1.good());
 
  vector<int> m_counts(dist*2+1, 0);
  vector<int> unm_counts(dist*2+1, 0);
 
  string line;
  vector<genome_coord> sites1;
  char tab_char = '\t';

  while(positions_ifstr1.good()){
    getline(positions_ifstr1, line);
    if(positions_ifstr1.good()){
      vector<string> line_fields = split(line, tab_char);
      if(line_fields.size() >= 2){
	genome_coord cur_site;
	cur_site.chrom = line_fields[0];
	cur_site.pos = atoi(line_fields[1].c_str());
	cur_site.record = line;
	//if(line_fields.size()>=3){
	//  cur_site.value = atof(line_fields[2].c_str());
	//  cur_site.using_value = true;
	//} else {
	//  cur_site.using_value = false;
	//}
	sites1.push_back(cur_site);
      }
    }
  }
  positions_ifstr1.close();

  cout<<"Read "<<sites1.size()<<" matches"<<endl;

  ifstream positions_ifstr2(positions_fname2.c_str());
  assert(positions_ifstr2.good());
  
  vector<genome_coord> sites2;

  while(positions_ifstr2.good()){
    getline(positions_ifstr2, line);
    if(positions_ifstr2.good()){
      vector<string> line_fields = split(line, tab_char);
      if(line_fields.size() >= 2){
	genome_coord cur_site;
	cur_site.chrom = line_fields[0];
	cur_site.pos = atoi(line_fields[1].c_str());
	cur_site.record = line;
	//if(line_fields.size()>=3){
	//  cur_site.value = atof(line_fields[2].c_str());
	//  cur_site.using_value = true;
	//} else {
	//  cur_site.using_value = false;
	//}
	sites2.push_back(cur_site);
	
      }
    }
  }
  positions_ifstr2.close();
  cout<<"Read "<<sites2.size()<<" genomic sites"<<endl;

  ifstream unm_ifstr1(unm_fname1.c_str());
  assert(unm_ifstr1.good());
  vector<genome_coord> unm1;
  while(unm_ifstr1.good()){
    if(unm_ifstr1.good()){
      getline(unm_ifstr1,line); 
      if(unm_ifstr1.good()){
	vector<string> fields = split(line, tab_char);
	genome_coord new_coord;
	new_coord.chrom = fields[0]; new_coord.pos = atoi(fields[1].c_str()); new_coord.value = atof(fields[2].c_str());
	unm1.push_back(new_coord);
      }
    }
  }
  unm_ifstr1.close();
  cout<<"Read "<<unm1.size()<<" unmatched dyads from file 1"<<endl;
  
  ifstream unm_ifstr2(unm_fname2.c_str());
  assert(unm_ifstr2.good());
  vector<genome_coord> unm2;
  while(unm_ifstr2.good()){
    if(unm_ifstr2.good()){
      getline(unm_ifstr2,line); 
      if(unm_ifstr2.good()){
	vector<string> fields = split(line, tab_char);
	if(fields.size()>=3){
	  genome_coord new_coord;
	  new_coord.chrom = fields[0]; new_coord.pos = atoi(fields[1].c_str()); new_coord.value = atof(fields[2].c_str());
	  unm2.push_back(new_coord);
	}
      }
    }
  }
  unm_ifstr2.close();
  cout<<"Read "<<unm2.size()<<" unmatched dyads from file 2"<<endl;
  
  cout<<"Read "<<sites2.size()<<" genomic positions"<<endl;

  for(unsigned int k=0; k<gt.contigs.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];
    cout<<"Processing "<<cur_chrom<<endl;

    vector<int> match_inds(cur_chrom_size, -1);
    for(unsigned int i=0; i<sites1.size(); i++){
      if(sites1[i].chrom == cur_chrom){
	match_inds[sites1[i].pos] = i;
      }
    }

    vector<int> unm_mask1(cur_chrom_size, -1);
    for(unsigned int i=0; i<unm1.size(); i++){
      if(unm1[i].chrom == cur_chrom && unm1[i].value>=0.4)
	unm_mask1[unm1[i].pos] = i;
    }

    vector<int> unm_mask2(cur_chrom_size, -1);
    for(unsigned int i=0; i<unm2.size(); i++){
      if(unm2[i].chrom == cur_chrom && unm2[i].value>=0.4)
	unm_mask2[unm2[i].pos] = i;
    }

    for(unsigned int i=0; i<sites2.size(); i++){
      if(sites2[i].chrom == cur_chrom){
	int pos = sites2[i].pos;
	int cur_matches = 0;
	int cur_unm1 = 0; int cur_unm2 = 0;
	
	for(int j=pos-dist; j<=pos+dist; j++){
	  int offset = j-pos;
	  if(match_inds[j] != -1){
	    cur_matches++;
	    m_counts[dist+offset]++;
	    //ofstr<<sites1[match_inds[j]].record<<endl;
	  }
	  if(unm_mask1[j] != -1){ cur_unm1++; unm_counts[dist+offset]++; }
	  if(unm_mask2[j] != -1){ cur_unm2++; unm_counts[dist+offset]++; }

	  //if(unm_mask1[j] != -1) ofstr<<unm1[unm_mask1[j]].value<<endl;
	  //if(unm_mask2[j] != -1) ofstr<<unm2[unm_mask2[j]].value<<endl;
	}
	
	ofstr<<sites2[i].chrom<<"\t"<<sites2[i].pos;
	ofstr<<"\t"<<cur_matches<<"\t"<<cur_unm1<<"\t"<<cur_unm2<<endl;
      }
    }

  }
  ofstr.close();

  for(int i=0; i<(int)m_counts.size(); i++){
    d_ofstr<<i-dist<<"\t"<<m_counts[i]<<"\t"<<unm_counts[i]<<endl;
  }
  d_ofstr.close();

  return 0;
}
