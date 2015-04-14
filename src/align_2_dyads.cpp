#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <exception>
#include <assert.h>
#include <stdio.h>

#include "params.h"
#include "string_utils.h"
#include "genome_table.h"
#include "math_functions.h"
#include "paths.h"

using std::ifstream;
using std::ofstream;

using std::cerr;
using std::endl;
using std::cout;

struct region{
  string chrom;
  int begin;
  int end;
};

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("analysis_path", "analysis_path", STRING_TYPE);

  if(!pars.enforce()){
    exit(1);
  }

  //string regions_fname = pars.get_string_value("QuEST_regions_file");
  string analysis_path = pars.get_string_value("analysis_path");
  
  if(analysis_path[analysis_path.size()-1] != '/') analysis_path = analysis_path + "/";
  string genome_table_fname = analysis_path + genome_table_suffix;
  string bin_count_file_prefix = analysis_path + align_bin_suffix;
  string dyad_bin_file_prefix = analysis_path + dyad_bin_suffix;
  string fragment_estimate_fname = analysis_path + dist_plots_suffix + "whole_genome/" + "fragment_estimate.txt";

  ifstream fragment_estimate_ifstr(fragment_estimate_fname.c_str());
  if(!fragment_estimate_ifstr.good()){
    cout<<"Fragment estimate file "<<fragment_estimate_fname<<" is missing"<<endl;
    cout<<"You need to calculate distance plots first (dist_plots)"<<endl;
    exit(1);
  }

  string cur_line;
  getline(fragment_estimate_ifstr, cur_line);
  char space_char = ' ';
  vector<string> cur_line_fields = split(cur_line, space_char);
  assert(cur_line_fields.size() == 2);
  assert(cur_line_fields[0] == "size_estimate:");
  
  int fragment_size = atoi(cur_line_fields[1].c_str());
  fragment_estimate_ifstr.close();

  cout<<"Using size estimate of "<<fragment_size<<endl;
  int shift_dist = (fragment_size-fragment_size%2)/2+1;

  genome_table gt(genome_table_fname);


  for(unsigned int k=0; k<gt.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];

    //vector<int> dyad_counts(cur_chrom_size, 0);
    
    string bin_count_fname = bin_count_file_prefix + "." + cur_chrom;    
    string dyad_bin_count_fname = dyad_bin_file_prefix + "." + cur_chrom;
    
    ifstream bin_count_ifstr(bin_count_fname.c_str());
    if(!bin_count_ifstr.good()){
      cerr<<"Bad file: "<<bin_count_fname<<endl;
      exit(1);
    }
    
    cout<<endl<<"processing "<<bin_count_fname<<endl;
    
    int chr_size;
    bin_count_ifstr.read((char*)(&chr_size), sizeof(chr_size));
    assert(chr_size == cur_chrom_size); 

    vector<int> start_counts_for(cur_chrom_size);
    vector<int> start_counts_rev(cur_chrom_size);

    vector<int> dyad_counts(cur_chrom_size, 0);
    
    bin_count_ifstr.read((char*)(&start_counts_for[0]),sizeof(start_counts_for[0])*chr_size);
    bin_count_ifstr.read((char*)(&start_counts_rev[0]),sizeof(start_counts_rev[0])*chr_size);
    for(int i=shift_dist; i<cur_chrom_size-shift_dist; i++){
      dyad_counts[i] = start_counts_for[i-shift_dist] + start_counts_rev[i+shift_dist];
    }
    
    bin_count_ifstr.close();
    
    ofstream dyad_bin_count_ofstr(dyad_bin_count_fname.c_str());
    dyad_bin_count_ofstr.write((char*)(&cur_chrom_size), sizeof(int));
    dyad_bin_count_ofstr.write((char*)(&dyad_counts[0]), sizeof(int)*cur_chrom_size);
    dyad_bin_count_ofstr.close();

    cout<<endl;
  }    

  return 0;
}
