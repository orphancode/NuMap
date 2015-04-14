#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <exception>
#include <stdio.h>
#include <assert.h>

#include <time.h>

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
};

int main(int argc, char* argv[]){
  srand( time(NULL) );

  params pars(argc, argv);

  pars.require("analysis_path", "analysis_path", STRING_TYPE);
  pars.optional("window_size", "window_size", "1000", INT_TYPE);
  
  if(!pars.enforce()){
    exit(1);
  }

  int window_size = pars.get_int_value("window_size");

  string analysis_path = pars.get_string_value("analysis_path");
  if(analysis_path[analysis_path.size()-1]!='/')
    analysis_path = analysis_path + "/";
  
  string dyad_bin_prefix = analysis_path + dyad_bin_suffix;
  string mock_dyad_bin_prefix = analysis_path + mock_dyad_bin_suffix;
  
  string gt_fname = analysis_path + genome_table_suffix;
  genome_table gt(gt_fname);
  
  for(unsigned k=0; k<gt.contigs.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];
    
    string dyad_bin_fname = dyad_bin_prefix + "." + cur_chrom;
    cout<<"processing "<<dyad_bin_fname<<endl;
    ifstream dyad_bin_ifstr(dyad_bin_fname.c_str());
    assert(dyad_bin_ifstr.good());
    
    int chr_size;
    dyad_bin_ifstr.read((char*)(&chr_size), sizeof(chr_size));
    assert(chr_size == cur_chrom_size);
    vector<int> dyad_counts(chr_size);
    dyad_bin_ifstr.read((char*)(&dyad_counts[0]), sizeof(dyad_counts[0])*chr_size);
    dyad_bin_ifstr.close();
    
    vector<int> mock_dyad_counts(chr_size, 0);
    string mock_dyad_bin_fname = mock_dyad_bin_prefix + "." + cur_chrom;
    ofstream mock_dyad_bin_ofstr(mock_dyad_bin_fname.c_str());
    assert(mock_dyad_bin_ofstr.good());
    
    for(int i=0; i<cur_chrom_size-window_size; i+= window_size){
      int loc_dyads = 0;
      for(int j=i; j<i+window_size; j++){
	loc_dyads += dyad_counts[j];
      }
      for(int j=0; j<loc_dyads; j++){
	int mock_dyad_pos = i + rand()%window_size;
	mock_dyad_counts[mock_dyad_pos]++;
      }
      int mock_loc_dyads = 0;
      for(int j=i; j<i+window_size; j++){
	mock_loc_dyads += mock_dyad_counts[j];
      }
      assert(loc_dyads == mock_loc_dyads); //check the work
    }

    mock_dyad_bin_ofstr.write((char*)(&chr_size), sizeof(chr_size));
    mock_dyad_bin_ofstr.write((char*)(&mock_dyad_counts[0]), chr_size*sizeof(mock_dyad_counts[0]));
    mock_dyad_bin_ofstr.close();
  }

  return 0;
}
