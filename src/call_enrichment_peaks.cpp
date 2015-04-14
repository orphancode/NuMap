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

struct region{
  string chrom;
  int begin;
  int end;
};

class precomputed_kernel{
public:
  int bw;
  string type;
  int zero_ind; //index of zero in the values profile                                             
  int limit; //kernel support [-limit, limit]                                                     

  vector<double> values;

  precomputed_kernel(double _bw, string _type);
  double at(int ind); //returns precomputed value at ind                                          
  double& operator[](int i);

  double sum(); //sanity check to make sure kernel adds up to 1                                   
};

precomputed_kernel::precomputed_kernel(double _bw, string _type){
  bw = _bw;
  type = _type;
  assert(type == "triweight");

  if(type == "triweight"){
    limit = bw;
    values.resize(2*bw + 1, 0);
    zero_ind = bw;

    for(int i=-bw; i<=bw; i++){
      (*this)[i] = triweight_kernel(i, bw);
    }
  }
}

double& precomputed_kernel::operator[](int i){
  int offset = i + zero_ind;
  assert(offset >= 0 && offset < (int) values.size());
  return values[offset];
}
double precomputed_kernel::at(int ind){
  int offset = ind + zero_ind;

  if(offset < 0 || offset >= (int) values.size()){
    return 0;
  } else {
    return values[offset];
  }
}
double precomputed_kernel::sum(){
  double res = 0;
  for(unsigned int i=0; i<values.size(); i++){
    res += values[i];
  }
  return res;
}

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("analysis_path", "analysis_path", STRING_TYPE);
  pars.require("output_path", "output_path", STRING_TYPE);
  //pars.optional("mock", "yes/no", "no", STRING_TYPE);
  pars.optional("peak_threshold", "peak_threshold", "3.0", DOUBLE_TYPE);
  pars.optional("drop", "drop", "0.1", DOUBLE_TYPE); //the value of the profile should drop by this of the main peak
  pars.optional("drop_offset", "drop_offset", "50", INT_TYPE);
  pars.optional("max_w", "maximum_window", "100", INT_TYPE);
  pars.optional("invert_profile", "yes/no", "no", STRING_TYPE);

  if(!pars.enforce()){
    exit(1);
  }  

  double peak_threshold = pars.get_double_value("peak_threshold");
  double drop = pars.get_double_value("drop");
  int max_w = pars.get_int_value("max_w");
  string invert_profile_string = pars.get_string_value("invert_profile");

  int max_hw = max_w/2;
  
  bool invert_profile;
  if(invert_profile_string == "yes"){ invert_profile=true; }
  else{ if(invert_profile_string == "no"){ invert_profile = false; } else {assert(false);}}
    

  //string track_name = pars.get_string_value("track_name");

  int drop_offset = pars.get_int_value("drop_offset");
  assert(drop_offset > 0 && drop_offset<1000);

  string analysis_path = pars.get_string_value("analysis_path");
  string output_path = pars.get_string_value("output_path");

  if(output_path[output_path.size()-1] != '/')
    output_path = output_path + "/";  

  if(analysis_path[analysis_path.size()-1]!='/')
    analysis_path = analysis_path + "/";


  string gt_fname = analysis_path + genome_table_suffix;
  genome_table gt(gt_fname);

  string dyad_score_prefix = analysis_path + dyad_score_suffix;

  string peak_positions_fname = output_path + "peaks.txt";
  create_file_path(peak_positions_fname);

  ofstream peak_ofstr(peak_positions_fname.c_str());
  assert(peak_ofstr.good());
  peak_ofstr<<"#chrom\tpos\tenrchiment"<<endl;

  string peak_bed_fname = output_path + "DNase-FLASH_TF.peaks.bed";
  ofstream peak_bed_ofstr(peak_bed_fname.c_str());
  assert(peak_bed_ofstr.good());
  string track_name = "DNase_FLASH_TF_peaks";
  peak_bed_ofstr<<"track name=\"";
  peak_bed_ofstr<<track_name;
  peak_bed_ofstr<<"\" description=\""<<track_name<<"\"visibility=2 colorByStrand=\"69,139,0 69,139,0\""<<endl;

  int total_peaks = 0;
  
  for(unsigned int k=0; k<gt.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];
  
    string score_fname = dyad_score_prefix + "." + cur_chrom;

    cout<<"processing "<<score_fname<<endl;
    ifstream score_ifstr(score_fname.c_str());
    assert(score_ifstr.good());

    int chr_size;

    cout<<"Reading scores... "; cout.flush();
    score_ifstr.read((char*)(&chr_size), sizeof(chr_size));
    vector<double> scores(chr_size);
    
    score_ifstr.read((char*)(&scores[0]), sizeof(double)*chr_size);
    score_ifstr.close();
    cout<<" done!"<<endl;

    if(invert_profile){
      cout<<"inverting_profile"<<endl;
      for(unsigned int i=0; i<scores.size(); i++){
	scores[i] = -scores[i];
      }
    }

    cout<<"calling peaks"<<endl;
    int peaks = 0;
    for(int i=1000; i<cur_chrom_size-1000; i++){
      if(scores[i] > scores[i-1] && scores[i] > scores[i+1] &&
	 scores[i] >= peak_threshold && 
	 scores[i-drop_offset] <= scores[i]*(1.0-drop) &&
	 scores[i+drop_offset] <= scores[i]*(1.0-drop) ){
	
	bool peak_accepted = true;
	for(int j=i-max_hw; j<=i+max_hw; j++){
	  if(j!=i && scores[j] > scores[i]) peak_accepted = false;
	}
	if(peak_accepted){
	  peaks++; total_peaks++;
	  
	  peak_ofstr<<cur_chrom<<"\t"<<i<<"\t"<<scores[i]<<endl;
	  peak_bed_ofstr<<cur_chrom<<"\t"<<i-2<<"\t"<<i+2<<endl;
	}
      }
    }
    cout<<"Peaks: "<<peaks<<endl;

    /*
    char buf[50], buf1[50], buf2[50];
    double count_kde = dyads_kde[i];
    double count_fe = dyads_kde[i] / norm;
    sprintf(buf, "%.4f", stringency_values[i]);
    sprintf(buf1, "%.4f", count_kde);
    sprintf(buf2, "%.4f", count_fe);
    if(count_fe > 0.1){
      nucleosome_bed_ofstr<<cur_chrom<<"\t"<<i-73+1<<"\t"<<i+74+1<<"\t"<<buf<<"\t0\t+"<<endl;
      dyad_ofstr<<cur_chrom<<"\t"<<i<<"\t"<<buf<<"\t"<<buf2<<"\t"<<buf1<<endl;
    }
    nucleosome_bed_ofstr.close();
    */
  }

  cout<<endl;
  cout<<"Found "<<total_peaks<<" peaks"<<endl;
  peak_ofstr.close();
  peak_bed_ofstr.close();
  
  
  return 0;
}
