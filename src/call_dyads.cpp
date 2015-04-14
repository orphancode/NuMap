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
  pars.optional("mock", "yes/no", "no", STRING_TYPE);
  pars.optional("peak_value", "peak_value", "mean", STRING_TYPE);// DOUBLE_TYPE);
  pars.optional("drop", "drop", "0.1", DOUBLE_TYPE); //the value of the profile should drop by this much to the left and right
  pars.optional("with_trough", "yes/no", "no", STRING_TYPE);
  pars.optional("drop_window", "drop_window", "150", INT_TYPE);
  pars.optional("track_name", "track_name", "nucleosomes", STRING_TYPE);

  pars.optional("bw", "bw", "100", INT_TYPE);

  if(!pars.enforce()){
    exit(1);
  }  

  int bw = pars.get_int_value("bw");
  string peak_value_str = pars.get_string_value("peak_value");
  double peak_value = -1;
  bool using_mean_peak_value = false;
  if(peak_value_str == "mean") using_mean_peak_value = true;
  else {
    using_mean_peak_value = false;
    peak_value = atof(peak_value_str.c_str());
    assert(peak_value > 0 && peak_value <= 1);
  }

  double trough_drop = pars.get_double_value("drop");

  string track_name = pars.get_string_value("track_name");

  string with_trough_string = pars.get_string_value("with_trough");
  assert(with_trough_string == "no" || with_trough_string == "yes");
  bool with_trough = false;
  if(with_trough_string == "no") with_trough = false;
  if(with_trough_string == "yes") with_trough = true;

  int dw = pars.get_int_value("drop_window");
  assert(dw > 0 && dw<1000);

  string analysis_path = pars.get_string_value("analysis_path");
  string output_path = pars.get_string_value("output_path");
  if(output_path[output_path.size()-1] != '/')
    output_path = output_path + "/";
  

  if(analysis_path[analysis_path.size()-1]!='/')
    analysis_path = analysis_path + "/";

  string mock = pars.get_string_value("mock");
  assert(mock == "yes" || mock == "no");

  create_file_path(output_path + "/dummy");

  string gap_fname = output_path + "gaps.txt";
  string gap_bed_fname = output_path + "gaps.bed";
  ofstream gap_ofstr(gap_fname.c_str());
  assert(gap_ofstr.good());

  ofstream gap_bed_ofstr(gap_bed_fname.c_str());
  assert(gap_bed_ofstr.good());

  gap_bed_ofstr<<"track name=\"";
  gap_bed_ofstr<<track_name<<"_coverage_gaps";
  gap_bed_ofstr<<"\" description=\""<<track_name<<"_coverage_gaps"<<"\"visibility=2 colorByStrand=\"255,0,0 0,0,255\""<<endl;


  string gt_fname = analysis_path + genome_table_suffix;
  genome_table gt(gt_fname);

  string dyad_bin_prefix;
  if(mock == "no") dyad_bin_prefix = analysis_path + dyad_bin_suffix;
  else dyad_bin_prefix = analysis_path + mock_dyad_bin_suffix;
  
  string dyad_stringency_prefix;
  if(mock == "no") dyad_stringency_prefix = analysis_path + dyad_stringency_suffix;
  else dyad_stringency_prefix = analysis_path + mock_dyad_stringency_suffix; 

  string dyad_positions = output_path + "dyad_positions.txt";

  ofstream dyad_ofstr(dyad_positions.c_str());
  assert(dyad_ofstr.good());
  dyad_ofstr<<"#chrom\tpos\tstringency\tfold_enrichment"<<endl;
  //count_enrichment\tcount"<<endl;
  
  //string phasogram_fname = output_path + "phasogram";
  //ofstream phasogram_ofstr(phasogram_fname.c_str());
  //assert(phasogram_ofstr.good());

  precomputed_kernel kernel(bw, "triweight");
  //cout<<"kernel(0):"<<bw*kernel.at(-2)/1.09<<endl;

  //int max_dist = 3000;
  //vector<int> peak_phase(max_dist+1,0);
  
  int neighborhood = 75;

  for(unsigned int k=0; k<gt.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];
  
    string dyad_fname = dyad_bin_prefix + "." + cur_chrom;
    string stringency_fname = dyad_stringency_prefix + "." + cur_chrom;
    cout<<"processing "<<dyad_fname<<endl;
    cout<<"processing "<<stringency_fname<<endl;
    ifstream stringency_ifstr(stringency_fname.c_str());
    assert(stringency_ifstr.good());

    string nucleosome_bed_fname = output_path + "nucleosome." + cur_chrom + ".bed";
    ofstream nucleosome_bed_ofstr(nucleosome_bed_fname.c_str());
    assert(nucleosome_bed_ofstr.good());

    nucleosome_bed_ofstr<<"track name=\"";
    nucleosome_bed_ofstr<<track_name<<"_"<<cur_chrom;
    nucleosome_bed_ofstr<<"\" description=\""<<track_name<<"\"visibility=2 colorByStrand=\"255,0,0 0,0,255\""<<endl;

    int chr_size;
    stringency_ifstr.read((char*)(&chr_size), sizeof(chr_size));
    vector<double> stringency_values(chr_size);
    
    stringency_ifstr.read((char*)(&stringency_values[0]), sizeof(double)*chr_size);
    stringency_ifstr.close();
    cout<<"read stringencies"<<endl;

    ifstream dyad_ifstr(dyad_fname.c_str());
    if(!dyad_ifstr.good()){
      cerr<<"Bad file: "<<dyad_fname<<endl;
      exit(1);
    }

    vector<int> dyads(cur_chrom_size, 0);
    vector<double> dyads_kde(cur_chrom_size, 0);

    dyad_ifstr.read((char*)(&chr_size), sizeof(chr_size));
    assert(chr_size == cur_chrom_size);

    dyad_ifstr.read((char*)(&dyads[0]), sizeof(int)*cur_chrom_size);
    dyad_ifstr.close();
    cout<<"read dyads"<<endl;

    //calculate dyad kde
    cout<<"smoothing dyads"<<endl;
    for(int i=0; i<(int) dyads.size(); i++){
      if(i%100000 == 0){
	if(i>1000) cout<<"\r"<<((double)i)/1000000.0<<" M   ";//<<dyads_kde[i-1000]<<"       "; //<<kernel.at(1000)<<"   ";
	cout.flush();
      }
      if(dyads[i] > 0){
        for(int j=i-kernel.limit; j<=i+kernel.limit; j++){
          if(j>=0 && j<(int) dyads.size()){
            dyads_kde[j] += dyads[i] * kernel.at(i-j)*bw/1.09375; //triweight_kernel((double)(i-j), bw);     
          }
        }
      }
    }
    cout<<endl;
    //dyads.resize(0);

    /*
    double sum = 0;
    for(int i=0; i<(int)dyads_kde.size(); i++){
      sum += dyads_kde[i];
    }
    double norm = sum / dyads_kde.size();
    */
    //cout<<"Norm: "<<norm<<endl;

    //find gaps of coverage
    vector<region> coverage_gaps;
    bool in_a_gap = false;
    int gap_begin = -1;
    for(int i=0; i<(int)dyads_kde.size(); i++){
      if(in_a_gap){
	if(dyads_kde[i] <= 1.1){
	  //do nothing, wait till the end of the gap
	} else {
	  if(i-gap_begin >= 300){
	    region new_region;
	    new_region.chrom = cur_chrom;
	    new_region.begin = gap_begin;
	    new_region.end = i+1;
	    coverage_gaps.push_back(new_region);
	  }
	  in_a_gap = false;
	}
      } else {
	if(dyads_kde[i] < 1.1){
	  in_a_gap = true;
	  gap_begin = i;
	}
      }
    }
    //process last gap
    if(in_a_gap){
      region new_region;
      new_region.chrom = cur_chrom;
      new_region.begin = gap_begin;
      new_region.end = dyads_kde.size();
      coverage_gaps.push_back(new_region);
      in_a_gap = false;
    }
    cout<<"Found "<<coverage_gaps.size()<<" gaps"<<endl;

    vector<bool> coverage_mask(cur_chrom_size, true);
    for(unsigned int i=0; i<coverage_gaps.size(); i++){
      gap_ofstr<<coverage_gaps[i].chrom<<"\t"<<coverage_gaps[i].begin<<"\t"<<coverage_gaps[i].end<<endl;
      gap_bed_ofstr<<cur_chrom<<"\t"<<coverage_gaps[i].begin+1<<"\t"<<coverage_gaps[i].end+1<<"\t0\t0\t+"<<endl;
      for(int j=coverage_gaps[i].begin; j<coverage_gaps[i].end; j++){
	coverage_mask[j] = false;
      }
    }

    double sum = 0;
    int denom = 0;
    for(int i=0; i<(int)dyads_kde.size(); i++){
      if(coverage_mask[i]){
	sum += dyads_kde[i];
	denom ++;
      }
    }
    double norm = sum / ((double) denom); //dyads_kde.size();
    cout<<"Norm: "<<norm<<endl;


    vector<bool> peak_mask(stringency_values.size(), false);
    vector<bool> trough_mask(stringency_values.size(), false);
    int troughs = 0;

    if(using_mean_peak_value){
      double str_sum = 0;
      for(int i=0; i<(int)stringency_values.size(); i++){
	if(coverage_mask[i]){
	  str_sum += stringency_values[i];
	}
      }
      peak_value = str_sum / denom;
      cout<<"Using peak value of "<<peak_value<<endl;
    }


    cout<<"calling dyad peaks"<<endl;

    //find all troughs
    for(unsigned int i=1000; i<stringency_values.size()-1000; i++){
      if(stringency_values[i] < stringency_values[i-1] &&
	 stringency_values[i] < stringency_values[i+1]){
	//local max
	//check if it's the maximum in the bigger window
	double min_value = stringency_values[i];
	bool is_min = true;
	for(unsigned int j=i-neighborhood; j<=i+neighborhood; j++){
	  if(stringency_values[j] < min_value){
	    is_min = false;
	  }
	}
	if(is_min) trough_mask[i] = true;
	troughs++;
      }
    }
    cout<<"found "<<commify(troughs)<<" troughs"<<endl;
    
    //search for peaks
    int peaks = 0;
    for(unsigned int i=1000; i<stringency_values.size()-1000; i++){
      if(stringency_values[i] > stringency_values[i-1] &&
	 stringency_values[i] > stringency_values[i+1]){
	//this is local maximum. No check if it's the maximum in the 
	//bigger neighborhood
	
	bool is_max = true;
	double max_value = stringency_values[i];
	for(unsigned int j=i-neighborhood; j<=i+neighborhood; j++){
	  if(stringency_values[j] > max_value){
	    is_max = false;
	  }
	}
	if(is_max){
	  //find the trough to the left and to the right
	  if(with_trough){
	    int left_troughs = 0;
	    int left_trough_pos = -1;
	    for(unsigned int j= i-75-75; j<=i; j++){
	      if(trough_mask[j] == true){
		left_troughs++;
		left_trough_pos = j;
	      }
	    }
	    int right_troughs = 0;
	    int right_trough_pos = -1;
	    for(unsigned int j=i; j<=i+150; j++){
	      if(trough_mask[j] == true){
		right_troughs++;
		right_trough_pos = j;
	      }
	    }
	    //left_trough_pos =  i -75;
	    //right_trough_pos = i+75;
	    
	    if(left_troughs == 1 && right_troughs == 1 &&
	       max_value - stringency_values[left_trough_pos] >= trough_drop &&
	       max_value - stringency_values[right_trough_pos] >= trough_drop &&
	       max_value >= peak_value){
	      peak_mask[i] = true;
	      peaks++;
	    }
	  } else {
	    double left_min = stringency_values[i];
	    for(unsigned int j=i-dw; j<i; j++){
	      if(stringency_values[j] < left_min) left_min = stringency_values[j];
	    }
	    double right_min = stringency_values[i];
	    for(unsigned int j=i; j<i+dw; j++){
	      if(stringency_values[j] < right_min) right_min = stringency_values[j];
	    }
	    if(max_value >= peak_value && 
	       max_value - left_min >= trough_drop &&
	       max_value - right_min >= trough_drop){
	      peak_mask[i] = true;
	      peaks++;
	    }
	  }
	}
      }
    }
    cout<<"found "<<commify(peaks)<<" peaks"<<endl;
    /*
    for(unsigned int i = 10000; i<peak_mask.size()-10000; i++){
      if(peak_mask[i] == true){
	for(unsigned int j=i+1; j<=i+max_dist; j++){
	  if(peak_mask[j] == true){
	    peak_phase.at(j-i)++;
	  }
	}
      }
    }
    */
    for(unsigned int i=1000; i<peak_mask.size()-1000; i++){
      if(peak_mask[i] == true){
	char buf[50], buf1[50], buf2[50];
	//double count_kde = dyads_kde[i];
	double count_fe = dyads_kde[i] / norm;
	sprintf(buf, "%.4f", stringency_values[i]);
	//sprintf(buf1, "%.4f", count_kde);
	sprintf(buf2, "%.4f", count_fe);
	if(count_fe > 0.1){
	  nucleosome_bed_ofstr<<cur_chrom<<"\t"<<i-73+1<<"\t"<<i+74+1<<"\t"<<buf<<"\t0\t+"<<endl;
	  dyad_ofstr<<cur_chrom<<"\t"<<i<<"\t"<<buf<<"\t"<<buf2<<"\t"<<buf1<<endl;
	}
      }
    }
    nucleosome_bed_ofstr.close();
  }

  /*
  for(unsigned int j=0; j<peak_phase.size(); j++){
    phasogram_ofstr<<j<<"\t"<<peak_phase[j]<<endl;
  }
  */
  
  //phasogram_ofstr.close();
  dyad_ofstr.close();
  gap_ofstr.close();
  gap_bed_ofstr.close();
  return 0;
}
