#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <exception>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "paths.h"
#include "params.h"
#include "string_utils.h"
#include "genome_table.h"
#include "math_functions.h"

using std::ifstream;
using std::ofstream;

using std::istringstream;
using std::cerr;
using std::endl;
using std::cout;


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

bool mkdir(string& path){
  int status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  return status; //true;
}

bool directory_exists(string& path){
  struct stat sb;
  if(stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) return true;
  return false;
}

bool create_file_path(string& file_path){
  char delim = '/';
  vector<string> path_fields = split(file_path, delim);
  if(path_fields.size() == 1){
    cout<<"Warning, only one field in the path, no action taken"<<endl;
    return true;
  }
  string cur_path = "";
  for(unsigned int i=0; i<path_fields.size()-1; i++){
    cur_path = cur_path + path_fields[i] + "/";
    if(!directory_exists(cur_path)){
      cout<<"Creating directory "<<cur_path<<endl;
      mkdir(cur_path);
      if(!directory_exists(cur_path)){
	cout<<"Failed to create directory "<<cur_path<<endl;
	exit(1);
      }
    }
  }
  return true;
}

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("analysis_path", "analysis_path", STRING_TYPE);
  pars.optional("mock", "yes/no", "no", STRING_TYPE);
  //pars.require("bin_count_file_prefix","bin_count_file_prefix",STRING_TYPE);
  //pars.require("output_prefix","output_prefix", STRING_TYPE);  
  //pars.require("genome_table","genome_table", STRING_TYPE);
  pars.optional("core_size","core_size","147",INT_TYPE);
  //pars.optional("fragment_size", "fragment_size", "153", INT_TYPE);
  pars.optional("bw", "bw", "100", INT_TYPE); //KDE bandwidth

  if(!pars.enforce()){
    exit(1);
  }

  string mock = pars.get_string_value("mock");
  assert(mock == "yes" || mock == "no"); 

  string analysis_path = pars.get_string_value("analysis_path");
  if(analysis_path[analysis_path.size()-1] != '/'){
    analysis_path = analysis_path + "/";
  }
  
  string gt_fname = analysis_path + genome_table_suffix;
  string dyad_bin_file_prefix;
  if(mock == "no")
    dyad_bin_file_prefix = analysis_path + dyad_bin_suffix;
  else 
    dyad_bin_file_prefix = analysis_path + mock_dyad_bin_suffix;
  string dyad_stringency_file_prefix;
  if(mock == "no")
    dyad_stringency_file_prefix = analysis_path + dyad_stringency_suffix;
  else
    dyad_stringency_file_prefix = analysis_path + mock_dyad_stringency_suffix;
  
  cout<<"Creating the directory path "<<dyad_stringency_file_prefix<<endl;
  create_file_path(dyad_stringency_file_prefix);
  
  //int core_size = pars.get_int_value("core_size");
  //int half_core = (core_size - core_size%2) / 2;

  //int fragment_size = pars.get_int_value("fragment_size");
  int bw = pars.get_int_value("bw");

  //int dyad_shift = (fragment_size - fragment_size%2)/2;
  //int half_core = (core_size -  core_size%2)/2;

  cout<<"initializing kernel"<<endl;
  precomputed_kernel kernel(bw, "triweight");


  genome_table gt(gt_fname);
  
  for(unsigned int k=0; k<gt.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];
    string dyad_fname = dyad_bin_file_prefix + "." + cur_chrom;
    
    ifstream dyad_ifstr(dyad_fname.c_str());
    if(!dyad_ifstr.good()){
      cerr<<"Bad file: "<<dyad_fname<<endl;
      exit(1);
    }
    
    cout<<endl<<"processing "<<dyad_fname<<endl;
    string output_fname = dyad_stringency_file_prefix + "." + cur_chrom; 
    ofstream stringency_ofstr(output_fname.c_str());
    
    int chr_size;
    dyad_ifstr.read((char*)(&chr_size), sizeof(chr_size));
    assert(chr_size == cur_chrom_size); 
    
    //vector<int> start_counts_for(gt.contig_sizes[k]);
    //vector<int> start_counts_rev(gt.contig_sizes[k]);
    
    vector<int> dyads(cur_chrom_size, 0);
    vector<double> dyads_kde(cur_chrom_size, 0);
    //vector<int> core_coverage(gt.contig_sizes[k], 0);

    dyad_ifstr.read((char*)(&dyads[0]), sizeof(int)*cur_chrom_size); 
    dyad_ifstr.close();
    
    /*
    cout<<"calculating core coverage"<<endl;
    for(int i=0; i<(int) dyads.size(); i++){
      if(dyads[i] > 0){
	for(int j=i-half_core; j<=i+half_core; j++){
	  if(j>=0 && j<(int)core_coverage.size()){
	    core_coverage[j] += dyads[i];
	  }
	}
      }
    }
    */
    cout<<"smoothing dyads"<<endl;
    //calculate positioning stringency
    
    for(int i=0; i<(int) dyads.size(); i++){
      if(dyads[i] > 0){
	for(int j=i-kernel.limit; j<=i+kernel.limit; j++){
	  if(j>=0 && j<(int) dyads.size()){
	    dyads_kde[j] += dyads[i] * kernel.at(i-j); //triweight_kernel((double)(i-j), bw);
	  }
	}
      }
    }
    dyads.resize(0);
    
    cout<<"calculating infringement score"<<endl;
    vector<double> infringing_dyads(cur_chrom_size, 0);
    vector<double> stringency(cur_chrom_size, 0);
    
    int infringing_boundary = 150;
    for(int i=0; i<=infringing_boundary; i++){
      for(int j=i-infringing_boundary; j<=i+infringing_boundary; j++){
	if(j>=0 && j<(int) dyads_kde.size()){
	  infringing_dyads[i] += dyads_kde[j];
	}
      }
    }
    for(int i=(int) dyads_kde.size() - infringing_boundary; i<(int) dyads_kde.size(); i++){
      for(int j=i-infringing_boundary; j<=i+infringing_boundary; j++){
	if(j>=0 && j<(int) dyads_kde.size()){
	  infringing_dyads[i] += dyads_kde[j];
	}	
      }
      if(infringing_dyads[i] < 0){
	infringing_dyads[i] = 0;
      }
    }

    for(unsigned int i=infringing_boundary+1; i<dyads_kde.size()-infringing_boundary; i++){
      infringing_dyads[i] = infringing_dyads[i-1] + 
	dyads_kde[i+infringing_boundary] - dyads_kde[i-infringing_boundary-1];
    }

    for(unsigned int i=0; i<dyads_kde.size(); i++){
      if(infringing_dyads[i] <= 0){
	stringency[i] = 0;
      } else{
	//stringency[i] = dyads_kde[i] * bw / infringing_dyads[i];
	stringency[i] = dyads_kde[i] * bw / (1.09375 * infringing_dyads[i]);

	if(i%100000 == 0){
	  //printf("\r%.2f  ", stringency[i]);
	  cout.flush();
	}
      }
    }
    cout<<endl;
    
    cout<<"writing the profile to disk"<<endl;
    //int chr_size = gt.contig_sizes[k];
    stringency_ofstr.write((char*)(&chr_size), sizeof(chr_size));
    stringency_ofstr.write((char*)(&stringency[0]), sizeof(stringency[0])*chr_size);
    //stringency_ofstr.write((char*)(&core_coverage[0]), sizeof(core_coverage[0])*chr_size);
    stringency_ofstr.close();
  }

  /*
  //vector<double> dyad_stringency(dyad_hist_for.size(), 0);
  vector<double> infringing_dyads(dyad_hist_for.size(), 0);
  
  for(int i=0; i<(int) dyad_hist_kde.size(); i++){
    for(int j=i-150; j<=i+150; j++){
      if(j>=0 && j<=(int) dyad_hist_for.size()){
	infringing_dyads[i] += dyad_hist_kde[j];
      }
    }
  }

  vector<double> stringency(dyad_hist_for.size(), 0);
  for(int i=0; i<(int) dyad_hist_kde.size(); i++){
    stringency[i] = dyad_hist_kde[i] * bw / infringing_dyads[i];
  }

  double stringency_ave = mean(stringency);
  vector<double> stringency_normalized(dyad_hist_for.size(), 0);
  for(int i=0; i<(int) dyad_hist_kde.size(); i++){
    stringency_normalized[i] = stringency[i] / stringency_ave;
  }


  dyad_density_ofstr<<"#dist density stringency stringency_norm"<<endl;
  for(int i=0; i<(int) dyad_hist_kde.size(); i++){
    dyad_density_ofstr<<i-max_dist<<" "<<" "<<dyad_hist_kde[i]<<" "<<stringency[i]<<" "<<stringency_normalized[i]<<endl;
  }

  dyad_density_ofstr.close();
  
  for(int i=0; i<(int)nucl_coverage_for.size(); i++){
    nucl_coverage_ofstr<<i-max_dist<<" "<<nucl_coverage_for[i]<<" "<<nucl_coverage_rev[i]<<" "<<nucl_coverage_for[i]+nucl_coverage_rev[i]<<endl;      
  }
  nucl_coverage_ofstr.close();

  for(int i=0; i<(int)dyad_hist_for.size(); i++){
    dyad_hist_ofstr<<i-max_dist<<" "<<dyad_hist_for[i]<<" "<<dyad_hist_rev[i]<<" "<<dyad_hist_for[i] + dyad_hist_rev[i]<<endl;
  }
  dyad_hist_ofstr.close();
  */
  return 0;
}
