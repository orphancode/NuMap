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
#include "file_utilities.h"

using std::ifstream;
using std::ofstream;

using std::istringstream;
using std::cerr;
using std::endl;
using std::cout;


/*
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
*/
int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("analysis_path", "analysis_path", STRING_TYPE);
  pars.optional("mock", "yes/no", "no", STRING_TYPE);
  pars.optional("core_size","core_size","147",INT_TYPE);
  pars.optional("spacing", "spacing", "193", INT_TYPE);
  //pars.optional("fragment_size", "fragment_size", "153", INT_TYPE);
  //pars.optional("bw", "bw", "100", INT_TYPE); //KDE bandwidth

  if(!pars.enforce()){
    exit(1);
  }

  string mock = pars.get_string_value("mock");
  assert(mock == "yes" || mock == "no"); 

  string analysis_path = pars.get_string_value("analysis_path");
  if(analysis_path[analysis_path.size()-1] != '/'){
    analysis_path = analysis_path + "/";
  }

  int core = pars.get_int_value("core_size");
  int spacing = pars.get_int_value("spacing");
  
  string gt_fname = analysis_path + genome_table_suffix;
  string dyad_bin_file_prefix;
  if(mock == "no")
    dyad_bin_file_prefix = analysis_path + dyad_bin_suffix;
  else 
    dyad_bin_file_prefix = analysis_path + mock_dyad_bin_suffix;
  string coverage_file_prefix;
  if(mock == "no")
    coverage_file_prefix = analysis_path + coverage_suffix;
  else
    coverage_file_prefix = analysis_path + mock_coverage_suffix;
  
  cout<<"Creating the directory path "<<coverage_file_prefix<<endl;
  create_file_path(coverage_file_prefix);
  
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
    string output_fname = coverage_file_prefix + "." + cur_chrom; 
    ofstream ofstr(output_fname.c_str());
    assert(ofstr.good());
    
    int chr_size;
    dyad_ifstr.read((char*)(&chr_size), sizeof(chr_size));
    assert(chr_size == cur_chrom_size); 
    
    vector<int> dyads(cur_chrom_size, 0);

    dyad_ifstr.read((char*)(&dyads[0]), sizeof(int)*cur_chrom_size); 
    dyad_ifstr.close();

    
    vector<double> coverage(cur_chrom_size, 0);
    int total_dyads = 0;
    for(unsigned int i=0; i<dyads.size(); i++){
      if(dyads[i] > 0) total_dyads+= dyads[i];
    }

    double norm = (((double)total_dyads)*core) / ((double)dyads.size());
    norm = norm / (((double)spacing) / ((double)core));
    
    for(int i=0; i<(int)dyads.size(); i++){
      if(dyads[i] > 0){
	for(int j=i-core/2; j<=i+core/2; j++){
	  if(j>=0 && j<cur_chrom_size){
	    coverage[j] += dyads[i];
	  }
	}
      }
    }
    cout<<"norm: "<<norm<<endl;
    
    for(unsigned int i=0; i<coverage.size(); i++){
      coverage[i] = (coverage[i] / norm); // * ((double)spacing) / ((double)core));
    }

    double total_coverage = 0;
    for(unsigned int i=0; i<coverage.size(); i++){
      total_coverage += coverage[i];
    }
    double mean_coverage = total_coverage / ((double) coverage.size());
    cout<<"mean coverage: "<<mean_coverage<<endl;

    cout<<"writing the profile to disk"<<endl;

    ofstr.write((char*)(&chr_size), sizeof(chr_size));
    ofstr.write((char*)(&coverage[0]), sizeof(coverage[0])*chr_size);
    ofstr.close();
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
