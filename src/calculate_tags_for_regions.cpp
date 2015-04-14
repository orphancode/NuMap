#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <exception>
#include <assert.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "params.h"
#include "string_utils.h"
#include "genome_table.h"
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

struct genomic_region{
  string chrom;
  int begin;
  int end;
  
  int helper1;
  int helper2;
};

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("analysis_path", "analysis_path", STRING_TYPE);
  pars.require("regions_file", "regions_file", STRING_TYPE);
  pars.require("output_file", "output_file", STRING_TYPE);

  if(!pars.enforce()){
    exit(1);
  }

  string analysis_path = pars.get_string_value("analysis_path");
  string regions_fname = pars.get_string_value("regions_file");
  
  if(analysis_path[analysis_path.size()-1] != '/') analysis_path = analysis_path + "/";
  string genome_table_fname = analysis_path + genome_table_suffix;
  string bin_count_file_prefix = analysis_path + align_bin_suffix;

  string output_fname = pars.get_string_value("output_file");
  ofstream ofstr(output_fname.c_str());
  assert(ofstr.good());

  genome_table gt(genome_table_fname);

  vector<genomic_region> regions;
  ifstream regions_ifstr(regions_fname.c_str());
  assert(regions_ifstr.good());

  string line;

  while(regions_ifstr.good()){
    getline(regions_ifstr, line);
    if(regions_ifstr.good()){
      vector<string> line_fields = split(line, '\t');
      if(line_fields.size() >= 3){
	genomic_region new_region;
	new_region.chrom = line_fields[0];
	new_region.begin = atoi(line_fields[1].c_str());
	new_region.end = atoi(line_fields[2].c_str());
	regions.push_back(new_region);
      }
    }
  }
  cout<<"Read "<<regions.size()<<" regions"<<endl;
  
  for(unsigned int k=0; k<gt.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];

    string bin_count_fname = bin_count_file_prefix + "." + cur_chrom;
    
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
    
    bin_count_ifstr.read((char*)(&start_counts_for[0]),sizeof(start_counts_for[0])*chr_size);
    bin_count_ifstr.read((char*)(&start_counts_rev[0]),sizeof(start_counts_rev[0])*chr_size);
    
    bin_count_ifstr.close();
    
    cout<<"read: "<<chr_size<<" entries"<<endl;

    for(unsigned int j=0; j<regions.size(); j++){
      if(regions[j].chrom == cur_chrom){
	int count_for=0; int count_rev=0;
	for(int i=regions[j].begin; i<regions[j].end; i++){
	  count_for += start_counts_for[i];
	  count_rev += start_counts_rev[i];
	}
	regions[j].helper1 = count_for;
	regions[j].helper2 = count_rev;
      }
    }

    cout<<endl;
  }
   
  for(unsigned int i=0; i<regions.size(); i++){
    ofstr<<regions[i].chrom<<"\t"<<regions[i].begin<<"\t"<<regions[i].end<<"\t";
    ofstr<<regions[i].helper1<<"\t"<<regions[i].helper2<<endl;
  }
  ofstr.close();

  return 0;
}
