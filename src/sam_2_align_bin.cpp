#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <exception>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "paths.h"
#include "params.h"
#include "string_utils.h"
#include "genome_table.h"
#include "psam_alignment_file.h"

using std::ifstream;
using std::ofstream;

using std::istringstream;
using std::cerr;
using std::endl;
using std::cout;

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
  pars.require("genome_table","genome_table", STRING_TYPE);
  pars.require("input_sam_file", "input_sam_file", STRING_TYPE);
  pars.require("analysis_path", "analysis_path", STRING_TYPE);

  if(!pars.enforce()){
    exit(1);
  }

  string input_sam_fname = pars.get_string_value("input_sam_file");
  string gt_fname = pars.get_string_value("genome_table");
  string analysis_path = pars.get_string_value("analysis_path");



  if(analysis_path[analysis_path.size()-1] != '/')
    analysis_path = analysis_path + '/';
  string output_prefix = analysis_path + align_bin_suffix;
  create_file_path(output_prefix);
  //string output_prefix = pars.get_string_value("output_prefix");
  //string same_str_fname = output_prefix + ".same_str";
  //string opp_str_fname = output_prefix + ".opp_str";

  string tmp_path = analysis_path + tmp_suffix;
  if(!directory_exists(tmp_path)){
    cout<<"Creating path "<<tmp_path<<endl;
    string dummy_path = tmp_path + "/none";
    create_file_path(dummy_path);
  }

  //copy genome table
  string _gt_fname = analysis_path + genome_table_suffix;
  create_file_path(_gt_fname);
  ofstream _gt_ofstr(_gt_fname.c_str());
  assert(_gt_ofstr.good());

  ifstream gt_ifstr(gt_fname.c_str());
  assert(gt_ifstr.good());
  while(gt_ifstr.good()){
    string line;
    getline(gt_ifstr, line);
    if(gt_ifstr.good()){
      _gt_ofstr<<line<<endl;
    }
  }
  gt_ifstr.close();
  _gt_ofstr.close();

  
  vector<string> fnames_to_del;

  
  
  genome_table gt(_gt_fname);

  int thread_num = 8;
  PSamAlignmentFile input_sam_file(input_sam_fname, 500000);
  input_sam_file.SetThreadNum(thread_num);


  ofstream* sam_ofstrs = new ofstream[gt.contigs.size()];

  for(unsigned int i=0; i<gt.contigs.size(); i++){
    string cur_chrom = gt.contigs[i];
    string cur_sam_fname = tmp_path + cur_chrom +".sam";
    cout<<"Opening "<<cur_sam_fname<<endl;

    sam_ofstrs[i].open(cur_sam_fname.c_str());
    if(!sam_ofstrs[i].good()){
      cout<<"Could not open "<<cur_sam_fname<<endl;
    }
    assert(sam_ofstrs[i].good());
    fnames_to_del.push_back(cur_sam_fname);
  }

  int counter=0;
  int valid = 0;

  while(input_sam_file.good()){
    if(counter%10000 == 0){
      cout<<"\r Lines: "<<((double)counter)/1000000.0<<" M; valid alignments: "<<((double)valid)/1000000.0<<" M    ";
      cout.flush();
    }
    SamAlignment& cur_al = input_sam_file.read_next_alignment();
    
    if(input_sam_file.good()){
      //cout<<"SAM: "<<cur_al.sam_string<<endl;
      //cout<<"valid: "<<cur_al.valid()<<endl;
      if(cur_al.valid()){
	string cur_chrom = cur_al.chrom.str();
	//cout<<"cur_chrom: "<<cur_chrom<<endl;
	int cur_ind = gt.contig_ind(cur_chrom);
	if(cur_ind == -1){
	  //conit_ind(string) returns an error
	  //cout<<"gt.contig_ind(string&) returns -1 for the SAM string"<<endl;
	  //cout<<cur_al.sam_string<<endl;
	  //cout<<"crhom: "<<cur_chrom<<endl;
	}
	if(cur_ind >= 0 && cur_ind <= (int) gt.contigs.size()){
	  sam_ofstrs[cur_ind]<<cur_al.sam_string<<endl;
	  valid++;
	}
      }
      //valid++;
    }
    counter++;
  }

  input_sam_file.close();
  cout<<endl<<"Found: "<<counter<<" lines and "<<valid<<" valid_alignments"<<endl; 
  //exit(0);

  /*
  for(unsigned int i=0; i<gt.contigs.size(); i++){
    string cur_chrom = gt.contigs[i];
    string cur_sam_fname = tmp_path + cur_chrom +".sam";
    cout<<"Opening "<<cur_sam_fname<<endl;

    sam_ofstrs[i].open(cur_sam_fname.c_str());
    assert(sam_ofstrs[i].good());

  }

  for(unsigned int i=0; i<gt.contigs.size(); i++){
    sam_ofstrs[i].close();
    }*/

  delete[] sam_ofstrs;  


  for(unsigned int i=0; i<gt.size(); i++){
    string cur_chrom = gt.contigs[i];
    int cur_chrom_size = gt.contig_sizes[i];

    string cur_sam_fname = tmp_path + cur_chrom + ".sam";
    
    PSamAlignmentFile input_sam_file(cur_sam_fname, 500000);
    input_sam_file.SetThreadNum(thread_num);

    string count_bin_fname = output_prefix + "." + cur_chrom; 
    ofstream count_bin_ofstr(count_bin_fname.c_str());
    if(!count_bin_ofstr.good()){
      cerr<<"Bad file: "<<count_bin_fname<<endl;
      exit(1);
    }
    
    cout<<"Opened "<<count_bin_fname<<endl;

    vector<int> start_counts_for(cur_chrom_size, 0);
    vector<int> start_counts_rev(cur_chrom_size, 0);

    for(unsigned j=0; j<start_counts_for.size(); j++){
      start_counts_for[j] = 0;
      start_counts_rev[j] = 0;
    }

    int total_reads = 0;
    int mapped_to_this_chr_reads = 0;
    while(input_sam_file.good()){
      if(total_reads % 100000 == 0){
	//cout<<"\rTotal: "<<commify(total_reads)<<" this chrom: "<<commify(mapped_to_this_chr_reads)<<" ";
	printf("\rTotal: %.2f M, this chrom: %.2f M ", ((double)total_reads)/1000000.0, 
	       ((double)mapped_to_this_chr_reads) / 1000000.0);
	cout.flush();
      }
      SamAlignment& cur_al = input_sam_file.read_next_alignment();
      if(input_sam_file.good()){
	total_reads++;
	if(cur_al.chrom == cur_chrom){
	  mapped_to_this_chr_reads++;
	  int coord_5p = -1;
	  
	  if(cur_al.strand == '+') coord_5p = cur_al.begin;
	  else coord_5p = cur_al.end - 1;
	  
	  if(coord_5p < 0 || coord_5p >= cur_chrom_size){
	    cout<<"Warning, coordinate out of bounds. skipping"<<endl;
	  } else {
	    if(cur_al.strand == '+') start_counts_for[coord_5p]++;
	    else start_counts_rev[coord_5p]++;
	  }
	}
      }
    }
    
    cout<<endl<<endl;

    input_sam_file.close();
    count_bin_ofstr.write((char*)(&cur_chrom_size), sizeof(cur_chrom_size));
    count_bin_ofstr.write((char*)(&start_counts_for[0]),sizeof(start_counts_for[0])*cur_chrom_size);
    count_bin_ofstr.write((char*)(&start_counts_rev[0]),sizeof(start_counts_rev[0])*cur_chrom_size);
    count_bin_ofstr.close();
  }

  for(unsigned int i=0; i<fnames_to_del.size(); i++){
    unlink(fnames_to_del[i].c_str());
  }

  return 0;
}
