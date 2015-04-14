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
  pars.optional("regions_file", "regions_file", "none", STRING_TYPE);
  pars.optional("regions_file_type", "regions_file_type", "none", STRING_TYPE);
  pars.optional("max_dist","max_dist","3000",INT_TYPE);
  //pars.optional("start_count_thresh", "start_count_thresh", "1", INT_TYPE);

  if(!pars.enforce()){
    exit(1);
  }

  //string regions_fname = pars.get_string_value("QuEST_regions_file");
  string analysis_path = pars.get_string_value("analysis_path");
  string regions_fname = pars.get_string_value("regions_file");
  int max_dist = pars.get_int_value("max_dist");
  
  if(analysis_path[analysis_path.size()-1] != '/') analysis_path = analysis_path + "/";
  string genome_table_fname = analysis_path + genome_table_suffix;
  string bin_count_file_prefix = analysis_path + align_bin_suffix;

  bool using_regions;
  if(regions_fname == "none") using_regions = false;
  else{
    using_regions = true;
    assert(false); //need to add the code here
  }

  string distogram_fname;
  string phasogram_fname;
  string frag_est_fname;

  if(!using_regions){
    distogram_fname = analysis_path + dist_plots_suffix + "whole_genome/" + "distogram.txt";
    phasogram_fname = analysis_path + dist_plots_suffix + "whole_genome/" + "phasogram.txt";
    frag_est_fname = analysis_path + dist_plots_suffix + "whole_genome/" + "fragment_estimate.txt";
    create_file_path(frag_est_fname);
    
  }
  //int start_count_thresh = pars.get_int_value("start_count_thresh");
  
  remove(distogram_fname.c_str());
  remove(phasogram_fname.c_str());
  //remove(frag_est_fname.c_str());
  
  //double bw = pars.get_double_value("bw");
  // cout<<"bw: "<<bw<<endl;
  

  genome_table gt(genome_table_fname);

  //int genome_size = 0;
  //for(unsigned int i=0; i<gt.size(); i++){
  //  genome_size += gt.chr_size[i];
  //}

  /*
  ifstream regions_str(regions_fname.c_str());
  if(!regions_str.good()){
    cerr<<"Bad regions file: "<<regions_fname<<endl;
    exit(1);
  }

  vector<region> regions;
  
  while(regions_str.good()){
    string cur_string;
    getline(regions_str, cur_string);
    if(regions_str.good() && cur_string.size() > 0){
      if(cur_string[0] == 'R'){
	char gap_char = ' ';
	char dash_char = '-';
	vector<string> cur_string_fields = split(cur_string, gap_char);
	if(cur_string_fields.size() >= 3){
	  string cur_chrom = cur_string_fields[1];
	  string range = cur_string_fields[2];
	  
	  vector<string> range_fields = split(range, dash_char);
	  if(range_fields.size() == 2){
	    int region_begin = atoi(range_fields[0].c_str());
	    int region_end = atoi(range_fields[1].c_str());
	    
	    region new_region;
	    new_region.chrom = cur_chrom;
	    new_region.begin = region_begin;
	    new_region.end = region_end;

	    regions.push_back(new_region);
	  }
	}
      }
    }
  }
  regions_str.close();

  cout<<"read "<<regions.size()<<" regions"<<endl;  
  */

  int max_pile = 10;

  vector<unsigned int> dummy_distogram(max_dist*2 + 1, 0);
  vector<unsigned int> dummy_phasogram(max_dist + 1, 0);
  
  vector< vector <unsigned int> > distograms(max_pile+1, dummy_distogram); 
  vector< vector <unsigned int> > phasograms(max_pile+1, dummy_phasogram); 

  //vector<unsigned int> opp_str_plot(max_dist*2+1); //center is at max_dist
  

  //unsigned long total_region_area = 0;
  //int positions_occupied=0; 
  //int dyads_within_regions = 0;
  
  for(unsigned int k=0; k<1; k++){//gt.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];

    vector<bool> chrom_mask(cur_chrom_size, true);

    /*
    for(unsigned int i=0; i<chrom_mask.size(); i++){
      chrom_mask[i] = false;
    }
    
    int regions_in_this_chrom = 0;
    for(unsigned int l=0; l<regions.size(); l++){
      if(regions[l].chrom == cur_chrom){
	regions_in_this_chrom++;
	//cout<<"found region"<<endl;
	for(int i=regions[l].begin; i<regions[l].end; i++){
	  chrom_mask[i] = true;
	}
      }
    }

    unsigned int bp_within_regions = 0;
    for(unsigned int i=0; i<chrom_mask.size(); i++){
      if(chrom_mask[i] == true){
	bp_within_regions++;
      }
    }
  
    total_region_area += bp_within_regions;

    double fraction_within_regions = (((double)bp_within_regions)/((double)chrom_mask.size()));
    cout<<endl;
    cout<<"chrom: "<<cur_chrom<<endl;
    cout<<"regions: "<<regions_in_this_chrom<<" occupy: "<<100*fraction_within_regions<<" percent of this chromosome"<<endl;
    */
    
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

    /*
    cout<<"masking..."<<endl;
    for(unsigned int i=0; i<start_counts_for.size(); i++){
      if(start_counts_for[i] > 0){
	if(chrom_mask[i] == false){ 
	  start_counts_for[i] = 0; 
	}
	else{
	  positions_occupied ++;
	  dyads_within_regions += start_counts_for[i];
	}
      }
      if(start_counts_rev[i] > 0){
	if(chrom_mask[i] == false){ 
	  start_counts_rev[i] = 0; 	  
	}
	else{
	  positions_occupied ++; 
	  dyads_within_regions += start_counts_rev[i];
	}
      }
    }
    */

    cout<<"calculating distance plots..."<<endl;
    
    for(int k=1; k<=max_pile; k++){
      int pile_threshold = k;
      cout<<"Pile "<<k<<endl;
      for(unsigned int i=0; i<start_counts_for.size(); i++){
	if(i%10000 == 0){ 
	  printf( "\rstep 1 of 3: %.2f M / %.2f M", ((double)i)/1000000.0, ((double)start_counts_for.size())/1000000.0); cout.flush();}
	
	if(start_counts_for[i] >= pile_threshold){
	  for(unsigned int j=i; j<=min(i+max_dist,start_counts_for.size()); j++){
	    if(start_counts_for[j] >= pile_threshold){
	      int dist = ((int) j) - ((int) i);
	      phasograms[pile_threshold][dist]++;
	    }
	  }
	}
      }
    
      cout<<endl;
      for(unsigned int i=0; i<start_counts_rev.size(); i++){
	if(i%10000 == 0){ printf( "\rstep 2 of 3: %.2f M / %.2f M", ((double)i)/1000000.0, ((double)start_counts_for.size())/1000000.0); cout.flush();}
	if(start_counts_rev[i] >= pile_threshold){
	  for(unsigned int j=i; j<=min(i+max_dist,start_counts_rev.size()-1); j++){
	    if(start_counts_rev[j] >= pile_threshold){
	      int dist = ((int) j) - ((int) i);
	      phasograms[pile_threshold][dist]++;
	    }
	  }
	}
      }
      cout<<endl;
      
      for(int i=0; i<(int)start_counts_for.size(); i++){
	if(i%10000 == 0){ printf( "\rstep 3 of 3: %.2f M / %.2f M", ((double)i)/1000000.0, ((double)start_counts_for.size())/1000000.0); cout.flush();}
	if(start_counts_for[i] >= pile_threshold){
	  for(int j=max(0,i-max_dist); 
	      j<=min(i+max_dist, (int) start_counts_rev.size()-1); j++){
	    if(start_counts_rev[j] >= pile_threshold){
	      int dist  = ((int)j) - ((int)i);
	      distograms[pile_threshold][dist + max_dist] ++;
	      //opp_str_plot[dist + max_dist] += start_counts_rev[j];
	    }
	  }
	}
      }
    }
    cout<<endl;
  }    
   
  //double same_str_normalizer = ((double)positions_occupied) * (((double)positions_occupied) / (2*((double)total_region_area)));
  //double opp_str_normalizer = ((double)positions_occupied) * (((double)positions_occupied) / (((double)total_region_area)));
  //double dyad_frequency = ((double)dyads_within_regions) / ((double)(total_region_area));
  
  //cout<<"positions_above_thres: "<<positions_occupied<<endl;
  //cout<<"region area: "<<total_region_area<<endl;
  //cout<<"norm: "<<same_str_normalizer<<endl;
  //cout<<"dyad_frequency: "<<dyad_frequency<<endl;
  
  
  ofstream phasogram_ofstr(phasogram_fname.c_str());
  if(!phasogram_ofstr.good()){
    cerr<<"Bad file name: "<<phasogram_fname<<endl;
  }
  
  
  for(unsigned int i=0; i<phasograms[0].size(); i++){
    phasogram_ofstr<<i;
    for(unsigned int k=1; k<= (unsigned int) max_pile; k++){
      phasogram_ofstr<<"\t"<<phasograms[k][i];
    }
    phasogram_ofstr<<endl;
  }
  phasogram_ofstr.close();

  ofstream distogram_ofstr(distogram_fname.c_str());
  if(!distogram_ofstr.good()){
    cerr<<"Bad file name: "<<distogram_fname<<endl;
  }
  
  for(unsigned int i=0; i<distograms[0].size(); i++){
    distogram_ofstr<<((int)i)-((int)max_dist);
    for(unsigned int k=1; k<=(unsigned int) max_pile; k++){
      distogram_ofstr<<"\t"<<distograms[k][i];//<<" "<<((double)opp_str_plot[i])/opp_str_normalizer<<endl;
    }
    distogram_ofstr<<endl;
  }
  distogram_ofstr.close();

  ofstream frag_est_ofstr(frag_est_fname.c_str());

  int min_pos = 110;
  int max_pos = 190;


  vector<int> max_inds;
  
  for(unsigned int k=0; k<=(unsigned int) max_pile; k++){
    int cur_max_ind = -1;
    int cur_max_value = -1;
    for(int i=min_pos; i<max_pos; i++){
      int dist_pos = max_dist+i;
      if((int)distograms[k][dist_pos] > cur_max_value){
	cur_max_ind = i;
	cur_max_value = distograms[k][dist_pos];
      }
    }
    max_inds.push_back(cur_max_ind);
  }

  frag_est_ofstr<<"size_estimate: "<<max_inds.at(3);
  frag_est_ofstr<<endl;


  for(unsigned int i=1; i<max_inds.size(); i++){
    frag_est_ofstr<<"pile="<<i<<" estimate: "<<max_inds[i]<<endl;
  }
  
  frag_est_ofstr.close();

  cout<<endl;
  //cout<<"reads_processed: "<<reads_processed<<endl<<endl;

  return 0;
}
