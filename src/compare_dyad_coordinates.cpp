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
};

int assign_matches(int dist, vector<bool>& dyad_mask1, 
		   vector<bool>& dyad_mask2, vector<int>& match_mask){ 
  int res = 0;
  int cur_chrom_size = (int) dyad_mask1.size();
  for(int i=0; i<(int)dyad_mask1.size(); i++){
    //if(dyad_mask1[i] == true && dyad_mask2[i] == true)
    //  cout<<"Bingo!"<<endl;
    if(dyad_mask1[i] == true){
      //cout<<"here"<<endl;
      int closest_peak = -1;
      bool _continue = true;
      int cur_offset = 0;
      while(_continue){
	if(i-cur_offset >= 0 && i-cur_offset < cur_chrom_size){
	  if(dyad_mask2[i - cur_offset] == true){
	    closest_peak = i-cur_offset; _continue = false;
	  }
	}
	if(i+cur_offset >= 0 && i+cur_offset < cur_chrom_size){
	  if(dyad_mask2[i + cur_offset] == true){
	    closest_peak = i+cur_offset; _continue = false;
	  }
	}
	cur_offset++;
	if(cur_offset > dist) _continue = false;
      }
      if(closest_peak != -1){ match_mask[i] = closest_peak; res++; }
    }
  }
  return res;
}

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("called_dyads_file1", "called_dyads_file1", STRING_TYPE);
  pars.require("called_dyads_file2", "called_dyads_file2", STRING_TYPE);
  pars.require("genome_table", "genome_table", STRING_TYPE);
  pars.require("output_file", "outptu_file", STRING_TYPE);
  pars.optional("max_dist", "max_dist", "100", INT_TYPE);
  
  if(!pars.enforce()){
    exit(1);
  }


  string dyads_fname1 = pars.get_string_value("called_dyads_file1");
  string dyads_fname2 = pars.get_string_value("called_dyads_file2");

  string output_fname = pars.get_string_value("output_file");
  string score_sample_fname = output_fname + ".score_sample";
  int max_dist = pars.get_int_value("max_dist");

  string gt_fname = pars.get_string_value("genome_table");//analysis_path + genome_table_suffix;
  genome_table gt(gt_fname);

  ofstream ofstr(output_fname.c_str());
  assert(ofstr.good());
  
  //read the positions
  vector<genome_coord> dyads1;
  ifstream dyads_ifstr1(dyads_fname1.c_str());
  assert(dyads_ifstr1.good());

  string cur_line;
  char delim = '\t';
  while(dyads_ifstr1.good()){
    getline(dyads_ifstr1, cur_line);
    if(dyads_ifstr1.good()){
      vector<string> cur_line_fields = split(cur_line, delim);
      assert(cur_line_fields.size() >= 2);
      genome_coord new_site;
      new_site.chrom = cur_line_fields[0];
      new_site.pos = atoi(cur_line_fields[1].c_str());
      new_site.value = atof(cur_line_fields[2].c_str());
      dyads1.push_back(new_site);
    }
  }
  dyads_ifstr1.close();
  cout<<"Read "<<dyads1.size()<<" dyads from file 1"<<endl;

  vector<genome_coord> dyads2;
  ifstream dyads_ifstr2(dyads_fname2.c_str());
  assert(dyads_ifstr2.good());

  while(dyads_ifstr2.good()){
    getline(dyads_ifstr2, cur_line);
    if(dyads_ifstr2.good()){
      vector<string> cur_line_fields = split(cur_line, delim);
      assert(cur_line_fields.size() >= 2);
      genome_coord new_site;
      new_site.chrom = cur_line_fields[0];
      new_site.pos = atoi(cur_line_fields[1].c_str());
      new_site.value = atof(cur_line_fields[2].c_str());
      dyads2.push_back(new_site);
    }
  }
  dyads_ifstr2.close();
  cout<<"Read "<<dyads2.size()<<" dyads from file 2"<<endl;

  vector< vector<double> > score_sample1;
  vector< vector<double> > score_sample2;
  
  
  vector<int> matched_dyads;
  //vector<int> mock_matched_dyads;
  vector<int> match_dist;
  int step = 5;
  vector<double> dummy_vec;
  for(int i=0; i<=max_dist; i+=step){
    matched_dyads.push_back(0);
    //mock_matched_dyads.push_back(0);
    match_dist.push_back(i);
    score_sample1.push_back(dummy_vec);
    score_sample2.push_back(dummy_vec);
  }
  
  for(unsigned int k=0; k<gt.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];
    cout<<"Processing "<<cur_chrom<<endl;
    
    vector<bool> dyad_mask1(cur_chrom_size, false);
    vector<double> dyad_scores1(cur_chrom_size, 0);
    int peaks1 = 0;
    for(unsigned int i=0; i<dyads1.size(); i++){
      if(dyads1[i].chrom == cur_chrom){
	dyad_mask1.at(dyads1[i].pos) = true;
	dyad_scores1.at(dyads1[i].pos) = dyads1[i].value;
	peaks1++;
      }
    }
    cout<<"found "<<commify(peaks1)<<" peak in file1"<<endl;
    
    vector<bool> dyad_mask2(cur_chrom_size, false);
    vector<double> dyad_scores2(cur_chrom_size, 0);
    int peaks2 = 0;
    for(unsigned int i=0; i<dyads2.size(); i++){
      if(dyads2[i].chrom == cur_chrom){
	dyad_mask2.at(dyads2[i].pos) = true;
	dyad_scores2.at(dyads2[i].pos) = dyads2[i].value;
	peaks2++;
      }
    }
    cout<<"found "<<commify(peaks2)<<" peaks in file2"<<endl;

    int dind = 0;
    for(int d=0; d<=max_dist; d+=step){
      cout<<"cur_dist: "<<d<<endl;
      vector<int> match_mask1(cur_chrom_size,-1);
      vector<int> match_mask2(cur_chrom_size,-1);
      //int matches1 = assign_matches(d, dyad_mask1, dyad_mask2, match_mask1);
      //int matches2 = assign_matches(d, dyad_mask1, dyad_mask1, match_mask2);
      assign_matches(d, dyad_mask1, dyad_mask2, match_mask1);
      assign_matches(d, dyad_mask2, dyad_mask1, match_mask2);
      //cout<<"matches_computed"<<endl;

      int matches = 0;
      for(int i=0; i<cur_chrom_size; i++){
	if(match_mask1[i] != -1){
	  int match_ind1 = match_mask1[i];
	  assert(match_ind1 >= 0 && match_ind1 < cur_chrom_size);
	  int match_ind2 = match_mask2[match_ind1];
	  assert(match_ind2 >= 0 && match_ind2 < cur_chrom_size);
	  if(match_ind2 == i) matches++;
	}
	//if(matches2[matches1[i]] == i) matches++;

      }
      matched_dyads[dind] += matches;
      cout<<"matches: "<<matches<<endl;

      if(matches > 1000){
	int divisor = matches/1000;
	int counter = 0;
	for(int i=0; i<cur_chrom_size; i++){
	  if(match_mask1[i] != -1){
	    int match_ind1 = match_mask1[i];
	    int match_ind2 = match_mask2[match_ind1];
	    if(match_ind2 == i){
	      if(counter%divisor==0){
		double cur_value1 = dyad_scores1[i];
		double cur_value2 = dyad_scores2[match_ind1];
		score_sample1[dind].push_back(cur_value1);
		score_sample2[dind].push_back(cur_value2);
		//cout<<"adding"<<endl;
	      }
	    }
	    counter++;
	  }
	}
      }
					  
      
      dind++;
      /*
      for(int i=0; i<cur_chrom_size; i++){
	if(dyad_mask1[i] == true){
	  int closest_peak = -1;
	  bool _continue = 0;
	  int cur_offset = 0;
	  while(_continue){
	    if(cur_offset > d) _continue = false;
	    if(i-cur_offset >= 0 && i-cur_offset < cur_chrom_size){
	      if(dyad_mask2[i - cur_offset] == true){
		closest_peak = i-cur_offset; _continue = false;
	      }
	    }
	    if(i+cur_offset >= 0 && i+cur_offset < cur_chrom_size){
	      if(dyad_mask2[i + cur_offset] = true){
		closest_peak = i+cur_offset; _continue = false;
	      }
	    }
	    cur_offset++;
	  }
	  if(closest_peak != -1) match_mask1[closest_peak] = i;
	}
	}*/
    }
    //cout<<"found "<<commify(peaks2)<<" peaks in file2"<<endl;
  }


  for(unsigned int i=0; i<matched_dyads.size(); i++){
    ofstr<<match_dist[i]<<"\t"<<matched_dyads[i]<<endl;
  }
  ofstr.close();
  
  ofstream score_sample_ofstr(score_sample_fname.c_str());
  for(unsigned int i=0; i<match_dist.size(); i++){
    if(i!=0) score_sample_ofstr<<"\t";
    else score_sample_ofstr<<"#";
    score_sample_ofstr<<"dist_"<<match_dist[i]<<"_1"<<"\t"<<"dist_"<<match_dist[i]<<"_2";
  }
  cout<<endl;
  
  int min_size = -1;
  for(unsigned int j=0; j<score_sample1.size(); j++){
    if(j==0) min_size = score_sample1[j].size();
    else{
      if(min_size > (int)score_sample1[j].size()) min_size = score_sample1[j].size();
    }
  }
  cout<<"min_size = "<<min_size<<endl;
  for(int i=0; i<min_size; i++){
    for(int j=0; j<(int)score_sample1.size(); j++){
      if(j!=0) score_sample_ofstr<<"\t";
      score_sample_ofstr<<score_sample1[j][i]<<"\t"<<score_sample2[j][i];
    }
    score_sample_ofstr<<endl;
  }

  
  score_sample_ofstr.close();
  return 0;
}
