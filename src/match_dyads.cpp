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
  double value1;
  double value2;
  string record;
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
  //pars.require("output_file", "output_file", STRING_TYPE);
  pars.optional("match_dist", "match_dist", "50", INT_TYPE);
  pars.optional("output_path", "output_path", "none", STRING_TYPE);
  
  if(!pars.enforce()){
    exit(1);
  }

  precomputed_kernel kernel(5000, "triweight");

  string dyads_fname1 = pars.get_string_value("called_dyads_file1");
  string dyads_fname2 = pars.get_string_value("called_dyads_file2");

  string gap_fname1 = get_file_path(dyads_fname1) + "/gaps.txt";
  string gap_fname2 = get_file_path(dyads_fname2) + "/gaps.txt";

  string output_path = pars.get_string_value("output_path");
  assert(output_path != "none");
  if(output_path[output_path.size()-1] != '/') output_path = output_path + "/";

  string output_fname = output_path + "dyad_matches"; //pars.get_string_value("output_file");
  string output_fname1 = output_path + "/dyad_positions1.txt";
  string output_fname2 = output_path + "/dyad_positions2.txt";
  string unm_fname1 = output_path + "/unmatched_dyads1";
  string unm_fname2 = output_path + "/unmatched_dyads2";
  //string score_sample_fname = output_path + ".score_sample";
  int match_dist = pars.get_int_value("match_dist");

  bool using_bed_output = true;

  string gt_fname = pars.get_string_value("genome_table");//analysis_path + genome_table_suffix;

  genome_table gt(gt_fname);

  create_file_path(output_path + "/dummy");
  ofstream ofstr(output_fname.c_str());
  if(!ofstr.good()){
    cout<<"Error. Wrong output file: "<<endl;
    cout<<output_fname<<endl;
    cout<<"Check the output_path variable"<<endl;
    exit(1);
  }
  //assert(ofstr.good());

  ofstream unm_ofstr1(unm_fname1.c_str());
  ofstream unm_ofstr2(unm_fname2.c_str());

  ofstream m_ofstr1(output_fname1.c_str());
  ofstream m_ofstr2(output_fname2.c_str());

  //read the positions
  vector<genome_coord> dyads1;
  cout<<"Processing "<<dyads_fname1<<endl;
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
      new_site.value1 = atof(cur_line_fields[2].c_str());
      new_site.value2 = atof(cur_line_fields[3].c_str());
      new_site.record = cur_line;
      dyads1.push_back(new_site);
    }
  }
  dyads_ifstr1.close();
  cout<<"Read "<<dyads1.size()<<" dyads from file 1"<<endl;

  cout<<"Reading gaps from "<<gap_fname1<<endl;
  ifstream gap_ifstr1(gap_fname1.c_str());
  if(!gap_ifstr1.good()){
    cout<<"Failed to locate "<<gap_fname1<<endl;
    exit(1);
  }

  vector<region> gaps1;
  while(gap_ifstr1.good()){
    getline(gap_ifstr1, cur_line);
    if(gap_ifstr1.good()){
      if(cur_line.length() > 0){
	if(cur_line[0] != '#'){
	  vector<string> cur_line_fields = split(cur_line, delim);
	  assert(cur_line_fields.size() >= 3);
	  region new_region;
	  new_region.chrom = cur_line_fields[0];
	  new_region.begin = atoi(cur_line_fields[1].c_str());
	  new_region.end = atoi(cur_line_fields[2].c_str());
	  gaps1.push_back(new_region);
	}
      }
    }
  }
  cout<<"Read "<<gaps1.size()<<" gap regions from file 1"<<endl;
  

  cout<<"Processing "<<dyads_fname2<<endl;
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
      new_site.value1 = atof(cur_line_fields[2].c_str());
      new_site.value2 = atof(cur_line_fields[3].c_str());
      new_site.record = cur_line;
      dyads2.push_back(new_site);
    }
  }
  dyads_ifstr2.close();
  cout<<"Read "<<dyads2.size()<<" dyads from file 2"<<endl;
  cout<<"Reading gaps from "<<gap_fname1<<endl;

  ifstream gap_ifstr2(gap_fname2.c_str());
  if(!gap_ifstr2.good()){
    cout<<"Failed to locate "<<gap_fname2<<endl;
    exit(1);
  }

  vector<region> gaps2;
  while(gap_ifstr2.good()){
    getline(gap_ifstr2, cur_line);
    if(gap_ifstr2.good()){
      if(cur_line.length() > 0){
	if(cur_line[0] != '#'){
	  vector<string> cur_line_fields = split(cur_line, delim);
	  assert(cur_line_fields.size() >= 3);
	  region new_region;
	  new_region.chrom = cur_line_fields[0];
	  new_region.begin = atoi(cur_line_fields[1].c_str());
	  new_region.end = atoi(cur_line_fields[2].c_str());
	  gaps2.push_back(new_region);
	}
      }
    }
  }
  cout<<"Read "<<gaps2.size()<<" gap regions from file 2"<<endl;

  
  for(unsigned int k=0; k<gt.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];
    cout<<"Processing "<<cur_chrom<<endl;

    ofstream bed_ofstr;
    if(using_bed_output){
      string bed_output_fname = output_path + "/matching_dyads." + cur_chrom + ".bed";
      bed_ofstr.open(bed_output_fname.c_str());
      assert(bed_ofstr.good());
      bed_ofstr<<"track name=\"matching_dyads";
      bed_ofstr<<"\" description=\"matching_dyads\"visibility=2 colorByStrand=\"0,0,255 0,0,255\""<<endl;
    }

    string diff_fname = output_path + "/dyad_diff.bin." + cur_chrom;
    ofstream diff_ofstr(diff_fname.c_str());
    assert(diff_ofstr.good());

    vector<int> inds1(cur_chrom_size, -1);
    vector<int> inds2(cur_chrom_size, -1);

    vector<bool> gap_mask1(cur_chrom_size, false);
    for(unsigned int i=0; i<gaps1.size(); i++){
      if(gaps1[i].chrom == cur_chrom){
	for(int j=gaps1[i].begin; j<gaps1[i].end; j++){
	  if(j>0 && j<cur_chrom_size)
	    gap_mask1[j] = true;
	}
      }
    }

    vector<bool> gap_mask2(cur_chrom_size, false);
    for(unsigned int i=0; i<gaps2.size(); i++){
      if(gaps2[i].chrom == cur_chrom){
	for(int j=gaps2[i].begin; j<gaps2[i].end; j++){
	  if(j>0 && j<cur_chrom_size)
	    gap_mask2[j] = true;
	}
      }
    }

    
    vector<bool> dyad_mask1(cur_chrom_size, false);
    //vector<double> dyad_scores1(cur_chrom_size, 0);
    int peaks1 = 0;
    for(unsigned int i=0; i<dyads1.size(); i++){
      if(dyads1[i].chrom == cur_chrom){
	dyad_mask1.at(dyads1[i].pos) = true;
	//dyad_scores1.at(dyads1[i].pos) = dyads1[i].value;
	peaks1++;
	inds1[dyads1[i].pos] = i;
      }
    }
    cout<<"found "<<commify(peaks1)<<" peak in file1"<<endl;
    
    vector<bool> dyad_mask2(cur_chrom_size, false);
    //vector<double> dyad_scores2(cur_chrom_size, 0);
    int peaks2 = 0;
    for(unsigned int i=0; i<dyads2.size(); i++){
      if(dyads2[i].chrom == cur_chrom){
	dyad_mask2.at(dyads2[i].pos) = true;
	//dyad_scores2.at(dyads2[i].pos) = dyads2[i].value;
	peaks2++;
	inds2[dyads2[i].pos] = i;
      }
    }
    cout<<"found "<<commify(peaks2)<<" peaks in file2"<<endl;

    vector<int> match_mask1(cur_chrom_size,-1);
    vector<int> match_mask2(cur_chrom_size,-1);
    assign_matches(match_dist, dyad_mask1, dyad_mask2, match_mask1);
    assign_matches(match_dist, dyad_mask2, dyad_mask1, match_mask2);
    
    //vector<double> match_scores(cur_chrom_size, 0);
    int matches = 0;
    for(int i=0; i<cur_chrom_size; i++){
      if(match_mask1[i] != -1){
	int match_ind1 = match_mask1[i];
	assert(match_ind1 >= 0 && match_ind1 < cur_chrom_size);
	int match_ind2 = match_mask2[match_ind1];
	assert(match_ind2 >= 0 && match_ind2 < cur_chrom_size);
	if(match_ind2 == i){
	  matches++;
	  int pos1 = i;
	  int pos2 = match_mask1[i];
	  int begin_pos = -1; 
	  int end_pos = -1;
	  if(pos1 < pos2){ begin_pos = pos1; end_pos = pos2; }
	  else { begin_pos = pos2; end_pos = pos1; }
	  if(using_bed_output){
	    int dist = end_pos - begin_pos;
	    double s1 = dyads1.at(inds1[pos1]).value1;
	    double s2 = dyads2.at(inds2[pos2]).value1;
	    double e1 = dyads1.at(inds1[pos1]).value2;
	    double e2 = dyads2.at(inds2[pos2]).value2;
	    bed_ofstr<<cur_chrom<<"\t"<<begin_pos+1<<"\t"<<end_pos+2<<"\t"<<dist<<"\t0\t+"<<endl;
	    ofstr<<cur_chrom<<"\t"<<pos1<<"\t"<<pos2<<"\t";
	    ofstr<<s1<<"\t"<<s2<<"\t"<<e1<<"\t"<<e2<<endl;
	    m_ofstr1<<dyads1[inds1[pos1]].record<<endl;
	    m_ofstr2<<dyads2[inds2[pos2]].record<<endl;
	  }
	}
      }
      if(match_mask1[i] == -1 && dyad_mask1[i] == true){ //unmatched dyad
	double s = dyads1.at(inds1[i]).value1; 
	double e = dyads1.at(inds1[i]).value2;
	unm_ofstr1<<cur_chrom<<"\t"<<i<<"\t"<<s<<"\t"<<e<<endl; //dyad_scores1[i]<<endl;
      }
      if(match_mask2[i] == -1 && dyad_mask2[i] == true){
	double s = dyads2.at(inds2[i]).value1; 
	double e = dyads2.at(inds2[i]).value2;
	unm_ofstr2<<cur_chrom<<"\t"<<i<<"\t"<<s<<"\t"<<e<<endl;//dyad_scores2[i]<<endl;
      }
      
    }
    cout<<"matches: "<<matches<<endl;

    bed_ofstr.close();

  }
  
  unm_ofstr1.close();
  unm_ofstr2.close();

  m_ofstr1.close();
  m_ofstr2.close();
  
  ofstr.close();
  return 0;
}
