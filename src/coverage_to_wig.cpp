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
};

int main(int argc, char* argv[]){

  params pars(argc, argv);

  pars.require("analysis_path", "analysis_path", STRING_TYPE);
  pars.require("positions_file", "positions_file", STRING_TYPE);
  pars.require("output_file", "output_file", STRING_TYPE);
  pars.optional("max_dist","max_dist","10000",INT_TYPE);
  pars.optional("track_name", "track_name", "dyad_stringency", STRING_TYPE);

  if(!pars.enforce()){
    exit(1);
  }

  string analysis_path = pars.get_string_value("analysis_path");
  string output_fname = pars.get_string_value("output_file");
  string positions_fname = pars.get_string_value("positions_file");
  int max_dist = pars.get_int_value("max_dist");
  string track_name = pars.get_string_value("track_name");

  if(analysis_path[analysis_path.size()-1]!='/')
    analysis_path = analysis_path + "/";
  
  string dyad_stringency_prefix = analysis_path + coverage_suffix;
  
  ofstream ofstr(output_fname.c_str());
  assert(ofstr.good());
  
  string gt_fname = analysis_path + genome_table_suffix;

  genome_table gt(gt_fname);

  ofstream wig_ofstr(output_fname.c_str());
  assert(wig_ofstr.good());

  ifstream sites_ifstr(positions_fname.c_str());
  assert(sites_ifstr.good());

  vector<genome_coord> genomic_sites;

  string cur_line;
  char delim = '\t';
  while(sites_ifstr.good()){
    getline(sites_ifstr, cur_line);
    if(sites_ifstr.good()){
      vector<string> cur_line_fields = split(cur_line, delim);
      if(!(cur_line_fields.size() >= 2)){ 
	cout<<"Error in "<<cur_line<<endl;
      }
      assert(cur_line_fields.size() >= 2);
      genome_coord new_site;
      new_site.chrom = cur_line_fields[0];
      new_site.pos = atoi(cur_line_fields[1].c_str());
      genomic_sites.push_back(new_site);
    }
  }
  
  //int region_counter = 0;
  cout<<"Read: "<<genomic_sites.size()<<" sites"<<endl;
  
  wig_ofstr<<"track type=wiggle_0 name=\"";
  wig_ofstr<<track_name<<"\" description=\"nucl_coverage\" visibility=full color=202,73,170 priority=35 maxHeightPixels=50:50:11"<<endl;
  for(unsigned int k=0; k<gt.contigs.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];
    cout<<"Processing "<<cur_chrom<<endl;

    string stringency_bin_file = dyad_stringency_prefix + "." + cur_chrom;
    ifstream stringency_ifstr(stringency_bin_file.c_str());
    if(!stringency_ifstr.good()){ cout<<"Could not open "<<stringency_bin_file<<endl; }
    assert(stringency_ifstr.good());
    int chr_size;

    stringency_ifstr.read((char*)(&chr_size), sizeof(chr_size));
    assert(chr_size == cur_chrom_size);
    vector<double> stringency_values(chr_size);
    stringency_ifstr.read((char*)(&stringency_values[0]), sizeof(double)*chr_size);
    stringency_ifstr.close();

    for(unsigned j=0; j<genomic_sites.size(); j++){
      if(genomic_sites[j].chrom == cur_chrom){
	int pos = genomic_sites[j].pos;
	//wig_ofstr<<"track type=wiggle_0 name=\"dyad_stringency\" description=\"dyad_stringency\" visibility=full color=202,73,170 priority=35 maxHeightPixels=50:50:11"<<endl;
	wig_ofstr<<"fixedStep chrom="<<cur_chrom<<" start="<<pos-max_dist<<" step=1 span=1"<<endl;
	for(int i=pos-max_dist; i<=pos+max_dist; i++){
	  char buf[50];
	  sprintf(buf, "%.2f", stringency_values.at(i));
	  //wig_ofstr<<i<<" "<<
	  wig_ofstr<<buf<<endl;
	}
	//if(region_counter==2){
	//  wig_ofstr.close();
	//  return 0;
	//}
	//region_counter++;
      }
    }
    
  }


  wig_ofstr.close();
  return 0;
}
