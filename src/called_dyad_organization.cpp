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

  pars.require("positions_file", "positions_file", STRING_TYPE);
  pars.require("called_dyads_file", "called_dyads_file", STRING_TYPE);
  pars.require("genome_table", "genome_table", STRING_TYPE);
  pars.require("output_file", "outptu_file", STRING_TYPE);
  pars.optional("max_dist", "max_dist", "3000", INT_TYPE);
  pars.optional("directional", "yes/no", "no", STRING_TYPE);
  
  if(!pars.enforce()){
    exit(1);
  }


  string genomic_sites_fname = pars.get_string_value("positions_file");
  string dyads_fname = pars.get_string_value("called_dyads_file");
  string output_fname = pars.get_string_value("output_file");
  int max_dist = pars.get_int_value("max_dist");

  string gt_fname = pars.get_string_value("genome_table");
  string directional_string = pars.get_string_value("directional");

  assert(directional_string == "yes" || directional_string == "no");

  genome_table gt(gt_fname);

  bool directional;
  if(directional_string == "yes") directional = true;
  else directional = false;

  ofstream ofstr(output_fname.c_str());
  assert(ofstr.good());

  vector<int> peak_counts(2*max_dist+1, 0);

  //read the positions
  vector<genome_coord> genomic_sites;
  ifstream genomic_sites_ifstr(genomic_sites_fname.c_str());
  assert(genomic_sites_ifstr.good());
  string cur_line;
  char delim = '\t';
  while(genomic_sites_ifstr.good()){
    getline(genomic_sites_ifstr, cur_line);
    if(genomic_sites_ifstr.good()){
      vector<string> cur_line_fields = split(cur_line, delim);
      assert(cur_line_fields.size() >= 2);
      genome_coord new_site;
      new_site.chrom = cur_line_fields[0];
      new_site.pos = atoi(cur_line_fields[1].c_str());
      genomic_sites.push_back(new_site);
    }
  }
  genomic_sites_ifstr.close();

  cout<<"Read "<<genomic_sites.size()<<" genomic sites"<<endl;

  vector<genome_coord> dyad_positions;
  ifstream dyads_ifstr(dyads_fname.c_str());
  assert(dyads_ifstr.good());
  while(dyads_ifstr.good()){
    getline(dyads_ifstr, cur_line);
    if(dyads_ifstr.good()){
      vector<string> cur_line_fields = split(cur_line, delim);
      assert(cur_line_fields.size() >= 2);
      genome_coord new_site;
      new_site.chrom = cur_line_fields[0];
      new_site.pos = atoi(cur_line_fields[1].c_str());
      dyad_positions.push_back(new_site);
    }
  }
  dyads_ifstr.close();
  cout<<"read "<<dyad_positions.size()<<" dyad positions"<<endl;
  
  for(unsigned int k=0; k<gt.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];
  
    vector<bool> peak_mask(cur_chrom_size, false);
    int peaks = 0;
    for(unsigned int i=0; i<dyad_positions.size(); i++){
      if(dyad_positions[i].chrom == cur_chrom){
	peak_mask.at(dyad_positions[i].pos) = true;
	peaks++;
      }
    }
    cout<<"found "<<commify(peaks)<<" peaks"<<endl;

    int peaks_included = 0;
    int sites = 0;
    for(unsigned int l=0; l<genomic_sites.size(); l++){
      sites++;
      if(genomic_sites[l].chrom == cur_chrom){
	int pos = genomic_sites[l].pos;

	bool flip = false;
	if(directional){
	  //find out if left or right peak is closer
	  int left_peak_offset = -1;
	  int right_peak_offset = -1;
	  for(int i=pos; i>=pos-max_dist && left_peak_offset == -1; i--){
	    if(i>=0){
	      if(peak_mask[i] == true) left_peak_offset = i-pos;
	    }
	  }
	  for(int i=pos; i<=pos+max_dist && right_peak_offset == -1; i++){
	    if(i<=cur_chrom_size){
	      if(peak_mask[i] == true) right_peak_offset = pos - i;
	    }
	  }
	  if(left_peak_offset < right_peak_offset && left_peak_offset !=-1) flip = true;
	}

	for(int i=pos; i<pos+max_dist; i++){
	  if(i>=0 && i<cur_chrom_size){
	    if(peak_mask[i] == true){
	      int offset;
	      if(flip) offset = max_dist + i-pos;
	      else offset = max_dist + (pos-i);
	      peak_counts.at(offset)++;
	      peaks_included++;
	    }
	  }
	}
	for(int i=pos; i>=pos-max_dist; i--){
	  if(i>=0 && i<cur_chrom_size){
	    if(peak_mask[i] == true){
	      int offset;
	      if(flip) offset = max_dist + (i-pos);
	      else offset = max_dist + (pos-i);
	      peak_counts.at(offset)++;
	      peaks_included++;
	    }
	  }
	}
      }
    }
    cout<<"sites: "<<sites<<endl;
    cout<<"peaks included: "<<peaks_included<<endl;
  }
    
  for(int j=0; j<(int) peak_counts.size(); j++){
    int offset = j - max_dist;
    ofstr<<offset<<"\t"<<peak_counts[j]<<endl;
  }

  //cout<<"under_pos_nucl_bases: "<<under_pos_nucl_bases<<endl;
  //cout<<"total_bases: "<<total_bases<<endl;
  //cout<<"% of bases under positioned nucleosomes: "<<100*(((double)under_pos_nucl_bases)/((double)total_bases))<<endl;
  
  ofstr.close();
  return 0;
}
