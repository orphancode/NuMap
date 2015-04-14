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
  pars.require("coordinates_file", "chrom\tcoord", STRING_TYPE);
  pars.require("output_file", "output_file", STRING_TYPE);
  pars.optional("max_dist","max_dist","3000",INT_TYPE);
  //pars.optional("core_size","core_size","147",INT_TYPE);
  //pars.optional("bw", "bw", "50", INT_TYPE); //KDE bandwidth

  if(!pars.enforce()){
    exit(1);
  }

  string analysis_path = pars.get_string_value("analysis_path");
  string output_fname = pars.get_string_value("output_file");
  string coordinates_fname = pars.get_string_value("coordinates_file"); //this is where TFBS coordinates are stored

  if(analysis_path[analysis_path.size()-1]!='/')
    analysis_path = analysis_path + "/";
  
  string dyad_bin_prefix = analysis_path + dyad_bin_suffix;
  string genome_table_fname = analysis_path + genome_table_suffix;
  
  ofstream dyad_hist_ofstr(output_fname.c_str());
  assert(dyad_hist_ofstr.good());
  
  //string dyad_density_fname = output_prefix + ".dyad_density.txt";
  //remove(dyad_density_fname.c_str());
  //ofstream dyad_density_ofstr(dyad_density_fname.c_str());
  
  int max_dist = pars.get_int_value("max_dist");
  //int core_size = pars.get_int_value("core_size");
  //int bw = pars.get_double_value("bw");

  genome_table gt(genome_table_fname);

  vector<vector<int> > peak_coords;

  for(unsigned int i=0; i<gt.size(); i++){
    vector<int> dummy_vec;
    peak_coords.push_back(dummy_vec);
  }

  ifstream coord_ifstr(coordinates_fname.c_str());
  assert(coord_ifstr.good());
  
  vector<genome_coord> genomic_sites;
  string cur_line;
  char delim = '\t';
  while(coord_ifstr.good()){
    getline(coord_ifstr, cur_line);
    if(coord_ifstr.good()){
      vector<string> cur_line_fields = split(cur_line, delim);
      assert(cur_line_fields.size() >= 2);
      genome_coord cur_coord;
      cur_coord.chrom = cur_line_fields[0];
      cur_coord.pos = atoi(cur_line_fields[1].c_str());
      genomic_sites.push_back(cur_coord);
    }
  }
  coord_ifstr.close();

  cout<<"Read "<<genomic_sites.size()<<" coordinates"<<endl;

  vector<int> dyad_hist(2*max_dist + 1, 0);

  for(unsigned int k=0; k<gt.size(); k++){
    string cur_chrom = gt.contigs[k];
    int cur_chrom_size = gt.contig_sizes[k];
    string dyad_bin_fname = dyad_bin_prefix + "." + cur_chrom;
    
    ifstream dyad_ifstr(dyad_bin_fname.c_str());
    if(!dyad_ifstr.good()){
      cerr<<"Bad file: "<<dyad_bin_fname<<endl;
      exit(1);
    }
    
    cout<<endl<<"processing "<<dyad_bin_fname<<endl;
    
    int chr_size;
    dyad_ifstr.read((char*)(&chr_size), sizeof(chr_size));
    assert(chr_size == cur_chrom_size);
    
    vector<int> dyad_counts(cur_chrom_size);

    dyad_ifstr.read((char*)(&dyad_counts[0]),sizeof(dyad_counts[0])*chr_size);
    
    dyad_ifstr.close();

    for(unsigned int j=0; j<genomic_sites.size(); j++){
      if(genomic_sites[j].chrom == cur_chrom){
	int cur_site_coord = genomic_sites[j].pos;
	for(int i=cur_site_coord - max_dist; i<=cur_site_coord + max_dist; i++){
	  if(i > 0 && i<(int) cur_chrom_size){
	    int cur_offset = i-cur_site_coord + max_dist;
	    dyad_hist[cur_offset] += dyad_counts[i];
	  }
	}
      }
    }
  }

  //ofstream dyad_hist_ofstr(output_fname.c_str());
  for(int i=0; i<(int)dyad_hist.size(); i++){
    dyad_hist_ofstr<<i-max_dist<<"\t"<<dyad_hist.at(i)<<endl;
  }
  dyad_hist_ofstr.close();

  return 0;
  
  /*
  //calculate positioning stringency
  vector<double> dyad_hist_kde(dyad_hist_for.size(),0);
  for(int i=0; i<(int) dyad_hist_for.size(); i++){
    for(int j=i-bw; j<=i+bw; j++){
      if(j>=0 && j<(int) dyad_hist_for.size()){
	dyad_hist_kde[j] += (dyad_hist_for[i] + dyad_hist_rev[i]) * triweight_kernel((double)(i-j), bw);
      }
    }
  }

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

  return 0;
  */
}
