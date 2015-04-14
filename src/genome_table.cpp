#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include <stdlib.h>
#include "assert.h"
#include "string_utils.h"
#include "genome_table.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::pair;
using std::map;
using std::string;
using std::ifstream;
using std::istringstream;


genome_table::genome_table(){
}

genome_table::genome_table(string gt_fname){
  load(gt_fname);
  for(unsigned int i=0; i<contigs.size(); i++){
    contig_map[contigs[i]] = i;
  }
}

unsigned int genome_table::size() const {
  return contigs.size();
}

bool genome_table::load(string gt_fname){
  ifstream gt_ifstr(gt_fname.c_str());
  //cout<<"gt_fname: "<<gt_fname<<endl;
  if(!gt_ifstr.good()){
    cerr<<"Error in bool genome_table::load"<<endl;
    cerr<<"Bad file name: "<<gt_fname<<endl;
    exit(1);
  }
  else{
    while(gt_ifstr.good()){
      string cur_line;
      getline(gt_ifstr, cur_line);
      if(gt_ifstr.good()){

	istringstream cur_line_str(cur_line, istringstream::in);

	string cur_chr;
	int cur_chr_size;

	assert(cur_line_str.good());
	if(cur_line.length() > 0){
	  if( cur_line.substr( 0, 1) != "#"){
	    cur_line_str>>cur_chr;
	    assert(cur_line_str.good());
	    cur_line_str>>cur_chr_size;

	    contigs.push_back(cur_chr);
	    contig_sizes.push_back(cur_chr_size);
	  }
	}
      }
    }
  }

  //for(unsigned int i=0; i<contigs.size(); i++){
  //  cout<<contigs[i]<<" -> "<<contig_sizes[i]<<endl;
  //}

  
  return true;
}

bool genome_table::compare_coordinates(const string& chrom1, const int& pos1, const char& strand1,
				       const string& chrom2, const int& pos2, const char& strand2){
  
  if(chrom1 == chrom2){
    if( pos1<pos2 ){ return true; }
    else if (pos1 > pos2){ return false;}
    else{ //pos1 = pos2
      if(strand1 == '+' && strand2 == '-') return true;
      return false;
    }
  }

  //chrom1 != chrom2

  int ind1 = contig_ind(chrom1);
  int ind2 = contig_ind(chrom2);

  if(ind1 != -1 && ind2 != -1) return ind1 < ind2;
  else if (ind1 != -1 && ind2 == -1) return true; 
  else if (ind1 == -1 && ind2 != -1) return false;
  else{
    //ind1 == -1 && ind2 == -1 && chrom1 != chrom2
    return chrom1 < chrom2;
  }
  assert(false); //you should not be here
  return false;
}

bool genome_table::compare_coordinates(const StringView& chrom1, const int& pos1, const char& strand1,
				       const StringView& chrom2, const int& pos2, const char& strand2){
  if(chrom1 == chrom2){
    if( pos1<pos2 ){ return true; }
    else if (pos1 > pos2){ return false;}
    else{ //pos1 = pos2
      if(strand1 == '+' && strand2 == '-') return true;
      return false;
    }
  }

  //chrom1 != chrom2

  int ind1 = contig_ind(chrom1.str());
  int ind2 = contig_ind(chrom2.str());

  if(ind1 != -1 && ind2 != -1) return ind1 < ind2;
  else if (ind1 != -1 && ind2 == -1) return true; 
  else if (ind1 == -1 && ind2 != -1) return false;
  else{
    //ind1 == -1 && ind2 == -1 && chrom1 != chrom2
    return chrom1 < chrom2;
  }
  assert(false); //you should not be here
  return false;
}

bool genome_table::compare_coordinates(const StringView& chrom1, const int& pos1, 
				       const StringView& chrom2, const int& pos2){
  if(chrom1 == chrom2){
    if( pos1<pos2 ){ return true; }
    else if (pos1 > pos2){ return false;}
    else{ //pos1 = pos2
      return false;
    }
  }

  //chrom1 != chrom2

  int ind1 = contig_ind(chrom1.str());
  int ind2 = contig_ind(chrom2.str());

  if(ind1 != -1 && ind2 != -1) return ind1 < ind2;
  else if (ind1 != -1 && ind2 == -1) return true; 
  else if (ind1 == -1 && ind2 != -1) return false;
  else{
    //ind1 == -1 && ind2 == -1 && chrom1 != chrom2
    return chrom1 < chrom2;
  }
  assert(false); //you should not be here
  return false;
}

/*
bool genome_table::compare_coordinates(const string& chrom1, int pos1,
				       const string& chrom2, int pos2){
  if(chrom1 == chrom2){
  return pos1 < pos2;
  }

  //chrom1 != chrom2

  int ind1 = contig_ind(chrom1);
  int ind2 = contig_ind(chrom2);

  if(ind1 != -1 && ind2 != -1) return ind1 < ind2;
  else if (ind1 != -1 && ind2 == -1) return true; 
  else if (ind1 == -1 && ind2 != -1) return false;
  else{
    //ind1 == -1 && ind2 == -1 && chrom1 != chrom2
    return chrom1 < chrom2;
  }
  assert(false); //you should not be here
  return false;
  }
*/

int genome_table::contig_ind(const string& contig){
  if(contig_map.size() == 0){
    cout<<"Error. Contig map empty in genome table."<<endl;
    exit(1);
  }
  map<string,int>::iterator it = contig_map.find(contig);
  if(it == contig_map.end()){
    //cout<<"Warning in genome_table::contig_ind(const string&)"<<endl;
    //cout<<"could not find  "<<contig<<endl;
    return -1;
  }

  //cout<<"printing genome table map:"<<endl;
  //for(map<string,int>::iterator it=contig_map.begin(); it!= contig_map.end(); it++){
  //  cout<<it->first<<" "<<it->second<<endl;
  //}
  
  if(it == contig_map.end()) return -1; //element not found
 
  int res = it->second; 
  assert(res >= 0 && res < (int) contig_map.size());

  return res;
  
  /*
  for(unsigned int i=0; i<contigs.size(); i++){
    if(contig == contigs[i]){
      return i;
    }
  }
  return -1;
  */
}

int genome_table::contig_size(const string& contig){
  int _contig_ind = contig_ind(contig);
  if(_contig_ind != -1) return contig_sizes[_contig_ind];
  else return -1;
  /*
  for(unsigned int i=0; i<contigs.size(); i++){
    if(contigs[i] == contig){
      return contig_sizes[i];
    }
  }
  return -1;
  */
}

bool chrom_name_less(const chrom_info& lhs, const chrom_info& rhs){
  return lhs.chrom < rhs.chrom;
}

genome_hash_table::genome_hash_table(genome_table gt, int _seed_size){
  cout<<"constructing hash table"<<endl;
  //gt = _gt;
  for(unsigned int i=0; i<gt.contigs.size(); i++){
    chrom_info cur_chrom;
    cur_chrom.chrom = gt.contigs.at(i);
    cur_chrom.size = gt.contig_sizes.at(i);

    chroms.push_back(cur_chrom);
  }

  sort(chroms.begin(), chroms.end(), chrom_name_less);
  seed_size = _seed_size;

  int offset = 0; 
  for(unsigned int i=0; i<chroms.size(); i++){
    chroms.at(i).offset = offset;
    int cur_buckets = chroms.at(i).size / seed_size + 1; //one extra bucket to hold incomplete last part
    chroms.at(i).buckets = cur_buckets;
    offset += cur_buckets;
  }
}

int genome_hash_table::hash_value(string chrom, int coord){
  //first find index of this chrom using a binary search
  chrom_info cur_chrom;
  cur_chrom.chrom = chrom;

  pair<vector<chrom_info>::iterator, vector<chrom_info>::iterator> bounds;
  bounds = equal_range(chroms.begin(), chroms.end(), cur_chrom, chrom_name_less);
  
  if(bounds.first != bounds.second - 1){
    //chrom is not on the list
    return -1; 
  }
  
  int chrom_ind = int(bounds.first - chroms.begin());
  int cur_offset = coord / seed_size;
  assert(cur_offset - chroms.at(chrom_ind).buckets);
  int cur_ind = chroms[chrom_ind].offset + cur_offset;
  return cur_ind;
}

int genome_hash_table::size(){
  int res = 0;
  for(unsigned int i=0; i<chroms.size(); i++){
    res += chroms.at(i).buckets;
  }
  return res;
}

int genome_hash_table::dist_to_left_boundary(string chrom, int coord){
  int cur_dist = coord%seed_size;
  return cur_dist;
}
int genome_hash_table::dist_to_right_boundary(string chrom, int coord){
  int cur_dist = seed_size - coord%seed_size;
  return cur_dist;
}
