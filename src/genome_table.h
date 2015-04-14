#ifndef GENOME_TABLE_H
#define GENOME_TABLE_H

#include <vector>
#include <string>
#include <map>

using std::vector;
using std::string;


class StringView;

class genome_table{
 public:
  vector< string > contigs;
  std::vector<unsigned int> contig_sizes;
  
  genome_table();
  genome_table(std::string gt_fname);
  bool load(std::string gt_fname); //load the genome table from file
  int contig_ind(const string& contig); //return contig index if exists or -1 otherwise;  
  int contig_size(const string& contig); //returns the size of correspoinding chromosome
  unsigned int size() const;

  bool compare_coordinates(const string& chrom1, const int& pos1, const char& strand1, 
			   const string& chrom2, const int& pos2, const char& strand2);
  bool compare_coordinates(const StringView& chrom1, const int& pos1, const char& strand1, 
			   const StringView& chrom2, const int& pos2, const char& strand2);
  bool compare_coordinates(const StringView& chrom1, const int& pos1, 
			   const StringView& chrom2, const int& pos2);
  //returns true if (chrom1,pos1) < (chrom2,pos2), false otherwise

 private:
  std::map<std::string, int> contig_map; //allows to quickly search for contigs
};


struct chrom_info{
  std::string chrom;
  int size;
  int offset; //global offset
  int buckets; //buckets in this chrom
};


class genome_hash_table{
  //this class defines a hash table that allows for the fast
  //coordinate search 

 private: 
  std::vector<chrom_info> chroms;

  int seed_size; //this is the coarsess of the hash table.
  //good size would be 10^6

 public:
  int size(); //returns the number of hash values
  genome_hash_table(genome_table gt, int _seed_size);
  int hash_value(std::string chrom, int coord);
  int dist_to_left_boundary(std::string chrom, int coord);
  int dist_to_right_boundary(std::string chrom, int coord);
};

#endif
