#ifndef SAM_ALIGNMENT_H
#define SAM_ALIGNMENT_H

#include "string_utils.h"

class StringViewCigar{
 public:
  std::vector<int> numbers;
  std::vector<char> chars;

  StringViewCigar(){};

  StringViewCigar(const StringView& strView);
  int init(const StringView& strView); //-1 for error, 0 for ok
  int reference_length(); //the size of the reference region spanned by left and right end of the mapped part  
  int mapped_length();  //the number of base pairs mapped to the reference                                     
  int get_cigar_element_length(int ind); //returns the reference length for the cigar element with index ind   
  StringViewCigar& operator=(const StringViewCigar& rhs){
    if(this == &rhs) return *this;
    numbers = rhs.numbers;
    chars = rhs.chars;
    return *this;
  }
};

class mdfield{
 public:
  std::vector<int> numbers;
  std::vector<char> ops; //operations: I - insertion, X - mismatch (i.e. SNP)
  std::vector<std::string> bases; //bases

  mdfield(){};
  //mdfield(const std::string& inp_string);
  mdfield(const StringView& inp_string);
  void init(const StringView& inp_string);
};

// this file provides definitions for the SAM alignment class

class SamAlignment{//: public Alignment{ 
 public:
  
  std::string sam_string;

  unsigned int flag; //described in the sam file document
  StringViewCigar _cigar;
  mdfield* ptr_mdfield;

  int begin; 
  int end;
  int mapq;

  StringView chrom;
  char strand;

  StringView qname;
  StringView seq;
  StringView qual;
  StringView cigar_string;
  StringView mdfield_string;

  bool _valid;
  bool _mapped;

 public:
  SamAlignment(std::string& _sam_string, std::string& error_msg);
  SamAlignment();
  SamAlignment(const SamAlignment& al);
  ~SamAlignment();
  SamAlignment& operator=(const SamAlignment& rhs); 

  int init_cigar(); //fills in values of _cigar
  void init_mdfield(); //fills in values of _mdfield

  bool output(std::ofstream& ofstr); //output sam string into this stream;
  bool mapped();
  bool valid();
};

/*
class SamMPAlignment{
public:
	unsigned flag;
	std::string cigar_string;
	std::string sam_string;
	
	std::string qname;
	int mapq;
	std::string qual;
	std::string seq;	

	//std::string mdfield;

	bool mapped1;
	bool mapped2;
		
	char strand1; 
	char strand2;
	
	std::string chrom1;
	std::string chrom2;

	int begin1;
	int begin2;

	SamMPAlignment(std::string& _sam_string);
	SamMPAlignment();
	~SamMPAlignment(){};

	//bool mapped();
	bool valid();
	
private:
	bool _valid; //valid mp entry
	//bool _mapped;
};
*/
#endif
