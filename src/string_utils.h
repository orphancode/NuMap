#ifndef UTILS_H
#define UTILS_H

class StringView{
 public:
  std::string const *ptrStr;
  int begin;
  int _length;
  StringView(){ ptrStr = NULL; begin = 0; _length = 0; }
  StringView(const std::string &str, int _begin, int __length) { ptrStr = &str; begin=_begin; _length=__length;}
  StringView(const StringView &strView, int _begin, int __length){
    ptrStr = strView.ptrStr; begin = strView.begin + _begin; _length = __length; };

  const char* c_str() const;
  unsigned int size() const{ return (unsigned int) _length; };
  unsigned int length() const{ return (unsigned int) _length; };
  //std::string get_string() const;
  std::string str() const;

  StringView& operator=(const StringView& rhs){
    ptrStr = rhs.ptrStr; begin = rhs.begin; _length = rhs._length;
    return *this;
  }
  //bool operator==(const char* str);
  //bool operator==(const StringView& rhs);
  //bool operator <(StringView& rhs) const; 

  const char operator[](int i) const{ return (*ptrStr)[begin+i]; }
  const char at(int i) const;//{ assert(i<_length); return (*ptrStr)[begin+i]; }
  friend std::ostream &operator<<(std::ostream& stream, const StringView& strView){
    if(strView._length==0) return stream;
    stream<<(strView.ptrStr->substr(strView.begin,strView._length));
    return stream;
  }

  std::string substr(int beg, int len) const{ return ptrStr->substr(begin + beg, len); }
};

std::vector<std::string> split(const std::string& inp_string, const std::string& regexp);
std::vector<std::string> split(const std::string& inp_string, char sep);
void split(const std::string& inp_string, char sep, std::vector<std::string>& res);
void split(std::string& inp_string, char sep, std::vector<StringView>& res);
void split(const StringView& inpStringView, char sep, std::vector<StringView>& res);

void chomp(std::string& str); //remove \n
std::ios::pos_type perc_file_length(std::ifstream& ifs);
unsigned int get_string_count(std::ifstream& ifs);

std::string int2string(int i);
std::string commify(int i); //for 12345678 it turns into more readable 12,345,678
std::string commify_long(long i);

bool operator==(const StringView& lhs, const StringView& rhs);
bool operator==(const StringView& lhs, const char* rhs);
bool operator==(const StringView& lhs, const std::string& rhs);
bool operator<(const StringView& lhs, const StringView& rhs);

#endif
