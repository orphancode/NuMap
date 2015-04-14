#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "assert.h"
#include "math_functions.h"
#include "kernels.h"


using std::vector;
using std::string;
using std::cerr;
using std::cout;
using std::endl;

discrete_kernel::discrete_kernel(double _bw, string _type){
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

double& discrete_kernel::operator[](int i){
  int offset = i + zero_ind;
  assert(offset >= 0 && offset < (int) values.size());
  return values[offset];
}

double discrete_kernel::at(int ind) const{
  int offset = ind + zero_ind;

  if(offset < 0 || offset >= (int) values.size()){
    return 0;
  } else {
    return values[offset];
  }
}
double discrete_kernel::sum(){
  double res = 0;
  for(unsigned int i=0; i<values.size(); i++){
    res += values[i];
  }
  return res;
}


