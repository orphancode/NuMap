// This class describes a kernel density object

#ifndef KERNELS_H
#define KERNELS_H

using std::vector;

class discrete_kernel{
 public:
  int bw;
  std::string type;
  int zero_ind; //index of zero in the values profile
  int limit; //kernel support [-limit, limit]

  vector<double> values;

  discrete_kernel(double _bw, std::string _type);
  double at(int ind) const; //returns precomputed value at ind
                                                                                            
  double& operator[](int i);

  double sum(); //sanity check to make sure kernel adds up to 1 
};


#endif
