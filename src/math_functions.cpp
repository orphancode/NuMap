#include <vector>
#include <algorithm>
#include "assert.h"
#include "math.h"

#include "math_functions.h"

using std::vector;

double sqr(double x){
  return x*x;
}

int max(int i, int j){
  if(i>j) return i;
  return j;
}

unsigned int max(unsigned int i, unsigned int j){
  if(i>j) return i;
  return j;
}

int min(int i, int j){
  if(i<j) return i;
  return j;
}

unsigned int min(unsigned int i, unsigned int j){
  if(i<j) return i;
  return j;
}

int max(const vector<int>& v){
  assert(v.size() > 0);
  int res = v[0];
  for(unsigned int i=1; i<v.size(); i++){
    if(res < v[i]){
      res = v[i];
    }
  }
  return res;
}

int sum(const std::vector<int>& v){
  int res = 0;
  for(unsigned int i=0; i<v.size(); i++){
    res += v[i];
  }
  return res;
}

double mean(const vector<double>& v){

  double res = 0;
  double _size = (double) v.size();

  for(unsigned int i=0; i<v.size(); i++){
    res += v[i] / _size;
  }
  return res;
}


double mean(const vector<int>& v){
  if(v.size() == 1){ return (double)v[0]; }

  double res = 0;
  double _size = (double) v.size();

  for(unsigned int i=0; i<v.size(); i++){
    res += ((double)v[i]) / _size;
  }
  return res;
}

double median(const vector<double>& v){
  assert(v.size() > 0);
  double res = 0;

  vector<double> copy_vec = v;
  sort(copy_vec.begin(), copy_vec.end());
  
  if(copy_vec.size() % 2 == 1){
    int mode_ind = (copy_vec.size()-1)/2;
    res = copy_vec[mode_ind];
  } else {
    res = 0.5 * (copy_vec[copy_vec.size()/2 - 1] + copy_vec[copy_vec.size()/2]);
  }
  
  return res;
}

double median(const vector<int>& v){
  assert(v.size() > 0);
  double res = 0;

  vector<int> copy_vec = v;
  sort(copy_vec.begin(), copy_vec.end());
  
  if(copy_vec.size() % 2 == 1){
    int mode_ind = (copy_vec.size()-1)/2;
    res = (double)copy_vec[mode_ind];
  } else {
    res = 0.5 * (((double)copy_vec[copy_vec.size()/2 - 1]) + ((double)copy_vec[copy_vec.size()/2]));
  }
  
  return res;
}

double stdev(const vector<double>& v){
  assert(v.size() > 0);
  double v_mean = mean(v);
  double var = 0;
  double _size = (double) v.size();
  for(unsigned int i=0; i<v.size(); i++){
    
    var += (v[i]-v_mean)*(v[i]-v_mean) / _size;
  }
  return sqrt(var);
}
double stdev(const vector<int>& v){
  assert(v.size() > 0);
  double v_mean = mean(v);
  double var = 0;
  double _size = (double) v.size();
  for(unsigned int i=0; i<v.size(); i++){
    
    var += ((double)v[i]-v_mean)*((double)v[i]-v_mean) / _size;
  }
  return sqrt(var);
}

double quantile(const vector<int>& v, double _quantile){
  assert(_quantile >= 0 && _quantile <= 100);
  assert(v.size() > 0);
  if(v.size() == 1){
    return v[0];
  }
  if(_quantile == 0.0) return v[0];
  if(_quantile == 100.0) return v[v.size()-1];
  
  double res = 0;

  vector<int> copy_vec = v;
  sort(copy_vec.begin(), copy_vec.end());
  
  int left_ind = floor((_quantile/100.0) * ((double)v.size()));
  int right_ind = left_ind + 1;

  if(right_ind >= (int)v.size()){
    assert(left_ind >= 0);
    return copy_vec[left_ind];
  }

  assert(left_ind >=0);
  assert(right_ind < (int)v.size());

  double remainder_frac = (_quantile / 100.0) * ((double) v.size()) - ((double)left_ind);
  assert(remainder_frac >=0 && remainder_frac <= 1);
  
  res = copy_vec[left_ind] + remainder_frac * ( ((double)copy_vec[right_ind]) - ((double)copy_vec[left_ind]) );
  return res;
}

double quantile(const vector<double>& v, double _quantile){
  assert(_quantile >= 0 && _quantile <= 100);
  assert(v.size() > 0);
  if(v.size() == 1){
    return v[0];
  }
  if(_quantile == 0.0) return v[0];
  if(_quantile == 100.0) return v[v.size()-1];
  
  double res = 0;

  vector<double> copy_vec = v;
  sort(copy_vec.begin(), copy_vec.end());
  
  int left_ind = floor((_quantile/100.0) * ((double)v.size()));
  int right_ind = left_ind + 1;

  if(right_ind >= (int)v.size()){
    assert(left_ind >= 0);
    return copy_vec[left_ind];
  }

  assert(left_ind >=0);
  assert(right_ind < (int)v.size());

  double remainder_frac = (_quantile / 100.0) * ((double) v.size()) - ((double)left_ind);
  assert(remainder_frac >=0 && remainder_frac <= 1);
  
  res = copy_vec[left_ind] + remainder_frac * ( ((double)copy_vec[right_ind]) - ((double)copy_vec[left_ind]) );
  return res;
}

/*
double covariance(const vector<double>& v1, const vector<double>& v2){
  assert(v1.size()  == v2.size());

  double v1_mean = mean(v1);
  double v2_mean = mean(v2);
  
  double _size = (double) v1.size();
  double res = 0;
  for(unsigned int i=0; i<v1.size(); i++){
    res += (v1[i] - v1_mean) * (v2[i] - v2_mean) / _size; 
  }
  return res;
}
*/


double covariance(const vector<double>& v1, const vector<double>& v2){
  assert(v1.size() == v2.size());
  double v1_mean = mean(v1);
  double v2_mean = mean(v2);

  vector<double> prod(v1.size(), 0);
  for(unsigned int i=0; i<prod.size(); i++){
    prod[i] = v1[i]*v2[i];    
  }

  double prod_mean = mean(prod);

  double res = prod_mean - v1_mean*v2_mean;
  return res;
}

double correlation(const vector<double>& v1, const vector<double>& v2){
  assert(v1.size() > 0);
  assert(v1.size() == v2.size());

  double v1_mean = mean(v1);
  double v2_mean = mean(v2);

  double v1_var = 0;
  double v2_var = 0;

  double _size = (double) v1.size();

  for(unsigned int i=0; i<v1.size(); i++){
    v1_var += (v1[i] - v1_mean) * (v1[i] - v1_mean) / _size; 
    v2_var += (v2[i] - v2_mean) * (v2[i] - v2_mean) / _size; 
  }

  double v1_stdev = sqrt(v1_var);
  double v2_stdev = sqrt(v2_var);

  double cor = 0;

  double den = _size * v1_stdev * v2_stdev;

  for(unsigned int i=0; i<v2.size(); i++){
    cor += (v1[i] - v1_mean) * (v2[i] - v2_mean) / den;
  }

  return cor;
}

double cube(double x){
  return x*x*x;
}

double triweight_kernel(double x, double bw){
  if(x<-bw || x>bw) return 0;
  return 1.09375 * cube(1-(x*x)/(bw*bw)) / bw;
}

double Epanechnikov_kernel(double x, double bw){

  if(x<-bw || x>bw) return 0;
  double u = x/bw;
  return 0.75 * (1-u*u) / bw;
}

