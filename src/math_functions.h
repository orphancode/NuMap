#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

int max(int i, int j);
unsigned int max(unsigned int i, unsigned int j);

int min(int i, int j);
unsigned int min(unsigned int i, unsigned int j);

double covariance(const std::vector<double> &v1, const std::vector<double> &v2); //covariance of the two vectors
double correlation(const std::vector<double>& v1, const std::vector<double>& v2);

double mean(const std::vector<double>& v);
double mean(const std::vector<int>& v);
double stdev(const std::vector<double>& v);
double stdev(const std::vector<int>& v);
double median(const std::vector<double>& v);
double median(const std::vector<int>& v);

double quantile(const std::vector<double>& v, double _quantile);
double quantile(const std::vector<int>& v, double _quantile);

int max(const std::vector<int>& v);
int sum(const std::vector<int>& v);

// KDE functionality
double triweight_kernel(double x, double bw);
double Epanechnikov_kernel(double x, double bw);

double sqr(double x);

#endif
