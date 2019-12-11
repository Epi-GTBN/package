#include <iostream>
//using namespace std;
#include <algorithm>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* These functions are defined in the entropy.cpp file */

double digamma(double z);
double entropy_empirical(std::map< std::vector<int>,int > frequencies, int nb_samples);
double entropy_dirichlet(std::map<  std::vector<int> ,int > frequencies, int nb_samples, double beta);
double entropy_miller_madow(std::map<  std::vector<int> ,int > frequencies, int nb_samples);
double entropy_shrink(std::map<  std::vector<int> ,int > frequencies, int nb_samples);
double entropy(const int *d, int nsamples, int nvars, int c, bool *v);

/* Entry points called from the R functions */
extern "C" 
{
    SEXP entropyR(SEXP data, SEXP nrows, SEXP ncols, SEXP choice);
}

