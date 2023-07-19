/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Igor Vasiljevic
_/    _/  _/_/_/  _/_/_/_/ email: igor.vasiljevic@studenti.unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Functions__
#define __Functions__

#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

// This function calculate statistical uncertainty estimation
double error(double, double, int);
// This function gives an armadillo's array filled with mean in each blocks
vec ArraySum(int, int,const vec& v);
// This function return the statistical uncertainty and the mean given number of blocks N and the number of throws in each block L and a armadillo's vector v with data blocking method
pair<vec,vec> ArraySumErr(int,int,double a, const vec& v1, const vec& v2);
pair<vec,vec> ArraySumErrConti(int N, int L,double a, const vec& v1, const vec& v2);

#endif // __Functions__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Igor Vasiljevic
_/    _/  _/_/_/  _/_/_/_/ email: igor.vasiljevic@studenti.unimi.it
*****************************************************************
*****************************************************************/
