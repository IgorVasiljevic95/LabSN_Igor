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

vec ArraySum(int N,int L, const vec& v);
// This function gives an armadillo's array filled with mean of uniform distribution in each blocks
vec ArraySumUnif(int N, int L, Random rnd);
// This function gives an armadillo's array filled with mean of exponential distribution in each blocks
vec ArraySumExp(int N, int L, double lambda, Random rnd);
// This function gives an armadillo's array filled with mean of cauchy distribution in each blocks
vec ArraySumCau(int N, int L, double mean, double delta, Random rnd);
//Generate M random number between [0,1)
void GenerateRandomNumber(int M, vec& values, Random rnd);
// This function returns the statistical uncertainty and the mean given number of blocks N and the number of throws in each block L and a armadillo's vector v with data blocking method
void DataBlocking(int N, int L, const vec& v, vec& mean, vec& err);
//Function to write results on a file named filename
void WriteToFile(const vec& mean, const vec& err, const string& filename);

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
