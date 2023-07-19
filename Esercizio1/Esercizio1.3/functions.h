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
//Generate M random number
void GenerateRandomNumber(int M, vec& values, Random rnd);
//Generate M random angle without the value pi
void GenerateRandomNumberTheta(int M, vec& values, Random rnd);
//Generate M random d times number
void GenerateRandomNumberTimes(int M, vec& values, double d, Random rnd);
// This function returns the statistical uncertainty and the mean of pi given number of blocks N,the number of throws in each block L, the length of a stick L1, the leight of gap lines and two armadillo's vector v1 (random position) and v2 (random angle) with data blocking method
void DataBlocking(int N, int L,double L1, double d, const vec& v,const vec& v2, vec& mean, vec& err);
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
