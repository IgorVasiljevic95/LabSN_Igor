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
// This function return the value of the function that is integrate
double Cosin(double x);
// This function return a random number [0,1) using a non uniform probability
double p(double x);
//Generate M random number with uniform distribution
void GenerateRandomNumber(int M, vec& values, Random rnd);
//Generate M random number with non-uniform distribution
void GenerateRandomNumberPx(int M, vec& values, Random rnd);
//This function return the mean and the statistical uncertainty for the integral using a uniform distribution
void DataBlocking(int N, int L, const vec& v, vec& mean, vec& err);
//This function return the mean and the statistical uncertainty for the integral using a non-uniform distribution
void DataBlockingPx(int N, int L, const vec& v, vec& mean, vec& err);
//Function that write data's inside a file
void WriteToFile(const vec& mean, const vec& err,const vec& meanPx, const vec& errPx, const string& filename);

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
