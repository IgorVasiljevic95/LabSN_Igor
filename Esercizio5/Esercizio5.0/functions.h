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

//This function calculate statistical uncertainty estimation
double error(double, double, int);
//Ground state wave function psi100 with n=1, l=0, m=0 for hydrogen atom in natural units
double psi100(double x, double y, double z);
//Ground state wave function psi210 with n=2, l=1, m=0  for hydrogen atom in natural units
double psi210(double x, double y, double z);
//Metropolis loop for Psi100 with random walk moves
vec metropolis100RW(int N, double x0,double y0,double z0, double a,Random rnd);
//Metropolis loop for Psi210 with random walk moves
vec metropolis100Gauss(int N, double x0,double y0,double z0, double sigma,Random rnd);
//Metropolis loop for Psi100 with Gaussian moves
vec metropolis210RW(int N, double x0,double y0,double z0, double a,Random rnd);
//Metropolis loop for Psi210 with Gaussian moves
vec metropolis210Gauss(int N, double x0,double y0,double z0, double sigma,Random rnd);
// This function returns the statistical uncertainty and the mean given number of blocks N and the number of throws in each block L and a armadillo's vector v with data blocking method
pair<vec,vec> ArraySumErr(int,int, const vec& v);


//Those 4 function are just for saving the coordination inside a matrix. I should add it inside a main loop of the metropolis algorithm, so the code will be a bit rusty. If you want just to do the metropolis algorithm then comment those function.
mat Coord100RW(int N, double x0,double y0,double z0, double stepsize,Random rnd);
mat Coord210RW(int N, double x0,double y0,double z0, double stepsize,Random rnd);
mat Coord100Gauss(int N, double x0,double y0,double z0, double sigma,Random rnd);
mat Coord210Gauss(int N, double x0,double y0,double z0, double sigma,Random rnd);
//Function to write results on a file named filename
void WriteToFile(const vec& avg100, const vec& err100,const vec& avg210, const vec& err210,const string& filename);
//Function to write results on a file named filename
void WriteMatrixToFile(const mat& matrix, int M, int N,const string& filename);
void WriteInstaToFile(const vec& avgr100, const vec& avgr210,const string& filename);
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
