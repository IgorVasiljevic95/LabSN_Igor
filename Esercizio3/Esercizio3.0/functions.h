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
double N(double);
// This function gives an armadillo's array filled with mean in each blocks
vec ArraySum(int, int,const vec& v);
// This function return the statistical uncertainty and the mean given number of blocks N and the number of throws in each block L and a armadillo's vector v with data blocking metho
pair<vec,vec> ArraySumErr(int,int, const vec& v);
//Direct and Discret simulation of call-option price
pair<vec,vec> sim_call(int N, int L, double S0, double K, double r, double sigma, double T,Random rnd);
//Direct and Discret simulation of out-option price
pair<vec,vec> sim_pull(int N, int L, double S0, double K, double r, double sigma, double T,Random rnd);

//Give me data of distribution of the price with both discret and direct method
vec simu_direct(int N, double S0, double K, double r, double sigma, double T,Random rnd);
vec simu_discret(int N, int L, double S0, double K, double r, double sigma, double T,Random rnd);
//Function to write results on a file named filename
void WriteToFile(const vec& avgCallDir, const vec& errCallDir,const vec& avgCallDis, const vec& errCallDis,const vec& avgPullDir, const vec& errPullDir,const vec& avgPullDis, const vec& errPullDis,const string& filename);
//Function to write prices on a file named filename
void WriteToFilePrice(const vec& price_direct, const vec& price_discret,const string& filename);


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
