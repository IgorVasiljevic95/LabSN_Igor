/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Igor Vasiljevic
_/    _/  _/_/_/  _/_/_/_/ email: igor.vasiljevic@studenti.unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "functions.h"
#include "random.h"
#include <armadillo>

using namespace std;
using namespace arma;

// This function calculate statistical uncertainty estimation
double error(double av, double av2, int n) {
   if (n==0) {
      return 0;
   }
   else {
      return sqrt((av2-av*av)/n);
   }
}

// This function gives an armadillo's array filled with mean in each blocks
vec ArraySum(int N,int L, const vec& v){
   vec avg = zeros<vec>(N);
   double sum=0.;
   for(unsigned int i=0; i<N;i++) {
      sum=0.;
      for(unsigned int j=0; j<L;j++) {
         sum+=v(j);
      }
      avg(i)=sum/L;
   }
   return avg;
}

// This function returns the statistical uncertainty and the mean given number of blocks N and the number of throws in each block L and a armadillo's vector v with data blocking method
pair<vec,vec> ArraySumErr(int N, int L,double a, const vec& v1, const vec& v2){
   vec avg = zeros<vec>(N);
   vec avg2 = zeros<vec>(N);
   vec sum_prog = zeros<vec>(N);
   vec sum2_prog = zeros<vec>(N);
   vec err = zeros<vec>(N);
   double x=0.;
   double y=0.;
   double z=0.;
   int k=0;
   for(unsigned int i=0; i<N;i++) {
      x=0.;
      y=0.;
      z=0.;
      for(unsigned int j=0;j<L;j++) {
         k = j+i*L;
         if(v1(k)==1) {
            x+=v2(k)*a;
         }
         else if (v1(k)==2) {
            y+=v2(k)*a;
         }
         else {
            z+=v2(k)*a;
         }
      }
      avg(i)=abs(x*x+y*y+z*z);
      avg2(i)=avg(i)*avg(i);
   }
   for(unsigned int i=0; i<N;i++) {
      for(unsigned int j=0;j<i+1;j++) {
         sum_prog(i) += avg(j);
         sum2_prog(i) += avg2(j);
      }
      sum_prog(i) /= (i+1);
      sum2_prog(i) /= (i+1);
      sum_prog(i)=sqrt(sum_prog(i));
      sum2_prog(i)=sqrt(sum2_prog(i));
      err(i)=error(sum_prog(i), sum2_prog(i),i);
   }
   return make_pair(err, sum_prog);
}

pair<vec,vec> ArraySumErrConti(int N, int L,double a, const vec& v1, const vec& v2){
   vec avg = zeros<vec>(N);
   vec avg2 = zeros<vec>(N);
   vec sum_prog = zeros<vec>(N);
   vec sum2_prog = zeros<vec>(N);
   vec err = zeros<vec>(N);
   double x=0.;
   double y=0.;
   double z=0.;
   int k=0;
   for(unsigned int i=0; i<N;i++) {
      x=0.;
      y=0.;
      z=0.;
      for(unsigned int j=0;j<L;j++) {
         k = j+i*L;
         x+=a*cos(v1(k))*sin(v2(k));
         y+=a*sin(v1(k))*sin(v2(k));
         z+=a*cos(v2(k));
         
      }
      avg(i)=abs(x*x+y*y+z*z);
      avg2(i)=avg(i)*avg(i);
   }
   for(unsigned int i=0; i<N;i++) {
      for(unsigned int j=0;j<i+1;j++) {
         sum_prog(i) += avg(j);
         sum2_prog(i) += avg2(j);
      }
      sum_prog(i) /= (i+1);
      sum2_prog(i) /= (i+1);
      sum_prog(i)=sqrt(sum_prog(i));
      sum2_prog(i)=sqrt(sum2_prog(i));
      err(i)=error(sum_prog(i), sum2_prog(i),i);
   }
   return make_pair(err, sum_prog);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Igor Vasiljevic
_/    _/  _/_/_/  _/_/_/_/ email: igor.vasiljevic@studenti.unimi.it
*****************************************************************
*****************************************************************/
