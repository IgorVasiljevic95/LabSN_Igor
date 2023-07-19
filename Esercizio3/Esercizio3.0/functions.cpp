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

double N(double x) {
	 return 0.5*(1.+erf(x/sqrt(2.)));
}

// This function returns the statistical uncertainty and the mean given number of blocks N and the number of throws in each block L and a armadillo's vector v with data blocking method
pair<vec,vec> ArraySumErr(int N, int L, const vec& v){
	 vec avg = zeros<vec>(N);
	 vec avg2 = zeros<vec>(N);
	 vec sum_prog = zeros<vec>(N);
	 vec sum2_prog = zeros<vec>(N);
	 vec err = zeros<vec>(N);
	 double sum1=0.;
	 int k=0;
	 for(unsigned int i=0; i<N;i++) {
			sum1=0.;
			for(unsigned int j=0;j<L;j++) {
				 k = j+i*L;
				 sum1 += v(k);
			}
			avg(i)=sum1/L;
			avg2(i)=avg(i)*avg(i);
	 }
	 for(unsigned int i=0; i<N;i++) {
			for(unsigned int j=0;j<i+1;j++) {
				 sum_prog(i) += avg(j);
				 sum2_prog(i) += avg2(j);
			}
			sum_prog(i) /= (i+1);
			sum2_prog(i) /= (i+1);
			err(i)=error(sum_prog(i), sum2_prog(i),i);
	 }
	 return make_pair(err, sum_prog);
}

vec simu_direct(int N, double S0, double K, double r, double sigma, double T,Random rnd) {
   vec ST=zeros<vec>(N);
   double v=0.;
   double STi=0.;
   for (unsigned int i = 0; i < N; i++) {
      v = rnd.Gauss(0.,1.);
      STi = S0 * exp((r - 0.5*sigma*sigma)*T + sigma*sqrt(T)*v);
      ST(i) = STi;
   }
   return ST;
}

vec simu_discret(int N, int L, double S0, double K, double r, double sigma, double T,Random rnd) {
   vec ST=zeros<vec>(N);
   double delta_t=0.;
   double v=0.;
   double Si=0.;
   delta_t = T/L;
   for (unsigned int i = 0; i < N; i++) {
      Si = S0;
      for (unsigned int j = 0; j < L; j++) {
         v = rnd.Gauss(0.,1.);
         Si = Si*exp((r - 0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*v);
      }
      ST(i) = Si;
   }
   return ST;
}

pair<vec,vec> sim_call(int N, int L, double S0, double K, double r, double sigma, double T,Random rnd) {
   vec C=zeros<vec>(N);
   vec Cd=zeros<vec>(N);
   vec St=zeros<vec>(N);
   vec Std=zeros<vec>(N);
   St=simu_direct(N,S0,K,r,sigma,T,rnd);
   Std=simu_discret(N,L,S0,K,r,sigma,T,rnd);
   for (unsigned int i = 0; i < N; i++) {
      C(i)=exp(-r*T)*max(0.,St(i)-K);
      Cd(i)=exp(-r*T)*max(0.,Std(i)-K);
   }
   return make_pair(C, Cd);
}

pair<vec,vec> sim_pull(int N, int L, double S0, double K, double r, double sigma, double T,Random rnd) {
   vec P=zeros<vec>(N);
   vec Pd=zeros<vec>(N);
   vec St=zeros<vec>(N);
   vec Std=zeros<vec>(N);
   St=simu_direct(N,S0,K,r,sigma,T,rnd);
   Std=simu_discret(N,L,S0,K,r,sigma,T,rnd);
   for (unsigned int i = 0; i < N; i++) {
      P(i)=exp(-r*T)*max(0.,K-St(i));
      Pd(i)=exp(-r*T)*max(0.,K-Std(i));
   }
   return make_pair(P, Pd);
}

//Function to write results on a file named filename
void WriteToFile(const vec& avgCallDir, const vec& errCallDir,const vec& avgCallDis, const vec& errCallDis,const vec& avgPullDir, const vec& errPullDir,const vec& avgPullDis, const vec& errPullDis,const string& filename) {
		ofstream file(filename);
		if (file.is_open()) {
			file << "C_Direct" << "\t" << "err_C_Direct" << "\t" << "C_Disc" << "\t" << "err_C_Disc"<<"\t" << "P_Direct" << "\t" << "err_P_Direct" << "\t" << "P_Disc" << "\t" << "err_P_Disc" <<endl;
			for(unsigned int i=0; i<avgCallDir.size(); i++) {
				file << avgCallDir(i) << "\t" << errCallDir(i) << "\t" << avgCallDis(i)<< "\t" << errCallDis(i)<<"\t" << avgPullDir(i) << "\t" << errPullDir(i) << "\t" << avgPullDis(i)<< "\t" << errPullDis(i) <<endl;
				}
				file.close();
		}
}

//Function to write prices on a file named filename
void WriteToFilePrice(const vec& price_direct, const vec& price_discret,const string& filename) {
		ofstream file(filename);
		if (file.is_open()) {
			file << "St_Direct" << "\t" << "St_Discret" <<endl;
			for(unsigned int i=0; i<price_direct.size(); i++) {
				file << price_direct(i) << "\t" << price_discret(i) <<endl;
			}
			file.close();
		}
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
