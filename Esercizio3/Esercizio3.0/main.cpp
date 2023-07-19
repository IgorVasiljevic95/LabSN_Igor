/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Igor Vasiljevic
_/    _/  _/_/_/  _/_/_/_/ email: igor.vasiljevic@studenti.unimi.it
*****************************************************************
*****************************************************************/
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <iomanip>
#include <armadillo>
#include "functions.h"

using namespace std;
using namespace arma;

//This program is to calcolate via Monte Carlo the European call-option price and put-option price, with two stregy, sampling directly the final assect price and sampling discretized path of the asset price.

int main (int argc, char *argv[]){
   
	//Define and Initialize the random seed and primes rnd.initializeRandom
	Random rnd;
	rnd.initializeRandom("Primes","seed.in");
	
	//Total number random values generated M
	const int M=100000;
	//Number of boxes N
	const int N=100;
	//Number of interval for divide the path of asset price
	const int Nt=100;
	//Initial values S(0) of the asset price
	double S0=100.;
	//Delivery time T=1.
	double T=1.;
	//Strike price:  𝐾=100
	double K=100.;
	//Risk-free interest rate
	double r=0.1;
	//Volatility
	double sigma=0.25;
	//Number of throw in each block
	int L=0.;
	L=int(M/N);
   
	//Inizialize armadillo's vector
	vec avgCallDir = zeros<vec>(N);
	vec avgPullDir = zeros<vec>(N);
	vec errCallDir = zeros<vec>(N);
	vec errPullDir = zeros<vec>(N);
	vec avgCallDis = zeros<vec>(N);
	vec avgPullDis = zeros<vec>(N);
	vec errCallDis = zeros<vec>(N);
	vec errPullDis = zeros<vec>(N);
	vec C_Dir = zeros<vec>(M);
	vec C_Dis = zeros<vec>(M);
	vec P_Dir=zeros<vec>(M);
	vec P_Dis=zeros<vec>(M);
	vec price_direct = zeros<vec>(M);
	vec price_discret = zeros<vec>(M);
	
	
	//Main loop for call/out-option price
	//Direct and Discret simulation of call-option price
	pair<vec,vec> Csum=sim_call(M,Nt,S0,K,r,sigma,T,rnd);
	C_Dir = Csum.first;
	C_Dis = Csum.second;
	//Direct and Discret simulation of out-option price
	pair<vec,vec> Psum=sim_pull(M,Nt,S0,K,r,sigma,T,rnd);
	P_Dir = Psum.first;
	P_Dis = Psum.second;
	
	
	//Data blocking call-option price
	pair<vec,vec> resultSumC=ArraySumErr(N,L,C_Dir);
	errCallDir = resultSumC.first;
	avgCallDir = resultSumC.second;
	pair<vec,vec> resultSumCs=ArraySumErr(N,L,C_Dis);
	errCallDis = resultSumCs.first;
	avgCallDis = resultSumCs.second;
	//Data blocking out-option price
	pair<vec,vec> resultSumP=ArraySumErr(N,L,P_Dir);
	errPullDir = resultSumP.first;
	avgPullDir = resultSumP.second;
	pair<vec,vec> resultSumPs=ArraySumErr(N,L,P_Dis);
	errPullDis = resultSumPs.first;
	avgPullDis = resultSumPs.second;
   
	//Distribution of the price
	price_direct=simu_direct(M,S0,K,r,sigma,T,rnd);
	price_discret=simu_discret(M,Nt,S0,K,r,sigma,T,rnd);
	
	//Save results into files
	WriteToFile(avgCallDir, errCallDir,avgCallDis,errCallDis, avgPullDir, errPullDir,avgPullDis,errPullDis,"dati.dat");
	WriteToFilePrice(price_direct, price_discret,"datiprice.dat");

	rnd.SaveSeed();
	return 0;
   
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
