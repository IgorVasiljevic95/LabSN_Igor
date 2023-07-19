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
 
int main (int argc, char *argv[]){
   	
	//Define and initialize the random seed and primes rnd.initializeRandom
	Random rnd;
	rnd.initializeRandom("Primes","seed.in");
	
	//Number of simulation
	const int M=1000000;
	//Number of boxes N
	const int N=100;
	//Number of throws in each block L
	int L=0.;
	L=int(M/N);
	
	//Starting point for Psi100, not fare from the ipotetic values
	double x0=1.5;
	double y0=1.5;
	double z0=1.5;
	//Starting point for Psi210, not fare from the ipotetic values (do not use (0,0,0)
	double x=5.;
	double y=5.;
	double z=5.;
	//Step for random walk for Psi100 to make metropolis acceptance ratio to be approximately 50%
	double k1=1.;
	//Step for random walk for Psi210 to make metropolis acceptance ratio to be approximately 50%
	double k2=2.;
	//Step for Gauss distribution for Psi100 to make metropolis acceptance ratio to be approximately 50%
	double sigma=0.7;
	//Step for Gauss distribution for Psi210 to make metropolis acceptance ratio to be approximately 50%
	double sigma210=1.3;
	
	//Initialization of armadillo's vector
	vec avgr100 = zeros<vec>(M);
	vec avgr210 = zeros<vec>(M);
	vec avgr100G = zeros<vec>(M);
	vec avgr210G = zeros<vec>(M);
	vec err100G = zeros<vec>(N);
  vec err210G = zeros<vec>(N);
  vec avg100G = zeros<vec>(N);
  vec avg210G = zeros<vec>(N);
	vec err100 = zeros<vec>(N);
	vec err210 = zeros<vec>(N);
	vec avg100 = zeros<vec>(N);
	vec avg210 = zeros<vec>(N);
	
	//Metropolis algorithm with approssimatly 50% acc/rej
	avgr100=metropolis100RW(M,x0,y0,z0,k1,rnd);
	avgr210=metropolis210RW(M,x,y,z,k2,rnd);
	avgr100G=metropolis100Gauss(M,x0,y0,z0,sigma,rnd);
	avgr210G=metropolis210Gauss(M,x,y,z,sigma210,rnd);
	
	//Data blocking
	pair<vec,vec> resultSumC=ArraySumErr(N,L,avgr100);
	err100 = resultSumC.first;
	avg100 = resultSumC.second;
	pair<vec,vec> resultSumCs=ArraySumErr(N,L,avgr210);
	err210 = resultSumCs.first;
	avg210 = resultSumCs.second;
	pair<vec,vec> resultSumCG=ArraySumErr(N,L,avgr100G);
	err100G = resultSumCG.first;
	avg100G = resultSumCG.second;
	pair<vec,vec> resultSumCsG=ArraySumErr(N,L,avgr210G);
	err210G = resultSumCsG.first;
	avg210G= resultSumCsG.second;

	//Printing the data for each simulation
	//Function to write results on a file named filename
	WriteToFile(avg100, err100, avg210,err210,"datiRW.dat");
	WriteToFile(avg100G,err100G,avg210G,err210G,"datiGauss.dat");

	//Initialization of armadillo's matrix which will cointain the coordinates of the point of all metropolis simulation
	mat matrix100RW(M, 3, fill::zeros);
	mat matrix210RW(M, 3, fill::zeros);
	mat matrix100Gauss(M, 3, fill::zeros);
	mat matrix210Gauss(M, 3, fill::zeros);

	//Saving coordinates inside a matrix
	matrix100RW=Coord100RW(M,x0,y0,z0,k1,rnd);
	matrix210RW=Coord210RW(M,x,y,z,k2,rnd);
	matrix100Gauss=Coord100Gauss(M,x0,y0,z0,sigma,rnd);
	matrix210Gauss=Coord210Gauss(M,x,y,z,sigma210,rnd);
	
	//Printing in a file the coordinates
	WriteMatrixToFile(matrix100RW, M, 3,"datiPsi100RW.dat");
	WriteMatrixToFile(matrix210RW, M, 3,"datiPsi210RW.dat");
	WriteMatrixToFile(matrix100Gauss, M, 3,"datiPsi100Gauss.dat");
	WriteMatrixToFile(matrix210Gauss, M, 3,"datiPsi210Gauss.dat");
	
	WriteInstaToFile(avgr100,avgr210,"datiInsta.dat");
	
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
