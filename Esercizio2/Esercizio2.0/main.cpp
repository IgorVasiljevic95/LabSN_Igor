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
	
	//Define and Initialize the random seed and primes rnd.initializeRandom
	Random rnd;
	rnd.initializeRandom("Primes","seed.in");
	
	// Total number random values generated M
	const int M=10000;
	// Values of boxes for each case N
	const int N=100;
	// Number of throws in each block L
	int L=0.;
	L=int(M/N);

	//Inizialize armadillo's vector
	vec avgUnif = zeros<vec>(N);
	vec errUnif = zeros<vec>(N);
	vec avgPx = zeros<vec>(N);
	vec errPx = zeros<vec>(N);
	vec randomgen = zeros<vec>(M);
	vec randomgenpx = zeros<vec>(M);
   
	//Generate M random numbers with uniform distribution and with non-distribution
	GenerateRandomNumber(M,randomgen,rnd);
	GenerateRandomNumberPx(M, randomgenpx, rnd);
	//Data blocking for mean and uncertainty of the integral with uniform distributio
	DataBlocking(N,L,randomgen, avgUnif, errUnif);
	DataBlockingPx(N,L,randomgenpx,avgPx, errPx);

	//Save results into files
	WriteToFile(avgUnif, errUnif,avgPx, errPx, "datiIngral.dat");

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
