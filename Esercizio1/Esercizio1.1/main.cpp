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
	
	// Total number random values generated M
	const int M=100000;
	// Number of boxes N
	const int N=100;
	// Number of throws in each block L
	int L=0.;
	L=int(M/N);
	
	// Inizialize armadillo's vector
	vec randomgen = zeros<vec>(M);
	vec mean = zeros<vec>(N);
	vec err = zeros<vec>(N);
   
	//Generate M random number
	GenerateRandomNumber(M,randomgen,rnd);
	DataBlocking(N,L,randomgen,mean,err);
	
	//Save results into file
	WriteToFile(mean,err,"dati.dat");
		
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
