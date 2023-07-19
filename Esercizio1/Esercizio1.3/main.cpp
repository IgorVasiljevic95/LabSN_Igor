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
   
	Random rnd;
	//Initialize the random seed and primes rnd.initializeRandom
	rnd.initializeRandom("Primes","seed.in");
	
	//Total number random values generated M
	const int M=100000;
	//Values of boxes for each case
	const int N=100;
  //Number of throws in each block
	int L=0.;
	L=int(M/N);
	//Length of space between lines
	double d=1.;
	//Length of the throwing stick
	double L1=0.8;

	//Inizialize armadillo's vector
	vec mean = zeros<vec>(N);
	vec err = zeros<vec>(N);
	vec randomtheta = zeros<vec>(M);
	vec randompos = zeros<vec>(M);
	
	//Generate M random angle theta and random position
	GenerateRandomNumberTheta(M,randomtheta,rnd);
	GenerateRandomNumberTimes(M,randompos,d,rnd);
	
	//Data blocking for mean and uncertainty of pi
	DataBlocking(N,L,L1,d,randompos,randomtheta,mean,err);

	//Save results into files
	WriteToFile(mean,err, "dati.dat");
   
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
