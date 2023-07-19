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
   
	ofstream uscita;
	//Path for the results of the coding for Python
	uscita.open("dati.dat");
	Random rnd;
	//Initialize the random seed and primes rnd.initializeRandom
	rnd.initializeRandom("Primes","seed.in");
	
	// Total number random values generated M
	const int M=10000;
	// Values of boxes for each case N
	const int N=100;
	// Number of throws in each block L
	int L=0.;
	L=int(M/N);
	
	// Cubic lattice with lattice constant a
	const double a = 1.0;

	// Inizialize armadillo's vector
	vec avg = zeros<vec>(N);
	vec err = zeros<vec>(N);
	vec avgCont = zeros<vec>(N);
	vec errCont = zeros<vec>(N);
	vec randomgen=zeros<vec>(M);
	vec randomgen2=zeros<vec>(M);
	vec randomgentheta=zeros<vec>(M);
	vec randomgenphi=zeros<vec>(M);
	
	// First fill randomgen with M random direction x=1, y=2 and z=3, then fill randomgen2 with -1 or +1 (backward or forward), then fill randomgentheta with random angle and fill rangomgenphi with random polar angle
	for (unsigned int i=0;i<M;i++) {
		randomgen(i)=rnd.RannyuInt(1,3);
		randomgen2(i)=rnd.RannyuNegInt();
		randomgentheta(i)=2*M_PI*rnd.Rannyu();
		randomgenphi(i)=M_PI*rnd.Rannyu();
	}
	
	// Data blocking for mean and uncertainty of the radius with discret walk
	pair<vec,vec> resultSum=ArraySumErr(N,L,a,randomgen,randomgen2);
	err = resultSum.first;
	avg = resultSum.second;
	// Data blocking for mean and uncertainty of the radius with continuum walk
	pair<vec,vec> resultSumCont=ArraySumErrConti(N,L,a,randomgentheta,randomgenphi);
	errCont = resultSumCont.first;
	avgCont = resultSumCont.second;
		
	// Save results into files
	uscita << "val" << "\t" << "err" << "\t" << "valCont" << "\t" << "errCont" <<endl;

	for(unsigned int i=0; i<N;i++) {
		uscita << avg(i) << "\t" << err(i) << "\t" << avgCont(i) << "\t" << errCont(i) <<endl;
	}
      
	uscita.close();
   
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
