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
	ofstream uscita1;
	ofstream uscita2;
	
	//Path for the results of the coding for Python
	uscita.open("datiUnif.dat");
	uscita1.open("datiExp.dat");
	uscita2.open("datiCauchy.dat");
   
	Random rnd;
	//Initialize the random seed and primes rnd.initializeRandom
	rnd.initializeRandom("Primes","seed.in");
	
	//Total number random values generated M
	const int M=10000;
	//Values of boxes for each case
	const int N1=1;
	const int N2=2;
	const int N10=10;
	const int N100=100;
	//Values for distributions
	double lambda=1.;
	double delta=1.;
	double mean=0.;

	//Inizialize armadillo's vector
	vec avgUnif1 = zeros<vec>(M);
	vec avgUnif2 = zeros<vec>(M);
	vec avgUnif10 = zeros<vec>(M);
	vec avgUnif100 = zeros<vec>(M);
	vec avgExp1 = zeros<vec>(M);
	vec avgExp2 = zeros<vec>(M);
	vec avgExp10 = zeros<vec>(M);
	vec avgExp100 = zeros<vec>(M);
	vec avgCa1 = zeros<vec>(M);
	vec avgCa2 = zeros<vec>(M);
	vec avgCa10 = zeros<vec>(M);
	vec avgCa100 = zeros<vec>(M);
      
	//Generate M random number for each N boxes with uniform distribution
	avgUnif1=ArraySumUnif(M,N1,rnd);
	avgUnif2=ArraySumUnif(M,N2,rnd);
	avgUnif10=ArraySumUnif(M,N10,rnd);
	avgUnif100=ArraySumUnif(M,N100,rnd);
	
	//Generate M random number for each N boxes with exponential distribution
	avgExp1=ArraySumExp(M,N1,lambda,rnd);
	avgExp2=ArraySumExp(M,N2,lambda,rnd);
	avgExp10=ArraySumExp(M,N10,lambda,rnd);
	avgExp100=ArraySumExp(M,N100,lambda,rnd);
	
	//Generate M random number for each N boxes with cauchy distribution
	avgCa1=ArraySumCau(M,N1,mean,delta,rnd);
	avgCa2=ArraySumCau(M,N2,mean,delta,rnd);
	avgCa10=ArraySumCau(M,N10,mean,delta,rnd);
	avgCa100=ArraySumCau(M,N100,mean,delta,rnd);
	
	//Save results into files
	uscita << "N=1" << "\t" << "N=2" << "\t" << "N=10" << "\t" << "N=100" << endl;
	uscita1 << "N=1" << "\t" << "N=2" << "\t" << "N=10" << "\t" << "N=100" << endl;
	uscita2 << "N=1" << "\t" << "N=2" << "\t" << "N=10" << "\t" << "N=100" << endl;

	for(unsigned int i=0; i<M;i++) {
		uscita << avgUnif1(i) << "\t" << avgUnif2(i) << "\t" << avgUnif10(i) << "\t" << avgUnif100(i) << endl;
		uscita1 << avgExp1(i) << "\t" << avgExp2(i) << "\t" << avgExp10(i) << "\t" << avgExp100(i) << endl;
		uscita2 << avgCa1(i) << "\t" << avgCa2(i) << "\t" << avgCa10(i) << "\t" << avgCa100(i) << endl;
	}
		
	uscita.close();
	uscita1.close();
	uscita2.close();
   
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
