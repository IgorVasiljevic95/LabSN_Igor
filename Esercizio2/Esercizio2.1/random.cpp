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
#include "random.h"

using namespace std;

// Default constructor, does not perform any action
Random :: Random(){}

// Default destructor, does not perform any action
Random :: ~Random(){}

//This function take seed and primes number from two files for initialize random generator
void Random::initializeRandom(const std::string& primesFile, const std::string& seedFile) {
	int seed[4];
	int p1,p2;
	ifstream primes(primesFile);
	if(primes.is_open()) {
		primes >> p1 >> p2;
	}
	else {
		cerr << "PROBLEM: Unable to open Primes" << endl;
	}
	primes.close();
	
	ifstream input(seedFile);
	string property;
	if (input.is_open()){
		 while ( !input.eof() ){
				input >> property;
				if( property == "RANDOMSEED" ){
					 input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
					 SetRandom(seed,p1,p2);
				}
		 }
		 input.close();
	 } else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

// This function sets the seed and parameters of the random number generator
void Random :: SetRandom(int * s, int p1, int p2){
	m1 = 502;
	m2 = 1521;
	m3 = 4071;
	m4 = 2107;
	l1 = s[0];
	l2 = s[1];
	l3 = s[2];
	l4 = s[3];
	n1 = 0;
	n2 = 0;
	n3 = p1;
	n4 = p2;

	return;
}

void Random :: SaveSeed(){
   // This function saves the current state of the random number generator to a file "seed.out"
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << "RANDOMSEED	" << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Rannyu(void){
	// This function generates a random number in the range [0,1)
	const double twom12=0.000244140625;
	int i1,i2,i3,i4;
	double r;

	i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
	i2 = l2*m4 + l3*m3 + l4*m2 + n2;
	i3 = l3*m4 + l4*m3 + n3;
	i4 = l4*m4 + n4;
	l4 = i4%4096;
	i3 = i3 + i4/4096;
	l3 = i3%4096;
	i2 = i2 + i3/4096;
	l2 = i2%4096;
	l1 = (i1 + i2/4096)%4096;
	r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

	return r;
}

double Random :: Rannyu(double min, double max){
	 // This function generates a random number in the range [min, max)
	 return min+(max-min)*Rannyu();
}

// This function generates a inter random number in the range [min, max)
int Random::RannyuInt(int min,int max) {
	 double s=Rannyu();
	 double x=min+((max+1)-min)*s;
	 return int(x);
}

//This function give randomly one number: 1 or -1
int Random::RannyuNegInt() {
	 double s=Rannyu();
	 if (s > 0.5) {
			return 1;
	 }
	 else {
			return -1;
	 }
}

// This function generates a random number from a Gaussian distribution with given mean and sigma
double Random :: Gauss(double mean, double sigma) {
	 double s=Rannyu();
	 double t=Rannyu();
	 double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
	 return mean + x * sigma;
}

// This function generates a random number from a Exponential distribution with given lambda
double Random :: Exp(double lambda) {
	 double s=Rannyu();
	 return -log(1.-s)/lambda;
}

// This function generates a random number from a Cauchy distribution with given mean and delta
double Random :: Cauchy(double mean, double delta) {
	 double s=Rannyu();
	 return delta*tan(M_PI*(s-0.5))+mean;
}

//This function give random angle with accept reject method
double Random::Pi_accept_reject() {
	double theta=0.;
	double x=0.;
	double y=0.;
	double r=0.;
	do {
		x=Rannyu();
		y=Rannyu();
		r=sqrt(x*x+y*y);
	} while(r>=1.0);
	theta = acos(x/r);
	return theta;
}

double Random::pdf(double x) {
   return sqrt(1-x*x);
}

double Random::Rannyu_accept_reject(double xmin, double xmax) {
   double ymax=pdf(xmin);
   double x=0.;
   double y=0.;
   do {
      x=Rannyu(xmin, xmax);
      y=ymax*Rannyu();
   } while(y>pdf(x));
   return x;
}

/****************************************************************
Cit. Prof. D.E. Galli, Universita' degli Studi di Milano, email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

/****************************************************************
*****************************************************************
		_/    _/  _/_/_/  _/       Numerical Simulation Laboratory
	 _/_/  _/ _/       _/       Physics Department
	_/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Igor Vasiljevic
_/    _/  _/_/_/  _/_/_/_/ email: igor.vasiljevic@studenti.unimi.it
*****************************************************************
*****************************************************************/





/*
double Random::Pi_accept_reject() {
	 double theta=0;
	 double x=Rannyu(-1,1);
	 double y=Rannyu(-1,1);
	 if(x*x+y*y<1) {
			if(y>=0.) {
				 theta=acos(x/(sqrt(x*x+y*y)));
			}
			else {
				 theta=2*M_PI-acos(x/(sqrt(x*x+y*y)));
			}
	 }
	 return theta;
}
*/
