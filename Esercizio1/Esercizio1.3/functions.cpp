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

//Generate M random number between [0,1)
void GenerateRandomNumber(int M, vec& values, Random rnd) {
	for(unsigned int i=0; i<M;i++) {
		values(i)=rnd.Rannyu();
		return;
	}
}
//Generate M random angle without pi
void GenerateRandomNumberTheta(int M, vec& values, Random rnd) {
	for(unsigned int i=0; i<M;i++) {
		values(i)=rnd.Pi_accept_reject();
	}
	return;
}
//Generate M random number between [0,d)
void GenerateRandomNumberTimes(int M, vec& values,double d, Random rnd) {
	for(unsigned int i=0; i<M;i++) {
		values(i)=d*rnd.Rannyu();
	}
	return;
}

// This function returns the statistical uncertainty and the mean of pi given number of blocks N,the number of throws in each block L, the length of a stick L1, the leight of gap lines and two armadillo's vector v1 (random position) and v2 (random angle) with data blocking method
void DataBlocking(int N, int L,double L1, double d, const vec& v1,const vec& v2, vec& mean, vec& err) {
	vec avg = zeros<vec>(N);
	vec avg2 = zeros<vec>(N);
	vec sum2_prog = zeros<vec>(N);
	double sum=0.;
	int Nhit=0;
	int k=0;
	for(unsigned int i=0; i<N;i++) {
		Nhit=0.;
		for(unsigned int j=0;j<L;j++) {
				k = j+i*L;
				if((v1(k) + L1*0.5*sin(v2(k))) >= d or (v1(k)-L1*0.5*sin(v2(k)))<=0.) {
					Nhit+=1.;
				}
		}
		avg(i)=2*L1*L/(Nhit*d);
		avg2(i)=avg(i)*avg(i);
	}
	for(unsigned int i=0; i<N;i++) {
		 for(unsigned int j=0;j<i+1;j++) {
				mean(i) += avg(j);
				sum2_prog(i) += avg2(j);
		 }
		 mean(i) /= (i+1);
		 sum2_prog(i) /= (i+1);
		 err(i)=error(mean(i), sum2_prog(i),i);
	}
	return;
}

//Function to write results on a file named filename
void WriteToFile(const vec& mean, const vec& err, const string& filename) {
		ofstream file(filename);
		if (file.is_open()) {
			file << "sum" << "\t" << "err" << endl;
			for(unsigned int i=0; i<mean.size(); i++) {
						file << mean(i) << "\t" << err(i) << endl;
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
