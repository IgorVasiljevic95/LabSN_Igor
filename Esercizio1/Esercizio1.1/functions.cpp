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
	}
}

// This function returns the statistical uncertainty and the mean given number of blocks N and the number of throws in each block L and a armadillo's vector v with data blocking method
void DataBlocking(int N, int L, const vec& v, vec& mean, vec& err) {
	vec avg = zeros<vec>(N);
	vec avg2 = zeros<vec>(N);
	vec sum2_prog = zeros<vec>(N);
	double sum=0.;
	int k=0;
	for(unsigned int i=0; i<N;i++) {
		 sum=0.;
		 for(unsigned int j=0;j<L;j++) {
				k = j+i*L;
				sum += (v(k)-0.5)*(v(k)-0.5);
		 }
		 avg(i)=sum/L;
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
