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
#include <algorithm>

using namespace std;
using namespace arma;

//This function calculate statistical uncertainty estimation
double error(double av, double av2, int n) {
	 if (n==0) {
			return 0;
	 }
	 else {
			return sqrt((av2-av*av)/n);
	 }
}

//Ground state wave function psi100 with n=1, l=0, m=0 for hydrogen atom in natural units
double psi100(double x, double y, double z) {
	return exp(-2.*sqrt(x*x+y*y+z*z));
}
//Ground state wave function psi210 with n=2, l=1, m=0  for hydrogen atom in natural units
double psi210(double x, double y, double z) {
	return (x*x+y*y+z*z)*exp(-sqrt(x*x+y*y+z*z))*cos(atan(y/x))*cos(atan(y/x));
}
//Metropolis loop for Psi100 with random walk moves
vec metropolis100RW(int N, double x0,double y0,double z0, double stepsize,Random rnd) {
	//Initialization of armadillo's vector
	vec mt=zeros<vec>(N);
	vec x=zeros<vec>(N+1);
	vec y=zeros<vec>(N+1);
	vec z=zeros<vec>(N+1);
	//Initialization of variable used
	double xnew=0., xold=0.;
	double znew=0., zold=0.;
	double ynew=0., yold=0.;
	//Initialization of parameter for acceptance of reject
	double A=0.;
	//Initialization of random number
	double r=0.;
	//Initialization of count acceptance
	int count=0;
	//Initialization of count of total simulation
	int counttot=0;
	//Initialization the acceptance values
	double accpr=0.;
	
	//Starting point
	x(0)=x0;
	y(0)=y0;
	z(0)=z0;
	
	//Main Metropolis loop
	for(unsigned int i=0; i<N; i++) {
		//Old values
		xold=x(i);
		yold=y(i);
		zold=z(i);
		//Random walk for new values of position
		xnew=xold+(2.*rnd.Rannyu()-1)*stepsize;
		ynew=yold+(2.*rnd.Rannyu()-1)*stepsize;
		znew=zold+(2.*rnd.Rannyu()-1)*stepsize;
		//Calculation of acceptance approssimatly 50%
		A=min(1.,psi100(xnew,ynew,znew)/psi100(xold,yold,zold));
		//Random values [0,1)
		r=rnd.Rannyu();
		//accept the move if r<=A then move
		if (r<=A) {
			x(i+1)=xnew;
			y(i+1)=ynew;
			z(i+1)=znew;
			//incremence count of accept moves
			count++;
		}
		//If reject, then stay and try again
		else {
			x(i+1)=xold;
			y(i+1)=yold;
			z(i+1)=zold;
		}
		//incremence count of total try
		counttot++;
	}
	//Calculate the values of r
	for(unsigned int i=0; i<N; i++) {
		mt(i)=sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i));
	}
	//Print the acceptance ratio
	accpr=static_cast<double>(count)/(counttot);
	cout << "Acceptance ratio 100 Random Walk: " << accpr << endl;
	return mt;
}

//Metropolis loop for Psi210 with random walk moves
vec metropolis210RW(int N, double x0,double y0,double z0, double stepsize,Random rnd) {
	//Initialization of armadillo's vector
	vec mt=zeros<vec>(N);
	vec x=zeros<vec>(N+1);
	vec y=zeros<vec>(N+1);
	vec z=zeros<vec>(N+1);
	//Initialization of variable used
	double xnew=0., xold=0.;
	double znew=0., zold=0.;
	double ynew=0., yold=0.;
	//Initialization of parameter for acceptance of reject
	double A=0.;
	//Initialization of random number
	double r=0.;
	//Initialization of count acceptance
	int count=0;
	//Initialization of count of total simulation
	int counttot=0;
	//Initialization the acceptance values
	double accpr=0.;
	
	//Starting position
	x(0)=x0;
	y(0)=y0;
	z(0)=z0;
	
	//Main Metropolis loop
	for(unsigned int i=0; i<N; i++) {
		//old position
		xold=x(i);
		yold=y(i);
		zold=z(i);
		//New position generate by random walk
		xnew=xold+(2.*rnd.Rannyu()-1)*stepsize;
		ynew=yold+(2.*rnd.Rannyu()-1)*stepsize;
		znew=zold+(2.*rnd.Rannyu()-1)*stepsize;
		//accept the move if r<=A then move
		A=min(1.,psi210(xnew,ynew,znew)/psi210(xold,yold,zold));
		//random number between [0,1)
		r=rnd.Rannyu();
		//accept the move if r<=A then move
		if (r<=A) {
			x(i+1)=xnew;
			y(i+1)=ynew;
			z(i+1)=znew;
			//increment count acceptance move
			count++;
		}
		//If reject, stay and try again
		else {
			x(i+1)=xold;
			y(i+1)=yold;
			z(i+1)=zold;
		}
		//incremence count of total try
		counttot++;
	}
	//calculate the values of r
	for(unsigned int i=0; i<N; i++) {
		mt(i)=sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i));
	}
	//print the acceptance ratio for 210 random walk
	accpr=static_cast<double>(count)/(counttot);
	cout << "Acceptance ratio 210 Random Walk: " << accpr << endl;
	return mt;
}

//Metropolis loop for Psi100 with Gaussian moves
vec metropolis100Gauss(int N, double x0,double y0,double z0, double sigma,Random rnd) {
	//Initialize armadillo's vector
	vec mt=zeros<vec>(N);
	vec x=zeros<vec>(N+1);
	vec y=zeros<vec>(N+1);
	vec z=zeros<vec>(N+1);
	//Initialize values
	double xnew=0., xold=0.;
	double znew=0., zold=0.;
	double ynew=0., yold=0.;
	//Initialization of parameter for acceptance of reject
	double A=0.;
	//Initialization of random number
	double r=0.;
	//Initialization of count acceptance
	int count=0;
	//Initialization of count of total simulation
	int counttot=0;
	//Initialization the acceptance values
	double accpr=0.;
	//Starting values
	x(0)=x0;
	y(0)=y0;
	z(0)=z0;
	//Main metropolis loop
	for(unsigned int i=0; i<N; i++) {
		//old positions
		xold=x(i);
		yold=y(i);
		zold=z(i);
		//New position generate by random values of Gaussian distribution centered in xold with sigma=sigma selected in order to have approssimatly 50%
		xnew=rnd.Gauss(xold,sigma);
		ynew=rnd.Gauss(yold,sigma);
		znew=rnd.Gauss(zold,sigma);
		//accept the move if r<=A then move
		A=min(1.,psi100(xnew,ynew,znew)/psi100(xold,yold,zold));
		//random number between [0,1)
		r=rnd.Rannyu();
		//accept the move if r<=A then move
		if (r<=A) {
			x(i+1)=xnew;
			y(i+1)=ynew;
			z(i+1)=znew;
			//increment count of accept moves
			count++;
		}
		else {
			//If reject, stay and try again
			x(i+1)=xold;
			y(i+1)=yold;
			z(i+1)=zold;
		}
		//increment total count of moves
		counttot++;
	}
	//calculate value of r
	for(unsigned int i=0; i<N; i++) {
		mt(i)=sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i));
	}
	//print acceptance ratio for 100 with Gauss
	accpr=static_cast<double>(count)/(counttot);
	cout << "Acceptance ratio 100 Gauss: " << accpr << endl;
	return mt;
}

//Metropolis loop for Psi210 with Gaussian moves
vec metropolis210Gauss(int N, double x0,double y0,double z0, double sigma,Random rnd) {
	//initialize armadillo's vector
	vec mt=zeros<vec>(N);
	vec x=zeros<vec>(N+1);
	vec y=zeros<vec>(N+1);
	vec z=zeros<vec>(N+1);
	//Initialize variables
	double xnew=0., xold=0.;
	double znew=0., zold=0.;
	double ynew=0., yold=0.;
	//Initialization of parameter for acceptance of reject
	double A=0.;
	//Initialization of random number
	double r=0.;
	//Initialization of count acceptance
	int count=0;
	//Initialization of count of total simulation
	int counttot=0;
	//Initialization the acceptance values
	double accpr=0.;
	//Starting position
	x(0)=x0;
	y(0)=y0;
	z(0)=z0;
	//Main metropolis loop
	for(unsigned int i=0; i<N; i++) {
		//old positions
		xold=x(i);
		yold=y(i);
		zold=z(i);
		//New position generate by random values of Gaussian distribution centered in xold with sigma=sigma selected in order to have approssimatly 50%
		xnew=rnd.Gauss(xold,sigma);
		ynew=rnd.Gauss(yold,sigma);
		znew=rnd.Gauss(zold,sigma);
		//accept the move if r<=A then move
		A=min(1.,psi210(xnew,ynew,znew)/psi210(xold,yold,zold));
		//random values [0,1)
		r=rnd.Rannyu();
		//accept the move if r<=A then move
		if (r<=A) {
			x(i+1)=xnew;
			y(i+1)=ynew;
			z(i+1)=znew;
			//increment count of accept moves
			count++;
		}
		else {
			//If reject, stay and try again
			x(i+1)=xold;
			y(i+1)=yold;
			z(i+1)=zold;
		}
		//increment total count of moves
		counttot++;
	}
	//calculate values of r
	for(unsigned int i=0; i<N; i++) {
		mt(i)=sqrt(x(i)*x(i)+y(i)*y(i)+z(i)*z(i));
	}
	//print acceptance ratio for 210 with Gauss
	accpr=static_cast<double>(count)/(counttot);
	cout << "Acceptance ratio 210 Gauss: " << accpr << endl;
	return mt;
}

// This function returns the statistical uncertainty and the mean given number of blocks N and the number of throws in each block L and a armadillo's vector v with data blocking method
pair<vec,vec> ArraySumErr(int N, int L, const vec& v){
	//initialize armadillo's vector
	vec avg = zeros<vec>(N);
	vec avg2 = zeros<vec>(N);
	vec sum_prog = zeros<vec>(N);
	vec sum2_prog = zeros<vec>(N);
	vec err = zeros<vec>(N);
	//initialize variables
	double sum1=0.;
	int k=0;
	//Loop summing values
	for(unsigned int i=0; i<N;i++) {
		sum1=0.;
		for(unsigned int j=0;j<L;j++) {
			k = j+i*L;
			sum1 += v(k);
		}
		avg(i)=sum1/L;
		avg2(i)=avg(i)*avg(i);
	}
	//Loop for mean ad err in each block
	for(unsigned int i=0; i<N;i++) {
		for(unsigned int j=0;j<i+1;j++) {
			sum_prog(i) += avg(j);
			sum2_prog(i) += avg2(j);
		}
		sum_prog(i) /= (i+1);
		sum2_prog(i) /= (i+1);
		err(i)=error(sum_prog(i), sum2_prog(i),i);
	}
	return make_pair(err, sum_prog);
}

//Those 4 function are just for saving the coordination inside a matrix. I should add it inside a main loop of the metropolis algorithm, so the code will be a bit rusty. If you want just to do the metropolis algorithm then comment those function.

mat Coord100RW(int N, double x0,double y0,double z0, double stepsize,Random rnd) {
	mat matrix=zeros<mat>(N,3);
	vec x=zeros<vec>(N+1);
	vec y=zeros<vec>(N+1);
	vec z=zeros<vec>(N+1);
	double xnew=0., xold=0.;
	double znew=0., zold=0.;
	double ynew=0., yold=0.;
	double A=0.;
	double r=0.;
	int count=0;
	int counttot=0;
	double accpr=0.;
	
	x(0)=x0;
	y(0)=y0;
	z(0)=z0;
	
	for(unsigned int i=0; i<N; i++) {
		xold=x(i);
		yold=y(i);
		zold=z(i);
		xnew=xold+(2.*rnd.Rannyu()-1)*stepsize;
		ynew=yold+(2.*rnd.Rannyu()-1)*stepsize;
		znew=zold+(2.*rnd.Rannyu()-1)*stepsize;
		A=min(1.,psi100(xnew,ynew,znew)/psi100(xold,yold,zold));
		r=rnd.Rannyu();
		if (r<=A) {
			x(i+1)=xnew;
			y(i+1)=ynew;
			z(i+1)=znew;
		}
		else {
			x(i+1)=xold;
			y(i+1)=yold;
			z(i+1)=zold;
		}
	}
	for (unsigned int i=0; i<N; i++) {
		matrix(i, 0) = x(i);
		matrix(i, 1) = y(i);
		matrix(i, 2) = z(i);
	}
	return matrix;
}

mat Coord210RW(int N, double x0,double y0,double z0, double stepsize,Random rnd) {
	mat matrix=zeros<mat>(N,3);
	vec x=zeros<vec>(N+1);
	vec y=zeros<vec>(N+1);
	vec z=zeros<vec>(N+1);
	double xnew=0., xold=0.;
	double znew=0., zold=0.;
	double ynew=0., yold=0.;
	double A=0.;
	double r=0.;
	int count=0;
	int counttot=0;
	double accpr=0.;
	
	x(0)=x0;
	y(0)=y0;
	z(0)=z0;
	
	for(unsigned int i=0; i<N; i++) {
		xold=x(i);
		yold=y(i);
		zold=z(i);
		xnew=xold+(2.*rnd.Rannyu()-1)*stepsize;
		ynew=yold+(2.*rnd.Rannyu()-1)*stepsize;
		znew=zold+(2.*rnd.Rannyu()-1)*stepsize;
		A=std::min(1.,psi210(xnew,ynew,znew)/psi210(xold,yold,zold));
		r=rnd.Rannyu();
		if (r<=A) {
			x(i+1)=xnew;
			y(i+1)=ynew;
			z(i+1)=znew;
		}
		else {
			x(i+1)=xold;
			y(i+1)=yold;
			z(i+1)=zold;
		}
	}
	for (unsigned int i=0; i<N; i++) {
		matrix(i, 0) = x(i);
		matrix(i, 1) = y(i);
		matrix(i, 2) = z(i);
	}
	return matrix;
}

mat Coord100Gauss(int N, double x0,double y0,double z0, double sigma,Random rnd) {
	mat matrix=zeros<mat>(N,3);
	vec x=zeros<vec>(N+1);
	vec y=zeros<vec>(N+1);
	vec z=zeros<vec>(N+1);
	double xnew=0., xold=0.;
	double znew=0., zold=0.;
	double ynew=0., yold=0.;
	double A=0.;
	double r=0.;
	int count=0;
	int counttot=0;
	double accpr=0.;
	
	x(0)=x0;
	y(0)=y0;
	z(0)=z0;
	
	for(unsigned int i=0; i<N; i++) {
		xold=x(i);
		yold=y(i);
		zold=z(i);
		xnew=rnd.Gauss(xold,sigma);
		ynew=rnd.Gauss(yold,sigma);
		znew=rnd.Gauss(zold,sigma);
		A=std::min(1.,psi100(xnew,ynew,znew)/psi100(xold,yold,zold));
		r=rnd.Rannyu();
		if (r<=A) {
			x(i+1)=xnew;
			y(i+1)=ynew;
			z(i+1)=znew;
		}
		else {
			x(i+1)=xold;
			y(i+1)=yold;
			z(i+1)=zold;
		}
	}
	for (unsigned int i=0; i<N; i++) {
		matrix(i, 0) = x(i);
		matrix(i, 1) = y(i);
		matrix(i, 2) = z(i);
	}
	
	return matrix;
}

mat Coord210Gauss(int N, double x0,double y0,double z0, double sigma,Random rnd) {
	mat matrix=zeros<mat>(N,3);
	vec x=zeros<vec>(N+1);
	vec y=zeros<vec>(N+1);
	vec z=zeros<vec>(N+1);
	double xnew=0., xold=0.;
	double znew=0., zold=0.;
	double ynew=0., yold=0.;
	double A=0.;
	double r=0.;
	int count=0;
	int counttot=0;
	double accpr=0.;
	
	x(0)=x0;
	y(0)=y0;
	z(0)=z0;
	
	for(unsigned int i=0; i<N; i++) {
		xold=x(i);
		yold=y(i);
		zold=z(i);
		xnew=rnd.Gauss(xold,sigma);
		ynew=rnd.Gauss(yold,sigma);
		znew=rnd.Gauss(zold,sigma);
		A=std::min(1.,psi210(xnew,ynew,znew)/psi210(xold,yold,zold));
		r=rnd.Rannyu();
		if (r<=A) {
			x(i+1)=xnew;
			y(i+1)=ynew;
			z(i+1)=znew;
		}
		else {
			x(i+1)=xold;
			y(i+1)=yold;
			z(i+1)=zold;
		}
	}

	for (unsigned int i=0; i<N; i++) {
		matrix(i, 0) = x(i);
		matrix(i, 1) = y(i);
		matrix(i, 2) = z(i);
	}
	
	return matrix;
}

//Function to write results on a file named filename
void WriteToFile(const vec& avg100, const vec& err100,const vec& avg210, const vec& err210,const string& filename) {
		ofstream file(filename);
		if (file.is_open()) {
			for(unsigned int i=0; i<avg100.size(); i++) {
				file << avg100(i) << "\t" << err100(i) << "\t" << avg210(i)<< "\t" << err210(i) <<endl;
				}
				file.close();
		}
}

//Function to write results on a file named filename
void WriteMatrixToFile(const mat& matrix, int M, int N,const string& filename) {
		ofstream file(filename);
		if (file.is_open()) {
			for(unsigned int i=0; i<N; i++) {
				for (unsigned int i = 0; i < M; ++i) {
					for (unsigned j = 0; j < N; ++j) {
						file << matrix(i, j) << "\t";
					}
					file << endl;
				}
			file.close();
			}
		}
}
	
//Function to write results on a file named filename
void WriteInstaToFile(const vec& avgr100, const vec& avgr210,const string& filename) {
	ofstream file(filename);
	if (file.is_open()) {
		for (unsigned int i = 0; i < avgr100.size(); ++i) {
			file << avgr100(i) << "\t" << avgr210(i) << endl;
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
