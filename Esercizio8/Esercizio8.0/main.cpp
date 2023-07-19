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
	
	//Initialize variables
	double x_init=0.; //Initial position
	double sigma=1.; //Starting values of sigma
	double mu=0.5;  //Starting values of mu
	double best_energy=0.; //Starting energy
	const int M=10000; //Number of Metropolis moves
	const int N=100; //Number of blocks
	const int equi_step=1000; //Number of equilibration step
	double stepsize=2.0; //Step random walk for metropolis
	int L=0.;  //Number of moves in each block
	L=int(M/N);
	int N_SA_steps=10000; //Number of SA step
	double starting_temp=25.; //Starting temperature
	double cooling_factor=0.9999; //Cooling factor
	
	//Initialize armadillo's vector
	vec sigma_traj = zeros<vec>(N_SA_steps);
	vec mu_traj = zeros<vec>(N_SA_steps);
	vec energy_traj=zeros<vec>(N_SA_steps);
	
	vec values=zeros<vec>(3);
	
	//Performing optimizeSA for sigma and mu
	optimizeSA(M,x_init,stepsize,equi_step,sigma,mu,N_SA_steps,starting_temp,cooling_factor, rnd,energy_traj,sigma_traj, mu_traj, values);
	
	//Saving best sigma, mu and energy
	sigma=values(0);
	mu=values(1);
	best_energy=values(2);
	
	//Initialize armadillo's vector
	vec energy_values = zeros<vec>(M);
	vec energy_err = zeros<vec>(N);
	vec energy_mean = zeros<vec>(N);
	
	vec x_pos=zeros<vec>(M);
	//Metropolis with best sigma and best mu. Adding the optional parameter for saving positioning values and output the accept ratio
	energy_values=metropolis(M,stepsize,x_init,sigma, mu,equi_step,rnd,x_pos,1.);
	
	//Data blocking statistic
	pair<vec,vec> resultSumC=ArraySumErr(N,L,energy_values);
	energy_err = resultSumC.first;
	energy_mean = resultSumC.second;
	
	//Write in a file the energy metropolis mean and error
	WriteEnergyToFile(energy_mean, energy_err,"datiEnergy.dat");
	
	//Writing the trajectory of energy, sigma and mu
	WriteTrajToFile(energy_traj, "datiEnergyTraj.dat");
	WriteTrajToFile(sigma_traj, "datiSigmaTraj.dat");
	WriteTrajToFile(mu_traj, "datiMuTraj.dat");
	
	double xmin=min(x_pos);
	double xmax=max(x_pos);
	
	//Histagram filled
	int nbins = 100;
	vec histogram = calculateHistogram(x_pos, nbins, xmin, xmax,sigma,mu);
	string histoFilename = "histo.dat";
	writeHistogramToFile(histogram, xmin,xmax, histoFilename);
	
	//Analyci solution of psi square
	string AnalyticFilename = "analyic.dat";
	calculateAnalytic(M,xmin,xmax,sigma,mu,AnalyticFilename);
	
	//Solving the matrix problem
	string MatFilename = "matValues.dat";
	calculateEquation(M, x_pos,MatFilename);

	//Print optimized parameter values and energy
	cout << "Best Sigma: " << sigma << endl;
	cout << "Best Mu: " << mu << endl;
	cout << "Best Energy: " << best_energy << endl;
	
	//Save random number generator seed
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
