/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Igor Vasiljevic
_/    _/  _/_/_/  _/_/_/_/ email: igor.vasiljevic@studenti.unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Functions__
#define __Functions__

#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

//Calculates the wavefunction at position x given: position, standard deviation and mean
double psi_prova(double x, double sigma, double mu);
//Evaluetes the second derivative of the wavefunction with respect to x
double psi_derivative(double x, double sigma, double mu);
//Calculates the statistical error
double error(double, double, int);
//Performs equilibration steps using the Metropolis algorithm and gives the position after equistep
double equilibration(int equistep, double stepsize, double x, double sigma, double mu, Random rnd);
//Calculates the energy of a particle Kin/psi + potential given position
double energy(double x, double sigma, double mu);
//Calculates the potential energy for each position in a vector
void potential(const vec& x, vec& v, int N);
//Performs N Metropolis Monte Carlo simulations using random walk moves with stepsize, starting with initial position x
vec metropolis(int N, double stepsize, double x, double sigma, double mu,int num_equilibration_steps,Random rnd, vec& x_pos, double opPar);
//Optimizes the parameters using simulated annealing, this function perform N_Sa_step number of SA, every SA perform N metropolis move
void optimizeSA(int N,double x_init,double stepsize,int equi_step,double sigma, double mu, int N_Sa_step, double t_start, double cooling_factor, Random rnd, vec& energyTraj, vec& sigmaTraj, vec& muTraj, vec& values);
//Calculates a histogram of particle positions and sampling the psi square
vec calculateHistogram(const vec& x_pos, int nbins, double xmin, double xmax, double sigma, double mu);
//Calculates the analytical solution of the wavefunction psi_prova, correctly normalized, between xmin and xmax.
void calculateAnalytic(int M, double xmin, double xmax, double sigma, double mu,const string& filename);
//Calculates the solution of the quantum equation using eig_sym of armadillo's library
void calculateEquation(int M, const vec& x_pos,const string& filename);
//Writes a histogram to a file with centered bins
void writeHistogramToFile(const vec& histogram, double xmin, double xmax, const string& filename);
//Performs statistical analysis on a vector of values in each blocks
pair<vec,vec> ArraySumErr(int,int, const vec& v);
//Function to write results on a file named filename
void WriteTrajToFile(const vec& traj, const string& filename);
void WriteEnergyToFile(const vec& mean, const vec& err, const string& filename);


#endif // __Functions__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Igor Vasiljevic
_/    _/  _/_/_/  _/_/_/_/ email: igor.vasiljevic@studenti.unimi.it
*****************************************************************
*****************************************************************/
