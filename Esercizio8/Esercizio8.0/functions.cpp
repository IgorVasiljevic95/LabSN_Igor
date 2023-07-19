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

//Calculates the wavefunction at position x given: position, standard deviation and mean
double psi_prova(double x, double sigma, double mu) {
	double exp1=0.;
	double exp2=0.;
	double psi=0.;
	exp1=exp(-pow(x-mu,2)/(2.*sigma*sigma));
	exp2=exp(-pow(x+mu,2)/(2.*sigma*sigma));
	psi=exp1+exp2;
	return psi;
}
//Evaluetes the second derivative explicitly calculated of the wavefunction with respect to x
double psi_derivative(double x, double sigma, double mu) {
	double exp1=0.;
	double exp2=0.;
	double d2Psi=0.;
	exp1=exp(-pow(x-mu,2)/(2.0*sigma*sigma));
	exp2=exp(-pow(x+mu,2)/(2.0*sigma*sigma));
	d2Psi = (1/pow(sigma,4))*((pow(x-mu,2)-sigma*sigma)*exp1+(pow(x+mu,2)-sigma*sigma)*exp2);
	return d2Psi;
}
//Calculates the energy of a particle Kin/psi + potential given position
double energy(double x, double sigma, double mu) {
	double pot=pow(x,4)-2.5*x*x;
	double kin=-0.5*psi_derivative(x, sigma, mu);
	return kin/psi_prova(x,sigma,mu)+pot;
}
//Calculates the potential energy for each position in a vector
void potential(const vec& x, vec& v, int N) {
	for(unsigned int i=0; i<N; i++) {
		v(i)=pow(x(i),4)-2.5*x(i)*x(i);
	}
	return;
}

//Calculates the statistical error
double error(double av, double av2, int n) {
   if (n==0) {
      return 0;
   }
   else {
      return sqrt((av2-av*av)/(n));
   }
}

//Performs equilibration steps using the Metropolis algorithm and gives the position after equistep
double equilibration(int equistep, double stepsize, double xinit, double sigma, double mu, Random rnd){
	vec x=zeros<vec>(equistep+1);
	double xold=0.;
	double xnew=0.;
	double A=0.;
	double r=0.;
	x(0)=xinit;
	for(int i=0; i < equistep; i++) {
		xold=x(i);
		xnew=xold+(2.*rnd.Rannyu()-1)*stepsize;;
		A=min(1.,(psi_prova(xnew,sigma,mu)*psi_prova(xnew,sigma,mu))/(psi_prova(xold,sigma,mu)*psi_prova(xold,sigma,mu)));
		r=rnd.Rannyu();
		if(r<A){
			x(i+1)=xnew;
		}
		else {
			x(i+1)=xold;
		}
	}
	return x(equistep);
}

//Performs N Metropolis Monte Carlo simulations using random walk moves with stepsize, starting with initial position x. The function tested is psi_prova with need mean mu and deviation standard sigma. Before starting the the function will do the equilibration with num_equilibration_steps steps. Takes x_pos armadillo's vector which is usend only if opPar is used and not 0.0. This is in because in optimizeSA function the function will not exstract the positions values of sampling and will not print the acceptance ratio.
vec metropolis(int N, double stepsize, double xinit, double sigma, double mu,int num_equilibration_steps,Random rnd, vec& x_pos, double opPar=0.0){
	vec x=zeros<vec>(N+1);
	vec energy_final=zeros<vec>(N);
	double xold=0.;
	double xnew=0.;
	int count=0;
	int counttot=0;
	double A=0.;
	double r=0.;
	double accpr=0.;
	
	x(0)=equilibration(num_equilibration_steps, stepsize, xinit, sigma, mu, rnd);
	
	for(unsigned int i=0; i<N; i++) {
		xold=x(i);
		xnew=xold+(2.*rnd.Rannyu()-1)*stepsize;
		A=min(1.,pow(psi_prova(xnew, sigma, mu),2) / pow(psi_prova(xold, sigma, mu),2));
		r=rnd.Rannyu();
		if(r<=A) {
			x(i+1)=xnew;
			count++;
		}
		else {
			x(i+1)=xold;
		}
		counttot++;
	}
	accpr=static_cast<double>(count)/(counttot);
	if (opPar != 0.0) {
		cout << "Acceptance ratio MC: " << accpr << endl;
		for(unsigned int i=0; i<N; i++) {
			x_pos(i)=x(i);
		}
	}
	
	for(unsigned int i=0; i<N; i++) {
		energy_final(i)=energy(x(i),sigma,mu);
	}
	return energy_final;
}

//Optimizes the parameters using simulated annealing, this function perform N_Sa_step number of SA, every SA perform N metropolis move. The trajectory of energy, sigma and mu are saved in energyTraj, sigmaTraj and muTraj respectly. Finally the best values of sigma, mu and energy will be saved in vec values. SA will change temperature with power law
void optimizeSA(int N,double x_init,double stepsize,int equi_step,double sigma, double mu, int N_Sa_step, double t_start, double cooling_factor, Random rnd, vec& energyTraj, vec& sigmaTraj, vec& muTraj,vec& values) {
	
	vec energy=zeros<vec>(N);
	vec sigmaTrajTemp=zeros<vec>(N_Sa_step+1);
	vec muTrajTemp=zeros<vec>(N_Sa_step+1);
	vec energyTrajTemp=zeros<vec>(N_Sa_step+1);

	double energySum=0.;
	double currentSigma=sigma;
	double currentMu=mu;
	double bestSigma = sigma;
	double bestMu = mu;
	double newSigma=0., oldSigma=0.;
	double newMu=0., oldMu=0.;
	double newEnergy=0., oldEnergy;
	double currentEnergy, bestEnergy;
	double temp;
	double accpr=0.;
	int count=0;
	int countot=0;
	vec temp1=zeros<vec>(N);
	
	energy=metropolis(N,stepsize,x_init,sigma,mu,equi_step,rnd,temp1);
	for(unsigned int i=0; i<N; i++) {
		energySum+=energy[i];
	}
	currentEnergy=energySum/N;
	bestEnergy=currentEnergy;
	
	sigmaTrajTemp(0) = currentSigma;
	muTrajTemp(0) = currentMu;
	energyTrajTemp(0) = currentEnergy;

	for(unsigned int i=0; i< N_Sa_step; i++) {
		temp=t_start*pow(cooling_factor,i);
		currentSigma=sigmaTrajTemp(i);
		currentMu=muTrajTemp(i);
		currentEnergy=energyTrajTemp(i);
		newSigma=currentSigma+(2.*rnd.Rannyu()-1)*stepsize;     //rnd.Gauss(currentSigma,stepsize);
		newMu=currentMu+(2.*rnd.Rannyu()-1)*stepsize;         //rnd.Gauss(currentMu,stepsize);
		energySum=0.;
		energy=metropolis(N,stepsize,x_init,newSigma,newMu,equi_step,rnd,temp1);
		for(unsigned int j=0; j<N; j++) {
			energySum+=energy[j];
		}
		newEnergy=energySum/N;
			//Metropolis acceptance test
		double r=rnd.Rannyu();
		double A=min(1.,exp((currentEnergy - newEnergy) / temp));
		if( r<=A) {
			sigmaTrajTemp(i+1)=newSigma;
			muTrajTemp(i+1)=newMu;
			energyTrajTemp(i+1)=newEnergy;
			count++;
		}
		else {
			sigmaTrajTemp(i+1)=currentSigma;
			muTrajTemp(i+1)=currentMu;
			energyTrajTemp(i+1)=currentEnergy;
		}
		// Update best parameters
		if (newEnergy < bestEnergy) {
			bestSigma = newSigma;
			bestMu = newMu;
			bestEnergy = newEnergy;
		}
		countot++;
	// Record parameters
	}
	for(unsigned int i=0; i<N_Sa_step; i++) {
		sigmaTraj(i)=sigmaTrajTemp(i);
		muTraj(i)=muTrajTemp(i);
		energyTraj(i)=energyTrajTemp(i);
	}
		
	accpr=static_cast<double>(count)/(countot);
	cout << "Acceptance ratio SA: " << accpr << endl;

	// Update final parameters
	values(0) = bestSigma;
	values(1) = bestMu;
	values(2)= bestEnergy;
	
	return;
}

//Calculates a histogram wiht nbins of particle positions and sampling the psi square
vec calculateHistogram(const vec& x_pos, int nbins, double xmin, double xmax, double sigma, double mu) {
		vec histogram = zeros<vec>(nbins);
		double dx = (xmax - xmin) / (nbins - 1);
		vec count = zeros<vec>(nbins);

		for (unsigned int i = 0; i < x_pos.n_elem; i++) {
				double psi_squared = pow(psi_prova(x_pos(i), sigma, mu) / (2 * pow(M_PI, 0.25) * sigma), 2);
				int bin_index = static_cast<int>((x_pos(i) - xmin) / dx);
				histogram(bin_index) += psi_squared;
				count(bin_index)++;
		}

		for (unsigned int i = 0; i < nbins; i++) {
				if (count(i) != 0) {
						histogram(i) /= count(i);
				}
		}

		return histogram;
}

//Calculates the analytical solution of the wavefunction psi_prova, correctly normalized, between xmin and xmax.
void calculateAnalytic(int M, double xmin, double xmax, double sigma, double mu,const string& filename) {
	ofstream Analytic(filename);
	vec x = linspace<vec>(xmin, xmax, M);
	
	if (!Analytic.is_open()) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	
	for(unsigned int i=0; i<M; i++) {
		double analytic_solut=pow(psi_prova(x(i),sigma,mu)/(2*pow(M_PI,0.25)*sigma),2);
		Analytic << x(i) << "\t"  << analytic_solut << endl;
	}
	Analytic.close();
}

void calculateEquation(int M, const vec& x_pos,const string& filename) {
	ofstream matValues(filename);
	if (!matValues.is_open()) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	vec pot=zeros<vec>(M);
	double xmin=min(x_pos);
	double xmax=max(x_pos);
	vec x = linspace<vec>(xmin, xmax, M);
	double dx2=x(1)-x(0);
	potential(x,pot,M);
	mat H(M,M,fill::zeros);
	mat V_diag(M, M, fill::zeros);
	// The central differences method: f" = (f_1 - 2*f_0 + f_-1) / dx^2
	mat CDiff = diagmat(ones<vec>(M - 1), -1) - 2 * diagmat(ones<vec>(M), 0) + diagmat(ones<vec>(M - 1), 1);
	H = (-(pow(1, 2)) * CDiff) / (2 *pow(dx2, 2)) + diagmat(pot);
	
	vec eigenvalues=zeros<vec>(M);
	mat eigenvectors(M, M,fill::zeros);
	
	eig_sym(eigenvalues, eigenvectors, H);
	
	eigenvectors=trans(eigenvectors);
	eigenvectors=eigenvectors/sqrt(dx2);
	cout << "Ground state energy: " << eigenvalues(0) << endl;
	cout << "1st excited state energy: " << eigenvalues(1) << std::endl;
	cout << "2nd excited state energy: " << eigenvalues(2) << std::endl;
	
	for (unsigned int i = 0; i < M; i++) {
		double waveFunctionSquared= pow(eigenvectors(0,i), 2);
		matValues << x(i) << "\t" << waveFunctionSquared << endl;
	}
	
	matValues.close();
}

//Writes a histogram to a file with centered bins
void writeHistogramToFile(const vec& histogram, double xmin, double xmax, const string& filename) {
	ofstream Histo(filename);
	if (!Histo.is_open()) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	vec bin_centers=zeros<vec>(histogram.n_elem);
	double dx = (xmax - xmin) / (histogram.n_elem- 1);
	for(unsigned int i = 0; i < histogram.n_elem; i++) {
		double bin_start = xmin + i * dx;
		double bin_end = bin_start + dx;
		//Calculate the middle point
		bin_centers(i) = (bin_start + bin_end) / 2.0;
	}

	for (unsigned int i = 0; i < histogram.n_elem; i++) {
		Histo << bin_centers(i) << "\t" << histogram(i) << endl;
	}
	Histo.close();
}

//Performs statistical analysis on a vector of values in each blocks
pair<vec,vec> ArraySumErr(int N, int L, const vec& v){
   vec avg = zeros<vec>(N);
   vec avg2 = zeros<vec>(N);
   vec sum_prog = zeros<vec>(N);
   vec sum2_prog = zeros<vec>(N);
   vec err = zeros<vec>(N);
   double sum1=0.;
   int k=0;
   for(unsigned int i=0; i<N;i++) {
      sum1=0.;
      for(unsigned int j=0;j<L;j++) {
         k = j+i*L;
         sum1 += v(k);
      }
      avg(i)=sum1/L;
      avg2(i)=avg(i)*avg(i);
   }
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

//Function to write results on a file named filename
void WriteTrajToFile(const vec& traj, const string& filename) {
		ofstream file(filename);
		if (file.is_open()) {
			for(unsigned int i=0; i<traj.size(); i++) {
						file << traj(i) << endl;
				}
				file.close();
		}
}

//Function to write results on a file named filename
void WriteEnergyToFile(const vec& mean, const vec& err, const string& filename) {
		ofstream file(filename);
		if (file.is_open()) {
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
