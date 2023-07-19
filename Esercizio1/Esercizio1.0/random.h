/****************************************************************
*****************************************************************
		_/    _/  _/_/_/  _/       Numerical Simulation Laboratory
	 _/_/  _/ _/       _/       Physics Department
	_/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Igor Vasiljevic
_/    _/  _/_/_/  _/_/_/_/ email: igor.vasiljevic@studenti.unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

// This class contains functions for generating random numbers using the RANNYU algorithm
class Random {

private:
	int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
	// Default constructor
	Random();
	// Destructor
	~Random();
	// Method for inizialize random
	void initializeRandom(const std::string& primesFile, const std::string& seedFile);
	// Method to set the seed for the RNG
	void SetRandom(int * , int, int);
	// Method to save the seed to a file
	void SaveSeed();
	// Method to generate a random number in the range [0,1)
	double Rannyu(void);
	// Method to generate a random number in the range [min,max)
	double Rannyu(double min, double max);
	// Method to generate a random int number in the range [min,max)
	int RannyuInt(int min, int max);
	// Method to generate a random number 1 or -1
	int RannyuNegInt();
	// Method to generate a random number with a Gaussian distribution
	double Gauss(double mean, double sigma);
	// Method to generate a random number with a Exponential distribution
	double Exp(double lambda);
	// Method to generate a random number with a Cauchy distribution
	double Cauchy(double mean, double delta);
	// Method to generate a random angle for pi
	double Pi_accept_reject();
	// Accept-reject method for any pdf
	double pdf(double x);
	double Rannyu_accept_reject(double xmin, double xmax);
};

#endif // __Random__

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
