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

#include <vector>
#include "random.h"

using namespace std;

//the City class is used to contain a coordinates of a city(x, y)
class City {
public:
		double x, y;
};

//represents a path consisting of a vector of city indices and its fitness value.
class Path {
public:
		vector<int> cities;
		double fitness;
	//function calculates the fitness of the path based on the total distance traveled between cities.
		void calculateFitness(const vector<City>& cityList);
private:
	//calculate the distance between two cities
		double calculateDistance(const City& city1, const City& city2);
};

//generates a random permutation of cities (shuffle does problems with armadillo's vector and the rnd generator, so i made the function)
vector<int> generateRandomPermutation(int size, Random rnd);
//creates the initial population of paths
vector<Path> createInitialPopulation(int populationSize, int numCities, Random rnd);
//checks whether a path satisfies the constraint of visiting each city only once
bool checkFunction(const Path& path);
//calculates the fitness of each path in the population based on the city positions
void EvFitness(vector<Path>& population, const vector<City>& cityList);
//performs tournament selection to choose the fittest paths for mating
Path tournamentSelection(const vector<Path>& population, int tournamentSize, Random& rnd);
//performs order crossover to create new paths from two parent paths during reproduction
Path Crossover(const Path& parent1, const Path& parent2, Random rnd);
//swap mutation is a mutation operator that randomly selects two positions in the path and swaps the cities at those positions
void swapMutation(Path& path, Random& rnd);
//pair permutation mutation selects a random position in the path and swaps the city at that position with the city at the next position
void pairPermutationMutation(Path& path, Random& rnd);
//shift mutation randomly selects a starting position in the path and shifts a m cities by n positions and it rotates a segment of the path
void shiftMutation(Path& path, int m, int n, Random& rnd);
//permutation mutation select a random starting position in the path and reverses the order of the m cities
void permutationMutation(Path& path, int m, Random& rnd);
//select a random starting position in the path and reverses the order of the next m cities
void inversionMutation(Path& path, int m, Random& rnd);
//return the best path (lowest fitness value) from the population
Path getBestPath(const vector<Path>& population);
//writes data (average path lengths) in a file
void writeDataToFile(const vector<double>& avgPathLengths, const string& filename);
//writes the best path in a file
void writeBestPathToFile(const Path& path, const vector<City>& cityList, const string& filename);

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
