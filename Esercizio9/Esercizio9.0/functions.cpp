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
#include <vector>
#include <algorithm>

using namespace std;
using namespace arma;

//calculates the fitness of a Path based on the total distance traveled
void Path::calculateFitness(const std::vector<City>& cityList) {
	fitness = 0.0;
	for (size_t i = 0; i < cities.size() - 1; i++) {
		int cityIndex1 = cities[i];
		int cityIndex2 = cities[i + 1];
		double distance = calculateDistance(cityList[cityIndex1], cityList[cityIndex2]);
		fitness += distance;
	}
	//Add distance from the last city back to the first city
	int lastIndex = cities[cities.size() - 1];
	int firstIndex = cities[0];
	double lastDistance = calculateDistance(cityList[lastIndex], cityList[firstIndex]);
	fitness += lastDistance;
}
//calculate the Euclidean distance between two cities
double Path::calculateDistance(const City& city1, const City& city2) {
		double dx = city1.x - city2.x;
		double dy = city1.y - city2.y;
		return sqrt(dx * dx + dy * dy);
}

//Function that generates a random permutation of cities, using swap with random int
vector<int> generateRandomPermutation(int size, Random rnd) {
	vector<int> permutation(size);
	for (int i = 0; i < size; i++) {
		permutation[i] = i;
	}
	for(int i=size-1; i>0; --i) {
		swap(permutation[i],permutation[rnd.RannyuInt(0,i+1)]);
	}
	return permutation;
}

//Function that create the starting population
vector<Path> createInitialPopulation(int populationSize, int numCities, Random rnd) {
	vector<Path> population(populationSize);
	for (int i = 0; i < populationSize; i++) {
		population[i].cities = generateRandomPermutation(numCities,rnd);
	}
	return population;
}

//Function that check if a path is satisfying the bond (visit ones) if the city is visited is marked with 1, else with 0
bool checkFunction(const Path& path) {
	vector<int> visited(path.cities.size(), 0);
	for (int city : path.cities) {
		if (visited[city] == 1) {
			return false;
		}
		visited[city] = 1;
	}
	return true;
}

//Function that evaluate the fitness of the population
void EvFitness(vector<Path>& population, const vector<City>& cityList) {
	for (Path& path : population) {
		path.calculateFitness(cityList);
	}
}

//Function for the tournament selection
//compare the function to determinate the order of the elements. return true if path1.fitness is less than path2.fitness
Path tournamentSelection(const vector<Path>& population, int tournamentSize, Random& rnd) {
	int populationSize = population.size();
	vector<int> tournamentIndices(tournamentSize);
	for (int i = 0; i < tournamentSize; i++) {
		tournamentIndices[i] = int(rnd.Rannyu(0, populationSize));
	}
	auto minElement = min_element(tournamentIndices.begin(), tournamentIndices.end(), [&](int i, int j) {
		return population[i].fitness < population[j].fitness; });
	return population[*minElement];
}

// Function for order crossover
Path Crossover(const Path& parent1, const Path& parent2, Random rnd) {
		int size = parent1.cities.size();
		int startPos = int(rnd.Rannyu(0, size - 1));
		int endPos = int(rnd.Rannyu(startPos + 1, size));

		vector<int> child(size, -1);

		//Copy the subsequence from parent1 to child
		for (int i = startPos; i < endPos; i++) {
				child[i] = parent1.cities[i];
		}

		//Fill the remaining positions with the cities from parent2
		int j = 0;
		for (int i = 0; i < size; i++) {
				if (j == startPos) {
						j = endPos;
				}
				if (find(child.begin(), child.end(), parent2.cities[i]) == child.end()) {
						child[j] = parent2.cities[i];
						j++;
				}
		}

		Path offspring;
		offspring.cities = child;
		return offspring;
}

//Function for swap mutation usign swap standard function
void swapMutation(Path& path, Random& rnd) {
	int size = path.cities.size();
	int pos1 = int(rnd.Rannyu(1, size - 1));
	int pos2 = int(rnd.Rannyu(1, size - 1));
	swap(path.cities[pos1], path.cities[pos2]);
}

//Function for pair permutation mutation usign swap standard function
void pairPermutationMutation(Path& path, Random& rnd) {
	int size = path.cities.size();
	int pos = int(rnd.Rannyu(1, size - 1));
	swap(path.cities[pos], path.cities[pos + 1]);
}

//Function for shift mutation using rotate stantard function
void shiftMutation(Path& path, int m, int n, Random& rnd) {
	int size = path.cities.size();
	int startPos = int(rnd.Rannyu(1, size - m));
	rotate(path.cities.begin() + startPos, path.cities.begin() + startPos + n, path.cities.begin() + startPos + m + 1);
}

//Function for permutation mutation using reverse standard function
void permutationMutation(Path& path, int m, Random& rnd) {
	int size = path.cities.size();
	int startPos1 = int(rnd.Rannyu(1, size - 2 * m));
	int startPos2 = startPos1 + m;
	reverse(path.cities.begin() + startPos1, path.cities.begin() + startPos1 + m);
	reverse(path.cities.begin() + startPos2, path.cities.begin() + startPos2 + m);
}

//Function for inversion mutation using reverse standard function
void inversionMutation(Path& path, int m, Random& rnd) {
	int size = path.cities.size();
	int startPos = int(rnd.Rannyu(1, size - m));
	reverse(path.cities.begin() + startPos, path.cities.begin() + startPos + m);
}

//Function for getting the best path from the population (lowest fitness value)
//The function comparison function to determinate the order of the elements and return true if path1.fitness is less than path2.fitness
Path getBestPath(const vector<Path>& population) {
	auto minElement = std::min_element(population.begin(), population.end(), [&](const Path& path1, const Path& path2) { return path1.fitness < path2.fitness; });
	return *minElement;
}

//Function to write path length and average path length to a file
void writeDataToFile(const vector<double>& avgPathLengths, const string& filename) {
	ofstream file(filename);
	if (file.is_open()) {
		for (size_t i = 0; i < avgPathLengths.size(); i++) {
			file << i << " "  << avgPathLengths[i] << endl;
		}
		file.close();
	}
}

//Function for writing the best path in a file take filename as a string
void writeBestPathToFile(const Path& path, const vector<City>& cityList, const string& filename) {
	ofstream file(filename);
	if (file.is_open()) {
		for (int city : path.cities) {
			file << cityList[city].x << " " << cityList[city].y << std::endl;
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
