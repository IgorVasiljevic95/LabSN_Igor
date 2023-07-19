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
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include "random.h"
#include "functions.h"

using namespace std;

//In this GA code to solving the Traveling Salesman Problem i calculated the L^2
int main() {

	//Define and Initialize the random seed and primes rnd.initializeRandom
	Random rnd;
	rnd.initializeRandom("Primes","seed.in");
	
	//Define parameters
	int numCities = 34;
	int populationSize = 500;
	int tournamentSize = 5;
	int numGenerations = 150;
	double SwapMutationRate = 0.09;
	double PairPermutationRate = 0.01;
	double ShiftMutationRate = 0.06;
	double PermutationRate = 0.02;
	double InversionRate = 0.08;
	int square_leight = 16.;

	//Generate random cities on a circumference
	vector<City> cityListCircumference(numCities);
	for (unsigned int i = 0; i < numCities; i++) {
		double angle = 2. * M_PI * i / numCities;
		cityListCircumference[i].x = cos(angle);
		cityListCircumference[i].y = sin(angle);
	}

	//Generate random cities inside a square
	vector<City> cityListSquare(numCities);
	for (unsigned int i = 0; i < numCities; i++) {
		cityListSquare[i].x = rnd.Rannyu(-square_leight, square_leight);
		cityListSquare[i].y = rnd.Rannyu(-square_leight, square_leight);
	}
	
	//Create the initial population for both circumference and square
	vector<Path> populationCircumference = createInitialPopulation(populationSize, numCities, rnd);
	vector<Path> populationSquare = createInitialPopulation(populationSize, numCities, rnd);

	//Evaluate the fitness of the initial population for both circumference and square
	EvFitness(populationCircumference, cityListCircumference);
	EvFitness(populationSquare, cityListSquare);

	//Store the best path and average path length for each generation
	vector<double> pathLengthsCircumference(numGenerations);
	vector<double> avgPathLengthsCircumference(numGenerations);
	vector<double> pathLengthsSquare(numGenerations);
	vector<double> avgPathLengthsSquare(numGenerations);

	//Main GA loop for both circumference and square
	//The main GA loop runs for numGenerations, where each generation involves the selection, crossover, and mutation of individuals in the population
	for (int generation = 0; generation < numGenerations; generation++) {
		//The population is sorted based on fitness, and the average path lengths of the best half of the population are calculated
		sort(populationCircumference.begin(), populationCircumference.end(), [&](const Path& path1, const Path& path2) { return path1.fitness < path2.fitness; });

		sort(populationSquare.begin(), populationSquare.end(), [&](const Path& path1, const Path& path2) { return path1.fitness < path2.fitness; });

		//Calculate average path length of the best half of the population for both circumference and square
		double avgPathLengthCircumference = 0.;
		double avgPathLengthSquare = 0.;
		for (int i = 0; i < populationSize / 2; i++) {
			avgPathLengthCircumference += populationCircumference[i].fitness;
			avgPathLengthSquare += populationSquare[i].fitness;
		}
		avgPathLengthCircumference /= (populationSize / 2);
		avgPathLengthSquare /= (populationSize / 2);

		pathLengthsCircumference[generation] = populationCircumference[0].fitness;
		avgPathLengthsCircumference[generation] = avgPathLengthCircumference;

		pathLengthsSquare[generation] = populationSquare[0].fitness;
		avgPathLengthsSquare[generation] = avgPathLengthSquare;

		//Tournament selection, crossover, and mutation for both scenarios
		//Tournament selection is performed to choose parents for crossover, and order crossover is used to create new offspring
		vector<Path> offspringCircumference;
		vector<Path> offspringSquare;
		while (offspringCircumference.size() < populationSize || offspringSquare.size() < populationSize) {
			//Selection for both circumference and square
			Path parent1Circumference = tournamentSelection(populationCircumference, tournamentSize, rnd);
			Path parent2Circumference = tournamentSelection(populationCircumference, tournamentSize, rnd);

			Path parent1Square = tournamentSelection(populationSquare, tournamentSize, rnd);
			Path parent2Square = tournamentSelection(populationSquare, tournamentSize, rnd);

			//Crossover for both circumference and square
			Path childCircumference = Crossover(parent1Circumference, parent2Circumference, rnd);
			Path childSquare = Crossover(parent1Square, parent2Square, rnd);

			//Mutation for both circumference and square
			if (rnd.Rannyu() < SwapMutationRate) {
				swapMutation(childCircumference, rnd);
			}
			if (rnd.Rannyu() < SwapMutationRate) {
				swapMutation(childSquare, rnd);
			}
			if (rnd.Rannyu() < PairPermutationRate) {
				pairPermutationMutation(childCircumference, rnd);
			}
			if (rnd.Rannyu() < PairPermutationRate) {
				pairPermutationMutation(childSquare, rnd);
			}
			if (rnd.Rannyu() < ShiftMutationRate) {
				shiftMutation(childCircumference, 2, 1, rnd);
			}
			if (rnd.Rannyu() < ShiftMutationRate) {
				shiftMutation(childSquare, 2, 1, rnd);
			}
			if (rnd.Rannyu() < PermutationRate) {
				permutationMutation(childCircumference, 2, rnd);
			}
			if (rnd.Rannyu() < PermutationRate) {
				permutationMutation(childSquare, 2, rnd);
			}
			if (rnd.Rannyu() < InversionRate) {
				inversionMutation(childCircumference, 3, rnd);
			}
			if (rnd.Rannyu() < InversionRate) {
				inversionMutation(childSquare, 3, rnd);
			}
			//The offspring replace the old population, and their fitness is evaluated again
			if (offspringCircumference.size() < populationSize) {
				offspringCircumference.push_back(childCircumference);
			}
			if (offspringSquare.size() < populationSize) {
				offspringSquare.push_back(childSquare);
			}
		}

		//Replace the old population with the offspring for both circumference and square
		populationCircumference = offspringCircumference;
		populationSquare = offspringSquare;

		//Evaluate the fitness of the new population for both circumference and square
		EvFitness(populationCircumference, cityListCircumference);
		EvFitness(populationSquare, cityListSquare);

		//Check if paths satisfy the bonds for both circumference and square
		for (Path& path : populationCircumference) {
			if (!checkFunction(path)) {
				cerr << "Invalid path found for circumference scenario!" << std::endl;
				return 1;
			}
		}
		for (Path& path : populationSquare) {
			if (!checkFunction(path)) {
				cerr << "Invalid path found for square scenario!" << std::endl;
				return 1;
			}
		}
	}

	//Get the best path from the final population for both circumference and square
	Path bestPathCircumference = getBestPath(populationCircumference);
	Path bestPathSquare = getBestPath(populationSquare);

	//Write data to files for both circumference and square
	writeDataToFile(avgPathLengthsCircumference, "L_values_circumference.dat");
	writeDataToFile(avgPathLengthsSquare, "L_values_square.dat");

	writeBestPathToFile(bestPathCircumference, cityListCircumference, "best_path_circumference.dat");
	writeBestPathToFile(bestPathSquare, cityListSquare, "best_path_square.dat");

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
