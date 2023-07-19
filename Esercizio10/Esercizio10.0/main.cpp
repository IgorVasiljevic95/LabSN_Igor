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
#include <string>
#include "random.h"
#include <sstream>
#include <mpi.h>

using namespace std;

// Defining the City class with coordinates x, y;
class City {
public:
    string name;
    string description;
    double x, y;
};

// Defining the ComputingNode class with city index and temperature
class ComputingNode {
public:
    int cityIndex;
    double temperature;

    bool operator==(const ComputingNode& other) const {
        return cityIndex == other.cityIndex && temperature == other.temperature;
    }
};

// Defining the Path class which contains the computingNodes vector and fitness
class Path {
public:
    vector<ComputingNode> computingNodes;
    double fitness;

    void calculateFitness(const std::vector<City>& cityList) {
        fitness = 0.0;
        for (size_t i = 0; i < computingNodes.size() - 1; i++) {
            int cityIndex1 = computingNodes[i].cityIndex;
            int cityIndex2 = computingNodes[i + 1].cityIndex;
            double distance = calculateDistance(cityList[cityIndex1], cityList[cityIndex2]);
            fitness += distance;
        }
        // Add distance from the last city back to the first city
        int lastIndex = computingNodes[computingNodes.size() - 1].cityIndex;
        int firstIndex = computingNodes[0].cityIndex;
        double lastDistance = calculateDistance(cityList[lastIndex], cityList[firstIndex]);
        fitness += lastDistance;
    }

private:
    // Calculate the distance between cities
    double calculateDistance(const City& city1, const City& city2) {
        double dx = city1.x - city2.x;
        double dy = city1.y - city2.y;
        return sqrt(dx * dx + dy * dy);
    }
};

//Function to generate a random permutation of computing nodes
vector<ComputingNode> generateRandomPermutation(int size, const vector<double>& temperatures, Random rnd) {
    vector<ComputingNode> permutation(size);
    for (int i = 0; i < size; i++) {
        permutation[i].cityIndex = i;
        permutation[i].temperature = temperatures[i % temperatures.size()];
    }
		for(int i=size-1; i>0; --i) {
			swap(permutation[i],permutation[rnd.RannyuInt(0,i+1)]);
		}
		return permutation;
}

//Function to create the starting population
vector<Path> createInitialPopulation(int populationSize, int numCities, const vector<double>& temperatures, Random rnd) {
    vector<Path> population(populationSize);
    for (int i = 0; i < populationSize; i++) {
        population[i].computingNodes = generateRandomPermutation(numCities, temperatures, rnd);
    }
    return population;
}

//Function to check if a path satisfies the bond
bool checkPath(const Path& path) {
    vector<int> visited(path.computingNodes.size(), 0);
    for (const auto& node : path.computingNodes) {
        if (visited[node.cityIndex] == 1) {
            return false;
        }
        visited[node.cityIndex] = 1;
    }
    return true;
}

//Function to evaluate the fitness of the population
void evaluateFitness(vector<Path>& population, const vector<City>& cityList) {
    for (Path& path : population) {
        path.calculateFitness(cityList);
    }
}

//Function for tournament selection
Path tournamentSelection(const vector<Path>& population, int tournamentSize, Random& rnd) {
    int populationSize = population.size();
    vector<int> tournamentIndices(tournamentSize);
    for (int i = 0; i < tournamentSize; i++) {
        tournamentIndices[i] = rnd.Rannyu(0, populationSize);
    }
    auto minElement = min_element(tournamentIndices.begin(), tournamentIndices.end(), [&](int i, int j) {
        return population[i].fitness < population[j].fitness;
    });
    return population[*minElement];
}

//Function for swap mutation
void swapMutation(Path& path, Random& rnd) {
    int size = path.computingNodes.size();
    int pos1 = rnd.Rannyu(1, size - 1);
    int pos2 = rnd.Rannyu(1, size - 1);
    swap(path.computingNodes[pos1], path.computingNodes[pos2]);
}

//Function for pair permutation mutation
void pairPermutationMutation(Path& path, Random& rnd) {
    int size = path.computingNodes.size();
    int pos = rnd.Rannyu(1, size - 1);
    swap(path.computingNodes[pos], path.computingNodes[pos + 1]);
}

// Function for shift mutation
void shiftMutation(Path& path, int m, int n, Random& rnd) {
    int size = path.computingNodes.size();
    int startPos = rnd.Rannyu(1, size - m);
    rotate(path.computingNodes.begin() + startPos, path.computingNodes.begin() + startPos + n, path.computingNodes.begin() + startPos + m + 1);
}

// Function for permutation mutation
void permutationMutation(Path& path, int m, Random& rnd) {
    int size = path.computingNodes.size();
    int startPos1 = rnd.Rannyu(1, size - 2 * m);
    int startPos2 = startPos1 + m;
    reverse(path.computingNodes.begin() + startPos1, path.computingNodes.begin() + startPos1 + m);
    reverse(path.computingNodes.begin() + startPos2, path.computingNodes.begin() + startPos2 + m);
}

// Function for inversion mutation
void inversionMutation(Path& path, int m, Random& rnd) {
    int size = path.computingNodes.size();
    int startPos = rnd.Rannyu(1, size - m);
    reverse(path.computingNodes.begin() + startPos, path.computingNodes.begin() + startPos + m);
}

// Function for getting the best path from the population
Path getBestPath(const vector<Path>& population) {
    auto minElement = min_element(population.begin(), population.end(), [&](const Path& path1, const Path& path2) {
        return path1.fitness < path2.fitness;
    });
    return *minElement;
}

// Function to write path length and average path length to a file
void writeDataToFile(const std::vector<double>& pathLengths, const std::vector<double>& avgPathLengths, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (size_t i = 0; i < pathLengths.size(); i++) {
            file << i << " " << pathLengths[i] << " " << avgPathLengths[i] << std::endl;
        }
        file.close();
    }
}

// Function to write the best path in a file
void writeBestPathToFile(const Path& path, const vector<City>& cityList, const string& filename) {
    ofstream file(filename);
    if (file.is_open()) {
        for (const auto& node : path.computingNodes) {
            int cityIndex = node.cityIndex;
            file << cityList[cityIndex].x << " " << cityList[cityIndex].y << endl;
        }
        file.close();
    }
}

// Function for Simulated Annealing with many temperatures
void simulatedAnnealing(vector<Path>& population, const vector<City>& cityList, const vector<double>& temperatures, int numGenerations, int tournamentSize, double mutationRate, Random& rnd) {
    int populationSize = population.size();
    int numTemperatures = temperatures.size();

    // Store the best path and average path length for each generation and temperature
    vector<vector<double>> pathLengths(numTemperatures, vector<double>(numGenerations));
    vector<vector<double>> avgPathLengths(numTemperatures, vector<double>(numGenerations));

    // Main GA loop
    for (int generation = 0; generation < numGenerations; generation++) {
        // Sort the population based on fitness
        sort(population.begin(), population.end(), [&](const Path& path1, const Path& path2) {
            return path1.fitness < path2.fitness;
        });

        // Calculate average path length of the best half of the population for each temperature
        for (int t = 0; t < numTemperatures; t++) {
            double avgPathLength = 0.0;
            for (int i = 0; i < populationSize / 2; i++) {
                if (population[i].computingNodes[0].temperature == temperatures[t]) {
                    avgPathLength += population[i].fitness;
                }
            }
            avgPathLength /= (populationSize / 2);

            pathLengths[t][generation] = population[0].fitness;
            avgPathLengths[t][generation] = avgPathLength;
        }

        // Tournament selection, mutation, and exchange trial move
        vector<Path> offspring;
        while (offspring.size() < populationSize) {
            // Selection
            Path parent1 = tournamentSelection(population, tournamentSize, rnd);
            Path parent2 = tournamentSelection(population, tournamentSize, rnd);

            // Crossover (order crossover is not used in Simulated Annealing)
            Path child = parent1;

            // Mutation
            if (rnd.Rannyu() < mutationRate) {
                swapMutation(child, rnd);
            }
            if (rnd.Rannyu() < mutationRate) {
                pairPermutationMutation(child, rnd);
            }
            if (rnd.Rannyu() < mutationRate) {
                shiftMutation(child, 2, 1, rnd);
            }
            if (rnd.Rannyu() < mutationRate) {
                permutationMutation(child, 2, rnd);
            }
            if (rnd.Rannyu() < mutationRate) {
                inversionMutation(child, 3, rnd);
            }

            // Apply the exchange trial move
            int nodeIndex = rnd.Rannyu(0, child.computingNodes.size() - 1);
            int nextNodeIndex = nodeIndex + 1;
            if (child.computingNodes[nodeIndex].temperature != child.computingNodes[nextNodeIndex].temperature) {
                swap(child.computingNodes[nodeIndex], child.computingNodes[nextNodeIndex]);
            }

            // Evaluate the fitness of the new path
            child.calculateFitness(cityList);

            // Acceptance criterion based on the Metropolis algorithm
            double deltaFitness = child.fitness - parent1.fitness;
            double temperature = parent1.computingNodes[0].temperature;
            double acceptanceProb = exp(-deltaFitness / temperature);
            if (rnd.Rannyu() < acceptanceProb) {
                offspring.push_back(child);
            }
            else {
                offspring.push_back(parent1);
            }
        }

        // Replace the old population with the offspring
        population = offspring;

        // Evaluate the fitness of the new population
        evaluateFitness(population, cityList);

        // Check if paths satisfy the bonds
        for (Path& path : population) {
            if (!checkPath(path)) {
                cerr << "Invalid path found!" << std::endl;
                return;
            }
        }

        // Exchange trial move
        if (generation % 10 == 0) {
            int numExchanges = numTemperatures - 1;
            for (int i = 0; i < numExchanges; i++) {
                int index1 = rnd.Rannyu(0, populationSize);
                int index2 = index1 + 1;

                // Swap the paths between adjacent temperatures
                if (population[index1].computingNodes[0].temperature != population[index2].computingNodes[0].temperature) {
                    swap(population[index1], population[index2]);
                }
            }
        }
			// Find the best path and temperature with the lowest fitness
							double minFitness = std::numeric_limits<double>::max();
							for (const Path& path : population) {
									if (path.fitness < minFitness) {
											minFitness = path.fitness;
											bestPath = path;
									}
							}

							// Update L_values and average path lengths for each temperature
							for (int t = 0; t < numTemperatures; t++) {
									double avgPathLength = 0.0;
									for (int i = 0; i < populationSize / 2; i++) {
											if (population[i].computingNodes[0].temperature == temperatures[t]) {
													avgPathLength += population[i].fitness;
											}
									}
									avgPathLength /= (populationSize / 2);

									pathLengths[t][generation] = population[0].fitness;
									avgPathLengths[t][generation] = avgPathLength;
							}
					}
					

					

					// Find the lowest temperature index and corresponding best path
					int lowestTemperatureIndex = -1;
					for (int t = 0; t < numTemperatures; t++) {
							if (temperatures[t] == bestPath.computingNodes[0].temperature) {
									lowestTemperatureIndex = t;
									break;
							}
					}
				
					// Write L_values for all temperatures to a file
					
						stringstream ss;
						ss << "L_values_" << 0 << ".dat";
						string filename = ss.str();
						writeDataToFile(avgPathLengths[0], filename);
					

					// Write the best path for the lowest temperature to a file
					if (lowestTemperatureIndex != -1) {
							stringstream ss;
							ss << "best_path_lowest_temp.dat";
							string bestPathFilename = ss.str();
							writeBestPathToFile(bestPath, cityList, bestPathFilename);
					}
			}

    }

    // Get the best path from each temperature's final population
    vector<Path> bestPaths(numTemperatures);
    for (int t = 0; t < numTemperatures; t++) {
        auto minElement = std::min_element(population.begin(), population.end(), [&](const Path& path1, const Path& path2) {
            return (path1.computingNodes[0].temperature == temperatures[t]) && (path1.fitness < path2.fitness);
        });
        bestPaths[t] = *minElement;
    }

    // Write data to files for each temperature
    for (int t = 0; t < numTemperatures; t++) {
        stringstream ss;
        ss << "L_values_" << t << ".dat";
        string filename = ss.str();
        writeDataToFile(pathLengths[t], avgPathLengths[t], filename);

        stringstream ss2;
        ss2 << "best_path_" << t << ".dat";
        string bestPathFilename = ss2.str();
        writeBestPathToFile(bestPaths[t], cityList, bestPathFilename);
    }
}

// MAIN
int main(int argc, char** argv) {
    // Initializing MPI
    MPI_Init(&argc, &argv);

    // Getting the number of processes and the rank of the current process
    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Initializing Random rnd *********
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    }
    else {
        cerr << "PROBLEM: Unable to open Primes" << endl;
        return 1;
    }
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    }
    else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
        return 1;
    }
    //*************

    // Read city coordinates from the "American_capitals.dat" file
    ifstream cityFile("American_capitals.dat");
    if (!cityFile.is_open()) {
        cerr << "Failed to open the city coordinates file." << endl;
        return 1;
    }

    // Skip the first row (header) in the file
    string line;
    getline(cityFile, line);

    // Read the city coordinates from the file and store them in the cityList vector
    vector<City> cityList;
    while (getline(cityFile, line)) {
        City city;
        stringstream ss(line);
        ss >> city.name >> city.description >> city.x >> city.y;
        cityList.push_back(city);
    }
    cityFile.close();

    // Define parameters
    int numCities = cityList.size();
    int populationSize = 100;
    int tournamentSize = 5;
    int numGenerations = 100;
    double mutationRate = 0.1;

    // Calculate the number of temperatures per process
    int numTemperatures = numProcesses > 10 ? 10 : numProcesses;
    int numTemperaturesPerProcess = numTemperatures / numProcesses;
    int remainingTemperatures = numTemperatures % numProcesses;

    // Calculate the start and end indices of temperatures for the current process
    int startTemperatureIndex = rank * numTemperaturesPerProcess;
    int endTemperatureIndex = startTemperatureIndex + numTemperaturesPerProcess;

    // Adjust the end temperature index for the last process if there are remaining temperatures
    if (rank == numProcesses - 1 && remainingTemperatures > 0) {
        endTemperatureIndex += remainingTemperatures;
    }

    // Create the vector of temperatures for the current process
    vector<double> temperatures;
    for (int i = startTemperatureIndex; i < endTemperatureIndex; i++) {
        double temperature = 0.1 * pow(2.0, i);
        temperatures.push_back(temperature);
    }

    // Create the initial population
    vector<Path> population = createInitialPopulation(populationSize, numCities, temperatures, rnd);

    // Evaluate the fitness of the initial population
    evaluateFitness(population, cityList);

    // Apply the Simulated Annealing algorithm with multiple temperatures
    simulatedAnnealing(population, cityList, temperatures, numGenerations, tournamentSize, mutationRate, rnd);

    // Finalize MPI
    MPI_Finalize();
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
