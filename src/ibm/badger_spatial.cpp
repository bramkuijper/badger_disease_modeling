// Cooperative breeding in variable environments
// Bram Kuijper 
// 2019
//
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>

#include <random>


// various functions, such as unique filename creation
#include "auxiliary.h"
#include "individual.h"

#define DEBUG

// function which clips values between min and max
double clamp(double val, double min, double max)
{
    return(val > max ? max : val < min ? min : val);
}

// standard namespace
using namespace std;

// set up the random number generator
// set random seed etc
int seed = get_nanoseconds();
mt19937 rng_r{static_cast<long unsigned int>(seed)};
uniform_real_distribution<> uniform(0.0,1.0);

// parameters & variables:

// number of individuals in population
int grid_width = 12;
int grid_height = 12;

// maximum number of inhabitants in a patch
int max_number_inhabitants = 10;


// define a cell in the grid
struct GridCell {

    // make a list of the inhabitants
    Individual inhabitants_L[max_number_inhabitants]; // latent infected
    Individual inhabitants_I[max_number_inhabitants]; // infectious
    Individual inhabitants_S[max_number_inhabitants]; // susceptible
   
    // number of individuals in each category
    int number_of_inhabitants_L; // latent number
    int number_of_inhabitants_I; // infectious
    int number_of_inhabitants_S; // susceptible


};

struct Individual {
    
    int infected_status 

}

// function or subroutine for intergroup transmission
void intergroup_transmission()
{


}

void dispersal()
{

}

void deaths()
{

}

void reproduce()
{

}

// the key part of the code
// accepting command line arguments
int main(int argc, char **argv)
{
    // subroutine to initialize all the parameter values
    init_arguments(argc, argv);
    
    // create output files with unique filenames, 
    // preventing overwiting by other simulations
    string filename = "sim_badger";
    create_filename(filename);
    ofstream DataFile(filename);
    
    // function that writes headers to my data file
    write_data_headers(DataFile);

    // initialize the population (give all the badgers
    // some values)
    init_population();

    // the key part of the code
    for (int generation = 0; generation < number_generations; ++generation)
    {
        reproduce();

        intergroup_transmission();
        
        dispersal();

        deaths();

        if (generation % skip == 0)
        {
            write_stats(generation, DataFile);
        }
    } //j5

    write_parameters(DataFile);
}



