// Disease modeling in populations of badgers - spatial model
// Bram Kuijper, Maria Wellbelove
// 2019
//
// roughly after the model by 
// White & Harris 1995 Phil Trans B
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

enum Sex {
    Female = 0,
    Male = 1
};

enum Age {
    Cub = 0,
    Yearling = 1,
    Adult = 2,
};

enum State {
    Susceptible = 0,
    Latent = 1,
    Infectious = 2
};

// define a cell in the grid
struct GridCell {

    // inhabitants in each patch
    // (2 Sexes) x (3 Age Classes) x (3 Infectious States) x (max N individuals)
    Individual inhabitants[2][3][3][max_number_inhabitants]; // latent infected
   
    // number of individuals in each category
    // (2 Sexes) x (3 Age Classes) x (3 Infectious States) 
    int N[2][3][3]; // latent number
};


// define the actual population of grid cells
GridCell Population[grid_width][grid_height];


struct Individual {
    
    int infected_status 
}

// function or subroutine for intergroup transmission
void intergroup_transmission()
{


}

void dispersal()
{
    for (int column_i = 0; column_i < grid_width; ++column_i)
    {
        for (int row_j = 0; row_j < grid_height; ++row_j)
        {

}

void mortality()
{
    for (int column_i = 0; column_i < grid_width; ++column_i)
    {
        for (int row_j = 0; row_j < grid_height; ++row_j)
        {
            // first calculate number of individuals for the mortality rate
            int Nindtot = 0;
            int NadF = 0;
            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
                {
                    for (int age_i = 0; age_i < 3; ++age_i)
                    {
                        Nindtot += Population[column_i][row_j].N[sex_i][inf_state_i][age_i];

                        if (age_i > 1 && sex_i == Female)
                        {
                            NadF += Population[column_i][row_j].N[sex_i][inf_state_i][age_i];
                        }
                    }
                }
            }

            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
                {
                    for (int age_i = 1; age_i < 3; ++age_i)
                    {
                        Nind = Population[column_i][row_j].N[sex_i][inf_state_i][age_i];
                        

                        mort_prob = mortality_rate;
                        
                        // now perform mortality on cubs
                        // more complicated
                        //
                        // 1. if there are no adult females in Jan - June 
                        // cub mortality rate = 90%
                        // 2. if group size = 2 * carrying capacity
                        // cub mortality is 70%
                        // 3. semiquadritic function on top of this all

                        if (age_i == 0)
                        {
                            mort_prob = mortality_rate;

                            if (NadF == 0)
                            {
                                mort_prob = mortality_rate_cubs_alone;
                            }
                            else if (Nindtot >= 2 * K)
                            {
                                mort_prob = mortality_rate_cubs_highK;
                            }
                        }

                        for (int individual_i = 0; individual_i < Nind: ++individual_i)
                        {
                            // individual is sampled as one that dies
                            //
                            // death in yearlings and beyond according to fixed mort rate
                            if (uniform(rng_r) < mort_prob)
                            {
                                // delete this individual
                                Population[column_i][row_j][sex_i].inhabitants[
                                    inf_state_i][age_i][individual_i] = 
                                Population[column_i][row_j][sex_i].inhabitants[
                                    inf_state_i][age_i][Nind - 1];
                            
                            --individual_i;
                            --Nind;
                            }
                        } // for (int i = 0; i < n_dead; ++i)
                    } // for (int age_i = 1; age_i < 3; ++age_i)
                }
            }
        }
    }
} // end mortality()


// the reproduction function
void reproduce(int const time)
{
    // go through columns of the grid
    // and make cubs
    //
    // at the same time, cubs from previous year become yearlings
    //

    // auxiliary variable counting females in the local patch
    int n_females_total = 0;
    for (int column_i = 0; column_i < grid_width; ++column_i)
    {
        for (int row_j = 0; row_j < grid_height; ++row_j)
        {
            n_females_total = 0; 

            // sum over total number of females
            for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
            {
                for (int age_i = 1; age_i < 3; ++age_i)
                {
                    n_females_total += Population[column_i][width_i].N[Female][inf_state_i][age_i];
                }
            }

            // first, for all sexes, states and ages,
            // move cubs to yearlings, yearlings to adults
            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
                {
                    // loop over the last two different ages
                    // but in reverse fashion, so that first the adults are
                    // supplemented with yearlings
                    //
                    // then yearlings with cubs
                    for (int age = 2; age >= 1; --age)
                    {


                        nyoung = Population[column_i][width_i].N[sex_i][age-1][inf_state_i];
                        nold = Population[column_i][width_i].N[sex_i][age][inf_state_i];

                        // bounds checking
                        assert(nyoung >= 0);
                        assert(nyoung <= max_number_inhabitants);

                        assert(nold + nyoung >= 0);
                        assert(nold + nyoung <= max_number_inhabitants);

                        for (int young_i = 0; young_i < nyoung; ++young_i)
                        {
                            // copy individual from age - 1 to class age
                            Population[column_i][width_i].inhabitants[sex_i][age][inf_state_i][
                                Population[column_i][width_i].N[sex_i][age][inf_state_i]++
                            ] = Population[column_i][width_i].inhabitants[sex_i][age-1][inf_state_i][young_i];
                        }

                        // all individuals from age-1 copied to age
                        // hence set counter of all age-1 individuals to 0
                        Population[column_i][width_i].N[sex_i][age-1][inf_state_i] = 0;
                               
                        // bounds checking again (trust nothing)
                        assert(Population[column_i][width_i].N[sex_i][age][inf_state_i] <= 
                                max_number_inhabitants);

                    } // end for int age
                } // end for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
            } // end for (int sex_i = 0; sex_i < 2; ++sex_i)

            // now reproduce
            ncubs = 

        }// end for (int row_j = 0; row_j < grid_height; ++row_j)
    } // end for (int column_i = 0; column_i < grid_width; ++column_i)
} // end  void reproduce(int const time)

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



