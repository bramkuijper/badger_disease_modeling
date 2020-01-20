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
#include <utility>

#include <random>


// various functions, such as unique filename creation
#include "auxiliary.hpp"
#include "individual.hpp"

//#define DEBUG

// function which clips values between (min and max)
void clamp(double &val, double min, double max)
{
    val = val > max ? max : val < min ? min : val;
}

// standard namespace
using namespace std;

// set up the random number generator
// set random seed etc
int seed =get_nanoseconds();
mt19937 rng_r{static_cast<unsigned int>(seed)};
uniform_real_distribution<> uniform(0.0,1.0);

// random direction of dispersal to neighbouring cells
uniform_int_distribution<> dispersal_direction(0, 3);

int skip = 1;
int space_skip = 1;

// parameters & variables:

// number of individuals in population
int const grid_width = 12;
int const grid_height = 12;

uniform_int_distribution<int> random_column(0, grid_width - 1);
uniform_int_distribution<int> random_row(0, grid_height - 1);

// maximum number of inhabitants in a patch
int max_number_inhabitants = 10;
int max_generations = 5;

// baseline mortality rate
double mortality_rate = 0.2;

// mortality rate when cubs are without adult female
double mortality_rate_cubs_alone = 0.9;
double mortality_rate_cubs_highK = 0.9;

double K = 5;
double tau = 1.0;
double phi = 0.5;

double mu = 0.01;
double sdmu = 0.01;

normal_distribution<> mutational_effect(0.0, sdmu);


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

// dispersal probability per year for females and males respectively
double dispersal_probs[2] = { 0.02,0.06 };

double init_v = 0.0;

// define a cell in the grid
struct GridCell {

    // inhabitants in each patch
    // (2 Sexes) x (3 Age Classes) x (3 Infectious States) x (max N individuals)
    vector <Individual> inhabitants[2][3][3];
    
    // inhabitants in each patch
    // (2 Sexes) x (3 Age Classes) x (3 Infectious States) x (max N individuals)
    vector <Individual> immigrants[2][3][3]; // latent infected
};


// define the actual population of grid cells
GridCell Population[grid_width][grid_height];


// initialize population
void init_population()
{
    // generate a standard individual
    Individual standard_individual(0.0);

    for (int column_i = 0; column_i < grid_width; ++column_i)
    {
        for (int row_j = 0; row_j < grid_height; ++row_j)
        {
            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
                {
                    // only make yearlings or older
                    for (int age_i = 1; age_i < 3; ++age_i)
                    {
                        for (int individual_i = 0; 
                                individual_i < K/2; 
                                ++individual_i)
                        {
                            // add new individual to this particular stack
                            Population[column_i][row_j].
                                inhabitants[sex_i][age_i][inf_state_i].
                                    push_back(standard_individual);
                        }
                    }
                }
            }
        }
    }
} // end void init_population()

// write statistics to file datafile
//
// parameters:
// -----------
//  int const generation: 
//      the current generation
//
//  ofstream &datafile:
//      reference to a ofstream object
//      to which the data will be written
void write_statistics(
        int const timestep,
        ofstream &datafile)
{
    // collect total numbers of individuals
    // in all classes
    int N[2][3][3];
    int Nss[2][3][3];

    // calculate virulence
    double meanv = 0.0;
    double ssv = 0.0;

    int Ntotal = 0;
    int Ntotal_infected = 0;

    // initialize counters to 0
    for (int sex_i = 0; sex_i < 2; ++sex_i)
    {
        for (int age_i = 0; age_i < 3; ++age_i)
        {
            for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
            {
                N[sex_i][age_i][inf_state_i] = 0;
                Nss[sex_i][age_i][inf_state_i] = 0;
            }
        }
    }

    // some temporary storage variables
    int the_size;
    double v;

    for (int column_i = 0; column_i < grid_width; ++column_i)
    {
        for (int row_j = 0; row_j < grid_height; ++row_j)
        {
            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int age_i = 0; age_i < 3; ++age_i)
                {
                    for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
                    {
                        the_size = Population[column_i][row_j].
                            inhabitants[sex_i][age_i][inf_state_i].size();

                        N[sex_i][age_i][inf_state_i] += the_size;
                        Ntotal += the_size;
                        
                        Nss[sex_i][age_i][inf_state_i] += the_size * the_size;

                        if ((State) inf_state_i != Susceptible)
                        {
                            Ntotal_infected += the_size;

                            for (int individual_i = 0; individual_i < the_size; ++individual_i)
                            {
                                v = Population[column_i][row_j].
                                    inhabitants[sex_i][age_i][inf_state_i][individual_i].v;

                                meanv += v;
                                ssv += v*v;
                            }
                        }
                    }
                }
            }
        }
    }

    int meanN;

    datafile << timestep << ";";

    for (int sex_i = 0; sex_i < 2; ++sex_i)
    {
        for (int age_i = 0; age_i < 3; ++age_i)
        {
            for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
            {
                meanN = N[sex_i][age_i][inf_state_i]/(grid_width * grid_height);

                datafile << N[sex_i][age_i][inf_state_i] << ";"
                            << meanN << ";"
                            << Nss[sex_i][age_i][inf_state_i]/(grid_width * grid_height) 
                                - meanN * meanN << ";";

            }
        }
    }

    meanv /= Ntotal_infected;
    ssv /= Ntotal_infected;

    datafile << meanv << ";"
            << ssv - meanv * meanv << ";" 
            << Ntotal << ";"
            << Ntotal_infected << ";"
            << endl;
} // write_statistics()

// write_spatial_snapshot
//
// writes a snapshot of the current spatial configuration
//
// parameters:
// -----------
//  int const generation: 
//      the current generation
//
//  ofstream &datafile:
//      reference to a ofstream object
//      to which the data will be written
void write_spatial_snapshot(
        int const timestep,
        ofstream &datafile)
{
    for (int column_i = 0; column_i < grid_width; ++column_i)
    {
        for (int row_j = 0; row_j < grid_height; ++row_j)
        {
            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
                {
                    // only make yearlings or older
                    for (int age_i = 1; age_i < 3; ++age_i)
                    {
                        datafile 
                            << timestep << ";"
                            << column_i << ";"
                            << row_j << ";"
                            << (sex_i == 0 ? "F" : "M") << ";"
                            << age_i << ";"
                            << inf_state_i << ";"
                            << Population[column_i][row_j].
                                inhabitants[sex_i][age_i][inf_state_i].size() << ";" << endl;
                    }
                }
            }
        }
    }
} // end void write_spatial_snapshot(



// write the headers of the file containing the data
//
// parameters:
// -----------
//  ofstream &datafile:
//      reference to an ofstream object
//      to which the 'normal' data will be written
//      (counts and averages)
//
//  ofstream &datafile_space
//      reference to an ofstream object
//      to which spatial snapshots will be written
//      see write_spatial_snapshot()
void write_data_headers(
        ofstream &datafile_counts
        ,ofstream &datafile_space)
{
    datafile_counts << "time;";
    string inf_state_chr[3] = {"S","L","I"};

    ostringstream oss;

    for (int sex_i = 0; sex_i < 2; ++sex_i)
    {
        for (int age_i = 0; age_i < 3; ++age_i)
        {
            for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
            {
                oss << (sex_i == 0 ? "F":"M") 
                    << age_i
                    << inf_state_chr[inf_state_i]; 

                datafile_counts 
                    << "N_" << oss.str() << ";"
                    << "meanN_" << oss.str() << ";"
                    << "varN_" << oss.str() << ";";
            }
        }
    }

    datafile_counts << "meanv;varv;Ntotal;Ntotal_infected;" << endl;

    datafile_space << 
        "time;column;row;sex;age;inf_state;N;" << endl;

} // end write_data_headers()

// get arguments from the command line
void init_arguments(int argc, char **argv)
{
    // command line parsing using getopt
    tau = atof(argv[1]);
    phi = atof(argv[2]);
    mortality_rate = atof(argv[3]);
    mortality_rate_cubs_alone = atof(argv[4]);
    mortality_rate_cubs_highK = atof(argv[5]);
    init_v = atof(argv[6]);
    K = atof(argv[7]);
    mu = atof(argv[8]);
    sdmu = atof(argv[9]);
}

// writes parameters used in the simulation to a file
//
// parameters:
// -----------
//  ofstream &datafile:
//      reference to a ofstream object
//      to which the data will be written
void write_parameters(
        ofstream &datafile)
{
    datafile
        << endl
        << endl // two line breaks to have a space between data and the parameters
        << "tau;" << tau << endl// the 
        << "phi;" << tau << endl;// the 
}

// overall product of contact rates and transmission probability
double beta(double const v)
{
    return(phi * v / (v + tau));
}

// function or subroutine for intragroup transmission
void intragroup_transmission()
{
    // we need to make a probability distribution of infection probabilities
    // where we have to remember the individual parasites
    //
    // hence we create a cumulative distribution of infection probabilities
    // draw a random number from this cumulative distribution 
    // and where the random number 'falls' within the distribution, this is 
    // the parasite that is going to infect matters 
    //
    // make a vector of pairs (first value cumulative value
    // second value the virulence trait of the parasite
    typedef pair <double, double> cumul_virulence_pair;
    vector < cumul_virulence_pair > prob_infect_cumul;

    // total cumulative value of infectivity
    double infectivity_cumul_total = 0.0;
    double random_cumul_sample = 0.0;

    // auxiliary variables
    double current_cumul_val; // current cumulative value
    double current_v; // current value of virulence
    double current_beta; // current value of virulence

    // auxiliary variable storing population sizes etc
    int Nind;
    
    int local_popsize;
    
    // first calculate infection probability
    for (int column_i = 0; column_i < grid_width; ++column_i)
    {
        for (int row_j = 0; row_j < grid_height; ++row_j)
        {
            // clear the vector containing (cumulative distribution, virulence) pairs
            prob_infect_cumul.clear();

            // make a cumulative distribution of per capita infectivity
            infectivity_cumul_total = 0.0;

            local_popsize = 0;
            
            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int age_i = 0; age_i < 3; ++age_i)
                {
                    for (int state_i = 0; state_i < 3; ++state_i)
                    {
                        local_popsize += Population[column_i][row_j].
                                            inhabitants[sex_i][age_i][state_i].size();
                    }
                }
            }


            // go over all infected individuals and make 
            // a per capita infectivity distribution
            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int age_i = 0; age_i < 3; ++age_i)
                {
                    Nind = Population[column_i][row_j].
                                inhabitants[sex_i][age_i][Infectious].size();

                    for (int individual_i = 0; individual_i < Nind; ++individual_i)
                    {
                        // remember current value of virulence
                        current_v = Population[column_i][row_j].
                                            inhabitants[sex_i][age_i][Infectious][individual_i].v;
                        assert(current_v >= 0.0);

                        // see Allen 1994 Math Biosci p 86
                        current_beta = beta(current_v) / local_popsize;

                        infectivity_cumul_total = current_cumul_val = 
                            infectivity_cumul_total + current_beta;

                        cumul_virulence_pair current_pair{current_cumul_val, current_v};
                        
                        prob_infect_cumul.push_back(current_pair);
                    }
                }
            } 
            
            // made cumulative distribution of infection
            // probabilities, now infect susceptibles

            assert(infectivity_cumul_total <= 1.0);

            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int age_i = 0; age_i < 3; ++age_i)
                {
                    // how many susceptibles in this category?
                    Nind = Population[column_i][row_j].
                                inhabitants[sex_i][age_i][Susceptible].size();

                    for (int individual_i = 0; individual_i < Nind; ++individual_i)
                    {
                        // this individual won't be infected, next
                        if (uniform(rng_r) > infectivity_cumul_total)
                        {
                            continue;
                        }
                        
                        // now who will infect this individual?
                        random_cumul_sample = uniform(rng_r) * infectivity_cumul_total;

                        // iterate over vector
                        for (vector<cumul_virulence_pair>::iterator it = 
                                    prob_infect_cumul.begin();
                                it != prob_infect_cumul.end();
                                ++it)
                        {
                            // found the individual?
                            if (random_cumul_sample <= it->first)
                            {
                                // transmit parasite's genetic virulence value
                                Population[column_i][row_j].
                                    inhabitants[sex_i][age_i][Susceptible][individual_i].v =
                                        it->second;

                                // add individual to stack of 
                                // infectious individuals
                                Population[column_i][row_j].
                                    inhabitants[sex_i][age_i][Infectious].push_back(
                                    Population[column_i][row_j].
                                        inhabitants[sex_i][age_i][Susceptible][individual_i]
                                    );

                                // now remove from susceptible stack
                                Population[column_i][row_j].
                                    inhabitants[sex_i][age_i][Susceptible].erase(
                                        Population[column_i][row_j].
                                            inhabitants[sex_i][age_i][Susceptible].begin()
                                            + individual_i
                                    );

                                --Nind;
                                --individual_i;
                                break; // break out of cumulative distribution for loop
                            } // end if(random_cumul_val)
                        } // end for (vector<T>::iterator it = prob_infect_cumul.begin();
                    } //end for int individual_i
                } // end for int age_i
            } //end for int sex_i
        } // end for (int row_j = 0; row_j < grid_height; ++row_j)
    } // end for (int column_i = 0; column_i < grid_width; ++column_i)
} // end void intragroup_transmission()

// dispersal of individuals among patches
void dispersal()
{
    // auxiliary variables to store destination location
    // for migrants
    int column_destination;
    int row_destination;

    // auxiliary variable to keep track of numbers 
    // of individuals in the current class
    int n_current_class;

    for (int column_i = 0; column_i < grid_width; ++column_i)
    {
        for (int row_j = 0; row_j < grid_height; ++row_j)
        {
            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
                {
                    for (int age_i = 0; age_i < 3; ++age_i)
                    {
                        n_current_class = Population[column_i][row_j].inhabitants[sex_i][age_i][inf_state_i].size();

                        for (int individual_i = 0;
                                individual_i < n_current_class; 
                                ++individual_i)
                        {
                            // see whether individual disperses or not
                            if (uniform(rng_r) < dispersal_probs[sex_i])
                            {
                                column_destination = column_i;
                                row_destination = row_j;

                                switch(dispersal_direction(rng_r))
                                {
                                    case 0: // north
                                        column_destination = column_i;
                                        row_destination = row_j + 1;
                                        break;
                                    case 1: // east
                                        column_destination = column_i + 1;
                                        row_destination = row_j;
                                        break;
                                    case 2: // south 
                                        column_destination = column_i;
                                        row_destination = row_j - 1;
                                        break;
                                    case 3: // west
                                        column_destination = column_i - 1;
                                        row_destination = row_j;
                                        break;
                                    default:
                                        cout << "direction fail" << endl;
                                        exit(1);
                                }

                                // Torus: 
                                // if coordinate negative
                                // move to top of grid (max_row/ max_col - 1)
                                //
                                // if coordinate exceeds rows/columns
                                // move to coordinate 0
                                if (row_destination < 0)
                                {
                                    row_destination = grid_height - 1;
                                }
                                else if (row_destination >= grid_height)
                                {
                                    row_destination = 0;
                                }

                                if (column_destination < 0)
                                {
                                    column_destination = grid_width - 1;
                                }
                                else if (column_destination >= grid_width)
                                {
                                    column_destination = 0;
                                }

                                assert(column_destination >= 0);
                                assert(column_destination < grid_width);
                                assert(row_destination >= 0);
                                assert(row_destination < grid_height);

                                // move the individual from original spot
                                // to new immigrant stack
                                Population[column_destination][row_destination].
                                    immigrants[sex_i][age_i][inf_state_i].push_back(
                                        Population[column_i][row_j].
                                            inhabitants[sex_i][age_i][inf_state_i][individual_i]
                                );

                                // remove individual from current habitat
                                Population[column_i][row_j].
                                    inhabitants[sex_i][age_i][inf_state_i].erase(
                                        Population[column_i][row_j].
                                            inhabitants[sex_i][age_i][inf_state_i].begin() 
                                            + individual_i
                                );

                                --n_current_class;
                                --individual_i;

                                assert(individual_i <= n_current_class);
                            } // end if if (uniform(rng_r) < dispersal[sex_i])

                        } // end for for (int individual_i = 0; individual_i < Nind; ++individual_i)

                    } // end for (int age_i = 0; age_i < 3; ++age_i)

                } // end for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)

            } // end for (int sex_i = 0; sex_i < 2; ++sex_i)

        } // end for (int row_j = 0; row_j < grid_height; ++row_j)

    } // end for (int column_i = 0; column_i < grid_width; ++column_i)
   
    // now add immigrants to inhabitants
    for (int column_i = 0; column_i < grid_width; ++column_i)
    {
        for (int row_j = 0; row_j < grid_height; ++row_j)
        {
            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
                {
                    for (int age_i = 0; age_i < 3; ++age_i)
                    {
                        // get total number of immigrants to add
                        int Nimm = Population[column_i][row_j].
                                    immigrants[sex_i][age_i][inf_state_i].size();

                        for (int immigrant_i = 0; immigrant_i < Nimm; ++immigrant_i)
                        {
                            // add arriving immigrants to the current population
                            // inhabitants
                            Population[column_i][row_j].
                                inhabitants[sex_i][age_i][inf_state_i].
                                    push_back(
                                            Population[column_i][row_j].
                                                immigrants[sex_i][age_i][inf_state_i][immigrant_i]
                            );
                        }

                        // all immigrants copied to the inhabitant stack, erase them
                        Population[column_i][row_j].
                            immigrants[sex_i][age_i][inf_state_i].clear();

                    } // end for (int age_i = 0; age_i < 3; ++age_i)

                } // end for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)

            } // end for (int sex_i = 0; sex_i < 2; ++sex_i)

        } // end for (int row_j = 0; row_j < grid_height; ++row_j)

    } // end for (int column_i = 0; column_i < grid_width; ++column_i)
} // end dispersal function
       

void mortality()
{
    // auxiliary variable to keep track of number
    // of individuals in current class
    int n_current_class;

    // auxiliary variable containing the mortality probabilty
    double mort_prob;

    double v;

    for (int column_i = 0; column_i < grid_width; ++column_i)
    {
        for (int row_j = 0; row_j < grid_height; ++row_j)
        {
            // first calculate number of individuals 
            // of all classes in local patch is the mortality rate
            // is density-dependent
            int Nindtot = 0;
            // moreover mortality rate is also dependent on number of adult females
            int NadF = 0;

            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
                {
                    for (int age_i = 0; age_i < 3; ++age_i)
                    {
                        Nindtot += Population[column_i][row_j].
                            inhabitants[sex_i][age_i][inf_state_i].size();

                        // calculate 
                        if (age_i > 1 && sex_i == Female)
                        {
                            NadF += Population[column_i][row_j].
                                inhabitants[sex_i][age_i][inf_state_i].size();
                        }
                    }
                }
            }

            for (int sex_i = 0; sex_i < 2; ++sex_i)
            {
                for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
                {
                    for (int age_i = 0; age_i < 3; ++age_i)
                    {
                        n_current_class = Population[column_i][row_j].
                            inhabitants[sex_i][age_i][inf_state_i].size();

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
                            if (NadF == 0)
                            {
                                mort_prob = mortality_rate_cubs_alone;
                            }
                            else if (Nindtot >= 2 * K)
                            {
                                mort_prob = mortality_rate_cubs_highK;
                            }
                        }

                        double mort_prob_i = 0.0;
                        for (int individual_i = 0; individual_i < n_current_class; ++individual_i)
                        {
                            mort_prob_i = mort_prob;

                            // if individual is infected the virulence of
                            // bovine TB will be important
                            if ((State) inf_state_i == Infectious)
                            {
                                v = Population[column_i][row_j].
                                        inhabitants[sex_i][age_i][inf_state_i][individual_i].v;

                                assert(v >= 0);

                                mort_prob_i += v;
                            }
                            // individual is sampled as one that dies
                            //
                            // death in yearlings and beyond according to fixed mort rate
                            if (uniform(rng_r) < mort_prob_i)
                            {
                                // delete this individual
                                Population[column_i][row_j].
                                    inhabitants[sex_i][age_i][inf_state_i].erase(
                                Population[column_i][row_j].
                                    inhabitants[sex_i][age_i][inf_state_i].begin() + individual_i
                                    );
                            
                                --individual_i;
                                --n_current_class;
                            }
                        } // for (int i = 0; i < n_dead; ++i)
                    } // for (int age_i = 1; age_i < 3; ++age_i)
                }
            }
        }
    }
} // end mortality()

// mutation function
void mutate(double &val)
{
    if (uniform(rng_r) < mu)
    {
        val += mutational_effect(rng_r);
    }
}

// create offspring
//
// Parameters:
//      Individual &mother  - reference to mother
//      Individual &father - reference to father
//      Individual &offspring - reference to offspring
//
void create_offspring(
        Individual &mother
        ,Individual &father
        ,Individual &offspring)
{
    offspring.v = uniform(rng_r) < 0.5 ? father.v : offspring.v;

    mutate(offspring.v);
    clamp(offspring.v, 0.0, 1.0);
}

// the reproduction function, see p393 White & Dodds
void reproduce()
{
    // go through columns of the grid
    // and make cubs
    //
    // at the same time, cubs from previous year become yearlings
    //

    // auxiliary variable counting individuals younger than the current
    // age class
    int nyoung = 0;

    // auxiliary individuals counting individuals in current age class
    int n_current_age= 0;

    // number of cubs to be produced as a floating point number
    double n_cubs_f = 0;
    // number of cubs to be produced as an integer
    double n_cubs_i = 0;

    // to randomly sample female we need to remember the number of
    // females present in 
    // age classes yearling and adult (cub females cannot
    // yet reproduce)
    int n_female_cumul_dist[2][3];
    int n_females_total = 0;

    int n_male_cumul_dist[2][3];
    int n_males_total = 0;

    // go over the different patches
    // and increase age of all individuals
    for (int column_i = 0; column_i < grid_width; ++column_i)
    {
        for (int row_j = 0; row_j < grid_height; ++row_j)
        {
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
                        nyoung = Population[column_i][row_j].
                            inhabitants[sex_i][age-1][inf_state_i].size();

                        n_current_age = Population[column_i][row_j].
                            inhabitants[sex_i][age][inf_state_i].size();

                        // bounds checking
                        assert(nyoung >= 0);
                        assert(nyoung <= max_number_inhabitants);

                        assert(n_current_age + nyoung >= 0);
                        assert(n_current_age + nyoung <= max_number_inhabitants);

                        for (int young_i = 0; young_i < nyoung; ++young_i)
                        {
                            assert(
                                    Population[column_i][row_j].
                                                inhabitants[sex_i][age-1][inf_state_i][young_i].v >= 0);
                            // copy individual from age - 1 to class age
                            Population[column_i][row_j].
                                inhabitants[sex_i][age][inf_state_i].
                                    push_back(
                                            Population[column_i][row_j].
                                                inhabitants[sex_i][age-1][inf_state_i][young_i]
                                            );
                            assert(
                                    Population[column_i][row_j].
                                    inhabitants[sex_i][age][inf_state_i].size() > 0);
                                    
                            assert(
                                    Population[column_i][row_j].
                                    inhabitants[sex_i][age][inf_state_i][
                                        Population[column_i][row_j].
                                        inhabitants[sex_i][age][inf_state_i].size() - 1
                                    ].v >= 0);

                        }


                        // all individuals from age-1 copied to age
                        // hence set counter of all age-1 individuals to 0
                        Population[column_i][row_j].
                            inhabitants[sex_i][age - 1][inf_state_i].clear();
                    } // end for int age
                } // end for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
            } // end for (int sex_i = 0; sex_i < 2; ++sex_i)
            
            n_females_total = 0; 
            n_males_total = 0; 

            // sum over total number of females
            for (int inf_state_i = 0; inf_state_i < 3; ++inf_state_i)
            {
                for (int age_i = 1; age_i < 3; ++age_i)
                {
                    assert(age_i - 1 >= 0);
                    n_female_cumul_dist[age_i - 1][inf_state_i] = n_females_total + 
                        Population[column_i][row_j].
                            inhabitants[Female][age_i][inf_state_i].size();

                    n_females_total = n_female_cumul_dist[age_i - 1][inf_state_i];
                    
                    n_male_cumul_dist[age_i - 1][inf_state_i] = n_males_total + 
                        Population[column_i][row_j].
                            inhabitants[Male][age_i][inf_state_i].size();

                    n_males_total = n_male_cumul_dist[age_i - 1][inf_state_i];
                }
            }
                        
            // now reproduce
            if (n_females_total < 1 || n_males_total < 1)
            {
                n_cubs_f = 0;
            }
            else
            {
                n_cubs_f = 0.60 + 0.63 * n_females_total;
            }
            
            // now round the number to an integer. 
            // normal rounding would just be using round() but 
            // this means that values of n_cubs_f of say, 1.4 
            // are always equal to 1. In reality, values of n_cubs_f of
            // 1.4 should be equal to 2 40% of the time and equal to 1
            // 60% of the time
            n_cubs_i = floor(n_cubs_f);

            if (uniform(rng_r) < n_cubs_f - n_cubs_i)
            {
                ++n_cubs_i;
            }

            // make random number generators for males and females
            uniform_int_distribution<> random_male_rng(0, n_males_total - 1);
            uniform_int_distribution<> random_female_rng(0, n_females_total - 1);

            // now make the offspring due to random mating in the patch
            for (int cub_i = 0; cub_i < n_cubs_i; ++cub_i)
            {
                // make variables to store the focal cub,
                // its mother and father
                Individual cub;
                Individual father;
                Individual mother;

                // get random number for mother
                // and father, find these numbers in the cumulative
                // distribution of males and females
                int random_male = random_male_rng(rng_r);
                int random_female = random_female_rng(rng_r);

                // indicator variables that are set to true 
                // when male and female are found so that search 
                // loop can be exited
                bool male_found = false;
                bool female_found = false;


                // start a search loop to find a random father and mother
                // we loop through a cumulative distribution of counts
                // of individuals. If our randomly chosen number lower
                // than the current count, pick this column
                //
                //
                // EXAMPLE
                // say, I have a total of 8 male badgers
                //
                // and three classes of male badgers (infected, susceptible, resistant)
                // containing, respectively 3 and 4 and 1 badger
                //
                // I generate a random number, say, random_male=5, which badger to choose?
                //
                // random_male > n_infected_males (as 5 > 3), moving on
                // random_male <= n_infected_males + n_susceptible_males (5 < 3 + 4; 
                // note cumulative count): stop & choose random susceptible male

                for (int infection_status_sample_i = 0; 
                        infection_status_sample_i < 3; 
                        ++infection_status_sample_i)
                {
                    for (int age_sample_i = 1; 
                            age_sample_i < 3; 
                            ++age_sample_i)
                    {
                        // is this male category 
                        // of infection_status_sample_i
                        // and age_sample_i the right one?
                        if (!male_found && 
                                random_male < 
                                    n_male_cumul_dist[age_sample_i - 1][infection_status_sample_i])
                        {
                            // check whether there are individuals here? 
                            assert(Population[column_i][row_j].
                                    inhabitants[Male][age_sample_i][infection_status_sample_i].size() > 0);
                            uniform_int_distribution<> male_sampler(0, 
                                    Population[column_i][row_j].
                                        inhabitants[Male][age_sample_i][infection_status_sample_i].size() - 1
                            );


                            father = Population[column_i][row_j].
                                inhabitants[Male][age_sample_i][infection_status_sample_i][
                                male_sampler(rng_r)
                                ];

                            male_found = true;
                        }

                        if (!female_found && random_female < n_female_cumul_dist[age_sample_i - 1][infection_status_sample_i])
                        {
                            assert(Population[column_i][row_j].
                                    inhabitants[Female][age_sample_i][infection_status_sample_i].size() > 0);
                            uniform_int_distribution<> female_sampler(0, 
                                    Population[column_i][row_j].
                                        inhabitants[Female][age_sample_i][infection_status_sample_i].size() - 1
                            );

                            mother = Population[column_i][row_j].
                                inhabitants[Female][age_sample_i][infection_status_sample_i][female_sampler(rng_r)];

                            female_found = true;
                        }

                        if (male_found && female_found)
                        {
                            break;
                        }
                    }
                    
                    if (male_found && female_found)
                    {
                        break;
                    }
                }

                assert(male_found);
                assert(female_found);

                assert(father.v >= 0.0);
                assert(mother.v >= 0.0);

                create_offspring(
                        mother
                        ,father
                        ,cub);

                assert(cub.v >= 0.0);

                // randomly determine sex of the cub
                Sex sex_cub = uniform(rng_r) < 0.5 ? Female : Male;

                // add cub to offsprign stack
                // all cubs are born as susceptible
                Population[column_i][row_j].
                    inhabitants[sex_cub][Cub][Susceptible].push_back(cub);
                            
            }
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
    ofstream SpatialFile(filename + "_space");
    
    // function that writes headers to my data file
    write_data_headers(DataFile, SpatialFile);

    // initialize the population (give all the badgers
    // some values)
    init_population();

    // the key part of the code
    for (int generation = 0; generation < max_generations; ++generation)
    {
        reproduce();

        intragroup_transmission();
        
        dispersal();

        mortality();

        if (generation % skip == 0)
        {
            write_statistics(generation, DataFile);
        }

        if (generation % space_skip == 0)
        {
            write_spatial_snapshot(generation, SpatialFile);
        }
    } 

    write_parameters(DataFile);
}



