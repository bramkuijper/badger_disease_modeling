#include "individual.h"


Individual::Individual():
    h{{0.0,0.0},{0.0,0.0}},
    d{{0.0,0.0},{0.0,0.0}}
{
}


Individual::Individual(
        double const hval1
        ,double const hval2
        ,double const dval1
        ,double const dval2):
    h{{hval1},{hval2}},
    d{{dval1},{dval2}}
{
}


// copy constructor
Individual::Individual(Individual const &other):
    h{{other.h[0][0],other.h[0][1]},{other.h[1][0],other.h[1][1]}},
    d{{other.d[0][0],other.d[0][1]},{other.d[1][0],other.d[1][1]}}
{
}

// overload the assignment operator 
void Individual::operator=(Individual const &other) 
{
    for (int envt_i = 0; envt_i < 2; ++envt_i)
    {
        for (int allele_i = 0; allele_i < 2; ++allele_i)
        {
            h[envt_i][allele_i] = other.h[envt_i][allele_i];
            d[envt_i][allele_i] = other.d[envt_i][allele_i];
        }
    }
}
