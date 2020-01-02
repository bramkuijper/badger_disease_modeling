#include "individual.hpp"


Individual::Individual():
    v{0.0}
{
}


Individual::Individual(
        double const vval
        ):
    v{vval}
{
}


// copy constructor
Individual::Individual(Individual const &other):
    v{other.v}
{
}

// overload the assignment operator 
void Individual::operator=(Individual const &other) 
{
    v = other.v;
}
