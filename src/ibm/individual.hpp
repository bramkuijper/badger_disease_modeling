#ifndef __INDIVIDUAL_INCLUDED__
#define __INDIVIDUAL_INCLUDED__

class Individual
{
    public:
        double v; // virulence

    Individual();

    Individual(
        double const vval
    );
    
    Individual(Individual const &other);

    void operator=(Individual const &other);
};


#endif
