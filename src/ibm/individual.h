#ifndef __INDIVIDUAL_INCLUDED__
#define __INDIVIDUAL_INCLUDED__

class Individual
{
    public:
        // diploid loci for helping in the two envts
        double h[2][2];
        // diploid loci for dispersal in the two envts
        double d[2][2];

    Individual();

    Individual(
        double const hval1
        ,double const hval2
        ,double const dval1
        ,double const dval2
    );
    
    Individual(Individual const &other);

    void operator=(Individual const &other);
};


#endif
