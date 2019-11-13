#ifndef __INDIVIDUAL_INCLUDED__
#define __INDIVIDUAL_INCLUDED__

class Individual
{
    public:

        int t_birth = 0;

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
