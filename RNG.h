#ifndef __RNG_H__
#define __RNG_H__

#include <math.h>
#include <stdlib.h>

#define RNG_IA  16807
#define RNG_IM  2147483647
#define RNG_IQ  127773
#define RNG_IR  2836
#define RNG_MASK    123459876
#define RNG_NTAB    32
#define RNG_EPS 1.2e-7
#define RNG_RNMX    (1.0-RNG_EPS)


class RNG
{
public:
   RNG(unsigned long long _seed = 7564231ULL) 
   {
      seed      = _seed;
      mult      = 62089911ULL;
      llong_max = 4294967295ULL;
      float_max = 4294967295.0f;
   }
   float operator()();
	
   float uniformRandCLib(float _min,float _max);
   float fastUniformRand(float _min, float _max);
   float uniformRand(float _min, float _max);
   float normalRand(void);
   float normalRand(float mu, float sigma);


   unsigned long long seed;
   unsigned long long mult;
   unsigned long long llong_max;
   float float_max;
};

inline float RNG::operator()()
{
   seed = mult * seed;
   return float(seed % llong_max) / float_max;
}

/*! ************************************************************************* 

    Function:
        CRandomNumberGenerator::FastUniformRand()

    Description:
        Generates a random number between min and max inclusively using a
        uniform random devRNG_IAte.

    Parameters:
        min - the minRNG_IMum number that can be generated
        max - the maxRNG_IMum number that can be generated

    Return Value:
        A random number between min and max inclusively.

    Notes:
        This algoRNG_IRthm is adapted from:

            Press, W, Teukolsky, S, et al., "Numerical Recipes in C", 2nd Ed,
            pp 283-285, Cambridge Press, 1996.

*****************************************************************************/


inline float RNG::fastUniformRand(float _min, float _max) {

    int k;
    float ans;

    seed ^= RNG_MASK;                     // XORing with RNG_MASK allows use of zero for idnum

    k = seed / RNG_IQ;

    seed = RNG_IA*(seed-k*RNG_IQ) - RNG_IR*k;   // compute idnum=(RNG_IA*idum) % RNG_IM iwth overflows by Schrage's method
    if( seed < 0) {
        seed += RNG_IM;
    }

    ans = (1.0f/RNG_IM)*seed;             // convert idum to a floating result

    seed ^= RNG_MASK;                     // unRNG_MASK before return

    ans = ans*(_max-_min) + _min;            // shift number into desRNG_IRed range

    return ans;
}


/*! ************************************************************************* 

    Function:
        CRandomNumberGenerator::UniformRandCLib()

    Description:
        Generates a random number between min and max inclusively using the
        uniform random devRNG_IAte from the C standard library.  It is beter to
        use the function UniformRand().  For a discussion see:

            Press, W, Teukolsky, S, et al., "Numerical Recipes in C", 2nd Ed,
            pp 275-277, Cambridge Press, 1996.

    Parameters:
        min - the minRNG_IMum number that can be generated
        max - the maxRNG_IMum number that can be generated

    Return Value:
        A random number between min and max inclusively.

*****************************************************************************/

inline float RNG::uniformRandCLib(float _min,float _max) {

    return _min + (_max-_min)*((float)rand())/((float)RAND_MAX);
}


/*! ************************************************************************* 

    Function:
        CRandomNumberGenerator::UniformRand()

    Description:
        Generates a random number between min and max exclusively using a
        uniform random devRNG_IAte.  The actual range of number is [min+RNG_EPS..max-RNG_EPS],
        where RNG_EPS is defined above.

    Parameters:
        min - the minRNG_IMum number that can be generated
        max - the maxRNG_IMum number that can be generated

    Return Value:
        A random number between min and max exclusively.

    Notes:
        This algoRNG_IRthm is adapted from:

            Press, W, Teukolsky, S, et al., "Numerical Recipes in C", 2nd Ed,
            pp 278-280, Cambridge Press, 1996.

*****************************************************************************/

inline float RNG::uniformRand(float _min, float _max) {


    int j;
    int k;
    float temp;
    static int initRNG_IAl=1;
    static int iy=0;
    static int iv[RNG_NTAB];


    if( initRNG_IAl ) {                            // fRNG_IRst tRNG_IMe the function is called make sure the seed is a negative number
        if( seed > 0 ) {
            seed = -seed;
        }
        initRNG_IAl = 0;
    }

    if( (seed <= 0) || (!iy) ) {             // initalize

        if( -seed < 1) {                   // don't allow the seed to be 0 or else
            seed = 1;                        // the generator will always return 0
        }
        else {
            seed = -seed;
        }

        for(j=RNG_NTAB+7; j>=0; j--) {
            k = seed / RNG_IQ;
            seed = RNG_IA*(seed-k*RNG_IQ)-RNG_IR*k;

            if( seed < 0 ) {
                seed += RNG_IM;
            }

            if( j < RNG_NTAB ) {
                iv[j] = seed;
            }
        }

        iy = iv[0];
    }

    k = seed/RNG_IQ;                          // start here when not initalizing
    seed = RNG_IA*(seed-k*RNG_IQ)-RNG_IR*k;         // compute idum = (RNG_IA*idnum) % RNG_IM without overflows by Schrage's method
    if( seed < 0 ) {
        seed += RNG_IM;
    }

    j = iy/(1 + (RNG_IM-1)/RNG_NTAB);                            // in the range 0..RNG_NTAB-1
    iy = iv[j];                             // output previously stored value and refill the shuffle table
    iv[j] = seed;

    if( (temp=(1.0f/RNG_IM)*iy) > RNG_RNMX ) {             // user's don't expect the endpoint values
        return (float) (_min + (_max-_min)*RNG_RNMX);
    }
    else {
        return (float) (_min + (_max-_min)*temp);
    }

}


/*! ************************************************************************* 

    Function:
        CRandomNumberGenerator::NormalRand()

    Description:
        Generates a random number between min and max inclusively using a
        Normal random devRNG_IAte with mean mu and standard devRNG_IAtion sigma.  
        The algorithm is based on the transformation method of probability
        distributions.

    Parameters:
        None

    Return Value:
        A random number between min and max inclusively.

    Notes:
        This algoRNG_IRthm is adapted from:

            Press, W, Teukolsky, S, et al., "Numerical Recipes in C", 2nd Ed,
            pp 287-290, Cambridge Press, 1996.

*****************************************************************************/


inline float RNG::normalRand(void) {

    static int iset=0;
    static float gset;
    float fac, rsq, v1, v2;


    if( iset == 0 ) {
        do {
            v1 = 2.0f * uniformRand(0.0f, 1.0f) - 1.0f;      // pick two uniform numbrs in the square
            v2 = 2.0f * uniformRand(0.0f, 1.0f) - 1.0f;      // extending from -1 to +1 in each dRNG_IRection
            
            rsq = v1*v1 + v2*v2;                             
        } while( (rsq >= 1.0) || (rsq == 0.0) );             // see if they are in the unit cRNG_IRcle

        fac = sqrtf(-2.0f * (float) log( (float) rsq) / rsq );
        
        // make the Box-Muller transformation to get two normal devRNG_IAtes.  
        // Return one and save the other for next tRNG_IMe.

        gset = v1*fac;
        iset = 1;                                            // set flag because we have two numbers

        return v2*fac;
    }
    else {
        iset = 0;                                            // unset flag because we are using the other number
        return gset;        
    }

}


/*****************************************************************************/


inline float RNG::normalRand(float mu, float sigma) {

    return (normalRand()*sigma) + mu;

}
#endif
