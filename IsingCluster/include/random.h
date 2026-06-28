#ifndef RANDOM_H
#define RANDOM_H

#include <math.h>
#include<stdint.h>

typedef struct { uint64_t state; uint64_t inc ;} pcg32_random_t;

typedef double(*function)(double);

uint32_t pcg32_random_r(pcg32_random_t *rng);

void pcg32_srandom_r(pcg32_random_t *rng, uint64_t initstate, uint64_t initseq);

void myrand_init(unsigned long int initstate, unsigned long int initseq);

double myrand(void);

pcg32_random_t pcg32_random_state;



#endif // end RANDOM_H