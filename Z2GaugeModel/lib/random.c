#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include"../include/random.h"


uint32_t pcg32_random_r(pcg32_random_t *rng)
{
  uint64_t oldstate = rng->state;
  // Advance internal state
  rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
  // Calculate output function (XSH RR), uses old state for max ILP
  uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
  uint32_t rot = (uint32_t)(oldstate >> 59u);
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

void pcg32_srandom_r(pcg32_random_t *rng, uint64_t initstate, uint64_t initseq)
{
  rng->state = 0U;
  rng->inc = (initseq << 1u) | 1u;
  pcg32_random_r(rng);
  rng->state += initstate;
  pcg32_random_r(rng);
}

void myrand_init(unsigned long int initstate, unsigned long int initseq)
{
  pcg32_srandom_r(&pcg32_random_state, (uint64_t)initstate, (uint64_t)initseq);
}


double myrand(void)
{
  return (double)pcg32_random_r(&pcg32_random_state) / (pow(2.0, 32.0));
}


