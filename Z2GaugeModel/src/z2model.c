/* Monte Carlo simulation of a Z2 lattice gauge theory in dim dimensions.

   Configurations are generated via ordered Metropolis sweeps on a dim 
   dimensional lattice with lexicographic indeces. Update attempts are 
   skipped with probability 'epsilon' to enforce aperiodicity and Boltzmann
   factors are precomputed in a lookup table for computational efficiency.

   Wilson loops W(Wt, Ws) are measured and averaged over all lattice sites
   and spatial directions. When MULTIHIT is enabled, link variables entering
   the loop are evaluated through  a variance-reduced link average, except
   for the first link of each segment of the path. For loops with Wt = 1 or
   Ws = 1, the multihit procedure is entirely disabled.

   The program outputs Wilson loop averages together with the Metropolis
   acceptance rate.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../include/random.h"
#include "../include/geometry.h"

#define MULTIHIT

#define dim 3
#define epsilon 0.05

#define STRING_LENGTH 50
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

const int maxStaples = 2 * (dim - 1);
const int dStaplesMax = 2 * maxStaples;

/* Initialize the lattice configuration.

   Allocates memory for a 'dim' dimensional lattice of linear size 'size'.
   Each site contains 'dim' links initialized randomly to ±1.

   The resulting lattice is stored in a dynamically allocated array:
   lattice[lex][dir], where 'lex' is the lexicographic site index and
   'dir' labels the direction.
*/
void init_lattice(int ***lattice, long int volume)
{
    int link;

    *lattice = (int **)malloc((unsigned long int)volume * sizeof(int *));
    if (!*lattice)
    {
        fprintf(stderr, "malloc failed for lattice");
        exit(EXIT_FAILURE);
    }
    for (long int i = 0; i < volume; i++)
    {
        (*lattice)[i] = (int *)malloc((unsigned int)dim * sizeof(int));
        if (!(*lattice)[i])
        {
            fprintf(stderr, "malloc failed for lattice[%li]", i);
            exit(EXIT_FAILURE);
        }
        for (int ii = 0; ii < dim; ii++)
        {
            link = 2 * (int)(2 * myrand()) - 1;
            (*lattice)[i][ii] = link;
        }
    }
}

/* Initialize nearest neighbour lookup tables.

   Allocates and fills the arrays 'nnp' and 'nnm', which store the forward
   and backward nearest neighbour indices for each lattice site.

   These tables are used to navigate the lattice efficiently in lexicographic
   representation.
*/
void init_neighbours(long int **nnp, long int **nnm, int size)
{
    long int volume = 1;
    for (int i = 0; i < dim; i++)
    {
        volume *= size;
    }
    (*nnp) = (long int *)malloc((unsigned long int)(dim * volume) * sizeof(long int));
    if (!(*nnp))
    {
        fprintf(stderr, "failed to initialize neighbours at (%s, %d)", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    (*nnm) = (long int *)malloc((unsigned long int)(dim * volume) * sizeof(long int));
    if (!(*nnm))
    {
        fprintf(stderr, "failed to initialize neighbours at (%s, %d)", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    init_geo(*nnp, *nnm, size, dim);
}

/* Initialize the exponential look-up table.

   Precomputes exp(beta * ΔStaples) for all possible variations in the
   Metropolis update  and stores them in 'boltzmannFactor' for faster access.

   The table indices are shifted by 'dStaplesMax' to allow negative deltas.

   Extra entries are allocated in the table (even if some are unused) to
   allow direct indexing from ΔStaples, avoiding the need to compute an
   index at runtime to improve efficiency.
*/
void init_expTable(double **boltzmannFactor, double beta)
{
    int range = 2 * dStaplesMax + 1;

    *boltzmannFactor = (double *)malloc((unsigned int)range * sizeof(double));
    for (int dStaples = -dStaplesMax; dStaples <= dStaplesMax; dStaples += 4)
    {
        (*boltzmannFactor)[dStaples + dStaplesMax] = exp(beta * dStaples);
    }
}

/* Initialize the exact multihit expectation table.

   For a Z2 link the conditional expectation value given the surrounding
   configuration admits the closed form tanh(beta*sumStaples). 
   
   This is precomputed for every achievable value of sumStaples
   and used by the Multihit algorithm in place of a stochastic estimate.
*/
void init_multihitMean(double **multihitMean, double beta)
{
    int range = 2 * maxStaples + 1;

    *multihitMean = (double *)malloc((unsigned int)range * sizeof(double));
    for (int staples = -maxStaples; staples <= maxStaples; staples += 2)
    {
        (*multihitMean)[staples + maxStaples] = tanh(beta * staples);
    }
}

/* Compute the sum of staples around a selected link.

   For the link at site 'lex' in direction 'dir' this function returns the sum
   of all staples of that link where each staple consists of a product of three
   neighbouring links forming the sides of the plaquette.
*/
int computeStaple(int **restrict Lattice,
                  long int const *const restrict nnp,
                  long int const *const restrict nnm,
                  long int lex, int dir, long int volume)
{
    int linkProd, sumStaples = 0;
    long int lex_minus_orth, lex_plus_dir, lex_plus_orth;

    /*                      ^ dir
                            |
                       lex_plus_dir
                   +--------+--------+
                   |        |        |
                   |        |        |
                   |        |        |    orth
                ---+--------+--------+--->
         lex_minus_orth    lex      lex_plus_orth
                            |
     */

    lex_plus_dir = nnp[dirgeo(lex, dir, volume)];

    for (int orth = 0; orth < dim; orth++)
    {
        if (orth == dir)
            continue;

        lex_plus_orth = nnp[dirgeo(lex, orth, volume)];
        lex_minus_orth = nnm[dirgeo(lex, orth, volume)];

        // ---------- forward staple ----------
        linkProd = Lattice[lex][orth];
        linkProd *= Lattice[lex_plus_dir][orth];
        linkProd *= Lattice[lex_plus_orth][dir];

        sumStaples += linkProd;

        // ---------- backward staple ----------
        linkProd = Lattice[lex_minus_orth][dir];
        linkProd *= Lattice[lex_minus_orth][orth];
        linkProd *= Lattice[nnp[dirgeo(lex_minus_orth, dir, volume)]][orth];

        sumStaples += linkProd;
    }
    return sumStaples;
}
/* Single Metropolis step for a selected link variable.

   Given the current link value and the corresponding staple sum, the
   function evaluates the local action variation and returns the updated
   flipped link according to the acceptance rule:

    - If ΔS > 0, the flip is always accepted.
    - If ΔS ≤ 0, the flip is accepted with probability exp(ΔS), looked up from the
      Boltzmann table previously defined.
*/
static inline int stepMetro(int link, int sumStaples,
     double const *const restrict boltzmannFactor)
{
    int dStaples = -2 * link * sumStaples;
    // to avoid out of range index in look-up Table
    if (dStaples < -dStaplesMax || dStaples > dStaplesMax)
    {
        fprintf(stderr, "Warning: out of range dStaples\n");
        return link;
    }
    if (dStaples > 0 || myrand() < boltzmannFactor[dStaples + dStaplesMax])
        return -link;

    return link;
}
/* Perform a full Metropolis sweep over the lattice.

   A single Metropolis update is applied to all links through an ordered
   sweep over sites and directions, where each update attempt is skipped with
   probability 'epsilon' to enforce aperiodicity.

   During the sweep the numbers of attempted and accepted updates are
   accumulated, and the corresponding acceptance rate is returned.
*/
double sweepMetro(int **restrict Lattice,
                  double const *const restrict boltzmannFactor,
                  long int const *const restrict nnp,
                  long int const *const restrict nnm,
                  long int volume)
{
    long int updateCounter = 0, trialCounter = 0;
    int sumStaples, oldLink, newLink;

    // ordered update of all lattice sites
    for (long int lex = 0; lex < volume; lex++)
    {
        // update in all possible direction for each site
        for (int dir = 0; dir < dim; dir++)
        {
            // probabilistic skip to enforce aperiodicity
            if (myrand() < epsilon)
            {
                continue;
            }
            else
            { // count total attempted link updates
                trialCounter += 1;

                // compute staples
                sumStaples = computeStaple(Lattice, nnp, nnm, lex, dir, volume);

                // count accepted updates
                oldLink = Lattice[lex][dir];
                newLink = stepMetro(oldLink,sumStaples,boltzmannFactor);

                if (oldLink != newLink)
                {
                    updateCounter += 1;

                    // update lattice
                    Lattice[lex][dir] = newLink;
                }
            }
        }
    }
    return (double)updateCounter / (double)(trialCounter);
}

/* Compute the multihit estimate of a link variable.

   The function performs a fixed number of successive Metropolis update
   attempts for the link at site 'lex' and direction 'dir', keeping the
   surrounding configuration fixed.

   The returned value is the average  over the resulting link values.
*/
double multihit(int **restrict Lattice,
                double const *const restrict multihitMean,
                long int lex, long int volume, int dir,
                long int const *const restrict nnp,
                long int const *const restrict nnm)
{
    double mean;

    int sumStaples = computeStaple(Lattice, nnp, nnm, lex, dir, volume);
    mean = multihitMean[sumStaples+maxStaples];
    
    return mean; 
}
/* Return the link value to be used in Wilson loop measurements.

   Depending on whether MULTIHIT is enabled, the function returns either
   the bare link variable Lattice[lex][dir] or its multihit estimate.

   The argument 'firstStep' forces the first link along each segment of
   the path to be evaluated without multihit, while 'disableMH' disables
   the multihit estimator entirely for the current Wilson loop.
*/
static inline double getLink(int **Lattice, long int lex, 
                             long int volume, int dir,
                             long int const *nnp, long int const *nnm,
                             int firstStep, int disableMH, 
                             double const *const restrict multihitMean)
{
    #ifdef MULTIHIT
        if (firstStep || disableMH)
            return (double)Lattice[lex][dir];
        return multihit(Lattice, multihitMean, lex, volume, dir, nnp, nnm);
    #else
        return (double)Lattice[lex][dir];
    #endif
}
/* Compute the Wilson loop W(Wt, Ws).

   The Wilson loop is defined as the product of link variables along a closed
   rectangular path in the (0, dir) plane, with temporal extent Wt and spatial extent Ws.

   Although a single loop would suffice in principle, the function averages over
   all lattice sites and spatial directions to exploit lattice symmetries and to
   improve statistics.
*/
double WilsonLoop(int **Lattice,
                  double const *const restrict multihitMean,
                  long int const *const restrict nnp,
                  long int const *const restrict nnm,
                  long int volume, int Wt, int Ws)
{
    double res = 0.0;
    long int check_lex;

    // disable MultiHit for entire Wilson loop
    int disableMH = (Wt == 1 || Ws == 1);

    // loop over all lattice sites
    for (long int lex = 0; lex < volume; lex++)
    {
        check_lex = lex;
        // loop over all spatial directions
        for (int dir = 1; dir < dim; dir++)
        {
            /*
             ^ 0-temporal direction
             | r1      r2
             +--------+
             |        |
             |        |
             |        | dir-spatial direction
          ---+--------+----->
             r0        r3
            */

            // initialize wilson loop variable
            double loop = 1.0;

            // starting point r0
            for (int i = 0; i < Wt; i++)
            {
                loop *= getLink(Lattice, lex, volume,
                     0, nnp, nnm, (i == 0),disableMH,multihitMean);
                lex = nnp[dirgeo(lex, 0, volume)];
            }
            // starting point r1
            for (int i = 0; i < Ws; i++)
            {
                loop *= getLink(Lattice, lex, volume,
                     dir, nnp, nnm, (i == 0),disableMH,multihitMean);
                lex = nnp[dirgeo(lex, dir, volume)];
            }
            // starting point r2
            for (int i = 0; i < Wt; i++)
            {
                lex = nnm[dirgeo(lex, 0, volume)];
                loop *= getLink(Lattice, lex, volume,
                     0, nnp, nnm, (i == 0),disableMH,multihitMean);
            }
            // starting point r3
            for (int i = 0; i < Ws; i++)
            {
                lex = nnm[dirgeo(lex, dir, volume)];
                loop *= getLink(Lattice, lex, volume,
                     dir, nnp, nnm, (i == 0),disableMH,multihitMean);
            }
            // check if loop has correctly been closed
            if (lex != check_lex)
            {
                fprintf(stderr, "Warning: loop not closed (lex=%ld, dir=%d)\n",
                        check_lex, dir);
                lex = check_lex;
                return NAN;
            }
            // at r0 again, loop closed
            res += (double)loop;
        }
    }
    // compute average
    res /= (double)volume;
    res /= (double)(dim - 1);

    return res;
}
int main(int argc, char **argv)
{
    FILE *fp;
    int **Lattice;
    long int *nnp, *nnm;
    double *boltzmannFactor;
    double *multihitMean;
    char datafile[STRING_LENGTH];

    int size, measevery;
    long int sample, therm;
    double beta, loop, accRate = 0.;

    const unsigned long int seed1 = (unsigned long int)time(NULL);
    const unsigned long int seed2 = seed1 + 127;

    if (argc != 7)
    {
        fprintf(stdout, "How to use this program:\n");
        fprintf(stdout, "  %s size beta sample therm measevery datafile\n\n", argv[0]);
        fprintf(stdout, "  size = temporal and spatial size of the lattice\n");
        fprintf(stdout, "  (space-time dimension defined by macro dim)\n");
        fprintf(stdout, "  beta = coupling\n");
        fprintf(stdout, "  sample = number of drawn to be extracted\n");
        fprintf(stdout, "  therm = number of thermalization sweeps\n");
        fprintf(stdout, "  measevery = MC sweep updates between consecutive measurements\n");
        fprintf(stdout, "  datafile = name of the file on which to write the data\n\n");
        fprintf(stdout, "Compiled for:\n");
        fprintf(stdout, "  dimensionality = %d\n\n", dim);
        fprintf(stdout, "Output:\n");
        fprintf(stdout, "  Wilson loops (Wt, Ws) using the following order\n");
        fprintf(stdout, "  for(Ws=1; Ws<=size/4; Ws++){\n");
        fprintf(stdout, "      for(Wt=1; Wt<=MIN(size/4,8); Wt++){\n");
        fprintf(stdout, "         Ws, Wt Wilson loop }}\n");

        return EXIT_SUCCESS;
    }
    else
    {
        // read and check input values
        size = atoi(argv[1]);
        beta = atof(argv[2]);
        sample = atol(argv[3]);
        therm = atol(argv[4]);
        measevery = atoi(argv[5]);

        if (strlen(argv[6]) >= STRING_LENGTH)
        {
            fprintf(stderr, "File name too long. Increse STRING_LENGTH or shorten the name (%s, %d)\n", __FILE__, __LINE__);
            return EXIT_FAILURE;
        }
        else
        {
            strcpy(datafile, argv[6]);
        }
        if (sample <= 0)
        {
            fprintf(stderr, "sample must be positive\n");
            return EXIT_FAILURE;
        }
        if (therm <= 0 || therm >= sample)
        {
            fprintf(stderr, "therm must be positive\n");
            return EXIT_FAILURE;
        }
        if (size <= 8)
        {
            fprintf(stderr, "'size' must be at least 8 to get useful data\n");
            return EXIT_FAILURE;
        }
    }
    // to get n°sample measures
    sample *= measevery;

    // initialize random number generator
    myrand_init(seed1, seed2);

    // compute lattice space time volume
    long int volume = 1;
    for (int dir = 0; dir < dim; dir++)
    {
        volume *= size;
    }

    // initialize data structures
    init_expTable(&boltzmannFactor, beta);
    init_multihitMean(&multihitMean, beta);
    init_lattice(&Lattice, volume);
    init_neighbours(&nnp, &nnm, size);


    // open file
    fp = fopen(datafile, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    // thermalization
    for (long int iter = 0; iter < therm; iter++)
    {
        sweepMetro(Lattice, boltzmannFactor, nnp, nnm, volume);
    }
    for (long int iter = 0; iter < sample; iter++)
    {
        accRate += sweepMetro(Lattice, boltzmannFactor, nnp, nnm, volume);
        if (iter % measevery == 0)
        {
            for (int Ws = 1; Ws <= MIN(size / 4, 7); Ws++)
            {
                for (int Wt = 1; Wt <= MIN(size / 4, 7); Wt++)
                {

                    loop = WilsonLoop(Lattice, multihitMean,
                                      nnp, nnm, volume, Wt, Ws);

                    fprintf(fp, "%.25f ", loop);
                }
            }
            fprintf(fp, "\n");
        }
    }
    accRate /= (double)sample;

    printf("\n");
    printf("Acceptance rate: %f\n", accRate);

    // close datafile
    fclose(fp);

    // free memory
    for (long int i = 0; i < volume; i++)
    {
        free(Lattice[i]);
    }
    free(Lattice);
    free(nnp);
    free(nnm);
    free(boltzmannFactor);
    free(multihitMean);

    return EXIT_SUCCESS;
}