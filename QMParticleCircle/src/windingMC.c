/*
   Monte Carlo simulation to sample the winding number of a 1D periodic lattice.

   The simulation uses the Metropolis algorithm to update lattice sites and runs
   multiple replicas of the system (if selected) at different simbeta values.
   To avoid topological freezing, configurations are occasionally exchanged to
   reduce correlation.

   The output is a sequence of winding numbers from the input simbeta replica.

   The HIGH_TEMPERATURE mode activates a reweighting potential that boosts the
   probability of sampling configurations with larger winding numbers
   (while keeping their magnitude bounded by QMAX).

   This mechanism improves statistical access to the tails of the
   winding number distribution in the high temperature regime, relevant for
   the interested observable, allowing rare topological sectors to be
   explored with a reasonable number of samples.

   When HIGH_TEMPERATURE is enabled, a post processing reweighting analysis is
   required to correctly restore the statistical weight of the sampled
   configurations.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "../include/random.h"

//#define DEBUG_MODE
//#define HIGH_TEMPERATURE

#define STRING_LENGTH 50

#define QMAX 5
#define CHI 1

typedef struct simulation
{
    double *lattice;
    long int windingNumber;

    double simbeta;
    long int nSites;

    double eta;

    int replicaID;

} simulation;

/* Limit winding number within the reweighting range */
inline long int clampToRange(long int x, long int bound)
{
    if (x > bound)
        return bound;
    if (x < -bound)
        return -bound;
    return x;
}

/* Compute the reweighting potential added in the HIGH_TEMPERATURE mode */
double reweightPotential(simulation *restrict system)
{
    double addedPotential;
    long int effQ = clampToRange(system->windingNumber, QMAX);

    addedPotential = (double)(effQ * effQ);
    addedPotential /= (2.0 * system->simbeta * CHI);

    return addedPotential;
}

/* Compute the minimal distance between two points on a circle */
double circleDistance(double y, double x)
{
    double diff = y - x;

    if (fabs(diff) <= 0.5)
        return diff;
    else if (diff > 0.5)
        return diff - 1.;
    else
        return diff + 1.;
}

/* Compute the winding number of the lattice configuration */
long int computeWindingNumber(simulation const *restrict system,
                              long int const *restrict nnp)
{
    double res = 0;
    for (long int i = 0; i < system->nSites; i++)
    {
        res += circleDistance(system->lattice[nnp[i]], system->lattice[i]);
    }
    // Approximate to the nearest integer
    return lround(res);
}

/* Compute the total energy of a lattice configuration */
double computeEnergy(simulation *restrict system, long int const *restrict nnp)
{
    double distance, totalEnergy = 0;
    for (long int i = 0; i < system->nSites; i++)
    {
        distance = circleDistance(system->lattice[nnp[i]], system->lattice[i]);
        totalEnergy += distance * distance;
    }
    totalEnergy *= 1. / (2. * system->eta);

    #ifdef HIGH_TEMPERATURE
    totalEnergy -= reweightPotential(system);
    #endif

    return totalEnergy;
}

/* Compute the change in winding number resulting from a local site trial update */
long int deltaWinding(const simulation *restrict system,
                      long int const *restrict nnp,
                      long int const *restrict nnm,
                      long int site, double trial)
{
    double distanceOld, distanceNew, delta;

    // Change from backward bond
    distanceOld = circleDistance(system->lattice[site], system->lattice[nnm[site]]);
    distanceNew = circleDistance(trial, system->lattice[nnm[site]]);

    delta = distanceNew - distanceOld;

    // Change from forward bond
    distanceOld = circleDistance(system->lattice[nnp[site]], system->lattice[site]);
    distanceNew = circleDistance(system->lattice[nnp[site]], trial);

    delta += distanceNew - distanceOld;

    // Round to get an integer output
    return lround(delta);
}

/* Compute the variation of the reweighting potential when the winding number changes */
double deltaReweightPotential(const simulation *restrict system, long int deltaQ)
{
    double delta;
    long int effQold, effQnew;

    if (!deltaQ)
        return 0;

    effQold = clampToRange(system->windingNumber, QMAX);
    effQnew = clampToRange((system->windingNumber + deltaQ), QMAX);

    if (effQold == effQnew)
        return 0;

    delta = (double)(effQnew * effQnew - effQold * effQold) / (2. * system->simbeta * CHI);

    return delta;
}

/* Compute the local contribution to the energy of a single site */
static inline double localEnergy(const simulation *restrict system, double siteValue,
                                 long int site, long int const *restrict nnp,
                                 long int const *restrict nnm)
{
    double distance, bondEnergy;

    // Backward bond
    distance = circleDistance(siteValue, system->lattice[nnm[site]]);
    bondEnergy = distance * distance;

    // Forward bond
    distance = circleDistance(system->lattice[nnp[site]], siteValue);
    bondEnergy += distance * distance;

    bondEnergy /= (2.0 * system->eta);

    return bondEnergy;
}

/* Perform a single Metropolis sweep over the lattice */
long int metroSweep(simulation *restrict system, long int const *const restrict nnp,
                    long int const *const restrict nnm)
{
    const double step = 0.5; 
    double trial, new, old, deltaE;

    long int accepted = 0, deltaQ;

    // for each site of lattice, in sequential order, do:
    for (long int site = 0; site < system->nSites; site++)
    {
        // Compute local energy for current state lattice
        old = localEnergy(system, system->lattice[site], site, nnp, nnm);

        // Propose random trial update
        trial = system->lattice[site] + (2 * myrand() - 1) * step;
        while (trial >= 1.0)
            trial -= 1.0;

        while (trial < 0.0)
            trial += 1.0;

        // Compute local energy for proposed trial
        new = localEnergy(system, trial, site, nnp, nnm);

        // energy variation
        deltaE = new - old;

        #ifdef HIGH_TEMPERATURE
        // Include reweighting correction
        deltaQ = deltaWinding(system, nnp, nnm, site, trial);
        deltaE -= deltaReweightPotential(system, deltaQ);
        #endif

        if (deltaE <= 0.0 || myrand() < exp(-deltaE))
        {
            system->lattice[site] = trial;
            accepted += 1;

            #ifdef HIGH_TEMPERATURE
            system->windingNumber = lround(system->windingNumber + deltaQ);
            #endif
        }
    }
    return accepted;
}

/* Initialize a lattice configuration */
void initializeSystem(simulation *system, long int nSites, double simbeta, int index,
                      long int const *restrict nnp)
{
    system->lattice = (double *)malloc((unsigned long int)nSites * sizeof(double));
    if (system->lattice == NULL)
    {
        fprintf(stderr, "allocation problem at (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    system->simbeta = simbeta;
    system->nSites = nSites;
    system->eta = simbeta / (double)nSites;
    system->replicaID = index;

    for (long int ii = 0; ii < nSites; ii++)
    {
        system->lattice[ii] = 0.5;
        #ifdef HIGH_TEMPERATURE
        system->lattice[ii] = myrand();
        #endif
    }
    #ifdef HIGH_TEMPERATURE
    // compute winding number
    system->windingNumber = computeWindingNumber(system, nnp);
    #endif
}

/* Attempt a replica exchange between two systems*/
int attemptReplicaExchange(simulation *restrict system1, simulation *restrict system2,
                           long int const *restrict nnp)
{
    double Eold1, Eold2, Enew1, Enew2, deltaE;
    double *tmp;

    // Compute total energy of both replicas
    Eold1 = computeEnergy(system1, nnp);
    Eold2 = computeEnergy(system2, nnp);

    // Compute energies after the exchange
    Enew1 = Eold2 * system2->eta / system1->eta;
    Enew2 = Eold1 * system1->eta / system2->eta;

    // energy variation
    deltaE = (Enew1 - Eold1) + (Enew2 - Eold2);

    if (deltaE <= 0.0 || myrand() < exp(-deltaE))
    {
        // Swap lattice configurations
        tmp = system1->lattice;
        system1->lattice = system2->lattice;
        system2->lattice = tmp;

        #ifdef HIGH_TEMPERATURE
        // Swap winding numbers
        long int tmpQ;
        tmpQ = system1->windingNumber;
        system1->windingNumber = system2->windingNumber;
        system2->windingNumber = tmpQ;
        #endif

        #ifdef DEBUG_MODE
        // Swap replica identifiers for tracking
        int tmpreplicaID;
        tmpreplicaID = system1->replicaID;
        system1->replicaID = system2->replicaID;
        system2->replicaID = tmpreplicaID;
        #endif

        return 1;
    }
    else
        return 0;
}
int main(int argc, char **argv)
{
    if (argc != 8)
    {

        fprintf(stdout, "How to use this program:\n");
        fprintf(stdout, "  %s simbeta Nt sample repnumber simbetamax therm datafile\n\n", argv[0]);
        fprintf(stdout, "  simbeta     : hbar^2/(4*pi^2*m*R^2*k_B*T)\n");
        fprintf(stdout, "  Nt          : number of temporal steps\n");
        fprintf(stdout, "  sample      : number of Monte Carlo samples\n");
        fprintf(stdout, "  repnumber   : number of Parallel Tempering  replicas\n");
        fprintf(stdout, "  simbetamax  : max value of simbeta in the parallel tempering\n");
        fprintf(stdout, "  therm       : number of thermalization steps\n");
        fprintf(stdout, "  datafile    : name of the output file\n\n");

        return EXIT_SUCCESS;
    }

    // init random number generator
    const unsigned long int seed1 = (unsigned long int)time(NULL);
    const unsigned long int seed2 = seed1 + 127;

    myrand_init(seed1, seed2);

    // define variables
    int repnumber;

    #ifdef HIGH_TEMPERATURE
    const int measevery = 10, swapevery = 20;
    #else
    const int measevery = 5, swapevery = 10;
    #endif

    char filename[STRING_LENGTH];
    double simbeta, simbetamax, repbeta, aux;
    long int sample, therm, Nt, *nnp, *nnm,
        acceptedUpdate, acceptedExchange, res;

    // check input
    simbeta = atof(argv[1]);
    if (simbeta <= 0)
    {
        fprintf(stderr,
                "simbeta must be a positive number"
                "(%s, %d)\n",
                __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    Nt = atol(argv[2]);
    if (Nt <= 1)
    {
        fprintf(stderr,
                "Nt must be greater than 1"
                "(%s, %d)\n",
                __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    sample = atol(argv[3]);
    if (sample <= 1)
    {
        fprintf(stderr,
                "sample must be greater than 1"
                "(%s, %d)\n",
                __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    repnumber = atoi(argv[4]);
    if (repnumber < 1)
    {
        fprintf(stderr,
                "repnumber must be positive and at least 1"
                "(%s, %d)\n",
                __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    simbetamax = atof(argv[5]);
    if (simbetamax < simbeta)
    {
        fprintf(stderr,
                "simbetamax must be at least equal to simbeta"
                "(%s, %d)\n",
                __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    therm = atol(argv[6]);
    if (therm > sample || therm < 0)
    {
        fprintf(stderr,
                "therm must be positive and smaller than sample"
                "(%s, %d)\n",
                __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    if (strlen(argv[7]) >= STRING_LENGTH)
    {
        fprintf(stderr,
                "File name too long(%s, %d)\n",
                __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    else
    {
        strcpy(filename, argv[7]);
    }

    // init data structures
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "problem opening file (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

    nnp = (long int *)malloc((unsigned long int)Nt * sizeof(long int));
    nnm = (long int *)malloc((unsigned long int)Nt * sizeof(long int));

    if (nnp == NULL)
    {
        fprintf(stderr, "nnp allocation problem at (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    if (nnm == NULL)
    {
        fprintf(stderr, "nnm allocation problem at (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    // set periodic boundary conditions
    for (long int i = 0; i < Nt; i++)
    {
        nnp[i] = i + 1;
        nnm[i] = i - 1;
    }
    nnp[Nt - 1] = 0;
    nnm[0] = Nt - 1;

    // Initialize repnumber replicas
    simulation *replica = (simulation *)malloc((unsigned int)repnumber * sizeof(simulation));
    if (replica == NULL)
    {
        fprintf(stderr, "Config allocation problem at (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    for (int n = 0; n < repnumber; n++)
    {
        if (repnumber == 1)
            repbeta = simbeta;
        else
        {
            // simbeta ladder to get same swapping probabilities among all replicas
            repbeta = simbeta * pow(simbetamax / simbeta, n / (double)(repnumber - 1));
        }
        initializeSystem(&replica[n], Nt, repbeta, n, nnp);
    }

    // thermalization
    for (long int iter = 0; iter < therm; iter++)
    {
        for (int n = 0; n < repnumber; n++)
            metroSweep(&replica[n], nnp, nnm);

        if (iter % swapevery == 0)
        {
            if (myrand() >= 0.5)
            {
                for (int n = 0; n < repnumber - 1; n++)
                {
                    (void)attemptReplicaExchange(&replica[n], &replica[n + 1], nnp);
                }
            }
            else
            {
                for (int n = repnumber - 1; n > 0; n--)
                {
                    (void)attemptReplicaExchange(&replica[n], &replica[n - 1], nnp);
                }
            }
        }
    }
    // metropolis update and exchange
    acceptedUpdate = 0;
    acceptedExchange = 0;

    for (long int iter = 0; iter < sample; iter++)
    {
        // Update each replica
        for (int n = 0; n < repnumber; n++)
        {
            if (n == 0)
            {
                acceptedUpdate += metroSweep(&replica[n], nnp, nnm);
            }
            else
            {
                (void)metroSweep(&replica[n], nnp, nnm);
            }
        }
        // Attempt replica exchanges every swapevery steps
        if (iter % swapevery == 0)
        {
            // Select swapping direction randomly
            if (myrand() >= 0.5)
            {
                // Forward direction: lower to higher simbeta
                for (int n = 0; n < repnumber - 1; n++)
                {
                    acceptedExchange += attemptReplicaExchange(&replica[n], &replica[n + 1], nnp);
                }
            }
            else
            {
                // Backward direction: higher to lower simbeta
                for (int n = repnumber - 1; n > 0; n--)
                {
                    acceptedExchange += attemptReplicaExchange(&replica[n], &replica[n - 1], nnp);
                }
            }
            #ifdef DEBUG_MODE
            for (int n = 0; n < repnumber; n++)
            {
                printf("%i", replica[n].replicaID);
            }
            printf("\n");
            #endif
        }
        // Measure winding number every measevery steps
        if (iter % measevery == 0)
        {
            res = computeWindingNumber(&replica[0], nnp);
            fprintf(fp, "%li \n", res);
        }
    }
    // print stats
    aux = (double)acceptedUpdate / (double)sample / (double)Nt;
    fprintf(stderr, "acceptance rate: %lf\n", aux);
    if (repnumber != 1)
    {
        aux = (double)acceptedExchange / ((double)sample * (repnumber - 1) / (double)swapevery);
        fprintf(stderr, "swapping rate: %lf\n", aux);
    }
    for (int n = 0; n < repnumber; n++)
    {
        free(replica[n].lattice);
    }

    fclose(fp);
    free(replica);
    free(nnp);
    free(nnm);
}
