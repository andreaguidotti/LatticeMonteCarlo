/* 
   Monte Carlo simulation of the Ising model on a periodic dim-dimensional 
   lattice of volume size^dim.

   Nearest neighbours in each spatial direction are precomputed and stored in
   lookup tables (nnp, nnm) to allow fast navigation of the lattice with
   lexicographic ordering.

   Configurations are updated using the Wolff cluster algorithm: starting from
   a random seed site, aligned neighbouring spins are added to the cluster with
   probability p = 1 - exp(-2*beta). The completed cluster is then flipped.

   After a thermalization phase of 'therm' cluster updates, the program
   performs further updates and measures the energy and magnetization per site
   at fixed intervals. Results are written to the specified output file.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "../include/random.h"
#include "../include/geometry.h"

#define STRING_LENGTH 50

/*
   Data structure for the Ising model simulation.
   Contains everything needed for Wolff algorithm
*/
typedef struct
{
    int *restrict lattice; // spin configuration

    int dim;         // lattice dimension
    int size;        // lattice size
    long int volume; // total sites in lattice

    long int *restrict nnp; // next neighbour plus
    long int *restrict nnm; // next neighbour minus

    long int *restrict cluster; // where sites of cluster are stored
    int *restrict keepTrack;    // Look up table to check if site is already in the cluster

    long int cluster_index; // last index of the cluster

} ising;

/*
   Initialize the ising system:

   Allocates memory for lattice, cluster, keepTrack, and neighbor arrays.
   Initializes spins randomly.
*/
void init_system(ising *sys, int size, int dim)
{
    int initHadProblems = 0;

    sys->size = size;
    sys->dim = dim;

    sys->volume = 1;
    for (int d = 0; d < sys->dim; d++)
    {
        sys->volume *= sys->size;
    }

    sys->lattice = (int *)malloc((unsigned long int)sys->volume * sizeof(int));
    if (sys->lattice == NULL)
    {
        fprintf(stderr, "lattice initialization failed at (%s, %d)\n", __FILE__, __LINE__);
        initHadProblems = 1;
    }
    for (long int i = 0; i < sys->volume; i++)
    {
        sys->lattice[i] = 2 * (int)(2 * myrand()) - 1;
    }

    sys->keepTrack = (int *)malloc((unsigned long int)sys->volume * sizeof(int));
    if (sys->keepTrack == NULL)
    {
        fprintf(stderr, "keepTrack initialization failed at (%s, %d)\n", __FILE__, __LINE__);
        initHadProblems = 1;
    }

    sys->cluster = (long int *)malloc((unsigned long int)sys->volume * sizeof(long int));
    if (sys->cluster == NULL)
    {
        fprintf(stderr, "cluster initialization failed at (%s, %d)\n", __FILE__, __LINE__);
        initHadProblems = 1;
    }

    sys->nnp = (long int *)malloc((unsigned long int)(dim * sys->volume) * sizeof(long int));
    if (sys->nnp == NULL)
    {
        fprintf(stderr, "nnp initialization failed at (%s, %d)\n", __FILE__, __LINE__);
        initHadProblems = 1;
    }

    sys->nnm = (long int *)malloc((unsigned long int)(dim * sys->volume) * sizeof(long int));
    if (sys->nnm == NULL)
    {
        fprintf(stderr, "nnm initialization failed at (%s, %d)\n", __FILE__, __LINE__);
        initHadProblems = 1;
    }

    if (initHadProblems)
        exit(EXIT_FAILURE);
    init_geo(sys->nnp, sys->nnm, size, dim);
}

/* compute configuration energy per unit volume */

static inline double energy(ising const *const restrict sys)
{
    double E, aux;

    aux = 0.0;

    for (long int i = 0; i < sys->volume; i++)
    {
        for (int ii = 0; ii < sys->dim; ii++)
        {
            aux -= sys->lattice[i] * sys->lattice[sys->nnp[dirgeo(i, ii, sys->volume)]];
        }
    }

    E = aux / (double)sys->volume;

    return E;
}

/* compute average magnetization of the system */

static inline double magn(ising const *const restrict sys)
{
    double m;

    m = 0.0;

    for (long int idx = 0; idx < sys->volume; idx++)
    {
        m += sys->lattice[idx];
    }

    m /= (double)sys->volume;

    return m;
}
/*
    Initialize a new cluster.

    Selects a random starting spin on the lattice, initializes the cluster
    with that site and marks it as visited in 'keepTrack'.
*/
void init_cluster(ising *restrict sys)
{
    memset(sys->keepTrack, 0, (unsigned long int)sys->volume * sizeof(int));
    sys->cluster_index = 0;

    long int start = (long int)(myrand() * sys->volume);

    sys->cluster[0] = start;
    sys->keepTrack[start] = 1;
}

/*
   Check and possibly add neighboring sites to the cluster.

   Examines the  nearest neighbors of given site (=cluster[idx]).
   If a neighbor has the same spin orientation and passes the probabilistic
   acceptance test (p < prob_add), it is added to the cluster.

   Returns:
     Number of new sites added to the cluster
*/
static inline int check_near(ising *restrict sys, int idx, double prob_add)

{
    int addcount = 0;

    long int lex = sys->cluster[idx];
    long int neighbour_plus, neighbour_minus;

    for (int ii = 0; ii < sys->dim; ii++)
    {
        neighbour_plus = sys->nnp[dirgeo(lex, ii, sys->volume)];

        if (sys->keepTrack[neighbour_plus] == 0 &&
            sys->lattice[lex] * sys->lattice[neighbour_plus] == 1)
        {
            if (myrand() < prob_add)
            {
                addcount++;
                sys->cluster[sys->cluster_index + addcount] = neighbour_plus;
                sys->keepTrack[neighbour_plus] = 1;
            }
        }
        neighbour_minus = sys->nnm[dirgeo(lex, ii, sys->volume)];

        if (sys->keepTrack[neighbour_minus] == 0 &&
            sys->lattice[lex] * sys->lattice[neighbour_minus] == 1)
        {
            if (myrand() < prob_add)
            {
                addcount++;
                sys->cluster[sys->cluster_index + addcount] = neighbour_minus;
                sys->keepTrack[neighbour_minus] = 1;
            }
        }
    }
    return addcount;
}
/*
   Construct a Wolff cluster.

   Expands the cluster starting from the seed site, adding neighboring spins
   that are aligned and accepted according to the probability 'prob_add'.

   Returns:
     Total number of sites in the constructed cluster
*/
int create_cluster(ising *restrict sys, double prob_add)
{
    init_cluster(sys);
    int n_old, n_new;

    n_old = 0;
    n_new = 1;

    while (n_old < n_new)
    {
        int addcount = 0;

        for (int idx = n_old; idx < n_new; idx++)
        {
            sys->cluster_index = n_new - 1 + addcount;
            addcount += check_near(sys, idx, prob_add);
        }
        n_old = n_new;
        n_new = addcount + n_old;
    }

    return n_new;
}
/*
   Perform a single Wolff cluster update.

   Constructs a cluster and flips all spins within it.
*/
void cluster_update(ising *restrict sys, double prob_add)
{
    int stop = create_cluster(sys, prob_add); // cluster lenght

    for (int j = 0; j < stop; j++)
    {
        sys->lattice[sys->cluster[j]] *= -1;
    }
}

int main(int argc, char **argv)
{
    if (argc != 7)
    {
        printf("How to use this program:\n\n");
        printf("%s sample beta size dim therm outfile\n\n", argv[0]);
        printf("Parameters:\n");
        printf("  sample   : number of drawn to be extracted\n");
        printf("  beta     : inverse temperature (1 / kT)\n");
        printf("  size     : linear lattice size\n");
        printf("  dim      : spatial dimension of the lattice\n");
        printf("  therm    : number of Wolff updates used for thermalization phase\n");
        printf("  outfile  : name of file where results are stored\n\n");
        printf("Notes:\n");
        printf("  The program first performs 'therm' Wolff cluster updates to reach thermal equilibrium.\n");
        printf("  It then carries out the simulation, performing cluster updates and recording the\n");
        printf("  energy and magnetization at regular intervals.\n\n");
        printf("  Output format: <Energy_per_site> <Magnetization_per_site>\n\n");

        return EXIT_SUCCESS;
    }

    // init variables
    char filename[STRING_LENGTH];

    double beta, prob_add;
    long int sample, therm;
    int size, dim, counter;

    const int measevery = 2;

    sample = (atol(argv[1]) * measevery);
    beta = atof(argv[2]);
    size = atoi(argv[3]);
    dim = atoi(argv[4]);
    therm = atol(argv[5]);

    prob_add = 1 - exp(-2 * beta);

    if (strlen(argv[6]) >= STRING_LENGTH)
    {
        fprintf(stderr, "File name too long (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    else
    {
        strcpy(filename, argv[6]);
    }

    if (size <= 0)
    {
        fprintf(stderr, "'size' must be positive\n");
        return EXIT_FAILURE;
    }

    if (sample <= 0)
    {
        fprintf(stderr, "'sample' must be positive\n");
        return EXIT_FAILURE;
    }

    if (therm <= 0)
    {
        fprintf(stderr, "'therm' must be positive\n");
        return EXIT_FAILURE;
    }

    if (beta <= 0)
    {
        fprintf(stderr, "'beta' must be positive\n");
        return EXIT_FAILURE;
    }

    if (dim <= 0)
    {
        fprintf(stderr, "'dim' must be at least 1\n");
        return EXIT_FAILURE;
    }

    // init data structures

    FILE *fp = fopen(filename, "w");

    if (fp == NULL)
    {
        fprintf(stderr, "Error opening file in %s,%d", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

    const unsigned long int seed1 = (unsigned long int)time(NULL);
    const unsigned long int seed2 = seed1 + 350;

    myrand_init(seed1, seed2);

    ising sys = {0};
    init_system(&sys, size, dim);

    // thermalization to reach equilibrium

    for (int i = 0; i < therm; i++)
    {
        cluster_update(&sys, prob_add);
    }

    // get measures

    counter = measevery;
    for (int j = 0; j < sample; j++)
    {
        counter--;
        cluster_update(&sys, prob_add);
        if (counter == 0)
        {
            fprintf(fp, "%lf %lf \n", energy(&sys), magn(&sys));
            counter = measevery;
        }
    }

    fclose(fp);
    free(sys.lattice);
    free(sys.keepTrack);
    free(sys.cluster);
    free(sys.nnp);
    free(sys.nnm);
}