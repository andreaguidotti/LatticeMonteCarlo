/*

  Program to compute observables from Monte Carlo data using
  blockking Jackknife analysis with optional reweighting flag
  for high temperature simulation.

  Note:
  Array of pointers, dataSets[], is defined as a simple wrapper
  mechanism to allow iteration over multiple datasets and avoid 
  code repetition.

*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define STRING_LENGTH 50
#define QMAX 5
#define CHI 1

// pointer to a function
typedef double (*f)(double value, double simbeta);

// Struct for storing data.
typedef struct data
{
    double *arr;       // Array of data
    double totSum;     // Total sum of array elements
    double jackSample; // Auxiliary variable for jackknife calculations

    double kahanCorrection; // Correction to be used in sumKahan

} data;

// Struct for storing observable results and jackknife statistics.
typedef struct obs
{
    double avg;         // Simple average of observable
    double std;         // Standard deviation
    double jackavg;     // Jackknife mean
    double jackavgSqrd; // Jackknife mean squared

    double kahanCorrection[2]; //  Corrections to be used in sumKahan

} obs;

/* Relevant functions to compute datasets */
double square(double x, double simbeta) { return x * x; }

/* Compute the reweighting factor used in high temperature simulation */
double reweighting(double windingNumber, double simbeta)
{
    double effQ;
    double reweight;

    // Clamp winding number to [-QMAX, QMAX]
    if (windingNumber > QMAX)
        effQ = QMAX;
    else if (windingNumber < -QMAX)
        effQ = -QMAX;
    else
        effQ = windingNumber;

    reweight = exp(-(effQ * effQ) / (2. * simbeta * CHI));

    return reweight;
}

/* Compute squared winding number times reweighting factor */
double reweightedSquare(double windingNumber, double simbeta)
{
    double res;
    res = pow(windingNumber, 2) * reweighting(windingNumber, simbeta);

    return res;
}

/* Kahan summation algorithm.

   Performs numerically stable summation by tracking a small
   correction term to compensate for floating-point errors.
*/
void sumKahan(double addend, double *sumtot, double *correction)
{
    double correctedAddend, tempSum;

    correctedAddend = addend + *correction;
    tempSum = *sumtot + correctedAddend;
    *correction = (*sumtot - tempSum) + correctedAddend;
    *sumtot = tempSum;
}

/*
  Extract data from the input file.

  Reads 'sampleEff' measurements from the file 'fp', skipping the first
  'therm' lines for thermalization.

  For each measurement, a set of functions 'fData' is applied to generate
  the relevant datasets for computing observables.

  Results are stored in the corresponding dataset arrays.
 */
void extractData(FILE *fp, data *datasets[], f fData[], double simbeta,
                 long int sampleEff, long int therm, int nDatasets)
{
    double value;

    // Skipping 'therm' rows for thermalization
    for (long int row = 0; row < therm; row++)
    {
        // integrity check of input file (must have 1 column)
        if (fscanf(fp, "%lf", &value) != 1)
        {
            fprintf(stderr, "Error: unexpected end of file"
                            "or read error at row %ld\n",
                    row);

            fclose(fp);
            exit(EXIT_FAILURE);
        }
    }
    // Load data and integrity check of input file
    for (long int row = 0; row < sampleEff; row++)
    {
        if (fscanf(fp, "%lf", &value) != 1)
        {
            fprintf(stderr, "Error: unexpected end of file"
                            "or read error at row %ld\n",
                    row);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        // Using datasets and fData arrays to avoid repetition of code
        for (int i = 0; i < nDatasets; i++)
        {
            datasets[i]->arr[row] = fData[i](value, simbeta);
            sumKahan(datasets[i]->arr[row], &(datasets[i]->totSum),
                     &(datasets[i]->kahanCorrection));
        }
    }
    // Add residual Kahan correction after the loop
    for (int i = 0; i < nDatasets; i++)
    {
        datasets[i]->totSum += datasets[i]->kahanCorrection;
    }
}

/* Jackknife leave-one-out.

   Computes the mean of the data excluding the i-th block.
   'blockdim' is the block size, 'nblocks' is the total number of blocks.
*/
static inline void jackLeaveOneOut(data *data, long int i,
                                   long int blockdim, long int nblocks)
{
    long int idx;

    data->jackSample = data->totSum;
    for (long int ii = 0; ii < blockdim; ii++)
    {
        idx = i * blockdim + ii;
        data->jackSample -= data->arr[idx];
    }
    data->jackSample /= (double)((nblocks - 1) * blockdim);
}

/* Accumulate jackknife mean and mean squared using Kahan summation */

static inline void accumulate(obs *o, double value)
{
    sumKahan(value, &o->jackavg, &o->kahanCorrection[0]);
    sumKahan(value * value, &o->jackavgSqrd, &o->kahanCorrection[1]);
}

/* Perform jackknife analysis.

   Computes the jackknife mean and standard deviation for each
   observables  over 'nblocks' blocks of size 'blockdim'.
*/
void jackknife(obs *obsv, data *datasets[], long int nblocks,
               long int blockdim, int useReweighting, int nDatasets)
{
    double value;

    // Get jacksamples for each dataset
    for (long int i = 0; i < nblocks; i++)
    {
        for (int ii = 0; ii < nDatasets; ii++)
        {
            jackLeaveOneOut(datasets[ii], i, blockdim, nblocks);
        }

        // Compute the observable value for this Jackknife sample
        if (useReweighting)
            value = datasets[0]->jackSample / datasets[1]->jackSample;

        else
            value = datasets[0]->jackSample;

        accumulate(obsv, value);
    }

    // Add the Kahan residual from the last addition
    obsv->jackavg += obsv->kahanCorrection[0];
    obsv->jackavgSqrd += obsv->kahanCorrection[1];

    // Normalization
    obsv->jackavg /= nblocks;
    obsv->jackavgSqrd /= nblocks;

    // Standard deviation
    obsv->std = sqrt((nblocks - 1) *
                     (obsv->jackavgSqrd - pow(obsv->jackavg, 2)));
}

/* Compute simple average of observable */
void meanObservable(obs *obsv, data *dataSets[], long int sampleEff, int useReweighting)
{
    if (useReweighting)
        obsv->avg = (dataSets[0]->totSum) / (dataSets[1]->totSum);

    else
        obsv->avg = dataSets[0]->totSum / sampleEff;
}
int main(int argc, char **argv)
{
    if (argc != 8)
    {
        printf("How to use this program: \n");
        printf("%s inputfile, outputfile, sample, blocksize, simbeta, therm, reweight\n\n", argv[0]);
        printf("  inputfile  : path to input data file\n");
        printf("  outputfile : path to output results file\n");
        printf("  sample     : number of measurements\n");
        printf("  blocksize  : size of each jackknife block (>=2)\n");
        printf("  simbeta    : inverse temperature parameter\n");
        printf("  therm      : thermalization steps (<sample)\n");
        printf("  reweight   : flag to use reweighting (0 = no, 1 = yes)\n");
        printf("\n");

        return EXIT_SUCCESS;
    }

    // Initialize variables 
    long int sample, sampleEff, blockdim, nblocks, therm;
    char infile[STRING_LENGTH];
    char outfile[STRING_LENGTH];
    double simbeta;
    int useReweighting, nDatasets;

    if (strlen(argv[1]) >= STRING_LENGTH || strlen(argv[2]) >= STRING_LENGTH)
    {
        fprintf(stderr, "File name too long (%s, %d)\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    else
    {
        strcpy(infile, argv[1]);
        strcpy(outfile, argv[2]);
    }

    sample = atol(argv[3]);
    blockdim = atol(argv[4]);

    simbeta = atof(argv[5]);
    therm = atol(argv[6]);
    useReweighting = atoi(argv[7]);

    // Check command line arguments
    if (sample < 0)
    {
        fprintf(stderr, "'sample' has to be positive\n");
        return EXIT_FAILURE;
    }
    if (blockdim < 2)
    {
        fprintf(stderr, "'blockdim' has to be at least 2\n");
        return EXIT_FAILURE;
    }
    if (simbeta < 0)
    {
        fprintf(stderr, "'beta' has to be positive\n");
        return EXIT_FAILURE;
    }
    if (therm > sample || therm < 0)
    {
        fprintf(stderr, "thermalization has to be less than sample and positive\n");
        return EXIT_FAILURE;
    }

    nblocks = (sample - therm) / blockdim;
    sampleEff = nblocks * blockdim;

    // initialize relevant structures
    data Q2 = {0}, reweight = {0};
    obs obsQ2 = {0};

    data **dataSets;
    f *fData;

    nDatasets = useReweighting ? 2 : 1;

    dataSets = malloc((unsigned)nDatasets * sizeof(data *));
    fData = malloc((unsigned)nDatasets * sizeof(f));

    if (dataSets == NULL || fData == NULL)
    {
        fprintf(stderr, "Error defining arrays at %s,%d\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    if (useReweighting)
    {
        dataSets[0] = &Q2;
        dataSets[1] = &reweight;

        fData[0] = reweightedSquare;
        fData[1] = reweighting;
    }
    else
    {
        dataSets[0] = &Q2;
        fData[0] = square;
    }

    // Open  input and output file
    FILE *fp = fopen(infile, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error opening file in %s,%d\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

    FILE *out = fopen(outfile, "a");
    if (out == NULL)
    {
        fprintf(stderr, "Error defining array %s,%d\n", __FILE__, __LINE__);
        fclose(fp);
        return EXIT_FAILURE;
    }

    // Allocate dynamic arrays
    for (int i = 0; i < nDatasets; i++)
    {
        dataSets[i]->arr = (double *)malloc((unsigned)sampleEff * sizeof(double));
        if (dataSets[i]->arr == NULL)
        {
            fprintf(stderr, "Error defining array %s,%d\n", __FILE__, __LINE__);
            fclose(fp);

            return EXIT_FAILURE;
        }
    }
    // Extract data from file
    extractData(fp, dataSets, fData, simbeta, sampleEff, therm, nDatasets);

    // Compute average value for observables
    meanObservable(&obsQ2, dataSets, sampleEff, useReweighting);

    // Perform jackknife analysis
    jackknife(&obsQ2, dataSets, nblocks, blockdim, useReweighting, nDatasets);

    // Export results
    fprintf(out, "%.100lf %.100lf\n", obsQ2.avg, obsQ2.std);

    // Save results in output file
    fprintf(stderr, "All valid data saved in %s\n", outfile);

    // Close files and free allocated memory
    fclose(fp);
    fclose(out);
    for (int i = 0; i < nDatasets; i++)
    {
        free(dataSets[i]->arr);
    }
    free(dataSets);
    free(fData);

    return EXIT_SUCCESS;
}