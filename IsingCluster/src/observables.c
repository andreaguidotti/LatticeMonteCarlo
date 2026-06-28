/*
  Note:

  Arrays of pointers, dataSets[] and obsv[], are defined
  as a simple wrapper mechanism to allow iteration over multiple
  datasets and observables. This avoids repeating the same code
  for each dataset or observable while still keeping explicit
  variables for each observable (e.g., binderU, Chi) to make the
  code clear and immediately understandable (hopefully).

*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define STRING_LENGTH 50
#define nDatasets 3
#define nObs 2

// pointer to a function
typedef double (*f)(double valueM);

// Struct for storing data.
typedef struct data
{
    double totSum; // Total sum of array elements
    double aux;    // Auxiliary variable for jackknife calculations
    double *arr;   // Array of data

    double kahanCorrection; // correction to be used in sumKahan

} data;

// Struct for storing observable results and jackknife statistics.
typedef struct obs
{
    double avg;         // Simple average of observable
    double std;         // Standard deviation
    double jackavg;     // Jackknife mean
    double jackavgSqrd; // Jackknife squared mean

    double kahanCorrection[2]; //  corrections to be used in sumKahan

} obs;

/* relevant functions to compute datasets */

double f_abs(double x) { return fabs(x); }
double square(double x) { return pow(x, 2); }
double quad(double x) { return pow(x, 4); }

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
void extractData(FILE *fp, data *datasets[], f fData[],
                 long int sampleEff, long int therm)
{
    double valueE, valueM;

    // skipping 'therm' rows for thermalization
    for (long int row = 0; row < therm; row++)
    {
        // integrity check of input file (must have 2 columns)
        if (fscanf(fp, "%lf %lf", &valueE, &valueM) != 2)
        {
            fprintf(stderr, "Error: unexpected end of file"
                            "or read error at row %ld\n",
                    row);

            fclose(fp);
            exit(EXIT_FAILURE);
        }
    }
    // load data and integrity check of input file
    for (long int row = 0; row < sampleEff; row++)
    {
        // valueE (Energy per unit volume) currently not used
        if (fscanf(fp, "%lf %lf", &valueE, &valueM) != 2)
        {
            fprintf(stderr, "Error: unexpected end of file"
                            "or read error at row %ld\n",
                    row);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        // using datasets[] and fData arrays to avoid repetition of code
        for (int i = 0; i < nDatasets; i++)
        {
            datasets[i]->arr[row] = fData[i](valueM);
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

/* compute Binder cumulant */

double Ubinder(double m4, double m2)
{
    return m4 / (m2 * m2);
}
/* compute magnetic susceptibility */

double chi(double m2, double absm, double size)
{
    return size * size * (m2 - absm * absm);
}

/* Jackknife leave-one-out.

   Computes the mean of the data excluding the i-th block.
   'blockdim' is the block size, 'nblocks' is the total number of blocks.
*/
static inline void jackLeaveOneOut(data *data, long int i,
                                   long int blockdim, long int nblocks)
{
    long int idx;

    data->aux = data->totSum;
    for (long int ii = 0; ii < blockdim; ii++)
    {
        idx = i * blockdim + ii;
        data->aux -= data->arr[idx];
    }
    data->aux /= (double)((nblocks - 1) * blockdim);
}

// Accumulate mean and mean squared using Kahan summation

static inline void accumulate(obs *o, double value)
{
    sumKahan(value, &o->jackavg, &o->kahanCorrection[0]);
    sumKahan(value * value, &o->jackavgSqrd, &o->kahanCorrection[1]);
}

/* Perform jackknife analysis.

   Computes the jackknife mean and standard deviation for each
   observables  over 'nblocks' blocks of size 'blockdim'.
*/
void jackknife(obs *binderU, obs *Chi, data *datasets[], obs *obsv[],
               long int nblocks, long int blockdim, int size)
{
    double valueChi, valueU;

    // get jacksamples for each dataset
    for (long int i = 0; i < nblocks; i++)
    {
        for (int ii = 0; ii < nDatasets; ii++)
        {
            jackLeaveOneOut(datasets[ii], i, blockdim, nblocks);
        }

        // mean and mean squared of jacksamples
        valueChi = chi(datasets[1]->aux, datasets[0]->aux, size);
        valueU = Ubinder(datasets[2]->aux, datasets[1]->aux);

        accumulate(binderU, valueU);
        accumulate(Chi, valueChi);
    }

    // using obsv[] as a wrapper to avoid code repetition
    for (int i = 0; i < nObs; i++)
    {
        // Add the Kahan residual from the last addition
        obsv[i]->jackavg += obsv[i]->kahanCorrection[0];
        obsv[i]->jackavgSqrd += obsv[i]->kahanCorrection[1];

        // normalization
        obsv[i]->jackavg /= nblocks;
        obsv[i]->jackavgSqrd /= nblocks;

        // Standard deviation
        obsv[i]->std = sqrt((nblocks - 1) *
                            (obsv[i]->jackavgSqrd - pow(obsv[i]->jackavg, 2)));
    }
}
int main(int argc, char **argv)
{
    if (argc != 8)
    {
        printf("How to use this program: \n");
        printf("%s inputfile, outputfile, sample, blocksize, beta, size, therm\n\n", argv[0]);
        printf("  inputfile  : path to input data file\n");
        printf("  outputfile : path to output results file\n");
        printf("  sample     : number of measurements\n");
        printf("  blocksize  : size of each jackknife block (>=2)\n");
        printf("  beta       : inverse temperature parameter\n");
        printf("  size       : lattice size (>=4)\n");
        printf("  therm      : thermalization steps (<sample)\n");
        printf("\n");

        return EXIT_SUCCESS;
    }

    // initialize variables and data structures
    long int sample, sampleEff, blockdim, nblocks, therm;
    char infile[STRING_LENGTH];
    char outfile[STRING_LENGTH];
    double beta;
    int size;

    data absM = {0}, m2 = {0}, m4 = {0};
    obs binderU = {0}, Chi = {0};

    data *dataSets[nDatasets] = {&absM, &m2, &m4};
    obs *obsv[nObs] = {&binderU, &Chi};

    // Function array to generate datasets from raw values using a loop structure
    f fData[nDatasets] = {f_abs, square, quad};

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

    beta = atof(argv[5]);
    size = atoi(argv[6]);
    therm = atol(argv[7]);

    // Check input
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
    if (beta < 0)
    {
        fprintf(stderr, "'beta' has to be positive\n");
        return EXIT_FAILURE;
    }
    if (size < 4)
    {
        fprintf(stderr, "size has to be at least 4\n");
        return EXIT_FAILURE;
    }
    if (therm > sample || therm < 0)
    {
        fprintf(stderr, "thermalization has to be less than sample and positive\n");
        return EXIT_FAILURE;
    }

    nblocks = (sample - therm) / blockdim;
    sampleEff = nblocks * blockdim;

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

    extractData(fp, dataSets, fData, sampleEff, therm);

    // compute average value for observables
    binderU.avg = Ubinder(m4.totSum / sampleEff, m2.totSum / sampleEff);
    Chi.avg = chi(m2.totSum / sampleEff, absM.totSum / sampleEff, size);

    // Perform jackknife analysis
    jackknife(&binderU, &Chi, dataSets, obsv, nblocks, blockdim, size);

    // Export results
    fprintf(out, "%lf ", beta);
    for (int i = 0; i < nObs; i++)
    {
        fprintf(out, "%lf %lf ", obsv[i]->avg, obsv[i]->std);
    }
    fprintf(out, "\n");
    fprintf(stderr, "All valid data saved in %s\n", outfile);

    // Close files and free allocated memory
    fclose(fp);
    fclose(out);
    for (int i = 0; i < nDatasets; i++)
    {
        free(dataSets[i]->arr);
    }
    return EXIT_SUCCESS;
}