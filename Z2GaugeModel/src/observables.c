
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define STRING_LENGTH 50
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

// Struct for storing (Wilson loop) data.
typedef struct data
{
    double *arr; // Array of data
    double totSum; // Total sum of array elements
    double jackSample; // Auxiliary variable for jackknife calculations

    double kahanCorrection; // Correction to be used in sumKahan

    int Wt; // Temporal size Wilson loop
    int Ws; // Spazial size Wilson loop

} data;

// Struct for storing observable results and jackknife statistics.
typedef struct obs
{
    double avg;         // Simple average of observable
    double std;         // Standard deviation
    double jackavg;     // Jackknife mean
    double jackavgSqrd; // Jackknife mean squared

    double kahanCorrection[2]; //  corrections to be used in sumKahan

} obs;

/* Validate input file format and data range.

   Checks that each of the 'sampleEff' lines from 'fp' (starting at 'pos')
   has at least 'ncolumns' valid floating point values and that each value is ≤ 1.

   If any line has fewer columns, or a value is out of range, or if an unexpected
   end of file is reached, an error message is printed and the program terminates.
*/
void checkFile(FILE *fp, long pos, long int sampleEff, int ncolumns)
{
    char buffer[4000];
    int fileHasErrors = 0;

    double value;
    int read;

    for (long int row = 0; row < sampleEff; row++)
    {
        // read one line from the file into 'buffer'
        if (fgets(buffer, sizeof(buffer), fp) == NULL)
        {
            fprintf(stderr, "Error: unexpected end of file or read error at row %ld\n", row);
            fclose(fp);
            exit(EXIT_FAILURE);
        }

        char *ptr = buffer;

        // check that each line has at least ncolumns and expected values
        for (int col = 0; col < ncolumns; col++)
        {
            if (sscanf(ptr, "%lf%n", &value, &read) != 1)
            {
                fprintf(stderr, "Line %ld has fewer columns "
                                "than expected (expected %d)\n",
                        row + 1, ncolumns);
                fileHasErrors = 1;
                break; // stop checking this line
            }
            if (value > 1)
            {
                fprintf(stderr, "invalid value %lf at line %ld column %d\n",
                        value, row + 1, col + 1);
                fileHasErrors = 1;
            }
            ptr += read;
        }
    }
    if (fileHasErrors)
    {
        fprintf(stderr, "\nError: invalid input file. Please fix and retry.\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }
    // reset pointer so data analysis can start from the top
    if (fseek(fp, pos, SEEK_SET) != 0)
    {
        fprintf(stderr, "Error: unable to reset file pointer after checkFile\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }
}

// Reset all fields in a data structure
void resetData(data *W)
{
    W->totSum = 0.0;
    W->jackSample = 0.0;
    W->kahanCorrection = 0.0;
    W->Wt = 0;
    W->Ws = 0;
}
// Reset all fields in a obs structure
void resetObs(obs *Potential)
{
    Potential->avg = 0.0;
    Potential->std = 0.0;
    Potential->jackavg = 0.0;
    Potential->jackavgSqrd = 0.0;
    Potential->kahanCorrection[0] = 0.0;
    Potential->kahanCorrection[1] = 0.0;
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

/* Extract data from a file for a specific column.

   Reads 'sampleEff' lines from file 'fp', extracts the value in
   column 'target_col' and stores it in the array 'W->arr'.

   Uses Kahan summation to accumulate the total sum.
   Updates the Wilson loop dimensions Wt and Ws.
*/
void extractData(FILE *fp, data *W, long int sampleEff, long pos,
                 int target_col, int Wt_counter, int Ws_counter)
{
    // buffer to store one line of text from the file
    char buffer[4000];

    double value;
    int read;

    for (long int row = 0; row < sampleEff; row++)
    {
        // read one line from the file into 'buffer'
        if (fgets(buffer, sizeof(buffer), fp) == NULL)
        {
            fprintf(stderr, "Error: unexpected end of file"
                            "or read error at row %ld\n",
                    row);
            free(W->arr);
            fclose(fp);
            exit(EXIT_FAILURE);
        }
        // pointer to current position in the line
        char *ptr = buffer;

        // reach desired column
        for (int col = 0; col <= target_col; col++)
        {
            // read a double and count chars read
            sscanf(ptr, "%lf%n", &value, &read);

            // advance pointer to next value in the line
            ptr += read;
        }
        // got to target column, now save data
        W->arr[row] = value;

        // Accumulate sum with Kahan correction for numerical accuracy
        sumKahan(value, &(W->totSum), &(W->kahanCorrection));

        W->Wt = Wt_counter;
        W->Ws = Ws_counter;
    }
    // Add residual Kahan correction after the loop
    W->totSum += W->kahanCorrection;

    // restore file pointer to position after thermalization ('pos)
    if (fseek(fp, pos, SEEK_SET) != 0)
    {
        fprintf(stderr,
                "Error: unable to reset file pointer to saved position, "
                "extraction of data in column %d failed\n",
                target_col);

        fclose(fp);
        free(W->arr);
        exit(EXIT_FAILURE);
    }
}

/* Compute String Potential.

   Given the average Wilson loop value and temporal size Wt,
   returns the string potential.
*/
double StringU(double Wilson, int Wt)
{
    double res;
    res = -log(Wilson);
    res /= (double)Wt;

    return res;
}

/* Jackknife leave-one-out.

   Computes the mean of the data excluding the i-th block.
   'blockdim' is the block size, 'nblocks' is the total number of blocks.
*/
static inline void jackLeaveOneOut(data *data, long int i, long int blockdim, long int nblocks)
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

// Accumulate mean and mean squared using Kahan summation
static inline void accumulate(obs *o, double value)
{
    sumKahan(value, &o->jackavg, &o->kahanCorrection[0]);
    sumKahan(value * value, &o->jackavgSqrd, &o->kahanCorrection[1]);
}

/* Perform jackknife analysis.

   Computes the jackknife mean and standard deviation for the Wilson 
   loop and String Potential over 'nblocks' blocks of size 'blockdim'.
*/
void jackknife(obs *Potential, obs *Wilson, data *W,
               long int nblocks, long int blockdim, bool skipStringU)
{
    double valuePotential;
    long int nblocks_valid = 0;

    for (long int i = 0; i < nblocks; i++)
    {
        jackLeaveOneOut(W, i, blockdim, nblocks);
        accumulate(Wilson, W->jackSample);

        // skip Static potential computation
        if (skipStringU)
            continue;

        if (W->jackSample <= 0)
        {
            fprintf(stderr, "Skipping block %ld for Wilson loop "
                            "(Wt=%d, Ws=%d) due to non-positive "
                            "mean value\n",
                    i, W->Wt, W->Ws);
            continue;
        }
        nblocks_valid += 1;

        // mean and mean squared of jacksamples

        valuePotential = StringU(W->jackSample, W->Wt);

        accumulate(Potential, valuePotential);
    }
    // Add the Kahan residual from the last addition

    Wilson->jackavg += Wilson->kahanCorrection[0];
    Wilson->jackavgSqrd += Wilson->kahanCorrection[1];

    // Normalization

    Wilson->jackavg /= nblocks;
    Wilson->jackavgSqrd /= nblocks;

    // Standard deviation

    Wilson->std = sqrt((nblocks - 1) * (Wilson->jackavgSqrd - pow(Wilson->jackavg, 2)));

    //skip Static potential computation
    if (skipStringU)
        return;

    // avoid unreliable jackknife analysis with too few valid blocks
    if (nblocks_valid < nblocks / 3)
    {
        fprintf(stderr,
                "Warning: too few valid blocks (%ld/%ld) for Wilson loop "
                "(Wt=%d, Ws=%d) — skipping this observable\n",
                nblocks_valid, nblocks, W->Wt, W->Ws);

        Potential->avg = NAN;
        Potential->std = NAN;
        return;
    }
    else
    {
        // Add the Kahan residual from the last addition

        Potential->jackavg += Potential->kahanCorrection[0];
        Potential->jackavgSqrd += Potential->kahanCorrection[1];

        // Normalization

        Potential->jackavg /= nblocks_valid;
        Potential->jackavgSqrd /= nblocks_valid;

        // Standard deviation

        Potential->std = sqrt((nblocks_valid - 1) * (Potential->jackavgSqrd - pow(Potential->jackavg, 2)));
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
        printf("  therm      : thermalization steps (>=0)\n");
        printf("\n");

        return EXIT_SUCCESS;
    }

    // initialize variables and data structures

    obs Potential = {0};
    obs Wilson = {0};
    data W = {0};

    long int sample, sampleEff, blockdim, nblocks, therm;
    int size, ncolumns;
    double beta;

    char infile[STRING_LENGTH];
    char outfile[STRING_LENGTH];

    sample = atol(argv[3]);
    blockdim = atol(argv[4]);
    beta = atof(argv[5]);
    size = atoi(argv[6]);
    therm = atol(argv[7]);

    // Check input

    if (blockdim < 2)
    {
        fprintf(stderr, "blockdim has to be at least 2\n");
        return EXIT_FAILURE;
    }
    if (size < 8)
    {
        fprintf(stderr, "size has to be at least 8\n");
        return EXIT_FAILURE;
    }
    if (therm < 0)
    {
        fprintf(stderr, "thermalization has to be at least 0\n");
        return EXIT_FAILURE;
    }
    nblocks = (sample - therm) / blockdim;
    sampleEff = nblocks * blockdim;
    ncolumns = MIN(size / 4, 8) * (size / 4);

    // Open  input and output file

    strcpy(infile, argv[1]);
    strcpy(outfile, argv[2]);

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

    // Discard thermalization

    char buffer[4000];
    for (long int i = 0; i < therm; i++)
    {
        if (fgets(buffer, sizeof(buffer), fp) == NULL)
            break;
    }

    // Save pointer position in FILE after thermalization
    long pos = ftell(fp);
    // Check file integrity
    checkFile(fp, pos, sampleEff, ncolumns);

    // Allocate dynamic arrays

    W.arr = (double *)malloc((unsigned)sampleEff * sizeof(double));
    if (W.arr == NULL)
    {
        fprintf(stderr, "Error defining array %s,%d\n", __FILE__, __LINE__);
        fclose(fp);

        return EXIT_FAILURE;
    }

    // Determine lattice loop parameters

    int Wt_max = MIN(size / 4, 8), Wt_counter = 1, Ws_counter = 1;

    // Loop over each Wilson loop column

    for (int target_col = 0; target_col < ncolumns; target_col++)
    {   
        bool skipStringU = false;

        resetData(&W);
        resetObs(&Wilson);
        resetObs(&Potential);

        extractData(fp, &W, sampleEff, pos, target_col, Wt_counter, Ws_counter);

        // Compute average potential for this column
        Wilson.avg = W.totSum / (double)sampleEff;
        if (Wilson.avg <= 0)
        {
            Potential.avg = NAN;
            Potential.std = NAN;
            skipStringU = true;

            fprintf(stderr,
                    "Warning: skipping column %d (Wt-%d Ws-%d) due "
                    "to non-positive average Wilson loop value %.12lf\n",
                    target_col, W.Wt, W.Ws, Wilson.avg);
        }
        else
            Potential.avg = StringU(Wilson.avg, W.Wt);

        // Perform jackknife analysis

        jackknife(&Potential, &Wilson, &W, nblocks, blockdim, skipStringU);

        // Export results

        fprintf(out, "%lf %d %d %.18lf %.18lf %.18lf %.18lf\n", 
                beta, W.Ws, W.Wt, Wilson.avg, Wilson.std, Potential.avg, Potential.std);

        // Update lattice counters

        Wt_counter += 1;
        if (Wt_counter > Wt_max)
        {
            Ws_counter += 1;
            Wt_counter = 1;
        }
    }

    fprintf(stderr, "\nAll valid data saved in %s\n", outfile);

    // Close files and free allocated memory

    fclose(fp);
    fclose(out);
    free(W.arr);

    return EXIT_SUCCESS;
}