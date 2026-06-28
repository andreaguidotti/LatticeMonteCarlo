#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include"../include/geometry.h"

/* Convert cartesian coordinates to lexicographic index.

   Given a lattice of size 'latticeSize' in each of 'dim' dimensions,
   this function converts an array of cartesian coordinates
   cartCoords[0..dim-1] into its corresponding lexicographic index.

   Example:

     For dim = 3 and latticeSize = L:
       cartCoords = (x0, x1, x2)
       => lex = x0 + x1*L + x2*L^2

   This mapping allows a multidimensional lattice to be represented
   as a one-dimensional array using lexicographic ordering.
*/

void cart_to_lex(long int *lex, int *cartcoord, int dim, int latticeSize)
{
    long int multiplier, lexIndex;

    multiplier = 1;
    lexIndex = 0;

    for (int i = 0; i < dim; i++)
    {
        lexIndex += cartcoord[i] * multiplier;
        multiplier *= latticeSize;
    }

    *lex = lexIndex;
}

/* Convert a lexicographic index to cartesian coordinates.

   Given a lattice of size 'latticeSize' in each of 'dim' dimensions,
   this function converts an integer index 'lex' (0 <= lex < latticeSize^dim)
   into its corresponding cartesian coordinates.

   Example:

     For dim = 3, latticeSize = L and lex = x0 + x1*L + x2*L^2: 
     cartCoords[2] = lex/L^2, cartCoords[1] = (lex - x2*L^2)/L^1, 
     cartCoords[1] = (lex - x2*L^2 -x1*L)/L^0, 
     => cartCoords = (x0, x1, x2)
*/

void lex_to_cart(long int lex,  int *cartcoord, int dim, int latticeSize)
{
    long int multiplier;
    multiplier = 1;

    for (int i = 0; i < dim-1; i++)
    {
        multiplier *= latticeSize;
    }

    for (int d = dim - 1; d >= 0; d--)
    {
        cartcoord[d] = (int)(lex / multiplier);
        lex -= cartcoord[d] * multiplier;
        multiplier /= latticeSize;
    }
}
/* Initialize the nearest neighbor lookup tables for a periodic lattice.

    For each site on a d dimensional periodic lattice, this function determines
    the lexicographic indices of its nearest neighbors in the + and − directions
    along each spatial dimension.

   INPUT:
     latticeSize  – number of sites along each lattice direction (L)
     dim          – number of spatial dimensions

   OUTPUT:
     nnp – array of lexicographic indices of forward neighbors (+ direction)
     nnm – array of lexicographic indices of backward neighbors (− direction)
            Both arrays have size (volume * dim), where volume = L^dim.

   HOW IT WORKS:
     Each site is represented by a unique lexicographic index "site".
     The function converts this index to its cartesian coordinates.
     Then, for each dimension:
       1. It increments or decrements the coordinate (applying periodic boundary conditions),
       2. Converts the modified coordinates back to lexicographic form,
       3. Stores the resulting indices in the neighbor arrays.

   Example:
     In 2D (L = 4), site (x,y) = (0,3) has:
       forward neighbors  → (1,3), (0,0)
       backward neighbors → (3,3), (0,2)
*/
void init_geo(long int *restrict nnp,
           long int *restrict nnm,
           int latticeSize,
           int dim)
{
    long int volume;
    int coordPlus, coordMinus, coordValue;
    long int lexFoward, lexBackward;

    int *cartcoord = ( int *)malloc((unsigned int)dim * sizeof(int));
    if (cartcoord == NULL)
    {
        fprintf(stderr, "Allocation problem at (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    volume = 1;
    for (int d = 0; d < dim; d++)
    {
        volume *= latticeSize;
    }
  //for each element i, in lexCoord, go to cartCoordinates and maps to lex its neighb
    for (long int i = 0; i < volume; i++)
    {
        lex_to_cart(i, cartcoord, dim, latticeSize);
        for (int ii = 0; ii < dim; ii++)
        {
            coordValue = cartcoord[ii];

            coordPlus = coordValue + 1;
            if (coordPlus >= latticeSize)
            {
                coordPlus -= latticeSize;
            }
            cartcoord[ii] = coordPlus;
            cart_to_lex(&lexFoward, cartcoord, dim, latticeSize);
            nnp[dirgeo(i, ii, volume)] = lexFoward;

            coordMinus = coordValue - 1;
            if (coordMinus < 0)
            {
                coordMinus += latticeSize;
            }
            cartcoord[ii] = coordMinus;
            cart_to_lex(&lexBackward, cartcoord, dim, latticeSize);
            nnm[dirgeo(i, ii, volume)] = lexBackward;

            cartcoord[ii] = coordValue;
        }
    }
    free(cartcoord);
}
