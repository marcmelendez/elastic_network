/****************************** elastic_network.c *****************************

 Given a list of coordinates, output an elastic newtork of bonds in the format
 "i j K r0" in rows, where i and j are the indices of the bonded points, K is
 the spring strength and r0 the equilibrium distance. Executable version of the
 function in elastic_network.h */

/* Standard library */
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "elastic_network.h"

int main(int argc, char * argv[])
{
  if(argc < 8) { /* Use message */
    printf("Usage: %s <N> <Rc> <K> <Lx> <Ly> <Lz> <file> [DIM] \nParameters:\n", argv[0]);
    printf("\tN:\tNumber of particles to read from file (set to -1 to read from the top of the input file).\n"
           "\tRc:\tCut-off radius for bonds.\n"
           "\tK:\tBond strength parameter.\n"
           "\tLx, Ly, Lz:\tBox dimensions "
           "(enter -1 for no periodic boundary conditions).\n"
           "\tfile:\tFilename of particle positions "
           "(Format: x y z ... by rows).\n"
           "\tDIM:\tDimensionality of space (1, 2 or 3).\n");
    return 0;
  }

  int nbonds = 0; /* Number of bonds generated */
  int DIM; /* Dimensionality */
  int k; /* Coordinate index */
  if(argc == 9) DIM = atoi(argv[8]); else DIM = 3;
  int N = atoi(argv[1]); /* Number of particles */
  FILE * data = fopen(argv[7], "r"); /* Pointer to data file */

  /* Check file pointer */
  if(data == NULL) {
    fprintf(stderr, "Error: unable to open file %s.\n", argv[7]);
    return -1;
  }

  float L[3]; /* Box dimensions */
  float Rc = atof(argv[2]); /* Cut-off radius */
  float K = atof(argv[3]); /* Bond strength parameter */

  /* Set dimensions of box */
  for(k = 0; k < DIM; k++) {
    if(atof(argv[4 + k]) > 0) {
      L[k] = atof(argv[4 + k]);
    }
  }

  /* Get box dimensions */
  for(k = 0; k < DIM; k++) L[k] = atof(argv[4 + k]);

  /* Direct to standard output */
  FILE * output = fopen("/dev/stdout", "w");

  /* Calculate bonds */
  nbonds =
    elastic_network(N, Rc, K, L, data, DIM, 0, output);

  /* Output number of bonds */
  fprintf(stderr, "Generated %d bonds.\n", nbonds);

  /* Clean up and exit */
  fclose(data);
  fclose(output);

  return 0;
}
