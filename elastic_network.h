/****************************** elastic_network.h *****************************

 Given a list of coordinates, output an elastic newtork of bonds in the format
 "i j K r0" in rows, where i and j are the indices of the bonded points, K is
 the spring strength and r0 the equilibrium distance. */

/* Standard library */
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

/* Print an elastic network given a coordinate file */
int elastic_network(int N,            /* Number of particles */
                    float Rc,         /* Cutoff radius */
                    float K,          /* Bond strength parameter */
                    float * L,        /* Box dimensions */
                    FILE * positions, /* Input file containing coordinates */
                    int dim,          /* Dimensionality */
                    int offset,       /* Offset in coordinate indices */
                    FILE * output)    /* Output file */
{
  /* Check file pointer */
  if(positions == NULL) {
    fprintf(stderr, "Error: unable to open file.\n");
    return -1;
  }

  char line[250]; /* Line buffer */
  int nread; /* Number of fields read */
  if(N < 0) {
    fgets(line, 250, positions);
    nread = sscanf(line, "%d", &N);
    fprintf(stderr, "Number of particles: %d\n", N);
  }
  float pos[N][dim]; /* Array of particle positions */
  int periodicL[3] = {0, 0, 0}; /* Flag indicating the periodic directions */
  int i, j; /* Particle indices */
  int k, l; /* Coordinate indices */
  int m, nm[3]; /* Cell indices */
  double r; /* Distance between particles */
  double r2; /* Distance squared */
  double rij[dim]; /* Displacement */
  int ncells = 1; /* Number of neighbour cells */
  int llist[N]; /* Linked list for particles */
  int n[dim]; /* Number of cells in each direction */
  float cellsize[dim]; /* Size of the cell in each direction */
  int cellcoord[dim], ncellcoord[dim]; /* Cell (integer) coordinates */
  int cellidx, tmpidx; /* Cell indices */
  int nbonds = 0; /* Total number of bonds output */

  /* Set dimensions of box */
  for(k = 0; k < dim; k++) {
    if(L[k] > 0) {
      periodicL[k] = 1;
    }
  }

  i = 0;
  while(fgets(line, 250, positions)) {
    /* Read in the position of particle i */
    switch(dim) {
      case 1:
        nread = sscanf(line, "%f", &pos[i][0]);
        break;
      case 2:
        nread = sscanf(line, "%f %f", &pos[i][0], &pos[i][1]);
        break;
      case 3:
        nread = sscanf(line, "%f %f %f", &pos[i][0], &pos[i][1], &pos[i][2]);
        break;
      default:
        fprintf(stderr, "Error: invalid dimensionality.\n");
        return -1;
    }
    if(nread == dim) {
      /* Set box size if not specified */
      for(k = 0; k < dim; k++) {
        if(periodicL[k] == 0) {
          if(fabs(pos[i][k]) > L[k]) L[k] = 2*fabs(pos[i][k]);
        }
      }
      i++; /* Move on to next particle */
      if(i >= N) break;
    }
  }

  for(k = 0; k < dim; k++) {
    if(L[k] < 3*Rc) {
      L[k] = 3*Rc;
      fprintf(stderr, "Warning: L[%d] box dimension too narrow, increased to "
                      "L[%d] = %f\n", k, k, L[k]);
    }
  }

  /* Have we read the positions of all the particles? */
  if(i < N) {
    fprintf(stderr, "Error: end of file reached prematurely.\n");
    return -1;
  }

  /* Calculate the number of cells and cell size */
  for(k = 0; k < dim; k++) {
    n[k] = (int) floorf(L[k]/Rc);
    cellsize[k] = L[k]/((float) n[k]);
    ncells *= n[k];
  }
  int hlist[ncells]; /* Head list for cells */

  /* Initialise head list */
  for(m = 0; m < ncells; m++) hlist[m] = -1;

  /* Build the linked list */
  for(i = 0; i < N; i++) {
    /* Cell coordinates of position */
    for(k = 0; k < dim; k++) {
      cellcoord[k] = (int) floorf((pos[i][k] + .5*L[k])/cellsize[k]);
      if(cellcoord[k] == n[k]) cellcoord[k] = 0;
    }

    /* Get cell index */
    cellidx = 0; /* Clear cell index */
    for(k = dim; k > 0; k--) {
      tmpidx = cellcoord[k - 1];
      for(l = k - 1; l > 0; l--)
        tmpidx *= n[l - 1];
      cellidx += tmpidx;
    }

    /* Add particle to linked list */
    llist[i] = hlist[cellidx];
    hlist[cellidx] = i;
  }

  /* Loop over all particle indices working out which ones are in bond range Rc */
  for(i = 0; i < N; i++) {
    /* Get position cell number */
    for(k = 0; k < dim; k++) {
      cellcoord[k] = (int) floorf((pos[i][k] + .5*L[k])/cellsize[k]);
    }
    /* Loop over neighbour cells */
    nm[0] = nm[1] = nm[2] = -1;
    for(m = 0; m < pow(3,dim); m++) {
      /* Neighbour cell coordinates */
      for(k = 0; k < dim; k++) ncellcoord[k] = cellcoord[k] + nm[k];

      /* Periodic boundary conditions */
      for(k = 0; k < dim; k++) {
        if(ncellcoord[k] < 0) ncellcoord[k] += n[k];
        if(ncellcoord[k] >= n[k]) ncellcoord[k] -= n[k];
      }

      /* Determine neighbour cell index */
      cellidx = 0; /* Clear cell index */
      for(k = dim; k > 0; k--) {
        tmpidx = ncellcoord[k - 1];
        for(l = k - 1; l > 0; l--)
          tmpidx *= n[l - 1];
        cellidx += tmpidx;
      }

      /* Go through linked lists checking distances */
      j = hlist[cellidx];
      while(j >= 0) {
        if(j >= i) {j = llist[j]; continue;} /* No bonds from a particle to itself */

        /* Displacement */
        for(k = 0; k < 3; k++) rij[k] = pos[j][k] - pos[i][k];

        /* Periodic boundary conditions */
        for(k = 0; k < dim; k++) {
          if(periodicL[k] == 1)
            rij[k] -= floorf(rij[k]/L[k] + .5)*L[k];
        }

        /* Distance */
        r2 = 0; for(k = 0; k < dim; k++) r2 += rij[k]*rij[k];
        r = sqrt(r2);

        /* Output bond */
        if(r <= Rc) {
          fprintf(output, "%d\t%d\t%f\t%f\n", offset + i, offset + j, K, r);
          nbonds++;
        }

        j = llist[j];
      }

      /* Next neighbour cell */
      nm[0]++;
      if(nm[0] > 1) {nm[0] = -1; nm[1]++;}
      if(nm[1] > 1) {nm[1] = -1; nm[2]++;}
    }
  }

  return nbonds;
}
