#ifdef PRINTSTRESS
/***************************************************************************
 * This program computes the stress at a point x in the domain
 * Sylvie Aubry, Apr 25 2008
 *
 **************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Matrix.h"

void PrintStress(Home_t *home)
{

  FILE *fpInf;
  char format[15];
  int i,j,k;
  int nx, ny;
  real8 r[3], s1[3][3],s1loc[3][3];
  double** GridPts[3];
  
  // Define the points where to print the stress
  real8 difX, difY, x, y;
  nx = 201;
  ny = 201;

  for (i=0;i<3;i++) 
    GridPts[i] = (double **)malloc(sizeof(double*)*nx);

  for (i=0;i<3;i++) 
    for (k=0;k<nx;k++) 
      {
	GridPts[i][k] = (double *)malloc(sizeof(double)*ny);
      }

  difX = home->param->Lx/(1.0*nx);
  difY = home->param->Ly/(1.0*ny);

  for (i=0; i<nx; i++) 
    {
      x = home->param->minSideX + difX*i;
      for (j=0; j<ny; j++) 
	{
	  y = home->param->minSideY + difY*j;
	  GridPts[0][i][j] = x;
	  GridPts[1][i][j] = y;
	  GridPts[2][i][j] = 100.0;
	}
    }

  
  if (home->myDomain == 0) 
    printf(" \n\n WRITING STRESS TO FILE %d x %d points \n\n",nx, ny);
  
  sprintf(format, "Stress%d.out",home->myDomain);
  printf("%s \n",format);
  fpInf = fopen(format,"w");

  for (i=0; i<nx; i++)
    {
      if (home->myDomain == 0)  printf("doing i=%d\n",i);
      for (j=0; j<ny; j++) 
	{
	  r[0] = GridPts[0][i][j];
	  r[1] = GridPts[1][i][j];
	  r[2] = GridPts[2][i][j]; 

	  /* infinite medium stress */
	  Matrix33_Zero(s1loc);
	  AllSegmentStress(home,r[0],r[1],r[2],s1loc);
 
#ifdef PARALLEL
	  MPI_Allreduce(s1loc, s1, 9, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
#else
        for (int ii = 0; ii < 3; ii++)
            for (int jj = 0; jj < 3; jj++)
                s1[ii][jj] = s1loc[ii][jj];
#endif

        /* Print stresses */ 
        fprintf(fpInf,"%.15e %.15e %.15e      %.15e %.15e %.15e %.15e %.15e %.15e\n",
               r[0],r[1],r[2], s1[0][0],s1[1][1],s1[2][2],s1[1][2],s1[2][0],s1[0][1]);
	}
    }

  fclose(fpInf);

 for (i=0;i<3;i++) 
    for (k=0;k<nx;k++) 
      {
	free(GridPts[i][k]);
      }

  for (i=0;i<3;i++) 
    {
      free(GridPts[i]);
    }  
}

#endif
