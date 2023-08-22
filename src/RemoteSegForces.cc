/****************************************************************************
 *
 *      Module:
 *      Description:   This module contains a number of functions used to
 *                     evaluate stress on dislocation segments from remote
 *                     sources.  For each cell, remote stress is evaluated
 *                     at some number of points, the stress filed defined
 *                     by these points is then used to interpolate the
 *                     the remote stress for dislocation segments within
 *                     the cell.
 *
 *      Included functions:
 *
 *          GaussQuadCoeffHalf()
 *          RemoteForceOneSeg()
 *          SegForceFromTaylorExp()
 *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"

#if 0 // Debug FMMs methods
static void N2Stress(Home_t *home, int cellX, int cellY, int cellZ, real8 R[3], 
		     real8 sigma[3][3], real8 sigmalocal[3][3], int lpbc)
{
	int xSkip1, xSkip2, xSkip3;
	int ySkip1, ySkip2, ySkip3;
	int zSkip1, zSkip2, zSkip3;

	int     ifar = 0;
        int     cellIndex, arm;
        int     i,j, cx, cy, cz;
        int     minXIndex, minYIndex, minZIndex;
        int     maxXIndex, maxYIndex, maxZIndex;
	int     ixpbc, iypbc, izpbc;
	double  pbcp1x, pbcp1y, pbcp1z;
	double  pbcp2x, pbcp2y, pbcp2z;	      
        real8   xc, yc, zc;
        real8   a, MU, NU;
        real8   bx, by, bz;
        real8   p1x, p1y, p1z;
        real8   p2x, p2y, p2z;
	real8   tin[3], t[3];
        real8   stress[6];
        Node_t  *node1, *node2;
        Cell_t  *cell;
        Param_t *param;

        param = home->param;

        a     = param->rc;
        MU    = param->shearModulus;
        NU    = param->pois;

	Matrix33_Zero(sigma);
	Matrix33_Zero(sigmalocal);


	/* Get right cells bounds */
         xSkip1 = cellX - 1 ;
         if (xSkip1 < 0) {
             if (param->xBoundType == Periodic)
                 xSkip1 = param->nXcells;
             else
                 xSkip1 = 0;
         }
         xSkip2 = cellX ;
         xSkip3 = cellX + 1 ;
         if (xSkip3 > param->nXcells) {
             if (param->xBoundType == Periodic)
                 xSkip3 = 0 ;
             else
                 xSkip3 = param->nXcells;
         }

         ySkip1 = cellY - 1 ;
         if (ySkip1 < 0) {
             if (param->yBoundType == Periodic)
                 ySkip1 = param->nYcells ;
             else
                 ySkip1 = 0;
         }
         ySkip2 = cellY ;
         ySkip3 = cellY + 1 ;
         if (ySkip3 > param->nYcells) {
             if (param->yBoundType == Periodic)
                 ySkip3 = 0 ;
             else
                 ySkip3 = param->nYcells;
         }

         zSkip1 = cellZ - 1 ;
         if (zSkip1 < 0) {
             if (param->zBoundType == Periodic)
                 zSkip1 = param->nZcells;
             else
                 zSkip1 = 0;
         }
         zSkip2 = cellZ ;
         zSkip3 = cellZ + 1 ;
         if (zSkip3 > param->nZcells) {
             if (param->zBoundType == Periodic)
                 zSkip3 = 0 ;
             else
                 zSkip3 = param->nZcells;
         }


	int start = 0;
	if ((param->xBoundType == Periodic) &&
	    (param->yBoundType == Periodic) &&
	    (param->zBoundType == Periodic)) 
	  {			      

             start = 10;//floor(0.5 + (pow(3,lpbc)-1)/2.);
	    if (start > 500) start = 500;
	    printf("For PBC stress calc, start=%d\n",start);
	    
	  }

/*
 *      Loop though all the cells : regular cells go from 1 to 4.
 *      Ghost cells at 0 and 5.
 */
        for (cx = 0; cx <= param->nXcells; cx++) {
	  for (cy = 0; cy <= param->nYcells; cy++) {
	    for (cz = 0; cz <= param->nZcells; cz++) {
	      
	      cellIndex = EncodeCellIdx(home, cx, cy, cz);
	      cell = LookupCell(home, cellIndex);
	      
	      if (cell == (Cell_t *)NULL) continue;

	      xc = cell->center[X];
	      yc = cell->center[Y];
	      zc = cell->center[Z];

	      PBCPOSITION(param, R[X], R[Y], R[Z], &xc, &yc, &zc);

	      ifar = 1;
	      
	      if ((cx==xSkip1 || cx==xSkip2 || cx==xSkip3) &&
		  (cy==ySkip1 || cy==ySkip2 || cy==ySkip3) &&
		  (cz==zSkip1 || cz==zSkip2 || cz==zSkip3)) 
		{
		  ifar = 0;
		}
		    
/*
 *            Loop over all nodes in this cell and over each segment
 *            attached to the node.  Skip any segment that is not
 *            owned by node1.
 */
	      node1 = cell->nodeQ;
	      
	      for (; node1 != (Node_t *)NULL; node1=node1->nextInCell) {
		for (arm = 0; arm < node1->numNbrs; arm++) {
		  
		  node2 = GetNeighborNode(home, node1, arm);
		  
		  if (NodeOwnsSeg(home, node1, node2) == 0) {
		    continue;
		  }

		  p1x = node1->x;
		  p1y = node1->y;
		  p1z = node1->z;
		  
		  p2x = node2->x;
		  p2y = node2->y;
		  p2z = node2->z;

		  PBCPOSITION(param, xc, yc, zc, &p1x, &p1y, &p1z);
		  PBCPOSITION(param, p1x, p1y, p1z, &p2x, &p2y, &p2z);

		  bx = node1->burgX[arm];
		  by = node1->burgY[arm];
		  bz = node1->burgZ[arm];

		  /* Add PBC corrections for N2 stress */
		  for (ixpbc=-start;ixpbc<start+1;ixpbc++)
		    for (iypbc=-start;iypbc<start+1;iypbc++)
		      for (izpbc=-start;izpbc<start+1;izpbc++)
			{
			  if ((param->xBoundType == Periodic) &&
			      (param->yBoundType == Periodic) &&
			      (param->zBoundType == Periodic)) 
			    {			      
			      // Periodic shift
			      pbcp1x = p1x + ixpbc*param->Lx;
			      pbcp1y = p1y + iypbc*param->Ly;
			      pbcp1z = p1z + izpbc*param->Lz;
			      
			      pbcp2x = p2x + ixpbc*param->Lx;
			      pbcp2y = p2y + iypbc*param->Ly;
			      pbcp2z = p2z + izpbc*param->Lz;
			    }
			  else
			    {
			      pbcp1x = p1x;
			      pbcp1y = p1y;
			      pbcp1z = p1z;
			      
			      pbcp2x = p2x;
			      pbcp2y = p2y;
			      pbcp2z = p2z;
			      
			    }

			      if (ifar == 1)
				{
				  /*if (ixpbc == 0 && iypbc == 0 && izpbc == 0)
				    {
				      printf("Far\n");
				      printf("p1=[%f %f %f];\n",p1x,p1y,p1z);
				      printf("p2=[%f %f %f];\n",p2x,p2y,p2z);
				      printf("burg=[%f %f %f];\n\n",bx,by,bz);
				      }*/
				  
				  StressDueToSeg(home, R[X], R[Y], R[Z],
						 pbcp1x, pbcp1y, pbcp1z,
						 pbcp2x, pbcp2y, pbcp2z, 
						 bx, by, bz,
						 a, MU, NU, stress);
				  
                                  sigma[0][0] += stress[0];
                                  sigma[1][1] += stress[1];
                                  sigma[2][2] += stress[2];
                                  
                                  sigma[0][1] += stress[5];
                                  sigma[1][2] += stress[3];
                                  sigma[0][2] += stress[4];
                                  
                                  sigma[1][0] = sigma[0][1];
                                  sigma[2][1] = sigma[1][2];
                                  sigma[2][0] = sigma[0][2];
				  
				}
			      else
				{
				  //printf("Close\n");
				  //printf("p1=[%f %f %f];\n",p1x,p1y,p1z);
				  //printf("p2=[%f %f %f];\n\n",p2x,p2y,p2z);
				  StressDueToSeg(home, R[X], R[Y], R[Z],
						 pbcp1x, pbcp1y, pbcp1z,
						 pbcp2x, pbcp2y, pbcp2z, 
						 bx, by, bz,
						 a, MU, NU, stress);
				  
                                  sigmalocal[0][0] += stress[0];
                                  sigmalocal[1][1] += stress[1];
                                  sigmalocal[2][2] += stress[2];
                                  
                                  sigmalocal[0][1] += stress[5];
                                  sigmalocal[1][2] += stress[3];
                                  sigmalocal[0][2] += stress[4];
                                  
                                  sigmalocal[1][0] += sigmalocal[0][1];
                                  sigmalocal[2][1] += sigmalocal[1][2];
                                  sigmalocal[2][0] += sigmalocal[0][2];
				}

			    } // end PBC images


			      
		} // arm
	      } // node 1
	      
	    }  // cz
	  }  // cy 
	}  // cx
	
        return;
}

#endif




#ifdef ESHELBY // Debugging routine

#if 0
static void N2EshStress(Home_t *home, real8 mu, real8 nu, real8 pos[3], real8 sigma[3][3])
{
   Matrix33_Zero(sigma);
   
   real8 sig[3][3];

   for (int m = 0; m < home->totInclusionCount; m++) 
   {
      EInclusion_t *inclusion = &home->eshelbyInclusions[m];
      int i,j;
      real8 vec[3];
      real8 Ev = inclusion->strainField[0] + inclusion->strainField[1] + inclusion->strainField[2];
      real8 radius = inclusion->radius[0];
   
      vec[0] = pos[0] - inclusion->position[0];
      vec[1] = pos[1] - inclusion->position[1];
      vec[2] = pos[2] - inclusion->position[2];
      
      real8 nvec  = Normal(vec);
      real8 nvec2 = nvec * nvec;
      real8 nvec3 = nvec * nvec2;
      real8 nvec5 = nvec2 * nvec3;
      real8 R3 = radius*radius*radius;
      
      real8 fac0 = mu*(1+nu)*Ev/(1-2*nu);
      real8 fac  = fac0 * R3;
      
      if (nvec > radius)
         for (i = 0; i < 3; i++) 
            for (j = 0; j < 3; j++)
               sig[i][j] = fac * ( (i==j)/nvec3 - 3 * vec[i]*vec[j]/nvec5 );
      else
         for (i = 0; i < 3; i++) 
            for (j = 0; j < 3; j++)
               sig[i][j] = -2*fac0 * (i==j);


           for (i = 0; i<3; i++)
              for (j = 0; j<3; j++)
                 sigma[i][j] += sig[i][j];
   }
}
#endif
#endif


#if 0 //def CALCENERGY // Debugging routine
static void N2Stress(Home_t *home, int cellX, int cellY, int cellZ, 
		     real8 *RemW, real8 *LocW, int lpbc)
{
	int xSkip1, xSkip2, xSkip3;
	int ySkip1, ySkip2, ySkip3;
	int zSkip1, zSkip2, zSkip3;

	int     ifar = 0;
        int     cellIndex, arm;
        int     i,j, c1x, c1y, c1z, c2x, c2y, c2z;
        int     minXIndex, minYIndex, minZIndex;
        int     maxXIndex, maxYIndex, maxZIndex;
	int     ixpbc, iypbc, izpbc;
	double  pbcp1x, pbcp1y, pbcp1z;
	double  pbcp2x, pbcp2y, pbcp2z;	      
        real8   xc, yc, zc;
        real8   a, MU, NU;
        real8   b1x, b1y, b1z, b2x, b2y, b2z;
        real8   p1x, p1y, p1z;
        real8   p2x, p2y, p2z;
        real8   p3x, p3y, p3z;
        real8   p4x, p4y, p4z;
	real8   tin[3], t[3];
        real8   stress[3][3];
        Node_t  *node1, *node2, *node3, *node4;
        Cell_t  *cell1;
        Param_t *param;	

        param = home->param;

        a     = param->rc;
        MU    = param->shearModulus;
        NU    = param->pois;

        *RemW = 0.0;
        *LocW = 0.0;

	/* Get right cells bounds */
         xSkip1 = cellX - 1 ;
         if (xSkip1 < 0) {
             if (param->xBoundType == Periodic)
                 xSkip1 = param->nXcells;
             else
                 xSkip1 = 0;
         }
         xSkip2 = cellX ;
         xSkip3 = cellX + 1 ;
         if (xSkip3 > param->nXcells) {
             if (param->xBoundType == Periodic)
                 xSkip3 = 0 ;
             else
                 xSkip3 = param->nXcells;
         }

         ySkip1 = cellY - 1 ;
         if (ySkip1 < 0) {
             if (param->yBoundType == Periodic)
                 ySkip1 = param->nYcells ;
             else
                 ySkip1 = 0;
         }
         ySkip2 = cellY ;
         ySkip3 = cellY + 1 ;
         if (ySkip3 > param->nYcells) {
             if (param->yBoundType == Periodic)
                 ySkip3 = 0 ;
             else
                 ySkip3 = param->nYcells;
         }

         zSkip1 = cellZ - 1 ;
         if (zSkip1 < 0) {
             if (param->zBoundType == Periodic)
                 zSkip1 = param->nZcells;
             else
                 zSkip1 = 0;
         }
         zSkip2 = cellZ ;
         zSkip3 = cellZ + 1 ;
         if (zSkip3 > param->nZcells) {
             if (param->zBoundType == Periodic)
                 zSkip3 = 0 ;
             else
                 zSkip3 = param->nZcells;
         }

	int startpbc = 0;
	if ((param->xBoundType == Periodic) &&
	    (param->yBoundType == Periodic) &&
	    (param->zBoundType == Periodic)) 
	  {			      
             startpbc = 4;
	  }

/*
 *      Loop though all the cells : regular cells go from 1 to 4.
 *      Ghost cells at 0 and 5.
 */
        for (c1x = 0; c1x <= param->nXcells; c1x++) {
	  for (c1y = 0; c1y <= param->nYcells; c1y++) {
	    for (c1z = 0; c1z <= param->nZcells; c1z++) {
	      
	      cellIndex = EncodeCellIdx(home, c1x, c1y, c1z);
	      cell1 = LookupCell(home, cellIndex);
	      
	      if (cell1 == (Cell_t *)NULL) continue;

	      xc = cell1->center[X];
	      yc = cell1->center[Y];
	      zc = cell1->center[Z];

	      ifar = 1;
	      
	      if ((c1x==xSkip1 || c1x==xSkip2 || c1x==xSkip3) &&
		  (c1y==ySkip1 || c1y==ySkip2 || c1y==ySkip3) &&
		  (c1z==zSkip1 || c1z==zSkip2 || c1z==zSkip3)) 
		{
		  ifar = 0;
		}
		    
/*
 *            Loop over all nodes in this cell and over each segment
 *            attached to the node.  Skip any segment that is not
 *            owned by node1.
 */
	      node1 = cell1->nodeQ;
	      
	      for (; node1 != (Node_t *)NULL; node1=node1->nextInCell) {
		for (arm = 0; arm < node1->numNbrs; arm++) {
		  
		  node2 = GetNeighborNode(home, node1, arm);
		  
		  if (node2 == (Node_t *)NULL) continue;

		  if (NodeOwnsSeg(home, node1, node2) == 0) {
		    continue;
		  }

		  p1x = node1->x;
		  p1y = node1->y;
		  p1z = node1->z;
		  
		  p2x = node2->x;
		  p2y = node2->y;
		  p2z = node2->z;

		  PBCPOSITION(param, xc, yc, zc, &p1x, &p1y, &p1z);
		  PBCPOSITION(param, p1x, p1y, p1z, &p2x, &p2y, &p2z);

		  b1x = node1->burgX[arm];
		  b1y = node1->burgY[arm];
		  b1z = node1->burgZ[arm];

		  /* Add PBC corrections for N2 stress */
		  for (ixpbc=-startpbc;ixpbc<startpbc+1;ixpbc++)
		    for (iypbc=-startpbc;iypbc<startpbc+1;iypbc++)
		      for (izpbc=-startpbc;izpbc<startpbc+1;izpbc++)
			{

			  if ((param->xBoundType == Periodic) &&
			      (param->yBoundType == Periodic) &&
			      (param->zBoundType == Periodic)) 
			    {			      
			      // Periodic shift
			      pbcp1x = p1x + ixpbc*param->Lx;
			      pbcp1y = p1y + iypbc*param->Ly;
			      pbcp1z = p1z + izpbc*param->Lz;
			      
			      pbcp2x = p2x + ixpbc*param->Lx;
			      pbcp2y = p2y + iypbc*param->Ly;
			      pbcp2z = p2z + izpbc*param->Lz;
			    }
			  else
			    {
			      pbcp1x = p1x;
			      pbcp1y = p1y;
			      pbcp1z = p1z;
			      
			      pbcp2x = p2x;
			      pbcp2y = p2y;
			      pbcp2z = p2z;
			      
			    }

			      if (ifar == 1)
				{				  
                                   SegSegEnergy(home, 
						 pbcp1x, pbcp1y, pbcp1z,
						 pbcp2x, pbcp2y, pbcp2z, 
						 b1x, b1y, b1z,
						 0.0, MU, NU, RemW);

				  real8 normStressFMM = 0.0;
				  for (i=0;i<3;i++)
				    for (j=0;j<3;j++)
				      {
					sigma[i][j] += stress[i][j];
					normStressFMM += stress[i][j]*stress[i][j];
				      }
				  fprintf(foutFMM,"%f %f %f %f\n",R[X],R[Y],R[Z],sqrt(normStressFMM));
				  

				  
				}
			      else
				{
				  start = clock();

				  int maxN = 10;
				  for (i=0;i<maxN;i++)
				    StressDueToSeg(home, R[X], R[Y], R[Z],
						   pbcp1x, pbcp1y, pbcp1z,
						   pbcp2x, pbcp2y, pbcp2z, 
						   b1x, b1y, b1z,
						   a, MU, NU, stress);
				  

				  end = clock();
				  mytime = ((double) (end - start)) / CLOCKS_PER_SEC ;
				  home->countN2+= maxN;
				  home->timeN2+=mytime;


				  real8 normStressN2 = 0.0;
				  for (i=0;i<3;i++)
				    for (j=0;j<3;j++)
				      {
					sigmalocal[i][j] += stress[i][j];
					normStressN2 += stress[i][j]*stress[i][j];
				      }

				  fprintf(foutN2,"%f %f %f %f\n",R[X],R[Y],R[Z],sqrt(normStressN2));


				}

			    } // end PBC images


			      
		} // arm
	      } // node 1
	      
	    }  // c1z
	  }  // c1y 
	}  // c1x
	

	fclose(foutN2);
	fclose(foutFMM);
        return;
}

#endif

/****************************************************************************
 *
 *      Function:
 *      Description:  Calculates the positions and weights for conducting
 *                    Gauss-Legendre integration along a normalized span
 *                    going from -1 to 1.  
 *
 *                    WARNING: This function halves the weights before
 *                             returning them to the caller.
 *
 *      Arguments:
 *          numPoints Number of positions/weights to generate.
 *          positions Pointer to array in which to return the positions.
 *                    This array must be at least <intOrder> elements 
 *                    in length.
 *          weights   pointer to array in which to return the weights.
 *                    This array must be at least <intOrder> elements 
 *                    in length.
 *
 ***************************************************************************/
void GaussQuadCoeffHalf(int numPoints, real8 *positions, real8 *weights)
{

        GaussQuadraturePointsAndWgts(numPoints, positions, weights);

        for (int i = 0; i < numPoints; i++) {
            weights[i] *= 0.5; 
        }

        return;
}



#if 0 // Debug FMMs methods
void SegForceFromTaylorExp(Home_t *home, int cellID,
			   real8 *positions, real8 *weights,
                           real8 *p1, real8 *p2, real8 *burg,
                           real8 p1f[3], real8 p2f[3])
{
        int        i, j, numExpansion, lpbc = 0;
	int        ix, iy, iz, numPoints;
	int        xIndex, yIndex, zIndex;
	real8      Lx, Ly, Lz;
	real8      pmidx, pmidy, pmidz;
        real8      pspanx, pspany, pspanz;
        real8      sigbx, sigby, sigbz;
        real8      fLinvx, fLinvy, fLinvz;
        real8      temp, mult1, mult2;
        real8      R[3], evalPos[3], sigmaN2[3][3], sigmalocal[3][3];
	real8      sigmaFMM[3][3], eshelbysigma[3][3];
        real8      cellXsize, cellYsize, cellZsize;

	real8 error, error1=0.0,error2=0.0, error3=0.0;
	real8 relerror, error4=0.0;

        FMCell_t   *cell;
        FMLayer_t  *layer;
        Param_t    *param;

	FILE *fpN2, *fpFMM, *fpError;
	char format[50];


	clock_t start, end;
	double mytime;
	
        param = home->param;
	numExpansion = param->fmExpansionOrder;

#ifdef UNIFORMFMM
	printf("\n\nUNIFORM GRID FMM method\n\n\n");
#endif

#ifdef BBFMM
	printf("\n\nBB FMM method\n\n\n");
#endif


#ifdef TAYLORFMM
	printf("\n\nTAYLOR FMM method\n\n\n");
#endif



	/* Compute the N2 stress to check the FMM */
#ifdef ANISOTROPIC	
	double A;
	AnisotropicVars_t *anisoVars;
	anisoVars = &home->anisoVars;

	A = 2.0*anisoVars->elasticConstantMatrix4D[0][1][0][1] / 
   	       (anisoVars->elasticConstantMatrix4D[0][0][0][0] - 
	        anisoVars->elasticConstantMatrix4D[0][0][1][1]);	

	printf("ANISOTROPIC : Comparison between stress from FMM and N2 interactions \n");
	printf("A=%f\n",A);
#else
	printf("ISOTROPIC : Comparison between stress from FMM and N2 interactions \n");
#endif

	lpbc = 10;
	if  (param->xBoundType == Periodic)
	  printf("PBC are on. For N^2, took lpbc=%d\n",lpbc);
	else
	  printf("PBC are off\n");


	if (lpbc != 0 && param->xBoundType == Free)
	  {
	    printf("\nPBC are off but lpbc is different than 0. N^2 stress will not be correct.\n");
	    printf("Setting lpbc to zero\n\n");
	    lpbc = 0;
	  }

	printf("Hack to compare with Chao's code .......\n");
        layer = &home->fmLayer[param->fmNumLayers-1];

	/* cellID is the cell of node 1 */
	evalPos[X] = -1.42;
	evalPos[Y] =  1.37;
	evalPos[Z] =  1.72;

        Lx = param->Lx;
        Ly = param->Ly;
        Lz = param->Lz;

        xIndex = (int)((evalPos[X] - param->minSideX) / (Lx / param->nXcells));
        yIndex = (int)((evalPos[Y] - param->minSideY) / (Ly / param->nYcells));
        zIndex = (int)((evalPos[Z] - param->minSideZ) / (Lz / param->nZcells));

        cellID = EncodeFMCellIndex(layer->lDim, xIndex, yIndex, zIndex);

        cell  = LookupFMCell(layer->cellTable, cellID);

	printf("cellID = %d\n",cellID);

	R[X] = evalPos[X] - cell->cellCtr[X];
	R[Y] = evalPos[Y] - cell->cellCtr[Y];
	R[Z] = evalPos[Z] - cell->cellCtr[Z];

	LocalToPoint(home, layer->cellSize, R, cell->expansionCoeff, sigmaFMM);
	
	/* Computes the N2 stress at position */
	Print3("Field Point",evalPos);

	printf("Calling N2Stress with xIndex=%d yIndex=%d zIndex=%d\n",xIndex, yIndex, zIndex);
	N2Stress(home, xIndex+1, yIndex+1, zIndex+1, evalPos, 
		 sigmaN2, sigmalocal, lpbc);
	
	Print3x3("sigmaN2",sigmaN2);	   
	Print3x3("sigmaFMM",sigmaFMM);
	
	printf("ratio = \n");
	printf("%.15e %.15e %.15e\n"  , sigmaN2[0][0]/sigmaFMM[0][0], sigmaN2[0][1]/sigmaFMM[0][1], 
	       sigmaN2[0][2]/sigmaFMM[0][2]);
	printf("%.15e %.15e %.15e\n"  , sigmaN2[1][0]/sigmaFMM[1][0], sigmaN2[1][1]/sigmaFMM[1][1], 
	       sigmaN2[1][2]/sigmaFMM[1][2]);
	printf("%.15e %.15e %.15e\n\n", sigmaN2[2][0]/sigmaFMM[2][0], sigmaN2[2][1]/sigmaFMM[2][1], 
	       sigmaN2[2][2]/sigmaFMM[2][2]);
	Print3x3("sigmalocal",sigmalocal);
	
	for (i=0;i<3;i++)
	  for (j=0;j<3;j++)
	    {
	      error1 += (sigmaN2[i][j] - sigmaFMM[i][j])*(sigmaN2[i][j] - sigmaFMM[i][j]);
	      error2 += sigmaN2[i][j]*sigmaN2[i][j];
	    }
	error1 = sqrt(error1);
	error2 = sqrt(error2);
	
	error = error1/error2;
	
	printf("Error =%.15e\n",error);

	printf("\n\n");
	exit(0);
}
#else


/****************************************************************************
 *
 *      Function:     SegForceFromTaylorExp
 *      Description:  This subroutine returns the forces fp1 and fp2 for a
 *                    dislocation segment starting a p1 and ending at p2 with
 *                    burgers vector b calculated from the taylor expansion
 *                    for the cell that encompasses the segment
 *
 *      Arguments:
 *          cell      pointer to the FM cell structure associated with
 *                    the cell owning the segment
 *          positions positions for conducting Gauss-Legendre integration.
 *          weights   weightings for conducting Gauss-Legendre integration.
 *          p1*       coordinates of first endpoint of segment (3 elements)
 *          p2*       coordinates of second endpoint of segment (3 elements)
 *          b*        burgers vector of segment (3 elements)
 *          p1f       Array in which to return forces at point p1 of segment
 *          p2f       Array in which to return forces at point p1 of segment
 *
 ***************************************************************************/
void SegForceFromTaylorExp(Home_t *home, int cellID,
                           real8 *positions, real8 *weights,
                           real8 *p1, real8 *p2, real8 *burg,
                           real8 p1f[3], real8 p2f[3])
{
        int        i, numPoints;
        real8      pmidx, pmidy, pmidz;
        real8      pspanx, pspany, pspanz;
        real8      sigbx, sigby, sigbz;
        real8      fLinvx, fLinvy, fLinvz;
        real8      temp, mult1, mult2;
        real8      R[3], evalPos[3], sigma[3][3];
        FMCell_t   *cell;
        FMLayer_t  *layer;
        Param_t    *param;


        param = home->param;

        layer = &home->fmLayer[param->fmNumLayers-1];
        cell  = LookupFMCell(layer->cellTable, cellID);

        p1f[0] = 0.0;
        p1f[1] = 0.0;
        p1f[2] = 0.0;

        p2f[0] = 0.0;
        p2f[1] = 0.0;
        p2f[2] = 0.0;

/*
 *      If PBC is enabled, the segment endpoints *may* have 
 *      moved outside the periodic boundaries and been folded
 *      back into the far side of the problem space.  In case
 *      this has happened, we need to adjust the coordinates
 *      to that of their periodic images closest to the cell
 *      center used for the taylor expansion.
 */
        PBCPOSITION(param, cell->cellCtr[X], cell->cellCtr[Y],
                    cell->cellCtr[Z], &p1[X], &p1[Y], &p1[Z]);
        PBCPOSITION(param, cell->cellCtr[X], cell->cellCtr[Y],
                    cell->cellCtr[Z], &p2[X], &p2[Y], &p2[Z]);

        numPoints = param->fmNumPoints;

        pmidx  = 0.5 * (p2[X]+p1[X]);
        pmidy  = 0.5 * (p2[Y]+p1[Y]);
        pmidz  = 0.5 * (p2[Z]+p1[Z]);

        pspanx = 0.5 * (p2[X]-p1[X]);
        pspany = 0.5 * (p2[Y]-p1[Y]);
        pspanz = 0.5 * (p2[Z]-p1[Z]);


        for (i = 0; i < numPoints; i++) {

            sigbx = 0.0;
            sigby = 0.0;
            sigbz = 0.0;

            evalPos[X] = pmidx+pspanx*positions[i];
            evalPos[Y] = pmidy+pspany*positions[i];
            evalPos[Z] = pmidz+pspanz*positions[i];

            if (param->fmEnabled) {

                R[X] = evalPos[X] - cell->cellCtr[X];
                R[Y] = evalPos[Y] - cell->cellCtr[Y];
                R[Z] = evalPos[Z] - cell->cellCtr[Z];

		LocalToPoint(home, layer->cellSize, R, cell->expansionCoeff, sigma);

                sigbx += sigma[0][0]*burg[X] +
                         sigma[0][1]*burg[Y] +
                         sigma[0][2]*burg[Z];

                sigby += sigma[1][0]*burg[X] +
                         sigma[1][1]*burg[Y] +
                         sigma[1][2]*burg[Z];

                sigbz += sigma[2][0]*burg[X] +
                         sigma[2][1]*burg[Y] +
                         sigma[2][2]*burg[Z];
            }

#ifdef ESHELBY
            if (param->eshelbyfmEnabled) 
            { 
                real8 eshelbysigma[3][3];

                EshelbyEvalTaylor(cell->eshelbytaylorCoeff, cell->cellCtr,
                                  evalPos, param->eshelbyfmTaylorOrder, eshelbysigma);

                sigbx += eshelbysigma[0][0]*burg[X] +
                         eshelbysigma[0][1]*burg[Y] +
                         eshelbysigma[0][2]*burg[Z];

                sigby += eshelbysigma[1][0]*burg[X] +
                         eshelbysigma[1][1]*burg[Y] +
                         eshelbysigma[1][2]*burg[Z];

                sigbz += eshelbysigma[2][0]*burg[X] +
                         eshelbysigma[2][1]*burg[Y] +
                         eshelbysigma[2][2]*burg[Z];
            }
#endif  // ESHELBY 

            fLinvx = (sigby*pspanz-sigbz*pspany);
            fLinvy = (sigbz*pspanx-sigbx*pspanz);
            fLinvz = (sigbx*pspany-sigby*pspanx);

            temp = weights[i]*positions[i];
            mult1 = weights[i]+temp;

            p2f[0] = p2f[0] + fLinvx*mult1;
            p2f[1] = p2f[1] + fLinvy*mult1;
            p2f[2] = p2f[2] + fLinvz*mult1;

            mult2 = weights[i]-temp;

            p1f[0] = p1f[0] + fLinvx*mult2;
            p1f[1] = p1f[1] + fLinvy*mult2;
            p1f[2] = p1f[2] + fLinvz*mult2;    
        }

        return;
}
#endif

/*---------------------------------------------------------------------------
 *
 *      Function:     RemoteForceOneSeg
 *      Description:  This subroutine calculates the force on a single
 *                    segment from all segments in remote cells.
 *                    IMPORTANT: this function requires that the <node1>
 *                    argument corresponds to the node owning the segment!
 *
 *      Arguments:
 *          node1   Pointer to node owning the segment
 *          node2   Pointer to 2nd node of segment
 *          f1,f2   Arrays in which to return the remote force componenets
 *                  for the node1--node2 segment at the respective nodes.
 *
 *-------------------------------------------------------------------------*/
void RemoteForceOneSeg(Home_t *home, Node_t *node1, Node_t *node2,
                       real8 f1[3], real8 f2[3])
{
        int       armID;
        int       cx, cy, cz, cellID;
        real8     p1[3], p2[3], burg[3];
        real8     p1f[3], p2f[3];
        Param_t   *param;
        FMLayer_t *layer;
        Cell_t    *cell;

        param = home->param;
        layer = &home->fmLayer[param->fmNumLayers-1];

        armID    = GetArmID(node1, node2);

        p1[X] = node1->x;
        p1[Y] = node1->y;
        p1[Z] = node1->z;

        p2[X] = node2->x;
        p2[Y] = node2->y;
        p2[Z] = node2->z;

        burg[X] = node1->burgX[armID];
        burg[Y] = node1->burgY[armID];
        burg[Z] = node1->burgZ[armID];

        PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p2[X], &p2[Y], &p2[Z]);

/*
 *      Find the indices (not shifted for ghost cells) of the cell
 *      owning <node> and convert to the corresponding cellID at
 *      the lowest FM layer.
 */
        cell = LookupCell(home, node1->cellIdx);

        cx = cell->xIndex;
        cy = cell->yIndex;
        cz = cell->zIndex;

        cx--; cy--; cz--;
        cellID = EncodeFMCellIndex(layer->lDim, cx, cy, cz);

        SegForceFromTaylorExp(home, cellID, home->glPositions,
                              home->glWeights, p1, p2, burg, p1f, p2f);

/*
 *      Return the arm-specific forces for both the nodes of the segment,
 *      the total nodal forces will be set to the sum of the arm forces later.
 */
        f1[X] = p1f[X];
        f1[Y] = p1f[Y];
        f1[Z] = p1f[Z];

        f2[X] = p2f[X];
        f2[Y] = p2f[Y];
        f2[Z] = p2f[Z];

        return;
}
