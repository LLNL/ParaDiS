/**************************************************************************
 *
 *      Module:  BlackBoxFMFunctions.c:
 *               Functions specific to the FMM using Black Box FMM.
 *               This file contains the equations of section 3 and 4 of Fong 
 *               and Darve's JCP paper but modified for the anisotropic
 *               kernel.
 *
 *      Includes functions:
 *       uniformfmm_S2M()
 *       uniformfmm_M2M()
 *       uniformfmm_L2L()
 *       uniformfmm_M2L()
 *       uniformfmm_L2P()
 *
 *      Note: The routines in this file comes from the code uniformfmm of Eric 
 *            Darve and Will Fong modified by Chao Chen so that it can handle 
 *            dislocation segments and compute the stress field in anisotropic
 *            elasticity.
 *
 *************************************************************************/
#define _ISOC99_SOURCE
#include <math.h>
#include "Home.h"

#ifdef UNIFORMFMM
#include "FM.h"
#include <fftw/fftw3.h>

/* 
 * FFT transform of multipole coefficients
 * In : mpCoeff of size 9 * n^3
 * Out: mpCoeff of size 9 * (2n-1)^3
 */
void FFTSmallToBig(real8 *mpCoeff)
{
  
  int n, n3, l1, l2, l3, l4;;
  int i, ii, padSize, halfSize,l1Size,l2Pad, l1Pad, shift, count,s;
  int FFTSize;
  
  n = uniformfmm->n;
  n3 = uniformfmm->n3;
  padSize = uniformfmm->pown3;
  halfSize = padSize/2;
  l1Size = n*n; 
  l2Pad = 2*n-1;
  l1Pad = l2Pad * l2Pad; 

  FFTSize     = 9 * padSize;

  real8 *padx    = (real8*)(fftw_malloc(FFTSize*sizeof(real8)));
  real8 *sourcefre = (real8*)(fftw_malloc(FFTSize*sizeof(real8)));

  for (i=0; i< FFTSize; i++) padx[i] = 0.0;

  for (i=0, count=0; i<n3; i++) {
    l3 = i % n;
    l1 = i / l1Size;
    l2 = i / n % n;
    shift = halfSize+l1*l1Pad+l2*l2Pad+l3;

    for (s=0; s<9; s++) {
      padx[shift] = mpCoeff[count++];
      shift += padSize; 
    }
  }

  fftw_plan p[9];
  for (s=0, shift=0; s<9; s++) {
    p[s] = fftw_plan_r2r_1d(padSize, padx + shift, sourcefre + shift, FFTW_R2HC, FFTW_FLAG);
    fftw_execute( p[s] );
    shift += padSize;
  }
  

  for (i = 0; i < FFTSize; i++) {
    mpCoeff[i] = sourcefre[i];
  }

  
  fftw_free(padx), padx=NULL;
  free(sourcefre);
  for (s=0; s<9; s++)
    fftw_destroy_plan(p[s]); 
}

/* 
 * FFT transform of multipole coefficients
 * In : mpCoeff of size 9 * (2n-1)^3
 * Out: mpCoeff of size 9 * n^3
 */
void FFTBigToSmall(real8 cellSize, real8 *expansionCoeff)
{
    int i, j, n, n3, padSize, l1Size,l2Pad, l1Pad;
    int shift, count;
    int f, l1, l2, l3;
    int FFTSize, InArraySize;
    real8 scale;
    real8 *Pf, *tmpexp; //, *res;
    
    
    n = uniformfmm->n;
    n3 = uniformfmm->n3;
    padSize =  uniformfmm->pown3;
    l1Size = n*n;
    l2Pad = 2*n-1;
    l1Pad = l2Pad * l2Pad;
    
    InArraySize = 6 * n3;
    FFTSize = 6 * padSize;  
    
    Pf = (real8 *)calloc(1, FFTSize * sizeof(real8));  
    //res = (real8 *)calloc(1,FFTSize * sizeof(real8));
    tmpexp = (real8 *)calloc(1,FFTSize * sizeof(real8));

    scale = pow((1.0/cellSize),uniformfmm->homogen);
    
    for (i=0; i<FFTSize; i++) { 
      tmpexp[i] = expansionCoeff[i];
    }   
   
    real8 *res = (real8*)(fftw_malloc( padSize*6  *sizeof(real8)));
    fftw_plan p[6];
    for (f=0; f<6; f++)
      {
	//rfftw_one(uniformfmm->p_c2r, tmpexp + f*padSize, res + f*padSize);
        p[f] = fftw_plan_r2r_1d(padSize, tmpexp + f*padSize, res + f*padSize, FFTW_HC2R, FFTW_FLAG);
	fftw_execute(p[f]);
     }
   
    for (f=0; f<6; f++)
      fftw_destroy_plan(p[f]);


    fftw_free(tmpexp), tmpexp=NULL;

    for (i=count=0; i<n*n*n; i++) {
      l3 = i % n;
      l1 = i / l1Size;
      l2 = i / n % n;
      shift = l1*l1Pad+l2*l2Pad+l3;
      
      for (f=0; f<6; f++, shift+=padSize) {
	Pf[count++] = res[shift]/padSize;
      }
    }
		  
    for (j=0; j<InArraySize; j++)
      {
	expansionCoeff[j] = scale *Pf[j];
      }
    
    
    free(Pf);
    fftw_free(res), res = NULL;
    
  }


/*---------------------------------------------------------------------------
 *
 *      Function:    uniformfmm_S2M
 *      Description: Segment to Moment.
 *                   Create a multipole expansion representing the stress
 *                   field from a straight segment [p1, p2]. The expansion 
 *                   is of order n, and the resulting multipole coefficients 
 *                   are stored in the array S, which is of length n*n*n*6
 *
 *      Arguments:
 *          n      Chebyshev order
 *          vec1   vector from the center of expansion, to the begining
 *                 of the line: R-p1
 *          vec2   vector from the beginning of the line to the end of
 *                 the line:  p2-p1
 *          b      Burgers vector.
 *          L      Length of the FMM box
 *          S      multipole coefficients. Its size is n*n*n*9
 *
 *-------------------------------------------------------------------------*/
void uniformfmm_S2M(real8 p1[3], real8 p2[3], real8 center[3], 
		    real8 burg[3], real8 L, real8* S)
{
  int n,n2,n3,i,j, NumGauss;
  int InArraySize;
  real8 ihalfL;
  real8 midpoint[3], vec[3], xi[3], LenAdjust;
  real8 *sourcet[3];
  real8 *Ss[3];

  n = uniformfmm->n;
  n2 = uniformfmm->n3; 
  n3 = uniformfmm->n3; 
  NumGauss = uniformfmm->NumGauss;

  // Adjust Length for segments outside the cell 
  L *= (1 + uniformfmm->AdjustLen);


  ihalfL = 2.0 / L;
  InArraySize = 9 * n3;

/*
 * Use Gauss quadrature to approximate the integral over the dislocation segment
 */

  for (i=0; i<3; i++)  
    sourcet[i] = (real8 *)malloc(NumGauss * sizeof(real8));

  midpoint[0] = (p1[0]+p2[0])*0.5;
  midpoint[1] = (p1[1]+p2[1])*0.5;
  midpoint[2] = (p1[2]+p2[2])*0.5;

  vec[0] = ihalfL*(midpoint[0] - center[0]);
  vec[1] = ihalfL*(midpoint[1] - center[1]);
  vec[2] = ihalfL*(midpoint[2] - center[2]);

  xi[0] = (p2[0] - p1[0])*0.5;
  xi[1] = (p2[1] - p1[1])*0.5;
  xi[2] = (p2[2] - p1[2])*0.5;

  // Map the source points to the box ([-1, 1])^3
  for (j=0;j<NumGauss;j++) {
    sourcet[0][j] = vec[0] + ihalfL*uniformfmm->Gausspoints[j] * xi[0];
    sourcet[1][j] = vec[1] + ihalfL*uniformfmm->Gausspoints[j] * xi[1];
    sourcet[2][j] = vec[2] + ihalfL*uniformfmm->Gausspoints[j] * xi[2];
  }

  // Compute Ss, the mapping function for the sources
  for (i=0; i<3; i++)  
    Ss[i] = (real8 *)malloc(n*NumGauss * sizeof(real8));
  ComputeSnUniform(sourcet,n,NumGauss,Ss); 

  int l = 0, m = 0;
  int k, w, r, l1, l2, l3, l4, l5, indx;
  real8 sum, tmp;
  real8 bxi[9];

  for (r=0;r<3;r++) 
    for (w=0;w<3;w++)  
      { 
	bxi[3*r+w] = burg[w] * xi[r];  
      }

  /* Three loops over 1D Uniform grid nodes */
  for (l1=0;l1<n;l1++) 
    for (l2=0;l2<n;l2++) 
      for (l3=0;l3<n;l3++) 
	{ 
	  /* Loop over 9  */
	  for (i=0;i<9;i++) 
	    { 
	      sum  = 0;
	      for (l5=0;l5<NumGauss;l5++) 
		{ /* Loop over number of Gauss points : page 8716 Eq. 1 */
		  sum += uniformfmm->Gaussweights[l5]*Ss[0][l1*NumGauss+l5]*
		    Ss[1][l2*NumGauss+l5]*Ss[2][l3*NumGauss+l5];  
		}
	      
	      // Multipole coefficients
	      S[l] = sum * bxi[i]; 
	      l++;
	    } // end w, r
	} // end l1, l2, l3
  
   for (i=0;i<3;i++) 
    {
      free(Ss[i]);
      free(sourcet[i]);
    }


}


/*---------------------------------------------------------------------------
 *
 *      Function:    uniformfmm_M2M
 *      Description: Moment2Moment.
 *                   Shift the center of the provided multipole expansion.
 *                   Used during the upward pass of the FM code.
 *
 *
 *      Arguments:
 *          n:  Uniform grid order
 *          r:  vector from the child box to the parent box
 *       Schild: Source values of the child. Its size is n*n*n*9
 *
 *      Output:
 *       SS: Source values of the parent. Its size is n*n*n*9 too.
 *
 *-------------------------------------------------------------------------*/
void uniformfmm_M2M(real8 *r, real8 *Schild, real8 *SS)
{
  int l, l1, l2, l3, l4, dofs;
  int count1, count2, count3, wstart;

  int n, incr, n2, n3, dofn, dofn2, InArraySize;
  
  real8 prefac;
  real8 *Cweights;

  incr = 1;

  n=uniformfmm->n;
  n2 = n*n;   
  n3 = uniformfmm->n3;
  dofs  = 9;
  dofn  = 9 * n;
  dofn2 = 9 * n2;
  InArraySize = 9 * n3;

  real8 *Sy, *Sz;
  Sy    = (real8 *)malloc(InArraySize * sizeof(real8));
  Sz    = (real8 *)malloc(InArraySize * sizeof(real8));

  Cweights = uniformfmm->Cweights;

  prefac  = 2./(real8)n;

  /* Recover the index of the child box 
     Careful here: assumes that r is the
     difference between child to parent */

  /* Gather the children source along the z-component */
  l = 0;
  if (r[2] < 0)
    wstart =  0;
  else
    wstart =  n2;

  for (l1=0;l1<n2;l1++) 
    {
      count1 = l1*n;
      for (l3=0;l3<n;l3++) 
	{
	  count2 = wstart + l3*n;
	  for (l4=0;l4<9;l4++) 
	    {
	      count3 = 9*count1 + l4;
	      Sz[l]  = ddot_(&n,&Schild[count3],&(dofs),&Cweights[count2],&incr);
	      l++;
	    }
	}
    }
  
  /* Gather the children sources along the y-component */
  l = 0;
  if (r[1] < 0)
    wstart =  0;
  else
    wstart =  n2;
  
  for (l1=0;l1<n;l1++) 
    {
      for (l2=0;l2<n;l2++) 
	{
	  count2 = wstart + l2*n;
	  for (l3=0;l3<n;l3++) 
	    {
	      count1 = l1*n2 + l3;
	      for (l4=0;l4<9;l4++) 
		{
		  count3 = 9*count1 + l4;
		  Sy[l]  = ddot_(&n,&Sz[count3],&dofn,&Cweights[count2],&incr);
		  l++;
		}
	    }
	}
    }
  
  /* Gather the children sources along the x-component and determine
     the parent sources */
  l = 0;
  if (r[0] < 0)
    wstart =  0;
  else
    wstart =  n2;
  
  for (l1=0;l1<n;l1++) 
    {
      count2 = wstart + l1*n;
      for (l2=0;l2<n;l2++) 
	{
	  for (l3=0;l3<n;l3++) 
	    {
	      count1 = l2*n + l3;
	      for (l4=0;l4<9;l4++) 
		{
		  count3  = 9*count1 + l4;
		  SS[l] = ddot_(&n,&Sy[count3],&dofn2,&Cweights[count2],&incr);
		  l++;
		}
	    }
	}
    }

  free(Sy);
  free(Sz);

}


/*---------------------------------------------------------------------------
 *      Function:    uniformfmm_L2L
 *      Description: Local2Local.
 *                   Shift the center of the provided taylor expansion.
 *                   Used during the downward pass of the FM code.
 *                    
 * Input:
 *       n: Number of Uniform grid nodes
 *       r[3]: Difference between child to parent vector
 *       F: Field values of the parent. Its size is n*n*n*6
 *
 * Output:
 *       Fchild: Field values of the child. Its size is n*n*n*6
 *
 *---------------------------------------------------------------------------*/
void uniformfmm_L2L(real8 *r, real8 *F, real8 *Fchild)
{
        int l, l1, l2, l3, l4, count1, count2, count3, wstart;
	int n, n2, n3, dofn, dofn2, InSizeArray,doff;
	real8 *Cweights;

	n=uniformfmm->n;
	n2 = n*n;                
	n3 = uniformfmm->n3;     
	doff  = 6;
	dofn  = 6 * n;
	dofn2 = 6 * n2;
	InSizeArray = 6 * n3;

	
	real8 *Fx, *Fy;
	Fx    = (real8 *)malloc(InSizeArray * sizeof(real8));
	Fy    = (real8 *)malloc(InSizeArray * sizeof(real8));
  

	Cweights = uniformfmm->Cweights;
		
/*
 *      Recover the index of the child box
 *      Assumes r is the difference between child to parent 
 */

	// Interpolate the parent field along the x-component
        l = 0;
	if (r[0] < 0)
	  wstart = 0;
	else
	  wstart = n2;
	
	for (l1=0;l1<n;l1++) {
	  count2 = wstart + l1;
	  for (l2=0;l2<n;l2++) {
	    for (l3=0;l3<n;l3++) {
	      count1 = l2*n + l3;
	      for (l4=0;l4<6;l4++) {
		count3 = 6*count1 + l4;
		Fx[l]  = ddot_(&n,&F[count3],&dofn2,&Cweights[count2],&n);
		l++;
	      }
	    }
	  }
	}
		
	// Interpolate the parent field along the y-component
	l = 0;
	if (r[1] < 0)
	  wstart = 0;
	else
	  wstart = n2;
	
	for (l1=0;l1<n;l1++) {
	  for (l2=0;l2<n;l2++) {
	    count2 = wstart + l2;
	    for (l3=0;l3<n;l3++) {
	      count1 = l1*n2 + l3;
	      for (l4=0;l4<6;l4++) {
		count3 = 6*count1 + l4;
		Fy[l] = ddot_(&n,&Fx[count3],&dofn,&Cweights[count2],&n);
		l++;
	      }
	    }
	  }
	}
	
	/* Interpolate the parent field along the z-component and add
	   to child field */
	//j = zindex;
	l = 0;
	if (r[2] < 0)
	  wstart = 0;
	else
	  wstart = n2;
	
	for (l1=0;l1<n2;l1++) {
	  count1 = l1*n;
	  for (l3=0;l3<n;l3++) {
	    count2 = wstart + l3;
	    for (l4=0;l4<6;l4++) {
	      count3 = 6*count1 + l4;
	      Fchild[l] = ddot_(&n,&Fy[count3],&(doff), &Cweights[count2],&n);
	      l++;
	    }
	  }
	}


	free(Fx);
	free(Fy);
}



/*---------------------------------------------------------------------------
 *      Function:     uniformfmm_L2P
 *      Description : Local2Particle
 *                    Compute stress at a point from a far field expansion.
 *                    Stress is a 6 dimentional vector.
 *
 *      Arguments:
 *
 *                   L : length of the current FMM cell 
 *                   R = field point - FieldCenter
 *                   F: Far Field values at the Uniform grid nodes (refer to the left hand side
 *                      of step 3c at page 8719). Its size is n*n*n*6
 *                   stress: size 3x3
 *-------------------------------------------------------------------------*/
void uniformfmm_L2P(real8 R[3], real8 L, real8 *F, real8 stress[3][3])
{
  int i, n, n3;
  int l, k, l1, l2, l3, l4;
  real8 tmp2, sum;
  real8 ihalfL, LenAdjust;
  real8 sigma[6];
  real8 *fieldt[3], *Sf[3];

  n  = uniformfmm->n;
  n3 = uniformfmm->n3;

  LenAdjust = (1 + uniformfmm->AdjustLen);

  ihalfL   = 2.0/L/LenAdjust; 

  for (i=0;i<3;i++) 
    {
      Sf[i]    = (real8 *)malloc(n * sizeof(real8));
      fieldt[i]= (real8 *)malloc(1 * sizeof(real8));
    }

  // Map all of the field points to the box ([-1 1])^3
  fieldt[0][0] = ihalfL*R[0];
  fieldt[1][0] = ihalfL*R[1];
  fieldt[2][0] = ihalfL*R[2];
  
  // Compute Sf, the mapping function for the field points
  ComputeSnUniform(fieldt,n,1,Sf);

  // Compute the values at the field point
  for (l4=0;l4<6;l4++) 
    {
      // Due to far field interactions
      sum = 0;
      l = l4;
      for (l1=0;l1<n;l1++) 
	for (l2=0;l2<n;l2++)
	  {
	    tmp2 = Sf[0][l1]*Sf[1][l2];
	    for (l3=0;l3<n;l3++) 
	      {
		sum += F[l]*tmp2*Sf[2][l3];
		l   += 6;
	      }
	  } // end l1, l2, l3
      sigma[l4] = sum;
    }
  
  for (i=0;i<3;i++) 
    {
      free(Sf[i]);
      free(fieldt[i]);
    }

    
  stress[0][0] = sigma[0];
  stress[1][1] = sigma[1];
  stress[2][2] = sigma[2];
  stress[0][1] = sigma[5];
  stress[1][2] = sigma[3];
  stress[0][2] = sigma[4];

  stress[1][0] = stress[0][1];
  stress[2][0] = stress[0][2];
  stress[2][1] = stress[1][2];
}

/*---------------------------------------------------------------------------
 *      Function:     uniformfmm_M2L
 *      Description : Moment2Local
 *                    Calculate the far field expansion contribution
 *                    for a cell from the multipole expansion of
 *                    a remote cell.
 *
 *      Arguments:
 *         cell_mpCoeff: Array of multipole expansion coefficients. Size (2n-1)^3*9
 *
 *	   K: (Precomputed) Pointer to the array of all the compressed M2L operators 
 *	      (refer to step 0 at the page of 8718)
 *
 *         Ktable: (Precomputed) Pointer to the array of indices of M2L operators, used to 
 *            find a specific compressed M2L operator in 'K' above
 *         L : length of the cell
 *         R : vector from the center of the 
 *             current cell to the center of the multipole
 *             expansion in the neighboring cell.
 *
 * Output:
 *        FFCoeff: Array in which the far field expansion coefficients are stored.
 *                 Size (2n-1)^3 * 6.
 *
 *-------------------------------------------------------------------------*/
void uniformfmm_M2L(real8 *R, real8 L, real8 *cell_mpCoeff, real8 *FFCoeff)
{
  int i,n, n3, pown3;
  int l, l1, l2;
  int k1, k2, k3;
  int FFTSize, InArraySize;
  int ninteract, count;
  int N, f, s, shift1=0, shift2=0, shift3=0;
  real8 iL;

  n=uniformfmm->n;
  n3 = uniformfmm->n3;
  pown3 = uniformfmm->pown3;

  //printf("n=%d n3=%d pown3=%d\n",n,n3,pown3);

  FFTSize = pown3 *  6;
  InArraySize = n3 *  6;
  
  // Not initialized in input. Do it here.
  for (i=0;i< FFTSize;i++) FFCoeff[i] = 0.0;
             
  // Final multipole expansion (Omega w) =  Sw ; v * Sw = Wl
  l = 0;
  iL=1.0/L;

  // Determine the corresponding index in the lookup table
  // R = cellSize * ( grid integer source cell i.e. neighboring cell - grid integer field cell i.e. current cell)
  k1 = (int)(-iL*R[0]) + 3;
  k2 = (int)(-iL*R[1]) + 3;
  k3 = (int)(-iL*R[2]) + 3;

  ninteract = uniformfmm->Ktable[49*k1+7*k2+k3];
  count = ninteract*pown3*54;
  
  // entry-wise product of cell_mpCoeff
  N = pown3;
  for (f=0; f< 6; f++, shift2+=N)
    for (s=shift1=0; s< 9; s++) {
      FrequencyProduct(N, &uniformfmm->K[count+shift3], &cell_mpCoeff[shift1],
		       &FFCoeff[shift2]);
      shift1 += N;
      shift3 += N;
    }

}

#endif  /* ifdef UNIFORMFMM */





