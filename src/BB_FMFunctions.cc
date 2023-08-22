/**************************************************************************
 *
 *      Module:  BlackBoxFMFunctions.c:
 *               Functions specific to the FMM using Black Box FMM.
 *               This file contains the equations of section 3 and 4 of Fong 
 *               and Darve's JCP paper but modified for the anisotropic
 *               kernel.
 *
 *      Includes functions:
 *       bbfmm_S2M()
 *       bbfmm_M2M()
 *       bbfmm_L2L()
 *       bbfmm_M2L()
 *       bbfmm_L2P()
 *
 *      Note: The routines in this file comes from the code bbfmm of Eric 
 *            Darve and Will Fong modified by Chao Chen so that it can handle 
 *            dislocation segments and compute the stress field in anisotropic
 *            elasticity.
 *
 *************************************************************************/
#include "Home.h"

#ifdef BBFMM
#include "FM.h"


void bbfmm_PostM2L(real8 cellSize, real8 *expansionCoeff)
{
  int i, index1, offSet;
  int n = bbfmm->n;
  int n3 = n*n*n;
  int incr = 1;
  int dofn3_f  = 6 * n3;
  int cutoff_f = bbfmm->cutofff;
  real8 beta=0;

  // Scaling factor for SVD 
  real8 scale = pow((1.0/cellSize),bbfmm->homogen);
  real8 alpha=1;
  char *trans = "n";
  real8 *temp;
  
  temp = (real8*) malloc(cutoff_f * sizeof(real8));

  for (i=0; i<cutoff_f; i++) {
    temp[i] = expansionCoeff[i];
  }
  
  dgemv_(trans, &dofn3_f, &cutoff_f, &scale, bbfmm->U,
	 &dofn3_f, temp, &incr, &beta,
	 expansionCoeff, &incr);

  free(temp);

/*
 * Adjust the field values by the appropriate weight
 */
  for (offSet = 0, index1 = 0; index1 < n3; index1++) {
    int   index2;
    real8 tmp;
    
    tmp = bbfmm->Kweights[index1];
    
    for (index2 = 0; index2 < 6; index2++) {
      
      expansionCoeff[offSet] *= tmp;
      offSet++;
    }
  }
}


/*---------------------------------------------------------------------------
 *
 *      Function:    bbfmm_S2M
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
void bbfmm_S2M(real8 p1[3], real8 p2[3], real8 center[3], 
	       real8 burg[3], real8 L, real8* S)
{
  int n,n2,n3,i,j, NumGauss, dofsn3;
  int cutoff_s = bbfmm->cutoffs;
  int l = 0, m = 0;
  int k, w, r, l1, l2, l3, l4, l5, indx;
  real8 sum, tmp, prefac3;
  real8 bxi[9];
  real8 ihalfL;
  real8 midpoint[3], vec[3], xi[3];
  real8 *sourcet[3];
  real8 *Ss[3];


  n = bbfmm->n;
  NumGauss = bbfmm->NumGauss;

  n2 = n*n;  
  n3 = n2*n; 
  prefac3 = 8.0/(real8) n3; 

 // Adjust Length for segments outside the cell 
  L *= (1 + bbfmm->AdjustLen);

  ihalfL = 2.0 / L;

  dofsn3= 9 * n3;

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
    sourcet[0][j] = vec[0] + ihalfL*bbfmm->Gausspoints[j] * xi[0];
    sourcet[1][j] = vec[1] + ihalfL*bbfmm->Gausspoints[j] * xi[1];
    sourcet[2][j] = vec[2] + ihalfL*bbfmm->Gausspoints[j] * xi[2];
  }

  // Compute Ss, the mapping function for the sources
  for (i=0; i<3; i++)  
    Ss[i] = (real8 *)malloc(n*NumGauss * sizeof(real8));

  ComputeSn(sourcet,bbfmm->Tkz,n,NumGauss,Ss); 

  for (r=0;r<3;r++) 
    for (w=0;w<3;w++)  
      { 
	bxi[3*r+w] = burg[w] * xi[r];  
      }

  /* Three loops over 1D Chebychev nodes */
  for (l1=0;l1<n;l1++) 
    for (l2=0;l2<n;l2++) 
      for (l3=0;l3<n;l3++) 
	{ 
	  /* Loop over dofs  */
	  for (i=0;i<9;i++) 
	    { 
	      sum  = 0;
	      for (l5=0;l5<NumGauss;l5++) 
		{ 
		  /* Loop over number of Gauss points : page 8716 Eq. 1 */
		  sum += bbfmm->Gaussweights[l5]*Ss[0][l1*NumGauss+l5]*
		    Ss[1][l2*NumGauss+l5]*Ss[2][l3*NumGauss+l5];  
		}
	      
	      // Multipole coefficients
	      S[l] = prefac3 * sum * bxi[i]; 
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
 *      Function:    bbfmm_M2M
 *      Description: Moment2Moment.
 *                   Shift the center of the provided multipole expansion.
 *                   Used during the upward pass of the FM code.
 *
 *
 *      Arguments:
 *          n:  Chebyshev order
 *          r:  vector from the child box to the parent box
 *       Schild: Source values of the child. Its size is n*n*n*9
 *
 *      Output:
 *       SS: Source values of the parent. Its size is n*n*n*9 too.
 *
 *-------------------------------------------------------------------------*/
void bbfmm_M2M(real8 *r, real8 *Schild, real8 *SS)
{
  int l, l1, l2, l3, l4;
  int count1, count2, count3, wstart;
  int n, incr,n2,n3,dofs,dofn,dofn2,dofn3;
  real8 prefac, prefac3;
  real8 *Cweights;

  n=bbfmm->n;
  n2 = n*n;                    
  n3 = n2*n;                   
  incr = 1;
  dofs  = 9;
  dofn  = 9 * n;
  dofn2 = 9 * n2;
  dofn3 = 9 * n3;
  
  real8 *Sy, *Sz;
  Sy = (real8*) malloc(dofn3 * sizeof(real8));
  Sz = (real8*) malloc(dofn3 * sizeof(real8));
  
  Cweights = bbfmm->Cweights;

  prefac  = 2./(real8)n;
  prefac3 = prefac*prefac*prefac;


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
	      Sz[l]  = ddot_(&n,&Schild[count3],&dofs,&Cweights[count2],&incr);
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
		  SS[l] = prefac3*ddot_(&n,&Sy[count3],&dofn2,&Cweights[count2],&incr);
		  l++;
		}
	    }
	}
    }

  free(Sy);
  free(Sz);
}


/*---------------------------------------------------------------------------
 *      Function:    bbfmm_L2L
 *      Description: Local2Local.
 *                   Shift the center of the provided taylor expansion.
 *                   Used during the downward pass of the FM code.
 *                    
 * Input:
 *       n: Number of Chebyshev nodes
 *       r[3]: Difference between child to parent vector
 *       F: Field values of the parent. Its size is n*n*n*6
 *
 * Output:
 *       Fchild: Field values of the child. Its size is n*n*n*6
 *
 *---------------------------------------------------------------------------*/
void bbfmm_L2L(real8 *r, real8 *F, real8 *Fchild)
{

  int n, n2, n3, doff, dofn, dofn2, dofn3;
  int l, l1, l2, l3, l4, count1, count2, count3, wstart;
  real8 prefac3;
  real8 *Cweights;
  real8 *Fx, *Fy;
  
  n=bbfmm->n;
  Cweights = bbfmm->Cweights;


  n2 = n*n;  
  n3 = n2*n;
  doff  = 6;
  dofn  = 6 * n;
  dofn2 = 6 * n2;
  dofn3 = 6 * n3;
  
  Fx = (real8*) malloc(dofn3 * sizeof(real8));
  Fy = (real8*) malloc(dofn3 * sizeof(real8));

  prefac3 = 8.0/(real8)n3;  
  
  // Recover the index of the child cell
  /* Assumes r is the difference between child to parent */
  
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
	Fchild[l] = prefac3*ddot_(&n,&Fy[count3],&(doff), &Cweights[count2],&n);
	l++;
      }
    }
  }

  free(Fx);
  free(Fy);
  
}



/*---------------------------------------------------------------------------
 *      Function:     bbfmm_L2P
 *      Description : Local2Particle
 *                    Compute stress at a point from a far field expansion.
 *                    Stress is a 6 dimentional vector.
 *
 *      Arguments:
 *
 *                   L : length of the current FMM cell 
 *                   R = field point - FieldCenter
 *                   F: Far Field values at the Chebyshev nodes (refer to the left hand side
 *                      of step 3c at page 8719). Its size is n*n*n*6
 *                   stress: size 3x3
 *-------------------------------------------------------------------------*/
void bbfmm_L2P(real8 R[3], real8 L, real8 *F, real8 stress[3][3])
{
  int i, n, n3, Nf;
  int l, k, l1, l2, l3, l4;
  real8 tmp2, sum, LenAdjust;
  real8 ihalfL, prefac3;
  real8 sigma[6];
  real8 *fieldt[3];
  real8 *Sf[3];

  n=bbfmm->n;
  n3 = n*n*n;


  LenAdjust = (1 + bbfmm->AdjustLen);


  ihalfL   = 2.0/L/LenAdjust;
  prefac3  = 8.0/(double)n3;

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
  ComputeSn(fieldt,bbfmm->Tkz,n,1,Sf);
  
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
      sigma[l4] = prefac3*sum;
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
 *      Function:     bbfmm_M2L
 *      Description : Moment2Local
 *                    Calculate the far field expansion contribution
 *                    for a cell from the multipole expansion of
 *                    a remote cell.
 *
 *      Arguments:
 *         cell_mpCoeff: Array of multipole expansion coefficients. Size n*n*n*9
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
 *                 Size n*n*n*6.
 *
 *-------------------------------------------------------------------------*/
void bbfmm_M2L(real8 *R, real8 L, real8 *cell_mpCoeff, real8 *FFCoeff)
{
  int n=bbfmm->n;
  int n3 = n*n*n, cutoff_f = bbfmm->cutofff, cutoff_s = bbfmm->cutoffs;
  int cutoff2  = cutoff_f * cutoff_s;
  int dofn3_s    = 9 * n*n*n;
                    
  
  real8 iL=1.0/L;

  // Final multipole expansion (Omega w) =  Sw ; v * Sw = Wl
  // 3a in page 8718
  int incr = 1;
  real8  alpha = 1, beta = 0;
  char trans = 'n';
  
  // Determine the corresponding index in the lookup table
  // R = cellSize * ( grid integer source cell i.e. neighboring cell - grid integer field cell i.e. current cell)
  int k1 = (int)(-iL*R[0]) + 3;
  int k2 = (int)(-iL*R[1]) + 3;
  int k3 = (int)(-iL*R[2]) + 3;

  int ninteract = bbfmm->Ktable[49*k1+7*k2+k3];
  int count = ninteract*cutoff2;
  
  /* Multiplication by U */
  // Compute the field proxy values: Matrix Vector multiplication in BLAS: dgemv - 3b
  // cutoff: cutoff on the number of SVD values,

  real8 *K;
  K = bbfmm->K;

  dgemv_(&trans,&cutoff_f,&cutoff_s,&alpha,K+count,&cutoff_f,cell_mpCoeff,&incr,
  	 &beta,FFCoeff,&incr); // 3b  

}

#endif  /* ifdef BBFMM */
