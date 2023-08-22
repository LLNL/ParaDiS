/**************************************************************************
 *
 *      Module:       FMSigma2NonRational.c
 *      Description:  This module contains the base FMSigma2NonRational
 *                    entry function used by the fast multipole code,
 *                    a generic function for computing the stress
 *                    values for any arbitrary expansion order, and
 *                    special functions for handling expansion
 *                    orders of 3 or less.
 *
 *                    The code currently provides two variants for
 *                    calculating forces. With the original, there
 *                    is a divergence of the stress field on an element by
 *                    element basis that cancels when the total segment
 *                    configuration closes on itself or is bounded by semi-
 *                    infinite segments.  The second method (referred to
 *                    as 'rational' segment forces introduces a correction
 *                    that makes the stress field due to a dislocation
 *                    segment divergence free on an element by element
 *                    basis so the system does not have to be closed or
 *                    terminated by semi-infinite segments in order to 
 *                    obtain the divergence-free stress field.
 *
 *      Includes public functions:
 *          FMSigma2NonRational()
 *
 *      Includes private functions:
 *          FMSigma2core()
 *          FMSigma2core0()
 *          FMSigma2core1()
 *          FMSigma2core2()
 *          FMSigma2core3()
 *
 *************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"

/*
 *      Defining PRINT_PE will force the FMSigma2core function
 *      to print out in C-code the operations it is performing.
 *      This is the mechanism that was used to initially 
 *      generate functions specific to the expansion orders with
 *      lots of unrolled loops and simple code the optimizers can
 *      handle well.
 */
#ifdef PRINT_PE
#undef PRINT_PE
#endif

#ifdef PRINT_PE
#define xprintf printf
#else
#define xprintf (void)
#endif

/* Unfortunately - defining xprintf() this way causese GCC to barf
 * out a tremendous volume of unused value warnings. The following
 * pragma diables those warnings.
 */

#pragma GCC diagnostic ignored "-Wunused-value"

/*
 *      Prototype the basic unoptimized function to handle
 *      any arbitrary expansion order, as well as any
 *      expansion-order specific optimized functions.
 */
static void FMSigma2core(int iorder, int pows[][3], real8 terms[],
                         real8 mu8pi, real8 two1nu, real8 Eeta[], matrix sigma);
static void FMSigma2core0(real8 terms[],real8 mu8pi, real8 two1nu,
                          real8 Eeta[],matrix sigma);
static void FMSigma2core1(real8 terms[],real8 mu8pi, real8 two1nu,
                          real8 Eeta[],matrix sigma);
static void FMSigma2core2(real8 terms[],real8 mu8pi, real8 two1nu,
                          real8 Eeta[],matrix sigma);
static void FMSigma2core3(real8 terms[],real8 mu8pi, real8 two1nu,
                          real8 Eeta[],matrix sigma);


/*---------------------------------------------------------------------------
 *
 *      Function:    FMSigma2NonRational
 *      Description: This routine is a little more than a dispatching
 *                   function which does some initialization and invokes
 *                   other special functions to do the real work.  The
 *                   special functions are called fmsigmacore<N> where <N>
 *                   is the value of norder.  Special routines are provided
 *                   for orders up to 5; orders higher than that will
 *                   call the generic function FMSigma2core() that
 *                   can handle higher orders.
 *
 *                   NOTE: The specialized functions can be regenerated
 *                   by defining the PRINT_PE macro at the start of this
 *                   moudule and invoking the FMSigma2core() function.
 *                   See code for more details.
 *
 *      Arguments:
 *          mu       Shear modulus
 *          nu       Poisson ratio
 *          norder   order of the multipole expansion
 *          Eeta     should be an array resulting from calls to makeeta, and
 *                   should have 9*(norder+1)*(noder+2)*(norder+3)/6 elements
 *                   (this is always an integer, and evaluated as written,
 *                   yields the correct result with integer arithmetic).
 *          sigmatot The resulting stress is written into sigmatot.
 *          pcterms
 *          pcpows
 *
 *-------------------------------------------------------------------------*/

void FMSigma2NonRational(real8 mu, real8 nu, int norder, real8 Eeta[],
                         matrix sigmatot, real8 pcterms[], int pcpows[][3])
{
        int    iorder, i, j, (*pows)[3], npows;
        real8  *terms;
        real8  pi = M_PI;
        real8  mu8pi = mu/(8.0*pi);
        real8  two1nu = 2.0/(1.0-nu);  
        matrix sigma;

        static real8 fact[NMAX+1],ifact[NMAX+1],dfact[2*NMAX-3+2];
        static int inited = 0;

/*
 *      Various arrays are static once initialized and shared among
 *      threads, but we have to make sure they are only initialized
 *      only once no matter how many threads are active.
 */
        if(inited == 0) {
#ifdef _OPENMP
#pragma omp critical (CRIT_INIT_FMSIGMA2)
#endif
            {
                if (inited == 0) {
                    makeftabs(fact,ifact,dfact);
                    inited = 1;
                }
            }
        }

        for(i = 0; i<3; i++) {
            for(j = 0; j<3; j++) {
                sigmatot[i][j] = 0.0;
            }
        }

        terms = pcterms;
        pows = pcpows;

        for(iorder=0; iorder<=norder; iorder++) {

            npows = (iorder+1+3)*(iorder+2+3)/2;

            switch(iorder) {
            case 0:
                FMSigma2core0(terms,mu8pi,two1nu,Eeta,sigma);
                break;
            case 1:
                FMSigma2core1(terms,mu8pi,two1nu,Eeta,sigma);
                break;
            case 2:
                FMSigma2core2(terms,mu8pi,two1nu,Eeta,sigma);
                break;
            case 3:
                FMSigma2core3(terms,mu8pi,two1nu,Eeta,sigma);
                break;
            default:
                FMSigma2core(iorder,pows,terms,mu8pi,two1nu,Eeta,sigma);
            }
        
            for(i = 0; i<3; i++) {
                for(j = 0; j<3; j++) {
                    sigmatot[i][j] = sigmatot[i][j] + sigma[i][j];
                }
            }

/*
 *          NB! Pointer arithmetic
 */
            terms += npows;
            pows += npows;
        }

        return;
}



static void FMSigma2core(int iorder, int pows[][3], real8 terms[],
		  real8 mu8pi, real8 two1nu, real8 Eeta[], matrix sigma) {
  int cyc[] = {0,1,2,0,1,2};
  int i,j,k,l,m,n,a,nx,ny,nz,npows,nidx;
  int etaoff,etadis0;
  matrix Gijmkl[3],Gmppkl[3] = {{{0.0}}};
  real8 t,t2,t11,t12,t3,tvec[NTMAX];

  static real8 fact[NMAX+1],ifact[NMAX+1],dfact[2*NMAX-3+2];
  static int inited = 0;

  if (inited == 0) {
#ifdef _OPENMP
#pragma omp critical (CRIT_INIT_FMSIGMA2)
#endif
    {
      if (inited == 0) {
        makeftabs(fact,ifact,dfact);
        inited = 1;
      }
    }
  }

  xprintf("#include \"real.h\"\n");
  xprintf("#ifndef NTMAX\n");
  xprintf("#define NTMAX  ( (NMAX+1)*(NMAX+2)/2 )\n");
  xprintf("#endif\n\n");

  xprintf("void FMSigma2core%d(real8 terms[],real8 mu8pi, "
	 "real8 two1nu,real8 Eeta[],matrix sigma) {\n",iorder);
  xprintf("  matrix Gijmkl[3],Gmppkl[3] = {{{0.0}}};\n");
  xprintf("  real8 t,t11,t12,t2,t3,tvec[NTMAX];\n");
  xprintf("  int npows,etaoff,a,nidx;\n\n");
  
  npows = (iorder+1+3)*(iorder+2+3)/2;
  xprintf("  npows = %d;\n",npows);

  etadis0 = (iorder+1)*(iorder+2) >> 1; /* (iorder+1)*(iorder+2)/2          */
  etaoff = 3*iorder*etadis0;            /* 9*iorder*(iorder+1)*(iorder+2)/6 */
  xprintf("  etaoff = %d;\n",etaoff);
  for(i=0; i<3; i++)
    for(j=0; j<3; j++) {
      for(m=0; m<3; m++) {
	/* From 'terms', reconstruct R(i,j,m,a1,...,aq) */
	int d[] = {0,0,0};
	d[i]++; d[j]++; d[m]++;
	nidx = 0;
	for(k = 0; k<npows; k++) {
	  nx = pows[k][0] - d[0];
	  ny = pows[k][1] - d[1];
	  nz = pows[k][2] - d[2];
	  if(nx>=0 && ny>=0 && nz>=0) {
	    tvec[nidx++] = terms[k] * fact[iorder] *
	      ifact[nx]*ifact[ny]*ifact[nz];
	    xprintf("  tvec[%3d] = terms[%3d] * %10.1f;\n",
		   nidx-1,k,fact[iorder]/fact[nx]/fact[ny]/fact[nz]);
	  }
	}
	xprintf("\n  nidx = %d;\n",nidx);

	etadis0 = etaoff;
	for(k = 0; k<3; k++)
	  for(l = 0; l<3; l++) {
	    t = 0.0;
	    xprintf("  t = 0.0;\n");
	    for(a = 0; a<nidx; a++) {
	      t = t + tvec[a]*Eeta[etadis0++];
	    }
	    xprintf("  for(a = 0; a<%d; a++)\n  ",nidx);
	    xprintf("  t = t + tvec[a]*Eeta[%3d+a];\n",etaoff+(3*k+l)*nidx);
	    
	    Gijmkl[m][k][l] = t;
	    xprintf("  Gijmkl[%3d][%3d][%3d] = t;\n",m,k,l);
	    if(i == j) {
	      Gmppkl[m][k][l] = Gmppkl[m][k][l] + t;
	      xprintf("  Gmppkl[%3d][%3d][%3d] = Gmppkl[%3d][%3d][%3d] + t;\n",
		     m,k,l,m,k,l);
	    }
	  }
      }
      
      t2 = 0.0;
      xprintf("\n  t2 = 0.0;\n");
      for(k = 0; k<3; k++) {
	/* Second term: 2/(1-nu) * ijk(k,m,n)G(i,j,m,n,k) */
	m = cyc[k+1]; n = cyc[k+2];
	t2 = t2 + (Gijmkl[m][n][k] - Gijmkl[n][m][k]);
	xprintf("  t2 = t2 + (Gijmkl[%3d][%3d][%3d] - Gijmkl[%3d][%3d][%3d]);\n",
	       m,n,k,n,m,k);
      }
      sigma[i][j] = two1nu * t2;
      xprintf("  sigma[%3d][%3d] = two1nu * t2;\n",i,j);
    }
  
  /* Term one: ijk(j,m,n)G(m,p,p,n,i) + ijk(i,m,n)G(m,p,p,n,j) */
  xprintf("\n");
  for(i = 0; i<3; i++)
    for(j = 0; j<3; j++) {
      m = cyc[j+1]; n = cyc[j+2];
      t11 = Gmppkl[m][n][i] - Gmppkl[n][m][i];
      xprintf("  t11 = Gmppkl[%3d][%3d][%3d] - Gmppkl[%3d][%3d][%3d];\n",
	     m,n,i,n,m,i);
      
      m = cyc[i+1]; n = cyc[i+2];
      t12 = Gmppkl[m][n][j] - Gmppkl[n][m][j];
      xprintf("  t12 = Gmppkl[%3d][%3d][%3d] - Gmppkl[%3d][%3d][%3d];\n",
	     m,n,j,n,m,j);
      
      sigma[i][j] = sigma[i][j] + (t11 + t12);
      xprintf("  sigma[%3d][%3d] = sigma[%3d][%3d] + (t11 + t12);\n",
	     i,j,i,j);
    }
  
  /* Third term: 2/(1-nu) * delta(i,j)ijk(km,n)G(p,p,m,n,k) */
  t3 = 0.0;
  xprintf("\n  t3 = 0.0;\n");
  for(k = 0; k<3; k++) {
    m = cyc[k+1]; n = cyc[k+2];
    t3 = t3 + Gmppkl[m][n][k] - Gmppkl[n][m][k];
    xprintf("  t3 = t3 + Gmppkl[%3d][%3d][%3d] - Gmppkl[%3d][%3d][%3d];\n",
	   m,n,k,n,m,k);
  }
  t3 = t3 * two1nu;
  xprintf("  t3 = t3 * two1nu;\n");
  for(i = 0; i<3; i++) {
    sigma[i][i] = sigma[i][i] - t3;
    xprintf("sigma[%3d][%3d] = sigma[%3d][%3d] - t3;\n",i,i,i,i);
  }
  t = ipow(-1.0,iorder);
  xprintf("\n  t = %30.20e * mu8pi;\n",ipow(-1.0,iorder) * ifact[iorder]);
  for(i = 0; i<3; i++)
    for(j = 0; j<3; j++) {
      sigma[i][j] = t * sigma[i][j] * mu8pi*ifact[iorder];
      xprintf("  sigma[%3d][%3d] = t * sigma[%3d][%3d];\n",i,j,i,j);
    }
  xprintf("}\n\n\n\n");
}


void FMSigma2core0(real8 terms[],real8 mu8pi, real8 two1nu,real8 Eeta[],matrix sigma) {
real8 H000,H001,H002,H010,H011,H012,H020,H021,H022;
real8 H100,H101,H102,H110,H111,H112,H120,H121,H122;
real8 H200,H201,H202,H210,H211,H212,H220,H221,H222;
real8 G000,G001,G002,G010,G011,G012,G020,G021,G022;
real8 G100,G101,G102,G110,G111,G112,G120,G121,G122;
real8 G200,G201,G202,G210,G211,G212,G220,G221,G222;
real8 t,t11,t12,t2,t3;
real8 tvec0;

  t2 = 0.0;
  tvec0 = terms[  0] *        1.0;

  G000 = tvec0*Eeta[  0];
  G001 = tvec0*Eeta[  1];
  G002 = tvec0*Eeta[  2];
  G010 = tvec0*Eeta[  3];
  G011 = tvec0*Eeta[  4];
  G012 = tvec0*Eeta[  5];
  G020 = tvec0*Eeta[  6];
  G021 = tvec0*Eeta[  7];
  G022 = tvec0*Eeta[  8];

  H000 = G000;
  H001 = G001;
  H002 = G002;
  H010 = G010;
  H011 = G011;
  H012 = G012;
  H020 = G020;
  H021 = G021;
  H022 = G022;

  tvec0 = terms[  1] *        1.0;

  G100 = tvec0*Eeta[  0];
  G101 = tvec0*Eeta[  1];
  G102 = tvec0*Eeta[  2];
  G110 = tvec0*Eeta[  3];
  G111 = tvec0*Eeta[  4];
  G112 = tvec0*Eeta[  5];
  t2 = t2 + (G012 - G102);
  G120 = tvec0*Eeta[  6];
  G121 = tvec0*Eeta[  7];
  G122 = tvec0*Eeta[  8];

  H100 = G100;
  H101 = G101;
  H102 = G102;
  H110 = G110;
  H111 = G111;
  H112 = G112;
  H120 = G120;
  H121 = G121;
  H122 = G122;

  tvec0 = terms[  4] *        1.0;

  G200 = tvec0*Eeta[  0];
  G201 = tvec0*Eeta[  1];
  G202 = tvec0*Eeta[  2];
  G210 = tvec0*Eeta[  3];
  t2 = t2 + (G120 - G210);
  G211 = tvec0*Eeta[  4];
  G212 = tvec0*Eeta[  5];
  G220 = tvec0*Eeta[  6];
  G221 = tvec0*Eeta[  7];
  t2 = t2 + (G201 - G021);
  G222 = tvec0*Eeta[  8];

  H200 = G200;
  H201 = G201;
  H202 = G202;
  H210 = G210;
  H211 = G211;
  H212 = G212;
  H220 = G220;
  H221 = G221;
  H222 = G222;

  sigma[  0][  0] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  1] *        1.0;

  G000 = tvec0*Eeta[  0];
  G001 = tvec0*Eeta[  1];
  G002 = tvec0*Eeta[  2];
  G010 = tvec0*Eeta[  3];
  G011 = tvec0*Eeta[  4];
  G012 = tvec0*Eeta[  5];
  G020 = tvec0*Eeta[  6];
  G021 = tvec0*Eeta[  7];
  G022 = tvec0*Eeta[  8];

  tvec0 = terms[  2] *        1.0;

  G100 = tvec0*Eeta[  0];
  G101 = tvec0*Eeta[  1];
  G102 = tvec0*Eeta[  2];
  G110 = tvec0*Eeta[  3];
  G111 = tvec0*Eeta[  4];
  G112 = tvec0*Eeta[  5];
  t2 = t2 + (G012 - G102);
  G120 = tvec0*Eeta[  6];
  G121 = tvec0*Eeta[  7];
  G122 = tvec0*Eeta[  8];

  tvec0 = terms[  5] *        1.0;

  G200 = tvec0*Eeta[  0];
  G201 = tvec0*Eeta[  1];
  G202 = tvec0*Eeta[  2];
  G210 = tvec0*Eeta[  3];
  t2 = t2 + (G120 - G210);
  G211 = tvec0*Eeta[  4];
  G212 = tvec0*Eeta[  5];
  G220 = tvec0*Eeta[  6];
  G221 = tvec0*Eeta[  7];
  t2 = t2 + (G201 - G021);
  G222 = tvec0*Eeta[  8];

  sigma[  0][  1] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  4] *        1.0;

  G000 = tvec0*Eeta[  0];
  G001 = tvec0*Eeta[  1];
  G002 = tvec0*Eeta[  2];
  G010 = tvec0*Eeta[  3];
  G011 = tvec0*Eeta[  4];
  G012 = tvec0*Eeta[  5];
  G020 = tvec0*Eeta[  6];
  G021 = tvec0*Eeta[  7];
  G022 = tvec0*Eeta[  8];

  tvec0 = terms[  5] *        1.0;

  G100 = tvec0*Eeta[  0];
  G101 = tvec0*Eeta[  1];
  G102 = tvec0*Eeta[  2];
  G110 = tvec0*Eeta[  3];
  G111 = tvec0*Eeta[  4];
  G112 = tvec0*Eeta[  5];
  t2 = t2 + (G012 - G102);
  G120 = tvec0*Eeta[  6];
  G121 = tvec0*Eeta[  7];
  G122 = tvec0*Eeta[  8];

  tvec0 = terms[  7] *        1.0;

  G200 = tvec0*Eeta[  0];
  G201 = tvec0*Eeta[  1];
  G202 = tvec0*Eeta[  2];
  G210 = tvec0*Eeta[  3];
  t2 = t2 + (G120 - G210);
  G211 = tvec0*Eeta[  4];
  G212 = tvec0*Eeta[  5];
  G220 = tvec0*Eeta[  6];
  G221 = tvec0*Eeta[  7];
  t2 = t2 + (G201 - G021);
  G222 = tvec0*Eeta[  8];

  sigma[  0][  2] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  2] *        1.0;

  G000 = tvec0*Eeta[  0];
  G001 = tvec0*Eeta[  1];
  G002 = tvec0*Eeta[  2];
  G010 = tvec0*Eeta[  3];
  G011 = tvec0*Eeta[  4];
  G012 = tvec0*Eeta[  5];
  G020 = tvec0*Eeta[  6];
  G021 = tvec0*Eeta[  7];
  G022 = tvec0*Eeta[  8];

  H000 = H000 + G000;
  H001 = H001 + G001;
  H002 = H002 + G002;
  H010 = H010 + G010;
  H011 = H011 + G011;
  H012 = H012 + G012;
  H020 = H020 + G020;
  H021 = H021 + G021;
  H022 = H022 + G022;

  tvec0 = terms[  3] *        1.0;

  G100 = tvec0*Eeta[  0];
  G101 = tvec0*Eeta[  1];
  G102 = tvec0*Eeta[  2];
  G110 = tvec0*Eeta[  3];
  G111 = tvec0*Eeta[  4];
  G112 = tvec0*Eeta[  5];
  t2 = t2 + (G012 - G102);
  G120 = tvec0*Eeta[  6];
  G121 = tvec0*Eeta[  7];
  G122 = tvec0*Eeta[  8];

  H100 = H100 + G100;
  H101 = H101 + G101;
  H102 = H102 + G102;
  H110 = H110 + G110;
  H111 = H111 + G111;
  H112 = H112 + G112;
  H120 = H120 + G120;
  H121 = H121 + G121;
  H122 = H122 + G122;

  tvec0 = terms[  6] *        1.0;

  G200 = tvec0*Eeta[  0];
  G201 = tvec0*Eeta[  1];
  G202 = tvec0*Eeta[  2];
  G210 = tvec0*Eeta[  3];
  t2 = t2 + (G120 - G210);
  G211 = tvec0*Eeta[  4];
  G212 = tvec0*Eeta[  5];
  G220 = tvec0*Eeta[  6];
  G221 = tvec0*Eeta[  7];
  t2 = t2 + (G201 - G021);
  G222 = tvec0*Eeta[  8];

  H200 = H200 + G200;
  H201 = H201 + G201;
  H202 = H202 + G202;
  H210 = H210 + G210;
  H211 = H211 + G211;
  H212 = H212 + G212;
  H220 = H220 + G220;
  H221 = H221 + G221;
  H222 = H222 + G222;

  sigma[  1][  1] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  5] *        1.0;

  G000 = tvec0*Eeta[  0];
  G001 = tvec0*Eeta[  1];
  G002 = tvec0*Eeta[  2];
  G010 = tvec0*Eeta[  3];
  G011 = tvec0*Eeta[  4];
  G012 = tvec0*Eeta[  5];
  G020 = tvec0*Eeta[  6];
  G021 = tvec0*Eeta[  7];
  G022 = tvec0*Eeta[  8];

  tvec0 = terms[  6] *        1.0;

  G100 = tvec0*Eeta[  0];
  G101 = tvec0*Eeta[  1];
  G102 = tvec0*Eeta[  2];
  G110 = tvec0*Eeta[  3];
  G111 = tvec0*Eeta[  4];
  G112 = tvec0*Eeta[  5];
  t2 = t2 + (G012 - G102);
  G120 = tvec0*Eeta[  6];
  G121 = tvec0*Eeta[  7];
  G122 = tvec0*Eeta[  8];

  tvec0 = terms[  8] *        1.0;

  G200 = tvec0*Eeta[  0];
  G201 = tvec0*Eeta[  1];
  G202 = tvec0*Eeta[  2];
  G210 = tvec0*Eeta[  3];
  t2 = t2 + (G120 - G210);
  G211 = tvec0*Eeta[  4];
  G212 = tvec0*Eeta[  5];
  G220 = tvec0*Eeta[  6];
  G221 = tvec0*Eeta[  7];
  t2 = t2 + (G201 - G021);
  G222 = tvec0*Eeta[  8];

  sigma[  1][  2] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  7] *        1.0;

  G000 = tvec0*Eeta[  0];
  G001 = tvec0*Eeta[  1];
  G002 = tvec0*Eeta[  2];
  G010 = tvec0*Eeta[  3];
  G011 = tvec0*Eeta[  4];
  G012 = tvec0*Eeta[  5];
  G020 = tvec0*Eeta[  6];
  G021 = tvec0*Eeta[  7];
  G022 = tvec0*Eeta[  8];

  H000 = H000 + G000;
  H001 = H001 + G001;
  H002 = H002 + G002;
  H010 = H010 + G010;
  H011 = H011 + G011;
  H012 = H012 + G012;
  H020 = H020 + G020;
  H021 = H021 + G021;
  H022 = H022 + G022;

  tvec0 = terms[  8] *        1.0;

  G100 = tvec0*Eeta[  0];
  G101 = tvec0*Eeta[  1];
  G102 = tvec0*Eeta[  2];
  G110 = tvec0*Eeta[  3];
  G111 = tvec0*Eeta[  4];
  G112 = tvec0*Eeta[  5];
  t2 = t2 + (G012 - G102);
  G120 = tvec0*Eeta[  6];
  G121 = tvec0*Eeta[  7];
  G122 = tvec0*Eeta[  8];

  H100 = H100 + G100;
  H101 = H101 + G101;
  H102 = H102 + G102;
  H110 = H110 + G110;
  H111 = H111 + G111;
  H112 = H112 + G112;
  H120 = H120 + G120;
  H121 = H121 + G121;
  H122 = H122 + G122;

  tvec0 = terms[  9] *        1.0;

  G200 = tvec0*Eeta[  0];
  G201 = tvec0*Eeta[  1];
  G202 = tvec0*Eeta[  2];
  G210 = tvec0*Eeta[  3];
  t2 = t2 + (G120 - G210);
  G211 = tvec0*Eeta[  4];
  G212 = tvec0*Eeta[  5];
  G220 = tvec0*Eeta[  6];
  G221 = tvec0*Eeta[  7];
  t2 = t2 + (G201 - G021);
  G222 = tvec0*Eeta[  8];

  H200 = H200 + G200;
  H201 = H201 + G201;
  H202 = H202 + G202;
  H210 = H210 + G210;
  H211 = H211 + G211;
  H212 = H212 + G212;
  H220 = H220 + G220;
  H221 = H221 + G221;
  H222 = H222 + G222;

  sigma[  2][  2] = two1nu * t2;

  t11 = H120 - H210;
  t12 = H120 - H210;
  sigma[  0][  0] = sigma[  0][  0] + (t11 + t12);
  t11 = H200 - H020;
  t12 = H121 - H211;
  sigma[  0][  1] = sigma[  0][  1] + (t11 + t12);
  t11 = H010 - H100;
  t12 = H122 - H212;
  sigma[  0][  2] = sigma[  0][  2] + (t11 + t12);
  t11 = H201 - H021;
  t12 = H201 - H021;
  sigma[  1][  1] = sigma[  1][  1] + (t11 + t12);
  t11 = H011 - H101;
  t12 = H202 - H022;
  sigma[  1][  2] = sigma[  1][  2] + (t11 + t12);
  t11 = H012 - H102;
  t12 = H012 - H102;
  sigma[  2][  2] = sigma[  2][  2] + (t11 + t12);

  t3 = 0.0;
  t3 = t3 + H120 - H210;
  t3 = t3 + H201 - H021;
  t3 = t3 + H012 - H102;
  t3 = t3 * two1nu;
sigma[  0][  0] = sigma[  0][  0] - t3;
sigma[  1][  1] = sigma[  1][  1] - t3;
sigma[  2][  2] = sigma[  2][  2] - t3;

  t =     1.00000000000000000000e+00 * mu8pi;
  sigma[  0][  0] = t * sigma[  0][  0];
  sigma[  0][  1] = t * sigma[  0][  1];
  sigma[  0][  2] = t * sigma[  0][  2];
  sigma[  1][  1] = t * sigma[  1][  1];
  sigma[  1][  2] = t * sigma[  1][  2];
  sigma[  2][  2] = t * sigma[  2][  2];
  sigma[1][0] = sigma[0][1];
  sigma[2][0] = sigma[0][2];
  sigma[2][1] = sigma[1][2];
}



void FMSigma2core1(real8 terms[],real8 mu8pi, real8 two1nu,real8 Eeta[],matrix sigma) {
real8 H000,H001,H002,H010,H011,H012,H020,H021,H022;
real8 H100,H101,H102,H110,H111,H112,H120,H121,H122;
real8 H200,H201,H202,H210,H211,H212,H220,H221,H222;
real8 G000,G001,G002,G010,G011,G012,G020,G021,G022;
real8 G100,G101,G102,G110,G111,G112,G120,G121,G122;
real8 G200,G201,G202,G210,G211,G212,G220,G221,G222;
real8 t,t11,t12,t2,t3;
real8 tvec0,tvec1,tvec2;

  t2 = 0.0;
  tvec0 = terms[  0] *        1.0;
  tvec1 = terms[  1] *        1.0;
  tvec2 = terms[  5] *        1.0;

  G000 = tvec0*Eeta[  9];
  G001 = tvec0*Eeta[ 12];
  G002 = tvec0*Eeta[ 15];
  G010 = tvec0*Eeta[ 18];
  G011 = tvec0*Eeta[ 21];
  G012 = tvec0*Eeta[ 24];
  G020 = tvec0*Eeta[ 27];
  G021 = tvec0*Eeta[ 30];
  G022 = tvec0*Eeta[ 33];

  G000 = G000 + tvec1*Eeta[ 10];
  G001 = G001 + tvec1*Eeta[ 13];
  G002 = G002 + tvec1*Eeta[ 16];
  G010 = G010 + tvec1*Eeta[ 19];
  G011 = G011 + tvec1*Eeta[ 22];
  G012 = G012 + tvec1*Eeta[ 25];
  G020 = G020 + tvec1*Eeta[ 28];
  G021 = G021 + tvec1*Eeta[ 31];
  G022 = G022 + tvec1*Eeta[ 34];

  G000 = G000 + tvec2*Eeta[ 11];
  G001 = G001 + tvec2*Eeta[ 14];
  G002 = G002 + tvec2*Eeta[ 17];
  G010 = G010 + tvec2*Eeta[ 20];
  G011 = G011 + tvec2*Eeta[ 23];
  G012 = G012 + tvec2*Eeta[ 26];
  G020 = G020 + tvec2*Eeta[ 29];
  G021 = G021 + tvec2*Eeta[ 32];
  G022 = G022 + tvec2*Eeta[ 35];

  H000 = G000;
  H001 = G001;
  H002 = G002;
  H010 = G010;
  H011 = G011;
  H012 = G012;
  H020 = G020;
  H021 = G021;
  H022 = G022;

  tvec0 = terms[  1] *        1.0;
  tvec1 = terms[  2] *        1.0;
  tvec2 = terms[  6] *        1.0;

  G100 = tvec0*Eeta[  9];
  G101 = tvec0*Eeta[ 12];
  G102 = tvec0*Eeta[ 15];
  G110 = tvec0*Eeta[ 18];
  G111 = tvec0*Eeta[ 21];
  G112 = tvec0*Eeta[ 24];
  G120 = tvec0*Eeta[ 27];
  G121 = tvec0*Eeta[ 30];
  G122 = tvec0*Eeta[ 33];

  G100 = G100 + tvec1*Eeta[ 10];
  G101 = G101 + tvec1*Eeta[ 13];
  G102 = G102 + tvec1*Eeta[ 16];
  G110 = G110 + tvec1*Eeta[ 19];
  G111 = G111 + tvec1*Eeta[ 22];
  G112 = G112 + tvec1*Eeta[ 25];
  G120 = G120 + tvec1*Eeta[ 28];
  G121 = G121 + tvec1*Eeta[ 31];
  G122 = G122 + tvec1*Eeta[ 34];

  G100 = G100 + tvec2*Eeta[ 11];
  G101 = G101 + tvec2*Eeta[ 14];
  G102 = G102 + tvec2*Eeta[ 17];
  G110 = G110 + tvec2*Eeta[ 20];
  G111 = G111 + tvec2*Eeta[ 23];
  G112 = G112 + tvec2*Eeta[ 26];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec2*Eeta[ 29];
  G121 = G121 + tvec2*Eeta[ 32];
  G122 = G122 + tvec2*Eeta[ 35];

  H100 = G100;
  H101 = G101;
  H102 = G102;
  H110 = G110;
  H111 = G111;
  H112 = G112;
  H120 = G120;
  H121 = G121;
  H122 = G122;

  tvec0 = terms[  5] *        1.0;
  tvec1 = terms[  6] *        1.0;
  tvec2 = terms[  9] *        1.0;

  G200 = tvec0*Eeta[  9];
  G201 = tvec0*Eeta[ 12];
  G202 = tvec0*Eeta[ 15];
  G210 = tvec0*Eeta[ 18];
  G211 = tvec0*Eeta[ 21];
  G212 = tvec0*Eeta[ 24];
  G220 = tvec0*Eeta[ 27];
  G221 = tvec0*Eeta[ 30];
  G222 = tvec0*Eeta[ 33];

  G200 = G200 + tvec1*Eeta[ 10];
  G201 = G201 + tvec1*Eeta[ 13];
  G202 = G202 + tvec1*Eeta[ 16];
  G210 = G210 + tvec1*Eeta[ 19];
  G211 = G211 + tvec1*Eeta[ 22];
  G212 = G212 + tvec1*Eeta[ 25];
  G220 = G220 + tvec1*Eeta[ 28];
  G221 = G221 + tvec1*Eeta[ 31];
  G222 = G222 + tvec1*Eeta[ 34];

  G200 = G200 + tvec2*Eeta[ 11];
  G201 = G201 + tvec2*Eeta[ 14];
  G202 = G202 + tvec2*Eeta[ 17];
  G210 = G210 + tvec2*Eeta[ 20];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec2*Eeta[ 23];
  G212 = G212 + tvec2*Eeta[ 26];
  G220 = G220 + tvec2*Eeta[ 29];
  G221 = G221 + tvec2*Eeta[ 32];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec2*Eeta[ 35];

  H200 = G200;
  H201 = G201;
  H202 = G202;
  H210 = G210;
  H211 = G211;
  H212 = G212;
  H220 = G220;
  H221 = G221;
  H222 = G222;

  sigma[  0][  0] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  1] *        1.0;
  tvec1 = terms[  2] *        1.0;
  tvec2 = terms[  6] *        1.0;

  G000 = tvec0*Eeta[  9];
  G001 = tvec0*Eeta[ 12];
  G002 = tvec0*Eeta[ 15];
  G010 = tvec0*Eeta[ 18];
  G011 = tvec0*Eeta[ 21];
  G012 = tvec0*Eeta[ 24];
  G020 = tvec0*Eeta[ 27];
  G021 = tvec0*Eeta[ 30];
  G022 = tvec0*Eeta[ 33];

  G000 = G000 + tvec1*Eeta[ 10];
  G001 = G001 + tvec1*Eeta[ 13];
  G002 = G002 + tvec1*Eeta[ 16];
  G010 = G010 + tvec1*Eeta[ 19];
  G011 = G011 + tvec1*Eeta[ 22];
  G012 = G012 + tvec1*Eeta[ 25];
  G020 = G020 + tvec1*Eeta[ 28];
  G021 = G021 + tvec1*Eeta[ 31];
  G022 = G022 + tvec1*Eeta[ 34];

  G000 = G000 + tvec2*Eeta[ 11];
  G001 = G001 + tvec2*Eeta[ 14];
  G002 = G002 + tvec2*Eeta[ 17];
  G010 = G010 + tvec2*Eeta[ 20];
  G011 = G011 + tvec2*Eeta[ 23];
  G012 = G012 + tvec2*Eeta[ 26];
  G020 = G020 + tvec2*Eeta[ 29];
  G021 = G021 + tvec2*Eeta[ 32];
  G022 = G022 + tvec2*Eeta[ 35];

  tvec0 = terms[  2] *        1.0;
  tvec1 = terms[  3] *        1.0;
  tvec2 = terms[  7] *        1.0;

  G100 = tvec0*Eeta[  9];
  G101 = tvec0*Eeta[ 12];
  G102 = tvec0*Eeta[ 15];
  G110 = tvec0*Eeta[ 18];
  G111 = tvec0*Eeta[ 21];
  G112 = tvec0*Eeta[ 24];
  G120 = tvec0*Eeta[ 27];
  G121 = tvec0*Eeta[ 30];
  G122 = tvec0*Eeta[ 33];

  G100 = G100 + tvec1*Eeta[ 10];
  G101 = G101 + tvec1*Eeta[ 13];
  G102 = G102 + tvec1*Eeta[ 16];
  G110 = G110 + tvec1*Eeta[ 19];
  G111 = G111 + tvec1*Eeta[ 22];
  G112 = G112 + tvec1*Eeta[ 25];
  G120 = G120 + tvec1*Eeta[ 28];
  G121 = G121 + tvec1*Eeta[ 31];
  G122 = G122 + tvec1*Eeta[ 34];

  G100 = G100 + tvec2*Eeta[ 11];
  G101 = G101 + tvec2*Eeta[ 14];
  G102 = G102 + tvec2*Eeta[ 17];
  G110 = G110 + tvec2*Eeta[ 20];
  G111 = G111 + tvec2*Eeta[ 23];
  G112 = G112 + tvec2*Eeta[ 26];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec2*Eeta[ 29];
  G121 = G121 + tvec2*Eeta[ 32];
  G122 = G122 + tvec2*Eeta[ 35];

  tvec0 = terms[  6] *        1.0;
  tvec1 = terms[  7] *        1.0;
  tvec2 = terms[ 10] *        1.0;

  G200 = tvec0*Eeta[  9];
  G201 = tvec0*Eeta[ 12];
  G202 = tvec0*Eeta[ 15];
  G210 = tvec0*Eeta[ 18];
  G211 = tvec0*Eeta[ 21];
  G212 = tvec0*Eeta[ 24];
  G220 = tvec0*Eeta[ 27];
  G221 = tvec0*Eeta[ 30];
  G222 = tvec0*Eeta[ 33];

  G200 = G200 + tvec1*Eeta[ 10];
  G201 = G201 + tvec1*Eeta[ 13];
  G202 = G202 + tvec1*Eeta[ 16];
  G210 = G210 + tvec1*Eeta[ 19];
  G211 = G211 + tvec1*Eeta[ 22];
  G212 = G212 + tvec1*Eeta[ 25];
  G220 = G220 + tvec1*Eeta[ 28];
  G221 = G221 + tvec1*Eeta[ 31];
  G222 = G222 + tvec1*Eeta[ 34];

  G200 = G200 + tvec2*Eeta[ 11];
  G201 = G201 + tvec2*Eeta[ 14];
  G202 = G202 + tvec2*Eeta[ 17];
  G210 = G210 + tvec2*Eeta[ 20];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec2*Eeta[ 23];
  G212 = G212 + tvec2*Eeta[ 26];
  G220 = G220 + tvec2*Eeta[ 29];
  G221 = G221 + tvec2*Eeta[ 32];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec2*Eeta[ 35];

  sigma[  0][  1] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  5] *        1.0;
  tvec1 = terms[  6] *        1.0;
  tvec2 = terms[  9] *        1.0;

  G000 = tvec0*Eeta[  9];
  G001 = tvec0*Eeta[ 12];
  G002 = tvec0*Eeta[ 15];
  G010 = tvec0*Eeta[ 18];
  G011 = tvec0*Eeta[ 21];
  G012 = tvec0*Eeta[ 24];
  G020 = tvec0*Eeta[ 27];
  G021 = tvec0*Eeta[ 30];
  G022 = tvec0*Eeta[ 33];

  G000 = G000 + tvec1*Eeta[ 10];
  G001 = G001 + tvec1*Eeta[ 13];
  G002 = G002 + tvec1*Eeta[ 16];
  G010 = G010 + tvec1*Eeta[ 19];
  G011 = G011 + tvec1*Eeta[ 22];
  G012 = G012 + tvec1*Eeta[ 25];
  G020 = G020 + tvec1*Eeta[ 28];
  G021 = G021 + tvec1*Eeta[ 31];
  G022 = G022 + tvec1*Eeta[ 34];

  G000 = G000 + tvec2*Eeta[ 11];
  G001 = G001 + tvec2*Eeta[ 14];
  G002 = G002 + tvec2*Eeta[ 17];
  G010 = G010 + tvec2*Eeta[ 20];
  G011 = G011 + tvec2*Eeta[ 23];
  G012 = G012 + tvec2*Eeta[ 26];
  G020 = G020 + tvec2*Eeta[ 29];
  G021 = G021 + tvec2*Eeta[ 32];
  G022 = G022 + tvec2*Eeta[ 35];

  tvec0 = terms[  6] *        1.0;
  tvec1 = terms[  7] *        1.0;
  tvec2 = terms[ 10] *        1.0;

  G100 = tvec0*Eeta[  9];
  G101 = tvec0*Eeta[ 12];
  G102 = tvec0*Eeta[ 15];
  G110 = tvec0*Eeta[ 18];
  G111 = tvec0*Eeta[ 21];
  G112 = tvec0*Eeta[ 24];
  G120 = tvec0*Eeta[ 27];
  G121 = tvec0*Eeta[ 30];
  G122 = tvec0*Eeta[ 33];

  G100 = G100 + tvec1*Eeta[ 10];
  G101 = G101 + tvec1*Eeta[ 13];
  G102 = G102 + tvec1*Eeta[ 16];
  G110 = G110 + tvec1*Eeta[ 19];
  G111 = G111 + tvec1*Eeta[ 22];
  G112 = G112 + tvec1*Eeta[ 25];
  G120 = G120 + tvec1*Eeta[ 28];
  G121 = G121 + tvec1*Eeta[ 31];
  G122 = G122 + tvec1*Eeta[ 34];

  G100 = G100 + tvec2*Eeta[ 11];
  G101 = G101 + tvec2*Eeta[ 14];
  G102 = G102 + tvec2*Eeta[ 17];
  G110 = G110 + tvec2*Eeta[ 20];
  G111 = G111 + tvec2*Eeta[ 23];
  G112 = G112 + tvec2*Eeta[ 26];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec2*Eeta[ 29];
  G121 = G121 + tvec2*Eeta[ 32];
  G122 = G122 + tvec2*Eeta[ 35];

  tvec0 = terms[  9] *        1.0;
  tvec1 = terms[ 10] *        1.0;
  tvec2 = terms[ 12] *        1.0;

  G200 = tvec0*Eeta[  9];
  G201 = tvec0*Eeta[ 12];
  G202 = tvec0*Eeta[ 15];
  G210 = tvec0*Eeta[ 18];
  G211 = tvec0*Eeta[ 21];
  G212 = tvec0*Eeta[ 24];
  G220 = tvec0*Eeta[ 27];
  G221 = tvec0*Eeta[ 30];
  G222 = tvec0*Eeta[ 33];

  G200 = G200 + tvec1*Eeta[ 10];
  G201 = G201 + tvec1*Eeta[ 13];
  G202 = G202 + tvec1*Eeta[ 16];
  G210 = G210 + tvec1*Eeta[ 19];
  G211 = G211 + tvec1*Eeta[ 22];
  G212 = G212 + tvec1*Eeta[ 25];
  G220 = G220 + tvec1*Eeta[ 28];
  G221 = G221 + tvec1*Eeta[ 31];
  G222 = G222 + tvec1*Eeta[ 34];

  G200 = G200 + tvec2*Eeta[ 11];
  G201 = G201 + tvec2*Eeta[ 14];
  G202 = G202 + tvec2*Eeta[ 17];
  G210 = G210 + tvec2*Eeta[ 20];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec2*Eeta[ 23];
  G212 = G212 + tvec2*Eeta[ 26];
  G220 = G220 + tvec2*Eeta[ 29];
  G221 = G221 + tvec2*Eeta[ 32];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec2*Eeta[ 35];

  sigma[  0][  2] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  2] *        1.0;
  tvec1 = terms[  3] *        1.0;
  tvec2 = terms[  7] *        1.0;

  G000 = tvec0*Eeta[  9];
  G001 = tvec0*Eeta[ 12];
  G002 = tvec0*Eeta[ 15];
  G010 = tvec0*Eeta[ 18];
  G011 = tvec0*Eeta[ 21];
  G012 = tvec0*Eeta[ 24];
  G020 = tvec0*Eeta[ 27];
  G021 = tvec0*Eeta[ 30];
  G022 = tvec0*Eeta[ 33];

  G000 = G000 + tvec1*Eeta[ 10];
  G001 = G001 + tvec1*Eeta[ 13];
  G002 = G002 + tvec1*Eeta[ 16];
  G010 = G010 + tvec1*Eeta[ 19];
  G011 = G011 + tvec1*Eeta[ 22];
  G012 = G012 + tvec1*Eeta[ 25];
  G020 = G020 + tvec1*Eeta[ 28];
  G021 = G021 + tvec1*Eeta[ 31];
  G022 = G022 + tvec1*Eeta[ 34];

  G000 = G000 + tvec2*Eeta[ 11];
  G001 = G001 + tvec2*Eeta[ 14];
  G002 = G002 + tvec2*Eeta[ 17];
  G010 = G010 + tvec2*Eeta[ 20];
  G011 = G011 + tvec2*Eeta[ 23];
  G012 = G012 + tvec2*Eeta[ 26];
  G020 = G020 + tvec2*Eeta[ 29];
  G021 = G021 + tvec2*Eeta[ 32];
  G022 = G022 + tvec2*Eeta[ 35];

  H000 = H000 + G000;
  H001 = H001 + G001;
  H002 = H002 + G002;
  H010 = H010 + G010;
  H011 = H011 + G011;
  H012 = H012 + G012;
  H020 = H020 + G020;
  H021 = H021 + G021;
  H022 = H022 + G022;

  tvec0 = terms[  3] *        1.0;
  tvec1 = terms[  4] *        1.0;
  tvec2 = terms[  8] *        1.0;

  G100 = tvec0*Eeta[  9];
  G101 = tvec0*Eeta[ 12];
  G102 = tvec0*Eeta[ 15];
  G110 = tvec0*Eeta[ 18];
  G111 = tvec0*Eeta[ 21];
  G112 = tvec0*Eeta[ 24];
  G120 = tvec0*Eeta[ 27];
  G121 = tvec0*Eeta[ 30];
  G122 = tvec0*Eeta[ 33];

  G100 = G100 + tvec1*Eeta[ 10];
  G101 = G101 + tvec1*Eeta[ 13];
  G102 = G102 + tvec1*Eeta[ 16];
  G110 = G110 + tvec1*Eeta[ 19];
  G111 = G111 + tvec1*Eeta[ 22];
  G112 = G112 + tvec1*Eeta[ 25];
  G120 = G120 + tvec1*Eeta[ 28];
  G121 = G121 + tvec1*Eeta[ 31];
  G122 = G122 + tvec1*Eeta[ 34];

  G100 = G100 + tvec2*Eeta[ 11];
  G101 = G101 + tvec2*Eeta[ 14];
  G102 = G102 + tvec2*Eeta[ 17];
  G110 = G110 + tvec2*Eeta[ 20];
  G111 = G111 + tvec2*Eeta[ 23];
  G112 = G112 + tvec2*Eeta[ 26];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec2*Eeta[ 29];
  G121 = G121 + tvec2*Eeta[ 32];
  G122 = G122 + tvec2*Eeta[ 35];

  H100 = H100 + G100;
  H101 = H101 + G101;
  H102 = H102 + G102;
  H110 = H110 + G110;
  H111 = H111 + G111;
  H112 = H112 + G112;
  H120 = H120 + G120;
  H121 = H121 + G121;
  H122 = H122 + G122;

  tvec0 = terms[  7] *        1.0;
  tvec1 = terms[  8] *        1.0;
  tvec2 = terms[ 11] *        1.0;

  G200 = tvec0*Eeta[  9];
  G201 = tvec0*Eeta[ 12];
  G202 = tvec0*Eeta[ 15];
  G210 = tvec0*Eeta[ 18];
  G211 = tvec0*Eeta[ 21];
  G212 = tvec0*Eeta[ 24];
  G220 = tvec0*Eeta[ 27];
  G221 = tvec0*Eeta[ 30];
  G222 = tvec0*Eeta[ 33];

  G200 = G200 + tvec1*Eeta[ 10];
  G201 = G201 + tvec1*Eeta[ 13];
  G202 = G202 + tvec1*Eeta[ 16];
  G210 = G210 + tvec1*Eeta[ 19];
  G211 = G211 + tvec1*Eeta[ 22];
  G212 = G212 + tvec1*Eeta[ 25];
  G220 = G220 + tvec1*Eeta[ 28];
  G221 = G221 + tvec1*Eeta[ 31];
  G222 = G222 + tvec1*Eeta[ 34];

  G200 = G200 + tvec2*Eeta[ 11];
  G201 = G201 + tvec2*Eeta[ 14];
  G202 = G202 + tvec2*Eeta[ 17];
  G210 = G210 + tvec2*Eeta[ 20];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec2*Eeta[ 23];
  G212 = G212 + tvec2*Eeta[ 26];
  G220 = G220 + tvec2*Eeta[ 29];
  G221 = G221 + tvec2*Eeta[ 32];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec2*Eeta[ 35];

  H200 = H200 + G200;
  H201 = H201 + G201;
  H202 = H202 + G202;
  H210 = H210 + G210;
  H211 = H211 + G211;
  H212 = H212 + G212;
  H220 = H220 + G220;
  H221 = H221 + G221;
  H222 = H222 + G222;

  sigma[  1][  1] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  6] *        1.0;
  tvec1 = terms[  7] *        1.0;
  tvec2 = terms[ 10] *        1.0;

  G000 = tvec0*Eeta[  9];
  G001 = tvec0*Eeta[ 12];
  G002 = tvec0*Eeta[ 15];
  G010 = tvec0*Eeta[ 18];
  G011 = tvec0*Eeta[ 21];
  G012 = tvec0*Eeta[ 24];
  G020 = tvec0*Eeta[ 27];
  G021 = tvec0*Eeta[ 30];
  G022 = tvec0*Eeta[ 33];

  G000 = G000 + tvec1*Eeta[ 10];
  G001 = G001 + tvec1*Eeta[ 13];
  G002 = G002 + tvec1*Eeta[ 16];
  G010 = G010 + tvec1*Eeta[ 19];
  G011 = G011 + tvec1*Eeta[ 22];
  G012 = G012 + tvec1*Eeta[ 25];
  G020 = G020 + tvec1*Eeta[ 28];
  G021 = G021 + tvec1*Eeta[ 31];
  G022 = G022 + tvec1*Eeta[ 34];

  G000 = G000 + tvec2*Eeta[ 11];
  G001 = G001 + tvec2*Eeta[ 14];
  G002 = G002 + tvec2*Eeta[ 17];
  G010 = G010 + tvec2*Eeta[ 20];
  G011 = G011 + tvec2*Eeta[ 23];
  G012 = G012 + tvec2*Eeta[ 26];
  G020 = G020 + tvec2*Eeta[ 29];
  G021 = G021 + tvec2*Eeta[ 32];
  G022 = G022 + tvec2*Eeta[ 35];

  tvec0 = terms[  7] *        1.0;
  tvec1 = terms[  8] *        1.0;
  tvec2 = terms[ 11] *        1.0;

  G100 = tvec0*Eeta[  9];
  G101 = tvec0*Eeta[ 12];
  G102 = tvec0*Eeta[ 15];
  G110 = tvec0*Eeta[ 18];
  G111 = tvec0*Eeta[ 21];
  G112 = tvec0*Eeta[ 24];
  G120 = tvec0*Eeta[ 27];
  G121 = tvec0*Eeta[ 30];
  G122 = tvec0*Eeta[ 33];

  G100 = G100 + tvec1*Eeta[ 10];
  G101 = G101 + tvec1*Eeta[ 13];
  G102 = G102 + tvec1*Eeta[ 16];
  G110 = G110 + tvec1*Eeta[ 19];
  G111 = G111 + tvec1*Eeta[ 22];
  G112 = G112 + tvec1*Eeta[ 25];
  G120 = G120 + tvec1*Eeta[ 28];
  G121 = G121 + tvec1*Eeta[ 31];
  G122 = G122 + tvec1*Eeta[ 34];

  G100 = G100 + tvec2*Eeta[ 11];
  G101 = G101 + tvec2*Eeta[ 14];
  G102 = G102 + tvec2*Eeta[ 17];
  G110 = G110 + tvec2*Eeta[ 20];
  G111 = G111 + tvec2*Eeta[ 23];
  G112 = G112 + tvec2*Eeta[ 26];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec2*Eeta[ 29];
  G121 = G121 + tvec2*Eeta[ 32];
  G122 = G122 + tvec2*Eeta[ 35];

  tvec0 = terms[ 10] *        1.0;
  tvec1 = terms[ 11] *        1.0;
  tvec2 = terms[ 13] *        1.0;

  G200 = tvec0*Eeta[  9];
  G201 = tvec0*Eeta[ 12];
  G202 = tvec0*Eeta[ 15];
  G210 = tvec0*Eeta[ 18];
  G211 = tvec0*Eeta[ 21];
  G212 = tvec0*Eeta[ 24];
  G220 = tvec0*Eeta[ 27];
  G221 = tvec0*Eeta[ 30];
  G222 = tvec0*Eeta[ 33];

  G200 = G200 + tvec1*Eeta[ 10];
  G201 = G201 + tvec1*Eeta[ 13];
  G202 = G202 + tvec1*Eeta[ 16];
  G210 = G210 + tvec1*Eeta[ 19];
  G211 = G211 + tvec1*Eeta[ 22];
  G212 = G212 + tvec1*Eeta[ 25];
  G220 = G220 + tvec1*Eeta[ 28];
  G221 = G221 + tvec1*Eeta[ 31];
  G222 = G222 + tvec1*Eeta[ 34];

  G200 = G200 + tvec2*Eeta[ 11];
  G201 = G201 + tvec2*Eeta[ 14];
  G202 = G202 + tvec2*Eeta[ 17];
  G210 = G210 + tvec2*Eeta[ 20];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec2*Eeta[ 23];
  G212 = G212 + tvec2*Eeta[ 26];
  G220 = G220 + tvec2*Eeta[ 29];
  G221 = G221 + tvec2*Eeta[ 32];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec2*Eeta[ 35];

  sigma[  1][  2] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  9] *        1.0;
  tvec1 = terms[ 10] *        1.0;
  tvec2 = terms[ 12] *        1.0;

  G000 = tvec0*Eeta[  9];
  G001 = tvec0*Eeta[ 12];
  G002 = tvec0*Eeta[ 15];
  G010 = tvec0*Eeta[ 18];
  G011 = tvec0*Eeta[ 21];
  G012 = tvec0*Eeta[ 24];
  G020 = tvec0*Eeta[ 27];
  G021 = tvec0*Eeta[ 30];
  G022 = tvec0*Eeta[ 33];

  G000 = G000 + tvec1*Eeta[ 10];
  G001 = G001 + tvec1*Eeta[ 13];
  G002 = G002 + tvec1*Eeta[ 16];
  G010 = G010 + tvec1*Eeta[ 19];
  G011 = G011 + tvec1*Eeta[ 22];
  G012 = G012 + tvec1*Eeta[ 25];
  G020 = G020 + tvec1*Eeta[ 28];
  G021 = G021 + tvec1*Eeta[ 31];
  G022 = G022 + tvec1*Eeta[ 34];

  G000 = G000 + tvec2*Eeta[ 11];
  G001 = G001 + tvec2*Eeta[ 14];
  G002 = G002 + tvec2*Eeta[ 17];
  G010 = G010 + tvec2*Eeta[ 20];
  G011 = G011 + tvec2*Eeta[ 23];
  G012 = G012 + tvec2*Eeta[ 26];
  G020 = G020 + tvec2*Eeta[ 29];
  G021 = G021 + tvec2*Eeta[ 32];
  G022 = G022 + tvec2*Eeta[ 35];

  H000 = H000 + G000;
  H001 = H001 + G001;
  H002 = H002 + G002;
  H010 = H010 + G010;
  H011 = H011 + G011;
  H012 = H012 + G012;
  H020 = H020 + G020;
  H021 = H021 + G021;
  H022 = H022 + G022;

  tvec0 = terms[ 10] *        1.0;
  tvec1 = terms[ 11] *        1.0;
  tvec2 = terms[ 13] *        1.0;

  G100 = tvec0*Eeta[  9];
  G101 = tvec0*Eeta[ 12];
  G102 = tvec0*Eeta[ 15];
  G110 = tvec0*Eeta[ 18];
  G111 = tvec0*Eeta[ 21];
  G112 = tvec0*Eeta[ 24];
  G120 = tvec0*Eeta[ 27];
  G121 = tvec0*Eeta[ 30];
  G122 = tvec0*Eeta[ 33];

  G100 = G100 + tvec1*Eeta[ 10];
  G101 = G101 + tvec1*Eeta[ 13];
  G102 = G102 + tvec1*Eeta[ 16];
  G110 = G110 + tvec1*Eeta[ 19];
  G111 = G111 + tvec1*Eeta[ 22];
  G112 = G112 + tvec1*Eeta[ 25];
  G120 = G120 + tvec1*Eeta[ 28];
  G121 = G121 + tvec1*Eeta[ 31];
  G122 = G122 + tvec1*Eeta[ 34];

  G100 = G100 + tvec2*Eeta[ 11];
  G101 = G101 + tvec2*Eeta[ 14];
  G102 = G102 + tvec2*Eeta[ 17];
  G110 = G110 + tvec2*Eeta[ 20];
  G111 = G111 + tvec2*Eeta[ 23];
  G112 = G112 + tvec2*Eeta[ 26];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec2*Eeta[ 29];
  G121 = G121 + tvec2*Eeta[ 32];
  G122 = G122 + tvec2*Eeta[ 35];

  H100 = H100 + G100;
  H101 = H101 + G101;
  H102 = H102 + G102;
  H110 = H110 + G110;
  H111 = H111 + G111;
  H112 = H112 + G112;
  H120 = H120 + G120;
  H121 = H121 + G121;
  H122 = H122 + G122;

  tvec0 = terms[ 12] *        1.0;
  tvec1 = terms[ 13] *        1.0;
  tvec2 = terms[ 14] *        1.0;

  G200 = tvec0*Eeta[  9];
  G201 = tvec0*Eeta[ 12];
  G202 = tvec0*Eeta[ 15];
  G210 = tvec0*Eeta[ 18];
  G211 = tvec0*Eeta[ 21];
  G212 = tvec0*Eeta[ 24];
  G220 = tvec0*Eeta[ 27];
  G221 = tvec0*Eeta[ 30];
  G222 = tvec0*Eeta[ 33];

  G200 = G200 + tvec1*Eeta[ 10];
  G201 = G201 + tvec1*Eeta[ 13];
  G202 = G202 + tvec1*Eeta[ 16];
  G210 = G210 + tvec1*Eeta[ 19];
  G211 = G211 + tvec1*Eeta[ 22];
  G212 = G212 + tvec1*Eeta[ 25];
  G220 = G220 + tvec1*Eeta[ 28];
  G221 = G221 + tvec1*Eeta[ 31];
  G222 = G222 + tvec1*Eeta[ 34];

  G200 = G200 + tvec2*Eeta[ 11];
  G201 = G201 + tvec2*Eeta[ 14];
  G202 = G202 + tvec2*Eeta[ 17];
  G210 = G210 + tvec2*Eeta[ 20];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec2*Eeta[ 23];
  G212 = G212 + tvec2*Eeta[ 26];
  G220 = G220 + tvec2*Eeta[ 29];
  G221 = G221 + tvec2*Eeta[ 32];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec2*Eeta[ 35];

  H200 = H200 + G200;
  H201 = H201 + G201;
  H202 = H202 + G202;
  H210 = H210 + G210;
  H211 = H211 + G211;
  H212 = H212 + G212;
  H220 = H220 + G220;
  H221 = H221 + G221;
  H222 = H222 + G222;

  sigma[  2][  2] = two1nu * t2;

  t11 = H120 - H210;
  t12 = H120 - H210;
  sigma[  0][  0] = sigma[  0][  0] + (t11 + t12);
  t11 = H200 - H020;
  t12 = H121 - H211;
  sigma[  0][  1] = sigma[  0][  1] + (t11 + t12);
  t11 = H010 - H100;
  t12 = H122 - H212;
  sigma[  0][  2] = sigma[  0][  2] + (t11 + t12);
  t11 = H201 - H021;
  t12 = H201 - H021;
  sigma[  1][  1] = sigma[  1][  1] + (t11 + t12);
  t11 = H011 - H101;
  t12 = H202 - H022;
  sigma[  1][  2] = sigma[  1][  2] + (t11 + t12);
  t11 = H012 - H102;
  t12 = H012 - H102;
  sigma[  2][  2] = sigma[  2][  2] + (t11 + t12);

  t3 = 0.0;
  t3 = t3 + H120 - H210;
  t3 = t3 + H201 - H021;
  t3 = t3 + H012 - H102;
  t3 = t3 * two1nu;
sigma[  0][  0] = sigma[  0][  0] - t3;
sigma[  1][  1] = sigma[  1][  1] - t3;
sigma[  2][  2] = sigma[  2][  2] - t3;

  t =    -1.00000000000000000000e+00 * mu8pi;
  sigma[  0][  0] = t * sigma[  0][  0];
  sigma[  0][  1] = t * sigma[  0][  1];
  sigma[  0][  2] = t * sigma[  0][  2];
  sigma[  1][  1] = t * sigma[  1][  1];
  sigma[  1][  2] = t * sigma[  1][  2];
  sigma[  2][  2] = t * sigma[  2][  2];
  sigma[1][0] = sigma[0][1];
  sigma[2][0] = sigma[0][2];
  sigma[2][1] = sigma[1][2];
}



void FMSigma2core2(real8 terms[],real8 mu8pi, real8 two1nu,real8 Eeta[],matrix sigma) {
real8 H000,H001,H002,H010,H011,H012,H020,H021,H022;
real8 H100,H101,H102,H110,H111,H112,H120,H121,H122;
real8 H200,H201,H202,H210,H211,H212,H220,H221,H222;
real8 G000,G001,G002,G010,G011,G012,G020,G021,G022;
real8 G100,G101,G102,G110,G111,G112,G120,G121,G122;
real8 G200,G201,G202,G210,G211,G212,G220,G221,G222;
real8 t,t11,t12,t2,t3;
real8 tvec0,tvec1,tvec2,tvec3,tvec4,tvec5;

  t2 = 0.0;
  tvec0 = terms[  0] *        1.0;
  tvec1 = terms[  1] *        2.0;
  tvec2 = terms[  2] *        1.0;
  tvec3 = terms[  6] *        2.0;
  tvec4 = terms[  7] *        2.0;
  tvec5 = terms[ 11] *        1.0;

  G000 = tvec0*Eeta[ 36];
  G001 = tvec0*Eeta[ 42];
  G002 = tvec0*Eeta[ 48];
  G010 = tvec0*Eeta[ 54];
  G011 = tvec0*Eeta[ 60];
  G012 = tvec0*Eeta[ 66];
  G020 = tvec0*Eeta[ 72];
  G021 = tvec0*Eeta[ 78];
  G022 = tvec0*Eeta[ 84];

  G000 = G000 + tvec1*Eeta[ 37];
  G001 = G001 + tvec1*Eeta[ 43];
  G002 = G002 + tvec1*Eeta[ 49];
  G010 = G010 + tvec1*Eeta[ 55];
  G011 = G011 + tvec1*Eeta[ 61];
  G012 = G012 + tvec1*Eeta[ 67];
  G020 = G020 + tvec1*Eeta[ 73];
  G021 = G021 + tvec1*Eeta[ 79];
  G022 = G022 + tvec1*Eeta[ 85];

  G000 = G000 + tvec2*Eeta[ 38];
  G001 = G001 + tvec2*Eeta[ 44];
  G002 = G002 + tvec2*Eeta[ 50];
  G010 = G010 + tvec2*Eeta[ 56];
  G011 = G011 + tvec2*Eeta[ 62];
  G012 = G012 + tvec2*Eeta[ 68];
  G020 = G020 + tvec2*Eeta[ 74];
  G021 = G021 + tvec2*Eeta[ 80];
  G022 = G022 + tvec2*Eeta[ 86];

  G000 = G000 + tvec3*Eeta[ 39];
  G001 = G001 + tvec3*Eeta[ 45];
  G002 = G002 + tvec3*Eeta[ 51];
  G010 = G010 + tvec3*Eeta[ 57];
  G011 = G011 + tvec3*Eeta[ 63];
  G012 = G012 + tvec3*Eeta[ 69];
  G020 = G020 + tvec3*Eeta[ 75];
  G021 = G021 + tvec3*Eeta[ 81];
  G022 = G022 + tvec3*Eeta[ 87];

  G000 = G000 + tvec4*Eeta[ 40];
  G001 = G001 + tvec4*Eeta[ 46];
  G002 = G002 + tvec4*Eeta[ 52];
  G010 = G010 + tvec4*Eeta[ 58];
  G011 = G011 + tvec4*Eeta[ 64];
  G012 = G012 + tvec4*Eeta[ 70];
  G020 = G020 + tvec4*Eeta[ 76];
  G021 = G021 + tvec4*Eeta[ 82];
  G022 = G022 + tvec4*Eeta[ 88];

  G000 = G000 + tvec5*Eeta[ 41];
  G001 = G001 + tvec5*Eeta[ 47];
  G002 = G002 + tvec5*Eeta[ 53];
  G010 = G010 + tvec5*Eeta[ 59];
  G011 = G011 + tvec5*Eeta[ 65];
  G012 = G012 + tvec5*Eeta[ 71];
  G020 = G020 + tvec5*Eeta[ 77];
  G021 = G021 + tvec5*Eeta[ 83];
  G022 = G022 + tvec5*Eeta[ 89];

  H000 = G000;
  H001 = G001;
  H002 = G002;
  H010 = G010;
  H011 = G011;
  H012 = G012;
  H020 = G020;
  H021 = G021;
  H022 = G022;

  tvec0 = terms[  1] *        1.0;
  tvec1 = terms[  2] *        2.0;
  tvec2 = terms[  3] *        1.0;
  tvec3 = terms[  7] *        2.0;
  tvec4 = terms[  8] *        2.0;
  tvec5 = terms[ 12] *        1.0;

  G100 = tvec0*Eeta[ 36];
  G101 = tvec0*Eeta[ 42];
  G102 = tvec0*Eeta[ 48];
  G110 = tvec0*Eeta[ 54];
  G111 = tvec0*Eeta[ 60];
  G112 = tvec0*Eeta[ 66];
  G120 = tvec0*Eeta[ 72];
  G121 = tvec0*Eeta[ 78];
  G122 = tvec0*Eeta[ 84];

  G100 = G100 + tvec1*Eeta[ 37];
  G101 = G101 + tvec1*Eeta[ 43];
  G102 = G102 + tvec1*Eeta[ 49];
  G110 = G110 + tvec1*Eeta[ 55];
  G111 = G111 + tvec1*Eeta[ 61];
  G112 = G112 + tvec1*Eeta[ 67];
  G120 = G120 + tvec1*Eeta[ 73];
  G121 = G121 + tvec1*Eeta[ 79];
  G122 = G122 + tvec1*Eeta[ 85];

  G100 = G100 + tvec2*Eeta[ 38];
  G101 = G101 + tvec2*Eeta[ 44];
  G102 = G102 + tvec2*Eeta[ 50];
  G110 = G110 + tvec2*Eeta[ 56];
  G111 = G111 + tvec2*Eeta[ 62];
  G112 = G112 + tvec2*Eeta[ 68];
  G120 = G120 + tvec2*Eeta[ 74];
  G121 = G121 + tvec2*Eeta[ 80];
  G122 = G122 + tvec2*Eeta[ 86];

  G100 = G100 + tvec3*Eeta[ 39];
  G101 = G101 + tvec3*Eeta[ 45];
  G102 = G102 + tvec3*Eeta[ 51];
  G110 = G110 + tvec3*Eeta[ 57];
  G111 = G111 + tvec3*Eeta[ 63];
  G112 = G112 + tvec3*Eeta[ 69];
  G120 = G120 + tvec3*Eeta[ 75];
  G121 = G121 + tvec3*Eeta[ 81];
  G122 = G122 + tvec3*Eeta[ 87];

  G100 = G100 + tvec4*Eeta[ 40];
  G101 = G101 + tvec4*Eeta[ 46];
  G102 = G102 + tvec4*Eeta[ 52];
  G110 = G110 + tvec4*Eeta[ 58];
  G111 = G111 + tvec4*Eeta[ 64];
  G112 = G112 + tvec4*Eeta[ 70];
  G120 = G120 + tvec4*Eeta[ 76];
  G121 = G121 + tvec4*Eeta[ 82];
  G122 = G122 + tvec4*Eeta[ 88];

  G100 = G100 + tvec5*Eeta[ 41];
  G101 = G101 + tvec5*Eeta[ 47];
  G102 = G102 + tvec5*Eeta[ 53];
  G110 = G110 + tvec5*Eeta[ 59];
  G111 = G111 + tvec5*Eeta[ 65];
  G112 = G112 + tvec5*Eeta[ 71];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec5*Eeta[ 77];
  G121 = G121 + tvec5*Eeta[ 83];
  G122 = G122 + tvec5*Eeta[ 89];

  H100 = G100;
  H101 = G101;
  H102 = G102;
  H110 = G110;
  H111 = G111;
  H112 = G112;
  H120 = G120;
  H121 = G121;
  H122 = G122;

  tvec0 = terms[  6] *        1.0;
  tvec1 = terms[  7] *        2.0;
  tvec2 = terms[  8] *        1.0;
  tvec3 = terms[ 11] *        2.0;
  tvec4 = terms[ 12] *        2.0;
  tvec5 = terms[ 15] *        1.0;

  G200 = tvec0*Eeta[ 36];
  G201 = tvec0*Eeta[ 42];
  G202 = tvec0*Eeta[ 48];
  G210 = tvec0*Eeta[ 54];
  G211 = tvec0*Eeta[ 60];
  G212 = tvec0*Eeta[ 66];
  G220 = tvec0*Eeta[ 72];
  G221 = tvec0*Eeta[ 78];
  G222 = tvec0*Eeta[ 84];

  G200 = G200 + tvec1*Eeta[ 37];
  G201 = G201 + tvec1*Eeta[ 43];
  G202 = G202 + tvec1*Eeta[ 49];
  G210 = G210 + tvec1*Eeta[ 55];
  G211 = G211 + tvec1*Eeta[ 61];
  G212 = G212 + tvec1*Eeta[ 67];
  G220 = G220 + tvec1*Eeta[ 73];
  G221 = G221 + tvec1*Eeta[ 79];
  G222 = G222 + tvec1*Eeta[ 85];

  G200 = G200 + tvec2*Eeta[ 38];
  G201 = G201 + tvec2*Eeta[ 44];
  G202 = G202 + tvec2*Eeta[ 50];
  G210 = G210 + tvec2*Eeta[ 56];
  G211 = G211 + tvec2*Eeta[ 62];
  G212 = G212 + tvec2*Eeta[ 68];
  G220 = G220 + tvec2*Eeta[ 74];
  G221 = G221 + tvec2*Eeta[ 80];
  G222 = G222 + tvec2*Eeta[ 86];

  G200 = G200 + tvec3*Eeta[ 39];
  G201 = G201 + tvec3*Eeta[ 45];
  G202 = G202 + tvec3*Eeta[ 51];
  G210 = G210 + tvec3*Eeta[ 57];
  G211 = G211 + tvec3*Eeta[ 63];
  G212 = G212 + tvec3*Eeta[ 69];
  G220 = G220 + tvec3*Eeta[ 75];
  G221 = G221 + tvec3*Eeta[ 81];
  G222 = G222 + tvec3*Eeta[ 87];

  G200 = G200 + tvec4*Eeta[ 40];
  G201 = G201 + tvec4*Eeta[ 46];
  G202 = G202 + tvec4*Eeta[ 52];
  G210 = G210 + tvec4*Eeta[ 58];
  G211 = G211 + tvec4*Eeta[ 64];
  G212 = G212 + tvec4*Eeta[ 70];
  G220 = G220 + tvec4*Eeta[ 76];
  G221 = G221 + tvec4*Eeta[ 82];
  G222 = G222 + tvec4*Eeta[ 88];

  G200 = G200 + tvec5*Eeta[ 41];
  G201 = G201 + tvec5*Eeta[ 47];
  G202 = G202 + tvec5*Eeta[ 53];
  G210 = G210 + tvec5*Eeta[ 59];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec5*Eeta[ 65];
  G212 = G212 + tvec5*Eeta[ 71];
  G220 = G220 + tvec5*Eeta[ 77];
  G221 = G221 + tvec5*Eeta[ 83];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec5*Eeta[ 89];

  H200 = G200;
  H201 = G201;
  H202 = G202;
  H210 = G210;
  H211 = G211;
  H212 = G212;
  H220 = G220;
  H221 = G221;
  H222 = G222;

  sigma[  0][  0] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  1] *        1.0;
  tvec1 = terms[  2] *        2.0;
  tvec2 = terms[  3] *        1.0;
  tvec3 = terms[  7] *        2.0;
  tvec4 = terms[  8] *        2.0;
  tvec5 = terms[ 12] *        1.0;

  G000 = tvec0*Eeta[ 36];
  G001 = tvec0*Eeta[ 42];
  G002 = tvec0*Eeta[ 48];
  G010 = tvec0*Eeta[ 54];
  G011 = tvec0*Eeta[ 60];
  G012 = tvec0*Eeta[ 66];
  G020 = tvec0*Eeta[ 72];
  G021 = tvec0*Eeta[ 78];
  G022 = tvec0*Eeta[ 84];

  G000 = G000 + tvec1*Eeta[ 37];
  G001 = G001 + tvec1*Eeta[ 43];
  G002 = G002 + tvec1*Eeta[ 49];
  G010 = G010 + tvec1*Eeta[ 55];
  G011 = G011 + tvec1*Eeta[ 61];
  G012 = G012 + tvec1*Eeta[ 67];
  G020 = G020 + tvec1*Eeta[ 73];
  G021 = G021 + tvec1*Eeta[ 79];
  G022 = G022 + tvec1*Eeta[ 85];

  G000 = G000 + tvec2*Eeta[ 38];
  G001 = G001 + tvec2*Eeta[ 44];
  G002 = G002 + tvec2*Eeta[ 50];
  G010 = G010 + tvec2*Eeta[ 56];
  G011 = G011 + tvec2*Eeta[ 62];
  G012 = G012 + tvec2*Eeta[ 68];
  G020 = G020 + tvec2*Eeta[ 74];
  G021 = G021 + tvec2*Eeta[ 80];
  G022 = G022 + tvec2*Eeta[ 86];

  G000 = G000 + tvec3*Eeta[ 39];
  G001 = G001 + tvec3*Eeta[ 45];
  G002 = G002 + tvec3*Eeta[ 51];
  G010 = G010 + tvec3*Eeta[ 57];
  G011 = G011 + tvec3*Eeta[ 63];
  G012 = G012 + tvec3*Eeta[ 69];
  G020 = G020 + tvec3*Eeta[ 75];
  G021 = G021 + tvec3*Eeta[ 81];
  G022 = G022 + tvec3*Eeta[ 87];

  G000 = G000 + tvec4*Eeta[ 40];
  G001 = G001 + tvec4*Eeta[ 46];
  G002 = G002 + tvec4*Eeta[ 52];
  G010 = G010 + tvec4*Eeta[ 58];
  G011 = G011 + tvec4*Eeta[ 64];
  G012 = G012 + tvec4*Eeta[ 70];
  G020 = G020 + tvec4*Eeta[ 76];
  G021 = G021 + tvec4*Eeta[ 82];
  G022 = G022 + tvec4*Eeta[ 88];

  G000 = G000 + tvec5*Eeta[ 41];
  G001 = G001 + tvec5*Eeta[ 47];
  G002 = G002 + tvec5*Eeta[ 53];
  G010 = G010 + tvec5*Eeta[ 59];
  G011 = G011 + tvec5*Eeta[ 65];
  G012 = G012 + tvec5*Eeta[ 71];
  G020 = G020 + tvec5*Eeta[ 77];
  G021 = G021 + tvec5*Eeta[ 83];
  G022 = G022 + tvec5*Eeta[ 89];

  tvec0 = terms[  2] *        1.0;
  tvec1 = terms[  3] *        2.0;
  tvec2 = terms[  4] *        1.0;
  tvec3 = terms[  8] *        2.0;
  tvec4 = terms[  9] *        2.0;
  tvec5 = terms[ 13] *        1.0;

  G100 = tvec0*Eeta[ 36];
  G101 = tvec0*Eeta[ 42];
  G102 = tvec0*Eeta[ 48];
  G110 = tvec0*Eeta[ 54];
  G111 = tvec0*Eeta[ 60];
  G112 = tvec0*Eeta[ 66];
  G120 = tvec0*Eeta[ 72];
  G121 = tvec0*Eeta[ 78];
  G122 = tvec0*Eeta[ 84];

  G100 = G100 + tvec1*Eeta[ 37];
  G101 = G101 + tvec1*Eeta[ 43];
  G102 = G102 + tvec1*Eeta[ 49];
  G110 = G110 + tvec1*Eeta[ 55];
  G111 = G111 + tvec1*Eeta[ 61];
  G112 = G112 + tvec1*Eeta[ 67];
  G120 = G120 + tvec1*Eeta[ 73];
  G121 = G121 + tvec1*Eeta[ 79];
  G122 = G122 + tvec1*Eeta[ 85];

  G100 = G100 + tvec2*Eeta[ 38];
  G101 = G101 + tvec2*Eeta[ 44];
  G102 = G102 + tvec2*Eeta[ 50];
  G110 = G110 + tvec2*Eeta[ 56];
  G111 = G111 + tvec2*Eeta[ 62];
  G112 = G112 + tvec2*Eeta[ 68];
  G120 = G120 + tvec2*Eeta[ 74];
  G121 = G121 + tvec2*Eeta[ 80];
  G122 = G122 + tvec2*Eeta[ 86];

  G100 = G100 + tvec3*Eeta[ 39];
  G101 = G101 + tvec3*Eeta[ 45];
  G102 = G102 + tvec3*Eeta[ 51];
  G110 = G110 + tvec3*Eeta[ 57];
  G111 = G111 + tvec3*Eeta[ 63];
  G112 = G112 + tvec3*Eeta[ 69];
  G120 = G120 + tvec3*Eeta[ 75];
  G121 = G121 + tvec3*Eeta[ 81];
  G122 = G122 + tvec3*Eeta[ 87];

  G100 = G100 + tvec4*Eeta[ 40];
  G101 = G101 + tvec4*Eeta[ 46];
  G102 = G102 + tvec4*Eeta[ 52];
  G110 = G110 + tvec4*Eeta[ 58];
  G111 = G111 + tvec4*Eeta[ 64];
  G112 = G112 + tvec4*Eeta[ 70];
  G120 = G120 + tvec4*Eeta[ 76];
  G121 = G121 + tvec4*Eeta[ 82];
  G122 = G122 + tvec4*Eeta[ 88];

  G100 = G100 + tvec5*Eeta[ 41];
  G101 = G101 + tvec5*Eeta[ 47];
  G102 = G102 + tvec5*Eeta[ 53];
  G110 = G110 + tvec5*Eeta[ 59];
  G111 = G111 + tvec5*Eeta[ 65];
  G112 = G112 + tvec5*Eeta[ 71];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec5*Eeta[ 77];
  G121 = G121 + tvec5*Eeta[ 83];
  G122 = G122 + tvec5*Eeta[ 89];

  tvec0 = terms[  7] *        1.0;
  tvec1 = terms[  8] *        2.0;
  tvec2 = terms[  9] *        1.0;
  tvec3 = terms[ 12] *        2.0;
  tvec4 = terms[ 13] *        2.0;
  tvec5 = terms[ 16] *        1.0;

  G200 = tvec0*Eeta[ 36];
  G201 = tvec0*Eeta[ 42];
  G202 = tvec0*Eeta[ 48];
  G210 = tvec0*Eeta[ 54];
  G211 = tvec0*Eeta[ 60];
  G212 = tvec0*Eeta[ 66];
  G220 = tvec0*Eeta[ 72];
  G221 = tvec0*Eeta[ 78];
  G222 = tvec0*Eeta[ 84];

  G200 = G200 + tvec1*Eeta[ 37];
  G201 = G201 + tvec1*Eeta[ 43];
  G202 = G202 + tvec1*Eeta[ 49];
  G210 = G210 + tvec1*Eeta[ 55];
  G211 = G211 + tvec1*Eeta[ 61];
  G212 = G212 + tvec1*Eeta[ 67];
  G220 = G220 + tvec1*Eeta[ 73];
  G221 = G221 + tvec1*Eeta[ 79];
  G222 = G222 + tvec1*Eeta[ 85];

  G200 = G200 + tvec2*Eeta[ 38];
  G201 = G201 + tvec2*Eeta[ 44];
  G202 = G202 + tvec2*Eeta[ 50];
  G210 = G210 + tvec2*Eeta[ 56];
  G211 = G211 + tvec2*Eeta[ 62];
  G212 = G212 + tvec2*Eeta[ 68];
  G220 = G220 + tvec2*Eeta[ 74];
  G221 = G221 + tvec2*Eeta[ 80];
  G222 = G222 + tvec2*Eeta[ 86];

  G200 = G200 + tvec3*Eeta[ 39];
  G201 = G201 + tvec3*Eeta[ 45];
  G202 = G202 + tvec3*Eeta[ 51];
  G210 = G210 + tvec3*Eeta[ 57];
  G211 = G211 + tvec3*Eeta[ 63];
  G212 = G212 + tvec3*Eeta[ 69];
  G220 = G220 + tvec3*Eeta[ 75];
  G221 = G221 + tvec3*Eeta[ 81];
  G222 = G222 + tvec3*Eeta[ 87];

  G200 = G200 + tvec4*Eeta[ 40];
  G201 = G201 + tvec4*Eeta[ 46];
  G202 = G202 + tvec4*Eeta[ 52];
  G210 = G210 + tvec4*Eeta[ 58];
  G211 = G211 + tvec4*Eeta[ 64];
  G212 = G212 + tvec4*Eeta[ 70];
  G220 = G220 + tvec4*Eeta[ 76];
  G221 = G221 + tvec4*Eeta[ 82];
  G222 = G222 + tvec4*Eeta[ 88];

  G200 = G200 + tvec5*Eeta[ 41];
  G201 = G201 + tvec5*Eeta[ 47];
  G202 = G202 + tvec5*Eeta[ 53];
  G210 = G210 + tvec5*Eeta[ 59];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec5*Eeta[ 65];
  G212 = G212 + tvec5*Eeta[ 71];
  G220 = G220 + tvec5*Eeta[ 77];
  G221 = G221 + tvec5*Eeta[ 83];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec5*Eeta[ 89];

  sigma[  0][  1] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  6] *        1.0;
  tvec1 = terms[  7] *        2.0;
  tvec2 = terms[  8] *        1.0;
  tvec3 = terms[ 11] *        2.0;
  tvec4 = terms[ 12] *        2.0;
  tvec5 = terms[ 15] *        1.0;

  G000 = tvec0*Eeta[ 36];
  G001 = tvec0*Eeta[ 42];
  G002 = tvec0*Eeta[ 48];
  G010 = tvec0*Eeta[ 54];
  G011 = tvec0*Eeta[ 60];
  G012 = tvec0*Eeta[ 66];
  G020 = tvec0*Eeta[ 72];
  G021 = tvec0*Eeta[ 78];
  G022 = tvec0*Eeta[ 84];

  G000 = G000 + tvec1*Eeta[ 37];
  G001 = G001 + tvec1*Eeta[ 43];
  G002 = G002 + tvec1*Eeta[ 49];
  G010 = G010 + tvec1*Eeta[ 55];
  G011 = G011 + tvec1*Eeta[ 61];
  G012 = G012 + tvec1*Eeta[ 67];
  G020 = G020 + tvec1*Eeta[ 73];
  G021 = G021 + tvec1*Eeta[ 79];
  G022 = G022 + tvec1*Eeta[ 85];

  G000 = G000 + tvec2*Eeta[ 38];
  G001 = G001 + tvec2*Eeta[ 44];
  G002 = G002 + tvec2*Eeta[ 50];
  G010 = G010 + tvec2*Eeta[ 56];
  G011 = G011 + tvec2*Eeta[ 62];
  G012 = G012 + tvec2*Eeta[ 68];
  G020 = G020 + tvec2*Eeta[ 74];
  G021 = G021 + tvec2*Eeta[ 80];
  G022 = G022 + tvec2*Eeta[ 86];

  G000 = G000 + tvec3*Eeta[ 39];
  G001 = G001 + tvec3*Eeta[ 45];
  G002 = G002 + tvec3*Eeta[ 51];
  G010 = G010 + tvec3*Eeta[ 57];
  G011 = G011 + tvec3*Eeta[ 63];
  G012 = G012 + tvec3*Eeta[ 69];
  G020 = G020 + tvec3*Eeta[ 75];
  G021 = G021 + tvec3*Eeta[ 81];
  G022 = G022 + tvec3*Eeta[ 87];

  G000 = G000 + tvec4*Eeta[ 40];
  G001 = G001 + tvec4*Eeta[ 46];
  G002 = G002 + tvec4*Eeta[ 52];
  G010 = G010 + tvec4*Eeta[ 58];
  G011 = G011 + tvec4*Eeta[ 64];
  G012 = G012 + tvec4*Eeta[ 70];
  G020 = G020 + tvec4*Eeta[ 76];
  G021 = G021 + tvec4*Eeta[ 82];
  G022 = G022 + tvec4*Eeta[ 88];

  G000 = G000 + tvec5*Eeta[ 41];
  G001 = G001 + tvec5*Eeta[ 47];
  G002 = G002 + tvec5*Eeta[ 53];
  G010 = G010 + tvec5*Eeta[ 59];
  G011 = G011 + tvec5*Eeta[ 65];
  G012 = G012 + tvec5*Eeta[ 71];
  G020 = G020 + tvec5*Eeta[ 77];
  G021 = G021 + tvec5*Eeta[ 83];
  G022 = G022 + tvec5*Eeta[ 89];

  tvec0 = terms[  7] *        1.0;
  tvec1 = terms[  8] *        2.0;
  tvec2 = terms[  9] *        1.0;
  tvec3 = terms[ 12] *        2.0;
  tvec4 = terms[ 13] *        2.0;
  tvec5 = terms[ 16] *        1.0;

  G100 = tvec0*Eeta[ 36];
  G101 = tvec0*Eeta[ 42];
  G102 = tvec0*Eeta[ 48];
  G110 = tvec0*Eeta[ 54];
  G111 = tvec0*Eeta[ 60];
  G112 = tvec0*Eeta[ 66];
  G120 = tvec0*Eeta[ 72];
  G121 = tvec0*Eeta[ 78];
  G122 = tvec0*Eeta[ 84];

  G100 = G100 + tvec1*Eeta[ 37];
  G101 = G101 + tvec1*Eeta[ 43];
  G102 = G102 + tvec1*Eeta[ 49];
  G110 = G110 + tvec1*Eeta[ 55];
  G111 = G111 + tvec1*Eeta[ 61];
  G112 = G112 + tvec1*Eeta[ 67];
  G120 = G120 + tvec1*Eeta[ 73];
  G121 = G121 + tvec1*Eeta[ 79];
  G122 = G122 + tvec1*Eeta[ 85];

  G100 = G100 + tvec2*Eeta[ 38];
  G101 = G101 + tvec2*Eeta[ 44];
  G102 = G102 + tvec2*Eeta[ 50];
  G110 = G110 + tvec2*Eeta[ 56];
  G111 = G111 + tvec2*Eeta[ 62];
  G112 = G112 + tvec2*Eeta[ 68];
  G120 = G120 + tvec2*Eeta[ 74];
  G121 = G121 + tvec2*Eeta[ 80];
  G122 = G122 + tvec2*Eeta[ 86];

  G100 = G100 + tvec3*Eeta[ 39];
  G101 = G101 + tvec3*Eeta[ 45];
  G102 = G102 + tvec3*Eeta[ 51];
  G110 = G110 + tvec3*Eeta[ 57];
  G111 = G111 + tvec3*Eeta[ 63];
  G112 = G112 + tvec3*Eeta[ 69];
  G120 = G120 + tvec3*Eeta[ 75];
  G121 = G121 + tvec3*Eeta[ 81];
  G122 = G122 + tvec3*Eeta[ 87];

  G100 = G100 + tvec4*Eeta[ 40];
  G101 = G101 + tvec4*Eeta[ 46];
  G102 = G102 + tvec4*Eeta[ 52];
  G110 = G110 + tvec4*Eeta[ 58];
  G111 = G111 + tvec4*Eeta[ 64];
  G112 = G112 + tvec4*Eeta[ 70];
  G120 = G120 + tvec4*Eeta[ 76];
  G121 = G121 + tvec4*Eeta[ 82];
  G122 = G122 + tvec4*Eeta[ 88];

  G100 = G100 + tvec5*Eeta[ 41];
  G101 = G101 + tvec5*Eeta[ 47];
  G102 = G102 + tvec5*Eeta[ 53];
  G110 = G110 + tvec5*Eeta[ 59];
  G111 = G111 + tvec5*Eeta[ 65];
  G112 = G112 + tvec5*Eeta[ 71];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec5*Eeta[ 77];
  G121 = G121 + tvec5*Eeta[ 83];
  G122 = G122 + tvec5*Eeta[ 89];

  tvec0 = terms[ 11] *        1.0;
  tvec1 = terms[ 12] *        2.0;
  tvec2 = terms[ 13] *        1.0;
  tvec3 = terms[ 15] *        2.0;
  tvec4 = terms[ 16] *        2.0;
  tvec5 = terms[ 18] *        1.0;

  G200 = tvec0*Eeta[ 36];
  G201 = tvec0*Eeta[ 42];
  G202 = tvec0*Eeta[ 48];
  G210 = tvec0*Eeta[ 54];
  G211 = tvec0*Eeta[ 60];
  G212 = tvec0*Eeta[ 66];
  G220 = tvec0*Eeta[ 72];
  G221 = tvec0*Eeta[ 78];
  G222 = tvec0*Eeta[ 84];

  G200 = G200 + tvec1*Eeta[ 37];
  G201 = G201 + tvec1*Eeta[ 43];
  G202 = G202 + tvec1*Eeta[ 49];
  G210 = G210 + tvec1*Eeta[ 55];
  G211 = G211 + tvec1*Eeta[ 61];
  G212 = G212 + tvec1*Eeta[ 67];
  G220 = G220 + tvec1*Eeta[ 73];
  G221 = G221 + tvec1*Eeta[ 79];
  G222 = G222 + tvec1*Eeta[ 85];

  G200 = G200 + tvec2*Eeta[ 38];
  G201 = G201 + tvec2*Eeta[ 44];
  G202 = G202 + tvec2*Eeta[ 50];
  G210 = G210 + tvec2*Eeta[ 56];
  G211 = G211 + tvec2*Eeta[ 62];
  G212 = G212 + tvec2*Eeta[ 68];
  G220 = G220 + tvec2*Eeta[ 74];
  G221 = G221 + tvec2*Eeta[ 80];
  G222 = G222 + tvec2*Eeta[ 86];

  G200 = G200 + tvec3*Eeta[ 39];
  G201 = G201 + tvec3*Eeta[ 45];
  G202 = G202 + tvec3*Eeta[ 51];
  G210 = G210 + tvec3*Eeta[ 57];
  G211 = G211 + tvec3*Eeta[ 63];
  G212 = G212 + tvec3*Eeta[ 69];
  G220 = G220 + tvec3*Eeta[ 75];
  G221 = G221 + tvec3*Eeta[ 81];
  G222 = G222 + tvec3*Eeta[ 87];

  G200 = G200 + tvec4*Eeta[ 40];
  G201 = G201 + tvec4*Eeta[ 46];
  G202 = G202 + tvec4*Eeta[ 52];
  G210 = G210 + tvec4*Eeta[ 58];
  G211 = G211 + tvec4*Eeta[ 64];
  G212 = G212 + tvec4*Eeta[ 70];
  G220 = G220 + tvec4*Eeta[ 76];
  G221 = G221 + tvec4*Eeta[ 82];
  G222 = G222 + tvec4*Eeta[ 88];

  G200 = G200 + tvec5*Eeta[ 41];
  G201 = G201 + tvec5*Eeta[ 47];
  G202 = G202 + tvec5*Eeta[ 53];
  G210 = G210 + tvec5*Eeta[ 59];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec5*Eeta[ 65];
  G212 = G212 + tvec5*Eeta[ 71];
  G220 = G220 + tvec5*Eeta[ 77];
  G221 = G221 + tvec5*Eeta[ 83];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec5*Eeta[ 89];

  sigma[  0][  2] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  2] *        1.0;
  tvec1 = terms[  3] *        2.0;
  tvec2 = terms[  4] *        1.0;
  tvec3 = terms[  8] *        2.0;
  tvec4 = terms[  9] *        2.0;
  tvec5 = terms[ 13] *        1.0;

  G000 = tvec0*Eeta[ 36];
  G001 = tvec0*Eeta[ 42];
  G002 = tvec0*Eeta[ 48];
  G010 = tvec0*Eeta[ 54];
  G011 = tvec0*Eeta[ 60];
  G012 = tvec0*Eeta[ 66];
  G020 = tvec0*Eeta[ 72];
  G021 = tvec0*Eeta[ 78];
  G022 = tvec0*Eeta[ 84];

  G000 = G000 + tvec1*Eeta[ 37];
  G001 = G001 + tvec1*Eeta[ 43];
  G002 = G002 + tvec1*Eeta[ 49];
  G010 = G010 + tvec1*Eeta[ 55];
  G011 = G011 + tvec1*Eeta[ 61];
  G012 = G012 + tvec1*Eeta[ 67];
  G020 = G020 + tvec1*Eeta[ 73];
  G021 = G021 + tvec1*Eeta[ 79];
  G022 = G022 + tvec1*Eeta[ 85];

  G000 = G000 + tvec2*Eeta[ 38];
  G001 = G001 + tvec2*Eeta[ 44];
  G002 = G002 + tvec2*Eeta[ 50];
  G010 = G010 + tvec2*Eeta[ 56];
  G011 = G011 + tvec2*Eeta[ 62];
  G012 = G012 + tvec2*Eeta[ 68];
  G020 = G020 + tvec2*Eeta[ 74];
  G021 = G021 + tvec2*Eeta[ 80];
  G022 = G022 + tvec2*Eeta[ 86];

  G000 = G000 + tvec3*Eeta[ 39];
  G001 = G001 + tvec3*Eeta[ 45];
  G002 = G002 + tvec3*Eeta[ 51];
  G010 = G010 + tvec3*Eeta[ 57];
  G011 = G011 + tvec3*Eeta[ 63];
  G012 = G012 + tvec3*Eeta[ 69];
  G020 = G020 + tvec3*Eeta[ 75];
  G021 = G021 + tvec3*Eeta[ 81];
  G022 = G022 + tvec3*Eeta[ 87];

  G000 = G000 + tvec4*Eeta[ 40];
  G001 = G001 + tvec4*Eeta[ 46];
  G002 = G002 + tvec4*Eeta[ 52];
  G010 = G010 + tvec4*Eeta[ 58];
  G011 = G011 + tvec4*Eeta[ 64];
  G012 = G012 + tvec4*Eeta[ 70];
  G020 = G020 + tvec4*Eeta[ 76];
  G021 = G021 + tvec4*Eeta[ 82];
  G022 = G022 + tvec4*Eeta[ 88];

  G000 = G000 + tvec5*Eeta[ 41];
  G001 = G001 + tvec5*Eeta[ 47];
  G002 = G002 + tvec5*Eeta[ 53];
  G010 = G010 + tvec5*Eeta[ 59];
  G011 = G011 + tvec5*Eeta[ 65];
  G012 = G012 + tvec5*Eeta[ 71];
  G020 = G020 + tvec5*Eeta[ 77];
  G021 = G021 + tvec5*Eeta[ 83];
  G022 = G022 + tvec5*Eeta[ 89];

  H000 = H000 + G000;
  H001 = H001 + G001;
  H002 = H002 + G002;
  H010 = H010 + G010;
  H011 = H011 + G011;
  H012 = H012 + G012;
  H020 = H020 + G020;
  H021 = H021 + G021;
  H022 = H022 + G022;

  tvec0 = terms[  3] *        1.0;
  tvec1 = terms[  4] *        2.0;
  tvec2 = terms[  5] *        1.0;
  tvec3 = terms[  9] *        2.0;
  tvec4 = terms[ 10] *        2.0;
  tvec5 = terms[ 14] *        1.0;

  G100 = tvec0*Eeta[ 36];
  G101 = tvec0*Eeta[ 42];
  G102 = tvec0*Eeta[ 48];
  G110 = tvec0*Eeta[ 54];
  G111 = tvec0*Eeta[ 60];
  G112 = tvec0*Eeta[ 66];
  G120 = tvec0*Eeta[ 72];
  G121 = tvec0*Eeta[ 78];
  G122 = tvec0*Eeta[ 84];

  G100 = G100 + tvec1*Eeta[ 37];
  G101 = G101 + tvec1*Eeta[ 43];
  G102 = G102 + tvec1*Eeta[ 49];
  G110 = G110 + tvec1*Eeta[ 55];
  G111 = G111 + tvec1*Eeta[ 61];
  G112 = G112 + tvec1*Eeta[ 67];
  G120 = G120 + tvec1*Eeta[ 73];
  G121 = G121 + tvec1*Eeta[ 79];
  G122 = G122 + tvec1*Eeta[ 85];

  G100 = G100 + tvec2*Eeta[ 38];
  G101 = G101 + tvec2*Eeta[ 44];
  G102 = G102 + tvec2*Eeta[ 50];
  G110 = G110 + tvec2*Eeta[ 56];
  G111 = G111 + tvec2*Eeta[ 62];
  G112 = G112 + tvec2*Eeta[ 68];
  G120 = G120 + tvec2*Eeta[ 74];
  G121 = G121 + tvec2*Eeta[ 80];
  G122 = G122 + tvec2*Eeta[ 86];

  G100 = G100 + tvec3*Eeta[ 39];
  G101 = G101 + tvec3*Eeta[ 45];
  G102 = G102 + tvec3*Eeta[ 51];
  G110 = G110 + tvec3*Eeta[ 57];
  G111 = G111 + tvec3*Eeta[ 63];
  G112 = G112 + tvec3*Eeta[ 69];
  G120 = G120 + tvec3*Eeta[ 75];
  G121 = G121 + tvec3*Eeta[ 81];
  G122 = G122 + tvec3*Eeta[ 87];

  G100 = G100 + tvec4*Eeta[ 40];
  G101 = G101 + tvec4*Eeta[ 46];
  G102 = G102 + tvec4*Eeta[ 52];
  G110 = G110 + tvec4*Eeta[ 58];
  G111 = G111 + tvec4*Eeta[ 64];
  G112 = G112 + tvec4*Eeta[ 70];
  G120 = G120 + tvec4*Eeta[ 76];
  G121 = G121 + tvec4*Eeta[ 82];
  G122 = G122 + tvec4*Eeta[ 88];

  G100 = G100 + tvec5*Eeta[ 41];
  G101 = G101 + tvec5*Eeta[ 47];
  G102 = G102 + tvec5*Eeta[ 53];
  G110 = G110 + tvec5*Eeta[ 59];
  G111 = G111 + tvec5*Eeta[ 65];
  G112 = G112 + tvec5*Eeta[ 71];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec5*Eeta[ 77];
  G121 = G121 + tvec5*Eeta[ 83];
  G122 = G122 + tvec5*Eeta[ 89];

  H100 = H100 + G100;
  H101 = H101 + G101;
  H102 = H102 + G102;
  H110 = H110 + G110;
  H111 = H111 + G111;
  H112 = H112 + G112;
  H120 = H120 + G120;
  H121 = H121 + G121;
  H122 = H122 + G122;

  tvec0 = terms[  8] *        1.0;
  tvec1 = terms[  9] *        2.0;
  tvec2 = terms[ 10] *        1.0;
  tvec3 = terms[ 13] *        2.0;
  tvec4 = terms[ 14] *        2.0;
  tvec5 = terms[ 17] *        1.0;

  G200 = tvec0*Eeta[ 36];
  G201 = tvec0*Eeta[ 42];
  G202 = tvec0*Eeta[ 48];
  G210 = tvec0*Eeta[ 54];
  G211 = tvec0*Eeta[ 60];
  G212 = tvec0*Eeta[ 66];
  G220 = tvec0*Eeta[ 72];
  G221 = tvec0*Eeta[ 78];
  G222 = tvec0*Eeta[ 84];

  G200 = G200 + tvec1*Eeta[ 37];
  G201 = G201 + tvec1*Eeta[ 43];
  G202 = G202 + tvec1*Eeta[ 49];
  G210 = G210 + tvec1*Eeta[ 55];
  G211 = G211 + tvec1*Eeta[ 61];
  G212 = G212 + tvec1*Eeta[ 67];
  G220 = G220 + tvec1*Eeta[ 73];
  G221 = G221 + tvec1*Eeta[ 79];
  G222 = G222 + tvec1*Eeta[ 85];

  G200 = G200 + tvec2*Eeta[ 38];
  G201 = G201 + tvec2*Eeta[ 44];
  G202 = G202 + tvec2*Eeta[ 50];
  G210 = G210 + tvec2*Eeta[ 56];
  G211 = G211 + tvec2*Eeta[ 62];
  G212 = G212 + tvec2*Eeta[ 68];
  G220 = G220 + tvec2*Eeta[ 74];
  G221 = G221 + tvec2*Eeta[ 80];
  G222 = G222 + tvec2*Eeta[ 86];

  G200 = G200 + tvec3*Eeta[ 39];
  G201 = G201 + tvec3*Eeta[ 45];
  G202 = G202 + tvec3*Eeta[ 51];
  G210 = G210 + tvec3*Eeta[ 57];
  G211 = G211 + tvec3*Eeta[ 63];
  G212 = G212 + tvec3*Eeta[ 69];
  G220 = G220 + tvec3*Eeta[ 75];
  G221 = G221 + tvec3*Eeta[ 81];
  G222 = G222 + tvec3*Eeta[ 87];

  G200 = G200 + tvec4*Eeta[ 40];
  G201 = G201 + tvec4*Eeta[ 46];
  G202 = G202 + tvec4*Eeta[ 52];
  G210 = G210 + tvec4*Eeta[ 58];
  G211 = G211 + tvec4*Eeta[ 64];
  G212 = G212 + tvec4*Eeta[ 70];
  G220 = G220 + tvec4*Eeta[ 76];
  G221 = G221 + tvec4*Eeta[ 82];
  G222 = G222 + tvec4*Eeta[ 88];

  G200 = G200 + tvec5*Eeta[ 41];
  G201 = G201 + tvec5*Eeta[ 47];
  G202 = G202 + tvec5*Eeta[ 53];
  G210 = G210 + tvec5*Eeta[ 59];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec5*Eeta[ 65];
  G212 = G212 + tvec5*Eeta[ 71];
  G220 = G220 + tvec5*Eeta[ 77];
  G221 = G221 + tvec5*Eeta[ 83];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec5*Eeta[ 89];

  H200 = H200 + G200;
  H201 = H201 + G201;
  H202 = H202 + G202;
  H210 = H210 + G210;
  H211 = H211 + G211;
  H212 = H212 + G212;
  H220 = H220 + G220;
  H221 = H221 + G221;
  H222 = H222 + G222;

  sigma[  1][  1] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  7] *        1.0;
  tvec1 = terms[  8] *        2.0;
  tvec2 = terms[  9] *        1.0;
  tvec3 = terms[ 12] *        2.0;
  tvec4 = terms[ 13] *        2.0;
  tvec5 = terms[ 16] *        1.0;

  G000 = tvec0*Eeta[ 36];
  G001 = tvec0*Eeta[ 42];
  G002 = tvec0*Eeta[ 48];
  G010 = tvec0*Eeta[ 54];
  G011 = tvec0*Eeta[ 60];
  G012 = tvec0*Eeta[ 66];
  G020 = tvec0*Eeta[ 72];
  G021 = tvec0*Eeta[ 78];
  G022 = tvec0*Eeta[ 84];

  G000 = G000 + tvec1*Eeta[ 37];
  G001 = G001 + tvec1*Eeta[ 43];
  G002 = G002 + tvec1*Eeta[ 49];
  G010 = G010 + tvec1*Eeta[ 55];
  G011 = G011 + tvec1*Eeta[ 61];
  G012 = G012 + tvec1*Eeta[ 67];
  G020 = G020 + tvec1*Eeta[ 73];
  G021 = G021 + tvec1*Eeta[ 79];
  G022 = G022 + tvec1*Eeta[ 85];

  G000 = G000 + tvec2*Eeta[ 38];
  G001 = G001 + tvec2*Eeta[ 44];
  G002 = G002 + tvec2*Eeta[ 50];
  G010 = G010 + tvec2*Eeta[ 56];
  G011 = G011 + tvec2*Eeta[ 62];
  G012 = G012 + tvec2*Eeta[ 68];
  G020 = G020 + tvec2*Eeta[ 74];
  G021 = G021 + tvec2*Eeta[ 80];
  G022 = G022 + tvec2*Eeta[ 86];

  G000 = G000 + tvec3*Eeta[ 39];
  G001 = G001 + tvec3*Eeta[ 45];
  G002 = G002 + tvec3*Eeta[ 51];
  G010 = G010 + tvec3*Eeta[ 57];
  G011 = G011 + tvec3*Eeta[ 63];
  G012 = G012 + tvec3*Eeta[ 69];
  G020 = G020 + tvec3*Eeta[ 75];
  G021 = G021 + tvec3*Eeta[ 81];
  G022 = G022 + tvec3*Eeta[ 87];

  G000 = G000 + tvec4*Eeta[ 40];
  G001 = G001 + tvec4*Eeta[ 46];
  G002 = G002 + tvec4*Eeta[ 52];
  G010 = G010 + tvec4*Eeta[ 58];
  G011 = G011 + tvec4*Eeta[ 64];
  G012 = G012 + tvec4*Eeta[ 70];
  G020 = G020 + tvec4*Eeta[ 76];
  G021 = G021 + tvec4*Eeta[ 82];
  G022 = G022 + tvec4*Eeta[ 88];

  G000 = G000 + tvec5*Eeta[ 41];
  G001 = G001 + tvec5*Eeta[ 47];
  G002 = G002 + tvec5*Eeta[ 53];
  G010 = G010 + tvec5*Eeta[ 59];
  G011 = G011 + tvec5*Eeta[ 65];
  G012 = G012 + tvec5*Eeta[ 71];
  G020 = G020 + tvec5*Eeta[ 77];
  G021 = G021 + tvec5*Eeta[ 83];
  G022 = G022 + tvec5*Eeta[ 89];

  tvec0 = terms[  8] *        1.0;
  tvec1 = terms[  9] *        2.0;
  tvec2 = terms[ 10] *        1.0;
  tvec3 = terms[ 13] *        2.0;
  tvec4 = terms[ 14] *        2.0;
  tvec5 = terms[ 17] *        1.0;

  G100 = tvec0*Eeta[ 36];
  G101 = tvec0*Eeta[ 42];
  G102 = tvec0*Eeta[ 48];
  G110 = tvec0*Eeta[ 54];
  G111 = tvec0*Eeta[ 60];
  G112 = tvec0*Eeta[ 66];
  G120 = tvec0*Eeta[ 72];
  G121 = tvec0*Eeta[ 78];
  G122 = tvec0*Eeta[ 84];

  G100 = G100 + tvec1*Eeta[ 37];
  G101 = G101 + tvec1*Eeta[ 43];
  G102 = G102 + tvec1*Eeta[ 49];
  G110 = G110 + tvec1*Eeta[ 55];
  G111 = G111 + tvec1*Eeta[ 61];
  G112 = G112 + tvec1*Eeta[ 67];
  G120 = G120 + tvec1*Eeta[ 73];
  G121 = G121 + tvec1*Eeta[ 79];
  G122 = G122 + tvec1*Eeta[ 85];

  G100 = G100 + tvec2*Eeta[ 38];
  G101 = G101 + tvec2*Eeta[ 44];
  G102 = G102 + tvec2*Eeta[ 50];
  G110 = G110 + tvec2*Eeta[ 56];
  G111 = G111 + tvec2*Eeta[ 62];
  G112 = G112 + tvec2*Eeta[ 68];
  G120 = G120 + tvec2*Eeta[ 74];
  G121 = G121 + tvec2*Eeta[ 80];
  G122 = G122 + tvec2*Eeta[ 86];

  G100 = G100 + tvec3*Eeta[ 39];
  G101 = G101 + tvec3*Eeta[ 45];
  G102 = G102 + tvec3*Eeta[ 51];
  G110 = G110 + tvec3*Eeta[ 57];
  G111 = G111 + tvec3*Eeta[ 63];
  G112 = G112 + tvec3*Eeta[ 69];
  G120 = G120 + tvec3*Eeta[ 75];
  G121 = G121 + tvec3*Eeta[ 81];
  G122 = G122 + tvec3*Eeta[ 87];

  G100 = G100 + tvec4*Eeta[ 40];
  G101 = G101 + tvec4*Eeta[ 46];
  G102 = G102 + tvec4*Eeta[ 52];
  G110 = G110 + tvec4*Eeta[ 58];
  G111 = G111 + tvec4*Eeta[ 64];
  G112 = G112 + tvec4*Eeta[ 70];
  G120 = G120 + tvec4*Eeta[ 76];
  G121 = G121 + tvec4*Eeta[ 82];
  G122 = G122 + tvec4*Eeta[ 88];

  G100 = G100 + tvec5*Eeta[ 41];
  G101 = G101 + tvec5*Eeta[ 47];
  G102 = G102 + tvec5*Eeta[ 53];
  G110 = G110 + tvec5*Eeta[ 59];
  G111 = G111 + tvec5*Eeta[ 65];
  G112 = G112 + tvec5*Eeta[ 71];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec5*Eeta[ 77];
  G121 = G121 + tvec5*Eeta[ 83];
  G122 = G122 + tvec5*Eeta[ 89];

  tvec0 = terms[ 12] *        1.0;
  tvec1 = terms[ 13] *        2.0;
  tvec2 = terms[ 14] *        1.0;
  tvec3 = terms[ 16] *        2.0;
  tvec4 = terms[ 17] *        2.0;
  tvec5 = terms[ 19] *        1.0;

  G200 = tvec0*Eeta[ 36];
  G201 = tvec0*Eeta[ 42];
  G202 = tvec0*Eeta[ 48];
  G210 = tvec0*Eeta[ 54];
  G211 = tvec0*Eeta[ 60];
  G212 = tvec0*Eeta[ 66];
  G220 = tvec0*Eeta[ 72];
  G221 = tvec0*Eeta[ 78];
  G222 = tvec0*Eeta[ 84];

  G200 = G200 + tvec1*Eeta[ 37];
  G201 = G201 + tvec1*Eeta[ 43];
  G202 = G202 + tvec1*Eeta[ 49];
  G210 = G210 + tvec1*Eeta[ 55];
  G211 = G211 + tvec1*Eeta[ 61];
  G212 = G212 + tvec1*Eeta[ 67];
  G220 = G220 + tvec1*Eeta[ 73];
  G221 = G221 + tvec1*Eeta[ 79];
  G222 = G222 + tvec1*Eeta[ 85];

  G200 = G200 + tvec2*Eeta[ 38];
  G201 = G201 + tvec2*Eeta[ 44];
  G202 = G202 + tvec2*Eeta[ 50];
  G210 = G210 + tvec2*Eeta[ 56];
  G211 = G211 + tvec2*Eeta[ 62];
  G212 = G212 + tvec2*Eeta[ 68];
  G220 = G220 + tvec2*Eeta[ 74];
  G221 = G221 + tvec2*Eeta[ 80];
  G222 = G222 + tvec2*Eeta[ 86];

  G200 = G200 + tvec3*Eeta[ 39];
  G201 = G201 + tvec3*Eeta[ 45];
  G202 = G202 + tvec3*Eeta[ 51];
  G210 = G210 + tvec3*Eeta[ 57];
  G211 = G211 + tvec3*Eeta[ 63];
  G212 = G212 + tvec3*Eeta[ 69];
  G220 = G220 + tvec3*Eeta[ 75];
  G221 = G221 + tvec3*Eeta[ 81];
  G222 = G222 + tvec3*Eeta[ 87];

  G200 = G200 + tvec4*Eeta[ 40];
  G201 = G201 + tvec4*Eeta[ 46];
  G202 = G202 + tvec4*Eeta[ 52];
  G210 = G210 + tvec4*Eeta[ 58];
  G211 = G211 + tvec4*Eeta[ 64];
  G212 = G212 + tvec4*Eeta[ 70];
  G220 = G220 + tvec4*Eeta[ 76];
  G221 = G221 + tvec4*Eeta[ 82];
  G222 = G222 + tvec4*Eeta[ 88];

  G200 = G200 + tvec5*Eeta[ 41];
  G201 = G201 + tvec5*Eeta[ 47];
  G202 = G202 + tvec5*Eeta[ 53];
  G210 = G210 + tvec5*Eeta[ 59];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec5*Eeta[ 65];
  G212 = G212 + tvec5*Eeta[ 71];
  G220 = G220 + tvec5*Eeta[ 77];
  G221 = G221 + tvec5*Eeta[ 83];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec5*Eeta[ 89];

  sigma[  1][  2] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[ 11] *        1.0;
  tvec1 = terms[ 12] *        2.0;
  tvec2 = terms[ 13] *        1.0;
  tvec3 = terms[ 15] *        2.0;
  tvec4 = terms[ 16] *        2.0;
  tvec5 = terms[ 18] *        1.0;

  G000 = tvec0*Eeta[ 36];
  G001 = tvec0*Eeta[ 42];
  G002 = tvec0*Eeta[ 48];
  G010 = tvec0*Eeta[ 54];
  G011 = tvec0*Eeta[ 60];
  G012 = tvec0*Eeta[ 66];
  G020 = tvec0*Eeta[ 72];
  G021 = tvec0*Eeta[ 78];
  G022 = tvec0*Eeta[ 84];

  G000 = G000 + tvec1*Eeta[ 37];
  G001 = G001 + tvec1*Eeta[ 43];
  G002 = G002 + tvec1*Eeta[ 49];
  G010 = G010 + tvec1*Eeta[ 55];
  G011 = G011 + tvec1*Eeta[ 61];
  G012 = G012 + tvec1*Eeta[ 67];
  G020 = G020 + tvec1*Eeta[ 73];
  G021 = G021 + tvec1*Eeta[ 79];
  G022 = G022 + tvec1*Eeta[ 85];

  G000 = G000 + tvec2*Eeta[ 38];
  G001 = G001 + tvec2*Eeta[ 44];
  G002 = G002 + tvec2*Eeta[ 50];
  G010 = G010 + tvec2*Eeta[ 56];
  G011 = G011 + tvec2*Eeta[ 62];
  G012 = G012 + tvec2*Eeta[ 68];
  G020 = G020 + tvec2*Eeta[ 74];
  G021 = G021 + tvec2*Eeta[ 80];
  G022 = G022 + tvec2*Eeta[ 86];

  G000 = G000 + tvec3*Eeta[ 39];
  G001 = G001 + tvec3*Eeta[ 45];
  G002 = G002 + tvec3*Eeta[ 51];
  G010 = G010 + tvec3*Eeta[ 57];
  G011 = G011 + tvec3*Eeta[ 63];
  G012 = G012 + tvec3*Eeta[ 69];
  G020 = G020 + tvec3*Eeta[ 75];
  G021 = G021 + tvec3*Eeta[ 81];
  G022 = G022 + tvec3*Eeta[ 87];

  G000 = G000 + tvec4*Eeta[ 40];
  G001 = G001 + tvec4*Eeta[ 46];
  G002 = G002 + tvec4*Eeta[ 52];
  G010 = G010 + tvec4*Eeta[ 58];
  G011 = G011 + tvec4*Eeta[ 64];
  G012 = G012 + tvec4*Eeta[ 70];
  G020 = G020 + tvec4*Eeta[ 76];
  G021 = G021 + tvec4*Eeta[ 82];
  G022 = G022 + tvec4*Eeta[ 88];

  G000 = G000 + tvec5*Eeta[ 41];
  G001 = G001 + tvec5*Eeta[ 47];
  G002 = G002 + tvec5*Eeta[ 53];
  G010 = G010 + tvec5*Eeta[ 59];
  G011 = G011 + tvec5*Eeta[ 65];
  G012 = G012 + tvec5*Eeta[ 71];
  G020 = G020 + tvec5*Eeta[ 77];
  G021 = G021 + tvec5*Eeta[ 83];
  G022 = G022 + tvec5*Eeta[ 89];

  H000 = H000 + G000;
  H001 = H001 + G001;
  H002 = H002 + G002;
  H010 = H010 + G010;
  H011 = H011 + G011;
  H012 = H012 + G012;
  H020 = H020 + G020;
  H021 = H021 + G021;
  H022 = H022 + G022;

  tvec0 = terms[ 12] *        1.0;
  tvec1 = terms[ 13] *        2.0;
  tvec2 = terms[ 14] *        1.0;
  tvec3 = terms[ 16] *        2.0;
  tvec4 = terms[ 17] *        2.0;
  tvec5 = terms[ 19] *        1.0;

  G100 = tvec0*Eeta[ 36];
  G101 = tvec0*Eeta[ 42];
  G102 = tvec0*Eeta[ 48];
  G110 = tvec0*Eeta[ 54];
  G111 = tvec0*Eeta[ 60];
  G112 = tvec0*Eeta[ 66];
  G120 = tvec0*Eeta[ 72];
  G121 = tvec0*Eeta[ 78];
  G122 = tvec0*Eeta[ 84];

  G100 = G100 + tvec1*Eeta[ 37];
  G101 = G101 + tvec1*Eeta[ 43];
  G102 = G102 + tvec1*Eeta[ 49];
  G110 = G110 + tvec1*Eeta[ 55];
  G111 = G111 + tvec1*Eeta[ 61];
  G112 = G112 + tvec1*Eeta[ 67];
  G120 = G120 + tvec1*Eeta[ 73];
  G121 = G121 + tvec1*Eeta[ 79];
  G122 = G122 + tvec1*Eeta[ 85];

  G100 = G100 + tvec2*Eeta[ 38];
  G101 = G101 + tvec2*Eeta[ 44];
  G102 = G102 + tvec2*Eeta[ 50];
  G110 = G110 + tvec2*Eeta[ 56];
  G111 = G111 + tvec2*Eeta[ 62];
  G112 = G112 + tvec2*Eeta[ 68];
  G120 = G120 + tvec2*Eeta[ 74];
  G121 = G121 + tvec2*Eeta[ 80];
  G122 = G122 + tvec2*Eeta[ 86];

  G100 = G100 + tvec3*Eeta[ 39];
  G101 = G101 + tvec3*Eeta[ 45];
  G102 = G102 + tvec3*Eeta[ 51];
  G110 = G110 + tvec3*Eeta[ 57];
  G111 = G111 + tvec3*Eeta[ 63];
  G112 = G112 + tvec3*Eeta[ 69];
  G120 = G120 + tvec3*Eeta[ 75];
  G121 = G121 + tvec3*Eeta[ 81];
  G122 = G122 + tvec3*Eeta[ 87];

  G100 = G100 + tvec4*Eeta[ 40];
  G101 = G101 + tvec4*Eeta[ 46];
  G102 = G102 + tvec4*Eeta[ 52];
  G110 = G110 + tvec4*Eeta[ 58];
  G111 = G111 + tvec4*Eeta[ 64];
  G112 = G112 + tvec4*Eeta[ 70];
  G120 = G120 + tvec4*Eeta[ 76];
  G121 = G121 + tvec4*Eeta[ 82];
  G122 = G122 + tvec4*Eeta[ 88];

  G100 = G100 + tvec5*Eeta[ 41];
  G101 = G101 + tvec5*Eeta[ 47];
  G102 = G102 + tvec5*Eeta[ 53];
  G110 = G110 + tvec5*Eeta[ 59];
  G111 = G111 + tvec5*Eeta[ 65];
  G112 = G112 + tvec5*Eeta[ 71];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec5*Eeta[ 77];
  G121 = G121 + tvec5*Eeta[ 83];
  G122 = G122 + tvec5*Eeta[ 89];

  H100 = H100 + G100;
  H101 = H101 + G101;
  H102 = H102 + G102;
  H110 = H110 + G110;
  H111 = H111 + G111;
  H112 = H112 + G112;
  H120 = H120 + G120;
  H121 = H121 + G121;
  H122 = H122 + G122;

  tvec0 = terms[ 15] *        1.0;
  tvec1 = terms[ 16] *        2.0;
  tvec2 = terms[ 17] *        1.0;
  tvec3 = terms[ 18] *        2.0;
  tvec4 = terms[ 19] *        2.0;
  tvec5 = terms[ 20] *        1.0;

  G200 = tvec0*Eeta[ 36];
  G201 = tvec0*Eeta[ 42];
  G202 = tvec0*Eeta[ 48];
  G210 = tvec0*Eeta[ 54];
  G211 = tvec0*Eeta[ 60];
  G212 = tvec0*Eeta[ 66];
  G220 = tvec0*Eeta[ 72];
  G221 = tvec0*Eeta[ 78];
  G222 = tvec0*Eeta[ 84];

  G200 = G200 + tvec1*Eeta[ 37];
  G201 = G201 + tvec1*Eeta[ 43];
  G202 = G202 + tvec1*Eeta[ 49];
  G210 = G210 + tvec1*Eeta[ 55];
  G211 = G211 + tvec1*Eeta[ 61];
  G212 = G212 + tvec1*Eeta[ 67];
  G220 = G220 + tvec1*Eeta[ 73];
  G221 = G221 + tvec1*Eeta[ 79];
  G222 = G222 + tvec1*Eeta[ 85];

  G200 = G200 + tvec2*Eeta[ 38];
  G201 = G201 + tvec2*Eeta[ 44];
  G202 = G202 + tvec2*Eeta[ 50];
  G210 = G210 + tvec2*Eeta[ 56];
  G211 = G211 + tvec2*Eeta[ 62];
  G212 = G212 + tvec2*Eeta[ 68];
  G220 = G220 + tvec2*Eeta[ 74];
  G221 = G221 + tvec2*Eeta[ 80];
  G222 = G222 + tvec2*Eeta[ 86];

  G200 = G200 + tvec3*Eeta[ 39];
  G201 = G201 + tvec3*Eeta[ 45];
  G202 = G202 + tvec3*Eeta[ 51];
  G210 = G210 + tvec3*Eeta[ 57];
  G211 = G211 + tvec3*Eeta[ 63];
  G212 = G212 + tvec3*Eeta[ 69];
  G220 = G220 + tvec3*Eeta[ 75];
  G221 = G221 + tvec3*Eeta[ 81];
  G222 = G222 + tvec3*Eeta[ 87];

  G200 = G200 + tvec4*Eeta[ 40];
  G201 = G201 + tvec4*Eeta[ 46];
  G202 = G202 + tvec4*Eeta[ 52];
  G210 = G210 + tvec4*Eeta[ 58];
  G211 = G211 + tvec4*Eeta[ 64];
  G212 = G212 + tvec4*Eeta[ 70];
  G220 = G220 + tvec4*Eeta[ 76];
  G221 = G221 + tvec4*Eeta[ 82];
  G222 = G222 + tvec4*Eeta[ 88];

  G200 = G200 + tvec5*Eeta[ 41];
  G201 = G201 + tvec5*Eeta[ 47];
  G202 = G202 + tvec5*Eeta[ 53];
  G210 = G210 + tvec5*Eeta[ 59];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec5*Eeta[ 65];
  G212 = G212 + tvec5*Eeta[ 71];
  G220 = G220 + tvec5*Eeta[ 77];
  G221 = G221 + tvec5*Eeta[ 83];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec5*Eeta[ 89];

  H200 = H200 + G200;
  H201 = H201 + G201;
  H202 = H202 + G202;
  H210 = H210 + G210;
  H211 = H211 + G211;
  H212 = H212 + G212;
  H220 = H220 + G220;
  H221 = H221 + G221;
  H222 = H222 + G222;

  sigma[  2][  2] = two1nu * t2;

  t11 = H120 - H210;
  t12 = H120 - H210;
  sigma[  0][  0] = sigma[  0][  0] + (t11 + t12);
  t11 = H200 - H020;
  t12 = H121 - H211;
  sigma[  0][  1] = sigma[  0][  1] + (t11 + t12);
  t11 = H010 - H100;
  t12 = H122 - H212;
  sigma[  0][  2] = sigma[  0][  2] + (t11 + t12);
  t11 = H201 - H021;
  t12 = H201 - H021;
  sigma[  1][  1] = sigma[  1][  1] + (t11 + t12);
  t11 = H011 - H101;
  t12 = H202 - H022;
  sigma[  1][  2] = sigma[  1][  2] + (t11 + t12);
  t11 = H012 - H102;
  t12 = H012 - H102;
  sigma[  2][  2] = sigma[  2][  2] + (t11 + t12);

  t3 = 0.0;
  t3 = t3 + H120 - H210;
  t3 = t3 + H201 - H021;
  t3 = t3 + H012 - H102;
  t3 = t3 * two1nu;
sigma[  0][  0] = sigma[  0][  0] - t3;
sigma[  1][  1] = sigma[  1][  1] - t3;
sigma[  2][  2] = sigma[  2][  2] - t3;

  t =     5.00000000000000000000e-01 * mu8pi;
  sigma[  0][  0] = t * sigma[  0][  0];
  sigma[  0][  1] = t * sigma[  0][  1];
  sigma[  0][  2] = t * sigma[  0][  2];
  sigma[  1][  1] = t * sigma[  1][  1];
  sigma[  1][  2] = t * sigma[  1][  2];
  sigma[  2][  2] = t * sigma[  2][  2];
  sigma[1][0] = sigma[0][1];
  sigma[2][0] = sigma[0][2];
  sigma[2][1] = sigma[1][2];
}



void FMSigma2core3(real8 terms[],real8 mu8pi, real8 two1nu,real8 Eeta[],matrix sigma) {
real8 H000,H001,H002,H010,H011,H012,H020,H021,H022;
real8 H100,H101,H102,H110,H111,H112,H120,H121,H122;
real8 H200,H201,H202,H210,H211,H212,H220,H221,H222;
real8 G000,G001,G002,G010,G011,G012,G020,G021,G022;
real8 G100,G101,G102,G110,G111,G112,G120,G121,G122;
real8 G200,G201,G202,G210,G211,G212,G220,G221,G222;
real8 t,t11,t12,t2,t3;
real8 tvec0,tvec1,tvec2,tvec3,tvec4,tvec5,tvec6,tvec7,tvec8,tvec9;

  t2 = 0.0;
  tvec0 = terms[  0] *        1.0;
  tvec1 = terms[  1] *        3.0;
  tvec2 = terms[  2] *        3.0;
  tvec3 = terms[  3] *        1.0;
  tvec4 = terms[  7] *        3.0;
  tvec5 = terms[  8] *        6.0;
  tvec6 = terms[  9] *        3.0;
  tvec7 = terms[ 13] *        3.0;
  tvec8 = terms[ 14] *        3.0;
  tvec9 = terms[ 18] *        1.0;

  G000 = tvec0*Eeta[ 90];
  G001 = tvec0*Eeta[100];
  G002 = tvec0*Eeta[110];
  G010 = tvec0*Eeta[120];
  G011 = tvec0*Eeta[130];
  G012 = tvec0*Eeta[140];
  G020 = tvec0*Eeta[150];
  G021 = tvec0*Eeta[160];
  G022 = tvec0*Eeta[170];

  G000 = G000 + tvec1*Eeta[ 91];
  G001 = G001 + tvec1*Eeta[101];
  G002 = G002 + tvec1*Eeta[111];
  G010 = G010 + tvec1*Eeta[121];
  G011 = G011 + tvec1*Eeta[131];
  G012 = G012 + tvec1*Eeta[141];
  G020 = G020 + tvec1*Eeta[151];
  G021 = G021 + tvec1*Eeta[161];
  G022 = G022 + tvec1*Eeta[171];

  G000 = G000 + tvec2*Eeta[ 92];
  G001 = G001 + tvec2*Eeta[102];
  G002 = G002 + tvec2*Eeta[112];
  G010 = G010 + tvec2*Eeta[122];
  G011 = G011 + tvec2*Eeta[132];
  G012 = G012 + tvec2*Eeta[142];
  G020 = G020 + tvec2*Eeta[152];
  G021 = G021 + tvec2*Eeta[162];
  G022 = G022 + tvec2*Eeta[172];

  G000 = G000 + tvec3*Eeta[ 93];
  G001 = G001 + tvec3*Eeta[103];
  G002 = G002 + tvec3*Eeta[113];
  G010 = G010 + tvec3*Eeta[123];
  G011 = G011 + tvec3*Eeta[133];
  G012 = G012 + tvec3*Eeta[143];
  G020 = G020 + tvec3*Eeta[153];
  G021 = G021 + tvec3*Eeta[163];
  G022 = G022 + tvec3*Eeta[173];

  G000 = G000 + tvec4*Eeta[ 94];
  G001 = G001 + tvec4*Eeta[104];
  G002 = G002 + tvec4*Eeta[114];
  G010 = G010 + tvec4*Eeta[124];
  G011 = G011 + tvec4*Eeta[134];
  G012 = G012 + tvec4*Eeta[144];
  G020 = G020 + tvec4*Eeta[154];
  G021 = G021 + tvec4*Eeta[164];
  G022 = G022 + tvec4*Eeta[174];

  G000 = G000 + tvec5*Eeta[ 95];
  G001 = G001 + tvec5*Eeta[105];
  G002 = G002 + tvec5*Eeta[115];
  G010 = G010 + tvec5*Eeta[125];
  G011 = G011 + tvec5*Eeta[135];
  G012 = G012 + tvec5*Eeta[145];
  G020 = G020 + tvec5*Eeta[155];
  G021 = G021 + tvec5*Eeta[165];
  G022 = G022 + tvec5*Eeta[175];

  G000 = G000 + tvec6*Eeta[ 96];
  G001 = G001 + tvec6*Eeta[106];
  G002 = G002 + tvec6*Eeta[116];
  G010 = G010 + tvec6*Eeta[126];
  G011 = G011 + tvec6*Eeta[136];
  G012 = G012 + tvec6*Eeta[146];
  G020 = G020 + tvec6*Eeta[156];
  G021 = G021 + tvec6*Eeta[166];
  G022 = G022 + tvec6*Eeta[176];

  G000 = G000 + tvec7*Eeta[ 97];
  G001 = G001 + tvec7*Eeta[107];
  G002 = G002 + tvec7*Eeta[117];
  G010 = G010 + tvec7*Eeta[127];
  G011 = G011 + tvec7*Eeta[137];
  G012 = G012 + tvec7*Eeta[147];
  G020 = G020 + tvec7*Eeta[157];
  G021 = G021 + tvec7*Eeta[167];
  G022 = G022 + tvec7*Eeta[177];

  G000 = G000 + tvec8*Eeta[ 98];
  G001 = G001 + tvec8*Eeta[108];
  G002 = G002 + tvec8*Eeta[118];
  G010 = G010 + tvec8*Eeta[128];
  G011 = G011 + tvec8*Eeta[138];
  G012 = G012 + tvec8*Eeta[148];
  G020 = G020 + tvec8*Eeta[158];
  G021 = G021 + tvec8*Eeta[168];
  G022 = G022 + tvec8*Eeta[178];

  G000 = G000 + tvec9*Eeta[ 99];
  G001 = G001 + tvec9*Eeta[109];
  G002 = G002 + tvec9*Eeta[119];
  G010 = G010 + tvec9*Eeta[129];
  G011 = G011 + tvec9*Eeta[139];
  G012 = G012 + tvec9*Eeta[149];
  G020 = G020 + tvec9*Eeta[159];
  G021 = G021 + tvec9*Eeta[169];
  G022 = G022 + tvec9*Eeta[179];

  H000 = G000;
  H001 = G001;
  H002 = G002;
  H010 = G010;
  H011 = G011;
  H012 = G012;
  H020 = G020;
  H021 = G021;
  H022 = G022;

  tvec0 = terms[  1] *        1.0;
  tvec1 = terms[  2] *        3.0;
  tvec2 = terms[  3] *        3.0;
  tvec3 = terms[  4] *        1.0;
  tvec4 = terms[  8] *        3.0;
  tvec5 = terms[  9] *        6.0;
  tvec6 = terms[ 10] *        3.0;
  tvec7 = terms[ 14] *        3.0;
  tvec8 = terms[ 15] *        3.0;
  tvec9 = terms[ 19] *        1.0;

  G100 = tvec0*Eeta[ 90];
  G101 = tvec0*Eeta[100];
  G102 = tvec0*Eeta[110];
  G110 = tvec0*Eeta[120];
  G111 = tvec0*Eeta[130];
  G112 = tvec0*Eeta[140];
  G120 = tvec0*Eeta[150];
  G121 = tvec0*Eeta[160];
  G122 = tvec0*Eeta[170];

  G100 = G100 + tvec1*Eeta[ 91];
  G101 = G101 + tvec1*Eeta[101];
  G102 = G102 + tvec1*Eeta[111];
  G110 = G110 + tvec1*Eeta[121];
  G111 = G111 + tvec1*Eeta[131];
  G112 = G112 + tvec1*Eeta[141];
  G120 = G120 + tvec1*Eeta[151];
  G121 = G121 + tvec1*Eeta[161];
  G122 = G122 + tvec1*Eeta[171];

  G100 = G100 + tvec2*Eeta[ 92];
  G101 = G101 + tvec2*Eeta[102];
  G102 = G102 + tvec2*Eeta[112];
  G110 = G110 + tvec2*Eeta[122];
  G111 = G111 + tvec2*Eeta[132];
  G112 = G112 + tvec2*Eeta[142];
  G120 = G120 + tvec2*Eeta[152];
  G121 = G121 + tvec2*Eeta[162];
  G122 = G122 + tvec2*Eeta[172];

  G100 = G100 + tvec3*Eeta[ 93];
  G101 = G101 + tvec3*Eeta[103];
  G102 = G102 + tvec3*Eeta[113];
  G110 = G110 + tvec3*Eeta[123];
  G111 = G111 + tvec3*Eeta[133];
  G112 = G112 + tvec3*Eeta[143];
  G120 = G120 + tvec3*Eeta[153];
  G121 = G121 + tvec3*Eeta[163];
  G122 = G122 + tvec3*Eeta[173];

  G100 = G100 + tvec4*Eeta[ 94];
  G101 = G101 + tvec4*Eeta[104];
  G102 = G102 + tvec4*Eeta[114];
  G110 = G110 + tvec4*Eeta[124];
  G111 = G111 + tvec4*Eeta[134];
  G112 = G112 + tvec4*Eeta[144];
  G120 = G120 + tvec4*Eeta[154];
  G121 = G121 + tvec4*Eeta[164];
  G122 = G122 + tvec4*Eeta[174];

  G100 = G100 + tvec5*Eeta[ 95];
  G101 = G101 + tvec5*Eeta[105];
  G102 = G102 + tvec5*Eeta[115];
  G110 = G110 + tvec5*Eeta[125];
  G111 = G111 + tvec5*Eeta[135];
  G112 = G112 + tvec5*Eeta[145];
  G120 = G120 + tvec5*Eeta[155];
  G121 = G121 + tvec5*Eeta[165];
  G122 = G122 + tvec5*Eeta[175];

  G100 = G100 + tvec6*Eeta[ 96];
  G101 = G101 + tvec6*Eeta[106];
  G102 = G102 + tvec6*Eeta[116];
  G110 = G110 + tvec6*Eeta[126];
  G111 = G111 + tvec6*Eeta[136];
  G112 = G112 + tvec6*Eeta[146];
  G120 = G120 + tvec6*Eeta[156];
  G121 = G121 + tvec6*Eeta[166];
  G122 = G122 + tvec6*Eeta[176];

  G100 = G100 + tvec7*Eeta[ 97];
  G101 = G101 + tvec7*Eeta[107];
  G102 = G102 + tvec7*Eeta[117];
  G110 = G110 + tvec7*Eeta[127];
  G111 = G111 + tvec7*Eeta[137];
  G112 = G112 + tvec7*Eeta[147];
  G120 = G120 + tvec7*Eeta[157];
  G121 = G121 + tvec7*Eeta[167];
  G122 = G122 + tvec7*Eeta[177];

  G100 = G100 + tvec8*Eeta[ 98];
  G101 = G101 + tvec8*Eeta[108];
  G102 = G102 + tvec8*Eeta[118];
  G110 = G110 + tvec8*Eeta[128];
  G111 = G111 + tvec8*Eeta[138];
  G112 = G112 + tvec8*Eeta[148];
  G120 = G120 + tvec8*Eeta[158];
  G121 = G121 + tvec8*Eeta[168];
  G122 = G122 + tvec8*Eeta[178];

  G100 = G100 + tvec9*Eeta[ 99];
  G101 = G101 + tvec9*Eeta[109];
  G102 = G102 + tvec9*Eeta[119];
  G110 = G110 + tvec9*Eeta[129];
  G111 = G111 + tvec9*Eeta[139];
  G112 = G112 + tvec9*Eeta[149];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec9*Eeta[159];
  G121 = G121 + tvec9*Eeta[169];
  G122 = G122 + tvec9*Eeta[179];

  H100 = G100;
  H101 = G101;
  H102 = G102;
  H110 = G110;
  H111 = G111;
  H112 = G112;
  H120 = G120;
  H121 = G121;
  H122 = G122;

  tvec0 = terms[  7] *        1.0;
  tvec1 = terms[  8] *        3.0;
  tvec2 = terms[  9] *        3.0;
  tvec3 = terms[ 10] *        1.0;
  tvec4 = terms[ 13] *        3.0;
  tvec5 = terms[ 14] *        6.0;
  tvec6 = terms[ 15] *        3.0;
  tvec7 = terms[ 18] *        3.0;
  tvec8 = terms[ 19] *        3.0;
  tvec9 = terms[ 22] *        1.0;

  G200 = tvec0*Eeta[ 90];
  G201 = tvec0*Eeta[100];
  G202 = tvec0*Eeta[110];
  G210 = tvec0*Eeta[120];
  G211 = tvec0*Eeta[130];
  G212 = tvec0*Eeta[140];
  G220 = tvec0*Eeta[150];
  G221 = tvec0*Eeta[160];
  G222 = tvec0*Eeta[170];

  G200 = G200 + tvec1*Eeta[ 91];
  G201 = G201 + tvec1*Eeta[101];
  G202 = G202 + tvec1*Eeta[111];
  G210 = G210 + tvec1*Eeta[121];
  G211 = G211 + tvec1*Eeta[131];
  G212 = G212 + tvec1*Eeta[141];
  G220 = G220 + tvec1*Eeta[151];
  G221 = G221 + tvec1*Eeta[161];
  G222 = G222 + tvec1*Eeta[171];

  G200 = G200 + tvec2*Eeta[ 92];
  G201 = G201 + tvec2*Eeta[102];
  G202 = G202 + tvec2*Eeta[112];
  G210 = G210 + tvec2*Eeta[122];
  G211 = G211 + tvec2*Eeta[132];
  G212 = G212 + tvec2*Eeta[142];
  G220 = G220 + tvec2*Eeta[152];
  G221 = G221 + tvec2*Eeta[162];
  G222 = G222 + tvec2*Eeta[172];

  G200 = G200 + tvec3*Eeta[ 93];
  G201 = G201 + tvec3*Eeta[103];
  G202 = G202 + tvec3*Eeta[113];
  G210 = G210 + tvec3*Eeta[123];
  G211 = G211 + tvec3*Eeta[133];
  G212 = G212 + tvec3*Eeta[143];
  G220 = G220 + tvec3*Eeta[153];
  G221 = G221 + tvec3*Eeta[163];
  G222 = G222 + tvec3*Eeta[173];

  G200 = G200 + tvec4*Eeta[ 94];
  G201 = G201 + tvec4*Eeta[104];
  G202 = G202 + tvec4*Eeta[114];
  G210 = G210 + tvec4*Eeta[124];
  G211 = G211 + tvec4*Eeta[134];
  G212 = G212 + tvec4*Eeta[144];
  G220 = G220 + tvec4*Eeta[154];
  G221 = G221 + tvec4*Eeta[164];
  G222 = G222 + tvec4*Eeta[174];

  G200 = G200 + tvec5*Eeta[ 95];
  G201 = G201 + tvec5*Eeta[105];
  G202 = G202 + tvec5*Eeta[115];
  G210 = G210 + tvec5*Eeta[125];
  G211 = G211 + tvec5*Eeta[135];
  G212 = G212 + tvec5*Eeta[145];
  G220 = G220 + tvec5*Eeta[155];
  G221 = G221 + tvec5*Eeta[165];
  G222 = G222 + tvec5*Eeta[175];

  G200 = G200 + tvec6*Eeta[ 96];
  G201 = G201 + tvec6*Eeta[106];
  G202 = G202 + tvec6*Eeta[116];
  G210 = G210 + tvec6*Eeta[126];
  G211 = G211 + tvec6*Eeta[136];
  G212 = G212 + tvec6*Eeta[146];
  G220 = G220 + tvec6*Eeta[156];
  G221 = G221 + tvec6*Eeta[166];
  G222 = G222 + tvec6*Eeta[176];

  G200 = G200 + tvec7*Eeta[ 97];
  G201 = G201 + tvec7*Eeta[107];
  G202 = G202 + tvec7*Eeta[117];
  G210 = G210 + tvec7*Eeta[127];
  G211 = G211 + tvec7*Eeta[137];
  G212 = G212 + tvec7*Eeta[147];
  G220 = G220 + tvec7*Eeta[157];
  G221 = G221 + tvec7*Eeta[167];
  G222 = G222 + tvec7*Eeta[177];

  G200 = G200 + tvec8*Eeta[ 98];
  G201 = G201 + tvec8*Eeta[108];
  G202 = G202 + tvec8*Eeta[118];
  G210 = G210 + tvec8*Eeta[128];
  G211 = G211 + tvec8*Eeta[138];
  G212 = G212 + tvec8*Eeta[148];
  G220 = G220 + tvec8*Eeta[158];
  G221 = G221 + tvec8*Eeta[168];
  G222 = G222 + tvec8*Eeta[178];

  G200 = G200 + tvec9*Eeta[ 99];
  G201 = G201 + tvec9*Eeta[109];
  G202 = G202 + tvec9*Eeta[119];
  G210 = G210 + tvec9*Eeta[129];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec9*Eeta[139];
  G212 = G212 + tvec9*Eeta[149];
  G220 = G220 + tvec9*Eeta[159];
  G221 = G221 + tvec9*Eeta[169];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec9*Eeta[179];

  H200 = G200;
  H201 = G201;
  H202 = G202;
  H210 = G210;
  H211 = G211;
  H212 = G212;
  H220 = G220;
  H221 = G221;
  H222 = G222;

  sigma[  0][  0] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  1] *        1.0;
  tvec1 = terms[  2] *        3.0;
  tvec2 = terms[  3] *        3.0;
  tvec3 = terms[  4] *        1.0;
  tvec4 = terms[  8] *        3.0;
  tvec5 = terms[  9] *        6.0;
  tvec6 = terms[ 10] *        3.0;
  tvec7 = terms[ 14] *        3.0;
  tvec8 = terms[ 15] *        3.0;
  tvec9 = terms[ 19] *        1.0;

  G000 = tvec0*Eeta[ 90];
  G001 = tvec0*Eeta[100];
  G002 = tvec0*Eeta[110];
  G010 = tvec0*Eeta[120];
  G011 = tvec0*Eeta[130];
  G012 = tvec0*Eeta[140];
  G020 = tvec0*Eeta[150];
  G021 = tvec0*Eeta[160];
  G022 = tvec0*Eeta[170];

  G000 = G000 + tvec1*Eeta[ 91];
  G001 = G001 + tvec1*Eeta[101];
  G002 = G002 + tvec1*Eeta[111];
  G010 = G010 + tvec1*Eeta[121];
  G011 = G011 + tvec1*Eeta[131];
  G012 = G012 + tvec1*Eeta[141];
  G020 = G020 + tvec1*Eeta[151];
  G021 = G021 + tvec1*Eeta[161];
  G022 = G022 + tvec1*Eeta[171];

  G000 = G000 + tvec2*Eeta[ 92];
  G001 = G001 + tvec2*Eeta[102];
  G002 = G002 + tvec2*Eeta[112];
  G010 = G010 + tvec2*Eeta[122];
  G011 = G011 + tvec2*Eeta[132];
  G012 = G012 + tvec2*Eeta[142];
  G020 = G020 + tvec2*Eeta[152];
  G021 = G021 + tvec2*Eeta[162];
  G022 = G022 + tvec2*Eeta[172];

  G000 = G000 + tvec3*Eeta[ 93];
  G001 = G001 + tvec3*Eeta[103];
  G002 = G002 + tvec3*Eeta[113];
  G010 = G010 + tvec3*Eeta[123];
  G011 = G011 + tvec3*Eeta[133];
  G012 = G012 + tvec3*Eeta[143];
  G020 = G020 + tvec3*Eeta[153];
  G021 = G021 + tvec3*Eeta[163];
  G022 = G022 + tvec3*Eeta[173];

  G000 = G000 + tvec4*Eeta[ 94];
  G001 = G001 + tvec4*Eeta[104];
  G002 = G002 + tvec4*Eeta[114];
  G010 = G010 + tvec4*Eeta[124];
  G011 = G011 + tvec4*Eeta[134];
  G012 = G012 + tvec4*Eeta[144];
  G020 = G020 + tvec4*Eeta[154];
  G021 = G021 + tvec4*Eeta[164];
  G022 = G022 + tvec4*Eeta[174];

  G000 = G000 + tvec5*Eeta[ 95];
  G001 = G001 + tvec5*Eeta[105];
  G002 = G002 + tvec5*Eeta[115];
  G010 = G010 + tvec5*Eeta[125];
  G011 = G011 + tvec5*Eeta[135];
  G012 = G012 + tvec5*Eeta[145];
  G020 = G020 + tvec5*Eeta[155];
  G021 = G021 + tvec5*Eeta[165];
  G022 = G022 + tvec5*Eeta[175];

  G000 = G000 + tvec6*Eeta[ 96];
  G001 = G001 + tvec6*Eeta[106];
  G002 = G002 + tvec6*Eeta[116];
  G010 = G010 + tvec6*Eeta[126];
  G011 = G011 + tvec6*Eeta[136];
  G012 = G012 + tvec6*Eeta[146];
  G020 = G020 + tvec6*Eeta[156];
  G021 = G021 + tvec6*Eeta[166];
  G022 = G022 + tvec6*Eeta[176];

  G000 = G000 + tvec7*Eeta[ 97];
  G001 = G001 + tvec7*Eeta[107];
  G002 = G002 + tvec7*Eeta[117];
  G010 = G010 + tvec7*Eeta[127];
  G011 = G011 + tvec7*Eeta[137];
  G012 = G012 + tvec7*Eeta[147];
  G020 = G020 + tvec7*Eeta[157];
  G021 = G021 + tvec7*Eeta[167];
  G022 = G022 + tvec7*Eeta[177];

  G000 = G000 + tvec8*Eeta[ 98];
  G001 = G001 + tvec8*Eeta[108];
  G002 = G002 + tvec8*Eeta[118];
  G010 = G010 + tvec8*Eeta[128];
  G011 = G011 + tvec8*Eeta[138];
  G012 = G012 + tvec8*Eeta[148];
  G020 = G020 + tvec8*Eeta[158];
  G021 = G021 + tvec8*Eeta[168];
  G022 = G022 + tvec8*Eeta[178];

  G000 = G000 + tvec9*Eeta[ 99];
  G001 = G001 + tvec9*Eeta[109];
  G002 = G002 + tvec9*Eeta[119];
  G010 = G010 + tvec9*Eeta[129];
  G011 = G011 + tvec9*Eeta[139];
  G012 = G012 + tvec9*Eeta[149];
  G020 = G020 + tvec9*Eeta[159];
  G021 = G021 + tvec9*Eeta[169];
  G022 = G022 + tvec9*Eeta[179];

  tvec0 = terms[  2] *        1.0;
  tvec1 = terms[  3] *        3.0;
  tvec2 = terms[  4] *        3.0;
  tvec3 = terms[  5] *        1.0;
  tvec4 = terms[  9] *        3.0;
  tvec5 = terms[ 10] *        6.0;
  tvec6 = terms[ 11] *        3.0;
  tvec7 = terms[ 15] *        3.0;
  tvec8 = terms[ 16] *        3.0;
  tvec9 = terms[ 20] *        1.0;

  G100 = tvec0*Eeta[ 90];
  G101 = tvec0*Eeta[100];
  G102 = tvec0*Eeta[110];
  G110 = tvec0*Eeta[120];
  G111 = tvec0*Eeta[130];
  G112 = tvec0*Eeta[140];
  G120 = tvec0*Eeta[150];
  G121 = tvec0*Eeta[160];
  G122 = tvec0*Eeta[170];

  G100 = G100 + tvec1*Eeta[ 91];
  G101 = G101 + tvec1*Eeta[101];
  G102 = G102 + tvec1*Eeta[111];
  G110 = G110 + tvec1*Eeta[121];
  G111 = G111 + tvec1*Eeta[131];
  G112 = G112 + tvec1*Eeta[141];
  G120 = G120 + tvec1*Eeta[151];
  G121 = G121 + tvec1*Eeta[161];
  G122 = G122 + tvec1*Eeta[171];

  G100 = G100 + tvec2*Eeta[ 92];
  G101 = G101 + tvec2*Eeta[102];
  G102 = G102 + tvec2*Eeta[112];
  G110 = G110 + tvec2*Eeta[122];
  G111 = G111 + tvec2*Eeta[132];
  G112 = G112 + tvec2*Eeta[142];
  G120 = G120 + tvec2*Eeta[152];
  G121 = G121 + tvec2*Eeta[162];
  G122 = G122 + tvec2*Eeta[172];

  G100 = G100 + tvec3*Eeta[ 93];
  G101 = G101 + tvec3*Eeta[103];
  G102 = G102 + tvec3*Eeta[113];
  G110 = G110 + tvec3*Eeta[123];
  G111 = G111 + tvec3*Eeta[133];
  G112 = G112 + tvec3*Eeta[143];
  G120 = G120 + tvec3*Eeta[153];
  G121 = G121 + tvec3*Eeta[163];
  G122 = G122 + tvec3*Eeta[173];

  G100 = G100 + tvec4*Eeta[ 94];
  G101 = G101 + tvec4*Eeta[104];
  G102 = G102 + tvec4*Eeta[114];
  G110 = G110 + tvec4*Eeta[124];
  G111 = G111 + tvec4*Eeta[134];
  G112 = G112 + tvec4*Eeta[144];
  G120 = G120 + tvec4*Eeta[154];
  G121 = G121 + tvec4*Eeta[164];
  G122 = G122 + tvec4*Eeta[174];

  G100 = G100 + tvec5*Eeta[ 95];
  G101 = G101 + tvec5*Eeta[105];
  G102 = G102 + tvec5*Eeta[115];
  G110 = G110 + tvec5*Eeta[125];
  G111 = G111 + tvec5*Eeta[135];
  G112 = G112 + tvec5*Eeta[145];
  G120 = G120 + tvec5*Eeta[155];
  G121 = G121 + tvec5*Eeta[165];
  G122 = G122 + tvec5*Eeta[175];

  G100 = G100 + tvec6*Eeta[ 96];
  G101 = G101 + tvec6*Eeta[106];
  G102 = G102 + tvec6*Eeta[116];
  G110 = G110 + tvec6*Eeta[126];
  G111 = G111 + tvec6*Eeta[136];
  G112 = G112 + tvec6*Eeta[146];
  G120 = G120 + tvec6*Eeta[156];
  G121 = G121 + tvec6*Eeta[166];
  G122 = G122 + tvec6*Eeta[176];

  G100 = G100 + tvec7*Eeta[ 97];
  G101 = G101 + tvec7*Eeta[107];
  G102 = G102 + tvec7*Eeta[117];
  G110 = G110 + tvec7*Eeta[127];
  G111 = G111 + tvec7*Eeta[137];
  G112 = G112 + tvec7*Eeta[147];
  G120 = G120 + tvec7*Eeta[157];
  G121 = G121 + tvec7*Eeta[167];
  G122 = G122 + tvec7*Eeta[177];

  G100 = G100 + tvec8*Eeta[ 98];
  G101 = G101 + tvec8*Eeta[108];
  G102 = G102 + tvec8*Eeta[118];
  G110 = G110 + tvec8*Eeta[128];
  G111 = G111 + tvec8*Eeta[138];
  G112 = G112 + tvec8*Eeta[148];
  G120 = G120 + tvec8*Eeta[158];
  G121 = G121 + tvec8*Eeta[168];
  G122 = G122 + tvec8*Eeta[178];

  G100 = G100 + tvec9*Eeta[ 99];
  G101 = G101 + tvec9*Eeta[109];
  G102 = G102 + tvec9*Eeta[119];
  G110 = G110 + tvec9*Eeta[129];
  G111 = G111 + tvec9*Eeta[139];
  G112 = G112 + tvec9*Eeta[149];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec9*Eeta[159];
  G121 = G121 + tvec9*Eeta[169];
  G122 = G122 + tvec9*Eeta[179];

  tvec0 = terms[  8] *        1.0;
  tvec1 = terms[  9] *        3.0;
  tvec2 = terms[ 10] *        3.0;
  tvec3 = terms[ 11] *        1.0;
  tvec4 = terms[ 14] *        3.0;
  tvec5 = terms[ 15] *        6.0;
  tvec6 = terms[ 16] *        3.0;
  tvec7 = terms[ 19] *        3.0;
  tvec8 = terms[ 20] *        3.0;
  tvec9 = terms[ 23] *        1.0;

  G200 = tvec0*Eeta[ 90];
  G201 = tvec0*Eeta[100];
  G202 = tvec0*Eeta[110];
  G210 = tvec0*Eeta[120];
  G211 = tvec0*Eeta[130];
  G212 = tvec0*Eeta[140];
  G220 = tvec0*Eeta[150];
  G221 = tvec0*Eeta[160];
  G222 = tvec0*Eeta[170];

  G200 = G200 + tvec1*Eeta[ 91];
  G201 = G201 + tvec1*Eeta[101];
  G202 = G202 + tvec1*Eeta[111];
  G210 = G210 + tvec1*Eeta[121];
  G211 = G211 + tvec1*Eeta[131];
  G212 = G212 + tvec1*Eeta[141];
  G220 = G220 + tvec1*Eeta[151];
  G221 = G221 + tvec1*Eeta[161];
  G222 = G222 + tvec1*Eeta[171];

  G200 = G200 + tvec2*Eeta[ 92];
  G201 = G201 + tvec2*Eeta[102];
  G202 = G202 + tvec2*Eeta[112];
  G210 = G210 + tvec2*Eeta[122];
  G211 = G211 + tvec2*Eeta[132];
  G212 = G212 + tvec2*Eeta[142];
  G220 = G220 + tvec2*Eeta[152];
  G221 = G221 + tvec2*Eeta[162];
  G222 = G222 + tvec2*Eeta[172];

  G200 = G200 + tvec3*Eeta[ 93];
  G201 = G201 + tvec3*Eeta[103];
  G202 = G202 + tvec3*Eeta[113];
  G210 = G210 + tvec3*Eeta[123];
  G211 = G211 + tvec3*Eeta[133];
  G212 = G212 + tvec3*Eeta[143];
  G220 = G220 + tvec3*Eeta[153];
  G221 = G221 + tvec3*Eeta[163];
  G222 = G222 + tvec3*Eeta[173];

  G200 = G200 + tvec4*Eeta[ 94];
  G201 = G201 + tvec4*Eeta[104];
  G202 = G202 + tvec4*Eeta[114];
  G210 = G210 + tvec4*Eeta[124];
  G211 = G211 + tvec4*Eeta[134];
  G212 = G212 + tvec4*Eeta[144];
  G220 = G220 + tvec4*Eeta[154];
  G221 = G221 + tvec4*Eeta[164];
  G222 = G222 + tvec4*Eeta[174];

  G200 = G200 + tvec5*Eeta[ 95];
  G201 = G201 + tvec5*Eeta[105];
  G202 = G202 + tvec5*Eeta[115];
  G210 = G210 + tvec5*Eeta[125];
  G211 = G211 + tvec5*Eeta[135];
  G212 = G212 + tvec5*Eeta[145];
  G220 = G220 + tvec5*Eeta[155];
  G221 = G221 + tvec5*Eeta[165];
  G222 = G222 + tvec5*Eeta[175];

  G200 = G200 + tvec6*Eeta[ 96];
  G201 = G201 + tvec6*Eeta[106];
  G202 = G202 + tvec6*Eeta[116];
  G210 = G210 + tvec6*Eeta[126];
  G211 = G211 + tvec6*Eeta[136];
  G212 = G212 + tvec6*Eeta[146];
  G220 = G220 + tvec6*Eeta[156];
  G221 = G221 + tvec6*Eeta[166];
  G222 = G222 + tvec6*Eeta[176];

  G200 = G200 + tvec7*Eeta[ 97];
  G201 = G201 + tvec7*Eeta[107];
  G202 = G202 + tvec7*Eeta[117];
  G210 = G210 + tvec7*Eeta[127];
  G211 = G211 + tvec7*Eeta[137];
  G212 = G212 + tvec7*Eeta[147];
  G220 = G220 + tvec7*Eeta[157];
  G221 = G221 + tvec7*Eeta[167];
  G222 = G222 + tvec7*Eeta[177];

  G200 = G200 + tvec8*Eeta[ 98];
  G201 = G201 + tvec8*Eeta[108];
  G202 = G202 + tvec8*Eeta[118];
  G210 = G210 + tvec8*Eeta[128];
  G211 = G211 + tvec8*Eeta[138];
  G212 = G212 + tvec8*Eeta[148];
  G220 = G220 + tvec8*Eeta[158];
  G221 = G221 + tvec8*Eeta[168];
  G222 = G222 + tvec8*Eeta[178];

  G200 = G200 + tvec9*Eeta[ 99];
  G201 = G201 + tvec9*Eeta[109];
  G202 = G202 + tvec9*Eeta[119];
  G210 = G210 + tvec9*Eeta[129];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec9*Eeta[139];
  G212 = G212 + tvec9*Eeta[149];
  G220 = G220 + tvec9*Eeta[159];
  G221 = G221 + tvec9*Eeta[169];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec9*Eeta[179];

  sigma[  0][  1] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  7] *        1.0;
  tvec1 = terms[  8] *        3.0;
  tvec2 = terms[  9] *        3.0;
  tvec3 = terms[ 10] *        1.0;
  tvec4 = terms[ 13] *        3.0;
  tvec5 = terms[ 14] *        6.0;
  tvec6 = terms[ 15] *        3.0;
  tvec7 = terms[ 18] *        3.0;
  tvec8 = terms[ 19] *        3.0;
  tvec9 = terms[ 22] *        1.0;

  G000 = tvec0*Eeta[ 90];
  G001 = tvec0*Eeta[100];
  G002 = tvec0*Eeta[110];
  G010 = tvec0*Eeta[120];
  G011 = tvec0*Eeta[130];
  G012 = tvec0*Eeta[140];
  G020 = tvec0*Eeta[150];
  G021 = tvec0*Eeta[160];
  G022 = tvec0*Eeta[170];

  G000 = G000 + tvec1*Eeta[ 91];
  G001 = G001 + tvec1*Eeta[101];
  G002 = G002 + tvec1*Eeta[111];
  G010 = G010 + tvec1*Eeta[121];
  G011 = G011 + tvec1*Eeta[131];
  G012 = G012 + tvec1*Eeta[141];
  G020 = G020 + tvec1*Eeta[151];
  G021 = G021 + tvec1*Eeta[161];
  G022 = G022 + tvec1*Eeta[171];

  G000 = G000 + tvec2*Eeta[ 92];
  G001 = G001 + tvec2*Eeta[102];
  G002 = G002 + tvec2*Eeta[112];
  G010 = G010 + tvec2*Eeta[122];
  G011 = G011 + tvec2*Eeta[132];
  G012 = G012 + tvec2*Eeta[142];
  G020 = G020 + tvec2*Eeta[152];
  G021 = G021 + tvec2*Eeta[162];
  G022 = G022 + tvec2*Eeta[172];

  G000 = G000 + tvec3*Eeta[ 93];
  G001 = G001 + tvec3*Eeta[103];
  G002 = G002 + tvec3*Eeta[113];
  G010 = G010 + tvec3*Eeta[123];
  G011 = G011 + tvec3*Eeta[133];
  G012 = G012 + tvec3*Eeta[143];
  G020 = G020 + tvec3*Eeta[153];
  G021 = G021 + tvec3*Eeta[163];
  G022 = G022 + tvec3*Eeta[173];

  G000 = G000 + tvec4*Eeta[ 94];
  G001 = G001 + tvec4*Eeta[104];
  G002 = G002 + tvec4*Eeta[114];
  G010 = G010 + tvec4*Eeta[124];
  G011 = G011 + tvec4*Eeta[134];
  G012 = G012 + tvec4*Eeta[144];
  G020 = G020 + tvec4*Eeta[154];
  G021 = G021 + tvec4*Eeta[164];
  G022 = G022 + tvec4*Eeta[174];

  G000 = G000 + tvec5*Eeta[ 95];
  G001 = G001 + tvec5*Eeta[105];
  G002 = G002 + tvec5*Eeta[115];
  G010 = G010 + tvec5*Eeta[125];
  G011 = G011 + tvec5*Eeta[135];
  G012 = G012 + tvec5*Eeta[145];
  G020 = G020 + tvec5*Eeta[155];
  G021 = G021 + tvec5*Eeta[165];
  G022 = G022 + tvec5*Eeta[175];

  G000 = G000 + tvec6*Eeta[ 96];
  G001 = G001 + tvec6*Eeta[106];
  G002 = G002 + tvec6*Eeta[116];
  G010 = G010 + tvec6*Eeta[126];
  G011 = G011 + tvec6*Eeta[136];
  G012 = G012 + tvec6*Eeta[146];
  G020 = G020 + tvec6*Eeta[156];
  G021 = G021 + tvec6*Eeta[166];
  G022 = G022 + tvec6*Eeta[176];

  G000 = G000 + tvec7*Eeta[ 97];
  G001 = G001 + tvec7*Eeta[107];
  G002 = G002 + tvec7*Eeta[117];
  G010 = G010 + tvec7*Eeta[127];
  G011 = G011 + tvec7*Eeta[137];
  G012 = G012 + tvec7*Eeta[147];
  G020 = G020 + tvec7*Eeta[157];
  G021 = G021 + tvec7*Eeta[167];
  G022 = G022 + tvec7*Eeta[177];

  G000 = G000 + tvec8*Eeta[ 98];
  G001 = G001 + tvec8*Eeta[108];
  G002 = G002 + tvec8*Eeta[118];
  G010 = G010 + tvec8*Eeta[128];
  G011 = G011 + tvec8*Eeta[138];
  G012 = G012 + tvec8*Eeta[148];
  G020 = G020 + tvec8*Eeta[158];
  G021 = G021 + tvec8*Eeta[168];
  G022 = G022 + tvec8*Eeta[178];

  G000 = G000 + tvec9*Eeta[ 99];
  G001 = G001 + tvec9*Eeta[109];
  G002 = G002 + tvec9*Eeta[119];
  G010 = G010 + tvec9*Eeta[129];
  G011 = G011 + tvec9*Eeta[139];
  G012 = G012 + tvec9*Eeta[149];
  G020 = G020 + tvec9*Eeta[159];
  G021 = G021 + tvec9*Eeta[169];
  G022 = G022 + tvec9*Eeta[179];

  tvec0 = terms[  8] *        1.0;
  tvec1 = terms[  9] *        3.0;
  tvec2 = terms[ 10] *        3.0;
  tvec3 = terms[ 11] *        1.0;
  tvec4 = terms[ 14] *        3.0;
  tvec5 = terms[ 15] *        6.0;
  tvec6 = terms[ 16] *        3.0;
  tvec7 = terms[ 19] *        3.0;
  tvec8 = terms[ 20] *        3.0;
  tvec9 = terms[ 23] *        1.0;

  G100 = tvec0*Eeta[ 90];
  G101 = tvec0*Eeta[100];
  G102 = tvec0*Eeta[110];
  G110 = tvec0*Eeta[120];
  G111 = tvec0*Eeta[130];
  G112 = tvec0*Eeta[140];
  G120 = tvec0*Eeta[150];
  G121 = tvec0*Eeta[160];
  G122 = tvec0*Eeta[170];

  G100 = G100 + tvec1*Eeta[ 91];
  G101 = G101 + tvec1*Eeta[101];
  G102 = G102 + tvec1*Eeta[111];
  G110 = G110 + tvec1*Eeta[121];
  G111 = G111 + tvec1*Eeta[131];
  G112 = G112 + tvec1*Eeta[141];
  G120 = G120 + tvec1*Eeta[151];
  G121 = G121 + tvec1*Eeta[161];
  G122 = G122 + tvec1*Eeta[171];

  G100 = G100 + tvec2*Eeta[ 92];
  G101 = G101 + tvec2*Eeta[102];
  G102 = G102 + tvec2*Eeta[112];
  G110 = G110 + tvec2*Eeta[122];
  G111 = G111 + tvec2*Eeta[132];
  G112 = G112 + tvec2*Eeta[142];
  G120 = G120 + tvec2*Eeta[152];
  G121 = G121 + tvec2*Eeta[162];
  G122 = G122 + tvec2*Eeta[172];

  G100 = G100 + tvec3*Eeta[ 93];
  G101 = G101 + tvec3*Eeta[103];
  G102 = G102 + tvec3*Eeta[113];
  G110 = G110 + tvec3*Eeta[123];
  G111 = G111 + tvec3*Eeta[133];
  G112 = G112 + tvec3*Eeta[143];
  G120 = G120 + tvec3*Eeta[153];
  G121 = G121 + tvec3*Eeta[163];
  G122 = G122 + tvec3*Eeta[173];

  G100 = G100 + tvec4*Eeta[ 94];
  G101 = G101 + tvec4*Eeta[104];
  G102 = G102 + tvec4*Eeta[114];
  G110 = G110 + tvec4*Eeta[124];
  G111 = G111 + tvec4*Eeta[134];
  G112 = G112 + tvec4*Eeta[144];
  G120 = G120 + tvec4*Eeta[154];
  G121 = G121 + tvec4*Eeta[164];
  G122 = G122 + tvec4*Eeta[174];

  G100 = G100 + tvec5*Eeta[ 95];
  G101 = G101 + tvec5*Eeta[105];
  G102 = G102 + tvec5*Eeta[115];
  G110 = G110 + tvec5*Eeta[125];
  G111 = G111 + tvec5*Eeta[135];
  G112 = G112 + tvec5*Eeta[145];
  G120 = G120 + tvec5*Eeta[155];
  G121 = G121 + tvec5*Eeta[165];
  G122 = G122 + tvec5*Eeta[175];

  G100 = G100 + tvec6*Eeta[ 96];
  G101 = G101 + tvec6*Eeta[106];
  G102 = G102 + tvec6*Eeta[116];
  G110 = G110 + tvec6*Eeta[126];
  G111 = G111 + tvec6*Eeta[136];
  G112 = G112 + tvec6*Eeta[146];
  G120 = G120 + tvec6*Eeta[156];
  G121 = G121 + tvec6*Eeta[166];
  G122 = G122 + tvec6*Eeta[176];

  G100 = G100 + tvec7*Eeta[ 97];
  G101 = G101 + tvec7*Eeta[107];
  G102 = G102 + tvec7*Eeta[117];
  G110 = G110 + tvec7*Eeta[127];
  G111 = G111 + tvec7*Eeta[137];
  G112 = G112 + tvec7*Eeta[147];
  G120 = G120 + tvec7*Eeta[157];
  G121 = G121 + tvec7*Eeta[167];
  G122 = G122 + tvec7*Eeta[177];

  G100 = G100 + tvec8*Eeta[ 98];
  G101 = G101 + tvec8*Eeta[108];
  G102 = G102 + tvec8*Eeta[118];
  G110 = G110 + tvec8*Eeta[128];
  G111 = G111 + tvec8*Eeta[138];
  G112 = G112 + tvec8*Eeta[148];
  G120 = G120 + tvec8*Eeta[158];
  G121 = G121 + tvec8*Eeta[168];
  G122 = G122 + tvec8*Eeta[178];

  G100 = G100 + tvec9*Eeta[ 99];
  G101 = G101 + tvec9*Eeta[109];
  G102 = G102 + tvec9*Eeta[119];
  G110 = G110 + tvec9*Eeta[129];
  G111 = G111 + tvec9*Eeta[139];
  G112 = G112 + tvec9*Eeta[149];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec9*Eeta[159];
  G121 = G121 + tvec9*Eeta[169];
  G122 = G122 + tvec9*Eeta[179];

  tvec0 = terms[ 13] *        1.0;
  tvec1 = terms[ 14] *        3.0;
  tvec2 = terms[ 15] *        3.0;
  tvec3 = terms[ 16] *        1.0;
  tvec4 = terms[ 18] *        3.0;
  tvec5 = terms[ 19] *        6.0;
  tvec6 = terms[ 20] *        3.0;
  tvec7 = terms[ 22] *        3.0;
  tvec8 = terms[ 23] *        3.0;
  tvec9 = terms[ 25] *        1.0;

  G200 = tvec0*Eeta[ 90];
  G201 = tvec0*Eeta[100];
  G202 = tvec0*Eeta[110];
  G210 = tvec0*Eeta[120];
  G211 = tvec0*Eeta[130];
  G212 = tvec0*Eeta[140];
  G220 = tvec0*Eeta[150];
  G221 = tvec0*Eeta[160];
  G222 = tvec0*Eeta[170];

  G200 = G200 + tvec1*Eeta[ 91];
  G201 = G201 + tvec1*Eeta[101];
  G202 = G202 + tvec1*Eeta[111];
  G210 = G210 + tvec1*Eeta[121];
  G211 = G211 + tvec1*Eeta[131];
  G212 = G212 + tvec1*Eeta[141];
  G220 = G220 + tvec1*Eeta[151];
  G221 = G221 + tvec1*Eeta[161];
  G222 = G222 + tvec1*Eeta[171];

  G200 = G200 + tvec2*Eeta[ 92];
  G201 = G201 + tvec2*Eeta[102];
  G202 = G202 + tvec2*Eeta[112];
  G210 = G210 + tvec2*Eeta[122];
  G211 = G211 + tvec2*Eeta[132];
  G212 = G212 + tvec2*Eeta[142];
  G220 = G220 + tvec2*Eeta[152];
  G221 = G221 + tvec2*Eeta[162];
  G222 = G222 + tvec2*Eeta[172];

  G200 = G200 + tvec3*Eeta[ 93];
  G201 = G201 + tvec3*Eeta[103];
  G202 = G202 + tvec3*Eeta[113];
  G210 = G210 + tvec3*Eeta[123];
  G211 = G211 + tvec3*Eeta[133];
  G212 = G212 + tvec3*Eeta[143];
  G220 = G220 + tvec3*Eeta[153];
  G221 = G221 + tvec3*Eeta[163];
  G222 = G222 + tvec3*Eeta[173];

  G200 = G200 + tvec4*Eeta[ 94];
  G201 = G201 + tvec4*Eeta[104];
  G202 = G202 + tvec4*Eeta[114];
  G210 = G210 + tvec4*Eeta[124];
  G211 = G211 + tvec4*Eeta[134];
  G212 = G212 + tvec4*Eeta[144];
  G220 = G220 + tvec4*Eeta[154];
  G221 = G221 + tvec4*Eeta[164];
  G222 = G222 + tvec4*Eeta[174];

  G200 = G200 + tvec5*Eeta[ 95];
  G201 = G201 + tvec5*Eeta[105];
  G202 = G202 + tvec5*Eeta[115];
  G210 = G210 + tvec5*Eeta[125];
  G211 = G211 + tvec5*Eeta[135];
  G212 = G212 + tvec5*Eeta[145];
  G220 = G220 + tvec5*Eeta[155];
  G221 = G221 + tvec5*Eeta[165];
  G222 = G222 + tvec5*Eeta[175];

  G200 = G200 + tvec6*Eeta[ 96];
  G201 = G201 + tvec6*Eeta[106];
  G202 = G202 + tvec6*Eeta[116];
  G210 = G210 + tvec6*Eeta[126];
  G211 = G211 + tvec6*Eeta[136];
  G212 = G212 + tvec6*Eeta[146];
  G220 = G220 + tvec6*Eeta[156];
  G221 = G221 + tvec6*Eeta[166];
  G222 = G222 + tvec6*Eeta[176];

  G200 = G200 + tvec7*Eeta[ 97];
  G201 = G201 + tvec7*Eeta[107];
  G202 = G202 + tvec7*Eeta[117];
  G210 = G210 + tvec7*Eeta[127];
  G211 = G211 + tvec7*Eeta[137];
  G212 = G212 + tvec7*Eeta[147];
  G220 = G220 + tvec7*Eeta[157];
  G221 = G221 + tvec7*Eeta[167];
  G222 = G222 + tvec7*Eeta[177];

  G200 = G200 + tvec8*Eeta[ 98];
  G201 = G201 + tvec8*Eeta[108];
  G202 = G202 + tvec8*Eeta[118];
  G210 = G210 + tvec8*Eeta[128];
  G211 = G211 + tvec8*Eeta[138];
  G212 = G212 + tvec8*Eeta[148];
  G220 = G220 + tvec8*Eeta[158];
  G221 = G221 + tvec8*Eeta[168];
  G222 = G222 + tvec8*Eeta[178];

  G200 = G200 + tvec9*Eeta[ 99];
  G201 = G201 + tvec9*Eeta[109];
  G202 = G202 + tvec9*Eeta[119];
  G210 = G210 + tvec9*Eeta[129];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec9*Eeta[139];
  G212 = G212 + tvec9*Eeta[149];
  G220 = G220 + tvec9*Eeta[159];
  G221 = G221 + tvec9*Eeta[169];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec9*Eeta[179];

  sigma[  0][  2] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  2] *        1.0;
  tvec1 = terms[  3] *        3.0;
  tvec2 = terms[  4] *        3.0;
  tvec3 = terms[  5] *        1.0;
  tvec4 = terms[  9] *        3.0;
  tvec5 = terms[ 10] *        6.0;
  tvec6 = terms[ 11] *        3.0;
  tvec7 = terms[ 15] *        3.0;
  tvec8 = terms[ 16] *        3.0;
  tvec9 = terms[ 20] *        1.0;

  G000 = tvec0*Eeta[ 90];
  G001 = tvec0*Eeta[100];
  G002 = tvec0*Eeta[110];
  G010 = tvec0*Eeta[120];
  G011 = tvec0*Eeta[130];
  G012 = tvec0*Eeta[140];
  G020 = tvec0*Eeta[150];
  G021 = tvec0*Eeta[160];
  G022 = tvec0*Eeta[170];

  G000 = G000 + tvec1*Eeta[ 91];
  G001 = G001 + tvec1*Eeta[101];
  G002 = G002 + tvec1*Eeta[111];
  G010 = G010 + tvec1*Eeta[121];
  G011 = G011 + tvec1*Eeta[131];
  G012 = G012 + tvec1*Eeta[141];
  G020 = G020 + tvec1*Eeta[151];
  G021 = G021 + tvec1*Eeta[161];
  G022 = G022 + tvec1*Eeta[171];

  G000 = G000 + tvec2*Eeta[ 92];
  G001 = G001 + tvec2*Eeta[102];
  G002 = G002 + tvec2*Eeta[112];
  G010 = G010 + tvec2*Eeta[122];
  G011 = G011 + tvec2*Eeta[132];
  G012 = G012 + tvec2*Eeta[142];
  G020 = G020 + tvec2*Eeta[152];
  G021 = G021 + tvec2*Eeta[162];
  G022 = G022 + tvec2*Eeta[172];

  G000 = G000 + tvec3*Eeta[ 93];
  G001 = G001 + tvec3*Eeta[103];
  G002 = G002 + tvec3*Eeta[113];
  G010 = G010 + tvec3*Eeta[123];
  G011 = G011 + tvec3*Eeta[133];
  G012 = G012 + tvec3*Eeta[143];
  G020 = G020 + tvec3*Eeta[153];
  G021 = G021 + tvec3*Eeta[163];
  G022 = G022 + tvec3*Eeta[173];

  G000 = G000 + tvec4*Eeta[ 94];
  G001 = G001 + tvec4*Eeta[104];
  G002 = G002 + tvec4*Eeta[114];
  G010 = G010 + tvec4*Eeta[124];
  G011 = G011 + tvec4*Eeta[134];
  G012 = G012 + tvec4*Eeta[144];
  G020 = G020 + tvec4*Eeta[154];
  G021 = G021 + tvec4*Eeta[164];
  G022 = G022 + tvec4*Eeta[174];

  G000 = G000 + tvec5*Eeta[ 95];
  G001 = G001 + tvec5*Eeta[105];
  G002 = G002 + tvec5*Eeta[115];
  G010 = G010 + tvec5*Eeta[125];
  G011 = G011 + tvec5*Eeta[135];
  G012 = G012 + tvec5*Eeta[145];
  G020 = G020 + tvec5*Eeta[155];
  G021 = G021 + tvec5*Eeta[165];
  G022 = G022 + tvec5*Eeta[175];

  G000 = G000 + tvec6*Eeta[ 96];
  G001 = G001 + tvec6*Eeta[106];
  G002 = G002 + tvec6*Eeta[116];
  G010 = G010 + tvec6*Eeta[126];
  G011 = G011 + tvec6*Eeta[136];
  G012 = G012 + tvec6*Eeta[146];
  G020 = G020 + tvec6*Eeta[156];
  G021 = G021 + tvec6*Eeta[166];
  G022 = G022 + tvec6*Eeta[176];

  G000 = G000 + tvec7*Eeta[ 97];
  G001 = G001 + tvec7*Eeta[107];
  G002 = G002 + tvec7*Eeta[117];
  G010 = G010 + tvec7*Eeta[127];
  G011 = G011 + tvec7*Eeta[137];
  G012 = G012 + tvec7*Eeta[147];
  G020 = G020 + tvec7*Eeta[157];
  G021 = G021 + tvec7*Eeta[167];
  G022 = G022 + tvec7*Eeta[177];

  G000 = G000 + tvec8*Eeta[ 98];
  G001 = G001 + tvec8*Eeta[108];
  G002 = G002 + tvec8*Eeta[118];
  G010 = G010 + tvec8*Eeta[128];
  G011 = G011 + tvec8*Eeta[138];
  G012 = G012 + tvec8*Eeta[148];
  G020 = G020 + tvec8*Eeta[158];
  G021 = G021 + tvec8*Eeta[168];
  G022 = G022 + tvec8*Eeta[178];

  G000 = G000 + tvec9*Eeta[ 99];
  G001 = G001 + tvec9*Eeta[109];
  G002 = G002 + tvec9*Eeta[119];
  G010 = G010 + tvec9*Eeta[129];
  G011 = G011 + tvec9*Eeta[139];
  G012 = G012 + tvec9*Eeta[149];
  G020 = G020 + tvec9*Eeta[159];
  G021 = G021 + tvec9*Eeta[169];
  G022 = G022 + tvec9*Eeta[179];

  H000 = H000 + G000;
  H001 = H001 + G001;
  H002 = H002 + G002;
  H010 = H010 + G010;
  H011 = H011 + G011;
  H012 = H012 + G012;
  H020 = H020 + G020;
  H021 = H021 + G021;
  H022 = H022 + G022;

  tvec0 = terms[  3] *        1.0;
  tvec1 = terms[  4] *        3.0;
  tvec2 = terms[  5] *        3.0;
  tvec3 = terms[  6] *        1.0;
  tvec4 = terms[ 10] *        3.0;
  tvec5 = terms[ 11] *        6.0;
  tvec6 = terms[ 12] *        3.0;
  tvec7 = terms[ 16] *        3.0;
  tvec8 = terms[ 17] *        3.0;
  tvec9 = terms[ 21] *        1.0;

  G100 = tvec0*Eeta[ 90];
  G101 = tvec0*Eeta[100];
  G102 = tvec0*Eeta[110];
  G110 = tvec0*Eeta[120];
  G111 = tvec0*Eeta[130];
  G112 = tvec0*Eeta[140];
  G120 = tvec0*Eeta[150];
  G121 = tvec0*Eeta[160];
  G122 = tvec0*Eeta[170];

  G100 = G100 + tvec1*Eeta[ 91];
  G101 = G101 + tvec1*Eeta[101];
  G102 = G102 + tvec1*Eeta[111];
  G110 = G110 + tvec1*Eeta[121];
  G111 = G111 + tvec1*Eeta[131];
  G112 = G112 + tvec1*Eeta[141];
  G120 = G120 + tvec1*Eeta[151];
  G121 = G121 + tvec1*Eeta[161];
  G122 = G122 + tvec1*Eeta[171];

  G100 = G100 + tvec2*Eeta[ 92];
  G101 = G101 + tvec2*Eeta[102];
  G102 = G102 + tvec2*Eeta[112];
  G110 = G110 + tvec2*Eeta[122];
  G111 = G111 + tvec2*Eeta[132];
  G112 = G112 + tvec2*Eeta[142];
  G120 = G120 + tvec2*Eeta[152];
  G121 = G121 + tvec2*Eeta[162];
  G122 = G122 + tvec2*Eeta[172];

  G100 = G100 + tvec3*Eeta[ 93];
  G101 = G101 + tvec3*Eeta[103];
  G102 = G102 + tvec3*Eeta[113];
  G110 = G110 + tvec3*Eeta[123];
  G111 = G111 + tvec3*Eeta[133];
  G112 = G112 + tvec3*Eeta[143];
  G120 = G120 + tvec3*Eeta[153];
  G121 = G121 + tvec3*Eeta[163];
  G122 = G122 + tvec3*Eeta[173];

  G100 = G100 + tvec4*Eeta[ 94];
  G101 = G101 + tvec4*Eeta[104];
  G102 = G102 + tvec4*Eeta[114];
  G110 = G110 + tvec4*Eeta[124];
  G111 = G111 + tvec4*Eeta[134];
  G112 = G112 + tvec4*Eeta[144];
  G120 = G120 + tvec4*Eeta[154];
  G121 = G121 + tvec4*Eeta[164];
  G122 = G122 + tvec4*Eeta[174];

  G100 = G100 + tvec5*Eeta[ 95];
  G101 = G101 + tvec5*Eeta[105];
  G102 = G102 + tvec5*Eeta[115];
  G110 = G110 + tvec5*Eeta[125];
  G111 = G111 + tvec5*Eeta[135];
  G112 = G112 + tvec5*Eeta[145];
  G120 = G120 + tvec5*Eeta[155];
  G121 = G121 + tvec5*Eeta[165];
  G122 = G122 + tvec5*Eeta[175];

  G100 = G100 + tvec6*Eeta[ 96];
  G101 = G101 + tvec6*Eeta[106];
  G102 = G102 + tvec6*Eeta[116];
  G110 = G110 + tvec6*Eeta[126];
  G111 = G111 + tvec6*Eeta[136];
  G112 = G112 + tvec6*Eeta[146];
  G120 = G120 + tvec6*Eeta[156];
  G121 = G121 + tvec6*Eeta[166];
  G122 = G122 + tvec6*Eeta[176];

  G100 = G100 + tvec7*Eeta[ 97];
  G101 = G101 + tvec7*Eeta[107];
  G102 = G102 + tvec7*Eeta[117];
  G110 = G110 + tvec7*Eeta[127];
  G111 = G111 + tvec7*Eeta[137];
  G112 = G112 + tvec7*Eeta[147];
  G120 = G120 + tvec7*Eeta[157];
  G121 = G121 + tvec7*Eeta[167];
  G122 = G122 + tvec7*Eeta[177];

  G100 = G100 + tvec8*Eeta[ 98];
  G101 = G101 + tvec8*Eeta[108];
  G102 = G102 + tvec8*Eeta[118];
  G110 = G110 + tvec8*Eeta[128];
  G111 = G111 + tvec8*Eeta[138];
  G112 = G112 + tvec8*Eeta[148];
  G120 = G120 + tvec8*Eeta[158];
  G121 = G121 + tvec8*Eeta[168];
  G122 = G122 + tvec8*Eeta[178];

  G100 = G100 + tvec9*Eeta[ 99];
  G101 = G101 + tvec9*Eeta[109];
  G102 = G102 + tvec9*Eeta[119];
  G110 = G110 + tvec9*Eeta[129];
  G111 = G111 + tvec9*Eeta[139];
  G112 = G112 + tvec9*Eeta[149];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec9*Eeta[159];
  G121 = G121 + tvec9*Eeta[169];
  G122 = G122 + tvec9*Eeta[179];

  H100 = H100 + G100;
  H101 = H101 + G101;
  H102 = H102 + G102;
  H110 = H110 + G110;
  H111 = H111 + G111;
  H112 = H112 + G112;
  H120 = H120 + G120;
  H121 = H121 + G121;
  H122 = H122 + G122;

  tvec0 = terms[  9] *        1.0;
  tvec1 = terms[ 10] *        3.0;
  tvec2 = terms[ 11] *        3.0;
  tvec3 = terms[ 12] *        1.0;
  tvec4 = terms[ 15] *        3.0;
  tvec5 = terms[ 16] *        6.0;
  tvec6 = terms[ 17] *        3.0;
  tvec7 = terms[ 20] *        3.0;
  tvec8 = terms[ 21] *        3.0;
  tvec9 = terms[ 24] *        1.0;

  G200 = tvec0*Eeta[ 90];
  G201 = tvec0*Eeta[100];
  G202 = tvec0*Eeta[110];
  G210 = tvec0*Eeta[120];
  G211 = tvec0*Eeta[130];
  G212 = tvec0*Eeta[140];
  G220 = tvec0*Eeta[150];
  G221 = tvec0*Eeta[160];
  G222 = tvec0*Eeta[170];

  G200 = G200 + tvec1*Eeta[ 91];
  G201 = G201 + tvec1*Eeta[101];
  G202 = G202 + tvec1*Eeta[111];
  G210 = G210 + tvec1*Eeta[121];
  G211 = G211 + tvec1*Eeta[131];
  G212 = G212 + tvec1*Eeta[141];
  G220 = G220 + tvec1*Eeta[151];
  G221 = G221 + tvec1*Eeta[161];
  G222 = G222 + tvec1*Eeta[171];

  G200 = G200 + tvec2*Eeta[ 92];
  G201 = G201 + tvec2*Eeta[102];
  G202 = G202 + tvec2*Eeta[112];
  G210 = G210 + tvec2*Eeta[122];
  G211 = G211 + tvec2*Eeta[132];
  G212 = G212 + tvec2*Eeta[142];
  G220 = G220 + tvec2*Eeta[152];
  G221 = G221 + tvec2*Eeta[162];
  G222 = G222 + tvec2*Eeta[172];

  G200 = G200 + tvec3*Eeta[ 93];
  G201 = G201 + tvec3*Eeta[103];
  G202 = G202 + tvec3*Eeta[113];
  G210 = G210 + tvec3*Eeta[123];
  G211 = G211 + tvec3*Eeta[133];
  G212 = G212 + tvec3*Eeta[143];
  G220 = G220 + tvec3*Eeta[153];
  G221 = G221 + tvec3*Eeta[163];
  G222 = G222 + tvec3*Eeta[173];

  G200 = G200 + tvec4*Eeta[ 94];
  G201 = G201 + tvec4*Eeta[104];
  G202 = G202 + tvec4*Eeta[114];
  G210 = G210 + tvec4*Eeta[124];
  G211 = G211 + tvec4*Eeta[134];
  G212 = G212 + tvec4*Eeta[144];
  G220 = G220 + tvec4*Eeta[154];
  G221 = G221 + tvec4*Eeta[164];
  G222 = G222 + tvec4*Eeta[174];

  G200 = G200 + tvec5*Eeta[ 95];
  G201 = G201 + tvec5*Eeta[105];
  G202 = G202 + tvec5*Eeta[115];
  G210 = G210 + tvec5*Eeta[125];
  G211 = G211 + tvec5*Eeta[135];
  G212 = G212 + tvec5*Eeta[145];
  G220 = G220 + tvec5*Eeta[155];
  G221 = G221 + tvec5*Eeta[165];
  G222 = G222 + tvec5*Eeta[175];

  G200 = G200 + tvec6*Eeta[ 96];
  G201 = G201 + tvec6*Eeta[106];
  G202 = G202 + tvec6*Eeta[116];
  G210 = G210 + tvec6*Eeta[126];
  G211 = G211 + tvec6*Eeta[136];
  G212 = G212 + tvec6*Eeta[146];
  G220 = G220 + tvec6*Eeta[156];
  G221 = G221 + tvec6*Eeta[166];
  G222 = G222 + tvec6*Eeta[176];

  G200 = G200 + tvec7*Eeta[ 97];
  G201 = G201 + tvec7*Eeta[107];
  G202 = G202 + tvec7*Eeta[117];
  G210 = G210 + tvec7*Eeta[127];
  G211 = G211 + tvec7*Eeta[137];
  G212 = G212 + tvec7*Eeta[147];
  G220 = G220 + tvec7*Eeta[157];
  G221 = G221 + tvec7*Eeta[167];
  G222 = G222 + tvec7*Eeta[177];

  G200 = G200 + tvec8*Eeta[ 98];
  G201 = G201 + tvec8*Eeta[108];
  G202 = G202 + tvec8*Eeta[118];
  G210 = G210 + tvec8*Eeta[128];
  G211 = G211 + tvec8*Eeta[138];
  G212 = G212 + tvec8*Eeta[148];
  G220 = G220 + tvec8*Eeta[158];
  G221 = G221 + tvec8*Eeta[168];
  G222 = G222 + tvec8*Eeta[178];

  G200 = G200 + tvec9*Eeta[ 99];
  G201 = G201 + tvec9*Eeta[109];
  G202 = G202 + tvec9*Eeta[119];
  G210 = G210 + tvec9*Eeta[129];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec9*Eeta[139];
  G212 = G212 + tvec9*Eeta[149];
  G220 = G220 + tvec9*Eeta[159];
  G221 = G221 + tvec9*Eeta[169];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec9*Eeta[179];

  H200 = H200 + G200;
  H201 = H201 + G201;
  H202 = H202 + G202;
  H210 = H210 + G210;
  H211 = H211 + G211;
  H212 = H212 + G212;
  H220 = H220 + G220;
  H221 = H221 + G221;
  H222 = H222 + G222;

  sigma[  1][  1] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[  8] *        1.0;
  tvec1 = terms[  9] *        3.0;
  tvec2 = terms[ 10] *        3.0;
  tvec3 = terms[ 11] *        1.0;
  tvec4 = terms[ 14] *        3.0;
  tvec5 = terms[ 15] *        6.0;
  tvec6 = terms[ 16] *        3.0;
  tvec7 = terms[ 19] *        3.0;
  tvec8 = terms[ 20] *        3.0;
  tvec9 = terms[ 23] *        1.0;

  G000 = tvec0*Eeta[ 90];
  G001 = tvec0*Eeta[100];
  G002 = tvec0*Eeta[110];
  G010 = tvec0*Eeta[120];
  G011 = tvec0*Eeta[130];
  G012 = tvec0*Eeta[140];
  G020 = tvec0*Eeta[150];
  G021 = tvec0*Eeta[160];
  G022 = tvec0*Eeta[170];

  G000 = G000 + tvec1*Eeta[ 91];
  G001 = G001 + tvec1*Eeta[101];
  G002 = G002 + tvec1*Eeta[111];
  G010 = G010 + tvec1*Eeta[121];
  G011 = G011 + tvec1*Eeta[131];
  G012 = G012 + tvec1*Eeta[141];
  G020 = G020 + tvec1*Eeta[151];
  G021 = G021 + tvec1*Eeta[161];
  G022 = G022 + tvec1*Eeta[171];

  G000 = G000 + tvec2*Eeta[ 92];
  G001 = G001 + tvec2*Eeta[102];
  G002 = G002 + tvec2*Eeta[112];
  G010 = G010 + tvec2*Eeta[122];
  G011 = G011 + tvec2*Eeta[132];
  G012 = G012 + tvec2*Eeta[142];
  G020 = G020 + tvec2*Eeta[152];
  G021 = G021 + tvec2*Eeta[162];
  G022 = G022 + tvec2*Eeta[172];

  G000 = G000 + tvec3*Eeta[ 93];
  G001 = G001 + tvec3*Eeta[103];
  G002 = G002 + tvec3*Eeta[113];
  G010 = G010 + tvec3*Eeta[123];
  G011 = G011 + tvec3*Eeta[133];
  G012 = G012 + tvec3*Eeta[143];
  G020 = G020 + tvec3*Eeta[153];
  G021 = G021 + tvec3*Eeta[163];
  G022 = G022 + tvec3*Eeta[173];

  G000 = G000 + tvec4*Eeta[ 94];
  G001 = G001 + tvec4*Eeta[104];
  G002 = G002 + tvec4*Eeta[114];
  G010 = G010 + tvec4*Eeta[124];
  G011 = G011 + tvec4*Eeta[134];
  G012 = G012 + tvec4*Eeta[144];
  G020 = G020 + tvec4*Eeta[154];
  G021 = G021 + tvec4*Eeta[164];
  G022 = G022 + tvec4*Eeta[174];

  G000 = G000 + tvec5*Eeta[ 95];
  G001 = G001 + tvec5*Eeta[105];
  G002 = G002 + tvec5*Eeta[115];
  G010 = G010 + tvec5*Eeta[125];
  G011 = G011 + tvec5*Eeta[135];
  G012 = G012 + tvec5*Eeta[145];
  G020 = G020 + tvec5*Eeta[155];
  G021 = G021 + tvec5*Eeta[165];
  G022 = G022 + tvec5*Eeta[175];

  G000 = G000 + tvec6*Eeta[ 96];
  G001 = G001 + tvec6*Eeta[106];
  G002 = G002 + tvec6*Eeta[116];
  G010 = G010 + tvec6*Eeta[126];
  G011 = G011 + tvec6*Eeta[136];
  G012 = G012 + tvec6*Eeta[146];
  G020 = G020 + tvec6*Eeta[156];
  G021 = G021 + tvec6*Eeta[166];
  G022 = G022 + tvec6*Eeta[176];

  G000 = G000 + tvec7*Eeta[ 97];
  G001 = G001 + tvec7*Eeta[107];
  G002 = G002 + tvec7*Eeta[117];
  G010 = G010 + tvec7*Eeta[127];
  G011 = G011 + tvec7*Eeta[137];
  G012 = G012 + tvec7*Eeta[147];
  G020 = G020 + tvec7*Eeta[157];
  G021 = G021 + tvec7*Eeta[167];
  G022 = G022 + tvec7*Eeta[177];

  G000 = G000 + tvec8*Eeta[ 98];
  G001 = G001 + tvec8*Eeta[108];
  G002 = G002 + tvec8*Eeta[118];
  G010 = G010 + tvec8*Eeta[128];
  G011 = G011 + tvec8*Eeta[138];
  G012 = G012 + tvec8*Eeta[148];
  G020 = G020 + tvec8*Eeta[158];
  G021 = G021 + tvec8*Eeta[168];
  G022 = G022 + tvec8*Eeta[178];

  G000 = G000 + tvec9*Eeta[ 99];
  G001 = G001 + tvec9*Eeta[109];
  G002 = G002 + tvec9*Eeta[119];
  G010 = G010 + tvec9*Eeta[129];
  G011 = G011 + tvec9*Eeta[139];
  G012 = G012 + tvec9*Eeta[149];
  G020 = G020 + tvec9*Eeta[159];
  G021 = G021 + tvec9*Eeta[169];
  G022 = G022 + tvec9*Eeta[179];

  tvec0 = terms[  9] *        1.0;
  tvec1 = terms[ 10] *        3.0;
  tvec2 = terms[ 11] *        3.0;
  tvec3 = terms[ 12] *        1.0;
  tvec4 = terms[ 15] *        3.0;
  tvec5 = terms[ 16] *        6.0;
  tvec6 = terms[ 17] *        3.0;
  tvec7 = terms[ 20] *        3.0;
  tvec8 = terms[ 21] *        3.0;
  tvec9 = terms[ 24] *        1.0;

  G100 = tvec0*Eeta[ 90];
  G101 = tvec0*Eeta[100];
  G102 = tvec0*Eeta[110];
  G110 = tvec0*Eeta[120];
  G111 = tvec0*Eeta[130];
  G112 = tvec0*Eeta[140];
  G120 = tvec0*Eeta[150];
  G121 = tvec0*Eeta[160];
  G122 = tvec0*Eeta[170];

  G100 = G100 + tvec1*Eeta[ 91];
  G101 = G101 + tvec1*Eeta[101];
  G102 = G102 + tvec1*Eeta[111];
  G110 = G110 + tvec1*Eeta[121];
  G111 = G111 + tvec1*Eeta[131];
  G112 = G112 + tvec1*Eeta[141];
  G120 = G120 + tvec1*Eeta[151];
  G121 = G121 + tvec1*Eeta[161];
  G122 = G122 + tvec1*Eeta[171];

  G100 = G100 + tvec2*Eeta[ 92];
  G101 = G101 + tvec2*Eeta[102];
  G102 = G102 + tvec2*Eeta[112];
  G110 = G110 + tvec2*Eeta[122];
  G111 = G111 + tvec2*Eeta[132];
  G112 = G112 + tvec2*Eeta[142];
  G120 = G120 + tvec2*Eeta[152];
  G121 = G121 + tvec2*Eeta[162];
  G122 = G122 + tvec2*Eeta[172];

  G100 = G100 + tvec3*Eeta[ 93];
  G101 = G101 + tvec3*Eeta[103];
  G102 = G102 + tvec3*Eeta[113];
  G110 = G110 + tvec3*Eeta[123];
  G111 = G111 + tvec3*Eeta[133];
  G112 = G112 + tvec3*Eeta[143];
  G120 = G120 + tvec3*Eeta[153];
  G121 = G121 + tvec3*Eeta[163];
  G122 = G122 + tvec3*Eeta[173];

  G100 = G100 + tvec4*Eeta[ 94];
  G101 = G101 + tvec4*Eeta[104];
  G102 = G102 + tvec4*Eeta[114];
  G110 = G110 + tvec4*Eeta[124];
  G111 = G111 + tvec4*Eeta[134];
  G112 = G112 + tvec4*Eeta[144];
  G120 = G120 + tvec4*Eeta[154];
  G121 = G121 + tvec4*Eeta[164];
  G122 = G122 + tvec4*Eeta[174];

  G100 = G100 + tvec5*Eeta[ 95];
  G101 = G101 + tvec5*Eeta[105];
  G102 = G102 + tvec5*Eeta[115];
  G110 = G110 + tvec5*Eeta[125];
  G111 = G111 + tvec5*Eeta[135];
  G112 = G112 + tvec5*Eeta[145];
  G120 = G120 + tvec5*Eeta[155];
  G121 = G121 + tvec5*Eeta[165];
  G122 = G122 + tvec5*Eeta[175];

  G100 = G100 + tvec6*Eeta[ 96];
  G101 = G101 + tvec6*Eeta[106];
  G102 = G102 + tvec6*Eeta[116];
  G110 = G110 + tvec6*Eeta[126];
  G111 = G111 + tvec6*Eeta[136];
  G112 = G112 + tvec6*Eeta[146];
  G120 = G120 + tvec6*Eeta[156];
  G121 = G121 + tvec6*Eeta[166];
  G122 = G122 + tvec6*Eeta[176];

  G100 = G100 + tvec7*Eeta[ 97];
  G101 = G101 + tvec7*Eeta[107];
  G102 = G102 + tvec7*Eeta[117];
  G110 = G110 + tvec7*Eeta[127];
  G111 = G111 + tvec7*Eeta[137];
  G112 = G112 + tvec7*Eeta[147];
  G120 = G120 + tvec7*Eeta[157];
  G121 = G121 + tvec7*Eeta[167];
  G122 = G122 + tvec7*Eeta[177];

  G100 = G100 + tvec8*Eeta[ 98];
  G101 = G101 + tvec8*Eeta[108];
  G102 = G102 + tvec8*Eeta[118];
  G110 = G110 + tvec8*Eeta[128];
  G111 = G111 + tvec8*Eeta[138];
  G112 = G112 + tvec8*Eeta[148];
  G120 = G120 + tvec8*Eeta[158];
  G121 = G121 + tvec8*Eeta[168];
  G122 = G122 + tvec8*Eeta[178];

  G100 = G100 + tvec9*Eeta[ 99];
  G101 = G101 + tvec9*Eeta[109];
  G102 = G102 + tvec9*Eeta[119];
  G110 = G110 + tvec9*Eeta[129];
  G111 = G111 + tvec9*Eeta[139];
  G112 = G112 + tvec9*Eeta[149];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec9*Eeta[159];
  G121 = G121 + tvec9*Eeta[169];
  G122 = G122 + tvec9*Eeta[179];

  tvec0 = terms[ 14] *        1.0;
  tvec1 = terms[ 15] *        3.0;
  tvec2 = terms[ 16] *        3.0;
  tvec3 = terms[ 17] *        1.0;
  tvec4 = terms[ 19] *        3.0;
  tvec5 = terms[ 20] *        6.0;
  tvec6 = terms[ 21] *        3.0;
  tvec7 = terms[ 23] *        3.0;
  tvec8 = terms[ 24] *        3.0;
  tvec9 = terms[ 26] *        1.0;

  G200 = tvec0*Eeta[ 90];
  G201 = tvec0*Eeta[100];
  G202 = tvec0*Eeta[110];
  G210 = tvec0*Eeta[120];
  G211 = tvec0*Eeta[130];
  G212 = tvec0*Eeta[140];
  G220 = tvec0*Eeta[150];
  G221 = tvec0*Eeta[160];
  G222 = tvec0*Eeta[170];

  G200 = G200 + tvec1*Eeta[ 91];
  G201 = G201 + tvec1*Eeta[101];
  G202 = G202 + tvec1*Eeta[111];
  G210 = G210 + tvec1*Eeta[121];
  G211 = G211 + tvec1*Eeta[131];
  G212 = G212 + tvec1*Eeta[141];
  G220 = G220 + tvec1*Eeta[151];
  G221 = G221 + tvec1*Eeta[161];
  G222 = G222 + tvec1*Eeta[171];

  G200 = G200 + tvec2*Eeta[ 92];
  G201 = G201 + tvec2*Eeta[102];
  G202 = G202 + tvec2*Eeta[112];
  G210 = G210 + tvec2*Eeta[122];
  G211 = G211 + tvec2*Eeta[132];
  G212 = G212 + tvec2*Eeta[142];
  G220 = G220 + tvec2*Eeta[152];
  G221 = G221 + tvec2*Eeta[162];
  G222 = G222 + tvec2*Eeta[172];

  G200 = G200 + tvec3*Eeta[ 93];
  G201 = G201 + tvec3*Eeta[103];
  G202 = G202 + tvec3*Eeta[113];
  G210 = G210 + tvec3*Eeta[123];
  G211 = G211 + tvec3*Eeta[133];
  G212 = G212 + tvec3*Eeta[143];
  G220 = G220 + tvec3*Eeta[153];
  G221 = G221 + tvec3*Eeta[163];
  G222 = G222 + tvec3*Eeta[173];

  G200 = G200 + tvec4*Eeta[ 94];
  G201 = G201 + tvec4*Eeta[104];
  G202 = G202 + tvec4*Eeta[114];
  G210 = G210 + tvec4*Eeta[124];
  G211 = G211 + tvec4*Eeta[134];
  G212 = G212 + tvec4*Eeta[144];
  G220 = G220 + tvec4*Eeta[154];
  G221 = G221 + tvec4*Eeta[164];
  G222 = G222 + tvec4*Eeta[174];

  G200 = G200 + tvec5*Eeta[ 95];
  G201 = G201 + tvec5*Eeta[105];
  G202 = G202 + tvec5*Eeta[115];
  G210 = G210 + tvec5*Eeta[125];
  G211 = G211 + tvec5*Eeta[135];
  G212 = G212 + tvec5*Eeta[145];
  G220 = G220 + tvec5*Eeta[155];
  G221 = G221 + tvec5*Eeta[165];
  G222 = G222 + tvec5*Eeta[175];

  G200 = G200 + tvec6*Eeta[ 96];
  G201 = G201 + tvec6*Eeta[106];
  G202 = G202 + tvec6*Eeta[116];
  G210 = G210 + tvec6*Eeta[126];
  G211 = G211 + tvec6*Eeta[136];
  G212 = G212 + tvec6*Eeta[146];
  G220 = G220 + tvec6*Eeta[156];
  G221 = G221 + tvec6*Eeta[166];
  G222 = G222 + tvec6*Eeta[176];

  G200 = G200 + tvec7*Eeta[ 97];
  G201 = G201 + tvec7*Eeta[107];
  G202 = G202 + tvec7*Eeta[117];
  G210 = G210 + tvec7*Eeta[127];
  G211 = G211 + tvec7*Eeta[137];
  G212 = G212 + tvec7*Eeta[147];
  G220 = G220 + tvec7*Eeta[157];
  G221 = G221 + tvec7*Eeta[167];
  G222 = G222 + tvec7*Eeta[177];

  G200 = G200 + tvec8*Eeta[ 98];
  G201 = G201 + tvec8*Eeta[108];
  G202 = G202 + tvec8*Eeta[118];
  G210 = G210 + tvec8*Eeta[128];
  G211 = G211 + tvec8*Eeta[138];
  G212 = G212 + tvec8*Eeta[148];
  G220 = G220 + tvec8*Eeta[158];
  G221 = G221 + tvec8*Eeta[168];
  G222 = G222 + tvec8*Eeta[178];

  G200 = G200 + tvec9*Eeta[ 99];
  G201 = G201 + tvec9*Eeta[109];
  G202 = G202 + tvec9*Eeta[119];
  G210 = G210 + tvec9*Eeta[129];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec9*Eeta[139];
  G212 = G212 + tvec9*Eeta[149];
  G220 = G220 + tvec9*Eeta[159];
  G221 = G221 + tvec9*Eeta[169];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec9*Eeta[179];

  sigma[  1][  2] = two1nu * t2;

  t2 = 0.0;
  tvec0 = terms[ 13] *        1.0;
  tvec1 = terms[ 14] *        3.0;
  tvec2 = terms[ 15] *        3.0;
  tvec3 = terms[ 16] *        1.0;
  tvec4 = terms[ 18] *        3.0;
  tvec5 = terms[ 19] *        6.0;
  tvec6 = terms[ 20] *        3.0;
  tvec7 = terms[ 22] *        3.0;
  tvec8 = terms[ 23] *        3.0;
  tvec9 = terms[ 25] *        1.0;

  G000 = tvec0*Eeta[ 90];
  G001 = tvec0*Eeta[100];
  G002 = tvec0*Eeta[110];
  G010 = tvec0*Eeta[120];
  G011 = tvec0*Eeta[130];
  G012 = tvec0*Eeta[140];
  G020 = tvec0*Eeta[150];
  G021 = tvec0*Eeta[160];
  G022 = tvec0*Eeta[170];

  G000 = G000 + tvec1*Eeta[ 91];
  G001 = G001 + tvec1*Eeta[101];
  G002 = G002 + tvec1*Eeta[111];
  G010 = G010 + tvec1*Eeta[121];
  G011 = G011 + tvec1*Eeta[131];
  G012 = G012 + tvec1*Eeta[141];
  G020 = G020 + tvec1*Eeta[151];
  G021 = G021 + tvec1*Eeta[161];
  G022 = G022 + tvec1*Eeta[171];

  G000 = G000 + tvec2*Eeta[ 92];
  G001 = G001 + tvec2*Eeta[102];
  G002 = G002 + tvec2*Eeta[112];
  G010 = G010 + tvec2*Eeta[122];
  G011 = G011 + tvec2*Eeta[132];
  G012 = G012 + tvec2*Eeta[142];
  G020 = G020 + tvec2*Eeta[152];
  G021 = G021 + tvec2*Eeta[162];
  G022 = G022 + tvec2*Eeta[172];

  G000 = G000 + tvec3*Eeta[ 93];
  G001 = G001 + tvec3*Eeta[103];
  G002 = G002 + tvec3*Eeta[113];
  G010 = G010 + tvec3*Eeta[123];
  G011 = G011 + tvec3*Eeta[133];
  G012 = G012 + tvec3*Eeta[143];
  G020 = G020 + tvec3*Eeta[153];
  G021 = G021 + tvec3*Eeta[163];
  G022 = G022 + tvec3*Eeta[173];

  G000 = G000 + tvec4*Eeta[ 94];
  G001 = G001 + tvec4*Eeta[104];
  G002 = G002 + tvec4*Eeta[114];
  G010 = G010 + tvec4*Eeta[124];
  G011 = G011 + tvec4*Eeta[134];
  G012 = G012 + tvec4*Eeta[144];
  G020 = G020 + tvec4*Eeta[154];
  G021 = G021 + tvec4*Eeta[164];
  G022 = G022 + tvec4*Eeta[174];

  G000 = G000 + tvec5*Eeta[ 95];
  G001 = G001 + tvec5*Eeta[105];
  G002 = G002 + tvec5*Eeta[115];
  G010 = G010 + tvec5*Eeta[125];
  G011 = G011 + tvec5*Eeta[135];
  G012 = G012 + tvec5*Eeta[145];
  G020 = G020 + tvec5*Eeta[155];
  G021 = G021 + tvec5*Eeta[165];
  G022 = G022 + tvec5*Eeta[175];

  G000 = G000 + tvec6*Eeta[ 96];
  G001 = G001 + tvec6*Eeta[106];
  G002 = G002 + tvec6*Eeta[116];
  G010 = G010 + tvec6*Eeta[126];
  G011 = G011 + tvec6*Eeta[136];
  G012 = G012 + tvec6*Eeta[146];
  G020 = G020 + tvec6*Eeta[156];
  G021 = G021 + tvec6*Eeta[166];
  G022 = G022 + tvec6*Eeta[176];

  G000 = G000 + tvec7*Eeta[ 97];
  G001 = G001 + tvec7*Eeta[107];
  G002 = G002 + tvec7*Eeta[117];
  G010 = G010 + tvec7*Eeta[127];
  G011 = G011 + tvec7*Eeta[137];
  G012 = G012 + tvec7*Eeta[147];
  G020 = G020 + tvec7*Eeta[157];
  G021 = G021 + tvec7*Eeta[167];
  G022 = G022 + tvec7*Eeta[177];

  G000 = G000 + tvec8*Eeta[ 98];
  G001 = G001 + tvec8*Eeta[108];
  G002 = G002 + tvec8*Eeta[118];
  G010 = G010 + tvec8*Eeta[128];
  G011 = G011 + tvec8*Eeta[138];
  G012 = G012 + tvec8*Eeta[148];
  G020 = G020 + tvec8*Eeta[158];
  G021 = G021 + tvec8*Eeta[168];
  G022 = G022 + tvec8*Eeta[178];

  G000 = G000 + tvec9*Eeta[ 99];
  G001 = G001 + tvec9*Eeta[109];
  G002 = G002 + tvec9*Eeta[119];
  G010 = G010 + tvec9*Eeta[129];
  G011 = G011 + tvec9*Eeta[139];
  G012 = G012 + tvec9*Eeta[149];
  G020 = G020 + tvec9*Eeta[159];
  G021 = G021 + tvec9*Eeta[169];
  G022 = G022 + tvec9*Eeta[179];

  H000 = H000 + G000;
  H001 = H001 + G001;
  H002 = H002 + G002;
  H010 = H010 + G010;
  H011 = H011 + G011;
  H012 = H012 + G012;
  H020 = H020 + G020;
  H021 = H021 + G021;
  H022 = H022 + G022;

  tvec0 = terms[ 14] *        1.0;
  tvec1 = terms[ 15] *        3.0;
  tvec2 = terms[ 16] *        3.0;
  tvec3 = terms[ 17] *        1.0;
  tvec4 = terms[ 19] *        3.0;
  tvec5 = terms[ 20] *        6.0;
  tvec6 = terms[ 21] *        3.0;
  tvec7 = terms[ 23] *        3.0;
  tvec8 = terms[ 24] *        3.0;
  tvec9 = terms[ 26] *        1.0;

  G100 = tvec0*Eeta[ 90];
  G101 = tvec0*Eeta[100];
  G102 = tvec0*Eeta[110];
  G110 = tvec0*Eeta[120];
  G111 = tvec0*Eeta[130];
  G112 = tvec0*Eeta[140];
  G120 = tvec0*Eeta[150];
  G121 = tvec0*Eeta[160];
  G122 = tvec0*Eeta[170];

  G100 = G100 + tvec1*Eeta[ 91];
  G101 = G101 + tvec1*Eeta[101];
  G102 = G102 + tvec1*Eeta[111];
  G110 = G110 + tvec1*Eeta[121];
  G111 = G111 + tvec1*Eeta[131];
  G112 = G112 + tvec1*Eeta[141];
  G120 = G120 + tvec1*Eeta[151];
  G121 = G121 + tvec1*Eeta[161];
  G122 = G122 + tvec1*Eeta[171];

  G100 = G100 + tvec2*Eeta[ 92];
  G101 = G101 + tvec2*Eeta[102];
  G102 = G102 + tvec2*Eeta[112];
  G110 = G110 + tvec2*Eeta[122];
  G111 = G111 + tvec2*Eeta[132];
  G112 = G112 + tvec2*Eeta[142];
  G120 = G120 + tvec2*Eeta[152];
  G121 = G121 + tvec2*Eeta[162];
  G122 = G122 + tvec2*Eeta[172];

  G100 = G100 + tvec3*Eeta[ 93];
  G101 = G101 + tvec3*Eeta[103];
  G102 = G102 + tvec3*Eeta[113];
  G110 = G110 + tvec3*Eeta[123];
  G111 = G111 + tvec3*Eeta[133];
  G112 = G112 + tvec3*Eeta[143];
  G120 = G120 + tvec3*Eeta[153];
  G121 = G121 + tvec3*Eeta[163];
  G122 = G122 + tvec3*Eeta[173];

  G100 = G100 + tvec4*Eeta[ 94];
  G101 = G101 + tvec4*Eeta[104];
  G102 = G102 + tvec4*Eeta[114];
  G110 = G110 + tvec4*Eeta[124];
  G111 = G111 + tvec4*Eeta[134];
  G112 = G112 + tvec4*Eeta[144];
  G120 = G120 + tvec4*Eeta[154];
  G121 = G121 + tvec4*Eeta[164];
  G122 = G122 + tvec4*Eeta[174];

  G100 = G100 + tvec5*Eeta[ 95];
  G101 = G101 + tvec5*Eeta[105];
  G102 = G102 + tvec5*Eeta[115];
  G110 = G110 + tvec5*Eeta[125];
  G111 = G111 + tvec5*Eeta[135];
  G112 = G112 + tvec5*Eeta[145];
  G120 = G120 + tvec5*Eeta[155];
  G121 = G121 + tvec5*Eeta[165];
  G122 = G122 + tvec5*Eeta[175];

  G100 = G100 + tvec6*Eeta[ 96];
  G101 = G101 + tvec6*Eeta[106];
  G102 = G102 + tvec6*Eeta[116];
  G110 = G110 + tvec6*Eeta[126];
  G111 = G111 + tvec6*Eeta[136];
  G112 = G112 + tvec6*Eeta[146];
  G120 = G120 + tvec6*Eeta[156];
  G121 = G121 + tvec6*Eeta[166];
  G122 = G122 + tvec6*Eeta[176];

  G100 = G100 + tvec7*Eeta[ 97];
  G101 = G101 + tvec7*Eeta[107];
  G102 = G102 + tvec7*Eeta[117];
  G110 = G110 + tvec7*Eeta[127];
  G111 = G111 + tvec7*Eeta[137];
  G112 = G112 + tvec7*Eeta[147];
  G120 = G120 + tvec7*Eeta[157];
  G121 = G121 + tvec7*Eeta[167];
  G122 = G122 + tvec7*Eeta[177];

  G100 = G100 + tvec8*Eeta[ 98];
  G101 = G101 + tvec8*Eeta[108];
  G102 = G102 + tvec8*Eeta[118];
  G110 = G110 + tvec8*Eeta[128];
  G111 = G111 + tvec8*Eeta[138];
  G112 = G112 + tvec8*Eeta[148];
  G120 = G120 + tvec8*Eeta[158];
  G121 = G121 + tvec8*Eeta[168];
  G122 = G122 + tvec8*Eeta[178];

  G100 = G100 + tvec9*Eeta[ 99];
  G101 = G101 + tvec9*Eeta[109];
  G102 = G102 + tvec9*Eeta[119];
  G110 = G110 + tvec9*Eeta[129];
  G111 = G111 + tvec9*Eeta[139];
  G112 = G112 + tvec9*Eeta[149];
  t2 = t2 + (G012 - G102);
  G120 = G120 + tvec9*Eeta[159];
  G121 = G121 + tvec9*Eeta[169];
  G122 = G122 + tvec9*Eeta[179];

  H100 = H100 + G100;
  H101 = H101 + G101;
  H102 = H102 + G102;
  H110 = H110 + G110;
  H111 = H111 + G111;
  H112 = H112 + G112;
  H120 = H120 + G120;
  H121 = H121 + G121;
  H122 = H122 + G122;

  tvec0 = terms[ 18] *        1.0;
  tvec1 = terms[ 19] *        3.0;
  tvec2 = terms[ 20] *        3.0;
  tvec3 = terms[ 21] *        1.0;
  tvec4 = terms[ 22] *        3.0;
  tvec5 = terms[ 23] *        6.0;
  tvec6 = terms[ 24] *        3.0;
  tvec7 = terms[ 25] *        3.0;
  tvec8 = terms[ 26] *        3.0;
  tvec9 = terms[ 27] *        1.0;

  G200 = tvec0*Eeta[ 90];
  G201 = tvec0*Eeta[100];
  G202 = tvec0*Eeta[110];
  G210 = tvec0*Eeta[120];
  G211 = tvec0*Eeta[130];
  G212 = tvec0*Eeta[140];
  G220 = tvec0*Eeta[150];
  G221 = tvec0*Eeta[160];
  G222 = tvec0*Eeta[170];

  G200 = G200 + tvec1*Eeta[ 91];
  G201 = G201 + tvec1*Eeta[101];
  G202 = G202 + tvec1*Eeta[111];
  G210 = G210 + tvec1*Eeta[121];
  G211 = G211 + tvec1*Eeta[131];
  G212 = G212 + tvec1*Eeta[141];
  G220 = G220 + tvec1*Eeta[151];
  G221 = G221 + tvec1*Eeta[161];
  G222 = G222 + tvec1*Eeta[171];

  G200 = G200 + tvec2*Eeta[ 92];
  G201 = G201 + tvec2*Eeta[102];
  G202 = G202 + tvec2*Eeta[112];
  G210 = G210 + tvec2*Eeta[122];
  G211 = G211 + tvec2*Eeta[132];
  G212 = G212 + tvec2*Eeta[142];
  G220 = G220 + tvec2*Eeta[152];
  G221 = G221 + tvec2*Eeta[162];
  G222 = G222 + tvec2*Eeta[172];

  G200 = G200 + tvec3*Eeta[ 93];
  G201 = G201 + tvec3*Eeta[103];
  G202 = G202 + tvec3*Eeta[113];
  G210 = G210 + tvec3*Eeta[123];
  G211 = G211 + tvec3*Eeta[133];
  G212 = G212 + tvec3*Eeta[143];
  G220 = G220 + tvec3*Eeta[153];
  G221 = G221 + tvec3*Eeta[163];
  G222 = G222 + tvec3*Eeta[173];

  G200 = G200 + tvec4*Eeta[ 94];
  G201 = G201 + tvec4*Eeta[104];
  G202 = G202 + tvec4*Eeta[114];
  G210 = G210 + tvec4*Eeta[124];
  G211 = G211 + tvec4*Eeta[134];
  G212 = G212 + tvec4*Eeta[144];
  G220 = G220 + tvec4*Eeta[154];
  G221 = G221 + tvec4*Eeta[164];
  G222 = G222 + tvec4*Eeta[174];

  G200 = G200 + tvec5*Eeta[ 95];
  G201 = G201 + tvec5*Eeta[105];
  G202 = G202 + tvec5*Eeta[115];
  G210 = G210 + tvec5*Eeta[125];
  G211 = G211 + tvec5*Eeta[135];
  G212 = G212 + tvec5*Eeta[145];
  G220 = G220 + tvec5*Eeta[155];
  G221 = G221 + tvec5*Eeta[165];
  G222 = G222 + tvec5*Eeta[175];

  G200 = G200 + tvec6*Eeta[ 96];
  G201 = G201 + tvec6*Eeta[106];
  G202 = G202 + tvec6*Eeta[116];
  G210 = G210 + tvec6*Eeta[126];
  G211 = G211 + tvec6*Eeta[136];
  G212 = G212 + tvec6*Eeta[146];
  G220 = G220 + tvec6*Eeta[156];
  G221 = G221 + tvec6*Eeta[166];
  G222 = G222 + tvec6*Eeta[176];

  G200 = G200 + tvec7*Eeta[ 97];
  G201 = G201 + tvec7*Eeta[107];
  G202 = G202 + tvec7*Eeta[117];
  G210 = G210 + tvec7*Eeta[127];
  G211 = G211 + tvec7*Eeta[137];
  G212 = G212 + tvec7*Eeta[147];
  G220 = G220 + tvec7*Eeta[157];
  G221 = G221 + tvec7*Eeta[167];
  G222 = G222 + tvec7*Eeta[177];

  G200 = G200 + tvec8*Eeta[ 98];
  G201 = G201 + tvec8*Eeta[108];
  G202 = G202 + tvec8*Eeta[118];
  G210 = G210 + tvec8*Eeta[128];
  G211 = G211 + tvec8*Eeta[138];
  G212 = G212 + tvec8*Eeta[148];
  G220 = G220 + tvec8*Eeta[158];
  G221 = G221 + tvec8*Eeta[168];
  G222 = G222 + tvec8*Eeta[178];

  G200 = G200 + tvec9*Eeta[ 99];
  G201 = G201 + tvec9*Eeta[109];
  G202 = G202 + tvec9*Eeta[119];
  G210 = G210 + tvec9*Eeta[129];
  t2 = t2 + (G120 - G210);
  G211 = G211 + tvec9*Eeta[139];
  G212 = G212 + tvec9*Eeta[149];
  G220 = G220 + tvec9*Eeta[159];
  G221 = G221 + tvec9*Eeta[169];
  t2 = t2 + (G201 - G021);
  G222 = G222 + tvec9*Eeta[179];

  H200 = H200 + G200;
  H201 = H201 + G201;
  H202 = H202 + G202;
  H210 = H210 + G210;
  H211 = H211 + G211;
  H212 = H212 + G212;
  H220 = H220 + G220;
  H221 = H221 + G221;
  H222 = H222 + G222;

  sigma[  2][  2] = two1nu * t2;

  t11 = H120 - H210;
  t12 = H120 - H210;
  sigma[  0][  0] = sigma[  0][  0] + (t11 + t12);
  t11 = H200 - H020;
  t12 = H121 - H211;
  sigma[  0][  1] = sigma[  0][  1] + (t11 + t12);
  t11 = H010 - H100;
  t12 = H122 - H212;
  sigma[  0][  2] = sigma[  0][  2] + (t11 + t12);
  t11 = H201 - H021;
  t12 = H201 - H021;
  sigma[  1][  1] = sigma[  1][  1] + (t11 + t12);
  t11 = H011 - H101;
  t12 = H202 - H022;
  sigma[  1][  2] = sigma[  1][  2] + (t11 + t12);
  t11 = H012 - H102;
  t12 = H012 - H102;
  sigma[  2][  2] = sigma[  2][  2] + (t11 + t12);

  t3 = 0.0;
  t3 = t3 + H120 - H210;
  t3 = t3 + H201 - H021;
  t3 = t3 + H012 - H102;
  t3 = t3 * two1nu;
sigma[  0][  0] = sigma[  0][  0] - t3;
sigma[  1][  1] = sigma[  1][  1] - t3;
sigma[  2][  2] = sigma[  2][  2] - t3;

  t =    -1.66666666666666657415e-01 * mu8pi;
  sigma[  0][  0] = t * sigma[  0][  0];
  sigma[  0][  1] = t * sigma[  0][  1];
  sigma[  0][  2] = t * sigma[  0][  2];
  sigma[  1][  1] = t * sigma[  1][  1];
  sigma[  1][  2] = t * sigma[  1][  2];
  sigma[  2][  2] = t * sigma[  2][  2];
  sigma[1][0] = sigma[0][1];
  sigma[2][0] = sigma[0][2];
  sigma[2][1] = sigma[1][2];
}



/*       3      180             0.007532             0.000000     1.978565308132407e-02*/
