/**************************************************************************
 *
 *      Module:  TAYLOR_FMFunctions.cc
 *
 *      Includes functions:
 *          DM()
 *          dmsym3()
 *          dmsym3int()
 *          EvalTaylor()
 *          FMShift()
 *          gcd0()
 *          ipow()
 *          makeeta()
 *          makeqtab()
 *          makeftabs()
 *          MkTaylor()
 *          MkTaylorAnisotropic()
 *          MkTaylorIsotropic()
 *          rdinit()
 *          TaylorShift()
 *
 *      NOTE: Functions are specific to the Taylor FMM method.
 *            Most of these functions are very nearly copied from code put
 *            together by Tomas Oppelstrup. 
 *
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "Home.h"
#include "FM.h"

#ifdef TAYLORFMM

#ifdef ANISOTROPIC
#include "G2sigma_tilde.h"
#endif

/*
static int rveclist[][3] = { {6,0,0} , {6,3,0} , {6,6,0} , {6,3,3} ,
                             {6,6,3} , {6,6,6} , {6,2,0} , {6,4,0} ,
                             {6,2,2} , {6,4,2} , {6,6,2} , {6,4,4} ,
                             {6,6,4} };
*/
static int rveclist[][3] = { {12, 0, 0} , {12,12, 0} , {12,12,12} ,
			     {12, 6, 0} , {12, 6, 6} , {12,12, 6} ,
			     {12, 4, 0} , {12, 4, 4} , {12, 8, 0} ,
			     {12, 8, 4} , {12, 8, 8} , {12, 12,4} ,
			     {12,12, 8} , {12, 3, 0} , {12, 3, 3} ,
			     {12, 6, 3} , {12, 9, 0} , {12, 9, 3} ,
			     {12, 9, 6} , {12, 9, 9} , {12,12, 3} ,
			     {12,12, 9} };


/*---------------------------------------------------------------------------
 *
 *      Function:     gcd0
 *      Description:  Greatest common divisor computation.  Returns
 *                    a non-negative greatest common divisor of a and b.
 *                    If the argument of smallest magnitude is zero, the
 *                    magnitude of the other argument is returned.
 *
 *--------------------------------------------------------------------------*/
int gcd0(int a1, int b1) {
    int a, b;

    if (a1 < 0) {
        a1 = -a1;
    }

    if (b1 < 0) {
        b1 = -b1;
    }

    if (a1 > b1) {
        a = a1;
        b = b1;
    } else {
        a = b1;
        b = a1;
    }

    if (b == 0) {
        return a;
    }

    while (1) {
        int m = a % b;

        if (m == 0) {
            break;
        }

        a = b;
        b = m;
    }

    return(b);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    ipow
 *      Description: Efficient and numerically sound integral power
 *                   computation
 *
 *-------------------------------------------------------------------------*/
real8 ipow(real8 x, int n)
{
        real8 xn = x, y = 1.0;

        if (n < 0) {
            xn = 1.0/xn;
            n = -n;
        }

        while (n > 0) {
            if (n & 0x01) y = y*xn;
            xn = xn*xn;
            n = n >> 1;
        }

        return(y);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    DM
 *      Description: 
 *
 *      Arguments:
 *          n       
 *          rvec    
 *          terms   
 *          npows   
 *
 *-------------------------------------------------------------------------*/
static void DM(int n, real8 *rvec, real8 *terms, int npows[][3])
{
        int   i, nx, ny, nz, a, b, c, m, p;
        real8 dfact[2*MAXORDER-3+2], fact[MAXORDER+1], gtab[MAXORDER+1][3];
        real8 mparity, f, ca, cb, cc, rfact, gx, gy, gz, r, x, y, z, xn, yn, zn;

/*
 *      Compute table for double factorial
 */
        dfact[0] = 1.0;

        for (i = 1; i <= 2*n-3; i += 2) {
            dfact[i+1] = dfact[i-1]*i;
        }
  
/*
 *      Compute table for factorial
 */
        fact[0] = 1.0;

        for (i = 1; i <= n; i++) {
            fact[i] = fact[i-1]*i;
        }

        r = 0.0;

        for (i = 0; i < 3; i++) {
            r = r + rvec[i]*rvec[i];
        }

        r = sqrt(r);
        rfact = ipow(-1/r, n-1);

        gx = rvec[0] / r;
        gy = rvec[1] / r;
        gz = rvec[2] / r;

/*
 *      This might be overkill, but has better accuracy than
 *          x(n) = gx*x(n-1)
 *
 *      This is just an unrolled version of ipow()
 */
        for (i = 0; i <= n; i++) {
            p = i;

            xn = gx;
            yn = gy;
            zn = gz;

            x = y = z = 1.0;

            while (p > 0) {

                if (p & 1) {
                    x = x*xn;
                    y = y*yn;
                    z = z*zn;
                }

                xn = xn*xn;
                yn = yn*yn;
                zn = zn*zn;
                p = p >> 1;
            }
            gtab[i][0] = x;
            gtab[i][1] = y;
            gtab[i][2] = z;
        }

        i = 0;

        for (nz = 0; nz <= n; nz++) {
            for (ny = 0; ny <= n-nz; ny++) {
                nx = n-nz-ny;

                npows[i][0] = nx;
                npows[i][1] = ny;
                npows[i][2] = nz;

                terms[i] = 0;
                ca = 1.0;

                for (a = 0; a <= nx; a += 2) {
                    cb = 1.0;
                    for (b = 0; b <= ny; b += 2) {
                        cc = 1.0;

                        m = (a+b)/2;
                        mparity = 1-2*(m%2);

                        for (c = 0; c <= nz; c += 2) {
                            f = mparity * dfact[2*(n-m)-3+1] * 
                                gtab[nx-a][0] *
                                gtab[ny-b][1] *
                                gtab[nz-c][2];

                            f = f * ca * cb * cc;
                            terms[i] = terms[i] + f;
                            cc = cc * (nz-c) * (nz-c-1) / (c+2);
                            m = m + 1;
                            mparity = -mparity;
                        }
                        cb = cb * (ny-b) * (ny-b-1) / (b+2);
                    }
                    ca = ca * (nx-a) * (nx-a-1) / (a+2);
                }
                terms[i] = terms[i] * fact[n] /
                           (fact[nx]*fact[ny]*fact[nz]) * rfact;
                i = i+1;
            }
        }

        return;
}


static void dmsym3init(int n,int M[][13][13],real8 drdata[])
{
  static int p[][3] = { {0,1,2} , {1,2,0} , {2,0,1} ,
                        {0,2,1} , {1,0,2} , {2,1,0} };
  int powvec[(MAXORDER+3+3)*(MAXORDER+3+2)*(MAXORDER+3+1)/6][3];
  int nold,i,i0,k,m,ir,nx,ny,nz;
  vector r;
  real8 fact[2*(MAXORDER+3)+1];

  nold = M[0][0][0]-3;
  if(nold < -1) nold = -1;

  fact[0] = 1.0;
/*
 * Bug fix...
 *
  for(i = 1; i<=n; i++) fact[i] = fact[i-1]*(real8) i;
 */
  for(i = 1; i<=n+3; i++) fact[i] = fact[i-1]*(real8) i;

  for(i = 0; i<13*13*13; i++) M[i/(13*13)][(i/13)%13][i%13] = -1;

  m = ( (NMAX+3+3)*(NMAX+3+2)*(NMAX+3+1) - (3+3)*(3+2)*(3+1) )/6;
  for(ir = 0; ir<sizeof(rveclist)/(3*sizeof(int)); ir++) {
    for(i = nold+1; i<=n; i++) {
      i0 = ( (i+3+2)*(i+3+1)*(i+3) - (3+2)*(3+1)*3 )/6;
      for(k = 0; k<3; k++) r[k] = rveclist[ir][k];
      /*printf("@ %s:%d: ir=%d i=%d n=%d i0=%d m=%d ir*m+i0=%d len(powvec)=%d\n",
	__FILE__,__LINE__,ir,i,n,i0,m,ir*m+i0,sizeof(powvec)/sizeof(int)/3);*/
      if(sizeof(powvec)/(3*sizeof(int))-i0 < (i+3+2)*(i+3+1)/2) {
	Fatal("@ %s:%d: powvec is too small!!!", __FILE__,__LINE__);
      }
      DM(i+3,r,&drdata[ir*m+i0],&powvec[i0]);

      /* Remove combinatorial factor now, since that needs to be
         done before use in mktaylor */
      k = 0;
      for(nz = 0; nz<=i+3; nz++)
        for(ny = 0; ny<=i+3-nz; ny++) {
          nx = i+3-nz-ny;

          /* I want to keep both divisions. The right hand side
             evaluates to an integer (although represented by a
             floating point number), so doing it this way makes
             the computation more accurate. This routine is only
             used for initialization anyway. */

          drdata[ir*m+i0+k] /= (fact[i+3]/(fact[nx]*fact[ny]*fact[nz]));
          k = k+1;
        }
    }
    for(i = 0; i<sizeof(p)/(3*sizeof(int)); i++)
      M [rveclist[ir][p[i][0]]]
        [rveclist[ir][p[i][1]]]
        [rveclist[ir][p[i][2]]] = ir*m;
  }
  M[0][0][0] = n+3;
}


static void dmsym3(int n,vector r,real8 *rvec/*,int npows[][3]*/)
{
  static real8 drdata[sizeof(rveclist)/(3*sizeof(int)) * 
		      (MAXORDER+3+3)*(MAXORDER+3+2)*(MAXORDER+3+1)/6];
  static int M[13][13][13] = {{{0}}};

  int i,j,m,jj[3],q[3],ridx,qi,qi0,nn[3],t;
  real8 rm,rmi,sgn[3][2],gn;
  vector ra;

/*
 *  The <M> and <drdata> arrays are static once initialized and shared
 *  among threads, but we have to make sure they are only initialized
 *  only once no matter how many threads are active.
 */
  if(M[0][0][0] < n+3) {
#ifdef _OPENMP
#pragma omp critical (CRIT_DMSYM3_INIT)
#endif
    {
      if(M[0][0][0] < n+3) {
        /*printf("len(drdata) = %d , sizeof(real8)=%d , sizeof(int)=%d\n",
          sizeof(drdata)/sizeof(real8),sizeof(real8),sizeof(int));*/
        if(sizeof(drdata)/sizeof(real8) < (n+3+3)*(n+3+2)*(n+3+1)/6)
          Fatal("@ %s:%d: drdata is too small!!",__FILE__,__LINE__);

        /* printf("Calling dmsym3init... sizeof(M)=%d  13*13*13=%d\n",
           (int) (sizeof(M)/sizeof(int)),13*13*13); */
        dmsym3init(n,M,drdata);
      }
    }
  }
  /* Collect signs and absolute values */
  rm = 0.0;
  for(i=0; i<3; i++) { /* I assume compiler unrolls this */
                       /* Replacing if-statements with max(), and fabs()
                          might be more efficient ? */
    sgn[i][0] = 1.0;
    if(r[i] >= 0.0) {
      ra[i] = r[i];
      sgn[i][1] = 1.0;
    } else {
      ra[i] = -r[i];
      sgn[i][1] = -1.0;
    }
    if(ra[i] > rm) rm = ra[i];
  }
  rmi = 12.0/rm;
  for(i = 0; i<3; i++)
    jj[i] = (int) (ra[i]*rmi + 0.1);

  for(i = 0; i<3; i++)
    if(fabs(jj[i] - ra[i]*rmi)>1e-10 || jj[i]>12 || jj[i]<0) {
      Fatal("dmsym vector error: %d,%d,%d  -- %15.5e,%15.5e,%15.5e",
	     jj[0],jj[1],jj[2],(double) ra[0],(double) ra[1],(double) ra[2]);
    }

  ridx = M[jj[0]][jj[1]][jj[2]];

  if(ridx < 0) {
    Fatal("Vector not in DM database, ridx=%d\n"
          "dmsym vector error: %d,%d,%d  -- %15.5e,%15.5e,%15.5e", ridx,
	  jj[0],jj[1],jj[2],(double) ra[0],(double) ra[1],(double) ra[2]);
  }

  /* Figure out permutation of coordinates */
  for(i = 0; i<3; i++) q[i] = i;
  for(i = 0; i<2; i++) /* This is a sorting routine, which should be
                          easily unrolled by the compiler */
    for(j = i+1; j<3; j++)
      if(jj[q[j]] > jj[q[i]]) {
        t = q[i];
        q[i] = q[j];
        q[j] = t;
      }

  /* Fill in dr vector using above permutation and inversion */

  /* This should be the most time consuming part by far. It may not
     be written in a smart way. By doing different versions for
     different n, all loops below can be totally unrolled, and all
     the integer arithmetic will disappear. I will make such a version
     and try it out in uBGL.

     This is essentially a copy/permute operation, where a vector
     from a table is copied to the output vector, but with different
     ordering of the element. Each element is also scaled and
     posibly has its sign changed.
     Below is suggested a way of optimizing the code by unrolling the
     entire loop below. A different approach might be to have six
     versions of the code below, one for each possibility of the
     permutation q. That way somee of the indexing, e.g. nn[q[2]] could
     be removed. The array nn could be replaced by three variables
     nx,ny,nx instead, which might make the life easier for the
     compiler.
  */

  m = 0;
  gn = rmi;
  for(i = 3; i<=n+3; i++) {
    gn *= rmi;
    qi0 = ( (i+2)*(i+1)*i - (3+2)*(3+1)*3 )/6;
    for(nn[2] = 0; nn[2]<=i; nn[2]++)
      for(nn[1] = 0; nn[1]<=i-nn[2]; nn[1]++) {
        nn[0] = i-nn[2]-nn[1];
        qi = nn[q[2]]*(2*i+3-nn[q[2]])/2 + nn[q[1]];

        rvec[m] = gn * sgn[0][nn[0]&1] * sgn[1][nn[1]&1] * sgn[2][nn[2]&1] *
          drdata[ridx+qi0+qi];

        /* Below is a generator for unrolled code. For the n used
           in the BGL simulations, the length of rvec is 276. Six
           versions of unrolled code would be needed, due to the
           different possibilities of the permutation q.

           A few sign/rescale variables, s0-s7, are needed by the
           unrolled and can be computed according to:
               for(i = 0; i<8; i++)
                 s_i = gn * sgn[0][(i/4)%2]*sgn[1][(i/2)%2]*sgn[2][i%2];

               k = 4*(nn[0]%2) + 2*(nn[1]%2) + (nn[0]%2);
               printf("rvec[%3d] = s%d*drdata[%3d];\n",m,k,ridx+qi0+qi);
        */

        m++;
      }
  }
}


/*---------------------------------------------------------------------------
 *
 *      Function:    makeqtab
 *      Description: Creates a table of quotients of the form
 *                   qtab[i][j] = i*(i-1)/(j+2)
 *
 *-------------------------------------------------------------------------*/
void makeqtab(real8 qtab[NMAX+1][NMAX+1])
{
        int   i, j;
        real8 cc;

        for (j = 0; j <= NMAX; j++) {
            qtab[0][j] = 1.0/(j+2);
        }

        for (i = 1; i <= NMAX; i++) {
            cc = i*(i-1);
            for (j = 0; j <= NMAX; j++) {
                qtab[i][j] = cc*qtab[0][j];
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    makeftabs
 *      Description: Builds tables of factorials fact[n] = n! = 1*2*3*...*n,
 *                   ifact[n] = 1/fact[n], and the double factorial
 *                   dfact[n] = n!! = 2*4*6*...*n if n even.
 *      Arguments:
 *          fact     Array of factorial values returned to caller. Must
 *                   be at least NMAX+1 elements in array!
 *          ifact    Array of inverse factorial values returned to caller.
 *                   Must be at least NMAX+1 elements in array!
 *          dfact    Array of double factorial values returned to caller. Must
 *                   be at least 2*NMAX-3+2 elements in array!
 *
 *-------------------------------------------------------------------------*/
void makeftabs(real8 *fact, real8 *ifact, real8 *dfact)
{
        int i;

/*
 *      Table for double factorial
 */
        dfact[0] = 1.0;

        for (i = 1; i <= 2*NMAX-3; i += 2) {
            dfact[i+1] = dfact[i-1]*i;
        }
  
/*
 *      Table for factorial and inverse factorial
 */
        fact[0] = 1.0;
        ifact[0] = 1.0;

        for (i = 1; i <= NMAX; i++) {
            fact[i] = fact[i-1]*i;
            ifact[i] = 1.0/fact[i];
        }

        return;
}


static void rdinit(int norder,int uorder,int maxorder,
		   int rdidx[],int pwvec[][3])
{
  int i,j,k,m,n,nn,nx,ny,nz,/*ix,*/iy,iz,idx,idx0;

  k = 0;
  idx0 = 0;
  for(i = 0; i<=uorder; i++) {
    for(nz = 0; nz<=i; nz++)
      for(ny = 0; ny<=i-nz; ny++) {
	nx = i-ny-nz;
		
	nn = norder+i;
	if(maxorder < nn) nn = maxorder;
	idx = idx0;

	for(j = i; j<=nn; j++) {
	  n = (j+2+3)*(j+1+3)/2;

	  for(iz = nz; iz<=j+3-nx; iz++)
	    for(iy = ny; iy<=j+3-iz-nx; iy++) {
	      /*ix = j+3-iz-iy;*/
	      m = iz*(2*(j+3) + 3 - iz)/2 + iy;
	      rdidx[k++] = idx+m;
	    }	
	  idx += n;
	}
      }
    idx0 += (i+2+3)*(i+1+3)/2;
  }

  idx0 = 0;
  for(i = 3; i<=norder+3; i++) {
    for(nz = 0; nz<=i; nz++)
      for(ny = 0; ny<=i-nz; ny++) {
	nx = i-nz-ny;
	pwvec[idx0][0] = nx;
	pwvec[idx0][1] = ny;
	pwvec[idx0][2] = nz;
	idx0 = idx0 + 1;
      }
  }
}


#ifdef ANISOTROPIC
/*---------------------------------------------------------------------------
 *
 *      Function:    MkTaylorAnisotropic
 *      Description: Version of MkTaylor for use with anisotropic elasticity.
 *                   See description of MkTaylor() for more details on
 *                   function and parameters.
 *
 *      Note: mu and nu are not used at all in this function, and any values
 *            can be passed. They are there just to make the calling sequence
 *            match the isotropic call.
 *-------------------------------------------------------------------------*/
void MkTaylorAnisotropic(Home_t *home, real8 mu, real8 nu, int norder,
                         int uorder, int maxorder, real8 r[], real8 eta[],
                         real8 *alpha)
{
    int            i, j, k;
    int            expc[3];
    int            kk[3], g;
    const int      neta = (norder+3)*(norder+2)*(norder+1)/6;
    real8          iscale, scale_factor;
    real8          ss[3], slen = 0.0;
    real8          fact[MAXORDER];
    real8          *etanew;
    real8          (*C)[3][3][3];
    real8          (*eta_tilde)[3][3];
    FMAnisoTable_t *anisoTable;
    static int     warnedOnOrderMismatch = 0;

    anisoTable = home->fmAnisoTable;
    C = (real8 (*)[3][3][3])home->anisoVars.elasticConstantMatrix4D;

/*
 *  Build factorial table
 */
    if (uorder+1 > MAXORDER) {
        Fatal("MkTaylorAnisotropic: uorder (%d) must be < than %d",
              uorder, MAXORDER);
    }

    fact[0] = 1.0;

    for (i = 1; i<=uorder; i++) {
        fact[i] = fact[i-1] * (double) i;
    }

/*
 *  Convert r to integral vector, and determine scaling factors.
 */
    for (i = 0; i < 3; i++) {
        ss[i] = fabs(r[i]/anisoTable->scaleVec[i]);
        if (ss[i] > slen) {
            slen = ss[i];
        }
    }

/*
 *  Deduce rational/integer coordinates from ss...
 */
    for (i = 0; i < 3; i++) {

        ss[i] *= (60.0/slen); /* All positive integers <= 6 divide 60... */
        kk[i] = (int) (ss[i] + 0.1);
    }

    for (i = 0; i < 3; i++) {
        if (kk[i] - ss[i] > 1e-5) {
            Fatal("MkTaylor_aniso: non-integral box vector\n"
                  "  {x,y,z}scale = (%15.5e , %15.5e , %15.5e)\n"
                  "  r-vector = (%15.5e , %15.5e , %15.5e)\n"
                  "  integer vector = (%18.8e , %18.8e , %18.8e)\n",
                  anisoTable->xScale, anisoTable->yScale, anisoTable->zScale,
                  r[0], r[1], r[2], ss[0], ss[1], ss[2]);
        }
    }

    g = gcd0(kk[0], gcd0(kk[1], kk[2]));

    if (g == 0) {
        Fatal("MkTaylor_aniso: zero returned from gcd(), kk = %d %d %d",
              kk[0],kk[1],kk[2]);
    }

/*
 *  expc will hold integer box translation vector / index into
 *  Greens function table
 */
    for (i = 0; i < 3; i++) {

        expc[i] = kk[i] / g;

        if (r[i] < 0.0) {
            expc[i] = -expc[i];
        }
    }

/*
 *  Check that Greens function derivative table contains expc
 */
    for (i = 0; i < 3; i++) {

        if (expc[i] < -5 || expc[i] > 5) {
            Fatal("MyTaylor_aniso: Integer translation vector outside of "
                  "bounds\nexpc = (%d,%d,%d) ; bounds are [-5,5] for each "
                  "coordinate.", expc[0], expc[1], expc[2]);
        }

        if (anisoTable->tagmap[expc[2]+5][expc[1]+5][expc[0]+5] == 0 ||
            anisoTable->Gderiv_map[expc[2]+5][expc[1]+5][expc[0]+5] == NULL) {
            Fatal("%s: Table missing data for integral translation vector\n"
                  "expc = (%d,%d,%d), tag = %d, ptr = %llx\n",
                  __func__, expc[0], expc[1], expc[2],
                  anisoTable->tagmap[expc[2]+5][expc[1]+5][expc[0]+5],
                  (unsigned long long int) anisoTable->Gderiv_map[expc[2]+5][expc[1]+5][expc[0]+5]);
        }
    }

/*
 *  Determine scaling...
 */
    for (j = 0, i = 1; i < 3; i++) {
        if (fabs(r[i]) > fabs(r[j])) {
            j = i;
        }
    }

    iscale = (expc[j]*anisoTable->scaleVec[j]) / r[j];

/*
 *  Compute a rescaled version of the multipole coefficients
 */
    etanew = (real8 *)malloc(sizeof(real8) * 9 * (norder+3) *
                             (norder+2) * (norder+1) / 6);
    {
        real8 scale2 = iscale;
        int k = 0;

        for (i = 0; i <= norder; i++) {
            int nz, ny, j;

            for (j = 0; j<9; j++) {
                for (nz = 0; nz <= i; nz++) {
                    for (ny = 0; ny <= i-nz; ny++) {
                        const int nx = i-nz-ny;

                        etanew[k] = eta[k] * scale2;
                        k++;
                    }
                }
            }
            scale2 *= iscale;
        }
    }

/*
 *  Convolve eta with C, to form:
 *      eta_tilde(a,b,c,p,g,q) = eps(n,g,r)*C(p,q,w,n)*eta(a,b,c,w,r)/(a!b!c!)
 */
    eta_tilde = (real8 (*)[3][3]) malloc(27 * neta * sizeof(real8));
    eta2eta_tilde(norder, C, etanew, eta_tilde);

/*
 *  Create Taylor expansion of Greens function
 */
    scale_factor = -1.0;

    for (k = 0, i = 0; i <=uorder; i++) {
        int nz,ny;
        int nn = norder+i;

        if (maxorder < nn) {
            nn = maxorder;
        }

        if (nn+1 > anisoTable->derivTblMaxOrder) {
            if (home->myDomain == 0) {
                if (warnedOnOrderMismatch == 0) {
#ifdef _OPENMP
#pragma omp critical (CRIT_WARN_ON_ORDER_MISMATCH)
#endif
                    {
                        if (warnedOnOrderMismatch == 0) {
#if 0
			    fprintf(stderr, "Warning in %s(): Requested "
                                    "norder = %d but Greens function "
                                    "derivative table only has up to "
                                    "order = %d. Truncating expansion. "
                                    "Loss of accuracy may occur.\n",
                                    __func__, nn, anisoTable->derivTblMaxOrder);
#endif

                            warnedOnOrderMismatch = 1;
                        }

                    }  /* end omp critical section */
                }
            }

            nn = anisoTable->derivTblMaxOrder - 1;
        }

/*
 *      In loop iterations where nn-i goes negative, we don't 
 *      don't want to do anything more, so just continue.
 */
        if (nn - i < 0) {
            continue;
        }

        scale_factor *= -iscale;

        for (nz = 0; nz <= i; nz++) {
            for (ny = 0; ny <= i-nz; ny++) {
                int         p, q;
                real8       dsigma[3][3];
                const int   nx = i-ny-nz;
                const real8 g = 1.0 / (fact[nx]*fact[ny]*fact[nz]);

                G2sigma_tilde_deriv(
                        norder, nn-i, nx, ny, nz,
                        anisoTable->Gderiv_map[expc[2]+5][expc[1]+5][expc[0]+5],
                        eta_tilde, dsigma);

                for (p = 0; p < 3; p++) {
                    for (q = 0; q < 3; q++) {
                        alpha[k*9+p*3+q] = g*dsigma[p][q] * scale_factor;
                    }
                }

                k++;
            }
        }
    }

    free(etanew);
    free(eta_tilde);

    return;
}
#endif  /* ifdef ANISOTROPIC */


/*---------------------------------------------------------------------------
 *
 *      Function:    MkTaylorIsotropic
 *      Description: Version of MkTaylor for use with isotropic elasticity.
 *                   See description of MkTaylor() for more details on
 *                   function and parameters.
 *
 *-------------------------------------------------------------------------*/
void MkTaylorIsotropic(real8 mu, real8 nu, int norder, int uorder,
                       int maxorder, real8 r[], real8 eta[], real8 *alpha)
{
  real8 rderiv[(MAXORDER+3)*(MAXORDER+2)*(MAXORDER+1)/6];
  real8  rdvec[(NMAX+3)*(NMAX+2)*(NMAX+1)/6];

  int i,j,k,m,p,q,/*n,idx,idx0,*/nx,ny,nz,k2,nn;
  real8 g;
  matrix dsigma;

  static int save_norder = -1,save_uorder = -1,save_maxorder = -1;
  static real8 fact[MAXORDER+4];
  static int rdidx[(NMAX+3)*(NMAX+2)*(NMAX+1)/6*(NMAX+3)*(NMAX+2)*(NMAX+1)/6];
  static int pwvec[(NMAX+3)*(NMAX+2)*(NMAX+1)/6][3];

  /* Table of factorials */
  if(maxorder > MAXORDER) {
    Fatal("MkTaylor: maxorder (%d) must be <= than %d",
          maxorder, MAXORDER);
  }

/*
 *  Various arrays are static once initialized and shared among threads,
 *  but we have to make sure they are only initialized only once no matter
 *  how many threads are active.
 */
  if(save_norder!=norder || save_uorder!=uorder || save_maxorder!=maxorder) {
#ifdef _OPENMP
#pragma omp critical (CRIT_INIT_MKTAYLOR)
#endif
    {
      if ((save_norder != norder) ||
          (save_uorder != uorder) ||
          (save_maxorder != maxorder)) {
        fact[0] = 1.0;
        for(i = 1; i<=maxorder+3; i++) fact[i] = fact[i-1] * (real8) i;
        rdinit(norder,uorder,maxorder,rdidx,pwvec);
        save_norder = norder;
        save_uorder = uorder;
        save_maxorder = maxorder;
      }
    }
  }


  /* Compute dervatives of r */
#if 1
  dmsym3(maxorder,r,rderiv);
#else
/*
 *  This block enables non-cubic FMM-cells. The previous code consists
 *  solely of the call to dmsym3, now commented out. This is an ugly
 *  hack, and will impact the FMM speed; I am guessing it will be
 *  30%-50% slower. -- Tomas
 */
  {
    static int powvec[(MAXORDER+3+3)*(MAXORDER+3+2)*(MAXORDER+3+1)/6][3];
    DM(maxorder,r,rderiv,powvec);
  }
#endif
  
  /*printf("Staring main computation loop\n");*/
  k = 0;
  k2 = 0;
  /*idx0 = 0;*/
  m = ((norder+3+3)*(norder+3+2)*(norder+3+1)-(2+3)*(2+2)*(2+1))/6;

  for(i = 0; i<=uorder; i++) {
    nn = norder+i;
    if(maxorder < nn) {
      nn = maxorder;
      m = ((nn-i+3+3)*(nn-i+3+2)*(nn-i+3+1)-(2+3)*(2+2)*(2+1))/6;
    }

    for(nz = 0; nz<=i; nz++)
      for(ny = 0; ny<=i-nz; ny++) {
	nx = i-ny-nz;
	
	for(j = 0; j<m; j++) {
	  rdvec[j] = rderiv[rdidx[k2++]];
	}

	/*
	k2 = 0;
	idx = idx0;
	for(j = i; j<=nn; j++) {
	  int ix,iy,iz;
	  n = (j+2+3)*(j+1+3)/2;
	  / *idx = n*i/3;* /

	  for(iz = nz; iz<=j+3-nx; iz++)
	    for(iy = ny; iy<=j+3-iz-nx; iy++) {
	      ix = j+3-iz-iy;
	      m = iz*(2*(j+3) + 3 - iz)/2 + iy;

	      rdvec[k2] = rderiv[idx+m];
	      k2++;
	    }	
	  idx += n;
	}
	*/

	FMSigma2(mu,nu,nn-i,eta,dsigma,rdvec,pwvec);
	
	g = 1.0 / (fact[nx]*fact[ny]*fact[nz]);
	for(p = 0; p<3; p++)
	  for(q = 0; q<3; q++)
            alpha[k*9+p*3+q] = g*dsigma[p][q];
	k++;
      }
    /*idx0 += (i+2+3)*(i+1+3)/2;*/
  }
}

/*---------------------------------------------------------------------------
 *
 *      Function:    MkTaylor
 *      Description: Calculate the taylor expansion contribution
 *                   for a cell from the multipole expansion of
 *                   a remote cell.
 *
 *      Arguments:
 *          MU       Shear modulus
 *          NU       Poisson ratio
 *          nOrder   order of the multipole expansion
 *          mOrder   order of the taylor expansion
 *          maxOrder
 *          r        vector from the center of taylor expansion
 *                   to the center of the multipole expansion
 *          eta      Array of multipole expansion coefficients.  Array
 *                   should have 9*(nOrder+1)*(nOrder+2)*(nOrder+3)/6
 *                   elements 
 *          alpha    Array in which the taylor expansion coefficients
 *                   will be returned to the caller.  Array should have 
 *                   9*(mOrder+1)*(mOrder+2)*(mOrder+3)/6 elements 
 *
 *-------------------------------------------------------------------------*/
void MkTaylor(Home_t *home, real8 mu, real8 nu, int norder, int uorder,
              int maxorder, real8 r[], real8 eta[], real8 *alpha)
{

/*
 *  Note: When general anisotropy is enabled via the ANISOTROPIC
 *        definition, the user may specify that the FMM code use
 *        either isotropy or anisotropy...
 */

#ifdef ANISOTROPIC
    if (home->param->fmEnableAnisotropy) {
/*
 *      WARNING!
 *
 *      The newer anisotropic version of MkTaylor was written with the
 *      opposite sign convention for the input vector <r>, so we have
 *      to negate the vector when using the anisotropic function.
 */
        real8 rNeg[3];

        rNeg[0] = -r[0];
        rNeg[1] = -r[1];
        rNeg[2] = -r[2];

        MkTaylorAnisotropic(home, mu, nu, norder, uorder, maxorder, rNeg,
                            eta, alpha);
        return;
    }
#endif

/*
 *  Either general anisotropy is not enabled, or the user does not
 *  want to use anisotropic FMM...
 */
    MkTaylorIsotropic(mu, nu, norder, uorder, maxorder, r, eta, alpha);

    return;
}



/*---------------------------------------------------------------------------
 *
 *      Function:    makeeta
 *      Description: Create a multipole expansion representing the stress
 *                   field from a straigth line. The expansion is of order
 *                   norder, and the resulting coefficients are stored in
 *                   the array Eeta, which is of length
 *                   9*norder*(norder+1)*(norder+2)/6.
 *
 *      Arguments:
 *          ksi0    vector from the center of expansion, to the begining
 *                  of the line.
 *          ksi     vector from the beginning of the line to the end of
 *                  the line.
 *          bn      Burgers vector.
 *
 *      To create a multipole expansion from several lines using
 *      the same expansion center, the coefficients from multiple calls
 *      to makeeta can be added, e.g.:
 *
 *          sz = 9*norder*(norder+1)*(norder+2)/6
 *          Eetatot(1..sz) = 0
 *          for i=1..nsegments
 *            ksi0 = startpoint(i) - center
 *            ksi  =   endpoint(i) - startpoint(i)
 *            makeeta(norderm,ksi0,ksi,b(i),Eetatemp)
 *            for j=1..sz
 *              Eetatot(j) = Eetatot(j) + Eetatemp(j)
 *            end
 *          end
 *
 *-------------------------------------------------------------------------*/
void makeeta(int norder, real8 ksi0[3], real8 ksi[3], real8 b[3], real8 *Eeta) {
        int   i, j, k, i0, i1, etaidx;
        int   n, nx, ny, nz, ngq, ixx, jfirst;
        real8 t, x, y, z;
        real8 ksinorm, ksinormi;
        real8 ksihat[3], E[9];
        real8 xp[NMAX+1][NMAX/2+1];
        real8 yp[NMAX+1][NMAX/2+1];
        real8 zp[NMAX+1][NMAX/2+1];

        static int   initlevel = 0;
        static real8 xx[(NMAX/2+1)*(NMAX/2+1+1)/2];
        static real8 ww[(NMAX/2+1)*(NMAX/2+1+1)/2];

        ngq = (norder >> 1) + 1;

        if (initlevel < ngq) {
            j = initlevel * (initlevel+1) >> 1;

            for (i = initlevel+1; i <= ngq; i++) {
                GaussQuadraturePointsAndWgts(i,&xx[j],&ww[j]);

                for (k = 0; k < i; k++) {
                    xx[j+k] = 0.5*xx[j+k] + 0.5;
                    ww[j+k] = 0.5*ww[j+k];
                }
                j = j + i;
            }
            initlevel = ngq;
        }

/*
 *      It should be a rare occurence, but just in case this function
 *      is called for a zero length segment, make sure we don't do any
 *      divide-by-zero.
 */
        ksinorm = sqrt(ksi[0]*ksi[0] + ksi[1]*ksi[1] + ksi[2]*ksi[2]);

        if (ksinorm > 0.0) {
            ksinormi = 1.0 / ksinorm;
        } else {
            ksinorm  = 0.0;
            ksinormi = 0.0;
        }

        ksihat[0] = ksi[0] * ksinormi;
        ksihat[1] = ksi[1] * ksinormi;
        ksihat[2] = ksi[2] * ksinormi;

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                E[3*i+j] = b[i] * ksihat[j];
            }
        }

        etaidx = 0;

        for (n = 0; n <= norder; n++) {
            i1 = (n+1)*(n+2) >> 1;
            i0 = 3*n*i1;
    
/*
 *         More heavy duty inlining...
 *         doeta(i,ksi0,ksi,&Eeta[i0+8*i1]);
 */
    
/*
 *          Gauss quadrature initialization
 */
            ngq = (n >> 1) + 1;
            ixx = ngq * (ngq-1) >> 1;

            if (n & 1) jfirst = n; /* ngq same as last time! */
            else jfirst = 1;

            for (i = 0; i < ngq; i++) {
                t = xx[ixx+i];
/*
 *              Power table again, no ipow
 */
                x = ksi0[0] + t*ksi[0];
                y = ksi0[1] + t*ksi[1];
                z = ksi0[2] + t*ksi[2];

                xp[0][i] = 1.0;
                yp[0][i] = 1.0;
                zp[0][i] = 1.0;

                for (j = jfirst; j <= n; j++) {
                    xp[j][i] = x * xp[j-1][i];
                    yp[j][i] = y * yp[j-1][i];
                    zp[j][i] = z * zp[j-1][i];
                }
            }
    
            i = i0 + 8*i1;

            for (nz=0; nz <= n; nz++) {
                for (ny = 0; ny <= n-nz; ny++) {
                    nx = n-nz-ny;
                    t = 0.0;
                    for (j = 0; j < ngq; j++)
                        t = t + ww[ixx+j] * xp[nx][j] * yp[ny][j] * zp[nz][j];
                    Eeta[i++] = t*ksinorm;
                }
           } /* End of inline */

        
            for (j = 0; j < 9; j++) {
                for (k = 0; k < i1; k++) {
                    /* etaidx == i0+j*i1+k */
                    Eeta[etaidx++] = E[j] * Eeta[i0+8*i1+k];
                }
            }
        }  /* for (n = 0; ...) */

        return;
}



/*---------------------------------------------------------------------------
 *
 *      Function:    FMShift
 *      Description: Shift the center of the provided multipole expansion.
 *                   Used during the upward pass of the FM code.
 *
 *      Arguments:
 *          norder   order of the multipole expansion
 *          r        vector from the center of the original expansion
 *                   to the center of the new expansion.
 *          eta      Array of multipole expansion coefficients.  Array
 *                   should have 9*(norder+1)*(norder+2)*(norder+3)/6
 *                   elements 
 *          neta     Array in which the shifted multipole expansion
 *                   coefficients are returned to the caller.  Sized the
 *                   same as eta array.
 *
 *-------------------------------------------------------------------------*/
void FMShift(int norder, real8 *r, real8 *eta, real8 *neta)
{
        int   i, j, k, nx, ny, nz, a, b, c, idx, n;
        real8 m, rp;
        real8 pw[NMAX+1][3], g[NMAX+1];

        static real8 A[NMAX+1][NMAX+1][NMAX+1][9];

/*
 *      FIX ME!  array <g> can be calculated one time and
 *      saved for future use since norder never varies during any 
 *      given execution.
 *
 *      g[] = array of factorials.
 */
        g[0] = 1.0;

        for (j = 0; j < 3; j++) pw[0][j] = 1.0;

        for (i = 1; i <= norder; i++) {
            g[i] = g[i-1] * (real8) i;
            for (j = 0; j < 3; j++) pw[i][j] = pw[i-1][j] * r[j];
        }

        for (i = 0; i<=norder; i++) {
            k = 0;
            n = (i+2) * (i+1) / 2;
            idx = 3 * n * i;
            for (nz = 0; nz <= i; nz++) {
                for (ny = 0; ny <= i-nz; ny++) {
                    nx = i-ny-nz;
                    for (j = 0; j < 9; j++) {
                        A[nz][ny][nx][j] = eta[idx+n*j+k];
                    }
                    k = k+1;
                }
            }
        }


        for (i = 0; i <= norder; i++) {
            k = 0;
            n = (i+2) * (i+1) / 2;
            idx = 3*n*i;

            for (nz = 0; nz <= i; nz++) {
                for (ny = 0; ny <= i-nz; ny++) {
                    nx = i-ny-nz;

                    for (j = 0; j < 9; j++) {
                        neta[idx+n*j+k] = 0.0;
                    }

                    for (c = 0; c <= nz; c++) {
                        for (b = 0; b <= ny; b++)
                            for (a = 0; a <= nx; a++) {
                                m = g[nx]/(g[a]*g[nx-a]) *
                                    g[ny]/(g[b]*g[ny-b]) *
                                    g[nz]/(g[c]*g[nz-c]);
                                rp = pw[a][0] * pw[b][1] * pw[c][2];

                                for (j=0; j < 9; j++) {
                                    neta[idx+n*j+k] += m*rp*A[nz-c][ny-b][nx-a][j];
                                }
                            }  /* for (a = 0; ...) */
                    }  /* for (c = 0; ...) */
                    k = k+1;
                }  /* for (ny = 0; ...) */
            }  /* for (nz = 0; ...) */
        }  /* for (i = 0; ...) */

       return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    TaylorShift
 *      Description: Shift the center of the provided taylor expansion.
 *                   Used during the downward pass of the FM code.
 *
 *      Arguments:
 *          norder   order of the taylor expansion
 *          r        vector from the center of the original expansion
 *                   to the center of the new expansion.
 *          alpha    Array of taylor expansion coefficients.  Array
 *                   should have 9*(norder+1)*(norder+2)*(norder+3)/6
 *                   elements 
 *          beta     Array in which the shifted taylor expansion
 *                   coefficients are returned to the caller.  Sized the
 *                   same as alpha array.
 *
 *-------------------------------------------------------------------------*/
void TaylorShift(int norder, real8 *r, real8 *alpha, real8 *beta)
{
        int   i, j, j2, k, nx, ny, nz, a, b, c;
        real8 pw[NMAX+1][3];
        real8 g[NMAX+1], m, rp;

        static real8 A[NMAX+1][NMAX+1][NMAX+1][3][3];

        g[0] = 1.0;

        for (j = 0; j < 3; j++) pw[0][j] = 1.0;

        for (i = 1; i <= norder; i++) {
            g[i] = g[i-1] * (real8) i;
            for (j = 0; j < 3; j++) pw[i][j] = pw[i-1][j] * r[j];
        }
  
        k = 0;

        for (i = 0; i <= norder; i++) {
            for (nz = 0; nz <= i; nz++) {
                for (ny = 0; ny <= i-nz; ny++) {
                    nx = i-ny-nz;
                    for (j = 0; j < 3; j++) {
                        for (j2 = 0; j2 < 3; j2++) {
                            A[nz][ny][nx][j][j2] = alpha[k*9+j*3+j2];
                        }
                    }
                    k = k+1;
                }
            }
        }
 
        k = 0;

        for (i = 0; i <= norder; i++) {
            for (nz = 0; nz <= i; nz++) {
                for (ny = 0; ny <= i-nz; ny++) {
                    nx = i-ny-nz;
                    for (j = 0; j < 3; j++) {
                        for (j2=0; j2 < 3; j2++) {
                            beta[k*9+j*3+j2] = 0.0;
                        }
                    }

                    for (c = nz; c <= norder; c++) {
                        for (b = ny; b <= norder-c; b++) {
                            for (a = nx; a <= norder-b-c; a++) {
                                m = g[a]/(g[nx]*g[a-nx]) *
                                    g[b]/(g[ny]*g[b-ny]) *
                                    g[c]/(g[nz]*g[c-nz]);
                                rp = pw[a-nx][0] * pw[b-ny][1] * pw[c-nz][2];
                                for (j = 0; j < 3; j++) {
                                    for (j2 = 0; j2 < 3; j2++) {
                                        beta[k*9+j*3+j2] +=
                                                m*rp*A[c][b][a][j][j2];
                                    }
                                }
                            }  /* for (a = nx; a <= norder-b-c ...) */
                        }  /* for (b = ny; b <= norder-c ...) */
                    }  /* for (c = nz; c <= norder ...) */
                    k = k+1;
                }  /* for (ny = 0; ny <= i-nz ...) */
            }  /* for (nz = 0; nz <= i ...) */
        }  /* for (i = 0; i <= norder ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    EvalTaylor
 *      Description: Compute stress at a point from a taylor expansion
 *
 *      Arguments:
 *          uorder   order of the taylor expansion
 *          r        vector from the center of taylor expansion
 *                   to the point at which stress is to be evaluated
 *          alpha    Array of taylor expansion coefficients.  Array
 *                   should have 9*(norder+1)*(norder+2)*(norder+3)/6
 *                   elements 
 *          sigma    3X3 matrix in which the calculated stress is returned
 *                   to the caller.
 *
 *-------------------------------------------------------------------------*/
void EvalTaylor(int uorder, real8 *r, real8 *alpha, real8 sigma[3][3])
{
        int   i, j, k, m, nx, ny, nz;
        real8 rp;
        real8 pw[NMAX+1][3];

        for (j = 0; j < 3; j++) pw[0][j] = 1.0;

        for (i = 1; i <= uorder; i++) {
            for (j = 0; j < 3; j++) {
                pw[i][j] = pw[i-1][j] * r[j];
            }
        }

        for (j = 0; j < 3; j++) {
            for (m = 0; m < 3; m++) {
                sigma[j][m] = 0.0;
            }
        }

        k = 0;

        for (i = 0; i <= uorder; i++) {
            for (nz = 0; nz <= i; nz++) {
                for (ny = 0; ny <= i-nz; ny++) {
                    nx = i-ny-nz;
                    rp = pw[nx][0] * pw[ny][1] * pw[nz][2];
                    for (j = 0; j < 3; j++) {
                        for (m = 0; m < 3; m++) {
                            sigma[j][m] += rp * alpha[k*9+j*3+m];
                        }
                    }
                    k = k+1;
                }
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Author:      Sylvie Aubry
 *      Function:    EvalMeanTaylor
 *      Description: Compute the average stress from a taylor expansion 
 *                   in the simulation cell to account for conditionally
 *                   convergence problem.
 *
 *      Arguments:
 *          avgStress  3X3 matrix in which the calculated average stress 
 *                     is returned to the caller.
 *          avgStressEshleby  3X3 matrix in which the calculated average
 *                     stress from eshelby inclusions is returned to the
 *                     caller.
 *
 *-------------------------------------------------------------------------*/
static void EvalMeanTaylor(Home_t *home, real8 avgStress[3][3], real8 avgStressEshelby[3][3])
{
        int       i, j, k, m, n, nx, ny, nz;
        int       numRows = 3;
        real8     rp;
        real8     pw[NMAX+1][3];
        real8     Min[3], Max[3], sigma[6][3], stress[6][3];

        Param_t   *param = home->param;
        FMLayer_t *layer = &home->fmLayer[param->fmNumLayers-1];

        /* Cell size bounds*/
        Min[0] = param->minSideX/layer->lDim[X];
        Max[0] = param->maxSideX/layer->lDim[X];

        Min[1] = param->minSideY/layer->lDim[Y];
        Max[1] = param->maxSideY/layer->lDim[Y];
        
        Min[2] = param->minSideZ/layer->lDim[Z];
        Max[2] = param->maxSideZ/layer->lDim[Z];

        int   Ncells = layer->lDim[X]*layer->lDim[Y]*layer->lDim[Z];
        real8 Volume = (Max[0]- Min[0])*(Max[1]- Min[1])*(Max[2]- Min[2])*Ncells;

#if 0 //def ESHELBY
        if (param->eshelbyfmEnabled)
            numRows += 3;
#endif

        for (j = 0; j < numRows; j++) {
            for (m = 0; m < 3; m++) {
                stress[j][m] = 0.0;
            }
        }

/*
 *      Compute int_Min ^ Max x^n:  First for standard FM, then for
 *      the FM for eshelby inclusions (if necessary)
 */
        int uorder = param->fmExpansionOrder;

        for (i = 0; i <= uorder; i++) {
            for (j = 0; j < 3; j++) {
                pw[i][j] = (pow(Max[j],i+1)-pow(Min[j],i+1))/((i+1)*1.0);
            }
        }

#if 0 //def ESHELBY
        int   uorderEshelby = param->eshelbyfmTaylorOrder;
        real8 pwEshelby[NMAX+1][3];

        for (i = 0; i <= uorderEshelby; i++) {
            for (j = 0; j < 3; j++) {
                pwEshelby[i][j] = (pow(Max[j],i+1)-pow(Min[j],i+1))/((i+1)*1.0);
            }
        }
#endif

/*
 *      Sum up the stress for all cells at the layer
 */
        for (n = 0; n < layer->ownedCnt; n++) 
        {
            int        cellID = layer->ownedMin + n;
            FMCell_t  *cell   = LookupFMCell(layer->cellTable, cellID);

            for (j = 0; j < 3; j++) {
                for (m = 0; m < 3; m++) {
                    sigma[j][m] = 0.0;
                }
            }

            k = 0;

            for (i = 0; i <= uorder; i++) {
                for (nz = 0; nz <= i; nz++) {
                    for (ny = 0; ny <= i-nz; ny++) {
                        nx = i-ny-nz;
                        rp = pw[nx][0] * pw[ny][1] * pw[nz][2];
                        for (j = 0; j < 3; j++) {
                            for (m = 0; m < 3; m++) {
                                sigma[j][m] += rp * cell->expansionCoeff[k*9+j*3+m];
                            }
                        }
                        k = k+1;
                    }
                }
            }

#if 0 //def ESHELBY
/*
 *          If FM is enabled for eshelby inclusions, also sum
 *          up stress from inclusions for all cells, but store it
 *          separately from the other sum.
 */
            if (param->eshelbyfmEnabled) 
            {
                for (j = 3; j < 6; j++) {
                    for (m = 0; m < 3; m++) {
                        sigma[j][m] = 0.0;
                    }
                }

                k = 0;

                for (i = 0; i <= uorderEshelby; i++) {
                    for (nz = 0; nz <= i; nz++) {
                        for (ny = 0; ny <= i-nz; ny++) {
                            nx = i-ny-nz;
                            rp = pwEshelby[nx][0] *
                                 pwEshelby[ny][1] *
                                 pwEshelby[nz][2];
                            for (j = 0; j < 3; j++) {
                                for (m = 0; m < 3; m++) {
                                    sigma[j+3][m] += rp * cell->eshelbyexpansionCoeff[k*9+j*3+m];
                                }
                            }
                            k = k+1;
                        }
                    }
                }
            }
#endif  // ESHELBY

            for (j = 0; j < numRows; j++) {
                for (m = 0; m < 3; m++) {
                    stress[j][m] += sigma[j][m];
                }
            }

        }  /* end for (n = 0; ... ) */

#ifdef PARALLEL
        real8  avgStressGlobal[6][3];

        MPI_Allreduce((real8 *) stress, (real8 *) avgStressGlobal, numRows * 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for (j = 0; j < 3; j++) {
            for (m = 0; m < 3; m++) {
                avgStress[j][m] = avgStressGlobal[j][m] / Volume;

#if 0//def ESHELBY
                if (param->eshelbyfmEnabled)
                    avgStressEshelby[j][m] = avgStressGlobal[j+3][m] / Volume;
#endif
            }
        }

#else
        for (m = 0; m < 3; m++) {
            for (k = 0; k < 3; k++) {
                avgStress[m][k] = stress[m][k]/Volume;

#if 0//def ESHELBY
                if (param->eshelbyfmEnabled) 
                    avgStressEshelby[m][k] = stress[m+3][k] / Volume;
#endif
            }
        }
#endif

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Author:      Sylvie Aubry
 *      Function:    MeanStressCorrection
 *      Description: Subtract the average stress computed in EvalMeanTaylor()
 *                   from the first order coefficients of the taylor
 *                   coefficients of each cell this domain intersects.
 *
 *-------------------------------------------------------------------------*/
void Taylor_MeanStressCorrection(Home_t *home)
{
        real8     avgStress[3][3], avgStressEshelby[3][3];
        Param_t   *param;
        FMLayer_t *layer;
  
        param = home->param;
  
        layer = &home->fmLayer[param->fmNumLayers-1];

        EvalMeanTaylor(home, avgStress, avgStressEshelby);

        for (int i = 0; i < layer->ownedCnt; i++) {

            int        cellID = layer->ownedMin + i;
            FMCell_t  *cell   = LookupFMCell(layer->cellTable, cellID);

            for (int j = 0; j < 3; j++) {
                for (int m = 0; m < 3; m++) {
                    cell->expansionCoeff[j*3+m] -= avgStress[j][m];

#if 0  // SA : FIX ME def ESHELBY
                    if (param->eshelbyfmEnabled) 
                        cell->eshelbyexpansionCoeff[j*3+m] -= avgStressEshelby[j][m];
#endif
                }
            }

        }  /* end for (i = 0; ... ) */

        return;
}
#endif


#ifdef CALCENERGY
// Energy corrections are specific to the Taylor FMM method.

static void DM2(int n, double *rvec, double rc, double *terms, int npows[][3])
{
        int   i, nx, ny, nz, a, b, c, m, p;
        double dfactx[2*MAXORDER-3+2 + 2], fact[MAXORDER+1], gtab[MAXORDER+1][3];
        double mparity, f, ca, cb, cc, rfact, gx, gy, gz, r, x, y, z, xn, yn, zn;
        double *dfact = dfactx+2;
/*
 *      Compute table for double factorial
 */
        dfact[-2] = -1.0;
        dfact[0] = 1.0;

        for (i = 1; i <= 2*n-3; i += 2) {
            dfact[i+1] = dfact[i-1]*i;
        }
  
/*
 *      Compute table for factorial
 */
        fact[0] = 1.0;

        for (i = 1; i <= n; i++) {
            fact[i] = fact[i-1]*i;
        }

        r = 0.0;

        for (i = 0; i < 3; i++) {
           r = r + rvec[i]*rvec[i];// + rc * rc;
        }

        r = sqrt(r);
        rfact = ipow(-1.0/r, n-1);

        gx = rvec[0] / r;
        gy = rvec[1] / r;
        gz = rvec[2] / r;

/*
 *      This might be overkill, but has better accuracy than
 *          x(n) = gx*x(n-1)
 *
 *      This is just an unrolled version of ipow()
 */
        for (i = 0; i <= n; i++) {
            p = i;

            xn = gx;
            yn = gy;
            zn = gz;

            x = y = z = 1.0;

            while (p > 0) {

                if (p & 1) {
                    x = x*xn;
                    y = y*yn;
                    z = z*zn;
                }

                xn = xn*xn;
                yn = yn*yn;
                zn = zn*zn;
                p = p >> 1;
            }
            gtab[i][0] = x;
            gtab[i][1] = y;
            gtab[i][2] = z;
        }

        i = 0;

        for (nz = 0; nz <= n; nz++) {
            for (ny = 0; ny <= n-nz; ny++) {
                nx = n-nz-ny;

                npows[i][0] = nx;
                npows[i][1] = ny;
                npows[i][2] = nz;

                terms[i] = 0;
                ca = 1.0;

                for (a = 0; a <= nx; a += 2) {
                    cb = 1.0;
                    for (b = 0; b <= ny; b += 2) {
                        cc = 1.0;

                        m = (a+b)/2;
                        mparity = 1-2*(m%2);

                        for (c = 0; c <= nz; c += 2) {
                            f = mparity * dfact[2*(n-m)-3+1] * 
                                gtab[nx-a][0] *
                                gtab[ny-b][1] *
                                gtab[nz-c][2];

                            f = f * ca * cb * cc;
                            terms[i] = terms[i] + f;
                            cc = cc * (nz-c) * (nz-c-1) / (c+2);
                            m = m + 1;
                            mparity = -mparity;
                        }
                        cb = cb * (ny-b) * (ny-b-1) / (b+2);
                    }
                    ca = ca * (nx-a) * (nx-a-1) / (a+2);
                }
                terms[i] = terms[i] * fact[n] /
                           (fact[nx]*fact[ny]*fact[nz]) * rfact;

                i = i+1;
            }
        }

        return;
}



static int ind(int nx,int ny,int nz)
{
   const int i = nx+ny+nz;
   const int idx0 = (i+2)*(i+1)*i/6;
   const int idx  = nz*(i+1) + ny - (nz-1)*nz/2;
   return idx0 + idx;
}

#endif
