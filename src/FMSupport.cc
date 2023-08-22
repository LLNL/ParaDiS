/**************************************************************************
 *
 *      Module:  FMSupport.c 
 *
 *      Includes functions:
 *          EncodeFMCellIndex()
 *          GaussQuadraturePointsAndWgts()
 *          IDToFMCellIndices()
 *
 *      NOTE: Functions in this file are used for all FMM methods.
 *            Most of these functions are very nearly copied from code put
 *            together by Tomas Oppelstrup. 
 *
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"

#ifdef ANISOTROPIC
#include "G2sigma_tilde.h"
#endif

/*
static int rveclist[][3] = { {6,0,0} , {6,3,0} , {6,6,0} , {6,3,3} ,
                             {6,6,3} , {6,6,6} , {6,2,0} , {6,4,0} ,
                             {6,2,2} , {6,4,2} , {6,6,2} , {6,4,4} ,
                             {6,6,4} };
*/

/*
static int rveclist[][3] = { {12, 0, 0} , {12,12, 0} , {12,12,12} ,
                             {12, 6, 0} , {12, 6, 6} , {12,12, 6} ,
                             {12, 4, 0} , {12, 4, 4} , {12, 8, 0} ,
                             {12, 8, 4} , {12, 8, 8} , {12, 12,4} ,
                             {12,12, 8} , {12, 3, 0} , {12, 3, 3} ,
                             {12, 6, 3} , {12, 9, 0} , {12, 9, 3} ,
                             {12, 9, 6} , {12, 9, 9} , {12,12, 3} ,
                             {12,12, 9} };
*/

/*---------------------------------------------------------------------------
 *
 *      Function:     IDToFMCellIndices
 *      Description:  Assuming a 3D block of cells with the dimensions
 *                    given by <cellsPerDim>, calculate the (x,y,z) 
 *                    coordinates in the block for the cell at position <id>
 *                    in the corresponding 1D array of the cells in the block.
 *                    Assumption is that indices change fastest in Z-direction
 *                    and slowest in X.
 *      Arguments:
 *          IN:   id          Index of cell in 1D representation of the 
 *                            cell block defined by <cellsPerDim>
 *          IN:   cellsPerDim Number of cells in each dimension of the
 *                            block of cells encompassing the cell in
 *                            question.
 *          OUT:  indices     Set to the x,y,z indices of the cell within
 *                            the block of cells indicated by <cellsPerDim>.
 *
 *--------------------------------------------------------------------------*/
void IDToFMCellIndices(int id, int cellsPerDim[3], int indices[3])
{
        int totalCount, perPlaneCount, perRowCount;

        totalCount = cellsPerDim[X] * cellsPerDim[Y] * cellsPerDim[Z];
        perPlaneCount = cellsPerDim[Y] * cellsPerDim[Z];
        perRowCount = cellsPerDim[Z];

        indices[X] = id / perPlaneCount;
        indices[Y] = (id - indices[X] * perPlaneCount) / perRowCount;
        indices[Z] = id - ((indices[X] * perPlaneCount) +
                           (indices[Y] * perRowCount));

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     EncodeFMCellIndex
 *      Description:  Assuming a 3D block of the dimensions
 *                    given by <dim>, calculate the index
 *                    in a corresponding 1D array of the element
 *                    at position (x,y,z).  Assumption is that
 *                    indices change fastest in Z-direction.
 *      Arguments:
 *          dim       3-element array
 *          x,y,z
 *
 *--------------------------------------------------------------------------*/
int EncodeFMCellIndex(int *dim, int x, int y, int z)
{
        return(z + (dim[Z]*y) + (dim[Z]*dim[Y]*x));
}


/*---------------------------------------------------------------------------
 *
 *      Function:     LookupFMCell
 *      Description:  Find a specified cell in the hash table and return
 *                    to the caller a pointer to the cell.  (Note: each
 *                    layer in the FM hierarchy has its own hash table
 *                    of cells.
 *      Arguments:
 *          cellTable  Pointer to the cell hash table in which to look
 *                     for the specified cell.
 *          cellID     ID number of the FM cell to locate.
 *
 *      Returns:  Pointer to the specified cell if found, NULL in all
 *                other cases.
 *
 *--------------------------------------------------------------------------*/
FMCell_t *LookupFMCell(FMCell_t *cellTable[], int cellID)
{
        int      hashVal;
        FMCell_t *cell = (FMCell_t *)NULL;

        hashVal = cellID % CELL_HASH_TABLE_SIZE;
        cell = cellTable[hashVal];

        while (cell != (FMCell_t *)NULL) {
            if (cell->cellID == cellID) {
                break;
            }
            cell = cell->next;
        }

        return(cell);
}



/*---------------------------------------------------------------------------
 *
 *      Function:    GaussQuadraturePointsAndWgts
 *      Description: Computes evaluation points xx and weights ww for n-point
 *                   Gaussian quadrature, so that
 *                         1             n-1
 *                        int f(x)dx  =  sum ww[i]*f(xx[i])
 *                        -1             i=0
 *
 *      Parameters:
 *          IN:   n        Number of points for gaussian quadrature
 *          OUT:  points   Array of <n> Gaussian quadrature points
 *          OUT:  weights  Array of <n> Gaussian quadrature weights
 *
 *  =========================================================================
 *  References:
 *
 *  http:
 *    www.mathworks.com/matlabcentral/fileexchange/loadCategory.do?objectId=16
 *
 *  Excerpt of relevant part of the above web-page:
 *    gaussq
 *    Numerically evaluates a integral using a Gauss quadrature.
 *    Author: Per A. Brodtkorb
 *    Category: Integration
 *
 *  GaussQuadraturePointsAndWgts is a C conversion of wfun=1 from
 *  qrule.m. Notes from qrule.m:
 *  %  wfun 1: copied from grule.m in NIT toolbox, see ref [2] 
 *  %  [2] Davis and Rabinowitz (1975) 'Methods of Numerical Integration',
 *  %  page 365, Academic Press.
 *
 *  Author note in qrule.m:
 *  % By Bryce Gardner, Purdue University, Spring 1993.
 *  % Modified by Per A. Brodtkorb 19.02.99 pab@marin.ntnu.no
 *  % to compute other quadratures  than the default
 *  =========================================================================
 */
void GaussQuadraturePointsAndWgts(int n, real8 *points, real8 *weights)
{
        int   ti, m, j, k;
        real8 e1, nn, t, xo, pkm1, pk, t1, pkp1;
        real8 den, d1, dpn, d2pn, d3pn, d4pn, u, v, h, p, dp, bp, fx, wf;

        m = (n+1)/2;
        e1 = n*(n+1);
        nn = 1.0 - (1.0 - 1.0/n)/(8*n*n);

        for (ti = 0; ti<m; ti++) {
            t = M_PI/(4*n+2)*(4*ti+3);
            xo = nn*cos(t);
            for (j = 0; j<2; j++) {
                pkm1 = 1.0;
                pk = xo;
                for (k = 2; k <= n; k++) {
                    t1 = xo*pk;
                    pkp1 = t1 - pkm1 - (t1-pkm1)/k + t1;
                    pkm1 = pk;
                    pk = pkp1;
                }
                den = 1.0 - xo*xo;
                d1 = n*(pkm1 - xo*pk);
                dpn = d1/den;
                d2pn = (2*xo*dpn - e1*pk)/den;
                d3pn = (4*xo*d2pn + (2.0-e1)*dpn)/den;
                d4pn = (6*xo*d3pn + (6.0-e1)*d2pn)/den;
                u = pk/dpn;
                v = d2pn/dpn;
                h = -u*(1.0 + (u/2)*(v + u*(v*v - u*d3pn/(3*dpn))));
                p = pk + h*(dpn + (h/2)*(d2pn + (h/3)*(d3pn + (h/4)*d4pn)));
                dp = dpn + h*(d2pn + (h/2)*(d3pn + h*d4pn/3));
                h = h - p/dp; xo = xo + h;
            }
            bp = -xo - h;
            fx = d1 - h*e1*(pk+(h/2)*(dpn+(h/3)*(d2pn + (h/4)*(d3pn + (h/5)*d4pn))));
            wf = 2*(1.0 - bp*bp)/(fx*fx);

            points[ti] = bp;
            points[n-1-ti] = -bp;
            weights[ti] = wf;
            weights[n-1-ti] = wf;
        }

        if (n%2 > 0) points[n/2] = 0.0;

        return;
}



#ifdef CALCENERGY

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

/*---------------------------------------------------------------------------
 *
 *      Function:    FMMEnergy
 *      Description: Energy contribution from two FMM cells
 *
 *      norder : linked to the size of etas
 *      R      : distance between the cell centers
 *      aeta1  : Multipole expansion array for first cell
 *      aeta2  : Multipole expansion array for second cell
 *
 *-------------------------------------------------------------------------*/

real8 FMMEnergy(Home_t *home, int norder, real8* R, 
                real8 *aeta1, real8* aeta2)
{
   real8 Energy = 0;

   real8 mu = home->param->shearModulus;
   real8 nu = home->param->pois;
   real8 rc  = home->param->rc;

   real8 mufac  = -mu/8.0/M_PI;
   real8 nufac2 = 2.0/(1-nu);
   real8 nufac1 = nu*nufac2;

   real8 r[3];
   r[0] = -R[0];
   r[1] = -R[1];
   r[2] = -R[2];

   // Factorials
   int i;
   real8 g[NMAX+1];   
   g[0] = 1.0;
   for (i = 1; i <= norder; i++) {
      g[i] = g[i-1] * (real8) i;
   }

   // Recover eta1 and eta2 in a form I can work with
   int k,n,idx,ii,jj;
   int nx,ny,nz;
   static real8 etaRa[NMAX+1][NMAX+1][NMAX+1][3][3];
   static real8 etaRb[NMAX+1][NMAX+1][NMAX+1][3][3];
   for (i = 0; i<=norder; i++) 
   {
      k = 0;
      n = (i+2) * (i+1) / 2;
      idx = 3 * n * i;
      for (nz = 0; nz <= i; nz++) 
         for (ny = 0; ny <= i-nz; ny++) 
         {
            nx = i-ny-nz;
            for (ii = 0; ii < 3; ii++) 
               for (jj = 0; jj < 3; jj++) 
               {      
                  // eta[3*ii+jj] = b[ii]t[jj] see makeeta
                  etaRa[nx][ny][nz][ii][jj] = aeta1[idx+n*(3*ii+jj)+k];
                  etaRb[nx][ny][nz][ii][jj] = aeta2[idx+n*(3*ii+jj)+k];
               }
            k = k+1;
         }
   }

   // Compute derivatives of R
   int maxorder = norder + 3;
   static int powvec[(MAXORDER+3+3)*(MAXORDER+3+2)*(MAXORDER+3+1)/6][3];
   double rderiv[(MAXORDER+3)*(MAXORDER+2)*(MAXORDER+1)/6];
   for (i = 0; i <= maxorder; i++) 
     DM2(i,r,rc,rderiv+(i+2)*(i+1)*i/6,powvec);

   // Calculate the energy from the two cells. 
   // According to LeSar and Rickman, formula (2), p. 144110-2.
   int kk,a,b,c;
   int i00,i11,i22,i01,i02,i12;
   double BinCs,d2R[3][3];
   real8 dkkR, valk;
   real8 sign=1;

   for (i = 0; i <= norder; i++) 
   {
      // First (nx, ny, nz) loop
      for (nz = 0; nz <= i; nz++) 
         for (ny = 0; ny <= i-nz; ny++) 
         {
            nx = i-ny-nz;

            // Second derivative of R
            i00 = ind(nx+2,ny  ,nz);
            i11 = ind(nx  ,ny+2,nz);
            i22 = ind(nx  ,ny  ,nz+2);
            
            i01 = ind(nx+1,ny+1,nz);
            i02 = ind(nx+1,ny  ,nz+1);
            i12 = ind(nx  ,ny+1,nz+1);           
            
            d2R[0][0] = rderiv[i00];
            d2R[1][1] = rderiv[i11];
            d2R[2][2] = rderiv[i22];
            
            d2R[0][1] = rderiv[i01]*0.5;
            d2R[0][2] = rderiv[i02]*0.5;
            d2R[1][2] = rderiv[i12]*0.5;
            
            d2R[1][0] = d2R[0][1];
            d2R[2][0] = d2R[0][2];
            d2R[2][1] = d2R[1][2];
            
            dkkR = d2R[0][0] + d2R[1][1] + d2R[2][2];

            // Second (a, b, c) loop
            for (c = 0; c <= nz; c++) 
               for (b = 0; b <= ny; b++)
                  for (a = 0; a <= nx; a++) 
                  {    
                     // Binomial coefficients
                     BinCs = 1.0/(g[a]*g[nx-a]) *
                             1.0/(g[b]*g[ny-b]) *
                             1.0/(g[c]*g[nz-c]);

                     // -ra
                     sign = ipow(-1,i-(a+b+c));

                     // Energy contribution
                     for (ii = 0; ii < 3; ii++)
                        for (jj = 0; jj < 3; jj++) 
                        {
                           Energy += sign*BinCs*dkkR*etaRb[a][b][c][jj][jj]*etaRa[nx-a][ny-b][nz-c][ii][ii] +
                              nufac1*sign*BinCs*dkkR*etaRb[a][b][c][ii][jj]*etaRa[nx-a][ny-b][nz-c][jj][ii];
                        }

                     for (ii = 0; ii < 3; ii++)
                        for (jj = 0; jj < 3; jj++) 
                        {
                           valk = 0.0;
                           for (kk = 0; kk < 3; kk++) 
                              valk += etaRb[a][b][c][ii][kk]*etaRa[nx-a][ny-b][nz-c][jj][kk];

                           Energy += sign*nufac2* BinCs*(d2R[ii][jj]-dkkR*(ii==jj))*valk;
                        }
                  }  
         } 
   }  

   return Energy*mufac;
}

#endif


