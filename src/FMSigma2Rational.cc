/**************************************************************************
 *
 *      Module:       FMSigma2Rational.c
 *      Description:  This module contains the base FMSigma2Rational
 *                    entry function used by the fast multipole code,
 *                    a generic function for computing the stress
 *                    values for any arbitrary expansion order
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
 *                    obtain the divergence-free stress field.  This
 *                    version of FMSigma2 has been modified
 *                    from the original version to include the same
 *                    stress symmetrization as the segment stress
 *                    routine for handling 'rational' segment.
 *
 *      Includes public functions:
 *          FMSigma2Rational()
 *
 *      Includes private functions:
 *          FMSigma2core()
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"

/*
 *      Prototype the basic unoptimized function to handle
 *      any arbitrary expansion order, as well as any
 *      expansion-order specific optimized functions.
 */
static void FMSigma2core(int iorder, int pows[][3], real8 terms[],
                         real8 mu8pi, real8 two1nu, real8 Eeta[], matrix sigma);

/*---------------------------------------------------------------------------
 *
 *      Function:    FMSigma2Rational
 *      Description: This routine is a little more than a dispatching
 *                   function which does some initialization and invokes
 *                   FMSigma2Core to do the real work.
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
void FMSigma2Rational(real8 mu, real8 nu, int norder, real8 Eeta[],
                  matrix sigmatot, real8 pcterms[], int pcpows[][3])
{
        int    iorder, i, j, (*pows)[3], npows;
        real8  *terms;
        real8  mu8pi = mu/(8.0*M_PI);
        real8  two1nu = 2.0/(1.0-nu);  
        matrix sigma;

        static real8 fact[NMAX+1],ifact[NMAX+1],dfact[2*NMAX-3+2];
        static int inited = 0;

/*
 *      Make sure only a single thread initializes the *fact tables
 */
        if (inited == 0) {
#ifdef _OPENMP
#pragma omp critical (CRIT_INIT_FMSIGMA2)
#endif
            {
                if (inited == 0) {
                    makeftabs(fact, ifact, dfact);
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

        for (iorder = 0; iorder <= norder; iorder++) {
            real8 t;

            npows = (iorder+1+3) * (iorder+2+3) / 2;

            FMSigma2core(iorder, pows, terms, mu8pi, two1nu, Eeta, sigma);

            t = ipow(-1.0, iorder);

            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    sigma[i][j] = t * sigma[i][j] * ifact[iorder];
                }
            }
        
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
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


static void FMSigma2core(int iorder, int pows[][3], real8 terms[], real8 mu8pi,
                         real8 two1nu, real8 Eeta[], matrix sigma) {
        int k, a, nx, ny, nz, npows, nidx;
        int etaoff, etadis0;
        real8 sig00, sig01, sig02, sig11, sig12, sig22;
        real8 t1, t2, t3, t4, t5, t6, t7, t8, t9;
        real8 u1, u2, u3, u4, u5, u6, u7, u8;
        real8 *Eeta1, *Eeta2, *Eeta3, *Eeta4;
        real8 *Eeta5, *Eeta6, *Eeta7, *Eeta8;
        real8 *Eeta9;
        real8 tvec[NTMAX];

        static real8 fact[NMAX+1],ifact[NMAX+1],dfact[2*NMAX-3+2];
        static int inited = 0;

        if (inited == 0) {
#ifdef _OPENMP
#pragma omp critical (CRIT_INIT_FMSIGMA2CORE)
#endif
            {
                if (inited == 0) {
                    makeftabs(fact,ifact,dfact);
                    inited = 1;
                }
            }
        }

        sig00 = 0.0;
        sig01 = 0.0;
        sig02 = 0.0;
        sig11 = 0.0;
        sig12 = 0.0;
        sig22 = 0.0;

        npows = (iorder+1+3)*(iorder+2+3)/2;

        etadis0 = (iorder+1)*(iorder+2) >> 1;  /* (iorder+1)*(iorder+2)/2 */
        etaoff = 3*iorder*etadis0;  /* 9*iorder*(iorder+1)*(iorder+2)/6 */

/*
 *      Term: i = 0, j = 0, m = 0
 */
        nidx = 0;

        for (k = 0; k < npows; k++) {

            nx = pows[k][0] - 3;
            ny = pows[k][1] - 0;
            nz = pows[k][2] - 0;

            if (nx >= 0 && ny >= 0 && nz >= 0) {
                tvec[nidx++] = terms[k] * fact[iorder] *
                               ifact[nx]*ifact[ny]*ifact[nz];
            }
        }

        etadis0 = etaoff + nidx * 4;

        t1 = 0.0; t2 = 0.0; t3 = 0.0; t4 = 0.0;

        Eeta1 = &Eeta[etadis0];
        Eeta2 = &Eeta[etadis0+nidx];
        Eeta3 = &Eeta[etadis0+nidx*3];
        Eeta4 = &Eeta[etadis0+nidx*4];

        for (a = 0; a < nidx; a++) {
            real8 tvecVal = tvec[a];

            t1 += tvecVal * *(Eeta1++);
            t2 += tvecVal * *(Eeta2++);
            t3 += tvecVal * *(Eeta3++);
            t4 += tvecVal * *(Eeta4++);
        }

        u2 = t2 * two1nu;
        u3 = t3 * two1nu;

        sig11 = sig11  - u2;
        sig11 = sig11  - 2*t3 + u3;

        sig12 = sig12  + t1;
        sig12 = sig12  - t4;

        sig22 = sig22  + 2*t2 - u2;
        sig22 = sig22  + u3;

/*
 *      Term: i = 0, j = 0, m = 1
 */
        nidx = 0;

        for (k = 0; k < npows; k++) {

            nx = pows[k][0] - 2;
            ny = pows[k][1] - 1;
            nz = pows[k][2] - 0;

            if (nx >= 0 && ny >= 0 && nz >= 0) {
                tvec[nidx++] = terms[k] * fact[iorder] *
                               ifact[nx]*ifact[ny]*ifact[nz];
            }
        }

        etadis0 = etaoff;
        etadis0 += nidx;

/*
 *      Term: i = 0, j = 0, m = 1, k = 0, l = 1
 */
        t1 = 0.0; t2 = 0.0; t3 = 0.0; t4 = 0.0;
        t5 = 0.0; t6 = 0.0; t7 = 0.0; t8 = 0.0;

        Eeta1 = &Eeta[etadis0];
        Eeta2 = &Eeta[etadis0+nidx];
        Eeta3 = &Eeta[etadis0+nidx*2];
        Eeta4 = &Eeta[etadis0+nidx*3];
        Eeta5 = &Eeta[etadis0+nidx*4];
        Eeta6 = &Eeta[etadis0+nidx*5];
        Eeta7 = &Eeta[etadis0+nidx*6];
        Eeta8 = &Eeta[etadis0+nidx*7];

        for (a = 0; a < nidx; a++) {
            real8 tvecVal = tvec[a];

            t1 += tvecVal * *(Eeta1++);
            t2 += tvecVal * *(Eeta2++);
            t3 += tvecVal * *(Eeta3++);
            t4 += tvecVal * *(Eeta4++);
            t5 += tvecVal * *(Eeta5++);
            t6 += tvecVal * *(Eeta6++);
            t7 += tvecVal * *(Eeta7++);
            t8 += tvecVal * *(Eeta8++);
        }

        u2 = t2 * two1nu;
        u5 = t5 * two1nu;
        u6 = t6 * two1nu;
        u7 = t7 * two1nu;

        sig01 = sig01  + u5;
        sig01 = sig01  + 2*t7 - u7;

        sig02 = sig02  - t4;
        sig02 = sig02  + t8;

        sig11 = sig11  + u2;
        sig11 = sig11  + 2*t6 - u6;

        sig12 = sig12  - t1;
        sig12 = sig12  - t3;

        sig22 = sig22  - 2*t2 + u2;
        sig22 = sig22  - u6;


/*
 *      Term: i = 0, j = 0, m = 2
 */
        nidx = 0;

        for (k = 0; k < npows; k++) {

            nx = pows[k][0] - 2;
            ny = pows[k][1] - 0;
            nz = pows[k][2] - 1;

            if (nx >= 0 && ny >= 0 && nz >= 0) {
                tvec[nidx++] = terms[k] * fact[iorder] *
                               ifact[nx]*ifact[ny]*ifact[nz];
            }
        }

        etadis0 = etaoff + nidx;

/*
 *      Term: i = 0, j = 0, m = 2, k = 0, l = 1
 */
        t1 = 0.0; t2 = 0.0; t3 = 0.0; t4 = 0.0;
        t5 = 0.0; t6 = 0.0; t7 = 0.0; t8 = 0.0;

        Eeta1 = &Eeta[etadis0];
        Eeta2 = &Eeta[etadis0+nidx];
        Eeta3 = &Eeta[etadis0+nidx*2];
        Eeta4 = &Eeta[etadis0+nidx*3];
        Eeta5 = &Eeta[etadis0+nidx*4];
        Eeta6 = &Eeta[etadis0+nidx*5];
        Eeta7 = &Eeta[etadis0+nidx*6];
        Eeta8 = &Eeta[etadis0+nidx*7];

        for (a = 0; a < nidx; a++) {
            real8 tvecVal = tvec[a];

            t1 += tvecVal * *(Eeta1++);
            t2 += tvecVal * *(Eeta2++);
            t3 += tvecVal * *(Eeta3++);
            t4 += tvecVal * *(Eeta4++);
            t5 += tvecVal * *(Eeta5++);
            t6 += tvecVal * *(Eeta6++);
            t7 += tvecVal * *(Eeta7++);
            t8 += tvecVal * *(Eeta8++);
        }

        u1 = t1 * two1nu;
        u3 = t3 * two1nu;
        u5 = t5 * two1nu;
        u7 = t7 * two1nu;

        sig01 = sig01  - t4;
        sig01 = sig01  + t8;

        sig02 = sig02  - 2*t5 + u5;
        sig02 = sig02  - u7;

        sig11 = sig11  + 2*t1 - u1;
        sig11 = sig11  + u3;

        sig12 = sig12  + t2;
        sig12 = sig12  + t6;

        sig22 = sig22  - u1;
        sig22 = sig22  - 2*t3 + u3;

/*
 *      Term: i = 0, j = 1, m = 1
 */
        nidx = 0;

        for (k = 0; k < npows; k++) {

            nx = pows[k][0] - 1;
            ny = pows[k][1] - 2;
            nz = pows[k][2] - 0;

            if (nx >= 0 && ny >= 0 && nz >= 0) {
                tvec[nidx++] = terms[k] * fact[iorder] *
                               ifact[nx]*ifact[ny]*ifact[nz];
            }
        }

        etadis0 = etaoff;

        t1 = 0.0; t2 = 0.0; t3 = 0.0; t4 = 0.0;
        t5 = 0.0; t6 = 0.0; t7 = 0.0; t8 = 0.0;

        Eeta1 = &Eeta[etadis0];
        Eeta2 = &Eeta[etadis0+nidx];
        Eeta3 = &Eeta[etadis0+nidx*2];
        Eeta4 = &Eeta[etadis0+nidx*3];
        Eeta5 = &Eeta[etadis0+nidx*5];
        Eeta6 = &Eeta[etadis0+nidx*6];
        Eeta7 = &Eeta[etadis0+nidx*7];
        Eeta8 = &Eeta[etadis0+nidx*8];

        for (a = 0; a < nidx; a++) {
            real8 tvecVal = tvec[a];

            t1 += tvecVal * *(Eeta1++);
            t2 += tvecVal * *(Eeta2++);
            t3 += tvecVal * *(Eeta3++);
            t4 += tvecVal * *(Eeta4++);
            t5 += tvecVal * *(Eeta5++);
            t6 += tvecVal * *(Eeta6++);
            t7 += tvecVal * *(Eeta7++);
            t8 += tvecVal * *(Eeta8++);
        }

        u3 = t3 * two1nu;
        u5 = t5 * two1nu;
        u6 = t6 * two1nu;
        u7 = t7 * two1nu;

        sig00 = sig00  - u5;
        sig00 = sig00  - 2*t7 + u7;

        sig01 = sig01  - u3;
        sig01 = sig01  - 2*t6 + u6;

        sig02 = sig02  + t2;
        sig02 = sig02  + t4;

        sig12 = sig12  + t1;
        sig12 = sig12  - t8;

        sig22 = sig22  + 2*t5 - u5;
        sig22 = sig22  + u7;

/*
 *      Term: i = 0, j = 1, m = 2
 */
        nidx = 0;

        for (k = 0; k < npows; k++) {

            nx = pows[k][0] - 1;
            ny = pows[k][1] - 1;
            nz = pows[k][2] - 1;

            if (nx >= 0 && ny >= 0 && nz >= 0) {
                tvec[nidx++] = terms[k] * fact[iorder] *
                               ifact[nx]*ifact[ny]*ifact[nz];
            }
        }

        etadis0 = etaoff;

        t1 = 0.0; t2 = 0.0; t3 = 0.0; t4 = 0.0;
        t5 = 0.0; t6 = 0.0; t7 = 0.0; t8 = 0.0;
        t9 = 0.0;

        Eeta1 = &Eeta[etadis0];
        Eeta2 = &Eeta[etadis0+nidx];
        Eeta3 = &Eeta[etadis0+nidx*2];
        Eeta4 = &Eeta[etadis0+nidx*3];
        Eeta5 = &Eeta[etadis0+nidx*4];
        Eeta6 = &Eeta[etadis0+nidx*5];
        Eeta7 = &Eeta[etadis0+nidx*6];
        Eeta8 = &Eeta[etadis0+nidx*7];
        Eeta9 = &Eeta[etadis0+nidx*8];

        for (a = 0; a < nidx; a++) {
            real8 tvecVal = tvec[a];

            t1 += tvecVal * *(Eeta1++);
            t2 += tvecVal * *(Eeta2++);
            t3 += tvecVal * *(Eeta3++);
            t4 += tvecVal * *(Eeta4++);
            t5 += tvecVal * *(Eeta5++);
            t6 += tvecVal * *(Eeta6++);
            t7 += tvecVal * *(Eeta7++);
            t8 += tvecVal * *(Eeta8++);
            t9 += tvecVal * *(Eeta9++);
        }

        u2 = t2 * two1nu;
        u3 = t3 * two1nu;
        u4 = t4 * two1nu;
        u6 = t6 * two1nu;
        u7 = t7 * two1nu;
        u8 = t8 * two1nu;

        sig00 = sig00  + 2*t5;
        sig00 = sig00  - 2*t9;

        sig01 = sig01  - t2 + u2;
        sig01 = sig01  + t4 - u4;

        sig02 = sig02  + t3 - u3;
        sig02 = sig02  - t7 + u7;

        sig11 = sig11  - 2*t1;
        sig11 = sig11  + 2*t9;

        sig12 = sig12  - t6 + u6;
        sig12 = sig12  + t8 - u8;

        sig22 = sig22  + 2*t1;
        sig22 = sig22  - 2*t5;

/*
 *      Term: i = 0, j = 2, m = 2
 */
        nidx = 0;

        for (k = 0; k < npows; k++) {

            nx = pows[k][0] - 1;
            ny = pows[k][1] - 0;
            nz = pows[k][2] - 2;

            if (nx >= 0 && ny >= 0 && nz >= 0) {
                tvec[nidx++] = terms[k] * fact[iorder] *
                               ifact[nx]*ifact[ny]*ifact[nz];
            }
        }

        etadis0 = etaoff;

        t1 = 0.0; t2 = 0.0; t3 = 0.0; t4 = 0.0;
        t5 = 0.0; t6 = 0.0; t7 = 0.0; t8 = 0.0;

        Eeta1 = &Eeta[etadis0];
        Eeta2 = &Eeta[etadis0+nidx];
        Eeta3 = &Eeta[etadis0+nidx*2];
        Eeta4 = &Eeta[etadis0+nidx*3];
        Eeta5 = &Eeta[etadis0+nidx*4];
        Eeta6 = &Eeta[etadis0+nidx*5];
        Eeta7 = &Eeta[etadis0+nidx*6];
        Eeta8 = &Eeta[etadis0+nidx*7];

        for (a = 0; a < nidx; a++) {
            real8 tvecVal = tvec[a];

            t1 += tvecVal * *(Eeta1++);
            t2 += tvecVal * *(Eeta2++);
            t3 += tvecVal * *(Eeta3++);
            t4 += tvecVal * *(Eeta4++);
            t5 += tvecVal * *(Eeta5++);
            t6 += tvecVal * *(Eeta6++);
            t7 += tvecVal * *(Eeta7++);
            t8 += tvecVal * *(Eeta8++);
        }

        u2 = t2 * two1nu;
        u4 = t4 * two1nu;
        u6 = t6 * two1nu;
        u8 = t8 * two1nu;

        sig00 = sig00  + 2*t6 - u6;
        sig00 = sig00  + u8;

        sig01 = sig01  - t7;
        sig01 = sig01  - t3;

        sig02 = sig02  + u2;
        sig02 = sig02  + 2*t4 - u4;

        sig11 = sig11  - u6;
        sig11 = sig11  - 2*t8 + u8;

        sig12 = sig12  - t1;
        sig12 = sig12  + t5;

/*
 *      Term: i = 1, j = 1, m = 1
 */
        nidx = 0;

        for (k = 0; k < npows; k++) {

            nx = pows[k][0] - 0;
            ny = pows[k][1] - 3;
            nz = pows[k][2] - 0;

            if (nx >= 0 && ny >= 0 && nz >= 0) {
                tvec[nidx++] = terms[k] * fact[iorder] *
                               ifact[nx]*ifact[ny]*ifact[nz];
            }
        }

        etadis0 = etaoff;

        t1 = 0.0; t2 = 0.0; t3 = 0.0; t4 = 0.0;

        Eeta1 = &Eeta[etadis0];
        Eeta2 = &Eeta[etadis0+nidx*2];
        Eeta3 = &Eeta[etadis0+nidx*6];
        Eeta4 = &Eeta[etadis0+nidx*8];

        for (a = 0; a < nidx; a++) {
            real8 tvecVal = tvec[a];

            t1 += tvecVal * *(Eeta1++);
            t2 += tvecVal * *(Eeta2++);
            t3 += tvecVal * *(Eeta3++);
            t4 += tvecVal * *(Eeta4++);
        }

        u2 = t2 * two1nu;
        u3 = t3 * two1nu;

        sig00 = sig00  + u2;
        sig00 = sig00  + 2*t3 - u3;

        sig02 = sig02  - t1;
        sig02 = sig02  + t4;

        sig22 = sig22  - 2*t2 + u2;
        sig22 = sig22  - u3;


/*
 *      Term: i = 1, j = 1, m = 2
 */
        nidx = 0;

        for (k = 0; k < npows; k++) {

            nx = pows[k][0] - 0;
            ny = pows[k][1] - 2;
            nz = pows[k][2] - 1;

            if (nx >= 0 && ny >= 0 && nz >= 0) {
                tvec[nidx++] = terms[k] * fact[iorder] *
                               ifact[nx]*ifact[ny]*ifact[nz];
            }
        }

        etadis0 = etaoff;

        t1 = 0.0; t2 = 0.0; t3 = 0.0; t4 = 0.0;
        t5 = 0.0; t6 = 0.0; t7 = 0.0; t8 = 0.0;

        Eeta1 = &Eeta[etadis0];
        Eeta2 = &Eeta[etadis0+nidx];
        Eeta3 = &Eeta[etadis0+nidx*2];
        Eeta4 = &Eeta[etadis0+nidx*3];
        Eeta5 = &Eeta[etadis0+nidx*5];
        Eeta6 = &Eeta[etadis0+nidx*6];
        Eeta7 = &Eeta[etadis0+nidx*7];
        Eeta8 = &Eeta[etadis0+nidx*8];

        for (a = 0; a < nidx; a++) {
            real8 tvecVal = tvec[a];

            t1 += tvecVal * *(Eeta1++);
            t2 += tvecVal * *(Eeta2++);
            t3 += tvecVal * *(Eeta3++);
            t4 += tvecVal * *(Eeta4++);
            t5 += tvecVal * *(Eeta5++);
            t6 += tvecVal * *(Eeta6++);
            t7 += tvecVal * *(Eeta7++);
            t8 += tvecVal * *(Eeta8++);
        }

        u2 = t2 * two1nu;
        u3 = t3 * two1nu;
        u4 = t4 * two1nu;
        u6 = t6 * two1nu;

        sig00 = sig00  - u2;
        sig00 = sig00  - 2*t4 + u4;

        sig01 = sig01  + t1;
        sig01 = sig01  - t8;

        sig02 = sig02  - t5;
        sig02 = sig02  - t7;

        sig12 = sig12  + 2*t3 - u3;
        sig12 = sig12  + u6;

        sig22 = sig22  + 2*t2 - u2;
        sig22 = sig22  + u4;

/*
 *      Term: i = 1, j = 2, m = 2
 */
        nidx = 0;

        for (k = 0; k < npows; k++) {

            nx = pows[k][0] - 0;
            ny = pows[k][1] - 1;
            nz = pows[k][2] - 2;

            if (nx >= 0 && ny >= 0 && nz >= 0) {
                tvec[nidx++] = terms[k] * fact[iorder] *
                               ifact[nx]*ifact[ny]*ifact[nz];
            }
        }

        etadis0 = etaoff;

        t1 = 0.0; t2 = 0.0; t3 = 0.0; t4 = 0.0;
        t5 = 0.0; t6 = 0.0; t7 = 0.0; t8 = 0.0;

        Eeta1 = &Eeta[etadis0];
        Eeta2 = &Eeta[etadis0+nidx];
        Eeta3 = &Eeta[etadis0+nidx*2];
        Eeta4 = &Eeta[etadis0+nidx*3];
        Eeta5 = &Eeta[etadis0+nidx*4];
        Eeta6 = &Eeta[etadis0+nidx*5];
        Eeta7 = &Eeta[etadis0+nidx*6];
        Eeta8 = &Eeta[etadis0+nidx*7];

        for (a = 0; a < nidx; a++) {
            real8 tvecVal = tvec[a];

            t1 += tvecVal * *(Eeta1++);
            t2 += tvecVal * *(Eeta2++);
            t3 += tvecVal * *(Eeta3++);
            t4 += tvecVal * *(Eeta4++);
            t5 += tvecVal * *(Eeta5++);
            t6 += tvecVal * *(Eeta6++);
            t7 += tvecVal * *(Eeta7++);
            t8 += tvecVal * *(Eeta8++);
        }

        u2 = t2 * two1nu;
        u3 = t3 * two1nu;
        u4 = t4 * two1nu;
        u7 = t7 * two1nu;

        sig00 = sig00  + u3;
        sig00 = sig00  + 2*t7 - u7;

        sig01 = sig01  + t6;
        sig01 = sig01  + t8;

        sig02 = sig02  - t1;
        sig02 = sig02  + t5;

        sig11 = sig11  - 2*t3 + u3;
        sig11 = sig11  - u7;

        sig12 = sig12  - 2*t2 + u2;
        sig12 = sig12  - u4;

/*
 *      Term: i = 2, j = 2, m = 2
 */
        nidx = 0;

        for (k = 0; k < npows; k++) {

            nx = pows[k][0] - 0;
            ny = pows[k][1] - 0;
            nz = pows[k][2] - 3;

            if (nx >= 0 && ny >= 0 && nz >= 0) {
                tvec[nidx++] = terms[k] * fact[iorder] *
                               ifact[nx]*ifact[ny]*ifact[nz];
            }
        }

        etadis0 = etaoff;

        t1 = 0.0; t2 = 0.0; t3 = 0.0; t4 = 0.0;

        Eeta1 = &Eeta[etadis0];
        Eeta2 = &Eeta[etadis0+nidx];
        Eeta3 = &Eeta[etadis0+nidx*3];
        Eeta4 = &Eeta[etadis0+nidx*4];

        for (a = 0; a < nidx; a++) {
            real8 tvecVal = tvec[a];

            t1 += tvecVal * *(Eeta1++);
            t2 += tvecVal * *(Eeta2++);
            t3 += tvecVal * *(Eeta3++);
            t4 += tvecVal * *(Eeta4++);
        }

        u2 = t2 * two1nu;
        u3 = t3 * two1nu;

        sig00 = sig00  - u2;
        sig00 = sig00  - 2*t3 + u3;

        sig01 = sig01  + t1;
        sig01 = sig01  - t4;

        sig11 = sig11  + 2*t2 - u2;
        sig11 = sig11  + u3;

        sigma[0][0] = mu8pi * sig00;
        sig01 = mu8pi * sig01;
        sigma[0][1] = sig01;
        sig02 = mu8pi * sig02;
        sigma[0][2] = sig02;
        sigma[1][0] = sig01;
        sigma[1][1] = mu8pi * sig11;
        sig12 = mu8pi * sig12;
        sigma[1][2] = sig12;
        sigma[2][0] = sig02;
        sigma[2][1] = sig12;
        sigma[2][2] = mu8pi * sig22;

        return;
}
