/**************************************************************************
 *
 *  Module:  G2sigma_tilde.cc
 *
 *  Includes functions:
 *      alpha_tilde2alpha()
 *      eta2eta_tilde()
 *      G2sigma_tilde_deriv()
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

#ifdef ANISOTROPIC
#include "G2sigma_tilde.h"

/*---------------------------------------------------------------------------
 *
 *  Function:     eta2eta_tilde
 *  Description:  Multiply the eta coefficient by C (the elastic
 *                constant matrix) to get eta_tilde
 *
 *--------------------------------------------------------------------------*/
void eta2eta_tilde(int norder,
                   real8 C[3][3][3][3],
                   real8 eta[],
                   real8 (*eta_tilde)[3][3])
{
    int       i, c, b;
    int       eta_idx0, eta_skip;
    real8     fact_buf[16], *ifact = fact_buf;
    const int nbuf = sizeof(fact_buf)/sizeof(*fact_buf);
    const int neta = (norder+3)*(norder+2)*(norder+1)/6;


    if (norder >= nbuf) {
        ifact = (real8 *) malloc(sizeof(real8) * (norder+1));
    }
  
/*
 *  Initialize factorial table
 */
    ifact[0] = 1.0;

    for (i = 1; i <= norder; i++) {      
        ifact[i] = ifact[i-1]*i;
        ifact[i-1] = 1.0/ifact[i-1];
    }

    ifact[norder] = 1.0/ifact[norder];

/*
 *  Formula:
 *     eta_tilde[q][a,b,c][g][p] = eps[n,g,r]*C[p,q,w,n]*eta[a,b,c][w][r];
 *
 *  NB! eta order is really
 *     eta[a,b,c][w][r] = eta[9*(i+2)*(i+1)*i/6 +
 *                        (3*w+r)*(i+1)*(i+2)/2 +
 *                        c*(2*n+3-c)/2 + b];
 */
    eta_idx0 = 0;
    eta_skip = 0;

    for (i = 0; i <= norder; i++) {
        int eta_idx1 = 0;

        eta_skip += i+1;

        for (c = 0; c <= i; c++) {
            for (b = 0; b <= i-c; b++) {
                int         p, q;
                const int   a = i-c-b;
                const real8 f = ifact[a]*ifact[b]*ifact[c];

                for (p = 0; p < 3; p++) {
                    for (q = 0; q < 3; q++) {
                        int       g;
                        const int eta_t_off = q*neta;
                        const int eta_t_idx = eta_t_off + eta_idx0 + eta_idx1;
    
                        for (g = 0; g < 3; g++) {
                            int       w;
                            const int cycshift = 1*1 + 2*4 + 0*16 + 1*64;
                            const int g2 = g+g;
                            const int r = (cycshift>>(g2))&3;
                            const int n = (cycshift>>(g2+2))&3;
                            real8     deta = 0.0;

                            for (w = 0; w < 3; w++) {
                                const int eta_wr_idx = 9*eta_idx0 + (3*w+r)*eta_skip + eta_idx1;
                                const int eta_wn_idx = 9*eta_idx0 + (3*w+n)*eta_skip + eta_idx1;

                                deta /*eta_tilde[eta_t_idx][g][p]*/ +=
                                    f*( C[p][q][w][n]*eta[eta_wr_idx] -
                                    C[p][q][w][r]*eta[eta_wn_idx] );

                            }

                            eta_tilde[eta_t_idx][g][p] = deta;

                        }  // end for (g = 0; g < 3; g++)

                    }  // end for (q = 0; q < 3; q++)

                }  // end for (p = 0; p < 3; p++)

                eta_idx1++;

            }  // end for (b = 0; b <= i-c; b++)

        }  // end for (c = 0; c <= i; c++)

        eta_idx0 += eta_idx1;

    }  // end for (i = 0; i <= norder; i++)
  
    if (norder >= nbuf) {
        free(ifact);
    }

    return;
}


/*---------------------------------------------------------------------------
 *
 *  Function:     alpha_tilde2alpha
 *  Description:  Multiply the alpha_tilde coefficients by C (the elastic
 *                constant matrix) returning alpha coefficients.
 *
 *--------------------------------------------------------------------------*/
void alpha_tilde2alpha(int uorder,
                       real8 C[3][3][3][3],
                       real8 alpha_tilde[][3][3],
                       real8 alpha[][3][3])
{
    int u, c, b, idx = 0;

/*
 *  Formula:
 *      alpha[a,b,c][j][s] = C[j,s,v,g]*alpha_tilde[a,b,c][v][g]
 */
    for (u = 0; u <= uorder; u++) {
        for (c = 0; c <= u; c++) {
            for (b = 0; b <= u-c; b++) {
                int j, s;

                for (j = 0; j < 3; j++) {
                    for (s = 0; s < 3; s++) {
                        int v, g;

                        alpha[idx][j][s] = 0.0;

                        for (v = 0; v < 3; v++) {
                            for (g = 0; g < 3; g++) {
                                alpha[idx][j][s] += C[j][s][v][g] *
                                                    alpha_tilde[idx][v][g];
                            }
                        }

                    }  // end for (s = 0; s < 3; s++)

                }  // end for (j = 0; j < 3; j++)

                idx++;

            }  // end for (b = 0; b <= u-c; b++)

        }  // end for (c = 0; c <= u; c++)

    }  // end for (u = 0; u <= uorder; u++)

    return;
}


/*---------------------------------------------------------------------------
 *
 *  Function:     G2sigma_tilde_deriv
 *  Description:  Compute the kernel, the derivative of the Green's function,
 *                for the FMM method.  This version does not contain any
 *                multiplication by the C (elastic constants) matrix.
 *                The multiplications are done outside this function.
 *
 *--------------------------------------------------------------------------*/
void G2sigma_tilde_deriv(int norder_array, int norder_eval,
                         int ndx, int ndy, int ndz,
                         real8 Gderiv[][3][3],
                         real8 eta_tilde[][3][3],
                         real8 sigma_tilde[3][3])
{
    real8 S00 = 0.0, S01 = 0.0, S02 = 0.0;
    real8 S10 = 0.0, S11 = 0.0, S12 = 0.0;
    real8 S20 = 0.0, S21 = 0.0, S22 = 0.0;

/*
 *  Build Greens function derivative from multipole coefficients
 */ 
    if (norder_eval >= 0) {
        const int neta = (norder_array+3)*(norder_array+2)*(norder_array+1)/6;
        int       q;
        int       qvec = 0*1 + 0*2 + 1*4;
        int       eta_off = 0;
        const int nq0 = 1 + ndx + ndy + ndz;
        const int Gq_off = (nq0+2)*(nq0+1)*nq0/6;

        for (q = 0; q < 3; q++) {
            int n;
            int eta_idx = eta_off;
            int nq = nq0;
            int Gq_idx0 = Gq_off;

            for (n = 0; n <= norder_eval; n++) {
                int c,b;
                int cq = (qvec&1) + ndz;
                int Gq_idx1 = ( ( cq*(2*nq+3-cq) )>>1 ) + ((qvec>>1)&1) + ndy;

                for (c = 0; c <= n; c++) {
                    int Gq_idx = Gq_idx0 + Gq_idx1;

                    for (b = 0; b <= n-c; b++) {

#define term(g,v,p) eta_tilde[eta_idx][g][p]*Gderiv[Gq_idx][v][p]
#define acc(g,v) S ## g ## v += term(g,v,0) + term(g,v,1) + term(g,v,2)

                        acc(0,0);  acc(0,1);  acc(0,2);
                        acc(1,0);  acc(1,1);  acc(1,2);
                        acc(2,0);  acc(2,1);  acc(2,2);

                        eta_idx++;
                        Gq_idx++;

                    } /* b-loop */

                    Gq_idx1 += nq+1 - cq;
                    cq++;

                } /* c-loop */

                nq++;
                Gq_idx0 += ((nq+1)*nq)>>1;

            } /* n-loop */

            qvec >>= 1;
            eta_off += neta;

        } /* q-loop */

    } /* norder-loop */

/*
 *  Return reduced stress derivative
 */
    sigma_tilde[0][0] = S00;
    sigma_tilde[0][1] = S10;
    sigma_tilde[0][2] = S20;

    sigma_tilde[1][0] = S01;
    sigma_tilde[1][1] = S11;
    sigma_tilde[1][2] = S21;

    sigma_tilde[2][0] = S02;
    sigma_tilde[2][1] = S12;
    sigma_tilde[2][2] = S22;

    return;
}
#endif  // ifdef ANISOTROPIC
