/***************************************************************************
 *
 *      Author:       Chao Chen, Eric Darve and Sylvie Aubry
 *      Module:       CorrectionTable.c
 *      Description:  Contains functions for creating a table used
 *                    in conjuction with the Fast Multipole functions
 *                    to adjust the calculated stress to allow for
 *                    periodic images of the problem space.
 *
 *      Includes public functions:
 *
 *          bbFMM_CorrectionTableInit()
 *          bbFMM_CreateCorrectionTable()
 *          bbFMM_DoTableCorrection()
 *          bbFMM_FreeCorrectionTable()
 *
 *      Includes private functions:
 *
 *          ComputeWeightsPBC
 *
***************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "Home.h"
#include "FM.h"

#ifdef PARALLEL
#include  <mpi.h>
#endif

#ifdef UNIFORMFMM
#include "FM.h"

static real8 *correctionTbl = NULL;

// Unif2ChebGrid()
//
// Transform uniform grid points to Chebychev grid points
//---------------------------------------------------------------------------------

static void Unif2ChebGrid(int dof, real8 *Sn[3], real8 *C)
{
    int    n  = uniformfmm->n;
    int    n3 = uniformfmm->n3;

    real8 *U  = (real8 *)calloc(1,n3*dof * sizeof(real8));

    for (int i=0; i<(dof*n3); i++)
    {
        U[i] = C[i];
        C[i] = 0.0;
    }

    int lhs_idx=0;
    for (int l1=0; l1<n; l1++)
    {
        int x_idx = l1*n3;
        for (int l2=0; l2<n; l2++)
        {
            int y_idx = l2*n3;
            for (int l3=0; l3<n; l3++)
            {
                int z_idx = l3*n3;
                for (int j=0; j<n3; j++)
                {
                    real8 weight  = Sn[0][x_idx+j] *Sn[1][y_idx+j] *Sn[2][z_idx+j];
                    int   rhs_idx = j * dof;

                    for (int l4=0; l4<dof; l4++)
                        C[lhs_idx+l4] += weight * U[rhs_idx++];

                }
                lhs_idx += dof;
            }
        }
    }

    free(U);
}

// Cheb2UnifGrid()
//
// Transform Chebychev grid points to uniform grid points
//---------------------------------------------------------------------------------

static void Cheb2UnifGrid(int dof, real8 *Sn[3], real8 *U)
{
    int    n  = uniformfmm->n;
    int    n3 = uniformfmm->n3;

    real8 *C  = (real8 *)calloc(1,n3*dof * sizeof(real8));

    for (int i=0; i<dof*n3; i++)
    {
        C[i] = U[i];
        U[i] = 0;
    }

    for (int i=0; i<n3; i++)
    {
        int lhs_idx = dof * i;
        int rhs_idx = 0;
        for (int l1=0; l1<n; l1++)
        {
            real8 w1 = Sn[0][l1*n3+i];
            for (int l2=0; l2<n; l2++)
            {
                real8 w2 = w1 * Sn[1][l2*n3+i];
                for (int l3=0; l3<n; l3++)
                {
                    real8 w3 = w2 * Sn[2][l3*n3+i];

                    for (int idof=0; idof<dof; idof++)
                        U[lhs_idx+idof] += w3 * C[rhs_idx++];
                }
            }
        }
    }

    free(C);
}

// uniformFMM_DoTableCorrection()
//
// For simplicity, use Chebychev PBC correction in the
// Uniform grid method. Need to do some pre and post
// processing in order to do that.
//---------------------------------------------------------------------------------

void uniformFMM_DoTableCorrection(Home_t *home)
{
    real8 *Grid, *Grid3D[3], *Sn[3];

    int   n       = uniformfmm->n;
    int   n3      = uniformfmm->n3;
    real8 homogen = uniformfmm->homogen;

    FMLayer_t *layer = &home->fmLayer[0];
    if (layer->ownedCnt == 0) return;

    FMCell_t  *cell = LookupFMCell(layer->cellTable, 0);

    // Create a uniform grid
    Grid      = (real8 *) calloc(1,n  * sizeof(real8));
    Grid3D[0] = (real8 *) calloc(1,n3 * sizeof(real8));
    Grid3D[1] = (real8 *) calloc(1,n3 * sizeof(real8));
    Grid3D[2] = (real8 *) calloc(1,n3 * sizeof(real8));

    if (Grid     ==NULL) Fatal("memory allocation failure\n");
    if (Grid3D[0]==NULL) Fatal("memory allocation failure\n");
    if (Grid3D[1]==NULL) Fatal("memory allocation failure\n");
    if (Grid3D[2]==NULL) Fatal("memory allocation failure\n");

    real8 deltan = 2.0/((real8)(n-1));
    for (int m=0; m<n; m++)
    {
        Grid[m] = (1.0 - deltan*m);
    }

    int count=0;
    for (int l1=0; l1<n; l1++)
    {
        for (int l2=0; l2<n; l2++)
        {
            for (int l3=0; l3<n; l3++)
            {
                Grid3D[0][count] = Grid[l1];
                Grid3D[1][count] = Grid[l2];
                Grid3D[2][count] = Grid[l3];
                count++;
            }
        }
    }

    free(Grid); Grid = NULL;

    // Compute Sn
    Sn[0] = (real8 *)calloc(1,n3*n * sizeof(real8));
    Sn[1] = (real8 *)calloc(1,n3*n * sizeof(real8));
    Sn[2] = (real8 *)calloc(1,n3*n * sizeof(real8));

    if (Sn[0]==NULL) Fatal("memory allocation failure\n");
    if (Sn[1]==NULL) Fatal("memory allocation failure\n");
    if (Sn[2]==NULL) Fatal("memory allocation failure\n");

    ComputeSn2(Grid3D, uniformfmm->Tkz, n, n3, Sn);

    free(Grid3D[0]); Grid3D[0]=0;
    free(Grid3D[1]); Grid3D[1]=0;
    free(Grid3D[2]); Grid3D[2]=0;

    // Uniform to Chebychev
    Unif2ChebGrid(9, Sn, cell->mpCoeff);

    Stanford_DoTableCorrection(home,n,homogen);

#ifndef PARALLEL
    // Determine the mean stress here
    if (    (home->param->xBoundType == Periodic) 
         && (home->param->yBoundType == Periodic) 
         && (home->param->zBoundType == Periodic) )
    {
        MeanStressCorrection(home); //IN PARALLEL THIS CALL HANGS....
    }
#endif

    // Chebychev to Uniform
    Cheb2UnifGrid(6, Sn, cell->expansionCoeff);

    free(Sn[0]); Sn[0]=0;
    free(Sn[1]); Sn[1]=0;
    free(Sn[2]); Sn[2]=0;
}

#endif
