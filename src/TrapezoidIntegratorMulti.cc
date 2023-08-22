/**************************************************************************
 *
 *      Module:      TrapezoidIntegratorMulti.c
 *      Description: Implements a numerical timestep integrator using
 *                   the Trapezoid integration method. The integrator
 *                   is called multiple times until the maxDT time step
 *                   is reached.
 *
 ***************************************************************************/
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"
#include "Node.h"
#include "Timestep.h"

/*---------------------------------------------------------------------------
 *
 *      Function:     TrapezoidIntegratorMulti
 *
 *-------------------------------------------------------------------------*/

void TrapezoidIntegratorMulti(Home_t *home)
{
    // count current native nodes (skipping over null entries)...

    int ncnt=0;

    if (home->nodeKeys)
    {
        for (int i=0; i<home->newNodeKeyPtr; i++)
            if (home->nodeKeys[i]) ncnt++;
    }

    // add ghost nodes...

    for (Node_t *node=home->ghostNodeQ; (node); node=node->next ) { ncnt++; }

    // allocate space and store the current node positions and velocities.

    real8 *oldP = ( (ncnt>0) ? (real8 *) malloc(ncnt*3*sizeof(real8)) : 0 );
    real8 *oldV = ( (ncnt>0) ? (real8 *) malloc(ncnt*3*sizeof(real8)) : 0 );

    if (oldP && oldV)
    {
        int j=0;
        for (int i=0; i<home->newNodeKeyPtr; i++) 
        {
            Node_t *node = home->nodeKeys[i];

            if (node)
            {
               oldV[j  ] = node->vX;   oldP[j  ] = node->x;
               oldV[j+1] = node->vY;   oldP[j+1] = node->y;
               oldV[j+2] = node->vZ;   oldP[j+2] = node->z;
               j+=3;
            }
        }

        for (Node_t *node=home->ghostNodeQ; (node); node=node->next ) 
        {
            oldV[j  ] = node->vX;   oldP[j  ] = node->x;
            oldV[j+1] = node->vY;   oldP[j+1] = node->y;
            oldV[j+2] = node->vZ;   oldP[j+2] = node->z;
            j+=3;
        }
    }

    // Time integrate using the trapezoid integrator.

    Param_t *param = home->param;

    param->deltaTT = param->nextDT;

    int   totalSteps = 0;
    real8 totalDT    = 0.0;

    while (totalDT < param->maxDT && totalSteps < 20) 
    {
        TrapezoidIntegrator(home);
        totalSteps++;
        totalDT += param->realdt;
    }

    // if (home->myDomain == 0)
    //     printf("Cycle = %d: totalSteps = %d\n",home->cycle,totalSteps);

    param->deltaTT   = totalDT;
    param->realdt    = totalDT;
    param->timeStart = param->timeNow;

    // Now restore the old positions and velocities

    if (oldP && oldV)
    {
        int j=0;
        for (int i=0; i<home->newNodeKeyPtr; i++) 
        {
            Node_t *node = home->nodeKeys[i];

            if (node)
            {
                node->oldvX = oldV[j  ];   node->oldx = oldP[j  ];
                node->oldvY = oldV[j+1];   node->oldy = oldP[j+1];
                node->oldvZ = oldV[j+2];   node->oldz = oldP[j+2];
                j+=3;
            }
        }

        for (Node_t *node=home->ghostNodeQ; (node); node=node->next ) 
        {
            node->oldvX = oldV[j  ];   node->oldx = oldP[j  ];
            node->oldvY = oldV[j+1];   node->oldy = oldP[j+1];
            node->oldvZ = oldV[j+2];   node->oldz = oldP[j+2];
            j+=3;
        }
    }

    if (oldP) { free(oldP); oldP=0; }
    if (oldV) { free(oldV); oldV=0; }
}
