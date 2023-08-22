
#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"

#ifdef USE_KINSOL
#include "ParadisSUNDIALS.h"
#include "Timestep.h"


/*------------------------------------------------------------------------
 *
 *      Function:    Newton_Fcn
 *      Description: Used by the KINSOL timestep integrator to computes
 *                   the nonliner function value for Newton solve.
 *
 *      Parameters:
 *          IN:  solutionVec    Vector of current solution iterate values
 *          OUT: fVal           Vector of values of F(x)
 *          IN:  vHome          The <home> structure
 *
 *-----------------------------------------------------------------------*/
int Newton_Fcn(ParaDiS_Vector solutionVec, ParaDiS_Vector fVal, void *vHome)
{
    int       i, j;
    int       newNodeKeyPtr, ghostNodeCount;
    int       doAll = 1;
    int       mobIterError = 0;
    int       globalIterError = 0;
    long int  active_length;
    real8     halfDT;
    real8     localVals, globalVals;
    realtype *fdata, *sdata;
    Home_t   *home;
    Node_t  **nodeKeys, **ghostNodeList; 
    Param_t  *param;
       
/* 
 *  Get pointer to data and active length 
 */
    sdata = PV_DATA(solutionVec);
    fdata = PV_DATA(fVal);
    active_length = PV_ACTLENGTH(fVal);

/* 
 *  Shortcuts to home data 
 */
    home           = (Home_t *)vHome;
    param          = home->param;
    nodeKeys       = home->nodeKeys;
    newNodeKeyPtr  = home->newNodeKeyPtr;
    ghostNodeList  = home->ghostNodeList;
    ghostNodeCount = home->ghostNodeCount;
       
    halfDT = param->deltaTT * 0.5;

/* 
 *  Copy node and ghostNode positions from solutionVec vector
 *  into home for velocity calculations 
 */
    ParaDiSVectorToHome_Positions(solutionVec, home, 0);

/* 
 *  Calculate nodal force and velocity at the new positions 
 */
    NodeForce(home, FULL);
    mobIterError = CalcNodeVelocities(home, 0, doAll);
    CommSendVelocity(home);

/*
 *  Communicate whether any node was not able to iterate the mobility
 *  function to convergence
 */
#ifdef PARALLEL
    localVals = (real8)mobIterError;
    globalVals = 0;
       
    MPI_Allreduce(&localVals, &globalVals, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
       
    globalIterError = globalVals;
#else
    globalIterError = mobIterError;
#endif
       
/*
 *  If any domain encountered an error iterating inside the mobility
 *  function, return to ARKODE with function failure flag
 */
    if (globalIterError) {
        return(-1);
    }
       
/*
 *  Evaluate Newton function
 */
    j = 0;

    for (i = 0; i < newNodeKeyPtr; i++) {
        real8   ox, oy, oz;
        Node_t *hlocalNode;

        if (nodeKeys[i] == (Node_t *)NULL) {
            continue;
        }

        hlocalNode = nodeKeys[i];
  
        ox = hlocalNode->oldx;
        oy = hlocalNode->oldy;
        oz = hlocalNode->oldz;

        PBCPOSITION(param, sdata[j], sdata[j+1], sdata[j+2], &ox, &oy, &oz);
    
        fdata[j]   = sdata[j]   - ox -
                     (hlocalNode->vX + hlocalNode->currvX) * halfDT;

        fdata[j+1] = sdata[j+1] - oy -
                     (hlocalNode->vY + hlocalNode->currvY) * halfDT;

        fdata[j+2] = sdata[j+2] - oz -
                     (hlocalNode->vZ + hlocalNode->currvZ) * halfDT;
        j += 3;
    }
  
    j = active_length;

    for (i = 0; i < ghostNodeCount; i++) {
        real8   ox, oy, oz;
        Node_t *hghostNode;
    
        hghostNode = ghostNodeList[i];
    
        ox = hghostNode->oldx;
        oy = hghostNode->oldy;
        oz = hghostNode->oldz;
    
        PBCPOSITION(param, sdata[j], sdata[j+1], sdata[j+2], &ox, &oy, &oz);
    
        fdata[j]   = sdata[j]   - ox -
                     (hghostNode->vX + hghostNode->currvX) * halfDT;

        fdata[j+1] = sdata[j+1] - oy -
                     (hghostNode->vY + hghostNode->currvY) * halfDT;

        fdata[j+2] = sdata[j+2] - oz -
                     (hghostNode->vZ + hghostNode->currvZ) * halfDT;    
        j += 3;
    }

    return(0);
}
#endif  /* ifdef USE_KINSOL */
