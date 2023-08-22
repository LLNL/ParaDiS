
#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"
#include "Timestep.h"

#ifdef USE_ARKODE

#include "ParadisSUNDIALS.h"

/*------------------------------------------------------------------------
 *
 *      Function:    ARKodeIntegrator_RhsFi
 *      Description: Implements the right-hand side function for 
 *                   the ARKode::ARKStep timestep integrator.
 *
 *      Parameters:
 *          IN:  t              Current value of independent variable
 *          IN:  solutionVec    Vector of current solution iterate values
 *          OUT: fVal           Vector of values of g(x)
 *          IN:  vHome          The <home> structure sent to ARKODE cast
 *                              as a (void *).
 *-----------------------------------------------------------------------*/
int ARKodeIntegrator_RhsFi(real8 t, ParaDiS_Vector solutionVec,
                           ParaDiS_Vector fVal, void *vHome)
{
    Home_t *home = (Home_t *) vHome;
       
/* 
 *  Copy node and ghostNode positions from solutionVec vector
 *  into home for velocity calculations 
 */
    ParaDiSVectorToHome_Positions(solutionVec, home, 0);
       
/* 
 *  Calculate nodal force and velocity at the new positions 
 */
    NodeForce(home, FULL);

    int doAll=1;
    int mobIterError = CalcNodeVelocities(home, 0, doAll);

    CommSendVelocity(home);
       
/*
 *  Communicate whether any node was not able to iterate the mobility
 *  function to convergence
 */
#ifdef PARALLEL
    real8  localVals  = (real8) mobIterError;
    real8  globalVals = 0;
       
    MPI_Allreduce(&localVals, &globalVals, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
       
    int globalIterError = globalVals;
#else
    int globalIterError = mobIterError;
#endif
       
/*
 *  If any domain encountered an error iterating inside the mobility
 *  function, return to ARKODE with function failure flag
 */
    if (globalIterError) { return(-1); }
       
/* 
 *  Pack new velocities into fVal vector for ARKode 
 */
    HomeToParaDiSVector_Velocities(home, fVal);
       
    return(0);
}
#endif  // USE_ARKODE
