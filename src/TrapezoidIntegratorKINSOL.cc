//---------------------------------------------------------------------------------------------------------
// Module:      TrapezoidIntegratorKINSOL.cc
// Description: Implements a numerical timestep integrator using
//              the Trapezoid integration method.  This version
//              of the integrator is implemented using KINSOL,
//              with either an Anderson accelerated fixed point
//              nonlinear solver, or a Newton-Krylov nonlinear solver.
//
// Includes private functions:
//    AdvanceAllNodes()
//    PreserveNodalData()
//---------------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Timestep.h"
#include "DeltaTime.h"

#ifdef USE_KINSOL
#include "ParadisSUNDIALS.h"
#include <kinsol/kinsol.h>
#include <sunlinsol/sunlinsol_spgmr.h>

#ifdef _ARLFEM
#include "FEM.h"
#endif

#ifdef DEBUG_TIMESTEP
//---------------------------------------------------------------------------------------------------------
// Function:    PrintLimitingNode
// Description: Find the local node (if any) that is responsible
//              for the tiemstep being cut and print out the
//              nodal info.
//
// Parameters:
//     IN: globalErrMax  Maximum positioning error norm from all nodes
//                       in the system.  This value is obtained from
//                       via KINGetFuncNorm().
//     IN: newDT         Timestep size which when adjusted by
//                       <dtDecrementFact> will be the next attempted
//                       timestep size.
//---------------------------------------------------------------------------------------------------------

static void PrintLimitingNode(Home_t *home, real8 globalErrMax, real8 newDT)
{
        real8   localMax  = 0.0;
        real8   globalMax = 0.0;

        real8   errMax    = 0.0;
        Node_t *errNode   = (Node_t *)NULL;

        // Loop through the local nodes and find the max positioning error norm.

        for (int i=0; i < home->newNodeKeyPtr; i++)
        {
            Node_t *node = home->nodeKeys[i];

            if (node)
            {
                real8 oldx = node->oldx;
                real8 oldy = node->oldy;
                real8 oldz = node->oldz;

                PBCPOSITION(home->param, node->x, node->y, node->z, &oldx, &oldy, &oldz);

                localMax = errMax;

                errMax = MAX(errMax, fabs(node->x - oldx - ((node->vX+node->currvX)*0.5*newDT)));
                errMax = MAX(errMax, fabs(node->y - oldy - ((node->vY+node->currvY)*0.5*newDT)));
                errMax = MAX(errMax, fabs(node->z - oldz - ((node->vZ+node->currvZ)*0.5*newDT)));

                if (errMax > localMax) { errNode=node; }
            }
        }

        // FIX ME!
        // For some reason, the max error returned from the KINSOL
        // call and passed to this function does not match the maximum
        // error calculated by this function even though they *should*
        // be doing the samee calculation.  So, for now, ignore the
        // globalErrMax passed into this function, use the locally
        // generated error values and do an Allreduce to find the
        // global max from among all tasks.

#ifdef PARALLEL
        localMax = errMax;

        MPI_Allreduce(&localMax, &globalMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
        globalMax = errMax;
#endif

        // If the local max matches the global max, the associated node
        // is the one limiting the timestep.

        if (localMax == globalMax)
        {
            printf("  Cut timestep for (%d,%d):  errMax %e, newDT %e\n",
                   errNode->myTag.domainID, errNode->myTag.index,
                   localMax, newDT * home->param->dtDecrementFact);
            PrintNode(errNode);
            Timestep_Cut_Save(home,errNode);
        }

        return;
}
#endif  // ifdef DEBUG_TIMESTEP


// Limit_Node_Velocity()
//
// It's possible that velocities computed inside the mobility functions can result in non-physical
// node velocities. This routine will limit the computed velocities within this module to the
// shear sound velocity provided via param.
//
// If the shear velocity is greater than zero and the magnitiude of the velocity vector exceeds
// the shear velocity, the velocity is rescaled to a shear velocity vector directed along the
// original vector.
//---------------------------------------------------------------------------------------------------------

#ifdef LIMIT_NODE_VELOCITY
static void Limit_Node_Velocity
(
         real8  & vx  ,   ///< input velocity (x) (b/s) (potentially limited)
         real8  & vy  ,   ///< input velocity (y) (b/s) (potentially limited)
         real8  & vz  ,   ///< input velocity (z) (b/s) (potentially limited)
   const real8    bmag,   ///< magnitude of the burgers vector (m)
   const real8    vmax    ///< shear sound velocity (m/s)
)
{
   if (vmax>0.0)
   {
      real8 vmag = sqrt(vx*vx + vy*vy + vz*vz) * bmag;   // vmag = magnitude of node velocity (m/s)
      real8 vs   = ( (vmag>vmax) ? (vmax/vmag) : 1.0 );  // vs   = scale factor

            vx *= vs;                                    // apply scale factor
            vy *= vs;                                    //
            vz *= vs;                                    //
   }
}
#else
#define Limit_Node_Velocity(a,b,c,d,e)
#endif

// Advance_Node()
//
// Will update the node's position given a timestep (dt). This method will average
// the current and previous node velocities for a given timestep.
//
// Additional parameters are provided to adjust for periodic boundary conditions (PBC)...
//
//    cx,cy,cz - the center of the simulation box <xyz>
//    lx,ly,lz - the size   of the simulation box <xyz>
//    sx,sy,sz - PBC adjustment reciprocals       <xyz>
//
//               If PBC for a particular axis is active, the PBC reciprocal will be
//               the reciprocal of the simulation size, otherwise zero.
//               e.g  sx = ( pbcx && (lx>0.0) ? (1.0/lx) : 0.0 );
//                    sy = ( pbcy && (ly>0.0) ? (1.0/ly) : 0.0 );
//                    sz = ( pbcz && (lz>0.0) ? (1.0/lz) : 0.0 );
//               By setting  sx,sy,sz to the reciprocals, the rint() correction will be applied.
//---------------------------------------------------------------------------------------------------------

static void Advance_Node
(
         Node_t *node,                                 ///< points to node to advance
   const real8   dt,                                   ///< next time step
   const real8   cx, const real8 cy, const real8 cz,   ///< simulation center <xyz>
   const real8   lx, const real8 ly, const real8 lz,   ///< simulation size   <xyz>
   const real8   sx, const real8 sy, const real8 sz,   ///< pbc adjustment reciprocals <xyz>
   const real8   bmag,                                 ///< magnitude of the burgers vector (m)
   const real8   vmax                                  ///< shear sound velocity            (m/s)
)
{
   if (node)
   {
      // If we don't have a value for the previous velocity, assume
      // previous velocity was same as current velocity.

      if (    (node->oldvX == 0.0)
           && (node->oldvY == 0.0)
           && (node->oldvZ == 0.0) )
      {
         node->oldvX = node->currvX;
         node->oldvY = node->currvY;
         node->oldvZ = node->currvZ;
      }

      // compute node velocity...

      real8 vx = (node->currvX + node->oldvX)/2.0;  // vx = average of current and previous velocity (x)
      real8 vy = (node->currvY + node->oldvY)/2.0;  // vy = average of current and previous velocity (y)
      real8 vz = (node->currvZ + node->oldvZ)/2.0;  // vz = average of current and previous velocity (z)

      Limit_Node_Velocity(vx,vy,vz, bmag,vmax);      // limit velocity to shear sound velocity

      real8 dx = (vx*dt);
      real8 dy = (vy*dt);
      real8 dz = (vz*dt);

      real8 x  = node->oldx + dx;
      real8 y  = node->oldy + dy;
      real8 z  = node->oldz + dz;

      // adjust the updated position for periodic boundary...

      x -= rint((x-cx)*sx)*lx;
      y -= rint((y-cy)*sy)*ly;
      z -= rint((z-cz)*sz)*lz;

      node->x = x;
      node->y = y;
      node->z = z;
   }
}

// Advance_Nodes()
//
// Advances an array of nodes or node pointers. Basically - sets up the parallel for-loop that
// drives the nodal advance.
//---------------------------------------------------------------------------------------------------------

static void Advance_Nodes
(
         Node_t **narr,                                 ///< array of node pointers (can be sparse)
   const int      ncnt,                                 ///< number of nodes
   const real8    dt,                                   ///< dt = next time step
   const real8    cx, const real8 cy, const real8 cz,   ///< simulation center <xyz>
   const real8    lx, const real8 ly, const real8 lz,   ///< simulation size   <xyz>
   const real8    sx, const real8 sy, const real8 sz,   ///< pbc adjustment reciprocals <xyz>
   const real8    bmag,                                 ///< magnitude of the burgers vector (m)
   const real8    vmax                                  ///< shear sound velocity            (m/s)
)
{
   if (narr && (ncnt>0))
   {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i=0; (i<ncnt); i++)
         Advance_Node(narr[i], dt, cx,cy,cz, lx,ly,lz, sx,sy,sz, bmag,vmax);
   }
}

//---------------------------------------------------------------------------------------------------------
//  Function:     Advance_Nodes
//  Description:  Reposition all nodes (local and ghost) based
//                on the old/current nodal velocities and time deltas.
//
//                This function assumes current nodal positions and velocities
//                have already been preserved in the old* variables by a call
//                to PreserveNodalData().
//
//  Given the following:
//
//      currDT = desired delta time this timestep
//      oldDT  = delta time used in the previous timestep
//      currV  = velocity on entry to timestep integrator
//      oldV   = velocity on entry to timestep integrator on previous step
//      currP  = current nodal position
//      newP   = initial guess at position node will end up in after this timestep.
//
//      The new positions are calculates as:
//
//          newP = currP + 0.5 * (currV + oldV) * currDT;
//---------------------------------------------------------------------------------------------------------

static void Advance_Nodes(Home_t *home, real8 oldDT)
{
   Param_t *param = home->param;

   param->realdt = param->deltaTT;

   const real8 dt    = param->realdt;
   const real8 bmag  = param->burgMag;           // bmag = magnitude of the burgers vector (m)
   const real8 vmax  = param->shearVelocity;     // vmax = shear sound velocity (m/s)

   // set simlulation center...

   const real8 cx    = (param->maxSideX + param->minSideX)/2.0;
   const real8 cy    = (param->maxSideY + param->minSideY)/2.0;
   const real8 cz    = (param->maxSideZ + param->minSideZ)/2.0;

   // set simulation size...

   const real8 lx    = (param->maxSideX - param->minSideX);
   const real8 ly    = (param->maxSideY - param->minSideY);
   const real8 lz    = (param->maxSideZ - param->minSideZ);

   // set simulation periodic boundary reciprocals...

   const real8 sx    = ( (param->xBoundType == Periodic) && (lx>0.0) ? (1.0/lx) : 0.0 );
   const real8 sy    = ( (param->yBoundType == Periodic) && (ly>0.0) ? (1.0/ly) : 0.0 );
   const real8 sz    = ( (param->zBoundType == Periodic) && (lz>0.0) ? (1.0/lz) : 0.0 );

   // advance the native and ghost node arrays...

   Advance_Nodes(home->nodeKeys     , home->newNodeKeyPtr , dt, cx,cy,cz, lx,ly,lz, sx,sy,sz, bmag,vmax);
   Advance_Nodes(home->ghostNodeList, home->ghostNodeCount, dt, cx,cy,cz, lx,ly,lz, sx,sy,sz, bmag,vmax);
}

// Preserve_Node()
//
// Will preserve various node data elements specified via the items argument.
//
// The items argument is a bitwise 'or' of the following constants...
//    NODE_POSITION : saves the current node position into the old position
//    NODE_CURR_VEL : saves the current velocity without wiping out the copy of the previous velocity
//    NODE_OLD_VEL  : saves the old velocity (normally done at the end of the timestep integrator
//---------------------------------------------------------------------------------------------------------

static void Preserve_Node (Node_t *node, const int items)
{
   if (node)
   {
      if (items & NODE_POSITION)
      {
         node->oldx   = node->x;
         node->oldy   = node->y;
         node->oldz   = node->z;
      }

      if (items & NODE_CURR_VEL)
      {
         node->currvX = node->vX;
         node->currvY = node->vY;
         node->currvZ = node->vZ;
      }

      if (items & NODE_OLD_VEL)
      {
         node->oldvX  = node->currvX;
         node->oldvY  = node->currvY;
         node->oldvZ  = node->currvZ;
      }
   }
}

// Preserve_Nodes()
//
// Preserves the elements of an array of node pointers.
//---------------------------------------------------------------------------------------------------------

static void Preserve_Nodes
(
         Node_t **narr ,  ///< array of node pointers (can be sparse)
   const int      ncnt ,  ///< number of nodes
   const int      items   ///< items to preserve (see above)
)
{
   if (narr && (ncnt>0))
   {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i=0; (i<ncnt); i++)
         Preserve_Node(narr[i],items);
   }
}

//---------------------------------------------------------------------------------------------------------
// Function:    TrapezoidIntegratorKINSOL
//
// Description: Implements a numerical timestep integrator using
//              the Trapezoid integration method.  This version
//              of the integrator is implemented using KINSOL,
//              a general-purpose nonlinear system solver based on
//              Newton-Krylov solver technology
//
//              Note: This function assumes that the nodal
//              force/velocity data is accurate for the current
//              positions of the nodes on entry to the routine.
//---------------------------------------------------------------------------------------------------------

void TrapezoidIntegratorKINSOL(Home_t *home)
{
        int      kinsolStatus   = 0;
        long     numSolverIters = 0;
        real8    globalErrMax   = 0.0;

        Param_t *param = home->param;

        // update time step values

        real8 oldDT = param->deltaTT;
        real8 newDT = MIN(param->maxDT, param->nextDT);

        if (newDT <= 0.0) { newDT = param->maxDT; }

        param->deltaTT = newDT;

        // Preserve certain items of nodal data for later use.

        Preserve_Nodes (home->nodeKeys     , home->newNodeKeyPtr , NODE_POSITION | NODE_CURR_VEL );
        Preserve_Nodes (home->ghostNodeList, home->ghostNodeCount, NODE_POSITION | NODE_CURR_VEL );

        // If the user provided the name of a log file for use with
        // KINSOL, open it now.  Note: the log is only written from
        // a single task...

        FILE  *kinsolLog=0;

        if ((home->myDomain == 0) && (param->KINSOL_LogFile[0] != 0)) {
            if ((kinsolLog = fopen(param->KINSOL_LogFile, "a")) == NULL) {
                printf("Warning: can't open kinsol log file %s.\n", param->KINSOL_LogFile);
            }
        }

        if ((home->newNodeKeyPtr + home->ghostNodeCount <= 0))
        {
            Fatal("TrapezoidIntegratorKINSOL: no nodes!\n");
        }

        // Create new vector of solution values to send to KINSOL.
        // Primary solution is a vector of nodal positions.  Values are
        // copied from the <home> structure into this solution vector, solutionVec

        ParaDiS_Vector solutionVec = NewParaDiSVector  (home);
        ParaDiS_Vector ones        = CloneParaDiSVector(solutionVec);

        P_VConst(1.0e0, ones);

        // General KINSOL setup.  See KINSOL User's Guide for details on the setup functions.

        real8  fNormTol  = param->rTol;
        void  *kinsolMem = KINCreate();

        // Set up the solver; either Anderson accelerated fixed point
        // nonlinear solver or a Newton-Krylov nonlinear solver.  The
        // Anderson accelerated fixed-point solver is the default.

        SUNLinearSolver kinsolLS=NULL;

        int flag=0;
        flag = KINSetErrFile    (kinsolMem, kinsolLog);
        flag = KINSetInfoFile   (kinsolMem, kinsolLog);
        flag = KINSetUserData   (kinsolMem, home);
        flag = KINSetNumMaxIters(kinsolMem, param->KINSOL_MaxIterations);
        flag = KINSetFuncNormTol(kinsolMem, fNormTol);
        flag = KINSetPrintLevel (kinsolMem, param->KINSOL_PrintLevel);

        if (param->KINSOL_UseNewton)
        {
            kinsolLS = SUNLinSol_SPGMR(solutionVec, PREC_NONE, param->KINSOL_MaxLinIterations);

            KINInit            (kinsolMem, Newton_Fcn, solutionVec);
            KINSetLinearSolver (kinsolMem, kinsolLS, NULL);
            KINSetEtaForm      (kinsolMem, KIN_ETACONSTANT);
            KINSetEtaConstValue(kinsolMem, param->KINSOL_EtaVal);
        }
        else
        {
            // (anderson accelerated)
            flag = KINSetMAA(kinsolMem, param->KINSOL_NumPriorResiduals);

            KINInit(kinsolMem, FixedPt_Fcn, solutionVec);
        }

        // Loop until we converge on a time step.  First step is to
        // use the current positions and velocities and reposition the
        // nodes to where they would be after the suggested delta time.

        int convergent = 0;
        int incrDelta  = 1;

        while (!convergent)
        {
            // Advance all nodes from their previous positions to new positions
            // based on their current velocities and the delta T being tested.
            // This includes all ghost nodes, so nodal velocities must previously
            // have been distributed to neighboring domains.

            Advance_Nodes(home, oldDT);

            // Copy node and ghostNode positions from <home> into solutionVec

            HomeToParaDiSVector_Positions(home,solutionVec,1);

            // Initialize error

            globalErrMax = 0.0;

            // Invoke the appropriate solver

            if (param->KINSOL_UseNewton)
                kinsolStatus = KINSol(kinsolMem, solutionVec, KIN_NONE, ones, ones);
            else
                kinsolStatus = KINSol(kinsolMem, solutionVec, KIN_FP, ones, ones);

            // Get error and iteration count

            flag = KINGetFuncNorm(kinsolMem, &globalErrMax);
            flag = KINGetNumNonlinSolvIters(kinsolMem, &numSolverIters);

            // Check if KINSOL was able to converge on a solution

            if (    (kinsolStatus == KIN_SUCCESS)
                 || (kinsolStatus == KIN_INITIAL_GUESS_OK) )
            {
                convergent = 1;

#ifdef DEBUG_PARADIS_VECTOR
                printf("Converged Solution Vector \n");
                PrintParaDiSVector(solutionVec);
#endif

                // Copy node and ghostNode positions from solutionVec vector into <home>

                ParaDiSVectorToHome_Positions(solutionVec,home,1);
            }

            if (numSolverIters > MIN((param->KINSOL_MaxIterations-1), 1))
            { incrDelta = 0; }

            // If there is convergence, we've got a good delta T, otherwise
            // cut the delta T by a configured factor and try again.

            if (!convergent)
            {
#ifdef DEBUG_TIMESTEP
                PrintLimitingNode(home, globalErrMax, newDT);
#endif
                newDT *= param->dtDecrementFact;
                param->deltaTT = newDT;

                if ((newDT < 1.0e-20) && (home->myDomain == 0)) {
                    Fatal("TrapezoidIntegratorKINSOL(): Timestep has dropped below minimal threshold to %e.  Aborting!", newDT);
                }

            }  // if (!convergent)

        }  // while (!convergent)


        // Automatically increment timestep if convergence was reached
        // on the very first iteration of the above loop.
        //
        // If variable timestep adjustments are enabled, calculate an
        // adjustment factor based on the maximum allowed timestep increment
        // and the maximum error found above.  If variable timestep
        // adjustments are not enabled, adjust the timestep by the
        // maximum permitted factor.

        param->deltaTT   = newDT;
        param->realdt    = newDT;
        param->timeStart = param->timeNow;

        if (incrDelta)
        {
            if (param->dtVariableAdjustment)
            {
                real8 tmp1   = pow(param->dtIncrementFact, param->dtExponent);
                real8 tmp2   = globalErrMax/param->rTol;
                real8 tmp3   = 1.0 / param->dtExponent;
                real8 tmp4   = pow(1.0/(1.0+(tmp1-1.0)*tmp2), tmp3);
                real8 factor = param->dtIncrementFact * tmp4;

                param->nextDT = MIN(param->maxDT, newDT*factor);
            }
            else
            {
                param->nextDT = MIN(param->maxDT, newDT*param->dtIncrementFact);
            }

        } else {  // incrDelta==0
            param->nextDT = newDT;
        }

        // Copy the nodal velocities that existed on entry to the timestep integrator.

        Preserve_Nodes (home->nodeKeys     , home->newNodeKeyPtr , NODE_OLD_VEL );
        Preserve_Nodes (home->ghostNodeList, home->ghostNodeCount, NODE_OLD_VEL );

#ifdef _FEM
        // If we're using the FEM code for simulating free surfaces, we
        // have to handle any nodes/segments that have moved onto or past
        // the simulation surfaces.
        //
        // FIX ME!  Do we need to again recalculate forces/velocities for
        //          the nodes that get modified in AdjustNodePosition(), or
        //          can we live with the results until the next cycle???

        AdjustNodePosition(home, 1);
#endif

        // Release any remaining temporary resources allocated for the
        // KINSOL solver and close the log file.

        FreeParaDiSVector(solutionVec);
        FreeParaDiSVector(ones);
        SUNLinSolFree    (kinsolLS);
        KINFree          (&kinsolMem);

        if (kinsolLog) { fclose(kinsolLog); }

        return;
}
#endif  // USE_KINSOL


#ifndef USE_KINSOL
// If we get here, KINSOL was not enabled during compilation so
// just abort with a fatal error since we can't use the selected
// timestep integrator.


void TrapezoidIntegratorKINSOL(Home_t *home)
{
    Fatal("TrapezoidIntegratorKINSOL(): KINSOL integrator not available.\n"
          "You must enable KINSOL by setting the KINSOL_MODE=ON flag in makefile.setup\n");
}

#endif  // USE_KINSOL
