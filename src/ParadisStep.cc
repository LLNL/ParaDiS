//-------------------------------------------------------------------------
// Function:     ParadisStep
// Description:  This function controls everything needed for a
//               single step of a ParaDiS simulation including
//               force calculations, ghost cell communications,
//               node migration, dynamic load balance, output
//               generation, etc.
//
// Last Modified: 04/08/2008 gh - Added support for Vanadium non-linear
//                                mobility function.
//-------------------------------------------------------------------------

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "DisplayC.h"
#include "Comm.h"
#include "Mobility.h"
#include "Decomp.h"
#include "QueueOps.h"
#include "SMN.h"
#include "CrossSlip.h"
#include "Remesh.h"
#include "Timestep.h"
#include "Spectral.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

// By default, there are no runtime checks to see if all of
// the dislocations have annihilated themselves.  To enable
// a check with an abort if it happens, simply define the
// DEBUG_CHECK_FOR_ZERO_SEG value below to 1 rather than zero.

#define DEBUG_CHECK_FOR_ZERO_SEG 0

//-------------------------------------------------------------------------
//
// Function:    ApplyDeltaStress
// Description: Increment the force/vel for all native nodes based
//              on the change in applied stress (deltaStress)
//              during a timestep.  Since the applied stress plays
//              no impact on the segment/segment interactions this
//              is much less expensive than recomputing all the n^2
//              segment interactions.
//
// Special note for Eshelby inclusion support...
//      If we are using a mobility function that handles segments intersecting
//      eshelby inclusions in a different manner than it handles other segments
//      the mobility function requires knowledge of all segment/particle
//      intersections, but that info is not available at this point.  This
//      means we have to skip this function because we cannot recalculate proper
//      velocities.
//-------------------------------------------------------------------------

#ifndef ESHELBY
static void ApplyDeltaStress(Home_t *home, real8 deltaStress[3][3])
{
    Param_t *param = home->param;

    // Loop over all native nodes

    int numNodes = home->newNodeKeyPtr;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int nodeIndex=0; nodeIndex < numNodes; nodeIndex++)
    {
        Node_t  *node = home->nodeKeys[nodeIndex];

        if (node == (Node_t *)NULL) {
            continue;
        }

        real8 x1 = node->x;
        real8 y1 = node->y;
        real8 z1 = node->z;

        real8 f1Total[3] = { 0.0 };
        VECTOR_ZERO(f1Total);

        // For each node, recalculate forces for all the node's arms
        // that are either owned by this node, or terminate non-locally.

        int numSegs = node->numNbrs;

        for (int segIndex=0; segIndex < numSegs; segIndex++)
        {
            Node_t *nbrNode = GetNeighborNode(home, node, segIndex);

            if (nbrNode == (Node_t *)NULL) {
                printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                continue;
            }

            int nbrIsLocal = (nbrNode->myTag.domainID == home->myDomain);

            if (nbrIsLocal) {
                if (OrderNodes(node, nbrNode) >= 0) {
                    continue;
                }
            }

            real8 bx = node->burgX[segIndex];
            real8 by = node->burgY[segIndex];
            real8 bz = node->burgZ[segIndex];

            real8 dx = nbrNode->x - x1;
            real8 dy = nbrNode->y - y1;
            real8 dz = nbrNode->z - z1;

            ZImage(param, &dx, &dy, &dz) ;

            real8 x2 = x1 + dx;
            real8 y2 = y1 + dy;
            real8 z2 = z1 + dz;

            real8 f1[3] = { 0.0 };
            real8 f2[3] = { 0.0 };

            ExtPKForce(deltaStress, bx, by, bz, x1, y1, z1, x2, y2, z2, f1, f2);

            // Note: the segment force updates (via AddtoArmForce())
            // are technically safe here since each segment's forces
            // will only be updated by a single thread.  The nodal
            // forces, though, may be adjusted by multiple threads,
            // so some locking is necessary.  The AddtoNodeForce()
            // function itself handles the locking, so no explicit
            // locking is needed here.

            VECTOR_ADD(f1Total, f1);
            AddtoArmForce(node, segIndex, f1);

            if (nbrIsLocal)
            {
                int nbrArmID = GetArmID(nbrNode, node);

                AddtoNodeForce(nbrNode, f2);
                AddtoArmForce (nbrNode, nbrArmID, f2);
            }
        }

        // We've accumulated all the forces for this node, so now
        // update the nodal force and recalculate the velocity.

        MobArgs_t mobArgs;

        AddtoNodeForce(node, f1Total);
        EvaluateMobility(home, node, &mobArgs);
    }

    return;
}
#endif   // !ESHELBY


//-------------------------------------------------------------------------
// Function:    ReevaluateForces
// Description: Look for any local nodes whose force/velocity data
//              have been flagged as obsolete.  For all such nodes
//              recompute the forces and velocities.
//-------------------------------------------------------------------------
void ReevaluateForces(Home_t *home)
{
    for (int i=0; i < home->newNodeKeyPtr; i++)
    {
        Node_t *node = home->nodeKeys[i];

        if (node && (node->flags & NODE_RESET_FORCES) )
        {
            MobArgs_t mobArgs;

            SetOneNodeForce(home, node);
            EvaluateMobility(home, node, &mobArgs);

            node->flags &= (~NODE_RESET_FORCES);
        }
    }

    return;
}

#ifdef FIX_PLASTIC_STRAIN
//-------------------------------------------------------------------------
// ResetPlasticStrain()
//
// Set the plastic strain correction to zero before handling
// the topological operations.
//-------------------------------------------------------------------------

static void ResetPlasticStrain(Home_t *home)
{
    Param_t *param = home->param;

    for (int i=0; i<6; i++)
    {
        param->delpStrainCorr[i] = 0.0;
        param->delpSpinCorr  [i] = 0.0;
    }
}

//-------------------------------------------------------------------------
// Function:    CorrectPlasticStrain
// Description: Correct the plastic strain increments for segments
//              that have undergone topological operations.
//-------------------------------------------------------------------------
static void CorrectPlasticStrain(Home_t *home)
{
    Param_t *param = home->param;

#ifdef PARALLEL
    // Correct the plastic strain increments for segments that have
    // undergone topoligical operations. We've calculated processor
    // specific values, now sum the delta strain from all processors,
    // and accumulate into net plastic strain.

    real8 gdstn[6] = { 0.0 };
    real8 gdspn[6] = { 0.0 };

    MPI_Allreduce(param->delpStrainCorr, gdstn, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(param->delpSpinCorr  , gdspn, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int i=0; i<6; i++)
    {
        param->delpStrainCorr[i] = gdstn[i];
        param->delpSpinCorr  [i] = gdspn[i];
    }
#endif
    for (int i=0; i<6; i++)
    {
        param->delpStrain[i] += param->delpStrainCorr[i];
        param->delpSpin  [i] += param->delpSpinCorr  [i];
    }
}
#endif   // FIX_PLASTIC_STRAIN


//-------------------------------------------------------------------------
// ParadisStep()
//
// This is the main ParaDiS cycle/step routine....
//-------------------------------------------------------------------------

void ParadisStep(Home_t *home)
{
    Param_t  *param = home->param;

    home->cycleForceCalcCount = 0;

    // If this step is an initial load-balance-only step just
    // perform the minimal work needed to estimate per-process
    // load, shift boundaries, and migrate nodes among processors.

    if (param->numDLBCycles > 0)
    {
        // Note:  When param->numDLBCycles > 0, NodeForce() assumes
        // the cycle is a DLB-only cycle and only counts the number
        // of force calcs that would be done without actually calculating
        // any forces.

        NodeForce(home, FULL);
        Rebalance(home, DLB_USE_FORCECALC_COUNT);

        // Any time the boundaries are changed, we need to migrate
        // nodes to their new owning domains and go through all the
        // ghost node communication stuff.

        Migrate          (home);
        RecycleGhostNodes(home);
        SortNativeNodes  (home);
        CommSendGhosts   (home);

        return;
    }

#ifndef NO_XWINDOW
    while ( WinIsPaused() ) { sleep(1); }
#endif

    // TestSpectral(home);

    // Calculate the net charge tensor for each cell (includes global comm)

    CellCharge(home);

    // Calculate new force and velocity data for either all nodes or a
    // selected subset, and distribute the new data out to neighboring
    // domains.  The first cycle we'll recalculate all forces and velocities
    // We do this to get an initial estimate of forces on the first cycle,
    // After that, we only need to recompute values for nodes that were
    // involved in topological changes the previous step.
    //
    // NOTE: When we only do force calcs for a subset of the nodes, and
    //       we're using the special mobility for segments that intersect
    //       inclusions, we can only update velocities for nodes flagged
    //       for force updates, NOT their neighbors as well.  This is
    //       because we may not have the full segment/particle intersection
    //       info for the neighboring nodes.  In all other cases, we
    //       want to update all velocities.

    static int firstTime       = 1;
           int nodeForceLevel  = 0;
           int doAllVelocities = 0;
           int zeroOnError     = 0;

    if (firstTime)
    {
        nodeForceLevel  = FULL;
        zeroOnError     = 1;
        doAllVelocities = 1;
        firstTime       = 0;
    }
    else
    {
        nodeForceLevel  = PARTIAL;
        zeroOnError     = 0;
#ifdef ESHELBY
        doAllVelocities = 0;
#else
        doAllVelocities = 1;
#endif
    }

    NodeForce         (home, nodeForceLevel);
    CalcNodeVelocities(home, zeroOnError, doAllVelocities);
    CommSendVelocity  (home);

    // Invoke the selected time step integration method.  The
    // selected method will calculate the time step as well as
    // move nodes to their correct locations and communicate
    // the new nodal force/velocity data to neighboring domains.
    //
    // The 'trapezoid' method used to be specified as 'backard-euler', so
    // if the integration method is unrecognized, just use the trapezoid
    // method as default

         if (StrEquiv(param->timestepIntegrator, "forward-euler"    )) { ForwardEulerIntegrator   (home); }
    else if (StrEquiv(param->timestepIntegrator, "trapezoid"        )) { TrapezoidIntegrator      (home); }
    else if (StrEquiv(param->timestepIntegrator, "trapezoid-multi"  )) { TrapezoidIntegratorMulti (home); }
    else if (StrEquiv(param->timestepIntegrator, "trapezoid-kinsol" )) { TrapezoidIntegratorKINSOL(home); }
    else if (StrEquiv(param->timestepIntegrator, "arkode"           )) { ARKodeIntegrator         (home); }
    else                                                               { TrapezoidIntegrator      (home); }    //< (default)


    // In some simulations, it is necessary to recalculate and distribute
    // the glide plane information for segments after the nodes have been
    // repositioned.  Do so now if needed.

    ResetGlidePlanes(home);

    // Increment the per-burgers vector density gain/loss with
    // changes for this cycle.  This must be done immediately
    // after timestep integration!
    //
    // Note: This is currently only applicable to BCC simulations.

    GetDensityDelta(home);

    // Calculate the new plastic strain.

    DeltaPlasticStrain(home);

    // The call to GenerateOutput will update the time and cycle counters,
    // determine if any output needs to be generated at this stage, and
    // call the appropriate I/O functions if necessary.

    GenerateOutput(home, STAGE_CYCLE);

    // Before doing topological changes, set flags indicating any
    // nodes exempt from topological changes.  These flags are used
    // in both splitting multi-arm nodes and collisions, so this
    // function should be invoked before either of those items are done.

    InitTopologyExemptions(home);

    // Now do all the topological changes from segment interactions
    // (collisions, multinode splitting)...  Clear the list of local
    // operations that will be sent to the remote domains for processsing,
    // then split any multi-arm nodes that need splitting, cross slip
    // nodes (as needed/allowed), handle all local collisions, then
    // send remote nodes the list of ops needed to keep their data in sync.

    ClearOpList(home);
    SortNodesForCollision(home);

#ifdef FIX_PLASTIC_STRAIN
    // Set the plastic strain correction to zero before handling
    // the topological operations.

    ResetPlasticStrain(home);
#endif

    if (param->useParallelSplitMultiNode) { SplitMultiNodesParallel(home); }
    else                                  { SplitMultiNodes        (home); }

#ifdef _ARLFEM
    SplitSurfaceNodes(home);
#endif

    // Call a generic cross-slip dispatch function that will invoke
    // (if necessary) the cross-slip function appropriate to the
    // type of material in use.

    CrossSlip(home);

    // Search for dislocation segments in close proximity to each other
    // and if necessary handle any collision between them.

    HandleCollisions(home);

    // Check for any segments that should be decomposed (splintered)
    // into multiple segments of a different burgers vector.

    SplinterSegments(home);

    // For HCP simulations, we may want to split segments with <c+a>
    // burgers into two segments wth a <c> and an <a> burgers vector
    // respectively, and cross-slip one to a new plane...

    HCP_CA_Split(home);

#ifdef ESHELBY
    // Handle collsions between dislocations and particles

    if (home->param->RefineSegCollision == 1)
        HandleParticleCollisions(home);
#endif

#ifdef PARALLEL
#ifdef SYNC_TIMERS
    TimerStart (home, POST_COLLISION_BARRIER);
    MPI_Barrier(MPI_COMM_WORLD);
    TimerStop  (home, POST_COLLISION_BARRIER);
#endif
#endif
    TimerStart    (home, COL_SEND_REMESH);
    CommSendRemesh(home);
    TimerStop     (home, COL_SEND_REMESH);

    TimerStart    (home, COL_FIX_REMESH );
    FixRemesh     (home);
    TimerStop     (home, COL_FIX_REMESH );

    // For debugging purposes only, the following check has been added
    // to force an immediate abort if we detect an unconstrained node
    // at which the burgers vector is not conserved.  This code can
    // be removed once all the above topology modifying functions
    // have been verified to play well together.

    CheckBurgConservation(home, "After FixRemesh()");

#ifdef _ARLFEM
    FixSurfaceTopology(home);
#endif

    // Under certain circumstances, parallel topological changes can
    // create double links between nodes; links which can not be detected
    // until after FixRemesh() is called... so, a quick check has to be
    // done to clean up these potential double-links here, or they will
    // cause problems later on.  Should only have to check nodes local
    // to this domain.

    for (int i=0; i < home->newNodeKeyPtr; i++)
    {
        Node_t *node = home->nodeKeys[i];

        if (node)
        {
            RemoveDoubleLinks(home, node, 0);
            node->flags &= ~NODE_CHK_DBL_LINK;
        }
    }

    // If memory debugging is enabled, run a consistency check on all
    // allocated memory blocks looking for corruption in any of the
    // block headers or trailers.

#ifdef DEBUG_MEM
    ParadisMemCheck();
#endif
    // Invoke mesh coarsen/refine

    Remesh(home);

#ifdef FIX_PLASTIC_STRAIN
    // Correct the plastic strain increments for segments that have
    // undergone topoligical operations.

    CorrectPlasticStrain(home);
#endif

    // Define load curve and calculate change in applied stress this cycle

    real8  deltaStress[3][3] = { 0.0 };

    LoadCurve(home, deltaStress);

    // This is only needed when we do force calcs for only
    // a subset of the nodes at the beginning of the timestep.  It will
    // adjust the nodal forces based on the current delta stress and
    // recalculate the nodal velocities so we have more accurate values
    // when we enter the timestep integrator at the beginning of the next
    // cycle.

#ifndef ESHELBY
    ApplyDeltaStress(home, deltaStress);
#endif

    // If necessary, use the current load data to generate a new
    // domain decomposition to rebalance the workload among the
    // processors.

    Rebalance(home, DLB_USE_WALLCLK_TIME);

    // Send any nodes that have moved beyond the domain's
    // boundaries to the domain the node now belongs to.

    Migrate(home);

    // Recycle all the ghost nodes: move them back to the free Queue

    RecycleGhostNodes(home);

    // Sort the native nodes and Eshelby inclusions into their proper
    // subcells.

    SortNativeNodes(home);

    // Communicate ghost cells to/from neighbor domains

    CommSendGhosts(home);

#ifdef NAN_CHECK
    // For debug only:  Abort if any of the nodes have position or
    // velocity values that are NaNs or infinites.  Be sure to do this
    // before we write the restart files so we don't get bad data
    // in the restart.

    CheckForNANS(home, "In ParadisStep() after ghost commm");
#endif

    // If memory debugging is enabled, run a consistency check on all
    // allocated memory blocks looking for corruption in any of the
    // block headers or trailers.

#ifdef DEBUG_MEM
    ParadisMemCheck();
#endif

    CheckMemUsage(home, "ParadisStep-complete");

    // Zero out the count of force calculations done this cycle
    // so the load-balancing in the next step is based on accurate
    // values.

    home->cycleForceCalcCount = 0;

    // For Debugging and testing only...

#if DEBUG_CHECK_FOR_ZERO_SEG
    CheckForEmptySimulation(home);
#endif

    return;
}
