/*-------------------------------------------------------------------------
 *
 *      Module:       DebugFunctions.c
 *
 *      Description:  This module contains a variety of functions
 *                    that are intended for debug purposes only and
 *                    are not normally invoked by the code.  During
 *                    the debugging process, calls to these functions
 *                    can be inserted in the code as necessary.
 *
 *
 *      Included functions:
 *
 *          CheckBurgConservation()
 *          CheckSegLengths()
 *          CheckForces()
 *          CheckForNANS()
 *          CheckForEmptySimulation()
 *          CheckForUndefinedPlanes()
 *          CheckPlanes_RhomboVa()
 *          Print2x2()
 *          Synchronize()
 *          VerifyConnections()
 *          VerifyLocalConnections()
 *
 *-----------------------------------------------------------------------*/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Segment.h"

/**************************************************************************
 *
 *      Function:     CheckNodeBurgConservation
 *      Description:  Check if the burgers vector is conserved at
 *                    a single node.
 *
 *      Arguments:
 *          IN:    node  pointer to node structure
 *
 *      Returns: 1 if burgers vector is conserved at the node
 *               2 if burgers vector is not conserved but the node
 *                 has either PINNED or SURFACE_NODE constraints
 *               0 if burgerss vector is not conserved and the node
 *                 has neither the PINNED or SURFACE_NODE constraints
 *
 *************************************************************************/
int CheckNodeBurgConservation(Home_t *home, Node_t *node)
{
    int   i;
    real8 burgSum[3];

    if (HAS_ANY_OF_CONSTRAINTS(node->constraint, PINNED_NODE | SURFACE_NODE)) {
        return(2);
    }

    VECTOR_ZERO(burgSum);

    for (i = 0; i < node->numNbrs; i++) {
        burgSum[X] += node->burgX[i];
        burgSum[Y] += node->burgY[i];
        burgSum[Z] += node->burgZ[i];
    }

    if ((fabs(burgSum[X]) > 1e-4) ||
        (fabs(burgSum[Y]) > 1e-4) ||
        (fabs(burgSum[Z]) > 1e-4)) {
        return(0);
    }

    return(1);
}


/**************************************************************************
 *
 *      Function:     CheckBurgConservation
 *      Description:  Verify that the burgers vector is conserved at all
 *                    nodes that are not pinned or surface nodes.  If
 *                    burgers vector is not conserved at an unconstrained
 *                    node, the function aborts with an error.
 *
 *      Arguments:
 *          msg  text to be displayed with the error message when
 *               burgers vector is not conserved.
 *
 *************************************************************************/
void CheckBurgConservation(Home_t *home, const char *msg)
{
        int    i;
        int    numNodes;

        numNodes = home->newNodeKeyPtr;

        for (i = 0; i < numNodes; i++) {
            Node_t *node;

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            if (CheckNodeBurgConservation(home, node) == 0) {

                PrintNode(node);
                Fatal("%s: Burgers vector not conserved on (%d,%d)\n",
                      msg, node->myTag.domainID, node->myTag.index);
            }
        }

        return;
}


/**************************************************************************
 *
 *      Function:     CheckForces
 *      Description:  Check the force and velocity components of each
 *                    node looking for any with values exceeding some
 *                    hard-code thresholds.  If such a node is found,
 *                    the code aborts with an error identifying the offending
 *                    node.
 *
 *      Arguments:
 *          msg  text to be displayed with the error message when
 *               a long segment is found.
 *
 *************************************************************************/
void CheckForces(Home_t *home, const char *msg)
{
        int    i;
        int    numNodes;
        real8  forceThreshold, velThreshold;
        Node_t *node;

        forceThreshold = 1.0e+14;
        velThreshold   = 1.0e+14;

        numNodes = home->newNodeKeyPtr;

        for (i = 0; i < numNodes; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            if ((node->fX > forceThreshold) ||
                (node->fY > forceThreshold) ||
                (node->fZ > forceThreshold)) {
#if 1
                    Fatal("%s:  Node (%d,%d) force = %e %e %e\n",
                           msg, node->myTag.domainID, node->myTag.index,
                           node->fX, node->fY, node->fZ);
#else
                    printf("%s:  Node (%d,%d) force = %e %e %e\n",
                           msg, node->myTag.domainID, node->myTag.index,
                           node->fX, node->fY, node->fZ);
#endif
            }

            if ((node->vX > velThreshold) ||
                (node->vY > velThreshold) ||
                (node->vZ > velThreshold)) {
#if 1
                    Fatal("%s:  Node (%d,%d) velocity = %e %e %e\n",
                           msg, node->myTag.domainID, node->myTag.index,
                           node->vX, node->vY, node->vZ);
#else
                    printf("%s:  Node (%d,%d) velocity = %e %e %e\n",
                           msg, node->myTag.domainID, node->myTag.index,
                           node->vX, node->vY, node->vZ);
#endif
            }

        }

        return;
}


/**************************************************************************
 *
 *      Function:     CheckSegLengths
 *      Description:  Check the length of every segment owned by the
 *                    current domain and print out an error message
 *                    any time a segment exceeding the maximum segment
 *                    lenth <maxSeg> by 10% or more is found.
 *
 *      Arguments:
 *          msg  text to be displayed with the error message when
 *               a long segment is found.
 *
 *************************************************************************/
void CheckSegLengths(Home_t *home, const char *msg)
{
        int    i, j;
        int    numNodes, numNbrs;
        real8  xLen, yLen, zLen, totLen;
        Node_t *node1, *node2;


        numNodes = home->newNodeKeyPtr;

        for (i = 0; i < numNodes; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            numNbrs = node1->numNbrs;

            for (j = 0; j < numNbrs; j++) {

                node2 = GetNeighborNode(home, node1, j);

                if (node2 == (Node_t *)NULL) {
                    printf("%s:  Segment (%d,%d)--(%d,%d) neighbor not found\n",
                           msg, node1->myTag.domainID, node1->myTag.index,
                           node1->nbrTag[j].domainID, node1->nbrTag[j].index);
                    fflush(NULL);
                    continue;
                }

                xLen = node2->x - node1->x;
                yLen = node2->y - node1->y;
                zLen = node2->z - node1->z;

                ZImage(home->param, &xLen, &yLen, &zLen);

                totLen = sqrt(xLen*xLen + yLen*yLen + zLen*zLen);

                if (totLen > (home->param->maxSeg * 1.1)) {
                    printf("%s:  Segment (%d,%d)--(%d,%d) length = %lf\n",
                           msg, node1->myTag.domainID, node1->myTag.index,
                           node2->myTag.domainID, node2->myTag.index, totLen);
                    fflush(NULL);
                }
            }
        }

        return;
}


/**************************************************************************
 *
 *      Function:     CheckForNaNs
 *      Description:  Check for NaNs (or Infs) in the position, force and
 *                    velocity components for each local node.  If a NaN
 *                    or Inf is found, the code will abort.
 *
 *************************************************************************/
void CheckForNANS(Home_t *home, const char *msg)
{
        int    i;
        Node_t *node;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            if (isnan(node->x) ||
                isnan(node->y) ||
                isnan(node->z)) {
                PrintNode(node);
                Fatal("%s: Node position contains NaNs.  Aborting!", msg);
            }

            if (isinf(node->x) ||
                isinf(node->y) ||
                isinf(node->z)) {
                PrintNode(node);
                Fatal("%s: Node position contains Infs.  Aborting!", msg);
            }

            if (isnan(node->vX) ||
                isnan(node->vY) ||
                isnan(node->vZ)) {
                PrintNode(node);
                Fatal("%s: Node velocity contains NaNs.  Aborting!", msg);
            }

            if (isinf(node->vX) ||
                isinf(node->vY) ||
                isinf(node->vZ)) {
                PrintNode(node);
                Fatal("%s: Node velocity contains Infs.  Aborting!", msg);
            }

            if (isnan(node->fX) ||
                isnan(node->fY) ||
                isnan(node->fZ)) {
                PrintNode(node);
                Fatal("%s: Node force contains NaNs.  Aborting!", msg);
            }

            if (isinf(node->fX) ||
                isinf(node->fY) ||
                isinf(node->fZ)) {
                PrintNode(node);
                Fatal("%s: Node force contains NaNs.  Aborting!", msg);
            }

        }

        return;
}


/**************************************************************************
 *
 *      Function:     CheckForEmptySimulation
 *      Description:  Just a sanity check to kill the simulation if we
 *                    end up in a situation with no more dislocations in
 *                    the problem.
 *
 *************************************************************************/
void CheckForEmptySimulation(Home_t *home)
{
        int    i, localNodeCnt, globalNodeCnt;
        Node_t *node;

        localNodeCnt = 0;
        globalNodeCnt = 0;

        for (i = 0; i < home->newNodeKeyPtr; i++) {
            node = home->nodeKeys[i];
            if (node == (Node_t *)NULL) {
                continue;
            } else {
                localNodeCnt++;
            }
        }

#ifdef PARALLEL
        MPI_Reduce(&localNodeCnt, &globalNodeCnt, 1, MPI_INT, MPI_SUM,
                    0, MPI_COMM_WORLD);
#else
        globalNodeCnt = localNodeCnt;
#endif
        if ((home->myDomain == 0) && (globalNodeCnt == 0)) {
            Fatal("All dislocations in simulation have been "
                  "annihilated.\nSimulation terminating NOW!");
        }
}


/**************************************************************************
 *
 *      Function:     CheckForUndefinedPlanes
 *      Description:  Check for any local segment for which the glide
 *                    plane is undefined.  If found, print the nodal
 *                    information and an error message.
 *
 *      Arguments:
 *          msg  text to be displayed with the error message when
 *               a segment with an undefined glide plane is found.
 *
 *************************************************************************/
void CheckForUndefinedPlanes(Home_t *home, const char *msg)
{
        int    i, j;
        int    numNodes, numNbrs;
        Node_t *node;


        numNodes = home->newNodeKeyPtr;

        for (i = 0; i < numNodes; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            numNbrs = node->numNbrs;

            for (j = 0; j < numNbrs; j++) {

                if ((fabs(node->nx[j]) < 1.0e-4) &&
                    (fabs(node->ny[j]) < 1.0e-4) &&
                    (fabs(node->nz[j]) < 1.0e-4)) {
                    PrintNode(node);
                    Fatal("%s: undefined glide plane for (%d,%d) segment %d",
                          msg, node->myTag.domainID, node->myTag.index, j);
                }
            }
        }

        return;
}


/**************************************************************************
 *
 *      Function:     CheckPlanes_RhomboVa
 *      Description:  Examine the glide plane for each local segment and
 *                    compare the plane to the known valid planes for
 *                    rhombohedral Vanadium.  If an unexpected glide plane
 *                    is found, an error message is printed.
 *
 *      Arguments:
 *          msg  text to be displayed with the error message when
 *               an unexpected glide plane is found
 *
 *************************************************************************/
void CheckPlanes_RhomboVa(Home_t *home, const char *msg)
{
        int    i, j;
        Node_t *node;
        static real8 glideRef[47][3] = {
               /* for burgRef[0][*] */
           { -7.0710678119e-01,  7.0710678119e-01,  0.0000000000e+00},
           {  7.0710678119e-01,  0.0000000000e+00, -7.0710678119e-01},
           {  0.0000000000e+00, -7.0710678119e-01,  7.0710678119e-01},
               /* for burgRef[1][*] */
           { -7.0710678119e-01,  7.0710678119e-01,  0.0000000000e+00},
           {  7.0348253195e-02, -7.0535491891e-01, -7.0535491891e-01},
           {  7.0535491891e-01, -7.0348253195e-02,  7.0535491891e-01},
               /* for burgRef[2][*] */
           {  7.0710678119e-01,  0.0000000000e+00, -7.0710678119e-01},
           {  7.0348253195e-02, -7.0535491891e-01, -7.0535491891e-01},
           { -7.0535491891e-01, -7.0535491891e-01,  7.0348253195e-02},
               /* for burgRef[3][*] */
           {  0.0000000000e+00, -7.0710678119e-01,  7.0710678119e-01},
           {  7.0535491891e-01, -7.0348253195e-02,  7.0535491891e-01},
           { -7.0535491891e-01, -7.0535491891e-01,  7.0348253195e-02},
               /* for burgRef[4][*] */
           {  7.0710678119e-01, -7.0710678119e-01,  0.0000000000e+00},
           {  9.9775148521e-01, -4.7391844069e-02, -4.7391844069e-02},
           { -4.7391844069e-02,  9.9775148521e-01, -4.7391844069e-02},
           {  7.0535491891e-01,  7.0535491891e-01, -7.0348253195e-02},
               /* for burgRef[5][*] */
           { -7.0710678119e-01,  0.0000000000e+00,  7.0710678119e-01},
           {  7.0535491891e-01, -7.0348253195e-02,  7.0535491891e-01},
           {  9.9775148521e-01, -4.7391844069e-02, -4.7391844069e-02},
           {  4.7391844069e-02,  4.7391844069e-02, -9.9775148521e-01},
               /* for burgRef[6][*] */
           {  0.0000000000e+00,  7.0710678119e-01, -7.0710678119e-01},
           {  7.0348253195e-02, -7.0535491891e-01, -7.0535491891e-01},
           {  4.7391844069e-02, -9.9775148521e-01,  4.7391844069e-02},
           {  4.7391844069e-02,  4.7391844069e-02, -9.9775148521e-01},
               /* for burgRef[7][*] */
           { -8.1649658093e-01,  4.0824829046e-01,  4.0824829046e-01},
           { -8.4357735225e-01, -3.7972177365e-01, -3.7972177365e-01},
           {  7.0348253195e-02, -7.0535491891e-01, -7.0535491891e-01},
           {  9.9775148521e-01, -4.7391844069e-02, -4.7391844069e-02},
               /* for burgRef[8][*] */
           {  4.0824829046e-01, -8.1649658093e-01,  4.0824829046e-01},
           { -4.7391844069e-02,  9.9775148521e-01, -4.7391844069e-02},
           { -7.0535491891e-01,  7.0348253195e-02, -7.0535491891e-01},
               /* for burgRef[9][*] */
           {  4.0824829046e-01,  4.0824829046e-01, -8.1649658093e-01},
           {  4.7391844069e-02,  4.7391844069e-02, -9.9775148521e-01},
           { -3.7972177365e-01, -3.7972177365e-01, -8.4357735225e-01},
           { -7.0535491891e-01, -7.0535491891e-01,  7.0348253195e-02},
               /* for burgRef[10][*] */
           {  7.0710678119e-01, -7.0710678119e-01,  0.0000000000e+00},
           { -4.5837351184e-01,  3.8214699675e-01,  8.0240725104e-01},
           {  4.7391844069e-02,  4.7391844069e-02, -9.9775148521e-01},
           {  3.8214699675e-01, -4.5837351184e-01,  8.0240725104e-01},
               /* for burgRef[11][*] */
           {  7.0710678119e-01,  0.0000000000e+00, -7.0710678119e-01},
           {  4.5837351184e-01, -8.0240725104e-01, -3.8214699675e-01},
           {  4.7391844069e-02, -9.9775148521e-01,  4.7391844069e-02},
           {  3.8214699675e-01,  8.0240725104e-01, -4.5837351184e-01},
               /* for burgRef[12][*] */
           {  0.0000000000e+00,  7.0710678119e-01, -7.0710678119e-01},
           { -8.0240725104e-01,  4.5837351184e-01, -3.8214699675e-01},
           {  9.9775148521e-01, -4.7391844069e-02, -4.7391844069e-02},
           { -8.0240725104e-01, -3.8214699675e-01,  4.5837351184e-01}
        };


        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            for (j = 0; j < node->numNbrs; j++) {
                int k, match = 0;

                for (k = 0; k < 47; k++) {
                    real8 diff[3];

                    diff[0] = node->nx[j]-glideRef[k][0];
                    diff[1] = node->ny[j]-glideRef[k][1];
                    diff[2] = node->nz[j]-glideRef[k][2];

                    if ((fabs(diff[0]) < 1.0e-4) &&
                        (fabs(diff[1]) < 1.0e-4) &&
                        (fabs(diff[2]) < 1.0e-4)) {
                        match = 1;
                        break;
                    }

                    diff[0] = node->nx[j]+glideRef[k][0];
                    diff[1] = node->ny[j]+glideRef[k][1];
                    diff[2] = node->nz[j]+glideRef[k][2];

                    if ((fabs(diff[0]) < 1.0e-4) &&
                        (fabs(diff[1]) < 1.0e-4) &&
                        (fabs(diff[2]) < 1.0e-4)) {
                        match = 1;
                        break;
                    }
                }

                if (!match) {
                    printf("%s: unknown glide plane!\n", msg);
                    PrintNode(node);
                }
            }
        }

        return;
}


/**************************************************************************
 *
 *      Function:     Print2x2
 *      Description:  Print the components of a 2 X 2 matrix
 *
 *      Arguments:
 *          msg  text to be displayed with the matrix data.
 *          A    2 X 2 matrix to be displayed
 *
 *************************************************************************/
void Print2x2(const char *msg, real8 A[2][2])
{
        printf("\n %s\n", msg);

        printf("%.15e %.15e\n"  , A[0][0], A[0][1]);
        printf("%.15e %.15e\n"  , A[1][0], A[1][1]);

        return;
}


/**************************************************************************
 *
 *      Function:     Synchronize
 *      Description:  Explicitly synchronize parallel execution via an MPI
 *                    barrier and print a message when all tasks have reached
 *                    the barrier.  For debug only.
 *
 *      Arguments:
 *          msg  text to be displayed by task zero after all tasks have
 *               hit the barrier.
 *
 *************************************************************************/
void Synchronize(Home_t *home, const char *msg)
{

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        if (home->myDomain == 0) {
            printf(" *** %s: All tasks synchronized\n", msg);
            fflush(NULL);
        }
}


/**************************************************************************
 *
 *      Function:     VerifyConnections
 *      Description:  Loop through all local nodes, verify that the
 *                    neighboring nodes exist as local or ghost nodes in
 *                    this task, and then verify tha the neighboring node
 *                    is connected back to the local node.
 *
 *      Arguments:
 *          msg  text to be displayed along with any error message
 *
 *************************************************************************/
void VerifyConnections(Home_t *home, const char *msg)
{

/*
 *  Loop through all the native nodes
 */
    for (int i = 0; i < home->newNodeKeyPtr; i++) {
        Node_t *node;

        if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
            continue;
        }

/*
 *      Loop over every neighbor of the native node
 */
        for (int j = 0; j < node->numNbrs; j++) {
            int     foundLink;
            Node_t *nbr;

/*
 *          First verify that the neighboring node is found
 */
            nbr = GetNeighborNode(home, node, j);

            if (nbr == (Node_t *)NULL) {
                PrintNode(node);
                Fatal("VerifyConnections: %s Missing neighbor %d", msg, j);
            }

/*
 *          Now verify that the neighbor is also connected to the native
 *          node.
 */
            foundLink = 0;

            for (int k = 0; k < nbr->numNbrs; k++) {
                Node_t *backLink;

                backLink = GetNeighborNode(home, nbr, k);

                if (backLink == node) {
                    foundLink = 1;
                    break;
                }
            }

            if (foundLink == 0) {
                PrintNode(node);
                PrintNode(nbr);
                Fatal("VerifyConnections: %s inconsistent linkage", msg);
            }
        }
    }

    return;
}


/**************************************************************************
 *
 *      Function:     VerifyLocalConnections
 *      Description:  Loop through all native nodes, verify that any
 *                    local neighboring nodes exist on this task, and
 *                    then verify that the neighboring node
 *                    is connected back to the native node.
 *
 *      Arguments:
 *          msg  text to be displayed along with any error message
 *
 *************************************************************************/
void VerifyLocalConnections(Home_t *home, const char *msg)
{
    for (int i = 0; i < home->newNodeKeyPtr; i++) {
        Node_t *node;

        if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
            continue;
        }

        for (int j = 0; j < node->numNbrs; j++) {
            int     foundLink;
            Node_t *nbr;

/*
 *          Only look at neighboring nodes that are in the same domain
 */
            if (node->nbrTag[j].domainID == home->myDomain) {
                nbr = GetNeighborNode(home, node, j);

/*
 *              First verify that the local neighboring node is found
 */
                if (nbr == (Node_t *)NULL) {
                    PrintNode(node);
                    Fatal("VerifyLocalConnections: %s Missing neighbor %d",
                          msg, j);
                }

                foundLink = 0;

/*
 *              Now make sure the neighbor node thinks it is linked back
 */
                for (int k = 0; k < nbr->numNbrs; k++) {
                    Node_t *backLink;

                    backLink = GetNeighborNode(home, nbr, k);

                    if (backLink == node) {
                        foundLink = 1;
                        break;
                    }
                }

                if (foundLink == 0) {
                    PrintNode(node);
                    PrintNode(nbr);
                    Fatal("VerifyLocalConnections: %s inconsistent linkage",
                          msg);
                }

            }  // end if nbr is local

        }
    }

    return;
}

/**************************************************************************
 *
 *      Function:     CheckSegmentLengths
 *      Description:  Will construct a complete sgment list and check
 *                    for any segments that may be unusually large or small.
 *                    This was created to diagnose secondary ghost errors.
 *
 *************************************************************************/

static real8 Min (const real8 x, const real8 y, const real8 z) { return ( (x<y) ? ( (x<z) ? x : z ) : ( (y<z) ? y : z ) ); }

void CheckSegmentLengths 
(
   const Home_t *home,    ///< points to Home
   const char   *msg      ///< additional message text to be included
)
{
   Node_t    **native = home->nodeKeys      ;   // native = sparse array of native node pointers
   int         ncnt   = home->newNodeKeyPtr ;   // ncnt   = number of native node pointers
   Node_t    **ghosts = home->ghostNodeList ;   // ghosts = sparse array of ghost node pointers
   int         gcnt   = home->ghostNodeCount;   // gcnt   = number of ghost node pointers

   // create a list of segments for this domain...

   int         scnt = 0;
   Segment_t  *segs = Segment_List_Create (scnt, (const Node_t **) native, ncnt,
                                                 (const Node_t **) ghosts, gcnt );

   if (segs)
   {
      Param_t  *param  = ( home ? home->param : 0 );

      // obtain the overall simulation dimensions to compute the smallest cell dimensions
      // and apply periodic repositions...

      real8 lx = param->maxSideX - param->minSideX;   // lx = simulation size (x)
      real8 ly = param->maxSideY - param->minSideY;   // ly = simulation size (y)
      real8 lz = param->maxSideZ - param->minSideZ;   // lz = simulation size (z)

      int   nx = param->nXcells;                      // nx = number of cells (x)
      int   ny = param->nYcells;                      // ny = number of cells (y)
      int   nz = param->nZcells;                      // nz = number of cells (z)

      real8 cx = ( (nx>0) ? (lx/nx) : 0.0 );          // cx = cell dimension (x)
      real8 cy = ( (ny>0) ? (ly/ny) : 0.0 );          // cy = cell dimension (y)
      real8 cz = ( (nz>0) ? (lz/nz) : 0.0 );          // cz = cell dimension (z)

      real8 cmin = Min(cx,cy,cz);                     // cmin is shortest cell dimension

      real8 rann = param->rann;                       // rann = simulation annihilation distance

      real8 smin = rann/2.0;                          // smin = shortest acceptible segment
      real8 smax = cmin;                              // smax = longest  acceptible segment

      int nfail = Segment_List_Check(segs,scnt, lx,ly,lz, smin,smax, msg);

      if (nfail>0) Fatal("Segment length check fail (dom=%d)\n", home->myDomain );
   }

   if (segs) { delete [] segs; segs=0; }
}
