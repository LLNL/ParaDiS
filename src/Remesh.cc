/*****************************************************************************
 *
 *      Module:         Remesh.c
 *      Description:    This module contains functions common to more
 *                      than 1 of the supported version of remesh, plus
 *                      a generic entry function that invokes the
 *                      proper remesh version.
 *
 *      Included functions:
 *              CutSurfaceSegments()
 *              EstCoarsenForces()
 *              EstRefinementForces()
 *              MeshCoarsen()
 *              ProceedWithCoarsen()
 *              Remesh()
 *
 *****************************************************************************/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"
#include "Topology.h"
#include "Remesh.h"

static int dbgDom;

/*-------------------------------------------------------------------------
 *
 *      Function:    ProceedWithCoarsen
 *      Description: When removing nodes that have been flagged to be
 *                   coarsened out, there is the possibility that
 *                   multiple threads could simultaneously attempt
 *                   to perform actions on overlapping sets of nodes.
 *                   To prevent that, this function will look through
 *                   the list of 'committed' coarsen events and determine
 *                   if the specified nodes have been involved in any
 *                   previous events.  If not, a coarsen event involving
 *                   the specified nodes is permitted.
 *
 *                   WARNING: The thread calling this function is responsible
 *                            for insuring that no other thread will be
 *                            accessing the 'coarsen list' while this
 *                            function is executing.
 *
 *      Arguments:
 *          list      Pointer to list of potential or committed coarsen events
 *          listEnts  number of elements on <list>
 *          nodeTag   Tag of the node the caller wishes to coarsen out
 *          nbr1Tag   Tag of the 1st neighbor of node with tag <nodeTag>
 *          nbr2Tag   Tag of the 2nd neighbor of node with tag <nodeTag>
 *
 *      Returns:  0 if the node may not be coarsened out
 *                1 if the node may be coarsened out.
 *
 *------------------------------------------------------------------------*/
int ProceedWithCoarsen(CoarsenList_t *list, int listEnts, Tag_t *nodeTag,
                       Tag_t *nbr1Tag, Tag_t *nbr2Tag)
{
        int i;

        for (i = 0; i < listEnts; i++) {

            if (list[i].status != REMESH_COMMITTED) {
                continue;
            }

            if (((nodeTag->domainID == list[i].nodeTag.domainID) &&
                 (nodeTag->index    == list[i].nodeTag.index   )) ||
                ((nodeTag->domainID == list[i].nbr1Tag.domainID) &&
                 (nodeTag->index    == list[i].nbr1Tag.index   )) ||
                ((nodeTag->domainID == list[i].nbr2Tag.domainID) &&
                 (nodeTag->index    == list[i].nbr2Tag.index   ))) {
                return(0);
            }

            if (((nbr1Tag->domainID == list[i].nodeTag.domainID) &&
                 (nbr1Tag->index    == list[i].nodeTag.index   )) ||
                ((nbr1Tag->domainID == list[i].nbr1Tag.domainID) &&
                 (nbr1Tag->index    == list[i].nbr1Tag.index   )) ||
                ((nbr1Tag->domainID == list[i].nbr2Tag.domainID) &&
                 (nbr1Tag->index    == list[i].nbr2Tag.index   ))) {
                return(0);
            }

            if (((nbr2Tag->domainID == list[i].nodeTag.domainID) &&
                 (nbr2Tag->index    == list[i].nodeTag.index   )) ||
                ((nbr2Tag->domainID == list[i].nbr1Tag.domainID) &&
                 (nbr2Tag->index    == list[i].nbr1Tag.index   )) ||
                ((nbr2Tag->domainID == list[i].nbr2Tag.domainID) &&
                 (nbr2Tag->index    == list[i].nbr2Tag.index   ))) {
                return(0);
            }

        }

        return(1);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    EstRefineMentForces
 *      Description: Estimate the forces on the resultant segments
 *                   when an existing segment is split during
 *                   mesh refinement.
 *
 *      Arguments:
 *          node1    pointer to first endpoint of the original segment
 *          node2    pointer to second endpoint of the original segment
 *          newPos   coordinates of the point at which the node1/node2
 *                   segment will be split
 *          vec      vector from <node1> to <node2>
 *          f0Seg1   3 element array in which to return estimated segment
 *                   force on <node1> from a segment from <node1> to the
 *                   location <newPos>
 *          f1Seg1   3 element array in which to return estimated segment
 *                   force at location <newPos> from a segment from <node1>
 *                   to the location <newPos>
 *          f0Seg2   3 element array in which to return estimated segment
 *                   force at location <newPos> from a segment from <node2>
 *                   to the location <newPos>
 *          f1Seg2   3 element array in which to return estimated segment
 *                   force on <node2> from a segment from <node2> to the
 *                   location <newPos>
 *
 *------------------------------------------------------------------------*/
void EstRefinementForces(Home_t *home, Node_t *node1, Node_t *node2,
                         real8 newPos[3], real8 vec[3],
                         real8 f0Seg1[3], real8 f1Seg1[3],
                         real8 f0Seg2[3], real8 f1Seg2[3])
{
        int     arm12, arm21;
        real8   p1[3], p2[3], oldfp1[3], oldfp2[3], burg[3];

        arm12 = GetArmID(node1, node2);
        arm21 = GetArmID(node2, node1);

        burg[X] = node1->burgX[arm12];
        burg[Y] = node1->burgY[arm12];
        burg[Z] = node1->burgZ[arm12];

        oldfp1[X] = node1->armfx[arm12];
        oldfp1[Y] = node1->armfy[arm12];
        oldfp1[Z] = node1->armfz[arm12];

        oldfp2[X] = node2->armfx[arm21];
        oldfp2[Y] = node2->armfy[arm21];
        oldfp2[Z] = node2->armfz[arm21];

        p1[X] = node1->x;
        p1[Y] = node1->y;
        p1[Z] = node1->z;

        p2[X] = p1[X] + vec[X];
        p2[Y] = p1[Y] + vec[Y];
        p2[Z] = p1[Z] + vec[Z];

        FindSubFSeg(home, p1, p2, burg, oldfp1, oldfp2, newPos, f0Seg1,
                    f1Seg1, f0Seg2, f1Seg2);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    EstCoarsenForces
 *      Description: Estimate the forces on the resultant segment
 *                   when an existing discretization node is removed
 *                   during mesh coarsening.
 *
 *      Arguments:
 *          node1    pointer to first neighbor of node to be coarsened
 *          node2    pointer to node to be coarsened out
 *          node3    pointer to second neighbor of node to be coarsened
 *          f0Seg1   3 element array in which to return estimated segment
 *                   force on <node1> from a segment from <node1> to <node2>
 *          f1Seg1   3 element array in which to return estimated segment
 *                   force on <node2> from a segment from <node1> to <node2>
 *
 *------------------------------------------------------------------------*/
void EstCoarsenForces(Home_t *home, Node_t *node1, Node_t *node2,
                      Node_t *node3, real8 f0Seg[3], real8 f1Seg[3])
{
        int     i, arm21, arm23, arm12, arm32;
        real8   burg1[3], burg2[3];
        real8   p1[3], p2[3], p3[3];
        real8   fp0[3], fp1[3], fp2[3], fp3[3];
        real8   len1[3], len2[3];
        Param_t *param;

        param = home->param;

        arm21 = GetArmID(node2, node1);
        arm23 = GetArmID(node2, node3);
        arm12 = GetArmID(node1, node2);
        arm32 = GetArmID(node3, node2);

        p1[X] = node1->x; p2[X] = node2->x; p3[X] = node3->x;
        p1[Y] = node1->y; p2[Y] = node2->y; p3[Y] = node3->y;
        p1[Z] = node1->z; p2[Z] = node2->z; p3[Z] = node3->z;

/*
 *      If periodic boundaries are enabled, the nodes may
 *      be on opposite side of the problem space, so adjust
 *      the positions accordingly.
 */
        PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p2[X], &p2[Y], &p2[Z]);
        PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p3[X], &p3[Y], &p3[Z]);

        burg1[X] = node1->burgX[arm12];
        burg1[Y] = node1->burgY[arm12];
        burg1[Z] = node1->burgZ[arm12];

        burg2[X] = node2->burgX[arm23];
        burg2[Y] = node2->burgY[arm23];
        burg2[Z] = node2->burgZ[arm23];

        fp0[X] = node1->armfx[arm12];
        fp0[Y] = node1->armfy[arm12];
        fp0[Z] = node1->armfz[arm12];

        fp1[X] = node2->armfx[arm21];
        fp1[Y] = node2->armfy[arm21];
        fp1[Z] = node2->armfz[arm21];

        fp2[X] = node2->armfx[arm23];
        fp2[Y] = node2->armfy[arm23];
        fp2[Z] = node2->armfz[arm23];

        fp3[X] = node3->armfx[arm32];
        fp3[Y] = node3->armfy[arm32];
        fp3[Z] = node3->armfz[arm32];

/*
 *      It is possible to we'll have to deal with extremely short
 *      segments here (zero length segments created during collisions and
 *      short segments being coarsened out).  Since the FindFSegComb()
 *      function does not handle these cases well, we do special handling
 *      of them here.
 *
 *      If we find a zero-length segment, the force on a segment between
 *      nodes 1 and 3 would equal the force on the non-zero length segment
 *      attached to node 2.  Otherwise, use FindFSegComb() to get the force.
 */
        for (i = 0; i < 3; i++) {
            len1[i] = p1[i] - p2[i];
            len2[i] = p2[i] - p3[i];
        }

/*
 *      Segments less than .01b in length are treated as zero length
 */
        if ((len1[X]*len1[X] + len1[Y]*len1[Y] + len1[Z]*len1[Z]) < 1.0e-4) {
            f0Seg[X] = node2->armfx[arm23];
            f0Seg[Y] = node2->armfy[arm23];
            f0Seg[Z] = node2->armfz[arm23];
            f1Seg[X] = node3->armfx[arm32];
            f1Seg[Y] = node3->armfy[arm32];
            f1Seg[Z] = node3->armfz[arm32];
        } else if ((len2[X]*len2[X]+len2[Y]*len2[Y]+len2[Z]*len2[Z])<1.0e-4) {
            f0Seg[X] = node1->armfx[arm12];
            f0Seg[Y] = node1->armfy[arm12];
            f0Seg[Z] = node1->armfz[arm12];
            f1Seg[X] = node2->armfx[arm21];
            f1Seg[Y] = node2->armfy[arm21];
            f1Seg[Z] = node2->armfz[arm21];
        } else {
            FindFSegComb(home, p1, p2, p3, burg1, burg2,
                         fp0, fp1, fp2, fp3, f0Seg, f1Seg);
        }

        return;
}


#ifdef _ARLFEM
/*-------------------------------------------------------------------------
 *
 *      Function:       CutSurfaceSegments
 *      Description:    Loop through all the local nodes looking for 
 *                      surface nodes with connections to other nodes
 *                      on the same surface.
 *
 *-----------------------------------------------------------------------*/
static void CutSurfaceSegments(Home_t *home)
{
        int    i, j, thisDomain, globalOp = 1;
        real8  eps, tmp3[3];
        Node_t *node, *nbr;


        eps = 1.0e-12;
        thisDomain = home->myDomain;

        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

/*
 *          If the node is not a surface node, skip it.
 */
            if (HAS_NONE_OF_CONSTRAINTS(node->constraint, SURFACE_NODE)) {
                continue;
            }

            j = 0;

            while (j < node->numNbrs) {

                nbr = GetNeighborNode(home, node, j);
                if (nbr == (Node_t *)NULL) {
                    j++;
                    continue;
                }

/*
 *              if the neighboring node is not on the surface, skip
 *              to the next neighbor
 */
                if (HAS_NONE_OF_CONSTRAINTS(nbr->constraint, SURFACE_NODE)) {
                    j++;
                    continue;
                }
/*
 *              If the surface-normal vectors for the two nodes do not
 *              match the nodes are on different surfaces so we can skip
 *              to the next neighbor.
 */
                cross(node->surfaceNorm, nbr->surfaceNorm, tmp3);

                if ((tmp3[0]*tmp3[0] +
                     tmp3[1]*tmp3[1] +
                     tmp3[2]*tmp3[2]) > eps) {
                    j++;
                    continue;
                }
/*
 *              We found a segment we need to cut.  If the neighbor is
 *              in a remote domain, cut the connection in the lower
 *              numbered domain only.
 */
                if (thisDomain > nbr->myTag.domainID) {
                    j++;
                    continue;
                } 

                ChangeArmBurg(home, node, &node->nbrTag[j], 0, 0, 0,
                              0, 0, 0, globalOp, DEL_SEG_NONE);
                ChangeArmBurg(home, nbr, &node->myTag, 0, 0, 0,
                              0, 0, 0, globalOp, DEL_SEG_NONE);

                if ((nbr->myTag.domainID == thisDomain) &&
                    (nbr->numNbrs == 0)) {
                    RemoveNode(home, nbr, globalOp);
                }
            }
        }

        return;
}
#endif  /* ifdef _ARLFEM */


/*-------------------------------------------------------------------------
 *
 *      Function:    MeshCoarsen
 *      Description: 
 *
 *------------------------------------------------------------------------*/
void MeshCoarsen(Home_t *home)
{
#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom = -1;
#endif

        int      thisDomain = home->myDomain;
        Param_t *param      = home->param;

        int      allowFuzzyPlanes = 0;

        real8 cellLength = home->param->Lx / home->param->nXcells;
              cellLength = MIN(cellLength, home->param->Ly / home->param->nYcells);
              cellLength = MIN(cellLength, home->param->Lz / home->param->nZcells);

        real8 cutoffLength1 = param->maxSeg;
        real8 cutoffLength2 = MIN(cutoffLength1, 0.45 * cellLength);

        real8 areaMin    = param->remeshAreaMin;
        real8 areaMin2   = areaMin * areaMin;
        real8 delta      = 1.0e-16;

        int   localCoarsenCnt=0;

/*
 *      Allocate a block of structures to hold info about rediscretization
 *      events.  We allocate enough structures to refine coarsen twice as
 *      many segments as there are nodes which is overkill, but it means we
 *      shouldn't run out of entries and have to reallocate the block,
 *      plus minimizes the number of operations we'll have to do within
 *      the thread-critical sections of code.
 */
        int            coarsenListSize = 2 * home->newNodeKeyPtr;
        int            coarsenListEnts = 0;
        CoarsenList_t *coarsenList     = (CoarsenList_t *)NULL;

        if (coarsenListSize > 0) {
            coarsenList = (CoarsenList_t *)malloc(sizeof(CoarsenList_t) *
                                                  coarsenListSize);
        }

/*
 *      First loop through the nodes making a list of candidates for
 *      coarsening.  This part part is easily multi-threaded.  The
 *      actual node removal requires more care to be thread-safe and
 *      will be done later in a separate loop.
 */
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int    nIndex;
            int    numNodes;
            int    threadID, threadIterStart, threadIterEnd;        

            numNodes = home->newNodeKeyPtr;

            GetThreadIterationIndices(numNodes, &threadID, &threadIterStart,
                                      &threadIterEnd);

            for (nIndex = threadIterStart; nIndex < threadIterEnd; nIndex++) {
                int     hasRemoteNbr;
                int     domOwnsSeg1, domOwnsSeg2;
                int     nbr1IsRemote, nbr2IsRemote;
                real8   r1, r2, r3;
                real8   s, area2;
                real8   dvec1xdt, dvec1ydt, dvec1zdt;
                real8   dvec2xdt, dvec2ydt, dvec2zdt;
                real8   dvec3xdt, dvec3ydt, dvec3zdt;
                real8   dr1dt, dr2dt, dr3dt, dsdt, darea2dt;
                real8   vec1x, vec1y, vec1z;
                real8   vec2x, vec2y, vec2z;
                real8   vec3x, vec3y, vec3z;
                real8   gp0[3], gp1[3];
                Node_t *node, *nbr1, *nbr2;

                if ((node = home->nodeKeys[nIndex]) == (Node_t *)NULL) {
                    continue;
                }

/*
 *              Check for various conditions that will exempt a node from
 *              removal:
 *
 *              1) does not have exactly 2 arms
 *              2) node is a 'pinned' node
 *              3) node is flagged as exempt from coarsen operations
 *              4) current domain does not 'own' at least one of the  segments
 *              5) If the node's arms are on different glide planes we might
 *                 not allow the node to be removed.
 *
 *              NOTE: For now, a node that is pinned in ANY dimension is
 *                    treated here as if it is pinned in ALL dimensions.
 *                    This may change in the future
 */
                if ((node->numNbrs != 2)                                  ||
                    HAS_ANY_OF_CONSTRAINTS(node->constraint, PINNED_NODE) ||
                    (node->flags & NO_MESH_COARSEN)) {
                    continue;
                }

                nbr1 = GetNeighborNode(home, node, 0);
                nbr2 = GetNeighborNode(home, node, 1);

                if ((nbr1 == (Node_t *)NULL) || (nbr2 == (Node_t *)NULL)) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                domOwnsSeg1 = DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nbr1->myTag);
                domOwnsSeg2 = DomainOwnsSeg(home, OPCLASS_REMESH, thisDomain, &nbr2->myTag);

                if ((!domOwnsSeg1) && (!domOwnsSeg2)) {
                    continue;
                }

                nbr1IsRemote = (node->myTag.domainID != nbr1->myTag.domainID);
                nbr2IsRemote = (node->myTag.domainID != nbr2->myTag.domainID);

                hasRemoteNbr = nbr1IsRemote | nbr2IsRemote;

/*
 *              Calculate the lengths of the node's 2 arms plus
 *              the distance between the two neighbor nodes.
 *
 *              If periodic boundaries are enabled, the nodes may
 *              be on opposite side of the problem space, so adjust
 *              the lengths/distances accordingly.
 */
                vec1x = nbr1->x - node->x;
                vec1y = nbr1->y - node->y;
                vec1z = nbr1->z - node->z;

                vec2x = nbr2->x - node->x;
                vec2y = nbr2->y - node->y;
                vec2z = nbr2->z - node->z;

                vec3x = vec2x - vec1x;
                vec3y = vec2y - vec1y;
                vec3z = vec2z - vec1z;

                ZImage(param, &vec1x, &vec1y, &vec1z);
                ZImage(param, &vec2x, &vec2y, &vec2z);
                ZImage(param, &vec3x, &vec3y, &vec3z);

                r1 = sqrt(vec1x*vec1x + vec1y*vec1y + vec1z*vec1z);
                r2 = sqrt(vec2x*vec2x + vec2y*vec2y + vec2z*vec2z);
                r3 = sqrt(vec3x*vec3x + vec3y*vec3y + vec3z*vec3z);

/*
 *              If we are enforcing the use of segment glide planes, we do
 *              not want to coarsen out a node whose arms are on different
 *              glide planes.
 *
 *              However, when glide planes constraints are being enforced, we
 *              tend to get a lot of debris (i.e. tiny triangles or quadrangles,
 *              really short segments, etc) that oscillate quickly and adversely
 *              affect the timestep.  So, under the following circumstances,
 *              we'll allow the code to violate glide plane constraints.
 *
 *              1) If the node has an attached segment less than 1b in length
 *              2) If the node has two segments less than 20% of the minimum
 *                 segment length AND both of the node's neighbors are connected
 *                 to each other (i.e. 3 nodes form a triangle)
 *              3) If the glide planes are allowed to be 'fuzzy' add
 *                 a few extra exceptions (see below)
 */
                if (param->enforceGlidePlanes) {
                    int   connectionID, violateGlidePlanesOK = 0;

                    if ((r1 < 1.0) ||
                        (r2 < 1.0) ||
                        ((MAX(r1, r2) < MAX(1.0, 0.20 * param->minSeg) &&
                         (Connected(nbr1, nbr2, &connectionID))))) {
                        violateGlidePlanesOK = 1;
                    }

                    gp0[0] = node->nx[0];
                    gp0[1] = node->ny[0];
                    gp0[2] = node->nz[0];

                    gp1[0] = node->nx[1];
                    gp1[1] = node->ny[1];
                    gp1[2] = node->nz[1];

/*
 *                  Glide planes are being used, but if we are allowing
 *                  some 'fuzziness' in the planes there are a couple
 *                  situations in which we allow glide plane constraints
 *                  to be violated.
 *
 *                    1)  the two segment glide planes are within a small
 *                        number of degrees
 *                    2)  If either segment is shorter than the annihilation
 *                        distance
 *                    3)  when each segment is mapped to the closest
 *                        precise glide plane, if the two precise planes are
 *                        the same
 */
                    if (param->allowFuzzyGlidePlanes) {

                        if (fabs(DotProduct(gp0, gp1)) > 0.9555) {
                            violateGlidePlanesOK = 1;
                        } else if ((r1 < param->rann) || (r2 <  param->rann)) {
                            violateGlidePlanesOK = 1;
                        } else {
                            real8 burg1[3], burg2[3];
                            real8 lineDir1[3], lineDir2[3];
                            real8 testPlane1[3], testPlane2[3];

                            burg1[X] = node->burgX[0];
                            burg1[Y] = node->burgY[0];
                            burg1[Z] = node->burgZ[0];

                            burg2[X] = node->burgX[1];
                            burg2[Y] = node->burgY[1];
                            burg2[Z] = node->burgZ[1];

                            lineDir1[X] = vec1x;
                            lineDir1[Y] = vec1y;
                            lineDir1[Z] = vec1z;

                            lineDir2[X] = vec2x;
                            lineDir2[Y] = vec2y;
                            lineDir2[Z] = vec2z;

/*
 *                          FindPreciseGlidePlanes() just uses l cross b if
 *                          fuzzy planes are allowed, in this case we've
 *                          set 'allowFuzzyPlanes' to false to explicitly
 *                          disallow the use of 'fuzzy' planes so we truly
 *                          find the closest precise plane.
 */
                            FindPreciseGlidePlane(home, burg1, lineDir1,
                                                  testPlane1, allowFuzzyPlanes);
                            FindPreciseGlidePlane(home, burg2, lineDir2,
                                                  testPlane2, allowFuzzyPlanes);

                            if (fabs(DotProduct(testPlane1,testPlane2))>0.99) {
                                violateGlidePlanesOK = 1;
                            }
                        }

                    }  /* end if (param->allowFuzzyGlidePlanes) */

                    if (!violateGlidePlanesOK) {
                        real8 tmp3[3];

                        cross(gp0, gp1, tmp3);

                        if (fabs(DotProduct(tmp3, tmp3)) > 1.0e-3) {
                            continue;
                        }
                    }

                }  /* end if (param->enforceGlidePlanes) */

/*
 *              If coarsening out a node would leave a segment longer
 *              than a defined length, the node should not be removed.
 *              This 'cutoff length' is the maximum segment length
 *              if all involved nodes are within the same domain, but
 *              if a remote node is involved, we need to set the
 *              cutoff length to (at most) 1/2 the cell length.  This
 *              is needed because the remote node may potentially be involved
 *              in a simultaneous mesh coarsening in the remote domain,
 *              and although the node would not be removed, it could
 *              be repositioned resulting in a segment that spanned
 *              more than 2 cells... this is a bad thing.
 */
                if (((hasRemoteNbr == 0) && (r3 > cutoffLength1)) ||
                    ((hasRemoteNbr == 1) && (r3 > cutoffLength2))) {
                    continue;
                }

/*
 *              Check if the area of the triangle defined by node
 *              and its two neighbors, plus determine if that area
 *              is increasing or decreasing.
 */

                s = 0.5 * (r1 + r2 + r3);
                area2 = (s * (s-r1) * (s-r2) * (s-r3));

                dvec1xdt = nbr1->vX - node->vX;
                dvec1ydt = nbr1->vY - node->vY;
                dvec1zdt = nbr1->vZ - node->vZ;

                dvec2xdt = nbr2->vX - node->vX;
                dvec2ydt = nbr2->vY - node->vY;
                dvec2zdt = nbr2->vZ - node->vZ;

                dvec3xdt = dvec2xdt - dvec1xdt;
                dvec3ydt = dvec2ydt - dvec1ydt;
                dvec3zdt = dvec2zdt - dvec1zdt;

                dr1dt = ((vec1x * dvec1xdt) + (vec1y * dvec1ydt) +
                         (vec1z * dvec1zdt)) / (r1 + delta);

                dr2dt = ((vec2x * dvec2xdt) + (vec2y * dvec2ydt) +
                         (vec2z * dvec2zdt)) / (r2 + delta);
        
                dr3dt = ((vec3x * dvec3xdt) + (vec3y * dvec3ydt) +
                         (vec3z * dvec3zdt)) / (r3 + delta);

                dsdt = 0.5 * (dr1dt + dr2dt + dr3dt);

                darea2dt = (dsdt * (s-r1) * (s-r2) * (s-r3));
                darea2dt += s * (dsdt-dr1dt) * (s-r2) * (s-r3);
                darea2dt += s * (s-r1) * (dsdt-dr2dt) * (s-r3);
                darea2dt += s * (s-r1) * (s-r2) * (dsdt-dr3dt);

/*
 *              If the area is less than the specified minimum and shrinking,
 *              or one of the arms is less than the minimum segment length, the
 *              node should be removed, so add the node tag (and its neighbor
 *              tags to the 'coarsen' list.
 */
                if (((area2 < areaMin2) && (darea2dt < 0.0)) ||
                    ((r1 < param->minSeg) || (r2 < param->minSeg))) {
/*
 *                  If the node to be removed has a segment owned by a remote
 *                  node, we can't coarsen it out.  Reason is:
 *
 *                  Given nodes A, B, C, D connected as A--B--C--D where
 *                  A and B are in domain 1 and C and D are in domain 2 with
 *                  segment B--C owned by domain 1.  Domain 1 owns B and
 *                  segment B--C so it removes node B.  Domain 2 cannot
 *                  remove node C because it is connected to B via a segment
 *                  owned by domain 1.  However, MergeNode() would permit
 *                  nodes C and D to be merged by getting rid of D and
 *                  moving C to D's original location.  This breaks no
 *                  ownership rules, but the resulting single segment spans
 *                  the original positions A to D.  If the original segments
 *                  were of sufficient length, the resulting segment would
 *                  have a length exceeding <maxSeg> which can cause
 *                  problems...
 */
                    if (((domOwnsSeg1 == 0) && (nbr1IsRemote)) ||
                        ((domOwnsSeg2 == 0) && (nbr2IsRemote))) {
                        continue;
                    }

#ifdef _OPENMP
#pragma omp critical (CRIT_COARSEN_A)
#endif
                    {
/*
 *                      Just a quick sanity check, but make sure there is
 *                      still space on the list...
 */
                        if (coarsenListEnts < coarsenListSize) {

                            coarsenList[coarsenListEnts].nodeTag = node->myTag;
                            coarsenList[coarsenListEnts].nbr1Tag = nbr1->myTag;
                            coarsenList[coarsenListEnts].nbr2Tag = nbr2->myTag;
                            coarsenList[coarsenListEnts].status  = REMESH_PENDING;

                            coarsenListEnts += 1;
                        }

                    }  /* end "omp critical" section */
                }

            }  /* end for (i = threadIterStart; ...) */

        }  /* end "omp parallel" section */

/*
 *      Now loop through the list of nodes to be coarsened and remove
 *      them (if possible)
 */
#ifdef _OPENMP
#pragma omp parallel for if (coarsenListEnts > 1)
#endif
        for (int i=0; i < coarsenListEnts; i++) 
        {
            int    q, proceed;
            int    mergeDone, mergeStatus, globalOp;
            real8  newPos[3], f0seg1[3], f1seg1[3];
            Tag_t  nbr1Tag, nbr2Tag;
            Node_t *node, *nbr1, *nbr2, *mergedNode;
#ifdef DEBUG_TOPOLOGY_CHANGES
            Tag_t  oldTag1, oldTag2, oldTag3;
#endif

/*
 *          Verify that the neither the node nor its neighbors have
 *          been involved in a "committed" coarsen event in any other
 *          thread.  This is needed because coarsening a node modifies
 *          the node structures of both neighboring nodes plus various
 *          global data structures.  We can't risk different threads
 *          trying to access these data structures simultaneously,
 *          so we don't allow threads to work on overlapping sets
 *          of nodes.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_COARSEN_B)
#endif
            {
                proceed = ProceedWithCoarsen(coarsenList, coarsenListEnts,
                                             &coarsenList[i].nodeTag,
                                             &coarsenList[i].nbr1Tag,
                                             &coarsenList[i].nbr2Tag);
                if (proceed) {
                    coarsenList[i].status = REMESH_COMMITTED;
                } else {
                    coarsenList[i].status = REMESH_REJECTED;
                }

            }  /* end "omp critical" section */
                                       
/*
 *          Not allowed to jump out of a 'critical' section, so had to
 *          wait until we were out of the section before checking if we
 *          need to cycle to the next loop iteration.
 */
            if (!proceed) {
                continue;
            }

/*
 *          Get the node and its neighbors from the 'coarsen' list
 */
            node = GetNodeFromTag(home, coarsenList[i].nodeTag);
            nbr1 = GetNodeFromTag(home, coarsenList[i].nbr1Tag);
            nbr2 = GetNodeFromTag(home, coarsenList[i].nbr2Tag);

/*
 *          It is possible that another coarsen event has resulted
 *          in the removal of the current node or one of its neighbors
 *          (due to annihilation).  Check for this before attempting
 *          to remove the current node.
 *
 *          Given that the call to ProceedWithCoarsen() succeeded, this
 *          check *may* not be required.
 */
            if ((node == (Node_t *)NULL) ||
                (nbr1 == (Node_t *)NULL) ||
                (nbr2 == (Node_t *)NULL)) {
                continue;
            }


            EstCoarsenForces(home, nbr1, node, nbr2, f0seg1, f1seg1);

            mergeDone = 0;

            nbr1Tag = nbr1->myTag;
            nbr2Tag = nbr2->myTag;

/*
 *          If either of the neighbor nodes (or any of their neighbors)
 *          is in a remote domain, the operation must be treated as global.
 *          This is necessary to prevent an inconsistent linkage problem.
 *          For example, given the segments A--B--C--D where nodes A, B
 *          and C are in domain 1, D is in domain 2, and segment C--D is
 *          owned by domain 1:  domain 1 could coarsen node B into A, then
 *          node C into D.  If the first operation was not communicated
 *          to the remote domain, an inconsitency would arise.
 *
 *          NOTE:  It is safe to not distribute the purely local coarsen
 *                 operations so long as no other topological operations
 *                 are done after Remesh() but before the ghost node
 *                 are redistributed.
 */
            globalOp = ((nbr1->myTag.domainID != thisDomain) ||
                        (nbr2->myTag.domainID != thisDomain));

            for (q = 0; q < nbr1->numNbrs; q++) {
                globalOp |= (nbr1->nbrTag[q].domainID != thisDomain);
            }

            for (q = 0; q < nbr2->numNbrs; q++) {
                globalOp |= (nbr2->nbrTag[q].domainID != thisDomain);
            }

#ifdef DEBUG_TOPOLOGY_CHANGES
            oldTag1 = nbr1->myTag;
            oldTag2 = node->myTag;
            oldTag3 = nbr2->myTag;
#endif

/*
 *          The next block of code results in modification to 
 *          multiple node stuctures and global data structures.
 *          Access to those data structures is not yet thread-safe,
 *          so for now we need to insure that only 1 thread at a time
 *          will be modifying those structures.  Note: this will
 *          result in deletion of 1 or more nodes, but other checks
 *          above should guarantee that no other thread will be attempting
 *          to access the same nodes while this is being done.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_COARSEN_C)
#endif
            {
/*
 *              If the first neighbor is not exempt from a coarsen
 *              operation, attempt to merge the nodes.
 */
                if ((nbr1->flags & NO_MESH_COARSEN) == 0) {

                    newPos[X] = nbr1->x;
                    newPos[Y] = nbr1->y;
                    newPos[Z] = nbr1->z;

                    MergeNode(home, OPCLASS_REMESH, node, nbr1, newPos,
                              &mergedNode, &mergeStatus, globalOp);

                    mergeDone = mergeStatus & MERGE_SUCCESS;

                    if (mergeDone) {
                        localCoarsenCnt++;
                    }
                }

/*
 *              If the merge could not be done, try using
 *              the other neighbor.
 */
                if (mergeDone == 0) {
                    if ((nbr2->flags & NO_MESH_COARSEN) == 0) {
                        newPos[X] = nbr2->x;
                        newPos[Y] = nbr2->y;
                        newPos[Z] = nbr2->z;

                        MergeNode(home, OPCLASS_REMESH, node, nbr2, newPos,
                                  &mergedNode, &mergeStatus, globalOp);

                        mergeDone = mergeStatus & MERGE_SUCCESS;

                        if (mergeDone) {
                            localCoarsenCnt++;
                        }
                    }
                }

            }  /* end "omp critical" section */

/*
 *          If the merge was successful, update the forces
 *          on the remaining nodes.   Otherwise go back and
 *          continue looking for more nodes to coarsen out.
 */
            if (mergeDone == 0) continue;

#ifdef DEBUG_TOPOLOGY_CHANGES
            if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                printf("Coarsen: (%d,%d)--(%d,%d)--(%d,%d)\n",
                       oldTag1.domainID, oldTag1.index,
                       oldTag2.domainID, oldTag2.index,
                       oldTag3.domainID, oldTag3.index);
            }
#endif

/*
 *          It is possible coarsening out the node resulted in the
 *          annihilation of one or more of the neighboring nodes.
 *          if that's the case, skip any more operations on the 
 *          nodes involved.
 */
            nbr1 = GetNodeFromTag(home, nbr1Tag);
            nbr2 = GetNodeFromTag(home, nbr2Tag);

            if ((nbr1 == (Node_t *)NULL) && (nbr2 == (Node_t *)NULL)) {
                continue;
            }

/*
 *          The merge will have placed the resultant node at the
 *          location of either nbr1 or nbr2, but the merge function
 *          determines which of the two specified nodes is deleted
 *          and which survives, which means that nbr1 or nbr2 may
 *          have been the deleted and the <node> repositioned to
 *          the correct location... so if one of the nbr nodes does
 *          not exist anymore, <mergedNode> (if it exists) should
 *          be the node that replaced the nbr.
 */
            if (nbr1 == (Node_t *)NULL) {
                nbr1 = mergedNode;
            } else if (nbr2 == (Node_t *)NULL) {
                nbr2 = mergedNode;
            }

/*
 *          At this point, if we don't have a node at the location
 *          of at least one of the original nbr nodes, looks like
 *          some nodes were orphaned and deleted, so the force
 *          estimates we made are not applicable.
 */
            if ((nbr1 == (Node_t *)NULL) || (nbr2 == (Node_t *)NULL)) {
                continue;
            }

            mergedNode->flags |= NO_MESH_COARSEN;

/*
 *          Reset force/velocity for the two remaining nodes
 */
            ResetSegForces(home, nbr1, &nbr2->myTag, f0seg1[X],
                           f0seg1[Y], f0seg1[Z], 1);

            ResetSegForces(home, nbr2, &nbr1->myTag, f1seg1[X],
                           f1seg1[Y], f1seg1[Z], 1);

        }  /* end "omp parallel for" section */

/*
 *      Last phase of mesh coarsening is loop once more through
 *      the coarsen list and for each 'committed' coarsen, mark
 *      the forces for the nodes (and their neighbors as well) at
 *      the ends of the newly formed segment as obsolete.  This
 *      is done because the force estimates above are not accurate
 *      enough so we need to recalculate more exact forces for
 *      these nodes before (or at the beginning of) the next timestep.
 *
 *      This is not done in the loop above because we need to 
 *      access all the neighbors of the nodes forming the new
 *      segment, and in the loop above, those other nodes may
 *      be involved in a coarsen event in another thread.
 *      We didn't want to loop over them until all coarsened nodes
 *      have been completely removed and the topology is stable again.
 */
#ifdef _OPENMP
#pragma omp parallel for if (coarsenListEnts > 1)
#endif
        for (int i=0; i < coarsenListEnts; i++) {
            int    q;
            Node_t *nbr1, *nbr2, *nbr;

/*
 *          Only need to do the updates for nodes which were 
 *          involved in a coarsen event that was committed.
 */
            if (coarsenList[i].status != REMESH_COMMITTED) {
                continue;
            }

            nbr1 = GetNodeFromTag(home, coarsenList[i].nbr1Tag);
            nbr2 = GetNodeFromTag(home, coarsenList[i].nbr2Tag);

/*
 *          Usually, the segment endpoints will be the neighbors of
 *          the original node, however, it is possible that the
 *          original node to be coarsened out was NOT removed, but
 *          repositioned to one of the neighbors' positions and the
 *          neighbor removed.  In that case, the segment is defined
 *          by the original node and one of its neighbors rather than
 *          the two original neighbors.  If one of the original neighbors
 *          is gone, use the original node instead.
 */
            if (nbr1 == (Node_t *)NULL) {
                nbr1 = GetNodeFromTag(home, coarsenList[i].nodeTag);
            } else if (nbr2 == (Node_t *)NULL) {
                nbr2 = GetNodeFromTag(home, coarsenList[i].nodeTag);
            }

/*
 *          One last safety check... it is possible that coarsening
 *          out the node led to annihilation of the node and the
 *          neighboring nodes (i.e. a triangle collapses in on itself).
 *          If that happens, we don't need to do any more.
 */
            if ((nbr1 == (Node_t *)NULL) || (nbr2 == (Node_t *)NULL)) {
                continue;
            }

            MarkNodeForceObsolete(home, nbr1);
            MarkNodeForceObsolete(home, nbr2);

            for (q = 0; q < nbr1->numNbrs; q++) {
                nbr = GetNodeFromTag(home, nbr1->nbrTag[q]);
                if (nbr == (Node_t *)NULL) continue;
                MarkNodeForceObsolete(home, nbr);
            }

            for (q = 0; q < nbr2->numNbrs; q++) {
                nbr = GetNodeFromTag(home, nbr2->nbrTag[q]);
                if (nbr == (Node_t *)NULL) continue;
                MarkNodeForceObsolete(home, nbr);
            }

/*
 *          If we are enforcing glide planes but allowing them to be
 *          slightly fuzzy, we need to recalculate the glide plane for
 *          the new segment.
 */
            if (param->enforceGlidePlanes &&
                param->allowFuzzyGlidePlanes) {

                RecalcSegGlidePlane(home, nbr1, nbr2, 1);
            }

        }  /* end "omp parallel for" section */

#ifdef DEBUG_LOG_MESH_COARSEN
#ifdef PARALLEL
        int globalCoarsenCnt = 0;
        MPI_Reduce(&localCoarsenCnt, &globalCoarsenCnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        int globalCoarsenCnt = localCoarsenCnt;
#endif
        if (home->myDomain == 0) {
            printf("  Remesh: coarsen count = %d\n", globalCoarsenCnt);
        }
#endif

/*
 *      Free the coarsen list now that we no longer need it.
 */
        if (coarsenList != (CoarsenList_t *)NULL) {
            free(coarsenList);
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       Remesh
 *      Description:    This is just a generic function to set up
 *                      for a remesh operation then select and execute
 *                      the proper remesh function.
 *
 *-----------------------------------------------------------------------*/
void Remesh(Home_t *home)
{
        int     i;
        Node_t  *node;
        Param_t *param;

        param = home->param;

        if (param->remeshRule < 0) return;

#ifdef PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, REMESH_START_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, REMESH_START_BARRIER);
#endif
#endif

        TimerStart(home, REMESH);
        ClearOpList(home);
        InitTopologyExemptions(home);

        switch(param->remeshRule) {
        case 2:
            RemeshRule_2(home);
            break;
        case 3:
            RemeshRule_3(home);
            break;
        default:
            Fatal("Remesh: undefined remesh rule %d", param->remeshRule);
            break;
        }
        TimerStop(home, REMESH);

#ifdef _ARLFEM
/*
 *      If free-surfaces are in use rather than PBC, it is possible that
 *      some of the remesh operations have resulted in formation of
 *      segments between two nodes on the same surface.  Any such 
 *      segments must be cut here... This may orphan nodes, but those
 *      should be taken care of by the RemoveOrphanedNodes() call further
 *      on.
 */
        CutSurfaceSegments(home);
#endif  /* ifdef _ARLFEM */

/*
 *      Send to the neighboring domains, a list of all local
 *      operations (add/delete nodes, relinking of arms,
 *      etc) that may affect nodes in remote domains, and process
 *      the remesh operations from the neighboring domains
 */
        TimerStart(home, SEND_REMESH);
        CommSendRemesh(home);
        TimerStop(home, SEND_REMESH);

        TimerStart(home, FIX_REMESH);
        FixRemesh(home);
        TimerStop(home, FIX_REMESH);

/*
 *      Under certain circumstances, parallel topological changes
 *      can create double links between nodes; links which can not be
 *      detected until after FixRemesh() is called... so, a quick
 *      check has to be done to clean up these potential double-
 *      links here, or they will cause problems later on.  Should
 *      only have to check nodes local to this domain.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            (void)RemoveDoubleLinks(home, node, 0);
            node->flags &= ~NODE_CHK_DBL_LINK;
        }

/*
 *      It is possible that remesh/coarsen has left an orphaned node.
 *      We need to get rid of any such nodes before the exchange of
 *      node information in CommSendGhosts().
 */
        RemoveOrphanedNodes(home);

#ifdef PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, REMESH_END_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, REMESH_END_BARRIER);
#endif
#endif

/*
 *      If memory debugging is enabled, run a consistency check on all
 *      allocated memory blocks looking for corruption in any of the
 *      block headers or trailers.
 */
#ifdef DEBUG_MEM
        ParadisMemCheck();
#endif

        return;
}
