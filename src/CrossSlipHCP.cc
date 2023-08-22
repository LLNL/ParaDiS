/****************************************************************************
 *
 *      Module:         CrossSlipHCP.c
 *
 *      Author:         Adapted from CrossSlipFCC.c
 *                      S. Aubry Feb. 18 2011
 *
 *      Description:    This module contains functions for allowing
 *                      dislocations in HCP materials to cross slip to
 *                      a glide plane other than its original plane
 *
 *                      NOTE:  See comments at the top of the CrossSlip.c
 *                             module for a general description of the
 *                             cross-slip mechanism.
 *
 *      Includes public functions:
 *          CrossSlipHCP()
 *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "CrossSlip.h"

#ifdef DEBUG_CROSSSLIP_EVENTS
static int dbgDom;
#endif


/*---------------------------------------------------------------------------
 *
 *      Function:       CrossSlipHCP
 *
 *      Description:    Examines all nodes local to the domain, determines
 *                      if the node should cross slip, if the node is
 *                      permitted to cross slip, and if so, adjusts node
 *                      positions and segment glide planes as necessary.
 *
 *-------------------------------------------------------------------------*/
void CrossSlipHCP(Home_t *home)
{
        int     i;
        int     opClass, thisDom, numNodes;
        int     eventListSize, zipEventCnt, slipEventIndex;
        int     numBurg;
        int     *numGlissilePlanesPerBurg, *burgFirstPlaneIndex;
        real8   shearModulus;
        real8   noiseFactor, weightFactor;
        real8   eps, thetacrit, sthetacrit, s2thetacrit, areamin;
        real8   (*burgList)[3], (*planeList)[3];
        Param_t *param;
        SlipEvent_t *eventList = (SlipEvent_t *)NULL;

/*
 *      Allocate a block of structures to hold info about cross-slip
 *      events.  We allocate enough structures to cross-slip every local
 *      node, which is overkill, but it means we won't run out of entries
 *      and have to reallocate the block, plus minimizes the number of
 *      operations we'll have to do within the thread-critical sections
 *      of code.
 */
        eventListSize = home->newNodeKeyPtr;
        slipEventIndex = eventListSize;
        zipEventCnt = 0;

        if (eventListSize > 0) {
            int blockSize;
            blockSize = eventListSize * sizeof(SlipEvent_t);
            eventList = (SlipEvent_t *)malloc(blockSize);
        }

#ifdef DEBUG_CROSSSLIP_EVENTS
#ifdef DEBUG_CROSSSLIP_DOMAIN
        dbgDom = DEBUG_CROSSSLIP_DOMAIN;
#else
        dbgDom = -1;
#endif
#endif

        param = home->param;
        thisDom = home->myDomain;
        numNodes = home->newNodeKeyPtr;

        burgList = home->burgData.burgList;
        planeList = home->burgData.planeList;
        numBurg = home->burgData.numBurgVectors;
        numGlissilePlanesPerBurg = home->burgData.numGlissilePlanesPerBurg;
        burgFirstPlaneIndex = home->burgData.burgFirstPlaneIndex;

        eps = 1.0e-06;
        thetacrit = 2.0 / 180.0 * M_PI;
        sthetacrit = sin(thetacrit);
        s2thetacrit = sthetacrit * sthetacrit;
        areamin = param->remeshAreaMin;
        shearModulus= param->shearModulus;


/*
 *      In order to cross-slip a node, the force on the cross-slip
 *      plane must exceed a value based on the force on the current
 *      plane plus an amount that should guarantee the force on
 *      the cross-slip plane is not just in the noise.  The next
 *      two values are used in calculating that needed force.
 */
        noiseFactor=1e-04;
        weightFactor=1.0;

/*
 *      For purposes of determining segment ownership, we need to treat
 *      cross slip operations the same as we do during collision handling.
 */
        opClass = OPCLASS_COLLISION;

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int    nIndex;
            int    threadID, threadIterStart, threadIterEnd;


            GetThreadIterationIndices(numNodes, &threadID, &threadIterStart,
                                      &threadIterEnd);

/*
 *          Each thread loops through a portion of the nodes evaluating
 *          any 2-nodes as possible cross-slip candidates.
 */
            for (nIndex = threadIterStart; nIndex < threadIterEnd; nIndex++) {
                int    burgIndex;
                int    nbr1SegID, nbr2SegID;
                int    seg1_is_screw, seg2_is_screw, bothseg_are_screw;
                real8  burgSize;
                real8  test1, test2, test3, testmax1, testmax2, testmax3;
                real8  burgLab[3], burgCrystal[3];
                real8  fLab[3], fCrystal[3];
                real8  nodep[3], nbr1p[3], nbr2p[3];
                real8  vec1[3], vec2[3], vec3[3];
                Node_t *node, *nbr1, *nbr2;


                if ((node = home->nodeKeys[nIndex]) == (Node_t *)NULL) {
                    continue;
                }

                if (node->numNbrs != 2) {
                    continue;
                }

                if (!HAS_NO_CONSTRAINTS(node->constraint)) {
                    continue;
                }

                burgLab[0] = node->burgX[0];
                burgLab[1] = node->burgY[0];
                burgLab[2] = node->burgZ[0];

                VECTOR_COPY(burgCrystal, burgLab);

/*
 *              If needed, rotate a copy of the burgers vector from the lab
 *              frame to the crystal frame.  Otherwise, lab and crystal frames
 *              are identical.
 */
                if (param->useLabFrame) {
                    Matrix33Vector3Multiply(param->rotMatrixInverse, burgLab,
                                            burgCrystal);
                }

                burgSize = sqrt(DotProduct(burgLab, burgLab));

                NormalizeVec(burgLab);
                NormalizeVec(burgCrystal);

/*
 *              Just a quick sanity check...
 */
                if ((fabs(burgCrystal[X]) < eps) &&
                    (fabs(burgCrystal[Y]) < eps) &&
                    (fabs(burgCrystal[Z]) < eps)) {
                    continue;
                }

/*
 *              Find the index of the current burgers vector in the reference
 *              burgers vector list.  However, if the burgers vector is one of
 *              the junction burgers vectors (i.e. burgIndex > 9) we don't
 *              need to consider this for cross-slip.
 */
                GetBurgIndexHCP(burgCrystal, 1, numBurg, burgList, &burgIndex);

                if (burgIndex > 9) {
                    continue;
                }

/*
 *              We must not change the glide plane for any segment
 *              that is not 'owned' by the node in the current domain.
 *              Otherwise, neighboring domains *could* change the
 *              glide plane for the same segment... this would be bad.
 */
                if ((!DomainOwnsSeg(home,opClass,thisDom,&node->nbrTag[0])) ||
                    (!DomainOwnsSeg(home,opClass,thisDom,&node->nbrTag[1]))) {
                    continue;
                }

                nbr1 = GetNeighborNode(home, node, 0);
                nbr2 = GetNeighborNode(home, node, 1);

                if ((nbr1 == (Node_t *)NULL) ||
                    (nbr2 == (Node_t *)NULL)) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                nbr1SegID = GetArmID(nbr1, node);
                nbr2SegID = GetArmID(nbr2, node);

/*
 *              Get the positions of all three nodes and convert the neighbor
 *              node positions to the PBC coordinates relative to the main
 *              node.  These adjusted positions will be used during the
 *              cross slip algorithm, and if updated during the process, will
 *              be copied back into the respective nodal data structures and
 *              shifted back into the primary image space.
 */
                nodep[X] = node->x; nodep[Y] = node->y; nodep[Z] = node->z;
                nbr1p[X] = nbr1->x; nbr1p[Y] = nbr1->y; nbr1p[Z] = nbr1->z;
                nbr2p[X] = nbr2->x; nbr2p[Y] = nbr2->y; nbr2p[Z] = nbr2->z;

                PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
                            &nbr1p[X], &nbr1p[Y], &nbr1p[Z]);
                PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
                            &nbr2p[X], &nbr2p[Y], &nbr2p[Z]);

/*
 *              If the node is a point on a long screw then we can consider
 *              it for possible cross slip.
 */
                vec1[X] = nbr1p[X] - nbr2p[X];
                vec1[Y] = nbr1p[Y] - nbr2p[Y];
                vec1[Z] = nbr1p[Z] - nbr2p[Z];

                vec2[X] = nodep[X] - nbr1p[X];
                vec2[Y] = nodep[Y] - nbr1p[Y];
                vec2[Z] = nodep[Z] - nbr1p[Z];

                vec3[X] = nodep[X] - nbr2p[X];
                vec3[Y] = nodep[Y] - nbr2p[Y];
                vec3[Z] = nodep[Z] - nbr2p[Z];

/*
 *              Force vector (initially in laboratory frame)
 */
                fLab[X] = node->fX;
                fLab[Y] = node->fY;
                fLab[Z] = node->fZ;

/*
 *              Calculate some test conditions
 */
                test1 = DotProduct(vec1, burgLab);
                test2 = DotProduct(vec2, burgLab);
                test3 = DotProduct(vec3, burgLab);

                test1 = test1 * test1;
                test2 = test2 * test2;
                test3 = test3 * test3;

                testmax1 = DotProduct(vec1, vec1);
                testmax2 = DotProduct(vec2, vec2);
                testmax3 = DotProduct(vec3, vec3);

/*
 *              Set up the tests to see if this dislocation is close enough to
 *              screw to be considered for cross slip.  For a segment to be
 *              close to screw it must be within 2*thetacrit defined above
 */
                seg1_is_screw = ((testmax2 - test2) < (testmax2 * s2thetacrit));
                seg2_is_screw = ((testmax3 - test3) < (testmax3 * s2thetacrit));
                bothseg_are_screw =
                        (((testmax2-test2) < (4.0 * testmax2 * s2thetacrit)) &&
                         ((testmax3-test3) < (4.0 * testmax3 * s2thetacrit)) &&
                         ((testmax1-test1) < (testmax1 * s2thetacrit)));

/*
 *              If either (or both) of the segments is screw, evaluate the
 *              node to see if it should be cross-slipped.
 */
                if (seg1_is_screw || seg2_is_screw || bothseg_are_screw) {
                    int    j, n;
                    int    numGlidePlanes, firstGlidePlaneIndex;
                    int    pinned1, pinned2;
                    int    plane1, plane2, fplane;
                    int    thisEventIndex;
                    real8  fnodeThreshold;
                    real8  L1, L2;
                    real8  tmp;
                    real8  f1dotplane1, f1dotplane2, f1dotplanef;
                    real8  vec1dotb, vec2dotb, vec3dotb, fdotglide;
                    real8  newSegForce[3];
                    real8  segplane1[3], segplane2[3], newplane[3];
                    real8  tmp3[4], tmp3B[4], tmp3C[4];
                    real8  glideDirLab[4][3], glideDirCrystal[4][3];
                    SlipEvent_t *thisEvent;

/*
 *                  Set the force threshold for noise level within the code
 */
                    L1 = sqrt(testmax2);
                    L2 = sqrt(testmax3);

                    fnodeThreshold = noiseFactor * shearModulus * burgSize *
                                     0.5 * (L1 + L2);

/*
 *                  Find which glide planes the segments are on.
 *              
 *                  Use burgers vectors in crystal frame to generate initial
 *                  glide planes in crystal frame.
 *
 */
                    numGlidePlanes  = numGlissilePlanesPerBurg[burgIndex];
                    firstGlidePlaneIndex = burgFirstPlaneIndex[burgIndex];

                    for (n = 0; n < numGlidePlanes; n++) {
                        NormalizedCrossVector(
                                    burgList[burgIndex],
                                    planeList[firstGlidePlaneIndex+n],
                                    glideDirCrystal[n]);
                    }

/*
 *                  In HCP, the first 9 burgers vectors have 4 glissile glide
 *                  planes, but the 10th has only three.  To simplify things,
 *                  zero out the 4th row of the glide plane matrix for burgers
 *                  vector 10 so we can use the same code below regardless of
 *                  whether there are 3 or 4 glide planes.
 */
                    if (numGlidePlanes == 3) {
                        VECTOR_ZERO(glideDirCrystal[3]);
                    }

                    segplane1[X] = node->nx[0];
                    segplane1[Y] = node->ny[0];
                    segplane1[Z] = node->nz[0];

                    segplane2[X] = node->nx[1];
                    segplane2[Y] = node->ny[1];
                    segplane2[Z] = node->nz[1];

/*
 *                  Need copies of the glideDir matrix and the force vector
 *                  in both the lab frame and the crystal frame for later use.
 *                  Also need to rotate the segment planes into the crystal
 *                  frame.
 */
                    if (param->useLabFrame) {
                        int   row;
                        real8 tmpPlane[3];
/*
 *                      Rotate initial glide dir from crystal to lab frame
 */
                        for (row = 0; row < numGlidePlanes; row++) {
                            Matrix33Vector3Multiply(param->rotMatrix,
                                                    glideDirCrystal[row],
                                                    glideDirLab[row]);
                        }
/*
 *                      Rotate segment planes to crystal frame
 */
                        Matrix33Vector3Multiply(param->rotMatrixInverse,
                                                segplane1, tmpPlane);
                        VECTOR_COPY(segplane1, tmpPlane);

                        Matrix33Vector3Multiply(param->rotMatrixInverse,
                                                segplane2, tmpPlane);
                        VECTOR_COPY(segplane2, tmpPlane);

/*
 *                      Rotate force vector to crystal frame
 */
                        Matrix33Vector3Multiply(param->rotMatrixInverse, fLab,
                                                fCrystal);
                    } else {
/*
 *                      Lab and crystal frames are identical, so just copy
 *                      the vectors
 */
                        for (n = 0; n < numGlidePlanes; n++) {
                            VECTOR_COPY(glideDirLab[n], glideDirCrystal[n]);
                        }

                        VECTOR_COPY(fCrystal, fLab);
                    }

                    Matrix43_Vmul(tmp3 , glideDirCrystal, segplane1 );   // tmp3  = glideDirCrystal * segplane1
                    Matrix43_Vmul(tmp3B, glideDirCrystal, segplane2 );   // tmp3B = glideDirCrystal * segplane2
                    Matrix43_Vmul(tmp3C, glideDirCrystal, fCrystal  );   // tmp3C = glideDirCrystal * fCrystal

                    plane1 = 0;
                    plane2 = 0;
                    fplane = 0;

/*
 *                  For HCP there can be up to 4 slip planes for screw
 *                  dislocations
 */
                    for (j = 1; j < numGlidePlanes; j++) {
                        plane1=(fabs(tmp3[j])<fabs(tmp3[plane1])) ? j:plane1;
                        plane2=(fabs(tmp3B[j])<fabs(tmp3B[plane2])) ? j:plane2;
                        fplane=(fabs(tmp3C[j])>fabs(tmp3C[fplane])) ? j:fplane;
                    }

/*
 *                  Calculate the new plane in the lab frame since it will
 *                  only be used from here on out to update the existing
 *                  nodal data which is in the lab frame.
 */
                    NormalizedCrossVector(burgLab, glideDirLab[fplane],
                                          newplane);

                    if ((bothseg_are_screw) && (plane1 == plane2) &&
                        (plane1 != fplane) &&
                        (fabs(tmp3C[fplane]) >
                         (weightFactor*fabs(tmp3C[plane1])+fnodeThreshold))) {
/*
 *                      Both segments are close to screw and the average
 *                      direction is close to screw.
 *
 *                      Determine if the neighbor nodes should be considered
 *                      immobile.
 */
                        pinned1 = NodePinned(home, nbr1, plane1, glideDirLab,
                                             numGlidePlanes);

                        pinned2 = NodePinned(home, nbr2, plane2, glideDirLab,
                                             numGlidePlanes);

                        if (pinned1) {
                            if ((!pinned2) ||
                                ((testmax1-test1)<(eps*eps*burgSize*burgSize))){

                                vec1dotb = DotProduct(vec1, burgLab);
                                vec2dotb = DotProduct(vec2, burgLab);

/*
 *                              If the second neighbor node is pinned, it is
 *                              already perfectly aligned with the first
 *                              neighbor in the screw direction so there is
 *                              no need to move it.  If it is not pinned, go
 *                              ahead and move it.
 */
                                if (!pinned2) {
                                    nbr2p[X] = nbr1p[X] - (vec1dotb*burgLab[X]);
                                    nbr2p[Y] = nbr1p[Y] - (vec1dotb*burgLab[Y]);
                                    nbr2p[Z] = nbr1p[Z] - (vec1dotb*burgLab[Z]);
                                }

                                nodep[X] = nbr1p[X] + (vec2dotb * burgLab[X]);
                                nodep[Y] = nbr1p[Y] + (vec2dotb * burgLab[Y]);
                                nodep[Z] = nbr1p[Z] + (vec2dotb * burgLab[Z]);

                                fdotglide = DotProduct(fLab,
                                                       glideDirLab[fplane]);

                                tmp = areamin / fabs(vec1dotb) * 2.0 *
                                      (1.0 + eps) * Sign(fdotglide);

                                nodep[X] += tmp * glideDirLab[fplane][X];
                                nodep[Y] += tmp * glideDirLab[fplane][Y];
                                nodep[Z] += tmp * glideDirLab[fplane][Z];

/*
 *                              Add the necessary info to the list of normal
 *                              cross-slip events to attempt.  These normal
 *                              (non-zipper) type cross-slip events are added
 *                              from the top of the event list downward.  The
 *                              zipper events are added from the bottome of
 *                              the event list upward.
 */
#ifdef _OPENMP
#pragma omp critical (CROSS_SLIP_EVENT)
#endif
                                {
                                    slipEventIndex -= 1;
                                    thisEventIndex = slipEventIndex;
                                }

                                thisEvent = &eventList[thisEventIndex];

                                thisEvent->node = node;
                                thisEvent->nbrNodes[0] = nbr1;
                                thisEvent->nbrNodes[1] = nbr2;

                                VECTOR_COPY(thisEvent->newNodePos, nodep);
                                VECTOR_COPY(thisEvent->newPlane, newplane);
                                VECTOR_COPY(thisEvent->glideDir,
                                            glideDirLab[fplane]);

                                thisEvent->fdotglide = fdotglide;
                                thisEvent->segIndex = 1;
                                thisEvent->status = CROSS_SLIP_PENDING;
                                thisEvent->flags = SLIP_REPOSITION_NODE;

                                if (!pinned2) {
                                    VECTOR_COPY(thisEvent->newNbrPos, nbr2p);
                                    thisEvent->flags |= SLIP_REPOSITION_NBR;
                                }

                            }

                        } else {
/*
 *                          Neighbor 1 can be moved, so go ahead and add this
 *                          information to the list of potential cross-slip
 *                          events
 */
                            vec1dotb = DotProduct(vec1, burgLab);

                            nbr1p[X] = nbr2p[X] + (vec1dotb * burgLab[X]);
                            nbr1p[Y] = nbr2p[Y] + (vec1dotb * burgLab[Y]);
                            nbr1p[Z] = nbr2p[Z] + (vec1dotb * burgLab[Z]);

                            vec3dotb = DotProduct(vec3, burgLab);

                            nodep[X] = nbr2p[X] + (vec3dotb * burgLab[X]);
                            nodep[Y] = nbr2p[Y] + (vec3dotb * burgLab[Y]);
                            nodep[Z] = nbr2p[Z] + (vec3dotb * burgLab[Z]);

                            fdotglide = DotProduct(fLab, glideDirLab[fplane]);
                            tmp = areamin / fabs(vec1dotb) * 2.0 *
                                   (1.0 + eps) * Sign(fdotglide);

                            nodep[X] += tmp * glideDirLab[fplane][X];
                            nodep[Y] += tmp * glideDirLab[fplane][Y];
                            nodep[Z] += tmp * glideDirLab[fplane][Z];

/*
 *                          Add the necessary info to the list of normal
 *                          cross-slip events to attempt.  These normal
 *                          (non-zipper) type cross-slip events are added
 *                          from the top of the event list downward.  The
 *                          zipper events are added from the bottome of
 *                          the event list upward.
 */
#ifdef _OPENMP
#pragma omp critical (CROSS_SLIP_EVENT)
#endif
                            {
                                slipEventIndex -= 1;
                                thisEventIndex = slipEventIndex;
                            }  /* end "omp critical" section */

                            thisEvent = &eventList[thisEventIndex];

                            thisEvent->node = node;
                            thisEvent->nbrNodes[0] = nbr1;
                            thisEvent->nbrNodes[1] = nbr2;

                            VECTOR_COPY(thisEvent->newNodePos, nodep);
                            VECTOR_COPY(thisEvent->newNbrPos, nbr1p);
                            VECTOR_COPY(thisEvent->newPlane, newplane);
                            VECTOR_COPY(thisEvent->glideDir,
                                        glideDirLab[fplane]);

                            thisEvent->fdotglide = fdotglide;
                            thisEvent->segIndex = 0;
                            thisEvent->status = CROSS_SLIP_PENDING;
                            thisEvent->flags = SLIP_REPOSITION_NODE |
                                               SLIP_REPOSITION_NBR;
                        }

                    } else if ((seg1_is_screw) &&
                               (plane1 != plane2) &&
                               (plane2 == fplane) &&
                               (fabs(tmp3C[fplane]) >
                                (weightFactor*fabs(tmp3C[plane1]) +
                                 fnodeThreshold))){

/*
 *                      Zipper condition met for first segment.  If the first
 *                      neighbor is either not pinned or pinned but already
 *                      sufficiently aligned, this is a potential cross-slip
 *                      event.
 */

                        pinned1 = NodePinned(home, nbr1, plane1, glideDirLab,
                                             numGlidePlanes);

                        if ((!pinned1) ||
                            ((testmax2-test2) < (eps*eps*burgSize*burgSize))) {
/*
 *                          Before 'zippering' a segment, try a quick sanity
 *                          check to see if it makes sense.  If the force on
 *                          the segment to be 'zippered' is not sufficiently
 *                          larger on the new plane than the old plane, leave
 *                          the segment alone.
 */
                            newSegForce[X] = node->armfx[0] +
                                             nbr1->armfx[nbr1SegID];

                            newSegForce[Y] = node->armfy[0] +
                                             nbr1->armfy[nbr1SegID];

                            newSegForce[Z] = node->armfz[0] +
                                             nbr1->armfz[nbr1SegID];

                            f1dotplane1 = fabs(DotProduct(newSegForce,
                                                          glideDirLab[plane1]));
                            f1dotplanef = fabs(DotProduct(newSegForce,
                                                          newplane));

/* FIX ME!  replace with a 'zipperThreshold' as in other cross-slip modules? */
                            if (f1dotplane1 > 0.95 * f1dotplanef) {
                                continue;
                            }

                            if (!pinned1) {
                                vec2dotb = DotProduct(vec2, burgLab);

                                nbr1p[X] = nodep[X] - (vec2dotb * burgLab[X]);
                                nbr1p[Y] = nodep[Y] - (vec2dotb * burgLab[Y]);
                                nbr1p[Z] = nodep[Z] - (vec2dotb * burgLab[Z]);
                            }
/*
 *                          Add the necessary info to the list of zipper
 *                          cross-slip events to attempt.  The zipper
 *                          type cross-slip events are added from the bottom
 *                          of the event list upward.  The normal (non-zipper)
 *                          events are added from the top of the event list
 *                          downward.
 */
#ifdef _OPENMP
#pragma omp critical (CROSS_SLIP_ZIP_EVENT)
#endif
                            {
                                thisEventIndex = zipEventCnt;
                                zipEventCnt += 1;
                            }  /* end "omp critical" section */

                            thisEvent = &eventList[thisEventIndex];

                            thisEvent->flags = SLIP_ZIPPER_EVENT;
                            thisEvent->status = CROSS_SLIP_PENDING;
                            thisEvent->segIndex = 0;
                            thisEvent->node = node;
                            thisEvent->nbrNodes[0] = nbr1;

                            VECTOR_COPY(thisEvent->newPlane, newplane);

                            if (!pinned1) {
                                VECTOR_COPY(thisEvent->newNbrPos, nbr1p);
                                thisEvent->flags |= SLIP_REPOSITION_NBR;
                            }

                        }

                    } else if ((seg2_is_screw) &&
                               (plane1 != plane2) &&
                               (plane1 == fplane) &&
                               (fabs(tmp3C[fplane]) >
                                (weightFactor*fabs(tmp3C[plane2]) +
                                 fnodeThreshold))){
/*
 *                      Zipper condition met for second segment.  If the
 *                      second neighbor is either not pinned or pinned but
 *                      already sufficiently aligned, this is a potential
 *                      cross-slip event.
 */
                        pinned2 = NodePinned(home, nbr2, plane2, glideDirLab,
                                             numGlidePlanes);

                        if ((!pinned2) ||
                            ((testmax2-test2) < (eps*eps*burgSize*burgSize))) {
/*
 *                          Before 'zippering' a segment, try a quick sanity
 *                          check to see if it makes sense.  If the force on
 *                          the segment to be 'zippered' is not sufficiently
 *                          larger on the new plane than the old plane, leave
 *                          the segment alone.
 */
                            newSegForce[X] = node->armfx[1] +
                                             nbr2->armfx[nbr2SegID];

                            newSegForce[Y] = node->armfy[1] +
                                             nbr2->armfy[nbr2SegID];

                            newSegForce[Z] = node->armfz[1] +
                                             nbr2->armfz[nbr2SegID];

                            f1dotplane2 = fabs(DotProduct(newSegForce,
                                                          glideDirLab[plane2]));
                            f1dotplanef = fabs(DotProduct(newSegForce,
                                                          newplane));

/* FIX ME!  replace with a 'zipperThreshold' as in other cross-slip modules? */
                            if (f1dotplane2 > 0.95 * f1dotplanef) {
                                continue;
                            }

                            if (!pinned2) {
                                vec3dotb = DotProduct(vec3, burgLab);

                                nbr2p[X] = nodep[X] - (vec3dotb * burgLab[X]);
                                nbr2p[Y] = nodep[Y] - (vec3dotb * burgLab[Y]);
                                nbr2p[Z] = nodep[Z] - (vec3dotb * burgLab[Z]);
                            }

/*
 *                          Add the necessary info to the list of zipper
 *                          cross-slip events to attempt.  The zipper
 *                          type cross-slip events are added from the bottom
 *                          of the event list upward.  The normal (non-zipper)
 *                          events are added from the top of the event list
 *                          downward.
 */
#ifdef _OPENMP
#pragma omp critical (CROSS_SLIP_ZIP_EVENT)
#endif
                            {
                                thisEventIndex = zipEventCnt;
                                zipEventCnt += 1;
                            }  /* end "omp critical" section */

                            thisEvent = &eventList[thisEventIndex];

                            thisEvent->flags = SLIP_ZIPPER_EVENT;
                            thisEvent->status = CROSS_SLIP_PENDING;
                            thisEvent->segIndex = 1;
                            thisEvent->node = node;
                            thisEvent->nbrNodes[1] = nbr2;

                            VECTOR_COPY(thisEvent->newPlane, newplane);

                            if (!pinned2) {
                                thisEvent->flags |= SLIP_REPOSITION_NBR;
                                VECTOR_COPY(thisEvent->newNbrPos, nbr2p);
                            }

                        }  /* end if (!pinned2) */

                    }  /* end if (seg2isScrew) */

                }  /* end of check for any screws */

            }  /* end for (nIndex = threadIterStart; ...) */

        }  /* end "omp parallel" section */

/*
 *      It may be that due to conflicts we can not perform all the
 *      desired cross-slip events, however, if the conflict is between
 *      a 'zipper' event and a normal cross-slip event, we want to
 *      preferentially perform the 'zipper' event. So, we first
 *      loop through the 'zipper' events.  Reminder, the zipper events
 *      are in the first <zipEventCnt> elements of the event list.
 *
 *      NOTE: These loops are not explicitly threaded here for two
 *            reasons.  First, so much of the code in the loops
 *            would have to be in 'critical' sections that
 *            the lock contention would probably slow things down.
 *            Second, the primary expense is the force calculations
 *            in the loop over the non-zipper events, and those
 *            force calculations are already threaded within the
 *            SetOneNodeForce() function.
 */
        for (i = 0; i < zipEventCnt; i++) {
            int segIndex;
            Node_t *node, *nbrNode;

            if (!ProceedWithZipEvent(eventList, zipEventCnt, i)) {
                continue;
            }

            segIndex = eventList[i].segIndex;
            node = eventList[i].node;
            nbrNode = eventList[i].nbrNodes[segIndex];

#ifdef DEBUG_CROSSSLIP_EVENTS
            if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                char *zip1 = "Zipping first segment";
                char *zip2 = "Zipping second segment";
                DumpCrossSlipEvent(node, eventList[i].newPlane,
                                   (segIndex == 0 ? zip1 : zip2));
            }
#endif

            ResetGlidePlane(home, eventList[i].newPlane, &node->myTag,
                            &nbrNode->myTag, 1);

/*
 *          If the neighboring node was repositioned, shift the new coordinates
 *          back into the primary image space and update the corresponding
 *          node structure.  The operation will also be sent to the remote
 *          domains for processing.
 */
            if (eventList[i].flags & SLIP_REPOSITION_NBR) {
                real8 newNbrPos[3];
                VECTOR_COPY(newNbrPos, eventList[i].newNbrPos);
                FoldBox(param, &newNbrPos[X], &newNbrPos[Y], &newNbrPos[Z]);
                RepositionNode(home, newNbrPos, &nbrNode->myTag, 1);
            }
        }

/*
 *      Now loop through the normal cross-slip events.  Reminder: these
 *      are contained within elements slipEventIndex thru eventListSize-1
 *      of the event list.
 */
        for (i = slipEventIndex; i < eventListSize; i++) {
            int    nbr1ArmID, nbr2ArmID, segIndex;
            real8  newfdotglide;
            real8  *nbrPosOrig;
            real8  newPlane[3];
            real8  newNodePos[3], newNbrPos[3];
            real8  nodePosOrig[3], nbr1PosOrig[3], nbr2PosOrig[3];
            real8  newForce[3];
            real8  segForceOrig[4][3];
            Node_t *node, *nbr1, *nbr2, *nbrNode;

/*
 *          Since a cross-slip event can modify node positions and
 *          segments data, we don't want to attempt a cross-slip that
 *          involves nodes/segments that were modified in a prior
 *          cross-slip event this cycle.
 */
            if (!ProceedWithCrossSlip(eventList, zipEventCnt, slipEventIndex,
                                      eventListSize, i)) {
                continue;
            }

            segIndex = eventList[i].segIndex;
            nbrNode = eventList[i].nbrNodes[segIndex];

            node = eventList[i].node;
            nbr1 = eventList[i].nbrNodes[0];
            nbr2 = eventList[i].nbrNodes[1];

/*
 *          It's possible that once we've shifted nodes and recalculated
 *          forces after the cross-slip event, it will turn out that the
 *          original configuration was preferable.  We need to preserve
 *          the original configuration so that it can be restored if this
 *          situation occurs.
 */
            nbr1ArmID = GetArmID(nbr1, node);
            nbr2ArmID = GetArmID(nbr2, node);

            SaveCrossSlipInfo(node, nbr1, nbr2, nbr1ArmID, nbr2ArmID,
                              segForceOrig, nodePosOrig, nbr1PosOrig,
                              nbr2PosOrig);

            nbrPosOrig = (segIndex == 0 ? nbr1PosOrig : nbr2PosOrig);

/*
 *          Perform the cross-slip event and recalculate the nodal
 *          force.  If it appears the node will not continue to move out
 *          on the new plane, skip the cross-slip event and
 *          restore the old configuration.
 */
            VECTOR_COPY(newNodePos, eventList[i].newNodePos);
            ResetPosition(param, node, newNodePos);

            if (eventList[i].flags & SLIP_REPOSITION_NBR) {
                VECTOR_COPY(newNbrPos, eventList[i].newNbrPos);
                ResetPosition(param, nbrNode, newNbrPos);
            }

            SetOneNodeForce(home, node);

            newForce[X] = node->fX;
            newForce[Y] = node->fY;
            newForce[Z] = node->fZ;

            newfdotglide = DotProduct(newForce, eventList[i].glideDir);

            if ((Sign(newfdotglide) * Sign(eventList[i].fdotglide)) < 0.0) {
                ResetPosition(param, node, nodePosOrig);
                if (eventList[i].flags & SLIP_REPOSITION_NBR) {
                    ResetPosition(param, nbrNode, nbrPosOrig);
                }
                RestoreCrossSlipForce(node, nbr1, nbr2,
                                      nbr1ArmID, nbr2ArmID,
                                      segForceOrig);
                eventList[i].status = CROSS_SLIP_BACKED_OFF;
                continue;
            }

            VECTOR_COPY(newPlane, eventList[i].newPlane);

#ifdef DEBUG_CROSSSLIP_EVENTS
            if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                char msg[128];

                if (eventList[i].flags & SLIP_REPOSITION_NBR)  {
                    sprintf(msg, "%s neighbor repositioned",
                            segIndex == 0 ? "first" : "second");
                } else {
                    sprintf(msg, "neither neighbor repositioned");
                }
                DumpCrossSlipEvent(node, newPlane, msg);
            }
#endif

            ResetGlidePlane(home, newPlane, &node->myTag, &nbr1->myTag, 1);
            ResetGlidePlane(home, newPlane, &node->myTag, &nbr2->myTag, 1);

            FoldBox(param, &newNodePos[X], &newNodePos[Y], &newNodePos[Z]);
            RepositionNode(home, newNodePos, &node->myTag, 1);

            if (eventList[i].flags & SLIP_REPOSITION_NBR) {
                FoldBox(param, &newNbrPos[X], &newNbrPos[Y], &newNbrPos[Z]);
                RepositionNode(home, newNbrPos, &nbrNode->myTag, 1);
            }

/*
 *          Update velocity of the moved nodes
 */
            int mobError=0;
            MobArgs_t  mobArgs;
            mobError   = EvaluateMobility(home, node, &mobArgs);
            mobError  |= EvaluateMobility(home, nbr1, &mobArgs);
            mobError  |= EvaluateMobility(home, nbr2, &mobArgs);
        }
/*
 *      and free the event list before returning to the caller.
 */
        if (eventList != (SlipEvent_t *)NULL) {
            free(eventList);
        }

        return;
}
