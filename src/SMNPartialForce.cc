/**************************************************************************
 *
 *      Module:       SMNPartialForce.c
 *      Description:  Contains functions needed for calculating the 
 *                    contribution to forces on a node from all locally
 *                    owned segments as well as self-force, pk force, etc
 *                    for segments of the node that are owned locally.
 *
 *                    NOTE: These partial forces are needed only when using
 *                    the parallel multi-node splitting mechanism.
 *
 *      Public functions:
 *          BuildSMNLocalSegList()
 *          SMNPartialForce()
 *
 *      Private functions:
 *          BuildCellNbrList()
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "SMN.h"

/*-------------------------------------------------------------------------
 *
 *      Function:    BuildSMNLocalSegList
 *      Description: Build a list of locally owned segments.
 *                   This list is for use with SplitMultiNodeParallel() in
 *                   which partial forces on the split multinodes are
 *                   calculated by multiple domains then communicated
 *                   to the domain owning the multi-node.
 *
 *      Parameters:
 *          segList:OUT    Location in which to return pointer to an array
 *                         of segments (node pairs) owned by the local
 *                         domain.  Each segment in this list will be owned
 *                         by the first of the segment's nodes.  NOTE: caller
 *                         is responsible for freeing any <segList> returned.
 *          numSegs:OUT    Location in which to return the number of segs
 *                         placed onto <segList>
 *
 *------------------------------------------------------------------------*/
void BuildSMNLocalSegList(Home_t *home, SegNodes_t **segList, int *numSegs)
{
        int        nodeIndex, numNodes;
        int        segListSize;
        SegNodes_t *list = (SegNodes_t *)NULL;

        *segList = (SegNodes_t *)NULL;
        *numSegs = 0;

        numNodes = home->newNodeKeyPtr;
        segListSize = 0;

/*
 *      Allocate an initial array to hold all the locally owned
 *      segments (node pairs).  If necessary the list size can be
 *      increased later.  If there are no local nodes, there's no
 *      need to continue;
 */
        if ((segListSize = home->newNodeKeyPtr * 1.2) == 0) {
            return;
        }

        list = (SegNodes_t *)malloc(segListSize * sizeof(SegNodes_t));

        for (nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {
            int    numNbrs, segID;
            Node_t *node;

            if ((node = home->nodeKeys[nodeIndex]) == (Node_t *)NULL) {
                continue;
            }

            if (home->myDomain != node->myTag.domainID) {
                continue;
            }

            numNbrs = node->numNbrs;

            for (segID = 0; segID < numNbrs; segID++) {
                Node_t *nbrNode;

                nbrNode = GetNeighborNode(home, node, segID);

                if (nbrNode == (Node_t *)NULL) {
                    continue;
                }

/*
 *              Skip the segment if the neighboring node owns it
 */
                if (!NodeOwnsSeg(home, node, nbrNode)) {
                    continue;
                }

/*
 *              Add the segment to the list (reallocate mem for list
 *              if we've run out of space)
 */
                if (*numSegs >= segListSize) {
                    int allocSize;

                    segListSize += home->newNodeKeyPtr;
                    allocSize = segListSize * sizeof(SegNodes_t);
                    list = (SegNodes_t *)realloc(list, allocSize);
                }

                list[*numSegs].node1 = node;
                list[*numSegs].node2 = nbrNode;

                *numSegs += 1;

            }  /* end for (segID = 0; segID < numNbrs; ...) */

        }  /* end for (nodeIndex = 0; nodeIndex < numNodes; ...) */

        *segList = list;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    BuildCellNbrList
 *      Description: Given the XYZ indices of a cell (these indices are
 *                   shifted to account for ghost cells), create a list
 *                   of all cells including or immediately neighboring
 *                   the specified cell.
 *
 *                   This list is for use with the parallel version of
 *                   SplitMultiNode() in which partial forces on the split
 *                   multinodes are calculated by multiple domains
 *
 *      Parameters:
 *          cell{XYZ}Index IN: Indices of a cell. These indices should
 *                             already be shifted up by 1 to account for
 *                             ghost cells.
 *          nbrCellList   OUT: Array of cell IDs.  If periodic boundaries 
 *                             are in use, there will always be 27 ID's
 *                             in this list.  If one or more boundaries are
 *                             free surfaces, the number of cells may be
 *                             less than 27.
 *          numNbrCells   OUT: Number of cell IDs returned in <nbrCellList>.
 *
 *------------------------------------------------------------------------*/
static void BuildCellNbrList(Home_t *home, int cellXIndex, int cellYIndex,
                             int cellZIndex, int nbrCellList[27],
                             int *numNbrCells)
{
        int     i, j, k;
        int     xMin, yMin, zMin;
        int     xMax, yMax, zMax;
        int     xIndex, yIndex, zIndex;
        Param_t *param;

        param = home->param;

        *numNbrCells = 0;

        xMin = cellXIndex - 1;
        yMin = cellYIndex - 1;
        zMin = cellZIndex - 1;

        xMax = cellXIndex + 1;
        yMax = cellYIndex + 1;
        zMax = cellZIndex + 1;

        if (param->xBoundType != Periodic) {
            xMin = MAX(xMin, 1);
            xMax = MIN(xMax, param->nXcells);
        }

        if (param->yBoundType != Periodic) {
            yMin = MAX(yMin, 1);
            yMax = MIN(yMax, param->nYcells);
        }

        if (param->xBoundType != Periodic) {
            zMin = MAX(zMin, 1);
            zMax = MIN(zMax, param->nZcells);
        }

        for (i = xMin; i <= xMax; i++) {

            xIndex = i;

            if (xIndex < 1) xIndex = param->nXcells;
            if (xIndex > param->nXcells) xIndex = 1;

            for (j = yMin; j <= yMax; j++) {

                yIndex = j;

                if (yIndex < 1) yIndex = param->nYcells;
                if (yIndex > param->nYcells) yIndex = 1;

                for (k = zMin; k <= zMax; k++) {
                    int baseCellID;

                    zIndex = k;

                    if (zIndex < 1) zIndex = param->nZcells;
                    if (zIndex > param->nZcells) zIndex = 1;

                    baseCellID = EncodeCellIdx(home, xIndex, yIndex, zIndex);
                    nbrCellList[*numNbrCells] = baseCellID;
                    *numNbrCells += 1;
                }
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    SMNPartialForce
 *      Description: Calculate the contribution to forces on the segments
 *                   of a specified node from all locally owned segments as
 *                   well as self-force, pk force, etc for segments of the
 *                   node that are owned locally.  Segment specific forces
 *                   are then returned to the caller in the provided 
 *                   structure.
 *
 *      Parameters:
 *          node       IN:  Pointer to the node
 *          numSegs    IN:  Number of locally owned segments
 *          segList    IN:  Array of <numSegs> structs specifying nodal
 *                          endpoints of segments
 *          nodeForces OUT: Pointer to struct in which to return the
 *                          calculated segment specific forces.
 *
 *------------------------------------------------------------------------*/
void SMNPartialForce(Home_t *home, Node_t *node, int numSegs,
                     SegNodes_t *segList, SMN_Force_t *nodeForces)
{
        int     segID, numNbrs;
        int     locallyOwnedSeg;
        real8   a, MU, NU, Ecore, eps;
        real8   extStress[3][3];
        Param_t *param;

        eps   = 1.0e-06;

        param = home->param;

        a     = param->rc;
        MU    = param->shearModulus;
        NU    = param->pois;
        Ecore = param->Ecore;

        extStress[0][0] = param->appliedStress[0];
        extStress[1][1] = param->appliedStress[1];
        extStress[2][2] = param->appliedStress[2];
        extStress[1][2] = param->appliedStress[3];
        extStress[2][0] = param->appliedStress[4];
        extStress[0][1] = param->appliedStress[5];
        extStress[2][1] = extStress[1][2];
        extStress[0][2] = extStress[2][0];
        extStress[1][0] = extStress[0][1];

/*
 *      Loop through all the node's segments calculating forces
 *      on the segments from all segments in the local segment list
 *      that are owned by either the cell owning the target segment or
 *      an immediate neighbor of that cell.
 *
 *      The domain owning the target node will be responsible for 
 *      calculating the segment self force, PK force, osmotic force
 *      and other non-seg/seg force components.
 */
        numNbrs = node->numNbrs;

        for (segID = 0; segID < numNbrs; segID++) {
            int    targetOwnsSeg;
            int    armID12, armID21;
            int    cellX, cellY, cellZ;
            int    nbrCellCount;
            int    nbrCellList[27];
            real8  dx, dy, dz;
            real8  x1, y1, z1;
            real8  x2, y2, z2;
            real8  xCenter, yCenter, zCenter;
            real8  burg1[3];
            real8  f1Seg[3], f2Seg[3];
            Cell_t *cell;
            Node_t *node1, *node2, *nbrNode;

            nbrNode = GetNodeFromTag(home, node->nbrTag[segID]);

            if (nbrNode == (Node_t *)NULL) {
                continue;
            }

/*
 *          For a segment node1--node2, alwys point node1 to the node
 *          owning the segment.
 */
            if (NodeOwnsSeg(home, node, nbrNode)) {
                targetOwnsSeg = 1;
                node1 = node;
                node2 = nbrNode;
                armID12 = segID;
                armID21 = GetArmID(nbrNode, node);
            } else {
                targetOwnsSeg = 0;
                node1 = nbrNode;
                node2 = node;
                armID12 = GetArmID(node1, node2);
                armID21 = segID;
            }

            locallyOwnedSeg = (home->myDomain == node1->myTag.domainID);

/*
 *          We need the indices of the cell containing node1 in several
 *          places further on, so get that now.  If for some reason, the
 *          node has not yet been associated with a cell, we can't
 *          calculate all the forces, so skip it.
 */
            if (node1->cellIdx < 0) {
                continue;
            }

            cell = LookupCell(home, node1->cellIdx);

/*
 *          If the current domain does not have info on this cell,
 *          calculate the indices from the cellID and calculate the
 *          center of the cell.
 */
            if (cell == (Cell_t *)NULL) {

                DecodeCellIdx(home, node1->cellIdx, &cellX, &cellY, &cellZ);

                FindCellCenter(param, (real8)(cellX-1), (real8)(cellY-1),
                           (real8)(cellZ-1), 2, &xCenter, &yCenter,
                           &zCenter);
            } else {

                cellX = cell->xIndex;
                cellY = cell->yIndex;
                cellZ = cell->zIndex;

                xCenter = cell->center[X];
                yCenter = cell->center[Y];
                zCenter = cell->center[Z];
            }

/*
 *          When we loop through the local segment list later, we need
 *          to weed out segments that are not owned by cells neighboring
 *          the cell owning the target segment, so build a full list of
 *          base cell ID's for cells neighboring the target cell.
 */
            BuildCellNbrList(home, cellX, cellY, cellZ, nbrCellList,
                             &nbrCellCount);

            burg1[X] = node1->burgX[armID12];
            burg1[Y] = node1->burgY[armID12];
            burg1[Z] = node1->burgZ[armID12];

            x1 = node1->x;
            y1 = node1->y;
            z1 = node1->z;

            dx = node2->x - x1;
            dy = node2->y - y1;
            dz = node2->z - z1;

            ZImage(param, &dx, &dy, &dz);

            x2 = x1 + dx;
            y2 = y1 + dy;
            z2 = z1 + dz;

/*
 *          If the segment is zero-length, we cannot calculate forces
 *          on it, so skip it.
 */
            if (((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) < eps) {
                continue;
            }

            VECTOR_ZERO(f1Seg);
            VECTOR_ZERO(f2Seg);

/*
 *          If the multinode is owned by this domain, calculate the self-force,
 *          osmotic force, PK force, and other segment specific forces.
 */
            if (locallyOwnedSeg) {
                real8 f1[3], f2[3];
                real8 sigb[3], totRemSig[3][3];

/*
 *              If elastic interaction is not enabled, use a simple line
 *              tension model for calculating forces on the segment. (i.e.
 *              no segment/segment forces, no remote forces, etc)
 */
                if (!param->elasticinteraction) {

                    LineTensionForce(home, node1->x, node1->y, node1->z,
                                     node2->x, node2->y, node2->z,
                                     burg1[X], burg1[Y], burg1[Z], f1, f2);

                    if (targetOwnsSeg) {
                        VECTOR_ADD(nodeForces->segForce[segID].f1, f1);
                        VECTOR_ADD(nodeForces->segForce[segID].f2, f2);
                    } else {
                        VECTOR_ADD(nodeForces->segForce[segID].f1, f2);
                        VECTOR_ADD(nodeForces->segForce[segID].f2, f1);
                    }

                    continue;
                }

/*
 *              Add in segment force due to self-stress
 */
                SelfForce(home, 0, MU, NU, burg1[X], burg1[Y], burg1[Z],
                          x1, y1, z1, x2, y2, z2, a, Ecore, f1, f2);

                VECTOR_ADD(f1Seg, f1);
                VECTOR_ADD(f2Seg, f2);

/*
 *              Now PK force from external stress
 */
                ExtPKForce(extStress, burg1[X], burg1[Y], burg1[Z],
                           x1, y1, z1, x2, y2, z2, f1, f2);

                VECTOR_ADD(f1Seg, f1);
                VECTOR_ADD(f2Seg, f2);

/*
 *              Include osmotic forces if they are enabled
 */
                if (param->vacancyConcEquilibrium > 0.0) {

                    OsmoticForce(home, x1, y1, z1, x2, y2, z2,
                                 burg1[X], burg1[Y], burg1[Z], f1, f2);

                    VECTOR_ADD(f1Seg, f1);
                    VECTOR_ADD(f2Seg, f2);
                }

#ifdef ESHELBY
/*
 *              If the simulation includes eshelby inclusions, we need to add
 *              the force due to the inclusions
 */
                if (param->enableInclusions) 
                {
                   real8  pos1[3], pos2[3];
                   
                   int  inclusionCnt      = 0;
                   int  inclusionListSize = 0;
                   int *inclusionList     = (int *)NULL;
                   
                   pos1[X] = x1;
                   pos1[Y] = y1;
                   pos1[Z] = z1;
                   
                   pos2[X] = x2;
                   pos2[Y] = y2;
                   pos2[Z] = z2;
                
/*
 *                  If this function is called for two connected nodes, we
 *                  have to be sure we don't search for intersecting particles
 *                  multiple times, so only add the segment/particle
 *                  intersections if the segment is not yet on the list.
 *                  (Note: list may be unsorted at this point, so we have
 *                  to use a lookup function that handles an unsorted list)
 */
                   SegPartIntersect_t *intersection = SegPartListUnsortedLookup(home, &node1->myTag, &node2->myTag);
                   int                 addToList    = (intersection==(SegPartIntersect_t *)NULL?1:0);
                   
/*
 *                  Loop through the cell and all it's immediate neighbors.
 *                  (<nbrCellCnt> and <nbrCellList> are a little misleading
 *                  since they include the neighbors and the local cell)
 */
                   for (int k=0; k < nbrCellCount; k++) {
                      int    inclusionIndex, inclusionCellID;
                      Cell_t *inclusionCell;
                      
/*
 *                      If the cell is a ghost cell, alter the cell pointer to
 *                      the corresponding base cell.  Needed in order to pick
 *                      up the proper inclusionQ.
 */
                      inclusionCellID = nbrCellList[k];
                      inclusionCell = LookupCell(home, inclusionCellID);
                   
                      if (inclusionCell == (Cell_t *)NULL) {
                         continue;
                      }
                      
                      if (inclusionCell->baseIdx >= 0) {
                         int inclusionCellID;
                         inclusionCellID = inclusionCell->baseIdx;
                         inclusionCell = LookupCell(home, inclusionCellID);
                      }
                      
/*
 *                      Add the inclusion indices for all inclusions in this
 *                      cell to the list of inclusions for which we have to
 *                      calculate interaction forces for the current segment
 */
                      inclusionIndex = inclusionCell->inclusionQ;
                      
                      while (inclusionIndex >= 0) {
                         EInclusion_t *inclusion;
                         
                         inclusion = &home->eshelbyInclusions[inclusionIndex];
                         
#ifdef ESHELBYFORCE
                         // Force from Eshelby is acting on all segments and
                         // not only the ones that touch the particles.
                         int ninc = 0, isIntersecting;
                         real8 newPos1[3], newPos2[3];
                         real8 ratio[2];
                         isIntersecting = IntersectionInformation(home, pos1, pos2, inclusion, 
                            newPos1, newPos2, &ninc,
                            ratio);
                         
                         
                         // Eshelby force applies on/effects all segments in the simulation and not on 
                         // only the segments intersecting the particles!     
                         EshelbyForce(home, inclusion->position, inclusion->radius, 
                            inclusion->strainField, newPos1, newPos2, 
                            burg1, f1, f2);
                         
                         VECTOR_ADD(f1Seg, f1);
                         VECTOR_ADD(f2Seg, f2);

#endif   //ESHELBYFORCE


                         if (inclusionCnt == inclusionListSize) {
                            inclusionListSize += 1000;
                            inclusionList = (int *)realloc(inclusionList,
                                                           sizeof(int) * inclusionListSize);
                         }
                         
                         inclusionList[inclusionCnt] = inclusionIndex;
                         inclusionIndex = inclusion->nextInCell;
                         
                         inclusionCnt++;
                      }
                      
                   }  /* end for (k = 0; k <= nbrCellCount; ...) */
                   
/*
 *                  Now we have the list of all inclusions which directly
 *                  interact with the node1/node2 segment, so calculate the
 *                  force on the segment from each inclusion.
 */
#ifdef _OPENMP
#pragma omp parallel
#endif
                   {
                      int   m, threadID, threadIterStart, threadIterEnd;
                      real8 f1Thread[3], f2Thread[3];
                      
                      VECTOR_ZERO(f1Thread);
                      VECTOR_ZERO(f2Thread);
                   
                      GetThreadIterationIndices(inclusionCnt, &threadID,
                                                &threadIterStart,
                                                &threadIterEnd);
                      
/*
 *                      Loop through all the inclusions with which this
 *                      segment directly interacts and compute the force
 *                      on the segment from each inclusion.
 */
                      for (m = threadIterStart; m < threadIterEnd; m++) {
                         int  isIntersecting;
                         EInclusion_t *inclusion;
                         
                         inclusion = &home->eshelbyInclusions[inclusionList[m]];
                         
                         int ninc = 0;
                         real8 newPos1[3], newPos2[3];
                         real8 ratio[2];
                         isIntersecting = IntersectionInformation(home, pos1, pos2, inclusion, 
                                                                  newPos1, newPos2, &ninc,
                                                                  ratio);
                         
                         if (isIntersecting)
                         {
#ifdef ESHELBYCORE
                            real8 f1Eshelby[3], f2Eshelby[3];
                            SegPartCoreForce(home,inclusion->position, inclusion->radius, 
                                             inclusion->rotation, newPos1, newPos2, 
                                             burg1, f1Eshelby, f2Eshelby);
                            
                            
                            AddtoArmForce(node1, armID12, f1Eshelby);
                            AddtoArmForce(node2, armID21, f2Eshelby);
#endif //ESHELBYCORE


#ifdef ESHELBYSTEP
                            real8 f1EshelbyStep[3], f2EshelbyStep[3];
                            InclusionStepForce(home->param->IncStepEnergy, inclusion->position, 
                                               inclusion->radius, inclusion->rotation, 
                                               newPos1, newPos2,burg1, f1EshelbyStep, f2EshelbyStep);
                            
                            
                            AddtoArmForce(node1, armID12, f1EshelbyStep);
                            AddtoArmForce(node2, armID21, f2EshelbyStep);
#endif //ESHELBYSTEP



/*
 *                          If the segment intersects the inclusion and we
 *                          need to track such things for the mobility function,
 *                          add the segment segment/particle intersection
 *                          list.
 */
                         if (addToList)
                            SegPartListUpdate(home, &intersection,
                                              &node1->myTag, &node2->myTag,
                                              inclusionList[m],
                                              inclusion->id);
                         
                         }
                         
                      }
                      
/*
 *                      Add the eshelby force from each thread to the segment
 *                      forces.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_SMN_UPDATE_FORCE_A)
#endif
                      {
                         VECTOR_ADD(f1Seg, f1Thread);
                         VECTOR_ADD(f2Seg, f2Thread);
                      }
                      
                   }  /* end omp parallel section */
                   
                   if (inclusionList != (int *)NULL) {
                      free(inclusionList);
                   }
                }  /* if (param->enableInclusions) */
#endif  // ESHELBY & ESHELBYFORCE


#ifdef ESHELBY
                if ((param->fmEnabled) || (param->eshelbyfmEnabled)) 
#else
                if  (param->fmEnabled)
#endif
                {
                    RemoteForceOneSeg(home, node1, node2, f1, f2);

                    VECTOR_ADD(f1Seg, f1);
                    VECTOR_ADD(f2Seg, f2);

                }

/*
 *              If we're not using the FMM, find the segment midpoint and use
 *              that as the point at which to calculate the stress from remote
 *              segments.
 */
                if (param->fmEnabled == 0) {
                    real8 xm, ym, zm;

                    xm = x1 + (dx * 0.5);
                    ym = y1 + (dy * 0.5);
                    zm = z1 + (dz * 0.5);

                    GetFieldPointStressRem(home, xm, ym, zm,
                                           cellX-1, cellY-1, cellZ-1,
                                           totRemSig);

                    sigb[0] = totRemSig[0][0] * burg1[X] +
                              totRemSig[0][1] * burg1[Y] +
                              totRemSig[0][2] * burg1[Z];

                    sigb[1] = totRemSig[1][0] * burg1[X] +
                              totRemSig[1][1] * burg1[Y] +
                              totRemSig[1][2] * burg1[Z];

                    sigb[2] = totRemSig[2][0] * burg1[X] +
                              totRemSig[2][1] * burg1[Y] +
                              totRemSig[2][2] * burg1[Z];

                    VECTOR_COPY(&node1->sigbRem[armID12*3], sigb);
                    VECTOR_COPY(&node2->sigbRem[armID21*3], sigb);
#ifndef _ARLFEM
/*
 *                  If we are not linked with FEM, calculate the PK force
 *                  from remote segments and add it in.
 */
                    PKForce(sigb, x1, y1, z1, x2, y2, z2, f1, f2);

                    VECTOR_ADD(f1Seg, f1);
                    VECTOR_ADD(f2Seg, f2);
#endif
                }  /* if not fmEnabled */

#ifdef _ARLFEM
/*
 *              With free surfaces, we also have to calculate the force from
 *              all "virtual" segments in the simulation space.  Virtual
 *              segments are imaginary segments extending outward from any
 *              segment intersecting the surface.  The segments are treated
 *              as if they extend from the surface to infinity, but are in
 *              reality long but finite segments.
 */
                {
                    real8 pos1[3], pos2[3];

                    pos1[X] = x1;
                    pos1[Y] = y1;
                    pos1[Z] = z1;

                    pos2[X] = x2;
                    pos2[Y] = y2;
                    pos2[Z] = z2;

                    VirtualSegForce(home, pos1, pos2, burg1, f1, f2);

                    VECTOR_ADD(f1Seg, f1);
                    VECTOR_ADD(f2Seg, f2);
                }
#endif

            }  /* end if (locallyOwnedSeg) */

#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                int   m, threadID, threadIterStart, threadIterEnd;
                real8 f1Thread[3], f2Thread[3];

                VECTOR_ZERO(f1Thread);
                VECTOR_ZERO(f2Thread);

                GetThreadIterationIndices(numSegs, &threadID,
                                          &threadIterStart,
                                          &threadIterEnd);

/*
 *              Loop through all the locally owned segments looking for any
 *              that are owned by the cell containing node1 (or any cell
 *              neighboring that cell).
 */
                for (m = threadIterStart; m < threadIterEnd; m++) {
                    int        k, segCellID;
                    int        armID34;
                    real8      x3, y3, z3;
                    real8      x4, y4, z4;
                    real8      fx1, fy1, fz1;
                    real8      fx2, fy2, fz2;
                    real8      fx3, fy3, fz3;
                    real8      fx4, fy4, fz4;
                    real8      burg2[3];
                    Node_t     *node3, *node4;
                    SegNodes_t *seg;

                    seg = &segList[m];
                
/*
 *                  Be sure not to calculate force from a segment with itself.
 */
                    if ((node1->myTag.domainID == seg->node1->myTag.domainID) &&
                        (node1->myTag.index    == seg->node1->myTag.index   ) &&
                        (node2->myTag.domainID == seg->node2->myTag.domainID) &&
                        (node2->myTag.index    == seg->node2->myTag.index   )) {
                        continue;
                    }

                    segCellID = segList[m].node1->cellIdx;

/*
 *                  We only want to calculate forces on the target segment 
 *                  from segments owned by the same or immediately neighboring
 *                  cells, so check for that here.
 */
                    for (k = 0; k < nbrCellCount; k++) {
                        if (segCellID == nbrCellList[k]) {
                            break;
                        }
                    }

                    if (k >= nbrCellCount) {
                        continue;
                    }

/*
 *                  Calculate the force on the target segment from the locally
 *                  owned segment and increment the local partial force sum on
 *                  the target segment.
 */
                    node3 = seg->node1;
                    node4 = seg->node2;

                    x3 = node3->x;
                    y3 = node3->y;
                    z3 = node3->z;

                    x4 = node4->x;
                    y4 = node4->y;
                    z4 = node4->z;

                    armID34 = GetArmID(node3, node4);

                    burg2[X] = node3->burgX[armID34];
                    burg2[Y] = node3->burgY[armID34];
                    burg2[Z] = node3->burgZ[armID34];

/*
 *                  Find coordinates of image of the third node closest to the
 *                  center of the cell containing the first node, and use the
 *                  image of the fourth node closest to that point.  We need
 *                  use the cell center as the base because all the Burgers
 *                  vectors of segments inside one cell (belonging to the same
 *                  PBC image) are summed together for 'remote' stress
 *                  calculations, so cannot separate them for 'local' stress
 *                  calculations.
 */
                    PBCPOSITION(param, xCenter, yCenter, zCenter,
                                &x3, &y3, &z3);
                    PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);

                    SegSegForce(home,
                                x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4,
                                burg1[X], burg1[Y], burg1[Z],
                                burg2[X], burg2[Y], burg2[Z],
                                a, MU, NU, 1, 0,
                                &fx1, &fy1, &fz1, &fx2, &fy2, &fz2,
                                &fx3, &fy3, &fz3, &fx4, &fy4, &fz4);

                    f1Thread[X] += fx1;
                    f1Thread[Y] += fy1;
                    f1Thread[Z] += fz1;

                    f2Thread[X] += fx2;
                    f2Thread[Y] += fy2;
                    f2Thread[Z] += fz2;

                }  /* end for (m = threadIterStart; ...) */

/*
 *              Add the sge/seg force from each thread to the segment
 *              forces.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_SMN_UPDATE_FORCE_B)
#endif
                {
                    VECTOR_ADD(f1Seg, f1Thread);
                    VECTOR_ADD(f2Seg, f2Thread);
                }

            }  /* end omp parallel section */

/*
 *          Add the locally calculated forces on the target segment to
 *          the total segment forces.  Remember that in the above calcs
 *          node1 was the node owning the segment, so we have to make sure
 *          that the forces at the segment endpoints are applied to
 *          the correct nodes...
 */
            if (targetOwnsSeg) {
                VECTOR_ADD(nodeForces->segForce[segID].f1, f1Seg);
                VECTOR_ADD(nodeForces->segForce[segID].f2, f2Seg);
            } else {
                VECTOR_ADD(nodeForces->segForce[segID].f1, f2Seg);
                VECTOR_ADD(nodeForces->segForce[segID].f2, f1Seg);
            }

        }  /* end for (segID = 0; segID < numNbrs; ...) */

        return;
}
