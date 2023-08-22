/**************************************************************************
 *
 *      Module:       SMNSupport.c
 *      Description:  Contains various support functions needed by
 *                    the multi-node splitting process.
 *
 *                    Note: some of the included functiosn are used
 *                          in both versions of the multi-node splitting
 *                          and some are unique to the parallel version.
 *
 *      Public functions:
 *          BuildGhostMultiNodeList()
 *          BuildLocalMultiNodeList()
 *          BuildRecvDomList()
 *          CopyForcesToNode()
 *          HasShortSeg()
 *          MergeGhostNode()
 *          SMNGetArmSets()
 *          SplitGhostNode()
 *
 *************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "SMN.h"
#include "QueueOps.h"

/*
 *      Use a predfined of list of the unique ways it is possible to
 *      split 2 or more arms from nodes of up to 15 arms.  We predefine
 *      these values since they are constant and there's no point in
 *      recomputing the values every time we need them.
 */
extern int POSSIBLE_SPLITS[16];


/**************************************************************************
 *
 *      Function:    SMNGetArmSets
 *
 *      Description: For each possible manner in which a node of the
 *                   specified number of segments can be split, create
 *                   an array of integer arrays of flags (1 per arm of
 *                   the node).
 *
 *                   These flags will be set to 1 or zero depending on
 *                   whether the corresponding arm is selected to be
 *                   split off the node in the corresponding split.
 *                   (actually, there's an additional integer at the
 *                   end of each array for the count of arms within
 *                   that list that have been selected)
 *
 *      Parameters:
 *          in:  numNbrs    Number of segments of the node to be split
 *          out: setCnt     Number of possible ways in which the node
 *                          may be split.
 *          out: armSetList Array of integer arrays as stated in the
 *                          function description above.  NOTE: The
 *                          caller is responsible for freeing the memory
 *                          allocated for armSetList.
 *
 *
 *************************************************************************/
void SMNGetArmSets(int numNbrs, int *setCnt, int ***armSetList)
{
        int j, splitCnt;
        int maxSets, maxSplitCnt;
        int level, totalSets;
        int *currSet, **armSets;


        if (numNbrs > MAX_NBRS) {
           Fatal("SplitLocalMultinodes found node with too many segs (%d)",
                  numNbrs);
        }

/*
 *      The possible number of ways to split n arms from m arms
 *      (where n >= 2) is given by:
 *
 *          m! / ((m-n)! * n!)   # dividing result by 2 if n*2==m
 *
 *      For efficiency, the values have been pre-computed for the
 *      maximum number of segments we allow for a node.
 */
        maxSets = POSSIBLE_SPLITS[numNbrs];
        maxSplitCnt = numNbrs >> 1;

        totalSets = 0;
        level = 0;

        currSet = (int *)malloc(sizeof(int) * maxSplitCnt);
        armSets = (int **)calloc(1, maxSets * sizeof(int *));

        for (j = 0; j < maxSets; j++) {
            armSets[j] = (int *)calloc(1, (numNbrs + 1) * sizeof(int));
        }

        for (splitCnt = 2; splitCnt <= maxSplitCnt; splitCnt++) {
            int maxStartArm, startArm;

            maxStartArm = (numNbrs - splitCnt) + 1;

            for (startArm = 0; startArm < maxStartArm; startArm++) {
                int numSets;

                currSet[0] = startArm;
                numSets = BuildSplitList(numNbrs, splitCnt,
                                         level, 0, currSet,
                                         &armSets[totalSets]);
                totalSets += numSets;

                if ((splitCnt << 1) == numNbrs) {
                    break;
                }
            }
        }

        if (totalSets != maxSets) {
                Fatal("%s: expected %d %s %d-node, but found %d",
                        "SplitMultiNode()", maxSets,
                        "split possibilities for", numNbrs, totalSets);
        }

        free(currSet);

        *setCnt = totalSets;
        *armSetList = armSets;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    CopyForcesToNode
 *
 *      Description: Overwite the node's segment forces with the
 *                   provided values which are then summed together
 *                   to reset the total nodal force.
 *
 *                   Note: This function updates both the node and all
 *                         of its neighbors!
 *
 *      Parameters:
 *          in:  node   Pointer to node struct.  Nodal and segment force
 *                      data in this struct will be updated inside this
 *                      function.
 *        
 *          in:  forces pointer to struct containing segment specific
 *                      forces for eahc segment of <node>
 *
 *-------------------------------------------------------------------------*/
void CopyForcesToNode(Home_t *home, Node_t *node, SMN_Force_t *forces)
{
        int    segID;

        node->fX = 0;
        node->fY = 0;
        node->fZ = 0;

/*
 *      Loop through the node's segments
 */
        for (segID = 0; segID < node->numNbrs; segID++) {
            int    nbrSegID;
            Node_t *nbrNode;

            nbrNode = GetNeighborNode(home, node, segID);
            nbrSegID = GetArmID(nbrNode, node);

/*
 *          Reset the segment's forces and increment the total
 *          nodal force.
 */
            node->armfx[segID] = forces->segForce[segID].f1[X];
            node->armfy[segID] = forces->segForce[segID].f1[Y];
            node->armfz[segID] = forces->segForce[segID].f1[Z];

            node->fX += node->armfx[segID];
            node->fY += node->armfy[segID];
            node->fZ += node->armfz[segID];

/*
 *          Reset the segment force at the other end of the segment
 *          as well, and update the neighboring node's total force.
 */
            nbrNode->fX = 0;
            nbrNode->fY = 0;
            nbrNode->fZ = 0;

            nbrNode->armfx[nbrSegID] = forces->segForce[segID].f2[X];
            nbrNode->armfy[nbrSegID] = forces->segForce[segID].f2[Y];
            nbrNode->armfz[nbrSegID] = forces->segForce[segID].f2[Z];

            for (nbrSegID = 0; nbrSegID < nbrNode->numNbrs; nbrSegID++) {
                nbrNode->fX += nbrNode->armfx[nbrSegID];
                nbrNode->fY += nbrNode->armfy[nbrSegID];
                nbrNode->fZ += nbrNode->armfz[nbrSegID];
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    BuildRecvDomList
 *
 *      Description: Build a list of all domains intersecting any cell
 *                   owning one of the multi-nodes in the lmnInfo array.
 *                   This list will encompass all domains from which the
 *                   local domain expects to receive force data during the
 *                   parallel multi-node splitting process.
 *
 *      Parameters:
 *          in:  lmnCount    Number of locally owned
 *          in:  lmnInfo     Array of <lmnCount> structs containing
 *                           information about each local multi-node to
 *                           be evaluated when doing parallel multi-node
 *                           splitting.
 *          in:  numRecvDoms Number of domains in <recvDomList>
 *          in:  recvDomList Array of domain ID's.
 *
 *-------------------------------------------------------------------------*/
void BuildRecvDomList(Home_t *home, int lmnCount, SMN_Info_t *lmnInfo,
                         int *numRecvDoms, int **recvDomList)
{
        int    i;
        int    domCount, *domList;

        domCount = 0;
        domList = (int *)NULL;

        for (i = 0; i < lmnCount; i++) {
            int    j, cellID, newSize;
            Cell_t *cell;

            if (lmnInfo[i].multiNode->numNbrs < 4) {
                continue;
            }

            cellID = lmnInfo[i].multiNode->cellIdx;

            if ((cell = LookupCell(home, cellID)) == (Cell_t *)NULL) {
                Fatal("BuildRecvDomList: cell not found");
            }
/*
 *          Add to the list, any remote domain intersecting the current cell
 */
            newSize = (domCount + cell->domCount) * sizeof(int);
            domList = (int *)realloc(domList, newSize);

            for (j = 0; j < cell->domCount; j++) {
                if (cell->domains[j] != home->myDomain) {
                    domList[domCount++] = cell->domains[j];
                }
            }

/*
 *          Now loop over all the neighboring cells and add to the
 *          list any remote domains intersecting those cells.
 */
            for (j = 0; j < cell->nbrCount; j++) {
                int    k, nbrCellID;
                Cell_t *nbrCell;

                nbrCellID = cell->nbrList[j];
                nbrCell = LookupCell(home, nbrCellID);

                if (nbrCell == (Cell_t *)NULL) {
                    printf("BuildRecvDomList: nbr cell not found\n");
                    continue;
                }

                if (nbrCell->baseIdx >= 0) {
                    nbrCellID = nbrCell->baseIdx;
                    nbrCell = LookupCell(home, nbrCellID);
                }

                if (nbrCell == (Cell_t *)NULL) {
                    printf("BuildRecvDomList: base nbr cell not found\n");
                    continue;
                }

                newSize = (domCount + nbrCell->domCount) * sizeof(int);
                domList = (int *)realloc(domList, newSize);

                for (k = 0; k < nbrCell->domCount; k++) {
                    if (nbrCell->domains[k] != home->myDomain) {
                        domList[domCount++] = nbrCell->domains[k];
                    }
                }
            }
        }

/*
 *      Domains intersecting multiple cells may have been added to
 *      the list multiple times, so remove any duplicate entries on
 *      the list before returning it to the caller.  We may have
 *      allocated a domain list but ended up with no domains on it.
 *      if so just free the memory now.
 */
        qsort(domList, domCount, sizeof(int), IntCompare);

        Uniq(domList, &domCount);

        if ((domCount == 0) && (domList != (int *)NULL)) {
            free(domList);
            domList = (int *)NULL;
        }

        *numRecvDoms = domCount;
        *recvDomList = domList;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    HasShortSeg
 *
 *      Description: Check all segments attached to the specified node
 *                   to see if any of the segments is considered 'short'
 *                   according to the provided criteria.
 *
 *      Parameters:
 *          in:  node          Pointer to node structure
 *          in:  shortSegLenSq Square of the length under which any segment
 *                             should be considered a short segment.
 *
 *      Returns: 1 if the node has any short segments
 *               0 in all other cases
 *
 *-------------------------------------------------------------------------*/
int HasShortSeg(Home_t *home, Node_t *node, real8 shortSegLenSq)
{
        int segID;

        for (segID = 0; segID < node->numNbrs; segID++) {
            Node_t *nbrNode;

            nbrNode = GetNeighborNode(home, node, segID);

            if (nbrNode != (Node_t *)NULL) {
                real8  segLenSq;
                real8  segVector[3];

                segVector[X] = node->x - nbrNode->x;
                segVector[Y] = node->y - nbrNode->y;
                segVector[Z] = node->z - nbrNode->z;

                ZImage(home->param, &segVector[X], &segVector[Y],
                       &segVector[Z]);

                segLenSq = DotProduct(segVector, segVector);

                if (segLenSq < shortSegLenSq) {
                    return(1);
                }
            }
        }

        return(0);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    BuildLocalMultiNodeList
 *
 *      Description: Check all locally owned nodes and return an array
 *                   of SMN_Info_t structs containing information about
 *                   each multi-node that should be considered for
 *                   splitting this cycle.
 *
 *      Parameters:
 *          out:  nodeCount  Number of multinodes returned in the
 *                           <nodeInfo> array.
 *          out:  nodeInfo   Array of <nodeCount> structs containing
 *                           information about each local multi-node to
 *                           be evaluated when doing parallel multi-node
 *                           splitting.
 *
 *-------------------------------------------------------------------------*/
void BuildLocalMultiNodeList(Home_t *home, int *nodeCount,
                             SMN_Info_t **nodeInfo)
{
        int         i;
        int         listSize, listCount;
        real8       shortSegSq;
        Param_t     *param;
        SMN_Info_t  *list;

        param = home->param;
        shortSegSq = MIN(5.0, param->minSeg * 0.1);

        listSize = 0;
        listCount = 0;
        list = (SMN_Info_t *)NULL;

/*
 *      Loop through all native nodes; any node with at least 4 arms
 *      is a candidate for being split.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            int    numSplits;
            Node_t *node;

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            if (node->numNbrs < 4) {
                continue;
            }

            numSplits = POSSIBLE_SPLITS[node->numNbrs];

/*
 *          If any of the segments attached to the multinode are less
 *          than <shortSeg> in length, don't split the node.  This is
 *          an attempt to prevent node splits that would leave extremely
 *          short segments which may oscillate at very high frequency.
 */
            if (HasShortSeg(home, node, shortSegSq)) {
                continue;
            }

            if (listCount >= listSize) {
                listSize += 50;
                list = (SMN_Info_t *)realloc(list,
                                             listSize * sizeof(SMN_Info_t));
            }

            list[listCount].status    = 0;
            list[listCount].numSplits = numSplits;
            list[listCount].multiNode = node;
            list[listCount].splitInfo =
                    (SMN_Split_t *)calloc(1, numSplits * sizeof(SMN_Split_t));

            listCount++;
        }

        *nodeCount = listCount;
        *nodeInfo = list;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    BuildGhostMultiNodeList
 *
 *      Description: Check all ghost nodes and return an array
 *                   of SMN_Info_t structs containing information about
 *                   each multi-node that this domain should considered
 *                   for splitting this cycle.
 *
 *      Parameters:
 *          out:  nodeCount  Number of multinodes returned in the
 *                           <nodeInfo> array.
 *          out:  nodeInfo   Array of <nodeCount> structs containing
 *                           information about each ghost multi-node to
 *                           be evaluated when doing parallel multi-node
 *                           splitting.
 *
 *-------------------------------------------------------------------------*/
void BuildGhostMultiNodeList(Home_t *home, int *nodeCount,
                             SMN_Info_t **nodeInfo)
{
        int         i, numGhostNodes;
        int         listSize, listCount;
        real8       shortSegSq;
        Param_t     *param;
        SMN_Info_t  *list;

        param = home->param;
        shortSegSq = MIN(5.0, param->minSeg * 0.1);

        listSize = 0;
        listCount = 0;
        list = (SMN_Info_t *)NULL;

/*
 *      Loop through all the ghost nodes looking for multinodes
 *      for which this domain should compute a portion of forces.
 */
        numGhostNodes = home->ghostNodeCount;

        for (i = 0; i < numGhostNodes; i++) {
            int    cellIndex;
            int    ghostIsNearby, ghostNodeCellID;
            Node_t *ghostNode;

            ghostNode = home->ghostNodeList[i];

            if (ghostNode->numNbrs < 4) {
                continue;
            }

/*
 *          If the ghost node is not nearby (not in a cell that is native
 *          to this domain, nor in a cell immediately neighboring one of the
 *          native cells) this domain should not do any portion of the
 *          force calcs when splitting the multinode.
 */
            ghostIsNearby = 0;
            ghostNodeCellID = ghostNode->cellIdx;

            for (cellIndex = 0; cellIndex < home->cellCount; cellIndex++) {
                if (ghostNodeCellID == home->cellList[cellIndex]) {
                    ghostIsNearby = 1;
                    break;
                }
            }

            if (!ghostIsNearby) {
                continue;
            }

/*
 *          If any of the segments attached to the multinode are less
 *          than <shortSeg> in length, don't split the node.  This is
 *          an attempt to prevent node splits that would leave extremely
 *          short segments which may oscillate at very high frequency.
 */
            if (HasShortSeg(home, ghostNode, shortSegSq)) {
                continue;
            }

            if (listCount >= listSize) {
                listSize += 50;
                list = (SMN_Info_t *)realloc(list,
                                             listSize * sizeof(SMN_Info_t));
            }

            list[listCount].status = 0;
            list[listCount].numSplits = 0;
            list[listCount].multiNode = ghostNode;
            list[listCount].splitInfo = (SMN_Split_t *)NULL;
            listCount++;
        }

        *nodeCount = listCount;
        *nodeInfo = list;

        return;
}


/*---------------------------------------------------------------------------
 *
 *	Function:	SplitGhostNode
 *	Description:	Create a new ghost node and transfer the specified
 *                      set of connections from an existing ghost node to
 *                      the new node.  If necessary, a new link will also
 *                      be created between the existing node and the new
 *                      node.
 * 
 *                      WARNING!  This should only be used with
 *                                SplitMultiNodeParallel() for testing
 *                                ways of splitting multinodes in remote
 *                                domains.
 *
 *	Arguments:
 *		node            Pointer to the node to be split
 *              pos1            coordinates at which splitNode1 will be
 *                              left after the split 
 *              pos2            coordinates at which splitNode2 will be
 *                              left after the split
 *              vel1            velocity assigned to splitNode1 after the split
 *              vel2            velocity assigned to splitNode2 after the split
 *		armCount	number of arms of the original node selected
 *                              to be split off.
 *		armList		pointer to array of integers indicating the
 *				arms of existing node that are to be split off
 *              splitNode1      ptr to ptr to node to which all unselected arms
 *                              of the original node will be attached after the
 *                              split.  Returned to caller.
 *              splitNode2      ptr to ptr to node to which all selected arms
 *                              of the original node will be attached after the
 *                              after the split.  Returned to caller.
 *              flags           Bit field with additional processing flags
 *
 *	Returns:		1 if the split was successful
 *                              0 in all other cases
 *
 *-------------------------------------------------------------------------*/
int SplitGhostNode(Home_t *home, int opClass, Node_t *node, real8 *pos1,
              real8 *pos2, real8 *vel1, real8 *vel2, int armCount,
              int *armList, Node_t **splitNode1,
              Node_t **splitNode2, int flags)
{
	int	i, j, k, tarm, *tarmList, *tarmList2;
        int     newIndex;
        int     globalOp;
        int     owningDomain;
	real8	bx, by, bz;
	real8	cx, cy, cz;
	real8	nx, ny, nz;
        real8   segPlaneX=0.0, segPlaneY=0.0, segPlaneZ=0.0;
	real8	eps;
	real8   ftmp[3];
        real8   burg[3], vel[3], dirVec[3], glidePlane[3];
        Cell_t  *owningCell;
	Node_t	*newNode, *nbrNode;
	Param_t	*param;
        int     ownsSelectedArms, switchedArms;


	param      = home->param;

        owningDomain = node->myTag.domainID;

        *splitNode1 = (Node_t *)NULL;
        *splitNode2 = (Node_t *)NULL;

        switchedArms = 0;
        globalOp = 0;
        eps = 1.0e-12;

/*
 *      If we are only moving a single segment from the original node to the
 *      new node, it (likely) means we're simply bisecting a segment during
 *      MeshRefine().  In that case, preserve the segment's glide plane
 *      so when we add a new segment between the original node and new node
 *      we can have the new segment properly inherit the glide plane.
 */
        if (armCount == 1) {
            segPlaneX = node->nx[armList[0]];
            segPlaneY = node->ny[armList[0]];
            segPlaneZ = node->nz[armList[0]];
        }

/*
 *	The list of arms MUST be sorted in descending order (by index) to 
 *	avoid issues related to arm renumbering as arms are removed from the
 *	existing node.  Make a temporary copy of the arm list and sort
 *	it appropriately here, just in case the caller didn't.
 */
	tarmList = (int *)calloc(1, sizeof(int) * armCount);

	for (i = 0; i < armCount; i++) tarmList[i] = armList[i];

	for (i = 0; i < armCount-1; i++) {
		for (j = 0; j < armCount - (i+1); j++) {
			if (tarmList[j] < tarmList[j+1]) {
				tarm = tarmList[j+1];
				tarmList[j+1] = tarmList[j];
				tarmList[j] = tarm;
			}
		}
	}

/*
 *      All arms to be split from the original and moved to the new node
 *      MUST be owned by the domain containing the node to be split
 */
        ownsSelectedArms = 1;

        for (i = 0; i < armCount; i++) {
            ownsSelectedArms &= DomainOwnsSeg(home, opClass, owningDomain,
                                              &node->nbrTag[tarmList[i]]);
        }

/*
 *      We cannot do a split that would move the endpoint of a segment
 *      not owned by the splitting node's domain, so if that domain does
 *      not own all arms that are to be moved to the new node,  check
 *      if the arms to remain with <node> are all owned by node's
 *      domain.  If they are, we can still do the split by switching
 *      the selected segments with the unselected segments.
 */
        if (!ownsSelectedArms) {

            j = 0;
            k = 0;

	    tarmList2 = (int *)calloc(1, (node->numNbrs-armCount)*sizeof(int));

            for (i = node->numNbrs-1; i >= 0; i--) {
                if ((j >= armCount) || (tarmList[j] != i))
                    tarmList2[k++] = i;
                else
                    j++;
            }

            free(tarmList);
            tarmList = tarmList2;
            tarmList2 = (int *)NULL;
            armCount = node->numNbrs - armCount;

            ownsSelectedArms = 1;

            for (i = 0; i < armCount; i++) {
                ownsSelectedArms &= DomainOwnsSeg(home, opClass, owningDomain,
                                                  &node->nbrTag[tarmList[i]]);
            }

            if (!ownsSelectedArms) {
                free(tarmList);
                return(SPLIT_FAILED);
            }

            switchedArms = 1;
        }

/*
 *      We temporarily need a new ghost node, but don't know what ghost
 *      nodes are available in the remote domain, so always add the new
 *      one as the last one.
 */
        newIndex = home->remoteDomainKeys[node->myTag.domainID]->maxTagIndex;

/*
 *	Add a new node.  Warning:  GetNewGhostNode() pops a node off
 *      the queue, but does not free the arm related arrays.  It does
 *      this on the assumption that the number of arms previously
 *      assigned to the node is likely to be the same that will
 *      be reattached to the node and that the code will overwrite
 *      the existing arm data.  The following code, however, assumes
 *	the new node has zero arms initially, so we MUST free those arrays
 *	and reset the arm count!
 */
	newNode = GetNewGhostNode(home, node->myTag.domainID, newIndex);
	FreeNodeArms(newNode);

	newNode->native = 1;
	SET_CONSTRAINTS(newNode->constraint, UNCONSTRAINED);

	newNode->oldvX = 0.0;
	newNode->oldvY = 0.0;
	newNode->oldvZ = 0.0;

#ifdef _ARLFEM
/*
 *      If necessary copy the surface constraint and properties
 *      from the original node to the new node.
 */
        if ((flags & SPLIT_DUP_SURFACE_PROP) &&
            HAS_ANY_OF_CONSTRAINTS(node->constraint, SURFACE_NODE)) {

            ADD_CONSTRAINTS(newNode->constraint, SURFACE_NODE);
            VECTOR_COPY(newNode->surfaceNorm, node->surfaceNorm);
        }
#endif

/*
 *      Reset the velocities and positions for the two nodes involved
 *      in the split.  If we split off the arms selected by the caller
 *      splitNode1 will be mapped to the original node. Otherwise, it
 *      will be mapped to the new node.
 */
        if (switchedArms) {
            *splitNode1 = newNode;
            *splitNode2 = node;
        } else {
            *splitNode1 = node;
            *splitNode2 = newNode;
        }

        (*splitNode1)->x = pos1[X];
        (*splitNode1)->y = pos1[Y];
        (*splitNode1)->z = pos1[Z];

        (*splitNode1)->vX = vel1[X];
        (*splitNode1)->vY = vel1[Y];
        (*splitNode1)->vZ = vel1[Z];

        (*splitNode2)->x = pos2[X];
        (*splitNode2)->y = pos2[Y];
        (*splitNode2)->z = pos2[Z];

        (*splitNode2)->vX = vel2[X];
        (*splitNode2)->vY = vel2[Y];
        (*splitNode2)->vZ = vel2[Z];

	newNode->oldx = newNode->x;
	newNode->oldy = newNode->y;
	newNode->oldz = newNode->z;

/*
 *      If the node being split was pinned in any dimensions, the
 *      node left at position <pos1> must be pinned in the same way, and
 *      the other node must be explicitly unpinned.
 *
 *      NOTE: For now, a node that is pinned in ANY dimension is treated
 *            here as if it is pinned in ALL dimensions.
 */
        if (HAS_ANY_OF_CONSTRAINTS(node->constraint, PINNED_NODE)) {
            ADD_CONSTRAINTS((*splitNode1)->constraint,
                            node->constraint & PINNED_NODE);
            REMOVE_CONSTRAINTS((*splitNode2)->constraint, PINNED_NODE);
        }

#ifdef _ARLFEM
/*
 *      If the node being split is a surface node, then the
 *      node left at <pos1> must retain the surface properties.
 *      Additionally, if the surface properties of the original
 *      node are NOT to be duplicated, we must explicitly zero
 *      out surface constraint in the node at <pos2>.
 */
        if (HAS_ANY_OF_CONSTRAINTS(node->constraint, SURFACE_NODE)) {

            ADD_CONSTRAINTS((*splitNode1)->constraint, SURFACE_NODE);
            VECTOR_COPY((*splitNode1)->surfaceNorm, node->surfaceNorm);

            if ((flags & SPLIT_DUP_SURFACE_PROP) == 0) {
                REMOVE_CONSTRAINTS((*splitNode2)->constraint, SURFACE_NODE);
            }
        }
#endif

/*
 *	For each arm of the original node in armList, get the neighbor node,
 *	change the nbrNode->origNode linkage to nbrNode->newNode, and
 *	remove the origNode->nbrNode linkage.  Note: the DEL_SEG_NONE flag
 *	is passed to ChangeArmBurg() to prevent the function from treating
 *	the segments as annihilated since they are simply having an endpoint
 *	shifted.
 */
	for (i = 0; i < armCount; i++) {

		nbrNode = GetNeighborNode(home, node, tarmList[i]);

                if (nbrNode == (Node_t *)NULL) {
                    return(SPLIT_FAILED);
                }

		ftmp[X] = node->armfx[tarmList[i]];
		ftmp[Y] = node->armfy[tarmList[i]];
		ftmp[Z] = node->armfz[tarmList[i]];

		SubtractSegForce(node, nbrNode);

		GetBurgersVectorNormal(node, nbrNode, &bx, &by, &bz, &nx, &ny, &nz);
		ChangeConnection(home, nbrNode, &node->myTag,
				 &newNode->myTag, globalOp);
		InsertArm(home, newNode, &node->nbrTag[tarmList[i]],
			  bx, by, bz, nx, ny, nz, globalOp);
		ChangeArmBurg(home, node, &nbrNode->myTag, 0, 0, 0, 0, 0, 0,
			      globalOp, DEL_SEG_NONE);
		ResetSegForces(home, newNode, &nbrNode->myTag, ftmp[X],
			       ftmp[Y], ftmp[Z], globalOp);
	}

	free(tarmList);


/*
 *      Need to explicitly assign the new node to the same cell as
 *      the original node.  Can't use AssignNodeToCell() since the
 *      node may be outside any cells native to this local domain
 *      and AssignNodeToCell() will not assign cell a node to a
 *      non-native cell.
 */
        owningCell = LookupCell(home, node->cellIdx);

        newNode->cellIdx = node->cellIdx;
        newNode->nextInCell = owningCell->nodeQ;

        owningCell->nodeQ = newNode;
        owningCell->nodeCount++;

/*
 *	If necessary, create a link between the original node and new node
 *	in order to conserve burgers vector.  First sum up the burgers vectors
 *	for the arms of the new node...
 *
 *      One exception, though.  When we're dealing with free surfaces, if
 *      *both* nodes are on the surface, we do not want to add a new segment
 *      between the nodes...
 */
	bx = 0.0;
	by = 0.0;
	bz = 0.0;

        if (HAS_NO_CONSTRAINTS((*splitNode1)->constraint) ||
            HAS_NO_CONSTRAINTS((*splitNode2)->constraint)) {

            for (i = 0; i < (*splitNode2)->numNbrs; i++) {
                bx += (*splitNode2)->burgX[i];
                by += (*splitNode2)->burgY[i];
                bz += (*splitNode2)->burgZ[i];
            }
        }

/*
 *      Under two circumstances we need to add an additional link between
 *      the two nodes:
 *      1) if the sum of burgers vectors on the new node is non-zero
 *      2) if the sum of burgers vectors on the new node is zero BUT
 *         the node is 'pinned';  Scenario is this:
 *         -  the single arm of a pinned node is being split
 *         -  the 'pinned' node is repositioned and unpinned
 *         -  the new node is positioned at the location of
 *            the original node and pinned.
 *         -  the original node maintains its single connection
 *         
 *         At this point, the new node is pinned, but has no arms
 *         hence the sum of its burgers vectors is zero, so no
 *         link would be created between the original node and the
 *         new node.  This results in the pinned node being orphaned
 *         and deleted, and the other node is no longer attached to
 *         a pinned node and this is bad.
 */
	if (fabs(bx*bx + by*by + bz*bz) > 0.001) {

/*
 *              If this was just a call to bisect a segment,
 *              the glide plane should just be inherited from
 *              the bisected segment.  Otherwise, we need to
 *              calculate the glide plane of the new link.
 */
                if (((*splitNode1)->numNbrs == 1) ||
                    ((*splitNode2)->numNbrs == 1)) {
                    cx = segPlaneX;
                    cy = segPlaneY;
                    cz = segPlaneZ;
                } else {
                    dirVec[X] = (*splitNode1)->x - (*splitNode2)->x;
                    dirVec[Y] = (*splitNode1)->y - (*splitNode2)->y;
                    dirVec[Z] = (*splitNode1)->z - (*splitNode2)->z;

                    ZImage(param, &dirVec[X], &dirVec[Y], &dirVec[Z]);

                    burg[X] = bx;
                    burg[Y] = by;
                    burg[Z] = bz;

                    NormalizeVec(dirVec);
                    NormalizeVec(burg);

                    FindPreciseGlidePlane(home, burg, dirVec, glidePlane,
                                          home->param->allowFuzzyGlidePlanes);

                    if (DotProduct(glidePlane,glidePlane) < eps) {
                        vel[X] = ((*splitNode1)->vX + (*splitNode2)->vX) * 0.5;
                        vel[Y] = ((*splitNode1)->vY + (*splitNode2)->vY) * 0.5;
                        vel[Z] = ((*splitNode1)->vZ + (*splitNode2)->vZ) * 0.5;

                        NormalizeVec(vel);
                        FindPreciseGlidePlane(home, vel, dirVec, glidePlane,
                                          home->param->allowFuzzyGlidePlanes);

                        if (DotProduct(glidePlane,glidePlane) < eps) {
                            PickScrewGlidePlane(home, burg, glidePlane);
                        }
                    }

                    NormalizeVec(glidePlane);

                    cx = glidePlane[X];
                    cy = glidePlane[Y];
                    cz = glidePlane[Z];
                }

/*
 *		Add links from the original node to the new node
 *		and vice-versa.  Forces on this new segment will be
 *		added outside this function.
 */
		InsertArm(home, (*splitNode1), &(*splitNode2)->myTag,
                          bx, by, bz, cx, cy, cz, globalOp);

		InsertArm(home, (*splitNode2), &(*splitNode1)->myTag,
                          -bx, -by, -bz, cx, cy, cz, globalOp);

	}

	return(SPLIT_SUCCESS);
}


/*---------------------------------------------------------------------------
 *
 *	Function:	MergeGhostNode
 *	Description:	This function 'merges' two ghost nodes by moving
 *                      all arms of <deadNode> to <targetNode>, and then
 *			completely removing <deadNode>.
 * 
 *                      WARNING!  This should only be used with
 *                                SplitMultiNodeParallel() when testing
 *                                ways of splitting multinodes in remote
 *                                domains, and as such does not handle
 *                                all situations the normal MergeNode()
 *                                function does.
 *
 *	Arguments:
 *              opClass         Flag indicating the class of topological
 *                              operation invoking this function.
 *                              initiated.  Valid types are:
 *
 *                                  OPCLASS_SEPARATION
 *                                  OPCLASS_COLLISION
 *                                  OPCLASS_REMESH
 *
 *              node1           pointer to first node to be merged
 *              node2           pointer to second node to be merged
 *              position        coordinates (x,y,z) at which final merged
 *                              node is to be placed.
 *              mergedNode      pointer to location in which to return
 *                              pointer to the node resulting from the merge.
 *                              A NULL pointer will be returned if the
 *                              the merge fails.
 *              status          pointer to location in which to return
 *                              a completion status to the caller.  Valid
 *                              statuses are the following, where all statuses
 *                              indicating success may be logically OR'ed
 *                              
 *                                  MERGE_SUCCESS
 *                                  MERGE_NO_REPOSITION
 *                                  MERGE_NODE_ORPHANED
 *                                  MERGE_NOT_PERMITTED
 *                                  MERGE_DOUBLE_LINK
 *
 *-------------------------------------------------------------------------*/
void MergeGhostNode(Home_t *home, int opClass, Node_t *node1, Node_t *node2,
               real8 *position, Node_t **mergedNode, int *status)
{
        int	i;
        int     targetDomain;
        int     armID; 
        int     node1Deletable, node2Deletable;
        int     globalOp = 0;
        real8	bx, by, bz;
        real8	nx, ny, nz;
        real8   ftmp[3];
        Node_t	*nbrNode, *targetNode, *deadNode;
        Tag_t	tmpTag;
        RemoteDomain_t *remDom;

        targetDomain = node1->myTag.domainID;

        targetNode  = node1;
        deadNode    = node1;
        *mergedNode = (Node_t *)NULL;
        *status     = MERGE_NOT_PERMITTED;

/*
 *      Check if we are permitted to delete node1.  In order
 *      for a node to be deletable, it must satisfy the following:
 *
 *      1) Node may not be a 'pinned' node 
 *      2) Node must be owned by the target's domain
 *      3) All arms of the node must be owned by the target's domain
 *      4) Node may not be exempt from (additional) collisions this cycle
 *      5) Node may not be exempt from (additional) remeshes this cycle
 *      6) If node is a surface node, *both* nodes must be on the surface
 *         for it to be deletable.  (May need to add restrictions
 *         that surface nodes be on the same surface)
 *
 *      NOTE: For now, a node that is pinned in ANY dimension is treated
 *            here as if it is pinned in ALL dimensions
 */
        node1Deletable  = (node1->myTag.domainID == targetDomain); 

        node1Deletable &= HAS_NONE_OF_CONSTRAINTS(node1->constraint,
                                                  PINNED_NODE);

        node1Deletable &= (HAS_NO_CONSTRAINTS(node1->constraint) ||
                           HAS_ANY_OF_CONSTRAINTS(node2->constraint,
                                                  SURFACE_NODE));

        node1Deletable &= ((node1->flags & NO_COLLISIONS) == 0);
        node1Deletable &= ((node1->flags & NO_MESH_COARSEN) == 0);

        for (armID = 0; armID < node1->numNbrs; armID++) {
            node1Deletable &= DomainOwnsSeg(home, opClass, targetDomain,
                                            &node1->nbrTag[armID]);
        }

/*
 *      If we can delete node1, use node2 as the target, but if node1
 *      cannot be removed, use node2 if we can. If node2 cannot
 *      be deleted either, return a failure to the caller.
 */
        if (node1Deletable) {
            targetNode = node2;
            deadNode = node1;
        } else {
            node2Deletable  = (node2->myTag.domainID == targetDomain); 
            node2Deletable &= HAS_NONE_OF_CONSTRAINTS(node2->constraint,
                                                      PINNED_NODE);
            node2Deletable &= (HAS_NO_CONSTRAINTS(node2->constraint) ||
                               HAS_ANY_OF_CONSTRAINTS(node1->constraint,
                                                  SURFACE_NODE));
            node2Deletable &= ((node2->flags & NO_COLLISIONS) == 0);
            node2Deletable &= ((node2->flags & NO_MESH_COARSEN) == 0);

            for (armID = 0; armID < node2->numNbrs; armID++) {
                node2Deletable &= DomainOwnsSeg(home, opClass, targetDomain,
                                                &node2->nbrTag[armID]);
            }

            if (!node2Deletable) {
                *status = MERGE_NOT_PERMITTED;
                *mergedNode = (Node_t *)NULL;
                return;
            }

            targetNode = node1;
            deadNode = node2;
        }

/*
 *      We've made it past all the checks that would have prevented the
 *      merge operation from being performed.  The target and dead nodes
 *      have been selected, so if the target node is local to this
 *      domain, go ahead and reposition it.  If we can't reposition the
 *      node, the merge will still be done, but be sure to set the
 *      status flag for the caller indicating the node has not been
 *      repositioned.
 */
        *mergedNode = targetNode;
        *status = MERGE_SUCCESS;

        if ((targetNode->myTag.domainID == targetDomain) &&
            (HAS_NO_CONSTRAINTS(targetNode->constraint) ||
             (HAS_ANY_OF_CONSTRAINTS(targetNode->constraint, SURFACE_NODE) &&
              HAS_ANY_OF_CONSTRAINTS(deadNode->constraint, SURFACE_NODE)))) {

            targetNode->x = position[X];
            targetNode->y = position[Y];
            targetNode->z = position[Z];

        } else {
            *status |= MERGE_NO_REPOSITION;
        }

/*
 *      If there are any connections from the target/destination node
 *      to the dead node, use ChangeArmBurg() to get rid of those
 *      connections, then move all connections from the dead node to
 *      the target node, and add a new connection from the target node
 *      to each of the dead node's neighbors.
 */
        ChangeArmBurg(home, targetNode, &deadNode->myTag, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);
        ChangeArmBurg(home, deadNode, &targetNode->myTag, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);

/*
 *      Move all connections from the dead node to the target node
 *      and add a new connection from the target node to each of the
 *      dead node's neighbors.
 */
        for (i = (deadNode->numNbrs - 1); i >= 0; i--) {
            tmpTag.domainID = deadNode->nbrTag[i].domainID;
            tmpTag.index    = deadNode->nbrTag[i].index;
            nbrNode = GetNeighborNode(home, deadNode, i);

            if (nbrNode == (Node_t *)NULL) {
                Fatal("Neighbor not found at %s line %d\n",
                       __FILE__, __LINE__);
            }

            ftmp[X] = deadNode->armfx[i];
            ftmp[Y] = deadNode->armfy[i];
            ftmp[Z] = deadNode->armfz[i];

            GetBurgersVectorNormal(deadNode, nbrNode, &bx, &by, &bz, &nx, &ny, &nz);
            ChangeConnection(home, nbrNode, &deadNode->myTag,
                             &targetNode->myTag, globalOp);
            ChangeArmBurg(home, deadNode, &nbrNode->myTag, 0.0, 0.0, 0.0,
                          0.0, 0.0, 0.0, globalOp, DEL_SEG_NONE);
            InsertArm(home, targetNode, &tmpTag, bx, by, bz,
                      nx, ny, nz, globalOp);

            ResetSegForces(home, targetNode, &tmpTag, ftmp[X], ftmp[Y],
                           ftmp[Z], globalOp);
        }

/*
 *      Any constraints that may have been associated with the node
 *      being removed should be transferred to the target node.
 */
        ADD_CONSTRAINTS(targetNode->constraint, deadNode->constraint);

/*
 *      The <deadNode> should have no more arms, so we can safely
 *      remove it now.  Since it's a ghost node and we really want
 *      to get rid of it now, we can't use RemoveNode(), but have
 *      to explicitly take care of a bunch of stuff.
 *
 *        - Remove node from cell Q
 *        - Remove node from ghost node Q
 *        - Remove node from remote domain nodeKeys
 *        - Decrement remote domain maxTagIndex
 *        - Add node to free node Q
 *
 *      Based on the manner in which we temporarily added a ghost
 *      node to test multinode splits in remote domains, the deleted
 *      node should always occupy the first spot on the ghost node
 *      queue (and final entry in the ghost node list) and be
 *      the final entry in the remote domain's nodeKeys array...
 */
        RemoveNodeFromCellQ(home, deadNode);

        home->ghostNodeQ = deadNode->next;
        home->ghostNodeCount -= 1;
        home->ghostNodeList[home->ghostNodeCount] = (Node_t *)NULL;

        remDom = home->remoteDomainKeys[targetDomain];
        remDom->nodeKeys[deadNode->myTag.index] = (Node_t *)NULL;
        remDom->maxTagIndex -= 1;

        PushFreeNodeQ(home, deadNode);
        
        return;
}
