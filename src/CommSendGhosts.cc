/***************************************************************************
 *
 *      Module:      CommSendGhosts
 *      Description: Contains functions necessary for communicating
 *                   new nodal data between domains that own the
 *                   nodes and domains that have the nodes exported
 *                   as ghost nodes.
 *
 *                   This function is the prototypical communication
 *                   mechanism for the code.  Messages are sent using
 *                   asynchronous communications in the following steps:
 *
 *                      1) issue receives of incoming buffer lengths
 *                      2) pack data into output buffers
 *                      3) send lengths of output buffers
 *                      4) wait for length sends and receives to complete
 *                      5) allocate input buffers
 *                      6) issue receives for incoming data
 *                      7) issue sends of outgoing data
 *                      8) wait for sends and receives to complete
 *                      9) unpack data into local structs
 *                     10) free buffers
 *
 *      Includes functions:
 *
 *          CommPackGhosts()
 *          CommSendGhosts()
 *          CommUnpackGhosts()
 *          CompressNodeKeyArray()
 *          ReconcileCompressedNodeKeys()
 *          
 *
 **************************************************************************/
#include "mpi_portability.h"

#include "Home.h"
#include "RemoteDomain.h"
#include "Cell.h"
#include "Node.h"
#include "Comm.h"
#include "QueueOps.h"

/*
 *      Define a structure to hold some per-domain info regarding the
 *      number of entities being imported from remote domains and the indices
 *      of those entities in some temporary local arrays.
 */
typedef struct {
    int nodeCount;
    int nodeIndex;
    int inclusionIndex;
} DomEntityInfo_t;


/*-------------------------------------------------------------------------
 *
 *      Function:    CompressNodeKeyArray
 *
 *      Description: The local nodeKey array may become sparsely populated
 *                   due to node migration and/or node deletion.  For code
 *                   segments which loop through the nodekeys array, this
 *                   means many unnecesary loop iterations which is
 *                   particularly bad when threading is enabled.  As such,
 *                   this function is used to compress the array so there
 *                   are no empty slots in the array.
 *
 *------------------------------------------------------------------------*/
static void CompressNodeKeyArray(Home_t *home)
{
        int    i, lastNodeIndex, nextAvailableIndex, newIndex;
        Tag_t  oldTag;
        Node_t *node;

/*
 *      Loop through the nodeKeys in reverse, moving the node with the
 *      highest index into the free spot with the lowest index.  This
 *      insures we don't move a node from spot B into spot A and then
 *      re-use spot B for another node.  This keeps it simple and
 *      makes it easier to reconcile the changes everywhere.
 */
        lastNodeIndex = home->newNodeKeyPtr - 1;

        for (i = lastNodeIndex; i > 0; i--) {

            if (home->nodeKeys[i] == (Node_t *)NULL) {
                continue;
            }

            nextAvailableIndex = HeapPeek(home->recycledNodeHeap,
                                          home->recycledNodeHeapEnts);

/*
 *          If the there are no more available tags on the recycled node heap,
 *          or the next available tag has a larger index than the current node
 *          it means the list cannot be compressed any further, so break
 *          out of the loop and sort the tagMap list.
 */
            if ((nextAvailableIndex < 0) ||
                (nextAvailableIndex > i)) {
                break;
            }

/*
 *          Move the current node to the lowest available slot in the
 *          nodeKeys array and free up the current slot.
 */
            newIndex = GetFreeNodeTag(home);
            RecycleNodeTag(home, i);

            home->nodeKeys[newIndex] = home->nodeKeys[i];
            home->nodeKeys[i] = (Node_t *)NULL;

            node = home->nodeKeys[newIndex];
            node->myTag.index = newIndex;

/*
 *          Add a mapping from the old tag to the new tag so we
 *          can reconcile the change everywhere later.
 */
            oldTag.domainID = home->myDomain;
            oldTag.index    = i;

            AddTagMapping(home, &oldTag, &node->myTag);

            home->newNodeKeyPtr -= 1;
        }

/*
 *      Sort the tag map based on the 'new' tags for now.
 */
        if (home->tagMapEnts > 0) {
            qsort(home->tagMap, home->tagMapEnts, sizeof(TagMap_t),
                  TagMapCompareNew);
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ReconcileCompressedNodeKeys
 *
 *      Description: After the local nodeKeys array has been compressed
 *                   and the local domain has received all ghost nodes
 *                   along with mappings between old and new tags for any
 *                   of those ghost nodes, this function is called to 
 *                   loop through all the local and ghost nodes and 
 *                   update their list of neighbors using the list of
 *                   all mappings between old and new node tags.
 *
 *                   IMPORTANT!  This function asumes the tagMap array
 *                               has already been sorted (using the 'old'
 *                               tags as the sort key)
 *
 *------------------------------------------------------------------------*/
static void ReconcileCompressedNodeKeys(Home_t *home)
{
        int      i;

/*
 *      Loop through every nbrTag of each local node.  Check the
 *      tag mapping list to see if the neighbor nodes have been
 *      retagged.  If so, update the local node's corresponding nbrTag
 *
 *      NOTE: The node tags for all local and ghost nodes were updated
 *            before the ghost node comm was done.  That's why we only
 *            have to reconcile the nodes' nbrTag arrays here.
 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            int      j, numNbrs;
            Node_t   *node;
            TagMap_t key;
            TagMap_t *mapping;

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            numNbrs = node->numNbrs;

            for (j = 0; j < numNbrs; j++) {

                key.oldTag.domainID = node->nbrTag[j].domainID;
                key.oldTag.index    = node->nbrTag[j].index;

                mapping = (TagMap_t *)bsearch(&key, home->tagMap,
                                              home->tagMapEnts,
                                              sizeof(TagMap_t),
                                              TagMapCompareOld);

/*
 *              If there is a mapping to a new tag for this nbr node,
 *              reset the tag index (domain ID will not have changed
 *              during a comprerssion of the nodeKeys array).
 */
                if (mapping != (TagMap_t *)NULL) {
                    node->nbrTag[j].index = mapping->newTag.index;
                }
            }

        }  /* end "omp parallel for" */

/*
 *      Now we have to loop through all the ghost nodes and do the
 *      same...
 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->ghostNodeCount; i++) {
            int      j, numNbrs;
            Node_t   *node;
            TagMap_t *mapping;

            node = home->ghostNodeList[i];

            numNbrs = node->numNbrs;

            for (j = 0; j < numNbrs; j++) {

                mapping = (TagMap_t *)bsearch(&node->nbrTag[j], home->tagMap,
                                              home->tagMapEnts,
                                              sizeof(TagMap_t),
                                              TagMapCompareOld);

/*
 *              If there is a mapping to a new tag for this nbr node,
 *              reset the tag index (domain ID will not have changed
 *              during a comprerssion of the nodeKeys array).
 */
                if (mapping != (TagMap_t *)NULL) {
                    node->nbrTag[j].index = mapping->newTag.index;
                }
            }

        }  /* end "omp parallel for" */

/*
 *      We're done with the tag mapping list now, so clear it out
 */
        if (home->tagMapEnts > 0) {
            free(home->tagMap);
            home->tagMap = (TagMap_t *)NULL;
            home->tagMapSize = 0;
            home->tagMapEnts = 0;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    CommPackGhosts
 *
 *      Description: For each neighbor domain, pack into a send buffer
 *                   all the entities in any cell that borders the neighbor
 *
 *------------------------------------------------------------------------*/
static void CommPackGhosts(Home_t *home) 
{
#ifdef PARALLEL

/*
 *      We need to loop through all the remote domains packing a buffer
 *      for each with the nodal data for all local entities in cells
 *      that are exported to those domains.  Since these loop iterations
 *      will not take the same amount of processing time, don't use the
 *      OpenMP "parallel for" construct which assumes that loop iterations
 *      do take the same amout of time.
 */
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int i, threadID, threadIterStart, threadIterEnd;

            GetThreadIterationIndices(home->remoteDomainCount, &threadID,
                                      &threadIterStart, &threadIterEnd);

/*
 *          Have each thread loop through a portion of the remote domains.
 */
            for (i = threadIterStart; i < threadIterEnd; i++) {
                int    j, valCount;
                int    domainIndex, cellIndex;
                int    armCount, totNodeCount, totInclusionCount;
                int    bufIndex;
                int    outBufLen;
                real8  *outBuf;
                Cell_t *cell;
                Node_t *node;
                RemoteDomain_t *remDom;

                domainIndex = home->remoteDomains[i];
                remDom = home->remoteDomainKeys[domainIndex];

/*
 *              Get total counts of the entities in all cells being
 *              exported to this remote domain.
 */
                armCount = 0;
                totNodeCount = 0;
                totInclusionCount = 0;

                for (j = 0 ; j < remDom->numExpCells ; j++) {

                    cellIndex = remDom->expCells[j];
                    cell = LookupCell(home, cellIndex);

                    totNodeCount += cell->nodeCount;
#ifdef ESHELBY
                    totInclusionCount += cell->inclusionCount;
#endif

                    node = cell->nodeQ;

                    while (node) {
                        armCount += node->numNbrs;
                        node = node->nextInCell;
                    }
                }

/*
 *              Calculate the buffer size needed to hold the data for
 *              this remote domain and allocate it.
 */
                valCount = FLTS_PER_GHOST_NODE * totNodeCount +
                           FLTS_PER_GHOST_ARM * armCount +
                           FLTS_PER_GHOST_CELL * remDom->numExpCells +
                           FLTS_PER_GHOST_INCLUSION * totInclusionCount +
                           EXTRA_GHOST_FLTS;

                outBufLen = (valCount * sizeof(real8));
                outBuf = (real8 *)malloc(outBufLen);

                bufIndex = 0;

/*
 *              Send the remote domain the numbers of entities it will receive, 
 *              the maximum possible tag.index it will receive, and the
 *              number of cells in the message
 */
                bufIndex++;    /* reserve space for node count */
                bufIndex++;    /* reserve space for inclusion count */

                outBuf[bufIndex++] = home->newNodeKeyPtr;
                outBuf[bufIndex++] = remDom->numExpCells;

                for (j = 0; j < remDom->numExpCells; j++) {
                    int  k, nodesPacked;

                    cellIndex = remDom->expCells[j];
                    cell = LookupCell(home, cellIndex);

/*
 *                  Supply the cell Index and the number of entities to follow
 *                  for this cell
 */
                    outBuf[bufIndex++] = cellIndex;
                    outBuf[bufIndex++] = cell->nodeCount;
                    outBuf[bufIndex++] = cell->inclusionCount;

/*
 *                  loop through all the nodes on the cell's node queue
 */
                    nodesPacked = 0;
                    node = cell->nodeQ;

                    while (node) {
                        int       nbrCount;
                        TagMap_t  key;
                        TagMap_t  *mapping;

                        nbrCount = node->numNbrs;

/*
 *                      The node's domainID is known implicitly, don't
 *                      need to send it
 */
                        outBuf[bufIndex++] = node->myTag.index;
                        outBuf[bufIndex++] = node->numNbrs;

/*
 *                      Use the current tag to check if there is a mapping to
 *                      a previous tag (due to compression of the nodeKeys
 *                      array).  If so, we need to send along the old node
 *                      index.
 */
                        key.newTag.domainID = node->myTag.domainID;
                        key.newTag.index    = node->myTag.index;

                        mapping = (TagMap_t *)bsearch(&key, home->tagMap,
                                                      home->tagMapEnts,
                                                      sizeof(TagMap_t),
                                                      TagMapCompareNew);

                        if (mapping != (TagMap_t *)NULL) {
                            outBuf[bufIndex++] = mapping->oldTag.index;
                        } else {
                            outBuf[bufIndex++] = -1;
                        }

                        for (k = 0; k < nbrCount; k++) {

/*
 *                          IMPORTANT!  DO NOT REMOVE the following NULL
 *                          pointer check!  On the BG/P system, this loop was
 *                          having issues (core dumps) accessing the nbrTag
 *                          array.  It appears to be a compiler bug causing
 *                          the problem.  Without the check for a NULL pointer
 *                          before accessing nbrTag, we can have problems.
 */
                            if (node->nbrTag == (Tag_t *)NULL) {
                                Fatal("nbrTag NULL in CommSendGhosts");
                            }

                            outBuf[bufIndex++] = node->nbrTag[k].domainID;
                            outBuf[bufIndex++] = node->nbrTag[k].index;
 
                            outBuf[bufIndex++] = node->burgX[k];
                            outBuf[bufIndex++] = node->burgY[k];
                            outBuf[bufIndex++] = node->burgZ[k];

                            outBuf[bufIndex++] = node->nx[k];
                            outBuf[bufIndex++] = node->ny[k];
                            outBuf[bufIndex++] = node->nz[k];

                            outBuf[bufIndex++] = node->armfx[k];
                            outBuf[bufIndex++] = node->armfy[k];
                            outBuf[bufIndex++] = node->armfz[k];
                        }

                        outBuf[bufIndex++] = node->x;
                        outBuf[bufIndex++] = node->y;
                        outBuf[bufIndex++] = node->z;

                        outBuf[bufIndex++] = node->fX;
                        outBuf[bufIndex++] = node->fY;
                        outBuf[bufIndex++] = node->fZ;

                        outBuf[bufIndex++] = node->vX;
                        outBuf[bufIndex++] = node->vY;
                        outBuf[bufIndex++] = node->vZ;

                        outBuf[bufIndex++] = node->oldvX;
                        outBuf[bufIndex++] = node->oldvY;
                        outBuf[bufIndex++] = node->oldvZ;

                        outBuf[bufIndex++] = (real8)node->constraint;
                        outBuf[bufIndex++] = (real8)node->flags;

#ifdef _ARLFEM
                        outBuf[bufIndex++] = node->surfaceNorm[0];
                        outBuf[bufIndex++] = node->surfaceNorm[1];
                        outBuf[bufIndex++] = node->surfaceNorm[2];
#endif  /* ifdef _ARLFEM */

                        node = node->nextInCell;
                        nodesPacked++;

                    }  /* end while (node) */

                    if (nodesPacked != cell->nodeCount) {
                        Fatal("CommPackGhosts: dom %d, remDom %d, cell %d, "
                              "cell->nodeCount %d doesn't match cell->nodeQ, "
                              "nodesPacked value %d", home->myDomain,
                              domainIndex, cellIndex, cell->nodeCount,
                              nodesPacked);
                    }

#ifdef ESHELBY
                    // Loop through all inclusions in this cell and pack
                    // the necessary info into the buffer for the remote domain

                    int incsPacked = 0;
                    int incIndex   = cell->inclusionQ;

                    while (incIndex >= 0) 
                    {
                        incsPacked++;

                        EInclusion_t *inclusion = &home->eshelbyInclusions[incIndex];

                        outBuf[bufIndex++] = (real8)inclusion->id;

                        outBuf[bufIndex++] = inclusion->radius[0];
                        outBuf[bufIndex++] = inclusion->radius[1];
                        outBuf[bufIndex++] = inclusion->radius[2];

                        outBuf[bufIndex++] = inclusion->position[0];
                        outBuf[bufIndex++] = inclusion->position[1];
                        outBuf[bufIndex++] = inclusion->position[2];

                        outBuf[bufIndex++] = inclusion->rotation[0][0];
                        outBuf[bufIndex++] = inclusion->rotation[0][1];
                        outBuf[bufIndex++] = inclusion->rotation[0][2];

                        outBuf[bufIndex++] = inclusion->rotation[1][0];
                        outBuf[bufIndex++] = inclusion->rotation[1][1];
                        outBuf[bufIndex++] = inclusion->rotation[1][2];

                        outBuf[bufIndex++] = inclusion->strainField[0];
                        outBuf[bufIndex++] = inclusion->strainField[1];
                        outBuf[bufIndex++] = inclusion->strainField[2];
                        outBuf[bufIndex++] = inclusion->strainField[3];
                        outBuf[bufIndex++] = inclusion->strainField[4];
                        outBuf[bufIndex++] = inclusion->strainField[5];

                        incIndex = inclusion->nextInCell;

                    }  /* end while (incIndex >= 0) */

                    if (incsPacked != cell->inclusionCount) {
                        Fatal("%s::%s() ln=%d : dom %d, remDom %d, cell %d, cell->inclusionCount %d doesn't match cell->inclusionQ, incsPacked value %d",
                              __FILE__, __func__, __LINE__, home->myDomain, domainIndex, cellIndex, cell->inclusionCount, incsPacked);
                    }
#endif  // ESHELBY

                }   /* end for (j = 0 ...) */

                outBuf[0] = totNodeCount;
                outBuf[1] = totInclusionCount;

                remDom->outBuf = (char *)outBuf;
                remDom->outBufLen = outBufLen;

            }  /* end for (i = threadIterStart; ... ) */

        }  /* end "omp parallel" section */

#endif
        return;
}


/*-------------------------------------------------------------------------
 *
 *  Function    : CommUnpackGhosts
 *  Description : For each remote domain, unpack the comm packet which was
 *                just received into nodes, and queue the nodes on the
 *                ghost node queue
 *
 *  Updates:   09/06/01 - add invoice stuff, to support velocity comm - t.p.
 *             09/14/01 - changed name from CommUnpackNodes, and moved into
 *                        CommSendGhosts.c file. Converted to arbitrary
 *                        number of arms for each node - t.p.
 *
 *-------------------------------------------------------------------------*/
static void CommUnpackGhosts(Home_t *home) 
{
#ifdef PARALLEL
        int    tmpIndex;
        Node_t **tmpNodePtrs = (Node_t **)NULL;;
        DomEntityInfo_t *domEntityInfo = (DomEntityInfo_t *)NULL;

        int tmpNodeIndex      = 0;
        int tmpInclusionIndex = 0;

#ifdef ESHELBY
        tmpInclusionIndex = home->totInclusionCount;
#endif

        int numRemoteDomains = home->remoteDomainCount;

        if (numRemoteDomains > 0) {
            domEntityInfo = (DomEntityInfo_t *)malloc(numRemoteDomains * sizeof(DomEntityInfo_t));
        }

/*
 *      First do a quick loop, pull out the number of entities being
 *      sent from the remote domains and set up some arrays we can
 *      use to minimize the amount of locking we'll have to do later
 */
        for (tmpIndex = 0; tmpIndex < numRemoteDomains; tmpIndex++) {
            int            tmpDomIndex;
            real8          *tmpBufPtr;
            RemoteDomain_t *tmpDomPtr;

            tmpDomIndex = home->remoteDomains[tmpIndex];
            tmpDomPtr = home->remoteDomainKeys[tmpDomIndex];

            tmpBufPtr = (real8 *)tmpDomPtr->inBuf;

            domEntityInfo[tmpIndex].nodeCount = (int)tmpBufPtr[0];
            domEntityInfo[tmpIndex].nodeIndex = tmpNodeIndex;
            domEntityInfo[tmpIndex].inclusionIndex = tmpInclusionIndex;

            tmpNodeIndex += (int)tmpBufPtr[0];
            tmpInclusionIndex += (int)tmpBufPtr[1];
        }

#ifdef ESHELBY
        // Increase the size of the inclusion array if needed.  Inclusions
        // from each remote domain will be inserted in different portions
        // of the array as defined by the remoteEntities[*].nodeIndex value
        // set above.

        if (tmpInclusionIndex != home->totInclusionCount) 
        {
            int tmpSize = tmpInclusionIndex * sizeof(EInclusion_t);
            home->eshelbyInclusions = (EInclusion_t *)realloc(home->eshelbyInclusions, tmpSize);
            home->totInclusionCount = tmpInclusionIndex;
        }
#endif // ESHELBY

/*
 *      Create a temporary array of node pointers.  Nodes from each remote
 *      domain will be placed in different portions of this array as
 *      defined by the remoteEntities[*].nodeIndex value set above.
 */
        if (tmpNodeIndex != 0) {
            tmpNodePtrs = (Node_t **)malloc(tmpNodeIndex * sizeof(Node_t *));
        }

/*
 *      We need to loop through all the remote domains unpacking a buffer
 *      from each.  Since these loop iterations will not take the same
 *      amount of processing time, don't use the OpenMP "parallel for"
 *      construct which assumes that loop iterations do take the same
 *      amout of time.
 */
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int i, threadID, threadIterStart, threadIterEnd;

            GetThreadIterationIndices(numRemoteDomains, &threadID,
                                      &threadIterStart, &threadIterEnd);

/*
 *          Have each thread loop through a portion of the remote domains
 *          unpacking the ghost entities from that domain.
 */
            for (i = threadIterStart; i < threadIterEnd; i++) {
                int            j, newSize;
                int            domNodeIndex, domNodeCount;
                int            domainIndex, bufIndex, maxTagIndex;
                int            numNodes;
                int            numExpCells;
                int            nextNodeIndex;
                int            domInclusionIndex;
                int            numInclusions, nextInclusionIndex;
                real8          *inBuf;
                RemoteDomain_t *remDom;

                domainIndex = home->remoteDomains[i];
                remDom = home->remoteDomainKeys[domainIndex];

/*
 *              Set indices into the various arrays into which we'll store
 *              entities from this remote domain.
 */
                domNodeCount = domEntityInfo[i].nodeCount;
                domNodeIndex = domEntityInfo[i].nodeIndex;

                domInclusionIndex = domEntityInfo[i].inclusionIndex;

/*
 *              Free the old remDom->nodeKeys table, if any
 */
                if (remDom->maxTagIndex) {
                    remDom->maxTagIndex = 0;
                    free(remDom->nodeKeys);
                    remDom->nodeKeys = (Node_t **)NULL;
                }

/*
 *              We know how many ghost nodes are being imported from this
 *              remote domain, so pull as many nodes as needed from the
 *              free node queue right now rather than repeatedly lock and
 *              unlock access to the free node queue over and over later
 *              on.
 */
#ifdef _OPENMP
#pragma omp critical (FREE_NODE_Q)
#endif
                {
                    int nIndex;

                    for (nIndex = 0; nIndex < domNodeCount; nIndex++) {
                        tmpNodePtrs[domNodeIndex+nIndex] = PopFreeNodeQ(home);
                    }

                }  /* end "omp critical" section */

/*
 *              First unpack the counts of entities in the buffer, the size
 *              of the maximum tag index, and the number of cells being sent.
 */
                inBuf = (real8 *)remDom->inBuf;

                bufIndex = 0;

                numNodes            = inBuf[bufIndex++];
                numInclusions       = inBuf[bufIndex++];
                maxTagIndex         = inBuf[bufIndex++];
                numExpCells         = inBuf[bufIndex++];

                remDom->maxTagIndex = maxTagIndex;

/*
 *              Allocate and initialize the new nodeKeys table for this
 *              remote domain
 */
                newSize = maxTagIndex * sizeof(Node_t *);
                remDom->nodeKeys = (Node_t **)calloc(1, newSize);

                nextNodeIndex = domNodeIndex;
                nextInclusionIndex = domInclusionIndex;

/*
 *              Loop through the cells exported from this remote domain
 */
                for (j = 0; j < numExpCells; j++) {
                    int    k;
                    int    cellIndex, cellNodeIndex, cellNodeCount;
                    int    cellInclusionCount;
                    Cell_t *cell;
                    Node_t *prevNode, *node;

                    cellIndex     = inBuf[bufIndex++];
                    cellNodeCount = inBuf[bufIndex++];
                    cellInclusionCount  = inBuf[bufIndex++];

                    cellNodeIndex = nextNodeIndex;
                    prevNode = (Node_t *)NULL;

                    cell = LookupCell(home, cellIndex);

                    if (!cell) {
                        int cX, cY, cZ;
                        DecodeCellIdx(home, cellIndex, &cX, &cY, &cZ);
                        Fatal("%s: received an unknown cell %d, (%d,%d,%d)",
                              "CommUnpackGhosts", cellIndex, cX, cY, cZ);
                    }

/*
 *                  Loop through the nodes provided by the remote domain for
 *                  this cell and unpack the data into the structure
 */
                    for (k = 0; k < cellNodeCount; k++) {
                        int    numNbrs, segID;
                        int    oldTagIndex;

                        node = tmpNodePtrs[nextNodeIndex++];

                        node->myTag.domainID   = domainIndex;
                        node->myTag.index      = inBuf[bufIndex++];
                        numNbrs                = inBuf[bufIndex++];
                        oldTagIndex            = inBuf[bufIndex++];

/*
 *                      If there is a mapping for this node to an older tag,
 *                      add the mapping to the tagMap.  All other nodes
 *                      with nbrTag pointing to the old tag for this node
 *                      will be fixed when we have pulled in all ghost nodes
 *                      and have the complete tag mapping list.
 */
                        if (oldTagIndex >= 0) {
                            Tag_t oldTag;
                            oldTag.domainID = domainIndex;
                            oldTag.index    = oldTagIndex;
#ifdef _OPENMP
#pragma omp critical (ADD_TAG_MAPPING)
#endif
                            {
                                AddTagMapping(home, &oldTag, &node->myTag);
                            }
                        }

                        AllocNodeArms(node, numNbrs);

                        for (segID = 0; segID < numNbrs; segID++) {

                            node->nbrTag[segID].domainID = inBuf[bufIndex++];
                            node->nbrTag[segID].index    = inBuf[bufIndex++];

                            node->burgX[segID] = inBuf[bufIndex++];
                            node->burgY[segID] = inBuf[bufIndex++];
                            node->burgZ[segID] = inBuf[bufIndex++];

                            node->nx[segID] = inBuf[bufIndex++];
                            node->ny[segID] = inBuf[bufIndex++];
                            node->nz[segID] = inBuf[bufIndex++];

                            node->armfx[segID] = inBuf[bufIndex++];
                            node->armfy[segID] = inBuf[bufIndex++];
                            node->armfz[segID] = inBuf[bufIndex++];
                        }

                        node->x = inBuf[bufIndex++];
                        node->y = inBuf[bufIndex++];
                        node->z = inBuf[bufIndex++];

                        node->fX = inBuf[bufIndex++];
                        node->fY = inBuf[bufIndex++];
                        node->fZ = inBuf[bufIndex++];

                        node->vX = inBuf[bufIndex++];
                        node->vY = inBuf[bufIndex++];
                        node->vZ = inBuf[bufIndex++];

                        node->oldvX = inBuf[bufIndex++];
                        node->oldvY = inBuf[bufIndex++];
                        node->oldvZ = inBuf[bufIndex++];

                        node->constraint = (int)inBuf[bufIndex++];
                        node->flags = (int)inBuf[bufIndex++];
#ifdef _ARLFEM
                        node->surfaceNorm[0] = inBuf[bufIndex++];
                        node->surfaceNorm[1] = inBuf[bufIndex++];
                        node->surfaceNorm[2] = inBuf[bufIndex++];
#endif
                        node->native = 0;
                        node->cellIdx = cellIndex;

/*
 *                      This block of nodes will be linked into the cell node
 *                      queue later, but we set up the linkage for it now.
 */
                        node->nextInCell = prevNode;
                        prevNode = node;

/*
 *                      Add the node to the remote domain's nodeKeys array.
 */
                        remDom->nodeKeys[node->myTag.index] = node;

                    }  /* end for (k = 0; ...)  */

/*
 *                  Add this group of nodes to the cell's node queue.
 *                  Linkage between the nodes in the block is already
 *                  set up, so we just need to prepend this block to the
 *                  cell queue by linking the first node in the block
 *                  to the previous cell queue, repoint the cell queue
 *                  to the last node in this block, and update the node
 *                  count.
 *
 *                  Note: we only need to lock a single cell's node queue
 *                        here so it would be nicer to use the cell lock
 *                        rather than a "critical" code section, but for
 *                        now locks have been initialized only for the
 *                        'native' cells.  If we initialize the locks for
 *                        ghost cells as well as native cells, we can
 *                        avoid the critical code section below.
 */
                    if (cellNodeCount > 0) {

#if 0
                        LOCK(&cell->cellLock);

                        tmpNodePtrs[cellNodeIndex]->nextInCell = cell->nodeQ;
                        cell->nodeQ = node;
                        cell->nodeCount += cellNodeCount;

                        UNLOCK(&cell->cellLock);
#else
#ifdef _OPENMP
#pragma omp critical (UPDATE_CELL_NODEQ)
#endif
                        {
                            tmpNodePtrs[cellNodeIndex]->nextInCell = cell->nodeQ;
                            cell->nodeQ = node;
                            cell->nodeCount += cellNodeCount;
                        }
#endif
                    }

#ifdef ESHELBY
                    // Loop through all the inclusions for this cell from the
                    // remote domain and fill in the inclusion array data

                    int prevInclusionIndex = -1;
                    int cellInclusionIndex = nextInclusionIndex;

                    EInclusion_t *incList = home->eshelbyInclusions;

                    for (int k=0; (k<cellInclusionCount); k++) 
                    {

                        // When the previous list of ghost inclusions was wiped
                        // clean, the array was only resized if the total
                        // count of inclusions dropped to zero, so it's likely
                        // this realloc will have to do nothing and hence is
                        // not much of an expense.

                        EInclusion_t *inclusion = &incList[nextInclusionIndex];

                        inclusion->id             = (int)inBuf[bufIndex++];

                        inclusion->radius[0]      = inBuf[bufIndex++];
                        inclusion->radius[1]      = inBuf[bufIndex++];
                        inclusion->radius[2]      = inBuf[bufIndex++];

                        inclusion->position[0]    = inBuf[bufIndex++];
                        inclusion->position[1]    = inBuf[bufIndex++];
                        inclusion->position[2]    = inBuf[bufIndex++];

                        inclusion->rotation[0][0] = inBuf[bufIndex++];
                        inclusion->rotation[0][1] = inBuf[bufIndex++];
                        inclusion->rotation[0][2] = inBuf[bufIndex++];

                        inclusion->rotation[1][0] = inBuf[bufIndex++];
                        inclusion->rotation[1][1] = inBuf[bufIndex++];
                        inclusion->rotation[1][2] = inBuf[bufIndex++];

                        inclusion->strainField[0] = inBuf[bufIndex++];
                        inclusion->strainField[1] = inBuf[bufIndex++];
                        inclusion->strainField[2] = inBuf[bufIndex++];
                        inclusion->strainField[3] = inBuf[bufIndex++];
                        inclusion->strainField[4] = inBuf[bufIndex++];
                        inclusion->strainField[5] = inBuf[bufIndex++];

                        inclusion->cellID     = cellIndex;
                        inclusion->nextInCell = prevInclusionIndex;

                        nextInclusionIndex++;
                    }

                    // Add this group of inclusions to the cell's inclusion
                    // queue.  Linkage between the inclusions in the block
                    // is already set up, so we just need to prepend this
                    // block to the cell queue by linking the first inclusion
                    // in the block to the previous cell queue and repoint the
                    // cell queue to the last inclusion in this block
                    // 
                    // Note: we only need to lock a single cell's node queue
                    //       here so it would be nicer to use the cell lock
                    //       rather than a "critical" code section, but for
                    //       now locks have been initialized only for the
                    //       'native' cells.  If we initialize the locks for
                    //       ghost cells as well as native cells, we can
                    //       avoid the critical code section below.

                    if (cellInclusionCount > 0) 
                    {
                        EInclusion_t *inclusion = &incList[cellInclusionIndex];

#if 0
                        LOCK(&cell->cellLock);

                        inclusion->nextInCell = cell->inclusionQ;
                        cell->inclusionQ = cellInclusionIndex;

                        UNLOCK(&cell->cellLock);
#else
#ifdef _OPENMP
#pragma omp critical (UPDATE_CELL_INCLUSIONQ)
#endif
                        {
                            inclusion->nextInCell = cell->inclusionQ;
                            cell->inclusionQ = cellInclusionIndex;
                        }
#endif
                    }
#endif  // ESHELBY

                }  /* end for (j = 0; ...) */

/*
 *              Add all the ghost nodes from this remote domain onto the
 *              local ghost node queue
 */
#ifdef _OPENMP
#pragma omp critical (GHOST_NODE_Q)
#endif
                {
                    for (j = 0; j < domNodeCount; j++) {
                        PushGhostNodeQ(home, tmpNodePtrs[domNodeIndex+j]);
                    }
                }

            }  /* end for (i = threadIterStart; ... ) */

        }  /* end "omp parallel" section */

/*
 *      Free up any temporary buffers that were allocated
 */
        if (tmpNodePtrs != (Node_t **)NULL) {
            free(tmpNodePtrs);
        }

        if (domEntityInfo != (DomEntityInfo_t *)NULL) {
            free(domEntityInfo);
        }

#endif
        return;

}


/*------------------------------------------------------------------------
 *
 *      Function:    CommSendGhosts
 *      Description: Driver function to send nodal data for local nodes
 *                   to neighboring domains that export the nodes as
 *                   as ghosts, and to receive similar data for nodes
 *                   this domain maintains as ghosts.
 *
 *-----------------------------------------------------------------------*/
void CommSendGhosts(Home_t *home) 
{
#ifdef PARALLEL
        int            i, isrc, domainIdx, idst, idom, valCount;
        int            localBuffers = 0;
        int            remDomID, totRemDomCount;
        RemoteDomain_t *remDom;


        TimerStart(home, COMM_SEND_GHOSTS);

/*
 *      Compress the nodeKeys array.  This is important when using threads
 *      because if the nodeKeys array grows parse, the array will likely
 *      be more sparse at the upper end (since low numbered node keys are
 *      re-used preferentially).  Threaded loops over the sparse list will
 *      then starve some threads while others are overloaded
 */
        CompressNodeKeyArray(home);

/*
 *      All ghost nodes (including secondary ghosts) have been
 *      recycled, and the nodeKeys arrays for the primary remote
 *      domains will be reallocated as necessary when the primary
 *      ghosts are received, but we have to free the nodeKeys arrays
 *      for the secondary remote domains here or we risk leaving
 *      pointers to structures that have already been freed.
 */
        totRemDomCount = home->remoteDomainCount +
                         home->secondaryRemoteDomainCount;

        for (i = home->remoteDomainCount; i < totRemDomCount; i++) {

            remDomID = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[remDomID];

            if (remDom == (RemoteDomain_t *)NULL) {
                Fatal("Missing rmeote domain struct!");
            }

            if (remDom->maxTagIndex) {
                free(remDom->nodeKeys);
                remDom->nodeKeys = (Node_t **)NULL;
                free(remDom);
                home->remoteDomainKeys[remDomID] = (RemoteDomain_t *)NULL;
            }

            home->secondaryRemoteDomainCount--;
        }

/*
 *      Pre-issue receives of message lengths from each neighbor
 */
        for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {

            domainIdx = home->remoteDomains[isrc];
            remDom = home->remoteDomainKeys[domainIdx];

            MPI_Irecv(&remDom->inBufLen, 1, MPI_INT, domainIdx, MSG_GHOST_LEN,
                      MPI_COMM_WORLD, &home->inRequests[isrc]);
        }

/*
 *      Package up nodal data for neighboring domains and send
 *      out the buffer lengths
 */
        CommPackGhosts(home);

        for (idst = 0; idst < home->remoteDomainCount; idst++) {

            domainIdx = home->remoteDomains[idst];
            remDom = home->remoteDomainKeys[domainIdx];

            MPI_Isend(&remDom->outBufLen, 1, MPI_INT, domainIdx, MSG_GHOST_LEN,
                      MPI_COMM_WORLD, &home->outRequests[idst]);

            localBuffers += remDom->outBufLen;
        }

/*
 *      Wait for the length sends/receives to complete
 */
        MPI_Waitall(home->remoteDomainCount, home->outRequests,home->outStatus);
        MPI_Waitall(home->remoteDomainCount, home->inRequests,home->inStatus);

/*
 *      Allocate appropriately sized buffers for the incoming messages
 *      and pre-issue receives for buffers from all neighboring domains
 */
        for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {

            domainIdx = home->remoteDomains[isrc];
            remDom = home->remoteDomainKeys[domainIdx];

            valCount = remDom->inBufLen / sizeof(real8);
            remDom->inBuf = (char *)malloc(remDom->inBufLen);
            MPI_Irecv( (real8 *) remDom->inBuf, valCount, MPI_DOUBLE, domainIdx, 
                      MSG_GHOST, MPI_COMM_WORLD, &home->inRequests[isrc]);

            localBuffers += remDom->inBufLen;
        }

/*
 *      Send local data out to all the neighboring domains
 */
        for (idst = 0; idst < home->remoteDomainCount; idst++) {

            domainIdx = home->remoteDomains[idst];
            remDom = home->remoteDomainKeys[domainIdx];

            valCount = remDom->outBufLen / sizeof(real8);
            MPI_Isend( (real8 *) remDom->outBuf, valCount, MPI_DOUBLE, domainIdx, 
                      MSG_GHOST, MPI_COMM_WORLD, &home->outRequests[idst]);
        }

/*
 *      Wait for all data buffer sends/receives to complete and unpack
 *      the data
 */
        MPI_Waitall(home->remoteDomainCount, home->outRequests,home->outStatus);
        MPI_Waitall(home->remoteDomainCount, home->inRequests,home->inStatus);

        CommUnpackGhosts(home);

/*
 *      Release the message buffers...
 */
        for (idom = 0; idom < home->remoteDomainCount; idom++) {

            domainIdx = home->remoteDomains[idom];
            remDom = home->remoteDomainKeys[domainIdx];

            free(remDom->inBuf);
            free(remDom->outBuf);

            remDom->inBuf = (char *)NULL;
            remDom->outBuf = (char *)NULL;

            remDom->inBufLen = 0;
            remDom->outBufLen = 0;
        }

/*
 *      Previously the tagMap was sorted for lookups using the 'new' node
 *      tags; now we're going to need to do lookups using the 'old' tags,
 *      so we need to resort the tag mapping list before secondary ghost
 *      communication.
 */
        if (home->tagMapEnts > 0) {
            qsort(home->tagMap, home->tagMapEnts, sizeof(TagMap_t),
                  TagMapCompareOld);
        }

/*
 *      Turns out, in order to be sure we do all the necessary
 *      direct segment-to-segment force interactions, each domain
 *      needs a layer of secondary ghosts which are nodes outside
 *      any native cell or cell immediately neighboring a native
 *      one, but connected to a primary ghost (i.e. a ghost within
 *      a native or immediately adjoining cell).
 *
 *      So, we have to do one more communication to get those
 *      secondary ghosts.
 */
        CommSendSecondaryGhosts(home);

/*
 *      In case any of the domains compressed their nodeKeys array above
 *      and retagged any nodes, we now have to reconcile those changes
 *      in the local and ghost node arrays.
 */
        ReconcileCompressedNodeKeys(home);


        TimerStop(home, COMM_SEND_GHOSTS);

#endif  /* ifdef PARALLEL */

/*
 *      Measure dead time after ghost node communications.
 */
#ifdef PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, GHOST_COMM_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, GHOST_COMM_BARRIER);
#endif
#endif
        return;
}
