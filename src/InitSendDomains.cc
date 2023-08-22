/*****************************************************************************
 *
 *      Module:      InitSendDomains.c
 *      Description: This module contains the functions necessary for
 *                   distributing the initial nodal data (read from a
 *                   text restart file) from processor zero to all
 *                   remote domains. 
 *
 *      Included functions:
 *          PackInitialNodeData()
 *          RecvInitialNodeData()
 *          SendInitialNodeData()
 *          UnpackInitialNodeData()
 *
 *****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "mpi_portability.h"

#include "Init.h"
#include "Home.h"
#include "InData.h"
#include "Comm.h"
#include "Util.h"
#include "Node.h"
#include "QueueOps.h"

/*---------------------------------------------------------------------------
 *
 *      Function:     PackInitialNodeData
 *      Description:  Pack all data from the current block of nodes that is
 *                    relevant to the specified domain into a buffer for
 *                    trasnmission to a remote domain.  Data includes generic
 *                    parameters, and node info (node tags of its neighbor
 *                    nodes, burgers vector for node, glide plane orientation,
 *                    coordinates, constraint type).
 *      Arguments:
 *          inData      structure with initial nodal data
 *          nodeList    Integer array of indices into the inData node
 *                      list for nodes to be packed into the buffer
 *          nodeCount   Number of node indices in the <nodeList> array.
 *          buf         pointer to location in which to return to the
 *                      caller a pointer to the buffer of data to
 *                      be sent to the remote domain
 *          bufEnts     pointer to the location in which to return
 *                      to the caller the number of values packed into
 *                      the buffer
 *
 *-------------------------------------------------------------------------*/
static void PackInitialNodeData(InData_t *inData, int *nodeList,
                                int nodeCount, real8 **buf, int *bufEnts)
{
        int    i, armCount, armID, numNbrs;
        int    bufIndex, maxIndex, nIndex;
        real8  *dataBuf;
   
/*
 *      To calculate the buffer size we first need a count of all
 *      the arms associated with all nodes to be sent to the remote
 *      domain.  
 */
        armCount = 0;

        for (i = 0; i < nodeCount; i++) {
            armCount += inData->node[nodeList[i]].numNbrs;
        }

        maxIndex = INIT_PARAMS + (nodeCount * INIT_VALS_PER_NODE) +
                   (armCount * INIT_VALS_PER_ARM);

        dataBuf = (real8 *) malloc(maxIndex * sizeof(real8));

/*
 *      Pack the node counts as the first items in the buffer
 *      then pack all the nodal data.
 *
 *      Note: The remote domain will retag the nodes it receives
 *      as necessary and reconcile the old/new tags with neighboring
 *      domains later.
 */
        bufIndex = 0;

        dataBuf[bufIndex++] = (real8)nodeCount;
   
        for (i = 0; i < nodeCount; i++) {

            nIndex = nodeList[i];   /* index into the inData node list */

            dataBuf[bufIndex++] = (real8)inData->node[nIndex].myTag.domainID;
            dataBuf[bufIndex++] = (real8)inData->node[nIndex].myTag.index;

            numNbrs = inData->node[nIndex].numNbrs;
            dataBuf[bufIndex++] = (real8)numNbrs;

            for (armID = 0; armID < numNbrs; armID++) {

/*
 *              Tag of the node at the terminating end of the arm.
 */
                dataBuf[bufIndex++] =
                        (real8)inData->node[nIndex].nbrTag[armID].domainID;
                dataBuf[bufIndex++] =
                        (real8)inData->node[nIndex].nbrTag[armID].index;
/*
 *              bx, by, bz
 */
                dataBuf[bufIndex++]=inData->node[nIndex].burgX[armID];
                dataBuf[bufIndex++]=inData->node[nIndex].burgY[armID];
                dataBuf[bufIndex++]=inData->node[nIndex].burgZ[armID];
/*
 *              nx, ny, nz
 */
                dataBuf[bufIndex++]=inData->node[nIndex].nx[armID];
                dataBuf[bufIndex++]=inData->node[nIndex].ny[armID];
                dataBuf[bufIndex++]=inData->node[nIndex].nz[armID];
            }

            dataBuf[bufIndex++] = (real8)inData->node[nIndex].constraint;

            dataBuf[bufIndex++] = inData->node[nIndex].x; 
            dataBuf[bufIndex++] = inData->node[nIndex].y; 
            dataBuf[bufIndex++] = inData->node[nIndex].z; 
        }
   
        if (bufIndex != maxIndex) {
            Fatal("%s: buffer size %d did not match expected size %d",
                  "PackInitialNodeData", bufIndex, maxIndex);
        }

/*
 *      Return the buffer and its length to the caller
 */
        *buf     = dataBuf;
        *bufEnts = bufIndex;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     AddTagMapping
 *      Description:  
 *
 *-------------------------------------------------------------------------*/
void AddTagMapping(Home_t *home, Tag_t *oldTag, Tag_t *newTag)
{
        int      newSize;
        TagMap_t *mapping;

        if (home->tagMapEnts >= home->tagMapSize) {
              home->tagMapSize += NEW_NODEKEY_INC;
              newSize = home->tagMapSize * sizeof(TagMap_t);
              home->tagMap = (TagMap_t *)realloc(home->tagMap, newSize);
        }

        mapping = &home->tagMap[home->tagMapEnts];

        mapping->oldTag.domainID = oldTag->domainID;
        mapping->oldTag.index    = oldTag->index;

        mapping->newTag.domainID = newTag->domainID;
        mapping->newTag.index    = newTag->index;

        home->tagMapEnts++;

        return;
}


void AddTagMappingSorted(Home_t *home, Tag_t *oldTag, Tag_t *newTag,
                         int (*cmpFunc)(const void *, const void *))
{
        int      i, newSize, cmpResult;
        TagMap_t *mapping;
        TagMap_t tmpMapping;

        if (home->tagMapEnts >= home->tagMapSize) {
              home->tagMapSize += NEW_NODEKEY_INC;
              newSize = home->tagMapSize * sizeof(TagMap_t);
              home->tagMap = (TagMap_t *)realloc(home->tagMap, newSize);
        }

        mapping = &home->tagMap[home->tagMapEnts];

        mapping->oldTag.domainID = oldTag->domainID;
        mapping->oldTag.index    = oldTag->index;

        mapping->newTag.domainID = newTag->domainID;
        mapping->newTag.index    = newTag->index;

/*
 *      Shift the new tag mapping down the sorted list until it is in the
 *      proper position.
 */
        for (i = home->tagMapEnts; i > 0; i--) {

            cmpResult = cmpFunc(&home->tagMap[i], &home->tagMap[i-1]);

            if (cmpResult < 0) {
                memcpy(&tmpMapping, &home->tagMap[i], sizeof(TagMap_t));
                memcpy(&home->tagMap[i], &home->tagMap[i-1], sizeof(TagMap_t));
                memcpy(&home->tagMap[i-1], &tmpMapping, sizeof(TagMap_t));
            } else {
                break;
            }
        }

        home->tagMapEnts++;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     ExtendNodeKeys
 *      Description:  
 *
 *-------------------------------------------------------------------------*/
static void ExtendNodeKeys(Home_t *home, int newLength)
{
        int i, oldLength;

        oldLength = home->newNodeKeyMax;
        home->newNodeKeyMax = newLength;

        home->nodeKeys = (Node_t **)realloc(home->nodeKeys,
                                            newLength * sizeof(Node_t *));

        for (i = oldLength; i < newLength; i++) {
            home->nodeKeys[i] = (Node_t *)NULL;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     GetNewTag
 *      Description:  
 *
 *-------------------------------------------------------------------------*/
static void GetNewTag(Home_t *home, Tag_t *oldTag, Tag_t *newTag,
                      int *nextAvailableTag)
{
        int     nextTag, thisDomain;

        nextTag    = *nextAvailableTag;
        thisDomain = home->myDomain;

/*
 *      If the old tag belonged to a different domain, we must
 *      give the node a new tag.
 */
        if (oldTag->domainID != thisDomain) {

            for ( ; ; nextTag++) {

/*
 *              Extend the nodekeys array if necessary.
 */
                if (nextTag >= home->newNodeKeyMax) {
                    ExtendNodeKeys(home, home->newNodeKeyMax + NEW_NODEKEY_INC);
                }

                if (home->nodeKeys[nextTag] == (Node_t *)NULL) {
                    newTag->domainID = thisDomain;
                    newTag->index = nextTag++;
                    AddTagMapping(home, oldTag, newTag);
                    break;
                }
            }

        } else {

/*
 *          The old tag belonged to this domain...  
 *
 *          If we are attempting to preserve the previous node tags
 *          then check if that tag is available for this run.
 */
            if (home->param->assignNewNodeTags == 0) {
       
/*
 *              Extend the nodekeys array if necessary.
 */
                if (oldTag->index >= home->newNodeKeyMax) {
                    ExtendNodeKeys(home, oldTag->index + NEW_NODEKEY_INC);
                }

/*
 *              If the old tag is still available, use it and return
 *              to the caller.  If the tag is not available, we just
 *              drop through and get the next available tag.
 */
                if (home->nodeKeys[oldTag->index] == (Node_t *)NULL) {
                    newTag->domainID = oldTag->domainID;
                    newTag->index    = oldTag->index;
                    if (newTag->index >= home->newNodeKeyPtr) {
                        home->newNodeKeyPtr = newTag->index + 1;
                    }
                    return;
                }
            }

/*
 *          The old tag is either no longer available, or we're not
 *          preserving the old node tags, so just use the next available
 *          tag.
 */
            for ( ; ; nextTag++) {

/*
 *              Extend node keys array if necessary
 */
                if (nextTag >= home->newNodeKeyMax) {
                    ExtendNodeKeys(home, home->newNodeKeyMax + NEW_NODEKEY_INC);
                }

                if (home->nodeKeys[nextTag] == (Node_t *)NULL) {
                    newTag->domainID = thisDomain;
                    newTag->index = nextTag++;
                    AddTagMapping(home, oldTag, newTag);
                    break;
                }
            }
        }

        if (newTag->index >= home->newNodeKeyPtr) {
            home->newNodeKeyPtr = newTag->index + 1;
        }

        *nextAvailableTag = nextTag;

        return;
}


/*---------------------------------------------------------------------------
 *
 *	Function:	UnpackInitialNodeData
 *	Description:	Pull the necessary nodal data from the given
 *			buffers sent from processor zero.  The nodal
 *			data may be transmitted in multiple buffers, so
 *			this function may be called multiple times
 *			before all this domain's nodes are received.
 *			The final buffer will contain the full domain
 *			decomposition for the problem space.
 *
 *	Arguments:
 *
 *-------------------------------------------------------------------------*/
static void UnpackInitialNodeData(Home_t *home, real8 *buf,
                                  int *nextAvailableTag)
{
        int        i, armID, numNbrs, nodesInBuf, bufIndex;
        Tag_t      oldTag, newTag;
        Node_t     *node;

/*
 *      Pull the node count out of the buffer then loop through
 *      all the nodal data provided.
 */
        bufIndex   = 0;

        nodesInBuf = (int)buf[bufIndex++];

        for (i = 0; i < nodesInBuf; i++) {

            node = PopFreeNodeQ(home);

            oldTag.domainID = (int)buf[bufIndex++];
            oldTag.index    = (int)buf[bufIndex++];

/*
 *          Determine what tag will be used to identify this node.
 *          If the user requested that node tags be preserved, the
 *          new tag returned will be the same as the old tag
 *          if at all possible.  If the new tag is in fact different
 *          from the old tag, the GetNewTag() function will add
 *          a mapping from the old to the new tag for later use
 *          when reconciling all nodes that have been retagged.
 */
            GetNewTag(home, &oldTag, &newTag, nextAvailableTag);

            node->myTag.domainID = newTag.domainID;
            node->myTag.index    = newTag.index;

            numNbrs = (int)buf[bufIndex++];

            AllocNodeArms(node, numNbrs);

            for (armID = 0; armID < numNbrs; armID++) {

                node->nbrTag[armID].domainID = (int)buf[bufIndex++];
                node->nbrTag[armID].index    = (int)buf[bufIndex++];

                node->burgX[armID] = buf[bufIndex++];
                node->burgY[armID] = buf[bufIndex++];
                node->burgZ[armID] = buf[bufIndex++];

                node->nx[armID] = buf[bufIndex++];
                node->ny[armID] = buf[bufIndex++];
                node->nz[armID] = buf[bufIndex++];
            }

            node->constraint = (int)buf[bufIndex++];

            node->x = buf[bufIndex++]; 
            node->y = buf[bufIndex++]; 
            node->z = buf[bufIndex++]; 

/*
 *          Register the node in the nodeKeys array and push the node
 *          onto the queue of nodes native to this domain.
 */
            home->nodeKeys[node->myTag.index] = node;
            PushNativeNodeQ(home, node);

        }  /* for (i = 0; i < nodesInBuf; ... ) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SendInitialNodeData
 *      Description:  Assigns all nodes currently on the inData node list
 *                    to remote domains based on the node's coordinates and
 *                    the current decomposition data, packages the data
 *                    up and ships it off to the appropriate domains.
 *                    Since the initial node data may be read in blocks,
 *                    this function may be called multiple times before
 *                    all nodal data has been distributed to the remote
 *                    domains.
 *
 *      Arguments:
 *          inData
 *          msgCount
 *          nodeLists
 *          listCounts
 *
 *-------------------------------------------------------------------------*/
void SendInitialNodeData(Home_t *home, InData_t *inData, int *msgCount,
                         int **nodeLists, int *listCounts,
                         int *nextAvailableTag)
{
        int         i, domID, numDomains, thisDomain;
        int         sendIndex, selfSendIndex;
        int         numSendBufs, numRecvBufs;
        int         *sendBufEnts=0, *sendBufDest=0, *recvBufEnts=0;
        real8       **sendBuf=0, **recvBuf=0;
#ifdef PARALLEL
        int         recvIndex;
#endif

        numDomains = home->numDomains;
        thisDomain = home->myDomain;

/*
 *      Normally, the following arrays are allocated in 
 *      InitRemoteDomains(), but that has not yet been
 *      called when we enter this, so for now, allocate
 *      the arrays temporarily
 */
#ifdef PARALLEL
        home->inRequests  = (MPI_Request *) malloc(home->numDomains *
                                                  sizeof(MPI_Request));
        home->outRequests = (MPI_Request *) malloc(home->numDomains *
                                                  sizeof(MPI_Request));
        home->inStatus    = (MPI_Status *) malloc(home->numDomains *
                                                  sizeof(MPI_Status));
        home->outStatus   = (MPI_Status *) malloc(home->numDomains *
                                                  sizeof(MPI_Status));
#endif

/*
 *      Determine the number of remote domains to  which the
 *      current domain must send data and allocate arrays
 *      for the send buffer pointers and the send buffer
 *      sizes.
 */
        numSendBufs   = 0;
        selfSendIndex = -1;

        if (listCounts != (int *)NULL) {
            for (domID = 0; domID < numDomains; domID++) {
                if (listCounts[domID] != 0) {
                    if (domID == thisDomain) {
                        selfSendIndex = numSendBufs;
                    }
                    numSendBufs++;
                }
            }

            if (numSendBufs > 0) {
                sendBuf     = (real8 **)calloc(1, numSendBufs * sizeof(real8 *));
                sendBufEnts = (int *)calloc(1, numSendBufs * sizeof(int));
                sendBufDest = (int *)calloc(1, numSendBufs * sizeof(int));
            }

/*
 *          Pack up nodal data for each remote domain to which this
 *          domain will be sending data
 */
            sendIndex = 0;

            for (domID = 0; domID < numDomains; domID++) {
                if (listCounts[domID] != 0) {
    
                    PackInitialNodeData(inData, nodeLists[domID],
                                        listCounts[domID], &sendBuf[sendIndex],
                                        &sendBufEnts[sendIndex]);

                    sendBufDest[sendIndex] = domID;
                    sendIndex++;

                    if (sendIndex >= numSendBufs) {
                        break;
                    }
                }
            }
        }  /* if (listCount != NULL) */

/*
 *      Allocate array in which to receive incoming buffer lengths.
 *      Do not include any buffer this domain has for itself; any
 *      such buffer will be handled separately.
 */
        numRecvBufs = msgCount[thisDomain] - (selfSendIndex >= 0);

        if (numRecvBufs > 0) {
            recvBuf = (real8 **)calloc(1, numRecvBufs * sizeof(real8 *));
            recvBufEnts = (int *)calloc(1, numRecvBufs * sizeof(int));
        }
        
/*
 *      Pre-issue receives of buffer lengths from any domain that will
 *      be sending data to the current domain.
 */
#ifdef PARALLEL
        for (i = 0; i < numRecvBufs; i++) {
            MPI_Irecv(&recvBufEnts[i], 1, MPI_INT, MPI_ANY_SOURCE,
                      MSG_INIT_LENS, MPI_COMM_WORLD, &home->inRequests[i]);
        }

/*
 *      Have the current domain send out the sizes of the buffers it
 *      will be transmitting.  Note: If the buffer is for itself, the
 *      data is not actually transmitted
 */
        for (i = 0; i < numSendBufs; i++) {
            if (i == selfSendIndex) {
                home->outRequests[i] = MPI_REQUEST_NULL;
            } else {
                MPI_Isend(&sendBufEnts[i], 1, MPI_INT, sendBufDest[i],
                          MSG_INIT_LENS, MPI_COMM_WORLD,
                          &home->outRequests[i]);
            }
        }
#endif

/*
 *      If this domain packed up a buffer for itself, unpack that buffer
 *      now.
 */
        if (selfSendIndex >= 0) {
            UnpackInitialNodeData(home, sendBuf[selfSendIndex],
                                  nextAvailableTag);
        }

/*
 *      Wait for all length send/receives to complete
 */
#ifdef PARALLEL
        MPI_Waitall(numSendBufs, home->outRequests, home->outStatus);
        MPI_Waitall(numRecvBufs, home->inRequests, home->inStatus);

/*
 *      Allocate the receive buffers and post the receives
 *      associated with those buffers.
 */
        for (i = 0; i < numRecvBufs; i++) {
            recvBuf[i] = (real8 *)malloc(recvBufEnts[i] * sizeof(real8));
            MPI_Irecv(recvBuf[i], recvBufEnts[i], MPI_DOUBLE,
                      home->inStatus[i].MPI_SOURCE, MSG_INIT_NODES,
                      MPI_COMM_WORLD, &home->inRequests[i]);
        }

/*
 *      Send out all the packed buffers to the appropriate
 *      remote domains.  -- Unless the destination is the current domain
 *      in which case the buffer has already been processed.
 */
        for (i = 0; i < numSendBufs; i++) {
            
            if (i == selfSendIndex) {
                home->outRequests[i] = MPI_REQUEST_NULL;
            } else {
                MPI_Isend(sendBuf[i], sendBufEnts[i], MPI_DOUBLE,
                          sendBufDest[i], MSG_INIT_NODES, MPI_COMM_WORLD,
                          &home->outRequests[i]);
            }
        }

/*
 *      Process the incoming buffers as soon as they arrive.
 */
        for (i = 0; i < numRecvBufs; i++) {
            MPI_Waitany(numRecvBufs, home->inRequests, &recvIndex,
                        home->inStatus);
            UnpackInitialNodeData(home, recvBuf[recvIndex], nextAvailableTag);
            free(recvBuf[recvIndex]);
        }

/*
 *      Wait for all buffer sends to complete.
 */
        MPI_Waitall(numSendBufs, home->outRequests, home->outStatus);

        free(home->inRequests);
        free(home->outRequests);
        free(home->inStatus);
        free(home->outStatus);

        home->inRequests  = (MPI_Request *)NULL;
        home->outRequests = (MPI_Request *)NULL;
        home->inStatus    = (MPI_Status *)NULL;
        home->outStatus   = (MPI_Status *)NULL;
#endif

        for (i = 0; i < numSendBufs; i++) {
            free(sendBuf[i]);
        }

        if (numSendBufs > 0) {
            free(sendBuf);
            free(sendBufEnts);
            free(sendBufDest);
        }

        if (numRecvBufs > 0) {
            free(recvBuf);
            free(recvBufEnts);
        }
        
	return;
}
