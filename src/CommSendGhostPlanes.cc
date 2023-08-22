/***************************************************************************
 *
 *      Module:      CommSendGhostPlanes
 *      Description: Contains functions necessary for communicating
 *                   updated segment glide plane information for locally
 *                   owned segments to remote domains that may see those
 *                   segments as ghosts.
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
 *          CommPackGhostPlanes()
 *          CommSendGhostPlanes()
 *          CommUnpackGhostPlanes()
 *
 **************************************************************************/
#include "mpi_portability.h"

#include "Home.h"
#include "RemoteDomain.h"
#include "Cell.h"
#include "Node.h"
#include "Comm.h"
#include "QueueOps.h"


#define MSG_GHOST_PLANE_LEN 750
#define MSG_GHOST_PLANE     751

#define FLTS_PER_GHOST_PLANE       7
#define FLTS_PER_GHOST_PLANE_EXTRA 1

/*-------------------------------------------------------------------------
 *
 *      Function:    CommPackGhostPlanes
 *
 *      Description: For each neighbor domain, pack into a send buffer
 *                   the glide planes for segments owned by nodes in any
 *                   cell that borders the neighbor's cells.
 *
 *------------------------------------------------------------------------*/
static void CommPackGhostPlanes(Home_t *home) 
{
#ifdef PARALLEL
        int thisDomain;

        thisDomain = home->myDomain;

/*
 *      We need to loop through all the remote domains packing a buffer
 *      for each with the glide plane data for all local nodes in cells
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
                int            j;
                int            domainIdx, cellIndex, bufIndex;
                int            segCount, segmentsPacked;
                int            outBufLen, maxBufElements;
                real8          *outBuf;
                Cell_t         *cell;
                Node_t         *node;
                RemoteDomain_t *remDom;

                domainIdx = home->remoteDomains[i];
                remDom = home->remoteDomainKeys[domainIdx];

/*
 *              Get an estimate of the total number of segments to send
 *              to this remote domain.  This estimate may be high, but
 *              never too low.
 */
                segCount = 0;

                for (j = 0; j < remDom->numExpCells; j++) {

                    cellIndex = remDom->expCells[j];
                    cell = LookupCell(home, cellIndex);

                    node = cell->nodeQ;

                    while (node) {

                        if (node->myTag.domainID == thisDomain) {
                            segCount += node->numNbrs;
                        }

                        node = node->nextInCell;
                    }
                }

/*
 *              Calculate the maximum buffer size needed to hold the data for
 *              this remote domain and allocate it.
 */
                maxBufElements = FLTS_PER_GHOST_PLANE * segCount +
                                 FLTS_PER_GHOST_PLANE_EXTRA;

                outBufLen = (maxBufElements * sizeof(real8));
                outBuf = (real8 *)malloc(outBufLen);

                bufIndex = 0;
                segmentsPacked = 0;
        
/*
 *              Send the remote domain the number of segment glide planes
 *              it will receive
 */
                bufIndex++;    /* reserve space for segment count */

                for (j = 0; j < remDom->numExpCells; j++) {
        
                    cellIndex = remDom->expCells[j];
                    cell = LookupCell(home, cellIndex);
        
/*
 *                  loop through all the nodes on the cell's node queue
 */
                    node = cell->nodeQ;

                    while (node) {
                        int k, nbrCount;

/*
 *                      Don't bother with any nodes not owned by this domain
 */
                        if (node->myTag.domainID != thisDomain) {
                            node = node->nextInCell;
                            continue;
                        }

                        nbrCount = node->numNbrs;

/*
 *                      Loop though the node's segments; send glide plane info
 *                      for any segment owned by the node.
 */
                        for (k = 0; k < nbrCount; k++) {

#ifdef _BGP
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
#endif

                            if (OrderTags(&node->myTag,&node->nbrTag[k])<0){
                                continue;
                            }

                            outBuf[bufIndex++]=(real8)node->myTag.domainID;
                            outBuf[bufIndex++]=(real8)node->myTag.index;

                            outBuf[bufIndex++]=(real8)node->nbrTag[k].domainID;
                            outBuf[bufIndex++]=(real8)node->nbrTag[k].index;
 
                            outBuf[bufIndex++]=node->nx[k];
                            outBuf[bufIndex++]=node->ny[k];
                            outBuf[bufIndex++]=node->nz[k];

                            segmentsPacked++;
                        }

                        node = node->nextInCell;
        
                    }  /* end while (node) */

                }   /* end for (j = 0 ...) */
        
                outBuf[0] = segmentsPacked;
        
                remDom->outBuf = (char *)outBuf;
                remDom->outBufLen = bufIndex * sizeof(real8);
        
            }  /* end for (i = threadIterStart; ...) */

        }  /* end "omp parallel" section */
        
#endif
        return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:    CommUnpackGhostPlanes
 *  Description: For each remote domain, unpack the comm packet which was
 *               just received, and update the glide plane info for each
 *               segment
 *
 *-------------------------------------------------------------------------*/
static void CommUnpackGhostPlanes(Home_t *home) 
{
#ifdef PARALLEL

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
            int    i, threadID, threadIterStart, threadIterEnd;

            GetThreadIterationIndices(home->remoteDomainCount, &threadID,
                                      &threadIterStart, &threadIterEnd);

/*
 *          Have each thread loop through a portion of the buffers from the
 *          remote domains and unpack the ghost segment glide plane data.
 */
            for (i = threadIterStart; i < threadIterEnd; i++) {
                int            numSegs, segIndex;
                int            bufIndex, domainIdx;
                real8          *inBuf;
                RemoteDomain_t *remDom;

                domainIdx = home->remoteDomains[i];
                remDom = home->remoteDomainKeys[domainIdx];
        
/*
 *              First item in the buffer is the number of segments for which
 *              glide plane info is being sent.
 */
                inBuf = (real8 *)remDom->inBuf;

                bufIndex = 0;

                numSegs  = inBuf[bufIndex++];
        
/*
 *              Loop though the segments in the buffer pulling off the glide
 *              plane info.
 */
                for (segIndex = 0; segIndex < numSegs; segIndex++) {
                    int    j;
                    real8  newPlane[3];
                    Tag_t  tag1, tag2;
                    Node_t *node;

                    tag1.domainID = (int)inBuf[bufIndex++];
                    tag1.index    = (int)inBuf[bufIndex++];

                    tag2.domainID = (int)inBuf[bufIndex++];
                    tag2.index    = (int)inBuf[bufIndex++];

                    newPlane[X]   = inBuf[bufIndex++];
                    newPlane[Y]   = inBuf[bufIndex++];
                    newPlane[Z]   = inBuf[bufIndex++];

/*
 *                  This domain must update the segment glide plane info
 *                  for any of the segment's endpoints for which it
 *                  has data.  Local domain may not have data on both
 *                  endpoints, so we do a litle extra work.
 */
                    node = GetNodeFromTag(home, tag1);

                    if (node != (Node_t *)NULL) {
                        for (j = 0; j < node->numNbrs; j++) {
                            if ((node->nbrTag[j].domainID == tag2.domainID) &&
                                (node->nbrTag[j].index    == tag2.index   )) {
                                break;
                            }
                        }

                        if (j < node->numNbrs) {
                            node->nx[j] = newPlane[X];
                            node->ny[j] = newPlane[Y];
                            node->nz[j] = newPlane[Z];
                        }
                    }

                    node = GetNodeFromTag(home, tag2);

                    if (node != (Node_t *)NULL) {
                        for (j = 0; j < node->numNbrs; j++) {
                            if ((node->nbrTag[j].domainID == tag1.domainID) &&
                                (node->nbrTag[j].index    == tag1.index   )) {
                                break;
                            }
                        }

                        if (j < node->numNbrs) {
                            node->nx[j] = newPlane[X];
                            node->ny[j] = newPlane[Y];
                            node->nz[j] = newPlane[Z];
                        }
                    }
                }

            }  /* end for (i = threadIterStart; ...) */

        }  /* end "omp parallel" section */
        
#endif  /* end ifdef PARALLEL */

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    CommSendGhostPlanes
 *      Description: Driver function to send glide plane data for local nodes
 *                   to neighboring domains that export the nodes as
 *                   as ghosts, and to receive similar data for nodes
 *                   this domain maintains as ghosts.
 *
 *-----------------------------------------------------------------------*/
void CommSendGhostPlanes(Home_t *home) 
{
#ifdef PARALLEL
        int            i, isrc, domainIdx, idst, idom, valCount;
        int            localBuffers = 0;
        int            remDomID, totRemDomCount;
        RemoteDomain_t *remDom;
        
        TimerStart(home, COMM_SEND_GHOSTS);

/*
 *      All ghost nodes (including secondary ghosts) have been
 *      recycled, and the nodeKeys arrays for the primary remote
 *      domains will be reallocated as necessary when the primary
 *      ghosts are received, but we have to free the nodeKeys arrays
 *      for the secondary remote domains here or we risk leaving
 *      pointers to structures that have already been freed.
 */
        totRemDomCount = home->remoteDomainCount;

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
        }

/*
 *      Pre-issue receives of message lengths from each neighbor
 */
        for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {
        
            domainIdx = home->remoteDomains[isrc];
            remDom = home->remoteDomainKeys[domainIdx];
        
            MPI_Irecv(&remDom->inBufLen, 1, MPI_INT, domainIdx, MSG_GHOST_PLANE_LEN,
                      MPI_COMM_WORLD, &home->inRequests[isrc]);
        }

/*
 *      Package up nodal data for neighboring domains and send
 *      out the buffer lengths
 */
        CommPackGhostPlanes(home);
        
        for (idst = 0; idst < home->remoteDomainCount; idst++) {
        
            domainIdx = home->remoteDomains[idst];
            remDom = home->remoteDomainKeys[domainIdx];
        
            MPI_Isend(&remDom->outBufLen, 1, MPI_INT, domainIdx, MSG_GHOST_PLANE_LEN,
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
                      MSG_GHOST_PLANE, MPI_COMM_WORLD, &home->inRequests[isrc]);

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
                      MSG_GHOST_PLANE, MPI_COMM_WORLD, &home->outRequests[idst]);
        }
        
/*
 *      Wait for all data buffer sends/receives to complete and unpack
 *      the data
 */
        MPI_Waitall(home->remoteDomainCount, home->outRequests,home->outStatus);
        MPI_Waitall(home->remoteDomainCount, home->inRequests,home->inStatus);
        
        CommUnpackGhostPlanes(home);
        
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
