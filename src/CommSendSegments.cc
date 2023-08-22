/*****************************************************************************
 *
 *      Module:      CommSendSegments.c
 *      Description: This module contains the functions necessary for
 *                   distributing partial segment forces between domains
 *                   after the full force calculations.  This is necessary
 *                   because only a single domain will calculate the
 *                   forces between any given pair of segments, and the
 *                   the domains owning nodes attached to segments whose
 *                   forces were calculated elsewhere must receive the
 *                   component of the forces computed elsewhere.
 *
 *      Included public functions:
 *          CommSendSegForces()
 *          CommRecvSegForces()
 *
 *      Includes private functions:
 *          PackSegmentData()
 *          UnpackSegmentData()
 *
 *****************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"

#ifdef PARALLEL

/*
 *      Number of values to be communicated per segment sent
 */
#define NUM_VALS_PER_SEGMENT 11

/*---------------------------------------------------------------------------
 *
 *      Function:     PackSegmentData
 *      Description:  Loop through all segments for which some portion of
 *                    the forces were computed by this domain, and pack the
 *                    segment and its forces into buffers for any remote
 *                    domain owning either node in the segment.
 *
 *      Arguments:
 *          numSendBufs    Number of domains to which this domain will
 *                         be sending segment data.
 *          segCommInfo    Array of <numSendBufs> structs containing info
 *                         (domain ID, segment count, count of segment/particle
 *                         intersect structures to send) for each remote
 *                         domain to which this domain will be sending seg
 *                         force data.
 *          cellSegLists   Array of segment arrays -- 1 segment array
 *                         per cell known to this domain.
 *          cellSegCounts  Number of segments in each of the segment arrays
 *                         in <cellSegLists>
 *          sendBufs       Array in which to return to the caller the 
 *                         pointers to the data buffers to be sent out.
 *                         This function handles allocation and packing
 *                         of individual buffers
 *          sendBufEnts    Array in which to return the count of values
 *                         packed into the corresponding send buffers.
 *
 *-------------------------------------------------------------------------*/
static void PackSegmentData(
   Home_t *home,
   int             numSendBufs ,
   SegCommInfo_t  *segCommInfo ,
   Segment_t     **cellSegLists,
   int            *cellSegCnts ,
   real8         **sendBufs    ,
   int            *sendBufEnts  )
{
    int  thisDomain = home->myDomain;
    int  numCells   = home->cellCount;
    int *currIndex  = (int *)calloc(1, numSendBufs * sizeof(int));
 
    // We know how many segments will be sent to each domain, so allocate
    // appropriately sized buffers to be sent to those domains.

    for (int i=0; (i<numSendBufs); i++) 
    {
        int segCnt           = segCommInfo[i].segCnt;
        int numIntersectData = segCommInfo[i].intersectDataCnt;

#ifdef ESHELBY
        sendBufEnts[i]    = 1 + (segCnt * NUM_VALS_PER_SEGMENT) + (numIntersectData * MAX_SEGPART_INTERSECTIONS);
#else
        sendBufEnts[i]    = 1 + (segCnt * NUM_VALS_PER_SEGMENT);
#endif

        sendBufs   [i]    = (real8 *) malloc(sendBufEnts[i] * sizeof(real8));
        sendBufs   [i][0] = (real8  ) segCnt;
        currIndex  [i]    = 1;
    }

    // We need to loop through all the segments in all the cells sending
    // the segment data out to the remote domains owning the nodal endpoints
    // each segment.  Since these loop iterations will not take the same
    // amount of processing time, don't use the OpenMP "parallel for" construct
    // which assumes that loop iterations do take the same amout of time.

#ifdef _OPENMP
#pragma omp parallel
#endif
    {

        // Each thread will be responsible for packing data buffers for
        // some subset of the remote domains.  Calculate the range of
        // remote domains for each thread now.
        // 
        // NOTE: Since each thread is responsible for packing buffers for
        // a unique subset remote domains, there will be no
        // contention in accessing the buffers, so no locking is
        // necessary.

        int threadID=0, threadIterStart=0, threadIterEnd=0;

        GetThreadIterationIndices(numSendBufs, &threadID, &threadIterStart, &threadIterEnd);

        // Loop through all the cells unless this thread is not responsible
        // for packing any buffers.  If that is the case, just set the
        // initial index for the loop over cells so the thread does nothing.

        int i=0;
        if (threadIterStart >= threadIterEnd) { i=numCells; } 

        for (/* i initialized above */; i < numCells; i++) 
        {
            Segment_t *segmentList = cellSegLists[i];

            if (segmentList == (Segment_t *)NULL) { continue; }

            int numSegs = cellSegCnts[i];

            // Loop over all segments in this cell

            for (int j=0; (j<numSegs); j++) 
            {
                Segment_t segment = segmentList[j];

                // If this domain did not compute any portion of the forces
                // for this segment, skip it.

                if (segment.forcesSet == 0) { continue; }

                int domID1 = segment.node1->myTag.domainID;
                int domID2 = segment.node2->myTag.domainID;

                if (domID1 != thisDomain) 
                {
                    // Find the buffer associated with the domain owning
                    // the first node in the segment.  If this thread is
                    // responsible for the buffer for that domain, pack
                    // the segment data into its buffer.

                    int k=0;
                    for (k=threadIterStart; (k<threadIterEnd); k++)
                        if (segCommInfo[k].domID == domID1) break;

                    if (k<threadIterEnd) 
                    {
                        real8 *buf   = sendBufs [k];
                        int    index = currIndex[k];

                        buf[index++] = (real8)segment.node1->myTag.domainID;
                        buf[index++] = (real8)segment.node1->myTag.index;
                        buf[index++] = (real8)segment.node2->myTag.domainID;
                        buf[index++] = (real8)segment.node2->myTag.index;

                        // (not sending particle intersect data)

                        buf[index++] = 0.0;
                        buf[index++] = segment.f1[0];
                        buf[index++] = segment.f1[1];
                        buf[index++] = segment.f1[2];

                        buf[index++] = segment.f2[0];
                        buf[index++] = segment.f2[1];
                        buf[index++] = segment.f2[2];

                        currIndex[k] = index;
                    }
                }

                if ((domID2 != thisDomain) && (domID2 != domID1)) 
                {
                    // Find the buffer associated with the domain owning
                    // the second node in the segment.  If this thread is
                    // responsible for the buffer for that domain, pack
                    // the segment data into its buffer.

                    int k=0;
                    for (k=threadIterStart; (k<threadIterEnd); k++)
                        if (segCommInfo[k].domID == domID2) break;

                    if (k<threadIterEnd) 
                    {
                        real8 *buf   = sendBufs [k];
                        int    index = currIndex[k];

                        buf[index++] = (real8)segment.node1->myTag.domainID;
                        buf[index++] = (real8)segment.node1->myTag.index;
                        buf[index++] = (real8)segment.node2->myTag.domainID;
                        buf[index++] = (real8)segment.node2->myTag.index;
                        buf[index++] = (real8)segment.sendParticleIntersects;
                        buf[index++] = segment.f1[0];
                        buf[index++] = segment.f1[1];
                        buf[index++] = segment.f1[2];
                        buf[index++] = segment.f2[0];
                        buf[index++] = segment.f2[1];
                        buf[index++] = segment.f2[2];

#ifdef ESHELBY

                        // If we need to send the list of particles
                        // intersected by this segment, send them along.
                        // Note: since the intersection data only needs to be
                        // sent to domains not owning a segment, and the
                        // segment is always owned by the first node in the
                        // segment data, we never need to send this info to
                        // the domain owning the first node in the above loop.

                        if (segment.sendParticleIntersects) 
                        {
                            SegPartIntersect_t *intersect = SegPartListLookup(home, &segment.node1->myTag, &segment.node2->myTag);

                            if (!intersect) 
                            {  Fatal("%s::%s() ln=%d : Error looking up particle intersect data", __FILE__, __func__, __LINE__ ); }
                           

                            // Send the whole array of intersected particles,
                            // but set unused values to -1 so the receiver
                            // knows how many are valid.

                            for (int m=0; (m<MAX_SEGPART_INTERSECTIONS); m++) 
                            {
                                if (m < intersect->numIntersections) { buf[index++] = (real8)intersect->inclusionID[m]; } 
                                else                                 { buf[index++] = -1.0; }
                            }

                        }  // if (segment.sendParticleIntersects)
#endif

                        currIndex[k] = index;

                    }  // end if (k < threadIterEnd)
                }

            }  // end for (j = 0; j < numSegs; ...)

        }  // end loop over cells

    }  // end "omp parallel" section

    // Do a quick sanity check to verify we did not overrun any
    // of the buffers.  Once the code has been sufficiently tested,
    // this check should be removed.

    for (int i=0; (i<numSendBufs); i++) 
    {
        if (currIndex[i] > sendBufEnts[i])
            Fatal("%s::%s() ln=%d : Overpacked send buffer %d.\n" "Packed %d values, expected only %d", 
                  __FILE__, __func__, __LINE__, i, currIndex[i], sendBufEnts[i]);
    }

    free(currIndex);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     UnpackSegmentData
 *      Description:  Read segment data from the provided buffer and
 *                    update the segment forces for any nodes in the
 *                    segment owned by this domain.
 *
 *                    IMPORTANT!  This function may be called from
 *                                multiple threads simultaneously, so it
 *                                MUST be thread-safe!
 *
 *      Arguments:
 *          buf    pointer to array packed with segment data
 *
 *-------------------------------------------------------------------------*/
static void UnpackSegmentData(Home_t *home, real8 *buf)
{
    int thisDomain = home->myDomain;

    // Pull the segment count out of the buffer then loop through
    // all the segment data provided.

    int bufIndex = 0;
    int numSegs  = (int)buf[bufIndex++];

    for (int i=0; (i<numSegs); i++) 
    {
        Tag_t  tag1, tag2;
        real8  f1[3], f2[3];

        tag1.domainID = (int)buf[bufIndex++];
        tag1.index    = (int)buf[bufIndex++];

        tag2.domainID = (int)buf[bufIndex++];
        tag2.index    = (int)buf[bufIndex++];

        int getParticleIntersects = (int) buf[bufIndex++];

        f1[0] = buf[bufIndex++];
        f1[1] = buf[bufIndex++];
        f1[2] = buf[bufIndex++];

        f2[0] = buf[bufIndex++];
        f2[1] = buf[bufIndex++];
        f2[2] = buf[bufIndex++];

#ifdef ESHELBY

        // If the data includes the segment/particle intersection data
        // read that in and process it now.  NOTE: The remote domain
        // sends the inclusion ID's but the local domain needs to 
        // reference the inclusions by index into the local inclusion list
        // so we have to look that up...

        if (getParticleIntersects) 
        {
            SegPartIntersect_t intersect;

            intersect.tag1.domainID    = tag1.domainID;
            intersect.tag1.index       = tag1.index;
            intersect.tag2.domainID    = tag2.domainID;
            intersect.tag2.index       = tag2.index;
            intersect.numIntersections = 0;

            for (int j=0; (j<MAX_SEGPART_INTERSECTIONS); j++) 
            {
                intersect.inclusionID[j] = (int)buf[bufIndex++];

                if (intersect.inclusionID[j] >= 0) {
                    intersect.numIntersections += 1;
                    // Below should work now.
                    intersect.inclusionIndex[j] = GetIndexFromInclusionID(home, intersect.inclusionID[j]);
                } else {
                    intersect.inclusionIndex[j] = -1;
                }
            }

            // Append the segment's intersect info to the list. NOTE: this
            // is a quick append which does not guarantee ordering in the
            // list.  The list should be resorted after the communication
            // is completed.
            // 
            // Note: Access to the global intersection list is synchronized
            //       within SegPartListAppend(), so no additional locking is
            //       needed here.

            SegPartListAppend(home, &intersect);
        }
#endif

        // If this domain owns the first node, update the
        // segment force at the that node.  Multiple remote domains
        // may have calculated partial forces for this segment, so
        // the segment forces must be locked while being updated.

        if (tag1.domainID == thisDomain) 
        {
            Node_t *node = home->nodeKeys[tag1.index];

            for (int arm=0; (arm<node->numNbrs); arm++) 
            {
                if (    (node->nbrTag[arm].domainID == tag2.domainID)
                     && (node->nbrTag[arm].index    == tag2.index   ) ) 
                {
                    LOCK(&node->nodeLock);

                    node->armfx[arm] += f1[0];
                    node->armfy[arm] += f1[1];
                    node->armfz[arm] += f1[2];

                    UNLOCK(&node->nodeLock);

                    break;
                }
            }
        }

        // If this domain owns the second node, update the
        // segment force at the that node.  Multiple remote domains
        // may have calculated partial forces for this segment, so
        // the segment forces must be locked while being updated.

        if (tag2.domainID == thisDomain) 
        {
            Node_t *node = home->nodeKeys[tag2.index];

            for (int arm=0; (arm<node->numNbrs); arm++) 
            {
                if (    (node->nbrTag[arm].domainID == tag1.domainID) 
                     && (node->nbrTag[arm].index    == tag1.index   ) ) 
                {
                    LOCK(&node->nodeLock);

                    node->armfx[arm] += f2[0];
                    node->armfy[arm] += f2[1];
                    node->armfz[arm] += f2[2];

                    UNLOCK(&node->nodeLock);

                    break;
                }
            }
        }

    }  // for (i = 0; i < numSegs; ...)

}

/*---------------------------------------------------------------------------
 *
 *      Function:     CommSendSegForces
 *
 *      Description:  This function will package up the segment forces
 *                    that were calculated locally and send them off
 *                    to the approporiate remote domains.  These
 *                    buffers are sent asynchronously and the function
 *                    does not wait until the buffers arrive at the
 *                    destination.  The associated buffers and request
 *                    structs will be freed by CommRecvSegForces()
 *                    after all incoming messages have been processed
 *                    and the routine has verified all sent buffers
 *                    have been received at their destinations.
 *
 *      Arguments:
 *          IN:  numSendBufs    Number of remote domains to which this
 *                              domain will be sending segment data
 *          IN:  segCommInfo    Array of <numSendBufs> structs containing info
 *                              (domain ID, segment count, count of
 *                              segment/particle intersect structures to send)
 *                              for each remote domain to which this domain
 *                              will be sending seg force data.
 *          IN:  cellSegLists   Array of segment arrays -- 1 segment array
 *                              per cell known to this domain.
 *          IN:  cellSegCounts  Number of segments in each of the segment
 *                              arrays in <cellSegLists>
 *          OUT: sendBufs       Pointer to location in which to return
 *                              an array of <numSendBufs> buffers
 *                              associate with buffers sent to the remote
 *                              domains.
 *          OUT: sendBufLengths Pointer to array of <numSendBufs> integers
 *                              indicating the number of 'double' values
 *                              in the corresponding send buffers.
 *          OUT: sendBufLenRequests Pointer to location in which to return
 *                              an array of <numSendBufs> MPI request buffers
 *                              associated with buffer length messages sent
 *                              to the remote domains.
 *                              
 *
 *-------------------------------------------------------------------------*/
void CommSendSegForces(Home_t *home, int numSendBufs,
                       SegCommInfo_t *segCommInfo, Segment_t **cellSegLists,
                      int *cellSegCnts, real8 ***sendBufs,
                      int **sendBufLengths, MPI_Request **sendBufLenRequests)
{
        int          *sendBufEnts =0;
        real8       **sendBuf     =0;
        MPI_Request  *sendRequests=0;

        // Allocate arrays for the send buffer pointers and the send buffer sizes.

        if (numSendBufs>0) 
        {
            sendBuf      = (real8      **)calloc(1, numSendBufs * sizeof(real8     *));
            sendBufEnts  = (int         *)calloc(1, numSendBufs * sizeof(int        ));
            sendRequests = (MPI_Request *)calloc(1, numSendBufs * sizeof(MPI_Request));
        }

        // Pack up nodal data for each remote domain to which this
        // domain will be sending segment data

        if (numSendBufs>0) 
            PackSegmentData(home, numSendBufs, segCommInfo, cellSegLists, cellSegCnts, sendBuf, sendBufEnts);

        // Have the current domain asynchronously send out the sizes of the
        // buffers it will be transmitting followed by the actual buffers.
        // MPI request structs for buffer sends are in <home>, while structs
        // for buffer length sends are dynamically allocated and returned
        // explicitly to the caller.
        // 
        // NOTE: We will not wait for the sends to complete right now.  That
        // will be handled when the local domain is ready to begin receiving data.

        for (int i=0; (i<numSendBufs); i++) 
            MPI_Isend(&sendBufEnts[i], 1, MPI_INT, segCommInfo[i].domID, MSG_SEGDATA_LEN, MPI_COMM_WORLD, &sendRequests[i]);

        for (int i=0; (i<numSendBufs); i++) 
            MPI_Isend(sendBuf[i], sendBufEnts[i], MPI_DOUBLE, segCommInfo[i].domID, MSG_SEGDATA, MPI_COMM_WORLD, &home->outRequests[i]);

        *sendBufs           = sendBuf;
        *sendBufLenRequests = sendRequests;
        *sendBufLengths     = sendBufEnts;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CommRecvSegForces
 *
 *      Description:  This function will accept incoming buffers of
 *                    segment force data, update the local nodes
 *                    as necessary, and release all communication
 *                    buffers/requests before returning to the caller.
 *                    Incoming buffers were sent via CommSendSegForces().
 *
 *      Arguments:
 *          IN:  numRecvBufs    Number of remote domains from which this
 *                              domain will be receiving segment data
 *          IN:  numSendBufs    Number of remote domains to which this
 *                              domain has sent segment data
 *          IN/OUT: sendBufs    Location of the pointer to array of
 *                              <numSendBufs> buffers associate with
 *                              buffers sent to the remote domains.
 *                              All buffers and array of pointers will
 *                              be freed and pointers NULL'ed before
 *                              returning to caller.
 *          IN/OUT: sendBufLengths Pointer to array of <numSendBufs> buffer
 *                              lengths associated with buffers sent to the
 *                              remote domains.  Array of will be freed
 *                              and pointer NULL'ed before returning to
 *                              caller.
 *          IN/OUT: sendBufLenRequests Pointer to array of <numSendBufs> 
 *                              MPI request buffers associated with buffer
 *                              length messages sent to the remote domains.
 *                              Array will be freed and pointer NULL'ed before
 *                              returning to caller.
 *                              
 *
 *-------------------------------------------------------------------------*/
void CommRecvSegForces(
   Home_t         *home              ,      ///< points to the home data structure
   int             numRecvBufs       ,      ///< number of remote domains receiving from
   int             numSendBufs       ,      ///< number of remote domains sending to
   real8        ***sendBufs          ,      ///< the send/recv buffers
   int           **sendBufLengths    ,      ///< lengths of the send/recv buffers
   MPI_Request   **sendBufLenRequests )     ///< mpi request buffers
{

    // Allocate arrays for incoming message lengths and pointer
    // to buffers for incoming messages.  Lengths specified as count
    // of real8 values being communicated

    real8       **recvBuf=0;
    int          *recvBufEnts=0;
    MPI_Request  *recvLenRequests=0;

    if (numRecvBufs>0) 
    {
        recvBuf         = (real8      **)calloc(1,numRecvBufs*sizeof(real8 *)    );
        recvBufEnts     = (int         *)calloc(1,numRecvBufs*sizeof(int)        );
        recvLenRequests = (MPI_Request *)calloc(1,numRecvBufs*sizeof(MPI_Request));
    }
    
    // Pre-issue receives of buffer lengths from any domain that will
    // be sending data to the current domain.

    for (int i=0; (i<numRecvBufs); i++) {
        MPI_Irecv(&recvBufEnts[i], 1, MPI_INT, MPI_ANY_SOURCE, MSG_SEGDATA_LEN, MPI_COMM_WORLD, &recvLenRequests[i]);
    }

    // Wait for the incoming buffer lengths, allocate buffers of the
    // appropriate size and post receives for the incoming buffers.

    for (int i=0; (i<numRecvBufs); i++) 
    {
        int         recvIndex=0;
        MPI_Status  recvLenStatus;

        MPI_Waitany(numRecvBufs, recvLenRequests, &recvIndex, &recvLenStatus);

        recvBuf[recvIndex] = (real8 *)malloc(recvBufEnts[recvIndex] * sizeof(real8));

        MPI_Irecv(recvBuf[recvIndex], recvBufEnts[recvIndex],
                  MPI_DOUBLE, recvLenStatus.MPI_SOURCE, MSG_SEGDATA,
                  MPI_COMM_WORLD, &home->inRequests[recvIndex]);
    }

    // If code is threaded, wait for all the buffers to arrive then
    // process them in parallel on different threads.  Otherwise, 
    // wait for any buffer to arrive and start processing the buffers
    // immediately upon receipt.

#ifdef _OPENMP
    if (numRecvBufs>0)
        MPI_Waitall(numRecvBufs, home->inRequests, home->inStatus);

#pragma omp parallel for
    for (int i=0; (i<numRecvBufs); i++) 
    {
        UnpackSegmentData(home, recvBuf[i]);
        free(recvBuf[i]);
    }
#else  // _OPENMP not defined

    for (int i=0; (i<numRecvBufs); i++) 
    {
        int  recvIndex=0;

        MPI_Waitany(numRecvBufs, home->inRequests, &recvIndex, home->inStatus);

        UnpackSegmentData(home, recvBuf[recvIndex]);

        free(recvBuf[recvIndex]);
    }
#endif

    // Verify that all the buffers (and buffer length messages) sent
    // by this domain have been received at their destinations before
    // we free the buffers.

    if (numSendBufs > 0) 
    {
        MPI_Waitall(numSendBufs,  home->outRequests , home->outStatus);
        MPI_Waitall(numSendBufs, *sendBufLenRequests, home->outStatus);
    }

    // Free the buffers...

    if (numSendBufs > 0) 
    {
        for (int i=0; (i<numSendBufs); i++) 
            if((*sendBufs)[i]) { free((*sendBufs)[i]); (*sendBufs)[i]=0; }

        if(*sendBufs          ) { free(*sendBufs          ); *sendBufs          =0; }
        if(*sendBufLenRequests) { free(*sendBufLenRequests); *sendBufLenRequests=0; }
        if(*sendBufLengths    ) { free(*sendBufLengths    ); *sendBufLengths    =0; }
    }

    if (numRecvBufs > 0) 
    {
        if(recvBuf        ) { free(recvBuf        ); recvBuf        =0; }
        if(recvBufEnts    ) { free(recvBufEnts    ); recvBufEnts    =0; }
        if(recvLenRequests) { free(recvLenRequests); recvLenRequests=0; }
    }
}

#endif  // end ifdef PARALLEL
