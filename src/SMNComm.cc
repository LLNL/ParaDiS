/**************************************************************************
 *
 *      Module:       SMNComm.c
 *      Description:  Contains functions for packing, sending, receiving and
 *                    unpacking the partial force data calculated during
 *                    the parallel multinode splitting process.
 *
 *      Public functions:
 *          SMNPackForceBuffers()
 *          SMNUnpackBuf()
 *          SMNForceComm()
 *
 *************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "SMN.h"

/*
 *      Define some tags for MPI messages
 */
#define SMN_LEN_TAG 1900
#define SMN_BUF_TAG 1901


/*------------------------------------------------------------------------
 *
 *      Function:     SMNPackForceBuffers
 *      Description:  For each domain owning a multinode for which this
 *                    domain (may have) calculated some force data during
 *                    a parallel multinode split, pack appropriate force
 *                    data into a buffer for the remote domain, and build
 *                    the list of domains to which buffers must be sent.
 *
 *      Arguments:
 *          in:  gmnCount    Number of ghost multi-nodes for which this
 *                           domain was responsible for some force calcs
 *          in:  gmnInfo     Array of <gmnCount> structs containing the
 *                           ghost multi-node info and the related force
 *                           data calculated locally
 *          out: numSendBufs Set to the number of remote domains to which
 *                           this domain needs to send force data
 *          out: sendDoms    Array of <numSendBufs> integers indicating the
 *                           remote domains to which buffers must be sent.
 *                           Caller is responsible for freeing this array.
 *          out: sendBufLens Array of <numSendBufs> integers indicating the
 *                           length of the corresponding buffer to be sent.
 *                           Caller is responsible for freeing this array.
 *          out: sendBufs    Array of <numSendBufs> pointers to buffers
 *                           to be sent.  Caller is responsible for freeing
 *                           each buffer in this array and the array itself.
 *
 *----------------------------------------------------------------------*/
void SMNPackForceBuffers(Home_t *home, int gmnCount, SMN_Info_t *gmnInfo,
                         int *numSendBufs, int **sendDoms, int **sendBufLens,
                         real8 ***sendBufs)
{
#ifdef PARALLEL
        int  i;
        int  sendDomCount;
        int  *sendDomList, *sendBufOffset, *sendLen;

/*
 *      Loop through the ghost multinode list and build array of
 *      remote domains to which data must be sent.
 */
        sendDomCount = gmnCount;

        if (sendDomCount > 0) {
            sendDomList = (int *)malloc(gmnCount * sizeof(int));
            sendLen = (int *)malloc(gmnCount * sizeof(int));
        } else {
            sendDomList = (int *)NULL;
            sendLen = (int *)NULL;
        }

        for (i = 0; i < gmnCount; i++) {
            sendDomList[i] = gmnInfo[i].multiNode->myTag.domainID;
            sendLen[i] = 0;
        }

        qsort(sendDomList, sendDomCount, sizeof(int), IntCompare);
        Uniq(sendDomList, &sendDomCount);

/*
 *      Loop through the ghost multinodes, calculate how much associated
 *      data needs to be sent to the remote domain, and increment the
 *      corresponding buffer length appropriately.
 */
        for (i = 0; i < gmnCount; i++) {
            int j, numVals, bufferID, remDomID;

            remDomID = gmnInfo[i].multiNode->myTag.domainID;

            for (bufferID = 0; bufferID < sendDomCount; bufferID++) {
                if (remDomID == sendDomList[bufferID]) {
                    break;
                }
            }

            numVals = 3;

/*
 *          We'll only pack segment force data if there were no errors
 *          errors splitting the ghost multinode or calculating forces.
 */
            if (gmnInfo[i].status == 0) {
                for (j = 0; j < gmnInfo[i].numSplits; j++) {
                    numVals += (gmnInfo[i].splitInfo[j].node1Segs * 6) +
                               (gmnInfo[i].splitInfo[j].node2Segs * 6) + 2;
                }
            }

            sendLen[bufferID] += numVals * sizeof(real8);
        }

/*
 *      Allocate the buffers to be sent to the remote domains then
 *      pack them with data.
 */
        if (sendDomCount > 0) {
            *sendBufs = (real8 **)malloc(sendDomCount * sizeof(real8 *));
            sendBufOffset = (int *)calloc(1, sendDomCount * sizeof(int));
        } else {
            *sendBufs = (real8 **)NULL;
            sendBufOffset = (int *)NULL;
        }

        for (i = 0; i < sendDomCount; i++) {
            (*sendBufs)[i] = (real8 *)malloc(sendLen[i] * sizeof(real8));
        }

        for (i = 0; i < gmnCount; i++) {
            int   j, offset;
            int   bufferID, remDomID;
            real8 *buf;

            remDomID = gmnInfo[i].multiNode->myTag.domainID;

            for (bufferID = 0; bufferID < sendDomCount; bufferID++) {
                if (remDomID == sendDomList[bufferID]) {
                    break;
                }
            }

            buf = (*sendBufs)[bufferID];
            offset = sendBufOffset[bufferID];

            buf[offset++] = (real8)gmnInfo[i].multiNode->myTag.index;
            buf[offset++] = (real8)gmnInfo[i].status;
            buf[offset++] = (real8)gmnInfo[i].numSplits;

/*
 *          Only pack segmenmt force data if there were no errors
 *          splitting the ghost multinode or calculating forces
 */
            if (gmnInfo[i].status != 0) {

                sendBufOffset[bufferID] = offset;
                continue;
            }

            for (j = 0; j < gmnInfo[i].numSplits; j++) {
                int         k;
                SMN_Split_t *splitInfo;

                splitInfo = &gmnInfo[i].splitInfo[j];

                buf[offset++] = (real8)splitInfo->node1Segs;

                for (k = 0; k < splitInfo->node1Segs; k++) {

                    buf[offset++] = splitInfo->node1Force.segForce[k].f1[0];
                    buf[offset++] = splitInfo->node1Force.segForce[k].f1[1];
                    buf[offset++] = splitInfo->node1Force.segForce[k].f1[2];

                    buf[offset++] = splitInfo->node1Force.segForce[k].f2[0];
                    buf[offset++] = splitInfo->node1Force.segForce[k].f2[1];
                    buf[offset++] = splitInfo->node1Force.segForce[k].f2[2];
                }

                buf[offset++] = (real8)splitInfo->node2Segs;

                for (k = 0; k < splitInfo->node2Segs; k++) {

                    buf[offset++] = splitInfo->node2Force.segForce[k].f1[0];
                    buf[offset++] = splitInfo->node2Force.segForce[k].f1[1];
                    buf[offset++] = splitInfo->node2Force.segForce[k].f1[2];

                    buf[offset++] = splitInfo->node2Force.segForce[k].f2[0];
                    buf[offset++] = splitInfo->node2Force.segForce[k].f2[1];
                    buf[offset++] = splitInfo->node2Force.segForce[k].f2[2];
                }

            }

            if ( (offset * (int) sizeof(real8)) > sendLen[bufferID]) {
                Fatal("Buffer overrun");
            }

            sendBufOffset[bufferID] = offset;

        }  /* end for (i = 0; i < gmnCount; ... ) */

/*
 *      Free up any temporary buffers and store values for caller
 */
        if (sendDomCount > 0) {
            free(sendBufOffset);
        }

        *numSendBufs = sendDomCount;
        *sendDoms = sendDomList;
        *sendBufLens = sendLen;

#endif /* ifdef PARALLEL */

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:     SMNUnpackBuf
 *      Description:  Unpack a buffer of partial forces calculated in
 *                    a remote domain.  These partial forces are
 *                    associated with various multi-node splitting
 *                    possibilities.  Forces from the buffer will be
 *                    summed into the local structures containing
 *                    the forces for each possible split of the the
 *                    local multi-nodes.
 *
 *      Arguments:
 *          in:     bufLen      Length of <buf> as a count of double precision
 *                              values.
 *          in:     buf         Buffer containing a portion of the force
 *                              data needed to complete the parallel multi-
 *                              node splitting process on the local 
 *                              multi-nodes.
 *          in:     lmnCount    Number of local multi-nodes 
 *          in/out: lmnInfo     Array of <lmnCount> structs containing the
 *                              local multi-node info.  These structs will
 *                              be updated with the force data in the buffer
 *                              from the remote domain.
 *
 *----------------------------------------------------------------------*/
void SMNUnpackBuf(Home_t *home, int bufLen, real8 *buf, int lmnCount,
                  SMN_Info_t *lmnInfo)
{
#ifdef PARALLEL
        int   i, numVals;
        int   tagIndex, status, numSplits;
        int   offset, mnIndex;


        offset = 0;
        numVals = bufLen / sizeof(double);

        while (offset < numVals) {

            tagIndex   = buf[offset++];
            status     = buf[offset++];
            numSplits  = buf[offset++];

/*
 *          Look up the multinode in the local multinode list.  If it's not
 *          found something is bad.
 */
            mnIndex = -1;

            for (i = 0; i < lmnCount; i++) {
                if (tagIndex == lmnInfo[i].multiNode->myTag.index) {
                    mnIndex = i;
                    break;
                }
            }

            if (mnIndex < 0) {
                printf("Task %d: received buffer for unknown node (%d,%d)\n",
                       home->myDomain, home->myDomain, tagIndex);
                Fatal("No local multinode matching tag from remote domain");
            }

/*
 *          If the status from the remote domain is non-zero, there was
 *          some sort of problem, so just set the local status flag for
 *          the multinode so we don't bother with it. 
 */
            if (status != 0) {
                lmnInfo[mnIndex].status = status;
                continue;
            }

/*
 *          Pull the force data from the buffer and update the local
 *          multinode force data.
 */
            if (numSplits != lmnInfo[mnIndex].numSplits) {
                Fatal("Mismatch on splitCount");
            }

            for (i = 0; i < numSplits; i++) {
                int         j, segCount;
                SMN_Split_t *splitInfo;
                SegForce_t  *fSeg;

                splitInfo = &lmnInfo[mnIndex].splitInfo[i];

                segCount = buf[offset++];

/* FIX ME! This may be an obsolete check.  I think with recent changes 
 *         we should always know the node1Segs and node2Segs values at
 *         this point?
 */
 
                if (splitInfo->node1Segs == 0) {
                    splitInfo->node1Segs = segCount;
                } else if (segCount != splitInfo->node1Segs) {
                    Fatal("Mismatch on node1 segment count");
                }

                for (j = 0; j < splitInfo->node1Segs; j++) {

                    fSeg = &splitInfo->node1Force.segForce[j];

                    fSeg->f1[0] += buf[offset++];
                    fSeg->f1[1] += buf[offset++];
                    fSeg->f1[2] += buf[offset++];

                    fSeg->f2[0] += buf[offset++];
                    fSeg->f2[1] += buf[offset++];
                    fSeg->f2[2] += buf[offset++];
                }

                segCount = buf[offset++];

                if (splitInfo->node2Segs == 0) {
                    splitInfo->node2Segs = segCount;
                } else if (segCount != splitInfo->node2Segs) {
                    Fatal("Mismatch on node2 segment count");
                }

                for (j = 0; j < splitInfo->node2Segs; j++) {

                    fSeg = &splitInfo->node2Force.segForce[j];

                    fSeg->f1[0] += buf[offset++];
                    fSeg->f1[1] += buf[offset++];
                    fSeg->f1[2] += buf[offset++];

                    fSeg->f2[0] += buf[offset++];
                    fSeg->f2[1] += buf[offset++];
                    fSeg->f2[2] += buf[offset++];
                }
            }
        }

#endif  /* end ifdef PARALLEL */

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:     SMNForceComm
 *      Description:  This function controls the sending/receiving of
 *                    partial force data needed during the parallel
 *                    split multi-node process.
 *
 *      Arguments:
 *          in:  sendDomCount Number of partial force buffers to be sent
 *          in:  sendDoms     Array of <sendDomCount> send buffers
 *          in:  sendBufLen   Array of <sendDomCount> send buffer lengths
 *          in:  sendBuf      Array of <sendDomCount> pointers to buffers
 *                            this domain must send to remote domains
 *          in:  recvDomCount Number of remote domains from which this 
 *                            domain expects to receive partial forces
 *          in:  recvDoms     Array of <recvDoms> domain ID's.
 *          in:  lmnCount     Number of local multi-nodes 
 *          in:  lmnInfo      Array of <lmnCount> structs containing the
 *                            local multi-node info and the related force
 *                            data
 *
 *----------------------------------------------------------------------*/
void SMNForceComm(Home_t *home, int sendDomCount, int *sendDoms,
                  int *sendBufLen, real8 **sendBuf, int recvDomCount,
                  int *recvDoms, int lmnCount, SMN_Info_t *lmnInfo)
{
#ifdef PARALLEL
        int          i, valCount;
        int          *recvBufLen=0;
        real8        **recvBuf=0;
        MPI_Request  *sendRequest=0, *recvRequest=0;
        MPI_Status   *sendStatus=0, recvStatus;


        if (sendDomCount > 0) {
            sendStatus = (MPI_Status *)malloc(sendDomCount*sizeof(MPI_Status));
            sendRequest = (MPI_Request *)malloc(sendDomCount *
                                                sizeof(MPI_Request));
        }

        if (recvDomCount > 0) {
            recvRequest = (MPI_Request *)malloc(recvDomCount *
                                                sizeof(MPI_Request));
            recvBufLen = (int *)malloc(recvDomCount * sizeof(int));
            recvBuf = (real8 **)malloc(recvDomCount * sizeof(real8 *));
        }

/*
 *      Post receives of message lengths from each remote domain
 *      that should be sending force data to this domain
 */
        for (i = 0; i < recvDomCount; i++) {
            MPI_Irecv(&recvBufLen[i], 1, MPI_INT, recvDoms[i], SMN_LEN_TAG,
                      MPI_COMM_WORLD, &recvRequest[i]);
        }

/*
 *      Send buffer lengths out to the remote domains
 */
        for (i = 0; i < sendDomCount; i++) {
            MPI_Isend(&sendBufLen[i], 1, MPI_INT, sendDoms[i], SMN_LEN_TAG,
                      MPI_COMM_WORLD, &sendRequest[i]);
        }

/*
 *      Read the buffer size from each remote domain, allocate an
 *      appropriately sized buffer, then post receives for the incoming
 *      buffers
 */
        for (i = 0; i < recvDomCount; i++) {
            int recvIndex;
            MPI_Waitany(recvDomCount, recvRequest, &recvIndex, &recvStatus);
            recvBuf[recvIndex] = (real8 *)malloc(recvBufLen[recvIndex]);
        }

        for (i = 0; i < recvDomCount; i++) {
            valCount = recvBufLen[i] / sizeof(real8);
            MPI_Irecv(recvBuf[i], valCount, MPI_DOUBLE, recvDoms[i],
                      SMN_BUF_TAG, MPI_COMM_WORLD, &recvRequest[i]);
        }

/*
 *      Wait for all the buffer length sends to complete, then send the
 *      data buffers themselves.
 */
        MPI_Waitall(sendDomCount, sendRequest, sendStatus);

        for (i = 0; i < sendDomCount; i++) {
            valCount = sendBufLen[i] / sizeof(real8);
            MPI_Isend(sendBuf[i], valCount, MPI_DOUBLE, sendDoms[i],
                      SMN_BUF_TAG, MPI_COMM_WORLD, &sendRequest[i]);
        }

/*
 *      Read and process the data buffers as they arrive
 */
        for (i = 0; i < recvDomCount; i++) {
            int recvIndex;
            MPI_Waitany(recvDomCount, recvRequest, &recvIndex, &recvStatus);
            SMNUnpackBuf(home, recvBufLen[recvIndex], recvBuf[recvIndex],
                         lmnCount, lmnInfo);
            free(recvBuf[recvIndex]);
        }

/*
 *      Wait for all the data buffer sends to complete then clean
 *      up and we're done.
 */
        MPI_Waitall(sendDomCount, sendRequest, sendStatus);

        if (sendDomCount > 0) {
            free(sendRequest);
            free(sendStatus);
        }

        if (recvDomCount > 0) {
            free(recvRequest);
            free(recvBufLen);
            free(recvBuf);
        }

#endif  /* end ifdef PARALLEL */

        return;
}
