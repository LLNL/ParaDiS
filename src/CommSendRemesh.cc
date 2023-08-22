/****************************************************************************
 *
 *      Function:       CommSendRemesh
 *      Description:    Send to remote domains the opList containing
 *                      all topological operations performed by the
 *                      local domain that the remote domains need
 *                      to know about, and receive opList's from the
 *                      neighboring domains.  This is currently used
 *                      twice, after collision handling, then after
 *                      remeshing.
 *
 *****************************************************************************/
#include <stdio.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"

void CommSendRemesh(Home_t *home) 
{
        int            isrc, domainIdx, idst, outMsgLen;
        int            numRemoteDomains;
        int            *sendBufLen=0;
        char           **sendBuf=0;
        RemoteDomain_t *remDom=0;


        numRemoteDomains = home->remoteDomainCount;
        outMsgLen = home->opListUsedLen;

#ifdef PARALLEL
/*
 *      pre-issue receives of message length from each neighbor
 */
        for (isrc = 0; isrc < home->remoteDomainCount; isrc++) {
            domainIdx = home->remoteDomains[isrc];
            remDom = home->remoteDomainKeys[domainIdx];

            MPI_Irecv(&remDom->inBufLen, 1, MPI_INT, domainIdx, MSG_REMESH_LEN,
                      MPI_COMM_WORLD, &home->inRequests[isrc]);
        }

/*
 *      We're sending the same buffer to all remote domains, but
 *      MPI spec says we can't re-use the same buffer for multiple
 *      asynch sends at the same time, so need to make multiple
 *      copies of the buffer.  (Use the existing buffer as the
 *      one to be sent to the first remote domain)
 */
        if (numRemoteDomains > 0) {
            int i;

            sendBuf = (char **)malloc(numRemoteDomains * sizeof(char *));
            sendBufLen = (int *)malloc(numRemoteDomains * sizeof(int));

            sendBuf[0] = (char *)home->opListBuf;
            sendBufLen[0] = outMsgLen;

/*
 *          For first remote domain the existing buffer (opList) is used,
 *          so just need to set up remaining send buffers here.
 */
            for (i = 1; i < numRemoteDomains; i++) {
                sendBuf[i] = (char *)malloc(outMsgLen);
                sendBufLen[i] = outMsgLen;
                memcpy(sendBuf[i], home->opListBuf, outMsgLen);
            }
        }

/*
 *      Send the message length to the receiving neighbors
 */

        for (idst = 0; idst < numRemoteDomains; idst++) {
            domainIdx = home->remoteDomains[idst];
            remDom = home->remoteDomainKeys[domainIdx];

            MPI_Isend(&sendBufLen[idst], 1, MPI_INT, domainIdx, MSG_REMESH_LEN,
                      MPI_COMM_WORLD, &home->outRequests[idst]);
        } /* end for (idst = 0; ...) */

/*
 *      Wait for the length sends/receives to complete
 */
        MPI_Waitall(numRemoteDomains, home->outRequests, home->outStatus);
        MPI_Waitall(numRemoteDomains, home->inRequests, home->inStatus);

/*
 *      Allocate input buffers and pre-issue receives
 */
        for (isrc = 0; isrc < numRemoteDomains; isrc++) {
            domainIdx = home->remoteDomains[isrc];
            remDom = home->remoteDomainKeys[domainIdx];

            remDom->inBuf = (char *) malloc(remDom->inBufLen);
            MPI_Irecv(remDom->inBuf, remDom->inBufLen, MPI_BYTE, domainIdx, 
                      MSG_REMESH, MPI_COMM_WORLD, &home->inRequests[isrc]);
        }

/*
 *      Send the data
 */
        for (idst = 0; idst < home->remoteDomainCount; idst++) {
            domainIdx = home->remoteDomains[idst];
            remDom = home->remoteDomainKeys[domainIdx];

            MPI_Isend(sendBuf[idst], sendBufLen[idst], MPI_BYTE, domainIdx, 
                      MSG_REMESH, MPI_COMM_WORLD, &home->outRequests[idst]);
        } /* end for (idst = 0; ...) */

/*
 *      Wait for all traffic to complete
 */
        MPI_Waitall(numRemoteDomains, home->outRequests, home->outStatus);
        MPI_Waitall(numRemoteDomains, home->inRequests, home->inStatus);

/*
 *      Free any temporary buffers.  Note: first sendBuf points to the 
 *      actual <opList> which gets cleared elsewhere in the code, so
 *      don't free that one here.
 */
        if (numRemoteDomains > 0) {
            int i;

            for (i = 1; i < numRemoteDomains; i++) {
                free(sendBuf[i]);
            }

            free(sendBuf);
            free(sendBufLen);
        }
#endif
        return;
}
