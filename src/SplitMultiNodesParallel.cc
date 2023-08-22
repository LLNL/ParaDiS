/**************************************************************************
 *
 *      Module:       SplitMultiNodesParallel.c
 *      Description:  Main function to parallelize multi-node
 *                    splitting force calculations across multiple domains.
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "SMN.h"


void SplitMultiNodesParallel(Home_t *home)
{
        int        i, bufID;
        int        gmnCount, lmnCount;
        int        sendDomCount;
        int        *sendDoms;
        int        *sendBufLens;
        int        recvDomCount;
        int        *recvDoms;
        real8      **sendBufs;
        SMN_Info_t *gmnInfo, *lmnInfo;


/*
 *      Only do multinode splits on cycles that are multiples
 *      of the <splitMultiNodeFreq> control file parameter.
 */
        if (home->cycle % home->param->splitMultiNodeFreq) {
            return;
        }

        TimerStart(home, SPLIT_MULTI_NODES);

/*
 *      Get the lists of ghost and local multinodes for which this
 *      domain needs to do some force calcs, split each of those nodes
 *      and get the partial forces
 */
        gmnCount = 0;
        gmnInfo = (SMN_Info_t *)NULL;

        lmnCount = 0;
        lmnInfo = (SMN_Info_t *)NULL;

        BuildGhostMultiNodeList(home, &gmnCount, &gmnInfo);
        SMNEvalSplitForces(home, gmnCount, gmnInfo);

        BuildLocalMultiNodeList(home, &lmnCount, &lmnInfo);
        SMNEvalSplitForces(home, lmnCount, lmnInfo);

/*
 *      Pack the force data up for the appropriate remote domains
 */
        SMNPackForceBuffers(home, gmnCount, gmnInfo, &sendDomCount,
                            &sendDoms, &sendBufLens, &sendBufs);

/*
 *      Determine which remote domains should be sending force data
 *      to this one, send force data for ghost multinodes to the
 *      domains owning them, and receive force from any domain that may
 *      have calculated partial forces for a multinode local to this domain.
 */
        BuildRecvDomList(home, lmnCount, lmnInfo, &recvDomCount,
                         &recvDoms);

        SMNForceComm(home, sendDomCount, sendDoms, sendBufLens, sendBufs,
                     recvDomCount, recvDoms, lmnCount, lmnInfo);

/*
 *      Free various buffers allocated for the communication step
 */
        if (gmnCount > 0) {

            for (i = 0; i < gmnCount; i++) {
                free(gmnInfo[i].splitInfo);
            }

            free(gmnInfo);
        }

        if (sendDomCount > 0) {

            free(sendDoms);
            free(sendBufLens);

            for (bufID = 0; bufID < sendDomCount; bufID++) {
                free(sendBufs[bufID]);
            }

            free(sendBufs);
        }

        if (recvDomCount > 0) {
            free(recvDoms);
        }

/*
 *      We should have full forces for all possible splits of
 *      each local multi-node, so now do the final evaluation to
 *      determine if/how each multi-node should be split.
 */
        SMNEvalMultiNodeSplits(home, lmnCount, lmnInfo);

/*
 *      Free up the multinode lists...
 */
        if (lmnCount > 0) {

            for (i = 0; i < lmnCount; i++) {
                free(lmnInfo[i].splitInfo);
            }

            free(lmnInfo);
        }

        TimerStop(home, SPLIT_MULTI_NODES);

        return;
}
