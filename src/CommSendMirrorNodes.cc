/*
 *    CommSendMirrorNodes.c:  This module contains the functions
 *        needed for downloading the node data from the
 *        remote domains (aka. mirror domains) to processor 0
 *        during X-window plotting.  In order to prevent the
 *        exhaustion of memory on processor zero, the data is
 *        downloaded one remote domain at a time.  The data for
 *        a remote domain is transferred asynchronously while
 *        processor 0 is processing the data from the previous
 *        domain to help hide the latency of the data transfers.
 *
 *    Included functions:
 *        CommSendMirrorNodes()
 *        GenDomainOutput()
 *        GetNbrCoords()
 *        PackMirrorNodes()
 *        UnpackMirrorNodes()
 */

#include "mpi_portability.h"

#include "Home.h"
#include "Node.h"
#include "Comm.h"
#include "Util.h"
#include "QueueOps.h"
#include "MirrorDomain.h"

static void PackMirrorNodes(Home_t *home, int *numNodes);
static void UnpackMirrorNodes(Home_t *home, int domIndex);

/*-------------------------------------------------------------------------
 * Function:    GetNbrCoord
 * Description: Retrieve the x,y,z coordinates for the neighbor
 *              node at the end of the specified arm of <node>.
 *              This subroutine should only be invoked during
 *              the output generation process for nodes in a
 *              domain whose data is currently in memory on processor
 *              0 because it requires arrays in the mirror domain
 *              data that will only be present while the data for
 *              that domain is loaded on procesor 0.
 * Args:
 *        node
 *        arm
 *        x,y,z
 *
 *------------------------------------------------------------------------*/
void GetNbrCoords(Home_t *home, Node_t *node, int arm, real8 *x, real8 *y, real8 *z)
{
    // If <node> is native to domain zero, the neighboring node's
    // data is either native, or accessible in as a ghost node, so
    // grab the coords from the node structure.
  
    // If <node> is not a native node, the necessary coords should
    // be in the mirror domain's arm coordinate arrays.

    if (node->myTag.domainID==0) 
    {
        Node_t *nbrNode = GetNodeFromTag(home, node->nbrTag[arm]);

        if (nbrNode)
        {
            *x = nbrNode->x;
            *y = nbrNode->y;
            *z = nbrNode->z;
        }
    } 
    else 
    {
        MirrorDomain_t *domain = home->mirrorDomainKeys[node->myTag.domainID];

        if (domain)
        {
            *x = domain->armX[node->armCoordIndex[arm]];
            *y = domain->armY[node->armCoordIndex[arm]];
            *z = domain->armZ[node->armCoordIndex[arm]];
        }
    }

    return;
}


/*-------------------------------------------------------------------------
 * Function:    GenDomainOutput
 * Description: Invokes X-Window plotting appropriate to generate
 *              output data from the node information for a 
 *              specified domain at the given stage.
 *
 * Arguments:
 *     stage     indicates the current stage of
 *               the program execution.  This is
 *               used as an additional check since
 *               not all output types are to be
 *               generated at all stages. Valid
 *               values are:
 *                        STAGE_INIT
 *                        STAGE_CYCLE
 *                        STAGE_TERM
 *
 *     domIndex  id of the domain whose data is
 *               to be processed
 *     blkFlag   flag to be passed to the output
 *               functions so they know whether
 *               they need to do the specialized
 *               processing necessary for either
 *               the first or last blocks of data
 *
 *------------------------------------------------------------------------*/
static void GenDomainOutput(Home_t *home, int stage, int domIndex, int blkFlag)
{
#ifndef NO_XWINDOW

    switch (stage) 
    {
        case STAGE_CYCLE:
            TimerStart(home, PLOT);
            Plot(home, domIndex, blkFlag);
            TimerStop(home, PLOT);
#ifdef NO_THREAD
            WinEvolve();
#endif
            break;

        case STAGE_INIT:
            Plot(home, domIndex, blkFlag);
            break;
    }

#endif
    return;
}


/*-------------------------------------------------------------------------
 *
 *    Function:     CommSendMirrorNodes
 *    Description:    Download the data from all remote (mirror) domains
 *            to processor zero, and make calls to generate any
 *            appropriate output data for each block of data.
 *
 *------------------------------------------------------------------------*/
void CommSendMirrorNodes(Home_t *home, int stage)
{
    int    token=5319009;

    // Have each process package up its node data and get count
    // of nodes, segments and arms in the domain.  (For process 0
    // this call does not package up any data, but will return the
    // correct counts.
    
    int numNodes=0;
    PackMirrorNodes(home, &numNodes);

    int maxNodes   = numNodes;
    int maxPackSiz = 0;

    // Get the maximum buffer size of all the mirror domains as well
    // as the maximum number of nodes in any single domain.

#ifdef PARALLEL
    int    outVals[3], inVals[3];

    outVals[0] = home->maxPackSiz;
    outVals[1] = maxNodes;

    MPI_Reduce(outVals, inVals, 2, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    maxPackSiz = inVals[0];
    maxNodes   = inVals[1];
#else
    maxPackSiz = 0;
    maxNodes   = maxNodes;
#endif

    // processor 0 now does most of the work; prompting the remote domains
    // for their data, unpacking the data, and calling functions to generate
    // appropriate outputs

    if (home->myDomain==0) 
    {
        // Initiate the first of the data transfers from the
        // remote domains (if there are any remote domains).

#ifdef PARALLEL
        home->maxPackSiz = maxPackSiz;

        if (home->numDomains > 1) 
        {
            home->inBuf = (char *)malloc(maxPackSiz);

            MPI_Irecv(home->inBuf, maxPackSiz, MPI_PACKED, 1, MSG_MIRRORS, MPI_COMM_WORLD, home->inRequests);
            MPI_Send (&token, 1, MPI_INT, 1, MSG_SEND_MIRROR, MPI_COMM_WORLD);
        }
#endif

        // Process nodes from domain zero while the next domain
        // is sending its data.

        int blkFlag = FIRST_BLOCK;
        if (home->numDomains==1) blkFlag |= LAST_BLOCK;
        GenDomainOutput(home, stage, 0, blkFlag);

        // Allocate and initialize the mirror domains key array
        // and start processing data from remote domains.

        home->mirrorDomainKeys = (MirrorDomain_t **) calloc(1, home->numDomains * sizeof(MirrorDomain_t *) );

#ifdef PARALLEL
        for (int domIndex = 1; domIndex < home->numDomains; domIndex++) 
        {
            // Wait for the previous transfer to complete.  (If we're
            // lucky, the transfer completed while the previous
            // domain's data was being processed, and there's no
            // delay.

            MPI_Wait(home->inRequests, home->inStatus);

            if  (domIndex==(home->numDomains-1)) { blkFlag = LAST_BLOCK; }
            else                                 { blkFlag = 0; }
            
            UnpackMirrorNodes(home, domIndex);

            // Notify the next domain to send its data so the
            // transfer takes place while we process the 
            // buffer we just got

            if (domIndex<(home->numDomains-1))
            {
                MPI_Irecv(home->inBuf, maxPackSiz, MPI_PACKED, domIndex+1, MSG_MIRRORS, MPI_COMM_WORLD, home->inRequests);
                MPI_Send (&token, 1, MPI_INT, domIndex+1, MSG_SEND_MIRROR, MPI_COMM_WORLD);
            }

            // Handle the buffer we just received

            GenDomainOutput  (home, stage, domIndex, blkFlag);
            MirrorDomain_Free(home, domIndex);
        }
#endif

        if (home->numDomains!=1) {
            free(home->inBuf);
        }

        free(home->mirrorDomainKeys);
        home->mirrorDomainKeys=0;
    } 
    else 
    {
        // All tasks other than task zero wait for a prompt from
        // task zero then send their data off.

#ifdef PARALLEL
        MPI_Status status;

        MPI_Recv (&token, 1, MPI_INT, 0, MSG_SEND_MIRROR, MPI_COMM_WORLD, &status);
        MPI_Isend(home->outBuf, home->maxPackSiz, MPI_PACKED, 0, MSG_MIRRORS, MPI_COMM_WORLD, &home->outRequests[0]);
        MPI_Wait (home->outRequests, home->outStatus);

        free(home->outBuf);
#endif
    }

    return;
}


/*-------------------------------------------------------------------------
 * Function:    PackMirrorNodes
 * Description: Collect all nodal data on this domain into an
 *              MPI_PACKED buffer to send to processor 0.  Also
 *              obtain node, segment and arm counts and store
 *              them in the caller provided variables.
 *------------------------------------------------------------------------*/

static void PackMirrorNodes(Home_t *home, int *numNodes)
{

    // Count the number of nodes and the total number of arms

    int nodeCount=0, armCount=0, incCount=0;

    for (int i=0; (i<home->newNodeKeyPtr); i++) 
    {
        Node_t *node = home->nodeKeys[i];

        if (node)
        {
           nodeCount++;
           armCount += node->numNbrs;
        }
    }

    *numNodes = nodeCount;

      // If this is domain 0, there's no need to do anything else.
      // Otherwise, allocate buffers and package up the data to be sent
      // to domain 0.

    if (home->myDomain == 0) 
    {
#ifdef PARALLEL
        home->maxPackSiz = 0;
#endif
        return;
    }

#ifdef PARALLEL

#ifdef ESHELBY
    incCount = home->locInclusionCount;
#endif

    int intCount = nodeCount * INTS_PER_MIRROR_NODE + 
                   armCount  * INTS_PER_MIRROR_ARM  + EXTRA_MIRROR_INTS;

    int fltCount = nodeCount * FLTS_PER_MIRROR_NODE +
                   armCount  * FLTS_PER_MIRROR_ARM  +
                   incCount  * FLTS_PER_MIRROR_INCLUSION;

    int   *intBuf = (int   *) malloc(intCount * sizeof(int  ));
    real8 *fltBuf = (real8 *) malloc(fltCount * sizeof(real8));

    int iIdx = 0;
    int fIdx = 0;

    intBuf[iIdx++] = nodeCount;
    intBuf[iIdx++] = home->newNodeKeyPtr;
    intBuf[iIdx++] = intCount;
    intBuf[iIdx++] = fltCount;
    intBuf[iIdx++] = armCount;
    intBuf[iIdx++] = incCount;

    for (int i=0; (i<home->newNodeKeyPtr); i++) 
    {
        Node_t *node = home->nodeKeys[i];

        if (!node) continue;

        intBuf[iIdx++] = node->myTag.index;
        intBuf[iIdx++] = node->constraint;
        intBuf[iIdx++] = node->numNbrs;

        for (int j=0; (j<node->numNbrs); j++) 
        {
            intBuf[iIdx++] = node->nbrTag[j].domainID;
            intBuf[iIdx++] = node->nbrTag[j].index;

            fltBuf[fIdx++] = node->burgX[j];
            fltBuf[fIdx++] = node->burgY[j];
            fltBuf[fIdx++] = node->burgZ[j];

            fltBuf[fIdx++] = node->nx[j];
            fltBuf[fIdx++] = node->ny[j];
            fltBuf[fIdx++] = node->nz[j];

            // Various plotting functions require the x,y,z coords of
            // nodes at the far end of arms.  Since we no longer
            // download the entire problem space to task 0 in one
            // shot, each domain now needs to send with the node
            // data, the x,y,z coords of each of the node's neighbors
            // so they will be available to the plotting functions.

            Node_t *nbrNode = GetNodeFromTag(home, node->nbrTag[j]);

            if (nbrNode == (Node_t *)NULL) {
                Fatal("Task %d: %s, line %d -- error looking up node (%d,%d)",
                      home->myDomain, __FILE__, __LINE__,
                      node->nbrTag[j].domainID, node->nbrTag[j].index);
            }

            fltBuf[fIdx++] = nbrNode->x;
            fltBuf[fIdx++] = nbrNode->y;
            fltBuf[fIdx++] = nbrNode->z;
        }

        fltBuf[fIdx++] = node->x;
        fltBuf[fIdx++] = node->y;
        fltBuf[fIdx++] = node->z;
        fltBuf[fIdx++] = node->vX;
        fltBuf[fIdx++] = node->vY;
        fltBuf[fIdx++] = node->vZ;
    }

#ifdef ESHELBY
    // If there are any inclusions to send, pack the info
    // needed for plotting them.

    for (int i=0; (i<incCount); i++) 
    {
        EInclusion_t *inclusion = &home->eshelbyInclusions[i];

        if (inclusion)
        {
            fltBuf[fIdx++] = (real8) inclusion->id;
            fltBuf[fIdx++] = inclusion->radius[X];
            fltBuf[fIdx++] = inclusion->radius[Y];
            fltBuf[fIdx++] = inclusion->radius[Z];
            fltBuf[fIdx++] = inclusion->position[X];
            fltBuf[fIdx++] = inclusion->position[Y];
            fltBuf[fIdx++] = inclusion->position[Z];
            fltBuf[fIdx++] = inclusion->rotation[0][X];
            fltBuf[fIdx++] = inclusion->rotation[0][Y];
            fltBuf[fIdx++] = inclusion->rotation[0][Z];
            fltBuf[fIdx++] = inclusion->rotation[1][X];
            fltBuf[fIdx++] = inclusion->rotation[1][Y];
            fltBuf[fIdx++] = inclusion->rotation[1][Z];
        }
    }
#endif  // ESHELBY

    // Do a couple quick sanity checks on the amount of data stored
    // then package it all up for processor 0.

    if (iIdx > intCount) Fatal("PackMirrorNodes: intBuf mismatch");
    if (fIdx > fltCount) Fatal("PackMirrorNodes: fltBuf mismatch");

    int packedIntSiz=0;
    int packedFltSiz=0;

    MPI_Pack_size(intCount, MPI_INT   , MPI_COMM_WORLD, &packedIntSiz);
    MPI_Pack_size(fltCount, MPI_DOUBLE, MPI_COMM_WORLD, &packedFltSiz);

    int packSiz = packedIntSiz + packedFltSiz;

    home->outBuf = (char *) malloc(packSiz);

    int ipos = 0;
    MPI_Pack(intBuf, intCount, MPI_INT   , home->outBuf, packSiz, &ipos, MPI_COMM_WORLD);
    MPI_Pack(fltBuf, fltCount, MPI_DOUBLE, home->outBuf, packSiz, &ipos, MPI_COMM_WORLD);

    free(intBuf);
    free(fltBuf);

    home->maxPackSiz = packSiz;

#endif  // PARALLEL 

}


#ifdef PARALLEL
/*-------------------------------------------------------------------------
 * Function:    UnpackMirrorNodes
 * Description: Unpack all nodal data from the specified domain
 *              into the corresponding mirrorDomains element.
 *------------------------------------------------------------------------*/

static void UnpackMirrorNodes(Home_t *home, int domIndex) 
{
    // Allocate buffers large enough for this domain's data,
    // unpack the integer and float arrays, and pull out the
    // control words that indicate exact counts of data sent.

    int   *intBuf = (int   *) malloc(home->maxPackSiz);
    real8 *fltBuf = (real8 *) malloc(home->maxPackSiz);

    int    ipos=0, control[6];
    MPI_Unpack(home->inBuf, home->maxPackSiz, &ipos, control, 6, MPI_INT, MPI_COMM_WORLD);

    int numNodes      = control[0];
    int newNodeKeyPtr = control[1];
    int intCount      = control[2];
    int fltCount      = control[3];
    int armCount      = control[4];
#ifdef ESHELBY
    int incCount      = control[5];
#endif

    home->mirrorDomainKeys[domIndex] = (MirrorDomain_t *) malloc( sizeof(MirrorDomain_t) );

    MirrorDomain_t *mirrorDomain = home->mirrorDomainKeys[domIndex];

#ifdef ESHELBY
    mirrorDomain->inclusionCount = incCount;
#endif

    // Allocate space for the x,y,z coordinates associated with the neighbor
    // node at the end of each arm.  Needed for various plotting functions

    mirrorDomain->armX = (real8 *)malloc(sizeof(real8) * armCount);
    mirrorDomain->armY = (real8 *)malloc(sizeof(real8) * armCount);
    mirrorDomain->armZ = (real8 *)malloc(sizeof(real8) * armCount);

    mirrorDomain->newNodeKeyPtr = newNodeKeyPtr;
    mirrorDomain->nodeKeys      = (Node_t **) malloc( newNodeKeyPtr * sizeof(Node_t *) );

#ifdef ESHELBY
    if (incCount > 0) 
    {
        mirrorDomain->inclusionList = (EInclusion_t *)calloc(1, incCount * sizeof(EInclusion_t));
    }
#endif  // ESHELBY

    for (int i=0; (i<mirrorDomain->newNodeKeyPtr); i++)
        mirrorDomain->nodeKeys[i] = 0;

    intCount -= 6;  /* decrement by control words aleady read */

    MPI_Unpack(home->inBuf, home->maxPackSiz, &ipos, intBuf, intCount, MPI_INT   , MPI_COMM_WORLD);
    MPI_Unpack(home->inBuf, home->maxPackSiz, &ipos, fltBuf, fltCount, MPI_DOUBLE, MPI_COMM_WORLD);

    // Loop over the specified node count pulling out all data from
    // each buffer for the node.

    int iIdx=0;
    int fIdx=0;
    int coordIndex=0;

    for(int i=0; (i<numNodes); i++) 
    {
        Node_t *node    = PopFreeNodeQ(home);
        int     numNbrs = 0;

        node->myTag.domainID     = domIndex;
        node->myTag.index        = intBuf[iIdx++];
        node->constraint         = intBuf[iIdx++];
        numNbrs                  = intBuf[iIdx++];

        AllocNodeArms(node, numNbrs);

        node->armCoordIndex = (int *)malloc(sizeof(int) * numNbrs);

        for (int j=0; (j<numNbrs); j++) 
        {
            node->nbrTag[j].domainID = intBuf[iIdx++];
            node->nbrTag[j].index    = intBuf[iIdx++];

            node->burgX[j]=fltBuf[fIdx++];
            node->burgY[j]=fltBuf[fIdx++];
            node->burgZ[j]=fltBuf[fIdx++];

            node->nx[j] = fltBuf[fIdx++];
            node->ny[j] = fltBuf[fIdx++];
            node->nz[j] = fltBuf[fIdx++];

            mirrorDomain->armX[coordIndex] = fltBuf[fIdx++];
            mirrorDomain->armY[coordIndex] = fltBuf[fIdx++];
            mirrorDomain->armZ[coordIndex] = fltBuf[fIdx++];

            node->armCoordIndex[j] = coordIndex++;
        }

        node->x  = fltBuf[fIdx++];
        node->y  = fltBuf[fIdx++];
        node->z  = fltBuf[fIdx++];
        node->vX = fltBuf[fIdx++];
        node->vY = fltBuf[fIdx++];
        node->vZ = fltBuf[fIdx++];

        node->native = 0;
            
        mirrorDomain->nodeKeys[node->myTag.index] = node;
    }
    
#ifdef ESHELBY
    // If any eshelby inclusions were sent, unpack them now.

    for (int i=0; (i<incCount); i++) 
    {
        EInclusion_t *inclusion = &mirrorDomain->inclusionList[i];

        if (inclusion)
        {
            inclusion->id             = (int) fltBuf[fIdx++];
            inclusion->radius[X]      =       fltBuf[fIdx++];
            inclusion->radius[Y]      =       fltBuf[fIdx++];
            inclusion->radius[Z]      =       fltBuf[fIdx++];
            inclusion->position[X]    =       fltBuf[fIdx++];
            inclusion->position[Y]    =       fltBuf[fIdx++];
            inclusion->position[Z]    =       fltBuf[fIdx++];
            inclusion->rotation[0][X] =       fltBuf[fIdx++];
            inclusion->rotation[0][Y] =       fltBuf[fIdx++];
            inclusion->rotation[0][Z] =       fltBuf[fIdx++];
            inclusion->rotation[1][X] =       fltBuf[fIdx++];
            inclusion->rotation[1][Y] =       fltBuf[fIdx++];
            inclusion->rotation[1][Z] =       fltBuf[fIdx++];
        }
    }
#endif  // ESHELBY

    free(intBuf);
    free(fltBuf);

    if (iIdx > intCount) Fatal("UnpackMirrorNodes: intBuf mismatch");
    if (fIdx > fltCount) Fatal("UnpackMirrorNodes: fltBuf mismatch");
}
#endif  // #ifdef PARALLEL


