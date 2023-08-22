/**************************************************************************
 *
 *      Module:       Migrate.c
 *      Description:  This module contains a driver function 
 *                    to control the process of identifying nodes
 *                    whose ownership must be transferred from one
 *                    domain to another and coordinating the transfer
 *                    of ownership.  Also contained are all
 *                    support functions unique to this capability.
 *
 *      Includes public functions
 *          Migrate()
 *
 *      Includes private functions
 *          BuildMigLists()
 *          CommSendMigrators()
 *          PackMigrators()
 *          UnpackMigrators()
 *
 **************************************************************************/
#include <stdio.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "QueueOps.h"
#include "Decomp.h"


/*
 *      Define some tags for MPI messages related to node migration
 */
#define MSG_MIG_LENGTH 3000
#define MSG_MIG_NODES  3001

/*
 *      Define some constants indicating the number of values that
 *      must be communicated for each node and arm in order to properly
 *      migrate a node.
 */
#ifdef _ARLFEM
#define VALS_PER_MIG_NODE  16
#else
#define VALS_PER_MIG_NODE  13
#endif  /* ifdef _ARLFEM */

// VALS_PER_MIG_ARM  =  node dom id, node id, 3 b + 3 n
#define VALS_PER_MIG_ARM       11


#define VALS_PER_MIG_EXTRA     2
#define VALS_PER_MIG_INCLUSION 20 // SA: 15->20 add new inclusion format.

/*
 *      Define a structure that will hold some info about entities
 *      being migrated to a given remote domain
 */
typedef struct {
        int migNodeCount;
        int *migNodeList;
        int migInclusionCount;
        int *migInclusionList;
} MigInfo_t;


#ifdef PARALLEL

#ifdef ESHELBY
/*---------------------------------------------------------------------------
 *
 *      Author:       Gregg Hommes
 *
 *      Function:     CleanInclusionArray
 *
 *      Description:  Removes all non-local inclusions from the inclusion
 *                    array.  Local inclusions that have been migrated have
 *                    been flagged, and ghost inclusions are known because
 *                    they follow all local inclusions in the array.
 *                    
 *      Arguments:
 *
 *      Last Modified:  02/22/2008 - original version
 *
 *-------------------------------------------------------------------------*/
static void CleanInclusionArray(Home_t *home)
{
/*
 *      Get rid of current ghost inclusions by simply dropping
 *      the total inclusion count to the count of local inclusions.
 *      The inclusion array will only be freed if the inclusion
 *      count has dropped to zero.
 */
        home->totInclusionCount = home->locInclusionCount;

        int numInclusions = home->locInclusionCount;
        int incIndex      = 0;

        while (1) {

/*
 *          Break out of loop when we've checked all the 
 *          previously local inclusions
 */
            if (incIndex >= numInclusions) {
                break;
            }

            EInclusion_t *inclusion = &home->eshelbyInclusions[incIndex];

/*
 *          A non-negative <id> indicates the inclusions is still local
 *          so skip to the next one.
 */
            if (inclusion->id >= 0) {
                incIndex++;
                continue;
            }

/*
 *          Gotta delete an inclusion, so copy the final inclusion
 *          in the list on top of this one and decrement the inclusion
 *          count.
 *
 *          NOTE: We don't worry about the keeping the cell inclusionQ's
 *          up to date here since they will be completely rebuilt after
 *          the migration and updated during the ghost communications.
 */
            numInclusions--;

            if (incIndex != numInclusions) {
                memcpy(inclusion, &home->eshelbyInclusions[numInclusions],
                       sizeof(EInclusion_t));
            }

        }

        home->locInclusionCount = numInclusions;
        home->totInclusionCount = numInclusions;

/*
 *      If there are no more inclusions in this domain, free
 *      the array of inclusions
 */
        if (numInclusions == 0) {
            if (home->eshelbyInclusions != (EInclusion_t *)NULL) {
                free(home->eshelbyInclusions);
                home->eshelbyInclusions = (EInclusion_t *)NULL;
            }
        }

        return;
}
#endif  // ESHELBY


/*---------------------------------------------------------------------------
 *
 *      Author:       Gregg Hommes
 *
 *      Function:     BuildMigLists
 *
 *      Description:  Identify all local entities that are outside the
 *                    local domains current boundaries (due to motion of
 *                    the entity or the domain boundaries themselves) and
 *                    require migration to remote domains.  Arrays (1 per
 *                    domain to which entities will be migrated) will be
 *                    created with indices of the entities to be migrated
 *                    to the remote domains along with counts of the entities
 *                    being migrated.
 *
 *      Arguments:
 *          OUT: migCommList Pointer to an array of integers (1 per domain) 
 *                    On return to the caller, each element of the array
 *                    corresponding to a domain to which entities are to
 *                    be migrated will be set to 1
 *
 *          OUT: migInfo  Pointer in which to return to caller an array of
 *                    structures containing information on types and
 *                    counts of entities ebing migrated.  Returned array
 *                    will contain <numSendBufs> elements corresponding
 *                    to the domain ID returned in <sendBufDest>.
 *
 *          OUT: sendBufDest Pointer in which to return to the caller an
 *                    array of domain IDs to which this domain will be 
 *                    migrating nodes.
 *
 *          OUT: numSendBufs Pointer to location in which to return to the 
 *                    caller the number of domains to which the local
 *                    domain will be migrating nodes.
 *
 *-------------------------------------------------------------------------*/
static void BuildMigLists(Home_t *home, int *migCommList, MigInfo_t **migInfo,
                          int **sendBufDest, int *numSendBufs)
{
        int          i;
        int          nodeCount, thisDomain, destDom;
        int          allocatedValues, sendBufCnt;
        int          sendIndex;
        int          *sendBufList;
        Node_t       *node;
        MigInfo_t    *migInfoList;

        nodeCount = home->newNodeKeyPtr;
        thisDomain = home->myDomain;

        sendBufCnt = 0;
        allocatedValues = 0;
        sendBufList = (int *)NULL;
        migInfoList = (MigInfo_t *)NULL;

        for (i = 0; i < nodeCount; i++) {
            int migNodeCount;

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

/*
 *          Search the domain decomposition for the domain which
 *          should own this node.  If the current domain should retain
 *          ownership, move on to the next node.
 */
            destDom = FindCoordDomain(home, 1, &node->x, &node->y, &node->z);

            if (destDom == thisDomain) {
                continue;
            }

/*
 *          This node must be migrated.  If the remote domain has not
 *          already been added to the list of domains to which items must
 *          be migrated, add it to the list and add a structure to track
 *          entities being migrated to that particular remote domain.
 */
            for (sendIndex = 0; sendIndex < sendBufCnt; sendIndex++) {

                if (destDom == sendBufList[sendIndex]) {
                    break;
                }
            }

            if (sendIndex == sendBufCnt) {

                sendBufCnt += 1;

                if (sendBufCnt > allocatedValues) {

                    allocatedValues += 25;

                    sendBufList = (int *)realloc(sendBufList, allocatedValues *
                                                 sizeof(int));

                    migInfoList = (MigInfo_t *)realloc(migInfoList,
                                                       allocatedValues *
                                                       sizeof(MigInfo_t));
                }

                sendBufList[sendBufCnt-1] = destDom;
                memset(&migInfoList[sendBufCnt-1], 0, sizeof(MigInfo_t));
            }

/*
 *          Add the node to the list of nodes being sent to that particular
 *          remote domain and increment the count of nodes beng sent there.
 *          Also set the flag indicating entities are being migrated to
 *          this remote domain 
 */
            migNodeCount = migInfoList[sendIndex].migNodeCount;

            migInfoList[sendIndex].migNodeList =
                    (int *)realloc(migInfoList[sendIndex].migNodeList,
                                   (migInfoList[sendIndex].migNodeCount+1) *
                                   sizeof(int));

            migInfoList[sendIndex].migNodeList[migNodeCount] = i;
            migInfoList[sendIndex].migNodeCount += 1;

            migCommList[destDom] = 1;
        }

#ifdef ESHELBY
/*
 *      If the simulation includes Eshelby inclusions, go through the
 *      local inclusion list and determine if any inclusions need to 
 *      be migrated.  
 */
        if (home->param->enableInclusions==1) 
        {
           for (i = 0; i < home->locInclusionCount; i++) {
              int migInclusionCount;
              
/*
 *              Search the domain decomposition for the domain which
 *              should own this inclusion.  If the current domain should
 *              retain ownership, move on to the next inclusion.
 */
              EInclusion_t *inclusion = &home->eshelbyInclusions[i];
              
              destDom = FindCoordDomain(home, 0, &inclusion->position[X],
                                        &inclusion->position[Y],
                                        &inclusion->position[Z]);
              
              if (destDom == thisDomain) {
                 continue;
              }
           
/*
 *              This inclusion must be migrated.  If the remote domain has
 *              not already been added to the list of domains to which items
 *              must be migrated, add it to the list and add a structure to
 *              track entities being migrated to that particular remote domain.
 */
              for (sendIndex = 0; sendIndex < sendBufCnt; sendIndex++) {
                 
                 if (destDom == sendBufList[sendIndex]) {
                    break;
                 }
              }
              
              if (sendIndex == sendBufCnt) {
                 
                 sendBufCnt += 1;
                 
                 if (sendBufCnt > allocatedValues) {
                    
                    allocatedValues += 25;
                    
                    sendBufList = (int *)realloc(sendBufList,
                                                 allocatedValues *
                                                 sizeof(int));
                    
                    migInfoList = (MigInfo_t *)realloc(migInfoList,
                                                       allocatedValues *
                                                       sizeof(MigInfo_t));
                 }
                 
                 sendBufList[sendBufCnt-1] = destDom;
                 memset(&migInfoList[sendBufCnt-1], 0, sizeof(MigInfo_t));
              }
              
/*
 *              Add the inclusion to the list of inclusions being sent to
 *              that particular remote domain and increment the count of
 *              inclusions being sent there.  Also set the flag indicating
 *              entities are being migrated to this remote domain 
 */
              migInclusionCount = migInfoList[sendIndex].migInclusionCount;
              
              migInfoList[sendIndex].migInclusionList =
                 (int *)realloc(migInfoList[sendIndex].migInclusionList,
                                (migInclusionCount+1) * sizeof(int));
              
              migInfoList[sendIndex].migInclusionList[migInclusionCount] = i;
              migInfoList[sendIndex].migInclusionCount += 1;
              
              migCommList[destDom] = 1;

           } // enableInclusions ends
        }
#endif // ESHELBY

        *migInfo     = migInfoList;
        *sendBufDest = sendBufList;
        *numSendBufs = sendBufCnt;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     PackMigrators
 *      Description:  Allocate and pack a buffer with all data necessary
 *                    to migrate the entities specified in <migInfo> to 
 *                    a remote domain.
 *
 *      Arguments:
 *          IN: migInfo    Pointer to structure containing information on
 *                         types and counts of entities being migrated.
 *          IN: sendBuf    Pointer in which to return the pointer to the
 *                         allocated buffer.  Caller will be responsible for
 *                         freeing the buffer.
 *          IN: sendBufLen Length (in bytes) of the buffer returned to the
 *                         caller.
 *
 *-------------------------------------------------------------------------*/
static void PackMigrators(Home_t *home, MigInfo_t *migInfo, real8 **sendBuf, int *sendBufLen)
{
/*
 *      Determine how large a buffer is required to hold all the
 *      specified entities and allocate an appropriately sized buffer.
 */
        int armCount     = 0;
        int bufIndex     = 0;
        int migNodeCount = migInfo->migNodeCount;

        for (int i=0; (i<migNodeCount); i++)
            armCount += home->nodeKeys[migInfo->migNodeList[i]]->numNbrs;

        int numVals =                  VALS_PER_MIG_EXTRA 
                      + migNodeCount * VALS_PER_MIG_NODE
                      + armCount     * VALS_PER_MIG_ARM  ;

/*
 *      If eshelby inclusions are being simulated, adjust the buffer
 *      size to accomodate the migrating inclusions.
 */
#ifdef ESHELBY
        int migInclusionCount = 0;
        if (home->param->enableInclusions == 1) 
        {
           migInclusionCount = migInfo->migInclusionCount;
           numVals += migInclusionCount * VALS_PER_MIG_INCLUSION;
        }
#endif

        real8 *buf = (real8 *)malloc(numVals * sizeof(real8));

        buf[bufIndex++] = (real8)migNodeCount;
#ifdef ESHELBY
        buf[bufIndex++] = (real8)migInclusionCount;
#else
        buf[bufIndex++] = 0;
#endif

/*
 *      Loop through all the nodes to be sent and pack the necessary
 *      nodal data into the buffer.
 */
        for (int i=0; (i<migNodeCount); i++) 
        {
            Node_t *node = home->nodeKeys[migInfo->migNodeList[i]];

            buf[bufIndex++] = (real8)node->myTag.index;
            buf[bufIndex++] = (real8)node->constraint;
            buf[bufIndex++] = (real8)node->numNbrs;
            buf[bufIndex++] = (real8)node->flags;

            buf[bufIndex++] = node->x;
            buf[bufIndex++] = node->y;
            buf[bufIndex++] = node->z;

            buf[bufIndex++] = node->vX;
            buf[bufIndex++] = node->vY;
            buf[bufIndex++] = node->vZ;

            buf[bufIndex++] = node->oldvX;
            buf[bufIndex++] = node->oldvY;
            buf[bufIndex++] = node->oldvZ;

#ifdef _ARLFEM
            buf[bufIndex++] = node->surfaceNorm[0];
            buf[bufIndex++] = node->surfaceNorm[1];
            buf[bufIndex++] = node->surfaceNorm[2];
#endif  /* ifdef _ARLFEM */

            for (int j=0; (j<node->numNbrs); j++) 
            {
                buf[bufIndex++] = (real8)node->nbrTag[j].domainID;
                buf[bufIndex++] = (real8)node->nbrTag[j].index;

                buf[bufIndex++] = node->burgX[j];
                buf[bufIndex++] = node->burgY[j];
                buf[bufIndex++] = node->burgZ[j];

                buf[bufIndex++] = node->nx[j];
                buf[bufIndex++] = node->ny[j];
                buf[bufIndex++] = node->nz[j];

                buf[bufIndex++] = node->armfx[j];
                buf[bufIndex++] = node->armfy[j];
                buf[bufIndex++] = node->armfz[j];
            }

/*
 *          Free the node and recyle the tag
 */
            FreeNode(home, migInfo->migNodeList[i]);

        }  /* for (i = 0; i < migCount; ...) */

#ifdef ESHELBY
/*
 *      If there are any inclusions to be migrated, pack them into
 *      the buffers for the remote domains now.
 */
        if (home->param->enableInclusions == 1) 
        {
           for (int i=0; (i<migInclusionCount); i++) 
           {
              int           inclusionIndex = migInfo->migInclusionList[i]; 
              EInclusion_t *inclusion      = &home->eshelbyInclusions[inclusionIndex];
              
              buf[bufIndex++] = (real8)inclusion->id;
              buf[bufIndex++] = (real8)inclusion->cellID;
              
              buf[bufIndex++] = inclusion->radius[0];
              buf[bufIndex++] = inclusion->radius[1];
              buf[bufIndex++] = inclusion->radius[2];
              
              buf[bufIndex++] = inclusion->position[0];
              buf[bufIndex++] = inclusion->position[1];
              buf[bufIndex++] = inclusion->position[2];
              
              buf[bufIndex++] = inclusion->rotation[0][0];
              buf[bufIndex++] = inclusion->rotation[0][1];
              buf[bufIndex++] = inclusion->rotation[0][2];
              
              buf[bufIndex++] = inclusion->rotation[1][0];
              buf[bufIndex++] = inclusion->rotation[1][1];
              buf[bufIndex++] = inclusion->rotation[1][2];
              
              buf[bufIndex++] = inclusion->strainField[0];
              buf[bufIndex++] = inclusion->strainField[1];
              buf[bufIndex++] = inclusion->strainField[2];
              buf[bufIndex++] = inclusion->strainField[3];
              buf[bufIndex++] = inclusion->strainField[4];
              buf[bufIndex++] = inclusion->strainField[5];
              
              inclusion->id = -1;
           }
        } // end enableInclusions

/*
 *          NOTE: Local inclusions that have been migrated out will
 *          be removed from the inclusion list after buffers have been
 *          packed for ALL remote domains
 */
#endif // ESHELBY

/*
 *      Return pointer to buffer and the buffer size to the caller.
 *      Caller is responsible for freeing the buffer.
 */
        *sendBuf    = buf;
        *sendBufLen = numVals * sizeof(real8);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     UnpackMigrators
 *      Description:  Process nodal data migrated from the specified
 *                    remote domain.
 *
 *      Arguments:
 *          buf      Array containing the nodal data
 *          remDomID ID of the remote domain from which this buffer
 *                   was received.
 *
 *-------------------------------------------------------------------------*/
static void UnpackMigrators(Home_t *home, real8 *buf, int remDomID)
{
        int thisDomain = home->myDomain;
        int bufIndex   = 0;

/*
 *      Get the node count from the first element of the buffer.
 *      The second element contains the number of inclusions being migrated
 */
        int nodeCount      = (int)buf[bufIndex++];
        int inclusionCount = (int)buf[bufIndex++];

        for (int i=0; (i<nodeCount); i++) 
        {
/*
 *          Add a new node to the list of local nodes and populate the
 *          node structure with data from the remote domain.  Also
 *          add a mapping between the original node tag from the
 *          remote domain and the new node tag in the local domain.
 */
            int     nodeIndex = GetFreeNodeTag(home);
            Node_t *node      = PopFreeNodeQ  (home);
            int     numNbrs   = 0;

            home->nodeKeys[nodeIndex] = node;

            node->fX = 0.0;
            node->fY = 0.0;
            node->fZ = 0.0;

            node->myTag.domainID = thisDomain;
            node->myTag.index    = nodeIndex;

            Tag_t oldTag;

            oldTag.domainID = remDomID;
            oldTag.index    = (int)buf[bufIndex++];

            AddTagMapping(home, &oldTag, &node->myTag);

            node->constraint  = (int)buf[bufIndex++];
            numNbrs           = (int)buf[bufIndex++];
            node->flags       = (int)buf[bufIndex++];

            node->x     = buf[bufIndex++];
            node->y     = buf[bufIndex++];
            node->z     = buf[bufIndex++];

            node->vX    = buf[bufIndex++];
            node->vY    = buf[bufIndex++];
            node->vZ    = buf[bufIndex++];

            node->oldvX = buf[bufIndex++];
            node->oldvY = buf[bufIndex++];
            node->oldvZ = buf[bufIndex++];

#ifdef _ARLFEM
            node->surfaceNorm[0] = buf[bufIndex++];
            node->surfaceNorm[1] = buf[bufIndex++];
            node->surfaceNorm[2] = buf[bufIndex++];
#endif

/*
 *          Set all the segment specific nodal data values.
 */
            AllocNodeArms(node,numNbrs);

            for (int j=0; (j<numNbrs); j++) 
            {
                node->nbrTag[j].domainID = (int)buf[bufIndex++];
                node->nbrTag[j].index    = (int)buf[bufIndex++];

                node->burgX[j] = buf[bufIndex++];
                node->burgY[j] = buf[bufIndex++];
                node->burgZ[j] = buf[bufIndex++];

                node->nx   [j] = buf[bufIndex++];
                node->ny   [j] = buf[bufIndex++];
                node->nz   [j] = buf[bufIndex++];

                node->armfx[j] = buf[bufIndex++];
                node->armfy[j] = buf[bufIndex++];
                node->armfz[j] = buf[bufIndex++];

                node->fX += node->armfx[j];
                node->fY += node->armfy[j];
                node->fZ += node->armfz[j];

            }
        }

        if (inclusionCount > 0) 
        {
#ifdef ESHELBY
/*
 *          Increase the size of the local array of Eshelby inclusions
 *          with enough space to hold all the incoming inclusions.
 */
            int newSize = (home->locInclusionCount + inclusionCount) * sizeof(EInclusion_t);

            home->eshelbyInclusions = (EInclusion_t *)realloc(home->eshelbyInclusions, newSize);

            for (int i=0; (i<inclusionCount); i++) 
            {
                int           newIndex  =  home->locInclusionCount;
                EInclusion_t *inclusion = &home->eshelbyInclusions[newIndex];

                inclusion->id =     (int)buf[bufIndex++];
                inclusion->cellID = (int)buf[bufIndex++];
                inclusion->nextInCell = -1;

                inclusion->radius[0] = buf[bufIndex++];
                inclusion->radius[1] = buf[bufIndex++];
                inclusion->radius[2] = buf[bufIndex++];

                inclusion->position[0] = buf[bufIndex++];
                inclusion->position[1] = buf[bufIndex++];
                inclusion->position[2] = buf[bufIndex++];

                inclusion->rotation[0][0] = buf[bufIndex++];
                inclusion->rotation[0][1] = buf[bufIndex++];
                inclusion->rotation[0][2] = buf[bufIndex++];
                inclusion->rotation[1][0] = buf[bufIndex++];
                inclusion->rotation[1][1] = buf[bufIndex++];
                inclusion->rotation[1][2] = buf[bufIndex++];

                inclusion->strainField[0] = buf[bufIndex++];
                inclusion->strainField[1] = buf[bufIndex++];
                inclusion->strainField[2] = buf[bufIndex++];
                inclusion->strainField[3] = buf[bufIndex++];
                inclusion->strainField[4] = buf[bufIndex++];
                inclusion->strainField[5] = buf[bufIndex++];

                home->locInclusionCount += 1;
                home->totInclusionCount += 1;
            }
#endif // ESHELBY
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CommSendMigrators
 *      Description:  Handles communications required to migrate
 *                    nodes to/from remote domains.
 *
 *      Arguments:
 *          IN: migCommList Array of integers (1 per domain) where each value
 *                          is set to 0 or 1 indicating if the current domain
 *                          will be migrating nodes to the corresponding
 *                          remote domain.
 *          IN: migInfo     Pointer to structure containing information on
 *                          types and counts of entities being migrated.
 *          IN: sendBufDest Array of of domain IDs to which this domain will
 *                          be migrating nodes.
 *          IN: numSendBufs The number of domains to which the local domain
 *                          will be migrating nodes.
 *
 *-------------------------------------------------------------------------*/
static void CommSendMigrators(Home_t *home, int *migCommList,
                              MigInfo_t *migInfo, int *sendBufDest,
                              int numSendBufs)
{
        int          *sendBufLen =0,  *recvBufLen =0;
        char        **sendBuf    =0, **recvBuf    =0;
        MPI_Request  *sendRequest=0,  *recvRequest=0;
        MPI_Status   *sendStatus =0,  *recvStatus =0;

        int numDomains = home->numDomains;
        int thisDomain = home->myDomain;

/*
 *      First we need to determine how many remote domains will be
 *      migrating nodes to this domain.  Each domain has already set
 *      the flag in the migCommList for all domains to which it will
 *      migrate nodes.  When the all-reduce is done, each domain
 *      will know how many other domains will be migrating nodes
 *      to it...
 */
        int *glblMigCommList = home->globalReduceBuf;

        MPI_Allreduce(migCommList, glblMigCommList, numDomains, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        int numRecvBufs = glblMigCommList[thisDomain];

/*
 *      Pack buffers for each remote domain to which this domain
 *      will be migrating nodes
 */
        if (numSendBufs > 0) 
        {
            sendBuf     = (char       **) calloc(1, numSendBufs * sizeof(real8 *    ));
            sendBufLen  = (int         *) calloc(1, numSendBufs * sizeof(int        ));
            sendRequest = (MPI_Request *) malloc(   numSendBufs * sizeof(MPI_Request));
            sendStatus  = (MPI_Status  *) malloc(   numSendBufs * sizeof(MPI_Status ));
        }

        for (int i=0; (i<numSendBufs); i++)
            PackMigrators(home, &migInfo[i], (real8 **)&sendBuf[i], &sendBufLen[i]);

#ifdef ESHELBY
/*
 *      If we are simulating Eshelby inclusions, we need to remove
 *      all non-local (migrated or ghost) inclusions from the array.
 */
        CleanInclusionArray(home);
#endif

/*
 *      Allocate arrays for handling incoming migrated nodes
 */
        if (numRecvBufs>0) 
        {
            recvBuf     = (char       **) calloc(1, numRecvBufs * sizeof(real8 *    ));
            recvBufLen  = (int         *) calloc(1, numRecvBufs * sizeof(int        ));
            recvRequest = (MPI_Request *) malloc(   numRecvBufs * sizeof(MPI_Request));
            recvStatus  = (MPI_Status  *) malloc(   numRecvBufs * sizeof(MPI_Status ));
        }

/*
 *      Pre-issue receives of buffer lengths from all domains that
 *      will be be migrating nodes to this domain.  Lengths are
 *      specified in units of bytes.
 */
        for (int i=0; (i<numRecvBufs); i++) 
            MPI_Irecv(&recvBufLen[i], 1, MPI_INT, MPI_ANY_SOURCE, MSG_MIG_LENGTH, MPI_COMM_WORLD, &recvRequest[i]);

/*
 *      Have the current domain send out the sizes of buffers it will
 *      be transmitting.
 */
        for (int i=0; (i<numSendBufs); i++) 
            MPI_Isend(&sendBufLen[i], 1, MPI_INT, sendBufDest[i], MSG_MIG_LENGTH, MPI_COMM_WORLD, &sendRequest[i]);

/*
 *      Wait for all length send/receives to complete
 */
        MPI_Waitall(numSendBufs, sendRequest, sendStatus);
        MPI_Waitall(numRecvBufs, recvRequest, recvStatus);

/*
 *      Allocate receive buffers of the appropriate sizes and post the
 *      receives associated with those buffers.
 */
        for (int i=0; (i<numRecvBufs); i++) 
        {
            recvBuf[i] = (char *)malloc(recvBufLen[i]);
            int numValues  = recvBufLen[i] / sizeof(real8);
            MPI_Irecv( (double *) recvBuf[i], numValues, MPI_DOUBLE, recvStatus[i].MPI_SOURCE, MSG_MIG_NODES, MPI_COMM_WORLD, &recvRequest[i]);
        }

/*
 *      Send all the migrating nodes to the appropriate remote domains.
 */
        for (int i=0; (i<numSendBufs); i++) 
        {
            int numValues = sendBufLen[i] / sizeof(real8);
            MPI_Isend( (double *) sendBuf[i], numValues, MPI_DOUBLE, sendBufDest[i], MSG_MIG_NODES, MPI_COMM_WORLD, &sendRequest[i]);
        }

/*
 *      Process the incoming buffers as soon as they arrive.  Status
 *      will be placed in first recvStatus struct each time.
 */
#define EXP_NO_WAITANY
#ifdef  EXP_NO_WAITANY
        // This is David Boehme's suggested replacement for the MPI_Waitany()
        // loop that has apparently better scalable performance.

        MPI_Waitall(numRecvBufs, recvRequest, recvStatus);

        for (int i=0; (i<numRecvBufs); i++) 
        {
            UnpackMigrators(home, (real8 *)recvBuf[i], recvStatus[i].MPI_SOURCE);
            if (recvBuf[i]) { free(recvBuf[i]); recvBuf[i]=0; }
        }
#else
        // (original solution)

        for (int i=0; (i<numRecvBufs); i++) 
        {
            MPI_Waitany    (numRecvBufs, recvRequest, &recvIndex, &recvStatus[0]);
            UnpackMigrators(home, (real8 *)recvBuf[recvIndex], recvStatus[0].MPI_SOURCE);

            if (recvBuf[recvIndex]) { free(recvBuf[recvIndex]); recvBuf[recvIndex]=0; }
        }
#endif

/*
 *      Wait for all buffer sends to complete.
 */
        MPI_Waitall(numSendBufs, sendRequest, sendStatus);

/*
 *      Free all the temporary buffers before returning to the caller.
 */
        if (numSendBufs > 0) 
        {
            for (int i=0; (i<numSendBufs); i++) 
            { 
               if (sendBuf[i]) { free(sendBuf[i]); sendBuf[i]=0; } 
            }
            if (sendBuf    ) { free(sendBuf    ); sendBuf    =0; }
            if (sendBufLen ) { free(sendBufLen ); sendBufLen =0; }
            if (sendRequest) { free(sendRequest); sendRequest=0; }
            if (sendStatus ) { free(sendStatus ); sendStatus =0; }
        }

        if (numRecvBufs > 0) 
        {
            if (recvBuf    ) { free(recvBuf    ); recvBuf    =0; }
            if (recvBufLen ) { free(recvBufLen ); recvBufLen =0; }
            if (recvRequest) { free(recvRequest); recvRequest=0; }
            if (recvStatus ) { free(recvStatus ); recvStatus =0; }
        }

        return;
}
#endif  /* ifdef PARALLEL */


/*---------------------------------------------------------------------------
 *
 *      Function:     Migrate
 *      Description:  Driver function to control the process of
 *                    identifying nodes whose ownership must be
 *                    transferred from one domain to another and
 *                    coordinating the transfer of ownership.
 *
 *-------------------------------------------------------------------------*/
void Migrate(Home_t *home)
{
#ifdef PARALLEL
        TimerStart(home, MIGRATION);
/*
 *      Nodes will only migrate if we're running in parallel
 */
        int        numDomains  = home->numDomains;
        int        numSendBufs = 0;
        int       *sendBufDest = 0;
        MigInfo_t *migInfo     = 0;

        int       *migCommList = home->localReduceBuf;

/*
 *      Since we'll be incrementing values in migCommList, we need to
 *      re-initialize it to be sure we have valid starting values.
 */
        for (int i=0; (i<numDomains); i++) { migCommList[i]=0; }

/*
 *      Look through all local nodes and determine which nodes need
 *      to be migrated and the domains to which those nodes must
 *      be migrated.  For each remote domain to which this domain
 *      will migrate one or more nodes, build a list of the nodes
 *      to be sent to that domain.
 */
        BuildMigLists(home, migCommList, &migInfo, &sendBufDest, &numSendBufs);

/*
 *      Send out all nodes (if any) that need to be migrated
 *      to remote domains and receive any nodes migrating from
 *      other domains.
 */
        CommSendMigrators(home, migCommList, migInfo, sendBufDest, numSendBufs);

/*
 *      All migrated nodes have been retagged.  Each domain now needs
 *      to communicate to its neighboring domains the mappings between
 *      the old and new tags for all nodes it received during the
 *      migration.  Once that is done, the local domains go through
 *      all their own nodes reconciling the node tag changes.
 */
        DistributeTagMaps(home);

/*
 *      Free up all temporary arrays before returning to the caller.
 */
        for (int i=0; (i<numSendBufs); i++) 
        {
            if (migInfo[i].migNodeList     ) { free(migInfo[i].migNodeList     ); migInfo[i].migNodeList     =0; }
            if (migInfo[i].migInclusionList) { free(migInfo[i].migInclusionList); migInfo[i].migInclusionList=0; }
        }

        if (migInfo    ) { free(migInfo    ); migInfo    =0; }
        if (sendBufDest) { free(sendBufDest); sendBufDest=0; }

        TimerStop(home, MIGRATION);

#ifdef SYNC_TIMERS
/*
 *      Measure dead time after node migration
 */
        TimerStart(home, MIGRATION_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, MIGRATION_BARRIER);
#endif  // SYNC_TIMERS

#endif  // PARALLEL
}
