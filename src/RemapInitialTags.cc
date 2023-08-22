/***************************************************************************
 *
 *      Module:       RemapInitialTags.c
 *      Description:  This module contains the functions needed to 
 *                    reconcile between the domains node that have
 *                    been retagged.
 *
 *                    These functions are only necessary at initialization
 *                    time when nodes are assigned and distributed to
 *                    (possibly) new domains and get retagged with new
 *                    ids.
 *
 *      Included functions:
 *
 *          DistributeTagMaps()
 *          PackTagMap()
 *          RemapArmTags()
 *          TagMapCompareNew()
 *          TagMapCompareOld()
 *          UnpackTagMap()
 *
 **************************************************************************/
#include <stdio.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"


/*-------------------------------------------------------------------------
 *
 *      Function:       TagMapCompareOld
 *      Description:    Compares two TagMap_t structures based on the
 *                      <oldTag> values.  This function is compatible
 *                      for use by system sorting and searching functions
 *                      such as qsort(), bsearch(), etc.
 *
 *      Returns:        -1 if  a <  b
 *                       0 if  a == b
 *                       1 if  a >  b
 *
 *------------------------------------------------------------------------*/
int TagMapCompareOld(const void *a, const void *b)
{
        TagMap_t *tagMap1 = (TagMap_t *)a;
        TagMap_t *tagMap2 = (TagMap_t *)b;

        if (tagMap1->oldTag.domainID < tagMap2->oldTag.domainID) return(-1);
        if (tagMap1->oldTag.domainID > tagMap2->oldTag.domainID) return(1);

        if (tagMap1->oldTag.index < tagMap2->oldTag.index) return(-1);
        if (tagMap1->oldTag.index > tagMap2->oldTag.index) return(1);

        return(0);
}


/*-------------------------------------------------------------------------
 *
 *      Function:       TagMapCompareNew
 *      Description:    Compares two TagMap_t structures based on the
 *                      <newTag> values.  This function is compatible
 *                      for use by system sorting and searching functions
 *                      such as qsort(), bsearch(), etc.
 *
 *      Returns:        -1 if  a <  b
 *                       0 if  a == b
 *                       1 if  a >  b
 *
 *------------------------------------------------------------------------*/
int TagMapCompareNew(const void *a, const void *b)
{
        TagMap_t *tagMap1 = (TagMap_t *)a;
        TagMap_t *tagMap2 = (TagMap_t *)b;

        if (tagMap1->newTag.domainID < tagMap2->newTag.domainID) return(-1);
        if (tagMap1->newTag.domainID > tagMap2->newTag.domainID) return(1);

        if (tagMap1->newTag.index < tagMap2->newTag.index) return(-1);
        if (tagMap1->newTag.index > tagMap2->newTag.index) return(1);

        return(0);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     PackTagMap
 *      Description:  Packs a buffer with old/new tag mappings for
 *                    all local nodes which were retagged during
 *                    the initial node distribution.
 *      Arguments:
 *          buf       Pointer to location in which to return to the
 *                    caller a pointer to the buffer created by
 *                    this function.
 *          bufSize   Pointer to location in which to return to the
 *                    caller the size (in bytes) of the buffer
 *                    returned in <buf>.
 *
 *------------------------------------------------------------------------*/
void PackTagMap(Home_t *home, int **buf, int *bufSize)
{
        int      i, mapBufEnts, mapBufSize, bufIndex;
        int      *mapBuf;
        TagMap_t *mapping;

        mapBufEnts = 1 + (home->tagMapEnts * 4);
        mapBufSize = mapBufEnts * sizeof(int);
        mapBuf = (int *)malloc(mapBufSize);

        bufIndex = 0;

        mapBuf[bufIndex++] = home->tagMapEnts;

        for (i = 0; i < home->tagMapEnts; i++) {
 
            mapping = &home->tagMap[i];

            mapBuf[bufIndex++] = mapping->oldTag.domainID;
            mapBuf[bufIndex++] = mapping->oldTag.index;
            mapBuf[bufIndex++] = mapping->newTag.domainID;
            mapBuf[bufIndex++] = mapping->newTag.index;
        }

        *buf     = mapBuf;
        *bufSize = mapBufSize;

        return;

}


/*-------------------------------------------------------------------------
 *
 *      Function:     UnpackTagMap
 *      Description:  Unpacks the provided tag map from a remote
 *                    domain and adds all old/new tag mappings
 *                    to the local tag map.
 *
 *      Arguments:
 *          buf       Pointer to buffer to be unpacked.  Buffer is 
 *                    assumed to be an array of integers.
 *
 *------------------------------------------------------------------------*/
void UnpackTagMap(Home_t *home, int *buf)
{
        int     i, numMappings, bufIndex;
        Tag_t   oldTag, newTag;

        bufIndex = 0;

        numMappings = buf[bufIndex++];

        for (i = 0; i < numMappings; i++) {

            oldTag.domainID = buf[bufIndex++];
            oldTag.index    = buf[bufIndex++];
            newTag.domainID = buf[bufIndex++];
            newTag.index    = buf[bufIndex++];

            AddTagMapping(home, &oldTag, &newTag);
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     RemapArmTags
 *      Description:  Loop through all arms of the local nodes retagging
 *                    the arms if the node at the far end of the arm
 *                    was retagged during initialization.  The full
 *                    tag map (home->tagMap) must have been populated
 *                    with new mappings for the local and all nearby
 *                    domains, and sorted before this function is called.
 *
 *                    Note:  We don't alter the local node tags here
 *                    since they were retagged (if necessary) during
 *                    the initial node distribution.
 *
 *------------------------------------------------------------------------*/
static void RemapArmTags(Home_t *home)
{
        int      i, armID, maxNodeKey;
        Node_t   *node;
        TagMap_t key;
        TagMap_t *mapping;

        maxNodeKey = home->newNodeKeyPtr;

        for (i = 0; i < maxNodeKey; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            for (armID = 0; armID < node->numNbrs; armID++) {

                key.oldTag.domainID = node->nbrTag[armID].domainID;
                key.oldTag.index    = node->nbrTag[armID].index;

/*
 *              Use the arm's tag to look for any new tag mapping.
 *              If no mapping was found, the tag is fine as it is.
 *              Otherwise, reset the arm's tag to the new value.
 */
                mapping = (TagMap_t *)bsearch(&key, home->tagMap,
                                              home->tagMapEnts,
                                              sizeof(TagMap_t),
                                              TagMapCompareOld);

                if (mapping != (TagMap_t *)NULL) {
                    node->nbrTag[armID].domainID = mapping->newTag.domainID;
                    node->nbrTag[armID].index    = mapping->newTag.index;
                }
            }
        }

        return;
}



/*-------------------------------------------------------------------------
 *
 *      Function:     DistributeTagMaps
 *      Description:  Have each domain send to all its near-neighbors
 *                    a mapping between old and new tags for all nodes
 *                    the domain retagged during initialization.  And
 *                    of course, receive tag mappings from the remote
 *                    domains.
 *
 *------------------------------------------------------------------------*/
void DistributeTagMaps(Home_t *home) 
{
#ifdef PARALLEL
        int             i, remDomIndex, reqIndex = 0, tagBufLen;
        int             numRemoteDomains;
        int             *tagBuf, *remBuf;
        int             **sendBuf, *sendBufLen;
        RemoteDomain_t  *remDom;
        MPI_Status      reqStatus;


        numRemoteDomains = home->remoteDomainCount;

/*
 *	Pre-issue receives of tag map buffer lengths from each neighbor
 */
        for (i = 0; i < numRemoteDomains; i++) {
            remDomIndex = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[remDomIndex];
        
            MPI_Irecv(&remDom->inBufLen, 1, MPI_INT, remDomIndex,
                      MSG_TAGREMAP_LEN, MPI_COMM_WORLD,
                      &home->inRequests[i]);
        }

/*
 *	Pack a buffer with mappings between old and new tags for any
 *      local nodes this domain has retagged.
 */
        PackTagMap(home, &tagBuf, &tagBufLen);

/*
 *      We're sending the same buffer to all remote domains, but 
 *      MPI spec says we can't re-use the same buffer for multiple
 *      asynch sends at the same time, so need to make multiple
 *      copies of the buffer.  (Use the existing buffer as the 
 *      one to be sent to the first remote domain)
 */
        if (numRemoteDomains > 0) {

            sendBuf = (int **)malloc(numRemoteDomains * sizeof(int *));
            sendBufLen = (int *)malloc(numRemoteDomains * sizeof(int));

            sendBuf[0] = tagBuf;
            sendBufLen[0] = tagBufLen;

/*
 *          For first remote domain the existing buffer (tagBuf) is used,
 *          so just need to set up remaining send buffers here.
 */
            for (i = 1; i < numRemoteDomains; i++) {
                sendBuf[i] = (int *)malloc(tagBufLen);
                sendBufLen[i] = tagBufLen;
                memcpy(sendBuf[i], tagBuf, tagBufLen);
            }
        }

/*
 *	Send the length length (in bytes) of the tag map buffer to
 *      the neighboring domains who will be receiving the buffer.
 */
        for (i = 0; i < numRemoteDomains; i++) {
            remDomIndex = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[remDomIndex];

            MPI_Isend(&sendBufLen[i], 1, MPI_INT, remDomIndex,
                      MSG_TAGREMAP_LEN, MPI_COMM_WORLD,
                      &home->outRequests[i]);
        }

/*
 *	Wait until we've received the tag map buffer lengths from all
 *	neighboring domains, allocate buffers for the maps and pre-
 *	issue receives for the data.
 */
        MPI_Waitall(numRemoteDomains, home->inRequests, home->inStatus);

        for (i = 0; i < numRemoteDomains; i++) {

            remDomIndex = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[remDomIndex];

            remDom->inBuf = (char *) malloc(remDom->inBufLen);

            MPI_Irecv( (int *) remDom->inBuf, remDom->inBufLen / sizeof(int), MPI_INT, remDomIndex, MSG_TAGREMAP, MPI_COMM_WORLD, &home->inRequests[i]);
        }

/*
 *	Wait for all length sends to complete
 */
        MPI_Waitall(numRemoteDomains, home->outRequests, home->outStatus);

/*
 *	Send the local tag-map to all neighboring domains
 */
        for (i = 0; i < numRemoteDomains; i++) {

            remDomIndex = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[remDomIndex];

            MPI_Isend(sendBuf[i], sendBufLen[i] / sizeof(int), MPI_INT,
                      remDomIndex, MSG_TAGREMAP, MPI_COMM_WORLD,
                      &home->outRequests[i]);
        } 

/*
 *	Wait for tag-map from any remote domain as long as there are any
 *	outstanding receives.
 */
        for (i = 0; i < numRemoteDomains; i++) {
            MPI_Waitany(numRemoteDomains, home->inRequests,
                        &reqIndex, &reqStatus);
            remDomIndex = home->remoteDomains[reqIndex];
            remBuf = (int *)home->remoteDomainKeys[remDomIndex]->inBuf;
            UnpackTagMap(home, remBuf);
            free(remBuf);
            home->remoteDomainKeys[remDomIndex]->inBuf = (char *)NULL;
        }

/*
 *	Now just wait for all sends to complete (should be long done by
 *	now) and free up any temporary buffers.
 */
        MPI_Waitall(numRemoteDomains, home->outRequests, home->outStatus);

        free(tagBuf);

        if (numRemoteDomains > 0) {

/*
 *          sendBuf[0] pointed to tagBuf which was freed above.  Now just
 *          free the rest.
 */
            for (i = 1; i < numRemoteDomains; i++) {
                free(sendBuf[i]);
            }

            free(sendBuf);
            free(sendBufLen);
        }

#endif

/*
 *      All new tag mappings are known, so sort the tag mappings
 *      and retag all arms that need it (if there are any in this
 *      domain).
 */
        if (home->tagMapEnts > 0) {
            qsort(home->tagMap, home->tagMapEnts, sizeof(TagMap_t),
                  TagMapCompareOld);

            RemapArmTags(home);
        }

/*
 *      No longer need the tagMap stuff in the home structure, so free
 *      the array and reinitialize associated values.
 */
        free(home->tagMap);

        home->tagMap = (TagMap_t *)NULL;
        home->tagMapSize = 0;
        home->tagMapEnts = 0;

        return;
}
