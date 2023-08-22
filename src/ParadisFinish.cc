/*-------------------------------------------------------------------------
 *
 *      Function:     ParadisFinish
 *      Description:  Handles pre-termination processing including
 *                    output generation, cleanup of any X-Window
 *                    display, release of any remaining dynamically
 *                    allocated memory, etc.
 *
 *-----------------------------------------------------------------------*/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Init.h"
#include "Util.h"
#include "DisplayC.h"
#include "Decomp.h"
#include "ParadisThread.h"
#include "Timestep.h"
#include "SharedMem.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif


/**************************************************************************
 *
 *      Function:     ReleaseMemory
 *
 *      Description:  This function is called during the program
 *                    termination procedure to release all remaining
 *                    dynamically allocated memory.  This is not truly
 *                    necessary, however, doing so facilitates the
 *                    process of hunting memory leaks by cleaning up
 *                    all the known allocated memory blocks.
 *
 *************************************************************************/
void ReleaseMemory(Home_t *home)
{
        int         i;
        Node_t      *node;
        Param_t     *param;
        NodeBlock_t *nodeBlock, *thisNodeBlock;

        if (home == (Home_t *)NULL) {
            return;
        }

        param = home->param;

/*
 *      Use existing functions (if available) to clean up some of the
 *      allocated storage.
 */
        DLBfreeOld(home);
        FreeCellTable(home);
        FreeRijm();
        FreeRijmPBC();
        FreeCellCenters();
        FreeCorrectionTable();
        FMFree(home);

#ifdef ESHELBY
/*
 *      Free the array of Eshelby inclusions if it exists
 */
        if (home->totInclusionCount > 0) { free(home->eshelbyInclusions); home->eshelbyInclusions=0; }
#endif


#ifdef PARALLEL
/*
 *      For multi-processor runs using the X-Window display capability
 *      we'll need to clean up the mirror domain structures
 */
        if (home->mirrorDomainKeys != (MirrorDomain_t **)NULL) 
        {
            for (int domID=0; domID < home->numDomains; domID++) 
            {
                if ( home->mirrorDomainKeys[domID] ) 
                {   MirrorDomain_Free(home, domID); }
            }
        }
#endif

/*
 *      Loop through all allocated blocks of node structures.  Free
 *      any arrays associated with individual nodes then free up the
 *      blocks of node structures.
 */
        nodeBlock = home->nodeBlockQ;

        while (nodeBlock) {

            node = nodeBlock->nodes;

            for (i = 0; i < NODE_BLOCK_COUNT; i++) {

                if (node->numNbrs) {

                    FreeNodeArms(node);

                    if (node->armCoordIndex != (int *)NULL) {
                        free(node->armCoordIndex);
                        node->armCoordIndex = NULL;
                    }
                }

                DESTROY_LOCK(&node->nodeLock);

                node++;
            }

            thisNodeBlock = nodeBlock;
            nodeBlock = nodeBlock->next;

            memset(thisNodeBlock->nodes, 0, NODE_BLOCK_COUNT * sizeof(Node_t));
            free(thisNodeBlock->nodes);
            free(thisNodeBlock);
        }

        home->nodeBlockQ = 0;

        if(home->nodeKeys) {
            free(home->nodeKeys);
            home->nodeKeys = NULL;
        }

/*
 *      Free the heap used for recycling nodes that have bee deleted
 */
        if (home->recycledNodeHeapSize > 0) {
            free(home->recycledNodeHeap);
            home->recycledNodeHeap = (int *)NULL;
        }

/*
 *      Free any memory associated with the list of topological
 *      operations distributed to remote domains.
 */
        FreeOpList(home);

/*
 *      Remove all cell2 related arrays needed during collision handling
 */
        if (home->cell2) {
            free(home->cell2);
            home->cell2 = NULL;
        }

        if (home->cell2QentArray) {
            free(home->cell2QentArray);
            home->cell2QentArray = NULL;
        }

/*
 *      Free the buffer used to hold the global cell charge tensor
 *      needed for calculating the far field stresses when the fast
 *      multipole code is not enabled.
 */
        if (home->cellCharge) {
            free(home->cellCharge);
            home->cellCharge = NULL;
        }

/*
 *      Free all memory associated with the domain decomposition
 */
        if (home->decomp != (void *)NULL) {
            FreeDecomp(home, home->decomp);
            home->decomp = (void *)NULL;
        }

/*
 *      Release any memory used for the builtin coarse-grain timers
 */
        if (home->timers) {
            free(home->timers);
            home->timers = NULL;
        }

/*
 *      Release memory used to store the burgers-vector specific
 *      dislocation density values.
 */
        if (home->param->partialDisloDensity) {
            free(home->param->partialDisloDensity);
            param->partialDisloDensity = NULL;
        }

/*
 *      Free some arrays used by the fast multipole algorithm.
 */
#ifdef ESHELBY
        if ((home->param->fmEnabled) || (home->param->eshelbyfmEnabled)) 
#else
        if  (home->param->fmEnabled)
#endif
        {
            if (home->glPositions) { free(home->glPositions); home->glPositions=0; }
            if (home->glWeights  ) { free(home->glWeights  ); home->glWeights  =0; }
        }

/*
 *      Free arrays used in control and data file parsing
 */
        if (home->ctrlParamList != (ParamList_t *)NULL) {
            if (home->ctrlParamList->paramCnt > 0) {
                free(home->ctrlParamList->varList);
            }
            free(home->ctrlParamList);
            home->ctrlParamList = (ParamList_t *)NULL;
        }

        if (home->dataParamList != (ParamList_t *)NULL) {
            if (home->dataParamList->paramCnt > 0) {
                free(home->dataParamList->varList);
            }
            free(home->dataParamList);
            home->dataParamList = (ParamList_t *)NULL;
        }

/*
 *      Remove the burgers vector list and associated plane list
 *      if they were allocated.
 */
        if (home->burgData.burgList != (real8 (*)[3])NULL) {
            free(home->burgData.burgList);
            home->burgData.burgList = (real8 (*)[3])NULL;
        }

        if (home->burgData.planeList != (real8 (*)[3])NULL) {
            free(home->burgData.planeList);
            home->burgData.planeList = (real8 (*)[3])NULL;
        }

        if (home->burgData.numPlanesPerBurg != (int *)NULL) {
            free(home->burgData.numPlanesPerBurg);
            home->burgData.numPlanesPerBurg = (int *)NULL;
        }

        if (home->burgData.numGlissilePlanesPerBurg != (int *)NULL) {
            free(home->burgData.numGlissilePlanesPerBurg);
            home->burgData.numGlissilePlanesPerBurg = (int *)NULL;
        }

        if (home->burgData.planeType != (int *)NULL) {
            free(home->burgData.planeType);
            home->burgData.planeType = (int *)NULL;
        }

        if (home->burgData.burgFirstPlaneIndex != (int *)NULL) {
            free(home->burgData.burgFirstPlaneIndex);
            home->burgData.burgFirstPlaneIndex = (int *)NULL;
        }

        if (home->burgData.splinterableBurgList != (SplinterableBurg_t *)NULL) {
            free(home->burgData.splinterableBurgList);
            home->burgData.splinterableBurgList = (SplinterableBurg_t *)NULL;
        }

/*
 *      Destroy the lock used for synchronizing access to the
 *      segment/particle intersection list
 */
        DESTROY_LOCK(&home->segpartIntersectListLock);

/*
 *      If an array of ghost node pointers had been allocated, free it.
 */
        if (home->ghostNodeListSize > 0) {
            free(home->ghostNodeList);
            home->ghostNodeList = (Node_t **)NULL;
            home->ghostNodeListSize = 0;
        }

#ifdef ANISOTROPIC
/*
 *      Free the <Fq> tables and associated indices lists
 */
        if (home->anisoVars.FqReal) {
            free(home->anisoVars.FqReal);
            home->anisoVars.FqReal = (real8 *)NULL;
        }

        if (home->anisoVars.sphericalHarmonicsRemKept) {
            free(home->anisoVars.sphericalHarmonicsRemKept);
            home->anisoVars.sphericalHarmonicsRemKept = (short *)NULL;
        }

        if (home->anisoVars.sphericalHarmonicsRelKept) {
            free(home->anisoVars.sphericalHarmonicsRelKept);
            home->anisoVars.sphericalHarmonicsRelKept = (short *)NULL;
        }

        if (home->anisoVars.sphericalHarmonicsRe18Kept) {
            free(home->anisoVars.sphericalHarmonicsRe18Kept);
            home->anisoVars.sphericalHarmonicsRe18Kept = (short *)NULL;
        }

        if (home->anisoVars.sphericalHarmonicsReiKept) {
            free(home->anisoVars.sphericalHarmonicsReiKept);
            home->anisoVars.sphericalHarmonicsReiKept = (short *)NULL;
        }

        if (home->anisoVars.FqImag) {
            free(home->anisoVars.FqImag);
            home->anisoVars.FqImag = (real8 *)NULL;
        }

        if (home->anisoVars.sphericalHarmonicsImmKept) {
            free(home->anisoVars.sphericalHarmonicsImmKept);
            home->anisoVars.sphericalHarmonicsImmKept = (short *)NULL;
        }

        if (home->anisoVars.sphericalHarmonicsImlKept) {
            free(home->anisoVars.sphericalHarmonicsImlKept);
            home->anisoVars.sphericalHarmonicsImlKept = (short *)NULL;
        }

        if (home->anisoVars.sphericalHarmonicsIm18Kept) {
            free(home->anisoVars.sphericalHarmonicsIm18Kept);
            home->anisoVars.sphericalHarmonicsIm18Kept = (short *)NULL;
        }

        if (home->anisoVars.sphericalHarmonicsImiKept) {
            free(home->anisoVars.sphericalHarmonicsImiKept);
            home->anisoVars.sphericalHarmonicsImiKept = (short *)NULL;
        }

        if (home->anisoVars.coreFqReal) {
            free(home->anisoVars.coreFqReal);
            home->anisoVars.coreFqReal = (real8 *)NULL;

            free(home->anisoVars.coreFqImag);
            home->anisoVars.coreFqImag = (real8 *)NULL;
        }
#endif

        free(home->param);
        home->param = NULL;

/*
 *      Free buffers allocated for MPI reduction operations.
 */
        if (home->localReduceBuf != (int *)NULL) {
            free(home->localReduceBuf);
            home->localReduceBuf = (int *)NULL;
        }

        if (home->globalReduceBuf != (int *)NULL) {
            free(home->globalReduceBuf);
            home->globalReduceBuf = (int *)NULL;
        }

        free(home);

        return;
}


void ParadisFinish(Home_t *home)
{
        int maxMem;

        if (home->myDomain == 0) printf("ParadisFinish\n");

/*
 *      If the ARKode timestep integrator is being used, we need to
 *      invoke a special cleanup function for ARKode.
 */
#ifdef USE_ARKODE
        if (StrEquiv(home->param->timestepIntegrator, "ARKODE")) {
/*
            StatsARKode();
*/
            FreeARKode(home);
        }
#endif

/*
 *      Generate all type of output appropriate at program completion
 *      and write out the timer values 
 */
        GenerateOutput(home, STAGE_TERM);
        TimerPrint(home);

#ifndef NO_XWINDOW
#ifndef NO_THREAD
        if (WinAlive()) {
            if (home->myDomain == 0) printf("sleeping ...\n");
            Sleep();
        }
#else
        if (WinAlive()) {
            if (home->myDomain == 0) {
                printf("sleeping ...\n");
                WinRoutine();
            }
        }
#endif
        if (home->myDomain == 0) {
            WinSemRemove();
        }
#endif /* NO_XWINDOW */
    
#ifdef _ARLFEM
        FEM_ReleaseMemory(home);
#endif

/*
 *      Free any communicators and groups that may have been created
 *      for collective intra-node operations and collective operations
 *      involving the 'first' task on each node.  This MUST be done
 *      before MPI_Finalize() is called.
 */
        FreeNodeTaskGroups(home);

/*
 *      Print out memory usage info for processor zero
 */
        if (home->myDomain == 0) {
            Meminfo(&maxMem);
            printf("Estimated memory usage (proc 0): %dk bytes\n", maxMem);
        }

#ifdef PARALLEL
#ifdef USE_SCR
        SCR_Finalize();
#endif
        MPI_Finalize();
#endif

        ReleaseMemory(home);

        return;
}
