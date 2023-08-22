/**************************************************************************
 *
 *      Module:       ParadisThread.c
 *      Description:  This modules contains functions needed to support
 *                    threading within a single domain.
 *
 *      Functions:
 *          GetThreadIterationIndices()
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

/*------------------------------------------------------------------------
 *
 *      Function:    GetThreadIterationIndices
 *      Description: Given the number of iterations for a loop,
 *                   determine the optimal way to distribute the
 *                   loop iterations among the current thread pool
 *                   and return to the caller the loop indices
 *                   that should be used by this thread.  If the
 *                   code is not in a threaded section, this function
 *                   behaves as if it is dealing with a thread pool
 *                   of 1 thread, and the returned indices will cover
 *                   all iterations of the loop.
 *
 *                   i.e. the indices returned are such that the the caller
 *                   should set up the 'for' loop like this:
 *
 *                       for (i = startIndex, i < endIndex; i++)
 *
 *      Arguments:
 *          loopCount  number of loop iterations to distribute among threads
 *          threadID   Location in which to return current thread ID
 *          startIndex Location in which to return the starting loop
 *                     index for this thread
 *          endIndex   Location in which to return the ending loop
 *                     index for this thread.
 *
 *----------------------------------------------------------------------*/
void GetThreadIterationIndices(int loopCount, int *threadID, int *startIndex,
                               int *endIndex)
{
#ifdef _OPENMP
        int numThreads, minIterPerThread, numAtMaxIter;

        *threadID   = omp_get_thread_num();
        numThreads  = omp_get_num_threads();

        minIterPerThread = loopCount / numThreads;
        numAtMaxIter = loopCount - (numThreads * minIterPerThread);

        *startIndex = *threadID * minIterPerThread +
                      MIN(*threadID, numAtMaxIter);
        *endIndex   = *startIndex + minIterPerThread +
                      (*threadID < numAtMaxIter ? 1 : 0);
#else
        *threadID   = 0;
        *startIndex = 0;
        *endIndex   = loopCount;
#endif

        return;
}
