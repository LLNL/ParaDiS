/**************************************************************************
 *
 *      Module:       MemCheck.c
 *      Description:  This module contains wrappers for the basic C
 *                    memory allocation functions.  When memory
 *                    debugging is enabled, these wrappers provide
 *                    a mechanism for adding guard blocks around
 *                    allocated buffers and tracking allocated buffers
 *                    for debugging purposes.
 *
 *      Included functions:
 *               ParadisCalloc()
 *               ParadisFree()
 *               ParadisMalloc()
 *               ParadisMemCheck()
 *               ParadisRealloc()
 *
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>

#include "ParadisThread.h"

/*
 *      If threading is enabled, define a lock that will be used to
 *      synchronize access to the data structures used in the local
 *      memory management functions.
 */
#ifdef _OPENMP
omp_lock_t memcheckLock;
#endif

/*
 *      Define the number of bytes that will be prepended
 *      and appended to the memory requested by the caller.
 *      these header and trailer bytes will be set to a
 *      known pattern in order to more easily detect memory
 *      corruption.
 */
#define PARADIS_MEM_HEADER_LEN  32
#define PARADIS_MEM_TRAILER_LEN 32

#define PARADIS_MAX_NAMELEN 32
#define PARADIS_MEMBLOCK_INCREMENT 100

/*
 *      The application headers are purposely NOT included
 *      in this module, so we explicitly prototype any
 *      functions that may be called.
 */
void  Fatal(const char *format, ...);

void  ParadisMemCheck(void);
void  ParadisFree    (const char *fname, const char *func, const int ln, void *ptr);
void *ParadisMalloc  (const char *fname, const char *func, const int ln, const size_t size);
void *ParadisCalloc  (const char *fname, const char *func, const int ln, const size_t numElem, const size_t size);
void *ParadisRealloc (const char *fname, const char *func, const int ln, void *ptr, const size_t size);

/*
 *      Define a structure of info stored for each allocated
 *      block of memory.
 */
typedef struct {
    char  fileName[PARADIS_MAX_NAMELEN];  ///< module from which the memory block was allocated
    char  funcName[PARADIS_MAX_NAMELEN];  ///< function from which the memory block was allocated
    int   lineNum ;                       ///< source code line number at which memory block was allocated
    int   userSize;                       ///< size (in bytes) of the portion of the allocated memory buffer available to the caller
    void *realAddr;                       ///< true address of the allocated memory block
    void *userAddr;                       ///< memory address returned to caller
    char *header  ;                       ///< pointer to header prepended to memory buffer
    char *trailer ;                       ///< pointer to trailer appended to the memory buffer
} MemBlock_t;


/*
 *      And define some variables used by all the memory debugging
 *      functions for managing the array of memory blocks.
 */
static MemBlock_t *memBlocks          = (MemBlock_t *) NULL;
static int         memBlocksUsed      = 0;

#ifdef DEBUG_MEM
static int         memBlocksAllocated = 0;
#endif

/*
 *      If set, the following flag forces all the following memory
 *      management functions to explicitly check the consistency
 *      of the entire memory block array every time the functions
 *      are invoked.
 *
 *      This is extremely expensive, and may be best to set via
 *      the debugger only when necessary.  Hence, the default is
 *      to leave the flag off.
 */
int        doMemCheck = 0;


/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisMemCheckInit
 *      Description:    This function is used to do any necessary
 *                      initialization for the local memory management
 *                      functions.  
 *
 *                      WARNING! This function MUST be called before the
 *                      application invokes any of the other local memory
 *                      management functions!
 *
 *-------------------------------------------------------------------------*/
void ParadisMemCheckInit(void)
{
    INIT_LOCK(&memcheckLock);

    return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisMemCheckTerminate
 *      Description:    This function is used to do any necessary cleanup
 *                      after the local memory management functions are
 *                      no longer needed.  
 *
 *                      WARNING! This function should be called only during
 *                      application termination and once called, the other
 *                      local memory management functions MUST NOT be invoked!
 *
 *-------------------------------------------------------------------------*/
void ParadisMemCheckTerminate(void)
{
    DESTROY_LOCK(&memcheckLock);

    return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisMemCheck
 *      Description:    This function goes through the entire array of
 *                      allocated memory blocks checking for data corruption
 *                      by verifying that all memory block headers and
 *                      trailers contain the expected known sequence of
 *                      bytes
 *
 *-------------------------------------------------------------------------*/
void ParadisMemCheck(void)
{
    char  hdr [PARADIS_MEM_HEADER_LEN ];
    char  trlr[PARADIS_MEM_TRAILER_LEN];

    memset(hdr , 'H', PARADIS_MEM_HEADER_LEN );
    memset(trlr, 'T', PARADIS_MEM_TRAILER_LEN);

    LOCK(&memcheckLock);

    for (int i=0; i<memBlocksUsed; i++) 
    {
        if (memcmp(memBlocks[i].header, hdr, PARADIS_MEM_HEADER_LEN) != 0)
            Fatal("ERROR! memBlock[%d] header corruption!\n", i);
        
        if (memcmp(memBlocks[i].trailer, trlr, PARADIS_MEM_TRAILER_LEN) != 0)
            Fatal("ERROR! memBlock[%d] trailer corruption!\n", i);
        
    }

    UNLOCK(&memcheckLock);

    return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisFree
 *      Description:    If memory debugging is not enabled, this function
 *                      simply frees the specified memory block and 
 *                      returns to the caller.  Otherwise, this function
 *                      removes the specified memory block from the array
 *                      of allocated memory blocks, then uses the system
 *                      functions to free the memory.
 * 
 *                      If the specified memory address is not found in
 *                      the array of allocated memory blocks, a warning
 *                      message is generated before the memory is freed.
 *
 *      Arguments:
 *          fileName    Name of the module from which this function was invoked.
 *          lineNum     Line number in <fileName> from which this function was invoked.
 *          ptr         pointer to memory block to be freed
 *
 *-------------------------------------------------------------------------*/
void ParadisFree(const char *fileName, const char *funcName, const int lineNum, void *ptr)
{
#ifndef DEBUG_MEM
        free(ptr);
#else

/*
 *      If the incoming pointer is NULL, no need to free anything
 */
        if (ptr == (void *)NULL) {
            return;
        }

        LOCK(&memcheckLock);

/*
 *      Locate the specified memory addess in the array
 */
        int i=0;
        int found=0;
        for (i=0; i<memBlocksUsed; i++) {
            if (ptr == memBlocks[i].userAddr) {
                found = 1;
                break;
            }
        }

        if (!found) {
            printf("WARNING! Unmatched address in free(%p) at %s::%s(), line %d\n", ptr, fileName, funcName, lineNum);
            free(ptr);
            UNLOCK(&memcheckLock);
            return;
        }

/*
 *      Free the memory, and remove the block from the array
 */
        free(memBlocks[i].realAddr);

        memset(&memBlocks[i], 0, sizeof(MemBlock_t));
        memBlocksUsed--;

        if ((memBlocksUsed > 0) && (memBlocksUsed != i)) {
            memcpy(&memBlocks[i], &memBlocks[memBlocksUsed], sizeof(MemBlock_t));
            memset(&memBlocks[memBlocksUsed], 0, sizeof(MemBlock_t));
        }

        UNLOCK(&memcheckLock);

/*
 *      If required, do a full consistency check on the array of
 *      memory blocks.
 */
        if (doMemCheck) {
            ParadisMemCheck();
        }

#endif  /* ifndef DEBUG_MEM */

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisMalloc
 *      Description:    If memory debugging is not enabled, this function
 *                      simply allocates a memory block of the specified
 *                      size and returns to the caller.  Otherwise,
 *                      this function allocates a memory buffer for the
 *                      caller, add the gurad blocks and saves information
 *                      about the allocation in the memory block array
 *                      before returning to the caller.
 * 
 *      Arguments:
 *          fileName    Name of the module from which this function was invoked.
 *          lineNum     Line number in <fileName> from which this function was invoked.
 *          size        size in bytes of the memory buffer needed by the caller.
 *
 *-------------------------------------------------------------------------*/
void *ParadisMalloc(const char *fileName, const char *funcName, const int lineNum, const size_t size)
{
        void   *userMemAddr = (void *)NULL;

/*
 *      If the specified size is zero, don't attempt an allocation
 */
        if (size <= 0) {
            return(userMemAddr);
        }

#ifndef DEBUG_MEM
        
        userMemAddr = (void *)malloc(size);

        if (userMemAddr == (void *)NULL) {
            Fatal("Malloc error at %s::%s() line %d, size = %d\n", fileName, funcName, lineNum, size);
        }
#else

        LOCK(&memcheckLock);

/*
 *      If we need more space in the memory block array, reallocate
 *      the array with more space.
 */
        if (memBlocksUsed == memBlocksAllocated) 
        {
            size_t newSize = (memBlocksAllocated + PARADIS_MEMBLOCK_INCREMENT) * sizeof(MemBlock_t);

            memBlocks = (MemBlock_t *)realloc(memBlocks, newSize);

            if (memBlocks == (MemBlock_t *)NULL) {
                Fatal("no memory for %d malloc at %s::%s(), line %d", size, fileName, funcName, lineNum);
            }
            memBlocksAllocated += PARADIS_MEMBLOCK_INCREMENT;
        }

/*
 *      Allocate memeory for the caller with enough extra storage for
 *      a header and trailer to be prepended/appended.  Save the
 *      module and line number at which this request was generated,
 *      fill in the header and trailer, save the true address of
 *      the allocated memory, and the adjusted address which will
 *      be returned to the caller.
 */
        memset(&memBlocks[memBlocksUsed], 0, sizeof(MemBlock_t));

        size_t adjustedSize = size + PARADIS_MEM_HEADER_LEN + PARADIS_MEM_TRAILER_LEN;

        strncpy(memBlocks[memBlocksUsed].fileName, fileName, PARADIS_MAX_NAMELEN-1);
        strncpy(memBlocks[memBlocksUsed].funcName, funcName, PARADIS_MAX_NAMELEN-1);
        memBlocks[memBlocksUsed].lineNum  = lineNum;
        memBlocks[memBlocksUsed].userSize = size;
        memBlocks[memBlocksUsed].realAddr = (void *)malloc(adjustedSize);

        if (memBlocks[memBlocksUsed].realAddr == (void *)NULL) {
            Fatal("Malloc error at %s::%s() line %d, adjusted size = %d\n", fileName, funcName, lineNum, adjustedSize);
        }

        memBlocks[memBlocksUsed].userAddr = (char *)memBlocks[memBlocksUsed].realAddr + PARADIS_MEM_HEADER_LEN;
        memBlocks[memBlocksUsed].header   = (char *)memBlocks[memBlocksUsed].realAddr;
        memBlocks[memBlocksUsed].trailer  = (char *)memBlocks[memBlocksUsed].realAddr + PARADIS_MEM_HEADER_LEN + size;

        memset(memBlocks[memBlocksUsed].header , 'H', PARADIS_MEM_HEADER_LEN );
        memset(memBlocks[memBlocksUsed].trailer, 'T', PARADIS_MEM_TRAILER_LEN);

        userMemAddr = memBlocks[memBlocksUsed].userAddr;

        memBlocksUsed++;

        UNLOCK(&memcheckLock);

/*
 *      If required, do a full consistency check on the array of
 *      memory blocks.
 */
        if (doMemCheck) {
            ParadisMemCheck();
        }

#endif  /* ifndef DEBUG_MEM */

        return(userMemAddr);
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisCalloc
 *      Description:    This function allocates a zeroed-out memory
 *                      buffer for the caller.  If memory debugging is
 *                      not enabled, this function then simply returns
 *                      to the caller.  Otherwise, it add the guard
 *                      blocks and saves information about the allocation
 *                      in the memory block array before returning.
 * 
 *      Arguments:
 *          fileName    Name of the module from which this function was invoked.
 *          lineNum     Line number in <fileName> from which this function was invoked.
 *          numElem     Number of elements (of size <size>) to be allocated
 *          size        size of each element to be allocated
 *
 *-------------------------------------------------------------------------*/
void *ParadisCalloc(const char *fileName, const char *funcName, const int lineNum, const size_t numElem, const size_t size)
{
        void *userMemAddr = (void *)NULL;

/*
 *      If the specified size is zero, don't attempt an allocation
 */
        if ((numElem * size) <= 0) {
            return((void *)NULL);
        }

#ifndef DEBUG_MEM
        
        userMemAddr = (void *)calloc(numElem, size);

        if (userMemAddr == (void *)NULL) {
            Fatal("Calloc error at %s::%s() line %d, numElem = %d, size = %d\n", fileName, funcName, lineNum, numElem, size);
        }
#else

        LOCK(&memcheckLock);

/*
 *      If we need more space in the memory block array, reallocate
 *      the array with more space.
 */
        if (memBlocksUsed == memBlocksAllocated) {
            size_t newSize;

            newSize = (memBlocksAllocated + PARADIS_MEMBLOCK_INCREMENT) *
                      sizeof(MemBlock_t);
            memBlocks = (MemBlock_t *)realloc(memBlocks, newSize);
            if (memBlocks == (MemBlock_t *)NULL) {
                Fatal("no memory for %d calloc at %s::%s(), line %d", numElem * size, fileName, funcName, lineNum);
            }
            memBlocksAllocated += PARADIS_MEMBLOCK_INCREMENT;
        }

/*
 *      Allocate memory for the caller with enough extra storage for
 *      a header and trailer to be prepended/appended.  Save the
 *      module and line number at which this request was generated,
 *      fill in the header and trailer, save the true address of
 *      the allocated memory, and the adjusted address which will
 *      be returned to the caller.
 */
        memset(&memBlocks[memBlocksUsed], 0, sizeof(MemBlock_t));

        size_t adjustedSize =   (size * numElem)
                              + PARADIS_MEM_HEADER_LEN
                              + PARADIS_MEM_TRAILER_LEN;

        strncpy(memBlocks[memBlocksUsed].fileName, fileName, PARADIS_MAX_NAMELEN-1);
        strncpy(memBlocks[memBlocksUsed].funcName, funcName, PARADIS_MAX_NAMELEN-1);
        memBlocks[memBlocksUsed].lineNum  = lineNum;
        memBlocks[memBlocksUsed].userSize = size * numElem;
        memBlocks[memBlocksUsed].realAddr = (void *)calloc(1, adjustedSize);

        if (memBlocks[memBlocksUsed].realAddr == (void *)NULL) {
            Fatal("Calloc error at %s::%s() line %d, adjusted size = %d\n", fileName, funcName, lineNum, adjustedSize);
        }

        memBlocks[memBlocksUsed].userAddr = (char *)memBlocks[memBlocksUsed].realAddr + PARADIS_MEM_HEADER_LEN;
        memBlocks[memBlocksUsed].header   = (char *)memBlocks[memBlocksUsed].realAddr;
        memBlocks[memBlocksUsed].trailer  = (char *)memBlocks[memBlocksUsed].realAddr + PARADIS_MEM_HEADER_LEN + size;

        memset(memBlocks[memBlocksUsed].header , 'H', PARADIS_MEM_HEADER_LEN );
        memset(memBlocks[memBlocksUsed].trailer, 'T', PARADIS_MEM_TRAILER_LEN);

        userMemAddr = memBlocks[memBlocksUsed].userAddr;

        memBlocksUsed++;

        UNLOCK(&memcheckLock);

/*
 *      If required, do a full consistency check on the array of
 *      memory blocks.
 */
        if (doMemCheck) {
            ParadisMemCheck();
        }

#endif  /* ifndef DEBUG_MEM */

        return(userMemAddr);
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ParadisRealloc
 *      Description:    This function will reallocate a memory buffer
 *                      to the specified size.  If memory debugging is not
 *                      enabled, it then simply returns to the caller.
 *                      Otherwise, it removes from the memory block array,
 *                      information about the old allocation, adds gurad
 *                      blocks, and adds information about the new allocation
 *                      to the mormy lock array.
 * 
 *      Arguments:
 *          fileName    Name of the module from which this function was invoked.
 *          lineNum     Line number in <fileName> from which this function was invoked.
 *          ptr         pointer to memory block to be reallocated
 *          size        size in bytes of the memory buffer needed by the caller.
 *
 *-------------------------------------------------------------------------*/
void *ParadisRealloc(const char *fileName, const char *funcName, const int lineNum, void *ptr, const size_t size)
{
        void *userAddr = (void *)NULL;

#ifndef DEBUG_MEM

        if (size == 0) {

            if (ptr != (void *)NULL) {
                free(ptr);
            }

            return((void *)NULL);
        }

        userAddr = (void *)realloc(ptr, size);

        if (userAddr == (void *)NULL) {
            Fatal("Realloc error at %s::%s() line %d, addr = %p, size = %d\n", fileName, funcName, lineNum, ptr, size);
        }

#else
        LOCK(&memcheckLock);

        int  i=0, found=0;
/*
 *      If the caller provided an address, locate the associated
 *      entry in the memory block array.
 */
        if (ptr != (void *)NULL) 
        {
            for (i=0; i < memBlocksUsed; i++) {
                if (ptr == memBlocks[i].userAddr) { found=1; break; }
            }

            if (!found && (size != 0)) 
            {
                void *newMem;

                printf("WARNING! Unmatched address in realloc(%p) " "at %s::%s(), line %d\n", ptr, fileName, funcName, lineNum);

                newMem = (void *)realloc(ptr, size);

                if (newMem == (void *)NULL) {
                    Fatal("no memory for %d realloc at %s::%s(), line %d", size, fileName, funcName, lineNum);
                }

                UNLOCK(&memcheckLock);

                return(newMem);
            }
        }

/*
 *      If the requested block size is zero, treat the request as 
 *      a simple free() and return.
 */
        if ((size == 0) && (found)) 
        {
            free(memBlocks[i].realAddr);
            memset(&memBlocks[i], 0, sizeof(MemBlock_t));
            memBlocksUsed--;

/*
 *          Keep the memory block array compact by moving the
 *          last entry into the one we just freed up.
 */
            if ((memBlocksUsed > 0) && (memBlocksUsed != i)) 
            {
                memcpy(&memBlocks[i], &memBlocks[memBlocksUsed], sizeof(MemBlock_t));
                memset(&memBlocks[memBlocksUsed], 0, sizeof(MemBlock_t));
            }

            UNLOCK(&memcheckLock);

            return((void *)NULL);
        }

        if (size == 0) {
            printf("WARNING! Zero length realloc(%p) at %s::%s(), line %d\n", ptr, fileName, funcName, lineNum);
        }

/*
 *      Allocate a new memory block.  If this is to replace an old
 *      block, copy the old data into the new buffer, and free the
 *      old one.
 *
 *      NOTE: Need to release the memcheck lock before calling
 *            ParadisMalloc() so it can acquire the lock, then
 *            reacquire the lock afterward.
 */
        UNLOCK(&memcheckLock);

        userAddr = ParadisMalloc(fileName, funcName, lineNum, size);

        LOCK(&memcheckLock);

        if (found) 
        {
            size_t copySize;

/*
 *          We had to release the lock and reaacquire it above, which means
 *          the <memBlock> array may have been changed by another thread.
 *          Have to re-locate the block associated with the user's realloc
 *          request.
 */
            found=0;

            for (i=0; i < memBlocksUsed; i++) {
                if (ptr == memBlocks[i].userAddr) {
                    found = 1;
                    break;
                }
            }

/*
 *          If we couldn't find the block this time, something is screwy!
 */
            if (!found) {
                Fatal("Couldn't relocate block for realloc(%p) at %s::%s(), line %d", ptr, fileName, funcName, lineNum);
            }

            if (size < (size_t) memBlocks[i].userSize) { copySize = size; } 
            else                                       { copySize = memBlocks[i].userSize; }

            memcpy(userAddr, memBlocks[i].userAddr, copySize);

            free(memBlocks[i].realAddr);
            memset(&memBlocks[i], 0, sizeof(MemBlock_t));
            memBlocksUsed--;

            if ((memBlocksUsed > 0) && (memBlocksUsed != i)) 
            {
                memcpy(&memBlocks[i], &memBlocks[memBlocksUsed], sizeof(MemBlock_t));
                memset(&memBlocks[memBlocksUsed], 0, sizeof(MemBlock_t));
            }
        }

        UNLOCK(&memcheckLock);

/*
 *      If required, do a full consistency check on the array of
 *      memory blocks.
 */
        if (doMemCheck) { ParadisMemCheck(); }

#endif  /* ifndef DEBUG_MEM */

        return(userAddr);
}
