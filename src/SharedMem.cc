/*---------------------------------------------------------------------------
 *
 *   Module:      SharedMem.cc
 *
 *   Description: Contains functions to allocate shared memory buffers,
 *                memory map them into local memory, and synchronize
 *                access via a barrier.
 *
 *   Includes public functions:
 *       SharedMemBarrier()
 *       FreeSharedMem()
 *       AllocSharedMem()
 *
 *-------------------------------------------------------------------------*/

#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"

#if defined PARALLEL && defined SHARED_MEM
#include <time.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#endif


/*---------------------------------------------------------------------------
 *
 *   Function:    SharedMemBarrier
 *
 *   Description: Simple function used to have all MPI tasks on a
 *                compute-node synchronize at a barrier.  Use in
 *                conjunction with share memory buffers to keep tasks
 *                that have mapped a shared memory object from accessing
 *                it until the creating process has initialized the
 *                buffer contents.
 *
 *-------------------------------------------------------------------------*/
void SharedMemBarrier(Home_t *home)
{

//
//  If we are not running in parallel or do not have shared memory enabled,
//  there's no need to do anything.
//
#if defined PARALLEL && defined SHARED_MEM
    MPI_Barrier(home->commNode);
#endif

    return;
}


/*---------------------------------------------------------------------------
 *
 *   Function:    FreeSharedMem
 *
 *   Description: This function attempts to unmap the shared memory buffer
 *                previously memory mapped via mmap(),
 *
 *                Since shared memory buffers are only allocated when the
 *                simulation is compiled for parallel execution with shared
 *                memory enabled, if either is not enabled, the associated
 *                memory will have been allocated with the standard memory
 *                allocation mechanism, so the memory will be freed in a
 *                corresponding manner.
 *
 *                NOTE: We need the proper buffer size to unmap the shared memory
 *                buffer.  However, typically, the shared buffers won't be
 *                deallocated until the simulation is ending.  So, if the
 *                size of the shared memory block is not readily available
 *                at the time of this call, we'll just skip the unmapping and
 *                NULL out the pointer.  This leaves a system reference to the
 *                shared memory buffer which means it won't go away immediately.
 *                When the simulation ends, however, the reference will be
 *                automatically deleted and the shared memory buffer will
 *                be freed up at that time.
 *
 *
 *   Parameters:
 *       IN:     bufferSize  Size of the associated buffer in bytes.  If 
 *                           not provided and shared memory is enabled, the
 *                           shared memory segment is NOT unmapped, but the
 *                           contents of bufferPtr are zeroed.
 *       IN/OUT: bufferPtr   Address of the pointer to the buffer to be freed.
 *                           The contents of bufferPtr will be zeroed before
 *                           before control is returned to the caller.
 *
 *-------------------------------------------------------------------------*/
void FreeSharedMem(int *bufferSize, void **bufferPtr)
{

#if defined PARALLEL && defined SHARED_MEM

//
//  If we know the size of the shared memory buffer, try to unmap it
//  now.
//
    if (bufferSize > 0) {
        munmap(*bufferPtr, *bufferSize);
    }

#else  // Either PARALLEL or SHARED_MEM is not defined
//
//  Use of shared memory buffers is not enabled, so simply free
//  the specified buffer.
//
    free(*bufferPtr);

#endif  // end if defined PARALLEL && defined SHARED_MEM

    *bufferPtr = (void *)NULL;

    return;
}


/*---------------------------------------------------------------------------
 *
 *   Function:    AllocSharedMem
 *   Description: Attempt to memory map a shared memory buffer into the
 *                task's local memory.  The first MPI task on the node
 *                allocates the shared memory segment before all other
 *                tasks open and memory map it.
 *
 *                If shared memory support was not enabled at compile
 *                time, a buffer will be allocated using the standard
 *                memory allocation mechanism.
 *
 *      WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!
 *      WARNING!                                                       WARNING!
 *      WARNING!      When a shared memory buffer is created/opened,   WARNING!
 *      WARNING!      only the creating task has write access.         WARNING!
 *      WARNING!      All other tasks have read-only access.  Should   WARNING!
 *      WARNING!      a task with read-only access attempt to write    WARNING!
 *      WARNING!      to the shared memory, a segmentation fault       WARNING!
 *      WARNING!      will occur!                                      WARNING!
 *      WARNING!                                                       WARNING!
 *      WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!
 *
 *   Parameters:
 *       IN:     bufferSize  Size (in bytes) of the shared memory buffer
 *                           to be allocated.
 *       IN/OUT: bufferPtr   Location in which to return to the caller a
 *                           pointer to shared memory buffer allocated.
 *
 *-------------------------------------------------------------------------*/
void AllocSharedMem(Home_t *home, int *bufferSize, void **bufferPtr)
{
#if defined PARALLEL && defined SHARED_MEM
    int   i;
    int   isFirstTaskOnNode;
    int   thisTask;
    int   sharedMemID;
    int   protection;
    char  sharedBufName[64];


    thisTask          = home->myDomain;
    isFirstTaskOnNode = home->isFirstTaskOnNode;

//
//  The first task on the node creates a shared memory buffer
//  in /dev/shm.
//
    if (isFirstTaskOnNode) {
        int seed = (int)clock();
        int gotBuffer = 0;
        int tryCount = 1;

//
//      Attempt to create a shared memory object with a unique
//      name.  A failure due to an already existing object with
//      the same name or an access error will result in a retry
//      using another name.  Any other error is a serious problem.
//
//      Only try a set number of times before giving up...
//
        while (tryCount < 100) {

//
//          Use a random number to create a unique name for the shared
//          memory object.
//
//          Note: Actual shared memory object will be created with
//          the specified name under /dev/shm 
//
            sprintf(sharedBufName, "/ParaDiS.%016.0f", 1.0e+16 * randm(&seed));

            sharedMemID = shm_open(sharedBufName, O_RDWR | O_CREAT | O_EXCL,
                                   S_IRUSR | S_IWUSR);

//
//          If we got a shared memory object, try to set it to the proper
//          length.  On any error, just abort.
//
            if (sharedMemID >= 0) {
                int rc;

                if (ftruncate(sharedMemID, *bufferSize) == -1) {
                    shm_unlink(sharedBufName);
                    Fatal("Task %d: ftruncate error %d on shared buffer %s",
                          thisTask, errno, sharedBufName);
                }

                gotBuffer = 1;
                break;
            }

//
//          If the shared memory object already existed or we had an
//          access problem, try again.  For any other error, we'll
//          stop trying.
//
            if ((errno == EEXIST) || (errno == EACCES)) {
                ++tryCount;
                continue;
            }

            break;

        }  // end while (tryCount < 100)

//
//      If we did not get a shared memory object of the proper
//      length, abort.
//
        if (gotBuffer == 0) {
            perror("shm_open");
            Fatal("Task %d: Error %d creating shared buffer %s",
                  thisTask, errno, sharedBufName);
        }

    }  // end if (isFirstTaskOnNode) 

//
//  Broadcast the name and size of the shared memory object to
//  the other tasks on the node.  All tasks (other than the first 
//  on the node) then open the existing shared memory object
//
    MPI_Bcast(bufferSize, 1, MPI_INT, 0, home->commNode);

    MPI_Bcast(sharedBufName, sizeof(sharedBufName), MPI_CHAR, 0,
              home->commNode);

    if (!isFirstTaskOnNode) {

        sharedMemID = shm_open(sharedBufName, O_RDONLY, 0);

        if (sharedMemID < 0) {
            Fatal("Task %d: Error %d opening shared buffer %s",
                  thisTask, errno, sharedBufName);
        }
    }

//
//  All tasks have opened the shared memory object.  Now map it
//  to a memory address.  Creating task always gets read/write
//  access, all other tasks on the node get read-only access
//
    if (isFirstTaskOnNode) {
        protection = PROT_READ | PROT_WRITE;
    } else {
        protection = PROT_READ;
    }

    *bufferPtr = mmap((void *)NULL, *bufferSize, protection, MAP_SHARED,
                      sharedMemID, 0);

    if (*bufferPtr == (void *)NULL) {
        Fatal("Task %d: Error %d memory mapping shared buffer %s",
              thisTask, errno, sharedBufName);
    }

//
//  Documentation is unclear whether ftruncate() will zero
//  memory in the shared memory buffer, so explicitly do
//  so now in the task that created it.
//
    if (isFirstTaskOnNode) {
        char *cPtr;

        for (cPtr = (char *)*bufferPtr, i = 0; i < *bufferSize; ++i) {
            *cPtr++ = 0;
        }
    }

//
//  Wait until all tasks on the node have opened the shared buffer
//  and memory mapped it.
//
    MPI_Barrier(home->commNode);

//
//  As long as the memory mappings are not removed, there will be
//  system references to the shared object.  This means that
//  we may now unlink the shared memory object.  This will remove
//  the name from /dev/shm, but the shared memory will exist until
//  the last process on the node unmaps the memory or exits.
//
//  This is done to prevent the shared memory object from persisting
//  after the application ends.
//
    if (isFirstTaskOnNode) {
        if (shm_unlink(sharedBufName) < 0) {
            Fatal("Task %d: error %d unlinking shared object %s", 
                  thisTask, errno, sharedBufName);
        }
    }

#else  // Either PARALLEL or SHARED_MEM is not defined
//
//  Use of shared memory buffers is not enabled, so just do a
//  normal memory allocation
//
    *bufferPtr = calloc(1, *bufferSize);
#endif

    return;
}
