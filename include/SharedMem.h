#pragma once

#ifndef _PDS_SHARED_MEM_H
#define _PDS_SHARED_MEM_H

/****************************************************************************
 *
 *      SharedMem.h  Contains definitions, prototypes, etc. needed
 *                   when using shared memory segments to share
 *                   memory buffers among tasks on the same compute-nodes
 *
 ***************************************************************************/

#include "Home.h"

//
//   Prototypes for functions to manage the shared memory 
//
void AllocSharedMem(Home_t *home, int *bufferSize, void **bufferPtr);
void FreeSharedMem(int *bufferSize, void **bufferPtr);
void SharedMemBarrier(Home_t *home);

//
//   Prototypes for functions to identify tasks co-located on compute-nodes
//   and manage MPI communicators for collective operations specific to
//   the co-located tasks.
//
void FreeNodeTaskGroups(Home_t *home);
void GetNodeTaskGroups(Home_t *home);

#endif  // end ifndef _SharedMem_h
