#pragma once

#ifndef _MPI_PORTABILITY_H
#define _MPI_PORTABILITY_H

#ifdef PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
// If MPI is not actively part of the current build, the following
// empty macros are defined to eliminate the need for placing them
// within conditionals.
//
// Note that currently we only define the abort and barrier calls as
// empty macros. I will activate the remaining macros as I can validate
// their use within the code base. -bmw-
//----------------------------------------------------------------------------------------

#ifndef PARALLEL

#define MPI_Abort(a,b)
#define MPI_Barrier(a)

// #define MPI_Allgather(a,b,c,d,e,f,g)
// #define MPI_Allreduce(a,b,c,d,e,f)
// #define MPI_Bcast(a,b,c,d,e)
// #define MPI_Comm_create(a,b,c)
// #define MPI_Comm_free(a)
// #define MPI_Comm_group(a,b)
// #define MPI_Comm_rank(a,b)
// #define MPI_Finalize()
// #define MPI_Gather(a,b,c,d,e,f,g,h)
// #define MPI_Group_free(a)
// #define MPI_Group_incl(a,b,c,d)
// #define MPI_Init(a,b)
// #define MPI_Irecv(a,b,c,d,e,f,g)
// #define MPI_Isend(a,b,c,d,e,f,g)
// #define MPI_Pack(a,b,c,d,e,f,g)
// #define MPI_Pack_size(a,b,c,d)
// #define MPI_Recv(a,b,c,d,e,f,g)
// #define MPI_Reduce(a,b,c,d,e,f,g)
// #define MPI_Send(a,b,c,d,e,f)
// #define MPI_Unpack(a,b,c,d,e,f,g)
// #define MPI_Wait(a,b)
// #define MPI_Waitall(a,b,c)
// #define MPI_Waitany(a,b,c,d)

// If an application is not built for MPI, we do not need to initialize or finalize
// the MPI communications subsystems...
//----------------------------------------------------------------------------------------

#define MPI_Init(a,b)
#define MPI_Finalize()

#endif  // !PARALLEL

#endif  // _MPI_PORTABILITY_H
