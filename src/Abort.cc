//-----------------------------------------------------------------------------------------------
// Module:  Abort() - implementation of the application abort modules
//-----------------------------------------------------------------------------------------------
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Abort.h"

//-----------------------------------------------------------------------------------------------
// 'paradis_abort' is a global variable that will be used to flag an abort condition to
// the rest of the application (and ultimately - the entire MPI process collection).
//
// The variable is declared via an extern in the "Abort.h" file so applications that 
// reference the global abort will need to include the Abort.cc in their builds.
//
// Note that part of the rationale for having this global variable here is that the 
// ParaDiS Fatal() function does not include 'home' as part of the arguments.  This 
// global eliminates the need for having an abort flag in the 'home' data structure.
//-----------------------------------------------------------------------------------------------

int paradis_abort = 0;

// ParadisAbort()
//
// Sets the current abort flag.
//-----------------------------------------------------------------------------------------------

void ParadisAbort_Set (const int abort)   { paradis_abort = abort; }

//-----------------------------------------------------------------------------------------------
// ParadisAbort()
//
// Will return the status of the global abort variable.
//
// If MPI is active, will logically OR all the local abort variables
// into the global abort. If any of the processes are aborting, forces
// ALL the processes to abort.
//-----------------------------------------------------------------------------------------------

int ParadisAbort (void)
{
#ifdef PARALLEL
   int local_abort = paradis_abort;
   MPI_Allreduce(&local_abort, &paradis_abort, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
#endif

   return(paradis_abort);
}

//-----------------------------------------------------------------------------------------------
// ParadisAbort()
//
// Same as above, but will also save the state of home to an output file if the 
// size of the MPI simulation system is small (<1024 processes).
//-----------------------------------------------------------------------------------------------

int ParadisAbort (Home_t *home)
{
#ifdef PARALLEL
   int local_abort = paradis_abort;
   MPI_Allreduce(&local_abort, &paradis_abort, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
#endif

   // If we're aborting and the number of MPI process count is small, save
   // the home data structure to a file for debug.

   if (paradis_abort && home && (home->mpi_size < 1024) )
   {
      char path[256]; 

      mkdir("./abort", S_IRWXU);
      sprintf(path,"./abort/t%04d.home", home->mpi_rank );
      
      //Home_Save(path,home);
   }

   return(paradis_abort);
}

