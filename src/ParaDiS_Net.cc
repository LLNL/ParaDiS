/***************************************************************************
 *  Function    : Main
 *  Description : Hack to original Main routine to handle analysis of the 
 *                dislocation microstructure
 **************************************************************************/

#include <stdio.h>
#include <sys/stat.h>
#include <time.h>

#ifdef FPES_ON
#include <fpcontrol.h>
#endif

#include "mpi_portability.h"

#include "Home.h"
#include "Init.h"
#include "DislocationNetwork.h"

// Signal_Abort()
// 
// The following function is included to intercept external signal handling
// to effect a graceful shutdown.  It's still not clear how graceful this will
// be, but it's an attempt.
//--------------------------------------------------------------------------------

#ifndef PARALLEL
#define MPI_Comm_rank(a,b)
#define MPI_Allreduce(a,b,c,d,e,f)
#endif

#ifdef SIGNALS_ENABLED

#include <signal.h>

void Signal_Abort (int indx)
{
   int rank=0; 

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if (rank==0) 
   { 
      switch(indx)
      {
         case (SIGINT ) : { printf("paradis-abort - received signal SIGINT\n" ); break; }
         case (SIGTERM) : { printf("paradis-abort - received signal SIGTERM\n"); break; }
         default        : { printf("paradis-abort - received signal %d\n",indx); break; }
      }
   }

   ParadisAbort_Set(1);
}

#endif

/*------------------------------------------------------------------------
 *      Function:     main
 *      Description:  
 *-----------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
        int     i, memSize, initialDLBCycles;
        time_t  tp;
        Home_t  *home;
        Param_t *param;

#ifdef SIGNALS_ENABLED
        if (signal(SIGINT ,Signal_Abort) == SIG_ERR) { printf("Warning - can't catch SIGINT\n" ); }
        if (signal(SIGTERM,Signal_Abort) == SIG_ERR) { printf("Warning - can't catch SIGTERM\n"); }
#endif

#ifdef DEBUG_MEM
        // If the memory debugging code is enabled, we need to initialize
        // the local memory management layer before doing ANYTHING else!

        ParadisMemCheckInit();
#endif

#ifdef _OPENMP
        // Explicitly turn off nested OpenMP parallelism and the dynamic
        // adjustment of the number of parallel threads.  This is to
        // avoid potential problems later.

        omp_set_nested(0);
        omp_set_dynamic(0);
#endif

        // On some systems, the getrusage() call made by Meminfo() to get
        // the memory resident set size does not work properly.  In those
        // cases, the function will try to return the current heap size 
        // instead.  This initial call allows meminfo() to get a copy of
        // the original heap pointer so subsequent calls can calculate the
        // heap size by taking the diference of the original and current
        // heap pointers.

        Meminfo(&memSize);

        // on linux systems (e.g. MCR) if built to have floating point exceptions
        // turned on, invoke macro to do so
   
#ifdef FPES_ON
        unmask_std_fpes();
#endif

        // Now do the basic initialization (i.e. init MPI, read restart
        // files, blah, blah, blah...)

        ParadisInit(&home,argc,argv);

        home->cycle      = home->param->cycleStart;

        param            = home->param;
        initialDLBCycles = param->numDLBCycles;

        // Perform the needed number (if any) of load-balance-only
        // steps before doing the main processing loop.  These steps
        // perform only the minimal amount of stuff needed to
        // estimate per-process load, move boundaries and migrate
        // nodes among processsors to get a good initial balance.

        TimerStart(home, INITIALIZE);

        if ((home->myDomain == 0) && (initialDLBCycles != 0)) {
            time(&tp);
            printf("  +++ Beginning %d load-balancing steps at %s",
                   initialDLBCycles, asctime(localtime(&tp)));
        }

        while (param->numDLBCycles > 0) {
            ParadisStep(home);
            home->cycle++;
            param->numDLBCycles--;
        }

        if ((home->myDomain == 0) && (initialDLBCycles != 0)) {
            time(&tp);
            printf("  +++ Completed load-balancing steps at %s",
                   asctime(localtime(&tp)));
        }

        TimerStop(home, INITIALIZE);

        // Any time spent doing the initial DLB-only steps should
        // just be attributed to initialization time, so be sure to
        // reset the other timers before going into the main
        // computational loop

        TimerInitDLBReset(home);

        // The cycle number may have been incremented during the initial
        // load-balance steps, so reset it to the proper starting
        // value before entering the main processing loop.

        home->cycle = home->param->cycleStart;

        // Post-process data : Dislocation network
        printf("\n\nInitial verification\n");
        VerifyBurgersVectors(home);
        printf("Initial verification is good\n\n");
        
        printf("Calling Remove2Nodes\n");
        Remove2Nodes(home);
        VerifyBurgersVectors(home);
        printf("Remove2Nodes verification is good\n\n");
        
        printf("Calling SplitZeroPairs\n");
        SplitZeroPairs(home);
        printf("Checking Burgers vectors ....\n");
        VerifyBurgersVectors(home);
        printf("SplitZeroPairs verification is good\n\n");
              
        printf("Calling Remove2Nodes\n");
        Remove2Nodes(home);
        VerifyBurgersVectors(home);
        printf("Remove2Nodes verification is good\n\n");
        
        int three=0, four=0, more=0;
        Node_t *node;
        
        for (i = 0; i < home->newNodeKeyPtr; i++) {
           
           if ((node = home->nodeKeys[i]) == (Node_t *)NULL) 
              continue;
           else
           {
              if (node->numNbrs == 2)
                 printf("Nodes with two connections remains....\n");
              else if (node->numNbrs == 3) 
                 three++;         
              else if (node->numNbrs == 4) 
                 four++;
              else
                 more++;
           }
        }
        
        printf("\n\n%d %d %d\n",three,four,more);
        JunctionsLength(home);


        // Write any final output files, free dynamically allocated
        // memory, do MPI cleanup, etc.

        ParadisFinish(home);

#ifdef DEBUG_MEM
        // If the memory debugging code is enabled, we need to terminate
        // the local memory management layer after doing everything else.
 
        ParadisMemCheckTerminate();
#endif

        exit(0); 
        return(0);
}
