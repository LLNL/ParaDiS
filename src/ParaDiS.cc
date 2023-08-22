/***************************************************************************
 *
 *  Function    : Main
 *  Description : main routine for ParaDiS simulation
 *
 **************************************************************************/
#include <stdio.h>
#include <time.h>

#ifdef FPES_ON
#include <fpcontrol.h>
#endif

#ifdef SHUTDOWN_SIGNAL
#include <signal.h>
#endif

#ifdef SCOREP_USER_ENABLE
#include <scorep/SCOREP_User.h>
#endif

#ifdef USE_CALIPER
#include <caliper/cali.h>
#endif

#include "mpi_portability.h"

#include "Home.h"
#include "Init.h"

#ifdef SHUTDOWN_SIGNAL
/*
 *      A shutdown signal was defined when the code was compiled, so we
 *      need to include some stuff needed for catching and processing
 *      the signal.
 */

static int receivedShutdownSig = 0;
static int doGlobalShutdown = 0;

/*
 *      Simple signal handler...
 */
void SigHandler(int sigNum)
{
        switch (sigNum) {
/*
 *          When the shutdown signal is received, just set a local flag
 */
            case SHUTDOWN_SIGNAL:
                receivedShutdownSig = 1;
                break;
        }
}


/*------------------------------------------------------------------------
 *
 *      Function:     DoForcedShutdown
 *      Description:  This function will broadcast a flag from task
 *                    zero indicating if the the shutdown signal has
 *                    been received on task zero and return that flag
 *                    to the caller.
 *
 *      Returns:  0 if task zero has not received the shutdown signal
 *                1 if task zero has received the shutdown signal
 *
 *-----------------------------------------------------------------------*/
static int DoForcedShutdown(Home_t *home)
{

        if (home->myDomain == 0) {
            doGlobalShutdown = receivedShutdownSig;
            if (doGlobalShutdown) {
                printf("Received shutdown sig, attempting clean shutdown!\n");
            }
        }

#ifdef PARALLEL
        MPI_Bcast(&doGlobalShutdown, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

        return(doGlobalShutdown);
}
#endif  /* ifdef SHUTDOWN_SIGNAL */


/*------------------------------------------------------------------------
 *
 *      Function:     main
 *      Description:  
 *
 *-----------------------------------------------------------------------*/
int main (int argc, char *argv[])
{
        int     cycleEnd, memSize, initialDLBCycles;
        time_t  tp;
        Home_t  *home;
        Param_t *param;

#ifdef USE_CALIPER
        cali_config_preset("CALI_LOG_VERBOSITY", "0");
        cali_config_preset("CALI_CALIPER_ATTRIBUTE_PROPERTIES", 
                           "loop=process_scope:nested,"
                           "iteration#paradis.cycle=process_scope:nested,"
                           "mpi.function=process_scope:nested");
#endif

#ifdef DEBUG_MEM
/*
 *      If the memory debugging code is enabled, we need to initialize
 *      the local memory management layer before doing ANYTHING else!
 */
        ParadisMemCheckInit();
#endif

#ifdef SHUTDOWN_SIGNAL
/*
 *      A shutdown signal has been defined for the code  so register the
 *      signal to be caught by a signal handler.  If the code fails to
 *      register for the signal, it is not considered a fatal error and
 *      will quietly fail, continuing on without the capability to
 *      accept an external shutdown signal
 */
        struct sigaction action;

        action.sa_flags = SA_RESTART;
        action.sa_handler = SigHandler;
        sigemptyset(&action.sa_mask);

        sigaction(SHUTDOWN_SIGNAL, &action, (struct sigaction *)NULL);

#endif  /* end ifdef SHUTDOWN_SIGNAL */

#ifdef _OPENMP
/*
 *      Explicitly turn off nested OpenMP parallelism and the dynamic
 *      adjustment of the number of parallel threads.  This is to
 *      avoid potential problems later.
 */
        omp_set_nested(0);
        omp_set_dynamic(0);
#endif

/*
 *      On some systems, the getrusage() call made by Meminfo() to get
 *      the memory resident set size does not work properly.  In those
 *      cases, the function will try to return the current heap size 
 *      instead.  This initial call allows meminfo() to get a copy of
 *      the original heap pointer so subsequent calls can calculate the
 *      heap size by taking the diference of the original and current
 *      heap pointers.
 */
        Meminfo(&memSize);

/*
 *      on linux systems (e.g. MCR) if built to have floating point exceptions
 *      turned on, invoke macro to do so
 */
   
#ifdef FPES_ON
        unmask_std_fpes();
#endif


/*
 *      Now do the basic initialization (i.e. init MPI, read restart
 *      files, blah, blah, blah...)
 */
        ParadisInit(&home,argc,argv);
   
/* 
 * Determine the critical stress of a FR source
 */

#ifdef FRCRIT
        if (home->param->fmEnabled) Fatal("Critical stress calcs not implemented for FMM yet.\n");
        Calc_FR_Critical_Stress(argc, argv, &home);
        ParadisFinish(home);
        exit(0);
#endif

        home->cycle      = home->param->cycleStart;

        param            = home->param;
        cycleEnd         = param->cycleStart + param->maxstep;
        initialDLBCycles = param->numDLBCycles;

/*
 *      Perform the needed number (if any) of load-balance-only
 *      steps before doing the main processing loop.  These steps
 *      perform only the minimal amount of stuff needed to
 *      estimate per-process load, move boundaries and migrate
 *      nodes among processsors to get a good initial balance.
 */
        TimerStart(home, INITIALIZE);

#ifdef SCOREP_USER_ENABLE
        SCOREP_RECORDING_OFF();
#endif

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

#ifdef SCOREP_USER_ENABLE
        SCOREP_RECORDING_ON();
#endif

        TimerStop(home, INITIALIZE);

/*
 *      Any time spent doing the initial DLB-only steps should
 *      just be attributed to initialization time, so be sure to
 *      reset the other timers before going into the main
 *      computational loop
 */
        TimerInitDLBReset(home);

/*
 *      The cycle number may have been incremented during the initial
 *      load-balance steps, so reset it to the proper starting
 *      value before entering the main processing loop.
 */
        home->cycle = home->param->cycleStart;

#ifdef USE_CALIPER
        CALI_CXX_MARK_LOOP_BEGIN(cycle, "paradis.cycle");
#endif

        while (home->cycle < cycleEnd) 
        {

#ifdef USE_CALIPER
            CALI_CXX_MARK_LOOP_ITERATION(cycle, home->cycle);
#endif

            ParadisStep(home);

            TimerClearAll(home);

#ifdef SHUTDOWN_SIGNAL
            if (DoForcedShutdown(home)) {
                break;
            }
#endif

/*
 *          If the user specified a simulation time at which to terminate
 *          the simulation and we've reached/exceeded that time, then
 *          terminate the simulation even if we have not run for the
 *          specified number of cycles.
 */
            if ((param->timeEnd > 0.0) && (param->timeNow >= param->timeEnd)) {
                if (home->myDomain == 0) {
                    printf("+++\n"
                           "+++ Requested simulation end time: %.10e seconds\n"
                           "+++ Current simulation time:       %.10e seconds\n"
                           "+++\n"
                           "+++ Terminating simulation now...\n"
                           "+++\n", param->timeEnd, param->timeNow);
                }

                break;
            }
        }

#ifdef USE_CALIPER
        CALI_CXX_MARK_LOOP_END(cycle);
#endif        

/*
 *      Write any final output files, free dynamically allocated
 *      memory, do MPI cleanup, etc.
 */
        ParadisFinish(home);

#ifdef DEBUG_MEM
/*
 *      If the memory debugging code is enabled, we need to terminate
 *      the local memory management layer after doing everything else.
 */
        ParadisMemCheckTerminate();
#endif

        exit(0);
        return(0);
}
