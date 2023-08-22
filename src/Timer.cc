/***************************************************************************
 *
 *	Module:		Timer.c
 *	Description:	This module contains the various functions needed
 *			to provided a simple timing facility for ParaDiS.
 *
 *	Included functions:
 *
 *		TimerRegister
 *		TimerInit
 *		TimerStart
 *		TimerStop
 *		TimerSave
 *		TimerPrint
 *		TimerAtRestart
 *
 **************************************************************************/

#include <stdio.h>
#include <time.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"

#ifdef SCOREP_USER_ENABLE
#include <scorep/SCOREP_User.h>
#endif

#ifdef USE_CALIPER
#include <caliper/cali.h>
#endif

/*-------------------------------------------------------------------------
 *
 *	Function:	TimerRegister
 *	Description:    Initializes the timing data for the indicated
 *			category and assigns a name to the category
 *	Arguments:
 *		index	index of the timing category's data in the
 *			timer array
 *		label	NULL-terminated character string to be used
 *			to identify this category of timing data in
 *			output files.
 *
 *------------------------------------------------------------------------*/
static void TimerRegister(Home_t *home, int index, const char *label) 
{
	if (index >= TIMER_BLOCK_SIZE) Fatal("TimerRegister: invalid timer");

	home->timers[index].startTime = 0.0;
	home->timers[index].accum     = 0.0;
	home->timers[index].incr      = 0.0;
	home->timers[index].started   = 0;
	home->timers[index].name      = (char *) label;

	return;
}


/*-------------------------------------------------------------------------
 *
 *	Function:	TimerInit
 *	Description:    Allocates storage for the various timing values
 *			to be generated, initializes the values for each
 *			timing category and registers the category with
 *			an identifying string.  The identifying string
 *			is displayed with the associated timing data in
 *			the tmNNNN and time_N files.
 *
 *------------------------------------------------------------------------*/
void TimerInit(Home_t *home)
{
	home->timers = (Timer_t *) malloc(TIMER_BLOCK_SIZE * sizeof(Timer_t));

	TimerRegister(home, TOTAL_TIME,            "total time           ");
	TimerRegister(home, INITIALIZE,            "initialization time  ");
	TimerRegister(home, SORT_NATIVE_NODES,     "sort native nodes    ");
	TimerRegister(home, COMM_SEND_GHOSTS,      "comm send ghosts     ");
	TimerRegister(home, GHOST_COMM_BARRIER,    "sync after ghost comm");
	TimerRegister(home, CELL_CHARGE,           "cell charge          ");
	TimerRegister(home, CELL_CHARGE_BARRIER,   "sync after cell charge");
	TimerRegister(home, CALC_FORCE,            "nodal force          ");
	TimerRegister(home, LOCAL_FORCE,           "  local seg force    ");
	TimerRegister(home, REMOTE_FORCE,          "  remote seg force   ");
	TimerRegister(home, CALC_FORCE_BARRIER,    "sync after force calc");
	TimerRegister(home, CALC_VELOCITY,         "nodal velocity       ");
	TimerRegister(home, CALC_VELOCITY_BARRIER, "sync after vel calc  ");
	TimerRegister(home, COMM_SEND_VELOCITY,    "comm send velocity   ");
	TimerRegister(home, SPLIT_MULTI_NODES,     "split multi-nodes    ");
	TimerRegister(home, COLLISION_HANDLING,    "handle collisions    ");
	TimerRegister(home, POST_COLLISION_BARRIER,"barrier after col    ");
	TimerRegister(home, COL_SEND_REMESH,       "send col topo changes");
	TimerRegister(home, COL_FIX_REMESH,        "apply col topo changes");
	TimerRegister(home, COL_FORCE_UPDATE,      "collision force adj. ");
	TimerRegister(home, GENERATE_IO,           "generate all output  ");
	TimerRegister(home, PLOT,                  "  X window plot      ");
	TimerRegister(home, IO_BARRIER,            "sync after I/O       ");
	TimerRegister(home, REMESH_START_BARRIER,  "sync before remesh   ");
	TimerRegister(home, REMESH,                "remesh dislocations  ");
	TimerRegister(home, SEND_REMESH,           "send remesh changes  ");
	TimerRegister(home, FIX_REMESH,            "apply remesh changes ");
	TimerRegister(home, FORCE_UPDATE_REMESH,   "remesh force adj.    ");
	TimerRegister(home, REMESH_END_BARRIER,    "sync after remesh    ");
	TimerRegister(home, MIGRATION,             "migration            ");
	TimerRegister(home, MIGRATION_BARRIER,     "sync after migration ");
	TimerRegister(home, LOADCURVE,             "load curve           ");
	TimerRegister(home, LOAD_BALANCE,          "load balance         ");
	TimerRegister(home, SEGFORCE_COMM,         "segment force comm   ");

	return;
}


/*-------------------------------------------------------------------------
 *
 *	Function:	TimerStart
 *	Description:    (Re)starts the timer for the specified timing category.
 *	Arguments:
 *		index	index of the timing category's data in the
 *			timer array
 *
 *------------------------------------------------------------------------*/
void TimerStart(Home_t *home, int index)
{
/*
 *	If the timing category is invalid, or the timer for this category
 *	has already been started, there's a bug somewhere, so just abort
 *	and let the user fix the code.
 */
	if (index >= TIMER_BLOCK_SIZE)
		Fatal("TimerStart: invalid timer");

	if (home->timers[index].started != 0)
		Fatal("TimerStart: timer already started");

#ifdef PARALLEL
	home->timers[index].startTime = MPI_Wtime();
#else
	home->timers[index].startTime = (real8)clock() / CLOCKS_PER_SEC;
#endif
	home->timers[index].started = 1;

#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_BY_NAME_BEGIN(home->timers[index].name, SCOREP_USER_REGION_TYPE_COMMON );
#endif

#ifdef USE_CALIPER
    static cali::Annotation ann("paradis.timer", CALI_ATTR_SCOPE_PROCESS | CALI_ATTR_NESTED);
    ann.begin(home->timers[index].name);
#endif

    return;
}


/*-------------------------------------------------------------------------
 *
 *	Function:	TimerStop
 *	Description:    Stops the timer for the specified timing
 *			category.
 *	Arguments:
 *		index	index of the timing category's data in the
 *			timer array
 *
 *------------------------------------------------------------------------*/
void TimerStop(Home_t *home, int index)
{

#ifdef USE_CALIPER
    static cali::Annotation ann_index("paradis.timer.index", 
                                       CALI_ATTR_NOMERGE     | 
                                       CALI_ATTR_SKIP_EVENTS |
                                       CALI_ATTR_SCOPE_PROCESS);

    // allow us to keep ParaDiS' ordering of timers in Caliper
    ann_index.set(index);

    static cali::Annotation ann("paradis.timer");
    ann.end();
#endif

/*
 *	If the timing category is invalid, or the timer for this category
 *	was not already started, there's a bug somewhere, so just abort
 *	and let the user fix the code.
 */
	if (index >= TIMER_BLOCK_SIZE)
		Fatal("TimerStop: invalid timer");

	if (home->timers[index].started != 1)
		Fatal("TimerStop: timer not running");

/*
 *	Stop the timer for this category, calculate the current timing value,
 *	add it to the accumulated total and stop the timer.
 */
#ifdef PARALLEL
	real8 now = MPI_Wtime();
#else
	real8 now = (real8)clock() / CLOCKS_PER_SEC;
#endif
	home->timers[index].incr += now - home->timers[index].startTime;
	home->timers[index].accum += now - home->timers[index].startTime;
	home->timers[index].started = 0;

#ifdef SCOREP_USER_ENABLE
    SCOREP_USER_REGION_BY_NAME_END( home->timers[index].name );
#endif

	return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:       TimerInitDLBReset
 *      Description:    Resets all timers except total and initialization
 *                      timers to zero.  The only real reason this exists
 *                      is so we don't count the time spent in initial
 *                      load-balancing steps to any timers but the
 *                      initializtaion timer and of course total overall
 *                      time.
 *
 *------------------------------------------------------------------------*/
void TimerInitDLBReset(Home_t *home)
{
        int index;

        for (index = 0; index < TIMER_BLOCK_SIZE; index++) {

                if ((index == TOTAL_TIME) ||
                    (index == INITIALIZE)) {
                    continue;
                }

                home->timers[index].incr  = 0.0;
                home->timers[index].accum = 0.0;
                home->timers[index].save  = 0.0;
        }

        return;
}

/*-------------------------------------------------------------------------
 *
 *	Function:	TimerClearAll
 *	Description:    Resets the latest increment for all timer
 *			categories except TOTAL_TIME to zero.
 *
 *------------------------------------------------------------------------*/
void TimerClearAll(Home_t *home)
{
	int index;

	for (index = 0; index < TIMER_BLOCK_SIZE; index++) {
		if (index == TOTAL_TIME) continue;
		home->timers[index].incr = 0.0;
	}
	return;
}

/*-------------------------------------------------------------------------
 *
 *	Function:	TimerReinitialize
 *	Description:    Reinitialize all timers to zero.
 *
 *------------------------------------------------------------------------*/
void TimerReinitialize(Home_t *home)
{
	int index;

	for (index = 0; index < TIMER_BLOCK_SIZE; index++) {

                home->timers[index].startTime = 0.0;
                home->timers[index].accum     = 0.0;
                home->timers[index].incr      = 0.0;

/*
 *              If the timer is currently running, reset the start
 *              time to the current time.
 */
                if (home->timers[index].started) {
#ifdef PARALLEL
                    home->timers[index].startTime = MPI_Wtime();
#else
                    home->timers[index].startTime = (real8)clock() /
                                                    CLOCKS_PER_SEC;
#endif
                }
	}

	return;
}

/*-------------------------------------------------------------------------
 *
 *	Function:	TimerClear
 *	Description:    Resets the latest increment for the specified
 *			timer category to zero.
 *	Arguments:
 *		index	index of the timing category's data in the
 *			timer array
 *
 *------------------------------------------------------------------------*/
void TimerClear(Home_t *home, int index)
{
	home->timers[index].incr = 0.0;

	return;
}


/*-------------------------------------------------------------------------
 *
 *	Function:	TimerSave
 *	Description:    Saves the latest increment for the specified
 *			timer category. (Currently only needed to preserve
 *			NodeForce timings for use in later load-balancing
 *			calculations.
 *	Arguments:
 *		index	index of the timing category's data in the
 *			timer array
 *
 *------------------------------------------------------------------------*/
void TimerSave(Home_t *home, int index)
{
	home->timers[index].save = home->timers[index].incr;

	return;
}


/*-------------------------------------------------------------------------
 *
 *	Function:	TimerPrint
 *	Description:    Creates task-specific timing data files at job
 *			termination.  This can play havoc with file systems
 *			when the job spans thousands of tasks, so this
 *			capability by default is disabled.  It will only
 *			come into play if WRITE_TASK_TIMES is defined
 *			to the preprocessor.
 *
 *------------------------------------------------------------------------*/
void TimerPrint(Home_t *home) 
{
#ifdef WRITE_TASK_TIMES
/*
 *      If user has disabled saving of timer data via the control file
 *      just return without doing anything.
 */
        if (home->param->savetimers != 1) { return; }

/*
 * Task-specific output files will be named time_NNN where
 * <NNN> is the task id.
 */
    char path[128];

    sprintf(path, "%s/time_%d", DIR_TIMERS, home->myDomain);

    FILE *file = fopen(path,"w");

    if (file)
    {
        fprintf(file, "domain : %d     total cycles : %d\n\n", home->myDomain, home->cycle);
   
        for (int i=0; (i<TIMER_BLOCK_SIZE); i++) 
            fprintf(file, "%s\t%12.6f\n", home->timers[i].name, home->timers[i].accum);
   
        fclose(file);
    }

#endif
    return;
}


/*-------------------------------------------------------------------------
 *
 *	Function:	TimerAtRestart
 *	Description:    Collects timing data for every timing category
 *			from all processes, calculates min/max/avg values
 *			for each category, and creates a single timing
 *			data file with the information.  File name will
 *			be tmNNNN where <NNNN> is the savetimerscounter.
 *			Unless otherwise specified, this value will be the
 *			same as the savecncounter used in naming the
 *			rsNNNN restart files.
 *
 *------------------------------------------------------------------------*/
void TimeAtRestart(Home_t *home, int stage) 
{
	Param_t *param   = home->param;
	int      count   = TIMER_BLOCK_SIZE;
	int      numDoms = home->numDomains;
/*
 *	Collect timing values from all processors
 */
	real8 *localTimes = (real8 *) malloc(count * sizeof(real8));
    real8 *allTimes   = 0;

	for (int i=0; (i<count); i++) 
		localTimes[i] = home->timers[i].accum;

	if (home->myDomain == 0)
		allTimes = (real8 *) malloc(count * numDoms * sizeof(real8));

#ifdef PARALLEL
	MPI_Gather(localTimes, count, MPI_DOUBLE, allTimes, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
	for (int i=0; i<count; i++) 
       allTimes[i] = localTimes[i];
#endif

/*
 *	Have processor zero aggregate the values for all timing categories
 *	from all other processors and compute minimum, maximum, and average
 *	times over all processors.
 */
	if (home->myDomain == 0) 
    {
		real8 *minTimes = (real8 *) malloc(count * sizeof(real8));
		real8 *maxTimes = (real8 *) malloc(count * sizeof(real8));
		real8 *avgTimes = (real8 *) malloc(count * sizeof(real8));

		for (int i=0 ; (i<count); i++) 
        {
			minTimes[i] = 99.9e99;
			maxTimes[i] = 0.0;
			avgTimes[i] = 0.0;
		}

		for (int i=0; (i<count); i++) 
        {
			for (int j=0; (j<numDoms); j++) 
            {
				int k = (j*count) + i;

				if (allTimes[k] < minTimes[i]) { minTimes[i] = allTimes[k]; }
				if (allTimes[k] > maxTimes[i]) { maxTimes[i] = allTimes[k]; }

				avgTimes[i] += allTimes[k];
			}
			avgTimes[i] /= (real8) numDoms;
		}

/*
 *		construct the name for the output timing file for this
 *		restart, write the cycle number, cell and domain
 *		geometry, then the min/max/avg timing values for each
 *		timing category.
 */
	    char path[128];

		if (stage == STAGE_TERM)
			sprintf(path, "%s/tm.final", DIR_TIMERS);
		else
			sprintf(path, "%s/tm%04d", DIR_TIMERS, home->param->savetimerscounter);
   
		FILE *file = fopen(path, "w");

        if (file)
        {
    		fprintf(file, "cycle number :    %d\n", home->cycle);
    		fprintf(file, "domain geometry : %d x %d x %d\n"  , param->nXdoms , param->nYdoms , param->nZdoms );
    		fprintf(file, "cell geometry :   %d x %d x %d\n\n", param->nXcells, param->nYcells, param->nZcells);
    
    		fprintf(file, "                             min         max        avg     avg/max\n\n");
       
    		for (int i=0; (i<count); i++) 
            {
	            real8 balance = 0.0;

    			if (maxTimes[i] > 0.0)
    				balance = avgTimes[i] / maxTimes[i];
    			else
    				balance = 1.0;
    
    			fprintf(file, "%s\t%9.3f   %9.3f   %9.3f   %6.3f\n", 
    				home->timers[i].name, minTimes[i], maxTimes[i], avgTimes[i], balance);
    		}
    		fclose(file);
        }

/*
 *		Free up temporary storage allocated on domain zero only
 */
        if (minTimes) { free(minTimes); minTimes=0; }
        if (maxTimes) { free(maxTimes); maxTimes=0; }
        if (avgTimes) { free(avgTimes); avgTimes=0; }

	}   /* end if (home->myDomain == 0) */

	if (localTimes) { free(localTimes); localTimes=0; }
    if (allTimes  ) { free(allTimes  ); allTimes  =0; }

	return;
}
