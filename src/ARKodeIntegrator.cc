/**************************************************************************
 *
 *      Module:      ARKodeIntegrator.cc
 *      Description: Numerical timestep integrator using Runge Kutta
 *                   methods implemented in the ARKStep package within
 *                   SUNDIALS.
 *
 *      Includes exported functions:
 *          CreateARKode()
 *          FreeARKode()
 *
 *      Includes private functions:
 *          InitARKode();
 *          JvEvals();
 *          LinConvFails();
 *          LinIters();
 *          ResizeARKode()
 *          PreserveNodalData()
 *          StatsARKode()
 *
 *      Depends on external functions:
 *          ARKodeIntegrator_RHSFi()
 *
 ***************************************************************************/
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"

#ifdef USE_ARKODE

#include "ParadisSUNDIALS.h"        /* Paradis Vector definition */
#include <arkode/arkode.h>          /* prototypes for ARKODE     */
#include <arkode/arkode_arkstep.h>  /* prototypes for ARKStep    */
#include <arkode/arkode_ls.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include "Timestep.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif


/*
 *  Prototype local functions
 */
static void InitARKode(ParaDiS_Vector v, Home_t *home);
static int  JvEvals(Home_t *home, int update);
static int  LinConvFails(Home_t *home, int update);
static int  LinIters(Home_t *home, int update);
static void ResizeARKode(ParaDiS_Vector v, Home_t *home);


#ifdef DEBUG_TIMESTEP
/*------------------------------------------------------------------------
 *
 *      Function:    PrintLimitingNode
 *      Description: Find the local node (if any) that is responsible
 *                   for the timestep being cut and print out the
 *                   nodal info.
 *                   
 *      Parameters:
 *          IN:  arkode_mem  Pointer to block of memory used internal
 *                           to the ARKStep functions.
 *          IN:  attemptedDT Timestep size requested on entry to the
 *                           timestep inetgrator this cycle.
 *          IN:  takenDT     Timestep size actually achieved this cycle
 *
 *-----------------------------------------------------------------------*/
static void PrintLimitingNode(Home_t *home, void *arkode_mem,
                              real8 attemptedDT, real8 takenDT)
{
        int             i, nodeIndex, status;
        int             indexOfMax;
        long            vecLength;
        real8           localErrMax, globalErrMax;
        real8          *vecData;
        ParaDiS_Vector  errVec;


        localErrMax = 0.0;
        globalErrMax = 0.0;

        errVec = NewParaDiSVector(home);

/*
 *      Get the error for all nodes local to this task from ARKStep,
 *      and determine the local maximum and associated node.
 */
        status = ARKStepGetEstLocalErrors(arkode_mem, errVec);

        vecLength = PV_ACTLENGTH(errVec);
        vecData = PV_DATA(errVec);

        localErrMax = 0.0;
        indexOfMax = 0;

        for (i = 0; i < vecLength; i++) {
            real8  absErr;

            absErr = fabs(vecData[i]);

            if (absErr >= localErrMax) {
                localErrMax = absErr;
                indexOfMax = i;
            }
        }

/*
 *      Determine the largest error among all the tasks.
 */
#ifdef PARALLEL
        MPI_Allreduce(&localErrMax, &globalErrMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

/*
 *      If the local max matches the global max, the associated node
 *      is the one limiting the timestep.
 */
        if (localErrMax == globalErrMax) {
            Node_t *node;

/*
 *          The vector list consists of 3 error components per node, so 
 *          assuming the node list is compacted with no null entries,
 *          the index of the node of interest is easily calculated.
 */
            nodeIndex = indexOfMax / 3;
            node = home->nodeKeys[nodeIndex];

            if (node != (Node_t *)NULL) {
                printf("  Cut timestep: (%d,%d), err %e, attempted %e, "
                       "achieved %e\n", node->myTag.domainID,
                       node->myTag.index, localErrMax, attemptedDT, takenDT);
                PrintNode(node);
            } else {
                printf("  Cut timestep: unknown node, err %e, attempted %e, "
                       "achieved %e\n", localErrMax, attemptedDT, takenDT);
            }
        }

        FreeParaDiSVector(errVec);

        return;
}
#endif  /* ifdef DEBUG_TIMESTEP */


/*------------------------------------------------------------------------
 *
 *      Function:    PreserveNodalData
 *      Description: Both old and new values for certain nodal
 *                   data items are required during timestep
 *                   integration and calculating plastic strain.
 *                   This function copies the values for specified
 *                   items into the appropriate variables.
 *                   
 *-----------------------------------------------------------------------*/
static void PreserveNodalData(Home_t *home, int items)
{
        int i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            Node_t *node;

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

/*
 *          Save previous nodal positions
 */
            if (items & NODE_POSITION) {
                node->oldx = node->x;
                node->oldy = node->y;
                node->oldz = node->z;
            }

/*
 *          Make a copy of the current velocity without wiping out the
 *          copy of the previous velocity.
 */
            if (items & NODE_CURR_VEL) {
                node->currvX = node->vX;
                node->currvY = node->vY;
                node->currvZ = node->vZ;
            }

/*
 *          Make the current velocity the previous velocity (done at the
 *          end of the timestep integrator)
 */
            if (items & NODE_OLD_VEL) {
                node->oldvX = node->currvX;
                node->oldvY = node->currvY;
                node->oldvZ = node->currvZ;
            }
        }

/*
 *      Now do the same as the above for all the ghost nodes.
 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->ghostNodeCount; i++) {
            Node_t *node;

            node = home->ghostNodeList[i];

            if (items & NODE_POSITION) {
                node->oldx = node->x;
                node->oldy = node->y;
                node->oldz = node->z;
            }

            if (items & NODE_CURR_VEL) {
                node->currvX = node->vX;
                node->currvY = node->vY;
                node->currvZ = node->vZ;
            }

            if (items & NODE_OLD_VEL) {
                node->oldvX = node->currvX;
                node->oldvY = node->currvY;
                node->oldvZ = node->currvZ;
            }
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    ARKodeIntegrator
 *
 *      Description: Numerical timestep integrator using Runge Kutta
 *                   methods implemented in the ARKODE package within
 *                   SUNDIALS.
 *
 *                   Note: This function assumes that the nodal
 *                   force/velocity data is accurate for the current
 *                   positions of the nodes on entry to the routine.
 *
 *-----------------------------------------------------------------------*/
void ARKodeIntegrator(Home_t *home)
{
        int      i;                     /* loop counter                     */
        int      ierr;                  /* return code from ARKStep calls   */
        int      mobIterError;          /* mobility iteration error         */
        int      doAll = 1;             /* flag for velocity computation    */
        int      newNodeKeyPtr;
        int      ghostNodeCount;
        real8    reltol, abstol;        /* relative and absolute tolerance  */
        real8    oldDT, newDT;          /* old and new time step size       */
        real8    attemptedDT;           /* time step size to attempt        */
        real8    t0, tout, t;           /* initial, output and final times  */
        void    *arkode_mem=NULL;       /* pointer to ARKStep mem structure */
        Param_t *param;                 /* problem parameters               */
        ParaDiS_Vector solutionVec;     /* position vector for ARKode       */
        ParaDiS_VectorContent scontent; /* position vector content          */

/*
 *      shortcuts to home
 */
        param = home->param;
        newNodeKeyPtr  = home->newNodeKeyPtr;
        ghostNodeCount = home->ghostNodeCount;

/*
 *      If there are no nodes, just return.
 */
        if ((newNodeKeyPtr + ghostNodeCount) <= 0)
        {
            Fatal("ARKodeIntegrator: no nodes!\n");
        }

/*
 *      Create new vector of solution values to send to ARKODE.  
 *      Primary solution is a vector of nodal positions.  Values are
 *      copied from the <home> structure into this solution vector, solutionVec
 */
        solutionVec = NewParaDiSVector(home);

/*
 *      update time step values
 */
        oldDT = param->deltaTT;
        newDT = param->nextDT;

        if (newDT <= 0.0) {
            newDT = param->maxDT;
        }

        param->deltaTT = newDT;
        attemptedDT = newDT;

/*
 *      Set current time
 */
        t0 = param->timeNow; 

/*
 *      Preserve certain items of nodal data for later use.
 */
        PreserveNodalData(home, NODE_POSITION | NODE_CURR_VEL);

/*
 *      General ARKODE setup.  See ARKODE User's Guide for details
 *      on the setup functions.
 */
        if (home->arkodeInitFlag) {
            InitARKode(solutionVec, home);
            home->arkodeInitFlag = 0;
        } else {
            ResizeARKode(solutionVec, home);
        }

        arkode_mem = home->arkodeMemBlock;

#ifdef DEBUG_ARKODE
        FILE *ARKODELog = (FILE *)NULL;

/*
 *      Open files for ARKODE diagnostics output
 */
        if (home->myDomain == 0) {
            if ((ARKODELog = fopen("ARKODE.log", "a")) == NULL) {
                printf("Warning: can't open ARKODE.log file %s.\n",
                       "ARKODE_Log");
            }
            ARKStepSetDiagnostics(arkode_mem, ARKODELog);
        }
#endif

/*
 *      Let ARKStep know the maximum permitted delta time
 */
        ARKStepSetMaxStep(arkode_mem, param->maxDT);

/*
 *      ARKStep timestep, note tout below is not used
 */
        tout = t0 + newDT;
        ierr = ARKStepEvolve(arkode_mem, tout, solutionVec, &t, ARK_ONE_STEP); 

        if (ierr < 0) {
            Fatal("ARKodeIntegrator(): ARKStep Error");
        }

/*
 *      Copy node and ghostNode positions from solutionVec vector
 *      into <home>
 */
        ParaDiSVectorToHome_Positions(solutionVec, home, 1);

/*
 *      Copy the nodal velocities that existed on entry to the timestep
 *      integrator.
 */
        PreserveNodalData(home, NODE_OLD_VEL);  

/*
 *      Compute velocity at new positions
 */
        NodeForce(home, FULL);
        mobIterError = CalcNodeVelocities(home, 0, doAll);
        CommSendVelocity(home);

        if (mobIterError != 0) {
            Fatal("ARKodeIntegrator(): Mobility Error");
        }

/*
 *      ARKStepGetLastStep() returns the timestep size successfully
 *      taken on the last internal step.
 */
        ARKStepGetLastStep(arkode_mem, &oldDT);

/*
 *      ARKStepGetCurrentStep() returns the timestep size to be attempted
 *      on the next internal step.
 */
        ARKStepGetCurrentStep(arkode_mem, &newDT);

#ifdef DEBUG_TIMESTEP
/*
 *      If the timestep size actually taken is smaller than the
 *      timestep size that we originally specified on entry to
 *      this routine, print print out the node that
 *      was *likely* responsible for the time cut.
 */
        if (oldDT < attemptedDT) {
            PrintLimitingNode(home, arkode_mem, attemptedDT, oldDT);
        }
#endif
        param->deltaTT   = oldDT;
        param->realdt    = oldDT;
        param->timeStart = param->timeNow;
        param->nextDT    = newDT;


#ifdef _ARLFEM
 /*
  *     If we're using the FEM code for simulating free surfaces, we
  *     have to handle any nodes/segments that have moved onto or past
  *     the simulation surfaces.
  *
  *     FIX ME!  Do we need to again recalculate forces/velocities for
  *              the nodes that get modified in AdjustNodePosition(), or
  *              can we live with the results until the next cycle???
  */
        FixSurfaceTopology(home);
#endif

/*
 *      Release any remaining temporary resources allocated for the
 *      ARKODE solver and close the log file.
 */
        FreeParaDiSVector(solutionVec);

#ifdef DEBUG_ARKODE
        if (ARKODELog != (FILE *)NULL) {
            fclose(ARKODELog);
        }
#endif
/*
 *      Update linear solver stats
 */
        if (param->ARKODE_FPsolver==0) { 
            LinConvFails(home, 1);
            LinIters(home, 1);
            JvEvals(home, 1);
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    CreateARKode
 *      Description: Signals that ARKode needs to be created.
 *
 *-----------------------------------------------------------------------*/
void CreateARKode(Home_t *home)
{
        home->arkodeInitFlag = 1;

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    FreeARKode
 *      Description: Releases ARKStep memory block
 *
 *-----------------------------------------------------------------------*/
void FreeARKode(Home_t *home)
{
        SUNNonlinSolFree((SUNNonlinearSolver)home->arkodeNonlinearSolver);
        home->arkodeNonlinearSolver = NULL;

        SUNLinSolFree((SUNLinearSolver)home->arkodeLinearSolver);
        home->arkodeLinearSolver = NULL;

        ARKStepFree(&home->arkodeMemBlock);
        home->arkodeMemBlock = (void *)NULL;

        home->arkodeInitFlag = 1;

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    InitARKode
 *      Description: Initializes ARKStep memory and parameters
 *
 *-----------------------------------------------------------------------*/
static void InitARKode(ParaDiS_Vector v, Home_t *home)
{         
        real8   t0;
        real8   reltol, abstol; 
        void    *arkode_mem=NULL;
        Param_t *param;

/*
 *      shortcut to home parameters
 */
        param = home->param;

/*
 *      Create ARKode memory block
 */
        t0 = param->timeNow;
        arkode_mem = ARKStepCreate(NULL, ARKodeIntegrator_RhsFi, t0, v);          
        home->arkodeMemBlock = arkode_mem;  

/*
 *      pass home to ARKStep
 */
        ARKStepSetUserData(arkode_mem, home);

/*
 *      Specify tolerances
 */ 
        reltol = 1e-30;
        abstol = param->rTol;
        ARKStepSStolerances(arkode_mem, reltol, abstol);  

/*
 *      IRK method, min step size, adaptivity method
 */
        // ARKStepSetIRKTableNum(arkode_mem, param->ARKODE_IRKtable);
        ARKStepSetTableNum(arkode_mem, param->ARKODE_IRKtable, -1);
        ARKStepSetMinStep(arkode_mem, 1.0e-20);
        ARKStepSetAdaptivityMethod(arkode_mem, param->ARKODE_AdaptivityMethod,
                                   1, 0, NULL);

/*
 *      If this is a restart set initial step
 */
        if (home->cycle > 0) {
            ARKStepSetInitStep(arkode_mem, param->nextDT);
        }   

/*
 *      Nonlinear solver settings
 */
        ARKStepSetMaxNonlinIters(arkode_mem, param->ARKODE_MaxNonLinIters);
        ARKStepSetPredictorMethod(arkode_mem, param->ARKODE_PredictorMethod);
        ARKStepSetMaxCFailGrowth(arkode_mem, param->dtDecrementFact);
        ARKStepSetNonlinConvCoef(arkode_mem, param->ARKODE_NonLinCoef);

        if (param->ARKODE_FPsolver == 1) { 
/*
 *          Accelerated fixed point solver
 */
            home->arkodeNonlinearSolver = (void *)SUNNonlinSol_FixedPoint(v, param->ARKODE_FPaccel);
            ARKStepSetNonlinearSolver(arkode_mem, (SUNNonlinearSolver)home->arkodeNonlinearSolver);
        } else {
/*
 *          Nonlinear and linear solver setting for Newton method
 */
            home->arkodeNonlinearSolver = (void *)SUNNonlinSol_Newton(v);
            ARKStepSetNonlinearSolver(arkode_mem, (SUNNonlinearSolver)home->arkodeNonlinearSolver);

            home->arkodeLinearSolver = (void *)SUNLinSol_SPGMR(v, PREC_NONE, param->ARKODE_MaxLinIters);
            ARKStepSetLinearSolver(arkode_mem, (SUNLinearSolver)home->arkodeLinearSolver, NULL);
            ARKStepSetEpsLin(arkode_mem, param->ARKODE_LinCoef);
        }

    return;
}


/*------------------------------------------------------------------------
 *
 *      Function:    ResizeARKode
 *      Description: Resizes ARKStep internal vectors
 *
 *-----------------------------------------------------------------------*/
static void ResizeARKode(ParaDiS_Vector v, Home_t *home)
{
        real8   t0;             /* initial, output and final times   */
        real8   reltol, abstol; /* relative and absolute tolerance   */ 
        Param_t *param;         /* problem parameters                */
        void    *arkode_mem;    /* pointer to ARKStep mem structure  */

/*
 *      shortcut to home parameters
 */
        arkode_mem = home->arkodeMemBlock;
        param = home->param;

/*
 *      Resize ARKODE memory
 */
        t0 = param->timeNow; 
        ARKStepResize(arkode_mem, v, 1.0, t0, NULL, home);


        if (param->ARKODE_FPsolver == 1) {
/*
 *      Resize nonlinear solver for accelerated fixed point
 */
            SUNNonlinSolFree((SUNNonlinearSolver)home->arkodeNonlinearSolver);
            home->arkodeNonlinearSolver = (void *)SUNNonlinSol_FixedPoint(v, param->ARKODE_FPaccel);
            ARKStepSetNonlinearSolver(arkode_mem, (SUNNonlinearSolver)home->arkodeNonlinearSolver);
        } else {
/*
 *      Resize nonlinear and linear solver for Newton
 */
            SUNNonlinSolFree((SUNNonlinearSolver)home->arkodeNonlinearSolver);
            home->arkodeNonlinearSolver = (void *)SUNNonlinSol_Newton(v);
            ARKStepSetNonlinearSolver(arkode_mem, (SUNNonlinearSolver)home->arkodeNonlinearSolver);

            SUNLinSolFree((SUNLinearSolver)home->arkodeLinearSolver);
            home->arkodeLinearSolver = (void *)SUNLinSol_SPGMR(v, PREC_NONE, param->ARKODE_MaxLinIters);
            ARKStepSetLinearSolver(arkode_mem, (SUNLinearSolver)home->arkodeLinearSolver, NULL);
            ARKStepSetEpsLin(arkode_mem, param->ARKODE_LinCoef);
        }      
}


/*------------------------------------------------------------------------
 *
 *      Function:    LinConvFails
 *      Description: Since resizing the ARKODE memory requires destroying 
 *                   and re-allocating the linear solver memory structure 
 *                   a workaround is needed to accumulate total values for
 *                   the StatsARKode function. This function updates the 
 *                   total number of linear convergence fails. 
 *
 *-----------------------------------------------------------------------*/
static int LinConvFails(Home_t *home, int update)
{
        long int   lin_fails = 0;
        static int total_lin_convfails = 0;

        if (update) {
            ARKStepGetNumLinConvFails(home->arkodeMemBlock, &lin_fails);
            total_lin_convfails += lin_fails;
            return(0);
        }

        return(total_lin_convfails);
}


/*------------------------------------------------------------------------
 *
 *      Function:    LinIters
 *      Description: Since resizing the ARKODE memory requires destroying 
 *                   and re-allocating the linear solver memory structure 
 *                   a workaround is needed to accumulate total values for
 *                   the StatsARKode function. This function updates the 
 *                   total number of linear iterations. 
 *
 *-----------------------------------------------------------------------*/
static int LinIters(Home_t *home, int update)
{
        long int   lin_iters = 0;
        static int total_lin_iters = 0;

        if (update) {
            ARKStepGetNumLinIters(home->arkodeMemBlock, &lin_iters);
            total_lin_iters += lin_iters;
            return(0);
        }

        return(total_lin_iters);
}


/*------------------------------------------------------------------------
 *
 *      Function:    JvEvals
 *      Description: Since resizing the ARKODE memory requires destroying 
 *                   and re-allocating the linear solver memory structure 
 *                   a workaround is needed to accumulate total values for
 *                   the StatsARKode function. This function updates the 
 *                   total number of Jacobian vector products.
 *
 *-----------------------------------------------------------------------*/
static int JvEvals(Home_t *home, int update)
{
        long int   jv_evals = 0;
        static int total_Jv_evals = 0;

        if (update) {
            ARKStepGetNumJtimesEvals(home->arkodeMemBlock, &jv_evals);
            total_Jv_evals += jv_evals;
            return(0);
        }

        return(total_Jv_evals);
}


/*------------------------------------------------------------------------
 *
 *      Function:    StatsARKode
 *      Description: Outputs ARKODE statistics to a file
 *
 *-----------------------------------------------------------------------*/
void StatsARKode(Home_t *home)
{
        int flag = 0;                   /* error flag                        */
        long int step_attempts = 0;     /* num step attempts                 */
        long int acc_fails     = 0;     /* accuracy failed time steps        */
        long int nonlin_fails  = 0;     /* nonlinear solve convergence fails */
        long int lin_fails     = 0;     /* linear solve convergence fails    */
        long int nonlin_iter   = 0;     /* nonlinear solve iterations        */
        long int lin_iter      = 0;     /* linear solve iterations           */
        long int fi_evals      = 0;     /* number of implicit rhs evals      */
        long int fe_evals      = 0;     /* number of explicit rhs evals      */
        long int jv_evals      = 0;     /* Jv product evaluations            */
        Param_t *param;                 /* problem parameters                */
        FILE *StatFile = (FILE *)NULL;  /* pointer to stat file              */
        void *arkode_mem;

/*
 *      Output ARKODE stats on task zero only
 */
        arkode_mem = home->arkodeMemBlock;

        if (home->myDomain == 0) {

            if ((StatFile = fopen("ARKODE.stats", "w")) == NULL) {
                printf("Warning: can't open ARKODE.stats file %s.\n",
                       "ARKODE_Stats");
            }

            param = home->param;

/*
 *          output diagnostics
 */
            if (StatFile != NULL) {

                TimerStart(home, GENERATE_IO);

                ARKStepGetNumStepAttempts(arkode_mem, &step_attempts);    
                ARKStepGetNumErrTestFails(arkode_mem, &acc_fails);
                ARKStepGetNumNonlinSolvConvFails(arkode_mem, &nonlin_fails);

                if (param->ARKODE_FPsolver == 0) { 
                    lin_fails = LinConvFails(home, 0);
                }    

                ARKStepGetNumNonlinSolvIters(arkode_mem, &nonlin_iter);

                if (param->ARKODE_FPsolver == 0) {
                    lin_iter = LinIters(home, 0);
                }

                ARKStepGetNumRhsEvals(arkode_mem, &fe_evals, &fi_evals);

                if (param->ARKODE_FPsolver == 0) { 
                    jv_evals = JvEvals(home, 0);
                }

                fprintf(StatFile,"Attempts=%-9ld  Steps=%-9ld", 
                        step_attempts, step_attempts-acc_fails-nonlin_fails);
                fprintf(StatFile,"AccFail=%-9ld  NLinFail=%-9ld  "
                                 "LinFail=%-9ld  ",
                        acc_fails, nonlin_fails, lin_fails);
                fprintf(StatFile,"NLinIter=%-9ld  LinIter=%-9ld  ",
                        nonlin_iter, lin_iter);
                fprintf(StatFile,"RhsEval=%-9ld  JvEval=%-9ld  \n",
                        fi_evals, jv_evals);

                TimerStop(home, GENERATE_IO);
            }

/*
 *          Close stat file
 */
            if (StatFile != (FILE *)NULL) {
                fclose(StatFile);
            }

        }  /* end if (home->myDomain == 0) */

        return;
}


#else  /* USE_ARKODE has not been defined */

/*
 *      If we get here, ARKStep was not enabled during compilation so
 *      just abort with a fatal error since we can't use the selected
 *      timestep integrator.
 */
void ARKodeIntegrator(Home_t *home)
{
        Fatal("ARKodeIntegrator(): You must enable ARKODE by\n"
              "setting the <ARKODE_MODE=ON> flag in makefile.setup to\n"
              "use the ARKStep timestep integrator\n");
}
#endif  /* ifdef USE_ARKODE */
