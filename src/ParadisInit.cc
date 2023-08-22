/**************************************************************************
 *
 *      Mopdule:     ParadisInit.c
 *      Description: Contains functions for doing one-time initializations
 *                   before entering the main loop of the application.
 *
 *      Includes functions:
 *          ParadisInit()
 *
 *************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Init.h"
#include "Util.h"
#include "MPI_Utils.h"
#include "SharedMem.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

// ParadisInit_Summary()
//
// Will print a condensed summary of the simulation control parameters
// to the log file.  In some respects - just a duplication of the 
// information encapsulated in the control file, but also reflects
// command line arguments that override parameters in the control file.
//-------------------------------------------------------------------------

static void ParadisInit_Summary
(
   const Home_t *home
)
{
    if (home && home->param && (home->myDomain==0) )
    {
        const Home_t  *h = home;        // pointe to the Home_t  structure
        const Param_t *p = home->param; // points to the Param_t structure
              real8   *v = 0;           // (used for printing various matrices and arrays)

        printf("+----------------------------------------------------------------------------------------\n" ); 
        printf("| ParaDiS Simulation and control initialization summary\n" );
        printf("|\n" );
        printf("| MPI and Cell configuration controls\n" ); 
        printf("|\n" );
        printf("|   MPI ranks (total)            : %d\n"                  , h->numDomains  );
        printf("|   MPI geometry                 : %d x %d x %d (xyz)\n"  , p->nXdoms    , p->nYdoms  , p->nZdoms  );
        printf("|   cell geometry                : %d x %d x %d (xyz)\n"  , p->nXcells   , p->nYcells , p->nZcells );
        printf("|\n" );
        printf("| Simulation box and periodic boundary conditions\n" );
        printf("|   simulation node count        : %d\n"       , p->nodeCount  );
        printf("|   simulation box (x)           : [ %+0.0lf,%+0.0lf ] (%0.0lfb, %0.1lfum)\n" , p->minSideX, p->maxSideX, p->Lx, (p->Lx * p->burgMag * 1e6) );
        printf("|   simulation box (y)           : [ %+0.0lf,%+0.0lf ] (%0.0lfb, %0.1lfum)\n" , p->minSideY, p->maxSideY, p->Ly, (p->Ly * p->burgMag * 1e6) );
        printf("|   simulation box (z)           : [ %+0.0lf,%+0.0lf ] (%0.0lfb, %0.1lfum)\n" , p->minSideZ, p->maxSideZ, p->Lz, (p->Lz * p->burgMag * 1e6) );
        printf("|   burgmag                      : %0.4le\n"   , p->burgMag    );
        printf("|\n" );
        printf("|   periodic boundary (x)        : %s\n"       , (p->xBoundType==Periodic) ? "periodic" : ( (p->xBoundType==Free) ? "free" : "reflecting" ) );
        printf("|   periodic boundary (y)        : %s\n"       , (p->yBoundType==Periodic) ? "periodic" : ( (p->yBoundType==Free) ? "free" : "reflecting" ) );
        printf("|   periodic boundary (z)        : %s\n"       , (p->zBoundType==Periodic) ? "periodic" : ( (p->zBoundType==Free) ? "free" : "reflecting" ) );

        printf("|\n" );
        printf("|   decompType                   : %d\n"       , p->decompType     );
        printf("|   dlb_freq                     : %d\n"       , p->DLBfreq        );
        printf("|   dlb_cycles                   : %d\n"       , p->numDLBCycles   );
        printf("|\n" );
        printf("| Integration and timestep controls\n" );
        printf("|   integrator                   : %s\n"       , p->timestepIntegrator     );

        if ( StrEquiv(p->timestepIntegrator, "trapezoid") )
        {
        printf("|   trapezoid max integrations   : %d\n"       , p->trapezoidMaxIterations );
        printf("|\n" );
        }

        if ( StrEquiv(p->timestepIntegrator, "trapezoid-kinsol") )
        {
        printf("|   KINSOL_UseNewton             : %d\n"       , p->KINSOL_UseNewton         );
        printf("|   KINSOL_MaxIterations         : %d\n"       , p->KINSOL_MaxIterations     );
        printf("|   KINSOL_MaxLinIterations      : %d\n"       , p->KINSOL_MaxLinIterations  );
        printf("|   KINSOL_NumPriorResiduals     : %d\n"       , p->KINSOL_NumPriorResiduals );
        printf("|   KINSOL_PrintLevel            : %d\n"       , p->KINSOL_PrintLevel        );
        printf("|   KINSOL_LogFile               : %s\n"       , p->KINSOL_LogFile           );
        printf("|   KINSOL_EtaVal                : %+lf\n"     , p->KINSOL_EtaVal            );
        printf("|\n" );
        }

        if ( StrEquiv(p->timestepIntegrator, "arkode") )
        {
        printf("|   ARKODE_IRKtable              : %d\n"       , p->ARKODE_IRKtable          );
        printf("|   ARKODE_FPsolver              : %d\n"       , p->ARKODE_FPsolver          );
        printf("|   ARKODE_MaxNonLinIters        : %d\n"       , p->ARKODE_MaxNonLinIters    );
        printf("|   ARKODE_MaxLinIters           : %d\n"       , p->ARKODE_MaxLinIters       );
        printf("|   ARKODE_FPaccel               : %d\n"       , p->ARKODE_FPaccel           );
        printf("|   ARKODE_PredictorMethod       : %d\n"       , p->ARKODE_PredictorMethod   );
        printf("|   ARKODE_AdaptivityMethod      : %d\n"       , p->ARKODE_AdaptivityMethod  );
        printf("|   ARKODE_NonLinCoef            : %+lf\n"     , p->ARKODE_NonLinCoef        );
        printf("|   ARKODE_LinCoef               : %+lf\n"     , p->ARKODE_LinCoef           );
        printf("|\n" );
        }

        printf("|   start cycle                  : %d\n"       , p->cycleStart        );
        printf("|   last  cycle                  : %d\n"       , p->maxstep           );
        printf("|\n" );
        printf("|   timeNow                      : %le\n"      , p->timeNow           );
        printf("|   timeStart                    : %le\n"      , p->timeStart         );
        printf("|   timeEnd                      : %le\n"      , p->timeEnd           );
        printf("|\n" );
        printf("|   deltaTT                      : %0.2le\n"   , p->deltaTT           );
        printf("|   maxDT                        : %0.2le\n"   , p->maxDT             );
        printf("|   dtIncrementFact              : %0.1lf\n"   , p->dtIncrementFact   );
        printf("|   dtDecrementFact              : %0.1lf\n"   , p->dtDecrementFact   );
        printf("|   dtExponent                   : %0.1lf\n"   , p->dtExponent        );
        printf("|\n" );
        printf("|   minseg                       : %0.1lf\n"   , p->minSeg            );
        printf("|   maxseg                       : %0.1lf\n"   , p->maxSeg            );
        printf("|   rtol                         : %0.1lf\n"   , p->rTol              );
        printf("|   rmax                         : %0.1lf\n"   , p->rmax              );
        printf("|   rc                           : %0.4lf\n"   , p->rc                );
        printf("|   rann                         : %0.1lf\n"   , p->rann              );
        printf("|   remesh area (min)            : %0.1lf\n"   , p->remeshAreaMin     );
        printf("|   remesh area (max)            : %0.1lf\n"   , p->remeshAreaMax     );
        printf("|   split multinode freq         : %d\n"       , p->splitMultiNodeFreq  );
        printf("|   split multinode alpha        : %0.3lf\n"   , p->splitMultiNodeAlpha );
        printf("|\n" );

        printf("| FMM control parameters\n" );
        printf("|   FMM enabled                  : %s\n"       , p->fmEnabled ? "yes" : "no"   );

        if (p->fmEnabled)
        {
        printf("|   FMM layers                   : %d\n"       , p->fmNumLayers        );
        printf("|   FMM mp order                 : %d\n"       , p->fmMPOrder          );
        printf("|   FMM expansion                : %d\n"       , p->fmExpansionOrder   );
        printf("|   FMM npoints                  : %d\n"       , p->fmNumPoints        );
        }
        printf("|\n" );
        printf("| Mobility\n" );
        printf("|   mobility law                 : %s\n"       , p->mobilityLaw        );
        printf("|   material type                : %s\n"       , p->materialTypeName   );
        printf("|\n" );
        printf("|   mobility (screw)             : %0.1le\n"   , p->MobScrew           );
        printf("|   mobility (edge )             : %0.1le\n"   , p->MobEdge            );
        printf("|   mobility (climb)             : %0.1le\n"   , p->MobClimb           );
        printf("|   mobility (glide)             : %0.1le\n"   , p->MobGlide           );
        printf("|   mobility (line )             : %0.1le\n"   , p->MobLine            );

        printf("|\n" );
        printf("|   mobility climb junc          : %0.1le\n"   , p->MobClimbJunc       );
        printf("|   mobility split               : %d\n"       , p->MobSplit           );
        printf("|   mobility climb split         : %0.1le\n"   , p->MobClimbSplit      );
        printf("|\n" );
        printf("|   shear velocity               : %0.1lf\n"   , p->shearVelocity      );
        printf("|\n" );
        printf("|   enforce glide planes         : %s\n"       , p->enforceGlidePlanes ? "yes" : "no" );
        printf("|   enable  cross slip           : %s\n"       , p->enableCrossSlip    ? "yes" : "no" );
        printf("|\n" );

 
        v = (real8 *) p->elasticConstantMatrix;
        printf("|   elasticConstantMatrix        : [ %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le  \n", v[0], v[1], v[2], v[3], v[4], v[5] ); v+=6;
        printf("|                                    %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le  \n", v[0], v[1], v[2], v[3], v[4], v[5] ); v+=6; 
        printf("|                                    %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le  \n", v[0], v[1], v[2], v[3], v[4], v[5] ); v+=6; 
        printf("|                                    %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le  \n", v[0], v[1], v[2], v[3], v[4], v[5] ); v+=6; 
        printf("|                                    %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le  \n", v[0], v[1], v[2], v[3], v[4], v[5] ); v+=6; 
        printf("|                                    %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le %+0.2le ]\n", v[0], v[1], v[2], v[3], v[4], v[5] ); v+=6; 

        printf("|\n" );
        printf("| Simulation parameters\n" );
        printf("|   shear modulus                : %le\n"       , p->shearModulus      );
        printf("|   Young's modulus              : %le\n"       , p->YoungsModulus     );
        printf("|   poisson ratio                : %lf\n"       , p->pois              );
        printf("|   ecrit                        : %lf\n"       , p->ecrit             );
        printf("|   c/a                          : %lf\n"       , p->cOVERa            );
        printf("|   temp                         : %0.2lf K\n"  , p->TempK             );
        printf("|   pressure                     : %0.4lf MPa\n", p->pressure / 1e6    );
        printf("|   load type                    : %d\n"        , p->loadType          );

        v = (real8 *) p->appliedStress;
        printf("|   applied stress               : [ %+0.1le %+0.1le %+0.1le %+0.1le %+0.1le %+0.1le ]\n" ,  v[0], v[1], v[2], v[3], v[4], v[5] );
        printf("|   eRate                        : %+lf\n"      , p->eRate             );

        v = (real8 *) p->edotdir;
        printf("|   edotdir                      : [ %+0.4lf %+0.4lf %+0.4lf ]\n" ,  v[0], v[1], v[2] );

        v = (real8 *) p->rotMatrix;;
        printf("|   rotMatrix                    : [ %+0.4lf %+0.4lf %+0.4lf  \n" , v[0], v[1], v[2] ); v+=3;
        printf("|                                    %+0.4lf %+0.4lf %+0.4lf  \n" , v[0], v[1], v[2] ); v+=3;
        printf("|                                    %+0.4lf %+0.4lf %+0.4lf ]\n" , v[0], v[1], v[2] ); v+=3;

        printf("|\n" );

#if defined ANISOTROPIC
        printf("| Anisotropic control\n" );
        printf("|   harmonics terms base         : %d\n"    , p->anisoHarmonicsNumTermsBase );
        printf("|   Ntheta                       : %d\n"    , p->anisoNumThetaPoints        );
        printf("|   Nphi                         : %d\n"    , p->anisoNumPhiPoints          );
        printf("|   C11                          : %lf\n"   , p->C11                        );
        printf("|   C12                          : %lf\n"   , p->C12                        );
        printf("|   C13                          : %lf\n"   , p->C13                        );
        printf("|   C33                          : %lf\n"   , p->C33                        );
        printf("|   C44                          : %lf\n"   , p->C44                        );
        printf("|\n" );
#endif
        printf("| GPU and OpenMP acceleration\n" );
        printf("|   gpu_enabled                  : %s\n"    , p->gpu_enabled ? "yes" : "no" );
 
        if (p->gpu_enabled)
        {
        printf("|   gpu_nstrm                    : %d\n"    , p->gpu_nstrm    );
        printf("|   gpu_nthrd                    : %d\n"    , p->gpu_nthrd    );
        printf("|   gpu_npblk                    : %d\n"    , p->gpu_npblk    );
        }
        printf("|\n" );
        printf("+----------------------------------------------------------------------------------------\n" ); 
    }
}

/*-------------------------------------------------------------------------
 *
 *      Function:    ParadisInit
 *      Description: Create the 'home' structure, setup timing categories
 *                   and initialize MPI.
 *
 *------------------------------------------------------------------------*/
void ParadisInit
(
   Home_t **homeptr ,   ///< address of pointer to Home data structure (allocate/initialized here)
   int      argc    ,   ///< application argument count
   char    *argv[]      ///< application argument vector
)
{

        Home_t *home    = InitHome();
               *homeptr = home;
    
        TimerInit(home);
    
#ifdef PARALLEL
/*
 *      If ParaDiS has been compiled with OpenMP, initialize the MPI
 *      library for threads, otherwise do standard MPI initialization
 *
 *      NOTE: The standard ParaDiS build only requires MPI thread support
 *            at the MPI_THREAD_FUNNELED level, but when compiled with
 *            the hooks for communicating with the ARL FEM code (fed3)
 *            the required MPI thread support level is MPI_THREAD_MULTIPLE.
 */
#ifdef _OPENMP

#ifdef _ARLFEM
        int requiredLevel = MPI_THREAD_MULTIPLE;
#else
        int requiredLevel = MPI_THREAD_FUNNELED;
#endif  /* ifdef _ARLFEM */

        int providedLevel = MPI_THREAD_SINGLE;                         // (default)
        MPI_Init_thread(&argc, &argv, requiredLevel, &providedLevel);  // (sets actual provided level)

#else  /* _OPENMP not defined */
        MPI_Init(&argc, &argv); 

#endif  /* end ifdef _OPENMP */

#ifdef USE_SCR
/*
 *      Initialize the Scalable Checkpoint/Restart library.  Must be done
 *      *after* MPI_Init() is called.
 */
        SCR_Init();
#endif

        MPI_Comm_rank(MPI_COMM_WORLD, &home->myDomain      );
        MPI_Comm_size(MPI_COMM_WORLD, &home->numDomains    );

        MPI_Comm_rank(MPI_COMM_WORLD, &home->mpi_rank      );
        MPI_Comm_size(MPI_COMM_WORLD, &home->mpi_size      );
        MPI_Node_Rank(home->mpi_node_rank, home->mpi_node_size);
        MPI_Hostname (home->mpi_hostname, sizeof(home->mpi_hostname));
#endif  // ifdef PARALLEL

#ifndef PARALLEL
        // defaults when MPI is not active...

        home->myDomain      = 0;
        home->numDomains    = 1;

        home->mpi_rank      = 0;
        home->mpi_size      = 1;
        home->mpi_node_rank = 0;
        home->mpi_node_size = 1;
        gethostname(home->mpi_hostname,sizeof(home->mpi_hostname));

#endif  // ifndef PARALLEL

#if defined PARALLEL & defined _OPENMP
/*
 *      If MPI and threads are both enabled, make sure the MPI library
 *      provides a sufficient level of threading support.
 */
        if (providedLevel != requiredLevel) 
        {
            char *levelStr=0;

            switch (providedLevel) 
            {
                case (MPI_THREAD_SINGLE    ) : { levelStr = (char *) "MPI_THREAD_SINGLE"    ; break; }
                case (MPI_THREAD_FUNNELED  ) : { levelStr = (char *) "MPI_THREAD_FUNNELED"  ; break; }
                case (MPI_THREAD_SERIALIZED) : { levelStr = (char *) "MPI_THREAD_SERIALIZED"; break; }
                case (MPI_THREAD_MULTIPLE  ) : { levelStr = (char *) "MPI_THREAD_MULTIPLE"  ; break; }
                default                      : { levelStr = (char *) "UNKNOWN"              ; break; }
            }

            if (home->myDomain == 0) 
                Fatal("ParaDiS requires MPI thread support level MPI_THREAD_FUNNELED \n"
                      "but the current MPI version only provides level %s.", levelStr);

        }
#endif  /* end of 'if defined PARALLEL & defined _OPENMP' */

        TimerStart(home,TOTAL_TIME);

        TimerStart(home,INITIALIZE);
        Initialize(home,argc,argv );  
        TimerStop (home,INITIALIZE);
    
/*
 *      For parallel runs with shared memory support enabled, we need
 *      to identify MPI tasks co-located on compute-nodes and set up
 *      MPI communicators for collective communications specific to
 *      those tasks.
 */
        GetNodeTaskGroups(home);

#ifdef _ARLFEM
        FEM_Init(home);
#endif

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        ParadisInit_Summary(home); 

        if (home && (home->myDomain==0) )
           printf("ParadisInit finished\n");

        return;
}
