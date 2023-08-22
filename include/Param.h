#pragma once

#ifndef _PDS_PARAM_H
#define _PDS_PARAM_H

/****************************************************************************
 *  Param.h  Define the Param struct, which holds or points to all control
 *           parameters for this simulation
 ***************************************************************************/

#include "Parse.h"
#include "Home.h"
#include "Mobility.h"

/*
 *      Define a couple strings related to the nodal data files
 *      written along with the control file containing the global
 *      parameter values.
 */
#define HDF_DATA_FILE_SUFFIX   ".hdf"
#define NODEDATA_FILE_SUFFIX   ".data"

/*
 *      Version number incremented to 5 when constraint numbers
 *      changed to allow multiple constraints on a single node.
 */
#define NODEDATA_FILE_VERSION  5

#define ANGLE_MOB_CNT 8

struct _param {
/*
 *      Defines the number of domains in each dimension of the
 *      entire problem space.
 */
        int nXdoms;
        int nYdoms;
        int nZdoms;
#if defined _BGP || defined _BGQ
        int taskMappingMode;  /* Only valid for BG/P systems on which a  */
                              /* 3D hardware partition is allocated to   */
                              /* the users job.                          */
                              /*                                         */
                              /* Determines action taken when the user-  */
                              /* supplied domain decomposition is not    */
                              /* consistent with the 3D geometry of the  */
                              /* underlying hardware partition.  Allowed */
                              /* values are:                             */
                              /*     0 -- domain decomposition will be   */
                              /*          reset to be consistent with    */
                              /*          the hardware partition (if     */
                              /*          possible)                      */
                              /*     1 -- the user-supplied domain       */
                              /*          decomposition will be used     */
                              /*     2 -- the code will abort if there   */
                              /*          is a mismatch.                 */

#ifdef _BGQ
        int optimalTaskGeometry[3]; /* Only vaid for BG/Q systems.  This */
                              /* array will be set to the optimal task   */
                              /* geometry for the hardware partition when*/
                              /* possible.  If the optimal geometry can  */
                              /* not be determined, the components of    */
                              /* this array will be set to -1.  This will*/
                              /* be set internally and cannot be speci-  */
                              /* fied by the user.                       */
#endif
#endif
        int assignNewNodeTags; /* During initialization this controls  */
                               /* whether we try to preserve node tags */
                               /* from the restart file or assign new  */
                               /* tags.  This is not a user-supplied   */
                               /* value.  When the current domain      */
                               /* geometry matches that from the       */
                               /* file, node tags *can* be preserved so*/
                               /* this is set to zero.  Otherwise this */
                               /* value is set to 1 so the nodeKeys    */
                               /* array is not sparsely populated.     */
/*
 *      Defines number of cells in each of the dimensions of the
 *      entire problem space.
 */
        int nXcells;
        int nYcells;
        int nZcells;

/*
 *      "natural" min and max cell indices for this domain. (i.e. indices
 *      before being incremented by 1 to allow for ghost cells)
 */
        int iCellNatMin;
        int iCellNatMax;
        int jCellNatMin;
        int jCellNatMax;
        int kCellNatMin;
        int kCellNatMax;

/*
 *      Specifies type of boundary in each dimension:
 *          0 == periodic
 *          1 == free surface
 */
        BoundType_t xBoundType;
        BoundType_t yBoundType;
        BoundType_t zBoundType;

/*
 *      When free surfaces are used, these define the upper and
 *      lower limits on dislocation coordinates.  (i.e. the planes
 *      defining the free surfaces)
 */
        real8 xBoundMin, xBoundMax;
        real8 yBoundMin, yBoundMax;
        real8 zBoundMin, zBoundMax;

/*
 *      Define min/max coordinate limits for each dimension of
 *      the problem space.  NOTE: These are redundant but kept
 *      until all code references to them have been removed.
 *      Use <minCoordinates> and <maxCoordinates> arrays instead.
 */
        real8 minSideX, maxSideX;
        real8 minSideY, maxSideY;
        real8 minSideZ, maxSideZ;

/*
 *      Domain decomposition and rebalance values
 */
        int   decompType;       /* Selects decomposition type */
        int   DLBfreq;          /* how often to load balance */
        int   numDLBCycles;     /* Number of initial load-balance-only */
                                /* cycles to be executed before main   */
                                /* processing loop is entered.  This   */
                                /* is a command line option value, not */
                                /* a control file parameter            */

/*
 *      Simulation time and timestepping controls
 */
        int   cycleStart;    /* Starting cycle number for the simulation */
        int   maxstep;       /* Cycles to execute before terminating */
        real8 timeStart;     /* Initial simulation time */
        real8 timeNow;       /* current simulation time */
        real8 timeEnd;       /* simulation stop time */

        char  timestepIntegrator[MAX_STRING_LEN];
        int   trapezoidMaxIterations;

        real8 deltaTT;       /* duration of previous timestep */
        real8 realdt;        /* *almost* obsolete */
        real8 nextDT;        /* timestep to attempt the next cycle */
        real8 maxDT;         /* Maximum timestep duration permitted*/
        real8 dtIncrementFact;  /* Maximum delta time increment allowed */
        real8 dtDecrementFact;  /* Factor by which delta time is        */
                                /* multiplied when cutting timestep.    */
                                /* value should be between 0 and 1      */
        real8 dtExponent;       /* Used in calculating variable delta   */
                                /* time adjustments                     */
        int   dtVariableAdjustment;  /* Flag indicating if timestep     */
                                     /* adjustments are of variable     */
                                     /* size or set percentages of the  */
                                     /* current dt.                     */
        real8 rTol;      /* Maximum error allowed in timestep */
        real8 rmax;      /* maximum migration distance per timestep */
                         /* for any node   */

        int   gpu_enabled;          /* Set to 1 to enable GPU force calculations     */
        int   gpu_nstrm  ;          /* number of active CUDA streams                 */
        int   gpu_nthrd  ;          /* number of threads per kernel launch           */
        int   gpu_npblk  ;          /* number of segment pairs per cuda stream block */

/*
 *      Define some parameters that apply only if the KINSOL general-
 *      purpose nonlinear system solver has been enabled for use in
 *      the timestep integrator.
 */
        int   KINSOL_UseNewton;     /* Flag indicating if KINSOl integrator */
                                    /* should use the Newton-Krylov non-    */
                                    /* linear solver.  If zero, Anderson    */
                                    /* accelerated fixed-point solve is used*/
        int   KINSOL_MaxIterations; /* Maximum number of non-linear solver  */
                                    /* iterations the KINSOL based timestep */
                                    /* integrators try  before cutting the  */
                                    /* timestep.                            */
        int   KINSOL_MaxLinIterations; /* Maximum number of linear          */
                                    /* iterations within KINSOL's Newton-   */
                                    /* Krylov solver                        */
        int   KINSOL_NumPriorResiduals; /* Number of prior residual vectors */
                                    /* are used with the Anderson           */
                                    /* acceleration when using KINSOL in    */
                                    /* the timestep integrator              */
        int   KINSOL_PrintLevel;    /* Specifies the level of verbosity */
                                    /* in the output from KINSOL        */
        char  KINSOL_LogFile[MAX_STRING_LEN];  /* Name of the file to   */
                                    /* which the KINSOL solver writes   */
                                    /* it's debug information.          */
        real8 KINSOL_EtaVal;        /* Constant in stopping criteria of */
                                    /* linear solve                     */

/*
 *      Define some parameters that apply only if the ARKODEL multirate ODE
 *      system solver has been enabled for use in the timestep integrator.
 */
        int   ARKODE_IRKtable;       /* Flag for desired IRK method         */
        int   ARKODE_FPsolver;       /* Flag to activate Fixed Point solver */
        int   ARKODE_MaxNonLinIters; /* Max number of nonlinear iterations  */
        int   ARKODE_MaxLinIters;    /* Max number of linear iterations     */
        int   ARKODE_FPaccel;        /* Number of acceleration vectors      */
        int   ARKODE_PredictorMethod;   /* Predictor for nonlinear solves   */
        int   ARKODE_AdaptivityMethod;  /* Method for selecting next dt     */
        real8 ARKODE_NonLinCoef;  /* Safety factor for nonlinear solve      */
        real8 ARKODE_LinCoef;     /* Ratio between nonlinear and linear tol */

/*
 *      Discretization parameters and controls for topological changes
 */
        real8 minSeg;    /* min allowable segment length, before */
                         /* removing a node */
        real8 maxSeg;    /* max allowable segment length, before*/
                         /* adding a node*/
        real8 rann;      /* Annihilation distance (units of b) used in     */
                         /* collision handling, rediscretization, multi-   */
                         /* node splitting, etc.  Defines closest distance */
                         /* before dislocations are considered in contact  */

        int remeshRule;
        int collisionMethod;
        real8 remeshAreaMax; /* This is calculated from remeshAreaMin and  */
                             /* remeshAreaRatio and hence not user-provided*/
        real8 remeshAreaMin; /* This values is based on the minSeg value  */
                             /* and hence not specified by the user.      */
        int splitMultiNodeFreq;  /* Code will attempt to split multi-nodes */
                                 /* every cycle that is a multiple of this */
                                 /* value. */

        real8 splitMultiNodeAlpha;     /* factor used to control the split multi-node */
                                       /* vNoise term */

        int useParallelSplitMultiNode; /*Toggle indicating if the code   */
                                       /* should use the parallel method */
                                       /* of multi-node splitting or not */

        int useUniformJunctionLen; /* Toggle indicating whether to create */
                                   /* junctions of a uniform length in    */
                                   /* SplitMultiNodes() or try to calc    */
                                   /* the optimal junction length.        */

/*
 *      Fast Multipole Method parameters
 */
        int fmEnabled;       /* Set to 1 remote forces are caculated */
                             /* using the Fast Multipole code        */

        int fmNumLayers;     /* Number of layers of cells defined for*/
                             /* the Fast Multipole remote force calcs*/

        int fmMPOrder;       /* Order used for multipole expansions  */

        int fmExpansionOrder;/* Order used for polynomial expansions */

        int fmNumPoints;  /* Number of points along each segment    */
                          /* at which remote stress will be         */
                          /* evaluated from the taylor expansion    */
                          /* coefficients. Automatically calculated */
                          /* not supplied by user.                  */

        real8 fmAdjustLen; /* Local increase in cell size to account for segment */
                           /* owned by the cell, not in the cell.                */

        char fmCorrectionTbl[MAX_STRING_LEN];
        int  fmFullCorrectionTbl; /* Flag indicating if the provided    */
                          /* correction table is a full table or has    */
                          /* been compressed by having cubic symmetries */
                          /* factored out.  1 indicates table is full,  */
                          /* 0 (the default) indicates compression.     */
                          /* This value is not user provided, but       */
                          /* dynamically set when the correction table  */
                          /* is read.                                   */

        int  fmEnableAnisotropy;  /* Flag indicating if the anisotropic      */
                          /* fmm is to be enabled.  Default is zero.         */
                          /* TAYLORFMM on and this flag on DOES NOT WORK YET */
                          /* Additionally, if the general anisotropy is      */
                          /* not enabled at runtime via the ANISOTROPIC     */
                          /* definition, this flag will explicitly be        */
                          /* turned off no matter what the user specifies    */

#if (defined ANISOTROPIC) && (defined TAYLORFMM)
        char fmAnisoDerivTbl[MAX_STRING_LEN]; /* Table containing the  */
                                   /* Greeens function derivatives for */
                                   /* anisotropic Fast Multipole code  */
#endif

#ifdef ESHELBY
        int eshelbyfmEnabled; /* Set to 1, remote forces from Eshelby */
                           /* inclusions are caculated using the      */
                           /* Fast Multipole code                     */

        int eshelbyfmMPOrder; /* Order used for multipole expansions   */
                          /* when calculating coefficients for eshelby */
                          /* inclusions.                               */

        int eshelbyfmTaylorOrder; /* Order used for taylor expansions  */
                          /* when calculating coefficients for eshelby */
                          /* inclusions.                               */
        int RefineSegCollision; /* When particles are small. Collision detection */
                          /* is more difficult. Local remesh of dislocation*/
                          /* segments is done is this option is turned on*/
                          /* Default = 0 */
        real8 Ecore2;     /* core energy in Pa for the second material */
        real8 shearModulus2; /* Elastic constant of the second material*/
        real8 pois2;
        real8 MobEshelbyResist; /* Parameter used in the Eshelby mobility
                                   laws to account for mobility changes
                                   inside particles */
        real8 IncStepEnergy;   /* Energy associated with a step when a
                                  dislocation meets an inclusion.
                                  Unit is J/m^2                        */

#endif

#ifdef CALCENERGY
         char fmECorrTbl[MAX_STRING_LEN];  // Correction table for energy
#endif

#ifdef FRCRIT
         real8 minBisection, maxBisection;
#endif

#ifdef SPECTRAL
         int    FFTgrid;
         int    FFTenabled;
         int    localRotation;
#endif


/*
 *      Names of tables for non-FMM far-field forces
 */
        char Rijmfile[MAX_STRING_LEN];
        char RijmPBCfile[MAX_STRING_LEN];

/*
 *      Loading condition parameters
 */
        real8 TempK;         /* Temperature in deg K */
        real8 pressure;      /* Pressure in units of Pa */
        int   loadType;      /* 0 Creep test */
                             /* 1 Constant strain rate test */
                             /* 2 Displacement-controlled test */
                             /* 3 Load-controlled, load vs. time curve */
                             /* 4 Cyclic loading condition */
        real8 appliedStress[6]; /* External stress in units of Pa  */
                                /* as [sigma11, sigma22, sigma33,  */
                                /* sigma23, sigma31, sigma12] when */
                                /*  <loadType> == 0.               */
        real8 eRate;         /* Strain rate. Used when loadType == 1 */
        int   indxErate;     /* to be compatible with micro3d */
        real8 edotdir[3];    /* Uniaxial loading direction accompanying */
                             /* eRate                                   */

        int   crystalRotation;      /* Allow for crystal rotation due to spin during loading */
        int   cubicModulus;         /* Toggle to compute stress/strain modulus for cubic */
                                    /* materials as a function of the loading direction */
        real8 cubicC11;             /* Cubic elastic constants for stress/strain modulus */
        real8 cubicC12;
        real8 cubicC44;

        real8 cTimeOld;      /* Timestep related to cyclic loading */
        real8 netCyclicStrain; /* Net accumulated strain under cyclic load */
        real8 dCyclicStrain; /* Incremental strain under cyclic load */
        int   numLoadCycle;  /* Number of cyclic cycles */
        real8 eAmp;          /* Strain amplitude used with cyclic loading */

/*
 *      For some simulations, the geometric information (burgers vectors,
 *      normal planes) may be provided in a user-specified laboratory frame.
 *      (See the <useLabFrame> control file parameter).  If this is being
 *      done, we'll need to rotate between the lab frame and the standard
 *      crystal frame during various portions of the simulation, so allow
 *      the user to provide a rotation matrix and compute the inverse once
 *      during initialization for subsequent use.
 *
 *      Rotation matrix rows specify the X, Y and Z axes respectively
 *      for the user-defined laboratory frame rather than the crystallographic
 *      frame.  Default matrix is:
 *          1 0 0
 *          0 1 0
 *          0 0 1
 *
 *      Each row of the provided rotation matrix will be normalized during
 *      initialization.
 */
        int   useLabFrame;   /* 0 if standard crystalographic frame is to */
                             /* be used, 1 if user-supplied laboratory    */
                             /* frame is used                             */

        real8 rotMatrix[3][3];  /* Ignored if <useLabFrame == 0 */

        real8 rotMatrixInverse[3][3]; /* This matrix is calculated dynamically */
                                      /* from rotMatrix and is not provided    */
                                      /* through the control file              */

/*
 *      Material and mobility parameters
 */
#ifdef ESHELBY
        int   enableInclusions;
        char  inclusionFile[MAX_STRING_LEN];
                             /* Name of file defining Eshelby inclusions */
                             /* to be inserted into the simulation. By   */
                             /* default this string is empty indicating  */
                             /* no Eshelby inclusions are to be simulated*/
#endif

        char  ecFile[MAX_STRING_LEN];
                             /* Base name for files containing the material */
                             /* specific elastic constants for various      */
                             /* temperature and pressure combinations.      */
                             /* If provided, the poisson ratio and shear    */
                             /* modulus will be calculated based on the     */
                             /* given temp and pressure parameters.         */

        real8 meltTemp;      /* Material melting temperature in degrees K */
        char  meltTempFile[MAX_STRING_LEN];

        char  MobFileNL[MAX_STRING_LEN];
        char  mobilityLaw[MAX_STRING_LEN];
        int   mobilityType;  /* Integer value corresponding to the */
                             /* specified mobility law.  Redundant */
                             /* info, but easier to use in the code*/

/*
 *     Full elastic constant matrix for the simulation material.  For
 *     isotropic simulations this matrix is automatically recalculated
 *     based on the provided shear modulus and poisson ratio, and as
 *     such does not have to be specified in the control parameter file.
 *
 *     For anisotropic simulations this matrix *may* not have to be provided
 *     if the individual constants (C11, C12, etc) are provided (See User's
 *     Guide for details)
 */
       real8  elasticConstantMatrix[6][6];

/*
 *     Inverse of the elastic constant matrix.  This is dynamically calculated
 *     during initialization and is not provided in the control parameter
 *     file.
 */
       real8  elasticConstantMatrixInv[6][6];

#if defined ANISOTROPIC

/*
 *     The following are specific elastic constants for anisotropic
 *     elasticity.  Not all constants apply to all material types
 *     and these *may* be unnecessary if the <ElasticConstants> array
 *     above was provided.  See User's Guide for details.
 *
 *     For cubic crystals constants are:
 *         C11, C12, C44,
 *     For hexagonal crystals constants are:
 *         C11, C12, C13, C33, C44
 *
 *     Note: Cpr will be caluclated internally for simulations with cubic
 *           crystals, but is not specified via a control file
 *           parameter.
 */
       real8  C11, C12, C13, C33, C44, Cpr;

/*
 *     Define the base from which the number of terms in the spherical
 *     harmonics is determined.  The actual number of terms will be
 *     2 * <anisoHarmonicsNumTermsBase> + 1;
 */
       int   anisoHarmonicsNumTermsBase;

/*
 *     Define the number of points into which the theta and phi angles
 *     are discretized when defining the double integral over theta and
 *     phi angles when calculating the glm expansion coefficients
 *     needed for computing the Fq table for the anisotropic elasticity
 */
       int   anisoNumThetaPoints;
       int   anisoNumPhiPoints;


#endif  /* end if defined ANISOTROPIC */



        char  materialTypeName[64];  /* Specifies the material type (i.e  */
                                     /* "BCC", "FCC", etc) to be used for */
                                     /* the simulation                    */

        int   materialType;  /* Type of crystal structure (i.e. BCC, FCC)   */
                             /* This value is set within the code base on   */
                             /* the selected mobility law and is not user-  */
                             /* supplied                                    */

        real8 vacancyConc;            /* Concentration of vacancies in the */
                                      /* crystal: for calculating osmotic  */
                                      /* force (units ?)                   */

        real8 vacancyConcEquilibrium; /* Thermal equilibrium vacacy concen-*/
                                      /* ration in defect free crystal: for*/
                                      /* calculating osmotic force         */

        real8 shearModulus;
        real8 pois;
        real8 burgMag;
        real8 ecrit;

        real8 cOVERa; /* For HCP simulations only.  See description of in */
                      /* user's guide                                     */

        real8 YoungsModulus;

        real8 rc;     /* core radius in elastic interaction calculation */

        real8 Ecore;  /* core energy (wrt. the choice of rc) in unit of Pa */
        int   MDcore; /* core model fitted to MD simulations */
        char  MDcoreMaterial[MAX_STRING_LEN];

        real8 HCPEcoreA;   /* For HCP, there are 3 Ecore values which each */
        real8 HCPEcoreC;   /* apply to specific types of burgers vectors   */
        real8 HCPEcoreCpA; /* (types <a>, <c> and <c+a>) and are supplied  */
                           /* through these three parameters.  These core  */
                           /* energy prefactors are dependent on the core  */
                           /* radius and are in units of the shear modulus */
                           /* (energy per unit length / b^2 units).        */
                           /* These values override the <Ecore> value      */
                           /* above in Initialize.cc.                      */


        int   HCP_CA_SplitFrequency; /* Indicates the frequency (in cycles) */
                                     /* at which the code should attempt to */
                                     /* split HCP <c+a> burgers vectors     */
                                     /* into a <c> and an <a>.  Range is    */
                                     /* 0 <= HCP_CA_SplitFrequency, where   */
                                     /* turns off the c+a splitting         */

        int   splinterSegFrequency;  /* Indicates the frequency (in cycles) */
                                     /* at which the code should attempt to */
                                     /* splinter certain segments into      */
                                     /* a pair of identical segments        */
                                     /* whose burgers vectors sum to the    */
                                     /* burgers vector of the original seg. */
                                     /* Range is 0 <= SplinterSegFrequency  */
                                     /* where 0 turns off splintering       */

        int   enforceGlidePlanes;    /* Toggle indicating if the simulation */
                                     /* should restrict dislocation motion  */
                                     /* to the associated glide plane.      */

        int   allowFuzzyGlidePlanes; /* Only applicable when the            */
                                     /* <enforceGlidePlanes> toggle is set. */
                                     /* Set internally for specific mobility*/
                                     /* modules. Not read from control file */
        int   enableCrossSlip;       /* Toggle indicating if cross-slip is */
                                     /* turned on.  Only used when a glide-*/
                                     /* restrictive mobility is in use     */

        int   crossSlipBCCIndex;     // if cross slip is enabled and BCC active...
                                     //   0 = 3 glideplanes
                                     //   2 = 6 glideplanes
                                     //   3 = 6 glideplanes +              power dissipation + area growth
                                     //   4 = 6 glideplanes + randomness + power dissipation + area growth (default)

        int   mobilityIndex;  /* Internal use only.  Index into the */
                              /* mobility attributes array of the   */
                              /* selected mobility module.          */

/*
 *      The next two function pointers are set during initialization to
 *      indicate the appropriate mobility function for the simulation and
 *      the corresponding mobility initialization function (the init
 *      function is not required and may be left as a NULL pointer)
 */
        void  (*mobilityInitFunc)(Home_t *home);
        int   (*mobilityFunc)(Home_t *home, Node_t *node, MobArgs_t *mobArgs);


/*
 *      The following mobilities are in units of 1/(Pa*s)
 */
        real8 MobScrew;
        real8 MobEdge;
        real8 MobClimb;
        real8 MobGlide;  /* floor on mobility for glide dislocations */
        real8 MobLine;

        real8 MobClimbJunc; /* Climb mobility for the junctions only */
        real8 MobLineJunc ; /* Line mobility for the junctions only */
        int   MobSplit;
        real8 MobClimbSplit;
        real8 shearVelocity; /* Shear sound velocity in units of m/s */

/*
 *      Angle dependent mobility law
 */
        int    MobTheta_cnt;
        real8  MobThetas[8];
        real8  MobGlide_B[8];
        real8  MobGlide_D[8];
        real8  MobGlide_V0[8];

/*
 *      Some variables needed for the BCC_nl mobility
 */
        real8 MobExpP;
        real8 MobExpQ;
        real8 MobPeierls;
        real8 MobDeltaH0;       /* Kink pair activation enthalpy at 0 stress */
        real8 MobAlpha;

/*
 *       At high stress values, atomistic simulations are performed
 *       to obtain the velocity at different shear stresses. And the
 *       data are fitted using the following functions:
 *
 *           v0 = B * C0 * (1 - C*theta/theta_m + D*P)
 *           V = Ml * Tao/Tao_p - v0
 *
 *       where V is the velocity, Tao is the applied shear stress and
 *       Tao_p is the Peierls stress. The fitting is done at different
 *       pressure (P) and temperatures (T). The fitted results are:
 *
 *       Ml = A *C0 *(1-C*theta/theta_m +D*P)
 *
 *       A, B, C0, C, D are fitted constants, theta_m is the melting
 *       temperature
 */
        real8 MobCoeffA;                /* unitless */
        real8 MobCoeffC;                /* units of meter/sec */
        real8 MobCoeffC0;               /* unitless */
        real8 MobCoeffD;                /* units of GPa */

/*
 *      Parameters used for Vanadium non-linear mobility only
 */
        real8 MobG0;

/*
 *      Some values that are specific to the 'faceted' mobility
 *      functions
 */
        int   MobNumTrenches;  /* Number of angular windows within which */
                               /* drag on dislocations is increased.     */
                               /* Defines the number of valid elements in*/
                               /* the trench related arrays below        */

        real8 MobTrenchAngle[MAX_MOB_TRENCHES];
                               /* Centers (specified in angles) of the */
                               /* corresponding windows (trenches) of  */
                               /* reduced mobility                     */

        real8 MobTrenchWidth[MAX_MOB_TRENCHES];
                               /* Half the angular period (in degrees) of */
                               /* the corresponding window of reduced     */
                               /* mobility                                */

        real8 MobTrench[MAX_MOB_TRENCHES];
                               /* Mobility reduction on dislocations along   */
                               /* within the corresponding window of reduced */
                               /* mobility (trench)                          */

        int   MobFacetedUseCusps; /* Toggle indicating whether to do the     */
                                  /* faceting using cusps rather than teh    */
                                  /* standard trenches                       */

/*
 *      Define some parameters used only with HCP materials.  Note: values are
 *      drag coefficients rather than mobilities.
 *
 *      The following drag coefficients are in units of Pa * s.
 */
        real8 HCP_A_Basal_ScrewDrag;
        real8 HCP_A_Prismatic_ScrewDrag;
        real8 HCP_A_1stPyramidal_ScrewDrag;
        real8 HCP_A_2ndPyramidal_ScrewDrag;
        real8 HCP_C_Prismatic_ScrewDrag;
        real8 HCP_CpA_Prismatic_ScrewDrag;
        real8 HCP_CpA_1stPyramidal_ScrewDrag;
        real8 HCP_CpA_2ndPyramidal_ScrewDrag;
        real8 HCP_Sessile_ScrewDrag;

        real8 HCP_A_Basal_EdgeDrag;
        real8 HCP_A_Prismatic_EdgeDrag;
        real8 HCP_A_1stPyramidal_EdgeDrag;
        real8 HCP_A_2ndPyramidal_EdgeDrag;
        real8 HCP_C_Prismatic_EdgeDrag;
        real8 HCP_CpA_Prismatic_EdgeDrag;
        real8 HCP_CpA_1stPyramidal_EdgeDrag;
        real8 HCP_CpA_2ndPyramidal_EdgeDrag;
        real8 HCP_Sessile_EdgeDrag;

        real8 HCP_LineDrag;

/*
 *      For the "relaxation" mobility laws only
 */
        real8 MobRelaxX;
        real8 MobRelaxY;
        real8 MobRelaxZ;

        int   MobRelaxScaleByLength;

/*
 *      For linear tantalum mobility law
 */
        real8 TaDragEdge110;
        real8 TaDragEdge112;
        real8 TaDragScrew112T;
        real8 TaDragScrew112AT;
        real8 TaDragScrew110;
        real8 TaDragLine;
        real8 TaDragClimbEdge;
        real8 TaDragClimbScrew;

/*
 *      For non-linear tantalum mobility law
 */
        real8 TaDragParam110K;
        real8 TaDragParam110B;
        real8 TaDragParam110alpha;

        real8 TaDragParam112TK;
        real8 TaDragParam112TB;
        real8 TaDragParam112Talpha;

        real8 TaDragParam112ATK;
        real8 TaDragParam112ATB;
        real8 TaDragParam112ATalpha;

/*
 *      Allow for specifying that dislocations with certain types
 *      of burgers vectors or line directions are immobile.
 */
        real8 sessileburgspec[30];
        real8 sessilelinespec[30];

/*
 *      Add some variables needed to include inertial terms into
 *      mobility functions.  If the selected mobility does not
 *      have support for inertial effects, these values will be
 *      quietly ignored.
 */
        int   includeInertia;  /* Toggle: 0 == no, 1 == yes */
        real8 massDensity;     /* Units are kg/m^3 */

/*
 *      Velocity statistics and parameters
 */
        real8 vAverage;        /* average nodal velocity */
        real8 vStDev;          /* St.Dev of nodal velocity */

/*
 *      I/O parameters
 */
        char dirname[MAX_STRING_LEN];

        int  writeBinRestart; /* if set, will write data portion of */
                              /* restart file in a binary format    */

        int  doBinRead;  /* If set, will attempt to read binary format */
                         /* restart file.  This flag is set internally */
                         /* and not specified by the user.             */

        int  numIOGroups;  /* number of groups into which to split tasks */
                           /* when doing parallel I/O                    */

        int  skipIO;    /* if set, all major IO is skipped regardless */
                        /* of what types of output generation would   */
                        /* be enabled by other control parameters     */

        int   armfile, armfilefreq, armfilecounter;
        real8 armfiledt, armfiletime;

        int   writeFlux;
        int   writeFluxEdgeDecomp;
        int   writeFluxFullDecomp;
        int   writeFluxFullDecompTotals;
        int   writeFluxSimpleTotals;
        int   writeFluxFreq;
        int   writeFluxCounter;
        real8 writeFluxDT;
        real8 writeFluxTime;

        int   gnuplot, gnuplotfreq, gnuplotcounter;
        real8 gnuplotdt, gnuplottime;

        int   polefigfile, polefigfreq, polefigcounter;
        real8 polefigdt, polefigtime;

        int   povray,  povrayfreq,  povraycounter;
        real8 povraydt, povraytime;

        int   psfile,  psfilefreq;
        real8 psfiledt, psfiletime;

        int   savecn,  savecnfreq,  savecncounter;
        real8 savecndt, savecntime;

        int   savedt,  savedtfreq;

        int   save_narms,  save_narms_freq;

        int   saveprop, savepropfreq;
        real8 savepropdt, saveproptime;

        int   savetimers, savetimersfreq, savetimerscounter;
        real8 savetimersdt, savetimerstime;

        int   savedensityspec[3];

        int   tecplot, tecplotfreq, tecplotcounter;
        real8 tecplotdt, tecplottime;

        int   velfile, velfilefreq, velfilecounter;
        real8 velfiledt, velfiletime;

        int   writeForce, writeForceFreq, writeForceCounter;
        real8 writeForceDT, writeForceTime;

        int   writeVisit, writeVisitFreq, writeVisitCounter;
        int   writeVisitSegments, writeVisitSegmentsAsText;
        int   writeVisitNodes, writeVisitNodesAsText;
        int   writeVisitBurgID, writeVisitForceVector;
        int   writeVisitVelocityVector;
        real8 writeVisitDT, writeVisitTime;

        char  winDefaultsFile[MAX_STRING_LEN];



/*
 *      Lengths of each side of the problem space box.
 */
        real8 Lx, Ly, Lz;


/*
 *      General stuff
 */
        real8 springConst;

        int  numBurgGroups;     /* Number of groups into which different  */
                                /* burgers vectors are organized in order */
                                /* to track dislocation density by burgers*/
                                /* vector.  This number is dependent on   */
                                /* the type of mobility used.             */

        real8 *partialDisloDensity;  /* Dynamically sized array to hold   */
                                     /* dislocation density for distinct  */
                                     /* burgers vector groupings.  This   */
                                     /* array will be dynamically al-     */
                                     /* located to the appropriate size   */
                                     /* depending on the type of mobility */
                                     /* being used.                       */
        real8 disloDensity;

        real8 delSegLength;  /* accumulated length of deleted segments    */
                             /* since most recent write to result_density */
                             /* property file                             */

#ifdef CALCENERGY
        real8 TotalEnergy;   /* total energy                              */
        real8 FMMEnergy;     /* energy contribution from remote segments  */
        real8 PBCEnergy;     /* energy contribution from remote segments  */
#endif

#ifdef _ARLFEM
        real8 fem_delSegLength; /* Used for tracking dislocation length  */
                                /* lost as dislocations move outside the */
                                /* simulation bounds                     */
#endif

        real8 densityChange[14];  /* For tracking density change by burgers */
                                  /* vector; currently only used for BCC.   */
                                  /* Values accumulated since most recent   */
                                  /* write of density change data to the    */
                                  /* property file.                         */

        real8 TensionFactor;
        int elasticinteraction;

        real8 delpStrain[6],delSig[6],totpStn[6];
        real8 delpSpin[6],totpSpn[6];
#ifdef FIX_PLASTIC_STRAIN
        real8 delpStrainCorr[6], delpSpinCorr[6];
#endif

/*
 *      Added for strain decomposition and density flux decomp.
 *
 *      Note: The <elasticStrainTensor> is included in the control
 *            parameter file but is an output parameter only.  It
 *            is recalculated by the code each cycle and is written
 *            purely for informational purposes.
 */
        real8 totstraintensor[6];
        real8 elasticStrainTensor[6];
        real8 totedgepStrain[6], totscrewpStrain[6];
        real8 dedgepStrain[6],   dscrewpStrain[6];

/*
 *      For strain decomposition and density flux decomposition
 *      in BCC materials.
 */
        int  bcc_DensityFluxDecomp;  /* Density and flux decomposition flag
                                        for BCC materials. If equal to 1,
                                        density and flux are decomposed along
                                        3 planes: (110). If it equals to 2,
                                        6 planes decomposition: (110), (112). Default
                                        is 1: 3 planes decomposition.*/


        real8 BCC_LtotDecomp[4][11], BCC_dLtotDecomp[4][11];
        real8 BCC_fluxtot[4][7],     BCC_dfluxtot[4][7];

        real8 BCC2_LtotDecomp[4][20],  BCC2_dLtotDecomp[4][20];
        real8 BCC2_fluxtot[4][13],     BCC2_dfluxtot[4][13];

/*
 *      For strain decomposition and density flux decomp in
 *      rhombohedral vanadium.
 */
        real8 rhombo_Ltot[4][3][2],  rhombo_fluxtot[4][3][3];
        real8 rhombo_dLtot[4][3][2], rhombo_dfluxtot[4][3][3];

/*
 *      For strain decomposition and density flux decomp in FCC
 */
        real8 FCC_LtotDecomp[6][11], FCC_dLtotDecomp[6][11];
        real8 FCC_fluxtot[6][7],     FCC_dfluxtot[6][7];

/*
 *      For strain decomposition and density flux decomp in HCP
 */
        real8 HCP_LtotDecomp[10][17], HCP_dLtotDecomp[10][17];
        real8 HCP_fluxtot[10][9],     HCP_dfluxtot[10][9];

        int imgstrgrid[6];

/*
 *      The following two variables are used only temporarily and as such
 *      should not be specifiable in the control file, and hence, should
 *      not be bound to identifying strings via bindvar() as the other
 *      elements of this structure are.
 */
        char node_data_file[MAX_STRING_LEN];
                                  /* name of the file containing the    */
                                  /* nodal data for the run.  This data */
                                  /* is in the newer nodal data format, */
                                  /* (i.e. not contained as part of the */
                                  /* control file itself)               */

/*
 *      Define the parameters used in the nodal data file
 */
        int   dataFileVersion;       /* Version number of the data file */
        int   numFileSegments;       /* Number of files the nodal data  */
                                     /* is spread over.                 */
        int   nodeCount;             /* Total number of nodes in the    */
                                     /* data file (all file segments)   */
        int   dataDecompType;        /* Type of decomposition used when */
                                     /* generating the current data file*/
        int   dataDecompGeometry[3]; /* domain geometry used when gen-  */
                                     /* erating the current data file   */
        real8 minCoordinates[3];     /* minimum XYZ coordinates of sim  */
        real8 maxCoordinates[3];     /* maximum XYZ coordinates of sim  */

/*
 *      Define a couple factors used when calculating dislocation density.
 *      These will be calculated during initialization and will not be
 *      specified in the control parameter file.
 */
        real8 simVol;         /* Total volume of simulation */
        real8 burgVolFactor;  /* Volume factor used to convert dislocation */
                              /* length to dislocation density             */

/*
 *      Define the number of threads to use per MPI task. Ignored if
 *      threading has not been enabled.
 */
        int   maxNumThreads;

/*
 *      The following flag is used ONLY for some site-specific post-
 *      processing.  If set it will limit calculation of the flux
 *      decomposition to a subset of the simulation segments
 *      defined by certain constraints that are normally not
 *      used in the ParaDiS simulations!
 */
        int ignoreNonLoopFlux;

#ifdef BBFMM
        int  numGaussPoints;
        real8 ChebEps;
/*
 *      File names for kernel interaction matrix K, matrix of
 *      singular vectors U, and matrix of singular vectors V
 */
        char KMatFileName[MAX_STRING_LEN];
        char UMatFileName[MAX_STRING_LEN];
        char VMatFileName[MAX_STRING_LEN];
#endif

#ifdef UNIFORMFMM
        int  numGaussPoints;
/*
 *      File names for kernel interaction matrix K
 */
        char KMatFileName[MAX_STRING_LEN];
#endif

};

//------------------------------------------------------------------------------------------

extern void CtrlParamInit     (Param_t     *param, ParamList_t *CPList);
extern void DataParamInit     (Param_t     *param, ParamList_t *DPList);

extern void Print_Param       (Param_t     *param, const char *msg);
extern void Print_Param_List  (ParamList_t *plist, const char *msg);

#endif /* _PARAM_H */
