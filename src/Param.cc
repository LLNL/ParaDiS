/**************************************************************************
 *
 *      Module:       Param.c
 *      Description:  Contains functions for binding control and data
 *                    file parameter names to code variables.
 *
 *      Public functions:
 *          CtrlParamInit()
 *          DataParamInit()
 *          MarkParamDisabled()
 *          MarkParamEnabled()
 *
 *************************************************************************/
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Parse.h"
#include "Param.h"

/*------------------------------------------------------------------------
 *
 *      Functions:    MarkParamDisabled
 *      Description:  Explicitly mark a control file parameter
 *                    to be disabled.  Primarily used to
 *                    prevent writing to the restart file any
 *                    control file parameters that are not
 *                    applicable or appropriate to the current
 *                    execution of the code.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure
 *                    associated with the control file parameters.
 *          name      String containing the name of the control
 *                    parameter to be explicitly disabled.
 *
 *----------------------------------------------------------------------*/
void MarkParamDisabled(ParamList_t *CPList, const char *name)
{
    int indx = LookupParam(CPList,name);

    if (indx >= 0) CPList->varList[indx].flags |= VFLAG_DISABLED;
    else           Fatal("MarkParamDisabled: unknown parameter %s\n", name);

    return;
}

/*------------------------------------------------------------------------
 *
 *      Functions:    MarkParamEnabled
 *      Description:  Explicitly mark a control file parameter
 *                    to be enabled.  Primarily used when
 *                    explicitly controlling whether a control
 *                    parameter is to be written to a restart
 *                    file.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure
 *                    associated with the control file parameters.
 *          name      String containing the name of the control
 *                    parameter to be explicitly enabled.
 *
 *----------------------------------------------------------------------*/
void MarkParamEnabled(ParamList_t *CPList, const char *name)
{
    int indx = LookupParam(CPList,name);

    if (indx >= 0)
        CPList->varList[indx].flags &= ~VFLAG_DISABLED;
    else
        Fatal("MarkParamEnabled: unknown parameter %s\n", name);

    return;
}

/*------------------------------------------------------------------------
 *
 *      Function:     CtrlParamInit
 *      Description:  Bind all valid control file parameters to
 *                    the associated code variables.
 *
 *      Arguments:
 *          CPList    Pointer to the parameter list structure
 *                    associated with the control file parameters.
 *
 *      Last Modified: 04/08/2008 gh - Added bindings for new parameters
 *                                     <MobGo>, <MobP>, and <inclusionFile>.
 *                                     binding for pre-existing variable
 *                                     <MobLinear>.
 *----------------------------------------------------------------------*/
void CtrlParamInit(Param_t *param, ParamList_t *CPList)
{
    // Note: Parameters need only be initialized if their
    //       default values are non-zero.

    CPList->paramCnt = 0;
    CPList->varList  = (VarData_t *) NULL;

    // Simulation cell and processor setup

    BindVar(CPList, "Simulation cell and processor setup", (void *) NULL, V_COMMENT, 1, VFLAG_NULL);

    BindVar(CPList, "numXdoms", &param->nXdoms, V_INT, 1, VFLAG_NULL); param->nXdoms = 1;
    BindVar(CPList, "numYdoms", &param->nYdoms, V_INT, 1, VFLAG_NULL); param->nYdoms = 1;
    BindVar(CPList, "numZdoms", &param->nZdoms, V_INT, 1, VFLAG_NULL); param->nZdoms = 1;

#if defined _BGP || defined _BGQ
    // On BlueGene, user-specified task distribution is used by default;
    // system does not try to reset it.
    BindVar(CPList, "taskMappingMode", &param->taskMappingMode, V_INT, 1, VFLAG_NULL); param->taskMappingMode = 1;
#endif

    BindVar(CPList, "numXcells" , &param->nXcells   , V_INT, 1, VFLAG_NULL); param->nXcells    = 3;   // Must be >= 3
    BindVar(CPList, "numYcells" , &param->nYcells   , V_INT, 1, VFLAG_NULL); param->nYcells    = 3;   // Must be >= 3
    BindVar(CPList, "numZcells" , &param->nZcells   , V_INT, 1, VFLAG_NULL); param->nZcells    = 3;   // Must be >= 3

    BindVar(CPList, "xBoundType", &param->xBoundType, V_INT, 1, VFLAG_NULL); param->xBoundType = Periodic;
    BindVar(CPList, "yBoundType", &param->yBoundType, V_INT, 1, VFLAG_NULL); param->yBoundType = Periodic;
    BindVar(CPList, "zBoundType", &param->zBoundType, V_INT, 1, VFLAG_NULL); param->zBoundType = Periodic;

    BindVar(CPList, "xBoundMin" , &param->xBoundMin , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "xBoundMax" , &param->xBoundMax , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "yBoundMin" , &param->yBoundMin , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "yBoundMax" , &param->yBoundMax , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "zBoundMin" , &param->zBoundMin , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "zBoundMax" , &param->zBoundMax , V_DBL, 1, VFLAG_NULL);

    BindVar(CPList, "decompType", &param->decompType, V_INT, 1, VFLAG_NULL); param->decompType = 2;
    BindVar(CPList, "DLBfreq"   , &param->DLBfreq   , V_INT, 1, VFLAG_NULL); param->DLBfreq    = 3;

    // Simulation time and timestepping controls

    BindVar(CPList, "Simulation time and timestepping controls", (void *) NULL, V_COMMENT, 1, VFLAG_NULL);

    BindVar(CPList, "cycleStart"            , &param->cycleStart            , V_INT   , 1, VFLAG_NULL);
    BindVar(CPList, "maxstep"               , &param->maxstep               , V_INT   , 1, VFLAG_NULL); param->maxstep = 100;
    BindVar(CPList, "timeNow"               , &param->timeNow               , V_DBL   , 1, VFLAG_NULL);
    BindVar(CPList, "timeStart"             , &param->timeStart             , V_DBL   , 1, VFLAG_NULL);
    BindVar(CPList, "timeEnd"               , &param->timeEnd               , V_DBL   , 1, VFLAG_NULL); param->timeEnd = -1.0;
    BindVar(CPList, "timestepIntegrator"    ,  param->timestepIntegrator    , V_STRING, 1, VFLAG_NULL); strcpy(param->timestepIntegrator, "trapezoid");
    BindVar(CPList, "trapezoidMaxIterations", &param->trapezoidMaxIterations, V_INT   , 1, VFLAG_NULL); param->trapezoidMaxIterations = 2;
    BindVar(CPList, "deltaTT"               , &param->deltaTT               , V_DBL   , 1, VFLAG_NULL);
    BindVar(CPList, "maxDT"                 , &param->maxDT                 , V_DBL   , 1, VFLAG_NULL); param->maxDT = 1.0e-07;
    BindVar(CPList, "nextDT"                , &param->nextDT                , V_DBL   , 1, VFLAG_NULL);
    BindVar(CPList, "dtIncrementFact"       , &param->dtIncrementFact       , V_DBL   , 1, VFLAG_NULL); param->dtIncrementFact = 1.2;
    BindVar(CPList, "dtDecrementFact"       , &param->dtDecrementFact       , V_DBL   , 1, VFLAG_NULL); param->dtDecrementFact = 0.5;
    BindVar(CPList, "dtExponent"            , &param->dtExponent            , V_DBL   , 1, VFLAG_NULL); param->dtExponent      = 4.0;
    BindVar(CPList, "dtVariableAdjustment"  , &param->dtVariableAdjustment  , V_INT   , 1, VFLAG_NULL);
    BindVar(CPList, "rTol"                  , &param->rTol                  , V_DBL   , 1, VFLAG_NULL); param->rTol = -1.0;
    BindVar(CPList, "rmax"                  , &param->rmax                  , V_DBL   , 1, VFLAG_NULL); param->rmax = 100;

    //  GPU controls...

    BindVar(CPList, "gpu_enabled"           , &param->gpu_enabled           , V_INT   , 1, VFLAG_NULL); param->gpu_enabled = 0;            // 0=disable, 1=enable GPU compute kernels
    BindVar(CPList, "gpu_nstrm"             , &param->gpu_nstrm             , V_INT   , 1, VFLAG_NULL); param->gpu_nstrm   = 4;            // default CUDA stream count=4
    BindVar(CPList, "gpu_nthrd"             , &param->gpu_nthrd             , V_INT   , 1, VFLAG_NULL); param->gpu_nthrd   = 256;          // number of threads per kernel launch
    BindVar(CPList, "gpu_npblk"             , &param->gpu_npblk             , V_INT   , 1, VFLAG_NULL); param->gpu_npblk   = 64*1024;      // number of segment pairs per CUDA stream block

    // Sundials/KINSOL controls...

    BindVar(CPList, "KINSOL_UseNewton"         , &param->KINSOL_UseNewton        , V_INT    , 1, VFLAG_NULL ); param->KINSOL_UseNewton         = 0;     // 0=disable, 1=enable newton
    BindVar(CPList, "KINSOL_MaxIterations"     , &param->KINSOL_MaxIterations    , V_INT    , 1, VFLAG_NULL ); param->KINSOL_MaxIterations     = 5;     //
    BindVar(CPList, "KINSOL_MaxLinIterations"  , &param->KINSOL_MaxLinIterations , V_INT    , 1, VFLAG_NULL ); param->KINSOL_MaxLinIterations  = 0;     //
    BindVar(CPList, "KINSOL_NumPriorResiduals" , &param->KINSOL_NumPriorResiduals, V_INT    , 1, VFLAG_NULL ); param->KINSOL_NumPriorResiduals = param->KINSOL_MaxIterations-1;
    BindVar(CPList, "KINSOL_PrintLevel"        , &param->KINSOL_PrintLevel       , V_INT    , 1, VFLAG_NULL ); param->KINSOL_PrintLevel        = 0;     // 0=disable
    BindVar(CPList, "KINSOL_LogFile"           ,  param->KINSOL_LogFile          , V_STRING , 1, VFLAG_NULL ); memset(param->KINSOL_LogFile,0,sizeof(param->KINSOL_LogFile));  // (null string, must be set within control file)
    BindVar(CPList, "KINSOL_EtaVal"            , &param->KINSOL_EtaVal           , V_DBL    , 1, VFLAG_NULL ); param->KINSOL_EtaVal            = 0.0;   // note: 0.0=forces the value of the nonlinear residual factor to the KINSOL default to 0.1

    // Sundials/ARKode controls...

    BindVar(CPList, "ARKODE_IRKtable"          , &param->ARKODE_IRKtable         , V_INT    , 1, VFLAG_NULL ); param->ARKODE_IRKtable          = 13;    // Billington-3-3-2 Method (3rd Order)
    BindVar(CPList, "ARKODE_FPsolver"          , &param->ARKODE_FPsolver         , V_INT    , 1, VFLAG_NULL ); param->ARKODE_FPsolver          = 1;     // 0=Newton solver, 1=fixed-point solver
    BindVar(CPList, "ARKODE_MaxNonLinIters"    , &param->ARKODE_MaxNonLinIters   , V_INT    , 1, VFLAG_NULL ); param->ARKODE_MaxNonLinIters    = 4;     // 0=use ARKode's default of 3 Newton or 10 FP iterations
    BindVar(CPList, "ARKODE_MaxLinIters"       , &param->ARKODE_MaxLinIters      , V_INT    , 1, VFLAG_NULL ); param->ARKODE_MaxLinIters       = 0;     // 0=use ARKode's default of 5 iterations
    BindVar(CPList, "ARKODE_FPaccel"           , &param->ARKODE_FPaccel          , V_INT    , 1, VFLAG_NULL ); param->ARKODE_FPaccel           = param->ARKODE_MaxNonLinIters-1;   // number of acceleration vectors
    BindVar(CPList, "ARKODE_PredictorMethod"   , &param->ARKODE_PredictorMethod  , V_INT    , 1, VFLAG_NULL ); param->ARKODE_PredictorMethod   = 0;     // (trivial predictor)
    BindVar(CPList, "ARKODE_AdaptivityMethod"  , &param->ARKODE_AdaptivityMethod , V_INT    , 1, VFLAG_NULL ); param->ARKODE_AdaptivityMethod  = 2;     // (I step controller)
    BindVar(CPList, "ARKODE_NonLinCoef"        , &param->ARKODE_NonLinCoef       , V_DBL    , 1, VFLAG_NULL ); param->ARKODE_NonLinCoef        = 1.0;   // (nonlinear tolerance factor)
    BindVar(CPList, "ARKODE_LinCoef"           , &param->ARKODE_LinCoef          , V_DBL    , 1, VFLAG_NULL ); param->ARKODE_LinCoef           = 0.9;   // (linear    tolerance factor)

    // Discretization controls and controls for topological changes

    BindVar(CPList, "Discretization and topological change controls", (void *) NULL, V_COMMENT, 1, VFLAG_NULL);
    BindVar(CPList, "maxSeg"                   , &param->maxSeg                   , V_DBL, 1, VFLAG_NULL); param->maxSeg              = -1.0;
    BindVar(CPList, "minSeg"                   , &param->minSeg                   , V_DBL, 1, VFLAG_NULL); param->minSeg              = -1.0;
    BindVar(CPList, "rann"                     , &param->rann                     , V_DBL, 1, VFLAG_NULL); param->rann                = -1.0;
    BindVar(CPList, "remeshRule"               , &param->remeshRule               , V_INT, 1, VFLAG_NULL); param->remeshRule          =  2;
    BindVar(CPList, "splitMultiNodeFreq"       , &param->splitMultiNodeFreq       , V_INT, 1, VFLAG_NULL); param->splitMultiNodeFreq  =  1;
    BindVar(CPList, "splitMultiNodeAlpha"      , &param->splitMultiNodeAlpha      , V_DBL, 1, VFLAG_NULL); param->splitMultiNodeAlpha =  1.0e-3;
    BindVar(CPList, "useParallelSplitMultiNode", &param->useParallelSplitMultiNode, V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "useUniformJunctionLen"    , &param->useUniformJunctionLen    , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "collisionMethod"          , &param->collisionMethod          , V_INT, 1, VFLAG_NULL); param->collisionMethod     =  2;

    // FMM controls

    BindVar(CPList, "Fast Multipole Method controls", (void *) NULL, V_COMMENT, 1, VFLAG_NULL);
    BindVar(CPList, "fmEnabled"           , &param->fmEnabled           , V_INT   , 1, VFLAG_NULL);
    BindVar(CPList, "fmMPOrder"           , &param->fmMPOrder           , V_INT   , 1, VFLAG_NULL); param->fmMPOrder     = 2;

    // Treat "fmTaylorOrder" as an alias for fmExpansionOrder
    BindVar(CPList, "fmTaylorOrder"       , &param->fmExpansionOrder    , V_INT   , 1, VFLAG_NULL);
    BindVar(CPList, "fmExpansionOrder"    , &param->fmExpansionOrder    , V_INT   , 1, VFLAG_NULL); param->fmExpansionOrder = 5;

#ifdef TAYLORFMM
    BindVar(CPList, "fmCorrectionTbl"     ,  param->fmCorrectionTbl     , V_STRING, 1, VFLAG_NULL); strcpy(param->fmCorrectionTbl,"inputs/fm-ctab.Ta.600K.0GPa.m2.t5.dat");
#else
    BindVar(CPList, "fmCorrectionTbl"     ,  param->fmCorrectionTbl     , V_STRING, 1, VFLAG_NULL); strcpy(param->fmCorrectionTbl,"inputs/Stanfordfm-ctab.c3l3.dat");
#endif

    BindVar(CPList, "fmEnableAnisotropy"  , &param->fmEnableAnisotropy  , V_INT   , 1, VFLAG_NULL);

#ifdef BBFMM
        BindVar(CPList, "numGaussPoints", &param->numGaussPoints, V_INT,1, VFLAG_NULL);param->numGaussPoints = 5;
        BindVar(CPList, "KMatFileName", param->KMatFileName, V_STRING,1, VFLAG_NULL);strcpy(param->KMatFileName,"inputs/bbfmmKmatrix.mat");
        BindVar(CPList, "UMatFileName", param->UMatFileName, V_STRING,1, VFLAG_NULL);strcpy(param->UMatFileName,"inputs/bbfmmUmatrix.mat");
        BindVar(CPList, "VMatFileName", param->VMatFileName, V_STRING,1, VFLAG_NULL);strcpy(param->VMatFileName,"inputs/bbfmmVmatrix.mat");
        BindVar(CPList, "ChebEps", &param->ChebEps, V_DBL, 1, VFLAG_NULL);param->ChebEps = 1e-9;
#endif

#ifdef UNIFORMFMM
        BindVar(CPList, "numGaussPoints", &param->numGaussPoints, V_INT, 1, VFLAG_NULL);param->numGaussPoints = 5;
        BindVar(CPList, "KMatFileName", param->KMatFileName, V_STRING,1, VFLAG_NULL);strcpy(param->KMatFileName,"inputs/uniformfmmKmatrix.mat");
#endif

#if (defined ANISOTROPIC) && (defined TAYLORFMM)
    BindVar(CPList, "fmAnisoDerivTbl"     , param->fmAnisoDerivTbl      , V_STRING, 1, VFLAG_NULL); strcpy(param->fmAnisoDerivTbl,"inputs/fm-aniso-deriv-table.dat");
#endif

#ifdef ESHELBY
   BindVar(CPList, "eshelbyfmEnabled"          , &param->eshelbyfmEnabled    , V_INT   , 1, VFLAG_NULL );
   BindVar(CPList, "eshelbyfmMPOrder"          , &param->eshelbyfmMPOrder    , V_INT   , 1, VFLAG_NULL );   param->eshelbyfmMPOrder     = 3;
   BindVar(CPList, "eshelbyfmTaylorOrder"      , &param->eshelbyfmTaylorOrder, V_INT   , 1, VFLAG_NULL );   param->eshelbyfmTaylorOrder = 3;
   BindVar(CPList, "RefineSegCollision"        , &param->RefineSegCollision  , V_INT   , 1, VFLAG_NULL );   param->RefineSegCollision   = 0;
   BindVar(CPList, "MobEshelbyResist"          , &param->MobEshelbyResist    , V_DBL   , 1, VFLAG_NULL );   param->MobEshelbyResist     =  1.0e7; // Default : max resistance to dislocation motion in the particles
   BindVar(CPList, "shearModulus2"             , &param->shearModulus2       , V_DBL   , 1, VFLAG_NULL );   param->shearModulus2        = -1.0;
   BindVar(CPList, "pois2"                     , &param->pois2               , V_DBL   , 1, VFLAG_NULL );   param->pois2                = -1.0;
   BindVar(CPList, "Ecore2"                    , &param->Ecore2              , V_DBL   , 1, VFLAG_NULL );   param->Ecore2               = -1.0;

   BindVar(CPList, "IncStepEnergy"             , &param->IncStepEnergy       , V_DBL   , 1, VFLAG_NULL );   param->IncStepEnergy        = -1.0;
#endif

#ifdef CALCENERGY
   BindVar(CPList, "fmECorrTbl"                ,  param->fmECorrTbl          , V_STRING, 1, VFLAG_NULL );   strcpy(param->fmECorrTbl,"inputs/fmEnergy.m2.t5.fmm");
#endif

#ifdef FRCRIT
   BindVar(CPList, "minBisection"              , &param->minBisection        , V_DBL   , 1, VFLAG_NULL );   param->minBisection         = 0.0;
   BindVar(CPList, "maxBisection"              , &param->maxBisection        , V_DBL   , 1, VFLAG_NULL );   param->maxBisection         = 0.0;
#endif

#ifdef SPECTRAL
   BindVar(CPList, "FFTgrid"                   , &param->FFTgrid             , V_INT   , 1, VFLAG_NULL );   param->FFTgrid       = 64;
   BindVar(CPList, "FFTenabled"                , &param->FFTenabled          , V_INT   , 1, VFLAG_NULL );   param->FFTenabled    = 0;
   BindVar(CPList, "localRotation"             , &param->localRotation       , V_INT   , 1, VFLAG_NULL );   param->localRotation = 0;
#endif

    // Identify tables needed for remote force calculations if FMM is not enabled

    BindVar(CPList, "Tables for non-FMM far-field force calcs", (void *) NULL, V_COMMENT, 1, VFLAG_NULL);
    BindVar(CPList, "Rijmfile"     , param->Rijmfile   , V_STRING, 1, VFLAG_NULL); strcpy(param->Rijmfile   ,"inputs/Rijm.cube.out"   );
    BindVar(CPList, "RijmPBCfile"  , param->RijmPBCfile, V_STRING, 1, VFLAG_NULL); strcpy(param->RijmPBCfile,"inputs/RijmPBC.cube.out");

    // Loading condition parameters

    BindVar(CPList, "Loading conditions", (void *) NULL, V_COMMENT, 1, VFLAG_NULL);
    BindVar(CPList, "TempK"          , &param->TempK          , V_DBL, 1, VFLAG_NULL); param->TempK = 600;
    BindVar(CPList, "pressure"       , &param->pressure       , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "loadType"       , &param->loadType       , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "appliedStress"  , param->appliedStress   , V_DBL, 6, VFLAG_NULL);
    BindVar(CPList, "eRate"          , &param->eRate          , V_DBL, 1, VFLAG_NULL); param->eRate           = 1.0;
    BindVar(CPList, "indxErate"      , &param->indxErate      , V_INT, 1, VFLAG_NULL); param->indxErate       = 0;
    BindVar(CPList, "edotdir"        , param->edotdir         , V_DBL, 3, VFLAG_NULL); param->edotdir[0]      = 1.0;
                                                                                       param->edotdir[1]      = 0.0;
                                                                                       param->edotdir[2]      = 0.0;

    BindVar(CPList, "crystalRotation", &param->crystalRotation, V_INT, 1, VFLAG_NULL); param->crystalRotation = 1;
    BindVar(CPList, "cubicModulus"   , &param->cubicModulus   , V_INT, 1, VFLAG_NULL); param->cubicModulus    = 0;
    BindVar(CPList, "cubicC11"       , &param->cubicC11       , V_DBL, 1, VFLAG_NULL); param->cubicC11        = -1.0;
    BindVar(CPList, "cubicC12"       , &param->cubicC12       , V_DBL, 1, VFLAG_NULL); param->cubicC12        = -1.0;
    BindVar(CPList, "cubicC44"       , &param->cubicC44       , V_DBL, 1, VFLAG_NULL); param->cubicC44        = -1.0;

    BindVar(CPList, "cTimeOld"       , &param->cTimeOld       , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "dCyclicStrain"  , &param->dCyclicStrain  , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "netCyclicStrain", &param->netCyclicStrain, V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "numLoadCycle"   , &param->numLoadCycle   , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "eAmp"           , &param->eAmp           , V_DBL, 1, VFLAG_NULL);

    // To specify input sample axes in laboratory frame

    BindVar(CPList, "useLabFrame"    , &param->useLabFrame    , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "rotationMatrix" ,  param->rotMatrix      , V_DBL, 9, VFLAG_NULL);

    // Initialize the rotation matrix to the identity matrix

    for (int row=0; row<3; ++row)
        for (int col=0; col<3; ++col)
            param->rotMatrix[row][col] = ( (row==col) ? 1.0 : 0.0 );

    // Parameters for material specific constants and mobility values

    BindVar(CPList, "Material and mobility parameters", (void *) NULL, V_COMMENT, 1, VFLAG_NULL);

    BindVar(CPList, "mobilityLaw"          , param->mobilityLaw          , V_STRING,  1, VFLAG_NULL); strcpy(param->mobilityLaw     , "BCC_0");
    BindVar(CPList, "materialTypeName"     , param->materialTypeName     , V_STRING,  1, VFLAG_NULL); strcpy(param->materialTypeName, MAT_TYPE_NAME_UNDEFINED);
    BindVar(CPList, "elasticConstantMatrix", param->elasticConstantMatrix, V_DBL   , 36, VFLAG_NULL);

#if defined ANISOTROPIC
    BindVar(CPList, "C11", &param->C11, V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "C12", &param->C12, V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "C13", &param->C13, V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "C33", &param->C33, V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "C44", &param->C44, V_DBL, 1, VFLAG_NULL);

    BindVar(CPList, "anisoHarmonicsNumTermsBase", &param->anisoHarmonicsNumTermsBase, V_INT, 1, VFLAG_NULL); param->anisoHarmonicsNumTermsBase = 3;
    BindVar(CPList, "anisoNumThetaPoints"       , &param->anisoNumThetaPoints       , V_INT, 1, VFLAG_NULL); param->anisoNumThetaPoints        = 150;
    BindVar(CPList, "anisoNumPhiPoints"         , &param->anisoNumPhiPoints         , V_INT, 1, VFLAG_NULL); param->anisoNumPhiPoints          = 150;
#endif

#ifdef ESHELBY
    BindVar(CPList, "inclusionFile"         ,  param->inclusionFile         , V_STRING, 1, VFLAG_NULL );
#endif

    BindVar(CPList, "ecFile"                ,  param->ecFile                , V_STRING, 1, VFLAG_NULL);
    BindVar(CPList, "meltTempFile"          ,  param->meltTempFile          , V_STRING, 1, VFLAG_NULL);
    BindVar(CPList, "meltTemp"              , &param->meltTemp              , V_DBL   , 1, VFLAG_NULL); param->meltTemp               =  3.170683e+03;  // Ta: @ temp=600K, pressure=0GPa
    BindVar(CPList, "shearModulus"          , &param->shearModulus          , V_DBL   , 1, VFLAG_NULL); param->shearModulus           =  6.488424e+10;  // Ta: units=Pa @t=600K, p=0Pa
    BindVar(CPList, "pois"                  , &param->pois                  , V_DBL   , 1, VFLAG_NULL); param->pois                   =  3.327533e-01;  // Ta: temp=600K, pressure=0Pa
    BindVar(CPList, "burgMag"               , &param->burgMag               , V_DBL   , 1, VFLAG_NULL); param->burgMag                =  2.875401e-10;  // Ta: units=m @t=600K, pressure=0Pa
    BindVar(CPList, "ecrit"                 , &param->ecrit                 , V_DBL   , 1, VFLAG_NULL); param->ecrit                  =  1.0e-4;        // Critical angle for segment pair parallelism test
    BindVar(CPList, "cOVERa"                , &param->cOVERa                , V_DBL   , 1, VFLAG_NULL); param->cOVERa                 =  1.5680;        // Beryllium
    BindVar(CPList, "YoungModulus"          , &param->YoungsModulus         , V_DBL   , 1, VFLAG_NULL); param->YoungsModulus          = 172.95e+09;     // Ta: units=Pa @t=600K, p=0Pa

    BindVar(CPList, "vacancyConc"           , &param->vacancyConc           , V_DBL   , 1, VFLAG_NULL); param->vacancyConc            = -1.0;           // signifies not used
    BindVar(CPList, "vacancyConcEquilibrium", &param->vacancyConcEquilibrium, V_DBL   , 1, VFLAG_NULL); param->vacancyConcEquilibrium = -1.0;           // signifies not used

#if defined ANISOTROPIC
    // If anisotropic elasticity is enabled, elastic modulii consistent
    // with the elastic constants above must either be provided or
    // calculated later from those constants.  Therefore, in this case,
    // no defaults are permitted for these values.

    param->shearModulus  = 0.0;
    param->pois          = 0.0;
    param->YoungsModulus = 0.0;
#endif

    BindVar(CPList, "rc"            , &param->rc           , V_DBL   , 1, VFLAG_NULL); param->rc     = -1.0;
    BindVar(CPList, "Ecore"         , &param->Ecore        , V_DBL   , 1, VFLAG_NULL); param->Ecore  = -1.0;
    BindVar(CPList, "MDcore"        , &param->MDcore       , V_INT   , 1, VFLAG_NULL); param->MDcore = 0;
    BindVar(CPList, "MDcoreMaterial", param->MDcoreMaterial, V_STRING, 1, VFLAG_NULL); strcpy(param->MDcoreMaterial, "W");

    // For HCP, there are 3 Ecore values which each apply to specific
    // types of burgers vectors (types <a>, <c> and <c+a>) and are
    // provided in the HCPEcoreA,  HCPEcoreC, HCPEcoreCpA control parameters.
    //
    // These core energy prefactors are dependent on the core radius entered,
    // and in units of the shear modulus (energy per unit length / b^2 units).
    //
    // They replace the default Ecore factor defined in Initialize.cc
    //
    // Default values are for Berylium:  a = 3 |bA| ~ 7.8 A
    //                 HCPEcoreA   = 18.06977 GPa/|b<a>|^2
    //                 HCPEcoreC   = 21.70031 GPa/|b<c>|^2
    //                 HCPEcoreCpA =  7.18520 GPa/|b<c+a>|^2
    //
    // In the SelfForce*() functions, these Ecore values are multiplied
    // by b^2 i.e. A^2.  Giving core forces the right units of Newton per
    // length or energy per unit length.

    BindVar(CPList, "HCPEcoreA"  , &param->HCPEcoreA   , V_DBL, 1, VFLAG_NULL); param->HCPEcoreA   = 18.06977e9; // units: Pa/|b<a>|^2
    BindVar(CPList, "HCPEcoreC"  , &param->HCPEcoreC   , V_DBL, 1, VFLAG_NULL); param->HCPEcoreC   = 21.70031e9; // units: Pa/|b<c>|^2
    BindVar(CPList, "HCPEcoreCpA", &param->HCPEcoreCpA , V_DBL, 1, VFLAG_NULL); param->HCPEcoreCpA =  7.18520e9; // units: Pa/|b<c+a>|^2
    BindVar(CPList, "MobScrew"   , &param->MobScrew    , V_DBL, 1, VFLAG_NULL); param->MobScrew    = 10.0;       // units: 1/(Pa*s)
    BindVar(CPList, "MobEdge"    , &param->MobEdge     , V_DBL, 1, VFLAG_NULL); param->MobEdge     = 10.0;       // units: 1/(Pa*s)
    BindVar(CPList, "MobClimb"   , &param->MobClimb    , V_DBL, 1, VFLAG_NULL); param->MobClimb    =  1.0e-02;   // units: 1/(Pa*s)

    BindVar(CPList, "MobClimbJunc"   , &param->MobClimbJunc    , V_DBL, 1, VFLAG_NULL); param->MobClimbJunc  =  1.0;    // units: 1/(Pa*s)
    BindVar(CPList, "MobLineJunc"    , &param->MobLineJunc     , V_DBL, 1, VFLAG_NULL); param->MobLineJunc   = -1.0;    // units: 1/(Pa*s)
    BindVar(CPList, "MobSplit"       , &param->MobSplit        , V_INT, 1, VFLAG_NULL); param->MobSplit      =  0;      //
    BindVar(CPList, "MobClimbSplit"  , &param->MobClimbSplit   , V_DBL, 1, VFLAG_NULL); param->MobClimbSplit = -1.0;    // units: 1/(Pa*s)
    BindVar(CPList, "shearVelocity"  , &param->shearVelocity   , V_DBL, 1, VFLAG_NULL); param->shearVelocity =  3400.0; // units: m/s

    // Angle-dependent mobility law
    BindVar(CPList, "MobTheta_cnt"   , &param->MobTheta_cnt    , V_INT, 1, VFLAG_NULL); param->MobTheta_cnt = 0;
    BindVar(CPList, "MobThetas"      , &param->MobThetas       , V_DBL, ANGLE_MOB_CNT, VFLAG_NULL);
    BindVar(CPList, "MobGlide_B"     , &param->MobGlide_B      , V_DBL, ANGLE_MOB_CNT, VFLAG_NULL);
    BindVar(CPList, "MobGlide_D"     , &param->MobGlide_D      , V_DBL, ANGLE_MOB_CNT, VFLAG_NULL);
    BindVar(CPList, "MobGlide_V0"    , &param->MobGlide_V0     , V_DBL, ANGLE_MOB_CNT, VFLAG_NULL);

    // Non-linear mobility only...

    BindVar(CPList, "MobFileNL"  ,  param->MobFileNL   , V_STRING, 1, VFLAG_NULL);
    BindVar(CPList, "MobGlide"   , &param->MobGlide    , V_DBL, 1, VFLAG_NULL); param->MobGlide    = 1.0e+02;       // In units of 1/(Pa*s)
    BindVar(CPList, "MobLine"    , &param->MobLine     , V_DBL, 1, VFLAG_NULL); param->MobLine     = 1.0e+04;       // In units of 1/(Pa*s)
    BindVar(CPList, "MobDeltaH0" , &param->MobDeltaH0  , V_DBL, 1, VFLAG_NULL); param->MobDeltaH0  = 1.728748e-19;  // Ta: @ pressure=0.0GPa, temp=600K
    BindVar(CPList, "MobPeierls" , &param->MobPeierls  , V_DBL, 1, VFLAG_NULL); param->MobPeierls  = 3.240000e+08;  // Ta: @ pressure=0.0Gpa, temp=600k
    BindVar(CPList, "MobPieirls" , &param->MobPeierls  , V_DBL, 1, VFLAG_ALIAS);
    BindVar(CPList, "MobAlpha"   , &param->MobAlpha    , V_DBL, 1, VFLAG_NULL); param->MobAlpha    = 9.661236e-01;  // Ta: @ pressure=0.0GPa, temp=600K
    BindVar(CPList, "MobExpP"    , &param->MobExpP     , V_DBL, 1, VFLAG_NULL); param->MobExpP     = 0.5;
    BindVar(CPList, "MobExpQ"    , &param->MobExpQ     , V_DBL, 1, VFLAG_NULL); param->MobExpQ     = 1.23;
    BindVar(CPList, "MobCoeffA"  , &param->MobCoeffA   , V_DBL, 1, VFLAG_NULL); param->MobCoeffA   = 1.525e+00;     // Tantalum
    BindVar(CPList, "MobCoeffC"  , &param->MobCoeffC   , V_DBL, 1, VFLAG_NULL); param->MobCoeffC   = 2.03237e-01;   // Tantalum
    BindVar(CPList, "MobCoeffC0" , &param->MobCoeffC0  , V_DBL, 1, VFLAG_NULL); param->MobCoeffC0  = 2.048e+03;     // Tantalum: units of meter/sec
    BindVar(CPList, "MobCoeffD"  , &param->MobCoeffD   , V_DBL, 1, VFLAG_NULL); param->MobCoeffD   = 1.62695e-04;   // Tantalum: units of GPa

    // For Vanadium non-linear mobility only

    BindVar(CPList, "MobP"       , &param->pressure    , V_DBL, 1, VFLAG_ALIAS);
    BindVar(CPList, "MobG0"      , &param->MobG0       , V_DBL, 1, VFLAG_NULL );
    BindVar(CPList, "MobGo"      , &param->MobG0       , V_DBL, 1, VFLAG_ALIAS); param->MobG0 = 49.0e+09;  // units = GPa

    // For 'faceted' mobility functions only

    BindVar(CPList, "MobNumTrenches"     , &param->MobNumTrenches    , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "MobTrench"          ,  param->MobTrench         , V_DBL, MAX_MOB_TRENCHES, VFLAG_NULL);
    BindVar(CPList, "MobTrenchAngle"     ,  param->MobTrenchAngle    , V_DBL, MAX_MOB_TRENCHES, VFLAG_NULL);
    BindVar(CPList, "MobTrenchWidth"     ,  param->MobTrenchWidth    , V_DBL, MAX_MOB_TRENCHES, VFLAG_NULL);
    BindVar(CPList, "MobFacetedUseCusps" , &param->MobFacetedUseCusps, V_INT, 1, VFLAG_NULL);

    // For HCP mobility law functions only.
    //
    // Drag coefficients defaults are in units of Pa*s and are currently
    // based on values for Beryllium from MD simulations.

    BindVar(CPList, "HCP_A_Basal_EdgeDrag"          , &param->HCP_A_Basal_EdgeDrag          , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_A_Basal_ScrewDrag"         , &param->HCP_A_Basal_ScrewDrag         , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_A_Prismatic_EdgeDrag"      , &param->HCP_A_Prismatic_EdgeDrag      , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_A_Prismatic_ScrewDrag"     , &param->HCP_A_Prismatic_ScrewDrag     , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_A_1stPyramidal_EdgeDrag"   , &param->HCP_A_1stPyramidal_EdgeDrag   , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_A_1stPyramidal_ScrewDrag"  , &param->HCP_A_1stPyramidal_ScrewDrag  , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_A_2ndPyramidal_EdgeDrag"   , &param->HCP_A_2ndPyramidal_EdgeDrag   , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_A_2ndPyramidal_ScrewDrag"  , &param->HCP_A_2ndPyramidal_ScrewDrag  , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_CpA_Prismatic_EdgeDrag"    , &param->HCP_CpA_Prismatic_EdgeDrag    , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_CpA_Prismatic_ScrewDrag"   , &param->HCP_CpA_Prismatic_ScrewDrag   , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_CpA_1stPyramidal_EdgeDrag" , &param->HCP_CpA_1stPyramidal_EdgeDrag , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_CpA_1stPyramidal_ScrewDrag", &param->HCP_CpA_1stPyramidal_ScrewDrag, V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_CpA_2ndPyramidal_EdgeDrag" , &param->HCP_CpA_2ndPyramidal_EdgeDrag , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_CpA_2ndPyramidal_ScrewDrag", &param->HCP_CpA_2ndPyramidal_ScrewDrag, V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_C_Prismatic_EdgeDrag"      , &param->HCP_C_Prismatic_EdgeDrag      , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_C_Prismatic_ScrewDrag"     , &param->HCP_C_Prismatic_ScrewDrag     , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_Sessile_EdgeDrag"          , &param->HCP_Sessile_EdgeDrag          , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_Sessile_ScrewDrag"         , &param->HCP_Sessile_ScrewDrag         , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "HCP_LineDrag"                  , &param->HCP_LineDrag                  , V_DBL, 1, VFLAG_NULL);

    // Drag coefficients are in units of Pa/s and currently based on
    // values from MD simulations.

    param->HCP_A_Basal_EdgeDrag           = 1.18434387898773e-04;
    param->HCP_A_Basal_ScrewDrag          = 1.65421659553579e-04;

    param->HCP_A_Prismatic_EdgeDrag       = 5.69513923115224e-05;
    param->HCP_A_Prismatic_ScrewDrag      = 7.15096201802962e-05;

    param->HCP_A_1stPyramidal_EdgeDrag    = 7.28268834977848e-05;
    param->HCP_A_1stPyramidal_ScrewDrag   = 1.05154866364405e-04;

    param->HCP_A_2ndPyramidal_EdgeDrag    = 7.28268834977848e-05;
    param->HCP_A_2ndPyramidal_ScrewDrag   = 1.05154866364405e-04;

    // From polyxtal literature, values are assumed 20 times basal edge

    param->HCP_CpA_Prismatic_EdgeDrag     = 2.36868775797546e-03;
    param->HCP_CpA_Prismatic_ScrewDrag    = 2.36868775797546e-03;

    param->HCP_CpA_1stPyramidal_EdgeDrag  = 2.36868775797546e-03;
    param->HCP_CpA_1stPyramidal_ScrewDrag = 2.36868775797546e-03;

    param->HCP_CpA_2ndPyramidal_EdgeDrag  = 2.36868775797546e-03;
    param->HCP_CpA_2ndPyramidal_ScrewDrag = 2.36868775797546e-03;

    param->HCP_C_Prismatic_EdgeDrag       = 2.36868775797546e-03;
    param->HCP_C_Prismatic_ScrewDrag      = 2.36868775797546e-03;

    param->HCP_Sessile_EdgeDrag           = 23.6868;  // 10000 * pyramidal
    param->HCP_Sessile_ScrewDrag          = 23.6868;  // 10000 * pyramidal

    param->HCP_LineDrag                   = 5.0e-6;   // ~10 times lowest

    // For the "relaxation" mobility laws only
    BindVar(CPList, "MobRelaxX"            , &param->MobRelaxX, V_DBL, 1, VFLAG_NULL); param->MobRelaxX = 1.0;
    BindVar(CPList, "MobRelaxY"            , &param->MobRelaxY, V_DBL, 1, VFLAG_NULL); param->MobRelaxY = 1.0;
    BindVar(CPList, "MobRelaxZ"            , &param->MobRelaxZ, V_DBL, 1, VFLAG_NULL); param->MobRelaxZ = 1.0;
    BindVar(CPList, "MobRelaxScaleByLength", &param->MobRelaxScaleByLength, V_INT, 1, VFLAG_NULL); param->MobRelaxScaleByLength = 1;

    // For the linear and non-linear tantalum mobility law
    BindVar(CPList, "TaDragEdge110"        , &param->TaDragEdge110,    V_DBL, 1, VFLAG_NULL); param->TaDragEdge110= 0.00028532;
    BindVar(CPList, "TaDragEdge112"        , &param->TaDragEdge112,    V_DBL, 1, VFLAG_NULL); param->TaDragEdge112= 0.00028532;
    BindVar(CPList, "TaDragScrew112T"      , &param->TaDragScrew112T,  V_DBL, 1, VFLAG_NULL); param->TaDragScrew112T=0.017119;
    BindVar(CPList, "TaDragScrew112AT"     , &param->TaDragScrew112AT, V_DBL, 1, VFLAG_NULL); param->TaDragScrew112AT=0.11641;
    BindVar(CPList, "TaDragScrew110"       , &param->TaDragScrew110,   V_DBL, 1, VFLAG_NULL); param->TaDragScrew110= 0.034239;
    BindVar(CPList, "TaDragLine"           , &param->TaDragLine,       V_DBL, 1, VFLAG_NULL); param->TaDragLine=2.8532e-06;
    BindVar(CPList, "TaDragClimbEdge"      , &param->TaDragClimbEdge,  V_DBL, 1, VFLAG_NULL); param->TaDragClimbEdge= 342.39;
    BindVar(CPList, "TaDragClimbScrew"     , &param->TaDragClimbScrew, V_DBL, 1, VFLAG_NULL); param->TaDragClimbScrew=34.239;


    // For the non-linear tantalum mobility law : fitting parameters depends on (P,T). Default for (P=0Mbar,T=300K).
    BindVar(CPList, "TaDragParam110K"     , &param->TaDragParam110K, V_DBL, 1, VFLAG_NULL);     param->TaDragParam110K=9000;
    BindVar(CPList, "TaDragParam110B"     , &param->TaDragParam110B, V_DBL, 1, VFLAG_NULL);     param->TaDragParam110B=0.15;
    BindVar(CPList, "TaDragParam110alpha" , &param->TaDragParam110alpha, V_DBL, 1, VFLAG_NULL); param->TaDragParam110alpha=1;

    BindVar(CPList, "TaDragParam112TK"     , &param->TaDragParam112TK, V_DBL, 1, VFLAG_NULL);     param->TaDragParam112TK=2200;
    BindVar(CPList, "TaDragParam112TB"     , &param->TaDragParam112TB, V_DBL, 1, VFLAG_NULL);     param->TaDragParam112TB=0.2;
    BindVar(CPList, "TaDragParam112Talpha" , &param->TaDragParam112Talpha, V_DBL, 1, VFLAG_NULL); param->TaDragParam112Talpha=1.5;

    BindVar(CPList, "TaDragParam112ATK"     , &param->TaDragParam112ATK, V_DBL, 1, VFLAG_NULL);     param->TaDragParam112ATK=10600;
    BindVar(CPList, "TaDragParam112ATB"     , &param->TaDragParam112ATB, V_DBL, 1, VFLAG_NULL);     param->TaDragParam112ATB=0.5;
    BindVar(CPList, "TaDragParam112ATalpha" , &param->TaDragParam112ATalpha, V_DBL, 1, VFLAG_NULL); param->TaDragParam112ATalpha=3;

    // List of burgers vectors/line directions to be considered sessile

    BindVar(CPList, "sessileburgspec"      ,  param->sessileburgspec, V_DBL, 30, VFLAG_NULL);
    BindVar(CPList, "sessilelinespec"      ,  param->sessilelinespec, V_DBL, 30, VFLAG_NULL);

    // A couple values for including (if available) inertial terms in
    // the mobility function.

    BindVar(CPList, "includeInertia"       , &param->includeInertia, V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "massDensity"          , &param->massDensity   , V_DBL, 1, VFLAG_NULL); param->massDensity = -1.0;

    // Flux decomposition

    BindVar(CPList, "Flux decomposition"   , (void *) NULL, V_COMMENT, 1, VFLAG_NULL);
    BindVar(CPList, "bcc_DensityFluxDecomp", &param->bcc_DensityFluxDecomp, V_INT,   1, VFLAG_NULL); param->bcc_DensityFluxDecomp = 2;
    BindVar(CPList, "totstraintensor"      ,  param->totstraintensor      , V_DBL,   6, VFLAG_NULL);
    BindVar(CPList, "elasticStrainTensor"  ,  param->elasticStrainTensor  , V_DBL,   6, VFLAG_NULL);
    BindVar(CPList, "totpStn"              ,  param->totpStn              , V_DBL,   6, VFLAG_NULL);
    BindVar(CPList, "totpSpn"              ,  param->totpSpn              , V_DBL,   6, VFLAG_NULL);
    BindVar(CPList, "BCC_LtotDecomp"       ,  param->BCC_LtotDecomp       , V_DBL,  44, VFLAG_NULL);
    BindVar(CPList, "BCC2_LtotDecomp"      ,  param->BCC2_LtotDecomp      , V_DBL,  80, VFLAG_NULL);
    BindVar(CPList, "FCC_LtotDecomp"       ,  param->FCC_LtotDecomp       , V_DBL,  66, VFLAG_NULL);
    BindVar(CPList, "HCP_LtotDecomp"       ,  param->HCP_LtotDecomp       , V_DBL, 170, VFLAG_NULL);
    BindVar(CPList, "rhombo_Ltot"          ,  param->rhombo_Ltot          , V_DBL,  24, VFLAG_NULL);
    BindVar(CPList, "BCC_fluxtot"          ,  param->BCC_fluxtot          , V_DBL,  28, VFLAG_NULL);
    BindVar(CPList, "BCC2_fluxtot"         ,  param->BCC2_fluxtot         , V_DBL,  52, VFLAG_NULL);
    BindVar(CPList, "FCC_fluxtot"          ,  param->FCC_fluxtot          , V_DBL,  42, VFLAG_NULL);
    BindVar(CPList, "HCP_fluxtot"          ,  param->HCP_fluxtot          , V_DBL,  90, VFLAG_NULL);
    BindVar(CPList, "rhombo_fluxtot"       ,  param->rhombo_fluxtot       , V_DBL,  36, VFLAG_NULL);

    // Total system density: this is informational only and recalculated
    // each time a restart file is written.

    BindVar(CPList, "Total density. Informational only; ignored on input", (void *) NULL, V_COMMENT, 1, VFLAG_NULL);
    BindVar(CPList, "disloDensity"         , &param->disloDensity         , V_DBL, 1, VFLAG_NULL);

    // Velocity statistics

    BindVar(CPList, "Velocity statistics"  , (void *) NULL, V_COMMENT, 1, VFLAG_NULL);
    BindVar(CPList, "vAverage"             , &param->vAverage             , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "vStDev"               , &param->vStDev               , V_DBL, 1, VFLAG_NULL);

    // I/O controls and options

    BindVar(CPList, "I/O controls and parameters", (void *) NULL, V_COMMENT, 1, VFLAG_NULL);
    BindVar(CPList, "dirname"              ,  param->dirname              , V_STRING, 1, VFLAG_NULL);
    BindVar(CPList, "writeBinRestart"      , &param->writeBinRestart      , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "skipIO"               , &param->skipIO               , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "numIOGroups"          , &param->numIOGroups          , V_INT, 1, VFLAG_NULL); param->numIOGroups = 1;

    // segment/arm files

    BindVar(CPList, "armfile"              , &param->armfile              , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "armfilefreq"          , &param->armfilefreq          , V_INT, 1, VFLAG_NULL); param->armfilefreq = 100;
    BindVar(CPList, "armfiledt"            , &param->armfiledt            , V_DBL, 1, VFLAG_NULL); param->armfiledt   = -1.0;
    BindVar(CPList, "armfiletime"          , &param->armfiletime          , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "armfilecounter"       , &param->armfilecounter       , V_INT, 1, VFLAG_NULL);


    // flux decomposition files

    // Note: the flux related parameters below were renamed, but
    //       we're providing aliases for the original parameter names
    //       so they will still be recognized from older restart files.

    BindVar(CPList, "writeFlux"                , &param->writeFlux                , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "writeFluxEdgeDecomp"      , &param->writeFluxEdgeDecomp      , V_INT, 1, VFLAG_NULL); param->writeFluxEdgeDecomp = 1;
    BindVar(CPList, "writeFluxFullDecomp"      , &param->writeFluxFullDecomp      , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "writeFluxFullDecompTotals", &param->writeFluxFullDecompTotals, V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "writeFluxSimpleTotals"    , &param->writeFluxSimpleTotals    , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "writeFluxFreq"            , &param->writeFluxFreq            , V_INT, 1, VFLAG_NULL); param->writeFluxFreq = 100;
    BindVar(CPList, "writeFluxDT"              , &param->writeFluxDT              , V_DBL, 1, VFLAG_NULL); param->writeFluxDT   = -1.0;
    BindVar(CPList, "writeFluxTime"            , &param->writeFluxTime            , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "writeFluxCounter"         , &param->writeFluxCounter         , V_INT, 1, VFLAG_NULL);

    // Aliases

    BindVar(CPList, "fluxfile"                , &param->writeFlux                 , V_INT, 1, VFLAG_ALIAS);
    BindVar(CPList, "fluxfreq"                , &param->writeFluxFreq             , V_INT, 1, VFLAG_ALIAS);
    BindVar(CPList, "fluxdt"                  , &param->writeFluxDT               , V_DBL, 1, VFLAG_ALIAS);
    BindVar(CPList, "fluxtime"                , &param->writeFluxTime             , V_DBL, 1, VFLAG_ALIAS);
    BindVar(CPList, "fluxcounter"             , &param->writeFluxCounter          , V_INT, 1, VFLAG_ALIAS);

    // gnuplot files

    BindVar(CPList, "gnuplot"                 , &param->gnuplot                   , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "gnuplotfreq"             , &param->gnuplotfreq               , V_INT, 1, VFLAG_NULL); param->gnuplotfreq = 100;
    BindVar(CPList, "gnuplotdt"               , &param->gnuplotdt                 , V_DBL, 1, VFLAG_NULL); param->gnuplotdt   = -1.0;
    BindVar(CPList, "gnuplottime"             , &param->gnuplottime               , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "gnuplotcounter"          , &param->gnuplotcounter            , V_INT, 1, VFLAG_NULL);

    // pole figure files

    BindVar(CPList, "polefigfile"             , &param->polefigfile               , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "polefigfreq"             , &param->polefigfreq               , V_INT, 1, VFLAG_NULL); param->polefigfreq = 100;
    BindVar(CPList, "polefigdt"               , &param->polefigdt                 , V_DBL, 1, VFLAG_NULL); param->polefigdt   = -1.0;
    BindVar(CPList, "polefigtime"             , &param->polefigtime               , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "polefilecounter"         , &param->polefigcounter            , V_INT, 1, VFLAG_NULL);

    // povray files

    BindVar(CPList, "povray"                  , &param->povray                    , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "povrayfreq"              , &param->povrayfreq                , V_INT, 1, VFLAG_NULL); param->povrayfreq = 100;
    BindVar(CPList, "povraydt"                , &param->povraydt                  , V_DBL, 1, VFLAG_NULL); param->povraydt   = -1.0;
    BindVar(CPList, "povraytime"              , &param->povraytime                , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "povraycounter"           , &param->povraycounter             , V_INT, 1, VFLAG_NULL);

    // postscript files

    BindVar(CPList, "psfile"                  , &param->psfile                    , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "psfilefreq"              , &param->psfilefreq                , V_INT, 1, VFLAG_NULL); param->psfilefreq = 100;
    BindVar(CPList, "psfiledt"                , &param->psfiledt                  , V_DBL, 1, VFLAG_NULL); param->psfiledt   = -1.0;
    BindVar(CPList, "psfiletime"              , &param->psfiletime                , V_DBL, 1, VFLAG_NULL);

    // restart files

    BindVar(CPList, "savecn"                  , &param->savecn                    , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "savecnfreq"              , &param->savecnfreq                , V_INT, 1, VFLAG_NULL); param->savecnfreq = 100;
    BindVar(CPList, "savecndt"                , &param->savecndt                  , V_DBL, 1, VFLAG_NULL); param->savecndt   = -1.0;
    BindVar(CPList, "savecntime"              , &param->savecntime                , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "savecncounter"           , &param->savecncounter             , V_INT, 1, VFLAG_NULL);

    // delta time debug output

    BindVar(CPList, "savedt"                  , &param->savedt                    , V_INT, 1, VFLAG_NULL); param->savedt     = 0;
    BindVar(CPList, "savedtfreq"              , &param->savedtfreq                , V_INT, 1, VFLAG_NULL); param->savedtfreq = 100;

    // nodal arm diagnostic output

    BindVar(CPList, "save_narms"              , &param->save_narms                , V_INT, 1, VFLAG_NULL); param->save_narms      = 0;
    BindVar(CPList, "save_narms_freq"         , &param->save_narms_freq           , V_INT, 1, VFLAG_NULL); param->save_narms_freq = 100;

    // properties files

    BindVar(CPList, "saveprop"                , &param->saveprop                  , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "savepropfreq"            , &param->savepropfreq              , V_INT, 1, VFLAG_NULL); param->savepropfreq = 100;
    BindVar(CPList, "savepropdt"              , &param->savepropdt                , V_DBL, 1, VFLAG_NULL); param->savepropdt   = -1.0;
    BindVar(CPList, "saveproptime"            , &param->saveproptime              , V_DBL, 1, VFLAG_NULL);

    // timer files

    BindVar(CPList, "savetimers"              , &param->savetimers                , V_INT, 1, VFLAG_NULL); param->savetimers     = 0;
    BindVar(CPList, "savetimersfreq"          , &param->savetimersfreq            , V_INT, 1, VFLAG_NULL); param->savetimersfreq = 100;
    BindVar(CPList, "savetimersdt"            , &param->savetimersdt              , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "savetimerstime"          , &param->savetimerstime            , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "savetimerscounter"       , &param->savetimerscounter         , V_INT, 1, VFLAG_NULL); param->savetimerscounter = 0;

    // tecplot files

    BindVar(CPList, "tecplot"                 , &param->tecplot                   , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "tecplotfreq"             , &param->tecplotfreq               , V_INT, 1, VFLAG_NULL); param->tecplotfreq = 100;
    BindVar(CPList, "tecplotdt"               , &param->tecplotdt                 , V_DBL, 1, VFLAG_NULL); param->tecplotdt   = -1.0;
    BindVar(CPList, "tecplottime"             , &param->tecplottime               , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "tecplotcounter"          , &param->tecplotcounter            , V_INT, 1, VFLAG_NULL);

    // nodal velocity data

    BindVar(CPList, "velfile"                 , &param->velfile                   , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "velfilefreq"             , &param->velfilefreq               , V_INT, 1, VFLAG_NULL); param->velfilefreq = 100;
    BindVar(CPList, "velfiledt"               , &param->velfiledt                 , V_DBL, 1, VFLAG_NULL); param->velfiledt   = -1.0;
    BindVar(CPList, "velfiletime"             , &param->velfiletime               , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "velfilecounter"          , &param->velfilecounter            , V_INT, 1, VFLAG_NULL);

    // nodal force data

    BindVar(CPList, "writeForce"              , &param->writeForce                , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "writeForceFreq"          , &param->writeForceFreq            , V_INT, 1, VFLAG_NULL); param->writeForceFreq = 100;
    BindVar(CPList, "writeForceDT"            , &param->writeForceDT              , V_DBL, 1, VFLAG_NULL); param->writeForceDT   = -1.0;
    BindVar(CPList, "writeForceTime"          , &param->writeForceTime            , V_DBL, 1, VFLAG_NULL);
    BindVar(CPList, "writeForceCounter"       , &param->writeForceCounter         , V_INT, 1, VFLAG_NULL);

    // Parameters related to creation of VisIt output files

    BindVar(CPList, "writeVisit"              , &param->writeVisit                , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "writeVisitFreq"          , &param->writeVisitFreq            , V_INT, 1, VFLAG_NULL); param->writeVisitFreq     = 100;
    BindVar(CPList, "writeVisitCounter"       , &param->writeVisitCounter         , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "writeVisitSegments"      , &param->writeVisitSegments        , V_INT, 1, VFLAG_NULL); param->writeVisitSegments = 1;
    BindVar(CPList, "writeVisitSegmentsAsText", &param->writeVisitSegmentsAsText  , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "writeVisitNodes"         , &param->writeVisitNodes           , V_INT, 1, VFLAG_NULL); param->writeVisitNodes    = 1;
    BindVar(CPList, "writeVisitNodesAsText"   , &param->writeVisitNodesAsText     , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "writeVisitBurgID"        , &param->writeVisitBurgID          , V_INT, 1, VFLAG_NULL); param->writeVisitBurgID   = 1;
    BindVar(CPList, "writeVisitForceVector"   , &param->writeVisitForceVector     , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "writeVisitVelocityVector", &param->writeVisitVelocityVector  , V_INT, 1, VFLAG_NULL);
    BindVar(CPList, "writeVisitDT"            , &param->writeVisitDT              , V_DBL, 1, VFLAG_NULL); param->writeVisitDT       = -1.0;
    BindVar(CPList, "writeVisitTime"          , &param->writeVisitTime            , V_DBL, 1, VFLAG_NULL);

    // X-window display

    BindVar(CPList, "winDefaultsFile", param->winDefaultsFile, V_STRING, 1, VFLAG_NULL);
    strcpy(param->winDefaultsFile, "inputs/paradis.xdefaults");

    // 3D-dislocation density data

    BindVar(CPList, "savedensityspec", param->savedensityspec, V_INT, 3, VFLAG_NULL);

    // Miscellaneous parameters

    BindVar(CPList, "Miscellaneous parameters", (void *) NULL, V_COMMENT, 1, VFLAG_NULL);
    BindVar(CPList, "enforceGlidePlanes"    , &param->enforceGlidePlanes   , V_INT   , 1, VFLAG_NULL);
    BindVar(CPList, "enableCrossSlip"       , &param->enableCrossSlip      , V_INT   , 1, VFLAG_NULL); param->enableCrossSlip       = -1;
    BindVar(CPList, "crossSlipBCCIndex"     , &param->crossSlipBCCIndex    , V_INT   , 1, VFLAG_NULL); param->crossSlipBCCIndex     =  1;
    BindVar(CPList, "HCP_CA_SplitFrequency" , &param->HCP_CA_SplitFrequency, V_INT   , 1, VFLAG_NULL); param->HCP_CA_SplitFrequency =  0;
    BindVar(CPList, "splinterSegFrequency"  , &param->splinterSegFrequency , V_INT   , 1, VFLAG_NULL); param->splinterSegFrequency  =  1;
    BindVar(CPList, "TensionFactor"         , &param->TensionFactor        , V_DBL   , 1, VFLAG_NULL); param->TensionFactor         =  1;
    BindVar(CPList, "elasticinteraction"    , &param->elasticinteraction   , V_INT   , 1, VFLAG_NULL); param->elasticinteraction    =  1;
}

void DataParamInit(Param_t *param, ParamList_t *DPList)
{
    // Note: Parameters need only be initialized if their
    // default values are non-zero.

    DPList->paramCnt = 0;
    DPList->varList  = (VarData_t *) NULL;

    BindVar(DPList, "dataFileVersion"   , &param->dataFileVersion   , V_INT, 1, VFLAG_NULL); param->dataFileVersion = NODEDATA_FILE_VERSION;
    BindVar(DPList, "numFileSegments"   , &param->numFileSegments   , V_INT, 1, VFLAG_NULL); param->numFileSegments = 1;
    BindVar(DPList, "minCoordinates"    ,  param->minCoordinates    , V_DBL, 3, VFLAG_NULL);
    BindVar(DPList, "maxCoordinates"    ,  param->maxCoordinates    , V_DBL, 3, VFLAG_NULL);
    BindVar(DPList, "nodeCount"         , &param->nodeCount         , V_INT, 1, VFLAG_NULL);
    BindVar(DPList, "dataDecompType"    , &param->dataDecompType    , V_INT, 1, VFLAG_NULL); param->dataDecompType        = 2;
    BindVar(DPList, "dataDecompGeometry",  param->dataDecompGeometry, V_INT, 3, VFLAG_NULL); param->dataDecompGeometry[0] = 1;
                                                                                             param->dataDecompGeometry[1] = 1;
                                                                                             param->dataDecompGeometry[2] = 1;
}

/*------------------------------------------------------------------------
 *
 *      Functions:    Print
 *      Description:  Prints a human-readable summary of the parameter data structure.
 *
 *      Arguments:
 *          p         Pointer to the parameter structure
 *
 *          name      simple title message added to the output
 *
 *----------------------------------------------------------------------*/

void Print_Param (Param_t *p, const char *msg)
{

   if (p == (Param_t *) NULL) return;

   printf("\nParaDiS Param_t data structure status %s\n", msg ? msg : "");

   printf("   nXdoms                               : %d\n"      , p->nXdoms );
   printf("   nYdoms                               : %d\n"      , p->nYdoms );
   printf("   nZdoms                               : %d\n\n"    , p->nZdoms );

#if defined _BGP || defined _BGQ
   printf("   taskMappingMode                      : %d\n\n"    , p->taskMappingMode );

#ifdef _BGQ
   printf("   optimalTaskGeometry                  : [%d %d %d]\n\n" , p->optimalTaskGeometry[0],
                                                                       p->optimalTaskGeometry[1],
                                                                       p->optimalTaskGeometry[2] );
#endif
#endif

   printf("   assignNewNodeTags                    : %d\n"      , p->assignNewNodeTags );
   printf("   nXcells                              : %d\n"      , p->nXcells );
   printf("   nYcells                              : %d\n"      , p->nYcells );
   printf("   nZcells                              : %d\n\n"    , p->nZcells );

   printf("   iCellNatMin                          : %d\n"      , p->iCellNatMin );
   printf("   iCellNatMax                          : %d\n"      , p->iCellNatMax );
   printf("   jCellNatMin                          : %d\n"      , p->jCellNatMin );
   printf("   jCellNatMax                          : %d\n"      , p->jCellNatMax );
   printf("   kCellNatMin                          : %d\n"      , p->kCellNatMin );
   printf("   kCellNatMax                          : %d\n\n"    , p->kCellNatMax );

   printf("   xBoundType                           : %d (%s)\n"      , (int) p->xBoundType, ( p->xBoundType==Periodic ? "Periodic" : p->xBoundType==Free ? "Free" : "Reflecting" ) );
   printf("   yBoundType                           : %d (%s)\n"      , (int) p->yBoundType, ( p->yBoundType==Periodic ? "Periodic" : p->yBoundType==Free ? "Free" : "Reflecting" ) );
   printf("   zBoundType                           : %d (%s)\n\n"    , (int) p->zBoundType, ( p->zBoundType==Periodic ? "Periodic" : p->zBoundType==Free ? "Free" : "Reflecting" ) );

   printf("   xBoundMin                            : [%+lf %+lf]\n"    , p->xBoundMin, p->xBoundMax );
   printf("   yBoundMin                            : [%+lf %+lf]\n"    , p->yBoundMin, p->yBoundMax );
   printf("   zBoundMin                            : [%+lf %+lf]\n\n"  , p->zBoundMin, p->zBoundMax );

   printf("   minSideX                             : [%+lf %+lf]\n"    , p->minSideX, p->maxSideX );
   printf("   minSideY                             : [%+lf %+lf]\n"    , p->minSideY, p->maxSideY );
   printf("   minSideZ                             : [%+lf %+lf]\n\n"  , p->minSideZ, p->maxSideZ );

   printf("   decompType                           : %d\n"       , p->decompType );
   printf("   DLBfreq                              : %d\n"       , p->DLBfreq );
   printf("   numDLBCycles                         : %d\n\n"     , p->numDLBCycles );

   printf("   cycleStart                           : %d\n"       , p->cycleStart );
   printf("   maxstep                              : %d\n"       , p->maxstep );
   printf("   timeStart                            : %+lf\n"     , p->timeStart );
   printf("   timeNow                              : %+lf\n"     , p->timeNow );
   printf("   timeEnd                              : %+lf\n\n"   , p->timeEnd );

   printf("   timestepIntegrator                   : %s\n"       , p->timestepIntegrator );
   printf("   trapezoidMaxIterations               : %d\n\n"     , p->trapezoidMaxIterations );

   printf("   deltaTT                              : %+lf\n"     , p->deltaTT );
   printf("   realdt                               : %+lf\n"     , p->realdt );
   printf("   nextDT                               : %+lf\n"     , p->nextDT );
   printf("   maxDT                                : %+lf\n"     , p->maxDT );
   printf("   dtIncrementFact                      : %+lf\n"     , p->dtIncrementFact );
   printf("   dtDecrementFact                      : %+lf\n"     , p->dtDecrementFact );
   printf("   dtExponent                           : %+lf\n"     , p->dtExponent );
   printf("   dtVariableAdjustment                 : %d\n"       , p->dtVariableAdjustment );
   printf("   rTol                                 : %+lf\n"     , p->rTol );
   printf("   rmax                                 : %+lf\n\n"   , p->rmax );

   printf("   KINSOL_UseNewton                     : %d\n"       , p->KINSOL_UseNewton );
   printf("   KINSOL_MaxIterations                 : %d\n"       , p->KINSOL_MaxIterations );
   printf("   KINSOL_MaxLinIterations              : %d\n"       , p->KINSOL_MaxLinIterations );
   printf("   KINSOL_NumPriorResiduals             : %d\n"       , p->KINSOL_NumPriorResiduals );
   printf("   KINSOL_PrintLevel                    : %d\n"       , p->KINSOL_PrintLevel );
   printf("   KINSOL_LogFile                       : %s\n"       , p->KINSOL_LogFile );
   printf("   KINSOL_EtaVal                        : %+lf\n\n"   , p->KINSOL_EtaVal );

   printf("   ARKODE_IRKtable                      : %d\n"       , p->ARKODE_IRKtable );
   printf("   ARKODE_FPsolver                      : %d\n"       , p->ARKODE_FPsolver );
   printf("   ARKODE_MaxNonLinIters                : %d\n"       , p->ARKODE_MaxNonLinIters );
   printf("   ARKODE_MaxLinIters                   : %d\n"       , p->ARKODE_MaxLinIters );
   printf("   ARKODE_FPaccel                       : %d\n"       , p->ARKODE_FPaccel );
   printf("   ARKODE_PredictorMethod               : %d\n"       , p->ARKODE_PredictorMethod );
   printf("   ARKODE_AdaptivityMethod              : %d\n"       , p->ARKODE_AdaptivityMethod );
   printf("   ARKODE_NonLinCoef                    : %+lf\n"     , p->ARKODE_NonLinCoef );
   printf("   ARKODE_LinCoef                       : %+lf\n\n"   , p->ARKODE_LinCoef );

   printf("   minSeg                               : %+lf\n"     , p->minSeg );
   printf("   maxSeg                               : %+lf\n"     , p->maxSeg );
   printf("   rann                                 : %+lf\n\n"   , p->rann );

   printf("   remeshRule                           : %d\n"       , p->remeshRule );
   printf("   collisionMethod                      : %d\n"       , p->collisionMethod );
   printf("   remeshAreaMax                        : %+lf\n"     , p->remeshAreaMax );
   printf("   remeshAreaMin                        : %+lf\n"     , p->remeshAreaMin );
   printf("   splitMultiNodeFreq                   : %d\n"       , p->splitMultiNodeFreq );
   printf("   splitMultiNodeAlpha                  : %+lf\n\n"   , p->splitMultiNodeAlpha );

   printf("   useParallelSplitMultiNode            : %d\n\n"     , p->useParallelSplitMultiNode );

   printf("   useUniformJunctionLen                : %d\n\n"     , p->useUniformJunctionLen );

   printf("   fmEnabled                            : %d\n"       , p->fmEnabled );
   printf("   fmNumLayers                          : %d\n"       , p->fmNumLayers );
   printf("   fmMPOrder                            : %d\n"       , p->fmMPOrder );
   printf("   fmExpansionOrder                        : %d\n"       , p->fmExpansionOrder );
   printf("   fmNumPoints                          : %d\n"       , p->fmNumPoints );
   printf("   fmCorrectionTbl                      : %s\n"       , p->fmCorrectionTbl );
   printf("   fmFullCorrectionTbl                  : %d\n"       , p->fmFullCorrectionTbl );
   printf("   fmEnableAnisotropy                   : %d\n\n"     , p->fmEnableAnisotropy );
#if defined (ANISOTROPIC) && defined (TAYLORFMM)
   printf("   fmAnisoDerivTbl                      : %s\n\n"     , p->fmAnisoDerivTbl );
#endif

#ifdef ESHELBY
   printf("   eshelbyfmEnabled                     : %d\n"       , p->eshelbyfmEnabled );
   printf("   eshelbyfmMPOrder                     : %d\n"       , p->eshelbyfmMPOrder );
   printf("   eshelbyfmTaylorOrder                 : %d\n\n"     , p->eshelbyfmTaylorOrder );
#endif

   printf("   Rijmfile                             : %s\n"       , p->Rijmfile );
   printf("   RijmPBCfile                          : %s\n\n"     , p->RijmPBCfile );

   printf("   TempK                                : %+lf\n"     , p->TempK );
   printf("   pressure                             : %+lf\n"     , p->pressure );
   printf("   loadType                             : %d\n"      , p->loadType );
   printf("   appliedStress                        : [%+lf %+lf %+lf %+lf %+lf %+lf]\n"  , p->appliedStress[0], p->appliedStress[1], p->appliedStress[2]
                                                                                         , p->appliedStress[3], p->appliedStress[4], p->appliedStress[5] );
   printf("   eRate                                : %+lf\n"     , p->eRate );
   printf("   edotdir                              : [%+lf %+lf %+lf]\n"     , p->edotdir[0], p->edotdir[1], p->edotdir[2] );
   printf("   cTimeOld                             : %+lf\n"     , p->cTimeOld );
   printf("   netCyclicStrain                      : %+lf\n"     , p->netCyclicStrain );
   printf("   dCyclicStrain                        : %+lf\n"     , p->dCyclicStrain );
   printf("   numLoadCycle                         : %d\n"       , p->numLoadCycle );
   printf("   eAmp                                 : %+lf\n\n"   , p->eAmp );

   printf("   useLabFrame                          : %d\n\n"     , p->useLabFrame );

   printf("   rotMatrix                            : [ %+lf %+lf %+lf ]\n"  , p->rotMatrix[0][0], p->rotMatrix[0][1], p->rotMatrix[0][2] );
   printf("                                        : [ %+lf %+lf %+lf ]\n"  , p->rotMatrix[1][0], p->rotMatrix[1][1], p->rotMatrix[1][2] );
   printf("                                        : [ %+lf %+lf %+lf ]\n\n", p->rotMatrix[2][0], p->rotMatrix[2][1], p->rotMatrix[2][2] );

   printf("   rotMatrixInverse                     : [ %+lf %+lf %+lf ]\n"  , p->rotMatrixInverse[0][0], p->rotMatrixInverse[0][1], p->rotMatrixInverse[0][2] );
   printf("                                        : [ %+lf %+lf %+lf ]\n"  , p->rotMatrixInverse[1][0], p->rotMatrixInverse[1][1], p->rotMatrixInverse[1][2] );
   printf("                                        : [ %+lf %+lf %+lf ]\n\n", p->rotMatrixInverse[2][0], p->rotMatrixInverse[2][1], p->rotMatrixInverse[2][2] );

   //printf("   rotMatrixInverse                     : %+lf\n\n"   , p->rotMatrixInverse [3][3] );

#ifdef ESHELBY
   printf("   inclusionFile                        : %s\n\n"     , p->inclusionFile );
   printf("   enableInclusions                     : %d\n"       , p->enableInclusions );
#endif

   printf("   ecFile                               : %s\n\n"     , p->ecFile );

   printf("   meltTemp                             : %+lf\n"     , p->meltTemp );
   printf("   meltTempFile                         : %s\n\n"     , p->meltTempFile );

   printf("   MobFileNL                            : %s\n"       , p->MobFileNL );
   printf("   mobilityLaw                          : %s\n"       , p->mobilityLaw );
   printf("   mobilityType                         : %d\n\n"     , p->mobilityType );

   //real8  elasticConstantMatrix [6][6] );
   //real8  elasticConstantMatrixInv [6][6] );

#if defined ANISOTROPIC
   //real8  C11, C12, C13, C33, C44, Cpr );

   printf("  anisoHarmonicsNumTermsBase            : %d\n"       , p->anisoHarmonicsNumTermsBase );

   printf("  anisoNumThetaPoints                   : %d\n"       , p-> anisoNumThetaPoints );
   printf("  anisoNumPhiPoints                     : %d\n"       , p-> anisoNumPhiPoints );
#endif

   //printf("   materialTypeName                     : %s\n\n"    , p->materialTypeName [64] );

   printf("   materialType                         : %d\n\n"     , p->materialType );

   printf("   vacancyConc                          : %+lf\n\n"   , p->vacancyConc );

   printf("   vacancyConcEquilibrium               : %+lf\n\n"   , p->vacancyConcEquilibrium );

   printf("   shearModulus                         : %+lf\n"     , p->shearModulus );
   printf("   pois                                 : %+lf\n"     , p->pois );
   printf("   burgMag                              : %+lf\n\n"   , p->burgMag );

   printf("   cOVERa                               : %+lf\n\n"   , p->cOVERa );

   printf("   YoungsModulus                        : %+lf\n\n"   , p->YoungsModulus );

   printf("   rc                                   : %+lf\n\n"   , p->rc );

   printf("   Ecore                                : %+lf\n\n"   , p->Ecore );

   printf("   HCPEcoreA                            : %+lf\n"     , p->HCPEcoreA );
   printf("   HCPEcoreC                            : %+lf\n"     , p->HCPEcoreC );
   printf("   HCPEcoreCpA                          : %+lf\n\n"   , p->HCPEcoreCpA );

   printf("   HCP_CA_SplitFrequency                : %d\n\n"     , p->HCP_CA_SplitFrequency );

   printf("   splinterSegFrequency                 : %d\n\n"     , p->splinterSegFrequency );

   printf("   enforceGlidePlanes                   : %d\n\n"     , p->enforceGlidePlanes );

   printf("   allowFuzzyGlidePlanes                : %d\n"       , p->allowFuzzyGlidePlanes );
   printf("   enableCrossSlip                      : %d\n\n"     , p->enableCrossSlip );

   printf("   mobilityIndex                        : %d\n"        ,                 p->mobilityIndex    );
   printf("   mobilityInitFunc                     : 0x%012lx\n\n", (unsigned long) p->mobilityInitFunc );

   printf("   MobScrew                             : %+lf\n"     , p->MobScrew );
   printf("   MobEdge                              : %+lf\n"     , p->MobEdge  );
   printf("   MobClimb                             : %+lf\n"     , p->MobClimb );
   printf("   MobGlide                             : %+lf\n"     , p->MobGlide );
   printf("   MobLine                              : %+lf\n\n"   , p->MobLine  );

   printf("   MobExpP                              : %+lf\n"     , p->MobExpP    );
   printf("   MobExpQ                              : %+lf\n"     , p->MobExpQ    );
   printf("   MobPeierls                           : %+lf\n"     , p->MobPeierls );
   printf("   MobDeltaH0                           : %+lf\n"     , p->MobDeltaH0 );
   printf("   MobAlpha                             : %+lf\n\n"   , p->MobAlpha   );

   printf("   MobCoeffA                            : %+lf\n"     , p->MobCoeffA  );
   printf("   MobCoeffC                            : %+lf\n"     , p->MobCoeffC  );
   printf("   MobCoeffC0                           : %+lf\n"     , p->MobCoeffC0 );
   printf("   MobCoeffD                            : %+lf\n\n"   , p->MobCoeffD  );

   printf("   MobG0                                : %+lf\n\n"   , p->MobG0 );

   printf("   MobNumTrenches                       : %d\n\n"     , p->MobNumTrenches );

   //printf("   MobTrenchAngle                       : %+lf\n\n"   , p->MobTrenchAngle [MAX_MOB_TRENCHES] );

   //printf("   MobTrenchWidth                       : %+lf\n\n"   , p->MobTrenchWidth [MAX_MOB_TRENCHES] );

   //printf("   MobTrench                            : %+lf\n\n"   , p->MobTrench [MAX_MOB_TRENCHES] );

   printf("   MobFacetedUseCusps                   : %d\n\n"     , p->MobFacetedUseCusps                );

   printf("   HCP_A_Basal_ScrewDrag                : %+lf\n"     , p->HCP_A_Basal_ScrewDrag             );
   printf("   HCP_A_Prismatic_ScrewDrag            : %+lf\n"     , p->HCP_A_Prismatic_ScrewDrag         );
   printf("   HCP_A_1stPyramidal_ScrewDrag         : %+lf\n"     , p->HCP_A_1stPyramidal_ScrewDrag      );
   printf("   HCP_A_2ndPyramidal_ScrewDrag         : %+lf\n"     , p->HCP_A_2ndPyramidal_ScrewDrag      );
   printf("   HCP_C_Prismatic_ScrewDrag            : %+lf\n"     , p->HCP_C_Prismatic_ScrewDrag         );
   printf("   HCP_CpA_Prismatic_ScrewDrag          : %+lf\n"     , p->HCP_CpA_Prismatic_ScrewDrag       );
   printf("   HCP_CpA_1stPyramidal_ScrewDrag       : %+lf\n"     , p->HCP_CpA_1stPyramidal_ScrewDrag    );
   printf("   HCP_CpA_2ndPyramidal_ScrewDrag       : %+lf\n"     , p->HCP_CpA_2ndPyramidal_ScrewDrag    );
   printf("   HCP_Sessile_ScrewDrag                : %+lf\n\n"   , p->HCP_Sessile_ScrewDrag             );

   printf("   HCP_A_Basal_EdgeDrag                 : %+lf\n"     , p->HCP_A_Basal_EdgeDrag              );
   printf("   HCP_A_Prismatic_EdgeDrag             : %+lf\n"     , p->HCP_A_Prismatic_EdgeDrag          );
   printf("   HCP_A_1stPyramidal_EdgeDrag          : %+lf\n"     , p->HCP_A_1stPyramidal_EdgeDrag       );
   printf("   HCP_A_2ndPyramidal_EdgeDrag          : %+lf\n"     , p->HCP_A_2ndPyramidal_EdgeDrag       );
   printf("   HCP_C_Prismatic_EdgeDrag             : %+lf\n"     , p->HCP_C_Prismatic_EdgeDrag          );
   printf("   HCP_CpA_Prismatic_EdgeDrag           : %+lf\n"     , p->HCP_CpA_Prismatic_EdgeDrag        );
   printf("   HCP_CpA_1stPyramidal_EdgeDrag        : %+lf\n"     , p->HCP_CpA_1stPyramidal_EdgeDrag     );
   printf("   HCP_CpA_2ndPyramidal_EdgeDrag        : %+lf\n"     , p->HCP_CpA_2ndPyramidal_EdgeDrag     );
   printf("   HCP_Sessile_EdgeDrag                 : %+lf\n\n"   , p->HCP_Sessile_EdgeDrag              );

   printf("   HCP_LineDrag                         : %+lf\n\n"   , p->HCP_LineDrag );

   printf("   MobRelaxX                            : %+lf\n"     , p->MobRelaxX );
   printf("   MobRelaxY                            : %+lf\n"     , p->MobRelaxY );
   printf("   MobRelaxZ                            : %+lf\n\n"   , p->MobRelaxZ );

   printf("   MobRelaxScaleByLength                : %d\n\n"    , p->MobRelaxScaleByLength );

   //printf("   sessileburgspec                      : %+lf\n"     , p->sessileburgspec [30] );
   //printf("   sessilelinespec                      : %+lf\n\n"   , p->sessilelinespec [30] );

   printf("   includeInertia                       : %d\n"       , p->includeInertia );
   printf("   massDensity                          : %+lf\n\n"   , p->massDensity );

   printf("   vAverage                             : %+lf\n"     , p->vAverage );
   printf("   vStDev                               : %+lf\n\n"   , p->vStDev );

   printf("   dirname                              : %s\n\n"     , p->dirname );

   printf("   writeBinRestart                      : %d\n\n"     , p->writeBinRestart );

   printf("   doBinRead                            : %d\n\n"     , p->doBinRead );

   printf("   numIOGroups                          : %d\n\n"     , p->numIOGroups );

   printf("   skipIO                               : %d\n\n"     , p->skipIO );

   printf("   armfile                              : %d %3d %d dt=[%+lf %+lf]\n\n" , p->armfile    , p->armfilefreq   , p->armfilecounter   , p->armfiledt   , p->armfiletime    );
   printf("   writeFlux                            : %d %3d %d dt=[%+lf %+lf]\n"   , p->writeFlux  , p->writeFluxFreq , p->writeFluxCounter , p->writeFluxDT , p->writeFluxTime  );
   printf("   writeFluxEdgeDecomp                  : %d\n"       , p->writeFluxEdgeDecomp );
   printf("   writeFluxFullDecomp                  : %d\n"       , p->writeFluxFullDecomp );
   printf("   writeFluxFullDecompTotals            : %d\n"       , p->writeFluxFullDecompTotals );
   printf("   writeFluxSimpleTotals                : %d\n\n"     , p->writeFluxSimpleTotals );

   printf("   gnuplot                              : %d %3d %d dt=[%+lf %+lf]\n"   , p->gnuplot    , p->gnuplotfreq   , p->gnuplotcounter   , p->gnuplotdt   , p->gnuplottime    );
   printf("   polefigfile                          : %d %3d %d dt=[%+lf %+lf]\n"   , p->polefigfile, p->polefigfreq   , p->polefigcounter   , p->polefigdt   , p->polefigtime    );
   printf("   povray                               : %d %3d %d dt=[%+lf %+lf]\n"   , p->povray     , p->povrayfreq    , p->povraycounter    , p->povraydt    , p->povraytime     );
   printf("   psfile                               : %d %3d   dt=[%+lf %+lf]\n"    , p->psfile     , p->psfilefreq                          , p->psfiledt    , p->psfiletime     );
   printf("   savecn                               : %d %3d %d dt=[%+lf %+lf]\n"   , p->savecn     , p->savecnfreq    , p->savecncounter    , p->savecndt    , p->savecntime     );
   printf("   saveprop                             : %d %3d   dt=[%+lf %+lf]\n"    , p->saveprop   , p->savepropfreq                        , p->savepropdt  , p->saveproptime   );
   printf("   savetimers                           : %d %3d %d dt=[%+lf %+lf]\n"   , p->savetimers , p->savetimersfreq, p->savetimerscounter, p->savetimersdt, p->savetimerstime );
   printf("   tecplot                              : %d %3d %d dt=[%+lf %+lf]\n"   , p->tecplot    , p->tecplotfreq   , p->tecplotcounter   , p->tecplotdt   , p->tecplottime    );
   printf("   velfile                              : %d %3d %d dt=[%+lf %+lf]\n"   , p->velfile    , p->velfilefreq   , p->velfilecounter   , p->velfiledt   , p->velfiletime    );
   printf("   writeForce                           : %d %3d %d dt=[%+lf %+lf]\n\n" , p->writeForce , p->writeForceFreq, p->writeForceCounter, p->writeForceDT, p->writeForceTime );

   printf("   savedensityspec                      : %d %d %d\n\n", p->savedensityspec[0], p->savedensityspec[1], p->savedensityspec[2] );

   printf("   writeVisit                           : %d %3d %d dt=[%+lf %+lf]\n"  , p->writeVisit , p->writeVisitFreq, p->writeVisitCounter, p->writeVisitDT, p->writeVisitTime );
   printf("   writeVisitSegments                   : %d %d\n"     , p->writeVisitSegments, p->writeVisitSegmentsAsText );
   printf("   writeVisitNodes                      : %d %d\n"     , p->writeVisitNodes, p->writeVisitNodesAsText );
   printf("   writeVisitBurgID                     : %d %d\n"     , p->writeVisitBurgID, p->writeVisitForceVector );
   printf("   writeVisitVelocityVector             : %d\n\n"      , p->writeVisitVelocityVector );

   printf("   winDefaultsFile                      : %s\n\n"      , p->winDefaultsFile );

   printf("   Lx Ly Lz                             : [ %+lf %+lf %+lf ]\n"     , p->Lx   , p->Ly   , p->Lz    );

   printf("   springConst                          : %+lf\n\n"    , p->springConst );

   printf("   numBurgGroups                        : %d\n\n"      , p->numBurgGroups );

   printf("   partialDisloDensity                  : 0x%012lx\n"  , (unsigned long)  p->partialDisloDensity );
   printf("   disloDensity                         : %+lf\n\n"    , p->disloDensity );

   printf("   delSegLength                         : %+lf\n\n"    , p->delSegLength );

#ifdef _ARLFEM
   printf("   fem_delSegLength                     : %+lf\n\n"    , p->fem_delSegLength );
#endif

   //printf("   densityChange                        : %+lf\n\n"   , p->densityChange [14] );

   printf("   TensionFactor                        : %+lf\n"      , p->TensionFactor );
   printf("   elasticinteraction                   : %d\n\n"      , p->elasticinteraction );

   //printf("   delpStrain                           : %+lf\n"     , p->delpStrain [6],delSig [6],totpStn [6] );
   //printf("   delpSpin                             : %+lf\n\n"   , p->delpSpin [6],totpSpn [6] );

   //printf("   totstraintensor                      : %+lf\n"     , p->totstraintensor [6] );
   //printf("   elasticStrainTensor                  : %+lf\n"     , p->elasticStrainTensor [6] );
   //printf("   totedgepStrain                       : %+lf\n"     , p->totedgepStrain [6], totscrewpStrain [6] );
   //printf("   dedgepStrain                         : %+lf\n\n"   , p->dedgepStrain [6],   dscrewpStrain [6] );

   //printf("   BCC_LtotDecomp                       : %+lf\n"     , p->BCC_LtotDecomp [4][11], BCC_dLtotDecomp [4][11] );
   //printf("   BCC_fluxtot                          : %+lf\n\n"   , p->BCC_fluxtot [4][7],     BCC_dfluxtot [4][7] );

   //printf("   rhombo_Ltot                          : %+lf\n"     , p->rhombo_Ltot [4][3][2],  rhombo_fluxtot [4][3][3] );
   //printf("   rhombo_dLtot                         : %+lf\n\n"   , p->rhombo_dLtot [4][3][2], rhombo_dfluxtot [4][3][3] );

   //printf("   FCC_LtotDecomp                       : %+lf\n"     , p->FCC_LtotDecomp [6][11], FCC_dLtotDecomp [6][11] );
   //printf("   FCC_fluxtot                          : %+lf\n\n"   , p->FCC_fluxtot [6][7],     FCC_dfluxtot [6][7] );

   //printf("   HCP_LtotDecomp                       : %+lf\n"     , p->HCP_LtotDecomp [10][17], HCP_dLtotDecomp [10][17] );
   //printf("   HCP_fluxtot                          : %+lf\n\n"   , p->HCP_fluxtot [10][9],     HCP_dfluxtot [10][9] );

   printf("   imgstrgrid                           : [%d %d %d %d %d %d]\n\n"    , p->imgstrgrid[0], p->imgstrgrid[1], p->imgstrgrid[2], p->imgstrgrid[3], p->imgstrgrid[4], p->imgstrgrid[5] );

   printf("   node_data_file                       : %s\n\n"      , p->node_data_file  );

   printf("   dataFileVersion                      : %d\n"        , p->dataFileVersion );
   printf("   numFileSegments                      : %d\n"        , p->numFileSegments );
   printf("   nodeCount                            : %d\n"        , p->nodeCount );
   printf("   dataDecompType                       : %d\n"        , p->dataDecompType );
   printf("   dataDecompGeometry                   : [ %d %d %d ]\n"           , p->dataDecompGeometry[0], p->dataDecompGeometry[1], p->dataDecompGeometry[2] );
   printf("   minCoordinates                       : [ %+lf %+lf %+lf ]\n"     , p->minCoordinates    [0], p->minCoordinates    [1], p->minCoordinates    [2] );
   printf("   maxCoordinates                       : [ %+lf %+lf %+lf ]\n\n"   , p->maxCoordinates    [0], p->maxCoordinates    [1], p->maxCoordinates    [2] );

   printf("   simVol                               : %+lf\n"      , p->simVol );
   printf("   burgVolFactor                        : %+lf\n\n"    , p->burgVolFactor );

   printf("   maxNumThreads                        : %d\n\n"      , p->maxNumThreads );

   printf("   ignoreNonLoopFlux                    : %d\n\n"      , p->ignoreNonLoopFlux );

   printf("\n\n-------end of paramters-------\n\n");
};

/*------------------------------------------------------------------------
 *
 *      Functions:    Print_Param_List
 *      Description:  Prints a human-readable summary of the parameter list structure.
 *
 *      Arguments:
 *          p         Pointer to the parameter list structure
 *
 *          name      simple title message added to the output
 *
 *----------------------------------------------------------------------*/

void Print_Param_List (ParamList_t *plist, const char *msg)
{
   if (plist == (ParamList_t *) NULL) return;

   printf("ParaDiS ParamList_t data structure contents %s\n", msg ? msg : "");

   printf("   paramCnt  : %d\n"         ,                 plist->paramCnt);
   printf("   varList   : 0x%012lx\n\n" , (unsigned long) plist->varList );

   if ( (plist->paramCnt>0) && (plist->varList) )
   {
      int        pn = plist->paramCnt;
      VarData_t *pv = plist->varList ;

      const char *type_strs [] = { "null", "float", "int", "string", "comment" };

      for (int i=0; (i<pn); i++, pv++)
      {
         real8  *pdbl = (real8 *) pv->valList;
         int    *pint = (int   *) pv->valList;
         char   *pstr = (char  *) pv->valList;

         printf("      %-60s (%-8s) n=%-3d flg=%2d", pv->varName, type_strs[pv->valType], pv->valCnt, pv->flags);

         if ( (pv->valCnt<1) || !(pv->valList) ) { printf("\n"); }

         else if ( pv->valCnt==1 )
         {
            if (pv->valType==V_DBL    ) { printf(" : %+e\n", pdbl[0] ); }
            if (pv->valType==V_INT    ) { printf(" : %d\n" , pint[0] ); }
            if (pv->valType==V_STRING ) { printf(" : %s\n" , pstr    ); }
            if (pv->valType==V_COMMENT) { printf(" : %s\n" , pstr    ); }
         }

         else if ( pv->valCnt<5 )
         {
            printf(" : [ " );
            for (int j=0; (j<pv->valCnt); j++)
            {
               if (pv->valType==V_DBL) { printf("%+e ", pdbl[j] ); }
               if (pv->valType==V_INT) { printf("%d " , pint[j] ); }
            }
            printf("]\n" );
         }

         else
         {
            for (int j=0; (j<pv->valCnt); j++)
            {
               if (j%10==0)             { printf("\n"); }
               if (pv->valType==V_DBL) { printf("%+e ", pdbl[j] ); }
               if (pv->valType==V_INT) { printf("%d " , pint[j] ); }
            }
            printf("\n\n");
         }
      }
   }
}
