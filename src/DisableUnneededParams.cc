/****************************************************************************
 *
 *      Module:       DisableUnneededParams.c
 *
 *      Description:  This function is called by task 0 during
 *                    initialization to explicitly mark certain
 *                    control parameters in order to prevent them
 *                    from being written out when creating restart
 *                    files.  This is done because not all the known
 *                    control parameters are applicable to every
 *                    simulation.  For instance, the set of mobility
 *                    parameters used during a given simulation
 *                    is dependent on the mobility law in use, and
 *                    will only be a subset of all the available
 *                    mobility parameters.
 *
 ****************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

/*
 *      Need only be called by task zero
 */
void DisableUnneededParams(Home_t *home)
{
        Param_t  *param;

        param = home->param;

/*
 *      Not all mobility parameters are needed by all mobility functions.
 *      The easiest thing to do is first disable all of them and then only
 *      re-enable the ones appropriate to the selected mobility.
 */
        MarkParamDisabled(home->ctrlParamList, "MobScrew");
        MarkParamDisabled(home->ctrlParamList, "MobEdge");
        MarkParamDisabled(home->ctrlParamList, "MobClimb");
        MarkParamDisabled(home->ctrlParamList, "MobGlide");
        MarkParamDisabled(home->ctrlParamList, "MobLine");
        MarkParamDisabled(home->ctrlParamList, "MobNumTrenches");
        MarkParamDisabled(home->ctrlParamList, "MobTrench");
        MarkParamDisabled(home->ctrlParamList, "MobTrenchAngle");
        MarkParamDisabled(home->ctrlParamList, "MobTrenchWidth");
        MarkParamDisabled(home->ctrlParamList, "MobG0");
        MarkParamDisabled(home->ctrlParamList, "MobFileNL");
        MarkParamDisabled(home->ctrlParamList, "TempK");
        MarkParamDisabled(home->ctrlParamList, "MobExpP");
        MarkParamDisabled(home->ctrlParamList, "MobExpQ");
        MarkParamDisabled(home->ctrlParamList, "MobAlpha");
        MarkParamDisabled(home->ctrlParamList, "MobDeltaH0");
        MarkParamDisabled(home->ctrlParamList, "MobCoeffA");
        MarkParamDisabled(home->ctrlParamList, "MobCoeffC");
        MarkParamDisabled(home->ctrlParamList, "MobCoeffC0");
        MarkParamDisabled(home->ctrlParamList, "MobCoeffD");
        MarkParamDisabled(home->ctrlParamList, "MobPeierls");
        MarkParamDisabled(home->ctrlParamList, "meltTemp");
        MarkParamDisabled(home->ctrlParamList, "meltTempFile");
        MarkParamDisabled(home->ctrlParamList, "ecFile");
        MarkParamDisabled(home->ctrlParamList, "pressure");
        MarkParamDisabled(home->ctrlParamList, "cOVERa");
        MarkParamDisabled(home->ctrlParamList, "HCPEcoreA");
        MarkParamDisabled(home->ctrlParamList, "HCPEcoreC");
        MarkParamDisabled(home->ctrlParamList, "HCPEcoreCpA");
        MarkParamDisabled(home->ctrlParamList, "HCP_CA_SplitFrequency");
        MarkParamDisabled(home->ctrlParamList, "HCP_A_Basal_ScrewDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_A_Prismatic_ScrewDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_A_1stPyramidal_ScrewDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_A_2ndPyramidal_ScrewDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_C_Prismatic_ScrewDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_CpA_Prismatic_ScrewDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_CpA_1stPyramidal_ScrewDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_CpA_2ndPyramidal_ScrewDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_Sessile_ScrewDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_A_Basal_EdgeDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_A_Prismatic_EdgeDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_A_1stPyramidal_EdgeDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_A_2ndPyramidal_EdgeDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_C_Prismatic_EdgeDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_CpA_Prismatic_EdgeDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_CpA_1stPyramidal_EdgeDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_CpA_2ndPyramidal_EdgeDrag");

        MarkParamDisabled(home->ctrlParamList, "HCP_Sessile_EdgeDrag");
        MarkParamDisabled(home->ctrlParamList, "HCP_LineDrag");
        MarkParamDisabled(home->ctrlParamList, "MobRelaxX");
        MarkParamDisabled(home->ctrlParamList, "MobRelaxY");
        MarkParamDisabled(home->ctrlParamList, "MobRelaxZ");
        MarkParamDisabled(home->ctrlParamList, "MobRelaxScaleByLength");
        MarkParamDisabled(home->ctrlParamList, "MobFacetedUseCusps");

        MarkParamDisabled(home->ctrlParamList, "MobAlpha"   );
        MarkParamDisabled(home->ctrlParamList, "MobExpP"    );
        MarkParamDisabled(home->ctrlParamList, "MobExpQ"    );
        MarkParamDisabled(home->ctrlParamList, "MobCoeffA"  );
        MarkParamDisabled(home->ctrlParamList, "MobCoeffC"  );
        MarkParamDisabled(home->ctrlParamList, "MobCoeffC0" );
        MarkParamDisabled(home->ctrlParamList, "MobCoeffD"  );

        switch(param->mobilityType) {
            case MOB_BCC_0:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                break;
            case MOB_BCC_0B:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                break;
            case MOB_BCC_GLIDE:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                break;
            case MOB_BCC_LINEAR:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                break;
            case MOB_FCC_0:
            case MOB_FCC_0B:
            case MOB_FCC_CLIMB:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                break;
            case MOB_FCC_LINEAR:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                break;
            case MOB_HCP_LINEAR:
                MarkParamEnabled(home->ctrlParamList, "cOVERa");
                MarkParamEnabled(home->ctrlParamList, "HCP_CA_SplitFrequency");
                MarkParamEnabled(home->ctrlParamList, "HCP_A_Basal_ScrewDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_A_Prismatic_ScrewDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_A_1stPyramidal_ScrewDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_A_2ndPyramidal_ScrewDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_C_Prismatic_ScrewDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_CpA_Prismatic_ScrewDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_CpA_1stPyramidal_ScrewDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_CpA_2ndPyramidal_ScrewDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_Sessile_ScrewDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_A_Basal_EdgeDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_A_Prismatic_EdgeDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_A_1stPyramidal_EdgeDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_A_2ndPyramidal_EdgeDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_C_Prismatic_EdgeDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_CpA_Prismatic_EdgeDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_CpA_1stPyramidal_EdgeDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_CpA_2ndPyramidal_EdgeDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_Sessile_EdgeDrag");
                MarkParamEnabled(home->ctrlParamList, "HCP_LineDrag");
                break;
#ifdef ESHELBY
            case MOB_FCC_0B_ESHELBY:
            case MOB_BCC_0_ESHELBY:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                break;
            case MOB_FCC_LINEAR_ESHELBY:
                MarkParamEnabled(home->ctrlParamList, "MobScrew");
                MarkParamEnabled(home->ctrlParamList, "MobEdge");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                break;
#endif
            case MOB_BCC_FACETED:
                MarkParamEnabled(home->ctrlParamList, "MobGlide");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                MarkParamEnabled(home->ctrlParamList, "MobNumTrenches");
                MarkParamEnabled(home->ctrlParamList, "MobTrench");
                MarkParamEnabled(home->ctrlParamList, "MobTrenchAngle");
                MarkParamEnabled(home->ctrlParamList, "MobFacetedUseCusps");
                if (!param->MobFacetedUseCusps) {
                    MarkParamEnabled(home->ctrlParamList, "MobTrenchWidth");
                }
                break;
            case MOB_BCC_FE_NL:
                MarkParamEnabled(home->ctrlParamList, "TempK");
                break;
            case MOB_BCC_FE_NL_A:
                MarkParamEnabled(home->ctrlParamList, "TempK");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                break;
            case MOB_BCC_TA_NL:
                MarkParamEnabled(home->ctrlParamList, "MobGlide");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                MarkParamEnabled(home->ctrlParamList, "TempK");
                break;
            case MOB_BCC_TA_NL_B:
                MarkParamEnabled(home->ctrlParamList, "MobGlide");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                MarkParamEnabled(home->ctrlParamList, "TempK");
                break;
            case MOB_BCC_TA_NL_B_PLANAR:
                MarkParamEnabled(home->ctrlParamList, "MobGlide");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                MarkParamEnabled(home->ctrlParamList, "TempK");
                break;
            case MOB_BCC_VA_NL:
                MarkParamEnabled(home->ctrlParamList, "MobGlide");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                MarkParamEnabled(home->ctrlParamList, "MobLine");
                MarkParamEnabled(home->ctrlParamList, "MobG0");
                MarkParamEnabled(home->ctrlParamList, "pressure");
                MarkParamEnabled(home->ctrlParamList, "TempK");
                break;
            case MOB_BCC_VA_NL_PLANAR:
                MarkParamEnabled(home->ctrlParamList, "MobGlide");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                MarkParamEnabled(home->ctrlParamList, "MobLine");
                MarkParamEnabled(home->ctrlParamList, "MobG0");
                MarkParamEnabled(home->ctrlParamList, "pressure");
                MarkParamEnabled(home->ctrlParamList, "TempK");
                break;
            case MOB_RHOMBOHEDRAL_VA:
                MarkParamEnabled(home->ctrlParamList, "MobGlide");
                MarkParamEnabled(home->ctrlParamList, "TempK");
                break;
            case MOB_BCC_NL:
                MarkParamEnabled(home->ctrlParamList, "MobExpP");
                MarkParamEnabled(home->ctrlParamList, "MobExpQ");
                MarkParamEnabled(home->ctrlParamList, "MobAlpha");
                MarkParamEnabled(home->ctrlParamList, "MobDeltaH0");
                MarkParamEnabled(home->ctrlParamList, "MobCoeffA");
                MarkParamEnabled(home->ctrlParamList, "MobCoeffC");
                MarkParamEnabled(home->ctrlParamList, "MobCoeffC0");
                MarkParamEnabled(home->ctrlParamList, "MobCoeffD");
                MarkParamEnabled(home->ctrlParamList, "meltTemp");
                MarkParamEnabled(home->ctrlParamList, "meltTempFile");
                MarkParamEnabled(home->ctrlParamList, "MobLine");
                MarkParamEnabled(home->ctrlParamList, "MobGlide");
                MarkParamEnabled(home->ctrlParamList, "MobClimb");
                MarkParamEnabled(home->ctrlParamList, "MobPeierls");
                MarkParamEnabled(home->ctrlParamList, "TempK");
                MarkParamEnabled(home->ctrlParamList, "MobFileNL");
                MarkParamEnabled(home->ctrlParamList, "ecFile");
                break;

            case MOB_RELAX:
            case MOB_RELAX_GLIDE:
                MarkParamEnabled(home->ctrlParamList, "MobRelaxX");
                MarkParamEnabled(home->ctrlParamList, "MobRelaxY");
                MarkParamEnabled(home->ctrlParamList, "MobRelaxZ");
                MarkParamEnabled(home->ctrlParamList, "MobRelaxScaleByLength");
                break;

            case MOB_TA :
                MarkParamEnabled(home->ctrlParamList, "MobAlpha"   );
                MarkParamEnabled(home->ctrlParamList, "MobExpP"    );
                MarkParamEnabled(home->ctrlParamList, "MobExpQ"    );
                MarkParamEnabled(home->ctrlParamList, "MobCoeffA"  );
                MarkParamEnabled(home->ctrlParamList, "MobCoeffC"  );
                MarkParamEnabled(home->ctrlParamList, "MobCoeffC0" );
                MarkParamEnabled(home->ctrlParamList, "MobCoeffD"  );
                break;

            case MOB_TA_LINEAR :
                MarkParamEnabled(home->ctrlParamList, "MobAlpha"   );
                MarkParamEnabled(home->ctrlParamList, "MobExpP"    );
                MarkParamEnabled(home->ctrlParamList, "MobExpQ"    );
                MarkParamEnabled(home->ctrlParamList, "MobCoeffA"  );
                MarkParamEnabled(home->ctrlParamList, "MobCoeffC"  );
                MarkParamEnabled(home->ctrlParamList, "MobCoeffC0" );
                MarkParamEnabled(home->ctrlParamList, "MobCoeffD"  );
                break;


        }  /* end switch(mobilityType) */

/*
 *      If no free surfaces are used, we can drop the free surface
 *      boundary values.
 */
        if ((param->xBoundType == Periodic) &&
            (param->xBoundType == Periodic) &&
            (param->xBoundType == Periodic)) {
            MarkParamDisabled(home->ctrlParamList, "xBoundMin");
            MarkParamDisabled(home->ctrlParamList, "yBoundMin");
            MarkParamDisabled(home->ctrlParamList, "zBoundMin");
            MarkParamDisabled(home->ctrlParamList, "xBoundMax");
            MarkParamDisabled(home->ctrlParamList, "yBoundMax");
            MarkParamDisabled(home->ctrlParamList, "zBoundMax");
        }

/*
 *      If the dynamic load balancing is disabled, turn off the decomposition
 *      type.
 */
        if (param->DLBfreq < 1) {
            MarkParamDisabled(home->ctrlParamList, "decompType");
        }

/*
 *      If osmotic forces are not enabled, disable/enable any params
 *      as needed.
 *
 *      Note: TempK was disabled along with mobility params above but
 *            it is needed to do osmotic forces
 */
        if (param->vacancyConcEquilibrium <= 0)
        {
            MarkParamDisabled(home->ctrlParamList, "vacancyConc");
            MarkParamDisabled(home->ctrlParamList, "vacancyConcEquilibrium");
        }
        else
        {
            MarkParamEnabled(home->ctrlParamList, "TempK");
        }

/*
 *      If FMM is not enabled, disable the corresponding parameters
 *      otherwise diable the parameters related to the Rijm* files
 */
        if (param->fmEnabled == 0) {
            MarkParamDisabled(home->ctrlParamList, "fmMPOrder");
            MarkParamDisabled(home->ctrlParamList, "fmTaylorOrder");
            MarkParamDisabled(home->ctrlParamList, "fmCorrectionTbl");
#ifdef CALCENERGY
            MarkParamDisabled(home->ctrlParamList, "fmECorrTbl");
#endif
        } else {
            MarkParamDisabled(home->ctrlParamList, "Rijmfile");
            MarkParamDisabled(home->ctrlParamList, "RijmPBCfile");
        }

#ifdef ESHELBY
            MarkParamDisabled(home->ctrlParamList, "eshelbyfmEnabled");
            MarkParamDisabled(home->ctrlParamList, "inclusionFile");
            MarkParamDisabled(home->ctrlParamList, "eshelbyfmMPOrder");
            MarkParamDisabled(home->ctrlParamList, "eshelbyfmTaylorOrder");
            MarkParamDisabled(home->ctrlParamList, "MobEshelbyResist");
            MarkParamDisabled(home->ctrlParamList, "shearModulus2");
            MarkParamDisabled(home->ctrlParamList, "pois2");
            MarkParamDisabled(home->ctrlParamList, "Ecore2");
            MarkParamDisabled(home->ctrlParamList, "IncStepEnergy");
#endif // ESHELBY

#ifdef FRCRIT
            MarkParamDisabled(home->ctrlParamList, "minBisection");
            MarkParamDisabled(home->ctrlParamList, "maxBisection");
#endif

/*
 *      Not all of the timestep integrator parameters are used in
 *      all the integrator functions...
 */
        if (StrEquiv(param->timestepIntegrator, "forward-euler"))
        {
            MarkParamDisabled(home->ctrlParamList, "trapezoidMaxIterations");
            MarkParamDisabled(home->ctrlParamList, "dtDecrementFact");
            MarkParamDisabled(home->ctrlParamList, "dtExponent");
            MarkParamDisabled(home->ctrlParamList, "dtIncrementFact");
            MarkParamDisabled(home->ctrlParamList, "dtVariableAdjustment");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_MaxIterations");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_MaxLinIterations");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_NumPriorResiduals");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_PrintLevel");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_LogFile");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_EtaVal");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_IRKtable");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_FPsolver");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_MaxNonLinIters");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_MaxLinIters");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_FPaccel");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_PredictorMethod");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_AdaptivityMethod");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_NonLinCoef");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_LinCoef");
        }
        else if (StrEquiv(param->timestepIntegrator, "trapezoid"))
        {
            MarkParamDisabled(home->ctrlParamList, "rmax");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_MaxIterations");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_MaxLinIterations");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_NumPriorResiduals");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_PrintLevel");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_LogFile");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_EtaVal");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_IRKtable");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_FPsolver");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_MaxNonLinIters");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_MaxLinIters");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_FPaccel");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_PredictorMethod");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_NonLinCoef");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_LinCoef");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_AdaptivityMethod");
        }
        else if (StrEquiv(param->timestepIntegrator, "trapezoid-KINSOL"))
        {
            MarkParamDisabled(home->ctrlParamList, "trapezoidMaxIterations");
            MarkParamDisabled(home->ctrlParamList, "rmax");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_IRKtable");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_FPsolver");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_MaxNonLinIters");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_MaxLinIters");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_FPaccel");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_PredictorMethod");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_AdaptivityMethod");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_NonLinCoef");
            MarkParamDisabled(home->ctrlParamList, "ARKODE_LinCoef");

/*
 *          KINSOL timestep integrator uses Anderson acceleration or
 *          Newton-Krylov solver: disable KINSOl params not appropriate
 *          to the current solver.
 */
            if (param->KINSOL_UseNewton)
            {
                MarkParamDisabled(home->ctrlParamList, "KINSOL_NumPriorResiduals");
            }
            else
            {
                MarkParamDisabled(home->ctrlParamList, "KINSOL_MaxLinIterations");
                MarkParamDisabled(home->ctrlParamList, "KINSOL_EtaVal");
            }

        }
        else if (StrEquiv(param->timestepIntegrator, "ARKODE"))
        {
            MarkParamDisabled(home->ctrlParamList, "rmax");
            MarkParamDisabled(home->ctrlParamList, "trapezoidMaxIterations");
            MarkParamDisabled(home->ctrlParamList, "dtExponent");
            MarkParamDisabled(home->ctrlParamList, "dtIncrementFact");
            MarkParamDisabled(home->ctrlParamList, "dtVariableAdjustment");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_MaxIterations");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_NumPriorResiduals");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_PrintLevel");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_LogFile");
            MarkParamDisabled(home->ctrlParamList, "KINSOL_EtaVal");
        }

/*
 *      The selected <loadType> affects which parameters are used.
 *      do that setup now.  As with the mobility parameters above,
 *      for some of the parameters, it's easier to disable a group
 *      of them and only re-enable them as needed.
 */
        MarkParamDisabled(home->ctrlParamList, "cTimeOld");
        MarkParamDisabled(home->ctrlParamList, "dCyclicStrain");
        MarkParamDisabled(home->ctrlParamList, "netCyclicStrain");
        MarkParamDisabled(home->ctrlParamList, "numLoadCycle");
        MarkParamDisabled(home->ctrlParamList, "eAmp");

        switch(param->loadType)
        {
            case 0: break;
            case 1: break;
            case 2: break;
            case 3: break;
/*
 *          For loadType's 4 and 5, we use re-enabled the same set
 *          of parameters
 */
            case 4:
            case 5: MarkParamEnabled(home->ctrlParamList, "cTimeOld");
                    MarkParamEnabled(home->ctrlParamList, "dCyclicStrain");
                    MarkParamEnabled(home->ctrlParamList, "netCyclicStrain");
                    MarkParamEnabled(home->ctrlParamList, "numLoadCycle");
                    MarkParamEnabled(home->ctrlParamList, "eAmp");
                    break;
        }

/*
 *      If inertial terms are not being included, don't write parameters
 *      specific to that capability.
 */
        if (param->includeInertia == 0) {
            MarkParamDisabled(home->ctrlParamList, "massDensity");
        }

/*
 *      If no sessile burgers vectors have been specified, don't write
 *      the related arrays out.
 */
        if (param->sessileburgspec[0] <= 0) {
            MarkParamDisabled(home->ctrlParamList, "sessileburgspec");
            MarkParamDisabled(home->ctrlParamList, "sessilelinespec");
        }

/*
 *      If the geometries are not specified in a user-defined laboratory
 *      frame but use the standard crystalographic frame, don't write
 *      out the rotation matrix for the lab frame
 */
        if (param->useLabFrame == 0) {
            MarkParamDisabled(home->ctrlParamList, "rotationMatrix");
        }

/*
 *      The flux data arrays always get dumped to the control file, but
 *      which arrays are written depends on the type of material, so
 *      only enable the ones we need.
 */
        MarkParamDisabled(home->ctrlParamList, "BCC_LtotDecomp");
        MarkParamDisabled(home->ctrlParamList, "BCC2_LtotDecomp");
        MarkParamDisabled(home->ctrlParamList, "BCC_fluxtot");
        MarkParamDisabled(home->ctrlParamList, "BCC2_fluxtot");
        MarkParamDisabled(home->ctrlParamList, "FCC_LtotDecomp");
        MarkParamDisabled(home->ctrlParamList, "FCC_fluxtot");
        MarkParamDisabled(home->ctrlParamList, "HCP_LtotDecomp");
        MarkParamDisabled(home->ctrlParamList, "HCP_fluxtot");
        MarkParamDisabled(home->ctrlParamList, "rhombo_Ltot");
        MarkParamDisabled(home->ctrlParamList, "rhombo_fluxtot");

        // And a few other parameters tied to material type...
        MarkParamDisabled(home->ctrlParamList, "splinterSegFrequency");
        MarkParamDisabled(home->ctrlParamList, "HCP_CA_SplitFrequency");

        switch (param->materialType) {
            case MAT_TYPE_BCC:
                if (param->bcc_DensityFluxDecomp==2)
                {
                   MarkParamEnabled(home->ctrlParamList, "BCC2_LtotDecomp");
                   MarkParamEnabled(home->ctrlParamList, "BCC2_fluxtot");
                }
                else
                {
                   MarkParamEnabled(home->ctrlParamList, "BCC_LtotDecomp");
                   MarkParamEnabled(home->ctrlParamList, "BCC_fluxtot");
                }
                break;
            case MAT_TYPE_FCC:
                MarkParamEnabled(home->ctrlParamList, "FCC_LtotDecomp");
                MarkParamEnabled(home->ctrlParamList, "FCC_fluxtot");
                break;
            case MAT_TYPE_HCP:
                MarkParamEnabled(home->ctrlParamList, "HCP_LtotDecomp");
                MarkParamEnabled(home->ctrlParamList, "HCP_fluxtot");
                MarkParamEnabled(home->ctrlParamList, "HCPEcoreA");
                MarkParamEnabled(home->ctrlParamList, "HCPEcoreC");
                MarkParamEnabled(home->ctrlParamList, "HCPEcoreCpA");
                MarkParamEnabled(home->ctrlParamList, "splinterSegFrequency");
                MarkParamEnabled(home->ctrlParamList, "HCP_CA_SplitFrequency");
                break;
            case MAT_TYPE_RHOMBOHEDRAL_VA:
                MarkParamEnabled(home->ctrlParamList, "rhombo_Ltot");
                MarkParamEnabled(home->ctrlParamList, "rhombo_fluxtot");
                break;
        }

/*
 *      Disable writing of the output-related parameters that do not
 *      apply to the types of output actually selected.
 */
        if (param->armfile == 0) {
            MarkParamDisabled(home->ctrlParamList, "armfilefreq");
            MarkParamDisabled(home->ctrlParamList, "armfiledt");
            MarkParamDisabled(home->ctrlParamList, "armfiletime");
            MarkParamDisabled(home->ctrlParamList, "armfilecounter");
        };

        if (param->armfiledt > 0.0) {
            MarkParamDisabled(home->ctrlParamList, "armfilefreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "armfiledt");
            MarkParamDisabled(home->ctrlParamList, "armfiletime");
        }


        if (param->writeFlux == 0) {
            MarkParamDisabled(home->ctrlParamList, "writeFluxEdgeDecomp");
            MarkParamDisabled(home->ctrlParamList, "writeFluxFullDecomp");
            MarkParamDisabled(home->ctrlParamList, "writeFluxFullDecompTotals");
            MarkParamDisabled(home->ctrlParamList, "writeFluxSimpleTotals");
            MarkParamDisabled(home->ctrlParamList, "writeFluxFreq");
            MarkParamDisabled(home->ctrlParamList, "writeFluxDT");
            MarkParamDisabled(home->ctrlParamList, "writeFluxTime");
            MarkParamDisabled(home->ctrlParamList, "writeFluxCounter");
        }


        if (param->writeFluxDT > 0.0) {
            MarkParamDisabled(home->ctrlParamList, "writeFluxFreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "writeFluxDT");
            MarkParamDisabled(home->ctrlParamList, "writeFluxTime");
        }


        if (param->gnuplot == 0) {
            MarkParamDisabled(home->ctrlParamList, "gnuplotfreq");
            MarkParamDisabled(home->ctrlParamList, "gnuplotdt");
            MarkParamDisabled(home->ctrlParamList, "gnuplottime");
            MarkParamDisabled(home->ctrlParamList, "gnuplotcounter");
        }

        if (param->gnuplotdt > 0.0) {
            MarkParamDisabled(home->ctrlParamList, "gnuplotfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "gnuplotdt");
            MarkParamDisabled(home->ctrlParamList, "gnuplottime");
        }


        if (param->polefigfile == 0) {
            MarkParamDisabled(home->ctrlParamList, "polefigfreq");
            MarkParamDisabled(home->ctrlParamList, "polefigdt");
            MarkParamDisabled(home->ctrlParamList, "polefigtime");
            MarkParamDisabled(home->ctrlParamList, "polefilecounter");
        }

        if (param->polefigdt > 0) {
            MarkParamDisabled(home->ctrlParamList, "polefigfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "polefigdt");
            MarkParamDisabled(home->ctrlParamList, "polefigtime");
        }


        if (param->povray == 0) {
            MarkParamDisabled(home->ctrlParamList, "povrayfreq");
            MarkParamDisabled(home->ctrlParamList, "povraydt");
            MarkParamDisabled(home->ctrlParamList, "povraytime");
            MarkParamDisabled(home->ctrlParamList, "povraycounter");
        }

        if (param->povraydt > 0) {
            MarkParamDisabled(home->ctrlParamList, "povrayfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "povraydt");
            MarkParamDisabled(home->ctrlParamList, "povraytime");
        }


        if (param->psfile == 0) {
            MarkParamDisabled(home->ctrlParamList, "psfilefreq");
            MarkParamDisabled(home->ctrlParamList, "psfiledt");
            MarkParamDisabled(home->ctrlParamList, "psfiletime");
        }

        if (param->psfile > 0) {
            MarkParamDisabled(home->ctrlParamList, "psfilefreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "psfiledt");
            MarkParamDisabled(home->ctrlParamList, "psfiletime");
        }


        if (param->savecn == 0) {
            MarkParamDisabled(home->ctrlParamList, "savecnfreq");
            MarkParamDisabled(home->ctrlParamList, "savecndt");
            MarkParamDisabled(home->ctrlParamList, "savecntime");
            MarkParamDisabled(home->ctrlParamList, "savecncounter");
        }

        if (param->savecndt > 0) {
            MarkParamDisabled(home->ctrlParamList, "savecnfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "savecndt");
            MarkParamDisabled(home->ctrlParamList, "savecntime");
        }


        if (param->saveprop == 0) {
            MarkParamDisabled(home->ctrlParamList, "savepropfreq");
            MarkParamDisabled(home->ctrlParamList, "savepropdt");
            MarkParamDisabled(home->ctrlParamList, "saveproptime");
        }

        if (param->savepropdt > 0) {
            MarkParamDisabled(home->ctrlParamList, "savepropfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "savepropdt");
            MarkParamDisabled(home->ctrlParamList, "saveproptime");
        }


        if (param->savetimers == 0) {
            MarkParamDisabled(home->ctrlParamList, "savetimersfreq");
            MarkParamDisabled(home->ctrlParamList, "savetimersdt");
            MarkParamDisabled(home->ctrlParamList, "savetimerstime");
            MarkParamDisabled(home->ctrlParamList, "savetimerscounter");
        }

        if (param->savetimersdt > 0) {
            MarkParamDisabled(home->ctrlParamList, "savetimersfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "savetimersdt");
            MarkParamDisabled(home->ctrlParamList, "savetimerstime");
        }


        if (param->tecplot == 0) {
            MarkParamDisabled(home->ctrlParamList, "tecplotfreq");
            MarkParamDisabled(home->ctrlParamList, "tecplotdt");
            MarkParamDisabled(home->ctrlParamList, "tecplottime");
            MarkParamDisabled(home->ctrlParamList, "tecplotcounter");
        }

        if (param->tecplotdt > 0) {
            MarkParamDisabled(home->ctrlParamList, "tecplotfreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "tecplotdt");
            MarkParamDisabled(home->ctrlParamList, "tecplottime");
        }


        if (param->velfile == 0) {
            MarkParamDisabled(home->ctrlParamList, "velfilefreq");
            MarkParamDisabled(home->ctrlParamList, "velfiledt");
            MarkParamDisabled(home->ctrlParamList, "velfiletime");
            MarkParamDisabled(home->ctrlParamList, "velfilecounter");
        }

        if (param->velfiledt > 0) {
            MarkParamDisabled(home->ctrlParamList, "velfilefreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "velfiledt");
            MarkParamDisabled(home->ctrlParamList, "velfiletime");
        }


        if (param->writeForce == 0) {
            MarkParamDisabled(home->ctrlParamList, "writeForceFreq");
            MarkParamDisabled(home->ctrlParamList, "writeForceDT");
            MarkParamDisabled(home->ctrlParamList, "writeForceTime");
            MarkParamDisabled(home->ctrlParamList, "writeForceCounter");
        }

        if (param->writeForceDT > 0) {
            MarkParamDisabled(home->ctrlParamList, "writeForceFreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "writeForceDT");
            MarkParamDisabled(home->ctrlParamList, "writeForceTime");
        }


        if ((param->savedensityspec[0] == 0) ||
            (param->savedensityspec[1] == 0) ||
            (param->savedensityspec[2] == 0)) {
            MarkParamDisabled(home->ctrlParamList, "savedensityspec");
        }

        if (param->writeVisit == 0) {
            MarkParamDisabled(home->ctrlParamList, "writeVisitFreq");
            MarkParamDisabled(home->ctrlParamList, "writeVisitDT");
            MarkParamDisabled(home->ctrlParamList, "writeVisitTime");
            MarkParamDisabled(home->ctrlParamList, "writeVisitCounter");
            MarkParamDisabled(home->ctrlParamList, "writeVisitNodes");
            MarkParamDisabled(home->ctrlParamList, "writeVisitNodesAsText");
            MarkParamDisabled(home->ctrlParamList, "writeVisitSegments");
            MarkParamDisabled(home->ctrlParamList, "writeVisitSegmentsAsText");
            MarkParamDisabled(home->ctrlParamList, "writeVisitBurgID");
            MarkParamDisabled(home->ctrlParamList, "writeVisitForceVector");
            MarkParamDisabled(home->ctrlParamList, "writeVisitVelocityVector");
        }

        if (param->writeVisitDT > 0) {
            MarkParamDisabled(home->ctrlParamList, "writeVisitFreq");
        } else {
            MarkParamDisabled(home->ctrlParamList, "writeVisitDT");
            MarkParamDisabled(home->ctrlParamList, "writeVisitTime");
        }

        if (param->writeVisitNodes == 0) {
            MarkParamDisabled(home->ctrlParamList, "writeVisitNodesAsText");
            MarkParamDisabled(home->ctrlParamList, "writeVisitForceVector");
            MarkParamDisabled(home->ctrlParamList, "writeVisitVelocityVector");
        }

        if (param->writeVisitSegments == 0) {
            MarkParamDisabled(home->ctrlParamList, "writeVisitSegmentsAsText");
            MarkParamDisabled(home->ctrlParamList, "writeVisitBurgID");
        }

        return;
}
