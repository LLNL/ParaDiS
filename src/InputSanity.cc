
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Param.h"
#include "Util.h"
#include "V3.h"

/***************************************************************************
 *
 *  Function    : InputSanity
 *  Description : Do some checks on the consistency of the input
 *
 *      Last Modified: 04/08/2008 gh - Added check that FMM is used if
 *                                     PBC is disabled.
 *
 **************************************************************************/

void InputSanity (Home_t *home)
{
    Param_t *param = home->param;

	int domainCount = param->nXdoms * param->nYdoms * param->nZdoms;

#ifdef PARALLEL
/*
 *	If we are building for parallel execution, verify that the specified 
 *  geometry matches the domain count...
 */
	if (domainCount != home->numDomains) 
    {
		Fatal("%s; geometry (%dx%dx%d) mismatch with domain count %d\n", 
              __func__, param->nXdoms, param->nYdoms, param->nZdoms, home->numDomains);
	}
#else
/*
 *  If we are building serially, ignore/override the control file domain settings...
 */
	if (domainCount != home->numDomains) 
    {
		Warning("%s() ; ignoring domain geometry (%dx%dx%d) for serial execution",
		      __func__, param->nXdoms, param->nYdoms, param->nZdoms );

        param->nXdoms = 1;
        param->nYdoms = 1;
        param->nZdoms = 1;
	}
#endif

/*
 *      If the user elects to do parallel IO, insure that the IO
 *      parallelism does not exceed the number of domains.
 */
        if ((param->numIOGroups < 0) || (param->numIOGroups > domainCount)) {
            printf( "WARNING: The <numIOGroups> control parameter value (%d) "
                    "exceeds the domain count.\n Resetting it to %d\n",
                    param->numIOGroups, domainCount);
            param->numIOGroups = domainCount;
        }

/*
 *	turn off dynamic load balance if uniprocessor run
 */
	if (param->nXdoms * param->nYdoms * param->nZdoms == 1) { param->DLBfreq = 0; }

/*
 *	By default we want to dump timing results into a file
 *	every time a restart file is created.  If the user did
 *	not explicitly set variables controlling dumps of timing
 *	data, go ahead and use the same control values as for
 *	dumping restarts.
 */
	if (param->savetimers < 0)
		param->savetimers = param->savecn;

	if (param->savetimers > 0) {
		if (param->savetimersfreq < 1) {
			if (param->savecnfreq > 0) 
				param->savetimersfreq = param->savecnfreq;
			else
				param->savetimersfreq = 100;
		}
		if (param->savetimerscounter < 1) {
			if (param->savecnfreq > 0) 
				param->savetimerscounter = param->savecncounter;
			else
				param->savetimerscounter = 1;
		}
		if (param->savetimersdt < 0.0) {
			if (param->savecndt > 0.0) {
				param->savetimersdt = param->savecndt;
				param->savetimerstime = param->savecntime;
			} else {
				param->savetimersdt = 0.0;
				param->savetimerstime = 0.0;
			}
		}
	}

/*
 *      If there are no Eshelby inclusions being simulated, explicitly turn
 *      off the eshelby FMM (whether it was on or not).
 */
#ifdef ESHELBY
        if (param->inclusionFile[0] == 0) {
            param->eshelbyfmEnabled = 0;
        }
#endif

/*
 *    The current implementation of the fast-multipole code to deal
 *    with remote seg forces requires that the number of cells
 *    be identical in each dimension as well as being a power of 2.
 *    Here we just enforce that.
 */
#ifdef ESHELBY
        if ((param->fmEnabled) || (param->eshelbyfmEnabled)) 
#else
        if  (param->fmEnabled)
#endif
        {
            if ((param->nXcells != param->nYcells) ||
                (param->nXcells != param->nZcells)) 
            {
                Fatal("Cell geometry (%dX%dX%d) invalid.  %s",
                      param->nXcells, param->nYcells, param->nZcells,
                      "Number of cells must be identical in each dimension");
            }

            for (int i=2; (i<10); i++) 
            {
                int cellCount = 1 << i;

                if (param->nXcells == cellCount) 
                {
                    param->fmNumLayers = i+1;
                    break;
                }
    
                if (param->nXcells < cellCount)
                    Fatal("Number of cells in a dimension must be at least 4 and a power of 2");
            }
        }

/*
 *	If the fast multipole code is enabled for remote force
 *	calculations, the number of points along each segment
 *      at which to evaluate stress from the expansions
 *      is a function of the expansion order... but don't
 *      go higher than 6.
 */
#ifdef ESHELBY
	if ((param->fmEnabled) || (param->eshelbyfmEnabled)) 
#else
	if  (param->fmEnabled)
#endif
    {
            param->fmNumPoints = (param->fmExpansionOrder / 2) + 1;
            param->fmNumPoints = MIN(param->fmNumPoints, 6);
	}

/*
 *      Verify the domain decomposition type specified is valid.
 */
        if ((param->decompType < 1) || (param->decompType > 2)) {
            Fatal("decompType=%d is invalid.  Type must be 1 or 2\n", param->decompType);
        }

/*
 *      Do some validation of the boundary conditions.  We don't
 *      really support anything other than periodic boundaries
 *      or free surfaces yet, and the code will not properly handle
 *      a mixture of periodic boundaries and free surfaces in the
 *      same simulation.
 */
        if (((param->xBoundType != Periodic)&&(param->xBoundType != Free)) ||
            ((param->yBoundType != Periodic)&&(param->yBoundType != Free)) ||
            ((param->zBoundType != Periodic)&&(param->zBoundType != Free))) {
            Fatal("The only boundary types currently supported are types\n"
                  "       0 (periodic) and 1 (free surfaces).  At least\n"
                  "       one of the specified *BoundType control \n"
                  "       parameters is invalid\n");
        }

#if 0 // (disabled for now - BMW)
#ifndef _ARLFEM  
        if (((param->xBoundType == Free) ||
             (param->yBoundType == Free) ||
             (param->zBoundType == Free)) &&
            (param->fmEnabled == 0)) {
            Fatal("If free surfaces are used, the far-field forces must\n"
                  "       be calculated using the Fast Multipole code.\n"
                  "       (i.e. 'fmEnabled' parameter must be set to 1).\n");
        }
#endif  // _ARLFEM
#endif  // if 0

/*
 *      If the user wants to write binary restart files, make sure
 *      HDF5 support has been enabled.
 */
        if (param->writeBinRestart) {
#ifndef USE_HDF
            Fatal("Program must be compiled with HDF support (see HDF_MODE\n"
                  "in makefile.setup) to use binary restart file capability!");
#endif
        }

/*
 *      If the user wants the mobility functions to include inertial
 *      terms (if available), the associated mass density MUST be
 *      supplied as well.
 */
        if ((param->includeInertia) && (param->massDensity < 0.0)) {
            Fatal("The <includeInertia> toggle is enabled, but the required\n"
                  "<massDensity> control parameter has not been specified.");
        }
            
/*
 *      If the <loadType> is zero, explicitly set <eRate> to 1 and
 *      if <edotdir> is zeroed, set it to the default so there are
 *      no problems in plastic strain calculations.
 */
        if ( (param->loadType==0) && (param->eRate!=1.0) )
        {
           printf("Warning: eRate=%le is inconsistent with loadType=0, setting eRate to 1.0\n",  param->eRate);
           param->eRate = 1.0;
        }

/*
 *      Make sure the frequency at which to do multi-node splits is > 0.
 */
        param->splitMultiNodeFreq = MAX(1, param->splitMultiNodeFreq);

/*
 *      If the user is defining 'trenches' for a faceted mobility
 *      function, do a couple quick sanity checks.
 */
        if (param->MobNumTrenches < 0) {
            Fatal("The <MobNumTrenches> for faceted mobilities is set to\n"
                  "%d: a negative number of trenches is not permitted.",
                  param->MobNumTrenches);
        }

        for (int i = 0; i < param->MobNumTrenches; i++) 
        {
            if ((param->MobTrenchAngle[i] < -90.0) ||
                (param->MobTrenchAngle[i] > 90.0)) {
                Fatal("The <MobTrenchAngle[%d]> parameter is set to %lf:"
                      "range is -90 <= angle <= 90.", i,
                      param->MobTrenchAngle[i]);
            }
            if (param->MobTrenchWidth[i] < 0.0) {
                Fatal("The <MobTrenchWidth[%d]> parameter is set to %lf:"
                      "a negtive width is not permitted.", i,
                      param->MobTrenchWidth[i]);
            }
        }

/*
 *      For HCP materials, the cOVERa coefficient must be in the range
 *      1.0 < cOVERa < 2.0 or the burgers vectors and glide planes change.
 *      We'll actually enforce a slightly narrower range here.
 */
        if (StrEquiv(param->mobilityLaw, "HCP_linear")) {
            if ((param->cOVERa <= 1.05) || (param->cOVERa >= 1.95)) {
                Fatal("cOVERa value (%e) is invalid.  Must be within the \n"
                      "range 1.05 < cOVERa < 1.95", param->cOVERa);
            }
        }

/*
 *      If VisIt file output has been enabled, disable it if neither of the
 *      specific VisIt output types has been selected.
 */
        if (param->writeVisit) {
            if ((param->writeVisitSegments == 0) &&
                (param->writeVisitNodes == 0)) {
                printf( "The <writeVisit> control parameter has been been "
                        "disabled\nbecause neither <writeVisitSegments> "
                        "nor <writeVisitNodes>\nhas been enabled.\n");
            }
        }

/*
 *      If the 'trapezoid' timestep integrator is being used, do a quick
 *      sanity check on the max iteration count.
 */
        if (StrEquiv(param->timestepIntegrator, "trapezoid")) {
            if (param->trapezoidMaxIterations < 1) {
                printf("WARNING: The <trapezoidMaxIterations> parameter (%d) "
                       "must be greater than zero.\nSetting value to default "
                       "of 2.\n", param->trapezoidMaxIterations); 
                param->trapezoidMaxIterations = 2;
            }
        }

/*
 *      If KINSOL has been enabled, do some quick checks on the related
 *      input parameters
 */
#ifdef USE_KINSOL
        if ((param->KINSOL_PrintLevel < 0) || (param->KINSOL_PrintLevel > 3)) {
            printf("The <KINSOL_PrintLevel> (%d) is outside the range 0-3.\n"
                   "Parameter is being disabled\n", param->KINSOL_PrintLevel); 
            param->KINSOL_PrintLevel = 0;
        }

        if (param->KINSOL_MaxIterations < 1) {
            printf("The <KINSOL_MaxIterations> (%d) is invalid and being\n"
                   "reset to 1\n", param->KINSOL_MaxIterations); 
            param->KINSOL_MaxIterations = 1;
        }

        if ((param->KINSOL_NumPriorResiduals < 0) ||
            (param->KINSOL_NumPriorResiduals >= param->KINSOL_MaxIterations)) {
            printf("The <KINSOL_NumPriorResiduals> (%d) must be in the range\n"
                   "0 <= KINSOL_NumPriorResiduals < KINSOL_MaxIterations.\n"
                   "Resetting it to %d\n", param->KINSOL_NumPriorResiduals,
                   param->KINSOL_MaxIterations - 1);
        }
#endif

/*
 *      Do any checks related to anisotropy
 */
#ifdef ANISOTROPIC
        if ((param->anisoHarmonicsNumTermsBase > 20) ||
            (param->anisoHarmonicsNumTermsBase <  0)) {
                Fatal("The <anisoHarmonicsNumTermsBase> parameter is set to "
                      "%d.\nPermitted range is >= 0 and <= 20.",
                      param->anisoHarmonicsNumTermsBase);
        }

#ifdef TAYLORFMM
        if (param->fmEnableAnisotropy) {
            printf("Warning: Code not fully debugged for Anisotropy and Taylor FMM on."
                   " Turning off fmEnableAnisotropy flag!\n");
            param->fmEnableAnisotropy = 0;
        }
#endif

#else
        if (param->fmEnableAnisotropy) {
            printf("Warning: Code has not been compiled with support for\n"
                   "anisotropy.  Turning off fmEnableAnisotropy flag!\n");
            param->fmEnableAnisotropy = 0;
        }
#endif

/*
 *      The HCP_CA_SplitFrequency must be >= zero.  Zero disables the
 *      c+a splitting, a positive integer indicates the cycle frequency
 *      at which splitting is attempted.
 */
        if (param->HCP_CA_SplitFrequency < 0) {
            printf("Warning: HCP_CA_SplitFrequency was set to a negative "
                   "value.\nResetting it to zero!\n");
            param->HCP_CA_SplitFrequency = 0;
        }

#if defined ANISOTROPIC && defined _ARLFEM
/*
 *      When both anisotropy and FEM are in use, issue a warning that
 *      the contribution to the force from the FEM code is isotropic
 *      in nature.
 */
        printf("Warning: Both ANISOTROPIC and _ARLFEM have been defined\n"
               "         but the FEM component of the nodal force is\n"
               "         currently isotropic in nature\n");
#endif 

	return;
}
