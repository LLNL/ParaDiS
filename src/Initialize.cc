//---------------------------------------------------------------------------
//      Author:       Gregg Hommes
//
//      Module:       Initialize.c
//
//      Description:  Contains the driver routine for the initialization
//                    of the application. Handles some of the general
//                    initializations directly and calls all the more
//                    specific initialization routines.
//
//      Last Modified: 01/09/08: Gregg Hommes - Added VerifyBurgersVectors()
//                               sanity check.
//                     04/08/2008 gh - Added support for Vanadium non-linear
//                                     mobility function.  Modified
//                                     SetMobilityParameters() to do
//                                     initialization appropriate to the
//                                     selected mobility function.
//---------------------------------------------------------------------------

// Define INIT_MOB_ATTR_LIST so when we include Mobility.h we
// initialize the mobility module attribute list. All other
// modules will declare the attribute list as an external.

#define INIT_MOB_ATTR_LIST

#include <stdio.h>
#include <memory.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/utsname.h>
#include <fcntl.h>
#include <ctype.h>
#include <time.h>
#include <pwd.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Home.h"
#include "Init.h"
#include "InData.h"
#include "DisplayC.h"
#include "FM.h"
#include "Mobility.h"
#include "Decomp.h"
#include "Parse.h"
#include "Restart.h"
#include "Comm.h"
#include "Timestep.h"
#include "SSF_Driver.h"

#ifndef ANISOTROPIC
//---------------------------------------------------------------------------
//
//      Function:     InitIsotropicElasticConstantMatrix
//      Description:  This function creates the elastic constant matrix
//                    for isotropic simulations based on the provided
//                    shear modulus and poisson ratio.
//
//                    NOTE: This function only populates specific
//                          elements of the matrix, therefore the caller
//                          MUST insure that the full matrix was zeroed
//                          out before this function is called!
//
//---------------------------------------------------------------------------

static void InitIsotropicElasticConstantMatrix(Home_t *home)
{

        real8 shr  = home->param->shearModulus;
        real8 pois = home->param->pois;

        real8 (*ecMatrix)[6] = home->param->elasticConstantMatrix;

        // NOTE: All components not explicitly initialized here are
        //       assumed to have been zeroed out prior to this call.

        ecMatrix[0][0] = 2.0 * shr * (1.0 - pois) / (1.0 - 2.0 * pois);
        ecMatrix[0][1] = 2.0 * shr *        pois  / (1.0 - 2.0 * pois);
        ecMatrix[0][2] = 2.0 * shr *        pois  / (1.0 - 2.0 * pois);

        ecMatrix[1][0] = 2.0 * shr *        pois  / (1.0 - 2.0 * pois);
        ecMatrix[1][1] = 2.0 * shr * (1.0 - pois) / (1.0 - 2.0 * pois);
        ecMatrix[1][2] = 2.0 * shr *        pois  / (1.0 - 2.0 * pois);

        ecMatrix[2][0] = 2.0 * shr *        pois  / (1.0 - 2.0 * pois);
        ecMatrix[2][1] = 2.0 * shr *        pois  / (1.0 - 2.0 * pois);
        ecMatrix[2][2] = 2.0 * shr * (1.0 - pois) / (1.0 - 2.0 * pois);

        ecMatrix[3][3] = shr;
        ecMatrix[4][4] = shr;
        ecMatrix[5][5] = shr;

        return;
}
#endif  // !ANISOTROPIC

static void Usage(const char *program)
{
        const char *usage =
        "\n"
        "Usage: %s [args] <ctrl_file>\n"
        "\n"
        "   where [args] = optional command line args...\n"
        "\n"
        "      [-d    <data_file>    ]  : explicit path to data file                              \n"
        "                                   note - if the path to the data file is not provided,  \n"
        "                                   ParaDiS will attempt to guess the path using the path \n"
        "                                   to the control file and other settings.               \n"
        "                                      e.g.  <data_file> == <ctrl_file>+'.data'           \n"
        "                                                                                         \n"
        "      [-b                   ]  : sets binary data reads                                  \n"
        "      [-s                   ]  : skips IO pass                                           \n"
        "                                                                                         \n"
        "      [-r    <dlb_cycles>   ]  : sets initial dynamic load balance cycles prior to start \n"
        "                                                                                         \n"
        "      [-n       <n>         ]  : sets OpenMP max thread count (only if OpenMP enabled)   \n"
        "                                                                                         \n"
        "      [-f    loops          ]  : sets flux selection to loops                            \n"
        "                                                                                         \n"
        "      [-maxstep <n>         ]  : sets the maximum number of cycles to run                \n"
        "                                                                                         \n"
        "      [-doms  <nx ny nz>    ]  : sets the domain geometry to (nx x ny x nz)              \n"
        "      [-cells <nx ny nz>    ]  : sets the cell   geometry to (nx x ny x nz)              \n"
        "      [-gpu                 ]  : enables GPU acceleration (if available)                 \n"
        "\n";

        printf(usage,program);
        Fatal("Incorrect command line argument list");
}

//---------------------------------------------------------------------------
//
//      Function:     InitRecycleNodeHeap
//      Description:  If the user requested preservation (if possible)
//                    of the node tags found in the restart file, the
//                    <home->nodeKeys> array may be sparsely populated right
//                    from the start.  In this case, we have to
//                    create an initial heap of recycled nodes containing
//                    the indices of all unused <home->nodeKeys> entries
//                    less than <home->newNodeKeyPtr>
//
//---------------------------------------------------------------------------

static void InitRecycleNodeHeap(Home_t *home)
{
        for (int i=0; (i<home->newNodeKeyPtr); i++) {
            if (home->nodeKeys[i] == (Node_t *) NULL) {
                RecycleNodeTag(home, i);
            }
        }

}

//---------------------------------------------------------------------------
//
//      Function:     OpenDir
//      Description:  This function will create (if they do not yet exist)
//                    the primary output directory for the run, plus all
//                    necessary subdirectories for specific output types.
//
//---------------------------------------------------------------------------

int OpenDir(Home_t *home)
{
        char *dirname = home->param->dirname;

        // Only domain zero creates the primary output directory; don't
        // want thousands of processors beating on the file system.

        if (home->myDomain==0)
        {
            if (mkdir(dirname, S_IRWXU) != 0)
            {
                if (errno == EEXIST) {
                    printf("Warning: %s already exists\n", dirname);
                } else {
                    Fatal("Open error %d on directory %s", errno, dirname);
                }
            }
        }

        // All processes wait for task zero to create the primary output
        // directory then cd into that directory.

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        if (chdir(dirname) != 0) {
            Fatal("Task %d: Unable to cd into directory %s", home->myDomain, dirname);
        }

#ifdef PARALLEL
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (home->myDomain == 0) printf("chdir successful on all tasks.\n");

        // Create all subdirectories needed for specific types of output.
        // Again, only domain zero need do these creates.

        // Note: The current working directory for all tasks is now the
        // user specified output directory, so when we create the
        // subdirectories for various output types we just create them
        // local to the current working directory.

        if (home->myDomain == 0)
        {
            char subdir[256];

            if (home->param->armfile     ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_ARMDATA   ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->writeFlux   ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_FLUXDATA  ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->writeForce  ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_FORCE     ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->gnuplot     ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_GNUPLOT   ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->polefigfile ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_POLEFIG   ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->povray      ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_POVRAY    ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->saveprop    ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_PROPERTIES); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->savedt      ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_DEBUG     ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->save_narms  ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_DEBUG     ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->savecn      ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_RESTART   ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->tecplot     ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_TECPLOT   ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->savetimers  ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_TIMERS    ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->velfile     ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_VELOCITY  ); (void) mkdir(subdir,S_IRWXU); }
            if (home->param->writeVisit  ) { snprintf(subdir, sizeof(subdir), "./%s", DIR_VISIT     ); (void) mkdir(subdir,S_IRWXU); }
        }

        return(0);
}

//---------------------------------------------------------------------------
//
//      Function:     GetMeltTemp
//      Description:  Read a table of material specific melting temperatures
//                    vs. pressures from the provided file and calculate the
//                    appropriate melting temperature for the specified
//                    pressure.
//
//      Arguments:
//          param    pointer to the parameter structure containing the
//                   user provided pressure and in which to return to
//                   the caller the melting temperature.
//
//---------------------------------------------------------------------------

void GetMeltTemp(Param_t *param)
{
        real8 **columnData=0, *pressure=0, *meltTemp=0, *meltTempD2=0;

        // Read the table of melting temperatures vs. pressures.
        // Data consists of two columns of data; first is melting
        // temperatures in units if degrees K and second is
        // pressures in units of Mbars.

        int numRows = 0;
        int numCols = 2;

        ReadTabulatedData(param->meltTempFile, numCols, &columnData, &numRows);

        meltTemp = columnData[0];
        pressure = columnData[1];

        // Use a cubic spline method for interpolating the melting
        // temperature for the given pressure based on the table data.

        // NOTE: need a unit conversion since input pressure is
        // in Pa, but table pressures are in MBar.

        real8 targetPressure = param->pressure / 1.0e+11;

        meltTempD2 = (real8 *)calloc(1, numRows * sizeof(real8));

        CSpline(pressure, meltTemp, meltTempD2, numRows);
        CSplint(pressure, meltTemp, meltTempD2, numRows, targetPressure, &param->meltTemp);

        // Be sure to clean up any temporary arrays before returning.

        free(columnData[0]);
        free(columnData[1]);
        free(columnData);
        free(meltTempD2);

        return;
}

//---------------------------------------------------------------------------
//
//      Function:     GetNonLinearMobConstants
//      Description:  Read in the table of material specific
//                    peierls, deltah_0, and alpha values (over
//                    a range of temperatures and pressures) used
//                    in the non-linear mobility calculations, and
//                    do some interpolation to get the proper
//                    values for each of those constants based
//                    on the provided pressure and temperature.
//
//            The data is formatted as:
//
//              <pressure1> <temp1> <deltaH_0> <peierls> <alpha>
//              <pressure2> <temp1> <deltaH_0> <peierls> <alpha>
//              ...
//              <pressureN> <temp1> <deltaH_0> <peierls> <alpha>
//              <pressure1> <temp2> <deltaH_0> <peierls> <alpha>
//              <pressure2> <temp2> <deltaH_0> <peierls> <alpha>
//              ...
//              <pressureN> <tempM> <deltaH_0> <peierls> <alpha>
//
//---------------------------------------------------------------------------

void GetNonLinearMobConstants(Param_t *param)
{
        real8   pressureInGPa = param->pressure / 1.0e+09;
        char   *fileName      = param->MobFileNL;

         // Read the entire table into some internal arrays.

        int     numRows = 0;
        int     numCols = 5;
        real8 **columnData = (real8 **)NULL;

        ReadTabulatedData(fileName, numCols, &columnData, &numRows);

        real8 *p        = columnData[0];   // pressures
        real8 *t        = columnData[1];   // temps
        real8 *deltah0  = columnData[2];   // deltah0s
        real8 *peierls  = columnData[3];   // peierls
        real8 *alpha    = columnData[4];   // alphas

        // Use a bilinear interpolator to calculate the necessary
        // constants for the provided temperature and pressure.

        int status = InterpolateBL(t, p, peierls, numRows, param->TempK, pressureInGPa, &param->MobPeierls);

        if (status != 1) {
            Fatal("Error interpolating peierls value from table\n"
                  "%s\n    for pressure = %e, temperature = %e",
                  fileName, pressureInGPa, param->TempK);
        }


        status = InterpolateBL(t, p, deltah0, numRows, param->TempK,
                               pressureInGPa, &param->MobDeltaH0);
        if (status != 1) {
            Fatal("Error interpolating deltaH0 value from table\n"
                  "%s\n    for pressure = %e, temperature = %e",
                  fileName, pressureInGPa, param->TempK);
        }


        status = InterpolateBL(t, p, alpha, numRows, param->TempK,
                               pressureInGPa, &param->MobAlpha);
        if (status != 1) {
            Fatal("Error interpolating alpha value from table\n"
                  "%s\n    for pressure = %e, temperature = %e",
                  fileName, pressureInGPa, param->TempK);
        }

        // Convert things to the proper units:
        //
        //       peierls from Gpa to Pa
        // and   deltah0 from ev to ?

        param->MobPeierls *= 1.0e+09;
        param->MobDeltaH0 *= 1.60217653e-19;

        printf(" **   Setting peierls to %e\n", param->MobPeierls);
        printf(" **   Setting deltaH0 to %e\n", param->MobDeltaH0);
        printf(" **   Setting alpha   to %e\n", param->MobAlpha  );

        // And don't forget to free up the temporary arrays...

        for (int i=0; (i<numCols); i++) {
            free(columnData[i]);
        }

        free(columnData);

        return;
}

//---------------------------------------------------------------------------
//
//      Function:     SetMaterialTypeFromMobility
//      Description:  Based on the selected mobility module, set up some
//                    basic mobility parameters, and set the material
//                    type appropriately.  Most mobility modules are
//                    tied to specific mterial types, but some (the
//                    relaxation mobilities for instance) are not material
//                    specific.
//
//---------------------------------------------------------------------------

void SetMaterialTypeFromMobility(Home_t *home)
{
        Param_t *param    = home->param;
        int      mobIndex = 0;

        // Look up the attributes of the user-specified mobility module
        // and set various mobility parameters appropriately.

        for (mobIndex=0; mobIndex < MOB_MAX_INDEX; mobIndex++)
        {
            int mat=0;

            if (StrEquiv(param->mobilityLaw, (char *)mobAttrList[mobIndex].mobName))
            {
                param->mobilityType     = mobAttrList[mobIndex].mobilityType;
                param->mobilityInitFunc = mobAttrList[mobIndex].mobInitFunc;
                param->mobilityFunc     = mobAttrList[mobIndex].mobFunc;

                // If the mobility module is planar and glide planes are
                // not being enforced, explicitly enable the enforcement
                // and warn the user.

                if ((mobAttrList[mobIndex].isPlanarMobility) &&
                    (param->enforceGlidePlanes == 0)) {

                    param->enforceGlidePlanes = 1;

                    if (home->myDomain == 0) {
                        printf("The specified mobility (%s) requires the "
                               "enforceGlidePlanes\ncontrol parameter "
                               "toggle to be set.  Enabling toggle now.\n",
                               param->mobilityLaw);
                    }
                }

                if ((mobAttrList[mobIndex].isPlanarMobility == 0) &&
                    (param->enforceGlidePlanes == 1)) {

                    param->enforceGlidePlanes = 0;

                    if (home->myDomain == 0) {
                        printf("The specified mobility (%s) does not allow the "
                               "enforceGlidePlanes feature.\n So the control parameter "
                               "should be set to be zero.\n",
                               param->mobilityLaw);
                    }
                }

                param->allowFuzzyGlidePlanes = mobAttrList[mobIndex].allowFuzzyPlanes;

                // Map the user-specified material type name (if any) to
                // the corresponding integer value.  If the user did not
                // specify a material type, use the default type for the
                // mobility module selected.

                if (StrEquiv(param->materialTypeName, MAT_TYPE_NAME_UNDEFINED)){
                    int matType = mobAttrList[mobIndex].defaultMatType;
                    for (mat = 0; mat < MAT_TYPE_MAX_INDEX; mat++) {
                        if (matType == matInfo[mat].materialType) {
                            param->materialType = matType;
                            param->numBurgGroups = matInfo[mat].numBurgGroups;
                            strcpy(param->materialTypeName,
                                   matInfo[mat].materialTypeName);
                            break;
                        }
                    }
                } else {
                    for (mat = 0; mat < MAT_TYPE_MAX_INDEX; mat++) {
                        if (StrEquiv(param->materialTypeName,
                                     (char *)matInfo[mat].materialTypeName)) {

                            param->materialType = matInfo[mat].materialType;
                            param->numBurgGroups = matInfo[mat].numBurgGroups;
                            break;
                        }
                    }
                }

                if (mat >= MAT_TYPE_MAX_INDEX) {
                    Fatal("Specified <materialTypeName> (%s) is invalid.", param->materialTypeName);
                }

                // If the user-specified material type does not match the
                // default material type for the mobility module, we *may*
                // want to warn the user.

                if (param->materialType!=mobAttrList[mobIndex].defaultMatType) {
                    if (mobAttrList[mobIndex].warnOnMatTypeMismatch) {
                        printf("Warning: <materialType=%s> does not match the"
                               "\n         default material type for "
                               "<mobilityLaw=%s>\n", matInfo[mat].materialTypeName, param->mobilityLaw);
                    }
                }

                param->mobilityIndex = mobIndex;
                break;
            }
        }

        if (mobIndex >= MOB_MAX_INDEX) {
            Fatal("Unknown mobility function %s", param->mobilityLaw);
        }

        return;
}

//---------------------------------------------------------------------------
//
//      Function:     SetMobilityParameters
//      Description:  Calculate the values for all mobility related
//                    functions that are dependent on other values and
//                    have not yet been set elsewhere.  Assumption is
//                    that all the independent mobilityt parameters
//                    have already been set prior to the invocation
//                    of this function.
//
//---------------------------------------------------------------------------

static void SetMobilityParameters(Home_t *home)
{
        // Some of the mobility functions use values derived from some
        // of the input parameters.  For these cases, invoke mobility-
        // specific functions to initialize these values for subsequent
        // calls to the mobility functions.

        if (home->param->mobilityInitFunc != NULL) {
            home->param->mobilityInitFunc(home);
        }

        return;
}

//---------------------------------------------------------------------------
//
//      Function:     SetRemainingDefaults
//      Description:  The default values of certain global parameters
//                    are special in that they depend on values of
//                    other global parameters.  If the user did not
//                    specify values for these special parameters,
//                    this function will calculate the necessary
//                    defaults (as well as do some additional sanity
//                    checks on some of the values).
//
//---------------------------------------------------------------------------

void SetRemainingDefaults(Home_t *home)
{
        real8   tmp, eps;
        real8   xCellSize, yCellSize, zCellSize, minCellSize;
        Param_t *param;

        param = home->param;

        param->delSegLength = 0.0;

        xCellSize = param->Lx / param->nXcells;
        yCellSize = param->Ly / param->nYcells;
        zCellSize = param->Lz / param->nZcells;

        minCellSize = MIN(xCellSize, yCellSize);
        minCellSize = MIN(minCellSize, zCellSize);

        eps = 1.0e-02;

        // The core radius and maximum segment length are required
        // inputs.  If the user did not provide both values, abort
        // now.

        if (home->myDomain == 0) {
            if (param->rc < 0.0) {
                Fatal("The <rc> parameter is required but was not \n"
                      "    provided in the control file");
            }

            if (param->maxSeg < 0.0) {
                Fatal("The <maxSeg> parameter is required but was not \n"
                      "    provided in the control file");
            }
        }

        // If not provided, set position error tolerance based on <rc>

        if (param->rTol <= 0.0) {
            param->rTol = 0.25 * param->rc;
        }

        // The deltaTT is set in the timestep integrator, but some
        // mobility functions now use the deltaTT value, so it must
        // be initialized before ParadisStep() is called since there
        // is an initial mobility calculation done *before* timestep
        // integration the first time into the function.

        param->deltaTT = MIN(param->maxDT, param->nextDT);

        if (param->deltaTT <= 0.0) {
            param->deltaTT = param->maxDT;
        }

        // If the annihilation distance has not been provided, set
        // it to a default based on <rc>

        if (param->rann <= 0.0) {
            param->rann = 2.0 * param->rTol;
        }

        // Minimum area criteria for remesh is dependent on maximum
        // and minumum segment lengths and position error tolerance.

        param->remeshAreaMin = 2.0 * param->rTol * param->maxSeg;

        if (param->minSeg > 0.0) {
            param->remeshAreaMin = MIN(param->remeshAreaMin,
                                       (param->minSeg * param->minSeg *
                                        sqrt(3.0) / 4));
        }

        // Maximum area criteria for remesh is dependent on minimum area,
        // and maximum segment length.

        param->remeshAreaMax = 0.5 * ((4.0 * param->remeshAreaMin) +
                                      (0.25 * sqrt(3)) *
                                      (param->maxSeg*param->maxSeg));

        // If the user did not provide a minSeg length, calculate one
        // based on the remesh minimum area criteria.

        if (param->minSeg <= 0.0) {
            param->minSeg = sqrt(param->remeshAreaMin * (4.0 / sqrt(3)));
        }

        // If the user did not provide an Ecore value, set the default
        // based on the shear modulus and rc values

        if (param->Ecore < 0.0) {
            param->Ecore = (param->shearModulus / (4*M_PI)) *
                           log(param->rc/0.1);
        }

#ifdef ESHELBY
        if (param->shearModulus2 < 0.0) param->shearModulus2 =  param->shearModulus;
        if (param->pois2         < 0.0) param->pois2         =  param->pois;
        if (param->Ecore2        < 0.0) param->Ecore2        = (param->shearModulus2/(4*M_PI)) * log(param->rc/0.1);
#endif
        // Now do some additional sanity checks.

        if (home->myDomain == 0) {

        //     First check for some fatal errors...

            if (param->maxSeg <= param->rTol * (32.0 / sqrt(3.0))) {
                Fatal("Maximum segment length must be > rTol * 32 / sqrt(3)\n"
                      "    Current maxSeg = %lf, rTol = %lf",
                      param->maxSeg, param->rTol);
            }

            if (param->minSeg > (0.5 * param->maxSeg)) {
                Fatal("Minimum segment length must be < (0.5 * maxSeg)\n"
                      "    Current minSeg = %lf, maxSeg = %lf",
                      param->minSeg, param->maxSeg);
            }

            if (param->maxSeg <= param->minSeg) {
                Fatal("Max segment length (%e) must be greater than the\n"
                      "    minimum segment length (%e)", param->maxSeg,
                      param->minSeg);
            }

            if (param->maxSeg > (minCellSize * 0.9)) {
                Fatal("The maxSeg length must be less than the "
                      "minimum cell size * 0.9.  Current values:\n"
                      "    maxSeg    = %.1f\n    cellWidth = %.1f",
                      param->maxSeg, minCellSize);
            }

            if (param->remeshAreaMin > (0.25 * param->remeshAreaMax)) {
                Fatal("remeshAreaMin must be less than 0.25*remeshAreaMax\n"
                      "    Current remeshAreaMin = %lf, remeshAreaMax = %lf",
                      param->remeshAreaMin, param->remeshAreaMax);
            }

            // Now check for conditions that although not fatal, may result
            // in undesired behaviour, and warn the user.

            if (param->rc < 0.1) {
                fprintf(stderr, "WARNING: specified rc value (%e) will yield a \nnegative core energy\n", param->rc);
            }

            tmp = (param->maxSeg * param->maxSeg * param->maxSeg);

            if (param->remeshAreaMax > (0.25 * sqrt(3) * tmp)) {
                fprintf(stderr, "WARNING: Area criteria will be unused in remesh operations!\n");
                fprintf(stderr, "         rmeshAreaMax = %lf, maxSeg = %lf\n", param->remeshAreaMax, param->maxSeg);
            }

            if (param->rann > (0.5 * param->rc + eps)) {
                fprintf(stderr, "WARNING: Separation distance is larger than the core radius!\n");
                fprintf(stderr, "         rann = %lf, rc = %lf\n", param->rann, param->rc);
            }

            if (param->rann > (2.0 * param->rTol)) {
                fprintf(stderr, "WARNING: Collision distance is outside the positioning error tolerance!\n");
                fprintf(stderr, "         rann = %lf, rTol = %lf\n", param->rann, param->rTol);
            }

#if 0
            tmp = param->remeshAreaMin - (2.0 * param->rTol * param->maxSeg);

            if (fabs(tmp) > eps) {
                fprintf(stderr, "WARNING: remesh minimum area != 2.0 * rTol * maxSeg\n");
                fprintf(stderr, "         remeshAreaMin = %lf, rTol = %lf maxSeg = %lf\n", param->remeshAreaMin, param->rTol, param->maxSeg);
            }
#endif

            // If free suraces are used but the specified surfaces are
            // not within the primary bounding box, it's a problem.

            if (((param->xBoundType == Free) &&
                 ((param->xBoundMin < param->minCoordinates[X]) ||
                  (param->xBoundMax > param->maxCoordinates[X]))) ||

                ((param->yBoundType == Free) &&
                 ((param->yBoundMin < param->minCoordinates[Y]) ||
                  (param->yBoundMax > param->maxCoordinates[Y]))) ||

                ((param->zBoundType == Free) &&
                 ((param->zBoundMin < param->minCoordinates[Z]) ||
                  (param->zBoundMax > param->maxCoordinates[Z])))) {
                Fatal("Free surfaces are not within main bounding box!\n"
                      "    Surface min coordinates (%lf %lf %lf)\n"
                      "    Surface max coordinates (%lf %lf %lf)\n",
                      param->xBoundMin, param->yBoundMin, param->zBoundMin,
                      param->xBoundMax, param->yBoundMax, param->zBoundMax);
            }

#ifndef _ARLFEM
            // If free surfaces are enabled but the finite element code
            // is not linked in, results will not be accurate, so print
            // a warning.

            if ((param->xBoundType == Free) ||
                (param->yBoundType == Free) ||
                (param->zBoundType == Free)) {
                printf("***\n"
                       "*** WARNING!  Use of free surfaces in ParaDiS without the\n"
                       "*** FEM/ParaDiS coupling is not fully supported!\n"
                       "***\n" );
            }
#endif  // ifndef _ARLFEM

        }  // if domain == 0


        // If there are a mix of free surfaces and periodic boundaries,
        // the boundary min/max values must default to the simulation
        // boundaries in the dimensions without free surfaces.

        if ((param->xBoundType == Free) ||
            (param->yBoundType == Free) ||
            (param->zBoundType == Free)) {

            if (param->xBoundType == Periodic) {
                param->xBoundMin = param->minSideX;
                param->xBoundMax = param->maxSideX;
            }
            if (param->yBoundType == Periodic) {
                param->yBoundMin = param->minSideY;
                param->yBoundMax = param->maxSideY;
            }
            if (param->zBoundType == Periodic) {
                param->zBoundMin = param->minSideZ;
                param->zBoundMax = param->maxSideZ;
            }
        }


        param->partialDisloDensity = (real8 *)malloc(param->numBurgGroups * sizeof(real8));

        // If type 1 domain decompositionis enabled, the DLBfreq
        // value MUST be a multiple of 3.

        if ((param->DLBfreq > 0) && (param->decompType == 1)) {
            param->DLBfreq = (param->DLBfreq + 2) / 3 * 3;
        }

        // Some of the mobility related parameters are dependent on others.
        // At this point all the independent parameters should be set, so
        // we can calculate the dependent values now.

        SetMobilityParameters(home);

#ifdef ESHELBY
        // Set a flag if the user provided a file of eshelby inclusions
        // Allows to run regular ParaDiS with ESHELBY flag on and no inclusions
        // are present.

        param->enableInclusions = (param->inclusionFile[0] != 0);
#endif

        // If the cross slip flag has not been explicitly defined, give
        // it a default setting based on the mobility being used.

        if (param->enableCrossSlip < 0) {
            if (mobAttrList[param->mobilityIndex].doCrossSlip) {
                    param->enableCrossSlip = 1;
            } else {
                    param->enableCrossSlip = 0;
            }
        }

        // Several portions of the code need to calculate the total
        // volume of the simulation and a volume factor used for
        // converting dislocation length to density, so set those
        // factors now.

        // If the FEM code is linked in, this is going to depend on the
        // volume of the finite domain which is computed in the FEM and
        // sent to ParaDiS in ParadisInit.  We use the values of
        // x,y,z-bound as a placeholder.

#ifdef _ARLFEM
        param->simVol = (param->xBoundMax-param->xBoundMin) *
                        (param->yBoundMax-param->yBoundMin) *
                        (param->zBoundMax-param->zBoundMin);
#else
        if ((param->zBoundType == Free) ||
            (param->yBoundType == Free) ||
            (param->xBoundType == Free)) {

            param->simVol = (param->xBoundMax-param->xBoundMin) *
                            (param->yBoundMax-param->yBoundMin) *
                            (param->zBoundMax-param->zBoundMin);
        } else {
            param->simVol = param->Lx * param->Ly * param->Lz;
        }
#endif  // ifdef _ARLFEM

        param->burgVolFactor = 1.0 / (param->burgMag * param->burgMag *
                                      param->simVol);

#ifdef FULL_N2_FORCES
        // To do full n^2 force calculations without remote forces, we need
        // to explicitly set some flags.

        param->fmEnabled        = 0;
#ifdef ESHELBY
        param->eshelbyfmEnabled = 0;
#endif
        param->numDLBCycles     = 0;
        param->DLBfreq          = 0;
#endif

#ifdef SPECTRAL
        if (param->FFTenabled) {
            param->fmEnabled = 0;
        }
#endif

        // For HCP simulations, the Ecore value is a function of the burgers
        // vector type.  These Ecore values can not be calculated dynamically
        // which meanss the caller MUST have provided them in the cotrol
        // parameter file.  If they have not been provided, we have a problem.

        if (param->materialType == MAT_TYPE_HCP) {
            if ((param->HCPEcoreA   < 0.0) ||
                (param->HCPEcoreC   < 0.0) ||
                (param->HCPEcoreCpA < 0.0)) {
                const char *errMsg = "\n\nFor HCP simulations the HCPEcoreA, HCPEcoreC, and HCPEcoreCpA"
                                     "\nvalues MUST be provided in the control parameter file, but at"
                                     "\nleast one them was not.\n";
                Fatal(errMsg);
            }
        }

        return;
}

//---------------------------------------------------------------------------
//
//      Function:     VerifyConstraints
//      Description:  Do some validity checking on the nodal constraint
//                    values and either print a warning or abort depending
//                    on the problem.
//
//---------------------------------------------------------------------------

static void VerifyConstraints(Home_t *home)
{
        int allValidConstraints;
        int localError[2], globalError[2];

        // localError represents two types of "errors".  The first
        // component represents an error due to a completely invalid
        // constraint.  The second component represents an attempt to
        // pin a node in at least 1 but less than 3 dimensions (which
        // generates a warning rather than a fatalerror).

        localError[0]  = 0;
        localError[1]  = 0;

        globalError[0] = 0;
        globalError[1] = 0;

        for (int i=0; i < home->newNodeKeyPtr; i++)
        {
            Node_t *node;

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            // If a node has a completely unrecognized constraint, it's
            // an error worthy of an abort.

            // If the node has any constraints that are not in the list
            // of all valid constraints, it's a problem.

            allValidConstraints = PINNED_X | PINNED_Y | PINNED_Z |
                                  SURFACE_NODE | MULTI_JUNCTION_NODE |
                                  LOOP_NODE | CROSS_NODE;

            if (HAS_ANY_OF_CONSTRAINTS(node->constraint,
                                       ~allValidConstraints))
            {
                localError[0] = node->constraint;
                break;
            }

            // Current code does not fully support pinning a node in
            // only specific dimensions, so if we find a node that is
            // pinned but not in all dimensions, issue a warning.

            if (HAS_ANY_OF_CONSTRAINTS(node->constraint, PINNED_NODE)) {
                if (!HAS_ALL_OF_CONSTRAINTS(node->constraint, PINNED_NODE)) {
                    localError[1] = 1;
                    break;
                }
            }

        }

        // Domain zero gathers info on any errors and advises or aborts
        // as necessary.

#ifdef PARALLEL
        MPI_Allreduce(localError, globalError, 2, MPI_INT,
                      MPI_MAX, MPI_COMM_WORLD);
#else
        globalError[0] = localError[0];
        globalError[1] = localError[1];
#endif

        if (home->myDomain == 0)
        {
            // An invalid constraint is an error...

            if (globalError[0] != 0) {
                Fatal("An invalid nodal constraint (%d) was found in"
                      "the nodal data file.", globalError[0]);
            }

            // If a node was pinned in at least 1 but not all dimensions
            // it just requires a warning...

            if (globalError[1] != 0) {
                fprintf(stderr,
                        "WARNING: Currently ParaDiS only provides full "
                        "support for pinned nodes\n         if they are "
                        "pinned in All dimensions.  At least 1 node has "
                        "a\n         constraint pinning it in some but "
                        "not all dimensions!\n\n");
            }
        }

        return;
}

//---------------------------------------------------------------------------
//
//      Function:     CheckForGlidePlanes
//      Description:  Some of the mobility modules enforce the motion
//                    of dislocation segments along specific glide planes.
//                    If such a mobility is being used, verify that a
//                    glide plane has been specified for each segment and
//                    abort with an error if there are any segments without.
//
//---------------------------------------------------------------------------

static void CheckForGlidePlanes(Home_t *home)
{
        int     haveGlidePlanesLocal  = 1;
        int     haveGlidePlanesGlobal = 0;
        real8   eps = 1.0e-03;

        Param_t *param = home->param;

        // Unless we're dealing with one of the mobility laws that
        // requires glide planes to be defined for all segments,
        // nothing needs to be done here so just return.

        if (param->enforceGlidePlanes == 0) {
            return;
        }

        // Loop through every node and every segment attached to the
        // node and look for zeroed glide planes.

        int maxNodeIndex = home->newNodeKeyPtr;

#ifdef _OPENMP
#pragma omp parallel for shared(haveGlidePlanesLocal, eps, home)
#endif
        for (int i=0; i<maxNodeIndex; i++)
        {
            Node_t  *node;

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            for (int j=0; j < node->numNbrs; j++)
            {
                real8 tmp = (node->nx[j] * node->nx[j]) +
                            (node->ny[j] * node->ny[j]) +
                            (node->nz[j] * node->nz[j]);

                // <haveGlidePlanesLocal> is shared among all threads, but since
                // the only reference is to turn the flag off on error, there's
                // no contention problems.

                if (tmp < eps) {
                    haveGlidePlanesLocal = 0;
                    break;
                }
            }
        }

        // Find out if any procesor located a segment without a glide plane

#ifdef PARALLEL
        MPI_Allreduce(&haveGlidePlanesLocal, &haveGlidePlanesGlobal, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#else
        haveGlidePlanesGlobal = haveGlidePlanesLocal;
#endif

        // If there were any segments with specified glide planes
        // have processor zero print an error message and force
        // an abort.

        if ((home->myDomain == 0) && (haveGlidePlanesGlobal == 0)) {
            Fatal("The selected mobility law (%s) requires glide\n"
                  "       planes to be defined for all segments.  One or\n"
                  "       more segments in the provided nodal data file\n"
                  "       do not have glide planes specified.",
                  param->mobilityLaw);
        }

        return;
}

//---------------------------------------------------------------------------
//
//      Author:       Gregg Hommes
//
//      Function:     VerifyBurgersVectors
//
//      Description:  This function does a simple sanity check during
//                    initialization to verify that the burgers vector
//                    at one end of a segment matches the burgers vector
//                    at the other end.  Note: This function assumes the
//                    local domain has nodal info for all nodes terminating
//                    segments attached to local nodes.  This means that
//                    the ghost node information must be available before
//                    this function is called.
//
//      Last Modified:  01/09/08: - original version
//
//---------------------------------------------------------------------------

void VerifyBurgersVectors(Home_t *home)
{
        int      nodeID, maxNodeIndex;
        int      errorType = 0;
        int      errNodeSegID, errNbrNodeSegID;
        real8    eps = 1.0e-03;
        Node_t  *errNode, *errNbrNode;
        Param_t *param = home->param;

        // Loop over all local nodes.  If at any point the error flag gets
        // set, we can't break out of the loop due to threading issues, but
        // we can just skip any more checking in each loop iteration.

        maxNodeIndex = home->newNodeKeyPtr;

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (nodeID = 0; nodeID < maxNodeIndex; nodeID++) {
            int    armID, nbrArmID;
            real8  burgSumX, burgSumY, burgSumZ;
            Node_t *node, *nbr;

            if (errorType != 0) {
                continue;
            }

            node = home->nodeKeys[nodeID];

            if (node == (Node_t *)NULL) {
                continue;
            }

            // Loop over every segment attached to the local node

            for (armID = 0; armID < node->numNbrs; armID++) {

                if (errorType != 0) {
                    continue;
                }

                nbr = GetNeighborNode(home, node, armID);

                if (nbr == (Node_t *)NULL) {
#ifdef _OPENMP
#pragma omp critical (VERIFY_BURG_ERR)
#endif
                    {
                        errorType = 1;
                        errNode = node;
                        errNodeSegID = armID;
                    }  // end omp critical section
                }

                // Find the neighbor's arm that connects back to the current
                // node and get its index

                for (nbrArmID = 0; nbrArmID < nbr->numNbrs; nbrArmID++) {
                    if ((nbr->nbrTag[nbrArmID].domainID == home->myDomain) &&
                        (nbr->nbrTag[nbrArmID].index == node->myTag.index)) {
                        break;
                    }
                }

                if (nbrArmID >= nbr->numNbrs) {
#ifdef _OPENMP
#pragma omp critical (VERIFY_BURG_ERR)
#endif
                    {
                        errorType = 2;
                        errNode = node;
                        errNbrNode = nbr;
                    }  // end omp critical section
                }

                // If the sum of any of the corresponding components of the
                // burgers vectors at the two ends of the segment are not
                // zero, we have a problem.

                burgSumX = node->burgX[armID] + nbr->burgX[nbrArmID];
                burgSumY = node->burgY[armID] + nbr->burgY[nbrArmID];
                burgSumZ = node->burgZ[armID] + nbr->burgZ[nbrArmID];

                if ((fabs(burgSumX) > eps) ||
                    (fabs(burgSumY) > eps) ||
                    (fabs(burgSumZ) > eps)) {
#ifdef _OPENMP
#pragma omp critical (VERIFY_BURG_ERR)
#endif
                    {
                        errorType = 3;
                        errNode = node;
                        errNbrNode = nbr;
                        errNodeSegID = armID;
                        errNbrNodeSegID = nbrArmID;
                    }  // end omp critical section
                }

                // For HCP dislocations, the burgers vectors must match
                // the non-normalized burgers vectors used internal to the
                // code (calculated dynamically in GetBurgListHCP()) or
                // things (like mobility functions) may not work properly.
                // We won't abort on a mismatch but warnings will be issued.

                // So, rotate the specified burgers vector into the crystal
                // frame (if needed) and verify that it matches one of the
                // known burgers vectors.

                if (param->materialType == MAT_TYPE_HCP) {
                    int    burgIndex;
                    real8  bCryst[3];

                    if (param->useLabFrame) {
                        real8 bLab[3];
                        bLab[X] = node->burgX[armID];
                        bLab[Y] = node->burgY[armID];
                        bLab[Z] = node->burgZ[armID];
                        Matrix33Vector3Multiply(param->rotMatrixInverse,
                                                bLab, bCryst);
                    } else {
                        bCryst[X] = node->burgX[armID];
                        bCryst[Y] = node->burgY[armID];
                        bCryst[Z] = node->burgZ[armID];
                    }

                    GetBurgIndexHCP(bCryst, 0, HCP_NUM_TOT_BURG,
                                    home->burgData.burgList, &burgIndex);

                    if (burgIndex < 0) {
                        printf("WARNING: Unexpected HCP burgers vector for "
                               "node (%d,%d) segment %d.\n         "
                               "Unexpected burgers vector = %.6e %.6e %.6e\n",
                               node->myTag.domainID, node->myTag.index, armID,
                               node->burgX[armID], node->burgY[armID],
                               node->burgZ[armID]);

                    }

                }  // end if materialType == HCP

            }  // end for (armID = 0; ... )

        }  // end for (nodeID = 0; ...)

        // If there were any errors found above, abort with an appropriate
        // error message.

        switch (errorType) {
            case 1:
                Fatal("VerifyBurgersVectors(): Lookup of node (%d,%d) failed!",
                      errNode->nbrTag[errNodeSegID].domainID,
                      errNode->nbrTag[errNodeSegID].index);
                break;

            case 2:
                Fatal("VerifyBurgersVectors(): neighbor node (%d,%d) "
                      "not linked back\n    to local node (%d,%d)",
                      errNbrNode->myTag.domainID, errNbrNode->myTag.index,
                      errNode->myTag.domainID, errNode->myTag.index);
                break;

            case 3:
                Fatal("VerifyBurgersVectors(): Burgers vector mismatch!\n"
                      "    Segment (%d,%d)--(%d,%d)\n"
                      "    burg at first node  = %e %e %e\n"
                      "    burg at second node = %e %e %e\n",
                      errNode->myTag.domainID, errNode->myTag.index,
                      errNbrNode->myTag.domainID, errNbrNode->myTag.index,
                      errNode->burgX[errNodeSegID],
                      errNode->burgY[errNodeSegID],
                      errNode->burgZ[errNodeSegID],
                      errNbrNode->burgX[errNbrNodeSegID],
                      errNbrNode->burgY[errNbrNodeSegID],
                      errNbrNode->burgZ[errNbrNodeSegID]);
                break;

            default:
                break;
        }

        return;
}

//---------------------------------------------------------------------------
//
//      Function:    PrintBanner
//      Description: Prints a banner to stdout with some basic info
//                   about the current execution.
//
//---------------------------------------------------------------------------

static void PrintBanner(Home_t *home, int argc, char *argv[])
{
        time_t         currTime;
        char           *separator=0, *execName=0, *execPath=0;
        char           workingDir [512]; memset(workingDir ,0,sizeof(workingDir ));
        char           tmpExecName[512]; memset(tmpExecName,0,sizeof(tmpExecName));
        char           currTimeStr[128]; memset(currTimeStr,0,sizeof(currTimeStr));
        char           cmdLine    [512]; memset(cmdLine    ,0,sizeof(cmdLine    ));
        char           userName   [ 32]; memset(userName   ,0,sizeof(userName   ));
        struct utsname utsInfo;
        struct passwd  *pwInfo;

        time(&currTime);
        strcpy(currTimeStr, ctime(&currTime));
        currTimeStr[strlen(currTimeStr)-1] = 0;

        strcpy(tmpExecName, argv[0]);
        separator = strrchr(tmpExecName, '/');
        if (separator == (char *)NULL)
        {
                execPath = (char *) ".";
                execName = (char *) tmpExecName;
        } else {
                *separator = 0;
                execName = (char *) separator + 1;
                execPath = (char *) tmpExecName;
        }

#ifndef _XT4
        getcwd(workingDir, sizeof(workingDir) - 1);
        uname(&utsInfo);

        uid_t uid = getuid();
#if !defined _BGL && !defined _BGP && !defined _BGQ
        pwInfo = getpwuid(uid);
        strcpy(userName, pwInfo->pw_name);
#else
        sprintf(userName, "%d", uid);
#endif
#endif  // _XT4

        cmdLine[0] = 0;
        for (int i=1; i<argc; i++) {
            strcat(cmdLine, argv[i]);
            strcat(cmdLine, " ");
        }

        printf("***********************************************************\n");
        printf("**** \n");
        printf("**** Time of run:     %s\n", currTimeStr);
        printf("**** Executable name: %s\n", execName);
        printf("**** Executable path: %s\n", execPath);
#ifndef _XT4
        printf("**** Working dir:     %s\n", workingDir);
        printf("**** Execution host:  %s (%s %s %s)\n", utsInfo.nodename, utsInfo.sysname,
                                                        utsInfo.release, utsInfo.machine);
        printf("**** User name:       %s\n", userName);
#endif
        printf("**** Number of tasks: %d\n", home->numDomains);
        printf("**** Command line:    %s\n", cmdLine);
#ifdef _OPENMP
        printf("**** Thread count:    %d\n", omp_get_max_threads());
#endif
        printf("**** \n");
        printf("**** Build options enabled...\n");

#ifdef NO_XWINDOW
        printf("****    -DNO_XWINDOW         : X-Windows disabled\n");
#endif

#ifdef ANISOTROPIC
        printf("****    -DANISOTROPIC        : anisotropic elasticity enabled \n");
#endif

#ifdef ANISO_QMAX
        printf("****    -DANISO_QMAX         : anisotropic elasticity qmax=%d\n", ANISO_QMAX );
#endif

#ifdef ESHELBY
        printf("****    -DESHELBY            : Eshelby enabled\n");
#endif

#ifdef ESHELBYFORCE
        printf("****    -DESHELBYFORCE       : Eshelby force enabled\n");
#endif

#ifdef ESHELBYCORE
        printf("****    -DESHELBYCORE        : Eshelby core enabled\n");
#endif

#ifdef ESHELBYSTEP
        printf("****    -DESHELBYSTEP        : Eshelby step enabled\n");
#endif

#ifdef SPECTRAL
        printf("****    -DSPECTRAL           : Spectral forces enabled\n");
#endif

#ifdef FIX_PLASTIC_STRAIN
        printf("****    -DFIX_PLASTIC_STRAIN : Plastic strain fix enabled\n");
#endif

#ifdef FRCRIT
        printf("****    -DFRCRIT             : Critical force enabled\n");
#endif

#ifdef CALCENERGY
        printf("****    -DCALCENERGY         : Calculate energy enabled\n");
#endif

#ifdef PRINTSTRESS
        printf("****    -DPRINTSTRESS        : Print stress enabled\n");
#endif

#ifdef TAYLORFMM
        printf("****    -DTAYLORFMM          : Taylor expansion FMM enabled\n");
#endif

#ifdef BBFMM
        printf("****    -DBBFMM              : Chebyshev polynomial expansion FMM enabled\n");
#endif

#ifdef UNIFORMFMM
        printf("****    -DUNIFORMFMM         : Uniform polynomial expansion FMM enabled\n");
#endif

#ifdef USE_KINSOL
        printf("****    -DUSE_KINSOL         : Sundials KINSOL solver available (enable via control params)\n");
#endif

#ifdef USE_ARKODE
        printf("****    -DUSE_ARKODE         : Sundials ARKODE solver available (enable via control params)\n");
#endif

#ifdef GPU_ENABLED
        printf("****\n");
        printf("****    -DGPU_ENABLED        : GPU support available (enable via control params)\n");

        // dump some GPU device info to the log...
        {
           int gpu_cnt = CUDA_Device_Count();

           if (gpu_cnt>0)
           {
              cudaDeviceProp  *props = CUDA_Device_Properties(0);

              if (props)
              {
                 printf("****\n");
                 printf("****        GPU DEVICE INFO...\n");
                 printf("****        --------------------------------------------------------\n"                       );
                 printf("****        gpu device type       : %s\n"           , props->name                             );
                 printf("****        gpu device count      : %d\n"           , gpu_cnt                                 );
                 printf("****        global memory         : %ldMB (%ldGB)\n", props->totalGlobalMem>>20,
                                                                               props->totalGlobalMem>>30               );
                 printf("****        compute mode          : %d%d\n"         , props->major, props->minor              );
                 printf("****        warp size             : %d \n"          , props->warpSize                         );
                 printf("****        max threads per block : %d \n"          , props->maxThreadsPerBlock               );
                 printf("****        unified addressing    : %s\n"           , props->unifiedAddressing ? "yes" : "no" );
                 printf("****        device overlap        : %s\n"           , props->deviceOverlap     ? "yes" : "no" );
                 printf("****        SMU count             : %d\n"           , props->multiProcessorCount              );
                 printf("****        SMU shared mem        : %luK\n"         , props->sharedMemPerMultiprocessor>>10   );
                 printf("****        SMU regs              : %dK\n"          , props->regsPerMultiprocessor     >>10   );
                 printf("****\n");
              }
           }
           else
           { printf("****                    NB! No GPUs are detected on this platform!\n"); }
        }
#endif

        printf("***********************************************************\n");


        return;
}

//---------------------------------------------------------------------------
//
//      Function:     Initialize
//      Description:  This is the driver routine for initialization,
//                    handling some of the general initializations and
//                    calling all the more specific initialization routines.
//
//      Last Modified:  01/09/08: Gregg Hommes - added call to
//                                VerifyBurgersVectors() as a sanity check.
//                      04/25/08: gh  - Added call to read any eshelby
//                                inclusions from an inclusion file.
//
//---------------------------------------------------------------------------

void Initialize
(
   Home_t *home,    ///< pointer to initialized home structure
   int     argc,    ///< application argument count
   char   *argv[]   ///< application argument vector
)
{
        char         tmpDataFile[256], testFile[256];

        char         *ctrlFile=0;
        char         *dataFile=0;
        InData_t     *inData  =0;

        // Allocate some buffers used for MPI reduction operations in
        // various portions of ParaDiS.  NOTE: most of the code that
        // uses these buffers need buffers of <numDomains> elements, but
        // at least 1 portion requires <numDomains>+1 elements, so
        // we allocate for the maximum size.

        home->localReduceBuf  = (int *)calloc(1, (home->numDomains+1) * sizeof(int));
        home->globalReduceBuf = (int *)calloc(1, (home->numDomains+1) * sizeof(int));

        // Initialize map between old and new node tags before
        // reading nodal data and distributing it to the remote
        // domains.

        home->tagMap     = (TagMap_t *)NULL;
        home->tagMapSize = 0;
        home->tagMapEnts = 0;

        // Init defaults for some command-line parameters...

        int   skipIO=0;                  // skip IO reading
        int   numDLBCycles=0;            // initial dynamic load balance cycles
        int   maxNumThreads=1;           // OpenMP thread max
        int   ndx=0, ndy=0, ndz=0;       // domain geometry
        int   ncx=0, ncy=0, ncz=0;       // cell geometry
        int   doBinRead=0;               // disable binary read
        int   gpu_enabled=0;             // disble gpu
        int   maxstep=0;                 // max number of cycles
        char *dirname=0;                 // output directory

        Param_t *param=0;

        if (home->myDomain != 0) {
            param = home->param = (Param_t *)calloc(1, sizeof(Param_t));
        }

        inData = (InData_t *) calloc(1, sizeof(InData_t));

#ifdef _OPENMP
        // Need to explicitly set the default number of threads per task to
        // one.  If the user specifies a different thread count via the
        // command line options, this default will be overridden further on.

        maxNumThreads = 1;
        omp_set_num_threads(maxNumThreads);
#endif

#ifdef _BGQ
        // On the BGQ systems, depending how the system has been configured,
        // querying the system "personality" info for the partition shape/size
        // may not give the numbers for exact set of nodes on which the job is
        // running, so we have to determine where each task is in the partition
        // and determine the size/shape and allocation order ourselves.  Have
        // to do this *before* the ResetTaskGeometry() call because all tasks
        // have to participate in this process.

        int optimalTaskGeometry[3] = {0};

        GetOptimal3DTaskGeometry(home, optimalTaskGeometry);
#endif

        // Verify the command line syntax and pull out the control and
        // data file names (if provided).  If no control file is specified,
        // use a default file name If no dataFile name is provided, use
        // the control file name with the file suffix changed to the
        // appropriate data file name suffix.  Only need to do this on
        // processor zero.

        if (home->myDomain == 0)
        {
            int  ignoreNonLoopFlux = 0;
            char *fluxSelection = (char *)NULL;

            dataFile = (char *) NULL;
            ctrlFile = (char *) "control.script";

            for (int i=1; i<argc; i++)
            {
                if ( strcmp(argv[i], "-r")==0 )
                {
                    if (i >= (argc - 1)) Usage(argv[0]);

                    numDLBCycles = atoi(argv[++i]);

                    if (numDLBCycles <= 0) Usage(argv[0]);
                }
                else if (    (strcmp(argv[i], "-n"    )==0)
                          || (strcmp(argv[i], "-nthrd")==0) )
                {
                    if (i >= (argc - 1)) Usage(argv[0]);
                    maxNumThreads = atoi(argv[++i]);
#ifndef _OPENMP
                    printf("WARNING: Not compiled with OpenMP support; '-n' option ignored.\n");
#else
                    omp_set_num_threads(maxNumThreads);
                    if ((maxNumThreads < 1) ||
                        (maxNumThreads != omp_get_max_threads()))
                    {
                        Fatal("Error: unable to set thread count to %d!\n", maxNumThreads);
                    }
#endif
                }
                else if (    (strcmp(argv[i], "-maxstep" )==0)
                          || (strcmp(argv[i], "-nstep"   )==0) )
                {
                    if (i >= (argc - 1)) Usage(argv[0]);
                    maxstep = atoi(argv[++i]);
                }
                else if (strcmp(argv[i], "-b")==0) { doBinRead = 1; }
                else if (strcmp(argv[i], "-s")==0) { skipIO    = 1; }
                else if (strcmp(argv[i], "-d")==0)
                {
                    if (i >= (argc - 1)) Usage(argv[0]);
                    dataFile = (char *) argv[++i];

                // This argument is provided for some site-specific
                // post-processing and as such is not otherwise
                // documented or supported.  If you don't know what
                // it is for, there's no reason to use it.

                }
                else if (strcmp(argv[i], "-f")==0)
                {
                    if (i >= (argc - 1)) { Fatal("Missing argument for '-f' option"); }
                    fluxSelection = (char *) argv[++i];
                }
                else if (strcmp(argv[i], "-doms")==0)
                {
                    // '-doms' is used to override the domain specification in the control file

                    if (i >= (argc - 1)) { Fatal("Missing arguments for '-dom' option"); }

                    ndx = atoi(argv[++i]);  // (number of domains (x)
                    ndy = atoi(argv[++i]);  // (number of domains (y)
                    ndz = atoi(argv[++i]);  // (number of domains (z)
                }
                else if (strcmp(argv[i], "-cells")==0)
                {
                    // '-cells' is used to override the cell specification in the control file

                    if (i >= (argc - 1)) { Fatal("Missing arguments for '-cell' option"); }

                    ncx = atoi(argv[++i]);  // (number of cells (x)
                    ncy = atoi(argv[++i]);  // (number of cells (y)
                    ncz = atoi(argv[++i]);  // (number of cells (z)
                }
                else if (strcmp(argv[i], "-dir")==0)
                {
                    // '-dir' is used to override the save directory specification in the control file

                    if (i >= (argc - 1)) { Fatal("Missing arguments for '-dir' option"); }

                    dirname = (char *) argv[++i];   // path to destination directory
                }
                else if (strcmp(argv[i], "-gpu")==0)
                {
#ifdef GPU_ENABLED
                    gpu_enabled=1;
#else
                    printf("WARNING: Not compiled with GPU support; '-gpu' option ignored.\n");
                    gpu_enabled=0;
#endif
                }
                else
                {
                    if (i < (argc - 1)) Usage(argv[0]);
                    ctrlFile = (char *) argv[i];
                }
            }

            // The fluxSelection is related to some site-specific
            // post-processing only and as such should never be set
            // in *normal* ParaDiS simulations

            if (fluxSelection != (char *)NULL)
            {
                if (StrEquiv(fluxSelection, "loops")) { ignoreNonLoopFlux=1; }
                else                                  { Fatal("Unrecognized argument for '-f' option"); }
            }

            // If the user did not specify a data file name, set
            // a default name based on the base control file name.

            if (dataFile == (char *)NULL)
            {
                strcpy(tmpDataFile, ctrlFile);

                char *start = strrchr(tmpDataFile, '/');

                if (start == (char *)NULL) start = tmpDataFile;

                char *sep = strrchr(start, '.');

                if ((sep != (char *)NULL) && (sep > start)) *sep = 0;

                // If the user specified that a binary data file is
                // used set default name accordingly.  If the user
                // did not specify, try to figure out whether to
                // use a binary restart or not.

                if (doBinRead)
                {
                    strcat(tmpDataFile, HDF_DATA_FILE_SUFFIX);  // append ".hdf"
                }
                else
                {
                    strcpy(testFile, tmpDataFile);
                    strcat(testFile, HDF_DATA_FILE_SUFFIX);     // append ".hdf"

                    int fd = open(testFile, O_RDONLY, 0);

                    if (fd<0)
                    {
                        strcat(testFile, ".0");              // append ".0"
                        fd = open(testFile, O_RDONLY, 0);    // (try again)
                    }

                    if (fd>=0)
                    {
                        doBinRead=1;
                        close(fd);
                        strcat(tmpDataFile, HDF_DATA_FILE_SUFFIX);  // append ".hdf"
                    }
                    else
                    {
                        strcat(tmpDataFile, NODEDATA_FILE_SUFFIX);  // append ".data"
                    }

                }

                dataFile = tmpDataFile;

            }
            else
            {
                // User provided a data file name, but didn't indicate
                // if it was a binary data file or not.  Make a guess
                // based on the file name.

                if (strstr(dataFile, HDF_DATA_FILE_SUFFIX) != (char *)NULL)
                {
                    doBinRead = 1;
                }
            }

            PrintBanner(home, argc, argv);

            home->ctrlParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
            home->dataParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
            home->param         = (Param_t     *)calloc(1,sizeof(Param_t    ));

            param = home->param;

            param->doBinRead = doBinRead;
            param->maxNumThreads = maxNumThreads;
            param->ignoreNonLoopFlux = ignoreNonLoopFlux;

            CtrlParamInit(param, home->ctrlParamList);
            DataParamInit(param, home->dataParamList);

#ifdef _BGQ
            // Now that the param struct has been allocated we need to
            // store the optimal task geometry for later use.

            VECTOR_COPY(param->optimalTaskGeometry, optimalTaskGeometry);
#endif

#ifdef USE_SCR
            // See if there is a valid SCR checkpoint/restart file
            // from which to restart.  If SCR is enabled and a valid
            // checkpoint file is found, that checkpoint overrides any
            // restart file name provided on the command line.

            if (SCR_Route_file("ctrl.0", testFile) == SCR_SUCCESS)
            {
                ctrlFile = testFile;
            }
#endif

            // Read runtime parameters from the control file.

            printf("Initialize: Parsing control file %s\n", ctrlFile);
            ReadControlFile(home, ctrlFile);
            printf("Initialize: Control file parsing complete\n");

            // Override control file specifications if provided on command line args...

            if (param && (maxstep>0))
            {
               param->maxstep = maxstep;
            }

            if (param && ((ndx*ndy*ndz)>0))
            {
               param->nXdoms = ndx;
               param->nYdoms = ndy;
               param->nZdoms = ndz;
            }

            if (param && ((ncx*ncy*ncz)>0))
            {
               param->nXcells = ncx;
               param->nYcells = ncy;
               param->nZcells = ncz;
            }

            if (param && dirname)
            {
               strncpy(param->dirname,dirname,sizeof(param->dirname));
            }

#ifdef GPU_ENABLED
            if (param && gpu_enabled) { param->gpu_enabled=1; }

            if (param)
            { printf("GPU-accelerated SegSeg Forces are %s\n", ( param->gpu_enabled ? "ENABLED" : "DISABLED" ) ); }
#endif
            // Set number of initial load-balancing-only cycles. Will
            // be zero unless user provided value via "-r" command line
            // option.

            param->numDLBCycles = numDLBCycles;

            // If user explicitly requested major IO to be skipped via
            // command line options, override whatever happened to be
            // in the control file.

            if (skipIO) param->skipIO = 1;

            // Some checks on input consistency

            InputSanity(home);

#if defined _BGP || defined _BGQ
            // On some systems (BlueGene for example) jobs are executed on a
            // physical hardware partition configured as a multi-dimensional
            // grid.  This partition may be dynamically allocated, hence the
            // partition geometry is not known until run time.  We need to
            // determine the geometry of the hardware partition and (if needed
            // and permitted) reset the user-specified domain geometry to match
            // the phsical hardware.

            ResetTaskGeometry(home);
#endif
            // Some of the control parameters (particularly the
            // material specific constants) *may* be interpolated
            // from user-provided tables.  Check if the user
            // has provided tables for these types of parameters
            // and handle them now.

            // First see if the user provided a base name for
            // a set of files with material specific elastic constants,
            // find the constants for the current temp & pressure, and
            // recalculate the poisson ratio, shear modulus and burgers
            // vector magnitude.

            if (param->ecFile[0] != 0)
            {
                (void) GetElasticConstants(param);

                printf(" ** For temp = %6.2fK, pressure = %ePa:\n", param->TempK, param->pressure);
                printf(" **   Setting poisson ratio to %e\n", param->pois);
                printf(" **   Setting shear modulus to %e\n", param->shearModulus);
                printf(" **   Setting burgers vector magnitude to %e\n", param->burgMag);
            }

            // Now check for a file containing the material-specific melting
            // temperature over a range of pressures, read the file and
            // calculate the melting temperature for the given pressure.

            if (param->meltTempFile[0] != 0) {
                GetMeltTemp(param);
            }

            // Check for a file containing material specific constants
            // pertaining to the non-linear mobility calculations over a
            // range of pressures and temperatures, read the file and
            // calculate the appropriate constants (peierls, deltaH0, alpha)
            // for the current temp/pressure

            if (param->MobFileNL[0] != 0) {
                GetNonLinearMobConstants(param);
            }

#ifndef PARALLEL
            // If the code has not been compiled for parallel execution,
            // explicitly disable the flag for doing parallel multi-node
            // split force calcs.

            param->useParallelSplitMultiNode = 0;
#endif

        }  // if (home->myDomain == 0)

        // All domains need the current Param_t structure that's
        // been populated by domain zero, so broadcast it out.  There
        // may be subsequent changes to the Param_t structure, but if
        // so, those updates will be distributed as needed.

#ifdef PARALLEL
        MPI_Bcast((char *)param, sizeof(Param_t), MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

#ifdef _OPENMP
        // Only domain 0 originally knew how many threads to use
        // per MPI task, now that the info has been distributed to
        // all domains, each domain needs to explicitly set the
        // max thread count.

        omp_set_num_threads(param->maxNumThreads);
#endif

        // All processors have the domain geometry but need to calculate
        // the maximum number of decomposition levels (for the Recursive
        // Bisection decomposition only) in each dimension *before* an
        // old decomposition is read in or a new decomposition created.

        if (param->decompType == 2)
        {
            for (home->xMaxLevel=0; param->nXdoms >> home->xMaxLevel > 1; home->xMaxLevel++);
            for (home->yMaxLevel=0; param->nYdoms >> home->yMaxLevel > 1; home->yMaxLevel++);
            for (home->zMaxLevel=0; param->nZdoms >> home->zMaxLevel > 1; home->zMaxLevel++);

            // If the number of domains in a dimension is not a power of
            // two, we'll need 1 extra level in the decomposition.  make
            // sure we don't exceed the configured limit, though.

            if ((1 << home->xMaxLevel) < param->nXdoms) home->xMaxLevel++;
            if ((1 << home->yMaxLevel) < param->nYdoms) home->yMaxLevel++;
            if ((1 << home->zMaxLevel) < param->nZdoms) home->zMaxLevel++;

            if (home->myDomain == 0) {
                if ((home->xMaxLevel >= MAX_DECOMP_LVLS) ||
                    (home->yMaxLevel >= MAX_DECOMP_LVLS) ||
                    (home->zMaxLevel >= MAX_DECOMP_LVLS)) {
                    Fatal("Decomp level exceeds max allowed value of %d.\n"
                          "    Increase MAX_DECOMP_LVLS value in Decomp.h\n"
                          "    and recompile.", MAX_DECOMP_LVLS);
                }
            }
        }

        // Each process needs to get some information about which IO group
        // it is in before the nodal data files are read in.

        GetParallelIOGroup(home);

#ifdef USE_SCR
        // See if there is a valid SCR checkpoint/restart data file
        // from which to restart.  If SCR is enabled and a valid
        // checkpoint file is found, that checkpoint overrides any
        // restart data file name provided by the user.

        char scrFileName[256]; sprintf(scrFileName, "data.%d", home->myDomain);
        if (SCR_Route_file(scrFileName, testFile) == SCR_SUCCESS)
        {
            dataFile = testFile;
            param->doBinRead = param->writeBinRestart;
        }
#endif

        // Read the nodal data (and associated parameters) from the
        // data file (which may consist of multiple file segments).

        if (param->doBinRead) { ReadBinDataFile (home, inData, dataFile); }
        else                  { ReadNodeDataFile(home, inData, dataFile); }

        // Some of the parameters used in creating the nodal data file
        // used for this run may not match the values to be used for
        // this run.  We've completed processing the data file at this
        // point, so update the data file parameters to match the values
        // desired for this particular run.

        param->dataDecompGeometry[X] = param->nXdoms;
        param->dataDecompGeometry[Y] = param->nYdoms;
        param->dataDecompGeometry[Z] = param->nZdoms;

        param->dataDecompType  = param->decompType;
        param->dataFileVersion = NODEDATA_FILE_VERSION;
        param->numFileSegments = param->numIOGroups;

        // Calculate the length of the problem space in each
        // of the dimensions

        param->Lx = param->maxSideX - param->minSideX;
        param->Ly = param->maxSideY - param->minSideY;
        param->Lz = param->maxSideZ - param->minSideZ;

        // Check which mobility module has been selected and set
        // some basic mobility parameters and select the material
        // type appropriate to the mobility.  This must be done
        // before AnisotropicInit() because that function requires
        // the proper material type for some initializations.

        SetMaterialTypeFromMobility(home);

        // If necessary, create the rotation matrices needed for rotating
        // from geometries defined in the users laboratory frame to the
        // standard crystalographic frame and back.

        // NOTE: This MUST be done before the call to VerifyBurgersVectors(),
        //       AnisotropicInit(), and any other functions that required
        //       rotation matrix used to translate between the laboratory
        //       and crystallographic frames.

        if (param->useLabFrame)
        {
            NormalizeVec(param->rotMatrix[0]);
            NormalizeVec(param->rotMatrix[1]);
            NormalizeVec(param->rotMatrix[2]);

            Matrix33_Inverse(param->rotMatrixInverse, param->rotMatrix);
        }

#if defined ANISOTROPIC
        // Now that all processors have the basic parameter info,
        // we need to initialize all variables, tables, etc. that
        // are needed for the anisotropic elasticity code.  This
        // must be called *before* SetRemainingDefaults()!

        AnisotropicInit(home);
#else
        // For anisotropic simulations, the elastic constant matrix was
        // created above in AnisotropicInit().  For isotropic simulations
        // we need to create it now.

        InitIsotropicElasticConstantMatrix(home);
#endif

#ifdef BBFMM
        bbFMInit(home);
#endif

#ifdef UNIFORMFMM
        uniformFMInit(home);
#endif

        // Calculate and save the inverse of the elastic constant matrix
        // for use later on in the code...

        Matrix66_Inverse( param->elasticConstantMatrixInv,
                          param->elasticConstantMatrix     );

        // Now that the param structure is fully populated, do any
        // remaining sanity checks or initializations that could not
        // or have not been done yet.

        SetRemainingDefaults(home);

        // Initialize the SegSeg Forces structures for
        // both GPU and non-GPU vectorized force calculations

        SSF_Initialize(home);

        // Some of the control file parameters are only used when
        // specific other parameters/toggles have been enabled for
        // the simulation.  Here we try to identify parameters that
        // are not used in the current simulation and flag them
        // so they will not be written into the restart files.  Helps
        // remove some of the clutter.

        // This only needs to be done on domain 0, but requires
        // SetRemainingDefaults() to have be called first to complete
        // initializations.

        if (home->myDomain == 0) {
            DisableUnneededParams(home);
        }

        // For some material types we need to dynamically calculate
        // the available burgers vectors (and associated glide planes).

        GetBurgList(home);

        // Some of the mobility modules require glides planes to be
        // defined for all segments.  Check for those here.

        CheckForGlidePlanes(home);

        // Do some validity checking on the constraints applied to
        // each node.

        VerifyConstraints(home);

        // Free up some of the temporary buffers that *may* have been allocated

        FreeInitArrays(home, inData);
        free(inData);

        // We attempted to preserve the node tags from the previous
        // run, so the nodeKeys array may be sparsely populated.  If
        // that is the case, we have to add all unused tags lower than
        // <newNodeKeyPtr> to the recycled node array or those tags
        // will never be used again.

        InitRecycleNodeHeap(home);

#ifdef ESHELBY
        // Read in the list of Eshelby inclusions (if any)

        ReadEshelbyInclusions(home);
#endif

        // Find out which cells intersect this domain (the native cells), and
        // their neighbors cells. Find out which domains intersect each of these
        // native and ghost cells.

        InitCellNatives  (home);
        InitCellNeighbors(home);
        InitCellDomains  (home);

        // For each neighbor domain, build a list of native cells to send to that
        // domain.

        InitRemoteDomains(home);

        // Each domain still needs to communicate with its neighboring domains
        // to map old arm tags to new ones for any nodes that were retagged
        // during initialization.

        DistributeTagMaps(home);

        // If the Fast Multipole code is not enabled, allocate an array
        // to store the cell charge tensor for each cell (but don't do
        // it if remote forces are being disabled by the FULL_N2_FORCES flag)

#ifndef FULL_N2_FORCES
        if (!param->fmEnabled)
        {
            home->cellCharge = (real8 *) malloc(param->nXcells *
                                                param->nYcells *
                                                param->nZcells *
                                                9 * sizeof(real8));
        }
#endif

        // Initialize operation list used for collisions and other topological changes.

        InitOpList(home);

#ifndef FULL_N2_FORCES
        // If the Fast Multipole code is enabled, initialize the image
        // correction table.  Otherwise (if necessary) read in PBC image
        // correction table and PBC stress tables BEFORE creation of and
        // cd to the output file directory.

        // NOTE: This call to CorrectionTableInit() MUST be done after
        // the first call to FMInit() (which is invoked from within the
        // InitCellNatives() function called above).

        if (param->fmEnabled)
        {
            // Only need the FMM correction table if PBC is enabled

            if ((param->xBoundType == Periodic) ||
                (param->yBoundType == Periodic) ||
                (param->zBoundType == Periodic)) {
                CorrectionTableInit(home);
            }
        }
        else
        {
#ifdef SPECTRAL
            if (param->elasticinteraction && param->FFTenabled == 0) {
#else
            if (param->elasticinteraction) {
#endif
                ReadRijm(home);
                ReadRijmPBC(home);
            }
        }
#endif

#ifndef NO_XWINDOW
        if (home->myDomain == 0) {
            ReadWindowSpec(param->winDefaultsFile);
        }
#endif

        // Create the output directory (and sub-directories) and reset the
        // current working directory to the top level output directory.

        if (home->param->dirname[0] != 0) {
            OpenDir(home);
        }

        // Do the initial sort of native nodes into their proper subcells
        // and send the initial ghost node data.  Previously this was at
        // the top of the ParadisStep loop, but modifications to the way
        // output (including X-windows plot data) is generated requires that
        // the cpus have accurate ghost node data before calling GenerateOutput.
        // Hence, we do the first ghost node communications here, and move
        // the subsequent ghost node comm from the beginning of the ParadisStep
        // loop to the end.

        SortNativeNodes(home);
        CommSendGhosts(home);

        // Have each processor look at all segments attached to its local
        // nodes and verify that the burgers vector is the same at both
        // ends of the node.  Just a sanity check to prevent people from
        // doing silly things like creating a nodal configuration by hand
        // and putting in inconsistent burgers vectors.

        VerifyBurgersVectors(home);

        // This call to GenerateOutput will initialize the X-Window display.

#ifndef NO_XWINDOW
        GenerateOutput(home, STAGE_INIT);
#endif

        // The code for calculating the forces from remote segments
        // needs a couple arrays for doing Gauss-Legendre integration.
        // Once created, these arrays will remain static, so allocate
        // and initialize them now.

#ifdef ESHELBY
        if ((param->fmEnabled) || (param->eshelbyfmEnabled))
#else
        if  (param->fmEnabled)
#endif
        {

            home->glPositions = (real8 *)malloc(param->fmNumPoints * sizeof(real8));
            home->glWeights   = (real8 *)malloc(param->fmNumPoints * sizeof(real8));

            GaussQuadCoeffHalf(param->fmNumPoints, home->glPositions, home->glWeights);
        }

        // Create a lock used for synchronizing access to the segment/particle
        // intersection list.

        INIT_LOCK(&home->segpartIntersectListLock);

        // We need to initialize the ARKODE integrator if it is being used
        // as the main time step integrator...

#ifdef USE_ARKODE
        if ( StrEquiv(param->timestepIntegrator  , "ARKODE"    )   ) { CreateARKode(home); }
#endif

#ifdef SPECTRAL
        InitSpectral(home);
#endif

        CheckMemUsage(home, "Initialize");

        return;
}
