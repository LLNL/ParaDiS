/****************************************************************************
 *
 *      Function:     WriteProp
 *      Description:  Depending on the input parameter, write a particular
 *                    time-dependent property out to that property's
 *                    diagnostic file, along with a timestamp.
 *
 *      NOTE: The contents of the file written below containing the
 *      dislocation density have changed over time.  The contents of
 *      the various versions are defined below.
 *
 *      Version 0 contents:
 *        - strain
 *        - dislocation density
 *        - deleted segment length
 *        - average velocity
 *        - std deviation of dislocation velocities
 *
 *      Version 1 contents: After the velocity std deviation, the
 *      following items were added:
 *        - version number
 *        - a set of values indicating the dislocation density for
 *          segments of specific groupings of burgers vectors.  These
 *          groupings differ for BCC and FCC mobility.
 *
 *          BCC groupings:
 *
 *          group #     burgers vector types
 *            0         [ 1 1 1] [-1-1-1]
 *            1         [-1 1 1] [ 1-1-1]
 *            2         [ 1-1 1] [-1 1-1]
 *            3         [ 1 1-1] [-1-1 1]
 *            4         [ 1 0 0] [-1 0 0]
 *                      [ 0 1 0] [ 0-1 0]
 *                      [ 0 0 1] [ 0 0-1]
 *
 *          FCC groupings:
 *
 *          group #     burgers vector types
 *            0         [ 1 1 0] [-1-1 0]
 *            1         [-1 1 0] [ 1-1 0]
 *            2         [ 1 0 1] [-1 0-1]
 *            3         [-1 0 1] [ 1 0-1]
 *            4         [ 0 1 1] [ 0-1-1]
 *            5         [ 0-1 1] [ 0 1-1]
 *            6         all others
 *
 *      Version 2 contents:
 *        - plastic strain added as first column of data
 *
 ****************************************************************************/
#include <stdio.h>
#include <time.h>

#include "mpi_portability.h"

#include "Home.h"
#include "WriteProp.h"
#include "Util.h"
#include "Mobility.h"

#define DENSITY_FILE_VERSION 2

void WriteProp(Home_t *home, int property)
{
        int     i;
        real8   al, am, an, amag, sigijk, pstnijk, dpstnijk;
        real8   totDensityChange[14];
        char    fileName[256];
        FILE    *fp;
        Param_t *param;


        param = home->param;

        al = param->edotdir[0];
        am = param->edotdir[1];
        an = param->edotdir[2];

        amag = sqrt(al*al+am*am+an*an);
        amag = ( (fabs(amag)>0.0) ? (1.0/amag) : 1.0 );    // (avoid potential divide by null edotdir)

        al *= amag;
        am *= amag;
        an *= amag;

        switch(property) {
        case DENSITY:
        {

#ifdef PARALLEL
            int numItems = param->numBurgGroups+1;
            int tmpOffset = 1;
#ifdef _ARLFEM
                numItems += 1;
                tmpOffset += 1;
#endif

            real8 *localDensityVals  = (real8 *)malloc((numItems)*sizeof(real8));
            real8 *globalDensityVals = (real8 *)malloc((numItems)*sizeof(real8));

            for (i = 0; i < numItems; i++) {
                globalDensityVals[i] = 0.0;
            }

            localDensityVals[0] = param->delSegLength;
#ifdef _ARLFEM
            localDensityVals[1] = param->fem_delSegLength;
#endif

            for (i = 0; i < param->numBurgGroups; i++) {
                localDensityVals[i+tmpOffset] = param->partialDisloDensity[i];
            }

            MPI_Allreduce(localDensityVals, globalDensityVals, numItems,
                          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            param->delSegLength = globalDensityVals[0];
#ifdef _ARLFEM
            param->fem_delSegLength = globalDensityVals[1];
            param->fem_delSegLength *= param->burgVolFactor;
#endif

            for (i = 0; i < param->numBurgGroups; i++) {
                param->partialDisloDensity[i] = globalDensityVals[i+tmpOffset];
            }

            free(localDensityVals);
            free(globalDensityVals);
#endif  //  PARALLEL

/*
 *          convert to density and write results to the file
 */
            param->delSegLength *= param->burgVolFactor;

            if (home->myDomain == 0) {
                snprintf(fileName, sizeof(fileName), "%s/density",
                         DIR_PROPERTIES);

                WriteDensityFileDesc(home, fileName);
                WriteDensityFileGnu (home, fileName);

                fp = fopen(fileName, "a");

                pstnijk =  param->totpStn[0]*al*al     +
                           param->totpStn[1]*am*am     +
                           param->totpStn[2]*an*an     +
                           2.0*param->totpStn[3]*am*an +
                           2.0*param->totpStn[4]*an*al +
                           2.0*param->totpStn[5]*al*am;

/*
 *              First print the standard stuff that's common
 *              regardless of which mobility (BCC, FCC, etc)
 *              is in use.
 */
                fprintf(fp, "%e %e %e %e %e %e ",
                        pstnijk,param->eRate*param->timeNow, param->disloDensity,
                        param->delSegLength, param->vAverage,
                        param->vStDev);

/*
 *              ****************************************************
 *              ***                IMPORTANT!                    ***
 *              ***   If you change the contents of the density  ***
 *              ***   file, increment the definition for         ***
 *              ***   DENSITY_FILE_VERSION above and document    ***
 *              ***   the new contents.  DO NOT REMOVE the       ***
 *              ***   description of the previous version!       ***
 *              ****************************************************
 *
 *              Now append to the line the file version number
 *              and the dislocation density for each of the
 *              groupings of burgers vectors appropriate to the
 *              mobility type.  Should be kept in sync with
 *              DeltaPlasticStrain().
 */
                fprintf(fp, " %d", DENSITY_FILE_VERSION);

                for (i = 0; i < param->numBurgGroups; i++) {
                    fprintf(fp, " %e", param->partialDisloDensity[i]);
                }

                fprintf(fp, "\n");
                fclose(fp);

#ifdef _ARLFEM
                snprintf(fileName, sizeof(fileName), "%s/fem_density",
                         DIR_PROPERTIES);
                fp = fopen(fileName, "a");
                fprintf(fp, "%d %e %e\n", home->cycle,
                        param->eRate*param->timeNow, param->fem_delSegLength);
                fclose(fp);
#endif
            }

/*
 *          Reinitialize accumulated length of deleted segments
 *          after writing
 */
            param->delSegLength = 0.0;

#ifdef _ARLFEM
            param->fem_delSegLength = 0.0;
#endif
            break;
        }
        case EPS:
        {
            if (home->myDomain == 0)
            {
                    sigijk =       param->appliedStress[0]*al*al +
                                   param->appliedStress[1]*am*am +
                                   param->appliedStress[2]*an*an +
                               2.0*param->appliedStress[3]*am*an +
                               2.0*param->appliedStress[4]*an*al +
                               2.0*param->appliedStress[5]*al*am  ;

                    pstnijk =      param->totpStn[0]*al*al +
                                   param->totpStn[1]*am*am +
                                   param->totpStn[2]*an*an +
                               2.0*param->totpStn[3]*am*an +
                               2.0*param->totpStn[4]*an*al +
                               2.0*param->totpStn[5]*al*am  ;

                snprintf(fileName, sizeof(fileName), "%s/time_Plastic_strain", DIR_PROPERTIES);

                WriteTimePlasticStrainFileDesc(home, fileName);
                WriteTimePlasticStrainFileGnu (home, fileName);

                fp = fopen(fileName, "a");
                fprintf(fp,"%e %e\n", param->timeNow, pstnijk);
                fclose(fp);

                if ((param->loadType == 1) || (param->loadType == 4)) {

                    snprintf(fileName, sizeof(fileName),
                             "%s/stress_Plastic_strain", DIR_PROPERTIES);

                    WriteStressPlasticStrainFileDesc(home, fileName);
                    WriteStressPlasticStrainFileGnu (home, fileName);

                    fp = fopen(fileName, "a");
                    fprintf(fp, "%e %e\n", pstnijk, sigijk);
                    fclose(fp);

                    snprintf(fileName, sizeof(fileName),
                             "%s/stress_Total_strain", DIR_PROPERTIES);

                    WriteStressTotalStrainFileDesc(home, fileName);

                    fp = fopen(fileName, "a");
                    switch (param->loadType) {
                    case 1:
                        fprintf(fp, "%e %e %e %d\n",
                                param->eRate*param->timeNow, sigijk,
                                param->timeNow, home->cycle);
                        break;
                    case 4:
                        fprintf(fp, "%e %e %e %d\n",
                                param->netCyclicStrain, sigijk,
                                param->timeNow,param->numLoadCycle);
                        break;
                    }
                    fclose(fp);

                    snprintf(fileName, sizeof(fileName),
                             "%s/total_stress", DIR_PROPERTIES);

                    WriteTotalStressFileDesc(home, fileName);

                    fp = fopen(fileName, "a");
                    fprintf(fp, "%d %e %e %e %e %e %e %e\n",
                            home->cycle, param->timeNow,
                            param->appliedStress[0], param->appliedStress[1],
                            param->appliedStress[2], param->appliedStress[3],
                            param->appliedStress[4], param->appliedStress[5]);
                    fclose(fp);

#if 1
                    snprintf(fileName, sizeof(fileName),
                             "%s/edir", DIR_PROPERTIES);
                    fp = fopen(fileName, "a");
                    fprintf(fp, "%d %e %e %e %e\n",
                            home->cycle, param->eRate*param->timeNow,
                            al, am, an);
                    fclose(fp);

                    if (param->savetimers) {
                    #ifdef PARALLEL
                        real8 now = MPI_Wtime();
                    #else
                        real8 now = (real8)clock() / CLOCKS_PER_SEC;
                    #endif
                        real8 tottime = home->timers[0].accum + now - home->timers[0].startTime;

                        snprintf(fileName, sizeof(fileName),
                                 "%s/time", DIR_PROPERTIES);
                        fp = fopen(fileName, "a");
                        fprintf(fp, "%d %e %e %e %e\n",
                                home->cycle, param->eRate*param->timeNow,
                                param->timeNow, param->realdt, tottime);
                        fclose(fp);
                    }
#endif

                }
            }
            break;
        }
        case ALL_EPS:
        {
            if (home->myDomain == 0) {
                snprintf(fileName, sizeof(fileName), "%s/alleps",
                         DIR_PROPERTIES);

                WriteAllepsFileDesc(home, fileName);

                fp=fopen(fileName, "a");

                fprintf(fp,"%d %e %e %e %e %e %e %e %e\n",
                        home->cycle, param->timeNow,
                        param->totpStn[0],param->totpStn[1],param->totpStn[2],
                        param->totpStn[3],param->totpStn[4],param->totpStn[5],
                        param->disloDensity );
                fclose(fp);
            }
            break;
        }
        case EPSDOT:
        {
            if (home->myDomain == 0) {
                if (param->realdt > 1.e-20){
                    dpstnijk =  param->delpStrain[0]*al*al     +
                                param->delpStrain[1]*am*am     +
                                param->delpStrain[2]*an*an     +
                                2.0*param->delpStrain[3]*am*an +
                                2.0*param->delpStrain[4]*an*al +
                                2.0*param->delpStrain[5]*al*am;

                    snprintf(fileName, sizeof(fileName), "%s/epsdot",
                             DIR_PROPERTIES);

                    WriteEpsdotFileDesc(home, fileName);

                    fp = fopen(fileName, "a");
                    fprintf(fp, "%e %e\n", param->timeNow,
                            fabs(dpstnijk/param->realdt));
                    fclose(fp);
                }
            }
            break;
        }
        case DENSITY_DELTA:
        {
/*
 *          Currently only supported for BCC mobilities
 */
            if (param->materialType == MAT_TYPE_BCC) {

/*
 *              We track 2 values (density gain and density loss) for
 *              each of 7 burgers vectors
 */
                int numBurgVectors = 7;
                int numItems       = numBurgVectors * 2;
#ifdef PARALLEL

                for (i = 0; i < numItems; i++) {
                    totDensityChange[i] = 0.0;
                }

                MPI_Reduce(param->densityChange, totDensityChange, numItems,
                           MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
                for (i = 0; i < numItems; i++) {
                    totDensityChange[i] = param->densityChange[i];
                }
#endif
                if (home->myDomain == 0) {
                    real8 totGain = 0.0;
                    real8 totLoss = 0.0;

                    for (i = 0; i < numBurgVectors; i++) {
                        totGain += totDensityChange[i];
                        totLoss += totDensityChange[numBurgVectors+i];
                    }

                    snprintf(fileName, sizeof(fileName), "%s/density_delta",
                             DIR_PROPERTIES);


                    WriteDensityDeltaFileDesc(home, fileName);

                    fp = fopen(fileName, "a");

                    fprintf(fp, "%e %e %e ", param->eRate*param->timeNow,
                            totGain, totLoss);

                    for (i = 0; i < numItems; i++) {
                        fprintf(fp, "%e ", totDensityChange[i]);
                    }
                    fprintf(fp, "\n");

                    fclose(fp);
                }
/*
 *              Zero out the accumulated values since we've written the
 *              data to disk.
 */
                for (i = 0; i < numItems; i++) {
                    param->densityChange[i] = 0.0;
                }
            }
            break;
        }

#ifdef CALCENERGY
        case ENERGY:
        {
            if (home->myDomain == 0) {

               real8 TotalEnergy = 0.0;

                TotalEnergy = param->TotalEnergy;

                snprintf(fileName, sizeof(fileName), "%s/energy", DIR_PROPERTIES);

                WriteTimeEnergyFileDesc(home, fileName);

                fp = fopen(fileName, "a");
                fprintf(fp,"%e %e\n", param->timeNow, TotalEnergy);
                fclose(fp);

                if ((param->loadType == 1) || (param->loadType == 4)) {

                    snprintf(fileName, sizeof(fileName),
                             "%s/energy_Plastic_strain", DIR_PROPERTIES);

                    WriteEnergyPlasticStrainFileDesc(home, fileName);

                    fp = fopen(fileName, "a");
                    fprintf(fp, "%e %e\n", pstnijk, TotalEnergy);
                    fclose(fp);

                    snprintf(fileName, sizeof(fileName),
                             "%s/energy_Total_strain", DIR_PROPERTIES);

                    WriteEnergyStrainFileDesc(home, fileName);

                    fp = fopen(fileName, "a");
                    switch (param->loadType) {
                    case 1:
                        fprintf(fp, "%e %e %e %d\n",
                                param->eRate*param->timeNow, TotalEnergy,
                                param->timeNow, home->cycle);
                        break;
                    case 4:
                        fprintf(fp, "%e %e %e %d\n",
                                param->netCyclicStrain, TotalEnergy,
                                param->timeNow,param->numLoadCycle);
                        break;
                    }
                    fclose(fp);
                }
            }
            break;
        }
#endif

        default:
        {
            Fatal("WriteProp: input parameter (property=%d) invalid", property);
            break;
        }

        }  /* switch() */

        return;
}
