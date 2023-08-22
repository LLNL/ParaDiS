/****************************************************************************
 *
 *      Function:    WriteDensFlux
 *      Description: Writes information required for Tom's 
 *                   dislocation density-based continuum model
 *                     - Density and flux for each slip system
 *                     - Amount of slip (in strain) for all slip systems
 *                
 *      Authors:      
 *               Moon Rhee 06/23/2004
 *               Updated Sylvie Aubry for HCP in 2013
 *               Updated Sylvie Aubry for BCC more general decomposition 05/18/2016
 *
 *
 *      Includes public functions:
 *          WriteDensFlux()
 *
 *      Includes private functions:
 *          WriteDensFlux_BCC2()
 *          WriteDensFlux_BCC()
 *          WriteDensFlux_FCC()
 *          WriteDensFlux_HCP()
 *          WriteDensFlux_RhomboVa()
 *
 ****************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"

// WriteDensFlux_BCC2 uses the decomposition of the density and flux into (110) and (112)
// types directions. A total of 6 possible planes.
static void WriteDensFlux_BCC2(Home_t *home, real8 pstnijk, real8 tmpstn)
{
        int        i, j, bIndex;
        char       fileName[512];
        Param_t    *param;
        FILE       *fp;
        static int numBurg   = 4;
        static int numPlanes = 6;
        static int sizeFlux  = 13;

        // Define indices into the BCC2_LtotDecomp array for various items
        static int indxTotEdge  = 0;
        static int indxEdges    = 1;
        static int indxTotScrew = 7;
        static int indxScrews   = 8;
        static int indxTotals   = 14;

//
//   Density data is contained in the BCC2_LtotDecomp array where
//   the per-burgers vector data components are organized as:
//
//   [burg][0]   total edge_density (based on total length - screw length)
//   [burg][1]   edge_density plane 1
//   [burg][2]   edge_density plane 2
//   [burg][3]   edge_density plane 3
//   [burg][4]   edge_density plane 4
//   [burg][5]   edge_density plane 5
//   [burg][6]   edge_density plane 6
//   [burg][7]   total screw_density (based on screw length)
//   [burg][8]   screw_density plane 1
//   [burg][9]   screw_density plane 2
//   [burg][10]  screw_density plane 3
//   [burg][11]  screw_density plane 4
//   [burg][12]  screw_density plane 5
//   [burg][13]  screw_density plane 6
//   [burg][14]  total_density plane 1 = sqrt(edge_density1^2 + screw_density1^2)
//   [burg][15]  total_density plane 2 = sqrt(edge_density2^2 + screw_density2^2)
//   [burg][16]  total_density plane 3 = sqrt(edge_density3^2 + screw_density3^2)
//   [burg][17]  total_density plane 4 = sqrt(edge_density1^2 + screw_density1^2)
//   [burg][18]  total_density plane 5 = sqrt(edge_density2^2 + screw_density2^2)
//   [burg][19]  total_density plane 6 = sqrt(edge_density2^2 + screw_density2^2)

//   Flux contents:
//   0-1:  flux due to climb
//   1-6:  flux from edge components on 6 planes
//   7-12: flux from screw components on 6 planes
//


        param = home->param;

        WriteDensFluxDesc_BCC2(home);


        for (bIndex = 0; bIndex < numBurg; bIndex++) {

/*
 *          Create the Ltot_bN file if requested
 */
            if (param->writeFluxEdgeDecomp) {
                real8 allEdge, allScrew;
                real8 totalDensityEdge[4];

                allEdge = 0.0;
                allScrew = 0.0;

/*
 *              Get some totals we'll need for later
 */
                for (i = 0; i < numBurg; i++) {

                    totalDensityEdge[i] = 0.0;

                    for (j = 0; j < numPlanes; j++) {
                        totalDensityEdge[i]  +=
                                param->BCC2_LtotDecomp[i][indxEdges+j];
                    }

                    allEdge  += totalDensityEdge[i];
                    allScrew += param->BCC2_LtotDecomp[i][indxTotScrew];
                }

                snprintf(fileName, sizeof(fileName), "%s/Ltot_b%d",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // Total screw density
                fprintf(fp, "%e ", param->BCC2_LtotDecomp[bIndex][indxTotScrew]);

                // The per-plane edge density
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ", param->BCC2_LtotDecomp[bIndex][indxEdges+i]);
                }

                // Some totals...
                fprintf(fp,"%e %e %e ", totalDensityEdge[bIndex],
                        allEdge, allScrew);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxEdgeDecomp)

/*
 *          Create the Ltot_bN_full_decomp file if requested
 */
            if (param->writeFluxFullDecomp) {
                snprintf(fileName, sizeof(fileName), "%s/Ltot_b%d_full_decomp",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The per-plane screw density values
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ",
                            param->BCC2_LtotDecomp[bIndex][indxScrews+i]);
                }

                // The per-plane edge density values
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ",
                            param->BCC2_LtotDecomp[bIndex][indxEdges+i]);
                }

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxFullDecomp)

/*
 *          Create the Ltot_bN_full_decomp_totals file if requested
 */
            if (param->writeFluxFullDecompTotals) {

                snprintf(fileName, sizeof(fileName),
                         "%s/Ltot_b%d_full_decomp_totals", DIR_FLUXDATA,
                         bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The per-plane density totals followed by the sum of
                // all those values.

                real8 sumTotal = 0.0;
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ", param->BCC2_LtotDecomp[bIndex][indxTotals+i]);
                    sumTotal += param->BCC2_LtotDecomp[bIndex][indxTotals+i];
                }

                fprintf(fp, "%e ", sumTotal);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxFullDecompTotals)


/*
 *          Create the Ltot_bN_simple_totals file if requested
 */
            if (param->writeFluxSimpleTotals) {

                snprintf(fileName, sizeof(fileName),
                         "%s/Ltot_b%d_simple_totals",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The total screw density
                fprintf(fp, "%e ", param->BCC2_LtotDecomp[bIndex][indxTotScrew]);

                // The total edge density
                fprintf(fp, "%e ", param->BCC2_LtotDecomp[bIndex][indxTotEdge]);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxSimpleTotals)

/*
 *          Create the fluxtot_bN file
 */
            snprintf(fileName, sizeof(fileName), "%s/fluxtot_b%d",
                     DIR_FLUXDATA, bIndex+1);

            fp = fopen(fileName,"a");

            fprintf(fp, "%e %e ", pstnijk, tmpstn);

            for (i = 0; i < sizeFlux; i++) {
                fprintf(fp, "%e ", param->BCC2_fluxtot[bIndex][i]);
            }
            fprintf(fp,"%e \n", param->realdt);

            fclose(fp);
        }

        return;
}

// WriteDensFlux_BCC uses the decomposition of the density and flux into (110)
// types directions. A total of 3 possible planes.
static void WriteDensFlux_BCC(Home_t *home, real8 pstnijk, real8 tmpstn)
{
        int        i, j, bIndex;
        char       fileName[512];
        Param_t    *param;
        FILE       *fp;
        static int numBurg = 4;
        static int numPlanes = 3;
        static int sizeFlux = 7;
        // Define indices into the FCC_LtotDecomp array for various items
        static int indxTotEdge  = 0;
        static int indxEdges    = 1;
        static int indxTotScrew = 4;
        static int indxScrews   = 5;
        static int indxTotals   = 8;

//
//      Density data is contained in the BCC_LtotDecomp array where
//      the per-burgers vector data components are organized as:
//
//        [burg][0]  total edge_density (based on total length - screw length)
//        [burg][1]  edge_density plane 1
//        [burg][2]  edge_density plane 2
//        [burg][3]  edge_density plane 3
//        [burg][4]  total screw_density (based on screw length)
//        [burg][5]  screw_density plane 1
//        [burg][6]  screw_density plane 2
//        [burg][7]  screw_density plane 3
//        [burg][8]  total_density plane 1 = sqrt(edge_density1^2+screw_density1^2)
//        [burg][9]  total_density plane 2 = sqrt(edge_density2^2+screw_density2^2)
//        [burg][10] total_density plane 3 = sqrt(edge_density3^2+screw_density3^2)
//

        param = home->param;

        WriteDensFluxDesc_BCC(home);


        for (bIndex = 0; bIndex < numBurg; bIndex++) {

/*
 *          Create the Ltot_bN file if requested
 */
            if (param->writeFluxEdgeDecomp) {
                real8 allEdge, allScrew;
                real8 totalDensityEdge[4];

                allEdge = 0.0;
                allScrew = 0.0;

/*
 *              Get some totals we'll need for later
 */
                for (i = 0; i < numBurg; i++) {

                    totalDensityEdge[i] = 0.0;

                    for (j = 0; j < numPlanes; j++) {
                        totalDensityEdge[i]  +=
                                param->BCC_LtotDecomp[i][indxEdges+j];
                    }

                    allEdge  += totalDensityEdge[i];
                    allScrew += param->BCC_LtotDecomp[i][indxTotScrew];
                }

                snprintf(fileName, sizeof(fileName), "%s/Ltot_b%d",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // Total screw density
                fprintf(fp, "%e ", param->BCC_LtotDecomp[bIndex][indxTotScrew]);

                // The per-plane edge density
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ", param->BCC_LtotDecomp[bIndex][indxEdges+i]);
                }

                // Some totals...
                fprintf(fp,"%e %e %e ", totalDensityEdge[bIndex],
                        allEdge, allScrew);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxEdgeDecomp)

/*
 *          Create the Ltot_bN_full_decomp file if requested
 */
            if (param->writeFluxFullDecomp) {
                snprintf(fileName, sizeof(fileName), "%s/Ltot_b%d_full_decomp",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The per-plane screw density values
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ",
                            param->BCC_LtotDecomp[bIndex][indxScrews+i]);
                }

                // The per-plane edge density values
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ",
                            param->BCC_LtotDecomp[bIndex][indxEdges+i]);
                }

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxFullDecomp)

/*
 *          Create the Ltot_bN_full_decomp_totals file if requested
 */
            if (param->writeFluxFullDecompTotals) {

                snprintf(fileName, sizeof(fileName),
                         "%s/Ltot_b%d_full_decomp_totals", DIR_FLUXDATA,
                         bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The per-plane density totals followed by the sum of
                // all those values.

                real8 sumTotal = 0.0;
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ", param->BCC_LtotDecomp[bIndex][indxTotals+i]);
                    sumTotal += param->BCC_LtotDecomp[bIndex][indxTotals+i];
                }

                fprintf(fp, "%e ", sumTotal);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxFullDecompTotals)


/*
 *          Create the Ltot_bN_simple_totals file if requested
 */
            if (param->writeFluxSimpleTotals) {

                snprintf(fileName, sizeof(fileName),
                         "%s/Ltot_b%d_simple_totals",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The total screw density
                fprintf(fp, "%e ", param->BCC_LtotDecomp[bIndex][indxTotScrew]);

                // The total edge density
                fprintf(fp, "%e ", param->BCC_LtotDecomp[bIndex][indxTotEdge]);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxSimpleTotals)

/*
 *          Create the fluxtot_bN file
 */
            snprintf(fileName, sizeof(fileName), "%s/fluxtot_b%d",
                     DIR_FLUXDATA, bIndex+1);

            fp = fopen(fileName,"a");

            fprintf(fp, "%e %e ", pstnijk, tmpstn);

            for (i = 0; i < sizeFlux; i++) {
                fprintf(fp, "%e ", param->BCC_fluxtot[bIndex][i]);
            }
            fprintf(fp,"%e \n", param->realdt);

            fclose(fp);
        }

        return;
}


static void WriteDensFlux_FCC(Home_t *home, real8 pstnijk, real8 tmpstn)
{
        int        i, j, bIndex;
        char       fileName[512];
        Param_t    *param;
        FILE       *fp;
        static int numBurg = 6;
        static int numPlanes = 3;
        static int sizeFlux = 7;
        // Define indices into the FCC_LtotDecomp array for various items
        static int indxTotEdge  = 0;
        static int indxEdges    = 1;
        static int indxTotScrew = 4;
        static int indxScrews   = 5;
        static int indxTotals   = 8;

//
//      Density data is contained in the FCC_LtotDecomp array where
//      the per-burgers vector data components are organized as:
//
//      [burg][0]  total edge_density (based on total length - screw length)
//      [burg][1]  edge_density plane 1
//      [burg][2]  edge_density plane 2
//      [burg][3]  edge_density plane 3
//      [burg][4]  total screw_density (based on screw length)
//      [burg][5]  screw_density plane 1
//      [burg][6]  screw_density plane 2
//      [burg][7]  screw_density plane 3
//      [burg][8]  total_dens plane 1 = sqrt(edge_density1^2+screw_density1^2)
//      [burg][9]  total_dens plane 2 = sqrt(edge_density2^2+screw_density2^2)
//      [burg][10] total_dens plane 3 = sqrt(edge_density3^2+screw_density3^2)
//

        param = home->param;

        WriteDensFluxDesc_FCC(home);

        for (bIndex = 0; bIndex < numBurg; bIndex++) {

/*
 *          Create the Ltot_bN file if requested
 */
            if (param->writeFluxEdgeDecomp) {
                real8 allEdge, allScrew;
                real8 totalDensityEdge[6];

                allEdge = 0.0;
                allScrew = 0.0;

/*
 *              Get some totals we'll need for later
 */
                for (i = 0; i < numBurg; i++) {

                    totalDensityEdge[i] = 0.0;

                    for (j = 0; j < numPlanes; j++) {
                        totalDensityEdge[i]  +=
                                param->FCC_LtotDecomp[i][indxEdges+j];
                    }

                    allEdge  += totalDensityEdge[i];
                    allScrew += param->FCC_LtotDecomp[i][indxTotScrew];
                }

                snprintf(fileName, sizeof(fileName), "%s/Ltot_b%d",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // Total screw density
                fprintf(fp, "%e ", param->FCC_LtotDecomp[bIndex][indxTotScrew]);

                // The per-plane edge density
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ",
                            param->FCC_LtotDecomp[bIndex][indxEdges+i]);
                }

                // Some totals...
                fprintf(fp,"%e %e %e ", totalDensityEdge[bIndex],
                        allEdge, allScrew);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);
            }  // end if (param->writeFluxEdgeDecomp)

/*
 *          Create the Ltot_bN_full_decomp file if requested
 */
            if (param->writeFluxFullDecomp) {
                snprintf(fileName, sizeof(fileName), "%s/Ltot_b%d_full_decomp",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The per-plane screw density values
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ",
                            param->FCC_LtotDecomp[bIndex][indxScrews+i]);
                }

                // The per-plane edge density values
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ",
                            param->FCC_LtotDecomp[bIndex][indxEdges+i]);
                }

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxFullDecomp)

/*
 *          Create the Ltot_bN_full_decomp_totals file if requested
 */
            if (param->writeFluxFullDecompTotals) {

                snprintf(fileName, sizeof(fileName),
                         "%s/Ltot_b%d_full_decomp_totals", DIR_FLUXDATA,
                         bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The per-plane density totals followed by the sum of
                // all those values.

                real8 sumTotal = 0.0;
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ", param->FCC_LtotDecomp[bIndex][indxTotals+i]);
                    sumTotal += param->FCC_LtotDecomp[bIndex][indxTotals+i];
                }

                fprintf(fp, "%e ", sumTotal);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxFullDecompTotals)

/*
 *          Create the Ltot_bN_simple_totals file if requested
 */
            if (param->writeFluxSimpleTotals) {

                snprintf(fileName, sizeof(fileName),
                         "%s/Ltot_b%d_simple_totals",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The total screw density
                fprintf(fp, "%e ", param->FCC_LtotDecomp[bIndex][indxTotScrew]);

                // The total edge density
                fprintf(fp, "%e ", param->FCC_LtotDecomp[bIndex][indxTotEdge]);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxSimpleTotals)

/*
 *          Create the fluxtot_bN file
 */
            snprintf(fileName, sizeof(fileName), "%s/fluxtot_b%d",
                     DIR_FLUXDATA, bIndex+1);

            fp = fopen(fileName,"a");

            fprintf(fp, "%e %e ", pstnijk, tmpstn);

            for (i = 0; i < sizeFlux; i++) {
                fprintf(fp, "%e ", param->FCC_fluxtot[bIndex][i]);
            }
            fprintf(fp,"%e \n", param->realdt);

            fclose(fp);
        }

        return;
}


static void WriteDensFlux_HCP(Home_t *home, real8 pstnijk, real8 tmpstn)
{
        int        i, j, bIndex;
        char       fileName[512];
        Param_t    *param;
        FILE       *fp;
        static int numBurg = 10;
        // numPlanes is set to account for the 4 glissile planes plus one
        // additional grouping for all other planes
        static int numPlanes = 5;
        static int sizeFlux = 9;
        // Define indices into the HCP_LtotDecomp array for various items
        static int indxTotEdge  = 0;
        static int indxEdges    = 1;
        static int indxTotScrew = 6;
        static int indxScrews   = 7;
        static int indxTotals   = 12;

//
//      Density data is contained in the HCP_LtotDecomp array where
//      the per-burgers vector data components are organized as:
//
//      [burg][0]  total edge_density (based on total length - screw length)
//      [burg][1]  edge_density glissile plane 1
//      [burg][2]  edge_density glissile plane 2
//      [burg][3]  edge_density glissile plane 3
//      [burg][4]  edge_density glissile plane 4
//      [burg][5]  edge_density all other planes
//      [burg][6]  total screw_density (based on screw length)
//      [burg][7]  screw_density glissile plane 1
//      [burg][8]  screw_density glissile plane 2
//      [burg][9]  screw_density glissile plane 3
//      [burg][10] screw_density glissile plane 4
//      [burg][11] screw_density all other planes
//      [burg][12] total_dens plane 1 = sqrt(edge_density1^2+screw_density1^2)
//      [burg][13] total_dens plane 2 = sqrt(edge_density2^2+screw_density2^2)
//      [burg][14] total_dens plane 3 = sqrt(edge_density3^2+screw_density3^2)
//      [burg][15] total_dens plane 4 = sqrt(edge_density3^2+screw_density3^2)
//      [burg][16] total_dens other =   sqrt(edge_other^2+screw_other^2)
//

        param = home->param;

        WriteDensFluxDesc_HCP(home);


        for (bIndex = 0; bIndex < numBurg; bIndex++) {

/*
 *          Create the Ltot_bN file if requested
 */
            if (param->writeFluxEdgeDecomp) {
                real8 allEdge, allScrew;
                real8 totalDensityEdge[10];

                allEdge = 0.0;
                allScrew = 0.0;

/*
 *              Get some totals we'll need for later
 */
                for (i = 0; i < numBurg; i++) {

                    totalDensityEdge[i]  = 0.0;

                    for (j = 0; j < numPlanes; j++) {
                        totalDensityEdge[i] +=
                                param->HCP_LtotDecomp[i][indxEdges+j];
                    }

                    allEdge  += totalDensityEdge[i];
                    allScrew += param->HCP_LtotDecomp[i][indxTotScrew];
                }

                snprintf(fileName, sizeof(fileName), "%s/Ltot_b%d",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // Total screw density
                fprintf(fp, "%e ",
                        param->HCP_LtotDecomp[bIndex][indxTotScrew]);

                // The per-plane edge density
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ",
                            param->HCP_LtotDecomp[bIndex][indxEdges+i]);
                }

                // Some totals...
                fprintf(fp,"%e %e %e ", totalDensityEdge[bIndex],
                        allEdge, allScrew);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxEdgeDecomp)

/*
 *          Create the Ltot_bN_full_decomp file if requested
 */
            if (param->writeFluxFullDecomp) {
                snprintf(fileName, sizeof(fileName), "%s/Ltot_b%d_full_decomp",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The per-plane screw density values
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ",
                            param->HCP_LtotDecomp[bIndex][indxScrews+i]);
                }

                // The per-plane edge density values
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ",
                            param->HCP_LtotDecomp[bIndex][indxEdges+i]);
                }

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxFullDecomp)

/*
 *          Create the Ltot_bN_full_decomp_totals file if requested
 */
            if (param->writeFluxFullDecompTotals) {

                snprintf(fileName, sizeof(fileName),
                         "%s/Ltot_b%d_full_decomp_totals", DIR_FLUXDATA,
                         bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The per-plane density totals followed by the sum of
                // all those values.

                real8 sumTotal = 0.0;
                for (i = 0; i < numPlanes; i++) {
                    fprintf(fp, "%e ", param->HCP_LtotDecomp[bIndex][indxTotals+i]);
                    sumTotal += param->HCP_LtotDecomp[bIndex][indxTotals+i];
                }

                fprintf(fp, "%e ", sumTotal);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxFullDecompTotals)

/*
 *          Create the Ltot_bN_simple_totals file if requested
 */
            if (param->writeFluxSimpleTotals) {

                snprintf(fileName, sizeof(fileName),
                         "%s/Ltot_b%d_simple_totals",
                         DIR_FLUXDATA, bIndex+1);

                fp = fopen(fileName, "a");

                // Current plastic strain and strain values
                fprintf(fp, "%e %e ", pstnijk, tmpstn);

                // The total screw density
                fprintf(fp, "%e ",
                        param->HCP_LtotDecomp[bIndex][indxTotScrew]);

                // The total edge density
                fprintf(fp, "%e ",
                        param->HCP_LtotDecomp[bIndex][indxTotEdge]);

                // And lastly put in the delta t value
                fprintf(fp,"%e \n", param->realdt);

                fclose(fp);

            }  // end if (param->writeFluxSimpleTotals)

/*
 *          Create the fluxtot_bN file
 */
            snprintf(fileName, sizeof(fileName), "%s/fluxtot_b%d",
                     DIR_FLUXDATA, bIndex+1);

            fp = fopen(fileName,"a");

            fprintf(fp, "%e %e ", pstnijk, tmpstn);

/*
 *          For the first 9 burgers vectors, the flux_tot array data are:
 *
 *              1:   flux due to climb
 *              2-5: flux from edge components on 4 planes
 *              6-9: flux from screw components on 4 planes
 *
 *          but the 10th burgers vector has only 3 planes and the flux_tot
 *          array data are:
 *
 *              1:   flux due to climb
 *              2-4: flux from edge components on 3 planes
 *              5-7: flux from screw components on 3 planes
 *              8-9: zero
 *
 *          In the output files, we want the format to be the same for
 *          ALL the files, so we have to explicitly write out the
 *          data for the 10th burgers vector so it looks like it also
 *          has 4 planes, but the data for the fourth plane will always
 *          be zeroes.
 */
            if (bIndex < numBurg-1) {
                for (i = 0; i < sizeFlux; i++) {
                    fprintf(fp, "%e ", param->HCP_fluxtot[bIndex][i]);
                }
            } else {
                fprintf(fp, "%e ", param->HCP_fluxtot[bIndex][0]);
                fprintf(fp, "%e ", param->HCP_fluxtot[bIndex][1]);
                fprintf(fp, "%e ", param->HCP_fluxtot[bIndex][2]);
                fprintf(fp, "%e ", param->HCP_fluxtot[bIndex][3]);
                fprintf(fp, "%e ", 0.0);  /* edge component for 4th plane */
                fprintf(fp, "%e ", param->HCP_fluxtot[bIndex][4]);
                fprintf(fp, "%e ", param->HCP_fluxtot[bIndex][5]);
                fprintf(fp, "%e ", param->HCP_fluxtot[bIndex][6]);
                fprintf(fp, "%e ", 0.0);  /* screw component for 4th plane */
            }

            fprintf(fp,"%e \n", param->realdt);

            fclose(fp);
        }

        return;
}


static void WriteDensFlux_RhomboVa(Home_t *home, real8 pstnijk, real8 tmpstn)
{
        int        bIndex;
        char       fileName[512];
        Param_t    *param;
        FILE       *fp;
        static int numBurg = 4;


        param = home->param;

        WriteDensFluxDesc_RhomboVa(home);

        for (bIndex = 0; bIndex < numBurg; bIndex++) {

            snprintf(fileName, sizeof(fileName), "%s/Ltot_b%d",
                     DIR_FLUXDATA, bIndex+1);

            fp = fopen(fileName,"a");

            fprintf(fp, "%e %e  %e %e  %e %e  %e %e  %e \n", pstnijk, tmpstn,
                       param->rhombo_Ltot[bIndex][0][0],
                       param->rhombo_Ltot[bIndex][0][1],
                       param->rhombo_Ltot[bIndex][1][0],
                       param->rhombo_Ltot[bIndex][1][1],
                       param->rhombo_Ltot[bIndex][2][0],
                       param->rhombo_Ltot[bIndex][2][1],
                       param->realdt);
            fclose(fp);

            snprintf(fileName, sizeof(fileName), "%s/fluxtot_b%d",
                     DIR_FLUXDATA, bIndex+1);

            fp = fopen(fileName,"a");

            fprintf(fp, "%e %e  %e %e %e  %e %e %e  %e %e %e %e \n",
                    pstnijk, tmpstn,
                    /* climb */
                    param->rhombo_fluxtot[bIndex][0][0],
                    param->rhombo_fluxtot[bIndex][0][1],
                    param->rhombo_fluxtot[bIndex][0][2],
                    /* screw */
                    param->rhombo_fluxtot[bIndex][1][0],
                    param->rhombo_fluxtot[bIndex][1][1],
                    param->rhombo_fluxtot[bIndex][1][2],
                    /* 75 degree */
                    param->rhombo_fluxtot[bIndex][2][0],
                    param->rhombo_fluxtot[bIndex][2][1],
                    param->rhombo_fluxtot[bIndex][2][2],
                    param->realdt);
/*
 *                  rhombo_fluxtot[A][B][C]:
 *                      A: burgers vector (4)
 *                      B: plane (3)
 *                      C: index 0 == climb
 *                         index 1 == screw
 *                         index 2 == 75 degree
 */
            fclose(fp);
        }

        return;
}


void WriteDensFlux(char *fluxname, Home_t *home)
{
        real8   al, am, an, amag, pstnijk;
        real8   tmpstn;
        Param_t *param;

        param = home->param;

        if (home->myDomain != 0) {
            return;
        }
 
        param = home->param;

        al = param->edotdir[0];
        am = param->edotdir[1];
        an = param->edotdir[2];

        amag = sqrt(al*al + am*am + an*an);
        amag = ( (fabs(amag)>0.0) ? (1.0/amag) : 1.0 );    // (avoid potential divide by null edotdir)

        al *= amag;
        am *= amag;
        an *= amag;

        tmpstn = param->eRate * param->timeNow;

        pstnijk =  param->totpStn[0]*al*al     +
                   param->totpStn[1]*am*am     +
                   param->totpStn[2]*an*an     +
                   2.0*param->totpStn[3]*am*an +
                   2.0*param->totpStn[4]*an*al +
                   2.0*param->totpStn[5]*al*am;

        switch(param->materialType) {

            case MAT_TYPE_BCC:

               switch(param->bcc_DensityFluxDecomp) {
                   case 1: 
                      WriteDensFlux_BCC(home, pstnijk, tmpstn);
                      break;

                    case 2:
                       WriteDensFlux_BCC2(home, pstnijk, tmpstn);                 
                       break;
                    
                    default:
                       WriteDensFlux_BCC(home, pstnijk, tmpstn);
                       break;    
               }
               break;

            case MAT_TYPE_FCC:
                WriteDensFlux_FCC(home, pstnijk, tmpstn);
                break;

            case MAT_TYPE_HCP:
                WriteDensFlux_HCP(home, pstnijk, tmpstn);
                break;

            case MAT_TYPE_RHOMBOHEDRAL_VA:
                WriteDensFlux_RhomboVa(home, pstnijk, tmpstn);
                break;
        }

        return;
}
