
/****************************************************************************
 *
 *      Function:    WriteDensFluxDescription
 *      Description: The flux output data consists of two sets of files.
 *                   This module contains functions that will write
 *                   corresponding description files (appropriate to the
 *                   material being simulated) for each of the sets of
 *                   flux data.  The description files briefly describe the
 *                   contents and format of the flux output files.
 *                
 *      Includes public functions:
 *          WriteDensFluxDesc_BCC()
 *          WriteDensFluxDesc_BCC2()
 *          WriteDensFluxDesc_FCC()
 *          WriteDensFluxDesc_HCP()
 *          WriteDensFluxDesc_RhomboVa()
 *
 ****************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"


/****************************************************************************
 *
 *      Function:     WriteDensFluxDesc_BCC2
 *      Description:  Create a file containing a description of the contents
 *                    of the two sets of flux decomposition output files
 *                    appropriate to BCC simulations.
 *
 ****************************************************************************/
void WriteDensFluxDesc_BCC2(Home_t *home)
{
    char        fileName[MAX_STRING_LEN];
    FILE       *fpDesc;
    Param_t    *param = home->param;
    static int  descFileWritten = 0;

/*
 *  Only one task (domain zero) writes this file, and only writes
 *  it one timer per simulation.  If the associated description
 *  file already exists, the contents are overwritten
 */
    if ((home->myDomain != 0) || (descFileWritten == 1)) {
        return;
    }

    descFileWritten = 1;

/*
 *  Open and write the description file associated with the Ltot_* 
 *  flux output files if requested.
 */
    if (param->writeFluxEdgeDecomp) {

        sprintf(fileName, "%s/Ltot_bx.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The BCC2 Ltot_bx files each correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1    Burgers vector 1/2 [ 1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b2    Burgers vector 1/2 [ 1  1 -1]\n");
        fprintf(fpDesc, "    Ltot_b3    Burgers vector 1/2 [ 1 -1  1]\n");
        fprintf(fpDesc, "    Ltot_b4    Burgers vector 1/2 [-1  1  1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_bx files are identical in format and contain\n");
        fprintf(fpDesc, "the following data for the associated Burgers vectors\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Total screw density\n");

        fprintf(fpDesc, "4         Edge density for glide plane p1\n");
        fprintf(fpDesc, "5         Edge density for glide plane p2\n");
        fprintf(fpDesc, "6         Edge density for glide plane p3\n");
        fprintf(fpDesc, "7         Edge density for glide plane p4\n");
        fprintf(fpDesc, "8         Edge density for glide plane p5\n");
        fprintf(fpDesc, "9         Edge density for glide plane p6\n");

        fprintf(fpDesc, "10        Sum of edge densities (columns 10 thru 15)\n");
        fprintf(fpDesc, "11        Total system edge density (from all Ltot_bx files)\n");
        fprintf(fpDesc, "12        Total system screw density (from all Ltot_bx files)\n");
        fprintf(fpDesc, "13        Simulation timestep duration\n");

/*      
 *      The following planes and Burgers vectors correspond to the arrays
 *      burgRef and planeRef in BCC_Util.cc.
 */

        fprintf(fpDesc, "\nGlide planes are\n");       
        fprintf(fpDesc, "For Burgers vector 1/2 [ 1  1  1]\n");
        fprintf(fpDesc, "p1=( 0, -1,  1):\n");
        fprintf(fpDesc, "p2=( 1,  0, -1):\n");
        fprintf(fpDesc, "p3=(-1,  1,  0):\n");
        fprintf(fpDesc, "p4=( 2, -1, -1):\n");
        fprintf(fpDesc, "p5=(-1,  2, -1):\n");
        fprintf(fpDesc, "p6=(-1, -1,  2):\n");

        fprintf(fpDesc, "For Burgers vector 1/2 [ 1  1 -1]\n");
        fprintf(fpDesc, "p1=( 0,  1,  1):\n");
        fprintf(fpDesc, "p2=(-1,  0, -1):\n");
        fprintf(fpDesc, "p3=( 1, -1,  0):\n");
        fprintf(fpDesc, "p4=( 2, -1,  1):\n");
        fprintf(fpDesc, "p5=(-1,  2,  1):\n");
        fprintf(fpDesc, "p6=(-1, -1, -2):\n");

        fprintf(fpDesc, "For Burgers vector 1/2 [ 1 -1  1]\n");
        fprintf(fpDesc, "p1=( 0, -1, -1):\n");
        fprintf(fpDesc, "p2=(-1,  0,  1):\n");
        fprintf(fpDesc, "p3=( 1,  1,  0):\n");
        fprintf(fpDesc, "p4=( 2,  1, -1):\n");
        fprintf(fpDesc, "p5=(-1, -2, -1):\n");
        fprintf(fpDesc, "p6=(-1,  1,  2):\n");

        fprintf(fpDesc, "For Burgers vector 1/2 [-1  1  1]\n");
        fprintf(fpDesc, "p1=( 0,  1, -1):\n");
        fprintf(fpDesc, "p2=( 1,  0,  1):\n");
        fprintf(fpDesc, "p3=(-1, -1,  0):\n");
        fprintf(fpDesc, "p4=(-2, -1, -1):\n");
        fprintf(fpDesc, "p5=( 1,  2, -1):\n");
        fprintf(fpDesc, "p6=( 1, -1,  2):\n");
        fclose(fpDesc);
    }

/*
 *  Open and write the description file associated with the
 *  Ltot_b*_full_decomp flux output files if that class of output
 *  has been selected.
 */
    if (param->writeFluxFullDecomp) {

        sprintf(fileName, "%s/Ltot_bx_full_decomp.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The BCC2 Ltot_b*_full_decomp files each correspond to a specific\n");
        fprintf(fpDesc, "Burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1_*   Burgers vector 1/2 [ 1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b2_*   Burgers vector 1/2 [ 1  1 -1]\n");
        fprintf(fpDesc, "    Ltot_b3_*   Burgers vector 1/2 [ 1 -1  1]\n");
        fprintf(fpDesc, "    Ltot_b4_*   Burgers vector 1/2 [-1  1  1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_full_decomp files are identical in "
                        "format and contain\n");
        fprintf(fpDesc, "the following data for their associated burgers "
                        "vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Screw density plane p1\n");
        fprintf(fpDesc, "4         Screw density plane p2\n");
        fprintf(fpDesc, "5         Screw density plane p3\n");
        fprintf(fpDesc, "6         Screw density plane p4\n");
        fprintf(fpDesc, "7         Screw density plane p5\n");
        fprintf(fpDesc, "8         Screw density plane p6\n");

        fprintf(fpDesc, "9         Edge density plane p1\n");
        fprintf(fpDesc, "10        Edge density plane p2\n");
        fprintf(fpDesc, "11        Edge density plane p3\n");
        fprintf(fpDesc, "12        Edge density plane p4\n");
        fprintf(fpDesc, "13        Edge density plane p5\n");
        fprintf(fpDesc, "14        Edge density plane p6\n");
        fprintf(fpDesc, "15        Simulation timestep duration\n");


        fprintf(fpDesc, "\nGlide planes are\n");
/*      
 *      The following planes and Burgers vectors correspond to the arrays
 *      burgRef and planeRef in BCC_Util.cc.
 */

        fprintf(fpDesc, "\nGlide planes are\n");       
        fprintf(fpDesc, "For Burgers vector 1/2 [ 1  1  1]\n");
        fprintf(fpDesc, "p1=( 0, -1,  1):\n");
        fprintf(fpDesc, "p2=( 1,  0, -1):\n");
        fprintf(fpDesc, "p3=(-1,  1,  0):\n");
        fprintf(fpDesc, "p4=( 2, -1, -1):\n");
        fprintf(fpDesc, "p5=(-1,  2, -1):\n");
        fprintf(fpDesc, "p6=(-1, -1,  2):\n");

        fprintf(fpDesc, "For Burgers vector 1/2 [ 1  1 -1]\n");
        fprintf(fpDesc, "p1=( 0,  1,  1):\n");
        fprintf(fpDesc, "p2=(-1,  0, -1):\n");
        fprintf(fpDesc, "p3=( 1, -1,  0):\n");
        fprintf(fpDesc, "p4=( 2, -1,  1):\n");
        fprintf(fpDesc, "p5=(-1,  2,  1):\n");
        fprintf(fpDesc, "p6=(-1, -1, -2):\n");

        fprintf(fpDesc, "For Burgers vector 1/2 [ 1 -1  1]\n");
        fprintf(fpDesc, "p1=( 0, -1, -1):\n");
        fprintf(fpDesc, "p2=(-1,  0,  1):\n");
        fprintf(fpDesc, "p3=( 1,  1,  0):\n");
        fprintf(fpDesc, "p4=( 2,  1, -1):\n");
        fprintf(fpDesc, "p5=(-1, -2, -1):\n");
        fprintf(fpDesc, "p6=(-1,  1,  2):\n");

        fprintf(fpDesc, "For Burgers vector 1/2 [-1  1  1]\n");
        fprintf(fpDesc, "p1=( 0,  1, -1):\n");
        fprintf(fpDesc, "p2=( 1,  0,  1):\n");
        fprintf(fpDesc, "p3=(-1, -1,  0):\n");
        fprintf(fpDesc, "p4=(-2, -1, -1):\n");
        fprintf(fpDesc, "p5=( 1,  2, -1):\n");
        fprintf(fpDesc, "p6=( 1, -1,  2):\n");
        fclose(fpDesc);

    }  // end if (param->writeFluxFullDecomp)

/*
 *  Open and write the description file associated with the
 *  Ltot_b*_full_decomp_totals flux output files if that class
 *  of output has been selected.
 */
    if (param->writeFluxFullDecompTotals) {

        sprintf(fileName, "%s/Ltot_bx_full_decomp_totals.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The BCC2 Ltot_b*_full_decomp_totals files each "
                        "correspond to a specific\n");
        fprintf(fpDesc, "Burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1_*   Burgers vector 1/2 [ 1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b2_*   Burgers vector 1/2 [ 1  1 -1]\n");
        fprintf(fpDesc, "    Ltot_b3_*   Burgers vector 1/2 [ 1 -1  1]\n");
        fprintf(fpDesc, "    Ltot_b4_*   Burgers vector 1/2 [-1  1  1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_full_decomp_totals files are "
                        "identical in format and\n");
        fprintf(fpDesc, "contain the following data for their associated "
                        "burgers vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Total density plane p1\n");
        fprintf(fpDesc, "4         Total density plane p2\n");
        fprintf(fpDesc, "5         Total density plane p3\n");
        fprintf(fpDesc, "6         Total density plane p4\n");
        fprintf(fpDesc, "7         Total density plane p5\n");
        fprintf(fpDesc, "8         Total density plane p6\n");
        fprintf(fpDesc, "9         Total density (sum of columns 3 - 8)\n");
        fprintf(fpDesc, "10        Simulation timestep duration\n");
        fclose(fpDesc);

    }  // end if (param->writeFluxFullDecompTotals)

/*
 *  Open and write the description file associated with the
 *  Ltot_b*_simple_totals flux output files if that class
 *  of output has been selected.
 */
    if (param->writeFluxSimpleTotals) {

        sprintf(fileName, "%s/Ltot_bx_simple_totals.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The BCC2 Ltot_b*_simple_totals files each correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1_*   Burgers vector 1/2 [ 1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b2_*   Burgers vector 1/2 [ 1  1 -1]\n");
        fprintf(fpDesc, "    Ltot_b3_*   Burgers vector 1/2 [ 1 -1  1]\n");
        fprintf(fpDesc, "    Ltot_b4_*   Burgers vector 1/2 [-1  1  1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_simple_totals files are "
                        "identical in format and\n");
        fprintf(fpDesc, "contain the following data for their associated "
                        "burgers vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Total screw density\n");
        fprintf(fpDesc, "4         Total edge density\n");
        fprintf(fpDesc, "5         Simulation timestep duration\n");
        fclose(fpDesc);

    }  // end if (param->writeFluxSimpleTotals)

/*
 *  Open and write the description file associated with the fluxtot_* 
 *  flux output files.
 */
    sprintf(fileName, "%s/fluxtot_bx.%s", DIR_FLUXDATA,
            DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write flux description file %s\n", fileName);
        return;
    }

    fprintf(fpDesc, "The BCC2 fluxtot* files each correspond to a specific\n");
    fprintf(fpDesc, "burgers vector as indicated below.\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    fluxtot_b1    Burgers vector 1/2 [ 1  1  1]\n");
    fprintf(fpDesc, "    fluxtot_b2    Burgers vector 1/2 [ 1  1 -1]\n");
    fprintf(fpDesc, "    fluxtot_b3    Burgers vector 1/2 [ 1 -1  1]\n");
    fprintf(fpDesc, "    fluxtot_b4    Burgers vector 1/2 [-1  1  1]\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "The format of the fluxtot* files is similar, however\n");
    fprintf(fpDesc, "the flux from edge components on the 3 planes and\n");
    fprintf(fpDesc, "the flux from screw components on the 3 planes\n");
    fprintf(fpDesc, "(columns 4 thru 6 and 7 thru 9 respectively) differ\n");
    fprintf(fpDesc, "as seen below.\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");
    fprintf(fpDesc, "1         Plastic strain\n");
    fprintf(fpDesc, "2         Strain\n");
    fprintf(fpDesc, "3         Flux due to climb\n");
    fprintf(fpDesc, "16        Simulation timstep duration\n");
    fprintf(fpDesc, "\n");

/*      
 *      The following planes and Burgers vectors correspond to the arrays
 *      burgRef and planeRef in BCC_Util.cc.
 */

        fprintf(fpDesc, "\nGlide planes are\n");       
        fprintf(fpDesc, "For fluxtot_b1 and Burgers vector 1/2 [ 1  1  1]\n");
        fprintf(fpDesc, "4         edge ( 0, -1,  1):\n");
        fprintf(fpDesc, "5         edge ( 1,  0, -1):\n");
        fprintf(fpDesc, "6         edge (-1,  1,  0):\n");
        fprintf(fpDesc, "7         edge ( 2, -1, -1):\n");
        fprintf(fpDesc, "8         edge (-1,  2, -1):\n");
        fprintf(fpDesc, "9         edge (-1, -1,  2):\n");

        fprintf(fpDesc, "10        screw ( 0, -1,  1):\n");
        fprintf(fpDesc, "11        screw ( 1,  0, -1):\n");
        fprintf(fpDesc, "12        screw (-1,  1,  0):\n");
        fprintf(fpDesc, "13        screw ( 2, -1, -1):\n");
        fprintf(fpDesc, "14        screw (-1,  2, -1):\n");
        fprintf(fpDesc, "15        screw (-1, -1,  2):\n");

        fprintf(fpDesc, "For fluxtot_b2 and Burgers vector 1/2 [ 1  1 -1]\n");
        fprintf(fpDesc, "4         edge ( 0,  1,  1):\n");
        fprintf(fpDesc, "5         edge (-1,  0, -1):\n");
        fprintf(fpDesc, "6         edge ( 1, -1,  0):\n");
        fprintf(fpDesc, "7         edge ( 2, -1,  1):\n");
        fprintf(fpDesc, "8         edge (-1,  2,  1):\n");
        fprintf(fpDesc, "9         edge (-1, -1, -2):\n");

        fprintf(fpDesc, "10        screw ( 0,  1,  1):\n");
        fprintf(fpDesc, "11        screw (-1,  0, -1):\n");
        fprintf(fpDesc, "12        screw ( 1, -1,  0):\n");
        fprintf(fpDesc, "13        screw ( 2, -1,  1):\n");
        fprintf(fpDesc, "14        screw (-1,  2,  1):\n");
        fprintf(fpDesc, "15        screw (-1, -1, -2):\n");

        fprintf(fpDesc, "For fluxtot_b3 and Burgers vector 1/2 [ 1 -1  1]\n");
        fprintf(fpDesc, "4         edge ( 0, -1, -1):\n");
        fprintf(fpDesc, "5         edge (-1,  0,  1):\n");
        fprintf(fpDesc, "6         edge ( 1,  1,  0):\n");
        fprintf(fpDesc, "7         edge ( 2,  1, -1):\n");
        fprintf(fpDesc, "8         edge (-1, -2, -1):\n");
        fprintf(fpDesc, "9         edge (-1,  1,  2):\n");

        fprintf(fpDesc, "10        screw ( 0, -1, -1):\n");
        fprintf(fpDesc, "11        screw (-1,  0,  1):\n");
        fprintf(fpDesc, "12        screw ( 1,  1,  0):\n");
        fprintf(fpDesc, "13        screw ( 2,  1, -1):\n");
        fprintf(fpDesc, "14        screw (-1, -2, -1):\n");
        fprintf(fpDesc, "15        screw (-1,  1,  2):\n");

        fprintf(fpDesc, "For fluxtot_b4 and Burgers vector 1/2 [-1  1  1]\n");
        fprintf(fpDesc, "4         edge ( 0,  1, -1):\n");
        fprintf(fpDesc, "5         edge ( 1,  0,  1):\n");
        fprintf(fpDesc, "6         edge (-1, -1,  0):\n");
        fprintf(fpDesc, "7         edge (-2, -1, -1):\n");
        fprintf(fpDesc, "8         edge ( 1,  2, -1):\n");
        fprintf(fpDesc, "9         edge ( 1, -1,  2):\n");

        fprintf(fpDesc, "10        screw ( 0,  1, -1):\n");
        fprintf(fpDesc, "11        screw ( 1,  0,  1):\n");
        fprintf(fpDesc, "12        screw (-1, -1,  0):\n");
        fprintf(fpDesc, "13        screw (-2, -1, -1):\n");
        fprintf(fpDesc, "14        screw ( 1,  2, -1):\n");
        fprintf(fpDesc, "15        screw ( 1, -1,  2):\n");
        fclose(fpDesc);

    return;
}



/****************************************************************************
 *
 *      Function:     WriteDensFluxDesc_BCC
 *      Description:  Create a file containing a description of the contents
 *                    of the two sets of flux decomposition output files
 *                    appropriate to BCC simulations.
 *
 ****************************************************************************/
void WriteDensFluxDesc_BCC(Home_t *home)
{
    char        fileName[MAX_STRING_LEN];
    FILE       *fpDesc;
    Param_t    *param = home->param;
    static int  descFileWritten = 0;


/*
 *  Only one task (domain zero) writes this file, and only writes
 *  it one timer per simulation.  If the associated description
 *  file already exists, the contents are overwritten
 */
    if ((home->myDomain != 0) || (descFileWritten == 1)) {
        return;
    }

    descFileWritten = 1;

/*
 *  Open and write the description file associated with the Ltot_* 
 *  flux output files if requested.
 */
    if (param->writeFluxEdgeDecomp) {

        sprintf(fileName, "%s/Ltot_bx.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The BCC Ltot_bx files each correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1    Burgers vector 1/2 [ 1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b2    Burgers vector 1/2 [-1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b3    Burgers vector 1/2 [ 1 -1  1]\n");
        fprintf(fpDesc, "    Ltot_b4    Burgers vector 1/2 [ 1  1 -1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_bx files are identical in format and contain\n");
        fprintf(fpDesc, "the following data for the associated burgers vectors\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Screw density\n");
        fprintf(fpDesc, "4         Edge density 1\n");
        fprintf(fpDesc, "5         Edge density 2\n");
        fprintf(fpDesc, "6         Edge density 3\n");
        fprintf(fpDesc, "7         Sum of edge densities (columns 4 thru 6)\n");
        fprintf(fpDesc, "8         Total system edge density (from all Ltot_bx files)\n");
        fprintf(fpDesc, "9         Total system screw density (from all Ltot_bx files)\n");
        fprintf(fpDesc, "10        Simulation timestep duration\n");

        fclose(fpDesc);
    }

/*
 *  Open and write the description file associated with the
 *  Ltot_b*_full_decomp flux output files if that class of output
 *  has been selected.
 */
    if (param->writeFluxFullDecomp) {

        sprintf(fileName, "%s/Ltot_bx_full_decomp.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The BCC Ltot_b*_full_decomp files each correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1_*   Burgers vector 1/2 [ 1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b2_*   Burgers vector 1/2 [-1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b3_*   Burgers vector 1/2 [ 1 -1  1]\n");
        fprintf(fpDesc, "    Ltot_b4_*   Burgers vector 1/2 [ 1  1 -1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_full_decomp files are identical in "
                        "format and contain\n");
        fprintf(fpDesc, "the following data for their associated burgers "
                        "vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Screw density plane 1\n");
        fprintf(fpDesc, "4         Screw density plane 2\n");
        fprintf(fpDesc, "5         Screw density plane 3\n");
        fprintf(fpDesc, "6         Edge density plane 1\n");
        fprintf(fpDesc, "7         Edge density plane 2\n");
        fprintf(fpDesc, "8         Edge density plane 3\n");
        fprintf(fpDesc, "9         Simulation timestep duration\n");
        fclose(fpDesc);

    }  // end if (param->writeFluxFullDecomp)

/*
 *  Open and write the description file associated with the
 *  Ltot_b*_full_decomp_totals flux output files if that class
 *  of output has been selected.
 */
    if (param->writeFluxFullDecompTotals) {

        sprintf(fileName, "%s/Ltot_bx_full_decomp_totals.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The BCC Ltot_b*_full_decomp_totals files each "
                        "correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1_*   Burgers vector 1/2 [ 1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b2_*   Burgers vector 1/2 [-1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b3_*   Burgers vector 1/2 [ 1 -1  1]\n");
        fprintf(fpDesc, "    Ltot_b4_*   Burgers vector 1/2 [ 1  1 -1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_full_decomp_totals files are "
                        "identical in format and\n");
        fprintf(fpDesc, "contain the following data for their associated "
                        "burgers vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Total density plane 1\n");
        fprintf(fpDesc, "4         Total density plane 2\n");
        fprintf(fpDesc, "5         Total density plane 3\n");
        fprintf(fpDesc, "6         Total density (sum of columns 3 - 5)\n");
        fprintf(fpDesc, "7         Simulation timestep duration\n");
        fclose(fpDesc);

    }  // end if (param->writeFluxFullDecompTotals)

/*
 *  Open and write the description file associated with the
 *  Ltot_b*_simple_totals flux output files if that class
 *  of output has been selected.
 */
    if (param->writeFluxSimpleTotals) {

        sprintf(fileName, "%s/Ltot_bx_simple_totals.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The BCC Ltot_b*_simple_totals files each correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1_*   Burgers vector 1/2 [ 1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b2_*   Burgers vector 1/2 [-1  1  1]\n");
        fprintf(fpDesc, "    Ltot_b3_*   Burgers vector 1/2 [ 1 -1  1]\n");
        fprintf(fpDesc, "    Ltot_b4_*   Burgers vector 1/2 [ 1  1 -1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_simple_totals files are "
                        "identical in format and\n");
        fprintf(fpDesc, "contain the following data for their associated "
                        "burgers vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Total screw density\n");
        fprintf(fpDesc, "4         Total edge density\n");
        fprintf(fpDesc, "5         Simulation timestep duration\n");
        fclose(fpDesc);

    }  // end if (param->writeFluxSimpleTotals)

/*
 *  Open and write the description file associated with the fluxtot_* 
 *  flux output files.
 */
    sprintf(fileName, "%s/fluxtot_bx.%s", DIR_FLUXDATA,
            DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write flux description file %s\n", fileName);
        return;
    }

    fprintf(fpDesc, "The BCC fluxtot* files each correspond to a specific\n");
    fprintf(fpDesc, "burgers vector as indicated below.\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    fluxtot_b1    Burgers vector 1/2 [ 1  1  1]\n");
    fprintf(fpDesc, "    fluxtot_b2    Burgers vector 1/2 [-1  1  1]\n");
    fprintf(fpDesc, "    fluxtot_b3    Burgers vector 1/2 [ 1 -1  1]\n");
    fprintf(fpDesc, "    fluxtot_b4    Burgers vector 1/2 [ 1  1 -1]\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "The format of the fluxtot* files is similar, however\n");
    fprintf(fpDesc, "the flux from edge components on the 3 planes and\n");
    fprintf(fpDesc, "the flux from screw components on the 3 planes\n");
    fprintf(fpDesc, "(columns 4 thru 6 and 7 thru 9 respectively) differ\n");
    fprintf(fpDesc, "as seen below.\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");
    fprintf(fpDesc, "1         Plastic strain\n");
    fprintf(fpDesc, "2         Strain\n");
    fprintf(fpDesc, "3         Flux due to climb\n");
    fprintf(fpDesc, "10        Simulation timstep duration\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "For fluxtot_b1: burgers vector 1/2 [ 1  1  1]\n");
    fprintf(fpDesc, "4         ( 0  1 -1), [-2  1  1]\n");
    fprintf(fpDesc, "5         (-1  0  1), [ 1 -2  1]\n");
    fprintf(fpDesc, "6         ( 1 -1  0), [ 1  1 -2]\n");
    fprintf(fpDesc, "7         ( 0  1 -1), [-2  1  1]\n");
    fprintf(fpDesc, "8         (-1  0  1), [ 1 -2  1]\n");
    fprintf(fpDesc, "9         ( 1 -1  0), [ 1  1 -2]\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b2: burgers vector 1/2 [-1  1  1]\n");
    fprintf(fpDesc, "4         ( 0  1 -1), [ 2  1  1]\n");
    fprintf(fpDesc, "5         ( 1  0  1), [ 1  2 -1]\n");
    fprintf(fpDesc, "6         ( 0  1 -1), [ 1 -1  2]\n");
    fprintf(fpDesc, "7         ( 0  1 -1), [ 2  1  1]\n");
    fprintf(fpDesc, "8         ( 1  0  1), [ 1  2 -1]\n");
    fprintf(fpDesc, "9         ( 0  1 -1), [ 1 -1  2]\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b3: burgers vector 1/2 [ 1 -1  1]\n");
    fprintf(fpDesc, "4         ( 0  1  1), [ 2  1 -1]\n");
    fprintf(fpDesc, "5         ( 1  0 -1), [ 1  2  1]\n");
    fprintf(fpDesc, "6         ( 1  1  0), [-1  1  2]\n");
    fprintf(fpDesc, "7         ( 0  1  1), [ 2  1 -1]\n");
    fprintf(fpDesc, "8         ( 1  0 -1), [ 1  2  1]\n");
    fprintf(fpDesc, "9         ( 1  1  0), [-1  1  2]\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b4: burgers vector 1/2 [ 1  1 -1]\n");
    fprintf(fpDesc, "4         ( 0  1  1), [ 2 -1  1]\n");
    fprintf(fpDesc, "5         ( 1  0  1), [-1  2  1]\n");
    fprintf(fpDesc, "6         ( 1 -1  0), [ 1  1  2]\n");
    fprintf(fpDesc, "7         ( 0  1  1), [ 2 -1  1]\n");
    fprintf(fpDesc, "8         ( 1  0  1), [-1  2  1]\n");
    fprintf(fpDesc, "9         ( 1 -1  0), [ 1  1  2]\n");

    fclose(fpDesc);

    return;
}


/****************************************************************************
 *
 *      Function:     WriteDensFluxDesc_FCC
 *      Description:  Create a file containing a description of the contents
 *                    of the two sets of flux decomposition output files
 *                    appropriate to FCC simulations.
 *
 ****************************************************************************/
void WriteDensFluxDesc_FCC(Home_t *home)
{
    char        fileName[MAX_STRING_LEN];
    FILE       *fpDesc;
    Param_t    *param = home->param;
    static int  descFileWritten = 0;


/*
 *  Only one task (domain zero) writes this file, and only writes
 *  it one timer per simulation.  If the associated description
 *  file already exists, the contents are overwritten
 */
    if ((home->myDomain != 0) || (descFileWritten == 1)) {
        return;
    }

    descFileWritten = 1;

/*
 *  Open and write the description file associated with the Ltot_* 
 *  flux output files if requested.
 */
    if (param->writeFluxEdgeDecomp) {
        sprintf(fileName, "%s/Ltot_bx.%s", DIR_FLUXDATA, DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n", fileName);
            return;
        }

        fprintf(fpDesc, "The FCC Ltot_bx files each correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1    Burgers vector 1/2 [ 1  1  0]\n");
        fprintf(fpDesc, "    Ltot_b2    Burgers vector 1/2 [ 1 -1  0]\n");
        fprintf(fpDesc, "    Ltot_b3    Burgers vector 1/2 [ 1  0  1]\n");
        fprintf(fpDesc, "    Ltot_b4    Burgers vector 1/2 [ 1  0 -1]\n");
        fprintf(fpDesc, "    Ltot_b5    Burgers vector 1/2 [ 0  1  1]\n");
        fprintf(fpDesc, "    Ltot_b6    Burgers vector 1/2 [ 0  1 -1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_bx files are identical in format and contain\n");
        fprintf(fpDesc, "the following data for their associated burgers "
                        "vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Screw density\n");
        fprintf(fpDesc, "4         Edge density 1\n");
        fprintf(fpDesc, "5         Edge density 2\n");
        fprintf(fpDesc, "6         Edge density 3\n");
        fprintf(fpDesc, "7         Sum of edge densities (columns 4 thru 6)\n");
        fprintf(fpDesc, "8         Total system edge density (from all Ltot_bx files)\n");
        fprintf(fpDesc, "9         Total system screw density (from all Ltot_bx files)\n");
        fprintf(fpDesc, "10        Simulation timestep duration\n");

        fclose(fpDesc);

    }  // end if (param->writeFluxEdgeDecomp) {

/*
 *  Open and write the description file associated with the Ltot_b*_full_decomp
 *  flux output files if that class of output has been selected.
 */
    if (home->param->writeFluxFullDecomp) {

        sprintf(fileName, "%s/Ltot_bx_full_decomp.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The FCC Ltot_b*_full_decomp files each correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1_*   Burgers vector 1/2 [ 1  1  0]\n");
        fprintf(fpDesc, "    Ltot_b2_*   Burgers vector 1/2 [ 1 -1  0]\n");
        fprintf(fpDesc, "    Ltot_b3_*   Burgers vector 1/2 [ 1  0  1]\n");
        fprintf(fpDesc, "    Ltot_b4_*   Burgers vector 1/2 [ 1  0 -1]\n");
        fprintf(fpDesc, "    Ltot_b5_*   Burgers vector 1/2 [ 0  1  1]\n");
        fprintf(fpDesc, "    Ltot_b6_*   Burgers vector 1/2 [ 0  1 -1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_full_decomp files are identical in "
                        "format and contain\n");
        fprintf(fpDesc, "the following data for their associated burgers "
                        "vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Screw density plane 1\n");
        fprintf(fpDesc, "4         Screw density plane 2\n");
        fprintf(fpDesc, "5         Screw density plane 3\n");
        fprintf(fpDesc, "6         Edge density plane 1\n");
        fprintf(fpDesc, "7         Edge density plane 2\n");
        fprintf(fpDesc, "8         Edge density plane 3\n");
        fprintf(fpDesc, "9         Simulation timestep duration\n");
        fclose(fpDesc);

    }  // end if (home->param->writeFluxFulLDecomp)

/*
 *  Open and write the description file associated with the
 *  Ltot_b*_full_decomp_totals flux output files if that class
 *  of output has been selected.
 */
    if (param->writeFluxFullDecompTotals) {

        sprintf(fileName, "%s/Ltot_bx_full_decomp_totals.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The FCC Ltot_b*_full_decomp_totals files each "
                        "correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1_*   Burgers vector 1/2 [ 1  1  0]\n");
        fprintf(fpDesc, "    Ltot_b2_*   Burgers vector 1/2 [ 1 -1  0]\n");
        fprintf(fpDesc, "    Ltot_b3_*   Burgers vector 1/2 [ 1  0  1]\n");
        fprintf(fpDesc, "    Ltot_b4_*   Burgers vector 1/2 [ 1  0 -1]\n");
        fprintf(fpDesc, "    Ltot_b5_*   Burgers vector 1/2 [ 0  1  1]\n");
        fprintf(fpDesc, "    Ltot_b6_*   Burgers vector 1/2 [ 0  1 -1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_full_decomp_totals files are "
                        "identical in format and\n");
        fprintf(fpDesc, "contain the following data for their associated "
                        "burgers vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Total density plane 1\n");
        fprintf(fpDesc, "4         Total density plane 2\n");
        fprintf(fpDesc, "5         Total density plane 3\n");
        fprintf(fpDesc, "6         Total density (sum of columns 3 - 5)\n");
        fprintf(fpDesc, "7         Simulation timestep duration\n");
        fclose(fpDesc);

    }  // end if (param->writeFluxFullDecompTotals)

/*
 *  Open and write the description file associated with the
 *  Ltot_b*_simple_totals flux output files if that class
 *  of output has been selected.
 */
    if (param->writeFluxSimpleTotals) {

        sprintf(fileName, "%s/Ltot_bx_simple_totals.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The FCC Ltot_b*_simple_totals files each correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "    Ltot_b1_*   Burgers vector 1/2 [ 1  1  0]\n");
        fprintf(fpDesc, "    Ltot_b2_*   Burgers vector 1/2 [ 1 -1  0]\n");
        fprintf(fpDesc, "    Ltot_b3_*   Burgers vector 1/2 [ 1  0  1]\n");
        fprintf(fpDesc, "    Ltot_b4_*   Burgers vector 1/2 [ 1  0 -1]\n");
        fprintf(fpDesc, "    Ltot_b5_*   Burgers vector 1/2 [ 0  1  1]\n");
        fprintf(fpDesc, "    Ltot_b6_*   Burgers vector 1/2 [ 0  1 -1]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_simple_totals files are "
                        "identical in format and\n");
        fprintf(fpDesc, "contain the following data for their associated "
                        "burgers vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Total screw density\n");
        fprintf(fpDesc, "4         Total edge density\n");
        fprintf(fpDesc, "5         Simulation timestep duration\n");
        fclose(fpDesc);

    }  // end if (param->writeFluxSimpleTotals)

/*
 *  Open and write the description file associated with the fluxtot_* 
 *  flux output files.
 */
    sprintf(fileName, "%s/fluxtot_bx.%s", DIR_FLUXDATA,
            DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write flux description file %s\n", fileName);
        return;
    }

    fprintf(fpDesc, "The FCC fluxtot* files each correspond to a specific\n");
    fprintf(fpDesc, "burgers vector as indicated below.\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    fluxtot_b1    Burgers vector 1/2 [ 1  1  1]\n");
    fprintf(fpDesc, "    fluxtot_b2    Burgers vector 1/2 [-1  1  1]\n");
    fprintf(fpDesc, "    fluxtot_b3    Burgers vector 1/2 [ 1 -1  1]\n");
    fprintf(fpDesc, "    fluxtot_b4    Burgers vector 1/2 [ 1  1 -1]\n");
    fprintf(fpDesc, "    fluxtot_b5    Burgers vector 1/2 [ 0  1  1]\n");
    fprintf(fpDesc, "    fluxtot_b6    Burgers vector 1/2 [ 0  1 -1]\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "The format of the fluxtot* files is similar, however\n");
    fprintf(fpDesc, "the flux from edge components on the 3 planes and\n");
    fprintf(fpDesc, "the flux from screw components on the 3 planes\n");
    fprintf(fpDesc, "(columns 4 thru 6 and 7 thru 9 respectively) differ\n");
    fprintf(fpDesc, "as seen below.\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");
    fprintf(fpDesc, "1         Plastic strain\n");
    fprintf(fpDesc, "2         Strain\n");
    fprintf(fpDesc, "3         Flux due to climb\n");
    fprintf(fpDesc, "10        Simulation timstep duration\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "For fluxtot_b1: burgers vector 1/2 [ 1  1  0]\n");
    fprintf(fpDesc, "4         ( 1 -1  1), [-1  1  2]\n");
    fprintf(fpDesc, "5         (-1  1  1), [-1  1 -2]\n");
    fprintf(fpDesc, "6         ( 0 -0  1), [ 1  1  0]\n");
    fprintf(fpDesc, "7         ( 1 -1  1), [-1  1  2]\n");
    fprintf(fpDesc, "8         (-1  1  1), [-1  1 -2]\n");
    fprintf(fpDesc, "9         ( 1  0  1), [ 1  1  0]\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b2: burgers vector 1/2 [ 1 -1  0]\n");
    fprintf(fpDesc, "4         ( 1  1  1), [ 1  1 -2]\n");
    fprintf(fpDesc, "5         ( 1 -1 -1), [-1 -1 -2]\n");
    fprintf(fpDesc, "6         ( 0  0  1), [ 1  1  0]\n");
    fprintf(fpDesc, "7         ( 1  1  1), [ 1  1 -2]\n");
    fprintf(fpDesc, "8         ( 1 -1 -1), [-1 -1 -2]\n");
    fprintf(fpDesc, "9         ( 0  0  1), [ 1  1  0]\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b3: burgers vector 1/2 [ 1  0  1]\n");
    fprintf(fpDesc, "4         ( 1  1 -1), [ 1 -2 -1]\n");
    fprintf(fpDesc, "5         (-1  1  1), [ 1  2 -1]\n");
    fprintf(fpDesc, "6         ( 0  1  0), [ 1  0 -1]\n");
    fprintf(fpDesc, "7         ( 1  1 -1), [ 1 -2 -1]\n");
    fprintf(fpDesc, "8         (-1  1  1), [ 1  2 -1]\n");
    fprintf(fpDesc, "9         ( 0  1  0), [ 1  0 -1]\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b4: burgers vector 1/2 [ 1  0 -1]\n");
    fprintf(fpDesc, "4         ( 1  1  1), [-1  2 -1]\n");
    fprintf(fpDesc, "5         ( 1 -1  1), [ 1  2  1]\n");
    fprintf(fpDesc, "6         ( 0  1  0), [-1  0 -1]\n");
    fprintf(fpDesc, "7         ( 1  1  1), [-1  2 -1]\n");
    fprintf(fpDesc, "8         ( 1 -1  1), [ 1  2  1]\n");
    fprintf(fpDesc, "9         ( 0  1  0), [-1  0 -1]\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b5: burgers vector 1/2 [ 0  1  1]\n");
    fprintf(fpDesc, "4         ( 1  1 -1), [ 2 -1  1]\n");
    fprintf(fpDesc, "5         ( 1 -1  1), [-2 -1  1]\n");
    fprintf(fpDesc, "6         ( 1  0  0), [ 0  1  1]\n");
    fprintf(fpDesc, "7         ( 1  1 -1), [ 2 -1  1]\n");
    fprintf(fpDesc, "8         ( 1 -1  1), [-2 -1  1]\n");
    fprintf(fpDesc, "9         ( 1  0  0), [ 0  1  1]\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b6: burgers vector 1/2 [ 0  1 -1]\n");
    fprintf(fpDesc, "4         ( 1  1  1), [-2  1  1]\n");
    fprintf(fpDesc, "5         (-1  1  1), [-2 -1 -1]\n");
    fprintf(fpDesc, "6         ( 1  0  0), [ 0  1  1]\n");
    fprintf(fpDesc, "7         ( 1  1  1), [-2  1  1]\n");
    fprintf(fpDesc, "8         (-1  1  1), [-2 -1 -1]\n");
    fprintf(fpDesc, "9         ( 1  0  0), [ 0  1  1]\n");
    fprintf(fpDesc, "\n");

    fclose(fpDesc);

    return;
}


/****************************************************************************
 *
 *      Function:     WriteDensFluxDesc_HCP
 *      Description:  Create a file containing a description of the contents
 *                    of the two sets of flux decomposition output files
 *                    appropriate to HCP simulations.
 *
 ****************************************************************************/
void WriteDensFluxDesc_HCP(Home_t *home)
{
    char        fileName[MAX_STRING_LEN];
    FILE       *fpDesc;
    Param_t    *param = home->param;
    static int  descFileWritten = 0;


/*
 *  Only one task (domain zero) writes this file, and only writes
 *  it one timer per simulation.  If the associated description
 *  file already exists, the contents are overwritten
 */
    if ((home->myDomain != 0) || (descFileWritten == 1)) {
        return;
    }

    descFileWritten = 1;

/*
 *  Open and write the description file associated with the Ltot_* 
 *  flux output files if requested.
 */
    if (param->writeFluxEdgeDecomp) {
        sprintf(fileName, "%s/Ltot_bx.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The HCP Ltot_bx files each correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "   Ltot_b1  Burgers vector  a/2[-1  sqrt(3)  0]\n");
        fprintf(fpDesc, "   Ltot_b2  Burgers vector     [ a  0  0]\n");
        fprintf(fpDesc, "   Ltot_b3  Burgers vector  a/2[-1 -sqrt(3)  0]\n");
        fprintf(fpDesc, "   Ltot_b4  Burgers vector  1/2[-a  a*sqrt(3)  2c]\n");
        fprintf(fpDesc, "   Ltot_b5  Burgers vector     [ a  0  c]\n");
        fprintf(fpDesc, "   Ltot_b6  Burgers vector  1/2[ a  a*sqrt(3) -2c]\n");
        fprintf(fpDesc, "   Ltot_b7  Burgers vector  1/2[-a  a*sqrt(3) -2c]\n");
        fprintf(fpDesc, "   Ltot_b8  Burgers vector     [-a  0  c]\n");
        fprintf(fpDesc, "   Ltot_b9  Burgers vector -1/2[ a  a*sqrt(3)  2c]\n");
        fprintf(fpDesc, "   Ltot_b10 Burgers vector     [ 0  0  c]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_bx files are identical in format and contain\n");
        fprintf(fpDesc, "the following data for their associated burgers vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Screw density\n");
        fprintf(fpDesc, "4         Edge density 1\n");
        fprintf(fpDesc, "5         Edge density 2\n");
        fprintf(fpDesc, "6         Edge density 3\n");
        fprintf(fpDesc, "7         Edge density 4\n");
        fprintf(fpDesc, "8         Sum of edge densities (columns 4 thru 7)\n");
        fprintf(fpDesc, "9         Total system edge density (from all Ltot_bx files)\n");
        fprintf(fpDesc, "10        Total system screw density (from all Ltot_bx files)\n");
        fprintf(fpDesc, "          files\n");
        fprintf(fpDesc, "11        Simulation timestep duration\n");

        fclose(fpDesc);

    }  // end if (param->writeFluxEdgeDecomp)

/*
 *  Open and write the description file associated with the
 *  Ltot_b*_full_decomp flux output files if that class of output
 *  has been selected.
 */
    if (home->param->writeFluxFullDecomp) {

        sprintf(fileName, "%s/Ltot_bx_full_decomp.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The HCP Ltot_b*_full_decomp files each correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "  Ltot_b1_*  Burgers vec  a/2[-1  sqrt(3)  0]\n");
        fprintf(fpDesc, "  Ltot_b2_*  Burgers vec     [ a  0  0]\n");
        fprintf(fpDesc, "  Ltot_b3_*  Burgers vec  a/2[-1 -sqrt(3)  0]\n");
        fprintf(fpDesc, "  Ltot_b4_*  Burgers vec  1/2[-a  a*sqrt(3)  2c]\n");
        fprintf(fpDesc, "  Ltot_b5_*  Burgers vec     [ a  0  c]\n");
        fprintf(fpDesc, "  Ltot_b6_*  Burgers vec  1/2[ a  a*sqrt(3) -2c]\n");
        fprintf(fpDesc, "  Ltot_b7_*  Burgers vec  1/2[-a  a*sqrt(3) -2c]\n");
        fprintf(fpDesc, "  Ltot_b8_*  Burgers vec     [-a  0  c]\n");
        fprintf(fpDesc, "  Ltot_b9_*  Burgers vec -1/2[ a  a*sqrt(3)  2c]\n");
        fprintf(fpDesc, "  Ltot_b10_* Burgers vec     [ 0  0  c]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_full_decomp files are identical in "
                        "format and contain\n");
        fprintf(fpDesc, "the following data for their associated burgers "
                        "vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Screw density plane 1\n");
        fprintf(fpDesc, "4         Screw density plane 2\n");
        fprintf(fpDesc, "5         Screw density plane 3\n");
        fprintf(fpDesc, "6         Screw density plane 4\n");
        fprintf(fpDesc, "7         Screw density for all other planes\n");
        fprintf(fpDesc, "8         Edge density plane 1\n");
        fprintf(fpDesc, "9         Edge density plane 2\n");
        fprintf(fpDesc, "10        Edge density plane 3\n");
        fprintf(fpDesc, "11        Edge density plane 4\n");
        fprintf(fpDesc, "12        Edge density for all other planes\n");
        fprintf(fpDesc, "13        Simulation timestep duration\n");
        fclose(fpDesc);

    }  // end if (home->param->writeFluxFulLDecomp)

/*
 *  Open and write the description file associated with the
 *  Ltot_b*_full_decomp_totals flux output files if that class
 *  of output has been selected.
 */
    if (param->writeFluxFullDecompTotals) {

        sprintf(fileName, "%s/Ltot_bx_full_decomp_totals.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The HCP Ltot_b*_full_decomp_totals files each "
                        "correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "  Ltot_b1_*  Burgers vec  a/2[-1  sqrt(3)  0]\n");
        fprintf(fpDesc, "  Ltot_b2_*  Burgers vec     [ a  0  0]\n");
        fprintf(fpDesc, "  Ltot_b3_*  Burgers vec  a/2[-1 -sqrt(3)  0]\n");
        fprintf(fpDesc, "  Ltot_b4_*  Burgers vec  1/2[-a  a*sqrt(3)  2c]\n");
        fprintf(fpDesc, "  Ltot_b5_*  Burgers vec     [ a  0  c]\n");
        fprintf(fpDesc, "  Ltot_b6_*  Burgers vec  1/2[ a  a*sqrt(3) -2c]\n");
        fprintf(fpDesc, "  Ltot_b7_*  Burgers vec  1/2[-a  a*sqrt(3) -2c]\n");
        fprintf(fpDesc, "  Ltot_b8_*  Burgers vec     [-a  0  c]\n");
        fprintf(fpDesc, "  Ltot_b9_*  Burgers vec -1/2[ a  a*sqrt(3)  2c]\n");
        fprintf(fpDesc, "  Ltot_b10_* Burgers vec     [ 0  0  c]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_full_decomp_totals files are "
                        "identical in format and\n");
        fprintf(fpDesc, "contain the following data for their associated "
                        "burgers vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Total density plane 1\n");
        fprintf(fpDesc, "4         Total density plane 2\n");
        fprintf(fpDesc, "5         Total density plane 3\n");
        fprintf(fpDesc, "6         Total density plane 4\n");
        fprintf(fpDesc, "7         Total density for all other planes\n");
        fprintf(fpDesc, "8         Total density (sum of columns 3 - 7)\n");
        fprintf(fpDesc, "9         Simulation timestep duration\n");
        fclose(fpDesc);

    }  // end if (param->writeFluxFullDecompTotals)

/*
 *  Open and write the description file associated with the
 *  Ltot_b*_simple_totals flux output files if that class of
 *  output has been selected.
 */
    if (param->writeFluxSimpleTotals) {

        sprintf(fileName, "%s/Ltot_bx_simple_totals.%s", DIR_FLUXDATA,
                DESCRIPTION_FILE_SUFFIX);

        if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
            printf("WARNING: Unable to write flux description file %s\n",
                   fileName);
            return;
        }

        fprintf(fpDesc, "The HCP Ltot_b*_simple_totals files each correspond to a specific\n");
        fprintf(fpDesc, "burgers vector as indicated below.\n");
        fprintf(fpDesc, "\n");
        fprintf(fpDesc, "  Ltot_b1_*  Burgers vec  a/2[-1  sqrt(3)  0]\n");
        fprintf(fpDesc, "  Ltot_b2_*  Burgers vec     [ a  0  0]\n");
        fprintf(fpDesc, "  Ltot_b3_*  Burgers vec  a/2[-1 -sqrt(3)  0]\n");
        fprintf(fpDesc, "  Ltot_b4_*  Burgers vec  1/2[-a  a*sqrt(3)  2c]\n");
        fprintf(fpDesc, "  Ltot_b5_*  Burgers vec     [ a  0  c]\n");
        fprintf(fpDesc, "  Ltot_b6_*  Burgers vec  1/2[ a  a*sqrt(3) -2c]\n");
        fprintf(fpDesc, "  Ltot_b7_*  Burgers vec  1/2[-a  a*sqrt(3) -2c]\n");
        fprintf(fpDesc, "  Ltot_b8_*  Burgers vec     [-a  0  c]\n");
        fprintf(fpDesc, "  Ltot_b9_*  Burgers vec -1/2[ a  a*sqrt(3)  2c]\n");
        fprintf(fpDesc, "  Ltot_b10_* Burgers vec     [ 0  0  c]\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "All the Ltot_b*_simple_totals files are "
                        "identical in format and\n");
        fprintf(fpDesc, "contain the following data for their associated "
                        "burgers vectors.\n");
        fprintf(fpDesc, "\n");

        fprintf(fpDesc, "Column    Description\n");
        fprintf(fpDesc, "------    -------------------------------------------\n");
        fprintf(fpDesc, "1         Plastic strain\n");
        fprintf(fpDesc, "2         Strain\n");
        fprintf(fpDesc, "3         Total screw density\n");
        fprintf(fpDesc, "4         Total edge density\n");
        fprintf(fpDesc, "5         Simulation timestep duration\n");
        fclose(fpDesc);

    }  // end if (param->writeFluxSimpleTotals)

/*
 *  Open and write the description file associated with the fluxtot_* 
 *  flux output files.
 */
    sprintf(fileName, "%s/fluxtot_bx.%s", DIR_FLUXDATA,
            DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write flux description file %s\n", fileName);
        return;
    }

    fprintf(fpDesc, "The HCP fluxtot* files each correspond to a specific\n");
    fprintf(fpDesc, "burgers vector as indicated below.\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    fluxtot_b1    Burgers vector  a/2[-1  sqrt(3)  0]\n");
    fprintf(fpDesc, "    fluxtot_b2    Burgers vector     [ a  0  0]\n");
    fprintf(fpDesc, "    fluxtot_b3    Burgers vector  a/2[-1 -sqrt(3)  0]\n");
    fprintf(fpDesc, "    fluxtot_b4    Burgers vector  1/2[-a  a*sqrt(3)  2c]\n");
    fprintf(fpDesc, "    fluxtot_b5    Burgers vector     [ a  0  c]\n");
    fprintf(fpDesc, "    fluxtot_b6    Burgers vector  1/2[ a  a*sqrt(3) -2c]\n");
    fprintf(fpDesc, "    fluxtot_b7    Burgers vector  1/2[-a  a*sqrt(3) -2c]\n");
    fprintf(fpDesc, "    fluxtot_b8    Burgers vector     [-a  0  c]\n");
    fprintf(fpDesc, "    fluxtot_b9    Burgers vector -1/2[ a  a*sqrt(3)  2c]\n");
    fprintf(fpDesc, "    fluxtot_b10   Burgers vector     [ 0  0  c]\n");

    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "The format of the fluxtot* files is similar, however\n");
    fprintf(fpDesc, "the flux from edge components on the 4 planes and\n");
    fprintf(fpDesc, "the flux from screw components on the 4 planes\n");
    fprintf(fpDesc, "(columns 4 thru 7 and 8 thru 11 respectively) differ\n");
    fprintf(fpDesc, "as seen below.\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");
    fprintf(fpDesc, "1         Plastic strain\n");
    fprintf(fpDesc, "2         Strain\n");
    fprintf(fpDesc, "3         Flux due to climb\n");
    fprintf(fpDesc, "12        Simulation timstep duration\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "For fluxtot_b1:  burgers vector  a/2[-1  sqrt(3)  0]\n");
    fprintf(fpDesc, "4         Plane ( 0          0  1)\n");
    fprintf(fpDesc, "5         Plane ( sqrt(3)    1  0)\n");
    fprintf(fpDesc, "6         Plane ( c*sqrt(3)  c -a*sqrt(3))\n");
    fprintf(fpDesc, "7         Plane ( c*sqrt(3)  c  a*sqrt(3))\n");
    fprintf(fpDesc, "8         Plane ( 0          0  c)\n");
    fprintf(fpDesc, "9         Plane ( sqrt(3)    1  0)\n");
    fprintf(fpDesc, "10        Plane ( c*sqrt(3)  c -a*sqrt(3))\n");
    fprintf(fpDesc, "11        Plane ( c*sqrt(3)  c  a*sqrt(3))\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b2:  burgers vector     [ a  0  0]\n");
    fprintf(fpDesc, "4         Plane ( 0   0  1)\n");
    fprintf(fpDesc, "5         Plane ( 0   1  0)\n");
    fprintf(fpDesc, "6         Plane ( 0  2c  a*sqrt(3))\n");
    fprintf(fpDesc, "7         Plane ( 0 -2c  a*sqrt(3))\n");
    fprintf(fpDesc, "8         Plane ( 0   0  1)\n");
    fprintf(fpDesc, "9         Plane ( 0   1  0)\n");
    fprintf(fpDesc, "10        Plane ( 0  2c  a*sqrt(3))\n");
    fprintf(fpDesc, "11        Plane ( 0 -2c  a*sqrt(3))\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b3:  burgers vector  a/2[-1 -sqrt(3)  0]\n");
    fprintf(fpDesc, "4         Plane ( 0           0  1)\n");
    fprintf(fpDesc, "5         Plane ( sqrt(3)    -1  0)\n");
    fprintf(fpDesc, "6         Plane ( c*sqrt(3)  -c  a*sqrt(3))\n");
    fprintf(fpDesc, "7         Plane (-c*sqrt(3)   c  a*sqrt(3))\n");
    fprintf(fpDesc, "8         Plane ( 0           0  1)\n");
    fprintf(fpDesc, "9         Plane ( sqrt(3)    -1  0)\n");
    fprintf(fpDesc, "10        Plane ( c*sqrt(3)  -c  a*sqrt(3))\n");
    fprintf(fpDesc, "11        Plane (-c*sqrt(3)   c  a*sqrt(3))\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b4:  burgers vector  1/2[-a  a*sqrt(3)  2c]\n");
    fprintf(fpDesc, "4         Plane ( c          -c*sqrt(3)  2a)\n");
    fprintf(fpDesc, "5         Plane ( sqrt(3)     1          0)\n");
    fprintf(fpDesc, "6         Plane ( 0          -2c         a*sqrt(3))\n");
    fprintf(fpDesc, "7         Plane ( c*sqrt(3)  -c          a*sqrt(3))\n");
    fprintf(fpDesc, "8         Plane ( c          -c*sqrt(3)  2a)\n");
    fprintf(fpDesc, "9         Plane ( sqrt(3)     1          0)\n");
    fprintf(fpDesc, "10        Plane ( 0          -2c         a*sqrt(3))\n");
    fprintf(fpDesc, "11        Plane ( c*sqrt(3)  -c          a*sqrt(3))\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b5:  burgers vector     [ a  0  c]\n");
    fprintf(fpDesc, "4         Plane ( c           0  a)\n");
    fprintf(fpDesc, "5         Plane ( 0           1  0)\n");
    fprintf(fpDesc, "6         Plane (-c*sqrt(3)   c  a*sqrt(3))\n");
    fprintf(fpDesc, "7         Plane ( c*sqrt(3)   c  a*sqrt(3))\n");
    fprintf(fpDesc, "8         Plane ( c           0  a)\n");
    fprintf(fpDesc, "9         Plane ( 0           1  0)\n");
    fprintf(fpDesc, "10        Plane (-c*sqrt(3)   c  a*sqrt(3))\n");
    fprintf(fpDesc, "11        Plane ( c*sqrt(3)   c  a*sqrt(3))\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b6:  burgers vector  1/2[ a  a*sqrt(3) -2c]\n");
    fprintf(fpDesc, "4         Plane ( c          c*sqrt(3) 2a)\n");
    fprintf(fpDesc, "5         Plane ( a*sqrt(3) -a         0)\n");
    fprintf(fpDesc, "6         Plane ( c*sqrt(3)  c         a*sqrt(3))\n");
    fprintf(fpDesc, "7         Plane ( 0          2c        a*sqrt(3))\n");
    fprintf(fpDesc, "8         Plane ( c          c*sqrt(3) 2a)\n");
    fprintf(fpDesc, "9         Plane ( a*sqrt(3) -a         0)\n");
    fprintf(fpDesc, "10        Plane ( c*sqrt(3)  c         a*sqrt(3))\n");
    fprintf(fpDesc, "11        Plane ( 0          2c        a*sqrt(3))\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b7:  burgers vector  1/2[-a  a*sqrt(3) -2c]\n");
    fprintf(fpDesc, "4         Plane (-c           c*sqrt(3) -2a)\n");
    fprintf(fpDesc, "5         Plane ( sqrt(3)     1          0)\n");
    fprintf(fpDesc, "6         Plane ( 0           2c         a*sqrt(3))\n");
    fprintf(fpDesc, "7         Plane ( 0           c          a*sqrt(3))\n");
    fprintf(fpDesc, "8         Plane ( -c*sqrt(3)  c*sqrt(3) -2a)\n");
    fprintf(fpDesc, "9         Plane ( -c          1          0)\n");
    fprintf(fpDesc, "10        Plane ( 0           2c         a*sqrt(3))\n");
    fprintf(fpDesc, "11        Plane ( -c*sqrt(3)  c          a*sqrt(3))\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b8:  burgers vector     [-a  0  c]\n");
    fprintf(fpDesc, "4         Plane ( c          0  a)\n");
    fprintf(fpDesc, "5         Plane ( 0          1  0)\n");
    fprintf(fpDesc, "6         Plane ( c*sqrt(3)  c  a*sqrt(3))\n");
    fprintf(fpDesc, "7         Plane ( c*sqrt(3) -c  a*sqrt(3))\n");
    fprintf(fpDesc, "8         Plane ( c          0  a)\n");
    fprintf(fpDesc, "9         Plane ( 0          1  0)\n");
    fprintf(fpDesc, "10        Plane ( c*sqrt(3)  c  a*sqrt(3))\n");
    fprintf(fpDesc, "11        Plane ( c*sqrt(3) -c  a*sqrt(3))\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b9:  burgers vector -1/2[ a  a*sqrt(3)  2c]\n");
    fprintf(fpDesc, "4         Plane ( c          c*sqrt(3) -2a)\n");
    fprintf(fpDesc, "5         Plane ( sqrt(3)   -1          0)\n");
    fprintf(fpDesc, "6         Plane ( 0         -2c         a*sqrt(3))\n");
    fprintf(fpDesc, "7         Plane ( c*sqrt(3) -c         -a*sqrt(3))\n");
    fprintf(fpDesc, "8         Plane ( c          c*sqrt(3) -2a)\n");
    fprintf(fpDesc, "9         Plane ( sqrt(3)   -1          0)\n");
    fprintf(fpDesc, "10        Plane ( 0         -2c         a*sqrt(3))\n");
    fprintf(fpDesc, "11        Plane ( c*sqrt(3)  c         -a*sqrt(3))\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "For fluxtot_b10: burgers vector     [ 0  0  c]\n");
    fprintf(fpDesc, "4         Plane ( sqrt(3)  1  0)\n");
    fprintf(fpDesc, "5         Plane ( sqrt(3) -1  0)\n");
    fprintf(fpDesc, "6         Plane ( 0        1  0)\n");
    fprintf(fpDesc, "7         ** Zero: only 3 planes for this burgers vec\n");
    fprintf(fpDesc, "8         Plane ( sqrt(3)  1  0)\n");
    fprintf(fpDesc, "9         Plane ( sqrt(3) -1  0)\n");
    fprintf(fpDesc, "10        Plane ( 0        1  0)\n");
    fprintf(fpDesc, "11        ** Zero: only 3 planes for this burgers vec\n");
    fprintf(fpDesc, "\n");

    fclose(fpDesc);

    return;
}


/****************************************************************************
 *
 *      Function:     WriteDensFluxDesc_RhomboVa
 *      Description:  Create a file containing a description of the contents
 *                    of the two sets of flux decomposition output files
 *                    appropriate to simulations of rhombohedral vanadium.
 *
 ****************************************************************************/
void WriteDensFluxDesc_RhomboVa(Home_t *home)
{
    char       fileName[MAX_STRING_LEN];
    FILE       *fpDesc;
    static int descFileWritten = 0;


/*
 *  Only one task (domain zero) writes this file, and only writes
 *  it one timer per simulation.  If the associated description
 *  file already exists, the contents are overwritten
 */
    if ((home->myDomain != 0) || (descFileWritten == 1)) {
        return;
    }

    descFileWritten = 1;

/*
 *  Open and write the description file associated with the Ltot_* 
 *  flux output files.
 */
    sprintf(fileName, "%s/Ltot_bx.%s", DIR_FLUXDATA, DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write flux description file %s\n", fileName);
        return;
    }

    fprintf(fpDesc, "The rhombohedral vanadium Ltot_bx files each correspond "
                    "to a specific burgers\n");
    fprintf(fpDesc, "vector as indicated below.\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    Ltot_b1    Burgers vector  [1.45 1.45 1.45]\n");
    fprintf(fpDesc, "    Ltot_b2    Burgers vector  [-1.1870 1.3185 1.3185]\n");
    fprintf(fpDesc, "    Ltot_b3    Burgers vector  [1.3185 -1.1870 1.3185]\n");
    fprintf(fpDesc, "    Ltot_b4    Burgers vector  [1.3185 1.3185 -1.1870]\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "Each of the burgers vectors has three associated glide\n");
    fprintf(fpDesc, "planes for which flux is defined.  The glide planes\n");
    fprintf(fpDesc, "are as follows:\n");
    fprintf(fpDesc, "    Burgers vector 1: [1.45 1.45 1.45]\n");
    fprintf(fpDesc, "        Plane 1: (-1  1  0) \n");
    fprintf(fpDesc, "        Plane 2: ( 1  0 -1) \n");
    fprintf(fpDesc, "        Plane 3: ( 0 -1  1) \n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    Burgers vector 2: [-1.1870 1.3185 1.3185]\n");
    fprintf(fpDesc, "        Plane 1: (  0 -1   1) \n");
    fprintf(fpDesc, "        Plane 2: ( 10 -1   10) \n");
    fprintf(fpDesc, "        Plane 3: (-10 -10  1) \n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    Burgers vector 3: [1.3185 -1.1870 1.3185]\n");
    fprintf(fpDesc, "        Plane 1: (  1  0  -1) \n");
    fprintf(fpDesc, "        Plane 2: (  1 -10 -10) \n");
    fprintf(fpDesc, "        Plane 3: (-10 -10  1) \n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    Burgers vector 4: [1.3185 1.3185 -1.1870]\n");
    fprintf(fpDesc, "        Plane 1: ( -1  1   0) \n");
    fprintf(fpDesc, "        Plane 2: (  1 -10 -10) \n");
    fprintf(fpDesc, "        Plane 3: ( 10 -1   10) \n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "All the Ltot_bx files are identical in format and contain\n");
    fprintf(fpDesc, "the following data for their associated burgers "
                    "vectors.\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");
    fprintf(fpDesc, "1         Plastic strain\n");
    fprintf(fpDesc, "2         Strain\n");
    fprintf(fpDesc, "3         Screw density for plane 1\n");
    fprintf(fpDesc, "4         75 degree (trench) density for plane 1\n");
    fprintf(fpDesc, "5         Screw density for plane 2\n");
    fprintf(fpDesc, "6         75 degree (trench) density for plane 2\n");
    fprintf(fpDesc, "7         Screw density for plane 3\n");
    fprintf(fpDesc, "8         75 degree (trench) density for plane 3\n");
    fprintf(fpDesc, "9         Simulation timestep duration\n");

    fclose(fpDesc);

/*
 *  Open and write the description file associated with the fluxtot_* 
 *  flux output files.
 */
    sprintf(fileName, "%s/fluxtot_bx.%s", DIR_FLUXDATA,
            DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write flux description file %s\n", fileName);
        return;
    }

    fprintf(fpDesc, "The rhombohedral vanadium fluxtot* files each correspond "
                    "to a specific\n");
    fprintf(fpDesc, "burgers vector as indicated below.\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    fluxtot_b1 Burgers vector [1.45 1.45 1.45]\n");
    fprintf(fpDesc, "    fluxtot_b1 Burgers vector [-1.1870 1.3185 1.3185]\n");
    fprintf(fpDesc, "    fluxtot_b1 Burgers vector [1.3185 -1.1870 1.3185]\n");
    fprintf(fpDesc, "    fluxtot_b1 Burgers vector [1.3185 1.3185 -1.1870]\n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "Each of the burgers vectors has three associated glide\n");
    fprintf(fpDesc, "planes for which flux is defined.  The glide planes\n");
    fprintf(fpDesc, "are as follows:\n");
    fprintf(fpDesc, "    Burgers vector 1: [1.45 1.45 1.45]\n");
    fprintf(fpDesc, "        Plane 1: (-1  1  0) \n");
    fprintf(fpDesc, "        Plane 2: ( 1  0 -1) \n");
    fprintf(fpDesc, "        Plane 3: ( 0 -1  1) \n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    Burgers vector 2: [-1.1870 1.3185 1.3185]\n");
    fprintf(fpDesc, "        Plane 1: (  0 -1   1) \n");
    fprintf(fpDesc, "        Plane 2: ( 10 -1   10) \n");
    fprintf(fpDesc, "        Plane 3: (-10 -10  1) \n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    Burgers vector 3: [1.3185 -1.1870 1.3185]\n");
    fprintf(fpDesc, "        Plane 1: (  1  0  -1) \n");
    fprintf(fpDesc, "        Plane 2: (  1 -10 -10) \n");
    fprintf(fpDesc, "        Plane 3: (-10 -10  1) \n");
    fprintf(fpDesc, "\n");
    fprintf(fpDesc, "    Burgers vector 4: [1.3185 1.3185 -1.1870]\n");
    fprintf(fpDesc, "        Plane 1: ( -1  1   0) \n");
    fprintf(fpDesc, "        Plane 2: (  1 -10 -10) \n");
    fprintf(fpDesc, "        Plane 3: ( 10 -1   10) \n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "The format of all fluxtot* files is identical and\n");
    fprintf(fpDesc, "contain the following data for their associated\n");
    fprintf(fpDesc, "burgers vectors.\n");
    fprintf(fpDesc, "\n");

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");
    fprintf(fpDesc, "1         Plastic strain\n");
    fprintf(fpDesc, "2         Strain\n");
    fprintf(fpDesc, "3         Climb flux for plane 1\n");
    fprintf(fpDesc, "4         Screw flux for plane 1\n");
    fprintf(fpDesc, "5         75 degree trench flux for plane 1\n");
    fprintf(fpDesc, "6         Climb flux for plane 2\n");
    fprintf(fpDesc, "7         Screw flux for plane 2\n");
    fprintf(fpDesc, "8         75 degree trench flux for plane 2\n");
    fprintf(fpDesc, "9         Climb flux for plane 3\n");
    fprintf(fpDesc, "10        Screw flux for plane 3\n");
    fprintf(fpDesc, "11        75 degree trench flux for plane 3\n");
    fprintf(fpDesc, "12        Simulation timesetp duration\n");

    fclose(fpDesc);

    return;
}
