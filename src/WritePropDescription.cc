/****************************************************************************
 *
 *      Function:     WritePropDescription
 *      Description:  For each of the property output files, we want
 *                    to create an associated file containing a description
 *                    of the file format (i.e. what each field of the
 *                    output represents).  This module contains the
 *                    functions to create those description files.
 *
 *     Includes public functions:
 *         WriteAllepsFileDesc()
 *         WriteDensityFileDesc()
 *         WriteDensityDeltaFileDesc()
 *         WriteEpsdotFileDesc()
 *         WriteStressPlasticStrainFileDesc()
 *         WriteStressTotalStrainFileDesc()
 *         WriteTimePlasticStrainFileDesc()
 *         WriteTotalStressFileDesc()
 *
 ****************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "WriteProp.h"
#include "Util.h"
#include "Mobility.h"


/****************************************************************************
 *
 *      Function:     WriteDensityFileDesc
 *      Description:  Create a file containing a description of the contents
 *                    of the 'density' output file.
 *
 *     Parameters:
 *         IN:     baseFileName  Base name of the output file for which
 *                               a description file is to be created.  The
 *                               description file name will be the base
 *                               name with the DESCRIPTION_FILE_SUFFIX
 *                               appended.
 *
 ****************************************************************************/
void WriteDensityFileDesc(Home_t *home, char *baseFileName)
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

    sprintf(fileName, "%s.%s", baseFileName, DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write output description file %s\n",
               fileName);
        return;
    }

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");
    fprintf(fpDesc, "1         Plastic strain\n");
    fprintf(fpDesc, "2         Strain\n");
    fprintf(fpDesc, "3         Dislocation density\n");
    fprintf(fpDesc, "4         Deleted dislocation density\n");
    fprintf(fpDesc, "5         Average dislocation velocity\n");
    fprintf(fpDesc, "6         Std deviation of dislocation velocity\n");
    fprintf(fpDesc, "7         File version number\n");
 
    switch (home->param->materialType) {
        case MAT_TYPE_BCC:
            fprintf(fpDesc, "8         Burgers vectors [1 1 1] [-1 -1 -1]\n");
            fprintf(fpDesc, "9         Burgers vectors [-1 1 1] [1 -1 -1]\n");
            fprintf(fpDesc, "10        Burgers vectors [1 -1 1] [-1 1 -1]\n");
            fprintf(fpDesc, "11        Burgers vectors [1 1 -1] [-1 -1 1]\n");
            fprintf(fpDesc, "12        Burgers vectors [1 0 0] [-1 0 0]\n");
            fprintf(fpDesc, "                          [0 1 0] [0 -1 0]\n");
            fprintf(fpDesc, "                          [0 0 1] [0 0 -1]\n");
            break;
        case MAT_TYPE_FCC:
            fprintf(fpDesc, "8         Burgers vectors [1 1 0] [-1 -1 0]\n");
            fprintf(fpDesc, "9         Burgers vectors [-1 1 0] [1 -1 0]\n");
            fprintf(fpDesc, "10        Burgers vectors [1 0 1] [-1 0 -1]\n");
            fprintf(fpDesc, "11        Burgers vectors [-1 0 1] [1 0 -1]\n");
            fprintf(fpDesc, "12        Burgers vectors [0 1 1] [0 -1 -1]\n");
            fprintf(fpDesc, "13        Burgers vectors [0 -1 1] [0 1 -1]\n");
            fprintf(fpDesc, "14        All other burgers vectors\n");
            break;
        case MAT_TYPE_HCP:
            fprintf(fpDesc, "8         Burgers vectors a/2[-1 sqrt(3) 0]\n");
            fprintf(fpDesc, "                          a/2[1 -sqrt(3) 0]\n");
            fprintf(fpDesc, "9         Burgers vectors [a 0 0] [-a 0 0]\n");
            fprintf(fpDesc, "10        Burgers vectors a/2[1 sqrt(3) 0]\n");
            fprintf(fpDesc, "                          a/2[-1 -sqrt(3) 0]\n");
            fprintf(fpDesc, "11        Burgers vectors 1/2[-a a*sqrt(3) 2c]\n");
            fprintf(fpDesc, "                          1/2[a -a*sqrt(3) -2c]\n");
            fprintf(fpDesc, "12        Burgers vectors 1/2[a 0 c]\n");
            fprintf(fpDesc, "                          1/2[-a 0 -c]\n");
            fprintf(fpDesc, "13        Burgers vectors 1/2[a a*sqrt(3) -2c]\n");
            fprintf(fpDesc, "                          1/2[-a -a*sqrt(3) 2c]\n");
            fprintf(fpDesc, "14        Burgers vectors 1/2[-a a*sqrt(3) -2c]\n");
            fprintf(fpDesc, "                          1/2[a -a*sqrt(3) 2c]\n");
            fprintf(fpDesc, "15        Burgers vectors 1/2[-a 0 c]\n");
            fprintf(fpDesc, "                          1/2[a 0 -c]\n");
            fprintf(fpDesc, "16        Burgers vectors 1/2[a a*sqrt(3) 2]\n");
            fprintf(fpDesc, "                          1/2[-a -a*sqrt(3) -2]\n");
            fprintf(fpDesc, "17        Burgers vectors [0 0 c] [0 0 -c]\n");
            fprintf(fpDesc, "18        All other burgers vectors\n");
            break;
    }  //  end switch()


    fclose(fpDesc);

    return;
}

//-----------------------------------------------------------------------------------

void WriteDensityFileGnu(Home_t *home, char *baseFileName)
{
    static int gnuFileWritten = 0;

/*
 *  Only one task (domain zero) writes this file, and only writes
 *  it one timer per simulation.  If the associated description
 *  file already exists, the contents are overwritten
 */
    if ((home->myDomain != 0) || (gnuFileWritten == 1))
        return;

    gnuFileWritten = 1;

    char  fname[MAX_STRING_LEN];

    sprintf(fname, "%s.%s", baseFileName, GNU_FILE_SUFFIX);

    FILE *fp = fopen(fname, "w");

    if (fp)
    {
        fprintf(fp,"#set term x11\n" );
        fprintf(fp,"set term png size 640,480\n" );
        fprintf(fp,"set output 'density.png'\n" );
        fprintf(fp,"\n" );
        fprintf(fp,"set title  'ParaDiS Dislocation Density'\n" );
        fprintf(fp,"set xlabel 'Strain'  \n" );
        fprintf(fp,"set ylabel 'Dislocation Density'  \n" );
        fprintf(fp,"unset key\n" );
        fprintf(fp,"\n" );
        fprintf(fp,"set format x '%%0.0e'\n" );
        fprintf(fp,"set format y '%%0.0e'\n" );
        fprintf(fp,"\n" );
        fprintf(fp,"plot 'density' u 2:3 w lp lc rgb '#ad0000' lt 1 lw 1 pt 7 ps 0.7\n" );
        fprintf(fp,"\n" );

       fclose(fp);
    }
    else
        printf("WARNING: Unable to write output gnu command file %s\n", fname);

    return;
}


/****************************************************************************
 *
 *      Function:     WriteTimePlasticStrainFileDesc
 *      Description:  Create a file containing a description of the contents
 *                    of the 'time_Plastic_strain' output file.
 *
 *     Parameters:
 *         IN:     baseFileName  Base name of the output file for which
 *                               a description file is to be created.  The
 *                               description file name will be the base
 *                               name with the DESCRIPTION_FILE_SUFFIX
 *                               appended.
 *
 ****************************************************************************/
void WriteTimePlasticStrainFileDesc(Home_t *home, char *baseFileName)
{
    static int descFileWritten = 0;

/*
 *  Only one task (domain zero) writes this file, and only writes
 *  it one timer per simulation.  If the associated description
 *  file already exists, the contents are overwritten
 */
    if ((home->myDomain != 0) || (descFileWritten == 1))
        return;

    descFileWritten = 1;

    char  fname[MAX_STRING_LEN];

    sprintf(fname, "%s.%s", baseFileName, DESCRIPTION_FILE_SUFFIX);

    FILE *fp = fopen(fname, "w");

    if (fp)
    {
       fprintf(fp, "Column    Description\n");
       fprintf(fp, "------    -------------------------------------------\n");
       fprintf(fp, "1         Elapsed simulation time\n");
       fprintf(fp, "2         plastic strain\n");

       fclose(fp);
    }
    else
        printf("WARNING: Unable to write output description file %s\n", fname);

    return;
}

#ifdef CALCENERGY
/****************************************************************************
 *
 *      Function:     WriteTimeEnergyFileDesc
 *      Description:  Create a file containing a description of the contents
 *                    of the 'time_Energy' output file.
 *
 *     Parameters:
 *         IN:     baseFileName  Base name of the output file for which
 *                               a description file is to be created.  The
 *                               description file name will be the base
 *                               name with the DESCRIPTION_FILE_SUFFIX
 *                               appended.
 *
 ****************************************************************************/
void WriteTimeEnergyFileDesc(Home_t *home, char *baseFileName)
{
    static int descFileWritten = 0;

/*
 *  Only one task (domain zero) writes this file, and only writes
 *  it one timer per simulation.  If the associated description
 *  file already exists, the contents are overwritten
 */
    if ((home->myDomain != 0) || (descFileWritten == 1))
        return;

    descFileWritten = 1;

    char  fname[MAX_STRING_LEN];

    sprintf(fname, "%s.%s", baseFileName, DESCRIPTION_FILE_SUFFIX);

    FILE *fp = fopen(fname, "w");

    if (fp)
    {
       fprintf(fp, "Column    Description\n");
       fprintf(fp, "------    -------------------------------------------\n");
       fprintf(fp, "1         Elapsed simulation time\n");
       fprintf(fp, "2         Total energy\n");

       fclose(fp);
    }
    else
        printf("WARNING: Unable to write output description file %s\n", fname);

    return;
}
#endif

//-----------------------------------------------------------------------------------
void WriteTimePlasticStrainFileGnu(Home_t *home, char *baseFileName)
{
    static int gnuFileWritten = 0;

/*
 *  Only one task (domain zero) writes this file, and only writes
 *  it one timer per simulation.  If the associated description
 *  file already exists, the contents are overwritten
 */
    if ((home->myDomain != 0) || (gnuFileWritten == 1))
        return;

    gnuFileWritten = 1;

    char  fname[MAX_STRING_LEN];

    sprintf(fname, "%s.%s", baseFileName, GNU_FILE_SUFFIX);

    FILE *fp = fopen(fname, "w");

    if (fp)
    {
        fprintf(fp,"#set term x11\n" );
        fprintf(fp,"set term png size 640,480\n" );
        fprintf(fp,"set output 'time_Plastic_strain.png'\n" );
        fprintf(fp,"\n" );
        fprintf(fp,"set title  'ParaDiS Simulation Time / Plastic Strain'\n" );
        fprintf(fp,"set xlabel 'Simulation Time'  \n" );
        fprintf(fp,"set ylabel 'Plastic Strain'  \n" );
        fprintf(fp,"unset key\n" );
        fprintf(fp,"\n" );
        fprintf(fp,"set format x '%%0.0e'\n" );
        fprintf(fp,"set format y '%%0.0e'\n" );
        fprintf(fp,"\n" );
        fprintf(fp,"plot 'time_Plastic_strain' u 1:2 w lp lc rgb '#ad0000' lt 1 lw 1 pt 7 ps 0.7\n" );
        fprintf(fp,"\n" );

       fclose(fp);
    }
    else
        printf("WARNING: Unable to write output gnu command file %s\n", fname);

    return;
}


/****************************************************************************
 *
 *      Function:     WriteStressPlasticStrainFileDesc
 *      Description:  Create a file containing a description of the contents
 *                    of the 'stress_Plastic_strain' output file.
 *
 *     Parameters:
 *         IN:     baseFileName  Base name of the output file for which
 *                               a description file is to be created.  The
 *                               description file name will be the base
 *                               name with the DESCRIPTION_FILE_SUFFIX
 *                               appended.
 *
 ****************************************************************************/
void WriteStressPlasticStrainFileDesc(Home_t *home, char *baseFileName)
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
 *  The stress/plastic strain output file is only written with
 *  constant strain rate (loadType=1) or cyclic loading (loadType=4)
 */
    if ((home->param->loadType != 1) && (home->param->loadType != 4)) {
        return;
    }

    sprintf(fileName, "%s.%s", baseFileName, DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write output description file %s\n",
               fileName);
        return;
    }

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");
    fprintf(fpDesc, "1         Plastic strain\n");
    fprintf(fpDesc, "2         Stress\n");

    fclose(fpDesc);

    return;
}

#ifdef CALCENERGY
/****************************************************************************
 *
 *      Function:     WriteEnergyPlasticStrainFileDesc
 *      Description:  Create a file containing a description of the contents
 *                    of the 'energy_Plastic_strain' output file.
 *
 *     Parameters:
 *         IN:     baseFileName  Base name of the output file for which
 *                               a description file is to be created.  The
 *                               description file name will be the base
 *                               name with the DESCRIPTION_FILE_SUFFIX
 *                               appended.
 *
 ****************************************************************************/
void WriteEnergyPlasticStrainFileDesc(Home_t *home, char *baseFileName)
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
 *  The stress/plastic strain output file is only written with
 *  constant strain rate (loadType=1) or cyclic loading (loadType=4)
 */
    if ((home->param->loadType != 1) && (home->param->loadType != 4)) {
        return;
    }

    sprintf(fileName, "%s.%s", baseFileName, DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write output description file %s\n",
               fileName);
        return;
    }

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");
    fprintf(fpDesc, "1         Plastic strain\n");
    fprintf(fpDesc, "2         Energy\n");

    fclose(fpDesc);

    return;
}
#endif

//-------------------------------------------------------------------------------------

void WriteStressPlasticStrainFileGnu(Home_t *home, char *baseFileName)
{
    static int gnuFileWritten = 0;

/*
 *  Only one task (domain zero) writes this file, and only writes
 *  it one timer per simulation.  If the associated description
 *  file already exists, the contents are overwritten
 */
    if ((home->myDomain != 0) || (gnuFileWritten == 1)) 
       return;

    gnuFileWritten = 1;

/*
 *  The stress/plastic strain output file is only written with
 *  constant strain rate (loadType=1) or cyclic loading (loadType=4)
 */
    if ((home->param->loadType != 1) && (home->param->loadType != 4))
        return;

    char  fname[MAX_STRING_LEN];
    sprintf(fname, "%s.%s", baseFileName, GNU_FILE_SUFFIX);

    FILE *fp = fopen(fname, "w");

    if (fp)
    {
        fprintf(fp, "#set term x11\n" );
        fprintf(fp, "set term png size 640,480\n" );
        fprintf(fp, "set output 'stress_Plastic_strain.png'\n" );
        fprintf(fp, "\n" );
        fprintf(fp, "set title  'ParaDiS Stress/Strain'\n" );
        fprintf(fp, "set xlabel 'Plastic Strain'  \n" );
        fprintf(fp, "set ylabel 'Stress'  \n" );
        fprintf(fp, "unset key\n" );
        fprintf(fp, "\n" );
        fprintf(fp, "set format x '%%0.1e'\n" );
        fprintf(fp, "set format y '%%0.1e'\n" );
        fprintf(fp, "\n" );
        fprintf(fp, "plot 'stress_Plastic_strain' u 1:2 w lp lc rgb '#ad0000' lt 1 lw 1 pt 7 ps 0.7\n" );
        fprintf(fp, "\n" );

        fclose(fp);
    }
    else
        printf("WARNING: Unable to write output gnu command file %s\n", fname);

    return;
}


/****************************************************************************
 *
 *      Function:     WriteStressTotalStrainFileDesc
 *      Description:  Create a file containing a description of the contents
 *                    of the 'stress_Total_strain' output file.
 *
 *     Parameters:
 *         IN:     baseFileName  Base name of the output file for which
 *                               a description file is to be created.  The
 *                               description file name will be the base
 *                               name with the DESCRIPTION_FILE_SUFFIX
 *                               appended.
 *
 ****************************************************************************/
void WriteStressTotalStrainFileDesc(Home_t *home, char *baseFileName)
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
 *  The stress/total strain output file is only written with
 *  constant strain rate (loadType=1) or cyclic loading (loadType=4)
 */
    if ((home->param->loadType != 1) && (home->param->loadType != 4)) {
        return;
    }

    sprintf(fileName, "%s.%s", baseFileName, DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write output description file %s\n", fileName);
        return;
    }

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");

    if (home->param->loadType == 1) {

        // Constant strain rate
        fprintf(fpDesc, "1         Strain\n");
        fprintf(fpDesc, "2         Stress\n");
        fprintf(fpDesc, "3         Elapsed simulation time\n");
        fprintf(fpDesc, "4         Curent simulation cycle number\n");

    } else if (home->param->loadType == 4) {
        // Cyclic loading
        fprintf(fpDesc, "1         Net accumulated strain\n");
        fprintf(fpDesc, "2         Stress\n");
        fprintf(fpDesc, "3         Elapsed simulation time\n");
        fprintf(fpDesc, "4         Number of loading cycles\n");
    }

    fclose(fpDesc);

    return;
}


/****************************************************************************
 *
 *      Function:     WriteTotalStressFileDesc
 *      Description:  Create a file containing a description of the contents
 *                    of the 'total_stress' output file.
 *
 *     Parameters:
 *         IN:     baseFileName  Base name of the output file for which
 *                               a description file is to be created.  The
 *                               description file name will be the base
 *                               name with the DESCRIPTION_FILE_SUFFIX
 *                               appended.
 *
 ****************************************************************************/
void WriteTotalStressFileDesc(Home_t *home, char *baseFileName)
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
 *  The total stress output file is only written with
 *  constant strain rate (loadType=1) or cyclic loading (loadType=4)
 */
    if ((home->param->loadType != 1) && (home->param->loadType != 4)) {
        return;
    }

    sprintf(fileName, "%s.%s", baseFileName, DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write output description file %s\n",
               fileName);
        return;
    }

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");

    fprintf(fpDesc, "1         Current simulation cycle number\n");
    fprintf(fpDesc, "2         Current simulation time\n");
    fprintf(fpDesc, "3         Applied stress: sigma11\n");
    fprintf(fpDesc, "4         Applied stress: sigma22\n");
    fprintf(fpDesc, "5         Applied stress: sigma33\n");
    fprintf(fpDesc, "6         Applied stress: sigma23\n");
    fprintf(fpDesc, "7         Applied stress: sigma31\n");
    fprintf(fpDesc, "8         Applied stress: sigma12\n");

    fclose(fpDesc);

    return;
}

#ifdef CALCENERGY
/****************************************************************************
 *
 *      Function:     WriteEnergyStrainFileDesc
 *      Description:  Create a file containing a description of the contents
 *                    of the 'total_energy' output file.
 *
 *     Parameters:
 *         IN:     baseFileName  Base name of the output file for which
 *                               a description file is to be created.  The
 *                               description file name will be the base
 *                               name with the DESCRIPTION_FILE_SUFFIX
 *                               appended.
 *
 ****************************************************************************/
void WriteEnergyStrainFileDesc(Home_t *home, char *baseFileName)
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
 *  The total stress output file is only written with
 *  constant strain rate (loadType=1) or cyclic loading (loadType=4)
 */
    if ((home->param->loadType != 1) && (home->param->loadType != 4)) {
        return;
    }

    sprintf(fileName, "%s.%s", baseFileName, DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write output description file %s\n",
               fileName);
        return;
    }

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");

    fprintf(fpDesc, "1         Total strain\n");
    fprintf(fpDesc, "2         Total energy\n");
    fprintf(fpDesc, "3         Current simulation time\n");
    fprintf(fpDesc, "4         Current cycle\n");

    fclose(fpDesc);

    return;
}
#endif

/****************************************************************************
 *
 *      Function:     WriteAllepsFileDesc
 *      Description:  Create a file containing a description of the contents
 *                    of the 'alleps' output file.
 *
 *     Parameters:
 *         IN:     baseFileName  Base name of the output file for which
 *                               a description file is to be created.  The
 *                               description file name will be the base
 *                               name with the DESCRIPTION_FILE_SUFFIX
 *                               appended.
 *
 ****************************************************************************/
void WriteAllepsFileDesc(Home_t *home, char *baseFileName)
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

    sprintf(fileName, "%s.%s", baseFileName, DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write output description file %s\n",
               fileName);
        return;
    }

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");

    fprintf(fpDesc, "1         Current simulation cycle number\n");
    fprintf(fpDesc, "2         Elapsed simulation time\n");
    fprintf(fpDesc, "3-8       Plastic strain tensor elements in sequence\n");
    fprintf(fpDesc, "          [0][0], [1][1], [2][2], [1][2], [0][1], [0][2]\n");
    fprintf(fpDesc, "9         Dislocation density\n");

    fclose(fpDesc);

    return;
}


/****************************************************************************
 *
 *      Function:     WriteEpsdotFileDesc
 *      Description:  Create a file containing a description of the contents
 *                    of the 'epsdot' output file.
 *
 *     Parameters:
 *         IN:     baseFileName  Base name of the output file for which
 *                               a description file is to be created.  The
 *                               description file name will be the base
 *                               name with the DESCRIPTION_FILE_SUFFIX
 *                               appended.
 *
 ****************************************************************************/
void WriteEpsdotFileDesc(Home_t *home, char *baseFileName)
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

    sprintf(fileName, "%s.%s", baseFileName, DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write output description file %s\n",
               fileName);
        return;
    }

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");

    fprintf(fpDesc, "1         Elapsed simulation time\n");
    fprintf(fpDesc, "2         Plastic strain rate\n");

    fclose(fpDesc);

    return;
}


/****************************************************************************
 *
 *      Function:     WriteDensityDeltaFileDesc
 *      Description:  Create a file containing a description of the contents
 *                    of the 'density_delta' output file.
 *
 *     Parameters:
 *         IN:     baseFileName  Base name of the output file for which
 *                               a description file is to be created.  The
 *                               description file name will be the base
 *                               name with the DESCRIPTION_FILE_SUFFIX
 *                               appended.
 *
 ****************************************************************************/
void WriteDensityDeltaFileDesc(Home_t *home, char *baseFileName)
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
 *  The density delta file is only written when one of the BCC
 *  mobilities is used.
 */
    if (home->param->materialType != MAT_TYPE_BCC) {
        return;
    }

    sprintf(fileName, "%s.%s", baseFileName, DESCRIPTION_FILE_SUFFIX);

    if ((fpDesc = fopen(fileName, "w")) == (FILE *)NULL) {
        printf("WARNING: Unable to write output description file %s\n",
               fileName);
        return;
    }

    fprintf(fpDesc, "Column    Description\n");
    fprintf(fpDesc, "------    -------------------------------------------\n");

    fprintf(fpDesc, "1         Strain\n");
    fprintf(fpDesc, "2         Total incremental gain (sum of columns\n");
    fprintf(fpDesc, "          4 thru 10\n");
    fprintf(fpDesc, "3         Total incremental loss (sum of columns\n");
    fprintf(fpDesc, "          11 thru 17\n");
    fprintf(fpDesc, "4         Density gain: Burgers vectors [-1  1  1]\n");
    fprintf(fpDesc, "                                        [ 1 -1 -1]\n");
    fprintf(fpDesc, "5         Density gain: Burgers vectors [ 1 -1  1]\n");
    fprintf(fpDesc, "                                        [-1  1 -1]\n");
    fprintf(fpDesc, "6         Density gain: Burgers vectors [ 1  1 -1]\n");
    fprintf(fpDesc, "                                        [-1 -1  1]\n");
    fprintf(fpDesc, "7         Density gain: Burgers vectors [ 1  1  1]\n");
    fprintf(fpDesc, "                                        [-1 -1 -1]\n");
    fprintf(fpDesc, "8         Density gain: Burgers vectors [ 1  0  0]\n");
    fprintf(fpDesc, "                                        [-1  0  0]\n");
    fprintf(fpDesc, "9         Density gain: Burgers vectors [ 0  1  0]\n");
    fprintf(fpDesc, "                                        [ 0 -1  0]\n");
    fprintf(fpDesc, "10        Density gain: Burgers vectors [ 0  0  1]\n");
    fprintf(fpDesc, "                                        [ 0  0 -1]\n");
    fprintf(fpDesc, "11        Density loss: Burgers vectors [-1  1  1]\n");
    fprintf(fpDesc, "                                        [ 1 -1 -1]\n");
    fprintf(fpDesc, "12        Density loss: Burgers vectors [ 1 -1  1]\n");
    fprintf(fpDesc, "                                        [-1  1 -1]\n");
    fprintf(fpDesc, "13        Density loss: Burgers vectors [ 1  1 -1]\n");
    fprintf(fpDesc, "                                        [-1 -1  1]\n");
    fprintf(fpDesc, "14        Density loss: Burgers vectors [ 1  1  1]\n");
    fprintf(fpDesc, "                                        [-1 -1 -1]\n");
    fprintf(fpDesc, "15        Density loss: Burgers vectors [ 1  0  0]\n");
    fprintf(fpDesc, "                                        [-1  0  0]\n");
    fprintf(fpDesc, "16        Density loss: Burgers vectors [ 0  1  0]\n");
    fprintf(fpDesc, "                                        [ 0 -1  0]\n");
    fprintf(fpDesc, "17        Density loss: Burgers vectors [ 0  0  1]\n");
    fprintf(fpDesc, "                                        [ 0  0 -1]\n");

    fclose(fpDesc);

    return;
}
