/****************************************************************************
 *
 *      Module:         CTableGen.c
 *      Description:    Contains the main routine for controlling
 *                      generation of a PBC image correction table
 *                      used in conjunction with the fast multipole
 *                      code.
 *
 *      Included functions:
 *          main()
 *          GetInArgs()
 *          InitValues()
 *          PrintHelp()
 *          Usage()
 *
 *      Usage:  ctablegen {-nu <poissonratio> -mu <shearmodulus> ||       \
 *                  -ecfile <basefile> -pressure <pressure> -temp <temp>} \
 *                  [-aniso]                                              \
 *                  [-cubesize <boxsize>]                                 \
 *                  [-derivfile <derivatives_file>]                       \
 *                  [-help]                                               \
 *                  [-levels <numlayers>]                                 \
 *                  [-outfile <outputfilename>]                           \
 *                  [-mporder <multipole_expansion_order>]                \
 *                  [-mattype {bcc|fcc|hcp|rhombohedral}]                 \
 *                  [-torder <taylor_expansion_order>]                    \
 *                  [-pbc <pbcVal>]
 *
 *      Overview:
 *
 *      This utility will create a file containing a table of
 *      values needed to do correct the primary image multipole
 *      expansion to allow for multiple PBC images of the problem
 *      space.
 *
 *      The contents of the table are determined by the user
 *      supplied command line arguments such as shear modulus,
 *      multipole expansion order, taylor expansion order, etc.
 *      
 *      All options may be abbreviated to the shortest non-ambiguous
 *      abbreviation of the option.  Most options take on default values
 *      if not specified on the command line.  
 *      
 ***************************************************************************/
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>

#include "mpi_portability.h"

#include "Home.h"

#ifdef TAYLORFMM
#include "FM.h"


/*
 *      Define an integer identifier to be associated with each
 *      posible command line argument.  To be used to index the
 *      option-specific data in the optList array below.  OPT_MAX
 *      must be the last value in the list
 */
enum {
    OPT_ANISO_ENABLED = 0,
    OPT_CUBESIZE,
    OPT_ECFILE,
    OPT_DERIV_FILE,
    OPT_HELP,
    OPT_LEVELS,
    OPT_MAT_TYPE,
    OPT_MPORDER,
    OPT_MU,
    OPT_NU,
    OPT_OUTFILE,
    OPT_PBC,
    OPT_PRESSURE,
    OPT_TEMP,
    OPT_TORDER,
    OPT_MAX       /* This MUST be the last element in the enumerated list */
};


/*
 *      Define a structure to hold a command line option's id (type),
 *      name, the shortest possible unique abbreviation of the option
 *      name, and a flag indicating if the option is paired with
 *      a value or not.
 */
typedef struct {
         int     optType;
   const char   *optName;
         int     optMinAbbrev;
         int     optPaired;
} Option_t;


/*
 *      Define and initialize an array of structures containing
 *      all the possible command line arguments, and some info
 *      determining how it is treated.
 *
 *      option            option       #characters     1 unless option
 *      type              name          in unique      has no associated
 *                                     abbreviation    value
 */
Option_t        optList[OPT_MAX] = {
        {OPT_ANISO_ENABLED, (const char *) "aniso",     1,              0},
        {OPT_CUBESIZE,      (const char *) "cubesize",  1,              1},
        {OPT_DERIV_FILE,    (const char *) "derivfile", 1,              1},
        {OPT_ECFILE,        (const char *) "ecfile",    1,              1},
        {OPT_HELP,          (const char *) "help",      1,              0},
        {OPT_LEVELS,        (const char *) "levels",    1,              1},
        {OPT_MAT_TYPE,      (const char *) "mattype",   2,              1},
        {OPT_MPORDER,       (const char *) "mporder",   2,              1},
        {OPT_MU,            (const char *) "mu",        2,              1},
        {OPT_NU,            (const char *) "nu",        1,              1},
        {OPT_OUTFILE,       (const char *) "outfile",   1,              1},
        {OPT_PBC,           (const char *) "pbc",       2,              1},
        {OPT_PRESSURE,      (const char *) "pressure",  2,              1},
        {OPT_TEMP,          (const char *) "temp",      2,              1},
        {OPT_TORDER,        (const char *) "torder",    2,              1}
};


/*---------------------------------------------------------------------------
 *
 *      Function:       Usage
 *      Description:    Print out a brief message indicating the possible
 *                      command line options and terminated.
 *
 *      Arguments:
 *          program   name of the program being executed
 *
 *-------------------------------------------------------------------------*/
static void Usage(char *program)
{

        printf("Usage:  %10s {-nu <poissonratio> -mu <shearmodulus> || \\\n",
               program);
        printf("              -ecfile <basefile> -pressure <pressure> -temp <temp>}\\\n");
        printf("                  [-aniso]                             \\\n");
        printf("                  [-derivfile <derivatives_file>]      \\\n");
        printf("                  [-cubesize <boxsize>]                \\\n");
        printf("                  [-help]                              \\\n");
        printf("                  [-levels <numlayers>]                \\\n");
        printf("                  [-outfile <outputfilename>]          \\\n");
        printf("                  [-mporder <multipole_exp_order>]     \\\n");
        printf("                  [-mattype {bcc|fcc|hcp|rhombohedral}]\\\n");
        printf("                  [-pbc <pbc_val>]                     \\\n");
        printf("                  [-torder <taylor_exp_order>]\n");

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       PrintHelp
 *      Description:    Print out a detailed description of the available
 *                      program options, what they represent, how they
 *                      relate to each other, which are interdependent, or
 *                      mutually exclusive, default values, etc.
 *
 *      Arguments:
 *              program         name of the program being executed
 *
 *-------------------------------------------------------------------------*/
static void PrintHelp(char *program)
{
    Usage(program);

    printf("    Options may be abbreviated to the shortest non-ambiguous\n");
    printf("    abbreviation of the option.\n\n");
    printf("Options:\n\n");
    printf("  -aniso    If the code was compiled with support for \n");
    printf("            anisotropy, this flag will enable the creation \n");
    printf("            of a correction table for use with aniostropic fmm.\n");
    printf("            This flag is disabled by default.\n\n");
    printf("  -cubesize Define the length (units of b) of a single side\n");
    printf("            of the cubic problem space for which this \n");
    printf("            correction table is being generated.\n\n");
    printf("  -derivfile Specifies the name of a file containing the\n");
    printf("            Greens function derivatives used with\n");
    printf("            anisotropic elasticity only.  Required with\n");
    printf("            when anisotropy is enabled.  Default is:\n");
    printf("            'fm-aniso-deriv-table.dat'\n\n");
    printf("  -ecfile   Base file name of material specific elastic\n");
    printf("            constants files.\n\n");
    printf("  -help     Prints this help information.\n\n");
    printf("  -mu       Shear modulus.\n\n");
    printf("  -nu       Poisson ratio.\n\n");
    printf("  -levels   Number of levels.  Used to determine the number\n");
    printf("            periodic images of the primary problem space\n");
    printf("            in each dimension. Number of periodic images\n");
    printf("            per dimension = 2^levels.\n\n");
    printf("  -mattype  Type of material structure.  Determines if the\n");
    printf("            full correction table is created or if redundant\n");
    printf("            coefficents due to crystal symmetries are factored\n");
    printf("            out of the table.  Valid types are: 'bcc', 'fcc',\n");
    printf("            'hcp', 'rhombohedral'.  Default = 'bcc'.\n\n");
    printf("  -mporder  Order of the multipole expansion for which this.\n");
    printf("            correction table is being built.  Must be integer\n");
    printf("            greater than or equal to zero.  Default = 2.\n\n");
    printf("  -torder   Order of the taylor expansion for which this.\n");
    printf("            correction table is being built.  Must be integer\n");
    printf("            greater than or equal to zero.  Default = 5.\n\n");
    printf("  -pbc      Integer value specifying dimensions in which \n");
    printf("            periodic boundaries are enabled.  Default is \n");
    printf("            all dimensions.  Values can be OR'ed together to\n");
    printf("            yield periodic boundaries in any combination of\n");
    printf("            dimensions.\n");
    printf("              1 == periodic in X dimensions.\n");
    printf("              2 == periodic in Y dimensions.\n");
    printf("              4 == periodic in Z dimensions.\n\n");
    printf("  -pressure Pressure in units of Pa.\n\n");
    printf("  -temp     Temperature in degress K.\n\n");
    printf("  -outfile  Name of the file into which to write the\n");
    printf("            image correction table data. \n\n");

    exit(0);
}


/*---------------------------------------------------------------------------
 *
 *      Function:  GetInArgs
 *      Description: Parse and process and user supplied command line
 *                   options.  Set appropriate default values (if
 *                   permitted) for any values not provided by the
 *                   user.
 *
 *      Arguments:
 *          argc       count of command line args
 *          argv       command line args
 *          param      Pointer to structure to be populated with
 *                     user supplied or default values.
 *          numLevels  Number of levels used in hierarchy when
 *                     creating the image correction table.
 *
 *-------------------------------------------------------------------------*/
static void GetInArgs(int argc, char *argv[], Param_t *param, int *numLevels,
                      int thisTask, int pbc[3], char *derivFile)
{
        int     i, j, pbcVal;
        real8   maxSide, minSide;
        char    *argName=0;
        char    *argValue=0;

        for (i = 1; i < argc; i++) {
/*
 *              If the option doesn't begin with a '-' something
 *              is wrong, so notify the user and terminate.
 */
        if (argv[i][0] != '-') {
            Usage(argv[0]);
            exit(1);
        }

        argName = &argv[i][1];

/*
 *      Scan the array of valid options for the user supplied
 *      option name.  (This may be any unique abbreviation of
 *      one of any of the options.  If we don't find the option
 *      notify the user and terminate.
 */
        for (j = 0; j < OPT_MAX; j++) {
            if (!strncmp(argName, optList[j].optName,
                         optList[j].optMinAbbrev)) {
                break;
            }
        }

        if (j == OPT_MAX) {
            Usage(argv[0]);
            exit(1);
        }

/*
 *      Verify that there is an associated value if the specified
 *      option is supposed to be paired with a value.
 */
        if (optList[j].optPaired) {
            if (i+1 >= argc) {
                Usage(argv[0]);
                exit(1);
            } else {
                argValue = argv[++i];
            }
        }

/*
 *      Do any option-specific processing...
 */
        switch (optList[j].optType)  {
            case OPT_ANISO_ENABLED:
                param->fmEnableAnisotropy = 1;
                break;
            case OPT_CUBESIZE:
                sscanf(argValue, "%le", &maxSide);
                maxSide  = maxSide / 2.0;
                minSide  = -maxSide;

                param->minSideX=minSide;
                param->minSideY=minSide;
                param->minSideZ=minSide;

                param->maxSideX=maxSide;
                param->maxSideY=maxSide;
                param->maxSideZ=maxSide;

                break;
            case OPT_DERIV_FILE:
                strcpy(derivFile, argValue);
                break;
            case OPT_ECFILE:
                strcpy(param->ecFile, argValue);
                break;
            case OPT_OUTFILE:
                strcpy(param->fmCorrectionTbl, argValue);
                break;
            case OPT_LEVELS:
                *numLevels = atoi(argValue);
                break;
            case OPT_MU:
                sscanf(argValue, "%le", &param->shearModulus);
                break;
            case OPT_NU:
                sscanf(argValue, "%le", &param->pois);
                break;
            case OPT_PBC:
                pbcVal = atoi(argValue);
                pbc[0] = ((pbcVal & 0x01) > 0);
                pbc[1] = ((pbcVal & 0x02) > 0);
                pbc[2] = ((pbcVal & 0x04) > 0);
                break;
            case OPT_PRESSURE:
                sscanf(argValue, "%le", &param->pressure);
                break;
            case OPT_TEMP:
                sscanf(argValue, "%le", &param->TempK);
                break;
            case OPT_MAT_TYPE:
                if (StrEquiv(argValue, "bcc")) {
                    param->materialType = MAT_TYPE_BCC;
                } else if (StrEquiv(argValue, "fcc")) {
                    param->materialType = MAT_TYPE_FCC;
                } else if (StrEquiv(argValue, "hcp")) {
                    param->materialType = MAT_TYPE_HCP;
                } else if (StrEquiv(argValue, "rhombohedral")) {
                    param->materialType = MAT_TYPE_RHOMBOHEDRAL_VA;
                } else {
                    Fatal("Invalid option '-mattype %s'", argValue);
                }
                break;
            case OPT_MPORDER:
                param->fmMPOrder = atoi(argValue);
                break;
            case OPT_TORDER:
                param->fmExpansionOrder = atoi(argValue);
                break;
            case OPT_HELP:
                PrintHelp(argv[0]);
                exit(0);
                break;
            }
        }

/*
 *      User must specify poisson ratio and shear modulus, OR
 *      temperature, pressure and an elastic constants file.
 */
        if ((param->pois < 0.0) || (param->shearModulus < 0.0)) {
/*
 *          If user did not provide either set of required parameters
 *          it's an error.
 */
            if ((param->TempK < 0.0) || (param->pressure < 0.0) ||
                (param->ecFile[0] == 0)) {
                Usage(argv[0]);
                exit(1);
            }

/*
 *          User provided temp and pressure, so calculate the
 *          shear modulus and poisson ratio
 */
            if (GetElasticConstants(param) == 0) {
                if (thisTask == 0) {
                    printf("Error calculating MU and NU for specified"
                           "temperature and pressure\n");
                }
                exit(1);
            }
 
            if (thisTask == 0) {
                printf(" ** For temp = %6.2fK, pressure = %ePa:\n",
                       param->TempK, param->pressure);
                printf(" **   Setting poisson ratio to %e\n", param->pois);
                printf(" **   Setting shear modulus to %e\n",
                       param->shearModulus);
                printf(" **   Setting burgers vector magnitude to %e\n",
                       param->burgMag);
            }
            
            
        } else {
/*
 *          If user provided both sets of parameters, just use
 *          the provided shear modulus and poisson ratio but
 *          warn him that is being done.
 */
            if ((param->TempK > 0.0) || (param->pressure > 0.0) ||
                (param->ecFile[0] != 0)) {
                if (thisTask == 0) {
                    printf("Warning: MU and NU specified; temperature, pressure\n"
                           "and elastic constants file will be ignored.\n");
                }
            }
        }

/*
 *      Set remaining default values if permitted
 */
        if (*numLevels < 3) {
            *numLevels = 3;
        }

        if (param->maxSideX <= 0.0) {
            real8  maxSide, minSide;

            maxSide  = 35000.0 / 2.0;
            minSide  = -maxSide;

            param->minSideX=minSide;
            param->minSideY=minSide;
            param->minSideZ=minSide;

            param->maxSideX=maxSide;
            param->maxSideY=maxSide;
            param->maxSideZ=maxSide;
        }

        if (param->fmMPOrder < 0) {
            param->fmMPOrder = 2;
            param->fmExpansionOrder = 5;
        }

        if (param->fmExpansionOrder < 0) {
            param->fmExpansionOrder = (param->fmMPOrder * 2) + 1;
        }

/*
 *      Default correction table name is dependent on whether
 *      anisotropy is enabled or disabled
 */
        if (param->fmCorrectionTbl[0] == 0) {
#ifdef ANISOTROPIC
            sprintf(param->fmCorrectionTbl, "fm-ctab-aniso.m%d.t%d.l%d.dat",
                    param->fmMPOrder, param->fmExpansionOrder, *numLevels);
#else
            sprintf(param->fmCorrectionTbl, "fm-ctab.m%d.t%d.l%d.dat",
                    param->fmMPOrder, param->fmExpansionOrder, *numLevels);
#endif
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    InitValues
 *      Description: Set some initial invalid values for the various
 *                   parameters controlling the correction table 
 *                   generation.  After input arguments have been
 *                   processed, these configuration values can be
 *                   examined and set to appropriate defaults if
 *                   no values were provided by the caller.
 *
 *-------------------------------------------------------------------------*/
static void InitValues(Param_t *param, int *numLayers, int pbc[3],
                       char *derivFile)
{
        param->fmMPOrder     = -1;
        param->fmExpansionOrder = -1;
        param->shearModulus  = -1.0;
        param->pois          = -1.0;
        param->pressure      = -1.0;
        param->TempK         = -1.0;
        param->materialType  = MAT_TYPE_BCC;

        param->minSideX = -1.0;
        param->minSideY = -1.0;
        param->minSideZ = -1.0;

        param->maxSideX = -1.0;
        param->maxSideY = -1.0;
        param->maxSideZ = -1.0;

        memset(param->fmCorrectionTbl, 0, sizeof(param->fmCorrectionTbl));

        *numLayers = 10;

/*
 *      Assume periodic boundaries in all dimensions as default
 */
        pbc[0] = 1;
        pbc[1] = 1;
        pbc[2] = 1;

/*
 *      Set the default file name for the Greens function
 *      derivatives files.  Only needed with anisotropic elasticity.
 */
        strcpy(derivFile, "fm-aniso-deriv-table.dat");

        return;
}


int main(int argc, char *argv[])
{
        int      numLevels, numTasks, thisTask;
        int      pbc[3];
        char     derivativesFile[1024];
        Home_t  *home;
        Param_t  param;
#ifdef ANISOTROPIC
        int      tblMaxOrder;
        real8    xScale, yScale, zScale;
        FILE    *fpTable;
        FMAnisoTable_t *anisoTable;
#endif

        home = (Home_t *)calloc(1, sizeof(Home_t));
        memset((void *)&param, 0, sizeof(Param_t));

#ifdef ANISOTROPIC
        home->fmAnisoTable =
               (FMAnisoTable_t *)calloc(1, sizeof(FMAnisoTable_t));
        anisoTable = home->fmAnisoTable;
#endif  /* ifdef ANISOTROPIC */

#ifdef PARALLEL
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &thisTask);
        MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
#else
        thisTask = 0;
        numTasks = 1;
#endif

/*
 *      Do some intialization stuff
 */
        home->myDomain   = thisTask;
        home->numDomains = numTasks;

        InitValues(&param, &numLevels, pbc, derivativesFile);

        GetInArgs(argc, argv, &param, &numLevels, thisTask, pbc,
                  derivativesFile);

        home->param = &param;

        param.minCoordinates[0] = param.minSideX;
        param.minCoordinates[1] = param.minSideY;
        param.minCoordinates[2] = param.minSideZ;

        param.maxCoordinates[0] = param.maxSideX;
        param.maxCoordinates[1] = param.maxSideY;
        param.maxCoordinates[2] = param.maxSideZ;


#ifdef ANISOTROPIC
/*
 *      If the user requested anisotropic fmm, read the metadata from
 *      the derivatives file and save some of the values for later use.
 */
        if (param.fmEnableAnisotropy) {
            real8    cellSize[3];

            strcpy(param.fmAnisoDerivTbl, derivativesFile);

            fpTable = fopen(param.fmAnisoDerivTbl, "r");

            if (fpTable == (FILE *)NULL) {
                Fatal("FMAnisotropicInit: Open error %d on anistropy "
                      "derivatives file %s", errno, param.fmAnisoDerivTbl);
            }

            ReadDerivTableMetadata(fpTable, param.fmAnisoDerivTbl,
                                   param.elasticConstantMatrix,
                                   &xScale, &yScale, &zScale, &tblMaxOrder);

            anisoTable->xScale = xScale;
            anisoTable->yScale = yScale;
            anisoTable->zScale = zScale;

            anisoTable->scaleVec[0] = xScale;
            anisoTable->scaleVec[1] = yScale;
            anisoTable->scaleVec[2] = zScale;

            anisoTable->derivTblMaxOrder = tblMaxOrder;

/*
 *         Convert the 6x6 elastic constants matrix into the 3x3x3x3 form
 *         and initialize the table of Greens function derivatives
 */
           SetElasticConstantMatrix4D(param.elasticConstantMatrix,
                                      home->anisoVars.elasticConstantMatrix4D);

/*
 *          Verify the correctness of the {x,y,z} scale...
 */
            cellSize[0] = (param.maxCoordinates[0] - param.minCoordinates[0]) /
                          param.nXcells;
            cellSize[1] = (param.maxCoordinates[1] - param.minCoordinates[1]) /
                          param.nYcells;
            cellSize[2] = (param.maxCoordinates[2] - param.minCoordinates[2]) /
                          param.nZcells;

            if ((fabs(xScale*cellSize[1]/(yScale*cellSize[0]) - 1.0) > 1e-6) ||
                (fabs(yScale*cellSize[2]/(zScale*cellSize[1]) - 1.0) > 1e-6)) {
                Fatal("FMAnisotropicInit: Coordinate scalings in Greens "
                      "function file %s\nare not commensurate with the "
                      "current simulation cell size", param.fmAnisoDerivTbl);
            }

            ReadDerivTable(home, fpTable);

            fclose(fpTable);
        }

#endif  /* ifdef ANISOTROPIC */

#ifdef PARALLEL
/*
 *      Distribute param structure to all processes
 */
        MPI_Bcast((char *)&param, sizeof(Param_t), MPI_CHAR, 0,
                  MPI_COMM_WORLD);

        if (thisTask == 0)
#endif
        {
            printf(" ** Using multipole expansion order %d\n"
                   " ** Using taylor expansion order    %d\n"
                   " ** Using number of levels          %d\n",
                   param.fmMPOrder, param.fmExpansionOrder, numLevels);
        }

        Taylor_CreateCorrectionTable(home, numLevels, pbc, numTasks, thisTask);

#ifdef PARALLEL
        MPI_Finalize();
#endif

        exit(0);
        return(0);
}

#else

int main(int argc, char *argv[])
{
        printf("\nTaylor FMM must be enabled in  makefile.setup and\n");
        printf("this utility recompiled before it may be used!\n\n");

        exit(0);
        return(0);
}

#endif

