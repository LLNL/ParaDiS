/****************************************************************************
 *
 *      Module:         StanfordFMsMTableGen.c
 *      Description:    Contains the main routines for controlling
 *                      generation of a PBC image correction table
 *                      used in conjunction with the BB and UNIFORM Fast 
 *                      Multipole code.
 *
 *      Included functions:
 *          main()
 *          GetInArgs()
 *          InitValues()
 *          PrintHelp()
 *          Usage()
 *
 *      Usage:  StanfordFMMsTableGen [-cubesize <boxsize>] [-C11 <c11>]      \
 *                  [-C12 <c12>] [-C13 <c13>] [-C33 <c33>] [-C44 <c44>]      \
 *                  [-cmatfile <filename>] [-help]  [-levels <numlayers>]    \
 *                  [-nphi <nphi>] [-ntheta <ntheta>] [-qmax <qMax>]         \
 *                  [-corder <chebyshev_order>] [-outfile <outputfilename>] 
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
#include "Home.h"
#include "Mobility.h"
#include "FM.h"

#if defined (BBFMM) || defined (UNIFORMFMM)

#if defined ANISOTROPIC
#include "Anisotropic.h"
#endif

#define MAX_LINE_LEN 512

typedef struct {
        int   FMMOrder;
        int   numLevels;
        char  cMatrixFile[512];
} InputArgs_t;

/*
 *      Define an integer identifier to be associated with each
 *      posible command line argument.  To be used to index the
 *      option-specific data in the optList array below.  OPT_MAX
 *      must be the last value in the list and 1 greater than
 */

#ifdef ANISOTROPIC
enum {
    OPT_C11 = 0,
    OPT_C12,
    OPT_C13,
    OPT_C33,
    OPT_C44,
    OPT_nphi,
    OPT_ntheta,
    OPT_qmax,
    OPT_CMATRIX_FILE,
    OPT_CORDER,
    OPT_CUBESIZE,
    OPT_HELP,
    OPT_LEVELS,
    OPT_MAT_TYPE,
    OPT_OUTFILE,
    OPT_MAX       /* This MUST be the last element in the enumerated list */
};
#else
enum {
    OPT_MU = 0,
    OPT_NU,
    OPT_CORDER,
    OPT_CUBESIZE,
    OPT_HELP,
    OPT_LEVELS,
    OPT_MAT_TYPE,
    OPT_OUTFILE,
    OPT_MAX       /* This MUST be the last element in the enumerated list */
};
#endif


/*
 *      Define a structure to hold a command line option's id (type),
 *      name, the shortest possible unique abbreviation of the option
 *      name, and a flag indicating if the option is paired with
 *      a value or not.
 */
typedef struct {
        int     optType;
  const char    *optName;
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
#ifdef ANISOTROPIC
Option_t        optList[OPT_MAX] = {
        {OPT_C11,        "c11",          3,             1},
        {OPT_C12,        "c12",          3,             1},
        {OPT_C13,        "c13",          3,             1},
        {OPT_C33,        "c33",          2,             1},
        {OPT_C44,        "c44",          2,             1},
        {OPT_nphi,       "nphi",         4,             1},
        {OPT_ntheta,     "ntheta",       5,             1},
        {OPT_qmax,       "qmax",         4,             1},
        {OPT_CMATRIX_FILE, "cmatfile",   2,             1},
        {OPT_CORDER,     "corder",       2,             1},
        {OPT_CUBESIZE,   "cubesize",     2,             1},
        {OPT_HELP,       "help",         1,             0},
        {OPT_LEVELS,     "levels",       1,             1},
        {OPT_MAT_TYPE,   "mattype",      1,             1},
        {OPT_OUTFILE,    "outfile",      1,             1} 
};
#else
Option_t        optList[OPT_MAX] = {
        {OPT_MU,        "mu",            2,             1},
        {OPT_NU,        "nu",            2,             1},
        {OPT_CORDER,     "corder",       2,             1},
        {OPT_CUBESIZE,   "cubesize",     2,             1},
        {OPT_HELP,       "help",         1,             0},
        {OPT_LEVELS,     "levels",       1,             1},
        {OPT_MAT_TYPE,   "mattype",      1,             1},
        {OPT_OUTFILE,    "outfile",      1,             1} 
};
#endif


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
#ifdef ANISOTROPIC
        printf("Usage:  %10s \\\n", program);
        printf("                [-cubesize <boxsize>]                \\\n");
        printf("                [-c11 <c11>] [-c12 <c12>] [-c13 <c13>]\\\n");
        printf("                [-c33 <c33>] [-c44 <c44>]            \\\n");
        printf("                [-cmatfile <C_matrix_filename>]      \\\n");
        printf("                [-corder <FMM_order>]          \\\n");
        printf("                [-mattype {bcc|fcc|hcp|rhombohedral} \\\n");
        printf("                [-help] [-levels <numlayers>]        \\\n");
        printf("                [-outfile <outputfilename>]          \\\n");

#else
        printf("Usage:  %10s \\\n", program);
        printf("                [-cubesize <boxsize>]                \\\n");
        printf("                [-mu <mu>] [-nu <nu>] \\\n");
        printf("                [-cmatfile <C_matrix_filename>]      \\\n");
        printf("                [-corder <FMM_order>]          \\\n");
        printf("                [-mattype {bcc|fcc|hcp|rhombohedral} \\\n");
        printf("                [-help] [-levels <numlayers>]        \\\n");
        printf("                [-outfile <outputfilename>]          \\\n");

#endif
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
    printf("  -cubesize Define the length (units of b) of a single side\n");
    printf("            of the cubic problem space for which this \n");
    printf("            correction table is being generated.\n\n");
#ifdef ANISOTROPIC
    printf("  -c11      Elastic constant\n\n");
    printf("  -c12      Elastic constant\n\n");
    printf("  -c13      Elastic constant\n\n");
    printf("  -c33      Elastic constant\n\n");
    printf("  -c44      Elastic constant\n\n");
    printf("  -qmax     Harmonics series cut-off parameter\n");
#else
    printf("  -mu      Elastic constant\n\n");
    printf("  -nu      Elastic constant\n\n");
#endif
    printf("  -cmatfile Name of file containing the 6x6 elastic constant\n");
    printf("            matrix.  Not needed if individual elastic \n");
    printf("            constants are provided: format: format is 6\n");
    printf("            lines, each with 1 row of the matrix. \n\n");
    printf("  -help     Prints this help information.\n\n");
    printf("  -levels   Number of levels.  Used to determine the number\n");
    printf("            periodic images of the primary problem space\n");
    printf("            in each dimension. Number of periodic images\n");
    printf("            per dimension = 2^levels. Integer greater than 3.\n\n");
    printf("  -corder   Order of the FMM expansion with which this.\n");
    printf("            correction table is being built.  Must be integer\n");
    printf("            greater than or equal to zero.  Default = 3.\n\n");
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
static void GetInArgs(int argc, char *argv[], InputArgs_t *inArgs,
                      Param_t *param)
{
        int     i, j, pbcVal;
        real8   maxSide, minSide;
        char    *argName;
        char    *argValue;

        for (i = 1; i < argc; i++) {
/*
 *          If the option doesn't begin with a '-' something
 *          is wrong, so notify the user and terminate.
 */
            if (argv[i][0] != '-') {
                Usage(argv[0]);
                exit(1);
            }

            argName = &argv[i][1];

/*
 *          Scan the array of valid options for the user supplied
 *          option name.  (This may be any unique abbreviation of
 *          one of any of the options.  If we don't find the option
 *          notify the user and terminate.
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
 *          Verify that there is an associated value if the specified
 *          option is supposed to be paired with a value.
 */
            if (optList[j].optPaired && (i+1 >= argc)) {
                Usage(argv[0]);
                exit(1);
            } else argValue = argv[++i];

/*
 *          Do any option-specific processing...
 */

#ifdef ANISOTROPIC
            switch (optList[j].optType)  {
                case OPT_CORDER:
                    inArgs->FMMOrder = atoi(argValue);
                    param->fmExpansionOrder = atoi(argValue);
                    break;
                case OPT_C11:
                    param->C11 = atof(argValue);
                    break;
                case OPT_C12:
                    param->C12 = atof(argValue);
                    break;
                case OPT_C13:
                    param->C13 = atof(argValue);
                    break;
                case OPT_C33:
                    param->C33 = atof(argValue);
                    break;
                case OPT_C44:
                    param->C44 = atof(argValue);
                    break;
                case OPT_nphi:
                    param->anisoNumPhiPoints = atof(argValue);
                    break;
                case OPT_ntheta:
                    param->anisoNumThetaPoints = atof(argValue);
                    break;
                case OPT_qmax:
                    param->anisoHarmonicsNumTermsBase = atof(argValue);
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
                case OPT_CMATRIX_FILE:
                    strcpy(inArgs->cMatrixFile, argValue);
                    break;
                case OPT_LEVELS:
                    inArgs->numLevels = atoi(argValue);
                    break;
                case OPT_MAT_TYPE:
                    if (StrEquiv(argValue, "bcc")) {
                        param->materialType = MAT_TYPE_BCC;
                        strcpy(param->materialTypeName, "BCC");
                    } else if (StrEquiv(argValue, "fcc")) {
                        param->materialType = MAT_TYPE_FCC;
                        strcpy(param->materialTypeName, "FCC");
                    } else if (StrEquiv(argValue, "hcp")) {
                        param->materialType = MAT_TYPE_HCP;
                        strcpy(param->materialTypeName, "HCP");
                    } else if (StrEquiv(argValue, "rhmobohedral")) {
                        param->materialType = MAT_TYPE_RHOMBOHEDRAL_VA;
                        strcpy(param->materialTypeName, "RHOMBOHEDRAL");
                    } else {
                        printf("Error: Unknown material type %s\n", argValue);
                        exit(1);
                    }
                    break;
                case OPT_OUTFILE:
                    strcpy(param->fmCorrectionTbl, argValue);
                    break;
                case OPT_HELP:
                    PrintHelp(argv[0]);
                    exit(0);
                    break;

            }  /* end switch (optList[j].optType) */

#else

            switch (optList[j].optType)  {
                case OPT_CORDER:
                    inArgs->FMMOrder = atoi(argValue);
                    param->fmExpansionOrder = atoi(argValue);
                    break;
                case OPT_MU:
		    param->shearModulus = atof(argValue);
                    break;
                case OPT_NU:
		    param->pois = atof(argValue);
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
                case OPT_LEVELS:
                    inArgs->numLevels = atoi(argValue);
                    break;
                case OPT_MAT_TYPE:
                    if (StrEquiv(argValue, "bcc")) {
                        param->materialType = MAT_TYPE_BCC;
                        strcpy(param->materialTypeName, "BCC");
                    } else if (StrEquiv(argValue, "fcc")) {
                        param->materialType = MAT_TYPE_FCC;
                        strcpy(param->materialTypeName, "FCC");
                    } else if (StrEquiv(argValue, "hcp")) {
                        param->materialType = MAT_TYPE_HCP;
                        strcpy(param->materialTypeName, "HCP");
                    } else if (StrEquiv(argValue, "rhmobohedral")) {
                        param->materialType = MAT_TYPE_RHOMBOHEDRAL_VA;
                        strcpy(param->materialTypeName, "RHOMBOHEDRAL");
                    } else {
                        printf("Error: Unknown material type %s\n", argValue);
                        exit(1);
                    }
                    break;
                case OPT_OUTFILE:
                    strcpy(param->fmCorrectionTbl, argValue);
                    break;
                case OPT_HELP:
                    PrintHelp(argv[0]);
                    exit(0);
                    break;

            }  /* end switch (optList[j].optType) */
#endif
	} // end for(i = 1; i < argc; i++) {


/*
 *      Set remaining default values if permitted
 */
        if (inArgs->numLevels < 3) {
            inArgs->numLevels = 3;
        }

        if (param->maxSideX <= 0.0) {
            real8  maxSide, minSide;

            maxSide  = 1.0;
            minSide  = -maxSide;

            param->minSideX=minSide;
            param->minSideY=minSide;
            param->minSideZ=minSide;

            param->maxSideX=maxSide;
            param->maxSideY=maxSide;
            param->maxSideZ=maxSide;
        }

        if (param->fmCorrectionTbl[0] == 0) {
            sprintf(param->fmCorrectionTbl, "Stanfordfm-ctab.c%d.l%d.dat",
                    inArgs->FMMOrder, inArgs->numLevels);
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
static void InitValues(InputArgs_t *inArgs, Param_t *param)
{
        inArgs->FMMOrder = 3;
        inArgs->numLevels = 5;
        inArgs->cMatrixFile[0] = 0;


        param->fmExpansionOrder = 3;

        param->minSideX = -1.0;
        param->minSideY = -1.0;
        param->minSideZ = -1.0;

        param->maxSideX = -1.0;
        param->maxSideY = -1.0;
        param->maxSideZ = -1.0;

#ifdef ANISOTROPIC
	param->anisoHarmonicsNumTermsBase = 4;
	param->anisoNumPhiPoints = 50;
	param->anisoNumThetaPoints = 50;
#endif

        param->materialType = MAT_TYPE_BCC;
        strcpy(param->materialTypeName, "BCC");

        memset(param->fmCorrectionTbl, 0, sizeof(param->fmCorrectionTbl));

        return;
}

int main(int argc, char *argv[])
{
        int         pbc[3];
        real8       epsC = 1.0e-05;
        Home_t      home;
        Param_t     param;
        InputArgs_t inArgs;

        memset(&home, 0, sizeof(Home_t));
        memset(&param, 0, sizeof(Param_t));
        memset(&inArgs, 0, sizeof(InputArgs_t));

        home.param = &param;

/*
 *      For now, we're using PBC in all dimensions.
 */
        pbc[0] = 1;
        pbc[1] = 1;
        pbc[2] = 1;

        InitValues(&inArgs, &param);
        GetInArgs(argc, argv, &inArgs, &param);

/*
 *      If the user provided a file with the full elastic constant
 *      matrix, read it in, otherwise we'll need to calculate it from
 *      the elastic constants;  if the elastic constants have not been
 *      provided, we have a problem.
 */

#ifdef ANISOTROPIC 
        if (inArgs.cMatrixFile[0] != 0) {
            int   i;
            char  inLine[MAX_LINE_LEN];
            FILE *fpMatrix;

            fpMatrix = fopen(inArgs.cMatrixFile, "r");

            if (fpMatrix == (FILE *)NULL) {
                printf("Error opening file %s\n", inArgs.cMatrixFile);
                exit(1);
            }

            for (i = 0; i < 6; i++) {
                fgets(inLine, MAX_LINE_LEN, fpMatrix);
                if (sscanf(inLine, "%lf %lf %lf %lf %lf %lf",
                    &param.elasticConstantMatrix[i][0],
                    &param.elasticConstantMatrix[i][1],
                    &param.elasticConstantMatrix[i][2],
                    &param.elasticConstantMatrix[i][3],
                    &param.elasticConstantMatrix[i][4],
                    &param.elasticConstantMatrix[i][5]) != 6) {
                    printf("Error reading C matrix file %s\n",
                           inArgs.cMatrixFile);
                    exit(1);
                }
            }

            fclose(fpMatrix);

        } else {
            switch(param.materialType) {
                case MAT_TYPE_BCC:
                case MAT_TYPE_FCC:
                    if (fabs(param.C11) * fabs(param.C12) *
                        fabs(param.C44) <= epsC) {
                        printf("Error: The needed elastic constants were "
                               "not provided\n");
                        exit(1);
                    }
                    break;

                case MAT_TYPE_HCP:
                    if (fabs(param.C11) * fabs(param.C12) *
                        fabs(param.C13) * fabs(param.C33) *
                        fabs(param.C44) <= epsC) {
                        printf("Error: The needed elastic constants were "
                               "not provided\n");
                        exit(1);
                    }
                    break;
                case MAT_TYPE_RHOMBOHEDRAL_VA:
                    printf("Error: The full elastic constant matrix must "
                           "be provided for rhombohedral materials\n");
                    exit(1);
            }

            SetElasticConstantMatrix(&home);
	}
#endif

        printf("\n");
        printf("Simulation size:  %d\n",
               (int)(param.maxSideX - param.minSideX));
        printf("Num levels:       %d\n", inArgs.numLevels);
        printf("FMM order:        %d\n", inArgs.FMMOrder);
        printf("Correction table: %s\n", param.fmCorrectionTbl);
        printf("Material type:    %s\n", param.materialTypeName);
	printf("Elastic Constants:\n");

#ifdef ANISOTROPIC
	printf("                   C11=%8.5f C12=%8.5f C44=%8.5f\n",param.C11, param.C12, param.C44);
	printf("qMax :            %d\n",param.anisoHarmonicsNumTermsBase);
	printf("nphi :            %d\n",param.anisoNumPhiPoints);
	printf("ntheta :          %d\n",param.anisoNumThetaPoints);
#else
	printf("                   mu = %8.5f nu=%8.5f\n",param.shearModulus,param.pois);
#endif

        printf("\n");

        Stanford_CreateTableCorrection(&home, inArgs.FMMOrder, inArgs.numLevels, pbc);

   return(0);
}

#else

int main(int argc, char *argv[])
{
   printf("Turn on BBFMM or UNIFORMFMM methods to obtain FMM tables for these methods\n");

   return(0);
}
#endif
