/***************************************************************************
 *
 *      Module:      StressTableGen.c
 *      Description: Contains all functions needed to compute and
 *                   write the two stress tables used for computing
 *                   farfield stress from within the primary image
 *                   and for periodic images if the fast-multipole
 *                   code is not being used.
 *
 *      Includes functions:
 *          CalRijm()
 *          CheckArgs()
 *          ComputeRijmTable()
 *          GetArgs()
 *          InitArgs()
 *          main()
 *          PrintHelp()
 *          Usage()
 *          WriteRijm()
 *          WriteRijmPBC()
 *
 *      Usage:  stresstablegen [-help] [-size <xlen[,ylen,zlen]>]      \
 *                             [-nimgx <numximages>] [-nimgy <numyimages>] \
 *                             [-nimgz <numzimages>] [-nx <nxval>          \
 *                             [-ny <nyval> [-nz <nzval> [-pbc|-nopbc] 
 *
 *      Overview:
 *
 *      ...
 *
 ***************************************************************************/
#include <string.h>
#include "Home.h"
#include <math.h>


/*
 *      Indices of 3rd derivatives
 */
static int DER3IND[10][3] = {{0,0,0}, {1,1,1}, {2,2,2},
                             {0,0,1}, {0,0,2}, {1,1,0},
                             {1,1,2}, {2,2,0}, {2,2,1},
                             {0,1,2} };

static int DER3INVIND[3][3][3]={{{0,3,4}, {3,5,9}, {4,9,7}},
                                {{3,5,9}, {5,1,6}, {9,6,8}},
                                {{4,9,7}, {9,6,8}, {7,8,2}}};
/*
 *      Stress tables.
 */
static real8 *rijmTable;
static real8 *rijmPBCTable;

/*
 *      Define an integer identifier to be associated with each
 *      posible command line argument.  To be used to index the
 *      option-specific data in the optList array below.
 */
#define OPT_HELP      0
#define OPT_SIZE      1
#define OPT_NIMGX     2
#define OPT_NIMGY     3
#define OPT_NIMGZ     4
#define OPT_NX        5
#define OPT_NY        6
#define OPT_NZ        7
#define OPT_NOPBC     8
#define OPT_OUTFILE   9
#define OPT_PBC      10
#define OPT_MAX      11

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


Option_t optList[OPT_MAX] = {
        { OPT_HELP,     "help",    1,  0 },
        { OPT_SIZE,     "size",    1,  1 },
        { OPT_NIMGX,    "nimgx",   5,  1 },
        { OPT_NIMGX,    "nimgy",   5,  1 },
        { OPT_NIMGX,    "nimgz",   5,  1 },
        { OPT_NX,       "nx",      2,  1 },
        { OPT_NY,       "ny",      2,  1 },
        { OPT_NZ,       "nz",      2,  1 },
        { OPT_NOPBC,    "nopbc",   2,  0 },
        { OPT_OUTFILE,  "outfile", 1,  1 },
        { OPT_PBC,      "pbc",     1,  0 }
};

/*
 *      Define a structure containing all items corresponding to
 *      command line options that have assoictaed values.  This
 *      gives us an easy way to pass the command line arg values
 *      around to various functions.
 */
typedef struct {
        int   pbc;
        int   nx;
        int   ny;
        int   nz;
        int   nimagex;
        int   nimagey;
        int   nimagez;
        real8 xLen, yLen, zLen;
        char  *rijmFileName;
        char  *rijmPBCFileName;
} InArgs_t;


/*---------------------------------------------------------------------------
 *
 *      Function:       Usage
 *      Description:    Print out a brief message indicating the possible
 *                      command line options and terminated.
 *
 *      Arguments:
 *              program         name of the program being executed
 *
 *-------------------------------------------------------------------------*/
static void Usage(char *program)
{
        printf("Usage:  %16s [-help] [-size <xlen[,ylen,zlen]>] \n",
               program);
        printf("                       [-nimgx <numximages>]\n");
        printf("                       [-nimgy <numyimages>]\n");
        printf("                       [-nimgz <numzimages>] [-nx <nxval>]\n");
        printf("                       [-ny <nyval>] [-nz <nzval>] \n");
        printf("                       [-pbc|-nopbc]\n");

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

        printf("    Options may be abbreviated to the shortest\n");
        printf("    non-ambiguous abbreviation of the option.  For\n");
        printf("    example, '-no' is a valid abbreviation for \n");
        printf("    '-nopbc' however, '-n' could refer to \n");
        printf("    numerous options and hence is invalid.\n\n");

        printf("Options:\n\n");
        printf("  -help             Prints this help information.\n\n");
        printf("  -size     Define the lengths (units of b) of the sides\n");
        printf("            of the problem space to be generated.  If\n");
        printf("            a single value is provided, the problem space\n");
        printf("            is assumed to be cubic, otherwise this \n");
        printf("            should contain 3 comma-delimited values \n");
        printf("            corresponding to the lengths of the problem\n");
        printf("            space in X, Y and Z dimensions respectively.\n\n");
        printf("  -nimgx    Number of periodic images of the problem\n");
        printf("            in the X dimension.\n\n");
        printf("  -nimgy    Number of periodic images of the problem\n");
        printf("            in the Y dimension.\n\n");
        printf("  -nimgz    Number of periodic images of the problem\n");
        printf("            in the Z dimension.\n\n");
        printf("  -nopbc    Indicates this invocation of the utilitiy is\n");
        printf("            to create the stress table used for farfield\n");
        printf("            stress within the primary image\n\n");
        printf("  -nx       For future use\n\n");
        printf("  -ny       For future use\n\n");
        printf("  -nz       For future use\n\n");
        printf("  -outfile  Name of the file into which to write the\n");
        printf("            table. \n\n");
        printf("  -pbc      Indicates this invocation of the utilitiy is\n");
        printf("            to create the stress table used for periodic\n");
        printf("            image stress\n\n");

        exit(0);
}


static void CalRijm(real8 rt[3], real8 Rijm[10])
{
        int i, j, m, k;
        real8 R, R2, R3;
        
        if ((fabs(rt[0]) < 1e-10) && (fabs(rt[1]) < 1e-10) &&
            (fabs(rt[2]) < 1e-10)) {

            for (k=0;k<10;k++) {
                Rijm[k]=0;
            }

            return;
        }
        
        R2 = rt[0]*rt[0] + rt[1]*rt[1] + rt[2]*rt[2];
        R  = sqrt(R2);
        R3 = R*R2;
            
        for (k=0;k<10;k++) {

            i = DER3IND[k][0];
            j = DER3IND[k][1];
            m = DER3IND[k][2];

            Rijm[k] = (-rt[m]*(i==j) - rt[j]*(i==m) - rt[i]*(j==m) +
                       3*rt[i]*rt[j]*rt[m]/R2) / R3;
        }

        return;
}


static void ComputeRijmTable(InArgs_t *inArgs)
{
        int   i, j, k, m;
        int   imgx, imgy, imgz;
        int   pbc, NX, NY, NZ, NIMGX, NIMGY, NIMGZ;
        int   offset1, offset2, offset3;
        real8 Lx, Ly, Lz;
        real8 rt[3],Rijmarray[10];
        real8 dRijmdx[10],dRijmdy[10],dRijmdz[10];
        
        pbc = inArgs->pbc;

        NX = inArgs->nx;
        NY = inArgs->ny;
        NZ = inArgs->nz;

        NIMGX = inArgs->nimagex;
        NIMGY = inArgs->nimagey;
        NIMGZ = inArgs->nimagez;

        Lx = inArgs->xLen;
        Ly = inArgs->yLen;
        Lz = inArgs->zLen;

        offset3 = 10;
        offset2 = 10 * inArgs->nz;
        offset1 = 10 * inArgs->ny * inArgs->nz;

        for (m=0;m<10;m++) {
            dRijmdx[m]=dRijmdy[m]=dRijmdz[m]=0;
        }

        for (imgx=-NIMGX;imgx<=NIMGX;imgx++) {
            for (imgy=-NIMGY;imgy<=NIMGY;imgy++) {
                for (imgz=-NIMGZ;imgz<=NIMGZ;imgz++) {

                    rt[0]=(-0.5+imgx)*Lx;
                    rt[1]=(     imgy)*Ly;
                    rt[2]=(     imgz)*Lz;
                    CalRijm(rt, Rijmarray);
                    for (m=0;m<10;m++) {
                        dRijmdx[m]-=Rijmarray[m];
                    }

                    rt[0]=( 0.5+imgx)*Lx;
                    rt[1]=(     imgy)*Ly;
                    rt[2]=(     imgz)*Lz;
                    CalRijm(rt, Rijmarray);
                    for (m=0;m<10;m++) {
                        dRijmdx[m]+=Rijmarray[m];
                    }
                
                    rt[0]=(     imgx)*Lx;
                    rt[1]=(-0.5+imgy)*Ly;
                    rt[2]=(     imgz)*Lz;
                    CalRijm(rt, Rijmarray);
                    for (m=0;m<10;m++) {
                        dRijmdy[m]-=Rijmarray[m];
                    }

                    rt[0]=(     imgx)*Lx;
                    rt[1]=( 0.5+imgy)*Ly;
                    rt[2]=(     imgz)*Lz;
                    CalRijm(rt, Rijmarray);
                    for (m=0;m<10;m++) {
                        dRijmdy[m]+=Rijmarray[m];
                    }
        
                    rt[0]=(     imgx)*Lx;
                    rt[1]=(     imgy)*Ly;
                    rt[2]=(-0.5+imgz)*Lz;
                    CalRijm(rt, Rijmarray);
                    for (m=0;m<10;m++) {
                        dRijmdz[m]-=Rijmarray[m];
                    }

                    rt[0]=(     imgx)*Lx;
                    rt[1]=(     imgy)*Ly;
                    rt[2]=( 0.5+imgz)*Lz;
                    CalRijm(rt, Rijmarray);
                    for (m=0;m<10;m++) {
                        dRijmdz[m]+=Rijmarray[m];
                    }
                }
            }
        }
        
            
        for (i=0;i<NX;i++) {
            for (j=0;j<NY;j++) {
                for (k=0;k<NZ;k++) {
                    if (!pbc) {
                        for (m=0;m<10;m++) {
                            rijmTable[i*offset1+j*offset2+k*offset3+m]=0;
                        }
                    } else {
                        for (m=0;m<10;m++) {
                            rijmPBCTable[i*offset1+j*offset2+k*offset3+m]=0;
                        }
                    }
            
                    for (imgx=-NIMGX;imgx<=NIMGX;imgx++) {
                        for (imgy=-NIMGY;imgy<=NIMGY;imgy++) {
                            for (imgz=-NIMGZ;imgz<=NIMGZ;imgz++) {
                                if (!pbc) {
                                    if ((imgx==0) &&
                                        (imgy==0) &&
                                        (imgz==0)) {
                                        continue;
                                    }
                                }

                                rt[0]=Lx*(-0.5+1.0*i/(NX-1)+imgx);
                                rt[1]=Ly*(-0.5+1.0*j/(NY-1)+imgy);
                                rt[2]=Lz*(-0.5+1.0*k/(NZ-1)+imgz);
                    
                                if (!((fabs(rt[0]) < 1e-10) &&
                                      (fabs(rt[1]) < 1e-10) &&
                                      (fabs(rt[2]) < 1e-10))) {
                                    CalRijm(rt, Rijmarray);
                                }
        
                                if (!pbc) {
                                    for (m=0;m<10;m++) {
                                        rijmTable[i*offset1+j*offset2+k*offset3+m]+=Rijmarray[m];
                                    }
                                } else {
                                    for (m=0;m<10;m++) {
                                        rijmPBCTable[i*offset1+j*offset2+k*offset3+m]+=Rijmarray[m];
                                    }
                                }
                            }
                        }
                    }
/*
 *                  Correct for slope
 */
                    if (!pbc) {
                        for (m=0;m<10;m++) {
                            rijmTable[i*offset1+j*offset2+k*offset3+m] -=
                                    dRijmdx[m]*(-0.5+1.0*i/(NX-1));
                            rijmTable[i*offset1+j*offset2+k*offset3+m] -=
                                    dRijmdy[m]*(-0.5+1.0*j/(NY-1));
                            rijmTable[i*offset1+j*offset2+k*offset3+m] -=
                                    dRijmdz[m]*(-0.5+1.0*k/(NZ-1));
                        }
                    } else {
                        for (m=0;m<10;m++) {
                            rijmPBCTable[i*offset1+j*offset2+k*offset3+m] -=
                                    dRijmdx[m]*(-0.5+1.0*i/(NX-1));
                            rijmPBCTable[i*offset1+j*offset2+k*offset3+m] -=
                                    dRijmdy[m]*(-0.5+1.0*j/(NY-1));
                            rijmPBCTable[i*offset1+j*offset2+k*offset3+m] -=
                                    dRijmdz[m]*(-0.5+1.0*k/(NZ-1));
                        }
                    }
                }

                printf("  <%06.3f%% complete>",
                       (double)(i*NY*NZ+j*NZ+NZ)/(double)(NX*NY*NZ)*100.0);
                printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                fflush(NULL);
            }
        }

        printf("                    \n");

        return;
}


static void WriteRijm(InArgs_t *inArgs)
{
        int   NX, NY, NZ;
        int   i, j, k, m;
        int   NIMGX, NIMGY, NIMGZ;
        int   offset1, offset2, offset3;
        real8 Lx, Ly, Lz;
        FILE  *fp;
        
        printf("Creating file %s...\n", inArgs->rijmFileName);
        
        NX = inArgs->nx;
        NY = inArgs->ny;
        NZ = inArgs->nz;
        
        NIMGX = inArgs->nimagex;
        NIMGY = inArgs->nimagey;
        NIMGZ = inArgs->nimagez;
        
        Lx = inArgs->xLen;
        Ly = inArgs->yLen;
        Lz = inArgs->zLen;
        
        offset3 = 10;
        offset2 = 10 * inArgs->nz;
        offset1 = 10 * inArgs->ny * inArgs->nz;

        fp=fopen(inArgs->rijmFileName,"w");
        
        fprintf(fp,"%d\n%d\n%d\n",NX,NY,NZ);
        fprintf(fp,"%e\n%e\n%e\n",Lx,Ly,Lz);
        fprintf(fp,"%d\n%d\n%d\n",NIMGX,NIMGY,NIMGZ);
            
        for (i=0;i<NX;i++) {
            for (j=0;j<NY;j++) {
                for (k=0;k<NZ;k++) {
                    for (m=0;m<10;m++) {
                        fprintf(fp, "%e\n",
                                rijmTable[i*offset1+j*offset2+k*offset3+m]);
                    }
                }
            }
        }

        fclose(fp);

        return;
}


static void WriteRijmPBC(InArgs_t *inArgs)
{
        int   NX, NY, NZ;
        int   i, j, k, m;
        int   NIMGX, NIMGY, NIMGZ;
        int   offset1, offset2, offset3;
        real8 Lx, Ly, Lz;
        FILE  *fp;
            
        printf("Creating file %s...\n", inArgs->rijmPBCFileName);
        
        NX = inArgs->nx;
        NY = inArgs->ny;
        NZ = inArgs->nz;
        
        NIMGX = inArgs->nimagex;
        NIMGY = inArgs->nimagey;
        NIMGZ = inArgs->nimagez;
        
        Lx = inArgs->xLen;
        Ly = inArgs->yLen;
        Lz = inArgs->zLen;
        
        offset3 = 10;
        offset2 = 10 * inArgs->nz;
        offset1 = 10 * inArgs->ny * inArgs->nz;

        fp = fopen(inArgs->rijmPBCFileName, "w");
        
        fprintf(fp, "%d\n%d\n%d\n", NX, NY, NZ);
        fprintf(fp, "%e\n%e\n%e\n", Lx, Ly, Lz);
        fprintf(fp, "%d\n%d\n%d\n", NIMGX, NIMGY, NIMGZ);
            
        for (i=0;i<NX;i++) {
            for (j=0;j<NY;j++) {
                for (k=0;k<NZ;k++) {
                    for (m=0;m<10;m++) {
                        fprintf(fp, "%e\n",
                                rijmPBCTable[i*offset1+j*offset2+k*offset3+m]);
                    }
                }
            }
        }

        fclose(fp);

        return;
}



static void InitArgs(InArgs_t *inArgs)
{
        inArgs->pbc =  1;

        inArgs->nx  = -1;
        inArgs->ny  = -1;
        inArgs->nz  = -1;

        inArgs->nimagex = -1;
        inArgs->nimagey = -1;
        inArgs->nimagez = -1;

        inArgs->xLen = -1.0;
        inArgs->yLen = -1.0;
        inArgs->zLen = -1.0;

        inArgs->rijmFileName    = (char *) "Rijm.tbl";
        inArgs->rijmPBCFileName = (char *) "RijmPBC.tbl";

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       GetArgs
 *      Description:    Parse the command line arguments.  Verify all options
 *                      are valid keywords and that all options requiring
 *                      associated values have them.  All user supplied values
 *                      will be stored in the inArgs structure, and a sanity
 *                      check on those values will be done elsewhere.
 *      Arguments
 *              argc    count of command line args
 *              argv    command line args
 *              inArgs  structure to hold user supplied values.  This
 *                      structure should have been populated with any
 *                      appropriate default values before this function
 *                      is called.
 *
 *-------------------------------------------------------------------------*/
static void GetArgs(int argc, char *argv[], InArgs_t *inArgs)
{
        int     i, j, k;
        char    *argName=0;
        char    *argValue=0;
        char    *token=0;

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
 *              Scan the array of valid options for the user supplied
 *              option name.  (This may be any unique abbreviation of
 *              one of any of the options.  If we don't find the option
 *              notify the user and terminate.
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
 *              Verify that there is an associated value if the specified
 *              option is supposed to be paired with a value.
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
 *              Do any option-specific processing...
 */
                switch (j)  {
                case OPT_HELP:
                        PrintHelp(argv[0]);
                        exit(0);
                        break;
                case OPT_SIZE:
/*
 *                      Allow a single value indicating a cubic problem
 *                      space with the specified dimensions, or 3
 *                      comma-delimited values specifying the lenth of
 *                      each dimension of the problem space explicitly
 */
                        token = strtok(argValue, ",");
                        inArgs->xLen = atoi(token);
                        token = strtok(NULL, ",");
                        if (token == (char *)NULL) {
                            inArgs->yLen = inArgs->xLen;
                            inArgs->zLen = inArgs->xLen;
                            break;
                        }
                        inArgs->yLen = atoi(token);
                        token = strtok(NULL, ",");
                        if (token == (char *)NULL) {
                            Usage(argv[0]);
                            exit(1);
                        }
                        inArgs->zLen = atoi(token);
                        break;
                case OPT_NX:
                        inArgs->nx = atoi(argValue);
                        break;
                case OPT_NY:
                        inArgs->ny = atoi(argValue);
                        break;
                case OPT_NZ:
                        inArgs->nz = atoi(argValue);
                        break;
                case OPT_NIMGX:
                        inArgs->nimagex = atoi(argValue);
                        break;
                case OPT_NIMGY:
                        inArgs->nimagey = atoi(argValue);
                        break;
                case OPT_NIMGZ:
                        inArgs->nimagez = atoi(argValue);
                        break;
                case OPT_NOPBC:
                        inArgs->pbc = 0;
                        break;
                case OPT_PBC:
                        inArgs->pbc = 1;
                        break;
                case OPT_OUTFILE:
                        inArgs->rijmFileName = argValue;
                        inArgs->rijmPBCFileName = argValue;
                        break;
                }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       CheckArgs
 *      Description:    Check the current values for the input arguments.
 *                      For each value not provided by the caller, apply
 *                      a default value if it has not yet been set.
 *
 *-------------------------------------------------------------------------*/
static void CheckArgs(InArgs_t *inArgs)
{
/*
 *      Set default to a cubic problem space of approximately 1 micron.
 */
        if (inArgs->xLen < 0) {
            inArgs->xLen = 3500;
            inArgs->yLen = 3500;
            inArgs->zLen = 3500;
        }


/*
 *      Default number of periodic images is 20 in X dimension and
 *      in Y and Z dimensions defaults to same value as in X dimension
 *      whether that be the default value or a user provided value.
 */
        if (inArgs->nimagex < 0) {
            inArgs->nimagex = 20;
        }

        if (inArgs->nimagey < 0) {
            inArgs->nimagey = inArgs->nimagex;
        }

        if (inArgs->nimagez < 0) {
            inArgs->nimagez = inArgs->nimagex;
        }

/*
 *      Default number of xxx is 51 in the X dimension and
 *      in Y and Z dimensions defaults to same value as in X dimension
 *      whether that be the default value or a user provided value.
 */
        if (inArgs->nx < 0) {
            inArgs->nx = 51;
        }

        if (inArgs->ny < 0) {
            inArgs->ny = inArgs->nx;
        }

        if (inArgs->nz < 0) {
            inArgs->nz = inArgs->nx;
        }

        return;
}


int main(int argc, char *argv[])
{
        int      tableSize;
        InArgs_t inArgs;

/*
 *      First initialize some stuff based on either user-provided
 *      command line arguments or hard-coded defaults.
 */
        InitArgs(&inArgs);
        GetArgs(argc, argv, &inArgs);
        CheckArgs(&inArgs);

/*
 *      Now allocate whichever table is being created
 */
        tableSize = (inArgs.nx * inArgs.ny * inArgs.nz * 10) * sizeof(real8);

        if (inArgs.pbc) {
            rijmPBCTable = (real8 *)calloc(1, tableSize);
        } else {
            rijmTable    = (real8 *)calloc(1, tableSize);
        }
        
/*
 *      And last compute the appropriate table and write it out to disk.
 */
        ComputeRijmTable(&inArgs);
 
        if (inArgs.pbc) {
            WriteRijmPBC(&inArgs);
        } else {
            WriteRijm(&inArgs);
        }

        exit(0);
        return(0);
}
