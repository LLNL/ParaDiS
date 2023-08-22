#pragma once

#ifndef _PDS_PARADIS_GEN_H
#define _PDS_PARADIS_GEN_H

/*
 *      For each type of nodal configuration that can be created
 *      define a name and integer id to be associated with the type
 */
#define FTYPE_SCREW		0
#define FTYPE_FCCSCREW		1
#define FTYPE_FINITE_MIXED	2
#define FTYPE_PRISMATIC_LOOP	3
#define FTYPE_FCCPRISMATIC_LOOP 4
#define FTYPE_FRANK_READ	5
#define FTYPE_FCC		6
#define FTYPE_FCC_IRRAD		7
#define FTYPE_EDGE		8
#define FTYPE_RHOMBOHEDRAL_VA   9
#define FTYPE_HCP		10
#define FTYPE_MAX		11

#define FNAME_SCREW		"screw"
#define FNAME_FCCSCREW		"fccscrew"
#define FNAME_EDGE		"edge"
#define FNAME_FINITE_MIXED	"finite-mixed"
#define FNAME_PRISMATIC_LOOP	"prismatic-loop"
#define FNAME_FCCPRISMATIC_LOOP "fccprismatic-loop"
#define FNAME_FRANK_READ	"frank-read-src"
#define FNAME_FCC		"fcc"
#define FNAME_FCC_IRRAD		"fcc-irrad"
#define FNAME_RHOMBOHEDRAL_VA 	"rhombohedral-va"
#define FNAME_HCP               "hcp"


/*
 *      Define a structure to hold a nodal configuration type, name, and
 *      a pointer to the function to invoke to create that type of nodal
 *      configuration.
 */
typedef struct {
        int   funcType;
        char  *funcName;
/*
        void  (* func)();
*/
} FuncData_t;


/*
 *      Define an integer identifier to be associated with each
 *      possible command line argument.  To be used to index the
 *      option-specific data in the optList array below.
 */
typedef enum {
    OPT_COVERA = 0,
    OPT_CUBEL,
    OPT_FRLEN,
    OPT_HELP,
    OPT_HEXSIZE,
    OPT_LOOPTYPE,
    OPT_MAXSEG,
    OPT_NCHAINS,
    OPT_NFRSRCS,
    OPT_NLOOPS,
    OPT_NODES_PER_FILE,
    OPT_OUTFILE,
    OPT_PBC,
    OPT_RADIUS,
    OPT_SEED,
    OPT_TYPE,
    OPT_VACANCY,
    OPT_XSURF,
    OPT_YSURF,
    OPT_ZSURF,
    OPT_MAX
} ParadisGenArgs_t;

/*
 *      Define a structure to hold a command line option's id (type),
 *      name, the shortest possible unique abbreviation of the option
 *      name, and a flag indicating if the option is paired with
 *      a value or not.
 */
typedef struct {
        int   optType;
        char  *optName;
        int   optMinAbbrev;
        int   optPaired;
} Option_t;


/*
 *      Define a structure containing all items corresponding to
 *      all command line options that have associated values.  This
 *      gives us an easy way to pass the command line arg values around
 *      to various functions.
 */
typedef struct {
        real8   cOVERa;
        int     cubeLength;
        int     frLenMin;
        int     frLenMax;
        int     hexl;
        int     interstitialLoops;
        int     loopType;
        real8   maxSegLen;
        int     numChains;
        int     numFRSrcs;
        int     numLoops;
        char    *outputFile;
        int     pbcVal;
        real8   radius;
        int     seed;
        int     type;
        int     targetNodesPerFile;
        real8   xSurf[2];
        real8   ySurf[2];
        real8   zSurf[2];
/*
 *      The following are not input variables, but it is useful to include
 *      them in this structure to pass them around with the values from
 *      the command line arguments.
 */
        int     useSegmentedDataFile;
        int     nodesWrittenToFile;
        int     currDataFileSegment;
} InArgs_t;

/*
 *      Prototype the various global functions needed during the initial
 *      problem generation.
 */
void    InitRemesh(InData_t *inData, int domValue, int startIndex);
Node_t *LookupInNode(InData_t *inData, Tag_t *tag, int *index);
real8   randm(int *seed);
void    WriteInitialNodeData(Home_t *home, InData_t *inData, InArgs_t *inArgs,
            int lastBlock);

void CreateFCCScrewConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                          real8 *totDislocLen);
void CreateEdges(Home_t *home, InData_t *inData, InArgs_t *inArgs,
            real8 *totDislocLen);
void CreateFCCConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
            real8 *totDislocLen);
void CreateFCCIrradConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
            real8 *totDislocLen);
void CreateFiniteMixedConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
            real8 *totDislocLen);
void CreateFRSource(Home_t *home, InData_t *inData, InArgs_t *inArgs,
            real8 *totDislocLen);
void CreateHCPLoop(Home_t *home, InData_t *inData, InArgs_t *inArgs,
            real8 *totDislocLen);
void CreatePrismaticLoop(Home_t *home, InData_t *inData, InArgs_t *inArgs,
            real8 *totDislocLen);
void CreateFCCPrismaticLoop(Home_t *home, InData_t *inData, InArgs_t *inArgs, 
                            real8 *totDislocLen);
void CreateScrewConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
            real8 *totDislocLen);
void RhombohedralVaConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
            real8 *totDislocLen);

#endif /* _ParadisGen_h */
