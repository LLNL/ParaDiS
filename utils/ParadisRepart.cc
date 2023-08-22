/*-------------------------------------------------------------------------
 *
 *      Module:       ParadisRepart.c
 *      Description:  This module contains the main function and various
 *                    support functions needed to replace the domain
 *                    decomposition/partitioning in an existing nodal
 *                    data file with a new decomposition.  The new
 *                    decomposition attempts to partition domains such
 *                    that all domains will require roughly the same
 *                    amount of processing time for local seg/seg force
 *                    calculations.
 * 
 *      Usage:  paradisrepart -infile inputfile -cells xcells[,ycells,zcells] \
 *                         -domains xdoms[,ydoms,zdoms] [-outfile outputfile] \
 *                         [-decompType <type>] [-help]
 *
 *      where
 *          
 *      -infile <inputFile>    Specifies the name of the nodal data file.
 *                             On success, the first nodal data file
 *                             segment will be copied to <inputFile>.bkup.
 *                             This command line argument is not optional.
 *
 *      -cells <xcells[,ycells,zcells]>  Specifies the number of
 *                             cells in each dimension.  If the
 *                             number of cells in the Y and Z dimensions
 *                             are not provided, they will default to the
 *                             same value as <xcells>. This command line
 *                             argument is not optional.
 *
 *      -domains <xdoms[,ydoms,zdoms]> Specifies the number of domains
 *                             in each dimension.  If the number of
 *                             domains in the Y and Z dimensions are
 *                             not provided, they will default to the
 *                             same value as <xdoms>.  This command line
 *                             argument is not optional.
 *
 *      -decompType <type>     Specifies type of domain decomposition:
 *                             If not specified, defaults to 1.  Valid
 *                             types are:
 *
 *                               1 == Recursive Sectioning
 *                               2 == Recursive Bisectioning
 *
 *      -help                  Causes the utility to display the command line
 *                             format and option descriptions then terminate.
 *
 *      -outfile <outputFile>  Specifies the name of the file into which
 *                             to write the new domain decomposition.  If
 *                             not specified, new domain decomposition
 *                             will be written into the file specified
 *                             by <inputFile>.
 *
 *
 *      All options may be abbreviated to the shortest non-ambiguous
 *      abbreviation of the option.  For example, '-c' is a valid abbreviation
 *      for '-cells'.  
 *
 *    ***********************************************************
 *    *                        FIX ME!                          *
 *    * Currently this utility assumes periodic boundaries are  *
 *    * used.  The code should probably be tweaked to allow the *
 *    * user to specify whether PBC is to be used or not.       *
 *    *                                                         *
 *    ***********************************************************
 *
 *-----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include "Home.h"
#include "Decomp.h"
#include "RBDecomp.h"
#include "Restart.h"


typedef struct {
        int   cellID;
        real8 coord[3];
} NodeData_t;


typedef struct {
        real8 cellWeight;
        real8 cellBlockWeight;
        int   cellNodeCount;
} CellData_t;

typedef enum {
        OPT_CELLS,
        OPT_DECOMPTYPE,
        OPT_DOMAINS,
        OPT_HELP,
        OPT_INFILE,
        OPT_OUTFILE,
        OPT_PBC,
        OPT_MAX
} RepartOpt_t;

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
 *      option          option          #characters     1 unless option
 *      type            name            in unique       has no associated
 *                                      abbreviation    value
 */
Option_t        optList[OPT_MAX] = {
        { OPT_CELLS,     "cells",       1, 1 },
        { OPT_DOMAINS,   "decompType",  2, 1 },
        { OPT_DOMAINS,   "domains",     2, 1 },
        { OPT_HELP,      "help",        1, 0 },
        { OPT_INFILE,    "infile",      1, 1 },
        { OPT_OUTFILE,   "outfile",     1, 1 },
        { OPT_PBC,       "pbc",         1, 1 }
};


static void Usage(char *prog)
{
        fprintf(stderr, "  Usage:  %s -infile inputfile        \\\n", prog);
        fprintf(stderr, "                -cells xcells[,ycells,zcells] \\\n");
        fprintf(stderr, "                -domains xdoms[,ydoms,zdoms]  \\\n");
        fprintf(stderr, "                [-decompType type]            \\\n");
        fprintf(stderr, "                [-outfile outputfile] [-help] \n");
        fprintf(stderr, "\n");

       return;
}

static void PrintHelp(char *prog)
{
    Usage(prog);


    fprintf(stderr, " where\n\n");
    fprintf(stderr, "    -infile <inputFile>    Specifies the name of the nodal data file.\n");
    fprintf(stderr, "                           On success, the first nodal data file\n");
    fprintf(stderr, "                           segment will be copied to <inputFile>.bkup.\n");
    fprintf(stderr, "                           This command line argument is not optional.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -cells <xcells[,ycells,zcells]>  Specifies the number of\n");
    fprintf(stderr, "                           cells in each dimension.  If the\n");
    fprintf(stderr, "                           number of cells in the Y and Z dimensions\n");
    fprintf(stderr, "                           are not provided, they will default to the\n");
    fprintf(stderr, "                           same value as <xcells>. This command line\n");
    fprintf(stderr, "                           argument is not optional.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -domains <xdoms[,ydoms,zdoms]> Specifies the number of domains\n");
    fprintf(stderr, "                           in each dimension.  If the number of\n");
    fprintf(stderr, "                           domains in the Y and Z dimensions are\n");
    fprintf(stderr, "                           not provided, they will default to the\n");
    fprintf(stderr, "                           same value as <xdoms>.  This command line\n");
    fprintf(stderr, "                           argument is not optional.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -decompType <type>     Specifies type of domain decomposition:\n");
    fprintf(stderr, "                           If not specified, defaults to 1.  Valid\n");
    fprintf(stderr, "                           types are:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "                             1 == Recursive Sectioning\n");
    fprintf(stderr, "                             2 == Recursive Bisectioning\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -help                  Causes the utility to display the command line\n");
    fprintf(stderr, "                           format and option descriptions then terminate.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    -outfile <outputFile>  Specifies the name of the file into which\n");
    fprintf(stderr, "                           to write the new domain decomposition.  If\n");
    fprintf(stderr, "                           not specified, new domain decomposition\n");
    fprintf(stderr, "                           will be written into the file specified\n");
    fprintf(stderr, "                           by <inputFile>.\n");


    return;
}


static void GetInArgs(int argc, char *argv[], Param_t *param,
                      char **baseFileName, char **outFileName)
{
        int  i, j;
        char *argName=0, *argValue=0, *token=0, *suffix=0;

        *baseFileName = (char *)NULL;
        *outFileName  = (char *)NULL;

        param->nXcells = -1;
        param->nYcells = -1;
        param->nZcells = -1;

        param->nXdoms  = -1;
        param->nYdoms  = -1;
        param->nZdoms  = -1;

        param->decompType = 1;

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
            if (optList[j].optPaired) {
                if (i+1 >= argc) {
                    Usage(argv[0]);
                    exit(1);
                } else {
                    argValue = argv[++i];
                }
            }
/*
 *          Do any option-specific processing...
 */
            switch (j)  {
            case OPT_CELLS:
                token = strtok(argValue, ",");
                param->nXcells = atoi(token);
                token = strtok(NULL, ",");
                if (token == (char *)NULL) {
                    param->nYcells = param->nXcells;
                    param->nZcells = param->nXcells;
                    break;
                }
                param->nYcells = atoi(token);
                token = strtok(NULL, ",");
                if (token == (char *)NULL) {
                    Usage(argv[0]);
                    exit(1);
                }
                param->nZcells = atoi(token);
                break;

            case OPT_DECOMPTYPE:
                param->decompType = atoi(argValue);
                break;

            case OPT_DOMAINS:
                token = strtok(argValue, ",");
                param->nXdoms = atoi(token);
                token = strtok(NULL, ",");
                if (token == (char *)NULL) {
                    param->nYdoms = param->nXdoms;
                    param->nZdoms = param->nXdoms;
                    break;
                }
                param->nYdoms = atoi(token);
                token = strtok(NULL, ",");
                if (token == (char *)NULL) {
                    Usage(argv[0]);
                    exit(1);
                }
                param->nZdoms = atoi(token);
                break;

            case OPT_HELP:
                PrintHelp(argv[0]);
                exit(0);
                break;

            case OPT_PBC:
                printf("Option not yet supported. PBC assumed for now\n");
                break;

            case OPT_INFILE:
                *baseFileName = argValue;
                suffix = (*baseFileName + strlen(*baseFileName)) - 2;
                if (strcmp(suffix, ".0") == 0) *suffix = 0;
                break;

            case OPT_OUTFILE:
                *outFileName = argValue;
                break;

            default:
                Usage(argv[0]);
                exit(1);
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    IndexToCellID
 *      Description: Given a cell's X, Y and Z indices, return the
 *                   integer cellID 
 *
 *      Arguments:
 *          param      Pointer to Param_t structure containing the
 *                     parameters defining the problem space, including
 *                     cell geometry.
 *          cellX      Cell's X index
 *          cellY      Cell's Y index
 *          cellZ      Cell's Z index
 *
 *      Returns:  Cell's integer ID
 *
 *------------------------------------------------------------------------*/
static int IndexToCellID(Param_t *param, int cellX, int cellY, int cellZ)
{
        int cellID;

        cellID = cellZ                +
                 param->nZcells*cellY +
                 param->nZcells*param->nYcells*cellX;

        return(cellID);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetNodesCellID
 *      Description: Given the specified x,y,z location return to the
 *                   caller the ID of the cell containing the point.
 *
 *      Arguments:
 *          param        Pointer to Param_t structure containing the
 *                       parameters defining the problem space, including
 *                       cell geometry.
 *          cellSizeX    Size of each cell in the X dimension
 *          cellSizeY    Size of each cell in the Y dimension
 *          cellSizeZ    Size of each cell in the Z dimension
 *          x,y,z        Coordinates of the point in question
 *
 *------------------------------------------------------------------------*/
static int GetNodesCellID(Param_t *param, real8 cellSizeX, real8 cellSizeY,
                          real8 cellSizeZ, real8 x, real8 y, real8 z)
{
        int cellID, cellX, cellY, cellZ;

	x -= param->minSideX;
	y -= param->minSideY;
	z -= param->minSideZ;

        cellX = (int)(x/cellSizeX);
        cellY = (int)(y/cellSizeY);
        cellZ = (int)(z/cellSizeZ);

/*
 *      We have to force any nodes exactly on the upper boundaries
 *      into the proper cell...
 */
        if (cellX >= param->nXcells) cellX--;
        if (cellY >= param->nYcells) cellY--;
        if (cellZ >= param->nZcells) cellZ--;

        cellID = IndexToCellID(param, cellX, cellY, cellZ);
        
        return(cellID);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetWeighting
 *      Description: For each cell, sets the cellblock (3X3X3 block
 *                   of cells centered around the cell in question)
 *                   weight.  Note: total weighting for the whole problem
 *                   is the product of the cell weight and cell block
 *                   weight summed over all cells.
 *
 *                   NOTE: the cellBlockWeight values for each cell are
 *                   assumed to have been zero'ed out before this function
 *                   is invoked, and cell weights should have been
 *                   calculated when reading in the nodal data.
 *
 *      Arguments:
 *          param        Pointer to Param_t structure containing the
 *                       parameters defining the problem space, including
 *                       cell geometry.
 *          cellList     Array of cell specific data
 *
 *------------------------------------------------------------------------*/
static void GetWeighting(Param_t *param, CellData_t *cellList)
{
        int i, j, k;
        int cx1, cy1, cz1;
        int cx2, cy2, cz2;
        int cellID1, cellID2;
        int numXCells, numYCells, numZCells;

        numXCells = param->nXcells;
        numYCells = param->nYcells;
        numZCells = param->nZcells;

/*
 *      Loop over every cell...
 */
        for (cx1 = 0; cx1 < numXCells; cx1++) {
            for (cy1 = 0; cy1 < numYCells; cy1++) {
                for (cz1 = 0; cz1 < numZCells; cz1++) {

                    cellID1 = IndexToCellID(param, cx1, cy1, cz1);

/*
 *                  The cell block weight is essentially an estimate of
 *                  the number of force calculations that will be needed
 *                  for a node in the cell.  The nature of the current
 *                  force calculations in the parallel code are such that
 *                  for any pair of segments, the force calc will only
 *                  be done once by the domain owning the dominant
 *                  segment, with results being communicated to remote
 *                  domains as necessary.  In order to emulate that
 *                  the cell block weight will be the sum of the cell weights
 *                  of all cells with a higher cell ID (taking periodic
 *                  boundaries into consideration, of course).
 *
 *                  Assume periodic bounaries for now.  Probably should
 *                  fix this to handle non-PBC also, though.
 */
 
                    for (i = cx1; i <= cx1+1; i++) {
                        cx2 = i;
                        cx2 = (cx2 >= numXCells) ? 0 : cx2;

                        for (j = cy1; j <= cy1+1; j++) {
                            cy2 = j;
                            cy2 = (cy2 >= numYCells) ? 0 : cy2;

                            for (k = cz1; k <= cz1+1; k++) {
                                cz2 = k;
                                cz2 = (cz2 >= numZCells) ? 0 : cz2;

                                cellID2 = IndexToCellID(param, cx2, cy2, cz2);
                                cellList[cellID1].cellBlockWeight +=
                                        cellList[cellID2].cellWeight;

                            }  /* for (k = cz1-1; ...) */
                        }  /* for (j = cy1-1; ...) */
                    }  /* for (i = cx1-1; ...) */

                }  /* for (cz1 = 0; ...) */
            }  /* for (cy1 = 0; ...) */
        }  /* for (cx1 = 0; ...) */

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Functions:   XCmp, YCmp, ZCmp
 *      Description: These function are simple comparison functions
 *                   which operate respectively on the X, Y and Z
 *                   coordinate components of the NodeData_t structure.
 *                   These functions are compatible with system sort
 *                   functions such as qsort().
 * 
 *      Returns:     -1  if (a) <  (b)
 *                    0  if (a) == (b)
 *                    1  if (a) >  (b)
 *
 *------------------------------------------------------------------------*/
static int XCmp(const void *a, const void *b)
{
	NodeData_t *nodeA = (NodeData_t *)a, *nodeB = (NodeData_t *)b;

	if (nodeA->coord[X] < nodeB->coord[X]) return(-1);
	return(nodeA->coord[X] > nodeB->coord[X]);
}

static int YCmp(const void *a, const void *b)
{
	NodeData_t *nodeA = (NodeData_t *)a, *nodeB = (NodeData_t *)b;

	if (nodeA->coord[Y] < nodeB->coord[Y]) return(-1);
	return(nodeA->coord[Y] > nodeB->coord[Y]);
}

static int ZCmp(const void *a, const void *b)
{
	NodeData_t *nodeA = (NodeData_t *)a, *nodeB = (NodeData_t *)b;

	if (nodeA->coord[Z] < nodeB->coord[Z]) return(-1);
	return(nodeA->coord[Z] > nodeB->coord[Z]);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetBounds
 *      Description: Given a set of nodal coordinates, the number of
 *                   slices (and number of nodes per slice) over which
 *                   to spread the nodes, calculate the boundaries.
 *
 *      Arguments:
 *          cellList   array of structures (1 per cell) containing
 *                     the cell and cell block weighting values.
 *          nodeList   array of structures containing all data specific
 *                     to the node that is required for the repartitioning
 *          dim        Indicates the element of the nodal coordinates
 *                     by which nodes in <nodeList> have been sorted.
 *                         0 == X coordinate
 *                         1 == Y coordinate
 *                         2 == Z coordinate
 *          numNodes   number of structures in array <nodeList>
 *          numSlices  Number of sections into which the specified area
 *                     should be sliced along the <dim> dimension
 *          minSide    Minimum coordinate in the dimension specified by
 *                     <dim> of the volume encompassing all nodes in
 *                     the <nodeList> array.
 *          maxSide    Maximum coordinate in the dimension specified by
 *                     <dim> of the volume encompassing all nodes in
 *                     the <nodeList> array.
 *          bounds     array in which to return to the caller the
 *                     calculated boundaries.  This array must have 
 *                     <numSlices>+1 elements.  On return, bounds[0] will
 *                     always be minSide and bounds[numSlices] will always
 *                     be maxSide.
 *          sliceNodeCount  array indicating the number of nodes that should
 *                     be contained in each of the <numSlices> slices.
 *
 *------------------------------------------------------------------------*/
static void GetBounds(CellData_t *cellList, NodeData_t *nodeList, int dim,
                      int numNodes, int numSlices, real8 minSide, real8 maxSide,
                      real8 *bounds, int *sliceNodeCount)
{
        int   i, j;
        int   cellID, nextNode, remainingParts, nodesInPart, remainingNodes;
        real8  minCoord, maxCoord, nodeCoord;
        real8  optimalWeight, totalWeight, partWeight, remainingWeight;
        real8 minBound, maxBound;

/*
 *      lower and upper boundaries will always be set to the 
 *      minimum and maximum boundaries provided by caller.
 */
        bounds[0] = minSide;
        bounds[numSlices] = maxSide;

/*
 *      Find the total weighting for all nodes in this section of
 *      the problem space.
 */
        totalWeight = 0.0;

        for (i = 0; i < numNodes; i++) {
            cellID = nodeList[i].cellID;
            totalWeight += cellList[cellID].cellBlockWeight;
        }

        remainingNodes = numNodes;
        remainingWeight = totalWeight;
        nextNode = 0;

/*
 *      First and last boundaries are already set, so loop to
 *      set remaining numSlices-1 boundaries.
 */
        for (i = 1; i < numSlices; i++) {

            minCoord = bounds[i-1];
            if (numNodes == 0) {
                maxCoord = bounds[numSlices];
            } else {
                maxCoord = nodeList[nextNode].coord[dim];
            }

/*
 *          Optimal weight for this partition is summed weight of
 *          all nodes in the portion of this slice not yet partitiond
 *          divided by the number of pieces into which the remaining
 *          portion is to be partitioned.
 */
            remainingParts = numSlices - (i - 1);
            optimalWeight  = remainingWeight / (real8)remainingParts;

/*
 *          Set minimum and maximum limits on the upper boundary
 *          for this partition.
 */
            minBound = bounds[i-1];
            maxBound = maxSide;

/*
 *          Start looping though the nodes to determine where
 *          the upper boundary for this partition should be.
 */
            partWeight = 0.0;
            nodesInPart = 0;

            for (j = nextNode; j < numNodes; j++) {

                nodeCoord = nodeList[j].coord[dim];
                cellID = nodeList[j].cellID;

/*
 *              If the next node's coordinate is below the minimum
 *              upper boundary for this partition, it must be included
 *              in the partition...
 */
                if (nodeCoord < minBound) {
                    partWeight += cellList[cellID].cellBlockWeight;
                    minCoord    = nodeCoord;
                    nextNode++;
                    nodesInPart++;
                    continue;
                }

/*
 *              If the next node's coordinate is above the maximum
 *              upper boundary for this partition, break out of the
 *              loop.  the boundary will be set between the coordinates
 *              of this node nd the previous node.
 */
                if (nodeCoord >= maxBound) {
                    break;
                }

/*
 *              As soon as the weighting for the current partition
 *              reaches the optimal, stop adding nodes, and break out
 *              to do the boundary calculation.
 */
                partWeight += cellList[cellID].cellBlockWeight;

                if (partWeight >= optimalWeight) {
                    minCoord = nodeCoord;
                    nextNode++;
                    nodesInPart++;
                    break;
                }

                minCoord = nodeCoord;
                nextNode++;
                nodesInPart++;

            }  /* for (j = nextNode; ...) */

/*
 *          Set boundary between minCoord and maxCoord, then
 *          apply limits?
 */
            if (nextNode == numNodes) {
                maxCoord = bounds[numSlices];
            } else {
                maxCoord = nodeList[nextNode].coord[dim];
            }

            bounds[i] = (minCoord + maxCoord) / 2.0;

            if (bounds[i] < minBound) {
                bounds[i] = minBound;
            } else if (bounds[i] > maxBound) {
                bounds[i] = maxBound;
            }
  
            remainingWeight -= partWeight;

            if (remainingWeight < 0) {
                remainingWeight = 0;
            }

            minCoord = bounds[i];
            sliceNodeCount[i-1] = nodesInPart;
            remainingNodes -= nodesInPart;

        }  /* for (i = 1; i < numSlices; ...) */

        sliceNodeCount[numSlices-1] = remainingNodes;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ReadDataFile
 *      Description: Read nodal data and ancillary data from the
 *                   nodal data file.
 *
 *      Arguments:
 *          fileName         Name of source data file.  NOTE: For
 *                           segmented data files, this may be either
 *                           the base file name, or the name of the 
 *                           first data file segment. (i.e. both
 *                           rsNNNN.data and rsNNNN.data.0 are
 *                           valid if rsNNNN.data consists of multiple
 *                           file segments.
 *          cellList         Array of cell specific data
 *          nodeListp        location in which to return to the caller a
 *                           pointer to the array of nodal coordinates read
 *                           from the data file
 *          nodeDataPos      location in which to return the position in
 *                           the file at which the nodal data begins.
 *
 *------------------------------------------------------------------------*/
static void ReadDataFile(Home_t *home, char *fileName, CellData_t *cellList,
                         NodeData_t **nodeListp, fpos_t *nodeDataPos)
{
        int     i, j, numArms;
        int     cellID;
        int     fileSeqNum = 0, binRead = 0;
        int     valType, numVals, tokenType, maxTokenLen, pIndex, okay;
        real8   cellSizeX, cellSizeY, cellSizeZ;
        char    token[256];
        char    inLine[512], baseFileName[256], segmentName[256];
        char    *suffix;
        NodeData_t *nodeList;
        Param_t *param;
        void    *valList;
        FILE    *fp;
        ParamList_t *dataParamList;

        param = home->param;
        dataParamList = home->dataParamList;
        maxTokenLen   = sizeof(token);

        memset(inLine, 0, sizeof(inLine));

/*
 *      First get the data file's base file name (i.e. no
 *      sequence number suffix)
 */
        memset(baseFileName, 0, sizeof(baseFileName));
        strncpy(baseFileName, fileName, sizeof(baseFileName)-1);

        suffix = (baseFileName + strlen(baseFileName)) - 2;
        if (strcmp(suffix, ".0") == 0) *suffix = 0;

        snprintf(segmentName, sizeof(segmentName), "%s", baseFileName);

/*
 *      Attempt to open the base data file.  If that doesn't work,
 *      try opening the first segment of that data file.  If we
 *      still can't open the file, commit suicide
 */
        if ((fp = fopen(segmentName, "r")) == (FILE *)NULL) {
            strcat(segmentName, ".0");
            if ((fp = fopen(segmentName, "r")) == (FILE *)NULL) {
                printf("Error %d opening file %s\n", errno, baseFileName);
                exit(1);
            }
        }

/*
 *      Get the first token.  This should either be a known
 *      parameter identifier, or a file version number.
 */
        tokenType = GetNextToken(fp, token, maxTokenLen);
        pIndex = LookupParam(dataParamList, token);

        if (pIndex < 0) {
/*
 *          If the token does not correspond to a known parameter, it
 *          should be the version number, so convert it to an integer
 *          and parse the rest of the file the old way (i.e. values
 *          are positional, not identified by name/value pairs)
 */
            param->dataFileVersion = atoi(token);

            if ((param->dataFileVersion < 1) ||
                (param->dataFileVersion > 3)) {
                Fatal("ReadNodeDatFile: Unsupported file version %d",
                      param->dataFileVersion);
            }

            ReadPreV4DataParams(home, fp, &home->decomp);

        } else {
/*
 *          Just go through the nodal data file reading all
 *          the associated parameters.  Need to do special
 *          processing of the domain decomposition, and when
 *          when we hit the 'nodaldata' identifier, just break
 *          out of this loop.
 */
            while ((tokenType != TOKEN_ERR) && (tokenType != TOKEN_NULL)) {

                if (pIndex >= 0) {
/*
 *                  Token represents a known parameter identifier, so
 *                  read the associated value(s).
 */
                    valType = dataParamList->varList[pIndex].valType;
                    numVals = dataParamList->varList[pIndex].valCnt;
                    valList = dataParamList->varList[pIndex].valList;
                    okay = GetParamVals(fp, valType, numVals, valList);
                    if (!okay) {
                        Fatal("Parsing Error obtaining values for "
                              "parameter %s\n",
                              dataParamList->varList[pIndex].varName);
                    }
                } else {
/*
 *                  Token does not represent one of the simple
 *                  parameters.  If it's not one of the identifiers
 *                  that needs special handling, skip it.
 */
                    if (strcmp(token, "domainDecomposition") == 0) {
/*
 *                      The minSide{XYZ} and maxSide{XYZ} values are
 *                      now redundant but until the rest of the code
 *                      is modified to remove references to these
 *                      values, we need to explicitly set them now.
 */
                        param->minSideX = param->minCoordinates[X];
                        param->minSideY = param->minCoordinates[Y];
                        param->minSideZ = param->minCoordinates[Z];

                        param->maxSideX = param->maxCoordinates[X];
                        param->maxSideY = param->maxCoordinates[Y];
                        param->maxSideZ = param->maxCoordinates[Z];

                        tokenType = GetNextToken(fp, token, maxTokenLen);
                        ReadDecompBounds(home, (void **)&fp, binRead,
                                         param->decompType, &home->decomp);

                   } else if (strcmp(token, "nodalData") == 0) {
/*
 *                      When we hit the nodal data, we can break
 *                      out of the loop because we are assuming
 *                      all other data file parameters have been
 *                      processed.  If they have not, we have a
 *                      problem since processing of the nodal data
 *                      requires the other parameters.
 *                      Note: Remainder of the file should just
 *                      contain " = " followed by the nodal data,
 *                      so be sure to skip the next token before
 *                      processing the nodal data.
 */
                        tokenType = GetNextToken(fp, token, maxTokenLen);
                        break;
                    } else {
/*
 *                      If the parameter is not recognized, skip the
 *                      parameter and any associated value(s).
 */
                        printf("Ignoring unknown data file parameter %s\n",
                               token);
                        valType = V_NULL;
                        numVals = 0;
                        valList = (void *)NULL;
                        okay = GetParamVals(fp, valType, numVals,
                                            valList);
                    }
                }
                tokenType = GetNextToken(fp, token, maxTokenLen);

                if ((tokenType == TOKEN_NULL)||(tokenType == TOKEN_ERR)) {
                    Fatal("Parsing error on file %s\n", segmentName);
                }
                pIndex = LookupParam(dataParamList, token);
            }
        }

/*
 *      Preserve the current value of the file pointer.  Later we'll
 *      be seeking back to this location in order to do a copy
 *      of the nodal data from this file to the new data file.
 */
        if (fgetpos(fp, nodeDataPos) != 0) {
            fprintf(stderr, "fgetpos() failed with error%d\n", errno);
            exit(1);
        }

/*
 *      When we read the nodal data in, it is possible some coordinates
 *      may be outside the primary image.  We need to set a few parameters
 *      related to the box size in order to properly fold the the nodes
 *      back into the primary image.
 */
        param->Lx = param->maxSideX - param->minSideX;
        param->Ly = param->maxSideY - param->minSideY;
        param->Lz = param->maxSideZ - param->minSideZ;

        cellSizeX = param->Lx / param->nXcells;
        cellSizeY = param->Ly / param->nYcells;
        cellSizeZ = param->Lz / param->nZcells;

/*
 *      Read only the nodal coordinates into memory and fold any
 *      nodes outside the problem space back into the primary image.
 *
 *      If we hit hit an EOF reading nodal data, (i.e. inLine
 *      is empty) and we're on the last data file segment to be
 *      read then there's a problem.  If we're not on the last
 *      file segment, close the current segment, open the next
 *      and continue reading.  As usual, abort on an open error.
 */
        nodeList = (NodeData_t *)malloc(param->nodeCount * sizeof(NodeData_t));

        for (i = 0; i < param->nodeCount; i++) {

            while (1) {
                Getline(inLine, sizeof(inLine), fp);
                if (inLine[0] == 0) {
                    fileSeqNum++;
                    if (fileSeqNum >= param->numFileSegments) {
                        Fatal("Unexpected EOF in nodal data file");
                    }
                    fclose(fp);
                    sprintf(segmentName, "%s.%d", baseFileName, fileSeqNum);
                    fp = fopen(segmentName, "r");
                    if (fp == (FILE *)NULL) {
                        Fatal("Error %d opening node data file %s",
                              errno, fileName);
                    }
                } else {
                    break;
                }
            }

/*
 *          Got information on the next node, so deal with it.
 */
            sscanf(inLine, "%*d,%*d %lf %lf %lf %d %*d",
                   &nodeList[i].coord[X], &nodeList[i].coord[Y],
                   &nodeList[i].coord[Z], &numArms);

            FoldBox(param, &nodeList[i].coord[X], &nodeList[i].coord[Y],
                    &nodeList[i].coord[Z]);

            for (j = 0; j < numArms; j++) {
                Getline(inLine, sizeof(inLine), fp);
                Getline(inLine, sizeof(inLine), fp);
            }
/*
 *          Determine the ID of the cell encompassing this node and
 *          increment that cell's node count.
 */
            cellID = GetNodesCellID(param, cellSizeX, cellSizeY, cellSizeZ,
                                    nodeList[i].coord[X], nodeList[i].coord[Y],
                                    nodeList[i].coord[Z]);

            cellList[cellID].cellWeight += (real8)(numArms-1);
            cellList[cellID].cellNodeCount++;
            nodeList[i].cellID = cellID;
        }

        *nodeListp = nodeList;

        fclose(fp);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    RepartRSDecomp
 *      Description: Generate a new domain decomposition using the
 *                   recusrive sectioning algorithm.  The repartioning
 *                   will be based on the provided list of nodes, the
 *                   cells encompassing those nodes, and the estimated
 *                   weight for each node.
 *
 *      Arguments:
 *          nodeList         array of nodal data read from the restart file
 *          numNodes         the number of structures contained in <nodeList> 
 *          cellList         Array of cell specific weighting data
 *
 *------------------------------------------------------------------------*/
static void RepartRSDecomp(Home_t *home, NodeData_t *nodeList,
                           CellData_t *cellList)
{
        int        xDom, yDom, zDom;
        int        numXdoms, numYdoms, numZdoms;
        int        xFirstIndex, yFirstIndex, zFirstIndex;
        int        xLastIndex, yLastIndex, zLastIndex;
        int        xNodeIndex, yNodeIndex, zNodeIndex;
        int        nodesInXSlice, nodesInYSlice;
        int        numNodes;
        int        *xSlices, *ySlices, *zSlices;
        real8      *xBound, *yBound, *zBound;
        Param_t    *param;
        RSDecomp_t *decomp;

        param = home->param;

        numXdoms = param->nXdoms;
        numYdoms = param->nYdoms;
        numZdoms = param->nZdoms;

        numNodes = param->nodeCount;

/*
 *      Allocate some temporary arrays
 */
        xBound = (real8 *)malloc((numXdoms+1) * sizeof(real8));
        yBound = (real8 *)malloc((numYdoms+1) * sizeof(real8));
        zBound = (real8 *)malloc((numZdoms+1) * sizeof(real8));

        xSlices = (int *)malloc(numXdoms * sizeof(int));
        ySlices = (int *)malloc(numYdoms * sizeof(int));
        zSlices = (int *)malloc(numZdoms * sizeof(int));

/*
 *      Allocate storage for the domain decomposition
 */
        home->decomp = (void *)calloc(1, sizeof(RSDecomp_t));
        decomp = (RSDecomp_t *)home->decomp;

        decomp->domBoundX = (real8 *)malloc((numXdoms+1) * sizeof(real8));
        decomp->domBoundY = (real8 **)malloc(numXdoms * sizeof(real8 *));
        decomp->domBoundZ = (real8 ***)malloc(numXdoms * sizeof(real8 **));

        for (xDom = 0; xDom < numXdoms; xDom++) {
            decomp->domBoundY[xDom] =
                    (real8 *)malloc((numYdoms+1) * sizeof(real8));
            decomp->domBoundZ[xDom] =
                    (real8 **)malloc(numYdoms * sizeof(real8 *));
            for (yDom = 0; yDom < numYdoms; yDom++) {
                decomp->domBoundZ[xDom][yDom] =
                        (real8 *)malloc((numZdoms+1) * sizeof(real8));
            }
        }

/*
 *      Sort the full list of nodal coordinates by the X component.
 *      Then carve up the problem space into slices along the X axis
 *      such that all slices contain roughly the same weighting 
 */
        xFirstIndex = 0;
        xLastIndex  = numNodes;
        xNodeIndex  = 0;

        qsort(nodeList, numNodes, sizeof(NodeData_t), XCmp);

        GetBounds(cellList, nodeList, X, numNodes, numXdoms, param->minSideX,
                  param->maxSideX, xBound, xSlices);

/*
 *      Now loop through each X slice, resorting coordinates for
 *      all nodes within the slab by the Y component, and carving
 *      the slice into rows along the Y axis such that each row
 *      in the X slice contain roughly the same weighting.
 */
        for (xDom = 0; xDom < numXdoms; xDom++) {

            nodesInXSlice = xSlices[xDom];

            yFirstIndex = xNodeIndex;
            yLastIndex  = xNodeIndex + nodesInXSlice;
            yNodeIndex  = yFirstIndex;

            xNodeIndex += nodesInXSlice;

            qsort(&nodeList[yFirstIndex], yLastIndex-yFirstIndex,
                  sizeof(NodeData_t), YCmp);

            GetBounds(cellList, &nodeList[yFirstIndex], Y,
                      yLastIndex-yFirstIndex, numYdoms, param->minSideY,
                      param->maxSideY, yBound, ySlices);

            decomp->domBoundX[xDom] = xBound[xDom];

/*
 *          Lastly loop through each row in the Y slice, resorting
 *          coordinates for all nodes within the row by the Z component,
 *          and carving the row into elements along the Z axis such that each
 *          element of the row contains roughly the same weighting.
 */
            for (yDom = 0; yDom < numYdoms; yDom++) {

                nodesInYSlice = ySlices[yDom];

                zFirstIndex = yNodeIndex;
                zLastIndex  = yNodeIndex + nodesInYSlice;
                zNodeIndex  = zFirstIndex;

                yNodeIndex += nodesInYSlice;

                qsort(&nodeList[zFirstIndex], zLastIndex-zFirstIndex,
                      sizeof(NodeData_t), ZCmp);

                GetBounds(cellList, &nodeList[zFirstIndex], Z,
                          zLastIndex-zFirstIndex, numZdoms, param->minSideZ,
                          param->maxSideZ, zBound, zSlices);

                decomp->domBoundY[xDom][yDom] = yBound[yDom];

/*
 *              Now just loop over the Z domains and set the boundaries
 *              for all domains in this row of the problem space.
 */
                for (zDom = 0; zDom < numZdoms; zDom++) {
                    decomp->domBoundZ[xDom][yDom][zDom] = zBound[zDom];
                }

                decomp->domBoundZ[xDom][yDom][numZdoms] = zBound[numZdoms];
            }  /* loop over Y domains */

            decomp->domBoundY[xDom][numYdoms] = yBound[numYdoms];
        }  /* loop over X domains */

        decomp->domBoundX[numXdoms] = xBound[numXdoms];


        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GenUGrid
 *      Description: Generate a uniform grid encompassing the problem
 *                   space.  Note: The caller is responsible for
 *                   freeing all memory buffers allocated by this
 *                   function.
 *
 *      Arguments:
 *          gridDim     3 element array specifying the number of sectors
 *                      per dimension to have in the uniform grid
 *          xBoundaries Location in which to return a pointer to the
 *                      array of coordinates defining the planes along
 *                      which the problems space is partitioned in the
 *                      X dimension.  Array contains <gridDim[X]+1> elements.
 *          yBoundaries Same description as <xBoundaries> for Y dimension.
 *          zBoundaries Same description as <xBoundaries> for Z dimension.
 *
 *------------------------------------------------------------------------*/
static void GenUGrid(Home_t *home, int gridDim[3], real8 **xBoundaries,
                     real8 **yBoundaries, real8 **zBoundaries)
{
        int   i;
        real8 xSectorLen, ySectorLen, zSectorLen;
        real8 xMinCoord, yMinCoord, zMinCoord;
        real8 *xBounds, *yBounds, *zBounds;

        xSectorLen = home->param->Lx / gridDim[X];
        ySectorLen = home->param->Ly / gridDim[Y];
        zSectorLen = home->param->Lz / gridDim[Z];

        xMinCoord = home->param->minSideX;
        yMinCoord = home->param->minSideY;
        zMinCoord = home->param->minSideZ;

/*
 *      Allocated arrays to store the coordinates of the planes
 *      along which each dimension of the problem space is to be
 *      uniformly sectioned. Then go ahead and calculate the plane
 *      coordinates.
 */
        xBounds = (real8 *)calloc(1, (gridDim[X]+1) * sizeof(real8));
        yBounds = (real8 *)calloc(1, (gridDim[Y]+1) * sizeof(real8));
        zBounds = (real8 *)calloc(1, (gridDim[Z]+1) * sizeof(real8));

        for (i = 0; i < gridDim[X]; i++) {
            xBounds[i] = xMinCoord + (i * xSectorLen);
        }

        for (i = 0; i < gridDim[Y]; i++) {
            yBounds[i] = yMinCoord + (i * ySectorLen);
        }

        for (i = 0; i < gridDim[Z]; i++) {
            zBounds[i] = zMinCoord + (i * zSectorLen);
        }

        xBounds[gridDim[X]] = home->param->maxSideX;
        yBounds[gridDim[Y]] = home->param->maxSideY;
        zBounds[gridDim[Z]] = home->param->maxSideZ;

        *xBoundaries = xBounds; 
        *yBoundaries = yBounds; 
        *zBoundaries = zBounds; 

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetRepartULoadData
 *      Description: From the cell weighting values, estimate the load
 *                   data for each sector of a uniform grid encompassing
 *                   the problem space. 
 *
 *                   NOTE: This functions assumes the number of sectors
 *                   in each dimension of the uniform grid matches the
 *                   number of cells in the corresponding dimension.
 *
 *      Arguments:
 *          numCells   Total number of cell structures in <cellList>
 *          cellList   Array of cell-specific load data
 *          uloadData  Location in which to return to caller the
 *                     array of load data values based on the 
 *                     cell weighting information.
 *
 *------------------------------------------------------------------------*/
static void GetRepartULoadData(int numCells, CellData_t *cellList, real8 **uloadData)
{
        int   cellIndex;
        real8 *cellLoadData;

        cellLoadData = (real8 *)calloc(1, numCells * sizeof(real8));

        for (cellIndex = 0; cellIndex < numCells; cellIndex++) {
            cellLoadData[cellIndex] = cellList[cellIndex].cellNodeCount *
                                      cellList[cellIndex].cellBlockWeight;
        }

        *uloadData = cellLoadData;

        return;
}


static void GetDecompCnt(int xNumDoms, int yNumDoms, int zNumDoms, int *totDecomps)
{
        int  i;
        int  xDoms[2], yDoms[2], zDoms[2];
        int  cut[3], subDecompUsed[8];

        *totDecomps += 1;

/*
 *      Figure out the dimensions (if any) in which this portion of the
 *      deocmposition must be further decomposed.
 */
        cut[X] = (xNumDoms > 1);
        cut[Y] = (yNumDoms > 1);
        cut[Z] = (zNumDoms > 1);

        if (cut[X] + cut[Y] + cut[Z] == 0) {
            return;
        }

        for (i = 0; i < 8; i++) {
            subDecompUsed[i] = 1;
        }
/*
 *      If we're not cutting in the X dimension, none of the
 *      right side octants will be used.  Also set the number of
 *      domains in the X dimension that will be encompassed by
 *      octants on each side of the X bisection (left and right
 *      side octants).
 *
 *      Then the same sort of thing for Y and Z dimensions.
 */
        if (!cut[X]) {
            subDecompUsed[LRF] = 0;
            subDecompUsed[URF] = 0;
            subDecompUsed[LRB] = 0;
            subDecompUsed[URB] = 0;
            xDoms[0] = xNumDoms;
        } else {
            xDoms[0] = xNumDoms >> 1;
            xDoms[1] = xNumDoms - xDoms[0];
        }

        if (!cut[Y]) {
            subDecompUsed[ULF] = 0;
            subDecompUsed[URF] = 0;
            subDecompUsed[ULB] = 0;
            subDecompUsed[URB] = 0;
            yDoms[0] = yNumDoms;
        } else {
            yDoms[0] = yNumDoms >> 1;
            yDoms[1] = yNumDoms - yDoms[0];
        }

        if (!cut[Z]) {
            subDecompUsed[LLB] = 0;
            subDecompUsed[LRB] = 0;
            subDecompUsed[ULB] = 0;
            subDecompUsed[URB] = 0;
            zDoms[0] = zNumDoms;
        } else {
            zDoms[0] = zNumDoms >> 1;
            zDoms[1] = zNumDoms - zDoms[0];
        }

/*
 *      Loop through all the octants that will be populated, incrementing
 *      the total decomp count with the number of portions into which the
 *      octant is sub-partitioned.
 */
        for (i = 0; i < 8; i++) {

            if (!subDecompUsed[i]) {
                continue;
            }

            switch (i) {
                case (LLF): { GetDecompCnt(xDoms[0], yDoms[0], zDoms[0], totDecomps); break; }
                case (LRF): { GetDecompCnt(xDoms[1], yDoms[0], zDoms[0], totDecomps); break; }
                case (ULF): { GetDecompCnt(xDoms[0], yDoms[1], zDoms[0], totDecomps); break; }
                case (URF): { GetDecompCnt(xDoms[1], yDoms[1], zDoms[0], totDecomps); break; }
                case (LLB): { GetDecompCnt(xDoms[0], yDoms[0], zDoms[1], totDecomps); break; }
                case (LRB): { GetDecompCnt(xDoms[1], yDoms[0], zDoms[1], totDecomps); break; }
                case (ULB): { GetDecompCnt(xDoms[0], yDoms[1], zDoms[1], totDecomps); break; }
                case (URB): { GetDecompCnt(xDoms[1], yDoms[1], zDoms[1], totDecomps); break; }
            }
        }
        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:    GetGridOverlap
 *      Description: Given a set of boundaries <cMin> and <cMax>
 *                   determine the minimum and maximum indices
 *                   in each dimension of the elements of the
 *                   uniform grid defined by <[xyz]Boundaries> and
 *                   return these indices to the caller along with
 *                   ratio of overlap of each slab (in each dimension)
 *                   with the specified block.
 *
 *      Arguments:
 *          cMin        3 element array containing minimum coordinates
 *                      in all three dimensions of the volume of space
 *                      associated with the load data <currLoad>
 *          cMax        3 element array containing maximum coordinates
 *                      in all three dimensions of the volume of space
 *                      associated with the load data <currLoad>
 *          gridDim     3 element array specifying the number of sectors
 *                      per dimension in the uniform load data grid 
 *          xBoundaries Array of coordinates defining the planes along
 *                      which the problems space is partitioned in the
 *                      X dimension.  Array contains <gridDim[X]+1> elements.
 *          yBoundaries Same description as <xBoundaries> for Y dimension.
 *          zBoundaries Same description as <xBoundaries> for Z dimension.
 *          xMinIndex   Minimum index in X of sectors of the uniform
 *                      grid intersecting the specified block
 *          yMinIndex   Same as xMinIndex for Y dimension
 *          zMinIndex   Same as xMinIndex for Z dimension
 *          xMaxIndex   Maximum index in X of sectors of the uniform
 *                      grid intersecting the specified block
 *          yMaxIndex   Same as xMaxIndex for Y dimension
 *          zMaxIndex   Same as xMaxIndex for Y dimension
 *          xOverlap    Array of values indicating what portion of the
 *                      slab of the uniform grid intersects the caller
 *                      specified volume.  Array contains 1 element for
 *                      each X slab.
 *          yOverlap    Same as xOverlap for Y dimension
 *          yOverlap    Same as xOverlap for Z dimension
 *
 *------------------------------------------------------------------------*/
static void GetGridOverlap(real8 *cMin, real8 *cMax, int gridDim[3],
                           real8 *xBoundaries, real8 *yBoundaries,
                           real8 *zBoundaries,
                           int *xMinIndex, int *yMinIndex, int *zMinIndex,
                           int *xMaxIndex, int *yMaxIndex, int *zMaxIndex,
                           real8 **xOverlap, real8 **yOverlap, real8 **zOverlap)
{
        int   i;
        int   xMin, yMin, zMin;
        int   xMax, yMax, zMax;
        real8 *xFactor, *yFactor, *zFactor, overLap;

        xFactor = (real8 *)calloc(1, gridDim[X] * sizeof(real8));
        yFactor = (real8 *)calloc(1, gridDim[Y] * sizeof(real8));
        zFactor = (real8 *)calloc(1, gridDim[Z] * sizeof(real8));

        xMin = 0;
        yMin = 0;
        zMin = 0;

        xMax = gridDim[X] - 1;
        yMax = gridDim[Y] - 1;
        zMax = gridDim[Z] - 1;

/*
 *      Find the minimum and maximum indices (in each dimension)
 *      of the group of sectors of the uniform grid that intersect
 *      the volume defined by the minimum and maximum coordinates
 *      in <cMin> and <cMax>.
 */
        for (i = 0; i < gridDim[X]; i++) {
            if (cMin[X] >= xBoundaries[i+1]) {
                xMin = i + 1;
            } else if (cMax[X] <= xBoundaries[i+1]) {
                xMax = i;
                break;
            }
        }

        for (i = 0; i < gridDim[Y]; i++) {
            if (cMin[Y] >= yBoundaries[i+1]) {
                yMin = i + 1;
            } else if (cMax[Y] <= yBoundaries[i+1]) {
                yMax = i;
                break;
            }
        }

        for (i = 0; i < gridDim[Z]; i++) {
            if (cMin[Z] >= zBoundaries[i+1]) {
                zMin = i + 1;
            } else if (cMax[Z] <= zBoundaries[i+1]) {
                zMax = i;
                break;
            }
        }

/*
 *      For each slab of the uniform grid (in each dimension)
 *      that intersects the specified volume, determine what
 *      portion of the volume overlaps the slab.
 */
        for (i = xMin; i <= xMax; i++) {
            overLap = MIN(cMax[X], xBoundaries[i+1]) -
                      MAX(cMin[X], xBoundaries[i]);
            xFactor[i] = overLap / (cMax[X] - cMin[X]);
        }

        for (i = yMin; i <= yMax; i++) {
            overLap = MIN(cMax[Y], yBoundaries[i+1]) -
                      MAX(cMin[Y], yBoundaries[i]);
            yFactor[i] = overLap / (cMax[Y] - cMin[Y]);
        }

        for (i = zMin; i <= zMax; i++) {
            overLap = MIN(cMax[Z], zBoundaries[i+1]) -
                      MAX(cMin[Z], zBoundaries[i]);
            zFactor[i] = overLap / (cMax[Z] - cMin[Z]);
        }

        *xOverlap = xFactor;
        *yOverlap = yFactor;
        *zOverlap = zFactor;

        *xMinIndex = xMin;
        *yMinIndex = yMin;
        *zMinIndex = zMin;

        *xMaxIndex = xMax;
        *yMaxIndex = yMax;
        *zMaxIndex = zMax;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       GetBisectionCoord
 *      Description:    Calculates the coordinate of the plane (in
 *                      dimension <dim> that bisects the volume 
 *                      given by <cMin> and <cMax> such that 1/2
 *                      the computational load for the volume will
 *                      fall on each side of the plane.
 *
 *      Arguments:
 *          cMin        3 element array containing minimum coordinates
 *                      in all three dimensions of the volume of space
 *                      to be bisected.
 *          cMax        3 element array containing maximum coordinates
 *                      in all three dimensions of the volume of space
 *                      to be bisected.
 *          dim         Dimension in which to bisect the specified volume.
 *                      0 == X, 1 == Y, 2 == Z.
 *          gridDim     3 element array specifying the number of sectors
 *                      per dimension in the uniform load data grid
 *          uloadData   Array of load/weight values for each sector of the
 *                      uniform grid of dimensions <gridDim> overlaid on the
 *                      problem space.
 *          xBoundaries Array of coordinates defining the planes along
 *                      which the problems space is partitioned in the
 *                      X dimension.  Array contains <gridDim[X]+1> elements.
 *          yBoundaries Same description as <xBoundaries> for Y dimension.
 *          zBoundaries Same description as <xBoundaries> for Z dimension.
 *          coord       Location in which to return the coordinate of the
 *                      bisecting plane.
 *          minSideRatio  portion of load to be accumulated on the minimum
 *                        side of the bisection coordinate.  Must be
 *                        between 0.0 and 1.0.
 * 
 *-----------------------------------------------------------------------*/
static void GetBisectionCoord(real8 *cMin, real8 *cMax, int dim,
                              int gridDim[3], real8 *uloadData,
                              real8 *xBoundaries, real8 *yBoundaries,
                              real8 *zBoundaries, real8 *coord,
                              real8 minSideRatio)
{
        int i, offset, numSlabs;
        int slabMinIndex=0, slabMaxIndex=0;
        int xMin, yMin, zMin;
        int xMax, yMax, zMax;
        int index[3]; 
        real8 totalLoad, neededLoad, overlapFactor;
        real8 *slabLoad=0, *boundaries=0;
        real8 *xOverlap=0, *yOverlap=0, *zOverlap=0;

/*
 *      Find the min and max indices (per dimension) of the sectors
 *      of the uniform grid overlapped by the specified volume.
 *      Also return the portion overlap of the volume with the slabs of
 *      the grid (in each dimension)
 */
        GetGridOverlap(cMin, cMax, gridDim, xBoundaries, yBoundaries,
                       zBoundaries, &xMin, &yMin, &zMin, &xMax, &yMax,
                       &zMax, &xOverlap, &yOverlap, &zOverlap);

        numSlabs = gridDim[dim];
        slabLoad = (real8 *)calloc(1, numSlabs * sizeof(real8));

        switch(dim) {
            case X:
                slabMinIndex = xMin;
                slabMaxIndex = xMax;
                boundaries = xBoundaries;
                break;
            case Y:
                slabMinIndex = yMin;
                slabMaxIndex = yMax;
                boundaries = yBoundaries;
                break;
            case Z:
                slabMinIndex = zMin;
                slabMaxIndex = zMax;
                boundaries = zBoundaries;
                break;
        }

/*
 *      Find the total weight/load for the specified volume, plus
 *      find the weight of each slab of the uniform grid (in the
 *      specified dimension) that interescts the specified volume.
 *      Note: slab weight is only for the portion of the slab
 *      that intersects the specified volume.
 */
        totalLoad = 0.0;

        for (index[X] = xMin; index[X] <= xMax; index[X]++) {
            for (index[Y] = yMin; index[Y] <= yMax; index[Y]++) {
                for (index[Z] = zMin; index[Z] <= zMax; index[Z]++) {

                    offset = (index[X] * gridDim[Y] * gridDim[Z]) +
                             (index[Y] * gridDim[Z]) + index[Z];

                    overlapFactor = xOverlap[index[X]] *
                                    yOverlap[index[Y]] *
                                    zOverlap[index[Z]];

                    slabLoad[index[dim]] += uloadData[offset] * overlapFactor;
                    totalLoad += uloadData[offset] * overlapFactor;
                }
            }
        }

/*
 *      We have the load-weight (assumed evenly distributed) for each slab
 *      in the specified dimension, the coordinates of each slab, and
 *      the coordinates of the volume we want bisected.  Now find the
 *      coordinate of the plane (in the specified dimension) that would
 *      result in the ratio (of load on the minimum side of the bisection
 *      plane to the total load) of <minSideRatio>.
 */
        neededLoad = totalLoad * minSideRatio;
        *coord = MAX(boundaries[slabMinIndex], cMin[dim]);

/*
 *      Just a quick sanity check...
 */
        if (neededLoad == 0.0) {
            *coord = cMin[dim] + 0.5 * (cMax[dim] - cMin[dim]);
            return;
        }

        for (i = slabMinIndex; i <= slabMaxIndex; i++) {
            if (neededLoad >= slabLoad[i]) {
                neededLoad -= slabLoad[i];
                *coord = boundaries[i+1];
            } else {
                *coord += (boundaries[i+1] - *coord) *
                           (neededLoad / slabLoad[i]);
                break;
            }
        }

        free(xOverlap);
        free(yOverlap);
        free(zOverlap);
        free(slabLoad);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       RepartRBDecomp
 *      Description:    Calculate the boundaries for the branch of
 *                      the decomposition tree terminating at the
 *                      local domain's decomposition based on the
 *                      current load (computation cost) distribution among
 *                      the processors.
 *
 *                      Once each domain has calculated its own boundaries
 *                      A call to ExchangeRBDecomp() can be used to
 *                      share those boundaries with all other processors
 *                      so all processors have the full decomposition.
 *
 *      Arguments:
 *          decomp       pointer to new decomposition tree which
 *                       has been pre-allocated and partially initialized
 *                       via AllocRBDecomp() and InitRBDecomp().
 *          oldDecomp    pointer to previous decomposition tree.
 *          domDecompID  ID of portion of decomposition associated with
 *                       the current domain
 *          currLevel    current level in the domain decomposition.
 *          startLevel   Level at which to start the decomposition.  The
 *                       caller may wish to do only a partial decomposition
 *                       in which case some number of levels of the
 *                       decompositiontree will simply have been copied
 *                       from a previous decomposition.  All levels less
 *                       than <startLevel> will simply be left as they are
 *          uloadData    Array of load/weight values for each sector of the
 *                       uniform grid of dimensions <gridDim> overlaid on the
 *                       problem space.
 *          gridDim      3 element array specifying the number of sectors
 *                       per dimension in the uloadData grid
 *          xBoundaries  Array of coordinates defining the planes along
 *                       which the load data grid is partitioned in the
 *                       X dimension.  Array contains <gridDim[X]+1> elements.
 *          yBoundaries  Same description as <xBoundaries> for Y dimension.
 *          zBoundaries  Same description as <xBoundaries> for Z dimension.
 *
 *-------------------------------------------------------------------------*/
static void RepartRBDecomp(Home_t *home, RBDecomp_t *decomp,
                           RBDecomp_t *oldDecomp, char *domDecompID,
                           int currLevel, int startLevel, real8 *uloadData,
                           int gridDim[3], real8 *xBoundaries,
                           real8 *yBoundaries, real8 *zBoundaries)
{
        int          i, octant;
        int          cut[3];
        real8        coord;
        real8        minRatio[3];
        Param_t      *param;
        RBDecomp_t   *subDecomp, **subDecompList;

        param = home->param;


        if (decomp == (RBDecomp_t *)NULL) {
            return;
        }

/*
 *      Figure out the dimensions in which this partition is to be
 *      sub-partitioned.  Sub-partitioning will be into octants, quadrants
 *      halves or no subpartitions at all.
 */
        cut[X] = (decomp->numDoms[X] > 1);
        cut[Y] = (decomp->numDoms[Y] > 1);
        cut[Z] = (decomp->numDoms[Z] > 1);

/*
 *      If there is no need to do further decomposition, we're
 *      down to an individual domain, so we can return to the caller.
 */
        if ((cut[X] + cut[Y] + cut[Z]) == 0) {
            return;
        }

/*
 *      When we cut the partion in a dimension, the two 'halves' may
 *      not be asigned the same number of procesors, so we need to
 *      calculate the portion of the load that will be associated with
 *      the lower side of the bisection plane in each dimension.
 */
        minRatio[X] = (real8)(decomp->numDoms[X] / 2) /
                      (real8)decomp->numDoms[X];
        minRatio[Y] = (real8)(decomp->numDoms[Y] / 2) /
                      (real8)decomp->numDoms[Y];
        minRatio[Z] = (real8)(decomp->numDoms[Z] / 2) /
                      (real8)decomp->numDoms[Z];

/*
 *      Set up the initial boundaries on all used subpartitions so
 *      they are identical to the current pratition boundaries.
 *      As we determine the proper planes along which to subpartition
 *      the current partition, the subpartition boundaries will be 
 *      updated appropriately.
 */
        for (i = 0; i < 8; i++) {

            subDecomp = decomp->subDecomp[i];

            if (subDecomp != (RBDecomp_t *)NULL) {
    
                subDecomp->cMin[X] = decomp->cMin[X];
                subDecomp->cMin[Y] = decomp->cMin[Y];
                subDecomp->cMin[Z] = decomp->cMin[Z];

                subDecomp->cMax[X] = decomp->cMax[X];
                subDecomp->cMax[Y] = decomp->cMax[Y];
                subDecomp->cMax[Z] = decomp->cMax[Z];
            }
        }

/*
 *      If the partition needs to be cut in the X dimension
 *      Find the coordinate of the plane along which the 
 *      partition needs to be cut in the X dimension and
 *      update the min/max coordinates for each subpartition
 *      that could be affected by slicing the partition in
 *      this dimension.
 */
        subDecompList = decomp->subDecomp;

        if (cut[X]) {

            GetBisectionCoord(subDecompList[LLF]->cMin,
                              subDecompList[LLF]->cMax, X, gridDim,
                              uloadData, xBoundaries, yBoundaries,
                              zBoundaries, &coord, minRatio[X]);

/*
 *          If the current partition is not being subpartitioned
 *          in all 3 dimensions, not all of the octants will be
 *          have been allocated... We know certain octants must
 *          exist, but we have to verify the existence of the others
 *          before we try accessing them.  The same situation also
 *          applies further on when we cut the Y and Z dimension.
 */
            subDecompList[LLF]->cMax[X] = coord;
            subDecompList[LRF]->cMin[X] = coord;

            if (subDecompList[ULF]) subDecompList[ULF]->cMax[X] = coord;
            if (subDecompList[URF]) subDecompList[URF]->cMin[X] = coord;
            if (subDecompList[LLB]) subDecompList[LLB]->cMax[X] = coord;
            if (subDecompList[LRB]) subDecompList[LRB]->cMin[X] = coord;
            if (subDecompList[ULB]) subDecompList[ULB]->cMax[X] = coord;
            if (subDecompList[URB]) subDecompList[URB]->cMin[X] = coord;
        }

/*
 *      If the partition needs to be cut in the Y dimension
 *      determine the Y-coordinates of the planes along which
 *      to slice the (one or more) subpartitions formed after
 *      the partition has been sliced (possibly) in the X dimension,
 *      then update the min/max coordinates for each subpartition
 *      that could be affected by slicing the partition in this
 *      dimension.
 */
        if (cut[Y]) {

/*
 *          Cut all the existing LLF, LRF halves in the Y dimension.
 */
            GetBisectionCoord(subDecompList[LLF]->cMin,
                              subDecompList[LLF]->cMax, Y, gridDim,
                              uloadData, xBoundaries, yBoundaries,
                              zBoundaries, &coord, minRatio[Y]);

            subDecompList[LLF]->cMax[Y] = coord;
            subDecompList[ULF]->cMin[Y] = coord;

            if (subDecompList[LLB]) subDecompList[LLB]->cMax[Y] = coord;
            if (subDecompList[ULB]) subDecompList[ULB]->cMin[Y] = coord;

            if (subDecompList[LRF]) {
                GetBisectionCoord(subDecompList[LRF]->cMin,
                                  subDecompList[LRF]->cMax, Y, gridDim,
                                  uloadData, xBoundaries, yBoundaries,
                                  zBoundaries, &coord, minRatio[Y]);

                subDecompList[LRF]->cMax[Y] = coord;
                subDecompList[URF]->cMin[Y] = coord;

                if (subDecompList[LRB]) subDecompList[LRB]->cMax[Y] = coord;
                if (subDecompList[URB]) subDecompList[URB]->cMin[Y] = coord;
            }

        }  /* if (cut[Y]) */

/*
 *      If the partition needs to be cut in the Z dimension
 *      determine the z-coordinates of the planes along which
 *      to slice the (one or more) subpartitions formed after
 *      the partition has been sliced (possibly) in the X and
 *      Y dimensions, then update the min/max coordinates for
 *      each subpartition that could be affected by slicing
 *      the partition in this dimension.
 */
        if (cut[Z]) {

/*
 *          Cut all the existing LLF, LRF, ULF and URF quadrants
 *          in the Z dimension.
 */
            GetBisectionCoord(subDecompList[LLF]->cMin,
                              subDecompList[LLF]->cMax, Z, gridDim,
                              uloadData, xBoundaries, yBoundaries,
                              zBoundaries, &coord, minRatio[Z]);

            subDecompList[LLF]->cMax[Z] = coord;
            subDecompList[LLB]->cMin[Z] = coord;

            if (subDecompList[LRF]) {
                GetBisectionCoord(subDecompList[LRF]->cMin,
                                  subDecompList[LRF]->cMax, Z, gridDim,
                                  uloadData, xBoundaries, yBoundaries,
                                  zBoundaries, &coord, minRatio[Z]);

                subDecompList[LRF]->cMax[Z] = coord;
                subDecompList[LRB]->cMin[Z] = coord;
            }

            if (subDecompList[ULF]) {
                GetBisectionCoord(subDecompList[ULF]->cMin,
                                  subDecompList[ULF]->cMax, Z, gridDim,
                                  uloadData, xBoundaries, yBoundaries,
                                  zBoundaries, &coord, minRatio[Z]);

                subDecompList[ULF]->cMax[Z] = coord;
                subDecompList[ULB]->cMin[Z] = coord;
            }

            if (subDecompList[URF]) {
                GetBisectionCoord(subDecompList[URF]->cMin,
                                  subDecompList[URF]->cMax, Z, gridDim,
                                  uloadData, xBoundaries, yBoundaries,
                                  zBoundaries, &coord, minRatio[Z]);
    
                subDecompList[URF]->cMax[Z] = coord;
                subDecompList[URB]->cMin[Z] = coord;
            }
    
        }  /* if (cut[Z]) */

/*
 *      Boundaries are now properly set for all necessary octants of
 *      the current portion of the decomposition.  Now continue with
 *      a full decomposition of all subpartitions.
 */
        for (octant = 0; octant < 8; octant++) {
            RepartRBDecomp(home, decomp->subDecomp[octant], oldDecomp,
                           domDecompID, currLevel+1, startLevel, uloadData,
                           gridDim, xBoundaries, yBoundaries, zBoundaries);
        }


        return;
}


int main (int argc, char *argv[])
{
        int     numXdoms, numYdoms, numZdoms;
        int     numXcells, numYcells, numZcells, numCells;
        int     gridDim[3];
        real8   *xBoundaries, *yBoundaries, *zBoundaries, *uloadData;
        char    backupFileName[256], textLine[512];
        char    oldFileName[256], newFileName[256];
        char    *baseFileName, *outFileName;
        fpos_t  nodeDataPos;
        FILE    *oldFP, *newFP;
        NodeData_t *nodeList;
        CellData_t *cellList;
        Param_t *param;
        Home_t  home;

/*
 *      Get the input arguments and do general initialization
 */
        param = (Param_t *)malloc(sizeof(Param_t));
        param->xBoundType = Periodic;
        param->yBoundType = Periodic;
        param->zBoundType = Periodic;

        home.param = param;
        home.dataParamList = (ParamList_t *)calloc(1, sizeof(ParamList_t));

        DataParamInit(param, home.dataParamList);

        GetInArgs(argc, argv, param, &baseFileName, &outFileName);

        numXcells = param->nXcells;
        numYcells = param->nYcells;
        numZcells = param->nZcells;

        numXdoms = param->nXdoms;
        numYdoms = param->nYdoms;
        numZdoms = param->nZdoms;

        home.numDomains = numXdoms * numYdoms * numZdoms;

/*
 *      And a couple quick sanity checks on the input args...
 */
        if ((numXdoms < 1) || (numYdoms < 1) || (numZdoms < 1)) {
            fprintf(stderr, "Domain geometry (%d X %d X %d) invalid.\n",
                    numXdoms, numYdoms, numZdoms);
            exit(1);
        }

        if ((numXcells < 3) || (numYcells < 3) || (numZcells < 3)) {
            fprintf(stderr, "Cell geometry (%d X %d X %d) invalid.\n",
                    numXcells, numYcells, numZcells);
            exit(1);
        }

        if (baseFileName == (char *)NULL) {
            fprintf(stderr, "No input file specified.\n");
            exit(1);
        }

/*
 *      For RB decompositions we need to initialize some stuff...
 */
        if (param->decompType == 2) {
            for (home.xMaxLevel = 0; param->nXdoms >> home.xMaxLevel > 1;
                 home.xMaxLevel++);

            for (home.yMaxLevel = 0; param->nYdoms >> home.yMaxLevel > 1;
                 home.yMaxLevel++);

            for (home.zMaxLevel = 0; param->nZdoms >> home.zMaxLevel > 1;
                 home.zMaxLevel++);

            if ((1 << home.xMaxLevel) < param->nXdoms) home.xMaxLevel++;
            if ((1 << home.yMaxLevel) < param->nYdoms) home.yMaxLevel++;
            if ((1 << home.zMaxLevel) < param->nZdoms) home.zMaxLevel++;

        }

/*
 *      Allocate array of structures for cell-specific data
 */
        numCells  = numXcells * numYcells * numZcells;

        cellList = (CellData_t *)calloc(1, numCells * sizeof(CellData_t));

/*
 *      Read the nodal coordinates and other ancillary data from the
 *      data file.
 */
        ReadDataFile(&home, baseFileName, cellList, &nodeList, &nodeDataPos);

#if 0
/*
 *      The repartitioning utility will not convert nodal data from
 *      one format to another.  Therefore, if the nodal data in the
 *      input file is not of the latest format, have the user convert
 *      the nodal data file to the current version before trying to
 *      repartition the data.
 *
 *      This is stubbed out for now since all versions of the data
 *      file use the same format for the nodal data.
 */
        if (param->dataFileVersion != NODEDATA_FILE_VERSION) {
            Fatal("The nodal data in the specified input file is\n"
                  "in an older format.  Please convert the input file\n"
                  "to the current format using 'paradisconvert' and then\n"
                  "try the repartitioning utility again.");
        }
#endif
        
/*
 *      Overwrite the domain geometry from the restart file with
 *      the user provided domain geometry.
 */
        param->dataDecompGeometry[X] = numXdoms;
        param->dataDecompGeometry[Y] = numYdoms;
        param->dataDecompGeometry[Z] = numZdoms;

/*
 *      Calculate the total weighting for the problem, and
 *      calculate the cell and cellblock weight for each cell
 */
        GetWeighting(param, cellList);

/*
 *      Call the decomposition function appropriate to the user
 *      specified decomposition type.
 */
        switch (param->decompType) {
        case 1:
            RepartRSDecomp(&home, nodeList, cellList);
            break;
        case 2:
/*
 *          The type 2 decomposition will be based on the estimated
 *          load for each sector of a uniform grid encompassing the
 *          entire problem space.  Before doing the decomposition
 *          we need to generate the grid and set the load data for
 *          each sector.
 */
            gridDim[X] = numXcells;
            gridDim[Y] = numYcells;
            gridDim[Z] = numZcells;

            GenUGrid(&home, gridDim, &xBoundaries, &yBoundaries, &zBoundaries);
            GetRepartULoadData(numCells, cellList, &uloadData);

/*
 *          Pre-allocate the decomposition tree structure then do the
 *          domain decomposition based on the estimated load
 */
            AllocRBDecomp(&home, (RBDecomp_t *)NULL,
                          (RBDecomp_t **)&home.decomp, ALLOC_NEW_DECOMP);

            RepartRBDecomp(&home, (RBDecomp_t *)home.decomp, (RBDecomp_t *)NULL,
                           (char *)NULL, 0, 0, uloadData, gridDim, xBoundaries,
                           yBoundaries, zBoundaries);
            break;
        default:
            printf("Error: Unsupported domain decomposition type %d\n",
                   param->decompType);
            exit(-1);
        }

/*
 *      Now we need to reopen the first segment of the original
 *      data file in order to copy out the nodal data, plus create
 *      a new file which will be a copy of the old file expect that
 *      the domain decomposition will be updated.
 */
        strncpy(oldFileName, baseFileName, sizeof(oldFileName));

        if ((oldFP = fopen(oldFileName, "r")) == (FILE *)NULL) {
            strcat(oldFileName, ".0");
            if ((oldFP = fopen(oldFileName, "r")) == (FILE *)NULL) {
                printf("Error %d opening file %s\n", errno, baseFileName);
                exit(1);
            }
        }

        snprintf(backupFileName, sizeof(backupFileName), "%s.bkup",
                 oldFileName);
        sprintf(newFileName, "%s.decomptmp", oldFileName);

/*
 *      Open the new data file and write the version number,
 *      domain boundaries, node count, etc.  Make sure the
 *      file is rewritten in the format applicable to the
 *      the version number of the original file.
 */
        newFP = fopen(newFileName, "w");

        if ((newFP = fopen(newFileName, "w")) == (FILE *)NULL) {
            fprintf(stderr, "fopen error %d for file %s\n",
                    errno, newFileName);
            exit(1);
        }

/*
 *      The format of the new nodal data file is forced to the
 *      default version supported by this release of the code
 *      regardless of the version of the user specified data file.
 */
        param->dataFileVersion = NODEDATA_FILE_VERSION;
        param->dataDecompType = param->decompType;
        

        WriteParam(home.dataParamList, -1, newFP);

/*
 *      Write the domain decomposition into nodal data file
 *      and then some comment lines describing the nodal
 *      data that will follow.
 */
        fprintf(newFP, "\n#\n#  END OF DATA FILE PARAMETERS\n#\n\n");
        fprintf(newFP, "domainDecomposition = \n");
        WriteDecompBounds(&home, newFP);

        fprintf(newFP, "nodalData = \n");

/*
 *      Last thing we need to do is copy the nodal data from the
 *      original file to the new file.  Reset the original file's
 *      file pointer to the location of the nodal information,
 *      read each remaining line from the original file and dump
 *      it into the new file.
 */
        fsetpos(oldFP, &nodeDataPos);

        while (fgets(textLine, sizeof(textLine)-1, oldFP) != NULL) {
            fprintf(newFP, "%s", textLine);
        }

        fclose(oldFP);
        fclose(newFP);

        if (outFileName == (char *)NULL) {
            outFileName = oldFileName;
        }

        rename(oldFileName, backupFileName);
        rename(newFileName, outFileName);

        printf("\n  Original file:  %s\n", backupFileName);
        printf("  New file:       %s\n\n", outFileName);

        exit(0);
        return(0);
}
