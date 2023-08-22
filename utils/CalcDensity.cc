/*****************************************************************************
 *
 *      Module:         CalcDensity.c
 *      Description:    This module contains functions needed to read
 *                      the specified restart file, calculate the
 *                      dislocation density at the center of each
 *                      cell of a 'grid' of the size defined by the command
 *                      line args, and write the data to a file that
 *                      can be visualized via some external tool.
 *
 *      Usage:  CalcDensity -b -g <gridx[,y,z]> -h -d <dataFile> \
 *                          -o <outFile> <ctrlFile>
 *
 *      Includes functions: 
 *          DeleteUnreferencedNodes()
 *          FindDensity()
 *          FindNode()
 *          FindNodeIndex()
 *          ProcessBinNodalData()
 *          ProcessNodalData()
 *          ReadDataParams()
 *          ReadBinDataParams()
 *          ReadDataParams()
 *          ReallocNodeInfo()
 *          ReallocSegList()
 *          SetDensity()
 *          Usage()
 *
 *      Output file format:
 *
 ****************************************************************************/
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "Home.h"
#include "Decomp.h"
#include "Restart.h"

#ifdef USE_HDF
#include "hdf5.h"
#endif

#define SURF_NONE  -1
#define SURF_NEG_X  0
#define SURF_POS_X  1
#define SURF_NEG_Y  2
#define SURF_POS_Y  3
#define SURF_NEG_Z  4
#define SURF_POS_Z  5

#define DENSITY_FILE_SUFFIX ".dens"

#define SEG_BLOCK_INC       50000
#define NODE_BLOCK_INC      50000
#define NODE_READ_BLOCK_CNT 50000

/*
 *      Define a structure to hold the tags of the two endpoints uniquely
 *      describing a segment.
 */
typedef struct {
    Tag_t tag1;
    Tag_t tag2;
} Seg_t;

/*
 *      Define a structure to hold just the nodal data we need for this
 *      utility.
 */
typedef struct {
    Tag_t tag;
    int   refCnt;
    real8 pos[3];
} NodeInfo_t;


/*---------------------------------------------------------------------------
 *
 *      Function:    SetDensity
 *      Description: Given a segment from <startPos> to <endPos>, determines
 *                   how much of the segment is contained within the density
 *                   grid cell specified by <cellIndex>, updates the
 *                   cell's density value appropriately, and (if necesary)
 *                   recursively calls SetDensity() with the portion
 *                   of the segment that was outside the cell.
 *
 *      Arguments:
 *          densityGrid  3-D density grid represented as a single
 *                       dimension array.
 *          gridSize     specifies the size of the density grid in
 *                       X, Y and Z dimensions.
 *          cellSize     length (units of b) of density grid cells
 *                       in X, Y and Z.
 *          startPos     Coordinates of the starting position of the
 *                       segment.
 *          endPos       Coordinates of the ending position of the
 *                       segment.
 *          lineDir      Indicates the direction (-1 or +1) from <startPos>
 *                       to <endPos> in X, Y and Z.  
 *          cellIndex    Coordinates of the density grid cell to update.
 *                       Note:  This *could* be outside the range specified
 *                       in <gridSize>, but will be adjusted to handle pbc.
 *          burgVolFactor  factor by which segment length is adjusted to
 *                       calculate density.
 *
 *-------------------------------------------------------------------------*/
static void SetDensity(Home_t *home, real8 *densGrid, int gridSize[3],
                       real8 cellSize[3], real8 startPos[3], real8 endPos[3],
                       real8 lineDir[3], int cellIndex[3], real8 burgVolFactor)
{
        int     i, offset;
        int     intersectSurface, intersectIndex;
        int     trueIndex[3], newIndex[3];
        real8   xLen, yLen, zLen, rho;
        real8   intersectPlane;
        real8   lenInDim, lenPartial, lenFrac, minFrac;
        real8   intersectPos[3];
        real8   minCellBound[3], maxCellBound[3];
        Param_t *param;


        param = home->param;

        for (i = 0; i < 3; i++) {
            trueIndex[i] = cellIndex[i];
            newIndex[i] = cellIndex[i];
        }

/*
 *      We know the indices for the density grid cell, so calculate
 *      the boundaries for the cell.
 */
        for (i = 0; i < 3; i++) {
            minCellBound[i] = (cellIndex[i] * cellSize[i]) +
                              param->minCoordinates[i];
            maxCellBound[i] = ((cellIndex[i]+1) * cellSize[i]) +
                              param->minCoordinates[i];
        }

/*
 *      Look for the intersection point of the segment and the
 *      cell boundaries (only need to check boundaries in the
 *      same direction as the lineDir).
 */
        intersectIndex = -1;
        intersectPlane = 0.0;
        intersectSurface = SURF_NONE;

        minFrac = 1.0;

        for (i = 0; i < 3; i++) {
            lenInDim = endPos[i] - startPos[i];
            if (lineDir[i] < 0) {
                if (endPos[i] < minCellBound[i]) {
                    lenPartial = minCellBound[i] - startPos[i];
                    lenFrac = lenPartial / lenInDim;
                    if (lenFrac < minFrac) {
                       intersectPlane = minCellBound[i];
                       intersectIndex = i;
                       intersectSurface = i;
                       minFrac = lenFrac;
                    }
                }
            } else {
                if (endPos[i] > maxCellBound[i]) {
                    lenPartial = maxCellBound[i] - startPos[i];
                    lenFrac = lenPartial / lenInDim;
                    if (lenFrac < minFrac) {
                       intersectPlane = maxCellBound[i];
                       intersectIndex = i;
                       intersectSurface = i+3;
                       minFrac = lenFrac;
                    }
                }
            }
        }

/*
 *      If the segment pierces the cell boundaries, calculate the
 *      point at which it occurs.
 */
        for (i = 0; i < 3; i++) {
            if (i == intersectIndex) {
                intersectPos[i] = intersectPlane;
            } else {
                intersectPos[i] = startPos[i] +
                                  (minFrac * (endPos[i]-startPos[i]));
            }
        }

/*
 *      Calculate rho for segment from <startPos> to <intersectPos>
 */
        xLen = startPos[X] - intersectPos[X];
        yLen = startPos[Y] - intersectPos[Y];
        zLen = startPos[Z] - intersectPos[Z];


        rho = sqrt(xLen*xLen + yLen*yLen + zLen*zLen);
        rho *= burgVolFactor;

/*
 *      Get cell indices adjusted for PBC (if necessary).  If periodic
 *      boundaries are not used and the cell indices are outside the
 *      primary image, use the closest density field cell.
 */ 
        
        if (trueIndex[X] < 0) {
            trueIndex[X] = (gridSize[X]-1) * (param->xBoundType == Periodic);
        } else if (trueIndex[X] >= gridSize[X]) {
            trueIndex[X] = (gridSize[X]-1) * (param->xBoundType != Periodic);
        }

        if (trueIndex[Y] < 0) {
            trueIndex[Y] = (gridSize[Y]-1) * (param->yBoundType == Periodic);
        } else if (trueIndex[Y] >= gridSize[Y]) {
            trueIndex[Y] = (gridSize[Y]-1) * (param->yBoundType != Periodic);
        }

        if (trueIndex[Z] < 0) {
            trueIndex[Z] = (gridSize[Z]-1) * (param->zBoundType == Periodic);
        } else if (trueIndex[Z] >= gridSize[Z]) {
            trueIndex[Z] = (gridSize[Z]-1) * (param->zBoundType != Periodic);
        }

/*
 *      Increment rho for the cell containing the segment
 */
        offset = trueIndex[X] * gridSize[Y] * gridSize[Z] +
                 trueIndex[Y] * gridSize[Z] + trueIndex[Z];

        densGrid[offset] += rho;

/*
 *      If the segment was not fully contained within this cell (i.e.
 *      it intersected one of the cell boundaries), recursively calculate
 *      rho for the remaining portion of the segment that is outside
 *      the current cell.  (set <newIndex> to the indices of the cell
 *      containing <intersectPos>)
 */
        if (intersectSurface != SURF_NONE) {
/*
 *          A surface ID > 2 indicates the segment pierced the cell
 *          in the positive direction, surface IDs less than 3
 *          represent the negative surfaces.
 */
            if (intersectSurface > 2) {
                newIndex[intersectSurface-3]++;
            } else {
                newIndex[intersectSurface]--;
            }

            SetDensity(home, densGrid, gridSize, cellSize, intersectPos,
                       endPos, lineDir, newIndex, burgVolFactor);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    FindNode
 *      Description: Comparison function (suitable for use with bsearch)
 *                   that compares tag <a> with the tag in the
 *                   NodeInfo_t structure at <b>.
 *
 *      Returns:
 *          -1 if tag <a> <  tag <b.tag>
 *           0 if tag <a> == tag <b.tag>
 *          +1 if tag <a> >  tag <b.tag>
 *
 *-------------------------------------------------------------------------*/
static int FindNode(const void *a, const void *b)
{
        Tag_t      *tag = (Tag_t *)a;
        NodeInfo_t *node = (NodeInfo_t *)b;

        if (tag->domainID < node->tag.domainID) return(-1);
        if (tag->domainID > node->tag.domainID) return(1);
        if (tag->index < node->tag.index) return(-1);

        return(tag->index > node->tag.index);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    FindNodeIndex
 *      Description: Search the <nodeInfo> array for a node with the
 *                   specified tag.
 *
 *      Arguments:
 *          tag       tag of the node to search for
 *          nodeInfo  array of NodeInfo_t structures
 *          numNodes  number of elements in array <nodeInfo>
 *
 *      Returns:  Index in <nodeInfo> of the element whose tag matches <tag>
 *                or -1 if not found.
 *
 *-------------------------------------------------------------------------*/
static int FindNodeIndex(Tag_t *tag, NodeInfo_t *nodeInfo, int numNodes)
{
        NodeInfo_t *found;

        found = (NodeInfo_t *)bsearch(tag, nodeInfo, numNodes,
                                      sizeof(NodeInfo_t), FindNode);

        if (found == (NodeInfo_t *)NULL) {
            return(-1);
        } else {
            return(found - nodeInfo);
        }
}


/*---------------------------------------------------------------------------
 *
 *      Function:    FindDensity
 *      Description: Given a list of segments, find the location of both
 *                   endpoints (if known) and calculate the contribution
 *                   to all appropriate density grid cells for each segment.
 *
 *      Arguments:
 *          segList       array of segments
 *          numAllocSegs  number of allocated elements in array <segList>
 *          numSegs       number of used elements in array <segList>
 *          nodeInfo      array of nodal data
 *          numNodes      number of nodes in array <nodeInfo>
 *          densityGrid   3-D density grid represented as a single
 *                        dimension array.
 *          gridSize      specifies the size of the density grid in
 *                        X, Y and Z dimensions.
 *          cellSize      length (units of b) of density grid cells
 *                        in X, Y and Z.
 *          burgVolFactor  factor by which segment length is adjusted to
 *                       calculate density.
 *
 *-------------------------------------------------------------------------*/
static void FindDensity(Home_t *home, Seg_t *segList, int *numAllocSegs,
                        int *numSegs, NodeInfo_t *nodeInfo, int numNodes,
                        real8 *densityGrid, int gridSize[3],
                        real8 cellSize[3], real8 burgVolFactor)
{
        int     i, j, index1, index2;
        int     cellIndex[3];
        real8   p1[3], p2[3], lineDir[3];
        Param_t *param;

        param = home->param;

/*
 *      Loop through all the segments currently on the list.  If
 *      we have the node info (position) of eahc endpoint loaded
 *      into memory, we can process the segment.  Otherwise skip
 *      it for now.
 */
        for (i = 0; i < *numSegs; i++) {

            index1 = FindNodeIndex(&segList[i].tag1, nodeInfo, numNodes);
            index2 = FindNodeIndex(&segList[i].tag2, nodeInfo, numNodes);

            if ((index1 < 0) || (index2 < 0)) {
                continue;
            }

            for (j = 0; j < 3; j++) {
               p1[j] = nodeInfo[index1].pos[j];
               p2[j] = nodeInfo[index2].pos[j];
            }

            PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p2[X], &p2[Y], &p2[Z]);
/*
 *          Get the direction of the segment from point p1 to p2.
 */
            for (j = 0; j < 3; j++) {
                lineDir[j] = (p2[j] >= p1[j] ? 1 : -1);
            }
/*
 *          Get the cell indices in the density grid for 
 *          point p1.  If the point is actually outside the
 *          primary simulation space, the cell indices may
 *          be less than zero or greater than the maximum
 *          index, but that is okay since the indices
 *          will be shifted for PBC (or to the closest primary
 *          cell if PBC is not enabled).
 */
            for (j = 0; j < 3; j++) {

                cellIndex[j] = (int) floor((p1[j] - param->minCoordinates[j]) /
                               cellSize[j]);
            }

/*
 *          Increment the density values for all cells in the
 *          density field containing any portion of this
 *          segment.
 */
            SetDensity(home, densityGrid, gridSize, cellSize,
                       p1, p2, lineDir, cellIndex, burgVolFactor);

/*
 *          Decrement the reference count for each of the segment
 *          endpoints, and remove this segment from the list by
 *          copying the last segment in the list over top of the
 *          current one.
 */
            nodeInfo[index1].refCnt--;
            nodeInfo[index2].refCnt--;

            if (i < (*numSegs - 1)) {
                memcpy(&segList[i], &segList[*numSegs-1], sizeof(Seg_t));
            }

            *numSegs -= 1;
/*
 *          And decrement i so we process the segment we just moved
 *          in the array.
 */
            i--;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    DeleteUnreferencedNodes
 *      Description: Go through all the nodes in <nodeInfo> and removed
 *                   any nodes whose reference count is zero.
 *
 *      Arguments:
 *          nodeInfo      array of nodal data
 *          numNodes      pointer to number of nodes in array <nodeInfo>
 *                        Will be updated to contain the remaining number
 *                        of nodes in <nodeInfo>.
 *
 *-------------------------------------------------------------------------*/
static void DeleteUnreferencedNodes(NodeInfo_t *nodeInfo, int *numNodes)
{
        int i;

        for (i = 0; i < *numNodes; i++) {

            if (nodeInfo[i].refCnt == 0) {

                if (i < (*numNodes - 1)) {
                    memcpy(&nodeInfo[i], &nodeInfo[*numNodes-1],
                           sizeof(NodeInfo_t));
                }

                *numNodes -= 1;
                i--;
            }
        }

        return;
}


static void Usage(char *progName)
{
        printf("%-15s:  [-b] \n", progName);
        printf("                  [-d <dataFile>]\n");
        printf("                  [-g <gridX[,gridY,gridZ]>]\n");
        printf("                  [-o <outFile>]\n");
        printf("                  <ctrlFile>\n");
        printf("\n");

        printf("  Where:\n");
        printf("      -b    indicates the use of an HDF5 restart file\n");
        printf("      -d    specifies the name of the nodal data file.  If\n");
        printf("            none is given default behavior is the same as\n");
        printf("            the ParaDiS default for data file name.\n");
        printf("      -g    specifies the density grid size.  If the the\n");
        printf("            size in Y and Z dimensions are not given they\n");
        printf("            default to the size in the X dimension.  Default\n");
        printf("            is a 25X25X25 density grid.\n");
        printf("      -o    specifies the name of the density output file.  If\n");
        printf("            none is given default behavior is the same as\n");
        printf("            the base control file name with a '.dens' suffix\n");
        printf("            appended\n");
        printf("\n");

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    ReallocNodeInfo
 *      Description: Increase the size of the <nodeInfo> array if necessary.
 *
 *      Arguments:
 *          nodeInfo      pointer to array of nodal data.  Will be updated
 *                        for caller if the node array is reallocated.
 *          numAllocNodes pointer to size (number of allocated elements) of
 *                        array <nodeInfo>.  Will be updated for caller
 *                        if the node array is reallocated.
 *          numNodes      number of nodes in array <nodeInfo>.
 *
 *-------------------------------------------------------------------------*/
static void ReallocNodeInfo(NodeInfo_t **nodeInfo, int *numAllocNodes, int numNodes)
{
        int   neededNodes;

        neededNodes = numNodes + NODE_BLOCK_INC;

        if (neededNodes > *numAllocNodes) {
            *nodeInfo = (NodeInfo_t *)realloc(*nodeInfo, neededNodes *
                                              sizeof(NodeInfo_t));
            if (*nodeInfo == (NodeInfo_t *)NULL) {
                Fatal("Error reallocating NodeInfo buffer");
            }
            *numAllocNodes = neededNodes;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    ReallocSegList
 *      Description: Increase the size of the <segList> array if necessary.
 *
 *      Arguments:
 *          segList       pointer to array of segment data.  Will be updated
 *                        for caller if the segment array is reallocated.
 *          numAllocSegs  pointer to size (number of allocated elements) of
 *                        array <segList>.  Will be updated for caller
 *                        if the segment array is reallocated.
 *          numSegs       number of segments in array <segList>.
 *
 *-------------------------------------------------------------------------*/
static void ReallocSegList(Seg_t **segList, int *numAllocSegs, int numSegs)
{
        int   neededSegs;

        neededSegs = numSegs + SEG_BLOCK_INC;

        if (neededSegs > *numAllocSegs) {
            *segList = (Seg_t *)realloc(*segList, neededSegs * sizeof(Seg_t));
            if (*segList == (Seg_t *)NULL) {
                Fatal("Error reallocating segList buffer");
            }
            *numAllocSegs = neededSegs;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    ReadDataParams
 *      Description: Read the data file parameters from a text version
 *                   nodal data file.
 *
 *      Arguments:
 *          fpData  file pointer to open nodal data file.
 *
 *-------------------------------------------------------------------------*/
static void ReadDataParams(Home_t *home, FILE *fpData)
{
        int         maxTokenLen, tokenType, pIndex, valType, numVals, okay;
        int         doBinRead;
        char        token[256];
        void        *valList;
        Param_t     *param;
        ParamList_t *dataParamList;

        doBinRead = 0;
        maxTokenLen = sizeof(token);

        param = home->param;
        dataParamList = home->dataParamList;

/*
 *      Get the first token.  This should either be a known
 *      parameter identifier, or a file version number.
 */
        tokenType = GetNextToken(fpData, token, maxTokenLen);
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

            ReadPreV4DataParams(home, fpData, &home->decomp);

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
                    okay = GetParamVals(fpData, valType, numVals, valList);
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

                        tokenType = GetNextToken(fpData, token, maxTokenLen);
/*
 *                      Do a quick verification of the decomposition type
 */
                        if ((param->dataDecompType < 1) ||
                            (param->dataDecompType > 2)) {
                            Fatal("dataDecompType=%d is invalid.  Type must be 1 or 2\n",
                                  param->dataDecompType);
                        }

                        ReadDecompBounds(home, (void **)&fpData, doBinRead,
                                         param->decompType,
                                         &home->decomp);
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
                        tokenType = GetNextToken(fpData, token, maxTokenLen);
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
                        okay = GetParamVals(fpData, valType, numVals,
                                            valList);
                    }
                }

                tokenType = GetNextToken(fpData, token, maxTokenLen);

                if ((tokenType == TOKEN_NULL)||(tokenType == TOKEN_ERR)) {
                    Fatal("Parsing error in data file\n");
                }
                pIndex = LookupParam(dataParamList, token);
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    ProcessNodalData
 *      Description: Loop through all segments of the text nodal data file.
 *                   Read node positions, create the segment list, and
 *                   calculate the density cell contribution from each segment
 *
 *      Arguments:
 *          fpData       file pointer to open nodal data file.
 *          baseFileName Name of the nodal data file stripped of any
 *                       sequence number suffix.
 *          densityGrid  3-D density grid represented as a single
 *                       dimension array.
 *          gridSize     specifies the size of the density grid in
 *                       X, Y and Z dimensions.
 *          cellSize     length (units of b) of density grid cells
 *                       in X, Y and Z.
 *          burgVolFactor  factor by which segment length is adjusted to
 *                       calculate density.
 *
 *-------------------------------------------------------------------------*/
static void ProcessNodalData(Home_t *home, FILE *fpData, char *baseFileName,
                             real8 *densityGrid, int gridSize[3],
                             real8 cellSize[3], real8 burgVolFactor)
{
        int        i, j;
        int        numArms, fileSeqNum;
        int        numSegs, numAllocSegs;
        int        numNodes, numAllocNodes;
        int        nodeIndex, segIndex, readThisBlock;
        char       inLine[256];
        char       segmentName[256];
        Tag_t      nbrTag;
        Seg_t      *segList;
        Param_t    *param;
        NodeInfo_t *nodeInfo;

        param = home->param;

        fileSeqNum    = 0;
        nodeIndex     = 0;
        segIndex      = 0;
        readThisBlock = 0;

        numSegs      = 0;
        numAllocSegs = 0;
        segList      = (Seg_t *)NULL;

        numNodes      = 0;
        numAllocNodes = 0;
        nodeInfo      = (NodeInfo_t *)NULL;

/*
 *      Set up some initial nodeInfo and segment arrays.
 */
        ReallocSegList(&segList, &numAllocSegs, numSegs);
        ReallocNodeInfo(&nodeInfo, &numAllocNodes, numNodes);

/*
 *      Loop until we've read every node in the system.   Don't
 *      assume we can read the whole thing into memory in one shot,
 *      though.  Read in blocks, and at the end of each block
 *      process any segments for which we have both endpoint
 *      coordinates.
 */
        for (i = 0; i < param->nodeCount; i++) {

            while (1) {
                Getline(inLine, sizeof(inLine), fpData);
/*
 *              If we've hit EOF, open up the next data file segment
 *              (if any).
 */
                if (inLine[0] == 0) {
                    fileSeqNum++;
                    if (fileSeqNum >= param->numFileSegments) {
                        Fatal("Unexpected EOF in nodal data file");
                    }
                    fclose(fpData);
                    sprintf(segmentName, "%s.%d", baseFileName, fileSeqNum);
                    fpData = fopen(segmentName, "r");
                    if (fpData == (FILE *)NULL) {
                        Fatal("Error %d opening node data file %s",
                              errno, segmentName);
                    }
                } else {
                    break;
                }
            }

/*
 *          Got info on the next node.  If we've read already read in
 *          the maximum number of nodes for this block, process any
 *          data we can from the buffers before doing anything else.
 */
            if (readThisBlock > NODE_READ_BLOCK_CNT) {

                FindDensity(home, segList, &numAllocSegs, &numSegs,
                            nodeInfo, numNodes, densityGrid,
                            gridSize, cellSize, burgVolFactor);

                DeleteUnreferencedNodes(nodeInfo, &numNodes);

                readThisBlock = 0;
                nodeIndex = numNodes;
                segIndex = numSegs;
            }

/*
 *          If the nodeInfo array is full, increase it's size.
 */
            if (numNodes == numAllocNodes) {
                ReallocNodeInfo(&nodeInfo, &numAllocNodes, numNodes);
            }

            sscanf(inLine, "%d,%d %lf %lf %lf %d %*d",
                   &nodeInfo[nodeIndex].tag.domainID,
                   &nodeInfo[nodeIndex].tag.index,
                   &nodeInfo[nodeIndex].pos[X], &nodeInfo[nodeIndex].pos[Y],
                   &nodeInfo[nodeIndex].pos[Z], &numArms);

            FoldBox(param, &nodeInfo[nodeIndex].pos[X], &nodeInfo[nodeIndex].pos[Y],
                    &nodeInfo[nodeIndex].pos[Z]);

            nodeInfo[nodeIndex].refCnt = numArms;

/*
 *          Read the node's segment data and add any necessary segments to the
 *          segment list.
 */
            for (j = 0; j < numArms; j++) {

                Getline(inLine, sizeof(inLine), fpData);
                sscanf(inLine, "%d,%d, %*f %*f %*f", &nbrTag.domainID, &nbrTag.index);
                Getline(inLine, sizeof(inLine), fpData);

                if ((nodeInfo[nodeIndex].tag.domainID < nbrTag.domainID) ||
                    ((nodeInfo[nodeIndex].tag.domainID == nbrTag.domainID) &&
                     (nodeInfo[nodeIndex].tag.index < nbrTag.index))) {
/*
 *                  Add the segment to the seg list.  Increase the size of
 *                  the segment list if there's no more room.
 */
                    if (numSegs == numAllocSegs) {
                        ReallocSegList(&segList, &numAllocSegs, numSegs);
                    }

                    segList[segIndex].tag1.domainID =
                            nodeInfo[nodeIndex].tag.domainID;
                    segList[segIndex].tag1.index =
                            nodeInfo[nodeIndex].tag.index;
                    segList[segIndex].tag2.domainID = nbrTag.domainID;
                    segList[segIndex].tag2.index = nbrTag.index;

                    segIndex++;
                    numSegs++;
                }
            }

            nodeIndex++;
            numNodes++;
            readThisBlock++;
        }

/*
 *      Process any remaining segments in the list
 */
        FindDensity(home, segList, &numAllocSegs, &numSegs,
                    nodeInfo, numNodes, densityGrid,
                    gridSize, cellSize, burgVolFactor);

        return;
}


#ifdef USE_HDF
/*---------------------------------------------------------------------------
 *
 *      Function:    ProcessBinNodalData
 *      Description: Loop through all segments of the HDF nodal data file.
 *                   Read node positions, create the segment list, and
 *                   calculate the density cell contribution from each segment
 *
 *      Arguments:
 *          fileID       file pointer to open HDF5 file.
 *          baseFileName Name of the nodal data file stripped of any
 *                       sequence number suffix.
 *          densityGrid  3-D density grid represented as a single
 *                       dimension array.
 *          gridSize     specifies the size of the density grid in
 *                       X, Y and Z dimensions.
 *          cellSize     length (units of b) of density grid cells
 *                       in X, Y and Z.
 *          burgVolFactor  factor by which segment length is adjusted to
 *                       calculate density.
 *
 *-------------------------------------------------------------------------*/
static void ProcessBinNodalData(Home_t *home, hid_t fileID, char *baseFileName,
                                real8 *densityGrid, int gridSize[3],
                                real8 cellSize[3], real8 burgVolFactor)
{
        int        i;
        int        totNodesRead;
        int        nextTaskID, nextFileSeg, maxFileSeg;
        int        taskNodeCount, taskSegCount;
        int        numSegs, numAllocSegs;
        int        numNodes, numAllocNodes;
        int        nodeIndex, segIndex, readThisBlock;
        int        *nodeIndices, *nodeNumSegs, *segTags;
        int        taskIDRange[2];
        real8      *nodePos;
        char       tmpFileName[256], itemName[64], taskDir[64];
        Seg_t      *segList;
        Param_t    *param;
        NodeInfo_t *nodeInfo;


        param = home->param;

        maxFileSeg = param->numFileSegments - 1;

        nextFileSeg = 1;
        nextTaskID = 0;
        totNodesRead = 0;

        nodeIndex     = 0;
        segIndex      = 0;
        readThisBlock = 0;

        numSegs      = 0;
        numAllocSegs = 0;
        segList      = (Seg_t *)NULL;

        numNodes      = 0;
        numAllocNodes = 0;
        nodeInfo      = (NodeInfo_t *)NULL;

/*
 *      Find the range of tasks for which data is included in this
 *      first data file segment.
 */
        sprintf(itemName, "/taskIDRange");
        int status = ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT, 2, taskIDRange);

        if (status != 0) {
            Fatal("ProcessBinNodalData: Error reading taskIDRange");
        }

/*
 *      Loop until we've read every node in the system.   Don't
 *      assume we can read the whole thing into memory in one shot,
 *      though.  Read in blocks, and at the end of each block
 *      process any segments for which we have both endpoint
 *      coordinates.
 */
        while (totNodesRead < param->nodeCount) {

/*
 *          If we don't have an open data file segment, open the next one
 *          now.
 */
            if ((fileID < 0) && (nextFileSeg <= maxFileSeg)) {

                sprintf(tmpFileName, "%s.%d", baseFileName, nextFileSeg);
                fileID = H5Fopen(tmpFileName, H5F_ACC_RDONLY, H5P_DEFAULT);

                if (fileID < 0) {
                    Fatal("Task %d: Error opening %s", home->myDomain,
                          tmpFileName);
                }
/*
 *              Find the min/max ID's of the tasks that wrote nodal
 *              data to this data file.
 */
                sprintf(itemName, "/taskIDRange");
                status = ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                        2, taskIDRange);

                if (status != 0) {
                    Fatal("ProcessBinNodalData: Error reading %s taskIDRange",
                          tmpFileName);
                }

                nextTaskID = taskIDRange[0];

            }

/*
 *          If we don't have an open file descriptor at this point, it's
 *          because we've read all the data, so break out of the loop.
 */
            if (fileID < 0) {
                break;
            }

/*
 *          Loop through all the task(i.e. domain) directories in this
 *          file.  Read the nodal data and process it as necesary.
 */
            while (nextTaskID <= taskIDRange[1]) {

/*
 *              Get the node count for this task.  If there is no
 *              nodal data in the task's directory, we don't need to
 *              read anything else here.
 */
                sprintf(taskDir, "/task%d", nextTaskID);
                sprintf(itemName, "%s/nodeCount", taskDir);

                status = ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                        1, &taskNodeCount);

                if (status != 0) {
                    Fatal("ReadBinDataFile: Error reading %s %s",
                          tmpFileName, itemName);
                }

                if (taskNodeCount <= 0) {
                    nextTaskID++;
                    continue;
                }

/*
 *              Find out how much segment data is stored
 */
                sprintf(itemName, "%s/segCount", taskDir);
                status = ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                            1, &taskSegCount);
                if (status != 0) {
                    Fatal("ReadBinDataFile: Error reading %s %s",
                          tmpFileName, itemName);
                }

/*
 *              If we've already read in a ssufficient number of nodes
 *              in this loop, process the data that's already in memory
 *              before reading more data.
 */
                if ((readThisBlock != 0) &&
                    ((readThisBlock + taskNodeCount) > NODE_READ_BLOCK_CNT)) {

                    FindDensity(home, segList, &numAllocSegs, &numSegs,
                                nodeInfo, numNodes, densityGrid,
                                gridSize, cellSize, burgVolFactor);

                    DeleteUnreferencedNodes(nodeInfo, &numNodes);

                    readThisBlock = 0;
                }

/*
 *              Allocate arrays for the binary node and segment data
 *              and read in the arrays from the file.
 */
                nodeIndices = (int *)malloc(taskNodeCount * sizeof(int));
                nodeNumSegs = (int *)malloc(taskNodeCount * sizeof(int));
                nodePos = (real8 *)malloc(taskNodeCount * 3 * sizeof(real8));
                segTags = (int *)malloc(taskSegCount * 3 * sizeof(int));


                sprintf(itemName, "%s/nodeIndex", taskDir);
                status |= ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT, taskNodeCount, nodeIndices);

                sprintf(itemName, "%s/nodeNumSegs", taskDir);
                status |= ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT, taskNodeCount, nodeNumSegs);

                sprintf(itemName, "%s/nodePos", taskDir);
                status |= ReadHDFDataset(fileID, itemName, H5T_NATIVE_DOUBLE, taskNodeCount * 3, nodePos);

                sprintf(itemName, "%s/segTags", taskDir);
                status |= ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT, taskSegCount * 3, segTags);

                if (status != 0) {
                    Fatal("ReadBinDataFile: Error reading %s in %s",
                          taskDir, tmpFileName);
                }

/*
 *              Loop through the nodes and store the nodal positions and
 *              reference counts in the nodeInfo array.  If necessary,
 *              increase the size of the nodeInfo array before adding in
 *              the new nodal data.
 */
                for (i = 0; i < taskNodeCount; i++) {

                    if (numNodes == numAllocNodes) {
                        ReallocNodeInfo(&nodeInfo, &numAllocNodes, numNodes);
                    }

                    nodeInfo[numNodes].pos[X] = nodePos[i*3  ];
                    nodeInfo[numNodes].pos[Y] = nodePos[i*3+1];
                    nodeInfo[numNodes].pos[Z] = nodePos[i*3+2];

                    nodeInfo[numNodes].refCnt = nodeNumSegs[i];

                    nodeInfo[numNodes].tag.domainID = nextTaskID;
                    nodeInfo[numNodes].tag.index    = nodeIndices[i];

                    FoldBox(param, &nodeInfo[numNodes].pos[X],
                            &nodeInfo[numNodes].pos[Y],
                            &nodeInfo[numNodes].pos[Z]);

                    numNodes++;
                }

/*
 *              Loop through all the segments just read.  Segment list
 *              in restart file contains all segments that cross domain
 *              boundaries, and all locally owned segments, so be sure
 *              to avoid adding duplicate segments to the segList array.
 */
                for (i = 0; i < taskSegCount; i++) {

                    if ((nextTaskID > segTags[i*3+1]) ||
                        ((nextTaskID == segTags[i*3+1]) &&
                         (segTags[i*3] > segTags[i*3+2]))) {
                        continue;
                    }
/*
 *                  If necessary, increase the size of the segList array
 */
                    if (numSegs == numAllocSegs) {
                        ReallocSegList(&segList, &numAllocSegs, numSegs);
                    }

                    segList[numSegs].tag1.domainID = nextTaskID;
                    segList[numSegs].tag1.index    = segTags[i*3  ];
                    segList[numSegs].tag2.domainID = segTags[i*3+1];
                    segList[numSegs].tag2.index    = segTags[i*3+2];

                    numSegs++;
                }

/*
 *              Free up the task/domain specific node and segment arrays
 */
                free(nodeIndices);
                free(nodeNumSegs);
                free(nodePos);
                free(segTags);

                totNodesRead  += taskNodeCount;
                readThisBlock += taskNodeCount;

                nextTaskID++;
            }

            H5Fclose(fileID);
            fileID = -1;
            nextFileSeg++;
        }


/*
 *      Process any remaining segments in the list
 */
        FindDensity(home, segList, &numAllocSegs, &numSegs,
                    nodeInfo, numNodes, densityGrid,
                    gridSize, cellSize, burgVolFactor);

        return;
}
#endif


int main(int argc, char *argv[])
{
        int         i, j, k, fd, offset;
        int         doBinRead, fileSeqNum;
        int         numGridElements;
        int         gridSize[3];
        real8       cellVol, burgVolFactor;
        real8       cellSize[3], centerCoord[3];
        real8       *densityGrid;
        char        *dataFile, *ctrlFile, *outFile;
        char        *sep, *start, *argValue, *argTok;
        char        tmpDataFile[256], tmpOutFile[256], testFile[256];
        char        tmpFileName[256], baseFileName[256];
        Home_t      home;
        Param_t     *param;
        FILE        *fpData, *fpOut;
#ifdef USE_HDF
        hid_t       fileID;

        fileID = -1;
#endif

        ctrlFile = (char *)NULL;
        dataFile = (char *)NULL;
        outFile  = (char *)NULL;

        doBinRead = 0;
        fileSeqNum = 0;

        gridSize[X] = 25;
        gridSize[Y] = 25;
        gridSize[Z] = 25;

        memset(&home, 0, sizeof(Home_t));
/*
 *      Get command line options
 */
        for (i = 1; i < argc; i++) {

            if (!strcmp(argv[i], "-b")) {
                doBinRead = 1;
            } else if (!strcmp(argv[i], "-d")) {
                if (i >= (argc - 1)) {
                    Usage(argv[0]);
                    exit(1);
                }
                dataFile = argv[++i];
            } else if (!strcmp(argv[i], "-h")) {
                Usage(argv[0]);
                exit(0);
            } else if (!strcmp(argv[i], "-g")) {
                if (i >= (argc - 1)) {
                    Usage(argv[0]);
                    exit(1);
                }
                argValue = argv[++i];
                argTok = strtok(argValue, ",");
                gridSize[X] = atoi(argTok);
                argTok = strtok(NULL, ",");
                if (argTok == (char *)NULL) {
                    gridSize[Y] = gridSize[X];
                    gridSize[Z] = gridSize[X];
                    continue;
                }
                gridSize[Y] = atoi(argTok);
                argTok = strtok(NULL, ",");
                if (argTok == (char *)NULL) {
                    Usage(argv[0]);
                    exit(1);
                }
                gridSize[Z] = atoi(argTok);
            } else if (!strcmp(argv[i], "-o")) {
                if (i >= (argc - 1)) {
                    Usage(argv[0]);
                    exit(1);
                }
                outFile = argv[++i];
            } else {
                ctrlFile = argv[i];
            }
        }

        if ((gridSize[X] < 1) ||
            (gridSize[Y] < 1) ||
            (gridSize[Z] < 1)) {
            Fatal("Invalid density grid size %dX%dX%d",
                  gridSize[X], gridSize[Y], gridSize[Z]);
        }

/*
 *      If no control file was specified, abort.
 */
        if (ctrlFile == (char *)NULL) {
            Usage(argv[0]);
            exit(1);
        }

/*
 *      Allocate the density grid
 */
        numGridElements = gridSize[X] * gridSize[Y] * gridSize[Z];

        densityGrid = (real8 *)calloc(1, numGridElements*sizeof(real8));

        if (densityGrid == (real8 *)NULL) {
            Fatal("calloc error getting density field of size %d elements",
                  numGridElements);
        }

/*
 *      If user did not provide an output file name, just use
 *      the control file name with a '.dens' suffix appended.
 */
        if (outFile == (char *)NULL) {
            sprintf(tmpOutFile, "%s%s", ctrlFile, DENSITY_FILE_SUFFIX);
            outFile = tmpOutFile;
        }

/*
 *      If the user did not specify a data file name, set
 *      a default name based on the base control file name.
 */
        if (dataFile == (char *)NULL) {
            strcpy(tmpDataFile, ctrlFile);
            start = strrchr(tmpDataFile, '/');
            if (start == (char *)NULL) start = tmpDataFile;
            sep = strrchr(start, '.');
            if ((sep != (char *)NULL) && (sep > start)) *sep = 0;

/*
 *          If the user specified that a binary data file is
 *          used set default name accordingly.  If the user
 *          did not specify, try to figure out whether to
 *          use a binary restart or not.
 */
            if (doBinRead) {
                strcat(tmpDataFile, HDF_DATA_FILE_SUFFIX);
            } else {
                strcpy(testFile, tmpDataFile);
                strcat(testFile, HDF_DATA_FILE_SUFFIX);

                if ((fd = open(testFile, O_RDONLY, 0)) < 0) {
                    strcat(testFile, ".0");
                    fd = open(testFile, O_RDONLY, 0);
                }
                if (fd >= 0) {
                    doBinRead = 1;
                    close(fd);
                    strcat(tmpDataFile, HDF_DATA_FILE_SUFFIX);
                } else {
                    strcat(tmpDataFile, NODEDATA_FILE_SUFFIX);
                }

            }

            dataFile = tmpDataFile;

        } else {
/*
 *          User provided a data file name, but didn't indicate
 *          if it was a binary data file or not.  Make a guess
 *          based on the file name.
 */
            if (strstr(dataFile, HDF_DATA_FILE_SUFFIX) != (char *)NULL) {
                doBinRead = 1;
            }

        }

#ifndef USE_HDF
        if (doBinRead) {
            printf("\nERROR:\n");
            printf("The specified restart file is apparently a binary\n");
            printf("restart file, but this utility has not been compiled\n");
            printf("with HDF support.  Please recompile with HDF support\n");
            printf("and try again. \n\n");
            exit(1);
        }
#endif

        home.param = (Param_t *)calloc(1, sizeof(Param_t));
        param = home.param;

        home.ctrlParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));
        home.dataParamList = (ParamList_t *)calloc(1,sizeof(ParamList_t));

/*
 *      Initialize the control and data parameter lists and read the
 *      control file.
 */
        CtrlParamInit(param, home.ctrlParamList);
        DataParamInit(param, home.dataParamList);

        ReadControlFile(&home, ctrlFile);

/*
 *      Open up the nodal data file (or first segment) for reading.
 *      If we can't open the file, exit with an error.
 */
        snprintf(tmpFileName, sizeof(tmpFileName), "%s", dataFile);
 
        if (!strcmp(&tmpFileName[strlen(tmpFileName)-2], ".0")) {
            tmpFileName[strlen(tmpFileName)-2] = 0;
        }

        snprintf(baseFileName, sizeof(baseFileName), "%s", tmpFileName);


        if (doBinRead) {
#ifdef USE_HDF
            struct stat statbuf;
            if (stat(dataFile, &statbuf) == 0) {
                fileID = H5Fopen(baseFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
                if (fileID < 0) {
                    Fatal("Error opening file %s to read nodal data",
                    dataFile);
                }
            } else {
                snprintf(tmpFileName, sizeof(tmpFileName), "%s.%d",
                         baseFileName, fileSeqNum);
                fileID = H5Fopen(tmpFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
                if (fileID < 0) {
                    Fatal("Error %d opening file %s to read nodal data",
                          dataFile);
                }
            }
#endif
        } else {
            if ((fpData = fopen(tmpFileName, "r")) == (FILE *)NULL) {
                snprintf(tmpFileName, sizeof(tmpFileName), "%s.%d",
                         baseFileName, fileSeqNum);
                if ((fpData = fopen(tmpFileName, "r")) == (FILE *)NULL) {
                    Fatal("Error %d opening file %s to read nodal data",
                          errno, dataFile);
                }
            }
        }

/*
 *      Read the data file parameters and do some initializations
 *      before reading the rest of the nodal data
 */
        if (doBinRead) {
#ifdef USE_HDF
            int status = ReadBinDataParams(&home, fileID);
            if (status != 0) {
                Fatal("ReadBinDataFile: Error reading basic params from %s",
                      dataFile);
            }
#endif
        } else {
            ReadDataParams(&home, fpData);
        }


        param->Lx = param->maxCoordinates[X] - param->minCoordinates[X];
        param->Ly = param->maxCoordinates[Y] - param->minCoordinates[Y];
        param->Lz = param->maxCoordinates[Z] - param->minCoordinates[Z];

        if ((param->xBoundType == Free) ||
            (param->yBoundType == Free) ||
            (param->zBoundType == Free)) {

            if (param->xBoundType == Periodic) {
                param->xBoundMin = param->minCoordinates[X];
                param->xBoundMax = param->maxCoordinates[X];
            }
            if (param->yBoundType == Periodic) {
                param->yBoundMin = param->minCoordinates[Y];
                param->yBoundMax = param->maxCoordinates[Y];
            }
            if (param->zBoundType == Periodic) {
                param->zBoundMin = param->minCoordinates[Z];
                param->zBoundMax = param->maxCoordinates[Z];
            }
        }

/*
 *      Calculate the volume of the full simulation, the size of
 *      each density grid cell, etc.
 */
#ifdef _FEM
        if (param->mesh_type == 1) {
            cellVol = (param->xBoundMax - param->xBoundMin) *
                      (param->yBoundMax - param->yBoundMin) *
                      (param->zBoundMax - param->zBoundMin);
        }

        if (param->mesh_type == 2) {
            cellVol = M_PI * param->fem_radius * param->fem_radius *
                      param->fem_height;
        }
#else
        if ((param->zBoundType == Free) ||
            (param->yBoundType == Free) ||
            (param->xBoundType == Free)) {

            cellVol = (param->xBoundMax-param->xBoundMin) *
                      (param->yBoundMax-param->yBoundMin) *
                      (param->zBoundMax-param->zBoundMin);
        } else {
            cellVol = param->Lx * param->Ly * param->Lz;
        }
#endif

        burgVolFactor = 1.0 / (param->burgMag * param->burgMag * cellVol);

        cellSize[X] = param->Lx / (real8)gridSize[X];
        cellSize[Y] = param->Ly / (real8)gridSize[Y];
        cellSize[Z] = param->Lz / (real8)gridSize[Z];

/*
 *      Now need to read the nodal data and process it.  There may
 *      be too much data to read all into memory, so read it in
 *      in blocks and process as much as we can before reaing the
 *      next block of data.
 *
 *      Note: file descriptor for open data file will be passed to
 *      the functions to proces the nodal data.  The function will
 *      close the data file when it is done.
 */
        if (doBinRead) {
#ifdef USE_HDF
            ProcessBinNodalData(&home, fileID, baseFileName, densityGrid,
                                gridSize, cellSize, burgVolFactor);
#endif
        } else {
            ProcessNodalData(&home, fpData, baseFileName, densityGrid,
                             gridSize, cellSize, burgVolFactor);
        }
        
/*
 *      Density field has been updated with contributions from
 *      all segments, so now just write out the density file.
 */
        if ((fpOut = fopen(outFile, "w")) == (FILE *)NULL) {
            Fatal("CalcDensity: Open error %d on %s",
                  errno, outFile);
        }

        fprintf(fpOut, "#\n");
        fprintf(fpOut, "#   This file contains the dislocation density\n");
        fprintf(fpOut, "#   specified at the center of each cell in the\n");
        fprintf(fpOut, "#   density grid. Each data line consists of the\n");
        fprintf(fpOut, "#   center coordinates (X, Y and Z) of a single\n");
        fprintf(fpOut, "#   density grid cell followed by the density\n");
        fprintf(fpOut, "#   associated with that cell.\n");
        fprintf(fpOut, "#\n");
        fprintf(fpOut, "#  Grid size (in X,Y,Z):  %d %d %d\n",
                gridSize[X], gridSize[Y], gridSize[Z]);
        fprintf(fpOut, "#\n");

        for (i = 0; i < gridSize[X]; i++) {
            centerCoord[X] = param->minCoordinates[X] +
                             ((0.5 + (real8)i) * cellSize[X]);
            for (j = 0; j < gridSize[Y]; j++) {
                centerCoord[Y] = param->minCoordinates[Y] +
                                 ((0.5 + (real8)j) * cellSize[Y]);
                for (k = 0; k < gridSize[Z]; k++) {
                    centerCoord[Z] = param->minCoordinates[Z] +
                                     ((0.5 + (real8)k) * cellSize[Z]);
                    offset = i * gridSize[Y] * gridSize[Z] +
                             j * gridSize[Z] + k;
                    fprintf(fpOut, "%13.5e %13.5e %13.5e %13.5e  "
                            "# cell = %d,%d,%d\n",
                            centerCoord[X], centerCoord[Y], centerCoord[Z],
                            densityGrid[offset], i, j, k);
                }
            }
        }

        fclose(fpOut);

        exit(0);
        return(0);
}
