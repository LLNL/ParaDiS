/****************************************************************************
 *
 *      Module:      ParadisConvert.c
 *      Description: Converts older style control parameter and
 *                   nodal data files to the current format
 *
 ***************************************************************************/
#include <stdarg.h>
#include "Home.h"
#include "Parse.h"
#include "Decomp.h"
#include "RSDecomp.h"
#include "Restart.h"


static void Usage(char *progName)
{
        printf("\nUsage:  %s  <controlFile> [<dataFile>]\n\n", progName);
        return;
}

/*-----------------------------------------------------------------------
 *
 *      Function:     GetConfigFromCtrlFile
 *      Description:  Reads the data associated with the "config="
 *                    specification of the original style control
 *                    files.  Data format for the specification is
 *                    config = [
 *                      <minimum boundary coordinates>
 *                      <maximum boundary coordinates>
 *                      <number of burgers vectors:==nb>
 *                      [1  bx by bz
 *                       2  bx by bz
 *                       ...
 *                       nb bx by bz]
 *                      <node count>
 *                      <nodal data>
 *                      <numXDoms numYDoms numZDoms>
 *                      [domain decomposition>
 *                    ]
 *
 *      Arguments:
 *          fp         File stream from which to read the old-style
 *                     configuration data 
 *          nodeConfig Location in which to return to the caller an
 *                     array of Node_t structures containing all nodal
 *                     data from the file.
 *
 *----------------------------------------------------------------------*/
static void GetConfigFromCtrlFile(Home_t *home, FILE *fp, Node_t **nodeConfig)
{
        int    i, j, tokenType;
        int    nodeCount, numNbrs, burgID, nodeWeight;
        int    numXDoms, numYDoms, numZDoms, saveDecomp;
        int    numBurgersVectors, *burgersVectors;
        int    binRead = 0;
        char   token[MAX_STRING_LEN];
        char   inLine[256], oldTagStr[32];
        char   *separator;
        real8  boundary[6];
        real8  *burgX=0, *burgY=0, *burgZ=0;
        Node_t *nodeList=0, *node=0;
        Param_t *param=0;
        ParamList_t *dataParamList;

        param = home->param;
        dataParamList = home->dataParamList;

/*
 *      The remainder of the current line should contain "= ["
 *      so just read that off before we pull in the real data.
 */
        Getline(inLine, sizeof(inLine), fp);

/*
 *      Get the boundary coordinates of the problem space and
 *      set some parameters dependent on the simulation boundary
 *      coordinates
 */
        Getline(inLine, sizeof(inLine), fp);
        sscanf(inLine, "%lf %lf %lf", &param->minCoordinates[0],
               &param->minCoordinates[1], &param->minCoordinates[2]);

        Getline(inLine, sizeof(inLine), fp);
        sscanf(inLine, "%lf %lf %lf", &param->maxCoordinates[0],
               &param->maxCoordinates[1], &param->maxCoordinates[2]);

        param->minSideX = param->minCoordinates[X];
        param->minSideY = param->minCoordinates[Y];
        param->minSideZ = param->minCoordinates[Z];

        param->maxSideX = param->maxCoordinates[X];
        param->maxSideY = param->maxCoordinates[Y];
        param->maxSideZ = param->maxCoordinates[Z];

        param->Lx = param->maxSideX - param->minSideX;
        param->Ly = param->maxSideY - param->minSideY;
        param->Lz = param->maxSideZ - param->minSideZ;

/*
 *      There may be an array of burgers vectors provided.  If so
 *      read in the array.
 */
        Getline(inLine, sizeof(inLine), fp);
        sscanf(inLine, "%d", &numBurgersVectors);

        if (numBurgersVectors > 0) {

            burgX = (real8 *)calloc(1, numBurgersVectors * sizeof(real8));
            burgY = (real8 *)calloc(1, numBurgersVectors * sizeof(real8));
            burgZ = (real8 *)calloc(1, numBurgersVectors * sizeof(real8));

            for (i = 0; i < numBurgersVectors; i++) {
                    Getline(inLine, sizeof(inLine), fp);
                    sscanf(inLine, "%*d %lf %lf %lf\n", 
                           &burgX[i], &burgY[i], &burgZ[i]);
            }
        }

/*
 *      Next is the node count followed by nodal data
 */
        Getline(inLine, sizeof(inLine), fp);
        sscanf(inLine, "%d", &nodeCount);

        if (nodeCount > 0) {
            nodeList = (Node_t *)calloc(1, nodeCount * sizeof(Node_t));
        }

        for (i = 0; i < nodeCount; i++) {

            node = &nodeList[i];

            Getline(inLine, sizeof(inLine), fp);
            
            sscanf(inLine, "%*1c %[^ ] %d %lf %lf %lf %d %d %d %d",
                       oldTagStr, &nodeWeight, &node->x, &node->y,
                       &node->z, &numNbrs, &node->constraint,
                       &node->myTag.domainID, &node->myTag.index);

/*
 *          With the version 5 nodal data file, the values for
 *          constraints changed some so we have to remap the constraint
 *          numbers.
 */
            MapV4ToV5Constraint(&node->constraint);
                   
/*
 *          Older format config files use an integer global id;
 *          new format files use a tag of the form
 *          <domain>,<index>.  For now, need to handle both formats.
 *          Use this to overwrite the node->tag since we no
 *          longer have a need for both old and new tags in
 *          the Node_t structure.
 */
            if ((separator = strchr(oldTagStr,',')) == (char *)NULL) {
                node->myTag.domainID = 0;
                node->myTag.index = atoi(oldTagStr);
            } else {
                node->myTag.domainID = atoi(oldTagStr);
                node->myTag.index = atoi(separator+1);
            }

/*
 *          Allocate structures for the node's segments and read all
 *          segment data in
 */
            AllocNodeArms(node, numNbrs);

            for (j = 0; j < node->numNbrs; j++) {
                if (numBurgersVectors > 0) {
                    Getline(inLine, sizeof(inLine), fp);
                    sscanf(inLine, "%*[ \t]%[^ ] %d %lf %lf %lf",
                           oldTagStr, &burgID, &node->nx[j],
                           &node->ny[j], &node->nz[j]);
                    node->burgX[j] = burgX[burgID-1];
                    node->burgY[j] = burgY[burgID-1];
                    node->burgZ[j] = burgZ[burgID-1];
                } else {
                    Getline(inLine, sizeof(inLine), fp);
                    sscanf(inLine, "%*[ \t]%[^ ] %lf %lf %lf",
                           oldTagStr, &node->burgX[j],
                           &node->burgY[j], &node->burgZ[j]);
                    Getline(inLine, sizeof(inLine), fp);
                    sscanf(inLine, "%lf %lf %lf", &node->nx[j],
                           &node->ny[j], &node->nz[j]);
                }

                if ((separator=strchr(oldTagStr,','))==(char *)NULL) {
                    node->nbrTag[j].domainID = 0;
                    node->nbrTag[j].index = atoi(oldTagStr);
                } else {
                    node->nbrTag[j].domainID = atoi(oldTagStr);
                    node->nbrTag[j].index = atoi(separator+1);
                }

            }  /* for (j = 0; j < unmNbrs; ...) */

            FoldBox(param, &node->x, &node->y, &node->z);

        }  /* for (i = 0; i < nodeCount; ...) */

/*
 *      Next comes the domain geometry and decomposition.
 */
        Getline(inLine, sizeof(inLine), fp);
        sscanf(inLine, "%d %d %d",
               &param->dataDecompGeometry[X],
               &param->dataDecompGeometry[Y],
               &param->dataDecompGeometry[Z]);

        saveDecomp = 1;
        
/*
 *      The input decomp type should always be one for these files, so the
 *      output decomp type must match it.
 */
        param->dataDecompType = 1;

        if ((param->dataDecompGeometry[X] != param->nXdoms) ||
            (param->dataDecompGeometry[Y] != param->nYdoms) ||
            (param->dataDecompGeometry[Z] != param->nZdoms) ||
            (param->decompType != 1)) {
            saveDecomp = 0;
        }

        ReadDecompBounds(home, (void **)&fp, binRead, param->decompType,
                         &home->decomp);

/*
 *      The next line should just contain a "]" terminating the
 *      "config=[" specification, so just read that in and discard
 *      it before continuing.
 */
        Getline(inLine, sizeof(inLine), fp);

        *nodeConfig = nodeList;
        param->nodeCount = nodeCount;

        return;
}


/*-----------------------------------------------------------------------
 *
 *      Function:     CreateDataFileFromCtrl
 *      Description:  Given an array of node structures and the
 *                    current problem size, domain geometry and
 *                    domain decomposition, create a new nodal
 *                    data file.
 *
 *      Arguments:
 *          dataFileName Name of the data file to be created
 *          nodeList     Array of structures containing the nodal
 *                       data to be written.
 *
 *----------------------------------------------------------------------*/
static int CreateDataFileFromCtrl(Home_t *home, char *dataFileName,
                                  Node_t *nodeList)
{
        int    i, iArm;
        FILE   *fpData;
        Node_t *node;
        Param_t *param;
        ParamList_t *dataParamList;
        char   newFileName[512];

        param = home->param;
        dataParamList = home->dataParamList;

        sprintf(newFileName, "%s.new", dataFileName);

        if ((fpData = fopen(newFileName, "w")) == (FILE *)NULL) {
            Fatal("Open error (%d) on file %s\n", errno, newFileName);
        }

/*
 *      Write the data file parameters
 */
        WriteParam(dataParamList, -1, fpData);

/*
 *      Write the domain decomposition into nodal data file
 *      and then some comment lines describing the nodal
 *      data that will follow.
 */
        fprintf(fpData, "\n#\n#  END OF DATA FILE PARAMETERS\n#\n\n");
        fprintf(fpData, "domainDecomposition = \n");
        WriteDecompBounds(home, fpData);

        fprintf(fpData, "nodalData = \n");
        fprintf(fpData, "#  Primary lines: node_tag, x, y, z, "
                "num_arms, constraint\n");
        fprintf(fpData, "#  Secondary lines: arm_tag, burgx, burgy, "
                "burgz, nx, ny, nz\n");

/*
 *      Dump the nodal data to the data file
 */
        for (i = 0; i < param->nodeCount; i++) {
            node = &nodeList[i];
 
            if ((node->numNbrs == 1) &&
                HAS_NONE_OF_CONSTRAINTS(node->constraint,
                                        PINNED_NODE | SURFACE_NODE))
            {
                Fatal("CreateDataFromCtrl: node (%d,%d) singly linked!",
                      node->myTag.domainID, node->myTag.index);
            }

            fprintf(fpData,
                    " %d,%d %.8f %.8f %.8f %d %d\n",
                    node->myTag.domainID, node->myTag.index,
                    node->x, node->y, node->z, node->numNbrs,
                    node->constraint);

            for (iArm = 0; iArm < node->numNbrs; iArm++) {
                fprintf(fpData, "   %d,%d %16.10f %16.10f %16.10f\n"
                        "       %16.10f %16.10f %16.10f\n",
                        node->nbrTag[iArm].domainID,
                        node->nbrTag[iArm].index, node->burgX[iArm],
                        node->burgY[iArm], node->burgZ[iArm],
                        node->nx[iArm], node->ny[iArm], node->nz[iArm]);
            }
        }

        fclose(fpData);

        return(0);
}


/*-----------------------------------------------------------------------
 *
 *      Function:     ConvertControlFile
 *      Description:  Converts the specified control parameter file
 *                    to the current format.  In the event the specified
 *                    file is the original format which also contains
 *                    the nodal data, a new nodal data file will be
 *                    created and the caller will be informed that
 *                    no additional data file conversion will be
 *                    necessary.
 *
 *      Arguments:
 *          ctrlFileName  Name of control parameter file to be converted
 *          dataFileName  Name of nodal data file.
 *          doDataConvert Location in which to return to the caller a flag
 *                        indicating if the nodal data file still needs to
 *                        be converted (0 == no, 1 == yes)
 *
 *----------------------------------------------------------------------*/
static int ConvertControlFile(Home_t *home, char *ctrlFileName,
                              char *dataFileName, int *doDataConvert)
{

        int  okay, maxTokenLen, tokenType, pIndex;
        int  numVals, valType;
        char newFileName[256], token[256];
        FILE *fpOld, *fpNew;
        void *valList;
        Node_t *nodeList;
        void *decomp;
        Param_t *param;
        ParamList_t *ctrlParamList, *dataParamList;

        param = home->param;

        ctrlParamList = home->ctrlParamList;
        dataParamList = home->dataParamList;

        *doDataConvert = 1;

        sprintf(newFileName, "%s.new", ctrlFileName);

        fpOld = fopen(ctrlFileName, "r");
        fpNew = fopen(newFileName, "w");

        maxTokenLen = sizeof(token);
        tokenType = TOKEN_GENERIC;

        while ((tokenType != TOKEN_ERR) && (tokenType != TOKEN_NULL)) {

/*
 *          Next token (if any) should be a parameter name...
 */
            tokenType = GetNextToken(fpOld, token, maxTokenLen);

            if (tokenType == TOKEN_NULL) {
                printf("Parsing complete!\n");
                break;
            }

            pIndex = LookupParam(ctrlParamList, token);

/*
 *          Obtain values associated with the parameter; don't
 *          save values for unknown/obsolete parameters.
 */
            if (pIndex < 0) {
                if (strcmp(token, "config") == 0) {
/*
 *                  Some very old format files contained the
 *                  nodal data within a 'config' specification
 *                  in the control file.  For those cases, we'll
 *                  read the nodal data and create a brand new
 *                  nodal data file.
 */
                    GetConfigFromCtrlFile(home, fpOld, &nodeList);
                    *doDataConvert = 0;

                    param->dataFileVersion = NODEDATA_FILE_VERSION;
                    param->numFileSegments = 1;
                    param->decompType = 1;

                    if (home->decomp == (void *)NULL) {
                        UniformDecomp(home, &home->decomp);
                    }

                    CreateDataFileFromCtrl(home, dataFileName, nodeList);

                    *doDataConvert = 0;
                    continue;
                } else {
                    printf("Ignoring unrecognized parameter (%s) type (%d)\n",
                           token, tokenType);
                    valType = V_NULL;
                    numVals = 0;
                    valList = (void *)NULL;
                }
            } else {
                valType = ctrlParamList->varList[pIndex].valType;
                numVals = ctrlParamList->varList[pIndex].valCnt;
                valList = ctrlParamList->varList[pIndex].valList;
            }

/*
 *          Obtain the value(s) associated with the parameter name
 *          and write the information back out to the new file in the
 *          current format.
 */
            okay = GetParamVals(fpOld, valType, numVals, valList);

            if (!okay) {
                printf("Error obtaining values for parameter %s\n",
                       ctrlParamList->varList[pIndex].varName);
                exit(1);
                
            }

            if (pIndex >= 0) {
                WriteParam(ctrlParamList, pIndex, fpNew);
            }

        }  /* end parsing loop */

        fclose(fpOld);
        fclose(fpNew);

        return(0);
}


/*-----------------------------------------------------------------------
 *
 *      Function:     GetOldStyleDataParams
 *      Description:  Read the old-style data file in which values
 *                    were dependent on their position in the file.
 *
 *      Arguments:
 *          fp   File stream from which to read the data
 *
 *----------------------------------------------------------------------*/
static void GetOldStyleDataParams(Home_t *home, FILE *fp)
{
        int        saveDecomp;
        int        binRead = 0;
        char       inLine[512];
        RSDecomp_t *newDecomp;
        Param_t    *param;

        param = home->param;

/*
 *      Get the count of file segments
 */
        if (param->dataFileVersion > 1) {
            Getline(inLine, sizeof(inLine), fp);
            sscanf(inLine, "%d", &param->numFileSegments);
        }

/*
 *      Get the boundary coordinates of the problem space
 */
        Getline(inLine, sizeof(inLine), fp);
        sscanf(inLine, "%lf %lf %lf", &param->minCoordinates[0],
               &param->minCoordinates[1], &param->minCoordinates[2]);

        Getline(inLine, sizeof(inLine), fp);
        sscanf(inLine, "%lf %lf %lf", &param->maxCoordinates[0],
               &param->maxCoordinates[1], &param->maxCoordinates[2]);

/*
 *      In version 2 and up, the total node count was relocated
 *      to a position in front of the domain decomposition...
 */
        if (param->dataFileVersion > 1) {
            Getline(inLine, sizeof(inLine), fp);
            sscanf(inLine, "%d", &param->nodeCount);
        }

/*
 *      Get the boundary coordinates of the problem space and
 *      set some parameters dependent on the simulation boundary
 *      coordinates
 */
        Getline(inLine, sizeof(inLine), fp);
        sscanf(inLine, "%d %d %d", &param->dataDecompGeometry[0],
               &param->dataDecompGeometry[1], &param->dataDecompGeometry[2]);

        param->minSideX = param->minCoordinates[X];
        param->minSideY = param->minCoordinates[Y];
        param->minSideZ = param->minCoordinates[Z];

        param->maxSideX = param->maxCoordinates[X];
        param->maxSideY = param->maxCoordinates[Y];
        param->maxSideZ = param->maxCoordinates[Z];

        param->Lx = param->maxSideX - param->minSideX;
        param->Ly = param->maxSideY - param->minSideY;
        param->Lz = param->maxSideZ - param->minSideZ;

/*
 *      Get domain decomposition (if present)
 *
 *      Explicitly set the domain geometry from the control file to match
 *      the domain geometry form the data file.  This forces us to read
 *      the decomposition in the data file and dump it back out as is.
 *      We can do this because the domain geometry from the old control
 *      file has already been written back out to the new control file,
 *      so messing with the value does no harm.
 */
        if (param->dataFileVersion > 2) {
            Getline(inLine, sizeof(inLine), fp);
            sscanf(inLine, "%d", &param->dataDecompType);
        } else {
            param->dataDecompType = param->decompType;
        }

        param->nXdoms = param->dataDecompGeometry[X];
        param->nYdoms = param->dataDecompGeometry[Y];
        param->nZdoms = param->dataDecompGeometry[Z];

        ReadDecompBounds(home, (void **)&fp, binRead, param->dataDecompType,
                         &home->decomp);

/*
 *      If no old domain decomposition was provided, generate
 *      a new single domain RS decomposition.
 */
        if (home->decomp == (void *)NULL) {

            param->decompType = param->dataDecompType = 1;

            param->nXdoms = param->dataDecompGeometry[X] = 1;
            param->nYdoms = param->dataDecompGeometry[Y] = 1;
            param->nZdoms = param->dataDecompGeometry[Z] = 1;

            UniformDecomp(home, &home->decomp);
        }

/*
 *      For version 1 files, the total node count is located
 *      next in the file, so read it in now.
 */
        if (param->dataFileVersion == 1) {
            Getline(inLine, sizeof(inLine), fp);
            sscanf(inLine, "%d", &param->nodeCount);
        }

        return;
}


/*-----------------------------------------------------------------------
 *
 *      Function:     ConvertDataFile
 *      Description:  Reads the contents of the specified nodal data
 *                    file(s) and creates a new set in the current format
 *
 *      Arguments:
 *          baseDataFileName Base file name of the set of nodal data
 *                           file segments
 *
 *----------------------------------------------------------------------*/
static void ConvertDataFile(Home_t *home, char *baseDataFileName)
{
        int  maxTokenLen, tokenType, pIndex;
        int  nextSeqNum, notDone;
        int  oldFileVersion;
        int  binRead = 0;
        char inLine[512], token[256], segmentName[512], newSegName[512];
        FILE *fpOld, *fpNew;
        void *oldDecomp;
        void *valList;
        int  valType, numVals, okay;
        ParamList_t *dataParamList;
        Param_t *param;

        param = home->param;
        dataParamList = home->dataParamList;

        maxTokenLen = sizeof(token);

/*      Figure out the name of the first data file segment */

        sprintf(segmentName, "%s", baseDataFileName);
        sprintf(newSegName, "%s.new", segmentName);
        fpOld = fopen(segmentName, "r");

        if (fpOld == (FILE *)NULL) {
            sprintf(segmentName, "%s.0", baseDataFileName);
            sprintf(newSegName, "%s.0.new", baseDataFileName);
            fpOld = fopen(segmentName, "r");
        }

        if (fpOld == (FILE *)NULL) {
            printf("Open error %d on file %s\n", errno, baseDataFileName);
            exit(1);
        }

        fpNew = fopen(newSegName, "w");
        if (fpNew == (FILE *)NULL) {
            printf("Open error %d on new file %s\n", errno, newSegName);
            exit(1);
        }

/*
 *      Get the first token.  If the token is a known token, use
 *      new parser, otherwise assume it is version number
 */
        tokenType = GetNextToken(fpOld, token, maxTokenLen);
        pIndex = LookupParam(dataParamList, token);

        if (pIndex < 0) {
/*
 *          If the token does not correspond to a know parameter, it
 *          should be the version number, so convert it to an integer
 *          and parse the rest of the file the old way (i.e. values
 *          are positional, not identified by name/value pairs.
 */
            param->dataFileVersion = atoi(token);
            if (param->dataFileVersion < 0) {
                printf("Unsupported version number for file %s\n", newSegName);
                exit(1);
            }

            GetOldStyleDataParams(home, fpOld);

        } else {

            while ((tokenType != TOKEN_ERR) && (tokenType != TOKEN_NULL)) {

                if (pIndex >= 0) {
/*
 *                  Token represents a known parameter identifier, so
 *                  read the associated value(s).
 */
                    valType = dataParamList->varList[pIndex].valType;
                    numVals = dataParamList->varList[pIndex].valCnt;
                    valList = dataParamList->varList[pIndex].valList;
                    okay = GetParamVals(fpOld, valType, numVals, valList);
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

/*
 *                      Explicitly set the domain geometry from the
 *                      control file to match the domain geometry
 *                      form the data file.  This forces us to read
 *                      the decomposition in the data file and dump
 *                      it back out as is.  We can do this because 
 *                      the domain geometry from the old control file
 *                      has already been written back out to the new
 *                      control file, so messing with the value does
 *                      no harm.
 */
                        param->nXdoms = param->dataDecompGeometry[X];
                        param->nYdoms = param->dataDecompGeometry[Y];
                        param->nZdoms = param->dataDecompGeometry[Z];

/*
 *                      For RB decompositions we need to initialize
 *                      some stuff before we read the decomposition data.
 */
                        if (param->decompType == 2) {
                            for (home->xMaxLevel = 0;
                                 param->nXdoms >> home->xMaxLevel > 1;
                                 home->xMaxLevel++);
                
                            for (home->yMaxLevel = 0;
                                 param->nYdoms >> home->yMaxLevel > 1;
                                 home->yMaxLevel++);

                            for (home->zMaxLevel = 0;
                                 param->nZdoms >> home->zMaxLevel > 1;
                                 home->zMaxLevel++);

                            if ((1 << home->xMaxLevel) < param->nXdoms) {
                                home->xMaxLevel++;
                            }

                            if ((1 << home->yMaxLevel) < param->nYdoms) {
                                home->yMaxLevel++;
                            }

                            if ((1 << home->zMaxLevel) < param->nZdoms) {
                                home->zMaxLevel++;
                            }
                        }

                        tokenType = GetNextToken(fpOld, token, maxTokenLen);
                        ReadDecompBounds(home, (void **)&fpOld, binRead,
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
                        tokenType = GetNextToken(fpOld, token, maxTokenLen);
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
                        okay = GetParamVals(fpOld, valType, numVals, valList);
                    }
                }

                tokenType = GetNextToken(fpOld, token, maxTokenLen);

                if ((tokenType == TOKEN_NULL)||(tokenType == TOKEN_ERR)) {
                    Fatal("Parsing error on file %s\n", segmentName);
                }
                pIndex = LookupParam(dataParamList, token);
            }
        }
/*
 *      The domain decomposition type in the control file may not
 *      match the type in the data file.  We've processed the control
 *      file parameters now, so explicitly set param->decompType
 *      to the value from the data file... so when we rewrite the
 *      domain decomposition we use the same decomp type that
 *      was read in from the data file.
 */
        param->decompType = param->dataDecompType;
        
/*
 *      If there is currently no domain decomposition then
 *      generate a uniform decomposition based on the 
 *      geometry in param->dataDecompGeometry
 */
        if (home->decomp == (void **)NULL) {
            UniformDecomp(home, &home->decomp);
        }

/*
 *      Write data file parameters in newest format
 */
        oldFileVersion = param->dataFileVersion;
        param->dataFileVersion = NODEDATA_FILE_VERSION;

        WriteParam(dataParamList, -1, fpNew);

/*
 *      Write domain decomposition
 */
        fprintf(fpNew, "\n#\n#  END OF DATA FILE PARAMETERS\n#\n\n");
        fprintf(fpNew, "domainDecomposition = \n");

        WriteDecompBounds(home, fpNew);

        fprintf(fpNew, "nodalData = \n");
        fprintf(fpNew, "#  Primary lines: node_tag, x, y, z, "
                "num_arms, constraint\n");
        fprintf(fpNew, "#  Secondary lines: arm_tag, burgx, burgy, "
                "burgz, nx, ny, nz\n");

/*
 *      Now just process all the nodal data left in this file and
 *      in the other segments of the data file (if multiple segments
 *      exist)
 */
        nextSeqNum = 1;
        notDone = 1;

        while (notDone) {

/*
 *          If necessary, open the next segment of the nodal data file
 */
            if ((fpOld == (FILE *)NULL)&&(nextSeqNum<param->numFileSegments)) {
                sprintf(segmentName, "%s.%d", baseDataFileName, nextSeqNum);
                sprintf(newSegName, "%s.%d.new", baseDataFileName, nextSeqNum);
                fpOld = fopen(segmentName, "r");
                fpNew = fopen(newSegName, "w");
                nextSeqNum++;
                if ((fpOld == (FILE *)NULL) || (fpNew == (FILE *)NULL)) { 
                    printf("Open error %d on file %s\n", errno, segmentName);
                    exit(1);
                }
            }

/*
 *          If we have an open data file segment, read lines of
 *          data associated with the next node.
 */
            if (fpOld != (FILE *)NULL) {
                int   i, numNbrs, constraint;
                real8 pos_x, pos_y, pos_z;
                Tag_t tag;

                Getline(inLine, sizeof(inLine), fpOld);

/*
 *              If there is no more data in this file/segment, close
 *              the file and continue
 */
                if (inLine[0] == 0) {
                    fclose(fpOld);
                    fclose(fpNew);
                    fpOld = (FILE *)NULL;
                    fpNew = (FILE *)NULL;
                    continue;
                }

                sscanf(inLine, "%d,%d %lf %lf %lf %d %d",
                       &tag.domainID, &tag.index,
                       &pos_x, &pos_y, &pos_z,
                       &numNbrs, &constraint);

/*
 *              With the version 5 nodal data file, the values for
 *              constraints changed some.  If the data file version
 *              is less than 5, we have to remap the constraint
 *              numbers.
 */
                if (oldFileVersion < 5) {
                    MapV4ToV5Constraint(&constraint);
                }

/*
 *              Write the primary node data line out to the new file.
 */
                fprintf(fpNew, " %d,%d %.8f %.8f %.8f %d %d\n",
                        tag.domainID, tag.index,
                        pos_x, pos_y, pos_z, numNbrs, constraint);

/*
 *              No conversions currently needed for the node's
 *              segment data, so just read in each line of that data
 *              and write it out.
 *
 *              Note: Assumption is that for segmented data files, all the
 *              data for a specific node will be contained in the same
 *              file.
 */
                for (i = 0; i < numNbrs; i++) {
                    Getline(inLine, sizeof(inLine), fpOld);
                    fprintf(fpNew, "%s", inLine);
                    Getline(inLine, sizeof(inLine), fpOld);
                    fprintf(fpNew, "%s", inLine);
                }

            } else {
/*
 *              We completed reading a file and there were no
 *              more data file segments left... so we're done.
 */
                notDone = 0;
            }
        }

        return;
}



/*-----------------------------------------------------------------------
 *
 *      Function:     RenameAllFiles
 *      Description:  On successful conversion, we'll have two copies
 *                    of control and data files, one with a ".new"
 *                    suffix.  Rename all the original files by
 *                    appending a ".bkup" suffix, and rename the
 *                    .new files to the original names.
 *
 *      Arguments:
 *          ctrlFileName     Name of the original control file
 *          baseDataFileName Base name used for the new nodal
 *                           data file(s).
 *          doDataConvert    Flag indicating if a nodal data file had
 *                           been converted and an orignal version exists.
 *
 *----------------------------------------------------------------------*/
static void RenameAllFiles(Home_t *home, char *ctrlFileName,
                           char *baseDataFileName, int doDataConvert)
{
        int  i;
        char bkupFileName[512], origFileName[512], newFileName[512];

/*
 *      Rename the original control file to be a backup and
 *      give the new file the original name.
 */
        sprintf(bkupFileName, "%s.bkup", ctrlFileName);
        sprintf(newFileName, "%s.new", ctrlFileName);

        if (rename(ctrlFileName, bkupFileName) < 0) {
            printf("Error %d renaming %s to %s.  Aborting\n",
                   errno, ctrlFileName, bkupFileName);
            exit(1);
        }
        if (rename(newFileName, ctrlFileName) < 0) {
            printf("Error %d renaming %s to %s.  Aborting\n",
                   errno, newFileName, ctrlFileName);
            exit(1);
        }

/*
 *      Rename each segment of the original nodal data file to be
 *      a backup and give the new segment the original name.
 */
        for (i = 0; i < home->param->numFileSegments; i++) {

            if (home->param->numFileSegments > 1) {
                sprintf(origFileName, "%s.%d", baseDataFileName, i);
                sprintf(bkupFileName, "%s.%d.bkup", baseDataFileName, i);
                sprintf(newFileName, "%s.%d.new", baseDataFileName, i);
            } else {
                sprintf(origFileName, "%s", baseDataFileName);
                sprintf(bkupFileName, "%s.bkup", baseDataFileName);
                sprintf(newFileName, "%s.new", baseDataFileName);
            }

/*
 *          If the control file was an old style which
 *          contained the nodal data within itself, there's no
 *          original data file to rename as a backup.
 */
            if (doDataConvert) {
                if (rename(origFileName, bkupFileName) < 0) {
                    printf("Error %d renaming %s to %s.  Aborting\n",
                           errno, origFileName, bkupFileName);
                    exit(1);
                }
            }

            if (rename(newFileName, origFileName) < 0) {
                printf("Error %d renaming %s to %s.  Aborting\n",
                       errno, newFileName, origFileName);
                exit(1);
            }
        }

        return;
}


int main(int argc, char *argv[])
{
        int         doDataConvert = 1;
        char        *ctrlFileName = (char *)NULL;
        char        *dataFileName = (char *)NULL;
        char        tmpFileName[256];
        char        *start, *sep;
        Param_t     param;
        Home_t      home;

        switch(argc) {
            case 1:
                Usage(argv[0]);
                exit(1);
                break;
            case 2:
                if (strcmp(argv[1], "-help") == 0) {
                    Usage(argv[0]);
                    exit(0);
                }
/*
 *              The data file name has not been provided, so use the
 *              control file name stripped of any suffix.
 */
                ctrlFileName = argv[1];
                strcpy(tmpFileName, ctrlFileName);
                start = strrchr(tmpFileName, '/');
                if (start == (char *)NULL) start = tmpFileName;
                sep = strrchr(start, '.');
                if ((sep != (char *)NULL) && (sep > start)) *sep = 0;
                strcat(tmpFileName, NODEDATA_FILE_SUFFIX);
                dataFileName = &tmpFileName[0];
                break;
            case 3:
                ctrlFileName = argv[1];
                dataFileName = argv[2];
                break;
            default:
                Usage(argv[0]);
                exit(1);
                break;
        }

        memset(&home, 0, sizeof(Home_t));

        home.param = &param;
        home.ctrlParamList = (ParamList_t *)calloc(1, sizeof(ParamList_t));
        home.dataParamList = (ParamList_t *)calloc(1, sizeof(ParamList_t));
        home.numDomains = 1;

/*
 *      If the specified data file name ends in with a
 *      suffix of ".0", assume it is the first of a
 *      set of data file segments and strip the suffix
 *      off to get the base file name.
 */
        start = strrchr(dataFileName, '/');
        if (start == (char *)NULL) start = dataFileName;
        sep = strrchr(start, '.');
        if ((sep != (char *)NULL) && (sep > start)) {
            if (strcmp(sep, ".0") == 0) *sep = 0;
        }

        printf("Control file:    %s\n", ctrlFileName);
        printf("Base data file:  %s\n", dataFileName);

        memset(&param, 0, sizeof(Param_t));

        CtrlParamInit(&param, home.ctrlParamList);
        DataParamInit(&param, home.dataParamList);

        ConvertControlFile(&home, ctrlFileName, dataFileName, &doDataConvert);

        if (doDataConvert) {
            ConvertDataFile(&home, dataFileName);
        }

        RenameAllFiles(&home, ctrlFileName, dataFileName, doDataConvert);

        exit(0);
        return(0);
}
