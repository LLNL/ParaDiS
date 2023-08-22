/**************************************************************************
 *
 *      Module:       ReadRestart.c
 *      Description:  This module contains the functions for reading
 *                    parameters and nodal data from the control and
 *                    nodal data files indicated (or implied) by the
 *                    command line arguments.
 *
 *                    Several versions of the nodal data file are
 *                    supported, so there are multiple functions
 *                    provided to deal with different formats.
 *
 *                    NOTE:  The paradisconvert utility can be used
 *                    to translate any old style control and nodal
 *                    data files to the most current version.
 *
 *      Included public functions:
 *          AssignNodesToDomains()
 *          FreeNodeLists()
 *          ReadControlFile()
 *          ReadNodeDataFile()
 *
 *      Included private functions:
 *          ReadPreV4DataParams()
 *
 *************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "InData.h"
#include "Tag.h"
#include "Util.h"
#include "Decomp.h"
#include "Restart.h"
#include "Parse.h"


/*---------------------------------------------------------------------------
 *
 *      Function:       ReadControlFile
 *      Description:    Read and parse the contents of the control
 *                      parameter file, saving the values associated
 *                      with the parameters in the corresponding
 *                      variables.
 *
 *      Arguments:
 *          ctrlFileName  Name of the control file to be read.
 *
 *-------------------------------------------------------------------------*/
void ReadControlFile(Home_t *home, char *ctrlFileName)
{
        int  maxTokenLen, tokenType, pIndex, okay;
        int  valType, numVals;
        char token[256];
        void *valList;
        FILE *fpCtrl;

        maxTokenLen = sizeof(token);
        tokenType = TOKEN_GENERIC;

        if ((fpCtrl = fopen(ctrlFileName, "r")) == (FILE *)NULL) {
            Fatal("ReadControlFile: Error %d opening file %s",
                  errno, ctrlFileName);
        }

        while (1) {

/*
 *          Next token (if any) should be a parameter name...
 */
            tokenType = GetNextToken(fpCtrl, token, maxTokenLen);

/*
 *          If we hit the end of the file (i.e. no more tokens)
 *          we're done...
 */
            if (tokenType == TOKEN_NULL) {
                break;
            }

/*
 *          If there were any parsing errors, just abort
 */
            if (tokenType == TOKEN_ERR) {
                Fatal("ReadControlFile: Error parsing file %s", ctrlFileName);
            }

/*
 *          Obtain values associated with the parameter; don't
 *          save values for unknown/obsolete parameters.
 */
            pIndex = LookupParam(home->ctrlParamList, token);

            if (pIndex < 0) {
                printf("Ignoring unrecognized parameter (%s)\n", token);
                valType = V_NULL;
                numVals = 0;
                valList = (void *)NULL;
            } else {
                valType = home->ctrlParamList->varList[pIndex].valType;
                numVals = home->ctrlParamList->varList[pIndex].valCnt;
                valList = home->ctrlParamList->varList[pIndex].valList;
                home->ctrlParamList->varList[pIndex].flags |= VFLAG_SET_BY_USER;
            }

            okay = GetParamVals(fpCtrl, valType, numVals, valList);
            if (!okay) {
                Fatal("Error obtaining values for parameter %s",
                       home->ctrlParamList->varList[pIndex].varName);
            }
        }

        fclose(fpCtrl);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       AssignNodesToDomains
 *      Description:    Using the provided domain decomposition data,
 *                      loop though all nodes in the inData node array
 *                      and assign each node to the domain which encompasses
 *                      the node's coordinates.
 *
 *      Arguments:
 *          nodeCount   number of nodes in the inData node array.
 *          nodeLists   Location in which to return to the caller the
 *                      array of lists of nodes to be sent to the
 *                      remote domains.  The nodeLists array has
 *                      a list pointer for every domain in the problem.
 *                      (each list is an array of indices into the
 *                      inData node array)
 *          listCounts  Location in which to return to the caller an
 *                      array of integers indicating the count of nodes
 *                      on the corresponding list in <nodeLists>.
 *
 *-------------------------------------------------------------------------*/
void AssignNodesToDomains(Home_t *home, InData_t *inData, int nodeCount,
                          int ***nodeLists, int **listCounts)
{
        int      i, j;
        int      nXdoms, nYdoms, nZdoms;
        int      domIndex, nDoms, len, nexti;
        int      *qHead, *qCount, *list;
        int      **listArray, *countArray;
        Param_t  *param;
        Node_t   *node;


        if (nodeCount == 0) {
            *nodeLists = (int **)NULL;
            *listCounts = (int *)NULL;
            return;
        }

        param = home->param;

        nXdoms = param->nXdoms;
        nYdoms = param->nYdoms;
        nZdoms = param->nZdoms;

        nDoms = nXdoms * nYdoms * nZdoms;

        listArray  = (int **) calloc(1, nDoms * sizeof(int *));
        countArray = (int *) calloc(1, nDoms * sizeof(int));

/*
 *      Allocate and initialize an array of integers (1 per domain)
 *      to be used as pointers to the head of the queue of nodes
 *      assigned to the domains.
 */
        qHead  = (int *)malloc(nDoms * sizeof(int));
        qCount = (int *)malloc(nDoms * sizeof(int));

        for (i = 0; i < nDoms; i++) {
            qHead[i] = -1;
            qCount[i] = 0;
        }

/*
 *      Loop through all the nodes on the current node list, find
 *      the proper domain for the node based on the node coordinates
 *      and add the node to the domain's queue.
 */
        for (i = 0; i < nodeCount; i++) {

            node = &inData->node[i];
            domIndex = FindCoordDomain(home, 0, &node->x, &node->y, &node->z);

/*
 *          Each node on a queue contains the index of the next node
 *          on the queue.  To add a node to the queue, we just set
 *          that node's 'next node' pointer to the current head of
 *          the queue, and then point the queue head to the new node.
 */
            nexti = qHead[domIndex];
            if (nexti < 0) {
                node->next = (Node_t *)NULL;
            } else {
                node->next = &inData->node[qHead[domIndex]];
            }
            qHead[domIndex] = i;
            qCount[domIndex]++;
        }

/*
 *      For each domain, generate the list of indices in the node array
 *      for nodes to be sent to the domain.
 */
        for (domIndex = 0; domIndex < nDoms; domIndex++) {

            list  = (int *)NULL;
            len   = qCount[domIndex];
            nexti = qHead[domIndex];

            if (len > 0) list = (int *)malloc(len * sizeof(int));

            for (j = 0; j < len; j++) {
                list[j] = nexti;
/*
 *              Do some pointer arithmetic to get convert pointer
 *              addresses into array indices.
 */
                if (inData->node[nexti].next == (Node_t *)NULL) {
                    nexti = -1;
                } else {
                    nexti = inData->node[nexti].next - inData->node;
                }
            }

            if (nexti != -1) {
                Fatal("Queue error. domain %d, queue len %d", domIndex, len);
            }

            listArray[domIndex] = list;
            countArray[domIndex] = len;
        }

/*
 *      Free up the now unneeded arrays.  (the lists of node indices
 *      will be freed elsewhere when they are no longer needed)
 */
        free(qHead);
        free(qCount);

        *nodeLists  = listArray;
        *listCounts = countArray;

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       FreeNodeLists
 *      Description:    Releases temporary storage associated with 
 *                      the arrays of node indices identifying nodes
 *                      to be sent to remote domains.
 *
 *      Arguments:
 *          nodeLists   Address of the the array of lists of nodes to
 *                      that were sent to the remote domains.  The nodeLists
 *                      array has a list pointer for every domain in the
 *                      problem.  On return, the contents of this adress
 *                      will be zeroed.
 *          listCounts  Address of the array of integers indicating the
 *                      count of nodes on the corresponding list in
 *                      <nodeLists>.  On return, the contents of this
 *                      adress will be zeroed.
 *
 *-------------------------------------------------------------------------*/
void FreeNodeLists(Home_t *home, int ***nodeLists, int **listCounts)
{
        int dom;

        if (*listCounts != (int *)NULL) {
            free(*listCounts);
            *listCounts = (int *)NULL;
        }

        if (*nodeLists != (int **)NULL) {
            for (dom = 0; dom < home->numDomains; dom++) {
                if ((*nodeLists)[dom] != (int *)NULL) {
                    free((*nodeLists)[dom]);
                    (*nodeLists)[dom] = (int *)NULL;
                }
            }
            free(*nodeLists);
            *nodeLists = (int **)NULL;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ReadPreV4DataParams
 *      Description:    Read all parameters and the domain decomposition
 *                      (if present) from the opened pre-version-4 nodal
 *                      data file.  On exit from this function, the file
 *                      pointer will be positioned at the start of the
 *                      raw nodal data.
 *
 *      Arguments:
 *          fp          File stream associated with the nodal data 
 *                      file to be read.
 *          dataDecomp  Pointer to location in which to return a pointer
 *                      to the domain decomposition (if any) read in
 *                      from the file.  
 *
 *-------------------------------------------------------------------------*/
void ReadPreV4DataParams(Home_t *home, FILE *fp, void **dataDecomp)
{
        int     binRead = 0;
        char    inLine[512];
        void    *oldDecomp;
        Param_t *param;

        param = home->param;
        oldDecomp = (void *)NULL;

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
 *      Get the domain geometry, and (if present) the domain
 *      decomposition type
 */
        Getline(inLine, sizeof(inLine), fp);
        sscanf(inLine, "%d %d %d", &param->dataDecompGeometry[0],
               &param->dataDecompGeometry[1], &param->dataDecompGeometry[2]);

        if (param->dataFileVersion > 2) {
            Getline(inLine, sizeof(inLine), fp);
            sscanf(inLine, "%d", &param->dataDecompType);
        } else {
            param->dataDecompType = 1;
        }

/*
 *      The minSide{XYZ} and maxSide{XYZ} values are now redundant
 *      but until the rest of the code is modified to remove references
 *      to these values, we need to explicitly set them here.
 */
        param->minSideX = param->minCoordinates[X];
        param->minSideY = param->minCoordinates[Y];
        param->minSideZ = param->minCoordinates[Z];

        param->maxSideX = param->maxCoordinates[X];
        param->maxSideY = param->maxCoordinates[Y];
        param->maxSideZ = param->maxCoordinates[Z];

/*
 *      Read domain decomposition (if present).  If the type of domain
 *      decomposition in the file does not match the domain decomp
 *      type for this run, the decomposition from the file will
 *      simply be read in and discarded and not returned by the
 *      following call.  In that case, a uniform decomposition
 *      of the proper type will have to be generated later.
 */
        ReadDecompBounds(home, (void **)&fp, binRead, param->decompType,
                         &oldDecomp);

/*
 *      For version 1 files, the total node count is located
 *      next in the file, so read it in now.
 */
        if (param->dataFileVersion == 1) {
            Getline(inLine, sizeof(inLine), fp);
            sscanf(inLine, "%d", &param->nodeCount);
        }

        *dataDecomp = oldDecomp;

        return;
}


/*------------------------------------------------------------------------
 *
 *	Function:	MapV4ToV5Constraint
 *	Description:	Nodal constraint values changed between versions
 *	                4 and 5 of the nodal data files, so for any
 *	                pre-v5 files, we need to remap nodal contraints
 *	                from the old values to the corresponding new values.
 *
 *-----------------------------------------------------------------------*/
void MapV4ToV5Constraint(int *constraint)
{

/*
 *      Old SURFACE_NODE constraint value was 1; remap to current value
 */
        if (*constraint == 1) {
            SET_CONSTRAINTS(*constraint, SURFACE_NODE);
        }

/*
 *      Old PINNED_NODE constraint value was 7; remap to current value
 */
        if (*constraint == 7) {
            SET_CONSTRAINTS(*constraint, PINNED_NODE);
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *	Function:	ReadNodeDataFile
 *	Description:	Read the nodal data from a single or multi-segment
 *                      data file, assign the nodes to appropriate domains
 *                      based on the node coordinates and the current domain
 *			decomposition, and distribute the nodal data
 *			to the appropriate remote domains.
 *			
 *			The data will be processed in blocks so as to 
 *			avoid loading the entire problem space onto
 *			processor zero and potentially exhausting the
 *			memory on small-memory systems.
 *
 *			NOTE: This function is only for use with newer
 *			data files, not for reading node data provided
 *			in the old format within the control file itself!
 *
 *      Arguments:
 *          inData    pointer to structure in which to temporarily
 *                    store domain decomposition, nodal data, etc
 *                    from the restart file.
 *          dataFile  Name of the nodal data file.  For segmented
 *                    restart files, this may be the base name
 *                    (i.e. no trailing sequence number) or the
 *                    name of the first file segment.
 *
 *-----------------------------------------------------------------------*/
void ReadNodeDataFile(Home_t *home, InData_t *inData, char *dataFile)
{
        int      i, dom, iNbr;
        int      maxTokenLen, tokenType, pIndex, okay;
        int      valType, numVals;
        int      numDomains, numReadTasks;
        int      nextFileSeg, maxFileSeg, segsPerReader, taskIsReader;
        int      distIncomplete, readCount, numNbrs;
        int      fileSegCount, fileSeqNum = 0;
        int      *globalMsgCnt, *localMsgCnt, **nodeLists, *listCounts;
        int      miscData[2]={0,0};
        int      localNodeCount, globalNodeCount;
        int      binRead = 0;
        int      nextAvailableTag = 0;
        real8    burgSumX, burgSumY, burgSumZ;
        void     *valList;
        char     inLine[500], token[256];
        char     baseFileName[256], tmpFileName[256];
        FILE     *fpSeg;
        Node_t   *node;
        Param_t  *param;
        ParamList_t *dataParamList;

        param      = home->param;
        numDomains = home->numDomains;
        fpSeg      = (FILE *)NULL;

        dataParamList = home->dataParamList;

        globalMsgCnt = home->globalReduceBuf;
        localMsgCnt  = home->localReduceBuf;

        memset(inLine, 0, sizeof(inLine));
        maxTokenLen = sizeof(token);

/*
 *      Only domain zero reads the initial stuff...
 */
        if (home->myDomain == 0) {

            if (dataFile == (char *)NULL) {
                Fatal("ReadNodeDataFile: No data file provided");
            }

            snprintf(tmpFileName, sizeof(tmpFileName), "%s", dataFile);

            if (!strcmp(&tmpFileName[strlen(tmpFileName)-2], ".0")) {
                    tmpFileName[strlen(tmpFileName)-2] = 0;
            }

            snprintf(baseFileName, sizeof(baseFileName), "%s", tmpFileName);

/*
 *          Try to open the nodal data file.  If the specified file
 *          can not be opened, try looking for a segmented data file
 *          (i.e. look for the first file segment <dataFile>.0). If
 *          that can't be opened either, then exit with an error.
 */

            if ((fpSeg = fopen(tmpFileName, "r")) == (FILE *)NULL) {
                snprintf(tmpFileName, sizeof(tmpFileName), "%s.%d",
                         baseFileName, fileSeqNum);
                if ((fpSeg = fopen(tmpFileName, "r")) == (FILE *)NULL) {
                    Fatal("Error %d opening file %s to read nodal data",
                          errno, dataFile);
                }
            }

/*
 *          Get the first token.  This should either be a known
 *          parameter identifier, or a file version number.
 */
            tokenType = GetNextToken(fpSeg, token, maxTokenLen);
            pIndex = LookupParam(dataParamList, token);

            if (pIndex < 0) {
/*
 *              If the token does not correspond to a known parameter, it
 *              should be the version number, so convert it to an integer
 *              and parse the rest of the file the old way (i.e. values
 *              are positional, not identified by name/value pairs)
 */
                param->dataFileVersion = atoi(token);

                if ((param->dataFileVersion < 1) ||
                    (param->dataFileVersion > 3)) {
                    Fatal("ReadNodeDatFile: Unsupported file version %d",
                          param->dataFileVersion);
                }

                ReadPreV4DataParams(home, fpSeg, &inData->decomp);

            } else {
/*
 *              Just go through the nodal data file reading all
 *              the associated parameters.  Need to do special
 *              processing of the domain decomposition, and when
 *              when we hit the 'nodaldata' identifier, just break
 *              out of this loop.
 */
                while ((tokenType != TOKEN_ERR) && (tokenType != TOKEN_NULL)) {

                    if (pIndex >= 0) {
/*
 *                      Token represents a known parameter identifier, so
 *                      read the associated value(s).
 */
                        valType = dataParamList->varList[pIndex].valType;
                        numVals = dataParamList->varList[pIndex].valCnt;
                        valList = dataParamList->varList[pIndex].valList;
                        okay = GetParamVals(fpSeg, valType, numVals, valList);
                        if (!okay) {
                            Fatal("Parsing Error obtaining values for "
                                  "parameter %s\n",
                                  dataParamList->varList[pIndex].varName);
                        }
    
                    } else {
/*
 *                      Token does not represent one of the simple 
 *                      parameters.  If it's not one of the identifiers
 *                      that needs special handling, skip it.
 */
                        if (strcmp(token, "domainDecomposition") == 0) {
/*
 *                          The minSide{XYZ} and maxSide{XYZ} values are
 *                          now redundant but until the rest of the code
 *                          is modified to remove references to these
 *                          values, we need to explicitly set them now.
 */
                            param->minSideX = param->minCoordinates[X];
                            param->minSideY = param->minCoordinates[Y];
                            param->minSideZ = param->minCoordinates[Z];

                            param->maxSideX = param->maxCoordinates[X];
                            param->maxSideY = param->maxCoordinates[Y];
                            param->maxSideZ = param->maxCoordinates[Z];

                            tokenType = GetNextToken(fpSeg, token, maxTokenLen);
/*
 *                          Do a quick verification of the decomposition type
 */
                            if ((param->dataDecompType < 1) ||
                                (param->dataDecompType > 2)) {
                                Fatal("dataDecompType=%d is invalid.  Type must be 1 or 2\n",
                                      param->dataDecompType);
                            }

                            ReadDecompBounds(home, (void **)&fpSeg, binRead,
                                             param->decompType,
                                             &inData->decomp);

                        } else if (strcmp(token, "nodalData") == 0) {
/*
 *                          When we hit the nodal data, we can break
 *                          out of the loop because we are assuming
 *                          all other data file parameters have been
 *                          processed.  If they have not, we have a
 *                          problem since processing of the nodal data
 *                          requires the other parameters.
 *                          Note: Remainder of the file should just
 *                          contain " = " followed by the nodal data,
 *                          so be sure to skip the next token before
 *                          processing the nodal data.
 */
                            tokenType = GetNextToken(fpSeg, token, maxTokenLen);
                            break;
                        } else {
/*
 *                          If the parameter is not recognized, skip the
 *                          parameter and any associated value(s).
 */
                            printf("Ignoring unknown data file parameter %s\n",
                                   token);
                            valType = V_NULL;
                            numVals = 0;
                            valList = (void *)NULL;
                            okay = GetParamVals(fpSeg, valType, numVals,
                                                valList);
                        }
                    }

                    tokenType = GetNextToken(fpSeg, token, maxTokenLen);

                    if ((tokenType == TOKEN_NULL)||(tokenType == TOKEN_ERR)) {
                        Fatal("Parsing error on file %s\n", tmpFileName);
                    }
                    pIndex = LookupParam(dataParamList, token);
                }
            }

/*
 *          Need to set some of the values that are dependent on
 *          the simulation size before we go any further.
 */
            param->Lx = param->maxSideX - param->minSideX;
            param->Ly = param->maxSideY - param->minSideY;
            param->Lz = param->maxSideZ - param->minSideZ;

/*
 *          If we did not get a domain decomposition from the restart
 *          file (whether because of a mismatch in domain geometry or
 *          domain decomposition type between the current run and
 *          that from the restart file) we'll need to create a new
 *          uniform decomposition to start off.
 */
            if (inData->decomp == (void *)NULL) {
                printf("Generating uniform domain decomposition.\n");
                UniformDecomp(home, &inData->decomp);
/*
 *              Since the current domain geometry does not match the
 *              one from the restart file, some of the nodes from the
 *              restart file will be assigned new node IDs.  Given
 *              that, we won't try to preserve *any* of the node IDs.
 *              This is so we don't end up with very sparsely populated
 *              nodeKeys arrays which are bad when the code is threaded.
 */
                param->assignNewNodeTags = 1;
            }

            if (numDomains >= param->numFileSegments) {
                numReadTasks = param->numFileSegments;
            } else {
                numReadTasks = MIN(((param->numFileSegments + (numDomains - 1)) /
                                   numDomains), numDomains);
            }

            miscData[0] = numReadTasks;
            miscData[1] = param->numFileSegments;

        }  /* if (home->myDomain == 0) */

#ifdef PARALLEL
/*
 *      Domain zero now needs to pass various pieces of data to
 *      the remote domains before any processes begin to read
 *      and distribute the nodal data... start with the param structure.
 */
        MPI_Bcast((char *)param, sizeof(Param_t), MPI_CHAR, 0, MPI_COMM_WORLD);

/*
 *      Now the base nodal data file name.
 */
        MPI_Bcast((char *)baseFileName, sizeof(baseFileName),
                  MPI_CHAR, 0, MPI_COMM_WORLD);

/*
 *      Also send off the count of domains that will be involved in
 *      the parallel read of nodal data and the number of file segments
 *      comprising the nodal data.
 */
        MPI_Bcast(miscData, 2, MPI_INT, 0, MPI_COMM_WORLD);
#endif

        numReadTasks = miscData[0];
        fileSegCount = miscData[1];

/*
 *      Lastly, only domain zero knows the current domain decomposition.
 *      Invoke a function to distribute that data to remote domains (if
 *      any).
 */
        BroadcastDecomp(home, inData->decomp);

/*
 *      Have each domain determine which (if any) of the
 *      nodal data file segments it is to read in.
 */
        segsPerReader = (fileSegCount + (numReadTasks - 1)) / numReadTasks;
        nextFileSeg   = home->myDomain * segsPerReader;
        maxFileSeg    = MIN((nextFileSeg+(segsPerReader-1)),(fileSegCount-1));
        taskIsReader  = (nextFileSeg <= maxFileSeg);

/*
 *      All processes loop until all nodal data has been read in
 *      and distributed to the appropriate domains.
 */
        distIncomplete = 1;
        readCount = 0;

        while (distIncomplete) {

/*
 *          Have each reader task allocate a buffer for the next
 *          block of nodes to be read in.  Then have the readers
 *          read their next blocks of nodes.
 */
            if (taskIsReader) {
                inData->node = (Node_t *)calloc(1, MAX_NODES_PER_BLOCK *
                                                sizeof(Node_t));
            }

            while (taskIsReader) {

/*
 *              If the current domain doesn't have a data file segment
 *              opened up yet, open the next segment now.
 */
                if ((fpSeg == (FILE *)NULL) && (nextFileSeg <= maxFileSeg)) {
               
                    sprintf(tmpFileName, "%s.%d", baseFileName, nextFileSeg);
                    fpSeg = fopen(tmpFileName, "r");

                    if (fpSeg == (FILE *)NULL) {
                        Fatal("Task %d: Error %d opening %s", home->myDomain,
                              errno, tmpFileName);
                    }
                    
                }

/*
 *              If we have a valid open file pointer now, continue
 *              reading nodal data.  If we don't, it's because the
 *              task has no more data to read so terminate its
 *              'reader' status.
 */
                if (fpSeg != (FILE *)NULL) {

                    Getline(inLine, sizeof(inLine), fpSeg);

/*
 *                  If we hit EOF, close the current segment and
 *                  loop back to open the next segment (if any)
 *                  this domain is responsible for reading.
 */
                    if (inLine[0] == 0) {
                        nextFileSeg++;
                        fclose(fpSeg);
                        fpSeg = (FILE *)NULL;
                        continue;  
                    }

/*
 *                  Got information on the next node, so deal with it.
 */
                    node = &inData->node[readCount++];

                    sscanf(inLine, "%d,%d %lf %lf %lf %d %d",
                           &node->myTag.domainID, &node->myTag.index,
                           &node->x, &node->y, &node->z,
                           &numNbrs, &node->constraint);

/*
 *                  With the version 5 nodal data file, the values for
 *                  constraints changed some.  If the data file version
 *                  is less than 5, we have to remap the constraint
 *                  numbers.
 */
                    if (param->dataFileVersion < 5) {
                        MapV4ToV5Constraint(&node->constraint);
                    }

                    AllocNodeArms(node, numNbrs); 
/*
 *                  Read in data for each arm of the node
 */
                    burgSumX = 0.0;
                    burgSumY = 0.0;
                    burgSumZ = 0.0;

                    for (iNbr = 0; iNbr < node->numNbrs; iNbr++) {

                        Getline(inLine, sizeof(inLine), fpSeg); 
                        sscanf(inLine, "%d,%d %lf %lf %lf",
                               &node->nbrTag[iNbr].domainID,
                               &node->nbrTag[iNbr].index,
                               &node->burgX[iNbr],
                               &node->burgY[iNbr],
                               &node->burgZ[iNbr]);

                        Getline(inLine, sizeof(inLine), fpSeg);
                        sscanf(inLine, "%lf %lf %lf", &node->nx[iNbr],
                               &node->ny[iNbr], &node->nz[iNbr]);

                        Normalize(&node->nx[iNbr], &node->ny[iNbr], &node->nz[iNbr]);

                        burgSumX += node->burgX[iNbr];
                        burgSumY += node->burgY[iNbr];
                        burgSumZ += node->burgZ[iNbr];
                    }

                    FoldBox(param, &node->x, &node->y, &node->z);

/*
 *                  Just a quick sanity check, to make sure burgers
 *                  vector is conserved for all nodes that are neither
 *                  pinned nor surface nodes.
 */
                    if (HAS_NONE_OF_CONSTRAINTS(node->constraint, PINNED_NODE |
                                                SURFACE_NODE)) {

                        if ((fabs(burgSumX) > 0.0001) ||
                            (fabs(burgSumY) > 0.0001) ||
                            (fabs(burgSumZ) > 0.0001)) {

                            printf("Error: node (%d,%d)\n",
                                   node->myTag.domainID, node->myTag.index);

                            for (iNbr=0; iNbr < node->numNbrs; iNbr++) {
                                printf("  arm[%d] burg = %e %e %e\n",
                                       iNbr, node->burgX[iNbr],
                                       node->burgY[iNbr],
                                       node->burgZ[iNbr]);
                            }

                            Fatal("Burger's vector not conserved!");
                        }
                    }
                } else {
/*
 *                  No more files to be read by this task
 */
                    taskIsReader = 0;
                }

/*
 *              We read the nodal data in blocks rather than all at once.
 *              If this domain has read in the maximum number of nodes
 *              allowed in a block, stop reading data and get ready
 *              to distribute it to the appropriate remote domains.
 */
                if (readCount >= MAX_NODES_PER_BLOCK) {
                    break;
                }

            }  /* while (taskIsReader) */

/*
 *          Determine the domains to which to send any nodes
 *          read in by this domain.
 */
            AssignNodesToDomains(home, inData, readCount, &nodeLists,
                                 &listCounts);

/*
 *          Set up an array (1 element per domain).  Each reader
 *          task sets to 1 the array entry for each remote domain
 *          to which it will be sending data.  When we do a global
 *          reduction to sum up the arrays, the resulting array
 *          contains the number of messages that will be sent to
 *          each domain during this communication.
 *
 *          Plus 1 extra element set to 1 if ANY process is sending
 */
#ifdef PARALLEL

            memset(localMsgCnt, 0, (numDomains+1) * sizeof(int));

            if (listCounts != (int *)NULL) {
                for (dom = 0; dom < numDomains; dom++) {
                    if (listCounts[dom] > 0) {
                       localMsgCnt[dom] = 1;
                       localMsgCnt[numDomains] = 1;
                    }
                }
            }

            MPI_Allreduce(localMsgCnt, globalMsgCnt, numDomains + 1,
                          MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#else
            memset(globalMsgCnt, 0, (numDomains+1) * sizeof(int));

            if (listCounts != (int *)NULL) {
                for (dom = 0; dom < numDomains; dom++) {
                    if (listCounts[dom] > 0) {
                       globalMsgCnt[dom] = 1;
                       globalMsgCnt[numDomains] = 1;
                    }
                }
            }
#endif

/*
 *          If no domains have any more data to send out, we're
 *          done reading/distributing the data.
 */
            if (globalMsgCnt[numDomains] == 0) {
                distIncomplete = 0;
                FreeInNodeArray(inData, readCount);
                FreeNodeLists(home, &nodeLists, &listCounts);
                continue;
            }

/*
 *          Do next send/receive of nodal data
 */
            SendInitialNodeData(home, inData, globalMsgCnt, nodeLists,
                                listCounts, &nextAvailableTag);

            FreeInNodeArray(inData, readCount);
            FreeNodeLists(home, &nodeLists, &listCounts);
            readCount = 0;

        }  /* if (distIncomplete) */

/*
 *      This is a good place for a quick sanity check that the sum of
 *      nodes on all domains equals the total node count from the
 *      data file.
 */
        localNodeCount = 0;
        globalNodeCount = 0;

        for (i = 0; i < home->newNodeKeyPtr; i++) {
            if (home->nodeKeys[i] != (Node_t *)NULL) {
                localNodeCount++;
            }
        }

#ifdef PARALLEL
        MPI_Reduce(&localNodeCount, &globalNodeCount, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);

#else
        globalNodeCount = localNodeCount;
#endif

        if ((home->myDomain == 0) && (param->nodeCount != globalNodeCount)) {
            Fatal("ReadNodedataFile: Read %d nodes, expected %d!",
                  globalNodeCount, param->nodeCount);
        }

        return;
}
