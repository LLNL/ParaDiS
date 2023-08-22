/**************************************************************************
 *
 *      Module:       ReadBinaryRestart.c
 *      Description:  This module contains the functions for reading
 *                    parameters and nodal data from the binary
 *                    format (HDF) nodal data files indicated (or
 *                    implied) by the command line arguments.
 *
 *      Included public functions:
 *          FindInNode()
 *          ReadBinDataFile()
 *          ReadBinDataParams()
 *          ReadHDFDataset()
 *
 *      Included public functions:
 *
 *************************************************************************/
#include <stdio.h>
#include <sys/stat.h>

#include "mpi_portability.h"

#include "Home.h"
#include "InData.h"
#include "Tag.h"
#include "Util.h"
#include "Decomp.h"
#include "Restart.h"
#include "Parse.h"

#ifdef USE_HDF
#include <hdf5.h>
#endif


#ifdef USE_HDF
/*---------------------------------------------------------------------------
 *
 *      Function:      FindInNode
 *
 *      Description:   Search the array of node structures for the one
 *                     identified by <tag>.
 *
 *                     NOTE: Assumes the <nodelist> array is sorted.
 *
 *      Arguments:
 *          nodelist   pointer to array of <numNodes> node structures
 *          tag        tag of the node for which to search
 *
 *      Last Modified: 09/10/08 - original version
 *
 *-------------------------------------------------------------------------*/
static Node_t *FindInNode(Node_t *nodeList, int numNodes, Tag_t *tag)
{
        Node_t *found;
        Node_t tmpNode;

        tmpNode.myTag = *tag;

        found = (Node_t *)bsearch(&tmpNode, nodeList, numNodes,
                                  sizeof(Node_t), OrderNodes);

        return(found);
}


/*---------------------------------------------------------------------------
 *
 *      Function:      ReadHDFDataset
 *
 *      Description:   Read the specified data item(s) from an HDF file
 *
 *      Arguments:
 *          fileID     File handle for previously opened HDF5 file
 *          datasetName Character string identifying the data item
 *                     to be read from the file.
 *          itemType   HDF5 typy associated with the data to be read
 *          numItems   The number of values the caller expects the
 *                     data item to contain (may be a scalar or an array).
 *          itemList   Location at which the item(s) read are to be stored
 *
 *      Returns:  0 on success
 *                1 if the data could not be read.
 *
 *      Last Modified: 09/10/08 - original version
 *
 *-------------------------------------------------------------------------*/
int ReadHDFDataset(hid_t fileID, char *datasetName, hid_t itemType,
                   hsize_t numItems, void *itemList)
{
        int     status;
        hid_t   dataspaceID, datasetID;
        hsize_t numFileItems;

/*
 *      If the expected item count is zero, don't bother trying
 *      to read the dataset.  (It's easier to have the check
 *      once in here rather than every place we might call the
 *      function.)
 */
        if (numItems < 1) {
            return(0);
        }

/*
 *      Verify the number of items requested is the same as the
 *      number of items in the dataset.
 */
        datasetID = H5Dopen(fileID, datasetName, H5P_DEFAULT);
        dataspaceID = H5Dget_space(datasetID);

        numFileItems = H5Sget_simple_extent_npoints(dataspaceID);

        if (numItems == numFileItems) {
            status = H5Dread(datasetID, itemType, H5S_ALL, H5S_ALL,
                                H5P_DEFAULT, itemList);
        } else {
            status = 1;
        }

        H5Sclose(dataspaceID);
        H5Dclose(datasetID);

        return(status != 0);
}


/*---------------------------------------------------------------------------
 *
 *      Function:      ReadBinDataParams
 *
 *      Description:   Reads the data file parameters (such as decompType,
 *                     nodeCount, etc) from the HDF file.
 *
 *      Arguments:
 *          fileID     File handle for previously opened HDF5 file
 *
 *      Returns:  0 on success
 *                1 on error
 *
 *      Last Modified: 09/10/08 - original version
 *
 *-------------------------------------------------------------------------*/
int ReadBinDataParams(Home_t *home, hid_t fileID)
{
        int         i, maxIndex;
        int         count, type, status, hdfStatus, hdfType;
        char        itemName[128];
        ParamList_t *list;

        list = home->dataParamList;
        maxIndex = list->paramCnt;
        status = 0;

        for (i = 0; i < maxIndex; i++) {
/*
 *          If this parameter is an alias for another, we don't
 *          need to do anything with it.
 */
            if ((list->varList[i].flags & VFLAG_ALIAS) != 0) {
                continue;
            }

            type = list->varList[i].valType;
            count = list->varList[i].valCnt;

/*
 *          For now we're only concerned with the numerical values
 *          stored in the data file.  Skip all others.
 */
            if ((type != V_INT) && (type != V_DBL)) {
                continue;
            }

            switch(type) {
                case V_INT:
                    hdfType = H5T_NATIVE_INT;
                    break;
                case V_DBL:
                    hdfType = H5T_NATIVE_DOUBLE;
                    break;
            }

            sprintf(itemName, "/%s", list->varList[i].varName);
            status |= ReadHDFDataset(fileID, itemName, hdfType,
                                     count, list->varList[i].valList);
        }

        return(status);
}
#endif  /* ifdef USE_HDF */


/*------------------------------------------------------------------------
 *
 *	Function:	ReadBinDataFile
 *	Description:	Read the nodal data from a single or multi-segment
 *                      binary data file, assign the nodes to appropriate
 *                      domains based on the node coordinates and the
 *                      current domain decomposition, and distribute the
 *                      nodal data to the appropriate remote domains.
 *			
 *			The data will be processed in blocks so as to 
 *			avoid loading the entire problem space onto
 *			processor zero and potentially exhausting the
 *			memory on small-memory systems.
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
void ReadBinDataFile(Home_t *home, InData_t *inData, char *dataFile)
{
#ifdef USE_HDF
        int      i, j, dom, iNbr;
        int      numVals, status;
        int      numDomains, numReadTasks;
        int      nextFileSeg, maxFileSeg, segsPerReader, taskIsReader;
        int      distIncomplete, readCount, numNbrs;
        int      binRead = 1;
        int      nextAvailableTag = 0;
        int      fileSegCount = 1, fileSeqNum = 0;
        int      *globalMsgCnt, *localMsgCnt, **nodeLists, *listCounts;
        int      miscData[2];
        int      localNodeCount, globalNodeCount;
        int      doDist, segID;
        int      nodeArrayCnt, firstNodeIndex, neededCnt, nextTaskID;
        int      taskNodeCount, taskSegCount;
        int      taskIDRange[2];
        int      *nodeIndex, *nodeConstraint, *nodeNumSegs, *segTags;
        real8    burgSumX, burgSumY, burgSumZ;
        real8    *nodePos, *burgersVec, *glidePlane;
        char     baseFileName[256], tmpFileName[256];
        char     itemName[64], taskDir[64];
        Node_t   *node;
        Param_t  *param;
        hid_t    fileID;
        Tag_t    tag1, tag2;

        param      = home->param;
        numDomains = home->numDomains;
        fileID     = -1;
        doDist     = 0;

        globalMsgCnt = home->globalReduceBuf;
        localMsgCnt  = home->localReduceBuf;

/*
 *      May need to handle different file versions in the future,
 *      but for now, we only support one version of binary data
 *      files.
 *
 *      Only domain zero reads the initial stuff...
 */
        if (home->myDomain == 0) {

            if (dataFile == (char *)NULL) {
                Fatal("ReadNodeDataFile: No data file provided");
            }

/*
 *          Get the base name of the data file (i.e. no segment
 *          ID appended) and try to open the nodal data file.
 *          If the base file can not be opened, try opening 
 *          a segmented data file (i.e. look for the first
 *          file segment <dataFile>.0). If that can't be opened
 *          either, then exit with an error.
 */
            snprintf(tmpFileName, sizeof(tmpFileName), "%s", dataFile);

            if (!strcmp(&tmpFileName[strlen(tmpFileName)-2], ".0")) {
                    tmpFileName[strlen(tmpFileName)-2] = 0;
            }

            snprintf(baseFileName, sizeof(baseFileName), "%s", tmpFileName);

            struct stat statbuf;
            if (stat(baseFileName, &statbuf) == 0) {
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

/*
 *          Read the basic data file parameters.  All the basic
 *          parameters should exist and be readable.  If not,
 *          something is wrong.
 */
            status = ReadBinDataParams(home, fileID);
            if (status != 0) {
                Fatal("ReadBinDataFile: Error reading basic params from %s",
                      dataFile);
            }

/*
 *          Do a quick sanity check on the decomposition type
 */
            if ((param->dataDecompType < 1) || (param->dataDecompType > 2)) {
                Fatal("dataDecompType=%d is invalid.  Type must be 1 or 2\n",
                      param->dataDecompType);
            }

/*
 *          The minSide{XYZ} and maxSide{XYZ} values are
 *          now redundant but until the rest of the code
 *          is modified to remove references to these
 *          values, we need to explicitly set them now.
 */
            param->minSideX = param->minCoordinates[X];
            param->minSideY = param->minCoordinates[Y];
            param->minSideZ = param->minCoordinates[Z];

            param->maxSideX = param->maxCoordinates[X];
            param->maxSideY = param->maxCoordinates[Y];
            param->maxSideZ = param->maxCoordinates[Z];

/*
 *          Read domain decomposition, node and segment data (or 
 *          whichever items are necessary).
 */
            ReadDecompBounds(home, (void **)&fileID, binRead, param->dataDecompType,
                             &inData->decomp);

/*
 *          Need to set some of the values that are dependent on
 *          the simulation size before we go any further.
 */
            param->Lx = param->maxSideX - param->minSideX;
            param->Ly = param->maxSideY - param->minSideY;
            param->Lz = param->maxSideZ - param->minSideZ;

/*
 *          If the domain decomposition used to create the data file
 *          does not match the domain decomposition requested for this
 *          run, we'll need to create a new uniform decomposition
 *          to start off.
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

/*
 *          Find the min/max ID's of the tasks that wrote nodal
 *          data to this data file.
 */
            sprintf(itemName, "/taskIDRange");
            status = ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                    2, taskIDRange);
                    
            if (status != 0) {
                Fatal("ReadBinDataFile: Error reading %s taskIDRange",
                      baseFileName);
            }

            nextTaskID = taskIDRange[0];

        }  /* if (home->myDomain == 0) */

/*
 *      If there are fewer tasks than data file segments, the 
 *      number of read tasks must be reset to the current task count.
 */
        if (numDomains >= param->numFileSegments) {
            numReadTasks = param->numFileSegments;
        } else {
            numReadTasks = MIN(((param->numFileSegments + (numDomains - 1)) /
                               numDomains), numDomains);
        }

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
 *      the parallel read of nodal data and the number of files segments
 *      comprising the nodal data.
 */
        miscData[0] = numReadTasks;
        miscData[1] = param->numFileSegments;

        MPI_Bcast(miscData, 2, MPI_INT, 0, MPI_COMM_WORLD);

        numReadTasks = miscData[0];
        fileSegCount = miscData[1];

#endif

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
                nodeArrayCnt = MAX_NODES_PER_BLOCK;
                inData->node = (Node_t *)calloc(1, nodeArrayCnt *
                                                sizeof(Node_t));
            }

            while (taskIsReader) {

/*
 *              If the current domain doesn't have a data file segment
 *              opened up yet, open the next segment now.
 */
                if ((fileID < 0) && (nextFileSeg <= maxFileSeg)) {
               
                    sprintf(tmpFileName, "%s.%d", baseFileName, nextFileSeg);
                    fileID = H5Fopen(tmpFileName, H5F_ACC_RDONLY, H5P_DEFAULT);

                    if (fileID < 0) {
                        Fatal("Task %d: Error opening %s", home->myDomain,
                              tmpFileName);
                    }
/*
 *                  Find the min/max ID's of the tasks that wrote nodal
 *                  data to this data file.
 */
                    sprintf(itemName, "/taskIDRange");
                    status = ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                            2, taskIDRange);
                    
                    if (status != 0) {
                        Fatal("ReadBinDataFile: Error reading %s taskIDRange",
                              tmpFileName);
                    }

                    nextTaskID = taskIDRange[0];
                }

/*
 *              If we don't have an open file at this point, it's because
 *              the task has no more data to read, so terminate its
 *              'reader' status.
 */
                if (fileID < 0) {
                    taskIsReader = 0;
                    break;
                }

/*
 *              Loop through the task directories in this file that
 *              we haven't read yet.
 */
                while (nextTaskID <= taskIDRange[1]) {

                    sprintf(taskDir, "/task%d", nextTaskID);
                    sprintf(itemName, "%s/nodeCount", taskDir);
                    status = ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                            1, &taskNodeCount);

                    if (status != 0) {
                        Fatal("ReadBinDataFile: Error reading %s %s",
                              tmpFileName, itemName);
                    }

/*
 *                  If there is no nodal data in the task's directory,
 *                  we don't need to read anything else here.
 */
                    if (taskNodeCount <= 0) {
                        nextTaskID++;
                        continue;
                    }

/*
 *                  Find out how much segment data is stored
 */
                    sprintf(itemName, "%s/segCount", taskDir);
                    status = ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                            1, &taskSegCount);
                    if (status != 0) {
                        Fatal("ReadBinDataFile: Error reading %s %s",
                              tmpFileName, itemName);
                    }

/*
 *                  If the node structure array is non-empty and we're
 *                  about to read in enough nodes to overflow the array,
 *                  break out of the loop now to distribute the nodal
 *                  data.  Otherwise, increase the size of the node
 *                  structure array if necessary, then read in the
 *                  block of nodes.
 */
                    neededCnt = taskNodeCount + readCount;

                    if (neededCnt > nodeArrayCnt) {
                        if (readCount > 0) {
                            doDist = 1;
                            break;
                        } else {
                            inData->node = (Node_t *)realloc(inData->node,
                                                   sizeof(Node_t) * neededCnt);
                            memset(&inData->node[nodeArrayCnt], 0,
                                   sizeof(Node_t) * (neededCnt - nodeArrayCnt));
                            nodeArrayCnt = neededCnt;
                        }
                    } 

/*
 *                  Allocate arrays and read the arrays of nodal and
 *                  segment information from the data file.
 */
                    nodeIndex =
                        (int *)malloc(taskNodeCount * sizeof(int));
                    nodeConstraint =
                        (int *)malloc(taskNodeCount * sizeof(int));
                    nodeNumSegs =
                        (int *)malloc(taskNodeCount * sizeof(int));
                    nodePos =
                        (real8 *)malloc(taskNodeCount * 3 * sizeof(real8));

                    segTags =
                        (int *)malloc(taskSegCount * 3 * sizeof(int));
                    burgersVec =
                        (real8 *)malloc(taskSegCount * 3 * sizeof(real8));
                    glidePlane =
                        (real8 *)malloc(taskSegCount * 3 * sizeof(real8));


                    sprintf(itemName, "%s/nodeIndex", taskDir);
                    status |= ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                             taskNodeCount, nodeIndex);

                    sprintf(itemName, "%s/nodeConstraint", taskDir);
                    status |= ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                             taskNodeCount, nodeConstraint);

                    sprintf(itemName, "%s/nodeNumSegs", taskDir);
                    status |= ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                             taskNodeCount, nodeNumSegs);

                    sprintf(itemName, "%s/nodePos", taskDir);
                    status |= ReadHDFDataset(fileID, itemName,
                                             H5T_NATIVE_DOUBLE,
                                             taskNodeCount * 3, nodePos);

                    sprintf(itemName, "%s/segTags", taskDir);
                    status |= ReadHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                             taskSegCount * 3, segTags);

                    sprintf(itemName, "%s/burgersVec", taskDir);
                    status |= ReadHDFDataset(fileID, itemName,
                                             H5T_NATIVE_DOUBLE,
                                             taskSegCount * 3, burgersVec);

                    sprintf(itemName, "%s/glidePlane", taskDir);
                    status |= ReadHDFDataset(fileID, itemName,
                                             H5T_NATIVE_DOUBLE,
                                             taskSegCount * 3, glidePlane);

                    if (status != 0) {
                        Fatal("ReadBinDataFile: Error reading %s in %s",
                              taskDir, tmpFileName);
                    }

/*
 *                  Loop through the nodes just read in, and populate
 *                  the node structures with the data.
 */
                    firstNodeIndex = readCount;

                    for (i = 0; i < taskNodeCount; i++) {

                        node = &inData->node[readCount++];

                        node->myTag.domainID = nextTaskID;
                        node->myTag.index = nodeIndex[i];

                        node->constraint = nodeConstraint[i];

/*
 *                      With the version 5 nodal data file, the values for
 *                      constraints changed some.  If the data file version
 *                      is less than 5, we have to remap the constraint
 *                      numbers.
 */
                        if (param->dataFileVersion < 5) {
                            MapV4ToV5Constraint(&node->constraint);
                        }

                        node->x = nodePos[i*3  ];
                        node->y = nodePos[i*3+1];
                        node->z = nodePos[i*3+2];

                        FoldBox(param, &node->x, &node->y, &node->z);
                    
                        numNbrs = nodeNumSegs[i];

                        AllocNodeArms(node, numNbrs);
/*
 *                      We've allocated the correct number of
 *                      arms for the node, but we're going to zero
 *                      the value out for now.  The value will
 *                      be incremented later as segments get
 *                      processed and the data is associated with
 *                      the proper nodes.  After all segments
 *                      have been processed, the node numNbrs
 *                      values should have been incremented back
 *                      to the correct value.
 */
                        node->numNbrs = 0;
                    }

/*
 *                  Loop though all the segments just read in, locate
 *                  the local nodes associated with these segments and
 *                  add the segment data to the node structures.
 */
                    for (i = 0; i < taskSegCount; i++) {

                        tag1.domainID = nextTaskID;
                        tag1.index    = segTags[i*3  ];
                        tag2.domainID = segTags[i*3+1];
                        tag2.index    = segTags[i*3+2];

/*
 *                      The endpoint specified by tag1 is local
 *                      to the task whose nodes were just read
 *                      so the corresponding node always exists.
 */
                        node = FindInNode(&inData->node[firstNodeIndex],
                                          taskNodeCount, &tag1);

                        segID = node->numNbrs;

                        node->nbrTag[segID].domainID = tag2.domainID;
                        node->nbrTag[segID].index    = tag2.index;

                        node->burgX[segID] = burgersVec[i*3  ];
                        node->burgY[segID] = burgersVec[i*3+1];
                        node->burgZ[segID] = burgersVec[i*3+2];

                        node->nx[segID] = glidePlane[i*3  ];
			node->ny[segID] = glidePlane[i*3+1];
                        node->nz[segID] = glidePlane[i*3+2];

                        Normalize(&node->nx[segID], &node->ny[segID],
                                  &node->nz[segID]);

                        node->numNbrs += 1;

/*
 *                      If the 2nd endpoint is in the same domain
 *                      it to must exist in the array, so update
 *                      its segment list also.  Note that the
 *                      burgers vector will be opposite from
 *                      this side of the segment!
 */
                        if (tag2.domainID != nextTaskID) {
                            continue;
                        }

                        node = FindInNode(&inData->node[firstNodeIndex],
                                          taskNodeCount, &tag2);

                        segID = node->numNbrs;

                        node->nbrTag[segID].domainID = tag1.domainID;
                        node->nbrTag[segID].index    = tag1.index;

                        node->burgX[segID] = -burgersVec[i*3  ];
                        node->burgY[segID] = -burgersVec[i*3+1];
                        node->burgZ[segID] = -burgersVec[i*3+2];

                        node->nx[segID] = glidePlane[i*3  ];
			node->ny[segID] = glidePlane[i*3+1];
                        node->nz[segID] = glidePlane[i*3+2];

                        Normalize(&node->nx[segID], &node->ny[segID],
                                  &node->nz[segID]);

                        node->numNbrs += 1;
                    }

/*
 *                  Free the array buffers since we don't know if we'll
 *                  need them anymore or what sizes they would be.
 */
                    free(nodeIndex);       nodeIndex = (int *)NULL;
                    free(nodeConstraint);  nodeConstraint = (int *)NULL;
                    free(nodeNumSegs);     nodeNumSegs = (int *)NULL;
                    free(nodePos);         nodePos = (real8 *)NULL;

                    free(segTags);         segTags = (int *)NULL;
                    free(burgersVec);      burgersVec = (real8 *)NULL;
                    free(glidePlane);      glidePlane = (real8 *)NULL;

                    nextTaskID++;

                }  /* end while (nextTaskID <= maxTaskID)  */

                if (nextTaskID > taskIDRange[1]) {
                    H5Fclose(fileID);
                    fileID = -1;
                    nextFileSeg++;
                }

/*
 *              If we need to distribute nodal data to remote domains
 *              break out of the current loop.
 */
                if (doDist) {
                    doDist = 0;
                    break;
                }

            }  /* end if (taskIsReader) */

/*
 *          Do a run through all the nodes just read in and verify that
 *          the burgers vector is conserved for all nodes that are not
 *          pinned or surface nodes.
 */
            for (i = 0; i < readCount; i++) {

                node = &inData->node[i];

                if (HAS_NONE_OF_CONSTRAINTS(node->constraint, PINNED_NODE |
                                            SURFACE_NODE)) {

                    burgSumX = 0.0;
                    burgSumY = 0.0;
                    burgSumZ = 0.0;

                    for (j = 0; j < node->numNbrs; j++) {
                        burgSumX += node->burgX[j];
                        burgSumY += node->burgY[j];
                        burgSumZ += node->burgZ[j];
                    }

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
            }

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

#endif  /* ifdef USE_HDF */

        return;
}
