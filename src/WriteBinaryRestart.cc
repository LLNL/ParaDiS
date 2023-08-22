/*---------------------------------------------------------------------------
 *
 *      Module:      WriteBinaryRestart.c
 *      Description: Contains functions needed to create the restart
 *                   control file and binary (HDF) nodal data file.
 *
 *      Includes functions:
 *
 *          CreateHDFGroup()
 *          CreateHDFDataset()
 *          WriteBinDataParams()
 *          WriteBinData()
 *          WriteBinaryRestart()
 *
 *      Structure of the HDF file is:
 *
 *       /dataFileVersion    int
 *       /numFileSegments    int
 *       /nodeCount          int
 *       /minCoordinates     double[3]
 *       /maxCoordinates     double[3]
 *       /dataDecompGeometry int[3]
 *       /dataDecompType     int
 *       /decomposition      double[*]
 *       /taskIDRange        int[2]
 *       /task0/
 *           nodeCount      int
 *           segCount       int
 *           nodeIndex      int[nodeCount]
 *           nodeConstraint int[nodeCount]
 *           nodeNumSegs    int[nodeCount]
 *           nodePos        double[nodeCount*3]
 *           segTags        int[segCount*3]
 *           burgersVec     double[segCount*3]
 *           glidePlane     double[segCount*3]
 *       /task1/
 *          .
 *          .
 *          .
 *
 *-------------------------------------------------------------------------*/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Restart.h"
#include "Decomp.h"
#include <sys/types.h>
#include <sys/stat.h>

#ifdef USE_HDF
#include <hdf5.h>
#endif


#ifdef USE_HDF
/*---------------------------------------------------------------------------
 *
 *      Function:    CreateHDFGroup
 *      Description: Create an HDF5 'group' of the specified name in the
 *                   file previously opened and associated with the 
 *                   provided file handle.
 *
 *      Arguments:
 *          fileID    File handle for an open HDF5 file.
 *          groupName Character string specifying the 'group' name
 *                    to be created in the file.
 *
 *      Returns:   0 on success
 *                 1 on error
 *
 *-------------------------------------------------------------------------*/
int CreateHDFGroup(hid_t fileID, char *groupName)
{
        hid_t  groupID;

        groupID = H5Gcreate(fileID, groupName, H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
        H5Gclose(groupID);

        return(groupID < 0);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    CreateHDFDataset
 *      Description: Create an HDF5 data item of the specified name in the
 *                   file previously opened and asscoiated with the 
 *                   provided file handle.
 *
 *                   Note: A request to create a dataset of zero elements
 *                   is silently ignored but treated as a successful request.
 *                   Attempting to actually create a zero element dataset
 *                   in the HDF5 file would cause an abort in the HDF library.
 *
 *      Arguments:
 *          fileID    File handle for an open HDF5 file.
 *          datasetName Character string specifying the name of the data
 *                    item to be created in the file.
 *          itemType  indicates the HDF5 data type
 *          numItems  Specifies the number of values of type <itemType>
 *                    to be written into the dataset.
 *          itemList  Pointer to data to be written.  Must point to
 *                    at least <numItems> values of the specified type
 *                    or the results are undefined.
 *
 *      Returns:   0 on success
 *                 1 on error
 *
 *-------------------------------------------------------------------------*/
int CreateHDFDataset(hid_t fileID, char *datasetName, hid_t itemType,
                     hsize_t numItems, void *itemList)
{
        int     numDims, hdfStatus;
        hid_t   dataspaceID, datasetID;

/*
 *      Don't create a zero-element dataset.  Could be problematic later.
 */
        if (numItems < 1) {
            return(0);
        }

/*
 *      We're only using 1 dimensional arrays for now.
 */
        numDims = 1;

        dataspaceID = H5Screate_simple(numDims, &numItems, NULL);
        if (dataspaceID < 0) {
            return(1);
        }

        datasetID = H5Dcreate(fileID, datasetName, itemType,
                              dataspaceID, H5P_DEFAULT, H5P_DEFAULT,
                              H5P_DEFAULT);
        if (datasetID < 0) {
            return(1);
        }

        if (numItems > 0) {
            hdfStatus = H5Dwrite(datasetID, itemType, H5S_ALL, H5S_ALL,
                                 H5P_DEFAULT, itemList);
        }

        H5Dclose(datasetID);
        H5Sclose(dataspaceID);

        return(hdfStatus != 0);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    WriteBinDataParams
 *      Description: Loop through the known data file parameters
 *                   and write each one out into its own HDF5 dataset.
 *
 *                   NOTE: We only write out the numeric parameters
 *                   in this function.  No other ones are needed
 *                   at this time.
 *
 *      Arguments:
 *          fileID    File handle for an open HDF5 file.
 *
 *      Returns:   0 on success
 *                 1 on error
 *
 *-------------------------------------------------------------------------*/
int WriteBinDataParams(Home_t *home, hid_t fileID)
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
            status |= CreateHDFDataset(fileID, itemName, hdfType,
                                       count, list->varList[i].valList);
        }

        return(status);
}
#endif  /* USE_HDF */


/*---------------------------------------------------------------------------
 *
 *      Function:    WriteBinData
 *      Description: Writes nodal data to the specified file.
 *
 *      Arguments:
 *          dataFile        Name of the data file.
 *          mode            1 to open file in write mode,
 *                          0 to open file in append mode.
 *          writePrologue   1 if this process should write all needed
 *                          headers and do any initialization associated
 *                          with the file, zero otherwise.
 *          writeEpilogue   1 if this process should write all needed
 *                          trailers and do any terminal processing
 *                          associated with the file, zero otherwise.
 *          binData         Pointer to a structure containing arrays
 *                          of node and segment data preformatted and
 *                          ready to be written to the HDF file.
 *
 *      Returns:   0 on success
 *                 non-zero on error
 *
 *-------------------------------------------------------------------------*/
int WriteBinData(Home_t *home, char *dataFile, int mode,
                 int writePrologue, int writeEpilogue,
                 BinFileData_t *binData)
{
        int      status = -1;
#ifdef USE_HDF
        int      numValues;
        int      taskIDRange[2];
        real8    *decompBounds;
        char     itemName[256], taskDir[128];
        hid_t    fileID;

        numValues = 0;
        decompBounds = (real8 *)NULL;

        if (mode) {
/*
 *          Create or overwrite the data file
 */
            fileID = H5Fcreate(dataFile, H5F_ACC_TRUNC, H5P_DEFAULT,
                               H5P_DEFAULT);
            if (fileID < 0) {
                printf("WriteBinData: H5Fcreate() error on %s", dataFile);
                return(-1);
            }

/*
 *          Store the ID range of tasks that will be writing to this file.
 */
            taskIDRange[0] = binData->firstInGroup;
            taskIDRange[1] = binData->lastInGroup;

            sprintf(itemName, "/taskIDRange");
            status = CreateHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                      2, taskIDRange);
            if (status != 0) {
                printf("WriteBinData: Error writing task ID range!");
                return(status);
            }

            if (writePrologue) {

/*
 *              Write the data file parameters
 */
                status = WriteBinDataParams(home, fileID);
                if (status != 0) {
                    printf("WriteBinData: Error writing basic parameters!");
                    return(status);
                }

/*
 *              Write the domain decomposition
 */
                GetAllDecompBounds(home, &decompBounds, &numValues);

                sprintf(itemName, "/decomposition");
                status = CreateHDFDataset(fileID, itemName, H5T_NATIVE_DOUBLE,
                                          numValues, decompBounds);
                if (status != 0) {
                    printf("WriteBinData: Error writing decomposition!");
                    return(status);
                }

                free(decompBounds);
                decompBounds = (real8 *)NULL;
            }
        } else {
/*
 *          Open an existing data file in a write mode
 */
            fileID = H5Fopen(dataFile, H5F_ACC_RDWR, H5P_DEFAULT);
            if (fileID < 0) {
                printf("WriteBinData: H5Fopen() error on %s", dataFile);
                return(-1);
            }
        }

/*
 *      Now create the task specific sub-directory structure in the
 *      file and dump the information into the data file.
 *
 *      NOTE: If the node count for this task is zero, the
 *      calls to create the arrays of nodal and segment data will
 *      silently return success so as not attempt to create zero-length
 *      datasets (which are a bad thing) in the hdf file.  The
 *      reader of the file should check the node count first to
 *      determine the necessity of reading additional data.
 */
        status = 0;

        sprintf(taskDir, "/task%d", home->myDomain);
        status |= CreateHDFGroup(fileID, taskDir);

        sprintf(itemName, "%s/nodeCount", taskDir);
        status |= CreateHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                   1, &binData->nodeCount);

        sprintf(itemName, "%s/nodeIndex", taskDir);
        status |= CreateHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                   binData->nodeCount, binData->nodeIndex);

        sprintf(itemName, "%s/nodePos", taskDir);
        status |= CreateHDFDataset(fileID, itemName, H5T_NATIVE_DOUBLE,
                                   binData->nodeCount * 3, binData->nodePos);

        sprintf(itemName, "%s/nodeConstraint", taskDir);
        status |= CreateHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                   binData->nodeCount, binData->nodeConstraint);

        sprintf(itemName, "%s/nodeNumSegs", taskDir);
        status |= CreateHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                   binData->nodeCount, binData->nodeNumSegs);

        sprintf(itemName, "%s/segCount", taskDir);
        status |= CreateHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                   1, &binData->segCount);

        sprintf(itemName, "%s/segTags", taskDir);
        status |= CreateHDFDataset(fileID, itemName, H5T_NATIVE_INT,
                                   binData->segCount * 3, binData->segTags);

        sprintf(itemName, "%s/burgersVec", taskDir);
        status |= CreateHDFDataset(fileID, itemName, H5T_NATIVE_DOUBLE,
                                   binData->segCount * 3, binData->burgersVec);

        sprintf(itemName, "%s/glidePlane", taskDir);
        status |= CreateHDFDataset(fileID, itemName, H5T_NATIVE_DOUBLE,
                                   binData->segCount * 3, binData->glidePlane);

        if (status != 0) {
            printf("WriteBinData: Error writing %s data\n", taskDir);
            return(status);
        }

        status = H5Fclose(fileID);
        if (status != 0) {
            printf("WriteBinData: H5Fclose() error on %s\n", dataFile);
        }
#endif  /* USE_HDF */

        return(status);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    WriteBinaryRestart
 *      Description: Writes control and nodal data for the current domain
 *                   to files with names derived from a base name.
 *
 *      Arguments:
 *          baseFileName     Base name of the restart file.  Restart data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          ioGroup          I/O group number associated with this domain
 *          firstInGroup     1 if this domain is the first processor in
 *                           its I/O group, zero otherwise.
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the file, zero otherwise
 *          binData          Pointer to a structure containing arrays
 *                           of node and segment data preformatted and
 *                           ready to be written to the HDF file.
 *
 *-------------------------------------------------------------------------*/
void WriteBinaryRestart(Home_t *home, char *baseFileName, int ioGroup,
                        int firstInGroup, int writePrologue, int writeEpilogue,
                        BinFileData_t *binData)
{
        int     status;
        char    *suffix, *start;
        char    fileName[256], ctrlFile[256], dataFile[256];

/*
 *      NOTE: The name of the nodal data file for this restart will
 *      be the same as the control file name with the exception
 *      that a ".hdf[.seqnum]" suffix will replace the file
 *      name suffix i.e. '.<anything>') of the control file name,
 *      or if the control file has no suffix, the new suffix will
 *      simply be appended. 
 */
        snprintf(ctrlFile, sizeof(ctrlFile), "%s/%s", DIR_RESTART,
                 baseFileName);

        strcpy(fileName, baseFileName);

        start = strrchr(fileName, '/');
        suffix = strrchr(fileName, '.');

        if (start == (char *)NULL) {
            start = fileName;
        }

        if ((suffix != (char *)NULL) && (suffix > start)) {
            *suffix = 0;
        }

/*
 *      Set control and data file names.  Only append a sequence
 *      number to the data file name if the data is to be spread
 *      across multiple files.
 */

        if (home->param->numIOGroups == 1) {
            snprintf(dataFile, sizeof(dataFile), "%s/%s%s",
                     DIR_RESTART, fileName, HDF_DATA_FILE_SUFFIX);
        } else {
            snprintf(dataFile, sizeof(dataFile), "%s/%s%s.%d",
                     DIR_RESTART, fileName, HDF_DATA_FILE_SUFFIX, ioGroup);
        }

#ifdef PARALLEL
#ifdef DO_IO_TO_NFS
/*
 *      It appears that when multiple processes on different hosts
 *      write to the same NFS-mounted file (even if access is
 *      sequentialized), consistency problems can arise resulting
 *      in corrupted data files.  Explicitly doing a 'stat' of the
 *      file on each process immediately before opening it is a kludge,
 *      but *seems* to force the NFS client on the local machine to 
 *      invalidate any cached data and acquire accurate info for
 *      the file and avoids the problem.
 */
        struct  stat statbuf;
        memset(&statbuf, 0, sizeof(statbuf));
        (void) stat(dataFile, &statbuf);
#endif
#endif

/*
 *      If this process is the first member of the first I/O group
 *      it needs to create the control file
 */
        if (firstInGroup && writePrologue) {
            status = WriteCtrl(home, ctrlFile);
            if (status != 0) {
                Fatal("WriteBinaryRestart: WriteCtrl() error on %s", ctrlFile);
            }
        }

/*
 *      First task in the I/O group must open the data file for writing
 *      to overwrite any existing file of the same name, otherwise append.
 */
        status = WriteBinData(home, dataFile, firstInGroup,
                              writePrologue, writeEpilogue,
                              binData);
        if (status != 0) {
            Fatal("WriteBinaryRestart: WriteBinData() error on %s", dataFile);
        }

        return;
}
