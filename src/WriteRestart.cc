/*---------------------------------------------------------------------------
 *
 *      Module:      WriteRestart.c
 *      Description: Contains functions needed to create the control and
 *                   nodal data files for restarts
 *
 *      Includes functions:
 *
 *          SetLatestRestart()
 *          WriteCtrlParams()
 *          WriteRestart()
 *
 *-------------------------------------------------------------------------*/
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Restart.h"
#include "Decomp.h"

/*-------------------------------------------------------------------------
 *
 *      Function:    SetLatestRestart
 *      Description: Writes the base name of the most recently written
 *                   restart file into the "latest_restart" file.
 *
 *------------------------------------------------------------------------*/
void SetLatestRestart(char *fileName)
{
        FILE *fp;

        fp = fopen("latest_restart", "w");
        fprintf(fp, "%s/%s\n", DIR_RESTART, fileName);
        fclose(fp);

        printf(" +++ Wrote restart file(s) %s\n", fileName);

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:    WriteCtrl
 *      Description: Write control parameter data to specified file.
 *
 *      Arguments:
 *          ctrlFile    Control file name.
 *
 *      Returns:   0 on success
 *                 non-zero on error
 *
 *-------------------------------------------------------------------------*/
int WriteCtrl(Home_t *home, char *ctrlFile)
{
        FILE    *fpCtrl;

        if ((fpCtrl = fopen(ctrlFile, "w")) == (FILE *)NULL) {
            printf("WriteCtrl: Open error %d on %s", errno, ctrlFile);
            return(errno);
        }

/*
 *      Restart at current cycle count
 */
        home->param->cycleStart = home->cycle;

/*
 *      Write out all the control parameters
 */
        fprintf(fpCtrl, "########################################\n");
        fprintf(fpCtrl, "###                                  ###\n");
        fprintf(fpCtrl, "###  ParaDiS control parameter file  ###\n");
        fprintf(fpCtrl, "###                                  ###\n");
        fprintf(fpCtrl, "########################################\n\n");
        WriteParam(home->ctrlParamList, -1, fpCtrl);

        return(fclose(fpCtrl));
}

/*---------------------------------------------------------------------------
 *
 *      Function:    WriteData
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
 *
 *      Returns:   0 on success
 *                 non-zero on error
 *
 *-------------------------------------------------------------------------*/
int WriteData(Home_t *home, char *dataFile, int mode,
              int writePrologue, int writeEpilogue)
{
        int     i, newNodeKeyPtr;
        int     iArm;
        FILE    *fp;
        Node_t  *node;

        if ((fp = fopen(dataFile, mode ? "w" : "a")) == (FILE *)NULL) {
            printf("WriteData: Open error %d on %s\n", errno, dataFile);
            return(errno);
        }

        if (writePrologue) {
/*
 *          Write the data file parameters
 */
            WriteParam(home->dataParamList, -1, fp);

/*
 *          Write the domain decomposition into nodal data file
 *          and then some comment lines describing the nodal
 *          data that will follow.
 */
            fprintf(fp, "\n#\n#  END OF DATA FILE PARAMETERS\n#\n\n");
            fprintf(fp, "domainDecomposition = \n");
            WriteDecompBounds(home, fp);

            fprintf(fp, "nodalData = \n");
            fprintf(fp, "#  Primary lines: node_tag, x, y, z, "
                    "num_arms, constraint\n");
            fprintf(fp, "#  Secondary lines: arm_tag, burgx, burgy, "
                         "burgz, nx, ny, nz\n");
        }

/*
 *      Now dump the data for all nodes in this block.
 */
        newNodeKeyPtr = home->newNodeKeyPtr;

        for (i = 0; i < newNodeKeyPtr; i++) {
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

/*
 *          For now, add a temporary sanity check. If we find any 
 *          unconstrained node (i.e. not a pinned node, surface
 *          node, etc) that is singly-linked we've got a problem,
 *          so abort without writing the restart file.
 *
 *          Once we identify the problem in the code that is leading to
 *          the singly-linked nodes, we can get rid of this check.
 *
 *          NOTE: For now, a node that is pinned in ANY dimension is
 *                treated here as if it is pinned in ALL dimensions.
 */
            if ((node->numNbrs == 1) &&
                HAS_NONE_OF_CONSTRAINTS(node->constraint,
                                        PINNED_NODE | SURFACE_NODE))
            {
                PrintNode(node);
                Fatal("WriteData: Node (%d,%d) singly linked!",
                      node->myTag.domainID, node->myTag.index);
            }

            fprintf(fp,
                    " %d,%d %.8f %.8f %.8f %d %d\n",
                    node->myTag.domainID, node->myTag.index,
                    node->x, node->y, node->z, node->numNbrs,
                    node->constraint);
        
/*
 *          Write the segment specific data
 */
            for (iArm = 0; iArm < node->numNbrs; iArm++) {
                fprintf(fp, "   %d,%d %16.10e %16.10e %16.10e\n"
                        "       %16.10e %16.10e %16.10e\n",
                        node->nbrTag[iArm].domainID,
                        node->nbrTag[iArm].index, node->burgX[iArm],
                        node->burgY[iArm], node->burgZ[iArm],
                        node->nx[iArm], node->ny[iArm], node->nz[iArm]);
            }
        }

        return(fclose(fp));
}

/*---------------------------------------------------------------------------
 *
 *      Function:    WriteRestart
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
 *                           with the restart file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the file, zero otherwise
 *
 *-------------------------------------------------------------------------*/
void WriteRestart(Home_t *home, char *baseFileName, int ioGroup,
                  int firstInGroup, int writePrologue, int writeEpilogue)
{
        int     status;
        char    *suffix, *start;
        char    fileName[256], ctrlFile[256], dataFile[256];

/*
 *      NOTE: The name of the nodal data file for this restart will
 *      be the same as the control file name with the exception
 *      that a ".data[.seqnum]" suffix will replace the file
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
                     DIR_RESTART, fileName, NODEDATA_FILE_SUFFIX);
        } else {
            snprintf(dataFile, sizeof(dataFile), "%s/%s%s.%d",
                     DIR_RESTART, fileName, NODEDATA_FILE_SUFFIX, ioGroup);
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
                Fatal("WriteRestart: WriteCtrl() error on %s", ctrlFile);
            }
        }
/*
 *      First task in the I/O group must open the data file for writing
 *      to overwrite any existing file of the same name, otherwise append.
 */
        status = WriteData(home, dataFile, firstInGroup,
                           writePrologue, writeEpilogue);
        if (status != 0) {
            Fatal("WriteRestart: WriteData() error on %s", dataFile);
        }

        return;
}
