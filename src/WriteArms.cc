#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"

/*---------------------------------------------------------------------------
 *
 *      Function:     WriteArms
 *      Description:  For all active nodes in a domain, write the arm
 *                    burgers vectors and line directions out to a file.
 *      Args:
 *          baseFileName     Base name of the plot file.  Plot data
 *                           will be written to 1 or more file segments
 *                           named <baseFileName>.n
 *          ioGroup          I/O group number associated with this domain
 *          firstInGroup     1 if this domain is the first processor in
 *                           its I/O group, zero otherwise.
 *          writePrologue    1 if this process should write all needed
 *                           headers and do any initialization associated
 *                           with the plot file, zero otherwise
 *          writeEpilogue    1 if this process should write all needed
 *                           trailers and do any terminal processing
 *                           associated with the plot file, zero otherwise
 *
 *-------------------------------------------------------------------------*/
void WriteArms(Home_t *home, char *baseFileName, int ioGroup, int firstInGroup,
               int writePrologue, int writeEpilogue)
{
        int      i, iarm, newNodeKeyPtr;
        real8    ux, uy, uz, disMag;
        real8    dx, dy, dz;
        real8    nbrx, nbry, nbrz;
        char     fileName[256];
        Node_t   *node, *nbrNode;
        Param_t  *param;
        FILE     *fp;


        param=home->param;

/*
 *      Set data file name.  Only append a sequence number to
 *      the data file name if the data is to be spread across
 *      multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s",
                     DIR_ARMDATA, baseFileName);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%d",
                     DIR_ARMDATA, baseFileName, ioGroup);
        }

#ifdef PARALLEL
#ifdef DO_IO_TO_NFS
/*
 *      It appears that when multiple processes on different hosts
 *      write to the same NFS-mounted file (even if access is
 *      sequentialized), consistency problems can arise resulting
 *      in corrupted data files.  Explicitly doing a 'stat' of the
 *      file on each process immediately before opening it *seems*
 *      to clear up the problem.
 */
        struct stat statbuf;
        memset(&statbuf, 0, sizeof(statbuf));
        (void) stat(fileName, &statbuf);
#endif
#endif

/*
 *      First task in the I/O group must open the data file for writing
 *      to overwrite any existing file of the same name, all other
 *      tasks in I/O group must open the data file in an append mode
 *      so everything gets added to the end of the file.
 */
        if (firstInGroup) {
            if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("WriteArms: Open error %d on %s\n", errno, fileName);
            }
            if (writePrologue) {
                printf(" +++ Writing arm data file(s) %s\n", baseFileName);
            }
        } else {
            if ((fp = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("WriteArms: Open error %d on %s\n", errno, fileName);
            }
        }


        newNodeKeyPtr = home->newNodeKeyPtr;

        for (i=0;i<newNodeKeyPtr;i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

/*
 *          Check each of the node's arms. If nbrTag is lower than myTag,
 *          print the arm info
 */
            for (iarm = 0; iarm < node->numNbrs; iarm++) {

                nbrNode = GetNeighborNode(home, node, iarm);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                if (OrderNodes(node, nbrNode) > 0) {

/*
 *                  Get the neighbor's coordinates and the line
 *                  direction vector
 */
                    nbrx = nbrNode->x;
                    nbry = nbrNode->y;
                    nbrz = nbrNode->z;

                    dx = nbrx - node->x;
                    dy = nbry - node->y;
                    dz = nbrz - node->z;

                    ZImage(param, &dx, &dy, &dz);

                    disMag = sqrt(dx*dx + dy*dy + dz*dz);

                    ux = dx / disMag;
                    uy = dy / disMag;
                    uz = dz / disMag;

/*
 *                  Print the burgers vector and line direction vector
 */
                    fprintf(fp,"%e %e %e %e %e %e %e %e %e %e\n",
                            node->burgX[iarm], node->burgY[iarm],
                            node->burgZ[iarm], ux, uy, uz, disMag,
                            node->x, node->y, node->z);

                }
            }
        }

        fclose(fp);

        return;
}
