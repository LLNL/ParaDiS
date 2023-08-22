/*****************************************************************************
 *
 *      Module:         WritePovray.c
 *      Description:    This module contains general functions needed for
 *                      plotting dislocation segments and domain boundaries
 *                      formatted for use with povary.
 *
 ****************************************************************************/
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "Decomp.h"


/*---------------------------------------------------------------------------
 *
 *      Function:    WritePovray
 *      Description: Write nodal data formatted for use with
 *                   the Povray tool. The filew written by this
 *                   function must be post-processed via the
 *                   genPovrayFrames tool before being visualized
 *                   with Povray.  See documentation on this output
 *                   format for details.
 *
 *      Arguments:
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
 *          numSegs          count of all segments in the problem space
 *
 *-------------------------------------------------------------------------*/
void WritePovray(Home_t *home, char *baseFileName, int ioGroup,
                 int firstInGroup, int writePrologue, int writeEpilogue)
{
        int     i, iarm;
        int     color, transmit, radius, newNodeKeyPtr;
        real8   bx, by, bz;
        real8   dx, dy, dz;
        real8   nbrx, nbry, nbrz;
        char    fileName[256];
        Param_t *param;
        Node_t  *node, *nbrNode;
        FILE    *fp; 
           
        
        param=home->param;
        

/*
 *      Set data file name.  Only append a sequence number to
 *      the data file name if the data is to be spread across
 *      multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s",
                     DIR_POVRAY, baseFileName);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%d",
                     DIR_POVRAY, baseFileName, ioGroup);
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
                Fatal("WritePovray: Open error %d on %s\n", errno, fileName);
            }
            if (writePrologue) {
                printf(" +++ Writing Povray file(s) %s\n", baseFileName);
            }
        } else {
            if ((fp = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("WritePovray: Open error %d on %s\n", errno, fileName);
            }
        }
        

        newNodeKeyPtr = home->newNodeKeyPtr;
        
        for (i = 0; i < newNodeKeyPtr; i++) {
        
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }
        
            if (node->numNbrs == 2) {
                color=1;
            } else {
                color=2;
            }
              
            transmit = color;
            radius = color;

            fprintf(fp, "sphere{ <%e,%e,%e>, radius%02d "
                    "pigment { color color%02d "
                    "transmit transmit%02d } "
                    "finish { phong 1 metallic }"
                    " }\n",
                    node->x, node->y, node->z,
                    radius,color,transmit);

/*
 *          Check each of the node's arms. If nbrTag is lower
 *          than myTag, print the arm info
 */
            for (iarm = 0; iarm < node->numNbrs; iarm++) {

                nbrNode = GetNeighborNode(home, node, iarm);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                if (OrderNodes(node, nbrNode) > 0) {
        
                    bx = node->burgX[iarm];
                    by = node->burgY[iarm];
                    bz = node->burgZ[iarm];
        
                    if ((fabs(fabs(bx)-fabs(by))<1e-3)
                        &&(fabs(fabs(bx)-fabs(bz))<1e-3)) {
                        color=1;
                    } else {
                        color=2;
                    }

                    if (node->myTag.domainID !=
                        node->nbrTag[iarm].domainID) {
                        color=3;
                    }
                         
                    transmit = color;
                    radius = color;

                    nbrx = nbrNode->x;
                    nbry = nbrNode->y;
                    nbrz = nbrNode->z;
        
                    dx = nbrx - node->x;
                    dy = nbry - node->y;
                    dz = nbrz - node->z;

                    ZImage(param, &dx, &dy, &dz);
                         
                    fprintf(fp, "cylinder{ <%e,%e,%e>,<%e,%e,%e>,radius%02d "
                            "pigment { color color%02d "
                            "transmit transmit%02d } "
                            "finish { phong 1 metallic }"
                            " }\n",
                            node->x, node->y, node->z,
                            node->x+dx, node->y+dy, node->z+dz,
                            radius,color,transmit);
                }
            }
        }

/*
 *      Deal with all the local Eshelby inclusions (if any)
 */
        color = 4;
        transmit = 4;

#ifdef ESHELBY
        for (i = 0; i < home->locInclusionCount; i++) 
        {

            EInclusion_t *inclusion = &home->eshelbyInclusions[i];

//SA : want an ellipsoid here instead of a sphere.....

            radius = MAX(1, (int)inclusion->radius[0]);

            fprintf(fp, "sphere{ <%e,%e,%e>, %d pigment "
                    "{ color color%02d transmit transmit%02d } "
                    "finish { phong 1 metallic } }\n",
                    inclusion->position[X], inclusion->position[Y],
                    inclusion->position[Z], radius, color, transmit);
        }
#endif

/*
 *      If necesary (i.e. this is the final processor in the
 *      last I/O group), handle anything epilogue stuff that
 *      needs doing.
 */
        if (writeEpilogue) {
        }
        
        fclose(fp);

        return; 
}
