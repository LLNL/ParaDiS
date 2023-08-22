#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"

/*---------------------------------------------------------------------------
 *
 *      Function:     Gnuplot
 *      Description:  Create plot file in format used for moviemaking
 *                    with gnuplot (M. Rhee)
 *
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
void Gnuplot(Home_t *home, char *baseFileName, int ioGroup, int firstInGroup,
             int writePrologue, int writeEpilogue)
{
	int		i, j, index, thisDomain;
	int		newNodeKeyPtr;
	real8		x, y, z, x2, y2, z2;
	real8		prx, pry, prz;
	real8	        Lx, Ly, Lz, x0, y0, z0, x1, y1, z1;
        char            fileName[256], tmpName[256];
	Node_t		*node, *nbrNode;
	Param_t		*param; /* Need this to plot box */
	FILE	        *fp;


        thisDomain = home->myDomain;
	param      = home->param;

        x0 = param->minSideX;
        y0 = param->minSideY;
        z0 = param->minSideZ;

        x1 = param->maxSideX;
        y1 = param->maxSideY;
        z1 = param->maxSideZ;

        Lx = x1 - x0;
        Ly = y1 - y0;
        Lz = z1 - z0;

/*
 *      Set data file name.  Only append a sequence number to
 *      the data file name if the data is to be spread across
 *      multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s",
                     DIR_GNUPLOT, baseFileName);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%d",
                     DIR_GNUPLOT, baseFileName, ioGroup);
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
 *      If this process is the first member of its I/O group
 *      it needs to create the output file
 */
        if (firstInGroup) {

/*
 *	    If necessary, do the basic initialization for the plot file
 */
            if (writePrologue) {

                snprintf(tmpName, sizeof(tmpName), "%s/box.in", DIR_GNUPLOT);
                fp=fopen(tmpName, "w");

                fprintf(fp,"%f %f %f\n",x0,y0,z0);
                fprintf(fp,"%f %f %f\n",x0,y0,z1);
                fprintf(fp,"\n\n\n");

                fprintf(fp,"%f %f %f\n",x0,y1,z0);
                fprintf(fp,"%f %f %f\n",x0,y1,z1);
                fprintf(fp,"\n\n\n");

                fprintf(fp,"%f %f %f\n",x1,y1,z0);
                fprintf(fp,"%f %f %f\n",x1,y1,z1);
                fprintf(fp,"\n\n\n");

                fprintf(fp,"%f %f %f\n",x1,y0,z0);
                fprintf(fp,"%f %f %f\n",x1,y0,z1);
                fprintf(fp,"\n\n\n");

                fprintf(fp,"%f %f %f\n",x0,y0,z0);
                fprintf(fp,"%f %f %f\n",x0,y1,z0);
                fprintf(fp,"\n\n\n");

                fprintf(fp,"%f %f %f\n",x0,y0,z1);
                fprintf(fp,"%f %f %f\n",x0,y1,z1);
                fprintf(fp,"\n\n\n");

                fprintf(fp,"%f %f %f\n",x1,y0,z1);
                fprintf(fp,"%f %f %f\n",x1,y1,z1);
                fprintf(fp,"\n\n\n");

                fprintf(fp,"%f %f %f\n",x1,y0,z0);
                fprintf(fp,"%f %f %f\n",x1,y1,z0);
                fprintf(fp,"\n\n\n");

                fprintf(fp,"%f %f %f\n",x0,y0,z0);
                fprintf(fp,"%f %f %f\n",x1,y0,z0);
                fprintf(fp,"\n\n\n");

                fprintf(fp,"%f %f %f\n",x0,y0,z1);
                fprintf(fp,"%f %f %f\n",x1,y0,z1);
                fprintf(fp,"\n\n\n");

                fprintf(fp,"%f %f %f\n",x0,y1,z1);
                fprintf(fp,"%f %f %f\n",x1,y1,z1);
                fprintf(fp,"\n\n\n");

                fprintf(fp,"%f %f %f\n",x0,y1,z0);
                fprintf(fp,"%f %f %f\n",x1,y1,z0);
                fprintf(fp,"\n\n\n");

                fclose(fp);

                printf(" +++ Writing Gnuplot file(s) %s\n", baseFileName);

            }  /* if (writePrologue) */

/*
 *          First task in the I/O group must open the data file for writing
 *          to overwrite any existing file of the same name.
 */
            if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("gnu_plot: Open error %d on %s\n", errno, fileName);
            }

        } else {
/*
 *          Any task NOT first in its I/O group must open the data file
 *          in an append mode so everything gets added to the end of
 *          the file.
 */
            if ((fp = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("gnu_plot: Open error %d on %s\n", errno, fileName);
            }
        }

/*
 *      Plot segments
 */
        newNodeKeyPtr = home->newNodeKeyPtr;

        for (i = 0; i < newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            x = node->x;
            y = node->y;
            z = node->z;

            x = x - Lx*rint(x/Lx);
            y = y - Ly*rint(y/Ly);
            z = z - Lz*rint(z/Lz);
            
            for (j = 0; j < node->numNbrs; j++) {

                index = node->nbrTag[j].index;

                if (index < 0) {
                    continue;
                }

                if ((node->nbrTag[j].domainID == thisDomain) && (index < i)) {
                    continue;
                }

                if (node->nbrTag[j].domainID < thisDomain) {
                    continue;
                }

                nbrNode = GetNeighborNode(home, node, j);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                x2 = nbrNode->x;
                y2 = nbrNode->y;
                z2 = nbrNode->z;

                prx = x2 - x;
                pry = y2 - y;
                prz = z2 - z;

                prx = prx - Lx*rint(prx/Lx);
                pry = pry - Ly*rint(pry/Ly);
                prz = prz - Lz*rint(prz/Lz);

                x2 = x + prx;
                y2 = y + pry;
                z2 = z + prz;

                fprintf(fp, "%f %f %f\n", x, y, z);
                fprintf(fp, "%f %f %f\n", x2, y2, z2);
                fprintf(fp, "# %d.%d %d\n", node->myTag.domainID,
                        node->myTag.index, 0);
                fprintf(fp, "\n");
                fprintf(fp, "\n");                       
            }
        }

/*
 *      No epilogue stuff needed, so just close the output file
 *      and we're done.
 */
        fclose(fp);

        return;
}
