#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"

/*---------------------------------------------------------------------------
 *
 *      Function:    Tecplot
 *      Description: Write segment data in tecplot format... Masato
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
 *          numSegs          count of all segments in the problem space
 *
 *-------------------------------------------------------------------------*/
void Tecplot(Home_t *home, char *baseFileName, int ioGroup, int firstInGroup,
             int writePrologue, int writeEpilogue, int numSegs)
{
        int      i, j, thisDomain;
        int      newNodeKeyPtr, narm;
        int      btype=0;  
        real8    x,y,z;
        real8    x2,y2,z2;
        real8    x2o,y2o,z2o;
        real8    prx,pry,prz;
        real8    vx,vy,vz;
        real8    bX, bY, bZ;
        real8    Lx, Ly, Lz;
        char     fileName[256];
        Node_t   *node, *nbrNode;
        Param_t  *param;
        FILE     *fp;


        param      = home->param;
        thisDomain = home->myDomain;

        Lx = param->Lx;
        Ly = param->Ly;
        Lz = param->Lz;

/*
 *      Set data file name.  Only append a sequence number to
 *      the data file name if the data is to be spread across
 *      multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s",
                     DIR_TECPLOT, baseFileName);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%d",
                     DIR_TECPLOT, baseFileName, ioGroup);
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
 *      it needs to create the output file and do any intializations
 *      needed.
 */
        if (firstInGroup) {
/*
 *          First task in the I/O group must open the data file for writing
 *          to overwrite any existing file of the same name.
 */
            if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("tec_plot: Open error %d on %s\n", errno, fileName);
            }

            if (writePrologue) {
                printf(" +++ Writing Tecplot file(s) %s\n", baseFileName);

                fprintf(fp, "variables = X,Y,Z,V1,V2,V3,V4,V5,V6,V7,V8\n");
                fprintf(fp, "zone i = %d  F=POINT\n", 2*numSegs); 
            }
        } else {
/*
 *          Any task NOT first in its I/O group must open the data file
 *          in an append mode so everything gets added to the end of
 *          the file.
 */
            if ((fp = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("tec_plot: Open error %d on %s\n", errno, fileName);
            }
        }
 

/*
 *      Generate the plot data for the segments associated with this
 *      domain's data
 */
        newNodeKeyPtr = home->newNodeKeyPtr;

        for (i = 0; i < newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {;
                continue;
            }
        
            x = node->x;
            y = node->y;
            z = node->z;
          
/*
 *          Don't assume PBC is enabled...
 */
            if (param->xBoundType == Periodic) {
                x = x - Lx*rint(x/Lx);
            }

            if (param->yBoundType == Periodic) {
                y = y - Ly*rint(y/Ly);
            }

            if (param->zBoundType == Periodic) {
                z = z - Lz*rint(z/Lz);
            }
        
            for (j = 0; j < node->numNbrs; j++) {
#if 0
/*
 *              This check should be obsolete; should never
 *              have index < zero at this point in execution
 */
                if ((index = node->nbrTag[j].index) < 0) {
                    continue; 
                }
#endif
        
/*
 *  FIX ME? This will result in segments crossing domain boundaries
 *  to be added to the file by both domains!  Is this what is wanted?
 */
                if ((node->nbrTag[j].domainID == thisDomain) && 
                    (node->nbrTag[j].index < i)) {
                    continue;
                }

                bX = node->burgX[j];
                bY = node->burgY[j];
                bZ = node->burgZ[j];

/*
 *              For the following burgers vector checks, convert the burgers
 *              vector to the crystalographic frame if necessary.
 */
                if (param->useLabFrame) {
                    real8 burgLab[3] = {bX, bY, bZ}, burgCrystal[3];

                    Matrix33Vector3Multiply(param->rotMatrixInverse, burgLab,
                                            burgCrystal);
                    bX = burgCrystal[X];
                    bY = burgCrystal[Y];
                    bZ = burgCrystal[Z];
                }

                if (bX*bY*bZ == 0.0) {
                    if (bX*bY != 0.0 ||
                        bY*bZ != 0.0 ||
                        bZ*bX != 0.0) {
/*
 *                      unconventional [110] junction (only 1 comp is zero)
 */
                        btype=10;
                    } else {
/*
 *                      conventional [100] junction (5)
 */
                        btype=0;
                    }
                } else {
/*
 *                  nonzero BVs
 */
                    if(bX*bY*bZ<0.0) {
                        bX*=-1; bY*=-1; bZ*=-1;
                    }

/*
 *                  bv magnitude is not equal
 */
                    if (fabs(fabs(bX)-fabs(bY)) > 1e-2 ||
                        fabs(fabs(bY)-fabs(bZ)) > 1e-2 ||
                        fabs(fabs(bZ)-fabs(bX)) > 1e-2) {
                        btype = 20;
                    } else {
/*
 *                      products are always positive
 *                      this is rudandant condition 
 */
                        if (bY < 0 && bZ < 0) btype = 1;
                        if (bZ < 0 && bX < 0) btype = 2;
                        if (bX < 0 && bY < 0) btype = 3;
                        if (bX > 0 && bY > 0 && bZ > 0) btype = 4;
                    }
                }
        
                nbrNode = GetNeighborNode(home, node, j);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                x2o = nbrNode->x;
                y2o = nbrNode->y;
                z2o = nbrNode->z;
        
                prx = x2o-x;
                pry = y2o-y;
                prz = z2o-z;
        
                if (param->xBoundType == Periodic) {
                    prx = prx - Lx*rint(prx/Lx);
                }

                if (param->yBoundType == Periodic) {
                    pry = pry - Ly*rint(pry/Ly);
                }

                if (param->zBoundType == Periodic) {
                    prz = prz - Lz*rint(prz/Lz);
                }
        
                x2 = x + prx;
                y2 = y + pry;
                z2 = z + prz;

                vx =  node->vX * param->burgMag;
                vy =  node->vY * param->burgMag;
                vz =  node->vZ * param->burgMag;
                  
                narm = node->numNbrs;
        
/*
 *              This part is for output of junction investigation 
 */
                fprintf(fp, "%7.1f %7.1f %7.1f %6.1f %6.1f %6.1f %7.4f %7.4f %7.4f %d %d\n", 
                        x,y,z,prx,pry,prz,vx,vy,vz,narm,btype);
                fprintf(fp, "%7.1f %7.1f %7.1f %6.1f %6.1f %6.1f %7.4f %7.4f %7.4f %d %d\n", 
                        x2,y2,z2,-prx,-pry,-prz,0.0,0.0,0.0,narm,btype);
         
            }
        }
        
/*
 *      Handle anything that needs to be done after all data has been
 *      written to the file
 */
        if (writeEpilogue) {
            fprintf(fp, "\n");
        }
      
        fclose(fp);
       
	return;
}
