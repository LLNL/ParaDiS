#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"

/*---------------------------------------------------------------------------
 *
 *      Function:    WritePoleFig
 *      Description: Creates <111> type burgers vector pole figures
 *                   (Moon Rhee)
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
void WritePoleFig(Home_t *home, char *baseFileName, int ioGroup,
                  int firstInGroup, int writePrologue, int writeEpilogue)
{
        int     i, ipole, iburg, newNodeKeyPtr;
        real8   ux, uy, uz, disMag;
        real8   dx, dy, dz;
        real8   nbrx, nbry, nbrz;
        real8   rmag,r00,r01,r02,r10,r11,r12,r20,r21,r22;
        real8   polx,poly,polz,den;
        real8   bx,by,bz,weight;
        char    fileName[256];
        FILE    *fp;
        Node_t  *node, *nbrNode;
        Param_t *param;
           
        
        param=home->param;
        
/*
 *      Set data file name.  Only append a sequence number to
 *      the data file name if the data is to be spread across
 *      multiple files.
 */
        if (param->numIOGroups == 1) {
            snprintf(fileName, sizeof(fileName), "%s/%s",
                     DIR_POLEFIG, baseFileName);
        } else {
            snprintf(fileName, sizeof(fileName), "%s/%s.%d",
                     DIR_POLEFIG, baseFileName, ioGroup);
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
                Fatal("WritePoleFig: Open error %d on %s\n", errno, fileName);
            }
            if (writePrologue) {
                printf(" +++ Writing PoleFig file(s) %s\n", baseFileName);
            }
        } else {
            if ((fp = fopen(fileName, "a")) == (FILE *)NULL) {
                Fatal("WritePoleFig: Open error %d on %s\n", errno, fileName);
            }
        }
        
        newNodeKeyPtr = home->newNodeKeyPtr;
        
        for(i=0;i<newNodeKeyPtr;i++) {
        
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }
        
            for (ipole = 0; ipole < node->numNbrs; ipole++) {
        
                bx=node->burgX[ipole];
                by=node->burgY[ipole];
                bz=node->burgZ[ipole];
        
                if ((bx > 0.0) && (by > 0.0) && (bz > 0.0)) {

                    rmag = 1.0 / sqrt(2.0);
                    r00 = -rmag;
                    r01 = rmag;
                    r02 = 0.0;

                    rmag = 1.0 / sqrt(3.0);
                    r20 = rmag;
                    r21 = rmag;
                    r22 = rmag;

                    iburg = 1;

                } else if ((bx < 0.0) && (by > 0.0) && (bz > 0.0)) {

                    rmag = 1.0 / sqrt(2.0);
                    r00 = rmag;
                    r01 = rmag;
                    r02 = 0.0;

                    rmag = 1.0 / sqrt(3.0);
                    r20 = -rmag;
                    r21 = rmag;
                    r22 = rmag;

                    iburg = 2; 

                } else if ((bx > 0.0) && (by < 0.0) && (bz > 0.0)) {

                    rmag = 1.0 / sqrt(2.0);
                    r00 = rmag;
                    r01 = rmag;
                    r02 = 0.0;

                    rmag = 1.0 / sqrt(3.0);
                    r20 = rmag;
                    r21 = -rmag;
                    r22 = rmag;

                    iburg = 3;

                } else if ((bx > 0.0) && (by > 0.0) && (bz < 0.0)) {

                    rmag = 1.0 / sqrt(2.0);
                    r00 = rmag;
                    r01 = -rmag;
                    r02 = 0.0;

                    rmag = 1.0 / sqrt(3.0);
                    r20 = rmag;
                    r21 = rmag;
                    r22 = -rmag;

                    iburg = 4;

                } else {
                    continue;
                }
            
                xvector(r00, r01, r02, r20, r21, r22, &r10, &r11, &r12);
        
                nbrNode = GetNeighborNode(home, node, ipole);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }
        
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
        
                polx = r00*ux + r01*uy + r02*uz;
                poly = r10*ux + r11*uy + r12*uz;
                polz = r20*ux + r21*uy + r22*uz;
         
                if (polz >= 0.) {
                    den = 1.0 + polz;
                } else {
                    den = 1.0 - polz;
                }
        
                polx /= den;
                poly /= den;
                weight = disMag / 1000.0;

                fprintf(fp, " %d %e %e %e\n", iburg, polx, poly, weight);
            }
        }
        
        fclose(fp);

        return;
}
        
