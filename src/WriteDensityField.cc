#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Tag.h"
#include "Util.h"

#define MAXDENSITYMESH 51

/*---------------------------------------------------------------------------
 *
 *      Function:    WriteDensityField
 *      Description: Write a file containing dislocation density field in
 *                   3-d array... Wei Cai
 *
 *      Args:
 *
 *-------------------------------------------------------------------------*/
void WriteDensityField(Home_t *home, char *fileName)
{
        int     nx, ny, nz, elemCount;
        int     i, j, k, im, jm, km, i0, j0, k0, di, dj, dk;
        int     newNodeKeyPtr;
        int     thisDomain;
        real8   x, y, z, x1, y1, z1, dx, dy, dz, xm, ym, zm;
        real8   dr2, rho;
        real8   (*totDensity)[MAXDENSITYMESH][MAXDENSITYMESH][MAXDENSITYMESH];
        Node_t  *node, *nbrNode;
        Param_t *param;
        static real8 xmin, xmax, ymin, ymax, zmin, zmax, Lx, Ly, Lz;
        static real8 den[MAXDENSITYMESH][MAXDENSITYMESH][MAXDENSITYMESH];
            
        
        thisDomain = home->myDomain;
        param      = home->param;
        
        nx = param->savedensityspec[0];
        ny = param->savedensityspec[1];
        nz = param->savedensityspec[2];
        
        if ((nx <= 1) || (ny <= 1) || (nz <= 1)) {
            return;
        }

        if ((nx > MAXDENSITYMESH) || (ny > MAXDENSITYMESH) ||
            (nz > MAXDENSITYMESH)) {
            if (thisDomain == 0) {
                printf("WriteDensityField: (%d,%d,%d) > limit (%d)\n",
                       nx, ny, nz, MAXDENSITYMESH);
            }
            return;
        }
        
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                for (k = 0; k < nz; k++) {
                    den[i][j][k] = 0;
                }
            }
        }

        xmin = param->minSideX;
        xmax = param->maxSideX;
        ymin = param->minSideY;
        ymax = param->maxSideY;
        zmin = param->minSideZ;
        zmax = param->maxSideZ;
        
        Lx = xmax - xmin;
        Ly = ymax - ymin;
        Lz = zmax - zmin;
         
        FILE *fp=0;

        if (thisDomain == 0) {

            if ((fp = fopen(fileName, "w")) == (FILE *)NULL) {
                Fatal("WriteDensityField: Open error %d on %s\n", errno, fileName);
            }
/*
 *          Write mesh dimensions parameters
 */
            fprintf(fp, "ParaDiS dislocation density file (viewable by vaspview)\n");
            fprintf(fp, "1.0\n");
            fprintf(fp, "%e 0 0\n", 1.0);
            fprintf(fp, "0 %e 0\n", Ly / Lx);
            fprintf(fp, "0 0 %e\n", Lz / Lx);
            fprintf(fp, "0\nDirect\n\n");
            fprintf(fp, "%d %d %d\n", nx, ny, nz);
        }

/*
 *      Go through all the nodes of this domain and update the density
 *      field array appropriately.
 */
        newNodeKeyPtr = home->newNodeKeyPtr;
        
        for (i=0;i<newNodeKeyPtr;i++) {
        
            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            x = node->x;
            y = node->y;
            z = node->z;
                
            for (j = 0; j < node->numNbrs; j++) {

                if (node->nbrTag[j].index < 0) {
                    continue;
                }

                nbrNode = GetNeighborNode(home, node, j);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }
        
                x1 = nbrNode->x;
                y1 = nbrNode->y;
                z1 = nbrNode->z;

                dx = x1 - x;
                dy = y1 - y;
                dz = z1 - z;

                ZImage(param, &dx, &dy, &dz);

                rho = sqrt(dx*dx + dy*dy + dz*dz);
                    
                xm = x + dx*0.5 - (xmin+xmax)*0.5;
                ym = y + dy*0.5 - (ymin+ymax)*0.5;
                zm = z + dz*0.5 - (zmin+zmax)*0.5;

                xm /= Lx;
                ym /= Ly;
                zm /= Lz;

                xm -= rint(xm);
                ym -= rint(ym);
                zm -= rint(zm);
                    
                i0 = (int)rint((xm+0.5)*(nx-1));
                j0 = (int)rint((ym+0.5)*(ny-1));
                k0 = (int)rint((zm+0.5)*(nz-1));
        
                for (di = -5; di <= 5; di++) {
                    for (dj = -5; dj <= 5; dj++) {
                        for (dk = -5; dk <= 5; dk++) {

                            im = (i0+di+nx-1) % (nx-1);
                            jm = (j0+dj+ny-1) % (ny-1);
                            km = (k0+dk+nz-1) % (nz-1);

                            dr2= (((xm+0.5)*(nx-1)-i0-di) *
                                  ((xm+0.5)*(nx-1)-i0-di) +
                                  ((ym+0.5)*(ny-1)-j0-dj) *
                                  ((ym+0.5)*(ny-1)-j0-dj) +
                                  ((zm+0.5)*(nz-1)-k0-dk) *
                                  ((zm+0.5)*(nz-1)-k0-dk));

                            den[im][jm][km] += rho * exp(-dr2);
                        }
                    }
                }
            }
        }
        
/*
 *      Handle periodicity
 */
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                den[i][j][nz-1] = den[i][j][0];
            }
        }
        
        for (i = 0; i < nx; i++) {
            for (k = 0; k < nz; k++) {
                den[i][ny-1][k] = den[i][0][k];
            }
        }
        
        for (j = 0; j < ny; j++) {
            for (k = 0; k < nz; k++) {
                den[nx-1][j][k] = den[0][j][k];
            }
        }
        
/*
 *      Sum the density field from all domains leaving the result
 *      on processor zero.
 */
        elemCount  = MAXDENSITYMESH*MAXDENSITYMESH*MAXDENSITYMESH;
        totDensity = (real8 (*)[MAXDENSITYMESH][MAXDENSITYMESH][MAXDENSITYMESH]) calloc(1, sizeof(real8) * elemCount);

#ifdef PARALLEL
        MPI_Reduce( &den[0][0][0], (real8 *) totDensity, elemCount, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                for (k = 0; k < nz; k++) {
                    totDensity[0][i][j][k] = den[i][j][k];
                }
            }
        }
#endif

        if (thisDomain == 0) {
            for (i = 0; i < nx; i++) {
                for (j = 0; j < ny; j++) {
                    for (k = 0; k < nz; k++) {
                        if (totDensity[0][i][j][k] != 0) {
                            fprintf(fp, "%20.16e\n", totDensity[0][i][j][k]);
                        } else {
                            fprintf(fp, "0\n");
                        }
                    }
                }
            }

            fclose(fp);
        }

        free(totDensity);


        return; 
}
