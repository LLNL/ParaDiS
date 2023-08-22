/***************************************************************************
 *
 *      Module:      FMAnisotropicInit.c
 *
 *      Description: This module contains functions to initialize certain
 *                   Fast Multipole tables needed for use with
 *                   anisotropic elasticity when the Taylor FMM method is on.
 *
 *      Includes public functions:
 *          FMAnisotropicInit()
 *
 *      Includes private functions:
 *          AnisotropicInitElasticConstants()
 *
 **************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

#if (defined ANISOTROPIC) && (defined TAYLORFMM)

void ReadDerivTableMetadata(FILE *fpTable, char *fileName, real8 C66[6][6],
                            real8 *xScale, real8 *yScale, real8 *zScale,
                            int *tblMaxOrder)
{
    int    i, j;
    int    maxLineLen;
    int    tag66[6][6];
    char   inLine[1024];

    maxLineLen = sizeof(inLine);

/*
 *  NOTE: Green's function derivatives file is assumed to be open already
 *
 *  Read the max order
 */
    fgets(inLine, maxLineLen, fpTable);

    if (sscanf(inLine, "%%%% maxorder = %d;",
               tblMaxOrder) != 1) {
        Fatal("FMAnisotropicInit: Error reading maxOrder from anistropy "
              "derivatives file %s", fileName);
    }

/*
 *  Read the cell size scale factors and store them for later
 */
    fgets(inLine, maxLineLen, fpTable);

    if (sscanf(inLine, "%%%% xscale = %lf; yscale = %lf; zscale = %lf;",
               xScale, yScale, zScale) != 3) {
        Fatal("FMAnisotropicInit: Error reading scale factors from "
              "anistropy derivatives file %s", fileName);
    }

/*
 *  Read the C matrix in 6x6 form and convert to the 3x3x3x3 form
 */
    for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
            tag66[i][j] = 0;
        }
    }

    for (i = 0; i < 6; i++) {
        int    k;
        int    ii[6], jj[6];
        real8  Cstrip[6];

        fgets(inLine, maxLineLen, fpTable);

        if (sscanf(inLine,
                   "%%%% C%1d%1d = %lf; C%1d%1d = %lf; C%1d%1d = %lf; "
                   "C%1d%1d = %lf; C%1d%1d = %lf; C%1d%1d = %lf;",
                   ii+0,jj+0,Cstrip+0,
                   ii+1,jj+1,Cstrip+1,
                   ii+2,jj+2,Cstrip+2,
                   ii+3,jj+3,Cstrip+3,
                   ii+4,jj+4,Cstrip+4,
                   ii+5,jj+5,Cstrip+5) != 3*6) {
            Fatal("FMAnisotropicInit: Error reading C matrix from "
                  "anistropy derivatives file %s", fileName);
        }

        for (k = 0; k < 6; k++) {

            if (ii[k] < 0 || ii[k] >= 6 || jj[k] < 0 || jj[k] >= 6) {
            
                Fatal("FMAnisotropicInit: Index error in C matrix from "
                      "anistropy derivatives file %s", fileName);
            }

            C66[ii[k]][jj[k]] = Cstrip[k];
            tag66[ii[k]][jj[k]]++;
        }
    }

    for (i = 0; i<6; i++) {
        for (j = 0; j<6; j++) {
            if (tag66[i][j] != 1) {
                Fatal("FMAnisotropicInit: Missing or redundant "
                      "elements  in C matrix from\nanistropy "
                      "derivatives file %s", fileName);
            }
        }
    }

    return;
}


void ReadDerivTable(Home_t *home, FILE *fpTable)
{
    int             i, j, k;
    int             tblMaxOrder, numDeriv;
    int             maxLineLen;
    char            inLine[1024];
    Param_t        *param;
    FMAnisoTable_t *anisoTable;


    param      = home->param;
    anisoTable = home->fmAnisoTable;

    maxLineLen = sizeof(inLine) / sizeof(char);
    tblMaxOrder = anisoTable->derivTblMaxOrder;
    numDeriv = (tblMaxOrder+3)*(tblMaxOrder+2)*(tblMaxOrder+1)/6;

/*
 *  Read the remainder of the derivatives files
 */
    while (fgets(inLine, maxLineLen, fpTable) != NULL) {
        int   ii, jj;
        int   n, idx;
        int   ix, iy, iz;
        int   ndx, ndy, ndz;
        real8 val;

        if (strlen(inLine) == 0) continue;

/*
 *      Read indices of the derivative table value and the derivative
 *      value from the file...
 *
 *      nd{xyz} => number of derivatives in the block of the
 *                 table at i{xyz}
 *      i{xyz}  => indices of the block of the derivatives table
 *                 containing the value.
 *      ii      => tensor index of the value in 3x3
 *      jj      =>
 */
        sscanf(inLine,"%d%d %d%d%d %d%d%d %lf",
               &ii, &jj,   &ndx, &ndy, &ndz,   &ix, &iy, &iz,   &val);

        if (ii<0 || ii>=3 || jj<0 || jj>=3) {
            Fatal("FMAnisotropicInit: ii or jj out of range in %s\n"
                  "ii = %d, jj = %d", param->fmAnisoDerivTbl, ii, jj);
        }

        if (ix*ix > 25 || iy*iy > 25 || iz*iz > 25) {
            Fatal("FMAnisotropicInit: ix, iy, or iz out of range in %s\n"
                  "ix = %d, iy = %d, iz = %d ",
                  param->fmAnisoDerivTbl, ii, jj);
        }

/*
 *      If we have not yet allocated the memory for this portion
 *      of the derivative table, go ahead and do that now and mark
 *      it as allocated.
 */
        if (anisoTable->tagmap[iz+5][iy+5][ix+5] == 0) {
            anisoTable->tagmap[iz+5][iy+5][ix+5] = 1;
            anisoTable->Gderiv_map[iz+5][iy+5][ix+5] =
                    (matrix *)malloc(sizeof(matrix) * numDeriv);
        }

/*
 *      Verify the number of derivative values matches the expected
 */
        n = ndx + ndy + ndz;
        idx = (n+2)*(n+1)*n/6 + ndz*(2*n+3-ndz)/2 + ndy;

        if ((idx < 0) || (idx >= numDeriv) || (n > tblMaxOrder) ||
	    (ndx < 0) || (ndy < 0) || (ndz < 0)) {
            Fatal("FMAnisotropicInit: number of derivatives is out of "
                  "range in %s\nexpected %d, found %d "
		  "(ndx=%d, ndy=%d, ndz=%d, n=%d, tblMaxOrder=%d)",
                  param->fmAnisoDerivTbl, numDeriv, idx,
		  ndx, ndy, ndz, n, tblMaxOrder);
        }

/*
 *      Store the derivative value at the appropriate spot in the table
 */
        anisoTable->Gderiv_map[iz+5][iy+5][ix+5][idx][ii][jj] = val;
    }

    return;
}


void DistributeDerivTable(Home_t *home)
{
    int            i, j, k;
    int            tblMaxOrder, numDeriv;
    FMAnisoTable_t *anisoTable;

    anisoTable = home->fmAnisoTable;

#ifdef PARALLEL
/*
 *  We've got the derivative table read in on task 0, now all the
 *  appropriate data has to be broadcast out to the remote domains...
 *
 *  First broadcast the taylor expansion order used when creating the
 *  derivatives file
 */
    MPI_Bcast(&anisoTable->derivTblMaxOrder, 1, MPI_INT, 0, MPI_COMM_WORLD);

/*
 *  Broadcast the cell scaling factors.
 */
    MPI_Bcast(anisoTable->scaleVec, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

/*
 *  Now the tagmap so that all processes will know which blocks of the
 *  derivative map have been populated
 */
    MPI_Bcast(&anisoTable->tagmap, 11*11*11, MPI_INT, 0, MPI_COMM_WORLD);

    anisoTable->xScale = anisoTable->scaleVec[0];
    anisoTable->yScale = anisoTable->scaleVec[1];
    anisoTable->zScale = anisoTable->scaleVec[2];

    tblMaxOrder = anisoTable->derivTblMaxOrder;
    numDeriv = (tblMaxOrder+3)*(tblMaxOrder+2)*(tblMaxOrder+1)/6;

/*
 *  Loop through all the blocks of derivatives.  Each block that
 *  has been populated gets broadcast from task 0 to all remote
 *  domains.
 */
    for (i = 0; i < 11; i++) {
        for (j = 0; j < 11; j++) {
            for (k = 0; k < 11; k++) {
/*
 *              If the block is populated, everybody participates in
 *              a broadcast to acquire the data.  All tasks other than
 *              task 0 must allocate the buffer into which the data will
 *              be received.  Block being broadcast consists of
 *              <numDeriv> tensors.
 */
                if (anisoTable->tagmap[i][j][k] == 1) {

                    if (home->myDomain != 0) {
                        anisoTable->Gderiv_map[i][j][k] =
                                (matrix *)malloc(sizeof(matrix) * numDeriv);
                    }

                    MPI_Bcast(anisoTable->Gderiv_map[i][j][k],
                              numDeriv*9, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                }
            }
        }
    }
#endif  /* ifdef PARALLEL */

    return;
}


void FMAnisotropicInit(Home_t *home)
{
    Param_t        *param = home->param;
    FMAnisoTable_t *anisoTable;

/*
 *  Have all tasks allocate the structure to hold the Greens function
 *  derivatives and associated values.  All value initially zero.
 */
    home->fmAnisoTable = (FMAnisoTable_t *)calloc(1, sizeof(FMAnisoTable_t));
    anisoTable = home->fmAnisoTable;

/*
 *  Only the master task reads the file of Greens function derivatives
 */
    if (home->myDomain == 0) {
        int    i, j;
        int    tblMaxOrder, numDeriv;
        real8  nrm, nrmDiff;
        real8  xScale, yScale, zScale;
        real8  relativeTol = 1.0e-5;
        real8  cellSize[3];
        real8  C66[6][6], C3333[3][3][3][3];
        FILE  *fpTable;

/*
 *      Open the Green's function derivatives file
 */
        fpTable = fopen(param->fmAnisoDerivTbl, "r");

        if (fpTable == (FILE *)NULL) {
            Fatal("FMAnisotropicInit: Open error %d on anistropy derivatives "
                  "file %s", errno, param->fmAnisoDerivTbl);
        }

/*
 *      Read the metadata from the derivatives file and save some of 
 *      the values for later use.
 */
        ReadDerivTableMetadata(fpTable, param->fmAnisoDerivTbl, C66,
                               &xScale, &yScale, &zScale, &tblMaxOrder);

        anisoTable->xScale = xScale;
        anisoTable->yScale = yScale;
        anisoTable->zScale = zScale;

        anisoTable->scaleVec[0] = xScale;
        anisoTable->scaleVec[1] = yScale;
        anisoTable->scaleVec[2] = zScale;

        anisoTable->derivTblMaxOrder = tblMaxOrder;
        numDeriv = (tblMaxOrder+3)*(tblMaxOrder+2)*(tblMaxOrder+1)/6;

/*
 *      Convert the 6x6 C matrix from the file into the 3x3x3x3 form
 */
        SetElasticConstantMatrix4D(C66, C3333);

/*
 *      Compare the 6x6 C matrix from the file with the C matrix
 *      being used in the current simulation to verify they are
 *      consistent.
 */
        nrm = 0.0;
        nrmDiff = 0.0;

        for (i = 0; i < 6; i++) {
            for (j = 0; j < 6; j++) {
                real8 simVal  = param->elasticConstantMatrix[i][j];
                real8 fileVal = C66[i][j] - simVal;

                nrm = nrm + simVal*simVal;
                nrmDiff = nrmDiff + fileVal*fileVal;
            }
        }

        nrm = sqrt(nrm);
        nrmDiff = sqrt(nrmDiff);

        if ((nrm == 0.0) || (fabs(nrmDiff) > (relativeTol * nrm))) {
            Fatal("FMAnisotropicInit: C matrix from %s does not match "
                  "C matrix for current simulation",
                  param->fmAnisoDerivTbl);
        }

/*
 *      Verify the correctness of the {x,y,z} scale...
 */
        cellSize[0] = (param->maxCoordinates[0] - param->minCoordinates[0]) /
                      param->nXcells;
        cellSize[1] = (param->maxCoordinates[1] - param->minCoordinates[1]) /
                      param->nYcells;
        cellSize[2] = (param->maxCoordinates[2] - param->minCoordinates[2]) /
                      param->nZcells;

        if ((fabs(xScale*cellSize[1]/(yScale*cellSize[0]) - 1.0) > 1e-6) ||
            (fabs(yScale*cellSize[2]/(zScale*cellSize[1]) - 1.0) > 1e-6)) {
            Fatal("FMAnisotropicInit: Coordinate scalings in Greens function "
                  "file %s\nare not commensurate with the current simulation "
                  "cell size", param->fmAnisoDerivTbl);
        }

        ReadDerivTable(home, fpTable);

        fclose(fpTable);

    }  /* end if (home->myDomain == 0) */

/*
 *  Task zero is the only one with the derivatives table, so it needs
 *  to distribute it to the rest of the tasks
 */
    DistributeDerivTable(home);

    return;
}
#endif  /* end ifdef ANISOTROPIC && TAYLORFMM */
