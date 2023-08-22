
#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"


#ifdef ESHELBY
/****************************************************************************
 *
 *      Function:    EshelbyMakeTaylor
 *      Description: Calculate the eshleby inclusion taylor expansion
 *                   contribution for a cell from the multipole expansion of
 *                   a remote cell.
 *
 *      Arguments:
 *          targetPos      Coordinates of the center of the taylor expansion
 *          sourcePos      Coordinates of the center of the multipole
 *                         expansion
 *          eta            Array of multipole expansion coefficients
 *          multipoleOrder Self-explanatory
 *          taylorOrder    Also self-explanatory
 *          MU             Shear modulus
 *          NU             Pois ratio
 *          taylorSeries   Array in wich to return the taylor coefficients
 *                         to the caller.  Must be at least
 *                         sequence[taylorOrder+3] - 4 elements in length
 *
 ***************************************************************************/
void EshelbyMakeTaylor(real8 R[3],
                       real8 *eta, int multipoleOrder,
                       int taylorOrder, real8 MU, real8 NU,
                       real8 *taylorSeries)
{
        int   i, j, m, n;
        int   rOffset, cOffset, offset, numVals, seriesLen;
        int   numRows, numCols, tRows, tCols, matrixSize, drListSize; 
        int   mRows, mCols, etaCnt;
        real8 factor;
        real8 *matrix, *tMatrix, *etav, *drList;

    
        etaCnt = sequence[multipoleOrder+1];

        drListSize = sequence[multipoleOrder+taylorOrder+3] - 4;
        drList = (real8 *)calloc(1, drListSize * sizeof(real8));

        EshelbyDerofR(R, multipoleOrder+taylorOrder, drList);

        numRows = sequence[taylorOrder+3] - 4;
        numCols = sequence[multipoleOrder+1];

        matrixSize = numRows * numCols * sizeof(real8);
        matrix = (real8 *)calloc(1, matrixSize);

        factor = MU * (1.0+NU) / (1.0-(2.0*NU));

/*
 *      Create a temporary matrix of the maximum size and reuse
 *      it rather than reallocate over and over.
 */
        mRows = sequence[taylorOrder + 3] - sequence[taylorOrder + 2];
        mCols = sequence[multipoleOrder + 1] - sequence[multipoleOrder];

        tMatrix = (real8 *)calloc(1, mRows * mCols * sizeof(real8));

        for (i = 0; i <= taylorOrder; i++) {

            for (j = 0; j <= multipoleOrder; j++) {

                tRows = sequence[i+3] - sequence[i+2];
                tCols = sequence[j+1] - sequence[j];

                EshelbyBuildMatrix(drList, i, j, tMatrix, mRows, mCols);

                rOffset = sequence[i+2] - 4;
                cOffset = sequence[j];

                for (m = 0; m < tRows; m++) {
                    for (n = 0; n < tCols; n++) {
                        matrix[numCols*(rOffset+m)+(cOffset+n)] =
                                tMatrix[m*mCols+n];
                    }
                }
            }
        }

        free(tMatrix);
        tMatrix = (real8 *)NULL;

        etav = (real8 *)malloc(etaCnt * sizeof(real8));

        for (i = 0; i < etaCnt; i++) {
            etav[i] = eta[i];
        }

        if (multipoleOrder > 1) {

            numVals = sequence[multipoleOrder+1] - sequence[2];
            offset = sequence[2];

            for (i = 0; i < numVals; i++) {
                etav[offset+i] *= multFactor[i];
            }
        }

        seriesLen = sequence[taylorOrder+3] - 4;

        MatrixMult(matrix, numRows, numCols, numCols,
                   etav, 1, 1, taylorSeries, 1);

        for (i = 0; i < seriesLen; i++) {
            taylorSeries[i] *= factor;
        }

/*
 *      Make sure we clean up...
 */
        free(matrix);
        free(etav);
        free(drList);

        return;
}
#endif // ESHELBY
