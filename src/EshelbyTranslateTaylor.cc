
#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"


#ifdef ESHELBY
/*---------------------------------------------------------------------------
 *
 *      Function:    EshelbyTranslateTaylor
 *      Description: Shift the center of the provided taylor expansion.
 *                   Used during the downward pass of the FM code.
 *
 *      Arguments:
 *          taylorSeries taylor expansion coefficients at old
 *                       expansion center
 *          oldCenter    Coordinates of old taylor expansion center
 *          newCenter    Coordinates of the new taylor expansion center
 *          taylorOrder  self-explanatory
 *          tPrime       Array in which to return the shifted taylor
 *                       coefficients to the caller
 *
 *-------------------------------------------------------------------------*/
void EshelbyTranslateTaylor(real8 *taylorSeries, real8 oldCenter[3],
                            real8 newCenter[3], int taylorOrder,
                            real8 *tPrime)
{
        int    i, j, m, n;
        int    rOffset, cOffset;
        int    numRows, numCols, tRows, tCols, matrixSize;
        int    mRows, mCols;
        int    submNumRows, submNumCols, submRowOffset;
        int    tpOffset, etaCnt, indexStart, indexEnd, seriesLen;
        real8  *matrix, *tMatrix, *etav;


        numRows = sequence[taylorOrder+3] - 4;
        numCols = sequence[taylorOrder+1];

        matrixSize = numRows * numCols * sizeof(real8);
        matrix = (real8 *)calloc(1, matrixSize);

        etaCnt = sequence[taylorOrder+1];
        etav = (real8 *)calloc(1, etaCnt * sizeof(real8));

        EshelbyMakeEta(1.0, 1.0, newCenter, oldCenter, taylorOrder, etav);

        if (taylorOrder > 1) {

            indexStart = sequence[2];
            indexEnd = sequence[taylorOrder+1] - 1;

            for (j = 0, i = indexStart; i < indexEnd; j++, i++) {
                etav[i] *= multFactor[j];
            }
        }

/*
 *      Create a temporary matrix of the maximum size and reuse
 *      it rather than reallocate over and over.
 */
        mRows = sequence[taylorOrder + 3] - sequence[taylorOrder + 2];
        mCols = sequence[taylorOrder + 1] - sequence[taylorOrder];

        tMatrix = (real8 *)calloc(1, mRows * mCols * sizeof(real8));

        for (i = 0; i <= taylorOrder; i++) {

            for (j = 0; j <= (taylorOrder - i); j++) {

                tRows = sequence[i+3] - sequence[i+2];
                tCols = sequence[j+1] - sequence[j];

                EshelbyBuildMatrix(taylorSeries, i, j, tMatrix, mRows, mCols);

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

        seriesLen = sequence[taylorOrder+3] - 4;

        for (i = 0; i < seriesLen; i++) {
            tPrime[i] = 0.0;
        }

/*
 *      Multiply the appropriate submatrix by etav and store the
 *      results in the correct portion of tPrime. (etav and
 *      tPrime are vectors, but for the multiply, we can treat
 *      them like N X 1 matrices where N == number of columns in
 *      the portion of matrix being multiplied by eta).
 */
        for (i = 0; i <= taylorOrder; i++) {

            submNumCols   = sequence[taylorOrder+1-i];
            submRowOffset = sequence[i+2] - 4;
            submNumRows   = (sequence[i+3] - 5) - submRowOffset + 1;
            tpOffset      = submRowOffset;

            MatrixMultArb(matrix, numCols, submRowOffset, 0,
                          submNumRows, submNumCols,
                          etav, 1, 0, 0, 1,
                          tPrime, 1, tpOffset, 0);

        }

/*
 *      Make sure we clean up...
 */
        free(etav);
        free(matrix);

        return;
}
#endif // ESHELBY
