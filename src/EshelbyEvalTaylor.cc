
#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"


#ifdef ESHELBY
/****************************************************************************
 *
 *      Function:     EshelbyEvalTaylor
 *      Description:  Compute the stress at a point due to the effects
 *                    of remote eshelby inclusions.
 *
 *      Arguments:
 *          taylorSeries  Array of taylor coefficients
 *          taylorCenter  Coordinates of the taylor expansion center
 *          position      Coordinates at which to evaluate the stress
 *          taylorOrder   order of the talor expansion
 *          stressMatrix  Array in which to return the stress components
 *                        to the caller
 *
 ***************************************************************************/
void EshelbyEvalTaylor(real8 *taylorSeries, real8 taylorCenter[3],
                       real8 position[3], int taylorOrder,
                       real8 stressMatrix[3][3])
{
        int    i, m, n, cIndex, offset, numVals, etaCnt;
        int    numRows, numCols, mRows, mCols, matrixSize;
        real8  *matrix, *etav, *tMatrix, *stress;

        numRows = sequence[taylorOrder+3] - 4;
        numCols = sequence[taylorOrder+1];

        matrixSize = numRows * numCols * sizeof(real8);
        matrix = (real8 *)calloc(1, matrixSize);

        etaCnt = sequence[taylorOrder+1];
        etav = (real8 *)calloc(1, etaCnt * sizeof(real8));

        EshelbyMakeEta(1.0, 1.0, taylorCenter, position, taylorOrder, etav);

/*
 *      Create a temporary matrix of the maximum size and reuse
 *      it rather than reallocate over and over.
 */
        mRows = sequence[3] - sequence[2];
        mCols = sequence[taylorOrder+1] - sequence[taylorOrder];

        tMatrix = (real8 *)calloc(1, mRows * mCols * sizeof(real8));

        for (i = 0; i <= taylorOrder; i++) {


            EshelbyBuildMatrix(taylorSeries, 0, i, tMatrix, mRows, mCols);

            for (m = 0; m < 6; m++) {
                cIndex = 0;
                for (n = sequence[i]; n < sequence[i+1]; n++) {
                    matrix[m*numCols+n] = tMatrix[m*mCols+cIndex];
                    cIndex++;
                }
            }
        }

        free(tMatrix);
        tMatrix = (real8 *)NULL;

        if (taylorOrder > 1) {
            numVals = sequence[taylorOrder+1] - sequence[2];
            offset = sequence[2];

            for (i = 0; i < numVals; i++) {
                etav[offset+i] *= multFactor[i];
            }
        }

        stress = (real8 *)malloc(numRows * sizeof(real8));
        MatrixMult(matrix, numRows, numCols, numCols,
                   etav, 1, 1, stress, 1);

        stressMatrix[0][0] = stress[0];
        stressMatrix[0][1] = stress[1];
        stressMatrix[0][2] = stress[2];
        stressMatrix[1][0] = stress[1];
        stressMatrix[1][1] = stress[3];
        stressMatrix[1][2] = stress[4];
        stressMatrix[2][0] = stress[2];
        stressMatrix[2][1] = stress[4];
        stressMatrix[2][2] = stress[5];

/*
 *      Make sure we clean up...
 */
        free(stress);
        free(etav);
        free(matrix);

        return;
}
#endif // ESHELBY
