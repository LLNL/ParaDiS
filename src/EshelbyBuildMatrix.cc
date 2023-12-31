
#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"

/****************************************************************************
 *
 *      Function:    BuildMatrix
 *
 *      Description: This function is specifically for use with certain
 *                   functions comprising the Esehlby FMM code.  It
 *                   generate a matrix from the provided coefficients
 *                   appropriate to the taylor and multipole expansion
 *                   orders specified.
 *
 *                   The size of the matrix data returned is:
 *                      rRows  = sequence[taylorOrder+3] -
 *                               sequence[taylorOrder+2]
 *                      rCols  = sequence[multipoleOrder+1] - 
 *                               sequence[multipoleOrder]
 *
 *      drList          Array of coefficients generated by previous call
 *                      to derofR().
 *      taylorOrder     taylor expansion order
 *      multipoleOrder  multipole expansion order
 *      matrix          pre-allocated storage in which to return
 *                      the generated matrix to the caller.  NOTE: the
 *                      incoming storage for <matrix> is permitted to be
 *                      of dimensions larger than needed, in which case
 *                      the function will set only rows 0 thru (rRows-1) and
 *                      columns 0 thru (rCols-1) where rRows and rCols
 *                      are defined as indicated above.
 *      mRows           Number of rows in the full <matrix> provided by
 *                      the caller
 *      mCols           Number of columns in the full <matrix> provided
 *                      by the caller.
 *
 ***************************************************************************/
void EshelbyBuildMatrix(double *drList, int taylorOrder, int multipoleOrder,
                        double *matrix, int mRows, int mCols)
{
        int  i, j, k, l;
        int  numRows, numCols;
        int  rowIndex, colIndex, seqIndex, drIndex, drOffset, matOffset;

        numRows = sequence[taylorOrder+3] - sequence[taylorOrder+2];
        numCols = sequence[multipoleOrder+1] - sequence[multipoleOrder];

        if ((numRows > mRows) || (numCols > mCols)) {
            Fatal("BuildMatrix: Provided storage too small!\n");
        }

/*
 *      Zero out the portion of the incoming matrix that will be
 *      populated with data
 */
        for (i = 0; i < numRows; i++) {
            for (j = 0; j < numCols; j++) {
                matOffset = i*mCols+j;
                matrix[matOffset] = 0.0;
            }
        }

        seqIndex = taylorOrder + multipoleOrder + 2;
        drIndex  = sequence[seqIndex] - 4;
        rowIndex = 0;

        for (i = 0; i < taylorOrder+3; i++) {
            for (j = 0; j <= i; j++) {
                colIndex = 0;
                drOffset = rowIndex;
                for (k = 0; k < multipoleOrder+1; k++) {
                    for (l = 0; l <= k; l++) {
                        matOffset = rowIndex*(mCols)+colIndex;
                        matrix[matOffset] = drList[drIndex+drOffset];
                        colIndex++;
                        drOffset++;
                    }
                    drOffset = drOffset+(i);
                }
                rowIndex++;
            }
        }

        return;
}
