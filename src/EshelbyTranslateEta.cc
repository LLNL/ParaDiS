
#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"

#define ELEM(a,b,c) ((a)*(c)+(b))


#ifdef ESHELBY
/****************************************************************************
 *
 *      Function:    EshelbyTranslateEta
 *      Description: Translate the multipole coefficients from one
 *                   expansion center to another.
 *
 *      Arguments:
 *          eta            Multipole coefficients associated with the
 *                         old expansion center
 *          newCenter      Coordiantes of the new expansion center
 *          oldCenter      Coordinates of the old expansion center
 *          multipoleOrder self-explanatory
 *          ep             Array in which to return the shifted multi-
 *                         pole coefficients to the caller
 * 
 ***************************************************************************/
void EshelbyTranslateEta(real8 *eta, real8 newCenter[3], real8 oldCenter[3],
                         int multipoleOrder, real8 *ep)
{
        int    i, nRows, nCols;
        int    etaSize, etaTransSize;
        real8  *etaTrans, *mat;


        etaSize = sequence[multipoleOrder+1];
        etaTransSize = sequence[multipoleOrder+1];
        etaTrans = (real8 *)calloc(1, etaTransSize * sizeof(real8));

        EshelbyMakeEta(1.0, 1.0, newCenter, oldCenter, multipoleOrder,
                       etaTrans);

/*
 *      construct matrix of intermingled terms
 *      set mat = eta*etatrans' which yields an
 *      etaSize X etaTransSize matrix.
 */
        mat = (real8 *)malloc(etaSize * etaTransSize * sizeof(real8));

        nRows = etaSize;
        nCols = etaTransSize;

        MatrixMult(eta, nRows, 1, 1,
                   etaTrans, nCols, nCols,
                   mat, nCols);

        if (multipoleOrder >= 0) {
            ep[0] = mat[0];
        }

        if (multipoleOrder > 0) {
            ep[1] = mat[ELEM(1,0,nCols)] + mat[ELEM(0,1,nCols)];
            ep[2] = mat[ELEM(2,0,nCols)] + mat[ELEM(0,2,nCols)];
            ep[3] = mat[ELEM(3,0,nCols)] + mat[ELEM(0,3,nCols)];
        }


        if (multipoleOrder > 1) {

            for (i = 4; i < 10; i++) {
                ep[i] = mat[ELEM(i,0,nCols)] + mat[ELEM(0,i,nCols)];
            }

            ep[4] += 2.0 * mat[ELEM(1,1,nCols)];
            ep[5] += mat[ELEM(1,2,nCols)] + mat[ELEM(2,1,nCols)];
            ep[6] += mat[ELEM(1,3,nCols)] + mat[ELEM(3,1,nCols)];
            ep[7] += 2.0 * mat[ELEM(2,2,nCols)];
            ep[8] += mat[ELEM(2,3,nCols)] + mat[ELEM(3,2,nCols)];
            ep[9] += 2.0 * mat[ELEM(3,3,nCols)];
        }


        if (multipoleOrder > 2) {

            for (i = 10; i < 20; i++) {
                ep[i] = mat[ELEM(i,0,nCols)] + mat[ELEM(0,i,nCols)];
            }

            /* 1 1 1 */
            ep[10] += 3.0 * (mat[ELEM(1,4,nCols)] + mat[ELEM(4,1,nCols)]);

            /* 1 1 2 */
            ep[11] += 2.0 * (mat[ELEM(1,5,nCols)] + mat[ELEM(5,1,nCols)]) +
                             mat[ELEM(2,4,nCols)] + mat[ELEM(4,2,nCols)];

            /* 1 1 3 */
            ep[12] += 2.0 * (mat[ELEM(1,6,nCols)] + mat[ELEM(6,1,nCols)]) +
                             mat[ELEM(3,4,nCols)] + mat[ELEM(4,3,nCols)];

            /* 1 2 2 */
            ep[13] +=        mat[ELEM(1,7,nCols)] + mat[ELEM(7,1,nCols)] +
                      2.0 * (mat[ELEM(2,5,nCols)] + mat[ELEM(5,2,nCols)]);

            /* 1 2 3 */
            ep[14] +=        mat[ELEM(1,8,nCols)] + mat[ELEM(8,1,nCols)] +
                             mat[ELEM(2,6,nCols)] + mat[ELEM(6,2,nCols)] +
                             mat[ELEM(3,5,nCols)] + mat[ELEM(5,3,nCols)];

            /* 1 3 3 */
            ep[15] +=        mat[ELEM(1,9,nCols)] + mat[ELEM(9,1,nCols)] +
                      2.0 * (mat[ELEM(3,6,nCols)] + mat[ELEM(6,3,nCols)]);

            /* 2 2 2 */
            ep[16] += 3.0 * (mat[ELEM(2,7,nCols)] + mat[ELEM(7,2,nCols)]);

            /* 2 2 3 */
            ep[17] += 2.0 * (mat[ELEM(2,8,nCols)] + mat[ELEM(8,2,nCols)]) +
                             mat[ELEM(3,7,nCols)] + mat[ELEM(7,3,nCols)];

            /* 2 3 3 */
            ep[18] +=        mat[ELEM(2,9,nCols)] + mat[ELEM(9,2,nCols)] +
                      2.0 * (mat[ELEM(3,8,nCols)] + mat[ELEM(8,3,nCols)]);

            /* 3 3 3 */
            ep[19] += 3.0 * (mat[ELEM(3,9,nCols)] + mat[ELEM(9,3,nCols)]);
        }




        if (multipoleOrder > 3) {

            for (i = 20; i < 35; i++) {
                ep[i] = mat[ELEM(i,0,nCols)] + mat[ELEM(0,i,nCols)];
            }

            /* 1 1 1 1 */
            ep[20] += 4.0 * (mat[ELEM(1,10,nCols)] + mat[ELEM(10,1,nCols)]) +
                      6.0 *  mat[ELEM(4,4,nCols)];

            /* 1 1 1 2 */
            ep[21] += 3.0 * (mat[ELEM(1,11,nCols)] + mat[ELEM(11,1,nCols)]) +
                             mat[ELEM(2,10,nCols)] + mat[ELEM(10,2,nCols)]  +
                      3.0 * (mat[ELEM(4,5,nCols)]  + mat[ELEM(5,4,nCols)]);

            /* 1 1 1 3 */
            ep[22] += 3.0 * (mat[ELEM(1,12,nCols)] + mat[ELEM(12,1,nCols)]) +
                             mat[ELEM(3,10,nCols)] + mat[ELEM(10,3,nCols)]  +
                      3.0 * (mat[ELEM(4,6,nCols)]  + mat[ELEM(6,4,nCols)]);

            /* 1 1 2 2 */
            ep[23] += 2.0 * (mat[ELEM(1,13,nCols)] + mat[ELEM(13,1,nCols)]) +
                      2.0 * (mat[ELEM(2,11,nCols)] + mat[ELEM(11,2,nCols)]) +
                             mat[ELEM(4,7,nCols)]  + mat[ELEM(7,4,nCols)]   +
                      4.0 *  mat[ELEM(5,5,nCols)];

            /* 1 1 2 3 */
            ep[24] += 2.0 * (mat[ELEM(1,14,nCols)] + mat[ELEM(14,1,nCols)]) +
                             mat[ELEM(2,12,nCols)] + mat[ELEM(12,2,nCols)]  +
                             mat[ELEM(3,11,nCols)] + mat[ELEM(11,3,nCols)]  +
                             mat[ELEM(4,8,nCols)]  + mat[ELEM(8,4,nCols)]   +
                      2.0 * (mat[ELEM(5,6,nCols)]  + mat[ELEM(6,5,nCols)]);

            /* 1 1 3 3 */
            ep[25] += 2.0 * (mat[ELEM(1,15,nCols)] + mat[ELEM(15,1,nCols)]) +
                      2.0 * (mat[ELEM(3,12,nCols)] + mat[ELEM(12,3,nCols)]) +
                             mat[ELEM(4,9,nCols)]  + mat[ELEM(9,4,nCols)]   +
                      4.0 *  mat[ELEM(6,6,nCols)];

            /* 1 2 2 2 */
            ep[26] +=        mat[ELEM(1,16,nCols)] + mat[ELEM(16,1,nCols)]  +
                      3.0 * (mat[ELEM(2,13,nCols)] + mat[ELEM(13,2,nCols)]) +
                      3.0 * (mat[ELEM(5,7,nCols)]  + mat[ELEM(7,5,nCols)]);

            /* 1 2 2 3 */
            ep[27] +=        mat[ELEM(1,17,nCols)] + mat[ELEM(17,1,nCols)]  +
                      2.0 * (mat[ELEM(2,14,nCols)] + mat[ELEM(14,2,nCols)]) +
                             mat[ELEM(3,13,nCols)] + mat[ELEM(13,3,nCols)]  +
                      2.0 * (mat[ELEM(5,8,nCols)]  + mat[ELEM(8,5,nCols)])  +
                             mat[ELEM(6,7,nCols)]  + mat[ELEM(7,6,nCols)];

            /* 1 2 3 3 */
            ep[28] +=        mat[ELEM(1,18,nCols)] + mat[ELEM(18,1,nCols)]  +
                             mat[ELEM(2,15,nCols)] + mat[ELEM(15,2,nCols)]  +
                      2.0 * (mat[ELEM(3,14,nCols)] + mat[ELEM(14,3,nCols)]) +
                             mat[ELEM(5,9,nCols)]  + mat[ELEM(9,5,nCols)]   +
                      2.0 * (mat[ELEM(6,8,nCols)]  + mat[ELEM(8,6,nCols)]);

            /* 1 3 3 3 */
            ep[29] +=        mat[ELEM(1,19,nCols)] + mat[ELEM(19,1,nCols)]  +
                      3.0 * (mat[ELEM(3,15,nCols)] + mat[ELEM(15,3,nCols)]) +
                      3.0 * (mat[ELEM(6,9,nCols)]  + mat[ELEM(9,6,nCols)]);

            /* 2 2 2 2 */
            ep[30] += 4.0 * (mat[ELEM(2,16,nCols)] + mat[ELEM(16,2,nCols)]) +
                      6.0 *  mat[ELEM(7,7,nCols)];

            /* 2 2 2 3 */
            ep[31] += 3.0 * (mat[ELEM(2,17,nCols)] + mat[ELEM(17,2,nCols)]) +
                             mat[ELEM(3,16,nCols)] + mat[ELEM(16,3,nCols)]  +
                      3.0 * (mat[ELEM(7,8,nCols)]  + mat[ELEM(8,7,nCols)]);

            /* 2 2 3 3 */
            ep[32] += 2.0 * (mat[ELEM(2,18,nCols)] + mat[ELEM(18,2,nCols)]) +
                      2.0 * (mat[ELEM(3,17,nCols)] + mat[ELEM(17,3,nCols)]) +
                             mat[ELEM(7,9,nCols)]  + mat[ELEM(9,7,nCols)]   +
                      4.0 *  mat[ELEM(8,8,nCols)];

            /* 2 3 3 3 */
            ep[33] +=        mat[ELEM(2,19,nCols)] + mat[ELEM(19,2,nCols)]  +
                      3.0 * (mat[ELEM(3,18,nCols)] + mat[ELEM(18,3,nCols)]) +
                      3.0 * (mat[ELEM(8,9,nCols)]  + mat[ELEM(9,8,nCols)]);

            /* 3 3 3 3 */
            ep[34] += 4.0 * (mat[ELEM(3,19,nCols)] + mat[ELEM(19,3,nCols)]) +
                      6.0 *  mat[ELEM(9,9,nCols)];

        }

 
        if (multipoleOrder > 4) {

            for (i = 35; i < 56; i++) {
                ep[i] = mat[ELEM(i,0,nCols)] + mat[ELEM(0,i,nCols)];
            }

            /* 1 1 1 1 1 */
            ep[35] +=  5.0 * (mat[ELEM(1,20,nCols)] + mat[ELEM(20,1,nCols)]) +
                      10.0 * (mat[ELEM(4,10,nCols)] + mat[ELEM(10,4,nCols)]);

            /* 1 1 1 1 2 */
            ep[36] +=  4.0 * (mat[ELEM(1,21,nCols)] + mat[ELEM(21,1,nCols)]) +
                              mat[ELEM(2,20,nCols)] + mat[ELEM(20,2,nCols)]  +
                       6.0 * (mat[ELEM(4,11,nCols)] + mat[ELEM(11,4,nCols)]) +
                       4.0 * (mat[ELEM(5,10,nCols)] + mat[ELEM(10,5,nCols)]);

            /* 1 1 1 1 3 */
            ep[37] +=  4.0 * (mat[ELEM(1,22,nCols)] + mat[ELEM(22,1,nCols)]) +
                              mat[ELEM(3,20,nCols)] + mat[ELEM(20,3,nCols)]  +
                       6.0 * (mat[ELEM(4,12,nCols)] + mat[ELEM(12,4,nCols)]) +
                       4.0 * (mat[ELEM(6,10,nCols)] + mat[ELEM(10,6,nCols)]);

            /* 1 1 1 2 2 */
            ep[38] +=  3.0 * (mat[ELEM(1,23,nCols)] + mat[ELEM(23,1,nCols)]) +
                       2.0 * (mat[ELEM(2,21,nCols)] + mat[ELEM(21,2,nCols)]) +
                       3.0 * (mat[ELEM(4,13,nCols)] + mat[ELEM(13,4,nCols)]) +
                       6.0 * (mat[ELEM(5,11,nCols)] + mat[ELEM(11,5,nCols)]) +
                              mat[ELEM(7,10,nCols)] + mat[ELEM(10,7,nCols)];

            /* 1 1 1 2 3 */
            ep[39] +=  3.0 * (mat[ELEM(1,24,nCols)] + mat[ELEM(24,1,nCols)]) +
                              mat[ELEM(2,22,nCols)] + mat[ELEM(22,2,nCols)]  +
                              mat[ELEM(3,21,nCols)] + mat[ELEM(21,3,nCols)]  +
                       3.0 * (mat[ELEM(4,14,nCols)] + mat[ELEM(14,4,nCols)]) +
                       3.0 * (mat[ELEM(5,12,nCols)] + mat[ELEM(12,5,nCols)]) +
                       3.0 * (mat[ELEM(6,11,nCols)] + mat[ELEM(11,6,nCols)]) +
                              mat[ELEM(8,10,nCols)] + mat[ELEM(10,8,nCols)];

            /* 1 1 1 3 3 */
            ep[40] +=  3.0 * (mat[ELEM(1,25,nCols)] + mat[ELEM(25,1,nCols)]) +
                       2.0 * (mat[ELEM(3,22,nCols)] + mat[ELEM(22,3,nCols)]) +
                       3.0 * (mat[ELEM(4,15,nCols)] + mat[ELEM(15,4,nCols)]) +
                       6.0 * (mat[ELEM(6,12,nCols)] + mat[ELEM(12,6,nCols)]) +
                              mat[ELEM(9,10,nCols)] + mat[ELEM(10,9,nCols)];

            /* 1 1 2 2 2 */
            ep[41] +=  2.0 * (mat[ELEM(1,26,nCols)] + mat[ELEM(26,1,nCols)]) +
                       3.0 * (mat[ELEM(2,23,nCols)] + mat[ELEM(23,2,nCols)]) +
                              mat[ELEM(4,16,nCols)] + mat[ELEM(16,6,nCols)]  +
                       6.0 * (mat[ELEM(5,13,nCols)] + mat[ELEM(13,5,nCols)]) +
                       3.0 * (mat[ELEM(7,11,nCols)] + mat[ELEM(11,7,nCols)]);

            /* 1 1 2 2 3 */
            ep[42] +=  2.0 * (mat[ELEM(1,27,nCols)] + mat[ELEM(27,1,nCols)]) +
                       2.0 * (mat[ELEM(2,24,nCols)] + mat[ELEM(24,2,nCols)]) +
                              mat[ELEM(3,23,nCols)] + mat[ELEM(23,3,nCols)]  +
                              mat[ELEM(4,17,nCols)] + mat[ELEM(17,4,nCols)]  +
                       4.0 * (mat[ELEM(5,14,nCols)] + mat[ELEM(14,5,nCols)]) +
                       2.0 * (mat[ELEM(6,13,nCols)] + mat[ELEM(13,6,nCols)]) +
                              mat[ELEM(7,12,nCols)] + mat[ELEM(12,7,nCols)]  +
                       2.0 * (mat[ELEM(8,11,nCols)] + mat[ELEM(11,8,nCols)]);

            /* 1 1 2 3 3 */
            ep[43] +=  2.0 * (mat[ELEM(1,28,nCols)] + mat[ELEM(28,1,nCols)]) +
                              mat[ELEM(2,25,nCols)] + mat[ELEM(25,2,nCols)]  +
                       2.0 * (mat[ELEM(3,24,nCols)] + mat[ELEM(24,3,nCols)]) +
                              mat[ELEM(4,18,nCols)] + mat[ELEM(18,4,nCols)]  +
                       2.0 * (mat[ELEM(5,15,nCols)] + mat[ELEM(15,5,nCols)]) +
                       4.0 * (mat[ELEM(6,14,nCols)] + mat[ELEM(14,6,nCols)]) +
                       2.0 * (mat[ELEM(8,12,nCols)] + mat[ELEM(12,8,nCols)]) +
                              mat[ELEM(9,11,nCols)] + mat[ELEM(11,9,nCols)];

            /* 1 1 3 3 3 */
            ep[44] +=  2.0 * (mat[ELEM(1,29,nCols)] + mat[ELEM(29,1,nCols)]) +
                       3.0 * (mat[ELEM(3,25,nCols)] + mat[ELEM(25,3,nCols)]) +
                              mat[ELEM(4,19,nCols)] + mat[ELEM(19,4,nCols)]  +
                       6.0 * (mat[ELEM(6,15,nCols)] + mat[ELEM(15,6,nCols)]) +
                       3.0 * (mat[ELEM(9,12,nCols)] + mat[ELEM(12,9,nCols)]);

            /* 1 2 2 2 2 */
            ep[45] +=         mat[ELEM(1,30,nCols)] + mat[ELEM(30,1,nCols)]  +
                       4.0 * (mat[ELEM(2,26,nCols)] + mat[ELEM(26,2,nCols)]) +
                       4.0 * (mat[ELEM(5,16,nCols)] + mat[ELEM(16,5,nCols)]) +
                       6.0 * (mat[ELEM(7,13,nCols)] + mat[ELEM(13,7,nCols)]);

            /* 1 2 2 2 3 */
            ep[46] +=         mat[ELEM(1,31,nCols)] + mat[ELEM(31,1,nCols)]  +
                       3.0 * (mat[ELEM(2,27,nCols)] + mat[ELEM(27,2,nCols)]) +
                              mat[ELEM(3,26,nCols)] + mat[ELEM(26,3,nCols)]  +
                       3.0 * (mat[ELEM(5,17,nCols)] + mat[ELEM(17,5,nCols)]) +
                              mat[ELEM(6,16,nCols)] + mat[ELEM(16,6,nCols)]  +
                       3.0 * (mat[ELEM(7,14,nCols)] + mat[ELEM(14,7,nCols)]) +
                       3.0 * (mat[ELEM(8,13,nCols)] + mat[ELEM(13,8,nCols)]);

            /* 1 2 2 3 3 */
            ep[47] +=         mat[ELEM(1,32,nCols)] + mat[ELEM(32,1,nCols)]  +
                       2.0 * (mat[ELEM(2,28,nCols)] + mat[ELEM(28,2,nCols)]) +
                       2.0 * (mat[ELEM(3,27,nCols)] + mat[ELEM(27,3,nCols)]) +
                       2.0 * (mat[ELEM(5,18,nCols)] + mat[ELEM(18,5,nCols)]) +
                       2.0 * (mat[ELEM(6,17,nCols)] + mat[ELEM(17,6,nCols)]) +
                              mat[ELEM(7,15,nCols)] + mat[ELEM(15,7,nCols)]  +
                       4.0 * (mat[ELEM(8,14,nCols)] + mat[ELEM(14,8,nCols)]) +
                              mat[ELEM(9,13,nCols)] + mat[ELEM(13,9,nCols)];

            /* 1 2 3 3 3 */
            ep[48] +=         mat[ELEM(1,33,nCols)] + mat[ELEM(33,1,nCols)]  +
                              mat[ELEM(2,29,nCols)] + mat[ELEM(29,2,nCols)]  +
                       3.0 * (mat[ELEM(3,28,nCols)] + mat[ELEM(28,3,nCols)]) +
                              mat[ELEM(5,19,nCols)] + mat[ELEM(19,5,nCols)]  +
                       3.0 * (mat[ELEM(6,18,nCols)] + mat[ELEM(18,6,nCols)]) +
                       3.0 * (mat[ELEM(8,15,nCols)] + mat[ELEM(15,8,nCols)]) +
                       3.0 * (mat[ELEM(9,14,nCols)] + mat[ELEM(14,9,nCols)]);

            /* 1 3 3 3 3 */
            ep[49] +=         mat[ELEM(1,34,nCols)] + mat[ELEM(34,1,nCols)]  +
                       4.0 * (mat[ELEM(3,29,nCols)] + mat[ELEM(29,3,nCols)]) +
                       4.0 * (mat[ELEM(6,19,nCols)] + mat[ELEM(19,6,nCols)]) +
                       6.0 * (mat[ELEM(9,15,nCols)] + mat[ELEM(15,9,nCols)]);

            /* 2 2 2 2 2 */
            ep[50] +=  5.0 * (mat[ELEM(2,30,nCols)] + mat[ELEM(30,2,nCols)]) +
                      10.0 * (mat[ELEM(7,16,nCols)] + mat[ELEM(16,7,nCols)]);

            /* 2 2 2 2 3 */
            ep[51] +=  4.0 * (mat[ELEM(2,31,nCols)] + mat[ELEM(31,2,nCols)]) +
                              mat[ELEM(3,30,nCols)] + mat[ELEM(30,3,nCols)]  +
                       6.0 * (mat[ELEM(7,17,nCols)] + mat[ELEM(17,7,nCols)]) +
                       4.0 * (mat[ELEM(8,16,nCols)] + mat[ELEM(16,8,nCols)]);

            /* 2 2 2 3 3 */
            ep[52] +=  3.0 * (mat[ELEM(2,32,nCols)] + mat[ELEM(32,2,nCols)]) +
                       2.0 * (mat[ELEM(3,31,nCols)] + mat[ELEM(31,3,nCols)]) +
                       3.0 * (mat[ELEM(7,18,nCols)] + mat[ELEM(18,7,nCols)]) +
                       6.0 * (mat[ELEM(8,17,nCols)] + mat[ELEM(17,8,nCols)]) +
                              mat[ELEM(9,16,nCols)] + mat[ELEM(16,9,nCols)];

            /* 2 2 3 3 3 */
            ep[53] +=  2.0 * (mat[ELEM(2,33,nCols)] + mat[ELEM(33,2,nCols)]) +
                       3.0 * (mat[ELEM(3,32,nCols)] + mat[ELEM(32,3,nCols)]) +
                              mat[ELEM(7,19,nCols)] + mat[ELEM(19,7,nCols)]  +
                       6.0 * (mat[ELEM(8,18,nCols)] + mat[ELEM(18,8,nCols)]) +
                       3.0 * (mat[ELEM(9,17,nCols)] + mat[ELEM(17,9,nCols)]);

            /* 2 3 3 3 3 */
            ep[54] +=         mat[ELEM(2,34,nCols)] + mat[ELEM(34,2,nCols)]  +
                       4.0 * (mat[ELEM(3,33,nCols)] + mat[ELEM(33,3,nCols)]) +
                       4.0 * (mat[ELEM(8,19,nCols)] + mat[ELEM(19,8,nCols)]) +
                       6.0 * (mat[ELEM(9,18,nCols)] + mat[ELEM(18,9,nCols)]);

            /* 3 3 3 3 3 */
            ep[55] +=  5.0 * (mat[ELEM(3,34,nCols)] + mat[ELEM(34,3,nCols)]) +
                      10.0 * (mat[ELEM(9,19,nCols)] + mat[ELEM(19,9,nCols)]);

        }


    
        if (multipoleOrder > 5) {

            for (i = 56; i < 84; i++) {
                ep[i] = mat[ELEM(i,0,nCols)] + mat[ELEM(0,i,nCols)];
            }

            /* 1 1 1 1 1 1 */
            ep[56] +=  6.0 * (mat[ELEM(1,35,nCols)] + mat[ELEM(35,1,nCols)])  +
                      15.0 * (mat[ELEM(4,20,nCols)] + mat[ELEM(20,4,nCols)])  +
                      20.0 *  mat[ELEM(10,10,nCols)];

            /* 1 1 1 1 1 2 */
            ep[57] +=  5.0 * (mat[ELEM(1,36,nCols)] + mat[ELEM(36,1,nCols)])  +
                              mat[ELEM(2,35,nCols)] + mat[ELEM(35,2,nCols)]   +
                      10.0 * (mat[ELEM(4,21,nCols)] + mat[ELEM(21,4,nCols)])  +
                       5.0 * (mat[ELEM(5,20,nCols)] + mat[ELEM(20,5,nCols)])  +
                      10.0 * (mat[ELEM(10,11,nCols)] + mat[ELEM(11,10,nCols)]);

            /* 1 1 1 1 1 3 */
            ep[58] +=  5.0 * (mat[ELEM(1,37,nCols)] + mat[ELEM(37,1,nCols)])  +
                              mat[ELEM(3,35,nCols)] + mat[ELEM(35,3,nCols)]   +
                      10.0 * (mat[ELEM(4,22,nCols)] + mat[ELEM(22,4,nCols)])  +
                       5.0 * (mat[ELEM(6,20,nCols)] + mat[ELEM(20,6,nCols)])  +
                      10.0 * (mat[ELEM(10,12,nCols)] + mat[ELEM(12,10,nCols)]);

            /* 1 1 1 1 2 2 */
            ep[59] +=  4.0 * (mat[ELEM(1,38,nCols)] + mat[ELEM(38,1,nCols)])  +
                       2.0 * (mat[ELEM(2,36,nCols)] + mat[ELEM(36,2,nCols)])  +
                       6.0 * (mat[ELEM(4,23,nCols)] + mat[ELEM(23,4,nCols)])  +
                       8.0 * (mat[ELEM(5,21,nCols)] + mat[ELEM(21,5,nCols)])  +
                              mat[ELEM(7,20,nCols)] + mat[ELEM(20,7,nCols)]   +
                       4.0 * (mat[ELEM(10,13,nCols)] + mat[ELEM(13,10,nCols)])+
                      12.0 *  mat[ELEM(11,11,nCols)];

            /* 1 1 1 1 2 3 */
            ep[60] +=  4.0 * (mat[ELEM(1,39,nCols)] + mat[ELEM(39,1,nCols)])  +
                              mat[ELEM(2,37,nCols)] + mat[ELEM(37,2,nCols)]   +
                              mat[ELEM(3,36,nCols)] + mat[ELEM(36,3,nCols)]   +
                       6.0 * (mat[ELEM(4,24,nCols)] + mat[ELEM(24,4,nCols)])  +
                       4.0 * (mat[ELEM(5,22,nCols)] + mat[ELEM(22,5,nCols)])  +
                       4.0 * (mat[ELEM(6,21,nCols)] + mat[ELEM(21,6,nCols)])  +
                              mat[ELEM(8,20,nCols)] + mat[ELEM(20,8,nCols)]   +
                       4.0 * (mat[ELEM(10,14,nCols)] + mat[ELEM(14,10,nCols)])+
                       6.0 * (mat[ELEM(11,12,nCols)] + mat[ELEM(12,11,nCols)]);

            /* 1 1 1 1 3 3 */
            ep[61] +=  4.0 * (mat[ELEM(1,40,nCols)] + mat[ELEM(40,1,nCols)])  +
                       2.0 * (mat[ELEM(3,37,nCols)] + mat[ELEM(37,3,nCols)])  +
                       6.0 * (mat[ELEM(4,25,nCols)] + mat[ELEM(25,4,nCols)])  +
                       8.0 * (mat[ELEM(6,22,nCols)] + mat[ELEM(22,6,nCols)])  +
                              mat[ELEM(9,20,nCols)] + mat[ELEM(20,9,nCols)]   +
                       4.0 * (mat[ELEM(10,15,nCols)] + mat[ELEM(15,10,nCols)])+
                      12.0 *  mat[ELEM(12,12,nCols)];

            /* 1 1 1 2 2 2 */
            ep[62] +=  3.0 * (mat[ELEM(1,41,nCols)] + mat[ELEM(41,1,nCols)])  +
                       3.0 * (mat[ELEM(2,38,nCols)] + mat[ELEM(38,2,nCols)])  +
                       3.0 * (mat[ELEM(4,26,nCols)] + mat[ELEM(26,4,nCols)])  +
                       9.0 * (mat[ELEM(5,23,nCols)] + mat[ELEM(23,5,nCols)])  +
                       3.0 * (mat[ELEM(7,21,nCols)] + mat[ELEM(21,7,nCols)])  +
                              mat[ELEM(10,16,nCols)] + mat[ELEM(16,10,nCols)] +
                       9.0 * (mat[ELEM(11,13,nCols)] + mat[ELEM(13,11,nCols)]);

            /* 1 1 1 2 2 3 */
            ep[63] +=  3.0 * (mat[ELEM(1,42,nCols)] + mat[ELEM(42,1,nCols)])  +
                       2.0 * (mat[ELEM(2,39,nCols)] + mat[ELEM(39,2,nCols)])  +
                              mat[ELEM(3,38,nCols)] + mat[ELEM(38,3,nCols)]   +
                       3.0 * (mat[ELEM(4,27,nCols)] + mat[ELEM(27,4,nCols)])  +
                       6.0 * (mat[ELEM(5,24,nCols)] + mat[ELEM(24,5,nCols)])  +
                       3.0 * (mat[ELEM(6,23,nCols)] + mat[ELEM(23,6,nCols)])  +
                              mat[ELEM(7,22,nCols)] + mat[ELEM(22,7,nCols)]   +
                       2.0 * (mat[ELEM(8,21,nCols)] + mat[ELEM(21,8,nCols)])  +
                              mat[ELEM(10,17,nCols)] + mat[ELEM(17,10,nCols)] +
                       6.0 * (mat[ELEM(11,14,nCols)] + mat[ELEM(14,11,nCols)])+
                       3.0 * (mat[ELEM(12,13,nCols)] + mat[ELEM(13,12,nCols)]);

            /* 1 1 1 2 3 3 */
            ep[64] +=  3.0 * (mat[ELEM(1,43,nCols)] + mat[ELEM(43,1,nCols)])  +
                              mat[ELEM(2,40,nCols)] + mat[ELEM(40,2,nCols)]   +
                       2.0 * (mat[ELEM(3,39,nCols)] + mat[ELEM(39,3,nCols)])  +
                       3.0 * (mat[ELEM(4,28,nCols)] + mat[ELEM(28,4,nCols)])  +
                       3.0 * (mat[ELEM(5,25,nCols)] + mat[ELEM(25,5,nCols)])  +
                       6.0 * (mat[ELEM(6,24,nCols)] + mat[ELEM(24,6,nCols)])  +
                       2.0 * (mat[ELEM(8,22,nCols)] + mat[ELEM(22,8,nCols)])  +
                              mat[ELEM(9,21,nCols)] + mat[ELEM(21,9,nCols)]   +
                              mat[ELEM(10,18,nCols)] + mat[ELEM(18,10,nCols)] +
                       3.0 * (mat[ELEM(11,15,nCols)] + mat[ELEM(15,11,nCols)])+
                       6.0 * (mat[ELEM(12,14,nCols)] + mat[ELEM(14,12,nCols)]);

            /* 1 1 1 3 3 3 */
            ep[65] +=  3.0 * (mat[ELEM(1,44,nCols)] + mat[ELEM(44,1,nCols)])  +
                       3.0 * (mat[ELEM(3,40,nCols)] + mat[ELEM(40,3,nCols)])  +
                       3.0 * (mat[ELEM(4,29,nCols)] + mat[ELEM(29,4,nCols)])  +
                       9.0 * (mat[ELEM(6,25,nCols)] + mat[ELEM(25,6,nCols)])  +
                       3.0 * (mat[ELEM(9,22,nCols)] + mat[ELEM(22,9,nCols)])  +
                              mat[ELEM(10,19,nCols)] + mat[ELEM(19,10,nCols)] +
                       9.0 * (mat[ELEM(12,15,nCols)] + mat[ELEM(15,12,nCols)]);

            /* 1 1 2 2 2 2 */
            ep[66] +=  2.0 * (mat[ELEM(1,45,nCols)] + mat[ELEM(45,1,nCols)])  +
                       4.0 * (mat[ELEM(2,41,nCols)] + mat[ELEM(41,2,nCols)])  +
                              mat[ELEM(4,30,nCols)] + mat[ELEM(30,4,nCols)]   +
                       8.0 * (mat[ELEM(5,26,nCols)] + mat[ELEM(26,5,nCols)])  +
                       6.0 * (mat[ELEM(7,23,nCols)] + mat[ELEM(23,7,nCols)])  +
                       4.0 * (mat[ELEM(11,16,nCols)] + mat[ELEM(16,11,nCols)])+
                      12.0 *  mat[ELEM(13,13,nCols)];

            /* 1 1 2 2 2 3 */
            ep[67] +=  2.0 * (mat[ELEM(1,46,nCols)] + mat[ELEM(46,1,nCols)])  +
                       3.0 * (mat[ELEM(2,42,nCols)] + mat[ELEM(42,2,nCols)])  +
                              mat[ELEM(3,41,nCols)] + mat[ELEM(41,3,nCols)]   +
                              mat[ELEM(4,31,nCols)] + mat[ELEM(31,4,nCols)]   +
                       6.0 * (mat[ELEM(5,27,nCols)] + mat[ELEM(27,5,nCols)])  +
                       2.0 * (mat[ELEM(6,26,nCols)] + mat[ELEM(26,6,nCols)])  +
                       3.0 * (mat[ELEM(7,24,nCols)] + mat[ELEM(24,7,nCols)])  +
                       3.0 * (mat[ELEM(8,23,nCols)] + mat[ELEM(23,8,nCols)])  +
                       3.0 * (mat[ELEM(11,17,nCols)] + mat[ELEM(17,11,nCols)])+
                              mat[ELEM(12,16,nCols)] + mat[ELEM(16,12,nCols)] +
                       6.0 * (mat[ELEM(13,14,nCols)] + mat[ELEM(14,13,nCols)]);

            /* 1 1 2 2 3 3 */
            ep[68] +=  2.0 * (mat[ELEM(1,47,nCols)] + mat[ELEM(47,1,nCols)])  +
                       2.0 * (mat[ELEM(2,43,nCols)] + mat[ELEM(43,2,nCols)])  +
                       2.0 * (mat[ELEM(3,42,nCols)] + mat[ELEM(42,3,nCols)])  +
                              mat[ELEM(4,32,nCols)] + mat[ELEM(32,4,nCols)]   +
                       4.0 * (mat[ELEM(5,28,nCols)] + mat[ELEM(28,5,nCols)])  +
                       4.0 * (mat[ELEM(6,27,nCols)] + mat[ELEM(27,6,nCols)])  +
                              mat[ELEM(7,25,nCols)] + mat[ELEM(25,7,nCols)]   +
                       4.0 * (mat[ELEM(8,24,nCols)] + mat[ELEM(24,8,nCols)])  +
                              mat[ELEM(9,23,nCols)] + mat[ELEM(23,9,nCols)]   +
                       2.0 * (mat[ELEM(11,18,nCols)] + mat[ELEM(18,11,nCols)])+
                       2.0 * (mat[ELEM(12,17,nCols)] + mat[ELEM(17,12,nCols)])+
                       2.0 * (mat[ELEM(13,15,nCols)] + mat[ELEM(15,13,nCols)])+
                       8.0 *  mat[ELEM(14,14,nCols)];

            /* 1 1 2 3 3 3 */
            ep[69] +=  2.0 * (mat[ELEM(1,48,nCols)] + mat[ELEM(48,1,nCols)])  +
                              mat[ELEM(2,44,nCols)] + mat[ELEM(44,2,nCols)]   +
                       3.0 * (mat[ELEM(3,43,nCols)] + mat[ELEM(43,3,nCols)])  +
                              mat[ELEM(4,33,nCols)] + mat[ELEM(33,4,nCols)]   +
                       2.0 * (mat[ELEM(5,29,nCols)] + mat[ELEM(29,5,nCols)])  +
                       6.0 * (mat[ELEM(6,28,nCols)] + mat[ELEM(28,6,nCols)])  +
                       3.0 * (mat[ELEM(8,25,nCols)] + mat[ELEM(25,8,nCols)])  +
                       3.0 * (mat[ELEM(9,24,nCols)] + mat[ELEM(24,9,nCols)])  +
                              mat[ELEM(11,19,nCols)] + mat[ELEM(19,11,nCols)] +
                       3.0 * (mat[ELEM(12,18,nCols)] + mat[ELEM(18,12,nCols)])+
                       6.0 * (mat[ELEM(14,15,nCols)] + mat[ELEM(15,14,nCols)]);

            /* 1 1 3 3 3 3 */
            ep[70] +=  2.0 * (mat[ELEM(1,49,nCols)] + mat[ELEM(49,1,nCols)])  +
                       4.0 * (mat[ELEM(3,44,nCols)] + mat[ELEM(44,3,nCols)])  +
                              mat[ELEM(4,34,nCols)] + mat[ELEM(34,4,nCols)]   +
                       8.0 * (mat[ELEM(6,29,nCols)] + mat[ELEM(29,6,nCols)])  +
                       6.0 * (mat[ELEM(9,25,nCols)] + mat[ELEM(25,9,nCols)])  +
                       4.0 * (mat[ELEM(12,19,nCols)] + mat[ELEM(19,12,nCols)])+
                      12.0 *  mat[ELEM(15,15,nCols)];

            /* 1 2 2 2 2 2 */
            ep[71] +=         mat[ELEM(1,50,nCols)] + mat[ELEM(50,1,nCols)]   +
                       5.0 * (mat[ELEM(2,45,nCols)] + mat[ELEM(45,2,nCols)])  +
                       5.0 * (mat[ELEM(5,30,nCols)] + mat[ELEM(30,5,nCols)])  +
                      10.0 * (mat[ELEM(7,26,nCols)] + mat[ELEM(26,7,nCols)])  +
                      10.0 * (mat[ELEM(13,16,nCols)] + mat[ELEM(16,13,nCols)]);

            /* 1 2 2 2 2 3 */
            ep[72] +=         mat[ELEM(1,51,nCols)] + mat[ELEM(51,1,nCols)]   +
                       4.0 * (mat[ELEM(2,46,nCols)] + mat[ELEM(46,2,nCols)])  +
                              mat[ELEM(3,45,nCols)] + mat[ELEM(45,3,nCols)]   +
                       4.0 * (mat[ELEM(5,31,nCols)] + mat[ELEM(31,5,nCols)])  +
                              mat[ELEM(6,30,nCols)] + mat[ELEM(30,6,nCols)]   +
                       6.0 * (mat[ELEM(7,27,nCols)] + mat[ELEM(27,7,nCols)])  +
                       4.0 * (mat[ELEM(8,26,nCols)] + mat[ELEM(26,8,nCols)])  +
                       6.0 * (mat[ELEM(13,17,nCols)] + mat[ELEM(17,13,nCols)])+
                       4.0 * (mat[ELEM(14,16,nCols)] + mat[ELEM(16,14,nCols)]);

            /* 1 2 2 2 3 3 */
            ep[73] +=         mat[ELEM(1,52,nCols)] + mat[ELEM(52,1,nCols)]   +
                       3.0 * (mat[ELEM(2,47,nCols)] + mat[ELEM(47,2,nCols)])  +
                       2.0 * (mat[ELEM(3,46,nCols)] + mat[ELEM(46,3,nCols)])  +
                       3.0 * (mat[ELEM(5,32,nCols)] + mat[ELEM(32,5,nCols)])  +
                       2.0 * (mat[ELEM(6,31,nCols)] + mat[ELEM(31,6,nCols)])  +
                       3.0 * (mat[ELEM(7,28,nCols)] + mat[ELEM(28,7,nCols)])  +
                       6.0 * (mat[ELEM(8,27,nCols)] + mat[ELEM(27,8,nCols)])  +
                              mat[ELEM(9,26,nCols)] + mat[ELEM(26,9,nCols)]   +
                       3.0 * (mat[ELEM(13,18,nCols)] + mat[ELEM(18,13,nCols)])+
                       6.0 * (mat[ELEM(14,17,nCols)] + mat[ELEM(17,14,nCols)])+
                              mat[ELEM(15,16,nCols)] + mat[ELEM(16,15,nCols)];

            /* 1 2 2 3 3 3 */
            ep[74] +=         mat[ELEM(1,53,nCols)] + mat[ELEM(53,1,nCols)]   +
                       2.0 * (mat[ELEM(2,48,nCols)] + mat[ELEM(48,2,nCols)])  +
                       3.0 * (mat[ELEM(3,47,nCols)] + mat[ELEM(47,3,nCols)])  +
                       2.0 * (mat[ELEM(5,33,nCols)] + mat[ELEM(33,5,nCols)])  +
                       3.0 * (mat[ELEM(6,32,nCols)] + mat[ELEM(32,6,nCols)])  +
                              mat[ELEM(7,29,nCols)] + mat[ELEM(29,7,nCols)]   +
                       6.0 * (mat[ELEM(8,28,nCols)] + mat[ELEM(28,8,nCols)])  +
                       3.0 * (mat[ELEM(9,27,nCols)] + mat[ELEM(27,9,nCols)])  +
                              mat[ELEM(13,19,nCols)] + mat[ELEM(19,13,nCols)] +
                       6.0 * (mat[ELEM(14,18,nCols)] + mat[ELEM(18,14,nCols)])+
                       3.0 * (mat[ELEM(15,17,nCols)] + mat[ELEM(17,15,nCols)]);

            /* 1 2 3 3 3 3 */
            ep[75] +=         mat[ELEM(1,54,nCols)] + mat[ELEM(54,1,nCols)]   +
                              mat[ELEM(2,49,nCols)] + mat[ELEM(49,2,nCols)]   +
                       4.0 * (mat[ELEM(3,48,nCols)] + mat[ELEM(48,3,nCols)])  +
                              mat[ELEM(5,34,nCols)] + mat[ELEM(34,5,nCols)]   +
                       4.0 * (mat[ELEM(6,33,nCols)] + mat[ELEM(33,6,nCols)])  +
                       4.0 * (mat[ELEM(8,29,nCols)] + mat[ELEM(29,8,nCols)])  +
                       6.0 * (mat[ELEM(9,28,nCols)] + mat[ELEM(28,9,nCols)])  +
                       4.0 * (mat[ELEM(14,19,nCols)] + mat[ELEM(19,14,nCols)])+
                       6.0 * (mat[ELEM(15,18,nCols)] + mat[ELEM(18,15,nCols)]);

            /* 1 3 3 3 3 3 */
            ep[76] +=         mat[ELEM(1,55,nCols)] + mat[ELEM(55,1,nCols)]   +
                       5.0 * (mat[ELEM(3,49,nCols)] + mat[ELEM(49,3,nCols)])  +
                       5.0 * (mat[ELEM(6,34,nCols)] + mat[ELEM(34,6,nCols)])  +
                      10.0 * (mat[ELEM(9,29,nCols)] + mat[ELEM(29,9,nCols)])  +
                      10.0 * (mat[ELEM(15,19,nCols)] + mat[ELEM(19,15,nCols)]);

            /* 2 2 2 2 2 2 */
            ep[77] +=  6.0 * (mat[ELEM(2,50,nCols)] + mat[ELEM(50,2,nCols)])  +
                      15.0 * (mat[ELEM(7,30,nCols)] + mat[ELEM(30,7,nCols)])  +
                      20.0 *  mat[ELEM(16,16,nCols)];

            /* 2 2 2 2 2 3 */
            ep[78] +=  5.0 * (mat[ELEM(2,51,nCols)] + mat[ELEM(51,2,nCols)])  +
                              mat[ELEM(3,50,nCols)] + mat[ELEM(50,3,nCols)]   +
                      10.0 * (mat[ELEM(7,31,nCols)] + mat[ELEM(31,7,nCols)])  +
                       5.0 * (mat[ELEM(8,30,nCols)] + mat[ELEM(30,8,nCols)])  +
                      10.0 * (mat[ELEM(16,17,nCols)] + mat[ELEM(17,16,nCols)]);

            /* 2 2 2 2 3 3 */
            ep[79] +=  4.0 * (mat[ELEM(2,52,nCols)] + mat[ELEM(52,2,nCols)])  +
                       2.0 * (mat[ELEM(3,51,nCols)] + mat[ELEM(51,3,nCols)])  +
                       6.0 * (mat[ELEM(7,32,nCols)] + mat[ELEM(32,7,nCols)])  +
                       8.0 * (mat[ELEM(8,31,nCols)] + mat[ELEM(31,8,nCols)])  +
                              mat[ELEM(9,30,nCols)] + mat[ELEM(30,9,nCols)]   +
                       4.0 * (mat[ELEM(16,18,nCols)] + mat[ELEM(18,16,nCols)])+
                      12.0 *  mat[ELEM(17,17,nCols)];

            /* 2 2 2 3 3 3 */
            ep[80] +=  3.0 * (mat[ELEM(2,53,nCols)] + mat[ELEM(53,2,nCols)])  +
                       3.0 * (mat[ELEM(3,52,nCols)] + mat[ELEM(52,3,nCols)])  +
                       3.0 * (mat[ELEM(7,33,nCols)] + mat[ELEM(33,7,nCols)])  +
                       9.0 * (mat[ELEM(8,32,nCols)] + mat[ELEM(32,8,nCols)])  +
                       3.0 * (mat[ELEM(9,31,nCols)] + mat[ELEM(31,9,nCols)])  +
                              mat[ELEM(16,19,nCols)] + mat[ELEM(19,16,nCols)] +
                       9.0 * (mat[ELEM(17,18,nCols)] + mat[ELEM(18,17,nCols)]);

            /* 2 2 3 3 3 3 */
            ep[81] +=  2.0 * (mat[ELEM(2,54,nCols)] + mat[ELEM(54,2,nCols)])  +
                       4.0 * (mat[ELEM(3,53,nCols)] + mat[ELEM(53,3,nCols)])  +
                              mat[ELEM(7,34,nCols)] + mat[ELEM(34,7,nCols)]   +
                       8.0 * (mat[ELEM(8,33,nCols)] + mat[ELEM(33,8,nCols)])  +
                       6.0 * (mat[ELEM(9,32,nCols)] + mat[ELEM(32,9,nCols)])  +
                       4.0 * (mat[ELEM(17,19,nCols)] + mat[ELEM(19,17,nCols)])+
                      12.0 *  mat[ELEM(18,18,nCols)];

            /* 2 3 3 3 3 3 */
            ep[82] +=        (mat[ELEM(2,55,nCols)] + mat[ELEM(55,2,nCols)])  +
                       5.0 * (mat[ELEM(3,54,nCols)] + mat[ELEM(54,3,nCols)])  +
                       5.0 * (mat[ELEM(8,34,nCols)] + mat[ELEM(34,8,nCols)])  +
                      10.0 * (mat[ELEM(9,33,nCols)] + mat[ELEM(33,9,nCols)])  +
                      10.0 * (mat[ELEM(18,19,nCols)] + mat[ELEM(19,18,nCols)]);

            /* 3 3 3 3 3 3 */
            ep[83] +=  6.0 * (mat[ELEM(3,55,nCols)] + mat[ELEM(55,3,nCols)])  +
                      15.0 * (mat[ELEM(9,34,nCols)] + mat[ELEM(34,9,nCols)])  +
                      20.0 *  mat[ELEM(19,19,nCols)];
        }
         
        return;
}
#endif // ESHELBY
