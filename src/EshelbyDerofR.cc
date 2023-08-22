
#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"
#include "Util.h"

/****************************************************************************
 *
 *      Function:     EshelbyDerofR
 *      Description:  
 *
 *      Arguments:
 *          position  coordinates of the target position
 *          origin    coordinates of the source position
 *          order     sum of taylor and multipole expansion orders
 *          drList    array in which to return the calculated coefficients
 *                    <drList> must be at least sequence[<order>+3]-4 elements
 *                    in length.
 *
 ***************************************************************************/
void EshelbyDerofR(real8 x[3], int order,real8 *drList)
            
{
        int   i, xvSize;
        int   index1, index2;
        real8 d, dinv2;
        real8 dinv[9];
        real8 *xv;

        for (i = 0; i < 9; i++) {
            dinv[i] = 0.0;
        }

        xvSize = sequence[order+3];
        xv     = (real8 *)calloc(1, xvSize * sizeof(real8));
         
        d = sqrt(DotProduct(x, x));

        dinv[0] = 1.0 / d;
        dinv2 = dinv[0] * dinv[0];

        dinv[1] = dinv[0] * dinv2;
        dinv[2] = dinv[1] * dinv2 * 3.0;

        xv[0] = 1.0;
        xv[1] = x[0];
        xv[2] = x[1];
        xv[3] = x[2];
        xv[4] = x[0] * x[0];
        xv[5] = x[0] * x[1];
        xv[6] = x[0] * x[2];
        xv[7] = x[1] * x[1];
        xv[8] = x[1] * x[2];
        xv[9] = x[2] * x[2];
        
        drList[0] = dinv[1] - xv[4] * dinv[2];  /* 1 1 */
        drList[1] =         - xv[5] * dinv[2];  /* 1 2 */
        drList[2] =         - xv[6] * dinv[2];  /* 1 3 */
        drList[3] = dinv[1] - xv[7] * dinv[2];  /* 2 2 */
        drList[4] =         - xv[8] * dinv[2];  /* 2 3 */
        drList[5] = dinv[1] - xv[9] * dinv[2];  /* 3 3 */

        if (order > 0) {

            dinv[3] = dinv[2] * dinv2 * 5.0;

            xv[10] = xv[4] * x[0];
            xv[11] = xv[5] * x[0];
            xv[12] = xv[6] * x[0];
            xv[13] = xv[7] * x[0];
            xv[14] = xv[8] * x[0];
            xv[15] = xv[9] * x[0];
            xv[16] = xv[7] * x[1];
            xv[17] = xv[8] * x[1];
            xv[18] = xv[9] * x[1];
            xv[19] = xv[9] * x[2];

            drList[ 6] = -3.0 * x[0] * dinv[2] + xv[10] * dinv[3]; /* 1 1 1 */
            drList[ 7] =       -x[1] * dinv[2] + xv[11] * dinv[3]; /* 1 1 2 */
            drList[ 8] =       -x[2] * dinv[2] + xv[12] * dinv[3]; /* 1 1 3 */
            drList[ 9] =       -x[0] * dinv[2] + xv[13] * dinv[3]; /* 1 2 2 */
            drList[10] =                         xv[14] * dinv[3]; /* 1 2 3 */
            drList[11] =       -x[0] * dinv[2] + xv[15] * dinv[3]; /* 1 3 3 */
            drList[12] = -3.0 * x[1] * dinv[2] + xv[16] * dinv[3]; /* 2 2 2 */
            drList[13] =       -x[2] * dinv[2] + xv[17] * dinv[3]; /* 2 2 3 */
            drList[14] =       -x[1] * dinv[2] + xv[18] * dinv[3]; /* 2 3 3 */
            drList[15] = -3.0 * x[2] * dinv[2] + xv[19] * dinv[3]; /* 3 3 3 */
        }


        if (order > 1) {

            dinv[4] = dinv[3] * dinv2 * 7.0;

            xv[20] = xv[10] * x[0];
            xv[21] = xv[11] * x[0];
            xv[22] = xv[12] * x[0];
            xv[23] = xv[13] * x[0];
            xv[24] = xv[14] * x[0];
            xv[25] = xv[15] * x[0];
            xv[26] = xv[16] * x[0];
            xv[27] = xv[17] * x[0];
            xv[28] = xv[18] * x[0];
            xv[29] = xv[19] * x[0];
            xv[30] = xv[16] * x[1];
            xv[31] = xv[17] * x[1];
            xv[32] = xv[18] * x[1];
            xv[33] = xv[19] * x[1];
            xv[34] = xv[19] * x[2];

            drList[16] =  dinv[2] * -3.0 +
                          dinv[3] *  6.0 * xv[4] -
                          dinv[4] * xv[20]; /* 1 1 1 1 */

            drList[17] =  dinv[3] * 3.0 * xv[5] -
                          dinv[4] * xv[21]; /* 1 1 1 2 */

            drList[18] =  dinv[3] * 3.0 * xv[6] -
                          dinv[4] * xv[22]; /* 1 1 1 3 */

            drList[19] = -dinv[2] +
                          dinv[3] * (xv[4]+xv[7]) -
                          dinv[4] *  xv[23]; /* 1 1 2 2 */

            drList[20] =  dinv[3] * xv[8] -
                          dinv[4] * xv[24]; /* 1 1 2 3 */

            drList[21] = -dinv[2] +
                          dinv[3] * (xv[4]+xv[9]) -
                          dinv[4] * xv[25]; /* 1 1 3 3 */

            drList[22] =  dinv[3] * 3.0 * xv[5] -
                          dinv[4] * xv[26]; /* 1 2 2 2 */

            drList[23] =  dinv[3] * xv[6] -
                          dinv[4] * xv[27]; /* 1 2 2 3 */

            drList[24] =  dinv[3] * xv[5] -
                          dinv[4] * xv[28]; /* 1 2 3 3 */

            drList[25] =  dinv[3] * 3.0 * xv[6] -
                          dinv[4] * xv[29]; /* 1 3 3 3 */

            drList[26] =  dinv[2] * -3.0 +
                          dinv[3] * 6.0 * xv[7] -
                          dinv[4] * xv[30]; /* 2 2 2 2 */

            drList[27] =  dinv[3] * 3.0 * xv[8] -
                          dinv[4] * xv[31]; /* 2 2 2 3 */

            drList[28] = -dinv[2] +
                          dinv[3] * (xv[7]+xv[9]) -
                          dinv[4] *  xv[32]; /* 2 2 3 3 */

            drList[29] =  dinv[3] * 3.0 * xv[8] -
                          dinv[4] * xv[33]; /* 2 3 3 3 */

            drList[30] =  dinv[2] * -3.0 +
                          dinv[3] *  6.0 * xv[9] -
                          dinv[4] *  xv[34]; /* 3 3 3 3 */
        }

        if (order > 2) {
 
            dinv[5] = dinv[4] * dinv2 * 9.0;

            for (index1 = 35, index2 = 20; index1 < 50; index1++, index2++) {
                xv[index1] = xv[index2] * x[0];
            }

            xv[50] = xv[30] * x[1];
            xv[51] = xv[31] * x[1];
            xv[52] = xv[32] * x[1];
            xv[53] = xv[33] * x[1];
            xv[54] = xv[34] * x[1];
            xv[55] = xv[34] * x[2];

            drList[31] = dinv[3] * 15.0 * x[0] -
                         dinv[4] * 10.0 * xv[10] +
                         dinv[5] * xv[35]; /* 1 1 1 1 1 */

            drList[32] = dinv[3] * 3.0 * x[1] -
                         dinv[4] * 6.0 * xv[11] +
                         dinv[5] * xv[36]; /* 1 1 1 1 2 */

            drList[33] = dinv[3] * 3.0 * x[2] -
                         dinv[4] * 6.0 * xv[12] +
                         dinv[5] * xv[37]; /* 1 1 1 1 3 */

            drList[34] = dinv[3] *  3.0 * x[0] -
                         dinv[4] * (xv[10] + 3.0 * xv[13]) +
                         dinv[5] *  xv[38]; /* 1 1 1 2 2 */

            drList[35] = dinv[4] * -3.0 * xv[14] +
                         dinv[5] * xv[39]; /* 1 1 1 2 3 */

            drList[36] = dinv[3] *  3.0 * x[0] -
                         dinv[4] * (xv[10] + 3.0 * xv[15]) +
                         dinv[5] *  xv[40]; /* 1 1 1 3 3 */

            drList[37] = dinv[3] *  3.0 * x[1] -
                         dinv[4] * (3.0 * xv[11] + xv[16]) +
                         dinv[5] * xv[41]; /* 1 1 2 2 2 */

            drList[38] = dinv[3] *  x[2] -
                         dinv[4] * (xv[12] + xv[17]) +
                         dinv[5] *  xv[42]; /* 1 1 2 2 3 */

            drList[39] = dinv[3] *  x[1] -
                         dinv[4] * (xv[11] + xv[18]) +
                         dinv[5] *  xv[43]; /* 1 1 2 3 3 */

            drList[40] = dinv[3] *  3.0 * x[2] -
                         dinv[4] * (3.0 * xv[12] + xv[19]) +
                         dinv[5] * xv[44]; /* 1 1 3 3 3 */

            drList[41] = dinv[3] *  3.0 * x[0] -
                         dinv[4] *  6.0 * xv[13] +
                         dinv[5] * xv[45]; /* 1 2 2 2 2 */

            drList[42] = dinv[4] * -3.0 * xv[14] +
                         dinv[5] * xv[46]; /* 1 2 2 2 3 */

            drList[43] = dinv[3] * x[0] -
                         dinv[4] * (xv[13] + xv[15]) +
                         dinv[5] * xv[47]; /* 1 2 2 3 3 */

            drList[44] = dinv[4] * -3.0 * xv[14] +
                         dinv[5] * xv[48]; /* 1 2 3 3 3 */

            drList[45] = dinv[3] * 3.0 * x[0] -
                         dinv[4] * 6.0 * xv[15] +
                         dinv[5] * xv[49]; /* 1 3 3 3 3 */

            drList[46] = dinv[3] * 15.0 * x[1] -
                         dinv[4] * 10.0 * xv[16] +
                         dinv[5] * xv[50]; /* 2 2 2 2 2 */

            drList[47] = dinv[3] * 3.0 * x[2] -
                         dinv[4] * 6.0 * xv[17] +
                         dinv[5] * xv[51]; /* 2 2 2 2 3 */

            drList[48] = dinv[3] *  3.0 * x[1] -
                         dinv[4] * (xv[16] + 3.0 * xv[18]) +
                         dinv[5] *  xv[52]; /* 2 2 2 3 3 */

            drList[49] = dinv[3] *  3.0 * x[2] -
                         dinv[4] * (xv[19] + 3.0 * xv[17]) +
                         dinv[5] *  xv[53]; /* 2 2 3 3 3 */

            drList[50] = dinv[3] * 3.0 * x[1] -
                         dinv[4] * 6.0 * xv[18] +
                         dinv[5] * xv[54]; /* 2 3 3 3 3 */

            drList[51] = dinv[3] * 15.0 * x[2] -
                         dinv[4] * 10.0 * xv[19] +
                         dinv[5] * xv[55]; /* 3 3 3 3 3 */
        }


        if (order > 3) {

            dinv[6] = dinv[5] * dinv2 * 11.0;

            for (index1 = 56, index2 = 35; index1 < 77; index1++, index2++) {
                xv[index1] = xv[index2] * x[0];
            }

            xv[77] = xv[50] * x[1];
            xv[78] = xv[51] * x[1];
            xv[79] = xv[52] * x[1];
            xv[80] = xv[53] * x[1];
            xv[81] = xv[54] * x[1];
            xv[82] = xv[55] * x[1];
            xv[83] = xv[55] * x[2];

            drList[52] = dinv[3] * 15.0 -
                         dinv[4] * 45.0 * xv[4] +
                         dinv[5] * 15.0 * xv[20] -
                         dinv[6] * xv[56]; /* 1 1 1 1 1 1 */

            drList[53] = dinv[4] * -15.0 * xv[5] +
                         dinv[5] *  10.0 * xv[21] -
                         dinv[6] *  xv[57]; /* 1 1 1 1 1 2 */

            drList[54] = dinv[4] * -15.0 * xv[6] +
                         dinv[5] *  10.0 * xv[22] -
                         dinv[6] *  xv[58]; /* 1 1 1 1 1 3 */

            drList[55] = dinv[3] *  3.0 -
                         dinv[4] *  3.0 * (2.0 * xv[4] + xv[7]) +
                         dinv[5] * (xv[20] + 6.0 * xv[23]) -
                         dinv[6] *  xv[59]; /* 1 1 1 1 2 2 */

            drList[56] = dinv[4] * -3.0 * xv[8] +
                         dinv[5] *  6.0 * xv[24] -
                         dinv[6] * xv[60]; /* 1 1 1 1 2 3 */
 
            drList[57] = dinv[3] *   3.0 -
                         dinv[4] *   3.0 * (2.0 * xv[4] + xv[9]) +
                         dinv[5] * (xv[20] + 6.0 * xv[25]) -
                         dinv[6] *  xv[61]; /* 1 1 1 1 3 3 */

            drList[58] = dinv[4] * -9.0 * xv[5] +
                         dinv[5] *  3.0 * (xv[21] + xv[26]) -
                         dinv[6] * xv[62]; /* 1 1 1 2 2 2 */

            drList[59] = dinv[4] *  -3.0 * xv[6] +
                         dinv[5] * (xv[22] + 3.0 * xv[27]) -
                         dinv[6] *  xv[63]; /* 1 1 1 2 2 3 */

            drList[60] = dinv[4] *  -3.0 * xv[5] +
                         dinv[5] * (xv[21] + 3.0 * xv[28]) -
                         dinv[6] *  xv[64]; /* 1 1 1 2 3 3 */

            drList[61] = dinv[4] * -9.0 * xv[6] +
                         dinv[5] *  3.0 * (xv[22] + xv[29]) -
                         dinv[6] * xv[65]; /* 1 1 1 3 3 3 */

            drList[62] = dinv[3] *  3.0 -
                         dinv[4] *  3.0 * (xv[4] + 2.0 * xv[7]) +
                         dinv[5] * (6.0 *  xv[23] + xv[30]) -
                         dinv[6] * xv[66]; /* 1 1 2 2 2 2 */

            drList[63] = dinv[4] * -3.0 * xv[8] +
                         dinv[5] * (3.0 * xv[24] + xv[31]) -
                         dinv[6] * xv[67]; /* 1 1 2 2 2 3 */

            drList[64] = dinv[3] -
                         dinv[4] * (xv[4] + xv[7] + xv[9]) +
                         dinv[5] * (xv[23] + xv[25] + xv[32]) -
                         dinv[6] *  xv[68]; /* 1 1 2 2 3 3 */

            drList[65] = dinv[4] * -3.0 * xv[8] +
                         dinv[5] * (3.0 * xv[24] + xv[33]) -
                         dinv[6] * xv[69]; /* 1 1 2 3 3 3 */

            drList[66] = dinv[3] *  3.0 -
                         dinv[4] *  3.0 * (xv[4] + 2.0 * xv[9]) +
                         dinv[5] * (6.0 * xv[25] + xv[34]) -
                         dinv[6] * xv[70]; /* 1 1 3 3 3 3 */

            drList[67] = dinv[4] * -15.0 * xv[5] +
                         dinv[5] *  10.0 * xv[26] -
                         dinv[6] *  xv[71]; /* 1 2 2 2 2 2 */

            drList[68] = dinv[4] * -3.0 * xv[6] +
                         dinv[5] *  6.0 * xv[27] -
                         dinv[6] * xv[72]; /* 1 2 2 2 2 3 */

            drList[69] = dinv[4] *  -3.0 * xv[5] +
                         dinv[5] * (xv[26] + 3.0 * xv[28]) -
                         dinv[6] *  xv[73]; /* 1 2 2 2 3 3 */

            drList[70] = dinv[4] * -3.0 * xv[6] +
                         dinv[5] * (3.0 * xv[27] + xv[29]) -
                         dinv[6] * xv[74]; /* 1 2 2 3 3 3 */

            drList[71] = dinv[4] * -3.0 * xv[5] +
                         dinv[5] *  6.0 * xv[28] -
                         dinv[6] * xv[75]; /* 1 2 3 3 3 3 */

            drList[72] = dinv[4] * -15.0 * xv[6] +
                         dinv[5] *  10.0 * xv[29] -
                         dinv[6] *  xv[76]; /* 1 3 3 3 3 3 */

            drList[73] = dinv[3] * 15.0 -
                         dinv[4] * 45.0 * xv[7] +
                         dinv[5] * 15.0 * xv[30] -
                         dinv[6] * xv[77]; /* 2 2 2 2 2 2 */

            drList[74] = dinv[4] * -15.0 * xv[8] +
                         dinv[5] *  10.0 * xv[31] -
                         dinv[6] *  xv[78]; /* 2 2 2 2 2 3 */

            drList[75] = dinv[3] *   3.0 -
                         dinv[4] *   3.0 * (2.0 * xv[7] + xv[9]) +
                         dinv[5] * (xv[30] + 6.0 *  xv[32]) -
                         dinv[6] *  xv[79]; /* 2 2 2 2 3 3 */

            drList[76] = dinv[4] * -9.0 * xv[8] +
                         dinv[5] *  3.0 * (xv[31] + xv[33]) -
                         dinv[6] * xv[80]; /* 2 2 2 3 3 3 */

            drList[77] = dinv[3] *  3.0 -
                         dinv[4] *  3.0 * (xv[7] + 2.0 * xv[9]) +
                         dinv[5] * (6.0 * xv[32] + xv[34]) -
                         dinv[6] * xv[81]; /* 2 2 3 3 3 3 */

            drList[78] = dinv[4] * -15.0 * xv[8] +
                         dinv[5] *  10.0 * xv[33] -
                         dinv[6] *  xv[82]; /* 2 3 3 3 3 3 */

            drList[79] = dinv[3] * 15.0 -
                         dinv[4] * 45.0 * xv[9] +
                         dinv[5] * 15.0 * xv[34] -
                         dinv[6] * xv[83]; /* 3 3 3 3 3 3 */
        }


        if (order > 4) {

            dinv[7] = dinv[6] * dinv2 * 13.0;

            for (index1 = 84, index2 = 56; index1 < 112; index1++, index2++) {
                xv[index1] = xv[index2] * x[0];
            }

            xv[112] = xv[77] * x[1];
            xv[113] = xv[78] * x[1];
            xv[114] = xv[79] * x[1];
            xv[115] = xv[80] * x[1];
            xv[116] = xv[81] * x[1];
            xv[117] = xv[82] * x[1];
            xv[118] = xv[83] * x[1];
            xv[119] = xv[83] * x[2];

            drList[80] = dinv[4] * -105.0 * x[0] +
                         dinv[5] *  105.0 * xv[10] -
                         dinv[6] *   21.0 * xv[35] +
                         dinv[7] *   xv[84]; /* 1 1 1 1 1 1 1 */

            drList[81] = dinv[4] * -15 * x[1] +
                         dinv[5] *  45 * xv[11] -
                         dinv[6] *  15 * xv[36] +
                         dinv[7] *  xv[85]; /* 1 1 1 1 1 1 2 */

            drList[82] = dinv[4] * -15 * x[2] +
                         dinv[5] *  45 * xv[12] -
                         dinv[6] *  15 * xv[37] +
                         dinv[7] *  xv[86]; /* 1 1 1 1 1 1 3 */

            drList[83] = dinv[4] * -15 * x[0] +
                         dinv[5] *   5 * (2 * xv[10] + 3 * xv[13]) -
                         dinv[6] * (xv[35] + 10 * xv[38]) +
                         dinv[7] *  xv[87]; /* 1 1 1 1 1 2 2 */

            drList[84] = dinv[5] * 15 * xv[14] -
                         dinv[6] * 10 * xv[39] +
                         dinv[7] * xv[88]; /* 1 1 1 1 1 2 3 */

            drList[85] = dinv[4] * -15 * x[0] +
                         dinv[5] *   5 * (2 * xv[10] + 3 * xv[15]) -
                         dinv[6] * (xv[35] + 10 * xv[40]) +
                         dinv[7] *  xv[89]; /* 1 1 1 1 1 3 3 */

            drList[86] = dinv[4] * -9 * x[1] +
                         dinv[5] *  3 * (6 * xv[11] + xv[16]) -
                         dinv[6] *  3 * (xv[36] + 2 * xv[41]) +
                         dinv[7] * xv[90]; /* 1 1 1 1 2 2 2 */

            drList[87] = dinv[4] *  -3 * x[2] +
                         dinv[5] *   3 * (2 * xv[12] + xv[17]) -
                         dinv[6] * (xv[37] + 6 * xv[42]) +
                         dinv[7] *  xv[91]; /* 1 1 1 1 2 2 3  */

            drList[88] = dinv[4] *  -3 * x[1] +
                         dinv[5] *   3 * (2 * xv[11] + xv[18]) -
                         dinv[6] * (xv[36] + 6 * xv[43]) +
                         dinv[7] *  xv[92]; /* 1 1 1 1 2 3 3  */

            drList[89] = dinv[4] * -9 * x[2] +
                         dinv[5] *  3 * (6 * xv[12] + xv[19]) -
                         dinv[6] *  3 * (xv[37] + 2 * xv[44]) +
                         dinv[7] * xv[93]; /* 1 1 1 1 3 3 3  */

            drList[90] = dinv[4] * -9 * x[0] +
                         dinv[5] *  3 * (xv[10] + 6 * xv[13]) -
                         dinv[6] *  3 * (2 * xv[38] + xv[45]) +
                         dinv[7] * xv[94]; /* 1 1 1 2 2 2 2  */

            drList[91] = dinv[5] *  9 * xv[14] -
                         dinv[6] *  3 * (xv[39] + xv[46]) +
                         dinv[7] * xv[95]; /* 1 1 1 2 2 2 3  */

            drList[92] = dinv[4] *  -3 * x[0] +
                         dinv[5] * (xv[10] + 3 * (xv[13] + xv[15])) -
                         dinv[6] * (xv[38] + xv[40] + 3 * xv[47]) +
                         dinv[7] *  xv[96]; /* 1 1 1 2 2 3 3  */

            drList[93] = dinv[5] *  9 * xv[14] -
                         dinv[6] *  3 * (xv[39] + xv[48]) +
                         dinv[7] * xv[97]; /* 1 1 1 2 3 3 3  */

            drList[94] = dinv[4] * -9 * x[0] +
                         dinv[5] *  3 * (xv[10] + 6 * xv[15]) -
                         dinv[6] *  3 * (2 * xv[40] + xv[49]) +
                         dinv[7] * xv[98]; /* 1 1 1 3 3 3 3  */

            drList[95] = dinv[4] * -15 * x[1] +
                         dinv[5] *   5 * (3 * xv[11] + 2 * xv[16]) -
                         dinv[6] * (10 * xv[41] + xv[50]) +
                         dinv[7] *  xv[99]; /* 1 1 2 2 2 2 2  */

            drList[96] = dinv[4] * -3 * x[2] +
                         dinv[5] *  3 * (xv[12] + 2 * xv[17]) -
                         dinv[6] * (6 * xv[42] + xv[51]) +
                         dinv[7] * xv[100]; /* 1 1 2 2 2 2 3  */

            drList[97] = dinv[4] *  -3 * x[1] +
                         dinv[5] * (xv[16] + 3 * (xv[11] + xv[18])) -
                         dinv[6] * (xv[41] + 3 * xv[43] + xv[52]) +
                         dinv[7] *  xv[101]; /* 1 1 2 2 2 3 3  */

            drList[98] = dinv[4] * -3 * x[2] +
                         dinv[5] * (3 * (xv[12] + xv[17]) + xv[19]) -
                         dinv[6] * (3 * xv[42] + xv[44] + xv[53]) +
                         dinv[7] * xv[102]; /* 1 1 2 2 3 3 3  */

            drList[99] = dinv[4] * -3 * x[1] +
                         dinv[5] *  3 * (xv[11] + 2 * xv[18]) -
                         dinv[6] * (6 * xv[43] + xv[54]) +
                         dinv[7] * xv[103]; /* 1 1 2 3 3 3 3  */

            drList[100] = dinv[4] * -15 * x[2] +
                          dinv[5] *   5 * (3 * xv[12] + 2 * xv[19]) -
                          dinv[6] * (10 * xv[44] + xv[55]) +
                          dinv[7] *  xv[104]; /* 1 1 3 3 3 3 3  */

            drList[101] = dinv[4] * -15 * x[0] +
                          dinv[5] *  45 * xv[13] -
                          dinv[6] *  15 * xv[45] +
                          dinv[7] *  xv[105]; /* 1 2 2 2 2 2 2  */

            drList[102] = dinv[5] * 15 * xv[14] -
                          dinv[6] * 10 * xv[46] +
                          dinv[7] * xv[106]; /* 1 2 2 2 2 2 3  */

            drList[103] = dinv[4] *  -3 * x[0] +
                          dinv[5] *   3 * (2 * xv[13] + xv[15]) -
                          dinv[6] * (xv[45] + 6 * xv[47]) +
                          dinv[7] *  xv[107]; /* 1 2 2 2 2 3 3  */

            drList[104] = dinv[5] *  9 * xv[14] -
                          dinv[6] *  3 * (xv[46] + xv[48]) +
                          dinv[7] * xv[108]; /* 1 2 2 2 3 3 3  */

            drList[105] = dinv[4] * -3 * x[0] +
                          dinv[5] *  3 * (xv[13] + 2 * xv[15]) -
                          dinv[6] * (6 * xv[47] + xv[49]) +
                          dinv[7] * xv[109]; /* 1 2 2 3 3 3 3  */

            drList[106] = dinv[5] * 15 * xv[14] -
                          dinv[6] * 10 * xv[48] +
                          dinv[7] * xv[110]; /* 1 2 3 3 3 3 3  */

            drList[107] = dinv[4] * -15 * x[0] +
                          dinv[5] *  45 * xv[15] -
                          dinv[6] *  15 * xv[49] +
                          dinv[7] *  xv[111]; /* 1 3 3 3 3 3 3  */

            drList[108] = dinv[4] * -105 * x[1] +
                          dinv[5] *  105 * xv[16] -
                          dinv[6] *   21 * xv[50] +
                          dinv[7] *   xv[112]; /* 2 2 2 2 2 2 2  */

            drList[109] = dinv[4] * -15 * x[2] +
                          dinv[5] *  45 * xv[17] -
                          dinv[6] *  15 * xv[51] +
                          dinv[7] *  xv[113]; /* 2 2 2 2 2 2 3  */

            drList[110] = dinv[4] * -15 * x[1] +
                          dinv[5] *   5 * (2 * xv[16] + 3 * xv[18]) -
                          dinv[6] * (xv[50] + 10 * xv[52]) +
                          dinv[7] *  xv[114]; /* 2 2 2 2 2 3 3  */

            drList[111] = dinv[4] * -9 * x[2] +
                          dinv[5] *  3 * (6 * xv[17] + xv[19]) -
                          dinv[6] *  3 * (xv[51] + 2 * xv[53]) +
                          dinv[7] *  xv[115]; /* 2 2 2 2 3 3 3  */

            drList[112] = dinv[4] * -9 * x[1] +
                          dinv[5] *  3 * (xv[16] + 6 * xv[18]) -
                          dinv[6] *  3 * (2 * xv[52] + xv[54]) +
                          dinv[7] * xv[116]; /* 2 2 2 3 3 3 3  */

            drList[113] = dinv[4] * -15 * x[2] +
                          dinv[5] *   5 * (3 * xv[17] + 2 * xv[19]) -
                          dinv[6] * (10 * xv[53] + xv[55]) +
                          dinv[7] *  xv[117]; /* 2 2 3 3 3 3 3  */

            drList[114] = dinv[4] * -15 * x[1] +
                          dinv[5] *  45 * xv[18] -
                          dinv[6] *  15 * xv[54] +
                          dinv[7] *  xv[118]; /* 2 3 3 3 3 3 3  */

            drList[115] = dinv[4] * -105 * x[2] +
                          dinv[5] *  105 * xv[19] -
                          dinv[6] *   21 * xv[55] +
                          dinv[7] *   xv[119]; /* 3 3 3 3 3 3 3  */
        }

        if (order > 5) {

            dinv[8] = dinv[7] * dinv2 * 15.0;

            for (index1 = 120, index2 = 84; index1 < 156; index1++, index2++) {
                xv[index1] = xv[index2] * x[0];
            }

            xv[156] = xv[112] * x[1];
            xv[157] = xv[113] * x[1];
            xv[158] = xv[114] * x[1];
            xv[159] = xv[115] * x[1];
            xv[160] = xv[116] * x[1];
            xv[161] = xv[117] * x[1];
            xv[162] = xv[118] * x[1];
            xv[163] = xv[119] * x[1];
            xv[164] = xv[119] * x[2];

            drList[116] = dinv[4] * -105 +
                          dinv[5] *  420 * xv[4] -
                          dinv[6] *  210 * xv[20] +
                          dinv[7] *   28 * xv[56] -
                          dinv[8] *   xv[120]; /*  1 1 1 1 1 1 1 1  */

            drList[117] = dinv[5] * 105 * xv[5] -
                          dinv[6] * 105 * xv[21] +
                          dinv[7] *  21 * xv[57] -
                          dinv[8] *  xv[121]; /*  1 1 1 1 1 1 1 2  */

            drList[118] = dinv[5] * 105 * xv[6] -
                          dinv[6] * 105 * xv[22] +
                          dinv[7] *  21 * xv[58] -
                          dinv[8] *  xv[122]; /*  1 1 1 1 1 1 1 3  */

            drList[119] = dinv[4] * -15 +
                          dinv[5] *  15 * (3 * xv[4] + xv[7]) -
                          dinv[6] *  15 * (xv[20] + 3 * xv[23]) +
                          dinv[7] * (xv[56] + 15 * xv[59]) -
                          dinv[8] *  xv[123]; /*  1 1 1 1 1 1 2 2  */

            drList[120] = dinv[5] *  15 * xv[8] -
                          dinv[6] *  45 * xv[24] +
                          dinv[7] *  15 * xv[60] -
                          dinv[8] *  xv[124]; /*  1 1 1 1 1 1 2 3  */

            drList[121] = dinv[4] * -15 +
                          dinv[5] *  15 * (3 * xv[4] + xv[9]) -
                          dinv[6] *  15 * (xv[20] + 3 * xv[25]) +
                          dinv[7] * (xv[56] + 15 * xv[61]) -
                          dinv[8] *  xv[125]; /*  1 1 1 1 1 1 3 3  */

            drList[122] = dinv[5] * 45 * xv[5] -
                          dinv[6] * 15 * (2 * xv[21] + xv[26]) +
                          dinv[7] * (3 * xv[57] + 10 * xv[62]) -
                          dinv[8] * xv[126]; /*  1 1 1 1 1 2 2 2  */
 
            drList[123] = dinv[5] *  15 * xv[6] -
                          dinv[6] *   5 * (2 * xv[22] + 3 * xv[27]) +
                          dinv[7] * (xv[58] + 10 * xv[63]) -
                          dinv[8] *  xv[127]; /*  1 1 1 1 1 2 2 3  */

            drList[124] = dinv[5] *  15 * xv[5] -
                          dinv[6] *   5 * (2 * xv[21] + 3 * xv[28]) +
                          dinv[7] * (xv[57] + 10 * xv[64]) -
                          dinv[8] *  xv[128]; /*  1 1 1 1 1 2 3 3  */

            drList[125] = dinv[5] * 45 * xv[6] -
                          dinv[6] * 15 * (2 * xv[22] + xv[29]) +
                          dinv[7] * (3 * xv[58] + 10 * xv[65]) -
                          dinv[8] * xv[129]; /*  1 1 1 1 1 3 3 3  */

            drList[126] = dinv[4] *  -9 +
                          dinv[5] *  18 * (xv[4] + xv[7]) - 
                          dinv[6] *   3 * (xv[20] + 12 * xv[23] + xv[30]) +
                          dinv[7] *   6 * (xv[59] + xv[66]) -
                          dinv[8] *  xv[130]; /*  1 1 1 1 2 2 2 2  */

            drList[127] = dinv[5] *  9 * xv[8] -
                          dinv[6] *  3 * (6 * xv[24] + xv[31]) +
                          dinv[7] *  3 * (xv[60] + 2 * xv[67]) -
                          dinv[8] * xv[131]; /*  1 1 1 1 2 2 2 3  */

            drList[128] = dinv[4] *  -3 +
                          dinv[5] *   3 * (2 * xv[4] + xv[7] + xv[9]) -
                          dinv[6] * (xv[20] + 6*xv[23] + 6*xv[25] + 3*xv[32]) +
                          dinv[7] * (xv[59] + xv[61] + 6 * xv[68]) -
                          dinv[8] *  xv[132]; /*  1 1 1 1 2 2 3 3  */

            drList[129] = dinv[5] *  9 * xv[8] -
                          dinv[6] *  3 * (6 * xv[24] + xv[33]) +
                          dinv[7] *  3 * (xv[60] + 2 * xv[69]) -
                          dinv[8] * xv[133]; /*  1 1 1 1 2 3 3 3  */

            drList[130] = dinv[4] * -9 +
                          dinv[5] * 18 * (xv[4] + xv[9]) -
                          dinv[6] *  3 * (xv[20] + 12 * xv[25] + xv[34]) +
                          dinv[7] *  6 * (xv[61] + xv[70]) -
                          dinv[8] * xv[134]; /*  1 1 1 1 3 3 3 3  */

            drList[131] = dinv[5] *  45 * xv[5] -
                          dinv[6] *  15 * (xv[21] + 2 * xv[26]) +
                          dinv[7] * (10 * xv[62] + 3 * xv[71]) -
                          dinv[8] *  xv[135]; /*  1 1 1 2 2 2 2 2  */

            drList[132] = dinv[5] * 9 * xv[6] -
                          dinv[6] * 3 * (xv[22] + 6 * xv[27]) +
                          dinv[7] * 3 * (2 * xv[63] + xv[72]) -
                          dinv[8] * xv[136]; /*  1 1 1 2 2 2 2 3  */

            drList[133] = dinv[5] *  9 * xv[5] -
                          dinv[6] *  3 * (xv[21] + xv[26] + 3 * xv[28]) +
                          dinv[7] * (xv[62] + 3 * (xv[64] + xv[73])) -
                          dinv[8] *  xv[137]; /*  1 1 1 2 2 2 3 3  */

            drList[134] = dinv[5] *  9 * xv[6] -
                          dinv[6] *  3 * (xv[22] + 3 * xv[27] + xv[29]) +
                          dinv[7] * (xv[65] + 3 * (xv[63] + xv[74])) -
                          dinv[8] *  xv[138]; /*  1 1 1 2 2 3 3 3  */

            drList[135] = dinv[5] * 9 * xv[5] -
                          dinv[6] * 3 * (xv[21] + 6 * xv[28]) +
                          dinv[7] * 3 * (2 * xv[64] + xv[75]) -
                          dinv[8] * xv[139]; /*  1 1 1 2 3 3 3 3  */

            drList[136] = dinv[5] *  45 * xv[6] -
                          dinv[6] *  15 * (xv[22] + 2 * xv[29]) +
                          dinv[7] * (10 * xv[65] + 3 * xv[76]) -
                          dinv[8] *  xv[140]; /*  1 1 1 3 3 3 3 3  */
 
            drList[137] = dinv[4] * -15 +
                          dinv[5] *  15 * (xv[4] + 3 * xv[7]) -
                          dinv[6] *  15 * (3 * xv[23] + xv[30]) +
                          dinv[7] * (15 * xv[66] + xv[77]) -
                          dinv[8] *  xv[141]; /*  1 1 2 2 2 2 2 2  */

            drList[138] = dinv[5] *  15 * xv[8] -
                          dinv[6] *   5 * (3 * xv[24] + 2 * xv[31]) +
                          dinv[7] * (10 * xv[67] + xv[78]) -
                          dinv[8] *  xv[142]; /*  1 1 2 2 2 2 2 3  */

            drList[139] = dinv[4] * -3 +
                          dinv[5] *  3 * (xv[4] + 2 * xv[7] + xv[9]) -
                          dinv[6] * (6 * xv[23] + 3*xv[25] + xv[30] + 6*xv[32])+
                          dinv[7] * (xv[66] + 6 * xv[68] + xv[79]) -
                          dinv[8] *  xv[143]; /*  1 1 2 2 2 2 3 3  */

            drList[140] = dinv[5] *  9 * xv[8] -
                          dinv[6] *  3 * (3 * xv[24] + xv[31] + xv[33]) +
                          dinv[7] * (3 * (xv[67] + xv[69]) + xv[80]) -
                          dinv[8] *  xv[144]; /*  1 1 2 2 2 3 3 3  */

            drList[141] = dinv[4] * -3 +
                          dinv[5] *  3 * (xv[4] + xv[7] + 2 * xv[9]) -
                          dinv[6] * (3 * xv[23] + 6*xv[25] + 6*xv[32] + xv[34])+
                          dinv[7] * (6 * xv[68] + xv[70] + xv[81]) -
                          dinv[8] *  xv[145]; /*  1 1 2 2 3 3 3 3  */

            drList[142] = dinv[5] *  15 * xv[8] -
                          dinv[6] *   5 * (3 * xv[24] + 2 * xv[33]) +
                          dinv[7] * (10 * xv[69] + xv[82]) -
                          dinv[8] *  xv[146]; /*  1 1 2 3 3 3 3 3  */

            drList[143] = dinv[4] * -15 +
                          dinv[5] *  15 * (xv[4] + 3 * xv[9]) -
                          dinv[6] *  15 * (3 * xv[25] + xv[34]) +
                          dinv[7] * (15 * xv[70] + xv[83]) -
                          dinv[8] *  xv[147]; /*  1 1 3 3 3 3 3 3  */

            drList[144] = dinv[5] * 105 * xv[5] -
                          dinv[6] * 105 * xv[26] +
                          dinv[7] * 21 * xv[71] -
                          dinv[8] * xv[148]; /*  1 2 2 2 2 2 2 2  */

            drList[145] = dinv[5] * 15 * xv[6] -
                          dinv[6] *  45 * xv[27] +
                          dinv[7] *  15 * xv[72] -
                          dinv[8] *  xv[149]; /*  1 2 2 2 2 2 2 3  */

            drList[146] = dinv[5] * 15 * xv[5] -
                          dinv[6] *  5 * (2 * xv[26] + 3 * xv[28]) +
                          dinv[7] * (xv[71] + 10 * xv[73]) -
                          dinv[8] *  xv[150]; /*  1 2 2 2 2 2 3 3  */

            drList[147] = dinv[5] * 9 * xv[6] -
                          dinv[6] * 3 * (6 * xv[27] + xv[29]) +
                          dinv[7] * 3 * (xv[72] + 2 * xv[74]) -
                          dinv[8] * xv[151]; /*  1 2 2 2 2 3 3 3  */

            drList[148] = dinv[5] * 9 * xv[5] -
                          dinv[6] * 3 * (xv[26] + 6 * xv[28]) +
                          dinv[7] * 3 * (2 * xv[73] + xv[75]) -
                          dinv[8] * xv[152]; /*  1 2 2 2 3 3 3 3  */

            drList[149] = dinv[5] *  15 * xv[6] -
                          dinv[6] *   5 * (3 * xv[27] + 2 * xv[29]) +
                          dinv[7] * (10 * xv[74] + xv[76]) -
                          dinv[8] * xv[153]; /*  1 2 2 3 3 3 3 3  */

            drList[150] = dinv[5] * 15 * xv[5] -
                          dinv[6] * 45 * xv[28] +
                          dinv[7] * 15 * xv[75] -
                          dinv[8] * xv[154]; /*  1 2 3 3 3 3 3 3  */

            drList[151] = dinv[5] * 105 * xv[6] -
                          dinv[6] * 105 * xv[29] +
                          dinv[7] *  21 * xv[76] -
                          dinv[8] * xv[155]; /*  1 3 3 3 3 3 3 3  */

            drList[152] = dinv[4] * -105 +
                          dinv[5] *  420 * xv[7] -
                          dinv[6] *  210 * xv[30] +
                          dinv[7] *   28 * xv[77] -
                          dinv[8] * xv[156]; /*  2 2 2 2 2 2 2 2  */

            drList[153] = dinv[5] * 105 * xv[8] -
                          dinv[6] * 105 * xv[31] +
                          dinv[7] * 21 * xv[78] -
                          dinv[8] * xv[157]; /*  2 2 2 2 2 2 2 3  */

            drList[154] = dinv[4] * -15 +
                          dinv[5] *  15 * (3 * xv[7] + xv[9]) -
                          dinv[6] *  15 * (xv[30] + 3 * xv[32]) +
                          dinv[7] * (xv[77] + 15 * xv[79]) -
                          dinv[8] *  xv[158]; /*  2 2 2 2 2 2 3 3  */

            drList[155] = dinv[5] * 45 * xv[8] -
                          dinv[6] * 15 * (2 * xv[31] + xv[33]) +
                          dinv[7] * (3 * xv[78] + 10 * xv[80]) -
                          dinv[8] * xv[159]; /*  2 2 2 2 2 3 3 3  */

            drList[156] = dinv[4] * -9 +
                          dinv[5] * 18 * (xv[7] + xv[9]) -
                          dinv[6] *  3 * (xv[30] + 12 * xv[32] + xv[34]) +
                          dinv[7] *  6 * (xv[79] + xv[81]) -
                          dinv[8] * xv[160]; /*  2 2 2 2 3 3 3 3  */

            drList[157] = dinv[5] *  45 * xv[8] -
                          dinv[6] *  15 * (xv[31] + 2 * xv[33]) +
                          dinv[7] * (10 * xv[80] + 3 * xv[82]) -
                          dinv[8] * xv[161]; /*  2 2 2 3 3 3 3 3  */

            drList[158] = dinv[4] * -15 +
                          dinv[5] *  15 * (xv[7] + 3 * xv[9]) -
                          dinv[6] *  15 * (3 * xv[32] + xv[34]) +
                          dinv[7] * (15 * xv[81] + xv[83]) -
                          dinv[8] * xv[162]; /*  2 2 3 3 3 3 3 3  */
 
            drList[159] = dinv[5] * 105 * xv[8] -
                          dinv[6] * 105 * xv[33] +
                          dinv[7] *  21 * xv[82] -
                          dinv[8] * xv[163]; /* 2 3 3 3 3 3 3 3  */
 
            drList[160] = dinv[4] * -105 +
                          dinv[5] *  420 * xv[9] -
                          dinv[6] *  210 * xv[34] +
                          dinv[7] *   28 * xv[83] -
                          dinv[8] * xv[164]; /* 3 3 3 3 3 3 3 3*/
        }

        free(xv);

        return;
}
