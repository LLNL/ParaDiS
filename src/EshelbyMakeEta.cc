
#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"


/****************************************************************************
 *
 *      Function:    EshelbyMakeEta
 *      Description: Create a multipole expansion representing the stress
 *                   field from an Eshelby inclusion.
 *
 *      Arguments:
 *          Eva      Volumetric strain associated with the inclusion
 *          R0       radius of the inclusion in units of b
 *          position Coordinates of the center of the inclusion
 *          center   Coordinates of the multipole expansion center
 *          order    Multipole expansion order.  Must be in the
 *                   range 0 <= order <= 6
 *          etav     Array in which to return the multipole coefficients.
 *                   This array must be at least sequence[<order>+1]
 *                   elements in length.
 *
 ***************************************************************************/

/*
 *      NOTE: multipoleOrder must be in the range:
 *          0 <= multipoleOrder <= 6
 *
 *      etav   Array in which to return the eta coefficients to the
 *             caller.  <etav> must be at least sequence[<order>+1]
 *             elements in length.
 */
void EshelbyMakeEta(real8 Ev, real8 R0, real8 position[3], real8 center[3],
                    int order, real8 *etav)
{
        int    i, etaCnt;
        real8  scalar;
        real8  r[3];

        r[0] = -(position[0] - center[0]);
        r[1] = -(position[1] - center[1]);
        r[2] = -(position[2] - center[2]);

        scalar = Ev * R0 * R0 * R0;

        etaCnt = sequence[order+1];

/*
 *      IMPORTANT!  The series of 'if' clauses below are set up
 *      so that the operations performed for a given multipole order
 *      include all operations for multipole orders less than the
 *      current value.
 */
        if (order >= 0) {
            etav[0] = 1.0;
        }

        if (order >= 1) {
            etav[1] = r[0];
            etav[2] = r[1];
            etav[3] = r[2];
        }

        if (order >= 2) {
            etav[4] = r[0] * etav[1];
            etav[5] = r[0] * etav[2];
            etav[6] = r[0] * etav[3];
            etav[7] = r[1] * etav[2];
            etav[8] = r[1] * etav[3];
            etav[9] = r[2] * etav[3];
        }

        if (order >= 3) {
            etav[10] = r[0] * etav[4];
            etav[11] = r[0] * etav[5];
            etav[12] = r[0] * etav[6];
            etav[13] = r[0] * etav[7];
            etav[14] = r[0] * etav[8];
            etav[15] = r[0] * etav[9];
            etav[16] = r[1] * etav[7];
            etav[17] = r[1] * etav[8];
            etav[18] = r[1] * etav[9];
            etav[19] = r[2] * etav[9];
        }

        if (order >= 4) {
            etav[20] = r[0] * etav[10];
            etav[21] = r[0] * etav[11];
            etav[22] = r[0] * etav[12];
            etav[23] = r[0] * etav[13];
            etav[24] = r[0] * etav[14];
            etav[25] = r[0] * etav[15];
            etav[26] = r[0] * etav[16];
            etav[27] = r[0] * etav[17];
            etav[28] = r[0] * etav[18];
            etav[29] = r[0] * etav[19];
            etav[30] = r[1] * etav[16];
            etav[31] = r[1] * etav[17];
            etav[32] = r[1] * etav[18];
            etav[33] = r[1] * etav[19];
            etav[34] = r[2] * etav[19];
        }

        if (order >= 5) {
            etav[35] = r[0] * etav[20];
            etav[36] = r[0] * etav[21];
            etav[37] = r[0] * etav[22];
            etav[38] = r[0] * etav[23];
            etav[39] = r[0] * etav[24];
            etav[40] = r[0] * etav[25];
            etav[41] = r[0] * etav[26];
            etav[42] = r[0] * etav[27];
            etav[43] = r[0] * etav[28];
            etav[44] = r[0] * etav[29];
            etav[45] = r[0] * etav[30];
            etav[46] = r[0] * etav[31];
            etav[47] = r[0] * etav[32];
            etav[48] = r[0] * etav[33];
            etav[49] = r[0] * etav[34];
            etav[50] = r[1] * etav[30];
            etav[51] = r[1] * etav[31];
            etav[52] = r[1] * etav[32];
            etav[53] = r[1] * etav[33];
            etav[54] = r[1] * etav[34];
            etav[55] = r[2] * etav[34];
        }

        if (order >= 6) {
            etav[56] = r[0] * etav[35];
            etav[57] = r[0] * etav[36];
            etav[58] = r[0] * etav[37];
            etav[59] = r[0] * etav[38];
            etav[60] = r[0] * etav[39];
            etav[61] = r[0] * etav[40];
            etav[62] = r[0] * etav[41];
            etav[63] = r[0] * etav[42];
            etav[64] = r[0] * etav[43];
            etav[65] = r[0] * etav[44];
            etav[66] = r[0] * etav[45];
            etav[67] = r[0] * etav[46];
            etav[68] = r[0] * etav[47];
            etav[69] = r[0] * etav[48];
            etav[70] = r[0] * etav[49];
            etav[71] = r[0] * etav[50];
            etav[72] = r[0] * etav[51];
            etav[73] = r[0] * etav[52];
            etav[74] = r[0] * etav[53];
            etav[75] = r[0] * etav[54];
            etav[76] = r[0] * etav[55];
            etav[77] = r[1] * etav[50];
            etav[78] = r[1] * etav[51];
            etav[79] = r[1] * etav[52];
            etav[80] = r[1] * etav[53];
            etav[81] = r[1] * etav[54];
            etav[82] = r[1] * etav[55];
            etav[83] = r[2] * etav[55];
        }

        for (i = 0; i < etaCnt; i++) {
            etav[i] *= scalar;
        }

        return;
}
