#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif


/**************************************************************************
 *
 *      Function:     TAEdgeDrag
 *      Description:  This function returns the drag function g_edge
 *                    along the line, the glide directions and the climb
 *                    directions as well as its derivative Dg_edge/Dv.
 *                    EdgeDrag is linear.
 *               
 *      Arguments:
 *
 *          IN:  vel           Test velocity for current iteration of
 *                             Newton-Raphson loop specified in the crystal frame.
 * 
 *          IN:  burg          Burgers vector specified in the crystal frame.
 *
 *          IN:  linedir       Direction of the edge part of the line direction.
 *
 *          IN:  dragClimb     Drag coefficient for the climb direction in Pa. s
 *
 *          IN:  dragLine      Drag coefficient for the line direction in Pa. s
 *
 *          IN:  dragGlide     Drag coefficient for glide associated with the
 *                             climbDir glide plane in Pa. s
 *
 *          OUT: fEdgeDrag     Drag vector for edge component
 *
 *          OUT: dfEdgeDragdv  Drag matrix, derivative of fEdgeDrag for
 *                             edge component.
 *
 *************************************************************************/
void TAEdgeDrag(real8 vel[3], real8 burg[3], real8 lineDir[3],
                real8 dragClimb, real8 dragLine, real8 dragGlide,
                real8 fEdgeDrag[3], real8 dfEdgeDragdv[3][3])
{
        int   m, n;
        real8 climbDir[3], glideDir[3];
        real8 glideDirMat[3][3], climbDirMat[3][3], lineDirMat[3][3];


	// glide direction is the Burgers vector normalized.
        VECTOR_COPY(glideDir, burg);
        NormalizeVec(glideDir);

	NormalizedCrossVector(lineDir, glideDir,climbDir);

        Vec3TransposeAndMult(glideDir, glideDirMat);
        Vec3TransposeAndMult(climbDir, climbDirMat);
        Vec3TransposeAndMult( lineDir,  lineDirMat);
       
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfEdgeDragdv[m][n] = dragClimb * climbDirMat[m][n] +
                                     dragLine  *  lineDirMat[m][n] +
                                     dragGlide * glideDirMat[m][n];
            }
        }

        Matrix33Vector3Multiply(dfEdgeDragdv, vel, fEdgeDrag);

        return;
}

/**************************************************************************
 *
 *      Function:     TAJunctionDrag
 *
 *      Description:  This function returns the drag for junction dislocation
 *                    segments.  Motion of junctions is restricted to the
 *                    corresponding line direction. It copies the mobility 
 *                    for BCC_0 for junctions.
 *
 *      Arguments:
 *
 *          IN:  vel               Test velocity for current iteration of
 *                                 Newton-Raphson loop specified in the
 *                                 crystal frame.
 *
 *          IN:  lineDir           Dislocation line direction
 *
 *          IN:  BJunctionLine     Drag coefficient for line direction
 *
 *          IN:  BJunctionGlide    Drag coefficient for glide
 *
 *          OUT: fJunctionDrag     Drag vector for junctions
 *
 *          OUT: dfJunctionDragdv  Drag matrix, derivative of fJunction
 *                                 for junctions.
 *
 *************************************************************************/
void TAJunctionDrag(real8 vel[3], real8 lineDir[3],
                  real8 BJunctionLine, real8 BJunctionGlide,
                  real8 fJunctionDrag[3], real8 dfJunctionDragdv[3][3])
{
        int    m, n;
        real8  glideDirMat[3][3], lineDirMat[3][3];

        glideDirMat[0][0] = 1.0 - lineDir[0] * lineDir[0];
        glideDirMat[1][1] = 1.0 - lineDir[1] * lineDir[1];
        glideDirMat[2][2] = 1.0 - lineDir[2] * lineDir[2];

        glideDirMat[0][1] =     - lineDir[0] * lineDir[1];
        glideDirMat[0][2] =     - lineDir[0] * lineDir[2];

        glideDirMat[1][0] =     - lineDir[1] * lineDir[0];
        glideDirMat[1][2] =     - lineDir[1] * lineDir[2];

        glideDirMat[2][0] =     - lineDir[2] * lineDir[0];
        glideDirMat[2][1] =     - lineDir[2] * lineDir[1];
 
        Vec3TransposeAndMult(lineDir, lineDirMat);

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfJunctionDragdv[m][n] = BJunctionGlide * glideDirMat[m][n]+
                                         BJunctionLine  * lineDirMat[m][n];
            }
        }

        Matrix33Vector3Multiply(dfJunctionDragdv, vel, fJunctionDrag);

        return;
}
