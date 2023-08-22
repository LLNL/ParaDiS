/**************************************************************************
 *
 *      Module:       LinearMobilityFunctions
 *      Author:       Sylvie Aubry
 *      Description:  Contains functions in common for different linear
 *                    mobilities.
 *
 ***************************************************************************/


#include "mpi_portability.h"

#include "Home.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

/**************************************************************************
 *
 *      Function:     EdgeDrag
 *      Description:  This function returns the drag function g_edge
 *                    along the line, the glide directions and the climb
 *                    directions as well as its derivative Dg_edge/Dv.
 *                    EdgeDrag is linear.
 *               
 *      Arguments:
 *
 *          IN:  vel       Test velocity for current iteration of
 *                         Newton-Raphson loop specified in the crystal frame.
 * 
 *          IN:  burg      Burgers vector specified in the crystal frame.
 *
 *          IN:  edgeDir   Line direction for one of the two edge components
 *                         into which the dislocation segment was decomposed.
 *                         Specified in the crystal frame.
 *
 *          IN:  dragClimb Drag coefficient for the climb direction
 *
 *          IN:  dragLine  Drag coefficient for the line direction
 *
 *          IN:  dragGlide Drag coefficient for glide associated with the
 *                         climbDir glide plane
 *
 *          OUT: fEdgeDrag Drag vector for edge component
 *
 *          OUT: dfEdgeDragdv  Drag matrix, derivative of fEdgeDrag for
 *                         edge component.
 *
 *************************************************************************/
void EdgeDrag(real8 vel[3], real8 burg[3], real8 edgeDir[3], 
              real8 dragClimb, real8 dragLine, real8 dragGlide,
              real8 fEdgeDrag[3], real8 dfEdgeDragdv[3][3])
{
        int   m, n;
        real8 lineDir[3], glideDir[3], climbDir[3];
        real8 glideDirMat[3][3], climbDirMat[3][3], lineDirMat[3][3];

        VECTOR_COPY(glideDir, burg);
        NormalizeVec(glideDir);

        VECTOR_COPY(lineDir, edgeDir);
        NormalizeVec(lineDir);

        cross(lineDir, glideDir, climbDir);

        Vec3TransposeAndMult(glideDir, glideDirMat);
        Vec3TransposeAndMult(climbDir, climbDirMat);
        Vec3TransposeAndMult( lineDir,  lineDirMat);
       
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfEdgeDragdv[m][n] = dragClimb * climbDirMat[m][n] +
                                     dragLine  *  lineDirMat[m][n]  +
                                     dragGlide * glideDirMat[m][n];
            }
        }

        Matrix33Vector3Multiply(dfEdgeDragdv, vel, fEdgeDrag);

        return;
}

/**************************************************************************
 *
 *      Function:     ScrewDrag
 *      Description:  This function returns the drag function g_screw
 *                    along the line, glide and climb directions as well
 *                    as its derivative Dg_screw/Dv.  ScrewDrag is linear.
 *
 *      Arguments:
 *
 *          IN:  vel       Test velocity for current iteration of
 *                         Newton-Raphson loop specified in the crystal frame.
 * 
 *          IN:  burg      Burgers vector specified in the crystal frame.
 *
 *          IN:  dragClimb Drag coefficient for the climb direction
 *
 *          IN:  dragLine  Drag coefficient for the line direction
 *
 *          IN:  dragGlide Drag coefficient for glide associated with the
 *                         <climbDir> glide plane
 *
 *          OUT: fScrewDrag Drag vector for screw component
 *
 *          OUT: dfScrewDragdv  Drag matrix, derivative of fScrewDrag for
 *                         screw component.
 *
 *************************************************************************/
void ScrewDrag(real8 vel[3], real8 burg[3],
               real8 dragClimb, real8 dragLine, real8 dragGlide,
               real8 fScrewDrag[3], real8 dfScrewDragdv[3][3])
{
        int   m, n;
        real8 lineDir[3];
        real8 glideDirMat[3][3], lineDirMat[3][3];

        VECTOR_COPY(lineDir, burg);
        NormalizeVec(lineDir);

        glideDirMat[0][0] = 1.0 - lineDir[0] * lineDir[0];
        glideDirMat[1][1] = 1.0 - lineDir[1] * lineDir[1];
        glideDirMat[2][2] = 1.0 - lineDir[2] * lineDir[2];

        glideDirMat[0][1] =     - lineDir[0] * lineDir[1];
        glideDirMat[0][2] =     - lineDir[0] * lineDir[2];

        glideDirMat[1][0] =     - lineDir[1] * lineDir[0];
        glideDirMat[1][2] =     - lineDir[1] * lineDir[2];

        glideDirMat[2][0] =     - lineDir[2] * lineDir[0];
        glideDirMat[2][1] =     - lineDir[2] * lineDir[1];

        Vec3TransposeAndMult(lineDir,  lineDirMat);
        
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfScrewDragdv[m][n] = dragLine  * lineDirMat[m][n]  +
                                      dragGlide * glideDirMat[m][n];
            }
        }

        Matrix33Vector3Multiply(dfScrewDragdv, vel, fScrewDrag);

        return;
}


/**************************************************************************
 *
 *      Function:     ScrewDrag2
 *      Description:  This function returns the drag function g_screw
 *                    along the line, glide and climb directions as well
 *                    as its derivative Dg_screw/Dv.  ScrewDrag is linear.
 *
 *      Arguments:
 *
 *          IN:  vel       Test velocity for current iteration of
 *                         Newton-Raphson loop specified in the crystal frame.
 * 
 *          IN:  burg      Burgers vector specified in the crystal frame.
 *
 *          IN:  climbDir  Glide plane for the dislocation segment specified
 *                         in the crystal frame.
 *
 *          IN:  dragClimb Drag coefficient for the climb direction
 *
 *          IN:  dragLine  Drag coefficient for the line direction
 *
 *          IN:  dragGlide Drag coefficient for glide associated with the
 *                         <climbDir> glide plane
 *
 *          OUT: fScrewDrag Drag vector for screw component
 *
 *          OUT: dfScrewDragdv  Drag matrix, derivative of fScrewDrag for
 *                         screw component.
 *
 *************************************************************************/
void ScrewDrag2(real8 vel[3], real8 burg[3],real8 climbDir[3], 
                real8 dragClimb, real8 dragLine, real8 dragGlide,
                real8 fScrewDrag[3], real8 dfScrewDragdv[3][3])
{
        int   m, n;
        real8 lineDir[3], glideDir[3];
        real8 glideDirMat[3][3], climbDirMat[3][3], lineDirMat[3][3];

        VECTOR_COPY(lineDir, burg);
        NormalizeVec(lineDir);

        cross(climbDir, lineDir, glideDir);
        NormalizeVec(glideDir);

        Vec3TransposeAndMult(glideDir, glideDirMat);
        Vec3TransposeAndMult(climbDir, climbDirMat);
        Vec3TransposeAndMult(lineDir,  lineDirMat);
        
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfScrewDragdv[m][n] = dragClimb * climbDirMat[m][n] +
                                      dragLine  * lineDirMat[m][n]  +
                                      dragGlide * glideDirMat[m][n];
            }
        }

        Matrix33Vector3Multiply(dfScrewDragdv, vel, fScrewDrag);

        return;
}

void AngleDrag(real8 vel[3], 
	       real8 lineDir[3],
	       real8 climbDir[3], 
	       real8 dragClimb, 
	       real8 dragLine, 
	       real8 dragGlide_B, 
	       real8 dragGlide_D, 
	       real8 dragGlide_v0, 
	       real8 fAngleDrag[3], 
	       real8 dfAngleDragdv[3][3])
{
        int   m, n;
        real8 glideDir[3], glideForce[3];
        real8 glideDirMat[3][3], climbDirMat[3][3], lineDirMat[3][3];
     // real8 dragGlide = 0.0;
        real8 velg = 0.0;

	cross(climbDir, lineDir, glideDir);
        NormalizeVec(glideDir);

        Vec3TransposeAndMult(climbDir, climbDirMat);
        Vec3TransposeAndMult(lineDir,  lineDirMat);
        
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfAngleDragdv[m][n] = dragClimb*climbDirMat[m][n]+dragLine*lineDirMat[m][n];
            }
        }

        Matrix33Vector3Multiply(dfAngleDragdv, vel, fAngleDrag);

/*
 *      Velocity magnitude in gliding direction
 */
	velg = DotProduct(vel, glideDir); 

/*
 *      Evaluate non-linear gliding force vector and its derivative matrix
 */	
	NonlinearGlideDrag(velg, vel, dragGlide_B, dragGlide_D, dragGlide_v0, 
			   glideDir, glideForce, glideDirMat);

/*
 *      Add gliding contributions
 */	
	for (m = 0; m < 3; m++) {
	  for (n = 0; n < 3; n++) {
	    dfAngleDragdv[m][n] += glideDirMat[m][n];
	  }
	}
	
	for (m = 0; m < 3; m++)
	  fAngleDrag[m] += glideForce[m];
  
        return;
}

void NonlinearGlideDrag(real8 Vg, real8 *V, 
			real8 B, real8 D, real8 V0, 
			real8 *g, real8 *F, real8 dFdv[3][3])
{
  real8 Fmag  = 0.0;
  real8 dFdvg = 0.0;
  real8 sign_v = Vg/fabs(Vg);
  if(fabs(Vg) > V0) {
    Fmag  = B*Vg + 1.0*sign_v*D*pow((fabs(Vg)-V0),1.5);
    dFdvg = B    + 1.5*D*pow((fabs(Vg)-V0),0.5);
  } else {
    Fmag = B*Vg;
    dFdvg= B;
  }
  F[0] = Fmag*g[0];
  F[1] = Fmag*g[1];
  F[2] = Fmag*g[2];

  Vec3TransposeAndMult(g, dFdv);
  dFdv[0][0] *= dFdvg;
  dFdv[1][0] *= dFdvg;
  dFdv[2][0] *= dFdvg;
  dFdv[0][1] *= dFdvg;
  dFdv[1][1] *= dFdvg;
  dFdv[2][1] *= dFdvg;
  dFdv[0][2] *= dFdvg;
  dFdv[1][2] *= dFdvg;
  dFdv[2][2] *= dFdvg;

  // dFdv[0][0] = F[0]*g[0];  dFdv[0][1] = F[0]*g[1];  dFdv[0][2] = F[0]*g[2];
  // dFdv[1][0] = F[1]*g[0];  dFdv[1][1] = F[1]*g[1];  dFdv[1][2] = F[1]*g[2];
  // dFdv[2][0] = F[2]*g[0];  dFdv[2][1] = F[2]*g[1];  dFdv[2][2] = F[2]*g[2];

  // dFdv[0][0] = Fmag/(V[0]+EPSILON)*g[0];  dFdv[0][1] = Fmag/(V[0]+EPSILON)*g[1];  dFdv[0][2] = Fmag/(V[0]+EPSILON)*g[2];
  // dFdv[1][0] = Fmag/(V[1]+EPSILON)*g[0];  dFdv[1][1] = Fmag/(V[1]+EPSILON)*g[1];  dFdv[1][2] = Fmag/(V[1]+EPSILON)*g[2];
  // dFdv[2][0] = Fmag/(V[2]+EPSILON)*g[0];  dFdv[2][1] = Fmag/(V[2]+EPSILON)*g[1];  dFdv[2][2] = Fmag/(V[2]+EPSILON)*g[2];
  return;
}


/**************************************************************************
 *
 *      Function:     JunctionDrag
 *
 *      Description:  This function returns the drag for junction dislocation
 *                    segments.  Motion of junctions is restricted to the
 *                    corresponding line direction.
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
 *          OUT: dfJunctionDragdv  Drag matrix, derivative of fScrewDrag
 *                                 for junctions.
 *
 *************************************************************************/
void JunctionDrag(real8 vel[3], real8 lineDir[3],
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
