
#include "mpi_portability.h"

#include "Home.h"

/*---------------------------------------------------------------------------
 *
 *      Function:     OsmoticForce
 *      Description:  This routine calculates the osmotic force on the
 *                    end nodes of a dislocation segment associated with
 *                    a super saturation of vacancies in the crystal. The
 *                    force is in the climb direction of the dislocation if
 *                    the dislocation has any edge character, and the force
 *                    vanishes for screw dislocations.  The model follows
 *                    the framework proposed by Mordehai et al.
 *
 *                    The units of the parameters of the force are tied to
 *                    the other units of the simulations but typical units
 *                    are:
 *                        Boltmann's constant -> Joules/Kelvin
 *                        Atomic Volume - > in cubic burgers vectors
 *                        cvEquilibrium - > number density of vacancies in
 *                                          thermal equilibrium per atomic
 *                                          volume
 *                        cv - > number density of vacancies per atomic volume
 *
 *                    When this routine is used, there is a balance between
 *                    the climb of dislocations changes and the vacancy
 *                    concentation in the simulation volume.
 *
 *      Arguments:
 *          x1,y1,z1  Coordinates of segment's first endpoint
 *          x2,y2,z2  Coordinates of segment's second endpoint
 *          bx,by,bz  Components of burgers vector from first endpoint
 *                    to second endpoint
 *          f1        Array in which to return osmotic force at pos1
 *          f2        Array in which to return osmotic force at pos2
 *
 *      Author:       Tuan Hoang
 *
 *-------------------------------------------------------------------------*/
void OsmoticForce(Home_t *home, real8 x1, real8 y1, real8 z1,
                  real8 x2, real8 y2, real8 z2, real8 bx, real8 by, real8 bz,
                  real8 f1[3], real8 f2[3])
{
        real8    eps, osForce, bDotrt;
        real8    rtNorm, checkScrewNorm, bEdgeNorm;
        real8    temperature, atomicVol, cv, cvFloor, cvEquilibrium;
        real8    burg[3];
        real8    rt[3], rtUnit[3], burgUnit[3], checkScrew[3];
        real8    bEdge[3], tmp3[3];
        Param_t *param;

        eps = 1.0e-12;

        param = home->param;
        cv = param->vacancyConc;
        cvEquilibrium = param->vacancyConcEquilibrium;
        temperature = param->TempK;
        atomicVol = param->burgMag * param->burgMag * param->burgMag;

        burg[0] = bx;
        burg[1] = by;
        burg[2] = bz;

/*
 *      Get directional vector of the dislocation and the corresponding
 *      unit vector
 */
        rt[0] = x2 - x1;
        rt[1] = y2 - y1;
        rt[2] = z2 - z1;

        rtNorm = Normal(rt);

        rtUnit[0] = rt[0] / rtNorm;
        rtUnit[1] = rt[1] / rtNorm;
        rtUnit[2] = rt[2] / rtNorm;

        burgUnit[0] = burg[0];
        burgUnit[1] = burg[1];
        burgUnit[2] = burg[2];

        Normalize(&burgUnit[0], &burgUnit[1], &burgUnit[2]);

/*
 *      Determine the characteristics of the segment and calculate
 *      edge component of the burgers vector in the case of a 
 *      mixed dislocation.
 */
        cross(burgUnit, rtUnit, checkScrew);
        checkScrewNorm = Normal(checkScrew);

        if (checkScrewNorm > eps) {
            bDotrt = DotProduct(burg, rt);
            bEdge[0] = burg[0] - bDotrt * rt[0] / (rtNorm * rtNorm);
            bEdge[1] = burg[1] - bDotrt * rt[1] / (rtNorm * rtNorm);
            bEdge[2] = burg[2] - bDotrt * rt[2] / (rtNorm * rtNorm);
            bEdgeNorm = Normal(bEdge);
/*
 *          Calculate the normal of the glide plane containing burg and rt
 *          and the osmotic force (per unit length of the dislocation)
 *
 *          Note: Set a floor on the minimum cv value so the argument to the
 *          log() call below does not got to zero or less.
 */
            cvFloor = MAX(cv, cvEquilibrium * 1.0e-05);
            cross(bEdge, rt, tmp3);
            Normalize(&tmp3[0], &tmp3[1], &tmp3[2]);
            osForce = (-(BOLTZMANNS_CONST*temperature/atomicVol) *
                       bEdgeNorm*log(cvFloor/cvEquilibrium) * rtNorm*0.5);
        } else {
            osForce = 0.0;
            tmp3[0] = tmp3[1] = tmp3[2] = 0.0;
        }

/*
 *      Osmotic force per unit length of the segment at each endpoint...
 */
        f1[0] = osForce * tmp3[0]; 
        f1[1] = osForce * tmp3[1]; 
        f1[2] = osForce * tmp3[2]; 

        f2[0] = f1[0];
        f2[1] = f1[1];
        f2[2] = f1[2];
                     
        return;
}
