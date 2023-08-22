/**************************************************************************
 *
 *                               ALERT!
 *
 *    THE FUNCTIONS IN THIS MODULE ARE NOT CURRENTLY INVOKED BY ANY
 *    OF THE PARADIS CODE.  THE FUNCTIONS HAVE BEEN DEFINED FOR FUTURE
 *    USE.
 *
 *************************************************************************/


/**************************************************************************
 *
 *      Module:  StressDueToSemiInfSeg.cc
 *
 *      Description: This module contains a simple wrapper function
 *               to invoke the semi-infinite stress function
 *               appropriate to the manner in which the code has 
 *               been configured.
 *
 *      Includes public functions:
 *              StressDueToSemiInfSeg()
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"


#ifdef _ARLFEM
/**************************************************************************
 *
 *      Function:    StressDueToSemiInfSeg
 *      Description: Wrapper function which invokes an appropriate
 *                   stress function based on whether the code has
 *                   been compiled for isotropic or anisotropic 
 *                   elasticity.
 *
 *      Arguments:
 *         px, py, pz     coordinates of field point at which stress is to
 *                        be evaluated
 *         p1x, p1y, p1z  starting position of the dislocation segment
 *         p2x, p2y, p2z  position on the dislocation segment in the direction
 *                        from p1 to infinity
 *         bx, by, bz     burgers vector associated with segment going
 *                        from p1 to infinity
 *         a              core value
 *         MU             shear modulus
 *         NU             poisson ratio
 *         stress         array of stresses form the indicated segment
 *                        at the field point requested
 *                            [0] = stressxx
 *                            [1] = stressyy
 *                            [2] = stresszz
 *                            [3] = stressyz
 *                            [4] = stressxz
 *                            [5] = stressxy
 *
 *************************************************************************/
void StressDueToSemiInfSeg(Home_t *home,
                           real8 px, real8 py, real8 pz,
                           real8 p1x, real8 p1y, real8 p1z,
                           real8 p2x, real8 p2y, real8 p2z,
                           real8 bx, real8 by, real8 bz,
                           real8 a, real8 MU, real8 NU, real8 *sigma)
{

#ifdef ANISOTROPIC
    StressDueToSemiInfSegAnisotropic(home, px, py, pz,
                                     p1x, p1y, p1z, p2x, p2y, p2z,
                                     bx, by, bz, a,
                                     home->param->anisoHarmonicsNumTermsBase,
                                     sigma);
#else
    StressDueToSemiInfSegIsotropic(px, py, pz,
                                   p1x, p1y, p1z, p2x, p2y, p2z,
                                   bx, by, bz, a, MU, NU, sigma);
#endif
}
#endif  // end ifdef _ARLFEM
