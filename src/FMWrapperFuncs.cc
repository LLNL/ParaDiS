/**************************************************************************
 *
 *      Module:  FMWrapperFuncs.c 
 *               Contains wrapper functions for the FMM.
 *               Calls either the Taylor or the Black Box FMM methods.
 *
 *      Includes functions:
 *               SegmentToMoment
 *               Moment2Moment
 *               Local2Local
 *               Moment2Local
 *               Local2Particle
 *
 *
 *
 *************************************************************************/
#include "Home.h"
#include "FM.h"


void SegmentToMoment(Home_t *home, real8 *cellSize, real8 *cellCenter,
                     real8 *vec1, real8 *vec2, real8 *burg,
                     real8 *moments)
{
#ifdef TAYLORFMM
  makeeta(home->param->fmMPOrder, vec1, vec2, burg, moments);
#endif

#ifdef BBFMM
/*
 *  Calculate segment endpoints from the provided vectors
 */
  real8 L;
  real8 p1[3], p2[3];
  p1[0] = cellCenter[0] + vec1[0];
  p1[1] = cellCenter[1] + vec1[1];
  p1[2] = cellCenter[2] + vec1[2];

  p2[0] = p1[0] + vec2[0];
  p2[1] = p1[1] + vec2[1];
  p2[2] = p1[2] + vec2[2];

  L = cellSize[0]; /* assuming cell is cubic, use cell size in X dimension */

  bbfmm_S2M(p1, p2, cellCenter, burg, L, moments); 

#endif  /* ifdef BBFMM */


#ifdef UNIFORMFMM
/*
 *  Calculate segment endpoints from the provided vectors
 */
  real8 L;
  real8 p1[3], p2[3];
  p1[0] = cellCenter[0] + vec1[0];
  p1[1] = cellCenter[1] + vec1[1];
  p1[2] = cellCenter[2] + vec1[2];

  p2[0] = p1[0] + vec2[0];
  p2[1] = p1[1] + vec2[1];
  p2[2] = p1[2] + vec2[2];

  L = cellSize[0]; /* assuming cell is cubic, use cell size in X dimension */

  uniformfmm_S2M(p1, p2, cellCenter, burg, L, moments); 

#endif  /* ifdef BBFMM */

}


void MomentToMoment(Home_t *home, real8 *r, real8 *eta, real8 *neta)
{
#ifdef TAYLORFMM
  FMShift(home->param->fmMPOrder, r, eta, neta);
#endif
  
#ifdef BBFMM
  bbfmm_M2M(r, eta, neta);
#endif


#ifdef UNIFORMFMM
  uniformfmm_M2M(r, eta, neta);
#endif
}


void LocalToLocal(Home_t *home, real8 *r, real8 *alpha, real8 *beta)
{
#ifdef TAYLORFMM
  TaylorShift(home->param->fmExpansionOrder, r, alpha, beta);
#endif
  
#ifdef BBFMM
  bbfmm_L2L(r, alpha, beta);
#endif

#ifdef UNIFORMFMM
  uniformfmm_L2L(r, alpha, beta);
#endif
}



void MomentToLocal(Home_t *home, real8 *cellSize, real8 *r, real8 *eta,
                  real8 *alpha)
{
#ifdef TAYLORFMM
  int     maxOrder;
  Param_t *param;

  param = home->param;
  maxOrder = MAX(param->fmMPOrder, param->fmExpansionOrder);

  MkTaylor(home, param->shearModulus, param->pois, param->fmMPOrder,
           param->fmExpansionOrder, maxOrder, r, eta, alpha);
#endif
  
#ifdef BBFMM
  real8 L;

  L = cellSize[0]; /* assuming cell is cubic, use cell size in X dimension */
  bbfmm_M2L(r, L, eta, alpha);
#endif

#ifdef UNIFORMFMM
  real8 L;

  L = cellSize[0]; /* assuming cell is cubic, use cell size in X dimension */
  uniformfmm_M2L(r, L, eta, alpha);
#endif

}


void LocalToPoint(Home_t *home, real8 *cellSize, real8 *r, real8 *alpha,
                 real8 sigma[3][3])
{
#ifdef TAYLORFMM
  /* sigma is 3x3 */
  EvalTaylor(home->param->fmExpansionOrder, r, alpha, sigma);
#endif
  
#ifdef BBFMM
  real8 L;
  L = cellSize[0]; /* assuming cell is cubic, use cell size in X dimension */
  /* sigma is 3x3 */
  bbfmm_L2P(r, L, alpha, sigma);
#endif

#ifdef UNIFORMFMM
  real8 L;
  L = cellSize[0]; /* assuming cell is cubic, use cell size in X dimension */
  /* sigma is 3x3 */
  uniformfmm_L2P(r, L, alpha, sigma);
#endif

}



/* Wrapper Functions for PBC */
void CorrectionTableInit(Home_t *home)
{
#ifdef TAYLORFMM
  Taylor_CorrectionTableInit(home);
#else
  Stanford_CorrectionTableInit(home);
#endif

}

void DoTableCorrection(Home_t *home)
{
#ifdef TAYLORFMM
  Taylor_DoTableCorrection(home);
#endif

#ifdef BBFMM
    bbFMM_DoTableCorrection(home);
#endif

#ifdef UNIFORMFMM
    uniformFMM_DoTableCorrection(home);
#endif
}


void FreeCorrectionTable(void )
{
#ifdef TAYLORFMM
  Taylor_FreeCorrectionTable();
#else
  Stanford_FreeCorrectionTable();
#endif
}


void MeanStressCorrection(Home_t *home)
{
#ifdef TAYLORFMM
  Taylor_MeanStressCorrection(home);
#else

  int n;
  real8 *Tkz;
#ifdef BBFMM
  n = bbfmm->n;
  Tkz = bbfmm->Tkz;
#else
  n = uniformfmm->n;
  Tkz = uniformfmm->Tkz;
#endif
  
  Stanford_MeanStressCorrection(home,n, Tkz);
#endif

}
