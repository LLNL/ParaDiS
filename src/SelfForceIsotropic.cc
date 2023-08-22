/***************************************************************************
 *
 *      Module:      SelfForceIsotropic.c
 *
 *      Description: This module contains the functions needed for
 *                   calculating isotropic self-force for a given
 *                   dislocation segment.
 *
 *      Includes public functions:
 *          SelfForceIsotropic()
 *          CoreForceandDerivativeAnisotropic()
 *
 ***************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

#if 0
/*--------------------------------------------------------------------------
 *
 *      Function:        CoreForceIsotropic2
 *      Description:     This function is for debug and validation only
 *
 *--------------------------------------------------------------------------*/
void CoreForceIsotropic2(Home_t *home, real8 x1, real8 y1, real8 z1,
			real8 x2, real8 y2, real8 z2,
			real8 b[3], real8 Ecore, real8 NU,
			real8 f2Core[3], real8 dfCoredx2[3][3])
{
		real8 L;
		real8 t[3];
        real8 Oneov1mNu = 1./(1-NU);
        Param_t *param = home->param;

		GetUnitVector(1, x1, y1, z1, x2, y2, z2, &(t[0]), &(t[1]), &(t[2]), &L);

		real8 b2 = DotProduct(b, b);
		real8 bs = b[0]*t[0] + b[1]*t[1] + b[2]*t[2];

		real8 be[3];
		be[0] = b[0]-bs*t[0];
		be[1] = b[1]-bs*t[1];
		be[2] = b[2]-bs*t[2];

		real8 cost = DotProduct(b, t) / sqrt(b2);
		real8 cos2t = cost*cost;
		real8 E = Ecore*b2*(cos2t + (1.0-cos2t)*Oneov1mNu);
		real8 dEdcos2 = Ecore*b2*(1.0-Oneov1mNu);

		f2Core[0] = -dEdcos2*2.0*bs/b2*be[0] -E*t[0];
		f2Core[1] = -dEdcos2*2.0*bs/b2*be[1] -E*t[1];
		f2Core[2] = -dEdcos2*2.0*bs/b2*be[2] -E*t[2];


		/* Derivative */
		real8 d2Edcos2 = 0.0;

		real8 bexbe[3][3], Itxt[3][3];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				bexbe[i][j] = be[i]*be[j];
				Itxt[i][j] = -t[i]*t[j];
			}
			Itxt[i][i] += 1.0;
		}

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				dfCoredx2[i][j] = -4.0*bs*bs/b2/b2/L*d2Edcos2*bexbe[i][j]
				                  +2.0/b2/L*dEdcos2*(-bexbe[i][j]+bs*bs*Itxt[i][j])
								  -E/L*Itxt[i][j];
			}
		}
}
#endif

/*--------------------------------------------------------------------------
 *
 *      Function:        CoreEnergyMD
 *      Description:     Core energy function fitted to MD data
 *
 *--------------------------------------------------------------------------*/
void CoreEnergyMD(Param_t *param, int bIndex, int pIndex, real8 theta,
	              real8 *E, real8 *dEdcos2, real8 *d2Edcos2)
{
		// Determine dislocation type
		int cType = -1;
		if (bIndex <= 3) {
			/* 1/2[111] glide dislocation */
			if ((pIndex >=  0 && pIndex <= 2 ) ||
			    (pIndex >=  6 && pIndex <= 8 ) ||
			    (pIndex >= 12 && pIndex <= 14) ||
			    (pIndex >= 18 && pIndex <= 20)) {
				/* 1/2[111](110) glide dislocation */
				cType = 0;
			} else if ((pIndex >=  3 && pIndex <= 5 ) ||
			           (pIndex >=  9 && pIndex <= 11) ||
			           (pIndex >= 15 && pIndex <= 17) ||
			           (pIndex >= 21 && pIndex <= 23)) {
				/* 1/2[111](112) glide dislocation */
				cType = 1;
			} else {
				/* Use 1/2[111](112) as default otherwise */
				cType = 1;
			}
		} else if (bIndex >= 4 && bIndex <= 6) {
			/* [100] binary junction */
			if (pIndex == 24 || pIndex == 25 || pIndex == 40 ||
			    pIndex == 41 || pIndex == 56 || pIndex == 57) {
				/* [100](110) dislocation junction */
				cType = 4;
			} else if (pIndex == 26 || pIndex == 27 || pIndex == 42 ||
			           pIndex == 43 || pIndex == 58 || pIndex == 59) {
				/* [100](010) dislocation junction */
				cType = 3;
			} else {
				/* [100](110) dislocation junction */
				cType = 4;
			}
		} else {
			/* Default */
			/* [100](110) dislocation junction */
			cType = 4;
		}


		real8 cost = cos(theta);
		real8 cos2t = cost*cost;

		/* Polynomial fitting */
		int   porder=0;
		real8 pcoeffs[10];

#if 0
		// Default model fitted (testing only - OK)
		if (theta >= 0.0) {
			porder = 1;
			pcoeffs[0] = -0.201519660755793;
			pcoeffs[1] = 0.594453276565760;
		} else {
			porder = 1;
			pcoeffs[0] = -0.201519660755793;
			pcoeffs[1] = 0.594453276565760;
		}
#else
		if (StrEquiv(param->MDcoreMaterial, "W")) {
			if (cType == 0) {
				/* W 1/2[111](110) */
				if (theta <= 0.0) {
					porder = 8;
					pcoeffs[0] = -1.804090167929550e+02;
					pcoeffs[1] = 6.477005814623024e+02;
					pcoeffs[2] = -9.267659610664427e+02;
					pcoeffs[3] = 6.700171014668820e+02;
					pcoeffs[4] = -2.548950234629997e+02;
					pcoeffs[5] = 47.477852915591185;
					pcoeffs[6] = -3.250603343349422;
					pcoeffs[7] = -0.386883184301567;
					pcoeffs[8] = 0.688220907173072;
				} else {
					porder = 8;
					pcoeffs[0] = -1.836627311347282e+02;
					pcoeffs[1] = 6.823472980373073e+02;
					pcoeffs[2] = -1.028479119971369e+03;
					pcoeffs[3] = 8.064740523651042e+02;
					pcoeffs[4] = -3.504435641840109e+02;
					pcoeffs[5] = 82.760570054439300;
					pcoeffs[6] = -9.846789932003803;
					pcoeffs[7] = 0.333934713189888;
					pcoeffs[8] = 0.691990059068584;
				}
			} else if (cType == 1) {
				/* W 1/2[111](112) */
				porder = 8;
				pcoeffs[0] = -1.446262220478938e+02;
				pcoeffs[1] = 5.499464725987633e+02;
				pcoeffs[2] = -8.519044267741900e+02;
				pcoeffs[3] = 6.911207142482818e+02;
				pcoeffs[4] = -3.144311414973766e+02;
				pcoeffs[5] = 79.560978962113580;
				pcoeffs[6] = -10.340027921845001;
				pcoeffs[7] = 0.114772126699434;
				pcoeffs[8] = 0.734459777834759;
			} else if (cType == 3) {
				/* W [100](010) */
				porder = 8;
				pcoeffs[0] = -3.180292783906065e+02;
				pcoeffs[1] = 1.321353671685343e+03;
				pcoeffs[2] = -2.262031059699275e+03;
				pcoeffs[3] = 2.060506026530759e+03;
				pcoeffs[4] = -1.077020397616367e+03;
				pcoeffs[5] = 3.247274043824662e+02;
				pcoeffs[6] = -53.822279841487685;
				pcoeffs[7] = 4.728706630743911;
				pcoeffs[8] = 0.538164809461604;
			} else if (cType == 4) {
				/* W [100](110) */
				if (fabs(theta) <= 0.9554 /* 54.74 deg */) {
					porder = 6;
					pcoeffs[0] = -1.707534975295302e+02;
					pcoeffs[1] = 6.879666960509496e+02;
					pcoeffs[2] = -1.129208362035896e+03;
					pcoeffs[3] = 9.660163918602670e+02;
					pcoeffs[4] = -4.546982764779750e+02;
					pcoeffs[5] = 1.122270219454686e+02;
					pcoeffs[6] = -10.603149328944712;
				} else {
					porder = 4;
					pcoeffs[0] = 27.410013057080224;
					pcoeffs[1] = -23.956509035087596;
					pcoeffs[2] = 6.193247421728476;
					pcoeffs[3] = -0.898251211463245;
					pcoeffs[4] = 0.878250608895093;
				}
			}
		} else if (StrEquiv(param->MDcoreMaterial, "Ta")) {
			if (cType == 0) {
				/* Ta 1/2[111](110) */
				if (theta <= 0.0) {
					porder = 8;
					pcoeffs[0] = -1.248039e+02;
					pcoeffs[1] = 4.736619e+02;
					pcoeffs[2] = -7.292838e+02;
					pcoeffs[3] = 5.896457e+02;
					pcoeffs[4] = -2.737382e+02;
					pcoeffs[5] = 7.590514e+01;
					pcoeffs[6] = -1.282221e+01;
					pcoeffs[7] = 7.389047e-01;
					pcoeffs[8] = 1.005656e+00;
				} else {
					porder = 8;
					pcoeffs[0] = -3.721121e+00;
					pcoeffs[1] = -4.449826e+01;
					pcoeffs[2] = 1.790668e+02;
					pcoeffs[3] = -2.536884e+02;
					pcoeffs[4] = 1.751634e+02;
					pcoeffs[5] = -6.271557e+01;
					pcoeffs[6] = 1.100602e+01;
					pcoeffs[7] = -1.318097e+00;
					pcoeffs[8] = 1.014785e+00;
				}
			} else if (cType == 1) {
				/* Ta 1/2[111](112) */
				porder = 8;
				pcoeffs[0] = -2.637448e+02;
				pcoeffs[1] = 1.050652e+03;
				pcoeffs[2] = -1.700436e+03;
				pcoeffs[3] = 1.434741e+03;
				pcoeffs[4] = -6.743630e+02;
				pcoeffs[5] = 1.756626e+02;
				pcoeffs[6] = -2.433798e+01;
				pcoeffs[7] = 1.111112e+00;
				pcoeffs[8] = 1.008781e+00;
			} else if (cType == 3) {
				/* Ta [100](010) */
				porder = 8;
				pcoeffs[0] = 1.346011e+00;
				pcoeffs[1] = -2.007501e+01;
				pcoeffs[2] = 3.449267e+01;
				pcoeffs[3] = -9.283937e-01;
				pcoeffs[4] = -3.771889e+01;
				pcoeffs[5] = 3.323840e+01;
				pcoeffs[6] = -1.203660e+01;
				pcoeffs[7] = 2.202739e+00;
				pcoeffs[8] = 7.474235e-01;
			} else if (cType == 4) {
				/* Ta [100](110) */
				if (fabs(theta) <= 0.9554 /* 54.74 deg */) {
					porder = 6;
					pcoeffs[0] = 4.328633e+01;
					pcoeffs[1] = -1.870063e+02;
					pcoeffs[2] = 3.307777e+02;
					pcoeffs[3] = -3.043766e+02;
					pcoeffs[4] = 1.529259e+02;
					pcoeffs[5] = -3.919007e+01;
					pcoeffs[6] = 4.910192e+00;
				} else {
					porder = 4;
					pcoeffs[0] = -1.250434e+01;
					pcoeffs[1] = 1.420643e+01;
					pcoeffs[2] = -5.607278e+00;
					pcoeffs[3] = 6.028495e-01;
					pcoeffs[4] = 9.890963e-01;
				}
			}
		} else {
			Fatal("CoreEnergyMD: unknown MDcoreMaterial = %s", param->MDcoreMaterial);
		}
#endif

		*E = pcoeffs[porder];
		*dEdcos2 = 0.0;
		real8 x = 1.0;
		for (int i = porder-1; i >= 0; i--) {
			*dEdcos2 += pcoeffs[i]*x;
			x *= cos2t;
			*E += pcoeffs[i]*x;
		}

		/* Second derivative */
		*d2Edcos2 = 0.0;
		if (porder > 1) {
			x = 1.0;
			for (int i = porder-2; i >= 0; i--) {
				*d2Edcos2 += pcoeffs[i]*x;
				x *= cos2t;
			}
		}

		/* Special cases */
		if (theta == 0.0 || (cType == 4 && fabs(theta) == 0.9554)) {
			*dEdcos2 = 0.0;
			*d2Edcos2 = 0.0;
		}

		// Units conversion from eV/A to J/b^2
		real8 convert = 1.6022e-19 * 1e10 / param->burgMag / param->burgMag;
		*E = (*E) * convert;
		*dEdcos2 = (*dEdcos2) * convert;
		*d2Edcos2 = (*d2Edcos2) * convert;
}

/*--------------------------------------------------------------------------
 *
 *      Function:        CoreForceandDerivativeMD
 *      Description:     Calculate the core force and derivative
 *		                 fitted to MD data
 *
 *--------------------------------------------------------------------------*/
void CoreForceandDerivativeMD(Home_t *home, real8 x1, real8 y1, real8 z1,
			real8 x2, real8 y2, real8 z2,
			real8 b[3], real8 Ecore, real8 NU,
			real8 f2Core[3], real8 dfCoredx2[3][3])
{
		real8 L;
		real8 t[3], burg[3], n[3];
		real8 E=0.0, dEdcos2=0.0, d2Edcos2=0.0;
        Param_t *param = home->param;

/*
 *      It is possible that this function gets called for a
 *		zero Burgers vector segment, e.g. during SplitMultiNodes.
 *		Just return here with a zero core energy.
 */
		if (fabs(DotProduct(b, b)) < 1.0e-10) {
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++)
					dfCoredx2[i][j] = 0.0;
				f2Core[i] = 0.0;
			}
			return;
		}


		GetUnitVector(1, x1, y1, z1, x2, y2, z2, &(t[0]), &(t[1]), &(t[2]), &L);

		real8 b2 = DotProduct(b, b);
		real8 bs = b[0]*t[0] + b[1]*t[1] + b[2]*t[2];

		real8 be[3];
		be[0] = b[0]-bs*t[0];
		be[1] = b[1]-bs*t[1];
		be[2] = b[2]-bs*t[2];

		burg[0] = b[0];
		burg[1] = b[1];
		burg[2] = b[2];

		FindPreciseGlidePlane(home, burg, t, n,
							  home->param->allowFuzzyGlidePlanes);

/*
 *      If the glide plane is sufficiently close to zero, the
 *      segment is screw, so explicitly select a default plane
 *		(the core energy of the screw is the same for all planes)
 */
		if (fabs(DotProduct(n, n)) < 1.0e-2) {

			int bIndex = -1;
            real8 (*burgList)[3];
            burgList = home->burgData.burgList;
            bIndex = GetBurgIndexBCC(b, BCC_NUM_TOT_BURG, burgList);

            if (bIndex < BCC_NUM_GLIDE_BURG) {
				/* [111] Burgers */
				n[0] = 2.0*b[0];
				n[1] = -b[1];
				n[2] = -b[2];
			} else {
				/* Junction Burgers */
				if (fabs(b[2]) > 1e-10) {
					n[0] = 1.0;
					n[1] = 1.0;
					n[2] = -(b[0]+b[1])/b[2];
		        } else if (fabs(b[1]) > 1e-10) {
		            n[0] = 1.0;
					n[1] = -(b[0]+b[2])/b[1];
					n[2] = 1.0;
		        } else {
		            n[0] = -(b[1]+b[2])/b[0];
					n[1] = 1.0;
					n[2] = 1.0;
		        }
			}
			NormalizeVec(n);
		}

/*
 *      If needed, rotate the line sense from the laboratory frame to
 *      the crystal frame.
 */
		if (param->useLabFrame) {
			double tmpRot[3];

			Matrix33Vector3Multiply(param->rotMatrixInverse, burg, tmpRot);
			burg[0] = tmpRot[0];
			burg[1] = tmpRot[1];
			burg[2] = tmpRot[2];

			Matrix33Vector3Multiply(param->rotMatrixInverse, n, tmpRot);
			n[0] = tmpRot[0];
			n[1] = tmpRot[1];
			n[2] = tmpRot[2];
		}

		/* Get signed character angle */
		if (param->materialType == MAT_TYPE_BCC && param->rc == 1.0) {

			int bIndex, pIndex, gIndex;
			double bref[3], nref[3], cref[3], tx, ty, theta;

			// Normalize Burgers and plane normals to get indices
			NormalizeVec(b);
			NormalizeVec(n);

			GetBCCAllIndices(home, b, n, &bIndex, &pIndex, &gIndex);

			bref[0] = home->burgData.burgList[bIndex][0];
			bref[1] = home->burgData.burgList[bIndex][1];
			bref[2] = home->burgData.burgList[bIndex][2];

			if (pIndex != -1) {
				nref[0] = home->burgData.planeList[pIndex][0];
				nref[1] = home->burgData.planeList[pIndex][1];
				nref[2] = home->burgData.planeList[pIndex][2];
			} else {
				// Use the original plane instaed for now on
				nref[0] = n[0];
				nref[1] = n[1];
				nref[2] = n[2];
			}

			cross(nref, b, cref);
			NormalizeVec(cref);

			tx = DotProduct(t, b);
			ty = DotProduct(t, cref);
			theta = atan2(ty, tx); // signed character angle

			// Fold theta in [-Pi/2,Pi/2]
			if (theta >  0.5*M_PI) theta = theta - M_PI;
			if (theta < -0.5*M_PI) theta = theta + M_PI;
			theta = MAX(MIN(theta, 0.5*M_PI), -0.5*M_PI);

			CoreEnergyMD(param, bIndex, pIndex, theta, &E, &dEdcos2, &d2Edcos2);

		} else {
			Fatal("CoreForceIsotropicMD function must be used with BCC crystal with rc = 1.0");
		}

		f2Core[0] = -dEdcos2*2.0*bs/b2*be[0] -E*t[0];
		f2Core[1] = -dEdcos2*2.0*bs/b2*be[1] -E*t[1];
		f2Core[2] = -dEdcos2*2.0*bs/b2*be[2] -E*t[2];

		/* Derivative */
		real8 bexbe[3][3], Itxt[3][3];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				bexbe[i][j] = be[i]*be[j];
				Itxt[i][j] = -t[i]*t[j];
			}
			Itxt[i][i] += 1.0;
		}

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				dfCoredx2[i][j] = -4.0*bs*bs/b2/b2/L*d2Edcos2*bexbe[i][j]
				                  +2.0/b2/L*dEdcos2*(-bexbe[i][j]+bs*bs*Itxt[i][j])
								  -E/L*Itxt[i][j];
			}
		}
}

/*--------------------------------------------------------------------------
 *
 *      Function:        CoreForceMD
 *      Description:     Calculate the core force fitted to MD data
 *
 *--------------------------------------------------------------------------*/
void CoreForceMD(Home_t *home, real8 x1, real8 y1, real8 z1,
				     real8 x2, real8 y2, real8 z2,
				     real8 b[3], real8 Ecore, real8 NU,
				     real8 f2Core[3])
{
		real8 dfCoredx2[3][3];
		CoreForceandDerivativeMD(home, x1, y1, z1, x2, y2, z2,
				                 b, Ecore, NU, f2Core, dfCoredx2);
}


/*---------------------------------------------------------------
 *
 *      Function:        CoreForceIsotropic
 *      Description:     Calculate the core force for
 *                       isotropic elasticity.
 *
 *      Parameters:
 *          IN:  x1,y1,z1, x2,y,2,z2  Dislocation segment
 *          IN:  b                    Burgers vector
 *          IN:  Ecore                Core energy
 *          OUT: fCore                Force coming from core effects at node 2
 *
 *--------------------------------------------------------------------------*/
void CoreForceIsotropic(Home_t *home, real8 x1, real8 y1, real8 z1,
			real8 x2, real8 y2, real8 z2,
			real8 b[3], real8 Ecore, real8 NU,
			real8 f2Core[3])
{
#if 1
	if (home->param->MDcore) {
		CoreForceMD(home, x1, y1, z1, x2, y2, z2,
				    b, Ecore, NU, f2Core);
		return;
	}
#endif

	real8 L;
	real8 t[3];
        real8 bs, bex, bey, bez, be2, bs2;
        real8 fL,ft;
        real8 Oneov1mNu = 1./(1-NU);
        Param_t *param = home->param;

/*
 *      For HCP Be, we have MD data that gives Ecore as a function of Burgers
 *      vector type, so replace the Ecore value as appropriate
 */
        if (param->materialType == MAT_TYPE_HCP) {
            real8 bCryst[3];

/*
 *          If needed, rotate the burgers vector from the laboratory frame
 *          to the crystal frame.
 */
            if (param->useLabFrame) {
                real8 bLab[3];
                bLab[X] = b[0];
                bLab[Y] = b[1];
                bLab[Z] = b[2];
                Matrix33Vector3Multiply(param->rotMatrixInverse, bLab, bCryst);
            } else {
                bCryst[X] = b[0];
                bCryst[Y] = b[1];
                bCryst[Z] = b[2];
            }

	    Ecore = HCPEcoreFactor(home, bCryst);

        }  /* end if (param->materialType == MAT_TYPE_HCP) */


        GetUnitVector(1, x1, y1, z1, x2, y2, z2, &(t[0]), &(t[1]), &(t[2]), &L);

        bs = b[0]*t[0] + b[1]*t[1] + b[2]*t[2];
        bex = b[0]-bs*t[0]; bey = b[1]-bs*t[1]; bez=b[2]-bs*t[2];
        be2 = (bex*bex+bey*bey+bez*bez);
        bs2 = bs*bs;

        fL = -Ecore*(bs2+be2*Oneov1mNu);
        ft =  Ecore*2*bs*NU*Oneov1mNu;

/*
 *      Deduce core force  F(x2) = -d (E . L) dx where L=x2-x1 is segment length
 */
        f2Core[0] = bex*ft + fL*t[0];
        f2Core[1] = bey*ft + fL*t[1];
        f2Core[2] = bez*ft + fL*t[2];
}


/*---------------------------------------------------------------
 *
 *      Function:        CoreForceandDerivativeAnisotropic
 *      Description:     Calculate the core force and its derivative for
 *                       isotropic elasticity.  The force derivative is
 *                       only computed when requested.
 *
 *      Parameters:
 *          IN:  t        Dislocation line direction (normalized)
 *          IN:  bx,by,bz Burgers vector
 *          IN:  C        Elastic constants provided in the 3x3x3x3
 *                        matrix form
 *          IN:  calcForceDeriv If 1, derivative of core force is
 *                        computed. Returns dFdt = 0 otherwise.
 *          OUT: fCore    Force coming from core effects
 *          OUT: dfCoredt Derivative of Force coming from core effects
 *
 *--------------------------------------------------------------------------*/

void CoreForceandDerivativeIsotropic(Home_t *home, real8 x1, real8 y1, real8 z1,
				     real8 x2, real8 y2, real8 z2,
				     real8 b[3], real8 Ecore, real8 NU,
				     real8 f2Core[3], real8 dfCoredx2[3][3])
{
#if 1
	if (home->param->MDcore) {
		CoreForceandDerivativeMD(home, x1, y1, z1, x2, y2, z2,
				                 b, Ecore, NU, f2Core, dfCoredx2);
		return;
	}
#endif

		int i,j;
	real8 t[3];
        real8 bs, bex, bey, bez, be2, bs2;
        real8 fL,ft;
        real8 Oneov1mNu = 1./(1-NU);
        Param_t *param = home->param;
	real8 dtdx2[3][3], dbsdx2[3], dbedx2[3][3], dbe2dx2[3];
	real8 Twobs, val, val2;
	real8 fac2 = 2*NU*Oneov1mNu;
	real8 be[3];
	real8 L, invL, tmp[3], fbe;


        GetUnitVector(1, x1, y1, z1, x2, y2, z2, &(t[0]), &(t[1]), &(t[2]), &L);

/*
 *      For HCP Be, we have MD data that gives Ecore as a function of Burgers
 *      vector type, so replace the Ecore value as appropriate
 */
        if (param->materialType == MAT_TYPE_HCP) {
            real8 bCryst[3];

/*
 *          If needed, rotate the burgers vector from the laboratory frame
 *          to the crystal frame.
 */
            if (param->useLabFrame) {
                real8 bLab[3];
                bLab[X] = b[0];
                bLab[Y] = b[1];
                bLab[Z] = b[2];
                Matrix33Vector3Multiply(param->rotMatrixInverse, bLab, bCryst);
            } else {
                bCryst[X] = b[0];
                bCryst[Y] = b[1];
                bCryst[Z] = b[2];
            }

	    Ecore = HCPEcoreFactor(home, bCryst);

        }  /* end if (param->materialType == MAT_TYPE_HCP) */

/*
 *      In case we have a zero-length vector...
 */
        if (L == 0.0) {
            invL = 0.0;
        } else {
            invL= 1./L;
        }

        bs = b[0]*t[0] + b[1]*t[1] + b[2]*t[2];
        bex = b[0]-bs*t[0]; bey = b[1]-bs*t[1]; bez=b[2]-bs*t[2];
        be2 = (bex*bex+bey*bey+bez*bez);
        bs2 = bs*bs;

        fL = -Ecore*(bs2+be2*Oneov1mNu);
        ft =  Ecore*2*bs*NU*Oneov1mNu;

/*
 *      Deduce core force  F(x1) = -d (E . L) dx where L=x2-x1 is segment length
 */
        f2Core[0] = bex*ft + fL*t[0];
        f2Core[1] = bey*ft + fL*t[1];
        f2Core[2] = bez*ft + fL*t[2];

/*
 *      Compute the derivative now
 */

	be[0]=bex;be[1]=bey;be[2]=bez;

	dtdx2[0][0] = (1.0 - t[0]*t[0]) * invL;
	dtdx2[1][1] = (1.0 - t[1]*t[1]) * invL;
	dtdx2[2][2] = (1.0 - t[2]*t[2]) * invL;

	dtdx2[0][1] = -t[0]*t[1] * invL;
	dtdx2[0][2] = -t[0]*t[2] * invL;

	dtdx2[1][0] = -t[1]*t[0] * invL;
	dtdx2[1][2] = -t[1]*t[2] * invL;

	dtdx2[2][0] = -t[2]*t[0] * invL;
	dtdx2[2][1] = -t[2]*t[1] * invL;

	Matrix33_Vmul(dbsdx2,b,dtdx2);

	Twobs = 2*bs;
	for (i=0;i<3;i++)
	  {
	    for (j=0;j<3;j++)
	      dbedx2[i][j]  = -t[i]*dbsdx2[j] - bs*dtdx2[i][j];

	    dbe2dx2[i] = -Twobs*dbsdx2[i];
	  }

	val = bs2+be2*Oneov1mNu;
	val2= bs*fac2;

	for (j=0;j<3;j++)
	  tmp[j] = Twobs*dbsdx2[j]+Oneov1mNu*dbe2dx2[j];

	/* This function computes the derivatives at node 2 */
	for (i=0; i<3;i++)
	  {
	    fbe = fac2*be[i];
	    for (j=0;j<3;j++)
	      {
		dfCoredx2[i][j]=-Ecore*(t[i]*tmp[j]+ val*dtdx2[i][j]-val2*dbedx2[i][j]-fbe*dbsdx2[j]);
	      }
	  }
}


void SelfForceIsotropic(Home_t *home, int coreOnly, real8 MU, real8 NU,
                        real8 bx, real8 by, real8 bz,
                        real8 x1, real8 y1, real8 z1,
                        real8 x2, real8 y2, real8 z2,
                        real8 a,  real8 Ecore,
                        real8 f1[3], real8 f2[3])
{
        real8 tx, ty, tz, L, La, S;
        real8 bs, bex, bey, bez;
        real8 Lainv, Linv, Scorr = 0.0;
	real8 b[3], f2Core[3];

/*
 *      If the segment was zero-length, there is no self-force, so
 *      zero the forces and return to the caller.
 */
        GetUnitVector(1, x1, y1, z1, x2, y2, z2, &tx, &ty, &tz, &L);

        if (L <= 0.0) {

            VECTOR_ZERO(f1);
            VECTOR_ZERO(f2);

            return;
        }

        bs = bx*tx + by*ty + bz*tz;
        bex = bx-bs*tx; bey = by-bs*ty; bez=bz-bs*tz;

        La=sqrt(L*L+a*a);
        Linv=1.0e0/L;
        Lainv=1.0e0/La;

/*
 *      NOTE: In the standard self-force calculation there
 *      is a divergence of the stress field on an element by
 *      element basis that cancels when the total segment
 *      configuration closes on itself or is bounded by semi-
 *      infinite segments.
 *
 *      If the USE_RATIONAL_SEG_FORCES flag has been set during
 *      compilation, the calculation includes a correction that makes
 *      the stress field due to a dislocation segment divergence free
 *      on an element by element basis so the system does not have
 *      to be closed or terminated by semi-infinite segments in order
 *      to obtain the divergence-free stress field.
 */
        if (coreOnly) {
            S = 0.0;
        } else {
            S = (-(2.0e0*NU*La+(1-NU)*a*a*Lainv-(1+NU)*a)*Linv +
                 (NU*log((La+L)/a)-(1-NU)*0.5e0*L*Lainv)) *
                MU*0.25e0/M_PI/(1-NU)*bs;
#ifdef USE_RATIONAL_SEG_FORCES
            Scorr = -(0.125*MU/M_PI)*bs*(log((La+L)/a)+L*Lainv-4*(La-a)*Linv);
#endif
            S = S + Scorr;
        }

/*
 *      Ecore = MU/(4*pi) log(a/a0)
 *
 *      Account for the selfforce due to the core energy change when
 *      the core radius changes from a0-->a, M. Tang, 7/19/2004
 */
	b[0] = bx;
        b[1] = by;
        b[2] = bz;

	/* To test aniso, set Ecore = MU/4/M_PI;*/

        CoreForceIsotropic(home, x1, y1, z1, x2, y2, z2, b,
			   Ecore, NU, f2Core);


	f2[0] = bex*S + f2Core[0];
        f2[1] = bey*S + f2Core[1];
        f2[2] = bez*S + f2Core[2];

        f1[0] = -f2[0];
        f1[1] = -f2[1];
        f1[2] = -f2[2];
        return;
}
