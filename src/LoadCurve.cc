/**************************************************************************
 *
 *      Author:  Moono Rhee
 *      Function: LoadCurve
 *
 *      Description: This subroutine defines the type of load curves.
 *                   Works only with the conventional x-y-z (global)
 *                   coordinate system.  If loading axis rotated, the
 *                   loading axis can be rotated, or One can rotate
 *                   the initial dislocation configuration to a
 *                   "laboratory" coordinate system.
 *
 *                   Types of load curves:
 *                      0  Creep
 *                      1  Constant strain test
 *                      2  Displacement-controlled
 *                      3  Junction unzipping jump test
 *                      4  Total strain controlled cyclic load
 *                      5  Plastic strain controlled cyclic load
 *                      6  Load-time curve
 *
 *      Last Modified:  01/03/2001 - original version
 *                      03/13/2003 - M. Rhee Removed anisotropic elastic
 *                                   constants.  Modified to include
 *                                   isotropic Hooke's law for arbitray
 *                                   loading.
 *                      11/11/2003 - MasatoH Implementation of loading axis
 *                                   rotation due to accumuration of
 *                                   material spin.  Instead of crystal
 *                                   system, lab frame is rotated in opposite
 *                                   way
 *                      06/23/2004 - M.Rhee Added strain decomposition and
 *                                   density flux decompostion.  Modified
 *                                   message passing calls for all decomposed
 *                                   strain/density info
 *                      07/12/2004 - Masato Strain contolled cyclic load
 *                                   is implemented.
 *                      01/07/2020 - N. Bertin added rotation to the applied 
 *                                   stress
 *
 ***************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include <stdio.h>
#include <math.h>


/*
 *      Ss(): Sine for small angle i.e. Ss~x
 *      Cs(): Cosine for small angle i.e. Cs~1-x^2/2
 */
#define Ss(a) ((a))
#define Cs(a) (1.0 - (0.5 * (a)*(a)))


/*
 *      Function:     SpinMatrix
 *      Description:  Small rotation matrix for accumulaed
 *                    rotations around axis 1, 2 and 3.
 *
 *                    Cs() = Cosine for small angle
 *                    Ss() = Sin for small angle
 */
static void SpinMatrix(real8 p1, real8 p2, real8 p3, real8 Rspin[3][3])
{
        Rspin[0][0] =  Cs(p3)*Cs(p2);
        Rspin[1][1] =  Cs(p3)*Cs(p1) + Ss(p3)*Ss(p2)*Ss(p1);
        Rspin[2][2] =  Cs(p2)*Cs(p1);
        Rspin[0][1] = -Ss(p3)*Cs(p1) + Cs(p3)*Ss(p1)*Ss(p2);
        Rspin[1][2] = -Cs(p3)*Ss(p1) + Ss(p3)*Ss(p2)*Cs(p1);
        Rspin[2][0] = -Ss(p2);
        Rspin[0][2] =  Ss(p3)*Ss(p1) + Cs(p3)*Cs(p1)*Ss(p2);
        Rspin[1][0] =  Ss(p3)*Cs(p2);
        Rspin[2][1] =  Cs(p2)*Ss(p1);

        return;
}


static real8 CubicYoungModulus(real8 C11, real8 C12, real8 C44, real8 hkl[3])
{
        real8 sd  = C11*C11+C11*C12-2.0*C12*C12;
        real8 s11 = (C11+C12)/sd;
        real8 s12 = -C12/sd;
        real8 s44 = 1.0/C44;
        real8 s0  = s11-s12-0.5*s44;

        real8 h2 = hkl[0]*hkl[0];
        real8 k2 = hkl[1]*hkl[1];
        real8 l2 = hkl[2]*hkl[2];
        real8 fact = (h2*k2+k2*l2+l2*h2)/(h2+k2+l2)/(h2+k2+l2);

        return 1.0/(s11-2.0*s0*fact);
}


void LoadCurve(Home_t *home, real8 deltaStress[3][3])
{
        int     i, j, loadtype;
        int     numLoadCycle, numLoadCycle2;
        real8   youngs, erate, dtt;
        real8   shr,pois;
        real8   modulus, dpl_stn, dStress, amag, al, am, an;
        real8   pstnijk, eAmp, timeNow, cTimeOld;
        real8   dCyclicStrain;
        real8   totCyclicStrain, netCyclicStrain;
        Param_t *param;


        TimerStart(home, LOADCURVE);

        param     = home->param;
        loadtype  = param->loadType;
        shr       = param->shearModulus;
        pois      = param->pois;
        youngs    = 2.0 * shr * (1.0+pois);
        modulus   = youngs;

        erate     = param->eRate;
        dtt       = param->realdt;

        totCyclicStrain = 0.0;

/*
 *      for cyclic load
 */
        eAmp            = param->eAmp;
        timeNow         = param->timeNow;
        cTimeOld        = param->cTimeOld;
        numLoadCycle    = param->numLoadCycle;
        netCyclicStrain = param->netCyclicStrain;
        dCyclicStrain   = param->dCyclicStrain;

        deltaStress[0][0] = param->appliedStress[0];
        deltaStress[1][1] = param->appliedStress[1];
        deltaStress[2][2] = param->appliedStress[2];
        deltaStress[1][2] = param->appliedStress[3];
        deltaStress[2][0] = param->appliedStress[4];
        deltaStress[0][1] = param->appliedStress[5];

#if debug
        printf("cosinesmall in LC %e\n",Cs(0.001));
        SpinMatrix(0.1,0.1,0.1,Rspin);

        gnu_plot(home);

        if (home->cycle == 2) {
            Fatal("Stopping at DeltaPlasticStrain to debug");
        }
#endif

/*
 *      If we're including osmotic forces on dislocation segments, we
 *      need to use the delta plastic strain to adjust the vacancy
 *      concentration.  Basically, sum the diagonal of the delta plastic
 *      strain and add to the vacancy concentration.
 */
        if (param->vacancyConcEquilibrium > 0.0) {
            param->vacancyConc += (param->delpStrain[0] +
                                   param->delpStrain[1] +
                                   param->delpStrain[2]);
        }

        for (i = 0; i < 6; i++) {
            param->totpStn[i] += param->delpStrain[i];
            param->totpSpn[i] += param->delpSpin[i];
            param->totedgepStrain[i] += param->dedgepStrain[i];
            param->totscrewpStrain[i] += param->dscrewpStrain[i];
        }


/*
 *      Part for Loading Axis Rotation due to Small Deformation Spin.
 *
 *      Some changes in rotation angles due to the deformation spins
 *      around x, y, and z axis
 */
        if (param->crystalRotation) {

            real8   phi1, phi2, phi3;
            real8   Rspin[3][3];
            real8   tempedot[3], temppassedot[3];
            real8   tmpS[3][3], tmpRS[3][3], RspinT[3][3];

            phi1 = -param->delpSpin[3];
            phi2 =  param->delpSpin[4];
            phi3 = -param->delpSpin[5];

/*
 *          Matrix for (combined) rotation around x,y,and z in the sequence.
 *          This sequential rotation is correct only for small changes
 *          in the angles since real material rotation occurs simultaneously.
 *          For counter-rotation, sign of phi is flipped.
 */
            SpinMatrix(phi1, phi2, phi3, Rspin);

/*
 *          Compute nodal velocity : 
 */
            tempedot[0] = param->edotdir[0];
            tempedot[1] = param->edotdir[1];
            tempedot[2] = param->edotdir[2];

            temppassedot[0] = 0.0;
            temppassedot[1] = 0.0;
            temppassedot[2] = 0.0;

            Matrix33Vector3Multiply(Rspin, tempedot, temppassedot);

            param->edotdir[0] = temppassedot[0];
            param->edotdir[1] = temppassedot[1];
            param->edotdir[2] = temppassedot[2];

/*
 *          Rotate the stress due to small spin
 */
            tmpS[0][0] = param->appliedStress[0];
            tmpS[1][1] = param->appliedStress[1];
            tmpS[2][2] = param->appliedStress[2];
            tmpS[1][2] = param->appliedStress[3];
            tmpS[2][0] = param->appliedStress[4];
            tmpS[0][1] = param->appliedStress[5];
            tmpS[2][1] = tmpS[1][2];
            tmpS[0][2] = tmpS[2][0];
            tmpS[1][0] = tmpS[0][1];

            Matrix33_Transpose(RspinT,Rspin);    // RspinT = Transpose(Rspin)

            Matrix33_Mul(tmpRS,Rspin, tmpS  );   // tmpRS = Rspin * tmpS
            Matrix33_Mul(tmpS ,tmpRS, RspinT);   // tmpS  = Rspin * tmpS * RspinT

            param->appliedStress[0] = tmpS[0][0];
            param->appliedStress[1] = tmpS[1][1];
            param->appliedStress[2] = tmpS[2][2];
            param->appliedStress[3] = tmpS[1][2];
            param->appliedStress[4] = tmpS[2][0];
            param->appliedStress[5] = tmpS[0][1];

/*
 *          Rotate the plastic strain due to small spin
 */
    	    tmpS[0][0] = param->totpStn[0];
            tmpS[1][1] = param->totpStn[1];
            tmpS[2][2] = param->totpStn[2];
            tmpS[1][2] = param->totpStn[3];
            tmpS[2][0] = param->totpStn[4];
            tmpS[0][1] = param->totpStn[5];
            tmpS[2][1] = tmpS[1][2];
            tmpS[0][2] = tmpS[2][0];
            tmpS[1][0] = tmpS[0][1];

            Matrix33_Transpose(RspinT,Rspin);    // RspinT = Transpose(Rspin)
	    
            Matrix33_Mul(tmpRS,Rspin, tmpS  );   // tmpRS = Rspin * tmpS
            Matrix33_Mul(tmpS ,tmpRS, RspinT);   // tmpS  = Rspin * tmpS * RspinT
	    
    	    param->totpStn[0] = tmpS[0][0];
            param->totpStn[1] = tmpS[1][1];
            param->totpStn[2] = tmpS[2][2];
            param->totpStn[3] = tmpS[1][2];
            param->totpStn[4] = tmpS[2][0];
            param->totpStn[5] = tmpS[0][1];

/*
 *          Rotate the edge plastic strain due to small spin
 */
    	    tmpS[0][0] = param->totedgepStrain[0];
            tmpS[1][1] = param->totedgepStrain[1];
            tmpS[2][2] = param->totedgepStrain[2];
            tmpS[1][2] = param->totedgepStrain[3];
            tmpS[2][0] = param->totedgepStrain[4];
            tmpS[0][1] = param->totedgepStrain[5];
            tmpS[2][1] = tmpS[1][2];
            tmpS[0][2] = tmpS[2][0];
            tmpS[1][0] = tmpS[0][1];

            Matrix33_Transpose(RspinT,Rspin);    // RspinT = Transpose(Rspin)
	    
            Matrix33_Mul(tmpRS,Rspin, tmpS  );   // tmpRS = Rspin * tmpS
            Matrix33_Mul(tmpS ,tmpRS, RspinT);   // tmpS  = Rspin * tmpS * RspinT
	    
    	    param->totedgepStrain[0] = tmpS[0][0];
            param->totedgepStrain[1] = tmpS[1][1];
            param->totedgepStrain[2] = tmpS[2][2];
            param->totedgepStrain[3] = tmpS[1][2];
            param->totedgepStrain[4] = tmpS[2][0];
            param->totedgepStrain[5] = tmpS[0][1];

/*
 *          Rotate the screw plastic strain due to small spin
 */
    	    tmpS[0][0] = param->totscrewpStrain[0];
            tmpS[1][1] = param->totscrewpStrain[1];
            tmpS[2][2] = param->totscrewpStrain[2];
            tmpS[1][2] = param->totscrewpStrain[3];
            tmpS[2][0] = param->totscrewpStrain[4];
            tmpS[0][1] = param->totscrewpStrain[5];
            tmpS[2][1] = tmpS[1][2];
            tmpS[0][2] = tmpS[2][0];
            tmpS[1][0] = tmpS[0][1];

            Matrix33_Transpose(RspinT,Rspin);    // RspinT = Transpose(Rspin)
	    
            Matrix33_Mul(tmpRS,Rspin,tmpS  );    // tmpRS = Rspin * tmpS
            Matrix33_Mul(tmpS ,tmpRS,RspinT);    // tmpS  = Rspin * tmpS * RspinT 
	    
    	    param->totscrewpStrain[0] = tmpS[0][0];
            param->totscrewpStrain[1] = tmpS[1][1];
            param->totscrewpStrain[2] = tmpS[2][2];
            param->totscrewpStrain[3] = tmpS[1][2];
            param->totscrewpStrain[4] = tmpS[2][0];
            param->totscrewpStrain[5] = tmpS[0][1];
	   
        }

/*
 *      Arbitrary loading direction but keep in lab frame
 */
        al = param->edotdir[0];
        am = param->edotdir[1];
        an = param->edotdir[2];

        amag = sqrt(al*al + am*am + an*an);
        amag = ( (fabs(amag)>0.0) ? (1.0/amag) : 1.0 );    // (avoid potential divide by null edotdir)

        al *= amag;
        am *= amag;
        an *= amag;

        switch(loadtype) {
/*
 *          creep - what we have been using so this should be default
 */
            case 0:
                break;
/*
 *          constant strain rate
 */
            case 1:

                dpl_stn=      param->delpStrain[0]*al*al +
                              param->delpStrain[1]*am*am +
                              param->delpStrain[2]*an*an +
                          2.0*param->delpStrain[3]*am*an +
                          2.0*param->delpStrain[4]*an*al +
                          2.0*param->delpStrain[5]*al*am;

#ifdef ANISOTROPIC

/*                In anisotropic elasticity, only in the case of cubic systems
 *                (BCC, FCC) can we define a modulus from C11, C12, C44. We cannot define 
 *                a modulus for HCP or rhombohedral, so use Young's modulus calculated from
 *                isotropic elasticity.
 */

                switch(home->param->materialType) 
                {
                case MAT_TYPE_BCC:
                case MAT_TYPE_FCC:
                   param->cubicModulus = 1;
                   param->cubicC11 = param->C11;
                   param->cubicC12 = param->C12;
                   param->cubicC44 = param->C44;
                   break;
                case MAT_TYPE_HCP:
                case MAT_TYPE_RHOMBOHEDRAL_VA:
                   param->cubicModulus = 0;
                   break;
                }
#endif

                if (param->cubicModulus) {
/*
 *                  Compute the stress/strain modulus for cubic crystals
 *                  as a function of the current loading direction
 */
                    modulus = CubicYoungModulus(param->cubicC11, param->cubicC12,
                                                param->cubicC44, param->edotdir);
                    param->YoungsModulus = modulus;

                } else {
                    modulus = youngs;
                }

/*
 *              local in the [l m n] frame
 */
                dStress = modulus * (erate*dtt - dpl_stn);
/*
 *              global (100)-(010)-(001) frame
 */
                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am;
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;

                param->totstraintensor[0] = erate * param->timeNow *  al*al;
                param->totstraintensor[1] = erate * param->timeNow *  am*am;
                param->totstraintensor[2] = erate * param->timeNow *  an*an;
                param->totstraintensor[3] = erate * param->timeNow *  an*am;
                param->totstraintensor[4] = erate * param->timeNow *  an*al;
                param->totstraintensor[5] = erate * param->timeNow *  al*am;
                break;

/*
 *          strain control cyclic load
 *
 *              stainCycle    = current loading cycle
 *              eAmp          = strain amplitude for each side
 *              dCyclicStrain = change in the strain for each side
 *              acumStrain    = accumulated strain
 *              sgnLoad       = sign of load
 */
            case 4:

/*
 *              Cover for specific loading direction also
 */
                dpl_stn =      param->delpStrain[0]*al*al +
                               param->delpStrain[1]*am*am +
                               param->delpStrain[2]*an*an +
                           2.0*param->delpStrain[3]*am*an +
                           2.0*param->delpStrain[4]*an*al +
                           2.0*param->delpStrain[5]*al*am;

                if (param->cubicModulus) {
/*
 *                  Compute the stress/strain modulus for cubic crystals
 *                  as a function of the current loading direction
 */
                    modulus = CubicYoungModulus(param->cubicC11, param->cubicC12,
                                                param->cubicC44, param->edotdir);
                    param->YoungsModulus = modulus;

                } else {
                    modulus = youngs;
                }

                dCyclicStrain = erate*dtt;
                param->dCyclicStrain = dCyclicStrain;

/*
 *              local in the [l m n] frame
 */
                dStress= modulus * (dCyclicStrain - dpl_stn);

                totCyclicStrain = fabs(erate*timeNow);
                numLoadCycle    = (int)rint(0.5*totCyclicStrain/eAmp);
                numLoadCycle2   = (int)rint(0.5*totCyclicStrain/eAmp-0.5);

                netCyclicStrain = fmod(totCyclicStrain, 2*eAmp);

                if (fabs(netCyclicStrain) > eAmp) {
                    netCyclicStrain = 2*eAmp - netCyclicStrain;
                }

                netCyclicStrain = pow(-1, numLoadCycle2) *
                                  fabs(netCyclicStrain);

                param->netCyclicStrain = netCyclicStrain;

#if 0
                printf("loading cycle 1&2 %d %d\n", numLoadCycle,
                       numLoadCycle2);
                printf("net strain %e\n", netCyclicStrain);
                printf("loading cycle %d\n", numLoadCycle);
                printf("Load Curve: dtt,totSt %e %e \n", dtt, totCyclicStrain);
                printf("Load Curve: dtt,dttOld dCyclic %e %e %e %e\n",
                       timeNow, cTimeOld, dCyclicStrain, eAmp);
#endif

                cTimeOld = timeNow;
                erate = fabs(erate)*pow(-1,numLoadCycle);
                dCyclicStrain = 0;
                param->cTimeOld = cTimeOld;
                param->numLoadCycle = numLoadCycle;
                param->eRate = erate;

/*
 *              global (100)-(010)-(001) frame
 */
                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am;
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;

                param->totstraintensor[0] = netCyclicStrain *  al*al;

                break;


/*
 *             Plastic strain control cyclic load
 *                 stainCycle    = current loading cycle
 *                 eAmp          = strain amplitude for each side
 *                 dCyclicStrain = change in the strain for each side
 *                 acumStrain    = accumulated strain
 *                 sgnLoad       = sign of load
 */
            case 5:

/*
 *              Cover for specific loading direction also
 */
                pstnijk = param->totpStn[0]*al*al +
                          param->totpStn[1]*am*am +
                          param->totpStn[2]*an*an +
                          2.0*param->totpStn[3]*am*an +
                          2.0*param->totpStn[4]*an*al +
                          2.0*param->totpStn[5]*al*am;

                dpl_stn = param->delpStrain[0]*al*al +
                          param->delpStrain[1]*am*am +
                          param->delpStrain[2]*an*an +
                          2.0*param->delpStrain[3]*am*an +
                          2.0*param->delpStrain[4]*an*al +
                          2.0*param->delpStrain[5]*al*am;

                if (param->cubicModulus) {
/*
 *                  Compute the stress/strain modulus for cubic crystals
 *                  as a function of the current loading direction
 */
                    modulus = CubicYoungModulus(param->cubicC11, param->cubicC12,
                                                param->cubicC44, param->edotdir);
                    param->YoungsModulus = modulus;

                } else {
                    modulus = youngs;
                }

                dCyclicStrain = erate*dtt;
                param->dCyclicStrain = dCyclicStrain;

/*
 *              local in the [l m n] frame
 */
                dStress= modulus * (dCyclicStrain - dpl_stn);

                totCyclicStrain += fabs(dpl_stn);
                numLoadCycle    = (int)rint(0.5*pstnijk/eAmp);
                numLoadCycle2   = (int)rint(0.5*pstnijk/eAmp-0.5);

                netCyclicStrain = fmod(pstnijk, 2*eAmp);

                if (fabs(netCyclicStrain) > eAmp ) {
                    netCyclicStrain = 2*eAmp - netCyclicStrain;
                }

                netCyclicStrain = pow(-1, numLoadCycle2) *
                                  fabs(netCyclicStrain);

                param->netCyclicStrain = netCyclicStrain;

#if 0
                printf("loading cycle 1&2 %d %d\n", numLoadCycle,
                       numLoadCycle2);
                printf("net strain %e\n", netCyclicStrain);
                printf("loading cycle %d\n", numLoadCycle);
                printf("Load Curve: dtt,totSt %e %e \n", dtt, totCyclicStrain);
                printf("Load Curve: dtt,dttOld dCyclic %e %e %e %e\n",
                       timeNow, cTimeOld, dCyclicStrain, eAmp);
#endif
                cTimeOld = timeNow;
                erate = fabs(erate) * pow(-1,numLoadCycle);
                dCyclicStrain = 0;
                param->cTimeOld = cTimeOld;
                param->numLoadCycle = numLoadCycle;
                param->eRate = erate;

/*
 *              global (100)-(010)-(001) frame
 */
                param->appliedStress[0] += dStress *al*al;
                param->appliedStress[1] += dStress *am*am;
                param->appliedStress[2] += dStress *an*an;
                param->appliedStress[3] += dStress *an*am;
                param->appliedStress[4] += dStress *an*al;
                param->appliedStress[5] += dStress *al*am;

                param->totstraintensor[0] = netCyclicStrain *  al*al;
                param->totstraintensor[1] = netCyclicStrain *  am*am;
                param->totstraintensor[2] = netCyclicStrain *  an*an;
                param->totstraintensor[3] = netCyclicStrain *  an*am;
                param->totstraintensor[4] = netCyclicStrain *  an*al;
                param->totstraintensor[5] = netCyclicStrain *  al*am;

                break;

/*
 *          User defined load-time curve
 */
            case 6:
                break;

            default:
                Fatal("Load curves not defined. Stopping the program. \n");
            break;

        }  /* end: switch(loadtype) */

        TimerStop(home, LOADCURVE);

        deltaStress[0][0] = param->appliedStress[0] - deltaStress[0][0];
        deltaStress[1][1] = param->appliedStress[1] - deltaStress[1][1];
        deltaStress[2][2] = param->appliedStress[2] - deltaStress[2][2];
        deltaStress[1][2] = param->appliedStress[3] - deltaStress[1][2];
        deltaStress[2][0] = param->appliedStress[4] - deltaStress[2][0];
        deltaStress[0][1] = param->appliedStress[5] - deltaStress[0][1];
        deltaStress[2][1] = deltaStress[1][2];
        deltaStress[0][2] = deltaStress[2][0];
        deltaStress[1][0] = deltaStress[0][1];

/*
 *      Update young's modulus in Param_t so the correct value is written
 *      into the control file.  Value is written to control file for
 *      informational purposes only.
 */
        param->YoungsModulus = modulus;
        return;
}
