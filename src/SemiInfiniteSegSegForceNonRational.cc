/**************************************************************************
 *
 *      Module:      SemiInfiniteSegSegForceNonRational.c
 *
 *      Authors:     Tom Arsenlis, Meijie Tang
 *
 *      Description: This module contains the functions needed for
 *                   calculating forces from interactions between a
 *                   non-rational finite length dislocation and a
 *                   semi-infinite dislocation.  This is only applicable
 *                   when the FEM code is hooked into ParaDiS for doing
 *                   simulations with free surfaces conditions.
 *
 *                   In SemiInfiniteSegSegForceNonRational(),
 *                   SpecialSemiInfiniteSegSegForceNonRational() and
 *                   SpecialSemiInfiniteSegSegForceHalfNonRationsl()
 *                   functions there is a divergence of the stress
 *                   field on an element by element basis that cancels
 *                   when the total segment configuration closes on itself
 *                   or is bounded by semi-infinite segments.
 *
 *                   Versions of these function have been modified from
 *                   these original versions to include a correction that
 *                   makes the stress field due to a dislocation segment
 *                   divergence free on an element by element basis so
 *                   the system does not have to be closed or terminated
 *                   by semi-infinite segments in order to obtain the
 *                   divergence-free stress field.  The modified versions
 *                   are contained in SemiInfiniteSegSegForceRational.c
 *                   and use the naming convention <funcName>Rational()
 *                   rather than <funcName>NonRational().
 *
 *                   See the User's Guide for more details on the method
 *                   method used to calculate the correction value.
 *
 *
 *      Includes public functions:
 *              SemiInfiniteSegSegForceNonRational()
 *
 *      Includes private functions:
 *              SpecialSemiInfiniteSegSegForceNonRational()
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"


#ifdef _ARLFEM
/*-------------------------------------------------------------------------
 *
 *      Function:     SpecialSemiInfiniteSegSegForceNonRational
 *      Description:  Special function for calculating forces between
 *                    dislocation segments too close to parallel to be
 *                    calculated via the function used for regular
 *                    segment/segment forces.
 *      Arguments:
 *          p1*          Coordinates of the surface node which is the
 *                       endpoint of the semi-infinite dislocation segment.
 *          p2*          Coordinates of a point which is outside the
 *                       free surface and on the semi-infinite dislocation
 *                       segment.
 *          p3*,p4*      endpoints for dislocation segment beginning
 *                       at point p3 and ending at point p4
 *          bpx,bpy,bpz  burgers vector for segment p1->p2
 *          bx,by,bz     burgers vector for segment p3->p4
 *          a            core parameter
 *          MU           shear modulus
 *          NU           poisson ratio
 *          fp1*         pointers to locations in which to return force
 *                       on the endpoint of the semi-infinite segment.
 *                       WARNING! Currently force at this point is
 *                       always zeroed out!
 *          fp3*, fp4*   pointers to arrays in which to return the
 *                       forces on nodes p3 and p4 respectively from
 *                       the semi-infinite segment.
 *
 *-----------------------------------------------------------------------*/
static void SpecialSemiInfiniteSegSegForceNonRational(Home_t *home,
                    real8 p1x, real8 p1y, real8 p1z,
                    real8 p2x, real8 p2y, real8 p2z,
                    real8 p3x, real8 p3y, real8 p3z,
                    real8 p4x, real8 p4y, real8 p4z,
                    real8 bpx, real8 bpy, real8 bpz,
                    real8 bx, real8 by, real8 bz,
                    real8 a, real8 MU, real8 NU, real8 ecrit,
                    real8 *fp1x, real8 *fp1y, real8 *fp1z,
                    real8 *fp3x, real8 *fp3y, real8 *fp3z,
                    real8 *fp4x, real8 *fp4y, real8 *fp4z)
{
        int i, j , alt1[3]={1,2,0}, alt2[3]={2,0,1};
        real8 eps, c, a2, d2, a2_d2, a2d2inv, flip;
        real8 x1[3], x2[3], x3[3], x4[3], b[3], bp[3];
        real8 f3[3], f4[3];
        real8 vec1[3], vec2[3], t[3], nd[3], dir[3];
        real8 temp1, temp2;
        real8 R[3], Rdt, x2mod[3], x2a[3], x2b[3];
        real8 oneoverL, oneoverLp;
        real8 y[2], z[2], yv[2], zv[2], ypz[2], ymz[2];
        real8 Ra[2], Rainv[2], Log_Ra_ypz[2];
        real8 temp, tmp[8];
        real8 common1[2], common2[3], common3[3];
        real8 magdiff, magdiff2;
        real8 wx, wy, wz;
        real8 fp3xcor, fp3ycor, fp3zcor;
        real8 fp4xcor, fp4ycor, fp4zcor;
        real8 f_003v[2], f_103v[2], f_113v[2], f_213v[2];
        real8 f_005v[2], f_105v[2], f_115v[2], f_215v[2];
        real8 f_003vinf[2], f_103vinf[2], f_113vinf[2], f_213vinf[2];
        real8 f_005vinf[2], f_105vinf[2], f_115vinf[2], f_215vinf[2];
        real8 f_003, f_103, f_113, f_213;
        real8 f_005, f_105, f_115, f_215;
        real8 Fint_003, Fint_113, Fint_005, Fint_115;
        real8 I_003[3], I_113[3], I_005[3], I_115[3];
        real8 m4p, m8p, m4pn, a2m4pn, a2m8p;
        real8 tdb, tdbp, nddb, bpctdb, bpctdnd;
        real8 bct[3], bpct[3], ndct[3], bpctct[3];
        real8 tanthetac;
        real8 tp[3]; 
        real8 shape,shapea,shapeb,shapemid; 


        tanthetac = sqrt((ecrit*1.01)/(1 - ecrit*1.01));

        eps    = 1e-12;
        a2     = a*a;
        m4p    = 0.25 * MU / M_PI;
        m8p    = 0.5 * m4p;
        m4pn   = m4p / ( 1 - NU );
        a2m4pn = a2 * m4pn;
        a2m8p  = a2 * m8p;
            
        *fp1x = 0.0;
        *fp1y = 0.0;
        *fp1z = 0.0;
            
        *fp3x = 0.0;
        *fp3y = 0.0;
        *fp3z = 0.0;
            
        *fp4x = 0.0;
        *fp4y = 0.0;
        *fp4z = 0.0;
        
        x1[0]=p1x;
        x1[1]=p1y;
        x1[2]=p1z;
        x2[0]=p2x;
        x2[1]=p2y;
        x2[2]=p2z;
        x3[0]=p3x;
        x3[1]=p3y;
        x3[2]=p3z;
        x4[0]=p4x;
        x4[1]=p4y;
        x4[2]=p4z;
        
        b[0]=bx;
        b[1]=by;
        b[2]=bz;
        bp[0]=bpx;
        bp[1]=bpy;
        bp[2]=bpz;
        
        for(i=0;i<3;i++) { 
            vec1[i]=x4[i]-x3[i];
            vec2[i]=x2[i]-x1[i];
        }

        temp1=0.0e0;
        temp2=0.0e0;

        for(i=0;i<3;i++) { 
            temp1+=vec1[i]*vec1[i];
            temp2+=vec2[i]*vec2[i];
        }

        oneoverL =1/sqrt(temp1);
        oneoverLp=1/sqrt(temp2);

        for(i=0;i<3;i++) { 
            t[i]=vec1[i]*oneoverL;
            tp[i]=vec2[i]*oneoverLp;
        }
        
        c=0.0e0;

        for(i=0;i<3;i++) { 
            c+=t[i]*tp[i];
        }

        flip=1.0e0;

        if (c < 0) {
            flip=-1.0e0;
            for(i=0;i<3;i++) { 
                bp[i]=-bp[i];
                tp[i]=-tp[i];
            }         
        } 
             
/*
 *      Find f3 and f4
 */
        temp1=0.0e0;

        for (i=0;i<3;i++) {
            temp1+=vec2[i]*t[i];
        }
            
        for (i=0;i<3;i++) {
            x2mod[i]=x1[i]+temp1*t[i];
        }
        
        for (i=0;i<3;i++) {
            vec2[i]=x2[i]-x2mod[i];
        }
        
        temp=0.0e0;

        for (i=0;i<3;i++) {
            temp+=vec2[i]*vec2[i];
        }
            
        magdiff=sqrt(temp);
        
        if (magdiff*oneoverLp >eps) { 
            for (i=0;i<3;i++) {
                dir[i]=vec2[i]/magdiff;
            }
        } else {
            for (i=0;i<3;i++) {
                vec2[i]=-t[0]*t[i];
            }

            vec2[0]+=1.0e0;
            magdiff2 = sqrt(DotProduct(vec2, vec2));

/*
 *          If segments are parallel in the [1 0 0] direction, change
 *          to the [0 1 0] direction so we don't seg fault.
 */
            if (magdiff2 < eps) {
                for (i = 0; i < 3; i++) {
                    vec2[i] = -t[0] * t[i];
                }

                vec2[1] += 1.0;
                magdiff2 = sqrt(DotProduct(vec2, vec2));
            }

            for (i=0;i<3;i++) {
                dir[i]=vec2[i]/magdiff2;
            }
        }

        shape=0.5e0*magdiff/(fabs(temp1)*tanthetac);

        for (i=0;i<3;i++) {
            temp=fabs(temp1)*tanthetac*dir[i];
            x2a[i]=x2mod[i]+temp;
            x2b[i]=x2mod[i]-temp;
        }

        
        for (i=0;i<3;i++) {
            R[i] = x3[i] - x1[i];
        }
        
        Rdt=0.0e0;

        for (i=0;i<3;i++) {
            Rdt+=R[i]*t[i];
        }
        
        for (i=0;i<3;i++) {
            nd[i]=R[i]-Rdt*t[i];
        }
        
        d2=0.0e0;

        for (i=0;i<3;i++) {
            d2+=nd[i]*nd[i];
        }    
        
        for (j=0;j<2;j++) {
            y[j]=0.0e0;
            z[j]=0.0e0;
        }  
        
        for (i=0;i<3;i++) {
            y[0]+=x3[i]*t[i];
            y[1]+=x4[i]*t[i];
            z[0]+=-x1[i]*t[i];
        } 

        yv[0]=y[0];
        yv[1]=y[1];
        zv[0]=z[0];
        zv[1]=z[0];
            
        a2_d2 = a2 + d2;   
        
        for (j=0;j<2;j++) {
            ypz[j] = yv[j] + zv[j];
            ymz[j] = yv[j] - zv[j];
        }
            
        for (j=0;j<2;j++) {
            tmp[j]=a2_d2 + ypz[j]*ypz[j];
        }
            
        for (j=0;j<2;j++) {
            Ra[j]=sqrt(tmp[j]);
        }
            
        for (j=0;j<2;j++) {
            Rainv[j]=1.0e0/Ra[j];
        }

        a2d2inv = 1.0e0 / a2_d2;
        
        for (j=0;j<2;j++) {
            tmp[j]=Ra[j] + ypz[j];
        }
            
        for (j=0;j<2;j++) {
            Log_Ra_ypz[j]=log(tmp[j]);
        }
        
        for (j=0;j<2;j++) {
            common1[j] = ymz[j] * Ra[j] * a2d2inv;
            f_115v[j] = -a2d2inv * ypz[j] * Rainv[j];
            f_115vinf[j]=0.0e0;
        }
        
        temp=2.0e0*a2d2inv;
        
        for (j=0;j<2;j++) {
            f_003v[j]    = Ra[j] * a2d2inv;
            f_003vinf[j] = yv[j] * a2d2inv;

            f_103v[j]    = Log_Ra_ypz[j] - common1[j];
            f_103vinf[j] = -yv[j] * yv[j] * a2d2inv;

            f_103v[j]    = -f_103v[j] / 2.0;
            f_103vinf[j] = -f_103vinf[j] / 2.0;

            f_113v[j]    = -Log_Ra_ypz[j];
            f_113vinf[j] = 0.0e0;

            f_213v[j]    = zv[j]*Log_Ra_ypz[j] - Ra[j];
            f_213vinf[j] =0.0e0;

            f_005v[j]    = a2d2inv * (2.0 * a2d2inv * Ra[j] - Rainv[j]);
            f_005vinf[j] = temp * yv[j] * a2d2inv;

            f_105v[j]    = (common1[j] - yv[j] * Rainv[j]) * a2d2inv;
            f_105vinf[j] = 0.5 * yv[j] * f_005vinf[j];

            f_215v[j]    =  Rainv[j] - zv[j] * f_115v[j];
            f_215vinf[j] = 0.0e0;
        }
        
        f_003 = flip*(f_003v[0]-f_003v[1])+(f_003vinf[0]-f_003vinf[1]);
        f_103 = flip*(f_103v[0]-f_103v[1])+(f_103vinf[0]-f_103vinf[1]);
        f_113 = flip*(f_113v[0]-f_113v[1])+(f_113vinf[0]-f_113vinf[1]);
        f_213 = flip*(f_213v[0]-f_213v[1])+(f_213vinf[0]-f_213vinf[1]);
        f_005 = flip*(f_005v[0]-f_005v[1])+(f_005vinf[0]-f_005vinf[1]);
        f_105 = flip*(f_105v[0]-f_105v[1])+(f_105vinf[0]-f_105vinf[1]);
        f_115 = flip*(f_115v[0]-f_115v[1])+(f_115vinf[0]-f_115vinf[1]);
        f_215 = flip*(f_215v[0]-f_215v[1])+(f_215vinf[0]-f_215vinf[1]);
        

        for (i=0;i<3;i++) {
            bct[i]=b[alt1[i]]*t[alt2[i]] - b[alt2[i]]*t[alt1[i]];
            bpct[i]=bp[alt1[i]]*t[alt2[i]] - bp[alt2[i]]*t[alt1[i]];
            ndct[i]=nd[alt1[i]]*t[alt2[i]] - nd[alt2[i]]*t[alt1[i]];
        }
        
        tdb=0.0e0;
        tdbp=0.0e0;
        nddb=0.0e0;
        bpctdb=0.0e0;
        bpctdnd=0.0e0;

        for (i=0;i<3;i++) {
            tdb += t[i]*b[i];
            tdbp+= t[i]*bp[i];
            nddb+= nd[i]*b[i];
            bpctdb += bpct[i]*b[i];
            bpctdnd += bpct[i]*nd[i];
            
        }
            
        temp = tdb*tdbp;
            
        for (i=0;i<3;i++) {
            bpctct[i] = tdbp*t[i] - bp[i];
            common2[i] = temp*nd[i];
            common3[i] = bpctdnd*bct[i];
        }   

        tmp[0]=(m4pn-m4p)*tdb;
        tmp[1]=m4pn*bpctdnd*nddb;
        tmp[2]=a2m8p*tdb;
        tmp[3]=m4pn*bpctdnd*tdb;
        
        for (i=0;i<3;i++) {
            I_003[i] = m4pn*(nddb*bpctct[i] + bpctdb*ndct[i] - common3[i]) -
                       m4p*common2[i]; 
            I_113[i] =  tmp[0]*bpctct[i];
            I_005[i] = -a2m8p*common2[i] - a2m4pn*common3[i] - tmp[1]*ndct[i];
            I_115[i] = -tmp[2]*bpctct[i] - tmp[3]*ndct[i];
        }
                     
        Fint_003 = f_103 - y[0]*f_003;
        Fint_113 = f_213 - y[0]*f_113;
        Fint_005 = f_105 - y[0]*f_005;
        Fint_115 = f_215 - y[0]*f_115;
        
        for (i=0;i<3;i++) {
            f4[i] = (I_003[i]*Fint_003 + I_113[i]*Fint_113 + I_005[i]*Fint_005 +
                     I_115[i]*Fint_115) * oneoverL;
        }

        Fint_003 = y[1]*f_003 - f_103;
        Fint_113 = y[1]*f_113 - f_213;
        Fint_005 = y[1]*f_005 - f_105;
        Fint_115 = y[1]*f_115 - f_215;

        for (i=0;i<3;i++) {
            f3[i] = (I_003[i]*Fint_003 + I_113[i]*Fint_113 + I_005[i]*Fint_005 +
                     I_115[i]*Fint_115) * oneoverL;
        }   

        shapea=shape*(1.0e0+2.0e0*shape);
        shapeb=-shape*(1.0e0-2.0e0*shape);
        shapemid=(1.0e0+2.0e0*shape)*(1.0e0-2.0e0*shape);
        *fp3x = shapemid*f3[0];
        *fp3y = shapemid*f3[1];
        *fp3z = shapemid*f3[2];
        *fp4x = shapemid*f4[0];
        *fp4y = shapemid*f4[1];
        *fp4z = shapemid*f4[2];
        
        SemiInfiniteSegSegForceNonRational(home,
                    x1[0], x1[1], x1[2], x2a[0], x2a[1], x2a[2],
                    x3[0], x3[1], x3[2], x4[0], x4[1], x4[2],
                    flip*bp[0], flip*bp[1], flip*bp[2],
                    b[0], b[1], b[2], a, MU, NU,
                    &wx, &wy, &wz,
                    &fp3xcor, &fp3ycor, &fp3zcor,
                    &fp4xcor, &fp4ycor, &fp4zcor);
        
        *fp3x += shapea*fp3xcor;
        *fp3y += shapea*fp3ycor;
        *fp3z += shapea*fp3zcor;
        *fp4x += shapea*fp4xcor;
        *fp4y += shapea*fp4ycor;
        *fp4z += shapea*fp4zcor;
        
        SemiInfiniteSegSegForceNonRational(home,
                    x1[0], x1[1], x1[2], x2b[0], x2b[1], x2b[2],
                    x3[0], x3[1], x3[2], x4[0], x4[1], x4[2],
                    flip*bp[0], flip*bp[1], flip*bp[2],
                    b[0], b[1], b[2], a, MU, NU,
                    &wx, &wy, &wz,
                    &fp3xcor, &fp3ycor, &fp3zcor,
                    &fp4xcor, &fp4ycor, &fp4zcor);
        
        *fp3x += shapeb*fp3xcor;
        *fp3y += shapeb*fp3ycor;
        *fp3z += shapeb*fp3zcor;
        *fp4x += shapeb*fp4xcor;
        *fp4y += shapeb*fp4ycor;
        *fp4z += shapeb*fp4zcor;

/*
 *      Caution: f1 is returned as zeros at this  time because it is
 *      believed that the result will never be used
 */
        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       SemiInfiniteSegSegForceNonRational
 *      Description:    Used to calculate the interaction forces between
 *                      a finite-length dislocation segment and a semi-
 *                      infinite-length dislocation segment.
 *
 *      Arguments:
 *              p1*          Coordinates of the surface node which is the
 *                           endpoint of the semi-infinite dislocation segment.
 *              p2*          Coordinates of a point which is outside the
 *                           free surface and on the semi-infinite dislocation
 *                           segment.
 *              p3*,p4*      endpoints for a finite-length dislocation segment
 *                           starting at p3x,p3y,p3z and ending at p4x,p4y,p4z
 *              bxp,byp,bzp  burgers vector for segment p1 to p2
 *              bx,by,bz     burgers vector for segment p3 to p4
 *              a            core parameter
 *              MU           shear modulus
 *              NU           poisson ratio
 *              fp1*         Location at which to return forces at the endpoint
 *                           of the semi-infinite segment.  WARNING! Currently
 *                           force at this point is always zeroed out!
 *              fp3*,fp4*    Locations in which to return forces on nodes
 *                           located at p3 and p4 respectively
 *                      
 *-----------------------------------------------------------------------*/
void SemiInfiniteSegSegForceNonRational(Home_t *home,
                 real8 p1x, real8 p1y, real8 p1z,
                 real8 p2x, real8 p2y, real8 p2z,
                 real8 p3x, real8 p3y, real8 p3z,
                 real8 p4x, real8 p4y, real8 p4z,
                 real8 bpx, real8 bpy, real8 bpz,
                 real8 bx, real8 by, real8 bz,
                 real8 a, real8 MU, real8 NU,
                 real8 *fp1x, real8 *fp1y, real8 *fp1z,
                 real8 *fp3x, real8 *fp3y, real8 *fp3z,
                 real8 *fp4x, real8 *fp4y, real8 *fp4z)
{
        real8 x1[3], x2[3], x3[3], x4[3], b[3], bp[3];
        real8 f1[3], f2[3], f3[3], f4[3];
        real8 vec1[3], vec2[3], t[3], tp[3], tctp[3];
        real8 R[2][3], tempa[2], tempb[2], y[2], z[2];
        int i, j , alt1[3]={1,2,0}, alt2[3]={2,0,1};
        real8 eps, d, c, c2, onemc2, onemc2inv, oneoverL, oneoverLp;
        real8 a2, m4p, m4pd, m8p, m8pd, m4pn, m4pnd, m4pnd2, m4pnd3;
        real8 a2m4pnd, a2m8pd, a2m4pn, a2m8p, a2_d2, a2_d2inv, denom;
        real8 temp1, temp2, temp3, temp4[4], tmp[10];
        real8 yv[2], zv[2], y2[2], z2[2], Ra[2], Rainv[2];
        real8 Ra_Rdot_tp[2], Ra_Rdot_t[2], log_Ra_Rdot_tp[2], log_Ra_Rdot_t[2];
        real8 Ra2_R_tpinv[2], Ra2_R_tinv[2];
        real8 ylog_Ra_Rdot_tp[2], zlog_Ra_Rdot_t[2];
        real8 yRa2_R_tpinv[2], zRa2_R_tinv[2];
        real8 y2Ra2_R_tpinv[2], z2Ra2_R_tinv[2];
        real8 adf_003[2], commonf223[2], commonf225[2];
        real8 commonf025[2], commonf205[2];
        real8 commonf305[2], commonf035[2], ycommonf025[2];
        real8 zcommonf205[2], zcommonf305[2];
        real8 Ra_Rdot_tpinf[2], log_Ra_Rdot_tpinf[2];
        real8 Ra2_R_tpinvinf[2], yRa2_R_tpinvinf[2];
        real8 y2Ra2_R_tpinvinf[2], y3Ra2_R_tpinvinf[2];
        real8 tf_113[2];
        real8 f_003v[2], f_103v[2], f_013v[2], f_113v[2];
        real8 f_203v[2], f_023v[2], f_005v[2], f_105v[2];
        real8 f_003,  f_103,  f_013,  f_113,  f_203,  f_023,  f_005,  f_105;
        real8 f_015v[2], f_115v[2], f_205v[2], f_025v[2];
        real8 f_215v[2], f_125v[2], f_225v[2], f_305v[2];
        real8 f_015,  f_115,  f_205,  f_025,  f_215,  f_125,  f_225,  f_305;
        real8 f_035v[2], f_315v[2], f_135v[2];
        real8 f_035,  f_315,  f_135;
        real8 tf_113inf[2], adf_003inf[2];
        real8 f_003vinf[2], f_103vinf[2], f_013vinf[2], f_113vinf[2];
        real8 f_203vinf[2], f_023vinf[2], f_005vinf[2], f_105vinf[2];
        real8 f_015vinf[2], f_115vinf[2], f_205vinf[2], f_025vinf[2];
        real8 f_215vinf[2], f_125vinf[2], f_225vinf[2], f_305vinf[2];
        real8 f_035vinf[2], f_315vinf[2], f_135vinf[2];
        real8 Fint_003, Fint_005, Fint_013, Fint_015, Fint_025, Fint_103;
        real8 Fint_105, Fint_115, Fint_125, Fint_205, Fint_215;
        real8 I_003[3], I_005[3], I_013[3], I_015[3], I_025[3], I_103[3];
        real8 I_105[3], I_115[3], I_125[3], I_205[3], I_215[3];
        real8 I00a[3], I01a[3], I10a[3], I00b[3], I01b[3], I10b[3];
        real8 bctctp[3], bct[3], bpctpct[3], bpctp[3], tcbpct[3];
        real8 bctdbp, bpctpdb, tcbpdb, tcbpdtp, tpcbdbp;
        real8 tcbp[3], tctpcbpct[3], tctpcbp[3], tctpct[3], tpct[3];
        real8 tctpcbpdb, tctpcbpdtp, tctpdb, tdb, tdbp;
        real8 tpcbctp[3], tpctcb[3], tpctctp[3], tpcb[3], tpctcbctp[3];
        real8 tpcbdt, tpctcbdbp, tpctcbdt, tpctdbp, tpdb, tpdbp;


        eps = 1e-4;            
        
        *fp1x = 0.0;
        *fp1y = 0.0;
        *fp1z = 0.0;

        *fp3x = 0.0;
        *fp3y = 0.0;
        *fp3z = 0.0;

        *fp4x = 0.0;
        *fp4y = 0.0;
        *fp4z = 0.0;

        x1[0]=p1x;
        x1[1]=p1y;
        x1[2]=p1z;
        x2[0]=p2x;
        x2[1]=p2y;
        x2[2]=p2z;
        x3[0]=p3x;
        x3[1]=p3y;
        x3[2]=p3z;
        x4[0]=p4x;
        x4[1]=p4y;
        x4[2]=p4z;
        
        b[0]=bx;
        b[1]=by;
        b[2]=bz;
        bp[0]=bpx;
        bp[1]=bpy;
        bp[2]=bpz;
        
        for(i=0;i<3;i++) { 
            vec1[i]=x4[i]-x3[i];
            vec2[i]=x2[i]-x1[i];
        }

        temp1=0.0e0;
        temp2=0.0e0;

        for(i=0;i<3;i++) { 
            temp1+=vec1[i]*vec1[i];
            temp2+=vec2[i]*vec2[i];
        }

        oneoverL =1/sqrt(temp1);
        oneoverLp=1/sqrt(temp2);
        
        for(i=0;i<3;i++) { 
            t[i]=vec1[i]*oneoverL;
            tp[i]=vec2[i]*oneoverLp;
        }
        
        c=0.0e0;

        for(i=0;i<3;i++) { 
            c+=t[i]*tp[i];
        }

        c2=c*c;
        onemc2=1-c2;
 
        if (onemc2 > eps) {
            for(i=0;i<3;i++) {
                tctp[i]=t[alt1[i]]*tp[alt2[i]]-t[alt2[i]]*tp[alt1[i]];
            }

            onemc2inv = 1/onemc2;
            
            for(i=0;i<3;i++) { 
                R[0][i]=x3[i]-x1[i];
                R[1][i]=x4[i]-x1[i];
            }

            d=0.0e0;

            for (j=0;j<2;j++) { 
                tempa[j]=0.0e0;
                tempb[j]=0.0e0;
            }

            for(i=0;i<3;i++) { 
                d+=0.5e0*(R[0][i]+R[1][i])*tctp[i];
                for (j=0;j<2;j++) { 
                    tempa[j]+=R[j][i]*t[i];
                    tempb[j]+=R[j][i]*tp[i];
                }
            }

            d*=onemc2inv;

            for (j=0;j<2;j++) { 
                y[j]=(tempa[j]-c*tempb[j])*onemc2inv;
                z[j]=(tempb[j]-c*tempa[j])*onemc2inv;
            }

/*
 *          Now we calculate the definite integrals of the force calculation
 */
            yv[0]=y[0];
            yv[1]=y[1];
            zv[0]=z[0];
            zv[1]=z[0];

            a2_d2 = a*a+d*d*onemc2;
                
            for (j=0;j<2;j++) {
                y2[j] = yv[j]*yv[j];
                z2[j] = zv[j]*zv[j];
                
            }

            for (j=0;j<2;j++) {
                temp4[j]=a2_d2 + y2[j] + z2[j] + 2.0e0*yv[j]*zv[j]*c;
            }

            temp1=onemc2*a2_d2;

            for (j=0;j<2;j++) {
                Ra[j]=sqrt(temp4[j]);
            }

            temp2=sqrt(temp1);

            for (j=0;j<2;j++) {
                Rainv[j]=1.0e0/Ra[j];
            }

            denom=1.0e0/temp2;
            a2_d2inv=1.0e0/a2_d2;

            for (j=0;j<2;j++) {
                Ra_Rdot_tp[j] = Ra[j]+(zv[j]+yv[j]*c);       
                Ra_Rdot_t[j]  = Ra[j]+(yv[j]+zv[j]*c);
                Ra_Rdot_tpinf[j] = onemc2*y2[j] + a2_d2;
            }
                
            for (j=0;j<2;j++) {
                log_Ra_Rdot_tp[j] =log(Ra_Rdot_tp[j]);
                log_Ra_Rdot_tpinf[j]=log(Ra_Rdot_tpinf[j]);
                log_Ra_Rdot_t[j]  =log(Ra_Rdot_t[j]);
            }
                
            for (j=0;j<2;j++) {
                Ra2_R_tpinv[j] = Rainv[j]/Ra_Rdot_tp[j];
                Ra2_R_tpinvinf[j] = 2/Ra_Rdot_tpinf[j];
                Ra2_R_tinv[j] =  Rainv[j]/Ra_Rdot_t[j];
            }
            
            for (j=0;j<2;j++) {
                ylog_Ra_Rdot_tp[j] = yv[j]*log_Ra_Rdot_tp[j];
                yRa2_R_tpinv[j]    = yv[j]*   Ra2_R_tpinv[j];
                yRa2_R_tpinvinf[j] = yv[j]*Ra2_R_tpinvinf[j];
                zlog_Ra_Rdot_t[j]  = zv[j]*log_Ra_Rdot_t[j];
                zRa2_R_tinv[j]     = zv[j]*   Ra2_R_tinv[j];
                
            }

            for (j=0;j<2;j++) {
                y2Ra2_R_tpinv[j] = yv[j]* yRa2_R_tpinv[j];
                y2Ra2_R_tpinvinf[j]=yv[j]* yRa2_R_tpinvinf[j];
                z2Ra2_R_tinv[j]  = zv[j]*  zRa2_R_tinv[j];
                y3Ra2_R_tpinvinf[j]=y2[j]* yRa2_R_tpinvinf[j];
            }

            temp1=denom*(1+c);

            for (j=0;j<2;j++) {
                temp4[j]=temp1*(Ra[j]+(yv[j]+zv[j]));
                temp4[j+2]=temp1*yv[j]*(1-c);
            }
            
            for (j=0;j<2;j++) {
                f_003v[j]=atan(temp4[j]);
                f_003vinf[j]=atan(temp4[j+2]);
            }
            
            temp1=-2.0e0*denom;

            for (j=0;j<2;j++) {
                f_003v[j]*=temp1;
                f_003vinf[j]*=-temp1;
            }

            for (j=0;j<2;j++) {
                adf_003[j]=f_003v[j]*a2_d2;
                adf_003inf[j]=f_003vinf[j]*a2_d2;
            }

            for (j=0;j<2;j++) {
                commonf223[j] = c*Ra[j] - adf_003[j];
                f_103v[j] = c*log_Ra_Rdot_t[j]  - log_Ra_Rdot_tp[j];
                f_103vinf[j] = log_Ra_Rdot_tpinf[j];
                f_013v[j] = c*log_Ra_Rdot_tp[j] - log_Ra_Rdot_t [j];
                f_013vinf[j] = -c*log_Ra_Rdot_tpinf[j];
                f_113v[j] = c*adf_003[j] - Ra[j];
                f_113vinf[j] = c*(adf_003inf[j]-yv[j]);
                f_203vinf[j] = yv[j]-adf_003inf[j];
            }
            
            for (j=0;j<2;j++) {
                commonf223[j] *=  onemc2inv;
                f_103v[j] *=      onemc2inv;
                f_013v[j] *=      onemc2inv;
                f_113v[j] *=      onemc2inv;
                f_103vinf[j] *=   onemc2inv;
                f_013vinf[j] *=   onemc2inv;
                f_113vinf[j] *=   onemc2inv;
                f_203vinf[j] *=   onemc2inv;
            }
        
            for (j=0;j<2;j++) {
                commonf225[j] = f_003v[j] - c*Rainv[j];
                commonf025[j] = c*yRa2_R_tpinv[j] - Rainv[j];
                commonf205[j] = c*zRa2_R_tinv[j]  - Rainv[j];
                commonf305[j] = log_Ra_Rdot_t[j]  - (yv[j]-c*zv[j])*Rainv[j] -
                                c2*z2Ra2_R_tinv[j];
                commonf035[j] = log_Ra_Rdot_tp[j] - (zv[j]-c*yv[j])*Rainv[j] -
                                c2*y2Ra2_R_tpinv[j]; 
                f_203v[j] =  zlog_Ra_Rdot_t[j]  + commonf223[j];
                f_023v[j] =  0.0e0;
                f_023vinf[j]=0.0e0;
                f_005v[j] = f_003v[j] - yRa2_R_tpinv[j] - zRa2_R_tinv[j];
                f_005vinf[j] = f_003vinf[j] + yRa2_R_tpinvinf[j];
                f_105v[j] = Ra2_R_tpinv[j] - c*Ra2_R_tinv[j];
                f_105vinf[j] = -Ra2_R_tpinvinf[j];
                f_015v[j] = Ra2_R_tinv[j]  - c*Ra2_R_tpinv[j];
                f_015vinf[j] = c*Ra2_R_tpinvinf[j];
                f_115v[j] = Rainv[j] - c*(yRa2_R_tpinv[j] + zRa2_R_tinv[j] +
                            f_003v[j]);
            }

            for (j=0;j<2;j++) {
                ycommonf025[j] = yv[j]*commonf025[j];
                zcommonf205[j] = zv[j]*commonf205[j];
                zcommonf305[j] = zv[j]*commonf305[j];
                tf_113[j]=2.0e0*f_113v[j];
                tf_113inf[j]=2.0e0*f_113vinf[j];
                f_205v[j] = yRa2_R_tpinv[j] + c2*zRa2_R_tinv[j] + commonf225[j];
                f_205vinf[j]= f_003vinf[j] - yRa2_R_tpinvinf[j];
                f_025v[j] = zRa2_R_tinv[j] + c2*yRa2_R_tpinv[j] + commonf225[j];
                f_025vinf[j] = f_003vinf[j] - c2*yRa2_R_tpinvinf[j];
                f_305v[j] = y2Ra2_R_tpinv[j] + c*commonf305[j] +
                            2.0e0*f_103v[j];
                f_305vinf[j] = -y2Ra2_R_tpinvinf[j] + 2.0e0*f_103vinf[j];
                f_035v[j] = z2Ra2_R_tinv[j] + c*commonf035[j] + 2.0e0*f_013v[j];
                f_035vinf[j] = c * (c2 * y2Ra2_R_tpinvinf[j] -
                                    log_Ra_Rdot_tpinf[j]) +
                               2.0e0 * f_013vinf[j];
            }
        
            for (j=0;j<2;j++) {
                f_115vinf[j] = -c * f_205vinf[j];
                f_215v[j]    = f_013v[j] - ycommonf025[j] +
                               c * (zcommonf205[j] - f_103v[j]); 
                f_215vinf[j] = f_013vinf[j] +
                               c * (y2Ra2_R_tpinvinf[j] - f_103vinf[j]);
                f_125v[j]    = f_103v[j] - zcommonf205[j] +
                               c * (ycommonf025[j] - f_013v[j]); 
                f_125vinf[j] = f_103vinf[j] - c2 * y2Ra2_R_tpinvinf[j] -
                               c * f_013vinf[j];
                f_225v[j]    = f_203v[j] - zcommonf305[j] +
                               c * (y2[j] * commonf025[j] - tf_113[j]);
                f_225vinf[j] = f_203vinf[j] - c * (c * y3Ra2_R_tpinvinf[j] +
                               tf_113inf[j]);
                f_315v[j]    = tf_113[j] - y2[j] * commonf025[j] +
                               c * (zcommonf305[j] - f_203v[j]);
                f_315vinf[j] = tf_113inf[j] +
                               c * (y3Ra2_R_tpinvinf[j] - f_203vinf[j]);
                f_135v[j]    = 0.0e0;
                f_135vinf[j] = 0.0e0;
            }
            
             
            f_003= (f_003v[0]+f_003vinf[0])-(f_003v[1]+f_003vinf[1]);
            f_013= (f_013v[0]+f_013vinf[0])-(f_013v[1]+f_013vinf[1]);
            f_103= (f_103v[0]+f_103vinf[0])-(f_103v[1]+f_103vinf[1]);
            f_113= (f_113v[0]+f_113vinf[0])-(f_113v[1]+f_113vinf[1]);
            f_023= (f_023v[0]+f_023vinf[0])-(f_023v[1]+f_023vinf[1]);
            f_203= (f_203v[0]+f_203vinf[0])-(f_203v[1]+f_203vinf[1]);
            f_005= (f_005v[0]+f_005vinf[0])-(f_005v[1]+f_005vinf[1]);
            f_015= (f_015v[0]+f_015vinf[0])-(f_015v[1]+f_015vinf[1]);
            f_105= (f_105v[0]+f_105vinf[0])-(f_105v[1]+f_105vinf[1]);
            f_115= (f_115v[0]+f_115vinf[0])-(f_115v[1]+f_115vinf[1]);
            f_025= (f_025v[0]+f_025vinf[0])-(f_025v[1]+f_025vinf[1]);
            f_205= (f_205v[0]+f_205vinf[0])-(f_205v[1]+f_205vinf[1]);
            f_215= (f_215v[0]+f_215vinf[0])-(f_215v[1]+f_215vinf[1]);
            f_125= (f_125v[0]+f_125vinf[0])-(f_125v[1]+f_125vinf[1]);
            f_035= (f_035v[0]+f_035vinf[0])-(f_035v[1]+f_035vinf[1]);
            f_305= (f_305v[0]+f_305vinf[0])-(f_305v[1]+f_305vinf[1]);
            f_225= (f_225v[0]+f_225vinf[0])-(f_225v[1]+f_225vinf[1]);
            f_135= (f_135v[0]+f_135vinf[0])-(f_135v[1]+f_135vinf[1]);
            f_315= (f_315v[0]+f_315vinf[0])-(f_315v[1]+f_315vinf[1]);
            
            
            f_005    *= a2_d2inv;
            f_105    *= onemc2inv;
            f_015    *= onemc2inv;
            f_115    *= onemc2inv;
            f_205    *= onemc2inv;
            f_025    *= onemc2inv;
            f_305    *= onemc2inv;
            f_035    *= onemc2inv; 
            f_215    *= onemc2inv; 
            f_125    *= onemc2inv; 
            f_225    *= onemc2inv;
            f_315    *= onemc2inv;
            f_135    *= onemc2inv;
            
      
/*
 *          Now construct the vector coefficients for the definite integrals
 */
            a2 = a*a;
            m4p = 0.25 * MU / M_PI;
            m4pd =  m4p * d;
            m8p = 0.5 * m4p;
            m8pd = m8p * d;
            m4pn = m4p / ( 1 - NU );
            m4pnd = m4pn * d;
            m4pnd2 = m4pnd * d;
            m4pnd3 = m4pnd2 * d;
            a2m4pnd = a2 * m4pnd;
            a2m8pd = a2 * m8pd;
            a2m4pn = a2 * m4pn;
            a2m8p = a2 * m8p;

            
            for (i=0;i<3;i++) {
                tpct[i]=-tctp[i];
                tcbp[i]=t[alt1[i]]*bp[alt2[i]]-t[alt2[i]]*bp[alt1[i]];
                tpcb[i]=tp[alt1[i]]*b[alt2[i]]-tp[alt2[i]]*b[alt1[i]];
                bct[i]=b[alt1[i]]*t[alt2[i]]-b[alt2[i]]*t[alt1[i]];
                bpctp[i]=bp[alt1[i]]*tp[alt2[i]]-bp[alt2[i]]*tp[alt1[i]];
                
            }

            tdb=0.0e0;
            tdbp=0.0e0;
            tpdb=0.0e0;
            tpdbp=0.0e0;
            tctpdb=0.0e0;
            tpctdbp=0.0e0;
            bpctpdb=0.0e0;
            bctdbp=0.0e0;
            
            for (i=0;i<3;i++) {
                tdb    +=t[i]*b[i];
                tdbp   +=t[i]*bp[i];
                tpdb   +=tp[i]*b[i];
                tpdbp  +=tp[i]*bp[i];
                tctpdb +=tctp[i]*b[i];
                tpctdbp+=tpct[i]*bp[i];
                bpctpdb+=bpctp[i]*b[i];
                bctdbp +=bct[i]*bp[i];
            }
            
            for (i=0;i<3;i++) {
                tctpct[i]    =        tp[i] -     c*t[i];
                tpctctp[i]   =         t[i] -    c*tp[i];
                tctpcbp[i]   =   tdbp*tp[i] - tpdbp*t[i];
                tpctcb[i]    =    tpdb*t[i] -  tdb*tp[i];
                tcbpct[i]    =        bp[i] -  tdbp*t[i];
                tpcbctp[i]   =         b[i] - tpdb*tp[i];
                bpctpct[i]   =   tdbp*tp[i] -    c*bp[i];
                bctctp[i]    =    tpdb*t[i] -     c*b[i];
                tctpcbpct[i] = tdbp*tpct[i];
                tpctcbctp[i] = tpdb*tctp[i];
            }
                
            
            tctpcbpdtp = tdbp - tpdbp*c;
            tpctcbdt = tpdb - tdb*c;
            tctpcbpdb =  tdbp*tpdb - tpdbp*tdb;
            tpctcbdbp = tctpcbpdb;
            tcbpdtp = tpctdbp; 
            tpcbdt = tctpdb;
            tcbpdb = bctdbp;
            tpcbdbp = bpctpdb;

/*
 *          Calculate the forces for segment p3->p4
 */
            temp1 = tdbp*tpdb + tctpcbpdb;

            for (i=0;i<3;i++) {
                I00a[i] = temp1 * tpct[i];
                I00b[i] = tctpcbpdtp * bct[i];
            }

            temp1 = (m4pnd * tctpdb);
            temp2 = (m4pnd * bpctpdb);
            temp3 = (m4pnd3 * tctpcbpdtp*tctpdb);
                
            for (i=0;i<3;i++) {
                I_003[i] = m4pd*I00a[i] - m4pnd*I00b[i] + temp1*bpctpct[i] +
                           temp2*tctpct[i]; 
                I_005[i] = a2m8pd*I00a[i] - a2m4pnd*I00b[i] -
                           temp3*tctpct[i];
                I10a[i] = tcbpct[i]*tpdb - tctp[i]*tcbpdb;
                I10b[i] = bct[i] * tcbpdtp;
            }

            temp1 = (m4pn * tdb);
            temp2 = m4pnd2 * (tcbpdtp*tctpdb + tctpcbpdtp*tdb);
                
            for (i=0;i<3;i++) {
                I_103[i] = temp1*bpctpct[i] + m4p*I10a[i] - m4pn*I10b[i];
                I_105[i] = a2m8p*I10a[i] - a2m4pn*I10b[i] - temp2*tctpct[i];
                I01a[i] = tctp[i]*bpctpdb - bpctpct[i]*tpdb;
            }

            tmp[0] = (m4pn * tpdb); 
            tmp[1] = (m4pn * bpctpdb);
            tmp[2] = (m4pnd2 * tctpcbpdtp * tpdb);
            tmp[3] = (m4pnd2 * tctpcbpdtp * tctpdb);
            tmp[4] = (m4pnd * tcbpdtp * tdb);
            tmp[5] = (m4pnd * tctpcbpdtp * tpdb) ;
            tmp[6] = (m4pnd * (tctpcbpdtp*tdb + tcbpdtp*tctpdb));
            tmp[7] = (m4pnd * tcbpdtp * tpdb);
            tmp[8] = (m4pn * tcbpdtp * tdb);
            tmp[9] = (m4pn * tcbpdtp * tpdb);
                
            for (i=0;i<3;i++) {
                I_013[i] = m4p*I01a[i] + tmp[0]*bpctpct[i] - tmp[1]*tctp[i];
                I_015[i] = a2m8p*I01a[i] - tmp[2]*tctpct[i] +
                           tmp[3]*tctp[i];
                I_205[i] = -tmp[4] * tctpct[i];
                I_025[i] = tmp[5] * tctp[i]; 
                I_115[i] = tmp[6]*tctp[i] - tmp[7]*tctpct[i];
                I_215[i] = tmp[8] * tctp[i];
                I_125[i] = tmp[9] * tctp[i];
            }
  
            Fint_003 = f_103 - y[0]*f_003;
            Fint_103 = f_203 - y[0]*f_103;
            Fint_013 = f_113 - y[0]*f_013;
            Fint_005 = f_105 - y[0]*f_005;
            Fint_105 = f_205 - y[0]*f_105;
            Fint_015 = f_115 - y[0]*f_015;
            Fint_115 = f_215 - y[0]*f_115;
            Fint_205 = f_305 - y[0]*f_205;
            Fint_025 = f_125 - y[0]*f_025;
            Fint_215 = f_315 - y[0]*f_215;
            Fint_125 = f_225 - y[0]*f_125;
                
            for (i=0;i<3;i++) {
                f4[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 +
                       I_013[i]*Fint_013 + I_005[i]*Fint_005 +
                       I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                       I_115[i]*Fint_115 + I_205[i]*Fint_205 +
                       I_025[i]*Fint_025 + I_215[i]*Fint_215 +
                       I_125[i]*Fint_125) * oneoverL;
            }

            Fint_003 = y[1]*f_003 - f_103;
            Fint_103 = y[1]*f_103 - f_203;
            Fint_013 = y[1]*f_013 - f_113;
            Fint_005 = y[1]*f_005 - f_105;
            Fint_105 = y[1]*f_105 - f_205;
            Fint_015 = y[1]*f_015 - f_115;
            Fint_115 = y[1]*f_115 - f_215;
            Fint_205 = y[1]*f_205 - f_305;
            Fint_025 = y[1]*f_025 - f_125;
            Fint_215 = y[1]*f_215 - f_315;
            Fint_125 = y[1]*f_125 - f_225;
                
            for (i=0;i<3;i++) {
                f3[i]=(I_003[i]*Fint_003 + I_103[i]*Fint_103 +
                       I_013[i]*Fint_013 + I_005[i]*Fint_005 +
                       I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                       I_115[i]*Fint_115 + I_205[i]*Fint_205 +
                       I_025[i]*Fint_025 + I_215[i]*Fint_215 +
                       I_125[i]*Fint_125) * oneoverL;
            }

            *fp3x=f3[0];
            *fp3y=f3[1];
            *fp3z=f3[2];
            *fp4x=f4[0];
            *fp4y=f4[1];
            *fp4z=f4[2];

#if 0
/*
 *          This portion will calculate the force from segment p3/p4 on
 *          the endpoint of the semi-infinite segment, but for now
 *          we do not need that value so we don't bother with it.
 */
            temp1 = tpdb*tdbp + tpctcbdbp;

            for (i=0;i<3;i++) {
                I00a[i] = temp1 * tctp[i];
                I00b[i] = bpctp[i] * tpctcbdt;
            }
                
            temp1 = m4pnd * tpctdbp;
            temp2 = m4pnd * bctdbp;
            temp3 = m4pnd3 * tpctcbdt * tpctdbp;
                
            for (i=0;i<3;i++) {
                I_003[i] = m4pd*I00a[i] - m4pnd*I00b[i] + temp1*bctctp[i] +
                           temp2*tpctctp[i];
                I_005[i] = a2m8pd*I00a[i] - a2m4pnd*I00b[i] -
                           temp3*tpctctp[i]; 
                I01a[i] = tpct[i]*tpcbdbp - tpcbctp[i]*tdbp;
                I01b[i] = -bpctp[i] * tpcbdt;
            }

            temp1 = m4pn * tpdbp;
            temp2 = m4pnd2 * (tpcbdt*tpctdbp + tpctcbdt*tpdbp);
                
            for (i=0;i<3;i++) {
                I_013[i] = -temp1 * bctctp[i] + m4p*I01a[i] - m4pn*I01b[i];
                I_015[i] = a2m8p*I01a[i] - a2m4pn*I01b[i] +
                           temp2*tpctctp[i];
                I10a[i] = bctctp[i]*tdbp - tpct[i]*bctdbp;
            }

            tmp[0] = m4pn * tdbp; 
            tmp[1] = m4pn * bctdbp;
            tmp[2] = m4pnd2 * tpctcbdt * tdbp;
            tmp[3] = m4pnd2 * tpctcbdt * tpctdbp;
            tmp[4] = (m4pnd * tpcbdt * tpdbp);
            tmp[5] = (m4pnd * tpctcbdt * tdbp);
            tmp[6] = m4pnd * (tpctcbdt*tpdbp + tpcbdt*tpctdbp);
            tmp[7] = m4pnd * tpcbdt * tdbp;
            tmp[8] = (m4pn * tpcbdt * tpdbp);
            tmp[9] = (m4pn * tpcbdt * tdbp);
                
            for (i=0;i<3;i++) {
                I_103[i] = m4p*I10a[i] - tmp[0]*bctctp[i] + tmp[1]*tpct[i];
                I_105[i] = a2m8p*I10a[i] + tmp[2]*tpctctp[i] -
                           tmp[3]*tpct[i];
                I_025[i] = -tmp[4] * tpctctp[i];
                I_205[i] = tmp[5] * tpct[i];
                I_115[i] = tmp[6]*tpct[i] - tmp[7]*tpctctp[i];
                I_125[i] = -tmp[8] * tpct[i];
                I_215[i] = -tmp[9] * tpct[i];
            }

            Fint_003 = f_003;
            Fint_103 = f_103;
            Fint_013 = f_013;
            Fint_005 = f_005;
            Fint_105 = f_105;
            Fint_015 = f_015;
            Fint_115 = f_115;
            Fint_205 = f_205;
            Fint_025 = f_025;
            Fint_215 = f_215;
            Fint_125 = f_125;
                
            for (i=0;i<3;i++) {
                f1[i]=I_003[i]*Fint_003 + I_103[i]*Fint_103 +
                      I_013[i]*Fint_013 + I_005[i]*Fint_005 +
                      I_105[i]*Fint_105 + I_015[i]*Fint_015 +
                      I_115[i]*Fint_115 + I_205[i]*Fint_205 +
                      I_025[i]*Fint_025 + I_215[i]*Fint_215 +
                      I_125[i]*Fint_125;
            }
                
            *fp1x=f1[0];
            *fp1y=f1[1];
            *fp1z=f1[2];
   
#endif  /* Code for calculating force on endpoint of semi-infinite seg */

        } else {
/*
 *          The two lines are parallel, so we have to use a special
 *          lower dimensional function
 */
            SpecialSemiInfiniteSegSegForceNonRational(home,
                    p1x, p1y, p1z, p2x, p2y, p2z,
                    p3x, p3y, p3z, p4x, p4y, p4z,
                    bpx, bpy, bpz, bx, by, bz, a, MU, NU,
                    eps, fp1x, fp1y, fp1z,
                    fp3x, fp3y, fp3z, fp4x, fp4y, fp4z);
       }

       return;
}
#endif  /* ifdef _ARLFEM */
