#ifdef CALCENERGY

/**************************************************************************
 *
 *      Module:  This module contains the functions needed for
 *               calculating interactions between dislocation
 *               segments. 
 *
 *      Author:  Sylvie Aubry based on some of Tom Arsenlis Matlab and C 
 *               routines 
 *
 *      Date :   Sept 2016
 *
 *      Includes public functions:
 *              SelfEnergy()
 *              ExtPKEnergy()
 *              SegSegEnergy()
 *              ComputeEnergy()
 *              LineTensionEnergy()
 *
 *      Includes private functions:
 *              SpecialSegSegEnergy()
 *
 *************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"
#include "ParadisThread.h"


/*---------------------------------------------------------------------------
 *
 *    Function:        SelfEnergy
 *    Description:     Computes the self energy of a given segment
 *
 *--------------------------------------------------------------------------*/
real8 SelfEnergy(int coreOnly, real8 MU, real8 NU,
		real8 bx, real8 by, real8 bz,
		real8 x1, real8 y1, real8 z1,
		real8 x2, real8 y2, real8 z2,
		real8 a,  real8 Ecore)
{
  real8 tx, ty, tz, L, La, S, S0;
  real8 bs, bs2, bex, bey, bez, be2;  

  real8 W = 0.0; 
 
  GetUnitVector(1, x1, y1, z1, x2, y2, z2, &tx, &ty, &tz, &L);
  bs = bx*tx + by*ty + bz*tz;
  bex = bx-bs*tx; bey = by-bs*ty; bez=bz-bs*tz;
  be2 = (bex*bex+bey*bey+bez*bez);
  bs2 = bs*bs;
  
  La=sqrt(L*L+a*a);
  
  if(coreOnly) 
    {
      S = 0.0;
    }
  else
    {
      S0 = 2*L*log((La+L)/a);
      S = (0.125*MU/M_PI)*(be2* S0/(1-NU) + bs2*(S0 -((3-NU)/(1-NU))*(La - a) ));
    }	

  W  = S +  Ecore*(bs2+be2/(1-NU))*L;
  return W;
}


real8 ExtPKEnergy(real8 s[3][3], real8 E, real8 NU,
		  real8 bx, real8 by, real8 bz,
		  real8 x1, real8 y1, real8 z1,
		  real8 x2, real8 y2, real8 z2)
{
  // Use Hirth and Lothe's formula p.45

  return (s[0][0] + s[1][1] + s[2][2])*(s[0][0] + s[1][1] + s[2][2])/(2*E) + 
    (s[1][2]*s[1][2] + s[2][0]*s[2][0] + s[0][1]*s[0][1] - 
     s[0][0]*s[1][1] - s[0][0]*s[2][2] - s[1][1]*s[2][2])*(1+NU)/E;
  
}

static real8 SpecialSegSegEnergy(real8 p1x, real8 p1y, real8 p1z,
				real8 p2x, real8 p2y, real8 p2z,
				real8 p3x, real8 p3y, real8 p3z,
				real8 p4x, real8 p4y, real8 p4z,
				real8 bpx, real8 bpy, real8 bpz,
				real8 bx, real8 by, real8 bz,
				real8 a, real8 MU, real8 NU, real8 ecrit)                               
{                                
  int i, j;
  real8 x1[3], x2[3], x3[3], x4[3], b[3], bp[3];
  real8 vec1[3], vec2[3], t[3], tp[3], nd[3];
  real8 R[3], Rdt;
  real8 d2, c, oneoverL, oneoverLp;
  real8 a2, a2_d2, a2d2inv;
  real8 Ra[4];
  real8 temp, tmp[8];
  real8 temp1, temp2, temp3, temp4;
  real8 W_001v[4], W_003v[4], W_113v[4], W_223v[4];
  real8 W_001, W_003, W_113, W_223;
  real8 I_001, I_003, I_113, I_223;
  real8 tdb, tdbp, bdbp;
  real8 omn, factor;
  real8 cotanthetac, eps;
  real8 magdiff, diffMag2, x1modMag2, x2modMag2;
  real8 x1mod[3], x2mod[3];
  real8 W12corr,y[2],z[2],yv[4],zv[4],ypz[4],Log_Ra_ypz[4];
  real8 opn,a2opnd2,nddb,nddbp;
  real8 W12=0.0;
  
  a2= a*a;
  
  cotanthetac = sqrt((1 - ecrit*1.01) / (ecrit*1.01));
  eps=1e-16;


  
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
  
  for(i=0;i<3;i++)
    { 
      vec1[i]=x4[i]-x3[i];
      vec2[i]=x2[i]-x1[i];
    }
  temp1=0.0e0;
  temp3=0.0e0;
  for(i=0;i<3;i++)
    { 
      temp1+=vec1[i]*vec1[i];
      temp3+=vec2[i]*vec2[i];
    }  
  
  temp2=1/temp1;
  temp4=1/temp3;
  
  oneoverL =sqrt(temp2);
  oneoverLp=sqrt(temp4);
  for(i=0;i<3;i++)
    { 
      t[i]=vec1[i]*oneoverL;
      tp[i]=vec2[i]*oneoverLp;
    }
  
  c=0.0e0;
  for(i=0;i<3;i++)
    { 
      c+=t[i]*tp[i];
    }

  if (c < 0) 
    {
      for(i=0;i<3;i++)
	{ 
	  temp=x2[i];
	  x2[i]=x1[i];
	  x1[i]=temp;
	  vec2[i]=-vec2[i];
	  bp[i]=-bp[i];
	  tp[i]=-tp[i];
	}
       c=-c;
    }

  temp=0.0e0;
  
  for (i=0;i<3;i++) {
    temp+=vec2[i]*t[i];
  }
  
  for (i=0;i<3;i++) {
    x2mod[i]=x1[i]+temp*t[i];
  }
  
  for (i=0;i<3;i++) {
    vec2[i]=x2[i]-x2mod[i];
  }
  
  temp=0.0e0;
  
  for (i=0;i<3;i++) {
    temp+=vec2[i]*vec2[i];
  }
  
  magdiff=sqrt(temp);
  temp=magdiff*0.5e0 * cotanthetac;
  
  for (i=0;i<3;i++) {
    vec1[i]=temp*t[i];
  }
  
  for (i=0;i<3;i++) {
    x1mod[i]=x1[i]+0.5e0*vec2[i]+vec1[i];
    x2mod[i]+=0.5e0*vec2[i]-vec1[i];
  }
  
  for (i=0;i<3;i++) {
    R[i]=x3[i]-x1mod[i];
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
    z[0]+=-x1mod[i]*t[i];
    z[1]+=-x2mod[i]*t[i];
  } 
  
  for (j=0;j<2;j++) {
    yv[2*j]=y[j];
    yv[2*j+1]=y[j];
    zv[j]=z[j];
    zv[j+2]=z[j];
  }    
  
  a2_d2 = a2 + d2;        
  a2d2inv = 1.0e0 / a2_d2;
  
  for (j=0;j<4;j++) {
    ypz[j] = yv[j] + zv[j];
  }
  
  for (j=0;j<4;j++) {
    tmp[j]=a2_d2 + ypz[j]*ypz[j];
  }
  
  for (j=0;j<4;j++) {
    Ra[j]=sqrt(tmp[j]);
  }
  
  
  for (j=0;j<4;j++) {
    tmp[j]=Ra[j] + ypz[j];
    tmp[j+4]=Ra[j]-ypz[j];
  }

  for (j=0;j<4;j++) {
    Log_Ra_ypz[j]=0.5e0*(log(tmp[j])-log(tmp[j+4]));
  }
  
  for (j=0;j<4;j++) {
    W_001v[j]=ypz[j]*Log_Ra_ypz[j]-Ra[j];
    W_003v[j]=Ra[j]*a2d2inv;
    W_113v[j]=-Log_Ra_ypz[j];
    W_223v[j]=ypz[j]*Log_Ra_ypz[j]-2*Ra[j];
  }
  
  W_001= W_001v[0]+W_001v[3]-W_001v[1]-W_001v[2]; 
  W_003= W_003v[0]+W_003v[3]-W_003v[1]-W_003v[2];
  W_113= W_113v[0]+W_113v[3]-W_113v[1]-W_113v[2];
  W_223= W_223v[0]+W_223v[3]-W_223v[1]-W_223v[2];
  
  opn=1+NU;
  omn=1-NU;
  a2opnd2=0.5*a2*opn;
  
  tdb=0.0e0;
  tdbp=0.0e0;
  bdbp=0.0e0;
  nddb=0.0e0;
  nddbp=0.0e0;

  for (i=0;i<3;i++) {
    tdb += t[i]*b[i];
    tdbp+= t[i]*bp[i];
    bdbp+= b[i]*bp[i];
    nddb+= nd[i]*b[i];
    nddbp+= nd[i]*bp[i];
  }

  I_001=opn*tdb*tdbp-bdbp;
  I_003=a2opnd2*tdb*tdbp - a2*bdbp - nddb*nddbp;
  I_113=-tdb*nddbp-tdbp*nddb;
  I_223=-tdb*tdbp;


  factor=(0.25*MU/(M_PI*(omn)));

  W12=I_001*W_001+I_003*W_003+I_113*W_113+I_223*W_223;
    
  W12 *=factor;

 x1modMag2 = 0.0e0;
  x2modMag2 = 0.0e0;
  
  for (i=0;i<3;i++) {
    x1modMag2 += x1mod[i]*x1mod[i];
    x2modMag2 += x2mod[i]*x2mod[i];
  }

  diffMag2 = magdiff*magdiff;

  if (diffMag2 > (eps * (x1modMag2+x2modMag2))) {
   W12corr =  SegSegEnergy(x1[0], x1[1], x1[2], x1mod[0], x1mod[1], x1mod[2],
			   x3[0], x3[1], x3[2], x4[0], x4[1], x4[2],
			   bp[0], bp[1], bp[2], b[0], b[1], b[2], a, MU, NU);

    W12 += W12corr;

    W12corr = SegSegEnergy(x2mod[0], x2mod[1], x2mod[2], x2[0], x2[1], x2[2],
			   x3[0], x3[1], x3[2], x4[0], x4[1], x4[2],
			   bp[0], bp[1], bp[2], b[0], b[1], b[2], a, MU, NU);
    W12 += W12corr;
  }

  return W12;
}



/*-------------------------------------------------------------------------
 *
 *      Function:       SegSegEnergy
 *      Description:    Used to calculate the interaction energy between
 *                      two dislocation segments analytically.
 *
 *      Arguments:
 *              p1*,p2*      endpoints for first dislocation segment starting
 *                           at p1x,p1y,p1z and ending at p2x,p2y,p2z
 *              p3*,p4*      endpoints for seond dislocation segment starting
 *                           at p3x,p3y,p3z and ending at p4x,p4y,p4z
 *              bxp,byp,bzp  burgers vector for segment p1 to p2
 *              bx,by,bz     burgers vector for segment p3 to p4
 *              a            core parameter
 *              MU           shear modulus
 *              NU           poisson ratio
 *              W            segment-segment interaction energy
 *                      
 *-----------------------------------------------------------------------*/
real8 SegSegEnergy(real8 p1x, real8 p1y, real8 p1z,
		  real8 p2x, real8 p2y, real8 p2z,
		  real8 p3x, real8 p3y, real8 p3z,
		  real8 p4x, real8 p4y, real8 p4z,
		  real8 bpx, real8 bpy, real8 bpz,
		  real8 bx, real8 by, real8 bz,
		  real8 a, real8 MU, real8 NU)
{
  real8 W12=0.0;
  real8 x1[3], x2[3], x3[3], x4[3], b[3], bp[3];
  real8 vec1[3], vec2[3], t[3], tp[3], tctp[3];
  real8 R[4][3], R2[4], Rdt[4], Rdtp[4];
  real8 eps, d, d2, c, c2, onemc2, onemc2inv, oneoverL, oneoverLp;
  real8 yv[4], zv[4], Ra[4], Rainv[4];
  real8 a2, a2_d2, denom;
  real8 Ra_Rdot_tp[4], Ra_Rdot_t[4], log_Ra_Rdot_tp[4], log_Ra_Rdot_t[4];
  real8 ylog_Ra_Rdot_tp[4], zlog_Ra_Rdot_t[4];
  real8 adW_003[4], commonW223[4];
  real8 W_001v[4], W_003v[4], W_103v[4], W_013v[4], W_113v[4], W_203v[4], W_023v[4];
  real8 W_001, W_003, W_103, W_013, W_113, W_203, W_023;
  real8 I_001, I_003, I_103, I_013, I_113, I_203, I_023;
  real8 tmp[8], temp1, temp2, temp4[4];
  real8 omn, a2omnd2, twonu, a2nu, factor;
  real8 tdb, tdbp, tpdb, tpdbp, bdbp, tctpdb, tctpdbp;
  real8 common1, common2;
  int i, j, alt1[3]={1,2,0}, alt2[3]={2,0,1};
  
  eps = 1e-4;            

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
  
  for(i=0;i<3;i++)
    { 
      vec1[i]=x4[i]-x3[i];
      vec2[i]=x2[i]-x1[i];
    }
  temp1=0.0e0;
  temp2=0.0e0;
  
  for(i=0;i<3;i++)
    { 
      temp1+=vec1[i]*vec1[i];
      temp2+=vec2[i]*vec2[i];
    }
  oneoverL =1/sqrt(temp1);
  oneoverLp=1/sqrt(temp2);
  
  for(i=0;i<3;i++)
    { 
      t[i]=vec1[i]*oneoverL;
      tp[i]=vec2[i]*oneoverLp;
    }
  
  c=0.0e0;
  for(i=0;i<3;i++)
    { 
      c+=t[i]*tp[i];
    }
  c2=c*c;
  onemc2=1-c2;
  
  if (onemc2 > eps)
    {
      onemc2inv = 1/onemc2;
      
      for(i=0;i<3;i++)
	{
	  tctp[i]=t[alt1[i]]*tp[alt2[i]]-t[alt2[i]]*tp[alt1[i]];
	}
       for (i=0;i<3;i++)
	{ 
	  R[0][i]=x3[i]-x1[i];
	  R[1][i]=x3[i]-x2[i];
	  R[2][i]=x4[i]-x1[i];  
	  R[3][i]=x4[i]-x2[i];
	}
      for (i=0;i<4;i++)
	{ 
	  R2[i]=0.0e0;
	  Rdt[i]=0.0e0;
	  Rdtp[i]=0.0e0;
	  for (j=0;j<3;j++)
	    {
	      R2[i]+=R[i][j]*R[i][j];
	      Rdt[i]+=R[i][j]*t[j];
	      Rdtp[i]+=R[i][j]*tp[j];
	    }
	}
      d=0.0e0;
      for(i=0;i<3;i++)
	{ 
	  d+=(R[0][i]+R[3][i])*tctp[i];
	}

      d*=0.5e0*onemc2inv;
      for (j=0;j<4;j++)
	{
	  yv[j]=Rdt[j]-c*Rdtp[j];
	  zv[j]=Rdtp[j]-c*Rdt[j];
	}
      
      d2=d*d;
      a2=a*a;
      a2_d2 = a2+d2*onemc2;

      for (j=0;j<4;j++)
	{	  
	  temp4[j]=R2[j]+a2;
	}
      temp1=onemc2*a2_d2;
      for (j=0;j<4;j++)
	{
	  Ra[j]=sqrt(temp4[j]);
	}
      temp2=sqrt(temp1);
      for (j=0;j<4;j++)
	{
	  Rainv[j]=1.0e0/Ra[j];
	}
      denom=1.0e0/temp2;
      
      for (j=0;j<4;j++)
	{
	  tmp[j]=0.5e0*Rainv[j];
	}
      
      for (j=0;j<4;j++)
	{
	  Ra_Rdot_tp[j] = Ra[j]+Rdtp[j];       
	  Ra_Rdot_t[j]  = Ra[j]+Rdt[j];
	}
      
      for (j=0;j<4;j++)
	{
	  log_Ra_Rdot_tp[j] =log(Ra_Rdot_tp[j]);
	  log_Ra_Rdot_t[j]  =log(Ra_Rdot_t[j]);
	}
      
      for (j=0;j<4;j++)
	{
	  ylog_Ra_Rdot_tp[j] = yv[j]*log_Ra_Rdot_tp[j];
	  zlog_Ra_Rdot_t[j]  = zv[j]*log_Ra_Rdot_t[j];
	}
     
      temp2=(1+c);
      for (j=0;j<4;j++)
	{
	  tmp[j]=temp2*Ra[j]+Rdt[j]+Rdtp[j];
	}
      for (j=0;j<4;j++)
	{
	  temp4[j]=denom*tmp[j];
	}
      
      for (j=0;j<4;j++)
	{
	  W_003v[j]=atan(temp4[j]);
	}
      
      for (j=0;j<4;j++)
	{
	  tmp[j]=1.0e0/tmp[j];
	}
      
      temp1=-2.0e0*denom;
      for (j=0;j<4;j++)
	{
	  W_003v[j]*=temp1;
	}
      for (j=0;j<4;j++)
	{
	  adW_003[j]=W_003v[j]*a2_d2;
	}
      
      for (j=0;j<4;j++)
	{
	   commonW223[j] = c*Ra[j] - adW_003[j];
	   W_113v[j] = c*adW_003[j] - Ra[j];
	}
      for (j=0;j<4;j++)
	{
	  W_003v[j] *= onemc2;
	  adW_003[j] *=onemc2;
	}
      
      for (j=0;j<4;j++)
	{
	  W_001v[j]= zlog_Ra_Rdot_t[j] + ylog_Ra_Rdot_tp[j] - adW_003[j];
	  W_103v[j] = c*log_Ra_Rdot_t[j]  - log_Ra_Rdot_tp[j];
	  W_013v[j] = c*log_Ra_Rdot_tp[j] - log_Ra_Rdot_t[j];
	}
      
      for (j=0;j<4;j++)
	{
	  W_203v[j] =  zlog_Ra_Rdot_t[j]  + commonW223[j];
	  W_023v[j] =  ylog_Ra_Rdot_tp[j] + commonW223[j];
	}
      
      W_001= W_001v[0]+W_001v[3]-W_001v[1]-W_001v[2]; 
      W_003= W_003v[0]+W_003v[3]-W_003v[1]-W_003v[2];
      W_013= W_013v[0]+W_013v[3]-W_013v[1]-W_013v[2];
      W_103= W_103v[0]+W_103v[3]-W_103v[1]-W_103v[2];
      W_113= W_113v[0]+W_113v[3]-W_113v[1]-W_113v[2];
      W_023= W_023v[0]+W_023v[3]-W_023v[1]-W_023v[2];
      W_203= W_203v[0]+W_203v[3]-W_203v[1]-W_203v[2];
      
      /* now construct the vector coefficients for the definite integrals */
      omn=1-NU;
      a2omnd2=0.5*a2*omn;
      twonu=2*NU;
      a2nu=a2*NU;
      
      
      tdb=0.0e0;
      tdbp=0.0e0;
      tpdb=0.0e0;
      tpdbp=0.0e0;
      bdbp=0.0e0;
      tctpdb=0.0e0;
      tctpdbp=0.0e0;
      
      for (i=0;i<3;i++)
	{
	  tdb    +=t[i]*b[i];
	  tdbp   +=t[i]*bp[i];
	  tpdb   +=tp[i]*b[i];
	  tpdbp  +=tp[i]*bp[i];
	  bdbp   +=bp[i]*b[i];
	  tctpdb +=tctp[i]*b[i];
	  tctpdbp+=tctp[i]*bp[i];
	}
      
      common1=tdb*tpdbp;
      common2=tdbp*tpdb;
      
      tmp[0]=a2*bdbp+d2*tctpdb*tctpdbp;
      tmp[1]=tdb*tctpdbp+tdbp*tctpdb;
      tmp[2]=d*c;
      tmp[3]=tpdb*tctpdbp+tpdbp*tctpdb;
      tmp[4]=common2+common1;
      tmp[5]=tdb*tdbp;
      tmp[6]=tpdb*tpdbp;
      
      I_001=omn*common1+twonu*common2-bdbp*c;
      I_003=a2omnd2*common1+a2nu*common2-tmp[0]*c;
      I_103=-tmp[2]*tmp[1];
      I_013=-tmp[2]*tmp[3];
      I_113=-c*tmp[4];
      I_203=-c*tmp[5];
      I_023=-c*tmp[6];
      
      temp2=(0.25*MU/(M_PI*(omn)));
      factor=temp2*onemc2inv;

      W12=I_001*W_001+I_003*W_003+I_013*W_013+I_103*W_103+I_113*W_113+I_023*W_023+I_203*W_203;
      
      W12 *=factor; 

    } else {
    /*
     *          The two lines are parallel, so we have to use a special
     *          lower dimensional function
     */
    W12 = SpecialSegSegEnergy(p1x, p1y, p1z, p2x, p2y, p2z,
                              p3x, p3y, p3z, p4x, p4y, p4z,
                              bpx, bpy, bpz, bx, by, bz, a, MU, NU, eps);
  }
 

#if 0 // For debugginf
  printf("x1 = [%f %f %f];\n",p1x,p1y,p1z);
  printf("x2 = [%f %f %f];\n",p2x,p2y,p2z);
  printf("x3 = [%f %f %f];\n",p3x,p3y,p3z);
  printf("x4 = [%f %f %f];\n",p4x,p4y,p4z);
  printf("bp = [%f %f %f];\n",bpx,bpy,bpz);
  printf("b  = [%f %f %f];\n",bx,by,bz);
  
  printf("W=%f;\n",W12);
#endif
  return W12;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     ComputeEnergy
 *      Description:  Obtains nodal information (i.e coordinates,
 *                    burgers vectors) of the endpoints defining
 *                    the segments, adjusts the coordinates (if
 *                    necessary) for boundary conditions, invokes
 *                    the function to calculate the interaction
 *                    between the two segments, and returns the
 *                    resulting energy at all four endpoints.
 *
 *      Arguments:
 *                    node1    first endpoint of segment node1<-->node2
 *                    node2    second endpoint of segment node1<-->node2
 *                    node3    first endpoint of segment node3<-->node4
 *                    node4    second endpoint of segment node3<-->node4
 *                    W        energy of the interactions from 
 *                             SegSegEnergy only.
 *
 *-----------------------------------------------------------------------*/
real8 ComputeEnergy(Home_t *home, Node_t *node1, Node_t *node2,
                   Node_t *node3, Node_t *node4)
{
        int     armID12, armID34;
        real8   xCenter, yCenter, zCenter;
        real8   a, MU, NU;
        real8   x1, x2, x3, x4;
        real8   y1, y2, y3, y4;
        real8   z1, z2, z3, z4;
        real8   dx, dy, dz;
        real8   b12x, b12y, b12z;
        real8   b34x, b34y, b34z;
	real8   W = 0.0;
        Param_t *param;
        Cell_t  *cell;

        param      = home->param;
/*
 *      This function is only used during the full energy calculations
 *      where energy for any given segment pair are calculated by only
 *      one domain in the problem.  This means that even if one of the
 *      segments is non-local, we still need its energy... so in this
 *      routine, we treat all segments as local so SegSegEnergy() sets
 *      energy for all nodes.
 */
        a  = param->rc;
        MU = param->shearModulus;
        NU = param->pois;

        x1 = node1->x;
        y1 = node1->y;
        z1 = node1->z;

        dx = node2->x - x1;
        dy = node2->y - y1;
        dz = node2->z - z1;

        ZImage(param, &dx, &dy, &dz);
        
/*
 *      It is possible to have a zero-length segment (created by collision
 *      handling).  If we find such a segment, there will be no seg/seg
 *      energy between the given segment pair, so just return zero energy.
 */
        if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) return 0.0;

        x2 = x1 + dx;
        y2 = y1 + dy;
        z2 = z1 + dz;

/*
 *      Convert the coordinates of node 3 to those of the image
 *      of that point nearest the center of the cell containing node1
 */
        cell = LookupCell(home, node1->cellIdx);


        xCenter = cell->center[X];
        yCenter = cell->center[Y];
        zCenter = cell->center[Z];

        x3 = node3->x;
        y3 = node3->y;
        z3 = node3->z;

        x4 = node4->x;
        y4 = node4->y;
        z4 = node4->z;

        PBCPOSITION(param, xCenter, yCenter, zCenter, &x3, &y3, &z3);
        PBCPOSITION(param, x3, y3, z3, &x4, &y4, &z4);

/*
 *      It is possible to have a zero-length segment (created by collision
 *      handling).  If we find such a segment, there will be no seg/seg
 *      energy between the given segment pair, so just return zero energy.
 */
        dx = x3 - x4;
        dy = y3 - y4;
        dz = z3 - z4;

        if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) return 0.0;

        armID12 = GetArmID(node1, node2);
        armID34 = GetArmID(node3, node4);

        b12x = node1->burgX[armID12];
        b12y = node1->burgY[armID12];
        b12z = node1->burgZ[armID12];

        b34x = node3->burgX[armID34];
        b34y = node3->burgY[armID34];
        b34z = node3->burgZ[armID34];

        W = SegSegEnergy(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4,
			 b12x, b12y, b12z, b34x, b34y, b34z,        
			 a, MU, NU);

        return W;
}

/*-------------------------------------------------------------------------
 *
 *      Function:       LineTensionEnergy
 *      Description:    Wrapper function which invokes an appropriate energy
 *                      function based on whether the code has been compiled
 *                      to do isotropic elasticity
 *
 *      Arguments:
 *
 *------------------------------------------------------------------------*/
real8 LineTensionEnergy(Home_t *home, real8 x1, real8 y1, real8 z1,
                       real8 x2, real8 y2, real8 z2,
                       real8 bx, real8 by, real8 bz)
{
        real8   E, MU, NU, TensionFactor, Ecore, a;
        real8   dx, dy, dz;
        real8   WSelf=0.0, WPK=0.0, W = 0.0;
        real8   extStress[3][3];
        Param_t *param;

        param = home->param;

        E = param->YoungsModulus;
        MU = param->shearModulus;
        NU = param->pois;
        a = param->rc;
        TensionFactor = param->TensionFactor;

        Ecore = 0.5 * TensionFactor * MU;

        extStress[0][0] = param->appliedStress[0];
        extStress[1][1] = param->appliedStress[1];
        extStress[2][2] = param->appliedStress[2];
        extStress[1][2] = param->appliedStress[3];
        extStress[2][0] = param->appliedStress[4];
        extStress[0][1] = param->appliedStress[5];
        extStress[2][1] = extStress[1][2];
        extStress[0][2] = extStress[2][0];
        extStress[1][0] = extStress[0][1];

/*
 *      If the segment is zero length (which can occur when testing whether
 *      multinodes should be split or not), set energy to zero and return.
 */
        dx = x2 - x1;
        dy = y2 - y1;
        dz = z2 - z1;

        if ((dx*dx + dy*dy + dz*dz) < 1.0e-06) {
           return 0.0;
        }
        
        ZImage(param, &dx, &dy, &dz);

        x2 = x1 + dx;
        y2 = y1 + dy;
        z2 = z1 + dz;

        WSelf = SelfEnergy(1, MU, NU, bx, by, bz, x1, y1, z1, x2, y2, z2,
			   a, Ecore);

        WPK = ExtPKEnergy(extStress, E, NU, bx, by, bz, x1, y1, z1, x2, y2, z2);

        W = WSelf + WPK;

        return W;
}


#endif
