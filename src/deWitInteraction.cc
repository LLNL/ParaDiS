/***************************************************************************
 * Function	: deWitInteraction
 * Description	: Calculates interaction stresses
 * 10/23/01 M.Rhee
 * 03/06/02 T.Pierce - modified the calling sequence, to pass in the
 *                     the coordinates, already periodically shifted as
 *                     necessary
 * 05/15/03 M.Hiratani - moved conditional statement out of Interpolation
 *                       removed unnecessary cal. for box size  
 * 09/16/03 T.Pierce - Unrolled a bunch of loops in dSegImgStress for 
 *                     optimization
 * 05/29/04 G.Hommes -	Removed obsolete functions calplanestress() and
 *			CalStress_fullimage().
***************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"

/* changed by Wei Cai, 10/16/2002 
   from
     void deWitInteraction(Home_t *home, ...)
   to
     void deWitInteraction(real8 MU, real8 NU, ...)
*/
void deWitInteraction(real8 MU, real8 NU,
                      real8  Sigma[][3],
                      real8  px, real8  py, real8  pz,
                      real8  tp1, real8  tp2, real8  tp3,
                      real8  burgX, real8  burgY, real8  burgZ,
                      real8  *cntx, real8  *cnty, real8  *cntz)
{
      real8  R1,R2,R3, bp1,bp2,bp3,Rmag,rho1,rho2,rho3 ;
      real8  Y1,Y2,Y3, Ymagsq;
      real8  Trpl, bY1,bY2,bY3, btp1,btp2,btp3 ;
      real8  M11,M22,M33,M12,M13,M23 ;
      real8  K11,K22,K33,K12,K13,K23 ; 
      real8  A11,A22,A33,A12,A13,A23 ; 
      real8  S11,S22,S33,S12,S13,S23 ; 
      real8  F1, F2, PP, onepois, Lp, shr, piinv, F0 ; 

      R1= *cntx -px; 
      R2= *cnty -py; 
      R3= *cntz -pz;

      bp1 = burgX ;
      bp2 = burgY ;
      bp3 = burgZ ;

      Rmag = sqrt( R1*R1 + R2*R2 + R3*R3 ) ;

#define FFACTOR_DEWITRMAG 1e-10
/*#define FFACTOR_DEWITRMAG 10*/
      
      if (Rmag <  FFACTOR_DEWITRMAG) Rmag = FFACTOR_DEWITRMAG ;

      Lp = R1*tp1 + R2*tp2 + R3*tp3 ;

      rho1 = R1 - Lp * tp1 ;
      rho2 = R2 - Lp * tp2 ;
      rho3 = R3 - Lp * tp3 ;

      Y1 = R1 + Rmag * tp1 ;
      Y2 = R2 + Rmag * tp2 ;
      Y3 = R3 + Rmag * tp3 ;

    
      Ymagsq = Y1*Y1 + Y2*Y2 + Y3*Y3 ;

       xvector(bp1,bp2,bp3,Y1,Y2,Y3,&bY1,&bY2,&bY3) ;

       Trpl = bY1*tp1 + bY2*tp2 + bY3*tp3 ;

       xvector(bp1,bp2,bp3,tp1,tp2,tp3,&btp1,&btp2,&btp3) ;

       M11 =  btp1 * Y1 ;
       M22 =  btp2 * Y2 ;
       M33 =  btp3 * Y3 ;
       M12 = 0.5*(btp1*Y2 + btp2*Y1) ;
       M13 = 0.5*(btp1*Y3 + btp3*Y1) ;
       M23 = 0.5*(btp2*Y3 + btp3*Y2) ;

       K11 =  bY1 * tp1 ;
       K22 =  bY2 * tp2 ;
       K33 =  bY3 * tp3 ;
       K12 = 0.5*(bY1*tp2 + bY2*tp1) ;
       K13 = 0.5*(bY1*tp3 + bY3*tp1) ;
       K23 = 0.5*(bY2*tp3 + bY3*tp2) ;

       shr = MU;
       onepois = 1./(1.-NU);
       
       piinv = 1./M_PI ;

       if (Ymagsq > FFACTOR_DEWITRMAG/2 ) 
       {
                F2 = 2./Ymagsq ;
                PP = shr * piinv / Ymagsq ;
        }
        else
        {
                F2 = 2./(Ymagsq+FFACTOR_DEWITRMAG/2) ;
                PP = shr * piinv / (Ymagsq+FFACTOR_DEWITRMAG/2) ;
        }

        if (Rmag > FFACTOR_DEWITRMAG/2) 
                F1 = Lp / Rmag ;
        else
                F1 = Lp / (Rmag+FFACTOR_DEWITRMAG/2) ;


        F0 = Trpl*0.5*onepois ;

        A11 = F2 * (2.*rho1*Y1 + Y1*Y1*F1) ;
        S11 = PP*(K11 - onepois*M11 - F0*(1. + tp1*tp1 + A11)) ;

        A22 = F2 * (2.*rho2*Y2 + Y2*Y2*F1) ;
        S22 = PP*(K22 - onepois*M22 - F0*(1. + tp2*tp2 + A22)) ;

        A33 = F2 * (2.*rho3*Y3 + Y3*Y3*F1) ;
        S33 = PP*(K33 - onepois*M33 - F0*(1. + tp3*tp3 + A33)) ;

        A12 = F2 * (rho1*Y2 + rho2*Y1 + Y1*Y2*F1) ;
        S12 = PP*(K12 - onepois * M12 - F0*(tp1*tp2 + A12)) ;

        A13 = F2 * (rho1*Y3 + rho3*Y1 + Y1*Y3*F1) ;
        S13 = PP*(K13 - onepois * M13 - F0*(tp1*tp3 + A13)) ;

        A23 = F2 * (rho2*Y3 + rho3*Y2 + Y2*Y3*F1) ;
        S23 = PP*(K23 - onepois * M23 - F0*(tp2*tp3 + A23)) ;

        Sigma[0][0] = S11 ;
        Sigma[1][1] = S22 ;
        Sigma[2][2] = S33 ;

        Sigma[0][1] = S12 ;
        Sigma[0][2] = S13 ;
        Sigma[1][2] = S23 ;

        Sigma[2][1] = S23 ;
        Sigma[2][0] = S13 ;
        Sigma[1][0] = S12 ;
}

void deWitStress(real8 MU, real8 NU,
                 real8 Sigma[][3],
                 real8 burgX, real8 burgY, real8 burgZ,
                 real8 xA, real8 yA, real8 zA,
                 real8 xB, real8 yB, real8 zB,
                 real8 x0, real8 y0, real8 z0)
{
    int kk, mm;
    real8 tp1, tp2, tp3, tpsiz;
    real8 sigA[3][3], sigB[3][3];
    
    GetUnitVector(1,xA,yA,zA,xB,yB,zB,&tp1,&tp2,&tp3,&tpsiz);
    
/*
 *  If the segment was zero-length, there is no associated stress
 *  so just zero out Sigma and return to the caller.
 */
    if (tpsiz  <= 0.0) {
        Matrix33_Zero(Sigma);
        return;
    }

    deWitInteraction(MU,NU,sigB,xB,yB,zB,
                     tp1,tp2,tp3, burgX,burgY,burgZ, &x0,&y0,&z0) ;
    deWitInteraction(MU,NU,sigA,xA,yA,zA,
                     tp1,tp2,tp3, burgX,burgY,burgZ, &x0,&y0,&z0) ;

    for (kk=0;kk<3;kk++)
        for (mm=0;mm<3;mm++)
            Sigma[mm][kk] = sigB[mm][kk]-sigA[mm][kk];
    
    return;
}

/* indices of 3rd derivatives */
#if 0
static int DER3IND[10][3] = {{0,0,0},{1,1,1},{2,2,2},
                             {0,0,1},{0,0,2},{1,1,0},
                             {1,1,2},{2,2,0},{2,2,1},
                             {0,1,2}};
static int DER3INVIND[3][3][3]={{{0,3,4},{3,5,9},{4,9,7}},
                                {{3,5,9},{5,1,6},{9,6,8}},
                                {{4,9,7},{9,6,8},{7,8,2}}};
#endif

/* PBC stress table */
#ifndef MAXSTRGRID
#define MAXSTRGRID 51
#endif

static int   Offset1, Offset2, Offset3;
static real8 *RIJMTABLE;
static real8 *RIJMPBCTABLE;
static real8 RIJMLX, RIJMLY, RIJMLZ;


void FreeRijm(void)
{
        if (RIJMTABLE) {
            free(RIJMTABLE);
            RIJMTABLE = NULL;
        }

        return;
}

void FreeRijmPBC(void)
{
        if (RIJMPBCTABLE) {
            free(RIJMPBCTABLE);
            RIJMPBCTABLE = NULL;
        }

        return;
}

void ReadRijm(Home_t *home)
{
	FILE	*fp;
	Param_t	*param;
	int	i, j, k, m;
	int	NX, NY, NZ;
	int	NIMGX, NIMGY, NIMGZ, RIJM_Table_Size;
	real8	Lx, Ly, Lz;
	real8	ratiox, ratioy, ratioz;

	real8	rijml[3] = { 0.0 };

	param=home->param;

	Offset1 = 10*MAXSTRGRID*MAXSTRGRID;
	Offset2 = 10*MAXSTRGRID;
	Offset3 = 10;

	RIJM_Table_Size = (MAXSTRGRID*MAXSTRGRID*MAXSTRGRID*10) * sizeof(real8);
	RIJMTABLE = (real8 *)calloc(1, RIJM_Table_Size);
/*
 *	Only processor zero reads in the table.  It will broadcast the
 *	table (and supporting data) to all other processors
 */
	if (home->myDomain == 0) 
    {
		printf("Reading file %s ...\n",param->Rijmfile);
		fp=fopen(param->Rijmfile,"r");
		if(fp==NULL) {
			fprintf(stderr, "ReadRijm file (%s) not found!\n",
				param->Rijmfile);
			exit (-1);
		}

		fscanf(fp,"%d %d %d",&NX,&NY,&NZ);
		fscanf(fp,"%le %le %le",&Lx,&Ly,&Lz);
		fscanf(fp,"%d %d %d",&NIMGX,&NIMGY,&NIMGZ);

		if ((NX>MAXSTRGRID)||(NY>MAXSTRGRID)||(NZ>MAXSTRGRID)) {
			fprintf(stderr,"NX=%d NY=%d NZ=%d MAXSTRGRID=%d\n",
				NX,NY,NZ,MAXSTRGRID);
			Fatal("table size larger than static array");
		}
    
		RIJMLX = Lx;
		RIJMLY = Ly;
		RIJMLZ = Lz;
    
		Lx=param->Lx;
		Ly=param->Ly;
		Lz=param->Lz;

		ratiox=RIJMLX/Lx;
		ratioy=RIJMLY/Ly;
		ratioz=RIJMLZ/Lz;

		if((fabs(ratiox-ratioy)>1e-3)||(fabs(ratioy-ratioz)>1e-3)) {
			printf("Lx=%e Ly=%e Lz=%e\nLX=%e LY=%e LZ=%e\n",
				Lx,Ly,Lz,RIJMLX,RIJMLY,RIJMLZ);
			Fatal("Wrong RIJM Table read in");
		}


		for(i=0;i<NX;i++)
		for(j=0;j<NY;j++)
		for(k=0;k<NZ;k++)
		for(m=0;m<10;m++) {
			fscanf(fp,"%le\n", &RIJMTABLE[i*Offset1+j*Offset2+k*Offset3+m]);
		}

		fclose(fp);

		param->imgstrgrid[0]=NX;
		param->imgstrgrid[1]=NY;
		param->imgstrgrid[2]=NZ;
    
		param->imgstrgrid[3]=NIMGX;
		param->imgstrgrid[4]=NIMGY;
		param->imgstrgrid[5]=NIMGZ;

		rijml[0] = RIJMLX;
		rijml[1] = RIJMLY;
		rijml[2] = RIJMLZ;

		printf("done.\n");
	}

#ifdef PARALLEL
	MPI_Bcast((char *)RIJMTABLE, RIJM_Table_Size, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast((char *)&param->imgstrgrid[0], sizeof(param->imgstrgrid), MPI_CHAR, 0,MPI_COMM_WORLD);
	MPI_Bcast(&rijml[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

	RIJMLX = rijml[0];
	RIJMLY = rijml[1];
	RIJMLZ = rijml[2];

	return;
}

void ReadRijmPBC(Home_t *home)
{
	FILE	*fp;
	Param_t	*param;
	int	i, j, k, m;
	int	NX, NY, NZ;
	int	NIMGX, NIMGY, NIMGZ, RIJM_Table_Size;
	real8	Lx, Ly, Lz;
	real8	ratiox, ratioy, ratioz;
	real8	rijml[3] = { 0.0 };

	param=home->param;

	Offset1 = 10*MAXSTRGRID*MAXSTRGRID;
	Offset2 = 10*MAXSTRGRID;
	Offset3 = 10;

	RIJM_Table_Size = (MAXSTRGRID*MAXSTRGRID*MAXSTRGRID*10) * sizeof(real8);
	RIJMPBCTABLE = (real8 *)calloc(1, RIJM_Table_Size);
/*
 *	Only processor zero reads in the table.  It will broadcast the
 *	table (and supporting data) to all other processors
 */
	if (home->myDomain == 0) 
    {
		printf("Reading file %s ...\n",param->RijmPBCfile);
		fp=fopen(param->RijmPBCfile,"r");
		if(fp==NULL) {
			fprintf(stderr, "ReadRijmPBC file (%s) not found!\n",
				param->RijmPBCfile);
			Fatal("ReadRijmPBC() error");
		}

		fscanf(fp,"%d %d %d",&NX,&NY,&NZ);
		fscanf(fp,"%le %le %le",&Lx,&Ly,&Lz);
		fscanf(fp,"%d %d %d",&NIMGX,&NIMGY,&NIMGZ);
 
		RIJMLX = Lx;
		RIJMLY = Ly;
		RIJMLZ = Lz;

		Lx=param->Lx;
		Ly=param->Ly;
		Lz=param->Lz;

		ratiox=RIJMLX/Lx;
		ratioy=RIJMLY/Ly;
		ratioz=RIJMLZ/Lz;

		if((fabs(ratiox-ratioy)>1e-3)||(fabs(ratioy-ratioz)>1e-3)) {
			printf("Lx=%e Ly=%e Lz=%e\n"
				"LX=%e LY=%e LZ=%e\n",
				Lx,Ly,Lz,RIJMLX,RIJMLY,RIJMLZ);
			Fatal("Wrong RIJMPBC Table read in");
		}
   
		for(i=0;i<NX;i++)
		for(j=0;j<NY;j++)
		for(k=0;k<NZ;k++)
		for(m=0;m<10;m++) {
			fscanf(fp,"%le\n", &RIJMPBCTABLE[i*Offset1+j*Offset2+k*Offset3+m]);
		}

		fclose(fp);

		param->imgstrgrid[0]=NX;
		param->imgstrgrid[1]=NY;
		param->imgstrgrid[2]=NZ;
		    
		param->imgstrgrid[3]=NIMGX;
		param->imgstrgrid[4]=NIMGY;
		param->imgstrgrid[5]=NIMGZ;

		rijml[0] = RIJMLX;
		rijml[1] = RIJMLY;
		rijml[2] = RIJMLZ;

		printf("done.\n");
	}

#ifdef PARALLEL
	MPI_Bcast((char *)RIJMPBCTABLE, RIJM_Table_Size, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast((char *)&param->imgstrgrid[0], sizeof(param->imgstrgrid), MPI_CHAR, 0,MPI_COMM_WORLD);
	MPI_Bcast(&rijml[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

	RIJMLX = rijml[0];
	RIJMLY = rijml[1];
	RIJMLZ = rijml[2];

	return;
}

void InterpolateRijm(Home_t *home, real8 rt[3], real8 Rijmarray[10], int pbc)
{
    Param_t *param;
    real8 Lx, Ly, Lz, alpha, beta, gamma, xmin, ymin, zmin;
    real8 ratiox, ratio2;
    int m, NX, NY, NZ, nx, ny, nz;

    param=home->param;
    NX=param->imgstrgrid[0];
    NY=param->imgstrgrid[1];
    NZ=param->imgstrgrid[2];

    xmin=param->minSideX;
    ymin=param->minSideY; 
    zmin=param->minSideZ; 
    
    Lx=param->Lx;
    Ly=param->Ly;
    Lz=param->Lz;
    ratiox=RIJMLX/Lx;
    ratio2=ratiox*ratiox;
    
    alpha=(rt[0]-xmin)/Lx*(NX-1);
    beta =(rt[1]-ymin)/Ly*(NY-1);
    gamma=(rt[2]-zmin)/Lz*(NZ-1);

    while(alpha<    0) alpha+=(NX-1);
    while(alpha>=NX-1) alpha-=(NX-1);
    while(beta <    0) beta +=(NY-1);
    while(beta >=NY-1) beta -=(NY-1);
    while(gamma<    0) gamma+=(NZ-1);
    while(gamma>=NZ-1) gamma-=(NZ-1);

#ifdef _CUBICINTERPOLATE
    nx=(int)floor(alpha);
    ny=(int)floor(beta );
    nz=(int)floor(gamma);

    alpha-=nx;
    beta -=ny;
    gamma-=nz;

    if (pbc)
    for(m=0;m<10;m++)
        Rijmarray[m]=
             (1-alpha)*(1-beta)*(1-gamma)*RIJMPBCTABLE[nx*Offset1+ny*Offset2+nz*Offset3+m]
            +(1-alpha)*(1-beta)*(  gamma)*RIJMPBCTABLE[nx*Offset1+ny*Offset2+(nz+1)*Offset3+m]   
            +(1-alpha)*(  beta)*(1-gamma)*RIJMPBCTABLE[nx*Offset1+(ny+1)*Offset2+nz*Offset3+m]    
            +(1-alpha)*(  beta)*(  gamma)*RIJMPBCTABLE[nx*Offset1+(ny+1)*Offset2+(nz+1)*Offset3+m]    
            +(  alpha)*(1-beta)*(1-gamma)*RIJMPBCTABLE[(nx+1)*Offset1+ny*Offset2+nz*Offset3+m]   
            +(  alpha)*(1-beta)*(  gamma)*RIJMPBCTABLE[(nx+1)*Offset1+ny*Offset2+(nz+1)*Offset3+m]
            +(  alpha)*(  beta)*(1-gamma)*RIJMPBCTABLE[(nx+1)*Offset1+(ny+1)*Offset2+nz*Offset3+m]    
            +(  alpha)*(  beta)*(  gamma)*RIJMPBCTABLE[(nx+1)*Offset1+(ny+1)*Offset2+(nz+1)*Offset3+m];
    else
    for(m=0;m<10;m++)
        Rijmarray[m]=
             (1-alpha)*(1-beta)*(1-gamma)*RIJMTABLE[nx*Offset1+ny*Offset2+nz*Offset3+m]
            +(1-alpha)*(1-beta)*(  gamma)*RIJMTABLE[nx*Offset1+ny*Offset2+(nz+1)*Offset3+m]   
            +(1-alpha)*(  beta)*(1-gamma)*RIJMTABLE[nx*Offset1+(ny+1)*Offset2+nz*Offset3+m]    
            +(1-alpha)*(  beta)*(  gamma)*RIJMTABLE[nx*Offset1+(ny+1)*Offset2+(nz+1)*Offset3+m]    
            +(  alpha)*(1-beta)*(1-gamma)*RIJMTABLE[(nx+1)*Offset1+ny*Offset2+nz*Offset3+m]   
            +(  alpha)*(1-beta)*(  gamma)*RIJMTABLE[(nx+1)*Offset1+ny*Offset2+(nz+1)*Offset3+m]
            +(  alpha)*(  beta)*(1-gamma)*RIJMTABLE[(nx+1)*Offset1+(ny+1)*Offset2+nz*Offset3+m]    
            +(  alpha)*(  beta)*(  gamma)*RIJMTABLE[(nx+1)*Offset1+(ny+1)*Offset2+(nz+1)*Offset3+m];

#else
    nx=(int)rint(alpha);
    ny=(int)rint(beta );
    nz=(int)rint(gamma);

    if (pbc)
    for(m=0;m<10;m++)
            Rijmarray[m]=RIJMPBCTABLE[nx*Offset1+ny*Offset2+nz*Offset3+m];
    else
    for(m=0;m<10;m++)
            Rijmarray[m]=RIJMTABLE[nx*Offset1+ny*Offset2+nz*Offset3+m];

#endif
	/* scale stress */
	for(m=0;m<10;m++)
	Rijmarray[m]*=ratio2;

}

#if 0
void dSegStress(Home_t *home,
                real8  Sigma[][3],
                real8  px,    real8 py,    real8 pz,
                real8  dlx,   real8 dly,   real8 dlz,
                real8  burgX, real8 burgY, real8 burgZ,
                real8  rx,    real8 ry,    real8 rz)
{
    /* calculate stress due to a differential dislocation segment
     * Sigma[3][3] return stress
     * (px, py, pz) location of segment (center)
     * (dlx, dly, dlz) vector of segment length
     * (burgX, burgY, burgZ) Burgers vector
     * (rx, ry, rz) location of field point
     *
     * Wei Cai, 7/6/2002
     */
    real8 rt[3], Rmpp[3], dl[3], bxdl[3], bxRmpp[3];
    real8 Rijmarray[10], Rijm[3][3][3], Rijmterm[3][3][3];
    real8 R, R2, R3;
    real8 MUover4pi, onepois;
    int i, j, m, k;

    MUover4pi = home->param->shearModulus/4/M_PI ;
    onepois = 1./(1.- home->param->pois) ;

    dl[0]=dlx; dl[1]=dly; dl[2]=dlz;
    rt[0]=rx-px; rt[1]=ry-py; rt[2]=rz-pz;
    R2=rt[0]*rt[0]+rt[1]*rt[1]+rt[2]*rt[2];
    R=sqrt(R2); R3=R*R2;

    CalRijm(rt, Rijmarray);
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            for(m=0;m<3;m++)
            {
                k=DER3INVIND[i][j][m];
                Rijm[i][j][m]=Rijmarray[k];
            }

    for(m=0;m<3;m++) Rmpp[m]=Rijm[0][0][m]+Rijm[1][1][m]+Rijm[2][2][m];

    bxdl[2] = burgX*dly - burgY*dlx;
    bxdl[0] = burgY*dlz - burgZ*dly;
    bxdl[1] = burgZ*dlx - burgX*dlz;
    
    bxRmpp[2] = burgX*Rmpp[1] - burgY*Rmpp[0];
    bxRmpp[0] = burgY*Rmpp[2] - burgZ*Rmpp[1];
    bxRmpp[1] = burgZ*Rmpp[0] - burgX*Rmpp[2];

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            Sigma[i][j] = -0.5*(bxRmpp[i]*dl[j]+bxRmpp[j]*dl[i]);

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            for(m=0;m<3;m++)
            {
                Rijmterm[i][j][m]=Rijm[i][j][m]-Rmpp[m]*(i==j);
            }
    
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
        {
            for(m=0;m<3;m++)
                Sigma[i][j]+=onepois*Rijmterm[i][j][m]*bxdl[m];
            Sigma[i][j]*=MUover4pi;
        }
}
#endif

void dSegImgStress(Home_t *home,
                   real8  Sigma[3][3],
                   real8  px,    real8 py,    real8 pz,
                   real8  dlx,   real8 dly,   real8 dlz,
                   real8  burgX, real8 burgY, real8 burgZ,
                   real8  rx,    real8 ry,    real8 rz,
                   int pbc)
{
    /* calculate image stress due to a differential dislocation segment
     * using the table Rijm
     * Sigma[3][3] return stress
     * (px, py, pz) location of segment (center)
     * (dlx, dly, dlz) vector of segment length
     * (burgX, burgY, burgZ) Burgers vector
     * (rx, ry, rz) location of field point
     *
     * Wei Cai, 7/6/2002
     */
    real8 rt[3], Rmpp[3], dl[3], bxdl[3], bxRmpp[3];
/*    real8 R, R2, R3;  */
    real8 MUover4pi, onepois;

    MUover4pi = home->param->shearModulus/4/M_PI ;
    onepois = 1./(1.- home->param->pois) ;

    dl[0]=dlx; dl[1]=dly; dl[2]=dlz;
    rt[0]=rx-px; rt[1]=ry-py; rt[2]=rz-pz;
/*    R2=rt[0]*rt[0]+rt[1]*rt[1]+rt[2]*rt[2];
      R=sqrt(R2); R3=R*R2;
*/
    real8 Rijmarray[10] = { 0.0 };

#if 0 /* instead of computing Rijmarray for the primary dislocation */
    CalRijm(rt, Rijmarray);
#else /* use the Rijm table for image sources */
    InterpolateRijm(home, rt, Rijmarray, pbc);
#endif
    
    real8 Rijm[3][3][3] = {{{ 0.0 }}};
/*
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++) {
            for(int m=0;m<3;m++)
            {
                k=DER3INVIND[i][j][m];
                Rijm[i][j][m]=Rijmarray[k];
            }
        }
*/

/* unroll above loop */

/*
    Rijm[0][0][0] = Rijmarray[0] ;
    Rijm[1][1][1] = Rijmarray[1] ;
    Rijm[2][2][2] = Rijmarray[2] ;
    Rijm[0][0][1] = Rijm[0][1][0] = Rijm[1][0][0] = Rijmarray[3] ;
    Rijm[0][0][2] = Rijm[0][2][0] = Rijm[2][0][0] = Rijmarray[4] ;
    Rijm[0][1][1] = Rijm[1][0][1] = Rijm[1][1][0] = Rijmarray[5] ;
    Rijm[1][1][2] = Rijm[1][2][1] = Rijm[2][1][1] = Rijmarray[6] ;
    Rijm[0][2][2] = Rijm[2][0][2] = Rijm[2][2][0] = Rijmarray[7] ;
    Rijm[1][2][2] = Rijm[2][1][2] = Rijm[2][2][1] = Rijmarray[8] ;
    Rijm[0][1][2] = Rijm[1][2][0] = Rijm[2][0][1] =
    Rijm[2][1][0] = Rijm[1][0][2] = Rijm[0][2][1] = Rijmarray[9] ;
*/

    Rijm[0][0][0] = Rijmarray[0] ;
    Rijm[1][1][1] = Rijmarray[1] ;
    Rijm[2][2][2] = Rijmarray[2] ;
    Rijm[0][0][1] = Rijmarray[3] ;
    Rijm[0][1][0] = Rijmarray[3] ;
    Rijm[1][0][0] = Rijmarray[3] ;
    Rijm[0][0][2] = Rijmarray[4] ;
    Rijm[0][2][0] = Rijmarray[4] ;
    Rijm[2][0][0] = Rijmarray[4] ;
    Rijm[0][1][1] = Rijmarray[5] ;
    Rijm[1][0][1] = Rijmarray[5] ;
    Rijm[1][1][0] = Rijmarray[5] ;
    Rijm[1][1][2] = Rijmarray[6] ;
    Rijm[1][2][1] = Rijmarray[6] ;
    Rijm[2][1][1] = Rijmarray[6] ;
    Rijm[0][2][2] = Rijmarray[7] ;
    Rijm[2][0][2] = Rijmarray[7] ;
    Rijm[2][2][0] = Rijmarray[7] ;
    Rijm[1][2][2] = Rijmarray[8] ;
    Rijm[2][1][2] = Rijmarray[8] ;
    Rijm[2][2][1] = Rijmarray[8] ;
    Rijm[0][1][2] = Rijmarray[9] ;
    Rijm[1][2][0] = Rijmarray[9] ;
    Rijm[2][0][1] = Rijmarray[9] ;
    Rijm[2][1][0] = Rijmarray[9] ;
    Rijm[1][0][2] = Rijmarray[9] ;
    Rijm[0][2][1] = Rijmarray[9] ;

/*     for(m=0;m<3;m++) Rmpp[m]=Rijm[0][0][m]+Rijm[1][1][m]+Rijm[2][2][m]; */

    Rmpp[0]=Rijm[0][0][0]+Rijm[1][1][0]+Rijm[2][2][0];
    Rmpp[1]=Rijm[0][0][1]+Rijm[1][1][1]+Rijm[2][2][1];
    Rmpp[2]=Rijm[0][0][2]+Rijm[1][1][2]+Rijm[2][2][2];

    bxdl[2] = burgX*dly - burgY*dlx;
    bxdl[0] = burgY*dlz - burgZ*dly;
    bxdl[1] = burgZ*dlx - burgX*dlz;
    
    bxRmpp[2] = burgX*Rmpp[1] - burgY*Rmpp[0];
    bxRmpp[0] = burgY*Rmpp[2] - burgZ*Rmpp[1];
    bxRmpp[1] = burgZ*Rmpp[0] - burgX*Rmpp[2];

/*
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            Sigma[i][j] = -0.5*(bxRmpp[i]*dl[j]+bxRmpp[j]*dl[i]);
*/

/* optimization: unroll */

    Sigma[0][0] =               -1.0*(bxRmpp[0]*dl[0]) ;
    Sigma[0][1] = Sigma[1][0] = -0.5*(bxRmpp[0]*dl[1]+bxRmpp[1]*dl[0]);
    Sigma[0][2] = Sigma[2][0] = -0.5*(bxRmpp[0]*dl[2]+bxRmpp[2]*dl[0]);
    Sigma[1][1] =               -1.0*(bxRmpp[1]*dl[1]) ;
    Sigma[1][2] = Sigma[2][1] = -0.5*(bxRmpp[1]*dl[2]+bxRmpp[2]*dl[1]);
    Sigma[2][2] =               -1.0*(bxRmpp[2]*dl[2]) ;


/*
    real8 Rijmterm[3][3][3];

    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            for(int m=0;m<3;m++)
            {
                Rijmterm[i][j][m]=Rijm[i][j][m]-Rmpp[m]*(i==j);
            }
*/

/* optimization change */

/*

    for(int i=0;i<3;i++)
        for(int m=0;m<3;m++)
        {
            Rijm[i][i][m] -= Rmpp[m] ;
        }
*/

    Rijm[0][0][0] -= Rmpp[0] ;
    Rijm[0][0][1] -= Rmpp[1] ;
    Rijm[0][0][2] -= Rmpp[2] ;
    
    Rijm[1][1][0] -= Rmpp[0] ;
    Rijm[1][1][1] -= Rmpp[1] ;
    Rijm[1][1][2] -= Rmpp[2] ;
    
    Rijm[2][2][0] -= Rmpp[0] ;
    Rijm[2][2][1] -= Rmpp[1] ;
    Rijm[2][2][2] -= Rmpp[2] ;
    
/*
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
        {
            for(m=0;m<3;m++)
**                Sigma[i][j]+=onepois*Rijmterm[i][j][m]*bxdl[m]; **
                Sigma[i][j]+=onepois*Rijm[i][j][m]*bxdl[m];
            Sigma[i][j]*=MUover4pi;
        }
*/

    Sigma[0][0] = MUover4pi*(Sigma[0][0]+onepois*(Rijm[0][0][0]*bxdl[0]+
                                                  Rijm[0][0][1]*bxdl[1]+
                                                  Rijm[0][0][2]*bxdl[2])) ;

    Sigma[0][1] = MUover4pi*(Sigma[0][1]+onepois*(Rijm[0][1][0]*bxdl[0]+
                                                  Rijm[0][1][1]*bxdl[1]+
                                                  Rijm[0][1][2]*bxdl[2])) ;

    Sigma[0][2] = MUover4pi*(Sigma[0][2]+onepois*(Rijm[0][2][0]*bxdl[0]+
                                                  Rijm[0][2][1]*bxdl[1]+
                                                  Rijm[0][2][2]*bxdl[2])) ;

    Sigma[1][1] = MUover4pi*(Sigma[1][1]+onepois*(Rijm[1][1][0]*bxdl[0]+
                                                  Rijm[1][1][1]*bxdl[1]+
                                                  Rijm[1][1][2]*bxdl[2])) ;

    Sigma[1][2] = MUover4pi*(Sigma[1][2]+onepois*(Rijm[1][2][0]*bxdl[0]+
                                                  Rijm[1][2][1]*bxdl[1]+
                                                  Rijm[1][2][2]*bxdl[2])) ;

    Sigma[2][2] = MUover4pi*(Sigma[2][2]+onepois*(Rijm[2][2][0]*bxdl[0]+
                                                  Rijm[2][2][1]*bxdl[1]+
                                                  Rijm[2][2][2]*bxdl[2])) ;

    Sigma[1][0] = Sigma[0][1] ;
    Sigma[2][0] = Sigma[0][2] ;
    Sigma[2][1] = Sigma[1][2] ;
}



#if 0
void testdeWitStress(Home_t *home)
{/* compute stress due to a finite segment on a plane
  * compare results with matlab script testcalstress.m
  * agreement reached, relative accuracy within 5e-7
  * Wei Cai, 7/5/2002
  */
    FILE *fp12, *fp23, *fp31;
    FILE *fp11, *fp22, *fp33;
    int NX, NY, i, j, ii, jj;
    real8 Lx, Ly, x, y, z;
    real8 burgX, burgY, burgZ;
    real8 sigB[3][3], sigA[3][3], str[3][3];
    real8 xA, yA, zA, xB, yB, zB, tpx, tpy, tpz, tpsize;
    real8 x0, y0, z0, dlx, dly, dlz;
    real8 MU, NU;

    MU=home->param->shearModulus;
    NU=home->param->pois;
    
    NX=41; NY=41;
    /* Lx=20; Ly=20; */
    Lx=80; Ly=19.5959;
    burgX=5; burgY=4; burgZ=3;
    xA=-1; yA=-2; zA=-4;
    xB= 1; yB= 2; zB= 4;
    x0=(xA+xB)/2; y0=(yA+yB)/2; z0=(zA+zB)/2;
    dlx=xB-xA; dly=yB-yA; dlz=zB-zA;
    tpsize=sqrt(dlx*dlx+dly*dly+dlz*dlz);
    tpx=dlx/tpsize;tpy=dly/tpsize;tpz=dlz/tpsize;
    
    fp12=fopen("dewitstr12.out","w");
    fp23=fopen("dewitstr23.out","w");
    fp31=fopen("dewitstr31.out","w");
    fp11=fopen("dewitstr11.out","w");
    fp22=fopen("dewitstr22.out","w");
    fp33=fopen("dewitstr33.out","w");

    for(i=0;i<NY;i++)
    {
        for(j=0;j<NX;j++)
        {
            x=-Lx/2+Lx*j/(NX-1);
            y=-Ly/2+Ly*i/(NY-1);
            z=5;
#if 1 /* if test deWitInteraction (done) */
#if 0
            deWitInteraction(MU,NU,sigB,xB,yB,zB,
                             tpx,tpy,tpz, burgX,burgY,burgZ,
                             &x,&y,&z) ;
            deWitInteraction(MU,NU,sigA,xA,yA,zA,
                             tpx,tpy,tpz, burgX,burgY,burgZ,
                             &x,&y,&z) ;
            for(ii=0;ii<3;ii++) for(jj=0;jj<3;jj++)
                str[ii][jj]=sigB[ii][jj]-sigA[ii][jj];
#else
            deWitStress(MU,NU,str,burgX,burgY,burgZ,
                        xA,yA,zA, xB,yB,zB, x,y,z);
#endif
#else /* else test dSegStress (stress of differential segment) */
            /*
            dSegStress(home,str,x0,y0,z0,dlx,dly,dlz,
                       burgX,burgY,burgZ,x,y,z);
            */
            dSegImgStress(home,str,x0,y0,z0,dlx,dly,dlz,
                          burgX,burgY,burgZ,x,y,z);
#endif
            fprintf(fp12,"%e ",str[0][1]);
            fprintf(fp23,"%e ",str[1][2]);
            fprintf(fp31,"%e ",str[2][0]);
            fprintf(fp11,"%e ",str[0][0]);
            fprintf(fp22,"%e ",str[1][1]);
            fprintf(fp33,"%e ",str[2][2]);
            /*
            printf("%e ",str[0][1]);
            printf("%e ",str[1][2]);
            printf("%e ",str[2][0]);
            printf("%e ",str[0][0]);
            printf("%e ",str[1][1]);
            printf("%e ",str[2][2]);
            */
        }
        fprintf(fp12,"\n");
        fprintf(fp23,"\n");
        fprintf(fp31,"\n");
        fprintf(fp11,"\n");
        fprintf(fp22,"\n");
        fprintf(fp33,"\n");
    }
    
    fclose(fp12);
    fclose(fp23);
    fclose(fp31);
    fclose(fp11);
    fclose(fp22);
    fclose(fp33);
}
#endif


#if 0
void testdeWitStress2()
{/* compute stress for benchmark
  * Wei Cai, 10/16/2002
  */

    real8 MU, NU;
    real8 xA, yA, zA, xB, yB, zB, bx, by, bz, x0, y0, z0;
    real8 stress[3][3];
    int i;
    
    real8 data[14][12]={
        /* xA, yA, zA,    xB, yB, zB,       bx, by, bz,      x0,  y0,  z0 */
        {243, 729, 486, 2673, 3159, 2916, 1.65, 1.65, 1.65, 2916, 3402, 3159},
        {243, 729, 486, 2673, 3159, 2916, 1.65, 1.65, 1.65, 2916, 3402, 2430},
        {243, 729, 486, 2673, 3159, 2916, 1.65, 1.65, -3.30, 2916, 3402, 3159},
        {243, 729, 486, 2673, 3159, 2916, 1.65, 1.65, -3.30, 2916, 3402, 2430},
        {243, 729, 486, 2673, 3159, 2916, 3.30, 3.30, -1.65, 2916, 3402, 3159},
        {243, 729, 486, 2673, 3159, 2916, 3.30, 3.30, -1.65, 2916, 3402, 2430},
        {243, 729, 486, 2673, 3159, 2916, 3.30, 3.30, -1.65, 1701, 2187, 1215},
        {243, 729, 486, 2673, 3159, 2916, 1.65, 1.65, 1.65, 3159, 3159, 3159}, 
        {243, 729, 486, 2673, 3159, 2916, 1.65, 1.65, -3.3, 3159, 3159, 3159},
        {243, 729, 486, 2916, 3402, 2430, 1.65, 1.65, -3.3, 2673, 3159, 2916},
        {0, 0, 0,  1,   0,    0,    3.30, 3.30, -1.65, 0,   0,   1.5},
        {0, 0, 0,  1,   0,    0,    3.30, 3.30, -1.65, 0,   0,   0.5},
        {1010, 1020, 1030, 10, 20, 30, 1.65, 1.65, -3.3, 1110, 1120, 1130},
        {1010, 1020, 1030, 10, 20, 30, 1.65, 1.65, -3.3, 1310, 1120, 930 }
    };    

    MU=6.85e10; /* Pa */
    NU=0.35;

    printf("        sig11,        sig22,        sig33,"
           "        sig12,        sig23,        sig31\n");
    
    for(i=0;i<14;i++)
    {
        xA=data[i][0];
        yA=data[i][1];
        zA=data[i][2];
        xB=data[i][3];
        yB=data[i][4];
        zB=data[i][5];
        bx=data[i][6];
        by=data[i][7];
        bz=data[i][8];
        x0=data[i][9];
        y0=data[i][10];
        z0=data[i][11];
#if 0 /* perturb it a bit */
        xA+=drand48()*0;
        yA+=drand48()*0;
        zA+=drand48()*0;
        xB+=drand48()*0;
        yB+=drand48()*0;
        zB+=drand48()*0;
        bx+=drand48()*0.01;
        by+=drand48()*0.01;
        bz+=drand48()*0.01;
        x0+=drand48()*0;
        y0+=drand48()*0;
        z0+=drand48()*0;
#endif
        deWitStress(MU,NU,stress,bx,by,bz,
                    xA,yA,zA, xB,yB,zB, x0,y0,z0);

        printf("%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
               stress[0][0],stress[1][1],stress[2][2],
               stress[0][1],stress[1][2],stress[2][0]);
    }
}
#endif


#if 0
void deWitStress_(Home_t *home,
                           real8 MU, real8 NU,
                           real8 Sigma[][3],
                           real8 burgX, real8 burgY, real8 burgZ,
                           real8 xA, real8 yA, real8 zA,
                           real8 xB, real8 yB, real8 zB,
                           real8 x0, real8 y0, real8 z0)
{
    Param_t *param;
    real8 Lx, Ly, Lz, xmin, ymin, zmin, xmax, ymax, zmax;
    int mm, kk;
    real8 rc, rc2, xm2, ym2, zm2, xm, ym, zm;
    real8 gsig[3][3];
    real8 dx, dy, dz, dx2, dy2, dz2, rMid2;
    
    param=home->param;
    Lx=param->Lx;
    Ly=param->Ly;
    Lz=param->Lz;
    
    rc=Lx; if(rc>Ly) rc=Ly; if(rc>Lz) rc=Lz; rc*=0.45; 
    rc2=rc*rc;
    
    dx = xB - xA;
    dy = yB - yA;
    dz = zB - zA;
    
    ZImage (param, &dx, &dy, &dz) ;
    xA = xB - dx ; yA = yB - dy ; zA = zB - dz ;
    
    /* calc the distance between midpoints of two segments, and apply cutoff */
    xm2 = xA + dx*0.5 ; ym2 = yA + dy*0.5 ; zm2 = zA + dz*0.5 ;

    dx2 = xm2-x0 ; dy2 = ym2-y0 ; dz2 = zm2-z0 ;
    ZImage (param, &dx2, &dy2, &dz2) ;

    xm = xm2 - dx2;  ym = ym2 - dy2;  zm = zm2 - dz2;

    rMid2 = dx2*dx2+dy2*dy2+dz2*dz2;
    if(rMid2<rc2)
    {
        deWitStress(MU, NU, Sigma, burgX, burgY, burgZ,
                    xA, yA, zA, xB, yB, zB, xm, ym, zm);
        
        /* tabulated image stress (image only) */
        /* need to double-check sign (!!) */
        dSegImgStress(home,gsig,xm2,ym2,zm2,dx,dy,dz,
                      burgX,burgY,burgZ, xm,ym,zm, 0);
        for (mm=0;mm<3;mm++) for (kk=0;kk<3;kk++)
            Sigma[mm][kk] += gsig[mm][kk];
    }
    else
    {
        /* tabulated image stress (periodic) */
        /* need to double-check sign (!!) */
        dSegImgStress(home,Sigma,xm2,ym2,zm2,dx,dy,dz,
                      burgX,burgY,burgZ, xm,ym,zm, 1);
    }
}
#endif


#if 0
void CalStress(Home_t *home,
               real8 xm,real8 ym,real8 zm,
               Node_t *node1, Node_t *node2,
               real8 totremSig[][3])
{/* calculate stress at point (xm, ym, zm) due to all dislocation segments
    and their images, excluding segments connecting to node1 and node2 */
   Param_t *param;
   real8 Lx, Ly, Lz, xmin, ymin, zmin, xmax, ymax, zmax;
   Node_t *rNodeA, *rNodeB ;
   real8 x1,y1,z1,x2,y2,z2,dx,dy,dz,xA,yA,zA,xB,yB,zB;
   real8 sigma[3][3];
   int i, mm, kk, nc2, ti2;
   real8 MU, NU, burgX, burgY, burgZ;
   
   param=home->param;
    Lx=param->Lx;
    Ly=param->Ly;
    Lz=param->Lz;

   MU=param->shearModulus;
   NU=param->pois;
   
   /* initialize stress to zero */
   for (mm=0;mm<3;mm++)
       for (kk=0;kk<3;kk++)
           totremSig[mm][kk] =0;
               
   /* loop through all segments */
   for(i=0; i<home->newNodeKeyPtr; i++)
   {
       rNodeB = home->nodeKeys[i];
       if(!rNodeB) continue;

       xB = rNodeB->x;
       yB = rNodeB->y;
       zB = rNodeB->z;

       nc2 = rNodeB->numNbrs ;

       for (ti2=0 ; ti2<nc2; ti2++) {
           
           rNodeA = GetNeighborNode(home,rNodeB,ti2);
           if (!rNodeA) continue ;
               
           if (OrderNodes (rNodeA,rNodeB)!=1) continue;

           dx = xB - rNodeA->x ;
           dy = yB - rNodeA->y ;
           dz = zB - rNodeA->z ;
           /* (PBC disabled for this function)
             ZImage (param, &dx, &dy, &dz) ;
           */
           xA = xB - dx ; yA = yB - dy ; zA = zB - dz ;
 
           burgX= rNodeB->burgX[ti2];
           burgY= rNodeB->burgY[ti2];
           burgZ= rNodeB->burgZ[ti2];
           
           deWitStress(MU, NU,
                       sigma, burgX, burgY, burgZ,
                       xB, yB, zB,
                       xA, yA, zA,
                       xm, ym, zm);
           for (mm=0;mm<3;mm++)
               for (kk=0;kk<3;kk++)
                   totremSig[mm][kk] += sigma[mm][kk];
           
       }/* end for(iNbr2...) */
   }/* end for(i=0;...) */
}
#endif
