/***************************************************************************
 *
 *      Author:       Tomas Oppelstrup
 *      Module:       CorrectionTable.c
 *      Description:  Contains functions for creating a table used
 *                    in conjuction with the Fast Multipole functions
 *                    to adjust the calculated stress to allow for
 *                    periodic images of the problem space.
 *
 *      Includes public functions:
 *
 *          fmsigma()
 *          CreateCorrectionTable()
 *          CorrectionTableInit()
 *          DoTableCorrection()
 *          FreeCorrectionTable()
 *
 *      Includes private functions:
 *
 *          countit()
 *          countit_full()
 *          paxby()
 *
 ***************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "Home.h"
#include "FM.h"

#ifdef PARALLEL
#include  <mpi.h>
#endif

#ifdef TAYLORFMM

#define CTABLE_VERSION_STR "version="
#define CTABLE_VERSION     2

static int   ctab_norder = 0;
static int   ctab_uorder = 0;
static int   ctab_levels = 0;
static real8 ctab_MU     = 0.0;
static real8 ctab_NU     = 0.0;
static real8 ctab_boxl   = 0.0;

static real8 *correctionTbl = NULL;

//#define DO_ALL_LEVELS 1

/*  EVALUATE STRESS RESULTING FROM MULTIPOLE EXPANSION
    Given material constans mu,nu, and a multipole expansion Eeta, or order norder,
    evaluate the stress tensor at point R. R is the vector from the center of expansion
    to the evaluation point. The resulting stress is written into sigmatot.

    Eeta should be an array resulting from calls to makeeta, and should have
    9*norder*(noder+1)*(norder+2)/6 elements (this is always an integer, and
    evaluated as written, yields the correct result with integer arithmetic).
*/
static void fmsigma(real8 mu, real8 nu, int norder,
	     real8 *R, real8 *Eeta, real8 sigmatot[3][3]) {
  int cyc[] = {0,1,2,0,1,2};
  real8 sigma[3][3],Gijmkl[3][3][3],Gmppkl[3][3][3];
  int iorder,i,j,k,l,m,n,a,nx,ny,nz,pows[NTMAX][3],npows,nidx,idx[NTMAX];
  int etaoff,etadis,etadis0;
  real8 t,t2,t11,t12,t3,terms[NTMAX],tvec[NTMAX];

  real8 pi = 3.1415926535897932385, mu8pi = mu/(8.0*pi), two1nu = 2.0/(1.0-nu);  

  /* dm stuff */
  int b,c,mparity;
  real8 f,ri,gx,gy,gz,rfact,ca,cb,cc,gtab[NMAX+1][3];

  static real8 fact[NMAX+1],ifact[NMAX+1],dfact[2*NMAX-3+2],qtab[NMAX+1][NMAX+1];
  static int inited = 0;
  if(inited == 0) {
    makeqtab(qtab);
    makeftabs(fact,ifact,dfact);
    inited = 1;
  }

  for(i = 0; i<3; i++)
    for(j = 0; j<3; j++)
      sigmatot[i][j] = 0.0;

  /* Set up for dm */
  ri = 0.0; for(i = 0; i<3; i++) ri = ri + R[i]*R[i];
  ri = 1.0/sqrt(ri);
  gx = R[0]*ri; gy = R[1]*ri; gz = R[2]*ri;
  rfact = -ri;
  
  gtab[0][0] = 1.0; gtab[0][1] = 1.0; gtab[0][2] = 1.0;
  for(i = 1; i<=norder+3; i++) {
    gtab[i][0] = gtab[i-1][0] * gx;
    gtab[i][1] = gtab[i-1][1] * gy;
    gtab[i][2] = gtab[i-1][2] * gz;
  }

  for(iorder=0; iorder<=norder; iorder++) {
    /* Here comes some hard core inlining of dm */
    /*dm(iorder+3,R,terms,pows);*/
    n = iorder+3;
    rfact = -rfact*ri;

    i = 0;
    for(nz=0; nz<=n; nz++) {
      for(ny=0; ny<=n-nz; ny++) {
	nx = n-nz-ny;
	pows[i][0] = nx;
	pows[i][1] = ny;
	pows[i][2] = nz;
	
	terms[i] = 0.0;
	ca = 1.0;
	for(a=0; a<=nx; a+=2) {
	  cb = 1.0;
	  for(b=0; b<=ny; b+=2) {
	    cc = 1.0;
	    m = (a+b)/2;
	    mparity = 1-((m&1)<<1);
	    for(c=0; c<=nz; c+=2) {
	      f = mparity*dfact[2*(n-m)-3+1] * 
		gtab[nx-a][0]*gtab[ny-b][1]*gtab[nz-c][2];
	      
	      f = f * ca*cb*cc;
	      terms[i] = terms[i] + f;
	      /*cc = cc * (nz-c)*(nz-c-1)/(c+2);*/
	      cc = cc * qtab[nz-c][c];
	      m = m+1;
	      mparity = -mparity;
	    }
	    /*cb = cb * (ny-b)*(ny-b-1)/(b+2);*/
	    cb = cb * qtab[ny-b][b];
	  }
	  /*ca = ca * (nx-a)*(nx-a-1)/(a+2);*/
	  ca = ca * qtab[nx-a][a];
	}
	terms[i] = terms[i] * fact[n]*ifact[nx]*ifact[ny]*ifact[nz] * rfact;
	i = i+1;
      }
    }
    /* End of inlining */

    npows = (iorder+1+3)*(iorder+2+3)/2;
    for(k = 0; k<npows; k++) {
      terms[k] = terms[k] * ifact[iorder+3] *
	fact[pows[k][0]]*fact[pows[k][1]]*fact[pows[k][2]];
    }

    for(m = 0; m<3; m++)
      for(k = 0; k<3; k++)
	for(l = 0; l<3; l++)
	  Gmppkl[m][k][l] = 0.0;
    
    etadis = (iorder+1)*(iorder+2) >> 1; /* (iorder+1)*(iorder+2)/2          */
    etaoff = 3*iorder*etadis;            /* 9*iorder*(iorder+1)*(iorder+2)/6 */
    for(i=0; i<3; i++)
      for(j=0; j<3; j++) {
	for(m=0; m<3; m++) {
	  /* From 'terms', reconstruct R(i,j,m,a1,...,aq) */
	  int d[] = {0,0,0};
	  d[i]++; d[j]++; d[m]++;
	  nidx = 0;
	  for(k = 0; k<npows; k++) {
	    nx = pows[k][0] - d[0];
	    ny = pows[k][1] - d[1];
	    nz = pows[k][2] - d[2];
	    if(nx>=0 && ny>=0 && nz>=0) {
	      tvec[nidx] = terms[k] * fact[iorder] *
		ifact[nx]*ifact[ny]*ifact[nz];
	      idx[nidx++] = k;
	    }
	  }

	  /* Gijmkl = R(i,j,m,a1,...,aq)*Eeta(k,l,a1,...,aq), sum(a1,...,aq)
	     Eeta(i,j,a1,...,aq) = E(i,j)*eta(a1,...,aq) */
	  etadis0 = etaoff;
	  for(k = 0; k<3; k++)
	    for(l = 0; l<3; l++) {
	      t = 0.0;
	      for(a = 0; a<nidx; a++) {
		/* t = t + tvec[a]*Eeta[iorder][k][l][a]; */
		/* t = t + tvec[a]*Eeta[etaoff + (3*k+l)*etadis + a]; */
		 t = t + tvec[a]*Eeta[etadis0++];
	      }
	      Gijmkl[m][k][l] = t;
	      if(i == j) Gmppkl[m][k][l] = Gmppkl[m][k][l] + t;
	    }
	}
  
	t2 = 0.0;
	for(k = 0; k<3; k++) {
	  /* Second term: 2/(1-nu) * ijk(k,m,n)G(i,j,m,n,k) */
	  m = cyc[k+1]; n = cyc[k+2];
	  t2 = t2 + (Gijmkl[m][n][k] - Gijmkl[n][m][k]);
	}
	sigma[i][j] = two1nu * t2;
      }

    /* Term one: ijk(j,m,n)G(m,p,p,n,i) + ijk(i,m,n)G(m,p,p,n,j) */
    for(i = 0; i<3; i++)
      for(j = 0; j<3; j++) {
	m = cyc[j+1]; n = cyc[j+2];
	t11 = Gmppkl[m][n][i] - Gmppkl[n][m][i];

	m = cyc[i+1]; n = cyc[i+2];
	t12 = Gmppkl[m][n][j] - Gmppkl[n][m][j];
      
	sigma[i][j] = sigma[i][j] + (t11 + t12);
      }

    /* Third term: 2/(1-nu) * delta(i,j)ijk(km,n)G(p,p,m,n,k) */
    t3 = 0.0;
    for(k = 0; k<3; k++) {
      m = cyc[k+1]; n = cyc[k+2];
      t3 = t3 + Gmppkl[m][n][k] - Gmppkl[n][m][k];
    }
    t3 = t3 * 2.0/(1.0-nu);
    for(i = 0; i<3; i++) sigma[i][i] = sigma[i][i] - t3;
    t = ipow(-1.0,iorder);
    for(i = 0; i<3; i++)
      for(j = 0; j<3; j++) {
	sigma[i][j] = sigma[i][j] * mu8pi*ifact[iorder];
	sigmatot[i][j] = sigmatot[i][j] + t*sigma[i][j];
      }
  }
}


static void countit(int norder, int *nlist_p, int **idxlist_p) {
  int i,n,idx,j,k,nalpha;
  int nx,ny,nz,m;
  int nlist = 0, *idxlist = NULL;

  nalpha = 9*(norder+3)*(norder+2)*(norder+1)/6;
  if(idxlist_p != NULL)
    idxlist = (int *) malloc(sizeof(int) * nalpha);
  
  for(i = 0; i<=norder; i++) {
    k = 0;
    n = (i+2)*(i+1)/2;
    idx = 3*n*i;
    for(nz = 0; nz<=i; nz++)
      for(ny = 0; ny<=i-nz; ny++) {
	nx = i-nz-ny;

	m = nz*(i+1) + ny + 1 - (nz-1)*nz/2;
	if(m != k+1)
	  fprintf(stderr,"m-k mismatch: m=%d, k=%d\n",m,k);

	if(nx>=ny && ny>=nz) {
	  for(j = 0; j<9; j++)
	    /*if(j%3 >= j/3)*/ {
	    if(idxlist != NULL) {
	      idxlist[nlist] = idx + j*n + k;
/*
	      printf("nlist=%3d nalpha=%4d listno=%4d %d,%d,%d\n",
		     nlist+1,nalpha,idx+j*n+k,nx,ny,nz);
*/
	    }
	    nlist = nlist+1;
	    }
	}
	k = k + 1;
      }
  }
  if(idxlist_p != NULL) *idxlist_p = idxlist;
  *nlist_p = nlist;
}


static void countit_full(int norder, int *nlist_p, int **idxlist_p) {
  int i,n,idx,j,k,nalpha;
  int nx,ny,nz,m;
  int nlist = 0, *idxlist = NULL;

  nalpha = 9*(norder+3)*(norder+2)*(norder+1)/6;
  if(idxlist_p != NULL)
    idxlist = (int *) malloc(sizeof(int) * nalpha);

  for(i = 0; i<=norder; i++) {
    k = 0;
    n = (i+2)*(i+1)/2;
    idx = 3*n*i;
    for(nz = 0; nz<=i; nz++)
      for(ny = 0; ny<=i-nz; ny++) {
        nx = i-nz-ny;

        m = nz*(i+1) + ny + 1 - (nz-1)*nz/2;
        if(m != k+1)
          fprintf(stderr,"m-k mismatch: m=%d, k=%d\n",m,k);

        /*if(nx>=ny && ny>=nz)*/ { /* Difference between this and compressed table is the commented-out if-statement on this line.*/
          for(j = 0; j<9; j++)
            /*if(j%3 >= j/3)*/ {
            if(idxlist != NULL) {
              idxlist[nlist] = idx + j*n + k;
/*
 *               printf("nlist=%3d nalpha=%4d listno=%4d %d,%d,%d\n",
 *                                    nlist+1,nalpha,idx+j*n+k,nx,ny,nz);
 *                                    */
            }
            nlist = nlist+1;
            }
        }
        k = k + 1;
      }
  }
  if(idxlist_p != NULL) *idxlist_p = idxlist;
  *nlist_p = nlist;
}


#define EPVEC_LEN 6

/*
 *  mpi_np   Number of mpi tasks involved in the calculation
 *  mpi_pid  MPI tasks number of current task
 *  pbc      3 element array of flags indicating whether periodic
 *           boundaries are enabled in each dimension.
 *              pbc[0] == 1 :  periodic in X dimension
 *              pbc[1] == 1 :  periodic in Y dimension
 *              pbc[2] == 1 :  periodic in Z dimension
 */
void Taylor_CreateCorrectionTable(Home_t *home, int numLevels, int pbc[3],
                           int mpi_np, int mpi_pid)
{
  Param_t *param;
  int   nep = EPVEC_LEN;
  int   anisoEnabled = 0;
  int   i,j,k,m,n,ix,iy,iz,nordertab,norder,uorder,nlevels,maxlevels,nalpha,neta;
  int   offset;
  int   compressSymmetries;
  real8 *etalist, *etanew;
  real8 *alpha, *alphanew, *leveldata = NULL;
  real8 r[3], ep[3];
  real8 epvec[EPVEC_LEN][3] = {
             {-.5, 0.0, 0.0}, {+.5, 0.0, 0.0},
             {0.0, -.5, 0.0}, {0.0, +.5, 0.0},
             {0.0, 0.0, -.5}, {0.0, 0.0, +.5}};
  real8 dsigma[3][3];
  real8 mu,nu,boxl,w;
  real8 (*tsigmadata)[EPVEC_LEN][3][3];
  real8 *taylordata;

  int i0,i1,nlist,*idxlist;

  FILE *fp;
#ifdef PARALLEL
  MPI_Status status;
#endif

  param = home->param;
  mu = param->shearModulus;
  nu = param->pois;
  nordertab = param->fmMPOrder;
  uorder = param->fmExpansionOrder;
  boxl = param->maxSideX - param->minSideX;
  nlevels = numLevels;
  maxlevels = nlevels;

  norder = MAX(nordertab, 13);

#ifdef ANISOTROPIC
  if (param->fmEnableAnisotropy) {
      anisoEnabled = 1;
  }
#endif

  /*MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_np);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_pid);*/

  nalpha = (uorder+3)*(uorder+2)*(uorder+1)/6;
  neta = 9*(norder+3)*(norder+2)*(norder+1)/6;;

/*
 * For CC and FCC type materials, coefficients that are redundant due
 * to cubic crystal symmetry will be compressed out of the table, for all
 * other crystal types, the full table will be created.
 */
  compressSymmetries = (param->materialType == MAT_TYPE_BCC) ||
                       (param->materialType == MAT_TYPE_FCC);

  if (compressSymmetries) {
    /* Make compressed table with all cubic symmetris factored out */
    countit(nordertab,&nlist,&idxlist);
  } else {
    /* Make full table, with all symmetries... */
    countit_full(nordertab,&nlist,&idxlist);
  }

  i0 = nlist/mpi_np * mpi_pid;
  if(mpi_pid < nlist%mpi_np) {
    i0 += mpi_pid;
    i1 = i0 + nlist/mpi_np + 1;
  } else {
    i0 += nlist%mpi_np;
    i1 = i0 + nlist/mpi_np;
  }

  taylordata = (real8 *) malloc((i1-i0) * nalpha * 6 * sizeof(real8));
  tsigmadata = (real8 (*)[6][3][3]) malloc(sizeof(real8)*9*nep*(i1-i0));

  etalist = (real8 *) malloc(nlevels*neta*sizeof(real8));
  etanew = (real8 *) malloc(neta*sizeof(real8));
  alpha = (real8 *) malloc(nalpha*9*sizeof(real8));
  alphanew = (real8 *) malloc(nalpha*9*sizeof(real8));

#ifdef DO_ALL_LEVELS
  leveldata = (real8 *) malloc((nlevels-2)*(i1-i0)*nalpha*9*sizeof(real8));
#endif

  if ((taylordata==NULL) ||
      (tsigmadata==NULL) ||
      (etalist==NULL)    ||
      (etanew==NULL)     ||
#if DO_ALL_LEVELS
      (leveldata==NULL)  ||
#endif
      (alpha==NULL)      ||
      (alphanew==NULL))  {
      Fatal("CreateCorrectiontable: Memory allocation error");
  }


  if(leveldata != NULL)
    for(i = 0; i<(maxlevels-2)*(i1-i0)*nalpha; i++)
      for(m = 0; m<3; m++) 
          for(n = 0; n<3; n++) leveldata[i*9+m*3+n] = 0.0;

  for(nlevels = maxlevels; nlevels>=3; nlevels--) {
    if(nlevels < maxlevels && leveldata == NULL) continue;

    for(i = 0; i<(i1-i0)*nalpha*6; i++) taylordata[i] = 0.0;


    for(i = i0; i<i1; i++) {
      printf("Task %3d  Coefficient %5d of %5d\n",mpi_pid,i+1,i1);

      for(j = 0; j<nlevels; j++)
        for(m = 0; m<neta; m++) etalist[j*neta+m] = 0.0;

      etalist[0*neta+idxlist[i]] = 1.0;

      for(k = 0; k<nalpha; k++)
        for(m = 0; m<3; m++)
          for(n = 0; n<3; n++)
            alpha[k*9+m*3+n] = 0.0;

      if(idxlist[i] > 8) {

        /* Upward pass */
        w = boxl;
        for(j = 1; j<nlevels; j++) {
          /*printf("  Level %2d of %2d\n",j+1,nlevels);*/
          for(iz=-1; iz<=1; iz++) if(pbc[2] || iz==0)
            for(iy=-1; iy<=1; iy++) if(pbc[1] || iy==0)
              for(ix=-1; ix<=1; ix++) if(pbc[0] || ix==0) {
                r[0] = ix*w;
                r[1] = iy*w;
                r[2] = iz*w;
                FMShift(norder,r,&etalist[(j-1)*neta],etanew);
                for(m = 0; m<neta; m++)
                  etalist[j*neta+m] += etanew[m];
              }
          w = 3*w;
        }

        for(k = 0; k<nalpha; k++)
          for(m = 0; m<3; m++)
            for(n = 0; n<3; n++)
                alpha[2*9+m*3+n] = 0.0;
        /*printf("  * Downward pass\n");*/
        w = w/3;
        for(j = nlevels-3; j>=0; j--) {
          /*printf("    At level %2d\n", j+1);*/
          w = w/3;
          for(iz=-4; iz<=4; iz++) if(pbc[2] || iz==0)
            for(iy=-4; iy<=4; iy++) if(pbc[1] || iy==0)
              for(ix=-4; ix<=4; ix++) if(pbc[0] || ix==0)
                if(ix*ix+iy*iy+iz*iz > 3) { /* Not close neighbor */
                  ep[0] = -w*ix;
                  ep[1] = -w*iy;
                  ep[2] = -w*iz;

                  if(nlevels == maxlevels) {

                    MkTaylor(home, mu,nu,norder,uorder,MAX(norder,uorder),
                             ep,&etalist[j*neta],alphanew);
                    if(leveldata == NULL) {
                      for(k = 0; k<nalpha; k++)
                        for(m = 0; m<3; m++)
                          for(n = 0; n<3; n++)
                            alpha[k*9+m*3+n] += alphanew[k*9+m*3+n];
                    } else {
                      for(k = 0; k<nalpha; k++) 
                        for(m = 0; m<3; m++)
                          for(n = 0; n<3; n++)
                            offset = (((i-i0)*(maxlevels-2)+j)*nalpha+k)*9;
                            leveldata[offset+m*3+n] += alphanew[k*9+m*3+n];
                    }
                  }
                }

          if(leveldata != NULL) {
            for(k = 0; k<nalpha; k++)
              for(m = 0; m<3; m++)
                for(n = 0; n<3; n++)
                  offset=(((i-i0)*(maxlevels-2)+j)*nalpha+k)*9;
                  alpha[k*9+m*3+n] +=
                    leveldata[offset+m*3+n];
          }

        }

        for(k = 0; k<nep; k++) {
          ep[0] = boxl*epvec[k][0];
          ep[1] = boxl*epvec[k][1];
          ep[2] = boxl*epvec[k][2];
          EvalTaylor(uorder,ep,alpha,dsigma);
          for(m = 0; m<9; m++)
            tsigmadata[i-i0][k][m/3][m%3] = dsigma[m/3][m%3];
        }

        if(fabs(w-boxl) > 1e-12) printf("Boxlength mismatch! w=%f boxl=%f\n",
                                        w,boxl);
        for(iz=-1; iz<2; iz++) if(pbc[2] || iz==0)
          for(iy=-1; iy<2; iy++) if(pbc[1] || iy==0)
            for(ix=-1; ix<2; ix++) if(pbc[0] || ix==0)
              for(k = 0; k<nep; k++) {
                ep[0] = boxl*epvec[k][0] - boxl*ix;
                ep[1] = boxl*epvec[k][1] - boxl*iy;
                ep[2] = boxl*epvec[k][2] - boxl*iz;

		MkTaylor(home,mu,nu,norder,0,norder,
			 ep,&etalist[0*neta],dsigma[0]);

                for(m = 0; m<9; m++)
                  tsigmadata[i-i0][k][m/3][m%3] += dsigma[m/3][m%3];

              }  /* end for(k = 0; k<nep; k++) */

      }  /* end if(idxlist[i] > 8) */

      for(k = 0; k<nalpha; k++) {
        j = 0;
        for(m = 0; m<9; m++)
          if(m%3 >= m/3)
            taylordata[((i-i0)*nalpha+k)*6 + j++] = alpha[k*9+(m/3)*3+(m%3)];
      }
    }

    if(mpi_pid == 0) {

      fp = fopen(param->fmCorrectionTbl,"w");

      if (fp == (FILE *)NULL) {
          Fatal("Error %d opening output file %s", errno,
                param->fmCorrectionTbl);
      }

      fprintf(fp, "%s%d\n", CTABLE_VERSION_STR, CTABLE_VERSION);
      fprintf(fp, "%d %d %d  # PBC flags in X, Y and Z respectively\n",
              pbc[0], pbc[1], pbc[2]);
      fprintf(fp, "%-25.15lf   # MU\n", mu);
      fprintf(fp, "%-25.15lf   # NU\n", nu);
      fprintf(fp, "%-25.15lf   # Simulation cube size\n", boxl);
      fprintf(fp, "%d   # Multipole expansion order\n", nordertab);
      fprintf(fp, "%d   # Taylor expansion order\n", uorder);
      fprintf(fp, "%d   # # of levels of PBC images\n", nlevels);
      fprintf(fp, "%d   # Compression flag (1 if cubic symmetries are "
                  "compresed out\n", compressSymmetries);
      fprintf(fp, "%d   # Anisotropy flag (1 == enabled)\n", anisoEnabled);


#ifdef ANISOTROPIC
/*
 *    When support for anisotropy is compiled in and enabled, we have to
 *    dump to the file the elastic constants matrix used when creating the
 *    correction table
 */
      if (anisoEnabled) {
          for (i = 0; i < 6; i++) {
              for (j = 0; j < 6; j++) {
                  fprintf(fp, "%-25.15lf   # C[%d][%d]\n",
                          param->elasticConstantMatrix[i][j], i, j);
              }
          }
      }
#endif

      for(i = 0; i<mpi_np; i++) {
        int ixx;
        ixx = nlist/mpi_np + (i < nlist%mpi_np);
        if(i > 0) {
#ifdef PARALLEL
          MPI_Send(&i,1,MPI_INT,i,999,MPI_COMM_WORLD);
          MPI_Recv(taylordata,ixx*nalpha*6,MPI_DOUBLE,i,1000,MPI_COMM_WORLD,
                   &status);
#endif
        }
        for(k = 0; k<6*ixx*nalpha; k++) {
          fprintf(fp,"%+.15e%s",taylordata[k],(k%6 < 5) ? " ":"\n");
          if(k%(6*nalpha) == 6*nalpha-1) fprintf(fp,"\n");
        }
      }
    } else {
#ifdef PARALLEL
      MPI_Recv(&i,1,MPI_INT,0,999,MPI_COMM_WORLD,&status);
      MPI_Send(taylordata,(i1-i0)*nalpha*6,MPI_DOUBLE,0,1000,MPI_COMM_WORLD);
#endif
    }
    if(mpi_pid == 0) fclose(fp);

  }

  free(taylordata);
  free(tsigmadata);

  free(etalist);
  free(etanew);
  free(alpha);
  free(alphanew);

  if(leveldata != NULL) free(leveldata);

  return;
}


void Taylor_FreeCorrectionTable(void)
{
        if (correctionTbl != (real8 *)NULL) {
            free(correctionTbl);
            correctionTbl = (real8 *)NULL;
        }

        return;
}


void Taylor_CorrectionTableInit(Home_t *home)
{
        int     n, k, j, i;
        int     version;
        int     pbc[3], ctab_pbc[3];
        int     anisoEnabled = 0;
        int     expectedCompression;
        int     tableCompressionFlag = 1; /* default assumes compression */
        real8   eps = 1.0e-06;
        real8   NU;
        char    *fileName;
        char    line[256];
        FILE    *fp;
        Param_t *param;
#ifdef ANISOTROPIC
        real8   relativeTol = 1.0e-5;
        real8   C[6][6];
#endif

        param = home->param;
        fileName = param->fmCorrectionTbl;
        NU = param->pois;

/*
 *      Just a sanity check to verify we have period boundaries
 *      in at least one dimension.
 */
        if ((param->xBoundType != Periodic) &&
            (param->yBoundType != Periodic) &&
            (param->zBoundType != Periodic)) {
            return;
        }

        pbc[0] = (param->xBoundType == Periodic);
        pbc[1] = (param->yBoundType == Periodic);
        pbc[2] = (param->zBoundType == Periodic);

/*
 *      Default correction table to PBC in all dimensions for
 *      older correction table files which assume PBC but do
 *      not specify it within the correction table file itself
 */
        ctab_pbc[0] = 1;
        ctab_pbc[1] = 1;
        ctab_pbc[2] = 1;
/*
 *      Only the domain owning the single FM cell at the coarsest
 *      layer has anything to do here.
 */
        if (home->fmLayer[0].ownedCnt == 0) return;

        if (fileName[0] == 0) {
            Fatal("Error: The file containing the PBC image correction\n"
                  "    table was not specified.  Please set the \n"
                  "    <fmCorrectionTbl> value in the control file\n"
                  "    to the name of the file containing the table.\n");
        }

        if ((fp = fopen(fileName, "r")) == (FILE *)NULL) {
            Fatal("Error %d opening PBC correction table file %s",
                  errno, fileName);
        }

/*
 *      Based on the material type set a flag indicating if we expect a
 *      full correction table or one with cubic symmetries compressed out.
 */
        expectedCompression = (param->materialType == MAT_TYPE_BCC) ||
                              (param->materialType == MAT_TYPE_FCC);

/*
 *      Get the version number of the correction table (if any)
 */
        Getline(line, sizeof(line)-1, fp);

        if (strstr(line, CTABLE_VERSION_STR) == (char *)NULL) {
            version = 0;
        } else {
            version = atoi(&line[strlen(CTABLE_VERSION_STR)]);
        }

/*
 *      Each correction table contains the specific constants
 *      such as MU, NU, multipole and taylor expansion orders, etc
 *      for which the table was built.  Read these constants 
 *      from the file and compare to the same values for this
 *      run.  If there is a mis-match, abort now.
 */
        switch (version) {
            case 0:
/*
 *              MU, NU and simulation box length are already in <line>.
 *              Other values must be read from the file.
 */
                sscanf(line,"%lf%lf%lf", &ctab_MU, &ctab_NU, &ctab_boxl);
                fscanf(fp,"%d%d%d", &ctab_norder, &ctab_uorder, &ctab_levels);
                break;

            case 1:
/*
 *              MU, NU, multipole and taylor expansion orders, etc. are
 *              written 1 per line possibly followed by comments
 */
                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%d%d%d", &ctab_pbc[0],
                       &ctab_pbc[1], &ctab_pbc[2]);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%lf", &ctab_MU);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%lf", &ctab_NU);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%lf", &ctab_boxl);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%d", &ctab_norder);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%d", &ctab_uorder);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%d", &ctab_levels);
                break;

            case 2:
/*
 *              MU, NU, multipole and taylor expansion orders, etc. are
 *              written 1 per line possibly followed by comments
 */
                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%d%d%d", &ctab_pbc[0],
                       &ctab_pbc[1], &ctab_pbc[2]);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%lf", &ctab_MU);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%lf", &ctab_NU);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%lf", &ctab_boxl);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%d", &ctab_norder);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%d", &ctab_uorder);

                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%d", &ctab_levels);

/*
 *              Next, we need to know if cubic crystal symmetries have
 *              been compressed out of the table.  Also check that the
 *              table compression state matches the state expected by
 *              the simulation.  (i.e. compressed for BCC or FCC, not
 *              compressed for any other material type)
 */
                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%d", &tableCompressionFlag);

                if (expectedCompression != tableCompressionFlag) {
                    Fatal("Expected correction table with cubic symmetries\n"
                          "%s but found %s table",
                          (expectedCompression == 1 ? "compressed" :
                                                      "uncompressed"),
                          (tableCompressionFlag == 1 ? "compressed" :
                                                       "uncompressed"));
                }

/*
 *              Next is anisotropy flag.  If it is set, we expect to find
 *              the 6x6 elastic constants matrix next.  If it is set and
 *              anisotropy has not been compiled in, we have a problem.
 */
                Getline(line, sizeof(line)-1, fp);
                sscanf(line,"%d", &anisoEnabled);

#ifdef ANISOTROPIC
                if ((param->fmEnableAnisotropy == 0) && (anisoEnabled)) {
                    Fatal("CorrectionTableInit: Correction table was "
                          "built with\nanisotropy enabled, but current "
                          "execution requests isotropic FMM.");
                }

                if ((param->fmEnableAnisotropy) && (!anisoEnabled)) {
                    Fatal("CorrectionTableInit: Correction table was not "
                          "built with\nanisotropy enabled, but current\n"
                          "execution requests anisotropic FMM.");
                }
                if (anisoEnabled) {
                    for (i = 0; i < 6; i++) {
                        for (j = 0; j < 6; j++) {
                            Getline(line, sizeof(line)-1, fp);
                            sscanf(line,"%lf", &C[i][j]);
                        }
                    }
                }
#else  /* ANISOTROPIC not defined */
                if (anisoEnabled) {
                    Fatal("CorrectionTableInit: Correction table built with "
                          "anisotropy enabled,\nbut current executable is "
                          "not built with anisotropy");
                }
#endif  /* ifdef ANISOTROPIC */

                break;

            default:
                Fatal("CorrectionTableInit: Unknown table version (%s)", line);
                break;
        }

#ifdef ANISOTROPIC
/*
 *      If anisotropic FMM has been requested, and anisotropy has been
 *      enabled at run time, and we expect a full correction table
 *      in which no symmetries have been compressed out, the correction
 *      table read in MUST be a full table. 
 *
 */
        if ((param->fmEnableAnisotropy) && (expectedCompression == 0) &&
            (tableCompressionFlag == 1)) {
            Fatal("Expected full correction table but found a table\n"
                  "with cubic symmetries compressed out");
        }

/*
 *      When the table was created with anisotropic constants we have to
 *      make sure the elastic constants used in creating the correction
 *      table match the elastic constants for the current simulation.
 */
        if ((!anisoEnabled) && (param->fmEnableAnisotropy)) {
            Fatal("CorrectionTableInit: Simulation is requesting anisotropic\n"
                  "FMM, but the correction table has not been built "
                  "with anisotropy");
        }
        if (anisoEnabled) {
            if (!Matrix66_Near(param->elasticConstantMatrix, C, relativeTol)) {
                Fatal("%s: C matrix from table does not match C matrix\n"
                      "for current simulation", __func__);
            }
        }
#endif

        for (i = 0; i < 3; i++) {
            if (pbc[i] != ctab_pbc[i]) {
                Fatal("CorrectionTableInit: Table created %s pbc in %s but "
                      "simulation run %s pbc in %s dimension",
                      (ctab_pbc[i] == 0 ? "without" : "with"),
                      (i == 0 ? "X" : (i == 1 ? "Y" : "Z")),
                      (pbc[i] == 0 ? "without" : "with"),
                      (i == 0 ? "X" : (i == 1 ? "Y" : "Z")));
            }
        }

        if (param->fmMPOrder != ctab_norder) {
            Fatal("Error: Current multipole expansion order "
                  "= %d, but\n    %s is built for multipole "
                  "expansion order = %d",
                  param->fmMPOrder, fileName, ctab_norder);
        }

        if (param->fmExpansionOrder != ctab_uorder) {
            Fatal("Error: Current taylor expansion order "
                  "= %d, but\n    %s is built for taylor "
                  "expansion order = %d",
                  param->fmExpansionOrder, fileName, ctab_uorder);
        }
 
        if ((NU < 0.0) || ((NU > 0.0) && fabs((NU-ctab_NU)/NU) > eps)) {
            Fatal("Mismatch between configured value "
                  "for pois\n    and value used to build FMM PBC correction "
                  "table:\n    Configured pois = %e    Table pois = %e"
                  "    Use ctablegen to create a new correction table",
                  NU, ctab_NU);
        }

/*
 *      The correction table appears to be consistent with the current
 *      simulation.  We do, however, have to preserve for later use the
 *      flag indicating if the table is full or if there were cubic
 *      symmetries that were compressed out.
 */
        param->fmFullCorrectionTbl = (tableCompressionFlag == 0);

/*
 *      Now we just read in the proper number of coefficients based
 *      on the multipole order, taylor order and compression flag.
 */
        n = (ctab_uorder+3)*(ctab_uorder+2)*(ctab_uorder+1)/6;

        if (tableCompressionFlag == 1) {
            countit(ctab_norder, &k, NULL);
        } else {
            countit_full(ctab_norder, &k, NULL);
        }

        correctionTbl = (real8 *)malloc(6 * k * n * sizeof(real8));

        for (i = 0; i < n*k; i++) {
            for (j = 0; j < 6; j++) {
                int numRead;
                numRead = fscanf(fp, "%lf", &correctionTbl[i*6+j]);
                if (numRead != 1) {
                    Fatal("EOF reading FMM correction table: i = %d, j = %d",
                          i, j);
                }
            }
        }

        fclose(fp);

        return;
}


static void paxby(int uorder, int p[], int q[], real8 weight,
                  real8 *corrtab, real8 *alpha) {
  int i,k,n,m,nx,ny,nz,ix,iy,iz,u,v,u1,v1,c[3];
  int offsetA1, offsetA2;
  int nts[3][3] = { {0,1,2} ,
                    {1,3,4} ,
                    {2,4,5} };
  k = 0;
  for(i = 0; i<=uorder; i++) {
    n = (i+2)*(i+1)*i/6;
    for(nz = 0; nz<=i; nz++)
      for(ny = 0; ny<=i-nz; ny++) {
        nx = i-nz-ny;
        c[0] = nx; c[1] = ny; c[2] = nz;
        ix = c[p[0]]; iy = c[p[1]]; iz = c[p[2]];

        m = iz*(i+1) + iy - (iz-1)*iz/2;

        offsetA1 = (n+m)*9;
        for(u = 0; u<3; u++) {
          u1 = q[u];
          offsetA2 = offsetA1 + u1*3;
          for(v = 0; v<3; v++) {
            v1 = q[v];
            alpha[offsetA2+v1] += weight*corrtab[(k*6)+nts[u][v]];
          }
        }
        k = k + 1;
      }
  }
}


void Taylor_DoTableCorrection(Home_t *home)
{
        int       i, j, k, m, n, u, v, ix, iy, iz, u1, v1, ip;
        int       nx, ny, nz, idx, norder, uorder, nalpha;
        int       offset1, offset2;
        int       isFullTable;
        int       *p, c[3], q[3];
        int       flags[(NMAX+2)*(NMAX+1)/2];
        int       ptab[][3] = {{0,1,2}, {1,2,0}, {2,0,1},
                               {2,1,0}, {1,0,2}, {0,2,1}};
        real8     s, MU;
        real8     *eta, *alpha;
        Param_t   *param;
        FMCell_t  *cell;
        FMLayer_t *layer;

        param = home->param;
        MU = param->shearModulus;

/*
 *      For cubic systems, some of the symmetries for cubic crystals 
 *      may have been compressed out of the correction table.  Set
 *      a flag indicating if the table is full or compressed.
 */
        isFullTable = param->fmFullCorrectionTbl;

/*
 *      Only the domain owning the single FM cell at the highest(coarsest)
 *      layer has anything to do here.
 */
        layer = &home->fmLayer[0];
        if (layer->ownedCnt == 0) return;

        norder = param->fmMPOrder;

        cell = LookupFMCell(layer->cellTable, 0);
        eta   = cell->mpCoeff;
 
        uorder = param->fmExpansionOrder;
        alpha = cell->expansionCoeff;
        nalpha = (uorder+3)*(uorder+2)*(uorder+1)/6;

        for (i = 0; i < nalpha; i++) {
            offset1 = i*9;
            for (u = 0; u < 3; u++) {
                offset2 = offset1 + (u*3);
                for (v = 0; v < 3; v++) {
                    alpha[offset2+v] = 0.0;
                }
            }
        }

/*
 *      Rescale eta to accomodate any difference between the dimensions
 *      of the current problem box length and the box length for which
 *      this correction table was created.
 */
        s = ctab_boxl/param->Lx;

        for (i = 0; i <= norder; i++) {
            idx = 9*(i+2)*(i+1)*i/6;
            n = 9*(i+2)*(i+1)/2;
            s = s * (ctab_boxl/param->Lx);
            for (j = 0; j < n; j++) {
                eta[idx+j] = eta[idx+j] * s;
            }
        }

        k = -1;

        for (i = 0; i<=norder; i++) {
            idx = 9*(i+2)*(i+1)*i/6;
            n = (i+2)*(i+1)/2;

            for (j = 0; j<n; j++) {
                flags[j] = 0;
            }

            for (nz = 0; nz<=i; nz++) {
                for (ny = 0; ny<=i-nz; ny++) {

                    nx = i-nz-ny;
		    
		    if (isFullTable == 0) {
                      /* cubic symmetries have been compressed out of table */
		      if (nx>=ny && ny>=nz) {
                          k = k+1;
                      }
		    } else {
		      k = k+1;
		    }

                    c[0] = nx;
                    c[1] = ny;
                    c[2] = nz;

/*                  ip=0 is the identity permutation...
 *
 *                  Need to check if the table is full or if cubic symmetries
 *                  were factored out...
 */
                    for (ip = 0; ip<( (isFullTable==0) ? 6:1); ip++) {
                        if (ip >= 3 && (nx!=ny && nx!=nz && ny!=nz)) s = -1.0;
                        else s = 1.0;
                        p = ptab[ip];
                        for (j = 0; j<3; j++) q[p[j]] = j;
                        ix = c[p[0]]; iy = c[p[1]]; iz = c[p[2]];
                        m = iz*(i+1) + iy - (iz-1)*iz/2;
        
                        if (flags[m] == 0) {
                            flags[m] = 1;
                            for (u = 0; u<3; u++) {
                                for (v = 0; v<3; v++) {
                                    u1 = q[u]; v1 = q[v];
                                    paxby(uorder,p,q,s*eta[idx+n*(3*u1+v1)+m],
                                        &correctionTbl[nalpha*(9*k+3*u+v)*6],alpha);
                                }
                            }
                        }
                    }
                }
            }
        }

/*
 *      Rescale alpha to match mu and boxlength with which the
 *      correction table was created.
 */
        s = 1.0;

        for (i = 0; i <= uorder; i++) {
            idx = (i+2)*(i+1)*i/6;
            n = (i+2)*(i+1)/2;
            for (j = 0; j < n; j++) {
                offset1 = (idx+j)*9;
                for (k = 0; k < 3; k++) {
                    offset2 = offset1 + k*3;
                    for (m = 0; m < 3; m++) {
                    alpha[offset2+m] *= s;
                    }
                }
            }
            s = s * (ctab_boxl/param->Lx);
        }

        return;
}

#endif
