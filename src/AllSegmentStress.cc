
//#ifdef PRINTSTRESS
#if defined(PRINTSTRESS) || defined(SPECTRAL)
/***************************************************************************
 *
 *      Functionf:     AllSegmentStress
 *      Description:   Calculate stress field at a point from all segments,
 *                     local and remote ones. PBC included.
 *
 *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"
#include "Util.h"
#include "Matrix.h"
#include "FM.h"

static int isRightDom(Home_t *home,real8 x, real8 y, real8 z)
{
  int DomID = -1;
  real8   xMin,yMin,zMin,xMax,yMax,zMax;

  xMin = home->domXmin;yMin = home->domYmin;zMin = home->domZmin;
  xMax = home->domXmax;yMax = home->domYmax;zMax = home->domZmax;

  if (((x >= xMin) && (x < xMax)) &&
      ((y >= yMin) && (y < yMax)) &&
      ((z >= zMin) && (z < zMax)))
    DomID = home->myDomain;

  return DomID;
}



static void xunitvector(real8 ax, real8 ay, real8 az,
                        real8 bx, real8 by, real8 bz,
                        real8 *cx, real8 *cy, real8 *cz)
{
     real8 L, dx, dy, dz;
     xvector(ax,ay,az, bx,by,bz, cx,cy,cz);
     L=(*cx)*(*cx)+(*cy)*(*cy)+(*cz)*(*cz);
     if(L<FFACTOR_LMIN2)
     {
         dx=bx+1; dy=by; dz=bz;
         xvector(ax,ay,az, dx,dy,dz, cx,cy,cz);
         L=(*cx)*(*cx)+(*cy)*(*cy)+(*cz)*(*cz);
         if(L<FFACTOR_LMIN2)
         {
             dx=bx; dy=by+1; dz=bz;
             xvector(ax,ay,az, dx,dy,dz, cx,cy,cz);
             L=(*cx)*(*cx)+(*cy)*(*cy)+(*cz)*(*cz);
         }
     }
     L=sqrt(L);
     *cx/=L; *cy/=L; *cz/=L;
}


static void SegmentStressHor(real8 MU, real8 NU,
                      real8 bx, real8 by, real8 bz,
                      real8 z1, real8 z2,
                      real8 x, real8 z,
                      real8 a,
                      real8 Sigma[3][3] )
{
/*
 *   segment lies horizontally (0,0,z1) to (0,0,z2)
 *   field point (x,0,z)  (i.e. y=0)
 *   Burgers vector (bx,by,bz)
 *   shear modulus MU, Poisson's ratio NU,
 *   regularization radius a (dislocation core width),
 *   return stress: Sigma[3][3]
 */
    real8 s[6];
    real8 L, Lp, L2, Lp2, Ra, Ra2, Ra3, Rap, Rap2, Rap3;
    real8 x2, a2, sunit, invral,invrapl, rhoa2;
    int form;
/*
 *    printf("%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
 *            bx,by,bz,z1,z2,x,z);
 */
    L=z1-z; Lp=z2-z;
    L2=L*L; Lp2=Lp*Lp; x2=x*x; a2=a*a;

    Ra2=x2+L2+a2;   Ra=sqrt(Ra2);   Ra3=Ra*Ra2;
    Rap2=x2+Lp2+a2; Rap=sqrt(Rap2); Rap3=Rap*Rap2;


    sunit=MU*(1.0/4.0/M_PI/(1-NU));

    if((L>0)&&(Lp>0))
        form=1;
    else if((L<0)&&(Lp<0))
        form=2;
    else
        form=3;

    /* for debug purposes (form 1,2,3 should yield same result) */
    /* form=3; */

    switch(form) {
    case(1):
        invral=1/Ra/(Ra+L);
        invrapl=1/Rap/(Rap+Lp);
        s[0]= by*x*( invrapl*(1-(x2+a2)/Rap2-(x2+a2)*invrapl)
                   - invral *(1-(x2+a2)/Ra2 -(x2+a2)*invral ) );
        s[1]=-bx*x*( invrapl - invral );
        s[2]= by*( (-NU/Rap+(x2+(1-NU)*a2/2.0)/Rap3)
                  -(-NU/Ra +(x2+(1-NU)*a2/2.0)/Ra3 ) );
        s[3]=-by*x*( invrapl*(1+a2/Rap2+a2*invrapl)
                    -invral *(1+a2/Ra2 +a2*invral ) );
        s[4]= (  bx*(NU/Rap-(1-NU)*(a2/2.0)/Rap3)
                -bz*x*(1-NU)*invrapl*(1+(a2/2.0)/Rap2+(a2/2.0)*invrapl) )
             -(  bx*(NU/Ra -(1-NU)*(a2/2.0)/Ra3 )
                -bz*x*(1-NU)*invral *(1+(a2/2.0)/Ra2 +(a2/2.0)*invral ) );
        s[5]= by*x*( (-2*NU*invrapl*(1+(a2/2.0)/Rap2+(a2/2.0)*invrapl)-Lp/Rap3)
                    -(-2*NU*invral *(1+(a2/2.0)/Ra2 +(a2/2.0)*invral )-L /Ra3 ) );
        break;
    case(2):
        invral=1/Ra/(Ra-L);
        invrapl=1/Rap/(Rap-Lp);
        s[0]=-by*x*( invrapl*(1-(x2+a2)/Rap2-(x2+a2)*invrapl)
                   - invral *(1-(x2+a2)/Ra2 -(x2+a2)*invral ) );
        s[1]= bx*x*( invrapl - invral );
        s[2]= by*( (-NU/Rap+(x2+(1-NU)*a2/2.0)/Rap3)
                  -(-NU/Ra +(x2+(1-NU)*a2/2.0)/Ra3 ) );
        s[3]= by*x*( invrapl*(1+a2/Rap2+a2*invrapl)
                    -invral *(1+a2/Ra2 +a2*invral ) );
        s[4]= (  bx*(NU/Rap-(1-NU)*(a2/2.0)/Rap3)
                +bz*x*(1-NU)*invrapl*(1+(a2/2.0)/Rap2+(a2/2.0)*invrapl) )
             -(  bx*(NU/Ra -(1-NU)*(a2/2.0)/Ra3 )
                +bz*x*(1-NU)*invral *(1+(a2/2.0)/Ra2 +(a2/2.0)*invral ) );
        s[5]= by*x*( ( 2*NU*invrapl*(1+(a2/2.0)/Rap2+(a2/2.0)*invrapl)-Lp/Rap3)
                    -( 2*NU*invral *(1+(a2/2.0)/Ra2 +(a2/2.0)*invral )-L /Ra3 ) );
        break;
    case(3):
        rhoa2=x2+a2;
        s[0]= by*x/rhoa2*( Lp/Rap*(1+rhoa2/Rap2)
                          -L /Ra *(1+rhoa2/Ra2 ) );
        s[1]= bx*x/rhoa2*( Lp/Rap - L/Ra );
        s[2]= by*( (-NU/Rap+x2/Rap3+(1-NU)*(a2/2.0)/Rap3)
                  -(-NU/Ra +x2/Ra3 +(1-NU)*(a2/2.0)/Ra3 ) );
        s[3]= by*x/rhoa2*( Lp/Rap*(1+(a2*2.0)/rhoa2+a2/Rap2)
                          -L /Ra *(1+(a2*2.0)/rhoa2+a2/Ra2 ) );
        s[4]= (  bx*(NU/Rap-(1-NU)*(a2/2.0)/Rap3)
                +bz*x/rhoa2*(1-NU)*Lp/Rap*(1+a2/rhoa2+(a2/2.0)/Rap2) )
             -(  bx*(NU/Ra -(1-NU)*(a2/2.0)/Ra3 )
                +bz*x/rhoa2*(1-NU)*L /Ra *(1+a2/rhoa2+(a2/2.0)/Ra2 ) );
        s[5]= by*x*( ( 2*NU*Lp/rhoa2/Rap*(1+a2/rhoa2+(a2/2.0)/Rap2)-Lp/Rap3)
                    -( 2*NU*L /rhoa2/Ra *(1+a2/rhoa2+(a2/2.0)/Ra2 )-L /Ra3 ) );
        break;
    default:
        Fatal("SegmentStressHor: unknown form!");
    }
    Sigma[0][0]=sunit*s[0];
    Sigma[0][1]=sunit*s[1];
    Sigma[0][2]=sunit*s[2];
    Sigma[1][1]=sunit*s[3];
    Sigma[1][2]=sunit*s[4];
    Sigma[2][2]=sunit*s[5];

    Sigma[1][0]=Sigma[0][1];
    Sigma[2][0]=Sigma[0][2];
    Sigma[2][1]=Sigma[1][2];
}


static void SegmentStress(real8 MU, real8 NU,
                   real8 bX, real8 bY, real8 bZ,
                   real8 xA, real8 yA, real8 zA,
                   real8 xB, real8 yB, real8 zB,
                   real8 x0, real8 y0, real8 z0,
                   real8 a,
                   real8 Sigma[3][3] )
{
/*
 *   Segment goes from (xA,yA,zA)->(xB,yB,zB)
 *   Burgers vector (bX,bY,bZ)
 *   field point at (x0,y0,z0)
 *   shear modulus MU, Poisson's ratio NU,
 *   regularization radius a (dislocation core width),
 *   return stress: Sigma[3][3]
 */
    real8 tmp[3][3], M[3][3], Mt[3][3];
    real8 xAB, yAB, zAB, xA0, yA0, zA0, rAB;
    real8 e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z;
    real8 z1, z2, xr0, zr0, brX, brY, brZ;

    xAB=xB-xA; yAB=yB-yA; zAB=zB-zA;
    xA0=x0-xA; yA0=y0-yA; zA0=z0-zA;

    rAB=sqrt(xAB*xAB+yAB*yAB+zAB*zAB);

    if(rAB<FFACTOR_LMIN)
        Fatal("SegmentStress: rAB < LMIN");

    e1x=xAB/rAB; e1y=yAB/rAB; e1z=zAB/rAB;

    xunitvector(e1x,e1y,e1z,xA0,yA0,zA0,&e3x,&e3y,&e3z);
    xvector(e3x,e3y,e3z,e1x,e1y,e1z,&e2x,&e2y,&e2z);

    /* e1,e2,e3 should all be normalized */

    M[0][0]=e2x; M[0][1]=e3x; M[0][2]=e1x;
    M[1][0]=e2y; M[1][1]=e3y; M[1][2]=e1y;
    M[2][0]=e2z; M[2][1]=e3z; M[2][2]=e1z;

    Mt[0][0]=e2x; Mt[0][1]=e2y; Mt[0][2]=e2z;
    Mt[1][0]=e3x; Mt[1][1]=e3y; Mt[1][2]=e3z;
    Mt[2][0]=e1x; Mt[2][1]=e1y; Mt[2][2]=e1z;

    /*
    printf("%20.10e %20.10e %20.10e\n"
           "%20.10e %20.10e %20.10e\n"
           "%20.10e %20.10e %20.10e\n",
           M[0][0], M[0][1], M[0][2],
           M[1][0], M[1][1], M[1][2],
           M[2][0], M[2][1], M[2][2] );
    */

    /* rr21=Mt*r21*/
    z1=0;
    z2=Mt[2][0]*xAB+Mt[2][1]*yAB+Mt[2][2]*zAB;
    /* rr=Mt*rt*/
    xr0=Mt[0][0]*xA0+Mt[0][1]*yA0+Mt[0][2]*zA0;
    zr0=Mt[2][0]*xA0+Mt[2][1]*yA0+Mt[2][2]*zA0;
    /* br=Mt*b*/
    brX=Mt[0][0]*bX+Mt[0][1]*bY+Mt[0][2]*bZ;
    brY=Mt[1][0]*bX+Mt[1][1]*bY+Mt[1][2]*bZ;
    brZ=Mt[2][0]*bX+Mt[2][1]*bY+Mt[2][2]*bZ;

    SegmentStressHor(MU,NU,brX,brY,brZ,z1,z2,xr0,zr0,a,Sigma);

    Matrix33_Mul(tmp  ,M  ,Sigma);   // tmp   = M*Sigma
    Matrix33_Mul(Sigma,tmp,Mt   );   // Sigma = M*Sigma*Mt

    return;
}



static void LocalStress(Home_t *home, real8 x, real8 y,real8 z,
                        int cellIndex[3],real8 locStress[3][3])
{
  int cellID, arm;
  int cx, cy, cz;
  real8 minXIndex, minYIndex, minZIndex;
  int  maxXIndex, maxYIndex, maxZIndex;
  real8 a, MU, NU;
  real8 bx, by, bz;
  real8 xc, yc, zc;
  real8 p1x, p1y, p1z;
  real8 p2x, p2y, p2z;
  Node_t  *node1, *node2;

  Cell_t  *cell;
  Param_t *param;

  param = home->param;

  a     = param->rc;
  MU    = param->shearModulus;
  NU    = param->pois;

/*
 *  Loop though all the cells in the block.
 *  min and max are in interval [0, param->ncells+1] when PBC are on.
 *  cells 0 and param->ncells+1 are added for PBC.
 */

  if (param->xBoundType == Periodic) {
    minXIndex = MAX(0, cellIndex[0]-1);
    maxXIndex = MIN(param->nXcells+1, cellIndex[0]+1);
  } else {
    minXIndex = MAX(1, cellIndex[0]-1);
    maxXIndex = MIN(param->nXcells, cellIndex[0]+1);
  }

  if (param->yBoundType == Periodic) {
    minYIndex = MAX(0, cellIndex[1]-1);
    maxYIndex = MIN(param->nYcells+1, cellIndex[1]+1);
  } else {
    minYIndex = MAX(1, cellIndex[1]-1);
    maxYIndex = MIN(param->nYcells, cellIndex[1]+1);
  }

  if (param->zBoundType == Periodic) {
    minZIndex = MAX(0, cellIndex[2]-1);
    maxZIndex = MIN(param->nZcells+1, cellIndex[2]+1);
  } else {
    minZIndex = MAX(1, cellIndex[2]-1);
    maxZIndex = MIN(param->nZcells, cellIndex[2]+1);
  }

/*
 *      Loop though all the cells in the block.
 */

  for (cx = (int)minXIndex; cx <= (int)maxXIndex; cx++) {
    for (cy = (int)minYIndex; cy <= (int)maxYIndex; cy++) {
      for (cz = (int)minZIndex; cz <= (int)maxZIndex; cz++) {

         cellID = EncodeCellIdx(home, cx, cy, cz);
         cell = LookupCell(home, cellID);

	if (cell == (Cell_t *)NULL) continue;


	/* Find center of this cell */
        xc = cell->center[X];
        yc = cell->center[Y];
        zc = cell->center[Z];

	/* put cell center to the nearest image of field point */
	PBCPOSITION(param, x, y, z, &xc, &yc, &zc);

/*
 *                  Loop over all nodes in this cell and over each segment
 *                  attached to the node.  Skip any segment that is not
 *                  owned by node1.
 */

	node1 = cell->nodeQ;

	for (; node1 != (Node_t *)NULL; node1=node1->nextInCell) {
	  for (arm = 0; arm < node1->numNbrs; arm++) {

	    node2 = GetNeighborNode(home, node1, arm);

	    if (node2 == (Node_t *)NULL ) continue;

            if (NodeOwnsSeg(home, node1, node2) == 0) {
               continue;
            }

	    p1x = node1->x;
	    p1y = node1->y;
	    p1z = node1->z;

	    p2x = node2->x;
	    p2y = node2->y;
	    p2z = node2->z;

	    PBCPOSITION(param, xc,  yc,  zc,  &p1x, &p1y, &p1z);
	    PBCPOSITION(param, p1x, p1y, p1z, &p2x, &p2y, &p2z);

	    bx = node1->burgX[arm];
	    by = node1->burgY[arm];
	    bz = node1->burgZ[arm];

#if 0
            real8 sigma[3][3];
            Matrix33_Zero(sigma);
            SegmentStress(MU, NU, bx, by, bz,
                          p1x, p1y, p1z,
                          p2x, p2y, p2z,
                          x, y, z, a, sigma);

            locStress[0][0] += sigma[0][0];
            locStress[1][1] += sigma[1][1];
            locStress[2][2] += sigma[2][2];
            locStress[0][1] += sigma[0][1];
            locStress[1][2] += sigma[1][2];
            locStress[0][2] += sigma[0][2];
#else
            /* Use StressDueToSeg function so that it is also
            valid for anisotropic elasticity - Nicolas */
            double sigma[6];
            StressDueToSeg(home, x, y, z, p1x, p1y, p1z,
                           p2x, p2y, p2z, bx, by, bz,
                           a, MU, NU, sigma);

            locStress[0][0] += sigma[0];
            locStress[1][1] += sigma[1];
            locStress[2][2] += sigma[2];
            locStress[0][1] += sigma[5];
            locStress[1][2] += sigma[3];
            locStress[0][2] += sigma[4];
#endif

 	  }
	}

      } /* end for(cz = 0; ...) */
    } /* end for(cy = 0; ...) */
  } /* end for(cx = 0; ...) */

  locStress[1][0] = locStress[0][1];
  locStress[2][0] = locStress[0][2];
  locStress[2][1] = locStress[1][2];
}

static void RemoteStressWithFMM(Home_t *home,real8 x, real8 y,real8 z,
		       int cellIndex[3],real8 remStress[3][3])
{
  int cellID;
  real8   R[3];
  FMLayer_t *layer;
  FMCell_t  *FMMcell;
  Param_t *param;

  param = home->param;

  cellIndex[0] --;
  cellIndex[1] --;
  cellIndex[2] --;

  layer = &home->fmLayer[param->fmNumLayers-1];

  cellID = EncodeFMCellIndex(layer->lDim, cellIndex[0],cellIndex[1],cellIndex[2]);

  FMMcell  = LookupFMCell(layer->cellTable, cellID);

  R[X] = x - FMMcell->cellCtr[X];
  R[Y] = y - FMMcell->cellCtr[Y];
  R[Z] = z - FMMcell->cellCtr[Z];

  ZImage(param, &R[X], &R[Y], &R[Z]);

  LocalToPoint(home, layer->cellSize, R, FMMcell->expansionCoeff, remStress);
}


#if 0

static void RemoteStressWithTable(Home_t *home,real8 x, real8 y,real8 z,
		       int cellIndex[3],real8 remStress[3][3])
{
  int i,j;
  int     xSkip1, xSkip2, xSkip3;
  int     ySkip1, ySkip2, ySkip3;
  int     zSkip1, zSkip2, zSkip3;
  int     cx, cy, cz;
  int     includePrimary,cellID;
  real8   cellXsize, cellYsize, cellZsize;
  real8   dx, dy, dz;
  real8   delSig[3][3];
  real8   burgX, burgY, burgZ;
  real8   xc, yc, zc;

  Param_t *param;

  param = home->param;


/*
 *      First time intop this function, allocate and initialize some
 *      static arrays
 */
  if (cellCterX == (real8 *)NULL) {

    real8   Lx, Ly, Lz;
    real8   cellXsize, cellYsize, cellZsize, xStart, yStart, zStart;

    Lx = param->Lx;
    Ly = param->Ly;
    Lz = param->Lz;

    cellXsize = Lx / param->nXcells;
    cellYsize = Ly / param->nYcells;
    cellZsize = Lz / param->nZcells;

    xStart = param->minSideX + cellXsize*0.5;
    yStart = param->minSideY + cellYsize*0.5;
    zStart = param->minSideZ + cellZsize*0.5;

    cellCterX = (real8 *) malloc(param->nXcells * sizeof(real8));
    cellCterY = (real8 *) malloc(param->nYcells * sizeof(real8));
    cellCterZ = (real8 *) malloc(param->nZcells * sizeof(real8));

    for (i = 0; i < param->nXcells; i++)
      cellCterX[i] = xStart + i*cellXsize;

    for (i = 0; i < param->nYcells; i++)
      cellCterY[i] = yStart + i*cellYsize;

    for (i = 0; i < param->nZcells; i++)
      cellCterZ[i] = zStart + i*cellZsize;
  }

/*
 *   Get the indices of the cell containing the field point

 *   Note: there is a change of cell index convention here
 *         The cellCharge array does not recognize padding of PBC image cells
 *         cellIndex[0] now goes from 0 to NCellX-1
 */

  //cellIndex[0] --;
  //cellIndex[1] --;
  //cellIndex[2] --;



/*
 *       Get the previous and next cell around cellIndex[X],Y and Z
 *       Wrap around if PBC are on.
 *
 *       Skips go from cell-1 to cell+1
 *       cells in interval [0, param->ncells - 1]
 */

  xSkip1 = cellIndex[0] - 1 ;
  if (xSkip1 < 0) {
    if (param->xBoundType == Periodic)
      xSkip1 = param->nXcells - 1 ;
    else
      xSkip1 = 0;
  }
  xSkip2 = cellIndex[0] ;
  xSkip3 = cellIndex[0] + 1 ;
  if (xSkip3 >= param->nXcells) {
    if (param->xBoundType == Periodic)
      xSkip3 = 0 ;
    else
      xSkip3 = param->nXcells - 1 ;
  }

  ySkip1 = cellIndex[1] - 1 ;
  if (ySkip1 < 0) {
    if (param->yBoundType == Periodic)
      ySkip1 = param->nYcells - 1 ;
    else
      ySkip1 = 0;
  }
  ySkip2 = cellIndex[1] ;
  ySkip3 = cellIndex[1] + 1 ;
  if (ySkip3 >= param->nYcells) {
    if (param->yBoundType == Periodic)
      ySkip3 = 0 ;
    else
      ySkip3 = param->nYcells - 1;
  }

  zSkip1 = cellIndex[2] - 1 ;
  if (zSkip1 < 0) {
    if (param->zBoundType == Periodic)
      zSkip1 = param->nZcells - 1 ;
    else
      zSkip1 = 0;
  }
  zSkip2 = cellIndex[2] ;
  zSkip3 = cellIndex[2] + 1 ;
  if (zSkip3 >= param->nZcells) {
    if (param->zBoundType == Periodic)
      zSkip3 = 0 ;
    else
      zSkip3 = param->nZcells - 1 ;
  }

#if deb
  printf("xSkip1=%d xSkip2=%d xSkip3=%d\n",xSkip1,xSkip2,xSkip3);
  printf("ySkip1=%d ySkip2=%d ySkip3=%d\n",ySkip1,ySkip2,ySkip3);
  printf("zSkip1=%d zSkip2=%d zSkip3=%d\n\n",zSkip1,zSkip2,zSkip3);
#endif

  for (cx = 0; cx < param->nXcells; cx++) {
    for (cy = 0; cy < param->nYcells; cy++) {
      for (cz = 0; cz < param->nZcells; cz++) {
	includePrimary = !(
			   (cx==xSkip1 || cx==xSkip2 || cx==xSkip3) &&
			   (cy==ySkip1 || cy==ySkip2 || cy==ySkip3) &&
			   (cz==zSkip1 || cz==zSkip2 || cz==zSkip3));


/*
 *              Get the center point of cell [cx, cy, cz]
 */
	xc = cellCterX[cx];
	yc = cellCterY[cy];
	zc = cellCterZ[cz];

/*
 *              Get the stress at the specified point caused
 *              by the net charge tensor of the current cell.
 */
	dx = xc - x;
	dy = yc - y;
	dz = zc - z;

	ZImage(param, &dx, &dy, &dz);

	xc = x + dx;
	yc = y + dy;
	zc = z + dz;


	cellID = cz + param->nZcells*cy +
	  param->nZcells*param->nYcells*cx;

 /*
  *              Stress (charge[.,1], [1,0,0])
  */
	burgX = home->cellCharge[9*cellID];
	burgY = home->cellCharge[9*cellID+3];
	burgZ = home->cellCharge[9*cellID+6];

	dx = 1;
	dy = 0;
	dz = 0;

	Matrix33_Zero(delSig);
	dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz,
		      burgX, burgY, burgZ, x, y, z,
		      includePrimary);

	for (i = 0; i < 3; i++)
	  for (j = 0; j < 3; j++)
	    remStress[i][j] += delSig[i][j];

	/*
	 *              Stress (charge[.,2], [0,1,0])
	 */
	burgX = home->cellCharge[9*cellID+1];
	burgY = home->cellCharge[9*cellID+4];
	burgZ = home->cellCharge[9*cellID+7];

	dx = 0;
	dy = 1;
	dz = 0;

	Matrix33_Zero(delSig);
	dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz,
		      burgX, burgY, burgZ, x, y, z,
		      includePrimary);

	for (i = 0; i < 3; i++)
	  for (j = 0; j < 3; j++)
	    remStress[i][j] += delSig[i][j];

/*
 *              Stress (charge[.,3], [0,0,1])
 */
	burgX = home->cellCharge[9*cellID+2];
	burgY = home->cellCharge[9*cellID+5];
	burgZ = home->cellCharge[9*cellID+8];

	dx = 0;
	dy = 0;
	dz = 1;

	Matrix33_Zero(delSig);
	dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz,
		      burgX, burgY, burgZ, x, y, z,
		      includePrimary);

	for (i = 0; i < 3; i++)
	  for (j = 0; j < 3; j++)
	    remStress[i][j] += delSig[i][j];

      } /* end for(cz = 0; ...) */
    } /* end for(cy = 0; ...) */
  } /* end for(cx = 0; ...) */

}
#endif

#ifndef FULL_N2_FORCES
/*      (newer more efficient version, use cells)
 *      Calculate stress at point (xm, ym, zm) due to
 *      local segments.
 */
void AllSegmentStress(Home_t *home, real8 x, real8 y, real8 z,
                      real8 Stress[3][3])
{
        int     i,j;
        int     cellIndex[3];
	real8   locStress[3][3],remStress[3][3];

	int DomID = -1;

	Param_t *param;

	param = home->param;

/*
 *      Initialize the different stresses components.
 */

        Matrix33_Zero(Stress);
        Matrix33_Zero(locStress);
        Matrix33_Zero(remStress);

        cellIndex[X] = (int)((x - param->minSideX) / (param->Lx / param->nXcells))+1;
        cellIndex[Y] = (int)((y - param->minSideY) / (param->Ly / param->nYcells))+1;
        cellIndex[Z] = (int)((z - param->minSideZ) / (param->Lz / param->nZcells))+1;

/*
 *   STRESS FROM LOCAL CELLS
 */
	LocalStress(home, x, y, z, cellIndex, locStress);

/*
 *   STRESS FROM REMOTE CELLS
 */

/*
 *       Only the processor encompassing a point should calculate the
 *       remote stress at that point.
 */

	DomID = isRightDom(home,x,y,z);

	if (DomID == home->myDomain)
	  {
	    if(param->fmEnabled)
	      {
		/* Compute the remote stress using FMM */
                 RemoteStressWithFMM(home,x,y,z,cellIndex,remStress);
	      }

	  } /*matches if Dom owns the seg */


/*
 * Sum up all stresses contributions
 *
 */
	for (i = 0; i < 3; i++)
	  for (j = 0; j < 3; j++)
	    {
	      Stress[i][j] = locStress[i][j] + remStress[i][j];
	    }
        return;
}

#else

void AllSegmentStress(Home_t *home,
                      real8 xm, real8 ym, real8 zm,
                      real8 totStress[3][3])
{
        int     i, mm, kk, nc2, ti2, includePrimary;
        real8   xA, yA, zA, xB, yB, zB;
        real8   xc, yc, zc;
        real8   MU, NU, bX, bY, bZ,a;
        real8   locstress[3][3], sigma[3][3], delSig[3][3];
	real8   rs[3],rm[3],b[3];
        Node_t  *node1, *node2;
        Param_t *param;

        param = home->param;

        MU = param->shearModulus;
        NU = param->pois;
        a  = param->rc;

/*
 *      Initialize stress to zero
 */
        for (mm = 0; mm < 3; mm++)
	  for (kk = 0; kk < 3; kk++)
	    locstress[mm][kk] = 0;


	if (param->fmEnabled)
	  Fatal("This part of AllSegmentStress cannot run with FMM enabled");


/*
 *      Loop through all segments owned by this domain
 */
        for (i = 0; i < home->newNodeKeyPtr; i++)
	  {
            node1 = home->nodeKeys[i];
            if (!node1) continue;

            nc2 = node1->numNbrs;

            for (ti2 = 0; ti2 < nc2; ti2++)
	      {
		node2 = GetNeighborNode(home, node1, ti2);
                if (!node2) continue;

                if (OrderNodes(node2, node1) != 1) continue;

                bX = node1->burgX[ti2];
                bY = node1->burgY[ti2];
                bZ = node1->burgZ[ti2];

                xB = node1->x;
                yB = node1->y;
                zB = node1->z;

		PBCPOSITION(param,xm,ym,zm,&xB,&yB,&zB);

                xA = node2->x;
                yA = node2->y;
                zA = node2->z;

                PBCPOSITION(param, xB,  yB,  zB,  &xA, &yA, &zA);

                dx = xA - xB;   dy = yA - yB;   dz = zA - zB;
                xc = (xA+xB)/2; yc = (yA+yB)/2; zc = (zA+zB)/2;

                PBCPOSITION(param, xm,  ym,  zm,  &xc, &yc, &zc);

                xA = xc+dx/2;   yA = yc+dy/2;   zA = zc+dz/2;
                xB = xc-dx/2;   yB = yc-dy/2;   zB = zc-dz/2;


   	        Matrix33_Zero(delSig);
#ifndef FULL_N2_FORCES
                /* add PBC image stress (virtual segment not counted) */
                includePrimary = 0;
	        dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz,
	                      bX, bY, bZ, xm, ym, zm,
	                      includePrimary);
#endif


		for (mm = 0; mm < 3; mm++)
		  for (kk = 0; kk < 3; kk++)
		    sigma[mm][kk] = 0.0;

		SegmentStress(MU, NU, bX, bY, bZ,
			      xB, yB, zB, xA, yA, zA,
			      xm, ym, zm, a, sigma);

                for (mm = 0; mm < 3; mm++)
                    for (kk = 0; kk < 3; kk++)
                        locstress[mm][kk] += sigma[mm][kk] +
                           delSig[mm][kk];


	      }  /* end for (ti2 = 0; ...) */
	  }  /* end for (i = 0;...) */

/*
 *      For serial runs, the local stress field is the complete
 *      stress field for the problem, but for parallel applications
 *      we have to sum the local stress fields from all processors
 *      into the total stress field.
 */
#ifdef PARALLEL
        MPI_Allreduce(locstress, totStress, 9, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
#else
        for (mm = 0; mm < 3; mm++)
            for (kk = 0; kk < 3; kk++)
                totStress[mm][kk] = locstress[mm][kk];
#endif
        return;
}

#endif


#endif
