/**************************************************************************
 *
 *  Function    : Mobility_FCC_0
 *  Author      : Wei Cai, Seok-Woo Lee (updated 07/14/09)
 *  Description : Generic Mobility Law of FCC metals
 *                Each line has a glide plane from input
 *                and it never changes
 *                If the plane normal is not of {111} type, dislocation
 *                motion is constrained along line direction
 *                If node flag == 7, node velocity is zero
 *
 *  Arguments:
 *      IN:  node        Pointer to the node for which to
 *                       calculate velocity
 *      IN/OUT: mobArgs  Structure containing additional
 *                       parameters for conveying information
 *                       to/from the mobility function.
 *
 *  Returns:  0 on success
 *            1 if velcoity could not be determined
 *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

#define _ENABLE_LINE_CONSTRAINT 1

int Mobility_FCC_0(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
    int numNonZeroLenSegs = 0;
    Param_t *param;
    real8 nVel[3];
    int i, j, nc, nlc;
    real8 lineX[MAX_NBRS], lineY[MAX_NBRS], lineZ[MAX_NBRS];
    real8 b;
    real8 dx, dy, dz, lx, ly, lz, lr, LtimesB;
    Node_t *nbr;
    real8 MobScrew, MobEdge, Mob;
    real8 bx, by, bz, br, dangle;
    real8 nForce[3];

    param = home->param;

    MobScrew = param->MobScrew;
    MobEdge  = param->MobEdge;
    
    nc = node->numNbrs;
    
    real8 (*invDragMatrix)[3]    = mobArgs->invDragMatrix;
 // real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
    int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

    Matrix33_Zero(invDragMatrix);
    *numGlideConstraints= 0;

/*
 *  If node is 'pinned' in ALL dimensions, or the node has any arms 
 *  with a burgers vector that has explicitly been set to be sessile (via
 *  control file inputs), the node may not be moved so just zero the velocity
 *  and return
 */
    if (HAS_ALL_OF_CONSTRAINTS(node->constraint, PINNED_NODE) ||
        NodeHasSessileBurg(home, node))
    {
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }

/*
 *  It's possible this function was called for a node which had only zero-
 *  length segments (during SplitSurfaceNodes() for example).  If that is
 *  the case, just set the velocity to zero and return.
 */
    for (i = 0; i < nc; i++) {
        if ((nbr = GetNeighborNode(home, node, i)) == (Node_t *)NULL) continue;
        dx = node->x - nbr->x;
        dy = node->y - nbr->y;
        dz = node->z - nbr->z;
        if ((dx*dx + dy*dy + dz*dz) > 1.0e-12) {
            numNonZeroLenSegs++;
        }
    }

    if (numNonZeroLenSegs == 0) {
        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;
        return(0);
    }


/*
 *  copy glide plane constraints and determine line constraints
 */
    for(i=0;i<nc;i++) {
        real8 norm[3];

        norm[X] = node->nx[i];
        norm[Y] = node->ny[i];
        norm[Z] = node->nz[i];

/*
 *      If needed, rotate the glide plane normal from the laboratory
 *      frame to the crystalographic frame.
 */
        if (param->useLabFrame) {
            real8 normRot[3];

            Matrix33Vector3Multiply(param->rotMatrixInverse, norm, normRot);
            VECTOR_COPY(norm, normRot);
        }

        if ( (fabs(fabs(norm[X]) - fabs(norm[Y])) > FFACTOR_NORMAL) ||
             (fabs(fabs(norm[Y]) - fabs(norm[Z])) > FFACTOR_NORMAL) )
	  { /* not {111} plane */
            if ((nbr=GetNeighborNode(home,node,i)) == (Node_t *)NULL) {
	      Fatal("Neighbor not found at %s line %d\n",__FILE__,__LINE__);
            }
            lineX[i] = nbr->x - node->x;
            lineY[i] = nbr->y - node->y; 
            lineZ[i] = nbr->z - node->z;
            ZImage (param, lineX+i, lineY+i, lineZ+i);
	    
/*
 *          If needed, rotate the line sense from the laboratory frame to
 *          the crystal frame.
 */
            if (param->useLabFrame) {
	      real8 lDir[3] = {lineX[i], lineY[i], lineZ[i]};
	      real8 lDirRot[3];
	      
	      Matrix33Vector3Multiply(param->rotMatrixInverse, lDir, lDirRot);
	      
	      lineX[i] = lDirRot[0];
	      lineY[i] = lDirRot[1];
	      lineZ[i] = lDirRot[2];
            }
	  }
	else
	  { /* no line constraint */
	    lineX[i] = lineY[i] = lineZ[i] = 0;
	  }
      }
    
    /* normalize lc line vectors*/
    for(i=0;i<nc;i++)
    {
	b=sqrt(lineX[i]*lineX[i]+lineY[i]*lineY[i]+lineZ[i]*lineZ[i]);

        if(b>0)
        {
            lineX[i]/=b;
            lineY[i]/=b;
            lineZ[i]/=b;
        }
    }

    /* Find independent line constraints */
    nlc = 0;
    for(i=0;i<nc;i++)
    {
        for(j=0;j<i;j++)
        {
            Orthogonalize(lineX+i,lineY+i,lineZ+i,lineX[j],lineY[j],lineZ[j]);
        }
        if((lineX[i]*lineX[i]+lineY[i]*lineY[i]+lineZ[i]*lineZ[i])<FFACTOR_ORTH)
        {
            lineX[i] = lineY[i] = lineZ[i] = 0;
        }
        else
        {
            nlc ++;
        }
    }

    /* find total dislocation length times drag coefficent (LtimesB)*/
    LtimesB=0;
    for(j=0;j<nc;j++)
    {
        if ((nbr=GetNeighborNode(home,node,j)) == (Node_t *)NULL) continue;
        dx=nbr->x - node->x;
        dy=nbr->y - node->y;
        dz=nbr->z - node->z;
        ZImage (param, &dx, &dy, &dz) ;

/*
 *      If needed, rotate the line sense from the laboratory frame to
 *      the crystal frame.
 */
        if (param->useLabFrame) {
            real8 dTmp[3] = {dx, dy, dz};
            real8 dRot[3];

            Matrix33Vector3Multiply(param->rotMatrixInverse, dTmp, dRot);

            dx = dRot[0]; dy = dRot[1]; dz = dRot[2];
        }

        lr=sqrt(dx*dx+dy*dy+dz*dz);
        
        if (lr==0)
        { /* zero arm segment can happen after node split 
           * it is OK to have plane normal vector == 0
           * Skip (do nothing)
           */
        }
        else 
        {
           if((node->nx[j]==0)&&(node->ny[j]==0)&&(node->nz[j]==0))
           {
              printf("Mobility_FCC_0: (%d,%d) glide plane norm = 0\n"
                     "for segment with nonzero length lr = %e!\n",
                     node->myTag.domainID, node->myTag.index, lr);
           }

           lx=dx/lr; ly=dy/lr; lz=dz/lr;

           bx = node->burgX[j];
           by = node->burgY[j];
           bz = node->burgZ[j];
/*
 *         If needed, rotate the burgers vector from the laboratory frame to
 *         the crystal frame.
 */
           if (param->useLabFrame) {
               real8 bTmp[3] = {bx, by, bz};
               real8 bRot[3];

               Matrix33Vector3Multiply(param->rotMatrixInverse, bTmp, bRot);

               bx = bRot[0]; by = bRot[1]; bz = bRot[2];
           }

           br = sqrt(bx*bx+by*by+bz*bz);
           bx/=br; by/=br; bz/=br; /* unit vector along Burgers vector */

           dangle = fabs(bx*lx+by*ly+bz*lz);
           Mob=MobEdge+(MobScrew-MobEdge)*dangle;

           LtimesB+=(lr/Mob);
	}
    }
    LtimesB/=2;

    nForce[0] = node->fX;
    nForce[1] = node->fY;
    nForce[2] = node->fZ;

/*
 *  If needed, rotate the force vector from the laboratory frame to the
 *  crystal frame
 */
    if (param->useLabFrame) {
        real8 rotForce[3];
        Matrix33Vector3Multiply(param->rotMatrixInverse, nForce, rotForce);
        VECTOR_COPY(nForce, rotForce);
    }
    
    /* Velocity is simply proportional to total force per unit length */
    nVel[0] = nForce[0]/LtimesB;
    nVel[1] = nForce[1]/LtimesB;
    nVel[2] = nForce[2]/LtimesB;
    
    invDragMatrix[0][0] = 1./ LtimesB;
    invDragMatrix[1][1] = 1./ LtimesB;
    invDragMatrix[2][2] = 1./ LtimesB;

/*
 *  If needed, rotate the velocity vector back to the laboratory frame
 *  from the crystal frame
 */
    if (param->useLabFrame) {
        real8 vTmp[3] = {nVel[0], nVel[1], nVel[2]};
        real8 vRot[3];

        Matrix33Vector3Multiply(param->rotMatrix, vTmp, vRot);
        VECTOR_COPY(nVel, vRot);
    }

    node->vX = nVel[0];
    node->vY = nVel[1];
    node->vZ = nVel[2];

    return(0);
}
