/**************************************************************************
 *
 *      Module:       Mobility_BCC_glide
 *      Description:  Contains functions for calculating mobility of nodes
 *                    in BCC metals.  Based on MobilityLaw_BCC_0b.c
 *                    This is really only slightly modified from
 *                    mobility BCC_0b by first obtaining the velocity assuming 
 *                    climbing exists, then project out the climb component 
 *                    from the velocity. 
 *
 *      Authors:      Meijie Tang         Nov. 1, 2007
 *
 *      Includes functions:
 *
 *            Mobility_BCC_glide()
 *                
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <sys/param.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

//#define M33_INVERSE(a,b) Matrix33_Inverse(a,b)
//#define M33_INVERSE(a,b) Matrix33_SVD_Inverse(a,b)
  #define M33_INVERSE(a,b) Matrix33_PseudoInverse(a,b)


/**************************************************************************
 *
 *      Function:     Mobility_BCC_glide
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Arguments:
 *          IN:  node          Pointer to the node for which to
 *                             calculate velocity
 *          IN/OUT: mobArgs    Structure containing additional
 *                             parameters for conveying information
 *                             to/from the mobility function.
 *
 *      Returns: 0 on success
 *               1 if velocity could not be determiend
 * 
 *************************************************************************/
int Mobility_BCC_glide(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int     i, nbrs;
        int     numNonZeroLenSegs = 0;
        real8   bx, by, bz;
        real8   dx, dy, dz;
        real8   mx, my, mz;
        real8   nx, ny, nz; 
        real8   mag, mag2, halfMag, invMag;
        real8   invbMag2, bMag2, costheta, costheta2;
        real8   Bline, Bscrew, Bedge, Bglide, Bclimb, Beclimb;
        real8   Bscrew2, Beclimb2;
        real8   invBscrew2, invBedge2;
        real8   BlmBsc, BclmBsc, BglmBsc, BlmBecl;
        real8   eps = 1.0e-12;
        real8   burgCryst[3];
        real8   nForce[3], nVel[3];
        real8   Btotal[3][3] = {{0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0},
                                {0.0, 0.0, 0.0}};

        Node_t  *nbrNode;
        Param_t *param;

        real8 (*invDragMatrix)[3]    =  mobArgs->invDragMatrix;
     // real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invDragMatrix);
        *numGlideConstraints = 0;

/*
 *      If the node is pinned in ALL dimensions it cannot be moved
 *      so just zero the velocity and return
 */
        if (HAS_ALL_OF_CONSTRAINTS(node->constraint, PINNED_NODE)) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      If any of the arms of the node has a burgers vector that
 *      has explicitly been set to be sessile (via control file
 *      inputs), the node may not be moved.
 */
        if (NodeHasSessileBurg(home, node)) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

        param = home->param;

        Bscrew     = 1.0 / param->MobScrew;
        Bedge      = 1.0 / param->MobEdge;

        /* Beclimb    = 1.0 / param->MobClimb; */ 
        Beclimb    = 1.0e10; 
        

        Bscrew2    = Bscrew * Bscrew;
        Beclimb2   = Beclimb * Beclimb;

        Bline      = 1.0e-2 * MIN(Bscrew, Bedge);
        BlmBsc     = Bline - Bscrew;
        BlmBecl    = Bline - Beclimb; 

        invBscrew2 = 1.0 / (Bscrew*Bscrew);
        invBedge2  = 1.0 / (Bedge*Bedge);

        nbrs = node->numNbrs;

        node->vX = 0.0e0;
        node->vY = 0.0e0;
        node->vZ = 0.0e0;

/*
 *      Loop over all arms of the node, adding each arm's contribution
 *      to the drag matrix.
 */
	for (i = 0; i < nbrs; i++) {
	  
	  bx = node->burgX[i];
	  by = node->burgY[i];
	  bz = node->burgZ[i];
	  
	  bMag2 = (bx*bx + by*by + bz*bz);
	  invbMag2 = 1.0 / bMag2;
	  
/*
 *        Calculate the length of the arm and its tangent line direction
 */
	  nbrNode = GetNeighborNode(home, node, i);
	  
	  if (nbrNode == (Node_t *)NULL) {
	    printf("WARNING: Neighbor not found at %s line %d\n",
		   __FILE__, __LINE__);
	    continue;
	  }
	  
	  dx = nbrNode->x - node->x;
	  dy = nbrNode->y - node->y;
	  dz = nbrNode->z - node->z;
	  
	  ZImage(param, &dx, &dy, &dz);

	  mag2    = dx*dx + dy*dy + dz*dz;
/*
 *        If the segment is zero length (which can happen when
 *        the mobility function is being called from SplitMultiNodes())
 *        just skip the segment.
 */
	  if (mag2 < eps) {
	    continue;
	  }
	  
	  numNonZeroLenSegs++;
	  
	  mag     = sqrt(mag2);
	  halfMag = mag/2.0;
	  invMag  = 1.0 / mag;
	  
	  dx *= invMag;
	  dy *= invMag;
	  dz *= invMag;

/*
 *        Calculate how close to screw the arm is
 */
	  costheta = (dx*bx + dy*by + dz*bz);
	  costheta2 = (costheta*costheta) * invbMag2;

/*
 *        [0 0 1] arms don't move as readily as other arms, so must be
 *        handled specially.
 *
 *        If needed, rotate a copy of the burgers vector from the
 *        laboratory frame to the crystal frame.
 */
	  if (param->useLabFrame) {
	    real8 bTmp[3] = {bx, by, bz};
	    Matrix33Vector3Multiply(param->rotMatrixInverse, bTmp,
				    burgCryst);
	  } else {
	    burgCryst[X] = bx;
	    burgCryst[Y] = by;
	    burgCryst[Z] = bz;
	  }
	  
	  if (fabs(burgCryst[X]*burgCryst[Y]*burgCryst[Z]) < eps) {
	    if (nbrs == 2) {
	      Btotal[0][0] += halfMag * Beclimb;
	      Btotal[1][1] += halfMag * Beclimb;
	      Btotal[2][2] += halfMag * Beclimb;
	    } else {
	      Btotal[0][0] += halfMag * (dx*dx * BlmBecl + Beclimb);
	      Btotal[0][1] += halfMag * (dx*dy * BlmBecl);
	      Btotal[0][2] += halfMag * (dx*dz * BlmBecl);
	      Btotal[1][1] += halfMag * (dy*dy * BlmBecl + Beclimb);
	      Btotal[1][2] += halfMag * (dy*dz * BlmBecl);
	      Btotal[2][2] += halfMag * (dz*dz * BlmBecl + Beclimb);
	    }
	  } else  {
/*
 *                  Arm is not [0 0 1], so build the drag matrix assuming the
 *                  dislocation is screw type
 */
	    Btotal[0][0] += halfMag * (dx*dx * BlmBsc + Bscrew);
	    Btotal[0][1] += halfMag * (dx*dy * BlmBsc);
	    Btotal[0][2] += halfMag * (dx*dz * BlmBsc);
	    Btotal[1][1] += halfMag * (dy*dy * BlmBsc + Bscrew);
	    Btotal[1][2] += halfMag * (dy*dz * BlmBsc);
	    Btotal[2][2] += halfMag * (dz*dz * BlmBsc + Bscrew);

/*
 *                  Now correct the drag matrix for dislocations that are
 *                  not screw
 */
	    if ((1.0 - costheta2) > eps) {
	      
	      
	      nx = node->nx[i];
	      ny = node->ny[i];
	      nz = node->nz[i];
	      
	      xvector(nx, ny, nz, dx, dy, dz, &mx, &my, &mz);
	      
	      Bglide = sqrt(invBedge2 + (invBscrew2-invBedge2) *
			    costheta2);
	      Bglide = 1.0 / Bglide;
	      Bclimb = sqrt(Beclimb2 + (Bscrew2 - Beclimb2) *
			    costheta2);

#ifdef NAN_CHECK
	      if (isnan(Bglide) != 0) {
		Fatal("Mobility_BCC_glide: Bglide = NaN\n"
		      "  Bglide = sqrt(invBedge2 + "
		      "(invBscrew2-invBedge2)*costheta2)\n"
		      "  where invBedge2 = %lf, invBscrew2 = %lf, "
		      "costheta2 = %lf", invBedge2, invBscrew2,
		      costheta2);
	      }
	      
	      if (isnan(Bclimb) != 0) {
		Fatal("Mobility_BCC_glide: Bclimb = NaN\n"
		      "  Bclimb = sqrt(Beclimb2 + "
		      "(Bscrew2-Beclimb2)*costheta2)\n"
		      "  where Beclimb2 = %lf, Bscrew2 = %lf, "
		      "costheta2 = %lf", Beclimb2, Bscrew2,
		      costheta2);
	      }
#endif
	      BclmBsc = Bclimb - Bscrew;
	      BglmBsc = Bglide - Bscrew;
	      
	      
	      Btotal[0][0] += halfMag * (nx*nx * BclmBsc +
					 mx*mx * BglmBsc);
	      Btotal[0][1] += halfMag * (nx*ny * BclmBsc +
					 mx*my * BglmBsc);
	      Btotal[0][2] += halfMag * (nx*nz * BclmBsc +
					 mx*mz * BglmBsc);
	      Btotal[1][1] += halfMag * (ny*ny * BclmBsc +
					 my*my * BglmBsc);
	      Btotal[1][2] += halfMag * (ny*nz * BclmBsc +
					 my*mz * BglmBsc);
	      Btotal[2][2] += halfMag * (nz*nz * BclmBsc +
					 mz*mz * BglmBsc);
	    }
	  }  /* End non-[0 0 1] arm */
	}  /* End loop over arms */
	
	Btotal[1][0] = Btotal[0][1];
	Btotal[2][0] = Btotal[0][2];
	Btotal[2][1] = Btotal[1][2];

/*
 *          It's possible this function was called for a node which only
 *          had zero length segments (during SplitSurfaceNodes() for example).
 *          If that is the case, just set the velocity to zero and return;
 */
	if (numNonZeroLenSegs == 0) {
	  node->vX = 0.0;
	  node->vY = 0.0;
	  node->vZ = 0.0;
	  return(0);
	}
/*
 *          At this point we should check if the matrix is invertable and
 *          if it isn't, find the eigen values and eigen vectors of the drag
 *          matrix, then invert the drag matrix keeping zero eigen values
 *          as zero.
 *
 *          FIX ME!  For now, we're assuming the matrix is invertible.
 */
        
        nForce[0] = node->fX;
        nForce[1] = node->fY;
        nForce[2] = node->fZ;
	
        if ( M33_INVERSE(invDragMatrix, Btotal) < 0 )
        { Fatal("%s::%s(%d) : Cannot invert drag matrix!", __FILE__, __func__, __LINE__ ); }

        Matrix33Vector3Multiply(invDragMatrix, nForce, nVel);
 
        node->vX = nVel[0];
        node->vY = nVel[1];
        node->vZ = nVel[2];

        return(0);
}
