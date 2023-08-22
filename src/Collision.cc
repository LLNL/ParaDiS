/*****************************************************************************
 *
 *      Module:         Collision.c
 *      Description:    This module contains a dispatch function to call the
 *                      user-selected collision handling mechanism plus
 *                      various support functions that are used by all the
 *                      collision handling mechanisms.
 *
 *                      Each domain handles local collisions and
 *                      distributes the necessary topological changes
 *                      to remote domains via the same mechanism
 *                      employed by remesh.
 *
 *                      NOTE: Certain strict rules govern what topological
 *                      changes are permitted based on noda and segment
 *                      ownership.  See comments at the beginning of
 *                      the module Topology.c for details of the
 *                      rule.  Additional restrictions may be implemented
 *                      for the collision handling; see code below
 *                      for details of any such restrictions.
 *
 *
 *      Included functions:
 *
 *          FindMergeNodeOnSegment()
 *          FindCollisionPoint()
 *          GetMinDist()
 *          HandleCollisions()
 *
 *****************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"

#ifdef DEBUG_TOPOLOGY_CHANGES
static int dbgDom = 0;
#endif

int FindMergeNodeOnSegment(Home_t *home, Node_t **node30, Node_t **node40, 
                           real8 x3, real8 y3, real8 z3,
                           real8 x4, real8 y4, real8 z4, int arm34,
                           real8 L2, Node_t **mergenode2)
{
   int globalOp = 1, arm43, splitStatus=0, q;
   Tag_t   oldTag1, oldTag2;
   real8   vx3, vy3, vz3;
   real8   vx4, vy4, vz4;
   real8   newx, newy, newz;
   real8   newvx, newvy, newvz;

   real8 burg1[3], oldfp0[3], oldfp1[3], f0seg1[3], f1seg1[3], f0seg2[3], f1seg2[3];
   real8 pnew[3], nodeVel[3],newNodeVel[3];

   real8   pos0[3], pos1[3];
   Node_t  *splitNode1, *splitNode2, *tmpNbr;

   Param_t *param;
   MobArgs_t mobArgs;


   Node_t *node3, *node4;

   node3 = *node30;
   node4 = *node40;

   param = home->param;

   *mergenode2 = (Node_t *)NULL;

   vx3 = node3->vX; vy3 = node3->vY; vz3 = node3->vZ;
   vx4 = node4->vX; vy4 = node4->vY; vz4 = node4->vZ;


   newx = x3 * (1.0-L2) + x4*L2;
   newy = y3 * (1.0-L2) + y4*L2;
   newz = z3 * (1.0-L2) + z4*L2;

   newvx = vx3 * (1.0-L2) + vx4*L2;
   newvy = vy3 * (1.0-L2) + vy4*L2;
   newvz = vz3 * (1.0-L2) + vz4*L2;

/*
 *                                   Estimate resulting forces on all segments
 *                                   involved in the split.
 */
   arm43 = GetArmID(node4, node3);
   
   burg1[X] = node3->burgX[arm34];
   burg1[Y] = node3->burgY[arm34];
   burg1[Z] = node3->burgZ[arm34];
   
   oldfp0[X] = node3->armfx[arm34];
   oldfp0[Y] = node3->armfy[arm34];
   oldfp0[Z] = node3->armfz[arm34];
   
   oldfp1[X] = node4->armfx[arm43];
   oldfp1[Y] = node4->armfy[arm43];
   oldfp1[Z] = node4->armfz[arm43];
   
   pos0[X] = x3;   pos1[X] = x4;
   pos0[Y] = y3;   pos1[Y] = y4;
   pos0[Z] = z3;   pos1[Z] = z4;
   
   pnew[X] = newx;
   pnew[Y] = newy;
   pnew[Z] = newz;
   
   newNodeVel[X] = newvx;
   newNodeVel[Y] = newvy;
   newNodeVel[Z] = newvz;
   
   nodeVel[X] = vx3;
   nodeVel[Y] = vy3;
   nodeVel[Z] = vz3;
   
   FindSubFSeg(home, pos0, pos1, burg1, oldfp0,
               oldfp1, pnew, f0seg1, f1seg1,
               f0seg2, f1seg2);
   
   oldTag1 = node3->myTag;
   oldTag2 = node4->myTag;
   
   FoldBox(param, &pnew[X], &pnew[Y], &pnew[Z]);
   
   splitStatus = SplitNode(home,
                           OPCLASS_COLLISION,
                           node3, pos0,
                           pnew, nodeVel,
                           newNodeVel, 1,
                           &arm34, globalOp,
                           &splitNode1,
                           &splitNode2, 0);

/*
 *                                   If we were unable to split the node
 *                                   go back to looking for more collision
 *                                   candidates.
 */
   if (splitStatus == SPLIT_FAILED) {
      return (splitStatus);
   }

/*
 *                                   The force estimates above are good enough
 *                                   for the remainder of this timestep, but
 *                                   mark the force and velocity data for some
 *                                   nodes as obsolete so that more accurate
 *                                   forces will be recalculated either at the
 *                                   end of this timestep, or the beginning of
 *                                   the next.
 */
   MarkNodeForceObsolete(home, splitNode2);
   
   for (q = 0; q < splitNode2->numNbrs; q++) {
      tmpNbr = GetNodeFromTag(home, splitNode2->nbrTag[q]);
      if (tmpNbr != (Node_t *)NULL) {
         tmpNbr->flags |= NODE_RESET_FORCES;
      }
   }
   
/*
 *                                   Reset nodal forces on nodes involved in the
 *                                   split.
 */
   ResetSegForces(home, splitNode1,
                  &splitNode2->myTag,
                  f0seg1[X], f0seg1[Y],
                  f0seg1[Z], 1);
   
   ResetSegForces(home, splitNode2,
                  &splitNode1->myTag,
                  f1seg1[X], f1seg1[Y],
                  f1seg1[Z], 1);
   
   ResetSegForces(home, splitNode2,
                  &node4->myTag,
                  f0seg2[X], f0seg2[Y],
                  f0seg2[Z], 1);
   
   ResetSegForces(home, node4,
                  &splitNode2->myTag,
                  f1seg2[X], f1seg2[Y],
                  f1seg2[Z], 1);
   
   (void)EvaluateMobility(home, splitNode1,
                          &mobArgs);
   
   (void)EvaluateMobility(home, splitNode2,
                          &mobArgs);
   
   (void)EvaluateMobility(home, node4,
                          &mobArgs);
   
/*
 *                                   When debugging, dump some info on
 *                                   topological changes taking place and
 *                                   the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
   if ((dbgDom < 0) ||
       (dbgDom == home->myDomain)) {
      printf("  Split-2(SegCollision): "
             "(%d,%d)--(%d,%d) ==> "
             "(%d,%d)--(%d,%d)--(%d,%d)\n",
             oldTag1.domainID,
             oldTag1.index,
             oldTag2.domainID,
             oldTag2.index,
             splitNode1->myTag.domainID,
             splitNode1->myTag.index,
             splitNode2->myTag.domainID,
             splitNode2->myTag.index,
             node4->myTag.domainID,
             node4->myTag.index);
      PrintNode(splitNode1);
      PrintNode(splitNode2);
      PrintNode(node4);
   }
#endif


   *node30 = node3;
   *node40 = node4;
   (*mergenode2) = splitNode2;
   return (splitStatus);
}



/*---------------------------------------------------------------------------
 *
 *      Function:       FindCollisionPoint
 *      Description:    This function attempts to select a collision
 *                      point on a plane common to the two nodes.
 *
 *      Arguments:
 *          node1, node2   pointers to the two node structures
 *          x, y, z        pointers to locations in which to return
 *                         the coordinates of the point at which the
 *                         two nodes should be collided.
 *
 *      Returns:  0 on error
 *                1 if a collision point was successfully found
 *
 *-------------------------------------------------------------------------*/
int FindCollisionPoint(Home_t *home, Node_t *node1, Node_t *node2,
                       real8 *x, real8 *y, real8 *z)
{
        int     i, j, m, n;
        real8   L, invL, tmp;
        real8   norm, invnorm;
        real8   n1mag2, n2mag2, eps;
        real8   dx, dy, dz;
        real8   n1x, n1y, n1z;
        real8   n2x, n2y, n2z;
        real8   dirx, diry, dirz;
        real8   p1[3], p2[3];
        real8   plane[3], vector[3];
        real8   tmp33[3][3];
        Node_t  *nbrNode;
        Param_t *param;

        param = home->param;

        eps = 1.0e-12;

        p1[0] = node1->x;
        p1[1] = node1->y;
        p1[2] = node1->z;

        p2[0] = node2->x;
        p2[1] = node2->y;
        p2[2] = node2->z;

        PBCPOSITION(param, p1[0], p1[1], p1[2], &p2[0], &p2[1], &p2[2]);

/*
 *      If a node is a 'pinned' node it can't be relocated, so use
 *      that node's coordinates as the collision point.
 *
 *      NOTE: For now, we treat a node that is pinned in ANY dimension
 *            as if it is pinned in ALL dimensions.  This may change
 *            in the future.
 */
        if (HAS_ANY_OF_CONSTRAINTS(node1->constraint, PINNED_NODE)) {
            *x = p1[0];
            *y = p1[1];
            *z = p1[2];
            return(1);
        } else if (HAS_ANY_OF_CONSTRAINTS(node2->constraint, PINNED_NODE)) {
            *x = p2[0];
            *y = p2[1];
            *z = p2[2];
            return(1);
        }

/*
 *      If strict glide planes are not being enforced, we
 *      can use 'fuzzy' glide planes so that we don't
 *      throw nodes long distances if the planes are
 *      close to parallel but don't quite meet the criteria
 *      for being parallel.
 */
        if (param->enforceGlidePlanes == 0) {

            real8   factor;
            real8   vector[3], result[3];
            real8   Matrix[3][3], invMatrix[3][3];

/*
 *          <factor> determines how 'fuzzy' planes are.  The smaller the
 *          value, the fuzzier the plane constraints.
 */
            factor = 5.0;

            vector[0] = (p1[0] + p2[0]) / 2.0;
            vector[1] = (p1[1] + p2[1]) / 2.0;
            vector[2] = (p1[2] + p2[2]) / 2.0;

            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    Matrix[m][n] = (real8)(m == n);
                }
            }

            for (i = 0; i < node1->numNbrs; i++) {

                nbrNode = GetNeighborNode(home, node1, i);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                dx = p1[0] - nbrNode->x;
                dy = p1[1] - nbrNode->y;
                dz = p1[2] - nbrNode->z;

                ZImage(param, &dx, &dy, &dz);

                L = sqrt(dx*dx + dy*dy + dz*dz);
                invL = 1.0 / L;
                dirx = dx * invL;
                diry = dy * invL;
                dirz = dz * invL;

                xvector(dirx, diry, dirz, node1->burgX[i], node1->burgY[i],
                        node1->burgZ[i], &n1x, &n1y, &n1z);

                xvector(dirx, diry, dirz, node1->vX, node1->vY, node1->vZ,
                        &n2x, &n2y, &n2z);

                n1mag2 = n1x*n1x + n1y*n1y + n1z*n1z;
                n2mag2 = n2x*n2x + n2y*n2y + n2z*n2z;

                if (n1mag2 > eps) {
                    norm = sqrt(n1mag2);
                    invnorm = 1.0 / norm;
                    plane[0] = n1x * invnorm;
                    plane[1] = n1y * invnorm;
                    plane[2] = n1z * invnorm;
                } else if (n2mag2 > eps) {
                    norm = sqrt(n2mag2);
                    invnorm = 1.0 / norm;
                    plane[0] = n2x * invnorm;
                    plane[1] = n2y * invnorm;
                    plane[2] = n2z * invnorm;
                } else continue;
                
                Matrix33_VVt_Mul(tmp33,plane, plane);

                for (m = 0; m < 3; m++) {
                    for (n = 0; n < 3; n++) {
                        Matrix[m][n] += factor * tmp33[m][n];
                    }
                }

                tmp = DotProduct(plane, p1);

                for (m = 0; m < 3; m++) {
                    vector[m] += factor * tmp * plane[m];
                }
            }

            for (i = 0; i < node2->numNbrs; i++) {

                nbrNode = GetNeighborNode(home, node2, i);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                dx = p2[0] - nbrNode->x;
                dy = p2[1] - nbrNode->y;
                dz = p2[2] - nbrNode->z;

                ZImage(param, &dx, &dy, &dz);

                L = sqrt(dx*dx + dy*dy + dz*dz);
                invL = 1.0 / L;
                dirx = dx * invL;
                diry = dy * invL;
                dirz = dz * invL;

                xvector(dirx, diry, dirz, node2->burgX[i], node2->burgY[i],
                        node2->burgZ[i], &n1x, &n1y, &n1z);

                xvector(dirx, diry, dirz, node2->vX, node2->vY, node2->vZ,
                        &n2x, &n2y, &n2z);

                n1mag2 = n1x*n1x + n1y*n1y + n1z*n1z;
                n2mag2 = n2x*n2x + n2y*n2y + n2z*n2z;

                if (n1mag2 > eps) {
                    norm = sqrt(n1mag2);
                    invnorm = 1.0 / norm;
                    plane[0] = n1x * invnorm;
                    plane[1] = n1y * invnorm;
                    plane[2] = n1z * invnorm;
                } else if (n2mag2 > eps) {
                    norm = sqrt(n2mag2);
                    invnorm = 1.0 / norm;
                    plane[0] = n2x * invnorm;
                    plane[1] = n2y * invnorm;
                    plane[2] = n2z * invnorm;
                } else continue;
                
                Matrix33_VVt_Mul(tmp33,plane, plane);

                for (m = 0; m < 3; m++) {
                    for (n = 0; n < 3; n++) {
                        Matrix[m][n] += factor * tmp33[m][n];
                    }
                }

                tmp = DotProduct(plane, p2);
    
                for (m = 0; m < 3; m++) {
                    vector[m] += factor * tmp * plane[m];
                }
            }

            if ( Matrix33_Inverse(invMatrix, Matrix) < 0 ) 
            {
                printf("%s::%s(%d) : Unable to invert 3x3 matrix!\n", __FILE__, __func__, __LINE__ );
#if 0
                printf("FindCollisionPoint(): node1 = (%d,%d), node2 = (%d,%d)\n",
                       node1->myTag.domainID, node1->myTag.index,
                       node2->myTag.domainID, node2->myTag.index);
                for (m = 0; m < 3; m++) {
                    for (n = 0; n < 3; n++) {
                        printf("    Matrix[%d][%d] = %e\n", m, n, Matrix[m][n]);
                    }
                }
                PrintNode(node1);
                PrintNode(node2);
#endif
                fflush(NULL);
                return(0);
            }

            Matrix33Vector3Multiply(invMatrix, vector, result);

            *x = result[0];
            *y = result[1];
            *z = result[2];

        } else {
/*
 *         enforceGlidePlanes == 1, so strict glide planes are being used...
 */
           int     conditionsmet, Nsize;
           real8   newplanecond, npc2, onemnpc4, detN;
           real8   Nmat[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
           real8   Matrix[6][6], invMatrix[6][6];
           real8   V[6], result[6];

           newplanecond = 0.875;
           npc2         = newplanecond * newplanecond;

           tmp          = 1.0 - newplanecond;
           onemnpc4     = tmp * tmp * tmp * tmp;

           vector[0] = 0.0;
           vector[1] = 0.0;
           vector[2] = 0.0;

           Nsize = 0;

           for (i = 0; i < node1->numNbrs; i++) {

               if (Nsize < 3) {

                   nbrNode = GetNeighborNode(home, node1, i);

                   if (nbrNode == (Node_t *)NULL) {
                       printf("WARNING: Neighbor not found at %s line %d\n",
                              __FILE__, __LINE__);
                       continue;
                   }

                   dx = p1[0] - nbrNode->x;
                   dy = p1[1] - nbrNode->y;
                   dz = p1[2] - nbrNode->z;

                   ZImage(param, &dx, &dy, &dz);

                   L = sqrt(dx*dx + dy*dy + dz*dz);
                   invL = 1.0 / L;
                   dirx = dx * invL;
                   diry = dy * invL;
                   dirz = dz * invL;

                   xvector(dirx, diry, dirz, node1->burgX[i], node1->burgY[i],
                           node1->burgZ[i], &n1x, &n1y, &n1z);

                   xvector(dirx, diry, dirz, node1->vX, node1->vY, node1->vZ,
                           &n2x, &n2y, &n2z);

                   n1mag2 = n1x*n1x + n1y*n1y + n1z*n1z;
                   n2mag2 = n2x*n2x + n2y*n2y + n2z*n2z;

                   if (n1mag2 > eps) {
/*
 *                 Preference for plane defined by l cross v
 */
                       norm = sqrt(n1mag2);
                       invnorm = 1.0 / norm;
                       plane[0] = n1x * invnorm;
                       plane[1] = n1y * invnorm;
                       plane[2] = n1z * invnorm;
                   } else if (n2mag2 > eps) {
/*
 *                 Preference for plane defined by l cross b
 */
                       norm = sqrt(n2mag2);
                       invnorm = 1.0 / norm;
                       plane[0] = n2x * invnorm;
                       plane[1] = n2y * invnorm;
                       plane[2] = n2z * invnorm;
                   } else {
/*
 *                 Use segment's defined glide plane if glide planes enforced
 */
                       plane[0] = node1->nx[i];
                       plane[1] = node1->ny[i];
                       plane[2] = node1->nz[i];
                   }


                   switch (Nsize) {
                   case 0:
                       conditionsmet = 1;
                       break;

                   case 1:
                       tmp = Nmat[0][0]*plane[0] +
                             Nmat[0][1]*plane[1] +
                             Nmat[0][2]*plane[2];
                       conditionsmet = (tmp*tmp < npc2);
                       break;

                   default:
                      Nmat[2][0] = plane[0];
                      Nmat[2][1] = plane[1];
                      Nmat[2][2] = plane[2];

                      detN = Matrix33_Det(Nmat);
                      conditionsmet = (detN*detN > onemnpc4);
                   }

                   if (conditionsmet) {
                       Nmat[Nsize][0] = plane[0];
                       Nmat[Nsize][1] = plane[1];
                       Nmat[Nsize][2] = plane[2];
                       vector[Nsize] = DotProduct(plane, p1);
                       Nsize++;
                   }
               }
           }

           for (i = 0; i < node2->numNbrs; i++) {

               if (Nsize < 3) {

                   nbrNode = GetNeighborNode(home, node2, i);

                   if (nbrNode == (Node_t *)NULL) {
                       printf("WARNING: Neighbor not found at %s line %d\n",
                              __FILE__, __LINE__);
                       continue;
                   }

                   dx = p2[0] - nbrNode->x;
                   dy = p2[1] - nbrNode->y;
                   dz = p2[2] - nbrNode->z;

                   ZImage(param, &dx, &dy, &dz);

                   L = sqrt(dx*dx + dy*dy + dz*dz);
                   invL = 1.0 / L;
                   dirx = dx * invL;
                   diry = dy * invL;
                   dirz = dz * invL;

                   xvector(dirx, diry, dirz, node2->burgX[i], node2->burgY[i],
                           node2->burgZ[i], &n1x, &n1y, &n1z);

                   xvector(dirx, diry, dirz, node2->vX, node2->vY, node2->vZ,
                           &n2x, &n2y, &n2z);

                   n1mag2 = n1x*n1x + n1y*n1y + n1z*n1z;
                   n2mag2 = n2x*n2x + n2y*n2y + n2z*n2z;

                   if (n1mag2 > eps) {
                       norm = sqrt(n1mag2);
                       invnorm = 1.0 / norm;
                       plane[0] = n1x * invnorm;
                       plane[1] = n1y * invnorm;
                       plane[2] = n1z * invnorm;
                   } else if (n2mag2 > eps) {
                       norm = sqrt(n2mag2);
                       invnorm = 1.0 / norm;
                       plane[0] = n2x * invnorm;
                       plane[1] = n2y * invnorm;
                       plane[2] = n2z * invnorm;
                   } else {
                       plane[0] = node2->nx[i];
                       plane[1] = node2->ny[i];
                       plane[2] = node2->nz[i];
                   }
 

                   switch (Nsize) {
                   case 1:
                       tmp = Nmat[0][0]*plane[0] +
                             Nmat[0][1]*plane[1] +
                             Nmat[0][2]*plane[2];
                       conditionsmet = (tmp*tmp < npc2);
                       break;
                   default:
                      Nmat[2][0] = plane[0];
                      Nmat[2][1] = plane[1];
                      Nmat[2][2] = plane[2];
                      detN = Matrix33_Det(Nmat);
                      conditionsmet = (detN*detN > onemnpc4);
                   }
 
                   if (conditionsmet) {
                       Nmat[Nsize][0] = plane[0];
                       Nmat[Nsize][1] = plane[1];
                       Nmat[Nsize][2] = plane[2];
                       vector[Nsize] = DotProduct(plane, p2);
                       Nsize++;
                   }
               }
           }

/*
 *         Upper left 3X3 of Matrix is identity matrix.
 *         Matrix rows 3 thru 3+(Nsize-1) colums 0 thru 2 are Nmat.
 *         Matrix columns 3 thru 3+(Nsize-1) rows 0 thru 2 are transpose of Nmat.
 *         All remaining elements are zeroed.
 */
           for (i = 0; i < 6; i++) {
               for (j = 0; j < 6; j++) {
                   Matrix[i][j] = 0.0;
               }
           }

           Matrix[0][0] = 1.0;
           Matrix[1][1] = 1.0;
           Matrix[2][2] = 1.0;

           for (i = 0; i < Nsize; i++) {
               for (j = 0; j < 3; j++) {
                   Matrix[3+i][j] = Nmat[i][j];
                   Matrix[j][3+i] = Nmat[i][j];
               }
           }
           V[0] = 0.5 * (p1[0] + p2[0]);
           V[1] = 0.5 * (p1[1] + p2[1]);
           V[2] = 0.5 * (p1[2] + p2[2]);
           V[3] = vector[0];
           V[4] = vector[1];
           V[5] = vector[2];

           Nsize += 3;

           if (!MatrixInvert((real8 *)Matrix, (real8 *)invMatrix, Nsize, 6)) {
               printf("FindCollisionPoint(): Unable to invert %dx%d matrix!\n",
                      Nsize, Nsize);
#if 0
               printf("FindCollisionPoint(): node1 = (%d,%d), node2 = (%d,%d)\n",
                      node1->myTag.domainID, node1->myTag.index,
                      node2->myTag.domainID, node2->myTag.index);
               for (m = 0; m < Nsize; m++) {
                   for (n = 0; n < Nsize; n++) {
                       printf("    Matrix[%d][%d] = %e\n", m, n, Matrix[m][n]);
                   }
               }
               PrintNode(node1);
               PrintNode(node2);
#endif
               fflush(NULL);
               return(0);
           }

           MatrixMult((real8 *)invMatrix, Nsize, Nsize, 6,
                      (real8 *)V, 1, 1,
                      (real8 *)result, 1);

           *x = result[0];
           *y = result[1];
           *z = result[2];

       }  /* enforceGlidePlanes != 0 */

       return(1);
}


/*---------------------------------------------------------------------------
 *
 *      Function:       GetMinDist
 *      Description:    Find the minimum distance between the two segments
 *                      p1->p2 andp3->p4.
 *
 *      Arguments:
 *          p1x, p1y, p1z  Coordinates of the first endpoint of segment 1
 *          v1x, v1y, v1z  Velocity of the node at point p1
 *          p2x, p2y, p2z  Coordinates of the second endpoint of segment 1
 *          v2x, v2y, v2z  Velocity of the node at point p2
 *          p3x, p3y, p3z  Coordinates of the first endpoint of segment 2
 *          v3x, v3y, v3z  Velocity of the node at point p3
 *          p4x, p4y, p4z  Coordinates of the second endpoint of segment 2
 *          v4x, v4y, v4z  Velocity of the node at point p4
 *          dist2          pointer to location in which to return the square
 *                         of the minimum distance between the two points
 *          ddist2dt       pointter to the location in which to return the
 *                         time rate of change of the distance between
 *                         L1 and L2
 *          L1             pointer to location at which to return the 
 *                         normalized position on seg 1 closest to segment 2
 *          L2             pointer to location at which to return the 
 *                         normalized position on seg 2 closest to segment 1
 *
 *-------------------------------------------------------------------------*/
void GetMinDist(real8 p1x, real8 p1y, real8 p1z,
                real8 v1x, real8 v1y, real8 v1z,
                real8 p2x, real8 p2y, real8 p2z,
                real8 v2x, real8 v2y, real8 v2z,
                real8 p3x, real8 p3y, real8 p3z,
                real8 v3x, real8 v3y, real8 v3z,
                real8 p4x, real8 p4y, real8 p4z,
                real8 v4x, real8 v4y, real8 v4z,
                real8 *dist2, real8 *ddist2dt, real8 *L1, real8 *L2)
{
        real8  A, B, C, D, E, F, G;
        real8  L1a=0, L1b=0;
        real8  L2a=0, L2b=0;
        real8  dist2a, dist2b;
        real8  eps = 1.0e-12;
        real8  ddistxdt, ddistydt, ddistzdt;
        real8  seg1Lx, seg1Ly, seg1Lz;
        real8  seg2Lx, seg2Ly, seg2Lz;
        real8  seg1vx, seg1vy, seg1vz;
        real8  seg2vx, seg2vy, seg2vz;

        seg1Lx = p2x - p1x;
        seg1Ly = p2y - p1y;
        seg1Lz = p2z - p1z;

        seg2Lx = p4x - p3x;
        seg2Ly = p4y - p3y;
        seg2Lz = p4z - p3z;

        seg1vx = v2x - v1x;
        seg1vy = v2y - v1y;
        seg1vz = v2z - v1z;

        seg2vx = v4x - v3x;
        seg2vy = v4y - v3y;
        seg2vz = v4z - v3z;

        A = seg1Lx*seg1Lx + seg1Ly*seg1Ly + seg1Lz*seg1Lz;
        B = ((seg1Lx*(p1x-p3x)) +
             (seg1Ly*(p1y-p3y)) +
             (seg1Lz*(p1z-p3z))) * 2.0;
        C = (seg1Lx*seg2Lx + seg1Ly*seg2Ly + seg1Lz*seg2Lz) * 2.0;
        D = ((seg2Lx*(p3x-p1x)) +
             (seg2Ly*(p3y-p1y)) +
             (seg2Lz*(p3z-p1z))) * 2.0;
        E = seg2Lx*seg2Lx + seg2Ly*seg2Ly + seg2Lz*seg2Lz;
        F = (p1x-p3x) * (p1x-p3x) +
            (p1y-p3y) * (p1y-p3y) +
            (p1z-p3z) * (p1z-p3z);
        G = C*C - 4*A*E;


        if (A < eps) {
/*
 *          Segment 1 is just a point...
 */
            *L1 = 0.0;
            if (E < eps) {
                *L2 = 0.0;
            } else {
                *L2 = -0.5 * D / E;
                *L2 = MIN(MAX(*L2, 0.0), 1.0);
            }

            *dist2 = F + (D + E * *L2) * *L2;

        } else if (E < eps) {
/*
 *          Segment 2 is just a point...
 */
            *L2 = 0.0;
            *L1 = -0.5 * B / A;
            *L1 = MIN(MAX(*L1, 0.0), 1.0);

            *dist2 = F + (B + A * *L1) * *L1;

        } else if (fabs(G) < eps) {
/*
 *          Segments are parallel.  Check ends of seg1
 */

            L2a = -0.5 * D / E;
            L2a = MIN(MAX(L2a, 0.0), 1.0);
            dist2a = F + (D + E * L2a) * L2a;

            L2b = 0.5 * (C - D) / E;
            L2b = MIN(MAX(L2b, 0.0), 1.0);
            dist2b = A + B + F + (D - C + E * L2b) * L2b;
/*
 *          Use the closest end of seg 1
 */
            if (dist2a < dist2b) {
                *L1 = 0.0;
                *L2 = L2a;
                *dist2 = dist2a;
            } else {
                *L1 = 1.0;
                *L2 = L2b;
                *dist2 = dist2b;
            }

        } else {
/*
 *          Non parallel segments
 */
            *L2 = (2.0 * A * D + B * C) / G;
            *L1 = 0.5 * (C * *L2 - B) / A;
/*
 *          Now make sure L1 and L2 are between 0.0 and 1.0
 */
            if ((*L1 >= 0.0) && (*L1 <= 1.0) &&
                (*L2 >= 0.0) && (*L2 <= 1.0)) {
                *dist2 = F +(B + A * *L1) * *L1 + (D + E * *L2) *
                        *L2 - C * *L1 * *L2;
            } else {
/*
 *              Special case where one or more of the endpoints is closest.
 *              First check L1=0 edge
 */
                if (*L1 < 0.0) {
                    L1a = 0.0;
                    L2a = -0.5 * D / E;
                    L2a = MIN(MAX(L2a, 0.0), 1.0);
                    dist2a = F + (D + E * L2a) * L2a;
/*
 *              Check L1=1 edge
 */
                } else if (*L1 > 1.0) {
                    L1a = 1.0;
                    L2a = 0.5 * (C - D) / E;
                    L2a = MIN(MAX(L2a, 0.0), 1.0);
                    dist2a = A + B + F + (D - C + E * L2a) * L2a;
                } else {
                    dist2a = -1.0;
                }

/*
 *              Check L2=0 edge
 */
                if (*L2 < 0.0) {
                    L2b = 0.0;
                    L1b = -0.5 * B / A;
                    L1b = MIN(MAX(L1b, 0.0), 1.0);
                    dist2b = F + (B + A * L1b) * L1b;
/*
 *              Check L2=1 edge
 */
                } else if (*L2 > 1.0) {
                    L2b = 1.0;
                    L1b = 0.5 * (C - B) / A;
                    L1b = MIN(MAX(L1b, 0.0), 1.0);
                    dist2b = D + E + F + (B - C + A * L1b) * L1b;
                } else {
                    dist2b = -1.0;
                }

/*
 *              Find the shortest distance
 */
                if ((dist2b < 0.0) ||
                    ((dist2a <= dist2b) && dist2a >= 0.0)) {
                    *L1 = L1a;
                    *L2 = L2a;
                    *dist2 = dist2a;
                } else {
                    *L1 = L1b;
                    *L2 = L2b;
                    *dist2 = dist2b;
                }
            }

        }  /* end code for non-parallel segments */

/*
 *      Lastly, the time rate of change between the points at L1 and L2.
 */
        ddistxdt = (v1x + (seg1vx * *L1) - v3x - (seg2vx * *L2)) *
                   (p1x + (seg1Lx * *L1) - p3x - (seg2Lx * *L2));
        ddistydt = (v1y + (seg1vy * *L1) - v3y - (seg2vy * *L2)) *
                   (p1y + (seg1Ly * *L1) - p3y - (seg2Ly * *L2));
        ddistzdt = (v1z + (seg1vz * *L1) - v3z - (seg2vz * *L2)) *
                   (p1z + (seg1Lz * *L1) - p3z - (seg2Lz * *L2));

        *ddist2dt = 2.0 * (ddistxdt + ddistydt + ddistzdt);

        return;
}

void HandleCollisions(Home_t *home)
{
    switch(home->param->collisionMethod)  
    {
        case 4:  { PredictiveCollisionsLG (home); break; }
        case 3:  { RetroactiveCollisions  (home); break; }
        case 2:  { PredictiveCollisions   (home); break; }
        case 1: 
        default: { ProximityCollisions    (home); break; }
    }

    return;
}
