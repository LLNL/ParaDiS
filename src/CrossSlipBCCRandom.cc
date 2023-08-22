//
//
// Module:         CrossSlipBCCRandom.c
//
//      Author:         Sylvie Aubry
//
//      Description:    Version of crossSlip for BCC crystals that uses power dissipation 
//                      and area growth criteria and applies on 6 glide planes + randomness
//                      in choice of cross-slip plane
//
//
//      Includes public functions:
//                      CrossSlipBCCRandom()
//
//      Date : 08/02/2017



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "CrossSlip.h"

#ifdef DEBUG_CROSSSLIP_EVENTS
static int dbgDom;
#endif


static real8 IsAreaGrowth(real8 P1[3], real8 P2[3], real8 Q[3],
                          real8 v1[3], real8 v2[3], real8 v[3])
{
   real8 delta = 1.0e-16, s, dsdt;
   real8 r1, dr1dt, r2, dr2dt, r3, dr3dt;
   real8 vecP1P2[3], vecQP1[3], vecQP2[3];
   real8 dvecP1P2dt[3], dvecQP1dt[3], dvecQP2dt[3];
   real8 dA2dt;

   vecQP1[X] = P1[X] - Q[X]; 
   vecQP1[Y] = P1[Y] - Q[Y];
   vecQP1[Z] = P1[Z] - Q[Z];
   
   vecQP2[X] = P2[X] - Q[X]; 
   vecQP2[Y] = P2[Y] - Q[Y];
   vecQP2[Z] = P2[Z] - Q[Z];
   
   vecP1P2[X] = P2[X] - P1[X]; 
   vecP1P2[Y] = P2[Y] - P1[Y];
   vecP1P2[Z] = P2[Z] - P1[Z];
   
   dvecQP1dt[X] = v1[X] - v[X];
   dvecQP1dt[Y] = v1[Y] - v[Y];
   dvecQP1dt[Z] = v1[Z] - v[Z];
   
   dvecQP2dt[X] = v2[X] - v[X];
   dvecQP2dt[Y] = v2[Y] - v[Y];
   dvecQP2dt[Z] = v2[Z] - v[Z];
   
   dvecP1P2dt[X] = v2[X] - v1[X];
   dvecP1P2dt[Y] = v2[Y] - v1[Y];
   dvecP1P2dt[Z] = v2[Z] - v1[Z];
   
   r1 = sqrt(DotProduct(vecQP1, vecQP1));
   r2 = sqrt(DotProduct(vecQP2, vecQP2));
   r3 = sqrt(DotProduct(vecP1P2, vecP1P2));
   
   s = 0.5 * (r1 + r2 + r3);
   
   dr1dt = (vecQP1[X]*dvecQP1dt[X] + vecQP1[Y]*dvecQP1dt[Y] +
            vecQP1[Z]*dvecQP1dt[Z] ) / (r1 + delta);
   
   dr2dt = (vecQP2[X]*dvecQP2dt[X] + vecQP2[Y]*dvecQP2dt[Y] +
            vecQP2[Z]*dvecQP2dt[Z] ) / (r2 + delta);
   
   dr3dt = (vecP1P2[X] * dvecP1P2dt[X] + vecP1P2[Y] *
            dvecP1P2dt[Y] + vecP1P2[Z] * dvecP1P2dt[Z]) /
           (r3 + delta);
   
   dsdt = 0.5 * (dr1dt + dr2dt + dr3dt);
   
   dA2dt = dsdt * (s-r1) * (s-r2) * (s-r3);
   dA2dt += s * (dsdt-dr1dt) * (s-r2) * (s-r3);
   dA2dt += s * (s-r1) * (dsdt-dr2dt) * (s-r3);
   dA2dt += s * (s-r1) * (s-r2) * (dsdt-dr3dt);
   return  dA2dt; //(dA2dt > 0) ? 1:0;
}

static void DefinePos(Param_t *param, Node_t *node, Node_t* nbr1, Node_t *nbr2,
                      real8 nodep[3], real8 nbr1p[3], real8 nbr2p[3],
                      real8 n1mn2[3], real8 nmn1[3], real8 nmn2[3])
{
   //  Get the positions of all three nodes and convert the neighbor
   //  node positions to the PBC coordinates relative to the main
   //  node.  These adjusted positions will be used during the
   //  cross slip algorithm, and if updated during the process, will
   //  be copied back into the respective nodal data structures and
   //  shifted back into the primary image space.
   nodep[X] = node->x; nodep[Y] = node->y; nodep[Z] = node->z;
   nbr1p[X] = nbr1->x; nbr1p[Y] = nbr1->y; nbr1p[Z] = nbr1->z;
   nbr2p[X] = nbr2->x; nbr2p[Y] = nbr2->y; nbr2p[Z] = nbr2->z;
   
   PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
               &nbr1p[X], &nbr1p[Y], &nbr1p[Z]);
   PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
               &nbr2p[X], &nbr2p[Y], &nbr2p[Z]);
   
   n1mn2[X] = nbr1p[X] - nbr2p[X];
   n1mn2[Y] = nbr1p[Y] - nbr2p[Y];
   n1mn2[Z] = nbr1p[Z] - nbr2p[Z];
   
   nmn1[X] = nodep[X] - nbr1p[X];
   nmn1[Y] = nodep[Y] - nbr1p[Y];
   nmn1[Z] = nodep[Z] - nbr1p[Z];
   
   nmn2[X] = nodep[X] - nbr2p[X];
   nmn2[Y] = nodep[Y] - nbr2p[Y];
   nmn2[Z] = nodep[Z] - nbr2p[Z];
}

//
// Takes in the power on 6 glide planes, return the index of the 
// chosen plane index chosen with some randomness
//
static int PowerChoose(real8 p[6], real8 units)
{
   static int   seed = 8917346;
   int iplane = 0, i;
   int pindexmin=0;
   real8 fp[6], sumfp, gamma[6], proba[6], alpha;
   real8 pmin=0.0;

#if 0
   {
      // For deterministic cross-slip
      int pindexmax=0;
      real8 pmax;
      FindAbsMax(p,6,&pmax,&pindexmax);      
      return pindexmax;
   }
#endif

  FindAbsMin(p,6,&pmin,&pindexmin);

   sumfp = 0.0;
   for (i = 0; i < 6; i++)  
   {
      fp[i] = pow(fabs(p[i]/pmin),8);
      sumfp += fp[i];
   }
         
   if (isnan(sumfp)) 
     {
        int pindexmax=0;
        real8 pmax=0.0;
        FindAbsMax(p,6,&pmax,&pindexmax);      
        return pindexmax;
     }

   for (i = 0; i < 6; i++)  
   {
      proba[i] = fp[i]/sumfp;
   }

   gamma[0] = proba[0];
   real8 sumgamma = gamma[0];
   for (i = 1; i < 6; i++)  
   {
      gamma[i] = gamma[i-1] + proba[i];
      sumgamma += gamma[i];
   }


   alpha = randm(&seed);

   if (alpha <= gamma[0]) iplane = 0;

   for (i = 1; i < 6; i++)  
      if (alpha > gamma[i-1] && alpha <= gamma[i])  iplane = i;

   return iplane;
}

//
// Function:       CrossSlipBCCRandom
//
// Description:    Examines all nodes local to the domain, determines
//                 if the node should cross slip, if the node is
//                 permitted to cross slip, adjusts node positions
//                 and segment glide planes as necessary.
void CrossSlipBCCRandom(Home_t *home)
{
        int    i;
        int    numNodes;
        int    opClass, thisDom;
        int    eventListSize, zipEventCnt, slipEventIndex;
        real8  eps, thetacrit, sthetacrit, s2thetacrit, areamin;
        Param_t *param;
        SlipEvent_t *eventList = (SlipEvent_t *)NULL;
        NodeState_t origNodeState, origNbr1State, origNbr2State;

        // Allocate a block of structures to hold info about cross-slip
        // events.  We allocate enough structures to cross-slip every local
        // node, which is overkill, but it means we won't run out of entries
        // and have to reallocate the block, plus minimizes the number of
        // operations we'll have to do within the thread-critical sections
        // of code.
        eventListSize = home->newNodeKeyPtr;
        slipEventIndex = eventListSize;
        zipEventCnt = 0;
        
        if (eventListSize > 0) 
        {
           int blockSize;
           blockSize = eventListSize * sizeof(SlipEvent_t);
           eventList = (SlipEvent_t *)malloc(blockSize);
        }
        
#ifdef DEBUG_CROSSSLIP_EVENTS
#ifdef DEBUG_CROSSSLIP_DOMAIN
        dbgDom = DEBUG_CROSSSLIP_DOMAIN;
#else
        dbgDom = -1;
#endif
#endif
        
        param = home->param;
        thisDom = home->myDomain;
        numNodes = home->newNodeKeyPtr;
        
        eps = 1.0e-06;
        thetacrit = 0.5 / 180.0 * M_PI;
        sthetacrit = sin(thetacrit);
        s2thetacrit = sthetacrit * sthetacrit;
        areamin = param->remeshAreaMin;
        

        // For purposes of determining segment ownership, we need to treat
        // cross slip operations the same as we do during collision handling.
        opClass = OPCLASS_COLLISION;
        
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
           int nIndex;
           int threadID, threadIterStart, threadIterEnd;
                   
           GetThreadIterationIndices(numNodes, &threadID, &threadIterStart,
                                     &threadIterEnd);
           
           // Each thread loops through a portion of the nodes evaluating
           // any 2-nodes as possible cross-slip candidates.
           for (nIndex = threadIterStart; nIndex < threadIterEnd; nIndex++) {
              int    seg1_is_screw, seg2_is_screw, bothseg_are_screw;
              real8  burgSize, tmp;
              real8  test1, test2, test3;
              real8  testmax1, testmax2, testmax3;
              real8  burgLab[3], burgCrystal[3];
              real8  nodep[3], nbr1p[3], nbr2p[3];
              real8  n1mn2[3], nmn1[3], nmn2[3];
              Node_t *node, *nbr1, *nbr2;
              real8  vAligned;
              
              if ((node = home->nodeKeys[nIndex]) == (Node_t *)NULL) {
                 continue;
              }
              
              // Discretization nodes only
              if (node->numNbrs != 2) {
                 continue;
              }
              
              if (!HAS_NO_CONSTRAINTS(node->constraint)) {
                 continue;
              }
              
              nbr1 = GetNeighborNode(home, node, 0);
              nbr2 = GetNeighborNode(home, node, 1);
              
              if ((nbr1 == (Node_t *)NULL) ||
                  (nbr2 == (Node_t *)NULL)) {
                 printf("WARNING: Neighbor not found at %s line %d\n",
                        __FILE__, __LINE__);
                 continue;
              }
              
              
              // We must not change the glide plane for any segment
              // that is not 'owned' by the node in the current domain.
              // Otherwise, neighboring domains *could* change the
              // glide plane for the same segment... this would be bad.
              if ((!DomainOwnsSeg(home,opClass,thisDom,&node->nbrTag[0])) ||
                  (!DomainOwnsSeg(home,opClass,thisDom,&node->nbrTag[1]))) {
                 continue;
              }

              // Recover the Burgers vector
              burgLab[X] = node->burgX[0];
              burgLab[Y] = node->burgY[0];
              burgLab[Z] = node->burgZ[0];
              burgSize = sqrt(DotProduct(burgLab,burgLab));
              vAligned = eps*eps*burgSize*burgSize;
              
              NormalizeVec(burgLab);

              VECTOR_COPY(burgCrystal, burgLab);
              NormalizeVec(burgCrystal);

             
              // If needed, rotate a copy of the burgers vector from the lab
              // frame to the crystal frame.  Otherwise, lab and crystal frames
              // are identical.
              if (param->useLabFrame) 
              {
                 Matrix33Vector3Multiply(param->rotMatrixInverse, burgLab,
                                         burgCrystal);
              }
              

              //  Check that we are dealing with a <111> type Burgers vector
              if (fabs(burgCrystal[X]*burgCrystal[Y]*burgCrystal[Z]) < eps)  continue;

              //  Get the reference Burgers vertor
              int m, bIndex=-1, numplanes = 6;
              int row;
              int  *burgFirstPlaneIndex;
              real8 burgRef[3];
              real8 (*burgList)[3], (*planeList)[3];
              burgList = home->burgData.burgList;
              planeList = home->burgData.planeList;
              burgFirstPlaneIndex = home->burgData.burgFirstPlaneIndex;
              
              bIndex = GetBurgIndexBCC(burgLab, BCC_NUM_TOT_BURG, burgList);
              
              VECTOR_COPY(burgRef,burgList[bIndex]);

               // Get the positions of all three nodes and convert the neighbor
              // node positions to the PBC coordinates relative to the main
              // node.  These adjusted positions will be used during the
              // cross slip algorithm, and if updated during the process, will
              // be copied back into the respective nodal data structures and
              // shifted back into the primary image space.
              DefinePos(param,node,nbr1,nbr2,nodep,nbr1p,nbr2p,n1mn2,nmn1,nmn2);

              // Calculate some cross-slip test conditions
              test1 = DotProduct(n1mn2, burgLab);
              test2 = DotProduct(nmn1,  burgLab);
              test3 = DotProduct(nmn2,  burgLab);
              
              test1 = test1 * test1;
              test2 = test2 * test2;
              test3 = test3 * test3;
              
              testmax1 = DotProduct(n1mn2, n1mn2);
              testmax2 = DotProduct(nmn1,  nmn1);
              testmax3 = DotProduct(nmn2,  nmn2);
              
              // Set up the tests to see if this dislocation is close enough to 
              // screw to be considered for cross slip.  For a segment to be
              // close to screw it must be within 2*thetacrit defined above
              seg1_is_screw = ((testmax2 - test2) < (testmax2 * s2thetacrit));
              seg2_is_screw = ((testmax3 - test3) < (testmax3 * s2thetacrit));

              bothseg_are_screw =
                 (((testmax2-test2) < (4.0 * testmax2 * s2thetacrit)) &&
                  ((testmax3-test3) < (4.0 * testmax3 * s2thetacrit)) &&
                  ((testmax1-test1) < (testmax1 * s2thetacrit)));
              
              // Neither seg1 nor seg2 is screw. Do nothing.
              if (!bothseg_are_screw && !seg1_is_screw && !seg2_is_screw) continue;
              
              // Initialize power dissipation for this initial configuration
              // composed of three perfectly aligned nodes along the screw 
              // configuration
              real8 powerMax =  (node->fX * node->vX) + 
                                (node->fY * node->vY) + 
                                (node->fZ * node->vZ);
                   

             // Either one or both segments are screw, examine cross-slip
              int    j, iglide, DoAddEvent=0, Keepglide=0;
              int    pinned1, pinned2;
              int    glide1, glide2;
              int    thisEventIndex;                
              MobArgs_t  mobArgs;
              real8  vNoise = param->splitMultiNodeAlpha * param->rTol / param->deltaTT;
              real8  coeff=areamin * 2.0 * (1.0 + eps);
              real8  Keepfdotglide;
              real8  n1mn2dotb, nmn1dotb, nmn2dotb, fdotglide;
              real8  powerTest[6], powerTestZip,units;
              real8  segplane1[3], segplane2[3];
              real8  tmp3[6], tmp3B[6];
              real8  glideDirCrystal[6][3], glideDirLab[6][3];
              real8  planeDirCrystal[6][3], planeDirLab[6][3];
              real8  KeepNodep[3], KeepNbrp[3], KeepNodepFor[3]; 
              SlipEvent_t *thisEvent;

              units = 1e2*fabs((node->fX * node->vX) + (node->fY * node->vY) + (node->fZ * node->vZ));
              units = 1/units;
           

              // Get the 6 plane and glide reference directions
              int start = burgFirstPlaneIndex[bIndex];
              int end   = burgFirstPlaneIndex[bIndex] + numplanes;
              for (m=0, j=start; j<end; j++, m++) 
              {
                 VECTOR_COPY(planeDirCrystal[m],planeList[j]);
                 NormalizedCrossVector(planeDirCrystal[m],burgRef,glideDirCrystal[m]);
              }

#if 0     // Check that the right planes and glides directions were set           
              for (row = 0; row < 6; row++) 
              {
                 Print3("glide",glideDirCrystal[row]);
                 Print3("plane",planeDirCrystal[row]);
                 printf("\n");
              }
              exit(0);
#endif

              // Planes associated with current node configuration
              segplane1[X] = node->nx[0];
              segplane1[Y] = node->ny[0];
              segplane1[Z] = node->nz[0];
              
              segplane2[X] = node->nx[1];
              segplane2[Y] = node->ny[1];
              segplane2[Z] = node->nz[1];
              
              // Need copies of the glideDir matrix and the force vector
              // in both the lab frame and the crystal frame for later use.
              // Also need to rotate the segment planes into the crystal
              // frame.
              if (param->useLabFrame)
              {
                 real8 tmpPlane[3];
                 // Rotate initial glide dir array from crystal to lab frame
                 for (row = 0; row < 6; row++) 
                 {
                    Matrix33Vector3Multiply(param->rotMatrix,
                                            glideDirCrystal[row],
                                            glideDirLab[row]);
                    
                    Matrix33Vector3Multiply(param->rotMatrix,
                                            planeDirCrystal[row],
                                            planeDirLab[row]);
                    
                 }
                   
                 // Rotate segment planes to crystal frame                 
                 Matrix33Vector3Multiply(param->rotMatrixInverse,
                                         segplane1, tmpPlane);
                 VECTOR_COPY(segplane1, tmpPlane);
                 
                 Matrix33Vector3Multiply(param->rotMatrixInverse,
                                         segplane2, tmpPlane);
                 VECTOR_COPY(segplane2, tmpPlane);                 
              } 
              else 
              {
                 // Lab and crystal frames are identical, so just copy
                 // the vectors
                 int row;
                 for (row = 0; row < 6; row++) 
                 {
                    VECTOR_COPY(glideDirLab[row], glideDirCrystal[row]);
                    VECTOR_COPY(planeDirLab[row], planeDirCrystal[row]);
                 }
              }


              // Find the glide direction associated with the current glide plane
              // set in the node structure. Used later on in zipper case              
              Matrix63_Vmul(tmp3 , glideDirCrystal, segplane1 );   // tmp3  =  glideDirCrystal * segplane1
              Matrix63_Vmul(tmp3B, glideDirCrystal, segplane2 );   // tmp3B =  glideDirCrystal * segplane2
              glide1 = 0;
              glide2 = 0;              
              for (j = 1; j < 6; j++) 
              {
                 glide1=(fabs(tmp3[j]) <fabs(tmp3[glide1]))  ? j:glide1;
                 glide2=(fabs(tmp3B[j])<fabs(tmp3B[glide2])) ? j:glide2;
                 
              }
              
              // Check if node's neighbors are pinned
              pinned1 = NodePinned(home, nbr1, glide1, glideDirLab, 6);                    
              pinned2 = NodePinned(home, nbr2, glide2, glideDirLab, 6);
             

              // Nodes position are going to change. Keep current positions.
              PreserveState(home, node, &origNodeState);
              PreserveState(home, nbr1, &origNbr1State);
              PreserveState(home, nbr2, &origNbr2State);
                            
              // Cross-slip of a node closely aligned to two neighbors
              if (bothseg_are_screw ) 
              {
                 if (pinned1) 
                   { 
                      // Neighbor 1 cannot be moved
                      DoAddEvent = 0;

                      // Align both segments exactly along screw by moving nbr2 and node.
                      n1mn2dotb = DotProduct(n1mn2, burgLab);
                      if  (!pinned2) 
                      {
                         nbr2p[X] = nbr1p[X]-(n1mn2dotb * burgLab[X]);
                         nbr2p[Y] = nbr1p[Y]-(n1mn2dotb * burgLab[Y]);
                         nbr2p[Z] = nbr1p[Z]-(n1mn2dotb * burgLab[Z]);
                      }
                      
                      nmn1dotb = DotProduct(nmn1, burgLab);                           
                      nodep[X] = nbr1p[X] + (nmn1dotb * burgLab[X]);
                      nodep[Y] = nbr1p[Y] + (nmn1dotb * burgLab[Y]);
                      nodep[Z] = nbr1p[Z] + (nmn1dotb * burgLab[Z]);

                      FoldBox(param, &nodep[X], &nodep[Y], &nodep[Z]);
                      RepositionNode(home, nodep, &node->myTag, 0);
                      
                      FoldBox(param, &nbr2p[X], &nbr2p[Y], &nbr2p[Z]);
                      RepositionNode(home, nbr2p, &nbr2->myTag, 0);

                      // Evaluate forces of node
                      SetOneNodeForce(home, node);  

                      // Force is the same for any of the 6 planes
                      KeepNodepFor[X] = node->fX;
                      KeepNodepFor[Y] = node->fY;
                      KeepNodepFor[Z] = node->fZ;
 
                      // Find the highest power configuration among the 6 glide planes
                      int mobError[6];
                      for (iglide = 0; iglide < 6; iglide++) 
                      {
                         ResetGlidePlane(home, planeDirLab[iglide], &node->myTag,
                                         &nbr1->myTag, 0);
                         ResetGlidePlane(home, planeDirLab[iglide], &node->myTag,
                                         &nbr2->myTag, 0);

                         mobError[iglide]  = EvaluateMobility(home, node, &mobArgs);
  
                          powerTest[iglide] = (node->fX * node->vX) +   (node->fY * node->vY) + (node->fZ * node->vZ) + 
                               vNoise * sqrt((node->fX * node->fX) +  (node->fY * node->fY) +  (node->fZ * node->fZ) );
                      }
                      // Use highest power dissipation criterion + randomness to choose the glide plane
                      // to cross-slip on
                      Keepglide = PowerChoose(powerTest,units);
                      // if (powerTest[Keepglide] > powerMax && mobError[Keepglide] == 0) 
                      if (mobError[Keepglide] == 0) 
                      {
                         DoAddEvent = 1;

                         // Node position and its neighbor position does not change with glide plane
                         VECTOR_COPY(KeepNodep,nodep);
                         VECTOR_COPY(KeepNbrp,nbr2p);
                      }
                     
                        
                      if (DoAddEvent)
                      {
                         // Cross-slip can happen
                         // Check that for this configuration, the area will grow.
                         // Move the node and evaluate the area growth criterion
 
                         n1mn2dotb = DotProduct(n1mn2, burgLab);
                         fdotglide = DotProduct(KeepNodepFor,glideDirLab[Keepglide]);                        
                         tmp = Sign(fdotglide) * coeff / fabs(n1mn2dotb) ;
                         
                         nodep[X] += tmp * glideDirLab[Keepglide][X];
                         nodep[Y] += tmp * glideDirLab[Keepglide][Y];
                         nodep[Z] += tmp * glideDirLab[Keepglide][Z];
                         
                         FoldBox(param, &nodep[X], &nodep[Y], &nodep[Z]);
                         RepositionNode(home, nodep, &node->myTag, 0);

                         SetOneNodeForce(home, node);
                         SetOneNodeForce(home, nbr1);
                         SetOneNodeForce(home, nbr2);

                         ResetGlidePlane(home, planeDirLab[Keepglide], &node->myTag,
                                         &nbr1->myTag, 0);
                         ResetGlidePlane(home, planeDirLab[Keepglide], &node->myTag,
                                         &nbr2->myTag, 0);

                         int mobError0;
                         mobError0  = EvaluateMobility(home, node, &mobArgs);
                         mobError0 |= EvaluateMobility(home, nbr1, &mobArgs);
                         mobError0 |= EvaluateMobility(home, nbr2, &mobArgs);           
                         
                         // Area growth criterion's test
                         real8 v1[3],v[3],v2[3];                         
                         v1[X] = nbr1->vX; v1[Y] = nbr1->vY; v1[Z]= nbr1->vZ;
                         v [X] = node->vX;  v[Y] = node->vY;  v[Z]= node->vZ;
                         v2[X] = nbr2->vX; v2[Y] = nbr2->vY; v2[Z]= nbr2->vZ;
                         real8 dAreadt = IsAreaGrowth(nbr1p,nbr2p,nodep,v1,v2,v);
                         if (dAreadt > 0.0 && mobError0 == 0)
                         {
                            DoAddEvent = 1;                           
                            Keepfdotglide = fdotglide;                          
                            VECTOR_COPY(KeepNodep,nodep);
                         }
                         else
                            DoAddEvent = 0;
                      }

                      // Restore initial node positions for next trial
                      RestoreState(home, node, &origNodeState);
                      RestoreState(home, nbr1, &origNbr1State);
                      RestoreState(home, nbr2, &origNbr2State);                                                 
                      DefinePos(param,node,nbr1,nbr2,nodep,nbr1p,nbr2p,n1mn2,nmn1,nmn2);

               
                      // Do cross-slip events if found
                      if (DoAddEvent) 
                      {
                         // Add the necessary info to the list of normal
                         // cross-slip events to attempt.  These normal
                         // (non-zipper) type cross-slip events are added
                         // from the top of the event list downward.  The
                         // zipper events are added from the bottome of
                         // the event list upward.
#ifdef _OPENMP
#pragma omp critical (CROSS_SLIP_EVENT)
#endif
                         {
                            slipEventIndex -= 1;
                            thisEventIndex = slipEventIndex;
                         }  // end "omp critical" section //
                         
                         thisEvent = &eventList[thisEventIndex];
                         
                         thisEvent->node = node;
                         thisEvent->nbrNodes[0] = nbr1;
                         thisEvent->nbrNodes[1] = nbr2;
                         
                         VECTOR_COPY(thisEvent->newNodePos, KeepNodep);
                         VECTOR_COPY(thisEvent->newPlane, planeDirLab[Keepglide]);
                         VECTOR_COPY(thisEvent->glideDir, glideDirLab[Keepglide]);
                         
                         thisEvent->fdotglide = Keepfdotglide;
                         thisEvent->segIndex = 1;
                         thisEvent->status = CROSS_SLIP_PENDING;
                         thisEvent->flags = SLIP_REPOSITION_NODE;
                         
                         if (!pinned2) {
                            VECTOR_COPY(thisEvent->newNbrPos, KeepNbrp);
                            thisEvent->flags |= SLIP_REPOSITION_NBR;
                         }
                      }   
                      
                   } 
                   else 
                   {  
                      //  Neighbor 1 can be moved
                      DoAddEvent = 0;

                      //  Align both segments exactly along screw by moving nbr1 and node.
                      n1mn2dotb = DotProduct(n1mn2, burgLab);
                      nbr1p[X] = nbr2p[X] + (n1mn2dotb * burgLab[X]);
                      nbr1p[Y] = nbr2p[Y] + (n1mn2dotb * burgLab[Y]);
                      nbr1p[Z] = nbr2p[Z] + (n1mn2dotb * burgLab[Z]);
                      
                      nmn2dotb = DotProduct(nmn2, burgLab);                         
                      nodep[X] = nbr2p[X] + (nmn2dotb * burgLab[X]);
                      nodep[Y] = nbr2p[Y] + (nmn2dotb * burgLab[Y]);
                      nodep[Z] = nbr2p[Z] + (nmn2dotb * burgLab[Z]);                         
                     
                      FoldBox(param, &nodep[X], &nodep[Y], &nodep[Z]);
                      RepositionNode(home, nodep, &node->myTag, 0);
                      
                      FoldBox(param, &nbr1p[X], &nbr1p[Y], &nbr1p[Z]);
                      RepositionNode(home, nbr1p, &nbr1->myTag, 0);

                      // Evaluate forces of node
                      SetOneNodeForce(home, node);

                      // Force is the same for any of the 6 planes
                      KeepNodepFor[X] = node->fX;
                      KeepNodepFor[Y] = node->fY;
                      KeepNodepFor[Z] = node->fZ;

                      // Calculate power configuration for the 6 glide planes
                      int mobError[6];
                      for (iglide = 0; iglide < 6; iglide++) 
                      {
                         ResetGlidePlane(home, planeDirLab[iglide], &node->myTag,
                                         &nbr1->myTag, 0);
                         ResetGlidePlane(home, planeDirLab[iglide], &node->myTag,
                                         &nbr2->myTag, 0);

                         mobError[iglide]  = EvaluateMobility(home, node, &mobArgs);
  
                         powerTest[iglide] = (node->fX * node->vX) + (node->fY * node->vY) + (node->fZ * node->vZ) + 
                               vNoise * sqrt((node->fX * node->fX) + (node->fY * node->fY) + (node->fZ * node->fZ) );
                      } 

                      // Use highest power dissipation criterion + randomness to choose the glide plane
                      // to cross-slip on
                      Keepglide = PowerChoose(powerTest,units);
                      //if (powerTest[Keepglide] > powerMax && mobError[Keepglide] == 0) 
                      if (mobError[Keepglide] == 0) 
                      {
                         DoAddEvent = 1;

                         // Node position and its neighbor position does not change with glide plane
                         VECTOR_COPY(KeepNodep,nodep);
                         VECTOR_COPY(KeepNbrp,nbr1p);
                      }

                      
                      if (DoAddEvent)
                      {
                         // Cross-slip can happen
                         // Check that for this configuration, the area will grow.
                         // Move the node and evaluate the area growth criterion
 
                         n1mn2dotb = DotProduct(n1mn2, burgLab);
                         fdotglide = DotProduct(KeepNodepFor,glideDirLab[Keepglide]);                        
                         tmp = Sign(fdotglide) * coeff / fabs(n1mn2dotb) ;
                         
                         nodep[X] += tmp * glideDirLab[Keepglide][X];
                         nodep[Y] += tmp * glideDirLab[Keepglide][Y];
                         nodep[Z] += tmp * glideDirLab[Keepglide][Z];
                         
                         FoldBox(param, &nodep[X], &nodep[Y], &nodep[Z]);
                         RepositionNode(home, nodep, &node->myTag, 0);

                         SetOneNodeForce(home, node);
                         SetOneNodeForce(home, nbr1);
                         SetOneNodeForce(home, nbr2);

                         ResetGlidePlane(home, planeDirLab[Keepglide], &node->myTag,
                                         &nbr1->myTag, 0);
                         ResetGlidePlane(home, planeDirLab[Keepglide], &node->myTag,
                                         &nbr2->myTag, 0);

                         int mobError0;
                         mobError0  = EvaluateMobility(home, node, &mobArgs);
                         mobError0 |= EvaluateMobility(home, nbr1, &mobArgs);                        
                         mobError0 |= EvaluateMobility(home, nbr2, &mobArgs);           

                         // Area growth criterion's test
                         real8 v1[3],v[3],v2[3];                         
                         v1[X] = nbr1->vX; v1[Y] = nbr1->vY; v1[Z]= nbr1->vZ;
                         v [X] = node->vX;  v[Y] = node->vY;  v[Z]= node->vZ;
                         v2[X] = nbr2->vX; v2[Y] = nbr2->vY; v2[Z]= nbr2->vZ;
                         real8 dAreadt = IsAreaGrowth(nbr1p,nbr2p,nodep,v1,v2,v);

                         if (dAreadt > 0.0 && mobError0 == 0)
                         {
                            DoAddEvent = 1;                           
                            Keepfdotglide = fdotglide;
                            VECTOR_COPY(KeepNodep,nodep);
                         }
                         else
                            DoAddEvent = 0;
                      }

                      // Restore initial node positions for next trial
                      RestoreState(home, node, &origNodeState);
                      RestoreState(home, nbr1, &origNbr1State);
                      RestoreState(home, nbr2, &origNbr2State);                                                 
                      DefinePos(param,node,nbr1,nbr2,nodep,nbr1p,nbr2p,n1mn2,nmn1,nmn2);

               
                      // Do cross-slip events if found
                      if (DoAddEvent) 
                      {
                         // Add the necessary info to the list of normal
                         // cross-slip events to attempt.  These normal
                         // (non-zipper) type cross-slip events are added
                         // from the top of the event list downward.  The
                         // zipper events are added from the bottome of
                         // the event list upward.
#ifdef _OPENMP
#pragma omp critical (CROSS_SLIP_EVENT)
#endif
                         {
                            slipEventIndex -= 1;
                            thisEventIndex = slipEventIndex;
                         }  // end "omp critical" section //
                         
                         thisEvent = &eventList[thisEventIndex];
                         
                         thisEvent->node = node;
                         thisEvent->nbrNodes[0] = nbr1;
                         thisEvent->nbrNodes[1] = nbr2;
                         
                         VECTOR_COPY(thisEvent->newNodePos, KeepNodep);
                         VECTOR_COPY(thisEvent->newPlane, planeDirLab[Keepglide]);
                         VECTOR_COPY(thisEvent->glideDir, glideDirLab[Keepglide]);
                         VECTOR_COPY(thisEvent->newNbrPos, KeepNbrp);
                         
                         thisEvent->fdotglide = Keepfdotglide;
                         thisEvent->segIndex = 0;
                         thisEvent->status = CROSS_SLIP_PENDING;
                         thisEvent->flags = SLIP_REPOSITION_NODE |
                            SLIP_REPOSITION_NBR;                         
                      }
                      
                   } 
              } // end if both_seg_are_screw
              
//
// ZIPPERS START HERE
//
              
              else if ((seg1_is_screw) && (glide1 != glide2) )
                {
                   real8 Keepnbrp[3];
                   VECTOR_ZERO(Keepnbrp);

                   DoAddEvent = 0;

                   if ((!pinned1) || ((testmax2-test2) < vAligned))
                   {
                      nmn1dotb = DotProduct(nmn1, burgLab);
                      if  (!pinned1) 
                      {
                         nbr1p[X] = nodep[X] - (nmn1dotb * burgLab[X]);
                         nbr1p[Y] = nodep[Y] - (nmn1dotb * burgLab[Y]);
                         nbr1p[Z] = nodep[Z] - (nmn1dotb * burgLab[Z]);
                      }
 
                      FoldBox(param, &nbr1p[X], &nbr1p[Y], &nbr1p[Z]);
                      RepositionNode(home, nbr1p, &nbr1->myTag, 0);
                      
                      // Recalculate force and velocity of node
                      SetOneNodeForce(home, node);
                      
                      ResetGlidePlane(home, planeDirLab[glide2], &node->myTag,
                                      &nbr1->myTag, 0);
                      
                      int mobError0  = EvaluateMobility(home, node, &mobArgs);
                      
                      powerTestZip = (node->fX * node->vX) + 
                                     (node->fY * node->vY) + 
                                     (node->fZ * node->vZ) +
                         vNoise * (sqrt((node->fX * node->fX) + 
                                        (node->fY * node->fY) + 
                                        (node->fZ * node->fZ)));
                      
                      
                      if ( (powerTestZip > powerMax) && (mobError0 == 0)) 
                      {                         
                         DoAddEvent = 1;                        
                         VECTOR_COPY(Keepnbrp,nbr1p);
                      }
                      
                      RestoreState(home, node, &origNodeState);
                      RestoreState(home, nbr1, &origNbr1State);
                      DefinePos(param,node,nbr1,nbr2,nodep,nbr1p,nbr2p,n1mn2,nmn1,nmn2);
                      
                      if (DoAddEvent)
                      {
                         // Add the necessary info to the list of zipper
                         // cross-slip events to attempt.  The zipper
                         // type cross-slip events are added from the bottom
                         // of the event list upward.  The normal (non-zipper)
                         // events are added from the top of the event list
                         // downward.
#ifdef _OPENMP
#pragma omp critical (CROSS_SLIP_ZIP_EVENT)
#endif
                         {
                            thisEventIndex = zipEventCnt;
                            zipEventCnt += 1;
                         }  // end "omp critical" section //
                         
                         thisEvent = &eventList[thisEventIndex];
                         
                         thisEvent->flags = SLIP_ZIPPER_EVENT;
                         thisEvent->status = CROSS_SLIP_PENDING;
                         thisEvent->segIndex = 0;
                         thisEvent->node = node;
                         thisEvent->nbrNodes[0] = nbr1;
                         
                         VECTOR_COPY(thisEvent->newPlane, planeDirLab[glide2]);
                         
                         if (!pinned1) 
                         {
                            thisEvent->flags |= SLIP_REPOSITION_NBR;
                            VECTOR_COPY(thisEvent->newNbrPos, Keepnbrp);
                         }
                      }
                   }


                } 
                else if ((seg2_is_screw) && (glide1 != glide2))
                {
                   real8 Keepnbrp[3];
                   VECTOR_ZERO(Keepnbrp);
                   DoAddEvent = 0;

                   if ((!pinned2) || (testmax2-test2) < vAligned)
                   {

                      nmn2dotb = DotProduct(nmn2, burgLab);
                      if  (!pinned2)
                      {
                         nbr2p[X] = nodep[X] - (nmn2dotb * burgLab[X]);
                         nbr2p[Y] = nodep[Y] - (nmn2dotb * burgLab[Y]);
                         nbr2p[Z] = nodep[Z] - (nmn2dotb * burgLab[Z]);
                      }
                                            
                         FoldBox(param, &nbr2p[X], &nbr2p[Y], &nbr2p[Z]);
                         RepositionNode(home, nbr2p, &nbr2->myTag, 0);
                         
                         // Recalculate force and velocity on the moved nodes
                         SetOneNodeForce(home, node);
                         
                         ResetGlidePlane(home, planeDirLab[glide1], &node->myTag,
                                         &nbr2->myTag, 0);
                         
                         int mobError0  = EvaluateMobility(home, node, &mobArgs);
                         
                         powerTestZip = (node->fX * node->vX) + 
                                        (node->fY * node->vY) + 
                                        (node->fZ * node->vZ) +
                            vNoise * (sqrt((node->fX * node->fX) + 
                                           (node->fY * node->fY) + 
                                           (node->fZ * node->fZ)));
                         
                         
                         if ( (powerTestZip > powerMax) && (mobError0==0)) 
                         {
                            DoAddEvent = 1;  
                            VECTOR_COPY(Keepnbrp,nbr2p);
                         }

                         // Restore initial node positions 
                         RestoreState(home, node, &origNodeState);
                         RestoreState(home, nbr2, &origNbr2State);
                         DefinePos(param,node,nbr1,nbr2,nodep,nbr1p,nbr2p,n1mn2,nmn1,nmn2);
                      

                      if (DoAddEvent)
                      {
                         // Add the necessary info to the list of zipper
                         // cross-slip events to attempt.  The zipper
                         // type cross-slip events are added from the bottom
                         // of the event list upward.  The normal (non-zipper)
                         // events are added from the top of the event list
                         // downward.
#ifdef _OPENMP
#pragma omp critical (CROSS_SLIP_ZIP_EVENT)
#endif
                         {
                            thisEventIndex = zipEventCnt;
                            zipEventCnt += 1;
                         }  // end "omp critical" section //
                         
                         thisEvent = &eventList[thisEventIndex];
                         
                         thisEvent->flags = SLIP_ZIPPER_EVENT;
                         thisEvent->status = CROSS_SLIP_PENDING;
                         thisEvent->segIndex = 1;
                         thisEvent->node = node;
                         thisEvent->nbrNodes[1] = nbr2;
                         
                         VECTOR_COPY(thisEvent->newPlane, planeDirLab[glide1]);
                         
                         if (!pinned2) {
                            thisEvent->flags |= SLIP_REPOSITION_NBR;
                            VECTOR_COPY(thisEvent->newNbrPos, Keepnbrp);
                         }
                         
                      }

                   }  // end if (!pinned2) //
                   
                }  // end if (seg2isScrew) //
             
                RestoreState(home, node, &origNodeState);
                RestoreState(home, nbr1, &origNbr1State);
                RestoreState(home, nbr2, &origNbr2State);   
                                              
		DiscardState(&origNodeState);
                DiscardState(&origNbr1State);
                DiscardState(&origNbr2State);                
		
   
	    }  // end for (nIndex = threadIterStart; ...) //
            
        }  // end "omp parallel" section // 

        
        // It may be that due to conflicts we can not perform all the
        // desired cross-slip events, however, if the conflict is between
        // a 'zipper' event and a normal cross-slip event, we want to 
        // preferentially perform the 'zipper' event. So, we first
        // loop through the 'zipper' events.  Reminder, the zipper events
        // are in the first <zipEventCnt> elements of the event list.
        //
        // NOTE: These loops are not explicitly threaded here for two
        //       reasons.  First, so much of the code in the loops
        //       would have to be in 'critical' sections that
        //       the lock contention would probably slow things down.
        //       Second, the primary expense is the force calculations
        //       in the loop over the non-zipper events, and those
        //       force calculations are already threaded within the
        //       SetOneNodeForce() function.
        //
        for (i = 0; i < zipEventCnt; i++) {
            int segIndex;
            Node_t *node, *nbrNode;

            if (!ProceedWithZipEvent(eventList, zipEventCnt, i)) {
                continue;
            }

            segIndex = eventList[i].segIndex;
            node = eventList[i].node;
            nbrNode = eventList[i].nbrNodes[segIndex];

#ifdef DEBUG_CROSSSLIP_EVENTS
            if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                char *zip1 = "Zipping first segment";
                char *zip2 = "Zipping second segment";
                DumpCrossSlipEvent(node, eventList[i].newPlane,
                                   (segIndex == 0 ? zip1 : zip2));
            }
#endif

            ResetGlidePlane(home, eventList[i].newPlane, &node->myTag,
                            &nbrNode->myTag, 1);

            // If the neighboring node was repositioned, shift the new coordinates
            // back into the primary image space and update the corresponding
            // node structure.  The operation will also be sent to the remote
            // domains for processing.
            if (eventList[i].flags & SLIP_REPOSITION_NBR) {
                real8 newNbrPos[3];
                VECTOR_COPY(newNbrPos, eventList[i].newNbrPos);
                FoldBox(param, &newNbrPos[X], &newNbrPos[Y], &newNbrPos[Z]);
                RepositionNode(home, newNbrPos, &nbrNode->myTag, 1);
            }
        }

        // Now loop through the normal cross-slip events.  Reminder: these
        // are contained within elements slipEventIndex thru eventListSize-1
        // of the event list.
        for (i = slipEventIndex; i < eventListSize; i++) {
            int    nbr1ArmID, nbr2ArmID, segIndex;
            real8  newfdotglide;
            real8  *nbrPosOrig;
            real8  newPlane[3];
            real8  newNodePos[3], newNbrPos[3];
            real8  nodePosOrig[3], nbr1PosOrig[3], nbr2PosOrig[3];
            real8  newForce[3];
            real8  segForceOrig[4][3];
            Node_t *node, *nbr1, *nbr2, *nbrNode;

            // Since a cross-slip event can modify node positions and
            // segments data, we don't want to attempt a cross-slip that
            // involves nodes/segments that were modified in a prior
            // cross-slip event this cycle.
            if (!ProceedWithCrossSlip(eventList, zipEventCnt, slipEventIndex,
                                      eventListSize, i)) {
                continue;
            }

            segIndex = eventList[i].segIndex;
            nbrNode = eventList[i].nbrNodes[segIndex];

            node = eventList[i].node;
            nbr1 = eventList[i].nbrNodes[0];
            nbr2 = eventList[i].nbrNodes[1];

            // It's possible that once we've shifted nodes and recalculated
            // forces after the cross-slip event, it will turn out that the
            // original configuration was preferable.  We need to preserve
            // the original configuration so that it can be restored if this
            // situation occurs.
            nbr1ArmID = GetArmID(nbr1, node);
            nbr2ArmID = GetArmID(nbr2, node);

            SaveCrossSlipInfo(node, nbr1, nbr2, nbr1ArmID, nbr2ArmID,
                              segForceOrig, nodePosOrig, nbr1PosOrig,
                              nbr2PosOrig);

            nbrPosOrig = (segIndex == 0 ? nbr1PosOrig : nbr2PosOrig);

            // Perform the cross-slip event and recalculate the nodal
            // force.  If it appears the node will not continue to move out
            // on the new plane, skip the cross-slip event and
            // restore the old configuration.
            VECTOR_COPY(newNodePos, eventList[i].newNodePos);
            ResetPosition(param, node, newNodePos);

            if (eventList[i].flags & SLIP_REPOSITION_NBR) {
                VECTOR_COPY(newNbrPos, eventList[i].newNbrPos);
                ResetPosition(param, nbrNode, newNbrPos);
            }

            SetOneNodeForce(home, node);

            newForce[X] = node->fX;
            newForce[Y] = node->fY;
            newForce[Z] = node->fZ;

            newfdotglide = DotProduct(newForce, eventList[i].glideDir);

            if ((Sign(newfdotglide) * Sign(eventList[i].fdotglide)) < 0.0) {
                ResetPosition(param, node, nodePosOrig);
                if (eventList[i].flags & SLIP_REPOSITION_NBR) {
                    ResetPosition(param, nbrNode, nbrPosOrig);

                }
                RestoreCrossSlipForce(node, nbr1, nbr2,
                                      nbr1ArmID, nbr2ArmID,
                                      segForceOrig);
                eventList[i].status = CROSS_SLIP_BACKED_OFF;
                continue;
            }

            VECTOR_COPY(newPlane, eventList[i].newPlane);

#ifdef DEBUG_CROSSSLIP_EVENTS
            if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                char msg[128];

                if (eventList[i].flags & SLIP_REPOSITION_NBR)  {
                    sprintf(msg, "%s neighbor repositioned",
                            segIndex == 0 ? "first" : "second");
                } else {
                    sprintf(msg, "neither neighbor repositioned");
                }
                DumpCrossSlipEvent(node, newPlane, msg);
            }
#endif
            ResetGlidePlane(home, newPlane, &node->myTag, &nbr1->myTag, 1);
            ResetGlidePlane(home, newPlane, &node->myTag, &nbr2->myTag, 1);

            FoldBox(param, &newNodePos[X], &newNodePos[Y], &newNodePos[Z]);
            RepositionNode(home, newNodePos, &node->myTag, 1);

            if (eventList[i].flags & SLIP_REPOSITION_NBR) {
                FoldBox(param, &newNbrPos[X], &newNbrPos[Y], &newNbrPos[Z]);
                RepositionNode(home, newNbrPos, &nbrNode->myTag, 1);
            }

/*
 *          Update velocity of the moved nodes
 */
            int mobError;
            MobArgs_t  mobArgs;
            mobError  = EvaluateMobility(home, node, &mobArgs);
            mobError |= EvaluateMobility(home, nbr1, &mobArgs);
            mobError |= EvaluateMobility(home, nbr2, &mobArgs);
        }

        if (eventList != (SlipEvent_t *)NULL) {
            free(eventList);
        }

        return;
}
