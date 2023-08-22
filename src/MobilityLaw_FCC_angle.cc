/**************************************************************************
 *
 *      Module:       Mobility_FCC_angle.c
 *      Description:  Contains functions for calculating generic isotropic
 *                    linear mobility for FCC materials in a non-linear
 *                    and character angles. Employed mobiity law is 
 *                    so-called Eshelby-Olmsted mobility law, and the 
 *                    law can be fitted based on MD data.
 *      Literatures:  1. Mobility law of dislocations with several
 *                       character angles and temperatures in FCC
 *                       aluminum (https://doi.org/10.1016/j.ijplas.2016.12.004)
 *                    2. Coupled 3D Dislocation Modeling at Nano-
 *                       and Micro-Scales (http://infoscience.epfl.ch/record/225968)
 *      Author   :    Jaehyun Cho
 *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Mobility.h"
#include "V3.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

//#define M33_INVERSE(a,b) Matrix33_Inverse(a,b)
//#define M33_INVERSE(a,b) Matrix33_SVD_Inverse(a,b)
  #define M33_INVERSE(a,b) Matrix33_PseudoInverse(a,b)

#if 0
static void FindThreeMaximum(real8 *seg, real8 *lineRef, 
			     real8 *Ar, int numElements, 
			     int *indexOfMax1, 
			     int *indexOfMax2, 
			     int *indexOfMax3)
{
        *indexOfMax1 = -1;
        *indexOfMax2 = -1;
        *indexOfMax3 = -1;

	real8 maxA1 = 0.0;
	real8 maxA2 = 0.0;
	real8 maxA3 = 0.0;
	
	real8 *absAr = new real8 [numElements];
	for (int i=0; i<numElements; i++)
	  absAr[i] = fabs(Ar[i]);
	
	FindMaxExcluding(absAr, numElements, &maxA1, indexOfMax1);

	FindMaxExcluding(absAr, numElements, &maxA2, indexOfMax2, 
			 *indexOfMax1);

	FindMaxExcluding(absAr, numElements, &maxA3, indexOfMax3, 
			 *indexOfMax1, *indexOfMax2);


	delete [] absAr;

        return;
}
#endif

static bool DecomposeArm(real8  seg[3], 
			 int    numThetasGlides, 
			 real8 *lineRef,
			 real8 *theta1Length, 
			 real8 *theta2Length,
			 real8 *theta3Length,
			 int   *theta1Index, 
			 int   *theta2Index,
			 int   *theta3Index)
{
	*theta1Index = -1; *theta1Length = 0.0;
	*theta2Index = -1; *theta2Length = 0.0;
	*theta3Index = -1; *theta3Length = 0.0;

        int    i;
	real8  thetaLine[3];
	real8 *LdotE     = new real8  [numThetasGlides];
	real8 *LdotEAbs  = new real8  [numThetasGlides];

	real8  Lthetas[3]; V3_ZERO(Lthetas);

/*
 *      Now get the three maximum projection lengths with the 
 *      given character angles and glide planes
 */
        for (i = 0; i < numThetasGlides; i++){
	  thetaLine[0] = lineRef[3*i+0];
	  thetaLine[1] = lineRef[3*i+1];
	  thetaLine[2] = lineRef[3*i+2];
	  LdotE[i] = DotProduct(thetaLine, seg);
	  LdotEAbs[i] = fabs(DotProduct(thetaLine, seg));
	}

	// fine two maximums
	real8 maxL1, maxL2;
	int   Index1, Index2;
	FindMaxExcluding(LdotEAbs, numThetasGlides, &maxL1, &Index1);
	FindMaxExcluding(LdotEAbs, numThetasGlides, &maxL2, &Index2, Index1);

	int   Index3 = -1;
	real8 EEMat[3][3];    Matrix33_Zero(EEMat);
	real8 EEMatInv[3][3]; Matrix33_Zero(EEMatInv);

	real8 e1[3] = {lineRef[3*Index1+0], lineRef[3*Index1+1], lineRef[3*Index1+2]};
	real8 e2[3] = {lineRef[3*Index2+0], lineRef[3*Index2+1], lineRef[3*Index2+2]};	

	for (i = 0; i < numThetasGlides; i++){
	  if( (i!=Index1) && (i!=Index2) ) {

	    Index3 = i;
	    real8 e3[3] = {lineRef[3*Index3+0], lineRef[3*Index3+1], lineRef[3*Index3+2]};

	    EEMat[0][0] = DotProduct(e1,e1);  EEMat[0][1] = DotProduct(e1,e2);  EEMat[0][2] = DotProduct(e1,e3);
	    EEMat[1][0] = DotProduct(e2,e1);  EEMat[1][1] = DotProduct(e2,e2);  EEMat[1][2] = DotProduct(e2,e3);
	    EEMat[2][0] = DotProduct(e3,e1);  EEMat[2][1] = DotProduct(e3,e2);  EEMat[2][2] = DotProduct(e3,e3);
	    real8 Le[3] = {DotProduct(e1,seg), DotProduct(e2,seg), DotProduct(e3,seg)};

	    if( M33_INVERSE(EEMatInv,EEMat) > 0 ) 
        {	
	      Matrix33Vector3Multiply(EEMatInv, Le, Lthetas);

	      real8 test_seg[3];
	      test_seg[0] = Lthetas[0]*e1[0] + Lthetas[1]*e2[0] + Lthetas[2]*e3[0];
	      test_seg[1] = Lthetas[0]*e1[1] + Lthetas[1]*e2[1] + Lthetas[2]*e3[1];
	      test_seg[2] = Lthetas[0]*e1[2] + Lthetas[1]*e2[2] + Lthetas[2]*e3[2];

          if( V3_NEAR(seg, test_seg, 1e-4)) 
          {
            *theta1Index  = Index1;       *theta2Index  = Index2;      *theta3Index  = Index3;
            *theta1Length = Lthetas[0];   *theta2Length = Lthetas[1];  *theta3Length = Lthetas[2];

            if( fabs(Lthetas[0]) < fabs(Lthetas[1]) ) 
            {
              *theta1Index  = Index2;
              *theta1Length = Lthetas[1];

              *theta2Index  = Index1;
              *theta2Length = Lthetas[0];
            } 

            delete [] LdotE;              
            delete [] LdotEAbs;
            return true;
          }
	    } 
	  }
	}	  
	Fatal("cannot decompose!");
    return(false);

	// if(maxL1 > maxL2) {
	//   *theta1Length = LdotE[Index1];
	//   *theta2Length = LdotE[Index2];
	//   *theta1Index = Index1;
	//   *theta2Index = Index2;
	// } else {
	//   *theta1Length = LdotE[Index2];
	//   *theta2Length = LdotE[Index1];
	//   *theta1Index = Index2;
	//   *theta2Index = Index1;
	// }

	// if(Index3 != -1) {
	//   *theta3Length = LdotE[Index3];
	//   *theta3Index = Index3;
	// } else {
	//   *theta3Length = 0.0;
	// }

	//IterativeInverseMatrix(numThetasGlides, EEmat, EEMatInv);

	//FindThreeMaximum(seg, lineRef, tempLength,  numThetasGlides, 
	//			 theta1Index,  theta2Index, theta3Index);

/*
 *      Three angles are used to decompose the seg
 */
	// real8 LdotE[3] = {tempLength[*theta1Index], 
	// 		  tempLength[*theta2Index], 
	// 		  tempLength[*theta3Index]};

	//real8 EEMat[3][3];
	//real8 EEMatInv[3][3];
	//Matrix33_Zero(EEMat);
	//Matrix33_Zero(EEMatInv);
	
	// real8 e1[3] = {lineRef[3*(*theta1Index)+0], 
	// 	       lineRef[3*(*theta1Index)+1],
	// 	       lineRef[3*(*theta1Index)+2]};
	// real8 e2[3] = {lineRef[3*(*theta2Index)+0], 
	// 	       lineRef[3*(*theta2Index)+1],
	// 	       lineRef[3*(*theta2Index)+2]};
	// real8 e3[3] = {lineRef[3*(*theta3Index)+0], 
	// 	       lineRef[3*(*theta3Index)+1],
	// 	       lineRef[3*(*theta3Index)+2]};
	
	// EEMat[0][0] = DotProduct(e1,e1);  EEMat[0][1] = DotProduct(e1,e2);  EEMat[0][2] = DotProduct(e1,e3);
	// EEMat[1][0] = DotProduct(e2,e1);  EEMat[1][1] = DotProduct(e2,e2);  EEMat[1][2] = DotProduct(e2,e3);
	// EEMat[2][0] = DotProduct(e3,e1);  EEMat[2][1] = DotProduct(e3,e2);  EEMat[2][2] = DotProduct(e3,e3);
	
	// //if ( M33_INVERSE(EEMatInv, EEMat) < 0 ) {
	// //  Fatal("cannot invert EEmat! probably the three maximum are not correct");
	// //}

	// real8 Lthetas[3];
	// //Matrix33Vector3Multiply(EEMatInv, LdotE, Lthetas);

	
	// real8 LthetasAbs[3];
	// V3_SET(LthetasAbs, fabs(Lthetas[0]), fabs(Lthetas[1]), fabs(Lthetas[2]));

	// int Indices[3];
	// V3_SET(Indices, *theta1Index, *theta2Index, *theta3Index);

	// real8 maxL1 = 0.0;	real8 maxL2 = 0.0;	real8 maxL3 = 0.0;
	// int Index1, Index2, Index3;
	// FindMaxExcluding(LthetasAbs, 3, &maxL1, &Index1);
	// FindMaxExcluding(LthetasAbs, 3, &maxL2, &Index2, Index1);
	// FindMaxExcluding(LthetasAbs, 3, &maxL3, &Index3, Index1, Index2);

	// // reconstruction of seg with line1, line2 and line3
	// real8 test_seg[3]; V3_ZERO(test_seg);

	// if(fabs(Lthetas[Index1]) < 1e-4){
	//   *theta1Length = 0.0;
	//   *theta1Index = -1;
	// } else {
	//   *theta1Length = Lthetas[Index1];
	//   *theta1Index = Indices[Index1];

	//   test_seg[0] += (*theta1Length)*lineRef[3*(*theta1Index)+0];
	//   test_seg[1] += (*theta1Length)*lineRef[3*(*theta1Index)+1];
	//   test_seg[2] += (*theta1Length)*lineRef[3*(*theta1Index)+2];
	// }

	// if(fabs(Lthetas[Index2]) < 1e-4){
	//   *theta2Length = 0.0;
	//   *theta2Index = -1;
	// } else {
	//   *theta2Length = Lthetas[Index2];
	//   *theta2Index = Indices[Index2];

	//   test_seg[0] += (*theta2Length)*lineRef[3*(*theta2Index)+0];
	//   test_seg[1] += (*theta2Length)*lineRef[3*(*theta2Index)+1];
	//   test_seg[2] += (*theta2Length)*lineRef[3*(*theta2Index)+2];
	// }

	// if(fabs(Lthetas[Index3]) < 1e-4){
	//   *theta3Length = 0.0;
	//   *theta3Index = -1;
	// } else {
	//   *theta3Length = Lthetas[Index3];
	//   *theta3Index = Indices[Index3];

	//   test_seg[0] += (*theta3Length)*lineRef[3*(*theta3Index)+0];
	//   test_seg[1] += (*theta3Length)*lineRef[3*(*theta3Index)+1];
	//   test_seg[2] += (*theta3Length)*lineRef[3*(*theta3Index)+2];
	// }
	
	// if (*theta1Index == -1) {
	//   Print3("seg",seg);
	//   printf("theta1Index  %i, theta2Index  %i, theta3Index  %i \n", 
	// 	 *theta1Index,  *theta2Index,  *theta3Index);
	//   printf("theta1Length %f, theta2Length %f, theta3Length %f \n", 
	// 	 *theta1Length, *theta2Length, *theta3Length);
	//   return false;
	// }

	// // release data
	// delete [] LdotE;
	// for (i = 0; i < numThetasGlides; i++){
	//   delete [] EEMat[i];
	//   delete [] EEMatInv[i];
	// }
	// delete [] EEMat;
	// delete [] EEMatInv;
	

	// if(V3_NEAR(test_seg, seg, 1e-4))
	//   return true;
	// Fatal("cannot decompose!");
}


#if 0
bool CheckLineDirection(int numNbrs,
			real8 *segTotLength,
			real8 angle1Dir[MAX_NBRS][3],
			real8 angle2Dir[MAX_NBRS][3],
			real8 angle3Dir[MAX_NBRS][3],
			char  *msg)
{
  printf("%s = \n", msg);

  real8 thetaLine1[3], thetaLine2[3], thetaLine3[3];
  for (int i = 0; i < numNbrs; i++) {
    if (fabs(segTotLength[i]) < 1e-12) {
      continue;
    }
    
    thetaLine1[0] = angle1Dir[i][0];
    thetaLine1[1] = angle1Dir[i][1];
    thetaLine1[2] = angle1Dir[i][2];

    printf("i=%i\n", i);
    Print3("thetaLine1", thetaLine1);
    Print3("thetaLine2", thetaLine2);
    Print3("thetaLine3", thetaLine3);

  }

  return(false);
}
#endif

/**************************************************************************
 *
 *      Function:     MobilityLaw_FCC_angle
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *************************************************************************/
int Mobility_FCC_angle(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int   i, m, n;
        int   count, maxIterCount;
        int   numNbrs, shortSegCutoff, typePlane1 = 0;
        int   numNonZeroLenSegs=0;
        int   firstPlaneIndex, numPlanes;
        int  *numPlanesPerBurg, *numGlissilePlanesPerBurg;
        int  *burgFirstPlaneIndex, *planeType;

        int   angle1Exists[MAX_NBRS];
        int   angle2Exists[MAX_NBRS];
        int   angle3Exists[MAX_NBRS];
        for (i=0; i < MAX_NBRS; i++){
        angle1Exists[i] = 0;	  angle2Exists[i] = 0;	  angle3Exists[i] = 0;
        }

        int   junctionExists[MAX_NBRS];
        real8 shortSegRatio = 1.0;
        real8 mu, maxSegLen, maxVal = 0;
        real8 massMult=0, massMatrix[3][3], eps;
        real8 vTest[3], vOld[3], fn[3], rt[3], plane[3];
        real8 (*burgList)[3], (*planeList)[3];

        real8 (*lineList)[3];
        int   lPlaneIndex[17];

        real8 segTotLength[MAX_NBRS];

        real8 burg[MAX_NBRS][3], lineDir[MAX_NBRS][3], glideDir[MAX_NBRS][3];

        real8 dragLine, dragClimb;
        real8 *dragGlide_B, *dragGlide_D, *dragGlide_V0;

        real8 *planeListXplane;

        real8 fError[3], dfErrordv[3][3];
        real8 fErrorTol;

        real8 forceError, totLen;
        real8 correction[3], inertia[3];
        Param_t *param;
        Node_t  *nbrNode;
        
        int    angle1Index, angle2Index, angle3Index;
        real8  angle1, angle2, angle3;
        real8  thetaLine1[3], thetaLine2[3], thetaLine3[3];

        int    angle1Indices[MAX_NBRS],   angle2Indices[MAX_NBRS],   angle3Indices[MAX_NBRS];
        real8  angle1Dir[MAX_NBRS][3],    angle2Dir[MAX_NBRS][3],    angle3Dir[MAX_NBRS][3];
        real8  segAngle1Length[MAX_NBRS], segAngle2Length[MAX_NBRS], segAngle3Length[MAX_NBRS];

        int    numThetas = home->param->MobTheta_cnt;
        real8 *thetas    = home->param->MobThetas;

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

        param = home->param;

        mu = param->shearModulus;
        maxSegLen = param->maxSeg;

        real8 (*invdfErrordv)[3]     =  mobArgs->invDragMatrix;
     // real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invdfErrordv);
        *numGlideConstraints = 0;

        eps = 1.0e-12;

/*     
 *      Define drag coefficients : chosen similar to FCC_0b.
 */
        
	dragGlide_B = new real8 [numThetas];
	dragGlide_D = new real8 [numThetas];
	dragGlide_V0= new real8 [numThetas];

	real8 minB = 1e+8;
	real8 maxB = 0.0;
	for (i=0; i < numThetas; i++){
	  dragGlide_B[i] = 1.0/(param->MobGlide_B[i]+1e-20);
	  dragGlide_D[i] = 1.0/(param->MobGlide_D[i]+1e-20);
	  dragGlide_V0[i]= param->MobGlide_V0[i];
	  if(dragGlide_B[i] < minB) {
	    minB = dragGlide_B[i];
	  }
	  if(dragGlide_B[i] > maxB) {
	    maxB = dragGlide_B[i];
	  }
	}
	dragLine  = minB * 1e-2; 
        dragClimb = maxB * 1e+2;

/*
 *      If we need to include inertial terms, need some more initializations
 */
        if (param->includeInertia != 0) {
            massMult = 0.25 * param->massDensity *
                       (param->burgMag * param->burgMag);
        }

	shortSegCutoff = 0.5 * param->minSeg;
        numNbrs = node->numNbrs;

/*
 *      This function saves various pieces of info for each neighbor
 *      node using static arrays for speed rather than doing memory
 *      allocations.  If the arrays are too  small, abort with an error
 *      to let the caller increase the array sizes.
 */
        if (numNbrs > MAX_NBRS) {
            Fatal("Mobility_FCC_linear: Node segment count (%d) exceeds\n"
                  "maximum allowed (%d) in static arrays.  Increase MAX_NBRS\n"
                  "and recompile\n", node->numNbrs, MAX_NBRS);
        }


/*
 *      Set some pointers to the reference burgers vector and glide plane
 *      lists that were created during initialization.
 */

        burgList = home->burgData.burgList;
        planeList = home->burgData.planeList;

        lineList = home->burgData.lineList;

     // real8 angleList[ANGLE_MOB_CNT];
     // for (i=0; i<ANGLE_MOB_CNT; i++) angleList[i]   = home->burgData.angleList[i];
        for (i=0; i<17;            i++) lPlaneIndex[i] = home->burgData.lPlaneIndex[i];

        numPlanesPerBurg = home->burgData.numPlanesPerBurg;
        numGlissilePlanesPerBurg = home->burgData.numGlissilePlanesPerBurg;
        burgFirstPlaneIndex = home->burgData.burgFirstPlaneIndex;
        planeType = home->burgData.planeType;

/*
 *      Loop over all arms attached to the node to get some information
 *      for the Newton-Raphson procedure
 */
        totLen = 0.0;

        for (i = 0; i < numNbrs; i++) {
            int   bIndex=-1, pIndex=-1;
            real8 rtSquared;

            nbrNode = GetNeighborNode(home, node, i);

            if (nbrNode == (Node_t *)NULL) {
                segTotLength[i] = 0.0;
                printf("WARNING: Neighbor not found at %s line %d\n",
                       __FILE__, __LINE__);
                continue;
            }

/*
 *          Get the segment line 
 */
            rt[0] = nbrNode->x - node->x;
            rt[1] = nbrNode->y - node->y;
            rt[2] = nbrNode->z - node->z;

            ZImage(param, &rt[0], &rt[1], &rt[2]);

/*
 *          If the segment is zero length (which may happen when the
 *          mobility function is called from SplitMultiNodes()), just
 *          skip the segment.
 */
	    rtSquared = DotProduct(rt, rt);

            if (rtSquared < eps) {
                segTotLength[i] = 0.0;
                continue;
            }

            numNonZeroLenSegs++;

/*
 *         If the node has a very short segment, we need to flag it
 *         for later
 */
	    real8 mag;
	    mag     = sqrt(rtSquared);

	    if (mag < shortSegCutoff) {
	      shortSegRatio = MIN(shortSegRatio, mag/param->minSeg);
	    }


/*
 *          Save the burgers vector and glide plane of the segment
 */
            burg[i][X] = node->burgX[i];
            burg[i][Y] = node->burgY[i];
            burg[i][Z] = node->burgZ[i];

            plane[X] = node->nx[i];
            plane[Y] = node->ny[i];
            plane[Z] = node->nz[i];
            NormalizeVec(plane);

            VECTOR_COPY(glideDir[i],plane);
/*
 *          If necessary, rotate the burgers vector and line sense from the
 *          laboratory frame to the crystal frame
 */
            if (param->useLabFrame) {
                real8 burgCryst[3], rtCryst[3], planeCryst[3];

                Matrix33Vector3Multiply(param->rotMatrixInverse, burg[i]    , burgCryst );
                Matrix33Vector3Multiply(param->rotMatrixInverse, rt         , rtCryst   );
                Matrix33Vector3Multiply(param->rotMatrixInverse, glideDir[i], planeCryst);
                VECTOR_COPY(burg[i], burgCryst);
                VECTOR_COPY(rt, rtCryst);
                VECTOR_COPY(glideDir[i], planeCryst);
            }

/*
 *          Calculate the segment total length
 */
            segTotLength[i] = sqrt(rtSquared);

            totLen += segTotLength[i];

/* 
 *          Identify if the dislocation is a junction
 *          by checking the index of its Burgers vector
 *          and plane in the list.
 *          typePlane 0 corresponds to a sessile plane
 */
            GetBurgIndexFCC(burg[i],FCC_NUM_TOT_BURG, burgList, &bIndex);

	    pIndex = -1;
	    typePlane1 = 0;

	    if (bIndex < 6) {
	      // Dislocation can be either a junction or not, depending on
	      // its glide plane
	      numPlanes  = numPlanesPerBurg[bIndex];
	      planeListXplane = (real8 *)calloc(1, numPlanes * sizeof(real8));
	      
	      firstPlaneIndex = burgFirstPlaneIndex[bIndex];
	      
	      MatrixMult(&planeList[firstPlaneIndex][0],
			 numPlanes, 3, 3,
			 plane, 1, 1,
			 planeListXplane, 1);
	      
	      FindAbsMax(planeListXplane, numPlanes, &maxVal,
			 &pIndex);
	      
	      
	      free(planeListXplane);
	      
	      pIndex += burgFirstPlaneIndex[bIndex];
	      typePlane1 = planeType[pIndex];
	    }

	    junctionExists[i] = 0;
	    angle1Exists[i] = 0;
	    angle2Exists[i] = 0;
	    angle3Exists[i] = 0;

/*
 *	    A dislocation segment is a junction if both its 
 *          Burgers vector and its planes are pre-defined as 
 *          junction Burgers vector and sessile plane.
 *          A bIndex < 6 can have a sessile plane too.
 */
            if (bIndex >= 6) junctionExists[i] = 1;
            else if  (typePlane1 == 0) junctionExists[i] = 1;
            else {
		
/*
 *              Not a junction
 */
        int start, end, j, m;
        int numGlidingPlanes;
        real8 burgRef[3];
     // real8 planeRef[2][3];
        
        start = burgFirstPlaneIndex[bIndex];
        end   = burgFirstPlaneIndex[bIndex] + numGlissilePlanesPerBurg[bIndex];

        VECTOR_COPY(burgRef, burgList[bIndex]);
        NormalizeVec(burgRef);

        numGlidingPlanes = numGlissilePlanesPerBurg[bIndex];
        
        real8 *lineRef  = new real8 [numThetas*numGlidingPlanes*3];
        real8 *angleRef = new real8 [numThetas*numGlidingPlanes];

        int    count=0;

        for (m = 0, j = start; j < end; j++) 
        {
            // VECTOR_COPY(planeRef[m], planeList[j]);
		    
		    // for each burg and gliding plane, evaluate all line
		    // directions of given mixed dislocations
		    for (n = 0; n < numThetas; n++){		      
		      real8 angle = thetas[n];
		      real8 line[3];  V3_ZERO(line);
		      int   lpindex = lPlaneIndex[j]+n;
		      if(lpindex == -1)
			Fatal("invalide index!");

		      V3_COPY(line, lineList[lpindex]);

		      //MixedDislocationLine(bIndex, m, n, lineList, line);

		      // check the line is already registered or not
		      // (to avoid double registrations of screw dislocation)
		      bool registered = false;
		      for (int l=0; l<count; l++) {
			real8 registered_line[3];
			V3_SET(registered_line, lineRef[3*l+0], lineRef[3*l+1], lineRef[3*l+2]);
			if(V3_NEAR(registered_line, line, 1e-8) == 1){
			  registered = true;
			  break;
			}
		      }

		      // add if it is not registered before
		      if(!registered) {
			angleRef[count] = angle;
			lineRef[3*count + 0] = line[0];
			lineRef[3*count + 1] = line[1];
			lineRef[3*count + 2] = line[2];
			count++;
		      } 
		    }

                    m++;
                }
		
/*
 *              Decompose segment into four angles.
 */
		bool success = false;
                success = DecomposeArm(rt, 
				       count,
				       lineRef,
				       &(segAngle1Length[i]),
				       &(segAngle2Length[i]),
				       &(segAngle3Length[i]),
				       &angle1Index,
				       &angle2Index,
				       &angle3Index
				       );
		if(!success){
		  PrintNode(node);
		  PrintNode(nbrNode);
		  Fatal("fail to decompose the segment (%.5e %.5e %.5e)", rt[0], rt[1], rt[2]);
		} 

		V3_ZERO(angle1Dir[i]);
		V3_ZERO(angle2Dir[i]);
		V3_ZERO(angle3Dir[i]);

		if(angle1Index != -1) {
		  angle1 = angleRef[angle1Index];
		  for (j=0; j<numThetas; j++) {
		    if(fabs(thetas[j]-angle1) < 1e-4) {
		      angle1Indices[i] = j;
		      break;
		    }
		  }

		  thetaLine1[0] = lineRef[3*angle1Index + 0];
		  thetaLine1[1] = lineRef[3*angle1Index + 1];
		  thetaLine1[2] = lineRef[3*angle1Index + 2];

		  VECTOR_COPY(angle1Dir[i], thetaLine1);
		  angle1Exists[i] = 1;
		} else {
		  angle1Exists[i] = 0;
		}
		  

		if(angle2Index != -1) {
		  angle2 = angleRef[angle2Index];
		  for (j=0; j<numThetas; j++) {
		    if(fabs(thetas[j]-angle2) < 1e-4) {
		      angle2Indices[i] = j;
		      break;
		    }
		  }

		  thetaLine2[0] = lineRef[3*angle2Index + 0];
		  thetaLine2[1] = lineRef[3*angle2Index + 1];
		  thetaLine2[2] = lineRef[3*angle2Index + 2];

		  VECTOR_COPY(angle2Dir[i], thetaLine2);
		  angle2Exists[i] = 1;
		} else {
		  angle2Exists[i] = 0;
		}

		if(angle3Index != -1) {
		  angle3 = angleRef[angle3Index];
		  for (j=0; j<numThetas; j++) {
		    if(fabs(thetas[j]-angle3) < 1e-4) {
		      angle3Indices[i] = j;
		      break;
		    }
		  }

		  thetaLine3[0] = lineRef[3*angle3Index + 0];
		  thetaLine3[1] = lineRef[3*angle3Index + 1];
		  thetaLine3[2] = lineRef[3*angle3Index + 2];

		  VECTOR_COPY(angle3Dir[i], thetaLine3);		  
		  angle3Exists[i] = 1;

		} else {
		  angle3Exists[i] = 0;
		}
		
		delete [] lineRef;
		delete [] angleRef;
		
            } /* end glide Burgers vector */

/*
 *          Calculate tangents and climb directions
 */
            lineDir[i][0] = rt[0] / (segTotLength[i]+1e-15);
            lineDir[i][1] = rt[1] / (segTotLength[i]+1e-15);
            lineDir[i][2] = rt[2] / (segTotLength[i]+1e-15);

        }  /* loop over neighbors */


	// CheckLineDirection(numNbrs,segTotLength,angle1Exists,angle2Exists,angle3Exists,angle1Dir,angle2Dir,angle3Dir,"loc1");
		
	/*
 *      It's possible this function was called for a node which only
 *      had zero length segments (during SplitSurfaceNodes() for example).
 *      If that is the case, just set the velocity to zero and return;
 */
        if (numNonZeroLenSegs == 0) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      Set up the mass matrix: if not needed, zero it out.
 */
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                massMatrix[m][n] = 0.0;
                if ((m == n) && (param->includeInertia)) {
                    massMatrix[m][n] = totLen * massMult / param->deltaTT;
                }
            }
        }

/* 
 *      Setup for the Newton-Raphson loop
 */
        count = 0;
        maxIterCount = 100000; 

        fn[0] = node->fX;
        fn[1] = node->fY;
        fn[2] = node->fZ;

        vTest[0] = node->vX;
        vTest[1] = node->vY;
        vTest[2] = node->vZ;

        VECTOR_COPY(vOld, vTest);

        fErrorTol = MAX((DotProduct(fn, fn) * 1.0e-16),
                        (mu * mu * maxSegLen * maxSegLen * 1.0e-24));

/*
 *      If necessary, rotate the force and velocity vectors from the
 *      laboratory frame to the crystal frame
 */
        if (param->useLabFrame) {
            real8 fnCryst[3], vTestCryst[3];
            Matrix33Vector3Multiply(param->rotMatrixInverse, fn, fnCryst);
            Matrix33Vector3Multiply(param->rotMatrixInverse, vTest, vTestCryst);
            VECTOR_COPY(fn, fnCryst);
            VECTOR_COPY(vTest, vTestCryst);
            VECTOR_COPY(vOld, vTest);
        }

/*
 *      If we need to include inertial terms...
 */
        if (param->includeInertia) {
            real8 vDiff[3];
            vDiff[0] = vTest[0] - vOld[0];
            vDiff[1] = vTest[1] - vOld[1];
            vDiff[2] = vTest[2] - vOld[2];
            Matrix33Vector3Multiply(massMatrix, vDiff, inertia);
        } else {
            VECTOR_ZERO(inertia);
        }

        forceError = fErrorTol + 1.0;

/*
 *      Begin the Newton-Raphson loop
 */
	// CheckLineDirection(numNbrs,segTotLength,angle1Exists,angle2Exists,angle3Exists,angle1Dir,angle2Dir,angle3Dir,"loc2");

        while ((count < maxIterCount) && (forceError > fErrorTol)) {

            for (m = 0; m < 3; m++) {
                fError[m] = fn[m] - inertia[m];
            }
    
            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    dfErrordv[m][n] = -massMatrix[m][n];
                }
            }

/*
 *          Loop over all segments (skipping zero length segments)
 *          adding up screw, edge and junction contributions to fError
 *          and dfErrordv.
 */
            for (i = 0; i < numNbrs; i++) {
/*
 *              Skip zero length segments or NULL neighbor
 */
                if (fabs(segTotLength[i]) < eps) {
                    continue;
                }

                real8 fJunction[3], dfJunctiondv[3][3];
                real8 fAngle1[3], dfAngle1dv[3][3];
                real8 fAngle2[3], dfAngle2dv[3][3];
                real8 fAngle3[3], dfAngle3dv[3][3];

      
/*
 *              Junction contribution
 */
                if (junctionExists[i]) {
                    JunctionDrag(vTest, lineDir[i], dragLine,
                                 dragClimb, fJunction, dfJunctiondv);

                    for (m = 0; m < 3; m++) {

 		        fError[m] -= 0.5 * fabs(segTotLength[i]) * fJunction[m];

                        for (n = 0; n < 3; n++) {
			  dfErrordv[m][n] -= 0.5 * fabs(segTotLength[i]) *
			                     dfJunctiondv[m][n];
                        }

                    }

                } else {  /* not a junction */

/*
 *                  Angle1 contribution
 */

                    if (angle1Exists[i]) {
		       thetaLine1[0] = angle1Dir[i][0];
		       thetaLine1[1] = angle1Dir[i][1];
		       thetaLine1[2] = angle1Dir[i][2];

                       AngleDrag(vTest, thetaLine1, 
				 glideDir[i], dragClimb, dragLine, 
				 dragGlide_B[angle1Indices[i]], // B 
				 dragGlide_D[angle1Indices[i]], // D
				 dragGlide_V0[angle1Indices[i]], // V0
				 fAngle1, 
				 dfAngle1dv);

		       for (m = 0; m < 3; m++) {
			 
			 fError[m] -= 0.5 * fabs(segAngle1Length[i]) * fAngle1[m];
			 
			 for (n = 0; n < 3; n++) {
			   dfErrordv[m][n] -= 0.5 * fabs(segAngle1Length[i]) *
			     dfAngle1dv[m][n];
			 }

		       }


                    }

/*
 *                  Angle2 contribution
 */

                    if (angle2Exists[i]) {
		       thetaLine2[0] = angle2Dir[i][0];
		       thetaLine2[1] = angle2Dir[i][1];
		       thetaLine2[2] = angle2Dir[i][2];
		       
                       AngleDrag(vTest, 
				 thetaLine2, 
				 glideDir[i], 
				 dragClimb, 
				 dragLine, 
				 dragGlide_B[angle2Indices[i]], // B 
				 dragGlide_D[angle2Indices[i]], // D
				 dragGlide_V0[angle2Indices[i]], // V0
				 fAngle2, 
				 dfAngle2dv);

		       for (m = 0; m < 3; m++) {
			 
			 fError[m] -= 0.5 * fabs(segAngle2Length[i]) * fAngle2[m];
			 
			 for (n = 0; n < 3; n++) {
			   dfErrordv[m][n] -= 0.5 * fabs(segAngle2Length[i]) *
			     dfAngle2dv[m][n];
			 }

		       }


                    }

/*
 *                  Angle3 contribution
 */
                    if (angle3Exists[i]) {

		       thetaLine3[0] = angle3Dir[i][0];
		       thetaLine3[1] = angle3Dir[i][1];
		       thetaLine3[2] = angle3Dir[1][2];

                       AngleDrag(vTest, 
				 thetaLine3, 
				 glideDir[i], 
				 dragClimb, 
				 dragLine, 
				 dragGlide_B[angle3Indices[i]], // B 
				 dragGlide_D[angle3Indices[i]], // D
				 dragGlide_V0[angle3Indices[i]], // V0
				 fAngle3, 
				 dfAngle3dv);

		       for (m = 0; m < 3; m++) {
			 
			 fError[m] -= 0.5 * fabs(segAngle3Length[i]) * fAngle3[m];
			 
			 for (n = 0; n < 3; n++) {
			   dfErrordv[m][n] -= 0.5 * fabs(segAngle3Length[i]) *
			     dfAngle3dv[m][n];
			 }			 
		       }
                    }
		    
                } /* fError and dfErrordv are computed now */

            } /* end loop over neighbors */

            if ( M33_INVERSE(invdfErrordv,dfErrordv) < 0 ) 
            {
                Print3x3("dfErrordv", dfErrordv);
                Fatal("Mobility_FCC_angle: Cannot invert dfErrordv!");
            }

            Matrix33Vector3Multiply(invdfErrordv, fError, correction);

	    vTest[0] -= correction[0];
	    vTest[1] -= correction[1];
	    vTest[2] -= correction[2];

            count++;
            forceError = DotProduct(fError, fError);

        } /* end Newton-Raphson while */

	delete [] dragGlide_B;
	delete [] dragGlide_D;
	delete [] dragGlide_V0;
	
        forceError = DotProduct(fError, fError);
/*
 *      If the mobility function failed to converge on a velocity,
 *      return an error.
 */
        if ((count > maxIterCount) || (forceError > fErrorTol)) {
	  PrintNode(node);
	  printf("[Warning] fail to converge velocity \n");
	  for (i = 0; i < numNbrs; i++) {
	    printf("i:%d angle1Exists:%d angle2Exists:%d angle3Exists:%d "\
		   "segAngle1Length:%.2f segAngl21Length:%.2f segAngle3Length:%.2f\n",
		   i, angle1Exists[i], angle2Exists[i],    angle3Exists[i],
		   segAngle1Length[i], segAngle2Length[i], segAngle3Length[i]);
	  }

	  return(1);
        }

/*
 *      We were able to converge on a velocity
 *
 *      If needed, rotate the velocity vector back to the laboratory frame
 *      from the crystal frame.
 */
        if (param->useLabFrame) {
           real8 vTestCryst[3];
           Matrix33Vector3Multiply(param->rotMatrix, vTest, vTestCryst);
	   VECTOR_COPY(vTest, vTestCryst);
        }
        
        node->vX = vTest[0];
        node->vY = vTest[1];
        node->vZ = vTest[2];

#ifdef NAN_CHECK
        if ((isnan(node->vX) || isinf(node->vX)) ||
            (isnan(node->vY) || isinf(node->vY)) ||
            (isnan(node->vZ) || isinf(node->vZ))) {
            PrintNode(node);
            Fatal("Mobility_FCC_linear: (%d, %d) velocity error.",
                  node->myTag.domainID, node->myTag.index);
        }
#endif

        return(0);
}
