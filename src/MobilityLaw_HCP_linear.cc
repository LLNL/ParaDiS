/**************************************************************************
 *
 *      Module:       Mobility_HCP_linear.c
 *      Author:       Sylvie Aubry
 *      Description:  Contains functions for calculating generic isotropic
 *                    linear mobility for HCP materials in a non-linear
 *                    fashion.
 *
 *                    NOTE: It is more energetically favorable for a screw
 *                    dislocation to move on the basal slip plane rather
 *                    than staying on the pyramidal slip plane. A dislocation
 *                    cannot extend on the pyramidal slip plane in a direction
 *                    perpendicular to the Burgers vector
 *
 *      ==============================================================
 *
 *      The HCP linear mobility law computes the velocity of a node
 *      solution of
 *
 *          F^{elastic} (v) + F^{drag}(v) = 0.
 *
 *      The elastic forces F^{elastic} (v) at the node are given as
 *      input. The function F^{drag}(v) is constructed using molecular
 *      dynamics input data. A Newton-Raphson procedure is used to
 *      calculate v solution of that equation. The mobility law for
 *      HCP is linear, this means that this system can be re-written
 *      as B . v = F and by inverting the 3x3 matrix B, v could be
 *      found directly.
 *
 *      It is assumed that the user has drag coefficients for the
 *      the HCP material. The drag coefficients for the different
 *      HCP planes should be provided for edge and screw components.
 *
 *      For each dislocation segments connected to the node, we compute:
 *
 *          F^{drag}(v) = 1/2 Ls Fs(v) + 1/2 Le1 Fe1(v) + 1/2 Le2 Fe2(v)
 *
 *      where Fs, Fe1 and Fe2 are the drag coefficients for screw and edge 
 *      directions and Ls, Le1 and Le2 the decomposition of the segment
 *      into screw and edge components.  or if the segment is a junction:
 *
 *          F^{drag}(v) = 1/2 L Fj(v)
 *
 *      where L is the total length of the segment and Fj is the drag 
 *      coefficient for junction.
 *
 *      The decomposition of Fs(v) supposes that the velocity v is
 *      decomposed into
 *
 *          v = (v. b/|b|) b/|b| + (v.n) n + (v.m) m
 *
 *      where m = (b x n)
 *           n is the glide plane of the segment considered
 *
 *      then, since the mobility law is linear
 *
 *         Fs(v) = Bsline * (v. b/|b|) b/|b| +
 *                 Bclimb * (v.n) n +
 *                 Bglide * (v.m) m
 *
 *      Bglide is the drag coefficient corresponding to the glide plane n,
 *      the dislocation segment is on. For instance, if n is a basal
 *      plane, Bsglide is Bbasalscrew.
 *
 *      The decomposition of Fe(v) supposed the velocity is decomposed into
 *
 *          v = (v. t) t + (v.n) n + (v. b/|b|) b/|b|
 *
 *      where t = (n x b/|b|)
 *            n is one of the two highest projections of the edge
 *              components for that segment
 *
 *      If Beglide is the corresponding drag coefficient for the edge
 *      component, then
 *
 *          Fe(v) = Beline * (v. t) t +
 *                  Beclimb * (v.n) n +
 *                  Bglide * (v. b/|b|) b/|b| 
 *
 *      Includes public functions:
 *                    MobilityLaw_HCP_linear()
 *
 *      Includes private functions:
 *                    DecomposeArm()
 *                
 ***************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

//#define M33_INVERSE(a,b) Matrix33_Inverse(a,b)
//#define M33_INVERSE(a,b) Matrix33_SVD_Inverse(a,b)
  #define M33_INVERSE(a,b) Matrix33_PseudoInverse(a,b)

/**************************************************************************
 *
 *      Function:     DecomposeArm
 *      Description:  This function decomposes a dislocation segment into
 *                    screw and edges components.  There are two edge
 *                    components which are the highest projections of the
 *                    dislocation segment onto the different edge directions.
 *                    This function assumes the first edge component index
 *                    is known. The second one is computed within this function.
 *
 *                    Decomposition of a segment is given by:
 *
 *                        seg = fabs(screwLength * burgNorm) +
 *                              fabs(edge1Length * refPlanes[edge1Index]) +
 *                              fabs(edge2Length * refPlanes[edge2Index])
 *      Arguments:
 *      
 *          IN:   seg        Dislocation segment vector specified in the
 *                           crystal frame.
 *
 *          IN:   numEdgePlanes Number of reference planes associated
 *                           with the burgers vector.
 *                          
 *          IN:   burgNorm   Normalized burgers vector specified in the
 *                           crystal frame.
 *
 *          IN:   refEdges   Edge reference planes for this burgers vector
 *                           specified in the crystal frame.
 *
 *          OUT:  screwLen   Length of the segment along the screw direction
 *
 *          OUT:  edge1Len   Length of the segment along the first edge
 *
 *          OUT:  edge2Len   Length of the segment along the second edge
 *
 *          IN:   edge1Index Index in <refPlanes> of the edge direction
 *                           with the highest projection
 *
 *          OUT:  edge1Index Index in <refPlanes> of the edge direction
 *                           with the 2nd highest projection
 *
 *************************************************************************/
static void DecomposeArm(real8 seg[3], int numEdgePlanes,
                         real8 burgNorm[3],
                         real8 refEdges[4][3],
                         real8 *screwLength, real8 *edge1Length,
                         real8 *edge2Length, int edge1Index, int *edge2Index)
{
        real8 segLeft[3], tempLength[4];

/*
 *      Calculate screw length
 */
        *screwLength = fabs(DotProduct(seg, burgNorm));

/*
 *      Now get the lengths for both edges
 */
        segLeft[0] = seg[0] - (*screwLength) * burgNorm[0];
        segLeft[1] = seg[1] - (*screwLength) * burgNorm[1];
        segLeft[2] = seg[2] - (*screwLength) * burgNorm[2];

        for (int i=0; i<numEdgePlanes; i++)
            tempLength[i] = fabs(DotProduct(refEdges[i], segLeft));

/*
 *      We know the index of the highest edge value, now find the index
 *      of the second highest.
 */
        *edge2Index = Get2ndMaxIndex(tempLength, numEdgePlanes, edge1Index);

        DecompVec(segLeft, refEdges[edge1Index], refEdges[*edge2Index], tempLength);

        *edge1Length = fabs(tempLength[0]);
        *edge2Length = fabs(tempLength[1]);

        return;
}


/**************************************************************************
 *
 *      Function:     MobilityLaw_HCP_linear
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Arguments:
 *          IN:  node                Pointer to the node for which to
 *                                   calculate velocity
 *          OUT: invdfErrordv        Array in which to return the inverse of
 *                                   the drag matrix used in the function.
 *          OUT: glideConstraints    Array in which to return the normals
 *                                   to the glide planes to which this node
 *                                   is restricted.
 *          OUT: numGlideConstraints Number of glide plane normals returned
 *
 *      Returns: 0 if successful
 *               1 if function could not converge on a velocity
 *
 *************************************************************************/
int Mobility_HCP_linear(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int   i, j, m, n;
        int   count, maxIterCount,numEdgePlanes;
        int   numNbrs, typePlane1, typePlane2; 
        int   numNonZeroLenSegs = 0;
        int   firstPlaneIndex, numPlanes;
        int  *numPlanesPerBurg, *numGlissilePlanesPerBurg;
        int  *burgFirstPlaneIndex, *planeType;
        int   edgeExists[MAX_NBRS];
        int   screwExists[MAX_NBRS];
        int   junctionExists[MAX_NBRS];
        real8 mu, maxSegLen;
        real8 dragJunctionLine, dragJunctionGlide, dragLine, dragClimb;
        real8 massMult=0.0, massMatrix[3][3], eps, maxVal;
        real8 vTest[3], vOld[3], fn[3], rt[3];
        real8 (*burgList)[3], (*planeList)[3];
        real8 segTotLength[MAX_NBRS],segScrewLength[MAX_NBRS];
        real8 segEdge1Length[MAX_NBRS], segEdge2Length[MAX_NBRS];
        real8 dragGlideScrew[MAX_NBRS], dragGlideEdge1[MAX_NBRS];
        real8 dragGlideEdge2[MAX_NBRS];
        real8 burg[MAX_NBRS][3], ndir[MAX_NBRS][3], lineDir[MAX_NBRS][3];
        real8 edgeDir1[MAX_NBRS][3], edgeDir2[MAX_NBRS][3];
        real8 edgeLength;
        real8 AdragScrew[5], AdragEdge[5];
        real8 CdragScrew[2], CdragEdge[2];
        real8 CpAdragScrew[4], CpAdragEdge[4];
        real8 *planeListXplane;
        real8 fError[3], dfErrordv[3][3];
        real8 fErrorTol;
        real8 forceError, totLen;
        real8 correction[3], inertia[3];
        Param_t *param;
        Node_t  *nbrNode;

/*
 *      If any of the arms of the node has a burgers vector that
 *      has explicitly been set to be sessile (via control file
 *      inputs), the node may not be moved.
 */
        if (NodeHasSessileBurg(home, node))
        {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      If the node is pinned in ALL dimensions it cannot be moved
 *      so just zero the velocity and return
 */
        if (HAS_ALL_OF_CONSTRAINTS(node->constraint, PINNED_NODE))
        {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

        param = home->param;

        mu = param->shearModulus;
        maxSegLen = param->maxSeg;

        real8 (*invdfErrordv)[3]     =  mobArgs->invDragMatrix;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invdfErrordv);
        *numGlideConstraints = 0;

        eps = 1.0e-12;

/* 
 *      Parameters for junctions 
 */
        dragJunctionLine  = param->HCP_LineDrag;
        dragJunctionGlide = param->HCP_Sessile_ScrewDrag;

/*
 *      Parameter for line mobility: fast mobility, small drag
 */
        dragLine = param->HCP_LineDrag;

/* 
 *      Parameter for climb mobility: slow mobility, large drag
 */
        dragClimb = param->HCP_Sessile_ScrewDrag;

/*     
 *      Set up some arrays with the various drag coefficients from the
 *      control parameters in order to get the proper plane type
 *      later.
 */
        AdragScrew[0] = param->HCP_A_Basal_ScrewDrag;
        AdragScrew[1] = param->HCP_A_Prismatic_ScrewDrag;
        AdragScrew[2] = param->HCP_A_1stPyramidal_ScrewDrag;
        AdragScrew[3] = param->HCP_A_2ndPyramidal_ScrewDrag;
        AdragScrew[4] = param->HCP_Sessile_ScrewDrag;

        CdragScrew[0] = param->HCP_C_Prismatic_ScrewDrag;
        CdragScrew[1] = param->HCP_Sessile_ScrewDrag;

        CpAdragScrew[0] = param->HCP_CpA_Prismatic_ScrewDrag;
        CpAdragScrew[1] = param->HCP_CpA_1stPyramidal_ScrewDrag;
        CpAdragScrew[2] = param->HCP_CpA_2ndPyramidal_ScrewDrag;
        CpAdragScrew[3] = param->HCP_Sessile_ScrewDrag;

        AdragEdge[0] = param->HCP_A_Basal_EdgeDrag;
        AdragEdge[1] = param->HCP_A_Prismatic_EdgeDrag;
        AdragEdge[2] = param->HCP_A_1stPyramidal_EdgeDrag;
        AdragEdge[3] = param->HCP_A_2ndPyramidal_EdgeDrag;
        AdragEdge[4] = param->HCP_Sessile_EdgeDrag;

        CdragEdge[0] = param->HCP_C_Prismatic_EdgeDrag;
        CdragEdge[1] = param->HCP_Sessile_EdgeDrag;

        CpAdragEdge[0] = param->HCP_CpA_Prismatic_EdgeDrag;
        CpAdragEdge[1] = param->HCP_CpA_1stPyramidal_EdgeDrag;
        CpAdragEdge[2] = param->HCP_CpA_2ndPyramidal_EdgeDrag;
        CpAdragEdge[3] = param->HCP_Sessile_EdgeDrag;

/*
 *      If we need to include inertial terms, need some more initializations
 */
        if (param->includeInertia != 0)
            massMult = 0.25 * param->massDensity * (param->burgMag * param->burgMag);

        numNbrs = node->numNbrs;

/*
 *      This function saves various pieces of info for each neighbor
 *      node using static arrays for speed rather than doing memory
 *      allocations.  If the arrays are too  small, abort with an error
 *      to let the caller increase the array sizes.
 */
        if (numNbrs > MAX_NBRS)
        {
            Fatal("Mobility_HCP_linear: Node segment count (%d) exceeds\n"
                  "maximum allowed (%d) in static arrays.  Increase MAX_NBRS\n"
                  "and recompile\n%s:%d", node->numNbrs, MAX_NBRS, __FILE__, __LINE__);
        }

/*
 *      Set some pointers to the reference burgers vector and glide plane
 *      lists that were created during initialization.
 */

        burgList = home->burgData.burgList;
        planeList = home->burgData.planeList;
        numPlanesPerBurg = home->burgData.numPlanesPerBurg;
        numGlissilePlanesPerBurg = home->burgData.numGlissilePlanesPerBurg;
        burgFirstPlaneIndex = home->burgData.burgFirstPlaneIndex;
        planeType = home->burgData.planeType;

/*
 *      Loop over all arms attached to the node to get some information
 *      for the Newton-Raphson procedure
 */
        totLen = 0.0;

        for (i=0; i<numNbrs; i++) 
        {
            real8 rtSquared;
            real8 plane[3];

            nbrNode = GetNeighborNode(home, node, i);

            if (nbrNode == (Node_t *)NULL)
            {
                printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                continue;
            }

/*
 *          Get the segment line direction
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
            if ((rtSquared = DotProduct(rt, rt)) < eps)
            {
                segTotLength[i] = 0.0;
                continue;
            }

            numNonZeroLenSegs++;

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

/*
 *          If necessary, rotate the burgers vector and line sense from the
 *          laboratory frame to the crystal frame
 */
            if (param->useLabFrame) 
            {
                real8 burgRot[3], rtRot[3], planeRot[3];

                Matrix33Vector3Multiply(param->rotMatrixInverse, burg[i], burgRot );
                Matrix33Vector3Multiply(param->rotMatrixInverse, rt     , rtRot   );
                Matrix33Vector3Multiply(param->rotMatrixInverse, plane  , planeRot);

                VECTOR_COPY(burg[i], burgRot);
                VECTOR_COPY(rt, rtRot);
                VECTOR_COPY(plane, planeRot);
            }

/*
 *          Calculate the segment total length
 */
            segTotLength[i] = sqrt(rtSquared);

            totLen += segTotLength[i];

/* 
 *          Get the index in the reference Burgers vectors list
 *          of current Burgers vector.
 */
            int bIndex = 0;
            int pIndex = 0;

            GetBurgIndexHCP(burg[i], 0, HCP_NUM_GLIDE_BURG , burgList, &bIndex);

/*
 *          If the Burgers vector is not a junction, find the 
 *          index of the reference plane.
 */
            if (bIndex >= 0)
            {
                numPlanes  = numPlanesPerBurg[bIndex];
                planeListXplane = (real8 *)calloc(1, numPlanes * sizeof(real8));
                
                firstPlaneIndex = burgFirstPlaneIndex[bIndex];
                
                MatrixMult(&planeList[firstPlaneIndex][0],
                           numPlanes, 3, 3,
                           plane    , 1, 1,
                           planeListXplane, 1);
                
                FindAbsMax(planeListXplane, numPlanes, &maxVal, &pIndex);
                
                free(planeListXplane);
      
                pIndex += burgFirstPlaneIndex[bIndex];
                typePlane1 = planeType[pIndex];
            }
            else
            {
                /* set type of plane to be sessile if bIndex<0. This should always be true.*/
                pIndex = -1;
                typePlane1 = 4;
            }

/* 
 *          Type plane is one of these numbers:
 *             0 if plane is basal
 *             1 if plane is prismatic
 *             2 if plane is first pyramidal
 *             3 if plane is second pyramidal
 *             4 if plane is sessile          
 */ 

/*
 *          A dislocation segment is a junction if both its 
 *          Burgers vector and its planes are pre-defined as 
 *          junction Burgers vector and sessile plane.
 *          A bIndex > 0 can have a sessile plane too.
 */
            if (bIndex < 0) 
            {
/*
 *              Junction Burgers vectors
 */
                junctionExists[i] = 1;
                edgeExists    [i] = 0;
                screwExists   [i] = 0;
            }
            else if (typePlane1==4) 
            {
/*
 *              Sessile plane : junction
 */
                junctionExists[i] = 1;
                edgeExists    [i] = 0;
                screwExists   [i] = 0;
            } 
            else 
            {
/*
 *              Not a junction
 */
                real8 burgRef[3];
                real8 edgeRef[4][3];
                real8 planeRef[4][3];

                junctionExists[i] = 0;

/*
 *              Find the drag coefficient associated with the glide plane
 *              The drag coefficient depends on the Burgers vector type too.
 * 
 *              The following assumes that we know where in the pre-defined
 *              Burgers table the <a>, <c> and <c+a> burgers vector are.
 *              If the pre-defined table of Burgers vector changes,
 *              the code below must change also. 
 */

                if (bIndex < 3) 
                {
                    // <a> type 
                    dragGlideScrew[i] = AdragScrew[typePlane1];
                } 
                else if (bIndex == 9) 
                {
                    // <c> type : Prismatic or sessile.
                    if (typePlane1 == 1) { dragGlideScrew[i] = CdragScrew[0]; } 
                    else                 { dragGlideScrew[i] = CdragScrew[1]; }
                } else 
                {
                    // <c+a> type : all the rest
                    dragGlideScrew[i] = CpAdragScrew[typePlane1-1];
                }
/* 
 *              Get the corresponding planes and edges directions.
 *              Can look among glide planes only. There should not
 *              be more than 4 glide planes.
 */
                int edgeIndex1=0;
                int edgeIndex2=0;
                int start = burgFirstPlaneIndex[bIndex];
                int end   = burgFirstPlaneIndex[bIndex] + numGlissilePlanesPerBurg[bIndex];

                VECTOR_COPY(burgRef, burgList[bIndex]);
                NormalizeVec(burgRef);

                for (m=0, j=start; (j<end); j++, m++) 
                {
                    VECTOR_COPY(planeRef[m], planeList[j]);

                    cross(planeRef[m], burgRef, edgeRef[m]);
                    NormalizeVec(edgeRef[m]);

/* 
 *                  The planeref corresponding to pIndex gives the 
 *                  plane for the first edge.  Set edgeIndex1.
 */
                    if (j==pIndex) { edgeIndex1 = (j-start); }
                }

/*
 *              The second edge direction is found
 *              in the group of glissile planes only. 
 */

                numEdgePlanes = numGlissilePlanesPerBurg[bIndex];

/*
 *              Decompose segment into screw and edges components.
 *              Compute length of screw and edges as well as their 
 *              directions.
 */
                DecomposeArm(rt, numEdgePlanes, burgRef,  edgeRef,
                             &(segScrewLength[i]), &(segEdge1Length[i]),
                             &(segEdge2Length[i]), edgeIndex1, &edgeIndex2);

/*
 *              Find the drag associated with the edge's plane
 *              The drag coefficient depends on the Burgers vector type too.
 */
                typePlane1 = planeType[start+edgeIndex1];
                typePlane2 = planeType[start+edgeIndex2];

                if (bIndex < 3) 
                {
                    // <a> type 
                    dragGlideEdge1[i] = AdragEdge[typePlane1];
                    dragGlideEdge2[i] = AdragEdge[typePlane2];
                } 
                else if (bIndex == 9) 
                {
                    // <c> type : Prismatic or sessile
                    if (typePlane1==1) { dragGlideEdge1[i] = CdragEdge[0]; } 
                    else               { dragGlideEdge1[i] = CdragEdge[1]; }

                    if (typePlane2==1) { dragGlideEdge2[i] = CdragEdge[0]; } 
                    else               { dragGlideEdge2[i] = CdragEdge[1]; }

                } 
                else 
                {
                    // <c+a> type : all the rest
                    dragGlideEdge1[i] = CpAdragEdge[typePlane1-1];
                    dragGlideEdge2[i] = CpAdragEdge[typePlane2-1];
                }

                VECTOR_COPY(edgeDir1[i], edgeRef[edgeIndex1]);
                VECTOR_COPY(edgeDir2[i], edgeRef[edgeIndex2]);

                edgeLength = sqrt((segEdge1Length[i] * segEdge1Length[i]) +
                                  (segEdge2Length[i] * segEdge2Length[i]));

                edgeExists [i] = ((edgeLength       /segTotLength[i])>1.0e-04);
                screwExists[i] = ((segScrewLength[i]/segTotLength[i])>1.0e-04);

            } /* end glide Burgers vector */

/*
 *          Calculate tangents and climb directions
 */
            lineDir[i][0] = rt[0] / segTotLength[i];
            lineDir[i][1] = rt[1] / segTotLength[i];
            lineDir[i][2] = rt[2] / segTotLength[i];

            VECTOR_COPY(ndir[i], plane);

        }  /* loop over neighbors */

/*
 *      It's possible this function was called for a node which only
 *      had zero length segments (during SplitSurfaceNodes() for example).
 *      If that is the case, just set the velocity to zero and return;
 */
        if (numNonZeroLenSegs == 0)
        {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      Set up the mass matrix: if not needed, zero it out.
 */
        for (m=0; m<3; m++)
        {
            for (n=0; n<3; n++) 
            {
                massMatrix[m][n] = 0.0;

                if ((m == n) && (param->includeInertia))
                    massMatrix[m][n] = totLen * massMult / param->deltaTT;
            }
        }

/* 
 *      Setup for the Newton-Raphson loop
 *
 *      NOTE: using this linear mobility, Newton-Raphson should converge
 *      in one iteration
 */
        count = 0;
        maxIterCount = 3; 

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
        if (param->useLabFrame)
        {
            real8 fnRot[3], vTestRot[3];
            Matrix33Vector3Multiply(param->rotMatrixInverse, fn, fnRot);
            Matrix33Vector3Multiply(param->rotMatrixInverse, vTest, vTestRot);
            VECTOR_COPY(fn, fnRot);
            VECTOR_COPY(vTest, vTestRot);
            VECTOR_COPY(vOld, vTest);
        }

/*
 *      If we need to include inertial terms...
 */
        if (param->includeInertia)
        {
            real8 vDiff[3];
            vDiff[0] = vTest[0] - vOld[0];
            vDiff[1] = vTest[1] - vOld[1];
            vDiff[2] = vTest[2] - vOld[2];
            Matrix33Vector3Multiply(massMatrix, vDiff, inertia);
        }
        else 
        {
            VECTOR_ZERO(inertia);
        }

        forceError = fErrorTol + 1;

/*
 *      Begin the Newton-Raphson loop
 */
        while ((count < maxIterCount) && (forceError > fErrorTol))
        {
            for (m=0; (m<3); m++)
                fError[m] = fn[m] - inertia[m];
    
            for (m=0; (m<3); m++)
                for (n=0; (n<3); n++)
                    dfErrordv[m][n] = -massMatrix[m][n];

/*
 *          Loop over all segments (skipping zero length segments)
 *          adding up screw, edge and junction contributions to fError
 *          and dfErrordv.
 */
            for (i=0; i<numNbrs; i++)
            {
                real8 fJunction[3], dfJunctiondv[3][3];
                real8 fScrew   [3], dfScrewdv   [3][3];
                real8 fEdge    [3], dfEdgedv    [3][3];

/*
 *              Skip zero length segments
 */
                if (fabs(segTotLength[i]) < eps) { continue; }
/*
 *              Junction contribution
 */
                if (junctionExists[i])
                {
                    JunctionDrag(vTest, lineDir[i], dragJunctionLine,
                                 dragJunctionGlide, fJunction, dfJunctiondv);

                    for (m=0; m<3; m++) 
                    {

                        fError[m] -= 0.5 * segTotLength[i] * fJunction[m];

                        for (n=0; (n<3); n++)
                            dfErrordv[m][n] -= 0.5 * segTotLength[i] * dfJunctiondv[m][n];
                    }

                } 
                else /* not a junction */
                {  
/*
 *                  Screw contribution
 */
                    if (screwExists[i])
                    {
                        ScrewDrag2(vTest, burg[i], ndir[i],
                                  dragClimb, dragLine, dragGlideScrew[i],
                                  fScrew, dfScrewdv);

                        for (m=0; m<3; m++)
                        {
                            fError[m] -= 0.5 * segScrewLength[i] * fScrew[m];

                            for (n=0; (n<3); n++)
                                dfErrordv[m][n] -= 0.5 * segScrewLength[i] * dfScrewdv[m][n];
                        }
                    }

/*
 *                  Edge contribution
 */
                    if (edgeExists[i])
                    {
/*
 *                      First edge
 */
                        EdgeDrag(vTest, burg[i], edgeDir1[i],
                                 dragClimb, dragLine, dragGlideEdge1[i],
                                 fEdge, dfEdgedv);

                        for (m=0; m<3; m++)
                        {
                            fError[m] -= 0.5 * segEdge1Length[i] * fEdge[m];

                                for (n=0; (n<3); n++)
                                {
                                    dfErrordv[m][n] -= 0.5 * segEdge1Length[i] * dfEdgedv[m][n];
                                }
                        }
/*
 *                      Second edge
 */
                        EdgeDrag(vTest, burg[i], edgeDir2[i],
                                 dragClimb, dragLine, dragGlideEdge2[i],
                                 fEdge, dfEdgedv);

                        for (m=0; m<3; m++)
                        {
                            fError[m] -= 0.5 * segEdge2Length[i] * fEdge[m];

                            for (n=0; (n<3); n++)
                                dfErrordv[m][n] -= 0.5 * segEdge2Length[i] * dfEdgedv[m][n];
                        }

                    }  /* end if (edgeExists[i]) */

                } /* fError and dfErrordv are computed now */

            } /* end loop over segments */

            if ( M33_INVERSE(invdfErrordv,dfErrordv) < 0 )
            {
                Print3x3("dfErrordv", dfErrordv);
                Fatal("%s::%s(%d) : Cannot invert dfErrordv!", __FILE__, __func__, __LINE__ );
            }

            Matrix33Vector3Multiply(invdfErrordv, fError, correction);

            vTest[0] -= correction[0];
            vTest[1] -= correction[1];
            vTest[2] -= correction[2];

            count++;
            forceError = DotProduct(fError, fError);

        } /* end Newton-Raphson while */


        forceError = DotProduct(fError, fError);

/*
 *      If the mobility function failed to converge on a velocity,
 *      return an error.
 */
        if ((count > maxIterCount) || (forceError > fErrorTol))
            return(1);

/*
 *      We were able to converge on a velocity
 *
 *      If needed, rotate the velocity vector back to the laboratory frame
 *      from the crystal frame.
 */
        if (param->useLabFrame)
        {
            real8 vTestRot[3];
            Matrix33Vector3Multiply(param->rotMatrix, vTest, vTestRot);
            VECTOR_COPY(vTest, vTestRot);
        }

        node->vX = vTest[0];
        node->vY = vTest[1];
        node->vZ = vTest[2];

#ifdef NAN_CHECK
        if ( (isnan(node->vX) || isinf(node->vX)) ||
             (isnan(node->vY) || isinf(node->vY)) ||
             (isnan(node->vZ) || isinf(node->vZ)))
        {
            PrintNode(node);
            Fatal("Mobility_HCP_linear: (%d, %d) velocity error.", node->myTag.domainID, node->myTag.index);
        }
#endif

        return(0);
}
