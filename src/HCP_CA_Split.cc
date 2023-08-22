/*-------------------------------------------------------------------------
 *
 *      Module:       HCP_CA_Split.cc
 *
 *      Description:  This module contains functions needed for splitting
 *                    HCP dislocations with burgers vector types <c+a>
 *                    into two segments with burgers vectors of types
 *                    <c> and <a> respectively and cross-slip one of
 *                    the segments to maximize power dissipation.
 *
 *      Includes public functions:
 *          HCP_CA_Split()
 *
 *      Included private functions:
 *          CreateCandidateNode()
 *          SetNode()
 *
 *-----------------------------------------------------------------------*/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"

/*-------------------------------------------------------------------------
 *
 *      Function:       SetNode()
 *
 *      Description:    Set burgers vectors and plane for a discretization
 *                      node and its neighbors.
 *
 * FIX ME!  document args: change arg list to use Node_t * only?
 *
 *      Parameters:
 *         IN/OUT: nodePtr  Pointer to pointer of node to be modified.
 *         IN:     burg     burgers vector node's segments should be given
 *         IN:     plane    plane node's segments should be placed on
 *         IN/OUT: node1Ptr Pointer to pointer to nodePtr's 1st neighbor
 *         IN/OUT: node2Ptr Pointer to pointer to nodePtr's 2nd neighbor
 *         IN:     notifyRemoteDoms  Toggle indicating if the remote domains
 *                                   should be notified of any action taken
 *
 *-----------------------------------------------------------------------*/
static void SetNode(Home_t *home, Node_t **nodePtr, real8 burg[3],
                    real8 plane[3], Node_t **nbr1Ptr, Node_t **nbr2Ptr,
                    int notifyRemoteDoms)
{
        Node_t *node, *nbr1, *nbr2;

        node = *nodePtr;
        nbr1 = *nbr1Ptr;
        nbr2 = *nbr2Ptr;

        node->oldvX = 0.0;
        node->oldvY = 0.0;
        node->oldvZ = 0.0;

        node->vX = 0.0;
        node->vY = 0.0;
        node->vZ = 0.0;

        ChangeArmBurg(home, node, &nbr1->myTag,
                      burg[X], burg[Y], burg[Z], 
                      plane[X], plane[Y], plane[Z],
                      notifyRemoteDoms, DEL_SEG_NONE);

        ChangeArmBurg(home, nbr1, &node->myTag,
                      -burg[X], -burg[Y], -burg[Z], 
                      plane[X], plane[Y], plane[Z],
                      notifyRemoteDoms, DEL_SEG_NONE);

        ChangeArmBurg(home, node, &nbr2->myTag,
                      -burg[X], -burg[Y], -burg[Z], 
                      plane[X], plane[Y], plane[Z],
                      notifyRemoteDoms, DEL_SEG_NONE);

        ChangeArmBurg(home, nbr2, &node->myTag,
                      burg[X], burg[Y], burg[Z], 
                      plane[X], plane[Y], plane[Z],
                      notifyRemoteDoms, DEL_SEG_NONE);

        *nodePtr = node;
        *nbr1Ptr = nbr1;
        *nbr2Ptr = nbr2;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       CreateCandidateNode()
 *
 *      Description:    Creates a new node called <newNode> at the same
 *                      position as <origNode>.  Both <newNode> and <origNode>
 *                      are connected to <nbrNode1> and <nbrNode2>.
 *
 *                      Note: the burgers vectors and planes for the new
 *                            node/segments are NOT set.
 *
 *      Parameters:
 *         IN:     origNodePtr
 *         OUT:    newNodePtr 
 *         IN:     nbrNode1Ptr
 *         IN:     nbrNode2Ptr
 *
 *-----------------------------------------------------------------------*/
static void CreateCandidateNode(Home_t *home,
                                Node_t **origNodePtr, Node_t **newNodePtr,
                                Node_t **nbrNode1Ptr, Node_t **nbrNode2Ptr)
{
/*
 *      Need to create a new node at the original location
 *      without telling remote domains to do the same we do this.
 */
        Node_t *newNode, *origNode, *nbrNode1, *nbrNode2;
        real8   dir[3];

        origNode = *origNodePtr;
        nbrNode1 = *nbrNode1Ptr;
        nbrNode2 = *nbrNode2Ptr;
  
        VECTOR_ZERO(dir);

        newNode = GetNewNativeNode(home);
        FreeNodeArms(newNode);

        newNode->native = 1;
        SET_CONSTRAINTS(newNode->constraint, UNCONSTRAINED);

        newNode->x = origNode->x;
        newNode->y = origNode->y;
        newNode->z = origNode->z;

/* FIX ME!  Can we use the cell ID of the original node here??? */
        AssignNodeToCell(home, newNode);

/*
 *      To insert new segments between the newly created node and
 *      the two neighbor nodes...
 *
 *      NOTE: InsertArm() must be called twice for each segment to
 *            make the necessary changes in all affected nodes.
 */

        InsertArm(home, newNode, &nbrNode1->myTag,
                  dir[X], dir[Y], dir[Z], 
                  dir[X], dir[Y], dir[Z], 0);

        InsertArm(home, nbrNode1, &newNode->myTag,
                  -dir[X], -dir[Y], -dir[Z], 
                   dir[X],  dir[Y],  dir[Z], 0);

/*
 *      For the second segment, burgers vector is just opposite of
 *      burgers vector for first segment, and plane is the same.
 */
        InsertArm(home, newNode, &nbrNode2->myTag,
                  -dir[X], -dir[Y], -dir[Z], 
                   dir[X],  dir[Y],  dir[Z], 0);

        InsertArm(home, nbrNode2, &newNode->myTag,
                  dir[X], dir[Y], dir[Z], 
                  dir[X], dir[Y], dir[Z], 0);

        *newNodePtr  = newNode;
        *origNodePtr = origNode;
        *nbrNode1Ptr = nbrNode1;
        *nbrNode2Ptr = nbrNode2;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       HCP_CA_Split() 
 *
 *      Description:    Attempt to split HCP segments with burgers vectors
 *                      of type <c+a> into two with burgers vectors of types
 *                      <c> and <a> respectively and to cross-slip one of
 *                      the segments to maximize power dissipation.
 *
 *-----------------------------------------------------------------------*/
void HCP_CA_Split(Home_t *home)
{
        int     inode, i, n;
        int     notifyRemoteDoms = 1;
        int     doOperation, pinned1, pinned2;
        int     nodeInd, newnodeInd;
        int     numNodes, burgIndex, planen1Index, planen2Index;
        int     numBurg, *numGlissilePlanesPerBurg;
        int     numGlidePlanes, firstGlidePlaneIndex;
        int     firstGlidePlaneaIndex, firstGlidePlanecIndex;
        int     *burgFirstPlaneIndex, *planeType;
        int     burgaIndex, burgcIndex;
        int     nbDir, FoundNode;
        int     mobError, planeaIndex, planecIndex;
        real8   minSegSquared, minSeg2Squared;
        real8   eps = 1.0e-12, delta = 1.0e-16, nu;
        real8   (*burgList)[3], (*planeList)[3];
        real8   burgLab[3], burgCrystal[3], n1[3], n2[3];
        real8   glidedir[3];
        real8   areamin, areaflat;
        real8   powerMax, powerTest;
        real8   tmp, tmp3[3];
        real8   burga[3], burgc[3];
        real8   planenode[3],planenewnode[3],burgnode[3],burgnewnode[3];
        real8   planea[3], planec[3];
        real8   fdotglide, vecdotdir;
        real8   nodep[3], nbr1p[3], nbr2p[3], newnodep[3], nbr3p[3];
        real8   savenodep[3], savenbr1p[3], savenbr2p[3], savenewnodep[3];
        real8   vecP1P2[3], vecP1P[3], vecP2P[3], vecP3P[3];
        real8   vecQP1[3], vecQP2[3];
        real8   planeRefa[3], planeRefc[3];
        real8   newNodePos[3], nodePos[3], nbr1Pos[3], nbr2Pos[3];
     // real8   newNodeF[3], NodeF[3];
        real8   Nbr1F[3]={0.0, 0.0, 0.0};
        real8   Nbr2F[3]={0.0, 0.0, 0.0};
        real8   forcetmp[3];
        real8   newNodeVeloc[3];
        real8   vNoise, normvec;
        real8   burgLabtest[3][3], burgCrystaltest[3][3], narray[3];
        Node_t  *node, *nbr1, *nbr2, *nbr3, *newnode;
        Param_t *param;
        MobArgs_t  mobArgs;
        NodeState_t origNodeState, origNbr1State, origNbr2State, origNbr3State;

        param = home->param;

        nu =  param->pois;


/* 
 *      This routine is only for HCP materials 
 */
        if (param->materialType != MAT_TYPE_HCP) {
            return;
        }

/*
 *      If c+a burgers vector splitting has not been enabled via the
 *      control file parameter, or it is not supposed to be done during
 *      the current cycle, no need to do anything.
 */
        if (param->HCP_CA_SplitFrequency == 0) {
            return;
        }

        if ((home->cycle % param->HCP_CA_SplitFrequency) != 0) {
            return;
        }

        numNodes       =   home->newNodeKeyPtr;
        minSegSquared  =   param->minSeg * param->minSeg;
        minSeg2Squared = 4*param->minSeg * param->minSeg;
        vNoise         =   param->splitMultiNodeAlpha * param->rTol / param->deltaTT;

/* 
 *      Get lists of HCP Burgers vectors and glide planes
 */
        burgList = home->burgData.burgList;
        planeList = home->burgData.planeList;
        numBurg = home->burgData.numBurgVectors;
        numGlissilePlanesPerBurg = home->burgData.numGlissilePlanesPerBurg;
        burgFirstPlaneIndex = home->burgData.burgFirstPlaneIndex;
        planeType = home->burgData.planeType;

/*
 *      Area criteria
 */ 
        areamin   = param->remeshAreaMin;

/* 
 *      Loop over nodes
 */
        for (inode = 0; inode < numNodes; inode++) 
        {
            int   n, ma, mtot;
            int   start, end;
            int   nodeplaneIndex[7], newnodeplaneIndex[7];
            real8 nodeplanes[7][3], newnodeplanes[7][3];
            real8 nodeburgs[7][3], newnodeburgs[7][3];
            real8 tstar[7][3];
            real8 t[3], Cat[3], normalCat;
            real8 glideDirCrystal[4][3],glideDirLab[4][3];
            real8 bat[3], bct[3];

            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) {
                continue;
            }

/*
 *          Add some restrictions on the segments we'll consider.
 *
 *          First, start with discretization nodes that are neither
 *          pinned nor surface nodes.
 */
            if ((node->numNbrs != 2) || 
                HAS_ANY_OF_CONSTRAINTS(node->constraint,
                                       PINNED_NODE | SURFACE_NODE)) {
                continue;
            }

/*
 *          Don't mess with any nodes that are flagged to prevent collisions
 *          or other topological changes this cycle.
 */
            if ((node->flags & NO_COLLISIONS) ||
                (node->flags & NO_MESH_COARSEN))  {
                continue;
            }

/*
 *          As with all topological operations, we have to make sure
 *          that we abide by node/segment ownership rules in order
 *          to perform these operations in parallel.  However, since
 *          this function will potentially reposition the nodes
 *          involved in the split/cross-slip, the local domain MUST
 *          own all the nodes involved!
 */
            doOperation = 1;

            if ((node->nbrTag[0].domainID != home->myDomain) ||
                (node->nbrTag[1].domainID != home->myDomain)) {
                continue;
            }

/*
 *          If needed, rotate a copy of the burgers vector from the lab
 *          frame to the crystal frame.  Otherwise, lab and crystal frames
 *          are identical.
 */
            burgLab[X] = node->burgX[0];
            burgLab[Y] = node->burgY[0];
            burgLab[Z] = node->burgZ[0];
            VECTOR_COPY(burgCrystal, burgLab);

            if (param->useLabFrame) {
                Matrix33Vector3Multiply(param->rotMatrixInverse, burgLab,
                                        burgCrystal);
            }

/*
 *          Just a quick sanity check...
 */
            if ((fabs(burgCrystal[X]) < eps) &&
                (fabs(burgCrystal[Y]) < eps) &&
                (fabs(burgCrystal[Z]) < eps)) {
                continue;
            }

/* 
 *          Only <c+a> Burgers vectors are considered.
 *          Burgers vector of type <c+a> are only for burgIndex 3 -> 8. 
 *
 *          Find the index of the current burgers vector in the reference
 *          burgers vector list. Keep only the <c+a> Burgers vectors 
 *          (i.e. index [3, 8])
 */
            GetBurgIndexHCP(burgCrystal, 1, numBurg, burgList, &burgIndex);

            if ((burgIndex < 3) || (burgIndex > 8)) {
                continue;
            }

            numGlidePlanes  = numGlissilePlanesPerBurg[burgIndex];
            firstGlidePlaneIndex = burgFirstPlaneIndex[burgIndex];

/*
 *          Check that vectors (nbr1, node, nbr2) are aligned i.e. the
 *          area their triangle forms is flat.  If not, loop on to the
 *          next discretization node.
 */
            nbr1 = GetNeighborNode(home, node, 0);
            nbr2 = GetNeighborNode(home, node, 1);

            if ((nbr1 == (Node_t *)NULL) ||
                (nbr2 == (Node_t *)NULL)) {
                printf("WARNING: Neighbor not found at %s line %d\n",
                        __FILE__, __LINE__);
                continue;
            }

/*
 *          Get the positions of all three nodes and convert the neighbor
 *          node positions to the PBC coordinates relative to the main
 *          node.  These adjusted positions will be used during the
 *          cross slip algorithm, and if updated during the process, will
 *          be copied back into the respective nodal data structures and
 *          shifted back into the primary image space.
 */
            nodep[X] = node->x; nodep[Y] = node->y; nodep[Z] = node->z;
            nbr1p[X] = nbr1->x; nbr1p[Y] = nbr1->y; nbr1p[Z] = nbr1->z;
            nbr2p[X] = nbr2->x; nbr2p[Y] = nbr2->y; nbr2p[Z] = nbr2->z;

            PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
                        &nbr1p[X], &nbr1p[Y], &nbr1p[Z]);
            PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
                        &nbr2p[X], &nbr2p[Y], &nbr2p[Z]);

            vecP1P2[X] = nbr2p[X] - nbr1p[X];
            vecP1P2[Y] = nbr2p[Y] - nbr1p[Y];
            vecP1P2[Z] = nbr2p[Z] - nbr1p[Z];

/*
 *          If vector vecP1P2 is less than twice the desired minimum
 *          segment length, don't bother trying to split this node.
 */
            if (DotProduct(vecP1P2, vecP1P2) < minSeg2Squared) {
                continue;
            }

            vecP1P[X] = nodep[X] - nbr1p[X];
            vecP1P[Y] = nodep[Y] - nbr1p[Y];
            vecP1P[Z] = nodep[Z] - nbr1p[Z];

            vecP2P[X] = nodep[X] - nbr2p[X];
            vecP2P[Y] = nodep[Y] - nbr2p[Y];
            vecP2P[Z] = nodep[Z] - nbr2p[Z];

/*
 *          Only keep pair of segments that are aligned: see if area
 *          formed by the triangle of the three nodes (nbr1, node, nbr2)
 *          small. CHECK WHAT SMALL MEANS HERE.
 */
            cross(vecP1P,vecP2P,tmp3);
            areaflat = 0.5*sqrt(DotProduct(tmp3,tmp3));

            if (areaflat > 2*areamin) {
                continue;
            }

/*
 *          Area is flat. This vector should be considered for possible
 *          split and cross-slip operation. 
 */

            VECTOR_COPY(t,vecP1P2);
            NormalizeVec(t);
/*
 *          Find candidates for cross-slip
 *
 *          Determine Burgers vectors <a> and <c> from Burgers vector <c+a>.
 *          Doing following decomposition to get <a> and <c> ensures
 *          Burgers vector conservation.
 *
 *          CHECK: Assigning b_a and b_c in the right lab/crystal frame here.
 */
            burga[X] = burgLab[X];
            burga[Y] = burgLab[Y];
            burga[Z] = 0.0;

            burgc[X] = 0.0;
            burgc[Y] = 0.0;
            burgc[Z] = burgLab[Z];

/*
 *          Do the Kroupa criterion (5-17 H&L): look if the interaction force between
 *          the two inifinitely long dislocations defined by (t, ba) and (t, bc) 
 *          are attractive or repulsive. 
 *          If they are attractive, do not attempt a splitting
 *          Note : F12 < 0 means attractive.
 */

            cross(burga,t,bat);
            cross(burgc,t,bct);

            if ( (DotProduct(burga,t) * DotProduct(burgc,t) + DotProduct(bat,bct)/(1.0-nu)) < 0.0 ) 
            {
                // Interaction force between the two dislocations is attractive. 
                continue;
            }

/*
 *          Get indices of these two Burgers vectors in the list
 *          For <a>, it has to be one of the 3 first Burgers vectors.
 *          For <c>, it is the 10th Burgers vector.
 */
            GetBurgIndexHCP(burga, 1, 3, burgList, &burgaIndex);
            burgcIndex = 9;

/* 
 *          Find information to get list of glide planes for <a> and <c> directions 
 */
            firstGlidePlaneaIndex = burgFirstPlaneIndex[burgaIndex];
            firstGlidePlanecIndex = burgFirstPlaneIndex[burgcIndex];

/*
 *          Find the two glide planes for the <c+a> segments: node-nbr1 and
 *          node-nbr2.  If plane is not specified, define one.
 */
            n1[X] = node->nx[0];
            n1[Y] = node->ny[0];
            n1[Z] = node->nz[0];

            n2[X] = node->nx[1];
            n2[Y] = node->ny[1];
            n2[Z] = node->nz[1];

            if (Normal(n1)<eps && Normal(n2)<eps) {
                FindPreciseGlidePlane(home, burgLab, vecP1P, n1,
                                      param->allowFuzzyGlidePlanes);
                VECTOR_COPY(n2, n1);
            } else if (Normal(n1)<eps && Normal(n2)>eps) {
                VECTOR_COPY(n1, n2);
            } else if (Normal(n2)<eps && Normal(n1)>eps) {
                VECTOR_COPY(n2, n1);
            }

            FindGlidePlaneHCP(home, burgIndex, n1, &planen1Index);
            FindGlidePlaneHCP(home, burgIndex, n2, &planen2Index);


            /* TEMPORARY: suppose for now that n1 = n2 */
            if (planen2Index != planen1Index) {
                continue;
            }

/* 
 *          Check that either nbr1 or nbr2 can be moved before attempting a
 *          split. 
 *
 *          node - nbrs is of type : <c+a>. It should not have more than
 *          4 glide planes.
 */
            for (n = 0; n < numGlidePlanes; n++) {
                NormalizedCrossVector(burgList[burgIndex],
                                      planeList[firstGlidePlaneIndex+n],
                                      glideDirCrystal[n]);
            }

            if (param->useLabFrame) {
                int   row;
/*
 *              Rotate initial glide dir from crystal to lab frame
 */
                for (row = 0; row < numGlidePlanes; row++) {
                    Matrix33Vector3Multiply(param->rotMatrix,
                                            glideDirCrystal[row],
                                            glideDirLab[row]);
                }
            } else {
/*
 *              Lab and crystal frames are identical, so just copy
 *              the vectors
 */
                for (n = 0; n < numGlidePlanes; n++) {
                    VECTOR_COPY(glideDirLab[n], glideDirCrystal[n]);
                }
            }

            pinned1 = NodePinned(home, nbr1, planen1Index, glideDirLab,
                                 numGlidePlanes);

            pinned2 = NodePinned(home, nbr2, planen2Index, glideDirLab,
                                 numGlidePlanes);

/*
 *          Check if the neighboring nodes can be repositioned.  If neither
 *          can be moved, don't try the split and cross-slip.
 */
            if ((pinned1 == 1) && (pinned2 == 1)) {
                continue;
            }


/* 
 *          Find all possible directions for which the cross-slip can be done
 *          Since we used glide planes only, there cannot be more than 4 for <a>
 *          and 3 planes for <c> found during the search.
 */

/*
 *          Find directions for <a>
 */
            start = firstGlidePlaneaIndex;
            end   = firstGlidePlaneaIndex +
                    numGlissilePlanesPerBurg[burgaIndex];
            ma = 0;

            for (n = start; n < end; n++) {

                if (ma > 4) {
                    Fatal("HCP_CA_Split(): Found more than 4 glissile "
                          "planes for a <a> direction\n");
                }

                /* t_cat orthogonal to n1 and n_a : newnode will be along <a> */
                cross(planeList[n],n1,Cat);
                normalCat = Normal(Cat);

                if (normalCat > eps) {
                    Cat[X] /= normalCat;
                    Cat[Y] /= normalCat;
                    Cat[Z] /= normalCat;

                    /* t parallel to t_cat */
                    if (fabs(DotProduct(t,Cat)) > 0.9995) {
                        VECTOR_COPY(planea,planeList[n]);
                        FindPreciseGlidePlane(home, burgc, Cat, planec,
                                              param->allowFuzzyGlidePlanes);

/*
 *                      Consider this couple (n_a and n_c) only if these
 *                      planes are different.  Both planes should be
 *                      normalized.
 */
                        if (fabs(DotProduct(planea,planec)) < 0.99995) 
                        {
                            FindGlidePlaneHCP(home, burgaIndex, planea,
                                              &(newnodeplaneIndex[ma]));
                            FindGlidePlaneHCP(home, burgcIndex, planec,
                                              &(nodeplaneIndex[ma]));

                            VECTOR_COPY(newnodeburgs[ma],burga);
                            VECTOR_COPY(newnodeplanes[ma],planea);

                            VECTOR_COPY(   nodeburgs[ma],burgc);
                            VECTOR_COPY(   nodeplanes[ma],planec);
                            VECTOR_COPY(tstar[ma],Cat);

                            ma++;
                        }


                    }
                }
            }

/*
 *          Add on directions for <c>
 */
            start = firstGlidePlanecIndex;
            end   = firstGlidePlanecIndex +
                    numGlissilePlanesPerBurg[burgcIndex];

            mtot = ma;

            for (n = start; n < end; n++) {

                /* t_cat orthogonal to n1 and n_c : newnode will be along <c> */
                cross(planeList[n],n1,Cat);
                normalCat = Normal(Cat);

                if (normalCat > eps) {
                    Cat[X] /= normalCat;
                    Cat[Y] /= normalCat;
                    Cat[Z] /= normalCat;
                    if (fabs(DotProduct(t,Cat)) > 0.9995) {
                        VECTOR_COPY(planec,planeList[n]);
                        FindPreciseGlidePlane(home, burga, Cat, planea,
                                              param->allowFuzzyGlidePlanes);

/*
 *                      Consider this couple (n_a and n_c) only if these
 *                      planes are different
 */
                        if (fabs(DotProduct(planea,planec)) < 0.99995) {

/*
 *                          Quick sanity check...
 */
                            if (mtot > 6) {
                                Fatal("HCP_CA_Split(): Found more than 3 "
                                      "glissile planes for a <c> direction\n");
                            }

                            FindGlidePlaneHCP(home, burgcIndex, planec,
                                              &(newnodeplaneIndex[mtot]));
                            FindGlidePlaneHCP(home, burgaIndex, planea,
                                              &(nodeplaneIndex[mtot]));
                            VECTOR_COPY(newnodeburgs[mtot],burgc);
                            VECTOR_COPY(newnodeplanes[mtot],planec);

                            VECTOR_COPY(   nodeburgs[mtot],burga);
                            VECTOR_COPY(   nodeplanes[mtot],planea);
                            VECTOR_COPY(tstar[mtot],Cat);
                            mtot++;
                        }
                    }
                }
            }

/*
 *          No directions were found. Stop here
 */
            if (mtot == 0) {
                continue;
            }

/*
 *          Preserve the original state of the three nodes before modifying
 *          nodes' connections.
 */
            PreserveState(home, node, &origNodeState);
            PreserveState(home, nbr1, &origNbr1State);
            PreserveState(home, nbr2, &origNbr2State);

/* 
 *          Calculate the power of the initial configuration composed of
 *          nbr1, node and nbr2.
 */
            powerMax = ((node->fX * node->vX) +
                        (node->fY * node->vY) +
                        (node->fZ * node->vZ) +

                        (nbr1->fX * nbr1->vX) +
                        (nbr1->fY * nbr1->vY) +
                        (nbr1->fZ * nbr1->vZ) +

                        (nbr2->fX * nbr2->vX) +
                        (nbr2->fY * nbr2->vY) +
                        (nbr2->fZ * nbr2->vZ));

/* 
 *          Create a new node <newnode> on top of node <node>. 
 *          The call to CreateCandidateNode does not set any Burgers vector
 *          or plane, just creates the new node with segments to the 
 *          proper neighbors.
 */
            CreateCandidateNode(home, &node, &newnode, &nbr1, &nbr2);

            newnodep[X] = newnode->x;
            newnodep[Y] = newnode->y;
            newnodep[Z] = newnode->z;

            PBCPOSITION(param, newnodep[X], newnodep[Y], newnodep[Z],
                        &nbr1p[X], &nbr1p[Y], &nbr1p[Z]);
            PBCPOSITION(param, newnodep[X], newnodep[Y], newnodep[Z],
                        &nbr2p[X], &nbr2p[Y], &nbr2p[Z]);

/*              
 *          Try splitting the node and cross-slipping one of the segments
 */
            FoundNode = 0;

            for (nbDir=0; nbDir < mtot; nbDir++) {
                real8 r1, r2, r3, s;
                real8 dr1dt, dr2dt, dr3dt, dsdt;
                real8 darea2dt=0;
                real8 dvecQP1dt[3], dvecQP2dt[3], dvecP1P2dt[3];
/*
 *              Re-initialize node, newnode, nbr1 and nbr2's positions
 */
                VECTOR_TO_SCALAR(node->x, node->y, node->z, nodep);
                VECTOR_TO_SCALAR(newnode->x, newnode->y,
                                 newnode->z, newnodep);
                VECTOR_TO_SCALAR(nbr1->x, nbr1->y, nbr1->z, nbr1p);
                VECTOR_TO_SCALAR(nbr2->x, nbr2->y, nbr2->z, nbr2p);

/*
 *              Set Burgers vectors and glide planes of node, newnode,
 *              nbr1 and nbr2
 */
                SetNode(home, &node, nodeburgs[nbDir], nodeplanes[nbDir],
                        &nbr1, &nbr2, 0);
                SetNode(home, &newnode, newnodeburgs[nbDir],
                        newnodeplanes[nbDir], &nbr1, &nbr2,0);

/*
 *              Compute glide direction of the newnode.  Need its force
 *              and velocity. node and newnode are on top of each other.
 */
                SetOneNodeForce(home, newnode); 

                VECTOR_FROM_SCALAR(forcetmp, newnode->fX, newnode->fY,
                                   newnode->fZ);
                mobError = EvaluateMobility(home, newnode, &mobArgs);

                if (mobError) {
                    doOperation = 0;
                }


                VECTOR_FROM_SCALAR(newNodeVeloc, newnode->vX, newnode->vY,
                                   newnode->vZ);

                NormalizedCrossVector(tstar[nbDir], newnodeplanes[nbDir],
                                      glidedir);

/*
 *              Move newnode along glide direction.  Adjust all other nodes
 *              as well.
 */
                vecdotdir = DotProduct(vecP1P2,tstar[nbDir]);

                savenbr1p[X] = nbr1p[X];
                savenbr1p[Y] = nbr1p[Y];
                savenbr1p[Z] = nbr1p[Z];

                savenbr2p[X] = nbr2p[X];
                savenbr2p[Y] = nbr2p[Y];
                savenbr2p[Z] = nbr2p[Z];

                if (pinned1) {
/*
 *                  nbr2 can be moved (since we wouldn't reach this point
 *                  if both nbrs were pinned). Adjust its position so that
 *                  it is exactly aligned with tstar
 */
                    savenbr2p[X] = nbr1p[X] + (vecdotdir * tstar[nbDir][X]);
                    savenbr2p[Y] = nbr1p[Y] + (vecdotdir * tstar[nbDir][Y]);
                    savenbr2p[Z] = nbr1p[Z] + (vecdotdir * tstar[nbDir][Z]);

                } else {
/*
 *                  nbr1 can be moved. Adjust its position so that
 *                  it is exactly aligned with tstar
 */
                    savenbr1p[X] = nbr2p[X] - (vecdotdir * tstar[nbDir][X]);
                    savenbr1p[Y] = nbr2p[Y] - (vecdotdir * tstar[nbDir][Y]);
                    savenbr1p[Z] = nbr2p[Z] - (vecdotdir * tstar[nbDir][Z]);
                }

/*
 *              Now adjust position of node
 */
                vecdotdir = DotProduct(vecP1P, tstar[nbDir]);
                savenodep[X] = savenbr1p[X] + (vecdotdir * tstar[nbDir][X]);
                savenodep[Y] = savenbr1p[Y] + (vecdotdir * tstar[nbDir][Y]);
                savenodep[Z] = savenbr1p[Z] + (vecdotdir * tstar[nbDir][Z]);

/*
 *              Finally, move newnode
 */
                fdotglide = DotProduct(forcetmp, glidedir);

                vecP1P2[X] = savenbr2p[X] - savenbr1p[X];
                vecP1P2[Y] = savenbr2p[Y] - savenbr1p[Y];
                vecP1P2[Z] = savenbr2p[Z] - savenbr1p[Z];

                normvec=  Normal(vecP1P2);

                if (normvec < eps) {
                    doOperation = 0;
                    tmp = 0.0;
                } else {
                    /* 4.0 because there are two triangles : 2 x areamin */
                    tmp = areamin / normvec * 4.0 * (1.0 + eps) * Sign(fdotglide);
                }

                savenewnodep[X] = savenodep[X] + tmp * glidedir[X];
                savenewnodep[Y] = savenodep[Y] + tmp * glidedir[Y];
                savenewnodep[Z] = savenodep[Z] + tmp * glidedir[Z];

/*
 *              Test 1: The area of the triangle formed by nbr1, newnode,
 *              nbr2 should be growing.
 *
 *              Need forces and velocities of all nodes involved in the
 *              calculations. Update positions for accurate calculations.
 */
                VECTOR_TO_SCALAR(newnode->x, newnode->y, newnode->z,
                                 savenewnodep);
                VECTOR_TO_SCALAR(node->x, node->y, node->z, savenodep);
                VECTOR_TO_SCALAR(nbr1->x, nbr1->y, nbr1->z, savenbr1p);
                VECTOR_TO_SCALAR(nbr2->x, nbr2->y, nbr2->z, savenbr2p);

                SetOneNodeForce(home, newnode);
                SetOneNodeForce(home, node);
                SetOneNodeForce(home, nbr1);
                SetOneNodeForce(home, nbr2);

                mobError  = EvaluateMobility(home, newnode, &mobArgs);
                mobError |= EvaluateMobility(home, node, &mobArgs);
                mobError |= EvaluateMobility(home, nbr1, &mobArgs);
                mobError |= EvaluateMobility(home, nbr2, &mobArgs);

                if (mobError) {
                    doOperation = 0;
                }

                if (mobError == 0) {
/* 
 *                  Area for triangle (nbr1, nbr2, newnode) using Heron's formula 
 */
                    vecQP1[X] = savenbr1p[X] - savenewnodep[X]; 
                    vecQP1[Y] = savenbr1p[Y] - savenewnodep[Y];
                    vecQP1[Z] = savenbr1p[Z] - savenewnodep[Z];

                    vecQP2[X] = savenbr2p[X] - savenewnodep[X]; 
                    vecQP2[Y] = savenbr2p[Y] - savenewnodep[Y];
                    vecQP2[Z] = savenbr2p[Z] - savenewnodep[Z];

                    vecP1P2[X] = savenbr2p[X] - savenbr1p[X]; 
                    vecP1P2[Y] = savenbr2p[Y] - savenbr1p[Y];
                    vecP1P2[Z] = savenbr2p[Z] - savenbr1p[Z];

                    dvecQP1dt[X] = nbr1->vX - newnode->vX;
                    dvecQP1dt[Y] = nbr1->vY - newnode->vY;
                    dvecQP1dt[Z] = nbr1->vZ - newnode->vZ;

                    dvecQP2dt[X] = nbr2->vX - newnode->vX;
                    dvecQP2dt[Y] = nbr2->vY - newnode->vY;
                    dvecQP2dt[Z] = nbr2->vZ - newnode->vZ;

                    dvecP1P2dt[X] = nbr2->vX - nbr1->vX;
                    dvecP1P2dt[Y] = nbr2->vY - nbr1->vY;
                    dvecP1P2dt[Z] = nbr2->vZ - nbr1->vZ;

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

                    darea2dt = dsdt * (s-r1) * (s-r2) * (s-r3);
                    darea2dt += s * (dsdt-dr1dt) * (s-r2) * (s-r3);
                    darea2dt += s * (s-r1) * (dsdt-dr2dt) * (s-r3);
                    darea2dt += s * (s-r1) * (s-r2) * (dsdt-dr3dt);

                }  // end if (!mobError) 

/*
 *              If area is growing, test 1 has passed.
 *              Test 2: Power dissipation criterion
 */
                if (darea2dt > 0.0 && mobError == 0) {

                    /* Power dissipation criterion */
                    powerTest = 
                        (node->fX    * node->vX) +
                        (node->fY    * node->vY) +
                        (node->fZ    * node->vZ) +
                        (nbr1->fX    * nbr1->vX) +
                        (nbr1->fY    * nbr1->vY) +
                        (nbr1->fZ    * nbr1->vZ) +
                        (nbr2->fX    * nbr2->vX) +
                        (nbr2->fY    * nbr2->vY) +
                        (nbr2->fZ    * nbr2->vZ) +
                        (newnode->fX * newnode->vX) +
                        (newnode->fY * newnode->vY) +
                        (newnode->fZ * newnode->vZ) -
                        vNoise * (sqrt((node->fX    * node->fX)    +
                                       (node->fY    * node->fY)    +
                                       (node->fZ    * node->fZ))   +
                                  sqrt((nbr1->fX    * nbr1->fX)    +
                                       (nbr1->fY    * nbr1->fY)    +
                                       (nbr1->fZ    * nbr1->fZ))   +
                                  sqrt((nbr2->fX    * nbr2->fX)    +
                                       (nbr2->fY    * nbr2->fY)    +
                                       (nbr2->fZ    * nbr2->fZ))   +
                                  sqrt((newnode->fX * newnode->fX) +
                                       (newnode->fY * newnode->fY) +
                                       (newnode->fZ * newnode->fZ)));

 /*
  *                 If this  node split would result in the
  *                 highest energy release, save enough info to 
  *                 perform this split later.
  */
                    if ((powerTest - powerMax) > eps) 
                    {
                        VECTOR_FROM_SCALAR(nodePos   , savenodep   [X], savenodep   [Y], savenodep   [Z]);
                        VECTOR_FROM_SCALAR(newNodePos, savenewnodep[X], savenewnodep[Y], savenewnodep[Z]);
                        VECTOR_FROM_SCALAR(nbr1Pos   , savenbr1p   [X], savenbr1p   [Y], savenbr1p   [Z]);
                        VECTOR_FROM_SCALAR(nbr2Pos   , savenbr2p   [X], savenbr2p   [Y], savenbr2p   [Z]);

                        /* Save Burgers vectors and planes */
                        VECTOR_COPY(planenode   , nodeplanes   [nbDir]);
                        VECTOR_COPY(planenewnode, newnodeplanes[nbDir]);
                        VECTOR_COPY(burgnode    , nodeburgs    [nbDir]);
                        VECTOR_COPY(burgnewnode , newnodeburgs [nbDir]);

                        /* Save velocity of newNode */
                        VECTOR_FROM_SCALAR(newNodeVeloc, newnode->vX, newnode->vY, newnode->vZ);

                        powerMax = powerTest;
                        FoundNode = 1;

                        /* Save force of all nodes */
                     // VECTOR_FROM_SCALAR(newNodeF, newnode->fX, newnode->fY, newnode->fZ);
                     // VECTOR_FROM_SCALAR(NodeF, node->fX, node->fY, node->fZ);
                        VECTOR_FROM_SCALAR(Nbr1F, nbr1->fX, nbr1->fY, nbr1->fZ);
                        VECTOR_FROM_SCALAR(Nbr2F, nbr2->fX, nbr2->fY, nbr2->fZ);
                    }

                } // if area growth

            } // nbDir
    

/*
 *          If none of the configurqations tested yielded a higher energy
 *          release than the orginal configuration, do nothing.
 */
            if (FoundNode == 0) {
                doOperation = 0;
            }

/*
 *          Get rid of the temporary node and segments
 *          added during the evaluation and restore the original state...
 */
            ChangeArmBurg(home, newnode, &nbr1->myTag   , 0, 0, 0, 0, 0, 0, 0, DEL_SEG_NONE);
            ChangeArmBurg(home, nbr1   , &newnode->myTag, 0, 0, 0, 0, 0, 0, 0, DEL_SEG_NONE);
            ChangeArmBurg(home, newnode, &nbr2->myTag   , 0, 0, 0, 0, 0, 0, 0, DEL_SEG_NONE);
            ChangeArmBurg(home, nbr2   , &newnode->myTag, 0, 0, 0, 0, 0, 0, 0, DEL_SEG_NONE);

            RemoveNode(home, newnode, 0);
            newnode = (Node_t *)NULL;

            RestoreState(home, node, &origNodeState);
            RestoreState(home, nbr1, &origNbr1State);
            RestoreState(home, nbr2, &origNbr2State);

/*
 *          We've finished the evaluation and it has been decided
 *          to do the operation for real, the configuration has already
 *          been restored to its original state, so now we redo the
 *          selected operation and make sure the remote domains
 *          find out about it.
 */
            if (doOperation) {
                Node_t *newNode;

/*
 *              First reset the burgers vector for the original segments.
 *              <burgnode> and <planenode> are the burgers
 *              vector and glide plane from the original node to the first
 *              neighbor node.
 */
                SetNode(home, &node, burgnode, planenode,
                        &nbr1, &nbr2, notifyRemoteDoms);

                FoldBox(param, &(nodePos[X]),&(nodePos[Y]),&(nodePos[Z]));
                RepositionNode(home, nodePos, &node->myTag, notifyRemoteDoms);

                FoldBox(param, &(nbr1Pos[X]),&(nbr1Pos[Y]),&(nbr1Pos[Z]));
                RepositionNode(home, nbr1Pos, &nbr1->myTag, notifyRemoteDoms);

                FoldBox(param, &(nbr2Pos[X]),&(nbr2Pos[Y]),&(nbr2Pos[Z]));
                RepositionNode(home, nbr2Pos, &nbr2->myTag, notifyRemoteDoms);

/*
 *              Allocate the new node...
 */
                newNode = GetNewNativeNode(home);
                FreeNodeArms(newNode);

                newNode->native = 1;
                newNode->constraint = UNCONSTRAINED;

/*
 *              Reposition the new node (if necessary)
 *              Reinstate its force and velocity.
 */
                FoldBox(param, &(newNodePos[X]), &(newNodePos[Y]), &(newNodePos[Z]));

                RepositionNode(home, newNodePos, &newNode->myTag, notifyRemoteDoms);

                newNode->oldvX = 0.0;
                newNode->oldvY = 0.0;
                newNode->oldvZ = 0.0;

                newNode->vX = newNodeVeloc[X];
                newNode->vY = newNodeVeloc[Y];
                newNode->vZ = newNodeVeloc[Z];

                AssignNodeToCell(home, newNode);

/*
 *              Have the remote domains add the new node in as a 
 *              ghost node by having them do a "split_node" operation.
 */
                AddOpSplitNode(home, &node->myTag, &newNode->myTag, newNodePos, newNodeVeloc);

                InsertArm(home, newNode, &nbr1->myTag,
                          burgnewnode[X], burgnewnode[Y], burgnewnode[Z], 
                          planenewnode[X], planenewnode[Y], planenewnode[Z], 1);

                InsertArm(home, nbr1, &newNode->myTag,
                          -burgnewnode[X], -burgnewnode[Y], -burgnewnode[Z], 
                          planenewnode[X], planenewnode[Y], planenewnode[Z], 1);

                InsertArm(home, newNode, &nbr2->myTag,
                          -burgnewnode[X], -burgnewnode[Y], -burgnewnode[Z], 
                          planenewnode[X], planenewnode[Y], planenewnode[Z], 1);

                InsertArm(home, nbr2, &newNode->myTag,
                          burgnewnode[X], burgnewnode[Y], burgnewnode[Z], 
                          planenewnode[X], planenewnode[Y], planenewnode[Z], 1);

/*
 *              Also mark all the nodes involved so we don't perform any
 *              collisions with these nodes for the remainder of this
 *              timestep. nbr1 and nbr2 are OK. No not flag them.
 */
                node->flags    |= NO_COLLISIONS;
                newNode->flags |= NO_COLLISIONS;

/*
 *              Mark the force and velocity data for some nodes as obsolete so
 *              that more accurate forces will be recalculated at the beginning
 *              of the next timestep. Forces of nbr1 and nbr2 are OK we've recalculated them.
 */
                MarkNodeForceObsolete(home, node);
                MarkNodeForceObsolete(home, newNode);

                /* Get force information for nbr1 and nbr2 */
                nbr1->fX = Nbr1F[X];
                nbr1->fY = Nbr1F[Y];
                nbr1->fZ = Nbr1F[Z];

                nbr2->fX = Nbr2F[X];
                nbr2->fY = Nbr2F[Y];
                nbr2->fZ = Nbr2F[Z];
            }  /* end if (doOperation) */

/*
 *          Before we end the loop iteration for this node, we need to
 *          discard any state information we may have saved for the
 *          original node.
 */
            DiscardState(&origNodeState);
            DiscardState(&origNbr1State);
            DiscardState(&origNbr2State);

        }  /* end for (inode = 0; inode < home->newNodeKeyPtr ...) */

/*
 *      ZIPPER CASE
 *
 *      Loop through all the nodes again, this time looking
 *      for 3-node to handle 'zipper' type events.
 */
        for (inode = 0; inode < numNodes; inode++) {
            int     m = 0;
            int     glidea = 0, glidec = 0;
            int     bSegIDcpa = -1, bSegIDc = -1, bSegIDa = -1;
            real8   glideDirCrystal[4][3], glideDirLab[4][3];
            real8   t[3], Cat[3], normalCat, tstarbis[3];
            Node_t *nbrPtr1 = (Node_t *)NULL, *nbrPtr2 = (Node_t *)NULL;

            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) {
                continue;
            }

/*
 *          Node has three neighbors
 */
            if ((node->numNbrs != 3) ||  
                HAS_ANY_OF_CONSTRAINTS(node->constraint,
                                       PINNED_NODE | SURFACE_NODE)) {
                continue;
            }

/*
 *          Don't look at any nodes that are flagged to prevent collisions
 *          or other topological changes this cycle.
 */
            if ((node->flags & NO_COLLISIONS) ||
                (node->flags & NO_MESH_COARSEN)) {
                continue;
            }


/*
 *          As with all topological operations, we have to make sure
 *          that we abide by node/segment ownership rules in order
 *          to perform these operations in parallel.  However, since
 *          this function will potentially reposition the nodes
 *          involved in the split/cross-slip, the local domain MUST
 *          own all the nodes involved!
 */
            doOperation = 1;

            if ((node->nbrTag[0].domainID != home->myDomain) ||
                (node->nbrTag[1].domainID != home->myDomain) ||
                (node->nbrTag[2].domainID != home->myDomain)) {
                continue;
            }

/*
 *          If needed, rotate a copy of the burgers vector from the lab
 *          frame to the crystal frame.  Otherwise, lab and crystal frames
 */
            burgLabtest[0][X] = node->burgX[0];
            burgLabtest[0][Y] = node->burgY[0];
            burgLabtest[0][Z] = node->burgZ[0];
            VECTOR_COPY(burgCrystaltest[0], burgLabtest[0]);

            if (param->useLabFrame) {
                Matrix33Vector3Multiply(param->rotMatrixInverse, burgLabtest[0],
                                        burgCrystaltest[0]);
            }

            burgLabtest[1][X] = node->burgX[1];
            burgLabtest[1][Y] = node->burgY[1];
            burgLabtest[1][Z] = node->burgZ[1];
            VECTOR_COPY(burgCrystaltest[1], burgLabtest[1]);

            if (param->useLabFrame) {
                Matrix33Vector3Multiply(param->rotMatrixInverse, burgLabtest[1],
                                        burgCrystaltest[1]);
            }

            burgLabtest[2][X] = node->burgX[2];
            burgLabtest[2][Y] = node->burgY[2];
            burgLabtest[2][Z] = node->burgZ[2];
            VECTOR_COPY(burgCrystaltest[2], burgLabtest[2]);

            if (param->useLabFrame) {
                Matrix33Vector3Multiply(param->rotMatrixInverse, burgLabtest[2],
                                        burgCrystaltest[2]);
            }

            narray[0] = Normal(burgCrystaltest[0]);
            narray[1] = Normal(burgCrystaltest[1]);
            narray[2] = Normal(burgCrystaltest[2]);

/*
 *         We are only interested in the following case : 
 *          - one neighbor is of type <c+a>
 *          - one neighbor is of type <c>
 *          - one neighbor is of type <a>
 */

/*
 *          Check if any of the segments has a <c+a> burgers vector.
 *          If none was found, loop to the next node.
 */
            FindAbsMax(narray, 3, &tmp, &bSegIDcpa);

            if (bSegIDcpa == -1) {
                continue;
            }

/*
 *          Found a node, check that is is a <c+a> and get its index.
 *          Get indices of <a> and <c> associated to <c+a>.
 */

            VECTOR_COPY(burgLab, burgLabtest[bSegIDcpa]);
            VECTOR_COPY(burgCrystal, burgCrystaltest[bSegIDcpa]);

            GetBurgIndexHCP(burgCrystal, 1, numBurg, burgList, &burgIndex);

            if (burgIndex < 3 || burgIndex > 8) {
                continue;
            }

            numGlidePlanes  = numGlissilePlanesPerBurg[burgIndex];
            firstGlidePlaneIndex = burgFirstPlaneIndex[burgIndex];

            burga[X] = burgLab[X];
            burga[Y] = burgLab[Y];
            burga[Z] = 0.0;

            burgc[X] = 0.0;
            burgc[Y] = 0.0;
            burgc[Z] = burgLab[Z];

/*
 *          Get indices of these two Burgers vectors in the list
 *          For <a>, it has to be one of the 3 first Burgers vectors.
 *          For <c>, it is the 9th Burgers vector.
 */

            GetBurgIndexHCP(burga, 1, 3, burgList, &burgaIndex);
            burgcIndex = 9;

/* 
 *          Find information to get list of planes for <a> and <c> directions 
 */

            firstGlidePlaneaIndex = burgFirstPlaneIndex[burgaIndex];
            firstGlidePlanecIndex = burgFirstPlaneIndex[burgcIndex];

/*
 *          Determine which segment has the <c> burgers vector
 *          and which has the <a>
 */
            for (i = 0; i < 3; i++) {
                if (i != bSegIDcpa) {
                    /* test for <c> */
                    if ((fabs(burgLabtest[i][X]) < eps) &&
                        (fabs(burgLabtest[i][Y]) < eps) ) {
                        bSegIDc = i;
                    } else {
                        bSegIDa = i;
                    }
                }
            }

/*
 *          If we don't have segments with all of <a>, <c> and <c+a>
 *          burgers vectors, loop to the next node.
 */
            if (bSegIDcpa ==-1 || bSegIDa ==-1 || bSegIDc==-1) {
                continue;
            }

/*
 *          Set the neighboring node pointers now:
 *
 *              First neighbor is <c+a>
 *              Second neighbor is <a>
 *              Third neighbor is <c>
 */
            nbr1 = GetNeighborNode(home, node, bSegIDcpa);
            nbr2 = GetNeighborNode(home, node, bSegIDa);
            nbr3 = GetNeighborNode(home, node, bSegIDc);

            PreserveState(home, node, &origNodeState);
            PreserveState(home, nbr1, &origNbr1State);
            PreserveState(home, nbr2, &origNbr2State);
            PreserveState(home, nbr3, &origNbr3State);

            nodep[X] = node->x; nodep[Y] = node->y; nodep[Z] = node->z;
            nbr1p[X] = nbr1->x; nbr1p[Y] = nbr1->y; nbr1p[Z] = nbr1->z;
            nbr2p[X] = nbr2->x; nbr2p[Y] = nbr2->y; nbr2p[Z] = nbr2->z;
            nbr3p[X] = nbr3->x; nbr3p[Y] = nbr3->y; nbr3p[Z] = nbr3->z;

            PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
                        &nbr1p[X], &nbr1p[Y], &nbr1p[Z]);
            PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
                        &nbr2p[X], &nbr2p[Y], &nbr2p[Z]);
            PBCPOSITION(param, nodep[X], nodep[Y], nodep[Z],
                        &nbr3p[X], &nbr3p[Y], &nbr3p[Z]);

            vecP1P[X] = nodep[X] - nbr1p[X];
            vecP1P[Y] = nodep[Y] - nbr1p[Y];
            vecP1P[Z] = nodep[Z] - nbr1p[Z];

/*
 *          If the vector vecP1P is less than the minimum desired segment
 *          length, don't bother trying to do anything with this node
 */
            if (DotProduct(vecP1P, vecP1P) < minSegSquared) {
                continue;
            }

            vecP2P[X] = nodep[X] - nbr2p[X];
            vecP2P[Y] = nodep[Y] - nbr2p[Y];
            vecP2P[Z] = nodep[Z] - nbr2p[Z];

            vecP3P[X] = nodep[X] - nbr3p[X];
            vecP3P[Y] = nodep[Y] - nbr3p[Y];
            vecP3P[Z] = nodep[Z] - nbr3p[Z];

/*
 *          Find the glide plane, if any, for this <c+a> segment
 *          If plane is not specified, define one.
 */
            n1[X] = node->nx[bSegIDcpa];
            n1[Y] = node->ny[bSegIDcpa];
            n1[Z] = node->nz[bSegIDcpa];

            if (Normal(n1) < eps) {
                FindPreciseGlidePlane(home, burgLab, vecP1P, n1,
                                      param->allowFuzzyGlidePlanes);
            }

            FindGlidePlaneHCP(home, burgIndex, n1, &planen1Index);

/* 
 *          Check that the node nbr1 can be moved before tempting a split
 */
            for (n = 0; n < numGlidePlanes; n++) {
                NormalizedCrossVector(burgList[burgIndex],
                                      planeList[firstGlidePlaneIndex+n],
                                      glideDirCrystal[n]);
            }

                  
            if (param->useLabFrame) {
                int   row;
/*
 *              Rotate initial glide dir from crystal to lab frame
 */
                for (row = 0; row < numGlidePlanes; row++) {
                    Matrix33Vector3Multiply(param->rotMatrix,
                                            glideDirCrystal[row],
                                            glideDirLab[row]);
                }
            } else {
/*
 *              Lab and crystal frames are identical, so just copy
 *              the vectors
 */
                for (n = 0; n < numGlidePlanes; n++) {
                    VECTOR_COPY(glideDirLab[n], glideDirCrystal[n]);
                }
            }

/*
 *          Check is nbr1 will be able to move
 */
            pinned1 = NodePinned(home, nbr1, planen1Index, glideDirLab,
                                 numGlidePlanes);

/*
 *          Need the planes for <a> and <c> too 
 */
            planeRefa[X] = node->nx[bSegIDa];
            planeRefa[Y] = node->ny[bSegIDa];
            planeRefa[Z] = node->nz[bSegIDa];

            if (Normal(planeRefa)<eps) {
                FindPreciseGlidePlane(home, burga, vecP2P, planeRefa,
                                      param->allowFuzzyGlidePlanes);
            }

            FindGlidePlaneHCP(home, burgaIndex, planeRefa, &planeaIndex);

            planeRefc[X] = node->nx[bSegIDc];
            planeRefc[Y] = node->ny[bSegIDc];
            planeRefc[Z] = node->nz[bSegIDc];

            if (Normal(planeRefc)<eps) {
                FindPreciseGlidePlane(home, burgc, vecP3P, planeRefc,
                                      param->allowFuzzyGlidePlanes);
            }

            FindGlidePlaneHCP(home, burgcIndex, planeRefc, &planecIndex);

/*
 *          Check that the planes <a> and <c> are different
 */
            if (fabs(DotProduct(planeRefa,planeRefc)) > 0.99995) {
                continue;
            }

/*
 *          Can move nodes on glide planes only.
 *          Check that the planes are glide planes
 */
    
            if (planeType[planeaIndex] != 4) { glidea = 1; }
            if (planeType[planecIndex] != 4) { glidec = 1; }

/*
 *          The possible directions for which the cross-plane can be done is either one
 *          or the other direction. Will move the glissile one. If both are glissile,
 *          any can move but that is symmetrical.
 */
            VECTOR_COPY(t,vecP1P);
            NormalizeVec(t);

            /* t_cat orthogonal to n1 and n_a : newnode will be along <a> */
            cross(planeRefa,n1,Cat);
            normalCat = Normal(Cat);

            if (normalCat > eps && glidea==1) {

                Cat[X] /= normalCat;
                Cat[Y] /= normalCat;
                Cat[Z] /= normalCat;

                if (fabs(DotProduct(t,Cat)) > 0.9995) {

                    /* node is <c> */
                    nodeInd = bSegIDc;
                    VECTOR_COPY(planenode, planeRefc);
                    VECTOR_COPY(burgnode, burgc);

                    /* newnode is <a> */
                    newnodeInd = bSegIDa;
                    VECTOR_COPY(planenewnode, planeRefa);
                    VECTOR_COPY(burgnewnode, burga);

                    VECTOR_COPY(tstarbis,Cat);

                    m++;
                }
            }

            /* t_cat orthogonal to n1 and n_c : newnode will be along <c> */
            cross(planeRefc,n1,Cat);  
            normalCat = Normal(Cat);

            if (normalCat > eps && glidec==1) {

                Cat[X] /= normalCat;
                Cat[Y] /= normalCat;
                Cat[Z] /= normalCat;

                if (fabs(DotProduct(t,Cat)) > 0.9995) {

                    /* node is <a> */
                    nodeInd = bSegIDa;

                    VECTOR_COPY(planenode,planeRefa);
                    VECTOR_COPY(burgnode,burga);

                    /* newnode is <c> */
                    newnodeInd = bSegIDc;

                    VECTOR_COPY(planenewnode,planeRefc);
                    VECTOR_COPY(burgnewnode,burgc);

                    VECTOR_COPY(tstarbis,Cat);
                    m++;
                }
            }

            if (m == 0) {
                doOperation = 0;
                continue;
            }

            nbrPtr1 = GetNeighborNode(home, node, newnodeInd);
            nbrPtr2 = GetNeighborNode(home, node, nodeInd);

/* 
 *          Power of initial configuration composed of P, P1, P2 and P3
 */
            powerMax = ((node->fX * node->vX) +
                        (node->fY * node->vY) +
                        (node->fZ * node->vZ) +

                        (nbr1->fX * nbr1->vX) +
                        (nbr1->fY * nbr1->vY) +
                        (nbr1->fZ * nbr1->vZ) +

                        (nbrPtr1->fX * nbrPtr1->vX) +
                        (nbrPtr1->fY * nbrPtr1->vY) +
                        (nbrPtr1->fZ * nbrPtr1->vZ) +

                        (nbrPtr2->fX * nbrPtr2->vX) +
                        (nbrPtr2->fY * nbrPtr2->vY) +
                        (nbrPtr2->fZ * nbrPtr2->vZ));



/*
 *          Create a new node <newnode> on top of node <node>
 *          Q is connected to nbr1 (<c+a>) and nbrPtr1 (<a> or <c>).
 */
            CreateCandidateNode(home, &node, &newnode, &nbr1, &nbrPtr1);
            VECTOR_COPY(newnodep,nodep);

/*
 *          Change <node> connections : node should not be connected to
 *          nbrPtr1 anymore.  node is connected to nbr1 (<c+a>) and nbrPtr2
 *          (<a> or <c>).
 */
            ChangeArmBurg(home, node, &nbrPtr1->myTag, 0, 0, 0, 0, 0, 0, 0, DEL_SEG_NONE);
            ChangeArmBurg(home, nbrPtr1, &node->myTag, 0, 0, 0, 0, 0, 0, 0, DEL_SEG_NONE);

            SetNode(home, &node, burgnode, planenode, &nbr1, &nbrPtr2, 0);
            SetNode(home, &newnode, burgnewnode, planenewnode, &nbr1, &nbrPtr1, 0);

/*
 *          Compute glide direction of the newnode.  Need its force
 *          and velocity.
 */
            SetOneNodeForce(home, newnode);
            VECTOR_FROM_SCALAR(forcetmp, newnode->fX, newnode->fY, newnode->fZ);
            NormalizedCrossVector(tstarbis, planenewnode, glidedir);

            vecdotdir = DotProduct(vecP1P,tstarbis);

/*
 *          nbr1 does not have to be moved. The alignment is done on node.
 *          Align node's position.
 */
            vecdotdir = DotProduct(vecP1P,tstarbis);
            nodep[X] = nbr1p[X] + (vecdotdir * tstarbis[X]);
            nodep[Y] = nbr1p[Y] + (vecdotdir * tstarbis[Y]);
            nodep[Z] = nbr1p[Z] + (vecdotdir * tstarbis[Z]);
/*
 *          Move newnode
 */
            fdotglide = DotProduct(forcetmp,glidedir);

            vecP1P[X] = nodep[X] - nbr1p[X];
            vecP1P[Y] = nodep[Y] - nbr1p[Y];
            vecP1P[Z] = nodep[Z] - nbr1p[Z];

            normvec = Normal(vecP1P);

            if (normvec < eps) {
                doOperation = 0;
                tmp = 0.0;
            } else {
                tmp = areamin / normvec * 2.0 * (1.0 + eps) * Sign(fdotglide);
            }

            newnodep[X] = nodep[X] + tmp * glidedir[X];
            newnodep[Y] = nodep[Y] + tmp * glidedir[Y];
            newnodep[Z] = nodep[Z] + tmp * glidedir[Z];

/*
 *          Compute forces and velocities for new configuration 
 *          Check power dissipation.
 *          nbrPtr1 and nbrPtr2 have not moved. Their neighbors 
 *          have changed
 */
            node->x = nodep[X];
            node->y = nodep[Y];
            node->z = nodep[Z];

            newnode->x = newnodep[X];
            newnode->y = newnodep[Y];
            newnode->z = newnodep[Z];

            nbr1->x = nbr1p[X];
            nbr1->y = nbr1p[Y];
            nbr1->z = nbr1p[Z];

            SetOneNodeForce(home, newnode); 
            SetOneNodeForce(home, node); 
            SetOneNodeForce(home, nbr1); 
            SetOneNodeForce(home, nbrPtr1); 
            SetOneNodeForce(home, nbrPtr2); 

            mobError  = EvaluateMobility(home, newnode, &mobArgs);
            mobError |= EvaluateMobility(home, node, &mobArgs);
            mobError |= EvaluateMobility(home, nbr1, &mobArgs);
            mobError |= EvaluateMobility(home, nbrPtr1, &mobArgs);
            mobError |= EvaluateMobility(home, nbrPtr2, &mobArgs);

            if (mobError) {
                doOperation = 0;
            }

            if (mobError == 0) {
/* 
 *              Power dissipation of the new configuration compared to the 
 *              power of initial configuration
 */
                powerTest = 
                    (node->fX    * node->vX) +
                    (node->fY    * node->vY) +
                    (node->fZ    * node->vZ) +
                    (nbr1->fX    * nbr1->vX) +
                    (nbr1->fY    * nbr1->vY) +
                    (nbr1->fZ    * nbr1->vZ) +
                    (nbrPtr1->fX    * nbrPtr1->vX) +
                    (nbrPtr1->fY    * nbrPtr1->vY) +
                    (nbrPtr1->fZ    * nbrPtr1->vZ) +
                    (nbrPtr2->fX    * nbrPtr2->vX) +
                    (nbrPtr2->fY    * nbrPtr2->vY) +
                    (nbrPtr2->fZ    * nbrPtr2->vZ) +
                    (newnode->fX * newnode->vX) +
                    (newnode->fY * newnode->vY) +
                    (newnode->fZ * newnode->vZ) -
                    vNoise * (sqrt((node->fX    * node->fX)  +
                                   (node->fY    * node->fY)  +
                                   (node->fZ    * node->fZ)) +
                              sqrt((nbr1->fX    * nbr1->fX)  +
                                   (nbr1->fY    * nbr1->fY)  +
                                   (nbr1->fZ    * nbr1->fZ)) +
                              sqrt((nbrPtr1->fX    * nbrPtr1->fX)  +
                                   (nbrPtr1->fY    * nbrPtr1->fY)  +
                                   (nbrPtr1->fZ    * nbrPtr1->fZ)) +
                              sqrt((nbrPtr2->fX    * nbrPtr2->fX)  +
                                   (nbrPtr2->fY    * nbrPtr2->fY)  +
                                   (nbrPtr2->fZ    * nbrPtr2->fZ)) +
                              sqrt((newnode->fX * newnode->fX) +
                                   (newnode->fY * newnode->fY) +
                                   (newnode->fZ * newnode->fZ)));

/*
 *              If this configuration has a higher power dissipation, save
 *              velocity of the newnode
 */
                if ((powerTest - powerMax) <= eps) {
                    doOperation = 0;
                } else {
                    VECTOR_FROM_SCALAR(newNodeVeloc, newnode->vX, newnode->vY,
                                       newnode->vZ);
                }

            }  // end if (!mobError)

/*
 *          Restore the connections of the original node
 */
            InsertArm(home, node, &nbrPtr1->myTag,
                      burgnode[X],burgnode[Y],burgnode[Z], 
                      planenode[X], planenode[Y], planenode[Z], 0);

            InsertArm(home, nbrPtr1, &node->myTag,
                      -burgnode[X],-burgnode[Y],-burgnode[Z], 
                       planenode[X],  planenode[Y],  planenode[Z], 0);
/*
 *          Get rid of the temporary node and segments added during the
 *          evaluation and restore the original state...
 */
            ChangeArmBurg(home, newnode, &nbr1->myTag   , 0, 0, 0, 0, 0, 0, 0, DEL_SEG_NONE);
            ChangeArmBurg(home, nbr1   , &newnode->myTag, 0, 0, 0, 0, 0, 0, 0, DEL_SEG_NONE);
            ChangeArmBurg(home, newnode, &nbrPtr1->myTag, 0, 0, 0, 0, 0, 0, 0, DEL_SEG_NONE);
            ChangeArmBurg(home, nbrPtr1, &newnode->myTag, 0, 0, 0, 0, 0, 0, 0, DEL_SEG_NONE);

            RemoveNode(home, newnode, 0);
            newnode = (Node_t *)NULL;

/*
 *          Need to restore node's connections....
 */
            RestoreState(home, node, &origNodeState);
            RestoreState(home, nbr1, &origNbr1State);
            RestoreState(home, nbr2, &origNbr2State);
            RestoreState(home, nbr3, &origNbr3State);

            nbrPtr1 = (Node_t *)NULL;
            nbrPtr2 = (Node_t *)NULL;

/*
 *          We've finished the evaluation and if it has been decided
 *          to do the operation for real, the configuration has already
 *          been restored to its original state, just redo the
 *          operation and make sure the remote domains find out about it.
 */
            if (doOperation) {
                Node_t *newNode;

                nbrPtr1 = GetNeighborNode(home, node, newnodeInd);
                nbrPtr2 = GetNeighborNode(home, node, nodeInd);

/*
 *              Remove the connection between node and nbrPtr1
 */
                ChangeArmBurg(home, node, &nbrPtr1->myTag,
                              0, 0, 0, 0, 0, 0, notifyRemoteDoms, 0);
                ChangeArmBurg(home, nbrPtr1, &node->myTag,
                              0, 0, 0, 0, 0, 0, notifyRemoteDoms, 0);

/*
 *              Change burgers and plane of node
 */
                SetNode(home, &node, burgnode, planenode, &nbr1, &nbrPtr2,
                        notifyRemoteDoms);

/*
 *              Allocate the new node...
 */
                newNode = GetNewNativeNode(home);
                FreeNodeArms(newNode);

                newNode->native = 1;
                SET_CONSTRAINTS(newNode->constraint, UNCONSTRAINED);

                newNode->x = newnodep[X];
                newNode->y = newnodep[Y];
                newNode->z = newnodep[Z];

                newNode->oldvX = 0.0;
                newNode->oldvY = 0.0;
                newNode->oldvZ = 0.0;

                newNode->vX = newNodeVeloc[X];
                newNode->vY = newNodeVeloc[Y];
                newNode->vZ = newNodeVeloc[Z];

                AssignNodeToCell(home, newNode);

/*
 *              Have the remote domains add the new node in as a 
 *              ghost node by having them do a "split_node" operation.
 */
                AddOpSplitNode(home, &node->myTag, &newNode->myTag,
                               newnodep, newNodeVeloc);

                InsertArm(home, newNode, &nbr1->myTag,
                          burgnewnode[X], burgnewnode[Y], burgnewnode[Z], 
                          planenewnode[X], planenewnode[Y], planenewnode[Z], 1);

                InsertArm(home, nbr1, &newNode->myTag,
                          -burgnewnode[X], -burgnewnode[Y], -burgnewnode[Z], 
                          planenewnode[X], planenewnode[Y], planenewnode[Z], 1);

                InsertArm(home, newNode, &nbrPtr1->myTag,
                          -burgnewnode[X], -burgnewnode[Y], -burgnewnode[Z], 
                          planenewnode[X], planenewnode[Y], planenewnode[Z], 1);

                InsertArm(home, nbrPtr1, &newNode->myTag,
                          burgnewnode[X], burgnewnode[Y], burgnewnode[Z], 
                          planenewnode[X], planenewnode[Y], planenewnode[Z], 1);

/*
 *            Also mark all the nodes involved so we don't perform any
 *            collisions with these nodes for the remainder of this
 *            timestep.
 */
                node->flags    |= NO_COLLISIONS;
                nbr1->flags    |= NO_COLLISIONS;
                nbr2->flags    |= NO_COLLISIONS;
                nbr3->flags    |= NO_COLLISIONS;
                newNode->flags |= NO_COLLISIONS;

/*
 *              Mark the force and velocity data for some nodes as obsolete so
 *              that more accurate forces will be recalculated at the beginning
 *              of the next timestep.
 */
                MarkNodeForceObsolete(home, node);
                MarkNodeForceObsolete(home, nbr1);
                MarkNodeForceObsolete(home, nbr2);
                MarkNodeForceObsolete(home, nbr3);
                MarkNodeForceObsolete(home, newNode);

            }  /* end if (doOperation) */

/*
 *          Before we end the loop iteration for this node, we need to
 *          discard any state information we may have saved for the
 *          original node.
 */
            DiscardState(&origNodeState);
            DiscardState(&origNbr1State);
            DiscardState(&origNbr2State);
            DiscardState(&origNbr3State);


        }  /* end for (inode = 0; inode < home->newNodeKeyPtr ...) */

        return;
}
