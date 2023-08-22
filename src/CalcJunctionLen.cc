/***************************************************************************
 *
 *      Module:       CalcJunctionLen.c
 *      Description:  This module contains functions needed in SplitMultiNode
 *                    to compute the optimal junction length based on line
 *                    tension approximations. The junction length is the
 *                    solution of the equation: sum of selforces at the
 *                    node equilibrates the force along the junction direction.
 *
 *      Includes public functions:
 *          CalcJunctionLen()
 *          AdjustJuncNodeForceAndVel()
 *
 *      Includes private functions:
 *          CalcSegLineTension()
 *          CalcNodeLineTension()
 *          VelocityWithoutConstraint()
 *
 **************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

#ifdef ANISOTROPIC
#include "Anisotropic.h"
#endif

/* Wrapper iso/aniso core force and derivative */

static void CoreForceandDerivative
(
   Home_t *home, real8 x1, real8 y1, real8 z1,
   real8 x2, real8 y2, real8 z2,
   real8 b[3], real8 Ecore, 
   real8 f2[3], real8 dfdx2[3][3]
)
{
#ifndef ANISOTROPIC
  real8 NU = home->param->pois;
  CoreForceandDerivativeIsotropic(home, x1, y1, z1, x2, y2, z2, b, Ecore, NU, f2, dfdx2);
#else
  int qMax = home->param->anisoHarmonicsNumTermsBase;
  CoreForceandDerivativeAnisotropic(home, x1, y1, z1, x2, y2, z2, b, Ecore, qMax, f2, dfdx2);
#endif
}

static void CoreForce
(
   Home_t *home, real8 x1, real8 y1, real8 z1,
   real8 x2, real8 y2, real8 z2,
   real8 b[3], real8 Ecore, real8 f2[3]
)
{
#ifndef ANISOTROPIC
  real8 NU = home->param->pois;
  CoreForceIsotropic(home, x1, y1, z1, x2, y2, z2, b, Ecore, NU, f2);
#else
  int qMax = home->param->anisoHarmonicsNumTermsBase;
  CoreForceAnisotropic(home, x1, y1, z1, x2, y2, z2, b, Ecore, qMax, f2);
#endif
}

/*-------------------------------------------------------------------------
 *
 *      Function:     CalcNodeLineTension
 *      Description:  Calculate the sum of forces on a node coming from
 *                    the subset of its arms defined in armList.  These
 *                    forces are line tension approximations or 
 *                    self-force + PK forces.
 *
 *      Arguments:
 *          position    IN: point along vector <dir> from node's original 
 *                          position at which the forces are to be evaluated
 *          node        IN: pointer to node for which forces are being
 *                          estimated
 *          armCount    IN: number of node's arms to be included in evaluation
 *          armList     IN: list (indices) of node's arms to be included in
 *                          evaluation
 *          A           IN: Arbitrary distant point along the junction
 *                          direction from node position.
 *          bA          IN: sum of the burgers vector at node
 *          vErr        OUT: vector error measuring how far from equilibrium
 *                           the sum of forces is
 *          dvErrDX     OUT: vector derivative error for the Newton-Raphson
 *                           method.
 *
 *------------------------------------------------------------------------*/
static void CalcNodeLineTension(Home_t *home, real8 position[3], Node_t *node,
                                int armCount, int *armList, int calcForceDeriv,
                                real8 A[3], real8 bA[3], real8 verr[3],
                                real8 dverrdx[3][3])
{
        int i, j;
        int     arm;
        real8   Ecore;
        Param_t *param;
        real8  F2[3];
        real8 dF2dx[3][3];

        param = home->param;
  
        Ecore  = param->Ecore;

        Matrix33_Zero(dverrdx);

/*
 *      Line tension forces and derivatives
 *      Need the force on node 1, CoreForceandDerivative returns the force on node 2
 */
        if (calcForceDeriv)
        {
            CoreForceandDerivative(home, position[0], position[1], position[2],
                                   A[0], A[1], A[2], bA, Ecore, F2, dF2dx);

            for (i = 0; i < 3; i++) {
                verr[i] = -F2[i];
                for (j = 0; j < 3; j++) {
                    dverrdx[i][j] = -dF2dx[i][j]; 
                }
            }
        }
        else
        {
            CoreForce(home, position[0], position[1], position[2],
                      A[0], A[1], A[2], bA, Ecore, F2);
            
            for (i = 0; i < 3; i++) 
                verr[i] = -F2[i];
        }
        
/*
 *      Line tension force on all arms
 */
        for (arm = 0; arm < armCount; arm++) 
        {
            int    i, j;
            real8  nbrPos[3], burg[3];
            Node_t *nbrNode;

            nbrNode = GetNeighborNode(home, node, armList[arm]);
      
            nbrPos[X] = nbrNode->x;
            nbrPos[Y] = nbrNode->y;
            nbrPos[Z] = nbrNode->z;

            FoldBox(param, &nbrPos[X], &nbrPos[Y], &nbrPos[Z]);
      
            burg[X] = -node->burgX[armList[arm]];
            burg[Y] = -node->burgY[armList[arm]];
            burg[Z] = -node->burgZ[armList[arm]];
      
/* 
 *          Add up forces/derivatives of forces coming from line tension
 *          for that arm 
 */

            if (calcForceDeriv)
            {
                CoreForceandDerivative(home, position[0], position[1], position[2],
                                       nbrPos[0], nbrPos[1], nbrPos[2],
                                       burg, Ecore, F2, dF2dx);
          
                for (i = 0; i < 3; i++) {
                    verr[i] -= F2[i];
                    for (j = 0; j < 3; j++) {
                        dverrdx[i][j] -= dF2dx[i][j]; 
                    }
                }
            }
            else
            {
                CoreForce(home, position[0], position[1], position[2],
                          nbrPos[0], nbrPos[1], nbrPos[2],
                          burg, Ecore, F2);
                
                for (i = 0; i < 3; i++)
                    verr[i] -= F2[i];
    
            }

        }  /* end for (arm = 0; ... ) */

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     CalcJunctionLen
 *      Description:  Use Newton-Raphson method to compute the junction
 *                    length where the function vanishes.  The function
 *                    resorts to a bisection method if Newton-Raphson fails.
 *
 *      Arguments:
 *          node        IN: pointer to node being evaluated
 *          armCount    IN: number of node's arms to be included in evaluation
 *          armList     IN: list (indices) of node's arms to be included in
 *                          evaluation
 *          dir{xyz}    IN: components of velocity direction
 *          splitDist   IN: minimum junction length.
 *
 *      Returns:  the calculated optimal junction length or <splitDist> if
 *                unable to calculate the length.
 *
 *------------------------------------------------------------------------*/
real8 CalcJunctionLen(Home_t *home, Node_t *node,
                      int armCount, int *armList,
                      real8 dirx, real8 diry, real8 dirz, 
                      real8 splitDist,int isOpposite)
{
        int   iterNR, maxIterNR, iterBS, maxIterBS, arm;
        real8 maxSegLen;
        real8 lenJunc, maxL;
        real8 error, errorTol;
        real8 position[3], dir[3];
        real8 A[3], bA[3], f[3], f1[3], f2[3];
        real8 verr[3], dverrdx[3][3];

        real8 positionfinal[3];
        real8 FrefInit[3],FrefFinal[3];
        real8 Finitdotdir, Ffinaldotdir;

/*
 *      Newton-Raphson max iter
 */
        iterNR = 1;
        maxIterNR = 10;

/*
 *      Bisection max iter
 */
        iterBS = 1;
        maxIterBS = 300;

        error = 1.0;
        errorTol = 1.0e-06;

        maxSegLen = home->param->maxSeg;
        maxL = 2*maxSegLen;

        VECTOR_ZERO(bA);

        dir[0] = dirx;
        dir[1] = diry;
        dir[2] = dirz;

/*
 *      Initialization : compute maxL and bA
 */
        for (arm = 0; arm < armCount; arm++) {
            real8  normL;
            real8  dx, dy, dz;
            Node_t *nbrNode; 

            nbrNode = GetNeighborNode(home, node, armList[arm]);
      
            dx = nbrNode->x - node->x;
            dy = nbrNode->y - node->y;
            dz = nbrNode->z - node->z;
      
            ZImage(home->param, &dx, &dy, &dz);
            normL = sqrt(dx*dx + dy*dy + dz*dz);

/*
 *          The length of the junction should be limited to
 *          the maximum segment length minus the length of the longest
 *          segment already associated with the node.  This will prevent
 *          any of the segments attached to the node from exceeding
 *          the param->maxSeg length when the node is moved during
 *          the topological operation.  However, never choose a junction
 *          length smaller than the base split distance.
 */
            maxL = MIN(maxL, (maxSegLen - normL));
            maxL = MAX(splitDist, maxL);

            bA[0] += node->burgX[armList[arm]];
            bA[1] += node->burgY[armList[arm]];
            bA[2] += node->burgZ[armList[arm]];
        }

        if (isOpposite) {
            maxL = MIN(maxL, 0.5*maxSegLen);
        } else {
            maxL = MIN(maxL, maxSegLen);
        }

        if (maxL <= splitDist) {
            return splitDist;
        }


/*
 *      A is a point on the segment dir direction where node sits
 */
        A[0] = node->x -  2*home->param->maxSeg * dir[0];
        A[1] = node->y -  2*home->param->maxSeg * dir[1];
        A[2] = node->z -  2*home->param->maxSeg * dir[2];


        lenJunc = splitDist;


/*
 *      Check that there is a zero in the interval 
 */
        position[0] = node->x + splitDist*dir[0];
        position[1] = node->y + splitDist*dir[1];
        position[2] = node->z + splitDist*dir[2];

        CalcNodeLineTension(home, position, node, 
                            armCount, armList, 0,
                            A, bA, FrefInit, dverrdx);
        //SA Do not need dverrdx here.

        positionfinal[0] = node->x + maxL *dir[0];
        positionfinal[1] = node->y + maxL *dir[1];
        positionfinal[2] = node->z + maxL *dir[2];

        CalcNodeLineTension(home, positionfinal, node, 
                            armCount, armList, 0,
                            A, bA, FrefFinal, dverrdx);

        //SA Do not need dverrdx here.

        Finitdotdir=DotProduct(FrefInit,dir);
        Ffinaldotdir=DotProduct(FrefFinal,dir);

/*
 *      Test if there is a zero in [splitDist, maxL].
 *      If no zero, return the lowest energy length
 */
        if (Finitdotdir * Ffinaldotdir > 0.0) {

            if (fabs(Finitdotdir) <= fabs(Ffinaldotdir)) {
                return(splitDist);
            } else {
                return(maxL);
            }
        }

/*
 *      One of splitDist or maxL is zero.
 */
        while ((fabs(error) > errorTol) && (iterNR < maxIterNR) && (lenJunc < maxL) ) 
        {
            int   i, j;
            real8 derrordl;
            real8 temp[3], juncVec[3];

            CalcNodeLineTension(home, position, node, 
                                armCount, armList, 1,
                                A, bA, verr, dverrdx);

/*
 *          Newton-Raphson
 */
            error = DotProduct(verr, dir);

            VECTOR_ZERO(temp);
      
            for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                    temp[i] -= dverrdx[i][j]*dir[j];
                    // SA: sign change here !
                }
            }
      
            derrordl = DotProduct(temp, dir);


            if (fabs(error) > errorTol) 
            {
/*
 *              If derrordl is zero, the Newton-Raphson method isn't going
 *              to work, so break out of the loop and try the bisection
 *              instead.
 */
                if (derrordl == 0.0) { break; }

                position[0] -= (error / derrordl) * dir[0];
                position[1] -= (error / derrordl) * dir[1];
                position[2] -= (error / derrordl) * dir[2];
            }

            juncVec[0] = position[0] - node->x;
            juncVec[1] = position[1] - node->y;
            juncVec[2] = position[2] - node->z;

            lenJunc = sqrt(DotProduct(juncVec, juncVec));

            iterNR++;

        }  /* end while(fabs(error) > errorTol ... ) */


/*
 *      if Newton-Raphson has failed, try bisection
 */
        if ((lenJunc >= maxL) || (lenJunc <= 0.0) || (iterNR >= maxIterNR)) 
        {
            real8 lengthMax;
            real8 fdotdir, f1dotdir, f2dotdir;
            real8 position1[3], position2[3];

            position1[0] = node->x + splitDist*dir[0];
            position1[1] = node->y + splitDist*dir[1];
            position1[2] = node->z + splitDist*dir[2];

            CalcNodeLineTension(home, position1, node, 
                                armCount, armList, 0,
                                A, bA, f1, dverrdx);

            //SA : Do not need dverrdx here

            f1dotdir = DotProduct(f1,dir);

            lengthMax = maxL + (maxL * 0.001);

            position2[0] = node->x + lengthMax * dir[0];
            position2[1] = node->y + lengthMax * dir[1];
            position2[2] = node->z + lengthMax * dir[2];

            CalcNodeLineTension(home, position2, node, 
                                armCount, armList, 0,
                                A, bA, f2, dverrdx);

            //SA : Do not need dverrdx here

            f2dotdir = DotProduct(f2,dir);

            if (f1dotdir * f2dotdir < 0.0) 
            { 
                if (fabs(f1dotdir) < errorTol) 
                {
                    VECTOR_COPY(position, position1);

                    lenJunc = sqrt((position[0]-node->x)*(position[0]-node->x) +
                                   (position[1]-node->y)*(position[1]-node->y) +
                                   (position[2]-node->z)*(position[2]-node->z));

                    return(lenJunc);
                }

                if (fabs(f2dotdir) < errorTol) 
                {
                    position[0] = position2[0];
                    position[1] = position2[1];
                    position[2] = position2[2];

                    lenJunc = sqrt((position[0]-node->x)*(position[0]-node->x) +
                                   (position[1]-node->y)*(position[1]-node->y) +
                                   (position[2]-node->z)*(position[2]-node->z));

                    return(lenJunc);
                }

                iterBS = 1;
                fdotdir = MIN(fabs(f1dotdir), fabs(f2dotdir));  

                while (fabs(fdotdir) >= errorTol && iterBS < maxIterBS) 
                {
                    position[0] = position1[0]+0.5*(position2[0]-position1[0]);
                    position[1] = position1[1]+0.5*(position2[1]-position1[1]);
                    position[2] = position1[2]+0.5*(position2[2]-position1[2]);

                    CalcNodeLineTension(home, position, node, 
                                        armCount, armList, 0,
                                        A, bA, f, dverrdx);

                    //SA : Do not need dverrdx here

                    fdotdir = DotProduct(f, dir);

                    iterBS ++;

                    if (fdotdir * DotProduct(f2, dir) < 0.0) {
                        position1[0] = position[0];
                        position1[1] = position[1];
                        position1[2] = position[2];
                    } else {
                        position2[0] = position[0];
                        position2[1] = position[1];
                        position2[2] = position[2];
                    }

                } /* end while */

/* 
 *              Bisection is a success
 */
                lenJunc = sqrt((position[0]-node->x)*(position[0]-node->x) + 
                               (position[1]-node->y)*(position[1]-node->y) + 
                               (position[2]-node->z)*(position[2]-node->z));

/* 
 *              Check that we have actually reduced forces 
 */
                CalcNodeLineTension(home, position, node, 
                                    armCount, armList, 0,
                                    A, bA, f, dverrdx);

                //SA : Do not need dverrdx here

                fdotdir = DotProduct(f,dir);

                if (fabs(fdotdir) > fabs(Finitdotdir) ) {
                    lenJunc  = splitDist;
                }

                return(lenJunc);
            } 
            else 
            {
                lenJunc = splitDist;

                return(lenJunc);
            }
        }

/* 
 *      Newton Raphson was successful, so check that we have
 *      actually reduced force 
 */
        CalcNodeLineTension(home, position, node, 
                            armCount, armList, 0,
                            A, bA, verr, dverrdx);

        //SA : Do not need dverrdx here

        Ffinaldotdir = DotProduct(verr,dir);

        if (fabs(Ffinaldotdir) > fabs(Finitdotdir)) {
            lenJunc  = splitDist;
        }

        return(lenJunc);
}

/*-------------------------------------------------------------------------
 *
 *      Function:     VelocityWithoutConstraint
 *      Description:  Computes f(v)=v(F)-v and df/dv for Newton-Raphson
 *                    procedure to get best velocity v and direction d=v/|v|
 *                    without accounting for any glide constraint(s)
 *                    associated with the node.
 *
 *      Arguments:
 *          IN:   node       Pointer to node of interest.  NOTE: The force
 *                           and velocity of this node will be updated within
 *                           this function!
 *          OUT:  verror     f(v)
 *          OUT:  deltav     Newton-Raphson update
 *          IN/OUT: mobArgs  Structure containing additional
 *                           parameters for conveying information
 *                           to/from the mobility function.
 *
 *      Returns: 0 on success
 *               1 on error
 *
 *------------------------------------------------------------------------*/
static int VelocityWithoutConstraint(Home_t *home, Node_t *node,
                                     real8 deltav[3], real8 verror[3],
                                     MobArgs_t *mobArgs)
{
        int     i, segID;
        int     mobError;
        real8   Ecore;
        real8   eps;
        real8   dx, dy, dz;
        real8   vold[3];
        real8   f2[3], df2dvold[3][3];
        real8   bA[3];
        real8   Dverror[3][3], DverrorInv[3][3];
        Node_t  *nbrNode;
        Param_t *param;

        param = home->param;

        Ecore = param->Ecore;

        eps = 1e-8;

        vold[0] = node->vX;
        vold[1] = node->vY;
        vold[2] = node->vZ;

        segID = -1;

        for (i = 0; i < node->numNbrs; i++) 
        {
            real8 normD ;

            nbrNode = GetNeighborNode(home, node, i);

            dx = nbrNode->x - node->x;
            dy = nbrNode->y - node->y;
            dz = nbrNode->z - node->z;

            normD = sqrt(dx*dx + dy*dy + dz*dz);

/*
 *          Look for the zero length segment and get its
 *          Burgers vector
 */
            if (normD < eps) {
                segID = i;
                break;
            }
        }

/*
 *      If we did not find a zero length segment, it means that
 *      burgers vector was conserved during the node split and
 *      no new segment was created, so there's nothing to do here.
 */
        if (segID == -1) {
            return(1);
        }

        bA[0] = node->burgX[segID];
        bA[1] = node->burgY[segID];
        bA[2] = node->burgZ[segID];


        CoreForceandDerivative(home, 0.0,0.0,0.0, vold[0], vold[1], vold[2], 
                               bA, Ecore, f2, df2dvold);

        /* The function CoreForceandDerivative returns f2 and df2dvold */

        node->armfx[segID] = -f2[0];
        node->armfy[segID] = -f2[1];
        node->armfz[segID] = -f2[2];

        Matrix33_Mul(Dverror, mobArgs->invDragMatrix, df2dvold);  // Dverror = invDragMatrix*df2dvold

        Dverror[0][0] -= 1.0;
        Dverror[1][1] -= 1.0;
        Dverror[2][2] -= 1.0;

        if ( Matrix33_Inverse(DverrorInv, Dverror) < 0 ) { return(1); }

/*
 *      F = sum (F) over its arms + F_{line tension} (bA, dir)
 */
        node->fX = 0.0;
        node->fY = 0.0;
        node->fZ = 0.0;

        for (i = 0; i < node->numNbrs; i++) {
            node->fX += node->armfx[i];
            node->fY += node->armfy[i];
            node->fZ += node->armfz[i];
        }

        mobError  = EvaluateMobility(home, node, mobArgs);

        if (mobError) {
            return(1);
        }

        verror[0] = node->vX - vold[0];
        verror[1] = node->vY - vold[1];
        verror[2] = node->vZ - vold[2];

        node->vX = vold[X];
        node->vY = vold[Y];
        node->vZ = vold[Z];

        Matrix33Vector3Multiply(DverrorInv, verror, deltav);    

        return(0);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     AdjustJuncNodeForceAndVel
 *      Description:  This function attempts to calculate the optimal
 *                    direction for the node's junction segment using
 *                    Newton-Raphson and adjusts the node's force
 *                    and velocity appropriately.
 *
 *      Arguments:
 *          IN:     node       Pointer to node of interest.
 *                             NOTE: The force and velocity of this
 *                             node will be updated within this
 *                             function!
 *          IN/OUT: mobArgs    Structure containing additional
 *                             parameters for conveying information
 *                             to/from the mobility function.
 *
 *------------------------------------------------------------------------*/
void AdjustJuncNodeForceAndVel(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int   iter, maxIter,ierr;
        real8 DotV;
        real8 errorTol;
        real8 fSave[3], vSave[3];
        real8 verror[3], verrorSave[3];
        real8 deltav[3];


        if (mobArgs->numGlideConstraints > 2) {
/*
 *          If there are more than two glide constraints, the velocity is
 *          explicitly zeroed.
 */
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;

            return;

        } else if (mobArgs->numGlideConstraints == 2) {
/*
 *          With exactly two glide constraints, the velocity direction
 *          is fixed.
 */
            DotV = node->vX * mobArgs->glideConstraints[2][X] + 
                   node->vY * mobArgs->glideConstraints[2][Y] + 
                   node->vZ * mobArgs->glideConstraints[2][Z];

            node->vX = DotV * mobArgs->glideConstraints[2][X];
            node->vY = DotV * mobArgs->glideConstraints[2][Y];
            node->vZ = DotV * mobArgs->glideConstraints[2][Z];

            return;
        }

/*
 *      There is either 0 or 1 glide constraint, so try Newton Raphson
 *
 *      Keep a copy of the node's intial force and velocity in case
 *      Newton-Raphson fails
 */
        fSave[0] = node->fX;
        fSave[1] = node->fY;
        fSave[2] = node->fZ;

        vSave[0] = node->vX;
        vSave[1] = node->vY;
        vSave[2] = node->vZ;

/*
 *      Initialize parameters
 */
        iter = 1;
        maxIter = 100;
        errorTol = 1.0e-8 * Normal(vSave);
        ierr = 0;

        ierr = VelocityWithoutConstraint(home, node, deltav, verror, mobArgs);

        if (mobArgs->numGlideConstraints == 1)
        {
           real8 S1;
           S1 = DotProduct(verror, mobArgs->glideConstraints[0]);
           verror[0] -= S1 * mobArgs->glideConstraints[0][0];
           verror[1] -= S1 * mobArgs->glideConstraints[0][1];
           verror[2] -= S1 * mobArgs->glideConstraints[0][2];
                      
           S1 = DotProduct(deltav, mobArgs->glideConstraints[0]);
           deltav[0] -= S1 * mobArgs->glideConstraints[0][0];
           deltav[1] -= S1 * mobArgs->glideConstraints[0][1];
           deltav[2] -= S1 * mobArgs->glideConstraints[0][2];
        }


        verrorSave[0] = verror[0];
        verrorSave[1] = verror[1];
        verrorSave[2] = verror[2];
   
/*
 *      If we did not (or could not) get velocity/force updates from
 *      the velocity functions above, just restore the original force
 *      and velocity and return to the caller.
 */
        if (ierr) {

            node->fX = fSave[0];
            node->fY = fSave[1];
            node->fZ = fSave[2];

            node->vX = vSave[0];
            node->vY = vSave[1];
            node->vZ = vSave[2];

            return;
        }

/*
 *      Start Newton-Raphson iterations
 */
        while (Normal(verror) > errorTol && iter < maxIter) {

            node->vX -= deltav[X];
            node->vY -= deltav[Y];
            node->vZ -= deltav[Z];

            ierr = VelocityWithoutConstraint(home, node, deltav, verror,
                                             mobArgs);

        if (mobArgs->numGlideConstraints == 1)
        {
           real8 S1;
           S1 = DotProduct(verror, mobArgs->glideConstraints[0]);
           verror[0] -= S1 * mobArgs->glideConstraints[0][0];
           verror[1] -= S1 * mobArgs->glideConstraints[0][1];
           verror[2] -= S1 * mobArgs->glideConstraints[0][2];
                      
           S1 = DotProduct(deltav, mobArgs->glideConstraints[0]);
           deltav[0] -= S1 * mobArgs->glideConstraints[0][0];
           deltav[1] -= S1 * mobArgs->glideConstraints[0][1];
           deltav[2] -= S1 * mobArgs->glideConstraints[0][2];
        }

/*
 *          If we did not (or could not) get velocity/force updates from
 *          the velocity functions above, just restore the original force
 *          and velocity and return to the caller.
 */
            if (ierr) {

                node->fX = fSave[0];
                node->fY = fSave[1];
                node->fZ = fSave[2];

                node->vX = vSave[0];
                node->vY = vSave[1];
                node->vZ = vSave[2];

                return;
            }

            iter++;

        }  /* end Newton Raphson loop */


        if (ierr == 1 || Normal(verror) > errorTol || iter > maxIter) {
/*
 *          Newton-Raphson has failed
 */
            if (Normal(verror) > Normal(verrorSave)) {

                node->fX = fSave[0];
                node->fY = fSave[1];
                node->fZ = fSave[2];

                node->vX = vSave[0];
                node->vY = vSave[1];
                node->vZ = vSave[2];
            }
        }

/*
 *      In all other cases, just return the current velocity and force
 */
        return;
}
