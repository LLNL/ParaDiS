/**************************************************************************
 *
 *      Module:  NodeForce -- This module contains various functions for
 *               calculating nodal forces
 *
 *      Includes functions:
 *               AddtoNodeForce()
 *               AddtoArmForce()
 *               ComputeSegSigbRem()
 *               ExtPKForce()
 *               NodeForce()
 *               PKForce()
 *               SelfForce()
 *
 *
 *      NOTE:  This module contains several blocks of code which are
 *             conditionally compiled based on the _ARLFEM definition.
 *             This pre-processor definition is only set when ParaDiS
 *             is being coupled with finite element modules which are
 *             not included as part of the ParaDiS release.
 *             Therefore, these blocks of code can be ignored by
 *             non-LLNL developers.
 *
 *************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "Mobility.h"
#include "V3.h"

#ifdef ANISOTROPIC
#include "Anisotropic.h"
#endif

#ifdef _ARLFEM
#include "FEM.h"
#include "DSMCommunicator.h"
#endif

/*
 *      Define a local structure in which to store some basic info we
 *      need when setting up a list of segment pairs for which forces
 *      need to be calculated
 */
typedef struct {
        Node_t *node1;
        Node_t *node2;
        Node_t *node3;
        Node_t *node4;
        real8  cellCenter[3];
} SegPair_t;

static real8 *cellCenterX = 0;
static real8 *cellCenterY = 0;
static real8 *cellCenterZ = 0;

/*---------------------------------------------------------------------------
 *    Function:        InitCellCenters
 *    Description:     Allocate and initialize the cell centers
 *-------------------------------------------------------------------------*/

void InitCellCenters(Home_t *home)
{
    const Param_t *param = home->param;

    const int   nx = param->nXcells;
    const int   ny = param->nYcells;
    const int   nz = param->nZcells;

    const real8 lx = param->Lx;
    const real8 ly = param->Ly;
    const real8 lz = param->Lz;

    const real8 dx = lx/nx;
    const real8 dy = ly/ny;
    const real8 dz = lz/nz;

    const real8 x0 = param->minSideX + (dx/2.0);
    const real8 y0 = param->minSideY + (dy/2.0);
    const real8 z0 = param->minSideZ + (dz/2.0);

    if (cellCenterX) { free( cellCenterX ); cellCenterX=0; }
    if (cellCenterY) { free( cellCenterY ); cellCenterY=0; }
    if (cellCenterZ) { free( cellCenterZ ); cellCenterZ=0; }

    cellCenterX = (real8 *) malloc(nx*sizeof(real8));
    cellCenterY = (real8 *) malloc(ny*sizeof(real8));
    cellCenterZ = (real8 *) malloc(nz*sizeof(real8));

    if (cellCenterX) { for (int i=0; (i<nx); i++) { cellCenterX[i] = x0 + i*dx; } }
    if (cellCenterY) { for (int i=0; (i<ny); i++) { cellCenterY[i] = y0 + i*dy; } }
    if (cellCenterZ) { for (int i=0; (i<nz); i++) { cellCenterZ[i] = z0 + i*dz; } }
}

/*---------------------------------------------------------------------------
 *
 *    Function:        FreeCellCenters
 *    Description:     Release memory allocated to arrays for storing
 *                     the coordinates of the cell centers.  These arrays
 *                     are used in the GetFieldPointStressRem() function.
 *
 *-------------------------------------------------------------------------*/
void FreeCellCenters(void)
{
    if (cellCenterX) { free(cellCenterX); cellCenterX=0; }
    if (cellCenterY) { free(cellCenterY); cellCenterY=0; }
    if (cellCenterZ) { free(cellCenterZ); cellCenterZ=0; }
}

/*---------------------------------------------------------------------------
 *
 *    Function:        AddtoNodeForce
 *    Description:     Increment the total nodal force for the indicated node
 *
 *-------------------------------------------------------------------------*/
void AddtoNodeForce(Node_t *node, real8 f[3])
{
    LOCK(&node->nodeLock);

    node->fX+=f[0];
    node->fY+=f[1];
    node->fZ+=f[2];

    UNLOCK(&node->nodeLock);
}

/*---------------------------------------------------------------------------
 *
 *    Function:        AddtoArmForce
 *    Description:     Increment the nodal force attributable to
 *                     a specific arm of a node by the specified
 *                     amount.
 *
 *-------------------------------------------------------------------------*/
void AddtoArmForce(Node_t *node, int arm, real8 f[3])
{

    LOCK(&node->nodeLock);

    node->armfx[arm] += f[0];
    node->armfy[arm] += f[1];
    node->armfz[arm] += f[2];

    UNLOCK(&node->nodeLock);

    return;
}

/*---------------------------------------------------------------------------
 *
 *    Function:        ZeroNodeForces
 *    Description:     Initialize all forces
 *
 *-------------------------------------------------------------------------*/

void ZeroNodeForces(Home_t *home, int reqType)
{
    ZeroNodeForces(home, home->nodeKeys     , home->newNodeKeyPtr,  reqType );
    ZeroNodeForces(home, home->ghostNodeList, home->ghostNodeCount, reqType );
}

void ZeroNodeForces
(
    Home_t    *home   ,   ///< points to home data structure
    Node_t   **nodes  ,   ///< array of node pointers
    int        ncnt   ,   ///< number of nodes
    int        reqType    ///< request type
)
{
    if (nodes && (ncnt>0))
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i=0; i<ncnt; i++)
        {
            Node_t *node = nodes[i];

            if (node)
            {
                // If we're doing a FULL zeroing of forces, or the node is flagged
                // to have forces reset, zero it's total forces.

                if ((reqType == FULL) || (node->flags & NODE_RESET_FORCES) != 0)
                {
                    node->fX = 0.0;
                    node->fY = 0.0;
                    node->fZ = 0.0;
                }

                // Loop through all the arms;  If we're doing a FULL zeroing of
                // forces OR either segment endpoint is flagged to have forces
                // reset, zero the arm's forces.

                for (int j=0; j < node->numNbrs; j++)
                {
                    Node_t *nbr = GetNeighborNode(home,node,j);

                    if (nbr)
                    {
                        if (     (reqType == FULL)
                             || ((node->flags & NODE_RESET_FORCES) != 0)
                             || ((nbr ->flags & NODE_RESET_FORCES) != 0) )
                        {
                            node->armfx[j] = 0.0;
                            node->armfy[j] = 0.0;
                            node->armfz[j] = 0.0;
                        }
                    }
                }

                // If we're doing a full reset of forces, we can zero the
                // 'reset forces' flag here; if we're only doing partial
                // forces, we still need the flag and will reset it elsewhere
                // when we no longer need it.

                if (reqType == FULL)
                {
                    node->flags &= (~NODE_RESET_FORCES);
                }
            }
        }
    }
}


void GetFieldPointStressRem(Home_t *home, real8 x, real8 y, real8 z,
                            int cellX, int cellY, int cellZ,
                            real8 totStress[3][3])
{
        int     xSkip1, xSkip2, xSkip3;
        int     ySkip1, ySkip2, ySkip3;
        int     zSkip1, zSkip2, zSkip3;
        real8   delSig[3][3];

        Param_t *param = home->param;

/*
 *      First time into this function, allocate and initialize some
 *      static arrays
 */
        if ( !cellCenterX || !cellCenterY || !cellCenterZ )
        {
            InitCellCenters(home);
        }

/*
 *      Initialize the stress at the field point
 */
        for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
           totStress[i][j] = 0.0;

/*
 *      Since we have to treat neighboring and remote cells
 *      differently, we need to identify which cells are which...
 *      If the neighbor cells are outside the primary image, wrap
 *      the cell indices around to the other side of the box ONLY
 *      if periodic boundaries are enabld!
 */
         xSkip1 = cellX - 1 ;
         if (xSkip1 < 0) {
             if (param->xBoundType == Periodic)
                 xSkip1 = param->nXcells - 1 ;
             else
                 xSkip1 = 0;
         }
         xSkip2 = cellX ;
         xSkip3 = cellX + 1 ;
         if (xSkip3 >= param->nXcells) {
             if (param->xBoundType == Periodic)
                 xSkip3 = 0 ;
             else
                 xSkip3 = param->nXcells - 1 ;
         }

         ySkip1 = cellY - 1 ;
         if (ySkip1 < 0) {
             if (param->yBoundType == Periodic)
                 ySkip1 = param->nYcells - 1 ;
             else
                 ySkip1 = 0;
         }
         ySkip2 = cellY ;
         ySkip3 = cellY + 1 ;
         if (ySkip3 >= param->nYcells) {
             if (param->yBoundType == Periodic)
                 ySkip3 = 0 ;
             else
                 ySkip3 = param->nYcells - 1;
         }

         zSkip1 = cellZ - 1 ;
         if (zSkip1 < 0) {
             if (param->zBoundType == Periodic)
                 zSkip1 = param->nZcells - 1 ;
             else
                 zSkip1 = 0;
         }
         zSkip2 = cellZ ;
         zSkip3 = cellZ + 1 ;
         if (zSkip3 >= param->nZcells) {
             if (param->zBoundType == Periodic)
                 zSkip3 = 0 ;
             else
                 zSkip3 = param->nZcells - 1 ;
         }


/*
 *      Loop through all cells, and add the stress contribution from
 *      the cell to stress at the segment midpoint.
 */
        for (int cx=0; cx < param->nXcells; cx++)
        {
          for (int cy=0; cy < param->nYcells; cy++)
          {
            for (int cz=0; cz < param->nZcells; cz++)
            {
                int includePrimary = !(    (cx==xSkip1 || cx==xSkip2 || cx==xSkip3)
                                        && (cy==ySkip1 || cy==ySkip2 || cy==ySkip3)
                                        && (cz==zSkip1 || cz==zSkip2 || cz==zSkip3) );

                // Get the center point of cell [cx, cy, cz]

                real8 xc = cellCenterX[cx];
                real8 yc = cellCenterY[cy];
                real8 zc = cellCenterZ[cz];

                // Get the stress at the specified point caused
                // by the net charge tensor of the current cell.

                real8 dx = xc - x;
                real8 dy = yc - y;
                real8 dz = zc - z;

                ZImage(param, &dx, &dy, &dz);

                xc = x + dx;
                yc = y + dy;
                zc = z + dz;

                int cellIndex = cz + param->nZcells*cy +
                                param->nZcells*param->nYcells*cx;

                // Stress (charge[.,1], [1,0,0])

                real8 burgX = home->cellCharge[9*cellIndex  ];
                real8 burgY = home->cellCharge[9*cellIndex+3];
                real8 burgZ = home->cellCharge[9*cellIndex+6];

                dx = 1;
                dy = 0;
                dz = 0;

                dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz,
                              burgX, burgY, burgZ, x, y, z, includePrimary);

                for (int i=0; i<3; i++)
                for (int j=0; j<3; j++)
                    totStress[i][j] += delSig[i][j];

                // Stress (charge[.,2], [0,1,0])

                burgX = home->cellCharge[9*cellIndex+1];
                burgY = home->cellCharge[9*cellIndex+4];
                burgZ = home->cellCharge[9*cellIndex+7];

                dx = 0;
                dy = 1;
                dz = 0;

                dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz,
                              burgX, burgY, burgZ, x, y, z, includePrimary);

                for (int i=0; i<3; i++)
                for (int j=0; j<3; j++)
                    totStress[i][j] += delSig[i][j];

                // Stress (charge[.,3], [0,0,1])

                burgX = home->cellCharge[9*cellIndex+2];
                burgY = home->cellCharge[9*cellIndex+5];
                burgZ = home->cellCharge[9*cellIndex+8];

                dx = 0;
                dy = 0;
                dz = 1;

                dSegImgStress(home, delSig, xc, yc, zc, dx, dy, dz,
                              burgX, burgY, burgZ, x, y, z, includePrimary);

                for (int i=0; i<3; i++)
                for (int j=0; j<3; j++)
                    totStress[i][j] += delSig[i][j];

            } /* end for(cz = 0; ...) */
          } /* end for(cy = 0; ...) */
        } /* end for(cx = 0; ...) */

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     ComputeSigSegbRem
 *      Description:  For each native segment, calculate the resultant
 *                    stress at the segment midpoint due to the cellCharges
 *                    (accumulated bXdl in a cell) from all cells. Store
 *                    the dot product of this tensor with the segment's
 *                    burger's vector (i.e sigma dot b) in the segment's
 *                    remote sigb.
 *
 *                    NOTE: For remote cells, the resultant stress on each
 *                          segment includes contributions from both the
 *                          primary and periodic images.  For the segment's
 *                          local and neighboring cells, this stress includes
 *                          contributions only from the periodic images.
 *                          Force from interactions with the primary image
 *                          have been computed elsewhere.
 *      Arguments:
 *          reqType   Indicates the type of calculation being requested.
 *                    Permitted values are:
 *
 *                        PARTIAL  Remote sigb is calculated only for nodes
 *                                 flagged to have forces recalculated.
 *                        FULL     Remote sigb is calculated for all nodes.
 *
 *-----------------------------------------------------------------------*/
void ComputeSegSigbRem (Home_t *home, int reqType)
{
        int     inode, ti, myXcell, myYcell, myZcell;
        int     armID1, armID2;
        real8   x1, y1, z1, dx, dy, dz, xm, ym, zm;
        real8   bx1, by1, bz1, sigb1, sigb2, sigb3;
        real8   totRemSig[3][3];
        Param_t *param;
        Node_t  *node, *nbr;
        Cell_t  *cell;

        param = home->param;

/*
 *      loop thru the native nodes
 */
        for (inode = 0; inode < home->newNodeKeyPtr; inode++)
        {
            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) { continue; }

/*
 *          loop thru the native node's arms.  If the segment is owned
 *          by the neighboring node, don't do the calculation in this
 *          iteration of the loop.
 */
            for (ti = 0; ti < node->numNbrs; ti++)
            {
                nbr = GetNeighborNode(home, node, ti);

                if (nbr == (Node_t *)NULL)
                {
                    printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                    continue;
                }

                if (NodeOwnsSeg(home, node, nbr) == 0) { continue; }

/*
 *              If we're only calculating the sigb for a portion of the
 *              nodes, skip this segment if neither node is flagged for
 *              a force update.
 */
                if ((reqType == PARTIAL) &&
                    (((node->flags & NODE_RESET_FORCES) == 0) &&
                     ((nbr->flags  & NODE_RESET_FORCES) == 0)))
                { continue; }

                armID1 = GetArmID(node, nbr);
                armID2 = GetArmID(nbr, node);

                x1 = node->x;
                y1 = node->y;
                z1 = node->z;

/*
 *              Initialize the sigb on the arm
 */
                node->sigbRem[3*armID1  ] = 0.0;
                node->sigbRem[3*armID1+1] = 0.0;
                node->sigbRem[3*armID1+2] = 0.0;

/*
 *              Get the midpoint of the segment
 */
                dx = nbr->x - x1;
                dy = nbr->y - y1;
                dz = nbr->z - z1;

                ZImage(param, &dx, &dy, &dz);

                xm = x1 + dx*0.5;
                ym = y1 + dy*0.5;
                zm = z1 + dz*0.5;

/*
 *              Get the cell indices for the cell containing the node
 *              and remove adjustments for periodic cells, then find
 *              the stress at the specied point from all segments in
 *              remote cells.
 */
                cell = LookupCell(home, node->cellIdx);

                myXcell = cell->xIndex;
                myYcell = cell->yIndex;
                myZcell = cell->zIndex;

                myXcell--;
                myYcell--;
                myZcell--;

                GetFieldPointStressRem(home, xm, ym, zm,
                                       myXcell, myYcell, myZcell, totRemSig);

/*
 *              Calculate Sigma dot b from remote stess on this segment, and
 *              store in the node arm
 */
                bx1= node->burgX[armID1];
                by1= node->burgY[armID1];
                bz1= node->burgZ[armID1];

                sigb1 = totRemSig[0][0]*bx1 +
                        totRemSig[0][1]*by1 +
                        totRemSig[0][2]*bz1;
                sigb2 = totRemSig[1][0]*bx1 +
                        totRemSig[1][1]*by1 +
                        totRemSig[1][2]*bz1;
                sigb3 = totRemSig[2][0]*bx1 +
                        totRemSig[2][1]*by1 +
                        totRemSig[2][2]*bz1;

                node->sigbRem[3*armID1  ]=sigb1;
                node->sigbRem[3*armID1+1]=sigb2;
                node->sigbRem[3*armID1+2]=sigb3;

                nbr->sigbRem[3*armID2  ]=sigb1;
                nbr->sigbRem[3*armID2+1]=sigb2;
                nbr->sigbRem[3*armID2+2]=sigb3;

            } /* for(ti=0;ti<nc;ti++) */

        } /* for(inode=0;...) */

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     NodeForce
 *      Description:  This function does no force calculations directly
 *                    but drives the calculations via lower level functions.
 *                    It generates the sig.b from remote segments, the
 *                    forces from direct segment-to-segment interactions,
 *                    then for each local segment, calculates the self-force,
 *                    the force from extrenal stress, and the force from
 *                    the remote sig.b.
 *
 *      Arguments:
 *          reqType   Specifies whether the function is to do a full set
 *                    or force calcs, or just for a subset of the nodes.
 *                    Valid values are: 1 (PARTIAL) or 2 (FULL)
 *
 *-----------------------------------------------------------------------*/
void NodeForce(Home_t *home, int reqType)
{
        int     i, nc, ti, nbrArm;
        real8   f1[3], f2[3];
        Node_t  *node, *nbr;

        Param_t *param      = home->param;
        int      elasticity = param->elasticinteraction;

        TimerStart(home, CALC_FORCE);

#ifdef ESHELBY
/*
 *      If we're simulating eshelby inclusions and using the
 *      segment/particle intersection list, we need to clear the
 *      list now.
 */
        SegPartListClear(home);
#endif

/*
 *      If this is a load-balance-only step, all we need to do
 *      is count the number of force calculations we *would*
 *      do if this were a real cycle.  If elastic interaction
 *      is enabled, a call to LocalSegForces() will set the
 *      count properly.  Otherwise, do a quick loop to figure
 *      out how many line tension calculations we'd do.
 *      Then, return without actually calculating new forces
 *      and updating nodal data.
 */
        if ((param->numDLBCycles > 0) && (reqType == FULL))
        {
            TimerStart(home, LOCAL_FORCE) ;

            if (elasticity)
            {
               LocalSegForces(home, FULL);
            }
            else
            {
                for (i = 0; i < home->newNodeKeyPtr; i++)
                {
                    node = home->nodeKeys[i];
                    if (!node) continue;
                    nc = node->numNbrs;
                    for (ti = 0; ti < nc; ti++)
                    {
                        nbr = GetNeighborNode(home, node, ti);
                        if (nbr == (Node_t *)NULL)
                        {
                            printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                            continue;
                        }
                        if ((nbr->myTag.domainID == home->myDomain) &&
                            (OrderNodes(node, nbr) >= 0)) continue;
                        home->cycleForceCalcCount++;
                    }
                }
            }

            TimerStop(home, LOCAL_FORCE) ;
            TimerStop(home, CALC_FORCE);

            return;
        }


/*
 *      Reset all node forces to zero (local and ghost nodes)
 */
        ZeroNodeForces(home, home->nodeKeys     , home->newNodeKeyPtr,  reqType );   // (zero local nodes)
        ZeroNodeForces(home, home->ghostNodeList, home->ghostNodeCount, reqType );   // (zero ghost nodes)

#ifdef CALCENERGY
        home->param->TotalEnergy = 0.0;   // set total energy to zero
#endif

/*
 *      If elastic interaction is not enabled, use a simple line
 *      tension model for calculating forces.  (Useful for doing
 *      initial testing of code with a quick run.)
 */
        if (!elasticity)
        {
            TimerStart(home, LOCAL_FORCE) ;

            for (i = 0; i < home->newNodeKeyPtr; i++)
            {
                if ((node = home->nodeKeys[i]) == (Node_t *)NULL) { continue; }

                for (ti = 0; ti < node->numNbrs; ti++)
                {
                    nbr = GetNeighborNode(home, node, ti);

                    if (nbr == (Node_t *)NULL)
                    {
                        printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                        continue;
                    }

/*
 *                  If we're only doing partial forces, skip any segments
 *                  for which neither nodal endpoint has been flagged for
 *                  force updates.
 */
                    if (reqType == PARTIAL)
                    {
                        if (((node->flags & NODE_RESET_FORCES) == 0) &&
                            ((nbr->flags & NODE_RESET_FORCES) == 0))
                        { continue; }
                    }

                    if ((nbr->myTag.domainID == home->myDomain) &&
                        (OrderNodes(node, nbr) >= 0)) continue;

                    LineTensionForce(home, node->x, node->y, node->z,
                                     nbr->x, nbr->y, nbr->z, node->burgX[ti],
                                     node->burgY[ti], node->burgZ[ti],
                                     f1, f2);

                    AddtoNodeForce(node, f1);
                    AddtoArmForce(node, ti, f1);


                    if (nbr->myTag.domainID == home->myDomain)
                    {
                        AddtoNodeForce(nbr, f2);
                        nbrArm = GetArmID(nbr, node);
                        AddtoArmForce(nbr, nbrArm, f2);
                    }

#ifdef CALCENERGY
                    real8 W = LineTensionEnergy(home,
                                                node->x, node->y, node->z,
                                                nbr ->x, nbr ->y, nbr ->z,
                                                node->burgX[ti], node->burgY[ti], node->burgZ[ti]);

                    home->param->TotalEnergy += W;
#endif
                }

            }

/*
 *          BUG FIXED: Recompute the nodes total force now.
 */
            for (i = 0; i < home->newNodeKeyPtr; i++)
            {
                if ((node = home->nodeKeys[i]) == (Node_t *)NULL)
                { continue; }

                node->fX = 0.0;
                node->fY = 0.0;
                node->fZ = 0.0;

                for (ti = 0; ti < node->numNbrs; ti++)
                {
                    node->fX += node->armfx[ti];
                    node->fY += node->armfy[ti];
                    node->fZ += node->armfz[ti];
                }
            }

            TimerStop(home, LOCAL_FORCE);
            TimerStop(home, CALC_FORCE);

            return;
        }

        TimerStart(home, REMOTE_FORCE);

/*
 *      If _ARLFEM is defined, we need to add the FEM image stress
 *      to each segment's sigbRem.  If FMM code is also enabled, we'll
 *      need to explicitly initialize the sigbRem values before
 *      computing the FEM contribution (if the FMM code is not enabled,
 *      ComputeSegSigbRem() will explicitly set the values, so no
 *      initialization is necessary).
 *
 *      NOTE: If full n^2 seg/seg forces are being calculated, we don't
 *      do any remote force calcs... unless the FEM code is hooked in,
 *      in which case we still need to factor in some Yoffe stress?
 */
#ifndef FULL_N2_FORCES
#ifdef SPECTRAL
        if ( (param->fmEnabled==0) && (param->FFTenabled==0) )
#else
        if   (param->fmEnabled==0)
#endif
        {
            ComputeSegSigbRem(home, reqType);
        }
#endif

        TimerStop(home, REMOTE_FORCE);

/*
 *      Now handle all the force calculations that must be done by the
 *      local domain.  This includes self-force, PK force, and far-field
 *      forces for all native segments plus any segment-to-segment
 *      interactions for segment pairs 'owned' by this domain.
 *
 *      All force calulations for a given segment will be calculated
 *      only once, and the calculating domain will distribute calculated
 *      forces to remote domains as necessary.
 */
        TimerStart    (home, LOCAL_FORCE);
        LocalSegForces(home, reqType);
        TimerStop     (home, LOCAL_FORCE);

#ifdef ESHELBY
/*
 *      If we're simulating eshelby inclusions and using the segment/particle
 *      intersection list, it should be fully built, but we need to sort
 *      it for later use.
 */
        SegPartListSort(home);
#endif

        TimerStop(home, CALC_FORCE);

#ifdef PARALLEL
#ifdef SYNC_TIMERS
        TimerStart(home, CALC_FORCE_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop(home, CALC_FORCE_BARRIER);
#endif
#endif

}


/*-------------------------------------------------------------------------
 *
 *      Function:       SetOneNodeForce
 *      Description:    Recalculate the force for all segments attached
 *                      to the specified node.  After all segment forces
 *                      are computed, total forces will be resummed for
 *                      for the nodes attached to the segments.
 *                      There are two version of this function.
 *
 *                      The first version is only used when calculating forces
 *                      for full N^2 segment-to-segment interactions.  It
 *                      assumes knowledge of all dislocation segments and no
 *                      remote force calculations are done.  This method is
 *                      only valid with serial compilation and is used
 *                      primarily for debugging or verification purposes.
 *
 *                      The second method is valid for serial or parallel
 *                      execution and computes all necessary local and
 *                      remote forces.
 *
 *      Arguments:
 *              node1           pointer to node for which force values
 *                              are to be updated
 *
 *------------------------------------------------------------------------*/
#ifdef FULL_N2_FORCES
/*
 *      WARNING! When FULL_N2_FORCES is defined, PARALLEL must NOT be
 *               defined.  The FEM code uses DSM which requires a parallel
 *               build, hence FEM cannot currently be used with
 *               FULL_N2_FORCES.
 */
void SetOneNodeForce(Home_t *home, Node_t *nodeA)
{
        int     i, j, k, l;
        int     armID12, armID21, armID34;
        int     isIntersecting;
        real8   dx, dy, dz;
        real8   pos1[3], pos2[3];
        real8   sigb[3], burg[3];
        real8   f1[3], f2[3], f3[3], f4[3];
        real8   extstress[3][3];
        Node_t  *node1, *node2, *node3, *node4;

        Param_t *param = home->param;

        real8    a     = param->rc;
        real8    E     = param->YoungsModulus;
        real8    MU    = param->shearModulus;
        real8    NU    = param->pois;
        real8    Ecore = param->Ecore;
        real8    eps   = 1.0e-6;

        extstress[0][0] = param->appliedStress[0];
        extstress[1][1] = param->appliedStress[1];
        extstress[2][2] = param->appliedStress[2];
        extstress[1][2] = param->appliedStress[3];
        extstress[2][0] = param->appliedStress[4];
        extstress[0][1] = param->appliedStress[5];
        extstress[2][1] = extstress[1][2];
        extstress[0][2] = extstress[2][0];
        extstress[1][0] = extstress[0][1];

/*
 *      Loop over all the segments attached to the node.  Recalculate
 *      the force on each segment (if possible) and reset the
 *      node/segment forces.
 */
        for (int armIndex = 0; armIndex < nodeA->numNbrs; armIndex++)
        {
            Node_t *nodeB = GetNodeFromTag(home, nodeA->nbrTag[armIndex]);

            if (nodeB == (Node_t *)NULL) { continue; }

/*
 *          Make sure node1 is always the node with the lower tag.
 */
            if (OrderNodes(nodeA, nodeB) < 0)
            {
                node1   = nodeA;
                node2   = nodeB;
                armID12 = armIndex;
                armID21 = GetArmID(node2, node1);
            }
            else
            {
                node1   = nodeB;
                node2   = nodeA;
                armID12 = GetArmID(node1, node2);
                armID21 = armIndex;
            }

/*
 *          Zero out the old segment forces.
 */
            node1->armfx[armID12] = 0.0;
            node1->armfy[armID12] = 0.0;
            node1->armfz[armID12] = 0.0;

            node2->armfx[armID21] = 0.0;
            node2->armfy[armID21] = 0.0;
            node2->armfz[armID21] = 0.0;

            pos1[0] = node1->x;
            pos1[1] = node1->y;
            pos1[2] = node1->z;

            dx = node2->x - pos1[0];
            dy = node2->y - pos1[1];
            dz = node2->z - pos1[2];

            ZImage(param, &dx, &dy, &dz);

            pos2[0] = pos1[0] + dx;
            pos2[1] = pos1[1] + dy;
            pos2[2] = pos1[2] + dz;

            burg[0] = node1->burgX[armID12];
            burg[1] = node1->burgY[armID12];
            burg[2] = node1->burgZ[armID12];

/*
 *          If the segment is zero-length (which can occur when we're
 *          testing whether multinodes should be split or not), ignore
 *          the segment.
 */
            if (((pos1[0]-pos2[0])*(pos1[0]-pos2[0])) +
                ((pos1[1]-pos2[1])*(pos1[1]-pos2[1])) +
                ((pos1[2]-pos2[2])*(pos1[2]-pos2[2])) < eps) { continue; }

/*
 *          If elastic interaction is not enabled, use a simple line
 *          tension model for calculating forces on the segment. (i.e.
 *          no segment/segment forces, no remote forces, etc)
 *
 *          Segment forcess are updated here; the nodal forces are
 *          resummed at the end of this routine.
 */
            if (!param->elasticinteraction)
            {
                LineTensionForce(home,
                                 node1->x, node1->y, node1->z,
                                 node2->x, node2->y, node2->z,
                                 burg[0], burg[1], burg[2], f1, f2);

                AddtoArmForce(node1, armID12, f1);
                AddtoArmForce(node2, armID21, f2);

#ifdef CALCENERGY
                real8 W = LineTensionEnergy(home,
                                            node1->x, node1->y, node1->z,
                                            node2->x, node2->y, node2->z,
                                            burg[0], burg[1], burg[2]);

                home->param->TotalEnergy += W;
#endif

/*
 *              We're done messing with the neighbor node, so recompute that
 *              node's total force now.
 */
                nodeB->fX = 0.0;
                nodeB->fY = 0.0;
                nodeB->fZ = 0.0;

                for (int nbrArmID = 0; nbrArmID < nodeB->numNbrs; nbrArmID++)
                {
                    nodeB->fX += nodeB->armfx[nbrArmID];
                    nodeB->fY += nodeB->armfy[nbrArmID];
                    nodeB->fZ += nodeB->armfz[nbrArmID];
                }
                continue;
            }

/*
 *          Calculate the segment self-force, force from external stress, etc.
 */
            SelfForce(home, 0, MU, NU,
                      burg[0], burg[1], burg[2],
                      pos1[0], pos1[1], pos1[2],
                      pos2[0], pos2[1], pos2[2], a, Ecore, f1, f2);

            AddtoArmForce(node1, armID12, f1);
            AddtoArmForce(node2, armID21, f2);


#ifdef CALCENERGY
            real8 W = SelfEnergy(0, MU, NU,
                                 burg[0], burg[1], burg[2],
                                 pos1[0], pos1[1], pos1[2],
                                 pos2[0], pos2[1], pos2[2], a, Ecore);

            home->param->TotalEnergy += W;
#endif


/*
 *          PK force from external stress
 */
            ExtPKForce(extstress,
                       burg[0], burg[1], burg[2],
                       pos1[0], pos1[1], pos1[2],
                       pos2[0], pos2[1], pos2[2], f1, f2);

            AddtoArmForce(node1, armID12, f1);
            AddtoArmForce(node2, armID21, f2);

#ifdef CALCENERGY
            W = ExtPKEnergy(extstress, E, NU,
                            burg[0], burg[1], burg[2],
                            pos1[0], pos1[1], pos1[2],
                            pos2[0], pos2[1], pos2[2]);

            home->param->TotalEnergy += W;
#endif

/*
 *          If we're including osmotic forces, add those in now
 */
            if (param->vacancyConcEquilibrium > 0.0)
            {
                OsmoticForce(home,
                             pos1[0], pos1[1], pos1[2],
                             pos2[0], pos2[1], pos2[2],
                             burg[0], burg[1], burg[2], f1, f2);

                AddtoArmForce(node1, armID12, f1);
                AddtoArmForce(node2, armID21, f2);
            }

#ifdef _ARLFEM
/*
 *          With free surfaces, we also have to calculate the force
 *          from all "virtual" segments in the simulation space.
 *          Virtual segments are imaginary segments extending
 *          outward from any segment intersecting the surface.
 *          The segments are treated as if they extend from
 *          the surface to infinity, but are in reality long
 *          but finite segments.
 */
            VirtualSegForce(home, pos1, pos2, burg, f1, f2);

            AddtoArmForce(node1, armID12, f1);
            AddtoArmForce(node2, armID21, f2);
#endif  /* ifdef _ARLFEM */


#ifdef ESHELBY
/*
 *          If this function is called for two connected nodes, we have to
 *          be sure we don't search for intersecting particles multiple
 *          times, so only add the segment/particle intersections if the
 *          segment is not yet on the list.  (Note: list may be unsorted at
 *          this point, so we have to use a lookup function that handles
 *          an unsorted list)
 */
            SegPartIntersect_t *intersection = SegPartListUnsortedLookup(home, &node1->myTag, &node2->myTag);
            int                 addToList    = (intersection == (SegPartIntersect_t *)NULL ? 1 : 0);

            for (j = 0; j < home->totInclusionCount; j++)
            {
               EInclusion_t *inclusion = &home->eshelbyInclusions[j];

               int ninc = 0;
               real8 newPos1[3], newPos2[3];
               real8 ratio[2];
               isIntersecting = IntersectionInformation(home, pos1, pos2, inclusion,
                                                        newPos1, newPos2, &ninc,
                                                        ratio);

#ifdef ESHELBYFORCE
               EshelbyForce(home, inclusion->position, inclusion->radius,
                            inclusion->strainField, newPos1, newPos2,
                            burg, f1, f2);

               AddtoArmForce(node1, armID12, f1);
               AddtoArmForce(node2, armID21, f2);
#endif  // ESHELBYFORCE

               if (isIntersecting)
               {

#ifdef ESHELBYCORE
                  SegPartCoreForce(home,inclusion->position, inclusion->radius,
                                   inclusion->rotation, newPos1, newPos2,
                                   burg, f1, f2);


                  AddtoArmForce(node1, armID12, f1);
                  AddtoArmForce(node2, armID21, f2);
#endif //ESHELBYCORE

#ifdef ESHELBYSTEP
                  InclusionStepForce(home->param->IncStepEnergy, inclusion->position, inclusion->radius,
                                     inclusion->rotation, newPos1, newPos2,
                                     burg, f1, f2);


                  AddtoArmForce(node1, armID12, f1);
                  AddtoArmForce(node2, armID21, f2);
#endif //ESHELBYSTEP

/*
 *              If the segment intersects the inclusion and we
 *              need to track such things for the mobility function,
 *              add the segment segment/particle intersection
 *              to the list.
 */
               if (addToList)
                  SegPartListUpdate(home, &intersection, &node1->myTag, &node2->myTag,
                                    j, inclusion->id);
               }
            }

#endif  // ESHELBY

/*
 *          Now loop through all the segments and do any necessary
 *          seg/seg force calcs.
 */
            for (k = 0; k < home->newNodeKeyPtr; k++)
            {
                if ((node3 = home->nodeKeys[k]) == (Node_t *)NULL)
                {
                    continue;
                }

                for (l = 0; l < node3->numNbrs; l++)
                {
                    node4 = GetNeighborNode(home, node3, l);

                    if (node4 == (Node_t *)NULL)
                    {
                        printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                        continue;
                    }

/*
 *                  Insures the node with the lower tag is the node3
 */
                    if (OrderNodes(node3, node4) >= 0) { continue; }

/*
 *                  Make sure we don't try to calculate seg/seg forces
 *                  on a segment with itself.
 */
                    if ((node1==node3) && (node2==node4)) { continue; }

/*
 *                  And after computing forces betweeen the two segments,
 *                  only update forces on segment (node1--node2).
 */
                    ComputeForces(home, node1,node2,node3,node4, f1,f2,f3,f4);

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);

#ifdef CALCENERGY
                    real8 W = 0.0;
                    W = ComputeEnergy(home, node1, node2, node3, node4);
                    home->param->TotalEnergy += W;
#endif

                }  /* for (l = 0; ...) */

            }  /* for (k = 0; ... ) */

/*
 *          We're not yet done modifying forces for the primary
 *          node, but we are done fiddling with the forces on the
 *          neighbor node (nodeB), so reset the total nodal force
 *          for the neighbor node to the sum of it's segment forces.
 */
            nodeB->fX = 0.0;
            nodeB->fY = 0.0;
            nodeB->fZ = 0.0;

            for (i = 0; i < nodeB->numNbrs; i++)
            {
                nodeB->fX += nodeB->armfx[i];
                nodeB->fY += nodeB->armfy[i];
                nodeB->fZ += nodeB->armfz[i];
            }
        }  /* for (armIndex = 0; ... ) */

/*
 *      All segment specific forces for the primary node (nodeA) have
 *      been updated, so now reset the total nodal force to the sum
 *      of the node's segment forces.
 */
        nodeA->fX = 0.0;
        nodeA->fY = 0.0;
        nodeA->fZ = 0.0;

        for (i = 0; i < nodeA->numNbrs; i++)
        {
            nodeA->fX += nodeA->armfx[i];
            nodeA->fY += nodeA->armfy[i];
            nodeA->fZ += nodeA->armfz[i];
        }
        return;
}
#else  /* FULL_N2_FORCES not defined */
/*
 *      WARNING! The FEM code reads and writes to the DSM (Distributed
 *               Shared Memory) which is a collective operation, however,
 *               SetOneNodeForce() is called by individual tasks, hence
 *               the FEM cannot currently be coupled with SetOneNodeForce().
 */
void SetOneNodeForce(Home_t *home, Node_t *nodeA)
{
        real8   sigb[3], totRemSig[3][3], extstress[3][3];
        real8   fseg1node1[3], fseg1node2[3];
        Node_t  *nodeB;

        Param_t   *param       = home->param;
        SegPair_t *segPairList = 0;

        real8 a     = param->rc;
        real8 MU    = param->shearModulus;
        real8 NU    = param->pois;
        real8 E     = param->YoungsModulus;
        real8 Ecore = param->Ecore;
        real8 eps   = 1.0e-6;

        extstress[0][0] = param->appliedStress[0];
        extstress[1][1] = param->appliedStress[1];
        extstress[2][2] = param->appliedStress[2];
        extstress[1][2] = param->appliedStress[3];
        extstress[2][0] = param->appliedStress[4];
        extstress[0][1] = param->appliedStress[5];
        extstress[2][1] = extstress[1][2];
        extstress[0][2] = extstress[2][0];
        extstress[1][0] = extstress[0][1];

        int  segPairCnt        = 0;
        int  segPairListSize   = 0;

#ifdef ESHELBY
        int  inclusionCnt      = 0;
        int  inclusionListSize = 0;
        int *inclusionList     = 0;
#endif

#ifdef SPECTRAL
        Spectral_t *spectral = home->spectral;
        double      rnei2    = 0.0;
        double      dist2    = 0.0;

        if (param->FFTenabled) { rnei2 = spectral->rcnei*spectral->rcnei; }
#endif

/*
 *      The following block of code is not yet supported!
 */
#if 0
#ifdef _ARLFEM
        std::vector<double> segmentData;

/*
 *      construct list of segments to send to FEM
 */
        for (int armIndex = 0; armIndex < nodeA->numNbrs; armIndex++)
        {
            nodeB = GetNodeFromTag(home, nodeA->nbrTag[armIndex]);

            if (nodeB == (Node_t *)NULL) { continue; }

            int armIDAB = GetArmID(nodeA, nodeB);

            segmentData.push_back(0.0); // surface node status does not matter
            segmentData.push_back(nodeA->x);
            segmentData.push_back(nodeA->y);
            segmentData.push_back(nodeA->z);
            segmentData.push_back(nodeB->x);
            segmentData.push_back(nodeB->y);
            segmentData.push_back(nodeB->z);
            segmentData.push_back(nodeA->burgX[armIDAB]);
            segmentData.push_back(nodeA->burgY[armIDAB]);
            segmentData.push_back(nodeA->burgZ[armIDAB]);

        }

/*
 *      Open hdf5 dsm buffer
 */
        hid_t hdf5FileHandle = H5Fcreate("dsm",
                                         H5F_ACC_TRUNC,
                                         H5P_DEFAULT,
                                         home->dsmFAPL);

/*
 *      Write segment data to buffer
 */
        DSMCollectiveWrite(&(segmentData[0]),
                           segmentData.size(),
                           "segment data",
                           hdf5FileHandle);

/*
 *      Send what type of force calculation, PARTIAL
 */
        DSMSingleTaskWrite(PARTIAL,
                           "force calculation type",
                           hdf5FileHandle);

/*
 *      Release control to FEM
 */
        herr_t hdf5Status = H5Fclose(hdf5FileHandle);
        std::cout << "Closed H5 file in SetOneNodeForce1 with status: "
                  << hdf5Status << std::endl;

        H5FD_dsm_server_update(home->dsmManager->GetDSMHandle());

#endif /* ifdef _ARLFEM */
#endif /* 0 */

/*
 *      Loop over all the segments attached to the node.  Build a list
 *      of all the segment pairs for which forces must be calculated.
 *
 *      Note: Since we're only resetting forces for segments attached
 *      to the specified node, only the forces on the first segment
 *      in each pair will be reset.
 */
        for (int armIndex = 0; armIndex < nodeA->numNbrs; armIndex++)
        {
            int     cx, cy, cz;
            int     cellX, cellY, cellZ;
            int     minXIndex, minYIndex, minZIndex;
            int     maxXIndex, maxYIndex, maxZIndex;
            real8   dx, dy, dz;
            real8   x1, y1, z1;
            real8   x2, y2, z2;
            real8   xCenter, yCenter, zCenter;
            real8   burg[3], f1[3], f2[3];
            Cell_t  *cell;
            Node_t  *node3, *node4;

            nodeB = GetNodeFromTag(home, nodeA->nbrTag[armIndex]);

            if (nodeB == (Node_t *)NULL) { continue; }

/*
 *          Always point node1 to the node to be treated as the segment
 *          owner.
 *
 *          The rationale behind this is: if the segment is owned by a
 *          node in a cell not native to the local domain, the domain
 *          probably does not have the FMM coefficients for the cell
 *          owning the segment and would be unable to properly calculate
 *          remote forces  Therefore, in this routine, we ignore the
 *          normal segment ownership rules and treat any segment attached
 *          to nodeA as if it is owned by nodeA.
 */
            Node_t  *node1   = nodeA;
            Node_t  *node2   = nodeB;
            int      armID12 = armIndex;
            int      armID21 = GetArmID(node2, node1);

/*
 *          If we're using the FMM code, and the node owning the segment is
 *          not local, we probably do not have the taylor expansion
 *          coefficients for the xcell owning the segment, and cannot
 *          calculate new remote forces, so leave the forces on this segment
 *          untouched.
 */
#ifdef ESHELBY
            if (((param->fmEnabled) || (param->eshelbyfmEnabled)) && (node1->myTag.domainID != home->myDomain))
                continue;
#else
            if (  param->fmEnabled                                && (node1->myTag.domainID != home->myDomain))
                continue;
#endif

/*
 *          We need the indices of the cell containing node1 in several
 *          places further on, so get that now.  If for some reason, the
 *          node has not yet been associated with a cell, we can't
 *          calculate all the forces, so skip it.
 */
            if (node1->cellIdx < 0) continue;

            cell = LookupCell(home, node1->cellIdx);

            cellX = cell->xIndex;
            cellY = cell->yIndex;
            cellZ = cell->zIndex;

/*
 *          We should be able to calculate new forces for the segment, so
 *          zero out the old segment forces.
 */
            node1->armfx[armID12] = 0.0;
            node1->armfy[armID12] = 0.0;
            node1->armfz[armID12] = 0.0;

            node2->armfx[armID21] = 0.0;
            node2->armfy[armID21] = 0.0;
            node2->armfz[armID21] = 0.0;

            V3_ZERO(fseg1node1);
            V3_ZERO(fseg1node2);

/*
 *          If we need to add in any FEM image stress, or we're not using
 *          the fmm code,  we'll need to explicitly zero out the segment
 *          sigbRem.
 */
#ifndef _ARLFEM
            if (param->fmEnabled == 0)
#endif
            {
                V3_ZERO(node1->sigbRem+(armID12*3));
                V3_ZERO(node2->sigbRem+(armID21*3));
            }

/*
 *          If the segment is zero-length, we cannot calculate forces on
 *          it, so skip it.
 */
            x1 = node1->x;
            y1 = node1->y;
            z1 = node1->z;

            dx = node2->x - x1;
            dy = node2->y - y1;
            dz = node2->z - z1;

            ZImage(param, &dx, &dy, &dz);

            x2 = x1 + dx;
            y2 = y1 + dy;
            z2 = z1 + dz;

            burg[0] = node1->burgX[armID12];
            burg[1] = node1->burgY[armID12];
            burg[2] = node1->burgZ[armID12];

            if (((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) < eps) { continue; }

/*
 *          If elastic interaction is not enabled, use a simple line
 *          tension model for calculating forces on the segment. (i.e.
 *          no segment/segment forces, no remote forces, etc)
 *
 *          Segment forcess are updated here; the nodal forces are
 *          resummed at the end of this routine.
 */
            if (!param->elasticinteraction)
            {
                LineTensionForce(home,
                                 node1->x, node1->y, node1->z,
                                 node2->x, node2->y, node2->z,
                                 burg[0], burg[1], burg[2], f1, f2);

                AddtoArmForce(node1, armID12, f1);
                AddtoArmForce(node2, armID21, f2);

#ifdef CALCENERGY
                real8 W = LineTensionEnergy(home,
                                            node1->x, node1->y, node1->z,
                                            node2->x, node2->y, node2->z,
                                            burg[0], burg[1], burg[2]);

                home->param->TotalEnergy += W;
#endif

                continue;
            }

/*
 *          Add the various forces on this segment which are not due to
 *          direct interactions with other segments.  Start with force
 *          due to self stress.
 */
            SelfForce(home, 0, MU, NU, burg[0], burg[1], burg[2], x1,y1,z1, x2,y2,z2, a,Ecore, f1,f2);

            V3_ACCUM(fseg1node1, f1);
            V3_ACCUM(fseg1node2, f2);

#ifdef CALCENERGY
            real8 W = SelfEnergy(0, MU, NU, burg[0], burg[1], burg[2], x1,y1,z1, x2,y2,z2, a, Ecore);

            home->param->TotalEnergy += W;
#endif
/*
 *          Now PK force from external stress
 */
            ExtPKForce(extstress, burg[0], burg[1], burg[2], x1,y1,z1, x2,y2,z2, f1,f2);

            V3_ACCUM(fseg1node1, f1);
            V3_ACCUM(fseg1node2, f2);

#ifdef CALCENERGY
            W = ExtPKEnergy(extstress, E, NU, burg[0], burg[1], burg[2], x1,y1,z1, x2,y2,z2);

            home->param->TotalEnergy += W;
#endif

/*
 *          If we're including osmotic forces, add those in
 */
            if (param->vacancyConcEquilibrium > 0.0)
            {
                OsmoticForce(home,
                             x1, y1, z1, x2, y2, z2,
                             burg[0], burg[1], burg[2], f1, f2);

                V3_ACCUM(fseg1node1, f1);
                V3_ACCUM(fseg1node2, f2);
            }

#ifdef ESHELBY
/*
 *          If the simulation includes Eshelby inclusions, we need to add
 *          the force due to the inclusions
 */
            if (param->enableInclusions)
            {
               int    k, numNbrCells;
               real8  pos1[3], pos2[3];

/*
 *              If the inclusion list was set in a previous iteration of the
 *              loop over the arms of the primary node, just zero out the
 *              inclusion count so we reuse the previously allocated array.
 */
               inclusionCnt = 0;

               pos1[0] = x1;
               pos1[1] = y1;
               pos1[2] = z1;

               pos2[0] = x2;
               pos2[1] = y2;
               pos2[2] = z2;

               numNbrCells = cell->nbrCount;

/*
 *              If this function is called for two connected nodes, we have to
 *              be sure we don't search for intersecting particles multiple
 *              times, so only add the segment/particle intersections if the
 *              segment is not yet on the list.  (Note: list may be unsorted at
 *              this point, so we have to use a lookup function that handles
 *              an unsorted list)
 */
               SegPartIntersect_t *intersection = SegPartListUnsortedLookup(home, &node1->myTag, &node2->myTag);
               int                 addToList    = (intersection == (SegPartIntersect_t *)NULL ? 1 : 0 );

               for (k = 0; k <= numNbrCells; k++)
               {
                  int inclusionCellID=0;
/*
 *                  First loop over all neighboring cells, then on
 *                  last loop iteration do the current cell.
 */
                  if (k < numNbrCells) { inclusionCellID = cell->nbrList[k]; }
                  else                 { inclusionCellID = node1->cellIdx;   }

/*
 *                  If the cell is a ghost cell, alter the cell pointer to
 *                  the corresponding base cell.  Needed in order to pick
 *                  up the proper inclusionQ.
 */
                  Cell_t *inclusionCell = LookupCell(home, inclusionCellID);

                  if (inclusionCell->baseIdx >= 0)
                  {
                     int inclusionCellID;
                     inclusionCellID = inclusionCell->baseIdx;
                     inclusionCell = LookupCell(home, inclusionCellID);
                  }

/*
 *                  Add the inclusion indices for all inclusions in this
 *                  cell to the list of inclusions for which we have to
 *                  calculate interaction forces for the current segment
 */
                  int inclusionIndex = inclusionCell->inclusionQ;

                  while (inclusionIndex >= 0)
                  {
                     EInclusion_t *inclusion;

                     inclusion = &home->eshelbyInclusions[inclusionIndex];

                     if (inclusionCnt == inclusionListSize)
                     {
                        inclusionListSize += 1000;
                        inclusionList = (int *)realloc(inclusionList,
                                                       sizeof(int) * inclusionListSize);
                     }

                     inclusionList[inclusionCnt] = inclusionIndex;
                     inclusionIndex = inclusion->nextInCell;

                     inclusionCnt++;
                  }

               }  /* end for (k = 0; k <= numNbrCells; ...) */

/*
 *              Now we have the list of all inclusions which directly
 *              interact with the node1/node2 segment, so calculate the
 *              force on the segment from each inclusion.
 */


#ifdef _OPENMP
#pragma omp parallel
#endif
               {
                  int   m, threadID, threadIterStart, threadIterEnd;
                  real8 f1Thread[3], f2Thread[3];

                  V3_ZERO(f1Thread);
                  V3_ZERO(f2Thread);

                  GetThreadIterationIndices(inclusionCnt, &threadID, &threadIterStart, &threadIterEnd);

/*
 *                  Loop through all the inclusions in the cell
 *                  and compute the force on the segment from
 *                  each inclusion.
 */
                  for (m = threadIterStart; m < threadIterEnd; m++)
                  {
                     int  isIntersecting;
                     EInclusion_t *inclusion;

                     inclusion = &home->eshelbyInclusions[inclusionList[m]];

                     int ninc=0;
                     real8 newPos1[3], newPos2[3];
                     real8 ratio[2];
                     isIntersecting = IntersectionInformation(home,pos1,pos2,inclusion,
                                                              newPos1, newPos2, &ninc,
                                                              ratio);

#ifdef ESHELBYFORCE
                     // Eshelby force applies on/effects all segments in the simulation and not on
                     // only the segments intersecting the particles!
                     real8 f1EshelbyForce[3], f2EshelbyForce[3];
                     EshelbyForce(home, inclusion->position, inclusion->radius,
                                  inclusion->strainField, newPos1, newPos2,
                                  burg, f1EshelbyForce, f2EshelbyForce);

                     V3_ACCUM(f1Thread, f1EshelbyForce);
                     V3_ACCUM(f2Thread, f2EshelbyForce);
#endif   //ESHELBYFORCE



                     if (isIntersecting)
                     {
#ifdef ESHELBYCORE
                        real8 f1EshelbyCore[3], f2EshelbyCore[3];
                        SegPartCoreForce(home, inclusion->position, inclusion->radius,
                                         inclusion->rotation, newPos1, newPos2,
                                         burg, f1EshelbyCore, f2EshelbyCore);


                        V3_ACCUM(f1Thread, f1EshelbyCore);
                        V3_ACCUM(f2Thread, f2EshelbyCore);
#endif //ESHELBYCORE


#ifdef ESHELBYSTEP
                        real8 f1EshelbyStep[3], f2EshelbyStep[3];
                        InclusionStepForce(home->param->IncStepEnergy, inclusion->position, inclusion->radius,
                                           inclusion->rotation, newPos1, newPos2,
                                           burg, f1EshelbyStep, f2EshelbyStep);


                        V3_ACCUM(f1Thread, f1EshelbyStep);
                        V3_ACCUM(f2Thread, f2EshelbyStep);
#endif //ESHELBYSTEP

/*
 *                      If the segment intersects the inclusion and we
 *                      need to track such things for the mobility function,
 *                      add the segment segment/particle intersection
 *                      list.
 */
                        if (addToList)
                        {
                            SegPartListUpdate(home, &intersection,
                                              &node1->myTag, &node2->myTag,
                                              inclusionList[m], inclusion->id);
                        }
                     }
                  }

/*
 *                  Add the eshelby force from each thread to the segment
 *                  forces.  Necessary locking is handled in AddtoArmForce().
 *                  f1 and f2Thread are initialized to zero. They are on-zero
 *                  when either ESHELBYFORCE or/and ESHELBYCORE are on.
 */
                  AddtoArmForce(node1, armID12, f1Thread);
                  AddtoArmForce(node2, armID21, f2Thread);

               }  /* end omp parallel section */
            }  /* if (param->enableInclusions) */
#endif  // ESHELBY

/*
 *          Last thing for this arm of the target node is to add force
 *          component from remote stress
 */
            V3_ZERO(f1);
            V3_ZERO(f2);

#ifdef ESHELBY
            if ((param->fmEnabled) || (param->eshelbyfmEnabled))
#else
            if  (param->fmEnabled)
#endif
            {

                RemoteForceOneSeg(home, node1, node2, f1, f2);

                V3_ACCUM(fseg1node1, f1);
                V3_ACCUM(fseg1node2, f2);
            }

#ifdef SPECTRAL
            if (param->FFTenabled)
            {
                SpectralForceOneSeg(home, node1, node2, f1, f2);

                V3_ACCUM(fseg1node1, f1);
                V3_ACCUM(fseg1node1, f2);
            }
#endif

/*
 *          If we're not using the FMM, find the segment midpoint and use
 *          that as the point at which to calculate the stress from remote
 *          segments.
 */
#ifdef SPECTRAL
            if ( (param->fmEnabled==0) && (param->FFTenabled==0) )
#else
            if   (param->fmEnabled==0)
#endif
            {
                real8 xm, ym, zm;

                xm = x1 + (dx * 0.5);
                ym = y1 + (dy * 0.5);
                zm = z1 + (dz * 0.5);

/*
 *              Indices for the cell containing node1 were obtained
 *              way up above... so we just need to subtract off the
 *              adjustment for periodic (ghost) cells and find the stress
 *              at the specified point from all segments in remote cells.
 */
                GetFieldPointStressRem(home, xm,ym,zm, cellX-1, cellY-1, cellZ-1, totRemSig);

                sigb[0] = totRemSig[0][0] * burg[0] +
                          totRemSig[0][1] * burg[1] +
                          totRemSig[0][2] * burg[2];

                sigb[1] = totRemSig[1][0] * burg[0] +
                          totRemSig[1][1] * burg[1] +
                          totRemSig[1][2] * burg[2];

                sigb[2] = totRemSig[2][0] * burg[0] +
                          totRemSig[2][1] * burg[1] +
                          totRemSig[2][2] * burg[2];

                V3_COPY(node1->sigbRem+(armID12*3), sigb);
                V3_COPY(node2->sigbRem+(armID21*3), sigb);
#ifndef _ARLFEM
/*
 *              If we are not linked with FEM, calculate the PK force
 *              from remote segments and add it in.
 */
                PKForce(sigb, x1, y1, z1, x2, y2, z2, f1, f2);

                V3_ACCUM(fseg1node1, f1);
                V3_ACCUM(fseg1node2, f2);
#endif


#if defined(_ARLFEM) && !defined(_ARLFEM_FULL_IMAGE_FORCES)
/* When the hack is on, need PKForces. */
/*
 *              If we are not linked with FEM, calculate the PK force
 *              from remote segments and add it in.
 */
                PKForce(sigb, x1, y1, z1, x2, y2, z2, f1, f2);

                V3_ACCUM(fseg1node1, f1);
                V3_ACCUM(fseg1node2, f2);

#endif
            }  /* if not fmEnabled */

#ifdef _ARLFEM
/*
 *          With free surfaces, we also have to calculate the force from all
 *          "virtual" segments in the simulation space.  Virtual segments are
 *          imaginary segments extending outward from any segment intersecting
 *          the surface.  The segments are treated as if they extend from the
 *          surface to infinity, but are in reality long but finite segments.
 */
            {
                real8 pos1[3] = { x1,y1,z1 };
                real8 pos2[3] = { x2,y2,z2 };

                VirtualSegForce(home, pos1, pos2, burg, f1, f2);

                V3_ACCUM(fseg1node1, f1);
                V3_ACCUM(fseg1node2, f2);
            }

#endif  /* ifdef _ARLFEM */

/*
 *          Add the forces accumulated so far into the segment's node forces
 */
            AddtoArmForce(node1, armID12, fseg1node1);
            AddtoArmForce(node2, armID21, fseg1node2);

/*
 *          We found the indices of the cell owning the first node
 *          earlier, but we also need to know the locations of the
 *          center of that cell.
 */
            xCenter = cell->center[0];
            yCenter = cell->center[1];
            zCenter = cell->center[2];

/*
 *          Now we need to find all the other segments for which force
 *          interactions must be calculated with segment node1/node2
 *          Need to compute the segment-to-segment force on this
 *          segment from all other segments in the same cell and all
 *          immediately neighboring cells.  So, we need to find the
 *          minimum and maximum cell indices of that block of cells
 *          allowing for periodic boundaries.
 */
            if (param->xBoundType == Periodic) { minXIndex = MAX(0, cellX-1); maxXIndex = MIN(param->nXcells+1, cellX+1); }
            else                               { minXIndex = MAX(1, cellX-1); maxXIndex = MIN(param->nXcells  , cellX+1); }

            if (param->yBoundType == Periodic) { minYIndex = MAX(0, cellY-1); maxYIndex = MIN(param->nYcells+1, cellY+1); }
            else                               { minYIndex = MAX(1, cellY-1); maxYIndex = MIN(param->nYcells  , cellY+1); }

            if (param->zBoundType == Periodic) { minZIndex = MAX(0, cellZ-1); maxZIndex = MIN(param->nZcells+1, cellZ+1); }
            else                               { minZIndex = MAX(1, cellZ-1); maxZIndex = MIN(param->nZcells  , cellZ+1); }

#ifdef SPECTRAL
            if ( param->FFTenabled && (rnei2 == 0.0) )  { continue; }
#endif

/*
 *          Now loop through all cells in that block.
 */
            for (cx = minXIndex; cx <= maxXIndex; cx++)
            {
                for (cy = minYIndex; cy <= maxYIndex; cy++)
                {
                    for (cz = minZIndex; cz <= maxZIndex; cz++)
                    {
                        int cellIndex = EncodeCellIdx(home, cx, cy, cz);
                        cell = LookupCell(home, cellIndex);

                        if (cell == (Cell_t *)NULL) continue;

/*
 *                      Loop though all nodes in this cell, and each arm
 *                      of each node.
 */
                        node3 = cell->nodeQ;

                        for (; node3!=(Node_t *)NULL; node3=node3->nextInCell)
                        {
                            for (int armID34=0; armID34<node3->numNbrs; armID34++)
                            {

                                node4 = GetNeighborNode(home, node3, armID34);

/*
 *                              Skip the node3/node4 segment if:
 *                              1) the current domain has no info on node4
 *                              2) node3 does not own the node3/node4 segment
 *                              3) the node3/node4 segment is the same as
 *                                 the node1/node2 segment.
 */
                                if ((node4 == (Node_t *)NULL)               ||
                                    (NodeOwnsSeg(home, node4, node3))       ||
                                    (((node1 == node3) && (node2 == node4)) ||
                                     ((node1 == node4) && (node2 == node3))))
                                {
                                    continue;
                                }

/*
 *                              Add segment pair to list for which direct
 *                              seg/seg forces need to be calculated.
 */
                                if (segPairCnt == segPairListSize)
                                {
                                    segPairListSize += 1000;
                                    segPairList = (SegPair_t *)realloc( segPairList, sizeof(SegPair_t) * segPairListSize);
                                }

                                segPairList[segPairCnt].node1 = node1;
                                segPairList[segPairCnt].node2 = node2;
                                segPairList[segPairCnt].node3 = node3;
                                segPairList[segPairCnt].node4 = node4;

                                segPairList[segPairCnt].cellCenter[0] = xCenter;
                                segPairList[segPairCnt].cellCenter[1] = yCenter;
                                segPairList[segPairCnt].cellCenter[2] = zCenter;

                                segPairCnt++;

                            }  /* for (armID34 = 0; ... ) */
                        }  /* loop over cell nodeQ */
                    }  /* for (cz = minZIndex; ... ) */
                }  /* for (cy = minYIndex; ...) */
            }  /* for (cx = minXIndex; ...) */
        }  /* for (armIndex = 0; ...) */

/*
 *      The following block of code is not yet functional
 */
#if 0
#ifdef _ARLFEM
/*
 *      Open distributed shared memory file and read node forces.
 */
        hdf5FileHandle = H5Fopen("dsm",
                                 H5F_ACC_RDONLY,
                                 home->dsmFAPL);

        int segmentForcesSize = (segmentData.size() * 3) / 5;
        real8 * segmentForces =
          (real8 *)malloc(segmentForcesSize * sizeof(real8));

        DSMCollectiveRead("segment forces",
                          hdf5FileHandle,
                          segmentForcesSize,
                          segmentForces);

        hdf5Status = H5Fclose(hdf5FileHandle);
        std::cout << "Closed H5 file in SetOneNodeForce2 with status: "
                  << hdf5Status << std::endl;

        // add image force contribution to nodes

        real8 * segmentForcesPtr = segmentForces;
        for (int armIndex = 0; armIndex < nodeA->numNbrs; armIndex++)
        {
            nodeB = GetNodeFromTag(home, nodeA->nbrTag[armIndex]);

            if (nodeB == (Node_t *)NULL) { continue; }

            int armIDAB = GetArmID(nodeA, nodeB);
            int armIDBA = GetArmID(nodeB, nodeA);

            AddtoArmForce(nodeA, armIDAB, segmentForcesPtr);
            AddtoArmForce(nodeB, armIDBA, segmentForcesPtr + 3);

            segmentForcesPtr += 6;
        }

        free(segmentForces);

#endif /* ifdef _ARLFEM */
#endif /* 0 */

/*
 *      Calculate the forces between every segment pair in the list
 *      created above, and update the forces on segment 1.
 */
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int     threadID=0, threadIterStart=0, threadIterEnd=0;
            GetThreadIterationIndices(segPairCnt, &threadID, &threadIterStart, &threadIterEnd);

            for (int pairID = threadIterStart; pairID < threadIterEnd; pairID++)
            {
                Node_t *node1 = segPairList[pairID].node1;
                Node_t *node2 = segPairList[pairID].node2;
                Node_t *node3 = segPairList[pairID].node3;
                Node_t *node4 = segPairList[pairID].node4;

#ifdef SPECTRAL
                if (param->FFTenabled)
                {
                    if (rnei2 == 0.0) continue;
                    MinSegSegDist(home, node1, node2, node3, node4, &dist2);
                    if ( (dist2 < 0) || (dist2 > rnei2) ) continue;
                }
#endif

                int armID12 = GetArmID(node1, node2);
                int armID21 = GetArmID(node2, node1);
                int armID34 = GetArmID(node3, node4);

                real8 b12[3]  = { node1->burgX[armID12],
                                  node1->burgY[armID12],
                                  node1->burgZ[armID12] };

                real8 b34[3]  = { node3->burgX[armID34],
                                  node3->burgY[armID34],
                                  node3->burgZ[armID34] };

/*
 *              Find coordinates of image of the third node closest to the
 *              center of the cell containing the first node, and use the
 *              image of the fourth node closest to that point.  We need
 *              use the cell center as the base because all the Burgers
 *              vectors of segments inside one cell (belonging to the same
 *              PBC image) are summed together for 'remote' stress
 *              calculations, so cannot separate them for 'local' stress
 *              calculations.
 */
                real8 p1[3] = { node1->x, node1->y, node1->z };        // p1 = position of node 1 <xyz>
                real8 p2[3] = { node2->x, node2->y, node2->z };        // p2 = position of node 2 <xyz>
                real8 p3[3] = { node3->x, node3->y, node3->z };        // p3 = position of node 3 <xyz>
                real8 p4[3] = { node4->x, node4->y, node4->z };        // p4 = position of node 4 <xyz>

                real8 dx = (p2[0]-p1[0]);
                real8 dy = (p2[1]-p1[1]);
                real8 dz = (p2[2]-p1[2]);

                ZImage(param,dx,dy,dz);

                p2[0] = (p1[0]+dx);
                p2[1] = (p1[1]+dy);
                p2[2] = (p1[2]+dz);

                real8 pc[3]   = { segPairList[pairID].cellCenter[0],      // pc = cell center <xyz>
                                  segPairList[pairID].cellCenter[1],      //
                                  segPairList[pairID].cellCenter[2] };    //

                PBCPOSITION(param, pc, p3);
                PBCPOSITION(param, p3, p4);

/*
 *              If either segment is zero-length, skip it
 */
                dx = (p1[0]-p2[0]);
                dy = (p1[1]-p2[1]);
                dz = (p1[2]-p2[2]);

                if ( ((dx*dx)+(dy*dy)+(dz*dz)) < eps ) { continue; }

                dx = (p3[0]-p4[0]);
                dy = (p3[1]-p4[1]);
                dz = (p3[2]-p4[2]);

                if ( ((dx*dx)+(dy*dy)+(dz*dz)) < eps ) { continue; }

/*
 *              Compute the forces on the 2 segments and add the force
 *              contribution from segment node3/node4 to segment node1/node2.
 */

                real8 f1[3] = { 0.0 };   // f1 = resulting force on node 1
                real8 f2[3] = { 0.0 };   // f2 = resulting force on node 2
                real8 f3[3] = { 0.0 };   // f3 = resulting force on node 3
                real8 f4[3] = { 0.0 };   // f4 = resulting force on node 4

                SegSegForce(home, f1,f2,f3,f4, p1,p2,p3,p4, b12,b34, a,MU,NU);

#ifdef SPECTRAL
                if (param->FFTenabled)
                {
                    real8  f1c[3] = { 0.0 };  // f1c = spectral correction force on node 1
                    real8  f2c[3] = { 0.0 };  // f2c = spectral correction force on node 2
                    real8  f3c[3] = { 0.0 };  // f3c = spectral correction force on node 3
                    real8  f4c[3] = { 0.0 };  // f4c = spectral correction force on node 4

                    SegSegForce(home, f1c,f2c,f3c,f4c, p1,p2,p3,p4, b12,b34, spectral->rcgrid,MU,NU);

                    V3_SUB(f1,f1,f1c);        // apply spectral correction
                    V3_SUB(f2,f2,f2c);        //
                    V3_SUB(f3,f3,f3c);        //
                    V3_SUB(f4,f4,f4c);        //
                }
#endif

/*
 *              Increment the accumulated forces on the node1/node2 segment
 */

                AddtoArmForce(node1, armID12, f1);
                AddtoArmForce(node2, armID21, f2);

#ifdef CALCENERGY
                real8 W = SegSegEnergy( p1[0],p1[1],p1[2],
                                        p2[0],p2[1],p2[2],
                                        p3[0],p3[1],p3[2],
                                        p4[0],p4[1],p4[2],
                                        b12[0], b12[1], b12[2],
                                        b34[0], b34[1], b34[2],
                                        a, MU, NU );

                home->param->TotalEnergy += W;
#endif

            }  /* end loop over segment pairs */

        }  /* end omp parallel section */

/*
 *      All segment specific forces for the primary node (nodeA) and
 *      all of its neighbors have been updated, so now reset the total
 *      nodal force for each of those nodes to the sum of the nodes'
 *      respective segment forces.
 */
        nodeA->fX = 0.0;
        nodeA->fY = 0.0;
        nodeA->fZ = 0.0;

        for (int armIndex = 0; armIndex < nodeA->numNbrs; armIndex++)
        {
            int nbrArmIndex;

            nodeA->fX += nodeA->armfx[armIndex];
            nodeA->fY += nodeA->armfy[armIndex];
            nodeA->fZ += nodeA->armfz[armIndex];

            nodeB = GetNodeFromTag(home, nodeA->nbrTag[armIndex]);

            nodeB->fX = 0.0;
            nodeB->fY = 0.0;
            nodeB->fZ = 0.0;

            for (nbrArmIndex = 0; nbrArmIndex < nodeB->numNbrs; nbrArmIndex++)
            {
                nodeB->fX += nodeB->armfx[nbrArmIndex];
                nodeB->fY += nodeB->armfy[nbrArmIndex];
                nodeB->fZ += nodeB->armfz[nbrArmIndex];
            }
        }

        if (segPairList  ) { free(segPairList  ); }

#ifdef ESHELBY
        if (inclusionList) { free(inclusionList); }
#endif

        return;
}
#endif  /* FULL_N2_FORCES not defined */


void ExtPKForce(real8 str[3][3],
                real8 bx, real8 by, real8 bz,
                real8 x1, real8 y1, real8 z1,
                real8 x2, real8 y2, real8 z2,
                real8 f1[3], real8 f2[3])
{
        real8 strb[3], xi[3], ft[3];
        int j;

        xi[0] = x2-x1;
        xi[1] = y2-y1;
        xi[2] = z2-z1;

        for (j = 0; j < 3; j++)
        {
            strb[j] = str[j][0]*bx + str[j][1]*by + str[j][2]*bz;
        }

        cross(strb, xi, ft);

        for (j = 0; j < 3; j++)
        {
            f1[j] = ft[j]*0.5;
            f2[j] = ft[j]*0.5;
        }
}


void PKForce(real8 sigb[3],
             real8 x1, real8 y1, real8 z1,
             real8 x2, real8 y2, real8 z2,
             real8 f1[3], real8 f2[3])
{
        real8 xi[3], ft[3];
        int j;

        xi[0] = x2-x1;
        xi[1] = y2-y1;
        xi[2] = z2-z1;

        cross(sigb, xi, ft);

        for (j = 0; j < 3; j++)
        {
            f1[j] = ft[j] * 0.5;
            f2[j] = ft[j] * 0.5;
        }
}



/*-------------------------------------------------------------------------
 *
 *      Function:       SelfForce
 *      Description:    Wrapper function which invokes an appropriate force
 *                      function based on whether the code has been compiled
 *                      to do isotropoic or anisotropic elasticity.
 *
 *                      NOTE: For HCP materials, Ecore is adjusted for
 *                            all <c+a> dislocations (see below).
 *
 *      Arguments:
 *
 *------------------------------------------------------------------------*/
void SelfForce(Home_t *home, int coreOnly, real8 MU, real8 NU,
               real8 bx, real8 by, real8 bz,
               real8 x1, real8 y1, real8 z1,
               real8 x2, real8 y2, real8 z2,
               real8 a,  real8 Ecore,
               real8 f1[3], real8 f2[3])
{

/*
 *      Invoke the appropriate self force function
 */
#ifdef ANISOTROPIC
        SelfForceAnisotropic(home, coreOnly,         bx,by,bz, x1,y1,z1, x2,y2,z2, a,Ecore, f1,f2);
#else
        SelfForceIsotropic  (home, coreOnly, MU, NU, bx,by,bz, x1,y1,z1, x2,y2,z2, a,Ecore, f1,f2);
#endif
        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:       LineTensionForce
 *      Description:    Wrapper function which invokes an appropriate force
 *                      function based on whether the code has been compiled
 *                      to do isotropoic or anisotropic elasticity
 *
 *      Arguments:
 *
 *------------------------------------------------------------------------*/
void LineTensionForce(Home_t *home,
                      real8 x1, real8 y1, real8 z1,
                      real8 x2, real8 y2, real8 z2,
                      real8 bx, real8 by, real8 bz,
                      real8 f1[3], real8 f2[3])
{
        Param_t *param = home->param;

        real8 MU = param->shearModulus;
        real8 NU = param->pois;
        real8 a  = param->rc;

/*
 *      For calculating line tension, the isotropic and anisotropic
 *      core energies are different.  For isotropic we reset it as
 *      below.  For anisotropic, we use the existing core energy which
 *      is either provided by the user or should default to log(param->rc/0.1).
 *
 *      WARNING: The function SelfForce redefines Ecore for HCP materials.
 *               The value of Ecore below is not used for HCP materials.
 *
 */
#ifdef ANISOTROPIC
        real8 Ecore = param->Ecore;
#else
        real8 Ecore = param->TensionFactor * MU / 2.0;
#endif

        real8   extStress[3][3];

        extStress[0][0] = param->appliedStress[0];
        extStress[1][1] = param->appliedStress[1];
        extStress[2][2] = param->appliedStress[2];
        extStress[1][2] = param->appliedStress[3];
        extStress[2][0] = param->appliedStress[4];
        extStress[0][1] = param->appliedStress[5];
        extStress[2][1] = extStress[1][2];
        extStress[0][2] = extStress[2][0];
        extStress[1][0] = extStress[0][1];

/*
 *      If the segment is zero length (which can occur when testing whether
 *      multinodes should be split or not), set forces to zero and return.
 */
        real8 dx = x2 - x1;
        real8 dy = y2 - y1;
        real8 dz = z2 - z1;

        if ((dx*dx + dy*dy + dz*dz) < 1.0e-06)
        {
            f1[0] = 0.0;
            f1[1] = 0.0;
            f1[2] = 0.0;

            f2[0] = 0.0;
            f2[1] = 0.0;
            f2[2] = 0.0;

            return;
        }

        ZImage(param, &dx, &dy, &dz);

        x2 = x1 + dx;
        y2 = y1 + dy;
        z2 = z1 + dz;

        real8   fSelf1[3], fSelf2[3];
        real8   fPK1  [3], fPK2  [3];

        SelfForce(home, 1, MU, NU, bx, by, bz, x1, y1, z1, x2, y2, z2, a, Ecore, fSelf1, fSelf2);

        ExtPKForce(extStress, bx, by, bz, x1, y1, z1, x2, y2, z2, fPK1, fPK2);

        f1[0] = fSelf1[0] + fPK1[0];
        f1[1] = fSelf1[1] + fPK1[1];
        f1[2] = fSelf1[2] + fPK1[2];

        f2[0] = fSelf2[0] + fPK2[0];
        f2[1] = fSelf2[1] + fPK2[1];
        f2[2] = fSelf2[2] + fPK2[2];

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:       GetFieldPointStress
 *      Description:
 *
 *      Arguments:
 *
 *------------------------------------------------------------------------*/
void GetFieldPointStress(Home_t *home, real8 x, real8 y, real8 z, real8 totStress[3][3])
{
        Param_t *param = home->param;

        real8    a     = param->rc;
        real8    MU    = param->shearModulus;
        real8    NU    = param->pois;

        totStress[0][0] = 0.0;
        totStress[0][1] = 0.0;
        totStress[0][2] = 0.0;
        totStress[1][0] = 0.0;
        totStress[1][1] = 0.0;
        totStress[1][2] = 0.0;
        totStress[2][0] = 0.0;
        totStress[2][1] = 0.0;
        totStress[2][2] = 0.0;

/*
 *      Get the indices of the cell containing the field point and
 *      determine the block of cells immediately neighboring that
 *      cell.
 */
        real8 Lx = param->Lx;
        real8 Ly = param->Ly;
        real8 Lz = param->Lz;

        int   cellX = (int)((x - param->minSideX) / (Lx / param->nXcells)) + 1;
        int   cellY = (int)((y - param->minSideY) / (Ly / param->nYcells)) + 1;
        int   cellZ = (int)((z - param->minSideZ) / (Lz / param->nZcells)) + 1;

/*
 *      Determine the minimum and maximum cell indices (in each
 *      dimension) of the block of cells encompassing the cell
 *      containing the field point, and all immediate neighbors
 *      of that cell.  Allow for periodic boundary conditions.
 */
        int  minXIndex=0, minYIndex=0, minZIndex=0;
        int  maxXIndex=0, maxYIndex=0, maxZIndex=0;

        if (param->xBoundType == Periodic) { minXIndex = MAX(0, cellX-1); maxXIndex = MIN(param->nXcells+1, cellX+1); }
        else                               { minXIndex = MAX(1, cellX-1); maxXIndex = MIN(param->nXcells  , cellX+1); }

        if (param->yBoundType == Periodic) { minYIndex = MAX(0, cellY-1); maxYIndex = MIN(param->nYcells+1, cellY+1); }
        else                               { minYIndex = MAX(1, cellY-1); maxYIndex = MIN(param->nYcells  , cellY+1); }

        if (param->zBoundType == Periodic) { minZIndex = MAX(0, cellZ-1); maxZIndex = MIN(param->nZcells+1, cellZ+1); }
        else                               { minZIndex = MAX(1, cellZ-1); maxZIndex = MIN(param->nZcells  , cellZ+1); }

/*
 *      Loop though all the cells in the block.
 */
        for (int cx = minXIndex; cx <= maxXIndex; cx++)
        {
            for (int cy = minYIndex; cy <= maxYIndex; cy++)
            {
                for (int cz = minZIndex; cz <= maxZIndex; cz++)
                {
                    int      cellIndex = EncodeCellIdx(home, cx, cy, cz);
                    Cell_t  *cell      = LookupCell(home, cellIndex);

                    if (cell == (Cell_t *)NULL) continue;

/*
 *                  Find the center of this cell and convert it to the
 *                  point in the image closest to the field point.
 *                  We need use the cell center as the base because
 *                  all the Burgers vectors of segments inside one cell
 *                  (belonging to the same PBC image) are summed
 *                  together for 'remote' stress calculations, so
 *                  cannot separate them for 'local' stress calculations.
 */
                    real8 xc = cell->center[0];
                    real8 yc = cell->center[1];
                    real8 zc = cell->center[2];

                    PBCPOSITION(param, x, y, z, &xc, &yc, &zc);

/*
 *                  Loop over all nodes in this cell and over each segment
 *                  attached to the node.  Skip any segment that is not
 *                  owned by node1.
 */
                    Node_t *node1 = cell->nodeQ;

                    for (; node1 != (Node_t *)NULL; node1=node1->nextInCell)
                    {
                        for (int arm = 0; arm < node1->numNbrs; arm++)
                        {
                            Node_t *node2 = GetNeighborNode(home, node1, arm);

                            if (node2 == (Node_t *)NULL) continue;
                            if (OrderNodes(node1, node2) >= 0) continue;

                            real8 p1x = node1->x;
                            real8 p1y = node1->y;
                            real8 p1z = node1->z;

                            real8 p2x = node2->x;
                            real8 p2y = node2->y;
                            real8 p2z = node2->z;

/*
 *                          Find coordinates of the image of the first node
 *                          closest to the center of the cell, and use the
 *                          image of the second node closest to that point.
 */
                            PBCPOSITION(param, xc, yc, zc, &p1x, &p1y, &p1z);
                            PBCPOSITION(param, p1x, p1y, p1z, &p2x, &p2y, &p2z);

                            real8 bx = node1->burgX[arm];
                            real8 by = node1->burgY[arm];
                            real8 bz = node1->burgZ[arm];

/*
 *                          The 6 components of stress returned from
 *                          the StressDueToSeg() call are in the order:
 *                              [sigma00
 *                               sigma11
 *                               sigma22
 *                               sigma12
 *                               sigma02
 *                               sigma01]
 *                          This order is different from the components of
 *                          <appliedStress> as given in the control parameter
 *                          file.
 */
                            real8 stress[6] = { 0.0 };

                            StressDueToSeg(home,
                                           x, y, z,
                                           p1x, p1y, p1z,
                                           p2x, p2y, p2z, bx, by, bz, a, MU, NU, stress);

                            totStress[0][0] += stress[0];
                            totStress[1][1] += stress[1];
                            totStress[2][2] += stress[2];
                            totStress[0][1] += stress[5];
                            totStress[1][2] += stress[3];
                            totStress[0][2] += stress[4];
                        }
                    }
                }
            }
        }

        totStress[1][0] = totStress[0][1];
        totStress[2][0] = totStress[0][2];
        totStress[2][1] = totStress[1][2];

        return;
}
