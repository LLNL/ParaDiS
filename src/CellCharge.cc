/***************************************************************************
 *
 *      Module:       CellCharge
 *      Description:  Sum the contributions from each segment in a cell to
 *                    the total charge tensor in that cell, and distribute
 *                    the sum to all processors. The result is that all
 *                    processors have the net charge tensor of each cell in
 *                    the problem.
 *
 *      Includes functions:
 *
 *          CellCharge()
 *          FMSetTaylorExpansions()
 *          FMCellCharge()
 *          MonopoleCellCharge()
 *
 *      Last Modified: 04/08/2008 - Added explicit initialization of
 *                                  taylor coefficients for highest layer
 *                                  fmm cell if PBC is disabled.
 *
 ***************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "FM.h"

#ifdef ANISOTROPIC
#include "G2sigma_tilde.h"
#endif


/*-------------------------------------------------------------------------
 *
 *      Function:     FMSetTaylorExpansions
 *      Description:  This is basically a control function to handle
 *                    the downward pass through the FM hierarchy
 *                    calculating calculating the taylor expansions
 *                    for each cell, and sending down fully aggregated
 *                    multipole expansions as necessary.
 *
 *      Last Modified: 04/08/2008 - Fixed call to FMFindNearNbrs()
 *
 *-----------------------------------------------------------------------*/
void FMSetTaylorExpansions(Home_t *home)
{
        int       layerID;
        int       cx, cy, cz;
        int       cellID;
        int       trimOverlap;
        int       *tMin, *tMax;
        int       skipX[3], skipY[3], skipZ[3];
        int       tmp1Min[3], tmp1Max[3], tmp2Min[3], tmp2Max[3];
        FMLayer_t *layer;
        FMCell_t  *cell;

        static real8     *tCoeff     = (real8 *)NULL;
        static real8     *tCoeffThrd = (real8 *)NULL;
#ifdef _OPENMP
#pragma omp threadprivate (tCoeff, tCoeffThrd)
#endif

#ifdef ESHELBY
        static real8     *eshelbytCoeff     = (real8 *)NULL;
        static real8     *eshelbytCoeffThrd = (real8 *)NULL;
#ifdef _OPENMP
#pragma omp threadprivate (eshelbytCoeff, eshelbytCoeffThrd)
#endif
#endif  // ESHELBY

        Param_t *param = home->param;

#ifdef ESHELBY
        real8 MU = param->shearModulus;
        real8 NU = param->pois;
#endif

        int   mpOrder=0, tOrder=0; //, maxOrder=0;
        int   numExpansionCoeff=0;
        int   numExpansionCoefflg=0;

        if (param->fmEnabled)
        {
            mpOrder        = param->fmMPOrder;
            tOrder         = param->fmExpansionOrder;
            //maxOrder       = MAX(mpOrder,tOrder);
            numExpansionCoeff = home->fmNumExpansionCoeff;
            numExpansionCoefflg = home->fmNumExpansionCoefflg;
       }

#ifdef ESHELBY
        int eshelbympOrder        = 0;
        int eshelbytOrder         = 0;
        int eshelbynumTaylorCoeff = 0;

        if ((param->eshelbyfmEnabled) && (!home->eshelbyFMInitialized))
        {
            eshelbympOrder        = param->eshelbyfmMPOrder;
            eshelbytOrder         = param->eshelbyfmTaylorOrder;
            eshelbynumTaylorCoeff = home->eshelbyfmNumTaylorCoeff;
        }
#endif  // ESHELBY


/*
 *      Loop through all FM layers except the lowest (coarsest),
 *      do the downward communication pass to shift taylor expansions
 *      from current layer cells to the next layer down as well
 *      distributing the accumulated lower layers' multipole
 *      expansions to the necessary cells at the lower layer.
 */
        for (layerID = 0; layerID < param->fmNumLayers-1; layerID++)
        {
            FMCommDownPass(home, layerID);
/*
 *          Each domain now adjusts the taylor expansions for cells
 *          it owns at the next lower layer based on the multipole
 *          expansions for all near neighbor cells excluding immediate
 *          neighbors.
 */
            layer = &home->fmLayer[layerID+1];

/*
 *          Loop through all cells owned by this domain at the next
 *          layer down (i.e. more refined) in the hierarchy.
 */
            for (int j=0; j < layer->ownedCnt; j++) {
                int totalNearbyCells, cellsPerPlane, cellsPerRow;
                int numNearbyCells[3];
                int cellIndices[3];

                cellID = layer->ownedMin + j;
                IDToFMCellIndices(cellID, layer->lDim, cellIndices);

                cx = cellIndices[X];
                cy = cellIndices[Y];
                cz = cellIndices[Z];

                cell = LookupFMCell(layer->cellTable, cellID);

                skipX[0] = cx-1; skipX[1] = cx; skipX[2] = cx+1;
                skipY[0] = cy-1; skipY[1] = cy; skipY[2] = cy+1;
                skipZ[0] = cz-1; skipZ[1] = cz; skipZ[2] = cz+1;

/*
 *              Now find all the "near" neighbors of the current cell
 *              and loop though each one.
 */
                tmp1Min[X] = cx;
                tmp1Min[Y] = cy;
                tmp1Min[Z] = cz;

                tmp1Max[X] = cx;
                tmp1Max[Y] = cy;
                tmp1Max[Z] = cz;

                trimOverlap = 0;

                FMFindNearNbrs(home, layerID+1, tmp1Min, tmp1Max,
                               tmp2Min, tmp2Max, trimOverlap);

                tMin = tmp2Min;
                tMax = tmp2Max;

                numNearbyCells[X] = tMax[X] - tMin[X] + 1;
                numNearbyCells[Y] = tMax[Y] - tMin[Y] + 1;
                numNearbyCells[Z] = tMax[Z] - tMin[Z] + 1;

                totalNearbyCells = numNearbyCells[X] *
                                   numNearbyCells[Y] *
                                   numNearbyCells[Z];

                cellsPerPlane = numNearbyCells[Y] * numNearbyCells[Z];
                cellsPerRow   = numNearbyCells[Z];

#ifdef _OPENMP
#pragma omp parallel
#endif
                {
                    int i;
                    int cellIter;
                    int threadID, threadIterStart, threadIterEnd;

                    GetThreadIterationIndices(totalNearbyCells,
                                              &threadID,
                                              &threadIterStart,
                                              &threadIterEnd);

/*
 *                  NOTE: The tCoeff* arrays are allocated one time
 *                  per thread and reused from then on.  They are
 *                  NEVER deallocated, so may show up as a mem leak
 */
                    if (param->fmEnabled)
                    {
                        int memSize = numExpansionCoefflg * sizeof(real8);

                        if (tCoeff == (real8 *)NULL) {
                            tCoeff     = (real8 *) malloc(   memSize);
                            tCoeffThrd = (real8 *) calloc(1, memSize);
                        } else {
                            memset(tCoeffThrd, 0, memSize);
                        }
                    }

#ifdef ESHELBY
/*
 *                  NOTE: The eshelbytCoeff* arrays are allocate one
 *                  time per thread and reused from then on.  They are
 *                  NEVER deallocated, so may show up as a mem leak
 *                  They need to be initialized to zero.
 */
                     if (param->eshelbyfmEnabled)
                     {
                          if ( eshelbytCoeff == (real8 *)NULL )
                          {
                             int memSize           = eshelbynumTaylorCoeff * sizeof(real8);
                                 eshelbytCoeff     = (real8 *) malloc(   memSize);
                                 eshelbytCoeffThrd = (real8 *) calloc(1, memSize);
                          }
                          else if(!home->eshelbyFMInitialized)
                          {
                             int memSize           = eshelbynumTaylorCoeff * sizeof(real8);
                             memset(eshelbytCoeffThrd, 0, memSize);
                          }
                     }
#endif  // ESHELBY

                    for (cellIter = threadIterStart; (cellIter < threadIterEnd); cellIter++)
                    {
                        int      tx, ty, tz;
                        int      cellX, cellY, cellZ;
                        int      nCellID;
                        real8    R[3];
                        FMCell_t *nCell;

                        cellX = cellIter / cellsPerPlane;
                        cellY = (cellIter - cellX * cellsPerPlane) /
                                cellsPerRow;
                        cellZ = cellIter -
                                ((cellX * cellsPerPlane) +
                                 (cellY * cellsPerRow));

                        cellX += tMin[X];
                        cellY += tMin[Y];
                        cellZ += tMin[Z];

/*
 *                      The initial indices may in fact represent periodic
 *                      cells outside the bounds of the problem, so
 *                      we fiddle with the cell{XYZ} indices to come up
 *                      with the t{xyz} indices matching the image of the
 *                      cell that is within the primary problem space
 */
                        tx = cellX % layer->lDim[X];
                        ty = cellY % layer->lDim[Y];
                        tz = cellZ % layer->lDim[Z];

                        tx = GETPBCINDEX(tx, layer->lDim[X]);
                        ty = GETPBCINDEX(ty, layer->lDim[Y]);
                        tz = GETPBCINDEX(tz, layer->lDim[Z]);

/*
 *                      If this neighbor is an immediate
 *                      neighbor of the owned cell (or the
 *                      cell itself), skip it.
 */
                        if ((cellX == skipX[0] || cellX == skipX[1] ||
                             cellX == skipX[2]) &&
                            (cellY == skipY[0] || cellY == skipY[1] ||
                             cellY == skipY[2]) &&
                            (cellZ == skipZ[0] || cellZ == skipZ[1] ||
                             cellZ == skipZ[2])) {
                            continue;
                        }

/*
 *                      Find the vector from the center of the
 *                      current cell to the center of the multipole
 *                      expansion in the neighboring cell.
 */
                        nCellID = EncodeFMCellIndex(layer->lDim, tx, ty, tz);
                        nCell   = LookupFMCell(layer->cellTable, nCellID);

                        R[X] = (layer->cellSize[X] * (cx - cellX));
                        R[Y] = (layer->cellSize[Y] * (cy - cellY));
                        R[Z] = (layer->cellSize[Z] * (cz - cellZ));

                        if (param->fmEnabled) {
/*
 *                          Calculate the contribution to this cell's
 *                          expansion coefficients from the neighbor's
 *                          multipole expansion
 */
                           memset(tCoeff, 0, numExpansionCoefflg * sizeof(real8));

                           MomentToLocal(home, layer->cellSize, R,
                                          nCell->mpCoeff, tCoeff);

                            for (i=0; i<numExpansionCoefflg; i++) {
                                tCoeffThrd[i] += tCoeff[i];
                            }


#ifdef CALCENERGY
                            real8 W = 0.0;
                            W = FMMEnergy(home, mpOrder, R, cell->mpCoeff,nCell->mpCoeff);
                            param->FMMEnergy += W;
#endif



                        }

#ifdef ESHELBY
                        if ( (param->eshelbyfmEnabled) && !home->eshelbyFMInitialized)
                        {
                           EshelbyMakeTaylor(R, nCell->eshelbympCoeff,
                                             eshelbympOrder, eshelbytOrder, MU, NU, eshelbytCoeff);

                            for (i=0; i<eshelbynumTaylorCoeff; i++)
                            {
                                eshelbytCoeffThrd[i] += eshelbytCoeff[i];
                            }
                        }
#endif  // ESHELBY

                    }  /* end for (cellIter = threadIterStart; ...) */

#ifdef BBFMM
/*
 *                  We have the contribution to the expansion coefficients
 *                  from all the nearby cells, but we must multiply them
 *                  by U before we go on.
 */
                    bbfmm_PostM2L(layer->cellSize[0],tCoeffThrd);
#endif

#ifdef UNIFORMFMM
                    // tCoeff's (2n-1)^3 in, n^3 out.
                    FFTBigToSmall(layer->cellSize[0], tCoeffThrd);
#endif

#ifdef ANISOTROPIC
#ifdef TAYLORFMM
/*
 *                  ATTENTION : This method is not fully debugged
 *                  If anisotropic Taylor FMM is being used, then post multiply the
 *                  alpha_tilde (tCoeffThrd array) returned from MkTaylor()
 *                  by C to get the standard alpha and copy it back into
 *                  tCoeffThrd.
 */

                    if (param->fmEnabled && param->fmEnableAnisotropy) {

                        alpha_tilde2alpha(param->fmExpansionOrder,
                                          (real8 (*)[3][3][3]) home->anisoVars.elasticConstantMatrix4D,
                                          (real8 (*)[3][3]   ) tCoeffThrd,
                                          (real8 (*)[3][3]   ) tCoeff );

                        for (i = 0; i < numExpansionCoeff; i++) {
                            tCoeffThrd[i] = tCoeff[i];
                        }
                    }
#endif
#endif

/*
 *                  Add the thread's contribution to the cell's
 *                  expansion coefficients.  Need the "critical" section
 *                  to insure only 1 thread updates the coefficients
 *                  at a time.
 */
                    if (param->fmEnabled)
                    {
#ifdef _OPENMP
#pragma omp critical(CRIT_SET_TAYLOR)
#endif
                        for (i = 0; i < numExpansionCoeff; i++) {
                            cell->expansionCoeff[i] += tCoeffThrd[i];
                        }
                    }

#ifdef ESHELBY
                    if ((param->eshelbyfmEnabled) && (!home->eshelbyFMInitialized))
                    {
#ifdef _OPENMP
#pragma omp critical(CRIT_SET_ESHELBY_TAYLOR)
#endif
                        for (i = 0; i < eshelbynumTaylorCoeff; i++) {
                            cell->eshelbytaylorCoeff[i] += eshelbytCoeffThrd[i];
                        }
                    }
#endif  // ESHELBY

                }  /* end omp parallel section */

            }  /* loop over cells */

        }  /* loop over layers */

#if TAYLORFMM
/*
 *      For simulations using periodic boundaries, the taylor coefficients
 *      of all cells at the lowest layer now need to be adjusted in order
 *      to account for conditional convergence issues.
 *      FOR THE OTHER FMM METHODS (BBFMM and UNIFORM GRID FMM), THE MEAN STRESS IS DONE
 *      IN DoTableCorrection.....
 */
        if ((param->xBoundType == Periodic) &&
            (param->yBoundType == Periodic) &&
            (param->zBoundType == Periodic)) {
           MeanStressCorrection(home);
        }
#endif

/*
 *      The downward pass is complete, but each domain still has to
 *      distribute the expansion expansion coeeficients for any cell
 *      it owns (at the most refined FM layer) to all domains
 *      intersecting that cell.
 */
        FMDistTaylorExp(home);

        return;
}


static void FMCellCharge(Home_t *home)
{
        int       i, inode, inbr, layerID;
        int       cx, cy, cz;
        int       cellID;
        int       *bMin, *bMax;
        real8     p1[3], p2[3];
        real8     vec1[3], vec2[3], burg[3];
        Node_t    *node, *nbr;
        FMLayer_t *layer;
        FMCell_t  *cell;
        Cell_t    *simCell;

        Param_t *param = home->param;

        int numMPCoeff        = 0;
        int numMPCoeffsm      = 0;
        int numExpansionCoeff = 0;
        real8 *etatemp        = 0;

        if (param->fmEnabled)
        {
            numMPCoeff        = home->fmNumMPCoeff;
            numMPCoeffsm      = home->fmNumMPCoeffsm;
            numExpansionCoeff = home->fmNumExpansionCoeff;
            etatemp = (real8 *)malloc(numMPCoeffsm * sizeof(real8));
        }

#ifdef ESHELBY
        int eshelbynumMPCoeff     = 0;
        int eshelbynumTaylorCoeff = 0;

        if ((param->eshelbyfmEnabled) && (!home->eshelbyFMInitialized))
        {
            eshelbynumMPCoeff     = home->eshelbyfmNumMPCoeff;
            eshelbynumTaylorCoeff = home->eshelbyfmNumTaylorCoeff;
        }
#endif  // ESHELBY

/*
 *      Zero out multipole expansion data for all cells intersecting this
 *      domain at the lowest FM layer.
 */
        layer = &home->fmLayer[param->fmNumLayers-1];

        bMin = layer->intersectMin;
        bMax = layer->intersectMax;

        for (cx = bMin[X]; cx <= bMax[X]; cx++) {
            for (cy = bMin[Y]; cy <= bMax[Y]; cy++) {
                for (cz = bMin[Z]; cz <= bMax[Z]; cz++)
                {
                    cellID = EncodeFMCellIndex(layer->lDim, cx, cy, cz);
                    cell   = LookupFMCell(layer->cellTable, cellID);

                    if (param->fmEnabled)
                    {
                        memset(cell->mpCoeff, 0, numMPCoeff * sizeof(real8));
                    }

#ifdef ESHELBY
                    if ((param->eshelbyfmEnabled) && (!home->eshelbyFMInitialized))
                    {
                        memset(cell->eshelbympCoeff, 0, eshelbynumMPCoeff * sizeof(real8));
                    }
#endif  // ESHELBY
                }
            }
        }

/*
 *      loop through the native nodes.  For each node look at
 *      its neighbors.  If neighbor has a higher tag than node,
 *      then it is a non-redundant segment.  For consistency with
 *      the way segments are chosen for local force calculation,
 *      a segment is considered to belong to the cell of its
 *      higher-tagged node.
 */
        inode = (param->fmEnabled ? 0 : home->newNodeKeyPtr);

        for ( ; inode < home->newNodeKeyPtr; inode++) {

            if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) continue;

            p1[X] = node->x;
            p1[Y] = node->y;
            p1[Z] = node->z;

/*
 *          Find out which simulation cell the current node is in and convert
 *          that cell's indices into the FM cell indices and ID for the
 *          corresponding FM cell at the bottom (most refined) FM layer.
 */
            simCell = LookupCell(home, node->cellIdx);
            cx = simCell->xIndex;
            cy = simCell->yIndex;
            cz = simCell->zIndex;

            cx--;
            cy--;
            cz--;

            layer  = &home->fmLayer[param->fmNumLayers-1];
            cellID = EncodeFMCellIndex(layer->lDim, cx, cy, cz);
            cell   = LookupFMCell(layer->cellTable, cellID);

            if (cell == (FMCell_t *)NULL) {
                Fatal("NULL FM cell ptr");
            }

            for (inbr = 0; inbr < node->numNbrs; inbr++) {

               nbr = GetNeighborNode (home, node, inbr);

               if (nbr == (Node_t *)NULL) {
                   printf("WARNING: Neighbor not found at %s line %d\n",
                          __FILE__, __LINE__);
                   continue;
               }

               if (NodeOwnsSeg(home, node, nbr) == 0) {
                   continue;
               }

               p2[X] = nbr->x;
               p2[Y] = nbr->y;
               p2[Z] = nbr->z;

               PBCPOSITION(param, p1[X], p1[Y], p1[Z], &p2[X], &p2[Y], &p2[Z]);

               burg[X] = node->burgX[inbr];
               burg[Y] = node->burgY[inbr];
               burg[Z] = node->burgZ[inbr];

/*
 *             Set vector from the segment starting point to segment
 *             end point
 */
               vec2[X] = p2[X] - p1[X];
               vec2[Y] = p2[Y] - p1[Y];
               vec2[Z] = p2[Z] - p1[Z];

/*
 *             Set vector from the segment starting point to the
 *             FM cell expansion center and increment the total
 *             multipole expansion for the lowest layer FM cell
 *             containing the segment with the contribution from this segment.
 */
               vec1[X] = p1[X] - cell->cellCtr[X];
               vec1[Y] = p1[Y] - cell->cellCtr[Y];
               vec1[Z] = p1[Z] - cell->cellCtr[Z];

               SegmentToMoment(home, layer->cellSize, cell->cellCtr,
                               vec1, vec2, burg, etatemp);


               for (i = 0; i < home->fmNumMPCoeffsm; i++) {
                   cell->mpCoeff[i] += etatemp[i];
               }

           }  /* end for (inbr = 0; ...)  */
        }  /* end for (inode = 0; ...)  */

        if (param->fmEnabled) {
            free(etatemp);
            etatemp = (real8 *)NULL;
        }

#ifdef ESHELBY
        if ((param->eshelbyfmEnabled) && (!home->eshelbyFMInitialized))
        {
/*
 *          Loop through all native inclusions and update the multipole
 *          expansion for the owning cell appropriately.
 */
            etatemp = (real8 *)malloc(eshelbynumMPCoeff * sizeof(real8));

            for (i = 0; i < home->locInclusionCount; i++)
            {
                EInclusion_t *inclusion = &home->eshelbyInclusions[i];

/*
 *              The inclusion's cellID is adjusted for ghost cells
 *              so must be converted to the FM cell ID
 */
                simCell = LookupCell(home, inclusion->cellID);
                cx = simCell->xIndex;
                cy = simCell->yIndex;
                cz = simCell->zIndex;

                cx--;
                cy--;
                cz--;

                layer  = &home->fmLayer[param->fmNumLayers-1];
                cellID = EncodeFMCellIndex(layer->lDim, cx, cy, cz);
                cell   = LookupFMCell(layer->cellTable, cellID);

                real8 Ev = inclusion->strainField[0] +
                           inclusion->strainField[1] +
                           inclusion->strainField[2];

                if (inclusion->radius[0] == inclusion->radius[1] &&
                    inclusion->radius[0] == inclusion->radius[2])
                {
                   EshelbyMakeEta(Ev, inclusion->radius[0], inclusion->position,
                                  cell->cellCtr, param->eshelbyfmMPOrder,
                                  etatemp);
                }
                else
                   Fatal("Calling Eshelby FMM with an ellipsoidal inclusion is not supported.\n");

                for (int j=0; j < eshelbynumMPCoeff; j++) {
                    cell->eshelbympCoeff[j] += etatemp[j];
                }
            }

#if 0
/*
 *          Need to explicitly zero out the taylor expansion at the
 *          highest layer before continuing.
 *
 *          NOTE:  At the time of this writing, the eshelby inclusions
 *                 are not permitted to move.  Hence, the multipole
 *                 and taylor expansion coefficients need only be
 *                 calculated once.  Since the memory for the taylor
 *                 coefficients was calloc'ed, we don't need to
 *                 re-zero it unless future changes allow inclusions
 *                 to move necessitating recalculation of the the
 *                 coefficients.
 */
            layer = &home->fmLayer[0];
            cell = LookupFMCell(layer->cellTable, 0);
            if (cell != (FMCell_t *)NULL) {
                memset(cell->eshelbytaylorCoeff, 0, eshelbynumTaylorCoeff *
                       sizeof(real8));
            }
#endif
        }
#endif  // ESHELBY


/*
 *      We have the multipole contribution from each native segment
 *      to the cells at the lowest FM layer, now we have to pass the
 *      data up the FM hierarchy.
 *      UPWARD PASS
 */
        for (layerID = param->fmNumLayers-1; layerID >= 0; layerID--) {
            FMCommUpPass(home, layerID);
        }

        if (param->fmEnabled) {
/*
 *          For the standard (non-eshelby) Fast Multipole stuff only,
 *          calculate the expansion expansion at the highest (coarsest) layer
 *          of the FM hierarchy.  This expansion expansion will account for
 *          all necessary periodic images of the full problem space except
 *          for the near neighbor images of the primary image.
 *          Add PBC
 */
            if ((param->xBoundType == Periodic) ||
                (param->yBoundType == Periodic) ||
                (param->zBoundType == Periodic)) {
               DoTableCorrection(home);
            } else {
/*
 *              If periodic boundaries are not enabled, we have to explicitly
 *              initialize the expansion coefficients at the highest FM layer.
 *              With periodic boundaries, this is taken care of in
 *              DoTableCorrection().
 */
                layer = &home->fmLayer[0];
                cell  = LookupFMCell(layer->cellTable, 0);
                if (cell != (FMCell_t *)NULL) {
                    memset(cell->expansionCoeff, 0,
                           home->fmNumExpansionCoeff * sizeof(real8));
                }
            }
        }

/*
 *      Now call the function that will handle passing all necessary
 *      data back down the FM hierarchy shifting the expansion expansions
 *      from ancestor to descendant cells as needed.
 *      DOWNWARD PASS
 */
        FMSetTaylorExpansions(home);

/*
 *      And since most of the eshelby FMM stuff need only be calculated
 *      one time, set a flag so we know we've done it.
 */
#ifdef ESHELBY
        home->eshelbyFMInitialized = 1;
#endif

        return;
}


static void MonopoleCellCharge(Home_t *home)
{
        int     i;
        int     nodeIndex, numNodes;
        int     numCells, numChargeVals;
        real8   *cellCharge;
        Param_t *param;


        param = home->param;

/*
 *      Allocate an array of 'charge' tensors and zero all the
 *      components.
 */
        numCells = param->nXcells * param->nYcells * param->nZcells;
        numChargeVals = numCells * 9;

        cellCharge = (real8 *)malloc(numChargeVals * sizeof(real8));

#ifdef _OPENMP
#pragma omp for
#endif
        for (i = 0 ; i < numChargeVals; i++) {
            cellCharge[i] = 0.0 ;
        }

/*
 *      Loop over each segment owned by each local node and add the
 *      segment's contribution to the charge tensor for the cell owning
 *      the segment.
 *
 *      For consistency with the way segments are chosen for local force
 *      calculation, a segment is considered to belong to the cell
 *      containing the node owning the segment
 */
        numNodes = home->newNodeKeyPtr;

#ifdef _OPENMP
#pragma omp for
#endif
        for (nodeIndex = 0; nodeIndex < numNodes; nodeIndex++) {
            int    m, n, offset;
            int    cellID, cellXIndex, cellYIndex, cellZIndex;
            int    numNbrs, nbrIndex;
            real8  x1, y1, z1;
            real8  cellChargeTensor[3][3];
            Cell_t *cell;
            Node_t *node;

            if ((node = home->nodeKeys[nodeIndex]) == (Node_t *)NULL) {
                continue;
            }

            Matrix33_Zero(cellChargeTensor);

            x1 = node->x;
            y1 = node->y;
            z1 = node->z;

/*
 *          Decrement the cell indices so we are dealing with the true cell
 *          indices here, not indices that have been shifted to account for
 *          ghost cells.
 */
            cell = LookupCell(home, node->cellIdx);

            cellXIndex = cell->xIndex - 1;
            cellYIndex = cell->yIndex - 1;
            cellZIndex = cell->zIndex - 1;

            cellID = cellZIndex +
                     param->nZcells * cellYIndex +
                     param->nZcells * param->nYcells * cellXIndex;

            numNbrs = node->numNbrs;

            for (nbrIndex = 0; nbrIndex < numNbrs; nbrIndex++) {
                real8  b[3], dl[3];
                Node_t *nbrNode;

                nbrNode = GetNeighborNode(home, node, nbrIndex);

                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }

                if (NodeOwnsSeg(home, node, nbrNode) == 0) {
                    continue;
                }

                dl[0] = nbrNode->x - x1;
                dl[1] = nbrNode->y - y1;
                dl[2] = nbrNode->z - z1;

                ZImage (param, &dl[0], &dl[1], &dl[2]);

                b[0] = node->burgX[nbrIndex];
                b[1] = node->burgY[nbrIndex];
                b[2] = node->burgZ[nbrIndex];

/*
 *              Increment a thread-local sum of the contributions to the
 *              cell's charge tensor from the node's segments.
 */
                for (m = 0; m < 3; m++) {
                    for (n = 0; n < 3; n++) {
                        cellChargeTensor[m][n] += b[m] * dl[n];
                    }
                }

            }  /* end for(nbrIndex = 0; ...) */

/*
 *          Calculate the offset to the charge tensor for the node's cell
 *          and increment the tensor with the summed contribution from
 *          this node's segments.
 *
 *          Use the cell's 'lock' to insure that only one thread is
 *          updating a given charge tensor at a time...
 */
            offset = cellID * 9;

            LOCK(&cell->cellLock);

            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    cellCharge[offset++] += cellChargeTensor[m][n];
                }
            }

            UNLOCK(&cell->cellLock);

        }  /* end for(nodeIndex = 0; ...) */

/*
 *      Sum the individual cell charges over all processors (only a few at
 *      most will contribute to any particular cell) and leave the sum for
 *      each cell on all processors.
 */
#ifdef PARALLEL
        MPI_Allreduce (cellCharge, home->cellCharge, numChargeVals, MPI_DOUBLE,
                       MPI_SUM, MPI_COMM_WORLD);
        free(cellCharge);
#else
        if (home->cellCharge) {
            free(home->cellCharge);
        }
        home->cellCharge = cellCharge;
#endif

        return;
}


void CellCharge(Home_t *home)
{
        TimerStart(home, CELL_CHARGE);

/*
 *      If we are not doing full n^2 force calculations, we need
 *      to prepare for the remote force calcs... but if we are
 *      doing full n^2 force calcs, this stuff is irrelevant.
 */
#ifndef FULL_N2_FORCES
#ifdef ESHELBY
        if ((home->param->fmEnabled) || (home->param->eshelbyfmEnabled))
#else
        if  (home->param->fmEnabled)
#endif
            FMCellCharge(home);
        else {
#ifdef SPECTRAL
            if (home->param->FFTenabled)
                FFTCellCharge(home);
            else
#endif
                MonopoleCellCharge(home);
        }
#endif

        TimerStop(home, CELL_CHARGE);

#ifdef PARALLEL
#ifdef SYNC_TIMERS
        TimerStart (home, CELL_CHARGE_BARRIER);
        MPI_Barrier(MPI_COMM_WORLD);
        TimerStop  (home, CELL_CHARGE_BARRIER);
#endif
#endif
        return;
}
