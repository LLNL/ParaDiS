/**************************************************************************
 *
 *      Module:  This module contains the functions needed for
 *               calculating interactions between dislocation
 *               segments.  See Tom Arsenlis for details on the
 *               method used to do the calculations.
 *
 *      Includes public functions:
 *              LocalSegForces()
 *              CellPriority()
 *              NodeOwnsSeg()
 *              FindSubFSeg()
 *              FindFSegComb()
 *
 *      Includes private functions:
 *              IncrDomSegCommCnts()
 *
 *************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Comm.h"
#include "ParadisThread.h"
#include "ComputeForces.h"
#include "V3.h"
#include "SegmentPairList.h"
#include "NativeSeg.h"
#include "Moment.h"
#include "SSF_Driver.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

#ifdef _ARLFEM
/*---------------------------------------------------------------------------
 *
 *      Function:    VirtualSegForce()
 *
 *      Description: Calculate the force from all the virtual segments
 *                   (i.e. semi-infinite segments extending from any
 *                   surface node outward from the simulation boundary)
 *                   on the segment beginning at <p1> and ending at <p2>.
 *
 *                   Note: Currently, we're only interested in the forces
 *                         on the endpoints of the finite segments.  If
 *                         necessary we can include the force at the surface
 *                         endpoint of the virtual segment.
 *
 *      Arguments:
 *          p1, p2  IN:  Positions of the two endpoints of the segment for
 *                       which forces are being calculated
 *          burg    IN:  burgers vector for the segment <p1> to <p2>.
 *          f1, f2 OUT:  Force from all virtual segments at points <p1>
 *                       and <p2> respectively
 *
 *-------------------------------------------------------------------------*/
void VirtualSegForce(Home_t *home, real8 p1[3], real8 p2[3], real8 burg[3],
                     real8 f1[3], real8 f2[3])
{
        real8 a, MU, NU;

        if (home->surfaceSegList == (real8 *)NULL) {
            return;
        }

        a  = home->param->rc;
        MU = home->param->shearModulus;
        NU = home->param->pois;

        V3_ZERO(f1);
        V3_ZERO(f2);

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int   i, bufIndex;
            int   numSegs;
            int   threadID, threadIterStart, threadIterEnd;
            real8 *segList;
            real8 f1ThreadTot[3], f2ThreadTot[3];

            bufIndex = 0;
            numSegs = home->surfaceSegCount;
            segList = home->surfaceSegList;

            V3_ZERO(f1ThreadTot);
            V3_ZERO(f2ThreadTot);

            GetThreadIterationIndices(numSegs, &threadID, &threadIterStart, &threadIterEnd);

/*
 *          Each thread loops through a portion of the surface segment list.
 *          The list contains all surface intersecting segments in the
 *          simulation space.  The first node of each node-pair for each
 *          segment on the list is the surface node.
 */
            for (i = threadIterStart; i < threadIterEnd; i++) {
                real8 burgSign, virtualSegLen;
                real8 p1S[3], p2S[3], pVirtual[3];
                real8 bVirtual[3], vec1[3];
                real8 f[3], f1[3], f2[3];
/*
 *              Get the positions of the endpoints of the surface
 *              intersecting segment.  <p1S> is the position of the
 *              internal node, <p2S> is the position of the surface
 *              node.
 */
                p1S[0] = segList[bufIndex++];
                p1S[1] = segList[bufIndex++];
                p1S[2] = segList[bufIndex++];

                p2S[0] = segList[bufIndex++];
                p2S[1] = segList[bufIndex++];
                p2S[2] = segList[bufIndex++];

/*
 *              Pull the burgers vector from the buffer
 */
                bVirtual[0] = segList[bufIndex++];
                bVirtual[1] = segList[bufIndex++];
                bVirtual[2] = segList[bufIndex++];

/*
 *              For the surface segment, calculate the vector from the
 *              internal node to the surface node and calculate the
 *              length of the virtual segment we want to use.
 */
                vec1[0] = p2S[0] - p1S[0];
                vec1[1] = p2S[1] - p1S[1];
                vec1[2] = p2S[2] - p1S[2];

                virtualSegLen = 10000 / sqrt(DotProduct(vec1, vec1));

/*
 *              Calculate the position of the "virtual" node that is
 *              outside the surface along the line direction from
 *              p2S to p1S.
 */
                pVirtual[0] = p2S[0] + vec1[0] * virtualSegLen;
                pVirtual[1] = p2S[1] + vec1[1] * virtualSegLen;
                pVirtual[2] = p2S[2] + vec1[2] * virtualSegLen;

                SemiInfiniteSegSegForce(home,
                                        p2S[0], p2S[1], p2S[2],
                                        pVirtual[0], pVirtual[1], pVirtual[2],
                                        p1[X], p1[Y], p1[Z],
                                        p2[X], p2[Y], p2[Z],
                                        bVirtual[0], bVirtual[1], bVirtual[2],
                                        burg[X],burg[Y],burg[Z],
                                        a, MU, NU,
                                        &f[0],  &f[1],  &f[2],
                                        &f1[0], &f1[1], &f1[2],
                                        &f2[0], &f2[1], &f2[2]);
/*
 *              Each thread keeps a running total of the forces it
 *              has calculated for the segment.
 */
                V3_ACCUM(f1ThreadTot, f1);
                V3_ACCUM(f2ThreadTot, f2);
            }

/*
 *          Thread has calculated its portion of the forces on the segment,
 *          so lock access to the total force and update it with this thread's
 *          contribution.
 */
#ifdef _OPENMP
#pragma omp critical (VIRTUAL_SEG_FORCE)
#endif
            {
                V3_ACCUM(f1, f1ThreadTot);
                V3_ACCUM(f2, f2ThreadTot);
            }
        }

        return;
}
#endif  /* end ifdef _ARLFEM */


/*---------------------------------------------------------------------------
 *
 *      Function:     NodeOwnsSeg
 *      Description:  This function examines the two provided nodes and
 *                    determines if the first node has priority over
 *                    the other node.  Priority is determined by
 *                    the following rules:
 *
 *                    WARNING! This function assumes that the nodes are
 *                    in the same or immediately neighboring cells (allowing
 *                    for PBC)!
 *
 *                    - If both nodes are in the same cell the node in the
 *                      higher domain has priority if the nodes are not
 *                      in the same domain, or the node with the higher
 *                      index if the nodes are in the same domain.
 *                    - If the nodes are in different cells, the node with the
 *                      higher cell index owns the segment unless the segment
 *                      crosses a periodic boundary, in which case the node
 *                      with the lower cell index owns the segment.
 *
 *      Arguments
 *          node1  pointer to the first node
 *          node2  pointer to the second node
 *
 *      Returns:   1 if <node1> has priority and owns the segment
 *                 0 if <node2> has priority and owns the segment
 *
 *-------------------------------------------------------------------------*/
int NodeOwnsSeg(Home_t *home, Node_t *node1, Node_t *node2)
{
        if (node1->cellIdx == node2->cellIdx)
        {
            if (node1->myTag.domainID == node2->myTag.domainID)
                return(node1->myTag.index    > node2->myTag.index   );
            else
                return(node1->myTag.domainID > node2->myTag.domainID);
        }
/*
 *      Since the max segment length is less than the cell width,
 *      no segment can span an entire cell.  So, if the cells
 *      containing the two nodes are not physically neighboring, the
 *      segment has crossed a periodic boundary...
 *      Note: The nodal cell IDs convert to indices for *real*
 *      cells, not the indices for periodic cells
 *
 *      For secondary ghost nodes, we may not have allocated a cell
 *      structure, so we just calculate the cell indices directly
 *      rather than do a cell lookup.
 */
        int cell1X=0, cell1Y=0, cell1Z=0;
        int cell2X=0, cell2Y=0, cell2Z=0;

        DecodeCellIdx(home, node1->cellIdx, &cell1X, &cell1Y, &cell1Z);
        DecodeCellIdx(home, node2->cellIdx, &cell2X, &cell2Y, &cell2Z);

        int diffX = abs(cell1X - cell2X);
        int diffY = abs(cell1Y - cell2Y);
        int diffZ = abs(cell1Z - cell2Z);

        if ((diffX > 1) || (diffY > 1) || (diffZ > 1))
            return(node1->cellIdx < node2->cellIdx);

        return(node1->cellIdx > node2->cellIdx);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CellPriority
 *      Description:  Given the IDs of two *neighboring* cells
 *                    (includes neighbors due to periodic boundaries)
 *                    determines the priority of cellID1 related to
 *                    cellID2
 *
 *      Returns:  -1 if cellID1 has lower priority than cellID2
 *                 0 if cellID1 has equal priority to cellID2
 *                 1 if cellID1 has higher priority than cellID2
 *
 *-------------------------------------------------------------------------*/
static int CellPriority(Home_t *home, int cellID1, int cellID2)
{
        int    cell1X, cell1Y, cell1Z;
        int    cell2X, cell2Y, cell2Z;
        int    diffX, diffY, diffZ;
        Cell_t *cell;

        if (cellID1 == cellID2) {
            return(0);
        }

/*
 *      This function assumes the cell structures have been allocated and
 *      initialized for any cell ID passed in.  For secondary ghosts, we may
 *      not have an allocated those structures yet, but since we *shouldn't*
 *      be calling this with cell IDs for secondary ghosts, we should be okay.
 */
        cell = LookupCell(home, cellID1);

        cell1X = cell->xIndex;
        cell1Y = cell->yIndex;
        cell1Z = cell->zIndex;

        cell = LookupCell(home, cellID2);

        cell2X = cell->xIndex;
        cell2Y = cell->yIndex;
        cell2Z = cell->zIndex;

/*
 *      Normally, the cell with the higher ID has a higher priority,
 *      however, if cell2 is more than 1 cell away from cell1 in any
 *      dimension, maximum cell in the current domain in any dimension
 *      the cells are neighbors due to periodic boundaries.  In that
 *      case, the cell with the lower ID has priority.
 */
        diffX = abs(cell1X - cell2X);
        diffY = abs(cell1Y - cell2Y);
        diffZ = abs(cell1Z - cell2Z);

        if ((diffX > 1) || (diffY > 1) || (diffZ > 1)) {
            return(cellID1 < cellID2 ? 1 : -1);
        }

        return(cellID1 > cellID2 ? 1 : -1);
}


/*---------------------------------------------------------------------------
 *
 *      Function:       FindSubFseg
 *      Description:    Given a segment p1-->p2 and the force at each endpoint
 *                      of the segment, estimate the resulting forces on the
 *                      segment pair created by bisecting p1-->p2 at newpos.
 *
 *      Arguments
 *          p1       Coordinates of point 1
 *          p2       Coordinates of point 2 (corresponding to the periodic
 *                   image of p2 closest to point p1)
 *          burg     burgers vector from p1 to p2
 *          oldfp1   force of segment p1-->p2 at point p1
 *          oldfp2   force of segment p1-->p2 at point p2
 *          newpos   coordinates of position along the p1-->p2
 *                   segment at which a new node is to be added.
 *          f0seg1   Resulting forces at p1 on the p1-->newpos segment
 *          f1seg1   Resulting forces at newpos on the p1-->newpos segment
 *          f0seg2   Resulting forces at newpos on the newpos-->p2 segment
 *          f1seg2   Resulting forces at p2 on the newpos-->p2 segment
 *
 *-------------------------------------------------------------------------*/
void FindSubFSeg(Home_t *home, real8 p1[3], real8 p2[3], real8 burg[3],
                 real8 oldfp1[3], real8 oldfp2[3], real8 newpos[3],
                 real8 f0seg1[3], real8 f1seg1[3], real8 f0seg2[3],
                 real8 f1seg2[3])
{
    real8    f1[3], f2[3], f3[3], f4[3];
    real8    fLinv1[3], fLinv2[3], fLinvPos[3];

    Param_t *param = home->param;

    real8 Ecore = param->Ecore;
    real8 MU    = param->shearModulus;
    real8 NU    = param->pois;
    real8 a     = param->rc;

    real8 L     = (p2[X]-p1[X])*(p2[X]-p1[X]) +
                  (p2[Y]-p1[Y])*(p2[Y]-p1[Y]) +
                  (p2[Z]-p1[Z])*(p2[Z]-p1[Z]);

    real8 Lsub1 = (newpos[X]-p1[X])*(newpos[X]-p1[X]) +
                  (newpos[Y]-p1[Y])*(newpos[Y]-p1[Y]) +
                  (newpos[Z]-p1[Z])*(newpos[Z]-p1[Z]);

    real8 Lsub2 = (p2[X]-newpos[X])*(p2[X]-newpos[X]) +
                  (p2[Y]-newpos[Y])*(p2[Y]-newpos[Y]) +
                  (p2[Z]-newpos[Z])*(p2[Z]-newpos[Z]);

/*
 *  It is possible to have a zero-length segment generated during
 *  a collision.  If we find such a case, the force on a segment
 *  between p1 and p2 would equal the force on the non-zero length
 *  segment between newpos and the more distant of p1 and p2. Otherwise,
 *  calc the forces as normal.
 */
    if (Lsub1 < 1.0e-20)
    {
        V3_ZERO(f0seg1);
        V3_ZERO(f1seg1);
        V3_COPY(f0seg2,oldfp1);
        V3_COPY(f1seg2,oldfp2);
        return;
    } else if (Lsub2 < 1.0e-20)
    {
        V3_COPY(f0seg1,oldfp1);
        V3_COPY(f1seg1,oldfp2);
        V3_ZERO(f0seg2);
        V3_ZERO(f1seg2);
        return;
    }


/*
 *  Find the self-force for this segment and subtract it off
 */
    SelfForce(home, 0, MU, NU, burg[X], burg[Y], burg[Z], p1[X], p1[Y], p1[Z],
              p2[X], p2[Y], p2[Z], a, Ecore, f1, f2);

    for (int i=0; (i<3); i++)
    {
        oldfp1[i] -= f1[i];
        oldfp2[i] -= f2[i];
    }

/*
 *  Use the shape function to find the fLinv at all three positions
 */
    for (int i=0; (i<3); i++)
    {
        fLinv1  [i] = ((4.0 * oldfp1[i]) - (2.0 * oldfp2[i])) / L;
        fLinv2  [i] = ((4.0 * oldfp2[i]) - (2.0 * oldfp1[i])) / L;
        fLinvPos[i] = fLinv1[i]*((L-Lsub1)/L) +
                      fLinv2[i]*((L-Lsub2)/L);
    }

/*
 *  Calculate the segment sub forces without self force
 */
    for (int i=0; (i<3); i++)
    {
        f0seg1[i] = ((2.0 * fLinv1  [i]) + fLinvPos[i]) * Lsub1 / 6.0;
        f1seg1[i] = ((2.0 * fLinvPos[i]) + fLinv1  [i]) * Lsub1 / 6.0;
        f0seg2[i] = ((2.0 * fLinvPos[i]) + fLinv2  [i]) * Lsub2 / 6.0;
        f1seg2[i] = ((2.0 * fLinv2  [i]) + fLinvPos[i]) * Lsub2 / 6.0;
    }

/*
 *  Add the self force back into the segment forces along with the
 *  new interaction between the two new segments that was not
 *  included in the previous segment calculation
 */
    SelfForce(home, 0, MU, NU, burg[X], burg[Y], burg[Z], p1[X], p1[Y], p1[Z],
              newpos[X], newpos[Y], newpos[Z], a, Ecore, f1, f2);

    for (int i=0; (i<3); i++)
    {
        f0seg1[i] += f1[i];
        f1seg1[i] += f2[i];
    }

    SelfForce(home, 0, MU, NU, burg[X], burg[Y], burg[Z], newpos[X], newpos[Y],
              newpos[Z], p2[X], p2[Y], p2[Z], a, Ecore, f1, f2);

    for (int i=0; (i<3); i++)
    {
        f0seg2[i] += f1[i];
        f1seg2[i] += f2[i];
    }

    SegSegForce(home, f1,f2,f3,f4, p1,newpos,newpos,p2, burg,burg, a, MU, NU); 

    for (int i=0; (i<3); i++)
    {
        f0seg1[i] += f1[i];
        f1seg1[i] += f2[i];
        f0seg2[i] += f3[i];
        f1seg2[i] += f4[i];
    }

    return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       FindFSegComb
 *      Description:
 *
 *      Arguments
 *          p1       Coordinates of point 1
 *          p2       Coordinates of point 2 (corresponding to the periodic
 *                   image of p2 closest to point p1)
 *          p3       Coordinates of point 3 (corresponding to the periodic
 *                   image of p3 closest to point p2)
 *          burg1    burgers vector from p1 to p2
 *          burg2    burgers vector from p2 to p3
 * WARNING the fp*seg* arrays are modified!
 *          fp1seg1  force of segment p1-->p2 at point p1
 *          fp2seg1  force of segment p1-->p2 at point p2
 *          fp2seg2  force of segment p2-->p3 at point p2
 *          fp3seg2  force of segment p2-->p3 at point p3
 *          fsegnew  resulting forces at p1 and p3 from the
 *                   segment p1-->p3
 *
 *-------------------------------------------------------------------------*/
void FindFSegComb(Home_t *home, real8 p0[3], real8 p1[3], real8 p2[3],
                  real8 burg1[3], real8 burg2[3], real8 fp0seg1[3],
                  real8 fp1seg1[3], real8 fp1seg2[3], real8 fp2seg2[3],
                  real8 f0new[3], real8 f1new[3])
{
    real8    L, Lsub1, Lsub2;
    real8    Ecore, MU, NU, a;
    real8    tmpfp0seg1[3], tmpfp1seg1[3], tmpfp1seg2[3], tmpfp2seg2[3];
    real8    fLinv1a[3], fLinv1b[3], fLinv2a[3], fLinv2b[3];
    real8    f1[3], f2[3], f3[3], f4[3];
    Param_t  *param;

    param = home->param;

    Ecore = param->Ecore;
    MU    = param->shearModulus;
    NU    = param->pois;
    a     = param->rc;

    Lsub1 = (p1[X]-p0[X])*(p1[X]-p0[X]) +
            (p1[Y]-p0[Y])*(p1[Y]-p0[Y]) +
            (p1[Z]-p0[Z])*(p1[Z]-p0[Z]);

    Lsub2 = (p1[X]-p2[X])*(p1[X]-p2[X]) +
            (p1[Y]-p2[Y])*(p1[Y]-p2[Y]) +
            (p1[Z]-p2[Z])*(p1[Z]-p2[Z]);

/*
 *  Find the self force of each segment at each point and subtract
 *  it off.
 */
    SelfForce(home, 0, MU, NU, burg1[X], burg1[Y], burg1[Z], p0[X], p0[Y], p0[Z],
              p1[X], p1[Y], p1[Z], a, Ecore, f1, f2);

    for (int i=0; (i<3); i++)
    {
        tmpfp0seg1[i] = fp0seg1[i] - f1[i];
        tmpfp1seg1[i] = fp1seg1[i] - f2[i];
    }

    SelfForce(home, 0, MU, NU, burg2[X], burg2[Y], burg2[Z], p1[X], p1[Y], p1[Z],
              p2[X], p2[Y], p2[Z], a, Ecore, f1, f2);

    for (int i=0; (i<3); i++)
    {
        tmpfp1seg2[i] = fp1seg2[i] - f1[i];
        tmpfp2seg2[i] = fp2seg2[i] - f2[i];
    }

/*
 *  Find the segment/segment forces between the two segments and subtract
 *  that off also.
 */
    SegSegForce(home, f1,f2,f3,f4, p0,p1,p1,p2, burg1,burg2, a, MU, NU); 

    for (int i=0; (i<3); i++)
    {
        tmpfp0seg1[i] -= f1[i];
        tmpfp1seg1[i] -= f2[i];
        tmpfp1seg2[i] -= f3[i];
        tmpfp2seg2[i] -= f4[i];
    }

    for (int i=0; (i<3); i++)
    {
        fLinv1a[i] = ((4.0 * tmpfp0seg1[i]) - (2.0 * tmpfp1seg1[i])) / Lsub1;
        fLinv1b[i] = ((4.0 * tmpfp1seg1[i]) - (2.0 * tmpfp0seg1[i])) / Lsub1;
        fLinv2a[i] = ((4.0 * tmpfp1seg2[i]) - (2.0 * tmpfp2seg2[i])) / Lsub2;
        fLinv2b[i] = ((4.0 * tmpfp2seg2[i]) - (2.0 * tmpfp1seg2[i])) / Lsub2;
    }

/*
 *  If the two segments are redundant links between two nodes, we have
 *  to do some special handling.
 */
    if ((p0[X] == p2[X]) && (p0[X] == p2[X]) && (p0[X] == p2[X]))
    {
        L = Lsub1;
        for (int i=0; (i<3); i++)
        {
            f0new[i] = (2.0 * (fLinv1a[i] + fLinv2b[i]) +
                        fLinv2a[i]+fLinv1b[i]) * L / 6.0;
            f1new[i] = (2.0 * (fLinv2a[i] + fLinv1b[i]) +
                        fLinv1a[i]+fLinv2b[i]) * L / 6.0;
        }
    }
    else
    {
        L = (p0[X]-p2[X]) * (p0[X]-p2[X]) +
            (p0[Y]-p2[Y]) * (p0[Y]-p2[Y]) +
            (p0[Z]-p2[Z]) * (p0[Z]-p2[Z]);

        for (int i=0; (i<3); i++)
        {
            f0new[i] = ((2.0 * fLinv1a[i]) + fLinv1b[i]) * L / 6.0;
            f1new[i] = (fLinv1a[i] + (2.0 * fLinv2b[i])) * L / 6.0;
        }
    }

    SelfForce(home, 0, MU, NU, burg1[X], burg1[Y], burg1[Z],
              p0[X], p0[Y], p0[Z], p2[X], p2[Y], p2[Z], a, Ecore, f1, f2);

    for (int i=0; (i<3); i++)
    {
        f0new[i] += f1[i];
        f1new[i] += f2[i];
    }

    return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     IncrDomSegCommCnts
 *      Description:  Given the 2 nodes defining a segment, check if
 *                    either node is not native to this domain.  If so,
 *                    add the remote domain(s) to the list of domains
 *                    to which this domain will be sending segment
 *                    force data, and increment the count of segments
 *                    being sent to any remote domains owning either
 *                    of the specified nodes.
 *      Arguments:
 *          node1: IN   Node pointer to segment's first endpoint
 *          node2: IN   Node pointer to segment's second endpoint
 *          sendParticleIntersects: IN Toggle set to 1 if the list of
 *                      particles intersected by this segment should
 *                      be sent to the remote domain, zero otherwise.
 *          sendDomCnt: IN/OUT Pointer to int containing the current count
 *                      of remote domains to which force updates will
 *                      be sent.  Contents of this value will be
 *                      updated as necessary for the caller.
 *          segCommInfo: IN/OUT Array of <numSendBufs> structs containing
 *                      info (domain ID, segment count, count of
 *                      segment/particle intersect structures to send) for
 *                      each remote domain to which this domain will be
 *                      sending seg force data.
 *          msgCnts: IN/OUT Array of ints (1 per domain) in which the
 *                      value for each remote domain will be set to 1 if the
 *                      current domain is sending force updates to the
 *                      remote domain, or zero if not.
 *
 *------------------------------------------------------------------------*/
void IncrDomSegCommCnts(Home_t *home, Node_t *node1, Node_t *node2,
                        int sendParticleIntersects, int *sendDomCnt,
                        SegCommInfo_t **segCommInfo, int *msgCnts)
{
        int            domID1, domID2;
        int            index, maxIndex;
        int            thisDomain, sendCnt;
        SegCommInfo_t  *list;

        thisDomain = home->myDomain;

        sendCnt = *sendDomCnt;
        maxIndex = sendCnt;

        list = *segCommInfo;

        domID1 = node1->myTag.domainID;
        domID2 = node2->myTag.domainID;

        if (domID1 != thisDomain) {

            msgCnts[domID1] = 1;

            for (index = 0; index < maxIndex; index++) {
                if (list[index].domID == domID1) break;
            }

/*
 *          If the remote domain was already on the send list, just increment
 *          the count of segment forces being sent.  Otherwise, add the domain
 *          to the send list.
 */
            if (index < maxIndex) {
                list[index].segCnt += 1;
            } else {
                list = (SegCommInfo_t *)realloc(list, (sendCnt+1) *
                                                      sizeof(SegCommInfo_t));
                list[index].domID = domID1;
                list[index].segCnt = 1;
/*
 *              When the segment is not native (i.e. first node owned by remote
 *              domain) we don't send out any particle intersection info for
 *              the segment.
 */
                list[index].intersectDataCnt = 0;
                sendCnt++;
                maxIndex++;
            }
        }

        if ((domID2 != thisDomain) && (domID2 != domID1)) {

            msgCnts[domID2] = 1;

            for (index = 0; index < maxIndex; index++) {
                if (list[index].domID == domID2) break;
            }

/*
 *          If the remote domain was already on the send list, just increment
 *          the count of segment forces being sent.  Otherwise, add the domain
 *          to the send list.
 */
            if (index < maxIndex) {
                list[index].segCnt += 1;
            } else {
                list = (SegCommInfo_t *)realloc(list, (sendCnt+1) *
                                                      sizeof(SegCommInfo_t));
                list[index].domID = domID2;
                list[index].segCnt = 1;
                list[index].intersectDataCnt = 0;
                sendCnt++;
            }
/*
 *          Bump up the count of segment/particle intersection data we'll
 *          be sending (if segment is so flagged).
 */
            list[index].intersectDataCnt += sendParticleIntersects;
        }

        *sendDomCnt = sendCnt;
        *segCommInfo = list;

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     LocalSegForces
 *      Description:  This is just a driver function for computing the
 *                    the interactions between pairs of nearby dislocation
 *                    segments.  There are two versions of this function.
 *
 *                    The first version is only used when calculating forces
 *                    for full N^2 segment-to-segment interactions.  It
 *                    assumes knowledge of all dislocation segments and no
 *                    remote force calculations are done.  This method is
 *                    only valid with serial compilation and is used
 *                    primarily for debugging or verification purposes.
 *
 *                    The second method is valid for serial or parallel
 *                    execution and builds lists of segments from native
 *                    and neighbor cells, handles selection of segment
 *                    pairs whose forces are to be calculated, and
 *                    communicates calculated segment forces to remote
 *                    domains as necessary.
 *
 *------------------------------------------------------------------------*/
#ifdef FULL_N2_FORCES
void LocalSegForces(Home_t *home, int reqType)
{
        int     i, j, k, l, m;
        int     setSeg1Forces, setSeg2Forces;
        int     armID12, armID21, armID34, armID43;
        int     isIntersecting;
        real8   a, E, MU, NU, Ecore, W=0;
        real8   pos1[3], pos2[3], burg[3];
        real8   dx, dy, dz;
        real8   f1[3], f2[3], f3[3], f4[3];
        real8   sigb[3], extStress[3][3];
        Node_t  *node1, *node2, *node3, *node4;
        Param_t *param;

        param = home->param;

        a     = param->rc;
        E     = param->YoungsModulus;
        MU    = param->shearModulus;
        NU    = param->pois;
        Ecore = param->Ecore;

        extStress[0][0] = param->appliedStress[0];
        extStress[1][1] = param->appliedStress[1];
        extStress[2][2] = param->appliedStress[2];
        extStress[1][2] = param->appliedStress[3];
        extStress[2][0] = param->appliedStress[4];
        extStress[0][1] = param->appliedStress[5];
        extStress[2][1] = extStress[1][2];
        extStress[0][2] = extStress[2][0];
        extStress[1][0] = extStress[0][1];

        int mval = 0;


        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            pos1[X] = node1->x;
            pos1[Y] = node1->y;
            pos1[Z] = node1->z;

            for (j = 0; j < node1->numNbrs; j++) {

                node2 = GetNeighborNode(home, node1, j);

                if (node2 == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                    continue;
                }

/*
 *              Ensure node1 is the node with the lower tag
 */
                if (OrderNodes(node1, node2) >= 0) {
                    continue;
                }

                dx = node2->x - pos1[X];
                dy = node2->y - pos1[Y];
                dz = node2->z - pos1[Z];

                ZImage(param, &dx, &dy, &dz);

/*
 *             It is possible to have a zero-length segment (created by
 *             collision handling).  If we find such a segment, there will
 *             be no forces on the segment, so just skip to the next segment.
 */
                if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
                    continue;
                }

                pos2[X] = pos1[X] + dx;
                pos2[Y] = pos1[Y] + dy;
                pos2[Z] = pos1[Z] + dz;

                burg[X] = node1->burgX[j];
                burg[Y] = node1->burgY[j];
                burg[Z] = node1->burgZ[j];

                armID12 = j;
                armID21 = GetArmID(node2, node1);

                setSeg1Forces = 1;

/*
 *              If we're doing a partial force calculation, only
 *              reset forces for this segment if one of the nodes
 *              is flagged for a force update.
 */
                if (reqType == PARTIAL) {
                    if (((node1->flags & NODE_RESET_FORCES) == 0) &&
                        ((node2->flags & NODE_RESET_FORCES) == 0)) {
                        setSeg1Forces = 0;
                    }
                }

/*
 *              Before calculating the force from other segments,
 *              calculate forces specific to this segment
 */
                if (setSeg1Forces) {
/*
 *                  Add in force due to self stress
 */
                    SelfForce(home, 0, MU, NU, burg[X], burg[Y], burg[Z],
                              pos1[X], pos1[Y], pos1[Z],
                              pos2[X], pos2[Y], pos2[Z],
                              a, Ecore, f1, f2);

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);

#ifdef CALCENERGY
                    W = SelfEnergy(0, MU, NU, burg[X], burg[Y], burg[Z],
				   pos1[X], pos1[Y], pos1[Z],
				   pos2[X], pos2[Y], pos2[Z],
				   a, Ecore);

                    home->param->TotalEnergy += W;
#endif


/*
 *                  Add in force due to external stress
 */
                    ExtPKForce(extStress, burg[X], burg[Y], burg[Z],
                               pos1[X], pos1[Y], pos1[Z],
                               pos2[X], pos2[Y], pos2[Z], f1, f2);

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);

#ifdef CALCENERGY
                    W = ExtPKEnergy(extStress, E, NU, burg[X], burg[Y], burg[Z],
				    pos1[X], pos1[Y], pos1[Z],
				    pos2[X], pos2[Y], pos2[Z]);

                    home->param->TotalEnergy += W;
#endif


#ifdef ESHELBY
/*
 *                  Collect dislocation segments intersecting particles
 *                  Add in force resulting from any eshelby inclusions
 */
                    SegPartIntersect_t *intersection = (SegPartIntersect_t *)NULL;
                    int                 addToList    = (intersection == (SegPartIntersect_t *)NULL ? 1 : 0);

                    for (m = 0; m < home->totInclusionCount; m++)
                    {
                        EInclusion_t *inclusion = &home->eshelbyInclusions[m];

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
#endif //ESHELBYFORCE


                        if (isIntersecting)
                        {
#ifdef ESHELBYCORE
                           SegPartCoreForce(home, inclusion->position, inclusion->radius,
                                            inclusion->rotation, newPos1, newPos2,
                                            burg, planenormal, f1, f2);


                           AddtoArmForce(node1, armID12, f1);
                           AddtoArmForce(node2, armID21, f2);
#endif //ESHELBYCORE



#ifdef ESHELBYSTEP
                           InclusionStepForce(home->param->IncStepEnergy, inclusion->position,
                                              inclusion->radius,inclusion->rotation,
                                              newPos1, newPos2,burg, f1, f2);


                           AddtoArmForce(node1, armID12, f1);
                           AddtoArmForce(node2, armID21, f2);
#endif //ESHELBYSTEP

                        }

/*
 *                      If the segment intersects the inclusion, and we need
 *                      to track such things for the mobility function,
 *                      add the segment/particle intersect info to the
 *                      list.
 */
                           if (isIntersecting && addToList)
                           {

                              SegPartListUpdate(home, &intersection, &node1->myTag,
                                                &node2->myTag, m, inclusion->id);

                           }
                   } // end loop on m totInclusionCount

#endif  // ESHELBY

#ifdef _ARLFEM
/*
 *                  With free surfaces, we also have to calculate the force
 *                  from all "virtual" segments in the simulation space.
 *                  Virtual segments are imaginary segments extending
 *                  outward from any segment intersecting the surface.
 *                  The segments are treated as if they extend from
 *                  the surface to infinity, but are in reality long
 *                  but finite segments.
 */
                    VirtualSegForce(home, pos1, pos2, burg, f1, f2)

                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);

#endif  /* end ifdef _ARLFEM */
                } // end setSeg1Forces

/*
 *              Now compute the force between segment (node1--node2)
 *              and all other segments in the simulation.
 */
                for (k = 0; k < home->newNodeKeyPtr; k++) {

                    if ((node3 = home->nodeKeys[k]) == (Node_t *)NULL) {
                        continue;
                    }

                    for (l = 0; l < node3->numNbrs; l++) {

                        node4 = GetNeighborNode(home, node3, l);

                        if (node4 == (Node_t *)NULL) {
                            printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                            continue;
                        }

/*
 *                      Ensure the node with the lower tag is the node3
 */
                        if (OrderNodes(node3, node4) >= 0) {
                            continue;
                        }

/*
 *                      Make sure we don't try to calculate seg/seg forces
 *                      on a segment with itself, and that we only do
 *                      forces between a given pair once.
 */
                        if ((node1 == node3) && (node2 == node4)) {
                            continue;
                        }

                        if (OrderNodes(node3, node1) < 0) {
                            continue;
                        }

                        if ((OrderNodes(node4, node2) < 0) &&
                            (node1 == node3)) {
                            continue;
                        }

                        setSeg2Forces = 1;
/*
 *                      If we're doing a partial force calculation, only
 *                      reset forces for this segment if one of the nodes
 *                      is flagged for a force update.
 */
                        if (reqType == PARTIAL) {
                            if (((node3->flags & NODE_RESET_FORCES) == 0) &&
                                ((node4->flags & NODE_RESET_FORCES) == 0)) {
                                setSeg2Forces = 0;
                            }
                        }

                        if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) {
                            continue;
                        }

/*
 *                      Calculate the forces between segment (node1--node2)
 *                      and segment (node3--node4).
 */
                        armID34 = l;
                        armID43 = GetArmID(node4, node3);

                        ComputeForces(home, node1,node2,node3,node4, f1,f2,f3,f4);

#ifdef CALCENERGY
                        W = ComputeEnergy(home, node1, node2, node3, node4);
                        home->param->TotalEnergy += W;
#endif

                        if (setSeg1Forces) {
                            AddtoArmForce(node1, armID12, f1);
                            AddtoArmForce(node2, armID21, f2);
                        }




                      if (setSeg2Forces) {
                            AddtoArmForce(node3, armID34, f3);
                            AddtoArmForce(node4, armID43, f4);
                        }

                    }  /* for (l = 0; ... ( */

                }  /* for (k = 0; ... ( */

            }  /* for (j = 0; ... ) */

        }  /* for (i = 0; ... ) */

/*
 *      Forces for all segments have been updated, so just
 *      sum up the segment forces to get the nodal forces.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            node1->fX = 0.0;
            node1->fY = 0.0;
            node1->fZ = 0.0;

            for (j = 0; j < node1->numNbrs; j++) {
                node1->fX += node1->armfx[j];
                node1->fY += node1->armfy[j];
                node1->fZ += node1->armfz[j];
            }

        }

        return;
}
#endif  // #ifdef FULL_N2_FORCES...

#ifndef FULL_N2_FORCES
void LocalSegForces(Home_t *home, int reqType)
{
        int        i, cellNum, cellID, cellSegCnt;
        int        arm;
        int        homeDomain, homeCells, homeNativeCells;
        int        totalNativeSegs;
        int        sendDomCnt, allocSize, numDomains;
        int        setSeg1Forces, setSeg2Forces;
        int        *nativeSegCounts, *totalSegCounts;
        int        nativeSegListCnt = 0;
        real8      MU, NU, a, Ecore, extstress[3][3];
        Node_t     *node1, *node2, *node3, *node4;
        Node_t     *node, *nbr;
        Cell_t     *cell;
        Param_t    *param;
        Segment_t  *segList, **cellSegLists;
        NativeSeg_t   *nativeSegList = NULL;
        SegCommInfo_t *segCommInfo;

        SegmentPairList_t  segPairList;

#ifdef PARALLEL
        int         *globalMsgCnts=0, *localMsgCnts=0;
        int         *sendBufLengths = (int *)NULL;
        real8       **sendBufs = (real8 **)NULL;
        MPI_Request *sendBufLenRequests = (MPI_Request *)NULL;
#endif
#ifdef _ARLFEM
        int         nonZeroGlobalSegCnt;
#endif


        homeCells       = home->cellCount;
        homeNativeCells = home->nativeCellCount;
        homeDomain      = home->myDomain;
        numDomains      = home->numDomains;
        param           = home->param;

        MU    = param->shearModulus;
        NU    = param->pois;
        a     = param->rc;
        Ecore = param->Ecore;

        extstress[0][0] = param->appliedStress[0];
        extstress[1][1] = param->appliedStress[1];
        extstress[2][2] = param->appliedStress[2];
        extstress[1][2] = param->appliedStress[3];
        extstress[2][0] = param->appliedStress[4];
        extstress[0][1] = param->appliedStress[5];
        extstress[2][1] = extstress[1][2];
        extstress[0][2] = extstress[2][0];
        extstress[1][0] = extstress[0][1];

        totalNativeSegs = 0;
        sendDomCnt = 0;
        segCommInfo = (SegCommInfo_t *)NULL;

/*
 *      Allocate and initialize some temporary arrays we'll be needing.
 *
 *      For each cell native to or neighboring the current domain, an
 *      array of unique segments is built.  These arrays contain two classes
 *      of segments; "native" and "ghost" (see descriptions in the code
 *      below).  Any "native" segments in a cell's segment list will
 *      preceed "ghost" segments.  The lists are set up this way to
 *      simplify the process of insuring that forces on each pair
 *      of segments are only evaluated one time.
 *
 *      The cell segment arrays are set up as arrays of Segment_t
 *      structures,  each segment represented by a pair of node
 *      pointers and force component for each node of the segment.
 *      These force values will only represent the force on the
 *      nodes from the seg-seg force calcs done by this domain.
 *      These values will be summed with values from remote domains
 *      to get the final forces on the nodes/segments after all
 *      local calculations are done.
 */
#ifdef PARALLEL
        if (param->numDLBCycles == 0) {
            localMsgCnts  = home->localReduceBuf;
            globalMsgCnts = home->globalReduceBuf;

            for (i = 0; i < numDomains; i++) {
                localMsgCnts[i] = 0;
            }
        }
#endif

        allocSize = sizeof(int) * homeCells;
        nativeSegCounts = (int *)calloc(1, allocSize);
        totalSegCounts  = (int *)calloc(1, allocSize);

        allocSize = sizeof(Segment_t *) * homeCells;
        cellSegLists = (Segment_t **)calloc(1, allocSize);

        for (cellNum=0; cellNum < homeCells; cellNum++) {

            cellID = home->cellList[cellNum];
            cell = LookupCell(home, cellID);

            if (cell->nodeCount == 0) continue;

/*
 *          Need to allocate a segment array large enough
 *          for all segments in the cell; Rather than loop through all
 *          the nodes in the cell counting arms, just multiply the cell
 *          node count by the maximum segments per node to get the largest
 *          possible seg count for the cell.  The array will be larger than
 *          needed, but should save some time.
 */
            allocSize = cell->nodeCount * MAX_NBRS * sizeof(Segment_t);
            cellSegLists[cellNum] = (Segment_t *)malloc(allocSize);
        }

/*
 *      Loop over all native cells adding "native" segments to the cell
 *      segment lists.  We only add native segments in this loop because
 *      all native segments in the array must preceed ghost segments (it
 *      makes things easier later on).  Ghost segments will be added
 *      to the arrays a little later.
 */
        for (cellNum=0; cellNum < homeNativeCells; cellNum++) {

            segList = cellSegLists[cellNum];
            cellID = home->cellList[cellNum];
            cell = LookupCell(home, cellID);
            cellSegCnt = 0;
            node = cell->nodeQ;

/*
 *          Cell "native" segments are segments for which the dominant
 *          (i.e. owning) node of the segment is in the current cell and
 *          domain.
 */
            for ( ; node != (Node_t *)NULL; node = node->nextInCell) {

                if (node->myTag.domainID != homeDomain) continue;

                for (arm = 0; arm < node->numNbrs; arm++) {

                    nbr = GetNeighborNode(home, node, arm);

                    if (nbr == (Node_t *)NULL) {
                        printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                        continue;
                    }

                    if (NodeOwnsSeg(home, node, nbr) == 0) {
                        continue;
                    }

/*
 *                  Initialize all components of the segList element
 *                  since the array was not calloc'ed.
 */
                    segList[cellSegCnt].node1                  = node;
                    segList[cellSegCnt].node2                  = nbr;
                    segList[cellSegCnt].forcesSet              = 0;
                    segList[cellSegCnt].sendParticleIntersects = 0;
                    V3_ZERO(segList[cellSegCnt].f1);
                    V3_ZERO(segList[cellSegCnt].f2);

/*
 *                  If needed, create a lock for each segment
 */
                    INIT_LOCK(&segList[cellSegCnt].segLock);

                    cellSegCnt++;
                }
            }

            nativeSegCounts[cellNum] = cellSegCnt;
            totalSegCounts [cellNum] = cellSegCnt;
            totalNativeSegs += cellSegCnt;
        }

/*
 *      Next add "ghost" segments to cell segment lists for
 *      all cells that are either native cells, or ghost cells
 *      which are NOT dominant over ALL cells native to the domain.
 *
 *      Note: A native cell may in fact only partially overlap
 *      the domain, hence the native cell may also contain ghost
 *      segments.
 *
 *      If there are NO native segments in this domain, there's
 *      no point in adding ghost segments to the lists because
 *      no force calcs will be done by this domain...
 */

        cellNum = ( (totalNativeSegs==0) ? homeCells : 0 );

        for (/* init cellNum above */ ; cellNum < homeCells; cellNum++) {

            cell = LookupCell(home, home->cellList[cellNum]);
            node = cell->nodeQ;
            cellSegCnt = totalSegCounts[cellNum];
            segList = cellSegLists[cellNum];

/*
 *          Cell "ghost" segments are comprised of segments owned by
 *          a "ghost" node.
 */
            for ( ; node != (Node_t *)NULL; node = node->nextInCell) {

                if (node->myTag.domainID == homeDomain) continue;

                for (arm = 0; arm < node->numNbrs; arm++) {

                    nbr = GetNeighborNode(home, node, arm);

                    if (nbr == (Node_t *)NULL) {
                        printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                        continue;
                    }

                    if (NodeOwnsSeg(home, node, nbr) == 0) {
                        continue;
                    }

/*
 *                  Initialize all components of the segList element
 *                  since the array was not calloc'ed.
 */
                    segList[cellSegCnt].node1 = node;
                    segList[cellSegCnt].node2 = nbr;
                    segList[cellSegCnt].forcesSet = 0;
                    segList[cellSegCnt].sendParticleIntersects = 0;
                    V3_ZERO(segList[cellSegCnt].f1);
                    V3_ZERO(segList[cellSegCnt].f2);

/*
 *                  Initialize the lock for each segment structure
 */
                    INIT_LOCK(&segList[cellSegCnt].segLock);

                    cellSegCnt++;
                }
            }

            totalSegCounts[cellNum]  = cellSegCnt;
        }

/*
 *      Okay, the cell segment lists are built; now go through and
 *      build a list of all the native segments for which we
 *      have to calculate forces, and a list of all the segment pairs
 *      for which forces need to be calculated.  (Note: in a segment
 *      pair, it's possible we don't need to get the forces for one
 *      of the segmenmts in the pair.)
 *
 *      Since the force calculation returns forces for all four nodes
 *      defining a pair of segments, we have to be sure we only
 *      compute the forces once for every pair of segments.  Additionally
 *      we don't need to compute forces between segment pairs if neither
 *      of the segments is native (i.e. has a node native to the current
 *      domain.)
 *
 *      The outer loop here only goes through the cells native to the
 *      domain because none of the other cells will have native segments.
 *
 */
        nativeSegList = (NativeSeg_t *)calloc(1, sizeof(NativeSeg_t) * totalNativeSegs);

        for (i = 0; i < homeNativeCells; i++) {
            int  j, k, l;
            int  numNbrCells, cellNativeSegs, cellTotalSegs;

            segList = cellSegLists[i];
            cellNativeSegs = nativeSegCounts[i];
            cellTotalSegs = totalSegCounts[i];

            cellID = home->cellList[i];
            cell = LookupCell(home, cellID);

/*
 *          We'll need seg/seg forces between every native
 *          segment in this cell and all other segments following
 *          the segment in the cell's segment list. (Unless the
 *          the 2nd segment is a ghost segment, and the node owning
 *          the second segment is lower priority (force calc will
 *          then be done by domain owning segment 2)).
 */
            for (j = 0; j < cellNativeSegs; j++) {

                setSeg1Forces = 1;

                node1 = segList[j].node1;
                node2 = segList[j].node2;

/*
 *              If we're only doing a partial force calc, we don't
 *              need to update the forces for a native segment if
 *              is not attached to a node marked for update.
 */
                if (reqType == PARTIAL) {

                    if (((node1->flags & NODE_RESET_FORCES) == 0) &&
                        ((node2->flags & NODE_RESET_FORCES) == 0)) {
                        setSeg1Forces = 0;
                    }
                }

/*
 *              If necessary, add native segment to list of native segs for
 *              which we need to calculate forces from interactions
 *              other than seg/seg interactions (i.e. self force,
 *              osmotic force, remote force, etc.)
 *
 *              Note: don't bother adding the segment to the list if we're
 *              doing a load-balance only step.
 */
                if (setSeg1Forces) {

                    if (param->numDLBCycles == 0) {
                        nativeSegList[nativeSegListCnt].seg = &segList[j];
                        nativeSegList[nativeSegListCnt].cell = cell;
                        nativeSegList[nativeSegListCnt].cellID = cellID;

                        nativeSegListCnt++;
                    }
                }

/*
 *              Now for segment pairs for which interactions must
 *              be computed.
 */
                for (k = j + 1; k < cellTotalSegs; k++) {

                    setSeg2Forces = 1;

                    node3 = segList[k].node1;
                    node4 = segList[k].node2;

/*
 *                  If we're only doing a partial force calc, we won't
 *                  need forces for segment 2 if it is not attached to a node
 *                  marked for update.
 */
                    if (reqType == PARTIAL) {

                        if (((node3->flags & NODE_RESET_FORCES) == 0) &&
                            ((node4->flags & NODE_RESET_FORCES) == 0)) {
                            setSeg2Forces = 0;
                        }
                    }

/*
 *                  If neither of the segments needs forces updated,
 *                  skip it.
 */
                    if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) {
                        continue;
                    }

/*
 *                  If segment 2 is a ghost, only do the force calc if the
 *                  node owning segment 1 has lower priority than the
 *                  node owning segment 2.
*/
                    if ((k >= nativeSegCounts[i]) &&
                        (NodeOwnsSeg(home, node1, node3))) {
                        continue;
                    }

/*
 *                  We need forces for this segment pair, so add the
 *                  segment pair to the seg pair list, and increment
 *                  the number of seg/seg force calculations being done
 *                  this step.
 *
 *                  Note: don't bother adding the segment pair to the
 *                  list if we're doing a load-balance only step
 */
                    home->cycleForceCalcCount++;

                    if (param->numDLBCycles == 0)
                    {
                        segPairList.Append(&segList[j], &segList[k], setSeg1Forces, setSeg2Forces);
                    }
                }

            }  /* Loop over native segments */

/*
 *          Next loop over all the neighbors of the current
 *          native cell.  If the current cell has priority
 *          over the the neighboring cell, we do need force
 *          calcs between these pairs in this loop; the segment
 *          pair will either be handled by one of the other
 *          iterations of this loop, or by the remote domain
 *          owning the segments in the other cell.
 *
 *          Note: Cell ids used here convert to ranges from
 *          zero -> num[XYZ]cells+1 allowing for periodic cells.
 *          But nbrCellID gets converted to the base cell index
 *          if it is a periodic cell.
 *
 */
            numNbrCells = cell->nbrCount;

            for (j = 0; j < numNbrCells; j++) 
            {
                int     nbrCellID = cell->nbrList[j];
                Cell_t *nbrCell   = LookupCell(home, nbrCellID);

                if (nbrCell->baseIdx >= 0) {
                    nbrCellID = nbrCell->baseIdx;
                }

/*
 *              If the neighbor cell has priority over the
 *              current cell we need to calculate seg/seg forces
 *              between native segs in the current cell and the
 *              segments in the neighboring cell.
 */
                if (CellPriority(home, cellID, nbrCellID) >= 0) {
                    continue;
                }

                for (k = 0; k < homeCells; k++) {
                    if (nbrCellID == home->cellList[k]) {
                        break;
                    }
                }

                Segment_t *nbrSegList          = cellSegLists   [k];
                int        nbrCellSegCnt       = totalSegCounts [k];
                int        nbrCellNativeSegCnt = nativeSegCounts[k];

/*
 *              If there are no segments in the neighboring cell, no
 *              need to do anything more with this neighbor cell.
 */
                if (nbrCellSegCnt < 1) {
                    continue;
                }

                for (k = 0; k < cellNativeSegs; k++) {

                    node1 = segList[k].node1;
                    node2 = segList[k].node2;

                    setSeg1Forces = 1;

/*
 *                  If we're only doing a partial force calc, we don't
 *                  need forces for segment 1 if it is not attached to a node
 *                  marked for update.
 */
                    if (reqType == PARTIAL) {

                        if (((node1->flags & NODE_RESET_FORCES) == 0) &&
                            ((node2->flags & NODE_RESET_FORCES) == 0)) {
                            setSeg1Forces = 0;
                        }
                    }

                    for (l = 0; l < nbrCellSegCnt; l++) {

                        node3 = nbrSegList[l].node1;
                        node4 = nbrSegList[l].node2;

                        setSeg2Forces = 1;

/*
 *                      If we're only doing a partial force calc, we don't
 *                      need forces for segment 2 if it is not attached to
 *                      a node marked for update.
 */
                        if (reqType == PARTIAL) {

                            if (((node3->flags & NODE_RESET_FORCES) == 0) &&
                                ((node4->flags & NODE_RESET_FORCES) == 0)) {
                                setSeg2Forces = 0;
                            }
                        }

/*
 *                      If the 2nd segment is native, we probably need the
 *                      forces, but if segment 2 is a ghost, only do the
 *                      calc if the node owning segment 1 has lower priority
 *                      than the node owning segment 2.
 */
                        if ( (l >= nbrCellNativeSegCnt) &&
                             (NodeOwnsSeg(home, node1, node3))) {
                            continue;
                        }

                        if ((setSeg1Forces == 0) && (setSeg2Forces == 0)) {
                            continue;
                        }

/*
 *                      We need forces for this segment pair, so add the
 *                      segment pair to the seg pair list, and increment
 *                      the number of seg/seg force calculations being done
 *                      this step.
 *
 *                      Note: don't bother adding the segment pair to the
 *                      list if we're doing a load-balance only step
 */
                        home->cycleForceCalcCount++;

                        if (param->numDLBCycles == 0)
                        {
                            segPairList.Append(&segList[k], &nbrSegList[l], setSeg1Forces, setSeg2Forces);
                        }
                    }
                }
            }  /* for (j = 0; j < numNbrCells...) */
        } /* for (i = 0; i < homeCells...) */

#ifdef _ARLFEM
/*
 *      If we're not doing a load-balancing-only step, send the segment
 *      list to the FEM code so the image forces can be calculated while
 *      ParaDiS is calculating other forces.
 *
 *      Since the send to FEM is a collective operation, don't include
 *      barrier time in timers.
 */
        if (param->numDLBCycles == 0) {
            TimerStop(home, LOCAL_FORCE);
            TimerStop(home, CALC_FORCE);
            nonZeroGlobalSegCnt = SendSegToFEM(home, nativeSegListCnt,
                                               nativeSegList, reqType);
            TimerStart(home, CALC_FORCE);
            TimerStart(home, LOCAL_FORCE);
        }

#endif  /* ifdef _ARLFEM */

/*
 *      Okay, we have explicit lists of all the native segments for
 *      which we need to calculate forces, plus a list of all the
 *      segment pairs for which we need to do seg/seg forces.
 *      Now loop over each list independently calculating the
 *      appropriate forces.  If openmp is being used, loops will
 *      be parallelized automatically.
 *
 *      NOTE: Since each iteration of this loop handles a unique
 *            segment, we technically don't need to do any locking
 *            of the nodal segment forces when we update them because
 *            there is no chance of a race condition.
 */


/*
 *      Accumulate forces on dislocation segments
 */

#ifdef _OPENMP
#pragma omp parallel for private(node1, node2, cell, cellID)
#endif
        for (i = 0; i < nativeSegListCnt; i++) {
            int   j;
            int   armID12, armID21;
            real8 x1, y1, z1;
            real8 x2, y2, z2;
            real8 bx1, by1, bz1;
            real8 dx, dy, dz;
            real8 f1[3], f2[3];
            real8 node1SegForce[3];
            real8 node2SegForce[3];

/*
 *          Zero out some arrays in which we'll accumulate nodal forces
 *          for the segment until we're done with the segment.  Since the
 *          AddtoArmForces() function locks the nodal info for update,
 *          we avoid calling that function until the end of the loop
 */
            V3_ZERO(node1SegForce);
            V3_ZERO(node2SegForce);

/*
 *          Before we calculate seg/seg forces, calculate all
 *          the forces affecting this single native segment
 *          (i.e. self force, pk force, etc) so those forces
 *          will be included when the seg/seg forces are
 *          communicated to remote domains.
 *
 *          Note: if fmm is enabled, calculate the remote force on
 *          each native segment here otherwise, the remote sigb
 *          has previously been computed and force will be calculated
 *          using that data.
 *
 *          This assumes node1 owns the segment!
 */
            node1 = nativeSegList[i].seg->node1;
            node2 = nativeSegList[i].seg->node2;

            x1 = node1->x;
            y1 = node1->y;
            z1 = node1->z;

            armID12 = GetArmID(node1, node2);
            armID21 = GetArmID(node2, node1);

            bx1 = node1->burgX[armID12];
            by1 = node1->burgY[armID12];
            bz1 = node1->burgZ[armID12];

            dx = node2->x - x1;
            dy = node2->y - y1;
            dz = node2->z - z1;

            ZImage(param, &dx, &dy, &dz);

/*
 *          It is possible to have a zero-length segment (created by
 *          collision handling).  If we find such a segment, there will
 *          be no forces on the segment, so just skip to the next segment.
 */
            if ((dx*dx + dy*dy + dz*dz) < 1.0e-20) {
                continue;
            }

            x2 = x1 + dx;
            y2 = y1 + dy;
            z2 = z1 + dz;

/*
 *          Add in force due to self stress
 */
            SelfForce(home, 0, MU, NU, bx1, by1, bz1, x1, y1, z1, x2, y2, z2, a, Ecore, f1, f2);
            V3_ACCUM(node1SegForce, f1);
            V3_ACCUM(node2SegForce, f2);
/*
 *          Shouldn't need to lock access while calculating non-seg/seg
 *          forces on native segments because each segment is unique
 *          and will only be handled in a single thread... so no race
 *          conditions should exist.
 *
 *          Still need to lock access when updating forces for individual
 *          nodes, but that will be handled in AddtoArmForce().
 */
            for (j = 0; j < 3; j++) {
                nativeSegList[i].seg->f1[j] += f1[j];
                nativeSegList[i].seg->f2[j] += f2[j];
            }


#ifdef CALCENERGY
            real8 W = SelfEnergy(0, MU, NU, bx1, by1, bz1, x1, y1, z1, x2, y2, z2, a, Ecore);
            home->param->TotalEnergy += W;
#endif


/*
 *          Add in PK force from external stress
 */
            ExtPKForce(extstress, bx1, by1, bz1, x1, y1, z1, x2, y2, z2, f1, f2);

            V3_ACCUM(node1SegForce, f1);
            V3_ACCUM(node2SegForce, f2);

            for (j = 0; j < 3; j++) {
                nativeSegList[i].seg->f1[j] += f1[j];
                nativeSegList[i].seg->f2[j] += f2[j];
            }

#ifdef CALCENERGY
            real8 E = param->YoungsModulus;
            W = ExtPKEnergy(extstress, E, NU, bx1, by1, bz1, x1, y1, z1, x2, y2, z2);
            home->param->TotalEnergy += W;
#endif

/*
 *          If we're including osmotic forces, add those in now
 */
            if (param->vacancyConcEquilibrium > 0.0) {

                OsmoticForce(home, x1, y1, z1, x2, y2, z2,
                             bx1, by1, bz1, f1, f2);

                V3_ACCUM(node1SegForce, f1);
                V3_ACCUM(node2SegForce, f2);

                for (j = 0; j < 3; j++) {
                    nativeSegList[i].seg->f1[j] += f1[j];
                    nativeSegList[i].seg->f2[j] += f2[j];
                }
            }

/*
 *          Add in remote force from fast-multipole method.
 */
#ifdef ESHELBY
            if ((param->fmEnabled) || (param->eshelbyfmEnabled))
#else
            if  (param->fmEnabled)
#endif
            {
                RemoteForceOneSeg(home, node1, node2, f1, f2);

                V3_ACCUM(node1SegForce, f1);
                V3_ACCUM(node2SegForce, f2);

                for (j = 0; j < 3; j++) {
                    nativeSegList[i].seg->f1[j] += f1[j];
                    nativeSegList[i].seg->f2[j] += f2[j];
                }
            }

#ifdef SPECTRAL
            if (param->FFTenabled)
            {
                SpectralForceOneSeg(home, node1, node2, f1, f2);

                V3_ACCUM(node1SegForce, f1);
                V3_ACCUM(node2SegForce, f2);

                for (j = 0; j < 3; j++) {
                    nativeSegList[i].seg->f1[j] += f1[j];
                    nativeSegList[i].seg->f2[j] += f2[j];
                }
            }
#endif

#ifndef _ARLFEM
#ifdef SPECTRAL
            if ( (param->fmEnabled==0) && (param->FFTenabled==0) )
#else
            if   (param->fmEnabled==0)
#endif
#endif
            {
                real8 sigb[3];

                sigb[0] = node1->sigbRem[armID12*3];
                sigb[1] = node1->sigbRem[armID12*3+1];
                sigb[2] = node1->sigbRem[armID12*3+2];
/*
 *              Add in PK force from non-fmm remote stress
 */
                PKForce(sigb, x1, y1, z1, x2, y2, z2, f1, f2);

                V3_ACCUM(node1SegForce, f1);
                V3_ACCUM(node2SegForce, f2);

                for (j = 0; j < 3; j++) {
                    nativeSegList[i].seg->f1[j] += f1[j];
                    nativeSegList[i].seg->f2[j] += f2[j];
                }
            }

#ifdef _ARLFEM
/*
 *          With free surfaces, we also have to calculate the force from
 *          all "virtual" segments in the simulation space.  Virtual segments
 *          are imaginary segments extending outward from any segment
 *          intersecting the surface.  The segments are treated as if they
 *          extend from the surface to infinity, but are in reality long
 *          but finite segments.
 */
            {
                real8 pos1[3], pos2[3], burg[3];

                pos1[0] = x1;
                pos1[1] = y1;
                pos1[2] = z1;

                pos2[0] = x2;
                pos2[1] = y2;
                pos2[2] = z2;

                burg[0] = bx1;
                burg[1] = by1;
                burg[2] = bz1;

                VirtualSegForce(home, pos1, pos2, burg, f1, f2);

                V3_ACCUM(node1SegForce, f1);
                V3_ACCUM(node2SegForce, f2);

                for (j = 0; j < 3; j++) {
                    nativeSegList[i].seg->f1[j] += f1[j];
                    nativeSegList[i].seg->f2[j] += f2[j];
                }
            }
#endif  /* end ifdef _ARLFEM */


#ifdef ESHELBY
/*
 *      Calculate interaction force between particles and dislocation segment
 */
        if (param->enableInclusions > 0)
        {
           int                inclusionIndex, inclusionCellID;
           int                numNbrCells;
           Cell_t             *inclusionCell;
           EInclusion_t       *inclusion;

           cell   = nativeSegList[i].cell;
           cellID = nativeSegList[i].cellID;

           numNbrCells = cell->nbrCount;

           SegPartIntersect_t *intersection = (SegPartIntersect_t *)NULL;
           int                 addToList    = (intersection == (SegPartIntersect_t *)NULL ? 1 : 0);


/*
 *         First loop over all neighboring cells, then on
 *         an extra loop iteration do the current cell.
 */
           for (j = 0; j <= numNbrCells; j++)
           {

              if (j < numNbrCells)
                 inclusionCellID = cell->nbrList[j];
              else
                 inclusionCellID = cellID;

/*
 *            Loop through all the inclusions in the cell
 *            and compute the force on the native segment from
 *            each inclusion.  NOTE: If the cell is a ghost
 *            cell, alter the cellID and cell pointer to the
 *            corresponding base cell.  Needed in order to pick
 *            up the proper inclusion.
 */
              inclusionCell = LookupCell(home, inclusionCellID);

              if (inclusionCell->baseIdx >= 0)
              {
                 inclusionCellID = inclusionCell->baseIdx;
                 inclusionCell = LookupCell(home, inclusionCellID);
              }

              inclusionIndex = inclusionCell->inclusionQ;

              int k;
              real8 pos1[3], pos2[3], burg[3];

              pos1[X] = node1->x;
              pos1[Y] = node1->y;
              pos1[Z] = node1->z;

              pos2[X] = node2->x;
              pos2[Y] = node2->y;
              pos2[Z] = node2->z;

              burg[X] = bx1;
              burg[Y] = by1;
              burg[Z] = bz1;


              while (inclusionIndex >= 0)
              {
                 inclusion = &home->eshelbyInclusions[inclusionIndex];

                 int ninc = 0, isIntersecting;
                 real8 newPos1[3], newPos2[3];
                 real8 ratio[2];
                 isIntersecting = IntersectionInformation(home, pos1, pos2, inclusion,
                                                          newPos1, newPos2, &ninc,
                                                          ratio);

#ifdef ESHELBYFORCE
                 // Eshelby force applies on/effects all segments in the simulation and not on
                 // only the segments intersecting the particles!
                 EshelbyForce(home, inclusion->position, inclusion->radius,
                              inclusion->strainField, newPos1, newPos2,
                              burg, f1, f2);

                 V3_ACCUM(node1SegForce,f1);
                 V3_ACCUM(node2SegForce,f2);

                 for (k = 0; k < 3; k++)
                 {
                    nativeSegList[i].seg->f1[k] += f1[k];
                    nativeSegList[i].seg->f2[k] += f2[k];
                 }
#endif   //ESHELBYFORCE


                 if (isIntersecting)
                 {
#ifdef ESHELBYCORE
                    // Core force or modulus mismatch force applies to segments intersecting the
                    // particles and segments completely inside the particle.
                    SegPartCoreForce(home, inclusion->position, inclusion->radius,
                                     inclusion->rotation, newPos1, newPos2,
                                     burg, f1, f2);

                    V3_ACCUM(node1SegForce,f1);
                    V3_ACCUM(node2SegForce,f2);

                    for (k = 0; k < 3; k++)
                    {
                       nativeSegList[i].seg->f1[k] += f1[k];
                       nativeSegList[i].seg->f2[k] += f2[k];
                    }

#endif //ESHELBYCORE


#ifdef ESHELBYSTEP
                    InclusionStepForce(home->param->IncStepEnergy, inclusion->position, inclusion->radius,
                                       inclusion->rotation, newPos1, newPos2,
                                       burg, f1, f2);


                    AddtoArmForce(node1, armID12, f1);
                    AddtoArmForce(node2, armID21, f2);

                    for (k = 0; k < 3; k++)
                    {
                       nativeSegList[i].seg->f1[k] += f1[k];
                       nativeSegList[i].seg->f2[k] += f2[k];
                    }

#endif //ESHELBYSTEP
                 }


/*
 *                  Store the segment that is intersecting the particle
 */
                   if (isIntersecting && addToList)
                   {
/*
 *                    If the segment intersects the inclusion
 *                    add the segment segment/particle intersection
 *                    list.
 */
                      SegPartListUpdate(home, &intersection,
                                        &node1->myTag, &node2->myTag,
                                        inclusionIndex, inclusion->id);

		      //printf("inclusionIndex=%d inclusion->id=%d\n",inclusionIndex,inclusion->id);
                   }

                   inclusionIndex = inclusion->nextInCell;

                 }  /* end while loop through cell's inclusions */

              }  /* end loop over neighbor cells */
        }  /* end if inclusions enabled */
#endif // ESHELBY


/*
 *          We've accumulated various forces on the segment, so now
 *          increment the segment forces (locking will be done within
 *          the update function).
 */
            AddtoArmForce(node1, armID12, node1SegForce);
            AddtoArmForce(node2, armID21, node2SegForce);

            nativeSegList[i].seg->forcesSet = 1;

        }  /* end for (i = 0; i < nativeSegCnt; i++) */

#ifdef ESHELBY
/*
 *          If we're tracking segment/particle intersections for later use in
 *          the mobility function, we need to sort the local list once now
 *          so we can do quick lookups on the list when we send the segment
 *          forces (and particle intersection data) to the remote domains.
 *          After the communication, the list will still have to be sorted
 *          again since data may have been received from remote domains
 *          and appended to the list.
 */
        if (param->numDLBCycles == 0)
        {
            SegPartListSort(home);
        }
#endif



/*
 *      Now loop through every segment pair for which forces must
 *      be calculated.  Note: it may be that we only update forces
 *      for one segment in a given pair because the forces on the
 *      other segment are not needed.
 *
 *      Note:  For threaded runs, updates to the forces in the segment
 *             structure must be explicitly synchronized via locking.
 *             since portions of the forces for a given segment may
 *             be calculated simultaneously in multiple threads.
 *             Updates to forces in the individual node structures
 *             will be synchronized in AddToArmForce().
 */
/*
 * FIX ME!  It might be better to use a pure parallel code block and
 *          calculate which interations each thread does.  Then we
 *          might be able to hold off locking and updating segment
 *          forces if the same thread will be calculating another
 *          set of forces on the same segment the next loop iteration...
 */

        {
            SegmentPair_t *ssf_sp   = segPairList.spl_pairs;
            int            ssf_np   = segPairList.spl_cnt  ;
#if defined(ANISOTROPIC)
            int            ssf_mode = 1;   // (anisotropic)
#else
            int            ssf_mode = 0;   // (isotropic  )
#endif

//#define USE_NEW_SSF_FORCES
#ifdef  USE_NEW_SSF_FORCES
            if ( ssf_sp && (ssf_np>0) ) { SSF_Compute_Forces (home,ssf_sp,ssf_np,ssf_mode); }
#else
            if ( ssf_sp && (ssf_np>0) ) {     ComputeForces  (home,ssf_sp,ssf_np); }
#endif
        }

/*
 *      If we're not doing a load-balance only step, we'll need to
 *      do some additional stuff including sending resulting forces
 *      out to other domains.
 */
        if (param->numDLBCycles == 0)
        {
#ifdef _ARLFEM
/*
 *      Read the force due to image stress from the FEM code
 *      if needed.
 *
 *      Since the receive from the FEM is a collective operation, don't
 *      include barrier time in timers.
 */
            if (param->numDLBCycles == 0)
            {
                if (nonZeroGlobalSegCnt)
                {
                    TimerStop(home, LOCAL_FORCE);
                    TimerStop(home, CALC_FORCE );

                    AddImageForce(home,nativeSegListCnt,nativeSegList);

                    TimerStart(home, CALC_FORCE );
                    TimerStart(home, LOCAL_FORCE);
                }
            }

#endif  /* ifdef _ARLFEM */

#ifdef PARALLEL
/*
 *          Bump up the count of segments that will be sent to
 *          any remote domain owning one of the nodes in any
 *          segment whose forces were updated by this domain.
 */
            for (cellNum = 0; cellNum < homeCells; cellNum++)
            {
                cellSegCnt = totalSegCounts[cellNum];
                segList    = cellSegLists  [cellNum];

                for (i = 0; i < cellSegCnt; i++)
                {
                    if (segList[i].forcesSet == 1)
                    {
                        IncrDomSegCommCnts(home, segList[i].node1, segList[i].node2,
                                           segList[i].sendParticleIntersects,
                                           &sendDomCnt, &segCommInfo,
                                           localMsgCnts);

                    }
                }
            }

/*
 *          Now we need to communicate the newly computed segment forces
 *          to all the appropriate remote domains.  We want the 'local'
 *          force timer to include the time spent packing up buffers and
 *          initiating the sends, but do not want to include the time spent
 *          in the all-reduce or waiting for buffers to arrive.  We do
 *          want to track the all-reduce/recv time separately, though.
 */
            CommSendSegForces(home, sendDomCnt, segCommInfo, cellSegLists,
                              totalSegCounts, &sendBufs, &sendBufLengths,
                              &sendBufLenRequests);

            TimerStop (home, LOCAL_FORCE  );
            TimerStop (home, CALC_FORCE   );
            TimerStart(home, SEGFORCE_COMM);

            MPI_Allreduce(localMsgCnts, globalMsgCnts, numDomains, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            CommRecvSegForces(home, globalMsgCnts[homeDomain], sendDomCnt, &sendBufs,
                              &sendBufLengths, &sendBufLenRequests);

            TimerStop (home, SEGFORCE_COMM);
            TimerStart(home, CALC_FORCE   );
            TimerStart(home, LOCAL_FORCE  );
#endif

/*
 *          We should now have updated forces for nodes/segments
 *          so now do a quick loop through all local nodes and set
 *          the nodes' total forces to the sum of their arms' forces.
 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (i=0; (i<home->newNodeKeyPtr); i++)
            {
                Node_t *node = home->nodeKeys[i];

                if (node)
                {
                    node->fX = 0.0;
                    node->fY = 0.0;
                    node->fZ = 0.0;

                    for (int j=0; (j<node->numNbrs); j++)
                    {
                        node->fX += node->armfx[j];
                        node->fY += node->armfy[j];
                        node->fZ += node->armfz[j];
                    }
                }
            }
        }  /* end if (param->numDLBCycles == 0) */

/*
 *      Free all temporary arrays
 */
        for (cellNum = 0; cellNum < homeCells; cellNum++)
        {
            segList = cellSegLists[cellNum];

            if (segList)
            {
#ifdef _OPENMP
                for (i=0; (i<totalSegCounts[cellNum]); i++)
                { DESTROY_LOCK(&segList[i].segLock); }
#endif
                free(segList);
            }
        }

        if (nativeSegList  ) free(nativeSegList  );

        if (cellSegLists   ) free(cellSegLists   );
        if (nativeSegCounts) free(nativeSegCounts);
        if (totalSegCounts ) free(totalSegCounts );

        if (segCommInfo    ) free(segCommInfo    );

        return;
}
#endif  /* FULL_N2_FORCES not defined */
