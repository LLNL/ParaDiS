/***************************************************************************
 *
 *      Module:       CrossSlip.c
 *      Description:  This module just provides a generic dispatch function
 *                    that will determine if cross-slip is enabled and
 *                    invoke a cross-slip function appropriate to the
 *                    material type (i.e. BCC, FCC, etc).  Also included 
 *                    are some support functions used by the material-specific
 *                    cross-slip functions.
 *
 *                    Although the exact implementation of the cross-slip
 *                    mechanism may not be the same for the different
 *                    material types, the different cross-slip mechanism
 *                    all share the same basic design.
 *
 *                    Cross slip is treated as a topological operation because
 *                    it modifies the connectivity data.  It does not introduce
 *                    new nodes and change any topology, but it can potentially
 *                    change the glide plane information in the connections
 *                    and effect slight changes in the nodal positions
 *
 *                    The cross slip routines consider the cross slip of
 *                    discretization nodes that are close to screw.  Closeness
 *                    to screw is determined by a critical angle (1 degree is
 *                    the value being used at the time this comment was added).
 *                    Both of the segments connected to the node being
 *                    considered must be within this critical angle as well
 *                    as the virtual line connecting the its two neighbor
 *                    segments.
 * 
 *                    If the segments connected to a node are considered close
 *                    to screw then the node is considered for a possible
 *                    cross slip operation.  A test is conducted to determine
 *                    which glide direction of the screw in its three glide
 *                    planes sees the greatest projection of force.  A
 *                    threshold is defined so a node's preference is to
 *                    remain on the primary (current) plane, so only if the
 *                    projection of the force is greatest on a glide plane
 *                    other than the primary plane, and the force on the
 *                    cross-slip plane exceeds the threshold, is a cross-slip
 *                    event attempted.
 *
 *                    There are two possibilities for cross-slip that
 *                    must be considered.
 *
 *                      a) both segments are on same plane (classic case)
 *                      b) the segments are on two different planes with one
 *                         plane being the intended cross slip plane (we
 *                         call this a zipper)
 *
 *                    For case a)  either one or both neighboring nodes are
 *                    moved into perfect screw alignment, the segments are
 *                    flipped to the cross slip plane and the cross slip node
 *                    is slightly moved into the cross slip plane and a small
 *                    areal nucleus is created.  For case b)  the node on the
 *                    primary plane is moved into screw alignment with the
 *                    cross slipping node, the segment is flipped into the
 *                    cross slip plane and the node is moved into the cross
 *                    slip plane. 
 *
 *                    There are conditions which will prevent the cross-slip
 *                    from occuring; in particular:
 *
 *                      1) The processor must own the segments whose data will
 *                         be modified
 *                      2) The processor must own the nodes whose positions
 *                         may be modified
 *                      3) No cross-slip events that would alter a 'pinned'
 *                         node (see the NodePinned() function) are permitted.
 *
 *                    Additionally when a node is physically cross-slipped,
 *                    a new force calculation is performed to see if the 
 *                    cross-slip nucleus continues to move outward into the
 *                    cross slip plane or whether it tends to move back and
 *                    cause a cross slip flicker.  If a flicker is detected
 *                    then the cross slip event is backed out and the node
 *                    is returned to it's original state.  Plus, when a 
 *                    'zipper' event is anticipated, the force on that
 *                    individual segment being affected is examined.  If the
 *                    force on that segment is higher in the cross-slip
 *                    plane than the primary plane, the zipper event is
 *                    performed, otherwise it will be skipped.
 *                    
 *                    NOTE: The original cross-slip mechanisms looped
 *                    through the nodes evaulating each node and performing
 *                    a cross-slip event if necessary before moving on to
 *                    the next node.  This behaviour was modified when the
 *                    cross-slip modules were threaded.  Now, all nodes
 *                    are evaluated in parallel, and when a cross-slip event
 *                    is deemed appropriate for a node, all information
 *                    required to perform that cross-slip is added to a
 *                    thread-shared list of potential cross-slip events.
 *                    This list is built such that zipper-type cross-slip
 *                    events and non-zipper cross-slips are separated.
 *                    Once the list is built, a single thread will loop 
 *                    first through the zipper-type events then the non-
 *                    zipper cross-slips, performing all cross-slip events
 *                    that do not conflict with a prior cross-slip event.
 *
 *      Includes functions:
 *          CrossSlip()
 *          DumpCrossSlipEvent()
 *          PrintCrossSlipEventList()
 *          ProceedWithCrossSlip()
 *          ProceedWithZipEvent()
 *          ResetPosition()
 *          RestoreCrossSlipForce()
 *          SaveCrossSlipInfo()
 *          
 *                    
 ***************************************************************************/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Mobility.h"
#include "CrossSlip.h"

#ifdef DEBUG_CROSSSLIP_EVENTS
void DumpCrossSlipEvent(Node_t *node, real8 newPlane[3], char *eventDesc)
{
        int i;

        printf("Cross-slip (%d,%d):  %s\n", node->myTag.domainID,
               node->myTag.index, eventDesc);

        printf("  node(%d,%d) arms %d, ", node->myTag.domainID,
               node->myTag.index, node->numNbrs);

        for (i = 0; i < node->numNbrs; i++) {
            printf("(%d,%d) ", node->nbrTag[i].domainID,
                   node->nbrTag[i].index);
        }
        printf("\n");

#if 1
/*
 *      Print the nodal position
 */
        printf("  node(%d,%d) position = (%.15e %.15e %.15e)\n",
               node->myTag.domainID, node->myTag.index,
               node->x, node->y, node->z);
#endif

#if 1
/*
 *      Print the burger's vector for each arm of the node
 */
        for (i = 0; i < node->numNbrs; i++) {
            printf("  node(%d %d) arm[%d]-> (%d %d) b = (%.15e %.15e %.15e)\n",
                   node->myTag.domainID, node->myTag.index, i,
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->burgX[i], node->burgY[i], node->burgZ[i]);
        }
#endif

#if 1
/*
 *      Print the old and new glide plane normals for each arm of the node
 *      New plane normal should arleady be in the lab frame, so no rotation
 *      is needed.
 */
        for (i = 0; i < node->numNbrs; i++) {
            printf("  segment(%d %d)--(%d %d) old glide plane = (%f %f %f)\n",
                   node->myTag.domainID, node->myTag.index,
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->nx[i], node->ny[i], node->nz[i]);
            printf("  segment(%d %d)--(%d %d) new glide plane = (%f %f %f)\n",
                   node->myTag.domainID, node->myTag.index,
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   newPlane[X], newPlane[Y], newPlane[Z]);
        }
#endif

        return;
}
#endif

/*---------------------------------------------------------------------------
 *
 *      Function:       PrintCrossSlipEventList
 *
 *      Description:    Prints out the current contents of the cross-slip
 *                      event list.  Prints the 'zipper' cross-slip events
 *                      first, then the normal cross-slip events.
 *                      This function is for debug purposes only.
 *
 *                      NOTE: The zipper-type cross-slip events are
 *                            contained in the first <zipEventCnt>
 *                            structures of <eventList>, while the
 *                            non-zipper cross-slip events are contained
 *                            in the last <eventListSize> - <slipEventIndex>
 *                            structures of <eventList>
 *
 *
 *      Arguments:
 *          eventList      Pointer to the cross-slip event list
 *          zipEventCnt    Number of zipper-type cross-slip events on
 *                         the <eventList>
 *          slipEventIndex Index in <eventList> of the first non-zipper
 *                         type cross-slip event
 *          eventListSize  Size of the <eventList> array as a count of
 *                         EventList_t structures.
 *
 *-------------------------------------------------------------------------*/
void PrintCrossSlipEventList(SlipEvent_t *eventList, int zipEventCnt,
                             int slipEventIndex, int eventListSize)
{
        int i;
        SlipEvent_t *event;

        for (i = 0; i < zipEventCnt; i++) {
            event = &eventList[i];
            printf("Event[%d]: Zip (%d,%d) -- (%d,%d) status = %d\n", i,
                   event->node->myTag.domainID,
                   event->node->myTag.index,
                   event->nbrNodes[event->segIndex]->myTag.domainID,
                   event->nbrNodes[event->segIndex]->myTag.index,
                   event->status);
        }

        for (i = slipEventIndex; i < eventListSize; i++) {
            event = &eventList[i];
            printf("Event[%d]: Non-zip (%d,%d) -- (%d,%d) -- (%d,%d) "
                   "status = %d\n", i,
                   event->nbrNodes[0]->myTag.domainID,
                   event->nbrNodes[0]->myTag.index,
                   event->node->myTag.domainID,
                   event->node->myTag.index,
                   event->nbrNodes[1]->myTag.domainID,
                   event->nbrNodes[1]->myTag.index,
                   event->status);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ProceedWithZipEvent
 *
 *      Description:    This function will determine if the current
 *                      zipper-type cross-slip event can be performed
 *                      without conflicting with a previously performed
 *                      zipper-type event.  Before returning to the caller
 *                      the status of the current cross-slip event is
 *                      updated to reflect whether the event has been
 *                      rejected or approved.
 *
 *                      Since all zipper-type cross-slips are performed
 *                      before the non-zipper cross-slips, we do not
 *                      check for conflicts with the non-zipper events here.
 *
 *                      IMPORTANT: This function does not use any locking
 *                                 mechanism to access the eventList and
 *                                 as such is not thread-safe!  This is
 *                                 not a problem so long as the function
 *                                 is not called from a threaded portion of
 *                                 code.
 *
 *                      NOTE: The zipper-type cross-slip events are
 *                            contained in the first <zipEventCnt>
 *                            structures of <eventList>, while the
 *                            non-zipper cross-slip events are contained
 *                            in the last <eventListSize> - <slipEventIndex>
 *                            structures of <eventList>
 *
 *
 *      Arguments:
 *          eventList      Pointer to the cross-slip event list
 *          zipEventCnt    Number of zipper-type cross-slip events on
 *                         the <eventList>
 *          currEventIndex Index in <eventList> of the zipper-type
 *                         cross-slip event we wish to perform.
 *
 *      Returns: 0 if the cross-slip event conflicts with a prior event
 *               1 if the corss-slip event may proceed.
 *
 *-------------------------------------------------------------------------*/
int ProceedWithZipEvent(SlipEvent_t *eventList, int zipEventCnt,
                         int currEventIndex)
{
        int         i;
        real8       cp[3];
        Node_t      *currNode, *currNbr, *prevNode, *prevNbr;
        SlipEvent_t *currEvent, *prevEvent;

        currEvent = &eventList[currEventIndex];
        currNode = currEvent->node;
        currNbr  = currEvent->nbrNodes[currEvent->segIndex];

        for (i = 0; i < currEventIndex; i++) {

            if (eventList[i].status != CROSS_SLIP_COMMITTED) {
                continue;
            }

            prevEvent = &eventList[i];
            prevNode = prevEvent->node;
            prevNbr  = prevEvent->nbrNodes[prevEvent->segIndex];

/*
 *          If the current cross-slip event affects a node that was
 *          involved in a previous cross-slip event, we will only allow
 *          permit the current event to proceeed if the glide planes
 *          for the two events match.
 */
            if ((currNode == prevNbr)  ||
                (currNbr  == prevNode) ||
                (currNbr  == prevNbr)) {

                cross(prevEvent->newPlane, currEvent->newPlane, cp);

                if (DotProduct(cp, cp) > 1.0e-3) {
                    currEvent->status = CROSS_SLIP_REJECTED;
                    return(0);
                }
            }
        }
        currEvent->status = CROSS_SLIP_COMMITTED;

        return(1);
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ProceedWithCrossSlip
 *
 *      Description:    This function will determine if the current
 *                      non-zipper cross-slip event can be performed
 *                      without conflicting with a previously performed
 *                      cross-slip event whether that be a zipper or
 *                      non-zipper event.  Before returning to the caller
 *                      the status of the current cross-slip event is
 *                      updated to reflect whether the event has been
 *                      rejected or approved.
 *
 *                      IMPORTANT: This function does not use any locking
 *                                 mechanism to access the eventList and
 *                                 as such is not thread-safe!  This is
 *                                 not a problem so long as the function
 *                                 is not called from a threaded portion of
 *                                 code.
 *
 *                      NOTE: The zipper-type cross-slip events are
 *                            contained in the first <zipEventCnt>
 *                            structures of <eventList>, while the
 *                            non-zipper cross-slip events are contained
 *                            in the last <eventListSize> - <slipEventIndex>
 *                            structures of <eventList>
 *
 *
 *      Arguments:
 *          eventList      Pointer to the cross-slip event list
 *          zipEventCnt    Number of zipper-type cross-slip events on
 *                         the <eventList>
 *          firstSlipEvent Index in <eventList> of the first non-zipper
 *                         type cross-slip event
 *          eventListSize  Size of the <eventList> array as a count of
 *                         EventList_t structures.
 *          currEventIndex Index in <eventList> of the zipper-type
 *                         cross-slip event we wish to perform.
 *
 *      Returns: 0 if the cross-slip event conflicts with a prior event
 *               1 if the corss-slip event may proceed.
 *
 *-------------------------------------------------------------------------*/
int ProceedWithCrossSlip(SlipEvent_t *eventList, int zipEventCnt,
                         int firstSlipIndex, int eventListSize,
                         int currEventIndex)
{
        int         i;
        real8       cp[3];
        Node_t      *currNode, *currNbr1, *currNbr2;
        Node_t      *prevNode, *prevNbr1, *prevNbr2, *prevNbr;
        SlipEvent_t *currEvent, *prevEvent;

/*
 *      If the current cross-slip event affects a node that was
 *      involved in a previous cross-slip event, we will only allow
 *      permit the current event to proceeed if the glide planes
 *      for the two events match.
 *
 *      First check that this normal cross-slip event does not conflict
 *      with any previous 'zipper' cross-slip events.
 */
        currEvent = &eventList[currEventIndex];
        currNode = currEvent->node;
        currNbr1 = currEvent->nbrNodes[0];
        currNbr2 = currEvent->nbrNodes[1];

        for (i = 0; i < zipEventCnt; i++) {

            if (eventList[i].status != CROSS_SLIP_COMMITTED) {
                continue;
            }

            prevEvent = &eventList[i];
            prevNode = prevEvent->node;
            prevNbr  = prevEvent->nbrNodes[prevEvent->segIndex];

            if ((currNode == prevNbr)  ||
                (currNbr1 == prevNode) ||
                (currNbr1 == prevNbr)  ||
                (currNbr2 == prevNode) ||
                (currNbr2 == prevNbr)) {

                cross(prevEvent->newPlane, currEvent->newPlane, cp);

                if (DotProduct(cp, cp) > 1.0e-3) {
                    currEvent->status = CROSS_SLIP_REJECTED;
                    return(0);
                }
            }
        }

/*
 *      Okay, no conflict with any of the 'zipper' events, so now check
 *      for any conflicts with previous non-zipper cross slip events.
 */
        for (i = firstSlipIndex; i < currEventIndex; i++) {

            if (eventList[i].status != CROSS_SLIP_COMMITTED) {
                continue;
            }

            prevEvent = &eventList[i];
            prevNode = prevEvent->node;
            prevNbr1 = prevEvent->nbrNodes[0];
            prevNbr2 = prevEvent->nbrNodes[1];

            if ((currNode == prevNbr1)  ||
                (currNode == prevNbr2)  ||
                (currNbr1 == prevNbr1)  ||
                (currNbr1 == prevNode)  ||
                (currNbr1 == prevNbr2)  ||
                (currNbr2 == prevNbr1)  ||
                (currNbr2 == prevNode)  ||
                (currNbr2 == prevNbr2)) {

                cross(prevEvent->newPlane, currEvent->newPlane, cp);

                if (DotProduct(cp, cp) > 1.0e-3) {
                    currEvent->status = CROSS_SLIP_REJECTED;
                    return(0);
                }
            }
        }

        currEvent->status = CROSS_SLIP_COMMITTED;

        return(1);
}


/*---------------------------------------------------------------------------
 *
 *      Function:       SaveCrossSlipInfo()
 *
 *      Description:    This function saves some of the original nodal
 *                      force/position information from a set of nodes
 *                      prior to attempting cross-slip events.  If after
 *                      adjusting the nodes involved, it turns out the
 *                      cross-slip event should not occur, this information
 *                      will be used to restored the configuration to
 *                      its original state.
 *
 *      Arguments:
 *          node         pointer to primary two-node
 *          nbr1         pointer to first neighboring node
 *          nbr2         pointer to second neighboring node
 *          nbr1ArmID    index for segment nbr1/node in the <nbr1> node struct
 *          nbr2ArmID    index for segment nbr2/node in the <nbr2> node struct
 *          segForceOrig array in whcih to store the original force on the two
 *                       segments.
 *
 *                       segForceOrig[0][*] force at <node> from node/nbr1 seg
 *                       segForceOrig[1][*] force at <nbr1> from node/nbr1 seg
 *                       segForceOrig[2][*] force at <node> from node/nbr2 seg
 *                       segForceOrig[3][*] force at <nbr2> from node/nbr2 seg
 *
 *          nodePosOrig  array in which to store coordinates for <node>
 *          nbr1PosOrig  array in which to store coordinates for <nbr1>
 *          nbr2PosOrig  array in which to store coordinates for <nbr2>
 *
 *-------------------------------------------------------------------------*/
void SaveCrossSlipInfo(Node_t *node, Node_t *nbr1, Node_t *nbr2,
                       int nbr1ArmID, int nbr2ArmID, real8 segForceOrig[4][3],
                       real8 nodePosOrig[3], real8 nbr1PosOrig[3],
                       real8 nbr2PosOrig[3])
{
        
       /* original force at both endpoints of the first segment */
       segForceOrig[0][X] = node->armfx[0];
       segForceOrig[0][Y] = node->armfy[0];
       segForceOrig[0][Z] = node->armfz[0];

       segForceOrig[1][X] = nbr1->armfx[nbr1ArmID];
       segForceOrig[1][Y] = nbr1->armfy[nbr1ArmID];
       segForceOrig[1][Z] = nbr1->armfz[nbr1ArmID];

       /* original force at both endpoints of the second segment */
       segForceOrig[2][X] = node->armfx[1];
       segForceOrig[2][Y] = node->armfy[1];
       segForceOrig[2][Z] = node->armfz[1];

       segForceOrig[3][X] = nbr2->armfx[nbr2ArmID];
       segForceOrig[3][Y] = nbr2->armfy[nbr2ArmID];
       segForceOrig[3][Z] = nbr2->armfz[nbr2ArmID];

       /* original position of all three nodes */
       nodePosOrig[X] = node->x;
       nodePosOrig[Y] = node->y;
       nodePosOrig[Z] = node->z;

       nbr1PosOrig[X] = nbr1->x;
       nbr1PosOrig[Y] = nbr1->y;
       nbr1PosOrig[Z] = nbr1->z;

       nbr2PosOrig[X] = nbr2->x;
       nbr2PosOrig[Y] = nbr2->y;
       nbr2PosOrig[Z] = nbr2->z;

       return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ResetPosition()
 *
 *      Description:    Reset the position of a node.  If periodic boundaries
 *                      are enabled and the new position is outside the
 *                      problem space, translate the coordinates to the 
 *                      appropriate location within the problem space.
 *
 *      Arguments:
 *          node  Pointer to node to be updated
 *          pos   new nodal coordinates
 *
 *-------------------------------------------------------------------------*/
void ResetPosition(Param_t *param, Node_t *node, real8 pos[3])
{
        node->x = pos[X];
        node->y = pos[Y];
        node->z = pos[Z];

        FoldBox(param, &node->x, &node->y, &node->z);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       RestoreCrossSlipForce()
 *
 *      Description:    Given a two-node and its neighbors, restore the
 *                      forces on the node's segments (at both endpoints
 *                      of each segment) from the preserved force values.
 *
 *      Arguments:
 *          node         pointer to primary two-node
 *          nbr1         pointer to first neighboring node
 *          nbr2         pointer to second neighboring node
 *          nbr1ArmID    index for segment nbr1/node in the <nbr1> node struct
 *          nbr2ArmID    index for segment nbr2/node in the <nbr2> node struct
 *          segForceOrig original force on the two segments.
 *                       segForceOrig[0][*] force at <node> from node/nbr1 seg
 *                       segForceOrig[1][*] force at <nbr1> from node/nbr1 seg
 *                       segForceOrig[2][*] force at <node> from node/nbr2 seg
 *                       segForceOrig[3][*] force at <nbr2> from node/nbr2 seg
 *
 *-------------------------------------------------------------------------*/
void RestoreCrossSlipForce(Node_t *node, Node_t *nbr1, Node_t *nbr2,
                           int nbr1ArmID, int nbr2ArmID,
                           real8 segForceOrig[4][3])
{
/*
 *      Reset force at <node> for both segments of and simply set the
 *      total nodal force to the sum of both.
 */
        node->armfx[0] = segForceOrig[0][X];
        node->armfy[0] = segForceOrig[0][Y];
        node->armfz[0] = segForceOrig[0][Z];

        node->armfx[1] = segForceOrig[2][X];
        node->armfy[1] = segForceOrig[2][Y];
        node->armfz[1] = segForceOrig[2][Z];

        node->fX = segForceOrig[0][X] + segForceOrig[2][X];
        node->fY = segForceOrig[0][Y] + segForceOrig[2][Y];
        node->fZ = segForceOrig[0][Z] + segForceOrig[2][Z];

/*
 *      The neighbor nodes are having forces reset for only a single segment
 *      each, so we have to subtract off the current portion of the total
 *      force for the specific segments, update the segment force, and add
 *      the new segment force into the total force for each node.
 */
        nbr1->fX -= nbr1->armfx[nbr1ArmID];
        nbr1->fY -= nbr1->armfy[nbr1ArmID];
        nbr1->fZ -= nbr1->armfz[nbr1ArmID];

        nbr2->fX -= nbr2->armfx[nbr2ArmID];
        nbr2->fY -= nbr2->armfy[nbr2ArmID];
        nbr2->fZ -= nbr2->armfz[nbr2ArmID];

        nbr1->armfx[nbr1ArmID] = segForceOrig[1][X];
        nbr1->armfy[nbr1ArmID] = segForceOrig[1][Y];
        nbr1->armfz[nbr1ArmID] = segForceOrig[1][Z];

        nbr2->armfx[nbr2ArmID] = segForceOrig[3][X];
        nbr2->armfy[nbr2ArmID] = segForceOrig[3][Y];
        nbr2->armfz[nbr2ArmID] = segForceOrig[3][Z];

        nbr1->fX += nbr1->armfx[nbr1ArmID];
        nbr1->fY += nbr1->armfy[nbr1ArmID];
        nbr1->fZ += nbr1->armfz[nbr1ArmID];

        nbr2->fX += nbr2->armfx[nbr2ArmID];
        nbr2->fY += nbr2->armfy[nbr2ArmID];
        nbr2->fZ += nbr2->armfz[nbr2ArmID];

        return;
}


void CrossSlip(Home_t *home)
{
        int enabled = home->param->enableCrossSlip;
        int matType = home->param->materialType;
        int bccIndx = home->param->crossSlipBCCIndex;

        if (enabled) 
        {
            if (matType==MAT_TYPE_BCC)
            {
               if (bccIndx==0) { CrossSlipBCC (home); }
               if (bccIndx==1) { CrossSlipBCCRandom(home); }

            }

            else if (matType==MAT_TYPE_FCC             ) { CrossSlipFCC(home); }
            else if (matType==MAT_TYPE_HCP             ) { CrossSlipHCP(home); }
            else if (matType==MAT_TYPE_RHOMBOHEDRAL_VA ) { CrossSlipRhombohedralVa(home); }
        }

        return;
}
