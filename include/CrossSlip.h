//--------------------------------------------------------------------------
// CrossSlip.h  Declare function prototypes and other relevant
//              cross-slip related data needed by modules performing
//              the cross-slip events.
//--------------------------------------------------------------------------

#ifndef PDS_CROSS_SLIP_H
#define PDS_CROSS_SLIP_H

#include "Home.h"

// Define a few flags used when storing information about
// potential cross-slip events

#define SLIP_REPOSITION_NODE 0x0001
#define SLIP_REPOSITION_NBR  0x0002
#define SLIP_ZIPPER_EVENT    0x0004

// Define a set of status values that are used to indicate
// the state of a potential cross-slip event that has been
// stored in the corss-slip event list

enum {
        CROSS_SLIP_PENDING     = 0,
        CROSS_SLIP_COMMITTED      ,
        CROSS_SLIP_BACKED_OFF     ,
        CROSS_SLIP_REJECTED
};

// Define a structure used by the cross-slip functions to store
// all the information needed to perform a cross-slip event.
//
// Note: not all cross-slip events require all data items.

typedef struct {
        int    flags;         // Contains various bit-flags pertaining
                              // to the specific cross-slip event     

        int    status;        // Status of the event

        int    segIndex;      // Index of the segment of <node> that is
                              // associated with this cross-slip event 

        real8  fdotglide;     // Original fdotglide value for standard
                              // cross-slip events, to be compared with
                              // a new fdotglide after new force calcs
                              // are done during the cross-slip attempt

        real8  newNodePos[3]; // Position to which <node> will be moved
                              // during the corss-slip event           

        real8  newNbrPos[3];  // Position to which <nbrNode[<segIndex>]>
                              // will be moved during the cross-slip (if
                              // the neighbor is to be repositioned     

        real8  newPlane[3];   // Components of the new glide plane

        real8  glideDir[3];   // Components of the original glide plane

        Node_t *node;         // Primary cross-slip node

        Node_t *nbrNodes[2];  // One or more neighboring of <node> that
                              // are affected by the cross-slip event  
} SlipEvent_t;

// Prototype the cross-slip related functions

extern void CrossSlip               (Home_t *home);
extern void CrossSlipBCC            (Home_t *home);
extern void CrossSlipBCCAllPlanes   (Home_t *home);
extern void CrossSlipBCCRandom      (Home_t *home);
extern void CrossSlipFCC            (Home_t *home);
extern void CrossSlipHCP            (Home_t *home);
extern void CrossSlipRhombohedralVa (Home_t *home);

#ifdef DEBUG_CROSSSLIP_EVENTS
extern void DumpCrossSlipEvent      (Node_t *node, real8 newPlane[3], char *eventDesc);
#endif

extern void PrintCrossSlipEventList (SlipEvent_t *eventList, int zipEventCnt, int slipEventIndex, int eventListSize);
extern int  ProceedWithZipEvent     (SlipEvent_t *eventList, int zipEventCnt, int currEventIndex);
extern int  ProceedWithCrossSlip    (SlipEvent_t *eventList, int zipEventCnt, int firstSlipIndex, int eventListSize, int currEventIndex);

extern void SaveCrossSlipInfo       (Node_t *node, Node_t *nbr1, Node_t *nbr2, int nbr1ArmID, int nbr2ArmID, 
                                     real8 segForceOrig[4][3], real8 nodePosOrig[3], real8 nbr1PosOrig[3], real8 nbr2PosOrig[3]);

extern void ResetPosition           (Param_t *param, Node_t *node, real8 pos[3]);

extern void RestoreCrossSlipForce   (Node_t *node, Node_t *nbr1, Node_t *nbr2, int nbr1ArmID, int nbr2ArmID, real8 segForceOrig[4][3]);

#endif   // PDS_CROSS_SLIP_H
