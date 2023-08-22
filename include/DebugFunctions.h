#pragma once

#ifndef _PDS_DEBUG_FUNCTIONS_H
#define _PDS_DEBUG_FUNCTIONS_H

/****************************************************************************
 *
 *      DebugFunctions.h  Contains prototypes for various debug-only functions
 *
 ***************************************************************************/

extern void   CheckBurgConservation     (Home_t *home, const char   *msg);
extern int    CheckNodeBurgConservation (Home_t *home, Node_t *node);
extern void   CheckForces               (Home_t *home, const char   *msg);
extern void   CheckSegLengths           (Home_t *home, const char   *msg);
extern void   CheckForNANS              (Home_t *home, const char   *msg);
extern void   CheckForEmptySimulation   (Home_t *home);
extern void   CheckForUndefinedPlanes   (Home_t *home, const char   *msg);
extern void   CheckPlanes_RhomboVa      (Home_t *home, const char   *msg);

extern void   CheckSegmentLengths       (const Home_t *home, const char *msg);

extern void   Print2x2                  (const char *msg, real8 A[2][2]);
extern void   Synchronize               (Home_t *home, const char   *msg);

extern void   VerifyConnections         (Home_t *home, const char   *msg);
extern void   VerifyLocalConnections    (Home_t *home, const char   *msg);

#ifdef CHECK_LARGE_SEGS
#define CHECK_SEGMENT_LENGTHS(a) \
{ \
   char msg[256];  sprintf(msg, "%s::%s(ln=%d)", __FILE__, __func__, __LINE__ ); \
   CheckSegmentLengths(a,msg); \
}
#else
#define CHECK_SEGMENT_LENGTHS(a)
#endif


#endif
