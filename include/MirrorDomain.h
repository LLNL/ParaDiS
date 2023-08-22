#pragma once

#ifndef _PDS_MIRROR_DOMAIN_H
#define _PDS_MIRROR_DOMAIN_H

//---------------------------------------------------------------------------
// MirrorDomain.h: Define the struct that holds all nodal position data
//                 from a remote (mirror) domain.  This data is used
//                 for X-window plotting and other output. Only domain 0
//                 populates and uses this structure.
//---------------------------------------------------------------------------

#include "Typedefs.h"
#include "Home.h"

typedef struct MirrorDomain_t
{
        Node_t **nodeKeys;              // indexed by node's tag.index, points
                                        // to Node_t
        int newNodeKeyPtr;              // elements of nodeKeys array

        real8 *armX;                    // arm[XYZ] are pointers to the corresponding
        real8 *armY;                    // X,Y and Z coordinates of nodes at the far
        real8 *armZ;                    // end of arms.  These arrays are only used
                                        // during the process of uploading the problem
                                        // to task 0 when generating output, and are
                                        // needed because some of the output functions
                                        // require the neighbors coordinates, but the
                                        // entire problem space may not fit into
                                        // task 0's memory.  So, the remote domains
                                        // send a node's neighbor node's coords along
                                        // with the arms data.

#ifdef ESHELBY
        int           inclusionCount;   // number of inclusions in remote domain
        EInclusion_t *inclusionList ;   // Array of inclusions:  NOTE: for these inclusions,
                                        // the only data provided is id, radius and the position!
#endif

} MirrorDomain_t;

// Forward declarations...

extern void MirrorDomain_Free (Home_t *home, int domIndex);

extern void MirrorDomain_Print(FILE *strm, MirrorDomain_t *mdom);

#endif
