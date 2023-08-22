#pragma once

#ifndef _PDS_MOBILITY_H
#define _PDS_MOBILITY_H

/***************************************************************************
 *
 *      Mobility.h      Declare function prototypes and other relevant
 *                      mobility related data needed by modules invoking
 *                      the mobility laws.
 *
 **************************************************************************/

#include "Home.h"
#include "MatType_t.h"

/*
 *      For the 'faceted' mobility functions, the user can specify
 *      a varying number of center directions (trenches) along which
 *      there is increased drag on dislocations.  This defines the
 *      maximum number of trenches permitted in the statically defined
 *      arrays of trench info
 */
#define MAX_MOB_TRENCHES 10

/*
 *      Define a set of integer values corresponding to the
 *      various types of available mobility laws.  During
 *      initialization the mobilityType parameter will be
 *      set to one of these values which can be used throughout
 *      the code to determine the current mobility in use.
 *      this is easier than using strcmp() everywhere we want
 *      to know make decisions based on the mobility law.
 *
 *      IMPORTANT:  When a new mobility type is added to the enumerated
 *                  list below, you must add an entry with the mobility
 *                  'attributes' to the <mobAttrList> array contained
 *                  in this file.
 */
typedef enum {
        MOB_BCC_0,
#ifdef ESHELBY
        MOB_BCC_0_ESHELBY,
#endif
        MOB_BCC_0B,
#ifdef ESHELBY
        MOB_BCC_0B_ESHELBY,
#endif
        MOB_BCC_GLIDE,
        MOB_BCC_LINEAR,
#ifdef ESHELBY
        MOB_BCC_LINEAR_ESHELBY,
#endif
        MOB_BCC_FACETED,
        MOB_BCC_NL,
        MOB_BCC_TA_NL,
        MOB_BCC_TA_NL_B,
        MOB_BCC_TA_NL_B_PLANAR,
        MOB_BCC_VA_NL,
        MOB_BCC_VA_NL_PLANAR,
        MOB_BCC_FE_NL,
        MOB_BCC_FE_NL_A,
        MOB_FCC_0,
        MOB_FCC_0B,
        MOB_FCC_CLIMB,
        MOB_FCC_LINEAR,
#ifdef ESHELBY
        MOB_FCC_0B_ESHELBY,
        MOB_FCC_LINEAR_ESHELBY,
#endif
        MOB_FCC_ANGLE,
        MOB_HCP_LINEAR,
#ifdef ESHELBY
        MOB_HCP_LINEAR_ESHELBY,
#endif
        MOB_TA,
        MOB_TA_LINEAR,
        MOB_TA_PENCIL,
        MOB_RHOMBOHEDRAL_VA,
        MOB_RELAX,
        MOB_RELAX_GLIDE,
        MOB_MAX_INDEX  /* MUST be the last item in the list */
} MobType_t;


/*
 *      Prototype the available mobility functions.  All mobility
 *      functions must have the same number and type of parameters
 *      and must return an integer error flag.  Flag is 1 if the
 *      nodal velocity could not be calculated, 0 if it could.
 */
int   Mobility_BCC_0                     (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_0b                    (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_faceted               (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_glide                 (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_linear                (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_nl                    (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_Ta_nl                 (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_Ta_nl_b               (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_Ta_nl_b_planar        (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_Va_nl                 (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_Va_nl_planar          (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_Fe_nl                 (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_Fe_nl_a               (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_FCC_0                     (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_FCC_0b                    (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_FCC_climb                 (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_FCC_linear                (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_FCC_angle                 (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_HCP_linear                (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_Ta                        (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_Ta_linear                 (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_Ta_pencil                 (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_Rhombohedral_Va_nl_planar (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_Relax                     (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_Relax_glide               (Home_t *home, Node_t *node, MobArgs_t *mobArgs);

/* Functions for Tantalum non-linear screw mobility law */
void  TAEdgeDrag     (real8 vel[3], real8 burg[3], real8 lineDir[3],
                      real8 dragClimb, real8 dragLine, real8 dragGlide,
                      real8 fEdgeDrag[3], real8 dfEdgeDragdv[3][3]);

void  TAJunctionDrag (real8 vel[3], real8 lineDir[3],
                      real8 BJunctionLine, real8 BJunctionGlide,
                      real8 fJunctionDrag[3], real8 dfJunctionDragdv[3][3]);

/*
 *      Some of the mobility functions use values that are dependent
 *      on the control parameters.  For those mobilities we provide
 *      an 'init' function to set those variables once during initialization
 */
void  InitMob_BCC_faceted                (Home_t *home);
void  InitMob_BCC_nl                     (Home_t *home);
void  InitMob_BCC_Ta_nl                  (Home_t *home);
void  InitMob_BCC_Ta_nl_b                (Home_t *home);
void  InitMob_BCC_Ta_nl_b_planar         (Home_t *home);
void  InitMob_BCC_Va_nl                  (Home_t *home);
void  InitMob_BCC_Va_nl_planar           (Home_t *home);
void  InitMob_BCC_Fe_nl                  (Home_t *home);
void  InitMob_BCC_Fe_nl_a                (Home_t *home);
void  InitMob_Rhombohedral_Va_nl_planar  (Home_t *home);

/*
 *  Functions used in several linear mobility laws
 */
void EdgeDrag     (real8 vel[3], real8 burg[3], real8 edgeDir[3],
                   real8 dragClimb, real8 dragLine, real8 dragGlide,
                   real8 fEdgeDrag[3], real8 dfEdgeDragdv[3][3]);

void ScrewDrag    (real8 vel[3], real8 burg[3],
                   real8 dragClimb, real8 dragLine, real8 dragGlide,
                   real8 fScrewDrag[3], real8 dfScrewDragdv[3][3]);

void ScrewDrag2   (real8 vel[3], real8 burg[3],real8 climbDir[3],
                   real8 dragClimb, real8 dragLine, real8 dragGlide,
                   real8 fScrewDrag[3], real8 dfScrewDragdv[3][3]);

void JunctionDrag (real8 vel[3], real8 lineDir[3],
                   real8 BJunctionLine, real8 BJunctionGlide,
                   real8 fJunctionDrag[3], real8 dfJunctionDragdv[3][3]);

void AngleDrag    (real8 vel[3],
                   real8 lineDir[3],
                   real8 climbDir[3],
                   real8 dragClimb,
                   real8 dragLine,
                   real8 dragGlide_B,
                   real8 dragGlide_D,
                   real8 dragGlide_v0,
                   real8 fAngleDrag[3],
                   real8 dfAngleDragdv[3][3]);

void NonlinearGlideDrag(real8 Vg, real8 *V, real8 B, real8 D, real8 V0, real8 *g, real8 *F, real8 dFdv[3][3]);

/*
 *      If we have material specific functions for dynamically calculating
 *      the available burgers vectors (and associated glide planes) prototype
 *      those functions here.
 */
void  GetBurgList(Home_t *home);

/*
 * For BCC
 */
void  GetBurgList_BCC     (real8 (**burgList)[3], real8 (**planeList)[3],
                           int **numPlanesPerBurg, int **burgFirstPlaneIndex,
                           int **numGlissilePlanes, int **planeType,
                           int *numBurg, int *numPlanes);

int  GetBurgIndexBCC     (real8 burgVec[3], int numBurg, real8 (*burgList)[3]);

void  FindGlidePlaneBCC   (Home_t *home, int bIndex, real8 plane[3], int *pIndex);


void GetBCCAllIndices(Home_t *home, real8 b[3], real8 p[3],
                      int *bIndex, int *pIndex, int *gIndex);

/*
 * For HCP
 */
void  GetBurgIndexHCP     (real8 burgVec[3], int useClosest, int numBurg, real8 (*burgList)[3], int *burgIndex);
void  GetBurgList_HCP     (real8 cOa, real8 (**burgList)[3], real8 (**planeList)[3],
                           SplinterableBurg_t **splinterableBurgList,
                           int **numPlanesPerBurg, int **burgFirstPlaneIndex,
                           int **numGlissilePlanes, int **planeType,
                           int *numBurg, int *numPlanes, int *numSplinterableBurgs);

void  FindGlidePlaneHCP   (Home_t *home, int bIndex, real8 plane[3], int *pIndex);
real8 AssignWeightToPlane (real8 Fac[5], int indexPlane);

/*
 * For FCC
 */
void  GetBurgList_FCC     (real8 (**burgList)[3], real8 (**planeList)[3],
                           int **numPlanesPerBurg, int **burgFirstPlaneIndex,
                           int **numGlissilePlanes, int **planeType,
                           int *numBurg, int *numPlanes);

void  GetBurgIndexFCC     (real8 burgVec[3], int numBurg, real8 (*burgList)[3], int *burgIndex);

void  FindGlidePlaneFCC   (Home_t *home, int bIndex, real8 plane[3], int *pIndex);

#ifdef ESHELBY
int   Mobility_BCC_0_Eshelby      (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_0b_Eshelby     (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_FCC_0b_Eshelby     (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_BCC_linear_Eshelby (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_FCC_linear_Eshelby (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
int   Mobility_HCP_linear_Eshelby (Home_t *home, Node_t *node, MobArgs_t *mobArgs);
#endif

/*
 *      Need to prototype some functions specific to the original non-linear
 *      mobility module.
 */
void GetNonLinearMobConstants(Param_t *param);

/*
 *      Define a structure containing various attributes of an individual
 *      mobility module.  By storing these mobility-specific data in
 *      this structure, we can do a simple lookup of the specified
 *      'mobilityLaw' and get all the data items at once.
 */
typedef struct {
        char mobName[64];  /* String used to select the mobility module */
                           /* via the <mobilityLaw> control parameter   */

        int  defaultMatType;   /* Default type of material structure for */
                               /* associated mobility module.  Must be   */
                               /* one of the enumerated types in         */
                               /* MatType_t above.                       */


        int  mobilityType; /* Integer value corresponding to the specific */
                           /* mobility law.  Redundant info, but easier to*/
                           /* use in the code than the string name        */

#define NO_WARN_ON_MAT_TYPE_MISMATCH 0
#define WARN_ON_MAT_TYPE_MISMATCH    1
        int  warnOnMatTypeMismatch; /* Set to 1 if a warning should be     */
                                    /* be issued if the user-specified     */
                                    /* material type does not match the    */
                                    /* default material type for this      */
                                    /* mobility module. Set to 0 otherwise */
#define NON_PLANAR 0
#define PLANAR     1
        int  isPlanarMobility; /* Set to 1 if the mobility module is a  */
                               /* planar mobility which requires that   */
                               /* the <enforceGlidePlanes> control      */
                               /* parameter be set.  Set to 0 otherwise */
#define NO_FUZZY_PLANES  0
#define FUZZY_PLANES     1
        int  allowFuzzyPlanes; /* Set to 1 if the mobility module is    */
                               /* planar but allows some climb yiedling */
                               /* some 'fuzziness' in the glide planes. */
                               /* Set to 0 in all other cases           */
#define NO_CROSS_SLIP  0
#define CROSS_SLIP     1
        int  doCrossSlip;  /* Set to 1 if cross-slip should be enabled */
                           /* by default for this mobility.  Set to 0  */
                           /* otherwise.                               */

        void (*mobInitFunc)(Home_t *home); /* Set to the initialization    */
                                           /* function (if applicable) for */
                                           /* the mobility module.  Set to */
                                           /* NULL if not applicable       */

        /* Indicates the function to be invoked for calculating mobilities */
        /* for this mobility module                                        */
        int  (*mobFunc)(Home_t *home, Node_t *node, MobArgs_t *mobArgs);
} MobAttr_t;

/*
 *      Since we now have mobility functions that are not material-specific
 *      the user can dynamically select the material type, which also
 *      changes the number of burgers vector groupings used for certain
 *      portions of the code.  So, rather than tie the number of burg groups
 *      to the mobility function, we need an array that ties the burg groups
 *      to the material type.
 */
typedef struct {
        char materialTypeName[64];
        int  materialType;
        int  numBurgGroups;
} MatInfo_t;

/*
 *      Now define/initialize the array of mobility function attributes,
 *      and the array mapping material type to number of burgers vector
 *      groupings.  During initialization we'll look up the mobility law
 *      specified for the simulation in this array and use it to set up
 *      some mobility-specific parameters.
 *
 *      IMPORTANT: Every mobility module defined above in the MobType_t
 *      enumerated list should have a corresponding entry in this array!
 *
 *      NOTE: Only 1 source module should define INIT_MOB_ATTR_LIST.
 *      This way the array is declared and initialized in only 1 module
 *      and referenced as an external in all others.
 */
#ifndef INIT_MOB_ATTR_LIST
extern const MatInfo_t matInfo[MAT_TYPE_MAX_INDEX];
extern const MobAttr_t mobAttrList[MOB_MAX_INDEX];
#else
const MatInfo_t matInfo[MAT_TYPE_MAX_INDEX] =
        {
            {MAT_TYPE_NAME_BCC,       MAT_TYPE_BCC,              5},
            {MAT_TYPE_NAME_FCC,       MAT_TYPE_FCC,              7},
            {MAT_TYPE_NAME_HCP,       MAT_TYPE_HCP,             11},
            {MAT_TYPE_NAME_RHOMBO_VA, MAT_TYPE_RHOMBOHEDRAL_VA,  0}
        };

const MobAttr_t mobAttrList[MOB_MAX_INDEX] =
        {
            {"BCC_0",                     // mobilityLaw name
              MAT_TYPE_BCC,               // default material type
              MOB_BCC_0,                  // mobility type
              WARN_ON_MAT_TYPE_MISMATCH,  // issue warning if user specified matType
                                          //    doesn't match default
              NON_PLANAR,                 // is a planar mobility?
              NO_FUZZY_PLANES,            // allows fuzzy planes?
              NO_CROSS_SLIP,              // use cross-slip as default?
              NULL,                       // mobility init function (if any)
              Mobility_BCC_0              // mobility function
            },

#ifdef ESHELBY
            {"BCC_0_Eshelby",
              MAT_TYPE_BCC,
              MOB_BCC_0_ESHELBY,
              WARN_ON_MAT_TYPE_MISMATCH,
              NON_PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              NULL,
              Mobility_BCC_0_Eshelby
            },
#endif

            {"BCC_0b",
              MAT_TYPE_BCC,
              MOB_BCC_0B,
              WARN_ON_MAT_TYPE_MISMATCH,
              NON_PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              NULL,
              Mobility_BCC_0b
            },

#ifdef ESHELBY
            {"BCC_0b_Eshelby",
              MAT_TYPE_BCC,
              MOB_BCC_0B_ESHELBY,
              WARN_ON_MAT_TYPE_MISMATCH,
              NON_PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              NULL,
              Mobility_BCC_0b_Eshelby
            },
#endif

            {"BCC_glide",
              MAT_TYPE_BCC,
              MOB_BCC_GLIDE,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_BCC_glide
            },

            {"BCC_linear",
              MAT_TYPE_BCC,
              MOB_BCC_LINEAR,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_BCC_linear
            },

#ifdef ESHELBY
             {"BCC_linear_Eshelby",
              MAT_TYPE_BCC,
              MOB_BCC_LINEAR_ESHELBY,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_BCC_linear_Eshelby
            },
#endif

           {"BCC_faceted",
              MAT_TYPE_BCC,
              MOB_BCC_FACETED,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              InitMob_BCC_faceted,
              Mobility_BCC_faceted
            },

            {"BCC_nl",
              MAT_TYPE_BCC,
              MOB_BCC_NL,
              WARN_ON_MAT_TYPE_MISMATCH,
              NON_PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              InitMob_BCC_nl,
              Mobility_BCC_nl
            },

            {"BCC_Ta_nl",
              MAT_TYPE_BCC,
              MOB_BCC_TA_NL,
              WARN_ON_MAT_TYPE_MISMATCH,
              NON_PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              InitMob_BCC_Ta_nl,
              Mobility_BCC_Ta_nl
            },

            {"BCC_Ta_nl_b",
              MAT_TYPE_BCC,
              MOB_BCC_TA_NL_B,
              WARN_ON_MAT_TYPE_MISMATCH,
              NON_PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              InitMob_BCC_Ta_nl_b,
              Mobility_BCC_Ta_nl_b
            },

            {"BCC_Ta_nl_b_planar",
              MAT_TYPE_BCC,
              MOB_BCC_TA_NL_B_PLANAR,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              InitMob_BCC_Ta_nl_b_planar,
              Mobility_BCC_Ta_nl_b_planar
            },

            {"BCC_Va_nl",
              MAT_TYPE_BCC,
              MOB_BCC_VA_NL,
              WARN_ON_MAT_TYPE_MISMATCH,
              NON_PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              InitMob_BCC_Va_nl,
              Mobility_BCC_Va_nl
            },

            {"BCC_Va_nl_planar",
              MAT_TYPE_BCC,
              MOB_BCC_VA_NL_PLANAR,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              InitMob_BCC_Va_nl_planar,
              Mobility_BCC_Va_nl_planar
            },

            {"BCC_Fe_nl",
              MAT_TYPE_BCC,
              MOB_BCC_FE_NL,
              WARN_ON_MAT_TYPE_MISMATCH,
              NON_PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              InitMob_BCC_Fe_nl,
              Mobility_BCC_Fe_nl
            },

            {"BCC_Fe_nl_a",
              MAT_TYPE_BCC,
              MOB_BCC_FE_NL_A,
              WARN_ON_MAT_TYPE_MISMATCH,
              NON_PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              InitMob_BCC_Fe_nl_a,
              Mobility_BCC_Fe_nl_a
            },

            {"FCC_0",
              MAT_TYPE_FCC,
              MOB_FCC_0,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_FCC_0
            },

            {"FCC_0b",
              MAT_TYPE_FCC,
              MOB_FCC_0B,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_FCC_0b
            },

            {"FCC_climb",
              MAT_TYPE_FCC,
              MOB_FCC_CLIMB,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_FCC_climb
            },

            {"FCC_linear",
              MAT_TYPE_FCC,
              MOB_FCC_LINEAR,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_FCC_linear
            },

#ifdef ESHELBY
            {"FCC_0b_Eshelby",
              MAT_TYPE_FCC,
              MOB_FCC_0B_ESHELBY,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_FCC_0b_Eshelby
            },

            {"FCC_linear_Eshelby",
              MAT_TYPE_FCC,
              MOB_FCC_LINEAR_ESHELBY,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_FCC_linear_Eshelby
            },
#endif
            {"FCC_angle",
              MAT_TYPE_FCC,
              MOB_FCC_ANGLE,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_FCC_angle
            },

            {"HCP_linear",
              MAT_TYPE_HCP,
              MOB_HCP_LINEAR,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,                 /* non-planar means : NO enforce glide planes */
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_HCP_linear
            },

#ifdef ESHELBY

            {"HCP_linear_Eshelby",
              MAT_TYPE_HCP,
              MOB_HCP_LINEAR_ESHELBY,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,                 /* no planar means : NO enforce glide planes */
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_HCP_linear_Eshelby
            },
#endif

            {"BCC_Ta",
             MAT_TYPE_BCC,
             MOB_TA,
             WARN_ON_MAT_TYPE_MISMATCH,
             PLANAR,
             NO_FUZZY_PLANES,
             CROSS_SLIP,
             NULL,
             Mobility_Ta
            },

            {"BCC_Ta_linear",
             MAT_TYPE_BCC,
             MOB_TA_LINEAR,
             WARN_ON_MAT_TYPE_MISMATCH,
             PLANAR,
             NO_FUZZY_PLANES,
             CROSS_SLIP,
             NULL,
             Mobility_Ta_linear
            },

            {"BCC_Ta_pencil",
             MAT_TYPE_BCC,
             MOB_TA_PENCIL,
             WARN_ON_MAT_TYPE_MISMATCH,
             NON_PLANAR,
             NO_FUZZY_PLANES,
             NO_CROSS_SLIP,
             NULL,
             Mobility_Ta_pencil
            },

            {"Rhombohedral_Va_nl_planar",
              MAT_TYPE_RHOMBOHEDRAL_VA,
              MOB_RHOMBOHEDRAL_VA,
              WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              InitMob_Rhombohedral_Va_nl_planar,
              Mobility_Rhombohedral_Va_nl_planar
            },

            {"Relax",
              MAT_TYPE_BCC,
              MOB_RELAX,
              NO_WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              NO_CROSS_SLIP,
              NULL,
              Mobility_Relax
            },

            {"Relax_Glide",
              MAT_TYPE_BCC,
              MOB_RELAX_GLIDE,
              NO_WARN_ON_MAT_TYPE_MISMATCH,
              PLANAR,
              NO_FUZZY_PLANES,
              CROSS_SLIP,
              NULL,
              Mobility_Relax_glide
            }

        };  /* end of mobAttrList */
#endif  /* if defined INIT_MOB_ATTR_LIST */

#endif /* _MOBILITY_H */
