#pragma once

#ifndef _PDS_CONSTANTS_H
#define _PDS_CONSTANTS_H

/****************************************************************************
 *
 *      Module:      Constants.h  
 *      Description: Contains definitions for some constants that are used
 *                   by other include files.  
 *
 ***************************************************************************/

/*
 *      Define the numbers of nodes allocated when obtaining new node blocks
 *      within portions of the code.
 */
#define NEW_NODEKEY_INC     1000
#define RECYC_NODESTACK_INC 100
#define NODE_BLOCK_COUNT    500

/*
 *      Some definitions used to select whether to do a full or
 *      partial set of force calculations
 */
#define PARTIAL   1
#define FULL      2

/*
 *      Define the subdirectories under which various types of
 *      output files will be created.  These directories are
 *      relative to the user-specified output directory in the
 *      control file.
 */
#define DIR_ARMDATA    "armdata"
#define DIR_DEBUG      "debug"
#define DIR_FLUXDATA   "fluxdata"
#define DIR_FORCE      "force"
#define DIR_GNUPLOT    "gnuplot"
#define DIR_POLEFIG    "polefig"
#define DIR_POVRAY     "povray"
#define DIR_PROPERTIES "properties"
#define DIR_RESTART    "restart"
#define DIR_TECPLOT    "tecplot"
#define DIR_TIMERS     "timers"
#define DIR_VELOCITY   "velocity"
#define DIR_VISIT      "visit"

#define DESCRIPTION_FILE_SUFFIX "desc"
#define GNU_FILE_SUFFIX         "gnu"

/*
 *      Some miscellaneous definitions used throughout the code
 */
#define MAX_NBRS            10 /* maximum number of segments for a node */
#define BOLTZMANNS_CONST    1.38e-23  /* joules/K */
#define FFACTOR_ORTH        0.05
#define FFACTOR_NORMAL      1.0e-4
#define FFACTOR_LMIN        1e-8  /* minimum allowed segment length */
#define FFACTOR_LMIN2       1e-16 /* minimum allowed segment length squared */
#define EPSILON             1e-15

#define MAX_STRING_LEN      512

/*
 *      Define some values that can be used by GenerateOutput()
 *      to indicate the current stage of the execution.  Needed
 *      because output generation differs depending on the point
 *      the program is at.
 */
#define STAGE_INIT      1
#define STAGE_CYCLE     2
#define STAGE_TERM      3

#define FIRST_BLOCK     1
#define LAST_BLOCK      2

/*
 *      Define a set of flags corresponding to the types of
 *      output that may be generated.  Note: these are intended
 *      as bit values that can be or'ed together.
 */
#define GEN_ARM_DATA            0x00001
#define GEN_DENSITY_DATA        0x00002
#define GEN_FLUX_DATA           0x00004
#define GEN_GNUPLOT_DATA        0x00008
#define GEN_POLEFIG_DATA        0x00010
#define GEN_POVRAY_DATA         0x00020
#define GEN_POSTSCRIPT_DATA     0x00040
#define GEN_PROPERTIES_DATA     0x00080
#define GEN_RESTART_DATA        0x00100
#define GEN_TECPLOT_DATA        0x00200
#define GEN_TIMER_DATA          0x00400
#define GEN_VISIT_DATA          0x00800
#define GEN_VELOCITY_DATA       0x01000
#define GEN_XWIN_DATA           0x02000
#define GEN_FORCE_DATA          0x04000
#define GEN_SCR_DATA            0x08000

/*
 *      The DEL_SEG_* definitions are to be used as the <del_seg_factor>
 *      argument in invocations of the ChangeArmBurg() function.  This
 *      value indicates what portion of length of any segment deleted by
 *      the function will be accumulated in the <delSegLength> variable.
 */
#define DEL_SEG_NONE    0.0
#define DEL_SEG_HALF    0.5

/*
 *      Some constants related to HCP materials
 */
#define HCP_MAX_PLANES_PER_BURG 20  /* Define the largest possible number */
                                    /* of glide planes that could be      */
                                    /* associated with an HCP burgers     */
                                    /* vector.  Used in sizing some temp  */
                                    /* arrays.                            */

#define HCP_NUM_GLIDE_BURG  10  /* Number of glide burgers vectors.  There */
                                /* are many more junction burgers vectors  */
                                /* that are not included in this number.   */

#define HCP_NUM_GLIDE_PLANES 4 /* Maximum number of planes for each of  */
                               /* the glide burgers vectors.  There is  */
                               /* one exception: the 10th glide burgers */
                               /* vector has only 3 associated planes.  */

#define HCP_NUM_TOT_BURG    41 /* Total number of HCP burgers vectors.     */
                               /* Includes both glide and junction burgers */
                               /* vectors.                                 */

#define HCP_NUM_TOT_PLANES  426 /* The summed total count of planes    */
                                /* associated with each of the         */
                                /* <HCP_NUM_TOT_BURG> burgers vectors. */

/*
 *      Some constants related to FCC materials
 */
#define FCC_NUM_GLIDE_BURG    6     /* Number of glide burgers vectors.  There  */
                                    /* are more junction burgers vectors  that  */
                                    /* are not included in this number.         */

#define FCC_NUM_GLIDE_PLANES  2     /* Maximum number of planes for each of     */
                                    /* the glide burgers vectors.               */

#define FCC_NUM_TOT_BURG    21      /* Total number of FCC burgers vectors.     */
                                    /* Includes both glide and junction burgers */
                                    /* vectors.                                 */

#define FCC_NUM_TOT_PLANES  54      /* The summed total count of planes         */
                                    /* associated with each of the              */
                                    /* burgers vectors.                         */

#define FCC_SIZE_DENS       10      /* See Ltot contents in DeltaPlasticStrain_FCC */
#define FCC_SIZE_FLUX        5      /* See AreaSwept contents in DeltaPlasticStrain_FCC */

/*
 *      Some constants related to BCC materials
 */
#define BCC_NUM_GLIDE_BURG    4     /* Number of glide burgers vectors.  There  */
                                    /* are more junction burgers vectors  that  */
                                    /* are not included in this number.         */

#define BCC_NUM_GLIDE_PLANES  6     /* Maximum number of planes for each of     */
                                    /* the glide burgers vectors.               */


#define BCC_NUM_TOT_BURG    13      /* Total number of FCC burgers vectors.     */
                                    /* Includes both glide and junction burgers */
                                    /* vectors.                                 */

#define BCC_NUM_TOT_PLANES 144      /* The summed total count of planes         */
                                    /* associated with each of the glide        */
                                    /* burgers vectors.                         */


#define FOUR_PI          1.256637061435917e+01
#define INV_4_PI_SQUARED 0.025330295910584  /* 1/(4*pi^2) */
#define INV_4_OVER_PI    0.079577471545948  /* 1/4/pi */

/*
 *      In order to cut down on memory usage, we created hash tables for
 *      both the standard and FMM cells rather than the old method of
 *      keeping an array of cell pointers for *all* cells.  This defines
 *      the number of base elements in the cell hash tables.
 */
#define CELL_HASH_TABLE_SIZE 197
#endif  /* _CONSTANTS_H */
