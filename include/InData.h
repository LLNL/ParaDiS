#pragma once

#ifndef _PDS_IN_DATA_H
#define _PDS_IN_DATA_H

/*****************************************************************************
 *
 *   InData.h  Define the data structures used to hold the input data
 *
 ****************************************************************************/

#include "Typedefs.h"
#include "Home.h"
#include "Tag.h"

/*
 *	During the parallel application initialization, we want to limit
 *	the number of nodes for which data will be held in memory at any
 *	given point.  When that number exceeds the below value, the
 *	current nodal data will be distributed to remote domains.
 */
#define	MAX_NODES_PER_BLOCK	50000

struct _indata {
	Param_t	*param;

	Node_t  *node;		/* array of node structs */
	int	nodeCount;	/* number of nodes in <node> array */

	real8	*burgX;		/* Array of XYZ burgers vector values */
	real8	*burgY;		/* which are now only used to support */
	real8	*burgZ;		/* old format ctrl files which include*/
				/* burger's vector arrays	      */

	int	nburg;		/* number of burgers vectors in the */
                                /* burgX, burgY and burgZ arrays    */

        void    *decomp;        /* pointer to memory allocated to hold   */
                                /* the domain decomposition read in from */
                                /* the restart file(s)                   */
};
#endif
