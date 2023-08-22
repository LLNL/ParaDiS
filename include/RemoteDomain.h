#pragma once

#ifndef _PDS_REMOTE_DOMAIN_H
#define _PDS_REMOTE_DOMAIN_H

/****************************************************************************
 *
 *	RemoteDomain.h	Define the struct that holds various data used in 
 *			communicating with a neighboring domain
 *
 ****************************************************************************/

#include "Typedefs.h"
#include "Tag.h"

struct _remotedomain {
	int	domainIdx;	/* encoded index of this domain */
	int	numExpCells;	/* number of native cells exported to */
				/* this domain                        */
	int	*expCells;	/* list of encoded indices of the */
				/* exported cells                 */

	int	maxTagIndex;	/* error if tag.index exceeds this value */
	Node_t	**nodeKeys;	/* indexed by node's tag.index, points */
				/* to Node_t                           */

	int	inBufLen;
	char	*inBuf;
	int	outBufLen;
	char	*outBuf;
};

#endif
