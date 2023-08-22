/****************************************************************************
 *
 *	Function:	InitRemesh
 *	Description:	If any segments in inData are bigger than
 *			param->maxSeg, break them up into smaller segments.
 *			
 *			Note: This function requires that the InData node
 *			      array be sorted in ascending order by the
 *			      nodes' tag values prior to entry.  The
 *			      array will be returned to the caller sorted
 *			      in the same manner.
 *
 ***************************************************************************/
#include "Home.h"
#include "InData.h"
#include "Param.h"
#include "Util.h"
#include "ParadisGen.h"


/*-------------------------------------------------------------------------
 *
 *      Function:      LookupInNode
 *      Description:   Search through the InData->nodes array of
 *                     Node_t structures for the node whose tag
 *                     value matches the provided tag
 *      Args:
 *          inData     pointer to structure containing the node array
 *                     to be searched
 *          tag        the <oldTag> value for the node to be located
 *          index      pointer to integer in which to return the node
 *                     array index of the node found.  If the desired
 *                     node is not found, a -1 will be returned to the
 *                     caller here.
 *
 *      Returns:    pointer to the node structure whose oldTag value
 *                  matched the provided tag, or NULL if no matching
 *                  node was found.
 *
 *------------------------------------------------------------------------*/
Node_t *LookupInNode(InData_t *inData, Tag_t *tag, int *index)
{
        Node_t *found;
        Node_t tmpNode;

        tmpNode.myTag = *tag;

        found = (Node_t *) bsearch(&tmpNode, inData->node, inData->nodeCount,
                                   sizeof(Node_t), OrderNodes);

        if (found != (Node_t *)NULL) {
            *index = found - inData->node;
        } else {
            *index = -1;
        }

        return(found);
}

void InitRemesh(InData_t *inData, int domValue, int startIndex)
{
	int		oldNodeCount, newNodeCount;
	int		idxA, idxB, idxBA, armAB, armBA;
	int		idxNew, ncount, found, armPrev, inew;
	real8		maxSeg2, r2, r, xa, ya, za, frac;
	real8		dx, dy, dz;
	Node_t          *nodeA, *nodeB, *prevNewNode, *newNode;
	Param_t		*param;


	param = inData->param;

	maxSeg2 = param->maxSeg * param->maxSeg;

/*
 *	Loop once thru the nodes in inData just to count the number of
 *	new nodes that will be added
 */
	oldNodeCount = inData->nodeCount;
	newNodeCount = 0;

	for (idxA = startIndex+1; idxA < oldNodeCount; idxA++) {

		nodeA = &inData->node[idxA];
		for (armAB = 0; armAB < nodeA->numNbrs; armAB++) {

			nodeB = LookupInNode(inData, &nodeA->nbrTag[armAB],
                                             &idxB);
			if (idxB < idxA) {

/*
 *				Count the number of nodes to stick in between
 *				nodeA and its neighbor
 */

				dx = nodeB->x - nodeA->x; 
				dy = nodeB->y - nodeA->y; 
				dz = nodeB->z - nodeA->z; 

				ZImage(param,&dx,&dy,&dz);

				r2 = dx*dx + dy*dy + dz*dz;
				if (r2 > maxSeg2) {
					r = sqrt(r2);
					newNodeCount +=
						(int)ceil(r/param->maxSeg) - 1;
				}

			}
		}
	}

/*
 *	If no new nodes are required, just return to the caller, otherwise
 *	expand the inData->node array to accomodate the new nodes.
 */
	if (newNodeCount == 0) return;

	inData->node = (Node_t *)realloc(inData->node,
			sizeof(Node_t) * (inData->nodeCount + newNodeCount));

        memset((char *)(&inData->node[inData->nodeCount]), 0,
               sizeof(Node_t) * newNodeCount);

/*
 *	Loop thru original nodes again, this time inserting new nodes
 *	between them as necessary, in order to reduce segment lengths
 *	to within param->maxSeg.
 */
	idxNew = oldNodeCount;  /* index of first added node */
	for (idxA = startIndex+1; idxA < oldNodeCount; idxA++) {

		nodeA = &inData->node[idxA];
		xa = nodeA->x; ya = nodeA->y; za = nodeA->z;

		for (armAB = 0; armAB < nodeA->numNbrs; armAB++) {

			nodeB = LookupInNode(inData, &nodeA->nbrTag[armAB], &idxB);
			if (idxB > idxA) continue; /* don't double-dip */

			dx = nodeB->x - xa;
			dy = nodeB->y - ya;
			dz = nodeB->z - za;

			ZImage(param,&dx,&dy,&dz);

			r2 = dx*dx + dy*dy + dz*dz;
			if (r2 > maxSeg2) {

/*
 *				calculate the number of new nodes needed
 *				between the two points.
 */
				r = sqrt(r2);
				ncount = (int) ceil(r/param->maxSeg) - 1;
				frac = 1.0 / (ncount+1);

/*
 *				Determine which arm of nodeB points to nodeA
 */
				found = 0;
				for (armBA = 0;armBA < nodeB->numNbrs;armBA++) 
                {
					LookupInNode(inData, &nodeB->nbrTag[armBA], &idxBA);

					if (idxBA==idxA) { found=1; break; }
				}

				if (!found) {
					Fatal("InitRemesh: broken Linkage");
				}

/*
 *				Add new nodes between nodeA and nodeB
 */
				prevNewNode = nodeA;
				armPrev = armAB;

				for (inew = 1; inew <= ncount; inew++) {

					newNode = &inData->node[idxNew];

					newNode->myTag.domainID = domValue;
					newNode->myTag.index =
                                                idxNew + param->nodeCount;

					newNode->x = xa + inew*frac*dx;
					newNode->y = ya + inew*frac*dy;
					newNode->z = za + inew*frac*dz;
					SET_CONSTRAINTS(newNode->constraint,
                                                        UNCONSTRAINED);

					AllocNodeArms(newNode, 2);
	
/*
 *					set up the linkage:
 *					prevNewNode <--> newNode <--> nodeB
 */
					newNode->nbrTag[0].domainID =
						prevNewNode->myTag.domainID;
					newNode->nbrTag[0].index =
						prevNewNode->myTag.index;
					newNode->nbrTag[1].domainID =
						nodeB->myTag.domainID;
					newNode->nbrTag[1].index =
						nodeB->myTag.index;

					prevNewNode->nbrTag[armPrev].domainID =
						domValue;
					prevNewNode->nbrTag[armPrev].index =
                                                idxNew + param->nodeCount;

					nodeB->nbrTag[armBA].domainID = domValue;
					nodeB->nbrTag[armBA].index =
                                                idxNew + param->nodeCount;
					armPrev = 1;

/*
 *					Initialize newNode's arms
 */
					newNode->burgX[0] = nodeB->burgX[armBA];
					newNode->burgY[0] = nodeB->burgY[armBA];
					newNode->burgZ[0] = nodeB->burgZ[armBA];
					newNode->nx[0] = nodeB->nx[armBA];
					newNode->ny[0] = nodeB->ny[armBA];
					newNode->nz[0] = nodeB->nz[armBA];

					newNode->burgX[1] = nodeA->burgX[armAB];
					newNode->burgY[1] = nodeA->burgY[armAB];
					newNode->burgZ[1] = nodeA->burgZ[armAB];
					newNode->nx[1] = nodeA->nx[armAB];
					newNode->ny[1] = nodeA->ny[armAB];
					newNode->nz[1] = nodeA->nz[armAB];

					prevNewNode = newNode;
					idxNew++;

				}  /* end for (inew = 1; ...)  */
			}  /* end if (r2 > maxSeg2)  */
		}  /* end for (armAB = 0; ...)  */
	}  /* end for idxA = 1; ...)  */

	inData->nodeCount += newNodeCount;

	if (idxNew != inData->nodeCount)
		Fatal("InitRemesh: miscount of added nodes!");

/*
 *	Must re-sort the node list in ascending order based on the
 *	tag values.  (This is required because LookupInNode() uses
 *	a binary search to locate nodes by their tags.)
 */
	qsort(&inData->node[startIndex], inData->nodeCount-startIndex,
              sizeof(Node_t), OrderNodes);

	return;
}
