/****************************************************************************
 * 
 *      Module:       SegPartIntersect.c
 *
 *      Description:  Contains routines used for maintaining and
 *                    searching the list of inclusion-intersecting
 *                    segments created during the force calculations
 *                    for use with mobility functions that differentiate
 *                    between mobilities inside and outside inclusions.
 *
 *      Includes public functions:
 *          GetIndexFromInclusionID()
 *          FindNodePartIntersects()
 *          SegPartListAppend()
 *          SegPartListClear()
 *          SegPartListInsert
 *          SegPartListLookup()
 *          SegPartListSort()
 *          SegPartListUnsortedLookup()
 *          SegPartListUpdate()
 *
 *      Includes local fucntions()
 *          FindSegPartIntersects()
 *          SegPartListCmp()
 *
 ****************************************************************************/

#ifdef ESHELBY

#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "Mobility.h"
#include "M33.h"
#include "Topology.h"

/*---------------------------------------------------------------------------
 *
 *      Function:     IntersectionSphereSegment
 *      Description:  Determine the number of intersecting points between
 *                    a sphere and a line segment, the fraction of the
 *                    segment inside the sphere, the the derivative of
 *                    these fractions.
 *
 *      Note : Before entering this routine, care must be made to ensure 
 *             that PBC is taken into account for pos1, pos2, and C before 
 *             rotation in the case of an ellipse.
 * 
 *      Arguments:
 *          p1 : node starting the segment
 *          p2 : node ending the segment
 *          C  : sphere center
 *          r  : sphere radius
 *          n  : number of intersection points
 *          ratio : fractions of the line segments where the intersections points
 *                   are. These fractions can be outside the range of [0, 1].
 *
 *-------------------------------------------------------------------------*/
static int IntersectionSphereSegment(real8 p1[3], real8 p2[3],real8 C[3], 
                                     real8 r, int *n, real8 ratio[2])
{
   real8 v[3], t[3], R[2][3], s[2], nd[3], sinter[2];  
   real8 invL, d2, test, sdiff;
   (*n) = 0;

   int isInter = 0;
   ratio[0] = 0.0;
   ratio[1] = 0.0;

   v[X] = p2[X]-p1[X];
   v[Y] = p2[Y]-p1[Y];
   v[Z] = p2[Z]-p1[Z];

   invL = 1./Normal(v);
   t[X] = v[X]*invL;
   t[Y] = v[Y]*invL;
   t[Z] = v[Z]*invL;

   R[0][X] = p1[X] - C[X];
   R[0][Y] = p1[Y] - C[Y];
   R[0][Z] = p1[Z] - C[Z];

   R[1][X] = p2[X] - C[X];
   R[1][Y] = p2[Y] - C[Y];
   R[1][Z] = p2[Z] - C[Z];

   s[0] = DotProduct(R[0],t);
   s[1] = DotProduct(R[1],t);

   nd[X] = R[0][X] - s[0]*t[X];
   nd[Y] = R[0][Y] - s[0]*t[Y];
   nd[Z] = R[0][Z] - s[0]*t[Z];

   d2 = DotProduct(nd,nd);
   test = r*r-d2;

   if (test < 0.0) return 0;
   
   sinter[0] = -sqrt(test);
   sinter[1] = -sinter[0];
   
   sdiff = s[1]-s[0];
   ratio[0] = (sinter[0] - s[0])/sdiff;
   ratio[1] = (sinter[1] - s[0])/sdiff;

   if (ratio[0]>0.0 && ratio[0]<1) (*n)++;
   if (ratio[1]>0.0 && ratio[1]<1) (*n)++;
   if ( (*n) > 0) 
      isInter = 1;
   else
      isInter = 0;

   // The whole segment could be inside the particle.
   // This counts as an intersection.
   if (ratio[0]<0.0 && ratio[1]>1) isInter = 1;

   ratio[0] = MIN(1.0, MAX(ratio[0], 0.0));
   ratio[1] = MIN(1.0, MAX(ratio[1], 0.0));
   
   return isInter;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     TransformEllips2Sphere
 *      Description:  Matrix defining the rotation and scaling of coordinate
 *                    system from an ellipsoid to a sphere. Used in the 
 *                    mobility laws called Eshelby for now.
 *                    The matrix B is defined as the product between the diagonal 
 *                    matrix of the semi-principal axis multiplied with the rotation
 *                    matrix of these axis
 * 
 *      Arguments:
 *            rotation[2][3] : two vectors defining the rotation matrix
 *            radius[3]      : semi-principal axis
 *             
 *            B[3][3]        : The matrix defining the transformation from ellipoidal
 *                             coordinates to spherical coordinates.
 *
 *---------------------------------------------------------------------------*/
static void TransformEllips2Sphere(real8 rotation[2][3], real8 radius[3], real8 B[3][3])
{
  real8 vec[3];
  NormalizedCrossVector(rotation[0],rotation[1],vec);

  B[0][0] = rotation[0][0]/radius[0];
  B[0][1] = rotation[0][1]/radius[0];
  B[0][2] = rotation[0][2]/radius[0];
  
  B[1][0] = rotation[1][0]/radius[1];
  B[1][1] = rotation[1][1]/radius[1];
  B[1][2] = rotation[1][2]/radius[1];
  
  B[2][0] = vec[0]/radius[2];
  B[2][1] = vec[1]/radius[2];
  B[2][2] = vec[2]/radius[2];

}

/*---------------------------------------------------------------------------
 *
 *      Function:     IntersectionEllipseSegment
 *      Description:  Determine the number of intersecting points between the sphere
 *                    and the line segment and the fraction of the segment inside the 
 *                    sphere.
 *
 *      Note : Before entering this routine, care must be made to ensure that PBC is 
 *             taken into account for pos1, pos2, and C before rotation in the case 
 *             of an ellipse.
 * 
 *      Arguments:
 *          C         : position of the ellipse
 *          Radius    : Sphere radius if ellispe is a sphere. Set to 1 in the case of 
 *                      an ellipse
 *          Rotation  : Two vectors defining the rotation matrix
 *          p1        : node starting the segment
 *          p2        : node ending the segment
 *          ninc      : number of intersection points
 *          ratio     : fractions of the line segments where the intersections points
 *                      are. These fractions can be outside the range of [0, 1].
 *
 *-------------------------------------------------------------------------*/
int IntersectionEllipseSegment(real8  C[3], real8 Radius[3], real8 Rotation[2][3],
                               real8 p1[3], real8 p2[3], real8 r, int *ninc, 
                               real8 ratio[2])
{
   int IsIntersect = 0;
   real8 B[3][3], Rp1[3], Rp2[3], RC[3];

   // Transform the ellipse into a sphere to get fractions of dislocation line segment
   // inside the ellipse
   TransformEllips2Sphere(Rotation, Radius, B);
   Matrix33Vector3Multiply(B,p1,Rp1);
   Matrix33Vector3Multiply(B,p2,Rp2);
   Matrix33Vector3Multiply(B,C,RC);
   IsIntersect = IntersectionSphereSegment(Rp1, Rp2, RC, r, ninc, ratio);
   return IsIntersect;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     IntersectionInformation
 *      Description:  Determine information between a given dislocation segment
 *                    and a given ellipse, take care of PBC conditions.
 *
 *      Arguments:
 *-------------------------------------------------------------------------*/
int IntersectionInformation(Home_t *home, real8 pos1[3], real8 pos2[3], EInclusion_t *inclusion, 
                            real8 newPos1[3], real8 newPos2[3], int *ninc, real8 ratio[2])
{
   int IsIntercept=0;
   real8 cellCtr[3];

   Param_t *param;
   Cell_t *cell;

   param = home->param;

   // Define the particle's info
   int cellID;
   real8 C[3], radius[3], Rotation[2][3];
   cellID = inclusion->cellID;
   VECTOR_COPY(C,inclusion->position);
   VECTOR_COPY(radius,inclusion->radius);
   VECTOR_COPY(Rotation[0],inclusion->rotation[0]);
   VECTOR_COPY(Rotation[1],inclusion->rotation[1]);
  
   // Ensure PBC 
   newPos1[X] = pos1[X]; newPos1[Y] = pos1[Y]; newPos1[Z] = pos1[Z];
   newPos2[X] = pos2[X]; newPos2[Y] = pos2[Y]; newPos2[Z] = pos2[Z];
   
   cell = LookupCell(home, cellID);
   
   if (cell == (Cell_t *)NULL) {
      FindCellCenter(param, C[X], C[Y], C[Z], 1,
                     &cellCtr[X], &cellCtr[Y], &cellCtr[Z]);
   } else {
      cellCtr[X] = cell->center[X];
      cellCtr[Y] = cell->center[Y];
      cellCtr[Z] = cell->center[Z];
   }
   
   PBCPOSITION(param, cellCtr[X], cellCtr[Y], cellCtr[Z],
               &newPos1[X], &newPos1[Y], &newPos1[Z]);

   PBCPOSITION(param, newPos1[X], newPos1[Y], newPos1[Z],
               &newPos2[X], &newPos2[Y], &newPos2[Z]);
   
   IsIntercept = IntersectionEllipseSegment(C, radius, Rotation, newPos1, newPos2, 
                                            1.0, ninc, ratio);

   return IsIntercept;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     GetIndexFromInclusionID()
 *      Description:  Locate the specified inclusion in the local domain's
 *                    inclusion list and return the index in the list
 *                    to the caller.  
 *
 *                    NOTE: The inclusions are NOT sorted, so unless we
 *                    modify the code to do so, we have to do a linear
 *                    search through the list for the inclusion.
 *
 *      Arguments:
 *          inclusionID  self-explanatory
 *
 *      Returns:  The index (0-based) in the local inclusion list of the
 *                inclusion with the provided ID.
 *
 *-------------------------------------------------------------------------*/
int GetIndexFromInclusionID(Home_t *home, int inclusionID)
{
        int                numEntries;
        SegPartIntersect_t *entry, *tmpList;

        // Access to segment/particle intersection list
        numEntries = home->segpartIntersectCnt;
        tmpList = home->segpartIntersectList;

        // Loop over all entries in the seg/part list
        for (int i = 0; i <numEntries ; i++) 
        {
           entry = &tmpList[i];           
           if (entry == (SegPartIntersect_t *)NULL)  return(-1);
           
           for (int j = 0; j < (entry)->numIntersections ; j++) 
           {
              if (inclusionID == entry->inclusionID[j])
              {
                 return (entry->inclusionIndex[j]);
              }
           }
        }

        return(-1);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListAppend()
 *      Description:  Append the provided segment/particle intersection
 *                    information to the segment/particle intersection
 *                    list.  Caller must resort the list before using
 *                    the SegPartListLookup() function.
 *
 *      Arguments:
 *          newInfo  Pointer to structure particle intersection data
 *                   for a segment.  Data is copied from this structure
 *                   into the segment/particle intersction list.
 *
 *-------------------------------------------------------------------------*/
void SegPartListAppend(Home_t *home, SegPartIntersect_t *newInfo)
{
        int                i, numEntries, allocedEntries;
        SegPartIntersect_t *intersectInfo, *tmpList;

/*
 *      The segment/particle intersection list is globally
 *      available to all threads, so we have to ensure only
 *      one thread accesses it at a time
 */
#ifdef _OPENMP
#pragma omp critical (ACCESS_SEG_PART_LIST)
#endif
        {
            numEntries = home->segpartIntersectCnt;
            allocedEntries = home->segpartIntersectListSize;
            tmpList = home->segpartIntersectList;

/*
 *          Make sure there is enough space on the list for a segment
 */
            if (numEntries >= allocedEntries) {
                allocedEntries += 100;
                tmpList = (SegPartIntersect_t *)realloc(tmpList,
                           allocedEntries * sizeof(SegPartIntersect_t));
            }

/*
 *          Get a pointer to the next free entry in the list and
 *          initialize it.  IMPORTANT: Make sure that the node tags
 *          identifying the segment are stored such that the first tag
 *          for a segment is "less than" the second as determined
 *          by OrderTags().  This is so we're consistent in the manner
 *          we sort the list and search for entries on the list later.
 */
            intersectInfo = &tmpList[numEntries];
            memset(intersectInfo, 0, sizeof(SegPartIntersect_t));

            if (OrderTags(&newInfo->tag1, &newInfo->tag2) < 0) {
                intersectInfo->tag1.domainID = newInfo->tag1.domainID;
                intersectInfo->tag1.index    = newInfo->tag1.index;
                intersectInfo->tag2.domainID = newInfo->tag2.domainID;
                intersectInfo->tag2.index    = newInfo->tag2.index;
            } else {
                intersectInfo->tag1.domainID = newInfo->tag2.domainID;
                intersectInfo->tag1.index    = newInfo->tag2.index;
                intersectInfo->tag2.domainID = newInfo->tag1.domainID;
                intersectInfo->tag2.index    = newInfo->tag1.index;
            }

/*
 *          Probably don't need to copy the inclusion IDs here since
 *          this function should only be called when reading the intersect
 *          data from a remote domain and will not be sending it any
 *          further, but...
 */
            for (i = 0; i < newInfo->numIntersections; i++) {
                intersectInfo->inclusionIndex[i] = newInfo->inclusionIndex[i];
                intersectInfo->inclusionID[i] = newInfo->inclusionID[i];
            }

            home->segpartIntersectCnt = numEntries + 1;
            home->segpartIntersectListSize = allocedEntries;
            home->segpartIntersectList = tmpList;

        }  /* end "omp critical" section */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FindSegPartIntersects
 *      Description:  Determine which inclusions (particles) intersect
 *                    the specified segment and add the info to the
 *                    seg/particle intersection list.  
 * 
 *                    NOTE:  This function does NOT verify that the
 *                           list does not already contain the intersection
 *                           data for the specified segment.  It is up to
 *                           the caller to ensure he does not call this
 *                           function for the same segment multiple times.
 *
 *      Arguments:
 *          tag1    Tag of one of the segment endpoints
 *          tag2    Tag for the other segment endpoint
 *          pos1    Coordinates of the node associated with tag1
 *          pos2    Coordinates of the node associated with tag2, adjusted
 *                  to the periodic image closest to pos1.
 *          cellID  ID of the cell containing pos1
 *          cellCtr Coordinates of the center of cell <cellID>
 *
 *-------------------------------------------------------------------------*/
static void FindSegPartIntersects(Home_t *home, Tag_t *tag1, Tag_t *tag2,
                                  real8 pos1[3], real8 pos2[3], int cellID,
                                  real8 cellCtr[3])
{
        int                j, numNbrCells, inclusionIndex, inclusionCellID;
        real8              sinter[2];
        real8              incPos[3];
        Cell_t             *cell, *inclusionCell;
        Param_t            *param;
        EInclusion_t       *inclusion;
        SegPartIntersect_t *intersection = (SegPartIntersect_t *)NULL;

        param = home->param;
        cell = LookupCell(home, cellID);
        numNbrCells = cell->nbrCount;
        intersection = (SegPartIntersect_t *)NULL;

/*
 *      First loop over all neighboring cells, then on last loop
 *      iteration do the current cell.
 */
        for (j = 0; j <= numNbrCells; j++) {

            if (j < numNbrCells) {
                inclusionCellID = cell->nbrList[j];
            } else {
                inclusionCellID = cellID;
            }

/*
 *          It *may* be possible that we don't have cell info for
 *          all the cells neighboring the provided <cellID>.  In
 *          that case we just do the best we can.
 */
            inclusionCell = LookupCell(home, inclusionCellID);
            if (inclusionCell == (Cell_t *)NULL) continue;

/*
 *          NOTE: If the cell is a ghost cell, alter the cellID and cell
 *          pointer to the corresponding base cell.  Needed in order to pick
 *          up the proper inclusionQ.
 */
            inclusionCell = LookupCell(home, inclusionCellID);

            if (inclusionCell->baseIdx >= 0) {
                inclusionCellID = inclusionCell->baseIdx;
                inclusionCell = LookupCell(home, inclusionCellID);
                if (inclusionCell == (Cell_t *)NULL) continue;
            }

            inclusionIndex = inclusionCell->inclusionQ;

/*
 *          Loop through all the inclusion in this cell
 */
            while (inclusionIndex >= 0) {

                int ninc;

                inclusion = &home->eshelbyInclusions[inclusionIndex];


/*
 *              Adjust the center of the inclusion to that of the periodic
 *              image closest to the center of the cell owning the first node.
 */
                incPos[X] = inclusion->position[X];
                incPos[Y] = inclusion->position[Y];
                incPos[Z] = inclusion->position[Z];

                PBCPOSITION(param, cellCtr[X], cellCtr[Y], cellCtr[Z],
                            &incPos[X], &incPos[Y], &incPos[Z]);

                int IsIntercept = 0;
                IsIntercept = IntersectionEllipseSegment(incPos, inclusion->radius, inclusion->rotation,
                                                         pos1, pos2, 1.0, &ninc, sinter);

/*
 *              If the segment does not intersect the inclusion, i.e number 
 *              of intersecting points is zero, skip it
 */
                if (IsIntercept)
                {
                   SegPartListUpdate(home, &intersection, tag1, tag2,
                                     inclusionIndex, inclusion->id);
                }

                inclusionIndex = inclusion->nextInCell;
            }

        }  /* loop over cells */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     FindNodePartIntersects
 *      Description:  Given a node, find out which inclusions are intersected
 *                    by any of the segments attached to the node.  The
 *                    segment/particle intersect list will be updated with
 *                    the info.
 *
 *                    NOTE:  It is up to the caller to sort the seg/particle
 *                    intersection list after he's done generating it.
 *
 *
 *      Returns:   Pointer to structure of segment/particle intersection
 *                 information, or NULL if segment does not intersect any
 *                 inclusions.
 *
 *-------------------------------------------------------------------------*/
void FindNodePartIntersects(Home_t *home, Node_t *node)
{
        int     i, cellID;
        real8   eps = 1.0e-6, dx, dy, dz;
        real8   pos1[3], pos2[3], cellCtr[3];
        Node_t  *nbrNode;
        Cell_t  *cell;
        SegPartIntersect_t *intersection;

        pos1[X] = node->x;
        pos1[Y] = node->y;
        pos1[Z] = node->z;

        cellID = node->cellIdx;
        cell = LookupCell(home, cellID);

        cellCtr[X] = cell->center[X];
        cellCtr[Y] = cell->center[Y];
        cellCtr[Z] = cell->center[Z];

        for (i = 0; i < node->numNbrs; i++) {

            if ((nbrNode = GetNeighborNode(home, node, i)) == (Node_t *)NULL) {
                continue;
            }

            dx = nbrNode->x - node->x;
            dy = nbrNode->y - node->y;
            dz = nbrNode->z - node->z;
            
            ZImage(home->param, &dx, &dy, &dz);


            pos2[X] = node->x + dx;
            pos2[Y] = node->y + dy;
            pos2[Z] = node->z + dz;


/*
 *          Zero length arms can occur when we're testing multinode splits.
 *          Ignore such segments.
 */
            if ((fabs(pos1[X]-pos2[X]) < eps) &&
                (fabs(pos1[Y]-pos2[Y]) < eps) &&
                (fabs(pos1[Z]-pos2[Z]) < eps)) {
                continue;
            }

/*
 *          If this function is called for two connected nodes, we have to
 *          be sure we don't search for intersecting particles multiple
 *          times, so only add the segment/particle intersections if the
 *          segment is not yet on the list.  (Note: list may be unsorted at
 *          this point, so we have to use a lookup function that handles
 *          an unsorted list)
 */
            intersection = SegPartListUnsortedLookup(home, &node->myTag,
                                                     &nbrNode->myTag);
            if (intersection != (SegPartIntersect_t *)NULL) {
                continue;
            }

            FindSegPartIntersects(home, &node->myTag, &nbrNode->myTag,
                                  pos1, pos2, cellID, cellCtr);
            
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListCmp()
 *      Description:  Compares two SegPartIntersect_t structures based
 *                    on the <tag1> and <tag2> values.  This function
 *                    is compatible with use by system sorting and
 *                    searching functions such as qsort(), bsearch(), etc.
 *
 *      Returns:        -1 if  a <  b
 *                       0 if  a == b
 *                       1 if  a >  b
 *
 *-------------------------------------------------------------------------*/
static int SegPartListCmp(const void *a, const void *b)
{
        SegPartIntersect_t *info1 = (SegPartIntersect_t *)a;
        SegPartIntersect_t *info2 = (SegPartIntersect_t *)b;

/*
 *      First compare first tag 
 */
        if (info1->tag1.domainID < info2->tag1.domainID) return(-1);
        if (info1->tag1.domainID > info2->tag1.domainID) return(1);

        if (info1->tag1.index < info2->tag1.index) return(-1);
        if (info1->tag1.index > info2->tag1.index) return(1);

/*
 *      Both entries have the same first tag, so compare the second
 *      tag of each.
 */
        if (info1->tag2.domainID < info2->tag2.domainID) return(-1);
        if (info1->tag2.domainID > info2->tag2.domainID) return(1);

        if (info1->tag2.index < info2->tag2.index) return(-1);
        if (info1->tag2.index > info2->tag2.index) return(1);

        return(0);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListSort()
 *      Description:  Sort the list of particle-intersecting segments
 *
 *-------------------------------------------------------------------------*/
void SegPartListSort(Home_t *home)
{
        if (home->segpartIntersectCnt > 1) {
            qsort(home->segpartIntersectList, home->segpartIntersectCnt,
                  sizeof(SegPartIntersect_t), SegPartListCmp);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListLookup()
 *      Description:  Given a pair of node tags defining a segment, return
 *                    to the caller the corresponding entry (if any) in the
 *                    segment/particle intersection list.
 *
 *      Arguments:
 *          tag1  Node tag of one of the segment endpoints
 *          tag2  Node tag of the other segment endpoint
 *
 *      Returns:  Pointer to the entry in the segment/particle
 *                intersection list with endpoints <tag1> and <tag2>
 *                or NULL if the segment is not found on the list
 *
 *-------------------------------------------------------------------------*/
SegPartIntersect_t *SegPartListLookup(Home_t *home, Tag_t *tag1, Tag_t *tag2)
{
        SegPartIntersect_t key;
        SegPartIntersect_t *entry = (SegPartIntersect_t *)NULL;

/*
 *      Node tags in the intersection list are stored such that the
 *      first tag for a segment is "less than" the second as determined
 *      by OrderTags(), so make sure we set up the tags for the search
 *      key the same way.
 */
        if (OrderTags(tag1, tag2) < 0) {
            key.tag1.domainID = tag1->domainID;
            key.tag1.index    = tag1->index;
            key.tag2.domainID = tag2->domainID;
            key.tag2.index    = tag2->index;
        } else {
            key.tag1.domainID = tag2->domainID;
            key.tag1.index    = tag2->index;
            key.tag2.domainID = tag1->domainID;
            key.tag2.index    = tag1->index;
        }

        entry = (SegPartIntersect_t *)bsearch(&key, home->segpartIntersectList,
                                              home->segpartIntersectCnt,
                                              sizeof(SegPartIntersect_t),
                                              SegPartListCmp);

        return(entry);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListUnsortedLookup()
 *      Description:  Given a pair of node tags defining a segment, return
 *                    to the caller the corresponding entry (if any) in the
 *                    segment/particle intersection list.  This function
 *                    does NOT assume the list is currently sorted and 
 *                    does a linear search through the array.  This is
 *                    only intended to be used when dealing with a very
 *                    limited set of segments, and so shouldn't be very
 *                    expensive.
 *
 *      Arguments:
 *          tag1  Node tag of one of the segment endpoints
 *          tag2  Node tag of the other segment endpoint
 *
 *      Returns:  Pointer to the entry in the segment/particle
 *                intersection list with endpoints <tag1> and <tag2>
 *                or NULL if the segment is not found on the list
 *
 *-------------------------------------------------------------------------*/
SegPartIntersect_t *SegPartListUnsortedLookup(Home_t *home, Tag_t *tag1,
                                              Tag_t *tag2)
{
        int                i;
        Tag_t              *tmpTag1, *tmpTag2;
        SegPartIntersect_t *entry = (SegPartIntersect_t *)NULL;

        if (OrderTags(tag1, tag2) < 0) {
           tmpTag1 = tag1;
           tmpTag2 = tag2;
        } else {
           tmpTag2 = tag1;
           tmpTag1 = tag2;
        }

        for (i = 0; i < home->segpartIntersectCnt; i++) {

            entry = &home->segpartIntersectList[i];

            if ((tmpTag1->domainID == entry->tag1.domainID) &&
                (tmpTag1->index    == entry->tag1.index   ) &&
                (tmpTag2->domainID == entry->tag2.domainID) &&
                (tmpTag2->index    == entry->tag1.index   )) {

                return(entry);
            }
        }

        return((SegPartIntersect_t *)NULL);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListUpdate()
 *      Description:  Allocate (if necessary) a new entry in the segment/
 *                    particle intersection list, initialize the
 *                    segment information, and return a pointer to the
 *                    entry to the caller.
 *
 *      Arguments:
 *          entry          Location containing pointer to the intersection
 *                         data to be updated.  If contents are NULL, a new
 *                         entry will be added to the segment/particle
 *                         intersection list and the address of the new
 *                         structure returned to the caller in <entry>
 *                         
 *          tag1           Node tag of one of the segment endpoints
 *          tag2           Node tag of the other segment endpoint
 *          inclusionIndex Index in the local inclusion list of the 
 *                         particle intersecting the specified segment.
 *          inclusionID    ID of the particle intersecting the specified
 *                         segment.
 *
 *-------------------------------------------------------------------------*/
void SegPartListUpdate(Home_t *home, SegPartIntersect_t **entry, Tag_t *tag1,
                       Tag_t *tag2, int inclusionIndex, int inclusionID)
{
        int                nextIndex, numEntries, allocedEntries;
        SegPartIntersect_t *intersectInfo, *tmpList;

/*
 *      For threaded processes, we have to synchronize updates to the
 *      seg/particle intersection list via locks.
 */
        LOCK(&home->segpartIntersectListLock);

/*
 *      if we don't have an entry to which to add the intersect info,
 *      insert a new entry into the list.
 */
        if (*entry == (SegPartIntersect_t *)NULL) {

            numEntries = home->segpartIntersectCnt;
            allocedEntries = home->segpartIntersectListSize;
            tmpList = home->segpartIntersectList;

/*
 *          Make sure there is enough space on the list for a segment
 */
            if (numEntries >= allocedEntries) {
                allocedEntries += 100;
                tmpList = (SegPartIntersect_t *)realloc(tmpList,
                          allocedEntries * sizeof(SegPartIntersect_t));
            }

/*
 *          Get a pointer to the next free entry in the list and
 *          initialize it.  IMPORTANT: Make sure that the node tags
 *          identifying the segment are stored such that the first tag
 *          for a segment is "less than" the second as determined
 *          by OrderTags().  This is so we're consistent in the manner
 *          we sort the list and search for entries on the list later.
 */
            intersectInfo = &tmpList[numEntries];
            memset(intersectInfo, 0, sizeof(SegPartIntersect_t));

            if (OrderTags(tag1, tag2) < 0) {
                intersectInfo->tag1.domainID = tag1->domainID;
                intersectInfo->tag1.index    = tag1->index;
                intersectInfo->tag2.domainID = tag2->domainID;
                intersectInfo->tag2.index    = tag2->index;
            } else {
                intersectInfo->tag1.domainID = tag2->domainID;
                intersectInfo->tag1.index    = tag2->index;
                intersectInfo->tag2.domainID = tag1->domainID;
                intersectInfo->tag2.index    = tag1->index;
            }

            home->segpartIntersectCnt = numEntries + 1;
            home->segpartIntersectListSize = allocedEntries;
            home->segpartIntersectList = tmpList;

            *entry = intersectInfo;
        }

/*
 *      We have a pointer to the segment/particle intersection info,
 *      so just update it now
 */
        nextIndex = (*entry)->numIntersections;

        if (nextIndex < MAX_SEGPART_INTERSECTIONS)
	  {
            (*entry)->inclusionIndex[nextIndex] = inclusionIndex;
            (*entry)->inclusionID[nextIndex] = inclusionID;
            (*entry)->numIntersections += 1;
	  }
	else 
	  {
            for (int i = 0; i < nextIndex; i++)
	      {
		printf("%d %d\n",inclusionID,inclusionIndex);
	    	    
		for (int j = 0; j < (*entry)->numIntersections; j++) 
		  printf("%d %d %d\n",(*entry)->numIntersections,
			 (*entry)->inclusionIndex[j],(*entry)->inclusionID[j]);
		
		Node_t *nodetmp1;
		nodetmp1 = GetNodeFromTag(home, *tag1);
		PrintNode(nodetmp1);
		Node_t *nodetmp2;
		nodetmp2 = GetNodeFromTag(home, *tag2);
		PrintNode(nodetmp2);
	      }

	    Fatal("Segment (%d,%d)--(%d,%d) intersects"
                   " more than %d inclusions\n", tag1->domainID, tag1->index,
                   tag2->domainID, tag2->index, MAX_SEGPART_INTERSECTIONS);
	  }

        UNLOCK(&home->segpartIntersectListLock);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SegPartListClear()
 *      Description:  Release any memory associated with the segment/
 *                    particle intersection list and reinitialize
 *                    the associated variable
 *
 *-------------------------------------------------------------------------*/
void SegPartListClear(Home_t *home)
{
        if (home->segpartIntersectListSize > 0) {
            free(home->segpartIntersectList);
        }

        home->segpartIntersectListSize = 0;
        home->segpartIntersectCnt = 0;
        home->segpartIntersectList = (SegPartIntersect_t *)NULL;

        return;
}



/*---------------------------------------------------------------------------
 *
 *      Function:     CollisionDislocPart()
 *      Description:  Locally remesh a dislocation segment when a particle
 *                    is in the vicinity of the dislocation
 *
 *-------------------------------------------------------------------------*/
void HandleParticleCollisions(Home_t *home)
{
   
   Node_t  *node1, *node2;
                
   int thisDomain;
   thisDomain = home->myDomain;
   
   Param_t *param;
   param      = home->param;

   real8 increaseRadius = 2;

/*
 * Start looping through native nodes looking for segments to collide with particles
 */
   for (int i = 0; i < home->newNodeKeyPtr; i++) 
   {
      if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
      
/*
 *    Loop with all arms of node 1
 *    Only if node 1 is a discretization node
 */
      int nbrs = node1->numNbrs;

/*
 *    Add some restrictions on the segments we'll consider.
 *
 *    Discretization nodes that are neither pinned nor surface nodes.
 */
      if ((nbrs != 2) || 
          HAS_ANY_OF_CONSTRAINTS(node1->constraint,
                                 PINNED_NODE | SURFACE_NODE)) {
         continue;
      }

      
      if ((node1->nbrTag[0].domainID != home->myDomain) ||
          (node1->nbrTag[1].domainID != home->myDomain)) {
         continue;
      }
      
      real8 pos1[3];
      pos1[X] = node1->x;
      pos1[Y] = node1->y;
      pos1[Z] = node1->z;

/*
 *    Start looping over particles now
 */
      for (int m = 0; m < home->totInclusionCount; m++) 
      {
         EInclusion_t *inclusion = &home->eshelbyInclusions[m];

/*  
 *       localMaxSeg value has to be greater than 2*minSeg+eps 
 *       so that the new created segment is not remeshed and
 *       slightly smaller than the diameter of the particle.
 */
         
	 real8 localMaxSeg = MAX(2*param->minSeg*1.01, 1.9*inclusion->radius[X]);
/*
 *       Loop over the node1's arms and determine how many arms to split.
 */
         int nbrSplit = 0, armkeep[2];
	 armkeep[0]=-1;armkeep[1]=-1;
	 
         for (int arm12 = 0; arm12 < nbrs; arm12++)
         {
	   Node_t *node2;
	   node2 = GetNodeFromTag(home, node1->nbrTag[arm12]);
	   if (node2 == (Node_t *)NULL) continue;
	   if (!DomainOwnsSeg(home, OPCLASS_COLLISION, thisDomain, &node2->myTag)) continue;
	   
	   real8 pos2[3];
	   pos2[X] = node2->x;
	   pos2[Y] = node2->y;
	   pos2[Z] = node2->z;
           
	   inclusion->radius[X] *= increaseRadius;
	   inclusion->radius[Y] *= increaseRadius;
	   inclusion->radius[Z] *= increaseRadius;
	   
	   int ninc = 0;
	   real8 newPos1[3], newPos2[3];
	   real8 ratio[2];                  
	   int isIntersecting = IntersectionInformation(home, pos1, pos2, inclusion, 
							newPos1, newPos2, &ninc, ratio);
	   
	   inclusion->radius[X] /= increaseRadius;
	   inclusion->radius[Y] /= increaseRadius;
	   inclusion->radius[Z] /= increaseRadius;
	   
	   real8 segLen = 
	     sqrt( (newPos2[0]-newPos1[0]) * (newPos2[0]-newPos1[0]) +
		   (newPos2[1]-newPos1[1]) * (newPos2[1]-newPos1[1]) +
		   (newPos2[2]-newPos1[2]) * (newPos2[2]-newPos1[2]) );	   
	   
	   if  ( (isIntersecting==1) && (segLen > localMaxSeg))
	     {
               armkeep[nbrSplit] = arm12;
               nbrSplit += 1;
	     }
         } // end loop over arm12
                 

	 for (int isplit=0; isplit<nbrSplit; isplit++) 
	   {
/* 
 *           Split segment  
 */
	     Node_t *node2;
	     node2 = GetNodeFromTag(home, node1->nbrTag[armkeep[isplit]]);
	   
	     real8 dx, dy, dz;
	     dx = node2->x - node1->x;
	     dy = node2->y - node1->y;
	     dz = node2->z - node1->z;

	     ZImage(param, &dx, &dy, &dz);
	   
	     real8 pos2[3];
	     pos2[X] = node1->x + dx;
	     pos2[Y] = node1->y + dy;
	     pos2[Z] = node1->z + dz;
	   
	     real8 newpos[3];
	     newpos[X] = (pos2[X]+pos1[X])*0.5;
	     newpos[Y] = (pos2[Y]+pos1[Y])*0.5;
	     newpos[Z] = (pos2[Z]+pos1[Z])*0.5;
	     
	     real8 nodeVel[3];
	     nodeVel[X] = node2->vX;
	     nodeVel[Y] = node2->vY;
	     nodeVel[Z] = node2->vZ;

	     
	     // Split node 2
	     int splitStatus;
	     Node_t *splitNode1, *splitNode2;
	     int arm21 = GetArmID(node2, node1);

	     splitStatus = SplitNode(home,
				     OPCLASS_COLLISION,
				     node2, pos2,
				     newpos, nodeVel,
				     nodeVel, 1,
				     &arm21, 1,
				     &splitNode1,
				     &splitNode2, 0);

	     if (splitStatus != SPLIT_SUCCESS) continue;

#ifdef FIX_PLASTIC_STRAIN
	     if (splitStatus == SPLIT_SUCCESS) {
	       real8 newPos[3];
	       
	       newPos[X] = splitNode1->x;
	       newPos[Y] = splitNode1->y;
	       newPos[Z] = splitNode1->z;
	       UpdateNodePlasticStrain(home, splitNode1, pos2, newPos);
	       
	       newPos[X] = splitNode2->x;
	       newPos[Y] = splitNode2->y;
	       newPos[Z] = splitNode2->z;
	       UpdateNodePlasticStrain(home, splitNode2, pos2, newPos);
	     }
#endif


	     MarkNodeForceObsolete(home, splitNode2);
	   }	 
	 
      } // end loop over inclusions
   } // loop over all nodes
}

#endif // ESHELBY
