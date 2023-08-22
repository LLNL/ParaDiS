#include <stdio.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "SegmentPair.h"

//-------------------------------------------------------------------------------------------------------------------
// Note - constructors and class inlined in header
//-------------------------------------------------------------------------------------------------------------------

static void PBC_Adjust
(
   const real8 *p1,   ///< position 1 (unchanged)
         real8 *p2,   ///< position 2 (updated)
   const real8  lx,   ///< simulation size <x>
   const real8  ly,   ///< simulation size <y>
   const real8  lz    ///< simulation size <z>
)
{
   const real8 lx2 = (lx/2.0);
   const real8 ly2 = (ly/2.0);
   const real8 lz2 = (lz/2.0);

   const real8 dx  = (p2[0]-p1[0]);
   const real8 dy  = (p2[1]-p1[1]);
   const real8 dz  = (p2[2]-p1[2]);

   p2[0] += ( (dx< -lx2) ? lx : ( (dx>lx2) ?  -lx : 0.0 ) );
   p2[1] += ( (dy< -ly2) ? ly : ( (dy>ly2) ?  -ly : 0.0 ) );
   p2[2] += ( (dz< -lz2) ? lz : ( (dz>lz2) ?  -lz : 0.0 ) );
}

// Print()
//
// Print diagnostic for the segment pair class
//-------------------------------------------------------------------------------------------------------------------

void SegmentPair_t::Print (FILE *fd)
{
   Node_t  *n1 = ( seg1 ? seg1->node1 : 0 );
   Node_t  *n2 = ( seg1 ? seg1->node2 : 0 );
   Node_t  *n3 = ( seg2 ? seg2->node1 : 0 );
   Node_t  *n4 = ( seg2 ? seg2->node2 : 0 );

   if (fd && n1 && n2 && n3 &&n4 )
   {
      char  t1[32]; sprintf(t1,"(%d,%d)", n1->myTag.domainID, n1->myTag.index );
      char  t2[32]; sprintf(t2,"(%d,%d)", n2->myTag.domainID, n2->myTag.index );
      char  t3[32]; sprintf(t3,"(%d,%d)", n3->myTag.domainID, n3->myTag.index );
      char  t4[32]; sprintf(t4,"(%d,%d)", n4->myTag.domainID, n4->myTag.index );

      real8 *bv1 = seg1->bv;
      real8 *bv2 = seg2->bv;
      real8 *nv1 = seg1->nv;
      real8 *nv2 = seg2->nv;

      real8  p1[3] = { n1->x, n1->y, n1->z };
      real8  p2[3] = { n2->x, n2->y, n2->z };
      real8  p3[3] = { n3->x, n3->y, n3->z };
      real8  p4[3] = { n4->x, n4->y, n4->z };

      fprintf(fd,"%12s %12s %12s %12s"     , t1, t2, t3, t4 );
      fprintf(fd," %12.6lf %12.6lf %12.6lf", p1 [0], p1 [1], p1 [2] );
      fprintf(fd," %12.6lf %12.6lf %12.6lf", p2 [0], p2 [1], p2 [2] );
      fprintf(fd," %12.6lf %12.6lf %12.6lf", p3 [0], p3 [1], p3 [2] );
      fprintf(fd," %12.6lf %12.6lf %12.6lf", p4 [0], p4 [1], p4 [2] );
      fprintf(fd," %12.8lf %12.8lf %12.8lf", bv1[0], bv1[1], bv1[2] );
      fprintf(fd," %12.8lf %12.8lf %12.8lf", nv1[0], nv1[1], nv1[2] );
      fprintf(fd," %12.8lf %12.8lf %12.8lf", bv2[0], bv2[1], bv2[2] );
      fprintf(fd," %12.8lf %12.8lf %12.8lf", nv2[0], nv2[1], nv2[2] );
      fprintf(fd,"\n");
   }
}

void SegmentPair_t::Print (FILE *fd, const real8 lx, const real8 ly, const real8 lz)
{
   Node_t  *n1 = ( seg1 ? seg1->node1 : 0 );
   Node_t  *n2 = ( seg1 ? seg1->node2 : 0 );
   Node_t  *n3 = ( seg2 ? seg2->node1 : 0 );
   Node_t  *n4 = ( seg2 ? seg2->node2 : 0 );

   if (fd && n1 && n2 && n3 &&n4 )
   {
      char  t1[32]; sprintf(t1,"(%d,%d)", n1->myTag.domainID, n1->myTag.index );
      char  t2[32]; sprintf(t2,"(%d,%d)", n2->myTag.domainID, n2->myTag.index );
      char  t3[32]; sprintf(t3,"(%d,%d)", n3->myTag.domainID, n3->myTag.index );
      char  t4[32]; sprintf(t4,"(%d,%d)", n4->myTag.domainID, n4->myTag.index );

      real8 *bv1 = seg1->bv;
      real8 *bv2 = seg2->bv;
      real8 *nv1 = seg1->nv;
      real8 *nv2 = seg2->nv;

      real8  p1[3] = { n1->x, n1->y, n1->z };
      real8  p2[3] = { n2->x, n2->y, n2->z };
      real8  p3[3] = { n3->x, n3->y, n3->z };
      real8  p4[3] = { n4->x, n4->y, n4->z };

      PBC_Adjust(p1,p2,lx,ly,lz);
      PBC_Adjust(p3,p4,lx,ly,lz);

      fprintf(fd,"%12s %12s %12s %12s"     , t1, t2, t3, t4 );
      fprintf(fd," %12.6lf %12.6lf %12.6lf", p1 [0], p1 [1], p1 [2] );
      fprintf(fd," %12.6lf %12.6lf %12.6lf", p2 [0], p2 [1], p2 [2] );
      fprintf(fd," %12.6lf %12.6lf %12.6lf", p3 [0], p3 [1], p3 [2] );
      fprintf(fd," %12.6lf %12.6lf %12.6lf", p4 [0], p4 [1], p4 [2] );
      fprintf(fd," %12.8lf %12.8lf %12.8lf", bv1[0], bv1[1], bv1[2] );
      fprintf(fd," %12.8lf %12.8lf %12.8lf", nv1[0], nv1[1], nv1[2] );
      fprintf(fd," %12.8lf %12.8lf %12.8lf", bv2[0], bv2[1], bv2[2] );
      fprintf(fd," %12.8lf %12.8lf %12.8lf", nv2[0], nv2[1], nv2[2] );
      fprintf(fd,"\n");
   }
}

//-------------------------------------------------------------------------------------------------------------------

static void Malloc_Check(void *p, const char *file, const char *func, const int ln)
{
   if (!p)
   {
      printf("%s::%s(ln=%d) memory allocation error\n", file, func, ln );
      MPI_Abort(MPI_COMM_WORLD,-1);
      exit(0);
   }
}

#define MALLOC_CHECK(p) Malloc_Check((p), __FILE__, __func__, __LINE__);

//-------------------------------------------------------------------------------------------------------------------
// SegmentPair_Append()
//
// Will append a segment pair to an existing segment pair array.  Reallocates and extends the 
// array as necessary. This could be more efficient by implementing the lists as linked lists, but 
// too much of the code base assumes the ability to iterate via array elements.
//
// Note - this routine is intended to be called via an assignment.
// 
// Example...
//    SegmentPair_t *segs=0; int scnt=0, smax=0;
//
//    for (int i=0; (i<...); i++)
//       segs = SegmentPair_Append(segs,scnt,smax, SegmentPair_t(...) );
//
//-------------------------------------------------------------------------------------------------------------------

SegmentPair_t *SegmentPairs_Append (SegmentPair_t *spairs, int & cnt, int & max, SegmentPair_t & spair)
{
   if ( !spairs || (cnt==max) )
   {
      max += ( (max==0) ? (10*1024) : (max/2) );

      SegmentPair_t *tmp = ( (max>0) ? new SegmentPair_t[max] : 0 );

      MALLOC_CHECK(tmp)

      if (tmp && spairs)
         for (int i=0; (i<cnt); i++) { tmp[i]=spairs[i]; }

      if (spairs) { delete [] spairs; }

      spairs = tmp;
   }

   if (spairs && (cnt<max) ) { spairs[cnt++]=spair; }

   return(spairs);
}

// Segment_Pairs_N2()
//
// Will construct and return a pure N-squared segment pair list from a source segment list.
// This routine was mainly provided for unit test support.
//-------------------------------------------------------------------------------------------------------------------

SegmentPair_t *Segment_Pairs_N2 
(
         int       & npairs,  ///< number of segment pairs (returned)
   const Segment_t  *segs  ,  ///< array of source segments
   const int         nsegs    ///< number of source segments
)
{
   int nmax = (nsegs*nsegs)/2;   // (assume (n*n)/2 segment pairs)

   SegmentPair_t *spairs = ( (nmax>0) ? new SegmentPair_t[nmax] : 0 );
                  npairs = 0;

   MALLOC_CHECK(spairs)

   if (segs && (nsegs>0))
   {
      for (int i=0;     (i<nsegs); i++ )
      for (int j=(i+1); (j<nsegs); j++ )
      {
         SegmentPair_t spair((segs+i),(segs+j),1,1);
      
         spairs = SegmentPairs_Append(spairs,npairs,nmax,spair);
      }
   }

   return(spairs);
}

