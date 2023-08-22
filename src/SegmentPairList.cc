#include <stdio.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "Home.h"
#include "Node.h"
#include "Segment.h"
#include "SegmentPair.h"
#include "SegmentPairList.h"
#include "V3.h"

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

SegmentPairList_t::SegmentPairList_t(void)
{
   spl_cnt   = 0;
   spl_max   = 0;
   spl_pairs = 0;
}

SegmentPairList_t::~SegmentPairList_t()
{
   if (spl_pairs) { delete [] spl_pairs; }

   spl_cnt   = 0;
   spl_max   = 0;
   spl_pairs = 0;
}

//-------------------------------------------------------------------------------------------------------------------

void SegmentPairList_t::Append (const SegmentPair_t & spair)
{
   if (!spl_pairs || (spl_cnt==spl_max))
   {
      spl_max += ( (spl_max==0) ? (10*1024) : (spl_max/2) );

      SegmentPair_t *tmp = ( (spl_max>0) ? new SegmentPair_t[spl_max] : 0 );

      MALLOC_CHECK(tmp)

      if (tmp && spl_pairs)
         for (int i=0; (i<spl_cnt); i++) { tmp[i]=spl_pairs[i]; }

      if (spl_pairs) { delete [] spl_pairs; }

      spl_pairs = tmp;
   }

   if (spl_pairs && (spl_cnt<spl_max) ) { spl_pairs[spl_cnt++]=spair; }
}

void SegmentPairList_t::Append
(
   const Segment_t *s1    ,     ///< points to segment 1
   const Segment_t *s2    ,     ///< points to segment 2
   const int        set_s1,     ///< set s1 forces (1=yes,0=no)
   const int        set_s2      ///< set s1 forces (1=yes,0=no)
)
{
   Append( SegmentPair_t(s1,s2,set_s1,set_s2) );
}

//-------------------------------------------------------------------------------------------------------------------

void SegmentPairList_t::Print (FILE *fd)
{
   if (fd && spl_pairs)
   {
      for (int i=0; (i<spl_cnt); i++)
         spl_pairs[i].Print(fd);
   }
}

// Save()
//
// Two forms. Given a path or open file descriptor will save the contents of the segment pair list
// to a binary output file.  Mainly used for debug and unit testing.
//-------------------------------------------------------------------------------------------------------------------

void SegmentPairList_t::Save(const char *path, Home_t *home)
{
   FILE *fd = (FILE *) ( path ? fopen(path,"w") : 0 );

   Save(fd,home);
}

void SegmentPairList_t::Save(FILE *fd, Home_t *home)
{
    if (fd)
    {
        SegmentPair_t *sp = spl_pairs;                               // sp points to source segment pair vector
        real8          q[7*3];                                       // q will be used to accumulate the pair info

        for (int i=0; (i<spl_cnt); i++, sp++ )
        {
            Node_t *n1   = sp->seg1->node1;                          // n1 points to node 1
            Node_t *n2   = sp->seg1->node2;                          // n2 points to node 2
            Node_t *n3   = sp->seg2->node1;                          // n3 points to node 3
            Node_t *n4   = sp->seg2->node2;                          // n4 points to node 4
            Cell_t *cell = LookupCell(home, n1->cellIdx);            // cell points to the cell containing node 1 (primary cell)

            real8  *x1   = q   ; V3_SET (x1, n1->x, n1->y, n1->z);   // copy position of node 1 <xyz>
            real8  *x2   = q+ 3; V3_SET (x2, n2->x, n2->y, n2->z);   // copy position of node 2 <xyz>
            real8  *x3   = q+ 6; V3_SET (x3, n3->x, n3->y, n3->z);   // copy position of node 3 <xyz>
            real8  *x4   = q+ 9; V3_SET (x4, n4->x, n4->y, n4->z);   // copy position of node 4 <xyz>
            real8  *c1   = q+12; V3_COPY(c1, cell->center);          // copy position of cell center <xyz>
            real8  *b12  = q+15; V3_ZERO(b12);                       // initialize segment 1 burgers vector to zero
            real8  *b34  = q+18; V3_ZERO(b34);                       // initialize segment 2 burgers vector to zero

            int     i12  = GetArmID(n1,n2);                          // i12 = arm index of segment 1
            int     i34  = GetArmID(n3,n4);                          // i34 = arm index of segment 2

            if (i12>=0) { V3_SET(b12, n1->burgX[i12], n1->burgY[i12], n1->burgZ[i12]); } // set segment 1 burgers vector
            if (i34>=0) { V3_SET(b34, n3->burgX[i34], n3->burgY[i34], n3->burgZ[i34]); } // set segment 2 burgers vector

            fwrite(q,sizeof(q),1,fd);
        }

        fclose(fd);
    }
}
