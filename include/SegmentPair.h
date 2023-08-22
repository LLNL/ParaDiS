#pragma once

#ifndef _PDS_SEGMENT_PAIR_H
#define _PDS_SEGMENT_PAIR_H

#include "Typedefs.h"
#include "Node.h"
#include "Segment.h"
#include "V3.h"

//-------------------------------------------------------------------------------------------------------------------

class SegmentPair_t
{
   public :
        Segment_t *seg1;
        Segment_t *seg2;
        int        setSeg1Forces;
        int        setSeg2Forces;

   public :
      SegmentPair_t(void)
      {
         seg1          = 0;
         seg2          = 0;
         setSeg1Forces = 0;
         setSeg2Forces = 0;
      }

      SegmentPair_t
      (
         const Segment_t  *s1    ,
         const Segment_t  *s2    ,
         const int         set_s1,
         const int         set_s2
      )
      {
         seg1          = (Segment_t *) s1;
         seg2          = (Segment_t *) s2;
         setSeg1Forces = set_s1;
         setSeg2Forces = set_s2;
      }

      SegmentPair_t(const SegmentPair_t & sp)
      {
         seg1          = sp.seg1;
         seg2          = sp.seg2;
         setSeg1Forces = sp.setSeg1Forces;
         setSeg2Forces = sp.setSeg2Forces;
      }

     ~SegmentPair_t() {} // (nothing allocated)

      const SegmentPair_t & operator =  (const SegmentPair_t & sp)
      {
         seg1          = sp.seg1;
         seg2          = sp.seg2;
         setSeg1Forces = sp.setSeg1Forces;
         setSeg2Forces = sp.setSeg2Forces;

         return(*this);
      }

      void Print (FILE *fd);
      void Print (FILE *fd, const real8 lx, const real8 ly, const real8 lz);
};

//-------------------------------------------------------------------------------------------------------------------

extern SegmentPair_t *SegmentPairs_Append (SegmentPair_t *spairs, int & cnt, int & max, SegmentPair_t & spair);

extern SegmentPair_t *Segment_Pairs_N2    (int & npairs, const Segment_t *segs, const int nsegs);

//-------------------------------------------------------------------------------------------------------------------

#endif  // (_PDS_SEGMENT_PAIR_H)
