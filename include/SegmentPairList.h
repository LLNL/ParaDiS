#pragma once

#ifndef _PDS_SEGMENT_PAIR_LIST_H
#define _PDS_SEGMENT_PAIR_LIST_H

#include <stdio.h>

#include "SegmentPair.h"

//-------------------------------------------------------------------------------------------------------------------

class SegmentPairList_t
{
   public:
      int            spl_cnt  ;   ///< current length of the segment pair list
      int            spl_max  ;   ///< maximum length of the segment pair list
      SegmentPair_t *spl_pairs;   ///< array of segment pairs

   public:
      SegmentPairList_t(void);
     ~SegmentPairList_t();

      void Append   (const SegmentPair_t & spair);

      void Append   (const Segment_t *s1    ,
                     const Segment_t *s2    ,
                     const int        set_s1,
                     const int        set_s2 );

#ifdef SUBCYCLING
      void Append   (const Segment_t *s1    ,
                     const Segment_t *s2    ,
                     const int        set_s1,
                     const int        set_s2,
                     const real8      freqDT,
                     const real8      nextDT );
#endif

      void Print    (FILE *fd);

      void Save     (      FILE *fd  , Home_t *home);
      void Save     (const char *path, Home_t *home);
};

#endif
