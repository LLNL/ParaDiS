#pragma once

#ifndef _PDS_NATIVE_SEGMENT_H
#define _PDS_NATIVE_SEGMENT_H

#include "Typedefs.h"
#include "Node.h"
#include "Segment.h"
#include "Cell.h"

//-------------------------------------------------------------------------------------------------------------------

class NativeSeg_t
{
   public :
      Segment_t *seg;
      Cell_t    *cell;
      int       cellID;

   public :
      NativeSeg_t (void)
      {
         seg    = 0;
         cell   = 0;
         cellID = 0;
      }

      NativeSeg_t(const NativeSeg_t & s) 
      { 
         seg    = s.seg;
         cell   = s.cell;
         cellID = s.cellID;
      }

     ~NativeSeg_t() {} // (nothing allocated)

      const NativeSeg_t & operator =  (const NativeSeg_t & s)
      { 
         seg    = s.seg;
         cell   = s.cell;
         cellID = s.cellID;

         return(*this);
      }
};

#endif  // (_PDS_NATIVE_SEGMENT_H)
