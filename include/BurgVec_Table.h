#pragma once

#ifndef _PDS_BURGVEC_TABLE_H
#define _PDS_BURGVEC_TABLE_H

#include "Typedefs.h"
#include "Node.h"
#include "Segment.h"

//------------------------------------------------------------------------------------------------------------
// ParaDiS Burgers Vector Table.
//------------------------------------------------------------------------------------------------------------

class BurgVec_Table
{
   public:
      int     bv_cnt  ;   ///< current   number of burgers vectors in the burgers vector table
      int     bv_max  ;   ///< allocated number of burgers vectors in the burgers vector table
      real8  *bv_tbl  ;   ///< burgers vector table, interleaved <xyz>
      int    *bv_histo;   ///< histogram of entries into the burgers vector table

   public:
      BurgVec_Table(void);
     ~BurgVec_Table();

      void   Append  (const real8  bx, const real8 by, const real8 bz);
      void   Append  (const real8 *bv);
      void   Append  (const Node_t    *nodes, const int ncnt);
      void   Append  (const Segment_t *segs , const int scnt);

      void   Rebuild (const Node_t    *nodes, const int ncnt);
      void   Rebuild (const Segment_t *segs , const int scnt);

      int    Index   (const real8  bx, const real8 by, const real8 bz) const;
      int    Index   (const real8 *bv) const;

      real8 *Burgers_Vector(const int indx) const;

      int   *Histogram (const Segment_t *segs, const int scnt);

      void   Print     (FILE *fd=stdout) const;

};

#endif  // _PDS_BV_TABLE_H
