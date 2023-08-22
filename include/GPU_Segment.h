#pragma once

#ifndef _PDS_GPU_SEGMENT_H
#define _PDS_GPU_SEGMENT_H

#include <stdlib.h>

#include "cuda_portability.h"

#include "Typedefs.h"
#include "Home.h"
#include "Node.h"
#include "Segment.h"

//------------------------------------------------------------------------------------------------------------

class GPU_Segment
{
   public :
      unsigned int     sndx  ;  ///< segment index

      unsigned long    n1_tag;  ///< node 1 tag 
      unsigned long    n2_tag;  ///< node 2 tag 

      real8            p1[3] ;  ///< node 1 position <xyz> (corrected for pbc)
      real8            p2[3] ;  ///< node 2 position <xyz> (corrected for pbc)
      real8            bv[3] ;  ///< burgers vector (node 1 -> node 2)
      real8            nv[3] ;  ///< glide plane    (node 1 -> node 2)
      real8            f1[3] ;  ///< accumulated force on node 1
      real8            f2[3] ;  ///< accumulated force on node 2

      GPU_Node        *n1    ;  ///< points to gpu node 1 (gpu)
      GPU_Node        *n2    ;  ///< points to gpu node 2 (gpu)

   public :
      __cuda_host__    GPU_Segment(void);
      __cuda_host__    GPU_Segment(const GPU_Segment & s);
      __cuda_host__   ~GPU_Segment() {}

      __cuda_host__    const GPU_Segment & operator =   (const GPU_Segment & s);
      __cuda_host__    const GPU_Segment & operator =   (const Segment_t   & s);

      __cuda_host__    int                 operator ==  (const Segment_t   & s) const;

      __cuda_hdev__    size_t           Bytes           (void) const;
      __cuda_hdev__    void             Reset_Forces    (void);

      __cuda_host__    unsigned char   *Serialize       (unsigned char *p) const;
      __cuda_host__    unsigned char   *Deserialize     (unsigned char *p);
   
      __cuda_host__    size_t           Pack            (unsigned char *pc, unsigned char *pg) const;

      __cuda_host__    void             Print           (FILE *fd=stdout) const;
};

//------------------------------------------------------------------------------------------------------------

__cuda_host__
size_t GPU_Segments_Bytes
(
   Segment_t  *segs,  ///< array of segments
   int         scnt   ///< number of segments
);

//------------------------------------------------------------------------------------------------------------

__cuda_host__
GPU_Segment *GPU_Segments
(
   Segment_t  *segs,  ///< array of segments
   int         scnt   ///< number of segments
);

//------------------------------------------------------------------------------------------------------------

__cuda_host__
size_t GPU_Segment_Pack
(
   unsigned char *pc    ,  ///< points to packing buffer allocated on host
   unsigned char *pg    ,  ///< points to packing buffer allocated on device (gpu)
   Segment_t     *seg      ///< points to segment to be packed (host)
);

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void GPU_Segments_Pack
(
   unsigned char *pc    ,  ///< points to packing buffer allocated on host
   unsigned char *pg    ,  ///< points to packing buffer allocated on device (gpu)
   size_t         pbytes,  ///< overall size of the packing buffers
   Segment_t     *segs  ,  ///< array  of segments (host)
   int            scnt     ///< number of segments
);

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void GPU_Segments_Unpack
(
   Segment_t     *segs  ,  ///< array  of segments (host)
   int            scnt  ,  ///< number of segments
   unsigned char *pc       ///< points to a buffer of packed gpu segments (host)
);

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void GPU_Segments_Reset
(
         GPU_Segment *segs, ///< array of gpu segments (gpu)
   const int          scnt  ///< number of segments
);

//------------------------------------------------------------------------------------------------------------

#endif   // _PDS_GPU_SEGMENT_H
