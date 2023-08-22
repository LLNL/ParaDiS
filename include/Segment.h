#pragma once

#ifndef _PDS_SEGMENT_H
#define _PDS_SEGMENT_H

#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef PTHREADS
#include <pthread.h>
#endif

#include "Typedefs.h"
#include "Node.h"
#include "PBC.h"

//-------------------------------------------------------------------------------------------------------------------

class Segment_t
{
   public:
      Node_t           *node1                 ;  ///< points to node 1
      Node_t           *node2                 ;  ///< points to node 2
      int               forcesSet             ;  ///< flag forces set (1=yes,0=no)
      int               sendParticleIntersects;  ///< flag send particle intersections (1=yes,0=no)
      real8             p1[3]                 ;  ///< position of node 1 <xyz> (corrected for pbc)
      real8             p2[3]                 ;  ///< position of node 2 <xyz> (corrected for pbc)
      real8             f1[3]                 ;  ///< force  on node 1
      real8             f2[3]                 ;  ///< force  on node 2
      real8             bv[3]                 ;  ///< burgers vector (node 1 -> node 2)
      real8             nv[3]                 ;  ///< glide plane    (node 1 -> node 2)
      int               cell_indx             ;  ///< encoded cell index
      int               native                ;  ///< is segment native (1) or ghost (0)

#ifdef _OPENMP
      omp_lock_t        segLock               ;  ///< mutex lock for openmp threads
#endif

#ifdef PTHREADS
      pthread_mutex_t   seg_plock             ;  ///< mutex lock for pthread threads
#endif

   public:
      Segment_t(void);
      Segment_t(const Segment_t & s);
      Segment_t(const Node_t *n1, const Node_t *n2);
      Segment_t(const Node_t *n1, const Node_t *n2, const real8 *b12, const real8 *n12);
      Segment_t(const Node_t *n1, const Node_t *n2, const real8 bx, const real8 by, const real8 bz,
                                                    const real8 nx, const real8 ny, const real8 nz );

     ~Segment_t();

      const Segment_t & operator =  (const Segment_t & s);

      real8 Segment_Length(void) const;

      void Reset_Forces (void);
      void Set_Forces   (const real8 *f1, const real8 *f2);
      void Add_Forces   (const real8 *f1, const real8 *f2);

      int  Parallel     (const Segment_t & s, const real8 ecrit) const;

      void PBC_Position (const real8 lx, const real8 ly, const real8 lz);

      void PBC_Position (const real8 lx, const real8 ly, const real8 lz,
                         const real8 sx, const real8 sy, const real8 sz);

      int  Set_Cell     (const int   nx, const int   ny, const int   nz,
                         const real8 lx, const real8 ly, const real8 lz );

      void Set_Native   (const int dom);
      void Set_Native   (const real8 x0, const real8 x1,
                         const real8 y0, const real8 y1,
                         const real8 z0, const real8 z1 );

      void Print        (FILE *fd=stdout) const;
      void Print        (FILE *fd, const real8 lx, const real8 ly, const real8 lz) const;
};

//-------------------------------------------------------------------------------------------------------------------

extern Segment_t       *Segment_Append             (Segment_t  *segs, int & scnt, int & smax, Segment_t & seg);

extern Segment_t       *Segment_List_Create        (int & scnt, const Node_t **nodes, const int ncnt);

extern Segment_t       *Segment_List_Create        (int & scnt,
                                                    const Node_t **native, const int native_cnt,
                                                    const Node_t **ghosts, const int ghost_cnt  );

extern Segment_t       *Segment_List_Crop          (Segment_t  *segs, const int scnt, int & nscnt,
                                                    const real8  x0, const real8 y0, const real8 z0,
                                                    const real8  x1, const real8 y1, const real8 z1,
                                                    const int    nx, const int   ny, const int   nz,
                                                    const real8  lx, const real8 ly, const real8 lz, real8 *msec=0 );

extern void             Segment_List_Set_Cells     (Segment_t  *segs, const int   scnt,
                                                    const int   nx  , const int   ny  , const int   nz  ,
                                                    const real8 lx  , const real8 ly  , const real8 lz  , real8 *msec=0 );

extern void             Segment_List_Set_Native    (Segment_t  *segs, const int   scnt,
                                                    const real8  x0, const real8 y0, const real8 z0,
                                                    const real8  x1, const real8 y1, const real8 z1, real8 *msec=0 );

extern int              Segment_List_Cnt           (const Segment_t  *segs, const int   scnt, const real8 *bv  );
extern real8            Segment_List_Sum           (const Segment_t  *segs, const int   scnt, const real8 *bv=0);
extern real8            Segment_List_Density       (const Segment_t  *segs, const int   scnt, const real8 *bv  , const real8 bmag, const real8 bvol);

extern void             Segment_List_Moment        (const Segment_t  *segs, const int   scnt, real8 & min, real8 & max, real8 & mean, real8 & sum );
extern void             Segment_List_Moment        (const Segment_t  *segs, const int   scnt, real8 & min, real8 & max, real8 & mean, real8 & sum, real8 & sdev);

extern void             Segment_List_Save          (const char *path, const Segment_t *segs, const int scnt, const real8 lx, const real8 ly, const real8 lz);

extern unsigned short  *Segment_Pairs              (int & npairs,
                                                    Segment_t *segs, const int scnt,
                                                    const real8  x0, const real8 y0, const real8 z0,
                                                    const real8  x1, const real8 y1, const real8 z1,
                                                    const int    nx, const int   ny, const int   nz,
                                                    const real8  lx, const real8 ly, const real8 lz, real8 *msec=0 );

extern void             Segment_Pairs_Save         (const char *path ,
                                                    const unsigned short      *pairs, const int   pcnt,
                                                    const          Segment_t  *segs , const int   scnt,
                                                    const          real8       lx   , const real8 ly  , const real8 lz );

extern int              Segment_List_Check         (Segment_t     *segs, const int   scnt,
                                                    const real8    lx  , const real8 ly  , const real8 lz ,
                                                    const real8    smin,
                                                    const real8    smax,
                                                    const char    *msg =0 );

extern int              Segment_List_Check         (const Node_t **native, const int native_cnt,
                                                    const Node_t **ghosts, const int ghost_cnt ,
                                                    const real8    lx, const real8 ly, const real8 lz,
                                                    const int      nx, const int   ny, const int   nz,
                                                          real8    smin=0.0,
                                                          real8    smax=0.0,
                                                    const char    *msg =0  );

//-------------------------------------------------------------------------------------------------------------------

#endif
