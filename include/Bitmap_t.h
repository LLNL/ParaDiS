#pragma once

#ifndef __PDS_BITMAP_H
#define __PDS_BITMAP_H

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mpi_portability.h"

//----------------------------------------------------------------------------------------

#include "Typedefs.h"

class Bitmap_t
{
   public :
      unsigned int *bmap   ;  ///< contiguous block of bit values
               int  bmap_n ;  ///< bitmap length (in bits)

   public :
       Bitmap_t(void);
       Bitmap_t(const int n);
       Bitmap_t(const Bitmap_t & b);
       Bitmap_t(const unsigned int  p );
       Bitmap_t(const unsigned int *p, const int pn );

      ~Bitmap_t();

      const Bitmap_t & operator =  (const Bitmap_t & b);
      const Bitmap_t & operator =  (const char      *b);

      Bitmap_t operator &  (const Bitmap_t & b);
      Bitmap_t operator |  (const Bitmap_t & b);
      Bitmap_t operator ^  (const Bitmap_t & b);

      void     operator &= (const Bitmap_t & b);
      void     operator |= (const Bitmap_t & b);
      void     operator ^= (const Bitmap_t & b);

      void     Clear     (void);
      void     Set       (void);
      void     Set       (const int j);
      void     Set       (const int *iv, const int icnt);

      int      Bit_Count (void);

      Bitmap_t AND       (const Bitmap_t & b);
      Bitmap_t OR        (const Bitmap_t & b);
      Bitmap_t XOR       (const Bitmap_t & b);

      int      Intersect (const Bitmap_t & b);

      void     Print     (FILE *fd=stdout) const;

};

// BITMAP_M
//
// The individual bits within the bitmap are stored within an unsigned int.
// BITMAP_M will be defined to hold the number of bits within an unsigned int.
//----------------------------------------------------------------------------------------

#define BITMAP_M (8*sizeof(unsigned int))

// inline constructors and destructor....
//----------------------------------------------------------------------------------------

inline Bitmap_t::Bitmap_t(void)               { bmap=0; bmap_n=0; }

inline Bitmap_t::Bitmap_t(const int n)        { bmap   = ( (n>0 ) ? new unsigned int[(n+BITMAP_M-1)/BITMAP_M] : 0 );
                                                bmap_n = ( (bmap) ? n : 0 );
                                                if (bmap) { memset(bmap,0,((n+BITMAP_M-1)/BITMAP_M)*sizeof(unsigned int)); }
                                              }

inline Bitmap_t::Bitmap_t(const Bitmap_t & b) { bmap   = ( (b.bmap) ? new unsigned int[(b.bmap_n+BITMAP_M-1)/BITMAP_M] : 0 );
                                                bmap_n = ( (  bmap) ? b.bmap_n : 0 );
                                                if (bmap && b.bmap) { memcpy(bmap,b.bmap,((bmap_n+BITMAP_M-1)/BITMAP_M)*sizeof(unsigned int)); }
                                              }

inline Bitmap_t::Bitmap_t(const unsigned int  p )               { bmap   = new unsigned int[1];
                                                                  bmap_n = ( (bmap) ? BITMAP_M : 0 );
                                                                  if (bmap) { bmap[0]=p; }
                                                                }

inline Bitmap_t::Bitmap_t(const unsigned int *p, const int pn ) { bmap   = ( (p && (pn>0)) ? new unsigned int[pn] : 0 );
                                                                  bmap_n = ( (bmap) ? pn*BITMAP_M : 0 );
                                                                  if (bmap) { memcpy(bmap,p,pn*sizeof(unsigned int)); }
                                                                }

inline Bitmap_t::~Bitmap_t()                  { if (bmap) { delete [] bmap; } bmap=0; bmap_n=0; }

#endif  // __PDS_BITMAP_H
