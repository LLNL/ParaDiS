#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "Bitmap_t.h"

const Bitmap_t & Bitmap_t::operator =  (const Bitmap_t & b)
{
   if (bmap) { delete [] bmap; }

   bmap_n = b.bmap_n;
   bmap   = ( (bmap_n>0) ? new unsigned int[(bmap_n+BITMAP_M-1)/BITMAP_M] : 0 );

   if (bmap && b.bmap) { memcpy(bmap,b.bmap,((bmap_n+BITMAP_M-1)/BITMAP_M)*sizeof(unsigned int)); }

   return(*this);
}

const Bitmap_t & Bitmap_t::operator =  (const char *b)
{
   if (bmap) { delete [] bmap; }

   bmap_n = 0;
   bmap   = 0;

   if (b)
   {
      char *p    = (char *) b;
      int   pn   = strlen(b);

      if ( (p[0]=='0') && ((p[1]=='b') || (p[1]=='B')) ) { p+=2; pn-=2; }

      bmap_n = ( (pn>0) ? pn : 0 );
      bmap   = ( (pn>0) ? new unsigned int[(pn+BITMAP_M-1)/BITMAP_M] : 0 );

      if (bmap)
      {
         memset(bmap,0x00,((bmap_n+BITMAP_M-1)/BITMAP_M)*sizeof(unsigned int));

         for (int i=0; (i<pn); i++)
         { if ( p[i]=='1') { bmap[i/BITMAP_M] |= ( (unsigned int) 1 ) << (BITMAP_M-(i%BITMAP_M)-1); } }
      }
   }

   return(*this);
}

Bitmap_t Bitmap_t::operator &  (const Bitmap_t & b) { Bitmap_t tmp(*this); tmp&=b; return(tmp); }
Bitmap_t Bitmap_t::operator |  (const Bitmap_t & b) { Bitmap_t tmp(*this); tmp|=b; return(tmp); }
Bitmap_t Bitmap_t::operator ^  (const Bitmap_t & b) { Bitmap_t tmp(*this); tmp^=b; return(tmp); }

Bitmap_t Bitmap_t::AND         (const Bitmap_t & b) { Bitmap_t tmp(*this); tmp&=b; return(tmp); }
Bitmap_t Bitmap_t::OR          (const Bitmap_t & b) { Bitmap_t tmp(*this); tmp|=b; return(tmp); }
Bitmap_t Bitmap_t::XOR         (const Bitmap_t & b) { Bitmap_t tmp(*this); tmp^=b; return(tmp); }

void Bitmap_t::operator &= (const Bitmap_t & b)
{
   int n = ( (bmap_n<b.bmap_n) ? bmap_n : b.bmap_n );
       n = ( (n+BITMAP_M-1)/BITMAP_M );

   unsigned int *ba =   bmap;
   unsigned int *bb = b.bmap;

   if ( ba && bb && (n>0) )
   {
      for (int i=0; (i<n); ++i) { ba[i]|=bb[i]; }
   }
}

void Bitmap_t::operator |= (const Bitmap_t & b)
{
   int n = ( (bmap_n<b.bmap_n) ? bmap_n : b.bmap_n );
       n = ( (n+BITMAP_M-1)/BITMAP_M );

   unsigned int *ba =   bmap;
   unsigned int *bb = b.bmap;

   if ( ba && bb && (n>0) )
   {
      for (int i=0; (i<n); ++i) { ba[i]|=bb[i]; }
   }
}

void Bitmap_t::operator ^= (const Bitmap_t & b)
{
   int n = ( (bmap_n<b.bmap_n) ? bmap_n : b.bmap_n );
       n = ( (n+BITMAP_M-1)/BITMAP_M );

   unsigned int *ba =   bmap;
   unsigned int *bb = b.bmap;

   if ( ba && bb && (n>0) )
   {
      for (int i=0; (i<n); ++i) { ba[i]^=bb[i]; }
   }
}

void Bitmap_t::Clear (void)        { if (bmap) memset(bmap,0x00,((bmap_n+BITMAP_M-1)/BITMAP_M)*sizeof(unsigned int)); }
void Bitmap_t::Set   (void)        { if (bmap) memset(bmap,0xff,((bmap_n+BITMAP_M-1)/BITMAP_M)*sizeof(unsigned int)); }

void Bitmap_t::Set   (const int j) { if (bmap && (0<=j) && (j<bmap_n)) bmap[j/BITMAP_M] |= ( (unsigned int) 1 ) << (BITMAP_M-(j%BITMAP_M)-1); }

void Bitmap_t::Set   (const int *iv, const int icnt)
{
   if (bmap) { memset(bmap,0x00,((bmap_n+BITMAP_M-1)/BITMAP_M)*sizeof(unsigned int)); }

   if (bmap && iv)
   {
      for (int i=0; (i<icnt); i++)
      {
         int j=iv[i];
         if ( (0<=j) && (j<bmap_n)) bmap[j/BITMAP_M] |= ( (unsigned int) 1 ) << (BITMAP_M-(j%BITMAP_M)-1);
      }
   }
}

int Bitmap_t::Intersect (const Bitmap_t & b)
{
   if (bmap && b.bmap && (bmap_n==b.bmap_n) )
   {
      int n = (bmap_n+BITMAP_M-1)/BITMAP_M;
      for (int i=0; (i<n); i++) { if (bmap[i]&&b.bmap[i]) return(1); }
   }
   return(0);
}

int Bitmap_t::Bit_Count (void)
{
   int cnt=0;

   if (bmap)
   {
      for (int i=0; (i<bmap_n); i++)
         if ( bmap[i/BITMAP_M] & ( ( (unsigned int) 1 ) << (BITMAP_M-(i%BITMAP_M)-1) ) )  cnt++;
   }

   return(cnt);
}

void Bitmap_t::Print (FILE *fd) const
{
   if (bmap)
   {
      int m=(bmap_n+BITMAP_M-1)/BITMAP_M;
      for (int i=0; (i<m); i++) { fprintf(fd,"%08x ", bmap[i]); }
   }
}

