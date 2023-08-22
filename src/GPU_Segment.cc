#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "cuda_portability.h"

#include "Typedefs.h"
#include "Home.h"
#include "Node.h"
#include "Segment.h"
#include "GPU_Node.h"
#include "GPU_Segment.h"
#include "V3.h"

//------------------------------------------------------------------------------------------------------------

#define VDELETE(a) if (a) { delete [] a; a=0; }

//------------------------------------------------------------------------------------------------------------

__cuda_host__
GPU_Segment::GPU_Segment(void)
{
   sndx   = 0;

   n1_tag = 0;
   n2_tag = 0;

   V3_ZERO(p1);
   V3_ZERO(p2);
   V3_ZERO(bv);
   V3_ZERO(nv);
   V3_ZERO(f1);
   V3_ZERO(f2);

   n1 = (GPU_Node *) 0;
   n2 = (GPU_Node *) 0;
}

__cuda_host__
GPU_Segment::GPU_Segment(const GPU_Segment & s)
{
   sndx   = s.sndx;

   n1_tag = s.n1_tag;
   n2_tag = s.n2_tag;

   V3_COPY(p1,s.p1);
   V3_COPY(p2,s.p2);
   V3_COPY(bv,s.bv);
   V3_COPY(nv,s.nv);
   V3_COPY(f1,s.f1);
   V3_COPY(f2,s.f2);

   n1 = s.n1;
   n2 = s.n2;
}

__cuda_host__
const GPU_Segment & GPU_Segment::operator = (const GPU_Segment & s)
{
   if (this!=&s)
   {
      sndx   = s.sndx;
   
      n1_tag = s.n1_tag;
      n2_tag = s.n2_tag;
   
      V3_COPY(p1,s.p1);
      V3_COPY(p2,s.p2);
      V3_COPY(bv,s.bv);
      V3_COPY(nv,s.nv);
      V3_COPY(f1,s.f1);
      V3_COPY(f2,s.f2);
   
      n1 = s.n1;
      n2 = s.n2;
   }

   return(*this);
}

__cuda_host__
const GPU_Segment & GPU_Segment::operator = (const Segment_t & s)
{
   sndx   = 0;

   n1_tag = ( s.node1 ? (unsigned long) s.node1->myTag : 0 );
   n2_tag = ( s.node2 ? (unsigned long) s.node2->myTag : 0 );

   V3_COPY(p1,s.p1);
   V3_COPY(p2,s.p2);
   V3_COPY(bv,s.bv);
   V3_COPY(nv,s.nv);
   V3_COPY(f1,s.f1);
   V3_COPY(f2,s.f2);

   n1 = (GPU_Node *) 0;
   n2 = (GPU_Node *) 0;

   return(*this);
}

__cuda_host__
int GPU_Segment::operator == (const Segment_t & s) const
{
   if ( n1_tag != ( s.node1 ? (unsigned long) s.node1->myTag : 0 ) ) return(0);
   if ( n2_tag != ( s.node2 ? (unsigned long) s.node2->myTag : 0 ) ) return(0);

   for (int i=0; (i<3); i++) { if ( fabs(p1[i]-s.p1[i]) > 1.0e-4 ) return(0); }
   for (int i=0; (i<3); i++) { if ( fabs(p2[i]-s.p2[i]) > 1.0e-4 ) return(0); }
   for (int i=0; (i<3); i++) { if ( fabs(bv[i]-s.bv[i]) > 1.0e-6 ) return(0); }
   for (int i=0; (i<3); i++) { if ( fabs(nv[i]-s.nv[i]) > 1.0e-6 ) return(0); }

   for (int i=0; (i<3); i++) { if ( fabs(f1[i]-s.f1[i]) > 1.0e-2 ) return(0); }
   for (int i=0; (i<3); i++) { if ( fabs(f2[i]-s.f2[i]) > 1.0e-2 ) return(0); }

   return(1);
}

size_t GPU_Segment::Bytes (void) const
{
   return(sizeof(GPU_Segment));
}

__cuda_hdev__
void GPU_Segment::Reset_Forces (void)
{
   V3_ZERO(f1);
   V3_ZERO(f2);
}

__cuda_host__
unsigned char   *GPU_Segment::Serialize (unsigned char *p) const
{
   if (p)
   {
      GPU_Segment *s = (GPU_Segment *) p;  p+=sizeof(GPU_Segment);

      s->sndx   = sndx;

      s->n1_tag = n1_tag;
      s->n2_tag = n2_tag;

      V3_COPY(s->p1,p1);
      V3_COPY(s->p2,p2);
      V3_COPY(s->bv,bv);
      V3_COPY(s->nv,nv);
      V3_COPY(s->f1,f1);
      V3_COPY(s->f2,f2);

      s->n1 = n1;
      s->n2 = n2;
   }

   return(p);
}

__cuda_host__
unsigned char   *GPU_Segment::Deserialize     (unsigned char *p)
{
   if (p)
   {
      GPU_Segment *s = (GPU_Segment *) p;  p+=sizeof(GPU_Segment);

      sndx   = s->sndx;

      n1_tag = s->n1_tag;
      n2_tag = s->n2_tag;

      V3_COPY(p1,s->p1);
      V3_COPY(p2,s->p2);
      V3_COPY(bv,s->bv);
      V3_COPY(nv,s->nv);
      V3_COPY(f1,s->f1);
      V3_COPY(f2,s->f2);

      n1 = s->n1;
      n2 = s->n2;
   }

   return(p);
}

__cuda_host__
size_t GPU_Segment::Pack (unsigned char *pc, unsigned char *pg) const
{

   return(0);
}

__cuda_host__
void GPU_Segment::Print (FILE *fd) const
{
   if (fd)
   {
      fprintf(fd,"seg=%5u "    , sndx );
      fprintf(fd,"t1=(%d,%d) " , (int) (n1_tag>>32), (int) (n1_tag & 0xffffffff) );
      fprintf(fd,"t2=(%d,%d)\n", (int) (n2_tag>>32), (int) (n2_tag & 0xffffffff) );

      fprintf(fd,"p1= %10.6lf %10.6lf %10.6lf\n", p1[0], p1[1], p1[2] );
      fprintf(fd,"p2= %10.6lf %10.6lf %10.6lf\n", p2[0], p2[1], p2[2] );
      fprintf(fd,"bv= %10.6lf %10.6lf %10.6lf\n", bv[0], bv[1], bv[2] );
      fprintf(fd,"nv= %10.6lf %10.6lf %10.6lf\n", nv[0], nv[1], nv[2] );
      fprintf(fd,"f1= %12.8le %12.8le %12.8le\n", f1[0], f1[1], f1[2] );
      fprintf(fd,"f2= %12.8le %12.8le %12.8le\n", f2[0], f2[1], f2[2] );
   }
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__
size_t GPU_Segments_Bytes
(
   Segment_t  *segs,  ///< array of segments
   int         scnt   ///< number of segments
)
{
   return( (segs && (scnt>0)) ? scnt*sizeof(GPU_Segment) : 0 );
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__
GPU_Segment *GPU_Segments
(
   Segment_t  *segs,  ///< array of segments
   int         scnt   ///< number of segments
)
{
   GPU_Segment *segs_gpu = ( (segs && (scnt>0)) ? new GPU_Segment[scnt] : 0 );

   if (segs_gpu)
   {
      for (int i=0; (i<scnt); i++)
      {
         segs_gpu[i]      = segs[i];
         segs_gpu[i].sndx = i;
      }
   }

   return(segs_gpu);
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void GPU_Segments_Reset
(
         GPU_Segment *segs, ///< array of gpu segments (gpu)
   const int          scnt  ///< number of segments
)
{
   if (segs && (scnt>0))
   { for (int i=0; (i<scnt); i++) segs[i].Reset_Forces(); }
}

//------------------------------------------------------------------------------------------------------------

