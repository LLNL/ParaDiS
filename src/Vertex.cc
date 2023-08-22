#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "Typedefs.h"
#include "V3.h"
#include "M44.h"
#include "Vertex.h"

#define MIN(a,b) ( (a)<(b) ? (a) : (b) )
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )

Vertex_List::Vertex_List(void)
{
   vtx  = 0;
   vmax = 0;
   vcnt = 0;
}

Vertex_List::Vertex_List(const int vn)
{
   vtx  = ( (vn>0) ? new real8[3*vn] : 0 );
   vmax = ( vtx ? vn : 0 );
   vcnt = 0;
}

Vertex_List::Vertex_List(const Vertex_List & v)
{
   vtx  = 0;
   vmax = 0;
   vcnt = 0;

   if (v.vtx && (v.vcnt>0))
   {
      vtx  = new real8[3*v.vcnt];
      vmax = v.vcnt;
      vcnt = v.vcnt;

      if (vtx) { memcpy(vtx,v.vtx,3*vcnt*sizeof(real8)); }
   }
}

Vertex_List::~Vertex_List()
{
   if (vtx) { delete [] vtx; vtx=0; }

   vcnt=0;
   vmax=0;
   vtx =0;
}

const Vertex_List & Vertex_List::operator =  (const Vertex_List & v)
{
   if (this!=&v)
   {
      Recycle();
   
      if (v.vtx && (v.vcnt>0))
      {
         vtx  = new real8[3*v.vcnt];
         vmax = v.vcnt;
         vcnt = v.vcnt;
   
         if (vtx) { memcpy(vtx,v.vtx,3*vcnt*sizeof(real8)); }
      }
   }

   return(*this);
}

void Vertex_List::Recycle(void)
{
   if (vtx) { delete [] vtx; vtx=0; }

   vcnt=0;
   vmax=0;
   vtx =0;
}

void Vertex_List::Append(const real8 vx, const real8 vy, const real8 vz)
{
   if ( !vtx || (vcnt==vmax) )
   {
      vmax += 200;

      real8 *tmp = new real8[3*vmax];

      if (tmp && vtx) { memcpy(tmp,vtx,3*vcnt*sizeof(real8)); }

      if (vtx) { delete [] vtx; }

      vtx = tmp;
   }

   if (vtx && (vcnt<vmax) ) 
   { 
      vtx[3*vcnt  ] = vx; 
      vtx[3*vcnt+1] = vy; 
      vtx[3*vcnt+2] = vz;  vcnt++;
   }
}

void Vertex_List::Scale  (const real8 sx, const real8 sy, const real8 sz)
{
   if (vtx && (vcnt>0))
   {
      real8 *v = vtx;
      for (int i=0; (i<vcnt); i++, v+=3)
      {
         v[0] *= sx;
         v[1] *= sy;
         v[2] *= sz;
      }
   }
}

void Vertex_List::Transform  (const real8 *m44)
{
   if (vtx && m44)
   {
      real8 *v = vtx;
      for (int i=0; (i<vcnt); i++, v+=3)
         V3_M44_V3_MUL(v,m44,v);
   }
}

void Vertex_List::Normalize (void)
{
   if (vtx && (vcnt>0))
   {
      real8 *v = vtx;
      for (int i=0; (i<vcnt); i++, v+=3)
         V3_NORMALIZE(v,v);
   }
}

void Vertex_List::Get_BBox  (real8 *bbox)
{
   if (bbox)
   {
      bbox[0]=0.0;
      bbox[1]=0.0;
      bbox[2]=0.0;
      bbox[3]=0.0;
      bbox[4]=0.0;
      bbox[5]=0.0;
   }

   if (bbox && vtx && (vcnt>0))
   {
      real8 *v = vtx;

      bbox[0] = v[0];
      bbox[1] = v[0];
      bbox[2] = v[1];
      bbox[3] = v[1];
      bbox[4] = v[2];
      bbox[5] = v[2];
      
      for (int i=0; (i<vcnt); i++, v+=3)
      {
         bbox[0] = MIN(bbox[0],v[0]);
         bbox[1] = MAX(bbox[1],v[0]);
         bbox[2] = MIN(bbox[2],v[1]);
         bbox[3] = MAX(bbox[3],v[1]);
         bbox[4] = MIN(bbox[4],v[2]);
         bbox[5] = MAX(bbox[5],v[2]);
      }
   }
}
