#pragma once

#ifndef _PDS_VERTEX_H
#define _PDS_VERTEX_H

/*--------------------------------------------------------------------------
 *	V3.h Declarations for 1x3 vector class
 *------------------------------------------------------------------------*/

#include <math.h>
#include "Typedefs.h"

//-----------------------------------------------------------------------------------------------

class Vertex_List
{
   public : 
      int    vcnt;
      int    vmax;
      real8 *vtx ;

   public :
      Vertex_List(void);
      Vertex_List(const int vn);
      Vertex_List(const Vertex_List & v);

     ~Vertex_List();

      real8 *Vertices  (void) const { return(vtx ); }
      int    Vertex_Cnt(void) const { return(vcnt); }

      const Vertex_List & operator =  (const Vertex_List & v);

      void Recycle   (void);
      void Append    (const real8 vx, const real8 vy, const real8 vz);

      void Scale     (const real8 sx, const real8 sy, const real8 sz);
      void Transform (const real8 *m44);
      void Normalize (void);

      void Get_BBox  (real8 *bbox);
};

#endif  // (_PDS_VERTEX_H)
