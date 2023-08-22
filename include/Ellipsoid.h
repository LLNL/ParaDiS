#pragma once

#ifndef _PDS_ELLIPSOID_H
#define _PDS_ELLIPSOID_H

#include <stdio.h>

#include "Typedefs.h"
#include "Vertex.h"
#include "Triangle.h"

//-----------------------------------------------------------------------------------------------

class Ellipsoid
{
   public :
      real8          pos  [3];
      real8          radii[3];
      real8          rot  [9];
      real8          bbox [6];
  
      Vertex_List    vlist; 
      Triangle_List  tlist; 
   
   public :
      Ellipsoid(void);
      Ellipsoid(const Ellipsoid & e);
     ~Ellipsoid() {}

      real8 *Vertices     (void) const { return(vlist.vtx ); }
      int    Vertex_Cnt   (void) const { return(vlist.vcnt); }
      int   *Triangles    (void) const { return(tlist.tndx); }
      int    Triangle_Cnt (void) const { return(tlist.tcnt); }

      const Ellipsoid & operator = (const Ellipsoid & p);

      void  Tesselate  (const int n);

      void  Scale      (const real8 sx, const real8 sy, const real8 sz);
      void  Transform  (const real8 *m44);

      int   Overlap    (Ellipsoid & e);

      void  Print      (FILE *fd);
      void  Save_GNU   (const char *path);

      void  Display    (                        unsigned long color, double r, unsigned long attr);
      void  Display    (real8 *m33, real8 *pos, unsigned long color, double r, unsigned long attr);
};

#endif // _PDS_ELLIPSOID_H

