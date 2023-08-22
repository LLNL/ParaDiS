#pragma once

#ifndef _PDS_TRIANGLE_H
#define _PDS_TRIANGLE_H

/*--------------------------------------------------------------------------
 *	V3.h Declarations for 1x3 vector class
 *------------------------------------------------------------------------*/

#include <math.h>
#include "Typedefs.h"

//-----------------------------------------------------------------------------------------------

class Triangle_List
{
   public :
      int    tcnt;
      int    tmax;
      int   *tndx;

   public :
      Triangle_List(void);
      Triangle_List(const int vn);
      Triangle_List(const Triangle_List & t);

     ~Triangle_List();

      int *Triangles   (void) const { return(tndx); }
      int  Triangle_Cnt(void) const { return(tcnt); }

      const Triangle_List & operator = (const Triangle_List & t);

      void Recycle(void);

      void Append(const int t0, const int t1, const int t2);
};

#endif   // (_PDS_TRIANGLE_H)
