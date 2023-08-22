#include <stdio.h>
#include <string.h>

#include "Typedefs.h"
#include "Triangle.h"

Triangle_List::Triangle_List(void)
{
   tcnt = 0;
   tmax = 0;
   tndx = 0;
}

Triangle_List::Triangle_List(const int tn)
{
   tndx = ( (tn>0) ? new int[3*tn] : 0 );
   tmax = ( tndx ? tn : 0 );
   tcnt = 0;
}

Triangle_List::Triangle_List(const Triangle_List & t)
{
   tcnt = 0;
   tmax = 0;
   tndx = 0;

   if (t.tndx && (t.tcnt>0))
   {
      tcnt = t.tcnt;
      tmax = t.tcnt;
      tndx = new int[3*t.tcnt];

      if (tndx) { memcpy(tndx,t.tndx,3*tcnt*sizeof(int)); }
   }
}

Triangle_List::~Triangle_List()
{
   if (tndx) { delete [] tndx; tndx=0; }

   tcnt=0;
   tmax=0;
   tndx=0;
}

void Triangle_List::Recycle(void)
{
   if (tndx) { delete [] tndx; tndx=0; }

   tcnt=0;
   tmax=0;
   tndx=0;
}

const Triangle_List & Triangle_List::operator = (const Triangle_List & t)
{
   if (this!=&t)
   {
      Recycle();
   
      if (t.tndx && (t.tcnt>0))
      {
         tndx = new int[3*t.tcnt];
         tmax = t.tcnt;
         tcnt = t.tcnt;
   
         if (tndx) { memcpy(tndx,t.tndx,3*tcnt*sizeof(int)); }
      }
   }

   return(*this);
}

void Triangle_List::Append(const int t0, const int t1, const int t2)
{
   if ( !tndx || (tcnt==tmax) )
   {
      tmax += 200;

      int *tmp = new int[3*tmax];

      if (tmp && tndx) { memcpy(tmp,tndx,3*tcnt*sizeof(int)); }

      if (tndx) { delete [] tndx; }

      tndx = tmp;
   }

   if (tndx && (tcnt<tmax) ) 
   { 
      tndx[3*tcnt  ] = t0; 
      tndx[3*tcnt+1] = t1; 
      tndx[3*tcnt+2] = t2;  tcnt++;
   }
}

