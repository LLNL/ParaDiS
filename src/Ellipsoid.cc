#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Ellipsoid.h"
#include "V3.h"
#include "M33.h"
#include "M44.h"
#include "Vertex.h"
#include "Triangle.h"
#include "DisplayC.h"

#define MIN(a,b) ( (a)<(b) ? (a) : (b) )
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )

//-----------------------------------------------------------------------------------------------

static void _init (Ellipsoid & e)
{
   e.pos  [0]=e.pos  [1]=e.pos  [2]=0.0;
   e.radii[0]=e.radii[1]=e.radii[2]=1.0;

   M33_IDENTITY(e.rot);

   e.bbox[0]  = -1.0; e.bbox[1]  = +1.0;
   e.bbox[2]  = -1.0; e.bbox[3]  = +1.0;
   e.bbox[4]  = -1.0; e.bbox[5]  = +1.0;
}

Ellipsoid::Ellipsoid (void)                { _init(*this); }
Ellipsoid::Ellipsoid (const Ellipsoid & e) { _init(*this);  *this=e; }

//-----------------------------------------------------------------------------------------------

const Ellipsoid & Ellipsoid::operator = (const Ellipsoid & e)
{
   if (this!=&e)
   {
      memcpy(pos  ,e.pos  ,sizeof(pos  ));
      memcpy(radii,e.radii,sizeof(radii));
      memcpy(rot  ,e.rot  ,sizeof(rot  ));
      memcpy(bbox ,e.bbox ,sizeof(bbox ));

      vlist = e.vlist; 
      tlist = e.tlist; 
   }

   return(*this);
}

//-----------------------------------------------------------------------------------------------

// Tesselate()
//
// Creates a tesselated ellipsoid.
// source : http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
//-----------------------------------------------------------------------------------------------

void Ellipsoid::Tesselate (const int rn)
{
   vlist.Recycle();
   tlist.Recycle();

   // initialize the base icosahedron vertices...
   {
      real8  t = (1.0+sqrt(5.0))/2.0;

      vlist.Append(-1.0,   t ,  0.0); 
      vlist.Append( 1.0,   t ,  0.0); 
      vlist.Append(-1.0,  -t ,  0.0); 
      vlist.Append( 1.0,  -t ,  0.0); 

      vlist.Append( 0.0, -1.0,   t ); 
      vlist.Append( 0.0,  1.0,   t ); 
      vlist.Append( 0.0, -1.0,  -t ); 
      vlist.Append( 0.0,  1.0,  -t ); 

      vlist.Append(  t ,  0.0, -1.0); 
      vlist.Append(  t ,  0.0,  1.0); 
      vlist.Append( -t ,  0.0, -1.0); 
      vlist.Append( -t ,  0.0,  1.0); 

      vlist.Normalize();
   }

   // initialize the icosohedran triangles...
   {
      // 5 faces around point 0
      tlist.Append( 0,11, 5);
      tlist.Append( 0, 5, 1);
      tlist.Append( 0, 1, 7);
      tlist.Append( 0, 7,10);
      tlist.Append( 0,10,11);
      
      // 5 adjacent faces
      tlist.Append( 1, 5, 9);
      tlist.Append( 5,11, 4);
      tlist.Append(11,10, 2);
      tlist.Append(10, 7, 6);
      tlist.Append( 7, 1, 8);
      
      // 5 faces around point 3
      tlist.Append( 3, 9, 4);
      tlist.Append( 3, 4, 2);
      tlist.Append( 3, 2, 6);
      tlist.Append( 3, 6, 8);
      tlist.Append( 3, 8, 9);
      
      // 5 adjacent faces
      tlist.Append( 4, 9, 5);
      tlist.Append( 2, 4,11);
      tlist.Append( 6, 2,10);
      tlist.Append( 8, 6, 7);
      tlist.Append( 9, 8, 1);
   }
   
   // refine the triangles via recursion...

   for (int i=0; (i<rn); i++)
   {
      Ellipsoid tmp;

      real8 *v0   = vlist.Vertices    ();
      int   *t0   = tlist.Triangles   (); 
      int    tcnt = tlist.Triangle_Cnt();

      // cycle through all the current triangles and subdivide them...

      for (int j=0; (j<(3*tcnt)); j+=3)
      {
         int i0 = 3*t0[j  ];
         int i1 = 3*t0[j+1];
         int i2 = 3*t0[j+2];

         real8 p0[3] = { v0[i0], v0[i0+1], v0[i0+2] };
         real8 p1[3] = { v0[i1], v0[i1+1], v0[i1+2] };
         real8 p2[3] = { v0[i2], v0[i2+1], v0[i2+2] };

         real8 p3[3] = { (p0[0]+p1[0])/2.0,
                         (p0[1]+p1[1])/2.0,
                         (p0[2]+p1[2])/2.0 };

         real8 p4[3] = { (p1[0]+p2[0])/2.0,
                         (p1[1]+p2[1])/2.0,
                         (p1[2]+p2[2])/2.0 };

         real8 p5[3] = { (p2[0]+p0[0])/2.0,
                         (p2[1]+p0[1])/2.0,
                         (p2[2]+p0[2])/2.0 };

         V3_NORMALIZE(p3,p3);
         V3_NORMALIZE(p4,p4);
         V3_NORMALIZE(p5,p5);

         int k=tmp.vlist.vcnt;
         
         tmp.vlist.Append(p0[0],p0[1],p0[2]);
         tmp.vlist.Append(p1[0],p1[1],p1[2]);
         tmp.vlist.Append(p2[0],p2[1],p2[2]);
         tmp.vlist.Append(p3[0],p3[1],p3[2]);
         tmp.vlist.Append(p4[0],p4[1],p4[2]);
         tmp.vlist.Append(p5[0],p5[1],p5[2]);
        
         tmp.tlist.Append( k+0, k+3, k+5 );
         tmp.tlist.Append( k+3, k+1, k+4 );
         tmp.tlist.Append( k+5, k+4, k+2 );
         tmp.tlist.Append( k+3, k+4, k+5 );
      }

      *this = tmp;
   }

   pos  [0] =  0.0; pos  [1] =  0.0; pos  [2] =  0.0;
   radii[0] = +1.0; radii[1] = +1.0; radii[2] = +1.0;

   rot  [0] = +1.0; rot  [1] =  0.0; rot  [2] =  0.0;
   rot  [3] =  0.0; rot  [4] = +1.0; rot  [5] =  0.0;
   rot  [6] =  0.0; rot  [7] =  0.0; rot  [8] = +1.0;

   bbox [0] = -1.0; bbox [1] = +1.0;
   bbox [2] = -1.0; bbox [3] = +1.0;
   bbox [4] = -1.0; bbox [5] = +1.0;
}

//-----------------------------------------------------------------------------------------------

void Ellipsoid::Scale (const real8 sx, const real8 sy, const real8 sz)
{
   vlist.Scale   (sx,sy,sz); 
   vlist.Get_BBox(bbox);
}

//-----------------------------------------------------------------------------------------------

void Ellipsoid::Transform (const real8 *m44)
{
   vlist.Transform(m44);
   vlist.Get_BBox (bbox);
}

//-----------------------------------------------------------------------------------------------

int Ellipsoid::Overlap (Ellipsoid & e)
{
   if ( e.bbox[1] < bbox[0] ) return(0);
   if ( e.bbox[0] > bbox[1] ) return(0);
   if ( e.bbox[3] < bbox[2] ) return(0);
   if ( e.bbox[2] > bbox[3] ) return(0);
   if ( e.bbox[5] < bbox[4] ) return(0);
   if ( e.bbox[4] > bbox[5] ) return(0);

   return(1);
}

//-----------------------------------------------------------------------------------------------

void Ellipsoid::Print (FILE *fd)
{
   real8 *vtx = vlist.Vertices    ();
   int   *tri = tlist.Triangles   ();
   int    tn  = tlist.Triangle_Cnt();

   if (fd && vtx && tri && (tn>0))
   {
      for (int i=0,j=0; (i<tn); i++,j+=3)
      {
         real8 *v0 = vtx + 3*tri[j  ];
         real8 *v1 = vtx + 3*tri[j+1];
         real8 *v2 = vtx + 3*tri[j+2];
 
         fprintf(fd,"%6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf\n"  , v0[0],v0[1],v0[2], v1[0],v1[1],v1[2] );
         fprintf(fd,"%6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf\n"  , v1[0],v1[1],v1[2], v2[0],v2[1],v2[2] );
         fprintf(fd,"%6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf\n"  , v2[0],v2[1],v2[2], v0[0],v0[1],v0[2] );
      }
   }
}

//-----------------------------------------------------------------------------------------------

static void Save_GNU_CMD (Ellipsoid & e, const char *path)
{
   char  fpath[256]; 
   FILE *fd = 0;

   if (path)
   {
      snprintf(fpath,sizeof(fpath),"%s.gnu", path);
      fd = fopen(fpath,"w");
   };

   if (fd)
   {
      real8 min = MIN(e.bbox[0], e.bbox[2]);
            min = MIN(min      , e.bbox[4]);
      real8 max = MAX(e.bbox[1], e.bbox[3]);
            max = MAX(max      , e.bbox[5]);

      fprintf(fd,"x0=%0.2lf\n", min );
      fprintf(fd,"y0=%0.2lf\n", min );
      fprintf(fd,"z0=%0.2lf\n", min );
      fprintf(fd,"x1=%0.2lf\n", max );
      fprintf(fd,"y1=%0.2lf\n", max );
      fprintf(fd,"z1=%0.2lf\n", max );
      fprintf(fd,"\n");
      fprintf(fd,"unset key\n");
      fprintf(fd,"\n");
      fprintf(fd,"set ticslevel 0\n");
      fprintf(fd,"unset xtics\n");
      fprintf(fd,"unset ytics\n");
      fprintf(fd,"unset ztics\n");
      fprintf(fd,"\n");
      fprintf(fd,"set xrange[x0:x1]\n");
      fprintf(fd,"set yrange[y0:y1]\n");
      fprintf(fd,"set zrange[z0:z1]\n");
      fprintf(fd,"\n");
      fprintf(fd,"set arrow from x0,y0,z0  to x1,y0,z0 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"set arrow from x0,y0,z0  to x0,y1,z0 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"set arrow from x0,y1,z0  to x1,y1,z0 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"set arrow from x1,y0,z0  to x1,y1,z0 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"set arrow from x0,y0,z1  to x1,y0,z1 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"set arrow from x0,y0,z1  to x0,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"set arrow from x0,y1,z1  to x1,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"set arrow from x1,y0,z1  to x1,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"set arrow from x0,y0,z0  to x0,y0,z1 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"set arrow from x1,y0,z0  to x1,y0,z1 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"set arrow from x0,y1,z0  to x0,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"set arrow from x1,y1,z0  to x1,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb 'black'\n");
      fprintf(fd,"\n");
      fprintf(fd,"splot '%s.dat' u 1:2:3:($4-$1):($5-$2):($6-$3) w vectors nohead notitle;\n", path);
      fprintf(fd,"pause -1;\n");
      fprintf(fd,"\n");
   }

   if (fd) fclose(fd);
}

static void Save_GNU_DAT (Ellipsoid & e, const char *path)
{
   char  fpath[256]; if (path) { snprintf(fpath,sizeof(fpath),"%s.dat", path); }

   FILE *fd = ( path ? fopen(fpath,"w") : 0 );

   if (fd)
   {
      real8 *vtx = e.Vertices    ();
      int   *tri = e.Triangles   ();
      int    tn  = e.Triangle_Cnt();

      if (vtx && tri && (tn>0))
      { 
         for (int i=0,j=0; (i<tn); i++,j+=3)
         {
            real8 *v0 = vtx + 3*tri[j  ];
            real8 *v1 = vtx + 3*tri[j+1];
            real8 *v2 = vtx + 3*tri[j+2];
    
            fprintf(fd,"%6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf\n"  , v0[0],v0[1],v0[2], v1[0],v1[1],v1[2] );
            fprintf(fd,"%6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf\n"  , v1[0],v1[1],v1[2], v2[0],v2[1],v2[2] );
            fprintf(fd,"%6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf\n"  , v2[0],v2[1],v2[2], v0[0],v0[1],v0[2] );
         }
      } 
   }

   if (fd) fclose(fd);
}

void Ellipsoid::Save_GNU (const char *path)
{
   Save_GNU_CMD(*this,path);
   Save_GNU_DAT(*this,path);
}

#ifndef NO_XWINDOW

// WinDrawLinePBC()
//
// Applies the periodic boundary to the enpoints of the line to be drawn and avoids
// drawing elements that cross the periodic boundary. This could be done better but 
// it works for now. Currently assumes the coordinates are in the normalized 
// reference frame [-1:+1].
//----------------------------------------------------------------------------------------------

static void WinDrawLinePBC(real8 *v0, real8 *v1, unsigned long color, double r, unsigned long attr)
{
    real8 p0[3] = { v0[0], v0[1], v0[2] };
    real8 p1[3] = { v1[0], v1[1], v1[2] };

    p0[0] += ( (p0[0] < -1.0) ? +2.0 : ( (p0[0]>1.0) ? -2.0 : 0.0 ) );
    p0[1] += ( (p0[1] < -1.0) ? +2.0 : ( (p0[1]>1.0) ? -2.0 : 0.0 ) );
    p0[2] += ( (p0[2] < -1.0) ? +2.0 : ( (p0[2]>1.0) ? -2.0 : 0.0 ) );
    p1[0] += ( (p1[0] < -1.0) ? +2.0 : ( (p1[0]>1.0) ? -2.0 : 0.0 ) );
    p1[1] += ( (p1[1] < -1.0) ? +2.0 : ( (p1[1]>1.0) ? -2.0 : 0.0 ) );
    p1[2] += ( (p1[2] < -1.0) ? +2.0 : ( (p1[2]>1.0) ? -2.0 : 0.0 ) );

    real8 dx = fabs(p1[0]-p0[0]);
    real8 dy = fabs(p1[1]-p0[1]);
    real8 dz = fabs(p1[2]-p0[2]);

    if ( (dx<1.0) && (dy<1.0) && (dz<1.0) )
       WinDrawLine(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], color, r, attr);
}

// Ellipsoid_Display()
//
// Two forms 
//   - traverse the triangle/vertex list and issue window line displays
//   - traverse the triangle/vertex list apply the rotation/translation transformation and
//        issue window line displays
//-----------------------------------------------------------------------------------------------

void Ellipsoid::Display (unsigned long color, double r, unsigned long attr)
{
   real8 *vtx = Vertices    ();
   int   *tri = Triangles   ();
   int    tn  = Triangle_Cnt();

   if (vtx && tri && (tn>0))
   {
      for (int i=0,j=0; (i<tn); i++,j+=3)
      {
         real8 *v0 = vtx + 3*tri[j  ];
         real8 *v1 = vtx + 3*tri[j+1];
         real8 *v2 = vtx + 3*tri[j+2];

         WinDrawLinePBC (v0,v1, color,r,attr);
         WinDrawLinePBC (v1,v2, color,r,attr);
         WinDrawLinePBC (v2,v0, color,r,attr);
      }
   }
}

void Ellipsoid::Display (real8 *m33, real8 *pos, unsigned long color, double r, unsigned long attr)
{
   real8 *vtx = Vertices    ();
   int   *tri = Triangles   ();
   int    tn  = Triangle_Cnt();

   if (vtx && tri && (tn>0))
   {
      for (int i=0,j=0; (i<tn); i++,j+=3)
      {
         real8 *v0 = vtx + 3*tri[j  ];                 // v0 = triangle vertex 0
         real8 *v1 = vtx + 3*tri[j+1];                 // v1 = triangle vertex 1
         real8 *v2 = vtx + 3*tri[j+2];                 // v2 = triangle vertex 2

         real8 p0[3] = { v0[0], v0[1], v0[2] };        // p0 = local copy of vertex 0
         real8 p1[3] = { v1[0], v1[1], v1[2] };        // p1 = local copy of vertex 1
         real8 p2[3] = { v2[0], v2[1], v2[2] };        // p2 = local copy of vertex 2

         V3_M33_V3_MUL(p0,m33,p0); V3_ADD(p0,p0,pos);  // rotate+translate vertex 0
         V3_M33_V3_MUL(p1,m33,p1); V3_ADD(p1,p1,pos);  // rotate+translate vertex 1
         V3_M33_V3_MUL(p2,m33,p2); V3_ADD(p2,p2,pos);  // rotate+translate vertex 2
 
         WinDrawLinePBC(p0,p1, color,r,attr);
         WinDrawLinePBC(p1,p2, color,r,attr);
         WinDrawLinePBC(p2,p0, color,r,attr);
      }
   }
}

#endif // !NO_XWINDOW

