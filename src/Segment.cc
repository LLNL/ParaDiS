#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi_portability.h"
#include "cuda_portability.h"
#include "omp_portability.h"

#include "Typedefs.h"
#include "CPU_Timer.h"
#include "Node.h"
#include "Segment.h"
#include "V3.h"
#include "PBC.h"

//-------------------------------------------------------------------------------------------------------------------

static void Malloc_Check(void *p, const char *file, const char *func, const int ln)
{
   if (!p)
   {
      printf("%s::%s(ln=%d) memory allocation error\n", file, func, ln );
      MPI_Abort(MPI_COMM_WORLD,-1);
      exit(0);
   }
}

#define MALLOC_CHECK(p) Malloc_Check((p), __FILE__, __func__, __LINE__);

//-------------------------------------------------------------------------------------------------------------------

static void _init(Segment_t &s)
{
   s.node1                  = 0;
   s.node2                  = 0;
   s.forcesSet              = 0;
   s.sendParticleIntersects = 0;
   s.cell_indx              = 0;
   s.native                 = 0;

   V3_ZERO(s.p1 );
   V3_ZERO(s.p2 );
   V3_ZERO(s.f1 );
   V3_ZERO(s.f2 );
   V3_ZERO(s.bv );
   V3_ZERO(s.nv );

   OMP_INIT_LOCK(&s.segLock);
}

//-------------------------------------------------------------------------------------------------------------------

Segment_t::Segment_t(void) { _init(*this); }

//-------------------------------------------------------------------------------------------------------------------

Segment_t::Segment_t(const Segment_t & s)
{
   _init(*this);

   node1                  = s.node1;
   node2                  = s.node2;
   forcesSet              = s.forcesSet;
   sendParticleIntersects = s.sendParticleIntersects;
   cell_indx              = s.cell_indx;
   native                 = s.native   ;

   V3_COPY(p1 ,s.p1 );
   V3_COPY(p2 ,s.p2 );
   V3_COPY(f1 ,s.f1 );
   V3_COPY(f2 ,s.f2 );
   V3_COPY(bv ,s.bv );
   V3_COPY(nv ,s.nv );
}

//-------------------------------------------------------------------------------------------------------------------

Segment_t::Segment_t(const Node_t *n1, const Node_t *n2)
{
   _init(*this);

   node1 = (Node_t *) n1;
   node2 = (Node_t *) n2;

   if (n1) { V3_SET(p1,n1->x, n1->y, n1->z); }
   if (n2) { V3_SET(p2,n2->x, n2->y, n2->z); }

   if (n1 && n2 && n1->nbrTag)
   {
      // search for burgers vector from node 1 to node 2...

      for (int i=0; (i<n1->numNbrs); i++)
      {
         if ( n1->nbrTag[i] == n2->myTag )
         {
            V3_SET(bv,n1->burgX[i], n1->burgY[i], n1->burgZ[i]);
            V3_SET(nv,n1->nx   [i], n1->ny   [i], n1->nz   [i]);
            break;
         }
      }
   }

   if (    (n1 && n1->native)
        || (n2 && n2->native) ) native=1;
}

//-------------------------------------------------------------------------------------------------------------------

Segment_t::Segment_t(const Node_t *n1, const Node_t *n2, const real8 *b12, const real8 *n12)
{
   _init(*this);

   node1 = (Node_t *) n1;
   node2 = (Node_t *) n2;

   if (b12) { V3_COPY(bv,b12); }
   if (n12) { V3_COPY(nv,n12); }

   if (n1 ) { V3_SET(p1, n1->x, n1->y, n1->z); }
   if (n2 ) { V3_SET(p2, n2->x, n2->y, n2->z); }

   if (    (n1 && n1->native)
        || (n2 && n2->native) ) native=1;
}

//-------------------------------------------------------------------------------------------------------------------

Segment_t::Segment_t
(
   const Node_t *n1, const Node_t *n2,
   const real8 bx, const real8 by, const real8 bz,
   const real8 nx, const real8 ny, const real8 nz
)
{
   _init(*this);

   node1 = (Node_t *) n1;
   node2 = (Node_t *) n2;

   V3_SET(bv,bx,by,bz);
   V3_SET(nv,nx,ny,nz);

   if (n1) { V3_SET(p1,n1->x, n1->y, n1->z); }
   if (n2) { V3_SET(p2,n2->x, n2->y, n2->z); }

   if (    (n1 && n1->native)
        || (n2 && n2->native) ) native=1;
}

//-------------------------------------------------------------------------------------------------------------------

const Segment_t & Segment_t::operator = (const Segment_t & s)
{
   if (this!=&s)
   {
      node1                  = s.node1;
      node2                  = s.node2;
      forcesSet              = s.forcesSet;
      sendParticleIntersects = s.sendParticleIntersects;
      cell_indx              = s.cell_indx;
      native                 = s.native   ;

      V3_COPY(p1 ,s.p1 );
      V3_COPY(p2 ,s.p2 );
      V3_COPY(f1 ,s.f1 );
      V3_COPY(f2 ,s.f2 );
      V3_COPY(bv ,s.bv );
      V3_COPY(nv ,s.nv );
   }

   return(*this);
}

//-------------------------------------------------------------------------------------------------------------------

Segment_t::~Segment_t()
{
   OMP_DELETE_LOCK(&segLock);
}

// Segment_Length()
//
// Returns the segment length.  Note - presumes p1,p2 already adjusted for PBC.
//-------------------------------------------------------------------------------------------------------------------

real8 Segment_t::Segment_Length(void) const
{
   real8 dx = (p1[0]-p2[0]);
   real8 dy = (p1[1]-p2[1]);
   real8 dz = (p1[2]-p2[2]);

   return ( sqrt(dx*dx + dy*dy + dz*dz) );
}

//-------------------------------------------------------------------------------------------------------------------

void Segment_t::Reset_Forces (void)
{
   OMP_LOCK(&segLock);

   V3_ZERO(f1 );
   V3_ZERO(f2 );

   OMP_UNLOCK(&segLock);
}

// Set_Forces()
//
// Initializes the segment force vectors to the values provided.
//-------------------------------------------------------------------------------------------------------------------

void Segment_t::Set_Forces (const real8 *pf1, const real8 *pf2)
{
   OMP_LOCK(&segLock);

   if (forcesSet && pf1) { f1[0]=pf1[0]; f1[1]=pf1[1]; f1[2]=pf1[2]; }
   if (forcesSet && pf2) { f2[0]=pf2[0]; f2[1]=pf2[1]; f2[2]=pf2[2]; }

   OMP_UNLOCK(&segLock);
}

// Add_Forces()
//
// Adds/accumulates the forces into the segment force vectors.
//-------------------------------------------------------------------------------------------------------------------

void Segment_t::Add_Forces (const real8 *pf1, const real8 *pf2)
{
   OMP_LOCK(&segLock);

   if (forcesSet && pf1) { f1[0]+=pf1[0]; f1[1]+=pf1[1]; f1[2]+=pf1[2]; }
   if (forcesSet && pf2) { f2[0]+=pf2[0]; f2[1]+=pf2[1]; f2[2]+=pf2[2]; }

   OMP_UNLOCK(&segLock);
}

// Parallel()
//
// Will determine if two segments are parallel based on the critical angle provided.
// Note - this is a duplication from the segment forces routines. This routine is used
// to build up the arrays of parallel/non-parallel segment pairs for the seg-seg force
// calculations.
//---------------------------------------------------------------------------------------------------------

int Segment_t::Parallel (const Segment_t & s, const real8 ecrit) const
{
   real8  x1[3]     = {   p1[0],   p1[1],   p1[2] };                          // x1 = position of node 1 <xyz>
   real8  x2[3]     = {   p2[0],   p2[1],   p2[2] };                          // x2 = position of node 2 <xyz>
   real8  x3[3]     = { s.p1[0], s.p1[1], s.p1[2] };                          // x3 = position of node 3 <xyz>
   real8  x4[3]     = { s.p2[0], s.p2[1], s.p2[2] };                          // x4 = position of node 4 <xyz>

   real8  v1[3]     = { x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2] };              // v1        = (x2-x1)
   real8  v3[3]     = { x4[0]-x3[0], x4[1]-x3[1], x4[2]-x3[2] };              // v3        = (x4-x3)
   real8  v1dot     = (v1[0]*v1[0]) + (v1[1]*v1[1]) + (v1[2]*v1[2]);          // v1dot     = dot(v1,v1)
   real8  v3dot     = (v3[0]*v3[0]) + (v3[1]*v3[1]) + (v3[2]*v3[2]);          // v3dot     = dot(v3,v3)
   real8  oneoverLp = ( (v1dot<1.0e-20) ? 0.0 : 1.0/sqrt(v1dot) );            // oneoverLp = 1.0/||x2-x1||
   real8  oneoverL  = ( (v3dot<1.0e-20) ? 0.0 : 1.0/sqrt(v3dot) );            // oneoverL  = 1.0/||x4-x3||
   real8  tp[3]     = { v1[0]*oneoverLp, v1[1]*oneoverLp, v1[2]*oneoverLp };  // tp        = (x2-x1)/||x2-x1|| (normalized direction vector of nodes x1->x2)
   real8  t [3]     = { v3[0]*oneoverL , v3[1]*oneoverL , v3[2]*oneoverL  };  // t         = (x4-x3)/||x4-x3|| (normalized direction vector of nodes x3->x4)
   real8  c         = t[0]*tp[0] + t[1]*tp[1] + t[2]*tp[2];                   // c         = dot(t,tp)
   real8  c2        = c*c;                                                    // c2        = c*c
          c2        = ( (c2<1.0) ? c2 : 1.0 );                                // c2        = min(c*c,1.0)
   real8  onemc2    = 1.0-c2;                                                 // onemc2    = 1.0-c*c

   return ( (onemc2<ecrit) ? 1 : 0 );
}

// PBC_Position()
//
// Will adjust the location of the second node for periodic boundary conditions. Note that this
// only affects the local copy of the node position stored in the segment (p2).  This does not change
// the actual position of the node in the node structure.
//
// Note there are two overloaded forms - one that includes the PBC reciprocal of simulation sizes (sx,sy,sz)
// and one that doesn't.  The version that doesn't include the PBC reciprocal's assumes PBC is always active
// (which is needed in some cases where printing the wrapped PBC segment is required (e.g. gnuplot output)).
//-------------------------------------------------------------------------------------------------------------------

void Segment_t::PBC_Position
(
   const real8 lx, const real8 ly, const real8 lz   ///< overall simulation size (xyz)
)
{
   real8 sx = ( (fabs(lx)>0.0) ? 1.0/fabs(lx) : 0.0 );   // sx = reciprocal of simulation size (x)
   real8 sy = ( (fabs(ly)>0.0) ? 1.0/fabs(ly) : 0.0 );   // sy = reciprocal of simulation size (y)
   real8 sz = ( (fabs(lz)>0.0) ? 1.0/fabs(lz) : 0.0 );   // sz = reciprocal of simulation size (z)

   ::PBC_Position(p1,p2, lx,ly,lz, sx,sy,sz);
}

void Segment_t::PBC_Position
(
   const real8 lx, const real8 ly, const real8 lz,  ///< overall simulation size (xyz)
   const real8 sx, const real8 sy, const real8 sz   ///< pbc reciprocal of simulation size (xyz)
)
{
   ::PBC_Position(p1,p2, lx,ly,lz, sx,sy,sz);
}

// Set_Cell()
//
// Computes the cell index based on the overall simulation size and cell geometry. Note that the assigned cell
// index is based on the segment midpoint.
//-------------------------------------------------------------------------------------------------------------------

int Segment_t::Set_Cell
(
   const int   nx,  ///< cell count <x>
   const int   ny,  ///< cell count <y>
   const int   nz,  ///< cell count <z>
   const real8 lx,  ///< simulation size <x>
   const real8 ly,  ///< simulation size <y>
   const real8 lz   ///< simulation size <z>
)
{
   real8 x1[3] = { p1[0], p1[1], p1[2] };                // x1 = local copy of node postion 1
   real8 x2[3] = { p2[0], p2[1], p2[2] };                // x2 = local copy of node postion 2

   real8 sx = ( (fabs(lx)>0.0) ? 1.0/fabs(lx) : 0.0 );   // sx = reciprocal of simulation size (x)
   real8 sy = ( (fabs(ly)>0.0) ? 1.0/fabs(ly) : 0.0 );   // sy = reciprocal of simulation size (y)
   real8 sz = ( (fabs(lz)>0.0) ? 1.0/fabs(lz) : 0.0 );   // sz = reciprocal of simulation size (z)

   ::PBC_Position(x1,x2, lx,ly,lz,  sx,sy,sz);           // update x2 for periodic boundary

   real8 cx = (x2[0]+x1[0])/2.0;                         // cx = midpoint of segment <x>
   real8 cy = (x2[1]+x1[1])/2.0;                         // cy = midpoint of segment <y>
   real8 cz = (x2[2]+x1[2])/2.0;                         // cz = midpoint of segment <z>

   int   ix = (int)((cx+(lx/2.0))/(lx/nx))%nx;           // ix = cell index <x>
   int   iy = (int)((cy+(ly/2.0))/(ly/ny))%ny;           // iy = cell index <y>
   int   iz = (int)((cz+(lz/2.0))/(lz/nz))%nz;           // iz = cell index <z>

   cell_indx = (iz*nx*ny)+(iy*nx)+ix;                    // cell_indx  = overall cell index

   return(cell_indx);
}

//-------------------------------------------------------------------------------------------------------------------

static inline int Node_Inside
(
   const real8 x , const real8 y , const real8 z ,    ///< node position <xyz>
   const real8 x0, const real8 y0, const real8 z0,    ///< lower bounds of domain
   const real8 x1, const real8 y1, const real8 z1     ///< upper bounds of domain
)
{
   if ( (x<x0) || (x1<=x) ) return(0);
   if ( (y<y0) || (y1<=y) ) return(0);
   if ( (z<z0) || (z1<=z) ) return(0);

   return(1);
}

//-------------------------------------------------------------------------------------------------------------------

void Segment_t::Set_Native (const int dom)
{
   native = (    (node1 && (node1->myTag.domainID==dom) )
              || (node2 && (node2->myTag.domainID==dom) ) ? 1 : 0 );
}

void Segment_t::Set_Native
(
   const real8 x0, const real8 x1,
   const real8 y0, const real8 y1,
   const real8 z0, const real8 z1
)
{
   native =    ( node1 ? Node_Inside(node1->x, node1->y, node1->z, x0,y0,z0, x1,y1,z1) : 0 )
            || ( node2 ? Node_Inside(node2->x, node2->y, node2->z, x0,y0,z0, x1,y1,z1) : 0 );
}

// Print()
//
// Will print the current segment to the output file descriptor.
//-------------------------------------------------------------------------------------------------------------------

void Segment_t::Print(FILE *fd) const
{
   if (fd && node1 && node2)
   {
      char t1[32]; sprintf(t1,"(%d,%d)", node1->myTag.domainID, node1->myTag.index );
      char t2[32]; sprintf(t2,"(%d,%d)", node2->myTag.domainID, node2->myTag.index );

      fprintf(fd,"%12s %12s %d %d"                                   , t1, t2, forcesSet, sendParticleIntersects );
      fprintf(fd," %14.8lf %14.8lf %14.8lf %14.8lf %14.8lf %14.8lf"  , p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]  );
      fprintf(fd," %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n", bv[0], bv[1], bv[2], nv[0], nv[1], nv[2]  );
   }
}

// Print()
//
// Will print the current segment to the output file descriptor.
// Will adjust the output for periodic boundaries on the simulation size.
//-------------------------------------------------------------------------------------------------------------------

void Segment_t::Print
(
         FILE  *fd,     ///< points to a currently open file descriptor
   const real8  lx,     ///< size of the simulation box (x) (needed for PBC corrections)
   const real8  ly,     ///< size of the simulation box (y) (needed for PBC corrections)
   const real8  lz      ///< size of the simulation box (z) (needed for PBC corrections)
) const
{
   if (fd && node1 && node2)
   {
      char t1[32]; sprintf(t1,"(%d,%d)", node1->myTag.domainID, node1->myTag.index);
      char t2[32]; sprintf(t2,"(%d,%d)", node2->myTag.domainID, node2->myTag.index);

      real8 x1[3] = { p1[0], p1[1], p1[2] };                // x1  = local copy of p1
      real8 x2[3] = { p2[0], p2[1], p2[2] };                // x2  = local copy of p2

      real8 sx = ( (fabs(lx)>0.0) ? 1.0/fabs(lx) : 0.0 );   // sx = reciprocal of simulation size (x)
      real8 sy = ( (fabs(ly)>0.0) ? 1.0/fabs(ly) : 0.0 );   // sy = reciprocal of simulation size (y)
      real8 sz = ( (fabs(lz)>0.0) ? 1.0/fabs(lz) : 0.0 );   // sz = reciprocal of simulation size (z)

      ::PBC_Position(x1,x2, lx,ly,lz,  sx,sy,sz);           // update x2 for periodic boundary

      fprintf(fd,"%-12s %-12s %d %d"                                 , t1, t2, forcesSet, sendParticleIntersects );
      fprintf(fd," %14.8lf %14.8lf %14.8lf %14.8lf %14.8lf %14.8lf"  , x1[0], x1[1], x1[2], x2[0], x2[1], x2[2]  );
      fprintf(fd," %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf\n", bv[0], bv[1], bv[2], nv[0], nv[1], nv[2]  );
   }
}

// Segment_Append()
//
// Will append a segment to an existing segment list (array).  Reallocates and extends the array as necessary.
// This could be more efficient by implementing the lists as linked lists, but too much of the code base
// assumes the ability to iterate via array elements.
//-------------------------------------------------------------------------------------------------------------------

Segment_t *Segment_Append (Segment_t *segs, int & cnt, int & max, Segment_t & seg)
{
   if ( !segs || (cnt==max) )
   {
      max += ( (max==0) ? 1000 : (max/2) );

      Segment_t *tmp = ( (max>0) ? new Segment_t[max] : 0 );

      MALLOC_CHECK(tmp);

      if (tmp && segs)
         for (int i=0; (i<cnt); i++) { tmp[i]=segs[i]; }

      if (segs) { delete [] segs; }

      segs = tmp;
   }

   if (segs && (cnt<max) ) { segs[cnt++]=seg; }

   return(segs);
}

// Find_Node()
//
// Performs a simple binary search of the source node array.
// Returns a pointer to the node if found, null if not found.
//---------------------------------------------------------------------------------------------------------

static Node_t *Find_Node (const Node_t **narr, const int ncnt, const Tag_t & tag)
{
   if (narr && (ncnt>0))
   {
      int i0=0;         // i0=lower index
      int i1=(ncnt-1);  // i1=upper index

      while (i0<=i1)
      {
         int i = (i0+i1)/2;   // i = midpoint index

         Node_t *node = (Node_t *) narr[i];
         int  cmp = ( (node->myTag<tag) ? -1 :
                      (node->myTag>tag) ? +1 : 0 );

         if (cmp==0) { return(node); }  // found!
         if (cmp< 0) { i0=i+1;       }  // adjust lower index
         else        { i1=i-1;       }  // adjust upper index
      }
   }

   return(0);  // not found!
}

// Segment_List_Create()
//
// Will construct and return a unique segment list from a given an array of sorted node pointers.
//
// A segment is considered unique if the tag for the first node is less than the tag of the
// second node. This insures that we don't wind up with duplicate segments e.g. (A->B) and (B->A).
//-------------------------------------------------------------------------------------------------------------------

Segment_t *Segment_List_Create
(
         int    & scnt,   ///< number of resulting segments (returned)
   const Node_t **narr,   ///< array  of node pointers
   const int      ncnt    ///< number of node pointers
)
{
   Segment_t *segs=0;  int smax=0; scnt=0;

   if (narr)
   {
      for (int i=0; (i<ncnt); i++)               // loop through all the nodes...
      {
         Node_t *n1 = (Node_t *) narr[i];        // n1 = points to first node

         if (n1)
         {
            int narms = n1->numNbrs;             // narms = number of arms

            for (int j=0; (j<narms); j++)        // loop through all the arms...
            {
               if (n1->myTag < n1->nbrTag[j])
               {
                  Node_t *n2 = Find_Node(narr,ncnt,n1->nbrTag[j]);

                  if (n1 && n2)
                  {
                     Segment_t seg = Segment_t(n1,n2,
                                               n1->burgX[j], n1->burgY[j], n1->burgZ[j],
                                               n1->nx   [j], n1->ny   [j], n1->nz   [j] );

                     segs = Segment_Append(segs,scnt,smax,seg);
                  }
               }
            }
         }
      }
   }

   return(segs);
}

// Segment_List_Create()
//
// Returns a list of segments given an array of native and ghost nodes.
//-------------------------------------------------------------------------------------------------------------------

Segment_t *Segment_List_Create
(
         int    & seg_cnt    ,   ///< number of resulting segments (returned)
   const Node_t **native     ,   ///< sparse array of native node pointers
   const int      native_cnt ,   ///< number of native node pointers
   const Node_t **ghosts     ,   ///< sparse array of ghost node pointers
   const int      ghost_cnt      ///< number of ghost node pointers
)
{
   Segment_t *segs=0; seg_cnt=0;

   // compose a dense array of node pointers (skips over null array entries)

   Node_t **narr = ( ((native_cnt+ghost_cnt)>0) ? new Node_t *[native_cnt+ghost_cnt] : 0 );

   if (narr)
   {
      int ncnt=0;
      if (native && (native_cnt>0)) for (int i=0; (i<native_cnt); ++i) { if (native[i]) narr[ncnt++] = (Node_t *) native[i]; }
      if (ghosts && (ghost_cnt >0)) for (int i=0; (i<ghost_cnt ); ++i) { if (ghosts[i]) narr[ncnt++] = (Node_t *) ghosts[i]; }

      narr = Nodes_Sort(narr,ncnt);  // (sort the array)

      segs = Segment_List_Create( seg_cnt, (const Node_t **) narr, ncnt);
   }

   if (narr) { delete [] narr; narr=0; }

   return(segs);
}

//-------------------------------------------------------------------------------------------------------------------

void Segment_List_Stats
(
   const Node_t **native     ,   ///< sparse array of native node pointers
   const int      native_cnt ,   ///< number of native node pointers
   const Node_t **ghosts     ,   ///< sparse array of ghost node pointers
   const int      ghost_cnt      ///< number of ghost node pointers
)
{
   Segment_t *segs=0; int scnt=0;

   real8      smin=0.0;  // shortest segment length
   real8      smax=0.0;  // longest  segment length
   real8      ssum=0.0;  // sum of  segment lengths
   int        sndx=0;    // index of longest segment

   segs = Segment_List_Create(scnt, native, native_cnt,
                                    ghosts, ghost_cnt );


   if (segs) { delete [] segs; segs=0; }
}

// Cell_Index()
//
// Converts a cell index triplet <xyz> to an absolute cell index.
//-------------------------------------------------------------------------------------------------------------------

static int Cell_Index
(
   int ix, int iy, int iz,  ///< cell index triplet <xyz>
   int nx, int ny, int nz   ///< cell counts        <xyz>
)
{
   ix = ( (ix<0) ? ((ix+nx)%nx) : (ix%nx) );  // constrain/wrap cell index <x>
   iy = ( (iy<0) ? ((iy+ny)%ny) : (iy%ny) );  // constrain/wrap cell index <y>
   iz = ( (iz<0) ? ((iz+nz)%nz) : (iz%nz) );  // constrain/wrap cell index <z>

   return((iz*nx*ny)+(iy*nx)+ix);
}

// Active_Cells()
//
// Returns an nx*ny*nz matrix of booleans that identify which cells are active within a
// given domain.  All cells that fall within the domain are set to 1.  This data structure
// can then be used to test whether a segment falls within the bounds after adjusting for
// periodic boundaries.
//-------------------------------------------------------------------------------------------------------------------

unsigned char *Active_Cells
(
   const real8  x0, const real8 y0, const real8 z0,   ///< lower bound of domain
   const real8  x1, const real8 y1, const real8 z1,   ///< upper bound of domain
   const int    nx, const int   ny, const int   nz,   ///< number of cells <xyz>
   const real8  lx, const real8 ly, const real8 lz    ///< simulation size <xyz>
)
{
   int ix0 = (int)((x0+(lx/2.0))/(lx/nx))-1;             // ix0 = lower bound of cells <x>
   int ix1 = (int)((x1+(lx/2.0))/(lx/nx))+1;             // ix1 = upper bound of cells <x>
   int iy0 = (int)((y0+(ly/2.0))/(ly/ny))-1;             // iy0 = lower bound of cells <y>
   int iy1 = (int)((y1+(ly/2.0))/(ly/ny))+1;             // iy1 = upper bound of cells <y>
   int iz0 = (int)((z0+(lz/2.0))/(lz/nz))-1;             // iz0 = lower bound of cells <z>
   int iz1 = (int)((z1+(lz/2.0))/(lz/nz))+1;             // iz1 = upper bound of cells <z>

   unsigned char *cells = ( (nx*ny*nz>0) ? new unsigned char[nx*ny*nz] : 0 );

   if (cells)
   {
      memset(cells,0,nx*ny*nz*sizeof(unsigned char));    // set all cells to inactive
                                                         //
      for (int k=iz0; (k<=iz1); k++)                     // loop over all the included cells
      for (int j=iy0; (j<=iy1); j++)                     //
      for (int i=ix0; (i<=ix1); i++)                     //
      {                                                  //
          int cndx = Cell_Index(i,j,k, nx,ny,nz);        // cndx = cell index of (i,j,k) triplet
                                                         //
          if ( (0<=cndx) && (cndx<(nx*ny*nz)) )          // if (cell index) legal
             cells[cndx] = 1;                            //    set cell to active
      }
   }

   return(cells);
}

// Segment_List_Crop()
//
// Returns a segment list containing only the segments that fall within the cells encapsulated by the
// domain (x0,y0,z0, x1,y1,z1).
//-------------------------------------------------------------------------------------------------------------------

Segment_t *Segment_List_Crop
(
         Segment_t   *segs,                           ///< source segment array
   const int          scnt,                           ///< number of segments in the source segment array
         int        & ncnt,                           ///< number of segments in cropped array (returned)
   const real8  x0, const real8 y0, const real8 z0,   ///< lower bound of domain
   const real8  x1, const real8 y1, const real8 z1,   ///< upper bound of domain
   const int    nx, const int   ny, const int   nz,   ///< number of cells <xyz>
   const real8  lx, const real8 ly, const real8 lz,   ///< simulation size <xyz>
         real8 *msec                                  ///< time in millisecs (if non-null)
)
{
   CPU_Timer tmr;
   CPU_TIMER_START(msec,tmr);

   Segment_t *nsegs=0; int nmax=0; ncnt=0;

   unsigned char *cells = Active_Cells(x0,y0,z0, x1,y1,z1, nx,ny,nz, lx,ly,lz);

   if (cells && segs && (scnt>0) )
   {
      for (int i=0; (i<scnt); i++)
      {
         int cndx = segs[i].Set_Cell(nx,ny,nz, lx,ly,lz);

         if ( (0<=cndx) && (cndx<(nx*ny*nz)) && (cells[cndx]) )
            nsegs = Segment_Append(nsegs,ncnt,nmax,segs[i]);
      }
   }

   if (cells) { delete [] cells; cells=0; }

   CPU_TIMER_STOP(msec,tmr);

   return(nsegs);
}

// Segment_List_Set_Cells()
//
// Sets the cell index for all the segments in a segment array
//-------------------------------------------------------------------------------------------------------------------

void Segment_List_Set_Cells
(
   Segment_t   *segs, const int   scnt,                     ///< segment array (and length)
   const int    nx  , const int   ny  , const int   nz,     ///< simulation cell geometry <xyz>
   const real8  lx  , const real8 ly  , const real8 lz,     ///< simulation size          <xyz>
         real8 *msec                                        ///< time in milliseconds (if non-null)
)
{
   CPU_Timer tmr;
   CPU_TIMER_START(msec,tmr);

   if (segs && (scnt>0))
   {
      for (int i=0; (i<scnt); i++)
          segs[i].Set_Cell(nx,ny,nz,lx,ly,lz);
   }

   CPU_TIMER_STOP(msec,tmr);
}

//-------------------------------------------------------------------------------------------------------------------

static inline int Node_Inside
(
   const real8 *p,                                     ///< node position <xyz>
   const real8  x0, const real8 y0, const real8 z0,    ///< lower bounds of domain
   const real8  x1, const real8 y1, const real8 z1     ///< upper bounds of domain
)
{
   if ( (p[0]<x0) || (x1<=p[0]) ) return(0);
   if ( (p[1]<y0) || (y1<=p[1]) ) return(0);
   if ( (p[2]<z0) || (z1<=p[2]) ) return(0);

   return(1);
}

// Segment_List_Set_Native()
//
// Will set all the segments within a segment list that have at least one node within the provided
// domain as native.
//-------------------------------------------------------------------------------------------------------------------

void Segment_List_Set_Native
(
   Segment_t   *segs,                                       ///< array of segments
   const int    scnt,                                       ///< number of segments
   const real8  x0, const real8 y0, const real8 z0,         ///< domain minimum
   const real8  x1, const real8 y1, const real8 z1,         ///< domain maximum
         real8 *msec                                        ///< time in milliseconds (if non-null)
)
{
   CPU_Timer tmr;
   CPU_TIMER_START(msec,tmr);

   if (segs && (scnt>0))
   {
      for (int i=0; (i<scnt); i++)
      {
          segs[i].native =    Node_Inside(segs[i].p1, x0,y0,z0, x1,y1,z1)
                           || Node_Inside(segs[i].p2, x0,y0,z0, x1,y1,z1);
      }
   }

   CPU_TIMER_STOP(msec,tmr);
}

// BV_Equiv()
//
// Returns 1 or 0 if two given burgers vectors are directionally equivalent.
//-------------------------------------------------------------------------------------------------------------------

static
int BV_Equiv (const real8 *b1, const real8 *b2)
{
   if (b1 && b2)
   {
      real8 t1[3]; V3_NORMALIZE(t1,b1);
      real8 t2[3]; V3_NORMALIZE(t2,b2);

      real8 pdot = V3_DOT(t1,t2);

      return ( ( fabs(fabs(pdot)-1.0) < 1.0e-4 ) ? 1 : 0 );
   }

   return(0);
}

// Segment_List_Cnt()
//
// Returns the number of segments that match a specific burgers vector in a segment list.
//-------------------------------------------------------------------------------------------------------------------

int Segment_List_Cnt
(
   const Segment_t  *segs,  ///< array of segments
   const int         scnt,  ///< number of segments
   const real8      *bv     ///< burgers vector (optional)
)
{
   int n=0;

   if (segs && (scnt>0) && bv )
   {
      for (int i=0; (i<scnt); i++)
      {
         if ( BV_Equiv(bv,segs[i].bv) ) { n++; }
      }
   }
   else { n=scnt; }

   return(n);
}

// Segment_List_Sum()
//
// Returns the sum of the segments that match a specific burgers vector in a segment list.
// If no burgers vector is provided, returns the sum of ALL segments.
//-------------------------------------------------------------------------------------------------------------------

real8 Segment_List_Sum
(
   const Segment_t  *segs,  ///< array of segments
   const int         scnt,  ///< number of segments
   const real8      *bv     ///< burgers vector (optional)
)
{
   real8 sum=0.0;

   if (segs && (scnt>0) )
   {
      // If the burgers vector was provided, sum all segments with the given bv...

      if (bv)
      {
         for (int i=0; (i<scnt); i++)
         {
            if ( BV_Equiv(bv,segs[i].bv) )
            { sum += segs[i].Segment_Length(); }
         }
      }

      // otherwise, sum ALL the segments...

      else
      {
         for (int i=0; (i<scnt); i++)
         { sum += segs[i].Segment_Length(); }
      }
   }

   return(sum);
}

// Segment_List_Density()
//
// Will return the density of segments that match a specific burgers vector in a segment list.
// If no burgers vector is provided, returns the density of ALL segments.
//-------------------------------------------------------------------------------------------------------------------

real8 Segment_List_Density
(
   const Segment_t  *segs,  ///< array of segments
   const int         scnt,  ///< number of segments
   const real8      *bv  ,  ///< burgers vector (optional)
   const real8       bmag,  ///< burgers vector magnitude
   const real8       bvol   ///< simulation volume (lx x ly x lz)
)
{
   real8 sum  = Segment_List_Sum(segs,scnt,bv);   // sum = sum lengths of segments
   real8 s    = fabs(bmag * bmag * bvol);         // s   = bmag^2*volume
   real8 rho  = ( (s>0.0) ? (sum/s) : 0.0 );      // density = sum/(bmag^2*volume)

   return (rho);
}

// Segment_List_Moment()
//
// Will travese a segment list compiling statistical data including min, max, and mean segment length.
//-------------------------------------------------------------------------------------------------------------------

void Segment_List_Moment
(
   const Segment_t  *segs,  ///< array of segments
   const int         scnt,  ///< segment count
         real8     & min ,  ///< minimum segment length  (returned)
         real8     & max ,  ///< maximum segment length  (returned)
         real8     & mean,  ///< mean    segment length  (returned)
         real8     & sum    ///< sum of  segment lengths (returned)
)
{
   min = max = mean = sum = 0.0;

   if (segs && (scnt>0) )
   {
      min = max = segs[0].Segment_Length();

      sum=0.0;
      for (int i=0; (i<scnt); i++)
      {
         real8 s = segs[i].Segment_Length();

         min  = ( (s<min) ? s : min );
         max  = ( (s>max) ? s : max );
         sum += s;
      }

      mean = (sum/scnt);
   }
}

void Segment_List_Moment
(
   const Segment_t   *segs,  ///< array of segments
   const int          scnt,  ///< segment count
         real8      & min ,  ///< minimum  segment length  (returned)
         real8      & max ,  ///< maximum  segment length  (returned)
         real8      & mean,  ///< mean     segment length  (returned)
         real8      & sum ,  ///< sum   of segment lengths (returned)
         real8      & sdev   ///< stdev of segment lengths (returned)
)
{
   min = max = mean = sum = sdev = 0.0;

   Segment_List_Moment(segs,scnt,min,max,mean,sum);

   if (segs && (scnt>0) )
   {
      real8 var=0.0;
      for (int i=0; (i<scnt); i++) { real8 s = (segs[i].Segment_Length()-mean); var += (s*s); }

      sdev = sqrt(var/scnt);
   }
}

// Segment_List_Save()
//
// Saves a segment list to an output file.
//-------------------------------------------------------------------------------------------------------------------

void Segment_List_Save
(
   const char      *path,   ///< path to the save file
   const Segment_t *segs,   ///< array of segments
   const int        scnt,   ///< number of segments
   const real8      lx  ,   ///< simulation extent <x>
   const real8      ly  ,   ///< simulation extent <y>
   const real8      lz      ///< simulation extent <z>
)
{
   char path_dat_n[256]; if (path) { sprintf(path_dat_n,"%s_n.dat", path); }
   char path_dat_g[256]; if (path) { sprintf(path_dat_g,"%s_g.dat", path); }
   char path_gnu  [256]; if (path) { sprintf(path_gnu  ,"%s.gnu"  , path); }

   if (segs && (scnt>0))
   {
      FILE *fd_n = ( path ? fopen(path_dat_n,"w") : 0 );
      FILE *fd_g = ( path ? fopen(path_dat_g,"w") : 0 );

      for (int i=0; (i<scnt); i++)
         segs[i].Print( (segs[i].native ? fd_n : fd_g), lx,ly,lz);

      fclose(fd_n);
      fclose(fd_g);
   }

   FILE *fd = ( path ? fopen(path_gnu,"w") : 0 );

   if (fd)
   {
      int x0 = (int) (-lx/2);  x0 = 1000*((x0/1000)-1);
      int y0 = (int) (-ly/2);  y0 = 1000*((y0/1000)-1);
      int z0 = (int) (-lz/2);  z0 = 1000*((z0/1000)-1);

      int x1 = (int) ( lx/2);  x1 = 1000*((x1/1000)+1);
      int y1 = (int) ( ly/2);  y1 = 1000*((y1/1000)+1);
      int z1 = (int) ( lz/2);  z1 = 1000*((z1/1000)+1);

      fprintf(fd,"#set term pngcairo size 800,640\n");
      fprintf(fd,"#set out '%s.png'\n", path );
      fprintf(fd,"\n");
      fprintf(fd,"reset\n");
      fprintf(fd,"\n");
      fprintf(fd,"set view 60,75,1.2 \n");
      fprintf(fd,"\n");
      fprintf(fd,"unset key \n");
      fprintf(fd,"\n");
      fprintf(fd,"set title  'Segments' \n");
      fprintf(fd,"set xlabel 'X' \n");
      fprintf(fd,"set ylabel 'Y' \n");
      fprintf(fd,"set zlabel 'Z' \n");
      fprintf(fd,"\n");
      fprintf(fd,"set ticslevel 0    \n");
      fprintf(fd,"set xrange [%d:%d] \n", x0, x1);
      fprintf(fd,"set yrange [%d:%d] \n", y0, y1);
      fprintf(fd,"set zrange [%d:%d] \n", z0, z1);
      fprintf(fd,"\n");
      fprintf(fd,"lx  = %d\n", (int) lx );
      fprintf(fd,"ly  = %d\n", (int) ly );
      fprintf(fd,"lz  = %d\n", (int) lz );
      fprintf(fd,"\n");
      fprintf(fd,"x0  = -lx/2 \n");
      fprintf(fd,"y0  = -ly/2 \n");
      fprintf(fd,"z0  = -lz/2 \n");
      fprintf(fd,"x1  =  lx/2 \n");
      fprintf(fd,"y1  =  ly/2 \n");
      fprintf(fd,"z1  =  lz/2 \n");
      fprintf(fd,"\n");
      fprintf(fd,"s   = lx/60 \n");
      fprintf(fd,"\n");
      fprintf(fd,"set arrow from x0,y0,z0  to x1,y0,z0 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y0,z0  to x0,y1,z0 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y1,z0  to x1,y1,z0 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x1,y0,z0  to x1,y1,z0 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y0,z1  to x1,y0,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y0,z1  to x0,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y1,z1  to x1,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x1,y0,z1  to x1,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y0,z0  to x0,y0,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x1,y0,z0  to x1,y0,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y1,z0  to x0,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x1,y1,z0  to x1,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"\n");
      fprintf(fd,"path_n = '%s_n.dat'\n", path);
      fprintf(fd,"path_g = '%s_g.dat'\n", path);
      fprintf(fd,"\n");
      fprintf(fd,"splot path_n u 5:6:7:($8-$5):($9-$6):($10-$7) w vectors nohead lt 1 lw 1   lc rgb '#0000ff' , \\\n");
      fprintf(fd,"      path_n u 5:6:7                          w points         pt 6 ps 0.5 lc rgb '#0000ff' , \\\n");
      fprintf(fd,"      path_n u 8:9:10                         w points         pt 6 ps 0.5 lc rgb '#0000ff' , \\\n");
      fprintf(fd,"      path_g u 5:6:7:($8-$5):($9-$6):($10-$7) w vectors nohead lt 1 lw 1   lc rgb '#ff0000' , \\\n");
      fprintf(fd,"      path_g u 5:6:7                          w points         pt 6 ps 0.5 lc rgb '#ff0000' , \\\n");
      fprintf(fd,"      path_g u 8:9:10                         w points         pt 6 ps 0.5 lc rgb '#ff0000' ;\n"   );
      fprintf(fd,"\n");

      fclose(fd);
   }
}

void Segment_Pairs_Save
(
   const          char       *path ,    ///< path to the save file
   const unsigned short      *pairs,    ///< array of segment pairs
   const          int         pcnt ,    ///< number of segment pairs
   const          Segment_t  *segs ,    ///< array of segments
   const          int         scnt ,    ///< number of segments
   const          real8       lx   ,    ///< simulation size (x)
   const          real8       ly   ,    ///< simulation size (y)
   const          real8       lz        ///< simulation size (z)
)
{
   char path_dat_n[256]; if (path) { sprintf(path_dat_n,"%s_n.dat", path); }
   char path_dat_g[256]; if (path) { sprintf(path_dat_g,"%s_g.dat", path); }
   char path_gnu  [256]; if (path) { sprintf(path_gnu  ,"%s.gnu"  , path); }

   FILE *fd_n = ( path ? fopen(path_dat_n,"w") : 0 );
   FILE *fd_g = ( path ? fopen(path_dat_g,"w") : 0 );

   if (pairs && (pcnt>0) && segs && (scnt>0))
   {
      for (int i=0,k=0; (i<pcnt); i++)
      {
         int i0=pairs[k++];
         int i1=pairs[k++];

         segs[i0].Print( (segs[i0].native ? fd_n : fd_g), lx,ly,lz);
         segs[i1].Print( (segs[i1].native ? fd_n : fd_g), lx,ly,lz);
      }
   }

   fclose(fd_n);
   fclose(fd_g);

   FILE *fd = ( path ? fopen(path_gnu,"w") : 0 );

   if (fd)
   {
      int x0 = (int) (-lx/2);  x0 = 1000*((x0/1000)-1);
      int y0 = (int) (-ly/2);  y0 = 1000*((y0/1000)-1);
      int z0 = (int) (-lz/2);  z0 = 1000*((z0/1000)-1);

      int x1 = (int) ( lx/2);  x1 = 1000*((x1/1000)+1);
      int y1 = (int) ( ly/2);  y1 = 1000*((y1/1000)+1);
      int z1 = (int) ( lz/2);  z1 = 1000*((z1/1000)+1);

      fprintf(fd,"#set term pngcairo size 800,640\n");
      fprintf(fd,"#set out '%s.png'\n", path );
      fprintf(fd,"\n");
      fprintf(fd,"reset\n");
      fprintf(fd,"\n");
      fprintf(fd,"set view 60,75,1.2 \n");
      fprintf(fd,"\n");
      fprintf(fd,"unset key \n");
      fprintf(fd,"\n");
      fprintf(fd,"set title  'Segments' \n");
      fprintf(fd,"set xlabel 'X' \n");
      fprintf(fd,"set ylabel 'Y' \n");
      fprintf(fd,"set zlabel 'Z' \n");
      fprintf(fd,"\n");
      fprintf(fd,"set ticslevel 0    \n");
      fprintf(fd,"set xrange [%d:%d] \n", x0, x1);
      fprintf(fd,"set yrange [%d:%d] \n", y0, y1);
      fprintf(fd,"set zrange [%d:%d] \n", z0, z1);
      fprintf(fd,"\n");
      fprintf(fd,"lx  = %d\n", (int) lx );
      fprintf(fd,"ly  = %d\n", (int) ly );
      fprintf(fd,"lz  = %d\n", (int) lz );
      fprintf(fd,"\n");
      fprintf(fd,"x0  = -lx/2 \n");
      fprintf(fd,"y0  = -ly/2 \n");
      fprintf(fd,"z0  = -lz/2 \n");
      fprintf(fd,"x1  =  lx/2 \n");
      fprintf(fd,"y1  =  ly/2 \n");
      fprintf(fd,"z1  =  lz/2 \n");
      fprintf(fd,"\n");
      fprintf(fd,"s   = lx/60 \n");
      fprintf(fd,"\n");
      fprintf(fd,"set arrow from x0,y0,z0  to x1,y0,z0 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y0,z0  to x0,y1,z0 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y1,z0  to x1,y1,z0 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x1,y0,z0  to x1,y1,z0 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y0,z1  to x1,y0,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y0,z1  to x0,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y1,z1  to x1,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x1,y0,z1  to x1,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y0,z0  to x0,y0,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x1,y0,z0  to x1,y0,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x0,y1,z0  to x0,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"set arrow from x1,y1,z0  to x1,y1,z1 nohead back nofilled lt 1 lw 1 lc rgb '#ff0080' \n");
      fprintf(fd,"\n");
      fprintf(fd,"path_n = '%s_n.dat'\n", path);
      fprintf(fd,"path_g = '%s_g.dat'\n", path);
      fprintf(fd,"\n");
      fprintf(fd,"splot path_n u 5:6:7:($8-$5):($9-$6):($10-$7) w vectors nohead lt 1 lw 1   lc rgb '#0000ff' , \\\n");
      fprintf(fd,"      path_g u 5:6:7:($8-$5):($9-$6):($10-$7) w vectors nohead lt 1 lw 1   lc rgb '#ff0000' ;\n"   );
      fprintf(fd,"\n");
      fprintf(fd,"#     path_n u 5:6:7                          w points         pt 6 ps 0.5 lc rgb '#0000ff' , \\\n");
      fprintf(fd,"#     path_n u 8:9:10                         w points         pt 6 ps 0.5 lc rgb '#0000ff' , \\\n");
      fprintf(fd,"#     path_g u 5:6:7                          w points         pt 6 ps 0.5 lc rgb '#ff0000' , \\\n");
      fprintf(fd,"#     path_g u 8:9:10                         w points         pt 6 ps 0.5 lc rgb '#ff0000' , \\\n");
      fprintf(fd,"\n");

      fclose(fd);
   }
}

// Append_Segment_Pair()
//
// Will append a segment pair to an existing segment pair list. Allocates and extends as necessary.
// Note that each entry in the segment pair list contains two indices (i,j). Adding a pair adds two
// indices to the pair list.
//-------------------------------------------------------------------------------------------------------------------

static unsigned short *Append_Segment_Pair
(
   unsigned short  *spairs,     ///< points to the segment pair array
            int   & npairs,     ///< current number of pairs in the segment pair array
            int   & nmax  ,     ///< maximum number of pairs in the segment pair array
   const    int     i     ,     ///< index of segment 1
   const    int     j           ///< index of segment 2
)
{
   // if needed, extend the size of the existing index array...

   if (!spairs || (npairs>=nmax))
   {
      nmax += ( (nmax==0) ? (1024*1024) : (nmax/2) );

      unsigned short *tmp = ( (nmax>0) ? new unsigned short[2*nmax] : 0 );

      MALLOC_CHECK(tmp);

      if (tmp && spairs) { memcpy(tmp,spairs,2*npairs*sizeof(unsigned short)); }

      if (spairs) { delete [] spairs; }

      spairs = tmp;
   }

   // append the index...

   if (spairs && (npairs<nmax))
   {
      spairs[(2*npairs)  ]=i;
      spairs[(2*npairs)+1]=j; npairs++;
   }

   return(spairs);
}

// Cell_Adjacent()
//
// Given two cell indices and the cell configuration (nx,ny,nz), will determine if the two
// cell indices are within a 3x3 grid of each other.
//---------------------------------------------------------------------------------------------------------

static int Cell_Adjacent
(
   const int c1,                               ///< encoded index of cell 1
   const int c2,                               ///< encoded index of cell 2
   const int nx, const int ny, const int nz    ///< cell configuration <xyz>
)
{
   if (c1==c2)  return(1);       // (trivial case)

   int nxy = (nx*ny);

   int k1  = (int)(c1/nxy);
   int k2  = (int)(c2/nxy);
   int dk  = abs(k2-k1);
       dk  = ( (dk<(nz-1)) ? dk : 1 );    if (dk>1) return(0);

   int j1  = (int)(c1%nxy)/nx;
   int j2  = (int)(c2%nxy)/nx;
   int dj  = abs(j2-j1);
       dj  = ( (dj<(ny-1)) ? dj : 1 );    if (dj>1) return(0);

   int i1  = (int)(c1%nxy)%nx;
   int i2  = (int)(c2%nxy)%nx;
   int di  = abs(i2-i1);
       di  = ( (di<(ny-1)) ? di : 1 );    if (di>1) return(0);

   return (1);
}

// Segment_Pairs()
//
// Given an segment list, will return the segment pair matchups that fall within the given domain
// and cell arrangement.  Note that prior to calling this routine, the segment list should be
// cropped to only those segments that fall within the cells encapsulated by the domain.  If not - the
// n-squared nature of the loops included below can get horrifically long.
//-------------------------------------------------------------------------------------------------------------------

unsigned short *Segment_Pairs
(
         int           & npairs ,                     ///< number of segment pairs (returned)
         Segment_t      *segs   ,                     ///< source segment array
   const int             scnt   ,                     ///< number of segments in the source segment array
   const real8  x0, const real8 y0, const real8 z0,   ///< lower bound of domain
   const real8  x1, const real8 y1, const real8 z1,   ///< upper bound of domain
   const int    nx, const int   ny, const int   nz,   ///< number of cells <xyz>
   const real8  lx, const real8 ly, const real8 lz,   ///< simulation size <xyz>
         real8 *msec                                  ///< time in millisecs (if non-null)
)
{
   CPU_Timer tmr;
   CPU_TIMER_START(msec,tmr);

   unsigned short *spairs=0;  npairs=0; int nmax=0;               // (set default return)

   if (segs && (scnt>0))
   {
      Segment_List_Set_Cells (segs,scnt, nx,ny,nz, lx,ly,lz);     // set cell indices    (side effect)
      Segment_List_Set_Native(segs,scnt, x0,y0,z0, x1,y1,z1);     // set native segments (side-effect)

      for (int i=0  ; (i<scnt); i++)  // N
      for (int j=i+1; (j<scnt); j++)  // N*(N-1)/2
      {
         // only add segment pairs where at least one of the segments is native
         // and the segments are located in the same or adjacent cells.

         if ( segs[i].native || segs[j].native )
         {
            int c1 = segs[i].cell_indx;   // c1 = cell index of segment 1
            int c2 = segs[j].cell_indx;   // c2 = cell index of segment 2

            if ( Cell_Adjacent(c1,c2,nx,ny,nz) )
               spairs = Append_Segment_Pair(spairs,npairs,nmax,i,j);
         }
      }
   }

   CPU_TIMER_STOP(msec,tmr);

   return(spairs);
}

// Segment_List_Check()
//
// Will create and check all the segments within the current domain and check that all the segments defined
// wthin the domain are all less than the size of the cells within the simulation.  This was written to
// diagnose the secondary ghost errors that were occasionally causing fatal errors within the simulation.
// The hypothesis was that either the reposition of the nodes or a topology event was causing segments to
// grow too large.
//
// Note that this routine will be expensive and impact performance. It routine should only be invoked
// as a diagnostic.
//-------------------------------------------------------------------------------------------------------------------

static inline real8 Min (const real8 x, const real8 y, const real8 z) { return ( (x<y) ? ( (x<z) ? x : z ) : ( (y<z) ? y : z ) ); }

int Segment_List_Check
(
   Segment_t   *segs,                                       ///< array  of segments
   const int    scnt,                                       ///< number of segments
   const real8  lx, const real8 ly, const real8 lz ,        ///< overall simulation size <xyz> (needed for PBC repositions)
   const real8  minseg,                                     ///< minimum allowed segment length (suggest (rann/2), 0.0 disables)
   const real8  maxseg,                                     ///< maximum allowed segment length (suggest min cell dimension)
   const char  *msg                                         ///< additional message text to be included in warning
)
{
   int nfail=0;

   if (segs && (scnt>0))
   {
      real8 sx = ( (fabs(lx)>0.0) ? 1.0/fabs(lx) : 0.0 );   // sx = reciprocal of simulation size (x) (pbc always on)
      real8 sy = ( (fabs(ly)>0.0) ? 1.0/fabs(ly) : 0.0 );   // sy = reciprocal of simulation size (y) (pbc always on)
      real8 sz = ( (fabs(lz)>0.0) ? 1.0/fabs(lz) : 0.0 );   // sz = reciprocal of simulation size (z) (pbc always on)

      for (int i=0; (i<scnt); ++i)                          // check all segments
      {
         real8 x1[3] = { segs[i].p1[0],                     // x1 = local copy of node 1 position
                         segs[i].p1[1],                     //
                         segs[i].p1[2] };                   //

         real8 x2[3] = { segs[i].p2[0],                     // x2 = local copy of node 2 position
                         segs[i].p2[1],                     //
                         segs[i].p2[2] };                   //

         ::PBC_Position(x1,x2, lx,ly,lz,  sx,sy,sz);        // reposition x2 to nearest image

         real8 s = fabs(V3_DIST(x1,x2));                    // s = segment length

         if ( (s<minseg) || (s>maxseg) )                    // if the segment length is too short (or too long)
         {                                                  //
            Tag_t t1 = segs[i].node1->myTag;                //    t1 = tag  1
            Tag_t t2 = segs[i].node2->myTag;                //    t2 = tag  2

            real8 dx = x2[0]-x1[0];
            real8 dy = x2[1]-x1[1];
            real8 dz = x2[2]-x1[2];

            if (s<minseg)
            {  printf("WARNING %s : short segment detected (%d,%d)->(%d,%d) s=[%0.2lf %0.2lf %0.2lf] slen=%0.2lf smin=%0.2lf\n",
                         (msg ? msg : ""), t1.domainID, t1.index, t2.domainID, t2.index, dx, dy, dz, s, minseg ); }

            if (s>maxseg)
            {  printf("WARNING %s : long  segment detected (%d,%d)->(%d,%d) s=[%0.2lf %0.2lf %0.2lf] slen=%0.2lf smax=%0.2lf\n",
                         (msg ? msg : ""), t1.domainID, t1.index, t2.domainID, t2.index, dx, dy, dz, s, maxseg ); }

            nfail++;
         }
      }
   }

   return(nfail);
}

int Segment_List_Check (
   const Node_t **native     ,                           ///< sparse array of native node pointers
   const int      native_cnt ,                           ///< number of native node pointers
   const Node_t **ghosts     ,                           ///< sparse array of ghost node pointers  (NIL okay)
   const int      ghost_cnt  ,                           ///< number of ghost node pointers        (NIL okay)
   const real8  lx, const real8 ly, const real8 lz ,     ///< overall simulation size       <xyz>
   const int    nx, const int   ny, const int   nz ,     ///< simulation cell configuration <xyz>
         real8  minseg,                                  ///< minimum allowed segment length (default=0.0, suggest (rann/2))
         real8  maxseg,                                  ///< maximum allowed segment length (default=shortest cell dimension)
   const char  *msg                                      ///< additional message text to be included in warning
)
{
   int        scnt = 0;
   Segment_t *segs = Segment_List_Create(scnt, native, native_cnt, ghosts, ghost_cnt);

   real8 cx = ( (nx>0) ? (lx/nx) : 0.0 );                // cx = cell dimension (x)
   real8 cy = ( (ny>0) ? (ly/ny) : 0.0 );                // cy = cell dimension (y)
   real8 cz = ( (nz>0) ? (lz/nz) : 0.0 );                // cz = cell dimension (z)

   real8 cmin   = Min(cx,cy,cz);                         // cmin = shortest cell dimension

         maxseg = ( (maxseg>0.0) ? maxseg : cmin );      // if maxseg is not provided, use shortest cell dimension

   int nfail = Segment_List_Check(segs,scnt, lx,ly,lz, minseg, maxseg, msg );

   if (segs) { delete [] segs; segs=0; }

   return(nfail);
}
