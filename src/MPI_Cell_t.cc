
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "Domain_t.h"
#include "MPI_Cell_t.h"

//------------------------------------------------------------------------------------------------------------

#define VDELETE(a)     if (a) { delete [] a; a=0; }

// (constructors)
//------------------------------------------------------------------------------------------------------------

MPI_Cell_t::MPI_Cell_t (void) 
{ 
   cndx=0;
   x0=x1=0.0; 
   y0=y1=0.0; 
   z0=z1=0.0; 

   ncnt=nmax=0; nodes=0;
   neighbors=0;
}

MPI_Cell_t::MPI_Cell_t (const MPI_Cell_t & cell)
{
   cndx = cell.cndx;
   x0   = cell.x0;
   x1   = cell.x1;
   y0   = cell.y0;
   y1   = cell.y1;
   z0   = cell.z0;
   z1   = cell.z1;

   ncnt=nmax=0; nodes=0;
   neighbors=0;
}

MPI_Cell_t::MPI_Cell_t 
(
   const int   _cndx, 
   const real8 _x0, const real8 _x1, 
   const real8 _y0, const real8 _y1, 
   const real8 _z0, const real8 _z1 
)
{
   cndx = _cndx;
   x0   = _x0;
   x1   = _x1;
   y0   = _y0;
   y1   = _y1;
   z0   = _z0;
   z1   = _z1;

   ncnt=nmax=0; nodes=0;
   neighbors=0;
}

// (destructor)
//------------------------------------------------------------------------------------------------------------

MPI_Cell_t::~MPI_Cell_t()
{
   Recycle();
}

void MPI_Cell_t::Recycle (void)
{
   VDELETE(nodes    );
   VDELETE(neighbors);

   ncnt=nmax=0; nodes=0;
   neighbors=0;
}

// (assignment)
//------------------------------------------------------------------------------------------------------------

const MPI_Cell_t & MPI_Cell_t::operator = (const MPI_Cell_t & c)
{
   if (this != &c)
   {
      cndx = c.cndx;
      x0   = c.x0;
      x1   = c.x1;
      y0   = c.y0;
      y1   = c.y1;
      z0   = c.z0;
      z1   = c.z1;

      VDELETE(nodes); ncnt=nmax=0; nodes=0;

      nodes = ( (c.ncnt>0) ? new Node_t *[c.ncnt] : 0 );

      ncnt=nmax= ( nodes ? c.ncnt : 0 );

      if (nodes && c.nodes) { memcpy(nodes, c.nodes, ncnt*sizeof(Node_t *)); }

      VDELETE(neighbors);
   }

   return(*this);
}

// Set_Neighbors()
//
// Will populate the array of neighbor indices based on the cell configuration provided.
//------------------------------------------------------------------------------------------------------------

void MPI_Cell_t::Set_Neighbors 
(
   const int   nx   , const int   ny  , const int   nz,   ///< nx,ny,nz    = number of cells (x,y,z)
   const int   pbx  , const int   pby , const int   pbz   ///< pbx,pby,pbz = periodic boundaries active (x,y,z)
)
{
   VDELETE(neighbors);

   neighbors = new int[3*3*3];

   if (neighbors)
   { for (int i=0; (i<(3*3*3)); i++) { neighbors[i]=-1; } }

   if (neighbors && (nx>0) && (ny>0) && (nz>0) )
   {
      for (int k=0,m=0; (k<3); k++)
      for (int j=0;     (j<3); j++)
      for (int i=0;     (i<3); i++, m++)
      {
         int indx = cndx;

         // convert cell index to (ix,iy,iz) index triple

         int iz = indx/(nx*ny);  indx%=(nx*ny);
         int iy = indx/(nx   );
         int ix = indx%(nx   );                   

         ix+=(i-1);
         iy+=(j-1);
         iz+=(k-1);

         // if we are outside of the simulation box and PBC not active...

         if ( ((ix<0) || (ix>=nx)) && (!pbx) ) { neighbors[m]=-1; continue; }
         if ( ((iy<0) || (iy>=ny)) && (!pby) ) { neighbors[m]=-1; continue; }
         if ( ((iz<0) || (iz>=nz)) && (!pbz) ) { neighbors[m]=-1; continue; }

         // otherwise adjust for periodic boundaries...

         ix += ( (ix<0) ? nx : 0 );  ix%=nx;
         iy += ( (iy<0) ? ny : 0 );  iy%=ny;
         iz += ( (iz<0) ? nz : 0 );  iz%=nz;

         indx = (iz*nx*ny)+(iy*nx)+ix;

         if ( (indx<0) || (indx>=(nx*ny*nz)) ) 
         { printf("%s::%s(ln=%d) cell indexing error\n", __FILE__, __func__, __LINE__ ); }

         neighbors[m] = ( (indx==cndx) ? -1 : indx ); 
      }
   }
}

// Intersect()
//
// Returns true/false (1/0) on whether two cells overlap each other. 
//------------------------------------------------------------------------------------------------------------

int MPI_Cell_t::Intersect (const MPI_Cell_t & cell) const
{
   return (    (x1<cell.x0) || (cell.x1<x0)
            || (y1<cell.y0) || (cell.y1<y0)
            || (z1<cell.z0) || (cell.z1<z0) ? 0 : 1 );
}

int MPI_Cell_t::Intersect (const Domain_t & dom) const
{
   return (    (x1<dom.x0) || (dom.x1<x0)
            || (y1<dom.y0) || (dom.y1<y0)
            || (z1<dom.z0) || (dom.z1<z0) ? 0 : 1 );
}

// Inside()
//
// Returns true/false (1/0) on whether a given position falls within a cell.
//------------------------------------------------------------------------------------------------------------

int MPI_Cell_t::Inside (const real8 *v) const
{
   return ( v && (    (v[0]<x0) || (x1<v[0]) 
                   || (v[1]<y0) || (x1<v[1]) 
                   || (v[2]<z0) || (x1<v[2]) ) ? 0 : 1 );
}

int MPI_Cell_t::Inside (const real8 x, const real8 y, const real8 z) const
{
   return (    (x<x0) || (x1<x) 
            || (y<y0) || (y1<y) 
            || (z<z0) || (z1<z) ? 0 : 1 );
}

// Outside()
//
// Returns true/false (1/0) on whether a given position falls outside a cell.
//------------------------------------------------------------------------------------------------------------

int MPI_Cell_t::Outside (const real8 *v) const
{
   return ( v && (    (v[0]<x0) || (x1<v[0]) 
                   || (v[1]<y0) || (y1<v[1]) 
                   || (v[2]<z0) || (z1<v[2]) ) ? 1 : 0 );
}

int MPI_Cell_t::Outside (const real8 x, const real8 y, const real8 z) const
{
   return (    (x<x0) || (x1<x) 
            || (y<y0) || (y1<y) 
            || (z<z0) || (z1<z) ? 1 : 0 );
}

// Append()
//
// Appends the node to the node pointer array.
//------------------------------------------------------------------------------------------------------------

void MPI_Cell_t::Append (const Node_t *node)
{
   if ( !nodes || (ncnt==nmax) )
   {
      nmax = ( (nmax==0) ? 32 : 2*nmax );

      Node_t **tmp = new Node_t *[nmax];

      if (tmp && nodes) { for (int i=0; (i<ncnt); i++) tmp[i]=nodes[i]; }

      if (nodes) { delete [] nodes; }

      nodes = tmp;
   }
   
   if (nodes && (ncnt<nmax) ) { nodes[ncnt++] = (Node_t *) node; }
}

// Print()
//
// Prints the current cell to the output file...
//------------------------------------------------------------------------------------------------------------

void MPI_Cell_t::Print (FILE *fd) const
{
   if (fd)
   {
      fprintf(fd,"%4d  %11.4f %11.4f %11.4f  %11.4f %11.4f %11.4f %3d\n", cndx, x0, x1, y0, y1, z0, z1, ncnt );
   }
}

// Print_GNU()
//
// Prints the cells in a gnuplot compatable format
//------------------------------------------------------------------------------------------------------------

void MPI_Cell_t::Print_GNU (FILE *fd) const
{
   if (fd)
   {
      real8 dx = (x1-x0);
      real8 dy = (y1-y0);
      real8 dz = (z1-z0);

      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x0, y0, z0, dx , 0.0, 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x0, y1, z0, dx , 0.0, 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x0, y0, z1, dx , 0.0, 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x0, y1, z1, dx , 0.0, 0.0 );

      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x0, y0, z0, 0.0, dy , 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x1, y0, z0, 0.0, dy , 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x0, y0, z1, 0.0, dy , 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x1, y0, z1, 0.0, dy , 0.0 );

      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x0, y0, z0, 0.0, 0.0, dz  );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x1, y0, z0, 0.0, 0.0, dz  );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x0, y1, z0, 0.0, 0.0, dz  );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", cndx, x1, y1, z0, 0.0, 0.0, dz  );

      fprintf(fd,"\n");
   }
}

// MPI_Cells_Append()
//
// Appends a cell to an existing array of cells. Returns the vector of cells.
//------------------------------------------------------------------------------------------------------------

MPI_Cell_t *MPI_Cells_Append 
(
   MPI_Cell_t  *cells,  ///< array of cells (allocated and extended if null)
   int        & cnt  ,  ///< current count of cells
   int        & max  ,  ///< current size of the cell array
   MPI_Cell_t & cell    ///< cell to append
)
{
   if ( !cells || (cnt==max) )
   {
      max = ( (max==0) ? 32 : 2*max );

      MPI_Cell_t *tmp = new MPI_Cell_t[max];

      if (tmp && cells) { for (int i=0; (i<cnt); i++) tmp[i]=cells[i]; }

      if (cells) { delete [] cells; }

      cells = tmp;
   }
   
   if (cells && (cnt<max) ) { cells[cnt++] = cell; }

   return(cells);
}

// MPI_Cells_Create()
//
// Will allocate and initialize an array of MPI communication cells.
//------------------------------------------------------------------------------------------------------------

MPI_Cell_t *MPI_Cells_Create 
(
   const int   nx  ,  ///< number of cells (x)
   const int   ny  ,  ///< number of cells (y)
   const int   nz  ,  ///< number of cells (z)
   const real8 lx  ,  ///< simulation size (x)
   const real8 ly  ,  ///< simulation size (y)
   const real8 lz  ,  ///< simulation size (z)
   const int   pbx ,  ///< periodic boundary active (x, 1=yes,0=no)?
   const int   pby ,  ///< periodic boundary active (y, 1=yes,0=no)?
   const int   pbz    ///< periodic boundary active (z, 1=yes,0=no)?
)
{
   MPI_Cell_t *cells = ( ((nx*ny*nz)>0) ? new MPI_Cell_t[nx*ny*nz] : 0 );

   if (cells)
   {
      const real8 dx =  ( (nx>0) ? (lx/nx) : 0.0 );
      const real8 dy =  ( (ny>0) ? (ly/ny) : 0.0 );
      const real8 dz =  ( (nz>0) ? (lz/nz) : 0.0 );
      const real8 x0 = -(lx/2.0);
      const real8 y0 = -(ly/2.0);
      const real8 z0 = -(lz/2.0);

      int m=0;
      int i,j,k; real8 x,y,z;
      for (k=0, z=z0; (k<nz); k++, z+=dz)
      for (j=0, y=y0; (j<ny); j++, y+=dy)
      for (i=0, x=x0; (i<nx); i++, x+=dx, m++)
      { cells[m] = MPI_Cell_t(m, x,x+dx, y,y+dy, z,z+dz); }

      // init the neighbor index lists for each cell...

      for (int i=0; (i<(nx*ny*nz)); i++)
      { cells[i].Set_Neighbors(nx,ny,nz, pbx,pby,pbz); } 
   }

   return(cells);
}

// MPI_Cells_Bin()
//
// Given an array of cells and the simulation geometry, will bin an array of nodes into the cells.
//------------------------------------------------------------------------------------------------------------

void MPI_Cells_Bin 
(
         MPI_Cell_t *cells,  ///< array  of cells
         Node_t     *nodes,  ///< array  of nodes
   const int         ncnt ,  ///< number of nodes
   const int         nx   ,  ///< number of cells (x)
   const int         ny   ,  ///< number of cells (y)
   const int         nz   ,  ///< number of cells (z)
   const real8       lx   ,  ///< simulation size (x)
   const real8       ly   ,  ///< simulation size (y)
   const real8       lz      ///< simulation size (z)
)
{
   if (cells && nodes && (ncnt>0) && (nx>0) && (ny>0) && (nz>0) )
   {
      const real8 dx =  (lx/nx);
      const real8 dy =  (ly/ny);
      const real8 dz =  (lz/nz);
      const real8 x0 = -(lx/2.0);
      const real8 y0 = -(ly/2.0);
      const real8 z0 = -(lz/2.0);

      for (int i=0; (i<(nx*ny*nz)); i++ ) { cells[i].Recycle(); }

      for (int i=0; (i<ncnt); i++)
      {
         int ix = (nodes[i].x-x0) / dx;
         int iy = (nodes[i].y-y0) / dy;
         int iz = (nodes[i].z-z0) / dz;

             ix += ( (ix<0) ? nx : 0 );  ix%=nx;
             iy += ( (iy<0) ? ny : 0 );  iy%=ny;
             iz += ( (iz<0) ? nz : 0 );  iz%=nz;

         int indx = (iz*nx*ny)+(iy*nx)+ix;

         if ( (indx<0) || (indx>=ncnt) ) { printf("%s::%s(ln=%d) invalid cell index\n", __FILE__, __func__, __LINE__ ); }

         if ( (0<=indx) && (indx<ncnt) )
         { cells[indx].Append( &nodes[i] ); }
      }
   }
}

