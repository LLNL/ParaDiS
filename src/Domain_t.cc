
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "Domain_t.h"

#define DELETE(a) if (a) { delete [] a; a=0; }

Domain_t::Domain_t (void)
{
   dndx=0;
   x0=x1=0.0;
   y0=y1=0.0;
   z0=z1=0.0;
}

Domain_t::Domain_t (const Domain_t & dom)
{
   dndx = dom.dndx;
   x0   = dom.x0;
   x1   = dom.x1;
   y0   = dom.y0;
   y1   = dom.y1;
   z0   = dom.z0;
   z1   = dom.z1;
}

Domain_t::Domain_t
(
   const int    i,
   const real8 *p1,
   const real8 *p2
)
{
   dndx=i;
   if (p1) { x0=p1[0]; y0=p1[1]; z0=p1[2]; } else { x0=y0=z0=0.0; }
   if (p2) { x1=p2[0]; y1=p2[1]; z1=p2[2]; } else { x1=y1=z1=0.0; }
}

Domain_t::Domain_t
(
   const real8 *p1,
   const real8 *p2
)
{
   dndx=0;
   if (p1) { x0=p1[0]; y0=p1[1]; z0=p1[2]; } else { x0=y0=z0=0.0; }
   if (p2) { x1=p2[0]; y1=p2[1]; z1=p2[2]; } else { x1=y1=z1=0.0; }
}

Domain_t::Domain_t
(
   const int   i,
   const real8 dx0, const real8 dx1,
   const real8 dy0, const real8 dy1,
   const real8 dz0, const real8 dz1
)
{
   dndx=i;
   x0=dx0; x1=dx1;
   y0=dy0; y1=dy1;
   z0=dz0; z1=dz1;
}

Domain_t::Domain_t
(
   const real8 dx0, const real8 dx1,
   const real8 dy0, const real8 dy1,
   const real8 dz0, const real8 dz1
)
{
   dndx=0;
   x0=dx0; x1=dx1;
   y0=dy0; y1=dy1;
   z0=dz0; z1=dz1;
}

const Domain_t & Domain_t::operator = (const Domain_t & dom)
{
   dndx = dom.dndx;
   x0   = dom.x0;
   x1   = dom.x1;
   y0   = dom.y0;
   y1   = dom.y1;
   z0   = dom.z0;
   z1   = dom.z1;

   return(*this);
}

const Domain_t & Domain_t::operator = (const RBDecomp_t & dom)
{
   dndx = dom.domID;
   x0   = dom.cMin[0];
   x1   = dom.cMax[0];
   y0   = dom.cMin[1];
   y1   = dom.cMax[1];
   z0   = dom.cMin[2];
   z1   = dom.cMax[2];

   return(*this);
}

int Domain_t::operator == (const Domain_t & dom) const
{
   real8 eps=1.0e-4;

   return (    (fabs(x0-dom.x0)>eps)
            || (fabs(x1-dom.x1)>eps)
            || (fabs(y0-dom.y0)>eps)
            || (fabs(y1-dom.y1)>eps)
            || (fabs(z0-dom.z0)>eps)
            || (fabs(z1-dom.z1)>eps) ? 0 : 1 );
}

// Center()
//
// Returns the position of the domain center (xyz).
//------------------------------------------------------------------------------------------------------------

void Domain_t::Center (real8 *c) const
{
   if (c)
   {
      c[0] = (x0+x1)/2.0;
      c[1] = (y0+y1)/2.0;
      c[2] = (z0+z1)/2.0;
   }
}

void Domain_t::Center (real8 &cx, real8 &cy, real8 &cz) const
{
   cx = (x0+x1)/2.0;
   cy = (y0+y1)/2.0;
   cz = (z0+z1)/2.0;
}

// Intersect()
//
// Returns true/false (1/0) on whether two domains overlap each other.
//------------------------------------------------------------------------------------------------------------

int Domain_t::Intersect (const Domain_t & dom) const
{
   return (    (x1<dom.x0) || (dom.x1<x0)
            || (y1<dom.y0) || (dom.y1<y0)
            || (z1<dom.z0) || (dom.z1<z0) ? 0 : 1 );
}

// Inside()
//
// Returns true/false (1/0) on whether two domains overlap each other.
//------------------------------------------------------------------------------------------------------------

int Domain_t::Inside (const real8 *v) const
{
   return ( v && (    (v[0]<x0) || (x1<v[0])
                   || (v[1]<y0) || (x1<v[1])
                   || (v[2]<z0) || (x1<v[2]) ) ? 0 : 1 );
}

int Domain_t::Inside (const real8 x, const real8 y, const real8 z) const
{
   return (    (x<x0) || (x1<x)
            || (y<y0) || (y1<y)
            || (z<z0) || (z1<z) ? 0 : 1 );
}

// Outside()
//
// Returns true/false (1/0) on whether two domains overlap each other.
//------------------------------------------------------------------------------------------------------------

int Domain_t::Outside (const real8 *v) const
{
   return ( v && (    (v[0]<x0) || (x1<v[0])
                   || (v[1]<y0) || (y1<v[1])
                   || (v[2]<z0) || (z1<v[2]) ) ? 1 : 0 );
}

int Domain_t::Outside (const real8 x, const real8 y, const real8 z) const
{
   return (    (x<x0) || (x1<x)
            || (y<y0) || (y1<y)
            || (z<z0) || (z1<z) ? 1 : 0 );
}

// Volume()
//
// Returns cubic volume of domain.
//------------------------------------------------------------------------------------------------------------

real8 Domain_t::Volume (void) const
{
   real8 dx = (x1-x0);
   real8 dy = (y1-y0);
   real8 dz = (z1-z0);

   return(dx*dy*dz);
}

// Expand()
//
// Will expand the size of the current domain to include the physical space of the cells that encapsulate
// the current domain. This method is commonly used to identify which MPI domains we need to communicate
// with. Note that this can result in a domain that can extend beyond the periodic boundaries of the current
// simulation.
//
// By default - will expand the domain to include just those cells that are included within the domain.
// Additional buffer cells can be obtained using the na argument. For example - if you want to include
// the 3x3 domain of cells that are included by the current domain, set nadd=1.
//------------------------------------------------------------------------------------------------------------

void Domain_t::Expand
(
   const real8 lx ,  ///< overall simulation size <x>
   const real8 ly ,  ///< overall simulation size <y>
   const real8 lz ,  ///< overall simulation size <z>
   const int   nx ,  ///< number of cells <x>
   const int   ny ,  ///< number of cells <y>
   const int   nz ,  ///< number of cells <z>
   const int   nadd  ///< number of additional cells on each boundary
)
{
   const real8 lx2=(lx/2),  dx=(lx/nx);
   const real8 ly2=(ly/2),  dy=(ly/ny);
   const real8 lz2=(lz/2),  dz=(lz/nz);

   x0 = (floor((x0+lx2)/dx)-nadd)*dx-lx2;
   x1 = (ceil ((x1+lx2)/dx)+nadd)*dx-lx2;
   y0 = (floor((y0+ly2)/dy)-nadd)*dy-ly2;
   y1 = (ceil ((y1+ly2)/dy)+nadd)*dy-ly2;
   z0 = (floor((z0+lz2)/dz)-nadd)*dz-lz2;
   z1 = (ceil ((z1+lz2)/dz)+nadd)*dz-lz2;
}

// Native_Cells()
//
// Given a simulation and cell geometry, will return an array of integers corresponding to the cell
// indices that intersect the current domain.
//
// Note assumes the current simulation space is [-lx/2:+lx/2].
//------------------------------------------------------------------------------------------------------------

#define Constrain(a,b,c) ( ((a)<(b)) ? (b) : ( ((a)<(c)) ? (a) : (c) ) )

void Domain_t::Native_Cells
(
   int & ix0, int & ix1,  ///< resulting cell min/max <x>
   int & iy0, int & iy1,  ///< resulting cell min/max <y>
   int & iz0, int & iz1,  ///< resulting cell min/max <z>
   const real8 lx ,       ///< simulation size <x>
   const real8 ly ,       ///< simulation size <y>
   const real8 lz ,       ///< simulation size <z>
   const int   nx ,       ///< number of cells <x>
   const int   ny ,       ///< number of cells <y>
   const int   nz         ///< number of cells <z>
)
{
   const real8 dx=(lx/nx);                  // dx  = cell size <x>
   const real8 dy=(ly/ny);                  // dy  = cell size <y>
   const real8 dz=(lz/nz);                  // dz  = cell size <z>

   ix0 = (int)((x0+(lx/2.0)   )/dx);        // ix0 = lower left cell index <x>
   iy0 = (int)((y0+(ly/2.0)   )/dy);        // iy0 = lower left cell index <y>
   iz0 = (int)((z0+(lz/2.0)   )/dz);        // iz0 = lower left cell index <z>

   ix1 = (int)((x1+(lx/2.0)-dx)/dx);        // ix1 = upper right cell index <x>
   iy1 = (int)((y1+(ly/2.0)-dy)/dy);        // iy1 = upper right cell index <y>
   iz1 = (int)((z1+(lz/2.0)-dz)/dz);        // iz1 = upper right cell index <z>

   if ( (dx*(ix1+1))-(lx/2.0) < x1 ) ix1++; // adjust for odd cell sizes <x>
   if ( (dy*(iy1+1))-(ly/2.0) < y1 ) iy1++; // adjust for odd cell sizes <y>
   if ( (dz*(iz1+1))-(lz/2.0) < z1 ) iz1++; // adjust for odd cell sizes <z>

   ix0 = Constrain(ix0, 0 ,(nx-1));         // constrain to valid indices
   ix1 = Constrain(ix1,ix0,(nx-1));         //
   iy0 = Constrain(iy0, 0 ,(ny-1));         //
   iy1 = Constrain(iy1,iy0,(ny-1));         //
   iz0 = Constrain(iz0, 0 ,(nz-1));         //
   iz1 = Constrain(iz1,iz0,(nz-1));         //
}

// Native_Cells()
//
// Given a simulation and cell geometry, will return an array of integers corresponding to the cell
// indices that intersect the current domain. Note that for native cells, there are no PBC issues to
// resolve. The returned index list will be sorted with gaps with the lowest cell index first.
// The number of native cells is also returned.
//------------------------------------------------------------------------------------------------------------

int *Domain_t::Native_Cells
(
         int   & ccnt,  ///< number of native cells (returned)
   const real8   lx  ,  ///< simulation size <x>
   const real8   ly  ,  ///< simulation size <y>
   const real8   lz  ,  ///< simulation size <z>
   const int     nx  ,  ///< number of cells <x>
   const int     ny  ,  ///< number of cells <y>
   const int     nz     ///< number of cells <z>
)
{
   int ix0=0, ix1=0;
   int iy0=0, iy1=0;
   int iz0=0, iz1=0;

   Native_Cells(ix0,ix1, iy0,iy1, iz0,iz1, lx,ly,lz, nx,ny,nz);

   int ncx = (ix1-ix0+1);
   int ncy = (iy1-iy0+1);
   int ncz = (iz1-iz0+1);

        ccnt  = (ncx*ncy*ncz);
   int *cells = ( (ccnt>0) ? new int [ccnt] : 0 );

   if (cells)
   {
      int m=0;
      for (int k=iz0; (k<=iz1); k++)
      for (int j=iy0; (j<=iy1); j++)
      for (int i=ix0; (i<=ix1); i++)
         cells[m++] = (k*nx*ny)+(j*nx)+i;
   }

   return(cells);
}

// Remote_Cells()
//
// Returns an array of cell indices that include the native cells and a 3x3 padded ghost block surrounding
// the current domain. Applies the periodic boundaries as appropriate.
//------------------------------------------------------------------------------------------------------------

int *Domain_t::Remote_Cells
(
         int   & ccnt,  ///< number of remote cells
   const real8   lx  ,  ///< simulation size <x>
   const real8   ly  ,  ///< simulation size <y>
   const real8   lz  ,  ///< simulation size <z>
   const int     nx  ,  ///< number of cells <x>
   const int     ny  ,  ///< number of cells <y>
   const int     nz  ,  ///< number of cells <z>
   const int     pbx ,  ///< periodic boundary active (1=yes,0=no) <x>
   const int     pby ,  ///< periodic boundary active (1=yes,0=no) <y>
   const int     pbz ,  ///< periodic boundary active (1=yes,0=no) <z>
   const int     rx  ,  ///< remote padding  <x>
   const int     ry  ,  ///< remote padding  <y>
   const int     rz     ///< remote padding  <z>
)
{
   int ix0=0, ix1=0;
   int iy0=0, iy1=0;
   int iz0=0, iz1=0;

   Native_Cells(ix0,ix1, iy0,iy1, iz0,iz1, lx,ly,lz, nx,ny,nz);

   // expand the region...
   // for the 3x3x3 region of ghost cells, set rx,ry,rz==1
   // for the 4x4x4 region of ghost cells, set rx,ry,rz==2

   ix0-=rx; ix1+=rx;
   iy0-=ry; iy1+=ry;
   iz0-=rz; iz1+=rz;

   // apply periodic boundaries...

   if (!pbx) { ix0=Constrain(ix0, 0 ,(nx-1));
               ix1=Constrain(ix1,ix0,(nx-1)); }

   if (!pby) { iy0=Constrain(iy0, 0 ,(ny-1));
               iy1=Constrain(iy1,iy0,(ny-1)); }

   if (!pbz) { iz0=Constrain(iz0, 0 ,(nz-1));
               iz1=Constrain(iz1,iz0,(nz-1)); }

   int  ncx    = (ix1-ix0+1);
   int  ncy    = (iy1-iy0+1);
   int  ncz    = (iz1-iz0+1);

        ccnt  = (ncx*ncy*ncz);
   int *cells = ( (ccnt>0) ? new int [ccnt] : 0 );

   if (cells)
   {
      int m=0;
      for (int iz=iz0; (iz<=iz1); iz++)
      for (int iy=iy0; (iy<=iy1); iy++)
      for (int ix=ix0; (ix<=ix1); ix++)
      {
         int i = ( (ix<0) ? (ix+nx) : (ix%nx) );
         int j = ( (iy<0) ? (iy+ny) : (iy%ny) );
         int k = ( (iz<0) ? (iz+nz) : (iz%nz) );

         cells[m++] = (k*nx*ny)+(j*nx)+i;
      }
   }

   return(cells);
}

// Print()
//
// Prints the current domain to an output stream.
//------------------------------------------------------------------------------------------------------------

void Domain_t::Print  (FILE *fd) const
{
   if (fd)
   {
      fprintf(fd,"%4d  %11.4f %11.4f %11.4f  %11.4f %11.4f %11.4f\n", dndx, x0, x1, y0, y1, z0, z1 );
   }
}

// Print_GNU()
//
// Prints the domain as a series of gnuplot vectors to an output stream.
//------------------------------------------------------------------------------------------------------------

void Domain_t::Print_GNU  (FILE *fd) const
{
   if (fd)
   {
      real8 dx = (x1-x0);
      real8 dy = (y1-y0);
      real8 dz = (z1-z0);

      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x0, y0, z0, dx , 0.0, 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x0, y1, z0, dx , 0.0, 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x0, y0, z1, dx , 0.0, 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x0, y1, z1, dx , 0.0, 0.0 );

      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x0, y0, z0, 0.0, dy , 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x1, y0, z0, 0.0, dy , 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x0, y0, z1, 0.0, dy , 0.0 );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x1, y0, z1, 0.0, dy , 0.0 );

      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x0, y0, z0, 0.0, 0.0, dz  );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x1, y0, z0, 0.0, 0.0, dz  );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x0, y1, z0, 0.0, 0.0, dz  );
      fprintf(fd,"%4d %11.4lf %11.4lf %11.4lf  %11.4lf %11.4lf %11.4lf\n", dndx, x1, y1, z0, 0.0, 0.0, dz  );

      fprintf(fd,"\n");
   }
}

// Append()
//
// Appends a domain to an existing array of domains. Returns the vector of domains.
//------------------------------------------------------------------------------------------------------------

Domain_t *Append
(
   Domain_t  *doms ,  ///< array of domains (allocated and extended if null)
   int      & cnt  ,  ///< current count of domains
   int      & max  ,  ///< current size of the domain array
   Domain_t & dom     ///< domain to append
)
{
   if ( !doms || (cnt==max) )
   {
      max = ( (max==0) ? 2000 : 2*max );

      Domain_t *tmp = new Domain_t[max];

      if (tmp && doms) { for (int i=0; (i<cnt); i++) tmp[i]=doms[i]; }

      if (doms) { delete [] doms; }

      doms = tmp;
   }

   if (doms && (cnt<max) ) { doms[cnt++] = dom; }

   return(doms);
}

// Create_Domains()
//
// Will create and return an array of domains for a given simulation size and cell configuration.
//------------------------------------------------------------------------------------------------------------

Domain_t *Create_Domains
(
   const real8 lx,   ///< simulation size <x>
   const real8 ly,   ///< simulation size <y>
   const real8 lz,   ///< simulation size <z>
   const int   nx,   ///< cell count <x>
   const int   ny,   ///< cell count <y>
   const int   nz    ///< cell count <z>
)
{
   real8 x0 = (-lx/2.0);
   real8 y0 = (-ly/2.0);
   real8 z0 = (-lz/2.0);

   real8 dx = ( (nx>0) ? (lx/nx) : 0.0 );
   real8 dy = ( (ny>0) ? (ly/ny) : 0.0 );
   real8 dz = ( (nz>0) ? (lz/nz) : 0.0 );

   Domain_t *doms = ( ((nx*ny*nz)>0) ? new Domain_t[nx*ny*nz] : 0 );

   real8 z=z0;
   for (int k=0,m=0; (k<nz); k++, z+=dz)
   {
      real8 y=y0;
      for (int j=0; (j<ny); j++, y+=dy)
      {
         real8 x=x0;
         for (int i=0; (i<nx); i++, m++, x+=dx)
            doms[m] = Domain_t(m, x,(x+dx), y,(y+dy), z,(z+dz) );
      }
   }

   return(doms);
}

// Create_Domains()
//
// Will return an array of domains initialized from an RS domain decomposition.
//------------------------------------------------------------------------------------------------------------

Domain_t *Create_Domains
(
   const RSDecomp_t *rs,   ///< points to an RS domain descriptor
   const int         nx,   ///< number of domains <x>
   const int         ny,   ///< number of domains <y>
   const int         nz    ///< number of domains <z>
)
{
   int         ndom = (nx*ny*nz);
   Domain_t   *doms = ( (ndom>0) ? new Domain_t[ndom] : 0 );

   real8      *rsx  = ( rs ? rs->domBoundX : 0 );
   real8     **rsy  = ( rs ? rs->domBoundY : 0 );
   real8    ***rsz  = ( rs ? rs->domBoundZ : 0 );

   if (doms && rsx && rsy && rsz)
   {
      int dndx=0;
      for (int i=0; (i<nx); i++)
      for (int j=0; (j<ny); j++)
      for (int k=0; (k<nz); k++, dndx++)
      {
         doms[dndx] = Domain_t( dndx,
                                rsx[i]      , rsx[i+1],
                                rsy[i][j]   , rsy[i  ][j+1],
                                rsz[i][j][k], rsz[i  ][j  ][k+1] );
      }
   }

   return(doms);
}

// Create_Domains()
//
// Will return an array of domains initialized from an RB domain decomposition.
//------------------------------------------------------------------------------------------------------------

Domain_t *Create_Domains
(
   const RBDecomp_t *rb,   ///< points to array of RB domains
   const int         nx,   ///< number of domains <x>
   const int         ny,   ///< number of domains <y>
   const int         nz    ///< number of domains <z>
)
{
   int       ndom = (nx*ny*nz);
   Domain_t *doms = ( (ndom>0) ? new Domain_t[ndom] : 0 );

   if (doms && rb)
   {
      for (int i=0,j=0; (j<ndom); i++)
         if (rb[i].domID>=0) { doms[j++] = rb[i]; }
   }

   return(doms);
}

