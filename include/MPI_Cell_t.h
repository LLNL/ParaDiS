#pragma once

#ifndef _PDS_MPI_CELL_H
#define _PDS_MPI_CELL_H

#include "mpi_portability.h"
#include "cuda_portability.h"
#include "Typedefs.h"
#include "Domain_t.h"
#include "Node.h"

//------------------------------------------------------------------------------------------------------------

class MPI_Cell_t
{
   public : 
      int          cndx ;      ///< cell index
      real8        x0,x1;      ///< cell min/max <x>
      real8        y0,y1;      ///< cell min/max <y>
      real8        z0,z1;      ///< cell min/max <z>

      int          ncnt ;      ///< current   node count
      int          nmax ;      ///< allocated node count
      Node_t     **nodes;      ///< array of node pointers that fall within this cell

      int         *neighbors;  ///< (3x3x3) array of cell neighbor indices 

   public : 
      MPI_Cell_t (void);
      MPI_Cell_t (const int   _cndx, 
                  const real8 _x0, const real8 _x1, 
                  const real8 _y0, const real8 _y1, 
                  const real8 _z0, const real8 _z1 );
      MPI_Cell_t (const MPI_Cell_t & cell);

     ~MPI_Cell_t();

     void  Recycle     (void);

     const MPI_Cell_t & operator = (const MPI_Cell_t & cell);

     void  Set_Neighbors (const int   nx   , const int   ny  , const int   nz,
                          const int   pbx  , const int   pby , const int   pbz );

     int   Intersect     (const MPI_Cell_t & cell) const;
     int   Intersect     (const Domain_t   & dom ) const;

     int   Inside        (const real8 *v) const;
     int   Inside        (const real8  x, const real8 y, const real8 z) const;

     int   Outside       (const real8 *v) const;
     int   Outside       (const real8  x, const real8 y, const real8 z) const;

     void  Append        (const Node_t *node);

     void  Print         (FILE *fd=stdout) const;
     void  Print_GNU     (FILE *fd=stdout) const;
};

//------------------------------------------------------------------------------------------------------------

MPI_Cell_t *MPI_Cells_Append        (MPI_Cell_t *cells, int & dcnt, int & dmax, MPI_Cell_t & cell);

MPI_Cell_t *MPI_Cells_Create        (const int   nx , const int   ny , const int   nz ,
                                     const real8 lx , const real8 ly , const real8 lz , 
                                     const int   pbx, const int   pby, const int   pbz );

void        MPI_Cells_Bin           (MPI_Cell_t *cells,
                                     Node_t     *nodes, const int   ncnt,
                                     const int   nx   , const int   ny  , const int   nz,
                                     const real8 lx   , const real8 ly  , const real8 lz );

#endif  //  _PDS_MPI_CELL_H
