#pragma once

#ifndef _PDS_DOMAIN_H
#define _PDS_DOMAIN_H

#include "mpi_portability.h"

#include "Typedefs.h"
#include "Decomp.h"

class Domain_t
{
   public :
      int    dndx     ;  ///< domain index
      real8  x0,x1    ;  ///< domain min/max <x>
      real8  y0,y1    ;  ///< domain min/max <y>
      real8  z0,z1    ;  ///< domain min/max <z>

   public :
      Domain_t (void);

      Domain_t (const Domain_t & dom);

      Domain_t (const int i, const real8 *p1, const real8 *p2);

      Domain_t (             const real8 *p1, const real8 *p2);

      Domain_t (const int i, const real8 dx0, const real8 dx1,
                             const real8 dy0, const real8 dy1,
                             const real8 dz0, const real8 dz1 );

      Domain_t (             const real8 dx0, const real8 dx1,
                             const real8 dy0, const real8 dy1,
                             const real8 dz0, const real8 dz1 );

     ~Domain_t (void) {}

     const Domain_t & operator = (const Domain_t   & dom);
     const Domain_t & operator = (const RBDecomp_t & dom);

     int   operator == (const Domain_t & dom) const;

     void  Center        (real8 *c) const;
     void  Center        (real8 &cx, real8 &cy, real8 &cz) const;

     int   Intersect     (const Domain_t & dom) const;

     int   Inside        (const real8 *v) const;
     int   Inside        (const real8  x, const real8 y, const real8 z) const;

     int   Outside       (const real8 *v) const;
     int   Outside       (const real8  x, const real8 y, const real8 z) const;

     real8 Volume        (void) const;

     void  Expand        (const real8 lx, const real8 ly, const real8 lz,
                          const int   nx, const int   ny, const int   nz, const int na=0 );

     void  Native_Cells  (int & ix0, int & ix1,
                          int & iy0, int & iy1,
                          int & iz0, int & iz1,
                          const real8 lx , const real8 ly , const real8 lz ,
                          const int   nx , const int   ny , const int   nz  );

     int  *Native_Cells  (      int  & ccnt,
                          const real8 lx , const real8 ly , const real8 lz ,
                          const int   nx , const int   ny , const int   nz  );

     int  *Remote_Cells  (      int  & ccnt,
                          const real8 lx , const real8 ly , const real8 lz ,
                          const int   nx , const int   ny , const int   nz ,
                          const int   pbx, const int   pby, const int   pbz,
                          const int   rx , const int   ry , const int   rz  );

     void  Print       (FILE *fd=stdout) const;
     void  Print_GNU   (FILE *fd=stdout) const;

};

extern Domain_t *Append (Domain_t *doms, int & dcnt, int & dmax, Domain_t & dom);

extern Domain_t *Create_Domains(const real8 lx, const real8 ly, const real8 lz,
                                const int   nx, const int   ny, const int   nz );

extern Domain_t *Create_Domains(const RSDecomp_t *rs, const int nx, const int ny, const int nz );
extern Domain_t *Create_Domains(const RBDecomp_t *rb, const int nx, const int ny, const int nz );

#endif  // _PDS_DOMAIN_H
