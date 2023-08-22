#pragma once

#ifndef _PDS_SSF_ISO_RATIONAL_H
#define _PDS_SSF_ISO_RATIONAL_H

#include "cuda_portability.h"
#include "Typedefs.h"

__cuda_hdev__
void SSF_Iso_Rational_Correction_Half
(
         real8    *f1   ,  ///< corrected force on node 1 <xyz> (updated)
         real8    *f2   ,  ///< corrected force on node 2 <xyz> (updated)
   const real8    *p1   ,  ///< position of node 1 <xyz>
   const real8    *p2   ,  ///< position of node 2 <xyz>
   const real8    *p3   ,  ///< position of node 3 <xyz>
   const real8    *p4   ,  ///< position of node 4 <xyz>
   const real8    *b1   ,  ///< burgers vector for segment p1->p2
   const real8    *b3   ,  ///< burgers vector for segment p3->p4
   const real8     a    ,  ///< core radius
   const real8     mu   ,  ///< shear modulus
   const real8     nu   ,  ///< poisson ratio
   const real8     ecrit   ///< critical angle for parallelism test
);

__cuda_hdev__
void SSF_Iso_Rational_Correction
(
         real8    *f1   ,  ///< corrected force on node 1 <xyz> (updated)
         real8    *f2   ,  ///< corrected force on node 2 <xyz> (updated)
         real8    *f3   ,  ///< corrected force on node 3 <xyz> (updated)
         real8    *f4   ,  ///< corrected force on node 4 <xyz> (updated)
   const real8    *p1   ,  ///< position of node 1 <xyz>
   const real8    *p2   ,  ///< position of node 2 <xyz>
   const real8    *p3   ,  ///< position of node 3 <xyz>
   const real8    *p4   ,  ///< position of node 4 <xyz>
   const real8    *b1   ,  ///< burgers vector for segment p1->p2
   const real8    *b3   ,  ///< burgers vector for segment p3->p4
   const real8     a    ,  ///< core radius
   const real8     mu   ,  ///< shear modulus
   const real8     nu   ,  ///< poisson ratio
   const real8     ecrit   ///< critical angle for parallelism test
);


#endif // _PDS_SSF_ISO_RATIONAL_H
