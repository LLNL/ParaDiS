#pragma once

#ifndef _PDS_ANISOTROPIC_VARS_H
#define _PDS_ANISOTROPIC_VARS_H

#include <stdio.h>
#include <stdlib.h>

#include "Typedefs.h"
#include "Complex_t.h"

// Define a structure containing various items that are specific 
// to the anisotropic elasticity but not provided as control
// parameters.

class AnisotropicVars_t
{
   public :
      int   qMax       ;  ///< qMax that was used to initialize the tables

      int   nz_elem_re0;  ///< number of non-zero real      elements retained (set in AnisotropicInit())
      int   nz_elem_re ;  ///< number of non-zero real      elements retained (set in AnisotropicInit())
      int   nz_elem_im0;  ///< number of non-zero imaginary elements retained (set in AnisotropicInit())
      int   nz_elem_im ;  ///< number of non-zero imaginary elements retained (set in AnisotropicInit())

      // Various portions of the code require the elastic constant matrix 
      // in either the 6x6 or 3X3X3X3 forms, so we make a local copy of the 
      // 6x6 matrix here and transform it to the 3x3x3x3 form here.

      real8 elasticConstantMatrix2D[6][6];
      real8 elasticConstantMatrix4D[3][3][3][3];

      // The anisotropic elasticity requires the <Fq> table.  This
      // table is calculated during initialization and is at first
      // a 3D table of complex doubles of dimensions [6][3][(qMax+2)*(qMax+1)].
      // However, many of the elements of the table may be zero, so the
      // table is flattened into a two single dimension arrays consisting
      // of only the non-zero real and imaginary components respectively.

      real8 *FqReal;
      real8 *FqImag;

      // The current GPU implementation of aniso forces uses an uncompressed
      // copy of the <Fq> table...

      real8 *FqReal_v3;
      real8 *FqImag_v3;

      // The core force calculation of a dislocation also requires an
      // <Fq> table, however, it is a different table than the table above.
      // The core force Fq table is also calculated during initialization
      // and consists of a pair of arrays containing the real and imaginary
      // portions of complex values of dimensions [(qMax+1)*(qMax+3)]

      real8 *coreFqReal;
      real8 *coreFqImag;

      // The anisotropic force routine is described in the paper entitled
      // "Use of spherical harmonics for disloction dynamics in anisotropic
      // elastic media" in the <ParadisDir>/docs directory.  The force
      // routine is an implementation begining with Equation 13 and ending
      // on page 10 of the paper.
      //
      // Within the force routine there are many zero entries in S(i,j,ii)
      // (which is defined in the above mentioned paper), and rather than do
      // many zero multiplies and have many zeros stored in the matrix the
      // array is stored in a column vector with the zeros removed.
      // Additional arrays are needed to specify which entries remain and to
      // make effective use of the analytical integrals used to compute forces.
      // In describingn these arrays it is more useful to use the expanded
      // notation of S rather than the collapsed into the column vector.
      //
      // Given the definition of S(i,j,ii), these additional arrays are
      // defined below.
      //
      // sphericalHarmonicsRemKept[i] is the number of j elements for
      // a given i that have at least one (j,ii) element that has
      // a non-zero 'real' portion.  sphericalHarmonicsImmKept[i] is
      // the corresponding array for the 'imaginary' portions of (j,ii).

      short   *sphericalHarmonicsRemKept;
      short   *sphericalHarmonicsImmKept;

      // sphericalHarmonicsRelKept contains the j indices associated
      // with the non-zero entries in sphericalHarmonicsRemKept[].
      // The number of elements in sphericalHarmonicsRelKept equals
      // the sum of all entries in sphericalHarmonicsRemKept. 
      // sphericalHarmonicsImlKept is the corresponding array for
      // sphericalHarmonicsImmKept.

      short   *sphericalHarmonicsRelKept;
      short   *sphericalHarmonicsImlKept;

      // sphericalHarmonicsRe18Kept[] contains the number of ii elements
      // associated with the an (i,j) that are non zero.  
      // sphericalHarmonicsIm18Kept[] is the corresponding array for
      // the 'imaginary' portions.

      short   *sphericalHarmonicsRe18Kept;
      short   *sphericalHarmonicsIm18Kept;

      // sphericalHarmonicsReiKept[] contains the indices (between 0
      // and 17) associated with the non-zero entries in
      // sphericalHarmonicsRe18Kept.  The number of elements in
      // sphericalHarmonicsRe18Kept equals the sum of all entries
      // in sphericalHarmonicsRem18Kept.  sphericalHarmonicsImiiKept
      // is the corresponding array for sphericalHarmonicsIm18Kept.

      short   *sphericalHarmonicsReiKept;
      short   *sphericalHarmonicsImiKept;

      // The following value is calculated from the elastic constant tensor
      // rotation matrix <anisoRotMatrix> and is defined as:
      //
      //    anisoRotMatrix[0][*] + I*anisoRotMatrix[1][*]

      complex8 anisoRotMatrix12[3];

      // Define the basis in which the elastic tensor C_ijkl is defined
      // for the spherical harmonics.

      real8 anisoRotMatrix[3][3];

   public :
      AnisotropicVars_t(void);
     ~AnisotropicVars_t() {} // (nothing allocated)

      static size_t Bytes           (const int qmax);
             size_t Bytes_Resized   (void);

             void   Allocate        (const int qmax);
             void   Recycle         (void);

             void   Set_EC_Matrices (const real8 ecMatrix[6][6]);

             size_t Serialize       (unsigned char *pc, unsigned char *pg) const;

             int    Validate        (void);

             void   Print           (FILE *fd=stdout) const;
};

#endif   // _PDS_ANISOTROPIC_VARS_H
