#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Complex_t.h"
#include "AnisotropicVars_t.h"

// align8() - returns a 64-bit aligned byte count
//------------------------------------------------------------------------------------------------------------

inline size_t align8 (const size_t bytes) { return( bytes + ( (bytes%8) ? (8-(bytes%8)) : 0 ) ); } 

// AnisotropicVars_t (constructor)
//------------------------------------------------------------------------------------------------------------

AnisotropicVars_t::AnisotropicVars_t(void)
{
   qMax        = 0;

   nz_elem_re0 = 0;
   nz_elem_re  = 0;
   nz_elem_im0 = 0;
   nz_elem_im  = 0;

   memset( elasticConstantMatrix2D,0,sizeof(elasticConstantMatrix2D) );
   memset( elasticConstantMatrix4D,0,sizeof(elasticConstantMatrix4D) );

   FqReal     = 0;
   FqImag     = 0;

   FqReal_v3  = 0;
   FqImag_v3  = 0;

   coreFqReal = 0;
   coreFqImag = 0;

   sphericalHarmonicsRemKept  = 0;
   sphericalHarmonicsImmKept  = 0;

   sphericalHarmonicsRelKept  = 0;
   sphericalHarmonicsImlKept  = 0;

   sphericalHarmonicsRe18Kept = 0;
   sphericalHarmonicsIm18Kept = 0;

   sphericalHarmonicsReiKept  = 0;
   sphericalHarmonicsImiKept  = 0;

   memset(anisoRotMatrix12, 0, sizeof(anisoRotMatrix12));
   memset(anisoRotMatrix  , 0, sizeof(anisoRotMatrix  ));

   // set default rotation to identity...

   anisoRotMatrix[0][0] = 1.0;
   anisoRotMatrix[0][1] = 0.0;
   anisoRotMatrix[0][2] = 0.0;

   anisoRotMatrix[1][0] = 0.0;
   anisoRotMatrix[1][1] = 1.0;
   anisoRotMatrix[1][2] = 0.0;

   anisoRotMatrix[2][0] = 0.0;
   anisoRotMatrix[2][1] = 0.0;
   anisoRotMatrix[2][2] = 1.0;

   anisoRotMatrix12[0] = complex8(anisoRotMatrix[0][0],
                                  anisoRotMatrix[1][0]);

   anisoRotMatrix12[1] = complex8(anisoRotMatrix[0][1],
                                  anisoRotMatrix[1][1]);

   anisoRotMatrix12[2] = complex8(anisoRotMatrix[0][2],
                                  anisoRotMatrix[1][2]);

}

// Bytes()
//
// Will return the space required to store an instance of the aniso vars structure.  This is needed
// for reserving and copying the structure to the GPU.
//------------------------------------------------------------------------------------------------------------

size_t AnisotropicVars_t::Bytes (const int qmax)
{
   size_t bytes=0;

   bytes += align8( sizeof(AnisotropicVars_t) );                  // reserve space for container

   bytes += align8(     18*(qmax+1)*(qmax+2)*sizeof(real8) );     // (FqReal)
   bytes += align8(     18*(qmax+1)*(qmax+2)*sizeof(real8) );     // (FqImag)

   bytes += align8(     18*(qmax+1)*(qmax+2)*sizeof(real8) );     // (FqReal_v3)
   bytes += align8(     18*(qmax+1)*(qmax+2)*sizeof(real8) );     // (FqImag_v3)

   bytes += align8(  6*((qmax+1)*(qmax+2)-1)*sizeof(real8) );     // (coreFqReal)
   bytes += align8(  6*((qmax+1)*(qmax+2)-1)*sizeof(real8) );     // (coreFqImag)

   bytes += align8(                 (qmax+2)*sizeof(short) );     // (sphericalHarmonicsRemKept)
   bytes += align8(                 (qmax+2)*sizeof(short) );     // (sphericalHarmonicsImmKept)

   bytes += align8(        (qmax+1)*(qmax+2)*sizeof(short) );     // (sphericalHarmonicsRelKept)
   bytes += align8(        (qmax+1)*(qmax+2)*sizeof(short) );     // (sphericalHarmonicsImlKept)

   bytes += align8(        (qmax+1)*(qmax+2)*sizeof(short) );     // (sphericalHarmonicsRe18Kept)
   bytes += align8(        (qmax+1)*(qmax+2)*sizeof(short) );     // (sphericalHarmonicsIm18Kept)

   bytes += align8(     18*(qmax+1)*(qmax+2)*sizeof(short) );     // (sphericalHarmonicsReiKept)
   bytes += align8(     18*(qmax+1)*(qmax+2)*sizeof(short) );     // (sphericalHarmonicsImiKept)

   return(bytes);
}

// Bytes_Resized()
//
// Will return the space required to store an instance of the aniso vars structure AFTER
// it has been initialized and the tables resized. 
//------------------------------------------------------------------------------------------------------------

size_t AnisotropicVars_t::Bytes_Resized (void)
{
   size_t bytes=0;

   bytes += align8( sizeof(AnisotropicVars_t) );                  // reserve space for container

   bytes += align8(             nz_elem_re  *sizeof(real8) );     // (FqReal)
   bytes += align8(             nz_elem_im  *sizeof(real8) );     // (FqImag)

   bytes += align8(     18*(qMax+1)*(qMax+2)*sizeof(real8) );     // (FqReal_v3)
   bytes += align8(     18*(qMax+1)*(qMax+2)*sizeof(real8) );     // (FqImag_v3)

   bytes += align8(  6*((qMax+1)*(qMax+2)-1)*sizeof(real8) );     // (coreFqReal)
   bytes += align8(  6*((qMax+1)*(qMax+2)-1)*sizeof(real8) );     // (coreFqImag)

   bytes += align8(                 (qMax+2)*sizeof(short) );     // (sphericalHarmonicsRemKept)
   bytes += align8(                 (qMax+2)*sizeof(short) );     // (sphericalHarmonicsImmKept)

   bytes += align8(             nz_elem_re0 *sizeof(short) );     // (sphericalHarmonicsRelKept)
   bytes += align8(             nz_elem_im0 *sizeof(short) );     // (sphericalHarmonicsImlKept)

   bytes += align8(             nz_elem_re0 *sizeof(short) );     // (sphericalHarmonicsRe18Kept)
   bytes += align8(             nz_elem_im0 *sizeof(short) );     // (sphericalHarmonicsIm18Kept)

   bytes += align8(             nz_elem_re  *sizeof(short) );     // (sphericalHarmonicsReiKept)
   bytes += align8(             nz_elem_im  *sizeof(short) );     // (sphericalHarmonicsImiKept)

   return(bytes);
}

// Allocate()
//
// Will allocate and zero the space needed based on qMax.
//------------------------------------------------------------------------------------------------------------

void AnisotropicVars_t::Allocate (const int qMax)
{
   // allocate...

   FqReal                     = (real8 *) malloc(     18*(qMax+1)*(qMax+2)*sizeof(real8) );
   FqImag                     = (real8 *) malloc(     18*(qMax+1)*(qMax+2)*sizeof(real8) );

   FqReal_v3                  = (real8 *) malloc(     18*(qMax+1)*(qMax+2)*sizeof(real8) );
   FqImag_v3                  = (real8 *) malloc(     18*(qMax+1)*(qMax+2)*sizeof(real8) );

   coreFqReal                 = (real8 *) malloc(  6*((qMax+1)*(qMax+2)-1)*sizeof(real8) );
   coreFqImag                 = (real8 *) malloc(  6*((qMax+1)*(qMax+2)-1)*sizeof(real8) );

   sphericalHarmonicsRemKept  = (short *) malloc(                 (qMax+2)*sizeof(short) );
   sphericalHarmonicsImmKept  = (short *) malloc(                 (qMax+2)*sizeof(short) );

   sphericalHarmonicsRelKept  = (short *) malloc(        (qMax+1)*(qMax+2)*sizeof(short) );
   sphericalHarmonicsImlKept  = (short *) malloc(        (qMax+1)*(qMax+2)*sizeof(short) );

   sphericalHarmonicsRe18Kept = (short *) malloc(        (qMax+1)*(qMax+2)*sizeof(short) );
   sphericalHarmonicsIm18Kept = (short *) malloc(        (qMax+1)*(qMax+2)*sizeof(short) );

   sphericalHarmonicsReiKept  = (short *) malloc(     18*(qMax+1)*(qMax+2)*sizeof(short) );
   sphericalHarmonicsImiKept  = (short *) malloc(     18*(qMax+1)*(qMax+2)*sizeof(short) );

   // initialize all the allocated buffers to zero...

   if (FqReal                    ) { memset( FqReal                     , 0, (     18*(qMax+1)*(qMax+2)*sizeof(real8) ) ); }
   if (FqImag                    ) { memset( FqImag                     , 0, (     18*(qMax+1)*(qMax+2)*sizeof(real8) ) ); }

   if (FqReal_v3                 ) { memset( FqReal_v3                  , 0, (     18*(qMax+1)*(qMax+2)*sizeof(real8) ) ); }
   if (FqImag_v3                 ) { memset( FqImag_v3                  , 0, (     18*(qMax+1)*(qMax+2)*sizeof(real8) ) ); }

   if (coreFqReal                ) { memset( coreFqReal                 , 0, (  6*((qMax+1)*(qMax+2)-1)*sizeof(real8) ) ); }
   if (coreFqImag                ) { memset( coreFqImag                 , 0, (  6*((qMax+1)*(qMax+2)-1)*sizeof(real8) ) ); }
  
   if (sphericalHarmonicsRemKept ) { memset( sphericalHarmonicsRemKept  , 0, (                 (qMax+2)*sizeof(short) ) ); }
   if (sphericalHarmonicsImmKept ) { memset( sphericalHarmonicsImmKept  , 0, (                 (qMax+2)*sizeof(short) ) ); }

   if (sphericalHarmonicsRelKept ) { memset( sphericalHarmonicsRelKept  , 0, (        (qMax+1)*(qMax+2)*sizeof(short) ) ); }
   if (sphericalHarmonicsImlKept ) { memset( sphericalHarmonicsImlKept  , 0, (        (qMax+1)*(qMax+2)*sizeof(short) ) ); }

   if (sphericalHarmonicsRe18Kept) { memset( sphericalHarmonicsRe18Kept , 0, (        (qMax+1)*(qMax+2)*sizeof(short) ) ); }
   if (sphericalHarmonicsIm18Kept) { memset( sphericalHarmonicsIm18Kept , 0, (        (qMax+1)*(qMax+2)*sizeof(short) ) ); }

   if (sphericalHarmonicsReiKept ) { memset( sphericalHarmonicsReiKept  , 0, (     18*(qMax+1)*(qMax+2)*sizeof(short) ) ); }
   if (sphericalHarmonicsImiKept ) { memset( sphericalHarmonicsImiKept  , 0, (     18*(qMax+1)*(qMax+2)*sizeof(short) ) ); }
}

// Recycle()
//
// Will free any memory allocated via Allocate()...
//------------------------------------------------------------------------------------------------------------

void AnisotropicVars_t::Recycle (void)
{
   if (FqReal                    ) { free( FqReal                     ); FqReal=0; }
   if (FqImag                    ) { free( FqImag                     ); FqImag=0; }

   if (FqReal_v3                 ) { free( FqReal_v3                  ); FqReal_v3=0; }
   if (FqImag_v3                 ) { free( FqImag_v3                  ); FqImag_v3=0; }

   if (coreFqReal                ) { free( coreFqReal                 ); coreFqReal=0; }
   if (coreFqImag                ) { free( coreFqImag                 ); coreFqImag=0; }

   if (sphericalHarmonicsRemKept ) { free( sphericalHarmonicsRemKept  ); sphericalHarmonicsRemKept =0; }
   if (sphericalHarmonicsImmKept ) { free( sphericalHarmonicsImmKept  ); sphericalHarmonicsImmKept =0; }

   if (sphericalHarmonicsRelKept ) { free( sphericalHarmonicsRelKept  ); sphericalHarmonicsRelKept =0; }
   if (sphericalHarmonicsImlKept ) { free( sphericalHarmonicsImlKept  ); sphericalHarmonicsImlKept =0; }

   if (sphericalHarmonicsRe18Kept) { free( sphericalHarmonicsRe18Kept ); sphericalHarmonicsRe18Kept=0; }
   if (sphericalHarmonicsIm18Kept) { free( sphericalHarmonicsIm18Kept ); sphericalHarmonicsIm18Kept=0; }

   if (sphericalHarmonicsReiKept ) { free( sphericalHarmonicsReiKept  ); sphericalHarmonicsReiKept =0; }
   if (sphericalHarmonicsImiKept ) { free( sphericalHarmonicsImiKept  ); sphericalHarmonicsImiKept =0; }
}

// Set_EC_Matrices
// 
// Will initialize the full 4D (3x3x3x3) elastic constant matrix from the source 6x6 matrix.
//------------------------------------------------------------------------------------------------------------

void AnisotropicVars_t::Set_EC_Matrices (const real8 ecMatrix[6][6])
{
   memcpy(elasticConstantMatrix2D, ecMatrix, sizeof(elasticConstantMatrix2D));

   int i3to6[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };

   for (int i=0; (i<3); i++)
   for (int j=0; (j<3); j++)
   for (int k=0; (k<3); k++)
   for (int l=0; (l<3); l++)
      elasticConstantMatrix4D[i][j][k][l] = ecMatrix[i3to6[i][j]][i3to6[k][l]];
}

// Serialize()
//
// Given a pointer to a location in memory on the host (cpu) and device (gpu) - serializes a copy
// of the data structure that can be used on the GPU. Note that the internal pointers of the serialized
// result are relative to the buffer location on the GPU. The resulting structure cannot be used
// on the host.
//
// Returns the total space required (bytes).
//------------------------------------------------------------------------------------------------------------

size_t AnisotropicVars_t::Serialize 
(
   unsigned char *pc,   ///< points to buffer on host   (cpu)
   unsigned char *pg    ///< points to buffer on device (gpu)
) const
{
   unsigned char *p0 = pc;

   if (pc && pg)
   {
      AnisotropicVars_t *av = (AnisotropicVars_t *) pc;

      // serialize the individual components into the host buffer...

      size_t bytes=0;

      bytes = align8( sizeof(AnisotropicVars_t) );               memcpy(pc, this  , sizeof(AnisotropicVars_t) );                                                 pc+=bytes; pg+=bytes;

      bytes = align8(            nz_elem_re  *sizeof(real8) );   memcpy(pc, FqReal    , bytes);                   av->FqReal                     = (real8 *) pg; pc+=bytes; pg+=bytes;
      bytes = align8(            nz_elem_im  *sizeof(real8) );   memcpy(pc, FqImag    , bytes);                   av->FqImag                     = (real8 *) pg; pc+=bytes; pg+=bytes;

      bytes = align8(    18*(qMax+1)*(qMax+2)*sizeof(real8) );   memcpy(pc, FqReal_v3 , bytes);                   av->FqReal_v3                  = (real8 *) pg; pc+=bytes; pg+=bytes;
      bytes = align8(    18*(qMax+1)*(qMax+2)*sizeof(real8) );   memcpy(pc, FqImag_v3 , bytes);                   av->FqImag_v3                  = (real8 *) pg; pc+=bytes; pg+=bytes;

      bytes = align8( 6*((qMax+1)*(qMax+2)-1)*sizeof(real8) );   memcpy(pc, coreFqReal, bytes);                   av->coreFqReal                 = (real8 *) pg; pc+=bytes; pg+=bytes;
      bytes = align8( 6*((qMax+1)*(qMax+2)-1)*sizeof(real8) );   memcpy(pc, coreFqImag, bytes);                   av->coreFqImag                 = (real8 *) pg; pc+=bytes; pg+=bytes;

      bytes = align8(                (qMax+2)*sizeof(short) );   memcpy(pc, sphericalHarmonicsRemKept , bytes);   av->sphericalHarmonicsRemKept  = (short *) pg; pc+=bytes; pg+=bytes;
      bytes = align8(                (qMax+2)*sizeof(short) );   memcpy(pc, sphericalHarmonicsImmKept , bytes);   av->sphericalHarmonicsImmKept  = (short *) pg; pc+=bytes; pg+=bytes;

      bytes = align8(            nz_elem_re0 *sizeof(short) );   memcpy(pc, sphericalHarmonicsRelKept , bytes);   av->sphericalHarmonicsRelKept  = (short *) pg; pc+=bytes; pg+=bytes;
      bytes = align8(            nz_elem_im0 *sizeof(short) );   memcpy(pc, sphericalHarmonicsImlKept , bytes);   av->sphericalHarmonicsImlKept  = (short *) pg; pc+=bytes; pg+=bytes;

      bytes = align8(            nz_elem_re0 *sizeof(short) );   memcpy(pc, sphericalHarmonicsRe18Kept, bytes);   av->sphericalHarmonicsRe18Kept = (short *) pg; pc+=bytes; pg+=bytes;
      bytes = align8(            nz_elem_im0 *sizeof(short) );   memcpy(pc, sphericalHarmonicsIm18Kept, bytes);   av->sphericalHarmonicsIm18Kept = (short *) pg; pc+=bytes; pg+=bytes;

      bytes = align8(            nz_elem_re  *sizeof(short) );   memcpy(pc, sphericalHarmonicsReiKept , bytes);   av->sphericalHarmonicsReiKept  = (short *) pg; pc+=bytes; pg+=bytes;
      bytes = align8(            nz_elem_im  *sizeof(short) );   memcpy(pc, sphericalHarmonicsImiKept , bytes);   av->sphericalHarmonicsImiKept  = (short *) pg; pc+=bytes; pg+=bytes;
   }

   return( (size_t) (pc-p0) );
}

// Validate()
//
// All of the indices within the aniso table should be within the run-time initialized qMax parameter.
// This routine will check to insure that all the parameters are within the bounds of the allocated tables.
//------------------------------------------------------------------------------------------------------------

int Valid ( const short *v, const int vn, const int vmin, const int vmax )
{
   if (v && (vn>0))
   {
      for (int i=0; (i<vn); i++) 
         if ( (v[i]<vmin) || (v[i]>=vmax) ) { return(0); }
   }

   return(1);
}

int AnisotropicVars_t::Validate(void)
{
   int vmin =  0;
   int vmax = 18*(qMax+1)*(qMax+2);

   if ( !Valid( sphericalHarmonicsRemKept  , qMax+2     , vmin, vmax ) ) return(0);
   if ( !Valid( sphericalHarmonicsImmKept  , qMax+2     , vmin, vmax ) ) return(0);

   if ( !Valid( sphericalHarmonicsRelKept  , nz_elem_re0, vmin, vmax ) ) return(0);
   if ( !Valid( sphericalHarmonicsImlKept  , nz_elem_im0, vmin, vmax ) ) return(0);

   if ( !Valid( sphericalHarmonicsRe18Kept , nz_elem_re0, vmin, vmax ) ) return(0);
   if ( !Valid( sphericalHarmonicsIm18Kept , nz_elem_im0, vmin, vmax ) ) return(0);

   if ( !Valid( sphericalHarmonicsReiKept  , nz_elem_re , vmin, vmax ) ) return(0);
   if ( !Valid( sphericalHarmonicsImiKept  , nz_elem_im , vmin, vmax ) ) return(0);

   return(1);
} 


static void Print_Vector(FILE *fd, const real8 *v, const int vn, const char *fmt, const char *indent, const int lcnt)
{
   if (fd && v && (vn>0))
   {
      fprintf(fd,"\n%s", indent);

      for (int i=0; (i<vn); )
      {
         int n = ( ((i+lcnt)<vn) ? lcnt : (vn-i) );

         for (int j=0; (j<n); j++) { fprintf(fd,fmt, v[i++]); }

         fprintf(fd,"\n%s", indent); 
      }
   }
}
      
static void Print_Vector(FILE *fd, const short *v, const int vn, const char *fmt, const char *indent, const int lcnt)
{
   if (fd && v && (vn>0))
   {
      fprintf(fd,"\n%s", indent);

      for (int i=0; (i<vn); )
      {
         int n = ( ((i+lcnt)<vn) ? lcnt : (vn-i) );

         for (int j=0; (j<n); j++) { fprintf(fd,fmt, v[i++]); }

         fprintf(fd,"\n%s", indent); 
      }
   }
}
      
// Print()
//
// Will print the aniso tables in a human-readable form...
//------------------------------------------------------------------------------------------------------------

void AnisotropicVars_t::Print (FILE *fd) const
{
   if (fd)
   {
      fprintf(fd,"AnisotropicVars_t::\n");
      fprintf(fd,"   qmax        = %d\n", qMax );
      fprintf(fd,"   nz_elem_re0 = %d\n", nz_elem_re0 );
      fprintf(fd,"   nz_elem_re  = %d\n", nz_elem_re  );
      fprintf(fd,"   nz_elem_im0 = %d\n", nz_elem_im0 );
      fprintf(fd,"   nz_elem_im  = %d\n", nz_elem_im  );
      fprintf(fd,"\n");
      fprintf(fd,"   FqReal                      = 0x%08lx\n", (unsigned long) FqReal                     );
      fprintf(fd,"   FqImag                      = 0x%08lx\n", (unsigned long) FqImag                     );
      fprintf(fd,"\n");
      fprintf(fd,"   FqReal_v3                   = 0x%08lx\n", (unsigned long) FqReal_v3                  );
      fprintf(fd,"   FqImag_v3                   = 0x%08lx\n", (unsigned long) FqImag_v3                  );
      fprintf(fd,"\n");
      fprintf(fd,"   coreFqReal                  = 0x%08lx\n", (unsigned long) coreFqReal                 );
      fprintf(fd,"   coreFqImag                  = 0x%08lx\n", (unsigned long) coreFqImag                 );
      fprintf(fd,"\n");
      fprintf(fd,"   sphericalHarmonicsRemKept   = 0x%08lx\n", (unsigned long) sphericalHarmonicsRemKept  );
      fprintf(fd,"   sphericalHarmonicsImmKept   = 0x%08lx\n", (unsigned long) sphericalHarmonicsImmKept  );
      fprintf(fd,"\n");
      fprintf(fd,"   sphericalHarmonicsRelKept   = 0x%08lx\n", (unsigned long) sphericalHarmonicsRelKept  );
      fprintf(fd,"   sphericalHarmonicsImlKept   = 0x%08lx\n", (unsigned long) sphericalHarmonicsImlKept  );
      fprintf(fd,"\n");
      fprintf(fd,"   sphericalHarmonicsRe18Kept  = 0x%08lx\n", (unsigned long) sphericalHarmonicsRe18Kept );
      fprintf(fd,"   sphericalHarmonicsIm18Kept  = 0x%08lx\n", (unsigned long) sphericalHarmonicsIm18Kept );
      fprintf(fd,"\n");
      fprintf(fd,"   sphericalHarmonicsReiKept   = 0x%08lx\n", (unsigned long) sphericalHarmonicsReiKept  );
      fprintf(fd,"   sphericalHarmonicsImiKept   = 0x%08lx\n", (unsigned long) sphericalHarmonicsImiKept  );
      fprintf(fd,"\n");
      fprintf(fd,"   elasticConstantMatrix2D[6][6] = \n" );
 
      for (int i=0; (i<6); i++)
      {
         fprintf(fd,"      %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf\n",
            elasticConstantMatrix2D[i][0],
            elasticConstantMatrix2D[i][1],
            elasticConstantMatrix2D[i][2],
            elasticConstantMatrix2D[i][3],
            elasticConstantMatrix2D[i][4],
            elasticConstantMatrix2D[i][5] );
      }

      fprintf(fd,"\n");
      fprintf(fd,"   elasticConstantMatrix4D[3][3][3][3] = \n" );
 
      for (int i=0; (i<3); i++)
      for (int j=0; (j<3); j++)
      {
         for (int k=0; (k<3); k++)
         {
            fprintf(fd,"      %8.4lf %8.4lf %8.4lf\n",  
               elasticConstantMatrix4D[i][j][k][0],
               elasticConstantMatrix4D[i][j][k][1],
               elasticConstantMatrix4D[i][j][k][2] );
         }
         fprintf(fd,"\n");
      }

      fprintf(fd,"\n");
      fprintf(fd,"   anisoRotMatrix[3][3] = \n" );
 
      for (int i=0; (i<3); i++)
      {
         fprintf(fd,"      %10.6lf %10.6lf %10.6lf\n", anisoRotMatrix[i][0], anisoRotMatrix[i][1], anisoRotMatrix[i][2] );
      }

      fprintf(fd,"\n   FqReal     (n=%d) = ", nz_elem_re ); Print_Vector(fd, FqReal, nz_elem_re, "%+0.6le ", "      ", 10 );
      fprintf(fd,"\n   FqImag     (n=%d) = ", nz_elem_im ); Print_Vector(fd, FqImag, nz_elem_im, "%+0.6le ", "      ", 10 );

      int n = 18*(qMax+1)*(qMax+2);
      fprintf(fd,"\n   FqReal_v3  (n=%d) = ", n ); Print_Vector(fd, FqReal_v3  , n, "%+0.6le ", "      ", 10 );
      fprintf(fd,"\n   FqImag_v3  (n=%d) = ", n ); Print_Vector(fd, FqImag_v3  , n, "%+0.6le ", "      ", 10 );

          n = 6*((qMax+1)*(qMax+2)-1);
      fprintf(fd,"\n   coreFqReal (n=%d) = ", n ); Print_Vector(fd, coreFqReal , n, "%+0.6le ", "      ", 10 );
      fprintf(fd,"\n   coreFqImag (n=%d) = ", n ); Print_Vector(fd, coreFqImag , n, "%+0.6le ", "      ", 10 );

      fprintf(fd,"\n   sphericalHarmonicsRemKept  (n=%d) = ", qMax+2      ); Print_Vector(fd, sphericalHarmonicsRemKept , qMax+2     , "%3d ", "      ", qMax+2 );
      fprintf(fd,"\n   sphericalHarmonicsImmKept  (n=%d) = ", qMax+2      ); Print_Vector(fd, sphericalHarmonicsImmKept , qMax+2     , "%3d ", "      ", qMax+2 );

      fprintf(fd,"\n   sphericalHarmonicsRelKept  (n=%d) = ", nz_elem_re0 ); Print_Vector(fd, sphericalHarmonicsRelKept , nz_elem_re0, "%3d ", "      ", 40 );
      fprintf(fd,"\n   sphericalHarmonicsImlKept  (n=%d) = ", nz_elem_im0 ); Print_Vector(fd, sphericalHarmonicsImlKept , nz_elem_im0, "%3d ", "      ", 40 );

      fprintf(fd,"\n   sphericalHarmonicsRe18Kept (n=%d) = ", nz_elem_re0 ); Print_Vector(fd, sphericalHarmonicsRe18Kept, nz_elem_re0, "%3d ", "      ", 40 );
      fprintf(fd,"\n   sphericalHarmonicsIm18Kept (n=%d) = ", nz_elem_im0 ); Print_Vector(fd, sphericalHarmonicsIm18Kept, nz_elem_im0, "%3d ", "      ", 40 );

      fprintf(fd,"\n   sphericalHarmonicsReiKept  (n=%d) = ", nz_elem_re  ); Print_Vector(fd, sphericalHarmonicsReiKept , nz_elem_re , "%3d ", "      ", 40 );
      fprintf(fd,"\n   sphericalHarmonicsImiKept  (n=%d) = ", nz_elem_im  ); Print_Vector(fd, sphericalHarmonicsImiKept , nz_elem_im , "%3d ", "      ", 40 );
      
      fprintf(fd,"\n");
   }
}

