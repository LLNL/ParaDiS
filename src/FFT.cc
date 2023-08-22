//-----------------------------------------------------------------------------------------------
// FFT - implementation of a simple fast Fourier transform.
//-----------------------------------------------------------------------------------------------

#include <math.h>

#include "FFT.h"

static int bfly_4[4]     = {  0,  2,  1,  3 };

static int bfly_8[8]     = {  0,  4,  2,  6,  1,  5,  3,  7 };

static int bfly_16[16]   = {  0,  8,  4, 12,  2, 10,  6, 14,  1,  9,  5, 13,  3, 11,  7, 15  };

static int bfly_32[32]   = {  0, 16,  8, 24,  4, 20, 12, 28,  2, 18, 10, 26,  6, 22, 14, 30,
                              1, 17,  9, 25,  5, 21, 13, 29,  3, 19, 11, 27,  7, 23, 15, 31  };

static int bfly_64[64]   = {  0, 32, 16, 48,  8, 40, 24, 56,  4, 36, 20, 52, 12, 44, 28, 60,
                              2, 34, 18, 50, 10, 42, 26, 58,  6, 38, 22, 54, 14, 46, 30, 62,
                              1, 33, 17, 49,  9, 41, 25, 57,  5, 37, 21, 53, 13, 45, 29, 61,
                              3, 35, 19, 51, 11, 43, 27, 59,  7, 39, 23, 55, 15, 47, 31, 63  };

static int bfly_128[128] = {  0, 64, 32, 96, 16, 80, 48,112,  8, 72, 40,104, 24, 88, 56,120,
                              4, 68, 36,100, 20, 84, 52,116, 12, 76, 44,108, 28, 92, 60,124,
                              2, 66, 34, 98, 18, 82, 50,114, 10, 74, 42,106, 26, 90, 58,122,
                              6, 70, 38,102, 22, 86, 54,118, 14, 78, 46,110, 30, 94, 62,126,
                              1, 65, 33, 97, 17, 81, 49,113,  9, 73, 41,105, 25, 89, 57,121,
                              5, 69, 37,101, 21, 85, 53,117, 13, 77, 45,109, 29, 93, 61,125,
                              3, 67, 35, 99, 19, 83, 51,115, 11, 75, 43,107, 27, 91, 59,123,
                              7, 71, 39,103, 23, 87, 55,119, 15, 79, 47,111, 31, 95, 63,127  };

static int bfly_256[256] = {  0,128, 64,192, 32,160, 96,224, 16,144, 80,208, 48,176,112,240,
                              8,136, 72,200, 40,168,104,232, 24,152, 88,216, 56,184,120,248,
                              4,132, 68,196, 36,164,100,228, 20,148, 84,212, 52,180,116,244,
                             12,140, 76,204, 44,172,108,236, 28,156, 92,220, 60,188,124,252,
                              2,130, 66,194, 34,162, 98,226, 18,146, 82,210, 50,178,114,242,
                             10,138, 74,202, 42,170,106,234, 26,154, 90,218, 58,186,122,250,
                              6,134, 70,198, 38,166,102,230, 22,150, 86,214, 54,182,118,246,
                             14,142, 78,206, 46,174,110,238, 30,158, 94,222, 62,190,126,254,
                              1,129, 65,193, 33,161, 97,225, 17,145, 81,209, 49,177,113,241,
                              9,137, 73,201, 41,169,105,233, 25,153, 89,217, 57,185,121,249,
                              5,133, 69,197, 37,165,101,229, 21,149, 85,213, 53,181,117,245,
                             13,141, 77,205, 45,173,109,237, 29,157, 93,221, 61,189,125,253,
                              3,131, 67,195, 35,163, 99,227, 19,147, 83,211, 51,179,115,243,
                             11,139, 75,203, 43,171,107,235, 27,155, 91,219, 59,187,123,251,
                              7,135, 71,199, 39,167,103,231, 23,151, 87,215, 55,183,119,247,
                             15,143, 79,207, 47,175,111,239, 31,159, 95,223, 63,191,127,255  };

// FFT_Butterfly()
//
// Will allocate and return a vector containing the radix-2 butterfly of length N
// Not a particularly fast implementation - should only be used once to generate 
// a vector.
//---------------------------------------------------------------------------------------------------------

int *FFT_Butterfly (const int n)
{
   int  m=1;                                       
   for (m=1; ((1<<m)<n); m++);                     
                                                   
   int *bfly = ( (n>0) ? new int[n] : 0 );         
                                                   
   if (bfly)                                       
   {                                               
      for (int i=0; (i<n); i++) { bfly[i]=i; }     
                                                   
      // --------  Bit Reversal  ----------        
                                                   
      int j  = 0;                                  
      int n2 = (n/2);                              
                                                   
      for (int i=1; i<(n-1); i++)                  
      {                                            
         int n1 = n2;                              
                                                   
         while(j>=n1)                              
         {                                         
            j  = j-n1;                             
            n1 = (n1/2);                           
         }                                         
                                                   
         j += n1;                                  
                                                   
         if (i<j)                                  
         {                                         
            int t   = bfly[i];                     
            bfly[i] = bfly[j];                     
            bfly[j] = t;                           
         }                                         
      }                                            
   }                                               
                                                   
   return (bfly);                                  
}

//---------------------------------------------------------------------------------------------------------

#define FFT_LOAD(q,p,b,n)                   \
if ( (q) && (p) && (b) && ((n)>=4) )        \
{                                           \
   for (int i=0,k=0; (i<(n)); ++i)          \
   {                                        \
      (q)[k++]=(p)[(b)[i]];                 \
      (q)[k++]=0.0;                         \
   }                                        \
}

#define FFT_LOAD_COMPLEX(q,p,b,n)           \
if ( (q) && (p) && (b) && ((n)>=4) )        \
{                                           \
   for (int i=0,k=0; (i<(n)); ++i)          \
   {                                        \
      int j = 2*(b)[i];                     \
      (q)[k++]=(p)[j  ];                    \
      (q)[k++]=(p)[j+1];                    \
   }                                        \
}

#define FFT_LOAD_COMPLEX_STRIDE(q,p,b,n,s)  \
if ( (q) && (p) && (b) && ((n)>=4) )        \
{                                           \
   for (int i=0,k=0; (i<(n)); ++i)          \
   {                                        \
      int j = s*(b)[i];                     \
      (q)[k++]=(p)[j  ];                    \
      (q)[k++]=(p)[j+1];                    \
   }                                        \
}

void FFT_Load (float  *dst, const unsigned char   *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }
void FFT_Load (float  *dst, const unsigned short  *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }
void FFT_Load (float  *dst, const          short  *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }
void FFT_Load (float  *dst, const unsigned int    *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }
void FFT_Load (float  *dst, const          int    *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }
void FFT_Load (float  *dst, const          float  *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }

void FFT_Load (double *dst, const unsigned char   *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }
void FFT_Load (double *dst, const unsigned short  *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }
void FFT_Load (double *dst, const          short  *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }
void FFT_Load (double *dst, const unsigned int    *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }
void FFT_Load (double *dst, const          int    *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }
void FFT_Load (double *dst, const          double *src, const int *bfly, const int n) { FFT_LOAD(dst,src,bfly,n) }

void FFT_Load_Complex (float  *dst, const float  *src, const int *bfly, const int n) { FFT_LOAD_COMPLEX(dst,src,bfly,n) }
void FFT_Load_Complex (double *dst, const double *src, const int *bfly, const int n) { FFT_LOAD_COMPLEX(dst,src,bfly,n) }

void FFT_Load_Complex (float  *dst, const float  *src, const int *bfly, const int n, const int m) { FFT_LOAD_COMPLEX_STRIDE(dst,src,bfly,n,m) }
void FFT_Load_Complex (double *dst, const double *src, const int *bfly, const int n, const int m) { FFT_LOAD_COMPLEX_STRIDE(dst,src,bfly,n,m) }

// FFT_Power_Spectrum()
//
// Converts a complex vector to power spectrum - q[i].real = sqrt(p[i].real*p[i].real + p[i].imag*p[i].imag)
// Note - result is packed vector containing only real's.
//---------------------------------------------------------------------------------------------------------

void  FFT_Power_Spectrum (float  *v, const int vn) { FFT_Power_Spectrum(v,v,vn); }
void  FFT_Power_Spectrum (double *v, const int vn) { FFT_Power_Spectrum(v,v,vn); }

void  FFT_Power_Spectrum (float  *dst, const float  *src, const int n)
{
   if (dst && src && (n>0))
   {
      for (int i=0,k=0; (i<n); ++i, k+=2)
      {
         float a=src[k  ];
         float b=src[k+1]; 
         dst[i]=sqrtf(a*a+b*b);
      }
   }
}

void  FFT_Power_Spectrum (double *dst, const double *src, const int n)
{
   if (dst && src && (n>0))
   {
      for (int i=0,k=0; (i<n); ++i, k+=2)
      {
         double a=src[k  ];
         double b=src[k+1]; 
         dst[i]=sqrt(a*a+b*b);
      }
   }
}

void FFT_Normalize (float  *v, const int n)
{
   float s = ( (n>0) ? (1.0/n) : 0.0 );

   for (int i=0; (i<(2*n)); ++i)  { v[i]*=s; }
}

void FFT_Normalize (double *v, const int n)
{
   double s = ( (n>0) ? (1.0/n) : 0.0 );

   for (int i=0; (i<(2*n)); ++i)  { v[i]*=s; }
}

// FFT_Transform()
//
// This is a revised implementation of Cooley Tukey FFT algorithm as extracted 
// from Numerical Recipes in C.  The revision basically uses the input complex 
// vector as a zero-indexed array rather than the original one-indexed implementation. 
// The bit reversal butterfly calculation has been pulled out of the transform to 
// improve the performance of repeated calls to the transform.
//---------------------------------------------------------------------------------------------------------

void FFT_Transform (float *v, const int vn, const int dir, const int norm)
{
   int    n     = (2*vn);
   int    istep = 0;

   float  isign = ( (dir<0) ? -1.0f : 1.0f );

   for (int mmax=2; (n>mmax); mmax=istep) 
   {
            istep = 2*mmax;

      float theta = isign*((2.0f*M_PI)/mmax);
      float wtemp = sinf(theta/2.0f);
      float wpr   = -2.0f*wtemp*wtemp;
      float wpi   = sinf(theta);
      float wr    = 1.0f;
      float wi    = 0.0f;

      for (int m=1; (m<mmax); m+=2) 
      {
         for (int i=m; (i<=n); i+=istep) 
         {
            int   j       = i+mmax;
            float tempr   = wr*v[j-1]-wi*v[j  ];
            float tempi   = wr*v[j  ]+wi*v[j-1];
                  v[j-1]  = v[i-1]-tempr;
                  v[j  ]  = v[i  ]-tempi;
                  v[i-1] += tempr;
                  v[i  ] += tempi;
         }

         wr = (wtemp=wr)*wpr-wi*wpi+wr;
         wi = wi*wpr+wtemp*wpi+wi;
      }

      mmax=istep;
   }

   if (norm>0)
      FFT_Normalize(v,vn);
}

void FFT_Transform (double *v, const int vn, const int dir, const int norm)
{
   int    n     = (2*vn);
   int    istep = 0;

   double isign = ( (dir<0) ? -1.0 : 1.0 );

   for (int mmax=2; (n>mmax); mmax=istep) 
   {
             istep = 2*mmax;

      double theta = isign*((2.0*M_PI)/mmax);
      double wtemp = sin(theta/2.0);
      double wpr   = -2.0*wtemp*wtemp;
      double wpi   = sin(theta);
      double wr    = 1.0;
      double wi    = 0.0;

      for (int m=1; (m<mmax); m+=2) 
      {
         for (int i=m; (i<=n); i+=istep) 
         {
            int    j       = i+mmax;
            double tempr   = wr*v[j-1]-wi*v[j  ];
            double tempi   = wr*v[j  ]+wi*v[j-1];
                   v[j-1]  = v[i-1]-tempr;
                   v[j  ]  = v[i  ]-tempi;
                   v[i-1] += tempr;
                   v[i  ] += tempi;
         }

         wr = (wtemp=wr)*wpr-wi*wpi+wr;
         wi = wi*wpr+wtemp*wpi+wi;
      }

      mmax=istep;
   }

   if (norm>0)
      FFT_Normalize(v,vn);
}

//---------------------------------------------------------------------------------------------------------

float *FFT (float *v, const int n, const int dir, const int norm)
{
   float *fft = ( (v && (n>0)) ? new float[2*n] : 0 );

   if (v && fft && (n>0))
   {
      int *bfly=0, bfly_alloc=0;

      switch(n)
      {
         case (  4) : { bfly=bfly_4;   break; }
         case (  8) : { bfly=bfly_8;   break; }
         case ( 16) : { bfly=bfly_16;  break; }
         case ( 32) : { bfly=bfly_32;  break; }
         case ( 64) : { bfly=bfly_64;  break; }
         case (128) : { bfly=bfly_128; break; }
         case (256) : { bfly=bfly_256; break; }
      }

      if (!bfly) { bfly = FFT_Butterfly(n); bfly_alloc=1; }

      if (bfly)
      {
         FFT_Load      (fft,v,bfly,n);
         FFT_Transform (fft,n,dir,norm);
      }

      if (bfly && bfly_alloc) { delete [] bfly; }
   }

   return(fft);
}

double *FFT (double *v, const int n, const int dir, const int norm)
{
   double *fft = ( (v && (n>0)) ? new double[2*n] : 0 );

   if (v && fft && (n>0))
   {
      int *bfly=0, bfly_alloc=0;

      switch(n)
      {
         case (  4) : { bfly=bfly_4;   break; }
         case (  8) : { bfly=bfly_8;   break; }
         case ( 16) : { bfly=bfly_16;  break; }
         case ( 32) : { bfly=bfly_32;  break; }
         case ( 64) : { bfly=bfly_64;  break; }
         case (128) : { bfly=bfly_128; break; }
         case (256) : { bfly=bfly_256; break; }
      }

      if (!bfly) { bfly = FFT_Butterfly(n); bfly_alloc=1; }

      if (bfly)
      {
         FFT_Load      (fft,v,bfly,n);
         FFT_Transform (fft,n,dir,norm);
      }

      if (bfly && bfly_alloc) { delete [] bfly; }
   }

   return(fft);
}

float *FFT (float *fft, float *v, const int n, const int dir, const int norm)
{
   if (!fft) fft = ( (v && (n>0)) ? new float[2*n] : 0 );

   if (v && fft && (n>0))
   {
      int *bfly=0, bfly_alloc=0;

      switch(n)
      {
         case (  4) : { bfly=bfly_4;   break; }
         case (  8) : { bfly=bfly_8;   break; }
         case ( 16) : { bfly=bfly_16;  break; }
         case ( 32) : { bfly=bfly_32;  break; }
         case ( 64) : { bfly=bfly_64;  break; }
         case (128) : { bfly=bfly_128; break; }
         case (256) : { bfly=bfly_256; break; }
      }

      if (!bfly) { bfly = FFT_Butterfly(n); bfly_alloc=1; }

      if (bfly)
      {
         FFT_Load      (fft,v,bfly,n);
         FFT_Transform (fft,n,dir,norm);
      }

      if (bfly && bfly_alloc) { delete [] bfly; }
   }

   return(fft);
}

double *FFT (double *fft, double *v, const int n, const int dir, const int norm)
{
   if (!fft) fft = ( (v && (n>0)) ? new double[2*n] : 0 );

   if (v && fft && (n>0))
   {
      int *bfly=0, bfly_alloc=0;

      switch(n)
      {
         case (  4) : { bfly=bfly_4;   break; }
         case (  8) : { bfly=bfly_8;   break; }
         case ( 16) : { bfly=bfly_16;  break; }
         case ( 32) : { bfly=bfly_32;  break; }
         case ( 64) : { bfly=bfly_64;  break; }
         case (128) : { bfly=bfly_128; break; }
         case (256) : { bfly=bfly_256; break; }
      }

      if (!bfly) { bfly = FFT_Butterfly(n); bfly_alloc=1; }

      if (bfly)
      {
         FFT_Load      (fft,v,bfly,n);
         FFT_Transform (fft,n,dir,norm);
      }

      if (bfly && bfly_alloc) { delete [] bfly; }
   }

   return(fft);
}

