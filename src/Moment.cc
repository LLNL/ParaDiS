#include <math.h>

#include "Moment.h"

//------------------------------------------------------------------------------------

#define MIN                                 \
   if (v && (vn>0))                         \
   {                                        \
      min = v[0];                           \
      for (int i=0; (i<vn); i++)            \
         min = ( (min<v[i]) ? min : v[i] ); \
   }
 
#define MAX                                 \
   if (v && (vn>0))                         \
   {                                        \
      max = v[0];                           \
      for (int i=0; (i<vn); i++)            \
         max = ( (max>v[i]) ? max : v[i] ); \
   }
 
float  Min (const float  *v, const int vn) { float min=0.0; MIN; return(min); }
double Min (const double *v, const int vn) { float min=0.0; MIN; return(min); }

float  Max (const float  *v, const int vn) { float max=0.0; MAX; return(max); }
double Max (const double *v, const int vn) { float max=0.0; MAX; return(max); }

#undef MIN
#undef MAX

//------------------------------------------------------------------------------------

#define MINMAX                              \
   if (v && (vn>0))                         \
   {                                        \
      min = max = v[0];                     \
      for (int i=0; (i<vn); i++)            \
      {                                     \
         min = ( (min<v[i]) ? min : v[i] ); \
         max = ( (max>v[i]) ? max : v[i] ); \
      }                                     \
   }
 
void MinMax (const float  *v, const int vn, float  & min, float  & max) { MINMAX; }
void MinMax (const double *v, const int vn, double & min, double & max) { MINMAX; }

#undef MINMAX

//------------------------------------------------------------------------------------

#define MEAN                       \
   if (v && (vn>0))                \
   {                               \
      for (int i=0; (i<vn); i++)   \
         mean+=v[i];               \
                                   \
      mean/=vn;                    \
   }
 

float  Mean (const float  *v, const int vn) { float  mean=0.0; MEAN; return(mean); }
double Mean (const double *v, const int vn) { double mean=0.0; MEAN; return(mean); }

#undef MEAN

//------------------------------------------------------------------------------------

#define SDEV                             \
   if (v && (vn>0))                      \
   {                                     \
      for (int i=0; (i<vn); i++)         \
      { s = (v[i]-mean); sdev+=(s*s); }  \
                                         \
      sdev = sqrt(sdev/vn);              \
   }

float  Sdev (const float  *v, const int vn) { float  mean=Mean(v,vn), sdev=0.0, s; SDEV; return(sdev); }
double Sdev (const double *v, const int vn) { double mean=Mean(v,vn), sdev=0.0, s; SDEV; return(sdev); }

float  Sdev (const float  *v, const int vn, const float  mean) { float  sdev=0.0,s; SDEV; return(sdev); }
double Sdev (const double *v, const int vn, const double mean) { double sdev=0.0,s; SDEV; return(sdev); }

#undef SDEV

//------------------------------------------------------------------------------------

int Sdev_Cnt (const float  *v, const int vn, const float  smax) 
{ 
   int cnt=0; 

   if (v && (vn>0))
   {
      for (int i=0; (i<vn); i++)
         if (fabsf(v[i])>smax) cnt++;
   }

   return(cnt); 
}

int Sdev_Cnt (const double *v, const int vn, const double smax)
{ 
   int cnt=0; 

   if (v && (vn>0))
   {
      for (int i=0; (i<vn); i++) 
         if (fabs(v[i])>smax) cnt++;
   }

   return(cnt); 
}

//------------------------------------------------------------------------------------

#define MOMENT                           \
{                                        \
   min = max = mean = range = 0.0;       \
                                         \
   if (v && (vn>0))                      \
   {                                     \
      sum = min = max = v[0];            \
                                         \
      for (int i=1; (i<vn); i++)         \
      {                                  \
         s    = v[i];                    \
         min  = ( (s<min) ? s : min );   \
         max  = ( (s>max) ? s : max );   \
         sum += s;                       \
      }                                  \
                                         \
      mean  = (sum/vn);                  \
      range = (max-min);                 \
   }                                     \
}

void Moment (const float  *v, const int vn, float  & min, float  & max, float  & mean, float  & range) { float  sum, s; MOMENT; }
void Moment (const double *v, const int vn, double & min, double & max, double & mean, double & range) { double sum, s; MOMENT; }

#undef MOMENT

//------------------------------------------------------------------------------------

#define MOMENT                                       \
{                                                    \
   min = max = mean = range = var = sdev = 0.0;      \
                                                     \
   Moment(v,vn,min,max,mean,range);                  \
                                                     \
   if (v && (vn>0))                                  \
   {                                                 \
      for (int i=0; (i<vn); i++)                     \
      {                                              \
         s    = (v[i]-mean);                         \
         var += (s*s);                               \
      }                                              \
                                                     \
      var  = (var/vn);                               \
      sdev = sqrt(var);                              \
   }                                                 \
}

void Moment (const float  *v, const int vn, float  & min, float  & max, float  & mean, float  & range, float  & var, float  & sdev) { float   s; MOMENT; }
void Moment (const double *v, const int vn, double & min, double & max, double & mean, double & range, double & var, double & sdev) { double  s; MOMENT; }

#undef MOMENT

void Moment (float  *m, const float  *v, const int vn) { if (m) Moment(v,vn, m[0], m[1], m[2], m[3], m[4], m[5] ); }
void Moment (double *m, const double *v, const int vn) { if (m) Moment(v,vn, m[0], m[1], m[2], m[3], m[4], m[5] ); }

