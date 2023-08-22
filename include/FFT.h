#pragma once

#ifndef __PDS_FFT_H
#define __PDS_FFT_H

extern int    *FFT_Butterfly      (const int n);

extern void    FFT_Load           (float  *dst, const unsigned char   *src, const int *bfly, const int n);
extern void    FFT_Load           (float  *dst, const unsigned short  *src, const int *bfly, const int n);
extern void    FFT_Load           (float  *dst, const          short  *src, const int *bfly, const int n);
extern void    FFT_Load           (float  *dst, const unsigned int    *src, const int *bfly, const int n);
extern void    FFT_Load           (float  *dst, const          int    *src, const int *bfly, const int n);
extern void    FFT_Load           (float  *dst, const          float  *src, const int *bfly, const int n);

extern void    FFT_Load           (double *dst, const unsigned char   *src, const int *bfly, const int n);
extern void    FFT_Load           (double *dst, const unsigned short  *src, const int *bfly, const int n);
extern void    FFT_Load           (double *dst, const          short  *src, const int *bfly, const int n);
extern void    FFT_Load           (double *dst, const unsigned int    *src, const int *bfly, const int n);
extern void    FFT_Load           (double *dst, const          int    *src, const int *bfly, const int n);
extern void    FFT_Load           (double *dst, const          double *src, const int *bfly, const int n);

extern void    FFT_Load_Complex   (float  *dst, const float  *src, const int *bfly, const int n);
extern void    FFT_Load_Complex   (double *dst, const double *src, const int *bfly, const int n);

extern void    FFT_Load_Complex   (float  *dst, const float  *src, const int *bfly, const int n, const int m);
extern void    FFT_Load_Complex   (double *dst, const double *src, const int *bfly, const int n, const int m);

extern void    FFT_Power_Spectrum (float  *vec, const int n);
extern void    FFT_Power_Spectrum (double *vec, const int n);
extern void    FFT_Power_Spectrum (float  *dst, const float  *src, const int n);
extern void    FFT_Power_Spectrum (double *dst, const double *src, const int n);

extern void    FFT_Normalize      (float  *vec, const int n);
extern void    FFT_Normalize      (double *vec, const int n);

extern void    FFT_Transform      (float  *vec, const int n, const int dir, const int norm);
extern void    FFT_Transform      (double *vec, const int n, const int dir, const int norm);

extern float  *FFT                (float  *vec, const int n, const int dir, const int norm);
extern double *FFT                (double *vec, const int n, const int dir, const int norm);

extern float  *FFT                (float  *fft, float  *v, const int n, const int dir, const int norm);
extern double *FFT                (double *fft, double *v, const int n, const int dir, const int norm);
#endif 
