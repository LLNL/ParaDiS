#pragma once

#ifndef _PDS_MOMENT_H
#define _PDS_MOMENT_H

//-----------------------------------------------------------------------------------------------

extern float  Min      (const float  *v, const int vn);
extern double Min      (const double *v, const int vn);

extern float  Max      (const float  *v, const int vn);
extern double Max      (const double *v, const int vn);

extern void   MinMax   (const float  *v, const int vn, float  & min, float  & max);
extern void   MinMax   (const double *v, const int vn, double & min, double & max);

extern float  Mean     (const float  *v, const int vn);
extern double Mean     (const double *v, const int vn);

extern float  Sdev     (const float  *v, const int vn);
extern double Sdev     (const double *v, const int vn);
extern float  Sdev     (const float  *v, const int vn, const float  mean);
extern double Sdev     (const double *v, const int vn, const double mean);

extern int    Sdev_Cnt (const float  *v, const int vn, const float  smax);
extern int    Sdev_Cnt (const double *v, const int vn, const double smax);

extern void   Moment   (const float  *v, const int vn, float  & min, float  & max, float  & mean, float  & range                             );
extern void   Moment   (const float  *v, const int vn, float  & min, float  & max, float  & mean, float  & range, float  & var, float  & sdev);
extern void   Moment   (const double *v, const int vn, double & min, double & max, double & mean, double & range                             );
extern void   Moment   (const double *v, const int vn, double & min, double & max, double & mean, double & range, double & var, double & sdev);

extern void   Moment   (float  *mmt, const float  *v, const int vn);
extern void   Moment   (double *mmt, const double *v, const int vn);

//-----------------------------------------------------------------------------------------------

#endif  // _PDS_MOMENT_H
