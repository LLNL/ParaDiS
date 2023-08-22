#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"

//---------------------------------------------------------------------------------------------------------
// Collection of matrix support routines available to ParaDiS...
//
// Note - the general conventions used within most of the new routines is consistent with
// most of the standard C/C++ libraray where the left hand side (LHS) arguments that receive
// output are first, followed by the arguments of the right hand side (RHS).
//
// e.g.
//     LHS    = RHS
//     result = function(rhs args);
//              function(lhs args, rhs args)
//
// Overall function summary...
//
// rslt   function(args...)                                                                                             brief description....
// ----   ------------------------------------------------------------------------------------------------------------- ---------------------------------------------------
// void   Matrix22_Zero(real8 a[2][2])                                                                                 : A = 0
// void   Matrix23_Zero(real8 a[2][3])                                                                                 : A = 0
// void   Matrix32_Zero(real8 a[3][2])                                                                                 : A = 0
// void   Matrix33_Zero(real8 a[3][3])                                                                                 : A = 0
// void   Matrix44_Zero(real8 a[4][4])                                                                                 : A = 0
// void   Matrix66_Zero(real8 a[6][6])                                                                                 : A = 0
// void   Matrix333_Zero(real8 a[3][3][3])                                                                             : A = 0
// void   MatrixMN_Zero(real8  *a, const int m , const int n)                                                          : A = 0
// void   MatrixMN_Zero(real8 **a, const int m , const int n)                                                          : A = 0
// void   MatrixMN_Zero(real8 **a, const int m0, const int m1, const int n0, const int n1)                             : A = 0
//
// void   MatrixMN_Set                                                                                                 : A = S
// void   Matrix22_Copy (real8 a[2][2], real8 b[2][2])                                                                 : A = B (2x2)
// void   Matrix33_Copy (real8 a[3][3], real8 b[3][3])                                                                 : A = B (3x3)
// void   Matrix44_Copy (real8 a[4][4], real8 b[4][4])                                                                 : A = B (4x4)
// void   Matrix66_Copy (real8 a[6][6], real8 b[6][6])                                                                 : A = B (6x6)
// void   MatrixMN_Copy                                                                                                : A = B (MxN)
// void   MatrixMN_Copy                                                                                                : A = B (MxN)
//
// (print utilities)
// void   Matrix22_Print (FILE *fd, const char *fmt, real8  a[2][2] )                                                  : fd << print(A) (2x2)
// void   Matrix33_Print (FILE *fd, const char *fmt, real8  a[3][3] )                                                  : fd << print(A) (3x3)
// void   Matrix44_Print (FILE *fd, const char *fmt, real8  a[4][4] )                                                  : fd << print(A) (4x4)
// void   Matrix66_Print (FILE *fd, const char *fmt, real8  a[6][6] )                                                  : fd << print(A) (6x6)
// void   MatrixMN_Print (FILE *fd, const char *fmt, real8  *a, const int m , const int n)                             : fd << print(A) (MxN) (zero indexed)
// void   MatrixMN_Print (FILE *fd, const char *fmt, real8 **a, const int m , const int n)                             : fd << print(A) (MxN) (zero indexed, iliffe)
// void   MatrixMN_Print (FILE *fd, const char *fmt, real8 **a, const int m0, const int m1,                            :
//                                                              const int n0, const int n1)                            : fd << print(A) (MxN) (general     . iliffe)
//
// void   Matrix22_Identity(real8 a[2][2])                                                                             : A = I (2x2)
// void   Matrix33_Identity(real8 a[3][3])                                                                             : A = I (3x3)
// void   Matrix44_Identity(real8 a[4][4])                                                                             : A = I (4x4)
// void   Matrix66_Identity(real8 a[6][6])                                                                             : A = I (6x6)
// void   MatrixMN_Identity(real8  *a, const int m , const int n )                                                     : A = I (MxN) (zero indexed)
// void   MatrixMN_Identity(real8 **a, const int m , const int n )                                                     : A = I (MxN) (zero indexed, iliffe)
// void   MatrixMN_Identity(real8 **a, const int m0, const int m1, const int n0, const int n1)                         : A = I (MxN) (general     , iliffe)
//
// void   Matrix22_Threshold  (real8   a[2][2],                           const real8 min);                            : A = ( fabs(a[i][j])<min ? 0.0 : a[i][j] )
// void   Matrix33_Threshold  (real8   a[3][3],                           const real8 min);                            : A = ( fabs(a[i][j])<min ? 0.0 : a[i][j] )
// void   Matrix44_Threshold  (real8   a[4][4],                           const real8 min);                            : A = ( fabs(a[i][j])<min ? 0.0 : a[i][j] )
// void   MatrixMN_Threshold  (real8  *a      , const int m, const int n, const real8 min);                            : A = ( fabs(a[i][j])<min ? 0.0 : a[i][j] )
// void   MatrixMN_Threshold  (real8 **a      , const int m, const int n, const real8 min);                            : A = ( fabs(a[i][j])<min ? 0.0 : a[i][j] )
//
// void   Matrix22_Cap        (real8   a[2][2],                           const real8 max);                            : A = ( a[i][j]<max ? a[i][j] : max )
// void   Matrix33_Cap        (real8   a[3][3],                           const real8 max);                            : A = ( a[i][j]<max ? a[i][j] : max )
// void   Matrix44_Cap        (real8   a[4][4],                           const real8 max);                            : A = ( a[i][j]<max ? a[i][j] : max )
// void   MatrixMN_Cap        (real8  *a      , const int m, const int n, const real8 max);                            : A = ( a[i][j]<max ? a[i][j] : max )
// void   MatrixMN_Cap        (real8 **a      , const int m, const int n, const real8 max);                            : A = ( a[i][j]<max ? a[i][j] : max )
//
// void   Matrix22_Constrain  (real8   a[2][2],                           const real8 min, const real8 max);           : A = ( a[i][j]<min ? min : ( a[i][j]>max ? max : a[i][j] ) )
// void   Matrix33_Constrain  (real8   a[3][3],                           const real8 min, const real8 max);           : A = ( a[i][j]<min ? min : ( a[i][j]>max ? max : a[i][j] ) )
// void   Matrix44_Constrain  (real8   a[4][4],                           const real8 min, const real8 max);           : A = ( a[i][j]<min ? min : ( a[i][j]>max ? max : a[i][j] ) )
// void   MatrixMN_Constrain  (real8  *a      , const int m, const int n, const real8 min, const real8 max);           : A = ( a[i][j]<min ? min : ( a[i][j]>max ? max : a[i][j] ) )
// void   MatrixMN_Constrain  (real8 **a      , const int m, const int n, const real8 min, const real8 max);           : A = ( a[i][j]<min ? min : ( a[i][j]>max ? max : a[i][j] ) )
//
// void   Matrix22_Abs        (real8   a[2][2]                          );                                             : A = Abs(A)    a[i][j] = fabs(a[i][j])
// void   Matrix33_Abs        (real8   a[3][3]                          );                                             : A = Abs(A)    a[i][j] = fabs(a[i][j])
// void   Matrix44_Abs        (real8   a[4][4]                          );                                             : A = Abs(A)    a[i][j] = fabs(a[i][j])
// void   MatrixMN_Abs        (real8  *a      , const int m, const int n);                                             : A = Abs(A)    a[i]    = fabs(a[i]   )
// void   MatrixMN_Abs        (real8 **a      , const int m, const int n);                                             : A = Abs(A)    a[i][j] = fabs(a[i][j])
//
// void   Matrix22_Sign       (real8   a[2][2]                          );                                             : A = Sign(A)   a[i][j] = ( a[i][j]<0 ? -1 : ( a[i][j]>0 ? +1 : 0 ) )
// void   Matrix33_Sign       (real8   a[3][3]                          );                                             : A = Sign(A)   a[i][j] = ( a[i][j]<0 ? -1 : ( a[i][j]>0 ? +1 : 0 ) )
// void   Matrix44_Sign       (real8   a[4][4]                          );                                             : A = Sign(A)   a[i][j] = ( a[i][j]<0 ? -1 : ( a[i][j]>0 ? +1 : 0 ) )
// void   MatrixMN_Sign       (real8  *a      , const int m, const int n);                                             : A = Sign(A)   a[i]    = ( a[i]   <0 ? -1 : ( a[i]   >0 ? +1 : 0 ) )
// void   MatrixMN_Sign       (real8 **a      , const int m, const int n);                                             : A = Sign(A)   a[i][j] = ( a[i][j]<0 ? -1 : ( a[i][j]>0 ? +1 : 0 ) )
//
// void   Matrix22_Residual   (real8   a[2][2], real8   b[2][2] );                                                     : A = Abs(A-B)
// void   Matrix33_Residual   (real8   a[3][3], real8   b[3][3] );                                                     : A = Abs(A-B)
// void   Matrix44_Residual   (real8   a[4][4], real8   b[4][4] );                                                     : A = Abs(A-B)
// void   MatrixMN_Residual   (real8  *a      , real8  *b      , const int m, const int n);                            : A = Abs(A-B)
// void   MatrixMN_Residual   (real8 **a      , real8 **b      , const int m, const int n);                            : A = Abs(A-B)
//
// real8  Matrix22_RMS_Error  (real8   a[2][2], real8   b[2][2] );                                                     : rms = RMS_Error(A-B)
// real8  Matrix33_RMS_Error  (real8   a[3][3], real8   b[3][3] );                                                     : rms = RMS_Error(A-B)
// real8  Matrix44_RMS_Error  (real8   a[4][4], real8   b[4][4] );                                                     : rms = RMS_Error(A-B)
// real8  MatrixMN_RMS_Error  (real8  *a      , real8  *b      , const int m, const int n);                            : rms = RMS_Error(A-B)
// real8  MatrixMN_RMS_Error  (real8 **a      , real8 **b      , const int m, const int n);                            : rms = RMS_Error(A-B)
//
// void   Matrix22_Rand       (real8   a[2][2]                          , const real8 min=0.0, const real8 max=1.0);   : A = a[i] = srand(min,max)
// void   Matrix33_Rand       (real8   a[3][3]                          , const real8 min=0.0, const real8 max=1.0);   : A = a[i] = srand(min,max)
// void   Matrix44_Rand       (real8   a[4][4]                          , const real8 min=0.0, const real8 max=1.0);   : A = a[i] = srand(min,max)
// void   MatrixMN_Rand       (real8  *a      , const int m, const int n, const real8 min=0.0, const real8 max=1.0);   : A = a[i] = srand(min,max)
// void   MatrixMN_Rand       (real8 **a      , const int m, const int n, const real8 min=0.0, const real8 max=1.0);   : A = a[i] = srand(min,max)
//
// (rotation)
// void   Matrix22_Rotate_Z(real8 a[2][2], const real8 rz)                                             : A = RotZ(rz)          (2x2)
//
// void   Matrix33_Rotate_X(real8 a[3][3], const real8 rx)                                             : A = RotX(rx)          (3x3)
// void   Matrix33_Rotate_Y(real8 a[3][3], const real8 ry)                                             : A = RotY(ry)          (3x3)
// void   Matrix33_Rotate_Z(real8 a[3][3], const real8 rz)                                             : A = RotZ(rz)          (3x3)
// void   Matrix33_Rotate  (real8 a[3][3], const real8 rx, const real8 ry, const real8 rz)             : A = Rot (rx,ry,rz)    (3x3) (euler angles)
//
// void   Matrix44_Rotate_X(real8 a[4][4], const real8 rx)                                             : A = RotX(rx)          (4x4)
// void   Matrix44_Rotate_Y(real8 a[4][4], const real8 ry)                                             : A = RotY(ry)          (4x4)
// void   Matrix44_Rotate_Z(real8 a[4][4], const real8 rz)                                             : A = RotZ(rz)          (4x4)
// void   Matrix44_Rotate  (real8 a[4][4], const real8 rx, const real8 ry, const real8 rz)             : A = Rot (rx,ry,rz,1)  (4x4) (euler angles)
//
// (scale)
// void   Matrix22_Diag(real8 a[2][2], const real8 s)                                                 : A = Scale(s ,s )      (2x2)
// void   Matrix22_Diag(real8 a[2][2], const real8 sx, const real8 sy)                                : A = Scale(sx,sy)      (2x2)
// void   Matrix33_Diag(real8 a[3][3], const real8 s)                                                 : A = Scale(s ,s ,s )   (3x3)
// void   Matrix33_Diag(real8 a[3][3], const real8 sx, const real8 sy, const real8 sz)                : A = Scale(sx,sy,sz)   (3x3)
// void   Matrix44_Diag(real8 a[4][4], const real8 s)                                                 : A = Scale(s ,s ,s ,1) (4x4)
// void   Matrix44_Diag(real8 a[4][4], const real8 sx, const real8 sy, const real8 sz)                : A = Scale(sx,sy,sz,1) (4x4)
//
// (add, subtract, multiply)
// void   Matrix22_Add (real8 a[2][2], real8 b[2][2], real8 c[2][2]);                                  : A = B+C   (2x2)
// void   Matrix22_Sub (real8 a[2][2], real8 b[2][2], real8 c[2][2]);                                  : A = B-C   (2x2)
// void   Matrix22_Mul (real8 a[2][2], real8 b[2][2], real8 c[2][2]);                                  : A = B*C   (2x2) (A cannot be either B or C)
// void   Matrix22_Mul (real8 a[2][2], real8 b[2][2], real8 c[2][2], real8 d[2][2]);                   : A = B*C*D (2x2) (A cannot be either B,C, or D)
// void   Matrix22_Mul (real8 a[2][2], real8 b[2][3], real8 c[3][2]);                                  : A = B*C   (2x2) (inner product of 2x3 and 3x2 matrices)
//
// void   Matrix33_Add (real8 a[3][3], real8 b[3][3], real8 c[3][3])                                   : A = B+C   (3x3)
// void   Matrix33_Sub (real8 a[3][3], real8 b[3][3], real8 c[3][3])                                   : A = B-C   (3x3)
// void   Matrix33_Mul (real8 a[3][3], real8 b[3][3], real8 c[3][3])                                   : A = B*C   (3x3) (A cannot be either B or C)
// void   Matrix33_Mul (real8 a[3][3], real8 b[3][3], real8 c[3][3], real8 d[3][3]);                   : A = B*C*D (3x3) (A cannot be either B,C, or D)
//
// void   Matrix44_Add (real8 a[4][4], real8 b[4][4], real8 c[4][4])                                   : A = B+C   (4x4)
// void   Matrix44_Sub (real8 a[4][4], real8 b[4][4], real8 c[4][4])                                   : A = B-C   (4x4)
// void   Matrix44_Mul (real8 a[4][4], real8 b[4][4], real8 c[4][4])                                   : A = B*C   (4x4) (A cannot be either B or C)
// void   Matrix44_Mul (real8 a[4][4], real8 b[4][4], real8 c[4][4], real8 d[4][4]);                   : A = B*C*D (4x4) (A cannot be either B,C, or D)
//
// void   Matrix66_Add (real8 a[6][6], real8 b[6][6], real8 c[6][6])                                   : A = B+C   (6x6)
// void   Matrix66_Sub (real8 a[6][6], real8 b[6][6], real8 c[6][6])                                   : A = B-C   (6x6)
// void   Matrix66_Mul (real8 a[6][6], real8 b[6][6], real8 c[6][6])                                   : A = B*C   (6x6) (A cannot be either B or C)
//
// void   MatrixMN_Add (real8  *a, real8  *b, real8  *c, const int m, const int n);                    : A = B+C   (MxN) (dense , zero-indexed)
// void   MatrixMN_Sub (real8  *a, real8  *b, real8  *c, const int m, const int n);                    : A = B-C   (MxN) (dense , zero-indexed)
// void   MatrixMN_Mul (real8  *a, real8  *b, real8  *c, const int m, const int n, const int p);       : A = B*C   (MxN) (dense , zero-indexed) (A cannot be either B or C)
//
// void   MatrixMN_Add (real8 **a, real8 **b, real8 **c, const int m, const int n);                    : A = B+C   (MxN) (iliffe, zero-indexed)
// void   MatrixMN_Sub (real8 **a, real8 **b, real8 **c, const int m, const int n);                    : A = B-C   (MxN) (iliffe, zero-indexed)
// void   MatrixMN_Mul (real8 **a, real8 **b, real8 **c, const int m, const int n, const int p);       : A = B*C   (MxN) (iliffe, zero-indexed) (A cannot be either B or C)
//
// (matrix vector operations)
// void   Matrix22_Vmul (real8 y[2]   , real8 a[2][2], real8 x[2]   );                                 : Y = A*X  (2x2 matrix vector post multiply)
// void   Matrix22_Vmul (real8 y[2]   , real8 x[2]   , real8 a[2][2]);                                 : Y = X*A  (2x2 vector matrix pre  multiply)
// void   Matrix23_Vmul (real8 y[2]   , real8 a[2][3], real8 x[3]   );                                 : Y = A*X  (2x3 matrix vector post multiply of 3x1 vector)
// void   Matrix32_Vmul (real8 y[3]   , real8 a[3][2], real8 x[2]   );                                 : Y = A*X  (3x2 matrix vector post multiply of 2x1 vector)
// void   Matrix33_Vmul (real8 y[3]   , real8 a[3][3], real8 x[3]   );                                 : Y = A*X  (3x3 matrix vector post multiply)
// void   Matrix33_Vmul (real8 y[3]   , real8 x[3]   , real8 a[3][3]);                                 : Y = X*A  (3x3 vector matrix pre  multiply)
// void   Matrix44_Vmul (real8 y[3]   , real8 a[4][4], real8 x[3]   );                                 : Y = A*X  (4x4 matrix vector post multiply)
// void   Matrix44_Vmul (real8 y[3]   , real8 x[3]   , real8 a[4][4]);                                 : Y = X*A  (4x4 vector matrix pre  multiply)
// void   Matrix43_Vmul (real8 y[4]   , real8 a[4][3], real8 x[3]   );                                 : Y = A*X  (4x4 matrix vector post multiply of 3x1 vector, assumed 1)
// void   Matrix63_Vmul (real8 y[6]   , real8 a[6][3], real8 x[3]   );                                 : Y = A*X  (6x3 matrix vector post multiply of 3x1 vector)
//
// void   Matrix33_VVt_Mul (real8 a[3][3], real8 v[3]);                                                : A = V*Vt (3x3 outer product of 3x1 vector and its transpose)
// void   Matrix33_VVt_Mul (real8 a[3][3], real8 v[3], real8 vt[3]);                                   : A = V*Vt (3x3 outer product of 3x1 vector and 1x3 vector)
//
// (transpose)
// void   Matrix22_Transpose (real8 a[2][2])                                                           : A = T(A)  (2x2) (in-place)
// void   Matrix33_Transpose (real8 a[3][3])                                                           : A = T(A)  (3x3) (in-place)
// void   Matrix44_Transpose (real8 a[4][4])                                                           : A = T(A)  (4x4) (in-place)
// void   Matrix66_Transpose (real8 a[6][6])                                                           : A = T(A)  (6x6) (in-place)
// void   Matrix22_Transpose (real8 a[2][2], real8 b[2][2])                                            : A = T(B)  (2x2) (A cannot be B)
// void   Matrix33_Transpose (real8 a[3][3], real8 b[3][3])                                            : A = T(B)  (3x3) (A cannot be B)
// void   Matrix44_Transpose (real8 a[4][4], real8 b[4][4])                                            : A = T(B)  (4x4) (A cannot be B)
// void   Matrix66_Transpose (real8 a[6][6], real8 b[6][6])                                            : A = T(B)  (6x6) (A cannot be B)
// void   Matrix23_Transpose (real8 a[3][2], real8 b[2][3] )                                           : A = T(B)  A(3x2), B(2x3) (A cannot be B)
// void   Matrix32_Transpose (real8 a[2][3], real8 b[3][2] )                                           : A = T(B)  A(2x3), B(3x2) (A cannot be B)
// void   MatrixMN_Transpose (real8  *a,            const int m, const int n)                          : A = T(A)  A(NxM), B(MxN) (in-place, zero-indexed)
// void   MatrixMN_Transpose (real8  *a,  real8 *b, const int m, const int n)                          : A = T(B)  A(NxM), B(MxN) (          zero-indexed,         A cannot be B)
// void   MatrixMN_Transpose (real8 **a, real8 **b, const int m, const int n)                          : A = T(B)  A(NxM), B(MxN) (          zero-indexed, iliffe, A cannot be B)
//
// (comparison)
// int    Matrix22_Near (const real8 a[2][2], const real8 s      , const real8 tol)                    : returns ( A==s ? 1 : 0 )  (2x2)
// int    Matrix22_Near (const real8 a[2][2], const real8 b[2][2], const real8 tol)                    : returns ( A==B ? 1 : 0 )  (2x2)
// int    Matrix33_Near (const real8 a[3][3], const real8 s      , const real8 tol)                    : returns ( A==s ? 1 : 0 )  (3x3)
// int    Matrix33_Near (const real8 a[3][3], const real8 b[3][3], const real8 tol)                    : returns ( A==B ? 1 : 0 )  (3x3)
// int    Matrix44_Near (const real8 a[4][4], const real8 s      , const real8 tol)                    : returns ( A==s ? 1 : 0 )  (4x4)
// int    Matrix44_Near (const real8 a[4][4], const real8 b[4][4], const real8 tol)                    : returns ( A==B ? 1 : 0 )  (4x4)
// int    Matrix66_Near (const real8 a[6][6], const real8 s      , const real8 tol)                    : returns ( A==s ? 1 : 0 )  (6x6)
// int    Matrix66_Near (const real8 a[6][6], const real8 b[6][6], const real8 tol)                    : returns ( A==B ? 1 : 0 )  (6x6)
// int    MatrixMN_Near (const real8  *a, const real8   s, const int m , const int n , real8 tol);     : returns ( A==s ? 1 : 0 )  (MxN,dense , zero-indexed)
// int    MatrixMN_Near (      real8 **a, const real8   s, const int m , const int n , real8 tol);     : returns ( A==s ? 1 : 0 )  (MxN,iliffe, zero-indexed)
// int    MatrixMN_Near (      real8 **a, const real8   s, const int m0, const int m1,
//                                                         const int n0, const int n1, real8 tol);     : returns ( A==s ? 1 : 0 )  (MxN,iliffe, general     )
// int    MatrixMN_Near (const real8  *a, const real8  *b, const int m , const int n , real8 tol);     : returns ( A==B ? 1 : 0 )  (MxN,dense , zero-indexed)
// int    MatrixMN_Near (      real8 **a,       real8 **b, const int m , const int n , real8 tol);     : returns ( A==B ? 1 : 0 )  (MxN,iliffe, zero-indexed)
// int    MatrixMN_Near (      real8 **a,       real8 **b, const int m0, const int m1,
//                                                         const int n0, const int n1, real8 tol);     : returns ( A==B ? 1 : 0 )  (MxN,iliffe, general     )
//
// (determinants, norms, condition, adjoints, inverse)
// real8  Matrix22_Det       (real8 a[2][2])                                                           : returns  |A|      (2x2)
// real8  Matrix22_Norm      (real8 a[2][2])                                                           : returns ||A||     (2x2)
// real8  Matrix22_Condition (real8 a[2][2])                                                           : returns Cond(A)   (2x2)
// void   Matrix22_Adjoint   (real8 a[2][2], real8 b[2][2])                                            : A = Adj(B)        (2x2)
// int    Matrix22_Inverse   (real8 a[2][2], real8 b[2][2])                                            : A = Inv(B)        (2x2) (returns 1 (success), -1 (error))
//
// real8  Matrix33_Det       (real8 a[3][3])                                                           : returns  |A|      (3x3)
// real8  Matrix33_Norm      (real8 a[3][3])                                                           : returns ||A||     (3x3)
// real8  Matrix33_Condition (real8 a[3][3])                                                           : returns Cond(A)   (3x3)
// void   Matrix33_Adjoint   (real8 a[3][3], real8 b[3][3])                                            : A = Adj(B)        (3x3)
// int    Matrix33_Inverse   (real8 a[3][3], real8 b[3][3])                                            : A = Inv(B)        (3x3) (returns 1 (success), -1 (error))
//
// real8  Matrix44_Det       (real8 a[4][4])                                                           : returns  |A|      (4x4)
// real8  Matrix44_Norm      (real8 a[4][4])                                                           : returns ||A||     (4x4)
// real8  Matrix44_Condition (real8 a[4][4])                                                           : returns Cond(A)   (4x4)
// void   Matrix44_Adjoint   (real8 a[4][4], real8 b[4][4])                                            : A = Adj(B)        (4x4)
// int    Matrix44_Inverse   (real8 a[4][4], real8 b[4][4])                                            : A = Inv(B)        (4x4) (returns 1 (success), -1 (error))
//
// (allocation, deallocation)
// real8   *MatrixMN_Alloc  (const int m , const int n )                                               : Allocates and returns matrix (MxN, dense , zero indexed)
// real8   *MatrixMN_Alloc  (const int m0, const int m1, const int n0, const int n1)                   : Allocates and returns matrix (MxN, dense , general     )
// real8  **MatrixMN_Iliffe (const int m , const int n )                                               : Allocates and returns matrix (MxN, iliffe, zero indexed)
// real8  **MatrixMN_Iliffe (const int m0, const int m1, const int n0, const int n1)                   : Allocates and returns matrix (MxN, iliffe, general     )
// real8  **MatrixMN_Iliffe (const real8 *mtx, const int m , const int n )                             : Allocates and returns iliffe vector, attaches to existing dense matrix  (MxN, iliffe, zero indexed)
// real8  **MatrixMN_Iliffe (const real8 *mtx, const int m0, const int m1,
//                                             const int n0, const int n1)                             : Allocates and returns iliffe vector, attaches to existing dense matrix  (MxN, iliffe, general     )
//
// void     MatrixMN_Free   (real8  *a)                                                                : releases matrix                   (MxN, dense , zero indexed)
// void     MatrixMN_Free   (real8 **a)                                                                : releases matrix and iliffe vector (MxN, iliffe, zero indexed)
// void     MatrixMN_Free   (real8 **a, real8 *am)                                                     : releases matrix and iliffe vector (MxN, iliffe, zero indexed)
// void     MatrixMN_Free   (real8 **a,            const int m0, const int m1,
//                                                 const int n0, const int n1)                         : releases matrix and iliffe vector (MxN, iliffe, general     )
// void     MatrixMN_Free   (real8 **a, real8 *am, const int m0, const int m1,
//                                                 const int n0, const int n1)                         : releases matrix and iliffe vector (MxN, iliffe, general     )
//
// (statistical)
// real8  MatrixMN_Min    (const real8 *a, const int m, const int n)                                   : returns Min  (A)                (MxN, dense, zero indexed)
// real8  MatrixMN_Max    (const real8 *a, const int m, const int n)                                   : returns Max  (A)                (MxN, dense, zero indexed)
// real8  MatrixMN_Mean   (const real8 *a, const int m, const int n)                                   : returns Mean (A)                (MxN, dense, zero indexed)
// real8  MatrixMN_Range  (const real8 *a, const int m, const int n)                                   : returns Range(A) (max-min)      (MxN, dense, zero indexed)
// real8  MatrixMN_Min    (int & iy, int & ix, const real8 *a, const int m, const int n)               : returns Min  (A) and location   (MxN, dense, zero indexed)
// real8  MatrixMN_Max    (int & iy, int & ix, const real8 *a, const int m, const int n)               : returns Max  (A) and location   (MxN, dense, zero indexed)
//
// void   MatrixMN_Moment (real8 & min, real8 & max, real8 & mean, real8 & range              , const real8 *a, const int m, const int n)  : returns statistical moment (min,max,mean,range     )(A) (MxN, dense, zero indexed)
// void   MatrixMN_Moment (real8 & min, real8 & max, real8 & mean, real8 & range, real8 & sdev, const real8 *a, const int m, const int n)  : returns statistical moment (min,max,mean,range,sdev)(A) (MxN, dense, zero indexed)
//
// (legacy routines)
// real8  MatrixMN_Norm(real8 **a, const int m0, const int m1, const int n0, const int n1)             :
// int    Matrix33Invert(real8 A[3][3], real8 B[3][3])                                                 :
// int    Matrix33Invert(real8 a[3][3], real8 b[3][3])                                                 :
// int    Matrix33Invert(real8 A[3][3], real8 B[3][3])                                                 :
// int    MatrixInvert(real8 *mat, real8 *invMat, int order, int lda)                                  :
// void   MatrixMult(real8 *a, int aRows, int aCols, int aLD,                                          :
// void   MatrixDblMultComplex(real8 *a, int aRows, int aCols, int aLD,                                :
// void   MatrixMultArb(real8 *a, int aCols, int aRowOffset,                                           :
// void   Matrix33Vector3Multiply(real8 A[3][3], real8 x[3], real8 y[3])                               :
// int    Matrix66_Inverse(real8 a[6][6], real8 b[6][6])                                               :
// void   Vec3TransposeAndMult(real8 vec[3], real8 result[3][3])                                       :
//
// (singular value decomposition)
// int    MatrixMN_SVD           (real8 **a, real8 **u, real8 *w, real8 *v, const int m, const int n);      : SVD(A) = U*W*V'
// int    MatrixMN_SVD           (real8  *a, real8  *u, real8 *w, real8 *v, const int m, const int n);      : SVD(A) = U*W*V'
//
// int    MatrixMN_SVD_Inverse   (real8  *a, const real8 *b, const int m, const int n);                     : A = SVD Inverse(B)
// int    MatrixMN_SVD_Inverse   (real8  *a,                 const int m, const int n);                     : A = SVD Inverse(A) (in-place)
//
// int    Matrix33_SVD           (real8 a[3][3] , real8 u [3][3], real8  w[3], real8 v[3][3] );             : SVD(A) = U*W*V'
// int    Matrix33_SVD_Inverse   (real8 a[3][3] );                                                          : A  = SVD Inverse(A) (in-place)
// int    Matrix33_SVD_Inverse   (real8 a[3][3] , real8 b[3][3] );                                          : A  = SVD Inverse(B)
//
// int    Matrix33_PseudoInverse (real8 a[3][3] , real8 b[3][3] );                                          : A  = PseudoInverse(B) (Nicolas Bertin version)
//
// real8  Matrix33_SVD_Condition (real8 a[3][3] );                                                          : returns SVD condition(A)
//
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

void Matrix22_Zero(real8 a[2][2])
{
    a[0][0]=0.0; a[0][1]=0.0;
    a[1][0]=0.0; a[1][1]=0.0;
}

void Matrix23_Zero(real8 a[2][3])
{
    a[0][0]=0.0; a[0][1]=0.0; a[0][2]=0.0;
    a[1][0]=0.0; a[1][1]=0.0; a[1][2]=0.0;
}

void Matrix32_Zero(real8 a[3][2])
{
    a[0][0]=0.0; a[0][1]=0.0;
    a[1][0]=0.0; a[1][1]=0.0;
    a[2][0]=0.0; a[2][1]=0.0;
}

void Matrix33_Zero(real8 a[3][3])
{
    a[0][0]=0.0; a[0][1]=0.0; a[0][2]=0.0;
    a[1][0]=0.0; a[1][1]=0.0; a[1][2]=0.0;
    a[2][0]=0.0; a[2][1]=0.0; a[2][2]=0.0;
}

void Matrix44_Zero(real8 a[4][4])
{
    a[0][0]=0.0; a[0][1]=0.0; a[0][2]=0.0; a[0][3]=0.0;
    a[1][0]=0.0; a[1][1]=0.0; a[1][2]=0.0; a[1][3]=0.0;
    a[2][0]=0.0; a[2][1]=0.0; a[2][2]=0.0; a[2][3]=0.0;
    a[3][0]=0.0; a[3][1]=0.0; a[3][2]=0.0; a[3][3]=0.0;
}

void Matrix66_Zero(real8 a[6][6])
{
    for (int i=0; (i<6); i++)
    {
       a[i][0]=0.0;
       a[i][1]=0.0;
       a[i][2]=0.0;
       a[i][3]=0.0;
       a[i][4]=0.0;
       a[i][5]=0.0;
    }
}

void Matrix333_Zero(real8 a[3][3][3])
{
    for(int i=0; (i<3); ++i)
    {
        a[i][0][0]=0.0; a[i][0][1]=0.0; a[i][0][2]=0.0;
        a[i][1][0]=0.0; a[i][1][1]=0.0; a[i][1][2]=0.0;
        a[i][2][0]=0.0; a[i][2][1]=0.0; a[i][2][2]=0.0;
    }
}

void MatrixMN_Zero(real8 *a, const int m, const int n)
{
    if (a && (m>0) && (n>0))
    {
        for(int i=0; (i<m*n); ++i) { a[i]=0.0; }
    }
}

void MatrixMN_Zero(real8 **a, const int m, const int n)
{
    if (a && (m>0) && (n>0))
    {
        for(int i=0; (i<m); ++i)
        for(int j=0; (j<n); ++j) { a[i][j]=0.0; }
    }
}

void MatrixMN_Zero(real8 **a, const int m0, const int m1, const int n0, const int n1)
{
    const int m = (m1-m0+1);
    const int n = (n1-n0+1);

    if (a && (m>0) && (n>0))
    {
        for(int i=m0; (i<=m1); ++i)
        for(int j=n0; (j<=n1); ++j) { a[i][j]=0.0; }
    }
}

//---------------------------------------------------------------------------------------------------------

void Matrix22_Set (real8 a[2][2], const real8 s )
{
   for (int i=0; i<2; ++i)
   for (int j=0; j<2; ++j) { a[i][j]=s; }
}

void Matrix23_Set (real8 a[2][3], const real8 s )
{
   for (int i=0; i<2; ++i)
   for (int j=0; j<3; ++j) { a[i][j]=s; }
}

void Matrix32_Set (real8 a[3][2], const real8 s )
{
   for (int i=0; i<3; ++i)
   for (int j=0; j<2; ++j) { a[i][j]=s; }
}

void Matrix33_Set (real8 a[3][3], const real8 s )
{
   for (int i=0; i<3; ++i)
   for (int j=0; j<3; ++j) { a[i][j]=s; }
}

void Matrix44_Set (real8 a[4][4], const real8 s )
{
   for (int i=0; i<4; ++i)
   for (int j=0; j<4; ++j) { a[i][j]=s; }
}

void Matrix66_Set (real8 a[6][6], const real8 s )
{
   for (int i=0; i<6; ++i)
   for (int j=0; j<6; ++j) { a[i][j]=s; }
}

void MatrixMN_Set
(
         real8  *a ,   ///< points to A matrix (dense)
   const real8   s ,   ///< set value
   const int     m ,   ///< number of rows
   const int     n     ///< number of columns
)
{
   if (a && (m>0) && (n>0))
   {
      for (int i=0; (i<m*n); ++i) { a[i]=s; }
   }
}

void MatrixMN_Set
(
         real8 **a ,   ///< points to A matrix (Iliffe matrix)
   const real8   s ,   ///< set value
   const int     m ,   ///< number of rows
   const int     n     ///< number of columns
)
{
   if (a && (m>0) && (n>0) )
   {
      for (int i=0; i<m; ++i)
      for (int j=0; j<n; ++j) { a[i][j]=s; }
   }
}

void MatrixMN_Set
(
         real8 **a ,   ///< points to A matrix (Iliffe matrix)
   const real8   s ,   ///< set value
   const int     m0,   ///< initial row index
   const int     m1,   ///< final   row index
   const int     n0,   ///< initial column index
   const int     n1    ///< final   column index
)
{
   const int m = (m1-m0+1);
   const int n = (n1-n0+1);

   if (a && (m>0) && (n>0) )
   {
      for (int i=m0; i<=m1; ++i)
      for (int j=n0; j<=n1; ++j) { a[i][j]=s; }
   }
}

void MatrixMN_Set
(
         real8  *a ,   ///< points to A matrix (dense)
   const real8  *s ,   ///< set values
   const int     m ,   ///< number of rows
   const int     n     ///< number of columns
)
{
   if (a && s && (m>0) && (n>0))
   {
      for (int i=0; (i<m*n); ++i) { a[i]=s[i]; }
   }
}

void MatrixMN_Set
(
         real8 **a ,   ///< points to A matrix (Iliffe matrix)
   const real8  *s ,   ///< set values
   const int     m0,   ///< initial row index
   const int     m1,   ///< final   row index
   const int     n0,   ///< initial column index
   const int     n1    ///< final   column index
)
{
   const int m = (m1-m0+1);
   const int n = (n1-n0+1);

   if (a && s && (m>0) && (n>0) )
   {
      for (int i=m0,k=0; i<=m1; ++i)
      for (int j=n0;     j<=n1; ++j,++k) { a[i][j]=s[k]; }
   }
}

void MatrixMN_Set
(
         real8 **a ,   ///< points to A matrix (Iliffe matrix)
         real8 **s ,   ///< points to A matrix (Iliffe matrix, assumed same geometry)
   const int     m0,   ///< initial row index
   const int     m1,   ///< final   row index
   const int     n0,   ///< initial column index
   const int     n1    ///< final   column index
)
{
   const int m = (m1-m0+1);
   const int n = (n1-n0+1);

   if (a && s && (m>0) && (n>0) )
   {
      for (int i=m0; i<=m1; ++i)
      for (int j=n0; j<=n1; ++j) { a[i][j]=s[i][j]; }
   }
}

//---------------------------------------------------------------------------------------------------------

void Matrix22_Copy (real8 a[2][2], real8 b[2][2])
{
   a[0][0]=b[0][0];
   a[0][1]=b[0][1];
   a[1][0]=b[1][0];
   a[1][1]=b[1][1];
}

void Matrix33_Copy (real8 a[3][3], real8 b[3][3])
{
   for (int i=0; (i<3); i++)
   {
      a[i][0]=b[i][0];
      a[i][1]=b[i][1];
      a[i][2]=b[i][2];
   }
}

void Matrix44_Copy (real8 a[4][4], real8 b[4][4])
{
   for (int i=0; (i<4); i++)
   {
      a[i][0]=b[i][0];
      a[i][1]=b[i][1];
      a[i][2]=b[i][2];
      a[i][3]=b[i][3];
   }
}

void Matrix66_Copy (real8 a[6][6], real8 b[6][6])
{
   for (int i=0; (i<6); i++)
   {
      a[i][0]=b[i][0];
      a[i][1]=b[i][1];
      a[i][2]=b[i][2];
      a[i][3]=b[i][3];
      a[i][4]=b[i][4];
      a[i][5]=b[i][5];
   }
}

// MatrixMN_Copy() : A=B
//
// Note - the dimensions of A and B are assumed identical.
//---------------------------------------------------------------------------------------------------------

void MatrixMN_Copy
(
         real8  *a,   ///< the destination matrix (dense)
   const real8  *b,   ///< the source      matrix (dense)
   const int     m,   ///< number of rows
   const int     n    ///< number of columns
)
{
   if (a && b && (a!=b) && (m>0) && (n>0) )
   {
      for (int i=0; i<(m*n); ++i) { a[i]=b[i]; }
   }
}

void MatrixMN_Copy
(
         real8 **a,   ///< the destination matrix (iliffe)
   const real8 **b,   ///< the source      matrix (iliffe)
   const int     m,   ///< number of rows
   const int     n    ///< number of columns
)
{
   if (a && b && (a!=b) && (m>0) && (n>0) )
   {
      for (int i=0; i<m; ++i)
      for (int j=0; j<n; ++j) { a[i][j]=b[i][j]; }
   }
}

// Matrix print utilities...
//---------------------------------------------------------------------------------------------------------

void Matrix22_Print (FILE  *fd, const char *fmt, real8  a[2][2] )
{
    if (fd && a)
    {
        for(int i=0; (i<2); ++i)
        {
           for(int j=0; (j<2); ++j)
              fprintf(fd,fmt,a[i][j]);

           fprintf(fd,"\n");
        }
    }
}

void Matrix33_Print (FILE  *fd, const char *fmt, real8  a[3][3] )
{
    if (fd && a)
    {
        for(int i=0; (i<3); ++i)
        {
           for(int j=0; (j<3); ++j)
              fprintf(fd,fmt,a[i][j]);

           fprintf(fd,"\n");
        }
    }
}

void Matrix44_Print (FILE  *fd, const char *fmt, real8  a[4][4] )
{
    if (fd && a)
    {
        for(int i=0; (i<4); ++i)
        {
           for(int j=0; (j<4); ++j)
              fprintf(fd,fmt,a[i][j]);

           fprintf(fd,"\n");
        }
    }
}

void Matrix66_Print (FILE  *fd, const char *fmt, real8  a[6][6] )
{
    if (fd && a)
    {
        for(int i=0; (i<6); ++i)
        {
           for(int j=0; (j<6); ++j)
              fprintf(fd,fmt,a[i][j]);

           fprintf(fd,"\n");
        }
    }
}

void MatrixMN_Print(FILE *fd, const char *fmt, real8 *a, const int m, const int n)
{
    if (fd && a && (m>0) && (n>0))
    {
        for(int i=0,k=0; (i<m); ++i)
        {
           for(int j=0;     (j<n); ++j, ++k)
              fprintf(fd,fmt,a[k]);

           fprintf(fd,"\n");
        }
    }
}

void MatrixMN_Print(FILE *fd, const char *fmt, real8 **a, const int m, const int n)
{
    if (fd && a && (m>0) && (n>0))
    {
        for(int i=0; (i<m); ++i)
        {
           for(int j=0; (j<n); ++j)
              fprintf(fd,fmt,a[i][j]);

           fprintf(fd,"\n");
        }
    }
}

void MatrixMN_Print(FILE *fd, const char *fmt, real8 **a, const int m0, const int m1, const int n0, const int n1)
{
    const int m = (m1-m0+1);
    const int n = (n1-n0+1);

    if (fd && a && (m>0) && (n>0))
    {
        for(int i=m0; (i<=m1); ++i)
        {
           for(int j=n0; (j<=n1); ++j)
              fprintf(fd,fmt,a[i][j]);

           fprintf(fd,"\n");
        }

    }
}

// Matrix22_Identity(m)     : m = | 1 0 |
//                                | 0 1 |
//
// Matrix33_Identity(m)     : m = | 1 0 0 |
//                                | 0 1 0 |
//                                | 0 0 1 |
//
// Matrix44_Identity(m)     : m = | 1 0 0 0 |
//                                | 0 1 0 0 |
//                                | 0 0 1 0 |
//                                | 0 0 0 1 |
//
// MatrixMN_Identity(m)     : m = | 1 0 0 . 0 |
//                                | 0 1 0 . 0 |
//                                | 0 0 1 . 0 |
//                                | . . . 1 0 |
//                                | 0 0 0 0 1 |
//
//---------------------------------------------------------------------------------------------------------

void Matrix22_Identity(real8 a[2][2])
{
    a[0][0]=1.0; a[0][1]=0.0;
    a[1][0]=0.0; a[1][1]=1.0;
}

void Matrix33_Identity(real8 a[3][3])
{
    a[0][0]=1.0; a[0][1]=0.0; a[0][2]=0.0;
    a[1][0]=0.0; a[1][1]=1.0; a[1][2]=0.0;
    a[2][0]=0.0; a[2][1]=0.0; a[2][2]=1.0;
}

void Matrix44_Identity(real8 a[4][4])
{
    a[0][0]=1.0; a[0][1]=0.0; a[0][2]=0.0; a[0][3]=0.0;
    a[1][0]=0.0; a[1][1]=1.0; a[1][2]=0.0; a[1][3]=0.0;
    a[2][0]=0.0; a[2][1]=0.0; a[2][2]=1.0; a[2][3]=0.0;
    a[3][0]=0.0; a[3][1]=0.0; a[3][2]=0.0; a[3][3]=1.0;
}

void Matrix66_Identity(real8 a[6][6])
{
   for (int i=0; (i<6); i++)
   for (int j=0; (j<6); j++) { a[i][j] = ( (i==j) ? 1.0 : 0.0 ); }
}

void MatrixMN_Identity(real8 *a, const int m, const int n)
{
   if (a && (m>0) && (n>0))
   {
      for(int i=0; (i<m*n); ++i)      { a[i]=0.0; }   // A=0
      for(int i=0; (i<m*n); i+=(n+1)) { a[i]=1.0; }   // (set diagonals to 1)
   }
}

void MatrixMN_Identity(real8 **a, const int m, const int n)
{
   if (a && (m>0) && (n>0))
   {
      for(int i=0; (i<m); ++i)
      for(int j=0; (j<n); ++j) { a[i][j] = ( (i==j) ? 1.0 : 0.0 ); }
   }
}

void MatrixMN_Identity(real8 **a, const int m0, const int m1, const int n0, const int n1)
{
   const int m = (m1-m0+1);
   const int n = (n1-n0+1);

   if (a && (m>0) && (n>0))
   {
      for(int i=0; (i<m); ++i)
      for(int j=0; (j<n); ++j) { a[m0+i][n0+j]= ( (i==j) ? 1.0 : 0.0 ); }
   }
}

//---------------------------------------------------------------------------------------------------------

void Matrix22_Threshold(real8 a[2][2], const real8 min)
{
   for (int i=0; (i<2); ++i)
   for (int j=0; (j<2); ++j)
   { a[i][j] = ( (fabs(a[i][j])<min) ? 0.0 : a[i][j] ); }
}

void Matrix33_Threshold(real8 a[3][3], const real8 min)
{
   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j)
   { a[i][j] = ( (fabs(a[i][j])<min) ? 0.0 : a[i][j] ); }
}

void Matrix44_Threshold(real8 a[4][4], const real8 min)
{
   for (int i=0; (i<4); ++i)
   for (int j=0; (j<4); ++j)
   { a[i][j] = ( (fabs(a[i][j])<min) ? 0.0 : a[i][j] ); }
}

void MatrixMN_Threshold (real8  *a, const int m, const int n, const real8 min)
{
   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<(m*n)); ++i)
      { a[i] = ( (fabs(a[i])<min) ? 0.0 : a[i] ); }
   }
}

void MatrixMN_Threshold (real8 **a, const int m, const int n, const real8 min)
{
   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<m); ++i)
      for (int j=0; (j<n); ++j)
      { a[i][j] = ( (fabs(a[i][j])<min) ? 0.0 : a[i][j] ); }
   }
}

//---------------------------------------------------------------------------------------------------------

void Matrix22_Cap(real8 a[2][2], const real8 max)
{
   for (int i=0; (i<2); ++i)
   for (int j=0; (j<2); ++j)
   { a[i][j] = ( (a[i][j]<max) ? a[i][j] : max ); }
}

void Matrix33_Cap(real8 a[3][3], const real8 max)
{
   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j)
   { a[i][j] = ( (a[i][j]<max) ? a[i][j] : max ); }
}

void Matrix44_Cap(real8 a[4][4], const real8 max)
{
   for (int i=0; (i<4); ++i)
   for (int j=0; (j<4); ++j)
   { a[i][j] = ( (a[i][j]<max) ? a[i][j] : max ); }
}

void MatrixMN_Cap (real8  *a, const int m, const int n, const real8 max)
{
   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<(m*n)); ++i)
      { a[i] = ( (a[i]<max) ? a[i] : max ); }
   }
}

void MatrixMN_Cap (real8 **a, const int m, const int n, const real8 max)
{
   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<m); ++i)
      for (int j=0; (j<n); ++j)
      { a[i][j] = ( (a[i][j]<max) ? a[i][j] : max ); }
   }
}

//---------------------------------------------------------------------------------------------------------

void Matrix22_Constrain(real8 a[2][2], const real8 min, const real8 max)
{
   for (int i=0; (i<2); ++i)
   for (int j=0; (j<2); ++j)
   { a[i][j] = ( (a[i][j]<min) ? min : ( (a[i][j]<max) ? a[i][j] : max ) ); }

}

void Matrix33_Constrain(real8 a[3][3], const real8 min, const real8 max)
{
   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j)
   { a[i][j] = ( (a[i][j]<min) ? min : ( (a[i][j]<max) ? a[i][j] : max ) ); }
}

void Matrix44_Constrain(real8 a[4][4], const real8 min, const real8 max)
{
   for (int i=0; (i<4); ++i)
   for (int j=0; (j<4); ++j)
   { a[i][j] = ( (a[i][j]<min) ? min : ( (a[i][j]<max) ? a[i][j] : max ) ); }
}

void MatrixMN_Constrain (real8  *a, const int m, const int n, const real8 min, const real8 max)
{
   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<(m*n)); ++i)
      { a[i] = ( (a[i]<min) ? min : ( (a[i]<max) ? a[i] : max ) ); }
   }
}

void MatrixMN_Constrain (real8 **a, const int m, const int n, const real8 min, const real8 max)
{
   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<m); ++i)
      for (int j=0; (j<n); ++j)
      { a[i][j] = ( (a[i][j]<min) ? min : ( (a[i][j]<max) ? a[i][j] : max ) ); }
   }
}

//---------------------------------------------------------------------------------------------------------

void Matrix22_Abs (real8   a[2][2])
{
   for (int i=0; (i<2); ++i)
   for (int j=0; (j<2); ++j)
   { a[i][j] = fabs(a[i][j]); }
}

void Matrix33_Abs (real8   a[3][3])
{
   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j)
   { a[i][j] = fabs(a[i][j]); }
}

void Matrix44_Abs (real8   a[4][4])
{
   for (int i=0; (i<4); ++i)
   for (int j=0; (j<4); ++j)
   { a[i][j] = fabs(a[i][j]); }
}

void MatrixMN_Abs (real8  *a, const int m, const int n)
{
   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<(m*n)); ++i)
      { a[i] = fabs(a[i]); }
   }
}

void MatrixMN_Abs (real8 **a, const int m, const int n)
{
   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<m); ++i)
      for (int j=0; (j<n); ++j)
      { a[i][j] = fabs(a[i][j]); }
   }
}

//---------------------------------------------------------------------------------------------------------

#define SIGN(a) ( ((a)<0.0) ? -1.0 : ((a)>0.0) ? 1.0 : 0.0 );

void Matrix22_Sign (real8   a[2][2])
{
   for (int i=0; (i<2); ++i)
   for (int j=0; (j<2); ++j)
   { a[i][j] = SIGN(a[i][j]); }
}

void Matrix33_Sign (real8   a[3][3])
{
   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j)
   { a[i][j] = SIGN(a[i][j]); }
}

void Matrix44_Sign (real8   a[4][4])
{
   for (int i=0; (i<4); ++i)
   for (int j=0; (j<4); ++j)
   { a[i][j] = SIGN(a[i][j]); }
}

void MatrixMN_Sign (real8  *a, const int m, const int n)
{
   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<(m*n)); ++i)
      { a[i] = SIGN(a[i]); }
   }
}

void MatrixMN_Sign (real8 **a, const int m, const int n)
{
   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<m); ++i)
      for (int j=0; (j<n); ++j)
      { a[i][j] = SIGN(a[i][j]); }
   }
}

#undef SIGN

//---------------------------------------------------------------------------------------------------------

void Matrix22_Residual (real8   a[2][2], real8   b[2][2] )
{
   for (int i=0; (i<2); ++i)
   for (int j=0; (j<2); ++j)
   { a[i][j] = fabs(a[i][j]-b[i][j]); }
}

void Matrix33_Residual (real8   a[3][3], real8   b[3][3] )
{
   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j)
   { a[i][j] = fabs(a[i][j]-b[i][j]); }
}

void Matrix44_Residual (real8   a[4][4], real8   b[4][4] )
{
   for (int i=0; (i<4); ++i)
   for (int j=0; (j<4); ++j)
   { a[i][j] = fabs(a[i][j]-b[i][j]); }
}

void MatrixMN_Residual (real8  *a, real8  *b, const int m, const int n)
{
   if ( a && b && (m>0) && (n>0) )
   {
      for (int i=0; (i<(m*n)); ++i)
      { a[i] = fabs(a[i]-b[i]); }
   }
}

void MatrixMN_Residual (real8 **a, real8 **b, const int m, const int n)
{
   if ( a && b && (m>0) && (n>0) )
   {
      for (int i=0; (i<m); ++i)
      for (int j=0; (j<n); ++j)
      { a[i][j] = fabs(a[i][j]-b[i][j]); }
   }
}

//---------------------------------------------------------------------------------------------------------

real8 Matrix22_RMS_Error (real8   a[2][2], real8   b[2][2] )
{
   real8 rms=0.0;

   for (int i=0; (i<2); ++i)
   for (int j=0; (j<2); ++j)
   { real8 s = fabs(a[i][j]-b[i][j]); rms += (s*s); }

   rms = sqrt(rms)/(2*2);

   return(rms);
}

real8 Matrix33_RMS_Error (real8   a[3][3], real8   b[3][3] )
{
   real8 rms=0.0;

   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j)
   { real8 s = fabs(a[i][j]-b[i][j]); rms += (s*s); }

   rms = sqrt(rms)/(3*3);

   return(rms);
}

real8 Matrix44_RMS_Error (real8   a[4][4], real8   b[4][4] )
{
   real8 rms=0.0;

   for (int i=0; (i<4); ++i)
   for (int j=0; (j<4); ++j)
   { real8 s = fabs(a[i][j]-b[i][j]); rms += (s*s); }

   rms = sqrt(rms)/(4*4);

   return(rms);
}

real8 MatrixMN_RMS_Error (real8  *a, real8  *b, const int m, const int n)
{
   real8 rms=0.0;

   if ( a && b && (m>0) && (n>0) )
   {
      for (int i=0; (i<(m*n)); ++i)
      { real8 s = fabs(a[i]-b[i]); rms += (s*s); }

      rms = sqrt(rms)/(m*n);
   }

   return(rms);
}

real8 MatrixMN_RMS_Error (real8 **a, real8 **b, const int m, const int n)
{
   real8 rms=0.0;

   if ( a && b && (m>0) && (n>0) )
   {
      for (int i=0; (i<m); ++i)
      for (int j=0; (j<n); ++j)
      { real8 s = fabs(a[i][j]-b[i][j]); rms += (s*s); }

      rms = sqrt(rms)/(m*n);
   }

   return(rms);
}

//---------------------------------------------------------------------------------------------------------

void Matrix22_Rand (real8   a[2][2], const real8 min, const real8 max)
{
   real8 s = ((real8)(max-min))/RAND_MAX;

   for (int i=0; (i<2); ++i)
   for (int j=0; (j<2); ++j)
   { a[i][j] = min+(s*rand()); }
}

void Matrix33_Rand (real8   a[3][3], const real8 min, const real8 max)
{
   real8 s = ((real8)(max-min))/RAND_MAX;

   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j)
   { a[i][j] = min+(s*rand()); }
}

void Matrix44_Rand (real8   a[4][4], const real8 min, const real8 max)
{
   real8 s = ((real8)(max-min))/RAND_MAX;

   for (int i=0; (i<4); ++i)
   for (int j=0; (j<4); ++j)
   { a[i][j] = min+(s*rand()); }
}

void MatrixMN_Rand (real8  *a, const int m, const int n, const real8 min, const real8 max)
{
   real8 s = ((real8)(max-min))/RAND_MAX;

   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<(m*n)); ++i)
      { a[i] = min+(s*rand()); }
   }
}

void MatrixMN_Rand (real8 **a, const int m, const int n, const real8 min, const real8 max)
{
   real8 s = ((real8)(max-min))/RAND_MAX;

   if ( a && (m>0) && (n>0) )
   {
      for (int i=0; (i<m); ++i)
      for (int j=0; (j<n); ++j)
      { a[i][j] = min+(s*rand()); }
   }
}


// ------------------------------------------------------------------------------------------
// 2x2 rotations...
//
// Matrix22_Rotate_Z (a,rz) : a = | cos(rz) -sin(rz) | (rz=rotation angle in radians)
//                                | sin(rz)  cos(rz) |
//
// Matrix33_Rotate_X (a,rx) : a = | 1   0        0     | (rx=rotation angle in radians)
//                                | 0 cos(rx) -sin(rx) |
//                                | 0 sin(rx)  cos(rx) |
//
// ------------------------------------------------------------------------------------------
// 3x3 rotations...
//
// Matrix33_Rotate_Y (a,ry) : a = |  cos(ry) 0 sin(ry) | (ry=rotation angle in radians)
//                                |    0     1   0     |
//                                | -sin(ry) 0 cos(ry) |
//
// Matrix33_Rotate_Z (a,rz) : a = | cos(rz) -sin(rz) 0 | (rz=rotation angle in radians)
//                                | sin(rz)  cos(rz) 0 |
//                                |   0        0     1 |
//
// ------------------------------------------------------------------------------------------
// 4x4 rotations...
//
// Matrix44_Rotate_X (a,rx) : a = | 1   0        0     0 | (rx=rotation angle in radians)
//                                | 0 cos(rx) -sin(rx) 0 |
//                                | 0 sin(rx)  cos(rx) 0 |
//                                | 0   0        0     1 |
//
// Matrix44_Rotate_Y (a,ry) : a = |  cos(ry) 0 sin(ry) 0 | (ry=rotation angle in radians)
//                                |    0     1   0     0 |
//                                | -sin(ry) 0 cos(ry) 0 |
//                                |    0     0   0     1 |
//
// Matrix44_Rotate_Z (a,rz) : a = | cos(rz) -sin(rz) 0 0 | (rz=rotation angle in radians)
//                                | sin(rz)  cos(rz) 0 0 |
//                                |   0        0     1 0 |
//                                |   0        0     0 1 |
//
//---------------------------------------------------------------------------------------------------------

void Matrix22_Rotate_Z(real8 a[2][2], const real8 rz)
{
   a[0][0] = cos(rz); a[0][1] = -sin(rz);
   a[1][0] = sin(rz); a[1][1] =  cos(rz);
}

void Matrix33_Rotate_X(real8 a[3][3], const real8 rx)
{
   a[0][0] = 1.0; a[0][1] =    0.0  ; a[0][2] =    0.0  ;
   a[1][0] = 0.0; a[1][1] =  cos(rx); a[1][2] = -sin(rx);
   a[2][0] = 0.0; a[2][1] =  sin(rx); a[2][2] =  cos(rx);
}

void Matrix33_Rotate_Y(real8 a[3][3], const real8 ry)
{
   a[0][0] =  cos(ry); a[0][1] = 0.0; a[0][2] = sin(ry);
   a[1][0] =    0.0  ; a[1][1] = 1.0; a[1][2] =   0.0  ;
   a[2][0] = -sin(ry); a[2][1] = 0.0; a[2][2] = cos(ry);
}

void Matrix33_Rotate_Z(real8 a[3][3], const real8 rz)
{
   a[0][0] = cos(rz); a[0][1] = -sin(rz); a[0][2] = 0.0;
   a[1][0] = sin(rz); a[1][1] =  cos(rz); a[1][2] = 0.0;
   a[2][0] =   0.0  ; a[2][1] =    0.0  ; a[2][2] = 1.0;
}

void Matrix33_Rotate (real8 a[3][3], const real8 rx, const real8 ry, const real8 rz)
{
   real8 mx[3][3];  Matrix33_Rotate_X(mx,rx);
   real8 my[3][3];  Matrix33_Rotate_Y(my,ry);
   real8 mz[3][3];  Matrix33_Rotate_Z(mz,rz);

   real8 mt[3][3];  Matrix33_Mul(mt,mx,my);
                    Matrix33_Mul(a ,mt,mz);
}

void Matrix44_Rotate_X(real8 a[4][4], const real8 rx)
{
   a[0][0] = 1.0; a[0][1] =    0.0  ; a[0][2] =    0.0  ; a[0][3] = 0.0;
   a[1][0] = 0.0; a[1][1] =  cos(rx); a[1][2] = -sin(rx); a[1][3] = 0.0;
   a[2][0] = 0.0; a[2][1] =  sin(rx); a[2][2] =  cos(rx); a[2][3] = 0.0;
   a[3][0] = 0.0; a[3][1] =    0.0  ; a[3][2] =    0.0  ; a[3][3] = 1.0;
}

void Matrix44_Rotate_Y(real8 a[4][4], const real8 ry)
{
   a[0][0] =  cos(ry); a[0][1] = 0.0; a[0][2] = sin(ry); a[0][3] = 0.0;
   a[1][0] =    0.0  ; a[1][1] = 1.0; a[1][2] =   0.0  ; a[1][3] = 0.0;
   a[2][0] = -sin(ry); a[2][1] = 0.0; a[2][2] = cos(ry); a[2][3] = 0.0;
   a[3][0] =    0.0  ; a[3][1] = 0.0; a[3][2] =   0.0  ; a[3][3] = 1.0;
}

void Matrix44_Rotate_Z(real8 a[4][4], const real8 rz)
{
   a[0][0] = cos(rz); a[0][1] = -sin(rz); a[0][2] = 0.0; a[0][3] = 0.0;
   a[1][0] = sin(rz); a[1][1] =  cos(rz); a[1][2] = 0.0; a[1][3] = 0.0;
   a[2][0] =   0.0  ; a[2][1] =    0.0  ; a[2][2] = 1.0; a[2][3] = 0.0;
   a[3][0] =   0.0  ; a[3][1] =    0.0  ; a[3][2] = 0.0; a[3][3] = 1.0;
}

void Matrix44_Rotate (real8 a[4][4], const real8 rx, const real8 ry, const real8 rz)
{
   real8 mx[4][4];  Matrix44_Rotate_X(mx,rx);
   real8 my[4][4];  Matrix44_Rotate_Y(my,ry);
   real8 mz[4][4];  Matrix44_Rotate_Z(mz,rz);

   real8 mt[4][4];  Matrix44_Mul(mt,mx,my);
                    Matrix44_Mul(a ,mt,mz);
}

// Matrix22_Diag  (a,s)      : a = | s 0 |
//                                 | 0 s |
//
// Matrix22_Diag  (a,x,y)    : a = | x 0 |
//                                 | 0 y |
//
// Matrix33_Diag  (a,s)      : a = | s 0 0 |
//                                 | 0 s 0 |
//                                 | 0 0 s |
//
// Matrix33_Diag  (a,x,y,z)  : a = | x 0 0 |
//                                 | 0 y 0 |
//                                 | 0 0 z |
//
// Matrix44_Diag  (a,s)      : a = | s 0 0 0 |
//                                 | 0 s 0 0 |
//                                 | 0 0 s 0 |
//                                 | 0 0 0 1 |
//
// Matrix44_Diag  (a,x,y,z)  : a = | x 0 0 0 |
//                                 | 0 y 0 0 |
//                                 | 0 0 z 0 |
//                                 | 0 0 0 1 |
//
//---------------------------------------------------------------------------------------------------------

void Matrix22_Diag(real8 a[2][2], const real8 s)
{
   a[0][0]=s; a[0][1]=0;
   a[1][0]=0; a[1][1]=s;
}

void Matrix22_Diag(real8 a[2][2], const real8 sx, const real8 sy)
{
   a[0][0]=sx; a[0][1]=0;
   a[1][0]=0;  a[1][1]=sy;
}

void Matrix33_Diag(real8 a[3][3], const real8 s)
{
   a[0][0]=s; a[0][1]=0; a[0][2]=0;
   a[1][0]=0; a[1][1]=s; a[1][2]=0;
   a[2][0]=0; a[2][1]=0; a[2][2]=s;
}

void Matrix33_Diag(real8 a[3][3], const real8 sx, const real8 sy, const real8 sz)
{
   a[0][0]=sx; a[0][1]=0;  a[0][2]=0;
   a[1][0]=0;  a[1][1]=sy; a[1][2]=0;
   a[2][0]=0;  a[2][1]=0;  a[2][2]=sz;
}

void Matrix44_Diag(real8 a[4][4], const real8 s)
{
   a[0][0]=s; a[0][1]=0; a[0][2]=0; a[0][3]=0;
   a[1][0]=0; a[1][1]=s; a[1][2]=0; a[1][3]=0;
   a[2][0]=0; a[2][1]=0; a[2][2]=s; a[2][3]=0;
   a[3][0]=0; a[3][1]=0; a[3][2]=0; a[3][3]=1;
}

void Matrix44_Diag(real8 a[4][4], const real8 sx, const real8 sy, const real8 sz)
{
   a[0][0]=sx; a[0][1]=0;  a[0][2]=0;  a[0][3]=0;
   a[1][0]=0;  a[1][1]=sy; a[1][2]=0;  a[1][3]=0;
   a[2][0]=0;  a[2][1]=0;  a[2][2]=sz; a[2][3]=0;
   a[3][0]=0;  a[3][1]=0;  a[3][2]=0;  a[3][3]=1;
}

void MatrixMN_Diag(real8  *a, const int m, const int n, const real8 *s)
{
  if (a && s && (m>0) && (n>0))
  {
     for (int i=0,k=0; (i<m); ++i)
     for (int j=0;     (j<n); ++j, k++)
        a[k] = ( (i==j) ? s[i] : 0.0 );
  }
}

void MatrixMN_Diag(real8 **a, const int m, const int n, const real8 *s)
{
  if (a && s && (m>0) && (n>0))
  {
     for (int i=0; (i<m); ++i)
     for (int j=0; (j<n); ++j)
        a[i][j] = ( (i==j) ? s[i] : 0.0 );
  }
}

// arithmetic ops....
//---------------------------------------------------------------------------------------------------------

void Matrix22_Add (real8 a[2][2], real8 b[2][2], real8 c[2][2])
{
   a[0][0] = b[0][0]+c[0][0];
   a[0][1] = b[0][1]+c[0][1];
   a[1][0] = b[1][0]+c[1][0];
   a[1][1] = b[1][1]+c[1][1];
}

void Matrix22_Sub (real8 a[2][2], real8 b[2][2], real8 c[2][2])
{
   a[0][0] = b[0][0]-c[0][0];
   a[0][1] = b[0][1]-c[0][1];
   a[1][0] = b[1][0]-c[1][0];
   a[1][1] = b[1][1]-c[1][1];
}

void Matrix22_Mul (real8 a[2][2], real8 b[2][2], real8 c[2][2])
{
   a[0][0] = b[0][0]*c[0][0] + b[0][1]*c[1][0];
   a[0][1] = b[0][0]*c[0][1] + b[0][1]*c[1][1];
   a[1][0] = b[1][0]*c[0][0] + b[1][1]*c[1][0];
   a[1][1] = b[1][0]*c[0][1] + b[1][1]*c[1][1];
}

void Matrix22_Mul (real8 a[2][2], real8 b[2][2], real8 c[2][2], real8 d[2][2])
{
   real8 t[2][2];

   Matrix22_Mul(t,b,c);  // T = B*C
   Matrix22_Mul(a,t,d);  // A = B*C*D
}

void Matrix22_Mul (real8 a[2][2], real8 b[2][3], real8 c[3][2] )
{
    a[0][0] =  b[0][0]*c[0][0] + b[0][1]*c[1][0] + b[0][2]*c[2][0];
    a[0][1] =  b[0][0]*c[0][1] + b[0][1]*c[1][1] + b[0][2]*c[2][1];
    a[1][0] =  b[1][0]*c[0][0] + b[1][1]*c[1][0] + b[1][2]*c[2][0];
    a[1][1] =  b[1][0]*c[0][1] + b[1][1]*c[1][1] + b[1][2]*c[2][1];
}

void Matrix33_Add (real8 a[3][3], real8 b[3][3], real8 c[3][3])
{
   for (int i=0; (i<3); ++i)
   {
      a[i][0]=b[i][0]+c[i][0];
      a[i][1]=b[i][1]+c[i][1];
      a[i][2]=b[i][2]+c[i][2];
   }
}

void Matrix33_Sub (real8 a[3][3], real8 b[3][3], real8 c[3][3])
{
   for (int i=0; (i<3); ++i)
   {
      a[i][0]=b[i][0]-c[i][0];
      a[i][1]=b[i][1]-c[i][1];
      a[i][2]=b[i][2]-c[i][2];
   }
}

void Matrix33_Mul (real8 a[3][3], real8 b[3][3], real8 c[3][3])
{
   for (int i=0; (i<3); ++i)
   {
      a[i][0] = b[i][0]*c[0][0] + b[i][1]*c[1][0] + b[i][2]*c[2][0];
      a[i][1] = b[i][0]*c[0][1] + b[i][1]*c[1][1] + b[i][2]*c[2][1];
      a[i][2] = b[i][0]*c[0][2] + b[i][1]*c[1][2] + b[i][2]*c[2][2];
   }
}

void Matrix33_Mul (real8 a[3][3], real8 b[3][3], real8 c[3][3], real8 d[3][3])
{
   real8 t[3][3];

   Matrix33_Mul(t,b,c);  // T = B*C
   Matrix33_Mul(a,t,d);  // A = B*C*D
}

void Matrix44_Add (real8 a[4][4], real8 b[4][4], real8 c[4][4])
{
   for (int i=0; (i<4); ++i)
   {
      a[i][0]=b[i][0]+c[i][0];
      a[i][1]=b[i][1]+c[i][1];
      a[i][2]=b[i][2]+c[i][2];
      a[i][3]=b[i][3]+c[i][3];
   }
}

void Matrix44_Sub (real8 a[4][4], real8 b[4][4], real8 c[4][4])
{
   for (int i=0; (i<4); ++i)
   {
      a[i][0]=b[i][0]-c[i][0];
      a[i][1]=b[i][1]-c[i][1];
      a[i][2]=b[i][2]-c[i][2];
      a[i][3]=b[i][3]-c[i][3];
   }
}

void Matrix44_Mul (real8 a[4][4], real8 b[4][4], real8 c[4][4])
{
   for (int i=0; (i<4); ++i)
   {
      a[i][0] = b[i][0]*c[0][0] + b[i][1]*c[1][0] + b[i][2]*c[2][0] + b[i][3]*c[3][0];
      a[i][1] = b[i][0]*c[0][1] + b[i][1]*c[1][1] + b[i][2]*c[2][1] + b[i][3]*c[3][1];
      a[i][2] = b[i][0]*c[0][2] + b[i][1]*c[1][2] + b[i][2]*c[2][2] + b[i][3]*c[3][2];
      a[i][3] = b[i][0]*c[0][3] + b[i][1]*c[1][3] + b[i][2]*c[2][3] + b[i][3]*c[3][3];
   }
}

void Matrix44_Mul (real8 a[4][4], real8 b[4][4], real8 c[4][4], real8 d[4][4])
{
   real8 t[4][4];

   Matrix44_Mul(t,b,c);  // T = B*C
   Matrix44_Mul(a,t,d);  // A = B*C*D
}

void Matrix66_Add (real8 a[6][6], real8 b[6][6], real8 c[6][6])
{
   for (int i=0; (i<6); ++i)
   {
      a[i][0]=b[i][0]+c[i][0];
      a[i][1]=b[i][1]+c[i][1];
      a[i][2]=b[i][2]+c[i][2];
      a[i][3]=b[i][3]+c[i][3];
      a[i][4]=b[i][4]+c[i][4];
      a[i][5]=b[i][5]+c[i][5];
   }
}

void Matrix66_Sub (real8 a[6][6], real8 b[6][6], real8 c[6][6])
{
   for (int i=0; (i<6); ++i)
   {
      a[i][0]=b[i][0]-c[i][0];
      a[i][1]=b[i][1]-c[i][1];
      a[i][2]=b[i][2]-c[i][2];
      a[i][3]=b[i][3]-c[i][3];
      a[i][4]=b[i][4]-c[i][4];
      a[i][5]=b[i][5]-c[i][5];
   }
}

void Matrix66_Mul (real8 a[6][6], real8 b[6][6], real8 c[6][6])
{
   for (int i=0; (i<6); ++i)
   {
      a[i][0] =   b[i][0]*c[0][0] + b[i][1]*c[1][0] + b[i][2]*c[2][0] + b[i][3]*c[3][0] + b[i][4]*c[4][0] + b[i][5]*c[5][0];
      a[i][1] =   b[i][0]*c[0][1] + b[i][1]*c[1][1] + b[i][2]*c[2][1] + b[i][3]*c[3][1] + b[i][4]*c[4][1] + b[i][5]*c[5][1];
      a[i][2] =   b[i][0]*c[0][2] + b[i][1]*c[1][2] + b[i][2]*c[2][2] + b[i][3]*c[3][2] + b[i][4]*c[4][2] + b[i][5]*c[5][2];
      a[i][3] =   b[i][0]*c[0][3] + b[i][1]*c[1][3] + b[i][2]*c[2][3] + b[i][3]*c[3][3] + b[i][4]*c[4][3] + b[i][5]*c[5][3];
      a[i][4] =   b[i][0]*c[0][4] + b[i][1]*c[1][4] + b[i][2]*c[2][4] + b[i][3]*c[3][4] + b[i][4]*c[4][4] + b[i][5]*c[5][4];
      a[i][5] =   b[i][0]*c[0][5] + b[i][1]*c[1][5] + b[i][2]*c[2][5] + b[i][3]*c[3][5] + b[i][4]*c[4][5] + b[i][5]*c[5][5];
   }
}

// Matrix-vector multiplies...
//---------------------------------------------------------------------------------------------------------

void Matrix22_Vmul (real8 y[2], real8 a[2][2], real8 x[2]   )
{
   y[0] = a[0][0]*x[0] + a[0][1]*x[1];
   y[1] = a[1][0]*x[0] + a[1][1]*x[1];
}

void Matrix22_Vmul (real8 y[2], real8 x[2]   , real8 a[2][2])
{
   y[0] = a[0][0]*x[0] + a[1][0]*x[1];
   y[1] = a[0][1]*x[0] + a[1][1]*x[1];
}

void Matrix23_Vmul(real8 y[2], real8 a[2][3], real8 x[3] )
{
   y[0] = a[0][0]*x[0] + a[0][1]*x[1] + a[0][2]*x[2];
   y[1] = a[1][0]*x[0] + a[1][1]*x[1] + a[1][2]*x[2];
}

void Matrix32_Vmul(real8 y[3], real8 a[3][2], real8 x[2])
{
    y[0] = a[0][0]*x[0] + a[0][1]*x[1];
    y[1] = a[1][0]*x[0] + a[1][1]*x[1];
    y[2] = a[2][0]*x[0] + a[2][1]*x[1];
}

void Matrix33_Vmul (real8 y[3], real8 a[3][3], real8 x[3]   )
{
   y[0] = a[0][0]*x[0] + a[0][1]*x[1] + a[0][2]*x[2];
   y[1] = a[1][0]*x[0] + a[1][1]*x[1] + a[1][2]*x[2];
   y[2] = a[2][0]*x[0] + a[2][1]*x[1] + a[2][2]*x[2];
}

void Matrix33_Vmul (real8 y[3], real8 x[3]   , real8 a[3][3])
{
   y[0] = a[0][0]*x[0] + a[1][0]*x[1] + a[2][0]*x[2];
   y[1] = a[0][1]*x[0] + a[1][1]*x[1] + a[2][1]*x[2];
   y[2] = a[0][2]*x[0] + a[1][2]*x[1] + a[2][2]*x[2];
}

void Matrix33_Vmul (real8 a[3][3], real8 v[3])
{
   a[0][0] = v[0]*v[0];
   a[0][1] = v[0]*v[1];
   a[0][2] = v[0]*v[2];

   a[1][0] = v[1]*v[0];
   a[1][1] = v[1]*v[1];
   a[1][2] = v[1]*v[2];

   a[2][0] = v[2]*v[0];
   a[2][1] = v[2]*v[1];
   a[2][2] = v[2]*v[2];
}

void Matrix33_Vmul (real8 a[3][3], real8 v[3], real8 vt[3] )
{
   a[0][0] = v[0]*vt[0];
   a[0][1] = v[0]*vt[1];
   a[0][2] = v[0]*vt[2];

   a[1][0] = v[1]*vt[0];
   a[1][1] = v[1]*vt[1];
   a[1][2] = v[1]*vt[2];

   a[2][0] = v[2]*vt[0];
   a[2][1] = v[2]*vt[1];
   a[2][2] = v[2]*vt[2];
}

void Matrix33_VVt_Mul (real8 a[3][3], real8 v[3])
{
   a[0][0] = v[0]*v[0];
   a[0][1] = v[0]*v[1];
   a[0][2] = v[0]*v[2];

   a[1][0] = v[1]*v[0];
   a[1][1] = v[1]*v[1];
   a[1][2] = v[1]*v[2];

   a[2][0] = v[2]*v[0];
   a[2][1] = v[2]*v[1];
   a[2][2] = v[2]*v[2];
}

void Matrix33_VVt_Mul (real8 a[3][3], real8 v[3], real8 vt[3])
{
   a[0][0] = v[0]*vt[0];
   a[0][1] = v[0]*vt[1];
   a[0][2] = v[0]*vt[2];

   a[1][0] = v[1]*vt[0];
   a[1][1] = v[1]*vt[1];
   a[1][2] = v[1]*vt[2];

   a[2][0] = v[2]*vt[0];
   a[2][1] = v[2]*vt[1];
   a[2][2] = v[2]*vt[2];
}

void Matrix44_Vmul (real8 y[3], real8 a[4][4], real8 x[3]   )
{
   y[0] = a[0][0]*x[0] + a[0][1]*x[1] + a[0][2]*x[2] + a[0][3];
   y[1] = a[1][0]*x[0] + a[1][1]*x[1] + a[1][2]*x[2] + a[1][3];
   y[2] = a[2][0]*x[0] + a[2][1]*x[1] + a[2][2]*x[2] + a[2][3];
}

void Matrix44_Vmul (real8 y[3], real8 x[3]   , real8 a[4][4])
{
   y[0] = a[0][0]*x[0] + a[1][0]*x[1] + a[2][0]*x[2] + a[3][0];
   y[1] = a[0][1]*x[0] + a[1][1]*x[1] + a[2][1]*x[2] + a[3][1];
   y[2] = a[0][2]*x[0] + a[1][2]*x[1] + a[2][2]*x[2] + a[3][2];
}

void Matrix43_Vmul (real8 y[4], real8 a[4][3], real8 x[3] )
{
   y[0] = a[0][0]*x[0] + a[0][1]*x[1] + a[0][2]*x[2];
   y[1] = a[1][0]*x[0] + a[1][1]*x[1] + a[1][2]*x[2];
   y[2] = a[2][0]*x[0] + a[2][1]*x[1] + a[2][2]*x[2];
   y[3] = a[3][0]*x[0] + a[3][1]*x[1] + a[3][2]*x[2];
}

void Matrix63_Vmul (real8 y[6], real8 a[6][3], real8 x[3] )
{
   y[0] = a[0][0]*x[0] + a[0][1]*x[1] + a[0][2]*x[2];
   y[1] = a[1][0]*x[0] + a[1][1]*x[1] + a[1][2]*x[2];
   y[2] = a[2][0]*x[0] + a[2][1]*x[1] + a[2][2]*x[2];
   y[3] = a[3][0]*x[0] + a[3][1]*x[1] + a[3][2]*x[2];
   y[4] = a[4][0]*x[0] + a[4][1]*x[1] + a[4][2]*x[2];
   y[5] = a[5][0]*x[0] + a[5][1]*x[1] + a[5][2]*x[2];
}

// MatrixMN_Add : A=B+C
//---------------------------------------------------------------------------------------------------------

void MatrixMN_Add
(
         real8  *a,   ///< matrix to receive the result  (dimension = MxN)
   const real8  *b,   ///< left -hand matrix of multiply (dimension = MxN)
   const real8  *c,   ///< right-hand matrix of multiply (dimension = MxN)
   const int     m,   ///< matrix dimension (M)
   const int     n    ///< matrix dimension (N)
)
{
  if (a && b && c && (m>0) && (n>0))
  {
     for (int i=0; (i<m*n); ++i)
        a[i]=b[i]+c[i];
  }
}

void MatrixMN_Add
(
         real8 **a,   ///< matrix to receive the result  (dimension = MxN)
         real8 **b,   ///< left -hand matrix of multiply (dimension = MxN)
         real8 **c,   ///< right-hand matrix of multiply (dimension = MxN)
   const int     m,   ///< matrix dimension (M)
   const int     n    ///< matrix dimension (N)
)
{
  if (a && b && c && (m>0) && (n>0))
  {
     for (int i=0; (i<m); ++i)
     for (int j=0; (j<n); ++j)
        a[i][j]=b[i][j]+c[i][j];
  }
}

// MatrixMN_Add : A=B-C
//---------------------------------------------------------------------------------------------------------

void MatrixMN_Sub
(
         real8  *a,   ///< matrix to receive the result  (dimension = MxN)
   const real8  *b,   ///< left -hand matrix of multiply (dimension = MxN)
   const real8  *c,   ///< right-hand matrix of multiply (dimension = MxN)
   const int     m,   ///< matrix dimension (M)
   const int     n    ///< matrix dimension (N)
)
{
  if (a && b && c && (m>0) && (n>0))
  {
     for (int i=0; (i<m*n); ++i)
        a[i]=b[i]-c[i];
  }
}

void MatrixMN_Sub
(
         real8 **a,   ///< matrix to receive the result  (dimension = MxN)
         real8 **b,   ///< left -hand matrix of multiply (dimension = MxN)
         real8 **c,   ///< right-hand matrix of multiply (dimension = MxN)
   const int     m,   ///< matrix dimension (M)
   const int     n    ///< matrix dimension (N)
)
{
  if (a && b && c && (m>0) && (n>0))
  {
     for (int i=0; (i<m); ++i)
     for (int j=0; (j<n); ++j)
        a[i][j]=b[i][j]-c[i][j];
  }
}

// MatrixMN_Mul() : matrix multiply  A=B*C
//
// Performs a general matrix-matrix multiply.
//
// Several assumptions...
//   - all matrices are dense, zero-indexed, row-major  memory spaces.
//   - the destination matrix (A) is not pointing to the same memory spaces
//     as either B or C (B and C can point to the same space).
//   - the dimensions of the matrices are compatible for the multiply
//       A = MxO matrix (result)
//       B = MxN matrix (source)
//       C = NxO matrix (source)
//---------------------------------------------------------------------------------------------------------

static real8 rowcol
(
   const real8  *r,   ///< points to matrix row
   const real8  *c,   ///< points to matrix column
   const int     m,   ///< number of elements to multiply
   const int     n    ///< stride of elements in column
)
{
   real8 s=0.0;

   for (int i=0,j=0; (i<m); i++, j+=n)
      s += r[i]*c[j];

   return(s);
}

void MatrixMN_Mul
(
         real8  *a,   ///< matrix to receive the result  (dimension = MxP)
   const real8  *b,   ///< left -hand matrix of multiply (dimension = MxN)
   const real8  *c,   ///< right-hand matrix of multiply (dimension = NxP)
   const int     m,   ///< matrix dimension (M)
   const int     n,   ///< matrix dimension (N)
   const int     p    ///< matrix dimension (P)
)
{
  if (a && b && c && (a!=b) && (a!=c) && (m>0) && (n>0) && (p>0) )
  {
     for (int i=0,k=0; (i<m); ++i)
     for (int j=0    ; (j<p); ++j, ++k)
     {
        a[k] = rowcol(b+i*n,c+j,n,p);
     }
  }
}

void MatrixMN_Mul
(
         real8 **a,   ///< matrix to receive the result  (dimension = MxP)
         real8 **b,   ///< left -hand matrix of multiply (dimension = MxN)
         real8 **c,   ///< right-hand matrix of multiply (dimension = NxP)
   const int     m,   ///< matrix dimension (M)
   const int     n,   ///< matrix dimension (N)
   const int     p    ///< matrix dimension (P)
)
{
  if (a && b && c && (a!=b) && (a!=c) && (m>0) && (n>0) && (p>0) )
  {
     for (int i=0; (i<m); ++i)
     for (int j=0; (j<p); ++j)
     {
        a[i][j] = 0.0;

        for (int k=0; (k<n); ++k)
           a[i][j] += b[i][k]*c[k][j];
     }
  }
}

// Matrix22_Transpose()
// Matrix33_Transpose()
// Matrix44_Transpose()
// MatrixMN_Transpose()
//    A(NxM) = Transpose(B(MxN))  (A and B cannot be the same matrix)
//    A(NxM) = Transpose(A(MxN))  (in-place methods provided)
//---------------------------------------------------------------------------------------------------------

#define SWAP(a,b) { real8 s=(a); (a)=(b); (b)=s; }

void Matrix22_Transpose (real8 a[2][2])
{
   SWAP( a[0][1], a[1][0] );
}

void Matrix33_Transpose (real8 a[3][3])
{
   SWAP( a[0][1], a[1][0] );
   SWAP( a[0][2], a[2][0] );

   SWAP( a[2][1], a[1][2] );
}

void Matrix44_Transpose (real8 a[4][4])
{
   SWAP( a[0][1], a[1][0] );
   SWAP( a[0][2], a[2][0] );
   SWAP( a[0][3], a[3][0] );

   SWAP( a[1][2], a[2][1] );
   SWAP( a[1][3], a[3][1] );

   SWAP( a[2][3], a[3][2] );
}

void Matrix66_Transpose (real8 a[6][6])
{
   SWAP( a[0][1], a[1][0] );
   SWAP( a[0][2], a[2][0] );
   SWAP( a[0][3], a[3][0] );
   SWAP( a[0][4], a[4][0] );
   SWAP( a[0][5], a[5][0] );

   SWAP( a[1][2], a[2][1] );
   SWAP( a[1][3], a[3][1] );
   SWAP( a[1][4], a[4][1] );
   SWAP( a[1][5], a[5][1] );

   SWAP( a[2][3], a[3][2] );
   SWAP( a[2][4], a[4][2] );
   SWAP( a[2][5], a[5][2] );

   SWAP( a[3][4], a[4][3] );
   SWAP( a[3][5], a[5][3] );

   SWAP( a[4][5], a[5][4] );
}

void Matrix22_Transpose (real8 a[2][2], real8 b[2][2])
{
   a[0][0] = b[0][0];
   a[0][1] = b[1][0];
   a[1][0] = b[0][1];
   a[1][1] = b[1][1];
}

void Matrix33_Transpose (real8 a[3][3], real8 b[3][3])
{
   for (int i=0; (i<3); ++i)
   {
      a[i][0] = b[0][i];
      a[i][1] = b[1][i];
      a[i][2] = b[2][i];
   }
}

void Matrix44_Transpose (real8 a[4][4], real8 b[4][4])
{
   for (int i=0; (i<4); ++i)
   {
      a[i][0] = b[0][i];
      a[i][1] = b[1][i];
      a[i][2] = b[2][i];
      a[i][3] = b[3][i];
   }
}

void Matrix66_Transpose (real8 a[6][6], real8 b[6][6])
{
   for (int i=0; (i<6); ++i)
   {
      a[i][0] = b[0][i];
      a[i][1] = b[1][i];
      a[i][2] = b[2][i];
      a[i][3] = b[3][i];
      a[i][4] = b[4][i];
      a[i][5] = b[5][i];
   }
}

// Matrix32_Transpose() : A(3x2) = B(2x3)
//---------------------------------------------------------------------------------------------------------

void Matrix23_Transpose (real8 a[3][2], real8 b[2][3] )
{
   a[0][0] = b[0][0];
   a[0][1] = b[1][0];
   a[1][0] = b[0][1];
   a[1][1] = b[1][1];
   a[2][0] = b[0][2];
   a[2][1] = b[1][2];
}

// Matrix32_Transpose() : A(2x3) = B(3x2)
//---------------------------------------------------------------------------------------------------------

void Matrix32_Transpose (real8 a[2][3], real8 b[3][2] )
{
   a[0][0] = b[0][0];
   a[0][1] = b[1][0];
   a[0][2] = b[2][0];
   a[1][0] = b[0][1];
   a[1][1] = b[1][1];
   a[1][2] = b[2][1];
}

void MatrixMN_Transpose (real8 *a, const int m, const int n)
{
   if (a && (m>0) && (n>0) )
   {
      real8 *b = new real8[m*n];

      if (b)
      {
         for (int i=0; (i<(m*n)); ++i) { b[i]=a[i]; }

         for (int i=0; (i<m); ++i)
         for (int j=0; (j<n); ++j) { a[j*m+i]=b[i*n+j]; }

         delete [] b;
      }
   }
}

void MatrixMN_Transpose (real8 *a, real8 *b, const int m, const int n)
{
   if (a && b && (m>0) && (n>0) )
   {
      for (int i=0; (i<m); ++i)
      for (int j=0; (j<n); ++j) { a[j*m+i]=b[i*n+j]; }
   }
}

void MatrixMN_Transpose (real8 **a, real8 **b, const int m, const int n)
{
   if (a && b && (a!=b) && (m>0) && (n>0) )
   {
      for (int i=0; (i<m); ++i)
      for (int j=0; (j<n); ++j) { a[j][i]=b[i][j]; }
   }
}

//---------------------------------------------------------------------------------------------------------
// comparison...
//---------------------------------------------------------------------------------------------------------

// nearly_eq()
//
// Local implementation of relative comparison between two floating point numbers.
// Note that this implementation deals with some of the comparison edge caese where
// one or both of the arguments is near zero, or there is a large difference in scale
// between the two arguments.
//
// Note that this routine is invoked by all of the matrix nearness tests.
//---------------------------------------------------------------------------------------------------------

static int nearly_eq(const real8 a, const real8 b, const real8 tol)
{

   if (a==b) return(1);                     // (bit-for-bit equal)

   // if either arg is near zero, can't use relative error...

   if (   (fabs(a)<1.0e-12)
       || (fabs(b)<1.0e-12) )
      return(  (fabs(a-b)<tol) ? 1 : 0 );

   // otherwise compute relative error by scaling by largest arg...

   real8 rerr = ( (fabs(a)>fabs(b)) ?  fabs((a-b)/a) : fabs((a-b)/b) );

   return( (rerr<tol) ? 1 : 0 );
}

// Matrix22_Near()
// Matrix33_Near()
// Matrix44_Near()
// Matrix66_Near()
//
// Returns 1 (true) or 0 (false) if ALL of the elements of A are less than tolerance.
// Overloaded to allow comparisons with matrices and scalars.
//---------------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------------

int Matrix22_Near (const real8 a[2][2], const real8 s      , const real8 tol)
{
   for (int i=0; (i<2); ++i)
   for (int j=0; (j<2); ++j)
   {
      if ( !nearly_eq(a[i][j],s,tol) ) return(0);
   }

   return(1);
}

int Matrix22_Near (const real8 a[2][2], const real8 b[2][2], const real8 tol)
{
   for (int i=0; (i<2); ++i)
   for (int j=0; (j<2); ++j)
   {
      if ( !nearly_eq(a[i][j],b[i][j],tol) ) return(0);
   }

   return(1);
}

int Matrix33_Near (const real8 a[3][3], const real8 s      , const real8 tol)
{
   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j)
   {
      if ( !nearly_eq(a[i][j],s,tol) ) return(0);
   }

   return(1);
}

int Matrix33_Near (const real8 a[3][3], const real8 b[3][3], const real8 tol)
{
   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j)
   {
      if ( !nearly_eq(a[i][j],b[i][j],tol) ) return(0);
   }

   return(1);
}

int Matrix44_Near (const real8 a[4][4], const real8 s      , const real8 tol)
{
   for (int i=0; (i<4); ++i)
   for (int j=0; (j<4); ++j)
   {
      if ( !nearly_eq(a[i][j],s,tol) ) return(0);
   }

   return(1);
}

int Matrix44_Near (const real8 a[4][4], const real8 b[4][4], const real8 tol)
{
   for (int i=0; (i<4); ++i)
   for (int j=0; (j<4); ++j)
   {
      if ( !nearly_eq(a[i][j],b[i][j],tol) ) return(0);
   }

   return(1);
}

int Matrix66_Near (const real8 a[6][6], const real8 s      , const real8 tol)
{
   for (int i=0; (i<6); ++i)
   for (int j=0; (j<6); ++j)
   {
      if ( !nearly_eq(a[i][j],s,tol) ) return(0);
   }

   return(1);
}

int Matrix66_Near (const real8 a[6][6], const real8 b[6][6], const real8 tol)
{
   for (int i=0; (i<6); ++i)
   for (int j=0; (j<6); ++j)
   {
      if ( !nearly_eq(a[i][j],b[i][j],tol) ) return(0);
   }

   return(1);
}

//---------------------------------------------------------------------------------------------------------

real8 Matrix22_Det (real8 a[2][2])
{
   return( a[0][0]*a[1][1] - a[0][1]*a[1][0] );
}

//---------------------------------------------------------------------------------------------------------

real8 Matrix22_Norm (real8 a[2][2])
{
   real8 norm=0.0;

   real8 sum  = fabs(a[0][0]);
         sum += fabs(a[0][1]);  norm = ( (sum>norm) ? sum : norm );

         sum  = fabs(a[1][0]);
         sum += fabs(a[1][1]);  norm = ( (sum>norm) ? sum : norm );

   return (norm);
}

//---------------------------------------------------------------------------------------------------------

real8 Matrix22_Condition (real8 a[2][2])
{
   real8 tmp[2][2];

   Matrix22_Inverse(tmp,a);

   return (   Matrix22_Norm(tmp)
            * Matrix22_Norm(a  ) );
}

//---------------------------------------------------------------------------------------------------------

void Matrix22_Adjoint (real8 a[2][2], real8 b[2][2])
{
   a[0][0] =  b[1][1];
   a[0][1] = -b[0][1];
   a[1][0] = -b[1][0];
   a[1][1] =  b[0][0];
}

//---------------------------------------------------------------------------------------------------------

int Matrix22_Inverse (real8 a[2][2], real8 b[2][2])
{
   real8 det = (b[0][0]*b[1][1]) - (b[0][1]*b[1][0]);  //  det = det(A)

   if ( fabs(det) < 1.e-20 )
   {
      printf("%s::%s(%d) - noninvertible matrix encountered det=%le\n", __FILE__, __func__, __LINE__, det );
      return(-1);
   }

   det = ( (fabs(det)>0.0) ? 1.0/det : 0.0 );    // det = 1/det(A)

   a[0][0] =  b[1][1]*det;                       // A = adj(B)/det(B)
   a[0][1] = -b[0][1]*det;                       //
   a[1][0] = -b[1][0]*det;                       //
   a[1][1] =  b[0][0]*det;                       //

   return(1);
}


// Matrix33_Det() : returns determinant of A
//---------------------------------------------------------------------------------------------------------

real8 Matrix33_Det (real8 a[3][3])
{
   return(   (a[0][0] * (a[1][1]*a[2][2] - a[2][1]*a[1][2]))
           - (a[0][1] * (a[1][0]*a[2][2] - a[2][0]*a[1][2]))
           + (a[0][2] * (a[1][0]*a[2][1] - a[2][0]*a[1][1])) );
}

// Matrix33_Norm() : returns norm of matrix :   norm = ||A||
//---------------------------------------------------------------------------------------------------------

real8 Matrix33_Norm (real8 a[3][3])
{
   real8 norm=0.0;

   real8 sum  = fabs(a[0][0]);
         sum += fabs(a[0][1]);
         sum += fabs(a[0][2]);  norm = ( (sum>norm) ? sum : norm );

         sum  = fabs(a[1][0]);
         sum += fabs(a[1][1]);
         sum += fabs(a[1][2]);  norm = ( (sum>norm) ? sum : norm );

         sum  = fabs(a[2][0]);
         sum += fabs(a[2][1]);
         sum += fabs(a[2][2]);  norm = ( (sum>norm) ? sum : norm );

   return (norm);
}

// Matrix33_Condition() : returns condition number of matrix :   cond = ||A|| * ||Inv(A)||
//---------------------------------------------------------------------------------------------------------

real8 Matrix33_Condition (real8 a[3][3])
{
   real8 tmp[3][3];

   Matrix33_Inverse(tmp,a);

   return (   Matrix33_Norm(tmp)
            * Matrix33_Norm(a  ) );
}

// Matrix33_Adjoint() : A = Adj(B)    (A and B cannot be the same matrices)
//---------------------------------------------------------------------------------------------------------

void Matrix33_Adjoint(real8 a[3][3], real8 b[3][3])
{
    a[0][0] = b[1][1]*b[2][2] - b[1][2]*b[2][1];
    a[0][1] = b[0][2]*b[2][1] - b[0][1]*b[2][2];
    a[0][2] = b[0][1]*b[1][2] - b[0][2]*b[1][1];
    a[1][0] = b[1][2]*b[2][0] - b[1][0]*b[2][2];
    a[1][1] = b[0][0]*b[2][2] - b[0][2]*b[2][0];
    a[1][2] = b[0][2]*b[1][0] - b[0][0]*b[1][2];
    a[2][0] = b[1][0]*b[2][1] - b[1][1]*b[2][0];
    a[2][1] = b[0][1]*b[2][0] - b[0][0]*b[2][1];
    a[2][2] = b[0][0]*b[1][1] - b[0][1]*b[1][0];
}

// Matrix33_Inverse() : A = Inv(B)
//
// Determines inverse using Kramer's rule
//---------------------------------------------------------------------------------------------------------

int Matrix33_Inverse (real8 a[3][3], real8 b[3][3])
{
   real8 det =    (b[0][0] * b[1][1] * b[2][2])           //  det = det(A)
                + (b[0][1] * b[1][2] * b[2][0])           //
                + (b[0][2] * b[1][0] * b[2][1])           //
                - (b[0][2] * b[1][1] * b[2][0])           //
                - (b[0][1] * b[1][0] * b[2][2])           //
                - (b[0][0] * b[1][2] * b[2][1]);          //

   if ( fabs(det) < 1.e-20 )
   {
      printf("%s::%s(%d) - noninvertible matrix encountered det=%le\n", __FILE__, __func__, __LINE__, det );
      return(-1);
   }

   det = ( (fabs(det)>0.0) ? 1.0/det : 0.0 );             // det = 1/det(A)

   a[0][0] = (b[1][1]*b[2][2] - b[1][2]*b[2][1]) * det;   // B = adj(A)/det(A)
   a[0][1] = (b[0][2]*b[2][1] - b[0][1]*b[2][2]) * det;   //
   a[0][2] = (b[0][1]*b[1][2] - b[0][2]*b[1][1]) * det;   //
   a[1][0] = (b[1][2]*b[2][0] - b[1][0]*b[2][2]) * det;   //
   a[1][1] = (b[0][0]*b[2][2] - b[0][2]*b[2][0]) * det;   //
   a[1][2] = (b[0][2]*b[1][0] - b[0][0]*b[1][2]) * det;   //
   a[2][0] = (b[1][0]*b[2][1] - b[1][1]*b[2][0]) * det;   //
   a[2][1] = (b[0][1]*b[2][0] - b[0][0]*b[2][1]) * det;   //
   a[2][2] = (b[0][0]*b[1][1] - b[0][1]*b[1][0]) * det;   //

   return(1);
}

//---------------------------------------------------------------------------------------------------------

real8 Matrix44_Det (real8 a[4][4])
{
   real8 tmp[4][4]; // temp matrix for pairs
   real8 mat[4][4]; // transpose of source matrix
   real8 cfx[4][4]; // matrix of cofactors

   Matrix44_Transpose(mat,a);

   // calculate pairs for first 8 cofactors

   tmp[0][0] = (mat[2][2]*mat[3][3]);
   tmp[0][1] = (mat[2][3]*mat[3][2]);
   tmp[0][2] = (mat[2][1]*mat[3][3]);
   tmp[0][3] = (mat[2][3]*mat[3][1]);
   tmp[1][0] = (mat[2][1]*mat[3][2]);
   tmp[1][1] = (mat[2][2]*mat[3][1]);
   tmp[1][2] = (mat[2][0]*mat[3][3]);
   tmp[1][3] = (mat[2][3]*mat[3][0]);
   tmp[2][0] = (mat[2][0]*mat[3][2]);
   tmp[2][1] = (mat[2][2]*mat[3][0]);
   tmp[2][2] = (mat[2][0]*mat[3][1]);
   tmp[2][3] = (mat[2][1]*mat[3][0]);

   // calculate first 8 cofactors

   cfx[0][0]  = (tmp[0][0]*mat[1][1]) + (tmp[0][3]*mat[1][2]) + (tmp[1][0]*mat[1][3]);
   cfx[0][0] -= (tmp[0][1]*mat[1][1]) + (tmp[0][2]*mat[1][2]) + (tmp[1][1]*mat[1][3]);
   cfx[0][1]  = (tmp[0][1]*mat[1][0]) + (tmp[1][2]*mat[1][2]) + (tmp[2][1]*mat[1][3]);
   cfx[0][1] -= (tmp[0][0]*mat[1][0]) + (tmp[1][3]*mat[1][2]) + (tmp[2][0]*mat[1][3]);
   cfx[0][2]  = (tmp[0][2]*mat[1][0]) + (tmp[1][3]*mat[1][1]) + (tmp[2][2]*mat[1][3]);
   cfx[0][2] -= (tmp[0][3]*mat[1][0]) + (tmp[1][2]*mat[1][1]) + (tmp[2][3]*mat[1][3]);
   cfx[0][3]  = (tmp[1][1]*mat[1][0]) + (tmp[2][0]*mat[1][1]) + (tmp[2][3]*mat[1][2]);
   cfx[0][3] -= (tmp[1][0]*mat[1][0]) + (tmp[2][1]*mat[1][1]) + (tmp[2][2]*mat[1][2]);
   cfx[1][0]  = (tmp[0][1]*mat[0][1]) + (tmp[0][2]*mat[0][2]) + (tmp[1][1]*mat[0][3]);
   cfx[1][0] -= (tmp[0][0]*mat[0][1]) + (tmp[0][3]*mat[0][2]) + (tmp[1][0]*mat[0][3]);
   cfx[1][1]  = (tmp[0][0]*mat[0][0]) + (tmp[1][3]*mat[0][2]) + (tmp[2][0]*mat[0][3]);
   cfx[1][1] -= (tmp[0][1]*mat[0][0]) + (tmp[1][2]*mat[0][2]) + (tmp[2][1]*mat[0][3]);
   cfx[1][2]  = (tmp[0][3]*mat[0][0]) + (tmp[1][2]*mat[0][1]) + (tmp[2][3]*mat[0][3]);
   cfx[1][2] -= (tmp[0][2]*mat[0][0]) + (tmp[1][3]*mat[0][1]) + (tmp[2][2]*mat[0][3]);
   cfx[1][3]  = (tmp[1][0]*mat[0][0]) + (tmp[2][1]*mat[0][1]) + (tmp[2][2]*mat[0][2]);
   cfx[1][3] -= (tmp[1][1]*mat[0][0]) + (tmp[2][0]*mat[0][1]) + (tmp[2][3]*mat[0][2]);

   // calculate pairs for second 8 cofactors

   tmp[0][0] = (mat[0][2]*mat[1][3]);
   tmp[0][1] = (mat[0][3]*mat[1][2]);
   tmp[0][2] = (mat[0][1]*mat[1][3]);
   tmp[0][3] = (mat[0][3]*mat[1][1]);
   tmp[1][0] = (mat[0][1]*mat[1][2]);
   tmp[1][1] = (mat[0][2]*mat[1][1]);
   tmp[1][2] = (mat[0][0]*mat[1][3]);
   tmp[1][3] = (mat[0][3]*mat[1][0]);
   tmp[2][0] = (mat[0][0]*mat[1][2]);
   tmp[2][1] = (mat[0][2]*mat[1][0]);
   tmp[2][2] = (mat[0][0]*mat[1][1]);
   tmp[2][3] = (mat[0][1]*mat[1][0]);

   // calculate second 8 cofactors

   cfx[2][0]  = (tmp[0][0]*mat[3][1]) + (tmp[0][3]*mat[3][2]) + (tmp[1][0]*mat[3][3]);
   cfx[2][0] -= (tmp[0][1]*mat[3][1]) + (tmp[0][2]*mat[3][2]) + (tmp[1][1]*mat[3][3]);
   cfx[2][1]  = (tmp[0][1]*mat[3][0]) + (tmp[1][2]*mat[3][2]) + (tmp[2][1]*mat[3][3]);
   cfx[2][1] -= (tmp[0][0]*mat[3][0]) + (tmp[1][3]*mat[3][2]) + (tmp[2][0]*mat[3][3]);
   cfx[2][2]  = (tmp[0][2]*mat[3][0]) + (tmp[1][3]*mat[3][1]) + (tmp[2][2]*mat[3][3]);
   cfx[2][2] -= (tmp[0][3]*mat[3][0]) + (tmp[1][2]*mat[3][1]) + (tmp[2][3]*mat[3][3]);
   cfx[2][3]  = (tmp[1][1]*mat[3][0]) + (tmp[2][0]*mat[3][1]) + (tmp[2][3]*mat[3][2]);
   cfx[2][3] -= (tmp[1][0]*mat[3][0]) + (tmp[2][1]*mat[3][1]) + (tmp[2][2]*mat[3][2]);
   cfx[3][0]  = (tmp[0][2]*mat[2][2]) + (tmp[1][1]*mat[2][3]) + (tmp[0][1]*mat[2][1]);
   cfx[3][0] -= (tmp[1][0]*mat[2][3]) + (tmp[0][0]*mat[2][1]) + (tmp[0][3]*mat[2][2]);
   cfx[3][1]  = (tmp[2][0]*mat[2][3]) + (tmp[0][0]*mat[2][0]) + (tmp[1][3]*mat[2][2]);
   cfx[3][1] -= (tmp[1][2]*mat[2][2]) + (tmp[2][1]*mat[2][3]) + (tmp[0][1]*mat[2][0]);
   cfx[3][2]  = (tmp[1][2]*mat[2][1]) + (tmp[2][3]*mat[2][3]) + (tmp[0][3]*mat[2][0]);
   cfx[3][2] -= (tmp[2][2]*mat[2][3]) + (tmp[0][2]*mat[2][0]) + (tmp[1][3]*mat[2][1]);
   cfx[3][3]  = (tmp[2][2]*mat[2][2]) + (tmp[1][0]*mat[2][0]) + (tmp[2][1]*mat[2][1]);
   cfx[3][3] -= (tmp[2][0]*mat[2][1]) + (tmp[2][3]*mat[2][2]) + (tmp[1][1]*mat[2][0]);

   // calculate determinant

   real8 det =   (mat[0][0]*cfx[0][0])
               + (mat[0][1]*cfx[0][1])
               + (mat[0][2]*cfx[0][2])
               + (mat[0][3]*cfx[0][3]);

   return (det);
}

//---------------------------------------------------------------------------------------------------------

real8 Matrix44_Norm (real8 a[4][4])
{
   real8 norm=0.0;

   for (int i=0; (i<4); i++)
   {
      real8 sum  = fabs(a[i][0]);
            sum += fabs(a[i][1]);
            sum += fabs(a[i][2]);
            sum += fabs(a[i][3]);  norm = ( (sum>norm) ? sum : norm );
   }

   return (norm);
}

//---------------------------------------------------------------------------------------------------------

real8 Matrix44_Condition (real8 a[4][4])
{
   real8 tmp[4][4];

   Matrix44_Inverse(tmp,a);

   return (   Matrix44_Norm(tmp)
            * Matrix44_Norm(a  ) );
}

//---------------------------------------------------------------------------------------------------------

void Matrix44_Adjoint (real8 a[4][4], real8 b[4][4])
{
   real8 tmp[4][4]; // temp matrix for pairs
   real8 mbt[4][4]; // transpose of source matrix

   Matrix44_Transpose(mbt,b);

   // calculate pairs for first 8 cofactors

   tmp[0][0] = (mbt[2][2]*mbt[3][3]);
   tmp[0][1] = (mbt[2][3]*mbt[3][2]);
   tmp[0][2] = (mbt[2][1]*mbt[3][3]);
   tmp[0][3] = (mbt[2][3]*mbt[3][1]);
   tmp[1][0] = (mbt[2][1]*mbt[3][2]);
   tmp[1][1] = (mbt[2][2]*mbt[3][1]);
   tmp[1][2] = (mbt[2][0]*mbt[3][3]);
   tmp[1][3] = (mbt[2][3]*mbt[3][0]);
   tmp[2][0] = (mbt[2][0]*mbt[3][2]);
   tmp[2][1] = (mbt[2][2]*mbt[3][0]);
   tmp[2][2] = (mbt[2][0]*mbt[3][1]);
   tmp[2][3] = (mbt[2][1]*mbt[3][0]);

   // calculate first 8 cofactors

   a[0][0]  = (tmp[0][0]*mbt[1][1]) + (tmp[0][3]*mbt[1][2]) + (tmp[1][0]*mbt[1][3]);
   a[0][0] -= (tmp[0][1]*mbt[1][1]) + (tmp[0][2]*mbt[1][2]) + (tmp[1][1]*mbt[1][3]);
   a[0][1]  = (tmp[0][1]*mbt[1][0]) + (tmp[1][2]*mbt[1][2]) + (tmp[2][1]*mbt[1][3]);
   a[0][1] -= (tmp[0][0]*mbt[1][0]) + (tmp[1][3]*mbt[1][2]) + (tmp[2][0]*mbt[1][3]);
   a[0][2]  = (tmp[0][2]*mbt[1][0]) + (tmp[1][3]*mbt[1][1]) + (tmp[2][2]*mbt[1][3]);
   a[0][2] -= (tmp[0][3]*mbt[1][0]) + (tmp[1][2]*mbt[1][1]) + (tmp[2][3]*mbt[1][3]);
   a[0][3]  = (tmp[1][1]*mbt[1][0]) + (tmp[2][0]*mbt[1][1]) + (tmp[2][3]*mbt[1][2]);
   a[0][3] -= (tmp[1][0]*mbt[1][0]) + (tmp[2][1]*mbt[1][1]) + (tmp[2][2]*mbt[1][2]);
   a[1][0]  = (tmp[0][1]*mbt[0][1]) + (tmp[0][2]*mbt[0][2]) + (tmp[1][1]*mbt[0][3]);
   a[1][0] -= (tmp[0][0]*mbt[0][1]) + (tmp[0][3]*mbt[0][2]) + (tmp[1][0]*mbt[0][3]);
   a[1][1]  = (tmp[0][0]*mbt[0][0]) + (tmp[1][3]*mbt[0][2]) + (tmp[2][0]*mbt[0][3]);
   a[1][1] -= (tmp[0][1]*mbt[0][0]) + (tmp[1][2]*mbt[0][2]) + (tmp[2][1]*mbt[0][3]);
   a[1][2]  = (tmp[0][3]*mbt[0][0]) + (tmp[1][2]*mbt[0][1]) + (tmp[2][3]*mbt[0][3]);
   a[1][2] -= (tmp[0][2]*mbt[0][0]) + (tmp[1][3]*mbt[0][1]) + (tmp[2][2]*mbt[0][3]);
   a[1][3]  = (tmp[1][0]*mbt[0][0]) + (tmp[2][1]*mbt[0][1]) + (tmp[2][2]*mbt[0][2]);
   a[1][3] -= (tmp[1][1]*mbt[0][0]) + (tmp[2][0]*mbt[0][1]) + (tmp[2][3]*mbt[0][2]);

   // calculate pairs for second 8 cofactors

   tmp[0][0] = (mbt[0][2]*mbt[1][3]);
   tmp[0][1] = (mbt[0][3]*mbt[1][2]);
   tmp[0][2] = (mbt[0][1]*mbt[1][3]);
   tmp[0][3] = (mbt[0][3]*mbt[1][1]);
   tmp[1][0] = (mbt[0][1]*mbt[1][2]);
   tmp[1][1] = (mbt[0][2]*mbt[1][1]);
   tmp[1][2] = (mbt[0][0]*mbt[1][3]);
   tmp[1][3] = (mbt[0][3]*mbt[1][0]);
   tmp[2][0] = (mbt[0][0]*mbt[1][2]);
   tmp[2][1] = (mbt[0][2]*mbt[1][0]);
   tmp[2][2] = (mbt[0][0]*mbt[1][1]);
   tmp[2][3] = (mbt[0][1]*mbt[1][0]);

   // calculate second 8 cofactors

   a[2][0]  = (tmp[0][0]*mbt[3][1]) + (tmp[0][3]*mbt[3][2]) + (tmp[1][0]*mbt[3][3]);
   a[2][0] -= (tmp[0][1]*mbt[3][1]) + (tmp[0][2]*mbt[3][2]) + (tmp[1][1]*mbt[3][3]);
   a[2][1]  = (tmp[0][1]*mbt[3][0]) + (tmp[1][2]*mbt[3][2]) + (tmp[2][1]*mbt[3][3]);
   a[2][1] -= (tmp[0][0]*mbt[3][0]) + (tmp[1][3]*mbt[3][2]) + (tmp[2][0]*mbt[3][3]);
   a[2][2]  = (tmp[0][2]*mbt[3][0]) + (tmp[1][3]*mbt[3][1]) + (tmp[2][2]*mbt[3][3]);
   a[2][2] -= (tmp[0][3]*mbt[3][0]) + (tmp[1][2]*mbt[3][1]) + (tmp[2][3]*mbt[3][3]);
   a[2][3]  = (tmp[1][1]*mbt[3][0]) + (tmp[2][0]*mbt[3][1]) + (tmp[2][3]*mbt[3][2]);
   a[2][3] -= (tmp[1][0]*mbt[3][0]) + (tmp[2][1]*mbt[3][1]) + (tmp[2][2]*mbt[3][2]);
   a[3][0]  = (tmp[0][2]*mbt[2][2]) + (tmp[1][1]*mbt[2][3]) + (tmp[0][1]*mbt[2][1]);
   a[3][0] -= (tmp[1][0]*mbt[2][3]) + (tmp[0][0]*mbt[2][1]) + (tmp[0][3]*mbt[2][2]);
   a[3][1]  = (tmp[2][0]*mbt[2][3]) + (tmp[0][0]*mbt[2][0]) + (tmp[1][3]*mbt[2][2]);
   a[3][1] -= (tmp[1][2]*mbt[2][2]) + (tmp[2][1]*mbt[2][3]) + (tmp[0][1]*mbt[2][0]);
   a[3][2]  = (tmp[1][2]*mbt[2][1]) + (tmp[2][3]*mbt[2][3]) + (tmp[0][3]*mbt[2][0]);
   a[3][2] -= (tmp[2][2]*mbt[2][3]) + (tmp[0][2]*mbt[2][0]) + (tmp[1][3]*mbt[2][1]);
   a[3][3]  = (tmp[2][2]*mbt[2][2]) + (tmp[1][0]*mbt[2][0]) + (tmp[2][1]*mbt[2][1]);
   a[3][3] -= (tmp[2][0]*mbt[2][1]) + (tmp[2][3]*mbt[2][2]) + (tmp[1][1]*mbt[2][0]);
}


// Matrix44_Inverse() : A = Inv(B)
//
// Determines inverse using Kramer's rule
//---------------------------------------------------------------------------------------------------------

int Matrix44_Inverse(real8 a[4][4], real8 b[4][4])
{
   real8 tmp[4][4]; // temp matrix for cofactor pairs
   real8 mbt[4][4]; // transpose of source matrix
   real8 cfx[4][4]; // matrix of cofactors

   Matrix44_Transpose(mbt,b);

   // calculate pairs for first 8 cofactors

   tmp[0][0] = (mbt[2][2]*mbt[3][3]);
   tmp[0][1] = (mbt[2][3]*mbt[3][2]);
   tmp[0][2] = (mbt[2][1]*mbt[3][3]);
   tmp[0][3] = (mbt[2][3]*mbt[3][1]);
   tmp[1][0] = (mbt[2][1]*mbt[3][2]);
   tmp[1][1] = (mbt[2][2]*mbt[3][1]);
   tmp[1][2] = (mbt[2][0]*mbt[3][3]);
   tmp[1][3] = (mbt[2][3]*mbt[3][0]);
   tmp[2][0] = (mbt[2][0]*mbt[3][2]);
   tmp[2][1] = (mbt[2][2]*mbt[3][0]);
   tmp[2][2] = (mbt[2][0]*mbt[3][1]);
   tmp[2][3] = (mbt[2][1]*mbt[3][0]);

   // calculate first 8 cofactors

   cfx[0][0]  = (tmp[0][0]*mbt[1][1]) + (tmp[0][3]*mbt[1][2]) + (tmp[1][0]*mbt[1][3]);
   cfx[0][0] -= (tmp[0][1]*mbt[1][1]) + (tmp[0][2]*mbt[1][2]) + (tmp[1][1]*mbt[1][3]);
   cfx[0][1]  = (tmp[0][1]*mbt[1][0]) + (tmp[1][2]*mbt[1][2]) + (tmp[2][1]*mbt[1][3]);
   cfx[0][1] -= (tmp[0][0]*mbt[1][0]) + (tmp[1][3]*mbt[1][2]) + (tmp[2][0]*mbt[1][3]);
   cfx[0][2]  = (tmp[0][2]*mbt[1][0]) + (tmp[1][3]*mbt[1][1]) + (tmp[2][2]*mbt[1][3]);
   cfx[0][2] -= (tmp[0][3]*mbt[1][0]) + (tmp[1][2]*mbt[1][1]) + (tmp[2][3]*mbt[1][3]);
   cfx[0][3]  = (tmp[1][1]*mbt[1][0]) + (tmp[2][0]*mbt[1][1]) + (tmp[2][3]*mbt[1][2]);
   cfx[0][3] -= (tmp[1][0]*mbt[1][0]) + (tmp[2][1]*mbt[1][1]) + (tmp[2][2]*mbt[1][2]);
   cfx[1][0]  = (tmp[0][1]*mbt[0][1]) + (tmp[0][2]*mbt[0][2]) + (tmp[1][1]*mbt[0][3]);
   cfx[1][0] -= (tmp[0][0]*mbt[0][1]) + (tmp[0][3]*mbt[0][2]) + (tmp[1][0]*mbt[0][3]);
   cfx[1][1]  = (tmp[0][0]*mbt[0][0]) + (tmp[1][3]*mbt[0][2]) + (tmp[2][0]*mbt[0][3]);
   cfx[1][1] -= (tmp[0][1]*mbt[0][0]) + (tmp[1][2]*mbt[0][2]) + (tmp[2][1]*mbt[0][3]);
   cfx[1][2]  = (tmp[0][3]*mbt[0][0]) + (tmp[1][2]*mbt[0][1]) + (tmp[2][3]*mbt[0][3]);
   cfx[1][2] -= (tmp[0][2]*mbt[0][0]) + (tmp[1][3]*mbt[0][1]) + (tmp[2][2]*mbt[0][3]);
   cfx[1][3]  = (tmp[1][0]*mbt[0][0]) + (tmp[2][1]*mbt[0][1]) + (tmp[2][2]*mbt[0][2]);
   cfx[1][3] -= (tmp[1][1]*mbt[0][0]) + (tmp[2][0]*mbt[0][1]) + (tmp[2][3]*mbt[0][2]);

   // calculate pairs for second 8 cofactors

   tmp[0][0] = (mbt[0][2]*mbt[1][3]);
   tmp[0][1] = (mbt[0][3]*mbt[1][2]);
   tmp[0][2] = (mbt[0][1]*mbt[1][3]);
   tmp[0][3] = (mbt[0][3]*mbt[1][1]);
   tmp[1][0] = (mbt[0][1]*mbt[1][2]);
   tmp[1][1] = (mbt[0][2]*mbt[1][1]);
   tmp[1][2] = (mbt[0][0]*mbt[1][3]);
   tmp[1][3] = (mbt[0][3]*mbt[1][0]);
   tmp[2][0] = (mbt[0][0]*mbt[1][2]);
   tmp[2][1] = (mbt[0][2]*mbt[1][0]);
   tmp[2][2] = (mbt[0][0]*mbt[1][1]);
   tmp[2][3] = (mbt[0][1]*mbt[1][0]);

   // calculate second 8 cofactors

   cfx[2][0]  = (tmp[0][0]*mbt[3][1]) + (tmp[0][3]*mbt[3][2]) + (tmp[1][0]*mbt[3][3]);
   cfx[2][0] -= (tmp[0][1]*mbt[3][1]) + (tmp[0][2]*mbt[3][2]) + (tmp[1][1]*mbt[3][3]);
   cfx[2][1]  = (tmp[0][1]*mbt[3][0]) + (tmp[1][2]*mbt[3][2]) + (tmp[2][1]*mbt[3][3]);
   cfx[2][1] -= (tmp[0][0]*mbt[3][0]) + (tmp[1][3]*mbt[3][2]) + (tmp[2][0]*mbt[3][3]);
   cfx[2][2]  = (tmp[0][2]*mbt[3][0]) + (tmp[1][3]*mbt[3][1]) + (tmp[2][2]*mbt[3][3]);
   cfx[2][2] -= (tmp[0][3]*mbt[3][0]) + (tmp[1][2]*mbt[3][1]) + (tmp[2][3]*mbt[3][3]);
   cfx[2][3]  = (tmp[1][1]*mbt[3][0]) + (tmp[2][0]*mbt[3][1]) + (tmp[2][3]*mbt[3][2]);
   cfx[2][3] -= (tmp[1][0]*mbt[3][0]) + (tmp[2][1]*mbt[3][1]) + (tmp[2][2]*mbt[3][2]);
   cfx[3][0]  = (tmp[0][2]*mbt[2][2]) + (tmp[1][1]*mbt[2][3]) + (tmp[0][1]*mbt[2][1]);
   cfx[3][0] -= (tmp[1][0]*mbt[2][3]) + (tmp[0][0]*mbt[2][1]) + (tmp[0][3]*mbt[2][2]);
   cfx[3][1]  = (tmp[2][0]*mbt[2][3]) + (tmp[0][0]*mbt[2][0]) + (tmp[1][3]*mbt[2][2]);
   cfx[3][1] -= (tmp[1][2]*mbt[2][2]) + (tmp[2][1]*mbt[2][3]) + (tmp[0][1]*mbt[2][0]);
   cfx[3][2]  = (tmp[1][2]*mbt[2][1]) + (tmp[2][3]*mbt[2][3]) + (tmp[0][3]*mbt[2][0]);
   cfx[3][2] -= (tmp[2][2]*mbt[2][3]) + (tmp[0][2]*mbt[2][0]) + (tmp[1][3]*mbt[2][1]);
   cfx[3][3]  = (tmp[2][2]*mbt[2][2]) + (tmp[1][0]*mbt[2][0]) + (tmp[2][1]*mbt[2][1]);
   cfx[3][3] -= (tmp[2][0]*mbt[2][1]) + (tmp[2][3]*mbt[2][2]) + (tmp[1][1]*mbt[2][0]);

   // calculate determinant

   real8 det =   (mbt[0][0]*cfx[0][0])
               + (mbt[0][1]*cfx[0][1])
               + (mbt[0][2]*cfx[0][2])
               + (mbt[0][3]*cfx[0][3]);

   if ( fabs(det) < 1.e-20 )
   {
      printf("%s::%s(%d) - noninvertible matrix encountered det=%le\n", __FILE__, __func__, __LINE__, det );
      return(-1);
   }

   det = ( (fabs(det)>0.0) ? 1.0/det : 0.0 );    // det = 1/det(A)

   // calculate matrix inverse

   a[0][0] = cfx[0][0]*det;
   a[0][1] = cfx[0][1]*det;
   a[0][2] = cfx[0][2]*det;
   a[0][3] = cfx[0][3]*det;
   a[1][0] = cfx[1][0]*det;
   a[1][1] = cfx[1][1]*det;
   a[1][2] = cfx[1][2]*det;
   a[1][3] = cfx[1][3]*det;
   a[2][0] = cfx[2][0]*det;
   a[2][1] = cfx[2][1]*det;
   a[2][2] = cfx[2][2]*det;
   a[2][3] = cfx[2][3]*det;
   a[3][0] = cfx[3][0]*det;
   a[3][1] = cfx[3][1]*det;
   a[3][2] = cfx[3][2]*det;
   a[3][3] = cfx[3][3]*det;

   return(1);
}

// MatrixMN_Alloc()
//
// Allocates a dense, zero-indexed, row-major  matrix of size M rows x N columns and initializes to zero.
// Note that you cannot access the resulting matrix using double-subscripts (e.g. a[i][j])
// Matrices allocated with this function can be free'd directly.
//---------------------------------------------------------------------------------------------------------

real8 *MatrixMN_Alloc(const int m, const int n)
{
   real8 *mtx = ( (m>0) && (n>0) ? new real8[m*n] : 0 );

   MatrixMN_Zero(mtx,m,n);

   return(mtx);
}

real8 *MatrixMN_Alloc(const int m0, const int m1, const int n0, const int n1)
{
   const int m = (m1-m0+1);
   const int n = (n1-n0+1);

   real8 *mtx = ( (m>0) && (n>0) ? new real8[m*n] : 0 );

   MatrixMN_Zero(mtx,m,n);

   return(mtx);
}

// MatrixMN_Iliffe()
//
// Given a pointer to an existing, dense MxN matrix, will return an Iliffe vector that allows the
// matrix to be accessed via double subscripts (e.g.  a[i][j]).  Note that this routine
// populates the returned vector such that in can be addressed via non-zero-indexed arguments.
//
// This routine is typically used to create subscript vectors that allow us to pass matrices
// that are interfaced to Fortran routines that expect 1-indexed arguments.
//
// If you use this routine to create non-zero-indexed Iliffe vector, you must offset that
// pointer when you delete the vector (e.g.  delete [] (vec+m0);
//
// Ultimately - I would like to replace all these routines with properly constructed instances
// of C++ matrix classes, until then - we need this as a crutch (bmw).
//---------------------------------------------------------------------------------------------------------

real8 **MatrixMN_Iliffe (const int m , const int n )
{
   real8   *mtx = ( (m>0) && (n>0) ? new real8  [m*n] : 0 );
   real8  **r   = (  mtx           ? new real8 *[m  ] : 0 );

   if (r)
   { for (int i=0,k=0; (i<m); i++, k+=n) { r[i]=(real8 *) mtx+k; } }

   return(r);
}

real8 **MatrixMN_Iliffe (const int m0, const int m1, const int n0, const int n1)
{
   const int m = ( (m1>=m0) ? (m1-m0+1) : 0 );   // m = number of rows
   const int n = ( (n1>=n0) ? (n1-n0+1) : 0 );   // n = number of columns

   real8   *mtx = ( (m>0) && (n>0) ? new real8  [m*n] : 0 );
   real8  **r   = (  mtx           ? new real8 *[m  ] : 0 );

   if (r)
   { for (int i=0,k=-n0; (i<m); i++, k+=n) { r[i]=(real8 *) mtx+k; } }

   return( r ? (r-m0) : 0 );
}

real8 **MatrixMN_Iliffe (const real8 *mtx, const int m , const int n )
{
   real8 **r = ( (mtx && (m>0)) ? new real8 *[m] : 0 );

   if (r)
   { for (int i=0,k=0; (i<m); i++, k+=n) { r[i]=(real8 *) mtx+k; } }

   return(r);
}

real8 **MatrixMN_Iliffe(const real8 *mtx, const int m0, const int m1, const int n0, const int n1)
{
   const int m = ( (m1>=m0) ? (m1-m0+1) : 0 );   // m = number of rows
   const int n = ( (n1>=n0) ? (n1-n0+1) : 0 );   // n = number of columns

   real8 **r = ( (mtx && (m>0)) ? new real8 *[m] : 0 );

   if (r)
   { for (int i=0,k=-n0; (i<m); i++, k+=n) { r[i]=(real8 *) mtx+k; } }

   return( r ? (r-m0) : 0 );
}

// MatrixMN_Free()
//
// Releases memory allocated vi MatrixMN_Alloc().
// Supports release of both iliffe and underlying matrix allocations.
//---------------------------------------------------------------------------------------------------------

void MatrixMN_Free (real8  *a)
{
   if (a) { delete [] a; }
}

void MatrixMN_Free (real8 **a)
{
   real8 *am = ( a ? a[0] : 0 );

   if (a ) { delete [] a ; }
   if (am) { delete [] am; }
}

void MatrixMN_Free (real8 **a, real8 *am)
{
   if (a ) { delete [] a ; }
   if (am) { delete [] am; }
}

void MatrixMN_Free (real8 **a,            const int m0, const int m1, const int n0, const int n1)
{
   real8 *am = ( a ? a[m0]+n0 : 0 );

   if (a ) { delete [] (a+m0); }
   if (am) { delete [] (am  ); }
}

void MatrixMN_Free (real8 **a, real8 *am, const int m0, const int m1, const int n0, const int n1)
{
   if (a ) { delete [] (a+m0); }
   if (am) { delete [] (am  ); }
}

// MatrixMN_Near()
//
// Returns 1 or 0 if two matrices are equal.  Assumes both matrices are the same dimension.
//---------------------------------------------------------------------------------------------------------

int MatrixMN_Near
(
   const real8 *a,     ///< points to A matrix (dense)
   const real8  s,     ///< scalar value
   const int    m,     ///< number of rows
   const int    n,     ///< number of columns
   const real8  tol    ///< comparison tolerance
)
{
   if (a && (m>0) && (n>0) )
   {
      for (int i=0,k=0; i<m; ++i)
      for (int j=0;     j<n; ++j,++k)
      {
        if ( !nearly_eq(a[k],s,tol) ) return(0);
      }

      return(1);
   }

   return(0);
}

int MatrixMN_Near
(
         real8 **a,    ///< points to A matrix (Iliffe matrix)
   const real8   s,    ///< scalar value
   const int     m,    ///< number of rows
   const int     n,    ///< number of columns
   const real8   tol   ///< comparison tolerance
)
{
   if (a && (m>0) && (n>0) )
   {
      for (int i=0; i<m; ++i)
      for (int j=0; j<n; ++j)
      {
        if ( !nearly_eq(a[i][j],s,tol) ) return(0);
      }

      return(1);
   }

   return(0);
}

int MatrixMN_Near
(
         real8 **a ,   ///< points to A matrix (Iliffe matrix)
   const real8   s,    ///< scalar value
   const int     m0,   ///< initial row index
   const int     m1,   ///< final   row index
   const int     n0,   ///< initial column index
   const int     n1,   ///< final   column index
   const real8   tol   ///< comparison tolerance
)
{
   const int m = (m1-m0+1);
   const int n = (n1-n0+1);

   if (a && (m>0) && (n>0) )
   {
      for (int i=m0; i<=m1; ++i)
      for (int j=n0; j<=n1; ++j)
      {
         if ( !nearly_eq(a[i][j],s,tol) ) return(0);
      }

      return(1);
   }

   return(0);
}

int MatrixMN_Near
(
   const real8 *a,     ///< points to A matrix (dense)
   const real8 *b,     ///< points to B matrix (dense)
   const int    m,     ///< number of rows
   const int    n,     ///< number of columns
   const real8  tol    ///< comparison tolerance
)
{
   if (a && b && (m>0) && (n>0) )
   {
      for (int i=0,k=0; i<m; ++i)
      for (int j=0;     j<n; ++j,++k)
      {
         if ( !nearly_eq(a[k],b[k],tol) ) return(0);
      }

      return(1);
   }

   return(0);
}

int MatrixMN_Near
(
         real8 **a ,   ///< points to A matrix (Iliffe matrix)
         real8 **b ,   ///< points to B matrix (Iliffe matrix)
   const int     m ,   ///< number of rows
   const int     n ,   ///< number of columns
   const real8   tol   ///< comparison tolerance
)
{
   if (a && b && (m>0) && (n>0) )
   {
      for (int i=0; i<m; ++i)
      for (int j=0; j<n; ++j)
      {
         if ( !nearly_eq(a[i][j],b[i][j],tol) ) return(0);
      }

      return(1);
   }

   return(0);
}

int MatrixMN_Near
(
         real8 **a ,   ///< points to A matrix (Iliffe matrix)
         real8 **b ,   ///< points to B matrix (Iliffe matrix)
   const int     m0,   ///< initial row index
   const int     m1,   ///< final   row index
   const int     n0,   ///< initial column index
   const int     n1,   ///< final   column index
   const real8   tol   ///< comparison tolerance
)
{
   const int m = (m1-m0+1);
   const int n = (n1-n0+1);

   if (a && b && (m>0) && (n>0) )
   {
      for (int i=m0; i<=m1; ++i)
      for (int j=n0; j<=n1; ++j)
      {
         if ( !nearly_eq(a[i][j],b[i][j],tol) ) return(0);
      }

      return(1);
   }

   return(0);
}

//---------------------------------------------------------------------------------------------------------

real8 MatrixMN_Min   (const real8 *a, const int m, const int n)
{
   real8 min=0.0;

   if (a && (m>0) && (n>0) )
   {
      min=a[0];

      for (int i=0; i<m*n; ++i) { min = ( (a[i]<min) ? a[i] : min ); }
   }

   return(min);
}

real8 MatrixMN_Max   (const real8 *a, const int m, const int n)
{
   real8 max=0.0;

   if (a && (m>0) && (n>0) )
   {
      max=a[0];

      for (int i=0; i<m*n; ++i) { max = ( (a[i]>max) ? a[i] : max ); }
   }

   return(max);
}

real8 MatrixMN_Mean  (const real8 *a, const int m, const int n)
{
   real8 mean=0.0;

   if (a && (m>0) && (n>0) )
   {
      for (int i=0; i<m*n; ++i) { mean += a[i]; }

      mean /= (m*n);
   }

   return(mean);
}

real8 MatrixMN_Range (const real8 *a, const int m, const int n)
{
   real8 min=0.0;
   real8 max=0.0;

   if (a && (m>0) && (n>0) )
   {
      min = max = a[0];

      for (int i=0; i<m*n; ++i) { min = ( (a[i]<min) ? a[i] : min );
                                  max = ( (a[i]>max) ? a[i] : max ); }
   }

   return(max-min);
}


real8 MatrixMN_Min   (int & iy, int & ix, const real8 *a, const int m, const int n)
{
   real8 min=0.0;

   if (a && (m>0) && (n>0) )
   {
      min=a[0]; iy=ix=0;

      for (int i=0,k=0; i<m; ++i)
      for (int j=0;     j<n; ++j, ++k)
      { if (a[k]<min) { min=a[k]; iy=i; ix=j; } }
   }

   return(min);
}

real8 MatrixMN_Max   (int & iy, int & ix, const real8 *a, const int m, const int n)
{
   real8 max=0.0;

   if (a && (m>0) && (n>0) )
   {
      max=a[0]; iy=ix=0;

      for (int i=0,k=0; i<m; ++i)
      for (int j=0;     j<n; ++j, ++k)
      { if (a[k]>max) { max=a[k]; iy=i; ix=j; } }
   }

   return(max);
}


void  MatrixMN_Moment (real8 & min, real8 & max, real8 & mean, real8 & range              , const real8 *a, const int m, const int n)
{
   min=max=mean=range=0.0;

   if (a && (m>0) && (n>0) )
   {
      min=max=a[0];

      for (int i=0; i<m*n; ++i)
      {
         min = ( (a[i]<min) ? a[i] : min );
         max = ( (a[i]>max) ? a[i] : max );
         mean += a[i];
      }

      mean /= (m*n);
      range = (max-min);
   }
}

void  MatrixMN_Moment (real8 & min, real8 & max, real8 & mean, real8 & range, real8 & sdev, const real8 *a, const int m, const int n)
{
   min=max=mean=range=0.0;

   if (a && (m>0) && (n>0) )
   {
      min=max=a[0];

      for (int i=0; i<m*n; ++i)
      {
         min = ( (a[i]<min) ? a[i] : min );
         max = ( (a[i]>max) ? a[i] : max );
         mean += a[i];
      }

      mean /= (m*n);
      range = (max-min);

      real8 var=0.0;
      for (int i=0; i<m*n; ++i)
      {
         real8 s    = (a[i]-mean);
               var += (s*s);
      }

      sdev = sqrt(var/m*n);
   }
}

//---------------------------------------------------------------------------------------------------------

/*-------------------------------------------------------------------------
 *
 *      Function:     MatrixMN_Norm
 *      Description:  Calculate the norm of an MxN matrix
 *
 *      Arguments:
 *          a    MxN array containing components of the matrix
 *
 *------------------------------------------------------------------------*/

real8 MatrixMN_Norm(real8 **a, const int m0, const int m1, const int n0, const int n1)
{
   real8 norm=0.0;

   for (int i=m0; (i<=m1); ++i)
   {
      real8 sum=0.0;
      for (int j=n0; (j<=n1); ++j) { sum += fabs(a[i][j]); }

      norm = ( (sum>norm) ? sum : norm );
   }

   return (norm);
}

#if 0
/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix33Invert
 *      Description:  Calculate the inverse of a 3x3 matrix
 *
 *      Arguments:
 *          a    3x3 array containing components of the matrix to
 *               be inverted.
 *          c    3x3 array in which to return to the caller the
 *               calculated inverse of martix <a>
 *
 *      Returns:  1 on success
 *                0 if the matrix was not invertible
 *
 *------------------------------------------------------------------------*/
int Matrix33Invert(real8 A[3][3], real8 B[3][3])
{
        int   k, j;
        real8 C[3][3], det;

/*
 *      C = adj (A)
 */
        C[0][0] = A[1][1]*A[2][2] - A[1][2]*A[2][1];
        C[1][1] = A[2][2]*A[0][0] - A[2][0]*A[0][2];
        C[2][2] = A[0][0]*A[1][1] - A[0][1]*A[1][0];
        C[0][1] = A[1][2]*A[2][0] - A[1][0]*A[2][2];
        C[1][2] = A[2][0]*A[0][1] - A[2][1]*A[0][0];
        C[2][0] = A[0][1]*A[1][2] - A[0][2]*A[1][1];
        C[0][2] = A[1][0]*A[2][1] - A[2][0]*A[1][1];
        C[1][0] = A[2][1]*A[0][2] - A[0][1]*A[2][2];
        C[2][1] = A[0][2]*A[1][0] - A[1][2]*A[0][0];

        det = A[0][0]*C[0][0] + A[0][1]*C[0][1] + A[0][2]*C[0][2];

        if (fabs(det) < 1e-20) {
            printf("Matrix33Invert: det==0\n");
            return(0);
        }

/*
 *      Transpose and divide by det(A)
 */
        for(j=0;j<3;j++) {
            for(k=0;k<3;k++) {
                B[j][k]=C[k][j]/det;
            }
        }

        return(1);
}
#endif

#if 0
/*
 *      Alternate matrix inverter used for testing.  Apparently
 *      there were some cases where version 9.0 of the intel compiler
 *      (on thunder) generated code for the other matrix inverted
 *      that resulted in a determinant that was significantly off.
 *      (When the numbers involved were on the order of 10e-12 or
 *      smaller)
 *
 *      Returns:  1 on success
 *                0 if the matrix was not invertible
 */
int Matrix33Invert(real8 a[3][3], real8 b[3][3])
{
        int    i, j, k;
        real8  p, fmax, fval, eps = 1.0e-20;
        real8  tmp[3][3];

/*
 *      Initialize the inverse to the identity matrix and create
 *      a temporary copy of the source matrix.
 */
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                b[i][j] = (real8)(i == j);
                tmp[i][j] = a[i][j];
            }
        }

        for (i = 0; i < 3; i++) {
            fmax = fabs(tmp[i][i]);
            for (j = i+1; j < 3; j++) {
                if (fabs(tmp[j][i]) > fmax) {
                    fmax = fabs(tmp[j][i]);
                    for (k = 0; k < 3; k++) {
                        p = tmp[i][k];
                        tmp[i][k] = tmp[j][k];
                        tmp[j][k] = p;
                        p = b[i][k];
                        b[i][k] = b[j][k];
                        b[j][k] = p;
                    }
                }
            }
/*
 *          If can't do the inversion, return 0
 */
            if (fmax < eps) {
#if 0
                printf("Matrix33Invert: fmax < eps, cannot invert!\n");
#endif
                return(0);
            }

            fval = 1.0 / tmp[i][i];

            for (j = 0; j < 3; j++)   {
                tmp[i][j] *= fval;
                b[i][j] *= fval;
            }

            for (k = 0; k < 3; k++) {
                if (k != i) {
                    fval = tmp[k][i];
                    for (j = 0;  j < 3;  j++) {
                        tmp[k][j] -= fval*tmp[i][j];
                        b[k][j] -= fval*b[i][j];
                    }
                }
            }
        }

        return(1);
}
#endif

/*-------------------------------------------------------------------------
 * Matrix33Invert()
 *
 * Computes the inverse of a 3x3 matrix using Cramer's rule.
 * This is the most efficient version of the 3x3 matrix invertion
 * routines.
 *------------------------------------------------------------------------*/

int Matrix33Invert(real8 A[3][3], real8 B[3][3])
{
    real8 det =    (A[0][0] * A[1][1] * A[2][2])       //  det = det(A)
                 + (A[0][1] * A[1][2] * A[2][0])       //
                 + (A[0][2] * A[1][0] * A[2][1])       //
                 - (A[0][2] * A[1][1] * A[2][0])       //
                 - (A[0][1] * A[1][0] * A[2][2])       //
                 - (A[0][0] * A[1][2] * A[2][1]);      //

    if ( fabs(det) < 1.e-20 )
    {
        printf("%s::%s(%d) - noninvertible matrix encountered det=%le\n", __FILE__, __func__, __LINE__, det );
        return(0);
    }

    det = ( (fabs(det)>0.0) ? 1.0/det : 0.0 );             // det = 1/det(A)

    B[0][0] = det * (A[1][1]*A[2][2] - A[1][2]*A[2][1]);   // B = adj(A)/det(A)
    B[0][1] = det * (A[0][2]*A[2][1] - A[0][1]*A[2][2]);   //
    B[0][2] = det * (A[0][1]*A[1][2] - A[0][2]*A[1][1]);   //
    B[1][0] = det * (A[1][2]*A[2][0] - A[1][0]*A[2][2]);   //
    B[1][1] = det * (A[0][0]*A[2][2] - A[0][2]*A[2][0]);   //
    B[1][2] = det * (A[0][2]*A[1][0] - A[0][0]*A[1][2]);   //
    B[2][0] = det * (A[1][0]*A[2][1] - A[1][1]*A[2][0]);   //
    B[2][1] = det * (A[0][1]*A[2][0] - A[0][0]*A[2][1]);   //
    B[2][2] = det * (A[0][0]*A[1][1] - A[0][1]*A[1][0]);   //

    return(1);
}

/*-------------------------------------------------------------------------
 *
 *      Function:     MatrixInvert
 *      Description:  A general matrix inverter
 *
 *      Arguments:
 *          mat     Memory containing matrix to be converted.  Assumes
 *                  components to be inverted reside in rows zero thru
 *                  order-1 and columns zero thru order-1
 *          invMat  Memory into which the inverted matrix will be
 *                  returned to the caller.
 *          lda     Specifies the leading dimension of the matrices <mat>
 *                  and <invMat>
 *          order   Specifies the order of the matrix being inverted.
 *
 *------------------------------------------------------------------------*/
int MatrixInvert(real8 *mat, real8 *invMat, int order, int lda)
{
        int    i, j, k, offset1, offset2, matSize;
        real8  tmpVal, tmpMax, fmax, fval, eps = 1.0e-12;
        real8  *tmpMat;


        matSize = lda * lda * sizeof(real8);
        tmpMat = (real8 *)calloc(1, matSize);

/*
 *      Allocate a temporary array to help form the augmented matrix
 *      for the system.  Initialize tmpMat to be a copy of the input
 *      matrix and the inverse to the identity matrix.
 */

        for (i = 0; i < order; i++) {
            for (j = 0; j < order; j++) {
                offset1 = i * lda + j;
                invMat[offset1] = (real8)(i == j);
                tmpMat[offset1] = mat[offset1];
            }
        }

        for (i = 0; i < order; i++) {

            fmax = fabs(tmpMat[i*lda+i]);
/*
 *          If tmpMat[i][i] is zero, find the next row with a non-zero
 *          entry in column i and switch that row with row i.
 */
            if (fmax < eps) {
                for (j = i+1; j < order; j++) {
                    if ((tmpMax = fabs(tmpMat[j*lda+i])) > fmax) {
                        fmax = tmpMax;
                        for (k = 0; k < order; k++) {
                            offset1 = i * lda + k;
                            offset2 = j * lda + k;
                            tmpVal = tmpMat[offset1];
                            tmpMat[offset1] = tmpMat[offset2];
                            tmpMat[offset2] = tmpVal;
                            tmpVal = invMat[offset1];
                            invMat[offset1] = invMat[offset2];
                            invMat[offset2] = tmpVal;
                        }
                        break;
                    }
                }
            }

/*
 *          If can't do the inversion, return 0
 */
            if (fmax < eps) {
                Fatal("MatrixInvert(): unable to invert matrix!");
            }

/*
 *          Multiply all elements in row i by the inverse of tmpMat[i][i]
 *          to obtain a 1 in tmpMat[i][i]
 */
            fval = 1.0 / tmpMat[i*lda+i];

            for (j = 0; j < order; j++)   {
                offset1 = i * lda + j;
                tmpMat[offset1] *= fval;
                invMat[offset1] *= fval;
            }

/*
 *          Insure that the only non-zero value in column i is in row i
 */
            for (k = 0; k < order; k++) {
                if (k != i) {
                    fval = tmpMat[k*lda+i];
                    for (j = 0; j < order;  j++) {
                        offset1 = k * lda + j;
                        offset2 = i * lda + j;
                        tmpMat[offset1] -= fval*tmpMat[offset2];
                        invMat[offset1] -= fval*invMat[offset2];
                    }
                }
            }

        }   /* for (i = 0; i < order; ...) */

        free(tmpMat);

        return(1);
}


/*------------------------------------------------------------------------
 *
 *      Function:     MatrixMult
 *      Description:  A generic (sort of) matrix multiplier that can
 *                    be used to multiply partially populated
 *                    matrices.
 *
 *          NOTE: The physical memory layout of the input and output
 *                matrices may be larger than the actual portions
 *                of the matrices being multiplied.  The assumption
 *                is made that the components of matrix <a> involved
 *                in the operation reside in rows zero thru <aRows>-1 and
 *                columns zero thru <aCols>-1, and the corresponding
 *                values for matrices <b> and <c>.
 *
 *          a         Pointer to row 0 column 0 of first matrix to be
 *                    multiplied
 *          aRows     Row order of matrix <a>
 *          aCols     Column order of matrix <a>
 *          aLD       Leading dimension of matrix <a>
 *          b         Pointer to row 0 column 0 of second matrix to be
 *                    multiplied
 *          bCols     Column order of matrix <b>
 *          bLD       Leading dimension of matrix <b>
 *          c         Pointer to location in which to store the
 *                    results of the matrix multiply.  Matrix <c> is
 *                    assumed to be of at last dimensions <aRows> X
 *                    <bCols>.
 *          cLD       Leading dimension of matrix <c>
 *
 *----------------------------------------------------------------------*/
void MatrixMult(real8 *a, int aRows, int aCols, int aLD,
                real8 *b, int bCols, int bLD, real8 *c, int cLD)
{
        int  k, m, n;
        int  aCol, bCol, cCol;
        int  aRow, bRow, cRow;
        int  aIndex, bIndex, cIndex;

        for (m = 0; m < aRows; m++) {

            aRow = m;
            cRow = m;

            for (n = 0; n < bCols; n++) {

                bCol = n;
                cCol = n;

                cIndex = cRow * cLD + cCol;
                c[cIndex] = 0.0;

                for (k = 0; k < aCols; k++) {

                    aCol = k;
                    bRow = k;

                    aIndex = aRow * aLD + aCol;
                    bIndex = bRow * bLD + bCol;

                    c[cIndex] += a[aIndex]*b[bIndex];
                }
            }
        }
}

#ifdef ANISOTROPIC
/*---------------------------------------------------------------------------
 *
 *      Function:    MatrixDblMultComplex
 *      Description: Multiplies a matrix of double values by a matrix of
 *                   complex values.
 *
 *          NOTE: The physical memory layout of the input and output
 *                matrices may be larger than the actual portions
 *                of the matrices being multiplied.  The assumption
 *                is made that the components of matrix <a> involved
 *                in the operation reside in rows zero thru <aRows>-1 and
 *                columns zero thru <aCols>-1, and the corresponding
 *                values for matrices <b> and <c>.
 *
 *          a         Pointer to row 0 column 0 of first matrix to be
 *                    multiplied
 *          aRows     Row order of matrix <a>
 *          aCols     Column order of matrix <a>
 *          aLD       Leading dimension of matrix <a>
 *          b         Pointer to row 0 column 0 of second matrix to be
 *                    multiplied
 *          bCols     Column order of matrix <b>
 *          bLD       Leading dimension of matrix <b>
 *
 *          c         Pointer to location in which to store the
 *                    results of the matrix multiply.  Matrix <c> is
 *                    assumed to be of at last dimensions <aRows> X
 *                    <bCols>.
 *          cLD       Leading dimension of matrix <c>
 *
 *-------------------------------------------------------------------------*/
void MatrixDblMultComplex(real8 *a, int aRows, int aCols, int aLD,
                          complex8 *b, int bCols, int bLD,
                          complex8 *c, int cLD)
{
        int  k, m, n;
        int  aCol, bCol, cCol;
        int  aRow, bRow, cRow;
        int  aIndex, bIndex, cIndex;

        for (m = 0; m < aRows; m++) {

            aRow = m;
            cRow = m;

            for (n = 0; n < bCols; n++) {

                bCol = n;
                cCol = n;

                cIndex = cRow * cLD + cCol;
                c[cIndex] = 0.0;

                for (k = 0; k < aCols; k++) {

                    aCol = k;
                    bRow = k;

                    aIndex = aRow * aLD + aCol;
                    bIndex = bRow * bLD + bCol;

                    c[cIndex] += a[aIndex]*b[bIndex];
                }
            }
        }
}
#endif  // (ANISOTROPIC)


/*------------------------------------------------------------------------
 *
 *      Function:     MatrixMultArb
 *      Description:  Another generic matrix multiplier for multiplying
 *                    partially populated matrices.  Unlike MatrixMult()
 *                    above, this function does not require that the
 *                    first component of each matrix being multiplied
 *                    be resident at row 0, column 0.
 *
 *                    Note: Assumes C-type memory allocation such that
 *                          column data varies faster than row data!
 *
 *      Arguments:
 *          a          first matrix for multiply
 *          aCols      number of columns in full <a> matrix
 *          aRowOffset row index of first element of submatrix of matrix <a>
 *                     to be used in the multiply
 *          aColOffset column index of first element of submatrix of
 *                     matrix <a> to be used in the multiply
 *          aMultRows  Number of rows in submatrix of <a> used in the multiply
 *          aMultRows  Number of columns in submatrix of <a> used in the
 *                     multiply
 *          b          second matrix for multiply
 *          bCols      number of columns in full <b> matrix
 *          bRowOffset row index of first element of submatrix of matrix <b>
 *                     to be used in the multiply
 *          bColOffset column index of first element of submatrix of
 *                     matrix <b> to be used in the multiply
 *          bMultCols  Number of columns of submatrix of <b> used in the
 *                     multiply
 *          c          result matrix
 *          cCols      number of columns in full <c> matrix
 *          cRowOffset row index of first element of submatrix of matrix <c>
 *                     in which to store results of the multiply
 *          cColOffset column index of first element of submatrix of
 *                     matrix <c> in which to store results of the multiply
 *
 *----------------------------------------------------------------------*/
void MatrixMultArb(real8 *a, int aCols, int aRowOffset,
                   int aColOffset, int aMultRows, int aMultCols,
                   real8 *b, int bCols, int bRowOffset,
                   int bColOffset, int bMultCols,
                   real8 *c, int cCols, int cRowOffset,
                   int cColOffset)
{
        int  k, m, n;
        int  aRow, bRow, cRow;
        int  aCol, bCol, cCol;
        int  aIndex, bIndex, cIndex;

        for (m = 0; m < aMultRows; m++) {

            aRow = m + aRowOffset;
            cRow = m + cRowOffset;

            for (n = 0; n < bMultCols; n++) {

                bCol = n + bColOffset;
                cCol = n + cColOffset;

                cIndex = cRow * cCols + cCol;
                c[cIndex] = 0.0;

                for (k = 0; k < aMultCols; k++) {

                    aCol = k + aColOffset;
                    bRow = k + bRowOffset;

                    aIndex = aRow * aCols + aCol;
                    bIndex = bRow * bCols + bCol;

                    c[cIndex] += a[aIndex]*b[bIndex];
                }
            }
        }
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix33Vector3Multiply
 *      Description:  Multiply a 3x3 matrix by a 3 element vector
 *
 *      Arguments:
 *          a    3x3 array containing components of the source matrix
 *          x    3 element array containing components of the vector by
 *               which to multiply matrix <a>
 *          y    3 element array in which to return to the caller the
 *               results of multiplying <a> by <x>.
 *
 *------------------------------------------------------------------------*/
void Matrix33Vector3Multiply(real8 A[3][3], real8 x[3], real8 y[3])
{
        y[0] = A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2];
        y[1] = A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2];
        y[2] = A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2];
}

/*-------------------------------------------------------------------------
 *      Function:    Matrix66_Inverse
 *
 *      Description:  Given a 6x6 matrix, calculate the inverse using the
 *                    Gauss-Jordan Method from Numerical Recipes
 *
 *      Arguments:
 *         OUT:  a     inverse of matrix <b>
 *         IN:   b     6x6 matrix to be inverted
 *
 *      Returns:
 *         1 if the matrix was inverted
 *        -1 if the matrix could not be inverted
 *------------------------------------------------------------------------*/
int Matrix66_Inverse(real8 a[6][6], real8 b[6][6])
{
    const int m=1;
    const int n=6;

    int   indxc[6] = { 0 };
    int   indxr[6] = { 0 };
    int   ipiv [6] = { 0 };

    real8 tmp66[6][6] = { 0.0 };

    for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) { tmp66[i][j] = b[i][j]; }      // TMP=B

    // loop over columns of the reduction.

    for (int i=0; i<n; i++)
    {
        real8 big=0.0;
        int   icol=0, irow=0;

        for (int j=0; j<n; j++)
        {
            if (ipiv[j] != 1)
            {
                for (int k=0; k<n; k++)
                {
                    if (ipiv[k] == 0)
                    {
                        if ( fabs(tmp66[j][k]) >= big)
                        {
                            big  = fabs(tmp66[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }  // for k...
            }  // if pivot...
        }

        ipiv[icol]++;

        if (irow != icol)
        {
            for (int l=0; l<n; l++) { SWAP( tmp66[irow][l], tmp66[icol][l] ); }
            for (int l=0; l<m; l++) { SWAP(     a[irow][l],     a[icol][l] ); }
        }

        indxr[i] = irow;
        indxc[i] = icol;

        // If the matrix cannot be inverted, return an error

        if ( tmp66[icol][icol] == 0.0) { return(-1); }  // (error)

        real8 pivinv = 1.0 / tmp66[icol][icol];
        tmp66[icol][icol] = 1.0;

        for (int l=0; l<n; l++) { tmp66[icol][l] *= pivinv; }
        for (int l=0; l<m; l++) {     a[icol][l] *= pivinv; }

        for (int ll=0; ll<n; ll++)
        {
            if (ll != icol)
            {
                real8 s = tmp66[ll][icol];
                          tmp66[ll][icol] = 0.0;

                for (int l=0; l<n; l++) { tmp66[ll][l] -= tmp66[icol][l]*s; }
                for (int l=0; l<m; l++) {     a[ll][l] -=     a[icol][l]*s; }
            }
        }

    } // (end of the main loop over columns of the reduction)

    for (int l = n-1; l >= 0; l--)
    {
        if (indxr[l] != indxc[l])
        {
            for (int k=0; k<n; k++)
            { SWAP( tmp66[k][indxr[l]], tmp66[k][indxc[l]] ); }
        }
    }

    for (int i=0; i<6; i++)
    for (int j=0; j<6; j++) { a[i][j] = tmp66[i][j]; }

    return(1);  // (success)
}

/*------------------------------------------------------------------------
 *
 *      Function:     Vec3TransposeAndMult
 *      Description:  Multiply a 3-element vector by it's transpose
 *
 *      Arguments:
 *          vec     3-element source vector
 *          result  3x3 matrix in which results are returned to the caller
 *
 *------------------------------------------------------------------------*/
void Vec3TransposeAndMult(real8 vec[3], real8 result[3][3])
{
    result[0][0] = vec[0]*vec[0];
    result[0][1] = vec[0]*vec[1];
    result[0][2] = vec[0]*vec[2];

    result[1][0] = vec[1]*vec[0];
    result[1][1] = vec[1]*vec[1];
    result[1][2] = vec[1]*vec[2];

    result[2][0] = vec[2]*vec[0];
    result[2][1] = vec[2]*vec[1];
    result[2][2] = vec[2]*vec[2];
}

//---------------------------------------------------------------------------------------------------------------------------------
// SVD Support...
//---------------------------------------------------------------------------------------------------------------------------------

#ifdef MIN
#undef MIN
#endif

#ifdef MAX
#undef MAX
#endif

static inline int    MIN   (const int   a, const int   b) { return( (a<b)    ? a : b ); }
static inline real8  MAX   (const real8 a, const real8 b) { return( (a>b)    ? a : b ); }
static inline real8  Sign  (const real8 a, const real8 b) { return( (b>=0.0) ? fabs(a) : -fabs(a) ); }

static        real8  Pythag(const real8 a, const real8 b)
{
	real8 absa = fabs(a);
	real8 absb = fabs(b);

	     if (absa>absb) { real8 x = (absb/absa); return( absa*sqrt(1.0+(x*x)) ); }
	else if (absb>0.0 ) { real8 x = (absa/absb); return( absb*sqrt(1.0+(x*x)) ); }

    return(0.0);
}

// MatrixMN_SVD()
//
// Computes the Singular Value Decomposition of the source matrix.
// Upon exit, populates three matrix and vector references that
// satisfy
//    SVD(A) = U*W*V'
//
// On entry:
//   A : MxN zero-indexed Iliffe matrix to be decomposed (input,updated)
//   V : NxN zero-indexed Iliffe matrix
//   W : 1xN zero-indexed vector
//   R : 1xN zero-indexed scratch vector
//
//   The source matrix A is provided as an Iliffe matrix, where
//   A[i] points to each row of the source matrix.  The source matrix
//   is passed by reference and modified by this routine.
//
// On successful exit:
//   a : MxN zero-indexed Iliffe matrix contains the decomposed result
//   v : NxN zero-indexed Iliffe matrix containing the transpose of V (e.g V')
//   w : 1xN vector containing the diagonal elements of the W matrix
//
// Returns...
//   0 : (success)
//  -1 : (failure - iteration count exceeded).
//
// Adapted from Numerical Recipes in C, Second Edition
//--------------------------------------------------------------------------------------

int MatrixMN_SVD
(
          real8 **a,  ///< MxN zero-indexed Iliffe matrix to be decomposed (input,updated)
          real8 **v,  ///< NxN zero-indexed Iliffe matrix
          real8  *w,  ///< 1xN allocated vector
          real8  *r,  ///< 1xN allocated scratch vector
   const  int     m,  ///< size of the source matrix (M-rows)
   const  int     n   ///< size of the source matrix (N-columns)
)
{
    int   i,its,j,jj,k,l,nm;
    real8 anorm,c,f,g,h,s,scale,x,y,z,*rv1;

    g=scale=anorm=0.0;

    for (i=0; (i<n); i++)
    {
        l=i+1;
        r[i]=scale*g;
        g=s=scale=0.0;
        if (i<m)
        {
            for (k=i; (k<m); k++) { scale += fabs(a[k][i]); }

            if (scale)
            {
                for (k=i; (k<m); k++)
                {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }

                f=a[i][i];
                g = -Sign(sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;

                for (j=l; (j<n); j++)
                {
                    for (s=0.0,k=i; (k<m); k++) { s += a[k][i]*a[k][j]; }
                    f=s/h;
                    for (      k=i; (k<m); k++) { a[k][j] += f*a[k][i]; }
                }
                for (k=i; (k<m); k++) { a[k][i] *= scale; }
            }
        }
        w[i]=scale *g;
        g=s=scale=0.0;
        if ( (i<m) && i != n-1)
        {
            for (k=l; (k<n); k++) scale += fabs(a[i][k]);
            if (scale)
            {
                for (k=l; (k<n); k++)
                {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f=a[i][l];
                g = -Sign(sqrt(s),f);
                h=f*g-s;
                a[i][l]=f-g;
                for (k=l; (k<n); k++) { r[k]=a[i][k]/h; }
                for (j=l; (j<m); j++)
                {
                    for (s=0.0,k=l; k<n; k++) { s += a[j][k]*a[i][k]; }
                    for (      k=l; k<n; k++) { a[j][k] += s*r[k];  }
                }
                for (k=l; k<n; k++) { a[i][k] *= scale; }
            }
        }
        anorm = MAX( anorm, (fabs(w[i])+fabs(r[i])) );
    }

    for (i=n-1; (i>=0); i--)
    {
        if (i<(n-1))
        {
            if (g)
            {
                for (j=l; (j<n); j++) { v[j][i]=(a[i][j]/a[i][l])/g; }
                for (j=l; (j<n); j++)
                {
                    for (s=0.0,k=l; (k<n); k++) { s += a[i][k]*v[k][j]; }
                    for (      k=l; (k<n); k++) { v[k][j] += s*v[k][i]; }
                }
            }
            for (j=l; (j<n); j++) { v[i][j]=v[j][i]=0.0; }
        }
        v[i][i]=1.0;
        g=r[i];
        l=i;
    }

    for (i=MIN(m,n)-1; (i>=0); i--)
    {
        l=i+1;
        g=w[i];

        for (j=l; (j<n); j++) { a[i][j]=0.0; }

        if (g)
        {
            g=1.0/g;

            for (j=l; (j<n); j++)
            {
                for (s=0.0,k=l; (k<m); k++) { s += a[k][i]*a[k][j]; }
                f=(s/a[i][i])*g;
                for (      k=i; (k<m); k++) { a[k][j] += f*a[k][i]; }
            }
            for (j=i; (j<m); j++) { a[j][i] *= g; }
        }
        else
            for (j=i; (j<m); j++) { a[j][i] = 0.0; }

        ++a[i][i];
    }

    for (k=(n-1); (k>=0); k--)
    {
        for (its=0; (its<=30); its++)
        {
            int flag=1;

            for (l=k; (l>=0); l--)
            {
                nm=l-1;
                if ( (fabs(r[l] )+anorm) == anorm) { flag=0; break; }
                if ( (fabs(w[nm])+anorm) == anorm) {         break; }
            }

            if (flag)
            {
                c=0.0;
                s=1.0;

                for (i=l; (i<k); i++)
                {
                    f=s*r[i];
                    r[i]=c*r[i];
                    if ( (fabs(f)+anorm) == anorm) break;
                    g=w[i];
                    h=Pythag(f,g);
                    w[i]=h;
                    h=1.0/h;
                    c=g*h;
                    s = -f*h;
                    for (j=0; (j<m); j++)
                    {
                        y=a[j][nm];
                        z=a[j][i ];
                        a[j][nm]=y*c+z*s;
                        a[j][i ]=z*c-y*s;
                    }
                }
            }

            z=w[k];

            if (l==k)
            {
                if (z<0.0)
                {
                    w[k] = -z;
                    for (j=0; (j<n); j++) { v[j][k] = -v[j][k]; }
                }
                break;
            }

            if (its>=30)
            {
               printf("%s::%s(%d) - maximum iteration count reached, solution invalid\n", __FILE__, __func__, __LINE__ );
               return(-1);
            }

            x=w[l];
            nm=k-1;
            y=w[nm];
            g=r[nm];
            h=r[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=Pythag(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+Sign(g,f)))-h))/x;
            c=s=1.0;

            for (j=l; (j<=nm); j++)
            {
                i=j+1;
                g=r[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=Pythag(f,h);
                r[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g = g*c-x*s;
                h=y*s;
                y *= c;

                for (jj=0; (jj<n); jj++)
                {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=Pythag(f,h);
                w[j]=z;

                if (z)
                {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=c*g+s*y;
                x=c*y-s*g;

                for (jj=0; (jj<m); jj++)
                {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }

            r[l]=0.0;
            r[k]=f;
            w[k]=x;
        }
    }

    return(0);
}

//--------------------------------------------------------------------------------------

int MatrixMN_SVD
(
         real8 *amat,    ///< A (dense, MxN matrix, zero-indexed) (unchanged)
         real8 *umat,    ///< U (dense, MxN matrix, zero-indexed) (returned )
         real8 *w   ,    ///< W (dense, 1xN vector, zero-indexed) (returned )
         real8 *vmat,    ///< V (dense, NxN matrix, zero-indexed) (returned )
   const int    m   ,    ///< number of rows    (M)
   const int    n        ///< number of columns (N)
)
{
   real8 **u = ( (umat && (m>0)) ? new real8 *[m] : 0 );
   real8 **v = ( (vmat && (n>0)) ? new real8 *[n] : 0 );
   real8  *r = (          (n>0)  ? new real8  [n] : 0 );

   if (umat && amat) for (int i=0; (i<m*n); ++i) { umat[i]=amat[i]; }   // U = A
   if (vmat        ) for (int i=0; (i<n*n); ++i) { vmat[i]=0.0;     }   // V = 0.0
   if (w           ) for (int i=0; (i<n  ); ++i) { w   [i]=0.0;     }   // W = 0.0
   if (r           ) for (int i=0; (i<n  ); ++i) { r   [i]=0.0;     }   // R = 0.0

   // create Iliffe matrices for SVD interface...

   for (int i=0,k=0; (i<m); i++, k+=n) { u[i] = umat+k; }   // a = A[0..(m-1)][0..(n-1)]
   for (int i=0,k=0; (i<n); i++, k+=n) { v[i] = vmat+k; }   // v = V[0..(n-1)][0..(n-1)]

   int err = MatrixMN_SVD(u,v,w,r,m,n);

   // cleanup and exit the routine...

   if (u) { delete [] u; }
   if (v) { delete [] v; }
   if (r) { delete [] r; }

   return (err);
}

// MatrixMN_SVD_Inverse()
//
// Computes AI=Inverse(A) using singular value decomposition.
// The A source matrix is a dense, zero-indexed, matrix and is unchanged.
//--------------------------------------------------------------------------------------

int MatrixMN_SVD_Inverse
(
         real8 *am,    ///< A (dense, MxN matrix, zero-indexed, row-major)
   const int    m ,    ///< number of rows    (M)
   const int    n      ///< number of columns (N)
)
{
   real8  *tmp = ( (m>0) && (n>0) ? new real8[2*m*n + n*n + n+n ] : 0 );
   real8  *p   = tmp;

   real8 **a   = ( am  ? MatrixMN_Iliffe(am,m,n) : 0 );
   real8 **u   = ( tmp ? MatrixMN_Iliffe(p ,m,n) : 0 ); p+=m*n;
   real8 **ut  = ( tmp ? MatrixMN_Iliffe(p ,n,m) : 0 ); p+=n*m;
   real8 **v   = ( tmp ? MatrixMN_Iliffe(p ,n,n) : 0 ); p+=n*n;
   real8  *w   = ( tmp ?                 p       : 0 ); p+=n;
   real8  *r   = ( tmp ?                 p       : 0 ); p+=n;

   if (u && a) { for(int i=0; (i<m); i++)
                 for(int j=0; (j<n); j++) { u[i][j]=a[i][j]; } }

   if (v)      { for(int i=0; (i<n); i++)
                 for(int j=0; (j<n); j++) { v[i][j]=0.0; } }

   if (w)      { for(int i=0; (i<n); i++) { w[i]=0.0; } }
   if (r)      { for(int i=0; (i<n); i++) { r[i]=0.0; } }

   int err = MatrixMN_SVD(u,v,w,r,m,n);

   if (!err)
   {
      for (int i=0; (i<n); i++)
      { w[i] = ( (fabs(w[i])>0.0) ? 1.0/w[i] : 0.0 ); }       //  W = 1/W

      MatrixMN_Transpose(ut,u,m,n);                           //  U' = Transpose(U)

      for (int i=0; (i<n); i++)                               //
      for (int j=0; (j<n); j++) { v[i][j] *= w[j]; }          //  V      = V x Diag(1/w)
                                                              //
      MatrixMN_Mul(a,v,ut,m,n,n);                             //  Inv(a) = V x Diag(1/w) x u'
   }

   if (tmp) { delete [] tmp; tmp=0; }
   if (a  ) { delete [] a  ; a  =0; }
   if (u  ) { delete [] u  ; u  =0; }
   if (ut ) { delete [] ut ; ut =0; }
   if (v  ) { delete [] v  ; v  =0; }

   return(err);
}

// Matrix33_SVD() : computes SVD(A) = U*W*V'
//
// This routine was specifically added for 3x3 matrices and avoids temporary
// allocations off the heap.
//--------------------------------------------------------------------------------------

int Matrix33_SVD
(
   real8 amat[3][3],    ///< A (3x3 matrix) (source, unchanged)
   real8 umat[3][3],    ///< U (3x3 matrix) (output)
   real8 wvec[3]   ,    ///< W (1x3 vector) (output)
   real8 vmat[3][3]     ///< V (3x3 matrix) (output)
)
{
   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j) { umat[i][j] = amat[i][j]; }        // U = A

   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j) { vmat[i][j] = 0.0; }               // V = 0.0

   for (int i=0; (i<3); ++i) { wvec[i]    = 0.0; }               // W = 0.0

   real8 rvec[3] = { 0.0 };  // (temporary)                      // R = 0.0

   // create Iliffe matrices and vectors for SVD interface...

   real8  *u[3] = { ((real8 *) umat)+0,                          // u = U
                    ((real8 *) umat)+3,                          //
                    ((real8 *) umat)+6 };                        //

   real8  *v[3] = { ((real8 *) vmat)+0,                          // v = A
                    ((real8 *) vmat)+3,                          //
                    ((real8 *) vmat)+6 };                        //

   real8  *w = (real8 *) wvec;                                   // w = W
   real8  *r = (real8 *) rvec;                                   // r = R

   int err = MatrixMN_SVD(u,v,w,r,3,3);

   return (err);
}

// Matrix33_SVD_Inverse()
//
// Will develop a pseudo-inverse via singular value decompostion.
//
// Essentially, performs the logic...
//
//   SVD [A]  = [U]*[W]*[V]'
//   INV([A]) = [V]*[diag(1/w)][U]'
//-------------------------------------------------------------------------------

int Matrix33_SVD_Inverse
(
   real8 a[3][3]   ///< A=Inv(A)
)
{
   real8  u[3][3] = { 0.0 };                           // U=0
   real8  v[3][3] = { 0.0 };                           // V=0
   real8  w[3]    = { 0.0 };                           // W=0

   if ( Matrix33_SVD(a,u,w,v) < 0 ) { return(0); }     // u,w,v = svd(a), on error=non-invertible
                                                       //
   for (int i=0; (i<3); i++)                           //
   { w[i] = ( (fabs(w[i])>0.0) ? 1.0/w[i] : 0.0 ); }   // w = diag(1/w)
                                                       //
   Matrix33_Transpose(u);                              // u = u'
                                                       //
   for (int i=0; (i<3); i++)                           //
   for (int j=0; (j<3); j++) { v[i][j] *= w[j]; }      //  V      = V x Diag(1/w)
                                                       //
   Matrix33_Mul(a,v,u);                                //  Inv(a) = V x Diag(1/w) x u'
                                                       //
   return(1);
}

// Matrix33_SVD_Inverse(alternate form) : A = SVD(B)
//
// Returns zero (failure) or one (success).
//-------------------------------------------------------------------------------

int Matrix33_SVD_Inverse
(
   real8 a[3][3],      ///< A = SVD(B) (result)
   real8 b[3][3]       ///< B          (unchanged)
)
{
   real8  u[3][3] = { 0.0 };                           // U=0
   real8  v[3][3] = { 0.0 };                           // V=0
   real8  w[3]    = { 0.0 };                           // W=0

   for (int i=0; (i<3); ++i)
   for (int j=0; (j<3); ++j) { a[i][j] = b[i][j]; }    // A = B  (copy)
                                                       //
   if ( Matrix33_SVD(a,u,w,v) < 0 ) { return (-1); }   // u,w,v = svd(a), on error=non-invertible
                                                       //
   for (int i=0; (i<3); i++)                           //
   { w[i] = ( (fabs(w[i])>0.0) ? 1.0/w[i] : 0.0 ); }   // w = diag(1/w)
                                                       //
   Matrix33_Transpose(u);                              // u = u'
                                                       //
   for (int i=0; (i<3); i++)                           //
   for (int j=0; (j<3); j++) { v[i][j] *= w[j]; }      //  V      = V x Diag(1/w)
                                                       //
   Matrix33_Mul(a,v,u);                                //  Inv(a) = V x Diag(1/w) x u'
                                                       //
   return (1);
}

// Matrix33_SVD_Condition()
//
// Returns condition number of matrix using pseudo inverse :  cond = ||A|| * ||Inv(A)||
//---------------------------------------------------------------------------------------------------------

real8 Matrix33_SVD_Condition (real8 a[3][3])
{
   real8 tmp[3][3];

   Matrix33_SVD_Inverse(tmp,a);

   return (   Matrix33_Norm(tmp)
            * Matrix33_Norm(a  ) );
}


// Matrix33_PseudoInverse
//
// This is Nicolas Bertin's version of the SVD-based pseudo inverse   A = PseudoInverse(B)
//---------------------------------------------------------------------------------------------------------

static inline real8 fmax3 (const real8 a, const real8 b, const real8 c) { return( (a>b) ? ( (a>c) ? a : c )  : ( (b>c) ? b : c ) ); }

int Matrix33_PseudoInverse
(
   real8 a[3][3],
   real8 b[3][3]
)
{
    real8 u [3][3] = { 0.0 };
    real8 w [3]    = { 0.0 };
    real8 v [3][3] = { 0.0 };
    real8 vs[3][3] = { 0.0 };

    int stat = Matrix33_SVD(b,u,w,v);

    real8 meps = __DBL_EPSILON__;                          // meps = machine precision (double)
    real8 tol  = 3.0*meps*fmax3( w[0], w[1], w[2] );       // tol  = 3*meps*Max(w)

    int rank=0;
    for (int i=0; (i<3); i++) { if (w[i]>tol) rank++; }

    Matrix33_Zero(a);  // B=0  (init)

    if (rank>0)
    {
       for (int i=0; (i<3);    i++)
       for (int j=0; (j<rank); j++) { vs[i][j] = v[i][j]/w[j]; }

       for (int i=0; (i<3); i++)
       for (int j=0; (j<3); j++)
       {
          a[i][j] = 0.0;
          for (int k=0; (k<rank); k++) { a[i][j] += vs[i][k]*u[j][k]; }
       }
    }

    return stat;
}
