#pragma once

#ifndef _PDS_MATRIX_H
#define _PDS_MATRIX_H

//---------------------------------------------------------------------------------------------------------
//  Matrix.h  Contains prototypes for various matrix operation functions
//---------------------------------------------------------------------------------------------------------

void    Matrix22_Zero         (real8 a[2][2] );                                                                          // A = 0
void    Matrix23_Zero         (real8 a[2][3] );                                                                          // A = 0
void    Matrix32_Zero         (real8 a[3][2] );                                                                          // A = 0
void    Matrix33_Zero         (real8 a[3][3] );                                                                          // A = 0
void    Matrix44_Zero         (real8 a[4][4] );                                                                          // A = 0
void    Matrix66_Zero         (real8 a[6][6] );                                                                          // A = 0
void    Matrix333_Zero        (real8 a[3][3][3] );                                                                       // A = 0

void    MatrixMN_Zero         (real8    *a, const int m , const int n );                                                 // A = 0
void    MatrixMN_Zero         (real8   **a, const int m , const int n );                                                 // A = 0
void    MatrixMN_Zero         (real8   **a, const int m0, const int m1, const int n0, const int n1 );                    // A = 0

void    Matrix22_Set          (real8 a[2][2], const real8 s );                                                           // A = s
void    Matrix23_Set          (real8 a[2][3], const real8 s );                                                           // A = s
void    Matrix32_Set          (real8 a[3][2], const real8 s );                                                           // A = s
void    Matrix33_Set          (real8 a[3][3], const real8 s );                                                           // A = s
void    Matrix44_Set          (real8 a[4][4], const real8 s );                                                           // A = s
void    Matrix66_Set          (real8 a[6][6], const real8 s );                                                           // A = s

void    MatrixMN_Set          (real8    *a, const real8   s, const int m , const int n );                                // A = s
void    MatrixMN_Set          (real8    *a, const real8  *s, const int m , const int n );                                // A = s
void    MatrixMN_Set          (real8   **a, const real8   s, const int m , const int n );                                // A = s
void    MatrixMN_Set          (real8   **a, const real8   s, const int m0, const int m1, const int n0, const int n1 );   // A = s
void    MatrixMN_Set          (real8   **a, const real8  *s, const int m0, const int m1, const int n0, const int n1 );   // A = s
void    MatrixMN_Set          (real8   **a,       real8 **s, const int m0, const int m1, const int n0, const int n1 );   // A = S

void    Matrix22_Copy         (real8 a[2][2], real8 b[2][2]);                                                            // A = B (2x2)
void    Matrix33_Copy         (real8 a[3][3], real8 b[3][3]);                                                            // A = B (3x3)
void    Matrix44_Copy         (real8 a[4][4], real8 b[4][4]);                                                            // A = B (4x4)
void    Matrix66_Copy         (real8 a[6][6], real8 b[6][6]);                                                            // A = B (6x6)

void    MatrixMN_Copy         (real8    *a, const real8  *b, const int m , const int n );                                // A = B
void    MatrixMN_Copy         (real8   **a, const real8 **b, const int m , const int n );                                // A = B

void    Matrix22_Identity     (real8 a[2][2] );                                                                          // A = I
void    Matrix33_Identity     (real8 a[3][3] );                                                                          // A = I
void    Matrix44_Identity     (real8 a[4][4] );                                                                          // A = I
void    Matrix66_Identity     (real8 a[6][6] );                                                                          // A = I

void    Matrix22_Threshold    (real8   a[2][2],                           const real8 min);                              // A = fabs(a[i][j])<min ? 0.0 : a[i][j]
void    Matrix33_Threshold    (real8   a[3][3],                           const real8 min);                              // A = fabs(a[i][j])<min ? 0.0 : a[i][j]
void    Matrix44_Threshold    (real8   a[4][4],                           const real8 min);                              // A = fabs(a[i][j])<min ? 0.0 : a[i][j]
void    MatrixMN_Threshold    (real8  *a      , const int m, const int n, const real8 min);                              // A = fabs(a[i]   )<min ? 0.0 : a[i]
void    MatrixMN_Threshold    (real8 **a      , const int m, const int n, const real8 min);                              // A = fabs(a[i][j])<min ? 0.0 : a[i][j]

void    Matrix22_Cap          (real8   a[2][2],                           const real8 max);                              // A = a[i][j]<max ? a[i][j] : max
void    Matrix33_Cap          (real8   a[3][3],                           const real8 max);                              // A = a[i][j]<max ? a[i][j] : max
void    Matrix44_Cap          (real8   a[4][4],                           const real8 max);                              // A = a[i][j]<max ? a[i][j] : max
void    MatrixMN_Cap          (real8  *a      , const int m, const int n, const real8 max);                              // A = a[i]   <max ? a[i]    : max
void    MatrixMN_Cap          (real8 **a      , const int m, const int n, const real8 max);                              // A = a[i][j]<max ? a[i][j] : max

void    Matrix22_Constrain    (real8   a[2][2],                           const real8 min, const real8 max);             // A = a[i][j]<min ? min : ( a[i][j]>max ? max : a[i][j] )
void    Matrix33_Constrain    (real8   a[3][3],                           const real8 min, const real8 max);             // A = a[i][j]<min ? min : ( a[i][j]>max ? max : a[i][j] )
void    Matrix44_Constrain    (real8   a[4][4],                           const real8 min, const real8 max);             // A = a[i][j]<min ? min : ( a[i][j]>max ? max : a[i][j] )
void    MatrixMN_Constrain    (real8  *a      , const int m, const int n, const real8 min, const real8 max);             // A = a[i]   <min ? min : ( a[i]   >max ? max : a[i]    )
void    MatrixMN_Constrain    (real8 **a      , const int m, const int n, const real8 min, const real8 max);             // A = a[i][j]<min ? min : ( a[i][j]>max ? max : a[i][j] )

void    Matrix22_Abs          (real8   a[2][2]                          );                                               // A = Abs(A)   a[i][j] = fabs(a[i][j])
void    Matrix33_Abs          (real8   a[3][3]                          );                                               // A = Abs(A)   a[i][j] = fabs(a[i][j])
void    Matrix44_Abs          (real8   a[4][4]                          );                                               // A = Abs(A)   a[i][j] = fabs(a[i][j])
void    MatrixMN_Abs          (real8  *a      , const int m, const int n);                                               // A = Abs(A)   a[i]    = fabs(a[i]   )
void    MatrixMN_Abs          (real8 **a      , const int m, const int n);                                               // A = Abs(A)   a[i][j] = fabs(a[i][j])

void    Matrix22_Sign         (real8   a[2][2]                          );                                               // A = Sign(A)  a[i][j] = ( a[i][j]<0 ? -1 : ( a[i][j]>0 ? +1 : 0 ) )
void    Matrix33_Sign         (real8   a[3][3]                          );                                               // A = Sign(A)  a[i][j] = ( a[i][j]<0 ? -1 : ( a[i][j]>0 ? +1 : 0 ) )
void    Matrix44_Sign         (real8   a[4][4]                          );                                               // A = Sign(A)  a[i][j] = ( a[i][j]<0 ? -1 : ( a[i][j]>0 ? +1 : 0 ) )
void    MatrixMN_Sign         (real8  *a      , const int m, const int n);                                               // A = Sign(A)  a[i]    = ( a[i]   <0 ? -1 : ( a[i]   >0 ? +1 : 0 ) )
void    MatrixMN_Sign         (real8 **a      , const int m, const int n);                                               // A = Sign(A)  a[i][j] = ( a[i][j]<0 ? -1 : ( a[i][j]>0 ? +1 : 0 ) )

void    Matrix22_Residual     (real8   a[2][2], real8   b[2][2] );                                                       // A = Abs(A-B)  a[i][j] = fabs(a[i][j]-b[i][j])
void    Matrix33_Residual     (real8   a[3][3], real8   b[3][3] );                                                       // A = Abs(A-B)  a[i][j] = fabs(a[i][j]-b[i][j])
void    Matrix44_Residual     (real8   a[4][4], real8   b[4][4] );                                                       // A = Abs(A-B)  a[i][j] = fabs(a[i][j]-b[i][j])
void    MatrixMN_Residual     (real8  *a      , real8  *b      , const int m, const int n);                              // A = Abs(A-B)  a[i]    = fabs(a[i]   -b[i]   )
void    MatrixMN_Residual     (real8 **a      , real8 **b      , const int m, const int n);                              // A = Abs(A-B)  a[i][j] = fabs(a[i][j]-b[i][j])

real8   Matrix22_RMS_Error    (real8   a[2][2], real8   b[2][2] );                                                       // rms = RMS_Error(A-B)
real8   Matrix33_RMS_Error    (real8   a[3][3], real8   b[3][3] );                                                       // rms = RMS_Error(A-B)
real8   Matrix44_RMS_Error    (real8   a[4][4], real8   b[4][4] );                                                       // rms = RMS_Error(A-B)
real8   MatrixMN_RMS_Error    (real8  *a      , real8  *b      , const int m, const int n);                              // rms = RMS_Error(A-B)
real8   MatrixMN_RMS_Error    (real8 **a      , real8 **b      , const int m, const int n);                              // rms = RMS_Error(A-B)

void    Matrix22_Rand         (real8   a[2][2]                          , const real8 min=0.0, const real8 max=1.0);     // A = a[i][j] = srand(min,max),  a[i][j]==[min,max]
void    Matrix33_Rand         (real8   a[3][3]                          , const real8 min=0.0, const real8 max=1.0);     // A = a[i][j] = srand(min,max),  a[i][j]==[min,max]
void    Matrix44_Rand         (real8   a[4][4]                          , const real8 min=0.0, const real8 max=1.0);     // A = a[i][j] = srand(min,max),  a[i][j]==[min,max]
void    MatrixMN_Rand         (real8  *a      , const int m, const int n, const real8 min=0.0, const real8 max=1.0);     // A = a[i]    = srand(min,max),  a[i][j]==[min,max]
void    MatrixMN_Rand         (real8 **a      , const int m, const int n, const real8 min=0.0, const real8 max=1.0);     // A = a[i][j] = srand(min,max),  a[i][j]==[min,max]

void    Matrix22_Rotate_Z     (real8 a[2][2], const real8 rz);                                                           // A = RZ(rz)

void    Matrix33_Rotate_X     (real8 a[3][3], const real8 rx);                                                           // A = RX(rx)
void    Matrix33_Rotate_Y     (real8 a[3][3], const real8 ry);                                                           // A = RY(ry)
void    Matrix33_Rotate_Z     (real8 a[3][3], const real8 rz);                                                           // A = RZ(rz)
void    Matrix33_Rotate       (real8 a[3][3], const real8 rx, const real8 ry, const real8 rz);                           // A = Rot(rx,ry,rz)  (Euler matrix)

void    Matrix44_Rotate_X     (real8 a[4][4], const real8 rx);                                                           // A = RX(rx)
void    Matrix44_Rotate_Y     (real8 a[4][4], const real8 ry);                                                           // A = RY(ry)
void    Matrix44_Rotate_Z     (real8 a[4][4], const real8 rz);                                                           // A = RZ(rz)
void    Matrix44_Rotate       (real8 a[4][4], const real8 rx, const real8 ry, const real8 rz);                           // A = Rot(rx,ry,rz)  (Euler matrix)

void    Matrix22_Diag         (real8 a[2][2], const real8 s );                                                           // A = Diag(s)
void    Matrix22_Diag         (real8 a[2][2], const real8 sx, const real8 sy);                                           // A = Diag(sx,sy,sz)

void    Matrix33_Diag         (real8 a[3][3], const real8 s );                                                           // A = Diag(s)
void    Matrix33_Diag         (real8 a[3][3], const real8 sx, const real8 sy, const real8 sz);                           // A = Diag(sx,sy,sz)

void    Matrix44_Diag         (real8 a[4][4], const real8 s );                                                           // A = Diag(s)
void    Matrix44_Diag         (real8 a[4][4], const real8 sx, const real8 sy, const real8 sz);                           // A = Diag(sx,sy,sz)

void    MatrixMN_Diag         (real8  *a    , const int m, const int n, const real8 *s);                                 // A = sets diagonal elements to s[i]
void    MatrixMN_Diag         (real8 **a    , const int m, const int n, const real8 *s);                                 // A = sets diagonal elements to s[i]

void    Matrix22_Add          (real8 a[2][2], real8 b[2][2], real8 c[2][2]);                                             // A = B+C
void    Matrix22_Sub          (real8 a[2][2], real8 b[2][2], real8 c[2][2]);                                             // A = B-C
void    Matrix22_Mul          (real8 a[2][2], real8 b[2][2], real8 c[2][2]);                                             // A = B*C   (A cannot be either B or C)
void    Matrix22_Mul          (real8 a[2][2], real8 b[2][2], real8 c[2][2], real8 d[2][2]);                              // A = B*C*D (A cannot be either B,C, or D)
void    Matrix22_Mul          (real8 a[2][2], real8 b[2][3], real8 c[3][2]);                                             // A = B*C   (inner product of 2x3 and 3x2 matrices)

void    Matrix33_Add          (real8 a[3][3], real8 b[3][3], real8 c[3][3]);                                             // A = B+C
void    Matrix33_Sub          (real8 a[3][3], real8 b[3][3], real8 c[3][3]);                                             // A = B-C
void    Matrix33_Mul          (real8 a[3][3], real8 b[3][3], real8 c[3][3]);                                             // A = B*C   (A cannot be either B or C)
void    Matrix33_Mul          (real8 a[3][3], real8 b[3][3], real8 c[3][3], real8 d[3][3]);                              // A = B*C*D (A cannot be either B,C, or D)

void    Matrix22_Vmul         (real8 y[2]   , real8 a[2][2], real8 x[2]   );                                             // Y = A*X  (2x2 matrix vector post multiply)
void    Matrix22_Vmul         (real8 y[2]   , real8 x[2]   , real8 a[2][2]);                                             // Y = X*A  (2x2 vector matrix pre  multiply)
void    Matrix23_Vmul         (real8 y[2]   , real8 a[2][3], real8 x[3]   );                                             // Y = A*X  (2x3 matrix vector post multiply of 3x1 vector)
void    Matrix32_Vmul         (real8 y[3]   , real8 a[3][2], real8 x[2]   );                                             // Y = A*X  (3x2 matrix vector post multiply of 2x1 vector)
void    Matrix33_Vmul         (real8 y[3]   , real8 a[3][3], real8 x[3]   );                                             // Y = A*X  (3x3 matrix vector post multiply)
void    Matrix33_Vmul         (real8 y[3]   , real8 x[3]   , real8 a[3][3]);                                             // Y = X*A  (3x3 vector matrix pre  multiply)
void    Matrix33_Vmul         (real8 a[3][3], real8 v[3]                  );                                             // A = V*Vt (3x3 outer product of 3x1 vector and its transpose)
void    Matrix33_Vmul         (real8 a[3][3], real8 v[3]   , real8 vt[3]  );                                             // A = V*Vt (3x3 outer product of 3x1 vector and 1x3 vector)
void    Matrix44_Vmul         (real8 y[3]   , real8 a[4][4], real8 x[3]   );                                             // Y = A*X  (4x4 matrix vector post multiply)
void    Matrix44_Vmul         (real8 y[3]   , real8 x[3]   , real8 a[4][4]);                                             // Y = X*A  (4x4 vector matrix pre  multiply)
void    Matrix43_Vmul         (real8 y[4]   , real8 a[4][3], real8 x[3]   );                                             // Y = A*X  (4x4 matrix vector post multiply of 3x1 vector, assumed 1)
void    Matrix63_Vmul         (real8 y[6]   , real8 a[6][3], real8 x[3]   );                                             // Y = A*X  (6x3 matrix vector post multiply of 3x1 vector)

void    Matrix33_VVt_Mul      (real8 a[3][3], real8 v[3]);                                                               // A = V*Vt (3x3 outer product of 3x1 vector and its transpose)
void    Matrix33_VVt_Mul      (real8 a[3][3], real8 v[3], real8 vt[3]);                                                  // A = V*Vt (3x3 outer product of 3x1 and 1x3 vectors)

void    Matrix44_Add          (real8 a[4][4], real8 b[4][4], real8 c[4][4]);                                             // A = B+C
void    Matrix44_Sub          (real8 a[4][4], real8 b[4][4], real8 c[4][4]);                                             // A = B-C
void    Matrix44_Mul          (real8 a[4][4], real8 b[4][4], real8 c[4][4]);                                             // A = B*C   (A cannot be either B or C)
void    Matrix44_Mul          (real8 a[4][4], real8 b[4][4], real8 c[4][4], real8 d[4][4]);                              // A = B*C*D (A cannot be either B,C, or D)

void    Matrix66_Add          (real8 a[6][6], real8 b[6][6], real8 c[6][6]);                                             // A = B+C
void    Matrix66_Sub          (real8 a[6][6], real8 b[6][6], real8 c[6][6]);                                             // A = B-C
void    Matrix66_Mul          (real8 a[6][6], real8 b[6][6], real8 c[6][6]);                                             // A = B*C  (A cannot be either B or C)

void    MatrixMN_Add          (real8  *a, const real8  *b, const real8  *c, const int m, const int n);                   // A = B+C
void    MatrixMN_Sub          (real8  *a, const real8  *b, const real8  *c, const int m, const int n);                   // A = B-C
void    MatrixMN_Mul          (real8  *a, const real8  *b, const real8  *c, const int m, const int n, const int p);      // A = B*C  (A cannot be either B or C)

void    MatrixMN_Add          (real8 **a,       real8 **b,       real8 **c, const int m, const int n);                   // A = B+C
void    MatrixMN_Sub          (real8 **a,       real8 **b,       real8 **c, const int m, const int n);                   // A = B-C
void    MatrixMN_Mul          (real8 **a,       real8 **b,       real8 **c, const int m, const int n, const int p);      // A = B*C  (A cannot be either B or C)

void    Matrix22_Transpose    (real8 a[2][2]);                                                                           // A = T(A)
void    Matrix33_Transpose    (real8 a[3][3]);                                                                           // A = T(A)
void    Matrix44_Transpose    (real8 a[4][4]);                                                                           // A = T(A)
void    Matrix66_Transpose    (real8 a[6][6]);                                                                           // A = T(A)
void    Matrix23_Transpose    (real8 a[3][2], real8 b[2][3] );                                                           // A(3x2) = T(B(2x3))
void    Matrix32_Transpose    (real8 a[2][3], real8 b[3][2] );                                                           // A(2x3) = T(B(3x2))
void    Matrix22_Transpose    (real8 a[2][2], real8 b[2][2]);                                                            // A = T(B)
void    Matrix33_Transpose    (real8 a[3][3], real8 b[3][3]);                                                            // A = T(B)
void    Matrix44_Transpose    (real8 a[4][4], real8 b[4][4]);                                                            // A = T(B)
void    Matrix66_Transpose    (real8 a[6][6], real8 b[6][6]);                                                            // A = T(B)
void    MatrixMN_Transpose    (real8  *a,            const int m, const int n);                                          // A = T(A)
void    MatrixMN_Transpose    (real8  *a, real8  *b, const int m, const int n);                                          // A = T(B)
void    MatrixMN_Transpose    (real8 **a, real8 **b, const int m, const int n);                                          // A = T(B)

int     Matrix22_Near         (const real8 a[2][2], const real8 s      , const real8 tol=1.0e-6);                        // (A==s) ? 1 : 0
int     Matrix22_Near         (const real8 a[2][2], const real8 b[2][2], const real8 tol=1.0e-6);                        // (A==B) ? 1 : 0
int     Matrix33_Near         (const real8 a[3][3], const real8 s      , const real8 tol=1.0e-6);                        // (A==s) ? 1 : 0
int     Matrix33_Near         (const real8 a[3][3], const real8 b[3][3], const real8 tol=1.0e-6);                        // (A==B) ? 1 : 0
int     Matrix44_Near         (const real8 a[4][4], const real8 s      , const real8 tol=1.0e-6);                        // (A==s) ? 1 : 0
int     Matrix44_Near         (const real8 a[4][4], const real8 b[4][4], const real8 tol=1.0e-6);                        // (A==B) ? 1 : 0
int     Matrix66_Near         (const real8 a[6][6], const real8 s      , const real8 tol=1.0e-6);                        // (A==s) ? 1 : 0
int     Matrix66_Near         (const real8 a[6][6], const real8 b[6][6], const real8 tol=1.0e-6);                        // (A==B) ? 1 : 0

real8   Matrix22_Det          (real8 a[2][2]);                                                                           // s = Det (A)
real8   Matrix22_Norm         (real8 a[2][2]);                                                                           // s = Norm(A)
real8   Matrix22_Condition    (real8 a[2][2]);                                                                           // s = Cond(A)
void    Matrix22_Adjoint      (real8 a[2][2], real8 b[2][2]);                                                            // A = Adj (B)
int     Matrix22_Inverse      (real8 a[2][2], real8 b[2][2]);                                                            // A = Inv (B)

real8   Matrix33_Det          (real8 a[3][3]);                                                                           // s = Det (A)
real8   Matrix33_Norm         (real8 a[3][3]);                                                                           // s = Norm(A)
real8   Matrix33_Condition    (real8 a[3][3]);                                                                           // s = Cond(A)
void    Matrix33_Adjoint      (real8 a[3][3], real8 b[3][3]);                                                            // A = Adj (B)
int     Matrix33_Inverse      (real8 a[3][3], real8 b[3][3]);                                                            // A = Inv (B)

real8   Matrix44_Det          (real8 a[4][4]);                                                                           // s = Det (A)
real8   Matrix44_Norm         (real8 a[4][4]);                                                                           // s = Norm(A)
real8   Matrix44_Condition    (real8 a[4][4]);                                                                           // s = Cond(A)
void    Matrix44_Adjoint      (real8 a[4][4], real8 b[4][4]);                                                            // A = Adj (B)
int     Matrix44_Inverse      (real8 a[4][4], real8 b[4][4]);                                                            // A = Inv (B)

int     Matrix66_Inverse      (real8 a[6][6], real8 b[6][6]);                                                            // A = Inv (B)

void    MatrixMN_Identity     (real8  *a, const int m , const int n );                                                   // A = I
void    MatrixMN_Identity     (real8 **a, const int m , const int n );                                                   // A = I
void    MatrixMN_Identity     (real8 **a, const int m0, const int m1, const int n0, const int n1 );                      // A = I


// (print utilities)

void    Matrix22_Print        (FILE  *fd, const char *fmt, real8  a[2][2] );                                                      // (print 2x2)
void    Matrix33_Print        (FILE  *fd, const char *fmt, real8  a[3][3] );                                                      // (print 3x3)
void    Matrix44_Print        (FILE  *fd, const char *fmt, real8  a[4][4] );                                                      // (print 4x4)
void    Matrix66_Print        (FILE  *fd, const char *fmt, real8  a[6][6] );                                                      // (print 6x6)

void    MatrixMN_Print        (FILE  *fd, const char *fmt, real8  *a, const int m , const int n  );                               // (print MxN, assumes zero-indexed)
void    MatrixMN_Print        (FILE  *fd, const char *fmt, real8 **a, const int m , const int n  );                               // (print MxN, assumes zero-indexed)
void    MatrixMN_Print        (FILE  *fd, const char *fmt, real8 **a, const int m0, const int m1, const int n0, const int n1 );   // (print MxN)

// (allocation/deallocation)

real8  *MatrixMN_Alloc        (const int m , const int n );                                                              // allocates and returns a dense, MxN matrix
real8  *MatrixMN_Alloc        (const int m0, const int m1, const int n0, const int n1);

real8 **MatrixMN_Iliffe       (                const int m , const int n );                                              // allocates and returns a  zero-indexed  MxN iliffe matrix,
real8 **MatrixMN_Iliffe       (                const int m0, const int m1, const int n0, const int n1);                  // allocates and returns an      indexed, MxN iliffe matrix,
real8 **MatrixMN_Iliffe       (const real8 *a, const int m , const int n );                                              // allocates and initializes an MxN iliffe vector pointing to and existing matrix (zero indexed)
real8 **MatrixMN_Iliffe       (const real8 *a, const int m0, const int m1, const int n0, const int n1);                  // allocates and initializes an MxN iliffe vector pointing to and existing matrix (     indexed)

void    MatrixMN_Free         (real8  *a);                                                                               // releases dense, zero-indexed matrix
void    MatrixMN_Free         (real8 **a);                                                                               // releases iliffe vector AND implied dense, zero-indexed matrix (pointed to by a[0]
void    MatrixMN_Free         (real8 **a, real8 *am);                                                                    // releases iliffe vector and named   dense, zero-indexed matrix
void    MatrixMN_Free         (real8 **a,            const int m0, const int m1, const int n0, const int n1);            // releases iliffe vector AND implied dense,      indexed matrix (pointed to by a[m0]
void    MatrixMN_Free         (real8 **a, real8 *am, const int m0, const int m1, const int n0, const int n1);            // releases iliffe vector and named   dense,      indexed matrix

// (comparison)

int     MatrixMN_Near         (const real8  *a, const real8   s, const int m , const int n ,                              const real8 tol=1.0e-6 );  // (A==s) ? 1 : 0
int     MatrixMN_Near         (      real8 **a, const real8   s, const int m , const int n ,                              const real8 tol=1.0e-6 );  // (A==s) ? 1 : 0
int     MatrixMN_Near         (      real8 **a, const real8   s, const int m0, const int m1, const int n0,  const int n1, const real8 tol=1.0e-6 );  // (A==s) ? 1 : 0

int     MatrixMN_Near         (const real8  *a, const real8  *b, const int m , const int n ,                              const real8 tol=1.0e-6 );  // (A==B) ? 1 : 0
int     MatrixMN_Near         (      real8 **a,       real8 **b, const int m , const int n ,                              const real8 tol=1.0e-6 );  // (A==B) ? 1 : 0
int     MatrixMN_Near         (      real8 **a,       real8 **b, const int m0, const int m1, const int n0,  const int n1, const real8 tol=1.0e-6 );  // (A==B) ? 1 : 0

// (statistical operations)

real8   MatrixMN_Min          (const real8 *a, const int m, const int n);                                                //  s = Min(A)
real8   MatrixMN_Max          (const real8 *a, const int m, const int n);                                                //  s = Max(A)
real8   MatrixMN_Mean         (const real8 *a, const int m, const int n);                                                //  s = Mean(A)
real8   MatrixMN_Range        (const real8 *a, const int m, const int n);                                                //  s = Max(A)-Min(A)

real8   MatrixMN_Min          (int & iy, int & ix, const real8 *a, const int m, const int n);                            //  s = Min(A), (iy,ix) identifies location (row,col)
real8   MatrixMN_Max          (int & iy, int & ix, const real8 *a, const int m, const int n);                            //  s = Max(A), (iy,ix) identifies location (row,col)

void    MatrixMN_Moment       (real8 & min, real8 & max, real8 & mean, real8 & range              , const real8 *a, const int m, const int n);  // (min,max,mean,range     )(A)  (statistical moment)
void    MatrixMN_Moment       (real8 & min, real8 & max, real8 & mean, real8 & range, real8 & sdev, const real8 *a, const int m, const int n);  // (min,max,mean,range,sdev)(A)  (statistical moment)

//---------------------------------------------------------------------------------------------------------

void  MatrixMult              (real8    *a, int aRows, int aCols, int aLD,
                               real8    *b, int bCols,            int bLD,
                               real8    *c,                       int cLD );

#ifdef ANISOTROPIC
void MatrixDblMultComplex     (real8    *a, int aRows, int aCols, int aLD,
                               complex8 *b, int bCols,            int bLD,
                               complex8 *c,                       int cLD );
#endif

void  MatrixMultArb           (real8 *a, int aCols, int aRowOffset, int aColOffset, int aMultRows, int aMultCols,
                               real8 *b, int bCols, int bRowOffset, int bColOffset, int bMultCols,
                               real8 *c, int cCols, int cRowOffset, int cColOffset );

int   MatrixInvert            (real8 *mat, real8 *invMat, int order, int lda );

void  Matrix33Vector3Multiply (real8 a[3][3], real8 x[3]   , real8 y[3]    );

void  Vec3TransposeAndMult    (real8 v[3]   , real8 result[3][3] );

// Pseudo-inverse support...
//---------------------------------------------------------------------------------------------------------

int   MatrixMN_SVD            (real8 **a, real8 **u, real8 *w, real8 *v, const int m, const int n);      // SVD(A) = U*W*V'
int   MatrixMN_SVD            (real8  *a, real8  *u, real8 *w, real8 *v, const int m, const int n);      // SVD(A) = U*W*V'

int   MatrixMN_SVD_Inverse    (real8 *ai, const real8 *a, const int m, const int n);                     // AI = SVD Inverse(A)
int   MatrixMN_SVD_Inverse    (                 real8 *a, const int m, const int n);                     // A  = SVD Inverse(A) (in-place)

int   Matrix33_SVD            (real8 a[3][3] , real8 u [3][3], real8  w[3], real8 v[3][3] );             // SVD(A) = U*W*V'
int   Matrix33_SVD_Inverse    (real8 a[3][3] );                                                          // A  = SVD Inverse(A) (in-place)
int   Matrix33_SVD_Inverse    (real8 a[3][3] , real8 b[3][3] );                                          // A  = SVD Inverse(B)

int   Matrix33_PseudoInverse  (real8 a[3][3] , real8 b[3][3] );                                          // A  = PseudoInverse(B) (Nicolas Bertin version)

real8 Matrix33_SVD_Condition  (real8 a[3][3] );                                                          // returns SVD condition

#endif
