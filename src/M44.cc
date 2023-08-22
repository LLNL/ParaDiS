/**************************************************************************
 *  Module:  M44 - implements a simple 3x3 transform matrix
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "V3.h"
#include "M44.h"

//-----------------------------------------------------------------------------------------------

real8 Matrix_4x4::Determinant (void)
{
   real8 mbt[16]; // array of transpose source matrix
   real8 tmp[16]; // array for cofactor pairs
   real8 cfx[16]; // array of cofactors

   M44_TRANSPOSE(mbt,mtx);

   // calculate pairs for first 8 cofactors

   tmp[ 0] = (mbt[10]*mbt[15]); 
   tmp[ 1] = (mbt[11]*mbt[14]); 
   tmp[ 2] = (mbt[ 9]*mbt[15]); 
   tmp[ 3] = (mbt[11]*mbt[13]); 
   tmp[ 4] = (mbt[ 9]*mbt[14]); 
   tmp[ 5] = (mbt[10]*mbt[13]); 
   tmp[ 6] = (mbt[ 8]*mbt[15]); 
   tmp[ 7] = (mbt[11]*mbt[12]); 
   tmp[ 8] = (mbt[ 8]*mbt[14]); 
   tmp[ 9] = (mbt[10]*mbt[12]); 
   tmp[10] = (mbt[ 8]*mbt[13]); 
   tmp[11] = (mbt[ 9]*mbt[12]); 

   // calculate first 8 cofactors

   cfx[ 0]  = (tmp[ 0]*mbt[ 5]) + (tmp[ 3]*mbt[ 6]) + (tmp[ 4]*mbt[ 7]); 
   cfx[ 0] -= (tmp[ 1]*mbt[ 5]) + (tmp[ 2]*mbt[ 6]) + (tmp[ 5]*mbt[ 7]); 
   cfx[ 1]  = (tmp[ 1]*mbt[ 4]) + (tmp[ 6]*mbt[ 6]) + (tmp[ 9]*mbt[ 7]); 
   cfx[ 1] -= (tmp[ 0]*mbt[ 4]) + (tmp[ 7]*mbt[ 6]) + (tmp[ 8]*mbt[ 7]); 
   cfx[ 2]  = (tmp[ 2]*mbt[ 4]) + (tmp[ 7]*mbt[ 5]) + (tmp[10]*mbt[ 7]); 
   cfx[ 2] -= (tmp[ 3]*mbt[ 4]) + (tmp[ 6]*mbt[ 5]) + (tmp[11]*mbt[ 7]); 
   cfx[ 3]  = (tmp[ 5]*mbt[ 4]) + (tmp[ 8]*mbt[ 5]) + (tmp[11]*mbt[ 6]); 
   cfx[ 3] -= (tmp[ 4]*mbt[ 4]) + (tmp[ 9]*mbt[ 5]) + (tmp[10]*mbt[ 6]); 
   cfx[ 4]  = (tmp[ 1]*mbt[ 1]) + (tmp[ 2]*mbt[ 2]) + (tmp[ 5]*mbt[ 3]); 
   cfx[ 4] -= (tmp[ 0]*mbt[ 1]) + (tmp[ 3]*mbt[ 2]) + (tmp[ 4]*mbt[ 3]); 
   cfx[ 5]  = (tmp[ 0]*mbt[ 0]) + (tmp[ 7]*mbt[ 2]) + (tmp[ 8]*mbt[ 3]); 
   cfx[ 5] -= (tmp[ 1]*mbt[ 0]) + (tmp[ 6]*mbt[ 2]) + (tmp[ 9]*mbt[ 3]); 
   cfx[ 6]  = (tmp[ 3]*mbt[ 0]) + (tmp[ 6]*mbt[ 1]) + (tmp[11]*mbt[ 3]); 
   cfx[ 6] -= (tmp[ 2]*mbt[ 0]) + (tmp[ 7]*mbt[ 1]) + (tmp[10]*mbt[ 3]); 
   cfx[ 7]  = (tmp[ 4]*mbt[ 0]) + (tmp[ 9]*mbt[ 1]) + (tmp[10]*mbt[ 2]); 
   cfx[ 7] -= (tmp[ 5]*mbt[ 0]) + (tmp[ 8]*mbt[ 1]) + (tmp[11]*mbt[ 2]); 

   // calculate pairs for second 8 cofactors

   tmp[ 0] = (mbt[ 2]*mbt[ 7]); 
   tmp[ 1] = (mbt[ 3]*mbt[ 6]); 
   tmp[ 2] = (mbt[ 1]*mbt[ 7]); 
   tmp[ 3] = (mbt[ 3]*mbt[ 5]); 
   tmp[ 4] = (mbt[ 1]*mbt[ 6]); 
   tmp[ 5] = (mbt[ 2]*mbt[ 5]);
   tmp[ 6] = (mbt[ 0]*mbt[ 7]); 
   tmp[ 7] = (mbt[ 3]*mbt[ 4]); 
   tmp[ 8] = (mbt[ 0]*mbt[ 6]); 
   tmp[ 9] = (mbt[ 2]*mbt[ 4]); 
   tmp[10] = (mbt[ 0]*mbt[ 5]); 
   tmp[11] = (mbt[ 1]*mbt[ 4]); 

   // calculate second 8 cofactors

   cfx[ 8]  = (tmp[ 0]*mbt[13]) + (tmp[ 3]*mbt[14]) + (tmp[ 4]*mbt[15]); 
   cfx[ 8] -= (tmp[ 1]*mbt[13]) + (tmp[ 2]*mbt[14]) + (tmp[ 5]*mbt[15]); 
   cfx[ 9]  = (tmp[ 1]*mbt[12]) + (tmp[ 6]*mbt[14]) + (tmp[ 9]*mbt[15]); 
   cfx[ 9] -= (tmp[ 0]*mbt[12]) + (tmp[ 7]*mbt[14]) + (tmp[ 8]*mbt[15]); 
   cfx[10]  = (tmp[ 2]*mbt[12]) + (tmp[ 7]*mbt[13]) + (tmp[10]*mbt[15]); 
   cfx[10] -= (tmp[ 3]*mbt[12]) + (tmp[ 6]*mbt[13]) + (tmp[11]*mbt[15]); 
   cfx[11]  = (tmp[ 5]*mbt[12]) + (tmp[ 8]*mbt[13]) + (tmp[11]*mbt[14]); 
   cfx[11] -= (tmp[ 4]*mbt[12]) + (tmp[ 9]*mbt[13]) + (tmp[10]*mbt[14]); 
   cfx[12]  = (tmp[ 2]*mbt[10]) + (tmp[ 5]*mbt[11]) + (tmp[ 1]*mbt[ 9]); 
   cfx[12] -= (tmp[ 4]*mbt[11]) + (tmp[ 0]*mbt[ 9]) + (tmp[ 3]*mbt[10]); 
   cfx[13]  = (tmp[ 8]*mbt[11]) + (tmp[ 0]*mbt[ 8]) + (tmp[ 7]*mbt[10]); 
   cfx[13] -= (tmp[ 6]*mbt[10]) + (tmp[ 9]*mbt[11]) + (tmp[ 1]*mbt[ 8]); 
   cfx[14]  = (tmp[ 6]*mbt[ 9]) + (tmp[11]*mbt[11]) + (tmp[ 3]*mbt[ 8]); 
   cfx[14] -= (tmp[10]*mbt[11]) + (tmp[ 2]*mbt[ 8]) + (tmp[ 7]*mbt[ 9]); 
   cfx[15]  = (tmp[10]*mbt[10]) + (tmp[ 4]*mbt[ 8]) + (tmp[ 9]*mbt[ 9]); 
   cfx[15] -= (tmp[ 8]*mbt[ 9]) + (tmp[11]*mbt[10]) + (tmp[ 5]*mbt[ 8]); 

   // calculate determinant

   real8 det = (mbt[ 0]*cfx[ 0]) + (mbt[ 1]*cfx[ 1]) + (mbt[ 2]*cfx[ 2]) + (mbt[ 3]*cfx[ 3]); 
         det = ( (det!=0.0) ? (1.0/det) : det ); 

   return(det);
}

//-----------------------------------------------------------------------------------------------

void Matrix_4x4::Inverse (void)
{
   Inverse(*this);
}

void Matrix_4x4::Inverse (const Matrix_4x4 & m)
{
   real8 tmp[16]; // temp array for pairs
   real8 mbt[16]; // array of transpose source matrix
   real8 cfx[16]; // array of cofactors

   M44_TRANSPOSE(mbt,m.mtx);

   // calculate pairs for first 8 cofactors

   tmp[ 0] = (mbt[10]*mbt[15]); 
   tmp[ 1] = (mbt[11]*mbt[14]); 
   tmp[ 2] = (mbt[ 9]*mbt[15]); 
   tmp[ 3] = (mbt[11]*mbt[13]); 
   tmp[ 4] = (mbt[ 9]*mbt[14]); 
   tmp[ 5] = (mbt[10]*mbt[13]); 
   tmp[ 6] = (mbt[ 8]*mbt[15]); 
   tmp[ 7] = (mbt[11]*mbt[12]); 
   tmp[ 8] = (mbt[ 8]*mbt[14]); 
   tmp[ 9] = (mbt[10]*mbt[12]); 
   tmp[10] = (mbt[ 8]*mbt[13]); 
   tmp[11] = (mbt[ 9]*mbt[12]); 

   // calculate first 8 cofactors

   cfx[ 0]  = (tmp[ 0]*mbt[ 5]) + (tmp[ 3]*mbt[ 6]) + (tmp[ 4]*mbt[ 7]); 
   cfx[ 0] -= (tmp[ 1]*mbt[ 5]) + (tmp[ 2]*mbt[ 6]) + (tmp[ 5]*mbt[ 7]); 
   cfx[ 1]  = (tmp[ 1]*mbt[ 4]) + (tmp[ 6]*mbt[ 6]) + (tmp[ 9]*mbt[ 7]); 
   cfx[ 1] -= (tmp[ 0]*mbt[ 4]) + (tmp[ 7]*mbt[ 6]) + (tmp[ 8]*mbt[ 7]); 
   cfx[ 2]  = (tmp[ 2]*mbt[ 4]) + (tmp[ 7]*mbt[ 5]) + (tmp[10]*mbt[ 7]); 
   cfx[ 2] -= (tmp[ 3]*mbt[ 4]) + (tmp[ 6]*mbt[ 5]) + (tmp[11]*mbt[ 7]); 
   cfx[ 3]  = (tmp[ 5]*mbt[ 4]) + (tmp[ 8]*mbt[ 5]) + (tmp[11]*mbt[ 6]); 
   cfx[ 3] -= (tmp[ 4]*mbt[ 4]) + (tmp[ 9]*mbt[ 5]) + (tmp[10]*mbt[ 6]); 
   cfx[ 4]  = (tmp[ 1]*mbt[ 1]) + (tmp[ 2]*mbt[ 2]) + (tmp[ 5]*mbt[ 3]); 
   cfx[ 4] -= (tmp[ 0]*mbt[ 1]) + (tmp[ 3]*mbt[ 2]) + (tmp[ 4]*mbt[ 3]); 
   cfx[ 5]  = (tmp[ 0]*mbt[ 0]) + (tmp[ 7]*mbt[ 2]) + (tmp[ 8]*mbt[ 3]); 
   cfx[ 5] -= (tmp[ 1]*mbt[ 0]) + (tmp[ 6]*mbt[ 2]) + (tmp[ 9]*mbt[ 3]); 
   cfx[ 6]  = (tmp[ 3]*mbt[ 0]) + (tmp[ 6]*mbt[ 1]) + (tmp[11]*mbt[ 3]); 
   cfx[ 6] -= (tmp[ 2]*mbt[ 0]) + (tmp[ 7]*mbt[ 1]) + (tmp[10]*mbt[ 3]); 
   cfx[ 7]  = (tmp[ 4]*mbt[ 0]) + (tmp[ 9]*mbt[ 1]) + (tmp[10]*mbt[ 2]); 
   cfx[ 7] -= (tmp[ 5]*mbt[ 0]) + (tmp[ 8]*mbt[ 1]) + (tmp[11]*mbt[ 2]); 

   // calculate pairs for second 8 cofactors

   tmp[ 0] = (mbt[ 2]*mbt[ 7]); 
   tmp[ 1] = (mbt[ 3]*mbt[ 6]); 
   tmp[ 2] = (mbt[ 1]*mbt[ 7]); 
   tmp[ 3] = (mbt[ 3]*mbt[ 5]); 
   tmp[ 4] = (mbt[ 1]*mbt[ 6]); 
   tmp[ 5] = (mbt[ 2]*mbt[ 5]);
   tmp[ 6] = (mbt[ 0]*mbt[ 7]); 
   tmp[ 7] = (mbt[ 3]*mbt[ 4]); 
   tmp[ 8] = (mbt[ 0]*mbt[ 6]); 
   tmp[ 9] = (mbt[ 2]*mbt[ 4]); 
   tmp[10] = (mbt[ 0]*mbt[ 5]); 
   tmp[11] = (mbt[ 1]*mbt[ 4]); 

   // calculate second 8 cofactors

   cfx[ 8]  = (tmp[ 0]*mbt[13]) + (tmp[ 3]*mbt[14]) + (tmp[ 4]*mbt[15]); 
   cfx[ 8] -= (tmp[ 1]*mbt[13]) + (tmp[ 2]*mbt[14]) + (tmp[ 5]*mbt[15]); 
   cfx[ 9]  = (tmp[ 1]*mbt[12]) + (tmp[ 6]*mbt[14]) + (tmp[ 9]*mbt[15]); 
   cfx[ 9] -= (tmp[ 0]*mbt[12]) + (tmp[ 7]*mbt[14]) + (tmp[ 8]*mbt[15]); 
   cfx[10]  = (tmp[ 2]*mbt[12]) + (tmp[ 7]*mbt[13]) + (tmp[10]*mbt[15]); 
   cfx[10] -= (tmp[ 3]*mbt[12]) + (tmp[ 6]*mbt[13]) + (tmp[11]*mbt[15]); 
   cfx[11]  = (tmp[ 5]*mbt[12]) + (tmp[ 8]*mbt[13]) + (tmp[11]*mbt[14]); 
   cfx[11] -= (tmp[ 4]*mbt[12]) + (tmp[ 9]*mbt[13]) + (tmp[10]*mbt[14]); 
   cfx[12]  = (tmp[ 2]*mbt[10]) + (tmp[ 5]*mbt[11]) + (tmp[ 1]*mbt[ 9]); 
   cfx[12] -= (tmp[ 4]*mbt[11]) + (tmp[ 0]*mbt[ 9]) + (tmp[ 3]*mbt[10]); 
   cfx[13]  = (tmp[ 8]*mbt[11]) + (tmp[ 0]*mbt[ 8]) + (tmp[ 7]*mbt[10]); 
   cfx[13] -= (tmp[ 6]*mbt[10]) + (tmp[ 9]*mbt[11]) + (tmp[ 1]*mbt[ 8]); 
   cfx[14]  = (tmp[ 6]*mbt[ 9]) + (tmp[11]*mbt[11]) + (tmp[ 3]*mbt[ 8]); 
   cfx[14] -= (tmp[10]*mbt[11]) + (tmp[ 2]*mbt[ 8]) + (tmp[ 7]*mbt[ 9]); 
   cfx[15]  = (tmp[10]*mbt[10]) + (tmp[ 4]*mbt[ 8]) + (tmp[ 9]*mbt[ 9]); 
   cfx[15] -= (tmp[ 8]*mbt[ 9]) + (tmp[11]*mbt[10]) + (tmp[ 5]*mbt[ 8]); 

   // calculate determinant

   real8 det = (mbt[ 0]*cfx[ 0]) + (mbt[ 1]*cfx[ 1]) + (mbt[ 2]*cfx[ 2]) + (mbt[ 3]*cfx[ 3]); 
         det = ( (det!=0.0) ? (1.0/det) : det ); 

   // calculate matrix inverse

   mtx[ 0] = (cfx[ 0]*det);
   mtx[ 1] = (cfx[ 1]*det);
   mtx[ 2] = (cfx[ 2]*det);
   mtx[ 3] = (cfx[ 3]*det);
   mtx[ 4] = (cfx[ 4]*det);
   mtx[ 5] = (cfx[ 5]*det);
   mtx[ 6] = (cfx[ 6]*det);
   mtx[ 7] = (cfx[ 7]*det);
   mtx[ 8] = (cfx[ 8]*det);
   mtx[ 9] = (cfx[ 9]*det);
   mtx[10] = (cfx[10]*det);
   mtx[11] = (cfx[11]*det);
   mtx[12] = (cfx[12]*det);
   mtx[13] = (cfx[13]*det);
   mtx[14] = (cfx[14]*det);
   mtx[15] = (cfx[15]*det);
}

//-----------------------------------------------------------------------------------------------

void Matrix_4x4::Transform_M44_V3N (real8 *v, const int vn) const 
{
   if (v && (vn>0)) for (int i=0; (i<vn); i++, v+=3) { V3_M44_V3_MUL(v,mtx,v); }
}

void Matrix_4x4::Transform_V3N_M44 (real8 *v, const int vn) const
{
   if (v && (vn>0)) for (int i=0; (i<vn); i++, v+=3) { V3_V3_M44_MUL(v,v,mtx); }
}

void Matrix_4x4::Transform_M44_V3N (real8 *va, real8 *vb, const int vn) const
{
   if (va && vb && (vn>0)) for (int i=0; (i<vn); i++, va+=3,vb+=3) { V3_M44_V3_MUL(va,mtx,vb); }
}

void Matrix_4x4::Transform_V3N_M44 (real8 *va, real8 *vb, const int vn) const
{
   if (va && vb && (vn>0)) for (int i=0; (i<vn); i++, va+=3,vb+=3) { V3_V3_M44_MUL(va,vb,mtx); }
}

/*===========================================================================================*/

void M44_INV(real8 *ma, const real8 *mb)
{
   real8 tmp[16]; // temp array for pairs
   real8 mbt[16]; // array of transpose source matrix
   real8 dst[16]; // 

   M44_TRANSPOSE(mbt,mb);

   // calculate pairs for first 8 cofactors

   tmp[ 0] = (mbt[10]*mbt[15]); 
   tmp[ 1] = (mbt[11]*mbt[14]); 
   tmp[ 2] = (mbt[ 9]*mbt[15]); 
   tmp[ 3] = (mbt[11]*mbt[13]); 
   tmp[ 4] = (mbt[ 9]*mbt[14]); 
   tmp[ 5] = (mbt[10]*mbt[13]); 
   tmp[ 6] = (mbt[ 8]*mbt[15]); 
   tmp[ 7] = (mbt[11]*mbt[12]); 
   tmp[ 8] = (mbt[ 8]*mbt[14]); 
   tmp[ 9] = (mbt[10]*mbt[12]); 
   tmp[10] = (mbt[ 8]*mbt[13]); 
   tmp[11] = (mbt[ 9]*mbt[12]); 

   // calculate first 8 cofactors

   dst[ 0]  = (tmp[ 0]*mbt[ 5]) + (tmp[ 3]*mbt[ 6]) + (tmp[ 4]*mbt[ 7]); 
   dst[ 0] -= (tmp[ 1]*mbt[ 5]) + (tmp[ 2]*mbt[ 6]) + (tmp[ 5]*mbt[ 7]); 
   dst[ 1]  = (tmp[ 1]*mbt[ 4]) + (tmp[ 6]*mbt[ 6]) + (tmp[ 9]*mbt[ 7]); 
   dst[ 1] -= (tmp[ 0]*mbt[ 4]) + (tmp[ 7]*mbt[ 6]) + (tmp[ 8]*mbt[ 7]); 
   dst[ 2]  = (tmp[ 2]*mbt[ 4]) + (tmp[ 7]*mbt[ 5]) + (tmp[10]*mbt[ 7]); 
   dst[ 2] -= (tmp[ 3]*mbt[ 4]) + (tmp[ 6]*mbt[ 5]) + (tmp[11]*mbt[ 7]); 
   dst[ 3]  = (tmp[ 5]*mbt[ 4]) + (tmp[ 8]*mbt[ 5]) + (tmp[11]*mbt[ 6]); 
   dst[ 3] -= (tmp[ 4]*mbt[ 4]) + (tmp[ 9]*mbt[ 5]) + (tmp[10]*mbt[ 6]); 
   dst[ 4]  = (tmp[ 1]*mbt[ 1]) + (tmp[ 2]*mbt[ 2]) + (tmp[ 5]*mbt[ 3]); 
   dst[ 4] -= (tmp[ 0]*mbt[ 1]) + (tmp[ 3]*mbt[ 2]) + (tmp[ 4]*mbt[ 3]); 
   dst[ 5]  = (tmp[ 0]*mbt[ 0]) + (tmp[ 7]*mbt[ 2]) + (tmp[ 8]*mbt[ 3]); 
   dst[ 5] -= (tmp[ 1]*mbt[ 0]) + (tmp[ 6]*mbt[ 2]) + (tmp[ 9]*mbt[ 3]); 
   dst[ 6]  = (tmp[ 3]*mbt[ 0]) + (tmp[ 6]*mbt[ 1]) + (tmp[11]*mbt[ 3]); 
   dst[ 6] -= (tmp[ 2]*mbt[ 0]) + (tmp[ 7]*mbt[ 1]) + (tmp[10]*mbt[ 3]); 
   dst[ 7]  = (tmp[ 4]*mbt[ 0]) + (tmp[ 9]*mbt[ 1]) + (tmp[10]*mbt[ 2]); 
   dst[ 7] -= (tmp[ 5]*mbt[ 0]) + (tmp[ 8]*mbt[ 1]) + (tmp[11]*mbt[ 2]); 

   // calculate pairs for second 8 cofactors

   tmp[ 0] = (mbt[ 2]*mbt[ 7]); 
   tmp[ 1] = (mbt[ 3]*mbt[ 6]); 
   tmp[ 2] = (mbt[ 1]*mbt[ 7]); 
   tmp[ 3] = (mbt[ 3]*mbt[ 5]); 
   tmp[ 4] = (mbt[ 1]*mbt[ 6]); 
   tmp[ 5] = (mbt[ 2]*mbt[ 5]);
   tmp[ 6] = (mbt[ 0]*mbt[ 7]); 
   tmp[ 7] = (mbt[ 3]*mbt[ 4]); 
   tmp[ 8] = (mbt[ 0]*mbt[ 6]); 
   tmp[ 9] = (mbt[ 2]*mbt[ 4]); 
   tmp[10] = (mbt[ 0]*mbt[ 5]); 
   tmp[11] = (mbt[ 1]*mbt[ 4]); 

   // calculate second 8 cofactors

   dst[ 8]  = (tmp[ 0]*mbt[13]) + (tmp[ 3]*mbt[14]) + (tmp[ 4]*mbt[15]); 
   dst[ 8] -= (tmp[ 1]*mbt[13]) + (tmp[ 2]*mbt[14]) + (tmp[ 5]*mbt[15]); 
   dst[ 9]  = (tmp[ 1]*mbt[12]) + (tmp[ 6]*mbt[14]) + (tmp[ 9]*mbt[15]); 
   dst[ 9] -= (tmp[ 0]*mbt[12]) + (tmp[ 7]*mbt[14]) + (tmp[ 8]*mbt[15]); 
   dst[10]  = (tmp[ 2]*mbt[12]) + (tmp[ 7]*mbt[13]) + (tmp[10]*mbt[15]); 
   dst[10] -= (tmp[ 3]*mbt[12]) + (tmp[ 6]*mbt[13]) + (tmp[11]*mbt[15]); 
   dst[11]  = (tmp[ 5]*mbt[12]) + (tmp[ 8]*mbt[13]) + (tmp[11]*mbt[14]); 
   dst[11] -= (tmp[ 4]*mbt[12]) + (tmp[ 9]*mbt[13]) + (tmp[10]*mbt[14]); 
   dst[12]  = (tmp[ 2]*mbt[10]) + (tmp[ 5]*mbt[11]) + (tmp[ 1]*mbt[ 9]); 
   dst[12] -= (tmp[ 4]*mbt[11]) + (tmp[ 0]*mbt[ 9]) + (tmp[ 3]*mbt[10]); 
   dst[13]  = (tmp[ 8]*mbt[11]) + (tmp[ 0]*mbt[ 8]) + (tmp[ 7]*mbt[10]); 
   dst[13] -= (tmp[ 6]*mbt[10]) + (tmp[ 9]*mbt[11]) + (tmp[ 1]*mbt[ 8]); 
   dst[14]  = (tmp[ 6]*mbt[ 9]) + (tmp[11]*mbt[11]) + (tmp[ 3]*mbt[ 8]); 
   dst[14] -= (tmp[10]*mbt[11]) + (tmp[ 2]*mbt[ 8]) + (tmp[ 7]*mbt[ 9]); 
   dst[15]  = (tmp[10]*mbt[10]) + (tmp[ 4]*mbt[ 8]) + (tmp[ 9]*mbt[ 9]); 
   dst[15] -= (tmp[ 8]*mbt[ 9]) + (tmp[11]*mbt[10]) + (tmp[ 5]*mbt[ 8]); 

   // calculate determinant

   real8 det = (mbt[ 0]*dst[ 0]) + (mbt[ 1]*dst[ 1]) + (mbt[ 2]*dst[ 2]) + (mbt[ 3]*dst[ 3]); 
         det = ( (det!=0.0) ? (1.0/det) : det ); 

   // calculate matrix inverse

   ma[ 0] = (dst[ 0]*det);
   ma[ 1] = (dst[ 1]*det);
   ma[ 2] = (dst[ 2]*det);
   ma[ 3] = (dst[ 3]*det);
   ma[ 4] = (dst[ 4]*det);
   ma[ 5] = (dst[ 5]*det);
   ma[ 6] = (dst[ 6]*det);
   ma[ 7] = (dst[ 7]*det);
   ma[ 8] = (dst[ 8]*det);
   ma[ 9] = (dst[ 9]*det);
   ma[10] = (dst[10]*det);
   ma[11] = (dst[11]*det);
   ma[12] = (dst[12]*det);
   ma[13] = (dst[13]*det);
   ma[14] = (dst[14]*det);
   ma[15] = (dst[15]*det);
}

int M44_Near (real8 *ma, const real8 *mb)
{
   return(M44_Near(ma,mb,1e-12));
}

int M44_Near (real8 *ma, const real8 *mb, const real8 eps)
{
   if (ma && mb)
   {
      if (fabs(ma[ 0]-mb[ 0])>eps) return(0);
      if (fabs(ma[ 1]-mb[ 1])>eps) return(0);
      if (fabs(ma[ 2]-mb[ 2])>eps) return(0);
      if (fabs(ma[ 3]-mb[ 3])>eps) return(0);
      if (fabs(ma[ 4]-mb[ 4])>eps) return(0);
      if (fabs(ma[ 5]-mb[ 5])>eps) return(0);
      if (fabs(ma[ 6]-mb[ 6])>eps) return(0);
      if (fabs(ma[ 7]-mb[ 7])>eps) return(0);
      if (fabs(ma[ 8]-mb[ 8])>eps) return(0);
      if (fabs(ma[ 9]-mb[ 9])>eps) return(0);
      if (fabs(ma[10]-mb[10])>eps) return(0);
      if (fabs(ma[11]-mb[11])>eps) return(0);
      if (fabs(ma[12]-mb[12])>eps) return(0);
      if (fabs(ma[13]-mb[13])>eps) return(0);
      if (fabs(ma[14]-mb[14])>eps) return(0);
      if (fabs(ma[15]-mb[15])>eps) return(0);

      return(1);
   }

   return(0);
}

int M44_Near (real8 *m, const real8 b, const real8 eps )
{
   if (m)
   {
      if (fabs(m[ 0]-b)>eps) return(0);
      if (fabs(m[ 1]-b)>eps) return(0);
      if (fabs(m[ 2]-b)>eps) return(0);
      if (fabs(m[ 3]-b)>eps) return(0);
      if (fabs(m[ 4]-b)>eps) return(0);
      if (fabs(m[ 5]-b)>eps) return(0);
      if (fabs(m[ 6]-b)>eps) return(0);
      if (fabs(m[ 7]-b)>eps) return(0);
      if (fabs(m[ 8]-b)>eps) return(0);
      if (fabs(m[ 9]-b)>eps) return(0);
      if (fabs(m[10]-b)>eps) return(0);
      if (fabs(m[11]-b)>eps) return(0);
      if (fabs(m[12]-b)>eps) return(0);
      if (fabs(m[13]-b)>eps) return(0);
      if (fabs(m[14]-b)>eps) return(0);
      if (fabs(m[15]-b)>eps) return(0);

      return(1);
   }

   return(0);
}

int M44_Near (
   real8 *m, 
   const real8 m00, const real8 m01, const real8 m02, const real8 m03,
   const real8 m10, const real8 m11, const real8 m12, const real8 m13,
   const real8 m20, const real8 m21, const real8 m22, const real8 m23,
   const real8 m30, const real8 m31, const real8 m32, const real8 m33, const real8 eps )
{
   if (m)
   {
      if (fabs(m[ 0]-m00)>eps) return(0);
      if (fabs(m[ 1]-m01)>eps) return(0);
      if (fabs(m[ 2]-m02)>eps) return(0);
      if (fabs(m[ 3]-m03)>eps) return(0);
      if (fabs(m[ 4]-m10)>eps) return(0);
      if (fabs(m[ 5]-m11)>eps) return(0);
      if (fabs(m[ 6]-m12)>eps) return(0);
      if (fabs(m[ 7]-m13)>eps) return(0);
      if (fabs(m[ 8]-m20)>eps) return(0);
      if (fabs(m[ 9]-m21)>eps) return(0);
      if (fabs(m[10]-m22)>eps) return(0);
      if (fabs(m[11]-m23)>eps) return(0);
      if (fabs(m[12]-m30)>eps) return(0);
      if (fabs(m[13]-m31)>eps) return(0);
      if (fabs(m[14]-m32)>eps) return(0);
      if (fabs(m[15]-m33)>eps) return(0);

      return(1);
   }

   return(0);
}

int M44_Near (
   real8 m[4][4], 
   const real8 m00, const real8 m01, const real8 m02, const real8 m03,
   const real8 m10, const real8 m11, const real8 m12, const real8 m13,
   const real8 m20, const real8 m21, const real8 m22, const real8 m23,
   const real8 m30, const real8 m31, const real8 m32, const real8 m33, const real8 eps )
{
   if (m)
   {
      if (fabs(m[0][0]-m00)>eps) return(0);
      if (fabs(m[0][1]-m01)>eps) return(0);
      if (fabs(m[0][2]-m02)>eps) return(0);
      if (fabs(m[0][3]-m03)>eps) return(0);
      if (fabs(m[1][0]-m10)>eps) return(0);
      if (fabs(m[1][1]-m11)>eps) return(0);
      if (fabs(m[1][2]-m12)>eps) return(0);
      if (fabs(m[1][3]-m13)>eps) return(0);
      if (fabs(m[2][0]-m20)>eps) return(0);
      if (fabs(m[2][1]-m21)>eps) return(0);
      if (fabs(m[2][2]-m22)>eps) return(0);
      if (fabs(m[2][3]-m23)>eps) return(0);
      if (fabs(m[3][0]-m30)>eps) return(0);
      if (fabs(m[3][1]-m31)>eps) return(0);
      if (fabs(m[3][2]-m32)>eps) return(0);
      if (fabs(m[3][3]-m33)>eps) return(0);

      return(1);
   }

   return(0);
}

void M44_Print (const real8 *m, const char *title)
{
   if (m)
   {
      if (title) printf("<m3x3 title=\"%s\">\n",title);
      else       printf("<m3x3>\n");

      printf("%7.4lf,%7.4lf,%7.4lf,%7.4lf," , m[ 0], m[ 1], m[ 2], m[ 3] );
      printf("%7.4lf,%7.4lf,%7.4lf,%7.4lf," , m[ 4], m[ 5], m[ 6], m[ 7] );
      printf("%7.4lf,%7.4lf,%7.4lf,%7.4lf," , m[ 8], m[ 9], m[10], m[11] );
      printf("%7.4lf,%7.4lf,%7.4lf,%7.4lf\n", m[12], m[13], m[14], m[15] );

      printf("</m3x3>\n");
   }
}

