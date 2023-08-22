#ifndef DGESVD_H
#define DGESVD_H

// The declarations below are included to provide an interface to the linear and SVD solvers within LAPACK.
//
// DGESV   - solves a real system of linear equations  A*X=B
// DGESVD  - computeis the singular value decomposition (SVD) of a real M-by-N matrix A, A=SVD(A)
//           optionally computing the left and/or right singular vectors.
//
// See http://performance.netlib.org/lapack for more details on using these routines.
//
// Note - the trailing underscore is needed on these routines because LAPACK is written
//        using FORTRAN. The FORTRAN compiler typically adds underscores to routines.
//------------------------------------------------------------------------------------------------------------

extern "C"
{
void dgesvd_(char JOBU[1], char JOBVT[1], int M[1], int N[1], double A[],
             int LDA[1], double S[], double U[], int LDU[1], double VT[],
             int LDVT[1], double WORK[], int LWORK[1], int INFO[1]);

void dgesv_ (int N[1], int NRHS[1], double A[], int LDA[1],
             int IPIV[], double B[], int LDB[1], int INFO[1]);
}

#endif   // DGESVD_H
