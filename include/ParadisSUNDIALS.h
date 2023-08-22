#pragma once

#ifndef _PDS_PARADIS_SUNDIALS_H
#define _PDS_PARADIS_SUNDIALS_H

/***************************************************************************
 *
 * Header: ParadisSUNDIALS.h 
 *
 * Description: SUNDIALS is a SUite of Nonlinear and DIfferential/ALgebraic
 *              equation Solvers. ParaDiS makes use of two components of
 *              this software for two of its timestep integrators, KINSOL
 *              and ARKODE.
 *
 *              This header defines the ParaDiS implementation of the
 *              NVECTOR module used for communication with the KINSOL and
 *              ARKODE. This implementation is based on the nvector_parallel
 *              implementation in SUNDIALS.
 *
 **************************************************************************/

// Pull in the Sundials includes first to avoid some ParaDiS macro 
// definitions that are  incompatible with Sundials.

#if defined(USE_KINSOL) || defined(USE_ARKODE)

#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>

#if defined(USE_ARKODE)
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_ls.h>
#endif

#include "Home.h"

/*
 *  define MPI data types
 */
#if defined(SUNDIALS_SINGLE_PRECISION)
#define PVEC_REAL_MPI_TYPE MPI_FLOAT

#elif defined(SUNDIALS_DOUBLE_PRECISION)
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define PVEC_REAL_MPI_TYPE MPI_LONG_DOUBLE

#endif

#define PVEC_INTEGER_MPI_TYPE MPI_LONG


/*
 *  ParaDiS implementation of the N_Vector 'content' structure
 */
typedef struct _ParaDiS_VectorContent {

    long int active_length; /* number of active vector components  */
    long int local_length;  /* vector length on local process active */
                            /* and ghosts                            */
    long int global_length; /* active vector length across all processes */
    booleantype own_data;   /* ownership of vector data */
    realtype *data;         /* pointer to vector data */ 

} *ParaDiS_VectorContent;

typedef N_Vector ParaDiS_Vector;
typedef N_Vector_Ops ParaDiS_Vector_Ops;


/*
 *  Macros to access ParaDiS vector content
 */
#define PV_CONTENT(v)    ( (ParaDiS_VectorContent)(v->content) )
#define PV_ACTLENGTH(v)  ( PV_CONTENT(v)->active_length )
#define PV_LOCLENGTH(v)  ( PV_CONTENT(v)->local_length )
#define PV_GLOBLENGTH(v) ( PV_CONTENT(v)->global_length )
#define PV_OWN_DATA(v)   ( PV_CONTENT(v)->own_data )
#define PV_DATA(v)       ( PV_CONTENT(v)->data )
#define PV_Ith(v,i)      ( PV_DATA(v)[i] )


/* -------------------------------------------------------------------------
 * Prototypes
 * -----------------------------------------------------------------------*/

/*
 *  Create vectors from the ParaDiS home structure.
 */
ParaDiS_Vector NewEmptyParaDiSVector();
ParaDiS_Vector NewParaDiSVector(Home_t *home);

/*
 *  Move ParaDiS_Vector data to the ParaDiS home structure.
 */
void ParaDiSVectorToHome_Positions(ParaDiS_Vector v, Home_t *home, int fold);
void HomeToParaDiSVector_Positions(Home_t *home, ParaDiS_Vector v, int fold);
void HomeToParaDiSVector_Velocities(Home_t *home, ParaDiS_Vector v);

/*
 *  Create new ParaDiS_Vectors from existing ParaDiS vectors.
 */
ParaDiS_Vector CloneEmptyParaDiSVector(ParaDiS_Vector w);
ParaDiS_Vector CloneParaDiSVector(ParaDiS_Vector w);

/*
 *  Free memory of ParaDiS_Vectors.
 */
void FreeParaDiSVector(ParaDiS_Vector v);

/*
 *  Print ParaDiS_Vector data
 */
void PrintParaDiSVector(ParaDiS_Vector v);

/*
 *  Vector operations
 */
ParaDiS_Vector P_VCloneEmpty(ParaDiS_Vector w); 
ParaDiS_Vector P_VClone(ParaDiS_Vector w); 

void P_VDestroy(ParaDiS_Vector v);
void P_VSpace(ParaDiS_Vector v, long int *lrw, long int *liw);
void P_VSetArrayPointer(realtype *v_data, ParaDiS_Vector v);
sunindextype P_VGetLength(ParaDiS_Vector v);
void* P_VGetCommunicator(ParaDiS_Vector v);
void P_VLinearSum(realtype a, ParaDiS_Vector x, realtype b, ParaDiS_Vector y,
                  ParaDiS_Vector z);
void P_VConst(realtype c, ParaDiS_Vector z);
void P_VProd(ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z);
void P_VDiv(ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z);
void P_VScale(realtype c, ParaDiS_Vector x, ParaDiS_Vector z);
void P_VAbs(ParaDiS_Vector x, ParaDiS_Vector z);
void P_VInv(ParaDiS_Vector x, ParaDiS_Vector z);
void P_VAddConst(ParaDiS_Vector x, realtype b, ParaDiS_Vector z);
void P_VCompare(realtype c, ParaDiS_Vector x, ParaDiS_Vector z);

realtype P_VDotProd(ParaDiS_Vector x, ParaDiS_Vector y);
realtype P_VMaxNorm(ParaDiS_Vector x);
realtype P_VWrmsNorm(ParaDiS_Vector x, ParaDiS_Vector w);
realtype P_VWrmsNormMask(ParaDiS_Vector x, ParaDiS_Vector w, ParaDiS_Vector id);
realtype P_VMin(ParaDiS_Vector x);
realtype P_VWL2Norm(ParaDiS_Vector x, ParaDiS_Vector w);
realtype P_VL1Norm(ParaDiS_Vector x);
realtype P_VMinQuotient(ParaDiS_Vector num, ParaDiS_Vector denom);
realtype *P_VGetArrayPointer(ParaDiS_Vector v);

booleantype P_VInvTest(ParaDiS_Vector x, ParaDiS_Vector z);
booleantype P_VConstrMask(ParaDiS_Vector c, ParaDiS_Vector x, ParaDiS_Vector m);

#endif  /* ifdef USE_KINSOL or USE_ARKODE */
#endif  /* ifndef _ParadisSUNDIALS_h */
