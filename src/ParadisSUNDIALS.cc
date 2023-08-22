/***************************************************************************
 *
 *      Module:       ParadisSUNDIALS.c
 *      Description:  SUNDIALS is a SUite of Nonlinear and
 *                    DIfferential/ALgebraic equation Solvers.  ParaDiS
 *                    makes use of one component of this software
 *                    suite (KINSOL) for one of its timestep integrators.
 *
 *                    This module contains the ParaDiS implementation of
 *                    an NVECTOR for the SUNDIALS Suite interface.
 *
 *      Includes private functions
 *          Vaxpy()
 *          VCopy()
 *          VDiff()
 *          VNeg()
 *          VScaleBy()
 *          VScaleDiff()
 *          VScaleSum()
 *          VLin1()
 *          VLin2()
 *          VSum()
 *
 *      Includes public functions
 *          NewEmptyParaDiSVector()
 *          NewParaDiSVector()
 *          CloneEmptyParaDiSVector()
 *          CloneParaDiSVector()
 *          FreeParaDiSVector()
 *          PrintParaDiSVector()
 *          P_VAbs()
 *          P_VAddConst()
 *          P_VClone()
 *          P_VCloneEmpty()
 *          P_VCompare()
 *          P_VConst()
 *          P_VConstrMask()
 *          P_VDestroy()
 *          P_VDiv()
 *          P_VDotProd()
 *          P_VInv()
 *          P_VInvTest()
 *          P_VL1Norm()
 *          P_VLinearSum()
 *          P_VMaxNorm()
 *          P_VMin()
 *          P_VMinQuotient()
 *          P_VProd()
 *          P_VScale()
 *          P_VSpace()
 *          P_VWL2Norm()
 *          P_VWrmsNorm()
 *          P_VWrmsNormMask()
 *
 *--------------------------------------------------------------------------*/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"

#if defined(USE_KINSOL) || defined(USE_ARKODE)
#include "ParadisSUNDIALS.h"

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/*
 *      Prototypes for the private helper functions
 */

/* z = x */
static void VCopy(ParaDiS_Vector x, ParaDiS_Vector z);
/* z = x + y */
static void VSum(ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z);
/* z = x - y */
static void VDiff(ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z);
/* z = -x */
static void VNeg(ParaDiS_Vector x, ParaDiS_Vector z);
/* z = c ( x + y ) */
static void VScaleSum(realtype c, ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z);
/* z = c ( x - y ) */
static void VScaleDiff(realtype c, ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z);
/* z = a x + y */
static void VLin1(realtype a, ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z);
/* z = a x - y */
static void VLin2(realtype a, ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z);
/* y <- a x + y */
static void Vaxpy(realtype a, ParaDiS_Vector x, ParaDiS_Vector y);
/* x <- a x */
static void VScaleBy(realtype a, ParaDiS_Vector x);

#ifdef PARALLEL
/* Reduction operations add/max/min over the processor group */
static realtype VAllReduce(realtype d, int op, MPI_Comm comm);
#else
#define MPI_COMM_WORLD 0
#define VAllReduce(a,b,c) (a)
#endif


/*--------------------------------------------------------------------------
 *      Functions exported
 *-------------------------------------------------------------------------*/

/***************************************************************************
 *
 *      Function:    NewEmptyParaDiSVector()
 *
 *      Description: Create a new ParaDiS_Vector with the operation
 *                   structure initialized and other contents zeroed.
 *
 *      Returns: Pointer to new ParaDiS_Vector or NULL on error
 *
 ***************************************************************************/
ParaDiS_Vector NewEmptyParaDiSVector()
{
    ParaDiS_Vector        v;
    ParaDiS_Vector_Ops    ops;
    ParaDiS_VectorContent content;

/*
 *  Create vector
 */
    v = (ParaDiS_Vector)N_VNewEmpty();

    if (v == (ParaDiS_Vector)NULL) {
        return((ParaDiS_Vector)NULL);
    }

/*
 *  Attach vector operations to ops structure
 */
    v->ops->nvclone           = P_VClone;
    v->ops->nvcloneempty      = P_VCloneEmpty;
    v->ops->nvdestroy         = P_VDestroy;
    v->ops->nvspace           = P_VSpace;
    v->ops->nvgetarraypointer = P_VGetArrayPointer;
    v->ops->nvsetarraypointer = P_VSetArrayPointer;
    v->ops->nvgetlength       = P_VGetLength;
    v->ops->nvgetcommunicator = P_VGetCommunicator;
    v->ops->nvlinearsum       = P_VLinearSum;
    v->ops->nvconst           = P_VConst;
    v->ops->nvprod            = P_VProd;
    v->ops->nvdiv             = P_VDiv;
    v->ops->nvscale           = P_VScale;
    v->ops->nvabs             = P_VAbs;
    v->ops->nvinv             = P_VInv;
    v->ops->nvaddconst        = P_VAddConst;
    v->ops->nvdotprod         = P_VDotProd;
    v->ops->nvmaxnorm         = P_VMaxNorm;
    v->ops->nvwrmsnormmask    = P_VWrmsNormMask;
    v->ops->nvwrmsnorm        = P_VWrmsNorm;
    v->ops->nvmin             = P_VMin;
    v->ops->nvwl2norm         = P_VWL2Norm;
    v->ops->nvl1norm          = P_VL1Norm;
    v->ops->nvcompare         = P_VCompare;
    v->ops->nvinvtest         = P_VInvTest;
    v->ops->nvconstrmask      = P_VConstrMask;
    v->ops->nvminquotient     = P_VMinQuotient;

/*
 *  Create content
 */
    content = (ParaDiS_VectorContent)malloc(sizeof(struct _ParaDiS_VectorContent));

    if (content == (ParaDiS_VectorContent)NULL) {
        N_VDestroy(v);
        return((ParaDiS_Vector)NULL);
    }

/*
 *  Attach lengths
 */
    content->active_length  = 0;
    content->local_length   = 0;
    content->global_length  = 0;
    content->own_data       = 0;
    content->data           = NULL;

/*
 *  Attach content
 */
    v->content = content;

    return(v);
}


/***************************************************************************
 *
 *      Function:    NewParaDiSVector()
 *
 *      Description: Create a new ParaDiS_Vector and initialize the
 *                   contents based on the current nodal data provided
 *                   in <home>.
 *
 *      Returns: Pointer to new ParaDiS_Vector or NULL on error
 *
 *      NOTE: This version of the ParaDiS Vector over allocates the array
 *            and sets the lengths based on the filled segment of the vector.
 *
 ***************************************************************************/
ParaDiS_Vector NewParaDiSVector(Home_t *home)
{
    int             newNodeKeyPtr, ghostNodeCount;
    long int        i, j;
    long int        allocate_length;
    long int        active_length, local_length, global_length;
    Node_t        **nodeKeys, **ghostNodeList;
    realtype       *data;
    ParaDiS_Vector  v;

/*
 *  shortcuts to home data
 */
    nodeKeys       = home->nodeKeys;
    newNodeKeyPtr  = home->newNodeKeyPtr;
    ghostNodeList  = home->ghostNodeList;
    ghostNodeCount = home->ghostNodeCount;

/*
 *  create new empty vector
 */
    v = NewEmptyParaDiSVector();

    if (v == (ParaDiS_Vector)NULL) {
        return((ParaDiS_Vector)NULL);
    }

/*
 *  Allocate memory
 */
    allocate_length = 3*(newNodeKeyPtr + ghostNodeCount);
    // if (allocate_length <= 0) { printf("allocate_length <= 0!\n"); allocate_length = 1; }
    data = (realtype *) malloc(allocate_length * sizeof(realtype));

    if (data == (realtype *)NULL) {
        P_VDestroy(v);
        return((ParaDiS_Vector)NULL);
    }

/*
 *  Copy data to vector
 */
    j = 0;
    active_length = 0;

    for (i = 0; i < newNodeKeyPtr; i++) {

        if (nodeKeys[i] == (Node_t *)NULL) {
            continue;
        }

        data[j]   = nodeKeys[i]->x;
        data[j+1] = nodeKeys[i]->y;
        data[j+2] = nodeKeys[i]->z;

        j += 3;
    }

    active_length = j;

    for (i = 0; i < ghostNodeCount; i++) {

        data[j]   = ghostNodeList[i]->x;
        data[j+1] = ghostNodeList[i]->y;
        data[j+2] = ghostNodeList[i]->z;

        j += 3;
    }

    local_length = j;

/*
 *  Set global length (active components only does not include ghosts)
 */
#ifdef PARALLEL
    MPI_Allreduce(&active_length, &global_length, 1, PVEC_INTEGER_MPI_TYPE,
                  MPI_SUM, MPI_COMM_WORLD);
#else
    global_length = active_length;
#endif

/*
 *  Attach lengths
 */
    PV_ACTLENGTH(v)  = active_length;
    PV_LOCLENGTH(v)  = local_length;
    PV_GLOBLENGTH(v) = global_length;

/*
 *  Attach data
 */
    PV_OWN_DATA(v) = 1;
    PV_DATA(v)     = data;

    return(v);
}


/***************************************************************************
 *
 *      Function:    CloneEmptyParaDiSVector()
 *
 *      Description: Exported wrapper to call P_VCloneEmpty which creates a
 *                   limited copy of an existing ParaDiS_Vector. Copy
 *                   includes the vector operation structure, but not the
 *                   other contents.
 *
 *      Parameters:
 *          IN: w   Original ParaDiS_Vector to be copied
 *
 *      Returns:  A new ParaDiS_Vector which is an empty clone of the
 *                vector <w>, or NULL on error.
 *
 ***************************************************************************/
ParaDiS_Vector CloneEmptyParaDiSVector(ParaDiS_Vector w)
{
    ParaDiS_Vector v;

    v = P_VCloneEmpty(w);

    return(v);
}


/***************************************************************************
 *
 *      Function:    CloneParaDiSVector()
 *
 *      Description: Exported wrapper to call P_VClone which Creates a copy
 *                   of an existing ParaDiS_Vector. Copy includes the vector
 *                   operation structure and the other contents.
 *
 *      Parameters:
 *          IN: w   Original ParaDiS_Vector to be copied
 *
 *      Returns:  A new ParaDiS_Vector which is a full clone of the
 *                vector <w>, or NULL on error.
 *
 ***************************************************************************/
ParaDiS_Vector CloneParaDiSVector(ParaDiS_Vector w)
{
    ParaDiS_Vector v;

    v = P_VClone(w);

    return(v);
}


/***************************************************************************
 *
 *      Function:    FreeParaDiSVector()
 *
 *      Description: Exported wrapper to call P_VDestroy which frees all
 *                   resources associated with the specified ParaDiS_Vector
 *                   <v>.
 *
 ***************************************************************************/
void FreeParaDiSVector(ParaDiS_Vector v)
{
    P_VDestroy(v);
}


/***************************************************************************
 *
 *      Function:    PrintParaDiSVector()
 *
 *      Description: Writes ParaDiSVector data to stdout.
 *
 ***************************************************************************/
void PrintParaDiSVector(ParaDiS_Vector x)
{
    long int  i, N;
    realtype *xd;

    N  = PV_ACTLENGTH(x);
    xd = PV_DATA(x);

    for (i = 0; i < N; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
        printf("%Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
        printf("%23.16e\n", xd[i]);
#else
        printf("%g\n", xd[i]);
#endif
    }
    printf("\n");

    return;
}


/*
 * -----------------------------------------------------------------
 * implementation of vector operations in table
 * -----------------------------------------------------------------
 */


/***************************************************************************
 *
 *      Function:    P_VCloneEmpty()
 *
 *      Description: Creates a limited copy of an existing ParaDiS_Vector.
 *                   Copy includes the vector operation structure, but
 *                   not the other contents.
 *
 *      Parameters:
 *          IN: w   Original ParaDiS_Vector to be copied
 *
 *      Returns:  A new ParaDiS_Vector which is an empty clone of the
 *                vector <w>, or NULL on error.
 *
 ***************************************************************************/
ParaDiS_Vector P_VCloneEmpty(ParaDiS_Vector w)
{
    ParaDiS_Vector        v;
    ParaDiS_Vector_Ops    ops;
    ParaDiS_VectorContent content;

    if (w == (ParaDiS_Vector)NULL) {
        return((ParaDiS_Vector)NULL);
    }

/*
 *  Create vector
 */
    v = (ParaDiS_Vector)N_VNewEmpty();

    if (v == (ParaDiS_Vector)NULL) {
        return((ParaDiS_Vector)NULL);
    }

/*
 *  Copy operations use for w to v
 */
    N_VCopyOps(w, v);

/*
 *  Create content
 */
    content = (ParaDiS_VectorContent)malloc(sizeof(struct _ParaDiS_VectorContent));

    if (content == (ParaDiS_VectorContent)NULL) {
        free(ops);
        free(v);
        return((ParaDiS_Vector)NULL);
    }

/*
 *  Attach lengths and data
 */
    content->active_length  = PV_ACTLENGTH(w);
    content->local_length   = PV_LOCLENGTH(w);
    content->global_length  = PV_GLOBLENGTH(w);
    content->own_data       = 0;
    content->data           = NULL;

/*
 *  Attach content
 */
    v->content = content;

    return(v);
}


/***************************************************************************
 *
 *      Function:    P_VClone()
 *
 *      Description: Creates a copy of an existing ParaDiS_Vector
 *                   Copy includes the vector operation structure and
 *                   the other contents.
 *
 *      Parameters:
 *          IN: w   Original ParaDiS_Vector to be copied
 *
 *      Returns:  A new ParaDiS_Vector which is a full clone of the
 *                vector <w>, or NULL on error.
 *
 ***************************************************************************/
ParaDiS_Vector P_VClone(ParaDiS_Vector w)
{
    long int        local_length;
    realtype       *data;
    ParaDiS_Vector  v;

    v = P_VCloneEmpty(w);

    if (v == (ParaDiS_Vector)NULL) {
       return((ParaDiS_Vector)NULL);
    }

    local_length = PV_LOCLENGTH(w);

/*
 *  Create data
 */
    if(local_length > 0) {

/*
 *      Allocate memory
 */
        data = (realtype *)malloc(local_length * sizeof(realtype));

        if (data == (realtype *)NULL) {
            P_VDestroy(v);
            return((ParaDiS_Vector)NULL);
        }

/*
 *      Attach data
 */
        PV_OWN_DATA(v) = 1;
        PV_DATA(v)     = data;
    }

    return(v);
}


/***************************************************************************
 *
 *      Function:    P_VDestroy()
 *
 *      Description: Frees the specified ParaDiS_Vector and all resources
 *                   associated with it.
 *
 *      Parameters:
 *          IN:  v  ParaDiS_Vector to be freed.
 *
 ***************************************************************************/
void P_VDestroy(ParaDiS_Vector v)
{
    if ((PV_OWN_DATA(v) == 1) && (PV_DATA(v) != NULL)) {
        free(PV_DATA(v));
        PV_DATA(v) = NULL;
    }

    free(v->content); v->content = NULL;
    free(v->ops); v->ops = NULL;
    free(v); v = NULL;

    return;
}


/***************************************************************************
 *
 *      Function:    P_VSpace()
 *
 *      Description: Number of realtype words, numer of ineteger words
 *
 ***************************************************************************/
void P_VSpace(ParaDiS_Vector v, long int *lrw, long int *liw)
{
    int npes;

#ifdef PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
#else
    npes = 1;
#endif

    *lrw = PV_GLOBLENGTH(v);
    *liw = 2*npes;

    return;
}


/***************************************************************************
 *
 *      Function:     P_VGetArrayPointer()
 *
 *      Description:  Returns pointer to realtype array from ParaDiS_Vector
 *                    <v>.
 *
 ***************************************************************************/
realtype *P_VGetArrayPointer(ParaDiS_Vector v)
{
    return((realtype *) PV_DATA(v));
}


/***************************************************************************
 *
 *      Function:     P_VSetArrayPointer()
 *
 *      Description:  Overwrites the data in ParaDiS_Vector <v> with given
 *                    array of realtype.
 *
 ***************************************************************************/
void P_VSetArrayPointer(realtype *v_data, ParaDiS_Vector v)
{
    if (PV_LOCLENGTH(v) > 0) {
        PV_DATA(v) = v_data;
    }

    return;
}

/***************************************************************************
 *
 *      Function:     P_VGetLength()
 *
 *      Description:  Returns the global length of the ParaDiS_Vector <v>.
 *
 ***************************************************************************/
sunindextype P_VGetLength(ParaDiS_Vector v)
{
    return PV_GLOBLENGTH(v);
}

/***************************************************************************
 *
 *      Function:     P_VGetCommunicator()
 *
 *      Description:  Returns the MPI communicator of the ParaDiS_Vector <v>.
 *
 ***************************************************************************/
void* P_VGetCommunicator(ParaDiS_Vector v)
{
    return ((void*)MPI_COMM_WORLD);
}


/***************************************************************************
 *
 *      Function:     P_VLinearSum()
 *
 *      Description:  Computes z = a x + b y
 *
 ***************************************************************************/
void P_VLinearSum(realtype a, ParaDiS_Vector x, realtype b, ParaDiS_Vector y, ParaDiS_Vector z)
{
    long int       i, N;
    realtype       c, *xd, *yd, *zd;
    booleantype    test;
    ParaDiS_Vector v1, v2;

    xd = yd = zd = NULL;

/*
 *  BLAS usage: axpy y <- ax+y
 */
    if ((b == ONE) && (z == y)) {
        Vaxpy(a, x, y);
        return;
    }

/*
 *  BLAS usage: axpy x <- by+x
 */
    if ((a == ONE) && (z == x)) {
        Vaxpy(b, y, x);
        return;
    }

/*
 *  Case: a == b == 1.0
 */
    if ((a == ONE) && (b == ONE)) {
        VSum(x, y, z);
        return;
    }

/*
 *  Cases: (1) a ==  1.0, b == -1.0
 *         (2) a == -1.0, b ==  1.0
 */
    if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
        v1 = test ? y : x;
        v2 = test ? x : y;
        VDiff(v2, v1, z);
        return;
    }

/*
 *  Cases:(1) a == 1.0, b == other or 0.0,
 *        (2) a == other or 0.0, b == 1.0
 *
 *  If a or b is 0.0, then user should have called N_VScale
 */
    if ((test = (a == ONE)) || (b == ONE)) {
        c = test ? b : a;
        v1 = test ? y : x;
        v2 = test ? x : y;
        VLin1(c, v1, v2, z);
        return;
    }

/*
 *  Cases: (1) a == -1.0, b != 1.0
 *         (2) a != 1.0, b == -1.0
 */
    if ((test = (a == -ONE)) || (b == -ONE)) {
        c = test ? b : a;
        v1 = test ? y : x;
        v2 = test ? x : y;
        VLin2(c, v1, v2, z);
        return;
    }

/*
 *  Case: a == b
 *
 *  Catches case both a and b are 0.0 - user should have
 *  called N_VConst
 */
    if (a == b) {
        VScaleSum(a, x, y, z);
        return;
    }

/*
 *  Case: a == -b
 */
    if (a == -b) {
        VScaleDiff(a, x, y, z);
        return;
    }

/*
 *  Do all cases not handled above:
 *      (1) a == other, b == 0.0 - user should have called N_VScale
 *      (2) a == 0.0, b == other - user should have called N_VScale
 *      (3) a,b == other, a !=b, a != -b
 */
    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    yd = PV_DATA(y);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,a,b,xd,yd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = (a*xd[i])+(b*yd[i]);
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VLinearSum\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    P_VConst()
 *
 *      Description: Computes z[i] = c for i = 0, 1, ..., N-1
 *
 ***************************************************************************/
void P_VConst(realtype c, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *zd;

    zd = NULL;

    N  = PV_LOCLENGTH(z);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,c,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = c;
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VConst\n");
    PrintParaDiSVector(z);
#endif

    return;
}

/***************************************************************************
 *
 *      Function:    P_VProd()
 *
 *      Description: Computes z[i] = x[i] * y[i] for i = 0, 1, ..., N-1
 *
 ***************************************************************************/
void P_VProd(ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *yd, *zd;

    xd = yd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    yd = PV_DATA(y);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,yd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = xd[i]*yd[i];
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VProd\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    P_VDiv()
 *
 *      Description: z[i] = x[i] / y[i] for i = 0, 1, ..., N-1
 *
 ***************************************************************************/
void P_VDiv(ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *yd, *zd;

    xd = yd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    yd = PV_DATA(y);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,yd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = xd[i]/yd[i];
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VDiv\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    P_VScale()
 *
 *      Description: Computes z = c * x
 *
 ***************************************************************************/
void P_VScale(realtype c, ParaDiS_Vector x, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *zd;

    xd = zd = NULL;

/*
 *  BLAS usage: scale x <- cx
 */
    if (z == x) {
        VScaleBy(c, x);
        return;
    }

    if (c == ONE) {
        VCopy(x, z);
    } else if (c == -ONE) {
        VNeg(x, z);
    } else {
        N  = PV_LOCLENGTH(x);
        xd = PV_DATA(x);
        zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,c,xd,zd) schedule(static)
#endif
        for (i = 0; i < N; i++) {
            zd[i] = c*xd[i];
        }

#ifdef DEBUG_PARADIS_VECTOR
        printf("P_VScale\n");
        PrintParaDiSVector(z);
#endif
    }

    return;
}


/***************************************************************************
 *
 *      Function:    P_VAbs()
 *
 *      Description: Computes z[i] = |x[i]|, for i = 0, 1, ..., N-1
 *
 ***************************************************************************/
void P_VAbs(ParaDiS_Vector x, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *zd;

    xd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = SUNRabs(xd[i]);
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VAbs\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    P_VInv()
 *
 *      Description: Computes z[i] = 1.0 / x[i], for i = 0, 1, ..., N-1
 *
 ***************************************************************************/
void P_VInv(ParaDiS_Vector x, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *zd;

    xd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = ONE/xd[i];
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VInv\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    P_VAddConst()
 *
 *      Description: Computes z[i] = x[i] + b, for i = 0, 1, ..., N-1
 *
 ***************************************************************************/
void P_VAddConst(ParaDiS_Vector x, realtype b, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *zd;

    xd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,b,xd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = xd[i]+b;
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VAddConst\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    P_VDotProd()
 *
 *      Description: Computes the ordinary dot product of x and y
 *                   (i.e. sum (i=0 to N-1) {x[i] * y[i]}) on the locally
 *                   owned data, then sums the results from all MPI tasks.
 *
 *      Returns: Dot product or 0.0 if N <= 0;
 *
 ***************************************************************************/
realtype P_VDotProd(ParaDiS_Vector x, ParaDiS_Vector y)
{
    long int i, N;
    realtype sum, *xd, *yd, gsum;

    sum = ZERO;
    xd = yd = NULL;

/*
 *  use active length
 */
    N  = PV_ACTLENGTH(x);
    xd = PV_DATA(x);
    yd = PV_DATA(y);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,yd) reduction(+:sum) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        sum += xd[i]*yd[i];
    }

    gsum = VAllReduce(sum, 1, MPI_COMM_WORLD);

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VDotProd: \n");
    printf("%23.16e \n\n",gsum);
#endif

    return(gsum);
}


/***************************************************************************
 *
 *      Function:    P_VMaxNorm()
 *
 *      Description: Computes the maximum norm of x
 *                   (i.e. max (i=0 to N-1) |x[i]|) on the locally owned
 *                   data, then performs a global reduction to find the
 *                   largest value among all the MPI tasks.
 *
 *      Returns: The global maximum norm or 0.0 if N <= 0;
 *
 ***************************************************************************/
realtype P_VMaxNorm(ParaDiS_Vector x)
{
    long int i, N;
    realtype tempmax, localMax, *xd, gmax;

    localMax = ZERO;
    xd = NULL;

/*
 *  Compute norm with active length
 */
    N  = PV_ACTLENGTH(x);
    xd = PV_DATA(x);

#ifdef _OPENMP
#pragma omp parallel default(none) private(i,tempmax) shared(N,localMax,xd)
#endif
    {
        tempmax = ZERO;

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (i = 0; i < N; i++) {
            if (SUNRabs(xd[i]) > tempmax) {
                tempmax = SUNRabs(xd[i]);
            }
        }
#ifdef _OPENMP
#pragma omp critical (P_VMAXNORM)
#endif
        {
            if (tempmax > localMax) {
                localMax = tempmax;
            }
        }
    }

    gmax = VAllReduce(localMax, 2, MPI_COMM_WORLD);

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VMaxNorm\n");
    printf("%13.16e\n\n",gmax);
#endif

    return(gmax);
}


/***************************************************************************
 *
 *      Function:    P_VWrmsNorm()
 *
 *      Description: Computes the weighted root mean square norm of x with
 *                   weight vector w:
 *
 *                       sqrt [(sum (i=0 to N-1) {(x[i] * w[i])^2}) / N]
 *
 *                   This calculation is done for the locally owned data,
 *                   then a global reduction is performed to sum the results
 *                   from all MPI tasks.
 *
 *      Returns: The global weighted root mean square norm
 *
 ***************************************************************************/
realtype P_VWrmsNorm(ParaDiS_Vector x, ParaDiS_Vector w)
{
//     long int i, N, N_global;
//     realtype sum, prodi, *xd, *wd, gsum;

//     sum = ZERO;
//     xd = wd = NULL;

// /*
//  * Compute norm with active length
//  */
//     N        = PV_ACTLENGTH(x);
//     N_global = PV_GLOBLENGTH(x);
//     xd       = PV_DATA(x);
//     wd       = PV_DATA(w);

// #ifdef _OPENMP
// #pragma omp parallel for default(none) private(i,prodi) shared(N,xd,wd) reduction(+:sum) schedule(static)
// #endif
//     for (i = 0; i < N; i++) {
//         prodi = xd[i]*wd[i];
//         sum += SUNSQR(prodi);
//     }

//     gsum = VAllReduce(sum, 1, MPI_COMM_WORLD);

// #ifdef DEBUG_PARADIS_VECTOR
//     printf("P_VWrmsNorm\n");
//     printf("%23.16e\n\n",SUNRsqrt(gsum/N_global));
// #endif

//     return(SUNRsqrt(gsum/N_global));

    return P_VMaxNorm(x);
}


/***************************************************************************
 *
 *      Function:    P_VWrmsNormMask()
 *
 *      Description: Computes the weighted root mean square norm of x with
 *                   weight vector w only for elements of x corresponding to
 *                   nonzero elements of <id>:
 *
 *                   sqrt [(sum (i=0 to N-1) {(sign(id[i]) * x[i] * w[i])^2}) / N]
 *
 *                   This calculation is done for the locally owned data,
 *                   then a global reduction is performed to sum the results
 *                   from all MPI tasks.
 *
 *      Returns: The global weighted root mean square norm
 *
 ***************************************************************************/
realtype P_VWrmsNormMask(ParaDiS_Vector x, ParaDiS_Vector w, ParaDiS_Vector id)
{
    long int i, N, N_global;
    realtype sum, prodi, *xd, *wd, *idd, gsum;

    sum = ZERO;
    xd = wd = idd = NULL;

/*
 *  Compute norm with active length
 */
    N        = PV_ACTLENGTH(x);
    N_global = PV_GLOBLENGTH(x);
    xd       = PV_DATA(x);
    wd       = PV_DATA(w);
    idd      = PV_DATA(id);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i,prodi) shared(N,xd,wd,idd) reduction(+:sum) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        if (idd[i] > ZERO) {
            prodi = xd[i]*wd[i];
            sum += SUNSQR(prodi);
        }
    }

    gsum = VAllReduce(sum, 1, MPI_COMM_WORLD);

    return(SUNRsqrt(gsum/N_global));
}


/***************************************************************************
 *
 *      Function:    P_VMin()
 *
 *      Description: Finds the minimum value in x. This calculation is
 *                   done on the locally owned data, then a global
 *                   reduction is performed to obtain the minimum
 *                   from among all the MPI tasks.
 *
 *      Returns: min x[i] if N > 0 or BIG_REAL if N <= 0
 *
 ***************************************************************************/
realtype P_VMin(ParaDiS_Vector x)
{
    long int i, N;
    realtype tempmin, localMin, *xd, gmin;

    localMin = BIG_REAL;
    xd  = NULL;

/*
 *  Compute with active length
 */
    N = PV_ACTLENGTH(x);

    if (N > 0) {
        xd = PV_DATA(x);

#ifdef _OPENMP
#pragma omp parallel default(none) private(i,tempmin) shared(N,localMin,xd)
#endif
        {
            tempmin = xd[0];

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (i = 1; i < N; i++) {
                if (xd[i] < tempmin) {
                    tempmin = xd[i];
                }
            }
#ifdef _OPENMP
#pragma omp critical (P_VMIN)
#endif
            {
                if (tempmin < localMin) {
                    localMin = tempmin;
                }
            }
        }

    }  /* end if (N > 0) */

    gmin = VAllReduce(localMin, 3, MPI_COMM_WORLD);

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VMin\n");
    printf("%23.16e\n\n",gmin);
#endif

    return(gmin);
}


/***************************************************************************
 *
 *      Function:    P_VWL2Norm()
 *
 *      Description: Computes the weighted L2 norm of x with
 *                   weight vector w:
 *
 *                       sqrt [(sum (i=0 to N-1) {(x[i] * w[i])^2}) ]
 *
 *                   This calculation is done the locally owned data, then
 *                   a global reduction is performed to sum the results
 *                   from all MPI tasks.
 *
 *      Returns: The global weighted L2 norm
 *
 ***************************************************************************/
realtype P_VWL2Norm(ParaDiS_Vector x, ParaDiS_Vector w)
{
    long int i, N;
    realtype sum, prodi, *xd, *wd, gsum;


    sum = ZERO;
    xd = wd = NULL;

/*
 *  Compute norm with active length
 */
    N  = PV_ACTLENGTH(x);
    xd = PV_DATA(x);
    wd = PV_DATA(w);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i,prodi) shared(N,xd,wd) reduction(+:sum) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        prodi = xd[i]*wd[i];
        sum += SUNSQR(prodi);
    }

    gsum = VAllReduce(sum, 1, MPI_COMM_WORLD);

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VWL2Norm\n");
    printf("%23.16e\n\n",SUNRsqrt(gsum));
#endif

    return(SUNRsqrt(gsum));
}


/***************************************************************************
 *
 *      Function:    P_VL1Norm()
 *
 *      Description: Computes the L1 norm of x:
 *
 *                   sqrt [(sum (i=0 to N-1) abs(x[i])
 *
 *                   This calculation is done the locally owned data, then
 *                   a global reduction is performed to sum the results
 *                   from all MPI tasks.
 *
 *      Returns: The global L1 norm
 *
 ***************************************************************************/
realtype P_VL1Norm(ParaDiS_Vector x)
{
    long int i, N;
    realtype sum, gsum, *xd;

    sum = ZERO;
    xd = NULL;

/*
 *  Compute norm with active length
 */
    N  = PV_ACTLENGTH(x);
    xd = PV_DATA(x);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd) reduction(+:sum) schedule(static)
#endif
    for (i = 0; i<N; i++) {
        sum += SUNRabs(xd[i]);
    }

    gsum = VAllReduce(sum, 1, MPI_COMM_WORLD);

    return(gsum);
}


/***************************************************************************
 *
 *      Function:    P_VCompare
 *
 *      Description: Compares the components of ParaDiS_Vector <x>
 *                   to a constant <c>
 *
 *      Returns: z[i] = 1.0 if |x[i]| >= c,
 *               z[i] = 0.0 if |x[i]| < c,
 *
 ***************************************************************************/
void P_VCompare(realtype c, ParaDiS_Vector x, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *zd;

    xd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,c,xd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO;
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VCompare\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    P_VInvTest
 *
 *      Description: Computes z[i] = 1.0 / x[i] for all non-zero x[i]
 *                   This calculation is done on al locally owned
 *                   data, then a global reduction is done to determine
 *                   if the inversion failed on *any* MPI task
 *
 *      Returns: 1 if all components of x (on *all* MPI tasks) were
 *               non-zero (successful inversion) or 0 otherwise.
 *
 ***************************************************************************/
booleantype P_VInvTest(ParaDiS_Vector x, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *zd, val, gval;

    xd = zd = NULL;

    N  = PV_ACTLENGTH(x);
    xd = PV_DATA(x);
    zd = PV_DATA(z);

    val = ONE;

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,val,xd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        if (xd[i] == ZERO) {
            val = ZERO;
        } else {
            zd[i] = ONE/xd[i];
        }
    }

    gval = VAllReduce(val, 3, MPI_COMM_WORLD);

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VInvTest\n");
    printf("gval = %23.16e\n",gval);
    PrintParaDiSVector(z);
#endif

    return(gval == ZERO ? 0 : 1);
}


/***************************************************************************
 *
 *      Function:    P_VConstrMask()
 *
 *      Description: This routine is used only for constraint checking
 *                   for the linear solve.  If the linear solve does
 *                   have constraints, this routine is unused.
 *
 *                   This function performs the following checks:
 *
 *                       x_i >  0 if c_i = 2,
 *                       x_i >= 0 if c_i = 1,
 *                       x_i <= 0 if c_i =-1,
 *                       x_i <  0 if c_i =-2.
 *
 *                   There is no constraint on x_i if c_i=0.
 *
 *      Returns: 0 if any element failed the constraint test, 1 if
 *               all passed.  It also sets a mask vector, m, with elements
 *               equal to 1.0 where the constraint test failed, and 0.0
 *               where the test passed.
 *
 ***************************************************************************/
booleantype P_VConstrMask(ParaDiS_Vector c, ParaDiS_Vector x, ParaDiS_Vector m)
{
    long int  i, N;
    realtype  temp;
    realtype *cd, *xd, *md;

    cd = xd = md = NULL;

    N  = PV_ACTLENGTH(x);
    xd = PV_DATA(x);
    cd = PV_DATA(c);
    md = PV_DATA(m);

    temp = ONE;

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,cd,md,temp) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        md[i] = ZERO;

        if (cd[i] == ZERO) {
            continue;
        }

        if (cd[i] > ONEPT5 || cd[i] < -ONEPT5) {

            if (xd[i]*cd[i] <= ZERO) {
                temp = ZERO;
                md[i] = ONE;
            }

            continue;
        }

        if (cd[i] > HALF || cd[i] < -HALF) {
            if (xd[i]*cd[i] < ZERO ) {
                temp = ZERO;
                md[i] = ONE;
            }
        }

    }  /* end for (i = 0; i < N; i++) */

    temp = VAllReduce(temp, 3, MPI_COMM_WORLD);

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VConstrMask\n");
    printf("%23.16e\n\n",temp);
#endif

    return(temp == ONE ? 1 : 0);
}


/***************************************************************************
 *
 *      Function:    P_VMinQuotient()
 *
 *      Description: This routine returns the minimum of the quotients
 *                   obtained by term-wise dividing num_i by denom_i.
 *                   A zero element in denom will be skipped.
 *                   This operation is performed on all local data
 *                   and then a global reduction is performed to find
 *                   the minimum from among all MPI tasks.
 *
 *      Returns: Minimum quotient as described above, or the large value
 *               BIG_REAL (defined in the header file sundials_types.h)
 *               if no such quotients are found.
 *
 ***************************************************************************/
realtype P_VMinQuotient(ParaDiS_Vector num, ParaDiS_Vector denom)
{
    long int     i, N;
    booleantype  notEvenOnce;
    realtype    *nd, *dd, tempmin, localMin;

    nd = dd = NULL;

    N  = PV_ACTLENGTH(num);
    nd = PV_DATA(num);
    dd = PV_DATA(denom);

    notEvenOnce = 1;
    localMin = BIG_REAL;

#ifdef _OPENMP
#pragma omp parallel default(none) private(i,tempmin,notEvenOnce) shared(N,localMin,nd,dd)
#endif
    {
        tempmin = BIG_REAL;

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (i = 0; i < N; i++) {

            if (dd[i] == ZERO) {
                continue;
            } else {
                if (!notEvenOnce) {
                    tempmin = MIN(tempmin, nd[i]/dd[i]);
                } else {
                    tempmin = nd[i]/dd[i];
                    notEvenOnce = 0;
                }
            }
        }
#ifdef _OPENMP
#pragma omp critical (P_VMINQUOTIENT)
#endif
        {
            if (tempmin < localMin) {
                localMin = tempmin;
            }
        }
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("P_VMinQuotient\n");
    printf("%23.16e\n\n",VAllReduce(localMin, 3, MPI_COMM_WORLD));
#endif

    return(VAllReduce(localMin, 3, MPI_COMM_WORLD));
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */


#ifdef PARALLEL
/***************************************************************************
 *
 *      Function:    VAllReduce()
 *
 *      Description: This function does a global reduction.  The operation is
 *                       sum if op = 1,
 *                       max if op = 2,
 *                       min if op = 3.
 *                   The operation is over all processors in the communicator
 *
 *      Returns: Global result of the opeartion <out>.
 *
 ***************************************************************************/
static realtype VAllReduce(realtype d, int op, MPI_Comm comm)
{
    realtype out = 0;;

    switch (op) {
        case 1:
            MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
            break;

        case 2:
            MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MAX, comm);
            break;

        case 3:
            MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MIN, comm);
            break;

        default:
            break;
    }

    return(out);
}
#endif


/***************************************************************************
 *
 *      Function:    VCopy()
 *
 *      Description: Copy the contents of ParaDiS_Vector <x> into <z>
 *
 ***************************************************************************/
static void VCopy(ParaDiS_Vector x, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *zd;

    xd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = xd[i];
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("VCopy\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    VSum()
 *
 *      Description: Computes z = x + y for the specified ParaDiS_Vectors
 *
 ***************************************************************************/
static void VSum(ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *yd, *zd;

    xd = yd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    yd = PV_DATA(y);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,yd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = xd[i]+yd[i];
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("VSum\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    VDiff()
 *
 *      Description: Computes the difference z = x - y for the specified
 *                   ParaDiS_Vector's <x> and <y>
 *
 ***************************************************************************/
static void VDiff(ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *yd, *zd;

    xd = yd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    yd = PV_DATA(y);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,yd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = xd[i]-yd[i];
    }
#ifdef DEBUG_PARADIS_VECTOR
    printf("VDiff\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    VNeg()
 *
 *      Description: Computes z = -x for all components of the
 *                   ParaDiS_Vector <x>.
 *
 ***************************************************************************/
static void VNeg(ParaDiS_Vector x, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *zd;

    xd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = -xd[i];
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("VNeg\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    VScaleSum()
 *
 *      Description: Computes z = c (x + y) for the specified
 *                   ParaDiS_Vectors <x> and <y>, and constant <c>.
 *
 ***************************************************************************/
static void VScaleSum(realtype c, ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *yd, *zd;

    xd = yd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    yd = PV_DATA(y);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,c,xd,yd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = c*(xd[i]+yd[i]);
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("VScaleSum\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    VScaleDiff()
 *
 *      Description: Computes z = c (x - y) for the specified
 *                   ParaDiS_Vectors <x> and <y>, and constant <c>.
 *
 ***************************************************************************/
static void VScaleDiff(realtype c, ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *yd, *zd;

    xd = yd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    yd = PV_DATA(y);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,c,xd,yd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = c*(xd[i]-yd[i]);
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("VScaleDiff\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    VLin1()
 *
 *      Description: Computes z = a x + y for the specified
 *                   ParaDiS_Vectors <x> and <y>, and scalar <a>.
 *
 ***************************************************************************/
static void VLin1(realtype a, ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *yd, *zd;

    xd = yd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    yd = PV_DATA(y);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,a,xd,yd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = (a*xd[i])+yd[i];
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("VLin1\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    VLin2()
 *
 *      Description: Computes z = a x - y for the specified
 *                   ParaDiS_Vectors <x> and <y>, and scalar <a>.
 *
 ***************************************************************************/
static void VLin2(realtype a, ParaDiS_Vector x, ParaDiS_Vector y, ParaDiS_Vector z)
{
    long int  i, N;
    realtype *xd, *yd, *zd;

    xd = yd = zd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    yd = PV_DATA(y);
    zd = PV_DATA(z);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,a,xd,yd,zd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        zd[i] = (a*xd[i])-yd[i];
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("VLin2\n");
    PrintParaDiSVector(z);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    Vaxpy()
 *
 *      Description: Computes y = a x + y for the specified
 *                   ParaDiS_Vectors <x> and <y>, and scalar <a>.
 *
 ***************************************************************************/
static void Vaxpy(realtype a, ParaDiS_Vector x, ParaDiS_Vector y)
{
    long int  i, N;
    realtype *xd, *yd;

    xd = yd = NULL;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);
    yd = PV_DATA(y);

    if (a == ONE) {
#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,yd) schedule(static)
#endif
        for (i = 0; i < N; i++) {
            yd[i] += xd[i];
        }

#ifdef DEBUG_PARADIS_VECTOR
        printf("Vaxpy1\n");
        PrintParaDiSVector(y);
#endif

        return;

    }  /* end if (a == ONE) */

    if (a == -ONE) {
#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,xd,yd) schedule(static)
#endif
        for (i = 0; i < N; i++) {
            yd[i] -= xd[i];
        }

#ifdef DEBUG_PARADIS_VECTOR
        printf("Vaxpy2\n");
        PrintParaDiSVector(y);
#endif

        return;

    }  /* end if (a == -ONE) */

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,a,xd,yd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        yd[i] += a*xd[i];
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("Vaxpy3\n");
    PrintParaDiSVector(y);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    VScaleBy()
 *
 *      Description: Computes x = a x for the specified
 *                   ParaDiS_Vector <x> and scalar <a>.
 *
 ***************************************************************************/
static void VScaleBy(realtype a, ParaDiS_Vector x)
{
    long int  i, N;
    realtype *xd;

    N  = PV_LOCLENGTH(x);
    xd = PV_DATA(x);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(i) shared(N,a,xd) schedule(static)
#endif
    for (i = 0; i < N; i++) {
        xd[i] *= a;
    }

#ifdef DEBUG_PARADIS_VECTOR
    printf("VScaleBy\n");
    PrintParaDiSVector(x);
#endif

    return;
}


/***************************************************************************
 *
 *      Function:    ParaDiSVectorToHome_Positions()
 *
 *      Description: Move ParaDiS_Vector data to the ParaDiS <home> structure.
 *
 ***************************************************************************/
void ParaDiSVectorToHome_Positions(ParaDiS_Vector v, Home_t *home, int foldbox)
{
    int        newNodeKeyPtr, ghostNodeCount;
    long int   i, j;
    realtype  *data;
    Node_t   **nodeKeys, **ghostNodeList;
    Param_t   *param;

/*
 *  Get pointer to vector data
 */
    data = PV_DATA(v);

/*
 *  Shortcuts to home data
 */
    param          = home->param;
    nodeKeys       = home->nodeKeys;
    newNodeKeyPtr  = home->newNodeKeyPtr;
    ghostNodeList  = home->ghostNodeList;
    ghostNodeCount = home->ghostNodeCount;

/*
 *  Copy data to home
 */
    for (j = 0, i = 0; i < newNodeKeyPtr; i++) {

        if (nodeKeys[i] == (Node_t *)NULL) {
            continue;
        }

        nodeKeys[i]->x = data[j];
        nodeKeys[i]->y = data[j+1];
        nodeKeys[i]->z = data[j+2];

        if (foldbox == 1) {
            FoldBox(param, &nodeKeys[i]->x, &nodeKeys[i]->y, &nodeKeys[i]->z);
        }

        j += 3;
    }

    for (i = 0; i < ghostNodeCount; i++) {

        ghostNodeList[i]->x = data[j];
        ghostNodeList[i]->y = data[j+1];
        ghostNodeList[i]->z = data[j+2];

        if (foldbox == 1) {
            FoldBox(param, &ghostNodeList[i]->x, &ghostNodeList[i]->y,
                    &ghostNodeList[i]->z);
        }

        j += 3;
    }

    return;
}


/***************************************************************************
 *
 *      Function:    HomeToParaDiSVector_Positions()
 *
 *      Description: Move data from ParaDiS <home> structure to ParaDiS_Vector.
 *
 ***************************************************************************/
void HomeToParaDiSVector_Positions(Home_t *home, ParaDiS_Vector v, int foldbox)
{
    int        newNodeKeyPtr, ghostNodeCount;
    long int   i, j;
    realtype  *data;
    Node_t   **nodeKeys, **ghostNodeList;
    Param_t   *param;

/*
 *  Get pointer to vector data
 */
    data = PV_DATA(v);

/*
 *  Shortcuts to home data
 */
    param          = home->param;
    nodeKeys       = home->nodeKeys;
    newNodeKeyPtr  = home->newNodeKeyPtr;
    ghostNodeList  = home->ghostNodeList;
    ghostNodeCount = home->ghostNodeCount;

/*
 *  Copy data to vector
 */
    for (j = 0, i = 0; i < newNodeKeyPtr; i++) {

        if (nodeKeys[i] == (Node_t *)NULL) {
            continue;
        }

        if (foldbox == 1) {
            FoldBox(param, &nodeKeys[i]->x, &nodeKeys[i]->y, &nodeKeys[i]->z);
        }

        data[j]   = nodeKeys[i]->x;
        data[j+1] = nodeKeys[i]->y;
        data[j+2] = nodeKeys[i]->z;

        j += 3;
    }

    for (i = 0; i < ghostNodeCount; i++) {

        if (foldbox == 1) {
            FoldBox(param, &ghostNodeList[i]->x, &ghostNodeList[i]->y,
                    &ghostNodeList[i]->z);
        }

        data[j]   = ghostNodeList[i]->x;
        data[j+1] = ghostNodeList[i]->y;
        data[j+2] = ghostNodeList[i]->z;

        j += 3;
    }

    return;
}


/***************************************************************************
 *
 *      Function:    HomeToParaDiSVector_Velocity()
 *
 *      Description: Move data from ParaDiS <home> structure to ParaDiS_Vector.
 *
 ***************************************************************************/
void HomeToParaDiSVector_Velocities(Home_t *home, ParaDiS_Vector v)
{
    int        newNodeKeyPtr, ghostNodeCount;
    long int   i, j;
    realtype  *data;
    Node_t   **nodeKeys, **ghostNodeList;
    Param_t   *param;

/*
 *  Get pointer to vector data
 */
    data = PV_DATA(v);

/*
 *  Shortcuts to home data
 */
    param          = home->param;
    nodeKeys       = home->nodeKeys;
    newNodeKeyPtr  = home->newNodeKeyPtr;
    ghostNodeList  = home->ghostNodeList;
    ghostNodeCount = home->ghostNodeCount;

/*
 *  Copy data to vector
 */
    for (j = 0, i = 0; i < newNodeKeyPtr; i++) {

        if (nodeKeys[i] == (Node_t *)NULL) {
            continue;
        }

        data[j]   = nodeKeys[i]->vX;
        data[j+1] = nodeKeys[i]->vY;
        data[j+2] = nodeKeys[i]->vZ;

        j += 3;
    }

    for (i = 0; i < ghostNodeCount; i++) {

        data[j]   = ghostNodeList[i]->vX;
        data[j+1] = ghostNodeList[i]->vY;
        data[j+2] = ghostNodeList[i]->vZ;

        j += 3;
    }

    return;
}

#endif  /* ifdef USE_KINSOL or USE_ARKODE */
