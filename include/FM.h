#pragma once

#ifndef _PDS_FM_H
#define _PDS_FM_H

#define X   0
#define Y   1
#define Z   2

#define NMAX        35    /* WARNING! MUST be at least param->fmOrder+4 */
#define NTMAX       (((NMAX+1)*(NMAX+2))>>1)
#define MAXORDER    (2*NMAX)

#define MAX_ORDER_ESHELBY     6

#define MP_COEFF     1
#define TAYLOR_COEFF 2


#define FFTW_FLAG FFTW_ESTIMATE
#define HOMO_THRESHOLD 1e-1 // threshold of homogeneous kernel

/*
 *  These macros are used in the FM code when calculating indices
 *  for all nearby neighbors of cells (where nearby is essentially
 *  all cells that are children of the parent (or children
 *  of the immediate neighbors of the parent) of cell <x>).
 *  Returned indices may be negative or larger than the dimensions
 *  of a block of cells, but adjustments for PBC are made elsewhere.
 */
#define SHIFTPOS(x)       (((((x) >> 1) + 1) << 1) + 1)
#define SHIFTNEG(x)       ((((x) >> 1) - 1) << 1)
/*
 *  Returns the index of <a> in a dimension of length <b> after
 *  correction for periodic boundary conditions.  If the initial
 *  condition is  0 <= a < b the returned value will be <a>,
 *  otherwise a will be adjusted.
 */
#define GETPBCINDEX(a,b)  (((a) % (b)) + ((b) * ((a) < 0)))


#ifdef BBFMM
typedef struct {

  int n; /* Chebyshev order */
  int n3;

  int NumGauss; /* Number of Gauss quadrature points */

  real8  homogen, AdjustLen;

  real8 *Kweights; /*Omega matrix */
  real8 *Cweights; /*Chebyshev interpolation coeffs. */

  real8 *Tkz; /* Chebyshev polynomials */

  int Ktable[343];

  /* Cut off from SVD for f and s */
  int cutofff;
  int cutoffs;

  /* Pre-computed matrices coming from SVD calculations */
  real8 *K;
  real8 *U;
  real8 *VT;

  /* Gauss weights and points */
  real8 *Gausspoints;
  real8 *Gaussweights;

  real8 ChebEps;

} bbFMM_t;

extern bbFMM_t *bbfmm;
#endif  /* ifdef BBFMM */

#ifdef UNIFORMFMM
typedef struct {

  int n; /* number of grid points */
  int n3;
  int pown3;

  int NumGauss; /* Number of Gauss quadrature points */

  real8  homogen, AdjustLen;

  real8 *Cweights; /*  interpolation coeffs. */

  real8 *Tkz; /* Chebyshev polynomials */

  int Ktable[343];

  /* Pre-computed matrices coming from SVD calculations */
  real8 *K;

  /* Gauss weights and points */
  real8 *Gausspoints;
  real8 *Gaussweights;

} uniformFMM_t;

extern uniformFMM_t *uniformfmm;

#endif  /* ifdef UNIFORMFMM */


typedef double matrix[3][3];
typedef double vector[3];

typedef struct _fmcell FMCell_t;
typedef struct _fmlayer FMLayer_t;

struct _fmcell {
        int   cellID;     /* ID of the cell in this FM layer */
        int   owningDom;  /* Index of domain owning this cell.  There is  */
                          /* a single unique owner for each cell at each  */
                          /* FM layer.                                    */
        int   domCnt;     /* number of domains intersecting this cell.    */
                          /* this cell.                                   */

        int   *domList;   /* list of all domains that intersect this cell.*/
                          /* This list will only be created for cells at  */
                          /* the most refined FM layer, and is recomputed */
                          /* each time domain boundaries shift.           */

        real8 cellCtr[3]; /* Coordinates of the center of the cell */

        real8 *mpCoeff;   /* Array of multipole expansion coefficients  */
                          /* for the cell.                              */
                          /* For Chebyshev polynomials, the  number of  */
                          /* coefficients is m^3 * 9 where m is the     */
                          /* order of the polynomial                    */
                          /* Otherwise, the  number of coefficients is  */
                          /* (m+3)*(m+2)*(m+1)/6*9 where m is the order */
                          /* multipole expansion.                       */

        real8 *expansionCoeff; /* Array of expansion coefficients           */
                               /* for the cell.                             */
                               /* For Taylor expansions, the  number of     */
                               /* coefficients is (m+3)*(m+2)*(m+1)/6*9     */
                               /* where m is the order of the taylor        */
                               /* expansion.                                */
                               /* For Chebyshev polynomials, the  number of */
                               /* coefficients is m^3 * 6 where m is the    */
                               /* order of the polynomial                   */

#ifdef ESHELBY
        real8 *eshelbympCoeff;
        real8 *eshelbytaylorCoeff;
#endif

        FMCell_t *next;
        FMCell_t *prev;
};


struct _fmlayer {
        int   lDim[3];     /* Number of cells in each dimension at this layer */

        real8 cellSize[3]; /* dimensions of each cell at this layer */
    
        int   ownedCnt;    /* Total number of cells owned by this domain */
                           /* at this layer                              */

        int   ownedMin;    /* Min and max IDs of cells owned by this   */
        int   ownedMax;    /* domain.  Each domain is assumed to own a */
                           /* contiguous range of <ownedCnt> cells.    */
                           /* This information is statically defined   */
                           /* during initialization.  If <ownedCnt> is */
                           /* zero, these values are undetermined.     */

        int   intersectCnt;     /* Number of cells intersecting this domain   */
        int   intersectMin[3];  /* Min and max indices of cells intersecting  */
        int   intersectMax[3];  /* this domain.  These are meaningless except */
                                /* at the most refined FM layer.  At the that */
                                /* layer, the intersection list will be       */
                                /* recomputed each cycle to handle shifting   */
                                /* domain boundries.                          */
        int        *domBuf;

/*
 *      Some lists of domains with which the current domain will 
 *      communicate during the updward and downward FM passes
 */
        int   fmUpPassSendDomCnt;
        int   *fmUpPassSendDomList;

        int   fmUpPassRecvDomCnt;
        int   *fmUpPassRecvDomList;

        int   fmDownPassSendDomCnt;
        int   *fmDownPassSendDomList;

        int   fmDownPassRecvDomCnt;
        int   *fmDownPassRecvDomList;

        int   *cellList; /* List of cellIDs that have been added to */
                         /* the cell table at this FM layer         */

        int   numCells;  /* Number of cells in the cell table at  */
                         /* this FM layer                         */

        FMCell_t *cellTable[CELL_HASH_TABLE_SIZE]; /* Hash table of  */
                         /* FM cell pointers. Contains all cells known  */
                         /* to the current domain at this FM layer      */
};

#if (defined ANISOTROPIC) && (defined TAYLORFMM)
/*
 *  Data structures for Green's function derivative table used with
 *  anisotropic elasticity
 */
typedef struct {
    int     tagmap[11][11][11];
    int     derivTblMaxOrder;

    real8   xScale, yScale, zScale;
    real8   scaleVec[3]; /* {x,y,z}scale in vector form for convenience */
    matrix *Gderiv_map[11][11][11];
} FMAnisoTable_t;
#endif

/*
 *  There are a couple global arrays that are needed by the
 *  Eshelby FMM stuff.  They're defined in FMComm.c, but we'll
 *  add the external declaration here...
 */
extern int    sequence[];
extern real8  multFactor[];

/*
 *      Prototypes for various FM related functions
 */


#ifdef BBFMM
/*
 *        Prototypes needed for the BB FMM implementation
 *  
 */
void bbFMInit(Home_t *home);
void bbfmm_PostM2L(real8 cellSize, real8 *tCoeff);
void bbfmm_S2M(real8 p1[3], real8 p2[3], real8 center[3],
               real8 burg[3], real8 L, real8* S);
void bbfmm_M2M(real8 *r, real8 *Schild, real8 *SS);
void bbfmm_L2L(real8 *r, real8 *F, real8 *Fchild);
void bbfmm_L2P(real8 R[3], real8 L, real8 *F, real8 stress[3][3]);
void bbfmm_M2L(real8 *R, real8 L, real8 *cell_mpCoeff, real8 *FFCoeff);

void ComputeWeights(real8 *Tkz, int *Ktable, real8 *Kweights,
                    real8 *Cweights, int n);
void ComputeKernelSVD(Home_t *home, real8 *Kweights, int n, real8 epsilon,
                      char *Kmat, char *Umat, char *Vmat);

void EvaluateKernelCell(Home_t *home, real8 *field[3], real8 *source[3], int Nf, int Ns,
                        real8 *kernel);

/* For PBC */
void bbFMM_DoTableCorrection(Home_t *home);
#endif

#ifdef UNIFORMFMM
/*
 *       Prototypes needed for the UNOFORM GRID FMM implementation
 *   
 */
void uniformFMInit(Home_t *home);
void uniformfmm_S2M(real8 p1[3], real8 p2[3], real8 center[3],
                    real8 burg[3], real8 L, real8* S);
void uniformfmm_M2M(real8 *r, real8 *Schild, real8 *SS);
void uniformfmm_L2L(real8 *r, real8 *F, real8 *Fchild);
void uniformfmm_L2P(real8 R[3], real8 L, real8 *F, real8 stress[3][3]);
void uniformfmm_M2L(real8 *R, real8 L, real8 *cell_mpCoeff, real8 *FFCoeff);

void FrequencyProduct(int N, real8 *Afre, real8 *xfre, real8 *res);
void ComputeSnUniform(real8 *point[3], int n, int N, real8 *Sn[3]);
void ComputeWeightsUniform(int *Ktable, real8 *Cweights, int n);
void ComputeKernelUniformGrid(Home_t *home, int n, char *Kmat);
void ComputeSn2(real8 *point[3], real8 *Tkz, int n, int N, real8 *Sn[3]);
/* For PBC */
void uniformFMM_DoTableCorrection(Home_t *home);
#endif



// Functions common to both the BBFMM and the UNIFORMFMM methods.
#if (defined BBFMM) || (defined UNIFORMFMM)
void DerivativeOfGreensFunction(Home_t *home,
				real8 p1x, real8 p1y, real8 p1z,
				real8 p2x, real8 p2y, real8 p2z,
				real8 dGdx[3][6]);

void ComputeTk(real8 x, int n, real8 *vec);
void ComputeSn(real8 *point[3], real8 *Tkz, int n, int N, real8 *Sn[3]);
void EvaluateKernel(Home_t *home,
		    real8 p1x, real8 p1y, real8 p1z,
		    real8 p2x, real8 p2y, real8 p2z, real8 *K);
void Stanford_CorrectionTableInit(Home_t *home);
void Stanford_CreateTableCorrection(Home_t *home, int n, int numLevels, int pbc[3]);
void Stanford_DoTableCorrection(Home_t *home, int n, real8 homogen);
void Stanford_FreeCorrectionTable(void);
void Stanford_MeanStressCorrection(Home_t *home, int n, real8* Tkz);


void DerivativeOfGreensFunction(Home_t *home,
                                real8 x, real8 y, real8 z,
                                real8 px, real8 py, real8 pz,
                                int qMax, real8 dGdx[3][6]);


#endif

#if defined (BBFMM) || defined (UNIFORMFMM)
extern "C" {

/*
 *      BLAS and LAPACK routines used by BB FMM
 */

// LAPACK's eigen routine
extern void dgeev_(char *jobvl, char *jobvr, int *n, real8 *A, int *lda,
                   real8 *wr, real8 *wi, real8 *vl, int *ldvl, real8 *vr,
                   int *ldvr, real8 *work, int *lwork, int *info);

// LAPACK's SVD routine
extern void dgesvd_(char *jobu, char *jobvt, int *m, int *n, real8 *A,
                    int *lda, real8 *S, real8 *U, int *ldu, real8 *VT,
                    int *ldvt, real8 *work, int *lwork, int *info);

// BLAS matrix-matrix multiply
extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
                   real8 *alpha, real8 *A, int *lda, real8 *B,
                   int *ldb, real8 *beta, real8 *C, int *ldc);

// BLAS matrix-vector multiply
extern void dgemv_(char *trans, int *m, int *n, real8 *alpha, real8 *A,
                   int *lda, real8 *x, int *incx, real8 *beta,
                   real8 *y, int *incy);

// BLAS dot product
extern real8 ddot_(int *n, real8 *dx, int *incx, real8 *dy, int *incy);

// BLAS daxpy
extern real8 daxpy_(int *n, real8 *da, real8 *dx, int *incx, real8 *dy,
                     int *incy);

// BLAS vector copy
extern real8 dcopy_(int *n, real8 *dx, int *incx, real8 *dy, int *incy);

}  // end extern "C" 
#endif  /* ifdef _BBFMM */


#ifdef UNIFORMFMM
// BLAS dot product
extern "C" {
extern real8 ddot_(int *n, double *dx, int *incx, double *dy, int *incy);
}

extern void FFTBigToSmall(real8 cellSize, real8 *expansionCoeff);
extern void FFTSmallToBig(real8 *mpCoeff);

extern void uniformfmm_CorrectionTableInit(Home_t *home);
extern void uniformfmm_CreateTableCorrection(Home_t *home, int numLevels,
        int pbc[3]);
extern void uniformfmm_DoTableCorrection(Home_t *home);
extern void uniformfmm_FreeCorrectionTable(void);
#endif  /* ifdef UNIFORMFMM */


#ifdef TAYLORFMM
void Taylor_CorrectionTableInit(Home_t *home);
void Taylor_DoTableCorrection(Home_t *home);
void Taylor_CreateCorrectionTable(Home_t *home, int numLevels, int pbc[3],
                                  int mpi_np, int mpi_pid);
void Taylor_FreeCorrectionTable(void);
void Taylor_MeanStressCorrection(Home_t *home);
#endif

//extern void      fmsigma                (real8 mu, real8 nu, int norder, real8 *R, real8 *Eeta, real8 sigmatot[3][3]);
extern void      CorrectionTableInit    (Home_t *home);
extern void      CreateCorrectionTable  (Home_t *home, int numLevels, int pbc[3], int mpi_np, int mpi_pid);

#ifdef CALCENERGY
extern void      CreateEnergyCorrectionTable  (Home_t *home, int numLevels, int pbc[3], int mpi_np, int mpi_pid);
#endif

extern void      DecodeFMCellIndex      (int dim[3], int cellID, int *x, int *y, int *z);
extern void      DoTableCorrection      (Home_t *home);
extern int       EncodeFMCellIndex      (int *dim, int x, int y, int z);
extern void      EvalTaylor             (int norder, real8 *r, real8 *alpha, real8 sigma[3][3]);
extern void      FMCommUpPass           (Home_t *home, int layerID);
extern void      FMCommDownPass         (Home_t *home, int layerID);
extern void      FMUnpackTaylorCoeff    (Home_t *home, int numCells, real8 *buf);
extern void      FMPackTaylorCoeff      (Home_t *home, int numCells, int *cellList, real8 *buf, int bufSize);
extern void      FMDistTaylorExp        (Home_t *home);
extern void      FMFree                 (Home_t *home);
extern void      FMInit                 (Home_t *home);
extern FMCell_t *LookupFMCell           (FMCell_t *cellTable[], int cellID);
extern void      GetCellDomainList      (Home_t *home, int cellIndex, int *domCount, int **domList);
extern void      MkTaylor               (Home_t *home, real8 MU, real8 NU, int norder, int uorder, int maxorder, real8 *r, real8 *eta, real8 *alpha);
extern void      MkTaylorIsotropic      (              real8 mu, real8 nu, int norder, int uorder, int maxorder, real8 *r, real8 *eta, real8 *alpha);
extern void      FMSetRemoteStress      (Home_t *home);
extern void      FMSetTaylorExpansions  (Home_t *home);
extern void      FMShift                (int norder, real8 *r, real8 *eta, real8 *neta);
extern void      FMSigma2               (real8 mu, real8 nu, int norder, real8 Eeta[], matrix sigmatot, real8 pcterms[], int pcpows[][3]);
extern void      FMSigma2NonRational    (real8 mu, real8 nu, int norder, real8 Eeta[], matrix sigmatot, real8 pcterms[], int pcpows[][3]);
extern void      FMSigma2Rational       (real8 mu, real8 nu, int norder, real8 Eeta[], matrix sigmatot, real8 pcterms[], int pcpows[][3]);
extern void      FMFindNearNbrs         (Home_t *home, int layerID, int *inMin, int *inMax, int *outMin, int *outMax, int trimOverlap);
extern void      GaussQuadCoeffHalf     (int numPoints, real8 *positions, real8 *weights);
extern real8     ipow                   (real8 x, int n);
extern int       gcd0                   (int a1, int b1);
extern void      IDToFMCellIndices      (int id, int cellsPerDim[3], int indices[3]);
extern void      makeeta                (int norder, real8 ksi0[3], real8 ksi[3], real8 b[3], real8 *Eeta);
extern void      makeftabs              (real8 *fact, real8 *ifact, real8 *dfact);
extern void      makeqtab               (real8 qtab[NMAX+1][NMAX+1]);
extern void      MeanStressCorrection   (Home_t *home);
extern void      RemoteForceOneSeg      (Home_t *home, Node_t *node1, Node_t *node2, real8 f1[3], real8 f2[3]);
extern void      SegForceFromTaylorExp  (Home_t *home, int cellID, real8 *positions, real8 *weights, real8 *p1, real8 *p2, real8 *burg, real8 p1f[3], real8 p2f[3]);
extern void      TaylorShift            (int norder, real8 *r, real8 *alpha, real8 *beta);
extern void      FMPrintMPCoeff         (Home_t *home, int layerID);
extern void      FMPrint                (Home_t *home);

#ifdef ESHELBY
/*
 *      Prototypes for the Eshelby FM for remote forces from inclusions
 */
extern void      EshelbyBuildMatrix     (real8 *drList, int taylorOrder, int multipoleOrder, real8 *matrix, int mRows, int mCols);
extern void      EshelbyDerofR          (real8 R[3], int order, real8 *drList);
extern void      EshelbyInitFMArray     ();
extern void      EshelbyMakeEta         (real8 Ev, real8 Ro, real8 position[3], real8 center[3], int order, real8 *etav);
extern void      EshelbyMakeTaylor      (real8 R[3], real8 *eta, int multipoleOrder, int taylorOrder, real8 MU, real8 NU, real8 *taylorSeries);
extern void      EshelbyEvalTaylor      (real8 *taylorSeries, real8 taylorCenter[3], real8 position[3], int taylorOrder, real8 stressMatrix[3][3]);
extern void      EshelbyTranslateEta    (real8 *eta, real8 newCenter[3], real8 oldCenter[3], int multipoleOrder, real8 *ep);
extern void      EshelbyTranslateTaylor (real8 *taylorSeries, real8 oldCenter[3], real8 newCenter[3], int taylorOrder, real8 *tPrime);
#endif

/*
 *      Prototype any support functions needed only with anisotropic FMM
 */
#if (defined ANISOTROPIC) && (defined TAYLORFMM)
extern void      FMAnisotropicInit      (Home_t *home);
extern void      G2sigma_deriv          (int norder, int ndx, int ndy, int ndz, real8 Gderiv[][3][3], real8 eta[], real8 C[3][3][3][3], real8 sigma[3][3]);
extern void      ReadDerivTableMetadata (FILE *fpTable, char *fileName, real8 C66[6][6], real8 *xScale, real8 *yScale, real8 *zScale, int *derivTblMaxOrder);
extern void      ReadDerivTable         (Home_t *home, FILE *fpTable);
#endif

#ifdef CALCENERGY
extern real8     FMMEnergy              (Home_t *home, int norder, real8 R[3], real8 *eta1, real8* eta2);
#endif


void SegmentToMoment(Home_t *home, real8 *cellSize, real8 *cellCenter,
                     real8 *vec1, real8 *vec2, real8 *burg,
                     real8 *moments);

void MomentToMoment(Home_t *home, real8 *r, real8 *eta, real8 *neta);
void LocalToLocal(Home_t *home, real8 *r, real8 *alpha, real8 *beta);
void MomentToLocal(Home_t *home, real8 *cellSize, real8 *r, real8 *eta,
                   real8 *alpha);
void LocalToPoint(Home_t *home, real8 *cellSize, real8 *r, real8 *alpha,
                  real8 sigma[3][3]);

#endif  /* FM_h */
