/**************************************************************************
 *
 *      Module:  BlackBoxFMMPrecompeFunction.c
 *               Functions specific to the FMM using Black Box FMM.
 *               Contains functions to do pre calculations. 
 *
 *      Includes functions:
 *
 *
 *      Note: The routines in this file comes from the code uniformfmm of Eric Darve
 *            modified by Chao Chen so that it can handle dislocation segments
 *            and compute the stress field.
 *
 *************************************************************************/
#define _ISOC99_SOURCE
#include <math.h>
#include "Home.h"

#ifdef UNIFORMFMM
#include "FM.h"
#include <fftw/fftw3.h>

uniformFMM_t *uniformfmm;


/*
 * Compute the entry-wise product of two frequencies from rfftw
 * Note 'res' has been initialized
 */
void FrequencyProduct(int N, real8 *Afre, real8 *xfre, real8 *res) 
{
  int i;
  res[0] += Afre[0]*xfre[0];
  for (i=1; i<N; i++) {
    if (i<(N+1)/2)
      res[i] += Afre[i]*xfre[i] + Afre[N-i]*xfre[N-i];
    else
      res[i] += Afre[N-i]*xfre[i] - Afre[i]*xfre[N-i]; 
  }
  if (N%2 == 0)
    res[N/2] += Afre[N/2]*xfre[N/2];
}


#if 0
/*
 * Function: ComputeSnUniform
 * ------------------------------------------------------------------
 * Computes S_n(x_m,x_i) for all Uniform grid node-point pairs using 
 * Clenshaw's recurrence relation.
 */
void ComputeSnUniform(real8 *point[3], int n, int N, real8 *Sn[3]) 
{
  int i, j, k, m;
  real8 vec[n], d[n+2], x, num, denom;
  

  for (m=0;m<n;m++) {
    
    k = m*N;
    for (i=0;i<N;i++) {
      x = point[0][i];
      num = 1.;
      denom = 1.;
      for (j=0;j<n;j++) {
	if(m!=j) {
	  num *= x - (1-2*(real8)j/(real8)(n-1));
	  denom *= 2*((real8)m-(real8)j)/(real8)(n-1);
	}
      }
      Sn[0][k+i] = num/denom;
      
      x = point[1][i];
      num = 1.;
      denom = 1.;
      for (j=0;j<n;j++) {
	if(m!=j) {
	  num *= x - (1-2*(real8)j/(real8)(n-1));
	  denom *= 2*((real8)m-(real8)j)/(real8)(n-1);
	}
      }
      Sn[1][k+i] = num/denom;
      
      x = point[2][i];
      num = 1.;
      denom = 1.;
      for (j=0;j<n;j++) {
	if(m!=j) {
	  num *= x - (1-2*(real8)j/(real8)(n-1));
	  denom *= 2*((real8)m-(real8)j)/(real8)(n-1);
	}
      }
      Sn[2][k+i] = num/denom;
    }
  }

}

#endif

static  double LagrangeWeight(int n, double x, int m) 
{
  int j;
  double num = 1., denom = 1.;
  for (j=0;j<n;j++) {
    if(m!=j) {
      num   *= x - (1-2*(double)j/(double)(n-1));
      denom *= -2*((double)m-(double)j)/(double)(n-1);
    }
  }
  return num/denom;
}


/* NEW FUNCTION */
void ComputeSnUniform(real8 *point[3], int n, int N, real8 *Sn[3]) 
{
  int i, m, k;

  for (m=0;m<n;m++) {
    k = m*N;
    for (i=0;i<N;i++) {
      Sn[0][k+i] = LagrangeWeight(n, point[0][i], m);
      Sn[1][k+i] = LagrangeWeight(n, point[1][i], m);
      Sn[2][k+i] = LagrangeWeight(n, point[2][i], m);
    }
  }
}


/*
 * Function: ComputeWeightsUniform
 * ------------------------------------------------------------------
 * Computes the weights for the Uniform grid nodes of all children cells
 * (identical for all cells and all levels so just compute once and
 * store in memory) and set up the lookup table.
 */
void ComputeWeightsUniform(int *Ktable, real8 *Cweights, int n) 
{
  int i, k, m, l1, l2, l3, count, ncell, ninteract, n2, n3;
  real8 tmp1, tmp2;
  real8 vtmp[3];
  real8 pi = M_PI; 
  real8 AdjustLen,scale;
  real8 *nodes;
  real8 *fieldt[3], *Sn[3];

  nodes = (real8 *)malloc(n * sizeof(real8));
  n2 = n*n;
  n3 = n2*n;

  AdjustLen = uniformfmm->AdjustLen;

  // Number of child Uniform grid nodes
  int Nc = 2*n3; 

  for (i=0;i<3;i++)
    fieldt[i] = (real8 *)malloc(Nc * sizeof(real8));                 

  // Uniform grid-transformed coordinates
  for (i=0;i<3;i++)
    Sn[i]     = (real8 *)malloc(n * Nc * sizeof(real8));

  // Define Tkz similar to Chebychev for PBC conditions
  uniformfmm->Tkz      = (real8 *)malloc(n2 * sizeof(real8)); 

  real8 vec[n];

  // Compute the n Chebychev nodes of T_n(x)
  for (m=0;m<n;m++)
    nodes[m] = cos(M_PI*((real8)m+0.5)/(real8)n);
  
  // Evaluate the Chebychev polynomials of degree 0 to n-1 at the nodes
  for (m=0;m<n;m++) {
    ComputeTk(nodes[m],n,vec);
    i = m*n;
    for (k=0;k<n;k++)
      {
	uniformfmm->Tkz[i+k] = vec[k];
      }
  }
  
  // Initialize lookup table
  for (i=0;i<343;i++)
    Ktable[i] = -1;
  
  // Create lookup table
  ncell = 0;
  ninteract = 0;
  for (l1=-3;l1<4;l1++) {
    for (l2=-3;l2<4;l2++) {
      for (l3=-3;l3<4;l3++) {
	if (abs(l1) > 1 || abs(l2) > 1 || abs(l3) > 1) {
	  Ktable[ncell] = ninteract;
	  ninteract++;
	}
	ncell++;
      }
    }
  }	
  
  // Compute the n uniform grid nodes of T_n(x)
  for (m=0;m<n;m++)
    nodes[m] = 1 - 2*(real8)m/((real8)(n-1));
  
  // Compute the weights for the kernel matrix K : Kweights were 1 : removed !
  
  // Map all Uniform grid nodes from the children cells to the parent
  k = 0;
  scale = 1 / (1 + AdjustLen);
  for (i=0;i<2;i++) {
    // Determine the mapping function for the specific child cell
    vtmp[0] = -0.5;//-1;
    vtmp[1] = -0.5;//-1;
    
    if (i == 0)
      vtmp[2] =-0.5;//-1;
    else
      vtmp[2] = 0.5;// 1;
    
    for (l1=0;l1<n;l1++) {
      for (l2=0;l2<n;l2++) {
	for (l3=0;l3<n;l3++) {
	  fieldt[0][k] = 0.5*nodes[l1] + vtmp[0]*scale;
	  fieldt[1][k] = 0.5*nodes[l2] + vtmp[1]*scale;
	  fieldt[2][k] = 0.5*nodes[l3] + vtmp[2]*scale;
	  k++;
	}
      }
    }
  }
    
  // Compute Sc, the mapping function for the field points
  ComputeSnUniform(fieldt,n,Nc,Sn);
  
  // Extract out the Uniform grid weights
  count = 0;
  for (l1=0;l1<n;l1++) {
    k = l1*Nc;
    for (l2=0;l2<n;l2++) {
      Cweights[count] = Sn[2][k+l2];
      count++;
    }
  }
  for (l1=0;l1<n;l1++) {
    k = l1*Nc;
    for (l2=0;l2<n;l2++) {
      Cweights[count] = Sn[2][k+n3+l2];
      count++;
    }
  }
  free(nodes);
  for (i=0;i<3;i++) 
    {
      free(fieldt[i]);
      free(Sn[i]);
    }
}


/*
 * Given n node positions and returns corresponding field and source
 * positions with respect to the index
 */
static void GetPosition(int n, int idx, real8 *fieldpos, real8 *sourcepos, real8 *nodepos) 
{
  
  if (idx < n) {
    *fieldpos  = nodepos[n-1]/2;
    *sourcepos = nodepos[idx]/2;
  } else {
    *fieldpos  = nodepos[2*(n-1)-idx]/2;
    *sourcepos = nodepos[n-1]/2;
  }
}


/*
 * Function: ComputeKernelUniformGrid
 * ---------------------------------------------------------------------
 * Computes the kernel for 316(2n-1)^3 interactions between Uniform
 * Grid nodes. Does not compute SVD.
 */
void ComputeKernelUniformGrid(Home_t *home, int n, char *Kmat)
{

   int i, m, k1, k2, k3, l1, l2, l3;
   int dof2 = 6 * 9;

    // Total size    
    // int dof2n6 = dof2 * (2*n-1)*(2*n-1)*(2*n-1); 

    real8 nodes[n], kernel[dof2];
    real8 fieldpos[3], sourcepos[3], scenter[3];
    
    
    // Compute Uniform grid nodes of T_n(x)
    for (m=0; m<n; m++)
      nodes[m] = 1 - 2*(real8)m/((real8)(n-1));


    // Create FFT plan
    int vecSize = 2*n-1, reducedMatSize = pow(vecSize, 3);
    int M2LSize = dof2 *reducedMatSize;

    real8 *MatM2L  = (real8 *)(fftw_malloc(M2LSize*sizeof(real8)));
    real8 *freqMat = (real8 *)(fftw_malloc(316*M2LSize*sizeof(real8)));
    fftw_plan p[316*dof2];
    
    // Compute the kernel values for interactions with all 316 cells
    int countM2L=0, count, count1, countPlan=0;
    int shift1, shift2, shiftGlo, shiftLoc;
    int reducedMatSizeDofs = reducedMatSize * 9;
    int f, s;


    for (i=0;i<M2LSize;i++) 
      MatM2L[i] = 0.0;

    for (i=0;i<316*M2LSize;i++) 
	freqMat[i] = 0.0;

    
    for (k1=-3;k1<4;k1++) {
      scenter[0] = (real8)k1;
      for (k2=-3;k2<4;k2++) {
	scenter[1] = (real8)k2;
	for (k3=-3;k3<4;k3++) {
	  scenter[2] = (real8)k3;
	  if (abs(k1) > 1 || abs(k2) > 1 || abs(k3) > 1) {
	    
	    for (count=0, l1=0; l1<vecSize; l1++) {
	      GetPosition(n, l1, &fieldpos[0], &sourcepos[0], nodes);
	      sourcepos[0] += scenter[0];
	      for (l2=0; l2<vecSize; l2++) {
		GetPosition(n, l2, &fieldpos[1], &sourcepos[1], nodes);
		sourcepos[1] += scenter[1];
		for (l3=0; l3<vecSize; l3++, count++) {
		  GetPosition(n, l3, &fieldpos[2], &sourcepos[2], nodes);
		  sourcepos[2] += scenter[2];
		  
		  EvaluateKernel(home,
				 fieldpos[0],fieldpos[1],fieldpos[2],
				 sourcepos[0],sourcepos[1],sourcepos[2],
				 kernel);
		  
		  // kernel[f x s] in column major
		  // MatM2L[(2n-1)^3 x s x f] in column
		  // major
		  count1 = 0;
		  shift1 = count;
		  for (s=0; s < 9; s++) {
		    for(f=0, shift2=0; f < 6; f++) {
		      MatM2L[shift1 + shift2] =  kernel[count1++];
		      shift2 += reducedMatSizeDofs;
		    }
		    shift1 += reducedMatSize;
		  }
		}
	      }
	    }
	    
	    // FFT
	    shiftGlo = countM2L *M2LSize;
	    
	    for (i=0, shiftLoc=0; i<dof2; i++) {
	      p[ countPlan ] = fftw_plan_r2r_1d(reducedMatSize, MatM2L + shiftLoc, freqMat + shiftGlo + shiftLoc, FFTW_R2HC, FFTW_FLAG);
	      fftw_execute( p[ countPlan ] );
	      countPlan++;
	      
	      shiftLoc += reducedMatSize;
	    }
	    
	    countM2L++;
	  }
	}
      }
    }

    FILE *ptr_file;
    ptr_file = fopen(Kmat, "wb");
    fwrite(freqMat, sizeof(real8),316 *M2LSize, ptr_file);
    fclose(ptr_file);

   fftw_free(MatM2L);
   fftw_free(freqMat);
   for (i=0; i<countPlan; i++)
     fftw_destroy_plan( p[i] );
}


/* Initialize the uniformfmm structure */
/* n : input parameter read in the input file */
void uniformFMInit(Home_t *home)
{
  int i,n,n2,n3;
  char *Kmat;
  Param_t *param;
  FILE *f;
  int domainID = home->myDomain;
  
  if (domainID == 0) printf("\n\nInitialization of uniform grid method\n");
  param = home->param;

#ifdef ANISOTROPIC
  AnisotropicVars_t *anisoVars;
  anisoVars = &home->anisoVars;
#endif

  n = param->fmExpansionOrder;
  n2 = n * n;
  n3 = n2 * n;

  Kmat = param->KMatFileName;

  if (domainID == 0) printf("Kmat name %s\n",Kmat);

  uniformfmm = (uniformFMM_t *)calloc(1, sizeof(uniformFMM_t));

  /* Number of Uniform grid nodes */
  uniformfmm->n = n;
  uniformfmm->n3 = n3;
  uniformfmm->pown3 = (2*n-1)*(2*n-1)*(2*n-1);

  if (domainID == 0) printf("Number of uniform grid points %d, n*n*n = %d, (2*n-1)*(2*n-1)*(2*n-1)=%d \n",
	 n, n3, uniformfmm->pown3);

  /* Kernel is homogeneous of order 2 */
  uniformfmm->homogen = 2.0; 


  uniformfmm->AdjustLen = param->fmAdjustLen;
  if (domainID == 0) printf("Cell size length adjustment is %f\n",uniformfmm->AdjustLen);

  /* Gauss quadrature  */
  uniformfmm->NumGauss = 5;
  if (domainID == 0) printf("Number of Gauss points : %d\n",uniformfmm->NumGauss);
  uniformfmm->Gausspoints =  (real8 *)malloc(uniformfmm->NumGauss * sizeof(real8));
  uniformfmm->Gaussweights = (real8 *)malloc(uniformfmm->NumGauss * sizeof(real8));
  GaussQuadraturePointsAndWgts(uniformfmm->NumGauss, uniformfmm->Gausspoints, uniformfmm->Gaussweights);
  
  int dofn3_s = 9 * n3;
  int dofn3_f = 6 * n3;

  // Uniform grid interpolation coefficients: Sn (page 8715)
  uniformfmm->Cweights = (real8 *)malloc(2 * n2 * sizeof(real8));

  // Compute the uniform weights and sets up the lookup table
  ComputeWeightsUniform(uniformfmm->Ktable,uniformfmm->Cweights,n);
  
  // Precompute the kernel interaction matrix (if necessary)
  if ((f = fopen(Kmat,"rb")) == NULL) {
    
    if (domainID == 0) printf("Computing Uniform grid kernel K\n");
    if (domainID == 0) {
      ComputeKernelUniformGrid(home, n, Kmat);
      if (domainID == 0) printf("Matrix K computed\n");
    }
  }
  else
    if (domainID == 0) printf("Matrix K already computed\n");
  
  /*
   *  Only domain zero creates the matrix files on disk, so make sure
   *  all processes wait until the files exist
   */
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* 
   * Initialize arrays, K is the compressed M2L operator C^{(i)}
   */
  
  int Ksize;
  Ksize = 316*uniformfmm->pown3*6*9;

  uniformfmm->K  = (real8 *)malloc(Ksize * sizeof(real8));
  
  if (uniformfmm->K == NULL)
    Fatal("K's size is too large\n");
  
  // Read kernel interaction matrix K
  FILE *ptr_file;

 // Read in kernel interaction matrix K
  ptr_file = fopen(Kmat,"rb");

  if(ptr_file == NULL)
    Fatal("Cannot read Kmat = %s\n",Kmat);

  fread(uniformfmm->K, sizeof(real8), Ksize, ptr_file);

 if (domainID == 0) printf("Initialization uniform grid method done\n\n");
}
	

void uniformfmmFinal(uniformFMM_t *uniformfmm)
{
  free(uniformfmm);
}
#endif /* ifdef UNIFORMFMM */
