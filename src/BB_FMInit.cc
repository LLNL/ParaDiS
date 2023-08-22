/**************************************************************************
 *
 *      Module:  BlackBoxFMMPrecompeFunction.c
 *               Functions specific to the FMM using Black Box FMM.
 *               Contains functions to do pre calculations. 
 *
 *      Includes functions:
 *
 *
 *      Note: The routines in this file comes from the code bbfmm of Eric Darve
 *            modified by Chao Chen so that it can handle dislocation segments
 *            and compute the stress field.
 *
 *************************************************************************/
#include "Home.h"

#ifdef BBFMM
#include "FM.h"

bbFMM_t *bbfmm;




/*
 * Function: ComputeWeights
 * ------------------------------------------------------------------
 * Computes the weights for the Chebyshev nodes of all children cells
 * (identical for all cells and all levels so just compute once and
 * store in memory) and set up the lookup table.
 */
void ComputeWeights(real8 *Tkz, int *Ktable, real8 *Kweights, 
		    real8 *Cweights, int n) 
{
  int i, k, m, l1, l2, l3, count, ncell, ninteract;
  real8 tmp1, tmp2;
  real8 vtmp[3];
  real8 AdjustLen,scale;
  real8 pi = M_PI;   
  real8 *nodes, *vec;
  real8 *fieldt[3], *Sn[3];

  nodes = (real8 *)malloc(n * sizeof(real8));
  if (nodes == NULL)
    {
      Fatal("nodes : Out of memory! %d\n",n);
    }

  vec   = (real8 *)malloc(n * sizeof(real8));
  if (vec == NULL)
    {
      Fatal("vec : Out of memory! %d\n",n);
    }

  int n3 = n*n*n;                     // n3 = n^3
  int Nc = 2*n3;                      // Number of child Chebyshev nodes

  AdjustLen = bbfmm->AdjustLen;


  for (i=0;i<3;i++)
    {
      fieldt[i] = (real8 *)malloc(Nc * sizeof(real8));                 

      if (fieldt[i] == NULL)
	{
	  Fatal("fieldt[%d] : Out of memory! %d\n",i, Nc);
	}
    }

  // Chebyshev-transformed coordinates
  for (i=0;i<3;i++)
    {
      Sn[i]     = (real8 *)malloc(n * Nc * sizeof(real8));

      if (Sn[i] == NULL)
	{
	  Fatal("Sn[%d] : Out of memory! %d\n",i, n*Nc);
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
  
  // Compute the n Chebyshev nodes of T_n(x)
  for (m=0;m<n;m++)
    nodes[m] = cos(pi*((real8)m+0.5)/(real8)n);
  
  // Evaluate the Chebyshev polynomials of degree 0 to n-1 at the nodes
  for (m=0;m<n;m++) {
    ComputeTk(nodes[m],n,vec);
    i = m*n;
    for (k=0;k<n;k++)
      Tkz[i+k] = vec[k];
  }
  
  // Compute the weights for the kernel matrix K
  count = 0;
  for (l1=0;l1<n;l1++) {
    tmp1 = 1./sqrt(1-nodes[l1]*nodes[l1]);
    for (l2=0;l2<n;l2++) {
      tmp2 = tmp1/sqrt(1-nodes[l2]*nodes[l2]);
      for (l3=0;l3<n;l3++) {
	Kweights[count] = tmp2/sqrt(1-nodes[l3]*nodes[l3]);
	count++;
      }
    }
  }
  
  // Map all Chebyshev nodes from the children cells to the parent
  k = 0;
  scale = 1/(1 + AdjustLen);
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

  // Compute Sn, the mapping function for the field points
  ComputeSn(fieldt,Tkz,n,Nc,Sn); 
  
  // Extract out the Chebyshev weights
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
  free(vec);
  for (i=0;i<3;i++) 
    {
      free(fieldt[i]);
      free(Sn[i]);
    }
}


/*
 * Function: ComputeKernelSVD
 * ---------------------------------------------------------------------
 * Computes the kernel for 316n^6 interactions between Chebyshev nodes
 * and then computes the SVD of the kernel matrix.
 */
void ComputeKernelSVD(Home_t *home, real8 *Kweights, int n, real8 epsilon, 
		      char *Kmat, char *Umat, char *Vmat)
{
  int i, j, l, m, k1, k2, k3, l1, l2, l3, m1, z;
  int count, count1, count2, count3;
  real8 scenter[3], vtmp[2];
  real8 sweight, fweight,scale;
  
  int n3 = n*n*n;            // n3 = n^3
  int dofn3_s = 9 * n3;
  int dofn3_f = 6 * n3;
  int dof2n6 = dofn3_s * dofn3_f; // Total size
  int Sigma_size;
  int cutoff[2];
  
  real8 *K0, *U0, *Sigma, *VT0;
  real8 *nodes, *kernel, *work;
  real8 *fieldpos[3], *sourcepos[3];

  K0 = (real8 *)malloc(316 * dof2n6 * sizeof(real8));
  kernel = (real8 *)malloc(dof2n6 * sizeof(real8));

  if (K0 == NULL)
    {
      Fatal("K0 : Out of memory! %d\n", 316*dof2n6);
    }

  if (kernel == NULL)
    {
      Fatal("kernel : Out of memory! %d\n", dof2n6);
    }

  for (z = 0; z < 3; ++z)
    {
      fieldpos[z]  = (real8 *)malloc(n3 * sizeof(real8));

      if (fieldpos[z] == NULL)
	{
	  Fatal("fieldpos[%d] : Out of memory! %d\n",z, n3);
	}


      sourcepos[z] = (real8 *)malloc(n3 * sizeof(real8));	

      if (sourcepos[z] == NULL)
	{
	  Fatal("sourcepos[%d] : Out of memory! %d\n",z, n3);
	}

    }

  // Compute Chebyshev nodes of T_n(x)
  nodes = (real8 *)malloc(n * sizeof(real8));

  scale = 1 + bbfmm->AdjustLen;
  for (m=0;m<n;m++)
    {
      nodes[m] = cos(M_PI*((real8)m+0.5)/(real8)n) * scale;
    }
  // Compute the locations of the field points in a unit cube
  count = 0;
  for (l1=0;l1<n;l1++) {
    vtmp[0] = 0.5*nodes[l1];
    for (l2=0;l2<n;l2++) {
      vtmp[1] = 0.5*nodes[l2];
      for (l3=0;l3<n;l3++) {
	fieldpos[0][count] = vtmp[0];
	fieldpos[1][count] = vtmp[1];
	fieldpos[2][count] = 0.5*nodes[l3];
	count++;
      }
    }
  }

  // Compute the kernel values for interactions with all 316 cells
  // Simplify the loops and use anti-symmetry property of the kernel
  count = 0;
  int countM2L=0;
  int col, row, idx1, idx2;

  for (k1=-3;k1<4;k1++) 
    {
      scenter[0] = (real8)k1;
      for (k2=-3;k2<4;k2++) 
	{
	  scenter[1] = (real8)k2;
	  for (k3=-3;k3<4;k3++) 
	    {
	      scenter[2] = (real8)k3;
	      if (abs(k1) > 1 || abs(k2) > 1 || abs(k3) > 1) 
		{
		  if (countM2L < 158) 
		    {
		      // Compute the locations of the source points in the cell
		      for (count1=0, l1=0;l1<n;l1++) {
			vtmp[0] = scenter[0] + 0.5*nodes[l1];
			for (l2=0;l2<n;l2++) {
			  vtmp[1] = scenter[1] + 0.5*nodes[l2];
			  for (l3=0; l3<n; l3++, count1++) {
			    sourcepos[0][count1] = vtmp[0];
			    sourcepos[1][count1] = vtmp[1];
			    sourcepos[2][count1] = scenter[2] + 0.5*nodes[l3];
			  }
			}
		      }
		      // Compute the kernel at each of the field Chebyshev nodes	  
		      EvaluateKernelCell(home, fieldpos, sourcepos, n3, n3, kernel);
	    
		      // Copy the kernel values to the appropriate location
		      count1=0;
		      count2=0;
		      for (l1=0; l1<n3; l1++, count2++) 
			{
			  sweight = Kweights[count2];
			  for (l=0;l<9;l++) 
			    { 
			      for (count3=0, m1=0; m1<n3; m1++, count3++) 
				{
				  fweight = Kweights[count3];
				  for (m=0; m<6; m++, count++, count1++) 
				    { 
				      K0[count] = kernel[count1]/(sweight*fweight);
				    }
				}
			    }
			}
		      countM2L++;
		    }
		  else 
		    { // For anti-symmetric kernel
		      for (col=0; col<dofn3_s; col++)
			for (row=0; row<dofn3_f; row++) 
			  {
			    idx1 = (315-countM2L)*dof2n6 + col*dofn3_f + row;
			    idx2 = countM2L*dof2n6 + (row/6*9 + col%9)*dofn3_f
			      + (col/9*6 + row%6);
			    K0[idx2] = - K0[idx1];
			}
		      countM2L++;
		    }
		}
	    }
	}
    }
  

  free(kernel);
 
  for (z = 0; z < 3; ++z)
    {
      free(fieldpos[z]);
      free(sourcepos[z]);
    } 

  // Extract the submatrix for each of the 316 cells
  // 316 M2L operators
  real8* Kcell[316];
  for (z = 0; z < 316; ++z)
    {
      Kcell[z] = (real8 *) malloc(dof2n6 * sizeof(real8));
  
      if (Kcell[z] == NULL)
	{
	  Fatal("Kcell[%d] : Out of memory! %d\n", z,dof2n6);
	}
    }
 
  count = 0;
  for (i=0;i<316;i++) {
    for (j=0;j<dof2n6;j++) {
      Kcell[i][j] = K0[count];
      count++;
    }
  }
  
  /****
   * Compute the SVD of K_fat
   ****/
  
  // Compute the SVD of K0
  char save[]="S", nosave[]="N";
  int nosavedim=1;
  int info, lwork;
  int cols_s = 316*dofn3_s;
  
  /* See dgesvd documentation:
   *          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
   *             - PATH 1  (M much larger than N, JOBU='N')  - our case for K_thin
   *             - PATH 1t (N much larger than M, JOBVT='N') - our case for K_fat
   */
  int max_dof = 9;
  lwork = 5*max_dof*n3; // Change
  work  = (real8 *) malloc(lwork * sizeof(real8));

  if (work == NULL)
    {
      Fatal("work : Out of memory! %d\n",lwork);
    }  

  int U0size;
  Sigma_size = dofn3_f;	
  U0size     = dofn3_f * dofn3_f;
  Sigma      = (real8 *)malloc(Sigma_size * sizeof(real8));

  if (Sigma == NULL)
    {
      Fatal("Sigma : Out of memory! %d\n",Sigma_size);
    } 

  U0         = (real8 *)malloc(U0size * sizeof(real8));	

  if (U0 == NULL)
    {
      Fatal("U0 : Out of memory! %d\n",U0size);
    } 

  VT0        = NULL;

  dgesvd_(save,nosave,&dofn3_f,&cols_s,K0,&dofn3_f,Sigma,U0,&dofn3_f,VT0,&nosavedim,
	  work,&lwork,&info);
  
  
  /* 
   * Determine the number of singular values to keep. We use epsilon for this.
   */
  
  real8 sumsqr = 0;
  for (i=Sigma_size-1; i>=0; --i)
    sumsqr += Sigma[i] * Sigma[i];
  
  real8 sum = 0, epsqr = sumsqr * epsilon * epsilon;
  
  cutoff[0] = Sigma_size;
  for (i=Sigma_size-1; i>=0; --i) {
    sum += Sigma[i] * Sigma[i];
    if (sum < epsqr)
      --cutoff[0];
    else
      break;
  }
  
  // Extract the needed columns from U0 and write out to a file
  FILE *ptr_file;
  ptr_file = fopen(Umat,"wb");
  fwrite(&cutoff[0], sizeof(int), 1, ptr_file);
  
  real8 truncatedSize = dofn3_f * cutoff[0];
  fwrite(U0, sizeof(real8), truncatedSize, ptr_file);
  fclose(ptr_file);
  
  free(Sigma); Sigma = NULL;
  
  
  /****
   * Compute the SVD of K_thin
   ****/
  
  // Form K_thin using all of 316 M2L operators stored in Kcell
  count = 0;
  for (j=0;j<dofn3_s;++j) {		
    for (i=0;i<316;i++) {
      for (l=0;l<dofn3_f;l++) {
	K0[count] = Kcell[i][l+j*dofn3_f];
	count++;
      }
    }
  }

  // save = "S"; nosave = "N"
  real8 *U1 = NULL;
  int rows_f = 316*dofn3_f;
  Sigma_size = dofn3_s;	
  Sigma      = (real8 *)malloc(Sigma_size * sizeof(real8));

  if (Sigma == NULL)
    {
      Fatal("Sigma : Out of memory! %d\n",Sigma_size);
    } 

  VT0        = (real8 *)malloc(dofn3_s * dofn3_s * sizeof(real8));

  if (VT0 == NULL)
    {
      Fatal("VT0 : Out of memory! %d\n",dofn3_s * dofn3_s);
    } 

  // Compute the SVD of the K_thin
  dgesvd_(nosave,save,&rows_f,&dofn3_s,K0,&rows_f,Sigma,U1,&nosavedim,VT0,&dofn3_s,
	  work,&lwork,&info);
  
  /* 
   * Determine the number of singular values to keep. We use epsilon for this.
   */

  sumsqr = 0;
  for (i=Sigma_size-1; i>=0; --i) 
    sumsqr += Sigma[i] * Sigma[i];
  
  sum = 0, epsqr = sumsqr * epsilon * epsilon;
  cutoff[1] = Sigma_size;
  for (i=Sigma_size-1; i>=0; --i) {
    sum += Sigma[i] * Sigma[i];
    if (sum < epsqr)
      --cutoff[1];
    else
      break;
  }
  
  // Extract truncated VT[cutoff_s][dofn3_s] from VT0[dofn3_s][dofn3_s] and write out to a file
  ptr_file = fopen(Vmat,"wb");
  fwrite(&cutoff[1], sizeof(int), 1, ptr_file);
  
  for (j=0;j<dofn3_s;j++) { // column
    count1 = j*dofn3_s;
    fwrite(VT0+count1, sizeof(real8), cutoff[1], ptr_file);
  }
  fclose(ptr_file);
  
  /** Computing the compressed kernel using the orthonormal basis U and VT.
   **/
  char *transa, *transb;
  int cutoff2 = cutoff[0] * cutoff[1];	
  real8 alpha=1, beta=0;
  real8 *Ecell, *KV;
  Ecell = (real8 *)malloc(cutoff2 * sizeof(real8));

  if (Ecell == NULL)
    {
      Fatal("Ecell : Out of memory! %d\n",cutoff2);
    } 
	
  KV    = (real8 *)malloc(dofn3_f * cutoff[1] * sizeof(real8));
 
  if (KV == NULL)
    {
      Fatal("KV : Out of memory! %d\n",dofn3_f * cutoff[1]);
    } 

  ptr_file = fopen(Kmat,"wb");
  for (i=0;i<316;i++) {       
    
    /* Compute K V:
     * K  is dofn3_f  x dofn3_s
     * VT is cutoff[1] x dofn3_s
     * V  is dofn3_s  x cutoff[1] 
     * KV is dofn3_f  x cutoff[1]
     * (Notice that VT is a submatrix of VT0)
     */
    transa = "n";
    transb = "t";
    dgemm_(transa, transb, &dofn3_f, &cutoff[1], &dofn3_s, &alpha, 
	   Kcell[i], &dofn3_f, VT0, &dofn3_s, &beta, KV, &dofn3_f);
    
    /* Compute U^T K V:
     * KV is dofn3_f x cutoff[1]
     * U  is dofn3_f x cutoff[0]
     * UT is cutoff[0] x dofn3_f
     * U^T K V is cutoff[0] x cutoff[1]
     * (Notice that U is a submatrix of U0)
     */
    transa = "t";
    transb = "n";
    dgemm_(transa, transb, &cutoff[0], &cutoff[1], &dofn3_f, &alpha, 
	   U0, &dofn3_f, KV, &dofn3_f, &beta, Ecell, &cutoff[0]);
    
    fwrite(Ecell, sizeof(real8), cutoff2, ptr_file);
  }
  fclose(ptr_file);
  
  free(K0);
  free(VT0);
  free(U0);
  free(Sigma);
  free(nodes);

  free(work);
  free(KV);
  free(Ecell);
  
  for (z = 0; z < 316; ++z)
    free(Kcell[z]);
}


/* Initialize the bbfmm structure */
/* n : input parameter read in the input file */
void bbFMInit(Home_t *home)
{
  int i,n,n2,n3;
  int dofn3_s, dofn3_f;
  int domainID = home->myDomain;

  FILE *f;
  real8 epsilon;
  char *Kmat, *Umat, *Vmat;
  Param_t *param;

  if (domainID == 0) printf("\n\nInitialization of Chebyshev method\n");
  param = home->param;


#ifdef ANISOTROPIC
  AnisotropicVars_t *anisoVars;
  anisoVars = &home->anisoVars;
  if (domainID == 0) printf("Anisotropic Elasticity qMax=%d\n", param->anisoHarmonicsNumTermsBase);
#endif

  n = param->fmExpansionOrder;
  n2 = n * n;
  n3 = n2 * n;

  dofn3_s = 9 * n3;
  dofn3_f = 6 * n3;

  Kmat = param->KMatFileName;
  Umat = param->UMatFileName;
  Vmat = param->VMatFileName;

  if (domainID == 0) printf("Kmat name %s\n",Kmat);
  if (domainID == 0) printf("Umat name %s\n",Umat);
  if (domainID == 0) printf("Vmat name %s\n",Vmat);

  epsilon = param->ChebEps;

  //printf("Overriding ChebEps from %f to %f\n",epsilon,1./(1.0*n));
  //epsilon = 1./(1.0*n); 
  //param->ChebEps = 1./(1.0*n); 

  bbfmm = (bbFMM_t *)calloc(1, sizeof(bbFMM_t));


  /* Number of Chebyshev nodes */
  bbfmm->n = n;
  bbfmm->n3 = n3;
  if (domainID == 0) printf("Chebyshev order %d\n",n);

  /* Kernel is homogeneous of order 2 */
  bbfmm->homogen = 2.0; 

  /* Gauss quadrature  */
  bbfmm->NumGauss = 5;
  if (domainID == 0) printf("Number of Gauss points : %d\n",bbfmm->NumGauss);
  bbfmm->Gausspoints =  (real8 *)malloc(bbfmm->NumGauss * sizeof(real8));
  bbfmm->Gaussweights = (real8 *)malloc(bbfmm->NumGauss * sizeof(real8));

  if (bbfmm->Gausspoints == NULL || bbfmm->Gaussweights == NULL)
    {
      Fatal("bbfmm->Gausspoint : Out of memory! %d\n",bbfmm->NumGauss);
    }

  GaussQuadraturePointsAndWgts(bbfmm->NumGauss, bbfmm->Gausspoints, bbfmm->Gaussweights);

  bbfmm->AdjustLen = param->fmAdjustLen;
  if (domainID == 0) printf("Cell size length adjustment is %f\n",bbfmm->AdjustLen);

  bbfmm->ChebEps = epsilon;

  // Omega matrix : weighting matrix for Chebyshev quadrature approximation
  bbfmm->Kweights = (real8 *)malloc(n3 * sizeof(real8));

  if (bbfmm->Kweights == NULL)
    {
      Fatal("bbfmm->Kweights : Out of memory! %d\n",n3);
    }

  
  // Chebyshev interpolation coefficients: Sn (page 8715)
  bbfmm->Cweights = (real8 *)malloc(2 * n2 * sizeof(real8));

  if (bbfmm->Kweights == NULL)
    {
      Fatal("bbfmm->Kweights : Out of memory! %d\n",2*n2);
    }


  // Evaluation of n Chebyshev nodes (of T_n) on 
  // n chebyshev polynomials from T_0 to T_{n-1}
  bbfmm->Tkz      = (real8 *)malloc(n2 * sizeof(real8)); 
  
  if (bbfmm->Tkz == NULL)
    {
      Fatal("bbfmm->Tkz: Out of memory! %d\n",n2);
    }


  // Set up for FMM
  int cutoff[2];

  // Compute the Chebyshev weights and sets up the lookup table
  ComputeWeights(bbfmm->Tkz,bbfmm->Ktable,bbfmm->Kweights,bbfmm->Cweights,n);

  // Precompute the SVD of the kernel interaction matrix (if necessary)
  if ((f = fopen(Kmat,"rb")) == NULL || (f = fopen(Umat,"rb")) == NULL || 
      (f = fopen(Vmat,"rb")) == NULL) {

    if (domainID == 0) printf("Computing matrices K, U, Vt\n");
    if (domainID == 0) {
      ComputeKernelSVD(home, bbfmm->Kweights, n, epsilon, Kmat, Umat, Vmat);
      if (domainID == 0) printf("Matrices computed\n");
    }
  }
  else
    if (domainID == 0) printf("Matrices already computed\n");

/*
 *  Only domain zero creates the matrix files on disk, so make sure
 *  all processes wait until the files exist
 */
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  f = fopen(Umat,"rb");
  fread(&(cutoff[0]), sizeof(int), 1, f);
  fclose(f);
  
  f = fopen(Vmat,"rb");
  fread(&(cutoff[1]), sizeof(int), 1, f);
  fclose(f);	

  /* 
   * Initialize arrays, K is the compressed M2L operator C^{(i)}
   * 
   * U: U^k_r p. 8719; downward pass; field
   * V: S^K_r p. 8718; upward pass; source
   */
  
  int Ksize, Usize, Vsize;
  Ksize = cutoff[0] * (316*cutoff[1]); 
  Usize = cutoff[0] * 6 *n3;
  Vsize = cutoff[1] * 9 *n3;

  bbfmm->K  = (real8 *)malloc(Ksize * sizeof(real8));
  if (bbfmm->K == NULL)
    {
      Fatal("bbfmm->K : Out of memory! %d\n",Ksize);
    }
  
  bbfmm->U  = (real8 *)malloc(Usize * sizeof(real8)); 
  if (bbfmm->U == NULL)
    {
      Fatal("bbfmm->U : Out of memory! %d\n",Usize);
    }

  bbfmm->VT = (real8 *)malloc(Vsize * sizeof(real8));
  if (bbfmm->VT == NULL)
    {
      Fatal("bbfmm->V : Out of memory! %d\n",Vsize);
    }

  // Read kernel interaction matrix K and singular vectors U and VT
  FILE *ptr_file;

  // Read in kernel interaction matrix K
  ptr_file = fopen(Kmat,"rb");
  fread(bbfmm->K, sizeof(real8), Ksize, ptr_file);
  fclose(ptr_file);

  // Read in matrix of singular vectors U
  ptr_file = fopen(Umat,"rb");
  fread(&i, sizeof(int), 1, ptr_file);
  fread(bbfmm->U, sizeof(real8), Usize, ptr_file);
  fclose(ptr_file);
  
  // Read in matrix of singular vectors VT
  ptr_file = fopen(Vmat,"rb");
  fread(&i, sizeof(int), 1, ptr_file);
  fread(bbfmm->VT, sizeof(real8), Vsize, ptr_file);
  fclose(ptr_file);	



  bbfmm->cutofff = cutoff[0];
  bbfmm->cutoffs = cutoff[1];

  if (domainID == 0) printf("Initialization Chebyshev done\n");
  if (domainID == 0) printf("Cutoff Info : epsilon = %g cutoff = %d dofn3 = %d\n",epsilon,bbfmm->cutoffs,dofn3_s);
  if (domainID == 0) printf("\n");

}
	


/*
 * Function: EvaluateKernelCell
 * -------------------------------------------------------------------
 * Evaluates the kernel for interactions between a pair of cells.
 * M2L operator initialization
 */
void EvaluateKernelCell(Home_t *home, real8 *field[3], real8 *source[3], 
			int Nf, int Ns, real8 *kernel) 
{
  int i, j, k, l, count, count_kernel;
  int dofNf = 6 * Nf;
  int dof2  = 54;
  real8 * Kij;
  
  int LDA_kernel = dofNf;
  
  Kij = (real8*) malloc(dof2 * sizeof(real8));

  if (Kij == NULL)
    {
      Fatal("Kij : Out of memory! %d\n",dof2);
    } 
  
  for (j=0;j<Ns;j++) {
    for (i=0;i<Nf;i++) {
      EvaluateKernel(home,  field[0][i], field[1][i],  field[2][i],
		           source[0][j],source[1][j], source[2][j], Kij);
      
      count_kernel = 6 * i + LDA_kernel * 9 * j;
      count = 0;			
      for (k=0;k<9;k++)
	for (l=0;l<6;l++, count++)
	  {
	    /* Column-major storage */
	    kernel[count_kernel + k * LDA_kernel + l] = Kij[count];
	  }
    }
  }
  
  free(Kij);
}


void bbFMMFinal(bbFMM_t *bbfmm)
{
  free(bbfmm);
}
#endif /* ifdef BBFMM */
