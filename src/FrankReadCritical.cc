#ifdef FRCRIT

/***************************************************************************
 *
 *  Function    : ANI_FrankRead.c
 *  Description : subroutines for computing critical stress of FR source
 *
 **************************************************************************/
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Init.h"
#include "Restart.h"
#include "Comm.h"
#include "QueueOps.h"
#include "WriteNodes.h"

static void RemoveAllNodes(Home_t *home)
{
   int    i, j, globalOp = 1;
   Node_t *node, *nbr;
   
   for (i = 0; i < home->newNodeKeyPtr; i++) 
   {
      node = home->nodeKeys[i];
      if (node == (Node_t *)NULL) continue;

      j = 0;      
      while (j < node->numNbrs) 
      {
         nbr = GetNeighborNode(home, node, j);
         if (nbr == (Node_t *)NULL) {
            j++;
            continue;
         }

         nbr->constraint = 0;
         
         ChangeArmBurg(home, node, &node->nbrTag[j], 0, 0, 0,
                       0, 0, 0, globalOp, DEL_SEG_NONE);
         ChangeArmBurg(home, nbr, &node->myTag, 0, 0, 0,
                       0, 0, 0, globalOp, DEL_SEG_NONE);

         if (nbr->numNbrs == 0)
         {
            RemoveNode(home, nbr, 0); 
            nbr == (Node_t *)NULL;
         }
      }
   } 

   for (i = 0; i < home->newNodeKeyPtr; i++) 
   {
      node = home->nodeKeys[i];
      if (node == (Node_t *)NULL) continue;
      node->constraint = 0;
         if (node->numNbrs == 0)
         {
            RemoveNode(home, node, 0);            
            node == (Node_t *)NULL;
         }
   }

   home->newNodeKeyPtr = 0;
}

static real8 MaxDisplacement(Home_t *home)
{
   Node_t *node, *nbrNode;

   Param_t *param;
   param = home->param;

   // Compute the maximum displacement of the nodes in the simulation
   real8 maxdisp = 0.0;
   
   for (int i = 0; i < home->newNodeKeyPtr; i++)
     {
       if (home->nodeKeys[i] != (Node_t *)NULL)
	 {
	   node = home->nodeKeys[i];

	   if (HAS_ALL_OF_CONSTRAINTS(node->constraint, PINNED_NODE)) continue;

	   // Get the displacement vector
	   real8 oldpos[3];
	   oldpos[X] = node->oldx;
	   oldpos[Y] = node->oldy;
	   oldpos[Z] = node->oldz;

	   real8 pos[3];
	   pos[X] = node->x;
	   pos[Y] = node->y;
	   pos[Z] = node->z;

	   real8 disp[3];
	   disp[X] = pos[X]-oldpos[X];
	   disp[Y] = pos[Y]-oldpos[Y];
	   disp[Z] = pos[Z]-oldpos[Z];

	   ZImage(param, &disp[X], &disp[Y], &disp[Z]);


	   int numNbrs = node->numNbrs;
	   for (int j = 0; j < numNbrs; j++)
	     {	    
	       nbrNode = GetNeighborNode(home, node, j);

	       if (nbrNode == (Node_t *)NULL) {
		 printf("WARNING: Neighbor not found at %s line %d\n",
			__FILE__, __LINE__);
		 continue;
	       }

	       // Get the segment line
	       real8 rt[3];
	       rt[0] = nbrNode->x - node->x;
	       rt[1] = nbrNode->y - node->y;
	       rt[2] = nbrNode->z - node->z;
	       
	       ZImage(param, &rt[0], &rt[1], &rt[2]);
	       
	       real8 pos2[3];
	       pos2[X] = node->x + rt[0];
	       pos2[Y] = node->y + rt[1];
	       pos2[Z] = node->z + rt[2];
	    
	       // Get the segment plane
	       real8 plane[3];
	       plane[X] = node->nx[j];
	       plane[Y] = node->ny[j];
	       plane[Z] = node->nz[j];

	       // Get the segment glide direction
   	       real8 g[3];
	       NormalizedCrossVector(rt,plane,g);

	       // Project the displacement along the glide direction
	       real8 dispdotg = fabs(disp[X]*g[X] + disp[Y]*g[Y] + disp[Z]*g[Z]);
	       
	       if (dispdotg > maxdisp) maxdisp = dispdotg;
	     }
	 }
     }

   return maxdisp;
}


static real8 BestFit(real8 L[10])
  {
    real8 meanX = 0, meanY = 0;

    for (int i = 1; i < 10; i++)
      {
	meanX += i;
	meanY += L[i];
      }
    meanX /= 9;
    meanY /= 9;

    real8 denom = 0.0;
    for (int i = 1; i < 10; i++)
      {
	denom += ((i*1.0-meanX)*(i*1.0-meanX));
      }

    real8 slope=0.0;
    for (int i = 1; i < 10; i++)
      {
	slope += (i-meanX)*(L[i]-meanY);
      }
    return slope/denom;
  }

static real8 DislocationLength(Home_t *home)
{
  int i,j;
  real8 length, rt[3], rtSquared;
  Node_t *node, *nbr;

  Param_t *param;
  param = home->param;

  length = 0;
  for (i = 0; i < home->newNodeKeyPtr; i++) 
  {
      node = home->nodeKeys[i];
      if (node == (Node_t *)NULL) continue;
      
      for (j = 0; j < node->numNbrs; j++) 
      {
          nbr = GetNeighborNode(home, node, j);
          if (nbr == (Node_t *)NULL) continue;
      
/*
 *        Ensure node is the node with the lower tag
 */
	  if (OrderNodes(node, nbr) != -1) {
	    continue;
	  }	  
	  
	  rt[0] = nbr->x - node->x;
          rt[1] = nbr->y - node->y;
          rt[2] = nbr->z - node->z;

          ZImage(param, &rt[0], &rt[1], &rt[2]);

          rtSquared = DotProduct(rt, rt);

          length += sqrt(rtSquared);
      }
  }
  return length;
}

static int FR_Status(Home_t *home, real8 alpha)
{
   int i, status;

   /* Check for stability */
   real8 maxdisp =  MaxDisplacement(home);
 
   status = 0 ;                        /* default: status = 0 (undecided)  */
   //if (maxdisp < 0.01)  status =-1;    /* conclude: FR source is stable    */
   //if (maxdisp > 0.5)   status = 1;    /* conclude: FR source is unstable  */


   if (maxdisp > 0.5)
     status = 1;
   else if (maxdisp<0.05)
     status =-1;
   else
     status = 0;
   

   printf("Doing alpha=%g  max displacement = %g\n",alpha,maxdisp);

   return status;
}


static int CalculateNumberNodes(Home_t *home)
{
  int numNodes = 0, i;
  Param_t *param;
  param = home->param;
  Node_t *node;
  
  for (i = 0; i < home->newNodeKeyPtr; i++) 
  {
      node = home->nodeKeys[i];
      if (node == (Node_t *)NULL) continue;

      numNodes++;
  }

  return numNodes;
}



void Calc_FR_Critical_Stress(int argc, char *argv[], Home_t **homeptr)
{
   Home_t  *home;
   Param_t *param;
   int cycleEnd;
   char ctrlFile[512], dataFile[512], rootFile[512];
   InData_t     *inData;
   int i, iter, Niter=15;
   real8 stress0[6], alpha, alpha_min=0, alpha_max=0, dalpha_eps=1e-3;
   FILE *fp;

#ifdef FPES_ON
   unmask_std_fpes();
#endif

   home  = *homeptr;
   param = home->param;

   int totinc = 0;
#ifdef ESHELBY	   
   totinc = home->totInclusionCount;
#endif	     
   
   
   printf("\nEntering Calc_FR_Critical_Stress\n\n");
   printf("Will do %d steps at each stress as set in the control file\n", param->maxstep);

   // Initial stress and dislocation length
   real8 normstress = 0;
   for (i = 0; i < 6; i++)
   {
      stress0[i] = param->appliedStress[i];
      normstress +=  param->appliedStress[i] * param->appliedStress[i];
   }
   normstress = sqrt(normstress);

   // Range for bisection
   FILE *ftmp;
   ftmp = fopen("alpharange.in","r");
   if (ftmp == NULL)
     {
       ftmp = fopen("alpharange.in","w");
       fprintf(ftmp,"%lf %lf",0.0,100.0);
       alpha_min =    0.0;
       alpha_max =  100.0;
     }
   else
     {
       fscanf(ftmp,"%lf %lf",&alpha_min,&alpha_max);
       if (alpha_min > alpha_max)
	 {
	   Fatal("alpha_min=%g > alpha_max=%g\n",alpha_min,alpha_max);
	 }
     }
   fclose(ftmp);

   printf("Starting interval for bisection is [%f, %f]\n",alpha_min,alpha_max);

   // Create restart file
   mkdir("restart",S_IRWXU);

   // Save init file
   char baseName[128];

   sprintf(baseName,"init.ctrl");
   WriteRestart(home, baseName, 0, 1, 1, 1);

   fp = fopen("Frank-Read.txt","a");
   
   // Before the bisection, make sure that the provided range is good
   alpha = alpha_max;
  
   // set stress according to current alpha
   normstress = 0;
   for (i = 0; i < 6; i++)
     {
       param->appliedStress[i] = stress0[i]*alpha;
       normstress +=  param->appliedStress[i] * param->appliedStress[i];
     }
   normstress = sqrt(normstress);

   // Don't do anything until startmcycle
   int startmcycle = 500;
   
   // Start regular ParaDiS
   home->cycle = param->cycleStart; 
   cycleEnd    = param->cycleStart + param->maxstep;
   
   // Check FR stability
   printf("Running alpha=%g with stress = %g\n",alpha,normstress);

   // Run startmcycle time step first
   while(home->cycle < startmcycle)
     {
       ParadisStep(home);
     }

   real8 L0 = DislocationLength(home);
   printf("L0=%g\n",L0);

   int status = -1;
   real8 Larray[10];
   for (int iY = 0; iY < 10; iY++) Larray[iY]=0;
   
   real8 L = L0;
   real8 slope = 0.0;
   real8 Linc = 0.05;
   while(home->cycle < cycleEnd)   
     {
       ParadisStep(home);

       L = DislocationLength(home);

       if ( fabs(L-L0)/L0 > Linc )
	 {
	   status = 1; // unstable
	   break;
	 }

       // Keep last 10 lengths
       Larray[home->cycle%10] = fabs(L-L0);
       
       if ( (home->cycle % 10) == 0)
	 {
	   slope = BestFit(Larray);
	   printf("cycle=%d L0=%10.5g L=%g (L-L0)/L0=%10.5g < %10.5g slope=%10.5g\n",
		  home->cycle,L0,L,fabs(L-L0)/L0,Linc,slope);
	 }
     }

   // If the length has consistently increased the last 10 times,
   // consider that the source is unstable
   if (status == -1 && slope > 0.6)
     {
       status = 1;
     }

   if (status == -1)
     {
       printf("alpha_max = %g is not a lower bound. Norm of stress is %g status =%d\n",
	      alpha_max,normstress,status);
       alpha_min =     alpha_max;
       alpha_max = 2 * alpha_max;
     }
   else
   {
      printf("alpha_max = %g is an upper bound.\n",alpha);
   }

   // reload config if stress changes
   RemoveAllNodes(home);
   strcpy(dataFile, home->param->dirname);
   
   //char cwd[1024];
   //getcwd(cwd, sizeof(cwd));
   //printf("Current working dir: %s\n", cwd);
   
   ParadisFinish(*homeptr);
   chdir("../..");
   
   //getcwd(cwd, sizeof(cwd));
   //printf("Current working dir: %s\n", cwd);
   
   
   strcat(dataFile, "/restart/init.ctrl");
   
   strcpy(argv[1],dataFile);
   
   
   ParadisInit(homeptr, argc, argv);
   home  = *homeptr;
   param = home->param;
 
   // Start bisection
   printf("Starting bisection. Will do %d iterations \n", Niter);
   printf("Stress range is  [%8g, %8g]\n\n",alpha_min*normstress, alpha_max*normstress);
   alpha = (alpha_max + alpha_min) * 0.5 ; 

   for(iter = 0; iter < Niter; iter ++)
   {
       printf("\nRunning ParaDiS at iter = %d/%d  alpha = %.6f\n\n", iter, Niter, alpha);
 
       // set stress according to current alpha
       for (i = 0; i < 6; i++)
           param->appliedStress[i] = stress0[i]*alpha;

       // Start regular ParaDiS
       home->cycle = param->cycleStart; 
       cycleEnd    = param->cycleStart + param->maxstep;

       // Check FR stability
       printf("Running alpha=%g with stress = %g\n",alpha,normstress);
       
       // Run startmcycle time step first
       while(home->cycle < startmcycle)
	 {
	   ParadisStep(home);
	 }

       real8 L0 = DislocationLength(home);
       printf("L0=%g\n",L0);

       status = -1;
       L = L0;

       real8 Larray[10];
       for (int iY = 0; iY < 10; iY++) Larray[iY]=0;
       
       while(home->cycle < cycleEnd)
       {
          ParadisStep(home);

	  L = DislocationLength(home);

	  if ( fabs(L-L0)/L0 > Linc)
	    {
	      status = 1;
	      break;
	    }

	  // Keep last 10 lengths
	  Larray[home->cycle%10] = fabs(L-L0);

       if ( (home->cycle % 10) == 0)
	 {
	   slope = BestFit(Larray);
	   printf("cycle=%d L0=%10.5g L=%g (L-L0)/L0=%10.5g < %10.5g slope=%10.5g\n",
		  home->cycle,L0,L,fabs(L-L0)/L0,Linc,slope);
	 }
	  
       }

       // If the length has consistently increased the last 10 times,
       // consider that the source is unstable
       if (status == -1 && slope > 0.6)
	 {
	   status = 1;
	 }

       fprintf(fp, "alpha = %15.10g  status = %2d   range=[%8g, %8g]\n", 
               alpha, status, alpha_min, alpha_max);
       printf("alpha = %15.10g  status = %2d   range=[%8g, %8g]\n", 
              alpha, status, alpha_min, alpha_max);
       fflush(fp);

       // Decide next alpha 
       if (status == -1)
       {
	 printf("FR source has reached a stable configuration : L0=%g L=%g\n",L0,L);
           alpha_min = alpha;
           alpha = (alpha_max + alpha_min) * 0.5 ;

           ftmp = fopen("alpharange.in","w");
           fprintf(ftmp,"%f %f",alpha_min,alpha_max);
           fclose(ftmp);
           
           if (  (alpha_max - alpha_min)  < dalpha_eps)
           {
              printf("desired accuracy reached (dalpha_eps = %g)\n", dalpha_eps);
              break;
           }

           printf("Niter=%d status = 0 : alpha = %f in [%f, %f]\n",
		  iter, alpha, alpha_min,alpha_max);
       }

       if (status == 1)
       {
           printf("FR source is unstable : L0=%g L=%g\n",L0,L);
           alpha_max = alpha;
           alpha = (alpha_max + alpha_min) * 0.5 ;

           ftmp = fopen("alpharange.in","w");
           fprintf(ftmp,"%f %f",alpha_min, alpha_max);
           fclose(ftmp);

           if (  (alpha_max - alpha_min)  < dalpha_eps)
           {
              printf("desired accuracy reached (dalpha_eps = %g)\n", dalpha_eps);
              break;
           }
           printf("Niter=%d status = 1 : alpha = %f in [%f, %f]\n",
		  iter, alpha, alpha_min, alpha_max);
       }

       FILE *fsemifinal;
       fsemifinal = fopen("semifinalalpha.txt","w");
       fprintf(fsemifinal,"%d %15.10g %15.10g %15.10g %d  %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",
	       totinc,alpha,alpha_min,alpha_max,iter,
	       param->appliedStress[0],param->appliedStress[1],param->appliedStress[2],
	       param->appliedStress[3],param->appliedStress[4],param->appliedStress[5]);
       fclose(fsemifinal);
       
       // reload config if stress changes
       if (iter < Niter - 1) 
       {
          RemoveAllNodes(home);
          strcpy(dataFile, home->param->dirname);           
          strcpy(rootFile, home->param->dirname);           

          ParadisFinish(*homeptr);
          chdir("../..");

          strcat(dataFile, "/restart/init.ctrl");
          strcpy(argv[1],dataFile);


          ParadisInit(homeptr, argc, argv);
          home  = *homeptr;
          param = home->param;
       }
   }  // end Niter

   fclose(fp);

   FILE *ffinal;
   ffinal = fopen("finalalpha.txt","w");
   fprintf(ffinal,"%d %15.10g %15.10g %15.10g  %15.10g %15.10g %15.10g %15.10g %15.10g %15.10g\n",
	   totinc,alpha,alpha_min,alpha_max,
	   param->appliedStress[0],param->appliedStress[1],param->appliedStress[2],
	   param->appliedStress[3],param->appliedStress[4],param->appliedStress[5]);
   fclose(ffinal);

   printf("Calc_FR_Critical_Stress finished\n");
}

#endif
