/***************************************************************************
 *
 *      Module:       DislocationNetwork.cc
 *      Description:  Contains functions to analyze the number of binary and 
 *                    tertiary junctions in a dislocation network.
 *
 *      Includes:
 *                Remove2Nodes()
 *                SplitZeroPairs()
 *
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "DislocationNetwork.h"

void Remove2Nodes(Home_t *home)
{
   int i;
   int nbrNodes, nbr2node;
   real8 newPos[3];
   Node_t *node, *nbr;
   
#if 0 // Debug
   nbrNodes = 0;
   for (i = 0; i < home->newNodeKeyPtr; i++) 
   {
      
      if ((node = home->nodeKeys[i]) == (Node_t *)NULL) 
         continue;
      else
         nbrNodes++;  
   }

   printf("Before - Number of nodes %d\n",nbrNodes);

#endif

   nbr2node = 0;
   for (i = 0; i < home->newNodeKeyPtr; i++) 
   {   
      if ((node = home->nodeKeys[i]) == (Node_t *)NULL) 
         continue;
      else
         if (node->numNbrs == 2) nbr2node++;
   }

   while (nbr2node > 0)
   {
      int  mergeStatus;
      Node_t *mergedNode;
      for (i = 0; i < home->newNodeKeyPtr; i++) 
      {      
         if ((node = home->nodeKeys[i]) == (Node_t *)NULL) continue;
         if (node->numNbrs != 2) continue;
         
         nbr = GetNeighborNode(home, node, 0);
         
         newPos[X] = nbr->x;
         newPos[Y] = nbr->y;
         newPos[Z] = nbr->z;
         
         // Merge node and its neighbor at neighbor location
         MergeNode(home,OPCLASS_REMESH, node, nbr, newPos,
                   &mergedNode, &mergeStatus, 1);
         
         if (mergeStatus != MERGE_SUCCESS) 
         {
            Fatal("Merge has failed\n");
         }
         

      } // end home->newNodeKeyPtr loop

      nbr2node = 0;
      for (i = 0; i < home->newNodeKeyPtr; i++) 
      {   
         if ((node = home->nodeKeys[i]) == (Node_t *)NULL) 
            continue;
         else
            if (node->numNbrs == 2) nbr2node++;
      }
   } // end while

#if 1 // Debug
   nbrNodes = 0;
   nbr2node = 0;
   for (i = 0; i < home->newNodeKeyPtr; i++) {
      
      if ((node = home->nodeKeys[i]) == (Node_t *)NULL) 
      {
         continue;
      }
      else
      {
         nbrNodes++;
         if (node->numNbrs == 2) nbr2node++;
      }
   }

   printf("After  - Number of nodes %d\n",nbrNodes);
   printf("number of 2-nodes=%d\n",nbr2node);
   //exit(0);
#endif
}




void SplitZeroPairs(Home_t *home)
{
   int i, j, ii, jj, numNbrs,remii,remjj;
   int armList[2];
   real8 eps = 1e-5;
   real8 origPos[3], origVel[3];
   real8 burg[MAX_NBRS][3];
   Node_t *node, *splitNode1, *splitNode2;

   for (i = 0; i < home->newNodeKeyPtr; i++) 
   {
      
      if ((node = home->nodeKeys[i]) == (Node_t *)NULL) 
      {
         continue;
      }
      
      origPos[X] = node->x;
      origPos[Y] = node->y;
      origPos[Z] = node->z;
      
      origVel[X] = node->vX;
      origVel[Y] = node->vY;
      origVel[Z] = node->vZ;
      
      numNbrs = node->numNbrs;
      
      // Consider nodes that have 4 or more arms
      if (numNbrs <= 3) continue;
      
      //printf("Before - Node\n");
      //PrintNode(node);


      for (j = 0; j < numNbrs; j++) {
         burg[j][X] = node->burgX[j];
         burg[j][Y] = node->burgY[j];
         burg[j][Z] = node->burgZ[j];
         NormalizeVec(burg[j]);
      }

      // When the node is split, its number of arms change.
      // Have to loop over the removal of 2 arms.
      
      numNbrs = node->numNbrs;
      while (numNbrs >=4)
      {
         //printf("node->numNbrs=%d\n",node->numNbrs);

         remii = -1;
         remjj = -1;
         for (ii = 0; ii < numNbrs; ii++) 
            for (jj = ii+1; jj < numNbrs; jj++) 
            {
               if ( fabs(burg[ii][X] + burg[jj][X]) < eps && 
                    fabs(burg[ii][Y] + burg[jj][Y]) < eps && 
                    fabs(burg[ii][Z] + burg[jj][Z]) < eps )
               {
                  remii = ii;
                  remjj = jj;
                  break;
               }
            }
         
         //printf("remii=%d remjj=%d\n",remii,remjj);
         
         if (remii != -1)
         {
            int splitStatus;

            // Create a new node on top of the existing one 
            // with these two Burgers vectors
            armList[0] = remii;
            armList[1] = remjj;
            
            //printf("Doing the split\n");
            splitStatus = SplitNode(home, OPCLASS_SEPARATION,
                                    node, origPos, origPos,
                                    origVel, origVel, 2,
                                    armList, 0, &splitNode1,
                                    &splitNode2, 1);
            
            if (splitStatus != SPLIT_SUCCESS) 
            {
               Fatal("Split has failed\n");
            }
            
#if 0
            printf("\nAfter - splitNode1\n");
            PrintNode(splitNode1);
            printf("\nAfter - splitNode2\n");
            PrintNode(splitNode2);
#endif

            numNbrs -= 2;
         }
         else
         {
            // No two Burgers vectors found. Stop the while loop.
            numNbrs = 0;
         }
      }
   }  // end home->newNodeKeyPtr



#if 0 // Debug

   int nbrNodes = 0;
   for (i = 0; i < home->newNodeKeyPtr; i++) {
      
      if ((node = home->nodeKeys[i]) == (Node_t *)NULL) 
      {
         continue;
      }
      else
      {
         nbrNodes++;
         printf("Number of neighbors = %d\n",node->numNbrs);

         if (node->numNbrs == 4) PrintNode(node);
      }
   }

   printf("After - Number of nodes %d\n",nbrNodes);
#endif


}

void JunctionsLength(Home_t *home)
{
   int i,j, numNbrs;
   int nb2=0, nb3=0, nb4=0, nbplus=0;
   real8 Dist, dx, dy, dz;
   real8 TotLen3=0.0, TotLen4=0.0;
   Node_t *node, *nbr;
   char name[80], name2[80];
   FILE *fp, *fp2;


   printf("cycle=%d\n",home->cycle/home->param->savecnfreq);
   sprintf(name,"%s%d%s","JunctionsLength",home->cycle/home->param->savecnfreq,".txt");
   sprintf(name2,"%s%d%s","JunctionsNumber",home->cycle/home->param->savecnfreq,".txt");

   fp = fopen(name,"w");
   fp2 = fopen(name2,"w");

   for (i = 0; i < home->newNodeKeyPtr; i++) {
      
      if ((node = home->nodeKeys[i]) == (Node_t *)NULL) 
      {
         continue;
      }
      else
      {
         numNbrs = node->numNbrs;

         if (numNbrs == 2) 
            nb2++;
         else if (numNbrs ==3) 
            nb3++;
         else if (numNbrs == 4)
            nb4++;
         else
            nbplus++;

        
         for (j = 0; j < node->numNbrs; j++) 
         {
            nbr = GetNeighborNode(home, node, j);
            
            if (nbr == (Node_t *)NULL) {
               printf("WARNING: Neighbor not found at %s line %d\n",
                      __FILE__, __LINE__);
               continue;
            }
            
            dx = nbr->x - node->x;
            dy = nbr->y - node->y;
            dz = nbr->z - node->z;
            
            ZImage(home->param, &dx, &dy, &dz);                     
            Dist = sqrt(dx*dx + dy*dy + dz*dz);

            if (numNbrs ==3 && nbr->numNbrs == 3)
               TotLen3 += Dist;
            else if (numNbrs ==4 && nbr->numNbrs == 4)
               TotLen4 += Dist;
            
            fprintf(fp,"%d %d %f\n",node->numNbrs,nbr->numNbrs,Dist);

         }
      }
   }

   fprintf(fp2,"%d %d %d %d %f %f\n",nb2,nb3,nb4,nbplus,TotLen3,TotLen4);


   fclose(fp);
   fclose(fp2);

}
