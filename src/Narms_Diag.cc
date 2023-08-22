#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Constants.h"
#include "MPI_Utils.h"
#include "Narms_Diag.h"

//------------------------------------------------------------------------------------------------------------

int  *Narms_Diag::narms_cnts  = 0;          //< array of accumulated nodal arm counts (across all domains)
int   Narms_Diag::narms_max   = MAX_NBRS;   //< maximum number of arms for a given node
int   Narms_Diag::narms_indx  = 0;          //< the current cycle index
int   Narms_Diag::narms_freq  = 100;        //< frequency of diagnostic saves

//------------------------------------------------------------------------------------------------------------
// Narms_Init() - initializes all the static members of this class
//------------------------------------------------------------------------------------------------------------

void Narms_Diag::Init (const Home_t *home)
{
   Param_t *param = ( home ? home->param : 0 );

   if (param && param->save_narms && !narms_cnts)
   {
      narms_indx = 0;
      narms_freq = ( (param->save_narms_freq>0) ? param->save_narms_freq : 0 );

      int n = (narms_max+2)*narms_freq;

      narms_cnts = ( (n>0) ? new int[n] : 0 ); 

      if (narms_cnts) { memset(narms_cnts,0,n*sizeof(int)); }
   }
}

// Cycle_Accum()
//
// At every cycle, will update the local arm counts. If MPI is active - the arm counts
// are then communicated to the root process via an MPI reduction.
//
// If the accumulation buffer is full, will trigger a save/append to the diagnostic output file.
//------------------------------------------------------------------------------------------------------------

void Narms_Diag::Cycle_Accum (const Home_t *home)
{
   Node_t **narr     = ( home ? home->nodeKeys      : 0 );   // narr     = sparse array of native node pointers
   int      ncnt     = ( home ? home->newNodeKeyPtr : 0 );   // ncnt     = number of native node pointers
   int      cycle    = ( home ? home->cycle         : 0 );   // cycle    = current cycle index
   int     *cnts     = ( narms_cnts );                       // cnts     = arm counts (all domains)
   int     *cnts_lcl = new int[narms_max];                   // cnts_lcl = arm counts (this domain)

   if (cnts_lcl) { memset(cnts_lcl,0,narms_max*sizeof(int)); }

   // accumulate the local arm counts for this domain...

   if (cnts_lcl && narr && (ncnt>0))
   {
      for (int i=0; (i<ncnt); i++)
      {
         if ( narr[i] )
         {
            int narms = narr[i]->numNbrs;

            if ( (0<=narms) && (narms<narms_max) ) { cnts_lcl[narms]++; }
         }
      }
   }

   // accumulate the arm counts across all the domains into the root domain...

   if (cnts && cnts_lcl)
   {
      cnts += (narms_max+2) * narms_indx;  

      *cnts = cycle; cnts++;
      *cnts =     0; cnts++;

      memset(cnts,0,narms_max*sizeof(int));

#ifdef PARALLEL
      MPI_Reduce ((void *) cnts_lcl, (void *) cnts, narms_max, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
#else
      memcpy(cnts,cnts_lcl,narms_max*sizeof(int));
#endif

      cnts--;
      for (int i=0; (i<narms_max); i++)
         cnts[0] += cnts[i+1]; 
   }

   if (cnts_lcl) { delete [] cnts_lcl; cnts_lcl=0; }

   narms_indx++;

   // If we've filled the buffer - save the accumulated arm counts to the output file...

   if (narms_indx==narms_freq) { Save(); }
}

//------------------------------------------------------------------------------------------------------------

void Narms_Diag::Save_GNU (void)
{
   static int saved=0;

   if ( !saved && ( MPI_Rank()==0 ) )
   {
      FILE *fd = fopen("debug/nodeArms.gnu","w");
   
      if (fd)
      {
         fprintf(fd,"#set term pngcairo size 1200,600\n");
         fprintf(fd,"#set title font 'helvetica,14'\n");
         fprintf(fd,"#set output 'nodeArms.png'\n");
         fprintf(fd,"\n");
         fprintf(fd,"reset\n");
         fprintf(fd,"\n");
         fprintf(fd,"set xlabel 'Simulation Cycle'\n");
         fprintf(fd,"set ylabel 'Arm Counts'\n");
         fprintf(fd,"#set yrange [0:50]\n");
         fprintf(fd,"\n");
         fprintf(fd,"set title  'ParaDiS Node Arm Distribution'\n");
         fprintf(fd,"\n");
         fprintf(fd,"plot  \\\n");
         for (int i=3; (i<narms_max); i++)
            fprintf(fd,"   'nodeArms.dat'   u 1:%2d  w l lw 1 title 'n=%-2d' %s", i+3, i, ( (i<(narms_max-1)) ? ", \\\n" : ";\n\n" ) );
         fprintf(fd,"pause -1;\n");
      }
   }

   saved=1;
}

//------------------------------------------------------------------------------------------------------------
// Save()
//
// Will open the diagnostic output file and flush the current accumulated counts to that file.
// When MPI is active, only the root process will save to the output file.
//------------------------------------------------------------------------------------------------------------

void Narms_Diag::Save (void)
{
   int mpi_rank = MPI_Rank();

   if ( (mpi_rank==0) && narms_cnts && (narms_indx>0) )
   {
      FILE *fd = fopen("debug/nodeArms.dat","a");   // (open in append mode)
   
      if (fd)
      {
         for (int i=0,k=0; (i<narms_indx); i++)
         {
            fprintf(fd," %7d", narms_cnts[k]); k++;  // save cycle
            fprintf(fd," %8d", narms_cnts[k]); k++;  // save total arm count

            for (int j=0; (j<narms_max); j++, k++)   // save arm counts histogram
               fprintf(fd, " %8d", narms_cnts[k]);   //

            fprintf(fd,"\n");
         }
          
         fclose(fd);
      }

      Save_GNU();
   }

   narms_indx=0;
}

//------------------------------------------------------------------------------------------------------------

void Narms_Diag_Save (const Home_t *home, const int flush)
{
   Narms_Diag narms;

   narms.Init(home);
   narms.Cycle_Accum(home);

   if (flush) { narms.Save(); }
}

