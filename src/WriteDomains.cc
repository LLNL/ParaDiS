//---------------------------------------------------------------------------------------------------
// Module:  Support module for printing diagnostic information about the state of the 
//          domain decomposition of the simulation.
//---------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <sys/stat.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Constants.h"
#include "Typedefs.h"
#include "Decomp.h"
#include "RSDecomp.h"
#include "RBDecomp.h"
#include "WriteDomains.h"

//---------------------------------------------------------------------------------------------------
// WriteDomains()
//
// Given a pointer to an RS decomposition array (of dimension Nx x Ny x Nz), will traverse
// the array and save the decomposion to an open file descriptor.
//---------------------------------------------------------------------------------------------------

void WriteDomains(FILE *fd, const RSDecomp_t *rsd, const int nx, const int ny, const int nz)
{
   if ( fd && rsd && rsd->domBoundX && rsd->domBoundY && rsd->domBoundZ )
   {
      for (int i=0,dom=0; (i<nx); i++)
         for (int j=0; (j<ny); j++)
            for (int k=0; (k<nz); k++, dom++)
            {
               fprintf(fd," %4d  %11.4f %11.4f %11.4f  %11.4f %11.4f %11.4f\n",
                       dom, 
                       rsd->domBoundX[i],
                       rsd->domBoundY[i][j],
                       rsd->domBoundZ[i][j][k],
                       rsd->domBoundX[i+1],
                       rsd->domBoundY[i][j+1],
                       rsd->domBoundZ[i][j][k+1] );
            }
   }
}

//---------------------------------------------------------------------------------------------------
// WriteDomains()
//
// Given a pointer to an RB decomposition array (of dimension Nx x Ny x Nz), will traverse
// the array and save the decomposion to an open file descriptor. 
//---------------------------------------------------------------------------------------------------

void WriteDomains(FILE *fd, const RBDecomp_t *rbd, const int nx, const int ny, const int nz)
{
   if (fd && rbd)
   {
      for (int i=0,j=0; j<(nx*ny*nz); i++) 
      {
         // Only print the data elements of the RB array that have positive domain ID's.

         if (rbd[i].domID>=0)
         {
            fprintf(fd," %4d  %11.4f %11.4f %11.4f  %11.4f %11.4f %11.4f\n",
                     rbd[i].domID, 
                     rbd[i].cMin[0], rbd[i].cMin[1], rbd[i].cMin[2],
                     rbd[i].cMax[0], rbd[i].cMax[1], rbd[i].cMax[2]);
            j++;
         }
      }
   }
}

//---------------------------------------------------------------------------------------------------
// WriteDomains()
//
// Will save the contents of the current domain decomposition to an open file descriptor.
//---------------------------------------------------------------------------------------------------

void WriteDomains (FILE *fd, Home_t *home)
{
   if (fd && home && home->decomp)
   {
      int nx = home->param->dataDecompGeometry[0];
      int ny = home->param->dataDecompGeometry[1];
      int nz = home->param->dataDecompGeometry[2];

      if (home->param->decompType==1) { WriteDomains(fd, (RSDecomp_t *) home->decomp, nx, ny, nz); }
      if (home->param->decompType==2) { WriteDomains(fd, (RBDecomp_t *) home->decomp, nx, ny, nz); }
   }
}

//---------------------------------------------------------------------------------------------------
// WriteDomains()
//
// Will save the contents of the current domain decomposition to an output file.
//---------------------------------------------------------------------------------------------------

void WriteDomains (const char *path, Home_t *home)
{
   FILE *fd = (FILE *) ( (path && home) ? fopen(path,"w") : 0 );

   if (fd)
   {
      WriteDomains(fd,home);
      fclose(fd);
   }
}

//---------------------------------------------------------------------------------------------------
// WriteDomains_GNU()
//
// Will save a gnuplot command file that can be used to visualize the domain decomposition.
// If you uncomment the terminal and output, will save the result as an animated GIF.
//---------------------------------------------------------------------------------------------------

void WriteDomains_GNU(const char *path)
{
   const char cmd[] = 
      "#set term gif animate size 640,480\n"
      "#set output 'doms.gif'\n"
      "\n"
      "set nokey\n"
      "\n"
      "list=system('ls d????.dat')\n"
      "\n"
      "do for [file in list] {\n"
      "   set title sprintf('%%s', file)\n"
      "\n"
      "   splot \\\n"
      "      file u 2:3:4:($5-$2):  (0)  :  (0)     w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  , \\\n"
      "      file u 2:3:4:  (0)  :($6-$3):  (0)     w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  , \\\n"
      "      file u 5:3:4:  (0)  :($6-$3):  (0)     w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  , \\\n"
      "      file u 2:6:4:($5-$2):  (0)  :  (0)     w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  , \\\n"
      "      file u 2:3:7:($5-$2):  (0)  :  (0)     w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  , \\\n"
      "      file u 2:3:7:  (0)  :($6-$3):  (0)     w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  , \\\n"
      "      file u 5:3:7:  (0)  :($6-$3):  (0)     w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  , \\\n"
      "      file u 2:6:7:($5-$2):  (0)  :  (0)     w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  , \\\n"
      "      file u 2:3:4:  (0)  :  (0)  :($7-$4)   w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  , \\\n"
      "      file u 5:3:4:  (0)  :  (0)  :($7-$4)   w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  , \\\n"
      "      file u 2:6:4:  (0)  :  (0)  :($7-$4)   w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  , \\\n"
      "      file u 5:6:4:  (0)  :  (0)  :($7-$4)   w vectors nohead lt 1 lw 1 lc rgb 'dark-violet'  ; \n"
      "\n"
      "   pause 0.1;\n"
      "}\n";

   FILE *fd = (FILE *) ( path ? fopen(path,"w") : 0 );

   if (fd)
   {
      fprintf(fd,cmd);
      fclose(fd);
   }
}

//---------------------------------------------------------------------------------------------------

void WriteDomains (Home_t *home)
{
    int   cycle    = home->cycle;       ///< cycle    = current simulation cycle
    int   mpi_rank = home->myDomain;    ///< mpi_rank = rank index of this process

    // first time through - create the subdirectory to save the domain files...

    static int init=1;

    if (init)
    {
        if (mpi_rank==0)   // (only the root process creates the output directory)
        {
            char path[512];
            snprintf(path, sizeof(path), "./%s"              , DIR_DEBUG); mkdir(path,S_IRWXU);
            snprintf(path, sizeof(path), "./%s/doms"         , DIR_DEBUG); mkdir(path,S_IRWXU);
            snprintf(path, sizeof(path), "./%s/doms/doms.gnu", DIR_DEBUG); WriteDomains_GNU(path);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        init=0;
    }

    // save the domain decomposition for each cycle (root only)

    if (mpi_rank==0)
    {
        char path[512];
        snprintf(path, sizeof(path), "./%s/doms/d%04d.dat", DIR_DEBUG,cycle); WriteDomains(path,home);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

