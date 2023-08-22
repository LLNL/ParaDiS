//---------------------------------------------------------------------------------------------------
// Module:  WriteNodes
//
// This module was created for diagnostic/debugging the ParaDiS simulation.  It was 
// originally created to identify a communications bug that was manifesting in the 
// ghost velocity communications.
//
// This module will enable the user to save the state of the simulation across all the 
// processes currently running. Each process will create it's own output file containing
// a snapshot of the current domain decomposition as well as saving all the nodes (native+ghost)
// present in each domain.
//
// This module can potentially generate HUGE amounts of output data. As such - it's not intended
// to be used for large-scale runs and/or extended periods.
//---------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "Home.h"
#include "Node.h"
#include "Restart.h"
#include "WriteNodes.h"

//-------------------------------------------------------------------------------
// max_domains - limits the number of processes that will save output. 
//-------------------------------------------------------------------------------

const int max_domains = 256;
const int max_history = 500;

//-------------------------------------------------------------------------------
// WriteNode()
//
// Will save a single node to a diagnostic output file.
//-------------------------------------------------------------------------------

char *PrintTag(const int dom, const int indx)
{
   static char tag[32];

   sprintf(tag,"(%d,%d)",dom,indx);

   return(tag);
}

void WriteNode (FILE *fd, Node_t *node)
{
   if (fd && node)
   {
      fprintf(fd,"%-11s %c %04x %2d %2d",
                 PrintTag(node->myTag.domainID, node->myTag.index),
                 ( node->native ? 'n' : 'g' ), node->flags, node->constraint, node->numNbrs );

      fprintf(fd," %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf ",
                 node->x    , node->y    , node->z    , 
                 node->oldx , node->oldy , node->oldz  );

      fprintf(fd," %13.6le %13.6le %13.6le %13.6le %13.6le %13.6le ",
                 node->vX   , node->vY   , node->vZ   , 
                 node->oldvX, node->oldvY, node->oldvZ );

      fprintf(fd," %13.6le %13.6le %13.6le ",
                 node->fX   , node->fY   , node->fZ    );

      fprintf(fd,"\n");

      for (int j=0; (j<node->numNbrs); j++)
      {
         fprintf(fd,"     %-11s ", PrintTag(node->nbrTag[j].domainID, node->nbrTag[j].index ) );
         fprintf(fd," %11.8lf  %11.8lf %11.8lf", node->burgX[j], node->burgY[j], node->burgZ[j] );
         fprintf(fd," %11.8lf  %11.8lf %11.8lf", node->nx   [j], node->ny   [j], node->nz   [j] );
         fprintf(fd,"\n");
      }
   }
}

//-------------------------------------------------------------------------------
// WriteNodes()
//
// Will traverse an array of nodes (native or ghosts) saving the contents of 
// those nodes to an open file descriptor.
//-------------------------------------------------------------------------------

void WriteNodes(FILE *fd, Node_t **nodes, const int n)
{
   if (fd && nodes && (n>0))
   {
      for (int i=0; (i<n); i++)
         WriteNode(fd,nodes[i]);
   }
}

void WriteNodes (const char *path, const char *mode, Node_t **nodes, const int n)
{
   if (path && mode && nodes && (n>0))
   {
      FILE *fd = fopen(path,mode);
 
      if (fd)
      {
         for (int i=0; (i<n); i++)
            WriteNode(fd,nodes[i]);

         fclose(fd);
      }
   }
}

//-------------------------------------------------------------------------------
// WriteNodes()
//
// Will save all the nodes (native and ghosts) for this process.
// If this is the root process, will create the output diagnostic
// directory and also save the current domain decomposition.
//
// Each MPI process will save an individual file to the output directory.
//-------------------------------------------------------------------------------

void WriteNodes (Home_t *home)
{
    char  path[512];                    ///< used to compose path names

    int   cycle    = home->cycle;       ///< cycle    = current simulation cycle
    int   mpi_rank = home->myDomain;    ///< mpi_rank = rank index of this process
    int   mpi_size = home->numDomains;  ///< mpi_size = total number of MPI processes

    // first time through - create the subdirectory to save the node files...

    static int init=1;

    if (init)
    {
        if (mpi_rank==0)   // (only the root process creates the output directory)
        {
            snprintf(path, sizeof(path), "./%s"      , DIR_DEBUG); mkdir(path,S_IRWXU);
            snprintf(path, sizeof(path), "./%s/nodes", DIR_DEBUG); mkdir(path,S_IRWXU);
        }

        MPI_Barrier(MPI_COMM_WORLD);  // (wait for root process to create flight directory)
        init=0;
    }

    // create the cycle subdirectory (root only)

    if (mpi_rank==0)
    {   snprintf(path, sizeof(path), "./%s/nodes/%04d", DIR_DEBUG,cycle); mkdir(path,S_IRWXU); }

    // if we are only saving recent history, delete the older directory

    if ( (mpi_rank==0) && (max_history>0) && (cycle>=max_history) )
    {   
       for (int i=0; (i<mpi_size); i++)
       { snprintf(path, sizeof(path), "./%s/nodes/%04d/n%04d.dat", DIR_DEBUG, (cycle-max_history), i); remove(path); }

       { snprintf(path, sizeof(path), "./%s/nodes/%04d"          , DIR_DEBUG, (cycle-max_history)   ); rmdir (path); }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Save the node diagnostics to the output file...

    if (mpi_rank<max_domains)
    {
        snprintf(path, sizeof(path), "./%s/nodes/%04d/n%04d.dat", DIR_DEBUG,cycle,mpi_rank);

        WriteNodes(path,"w",home->nodeKeys     , home->newNodeKeyPtr );  // save the native nodes
        WriteNodes(path,"a",home->ghostNodeList, home->ghostNodeCount);  // append the ghost nodes
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

//-------------------------------------------------------------------------------
// ReadNode()
// 
// Will create, read, and return a node from an input stream that was generated
// using the WriteNode() function. Presumably this routine is reading a file 
// of diagnostics. If the previous node is given, this node will be appended to
// the list being read.
//-------------------------------------------------------------------------------

Node_t *ReadNode (FILE *fd, Node_t *prev)
{
   Node_t *node=0;

   if (fd)
   {
      size_t  buf_size = 256;
      char   *buf      = (char *) malloc(buf_size);
      ssize_t bytes    = getline(&buf,&buf_size,fd);

      if (bytes>0)
      {
         node = (Node_t *) calloc(1,sizeof(Node_t));

         char nkey; 
         sscanf(buf," (%d,%d) %c %x %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", 
                 &node->myTag.domainID, &node->myTag.index, &nkey, 
                 &node->flags, &node->constraint, &node->numNbrs,
                 &node->x    , &node->y    , &node->z     ,
                 &node->oldx , &node->oldy , &node->oldz  ,
                 &node->vX   , &node->vY   , &node->vZ    ,
                 &node->oldvX, &node->oldvY, &node->oldvZ ,
                 &node->fX   , &node->fY   , &node->fZ     );
   
         node->native = ( nkey=='n' ? 1 : 0 );
   
         int n = node->numNbrs;
    
         if (n>0)
         { 
            node->nbrTag        = (Tag_t *) calloc(  n,sizeof(Tag_t));   
            node->burgX         = (real8 *) calloc(  n,sizeof(real8));
            node->burgY         = (real8 *) calloc(  n,sizeof(real8));
            node->burgZ         = (real8 *) calloc(  n,sizeof(real8));
            node->nx            = (real8 *) calloc(  n,sizeof(real8));
            node->ny            = (real8 *) calloc(  n,sizeof(real8));
            node->nz            = (real8 *) calloc(  n,sizeof(real8));
            node->armfx         = (real8 *) calloc(  n,sizeof(real8));
            node->armfy         = (real8 *) calloc(  n,sizeof(real8));
            node->armfz         = (real8 *) calloc(  n,sizeof(real8));
            node->sigbLoc       = (real8 *) calloc(3*n,sizeof(real8));
            node->sigbRem       = (real8 *) calloc(3*n,sizeof(real8));
            node->armCoordIndex = (int   *) calloc(  n,sizeof(int  ));
   
            for (int j=0; (j<n); j++)
            {
               bytes = getline(&buf,&buf_size,fd);
  
               if (bytes>0)
               {
                  sscanf(buf," (%d,%d) %lf %lf %lf %lf %lf %lf", 
                         &node->nbrTag[j].domainID, &node->nbrTag[j].index,
                         &node->burgX[j], &node->burgY[j], &node->burgZ[j],
                         &node->nx   [j], &node->ny   [j], &node->nz   [j] );
               }
            }
         }
      }

      if (buf) { free(buf); buf=0; }
   }

   if (node && prev) { prev->next = node; }

   return(node);
}

//-------------------------------------------------------------------------------
// ReadNodes()
// 
// Will ingest an entire input file of node diagnostics.
// Will return a linked list of nodes.
// Note that you will need to recursively free the linked list to release
//-------------------------------------------------------------------------------

Node_t *ReadNodes (const char *path, int & cnt)
{
   cnt=0;
   Node_t *root=0, *node=0, *prev=0;

   FILE *fd = (FILE *) ( path ? fopen(path,"r") : 0 );

   if (fd)
   {
      while (!feof(fd))
      {
         node = ReadNode(fd,prev);

         if (node)
         {
            if (!root) { root=node; }

            prev = node;
            cnt++;
         } 
         else  break;
      }
  
      fclose(fd); 
   }

   return(root);
}

//-------------------------------------------------------------------------------
// LoadNodes()
//
// Will ingest an entire directory of diagnostic node files.
// Assumes all the node files in a directory end with ".dat" 
//-------------------------------------------------------------------------------

Node_t *LoadNodes (const char *path, int & cnt)
{
   cnt=0;

   int     dcnt=0;
   Node_t *root=0, *last=0, *node=0;

   DIR *dir=opendir(path);

   if (dir)
   {
      struct dirent *dirent;

      for (dirent=readdir(dir); (dirent); dirent=readdir(dir))
      {
         if (strcmp(dirent->d_name, "." )==0) continue;     // ignore "."
         if (strcmp(dirent->d_name, "..")==0) continue;     // ignore ".."

         if (dirent->d_type==DT_REG)
         {
            char *p = strrchr(dirent->d_name,'.');

            node=0; dcnt=0;
            if (p && (strcmp(p,".dat")==0))
            {
               char dpath[256]; sprintf(dpath,"%s/%s", path, dirent->d_name);

               node = ReadNodes(dpath,dcnt);
            }

            if (node)
            {
               if (!root) { root = last = node; }
               else       { last->next  = node; }

               // point last to the final node in the list (for appending the next file).

               while (last && last->next) { last=last->next; }
            }

            cnt+=dcnt;
         }
      }

      closedir(dir);
   }

   return(root);
}

//-------------------------------------------------------------------------------
// DeleteNodes()
//
// Will safely release all memory associated with a node list created 
// using these utilities.
//-------------------------------------------------------------------------------

void DeleteNodes (Node_t *nodes)
{
   if (nodes)
   {
      Node_t *node=0, *next=0;

      for (node=nodes; (node); node=next)
      {
         if (node->nbrTag       ) { free(node->nbrTag       ); }
         if (node->burgX        ) { free(node->burgX        ); }
         if (node->burgY        ) { free(node->burgY        ); }
         if (node->burgZ        ) { free(node->burgZ        ); }
         if (node->nx           ) { free(node->nx           ); }
         if (node->ny           ) { free(node->ny           ); }
         if (node->nz           ) { free(node->nz           ); }
         if (node->armfx        ) { free(node->armfx        ); }
         if (node->armfy        ) { free(node->armfy        ); }
         if (node->armfz        ) { free(node->armfz        ); }
         if (node->sigbLoc      ) { free(node->sigbLoc      ); }
         if (node->sigbRem      ) { free(node->sigbRem      ); }
         if (node->armCoordIndex) { free(node->armCoordIndex); }

         next = node->next;

         free(node);
      }
   }
}

