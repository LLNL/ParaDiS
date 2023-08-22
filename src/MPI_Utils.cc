
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "mpi_portability.h"

#include "MPI_Utils.h"

#if defined(PARALLEL)
int MPI_Rank (MPI_Comm comm) { int rank=0;                               MPI_Comm_rank(comm, &rank); return(rank); } 
int MPI_Size (MPI_Comm comm) { int size=0;                               MPI_Comm_size(comm, &size); return(size); } 
int MPI_Rank (void)          { int rank=0; MPI_Comm comm=MPI_COMM_WORLD; MPI_Comm_rank(comm, &rank); return(rank); } 
int MPI_Size (void)          { int size=0; MPI_Comm comm=MPI_COMM_WORLD; MPI_Comm_size(comm, &size); return(size); } 
#else
int MPI_Rank (void) { return(0); }
int MPI_Size (void) { return(1); }
#endif

// MPI_Hostname()
//
// This is really just an alias to the standard gethostname() call with an error check.
// I've kept it here because the MPI_Node_Rank() utility uses it to generate a color
// key for the node rank utility.
//-----------------------------------------------------------------------------------------

void MPI_Hostname (char *hostname, const int bytes)
{
   int  err = ( hostname ? gethostname(hostname,bytes) : -1 );
 
   if (err==-1) { printf("%s:%s() - error obtaining MPI node host name\n", __FILE__, __func__ ); exit(0); }
}

#if defined(PARALLEL)

// CRC_Adler32()
//
// Returns a 32-bit Adler CRC code by traversing the values in the source buffer.
// The Adler CRC is used to create a unique integer index based on the machine
// hostname. That CRC is used to create a unique comm-split between all the MPI 
// processes occupying the same cluster node.
//-----------------------------------------------------------------------------------------

static unsigned int CRC_Adler32(const void *pbuf, const size_t plen)
{
   const unsigned char *p  = (const unsigned char *) pbuf;

   unsigned int s1 = 1;
   unsigned int s2 = 0;
  
   if (p && (plen>0))
   { 
      for (size_t i=0; (i<plen); i++) 
      {
         s1 = (s1 + p[i]) % 65521;
         s2 = (s2 + s1  ) % 65521;
      }
   }

   return( (s2<<16)|s1 );
}

// MPI_Node_Rank()
//
// Will return the local MPI rank of the processes running on a particular node.
// Uses the node's host name and a local comm block to determine the local node rank.
//
// Background - many of the new heterogeneous compute architectures include devices
// like GPUs and many-in-core (MIC) type devices (e.g. Intel Knight's Corner, Phi, etc.).
// Generally, an individual compute node will contain more CPU-cores than there are 
// available GPUs or MIC devices.  Ultimately, MPI processes that are deployed across
// the cluster will need to share/coordinate access to the GPUs/MICs that exist on each node.
//
// The MPI_Node_Rank utility was developed to determine how many processes are currently
// running on an individual cluster node. It also sets a unique index for each MPI process
// running on that node. These indices can be used to load balance access to the 
// heterogeneous compute devices.
//-----------------------------------------------------------------------------------------

void MPI_Node_Rank (int & node_rank, int & node_size)
{
   node_rank=0;
   node_size=1;

   char hostname[256]; memset(hostname,0,sizeof(hostname));

   MPI_Hostname(hostname,sizeof(hostname));                            // get the hostname string (which will be used to generate the color key for the comm-split)
 
   unsigned int crc = CRC_Adler32(hostname,strlen(hostname));          // crc = 32-bit Adler CRC of hostname (used as MPI comm-split index below)

   int      mpi_rank = 0;
   MPI_Comm comm     = MPI_COMM_NULL;

   MPI_Comm_rank (MPI_COMM_WORLD,     &mpi_rank       );               // mpi_rank  = global rank index 
   MPI_Comm_split(MPI_COMM_WORLD, crc, mpi_rank, &comm);               // comm      = local communication split
   MPI_Comm_rank (comm, &node_rank);                                   // node_rank = local rank index on this node
   MPI_Comm_size (comm, &node_size);                                   // node_size = number of mpi processes executing on this node

   int   bytes = sizeof(hostname);                                     // bytes = size of each hostname buffer
   char *send  = (char *) calloc(1,          bytes*sizeof(char));      // send  = will contain this hostname
   char *recv  = (char *) calloc(1,node_size*bytes*sizeof(char));      // recv  = will contain all the hostnames included in this comm split

   strncpy(send,hostname,sizeof(hostname));                            // copy hostname to the send buffer.

   MPI_Allgather(send, bytes, MPI_CHAR, recv, bytes, MPI_CHAR, comm);  // send this hostname to the local nodes, gather the names of all the nodes in this comm split

   // check for CRC collisions...

   int local_rank=0;
   for (int i=0; (i<node_size); i++)
   {
      if (strcmp(send,recv+(i*bytes))==0)
      {
         if (i<node_rank) { ++local_rank; }
         else             { break; }
      }
      else { /* hash collision */ }
   }
  
   // if a collision occurred - send a message to the console and adjust...

   if (node_rank!=local_rank) 
   {
      printf("%s:%s() - warning : collisions occured during node rank determinaton: "
             "mpi_rank=%5d node_rank=%5d local_rank=%5d, host=%s\n", __FILE__, __func__, mpi_rank, node_rank, local_rank, hostname);
      node_rank = local_rank;
   }

   MPI_Comm_free(&comm);
 
   if (send) { free(send); send=0; }
   if (recv) { free(recv); recv=0; }
}

int MPI_Node_Rank (void) { int nr=0,ns=0; MPI_Node_Rank(nr,ns); return(nr); }
int MPI_Node_Size (void) { int nr=0,ns=0; MPI_Node_Rank(nr,ns); return(ns); }

#endif

//-----------------------------------------------------------------------------------------
// If MPI is not currently active, MPI_Node_Rank() will return proper rank/size 
// for single process execution.
//-----------------------------------------------------------------------------------------

#if !defined(PARALLEL)

void MPI_Node_Rank (int & node_rank, int & node_size) { node_rank=0; node_size=1; } 
int  MPI_Node_Rank (void) { return(0); }
int  MPI_Node_Size (void) { return(1); }

#endif


// MPI_Send_Write_Token()
// MPI_Recv_Write_Token()
//
// There are numerous cases in ParaDiS where the production of output needs to be
// serialized across the various ranks. The send/receive write token basically
// enables MPI ranks to pass a baton from process to process (thus serializing)
// the output to a single file.
//
// For example, if the simulation used 4 MPI processes and wanted to save
// all the native nodes across all processes to a restart file, the 
// sequence occurs as follows...
//    rank 0 creates the file, writes  nodes, sends token to rank 1
//    rank 1 opens the file  , appends nodes, sends token to rank 2
//    rank 2 opens the file  , appends nodes, sends token to rank 3 
//    rank 3 opens the file  , appends nodes
//    note that ranks 1-3 are blocked, waiting to receive the token.
//
// Note that this code works for both parallel and non-parallel builds.
//-----------------------------------------------------------------------------------------

#ifdef PARALLEL
void MPI_Send_Write_Token(const int rank)
{
   int token=0;
   MPI_Send(&token, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
}

void MPI_Recv_Write_Token(const int rank)
{
   int        token = 0;
   MPI_Status mpi_status;

   MPI_Recv(&token, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, &mpi_status);
}
#else
void MPI_Send_Write_Token(const int rank) {}   // stub for non-parallel builds
void MPI_Recv_Write_Token(const int rank) {}   // stub for non-parallel builds
#endif

