#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"
#include "Home.h"
#include "Domain_t.h"
#include "MPI_Cell_t.h"
#include "MPI_Comms_t.h"

//------------------------------------------------------------------------------------------------------------

#define VDELETE(a)     if (a) { delete [] a; a=0; }
#define MDELETE(a,n)   if (a) { for (int i=0; (i<n); i++) { VDELETE((a)[i]); }  delete [] a; a=0; }

// MPI_Comms_t class statics...
//------------------------------------------------------------------------------------------------------------

         real8       MPI_Comms_t::sim_lx           = 0.0;   ///< simulation box size (x)
         real8       MPI_Comms_t::sim_ly           = 0.0;   ///< simulation box size (y)
         real8       MPI_Comms_t::sim_lz           = 0.0;   ///< simulation box size (z)

         real8       MPI_Comms_t::sim_sx           = 0.0;   ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
         real8       MPI_Comms_t::sim_sy           = 0.0;   ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
         real8       MPI_Comms_t::sim_sz           = 0.0;   ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)

         int         MPI_Comms_t::sim_pbx          = 0;     ///< periodic boundaries active (x) (1=yes,0=no)
         int         MPI_Comms_t::sim_pby          = 0;     ///< periodic boundaries active (y) (1=yes,0=no)
         int         MPI_Comms_t::sim_pbz          = 0;     ///< periodic boundaries active (z) (1=yes,0=no)

         int         MPI_Comms_t::doms_nx          = 0;     ///< number of domains (x)
         int         MPI_Comms_t::doms_ny          = 0;     ///< number of domains (y)
         int         MPI_Comms_t::doms_nz          = 0;     ///< number of domains (z)

         Domain_t   *MPI_Comms_t::domains          = 0;     ///< array of domain structures
         int         MPI_Comms_t::domain_cnt       = 0;     ///< number of domains

         int         MPI_Comms_t::cells_nx         = 0;     ///< number of MPI comm cells (x)
         int         MPI_Comms_t::cells_ny         = 0;     ///< number of MPI comm cells (y)
         int         MPI_Comms_t::cells_nz         = 0;     ///< number of MPI comm cells (z)
         MPI_Cell_t *MPI_Comms_t::mpi_cells        = 0;     ///< array  of MPI comm cells
         int         MPI_Comms_t::mpi_cell_cnt     = 0;     ///< number of MPI comm cells

         int         MPI_Comms_t::mpi_rank         = 0;     ///< MPI rank of this process

         int        *MPI_Comms_t::mpi_native_cells = 0;     ///< array of cell indices that overlap my domain
         int         MPI_Comms_t::mpi_native_ccnt  = 0;     ///< number of native cells that overlap my domain

         int         MPI_Comms_t::mpi_send_cnt     = 0;     ///< number of domains we are sending data
         int        *MPI_Comms_t::mpi_send_domains = 0;     ///< array of domain indices we are sending data
         int         MPI_Comms_t::mpi_recv_cnt     = 0;     ///< number of domains we are receiving data
         int        *MPI_Comms_t::mpi_recv_domains = 0;     ///< array of domain indices we are receiving data

         int        *MPI_Comms_t::mpi_send_cnts    = 0;     ///< array of domain send buffer sizes
unsigned char      **MPI_Comms_t::mpi_send_bufs    = 0;     ///< array of domain send buffers

         int        *MPI_Comms_t::mpi_recv_cnts    = 0;     ///< array of domain receive buffer sizes
unsigned char      **MPI_Comms_t::mpi_recv_bufs    = 0;     ///< array of domain receive buffers

// constructor...
//------------------------------------------------------------------------------------------------------------

MPI_Comms_t:: MPI_Comms_t(void) {}
MPI_Comms_t::~MPI_Comms_t()     {}

void MPI_Comms_t::Recycle (void)
{
   VDELETE(domains         );
   VDELETE(mpi_cells       );
   VDELETE(mpi_native_cells);
   VDELETE(mpi_send_domains);
   VDELETE(mpi_recv_domains);
   VDELETE(mpi_send_cnts   );
   VDELETE(mpi_recv_cnts   );

   int n = (doms_nx*doms_ny*doms_nz);

   MDELETE(mpi_send_bufs,n);
   MDELETE(mpi_recv_bufs,n);

   sim_lx          = 0.0;
   sim_ly          = 0.0;
   sim_lz          = 0.0;

   sim_sx          = 0.0;
   sim_sy          = 0.0;
   sim_sz          = 0.0;

   sim_pbx         = 0;
   sim_pby         = 0;
   sim_pbz         = 0;

   doms_nx         = 0;
   doms_ny         = 0;
   doms_nz         = 0;

   domain_cnt      = 0;

   cells_nx        = 0;
   cells_ny        = 0;
   cells_nz        = 0;
   mpi_cell_cnt    = 0;

   mpi_native_ccnt = 0;

   mpi_send_cnt    = 0;
   mpi_recv_cnt    = 0;
}

//------------------------------------------------------------------------------------------------------------

void MPI_Comms_t::Initialize
(
   Home_t *home
)
{
   Recycle();

   Param_t *param   = ( home ? home->param : 0 );          // param = points to home param structure

   doms_nx          = ( param ? param->nXdoms : 0 );
   doms_ny          = ( param ? param->nYdoms : 0 );
   doms_nz          = ( param ? param->nZdoms : 0 );
   domain_cnt       = (doms_nx * doms_ny * doms_nz);

   sim_pbx          = ( param->xBoundType==Periodic ? 1 : 0 );
   sim_pby          = ( param->yBoundType==Periodic ? 1 : 0 );
   sim_pbz          = ( param->zBoundType==Periodic ? 1 : 0 );

   mpi_recv_domains = ( (domain_cnt>0) ? new int[domain_cnt] : 0 );
   mpi_send_domains = ( (domain_cnt>0) ? new int[domain_cnt] : 0 );

   domains          = ( (domain_cnt>0) ? new Domain_t[domain_cnt] : 0 );

   sim_lx           = param->Lx;
   sim_ly           = param->Ly;
   sim_lz           = param->Lz;

   sim_sx           = ( (sim_pbx && (sim_lx>0.0)) ? (1.0/sim_lx) : 0.0 );
   sim_sy           = ( (sim_pby && (sim_ly>0.0)) ? (1.0/sim_ly) : 0.0 );
   sim_sz           = ( (sim_pbz && (sim_lz>0.0)) ? (1.0/sim_lz) : 0.0 );

   int nx = cells_nx = param->nXcells;
   int ny = cells_ny = param->nYcells;
   int nz = cells_nz = param->nZcells;

   mpi_cells         = MPI_Cells_Create(nx,ny,nz, sim_lx,sim_ly,sim_lz, sim_pbx,sim_pby,sim_pbz);
   mpi_cell_cnt      = (nx*ny*nz);
}

// Set_Simulation_Size()
//
// Sets the simulation extents and periodic boundary parameters used to compute the domain and
// cell adjacency rules.
//------------------------------------------------------------------------------------------------------------

void MPI_Comms_t::Set_Simulation_Size
(
   const real8 lx ,   ///< simulation size (x)
   const real8 ly ,   ///< simulation size (y)
   const real8 lz ,   ///< simulation size (z)
   const int   pbx,   ///< periodic boundaries active (x, 1=yes, 0=no)
   const int   pby,   ///< periodic boundaries active (y, 1=yes, 0=no)
   const int   pbz    ///< periodic boundaries active (z, 1=yes, 0=no)
)
{
   sim_pbx = pbx;
   sim_pby = pby;
   sim_pbz = pbz;

   sim_lx  = lx;
   sim_ly  = ly;
   sim_lz  = lz;

   sim_sx  = ( (pbx && (lx>0.0)) ? (1.0/lx) : 0.0 );
   sim_sy  = ( (pby && (ly>0.0)) ? (1.0/ly) : 0.0 );
   sim_sz  = ( (pbz && (lz>0.0)) ? (1.0/lz) : 0.0 );
}

// Set_Rank()
//
// Sets the MPI rank for this process.
//------------------------------------------------------------------------------------------------------------

void MPI_Comms_t::Set_Rank (const int rank)
{
   mpi_rank = rank;
}

// Set_Domains()
//
// Allocates and initializes the MPI comms domain-related parameters.
//------------------------------------------------------------------------------------------------------------

void MPI_Comms_t::Set_Domains
(
   const Domain_t *doms,   ///< array initialized domains
   const int       nx  ,   ///< number of domains (x)
   const int       ny  ,   ///< number of domains (y)
   const int       nz      ///< number of domains (z)
)
{
   VDELETE(domains);

   doms_nx    = nx;
   doms_ny    = ny;
   doms_nz    = nz;
   domain_cnt = (nx*ny*nz);

   domains    = ( (doms && (domain_cnt>0)) ? new Domain_t[domain_cnt] : 0 );
   domain_cnt = ( domains ? domain_cnt : 0 );

   if (doms && domains)
   { for (int i=0; (i<domain_cnt); i++) { domains[i]=doms[i]; } }
}

// Set_Cells()
//
// Allocates and initializes the communication cells.
// Note - the simulation size must be established prior to calling this method!
//------------------------------------------------------------------------------------------------------------

void MPI_Comms_t::Set_Cells
(
   const int nx,  ///< number of communication cells (x)
   const int ny,  ///< number of communication cells (y)
   const int nz   ///< number of communication cells (z)
)
{
   VDELETE(mpi_cells);

   cells_nx = nx;
   cells_ny = ny;
   cells_nz = nz;

   mpi_cells     = MPI_Cells_Create(nx,ny,nz, sim_lx,sim_ly,sim_lz, sim_pbx,sim_pby,sim_pbz);
   mpi_cell_cnt  = ( mpi_cells ? (nx*ny*nz) : 0 );
}

//
//------------------------------------------------------------------------------------------------------------

void MPI_Comms_t::Bin_Nodes
(
         Node_t *nodes,   ///< array  of nodes
   const int     ncnt     ///< number of nodes
)
{
   MPI_Cells_Bin(mpi_cells, nodes,ncnt, cells_nx,cells_ny,cells_nz, sim_lx,sim_ly,sim_lz );
}

// Pack_Node()
//
// Will serialize all the elements of a node into an output buffer suitable for MPI exchange.
// Note - all the integer values are packed as real8's to avoid sending mixed data in an MPI
// communication call.
//------------------------------------------------------------------------------------------------------------

static real8 *Pack_Node
(
   real8   *p,   ///< buffer to received packed node
   Node_t & n    ///< source node to be packed
)
{
   if (p)
   {
      *p++ = (real8) n.myTag.domainID;
      *p++ = (real8) n.myTag.index   ;

      *p++ = (real8) n.flags;

      *p++ = n.x     ; *p++ = n.y     ; *p++ = n.z     ;
      *p++ = n.fX    ; *p++ = n.fY    ; *p++ = n.fZ    ;
      *p++ = n.vX    ; *p++ = n.vY    ; *p++ = n.vZ    ;

      *p++ = (real8) n.numNbrs;

      int narms = n.numNbrs;

      if (n.nbrTag  ) for (int i=0; (i<narms); ++i) { *p++ = (real8) n.nbrTag[i].domainID;
                                                      *p++ = (real8) n.nbrTag[i].index   ; }

      if (n.burgX   ) for (int i=0; (i<narms); ++i) { *p++ =         n.burgX [i]; }
      if (n.burgY   ) for (int i=0; (i<narms); ++i) { *p++ =         n.burgY [i]; }
      if (n.burgZ   ) for (int i=0; (i<narms); ++i) { *p++ =         n.burgZ [i]; }

      if (n.nx      ) for (int i=0; (i<narms); ++i) { *p++ =         n.nx    [i]; }
      if (n.ny      ) for (int i=0; (i<narms); ++i) { *p++ =         n.ny    [i]; }
      if (n.nz      ) for (int i=0; (i<narms); ++i) { *p++ =         n.nz    [i]; }

      if (n.armfx   ) for (int i=0; (i<narms); ++i) { *p++ =         n.armfx [i]; }
      if (n.armfy   ) for (int i=0; (i<narms); ++i) { *p++ =         n.armfy [i]; }
      if (n.armfz   ) for (int i=0; (i<narms); ++i) { *p++ =         n.armfz [i]; }
   }

   return(p);
}

//------------------------------------------------------------------------------------------------------------

static real8 *Unpack_Node
(
   Node_t & n,   ///< node to receive unpack result
   real8   *p    ///< buffer containing packed node
)
{
   if (p)
   {
      n.myTag.domainID = (int) *p++;
      n.myTag.index    = (int) *p++;

      n.flags          = (int) *p++;

      n.x  = *p++; n.y  = *p++; n.z  = *p++;
      n.fX = *p++; n.fY = *p++; n.fZ = *p++;
      n.vX = *p++; n.vY = *p++; n.vZ = *p++;

      int narms = (int) *p++;

      n.Allocate_Arms(narms);

      if (n.nbrTag) for (int i=0; (i<narms); ++i) { n.nbrTag[i].domainID = (int) *p++;
                                                    n.nbrTag[i].index    = (int) *p++; }

      if (n.burgX ) for (int i=0; (i<narms); ++i) { n.burgX [i] = *p++; }
      if (n.burgY ) for (int i=0; (i<narms); ++i) { n.burgY [i] = *p++; }
      if (n.burgZ ) for (int i=0; (i<narms); ++i) { n.burgZ [i] = *p++; }

      if (n.nx    ) for (int i=0; (i<narms); ++i) { n.nx    [i] = *p++; }
      if (n.ny    ) for (int i=0; (i<narms); ++i) { n.ny    [i] = *p++; }
      if (n.nz    ) for (int i=0; (i<narms); ++i) { n.nz    [i] = *p++; }

      if (n.armfx ) for (int i=0; (i<narms); ++i) { n.armfx [i] = *p++; }
      if (n.armfy ) for (int i=0; (i<narms); ++i) { n.armfy [i] = *p++; }
      if (n.armfz ) for (int i=0; (i<narms); ++i) { n.armfz [i] = *p++; }
   }

   return(p);
}

// Pack_Velocity()
//------------------------------------------------------------------------------------------------------------

static real8 *Pack_Velocity
(
   real8   *p,   ///< buffer to received packed node velocity
   Node_t & n    ///< source node to be packed
)
{
   if (p)
   {
      *p++ = n.myTag.domainID;
      *p++ = n.myTag.index   ;
      *p++ = n.vX;
      *p++ = n.vY;
      *p++ = n.vZ;
   }

   return(p);
}

// Unpack_Velocity()
//------------------------------------------------------------------------------------------------------------

static real8 *Unpack_Velocity
(
   Tag_t & t,   ///< unpacked velocity
   real8  *v,   ///< unpacked velocity
   real8  *p    ///< buffer to received packed node velocity
)
{
   if (p && v)
   {
      t.domainID = *p++;
      t.index    = *p++;
      v[0]       = *p++;
      v[1]       = *p++;
      v[2]       = *p++;
   }

   return(p);
}

// Pack_Nodes()
//
// Will serialize an array of nodes into an output buffer suitable for MPI exchange.
//------------------------------------------------------------------------------------------------------------

static real8 *Pack_Nodes
(
         int    &pcnt ,   ///< total count of packed values (returned)
         Node_t *nodes,   ///< array  of native nodes
   const int     ncnt     ///< number of native nodes
)
{
   int    n0    = 0;  // n0 = number of items packed for each node
          n0   += 2;  // tag
          n0   += 1;  // flags
          n0   += 3;  // position <xyz>
          n0   += 3;  // force    <xyz>
          n0   += 3;  // velocity <xyz>
          n0   += 1;  // narms

   int    n1    = 0;  // n1 = number of items packed for each arm
          n1   += 2;  // neighbor tag
          n1   += 3;  // burgers vector <xyz>
          n1   += 3;  // normal  vector <xyz>
          n1   += 3;  // arm force      <xyz>

   // count up the arms...

   int    narms = 0;  // narms = total number of arms across ALL nodes

   if (nodes && (ncnt>0))
   {
      for (int i=0; (i<ncnt); i++) { narms += nodes[i].numNbrs; }
   }

   pcnt  = 0;
   pcnt += 1;          // reserve space for node count
   pcnt += 1;          // reserve space for packed values count
   pcnt += (n0*ncnt ); // reserve space for packed nodes
   pcnt += (n1*narms); // reserve space for packed arms

   real8 *pbuf = ( (pcnt>0) ? new real8[pcnt] : 0 );

   if (nodes && pbuf)
   {
      real8 *p =pbuf;
      real8 *pe=pbuf+pcnt;

      *p++ = ncnt;
      *p++ = pcnt;
      for (int i=0; (i<ncnt) && (p<pe); i++) { p = Pack_Node(p,nodes[i]); }
   }

   return(pbuf);
}

//------------------------------------------------------------------------------------------------------------

static Node_t *Unpack_Nodes
(
   int   & ncnt ,   ///< number of packed nodes
   real8  *pbuf     ///< source packed buffer
)
{
   ncnt=0;

   Node_t *nodes=0;

   if (pbuf)
   {
      real8 *p = pbuf;

          ncnt = (int) *p++; // pick up the node count
      int pcnt = (int) *p++; // pick up the packed values count

      real8 *pe = pbuf+pcnt;

      nodes = ( (ncnt>0) ? new Node_t[ncnt] : 0 );

      if (nodes)
      { for (int i=0; (i<ncnt) && (p<pe); i++) { p = Unpack_Node(nodes[i],p); } }
   }

   return(nodes);
}

// Send_Ghosts()
//
// Will send a buffer of packed nodes to all the adjacent domains.
//------------------------------------------------------------------------------------------------------------

void MPI_Comms_t::Send_Ghosts
(
         Node_t *nodes,   ///< array  of native nodes
   const int     ncnt     ///< number of native nodes
)
{
   VDELETE(mpi_send_cnts);

   const int n=doms_nx*doms_ny*doms_nz;

   mpi_send_cnts = ( (n>0) ? new int[n] : 0 );

   if (nodes && mpi_send_cnts)
   {
      memset(mpi_send_cnts,0,n*sizeof(int)); // (reset send counts to zero)

      //int *scnts = mpi_send_cnts+(mpi_rank);

      // Serialize the native nodes into a send buffer...

      int pcnt=0;
      for (int i=0; (i<ncnt); i++) { pcnt += nodes[i].MPI_Packed_Count(); }

      real8 *pbuf = ( (pcnt>0) ? new real8[pcnt] : 0 );

      if (pbuf)
      {
         real8 *p = pbuf;
         for (int i=0; (i<ncnt); i++) { p = nodes[i].MPI_Pack(p); }
      }
   }

}

// Send_Velocities()
//------------------------------------------------------------------------------------------------------------

void MPI_Comms_t::Send_Velocities
(
         Node_t *nodes,   ///< array  of native nodes
   const int     ncnt     ///< number of native nodes
)
{
   // (tbd)
}

// Cell_Grid()
//
// Given a domain, returns a populated grid of pixels
//------------------------------------------------------------------------------------------------------------

static unsigned char *Cell_Grid
(
   Domain_t & dom
)
{
            int   nx   = MPI_Comms_t::cells_nx;
            int   ny   = MPI_Comms_t::cells_ny;
            int   nz   = MPI_Comms_t::cells_nz;
            int   n    = nx*ny*nz;

   unsigned char *grid = ( (n>0) ? new unsigned char[n] : 0 );

   if (grid)
   {
      memset(grid,0,n*sizeof(unsigned char));
   }

   return(grid);
}

// Cell bitmap support.
//
// Identifying which domains need data exchange is done by determining which domains have overlapping cells.
// We construct bitmaps of the native cells for this MPI domain and intersect them with bitmaps of the 
// remote cells for each of the remote domains.
// 
// Bitmap_Allocate() - allocates space to hold a bitmap
//------------------------------------------------------------------------------------------------------------

static unsigned int *Bitmap_Allocate 
(
   const int     ncx  ,  ///< number of cells (x)
   const int     ncy  ,  ///< number of cells (y)
   const int     ncz     ///< number of cells (z)
)
{
   const int n = (ncx*ncy*ncz);                // n = number of cells
   const int M = 8*sizeof(unsigned int);       // M = number of bits in an unsigned int (32 or 64 bits)
   const int m = (n+M-1)/M;                    // m = total number of unsigned int's to hold the bitmap

   return( (m>0) ? new unsigned int[m] : 0 );  // allocate and return the space to hold the bitmap
}
   
// Bitmap_Set() - sets and returns a bitmap of cells indexed by an array of cell indices.
//------------------------------------------------------------------------------------------------------------

static void Bitmap_Set 
(
   unsigned int *bv   ,  ///< bit vector to set
   const int    *cells,  ///< array of cell indices
   const int     ccnt ,  ///< number of cell indices
   const int     ncx  ,  ///< number of cells (x)
   const int     ncy  ,  ///< number of cells (y)
   const int     ncz     ///< number of cells (z)
)
{
   const int M = 8*sizeof(unsigned int);
   const int m = (ncx*ncy*ncz+M-1)/M;

   if (bv) { memset(bv,0,m*sizeof(unsigned int)); }
   
   if (bv && cells)
   {
      for (int i=0; (i<ccnt); i++)
      {
         unsigned int k = (unsigned int) cells[i];
         bv[k/M] |= ((unsigned int) 1) << (M-(k%M)-1);
      }
   }
}

static int Bitmap_Intersect 
(
   unsigned int *bv0  ,  ///< bit vector 
   unsigned int *bv1  ,  ///< bit vector 
   const int     ncx  ,  ///< number of cells (x)
   const int     ncy  ,  ///< number of cells (y)
   const int     ncz     ///< number of cells (z)
)
{
   const int M = 8*sizeof(unsigned int);
   const int m = (ncx*ncy*ncz+M-1)/M;

   if (bv0 && bv1)
   {
      for (int i=0; (i<m); i++)
         if ( bv0[i] & bv1[i] ) { return(1); }
   }

   return(0);
}

// Append()
//
// Appends an int to an array of int's. Simple implementation of dynamic arrays.
//------------------------------------------------------------------------------------------------------------

static int *Append (int *v, int & vcnt, int & vmax, int ival)
{
   if ( !v || (vcnt==vmax) )
   {
      vmax += ( (vmax==0) ? 32 : vmax );

      int *tmp = ( (vmax>0) ? new int[vmax] : 0 );

      if (tmp     ) { memset(tmp,0,vmax*sizeof(int)); }
      if (tmp && v) { memcpy(tmp,v,vcnt*sizeof(int)); }

      if (v) { delete [] v; }

      v = tmp;
   }

   if (v && (vcnt<vmax) ) { v[vcnt++]=ival; }

   return(v);
}


// Set_Send_Domains()
//
// Will identify all the domains with which this domain needs to share data with.
// Note - assumes overall initialization has already occurred.
//------------------------------------------------------------------------------------------------------------

void MPI_Comms_t::Set_Send_Domains (void)
{
            real8 lx  = sim_lx;
            real8 ly  = sim_ly;
            real8 lz  = sim_lz;

            int   ncx = cells_nx;
            int   ncy = cells_ny;
            int   ncz = cells_nz;

            int   pbx = sim_pbx;
            int   pby = sim_pby;
            int   pbz = sim_pbz;
  
   unsigned int  *bmap_native = Bitmap_Allocate(ncx,ncy,ncz);  // space to hold the native cell bitmap
   unsigned int  *bmap_remote = Bitmap_Allocate(ncx,ncy,ncz);  // space to hold the remote cell bitmap's

   // create a bitmap of the native cells within this rank...

   int   ccnt  = 0;                                                                             // ccnt  = cell count (init to zero)
   int  *cells = ( domains ? domains[mpi_rank].Native_Cells(ccnt,lx,ly,lz,ncx,ncy,ncz) : 0 );   // cells = array of native cell indices for this rank

   Bitmap_Set(bmap_native,cells,ccnt,ncx,ncy,ncz);

   VDELETE(cells);

   int  *sv    =0;  // local copy of domain send vector
   int   sv_cnt=0;  // number of dommains 
   int   sv_max=0;  // allocated size of the send vector

   if (bmap_native && bmap_remote && domains)
   {
      for (int j=0; (j<(doms_nx*doms_ny*doms_nz)); j++)
      {
         if (j!=mpi_rank)
         {
            int  ccnt  = 0;
            int *cells = domains[j].Remote_Cells(ccnt, lx,ly,lz, ncx,ncy,ncz, pbx,pby,pbz, 1,1,1 );
 
            Bitmap_Set(bmap_remote,cells,ccnt,ncx,ncy,ncz);

            // If there are common cell indices shared by the native cells in this domain
            // and the remote cells of the remote domain, then we need to send ghost nodes
            // to that domain....

            if ( Bitmap_Intersect(bmap_native,bmap_remote,ncx,ncy,ncz) )
               sv = Append(sv,sv_cnt,sv_max,j);

            VDELETE(cells);
         }
      }
   }

   // cleanup...

   VDELETE(bmap_native); 
   VDELETE(bmap_remote); 

   VDELETE(mpi_send_domains); 

   mpi_send_domains = sv;
   mpi_send_cnt     = sv_cnt;
}

// Print()
//------------------------------------------------------------------------------------------------------------

void MPI_Comms_t::Print (FILE *fd) const
{
   if (!fd) return;

   fprintf(fd,"mpi_comms_t\n");
   fprintf(fd,"{\n");

   fprintf(fd,"   sim_lx            : %-8.2lf \n" ,                 sim_lx            );
   fprintf(fd,"   sim_ly            : %-8.2lf \n" ,                 sim_ly            );
   fprintf(fd,"   sim_lz            : %-8.2lf \n" ,                 sim_lz            );

   fprintf(fd,"   sim_sx            : %12.6le \n" ,                 sim_sx            );
   fprintf(fd,"   sim_sy            : %12.6le \n" ,                 sim_sy            );
   fprintf(fd,"   sim_sz            : %12.6le \n" ,                 sim_sz            );

   fprintf(fd,"   sim_pbx           : %d \n"      ,                 sim_pbx           );
   fprintf(fd,"   sim_pby           : %d \n"      ,                 sim_pby           );
   fprintf(fd,"   sim_pbz           : %d \n"      ,                 sim_pbz           );

   fprintf(fd,"   doms_nx           : %d \n"      ,                 doms_nx           );
   fprintf(fd,"   doms_ny           : %d \n"      ,                 doms_ny           );
   fprintf(fd,"   doms_nz           : %d \n"      ,                 doms_nz           );

   fprintf(fd,"   domains           : 0x%016lx \n", (unsigned long) domains           );
   fprintf(fd,"   domain_cnt        : %d \n"      ,                 domain_cnt        );

   fprintf(fd,"   cells_nx          : %d \n"      ,                 cells_nx          );
   fprintf(fd,"   cells_ny          : %d \n"      ,                 cells_ny          );
   fprintf(fd,"   cells_nz          : %d \n"      ,                 cells_nz          );
   fprintf(fd,"   mpi_cells         : 0x%016lx \n", (unsigned long) mpi_cells         );
   fprintf(fd,"   mpi_cell_cnt      : %d \n"      ,                 mpi_cell_cnt      );

   fprintf(fd,"   mpi_rank          : %d \n"      ,                 mpi_rank          );

   fprintf(fd,"   mpi_native_cells  : 0x%016lx \n", (unsigned long) mpi_native_cells  );
   fprintf(fd,"   mpi_native_ccnt   : %d \n"      ,                 mpi_native_ccnt   );

   fprintf(fd,"   mpi_send_cnt      : %d \n"      ,                 mpi_send_cnt      );
   fprintf(fd,"   mpi_send_domains  : 0x%016lx \n", (unsigned long) mpi_send_domains  );
   fprintf(fd,"   mpi_recv_cnt      : %d \n"      ,                 mpi_recv_cnt      );
   fprintf(fd,"   mpi_recv_domains  : 0x%016lx \n", (unsigned long) mpi_recv_domains  );

   fprintf(fd,"   mpi_send_cnts     : 0x%016lx \n", (unsigned long) mpi_send_cnts     );
   fprintf(fd,"   mpi_send_bufs     : 0x%016lx \n", (unsigned long) mpi_send_bufs     );

   fprintf(fd,"   mpi_recv_cnts     : 0x%016lx \n", (unsigned long) mpi_recv_cnts     );
   fprintf(fd,"   mpi_recv_bufs     : 0x%016lx \n", (unsigned long) mpi_recv_bufs     );

   fprintf(fd,"}\n");
}
