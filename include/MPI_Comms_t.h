#pragma once

#ifndef _PDS_MPI_COMMS_H
#define _PDS_MPI_COMMS_H

#include "mpi_portability.h"
#include "cuda_portability.h"
#include "Typedefs.h"
#include "Domain_t.h"
#include "Node.h"

//------------------------------------------------------------------------------------------------------------

class MPI_Comms_t
{
   public :
      static          real8       sim_lx           ;   ///< simulation box size (x)
      static          real8       sim_ly           ;   ///< simulation box size (y)
      static          real8       sim_lz           ;   ///< simulation box size (z)

      static          real8       sim_sx           ;   ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
      static          real8       sim_sy           ;   ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
      static          real8       sim_sz           ;   ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)

      static          int         sim_pbx          ;     ///< periodic boundaries active (x) (1=yes,0=no)
      static          int         sim_pby          ;     ///< periodic boundaries active (y) (1=yes,0=no)
      static          int         sim_pbz          ;     ///< periodic boundaries active (z) (1=yes,0=no)

      static          int         doms_nx          ;     ///< number of domains (x)
      static          int         doms_ny          ;     ///< number of domains (y)
      static          int         doms_nz          ;     ///< number of domains (z)

      static          Domain_t   *domains          ;     ///< array of domain structures
      static          int         domain_cnt       ;     ///< number of domains

      static          int         cells_nx         ;     ///< number of MPI comm cells (x)
      static          int         cells_ny         ;     ///< number of MPI comm cells (y)
      static          int         cells_nz         ;     ///< number of MPI comm cells (z)

      static          MPI_Cell_t *mpi_cells        ;     ///< array  of MPI comm cells
      static          int         mpi_cell_cnt     ;     ///< number of MPI comm cells

      static          int         mpi_rank         ;     ///< MPI rank of this process

      static          int        *mpi_native_cells ;     ///< array of cell indices that overlap my domain
      static          int         mpi_native_ccnt  ;     ///< number of native cells that overlap my domain

      static          int         mpi_send_cnt     ;     ///< number of domains we are sending data
      static          int        *mpi_send_domains ;     ///< array of domain indices we are sending data
      static          int         mpi_recv_cnt     ;     ///< number of domains we are receiving data
      static          int        *mpi_recv_domains ;     ///< array of domain indices we are receiving data

      static          int        *mpi_send_cnts    ;     ///< array of domain send buffer sizes
      static unsigned char      **mpi_send_bufs    ;     ///< array of domain send buffers

      static          int        *mpi_recv_cnts    ;     ///< array of domain receive buffer sizes
      static unsigned char      **mpi_recv_bufs    ;     ///< array of domain receive buffers

   public :
      MPI_Comms_t(void);
     ~MPI_Comms_t();

     void  Recycle             (void);
     void  Initialize          (Home_t *home);

     void  Set_Simulation_Size (const real8 lx , const real8 ly , const real8 lz ,
                                const int   pbx, const int   pby, const int   pbz );

     void  Set_Rank            (const int       rank);
     void  Set_Domains         (const Domain_t *doms, const int nx, const int ny, const int nz);
     void  Set_Cells           (                      const int nx, const int ny, const int nz);

     void  Bin_Nodes           (      Node_t   *nodes, const int   ncnt);

     void  Set_Send_Domains    (void);

     void  Send_Ghosts         (      Node_t   *nodes, const int   ncnt);
     void  Send_Velocities     (      Node_t   *nodes, const int   ncnt);

     void  Print               (FILE *fd=stdout) const;
};

//------------------------------------------------------------------------------------------------------------

#endif  //  _PDS_MPI_COMMS_H
