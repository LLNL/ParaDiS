#pragma once

#ifndef _PDS_MPI_UTILS_H
#define _PDS_MPI_UTILS_H

#if defined(PARALLEL)
extern int  MPI_Rank             (MPI_Comm comm);
extern int  MPI_Size             (MPI_Comm comm);
#endif

extern int  MPI_Rank             (void);
extern int  MPI_Size             (void);

extern void MPI_Hostname         (char *hostname, const int bytes);
extern void MPI_Node_Rank        (int & rank, int & size);
extern int  MPI_Node_Rank        (void);
extern int  MPI_Node_Size        (void);

extern void MPI_Send_Write_Token (const int rank);
extern void MPI_Recv_Write_Token (const int rank);

#endif

