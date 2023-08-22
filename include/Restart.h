#pragma once

#ifndef _PDS_RESTART_H
#define _PDS_RESTART_H

#include "Typedefs.h"
#include "Node.h"
#include "Segment.h"
#include "Domain_t.h"

//------------------------------------------------------------------------------------------------------------
// ParaDiS Restart file support....
//------------------------------------------------------------------------------------------------------------

class Paradis_Restart
{
   public:
      char       *rs_path         ;   ///< path to the source restart file

      int         rs_version      ;   ///< restart file version ID
      int         rs_nsegs        ;   ///< number of restart file segments (if saved in parallel)
      real8       rs_min_coords[3];   ///< extents min <x,y,z>
      real8       rs_max_coords[3];   ///< extents max <x,y,z>
   
      int         rs_decomp_type  ;   ///< domain decomposition type
      int         rs_decomp_geo[3];   ///< domain decomposition geometry <xyz>
   
      Domain_t   *rs_domains      ;   ///< array  of domains
      int         rs_domain_cnt   ;   ///< number of domains

      Node_t     *rs_nodes        ;   ///< array  of nodes
      int         rs_node_cnt     ;   ///< number of nodes

      int         rs_bv_cnt       ;   ///< number of burgers vectors in the burgers vector table
      real8      *rs_bv_tbl       ;   ///< burgers vector table

   public:
      Paradis_Restart(const char *path=0);
     ~Paradis_Restart();

      Node_t    *Nodes        (void) const { return(rs_nodes     ); }
      int        Node_Count   (void) const { return(rs_node_cnt  ); }

      Domain_t  *Domains      (void) const { return(rs_domains   ); }
      int        Domain_Count (void) const { return(rs_domain_cnt); }

      real8   X0   (void) const { return(rs_min_coords[0]); }
      real8   Y0   (void) const { return(rs_min_coords[1]); }
      real8   Z0   (void) const { return(rs_min_coords[2]); }

      real8   X1   (void) const { return(rs_max_coords[0]); }
      real8   Y1   (void) const { return(rs_max_coords[1]); }
      real8   Z1   (void) const { return(rs_max_coords[2]); }

      real8   Lx   (void) const { return(rs_max_coords[0]-rs_min_coords[0]); }
      real8   Ly   (void) const { return(rs_max_coords[1]-rs_min_coords[1]); }
      real8   Lz   (void) const { return(rs_max_coords[2]-rs_min_coords[2]); }

      int     Nx   (void) const { return( rs_decomp_geo[0]); }
      int     Ny   (void) const { return( rs_decomp_geo[1]); }
      int     Nz   (void) const { return( rs_decomp_geo[2]); }
      int     Ndom (void) const { return( rs_decomp_geo[0]
                                         *rs_decomp_geo[1]
                                         *rs_decomp_geo[2]); }

      void        Recycle         (void);

      void        Sort_Domains    (void);
      void        Sort_Nodes      (void);
      Node_t     *Find_Node       (const Tag_t & tag);

      void        Load            (const char *path);
      void        Save            (const char *path);

      void        Parse           (char *p, char *pe);
      char       *Parse_Domains   (char *p, char *pe);
      char       *Parse_Nodes     (char *p, char *pe);

      void        Print           (FILE *fd=stdout);
      void        Print           (const char *path);

      void        Print_GNU       (const char *path);
      void        Print_GNU_Doms  (const char *path);
      void        Print_GNU_Nodes (const char *path);
      void        Print_GNU_Cmd   (const char *path);

      int         Node_Count      (const int        dndx);
      Node_t     *Get_Nodes       (const int        dndx, int & ncnt);
      Node_t     *Get_Nodes       (const Domain_t & dom , int & ncnt);
      Node_t    **Get_Nodes       (int **ncnts);

      void        Index_Nodes     (void);

      Node_t    **Node_Array      (void);
      Node_t    **Node_Array      (int & ncnt);
      Node_t    **Node_Array      (const int   dndx, int & ncnt);

      Segment_t  *Segments        (int & nseg);
      Segment_t  *Segments        (int & nseg, const int dndx);

      void        ReIndex         (const int dndx=0);
      void        Refine          (const real8 smax);

      int         Validate        (const real8 eps=1.0e-6);

      void        Reset_Forces    (void);

      int         Max_Arms        (void);
      int        *Arm_Histogram   (int & max_arms);

      int         Burgers_Count   (void) const { return(rs_bv_cnt  ); }
      real8      *Burgers_Vectors (void);
      real8      *Burgers_Lengths (void);

      void        Extract         (const Paradis_Restart & rs, const real8 *p, const real8 *pw);
};


// Prototypes for functions involved in reading the restart files...

extern void AssignNodesToDomains (Home_t *home, InData_t *inData, int nodeCount, int ***nodeLists, int **listCounts);
extern void FreeNodeLists        (Home_t *home, int ***nodeLists, int **listCounts);
extern void MapV4ToV5Constraint  (int *constraint);
extern void ReadPreV4DataParams  (Home_t *home, FILE *fp, void **dataDecomp);
extern void ReadControlFile      (Home_t *home, char *ctrlFileName);
extern void ReadBinDataFile      (Home_t *home, InData_t *inData, char *dataFile);
extern void ReadNodeDataFile     (Home_t *home, InData_t *inData, char *dataFile);

#ifdef USE_HDF
extern int  ReadBinDataParams    (Home_t *home, hid_t fileID);
extern int  ReadHDFDataset       (hid_t fileID, char *datasetName, hid_t itemType, hsize_t numItems, void *itemList);
#endif


// Prototypes for functions involved in writing the restart files...

extern void SetLatestRestart     (char *fileName);
extern int  WriteCtrl            (Home_t *home, char *ctrlFile);
extern int  WriteData            (Home_t *home, char *dataFile    , int mode   , int writePrologue, int writeEpilogue );
extern void WriteRestart         (Home_t *home, char *baseFileName, int ioGroup, int firstInGroup , int writePrologue, int writeEpilogue);
extern int  WriteBinData         (Home_t *home, char *dataFile    , int mode   , int writePrologue, int writeEpilogue, BinFileData_t *binData);
extern void WriteBinaryRestart   (Home_t *home, char *baseFileName, int ioGroup, int firstInGroup , int writePrologue, int writeEpilogue, BinFileData_t *binData);

#endif
