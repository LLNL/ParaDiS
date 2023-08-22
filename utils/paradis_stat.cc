#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>

#include "Args.h"
#include "File.h"
#include "Directory.h"
#include "Restart.h"
#include "V3.h"

//----------------------------------------------------------------------------------------------------------

#define VDELETE(a) if (a) { delete [] a; a=0; }

//----------------------------------------------------------------------------------------------------------

char  *src_path     = 0;
char  *dst_path     = 0;

char **src_paths    = 0;
int    src_path_cnt = 0;

char  *density_path = 0;

real8  bv_scale     = 1.0;

real8  bmag         = 2.725e-10;

#define PRINT_BV_SUMS
#ifdef  PRINT_BV_SUMS
const real8 bv_ptbl [] = 
{
    0.577350,  0.577350,  0.577350,
    0.577350,  0.577350, -0.577350,
    0.577350, -0.577350,  0.577350,
    0.577350, -0.577350, -0.577350,
    1.154701,  0.000000,  0.000000,
    0.000000,  1.154701,  0.000000,
    0.000000,  0.000000,  1.154701
};
         
const int   bv_ptbl_cnt = sizeof(bv_ptbl)/(3*sizeof(real8));

      int   bv_ptbl_cnts[bv_ptbl_cnt] = { 0 };
      real8 bv_ptbl_sums[bv_ptbl_cnt] = { 0.0 };
#endif

//----------------------------------------------------------------------------------------------------------

const char  usage[] =
   "                                                                                              \n"
   " The ParaDiS stat application will produce statistics on a single restart file or collection  \n"
   " of restart files within a directory.                                                         \n"
   "                                                                                              \n"
   " usage: paradis_stat <path> -bs <scale>                                                       \n"
   "                                                                                              \n"
   " where                                                                                        \n"
   "        <path>    identifies a single restart file or directory containing a collection       \n"
   "                  collection of restart files (required).                                     \n"
   "                                                                                              \n"
   "    -bs   <scale> burgers vector scaling factor                                               \n"
   "                                                                                              \n"
   "    -bmag <val>   magnitude of the burgers vector                                             \n"
   "                                                                                              \n"
   "    -density <path>   path to the output density file                                         \n"
   "                                                                                              \n"
   "\n";

//----------------------------------------------------------------------------------------------------------

static void Process_Args (int argc, char *argv[])
{
   int indx=0;

   if ( Get_Flag("-help"     ,argc,argv) ) { fprintf(stderr,"%s", usage ); exit(0); }
   if ( Get_Flag("--help"    ,argc,argv) ) { fprintf(stderr,"%s", usage ); exit(0); }
   if ( argc==1 )                          { fprintf(stderr,"%s", usage ); exit(0); }

   src_path = strdup(argv[1]);

   if ( (indx=Find_Arg("-bs"      ,argc,argv))>0 ) { bv_scale     = (real8) strtod(argv[indx+1],0); }
   if ( (indx=Find_Arg("-bmag"    ,argc,argv))>0 ) { bmag         = (real8) strtod(argv[indx+1],0); }
   if ( (indx=Find_Arg("-gnu"     ,argc,argv))>0 ) { dst_path     = strdup(argv[indx+1]); }
   if ( (indx=Find_Arg("-density" ,argc,argv))>0 ) { density_path = strdup(argv[indx+1]); }
}

//----------------------------------------------------------------------------------------------------------

class Paradis_Stat
{
   public:
      char  *pds_path     ;  ///< path to source restart file
      int    node_cnt     ;  ///< number of nodes

      real8  pds_lx       ;  ///< simlation size (x)
      real8  pds_ly       ;  ///< simlation size (y)
      real8  pds_lz       ;  ///< simlation size (z)

      int    seg_cnt      ;  ///< number of unique segments
      real8  seg_min      ;  ///< minimum segment length (b)
      real8  seg_max      ;  ///< maximum segment length (b)
      real8  seg_mean     ;  ///< mean    segment length (b)
      real8  seg_sdev     ;  ///< std. deviation of segment lengths
      real8  seg_sum      ;  ///< total   segment length (b)
      real8  seg_density  ;  ///< overall segment density (b)/volume

      int    nv_arms_max  ;  ///< max node arms
      int   *nv_arms_histo;  ///< node arms histogram

      int    bv_tbl_cnt   ;  ///< number of burgers vectors
      real8 *bv_tbl       ;  ///< burgers vector table (x,y,z)
      int   *bv_tbl_cnts  ;  ///< burgers vector segment counts
      real8 *bv_tbl_sums  ;  ///< burgers vector segment length sums

   public:
      Paradis_Stat (void);
     ~Paradis_Stat ();

      void Recycle (void);

      real8 BV_Sum   (int & cnt, const real8 *bv);
      void  BV_Sums  (real8 *bv_sums, int *bv_cnts, const real8 *bv, const int bv_cnt );

      void Load      (Paradis_Restart & rs);
      void Print     (FILE *fd=stdout);
};

void Paradis_Stat::Recycle (void)
{
   VDELETE(pds_path     );
   VDELETE(nv_arms_histo); nv_arms_max = 0;
   VDELETE(bv_tbl       ); bv_tbl_cnt  = 0;
   VDELETE(bv_tbl_cnts  ); 
   VDELETE(bv_tbl_sums  ); 
}

Paradis_Stat:: Paradis_Stat (void) 
{
   pds_path      = 0;
   node_cnt      = 0;

   pds_lx        = 0.0;
   pds_ly        = 0.0;
   pds_lz        = 0.0;

   seg_cnt       = 0;
   seg_min       = 0.0;
   seg_max       = 0.0;
   seg_mean      = 0.0;
   seg_sdev      = 0.0;
   seg_sum       = 0.0;
   seg_density   = 0.0;

   nv_arms_max   = 0;
   nv_arms_histo = 0;

   bv_tbl_cnt    = 0;
   bv_tbl        = 0;
   bv_tbl_cnts   = 0;
   bv_tbl_sums   = 0;
}

Paradis_Stat::~Paradis_Stat ()     { Recycle(); }

static int BV_Index 
(
   const real8 *bv  ,   ///< source burgers vector
   const real8 *btbl,   ///< burgers table
   const int    bcnt    ///< length of the burgers table
)
{
   if (bv && btbl && (bcnt>0))
   {
      for (int i=0,j=0; (i<bcnt); i++, j+=3)
      {
         const real8 eps = 1.0e-6;

         if (    (fabs(bv[0]-btbl[j  ])<eps)
              && (fabs(bv[1]-btbl[j+1])<eps)
              && (fabs(bv[2]-btbl[j+2])<eps) ) return(i); // (found)
      }
   }

   return(-1); // (not found)
}

int BV_Equiv (const real8 *p1, const real8 *p2)
{
   if (p1 && p2)
   {
      real8 t1[3]; V3_NORMALIZE(t1,p1);
      real8 t2[3]; V3_NORMALIZE(t2,p2);

      real8 pdot = V3_DOT(t1,t2);

      return ( (fabs((fabs(pdot)-1.0)) < 1.0e-4 ) ? 1 : 0 );
   }

   return(0);
}

void Paradis_Stat::Load (Paradis_Restart & rs)
{
   Recycle();

   pds_path   = ( rs.rs_path ? strdup(rs.rs_path) : 0 );

   node_cnt   = rs.rs_node_cnt;

   pds_lx     = rs.Lx();
   pds_ly     = rs.Ly();
   pds_lz     = rs.Lz();

   seg_cnt    = 0;
   Segment_t  *segs = rs.Segments(seg_cnt);

   Segment_List_Moment (segs, seg_cnt, seg_min, seg_max, seg_mean, seg_sum, seg_sdev);

   real8 bMag2   = bmag * bmag;
   real8 pds_vol = pds_lx * pds_ly * pds_lz;

   seg_density   = seg_sum/(bMag2*pds_vol);

   // populate the arms histogram

   nv_arms_max   = rs.Max_Arms();
   nv_arms_histo = ( (nv_arms_max>0) ? new int[nv_arms_max+1] : 0 );

   if (nv_arms_histo) 
   { 
      memset(nv_arms_histo,0, (nv_arms_max+1)*sizeof(int)); 

      for (int i=0; (i<rs.rs_node_cnt); i++)
      {
         int narms = rs.rs_nodes[i].numNbrs;
   
         if ( (0<=narms) && (narms<=nv_arms_max) )  { nv_arms_histo[narms]++; }
      }
   }

   real8 *bv = rs.Burgers_Vectors();

   bv_tbl_cnt  = rs.rs_bv_cnt;
   bv_tbl      = ( (bv && (bv_tbl_cnt>0)) ? new real8[3*bv_tbl_cnt] : 0 );

   if (bv_tbl) { memcpy(bv_tbl,bv,3*bv_tbl_cnt*sizeof(real8)); }  // (copy from restart class)

   bv_tbl_sums = rs.Burgers_Lengths();

   bv_tbl_cnts = ( (bv && (bv_tbl_cnt>0)) ? new int[bv_tbl_cnt] : 0 );

   if (bv_tbl_cnts) { memset(bv_tbl_cnts,0,bv_tbl_cnt*sizeof(int)); }

   if (segs && bv_tbl && bv_tbl_cnts)
   {
      for (int i=0; (i<seg_cnt); i++)
      {
         int bndx = BV_Index( segs[i].bv, bv_tbl, bv_tbl_cnt);
         if (bndx>=0) bv_tbl_cnts[bndx]++;
      }
   }

   //VDELETE(segs);
}


// BV_Sum()
//
// Returns the total length and count of all the segments with a given burgers vector.
//-----------------------------------------------------------------------------------------------

real8 Paradis_Stat::BV_Sum (int & cnt, const real8 *bv)
{
   real8 sum = 0.0;
         cnt = 0;

   if (bv_tbl && bv)
   {
      for (int i=0,j=0; (i<bv_tbl_cnt); i++,j+=3)
      {
         if (BV_Equiv(bv,bv_tbl+j))
         {
            cnt += bv_tbl_cnts[i];
            sum += bv_tbl_sums[i];
         }
   
      }
   }

   return(sum);
}

// BV_Sums()
//
// Given a table of burgers vectors, returns the sums and counts of all the equivalent
// burgers vector table entries.
//-----------------------------------------------------------------------------------------------

void Paradis_Stat::BV_Sums 
(
         real8 *bv_sums,   ///< returned sums of burgers vector lengths
         int   *bv_cnts,   ///< returned counts of burgers vector entries
   const real8 *bv     ,   ///< input table of burgers vectors
   const int    bv_n       ///< number of burgers vectors
)
{
   if (bv_sums && bv_cnts && bv)
   {
      memset(bv_cnts,0,bv_n*sizeof(int  ));
      memset(bv_sums,0,bv_n*sizeof(real8));

      for (int i=0,j=0; (i<bv_n); i++, j+=3)
         bv_sums[i] = BV_Sum(bv_cnts[i],bv+j);
   }
}

void Paradis_Stat::Print (FILE *fd)
{
   if (fd)
   {
      fprintf(fd,"ParaDiS Restart Statistics...\n");
      fprintf(fd,"\n");
      fprintf(fd,"   pds_path    = %s\n"    ,   pds_path    );
      fprintf(fd,"\n");
      fprintf(fd,"   node_cnt    = %d\n"    ,   node_cnt    );
      fprintf(fd,"\n");
      fprintf(fd,"   pds_lx      = %0.1lf\n",   pds_lx      );
      fprintf(fd,"   pds_ly      = %0.1lf\n",   pds_ly      );
      fprintf(fd,"   pds_lz      = %0.1lf\n",   pds_lz      );
      fprintf(fd,"\n");
      fprintf(fd,"   seg_cnt     = %d\n"    ,   seg_cnt     );
      fprintf(fd,"   seg_min     = %0.2lf\n",   seg_min     );
      fprintf(fd,"   seg_max     = %0.2lf\n",   seg_max     );
      fprintf(fd,"   seg_mean    = %0.2lf\n",   seg_mean    );
      fprintf(fd,"   seg_sdev    = %0.2lf\n",   seg_sdev    );
      fprintf(fd,"   seg_sum     = %0.2lf\n",   seg_sum     );
      fprintf(fd,"   seg_density = %0.6le\n",   seg_density );

      if (nv_arms_histo) 
      { 
         fprintf(fd,"\n");
         fprintf(fd,"  node arm distribution...\n");
     
         for (int i=0; (i<=nv_arms_max); i++)
            fprintf(fd,"    %2d : %d\n", i, nv_arms_histo[i] );
      } 

      if (bv_tbl)
      {
         fprintf(fd,"\n");
         fprintf(fd,"  burgers vector table...\n");

         real8 s = ( (bv_scale>0.0) ? 1.0/bv_scale : 1.0 );
         real8 t = (bmag*bmag)*(pds_lx*pds_ly*pds_lz);
               t = ( (t>0.0) ? 1.0/t : 1.0 );
     
         for (int i=0,j=0; (i<bv_tbl_cnt); i++, j+=3)
            fprintf(fd,"    %2d : %9.6lf %9.6lf %9.6lf  %8d %10.2lf   %0.6le\n", i, s*bv_tbl[j], s*bv_tbl[j+1], s*bv_tbl[j+2], bv_tbl_cnts[i], bv_tbl_sums[i], t*bv_tbl_sums[i] );
      }

      fprintf(fd,"\n");
   }
}

static void Save_PDS_Info 
(
   FILE         *fd      ,  ///< output file descriptor
   Paradis_Stat *pds     ,  ///< array of paradis stat structs
   int           pds_cnt    ///< number of stats
)
{
   if (fd && pds && (pds_cnt>0) )
   {
      fprintf(fd,"# paradis file statistics...\n" ); 
      fprintf(fd,"#\n");
      fprintf(fd,"# columns : fname lx ly lz nodes segs min max mean sdev sum density\n");
      fprintf(fd,"#--------------------------------------------------------------------------------------------\n" ); 

      for (int i=0; (i<pds_cnt); i++)
      {
         fprintf(fd,"%s "                                       , pds[i].pds_path );
         fprintf(fd,"%8.1lf %8.1lf %8.1lf "                     , pds[i].pds_lx, pds[i].pds_ly, pds[i].pds_lz );
         fprintf(fd,"%6d %6d "                                  , pds[i].node_cnt, pds[i].seg_cnt );
         fprintf(fd,"%6.2lf %6.2lf %6.2lf %6.2lf %10.2lf %0.6le", pds[i].seg_min, pds[i].seg_max, pds[i].seg_mean, pds[i].seg_sdev, pds[i].seg_sum, pds[i].seg_density );
         fprintf(fd,"\n");
      }
   }
}

static void Save_PDS_GNU 
(
   FILE         *fd      ,  ///< output file descriptor
   char         *path    ,  ///< base path
   Paradis_Stat *pds     ,  ///< array of paradis stat structs
   int           pds_cnt    ///< number of stats
)
{
   if (fd && pds && (pds_cnt>0) )
   {
      fprintf(fd,"#set term x11\n"  );
      fprintf(fd,"#set term pngcairo size 1000,800\n");
      fprintf(fd,"reset\n" );
      fprintf(fd,"set title 'ParaDiS Node/Segment Counts'\n");
      fprintf(fd,"plot '%s.dat' u  5  w l lt 1 lw 1 lc rgb '#0000ff' title 'nodes'    , \\\n", path );
      fprintf(fd,"     '%s.dat' u  6  w l lt 1 lw 1 lc rgb '#ff00ff' title 'segments' ; \n"  , path );
      fprintf(fd,"pause -1;\n" );
      fprintf(fd,"set ylabel '(b)' \n" );
      fprintf(fd,"set title 'ParaDiS Minimum Segment Length (b)'\n");
      fprintf(fd,"plot '%s.dat' u  7  w l lt 1 lw 1 lc rgb '#ff0000' notitle         ; \n"  , path );
      fprintf(fd,"pause -1;\n" );
      fprintf(fd,"set ylabel '(b)' \n" );
      fprintf(fd,"set title 'ParaDiS Maximum Segment Length (b)'\n");
      fprintf(fd,"plot '%s.dat' u  8  w l lt 1 lw 1 lc rgb '#ff0000' notitle         ; \n"  , path );
      fprintf(fd,"pause -1;\n" );
      fprintf(fd,"set ylabel '(b)' \n" );
      fprintf(fd,"set title 'ParaDiS Mean Segment Length (b)'\n");
      fprintf(fd,"plot '%s.dat' u  9  w l lt 1 lw 1 lc rgb '#ff0000' notitle         ; \n"  , path );
      fprintf(fd,"pause -1;\n" );
      fprintf(fd,"set ylabel '(b)' \n" );
      fprintf(fd,"set title 'ParaDiS Total Segment Length (b)'\n");
      fprintf(fd,"plot '%s.dat' u 11  w l lt 1 lw 1 lc rgb '#ff0000' notitle         ; \n"  , path );
      fprintf(fd,"pause -1;\n" );
   }
}

//----------------------------------------------------------------------------------------------------------

int main (int argc, char *argv[])
{
   Process_Args(argc,argv);

   FILE *fd_density = ( density_path ? fopen(density_path,"w")  : 0 );

   if ( File_Type(src_path) == 1)
   {
      Paradis_Restart rs; 
      rs.Load(src_path);

      Paradis_Stat pds;
      pds.Load (rs);
      pds.Print();

#ifdef  PRINT_BV_SUMS
      pds.BV_Sums(bv_ptbl_sums, bv_ptbl_cnts, bv_ptbl, bv_ptbl_cnt);

      real8 s = ( (bv_scale>0.0) ? 1.0/bv_scale : 1.0 );
      real8 t = (bmag*bmag)*(pds.pds_lx*pds.pds_ly*pds.pds_lz);
            t = ( (t>0.0) ? 1.0/t : 1.0 );
  
      printf("   ptable densities...\n");

      for (int i=0,j=0; (i<bv_ptbl_cnt); i++,j+=3)
      {
         printf("    %2d : %9.6lf %9.6lf %9.6lf  %8d %10.2lf   %0.6le\n", i, s*bv_ptbl[j], s*bv_ptbl[j+1], s*bv_ptbl[j+2], bv_ptbl_cnts[i], bv_ptbl_sums[i], t*bv_ptbl_sums[i] );
      }
#endif
   }

   if ( File_Type(src_path) == 2)
   {
      Directory dir;
   
      src_paths    = dir.Find(src_path,".data");
      src_path_cnt = dir.dir_cnt;

      if ( src_paths && (src_path_cnt>0) )
      {
         Paradis_Stat *pds = new Paradis_Stat[src_path_cnt];

         for (int i=0; (i<src_path_cnt); i++)
         {
            printf("processing %s (%2d%%)\n", src_paths[i], (i*100)/src_path_cnt ); 

            Paradis_Restart rs; 
            rs.Load(src_paths[i]);

            pds[i].Load(rs);

#ifdef  PRINT_BV_SUMS
            if (fd_density)
            {
               pds[i].BV_Sums(bv_ptbl_sums, bv_ptbl_cnts, bv_ptbl, bv_ptbl_cnt);

               real8 t = (bmag*bmag)*(pds[i].pds_lx*pds[i].pds_ly*pds[i].pds_lz);
                     t = ( (t>0.0) ? 1.0/t : 1.0 );
  
               fprintf(fd_density, "%20s :  ", src_paths[i] );

               for (int i=0; (i<bv_ptbl_cnt); i++)
               {
                  fprintf(fd_density," %14.6le",  t*bv_ptbl_sums[i] );
               }

               fprintf(fd_density, "\n");
            }
#endif
         }

         if (dst_path)
         {
            char  path[256]; 
            FILE *fd=0;

            sprintf(path,"%s.dat", dst_path); if ( (fd=fopen(path,"w")) ) { Save_PDS_Info(fd,         pds,src_path_cnt);  fclose(fd); } 
            sprintf(path,"%s.gnu", dst_path); if ( (fd=fopen(path,"w")) ) { Save_PDS_GNU (fd,dst_path,pds,src_path_cnt);  fclose(fd); }
         }
         else
            Save_PDS_Info(stdout,pds,src_path_cnt);
      }
      else
      {
         printf("no restart files found in %s\n", src_path);
      }
   }

   if (fd_density) fclose(fd_density);

   exit(0);
   return(0);
}
