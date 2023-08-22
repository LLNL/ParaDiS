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

real8  bv_scale     = 1.0;

//----------------------------------------------------------------------------------------------------------

const char  usage[] =
   "                                                                                              \n"
   " The ParaDiS bvdump application will process one or more restart files dumping out a list of  \n"
   " segment burgers vectors and normalized line directions.                                      \n"
   "                                                                                              \n"
   " usage: paradis_bvdump <path> -bs <scale> -dst <path>                                         \n"
   "                                                                                              \n"
   " where                                                                                        \n"
   "    <path>         identifies a single restart file or directory containing a collection      \n"
   "                   collection of restart files (required).                                    \n"
   "                                                                                              \n"
   "    -bs   <scale>  burgers vector scaling factor (default=1.0)                                \n"
   "                                                                                              \n"
   "    -dst  <path>   path to the output file (or directory)                                     \n"
   "                    - if no output path is specified, output will be stdout                   \n"
   "                    - if the source path is a file, output will be a file                     \n"
   "                    - if the source path is a directory, output will be a sequence of files   \n"
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

   if ( (indx=Find_Arg("-bs" ,argc,argv))>0 ) { bv_scale = (real8) strtod(argv[indx+1],0); }
   if ( (indx=Find_Arg("-dst",argc,argv))>0 ) { dst_path =         strdup(argv[indx+1]);   }
}

//----------------------------------------------------------------------------------------------------------

static void Segment_BV_Dump (FILE *fd, Paradis_Restart & rs)
{
   int        nseg = 0;
   Segment_t *segs = rs.Segments(nseg);

   if (fd && segs && (nseg>0))
   {
      fprintf(fd,  "# t1        t2              bv(x)        bv(y)        bv(z)        t(x)         t(y)         t(z)       dot(bv,t)\n" );

      Segment_t *seg = segs;

      for (int i=0; (i<nseg); i++,seg++)
      {
         char t1_str[32]; seg->node1->myTag.Print(t1_str);
         char t2_str[32]; seg->node2->myTag.Print(t2_str);

         real8 p1[3] = { seg->p1[0], seg->p1[1], seg->p1[2] };
         real8 p2[3] = { seg->p2[0], seg->p2[1], seg->p2[2] };
         real8 bv[3] = { seg->bv[0], seg->bv[1], seg->bv[2] };

         real8 pv[3] = { p2[0]-p1[0], 
                         p2[1]-p1[1], 
                         p2[2]-p1[2] };

         V3_NORMALIZE(bv,bv);
         V3_NORMALIZE(pv,pv);

         real8 dot = V3_DOT(bv,pv);

         V3_SCALE(bv,bv,bv_scale);

         fprintf(fd," %-10s %-10s"            , t1_str, t2_str      );
         fprintf(fd," %12.8lf %12.8lf %12.8lf", bv[0], bv[1], bv[2] );
         fprintf(fd," %12.8lf %12.8lf %12.8lf", pv[0], pv[1], pv[2] );
         fprintf(fd," %12.8lf"                , dot                 );
         fprintf(fd,"\n");
      }
   }

   if (segs) { delete [] segs; segs=0; }
}

static void Segment_BV_Dump (const char *path, Paradis_Restart & rs)
{
   FILE *fd = ( path ? fopen(path,"w") : 0 );

   if (fd) { Segment_BV_Dump(fd,rs); fclose(fd); }
}

//----------------------------------------------------------------------------------------------------------

static char *File_Basename (const char *path)
{
   static char base[256];

   memset(base,0,sizeof(base));

   if (path)
   {
      char *p = strrchr((char *) path,'/');
            p = ( p ? (p+1) : (char *) path );

      if (p) { strcpy(base,p); }

      p = strrchr(base,'.');

      if (p) { *p=0; }
   }

   return(base);
}

//----------------------------------------------------------------------------------------------------------

int main (int argc, char *argv[])
{
   Process_Args(argc,argv);

   if ( File_Type(src_path) == 1)
   {
      Paradis_Restart rs(src_path); 

      if (dst_path) { Segment_BV_Dump(dst_path,rs); }
      else          { Segment_BV_Dump(stdout  ,rs); }
   }

   if ( File_Type(src_path) == 2)
   {
      Directory dir;
   
      char **src_paths    = dir.Find(src_path,".data");
      int    src_path_cnt = dir.dir_cnt;

      if ( src_paths && (src_path_cnt>0) )
      {
         for (int i=0; (i<src_path_cnt); i++)
         {
            printf("processing %s (%2d%%)\n", src_paths[i], (i*100)/src_path_cnt ); 

            Paradis_Restart rs(src_paths[i]);

            if (dst_path)
            {
               char *fname = File_Basename(src_paths[i]);

               char  path[256];
               sprintf(path,"%s/%s.bv", dst_path, fname); 
   
               Segment_BV_Dump(path,rs);
            }
            else
               Segment_BV_Dump(stdout,rs);
         }
      }
      else
      {
         printf("no restart files found in %s\n", src_path);
      }
   }

   exit(0);
   return(0);
}
