#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>

#include "Args.h"
#include "File.h"
#include "Directory.h"
#include "Restart.h"

//----------------------------------------------------------------------------------------------------------

char *src_path = 0;
char *gnu_path = 0;

//----------------------------------------------------------------------------------------------------------

static void Help(void)
{
   const char  msg[] =
      "                                                                                     \n"
      " usage: paradis_gnu  -src inputfile  -gnu outfile(s)                                 \n"
      "                                                                                     \n"
      " where                                                                               \n"
      "    -src     <path>        Specifies the path of the restart data file (required).   \n"
      "                                                                                     \n"
      "    -gnu     <path>        Specifies the path of the output gnuplot files            \n"
      "                             <path>.gnu (command)                                    \n" 
      "                             <path>.dat (data)                                       \n" 
      "                                                                                     \n"
      "\n";

   fprintf(stderr,"%s",msg);

   return;
}

//----------------------------------------------------------------------------------------------------------

static void Process_Args (int argc, char *argv[])
{
   int indx=0;

   if (      Get_Flag("-help"  ,argc,argv) ) { Help(); exit(0); }
   if (      Get_Flag("--help" ,argc,argv) ) { Help(); exit(0); }

   if ((indx=Find_Arg("-src"   ,argc,argv))>0 ) { src_path = strdup(argv[indx+1]); }
   if ((indx=Find_Arg("-gnu"   ,argc,argv))>0 ) { gnu_path = strdup(argv[indx+1]); }

   if (!src_path) { printf("error - you must provide a source restart file/dir\n"); exit(0); }

   if (!gnu_path) { gnu_path = strdup(src_path); }
}

//----------------------------------------------------------------------------------------------------------

void Save_Animation_Script (const char *pbase, int i0, int i1)
{
   char  path[256]; sprintf(path,"%s_anim.gnu", pbase);

   FILE *fd = fopen(path,"w");

   i0 = 1000*( (i0/1000) - ( (i0%1000) ? 1 : 0) );
   i1 = 1000*( (i1/1000) + ( (i1%1000) ? 1 : 0) );
  
   if (fd)
   {
      fprintf(fd,"#set term x11\n");
      fprintf(fd,"#set term pngcairo size 1000,800\n");
      fprintf(fd,"\n");
      fprintf(fd,"#reset\n");
      fprintf(fd,"\n");
      fprintf(fd,"#set view 60,30, 1.2, 1.2\n");
      fprintf(fd,"#set view 90, 0, 1.2, 1.2\n");
      fprintf(fd,"\n");
      fprintf(fd,"i0=%d\n", i0);
      fprintf(fd,"i1=%d\n", i1);
      fprintf(fd,"\n");
      fprintf(fd,"set xrange [i0:i1]\n");
      fprintf(fd,"set yrange [i0:i1]\n");
      fprintf(fd,"set zrange [i0:i1]\n");
      fprintf(fd,"\n");
      fprintf(fd,"unset key\n");
      fprintf(fd,"\n");
      fprintf(fd,"set title 'ParaDiS Simulation Animation'\n");
      fprintf(fd,"\n");
      fprintf(fd,"set xlabel 'X'  \n");
      fprintf(fd,"set ylabel 'Y'  \n");
      fprintf(fd,"set zlabel 'Z'  \n");
      fprintf(fd,"set ticslevel 0 \n");
      fprintf(fd,"set xtics autofreq i0, 1000, i1\n");
      fprintf(fd,"set ytics autofreq i0, 1000, i1\n");
      fprintf(fd,"set ztics autofreq i0, 1000, i1\n");
      fprintf(fd,"\n");
      fprintf(fd,"files = system(\"ls -1 %s????.dat\")\n", pbase);
      fprintf(fd,"\n");
      fprintf(fd,"do for [i=0:words(files)-1]           \\\n");
      fprintf(fd,"{                                     \\\n");
      fprintf(fd,"   print sprintf('%s%%04d',i);         \\\n", pbase);
      fprintf(fd,"                                      \\\n");
      fprintf(fd,"   dat = sprintf('%s%%04d.dat', i);    \\\n", pbase);
      fprintf(fd,"   dom = sprintf('%s%%04d.dom', i);    \\\n", pbase);
      fprintf(fd,"                                      \\\n");
      fprintf(fd,"   splot dat  u 5:6:7:11:12:13     w vec nohead lt 1 lw 1 lc rgb '#0000ff' ; \\\n");
      fprintf(fd,"               \\\n");
      fprintf(fd,"   pause 0.1;  \\\n");
      fprintf(fd,"}\n");
      fprintf(fd,"# you can add these lines to the splot to add node dots and domain...\n");
      fprintf(fd,"#        dat  u 5:6:(i0):11:12:(0) w vec nohead lt 1 lw 1 lc rgb '#00aaaa' , \\\n");
      fprintf(fd,"#        dat  u 5:6:7              w points pt 6 ps 0.5   lc rgb '#0000ff' , \\\n");
      fprintf(fd,"#        dat  u 8:9:10             w points pt 6 ps 0.5   lc rgb '#0000ff' , \\\n");
      fprintf(fd,"#        dom  u 2:3:4:5:6:7        w vec nohead lt 1 lw 1 lc rgb '#aa0000' ; \\\n");
      fprintf(fd,"\n");
      fprintf(fd,"#pause -1\n");

      fclose(fd);
   }
}

//----------------------------------------------------------------------------------------------------------

void Process_Directory (void)
{
   Directory dir;

   char **paths = dir.Load(src_path,".data");
   int    cnt   = dir.dir_cnt;

   if (cnt>0)
   {
      for (int i=0; (i<cnt); i++)
      {
         Paradis_Restart rs;

         char path[256]; sprintf(path,"%s%04d", gnu_path, i);

         printf("processing %s\n", paths[i]);

         rs.Load     (paths[i]);
         rs.Print_GNU(path);

         if (i==0)
         {
            int i0 = rs.rs_min_coords[0];
            int i1 = rs.rs_max_coords[0];

            Save_Animation_Script(gnu_path,i0,i1);
         }
      }
   }
   else
      printf("no restart files found in %s\n", src_path);
   
}

//----------------------------------------------------------------------------------------------------------

int main (int argc, char *argv[])
{
   Process_Args(argc,argv);

   if ( File_Type(src_path) == 1)
   {
      Paradis_Restart rs;

      rs.Load     (src_path);
      rs.Print_GNU(gnu_path);
   }

   if ( File_Type(src_path) == 2)
   {
      Process_Directory();
   }

   exit(0);
   return(0);
}
