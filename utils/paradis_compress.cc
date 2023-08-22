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
int   nthrd    = 1;

char *cmpr_cmd = strdup("bzip2");

//----------------------------------------------------------------------------------------------------------

const char  usage[] =
   "                                                                                              \n"
   " usage: paradis_compress <path> -d -n <cnt>                                                   \n"
   "                                                                                              \n"
   " where                                                                                        \n"
   "        <path>  identifies the directory containing a collection of restart files (required). \n"
   "                                                                                              \n"
   "    -n  <cnt>   process files in parallel using <cnt> threads.                                \n"
   "                   note - currently generates shell script output                             \n"
   "                                                                                              \n"
   "    -bzip2      compress files using bzip2    (default)                                       \n"
   "    -compress   compress files using compress                                                 \n"
   "                                                                                              \n"
   "\n";

//----------------------------------------------------------------------------------------------------------

static void Process_Args (int argc, char *argv[])
{
   if ( Get_Flag("-help"     ,argc,argv) ) { fprintf(stderr,"%s", usage ); exit(0); }
   if ( Get_Flag("--help"    ,argc,argv) ) { fprintf(stderr,"%s", usage ); exit(0); }
   if ( argc==1 )                          { fprintf(stderr,"%s", usage ); exit(0); }

   if ( Get_Arg ("-n"        ,argc,argv, nthrd) ) { nthrd = ( (nthrd< 1) ?  1 : nthrd );
                                                    nthrd = ( (nthrd>32) ? 32 : nthrd ); }

   if ( Get_Flag("-compress" , argc,argv) ) { cmpr_cmd = strdup("compress"); }

   src_path = strdup(argv[1]);
}

//----------------------------------------------------------------------------------------------------------

int main (int argc, char *argv[])
{
   Process_Args(argc,argv);

   if ( File_Type(src_path) == 1)
   {
      char cmd[256];  

      sprintf(cmd,"%s %s", cmpr_cmd, src_path );
      printf("%s\n", cmd); 

      system(cmd);
   }

   if ( File_Type(src_path) == 2)
   {
      Directory dir;
   
      char **paths = dir.Find(src_path,".data");
      int    cnt   = dir.dir_cnt;

      if (cnt==0)
      {
         printf("no restart files found in %s\n", src_path);
      }

      // if single threaded...

      if ( (cnt>0) && (nthrd==1) )
      {
         for (int i=0; (i<cnt); i++)
         {
            printf("compressing %s (%2d%%)\n", paths[i], (i*100)/cnt ); 

            char cmd[256];  sprintf(cmd,"%s %s", cmpr_cmd, paths[i]); system(cmd);
         }
      }

      // if multi threaded...
      // note that until I can get the waitpid() stuff to work, I just generate shell commands...
      
      if ( (cnt>0) && (nthrd>1) )
      {
         printf("#!/bin/sh -x\n");
         printf("#\n");

         for (int i=0; (i<cnt); i++)
         {
            char cmd[256];  
            sprintf(cmd,"%s %s &", cmpr_cmd, paths[i]);
           
            printf("%s\n", cmd); 

            if ( ((i+1)%nthrd)==0) printf("wait\n");
         }

         printf("wait\n");
      }
   }

   exit(0);
   return(0);
}
