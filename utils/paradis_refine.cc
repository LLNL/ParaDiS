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

char  *src_path  = 0;
char  *dst_path  = 0;

real8  smax      = 0.0;

//----------------------------------------------------------------------------------------------------------

const char  usage[] =
   "                                                                                         \n"
   " The ParaDiS refine application will take an existing ParaDiS node configuration file    \n"
   " and refine the mesh such that no segments exceed the specified maximum segment length.  \n"
   "                                                                                         \n"
   " usage: paradis_refine <path> -max <len>                                                 \n"
   "                                                                                         \n"
   " where                                                                                   \n"
   "    -src  <path>    identifies the path to the source    ParaDiS restart file [required] \n"
   "    -dst  <path>    identifies the path to the resulting ParaDiS restart file [required] \n"
   "    -max  <len>     specifies the maximem segment length (b)                  [required] \n"
   "                                                                                         \n"
   "\n";

//----------------------------------------------------------------------------------------------------------

static void Process_Args (int argc, char *argv[])
{
   int indx=0;

   if ( Get_Flag("-help"     ,argc,argv) )    { fprintf(stderr,"%s", usage ); exit(0); }
   if ( Get_Flag("--help"    ,argc,argv) )    { fprintf(stderr,"%s", usage ); exit(0); }
   if ( argc==1 )                             { fprintf(stderr,"%s", usage ); exit(0); }

   if ( (indx=Find_Arg("-src",argc,argv))>0 ) { src_path = (char *) strdup(argv[indx+1]); }
   if ( (indx=Find_Arg("-dst",argc,argv))>0 ) { dst_path = (char *) strdup(argv[indx+1]); }
   if ( (indx=Find_Arg("-max",argc,argv))>0 ) { smax = (real8) strtod(argv[indx+1],0); }

   if ( !src_path || !dst_path || (smax<=0.0) ) { fprintf(stderr,"%s", usage ); exit(0); }
}

//----------------------------------------------------------------------------------------------------------

int main (int argc, char *argv[])
{
   Process_Args(argc,argv);

   Paradis_Restart rs(src_path);

   rs.Refine(smax);
   rs.Save(dst_path);

   exit(0);
   return(0);
}
