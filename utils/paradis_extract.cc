#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Args.h"
#include "Restart.h"

/*-----------------------------------------------------------------------------------------------------------------
 * Globals...
 *-----------------------------------------------------------------------------------------------------------------*/

char *src_path = 0;
char *dst_path = 0;

double xpos[3] = { 0.0, 0.0, 0.0 };
double xwin[3] = { 0.0, 0.0, 0.0 };

/*-----------------------------------------------------------------------------------------------------------------*/

static void Help(void)
{
   const char  msg[] =
      "                                                                                            \n"
      " This application will extract the segments from a fixed region of an existing              \n"
      " ParaDiS restart file and generate a restart file containing only those nodes               \n"
      " within the specified region.                                                               \n"
      "                                                                                            \n"
      " usage: paradis_extract  -src inputfile  -dst outfile                                       \n"
      "                         -c x y z  -w wx wy wz                                              \n"
      "                                                                                            \n"
      " where                                                                                      \n"
      "    -src     <path>        Specifies the path of the restart data file (required).          \n"
      "    -dst     <path>        Specifies the path of the extracted restart file                 \n"
      "    -c       x y z         Specifies the center of the extraction window                    \n"
      "    -w       wx wy wz      Specifies the width, height, and depth of the extraction window  \n"
      "                                                                                            \n"
      "\n";

   fprintf(stderr,"%s",msg);

   return;
}

/*-----------------------------------------------------------------------------------------------------------------*/

static void Process_Args (int argc, char *argv[])
{
   int indx=0;

   if (      Get_Flag("-help"  ,argc,argv) ) { Help(); exit(0); }
   if (      Get_Flag("--help" ,argc,argv) ) { Help(); exit(0); }

   if ((indx=Find_Arg("-src"   ,argc,argv))>0 ) { src_path = strdup(argv[indx+1]); }
   if ((indx=Find_Arg("-dst"   ,argc,argv))>0 ) { dst_path = strdup(argv[indx+1]); }

   if (!src_path) { printf("error - you must provide a source restart file\n"); exit(0); }
   if (!dst_path) { printf("error - you must provide a destination file\n"   ); exit(0); }

   if ( !(Get_Args("-c",argc,argv,xpos,3)) ) { xpos[0]=xpos[1]=xpos[2]=   0.0; }
   if ( !(Get_Args("-w",argc,argv,xwin,3)) ) { xwin[0]=xwin[1]=xwin[2]=1000.0; }
}

/*-----------------------------------------------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
   Process_Args(argc,argv);

   Paradis_Restart rs(src_path);
   Paradis_Restart rx;

   rx.Extract(rs,xpos,xwin);
   rx.Save   (dst_path);

   exit(0);
   return(0);
}
