//--------------------------------------------------------------------------
// Generic application argument support
//--------------------------------------------------------------------------

#include <string.h>
#include <stdlib.h>

#include "Args.h"

//--------------------------------------------------------------------------

int Find_Arg (const char *str, int argc, char *argv[])
{
   if (str && argv && (argc>0) )
   {
      for (int i=0; (i<argc); i++)
         if (strcmp(str,argv[i])==0) return(i);
   }

   return(0);
}

//--------------------------------------------------------------------------

int Get_Flag (const char *str, int argc, char *argv[])
{
   if (str && argv && (argc>0) )
   {
      for (int i=0; (i<argc); i++)
         if (strcmp(str,argv[i])==0) return(1);
   }

   return(0);
}

//--------------------------------------------------------------------------

int Get_Arg (const char *str, int argc, char *argv[], int & val)
{
   int indx = Find_Arg(str,argc,argv);

   if ( (indx>0) && (indx<(argc-1)) )
   {
      val = (int) strtol(argv[indx+1],0,10);
      return(1);
   }

   return(0);
}

//--------------------------------------------------------------------------

int Get_Arg (const char *str, int argc, char *argv[], float & val)
{
   int indx = Find_Arg(str,argc,argv);

   if ( (indx>0) && (indx<(argc-1)) )
   {
      val = strtof(argv[indx+1],0);
      return(1);
   }

   return(0);
}

//--------------------------------------------------------------------------

int Get_Arg (const char *str, int argc, char *argv[], double & val)
{
   int indx = Find_Arg(str,argc,argv);

   if ( (indx>0) && (indx<(argc-1)) )
   {
      val = strtod(argv[indx+1],0);
      return(1);
   }

   return(0);
}

//--------------------------------------------------------------------------

char *Get_Arg (const char *str, int argc, char *argv[])
{
   int indx = Find_Arg(str,argc,argv);

   if ( (indx>0) && (indx<(argc-1)) )
      return( strdup(argv[indx+1]) );

   return(0);
}

//--------------------------------------------------------------------------

int Get_Args (const char *str, int argc, char *argv[], int     *v, int vn)
{
   int indx = Find_Arg(str,argc,argv);

   if ( (indx>0) && (indx<(argc-vn)) && v)
   {
      for (int i=0; (i<vn); i++)
         v[i] = (int) strtol(argv[indx+1+i],0,10);

      return(1);
   }

   return(0);
}

//--------------------------------------------------------------------------

int Get_Args (const char *str, int argc, char *argv[], float   *v, int vn)
{
   int indx = Find_Arg(str,argc,argv);

   if ( (indx>0) && (indx<(argc-vn)) && v)
   {
      for (int i=0; (i<vn); i++)
         v[i] = strtof(argv[indx+1+i],0);

      return(1);
   }

   return(0);
}

//--------------------------------------------------------------------------

int Get_Args (const char *str, int argc, char *argv[], double  *v, int vn)
{
   int indx = Find_Arg(str,argc,argv);

   if ( (indx>0) && (indx<(argc-vn)) && v)
   {
      for (int i=0; (i<vn); i++)
         v[i] = strtod(argv[indx+1+i],0);

      return(1);
   }

   return(0);
}

