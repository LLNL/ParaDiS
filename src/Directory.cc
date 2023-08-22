
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>

#include "Directory.h"

static char **append (char **paths, int & cnt, int & max, const char *path)
{
   if ( !paths || (cnt==max) )
   {
      max = ( (max==0) ? 64 : 2*max );

      char **tmp = new char *[max];

      if (tmp && paths) { for (int i=0; (i<cnt); i++) tmp[i]=paths[i]; }

      if (paths) { delete [] paths; }

      paths = tmp;
   }
   
   if (paths && (cnt<max) ) { paths[cnt++] = strdup(path); }

   return(paths);
}

Directory::Directory(void)
{
   dir_cnt  =0;
   dir_max  =0;
   dir_paths=0;
}

Directory::~Directory()
{
   Recycle();
}

void Directory::Recycle (void)
{
   if (dir_paths)
   { 
      for (int i=0; (i<dir_cnt); i++)
         if (dir_paths[i]) { delete [] dir_paths[i]; dir_paths[i]=0; }

      delete [] dir_paths; dir_paths=0;
   } 

   dir_cnt  =0;
   dir_max  =0;
   dir_paths=0;
}

// Load()
//
// Will return an array of paths to regular files in a given directory
//----------------------------------------------------------------------------------------------------------

char **Directory::Load (const char *path)
{
   Recycle();

   DIR *dir = opendir(path); 

   if (dir)
   {
      for (struct dirent *dirent=readdir(dir); (dirent); dirent=readdir(dir) )
      {
         if (strcmp(dirent->d_name, "." )==0) continue;  // (ignore current dir entry)
         if (strcmp(dirent->d_name, "..")==0) continue;  // (ignore parent  dir entry)

         // if regular file - compose and append the path to the directory entries...

         if ( dirent->d_type==DT_REG )
         {
            char tmp[256];  sprintf(tmp,"%s/%s", path, dirent->d_name);

            dir_paths = append(dir_paths,dir_cnt,dir_max,tmp);
         }
      }

      closedir (dir);
   }

   Sort();

   return(dir_paths);
}

// Load()
//
// Will return an array of paths to regular files in a given directory with file names that match the 
// given file suffix.
//----------------------------------------------------------------------------------------------------------

char **Directory::Load (const char *path, const char *sfx)
{
   Recycle();

   DIR *dir = opendir(path); 

   if (dir)
   {
      for (struct dirent *dirent=readdir(dir); (dirent); dirent=readdir(dir) )
      {
         if (strcmp(dirent->d_name, "." )==0) continue;  // (ignore current dir entry)
         if (strcmp(dirent->d_name, "..")==0) continue;  // (ignore parent  dir entry)

         // if regular file - compose and append the path to the directory entries...

         char *p = strrchr(dirent->d_name,'.');

         if ( (dirent->d_type==DT_REG) && p && (strcmp(p,sfx)==0) ) 
         {
            char tmp[256];  sprintf(tmp,"%s/%s", path, dirent->d_name);

            dir_paths = append(dir_paths,dir_cnt,dir_max,tmp);
         }
      }

      closedir (dir);
   }

   Sort();

   return(dir_paths);
}

char **Directory::Load  (const char *path, int (*cb)(const struct dirent *d) )
{
   Recycle();

   DIR *dir = opendir(path); 

   if (dir)
   {
      for (struct dirent *dirent=readdir(dir); (dirent); dirent=readdir(dir) )
      {
         if (strcmp(dirent->d_name, "." )==0) continue;  // (ignore current dir entry)
         if (strcmp(dirent->d_name, "..")==0) continue;  // (ignore parent  dir entry)
         // if regular file - compose and append the path to the directory entries...

         if ( (dirent->d_type==DT_REG) && (*cb)(dirent) )
         {
            char tmp[256];  sprintf(tmp,"%s/%s", path, dirent->d_name);

            dir_paths = append(dir_paths,dir_cnt,dir_max,tmp);
         }
      }

      closedir (dir);
   }

   Sort();

   return(dir_paths);
}

// Find()
//
// Will recursively load the paths of all regular files in the given directory.
// Note - because of the recursive nature - you must recycle() this class before
// executing this method.
//----------------------------------------------------------------------------------------------------------

char **Directory::Find (const char *path)
{
   DIR *dir = opendir(path); 

   if (dir)
   {
      for (struct dirent *dirent=readdir(dir); (dirent); dirent=readdir(dir) )
      {
         if (strcmp(dirent->d_name, "." )==0) continue;  // (ignore current dir entry)
         if (strcmp(dirent->d_name, "..")==0) continue;  // (ignore parent  dir entry)

         // if regular file - compose and append the path to the directory entries...

         if ( dirent->d_type==DT_REG )
         {
            char tmp[256];  sprintf(tmp,"%s/%s", path, dirent->d_name);

            dir_paths = append(dir_paths,dir_cnt,dir_max,tmp);
         }

         if ( dirent->d_type==DT_DIR )
         {
            char tmp[256];  sprintf(tmp,"%s/%s", path, dirent->d_name);

            Find(tmp); 
         }

      }

      closedir (dir);
   }

   Sort();

   return(dir_paths);
}

char **Directory::Find (const char *path, const char *sfx)
{
   DIR *dir = opendir(path); 

   if (dir)
   {
      for (struct dirent *dirent=readdir(dir); (dirent); dirent=readdir(dir) )
      {
         if (strcmp(dirent->d_name, "." )==0) continue;  // (ignore current dir entry)
         if (strcmp(dirent->d_name, "..")==0) continue;  // (ignore parent  dir entry)

         // if regular file - compose and append the path to the directory entries...

         char *p = strrchr(dirent->d_name,'.');

         if ( (dirent->d_type==DT_REG) && p && (strcmp(p,sfx)==0) ) 
         {
            char tmp[256];  sprintf(tmp,"%s/%s", path, dirent->d_name);

            dir_paths = append(dir_paths,dir_cnt,dir_max,tmp);
         }

         if ( dirent->d_type==DT_DIR )
         {
            char tmp[256];  sprintf(tmp,"%s/%s", path, dirent->d_name);

            Find(tmp,sfx);  // (recurse)
         }

      }

      closedir (dir);
   }

   Sort();


   return(dir_paths);
}

char **Directory::Find (const char *path, int (*cb)(const struct dirent *d) )
{
   DIR *dir = opendir(path); 

   if (dir)
   {
      for (struct dirent *dirent=readdir(dir); (dirent); dirent=readdir(dir) )
      {
         if (strcmp(dirent->d_name, "." )==0) continue;  // (ignore current dir entry)
         if (strcmp(dirent->d_name, "..")==0) continue;  // (ignore parent  dir entry)

         // if regular file - compose and append the path to the directory entries...

         if ( (dirent->d_type==DT_REG) && (*cb)(dirent) )
         {
            char tmp[256];  sprintf(tmp,"%s/%s", path, dirent->d_name);

            dir_paths = append(dir_paths,dir_cnt,dir_max,tmp);
         }

         if ( dirent->d_type==DT_DIR )
         {
            char tmp[256];  sprintf(tmp,"%s/%s", path, dirent->d_name);

            Find(tmp);  // (recurse)
         }

      }

      closedir (dir);
   }

   Sort();


   return(dir_paths);
}


// Sort()
//
// Will lexically sort the paths in the current directory.
//----------------------------------------------------------------------------------------------------------

static int pcmp (const void *path0, const void *path1)
{
   char *p0 = ( path0 ? *((char **) path0) : 0 );
   char *p1 = ( path1 ? *((char **) path1) : 0 );

   return( (p0 && p1) ? strcmp(p0,p1) : 0 );
}

void Directory::Sort (void)
{
   if (dir_paths && (dir_cnt>0))
      qsort(dir_paths, dir_cnt, sizeof(char *), pcmp);
}

// Print()
//
// Will print all the files found in the current directory.
//----------------------------------------------------------------------------------------------------------

void Directory::Print (FILE *fd)
{
   if (fd && dir_paths && (dir_cnt>0) )
   {
      for (int i=0; (i<dir_cnt); ++i)
      {
         if (dir_paths[i]) { fprintf(fd,"%s\n", dir_paths[i]); }
      }
   }
}

void Directory::Print (const char *path)
{
   FILE *fd = ( path ? fopen(path,"w") : 0 );

   if (fd) { Print(fd); fclose(fd); }
}

