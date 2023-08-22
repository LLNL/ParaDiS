#pragma once

#ifndef _PDS_DIRECTORY_H
#define _PDS_DIRECTORY_H

#include <dirent.h>

//----------------------------------------------------------------------------------------------------------

class Directory
{
   public :
      int     dir_cnt   ;   ///< current   number of paths in the directory
      int     dir_max   ;   ///< allocated number of paths in the directory
      char  **dir_paths ;   ///< current array of paths to individual files

   public :
      Directory(void);
     ~Directory();

     void   Recycle (void);

     int    Count   (void) const { return(dir_cnt  ); }
     char **Paths   (void) const { return(dir_paths); }

     char **Load    (const char *path);
     char **Load    (const char *path, const char *sfx);
     char **Load    (const char *path, int (*cb)(const struct dirent *d) );

     char **Find    (const char *path);
     char **Find    (const char *path, const char *sfx);
     char **Find    (const char *path, int (*cb)(const struct dirent *d) );

     void   Sort    (void);

     void   Print   (FILE *fd=stdout);
     void   Print   (const char *path);
};

#endif  // _PDS_DIRECTORY_H
