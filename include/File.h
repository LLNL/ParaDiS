#pragma once

#ifndef _PDS_FILE_H
#define _PDS_FILE_H

//----------------------------------------------------------------------------------------------------------

class File
{
   public :
      char    *fbuf  ;   ///< points to the loaded file
      size_t   fbytes;   ///< size of the loaded file

   public :
      File(void);
     ~File();

     size_t Bytes (void) { return( fbuf ? fbytes : 0 ); }

     char  *Load   (const char *path);

     void   Save   (const char *path, const char *buf, const size_t bytes);
     void   Save   (const char *path, const char *buf);

     void   Append (const char *path, const char *buf, const size_t bytes);
     void   Append (const char *path, const char *buf);
};

//----------------------------------------------------------------------------------------------------------

extern size_t   File_Size      (const char *path);

extern void     File_Save      (const char *path, const char *buf, const size_t bytes);
extern void     File_Save      (const char *path, const char *str);
extern void     File_Append    (const char *path, const char *buf, const size_t bytes);
extern void     File_Append    (const char *path, const char *str);

extern int      File_Type      (const char *path);
extern int      File_Exists    (const char *path);
extern int      File_Delete    (const char *path);

extern int      IsA_File       (const char *path);
extern int      IsA_Link       (const char *path);
extern int      IsA_Directory  (const char *path);

#endif
