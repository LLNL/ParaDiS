
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>

#include "File.h"

File::File(void)
{
   fbuf   = 0;
   fbytes = 0;
}

File::~File()
{
   if (fbuf) { delete [] fbuf; }

   fbuf   = 0;
   fbytes = 0;
}

char *File::Load (const char *path)
{
   if (fbuf) { delete [] fbuf; }

   fbuf   = 0;
   fbytes = 0;

   FILE *fp = ( path ? fopen(path,"r") : 0 );

   if (!fp) { printf("%s::%s() ln=%d error : failed to open file \"%s\"\n" , __FILE__, __func__, __LINE__, ( path ? path : "NIL" ) ); }

   if (fp)
   {
               fseek(fp,0,SEEK_END);
      fbytes = ftell(fp);
               fseek(fp,0,SEEK_SET);

      fbuf = ( (fbytes>0) ? new char[fbytes] : 0 );

      if (fbuf)
         fread(fbuf,fbytes,1,fp);

      fclose(fp);
   } 

   return(fbuf);
}

void File::Save (const char *path, const char *buf, const size_t bytes)
{
   if (!path ) { printf("%s::%s() ln=%d warning : null path provided to file save\n"       , __FILE__, __func__, __LINE__);  return; }
   if (!buf  ) { printf("%s::%s() ln=%d warning : null buffer provided to file save\n"     , __FILE__, __func__, __LINE__);  return; }
   if (!bytes) { printf("%s::%s() ln=%d warning : zero byte count provided to file save\n" , __FILE__, __func__, __LINE__);  return; }

   if (path && buf && (bytes>0))
   {
      FILE *fd = (FILE *) ( path ? fopen(path,"w") : 0 );

      if (!fd) { printf("%s::%s() ln=%d error : failed to open file \"%s\"\n" , __FILE__, __func__, __LINE__, path ); }

      if (fd)
      {
         fwrite(buf,bytes,1,fd);
         fclose(fd);
      }
   }
}

void File::Save (const char *path, const char *buf)
{
   Save(path,buf, buf ? strlen(buf) : 0 );
}

void File::Append (const char *path, const char *buf, const size_t bytes)
{
   if (!path ) { printf("%s::%s() ln=%d warning : null path provided to file save\n"       , __FILE__, __func__, __LINE__);  return; }
   if (!buf  ) { printf("%s::%s() ln=%d warning : null buffer provided to file save\n"     , __FILE__, __func__, __LINE__);  return; }
   if (!bytes) { printf("%s::%s() ln=%d warning : zero byte count provided to file save\n" , __FILE__, __func__, __LINE__);  return; }

   if (path && buf && (bytes>0))
   {
      FILE *fd = (FILE *) ( path ? fopen(path,"a") : 0 );

      if (!fd) { printf("%s::%s() ln=%d error : failed to open file \"%s\"\n" , __FILE__, __func__, __LINE__, path ); }

      if (fd)
      {
         fwrite(buf,bytes,1,fd);
         fclose(fd);
      }
   }
}

void File::Append (const char *path, const char *buf)
{
   Append(path,buf, buf ? strlen(buf) : 0 );
}

// File_Size()
//
// Given a path to a file, returns the size of that file in bytes.
//--------------------------------------------------------------------------------------

size_t File_Size(const char *path)
{
   size_t bytes = 0;
   FILE  *fp    = ( path ? fopen(path,"r") : 0 );

   if (!fp) { printf("%s::%s() ln=%d error : failed to open file \"%s\"\n" , __FILE__, __func__, __LINE__, ( path ? path : "NIL" ) ); }

   if ( fp)
   {
              fseek (fp,0,SEEK_END);    // seek to the end of the file
      bytes = ftell (fp);               // bytes = current file position
              fclose(fp);               // close the file
   } 

   return(bytes);
}

// File_Save()
//
// Given path to a file, a pointer to a buffer, and the size of the 
// buffer in bytes, will save the contents of that buffer to the output
// file.
//--------------------------------------------------------------------------------------

void File_Save (const char *path, const char *buf, const size_t bytes)
{
   if (!path ) { printf("%s::%s() ln=%d warning : null path provided to file save\n"       , __FILE__, __func__, __LINE__);  return; }
   if (!buf  ) { printf("%s::%s() ln=%d warning : null buffer provided to file save\n"     , __FILE__, __func__, __LINE__);  return; }
   if (!bytes) { printf("%s::%s() ln=%d warning : zero byte count provided to file save\n" , __FILE__, __func__, __LINE__);  return; }

   if (path && buf && (bytes>0))
   {
      FILE *fd = (FILE *) ( path ? fopen(path,"w") : 0 );

      if (!fd) { printf("%s::%s() ln=%d error : failed to open file \"%s\"\n" , __FILE__, __func__, __LINE__, path ); }

      if (fd)
      {
         fwrite(buf,bytes,1,fd);
         fclose(fd);
      }
   }
}

// (alternate form)

void File_Save (const char *path, const char *str) { File_Save(path, str, strlen(str) ); }

// File_Append()
//
// 
// Given path to a file, a pointer to a buffer, and the size of the 
// buffer in bytes, will append the contents of that buffer to the output
// file.
//--------------------------------------------------------------------------------------

void File_Append (const char *path, const char *buf, const size_t bytes)
{
   if (!path ) { printf("%s::%s() ln=%d warning : null path provided to file save\n"       , __FILE__, __func__, __LINE__);  return; }
   if (!buf  ) { printf("%s::%s() ln=%d warning : null buffer provided to file save\n"     , __FILE__, __func__, __LINE__);  return; }
   if (!bytes) { printf("%s::%s() ln=%d warning : zero byte count provided to file save\n" , __FILE__, __func__, __LINE__);  return; }

   if (path && buf && (bytes>0))
   {
      FILE *fd = (FILE *) ( path ? fopen(path,"a") : 0 );

      if (!fd) { printf("%s::%s() ln=%d error : failed to open file \"%s\"\n" , __FILE__, __func__, __LINE__, path ); }

      if (fd)
      {
         fwrite(buf,bytes,1,fd);
         fclose(fd);
      }
   }
}

// (alternate form)

void File_Append (const char *path, const char *str) { File_Append(path, (char *) str, strlen(str) ); }

//--------------------------------------------------------------------------------------

#if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS==64)
#define lstat lstat64
#define  stat  stat64
#endif

int File_Type (const char *path)
{
   struct stat fstat;

   if (path)
   {
      lstat(path, &fstat);

      if (S_ISLNK(fstat.st_mode))  { return(3); }
      else
      {
         stat(path, &fstat);
   
         if (S_ISREG(fstat.st_mode)) { return(1); }
         if (S_ISDIR(fstat.st_mode)) { return(2); }
      }
   }

   return (-1);
}

int IsA_File      (const char *path) { return( (File_Type(path)==1) ? 1 : 0 ); }
int IsA_Link      (const char *path) { return( (File_Type(path)==3) ? 1 : 0 ); }
int IsA_Directory (const char *path) { return( (File_Type(path)==2) ? 1 : 0 ); }

int File_Exists (const char *path)
{
   struct stat fstat;

   return ( path && (stat(path, &fstat)==0) ? 1 : 0 );
}

int File_Delete (const char *path)
{
   int err = ( path ? remove(path) : 0 );

   return(err);
}

//--------------------------------------------------------------------------------------

