#pragma once

#ifndef _PDS_ARGS_H
#define _PDS_ARGS_H

//--------------------------------------------------------------------------
// Generic application argument support
//--------------------------------------------------------------------------

extern int    Find_Arg (const char *str, int argc, char *argv[]);

extern int    Get_Flag (const char *str, int argc, char *argv[]);

extern int    Get_Arg  (const char *str, int argc, char *argv[], int    & val);
extern int    Get_Arg  (const char *str, int argc, char *argv[], float  & val);
extern int    Get_Arg  (const char *str, int argc, char *argv[], double & val);
extern char  *Get_Arg  (const char *str, int argc, char *argv[]);

extern int    Get_Args (const char *str, int argc, char *argv[], int     *v, int vn);
extern int    Get_Args (const char *str, int argc, char *argv[], float   *v, int vn);
extern int    Get_Args (const char *str, int argc, char *argv[], double  *v, int vn);

#endif
