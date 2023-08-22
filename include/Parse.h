#pragma once

#ifndef _PDS_PARSE_H
#define _PDS_PARSE_H

//---------------------------------------------------------------------------
// Module:      Parse.h  
// Description: Contains definitions, structures and prototypes
//              for code used in parsing the parameter files.
//---------------------------------------------------------------------------

#include "Home.h"

// Define a set values that may be returned by
// functions used in parsing the user supplied
// values in the control file.

#define TOKEN_ERR            -1
#define TOKEN_NULL            0
#define TOKEN_GENERIC         1
#define TOKEN_EQUAL           2
#define TOKEN_BEGIN_VAL_LIST  3
#define TOKEN_END_VAL_LIST    4

// Define the variable types that may be associated with input
// parameters

#define V_NULL    0
#define V_DBL     1
#define V_INT     2
#define V_STRING  3
#define V_COMMENT 4

#define VFLAG_NULL        0x00
#define VFLAG_ALIAS       0x01
#define VFLAG_INITIALIZED 0x02
#define VFLAG_DISABLED    0x04
#define VFLAG_SET_BY_USER 0x08

//-----------------------------------------------------------------------------------------------------------

extern void  BindVar               (ParamList_t *list, const char *name, const void *addr, const int type, const int cnt, const int flags);
extern void  DisableUnneededParams (Home_t      *home);
extern int   GetNextToken          (FILE        *fd  , char *token, int maxTokenSize);
extern int   GetParamVals          (FILE        *fd  , int   vType, int vExpected, void *vList);
extern int   LookupParam           (ParamList_t *list, const char *token);
extern void  MarkParamDisabled     (ParamList_t *list, const char *name );
extern void  MarkParamEnabled      (ParamList_t *list, const char *name );
extern void  WriteParam            (ParamList_t *list, const int   index, FILE *fd);

#endif // _PDS_PARSE_H
