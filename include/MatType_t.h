#pragma once

#ifndef _PDS_MATTYPE_T_H
#define _PDS_MATTYPE_T_H

//----------------------------------------------------------------------
// Define some values for indicating the type of crystal
// structure the material has (determined by the mobility
// function selected).  Used mainly in sections of the
// code that behave differently for BCC and FCC materials.
// By grouping things this way we don't have to check for every
// single BCC or FCC mobility in those areas, and in
// particular don't have to modify a dozen pieces of
// code to add a new mobility function.
//----------------------------------------------------------------------

typedef enum 
{
   MAT_TYPE_UNDEFINED        = -1,
   MAT_TYPE_BCC              =  0,
   MAT_TYPE_FCC                  ,
   MAT_TYPE_HCP                  ,
   MAT_TYPE_RHOMBOHEDRAL_VA      ,
   MAT_TYPE_MAX_INDEX              // note - MUST be the last item in this list
} MatType_t;

#define MAT_TYPE_NAME_UNDEFINED "UNDEFINED"
#define MAT_TYPE_NAME_BCC       "BCC"
#define MAT_TYPE_NAME_FCC       "FCC"
#define MAT_TYPE_NAME_HCP       "HCP"
#define MAT_TYPE_NAME_RHOMBO_VA "RHOMBO_VA"

#endif  // _PDS_MATTYPE_T_H
