#pragma once

#ifndef _PDS_TAG_H
#define _PDS_TAG_H

#include "cuda_portability.h"

//-------------------------------------------------------------------------------------------------------------------
// Defines a structure the is used to uniquely identify a node as a combination 
// of MPI domain and local index within a domain.
//-------------------------------------------------------------------------------------------------------------------

class Tag_t
{
   public :
      int  domainID;    ///< domain which owns this node/tag
      int  index   ;    ///< index within the domain

   public :
      __cuda_hdev__   Tag_t(void)                          { domainID=0;          index=0;       } 
      __cuda_hdev__   Tag_t(const int dom, const int indx) { domainID=dom;        index=indx;    } 
      __cuda_hdev__   Tag_t(const Tag_t & t)               { domainID=t.domainID; index=t.index; }

      __cuda_hdev__  ~Tag_t() {}  // (nothing allocated)

      __cuda_hdev__  int Domain(void) const { return(domainID); }
      __cuda_hdev__  int Index (void) const { return(index   ); }

      __cuda_hdev__  const Tag_t & operator = (const Tag_t & t) { domainID=t.domainID; index=t.index; return(*this); }

      __cuda_hdev__  operator unsigned long () const            { return( (((unsigned long) domainID) << 32) | (unsigned long) index ); }

      __cuda_hdev__  int operator == (const Tag_t & t) const { return ( (((unsigned long) *this) == ((unsigned long) t) ) ? 1 : 0 ); }
      __cuda_hdev__  int operator != (const Tag_t & t) const { return ( (((unsigned long) *this) != ((unsigned long) t) ) ? 1 : 0 ); }
      __cuda_hdev__  int operator <  (const Tag_t & t) const { return ( (((unsigned long) *this) <  ((unsigned long) t) ) ? 1 : 0 ); }
      __cuda_hdev__  int operator <= (const Tag_t & t) const { return ( (((unsigned long) *this) <= ((unsigned long) t) ) ? 1 : 0 ); }
      __cuda_hdev__  int operator >  (const Tag_t & t) const { return ( (((unsigned long) *this) >  ((unsigned long) t) ) ? 1 : 0 ); }
      __cuda_hdev__  int operator >= (const Tag_t & t) const { return ( (((unsigned long) *this) >= ((unsigned long) t) ) ? 1 : 0 ); }

      __cuda_host__  char *Print(char *sbuf=0, const char *sfmt=0) const;
};

//-------------------------------------------------------------------------------------------------------------------

__cuda_host__ Tag_t *Tag_Find   (const Tag_t *tags, const int ntags, const Tag_t & tag );
__cuda_host__ int    Tag_Index  (const Tag_t *tags, const int ntags, const Tag_t & tag );

//-------------------------------------------------------------------------------------------------------------------

#endif  // (_PDS_TAG_H)
