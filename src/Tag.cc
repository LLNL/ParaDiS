#include <stdio.h>

#include "Tag.h"
#include "TagMap.h"

// Print()
//
// Will convert the tag to a human-readable, printable string. Note that if the buffer is provided to 
// the method, will print into the buffer.  If the buffer is NOT provided, the string is printed into
// a locally static buffer. This routine is only thread-safe if a buffer is provided as an argument to
// the method.
//-------------------------------------------------------------------------------------------------------------------

__cuda_host__
char *Tag_t::Print(char *sbuf, const char *sfmt) const
{
   static char  tmp[32];
          char *str = ( sbuf ? (char *) sbuf : tmp );
          char *fmt = ( sfmt ? (char *) sfmt : (char *) "(%d,%d)" );

   sprintf(str,fmt,domainID,index);

   return(str); 
}

// Tag_Find()
//
// Given an array of sorted tags, returns a pointer to the tag containing the key.
// Returns NIL if the tag is not found.
//-------------------------------------------------------------------------------------------------------------------

__cuda_host__
Tag_t *Tag_Find
(
   const Tag_t  *tags ,  ///< array of tags
   const int     ntags,  ///< number of tags
   const Tag_t & key     ///< the search tag
)
{
   if (tags && (ntags>0))
   {
      int i0=0;
      int i1=(ntags-1);
   
      while (i0<=i1)
      {   
         int i = (i0+i1)/2;

         int  cmp  = ( (tags[i] < key) ? -1 : 
                       (tags[i] > key) ? +1 : 0 );

         if (cmp==0) { return( (Tag_t *)(tags+i) ); }
         if (cmp< 0) { i0=(i+1); }
         else        { i1=(i-1); }
      }
   }

   return(0);  // (not found)
}

// Tag_Index()
//
// Given an array of sorted tags, returns index of key tag.
// Returns -1 if the tag is not found.
//-------------------------------------------------------------------------------------------------------------------

__cuda_host__
int Tag_Index
(
   const Tag_t  *tags ,  ///< array of tags
   const int     ntags,  ///< number of tags
   const Tag_t & key     ///< the search tag
)
{
   Tag_t *ptag = Tag_Find(tags,ntags,key);

   return ( ptag ? (int)(ptag-tags) : -1 );
}

