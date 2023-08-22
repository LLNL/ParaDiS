
#include "TagMap.h"

//-------------------------------------------------------------------------------------------------------------------
// Note - constructors and class inlined in header
//-------------------------------------------------------------------------------------------------------------------

// TagMap_Append()
//
// Will append a tagmap pair to an existing tag map list (array).  
// Reallocates and extends the array as necessary. 
//
// Note - this routine is intended to be called via an assignment.
// 
// Example...
//    TagMap_t *tmap=0; int cnt=0, max=0;
//
//    for (int i=0; (i<...); i++)
//    {
//       TagMap_t t = TagMap_t(...);
//       segs = TagMap_Append(tmap,cnt,max,t);
//    }
//
//-------------------------------------------------------------------------------------------------------------------

TagMap_t *TagMap_Append (TagMap_t *tmap, int & cnt, int & max, TagMap_t & tm)
{
   if ( !tmap || (cnt==max) )
   {
      max = ( (max==0) ? 1000 : max+(max/2) );

      TagMap_t *tmp = new TagMap_t[max];

      if (tmp && tmap)
         for (int i=0; (i<cnt); i++) { tmp[i]=tmap[i]; }

      if (tmap) { delete [] tmap; }

      tmap = tmp;
   }

   if (tmap && (cnt<max) ) { tmap[cnt++]=tm; }

   return(tmap);
}

//-------------------------------------------------------------------------------------------------------------------

