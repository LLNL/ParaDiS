#pragma once

#ifndef _PDS_TAGMAP_H
#define _PDS_TAGMAP_H

#include "Tag.h"

class TagMap_t
{
   public : 
      Tag_t  oldTag;  ///< old tag
      Tag_t  newTag;  ///< new/renamed tag

   public :
      TagMap_t(void);
      TagMap_t(const TagMap_t & t);
   
     ~TagMap_t();

      const TagMap_t & operator =  (const TagMap_t & t);

      int operator == (const TagMap_t & t) const;
      int operator != (const TagMap_t & t) const;

};

inline TagMap_t::TagMap_t (void)                 { oldTag=Tag_t();  newTag=Tag_t();  }
inline TagMap_t::TagMap_t (const TagMap_t & t)   { oldTag=t.oldTag; newTag=t.newTag; }
   
inline TagMap_t::~TagMap_t () {}  // (nothing allocated by class)

inline const TagMap_t & TagMap_t::operator =  (const TagMap_t & t) { oldTag=t.oldTag; newTag=t.newTag; return(*this); }

inline int TagMap_t::operator == (const TagMap_t & t) const { return ( (oldTag==t.oldTag) && ( newTag==t.newTag) ? 1 : 0 ); }
inline int TagMap_t::operator != (const TagMap_t & t) const { return ( (oldTag==t.oldTag) && ( newTag==t.newTag) ? 0 : 1 ); }

extern TagMap_t *TagMap_Append (TagMap_t *tmap, int & cnt, int & max, TagMap_t & tm);

#endif
