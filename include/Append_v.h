#pragma once

#ifndef __PDS_APPEND_V_H
#define __PDS_APPEND_V_H

#define APPEND_CDECL(vtype) \
extern vtype *Append (vtype *v, int & vcnt, int & vmax, vtype val, int tmin=32, int t0=1, int t1=1);

APPEND_CDECL(         char )
APPEND_CDECL(unsigned char )
APPEND_CDECL(         short)
APPEND_CDECL(unsigned short)
APPEND_CDECL(         int  )
APPEND_CDECL(unsigned int  )
APPEND_CDECL(         long )
APPEND_CDECL(unsigned long )
APPEND_CDECL(         long long   )
APPEND_CDECL(unsigned long long   )
APPEND_CDECL(         float       )
APPEND_CDECL(         double      )
APPEND_CDECL(         long double )
APPEND_CDECL(         void *      )

#undef APPEND_CDECL

#define APPEND_CDECL(vtype) \
extern vtype *Append (vtype *v, int & vcnt, int & vmax, vtype *vsrc, int vn, int tmin=32, int t0=1, int t1=1);

APPEND_CDECL(         char )
APPEND_CDECL(unsigned char )
APPEND_CDECL(         short)
APPEND_CDECL(unsigned short)
APPEND_CDECL(         int  )
APPEND_CDECL(unsigned int  )
APPEND_CDECL(         long )
APPEND_CDECL(unsigned long )
APPEND_CDECL(         long long   )
APPEND_CDECL(unsigned long long   )
APPEND_CDECL(         float       )
APPEND_CDECL(         double      )
APPEND_CDECL(         long double )
APPEND_CDECL(         void *      )

#undef APPEND_CDECL

extern void *Append (void *v, int & vcnt, int & vmax, void *val, size_t vsize, int tmin=32, int t0=1, int t1=1);

#endif  // __PDS_APPEND_V_H
