#pragma once

#ifndef _PDS_MALLOC_H
#define _PDS_MALLOC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/resource.h>

#include "omp_portability.h"

//------------------------------------------------------------------------------------------------------------
// ParaDiS Memory Diagnostic Levels...
//
// Historically, the <stdlib.h> implementations of the dynamic memory allocation
// routines (e.g. malloc(), calloc(), free(), etc) - when space was allocated,
// an additional 4 or 8-bytes was added to the requested byte count and added to
// the start of an allocation. That byte-count prefix was used to identify the
// actual size of the allocated space. It was also used by the free() routine
// to determine how much memory was returned to the heap when an allocation
// was given back to the heap manager.
//
// While not part of the standard, many applications used that extra info to
// identify when memory leaks occurred.  ParaDiS took advantage of it as well.
//
// Unfortunately, more recent distributions of Linux have abandoned that support
// and instead, incorporate that meta data into the paging information retained
// within the heap manager.
//
// Within this module, I add additional information to all the dynamic allocations
// to assist with dealing with memory problems within ParaDiS.
//
// If the memory diagnostics are enabled, the conventional memory allocation
// routines (e.g. malloc(), calloc(), free(), etc) are replaced using a C/C++
// preprocessor macro that redirects all those allocations to replacement
// routines within this module...
//
// e.g.
//    malloc()  == PDS_Malloc()
//    calloc()  == PDS_Calloc()
//    realloc() == PDS_Realloc()
//    free()    == PDS_Free()
//
// Depending on the diagnostic level set during the build, the replacement routines
// will allocate additional space to retain things such as the size of the allocation,
// frame markers, and the file name and line number responsible for the allocation.
// The replacement routines will also force 64-bit alignment and total allocation
// tracking.
//
// Several levels are supported (with increasing impact on performance)...
//    0 : memory diagnostic support disabled (default)
//    1 : adds allocation sizes to each allocation (low impact)
//           - also adds 64-bit (8-bit) pad bytes for 64-bit alignment
//           - adds diagnostic allocation counts and high-water mark output
//    2 : adds frame markers before and after the user-addressable space
//           - also adds frame marker validation checks on delete
//           - can be used to detect out-of-bounds writes
//    3 : full diagnostic support, including block diagnostics
//           - adds block tracking (allocations, deallocations)
//           - significant performance impact
//
// Depending on the diagnostic level set, meta data is added to each
// allocation as follows...
//
//  +-------------------------------------+  <<<---- start of allocation
//  | allocation size (in bytes)          |  (level 1)
//  +-------------------------------------+
//  | frame marker (prefix)               |  (level 2)
//  +-------------------------------------+
//  | file and line number information    |  (level 3)
//  |    (e.g. __FILE__)                  |  (level 3)
//  |    (e.g. __LINE__)                  |  (level 3)
//  +-------------------------------------+  <<<---- start of user-addressable space
//  | user addressable space              |  (all levels)
//  |                                     |
//  |                                     |
//  |                                     |
//  +-------------------------------------+
//  | 64-bit alignment pad                |  (level 1)
//  +-------------------------------------+
//  | frame marker (suffix)               |  (level 2)
//  +-------------------------------------+
//
// Also note that if the diagnostic level is 3 - a separate array of allocation pointers
// is maintained that can be traversed independently to check for frame-marker checksum errors.
//
// Lastly - whatever diagnostic level is set - the allocation and free routines always assume
// to return pointers to the start of the user-addressable space. The diagnostics can be
// enabled and disabled with impunity and the upstream users will never know the difference.
//------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------
// MALLOC_DIAGS - you can use this #define to override the setting in the makefile.setup.
//------------------------------------------------------------------------------------------------------------

//#define MALLOC_DIAGS 3

#if defined(MALLOC_DIAGS) && ( (MALLOC_DIAGS<0) || (MALLOC_DIAGS>3) )
#error "MALLOC_DIAGS level not supported"
#endif

//------------------------------------------------------------------------------------------------------------
// The following macros will be used to intercept all the calls to the <stdlib> routines with their
// corresponding ParaDiS routines implemented within this module. Note that each of these macros will
// also insert C-compiler __FILE__ and __LINE__ macros to enabled better diagnostic/error output to the
// console.
//------------------------------------------------------------------------------------------------------------

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)

#define malloc(a)         PDS_Malloc (a,  __FILE__,__LINE__)
#define calloc(a,b)       PDS_Calloc (a,b,__FILE__,__LINE__)
#define realloc(a,b)      PDS_Realloc(a,b,__FILE__,__LINE__)
#define free(a)           PDS_Free   (a); (a)=0;

#endif

#define PDS_MALLOC(a)     PDS_Malloc (a,  __FILE__,__LINE__)
#define PDS_CALLOC(a,b)   PDS_Calloc (a,b,__FILE__,__LINE__)
#define PDS_REALLOC(a,b)  PDS_Realloc(a,b,__FILE__,__LINE__)
#define PDS_FREE(a)       PDS_Free   (a); (a)=0;

//------------------------------------------------------------------------------------------------------------
// If frame markers are enabled, each allocation includes a 64-bit (8-byte) prefix and suffix pattern
// that is used to help determine if an out-of-bound write has occurred.  They are checked when the
// memory is free'd.  They can also be manually checked using the PDS_Validate() routine within this
// module.
//------------------------------------------------------------------------------------------------------------

#define PDS_MBLK_PFX 0xf0f0f0f0f0f0f0f0ul
#define PDS_MBLK_SFX 0xf1f1f1f1f1f1f1f1ul

//------------------------------------------------------------------------------------------------------------

extern          size_t  PDS_Bytes    (size_t bytes);
extern          size_t  PDS_Bytes    (void  *mbuf );
extern unsigned long   *PDS_Prefix   (void  *mbuf );
extern unsigned long   *PDS_Suffix   (void  *mbuf );
extern unsigned char   *PDS_Base     (void  *mbuf );
extern          char   *PDS_Fname    (void  *mbuf );
extern          int     PDS_Fline    (void  *mbuf );

extern          void   *PDS_Malloc   (            size_t bytes, const char *file, const int ln);
extern          void   *PDS_Malloc   (size_t cnt, size_t bytes, const char *file, const int ln);

extern          void   *PDS_Calloc   (            size_t bytes, const char *file, const int ln);
extern          void   *PDS_Calloc   (size_t cnt, size_t bytes, const char *file, const int ln);

extern          void   *PDS_Realloc  (void *mbuf,             size_t bytes, const char *file, const int ln);
extern          void   *PDS_Realloc  (void *mbuf, size_t cnt, size_t bytes, const char *file, const int ln);

extern          void    PDS_Free     (void *mbuf);

extern          int     PDS_Validate (void *mbuf);

//------------------------------------------------------------------------------------------------------------

class PDS_Memory_Diagnostics
{
   public :
      static          int    pds_init       ;  ///< diags initialized?
      static unsigned long   pds_rcnt       ;  ///< current number of memory allocations
      static unsigned long   pds_rmax       ;  ///< maximum number of memory allocations
      static unsigned long   pds_alloc      ;  ///< current allocated heap bytes
      static unsigned long   pds_hwm        ;  ///< high water mark of total memory allocations

      static unsigned long   pds_rss_curr   ;  ///< current  memory allocation (not supported on all systems)
      static unsigned long   pds_rss_peak   ;  ///< peak     memory allocation (not supported on all systems)
      static unsigned long   pds_rss_phys   ;  ///< physical memory allocation (not supported on all systems)

      static          int    pds_pcnt       ;  ///< current number of PDS allocation pointers
      static          int    pds_pmax       ;  ///< maximum number of PDS allocation pointers (allocated)
      static unsigned char **pds_ptrs       ;  ///< array of PDS allocation pointers

#ifdef _OPENMP
      static omp_lock_t      pds_omp_lock   ;  ///< mutex lock when OpenMP active
#endif

   public :
      PDS_Memory_Diagnostics(void);
     ~PDS_Memory_Diagnostics() {}

      static unsigned long  Get_RSS_Current  (void);
      static unsigned long  Get_RSS_Peak     (void);
      static unsigned long  Get_RSS_Physical (void);

      static void Print (FILE *fd=stdout);

      static void Memblock_Append (const unsigned char *mbuf);
      static void Memblocks_Sort  (void);
      static int  Memblock_Index  (const unsigned char *mbuf);
      static void Memblock_Delete (const unsigned char *mbuf);

      static void Allocation_Push (const size_t bytes);
      static void Allocation_Pop  (const size_t bytes);

      static int  Validate        (void);
};

#endif
