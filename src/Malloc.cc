#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>
#endif

#include <sys/resource.h>

#include "omp_portability.h"

#include "Malloc.h"

//---------------------------------------------------------------------------------------------------------
// (note - malloc, calloc, realloc, and free cannot be overridden within this module)
//---------------------------------------------------------------------------------------------------------

#ifdef malloc
#undef malloc
#endif

#ifdef calloc
#undef calloc
#endif

#ifdef realloc
#undef realloc
#endif

#ifdef free
#undef free
#endif

//---------------------------------------------------------------------------------------------------------
// PDS_HDR_BYTES - defined for various diagnostic levels...
//---------------------------------------------------------------------------------------------------------

#if !defined(MALLOC_DIAGS)
#define PDS_HDR_BYTES 0
#endif

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS==0)
#define PDS_HDR_BYTES 0
#endif

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS==1)
#define PDS_HDR_BYTES (sizeof(size_t))
#endif

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS==2)
#define PDS_HDR_BYTES (sizeof(size_t)+sizeof(unsigned long))
#endif

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS==3)
#define PDS_HDR_BYTES (sizeof(size_t)+sizeof(unsigned long)+sizeof(char[64])+sizeof(unsigned long))
#endif

//------------------------------------------------------------------------------------------------------------
// ATTR_UNUSED is used to quiet compiler warnings...
//------------------------------------------------------------------------------------------------------------

#ifdef __GNUC__
#define  ATTR_UNUSED __attribute__ ((unused))
#else
#define  ATTR_UNUSED 
#endif

//---------------------------------------------------------------------------------------------------------
// Memory diagnostic statics....
//---------------------------------------------------------------------------------------------------------

int             PDS_Memory_Diagnostics::pds_init        = 0;  ///< diags initialized?

unsigned long   PDS_Memory_Diagnostics::pds_rcnt        = 0;  ///< current number of memory allocations
unsigned long   PDS_Memory_Diagnostics::pds_rmax        = 0;  ///< maximum number of memory allocations
unsigned long   PDS_Memory_Diagnostics::pds_alloc       = 0;  ///< current allocated heap bytes
unsigned long   PDS_Memory_Diagnostics::pds_hwm         = 0;  ///< current high water mark of memory allocations

unsigned long   PDS_Memory_Diagnostics::pds_rss_curr    = 0;  ///< current  memory allocation (not supported on all systems)
unsigned long   PDS_Memory_Diagnostics::pds_rss_peak    = 0;  ///< peak     memory allocation (not supported on all systems)
unsigned long   PDS_Memory_Diagnostics::pds_rss_phys    = 0;  ///< physical memory allocation (not supported on all systems)

int             PDS_Memory_Diagnostics::pds_pcnt        = 0;  ///< current number of PDS pointers
int             PDS_Memory_Diagnostics::pds_pmax        = 0;  ///< maximum number of PDS pointers available
unsigned char **PDS_Memory_Diagnostics::pds_ptrs        = 0;  ///< array of PDS pointers

#ifdef _OPENMP
omp_lock_t      PDS_Memory_Diagnostics::pds_omp_lock       ;  ///< mutex lock when OpenMP active
#endif

//---------------------------------------------------------------------------------------------------------

PDS_Memory_Diagnostics::PDS_Memory_Diagnostics(void)
{
   if (!pds_init)
   {
      OMP_INIT_LOCK(&pds_omp_lock); 
      pds_init=1;
   }
}

// Get_RSS_Current()
//
// Returns the current memory use.
//---------------------------------------------------------------------------------------------------------

unsigned long  PDS_Memory_Diagnostics::Get_RSS_Current (void)
{

#if defined(__APPLE__) && defined(__MACH__)
   // Apple (OSX)...
   {
      struct mach_task_basic_info info;

      mach_msg_type_number_t infoCount = { MACH_TASK_BASIC_INFO_COUNT };

      if ( task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) == KERN_SUCCESS) 
      {
         pds_rss_curr = (unsigned long) info.resident_size;
      }
   }
#endif

#if defined(__linux) || defined(__linux__) || defined(linux) || defined(__gnu_linux__)
   // Various flavors of Linux...
   {
      FILE *fd = fopen("/proc/self/statm","r");

      if (fd)
      {
         struct { unsigned long size,rss,share,text,lib,data,dt; } statm;
     
         fscanf(fd,"%ld %ld %ld %ld %ld %ld %ld",  
            &statm.size, &statm.rss, &statm.share, &statm.text, &statm.lib, &statm.data, &statm.dt );

         pds_rss_curr = statm.rss * (unsigned long) sysconf(_SC_PAGESIZE);
         pds_rss_peak = ( (pds_rss_curr>pds_rss_peak) ? pds_rss_curr : pds_rss_peak );

         fclose(fd);
      }
   }
#endif

   return(pds_rss_curr);
}

// Get_RSS_Peak()
//
// Returns the peak memory use.
//---------------------------------------------------------------------------------------------------------

unsigned long  PDS_Memory_Diagnostics::Get_RSS_Peak (void)
{
#if defined(__unix__) || defined(__unix) || defined(unix) || ( defined(__APPLE__) && defined(__MACH__) )
   // UNIX and Apple (OSX)...
   {
      struct rusage rusage;
      getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__)
      pds_rss_peak = (unsigned long) rusage.ru_maxrss;            // (osx  returns usage in bytes)
#else
      pds_rss_peak = (unsigned long) rusage.ru_maxrss * 1024L;    // (unix returns usage in kbytes)
#endif
   }
#endif

#if defined(__linux) || defined(__linux__) || defined(linux) || defined(__gnu_linux__)
   // Various flavors of Linux...
   {
      FILE *fd = fopen("/proc/self/statm","r");

      if (fd)
      {
         struct { unsigned long size,rss,share,text,lib,data,dt; } statm;
     
         fscanf(fd,"%ld %ld %ld %ld %ld %ld %ld",  
            &statm.size, &statm.rss, &statm.share, &statm.text, &statm.lib, &statm.data, &statm.dt );

         pds_rss_curr = statm.rss * (unsigned long) sysconf(_SC_PAGESIZE);
         pds_rss_peak = ( (pds_rss_curr>pds_rss_peak) ? pds_rss_curr : pds_rss_peak );

         fclose(fd);
      }
   }
#endif

   return(pds_rss_peak);
}

// Get_RSS_Physical()
//
// Returns physical memory size (currently only works on OSX)...
//---------------------------------------------------------------------------------------------------------

unsigned long  PDS_Memory_Diagnostics::Get_RSS_Physical (void)
{
#if defined(__unix__) || defined(__unix) || defined(unix) || ( defined(__APPLE__) && defined(__MACH__) )
   // Apple (OSX)...
   {
      unsigned long pages     = { (unsigned long) sysconf(_SC_PHYS_PAGES) };
      unsigned long page_size = { (unsigned long) sysconf(_SC_PAGE_SIZE ) };
      pds_rss_phys = pages * page_size;
   }
#endif

   return(pds_rss_phys);
}

// PDS_Memory_Diagnostics::Print()
//
// Simple formatted print of the current memory diagnostics.
//---------------------------------------------------------------------------------------------------------

void PDS_Memory_Diagnostics::Print (FILE *fd)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)
   if (fd)
   {
      fprintf(fd, "pds_memory_diagnostics...\n");
      fprintf(fd, "   pds_rcnt  = %lu\n"   ,                 pds_rcnt  );
      fprintf(fd, "   pds_rmax  = %lu\n"   ,                 pds_rmax  );
      fprintf(fd, "   pds_alloc = %lu\n"   , (unsigned long) pds_alloc );
      fprintf(fd, "   pds_hwm   = %lu\n"   , (unsigned long) pds_hwm   );

      fprintf(fd, "   pds_pcnt  = %d\n"    ,                 pds_pcnt  );
      fprintf(fd, "   pds_pmax  = %d\n"    ,                 pds_pmax  );
      fprintf(fd, "   pds_ptrs  = %-16p\n" , (void *)        pds_ptrs  );

#if (MALLOC_DIAGS>2)
      if (pds_ptrs)
      {
         for (int i=0; (i<pds_pcnt); i++)
         {
            if (pds_ptrs[i])
            {
               unsigned char   *p      = pds_ptrs[i];
                        size_t  pbytes = *((         size_t *) p);  p += sizeof(size_t);
               unsigned long   *pfx    =   (unsigned long   *) p ;  p += sizeof(unsigned long);
               unsigned long   *sfx    =   (unsigned long   *) (pds_ptrs[i] + pbytes - sizeof(unsigned long));
                        char   *file   =   (         char   *) p ;  p += sizeof(char[64]);
               unsigned long    ln     = *((unsigned long   *) p);  p += sizeof(unsigned long);

               fprintf(fd,"   %6d : %-16p %8lu %s(%lu)\n", i, (void *) pds_ptrs[i], pbytes, file, ln );
            }
         }
      }
#endif
   }
#endif
}

// PDS_Memory_Diagnostics::Validate()
//
// If the PDS block pointers are active, traverses the block pointers and checks the
// PDS frame markers. Returns 0 (false) upon encountering an invalid frame marker.
//---------------------------------------------------------------------------------------------------------

int PDS_Memory_Diagnostics::Validate (void)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>2)
   if (pds_ptrs)
   {
      for (int i=0; (i<pds_pcnt); i++)
      {
         unsigned char *p = pds_ptrs[i];

                  size_t  pbytes = *((         size_t *) p);  p += sizeof(size_t);
         unsigned long   *pfx    =   (unsigned long   *) p ;  p += sizeof(unsigned long);
         unsigned long   *sfx    =   (unsigned long   *) (pds_ptrs[i] + pbytes - sizeof(unsigned long));
                  char   *file   =   (         char   *) p ;  p += sizeof(char[64]);
         unsigned long    ln     = *((unsigned long   *) p);  p += sizeof(unsigned long);

         if ( !(   pfx && (*pfx==PDS_MBLK_PFX)
                && sfx && (*sfx==PDS_MBLK_SFX)) )
         {
                     size_t  pbytes = *((         size_t *) p);  p += sizeof(size_t);
            unsigned long   *pfx    =   (unsigned long   *) p ;  p += sizeof(unsigned long);
            unsigned long   *sfx    =   (unsigned long   *) (pds_ptrs[i] + pbytes - sizeof(unsigned long));

            printf("PDS memory validation error! %s(ln=%lu)\n", file, ln );
            return(0);
         }
      }
   }
#endif

   return(1);
}

//---------------------------------------------------------------------------------------------------------
// pblk_cmp()
//
// Custom comparison routine for comparing two PDS memory blocks.
//---------------------------------------------------------------------------------------------------------

ATTR_UNUSED
static int pblk_cmp (const void *p0, const void *p1)
{
   unsigned char *pds0 = ( p0 ? *((unsigned char **) p0) : 0 );
   unsigned char *pds1 = ( p1 ? *((unsigned char **) p1) : 0 );

   if (pds0 && pds1)
   {
      if ( pds0 < pds1) return(-1);
      if ( pds0 > pds1) return(+1);
   }

   return(0);
}

//---------------------------------------------------------------------------------------------------------
// Memblock_Append()
//
// Will append a memory block pointer to the memory diagnostics array.
//
// Note - this can be a huge performance killer because every time a diagnostic block is
// added to the diagnostics array - the array is sorted for later lookups and deletions.
//---------------------------------------------------------------------------------------------------------

void PDS_Memory_Diagnostics::Memblock_Append(const unsigned char *pblk)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>2)
   OMP_LOCK(&pds_omp_lock);

   if ( !pds_ptrs || (pds_pcnt==pds_pmax) )
   {
      pds_pmax = ( (pds_pmax==0) ? 1000 : (2*pds_pmax) );

      unsigned char **tmp = (unsigned char **) malloc(pds_pmax*sizeof(unsigned char *));

      if (tmp && pds_ptrs)
         for (int i=0; (i<pds_pcnt); i++) { tmp[i]=pds_ptrs[i]; }

      if (pds_ptrs) { free(pds_ptrs); }

      pds_ptrs = tmp;
   }

   if (pds_ptrs && (pds_pcnt<pds_pmax) )
   { pds_ptrs[pds_pcnt++]=(unsigned char *) pblk; }

   // Sort the addresses...

   if (pds_ptrs && (pds_pcnt>0))
   { qsort(pds_ptrs, pds_pcnt, sizeof(unsigned char *), pblk_cmp); }

   OMP_UNLOCK(&pds_omp_lock);
#endif
}

//---------------------------------------------------------------------------------------------------------
// Memblock_Index()
//
// Performs a binary search for a given memory block. Returns the index of the memory block
// entry or -1 (if not found).  Assumes the PDS memory block pointers are sorted above.
//---------------------------------------------------------------------------------------------------------

int PDS_Memory_Diagnostics::Memblock_Index (const unsigned char *pblk)
{
   if (!pblk    ) { printf("pds memory diags - null memory pointer\n");              return(-1); }
   if (!pds_ptrs) { printf("pds memory diags - cannot search empty pointer list\n"); return(-1); }

   if (pds_ptrs && (pds_pcnt>0) && pblk)
   {
      int i0=0;
      int i1=(pds_pcnt-1);

      while (i0<=i1)
      {
         int i = (i0+i1)/2;

         unsigned char *p = pds_ptrs[i];

         if ( p == pblk ) { return(i); }
         if ( p <  pblk ) { i0=i+1;    }
         if ( p >  pblk ) { i1=i-1;    }
      }

      printf("warning - pds memory error - cannot locate memory block\n");
   }

   return(-1);
}

//---------------------------------------------------------------------------------------------------------
// Memblock_Delete()
//
// Deletes a memory block from the PDS memory block pointer array.
//---------------------------------------------------------------------------------------------------------

void PDS_Memory_Diagnostics::Memblock_Delete(const unsigned char *pblk)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>2)
   OMP_LOCK(&pds_omp_lock);

   if (!pblk    ) { printf("warning - pds memory diags - null memory pointer\n");                   }
   if (!pds_ptrs) { printf("warning - pds memory diags - cannot delete from empty pointer list\n"); }

   if (pds_ptrs)
   {
      int pndx = Memblock_Index(pblk);

      if (pndx>=0)
      {
         for (int i=pndx; (i<(pds_pcnt-1)); i++)
            pds_ptrs[i]=pds_ptrs[i+1];

         pds_ptrs[pds_pcnt-1] = 0;
         pds_pcnt--;
      }
   }

   OMP_UNLOCK(&pds_omp_lock);
#endif
}

//---------------------------------------------------------------------------------------------------------

void PDS_Memory_Diagnostics::Allocation_Push (const size_t bytes)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)
   OMP_LOCK(&pds_omp_lock);

   pds_rcnt++;
   pds_rmax   = ( (pds_rcnt>pds_rmax) ? pds_rcnt : pds_rmax );
   pds_alloc += bytes;
   pds_hwm    = ( (pds_alloc>pds_hwm) ? pds_alloc : pds_hwm );

   OMP_UNLOCK(&pds_omp_lock);
#endif
}

void PDS_Memory_Diagnostics::Allocation_Pop (const size_t bytes)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)
   OMP_LOCK(&pds_omp_lock);

   if (pds_rcnt>0)       pds_rcnt--;
   if (pds_alloc>=bytes) pds_alloc -= bytes;

   OMP_UNLOCK(&pds_omp_lock);
#endif
}

//---------------------------------------------------------------------------------------------------------
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)
static inline size_t align8 (const size_t bytes) { return( bytes + ( (bytes%8) ? (8-(bytes%8)) : 0 ) ); }
#endif
//---------------------------------------------------------------------------------------------------------

// PDS_Bytes()
//
// Given a requested buffer size - returns the total byte count needed to store the buffer
// and any given diagnostic info currently included in the build.
//---------------------------------------------------------------------------------------------------------

size_t PDS_Bytes (size_t bytes)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS==1)
   bytes = ( (bytes>0) ? align8(bytes)+sizeof(size_t) : 0 );
#endif

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS==2)
   bytes = ( (bytes>0) ? align8(bytes)+sizeof(size_t)+2*sizeof(unsigned long) : 0 );
#endif

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS==3)
   bytes = ( (bytes>0) ? align8(bytes)+sizeof(size_t)+2*sizeof(unsigned long)+sizeof(char[64])+sizeof(unsigned long) : 0 );
#endif

   return(bytes);
}

// PDS_Bytes()
//
// Given a pointer to an existing buffer, returns the allocated size of the buffer using
// the size embedded in the diagnostics.
//---------------------------------------------------------------------------------------------------------

size_t PDS_Bytes (void *mbuf)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)
#if (MALLOC_DIAGS==1)
   return ( mbuf ? *((size_t *)(((unsigned char *) mbuf) - sizeof(size_t)) ) : 0 );
#endif

#if (MALLOC_DIAGS==2)
   return ( mbuf ? *((size_t *)(((unsigned char *) mbuf) - sizeof(size_t) - sizeof(unsigned long)) ) : 0 );
#endif

#if (MALLOC_DIAGS==3)
   return ( mbuf ? *((size_t *)(((unsigned char *) mbuf) - sizeof(size_t) - sizeof(unsigned long) - sizeof(char[64]) - sizeof(unsigned long)) ) : 0 );
#endif

#else
   return(0);
#endif
}

// PDS_Prefix()
//
// Returns a pointer to the PDS frame marker prefix.
//---------------------------------------------------------------------------------------------------------

unsigned long *PDS_Prefix (void *mbuf)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>1)
   if (mbuf)
   {
      return( (unsigned long *) (((unsigned char *) mbuf) - PDS_HDR_BYTES + sizeof(size_t)) );
   }
#endif

   return(0);
}

// PDS_Suffix()
//
// Returns a pointer to the PDS frame marker suffix.
//---------------------------------------------------------------------------------------------------------

unsigned long *PDS_Suffix (void *mbuf)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>1)
   if (mbuf)
   {
      unsigned char *p = ((unsigned char *) mbuf) - PDS_HDR_BYTES;

      size_t bytes = *((size_t *) p );

      p += bytes;
      p -= sizeof(unsigned long);

      return( (unsigned long *) p );
   }
#endif

   return(0);
}

// PDS_Base()
//
// Returns a pointer to the base address of a block of memory allocated via the PDS memory
// management routines.
//---------------------------------------------------------------------------------------------------------

unsigned char *PDS_Base (void *mbuf)
{
   return( mbuf ? (((unsigned char *) mbuf) - PDS_HDR_BYTES) : 0 );
}

// PDS_Fname()
//
// Returns a pointer to the file name that allocated the PDS memory block
//---------------------------------------------------------------------------------------------------------

char *PDS_Fname (void *mbuf)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>2)
   return( mbuf ? (char *) (((unsigned char *) mbuf) - sizeof(unsigned long) - sizeof(char[64])) : 0 );
#else
   return(0);
#endif
}

// PDS_Fline()
//
// Returns the line number of the file that allocated a PDS memory block.
//---------------------------------------------------------------------------------------------------------

int PDS_Fline (void *mbuf)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>2)
   return( mbuf ? (int) (*((unsigned long *) (((unsigned char *) mbuf)-sizeof(unsigned long)))) : 0 );
#else
   return(0);
#endif
}

ATTR_UNUSED
static void PDS_Diags_Set
(
          unsigned char   *mbuf ,   ///< points to the base of the allocation
                   size_t  bytes,   ///< the size of the allocation
   const           char   *file ,   ///< points to the source file name        (e.g. __FILE__)
   const           int     ln       ///< points to the source file line number (e.g. __LINE__)
)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)
   if (mbuf)
   {
      unsigned char *p = mbuf;

      *((size_t *) p) = bytes;  p += sizeof(size_t);

#if (MALLOC_DIAGS>1)
      unsigned long *pfx = (unsigned long *) p;   p += sizeof(unsigned long);
      unsigned long *sfx = (unsigned long *) (mbuf + bytes - sizeof(unsigned long));

      *pfx = PDS_MBLK_PFX;
      *sfx = PDS_MBLK_SFX;
#endif

#if (MALLOC_DIAGS>2)
      strncpy((char *) p, file, sizeof(char[64]));  p += sizeof(char[64]);
      *((unsigned long *) p) = ln;    p += sizeof(unsigned long);
#endif
   }
#endif
}

// PDS_Malloc()
//
// Allocates space for the request including any necessary space to hold the block
// diagnostics and frame info. The returned pointer is adjusted to point to the
// first available, user addressable byte of memory.
//
// Note - any memory allocated via PDS_Malloc() must be freed via PDS_Free().
//---------------------------------------------------------------------------------------------------------

void *PDS_Malloc(size_t bytes, const char *file, const int ln)
{
   unsigned char *mbuf=0;

   if (bytes>0)
   {
      bytes = PDS_Bytes(bytes);                  // (adjust request to hold space for the diagnostics)
      mbuf  = (unsigned char *) malloc(bytes);   // (allocate via malloc)

      if (!mbuf) { printf("%s(%d): pds memory allocation fail (req=%lu bytes)\n", file, ln, bytes); }

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)
      PDS_Diags_Set(mbuf,bytes,file,ln);

      PDS_Memory_Diagnostics::Allocation_Push(bytes);
      PDS_Memory_Diagnostics::Memblock_Append(mbuf );

      mbuf += PDS_HDR_BYTES;
#endif
   }

   return((void *) mbuf);
}

void *PDS_Malloc(size_t cnt, size_t bytes, const char *file, const int ln)
{
   return( PDS_Malloc( (cnt*bytes), file, ln ) );
}

//---------------------------------------------------------------------------------------------------------

void *PDS_Calloc(size_t bytes, const char *file, const int ln )
{
   unsigned char *mbuf=0;

   if (bytes>0)
   {
      bytes = PDS_Bytes(bytes);                    // (adjust request to hold space for the diagnostics)
      mbuf  = (unsigned char *) calloc(1,bytes);   // (allocate via calloc)

      if (!mbuf) { printf("%s(%d): pds memory allocation fail (req=%lu bytes)\n", file, ln, bytes); }

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)
      PDS_Diags_Set(mbuf,bytes,file,ln);

      PDS_Memory_Diagnostics::Allocation_Push(bytes);
      PDS_Memory_Diagnostics::Memblock_Append(mbuf );

      mbuf += PDS_HDR_BYTES;
#endif
   }

   return((void *) mbuf);
}

void *PDS_Calloc(size_t cnt, size_t bytes, const char *file, const int ln)
{
   return( PDS_Calloc((cnt*bytes),file,ln) );
}

//---------------------------------------------------------------------------------------------------------

void *PDS_Realloc(void *mbuf, size_t bytes, const char *file, const int ln)
{
   if (!mbuf)               { return( PDS_Malloc(bytes,file,ln) ); }
   if ( mbuf && (bytes==0)) {         PDS_Free(mbuf); return(0); }

#if !defined(MALLOC_DIAGS) || ( defined(MALLOC_DIAGS) && (MALLOC_DIAGS==0) )
   mbuf = realloc(mbuf,bytes);
#endif

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)
   size_t pbytes = PDS_Bytes(mbuf );  // pbytes = old buffer size (including diagnostics)
   size_t nbytes = PDS_Bytes(bytes);  // nbytes = new buffer size (including diagnostics)

   if (mbuf && (nbytes>0) && (nbytes>pbytes))
   {
      unsigned char *pbase = PDS_Base(mbuf);
      unsigned char *nbase = (unsigned char *) realloc(pbase,nbytes);

      if (!nbase) { printf("%s(%d): pds memory reallocation fail (req=%lu bytes)\n", file, ln, nbytes); }

      PDS_Diags_Set(nbase,nbytes,file,ln);

      PDS_Memory_Diagnostics::Allocation_Push(nbytes);
      PDS_Memory_Diagnostics::Allocation_Pop (pbytes);

      if (nbase!=pbase)
      {
         PDS_Memory_Diagnostics::Memblock_Append(nbase );
         PDS_Memory_Diagnostics::Memblock_Delete(pbase);
      }

      mbuf = (void *) (nbase + PDS_HDR_BYTES);
   }
#endif

   if (!mbuf) { printf("%s(%d): pds memory reallocation fail (req=%lu bytes)\n", file, ln, bytes); }

   return(mbuf);
}

void *PDS_Realloc(void *mbuf, size_t cnt, size_t bytes, const char *file, const int ln)
{
   return( PDS_Realloc(mbuf,(cnt*bytes),file,ln) );
}

//---------------------------------------------------------------------------------------------------------

int PDS_Validate (void *mbuf)
{
#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>1)
   if (mbuf)
   {
      unsigned char *p = PDS_Base(mbuf);

      size_t bytes = *((size_t *) p);

      unsigned long *pfx = (unsigned long *) (p+      sizeof(size_t)       );
      unsigned long *sfx = (unsigned long *) (p+bytes-sizeof(unsigned long));

      return(    pfx && (*pfx==PDS_MBLK_PFX)
              && sfx && (*sfx==PDS_MBLK_SFX) ? 1 : 0 );
   }
#endif

   return(1);
}

//---------------------------------------------------------------------------------------------------------

void PDS_Free(void *mbuf)
{
   if (mbuf)
   {
      unsigned char *p = PDS_Base(mbuf);

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)
      size_t bytes = *((size_t *) p);
#endif

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>1)
      unsigned long *pfx = (unsigned long *) (p+      sizeof(size_t)       );
      unsigned long *sfx = (unsigned long *) (p+bytes-sizeof(unsigned long));

      if ( *pfx != PDS_MBLK_PFX )  { printf("warning - pds memory prefix checksum error\n"); }
      if ( *sfx != PDS_MBLK_SFX )  { printf("warning - pds memory suffix checksum error\n"); }
#endif

#if defined(MALLOC_DIAGS) && (MALLOC_DIAGS>0)
      PDS_Memory_Diagnostics::Allocation_Pop(bytes);
      PDS_Memory_Diagnostics::Memblock_Delete(p);
#endif

      free(p);
   }
}

//---------------------------------------------------------------------------------------------------------

int PDS_Memcheck (void) { return( PDS_Memory_Diagnostics::Validate() ); }

//---------------------------------------------------------------------------------------------------------
