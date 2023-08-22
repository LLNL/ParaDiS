#pragma once

#ifndef _PDS_THREAD_H
#define _PDS_THREAD_H

//---------------------------------------------------------------------------------------------------------

#if defined(PTHREADS)
#include <pthread.h>
#endif

typedef enum 
{
   THREAD_SPAWNED,  ///< thread is currently active and running
   THREAD_JOINED    ///< thread has been joined and is idle
} Thread_Status;

class Thread
{
   public :
#ifdef PTHREADS
      pthread_t       thrd             ;   ///< thread ID
      pthread_attr_t  thrd_attr        ;   ///< attributes for this thread
#endif

      int             thrd_scope_id    ;   ///< thread scope ID
      int             thrd_sched_id    ;   ///< thread sceduler ID

      int             thrd_priority    ;   ///< the current thread priority
      int             thrd_priority_min;   ///< sets the minimum thread priority
      int             thrd_priority_max;   ///< sets the maximum thread priority

      Thread_Status   thrd_status      ;   ///< current thread run status (running or idle)

   public :
      Thread(void);
      Thread(void *(*cb)(void *), void *args=0);
     ~Thread();

      void    Spawn               (void *(*cb)(void *), void *args=0);

      void    Join                (void);

      void    Set_State_Detached  (void);
      void    Set_State_Joinable  (void);

      void    Set_Scope_System    (void);
      void    Set_Scope_Process   (void);
  
      void    Set_Sched_OTHER     (void);
      void    Set_Sched_FIFO      (void);
      void    Set_Sched_RR        (void);

      void    Set_Priority        (const int pval); 

      int     Spawned             (void) { return( (thrd_status==THREAD_SPAWNED) ? 1 : 0 ); }
      int     Joined              (void) { return( (thrd_status==THREAD_JOINED ) ? 1 : 0 ); }
};

//---------------------------------------------------------------------------------------------------------

#ifdef PTHREADS

#define PTHRD_LOCK_INIT(a)     { if (a) { pthread_mutex_init   ((a),NULL); } }
#define PTHRD_LOCK_DESTROY(a)  { if (a) { pthread_mutex_destroy((a))     ; } }
#define PTHRD_TRY_LOCK(a)      { if (a) { pthread_mutex_trylock((a))     ; } }
#define PTHRD_LOCK(a)          { if (a) { pthread_mutex_lock   ((a))     ; } }
#define PTHRD_UNLOCK(a)        { if (a) { pthread_mutex_unlock ((a))     ; } }

#else  
// empty macros when pthreads not active...

#define PTHRD_LOCK_INIT(a)
#define PTHRD_LOCK_DESTROY(a)
#define PTHRD_TRY_LOCK(a)
#define PTHRD_LOCK(a)
#define PTHRD_UNLOCK(a)

#endif  // PTHREADS

#endif  // _PDS_THREAD_H
