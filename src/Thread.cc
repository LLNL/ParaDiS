#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef PTHREADS
#include <pthread.h>
#endif

#include "Thread.h"

static void _init(Thread & t)
{
   t.thrd_scope_id     = 0;
   t.thrd_sched_id     = 0;

   t.thrd_priority     = 0;
   t.thrd_priority_min = 0;
   t.thrd_priority_max = 0;

   t.thrd_status       = THREAD_JOINED;

#ifdef PTHREADS
   t.thrd              = 0;
   pthread_attr_init(&t.thrd_attr); 

   t.thrd_priority     = SCHED_OTHER;
   t.thrd_priority_min = sched_get_priority_min(t.thrd_sched_id);
   t.thrd_priority_max = sched_get_priority_max(t.thrd_sched_id);
#endif
}

Thread::Thread(void)                            { _init(*this); }
Thread::Thread(void *(*cb)(void *), void *args) { _init(*this); Spawn(cb,args); }

Thread::~Thread()
{
   if (thrd_status==THREAD_SPAWNED) 
   { printf("Thread::~Thread() : warning - attempt to delete an active thread!\n"); }
}

#ifdef PTHREADS

void Thread::Spawn (void *(*cb)(void *), void *args)
{
   if (cb) 
   { 
      thrd_status=THREAD_SPAWNED; 
      pthread_create(&thrd, &thrd_attr, cb, (void *) args ); 
   }
}

void Thread::Join               (void) { pthread_join(thrd,NULL); 
                                         thrd_status=THREAD_JOINED; } 

void Thread::Set_State_Detached (void) { pthread_attr_setdetachstate(&thrd_attr, PTHREAD_CREATE_DETACHED); } 
void Thread::Set_State_Joinable (void) { pthread_attr_setdetachstate(&thrd_attr, PTHREAD_CREATE_JOINABLE); } 

void Thread::Set_Scope_System   (void) { thrd_scope_id=PTHREAD_SCOPE_SYSTEM ; pthread_attr_setscope(&thrd_attr,thrd_scope_id); } 
void Thread::Set_Scope_Process  (void) { thrd_scope_id=PTHREAD_SCOPE_PROCESS; pthread_attr_setscope(&thrd_attr,thrd_scope_id); } 

void Thread::Set_Sched_OTHER (void)
{ 
   thrd_sched_id=SCHED_OTHER; 
   pthread_attr_setschedpolicy ( &thrd_attr, thrd_sched_id); 
   thrd_priority_min = sched_get_priority_min(thrd_sched_id);
   thrd_priority_max = sched_get_priority_max(thrd_sched_id); 
}

void Thread::Set_Sched_FIFO (void)
{ 
   thrd_sched_id=SCHED_FIFO; 
   pthread_attr_setschedpolicy ( &thrd_attr, thrd_sched_id); 
   thrd_priority_min = sched_get_priority_min(thrd_sched_id);
   thrd_priority_max = sched_get_priority_max(thrd_sched_id); 
}

void Thread::Set_Sched_RR (void)
{
   thrd_sched_id=SCHED_RR; 
   pthread_attr_setschedpolicy ( &thrd_attr, thrd_sched_id); 
   thrd_priority_min = sched_get_priority_min(thrd_sched_id);
   thrd_priority_max = sched_get_priority_max(thrd_sched_id); 
}

void Thread::Set_Priority (const int pval)
{
   struct sched_param sp;

   int pv = pval;
       pv = ( (pv>thrd_priority_min) ? pv : thrd_priority_min );
       pv = ( (pv<thrd_priority_max) ? pv : thrd_priority_max );

   sp.sched_priority = pv; 

   int err = pthread_attr_setschedparam(&thrd_attr, &sp); 

   if (err!=0) printf("Thread_Set_Priority() - error setting scheduler priority\n");

   thrd_priority = pv;
}

#else

// stubs when pthreads not enabled...

void Thread::Spawn               (void *(*cb)(void *), void *args) {}
void Thread::Join                (void) {}

void Thread::Set_State_Detached  (void) {}
void Thread::Set_State_Joinable  (void) {}

void Thread::Set_Scope_System    (void) {}
void Thread::Set_Scope_Process   (void) {}
  
void Thread::Set_Sched_OTHER     (void) {}
void Thread::Set_Sched_FIFO      (void) {}
void Thread::Set_Sched_RR        (void) {}

void Thread::Set_Priority        (const int pval) {} 

#endif // PTHREADS

