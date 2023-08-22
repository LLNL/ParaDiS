#pragma once

#ifndef _PDS_TIMER_H
#define _PDS_TIMER_H

/*************************************************************************
 *
 *  Timer.h - define the timing structures
 *
 ************************************************************************/

#include "Typedefs.h"

class Timer_t 
{
   public :
      real8    startTime;  ///< time at which most recent TimerStart called
      real8    incr     ;  ///< time of most recent event for this event type
      real8    accum    ;  ///< accumulated time for this event type
      real8    save     ;  ///< to save full force update times till next WriteStep
      int      started  ;  ///< event active (1=yes, 0=no)
      char    *name     ;  ///< label used during TimerPrint

   public :
      Timer_t(void)
      {
         startTime = 0.0;
         incr      = 0.0;
         accum     = 0.0;
         save      = 0.0;
         started   = 0;
         name      = 0;
      }

     ~Timer_t() {}

};


enum {
    TOTAL_TIME=0           ,
    INITIALIZE             ,
    SORT_NATIVE_NODES      ,
    COMM_SEND_GHOSTS       ,
    GHOST_COMM_BARRIER     ,
    CELL_CHARGE            ,
    CELL_CHARGE_BARRIER    ,
    CALC_FORCE             ,
    LOCAL_FORCE            ,
    REMOTE_FORCE           ,
    CALC_FORCE_BARRIER     ,
    CALC_VELOCITY          ,
    CALC_VELOCITY_BARRIER  ,
    COMM_SEND_VELOCITY     ,
    SPLIT_MULTI_NODES      ,
    COLLISION_HANDLING     ,
    POST_COLLISION_BARRIER ,
    COL_SEND_REMESH        ,
    COL_FIX_REMESH         ,
    COL_FORCE_UPDATE       ,
    GENERATE_IO            ,
    PLOT                   ,
    IO_BARRIER             ,
    REMESH_START_BARRIER   ,
    REMESH                 ,
    SEND_REMESH            ,
    FIX_REMESH             ,
    FORCE_UPDATE_REMESH    ,
    REMESH_END_BARRIER     ,
    MIGRATION              ,
    MIGRATION_BARRIER      ,
    LOADCURVE              ,
    LOAD_BALANCE           ,
    SEGFORCE_COMM          ,
    TIMER_BLOCK_SIZE         ///< (must be last in the list)
};

// timer function prototypes....

extern void TimeAtRestart     (Home_t *home, int stage);
extern void TimerClear        (Home_t *home, int index);
extern void TimerClearAll     (Home_t *home);
extern void TimerInit         (Home_t *home);
extern void TimerInitDLBReset (Home_t *home);
extern void TimerPrint        (Home_t *home);
extern void TimerReinitialize (Home_t *home);
extern void TimerStart        (Home_t *home, int index);
extern void TimerStop         (Home_t *home, int index);
extern void TimerSave         (Home_t *home, int index);

#endif  // _PDS_TIMER_H
