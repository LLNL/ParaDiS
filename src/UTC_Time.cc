//------------------------------------------------------------------------------------
//
//   Module:        UTC_Time.cc
//   Description:   This module implements several UTC time utilities
//
//   Included functions:
//      Current_UTC_Time
//      UTC_Seconds
//      UTC_Microsecs
//      Elapsed_UTC_Microsecs
//      Elapsed_UTC_Seconds
//      UTC_To_Calendar
//      Calendar_To_UTC
//------------------------------------------------------------------------------------

#include <time.h>
#include <stdlib.h>
#include <sys/time.h>

#include "UTC_Time.h"

// Current_UTC_Time()
//
// Sets the current instance to the current UTC time.
//------------------------------------------------------------------------------------

long long UTC_Time::Current_UTC_Time (void)
{
   timeval tval;  gettimeofday(&tval,0);

   utc_usecs = (((long long) tval.tv_sec) * 1000000LL) + tval.tv_usec;

   return(utc_usecs);
}

// Current_UTC_Time()
//
// Sets the current instance to the current UTC time. 
// Returns the UTC time converted to conventional calendar units.
//------------------------------------------------------------------------------------

void UTC_Time::Current_UTC_Time (int &yy, int &yday, int &mo, int &dd, int &hh, int &mm, double &ss)
{
   timeval tval;  gettimeofday(&tval,0);

   utc_usecs = (((long long) tval.tv_sec) * 1000000LL) + tval.tv_usec;

   UTC_To_Calendar (yy,yday,mo,dd,hh,mm,ss);
}

// UTC_To_HHMMSS()
//
// Converts the current UTC clock time to elapsed hours, minutes, seconds.
//------------------------------------------------------------------------------------

void UTC_Time::UTC_To_HHMMSS 
(
   int    & hh,  ///< hours           (returned)
   int    & mm,  ///< minutes (0-59)  (returned)
   double & ss   ///< seconds (0-59)  (returned)
) const
{
   long long tmp = utc_usecs;

   hh = (int   ) (tmp/(3600*1000000LL));  tmp = tmp%(3600*1000000LL);
   mm = (int   ) (tmp/(  60*1000000LL));  tmp = tmp%(  60*1000000LL);
   ss = (double)  tmp/1000000.0;
}

// UTC_To_Calendar()
//
// Converts the current UTC clock time to conventional calendar time.
//------------------------------------------------------------------------------------

void UTC_Time::UTC_To_Calendar 
(
   int    &yy  ,   ///< years                (returned)
   int    &yday,   ///< year day (1-366)     (returned)
   int    &mo  ,   ///< months   (1-12)      (returned)
   int    &dd  ,   ///< days     (1-31)      (returned)
   int    &hh  ,   ///< hours    (0-23)      (returned)
   int    &mm  ,   ///< minutes  (0-59)      (returned)
   double &ss      ///< seconds  (0-59);     (returned)
) const
{
   time_t tt = (time_t) (utc_usecs/1000000LL);

   struct tm tms;

   gmtime_r(&tt,&tms);

   yy   =          tms.tm_year + 1900;
   yday =          tms.tm_yday;
   mo   =          tms.tm_mon  + 1;   
   dd   =          tms.tm_mday;
   hh   =          tms.tm_hour;
   mm   =          tms.tm_min ;
   ss   = (double) tms.tm_sec ;
   ss  += (utc_usecs%1000000LL)/1000000.0;
}

// Calendar_To_UTC()
//
// Converts a conventional calendar time to UTC time (microseconds)
//------------------------------------------------------------------------------------

void UTC_Time::Calendar_To_UTC 
(
   const int    yy,    ///< year
   const int    mo,    ///< month   (1-12)
   const int    dd,    ///< day     (1-31)
   const int    hh,    ///< hour    (00-23)
   const int    mm,    ///< minutes (00-59)
   const double ss)    ///< seconds (00-59)
{
   struct tm tms;

   tms.tm_year   = yy - 1900;
   tms.tm_mon    = mo - 1;     // (convert to zero indexed)
   tms.tm_mday   = dd;
   tms.tm_hour   = hh;
   tms.tm_min    = mm;
   tms.tm_sec    = (int) ss;

   utc_usecs  = (long long)(timegm(&tms))*1000000LL;
   utc_usecs += (long long)(ss*1000000.0)%1000000LL;
}

//------------------------------------------------------------------------------------
// Function:      Current_UTC_Time
// Description:   returns a 64-bit UTC time value in microseconds
//
// Arguments:     (none)
//------------------------------------------------------------------------------------

long long  Current_UTC_Time (void)
{
   timeval tval;  gettimeofday (&tval,0);

   return( (((long long) tval.tv_sec) * 1000000LL) + tval.tv_usec );
}

// UTC_To_HHMMSS
//
// converts a 64-bit UTC time value in microseconds to hours, minutes, seconds
//------------------------------------------------------------------------------------

void UTC_To_HHMMSS(const long long utc, int & hh, int & mm, double & ss)
{
   long long tmp = utc;

   hh = (int) (tmp/(3600*1000000LL));  tmp = tmp%(3600*1000000LL);
   mm = (int) (tmp/(  60*1000000LL));  tmp = tmp%(  60*1000000LL);
   ss = (double) tmp / 1000000.0;
}

