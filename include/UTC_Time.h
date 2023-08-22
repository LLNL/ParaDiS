#pragma once

#ifndef _PDS_UTC_TIME_H
#define _PDS_UTC_TIME_H

//------------------------------------------------------------------------------------
//  UTC_Time.h - Coordinated Universal Time (UTC) support
//------------------------------------------------------------------------------------

class UTC_Time 
{
   public :
      long long   utc_usecs;     // UTC time in microseconds

   public :
      UTC_Time(void)                  { utc_usecs=0LL; }
      UTC_Time(const UTC_Time  & utc) { utc_usecs=0LL; utc_usecs=utc.utc_usecs; }
      UTC_Time(const long long & utc) { utc_usecs=utc; }
     ~UTC_Time() {}

      const UTC_Time & operator = (const UTC_Time  & utc) { utc_usecs=utc.utc_usecs; return(*this); }
      const UTC_Time & operator = (const long long & utc) { utc_usecs=utc;           return(*this); }

      UTC_Time operator + (const UTC_Time & utc) { UTC_Time tmp(*this); tmp.utc_usecs+=utc.utc_usecs; return(tmp); }
      UTC_Time operator - (const UTC_Time & utc) { UTC_Time tmp(*this); tmp.utc_usecs-=utc.utc_usecs; return(tmp); }

      operator long long () const { return(utc_usecs); }

      int operator == (const UTC_Time & utc) const { return ( utc_usecs==utc.utc_usecs ? 1 : 0 ); }
      int operator != (const UTC_Time & utc) const { return ( utc_usecs!=utc.utc_usecs ? 1 : 0 ); }
      int operator <  (const UTC_Time & utc) const { return ( utc_usecs< utc.utc_usecs ? 1 : 0 ); }
      int operator <= (const UTC_Time & utc) const { return ( utc_usecs<=utc.utc_usecs ? 1 : 0 ); }
      int operator >  (const UTC_Time & utc) const { return ( utc_usecs> utc.utc_usecs ? 1 : 0 ); }
      int operator >= (const UTC_Time & utc) const { return ( utc_usecs>=utc.utc_usecs ? 1 : 0 ); }

      long long  Current_UTC_Time (void);
      void       Current_UTC_Time (int &yy, int &yday, int &mo, int &dd, int &hh, int &mm, double &ss);
      
      double     UTC_Seconds      (void) const { return(((double) utc_usecs)/1000000.0); }
      long long  UTC_Microsecs    (void) const { return(utc_usecs); }

      long long  Elapsed_Usecs    (const UTC_Time &utc) const { return(         (utc.utc_usecs-utc_usecs)           ); }
      double     Elapsed_Msecs    (const UTC_Time &utc) const { return( (double)(utc.utc_usecs-utc_usecs)/1000.0    ); }
      double     Elapsed_Secs     (const UTC_Time &utc) const { return( (double)(utc.utc_usecs-utc_usecs)/1000000.0 ); }
      
      void       UTC_To_HHMMSS    (int & hh, int & mm, double & ss) const;
      void       UTC_To_Calendar  (int &yy, int &yday, int &mo, int &dd, int &hh, int &mm, double &ss) const;
      void       Calendar_To_UTC  (const int yy, const int mo, const int dd, const int hh, const int mm, const double ss);
};

//------------------------------------------------------------------------------------
// Function prototypes...
//------------------------------------------------------------------------------------

extern long long  Current_UTC_Time (void);
extern void       UTC_To_HHMMSS    (const long long utc, int & hh, int & mm, double & ss);

#endif
