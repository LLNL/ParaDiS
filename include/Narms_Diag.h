#pragma once

#ifndef _PDS_NARMS_DIAG_H
#define _PDS_NARMS_DIAG_H

class Narms_Diag
{
   public : 
      static int  *narms_cnts;   //< array of accumulated nodal arm counts (across all domains)
      static int   narms_max ;   //< maximum number of arms for a given node
      static int   narms_indx;   //< the current cycle index
      static int   narms_freq;   //< frequency of diagnostic saves

   public : 
      Narms_Diag(void) {}
     ~Narms_Diag()     {}

      void Init        (const Home_t *home);
      void Cycle_Accum (const Home_t *home);
      void Save        (void);
      void Save_GNU    (void);
};

//------------------------------------------------------------------------------------------------------------

extern void Narms_Diag_Save (const Home_t *home, const int flush=0);

#endif  // _PDS_NARMS_DIAG_H
