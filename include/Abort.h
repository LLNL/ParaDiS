#pragma once

#ifndef _PDS_PARADIS_ABORT_H
#define _PDS_PARADIS_ABORT_H

#include "Home.h"

extern  int  paradis_abort;

extern  void ParadisAbort_Set (const int abort);

extern  int  ParadisAbort     (void);
extern  int  ParadisAbort     (Home_t *home);

#endif
