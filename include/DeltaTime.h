#pragma once

#ifndef _PDS_DELTA_TIME_H
#define _PDS_DELTA_TIME_H

#include "Home.h"

//-----------------------------------------------------------------------------------------------

extern void DeltaTime_Save    (const Home_t *home);
extern void DeltaTime_Close   (const Home_t *home);

extern void Timestep_Cut_Save (const Home_t *home, Node_t *node);

#endif
