#pragma once

#ifndef _PDS_INIT_H
#define _PDS_INIT_H

/*****************************************************************************
 *
 *	Init.h		Prototypes for the initialization routines
 *
 ****************************************************************************/

#include "Home.h"
#include "InData.h"

void    GetMeltTemp          (Param_t *param);
void    InitCellDomains      (Home_t  *home);
void    InitCellNatives      (Home_t  *home);
void    InitCellNeighbors    (Home_t  *home);
Home_t *InitHome             (void);
void    Initialize           (Home_t  *home   , int argc, char *argv[] );
int     OpenDir              (Home_t  *home);
void    ParadisInit          (Home_t **homeptr, int argc, char *argv[] );
void    RecvInitialNodeData  (Home_t  *home);
void    SendInitialNodeData  (Home_t  *home, InData_t *inData, int *msgCount, int **nodeLists, int *listCounts, int *nextAvailableTag);
void    SetRemainingDefaults (Home_t  *home);
void    VerifyBurgersVectors (Home_t  *home);

#endif  // _PDS_INIT_H
