#pragma once

#ifndef _PDS_WRITEPROP_H
#define _PDS_WRITEPROP_H

/*************************************************************************
 *
 *  WriteProp.h - define parameters for writing property-time curve files
 *
 *************************************************************************/

#include "Home.h"

extern void WriteProp                         (Home_t *home, int   property    );

extern void WriteAllepsFileDesc               (Home_t *home, char *baseFileName);
extern void WriteDensityDeltaFileDesc         (Home_t *home, char *baseFileName);
extern void WriteDensityFileDesc              (Home_t *home, char *baseFileName);
extern void WriteDensityFileGnu               (Home_t *home, char *baseFileName);
extern void WriteEpsdotFileDesc               (Home_t *home, char *baseFileName);
extern void WriteStressPlasticStrainFileDesc  (Home_t *home, char *baseFileName);
extern void WriteStressPlasticStrainFileGnu   (Home_t *home, char *baseFileName);
extern void WriteStressTotalStrainFileDesc    (Home_t *home, char *baseFileName);
extern void WriteTimePlasticStrainFileDesc    (Home_t *home, char *baseFileName);
extern void WriteTimePlasticStrainFileGnu     (Home_t *home, char *baseFileName);
extern void WriteTotalStressFileDesc          (Home_t *home, char *baseFileName);

#ifdef CALCENERGY
extern void WriteTimeEnergyFileDesc           (Home_t *home, char *baseFileName);
extern void WriteEnergyPlasticStrainFileDesc  (Home_t *home, char *baseFileName);
extern void WriteEnergyStrainFileDesc         (Home_t *home, char *baseFileName);
#endif

#define DENSITY         1
#define EPS             2
#define ALL_EPS         3
#define EPSDOT          4
#define DENSITY_DELTA   5
#define ENERGY          6

#endif
