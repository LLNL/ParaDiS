#pragma once

#ifndef _PDS_TIMESTEP_H
#define _PDS_TIMESTEP_H

//--------------------------------------------------------------------------
//  Module:       Timestep.h
//  Description:  This header is mainly just a dumping ground for
//                the miscellaneous function prototypes and definitions
//                related to the timestep intergrators.
//--------------------------------------------------------------------------

#if defined(USE_KINSOL) || defined(USE_ARKODE)
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>
#endif

#if defined(USE_ARKODE)
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_ls.h>
#endif

#include "Home.h"

void ForwardEulerIntegrator   (Home_t *home);
void TrapezoidIntegrator      (Home_t *home);
void TrapezoidIntegratorMulti (Home_t *home);
void TrapezoidIntegratorKINSOL(Home_t *home);
void ARKodeIntegrator         (Home_t *home);

#ifdef USE_KINSOL
#include "ParadisSUNDIALS.h"
int FixedPt_Fcn(ParaDiS_Vector solutionVec, ParaDiS_Vector fVal, void *vHome);
int Newton_Fcn (ParaDiS_Vector solutionVec, ParaDiS_Vector fVal, void *vHome);
#endif

#ifdef USE_ARKODE
#include "ParadisSUNDIALS.h"
int  ARKodeIntegrator_RhsFi(real8 t, ParaDiS_Vector solutionVec, ParaDiS_Vector fVal, void *vHome);
void CreateARKode          (Home_t *home);
void FreeARKode            (Home_t *home);
void StatsARKode           (Home_t *home);
#endif

#endif  // _PDS_TIMESTEP_H
