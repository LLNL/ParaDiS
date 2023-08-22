#pragma once

#ifndef _OPENMP_PORTABILITY_H
#define _OPENMP_PORTABILITY_H

#ifdef _OPENMP
//--------------------------------------------------
// OpenMP is enabled - activate the OMP macros...
//--------------------------------------------------

#include <omp.h>

#define OMP_INIT_LOCK(a)      omp_init_lock    (a)
#define OMP_DELETE_LOCK(a)    omp_destroy_lock (a)
#define OMP_LOCK(a)           omp_set_lock     (a)
#define OMP_UNLOCK(a)         omp_unset_lock   (a)

#else
//--------------------------------------------------
// OpenMP is disabled - define stubs...
//--------------------------------------------------

#define OMP_INIT_LOCK(a)
#define OMP_DELETE_LOCK(a)
#define OMP_LOCK(a)
#define OMP_UNLOCK(a)

#endif  // _OPENMP

#endif  // _OPENMP_PORTABILITY_H
