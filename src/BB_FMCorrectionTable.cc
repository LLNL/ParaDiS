/***************************************************************************
 *
 *      Author:       Chao Chen, Eric Darve and Sylvie Aubry
 *      Module:       CorrectionTable.c
 *      Description:  Contains functions for creating a table used
 *                    in conjuction with the Fast Multipole functions
 *                    to adjust the calculated stress to allow for
 *                    periodic images of the problem space.
 *
 *      Includes public functions:
 *
 *          bbFMM_CorrectionTableInit()
 *          bbFMM_CreateCorrectionTable()
 *          bbFMM_DoTableCorrection()
 *          bbFMM_FreeCorrectionTable()
 * 
 *      Includes private functions:
 *
 *          ComputeWeightsPBC
 *
***************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "Home.h"
#include "FM.h"

#ifdef PARALLEL
#include  <mpi.h>
#endif

#ifdef BBFMM
#include "FM.h"

static real8 *correctionTbl = NULL;

void bbFMM_DoTableCorrection(Home_t *home)
{
  int n;
  real8 homogen;

  n = bbfmm->n;
  homogen = bbfmm->homogen;
  Stanford_DoTableCorrection(home,n,homogen);
}

#endif
