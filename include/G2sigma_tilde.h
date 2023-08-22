#pragma once

#ifndef _PDS_G2_SIGMA_TILDE_H
#define _PDS_G2_SIGMA_TILDE_H

/***************************************************************************
 *
 *      Module:       G2sigma_tilde.h
 *      Description:  This header contains prototypes for functions
 *                    used only in support of the anisotropic FMM.
 *
 ***************************************************************************/

#ifdef ANISOTROPIC
#include "Home.h"

void eta2eta_tilde(int norder,
		   real8 C[3][3][3][3],
		   real8 eta[],
		   real8 (*eta_tilde)[3][3]);

void alpha_tilde2alpha(int uorder,
		       real8 C[3][3][3][3],
		       real8 alpha_tilde[][3][3],
		       real8 alpha[][3][3]);

void G2sigma_tilde_deriv(int norder_array, int norder_eval,
                         int ndx, int ndy, int ndz,
			 real8 Gderiv[][3][3],
			 real8 eta_tilde[][3][3],
			 real8 sigma_tilde[3][3]);

#endif  // ifdef ANISOTROPIC
#endif  // ifndef _G2sigma_tilde_h
