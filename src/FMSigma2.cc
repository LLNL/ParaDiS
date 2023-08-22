/**************************************************************************
 *
 *      Module:       FMSigma2.c
 *      Description:  This module contains the base FMSigma2 entry
 *                    function used by the fast multipole code.
 *
 *      Includes functions:
 *          FMSigma2()
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "FM.h"


/**************************************************************************
 *
 *      Function;     FMSigma2
 *
 *      Description:  Simple dispatch function to call the appropriate
 *                    function based on the compile-time flags.
 *
 *      Parameters:   See descriptions for functions being invoked
 *
 *************************************************************************/
void FMSigma2(real8 mu, real8 nu, int norder, real8 Eeta[],
              matrix sigmatot, real8 pcterms[], int pcpows[][3])
{
#ifdef USE_RATIONAL_SEG_FORCES
        FMSigma2Rational(mu, nu, norder, Eeta, sigmatot, pcterms, pcpows);
#else
        FMSigma2NonRational(mu, nu, norder, Eeta, sigmatot, pcterms, pcpows);
#endif
        return;
}
