/*************************************************************************
 *
 *      Module:       ElasticConstants.c
 *      Description:  Contains functions to interpolate various
 *                    material specific isothermal elastic constants
 *                    for specific temperature and pressure
 *                    combinations based on data files containing
 *                    the elastic constants sampled at various
 *                    pressures and temperatures.
 *
 *      Data file format:
 *          <temp1> <pressure1> <const value>
 *          <temp1> <pressure2> <const value>
 *          ...
 *          <temp1> <pressureN> <const value>
 *          <temp2> <pressure1> <const value>
 *          <temp2> <pressure2> <const value>
 *          ...
 *          <tempN> <pressureN> <const value>
 *
 ************************************************************************/
#include <stdio.h>
#include <sys/param.h>

#include "mpi_portability.h"

#include "Home.h"

#define C11_FILE_SUFFIX  "c11.dat"
#define C12_FILE_SUFFIX  "c12.dat"
#define C44_FILE_SUFFIX  "c44.dat"
#define BURG_FILE_SUFFIX "burg.dat"

#define BOHR_RADIUS  0.529   /* Angstroms */


/*------------------------------------------------------------------------
 *
 *      Function:     LookupElasticConst
 *      Description:  Use the provided temperature and pressure to
 *                    look up one of the elastic constants in the
 *                    specified data file.  If the constant is not
 *                    defined for the exact temp/pressure combination
 *                    the value will be interpolated from nearby values.
 *
 *      Arguments:
 *          fileName        name of the file containing the elastic constant
 *                          at various temperature/pressure combinations.
 *          targetTemp      material temperature (in degress Kelvin)
 *          targetPressure  external pressure (in Mbar)
 *          value           pointer to location in which constant is
 *                          returned to the caller.
 *
 *----------------------------------------------------------------------*/
static int LookupElasticConst(char *fileName, real8 targetTemp,
                              real8 targetPressure, real8 *value)
{
        int    i, numRows, numCols, status;
        real8  **columnData;
        real8  *p, *t, *c;

        numRows = 0;
        numCols = 3;
        columnData = (real8 **)NULL;

        ReadTabulatedData(fileName, numCols, &columnData, &numRows);

        t = columnData[0];   /* temps     */
        p = columnData[1];   /* pressures */
        c = columnData[2];   /* constants */

/*
 *      Use a bilinear interpolator to calculate the necessary
 *      constant value for the provided pressure and temperature.
 */
        status = InterpolateBL(t, p, c, numRows, targetTemp,
                               targetPressure, value);

        if (status != 1) {
            Fatal("LookupElasticConst(): Error interpolating constant\n"
                  "from table %s\nfor pressure = %e, temperature = %e",
                  fileName, targetPressure, targetTemp);
        }
/*
 *      And just free up the temporary arrays...
 */
        for (i = 0; i < numCols; i++) {
            free(columnData[i]);
        }

        free(columnData);

        return(1);
}


/*------------------------------------------------------------------------
 *
 *      Function:     GetElasticConstants
 *      Description:  If the caller provided the name of a set of
 *                    files with material specific elastic constants,
 *                    calculate the elastic constants, shear modulus,
 *                    poisson ratio and burgers vector magnitude for
 *                    the specified temperature and pressure.
 *
 *      Returns:   1 on success
 *                 0 in all other cases
 *
 *----------------------------------------------------------------------*/
int GetElasticConstants(Param_t *param)
{
        real8 c11, c12, c44;
        real8 targetTemp, targetPressure;
        char  *baseFile;
        char  fileName[MAXPATHLEN];

/*
 *      If the caller didn't provide files with the elastic constants,
 *      just return and use default values for the constants.  If a
 *      base file name has been provided, data file names will be of
 *      the form:
 *            <baseFile>.c11.dat
 *            <baseFile>.c12.dat
 *            <baseFile>.c44.dat
 *            <baseFile>.burg.dat
 */
        baseFile = param->ecFile;

        if (baseFile[0] == 0) {
            return(0);
        }

/*
 *      Convert pressure from units of Pa to Mbar
 */
        targetPressure = param->pressure / 1.0e+11;
        targetTemp = param->TempK;

/*
 *      C11 elastic constant is returned in units of Mbars, convert
 *      it to units of Pa for consistency with the rest of the code.
 */
        sprintf(fileName, "%s.%s", baseFile, C11_FILE_SUFFIX);
        LookupElasticConst(fileName, targetTemp, targetPressure, &c11);

        c11 *= 1.0e+11;

/*
 *      C12 elastic constant is returned in units of Mbars, convert
 *      it to units of Pa for consistency with the rest of the code.
 */
        sprintf(fileName, "%s.%s", baseFile, C12_FILE_SUFFIX);
        LookupElasticConst(fileName, targetTemp, targetPressure, &c12);

        c12 *= 1.0e+11;

/*
 *      C44 elastic constant is returned in units of Mbars, convert
 *      it to units of Pa for consistency with the rest of the code.
 */
        sprintf(fileName, "%s.%s", baseFile, C44_FILE_SUFFIX);
        LookupElasticConst(fileName, targetTemp, targetPressure, &c44);

        c44 *= 1.0e+11;

/*
 *      Find the atomic volume per atom, use it to get the burgers
 *      magnitude in angstroms and convert to meters.
 */
        sprintf(fileName, "%s.%s", baseFile, BURG_FILE_SUFFIX);
        LookupElasticConst(fileName, targetTemp, targetPressure,
                           &param->burgMag);

        param->burgMag  = pow(2.0*param->burgMag, 1.0/3.0) *
                          BOHR_RADIUS * sqrt(3.0) / 2.0;
        param->burgMag *= 1.0e-10;

/*
 *      Now set the appropriate poisson ratio and shear modulus
 *      for the given temp/pressure.
 */
        param->shearModulus = (c11 - c12 + 3.0*c44) / 5.0;
        param->pois         = c12 / (2.0 * (c12 + c44));

        return(1);
}
