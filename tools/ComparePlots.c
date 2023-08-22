/*
 *      Module:  ComparePlots.c
 *
 *      Description:  Compares data points from a 2-D plot to a
 *                    baseline plot.  If any data points in the plot
 *                    differ from the baseline plot by more than the
 *                    specified tolerance, the offending data points
 *                    are written into an error file and the utility
 *                    will exit with an error code indicating the
 *                    plots did not match.
 *
 *      Command line options (required, not really optional!)
 *
 *          -b <baselinePlot>    baseline plot file
 *          -c <currentPlot>     current plot file
 *          -e <errorFile>       error plot file
 *          -f1 <index1>         index (1 offset) of field defining the first
 *                               plot axis along which values will be
 *                               monotonically increasing
 *          -f2 <index2>         index (1 offset) of field defining the values
 *                               in the second plot axis
 *          -p <errorTol>        percent error tolerance between plot and
 *                               baseline points
 *          -t <numColumns>      total number of fields in plot file
 *
 *      Exit codes:
 *
 *          0 == No problems
 *          1 == processing error (i.e. missing files, command lines opts, etc
 *          2 == plots do not match within error tolerance
 */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/param.h>


void Usage(char *progname)
{
        fprintf(stdout, "\n");
        fprintf(stdout, "Usage:  %s -b <basePlotFile> -c <currPlotFile>\n",
                progname);
        fprintf(stdout, "                -e <errPlotFile> -f1 <fieldNum>\n");
        fprintf(stdout, "                -f2 <fieldNum> -p <prcntErr> \n");
        fprintf(stdout, "                -t <numFields>\n");
        fprintf(stdout, "\n");
        fprintf(stdout, "where:\n");
        fprintf(stdout, "\n");
        fprintf(stdout, "    -b <basePlotFile>  specifies the name of the\n");
        fprintf(stdout, "                       baseline plot file\n");
        fprintf(stdout, "    -c <currPlotFile>  specifies the name of the\n");
        fprintf(stdout, "                       plot file to be compared to\n");
        fprintf(stdout, "                       the baseline plot\n");
        fprintf(stdout, "    -e <errPlotFile>   specifies the name of the \n");
        fprintf(stdout, "                       file into which to dump the\n");
        fprintf(stdout, "                       lines from the current plot\n");
        fprintf(stdout, "                       that do not match the base\n");
        fprintf(stdout, "                       plot within tolerance\n");
        fprintf(stdout, "    -f1 <fieldNum>     index (1 offset) of the \n");
        fprintf(stdout, "                       plot file field correspond-\n");
        fprintf(stdout, "                       ing to the primary axis\n");
        fprintf(stdout, "                       along which values will be\n");
        fprintf(stdout, "                       monotonically increasing\n");
        fprintf(stdout, "    -f2 <fieldNum>     index (1 offset) of the \n");
        fprintf(stdout, "                       plot file field correspond-\n");
        fprintf(stdout, "                       ing to the secondary axis\n");
        fprintf(stdout, "    -p <prcntErr>      error tolerance (as %)\n");
        fprintf(stdout, "    -t <numFields>     total number of fields(i.e.\n");
        fprintf(stdout, "                       columns) in the plot files\n");
        fprintf(stdout, "\n");

        exit(1);
}


/*
 *      Function:    GetPlotFields()
 *
 *      Description: Read the next line of data from a plot file and copy
 *                   all of the values into an array to be returned to
 *                   the caller.  If a complete data line is not read in
 *                   the array of values from the caller is not updated.
 * 
 *                   Note: The values in the field defining the primary axis
 *                         of the plot should always be increasing, however,
 *                         if the code was stopped and restarted from a
 *                         previous restart file, there may be temporary
 *                         decreases in those values.  If so, just skip them.
 *      Arguments:
 *
 *          IN:  fp        File pointer from which to read the plot data
 *
 *          IN:  numFields Number of expected columns of data per line
 *
 *          IN:  axisField Index (0 offset) of the data column representing
 *                         the primary axis of the plot.
 *
 *          OUT: fieldList Array inwhich <numField> values read from the
 *                         next data line in the file are returned to the
 *                         caller.
 *
 *      Returns: 0  if no data read (probably EOF)
 *               1  if values returned as requested
 */
int GetPlotFields(FILE *fp, int numFields, int axisField, double *fieldList)
{
        int    index;
        double *newFieldList;
        char   *token, *savePtr, *s;
        char   line[512];

        newFieldList = calloc(1, numFields * sizeof(double));

        while (1) {

            line[0] = 0;

            if ((s = fgets(line, sizeof(line)-1, fp)) == (char *)NULL) {
                free(newFieldList);
                return(0);
            }

            line[strlen(line)-1] = 0;

/*
 *          Skip comment lines and blank lines
 */
            if (*s == '#') continue;

            for ( ; *s != 0; s++) {
                if ((*s != ' ') && (*s != '\t') && (*s != '\n')) break;
            }

            if (*s == 0) continue;

/*
 *          Read in all the fields of the input line.  If the data line
 *          did not contain at least the expected number of columns
 *          just print an error and return a failure.
 */
            savePtr = (char *)NULL;
            s = line;

            for (index = 0; index < numFields; index++) {
                token = strtok_r(s, " 	", &savePtr);
                s = (char *)NULL;
                if (token == (char *)NULL) {
                    fprintf(stderr,
                            "Error parsing file: too few columns\n");
                    free(newFieldList);
                    return(0);
                }
                newFieldList[index] = atof(token);
            }

/*
 *          If the new value in the field defining the primary axis
 *          is larger than the previous value, it's usable, so we
 *          can return the data to the caller.
 */
            if (newFieldList[axisField] > fieldList[axisField]) {
                break;
            }
        }

/*
 *      Copy out the data values into the caller's buffer
 */
        for (index = 0; index < numFields; index++) {
            fieldList[index] = newFieldList[index];
        }

        free(newFieldList);

        return(1);
}


main(int argc, char *argv[])
{
        int    status;
        int    opt, index;
        int    numFields, primaryAxis, secondaryAxis;
        double cmpVal, factor, prcntError;
        double maxErrTol, maxError;
        double basePrimaryMin, currPrimaryMin;
        double primaryRangeMin, primaryRangeMax;
        double currPrimary, currSecondary;
        double basePrimaryValMax, currPrimaryValMax;
        double primaryAxisErr;
        double *plotValues, *baseValuesMin, *baseValuesMax;
        char   *baselinePlotFile, *currPlotFile, *errPlotFile;
        FILE   *fpBaseline, *fpCurrent, *fpErrPlot;


        status = 0;

        fpErrPlot = (FILE *)NULL;

        baselinePlotFile = (char *)NULL;
        currPlotFile = (char *)NULL;
        errPlotFile = (char *)NULL;

        numFields      = -1;
        primaryAxis    = -1;
        secondaryAxis  = -1;

        maxErrTol      = -1.0;

/*
 *      Parse the command line
 */
        for (opt = 1; opt < argc; opt++) {
            if (!strcmp(argv[opt], "-f1")) {
                if (opt >= (argc - 1)) Usage(argv[0]);
                primaryAxis = atoi(argv[++opt]) - 1;
            } else if (!strcmp(argv[opt], "-f2")) {
                if (opt >= (argc - 1)) Usage(argv[0]);
                secondaryAxis = atoi(argv[++opt]) - 1;
            } else if (!strcmp(argv[opt], "-p")) {
                if (opt >= (argc - 1)) Usage(argv[0]);
                maxErrTol = atof(argv[++opt]);
            } else if (!strcmp(argv[opt], "-t")) {
                if (opt >= (argc - 1)) Usage(argv[0]);
                numFields = atoi(argv[++opt]);
            } else if (!strcmp(argv[opt], "-b")) {
                if (opt >= (argc - 1)) Usage(argv[0]);
                baselinePlotFile = argv[++opt];
            } else if (!strcmp(argv[opt], "-c")) {
                if (opt >= (argc - 1)) Usage(argv[0]);
                currPlotFile = argv[++opt];
            } else if (!strcmp(argv[opt], "-e")) {
                if (opt >= (argc - 1)) Usage(argv[0]);
                errPlotFile = argv[++opt];
            } else {
                Usage(argv[0]);
            }
        }

/*
 *      If file names were not provided return an error
 */
        if ((baselinePlotFile == (char *)NULL) ||
            (currPlotFile == (char *)NULL) ||
            (errPlotFile == (char *)NULL)) {
            fprintf(stdout, "Missing input or output file names!\n");
            Usage(argv[0]);
        }

/*
 *      If the caller did not provide the other needed command line
 *      info, return an error
 */
        if ((numFields < 0)     ||
            (primaryAxis < 0)   ||
            (secondaryAxis < 0) ||
            (maxErrTol < 0.0)) {
            fprintf(stdout, "Missing command line parameters!\n");
            Usage(argv[0]);
        }

/*
 *      Open the input files
 */
        if ((fpBaseline = fopen(baselinePlotFile, "r")) == (FILE *)NULL) {
            fprintf(stdout, "Error %d opening %s\n", errno, baselinePlotFile);
            exit(1);
        }

        if ((fpCurrent = fopen(currPlotFile, "r")) == (FILE *)NULL) {
            fprintf(stdout, "Error %d opening %s\n", errno, currPlotFile);
            exit(1);
        }

        plotValues = (double *)calloc(1, numFields * sizeof(double));
        baseValuesMin = (double *)calloc(1, numFields * sizeof(double));
        baseValuesMax = (double *)calloc(1, numFields * sizeof(double));

/*
 *      Read the first data point from the current plot file and
 *      the first two points from the baseline plot file.
 *
 *      Basically, we'll try to bracket each data point from the current
 *      plot with data points from the baseline then estimate the baseline
 *      plot value at the same position on the primary axis as the current
 *      data point.
 */
        GetPlotFields(fpCurrent, numFields, primaryAxis, plotValues);
        GetPlotFields(fpBaseline, numFields, primaryAxis, baseValuesMin);
        GetPlotFields(fpBaseline, numFields, primaryAxis, baseValuesMax);

/*
 *      Save the initial values on the primary axis from the baseline
 *      and current plots.
 */
        basePrimaryMin = baseValuesMin[primaryAxis];
        currPrimaryMin = plotValues[primaryAxis];

        maxError = 0.0;

        primaryRangeMin = baseValuesMin[primaryAxis];
        primaryRangeMax = baseValuesMax[primaryAxis];

        currPrimary   = plotValues[primaryAxis];
        currSecondary = plotValues[secondaryAxis];

        while (1) {

/*
 *          If the current point on the primary axis is less than the 
 *          lower bounding point on that axis from the baseline plot
 *          we don't want to try dealing with that, so skip the data point.
 */
            if (currPrimary < primaryRangeMin) {

                if (!GetPlotFields(fpCurrent, numFields, primaryAxis,
                                   plotValues)) {
                    break;
                }

                currPrimary   = plotValues[primaryAxis];
                currSecondary = plotValues[secondaryAxis];

                continue;
            }

/*
 *          If the current point on the primary axis exactly matches the 
 *          lower bounding point on that axis from the baseline plot
 *          just use the exact secondary axis value from the baseline
 */
            else if (currPrimary == primaryRangeMin) {
                cmpVal = baseValuesMin[secondaryAxis];

            }

/*
 *          If the current point on the primary axis exceeds the 
 *          uper bounding point on that axis from the baseline plot
 *          so far, read the next data line from the baseline.
 */
            else if (currPrimary >= primaryRangeMax) {

                for (index = 0; index < numFields; index++) {
                    baseValuesMin[index] = baseValuesMax[index];
                }

                primaryRangeMin = baseValuesMin[primaryAxis];

                if (!GetPlotFields(fpBaseline, numFields, primaryAxis,
                                   baseValuesMax)) {
                    break;
                }

                primaryRangeMax = baseValuesMax[primaryAxis];

                continue;
            }

/*
 *          The current plot value is bracketed by data points in the baseline
 *          so find the value for the secondary axis along the baseline plot
 *          that corresponds to the current position along the primary axis.
 *          (assuming linear change between points on the secondary axis in
 *          the baseline data).
 */
            else {

                factor = (currPrimary - primaryRangeMin) /
                         (primaryRangeMax - primaryRangeMin);

                cmpVal = baseValuesMin[secondaryAxis] +
                         factor * (baseValuesMax[secondaryAxis] -
                                   baseValuesMin[secondaryAxis]);
            }

            prcntError = fabs((cmpVal - currSecondary) / cmpVal * 100.0);
            maxError = MAX(maxError, prcntError);

/*
 *          If the data point is too far from the baseline value
 *          write the data point to a file we can use later to
 *          visualize the results.
 */
            if (prcntError > maxErrTol) {

                if (fpErrPlot == (FILE *)NULL) {
                    fpErrPlot = fopen(errPlotFile, "w");
                }

                for (index = 0; index < numFields; index++) {
                    fprintf(fpErrPlot, "%e ", plotValues[index]);
                }

                fprintf(fpErrPlot, "\n");

                status = 2;
            }

/*
 *          Read the next data point from the results file.  When
 *          there's no more data, we're done.
 */
            if (!GetPlotFields(fpCurrent,numFields,primaryAxis,plotValues)) {
                break;
            }

            currPrimary   = plotValues[primaryAxis];
            currSecondary = plotValues[secondaryAxis];

        }

/*
 *      One of the files may have ended before the other, so read
 *      to the end of each file in order to figure out how much difference
 *      there was 
 */
        while (1) {
            if (!GetPlotFields(fpCurrent,numFields,primaryAxis,plotValues)) {
                break;
            }
        }

        while (1) {
            if (!GetPlotFields(fpBaseline, numFields, primaryAxis,
                               baseValuesMax)) {
                break;
            }
        }

        fclose(fpBaseline);
        fclose(fpCurrent);

/*
 *      Check if the plot extends further or less than the baseline
 *      along the primary axis
 */
        basePrimaryValMax = baseValuesMax[primaryAxis];;
        currPrimaryValMax = plotValues[primaryAxis];;

        if (basePrimaryValMax > currPrimaryValMax) {
            primaryAxisErr = 100 * (basePrimaryValMax - currPrimaryValMax) /
                          basePrimaryValMax;
        } else {
            primaryAxisErr = 100 * (currPrimaryValMax - basePrimaryValMax) /
                          currPrimaryValMax;
        }

/*
 *      Print out some info from the comparison...
 */
        fprintf(stdout, "\n");
        fprintf(stdout, "Baseline file: %s\n", baselinePlotFile);
        fprintf(stdout, "Results file:  %s\n", currPlotFile);
        fprintf(stdout, "Maximum percent error: %5.2lf%%\n", maxError);
        fprintf(stdout, "Baseline max along primary axis %s results "
                "by %5.2lf%%\n", (basePrimaryValMax > currPrimaryValMax) ?
                "exceeds" : "is smaller than", primaryAxisErr);
        fprintf(stdout, "\n");
        fprintf(stdout, "Baseline primary axis: minimum %e, maximum %e\n",
                basePrimaryMin, basePrimaryValMax);
        fprintf(stdout, "Results primary axis:  minimum %e, maximum %e\n",
                currPrimaryMin, currPrimaryValMax);
        fprintf(stdout, "\n");

/*
 *      If the plots did not sufficiently match, let the user know.
 */

        if (fpErrPlot != (FILE *)NULL) {
            fclose(fpErrPlot);
            fprintf(stdout, "Data points exceeding threshold written to %s\n",
                    errPlotFile);
            fprintf(stdout, "\n");
        }

        exit(status);
}
