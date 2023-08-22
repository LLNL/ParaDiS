/****************************************************************************
 *
 *      Module:       InclusionGen.c
 *
 *      Description:  This utility is used for creating a file defining
 *                    eshelby inclusions in a format that can be understood
 *                    by ParaDiS.  The utility assumes the problem space
 *                    to be defined as a rectangular prism with sides of
 *                    arbitrary lengths and permits any combination of
 *                    free surfaces and periodic boundaries among the
 *                    three dimensions.
 *
 *
 *      For a description of the comand line options and default
 *      configuration values, execute:
 *
 *          inclusiongen -help
 *
 ****************************************************************************/
typedef double real8;


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/param.h>


#define DotProduct(vec1,vec2)       \
            ((vec1[0])*(vec2[0]) +  \
             (vec1[1])*(vec2[1]) +  \
             (vec1[2])*(vec2[2]))


/*
 *      Define globals and default values to use if not provided by the user
 */
               int kodd = 0;

static int    allowOverlap = 0;
static int    numToCreate = 100;
static int    isFreeSurf[3] = {0, 0, 0};
static int    pbcVal = 7;
static int    printDebug = 0;
static real8 radiusMin = 20.0;
static real8 radiusMax = 50.0;
static real8 VolumeFraction = 0.0;
static real8 simSize[3] = {35000.0, 35000.0, 35000.0};
static real8 simSizeInv[3];
static real8 simBoundsMin[3] = {-17500.0, -17500.0, -17500.0};
static real8 simBoundsMax[3] = { 17500.0,  17500.0,  17500.0};
static real8 strainField[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
static char   *inputFile = (char *)NULL;
static char   *outputFile = (char *)NULL;


#define X 0
#define Y 1
#define Z 2

/*
 *      Define structure of inclusion-related info we'll need
 */
typedef struct {
        int    inclusionID;
        real8 radius;
        real8 coordinates[3];
        real8 strainField[6];
} Inclusion_t;


/*
 *      Define the structure of info associed with each cell
 *      in the cellular grid which gets used if we need to
 *      check for overlapping inclusuions.
 */
typedef struct {
        int numInclusions;
        int arraySize;
        int *IDs;
} Cell_t;


/*-------------------------------------------------------------------------
 *
 *      Function:       Usage
 *      Description:    Print a brief synopsis of the acceptable
 *                      command line format enad exit.
 *
 *      Arguments:
 *          exeName   Name of the executable
 *
 *-----------------------------------------------------------------------*/
static void Usage(char *exeName)
{
        fprintf(stderr, "Usage:  %s [-a] [-d]\n", exeName);
        fprintf(stderr, "        [-i inputFile] [-h|-help]\n");
        fprintf(stderr, "        [-l lengthX[,lengthY,lengthZ]]\n");
        fprintf(stderr, "        [-n numInclusions]\n");
        fprintf(stderr, "        [-o outputFile] [-p pbcSpec]\n");
        fprintf(stderr, "        [-r radiusMin[,radiusMax]]\n");
        fprintf(stderr, "        [-s s00,s11,s22,s12,s02,s01]\n");

        exit(1);
}


/*-------------------------------------------------------------------------
 *
 *      Function:       PrintHelp
 *      Description:    Print a brief synopsis of the command line
 *                      format, a description of each option, and
 *                      print the values that would be used if
 *                      the utility was executed with the specified
 *                      options.
 *
 *      Arguments:
 *          exeName   Name of the executable
 *
 *-----------------------------------------------------------------------*/
static void PrintHelp(char *exeName)
{
        printf("Usage:  %s [-a] [-d]\n", exeName);
        printf("        [-i inputFile] [-h|-help]\n");
        printf("        [-l lengthX[,lengthY,lengthZ]]\n");
        printf("        [-n numInclusions]\n");
        printf("        [-o outputFile] [-p pbcSpec]\n");
        printf("        [-r radiusMin[,radiusMax]]\n");
        printf("        [-s s00,s11,s22,s12,s02,s01]\n");

        printf("\nwhere:\n");
        printf("    -a        Allows creation of overlapping inclusions.\n"
               "              By default, no overlap will be permitted\n");
        printf("    -d        Permits the utility to print some debug info\n"
               "              such as number of inclusions created and\n"
               "              overlaps encountered, etc. while the utility\n"
               "              is executing.\n");
        printf("    -i        Specifies the name of an input file containing\n"
               "              a set of inclusions to which the new inclusions\n"
               "              will be added\n");
        printf("    -h, -help Prints this help info\n");
        printf("    -l        Specifies the length of the simulation space.  \n");
        printf("    -n        Gives the number of new inclusions to create\n");
        printf("    -o        Gives the name of an output file\n");
        printf("    -p        Integer value specifying dimensions in which \n");
        printf("              periodic boundaries are enabled.  Values can \n");
        printf("              be bitwise OR'ed together to yield periodic \n");
        printf("              boundaries in any combination of dimensions.\n");
        printf("              Only applies if configuration type is\n");
        printf("              finitemixed\n");
        printf("                1 == periodic in X dimensions.\n");
        printf("                2 == periodic in Y dimensions.\n");
        printf("                3 == periodic in X and Y dimensions.\n");
        printf("                4 == periodic in Z dimensions.\n");
        printf("                5 == periodic in X and Z dimensions.\n");
        printf("                6 == periodic in Y and Z dimensions.\n");
        printf("                7 == periodic in all dimensions.\n");
        printf("    -r        Indicates the range of radii allowed for\n");
        printf("              new inclusions.  If maxRadius is not given\n");
        printf("              and all inclusions will be created with the\n");
        printf("              a radius equal to the minimum radius\n");
        printf("    -s        Specifies the 6 strain field components to\n");
        printf("              assign to the newly created inclusions\n");

        printf("\nDefault values:\n\n");
        printf("    Inclusions to create %d\n", numToCreate);
        printf("    Radius range:        %.1lf - %.1lf\n",
                radiusMin, radiusMax);
        printf("    Strain field used:   "
                "%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",
                strainField[0], strainField[1], strainField[2],
                strainField[3], strainField[4], strainField[5]);
        printf("    Overlap allowed:     %s\n",
                allowOverlap ? "YES" : "NO");
        printf("    Sim boundaries:       minBound / maxBound\n");
        printf("        X dimension      %.3lf / %.3lf\n", 
                simBoundsMin[X], simBoundsMax[X]);
        printf("        Y dimension      %.3lf / %.3lf\n", 
                simBoundsMin[Y], simBoundsMax[Y]);
        printf("        Z dimension      %.3lf / %.3lf\n", 
                simBoundsMin[Z], simBoundsMax[Z]);
        printf("    Boundary conditions:\n");
        printf("        X dimension      %s\n", 
                isFreeSurf[X] ? "Free surface" : "Periodic");
        printf("        Y dimension      %s\n", 
                isFreeSurf[Y] ? "Free surface" : "Periodic");
        printf("        Z dimension      %s\n", 
                isFreeSurf[Z] ? "Free surface" : "Periodic");
        printf("    Input file:          %s\n",
                inputFile == (char *)NULL ? "(NULL)" : inputFile);
        printf("    Output file:         %s\n",
                outputFile == (char *)NULL ? "(NULL)" : outputFile);
        printf("\n");

        exit(0);
}


/*-------------------------------------------------------------------------
 *
 *      Function:       IndicesToCellID
 *      Description:    Convert the X, Y and Z indices of a cell into
 *                      a unique integer ID.  Indices per dimension
 *                      are in the range:
 *                         -1 <= index <= number of cells in dimension
 *                      Cell indices less than zero or more than the
 *                      number of cell in the dimension are adjusted
 *                      to account for periodic boundaries.
 *
 *      Arguments:
 *          numCells  Three element array of number of cells in X, Y and
 *                    Z dimensions respectively.
 *          indices   Three element array of X, Y and Z indices of a cell
 *                    repsectively.  To account for periodic boundaries,
 *                    indices may span the range:
 *                       -1 <= indices[N] <= numCells[N]
 *
 *-----------------------------------------------------------------------*/
static int IndicesToCellID(int numCells[3], int indices[3])
{
        int x, y, z;

        x = indices[X];
        y = indices[Y];
        z = indices[Z];

/*
 *      Adjust indices that are outside the primary cell grid
 *      because of periodic boundaries.
 */
        if (x < 0) x = numCells[X] - 1;
        if (y < 0) y = numCells[Y] - 1;
        if (z < 0) z = numCells[Z] - 1;

        if (x >= numCells[X]) x = 0;
        if (y >= numCells[Y]) y = 0;
        if (z >= numCells[Z]) z = 0;

/*
 *      The cellIDs are calculated assuming cell indices vary fastest
 *      in the X dimension and slowest in the Z.
 */
        return(x + y * numCells[X] + z * numCells[X] * numCells[Y]);
}


/*-------------------------------------------------------------------------
 *
 *      Function:       AddToCell
 *      Description:    Add an inclusion to the list of inclusions
 *                      encompassed by a given cell.
 *
 *      Arguments:
 *          inclusionIndex Index in the inclusion list of the inclusion
 *                         to be added to the cell.
 *          indices        Three element array of X, Y and Z indices of a cell
 *                         respectively.  To account for periodic boundaries,
 *                         indices may span the range:
 *                            -1 <= indices[N] <= numCells[N]
 *          numCells       Three element array of number of cells in X, Y and
 *                         Z dimensions respectively.
 *          cellList       Array of cell structures: 1 per cell in grid.
 *
 *-----------------------------------------------------------------------*/
static void AddToCell(int inclusionIndex, int indices[3], int numCells[3],
                      Cell_t *cellList)
{
/*
 *      Determine the ID of the cell to which we'll add the inclusion.
 */
        int     cellID = IndicesToCellID(numCells, indices);
        Cell_t *cell   = &cellList[cellID];

/*
 *      If the cell's inclusion list does not have sufficient space
 *      for another inclusion, increase the size of the cell's
 *      inclusion list before appending the inclusion to the list.
 */
        if (cell->numInclusions == cell->arraySize) {
            cell->arraySize += 10;
            cell->IDs = (int *)realloc(cell->IDs, cell->arraySize *
                                       sizeof(int));
        }

        cell->IDs[cell->numInclusions] = inclusionIndex;
        cell->numInclusions += 1;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       CheckOverlap
 *      Description:    Read the next non-empty, non-comment line
 *
 *      Arguments:
 *          inc1,inc2  Pointers to 2 inclusions that need to be checked
 *                     for any overlapping volume.
 *
 *-----------------------------------------------------------------------*/
static int CheckOverlap(Inclusion_t *inc1, Inclusion_t *inc2)
{
        real8 x1, y1, z1;
        real8 x2, y2, z2;
        real8 distSquared, minDistSquared;

        x1 = inc1->coordinates[X];
        y1 = inc1->coordinates[Y];
        z1 = inc1->coordinates[Z];

        x2 = inc2->coordinates[X];
        y2 = inc2->coordinates[Y];
        z2 = inc2->coordinates[Z];

        minDistSquared = (inc1->radius + inc2->radius) *
                         (inc1->radius + inc2->radius);

        if (!isFreeSurf[X]) {
            x2 -= rint((x2-x1)*simSizeInv[X]) * simSize[X];
        }

        if (!isFreeSurf[Y]) {
            y2 -= rint((y2-y1)*simSizeInv[Y]) * simSize[Y];
        }

        if (!isFreeSurf[Z]) {
            z2 -= rint((z2-z1)*simSizeInv[Z]) * simSize[Z];
        }

        distSquared = (((x2-x1) * (x2-x1)) +
                       ((y2-y1) * (y2-y1)) +
                       ((z2-z1) * (z2-z1)));

        return(distSquared < minDistSquared);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     randm
 *      Description:  This is a special function for random number
 *                    generation on 32-bit machines that do not
 *                    support long integer multiplication and truncation.
 *                    The technique used is to do the multiplication and
 *                    addition in parts, by splitting all integers into a
 *                    'high' and a 'low' part.  The algorithm is exact,
 *                    and should give machine-independent results.
 *
 *                    The algorithm implemented is (following D.E. Knuth):
 *                    seed = seed*1592653589 + 453816691
 *                    if (seed.lt.0) seed = seed + 1 + 2147483647
 *                    Note that 1592653589 = 48603*2**15 + 30485
 *
 *      Returns:      real8 value in the range  0.0 <= value <= 1.0
 *
 *-------------------------------------------------------------------------*/
static real8 randm(int *seed)
{
        int ia, ib, i1, i2, i3;

        ia = *seed / 32768;
        ib = (*seed % 32768);
        i1 = ia * 30485;
        i2 = ib * 30485;
        i3 = ib * 48603;
        i1 = (i1 % 65536);
        i3 = (i3 % 65536);
        i1 = i1 + i3 + 13849 + i2 / 32768 + (ia % 2) * 32768;
        i2 = (i2 % 32768) + 12659;
        i1 = i1 + i2 / 32768;
        i2 = (i2 % 32768);
        i1 = (i1 % 65536);

        *seed = i1 * 32768 + i2;

        return(*seed * 4.65661287308E-10);
}


/*-------------------------------------------------------------------------
 *
 *      Function:       Getline
 *      Description:    Read the next non-empty, non-comment line
 *                      from the specified file stream and return
 *                      it to the caller.
 *
 *      Arguments:
 *              string          pointer to location at which to return
 *                              the line to the caller
 *              len             size (in bytes) of <string>
 *              fp              stream from which to read the data
 *
 *-----------------------------------------------------------------------*/
static void Getline(char *string, int len, FILE *fp)
{
        char    *s;

        while (1) {

                string[0] = 0;
                if ((s = fgets(string, len, fp)) == (char *)NULL) break;
                string[len-1] = 0;

/*
 *              If the line is a coment line, skip it. If
 *              it has any non-whitespace characters, return it to
 *              the caller, otherwise move to the next line.
 */
                if (*s =='#') continue;

                for ( ; *s != 0; s++) {
                        if ((*s != ' ') && (*s != '\t') && (*s != '\n')) break;
                }

                if (*s != 0) break;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       ReadInclusionFile
 *      Description:    Read inclusions from an existing inclusion
 *                      specification file and store all necessary
 *                      inclusion-related data in an array to be
 *                      returned to the caller.
 *
 *      Arguments:
 *          inputFile     Name of the inclusion specification file
 *                        to be read in.
 *          inclusionList Location in which to return to the caller a
 *                        pointer to the list of inclusions read in
 *          numInclusions Location in whcih to return to caller the
 *                        number of inclusions read in.
 *          largestRadius Pointer to the largest inclusion radius. If
 *                        any inclusions read have a radius that exceeds
 *                        this value, it will be updated before return
 *                        to the caller to contain the radius of the
 *                        largest inclusion in the input file.
 *
 *-----------------------------------------------------------------------*/
static void ReadInclusionFile(char *inputFile, Inclusion_t **inclusionList,
                       int *numInclusions, real8 *largestRadius)
{
        int         numAllocated = 0;
        int         numUsed = 0;
        real8       maxSize;
        char        *idStr, *radiusStr, *coordStr[3], *strainStr[6];
        char        inputLine[512];
        FILE        *fp;
        Inclusion_t *tmpList = (Inclusion_t *)NULL;

        fp = fopen(inputFile, "r");

        if (fp == (FILE *)NULL) {
            fprintf(stderr, "Open error on file %s!\n", inputFile);
            exit(1);
        }

        maxSize = *largestRadius;

/*
 *      Loop through each line of the input file.  Any comment
 *      lines or blank lines are ignored by Getline().  All other
 *      lines are considered to contain the specifications for
 *      a single inclusion and are parsed as such.
 */
        while (1) {

            Getline(inputLine, sizeof(inputLine) - 1, fp);

            if (inputLine[0] == 0) {
                break;
            }

/*
 *          Parse the line for a single inclusion's data
 */
            idStr        = strtok(inputLine, " \t");
            coordStr[X]  = strtok((char *)NULL, " \t");
            coordStr[Y]  = strtok((char *)NULL, " \t");
            coordStr[Z]  = strtok((char *)NULL, " \t");
            radiusStr    = strtok((char *)NULL, " \t");
            strainStr[0] = strtok((char *)NULL, " \t");
            strainStr[1] = strtok((char *)NULL, " \t");
            strainStr[2] = strtok((char *)NULL, " \t");
            strainStr[3] = strtok((char *)NULL, " \t");
            strainStr[4] = strtok((char *)NULL, " \t");
            strainStr[5] = strtok((char *)NULL, " \t");

/*
 *          If we did not get the expected number of data items
 *          for a single inclusion, tehre's a problem, so abort.
 */
            if ((idStr == (char *)NULL)        ||
                (radiusStr == (char *)NULL)    ||
                (coordStr[X] == (char *)NULL)  ||
                (coordStr[Y] == (char *)NULL)  ||
                (coordStr[Z] == (char *)NULL)  ||
                (strainStr[0] == (char *)NULL) ||
                (strainStr[1] == (char *)NULL) ||
                (strainStr[2] == (char *)NULL) ||
                (strainStr[3] == (char *)NULL) ||
                (strainStr[4] == (char *)NULL) ||
                (strainStr[5] == (char *)NULL)) {
                fprintf(stderr, "Parsing error on file %s!\n", inputFile);
                exit(1);
            }

/*
 *          Looks like we got a valid inclusion specification, so 
 *          insure the inclusion array has sufficient space, then
 *          add the inclusion to the list.
 */
            if (numUsed == numAllocated) {
                numAllocated += 500;
                tmpList = (Inclusion_t *)realloc(tmpList, numAllocated *
                                                 sizeof(Inclusion_t));
            }

            tmpList[numUsed].inclusionID = atoi(idStr);
            tmpList[numUsed].radius = atof(radiusStr);
            tmpList[numUsed].coordinates[X] = atof(coordStr[X]);
            tmpList[numUsed].coordinates[Y] = atof(coordStr[Y]);
            tmpList[numUsed].coordinates[Z] = atof(coordStr[Z]);
            tmpList[numUsed].strainField[0] = atof(strainStr[0]);
            tmpList[numUsed].strainField[1] = atof(strainStr[1]);
            tmpList[numUsed].strainField[2] = atof(strainStr[2]);
            tmpList[numUsed].strainField[3] = atof(strainStr[3]);
            tmpList[numUsed].strainField[4] = atof(strainStr[4]);
            tmpList[numUsed].strainField[5] = atof(strainStr[5]);
         
            maxSize = MAX(maxSize, tmpList[numUsed].radius);

            numUsed++;
        }

        fclose(fp);

/*
 *      Strip off any unused inclusion structures from the inclusion array
 *      before returning it to the caller.
 */
        if (numUsed != numAllocated) {
            tmpList = (Inclusion_t *)realloc(tmpList, numUsed *
                                             sizeof(Inclusion_t));
        }

        *numInclusions = numUsed;
        *inclusionList = tmpList;
        *largestRadius = maxSize;

        return;
}


int main(int argc, char *argv[])
{
        int          i, j, k, x, y, z;
        int          seed, argIndex, totInclusionCnt, overlapping;
        int          cellID, totCellCnt, numCellInclusions;
        int          numCells[3], cellSize[3];
        int          cellIndices[3], nbrIndices[3];
        int          numInclusions = 0, overlapCnt = 0;
        unsigned int overlapCheckCnt = 0;
        real8       radius, largestRadius;
        char         *token;
        Inclusion_t  *inclusion, *inclusionList = (Inclusion_t *)NULL;
        Cell_t       *cellList=0;
        FILE         *fpOut;

/*
 *      Get the command line arguments
 */
        for (argIndex = 1; argIndex < argc; ++argIndex) {
            if (!strcmp(argv[argIndex], "-a")) {
                allowOverlap = 1;
            } else if (!strcmp(argv[argIndex], "-d")) {
                printDebug = 1;
            } else if (!strcmp(argv[argIndex], "-i")) {
                if (argIndex >= (argc - 1)) {
                    Usage(argv[0]);
                }
                inputFile = argv[++argIndex];
            } else if (!strcmp(argv[argIndex], "-h")) {
                PrintHelp(argv[0]);
            } else if (!strcmp(argv[argIndex], "-help")) {
                PrintHelp(argv[0]);
            } else if (!strcmp(argv[argIndex], "-l")) {
                if (argIndex >= (argc - 1)) {
                    Usage(argv[0]);
                }
                argIndex++;
                token = strtok(argv[argIndex], ",");
                simSize[X] = atof(token);
                token = strtok(NULL, ",");
                if (token != (char *)NULL) {
                    simSize[Y] = atof(token);
                    token = strtok(NULL, ",");
                    if (token == (char *)NULL) Usage(argv[0]);
                    simSize[Z] = atof(token);
                } else {
                    simSize[Y] = simSize[X];
                    simSize[Z] = simSize[X];
                }
                simBoundsMin[X] = simSize[X] * -0.5;
                simBoundsMin[Y] = simSize[Y] * -0.5;
                simBoundsMin[Z] = simSize[Z] * -0.5;
                simBoundsMax[X] = simSize[X] *  0.5;
                simBoundsMax[Y] = simSize[Y] *  0.5;
                simBoundsMax[Z] = simSize[Z] *  0.5;
            } else if (!strcmp(argv[argIndex], "-n")) {
                if (argIndex >= (argc - 1)) {
                    Usage(argv[0]);
                }
                numToCreate = atoi(argv[++argIndex]);
            } else if (!strcmp(argv[argIndex], "-o")) {
                if (argIndex >= (argc - 1)) {
                    Usage(argv[0]);
                }
                outputFile = argv[++argIndex];
            } else if (!strcmp(argv[argIndex], "-p")) {
                if (argIndex >= (argc - 1)) {
                    Usage(argv[0]);
                }
                pbcVal = atoi(argv[++argIndex]);
            } else if (!strcmp(argv[argIndex], "-r")) {
                if (argIndex >= (argc - 1)) {
                    Usage(argv[0]);
                }
                argIndex++;
                token = strtok(argv[argIndex], ",");
                radiusMin = atof(token);
                token = strtok(NULL, ",");
                if (token != (char *)NULL) {
                    radiusMax = atof(token);
                } else {
                    radiusMax = radiusMin;
                }
                if (radiusMin > radiusMax) {
                    printf("Error: minimum radius > maximum radius!\n");
                    exit(1);
                }
            } else if (!strcmp(argv[argIndex], "-s")) {
                if (argIndex >= (argc - 1)) {
                    Usage(argv[0]);
                }
                argIndex++;
                token = strtok(argv[argIndex], ",");
                strainField[0] = atof(token);
                for (i = 1; i < 6; i++) {
                    token = strtok(NULL, ",");
                    if (token == (char *)NULL) {
                        Usage(argv[0]);
                    }
                    strainField[i] = atof(token);
                }
            } else {
                Usage(argv[0]);
            }
        }

/*
 *      A couple quick sanity checks... 
 */
        if ((simBoundsMin[X] > simBoundsMax[X]) ||
            (simBoundsMin[Y] > simBoundsMax[Y]) ||
            (simBoundsMin[Z] > simBoundsMax[Z])) {
            fprintf(stderr, "Error: minimum boundaries must be < "
                    "maximum boundaries!\n");
            exit(1);
        }

/*
 *      Now that we've reset defaults based on command-line
 *      args, do a few more initializations.
 */
        for (i = 0; i < 3; i++) {
            simSize[i] = simBoundsMax[i] - simBoundsMin[i];
            simSizeInv[i] = 1.0 / simSize[i];
        }

        largestRadius = radiusMax;

        isFreeSurf[X] = ((pbcVal & 0x01) == 0);
        isFreeSurf[Y] = ((pbcVal & 0x02) == 0);
        isFreeSurf[Z] = ((pbcVal & 0x04) == 0);

/*
 *      If the caller provided an input file of inclusions, read
 *      those inclusions into the local array before creating
 *      any new inclusions.
 */
        if (inputFile != (char *)NULL) {
            ReadInclusionFile(inputFile, &inclusionList, &numInclusions,
                              &largestRadius);
        }

/*
 *      Print out the configuration to be used for the run
 */
        fprintf(stderr, "+++\n");
        if (numInclusions > 0) {
        fprintf(stderr, "+++ Existing inclusions:  %d\n", numInclusions);
            
        }
        fprintf(stderr, "+++ Inclusions to create: %d\n", numToCreate);
        fprintf(stderr, "+++ Radius range:         %.1lf - %.1lf\n",
                radiusMin, radiusMax);
        fprintf(stderr, "+++ Strain field used:    "
                "%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",
                strainField[0], strainField[1], strainField[2],
                strainField[3], strainField[4], strainField[5]);
        fprintf(stderr, "+++ Overlap allowed:      %s\n",
                allowOverlap ? "YES" : "NO");
        fprintf(stderr, "+++ Sim boundaries:        minBound / maxBound\n");
        fprintf(stderr, "+++     X dimension       %.3lf / %.3lf\n", 
                simBoundsMin[X], simBoundsMax[X]);
        fprintf(stderr, "+++     Y dimension       %.3lf / %.3lf\n", 
                simBoundsMin[Y], simBoundsMax[Y]);
        fprintf(stderr, "+++     Z dimension       %.3lf / %.3lf\n", 
                simBoundsMin[Z], simBoundsMax[Z]);
        fprintf(stderr, "+++ Boundary conditions:\n");
        fprintf(stderr, "+++     X dimension       %s\n", 
                isFreeSurf[X] ? "Free surface" : "Periodic");
        fprintf(stderr, "+++     Y dimension       %s\n", 
                isFreeSurf[Y] ? "Free surface" : "Periodic");
        fprintf(stderr, "+++     Z dimension       %s\n", 
                isFreeSurf[Z] ? "Free surface" : "Periodic");
        fprintf(stderr, "+++ Input file:           %s\n",
                inputFile == (char *)NULL ? "(NULL)" : inputFile);
        fprintf(stderr, "+++ Output file:          %s\n",
                outputFile == (char *)NULL ? "(NULL)" : outputFile);
        fprintf(stderr, "+++\n");

/*
 *      If overlapping inclusions are not allowed, determine the dimensions
 *      of the cellular grid to be overlaid on the problem space, and
 *      initialize the cell structures.
 *
 *      Note: Memory usage is tied both to the number of inclusions and
 *            the size of the cellualr grid used.  If memory becomes an
 *            issue, the minimum number of cells can be dropped from
 *            75 to something smaller in the "MIN(numCells[i], 75)"
 *            statement below.  This will save some memory, but comes
 *            at the expense of increased computational cost.
 *  
 */
        if (!allowOverlap) {

            for (i = 0; i < 3; i++) {
                numCells[i] = (int)ceil(simSize[i] / (largestRadius * 2.0));
                numCells[i] = MIN(numCells[i], 75);

                cellSize[i] = simSize[i] / (real8)numCells[i];
            }

            totCellCnt = numCells[X] * numCells[Y] * numCells[Z];
            cellList = (Cell_t *)calloc(1, totCellCnt * sizeof(Cell_t));
        }

/*
 *      Resize the inclusion array so it is large enought to hold all the
 *      inclusions about to be created.
 */
        totInclusionCnt = numInclusions + numToCreate;
        inclusionList = (Inclusion_t *)realloc(inclusionList, totInclusionCnt *
                                               sizeof(Inclusion_t));

/*
 *      Loop creating new inclusions
 */

        real8 IncVolume = 0;

        seed = (int)time((time_t *)NULL);

        for (i = numInclusions; i < totInclusionCnt; i++) {

            inclusion = &inclusionList[i];

/*
 *          Get the radius for the next inclusion.
 */
            radius = radiusMin + (radiusMax-radiusMin) * randm(&seed);

            inclusion->radius = radius;

/*
 *          Get coordinates of inclusion.
 *
 *          Note:  If free surfaces are enabled in a dimension, inclusion must
 *                 be at least the inclusion's radius away from the free
 *                 surface to avoid extending through the surface.  (Don't
 *                 know whether this is really needed or not)
 */
            for (j = 0; j < 3; j++) {
                inclusion->coordinates[j] = (randm(&seed) - 0.5) *
                                      (simSize[j] - (isFreeSurf[j]*2.0*radius));
            }

#if 0 // HACK HACK
      // For inclusions to be clustered 
            for (j = 0; j < 3; j++) {
                inclusion->coordinates[j] = (randm(&seed) - 0.5) *
                                      (200 - (2.0*radius));
            }
#endif

#if 0
/*
 *        HACK HACK
 *        for the inclusions to be on specific planes
 */
            {
               double dot;
               double coord[3],dir[3];

               coord[0] = inclusion->coordinates[0];
               coord[1] = inclusion->coordinates[1];
               coord[2] = inclusion->coordinates[2];
                  

#if 0 // For screw
               if (kodd%2 == 0)
               {
                  // Even : (11-2) plane
                  dir[0] = 1./sqrt(6);
                  dir[1] = 1./sqrt(6);
                  dir[2] =-2./sqrt(6);
                  printf("112 dir");
               }
               else
               {
                  // Odd : (-101) plane
                  dir[0] =-1./sqrt(2);
                  dir[1] = 0;
                  dir[2] = 1./sqrt(2);
                  printf("101 dir");
               }
#else // For edge
               // (1 -1 0)
               dir[0] = 1./sqrt(2);
               dir[2] = 1./sqrt(2);
               dir[1] = 0.0;
#endif
               
               dot = DotProduct(coord,dir);
               
               coord[0] = dot * dir[0];
               coord[1] = dot * dir[1];
               coord[2] = dot * dir[2];
               
               inclusion->coordinates[0] -= coord[0];
               inclusion->coordinates[1] -= coord[1];
               inclusion->coordinates[2] -= coord[2];
               
               kodd++;
            }
#endif




/*
 *          If we are not permitting overlapped inclusions we have
 *          to verify that this inclusion does not intersect any
 *          of the previously created inclusions.
 */
            if (!allowOverlap) {

/*
 *              Determine which cell contains the current inclusion
 */
                for (j = 0; j < 3; j++) {
                    cellIndices[j] = (int)floor((inclusion->coordinates[j] +
                                                 (simSize[j]* 0.5)) /
                                                cellSize[j]);
                }

/*
 *              Given the manner in which the cellular grid is defined,
 *              an inclusion can only overlap another inclusion that
 *              is in the same or an immediately neighboring cell.  So,
 *              loop through the cells of the 27-cell cube encompassing
 *              the inclusion, and check if this inclusion overlaps
 *              any of the inclusions on those cell's inclusion lists.
 *              If any overlap is detected, we discard this inclusion and
 *              create another on to try.
 *
 *              Note:  If periodic boundaries are not enabled in a given
 *                     dimension, any 'neighboring' cells that cross a
 *                     periodic boundary in that dimension, are ignored.
 */
                overlapping = 0;

                for (x = cellIndices[X]-1; x <= cellIndices[X]+1; x++) {
                    if (((x < 0) || (x > numCells[X])) && isFreeSurf[X]) {
                        continue;
                    }
                    for (y = cellIndices[Y]-1; y <= cellIndices[Y]+1; y++) {
                        if (((y < 0) || (y > numCells[Y])) && isFreeSurf[Y]) {
                            continue;
                        }
                        for (z = cellIndices[Z]-1; z <= cellIndices[Z]+1; z++) {
                            if (((z < 0) || (z > numCells[Z])) &&
                                 isFreeSurf[Z]) {
                                continue;
                            }
                            nbrIndices[X] = x;
                            nbrIndices[Y] = y;
                            nbrIndices[Z] = z;

                            cellID = IndicesToCellID(numCells, nbrIndices);
                            numCellInclusions = cellList[cellID].numInclusions;
/*
 *                          Loop through all the inclusions contained
 *                          within this neighbor cell, looking for any
 *                          inclusions that overlap the current inclusion.
 */
                            for (k = 0; k < numCellInclusions; k++) {
                                overlapCheckCnt++;
                                overlapping = CheckOverlap(&inclusionList[i],
                                                           &inclusionList[cellList[cellID].IDs[k]]);
/*
 *                              If the inclusions overlap, there's no need
 *                              to continue, so set the loop indices to jump
 *                              out of the loop over neighboring cells.
 */
                                if (overlapping) {
                                    x += 3;
                                    y += 3;
                                    z += 3;
                                    break;
                                }
                            }
                        }
                    }
                }

/*
 *              If this inclusion was found to overlap another, discard
 *              it, and repeat this iteration of the inclusion generation
 *              loop.  Otherwise, add the inclusion to it's cell's
 *              inclusion list.
 */
                if (overlapping) {
                    i--;
                    overlapCnt++;
                    continue;
                }

                AddToCell(i, cellIndices, numCells, cellList);

            }  /* if (!allowOverlap) */

/*
 *          Finish addiing this inclusion to the list
 */
            inclusion->strainField[0] = strainField[0];
            inclusion->strainField[1] = strainField[1];
            inclusion->strainField[2] = strainField[2];
            inclusion->strainField[3] = strainField[3];
            inclusion->strainField[4] = strainField[4];
            inclusion->strainField[5] = strainField[5];

/*
 *          If we're checking for overlapping inclusions, the utility
 *          can take significantly longer.  In that case, periodically
 *          print out a status indicating how many inclusions have
 *          been generated, and just for kicks, the number of
 *          intersections encountered in the process.
 */
            if ((printDebug) && ((i+1) % 100 == 0)) {
                fprintf(stderr, "inclusions generated: %d "
                        "overlaps encountered: %d "
                        "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                        "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                        "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                        "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
                        i+1, overlapCnt);
            }

/*
 *     Add on volume of that inclusion
 */

            IncVolume += radius * radius * radius;


        }   // loops on number of inclusions


/*
 *    Determine the volume fraction of inclusions
 */
        real8 Volume = (simBoundsMax[0]-simBoundsMin[0])*
                       (simBoundsMax[1]-simBoundsMin[1])*
                       (simBoundsMax[2]-simBoundsMin[2]);
        
        VolumeFraction = 4.0*M_PI*IncVolume / Volume;
        printf("+++ Volume of inclusions = %4.2g Volume = %g Volume fraction = %10.5f \n",
               4.0*M_PI*IncVolume/3.0,Volume,VolumeFraction);
        printf("+++ \n");
        

        if (printDebug) {
            fprintf(stderr, "                                       "
                            "                                       \n");
            if (!allowOverlap) {
                printf("Overlap checks performed: %11u\n"
                       "Overlaps encountered:     %11d\n",
                       overlapCheckCnt, overlapCnt);
            }
        }

/*
 *      If the user specified an output file, open it up.  Otherwise
 *      just set up to write to stdout.
 */
        if (outputFile != (char *)NULL) {
            fpOut = fopen(outputFile, "w");
            if (fpOut == (FILE *)NULL) {
                fprintf(stderr, "Open error on output file %s!\n", outputFile);
                exit(1);
            }
        } else {
            fpOut = stdout;
        }

/*
 *      Before we write the inclusion data, put some comment lines
 *      in the output file
 */
   fprintf(fpOut, "#\n");
   fprintf(fpOut, "# Data for a precipitate consists of the following items separated with white space\n");
   fprintf(fpOut, "#\n");
   fprintf(fpOut, "# Precipitate ID     : 1 integer identifying the precipitate regardless of its position in the simulation or the domain encompassing it. (Should be sequential values)\n");
   fprintf(fpOut, "# Position           : 3 values specifying the XYZ coordinates of the center of the precipitate\n");
   fprintf(fpOut, "# Semi-principal axes: 3 values defining the three semi-principal axis of the ellipsoidal precipitate in units of b\n");
   fprintf(fpOut, "# Rotation matrix    : 2 vectors defining the rotation of the ellipsoidal precipitates\n");
   fprintf(fpOut, "# Strain field       : 6 components of the strain field for the particleStrain field is a symmetric matrix so only six components are specified\n");
   fprintf(fpOut, "#\n");
/*
 *      Loop through all the inclusions and write the data out.  This
 *      will also include all the inclusions read from the input file.
 *
 *      Note:  inclusions will be assigned new sequential ID's based on their
 *             positions in the array, just to insure no duplicate ID's are
 *             used.
 */
        for (i = 0; i < totInclusionCnt; i++) 
        {
           fprintf(fpOut, "%d  %11.3lf %11.3lf %11.3lf  %6.2lf %6.2lf %6.2lf  %lf %lf %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf\n",
                   i+1, /* array index determiens the inclusion ID */
                   inclusionList[i].coordinates[X],
                   inclusionList[i].coordinates[Y],
                   inclusionList[i].coordinates[Z],
                   inclusionList[i].radius, inclusionList[i].radius, inclusionList[i].radius,
                   1.0, 0.0, 0.0, 0.0, 1.0, 0.0,  // for compatibility with recent changes from spherical to ellipsoidal prec.
                   inclusionList[i].strainField[0],
                   inclusionList[i].strainField[1],
                   inclusionList[i].strainField[2],
                   inclusionList[i].strainField[3],
                   inclusionList[i].strainField[4],
                   inclusionList[i].strainField[5]);
        }

        fclose(fpOut);

        exit(0);
        return(0);
}
