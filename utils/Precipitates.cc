// Author:       S. Aubry
//
// Module:       Precipitates.c
//
// Description:  This utility is used for creating a file defining
//               spherical, theta prime and T precipitates in a format that can 
//               be understood by ParaDiS.  

// Attention:    The problem is more symetric if the number of precipitate is
//               a multiple of 15.
//
//
// For a description of the comand line options and default
// configuration values, execute:
//
//     Precipitates -help
//

typedef double real8;

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/param.h>

// Define globals and default values to use if not provided by the user
static int   numTPrec = 0, numThetaPrimePrec = 0, numSphePrec=0;
static int   printDebug  = 0;
static real8 axisMin   = 20.0;
static real8 axisMax   = 50.0;
static real8 VolumeFraction = 0.0;
static real8 simSize[3] = {35000.0, 35000.0, 35000.0};
static real8 simSizeInv[3];
static real8 simBoundsMin[3] = {-17500.0, -17500.0, -17500.0};
static real8 simBoundsMax[3] = { 17500.0,  17500.0,  17500.0};
static real8 strainField[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
static char   *outputFile = (char *)NULL;

// default is Aluminum in angstrom. 
// Lattice spacing is a = 4.050 angstroms. b = a/2 [110]. |b| = a/sqrt(2) = 2.86 angstroms.
static real8 burgMag = 2.86; 


#define X 0
#define Y 1
#define Z 2

// Define structure of precipitates-related info we'll need
typedef struct {
   int    precipitateID;
   real8 axis[3];
   real8 coordinates[3];
   real8 Rotation[2][3];
   real8 strainField[6];
} Precipitates_t;


//
// Define the structure of info associed with each cell
// in the cellular grid which gets used if we need to
// check for overlapping precipitates.

typedef struct {
   int numPrecipitates;
   int arraySize;
   int *IDs;
} Cell_t;


//
//
// Function:       Usage
// Description:    Print a brief synopsis of the acceptable
//                 command line format enad exit.
//
// Arguments:
//     exeName   Name of the executable
//
static void Usage(char *exeName)
{
   fprintf(stderr, "Usage:  %s [-a] [-d]\n", exeName);
   fprintf(stderr, "        [-h|-help]\n");
   fprintf(stderr, "        [-l sim box length\n");
   fprintf(stderr, "        [-nT numTPrec]\n");
   fprintf(stderr, "        [-nTheta ThetaPrimenumPrec]\n");
   fprintf(stderr, "        [-nSphe Spherical]\n");
   fprintf(stderr, "        [-o outputFile] \n");
   fprintf(stderr, "        [-b burgMag in angstroms\n");
   fprintf(stderr, "        [-s s00,s11,s22,s12,s02,s01]\n");

   exit(1);
}


//
//
// Function:       PrintHelp
// Description:    Print a brief synopsis of the command line
//                 format, a description of each option, and
//                 print the values that would be used if
//                 the utility was executed with the specified
//                 options.
//
// Arguments:
//     exeName   Name of the executable
//
static void PrintHelp(char *exeName)
{
   printf("Usage:  %s [-a] [-d]\n", exeName);
   printf("        [-l sim length]\n");
   printf("        [-h|-help]\n");
   printf("        [-l lengthX[,lengthY,lengthZ]]\n");
   printf("        [-nT numTPrec]\n");
   printf("        [-nThetaPrime ThetaPrimenumPrec]\n");
   printf("        [-nSphe SphericalPrec]\n");
   printf("        [-o outputFile] \n");

   printf("\nwhere:\n");
   printf("    -h, -help Prints this help info\n");
   printf("    -l        Specifies the length of the simulation space.  \n"
          "              The Y and Z lengths default to X length if not\n"
          "              specified. This option is mutually exclusive\n"
          "              with the bmin and bmax options\n");
   printf("    -nT           Gives the number of new T type precipitates to create\n");
   printf("    -nThetaPrime  Gives the number of new ThetaPrime precipitates to create\n");
   printf("    -nSphe        Gives the number of new Spherical precipitates to create\n");
   printf("    -o            Gives the name of an output file\n");

   printf("\nDefault values:\n\n");
   printf("    Precipitates to create - type T = %d type ThetaPrime = %d - type Sphe = %d\n", numTPrec, numThetaPrimePrec, numSphePrec);
   printf("    Strain field used:   "
          "%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",
          strainField[0], strainField[1], strainField[2],
          strainField[3], strainField[4], strainField[5]);
   printf("    Sim boundaries:       minBound / maxBound\n");
   printf("        X dimension      %.3lf / %.3lf\n", 
          simBoundsMin[X], simBoundsMax[X]);
   printf("        Y dimension      %.3lf / %.3lf\n", 
          simBoundsMin[Y], simBoundsMax[Y]);
   printf("        Z dimension      %.3lf / %.3lf\n", 
          simBoundsMin[Z], simBoundsMax[Z]);
   printf("    Output file:         %s\n",
          outputFile == (char *)NULL ? "(NULL)" : outputFile);
   printf("    burgMag =  %.3lf angstrom\n", burgMag);
   printf("\n");

   exit(0);
}

// Function:       IndicesToCellID
// Description:    Convert the X, Y and Z indices of a cell into
//                 a unique integer ID.  Indices per dimension
//                 are in the range:
//                    -1 <= index <= number of cells in dimension
//                 Cell indices less than zero or more than the
//                 number of cell in the dimension are adjusted
//                 to account for periodic boundaries.
//
// Arguments:
//     numCells  Three element array of number of cells in X, Y and
//               Z dimensions respectively.
//     indices   Three element array of X, Y and Z indices of a cell
//               repsectively.  To account for periodic boundaries,
//               indices may span the range:
//                  -1 <= indices[N] <= numCells[N]
//
static int IndicesToCellID(int numCells[3], int indices[3])
{
   int x, y, z;
   
   x = indices[X];
   y = indices[Y];
   z = indices[Z];
   
   // Adjust indices that are outside the primary cell grid
   // because of periodic boundaries.
   
   if (x < 0) x = numCells[X] - 1;
   if (y < 0) y = numCells[Y] - 1;
   if (z < 0) z = numCells[Z] - 1;
   
   if (x >= numCells[X]) x = 0;
   if (y >= numCells[Y]) y = 0;
   if (z >= numCells[Z]) z = 0;

   
   // The cellIDs are calculated assuming cell indices vary fastest
   // in the X dimension and slowest in the Z.
 
   return(x + y * numCells[X] + z * numCells[X] * numCells[Y]);
}


//
//
// Function:       AddToCell
// Description:    Add an precipitate to the list of precipitates
//                 encompassed by a given cell.
//
// Arguments:
//     precipitateIndex Index in the precipitate list of the precipitate
//                    to be added to the cell.
//     indices        Three element array of X, Y and Z indices of a cell
//                    respectively.  To account for periodic boundaries,
//                    indices may span the range:
//                       -1 <= indices[N] <= numCells[N]
//     numCells       Three element array of number of cells in X, Y and
//                    Z dimensions respectively.
//     cellList       Array of cell structures: 1 per cell in grid.
//
//
static void AddToCell(int precipitateIndex, int indices[3], int numCells[3], Cell_t *cellList)
{
   int    cellID;
   Cell_t *cell;

   // Determine the ID of the cell to which we'll add the precipitate. 
   cellID = IndicesToCellID(numCells, indices);
   cell = &cellList[cellID];

   // If the cell's precipitate list does not have sufficient space
   // for another precipitate, increase the size of the cell's
   // precipitate list before appending the precipitate to the list.
 
   if (cell->numPrecipitates == cell->arraySize) {
      cell->arraySize += 10;
      cell->IDs = (int *)realloc(cell->IDs, cell->arraySize *
                                 sizeof(int));
   }

   cell->IDs[cell->numPrecipitates] = precipitateIndex;
   cell->numPrecipitates += 1;

   return;
}

//
//
// Function:     randm
// Description:  This is a special function for random number
//               generation on 32-bit machines that do not
//               support long integer multiplication and truncation.
//               The technique used is to do the multiplication and
//               addition in parts, by splitting all integers into a
//               'high' and a 'low' part.  The algorithm is exact,
//               and should give machine-independent results.
//
//               The algorithm implemented is (following D.E. Knuth):
//               seed = seed//1592653589 + 453816691
//               if (seed.lt.0) seed = seed + 1 + 2147483647
//               Note that 1592653589 = 48603*2**15 + 30485
//
// Returns:      real8 value in the range  0.0 <= value <= 1.0
//
//
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


//
// Function:       Getline
// Description:    Read the next non-empty, non-comment line
//                 from the specified file stream and return
//                 it to the caller.
//
// Arguments:
//         string          pointer to location at which to return
//                         the line to the caller
//         len             size (in bytes) of <string>
//         fp              stream from which to read the data
//

static void Getline(char *string, int len, FILE *fp)
{
   char    *s;

   while (1) {

      string[0] = 0;
      if ((s = fgets(string, len, fp)) == (char *)NULL) break;
      string[len-1] = 0;

      //         If the line is a coment line, skip it. If
      //         it has any non-whitespace characters, return it to
      //         the caller, otherwise move to the next line.
 
      if (*s =='#') continue;

      for ( ; *s != 0; s++) {
         if ((*s != ' ') && (*s != '\t') && (*s != '\n')) break;
      }

      if (*s != 0) break;
   }

   return;
}



int main (int argc, char *argv[])
{
   int          i, j;
   int          seed, argIndex, totPrecipitateCnt;
   int          numPrecipitates=0;
   real8       axis[3][15];

   // 3 theta prime precipitates + 12 T precipitates + 1 Spherical precipitate
   real8          Rotation[2][3][15];

   char         *token;
   Precipitates_t  *precipitate, *precipitateList = (Precipitates_t *)NULL;
   FILE         *fpOut;

   static real8 Xdir[4][3] = {
     { 1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0)},
     { 1.0/sqrt(3.0),  1.0/sqrt(3.0), -1.0/sqrt(3.0)},
     { 1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0)},
     {-1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0)} };


   static real8 Ydir[12][3] = {
      // for Xdir[0][*] 
      { 0.0,  1.0/sqrt(2), -1.0/sqrt(2)},
      { 1.0/sqrt(2),  0.0, -1.0/sqrt(2)},
      { 1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0},
      // for Xdir[1][*] 
      { 0.0,  1.0/sqrt(2.0),  1.0/sqrt(2.0)},
      { 1.0/sqrt(2.0),  0.0,  1.0/sqrt(2.0)},
      { 1.0/sqrt(2.0), -1.0/sqrt(2.0),  0.0},
      // for Xdir[2][*] 
      { 0.0,  1.0/sqrt(2.0),  1.0/sqrt(2.0)},
      { 1.0/sqrt(2.0),  0.0, -1.0/sqrt(2.0)},
      { 1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0},
      // for Xdir[3][*] 
      { 0.0,  1.0/sqrt(2.0), -1.0/sqrt(2.0)},
      { 1.0/sqrt(2.0),  0.0,  1.0/sqrt(2.0)},
      { 1.0/sqrt(2.0),  1.0/sqrt(2.0),  0.0}};



   // Get the command line arguments
 
   for (argIndex = 1; argIndex < argc; ++argIndex) 
   {
      if (!strcmp(argv[argIndex], "-d")) 
      {
         // Print debug info
         printDebug = 1;
      } 
      else if (!strcmp(argv[argIndex], "-h")) 
      {
         // Print help
         PrintHelp(argv[0]);
      } 
      else if (!strcmp(argv[argIndex], "-help")) 
      {
         // Print help
         PrintHelp(argv[0]);
      } 
      else if (!strcmp(argv[argIndex], "-b")) 
      {
         // Read in burgMag
         if (argIndex >= (argc - 1)) 
         {
            Usage(argv[0]);
         }
         burgMag = atoi(argv[++argIndex]);
      } 
      else if (!strcmp(argv[argIndex], "-l")) 
      {
         // Read in length of simulation box
         if (argIndex >= (argc - 1)) 
         {
            Usage(argv[0]);
         }
         argIndex++;
         token = strtok(argv[argIndex], ",");
         simSize[X] = atof(token);
         token = strtok(NULL, ",");
         if (token != (char *)NULL) 
         {
            simSize[Y] = atof(token);
            token = strtok(NULL, ",");
            if (token == (char *)NULL) Usage(argv[0]);
            simSize[Z] = atof(token);
         } 
         else 
         {
            simSize[Y] = simSize[X];
            simSize[Z] = simSize[X];
         }
         simBoundsMin[X] = simSize[X] * -0.5;
         simBoundsMin[Y] = simSize[Y] * -0.5;
         simBoundsMin[Z] = simSize[Z] * -0.5;
         simBoundsMax[X] = simSize[X] *  0.5;
         simBoundsMax[Y] = simSize[Y] *  0.5;
         simBoundsMax[Z] = simSize[Z] *  0.5;
      } 
      else if (!strcmp(argv[argIndex], "-nT")) 
      {
         // Read in number of T precipitates
         if (argIndex >= (argc - 1)) {
            Usage(argv[0]);
         }
         numTPrec = atoi(argv[++argIndex]);
      } 
      else if (!strcmp(argv[argIndex], "-nThetaPrime")) 
      {
         // Read in number of theta prime precipitates
         if (argIndex >= (argc - 1)) {
            Usage(argv[0]);
         }
         numThetaPrimePrec = atoi(argv[++argIndex]);
      }
      else if (!strcmp(argv[argIndex], "-nSphe")) 
      {
         // Read in number of spherical precipitates
         if (argIndex >= (argc - 1)) {
            Usage(argv[0]);
         }
         numSphePrec = atoi(argv[++argIndex]);
      } 
      else if (!strcmp(argv[argIndex], "-o")) 
      {
         // Read in name of output file
         if (argIndex >= (argc - 1)) {
            Usage(argv[0]);
         }
         outputFile = argv[++argIndex];
      } 
      else if (!strcmp(argv[argIndex], "-s")) 
      {
         // Read in strain field
         if (argIndex >= (argc - 1)) {
            Usage(argv[0]);
         }
         argIndex++;
         token = strtok(argv[argIndex], ",");
         strainField[0] = atof(token);
         for (i = 1; i < 6; i++) {
            token = strtok(NULL, ",");
            if (token == (char *)NULL) 
            {
               Usage(argv[0]);
            }
            strainField[i] = atof(token);
         }
      } else 
      {
         Usage(argv[0]);
      }
   }

   // A couple quick sanity checks... 
   if ((simBoundsMin[X] > simBoundsMax[X]) ||
       (simBoundsMin[Y] > simBoundsMax[Y]) ||
       (simBoundsMin[Z] > simBoundsMax[Z])) {
      fprintf(stderr, "Error: minimum boundaries must be < "
              "maximum boundaries!\n");
      exit(1);
   }

   if (numTPrec + numThetaPrimePrec + numSphePrec == 0)  
   {
      fprintf(stderr, "Error: no precipitates to create\n");
      exit(1);     
   }

   // Now that we've reset defaults based on command-line
   // args, do a few more initializations.
 
   for (i = 0; i < 3; i++) {
      simSize[i] = simBoundsMax[i] - simBoundsMin[i];
      simSizeInv[i] = 1.0 / simSize[i];
   }

   // Print out the configuration to be used for the run
 
   fprintf(stderr, "+++\n");
   if (numPrecipitates > 0) {
      fprintf(stderr, "+++ Existing precipitates:  %d\n", numPrecipitates);
            
   }
   fprintf(stderr, "+++ Precipitates to create: type T = %d, type ThetaPrime = %d, type Spherical = %d\n",
           numTPrec,numThetaPrimePrec,numSphePrec);
   fprintf(stderr, "+++ burgMag: %3lf angstroms\n", burgMag);
   fprintf(stderr, "+++ Strain field used:    "
           "%.3lf %.3lf %.3lf %.3lf %.3lf %.3lf\n",
           strainField[0], strainField[1], strainField[2],
           strainField[3], strainField[4], strainField[5]);
   fprintf(stderr, "+++ Sim boundaries:        minBound / maxBound\n");
   fprintf(stderr, "+++     X dimension       %.3lf / %.3lf\n", simBoundsMin[X], simBoundsMax[X]);
   fprintf(stderr, "+++     Y dimension       %.3lf / %.3lf\n", simBoundsMin[Y], simBoundsMax[Y]);
   fprintf(stderr, "+++     Z dimension       %.3lf / %.3lf\n", simBoundsMin[Z], simBoundsMax[Z]);
   fprintf(stderr, "+++ Output file:          %s\n", outputFile == (char *)NULL ? "(NULL)" : outputFile);
   fprintf(stderr, "+++\n");

   // Resize the precipitate array so it is large enough to hold all the
   // precipitates about to be created.
 
   int numToCreate = numSphePrec + numTPrec + numThetaPrimePrec;
   totPrecipitateCnt = numPrecipitates + numToCreate;
   precipitateList = (Precipitates_t *)realloc(precipitateList, totPrecipitateCnt *
                                               sizeof(Precipitates_t));

   // To determine the volume fraction of particle,
   // the volume of the current ellipsoid needs to
   // be calculated
   real8 EllipsoidsVolume = 0;

   // Different precipitates
  // 
   // Spherical precipitates : 30nm or 85.8 |b|
   //
   int typeSpherical = 1;
   int precType = -1, nPrectype = typeSpherical;
   real8 rdaxis;
   seed = (int)time((time_t *)NULL);
   for (i = 0; i < numSphePrec; i++) 
   {    
      precipitate = &precipitateList[i];
        
      rdaxis = 85.8;

      precipitate->axis[X] = rdaxis;
      precipitate->axis[Y] = rdaxis;
      precipitate->axis[Z] = rdaxis;
      
      // Get the rotation matrix for the next precipitate.
      precipitate->Rotation[0][X] = 1;
      precipitate->Rotation[0][Y] = 0;
      precipitate->Rotation[0][Z] = 0;
      
      precipitate->Rotation[1][X] = 0;
      precipitate->Rotation[1][Y] = 1;
      precipitate->Rotation[1][Z] = 0;
      
      // Get coordinates of precipitate.
      for (j = 0; j < 3; j++) {
         precipitate->coordinates[j] = (randm(&seed) - 0.5) *
            (simSize[j]);
      }

      // Add on volume of that ellipsoid
      EllipsoidsVolume += precipitate->axis[0]*precipitate->axis[1]*precipitate->axis[2];
            
      // Finish adding this precipitate to the list      
      precipitate->strainField[0] = strainField[0];
      precipitate->strainField[1] = strainField[1];
      precipitate->strainField[2] = strainField[2];
      precipitate->strainField[3] = strainField[3];
      precipitate->strainField[4] = strainField[4];
      precipitate->strainField[5] = strainField[5];
   }



   // 
   // Theta Prime precipitates
   // Length range along [100] is 30-55nm. In units of b it's 104b - 195b. 1b = 2.86 angstroms
   // Length range along [0-10] is 28-82nm. In units of b it's 97b - 286b
   // Theta prime precipitates: [100], [010], [001].
   //
   int typeThetaPrime = 3;
   axisMin =  97;
   axisMax = 286;


   nPrectype = typeThetaPrime;
   real8 diffaxis = (axisMax-axisMin);
   seed = (int)time((time_t *)NULL);
   for (i = numSphePrec; i < numThetaPrimePrec+numSphePrec; i++) 
   {    
      if (++precType+1 > nPrectype) precType = 0;

      precipitate = &precipitateList[i];
        
      rdaxis = axisMin + diffaxis * randm(&seed);
      axis[X][0] = rdaxis;
      axis[Y][0] = rdaxis;
      axis[Z][0] = 3;
      
      rdaxis = axisMin + diffaxis * randm(&seed);
      axis[X][1] = 3;
      axis[Y][1] = rdaxis;
      axis[Z][1] = rdaxis;
      
      rdaxis = axisMin + diffaxis * randm(&seed);
      axis[X][2] = rdaxis;
      axis[Y][2] = 3;
      axis[Z][2] = rdaxis;

      precipitate->axis[X] = axis[X][precType];
      precipitate->axis[Y] = axis[Y][precType];
      precipitate->axis[Z] = axis[Z][precType];
      
      // Get the rotation matrix for the next precipitate.
      precipitate->Rotation[0][X] = 1;
      precipitate->Rotation[0][Y] = 0;
      precipitate->Rotation[0][Z] = 0;
      
      precipitate->Rotation[1][X] = 0;
      precipitate->Rotation[1][Y] = 1;
      precipitate->Rotation[1][Z] = 0;
      
      // Get coordinates of precipitate.
      for (j = 0; j < 3; j++) {
         precipitate->coordinates[j] = (randm(&seed) - 0.5) *
            (simSize[j]);
      }

      // Add on volume of that ellipsoid
      EllipsoidsVolume += precipitate->axis[0]*precipitate->axis[1]*precipitate->axis[2];
            
      // Finish adding this precipitate to the list      
      precipitate->strainField[0] = strainField[0];
      precipitate->strainField[1] = strainField[1];
      precipitate->strainField[2] = strainField[2];
      precipitate->strainField[3] = strainField[3];
      precipitate->strainField[4] = strainField[4];
      precipitate->strainField[5] = strainField[5];
   }


   // 
   // T precipitates
   // Length range along [-111] is 12-55nm. In units of b it's 34.32b - 212.3b
   // Length range along [1-11] is 36-54nm. In units of b it's 103b - 154.5b
   // Rotation[2x3][typeT][T]
   //
   int typeT = 12;
   axisMin = 34;
   axisMax = 213;
   diffaxis = (axisMax-axisMin);
   
   int Tindex = 0;
   for (i = 0; i < 4; i++)
      for (j = 0; j < 3; j++)
      {
         Rotation[0][X][Tindex] = Xdir[i][X];
         Rotation[0][Y][Tindex] = Xdir[i][Y];
         Rotation[0][Z][Tindex] = Xdir[i][Z];
         
         Rotation[1][X][Tindex] = Ydir[3*i + j][X];
         Rotation[1][Y][Tindex] = Ydir[3*i + j][Y];
         Rotation[1][Z][Tindex] = Ydir[3*i + j][Z];  

         Tindex++;
      }


   seed = (int)time((time_t *)NULL);
   precType = -1, nPrectype = typeT;
   // Loop creating T precipitates 
   for (i = numThetaPrimePrec+numSphePrec; i < numTPrec + numThetaPrimePrec+numSphePrec; i++) 
   {    
      if (++precType+1 > nPrectype) precType = 0;

      precipitate = &precipitateList[i];
         
      rdaxis = axisMin + diffaxis * randm(&seed);
      precipitate->axis[X] = rdaxis;
      precipitate->axis[Y] = 8.0;
      precipitate->axis[Z] = rdaxis;
      
      // Get the rotation matrix for the next precipitate.
      precipitate->Rotation[0][X] = Rotation[0][X][precType];
      precipitate->Rotation[0][Y] = Rotation[0][Y][precType];
      precipitate->Rotation[0][Z] = Rotation[0][Z][precType];
      
      precipitate->Rotation[1][X] = Rotation[1][X][precType];
      precipitate->Rotation[1][Y] = Rotation[1][Y][precType];
      precipitate->Rotation[1][Z] = Rotation[1][Z][precType];
      
      // Add on volume of that ellipsoid
      EllipsoidsVolume += precipitate->axis[0]*precipitate->axis[1]*precipitate->axis[2];

      // Get coordinates of precipitate.
      for (j = 0; j < 3; j++) {
         precipitate->coordinates[j] = (randm(&seed) - 0.5) *
            (simSize[j]);
      }
            
      // Finish adding this precipitate to the list      
      precipitate->strainField[0] = strainField[0];
      precipitate->strainField[1] = strainField[1];
      precipitate->strainField[2] = strainField[2];
      precipitate->strainField[3] = strainField[3];
      precipitate->strainField[4] = strainField[4];
      precipitate->strainField[5] = strainField[5];
   }

   real8 Volume = (simBoundsMax[0]-simBoundsMin[0])*
                  (simBoundsMax[1]-simBoundsMin[1])*
                  (simBoundsMax[2]-simBoundsMin[2]);
     
   VolumeFraction = 4.0*M_PI*EllipsoidsVolume/3.0 / Volume;
   printf("+++ Volume of Ellipsoids = %4.2g Volume = %g Volume fraction = %10.5f \n", 4.0*M_PI*EllipsoidsVolume/3.0,Volume,VolumeFraction);
   printf("+++ \n");

   // If the user specified an output file, open it up.  Otherwise
   // just set up to write to stdout.
   
   if (outputFile != (char *)NULL) {
      fpOut = fopen(outputFile, "w");
      if (fpOut == (FILE *)NULL) {
         fprintf(stderr, "Open error on output file %s!\n", outputFile);
         exit(1);
      }
   } else {
      fpOut = stdout;
   }
   
   
   // Before we write the inclusion data, put some comment lines
   // in the output file
   fprintf(fpOut, "#\n");
   fprintf(fpOut, "# Data for a precipitate consists of the following items separated with white space\n");
   fprintf(fpOut, "#\n");
   fprintf(fpOut, "# Precipitate ID     : 1 integer identifying the precipitate regardless of its position in the simulation or the domain encompassing it. (Should be sequential values)\n");
   fprintf(fpOut, "# Position           : 3 values specifying the XYZ coordinates of the center of the precipitate\n");
   fprintf(fpOut, "# Semi-principal axes: 3 values defining the three semi-principal axis of the ellipsoidal precipitate in units of b\n");
   fprintf(fpOut, "# Rotation matrix    : 2 vectors defining the rotation of the ellipsoidal precipitates\n");
   fprintf(fpOut, "# Strain field       : 6 components of the strain field for the particleStrain field is a symmetric matrix so only six components are specified\n");
   fprintf(fpOut, "#\n");
   
   // Loop through all the precipitates and write the data out.  
   //
   // Note:  precipitates will be assigned new sequential ID's based on their
   //        positions in the array, just to insure no duplicate ID's are
   //        used.
  for (i = 0; i < totPrecipitateCnt; i++) {
      fprintf(fpOut, "%d %11.3lf %11.3lf %11.3lf %11.3lf %11.3lf %11.3lf"
              " %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf  %lf %lf %lf %lf %lf %lf\n",
              i+1, // array index determines the precipitate ID 

              precipitateList[i].coordinates[X],
              precipitateList[i].coordinates[Y],
              precipitateList[i].coordinates[Z],

              precipitateList[i].axis[X],
              precipitateList[i].axis[Y],
              precipitateList[i].axis[Z],

              precipitateList[i].Rotation[0][X],
              precipitateList[i].Rotation[0][Y],
              precipitateList[i].Rotation[0][Z],

              precipitateList[i].Rotation[1][X],
              precipitateList[i].Rotation[1][Y],
              precipitateList[i].Rotation[1][Z],
              
              
              precipitateList[i].strainField[0],
              precipitateList[i].strainField[1],
              precipitateList[i].strainField[2],
              precipitateList[i].strainField[3],
              precipitateList[i].strainField[4],
              precipitateList[i].strainField[5]);  
   }
   
   fclose(fpOut);
   
   exit(0);
   return(0);
}
