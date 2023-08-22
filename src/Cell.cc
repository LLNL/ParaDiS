#include <stdio.h>
#include <math.h>
#include <float.h>

#include "Typedefs.h"
#include "Cell.h"

// _init()
//
// Initializes a node structure.
//---------------------------------------------------------------------------------------------------------

Cell_t::Cell_t(void)
{
   cellID         = 0;
   inUse          = 0;

   next           = 0;

   nodeQ          = 0;
   nodeCount      = 0;

   inclusionQ     = 0;
   inclusionCount = 0;

   nbrList        = 0;
   nbrCount       = 0;

   domains        = 0;
   domCount       = 0;

   baseIdx        = 0;

   xShift         = 0.0;
   yShift         = 0.0;
   zShift         = 0.0;

   xIndex         = 0;
   yIndex         = 0;
   zIndex         = 0;

   center[0]      = DBL_MAX;
   center[1]      = DBL_MAX;
   center[2]      = DBL_MAX;

#ifdef _OPENMP
   omp_init_lock(&cellLock);
#endif
}

Cell_t::~Cell_t()
{
   // (to be resolved - bmw)
}

void Cell_t::Print (FILE *fd) const
{
   if (fd)
   {
      fprintf(fd,"%d %d base=%d ncnt=%d icnt=%d\n", cellID, inUse, baseIdx, nodeCount, inclusionCount);

      if (center[0]<1.0e+30)  // (cell centers are initialized with DBL_MAX)
         fprintf(fd,"%3d center=[%14.6lf %14.6lf %14.6lf]\n", cellID, center[0], center[1], center[2] );

      fprintf(fd,"   shift =[%14.6lf %14.6lf %14.6lf]\n", xShift, yShift, zShift );
      fprintf(fd,"   index =[%d %d %d]\n"               , xIndex, yIndex, zIndex );

      if (nbrList)
      {
         fprintf(fd,"   ncnt=%d nbrList=[", nbrCount );
         for (int i=0; (i<nbrCount); i++) { fprintf(fd,"%3d ", nbrList[i] ); }
         fprintf(fd,"]\n");
      }

      if (domains)
      {
         fprintf(fd,"   dcnt=%d domains=[", domCount );
         for (int i=0; (i<domCount); i++) { fprintf(fd,"%3d ", domains[i] ); }
         fprintf(fd,"]\n");
      }
   }
}

//---------------------------------------------------------------------------
//  Function:     Cell_Table_Print
//  Description:  Simple diagnostic print of the current cell table
//
//  Arguments:
//      fd        points to a currently open file descriptor 
//      ctbl      points to a cell table (array of pointers to cells)
//      n         length of the cell table 
//---------------------------------------------------------------------------

void Cell_Table_Print (FILE *fd, const Cell_t **ctbl, const int n)
{
   if (fd && ctbl)
   {
      for (int i=0; (i<n); i++)
         if (ctbl[i]) ctbl[i]->Print(fd);
   }
}

