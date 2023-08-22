
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "BurgVec_Table.h"

//------------------------------------------------------------------------------------------------------------

#define VDELETE(a) if (a) { delete [] a; a=0; }

//------------------------------------------------------------------------------------------------------------

BurgVec_Table::BurgVec_Table(void)
{
   bv_cnt   = 0;
   bv_max   = 0;
   bv_tbl   = 0;
   bv_histo = 0;
}

//------------------------------------------------------------------------------------------------------------

BurgVec_Table::~BurgVec_Table()
{
   VDELETE(bv_tbl  );  bv_cnt=0; bv_max=0;
   VDELETE(bv_histo);
}

// Append()
//
// If the given vector is not already in the table, appends both positive and negative
// burgers vector the the existing table .
//------------------------------------------------------------------------------------------------------------

void BurgVec_Table::Append (const real8  bx, const real8 by, const real8 bz)
{
   VDELETE(bv_histo);

   if ( Index(bx,by,bz) < 0 )
   {
      if ( !bv_tbl || (bv_cnt==bv_max) )
      {
         bv_max = ( (bv_max==0) ? 32 : 2*bv_max );
   
         real8 *bv_tmp = new real8[3*bv_max];
   
         if (bv_tmp && bv_tbl) { memcpy(bv_tmp,bv_tbl,3*bv_cnt*sizeof(real8)); }
   
         VDELETE(bv_tbl);
   
         bv_tbl = bv_tmp;
      }
      
      if (bv_tbl && (bv_cnt<bv_max)) 
      { 
         int j = (3*bv_cnt);

         bv_tbl[j  ] =  bx;
         bv_tbl[j+1] =  by;
         bv_tbl[j+2] =  bz; 

         bv_tbl[j+3] = -bx;
         bv_tbl[j+4] = -by;
         bv_tbl[j+5] = -bz; 

         bv_cnt += 2;
      }
   }
}

void BurgVec_Table::Append (const real8 *bv)
{
   if (bv) { Append(bv[0],bv[1],bv[2]); }
}

//------------------------------------------------------------------------------------------------------------

void BurgVec_Table::Append (const Node_t *nodes, const int ncnt)
{
   VDELETE(bv_histo);

   if (nodes && (ncnt>0))
   {
      for (int i=0; (i<ncnt); i++)
      {
         int narms = nodes[i].numNbrs;
   
         for (int j=0; (j<narms); j++)
         {
            real8 bx = nodes[i].burgX[j];
            real8 by = nodes[i].burgY[j];
            real8 bz = nodes[i].burgZ[j];
   
            Append(bx,by,bz);  
         }
      }
   }
}

void BurgVec_Table::Append (const Segment_t *segs , const int scnt)
{
   VDELETE(bv_histo);

   if (segs && (scnt>0))
   {
      for (int i=0; (i<scnt); i++) 
         Append(segs[i].bv); 
   }
}

// Rebuild()
//
// Will delete the current table and rebuild the table using an array of nodes or segments.
//------------------------------------------------------------------------------------------------------------

void BurgVec_Table::Rebuild (const Node_t *nodes, const int ncnt)
{
   VDELETE(bv_tbl  );  bv_cnt=0; bv_max=0;
   VDELETE(bv_histo);

   Append(nodes,ncnt);
}

void BurgVec_Table::Rebuild (const Segment_t *segs , const int scnt)
{
   VDELETE(bv_tbl  );  bv_cnt=0; bv_max=0;
   VDELETE(bv_histo);

   Append(segs,scnt);
}

// Index()
//
// Returns the index of a particular burgers vector within the table. 
// Returns -1 if not found
// Note - currently a linear search through all burgers vectors.  Definitely needs a better search strategy.
//------------------------------------------------------------------------------------------------------------

int BurgVec_Table::Index (const real8  bx, const real8 by, const real8 bz) const
{
   if (bv_tbl && (bv_cnt>0) )
   {
      for (int i=0,j=0; (i<bv_cnt); ++i, j+=3)
      {
         const real8 eps=1.0e-6;

         if (    (fabs(bv_tbl[j  ]-bx)<eps)
              && (fabs(bv_tbl[j+1]-by)<eps)
              && (fabs(bv_tbl[j+2]-bz)<eps) ) return(i);
      }
   }

   return(-1);
}

int BurgVec_Table::Index  (const real8 *bv) const
{
   if (bv) { return( Index(bv[0],bv[1],bv[2]) ); }

   return(-1);
}

//------------------------------------------------------------------------------------------------------------

real8 *BurgVec_Table::Burgers_Vector(const int indx) const
{
   if ( bv_tbl && (bv_cnt>0) && (0<=indx) && (indx<bv_cnt) )
      return( bv_tbl+3*indx );

   return(0);
}

//------------------------------------------------------------------------------------------------------------

int *BurgVec_Table::Histogram (const Segment_t *segs, const int scnt)
{
   VDELETE(bv_histo);

   bv_histo = ( (bv_cnt>0) ? new int[bv_cnt+1] : 0 );

   if (bv_histo) { memset(bv_histo,0,(bv_cnt+1)*sizeof(int)); }

   if ( bv_histo && bv_tbl && (bv_cnt>0) && segs && (scnt>0) )
   {
      for (int i=0; (i<scnt); i++)
      {
         int indx = Index(segs[i].bv);  
            
         if ( (0<=indx) && (indx<bv_cnt) )
            bv_histo[indx]++; 
      }
   }

   return(bv_histo);
}

//------------------------------------------------------------------------------------------------------------

void BurgVec_Table::Print (FILE *fd) const
{
   if (fd && bv_tbl && (bv_cnt>0))
   {
      fprintf(fd,"# burgers vector table \n");
      fprintf(fd,"# bv_cnt=%d\n", bv_cnt);

      for (int i=0,j=0; (i<bv_cnt); i++, j+=3)
         fprintf(fd," %2d : %6d %12.8lf %12.8lf %12.8lf\n", i, (bv_histo ? bv_histo[i] : 0), bv_tbl[j], bv_tbl[j+1], bv_tbl[j+2] );
   }
}

