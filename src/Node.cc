//---------------------------------------------------------------------------
// Module:  Node.cc - implementations of various Node-related functions.
//
// Includes functions:
//    Node_New ()          : creates a new node object
//    Node_Reset()         : resets and existing node
//    Node_Recycle()       : releases and recyles node-related dynamic memory
//    Node_Allocate_Arms() : allocates and initializes node-related dynamic memory
//    Node_Copy()          : deep copy of static and dynamic elements of a node
//    Node_Print()         : prints elements of node in human-readable text
//    Node_Block_New()     : allocates and initializes a new node block
//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"
#include "Home.h"
#include "ParadisThread.h"
#include "Tag.h"
#include "Node.h"
#include "CPU_Timer.h"

#define MFREE(a)          if (a) { free(a); a=0; }
#define MZERO(a,bytes)    if (a) { memset(a,0,bytes); }
#define MCOPY(a,b,bytes)  if (a && b) { memcpy(a,b,bytes); }

// _init()
//
// Initializes a node structure.
//---------------------------------------------------------------------------------------------------------

static void _init(Node_t & n)
{
   n.flags     = 0;
        
   n.x  = n.y  = n.z  = 0.0;
   n.fX = n.fY = n.fZ = 0.0;
   n.vX = n.vY = n.vZ = 0.0;
        
   n.oldx   = n.oldy   = n.oldz   = 0.0;
   n.oldvX  = n.oldvY  = n.oldvZ  = 0.0;
   n.currvX = n.currvY = n.currvZ = 0.0;
 
   n.myTag   = Tag_t();
   n.myIndex = 0;

   n.numNbrs  = 0;                     
   n.nbrTag   = 0;                      
   n.nbrIndex = 0;

   n.armfx = n.armfy = n.armfz = 0;
   n.burgX = n.burgY = n.burgZ = 0;
   n.nx    = n.ny    = n.nz    = 0;

   n.sigbLoc       = 0;                     
   n.sigbRem       = 0;                     

   n.armCoordIndex = 0;               
   n.armNodeIndex  = 0;               

   n.constraint    = 0;                  

   n.cellIdx       = 0;                        
   n.cell2Idx      = 0;                       
   n.cell2QentIdx  = 0;                   

   n.native        = 0;                         

   n.next          = 0;                       
   n.nextInCell    = 0;                 

   n.tmark         = 0;                           

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
   n.multiNodeLife = 0;
#endif

#ifdef _ARLFEM
   n.surfaceNorm[0] = 0.0;
   n.surfaceNorm[1] = 0.0;
   n.surfaceNorm[2] = 0.0;
#endif

#ifdef _OPENMP
   INIT_LOCK(&n.nodeLock);
#endif
}

// Constructors and destructor.
//---------------------------------------------------------------------------------------------------------

Node_t::Node_t(void)                                  { _init(*this);                                  }
Node_t::Node_t(const Tag_t & tag)                     { _init(*this); myTag=tag;                       }
Node_t::Node_t(const Tag_t & tag, const int narms)    { _init(*this); myTag=tag; Allocate_Arms(narms); }
Node_t::Node_t(const Node_t & n)                      { _init(*this); *this=n;                         }

Node_t::~Node_t() { Recycle(); }

// assignment/copy operator
//---------------------------------------------------------------------------------------------------------

const Node_t & Node_t::operator =  (const Node_t & n)
{
   if (this==&n) { return(*this); }

   flags         = n.flags        ;
        
   x             = n.x            ;
   y             = n.y            ;
   z             = n.z            ;
   
   fX            = n.fX           ;
   fY            = n.fY           ;
   fZ            = n.fZ           ;
   
   vX            = n.vX           ;
   vY            = n.vY           ;
   vZ            = n.vZ           ;
   
   oldx          = n.oldx         ;
   oldy          = n.oldy         ;
   oldz          = n.oldz         ;
   
   oldvX         = n.oldvX        ;
   oldvY         = n.oldvY        ;
   oldvZ         = n.oldvZ        ;
   
   currvX        = n.currvX       ;
   currvY        = n.currvY       ;
   currvZ        = n.currvZ       ;
 
   myTag         = n.myTag        ;
   myIndex       = n.myIndex      ;

   Allocate_Arms(n.numNbrs);

   MCOPY(nbrTag        ,n.nbrTag        ,   numNbrs*sizeof(Tag_t) );
   MCOPY(nbrIndex      ,n.nbrIndex      ,   numNbrs*sizeof(int  ) );
   MCOPY(armfx         ,n.armfx         ,   numNbrs*sizeof(real8) );
   MCOPY(armfy         ,n.armfy         ,   numNbrs*sizeof(real8) );
   MCOPY(armfz         ,n.armfz         ,   numNbrs*sizeof(real8) );
   MCOPY(burgX         ,n.burgX         ,   numNbrs*sizeof(real8) );
   MCOPY(burgY         ,n.burgY         ,   numNbrs*sizeof(real8) );
   MCOPY(burgZ         ,n.burgZ         ,   numNbrs*sizeof(real8) );
   MCOPY(nx            ,n.nx            ,   numNbrs*sizeof(real8) );
   MCOPY(ny            ,n.ny            ,   numNbrs*sizeof(real8) );
   MCOPY(nz            ,n.nz            ,   numNbrs*sizeof(real8) );
   MCOPY(sigbLoc       ,n.sigbLoc       , 3*numNbrs*sizeof(real8) );
   MCOPY(sigbRem       ,n.sigbRem       , 3*numNbrs*sizeof(real8) );
   MCOPY(armCoordIndex ,n.armCoordIndex ,   numNbrs*sizeof(int  ) );
   MCOPY(armNodeIndex  ,n.armNodeIndex  ,   numNbrs*sizeof(int  ) );

   constraint    = n.constraint   ;

   cellIdx       = n.cellIdx      ;
   cell2Idx      = n.cell2Idx     ;
   cell2QentIdx  = n.cell2QentIdx ;

   native        = n.native       ;

   next          = n.next         ;
   nextInCell    = n.nextInCell   ;

   tmark         = n.tmark        ;

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
   multiNodeLife = n.multiNodeLife;
#endif

#ifdef _ARLFEM
   surfaceNorm[0] = n.surfaceNorm[0];
   surfaceNorm[1] = n.surfaceNorm[1];
   surfaceNorm[2] = n.surfaceNorm[2];
#endif

   return(*this);
}

// Comparison operations
//
// Two nodes are considered less than or greater than by comparing tags only. 
//---------------------------------------------------------------------------------------------------------

int Node_t::operator <  (const Node_t & n) const { return(myTag<n.myTag); }
int Node_t::operator >  (const Node_t & n) const { return(myTag>n.myTag); }
int Node_t::operator != (const Node_t & n) const { return( (*this==n) ? 0 : 1 ); }

// Node equality is based on a subset of the the node state parameters.
// Essentially scan through the relevant fields and return non-equal if any disagree.
//---------------------------------------------------------------------------------------------------------

int Node_t::operator == (const Node_t & n) const 
{ 
   if (this==&n) return(1);               // deal with (a==a)

   if (myTag!=n.myTag) return(0);

   real8 eps = 1.0e-4;

   if (fabs(x-n.x)>eps) return(0);
   if (fabs(y-n.y)>eps) return(0);
   if (fabs(z-n.z)>eps) return(0);
   
   if (flags!=n.flags) return(0);
           
   if (numNbrs!=n.numNbrs) return(0);

   for (int i=0; (i<numNbrs); ++i)
   {
      if (nbrTag  [i]!=n.nbrTag  [i]) return(0);
      if (nbrIndex[i]!=n.nbrIndex[i]) return(0);

      eps = 1.0e-1;
      if (fabs(armfx[i]-n.armfx[i])>eps) return(0);
      if (fabs(armfy[i]-n.armfy[i])>eps) return(0);
      if (fabs(armfz[i]-n.armfz[i])>eps) return(0);

      eps = 1.0e-6;
      if (fabs(burgX[i]-n.burgX[i])>eps) return(0);
      if (fabs(burgY[i]-n.burgY[i])>eps) return(0);
      if (fabs(burgZ[i]-n.burgZ[i])>eps) return(0);

      eps = 1.0e-6;
      if (fabs(nx   [i]-n.nx   [i])>eps) return(0);
      if (fabs(ny   [i]-n.ny   [i])>eps) return(0);
      if (fabs(nz   [i]-n.nz   [i])>eps) return(0);

   } 

   return(1); 
}

//---------------------------------------------------------------------------------------------------------

void Node_t::Allocate_Arms (const int narms)
{
   if (numNbrs!=narms)
   {
      Recycle();

      if (narms>0)
      {
         nbrTag        = (Tag_t *) malloc(  narms*sizeof(Tag_t));
         nbrIndex      = (int   *) malloc(  narms*sizeof(int  ));
         burgX         = (real8 *) malloc(  narms*sizeof(real8));
         burgY         = (real8 *) malloc(  narms*sizeof(real8));
         burgZ         = (real8 *) malloc(  narms*sizeof(real8));
         armfx         = (real8 *) malloc(  narms*sizeof(real8));
         armfy         = (real8 *) malloc(  narms*sizeof(real8));
         armfz         = (real8 *) malloc(  narms*sizeof(real8));
         nx            = (real8 *) malloc(  narms*sizeof(real8));
         ny            = (real8 *) malloc(  narms*sizeof(real8));
         nz            = (real8 *) malloc(  narms*sizeof(real8));
         sigbLoc       = (real8 *) malloc(3*narms*sizeof(real8)); 
         sigbRem       = (real8 *) malloc(3*narms*sizeof(real8)); 
         armCoordIndex = (int   *) malloc(  narms*sizeof(int  ));
         armNodeIndex  = (int   *) malloc(  narms*sizeof(int  ));
      }

      numNbrs = narms;
   }

   MZERO(nbrTag       ,   numNbrs*sizeof(Tag_t) );
   MZERO(nbrIndex     ,   numNbrs*sizeof(int  ) );
   MZERO(burgX        ,   numNbrs*sizeof(real8) );
   MZERO(burgY        ,   numNbrs*sizeof(real8) );
   MZERO(burgZ        ,   numNbrs*sizeof(real8) );
   MZERO(armfx        ,   numNbrs*sizeof(real8) );
   MZERO(armfy        ,   numNbrs*sizeof(real8) );
   MZERO(armfz        ,   numNbrs*sizeof(real8) );
   MZERO(nx           ,   numNbrs*sizeof(real8) );
   MZERO(ny           ,   numNbrs*sizeof(real8) );
   MZERO(nz           ,   numNbrs*sizeof(real8) );
   MZERO(sigbLoc      , 3*numNbrs*sizeof(real8) );
   MZERO(sigbRem      , 3*numNbrs*sizeof(real8) );
   MZERO(armCoordIndex,   numNbrs*sizeof(int  ) );
   MZERO(armNodeIndex ,   numNbrs*sizeof(int  ) );
}

void Node_t::Free_Arms (void)
{
   MFREE(nbrTag       );
   MFREE(nbrIndex     );
   MFREE(armfx        );
   MFREE(armfy        );
   MFREE(armfz        );
   MFREE(burgX        );
   MFREE(burgY        );
   MFREE(burgZ        );
   MFREE(nx           );
   MFREE(ny           );
   MFREE(nz           );
   MFREE(sigbLoc      );
   MFREE(sigbRem      );
   MFREE(armCoordIndex);
   MFREE(armNodeIndex );

   numNbrs = 0;
}

void Node_t::Recycle  (void)
{
   MFREE(nbrTag       );
   MFREE(nbrIndex     );
   MFREE(armfx        );
   MFREE(armfy        );
   MFREE(armfz        );
   MFREE(burgX        );
   MFREE(burgY        );
   MFREE(burgZ        );
   MFREE(nx           );
   MFREE(ny           );
   MFREE(nz           );
   MFREE(sigbLoc      );
   MFREE(sigbRem      );
   MFREE(armCoordIndex);
   MFREE(armNodeIndex );

   numNbrs = 0;
}

size_t Node_t::Bytes (void) const
{
   int n = numNbrs;

   size_t bytes  = sizeof(Node_t);
          bytes += ( nbrTag        ?   n*sizeof(Tag_t) : 0 );
          bytes += ( nbrIndex      ?   n*sizeof(int  ) : 0 );

          bytes += ( armfx         ?   n*sizeof(real8) : 0 );
          bytes += ( armfy         ?   n*sizeof(real8) : 0 );
          bytes += ( armfz         ?   n*sizeof(real8) : 0 );

          bytes += ( burgX         ?   n*sizeof(real8) : 0 );
          bytes += ( burgY         ?   n*sizeof(real8) : 0 );
          bytes += ( burgZ         ?   n*sizeof(real8) : 0 );

          bytes += ( nx            ?   n*sizeof(real8) : 0 );
          bytes += ( ny            ?   n*sizeof(real8) : 0 );
          bytes += ( nz            ?   n*sizeof(real8) : 0 );

          bytes += ( sigbLoc       ? 3*n*sizeof(real8) : 0 );
          bytes += ( sigbRem       ? 3*n*sizeof(real8) : 0 );
   
          bytes += ( armCoordIndex ?   n*sizeof(int  ) : 0 );
          bytes += ( armNodeIndex  ?   n*sizeof(int  ) : 0 );

   return(bytes);
}

// IsA_M_Node()
//
// An M-node is a node with more than two arms where the arms all have different burgers vectors.
//---------------------------------------------------------------------------------------------------------

int Node_t::IsA_M_Node (void) const
{
   if ( (numNbrs>2) && burgX && burgY && burgZ )
   {
      for (int i=0; (i<(numNbrs-1)); i++)
      {
         for (int j=(i+1); (j<numNbrs); j++)
         {
            if (    (fabs(burgX[i]-burgX[j])<1.0e-6) 
                 && (fabs(burgY[i]-burgY[j])<1.0e-6) 
                 && (fabs(burgZ[i]-burgZ[j])<1.0e-6) ) return(0);
         }
      }

      return(1);
   }

   return(0);
}

// IsA_N_Node()
//
// An N-node is a node with four arms where the burgers vectors represent two crossing lines.
//---------------------------------------------------------------------------------------------------------

int Node_t::IsA_N_Node (void) const
{
   if ( (numNbrs==4) && burgX && burgY && burgZ )
   {
      if (    (fabs(burgX[0]+burgX[1])<1.0e-6) 
           && (fabs(burgY[0]+burgY[1])<1.0e-6) 
           && (fabs(burgZ[0]+burgZ[1])<1.0e-6)
           && (fabs(burgX[2]+burgX[3])<1.0e-6) 
           && (fabs(burgY[2]+burgY[3])<1.0e-6) 
           && (fabs(burgZ[2]+burgZ[3])<1.0e-6) )  return(1);     // (b0==-b1) && (b2==-b3)

      if (    (fabs(burgX[0]+burgX[2])<1.0e-6) 
           && (fabs(burgY[0]+burgY[2])<1.0e-6) 
           && (fabs(burgZ[0]+burgZ[2])<1.0e-6)
           && (fabs(burgX[1]+burgX[3])<1.0e-6) 
           && (fabs(burgY[1]+burgY[3])<1.0e-6) 
           && (fabs(burgZ[1]+burgZ[3])<1.0e-6) )  return(1);     // (b0==-b2) && (b1==-b3)

      if (    (fabs(burgX[0]+burgX[3])<1.0e-6) 
           && (fabs(burgY[0]+burgY[3])<1.0e-6) 
           && (fabs(burgZ[0]+burgZ[3])<1.0e-6)
           && (fabs(burgX[1]+burgX[2])<1.0e-6) 
           && (fabs(burgY[1]+burgY[2])<1.0e-6) 
           && (fabs(burgZ[1]+burgZ[2])<1.0e-6) )  return(1);     // (b0==-b3) && (b1==-b2)
   }

   return(0);
}

// GPU_Bytes()
//
// Returs the number of bytes needed to store a gpu-packed node.
//---------------------------------------------------------------------------------------------------------

size_t Node_t::GPU_Bytes (void)
{
   int n = numNbrs;

   size_t bytes  =     sizeof(int  );  // tag domain id
          bytes +=     sizeof(int  );  // tag index   
          bytes +=     sizeof(int  );  // narms       
          bytes +=     sizeof(int  );  // native      
          bytes +=     sizeof(real8);  // position x
          bytes +=     sizeof(real8);  // position y
          bytes +=     sizeof(real8);  // position z
          bytes += 2*n*sizeof(int  );  // arm tags[]
          bytes +=   n*sizeof(int  );  // arm indices[]
          bytes +=   n*sizeof(real8);  // armfx[]
          bytes +=   n*sizeof(real8);  // armfy[]
          bytes +=   n*sizeof(real8);  // armfz[]

   return(bytes);
}

// GPU_Pack()
//
// Packs a gpu-packed node into a buffer. Returns the next available position in the buffer.
//---------------------------------------------------------------------------------------------------------

unsigned char *Node_t::GPU_Pack (unsigned char *q)
{
   int n = numNbrs;

   if (q)
   {
      *((int   *) q) = myTag.domainID;  q+=sizeof(int  );
      *((int   *) q) = myTag.index   ;  q+=sizeof(int  );
      *((int   *) q) = numNbrs       ;  q+=sizeof(int  );
      *((int   *) q) = native        ;  q+=sizeof(int  );
      *((real8 *) q) = x             ;  q+=sizeof(real8); 
      *((real8 *) q) = y             ;  q+=sizeof(real8);
      *((real8 *) q) = z             ;  q+=sizeof(real8);
   
      if (n>0)
      { 
         for (int i=0; (i<n); ++i)
         {
            *((int *) q) = nbrTag[i].domainID;  q+=sizeof(int);
            *((int *) q) = nbrTag[i].index   ;  q+=sizeof(int);
         }

         MCOPY(q,nbrIndex,(n*sizeof(int  ))); q+=(n*sizeof(int  ));
         MCOPY(q,armfx   ,(n*sizeof(real8))); q+=(n*sizeof(real8));
         MCOPY(q,armfy   ,(n*sizeof(real8))); q+=(n*sizeof(real8));
         MCOPY(q,armfz   ,(n*sizeof(real8))); q+=(n*sizeof(real8));
      }
   }

   return(q);
}

// GPU_Unpack()
//
// Unpacks a gpu-packed node from a buffer. Returns the next available position in the buffer.
//---------------------------------------------------------------------------------------------------------

unsigned char *Node_t::GPU_Unpack (unsigned char *p)
{
   int n = numNbrs;

   if (p)
   {
      myTag.domainID = *((int   *) p);  p+=sizeof(int  );
      myTag.index    = *((int   *) p);  p+=sizeof(int  );
      numNbrs        = *((int   *) p);  p+=sizeof(int  );
      native         = *((int   *) p);  p+=sizeof(int  );
      x              = *((real8 *) p);  p+=sizeof(real8);
      y              = *((real8 *) p);  p+=sizeof(real8);
      z              = *((real8 *) p);  p+=sizeof(real8);

      if (n>0)
      { 
         for (int i=0; (i<n); ++i)
         {
            nbrTag[i].domainID = *((int *) p);  p+=sizeof(int);
            nbrTag[i].index    = *((int *) p);  p+=sizeof(int);
         }
            
         MCOPY(nbrIndex,p,(n*sizeof(int  ))); p+=(n*sizeof(int  ));
         MCOPY(armfx   ,p,(n*sizeof(real8))); p+=(n*sizeof(real8));
         MCOPY(armfy   ,p,(n*sizeof(real8))); p+=(n*sizeof(real8));
         MCOPY(armfz   ,p,(n*sizeof(real8))); p+=(n*sizeof(real8));
      }
   }

   return(p);
}

// Validate()
//
// A node with 2 or more arms is considered valid if the burgers vectors of all the arms sum to zero.
//---------------------------------------------------------------------------------------------------------

int Node_t::Validate (void) const
{
   int narms = numNbrs;
 
   if ( (narms>1) && (constraint!=PINNED_NODE) )
   {
      real8 sum[3] = { 0.0, 0.0, 0.0 };

      for (int i=0; (i<narms); ++i)
      {
         sum[0] += burgX[i];
         sum[1] += burgY[i];
         sum[2] += burgZ[i];
      }

      real8 eps=1.0e-6;

      if (fabs(sum[0])>eps) return(0);
      if (fabs(sum[1])>eps) return(0);
      if (fabs(sum[2])>eps) return(0);
   }

   return(1);
}

// Inside()
//
// Will return if the node is physically located inside the domain provided.
//---------------------------------------------------------------------------------------------------------

int Node_t::Inside 
(
   const real8 x0, const real8 x1, 
   const real8 y0, const real8 y1, 
   const real8 z0, const real8 z1 ) const
{

   if ( (x<x0) || (x1<=x) ) return(0);
   if ( (y<y0) || (y1<=y) ) return(0);
   if ( (z<z0) || (z1<=z) ) return(0);

   return(1);
}

// Arm_Index()
//
// Scans through the node arm tags and returns the index of the arm associated with the tag.
// Returns -1 if tag not found.  Note - linear search - not particularly efficient.
//---------------------------------------------------------------------------------------------------------

int Node_t::Arm_Index (const Tag_t & tag) const
{
   if (nbrTag && (numNbrs>0))
   {
      for (int i=0; (i<numNbrs); ++i)
         if ( tag==nbrTag[i] ) return(i);
   }

   return(-1);
}

// Arm_Index()
//
// Scans through the node arms and returns the index of the arm that matches the burgers
// vector whose elements match (bx,by,bz);
//---------------------------------------------------------------------------------------------------------

int Node_t::Arm_Index (const real8 bx, const real8 by, const real8 bz) const
{
   if ( (numNbrs>0) && burgX && burgY && burgZ )
   {
      for (int i=0; (i<numNbrs); i++)
      {
         if (    (fabs(bx-burgX[i])<1.0e-6) 
              && (fabs(by-burgY[i])<1.0e-6) 
              && (fabs(bz-burgZ[i])<1.0e-6) ) return(i);
      }
   }

   return(-1);
}

void Node_t::Delete_Arm (const Tag_t & tag)
{
   if (numNbrs>0)
   {
      int indx = Arm_Index(tag);
   
      if (indx>=0)
      {
         for (int j=indx; (j<(numNbrs-1)); j++)
         {
            nbrTag       [j] = nbrTag       [j+1];
            nbrIndex     [j] = nbrIndex     [j+1];
            armfx        [j] = armfx        [j+1];
            armfy        [j] = armfy        [j+1];
            armfz        [j] = armfz        [j+1];
            burgX        [j] = burgX        [j+1];
            burgY        [j] = burgY        [j+1];
            burgZ        [j] = burgZ        [j+1];
            nx           [j] = nx           [j+1];
            ny           [j] = ny           [j+1];
            nz           [j] = nz           [j+1];
            sigbLoc      [j] = sigbLoc      [j+1];
            sigbRem      [j] = sigbRem      [j+1];
            armCoordIndex[j] = armCoordIndex[j+1];
            armNodeIndex [j] = armNodeIndex [j+1];
         }

         numNbrs--;
      }
   }

   if (numNbrs==0) { Free_Arms(); }
}

// Reset_Forces()
//
// Will reset the accumulated forces AND the arm forces to zero.
// Note - this method is NOT thread safe - any thread can reset the 
// forces to zero and override any other thread.
//---------------------------------------------------------------------------------------------------------

void Node_t::Forces_Reset (void) 
{ 
   LOCK(&nodeLock);

   fX = fY = fZ = 0.0; 

   MZERO(armfx,numNbrs*sizeof(real8));
   MZERO(armfy,numNbrs*sizeof(real8));
   MZERO(armfz,numNbrs*sizeof(real8));

   UNLOCK(&nodeLock);
}

// Add_Arm_Force()
//
// Adds the forces to a given arm.
//---------------------------------------------------------------------------------------------------------

void Node_t::Add_Arm_Force (const Tag_t & tag, const real8 *fv) 
{
   int i = Arm_Index(tag);

   if (i>=0)
   {
      LOCK(&nodeLock);
   
      armfx[i] += fv[0];
      armfy[i] += fv[1];
      armfz[i] += fv[2];
   
      UNLOCK(&nodeLock);
   }
}

// Add_Arm_Forces()
//
// This method will sum the individiual arm force controbutions into the main node force vector
//---------------------------------------------------------------------------------------------------------

void Node_t::Add_Arm_Forces (void) 
{ 
   LOCK(&nodeLock);

   fX = fY = fZ = 0.0; 

   if (numNbrs>0)
   {
      for (int i=0; (i<numNbrs); ++i)
      {
         fX += armfx[i];
         fY += armfy[i];
         fZ += armfz[i];
      }
   }

   UNLOCK(&nodeLock);
}

// Arm_Sum()
//
// Returns the sum of the arm forces
//---------------------------------------------------------------------------------------------------------

void Node_t::Arm_Sum (real8 *f) const
{
   if (f)
   {
      f[0]=f[1]=f[2]=0.0;

      if (armfx && armfy && armfz)
      {
         for (int i=0; (i<numNbrs); i++)
         {
            f[0] += armfx[i];
            f[1] += armfy[i];
            f[2] += armfz[i];
         }
      }
   }
}

// Advance()
//
// Will update the node's position given a timestep (dt). This method will average
// the current and previous node velocities for a given timestep.
//
// Additional parameters are provided to adjust for periodic boundary conditions (PBC)...
//
//    cx,cy,cz - the center of the simulation box <xyz>
//    lx,ly,lz - the size   of the simulation box <xyz>
//    sx,sy,sz - PBC adjustment reciprocals       <xyz>
//                 If PBC for a particular axis is active, the PBC reciprocal will be 
//                 the reciprocal of the simulation size, otherwise zero.
//                 e.g  sx = ( pbcx && (lx>0.0) ? (1.0/lx) : 0.0 );
//                      sy = ( pbcy && (ly>0.0) ? (1.0/ly) : 0.0 );
//                      sz = ( pbcz && (lz>0.0) ? (1.0/lz) : 0.0 );
//                 By setting  sx,sy,sz to the reciprocals, the rint() correction will be applied.
//---------------------------------------------------------------------------------------------------------

void Node_t::Advance 
(
   const real8 dt,                                   ///< dt = next time step
   const real8 cx, const real8 cy, const real8 cz,   ///< simulation center <xyz>
   const real8 lx, const real8 ly, const real8 lz,   ///< simulation size   <xyz>
   const real8 sx, const real8 sy, const real8 sz    ///< pbc adjustment reciprocals <xyz>
)
{
   // If we don't have a value for the previous velocity, assume
   // previous velocity was same as current velocity.

   if (    (oldvX == 0.0) 
        && (oldvY == 0.0)
        && (oldvZ == 0.0) )
   {
      oldvX = currvX;
      oldvY = currvY;
      oldvZ = currvZ;
   }

   // update the position...

   x = oldx + 0.5 * (currvX + oldvX) * dt;
   y = oldy + 0.5 * (currvY + oldvY) * dt;
   z = oldz + 0.5 * (currvZ + oldvZ) * dt;

   // adjust the updated position for periodic boundary...

   x -= rint((x-cx)*sx)*lx;
   y -= rint((y-cy)*sy)*ly;
   z -= rint((z-cz)*sz)*lz;
}

// Nodes_Advance()
//
// Advances an array of nodes or node pointers.
//---------------------------------------------------------------------------------------------------------

void Nodes_Advance 
(
         Node_t  *narr,                                 ///< array of nodes
   const int      ncnt,                                 ///< number of nodes
   const real8    dt,                                   ///< dt = next time step
   const real8    cx, const real8 cy, const real8 cz,   ///< simulation center <xyz>
   const real8    lx, const real8 ly, const real8 lz,   ///< simulation size   <xyz>
   const real8    sx, const real8 sy, const real8 sz    ///< pbc adjustment reciprocals <xyz>
)
{
   if (narr && (ncnt>0))
   {
      for (int i=0; (i<ncnt); i++)
         narr[i].Advance(dt, cx,cy,cz, lx,ly,lz, sx,sy,sz);
   }
}

void Nodes_Advance
(
         Node_t **narr,                                 ///< array of node pointers (can be sparse)
   const int      ncnt,                                 ///< number of nodes
   const real8    dt,                                   ///< dt = next time step
   const real8    cx, const real8 cy, const real8 cz,   ///< simulation center <xyz>
   const real8    lx, const real8 ly, const real8 lz,   ///< simulation size   <xyz>
   const real8    sx, const real8 sy, const real8 sz    ///< pbc adjustment reciprocals <xyz>
)
{
   if (narr && (ncnt>0))
   {
      for (int i=0; (i<ncnt); i++)
         if (narr[i]) narr[i]->Advance(dt, cx,cy,cz, lx,ly,lz, sx,sy,sz);
   }
}

// Preserve()
//
// Will preserve various node data elements specified via the items argument.
//
// The items argument is a bitwise 'or' of the following constants...
//    NODE_POSITION : saves the current node position into the old position
//    NODE_CURR_VEL : saves the current velocity without wiping out the copy of the previous velocity
//    NODE_OLD_VEL  : saves the old velocity (normally done at the end of the timestep integrator
//---------------------------------------------------------------------------------------------------------

void Node_t::Preserve (const int items)
{
   if (items & NODE_POSITION) 
   {
      oldx = x;
      oldy = y;
      oldz = z;
   }

   if (items & NODE_CURR_VEL) 
   {
      currvX = vX;
      currvY = vY;
      currvZ = vZ;
   }

   if (items & NODE_OLD_VEL) 
   {
      oldvX = currvX;
      oldvY = currvY;
      oldvZ = currvZ;
   }
}

// Nodes_Preserve()
//
// Preserves an array of nodes or node pointers.
//---------------------------------------------------------------------------------------------------------

void Nodes_Preserve (Node_t  *narr , const int ncnt, const int items)
{
   if (narr && (ncnt>0))
   {
      for (int i=0; (i<ncnt); i++)
         narr[i].Preserve(items);
   }
}

void Nodes_Preserve (Node_t **narr , const int ncnt, const int items)
{
   if (narr && (ncnt>0))
   {
      for (int i=0; (i<ncnt); i++)
         if (narr[i]) narr[i]->Preserve(items);
   }
}

//---------------------------------------------------------------------------------------------------------
// MPI data packing and unpacking.
//
// While MPI allows for user-defined data types, getting them to reliably function can be painful.
// Instread of dealing with that, we simply serialize all the data elements of the node into a 
// vector of double precision values that are sent/received using MPI mechanisms for sending basic
// data types (e.g. int, float, double, etc.)
//
// Note that most of the data elements of a node are already double precision values so packing the 
// non-DP elements only represents a small fraction of the overall space required.
//---------------------------------------------------------------------------------------------------------

// MPI_Packed_Count()
//
// Returns the number of double precision values needed to serialize the node.  Note that this can 
// be variable due to the dynamic nature of the number of arms currently linked to this node.
//---------------------------------------------------------------------------------------------------------

int Node_t::MPI_Packed_Count (void) const
{
   int cnt=0;
   int n=numNbrs;

   cnt++;       // flags
        
   cnt+=3;      // x , y , z 
   cnt+=3;      // fX, fY, fZ
   cnt+=3;      // vX, vY, vZ
        
   cnt+=3;      // oldx  , oldy  , oldz  
   cnt+=3;      // oldvX , oldvY , oldvZ 
   cnt+=3;      // currvX, currvY, currvZ
 
   cnt+=2;      // myTag
   cnt++ ;      // myIndex

   cnt++;       // numNbrs
   cnt+=2*n;    // nbrTag
   cnt+=  n;    // nbrIndex

   cnt+=3*n;    // armfx, armfy, armfz
   cnt+=3*n;    // burgX, burgY, burgZ
   cnt+=3*n;    // nx   , ny   , nz   

   cnt+=3*n;    // sigbLoc
   cnt+=3*n;    // sigbRem

   cnt+=  n;    // armCoordIndex
   cnt+=  n;    // armNodeIndex

   cnt++;       // constraint

   cnt++;       // cellIdx
   cnt++;       // cell2Idx
   cnt++;       // cell2QentIdx

   cnt++;       // native
   cnt++;       // tmark

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
   cnt++;       // multiNodeLife;
#endif
   
#ifdef _ARLFEM
   cnt+=3;      // surfaceNorm[3];
#endif

   return(cnt);
}

// MPI_Pack()
//
// Will serialize the node as a vector of doubles that can be exchanged between MPI processes.
// Returns a pointer to the next available postion in the serialization buffer.
//---------------------------------------------------------------------------------------------------------

real8 *Node_t::MPI_Pack 
(
   real8 *q   ///< the destination buffer (current position)
) const
{
   if (q)
   {
      *q++ = (real8) flags;
           
      *q++ = x     ; *q++ = y     ; *q++ = z     ;
      *q++ = fX    ; *q++ = fY    ; *q++ = fZ    ;
      *q++ = vX    ; *q++ = vY    ; *q++ = vZ    ;
           
      *q++ = oldx  ; *q++ = oldy  ; *q++ = oldz  ;
      *q++ = oldvX ; *q++ = oldvY ; *q++ = oldvZ ;
      *q++ = currvX; *q++ = currvY; *q++ = currvZ;
    
      *q++ = (real8) myTag.domainID;
      *q++ = (real8) myTag.index   ;

      *q++ = (real8) myIndex       ;
   
      int n = numNbrs;

      *q++ = (real8) n;

      if (nbrTag  ) for (int i=0; (i<n); ++i) { *q++ = (real8) nbrTag[i].domainID;
                                                *q++ = (real8) nbrTag[i].index   ; }

      if (nbrIndex) for (int i=0; (i<n); ++i) { *q++ = (real8) nbrIndex[i]; }
      if (armfx   ) for (int i=0; (i<n); ++i) { *q++ =         armfx   [i]; }
      if (armfy   ) for (int i=0; (i<n); ++i) { *q++ =         armfy   [i]; }
      if (armfz   ) for (int i=0; (i<n); ++i) { *q++ =         armfz   [i]; }

      if (burgX   ) for (int i=0; (i<n); ++i) { *q++ =         burgX   [i]; }
      if (burgY   ) for (int i=0; (i<n); ++i) { *q++ =         burgY   [i]; }
      if (burgZ   ) for (int i=0; (i<n); ++i) { *q++ =         burgZ   [i]; }

      if (nx      ) for (int i=0; (i<n); ++i) { *q++ =         nx      [i]; }
      if (ny      ) for (int i=0; (i<n); ++i) { *q++ =         ny      [i]; }
      if (nz      ) for (int i=0; (i<n); ++i) { *q++ =         nz      [i]; }

      if (sigbLoc)
      for (int i=0,k=0; (i<n); ++i) { *q++ = sigbLoc[k++];    // <x>
                                      *q++ = sigbLoc[k++];    // <y>
                                      *q++ = sigbLoc[k++]; }  // <z>

      if (sigbRem)
      for (int i=0,k=0; (i<n); ++i) { *q++ = sigbRem[k++];    // <x>
                                      *q++ = sigbRem[k++];    // <y>
                                      *q++ = sigbRem[k++]; }  // <z>
   
      if (armCoordIndex)
      for (int i=0; (i<n); ++i) { *q++ = (real8) armCoordIndex[i]; }
   
      if (armNodeIndex )
      for (int i=0; (i<n); ++i) { *q++ = (real8) armNodeIndex [i]; }
   
      *q++ = (real8) constraint  ;
   
      *q++ = (real8) cellIdx     ; 
      *q++ = (real8) cell2Idx    ;
      *q++ = (real8) cell2QentIdx;
   
      *q++ = (real8) native      ;
      *q++ = (real8) tmark       ;
   
#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
      *q++ = (real8) multiNodeLife;
#endif
   
#ifdef _ARLFEM
      *q++ = surfaceNorm[0];
      *q++ = surfaceNorm[1];
      *q++ = surfaceNorm[2];
#endif
   }

   return(q);  // return next available position in the serialization buffer
}

// MPI_Unpack()
//
// Will deserialize a node from a serialization buffer.
//---------------------------------------------------------------------------------------------------------

real8 *Node_t::MPI_Unpack 
(
   real8 *p    ///< points to a source serialized node
)
{
   if (p)
   {
      flags  = (int) *p++;
           
      x      = *p++; y      = *p++; z      = *p++;
      fX     = *p++; fY     = *p++; fZ     = *p++;
      vX     = *p++; vY     = *p++; vZ     = *p++;
           
      oldx   = *p++; oldy   = *p++; oldz   = *p++;
      oldvX  = *p++; oldvY  = *p++; oldvZ  = *p++;
      currvX = *p++; currvY = *p++; currvZ = *p++;
    
      myTag.domainID = (int) *p++;
      myTag.index    = (int) *p++;

      myIndex        = (int) *p++;
   
      int n = (int) *p++;

      Allocate_Arms(n);

      if (nbrTag  ) for (int i=0; (i<n); ++i) { nbrTag[i].domainID = (int) *p++;
                                                nbrTag[i].index    = (int) *p++; }

      if (nbrIndex) for (int i=0; (i<n); ++i) { nbrIndex[i] = (int) *p++; }
      if (armfx   ) for (int i=0; (i<n); ++i) { armfx   [i] =       *p++; }
      if (armfy   ) for (int i=0; (i<n); ++i) { armfy   [i] =       *p++; }
      if (armfz   ) for (int i=0; (i<n); ++i) { armfz   [i] =       *p++; }

      if (burgX   ) for (int i=0; (i<n); ++i) { burgX   [i] =       *p++; }
      if (burgY   ) for (int i=0; (i<n); ++i) { burgY   [i] =       *p++; }
      if (burgZ   ) for (int i=0; (i<n); ++i) { burgZ   [i] =       *p++; }

      if (nx      ) for (int i=0; (i<n); ++i) { nx      [i] =       *p++; }
      if (ny      ) for (int i=0; (i<n); ++i) { ny      [i] =       *p++; }
      if (nz      ) for (int i=0; (i<n); ++i) { nz      [i] =       *p++; }

      if (sigbLoc) 
      for (int i=0,k=0; (i<n); ++i) { sigbLoc[k++] = *p++;    // <x>
                                      sigbLoc[k++] = *p++;    // <y>
                                      sigbLoc[k++] = *p++; }  // <z>

      if (sigbRem) 
      for (int i=0,k=0; (i<n); ++i) { sigbRem[k++] = *p++;    // <x>
                                      sigbRem[k++] = *p++;    // <y>
                                      sigbRem[k++] = *p++; }  // <z>
  
      if (armCoordIndex) for (int i=0; (i<n); ++i) { armCoordIndex[i] = (int) *p++; } 
      if (armNodeIndex ) for (int i=0; (i<n); ++i) { armNodeIndex [i] = (int) *p++; }
   
      constraint   = (int) *p++;
   
      cellIdx      = (int) *p++;
      cell2Idx     = (int) *p++;
      cell2QentIdx = (int) *p++;
   
      native       = (int) *p++;
      tmark        = (int) *p++;
   
#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
      multiNodeLife = (int) *p++;
#endif
   
#ifdef _ARLFEM
      surfaceNorm[0] = *p++;
      surfaceNorm[1] = *p++;
      surfaceNorm[2] = *p++;
#endif
   
   }

   return(p);  // return pointer to next node
}

//---------------------------------------------------------------------------------------------------------
// MPI node forces serialization support.
//
// Note that when node forces are exchanged, we assume no topology updates, so the pack/unpack assume
// no geometry changes need to be communicated.
//
// We assume that when forces are updated, each process will have an array of ghosts.  We therefor
// pack the tag ID to enable the comm routines to look up the tag on the update receive.
//---------------------------------------------------------------------------------------------------------

int Node_t::MPI_Pack_Forces_Count (void) const
{
   int cnt=0;
   int n=numNbrs;

   cnt+=2;      // myTag
   cnt+=3;      // fX, fY, fZ
   cnt+=3*n;    // armfx, armfy, armfz

   return(cnt);
}

real8 *Node_t::MPI_Pack_Forces (real8 *q) const
{
   if (q)
   {
      *q++ = (real8) myTag.domainID;
      *q++ = (real8) myTag.index   ;
   
      *q++ = fX;
      *q++ = fY;
      *q++ = fZ;

      int n = numNbrs;

      if (armfx) { for (int i=0; (i<n); ++i) { *q++ = armfx[i]; } }
      if (armfy) { for (int i=0; (i<n); ++i) { *q++ = armfy[i]; } }
      if (armfz) { for (int i=0; (i<n); ++i) { *q++ = armfz[i]; } }
   }

   return(q);
}

real8 *Node_t::MPI_Unpack_Forces (real8 *p)
{
   if (p)
   {
      myTag.domainID = (int) *p++;
      myTag.index    = (int) *p++;
  
      fX = *p++;
      fY = *p++;
      fZ = *p++;

      int n = numNbrs;

      if (armfx) { for (int i=0; (i<n); ++i) { armfx[i] = *p++; } }
      if (armfy) { for (int i=0; (i<n); ++i) { armfy[i] = *p++; } }
      if (armfz) { for (int i=0; (i<n); ++i) { armfz[i] = *p++; } }
   }

   return(p);
}

//---------------------------------------------------------------------------------------------------------

static char *str_find(char *p, const char *pe, const char c)
{
   if ( p && pe && (p<pe))
   {
      while( (p<pe) && (*p!=c))  p++;
   }

   return( (p<pe) ? p : 0);
}

static char *str_tok(char *p, const char *pe)
{
   if (p && pe && (p<pe))
   {
      while( (p<pe) && ( isspace(*p)) ) p++;  // skip over leading whitespace
      while( (p<pe) && (!isspace(*p)) ) p++;  // skip over current token
      while( (p<pe) && ( isspace(*p)) ) p++;  // skip over trailing whitespace
   }

   return( (p<pe) ? p : 0);
}

static char *str_getln(char *p, const char *pe)
{
   if (p && pe && (p<pe))
   {
      while( (p<pe) && (*p!='\r') && (*p!='\n') ) p++; 

      if ( (p<pe) && (*p=='\r') ) p++;
      if ( (p<pe) && (*p=='\n') ) p++;
   }

   if (p && (p<pe) && (*p=='#') ) p=str_getln(p,pe);

   return( p && (p<pe) ? p : 0 );
}

// Parse()
//
// Will initialize a node from a source text string.
//---------------------------------------------------------------------------------------------------------

char *Node_t::Parse (char *p, const char *pe)
{
   if (p && pe && (p<pe))
   {
      int narms=0;

      myTag.domainID = (int  ) strtol(p,0,10); p = str_find (p,pe,',')+1;
      myTag.index    = (int  ) strtol(p,0,10); p = str_tok  (p,pe);
      x              = (real8) strtod(p,0);    p = str_tok  (p,pe);
      y              = (real8) strtod(p,0);    p = str_tok  (p,pe);
      z              = (real8) strtod(p,0);    p = str_tok  (p,pe);
      narms          = (int  ) strtol(p,0,10); p = str_tok  (p,pe);
      constraint     = (int  ) strtol(p,0,10); p = str_getln(p,pe);

      Allocate_Arms(narms);

      if (numNbrs>0)
      {
         for (int i=0; (i<numNbrs); i++)
         {
            nbrTag[i].domainID = (int  ) strtol(p,0,10); p = str_find (p,pe,',')+1;  
            nbrTag[i].index    = (int  ) strtol(p,0,10); p = str_tok  (p,pe);
            burgX [i]          = (real8) strtod(p,0);    p = str_tok  (p,pe);
            burgY [i]          = (real8) strtod(p,0);    p = str_tok  (p,pe);
            burgZ [i]          = (real8) strtod(p,0);    p = str_tok  (p,pe);
            nx    [i]          = (real8) strtod(p,0);    p = str_tok  (p,pe);
            ny    [i]          = (real8) strtod(p,0);    p = str_tok  (p,pe);
            nz    [i]          = (real8) strtod(p,0);    p = str_getln(p,pe);
         }
      }
   }

   return(p);
}

// Load()
//
// Will initialize a node from an open file descriptor
//---------------------------------------------------------------------------------------------------------

void Node_t::Load (FILE *fd, const int fmt)
{
   if (fd)
   {
      int  n=0;
      fscanf(fd," %d,%d %lf %lf %lf %d %d\n", &myTag.domainID, &myTag.index, &x, &y, &z, &n, &constraint );

      Allocate_Arms(n);

      if (numNbrs>0)
      {
         for (int i=0; (i<numNbrs); i++)
         {
            if (fmt==0) fscanf(fd, " %d,%d %lf %lf %lf %lf %lf %lf\n",           // (single line format)
                               &nbrTag[i].domainID, &nbrTag[i].index,            //
                               burgX+i, burgY+i, burgZ+i, nx+i, ny+i, nz+i );    //

            if (fmt==1) fscanf(fd, " %d,%d %lf %lf %lf\n %lf %lf %lf\n",         // (double-line format)
                               &nbrTag[i].domainID, &nbrTag[i].index,            //
                               burgX+i, burgY+i, burgZ+i, nx+i, ny+i, nz+i );    //
         }
      }
   }
}

//---------------------------------------------------------------------------------------------------------

void Node_t::Print (FILE *fd) const
{
   if (fd)
   {
      fprintf(fd," %4d,%-5d %14.6lf %14.6lf %14.6lf %d %d\n", 
                 myTag.domainID, myTag.index, x, y, z, numNbrs, constraint );

      int n = numNbrs;

      for (int i=0; (i<n); i++)
      {
         fprintf(fd,"    %4d,%-5d %11.8lf %11.8lf %11.8lf %11.8lf %11.8lf %11.8lf\n", 
                     nbrTag[i].domainID, nbrTag[i].index, 
                     burgX[i], burgY[i], burgZ[i],
                     nx[i]   , ny[i]   , nz[i]      );
      }
   }
}

//---------------------------------------------------------------------------------------------------------

void Node_t::Print (FILE *fd, const int info, const int pos, const int vel, const int arm, const int bvec, const int glide) const
{
   if (fd)
   {
      if (info)
      {
         fprintf(fd,"  node(%d,%d) bytes=%d narms=%d, ", myTag.domainID, myTag.index, (int) Bytes(), numNbrs);
   
         for (int i=0; i<numNbrs; i++) 
         { fprintf(fd,"(%d,%d) ", nbrTag[i].domainID, nbrTag[i].index); }

         fprintf(fd,"\n");
      }
       
      if (pos) // (print the nodal position) 
      { fprintf(fd,"  node(%d,%d) pos = (%.15e %.15e %.15e)\n", myTag.domainID, myTag.index, x, y, z); }
   
      if (vel) // (print the nodal velocity and total node force) 
      {
         fprintf(fd,"  node(%d,%d) v = (%.15e %.15e %.15e)\n", myTag.domainID, myTag.index, vX, vY, vZ);
         fprintf(fd,"  node(%d,%d) f = (%.15e %.15e %.15e)\n", myTag.domainID, myTag.index, fX, fY, fZ);
      }
   
      if (arm) // (print the arm specific forces) 
      {
         for (int i=0; i<numNbrs; i++) 
         {
            fprintf(fd,"  node(%d,%d) arm[%d]-> (%d %d) f = (%.15e %.15e %.15e)\n",
               myTag.domainID, myTag.index, i,       
               nbrTag[i].domainID, nbrTag[i].index,
               armfx[i], armfy[i], armfz[i]);
         }
      }
   
      if (bvec) // (print the Burger's vector for each arm of the node) 
      {
         for (int i=0; i<numNbrs; i++) 
         {
            fprintf(fd,"  node(%d,%d) arm[%d]-> (%d %d) b = (%.15e %.15e %.15e)\n",
               myTag.domainID, myTag.index, i,       
               nbrTag[i].domainID, nbrTag[i].index,
               burgX[i], burgY[i], burgZ[i]);
         }
      }
   
      if (glide) // (print the glide plane normal for each arm of the node) 
      {
         for (int i=0; i<numNbrs; i++) 
         {
            fprintf(fd,"  node(%d,%d) arm[%d]-> (%d %d) n = (%.15e %.15e %.15e)\n",
               myTag.domainID,myTag.index, i,       
               nbrTag[i].domainID,nbrTag[i].index,
               nx[i],ny[i],nz[i]);
         }
      }
   }
}

//---------------------------------------------------------------------------
//    Function:        Node_Count
//    Description:     Returns the number of nodes in a linked list of nodes.
//
//    Arguments:
//      node           points to the root node of a linked node list
//---------------------------------------------------------------------------

int Node_Count (Node_t *node)
{
   int cnt=0;

   while (node) { node=node->next; cnt++; }

   return(cnt);
}

// Node_Count()
//
// Returns the actual number of nodes in a sparse node array.
//---------------------------------------------------------------------------

int Node_Count 
(
   Node_t **narr,  ///< array of node pointers
   int      cnt    ///< number of node pointers
)
{
   int ncnt=0;

   if (narr && (cnt>0))
      for (int i=0; (i<cnt); i++) 
         if (narr[i]) ncnt++;

   return(ncnt);
}

//---------------------------------------------------------------------------
//    Function:        Node_Array
//    Description:     Returns an array of node pointers corresponding to a 
//                     linked node list. Output can be used to sort, compare, 
//                     search a linked node list as an array.
//                     Must be free'd by caller.
//
//    Arguments:
//      node           points to the root node of a linked node list
//      cnt            count of items in the array (returned)
//---------------------------------------------------------------------------

Node_t **Node_Array (Node_t *node, int & cnt)
{
   cnt = Node_Count(node);

   Node_t **narr = (Node_t **) ( (cnt>0) ? malloc(cnt*sizeof(Node_t *)) : 0 );

   MZERO(narr,cnt*sizeof(Node_t *));

   if (narr)
   {
      for (int i=0; (i<cnt); i++, node=node->next )
         narr[i]=node; 
   }

   return(narr);
}

// Nodes_Vector()
//
// Somewhat different than above. Given an array of populated nodes, creates
// a vector of pointers to those nodes.  This is typically done to sort and
// search a node array without all the overhead of actually copying physical
// nodes.
//---------------------------------------------------------------------------

Node_t **Nodes_Vector (Node_t *nodes, int ncnt)
{
   Node_t **narr = (Node_t **) ( nodes && (ncnt>0) ? new Node_t *[ncnt] : 0 );

   if (narr)
      for (int i=0; (i<ncnt); ++i ) { narr[i]=nodes+i; }

   return(narr);
}

//---------------------------------------------------------------------------
//    Function:        Nodes_Sort
//    Description:     Given an array of node pointers, will sort the elements 
//                     of the array. Uses the standard qsort() function provided
//                     in stdlib. Returns the array of pointers sorted.
//    Arguments:
//      narr           points to an array of node pointers
//      cnt            number of nodes in the array
//      cmp            a comparison function that takes pointers to two
//                     nodes for comparison. 
//                         cmp(n0,n1) returns...
//                             -1 : node 0 less than    node 1 
//                              0 : node 0 equal to     node 1 
//                             +1 : node 0 greater than node 1 
//
//                     An example comparison function is provided below.
//---------------------------------------------------------------------------

static int ncmp (const void *p0, const void *p1)
{
   Node_t *n0 = ( p0 ? *((Node_t **) p0) : 0 );
   Node_t *n1 = ( p1 ? *((Node_t **) p1) : 0 );

   if (n0 && n1)
   {
      if (n0->myTag.domainID < n1->myTag.domainID) return(-1);   // (n0<n1) (compares domains)
      if (n0->myTag.domainID > n1->myTag.domainID) return(+1);   // (n0>n1) (compares domains)
      if (n0->myTag.index    < n1->myTag.index   ) return(-1);   // (n0<n1) (compares indices)
      if (n0->myTag.index    > n1->myTag.index   ) return(+1);   // (n0>n1) (compares indices)
   }

   // fall through (n0==n1)

   return(0);
}

Node_t **Nodes_Sort (Node_t **narr, int ncnt, int (*cmp)(const void *, const void *) )
{
   if (!cmp) cmp=ncmp;

   if (narr && (ncnt>1)) { qsort(narr, ncnt, sizeof(Node_t *), cmp); }

   return(narr);
}

Node_t *Nodes_Sort (Node_t *narr, int ncnt, int (*cmp)(const void *, const void *) )
{
   if (!cmp) cmp=ncmp;

   if (narr && (ncnt>1))
   {
      Node_t **nv  = new Node_t *[ncnt];
      Node_t  *tmp = new Node_t  [ncnt];

      if (nv && tmp)
      {
         for (int i=0; (i<ncnt); i++) { nv[i] = narr+i; }

         qsort(nv,ncnt,sizeof(Node_t *),cmp);

         for (int i=0; (i<ncnt); i++) { tmp[i] = *nv[i]; }
         for (int i=0; (i<ncnt); i++) { tmp[i].next = tmp+i+1; }
                                        tmp[ncnt-1].next = 0;
      }

      delete [] narr; narr=tmp;

      if (nv) { delete [] nv; nv=0; }
   }

   return(narr);
}

// Nodes_Find()
//
// Performs a binary search of the source vector array. Assumes the source array is sorted (by Tag)
// and fully populated.
//---------------------------------------------------------------------------------------------------------

Node_t *Nodes_Find (Node_t *narr, const int ncnt, const Tag_t & tag)
{
   if (narr && (ncnt>0))
   {
      int i0=0;
      int i1=(ncnt-1);
   
      while (i0<=i1)
      {   
         int i = (i0+i1)/2;

         Node_t *node = narr+i;
         int  cmp = ( (node->myTag < tag) ? -1 : 
                      (node->myTag > tag) ? +1 : 0 );

         if (cmp==0) { return(node); }
         if (cmp< 0) { i0=i+1; }
         else        { i1=i-1; }
      }   
   }

   return(0);
}

Node_t *Nodes_Find (Node_t **narr, const int ncnt, const Tag_t & tag)
{
   if (narr && (ncnt>0))
   {
      int i0=0;
      int i1=(ncnt-1);
   
      while (i0<=i1)
      {   
         int i = (i0+i1)/2;

         Node_t *node = narr[i];
         int  cmp = ( (node->myTag < tag) ? -1 : 
                      (node->myTag > tag) ? +1 : 0 );

         if (cmp==0) { return(node); }
         if (cmp< 0) { i0=i+1; }
         else        { i1=i-1; }
      }   
   }

   return(0);
}

// Node_Append()
// 
// Appends a single node pointer to an array of node pointers. 
// Extends the array as necessary. Returns the updated array.
// 
// Important notes...
//
// There are a couple of ways to call this routine. Ultimately,
// the expense of this routine lies in the duplicative copies
// that are made as the array is extended.  When a node is appended
// to a full (or empty array), a temporary array is allocated to
// include space for the existing entries, plus space for new entries.
// The old entries are then copied to the temporary array and the
// original array is deleted.
//
// This process would be more efficient if the user knows and allocates
// the array prior to appending individual elements. 
//
// Generally - here is how the routine should be used...
//
//    Node_t **narr=0; int cnt=0, max=0;
//
//    for (int i=0; (i<...); i++)
//       narr = Node_Append(narr,cnt,max,nodes[i]);
//
//---------------------------------------------------------------------------

Node_t *Node_Append (Node_t *narr, int & cnt, int & max, Node_t &node )
{
   if ( !narr || (cnt==max) )
   {
      max = ( (max==0) ? 2000 : 2*max );

      Node_t *tmp = new Node_t[max];

      if (tmp && narr) { for (int i=0; (i<cnt); i++) tmp[i]=narr[i]; }

      if (narr) { delete [] narr; }

      narr = tmp;
   }
   
   if (narr && (cnt<max) ) { narr[cnt++] = node; }

   return(narr);
}

Node_t **Node_Append (Node_t **narr, int & cnt, int & max, Node_t  *node)
{
   if (node)
   {
      if ( !narr || (cnt==max) )
      {
         max = ( (max==0) ? 2000 : 2*max );
   
         Node_t **tmp = (Node_t **) malloc(max*sizeof(Node_t *));
  
         MZERO(tmp,     max*sizeof(Node_t *));
         MCOPY(tmp,narr,cnt*sizeof(Node_t *));
         MFREE(narr);
   
         narr = tmp;
      }
   
      if (narr && (cnt<max) ) { narr[cnt++] = (Node_t *) node; }
   }

   return(narr);
}

// Nodes_Append()
//
// Somewhat more efficient form. Appends an array of nodes to a node array.
// Extends the array as necessary. Returns the updated array.
//---------------------------------------------------------------------------

Node_t **Nodes_Append (Node_t **narr, int & cnt, int & max, Node_t **nodes, const int ncnt)
{
   if (nodes && (ncnt>0) && (narr!=nodes))
   {
      if ( !narr || ((cnt+ncnt)>=max) )
      {
         max = ( (max==0) ? ncnt : (max+ncnt) );
   
         Node_t **tmp = (Node_t **) malloc(max*sizeof(Node_t *));
   
         MCOPY(tmp,narr,cnt*sizeof(Node_t *));
         
         if (narr) { free(narr); }
   
         narr = tmp;
      }

      // Append all the non-null nodes from the source array...

      for (int i=0; (i<ncnt); i++)
      {  if (nodes[i]) { narr[cnt++] = (Node_t *) nodes[i]; } }
   }

   return(narr);
}

// Nodes_Set_Native()
//
// Will mark all the nodes within the provided domain as native or ghost.
//---------------------------------------------------------------------------

inline int Node_Inside(const real8 *p, const real8 *dom )
{
   return ( (    (p[0]<dom[0]) || (dom[1]<p[0]) 
              || (p[1]<dom[2]) || (dom[3]<p[1]) 
              || (p[2]<dom[4]) || (dom[5]<p[2]) ) ? 0 : 1 );
}

void Nodes_Set_Native 
(
         Node_t  *nodes,                ///< array of nodes
         int      ncnt,                 ///< number of nodes
   const real8 x0, const real8 x1,      ///< domain min/max (x)
   const real8 y0, const real8 y1,      ///< domain min/max (y)
   const real8 z0, const real8 z1       ///< domain min/max (z)
)
{
   if (nodes && (ncnt>0))
   {
      real8 dom[6] = { x0,x1, y0,y1, z0,z1 };

      for (int i=0; (i<ncnt); i++) 
      {
         real8 p[3] = { nodes[i].x, 
                        nodes[i].y, 
                        nodes[i].z };

         nodes[i].native = Node_Inside(p,dom);
      }
   }
}

void Nodes_Set_Native 
(
         Node_t **nodes,                ///< array of node pointers
         int      ncnt,                 ///< number of nodes
   const real8 x0, const real8 x1,      ///< domain min/max (x)
   const real8 y0, const real8 y1,      ///< domain min/max (y)
   const real8 z0, const real8 z1       ///< domain min/max (z)
)
{
   if (nodes && (ncnt>0))
   {
      real8 dom[6] = { x0,x1, y0,y1, z0,z1 };

      for (int i=0; (i<ncnt); i++) 
      {
         if (nodes[i])
         {
            real8 p[3] = { nodes[i]->x, 
                           nodes[i]->y, 
                           nodes[i]->z };

            nodes[i]->native = Node_Inside(p,dom);
         }
      }
   }
}

// Nodes_Validate()
//
// Returns whether an array of nodes is valid/invalid based on the 
// sum of the burgers vectors.
//---------------------------------------------------------------------------

int Nodes_Validate (Node_t  *narr , int ncnt)
{
   if (narr && (ncnt>0))
   {
      for (int i=0; (i<ncnt); i++)
         if (!narr[i].Validate()) return(0);
   }

   return(1);
}

int Nodes_Validate (Node_t **narr , int ncnt)
{
   if (narr && (ncnt>0))
   {
      for (int i=0; (i<ncnt); i++)
         if (narr[i] && !narr[i]->Validate() ) return(0);
   }

   return(1);
}

// Node_Forces_Reset()
//
// Will iterate through all the nodes in the node array and reset
// all the node-specific force vectors.
//---------------------------------------------------------------------------

void Node_Forces_Reset (Node_t *nodes, int ncnt)
{
   if (nodes && (ncnt>0))
      for (int i=0; (i<ncnt); i++) { nodes[i].Forces_Reset(); }
}

void Node_Forces_Reset (Node_t **narr , int ncnt)
{
   if (narr && (ncnt>0))
      for (int i=0; (i<ncnt); i++) { if (narr[i]) narr[i]->Forces_Reset(); }
}

// Node_Forces_Sum()
//
// Will iterate through all the nodes in the node array and sum all
// all the arm forces into the aggregate node force.
//---------------------------------------------------------------------------

void Node_Forces_Sum (Node_t *nodes, int ncnt)
{
   if (nodes && (ncnt>0))
      for (int i=0; (i<ncnt); i++) { nodes[i].Add_Arm_Forces(); }
}

void Node_Forces_Sum (Node_t **narr , int ncnt)
{
   if (narr && (ncnt>0))
      for (int i=0; (i<ncnt); i++) { if (narr[i]) narr[i]->Add_Arm_Forces(); }
}

//---------------------------------------------------------------------------

static char *file_load(const char *path, size_t & bytes)
{
   char *p    =0;
         bytes=0;

   FILE *fd = ( path ? fopen(path,"r") : 0 );

   if (fd)
   {
               fseek(fd,0,SEEK_END);                // seek to the end of the file
       bytes = ftell(fd);                           // bytes = current file position (tells us how much to allocate)
               fseek(fd,0,SEEK_SET);                // return to start of file
 
       p = (char *) ( (bytes>0) ? new char[bytes] : 0 );   // allocate space to receive the contents of the file
 
       if (p && (bytes>0))                                 // if (space allocated)
          fread(p,bytes,1,fd);                             //    read the entire file into memory
 
       fclose(fd);                                        // close the file
   }

   return(p);
}

// Nodes_Load()
//
// Will load an array of nodes from a source input file.
//---------------------------------------------------------------------------

Node_t *Nodes_Load (const char *path, int & ncnt)
{
   Node_t *nodes=0;  int nmax=0; ncnt=0;

   size_t  bytes = 0;
   char   *p     = file_load(path,bytes);
   char   *pe    = ( p ? (p+bytes) : 0 );
   char   *p0    = p;

   if (p && pe && (p<pe))
   {
      while (p && (p<pe))
      {
         Node_t n;  p=n.Parse(p,pe);

         nodes = Node_Append(nodes,ncnt,nmax,n);
      }
   }

   if (p0) { delete [] p0; p0=0; }

   return(nodes);
}

// Nodes_Save()
//
// Will save an array of nodes to an output file.
//---------------------------------------------------------------------------

void Nodes_Save (const char *path, Node_t *nodes, int ncnt)
{
   FILE *fd = ( path ? fopen(path,"w") : 0 );

   if (fd && nodes && (ncnt>0))
   {
      for (int i=0; (i<ncnt); i++) { nodes[i].Print(fd); }
   }

   if (fd) fclose(fd);
}

// Nodes_Save()
//
// Will save an array of nodes to an gnuplot compatible visualization file.
//---------------------------------------------------------------------------

void Nodes_Save_GNU (const char *path, Node_t *nodes, int ncnt, const real8 lx, const real8 ly, const real8 lz)
{
   real8 lx2=(lx/2.0);
   real8 ly2=(ly/2.0);
   real8 lz2=(lz/2.0);

   char path_dat[256]; if (path) { sprintf(path_dat,"%s.dat", path); }
   char path_gnu[256]; if (path) { sprintf(path_gnu,"%s.gnu", path); }

   FILE *fd = ( path ? fopen(path_dat,"w") : 0 );

   if (fd && nodes && (ncnt>0))
   {
      Node_t *n1 = (Node_t *) nodes;
     
      for (int i=0; (i<ncnt); i++, n1++)
      { 
         Tag_t *tags = n1->nbrTag;

         for (int j=0; (j<n1->numNbrs);  j++)
         {
            Node_t *n2 = Nodes_Find(nodes,ncnt,tags[j]);

            if (n1 && n2 && (*n1 < *n2) )
            {
               real8     p1 [3] = { n1->x, n1->y, n1->z };
               real8     p2 [3] = { n2->x, n2->y, n2->z };
               real8     p12[3] = { p2[0]-p1[0], 
                                    p2[1]-p1[1], 
                                    p2[2]-p1[2] };

               p12[0] += ( (p12[0] < -lx2) ? lx : ( (p12[0] > lx2) ?  -lx : 0.0 ) );    // correct for PBC (x)
               p12[1] += ( (p12[1] < -ly2) ? ly : ( (p12[1] > ly2) ?  -ly : 0.0 ) );    // correct for PBC (y)
               p12[2] += ( (p12[2] < -lz2) ? lz : ( (p12[2] > lz2) ?  -lz : 0.0 ) );    // correct for PBC (z)
     
               p2[0] = p1[0]+p12[0]; 
               p2[1] = p1[1]+p12[1]; 
               p2[2] = p1[2]+p12[2]; 

               fprintf(fd," %2d %4d %14.6lf %14.6lf %14.6lf"  , n1->myTag.domainID, n1->myTag.index, p1[0], p1[1], p1[2] );
               fprintf(fd," %2d %4d %14.6lf %14.6lf %14.6lf\n", n2->myTag.domainID, n2->myTag.index, p2[0], p2[1], p2[2] );
            }
         }
      }
   }

   if (fd) fclose(fd);

   fd = ( path ? fopen(path_gnu,"w") : 0 );

   if (fd)
   {
      fprintf(fd,"reset               \n");
      fprintf(fd,"\n"                   );
      fprintf(fd,"set view 60,30, 1.1 \n");
      fprintf(fd,"unset key           \n");
      fprintf(fd,"set xlabel 'X'      \n");
      fprintf(fd,"set ylabel 'Y'      \n");
      fprintf(fd,"set zlabel 'Z'      \n");
      fprintf(fd,"set ticslevel 0     \n");
      fprintf(fd,"set xrange [%d:%d]  \n", (int)(-lx2), (int)(lx2) );
      fprintf(fd,"set yrange [%d:%d]  \n", (int)(-ly2), (int)(ly2) );
      fprintf(fd,"set zrange [%d:%d]  \n", (int)(-lz2), (int)(lz2) );
      fprintf(fd,"\n");
      fprintf(fd,"path = '%s'\n", path_dat);
      fprintf(fd,"\n");
      fprintf(fd,"splot path u 3:4:5:($8-$3):($9-$4):($10-$5) with vectors nohead lt 1 lw 1   lc rgb '#0000ff' , \\\n");
      fprintf(fd,"      path u 3:4:5                          with points         pt 6 ps 0.5 lc rgb '#0000ff' , \\\n");
      fprintf(fd,"      path u 8:9:10                         with points         pt 6 ps 0.5 lc rgb '#0000ff' ;\n"   );
      fprintf(fd,"\n");

      fclose(fd);
   }

}

//---------------------------------------------------------------------------------------------------------
// MPI communications support to send/receive arrays of packed nodes.
//
// Node that this current API assumes a sparse pointer array.  These routines return a densely-packed 
// serialization of the nodes in the sparse node array by ignoring NIL entries in the current node list.
//---------------------------------------------------------------------------------------------------------

// MPI_Packed_Count()
//
// Returns the overall number of values needed to pack the nodes into an MPI serialization buffer.
//---------------------------------------------------------------------------------------------------------

int MPI_Packed_Count 
(
         Node_t **narr,   ///< points to a sparse array of node pointers
   const int      ncnt    ///< number of node pointers
)
{
   int cnt=0;

   if (narr && (ncnt>0))
   {
      for (int i=0; (i<ncnt); ++i)
         if (narr[i]) cnt += narr[i]->MPI_Packed_Count();
   }

   return(cnt);
}

// MPI_Packed_Count()
//
// Returns the overall number of values needed to pack the nodes into an MPI serialization buffer.
// Note - only counts nodes in the sparse array with domains that match the provided rank.
//---------------------------------------------------------------------------------------------------------

int MPI_Packed_Count 
(
         Node_t **narr,   ///< points to a sparse array of node pointers
   const int      ncnt,   ///< number of node pointers
   const int      rank    ///< mpi rank to match
)
{
   int cnt=0;

   if (narr && (ncnt>0))
   {
      for (int i=0; (i<ncnt); ++i)
         if (narr[i] && (narr[i]->myTag.domainID==rank)) 
            cnt += narr[i]->MPI_Packed_Count();
   }

   return(cnt);
}

// MPI_Packed_Count()
//
// Returns the overall number of values needed to pack the nodes into an MPI serialization buffer.
//---------------------------------------------------------------------------------------------------------

int MPI_Packed_Count 
(
   const Node_t *nodes,   ///< points to an array of nodes
   const int     ncnt     ///< number of nodes
)
{
   int cnt=0;

   if (nodes && (ncnt>0))
   {
      for (int i=0; (i<ncnt); ++i)
         cnt += nodes[i].MPI_Packed_Count();
   }

   return(cnt);
}

// MPI_Packed_Count()
//
// Returns the overall number of values needed to pack the nodes into an MPI serialization buffer.
//---------------------------------------------------------------------------------------------------------

int MPI_Packed_Count 
(
   const Node_t  *nodes,   ///< points to an array of nodes
   const int      ncnt ,   ///< number of nodes
   const int      rank     ///< mpi rank to match
)
{
   int cnt=0;

   if (nodes && (ncnt>0))
   {
      for (int i=0; (i<ncnt); ++i)
         if (nodes[i].myTag.domainID==rank) 
            cnt += nodes[i].MPI_Packed_Count();
   }

   return(cnt);
}

// MPI_Pack()
//
// Will serialize all the nodes in a node list into a serialization buffer.
// Note that if the buffer is nil - it will be allocated here.
//---------------------------------------------------------------------------------------------------------

real8 *MPI_Pack_Nodes         ///< returns a pointer to the serialized nodes.
(
         Node_t **narr    ,   ///< points to an array of source node pointers
   const int      narr_cnt    ///< number of nodes to be serialized
)
{
   int    pcnt = MPI_Packed_Count(narr,narr_cnt);
   real8 *pbuf = ( (pcnt>0) ? new real8[pcnt] : 0 );
   
   if (narr && (narr_cnt>0) && pbuf)
   {
      real8 *q = pbuf;

      for (int i=0; (i<narr_cnt); i++)
         if (narr[i]) { q = narr[i]->MPI_Pack(q); }
   }

   return(pbuf);
}

real8 *MPI_Pack_Nodes 
(
   const Node_t *nodes,  ///< array of source nodes to pack
   const int     ncnt    ///< number of nodes
)
{
   int pcnt = MPI_Packed_Count(nodes,ncnt);

   real8 *pbuf = ( (pcnt>0) ? new real8[pcnt] : 0 );
   
   if (nodes && (ncnt>0) && pbuf)
   {
      real8 *q = pbuf;

      for (int i=0; (i<ncnt); i++)
         q = nodes[i].MPI_Pack(q);
   }

   return(pbuf);
}

// MPI_Unpack()
//
// Returns an array of nodes that were previously serialized using MPI_Pack().
//---------------------------------------------------------------------------------------------------------

void MPI_Unpack_Nodes 
(
         Node_t  *nodes,  ///< array of nodes
   const int      ncnt ,  ///< number of nodes to unpack
   const real8   *pbuf    ///< source packed nodes
)
{
   if (nodes && (ncnt>0) && pbuf)
   {
      real8 *p = (real8 *) pbuf;

      for (int i=0; (i<ncnt); i++)
         p = nodes[i].MPI_Unpack(p);
   }
}

//---------------------------------------------------------------------------------------------------------
// GPU enabled support...
//---------------------------------------------------------------------------------------------------------

#ifdef GPU_ENABLED

//---------------------------------------------------------------------------------------------------------

static          size_t  mbuf_bytes_alloc = 0;   ///< bytes currently allocated on host/device
static          size_t  mbuf_bytes       = 0;   ///< current size of host/device buffers
static unsigned char   *mbuf_cpu         = 0;   ///< points to current host node buffer
static unsigned char   *mbuf_gpu         = 0;   ///< points to current gpu  node buffer

//---------------------------------------------------------------------------------------------------------

unsigned char    *Nodes_CPU_Mbuf  (void) { return(mbuf_cpu  ); }    // returns pointer to current host   node buffer
unsigned char    *Nodes_GPU_Mbuf  (void) { return(mbuf_gpu  ); }    // returns pointer to current device node buffer
         size_t   Nodes_GPU_Bytes (void) { return(mbuf_bytes); }    // returns current size of host/device buffers

// Nodes_GPU_Bytes()
//
// Given an array of nodes, will return the size of the buffer needed to gpu-pack all the nodes into
// a destination buffer.
//---------------------------------------------------------------------------------------------------------

size_t Nodes_GPU_Bytes (Node_t *nodes, int cnt)
{
   size_t bytes=0;

   if (nodes && (cnt>0))
   {
     bytes += cnt*sizeof(size_t  );      // reserve space for node offsets
     bytes += cnt*sizeof(Node_t *);      // reserve space for node pointers

     for (int i=0; (i<cnt); i++)         // reserve space for the packed nodes
        bytes += nodes[i].GPU_Bytes();   //
   }

   return(bytes);
}

// Nodes_GPU_Allocate()
//
// Given an array of nodes, manages the host/device memory buffer allocations needed to gpu-pack all
// the nodes in the buffers.
//---------------------------------------------------------------------------------------------------------

void Nodes_GPU_Allocate (Node_t *nodes, int cnt)
{
   size_t bytes=Nodes_GPU_Bytes   (nodes,cnt);
                Nodes_GPU_Allocate(bytes);
}

void Nodes_GPU_Allocate (size_t bytes)
{
   if ( (bytes>0) && (bytes>mbuf_bytes_alloc) )
   {
      if (mbuf_cpu) { CUDART_CHECK(cudaFreeHost((void *) mbuf_cpu));  mbuf_cpu=0; }
      if (mbuf_gpu) { CUDART_CHECK(cudaFree    ((void *) mbuf_gpu));  mbuf_gpu=0; }

      mbuf_bytes = bytes;

      bytes += (bytes/4);  // add a 25% growth fudge

      CUDART_CHECK(cudaHostAlloc((void **) &mbuf_cpu, bytes, cudaHostAllocPortable));
      CUDART_CHECK(cudaMalloc   ((void **) &mbuf_gpu, bytes));

      mbuf_bytes_alloc = bytes;
   }
}

// Nodes_GPU_Free()
//
// Releases all host/device memory buffer allocations for exchanging node data with the GPU.
//---------------------------------------------------------------------------------------------------------

void Nodes_GPU_Free (void)
{
   if (mbuf_cpu) { CUDART_CHECK(cudaFreeHost((void *) mbuf_cpu)); mbuf_cpu=0; }
   if (mbuf_gpu) { CUDART_CHECK(cudaFree    ((void *) mbuf_gpu)); mbuf_gpu=0; }

   mbuf_bytes_alloc=0; 
   mbuf_bytes      =0; 
   mbuf_cpu        =0;
   mbuf_gpu        =0;
}

// Nodes_GPU_Pack()
//
// Packs an array of nodes into a host buffer.
//---------------------------------------------------------------------------------------------------------

void Nodes_GPU_Pack (Node_t *nodes , int cnt, float *msec)
{
   CPU_Timer tmr; 

   if (msec) { tmr.Start(); }

   size_t bytes = Nodes_GPU_Bytes   (nodes,cnt);
                  Nodes_GPU_Allocate(bytes);

   if (nodes && (cnt>0) && mbuf_cpu && mbuf_gpu)
   {
      MZERO(mbuf_cpu,bytes);

      unsigned char    *q    = mbuf_cpu;
      unsigned char    *q0   = mbuf_cpu;
               size_t  *poff = (size_t  *) q;  q+=cnt*sizeof(size_t  );
               Node_t **pn   = (Node_t **) q;  q+=cnt*sizeof(Node_t *);

      for (int i=0; (i<cnt); i++)
      {
         q       = nodes[i].GPU_Pack(q); 
         poff[i] =           (         (size_t)(q-q0));
         pn  [i] = (Node_t *)(mbuf_gpu+(size_t)(q-q0));
      }
   }

   if (msec) { tmr.Stop(); *msec=(float) tmr.Msecs(); }
}

// Nodes_GPU_Unpack()
//
// Unpacks the gpu data from a host buffer.
//---------------------------------------------------------------------------------------------------------

void Nodes_GPU_Unpack (Node_t *nodes, int cnt, float *msec)
{
   CPU_Timer tmr; 

   if (msec) { tmr.Start(); }

   if (nodes && (cnt>0) && mbuf_cpu)
   {
      unsigned char  *p  = mbuf_cpu;
                      p += cnt*sizeof(size_t  );
                      p += cnt*sizeof(Node_t *);

      for (int i=0; (i<cnt); i++)
         p = nodes[i].GPU_Unpack(p); 
   }

   if (msec) { tmr.Stop(); *msec=(float) tmr.Msecs(); }
}

// Nodes_GPU_Send()
//
// Sends a packed array of nodes to the GPU.
//---------------------------------------------------------------------------------------------------------

void Nodes_GPU_Send (float *msec)
{
   CPU_Timer tmr; 

   if (msec) { tmr.Start(); }

   if (mbuf_cpu && mbuf_gpu && (mbuf_bytes>0) )
   {
      CUDART_CHECK(cudaMemcpy((void *) mbuf_gpu, (void *) mbuf_cpu, mbuf_bytes, cudaMemcpyHostToDevice));
      CUDART_CHECK(cudaDeviceSynchronize());
   }

   if (msec) { tmr.Stop(); *msec=(float) tmr.Msecs(); }
}

// Nodes_GPU_Recv()
//
// Receives a packed array of nodes to the GPU.
//---------------------------------------------------------------------------------------------------------

void Nodes_GPU_Recv (float *msec)
{
   CPU_Timer tmr; 

   if (msec) { tmr.Start(); }

   if (mbuf_cpu && mbuf_gpu && (mbuf_bytes>0) )
   {
      CUDART_CHECK(cudaMemcpy((void *) mbuf_cpu, (void *) mbuf_gpu, mbuf_bytes, cudaMemcpyDeviceToHost));
      CUDART_CHECK(cudaDeviceSynchronize());
   }

   if (msec) { tmr.Stop(); *msec=(float) tmr.Msecs(); }
}

#endif // GPU_ENABLED

//---------------------------------------------------------------------------------------------------------
// Stubs for non-GPU builds...
//---------------------------------------------------------------------------------------------------------

#ifndef GPU_ENABLED

unsigned char    *Nodes_CPU_Mbuf     (void)                                  { return(0); }
unsigned char    *Nodes_GPU_Mbuf     (void)                                  { return(0); }
         size_t   Nodes_GPU_Bytes    (void)                                  { return(0); }

         size_t   Nodes_GPU_Bytes    (Node_t *nodes, int cnt)                { return(0); } 
         void     Nodes_GPU_Allocate (Node_t *nodes, int cnt)                {} 
         void     Nodes_GPU_Allocate (size_t bytes)                          {} 
         void     Nodes_GPU_Free     (void)                                  {} 
         void     Nodes_GPU_Pack     (Node_t *nodes , int cnt, float *msec)  { if (msec) *msec=0.0; } 
         void     Nodes_GPU_Unpack   (Node_t *nodes , int cnt, float *msec)  { if (msec) *msec=0.0; } 
         void     Nodes_GPU_Send     (float *msec)                           { if (msec) *msec=0.0; } 
         void     Nodes_GPU_Recv     (float *msec)                           { if (msec) *msec=0.0; } 

#endif // !GPU_ENABLED
