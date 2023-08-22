
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "Typedefs.h"
#include "Tag.h"
#include "Node.h"
#include "NodeBlock.h"

//---------------------------------------------------------------------------

NodeBlock_t::NodeBlock_t(const int n)
{
   nodes = ( (n>0) ? new Node_t[n] : 0 );
   next  = 0;

   if (nodes)
   {
      Node_t *node=nodes;

      for (int i=0; (i<(n-1)); i++, node++)
         node->next = (node+1);

      node->next=0;
   }
}

NodeBlock_t::~NodeBlock_t()
{
   if (nodes) { delete [] nodes; nodes=0; }

   nodes=0;
   next =0;
}

#if 0
//---------------------------------------------------------------------------
//    Function:        Node_Block_New
//    Description:     Creates and initializes a new block of nodes
//
//    Arguments:
//      n              number of nodes to create in the block
//---------------------------------------------------------------------------

NodeBlock_t  *Node_Block_New (const int n)
{
   // allocate space for the node block structure ...

   NodeBlock_t *nb = (NodeBlock_t *) malloc(sizeof(NodeBlock_t));

   // allocate space for a block of new nodes and initialize...
    
   nb->next  = 0;
   nb->nodes = (Node_t *)calloc(1, n*sizeof(Node_t));
    
   // cycle through all the newly allocated nodes and assign a unique index 
   // and initialize the node mutex lock (for OpenMP access control)
 
   Node_t *node = nb->nodes;
   for (int i=0; (i<n); i++, node++) 
   {
      INIT_LOCK(&node->nodeLock);
      node->next = ( (i<(n-1)) ? (node+1) : 0 );
   }

   return(nb); 
}
#endif

