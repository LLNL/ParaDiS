#pragma once

#ifndef _PDS_NODE_BLOCK_H
#define _PDS_NODE_BLOCK_H

#include "Typedefs.h"
#include "Node.h"

class NodeBlock_t
{
   public:
      Node_t      *nodes;  ///< array of nodes allocated to this block
      NodeBlock_t *next ;  ///< points to next block of nodes

   public:
      NodeBlock_t(void) { nodes=0; next=0; } 
      NodeBlock_t(const int n);

     ~NodeBlock_t();
};

#endif  // _PDS_NODE_BLOCK_H
