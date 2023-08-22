#pragma once

#ifndef _PDS_QUEUE_OPS_H
#define _PDS_QUEUE_OPS_H

/*****************************************************************************
 *
 *  QueueOps.h   Define the Queue operation prototypes
 *
 ****************************************************************************/

Node_t *PopFreeNodeQ(Home_t *home);
void PushFreeNodeQ(Home_t *home, Node_t *node);
void PushNativeNodeQ(Home_t *home, Node_t *node);
void PushGhostNodeQ(Home_t *home, Node_t *node);
void RecycleGhostNodes(Home_t *home);
void RemoveNodeFromCellQ(Home_t *home, Node_t *node);
void RemoveNodeFromCell2Q(Home_t *home, Node_t *node);

#endif
