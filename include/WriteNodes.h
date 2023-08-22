#pragma once

#ifndef _PDS_WRITENODES_H
#define _PDS_WRITENODES_H

#include <stdio.h>

//#include "Typedefs.h"
//#include "Constants.h"
//#include "Home.h"
//#include "Node.h"
//#include "Restart.h"

extern void     WriteNode    (FILE *fd, Node_t  *node);
extern void     WriteNodes   (FILE *fd, Node_t **nodes, const int n);
extern void     WriteNodes   (const char *path, const char *mode, Node_t **nodes, const int n);
extern void     WriteNodes   (Home_t *home);

extern Node_t  *ReadNodes    (const char *path, int & cnt);
extern Node_t  *LoadNodes    (const char *path, int & cnt);

extern void     DeleteNodes  (Node_t *nodes);

#endif  // _PDS_WRITENODES_H
