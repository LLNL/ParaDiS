#pragma once

#ifndef _PDS_WRITEDOMAINS_H
#define _PDS_WRITEDOMAINS_H

#include <stdio.h>

#include "Home.h"
#include "Decomp.h"

extern void WriteDomains     (FILE *fd, const RSDecomp_t *rsd, const int nx, const int ny, const int nz);
extern void WriteDomains     (FILE *fd, const RBDecomp_t *rbd, const int nx, const int ny, const int nz);
extern void WriteDomains     (FILE *fd, Home_t *home);
extern void WriteDomains     (const char *path, Home_t *home);
extern void WriteDomains_GNU (const char *path);
extern void WriteDomains     (Home_t *home);

#endif  // _PDS_WRITEDOMAINS_H
