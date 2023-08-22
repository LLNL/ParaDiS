#pragma once

#ifndef _PDS_RS_DECOMP_H
#define _PDS_RS_DECOMP_H

/***************************************************************************
 *
 *      Module:         RSDecomp.h
 *      Description:    Contains prototypes for the public functions
 *                      used to create, access and manipulate a domain
 *                      decomposition based on the recursive sectioning
 *                      algorithm.
 *
 ***************************************************************************/

void BroadcastRSDecomp(Home_t *home, RSDecomp_t *decomp);
void DLBalanceX(Home_t *home, real8 *loadData);
void DLBalanceY(Home_t *home, real8 *loadData);
void DLBalanceZ(Home_t *home, real8 *loadData);
int  FindRSDecompCoordDomain(Home_t *home, RSDecomp_t *decomp, real8 x,
         real8 y, real8 z);
void FreeRSDecomp(Home_t *home, RSDecomp_t *decomp);
void GetAllRSDecompBounds(Home_t *home, RSDecomp_t *decomp, real8 *bounds);
void GetRSDecompCellDomainList(Home_t *home, int cellID, int *domCount,
         int **domList);
void GetRSDecompLocalDomainBounds(Home_t *home, RSDecomp_t *decomp);
void ReadBinRSDecompBounds(void *filePtr, int numXDoms, int numYDoms,
         int numZDoms, RSDecomp_t **oldDecomp);
void ReadRSDecompBounds(void **filePtr, int numXDoms, int numYDoms,
         int numZDoms, int saveDecomp, RSDecomp_t **oldDecomp);
void UniformRSDecomp(Param_t *param, RSDecomp_t **uniDecomp);
void WriteRSDecompBounds(Home_t *home, FILE *fp, RSDecomp_t *decomp);
void XPlotRSDecomp(Home_t *home, real8 xMin, real8 yMin, real8 zMin,
         real8 lMax, int color, real8 lineWidth);

#endif
