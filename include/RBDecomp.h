#pragma once

#ifndef _PDS_RB_DECOMP_H
#define _PDS_RB_DECOMP_H

/***************************************************************************
 *
 *      Module:         RBDecomp.h
 *      Description:    Contains prototypes for the public functions
 *                      used to create, access and manipulate a domain
 *                      decomposition based on the recursive bisectioning
 *                      algorithm.
 *
 ***************************************************************************/

/*
 *      Define the ordering for the octants in this type of
 *      domain decomposition
 */
#define LLF 0    /*  Lower Left Front  */
#define LRF 1    /*  Lower Right Front */
#define ULF 2    /*  Upper Left Front  */
#define URF 3    /*  Upper Right Front */
#define LLB 4    /*  Lower Left Back   */
#define LRB 5    /*  Lower Right Back  */
#define ULB 6    /*  Upper Left Back   */
#define URB 7    /*  Upper Right Back  */


void BroadcastRBDecomp(Home_t *home, RBDecomp_t *decomp);
void DomID2DecompID(Home_t *home, int domID, char decompID[MAX_DECOMP_LVLS]);
int  FindRBDecompCoordDomain(RBDecomp_t *decomp, real8 x, real8 y, real8 z);
void FreeRBDecomp(RBDecomp_t *decomp);
void GetAllRBDecompBounds(Home_t *home, RBDecomp_t *decomp, real8 *bounds);
void GetRBDecompCellDomainList(Home_t *home, int cellID, int *domainCnt,
         int **domainList);
void GetRBDecompLocalDomainBounds(Home_t *home, RBDecomp_t *decomp);
void AllocRBDecomp(Home_t *home, RBDecomp_t *oldDecomp, RBDecomp_t **newDecomp,
         int allocType);
real8 RBCheckLoadBalance(Home_t *home, RBDecomp_t *decomp, real8 *loadData,
                         real8 imbalanceThreshold, int currLevel,
                         int *rebalanceLevel);
void RBDecomp(Home_t *home, RBDecomp_t *decomp, RBDecomp_t *oldDecomp,
         char *domSubpartID, int currLevel, int startLevel);
void ExchangeRBDecomp(Home_t *home, RBDecomp_t *decomp);
void ReadBinRBDecompBounds(Home_t *home, void *filePtr, int numXDoms,
         int numYDoms, int numZDoms, RBDecomp_t **oldDecomp);
void ReadRBDecompBounds(Home_t *home, void **filePtr, int numXDoms,
         int numYDoms, int numZDoms, int saveDecomp, RBDecomp_t **oldDecomp);
void UniformRBDecomp(Param_t *param, RBDecomp_t *decomp, int level);
void WriteRBDecompBounds(Home_t *home, FILE *fp, RBDecomp_t *decomp,
         int level);
void XPlotRBDecomp(Home_t *home, RBDecomp_t *decomp, real8 xMin,
         real8 yMin, real8 zMin, real8 lMax, int color, real8 lineWidth);

#endif
