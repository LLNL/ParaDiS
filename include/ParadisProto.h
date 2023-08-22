#pragma once

#ifndef _PDS_PARADIS_PROTO_H
#define _PDS_PARADIS_PROTO_H

/***************************************************************************
 *
 *      Module:       ParadisProto.h
 *      Description:  This header is mainly just a dumping ground for
 *                    the miscellaneous funtion prototypes that have
 *                    not been included elsewhere.  This helps eliminate
 *                    some of the compiler whines...
 *
 ***************************************************************************/

#include "stdio.h"
#include "Tag.h"
#include "Home.h"

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

// Define some inline operations for handling 3-component vectors
//
// VECTOR_ADD         : add  the components of the second vector to the first
// VECTOR_COPY        : copy the components of the second vector to the first
// VECTOR_ZERO        : zero out the components of the vector
// VECTOR_FROM_SCALAR : copy the contents of <v> to the scalar values <x>, <y> and <z>.
// VECTOR_TO_SCALAR   : copy the contents of the three scalar values <x>, <y> and <z> into <v>

#define VECTOR_ADD(a,b)             { (a)[0]+=(b)[0]; (a)[1]+=(b)[1]; (a)[2]+=(b)[2]; }
#define VECTOR_COPY(a,b)            { (a)[0] =(b)[0]; (a)[1] =(b)[1]; (a)[2] =(b)[2]; }
#define VECTOR_ZERO(a)              { (a)[0] =0;      (a)[1] =0;      (a)[2] =0; }
#define VECTOR_FROM_SCALAR(v,x,y,z) { (v)[0] =(x); (v)[1]=(y); (v)[2]=(z); }
#define VECTOR_TO_SCALAR(x,y,z,v)   { (x)=(v)[0]; (y)=(v)[1]; (z)=(v)[2]; }

void AddTagMapping(Home_t *home, Tag_t *oldTag, Tag_t *newTag);
void AddTagMappingSorted(Home_t *home, Tag_t *oldTag, Tag_t *newTag, int (*cmpFunc)(const void *, const void *));
void ApplyNodeConstraintsToVelocity(Home_t *home, Node_t *node);
void AssignNodeToCell(Home_t *home, Node_t *node);
void Getline(char *string, int len, FILE *fp);
void GetVelocityStatistics(Home_t *home);
void BroadcastDecomp(Home_t *home, void *decomp);
int  CalcNodeVelocities(Home_t *home, int zeroOnErr, int doAll);
void CellCharge(Home_t *home);

void DeltaPlasticStrain(Home_t *home);
void DeltaPlasticStrain_BCC(Home_t *home);
void DeltaPlasticStrain_BCC2(Home_t *home);
void DeltaPlasticStrain_FCC(Home_t *home);
void DeltaPlasticStrain_HCP(Home_t *home);
void InitDeltaPlasticStrain_HCP(Home_t *home,
        real8 nanv[HCP_NUM_GLIDE_PLANES][3][HCP_NUM_GLIDE_BURG],
        real8 tanv[HCP_NUM_GLIDE_PLANES][3][HCP_NUM_GLIDE_BURG],
        real8 bnorm[HCP_NUM_GLIDE_BURG][3]);
void DeltaPlasticStrain_rhomboVa(Home_t *home);

void DistributeTagMaps(Home_t *home);
void FindPreciseGlidePlane(Home_t *home, real8 burgVecIn[3], real8 dirIn[3], real8 glidePlane[3], int allowFuzzyPlanes);
void FixRemesh(Home_t *home);
void InitCellCenters(Home_t *home);
void FreeCellCenters(void);
void FreeCorrectionTable(void);
#ifdef CALCENERGY
void FreeEnergyCorrectionTable(void);
#endif
void FreeInitArrays(Home_t *home, InData_t *inData);
void FreeInNodeArray(InData_t *inData, int numNodes);
void FreeRijm(void);
void FreeRijmPBC(void);
void GenerateOutput(Home_t *home, int stage);
void GetDensityDelta(Home_t *home);
int  GetElasticConstants(Param_t *param);

void GetMinDist      (real8 p1x, real8 p1y, real8 p1z, real8 v1x, real8 v1y, real8 v1z,
                      real8 p2x, real8 p2y, real8 p2z, real8 v2x, real8 v2y, real8 v2z,
                      real8 p3x, real8 p3y, real8 p3z, real8 v3x, real8 v3y, real8 v3z,
                      real8 p4x, real8 p4y, real8 p4z, real8 v4x, real8 v4y, real8 v4z,
                      real8 *dist2, real8 *ddist2dt, real8 *L1, real8 *L2);

void MinSegSegDist   (real8  *MinDist2, real8 *L1, real8 *L2, const real8 *x0, const real8 *x1, const real8 *y0, const real8 *y1);
void MinSegSegDist   (Home_t *home, Node_t *node1, Node_t *node2, Node_t *node3, Node_t *node4, real8 *dist2);
void GetNbrCoords    (Home_t *home, Node_t *node, int arm, real8 *x, real8 *y, real8 *z);


void GetParallelIOGroup(Home_t *home);
void GaussQuadraturePointsAndWgts(int numPoints, real8 *points, real8 *weights);
void HandleCollisions(Home_t *home);
void HCP_CA_Split(Home_t *home);

int FindMergeNodeOnSegment(Home_t *home, Node_t **node3, Node_t **node4,
                           real8 x3, real8 y3, real8 z3,
                           real8 x4, real8 y4, real8 z4, int arm34,
                           real8 L1, Node_t **mergenode2);

int FindCollisionPoint(Home_t *home, Node_t *node1, Node_t *node2,
                       real8 *x, real8 *y, real8 *z);
void RetroactiveCollisions(Home_t *home);
void PredictiveCollisions(Home_t *home);
void PredictiveCollisionsLG(Home_t *home);
void ProximityCollisions(Home_t *home);

void VacancyClimb(Home_t *home);

void HeapAdd(int **heap, int *heapSize, int *heapCnt, int value);
int  HeapPeek(int *heap, int heapCnt);
int  HeapRemove(int *heap, int *heapCnt);
void InitRemoteDomains(Home_t *home);
void InputSanity(Home_t *home);
void LoadCurve(Home_t *home, real8 deltaStress[3][3]);
void Migrate(Home_t *home);
int  NodeOwnsSeg(Home_t *home, Node_t *node1, Node_t *node2);
void ParadisStep(Home_t *home);
void ParadisFinish(Home_t *home);
void PickScrewGlidePlane(Home_t *home, real8 burgVec[3],
        real8 glidePlane[3]);
void ReadNodeDataFile(Home_t *home, InData_t *inData, char *dataFile);
void ReleaseMemory(Home_t *home);
void RemapArmTag(Home_t *home, Tag_t *oldTag, Tag_t *newTag);
void ResetGlidePlanes(Home_t *home);
void ResetTaskGeometry(Home_t *home);
void SetLatestRestart(char *fileName);
void SortNodesForCollision(Home_t *home);
void SplinterSegments(Home_t *home);
int  TagMapCompareNew(const void *a, const void *b);
int  TagMapCompareOld(const void *a, const void *b);
void PackTagMap      (Home_t *home, int **buf, int *bufSize);
void UnpackTagMap    (Home_t *home, int *buf);
void Tecplot(Home_t *home, char *baseFileName, int ioGroup, int firstInGroup,
        int writePrologue, int writeEpilogue, int numSegs);
void UniformDecomp(Home_t *home, void **decomp);
void WriteVelocity(Home_t *home, char *baseFileName, int ioGroup,
        int firstInGroup, int writePrologue, int writeEpilogue);
void WriteForce(Home_t *home, char *baseFileName, int ioGroup,
        int firstInGroup, int writePrologue, int writeEpilogue);
void WriteVisit(Home_t *home, char *baseFileName, int writePrologue,
        int writeEpilogue, int *nodesWritten, int *segsWritten);
void WriteVisitMetaDataFile(Home_t *home, char *baseFileName,
        int *groupDataCounts);

#ifdef ESHELBY
/*
 *      Prototypes for some functions dealing with eshelby inclusions
 */
int IntersectionEllipseSegment(real8  C[3], real8 Radius[3], real8 Rotation[2][3],
                               real8 p1[3], real8 p2[3], real8 r, int *ninc,
                               real8 ratio[2]);
int IntersectionInformation(Home_t *home, real8 pos1[3], real8 pos2[3], EInclusion_t *inclusion,
                            real8 newPos1[3], real8 newPos2[3], int *ninc, real8 ratio[2]);

void HandleParticleCollisions(Home_t *home);

void ReadEshelbyInclusions(Home_t *home);
void SegPartListSort(Home_t *home);
void FindNodePartIntersects(Home_t *home, Node_t *node);
int  GetBurgID(Home_t *home, Node_t *node, int segID);
int  GetIndexFromInclusionID(Home_t *home, int inclusionID);
void SegPartListAppend(Home_t *home, SegPartIntersect_t *newInfo);
void SegPartListUpdate(Home_t *home, SegPartIntersect_t **entry, Tag_t *tag1, Tag_t *tag2,
                       int inclusionIndex, int inclusionID);
SegPartIntersect_t *SegPartListLookup(Home_t *home, Tag_t *tag1, Tag_t *tag2);
SegPartIntersect_t *SegPartListUnsortedLookup(Home_t *home, Tag_t *tag1, Tag_t *tag2);
void SegPartListClear(Home_t *home);
#endif

/*
 *      Prototype functions used to create the flux out files and
 *      the associated description files.
 */
void WriteDensFlux(char *fname, Home_t *home);
void WriteDensFluxDesc_RhomboVa(Home_t *home);
void WriteDensFluxDesc_HCP(Home_t *home);
void WriteDensFluxDesc_FCC(Home_t *home);
void WriteDensFluxDesc_BCC(Home_t *home);
void WriteDensFluxDesc_BCC2(Home_t *home);

#ifdef CHECK_MEM_USAGE
void _CheckMemUsage(Home_t *home, char *msg);
#define CheckMemUsage(a,b) _CheckMemUsage((a),(b))
#else
#define CheckMemUsage(a,b) {}
#endif

#ifdef _BGQ
void  GetOptimal3DTaskGeometry(Home_t *home, int optimalTaskGeometry[3]);
#endif


#ifdef CALCENERGY
real8 SelfEnergy  (int coreOnly, real8 MU, real8 NU,
                   real8 bx, real8 by, real8 bz,
                   real8 x1, real8 y1, real8 z1,
                   real8 x2, real8 y2, real8 z2,
                   real8 a,  real8 Ecore);

real8 ExtPKEnergy (real8 s[3][3], real8 E, real8 NU,
                   real8 bx, real8 by, real8 bz,
                   real8 x1, real8 y1, real8 z1,
                   real8 x2, real8 y2, real8 z2);

real8 SegSegEnergy(real8 p1x, real8 p1y, real8 p1z,
                   real8 p2x, real8 p2y, real8 p2z,
                   real8 p3x, real8 p3y, real8 p3z,
                   real8 p4x, real8 p4y, real8 p4z,
                   real8 bpx, real8 bpy, real8 bpz,
                   real8 bx , real8 by , real8 bz ,
                   real8 a  , real8 MU , real8 NU  );

real8 ComputeEnergy(Home_t *home, Node_t *node1, Node_t *node2, Node_t *node3, Node_t *node4);

real8 LineTensionEnergy(Home_t *home,
                        real8 x1, real8 y1, real8 z1,
                        real8 x2, real8 y2, real8 z2,
                        real8 bx, real8 by, real8 bz );
#endif

#ifdef FRCRIT
void Calc_FR_Critical_Stress(int argc, char *argv[], Home_t **homeptr);
#endif


#endif  // _PDS_PARADIS_PROTO_H
