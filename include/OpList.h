#pragma once

#ifndef _PDS_OPLIST_H
#define _PDS_OPLIST_H

/***************************************************************************
 *
 *  OpList.h Define the structures that holds all data needed
 *           to handle a cross-domain topology changes, and prototype
 *           the needed functions.
 *
 **************************************************************************/

#define REM_OP_INCR_LEN 4096

typedef enum {
        REM_OP_CHANGE_CONN,
        REM_OP_INSERT_ARM,
        REM_OP_REMOVE_NODE,
        REM_OP_CHANGE_BURG,
        REM_OP_SPLIT_NODE,
        REM_OP_RESET_COORD,
        REM_OP_RESET_SEG_FORCE1,
        REM_OP_RESET_SEG_FORCE2,
        REM_OP_MARK_FORCE_OBSOLETE,
        REM_OP_RESET_PLANE
} OpType_t;

typedef struct {
    int   opType;
    Tag_t tag1;
    Tag_t tag2;
    Tag_t tag3;
} RemOpChangeConn_t;

typedef struct {
    int   opType;
    Tag_t tag1;
    Tag_t tag2;
    real8 burg[3];
    real8 plane[3];
} RemOpInsertArm_t;

typedef struct {
    int   opType;
    Tag_t tag1;
} RemOpRemoveNode_t;

typedef struct {
    int   opType;
    Tag_t tag1;
    Tag_t tag2;
    real8 newBurg[3];
    real8 newPlane[3];
} RemOpChangeBurg_t;

typedef struct {
    int   opType;
    Tag_t tag1;
    Tag_t tag2;
    real8 newCoord[3];
    real8 newVelocity[3];
} RemOpSplitNode_t;

typedef struct {
    int   opType;
    Tag_t tag1;
    real8 newCoord[3];
} RemOpResetCoord_t;

typedef struct {
    int   opType;
    Tag_t tag1;
    Tag_t tag2;
    real8 f1[3];
} RemOpResetSegForce1_t;

typedef struct {
    int   opType;
    Tag_t tag1;
    Tag_t tag2;
    real8 f1[3];
    real8 f2[3];
} RemOpResetSegForce2_t;

typedef struct {
    int   opType;
    Tag_t tag1;
} RemOpMarkForceObsolete_t;

typedef struct {
    int   opType;
    Tag_t tag1;
    Tag_t tag2;
    real8 newPlane[3];
} RemOpResetPlane_t;


/*
 *      Prototype functions related to managing the remote operation list
 */
void AddOpChangeBurg(Home_t *home, Tag_t *tag1, Tag_t *tag2, real8 newBurg[3],
        real8 newPlane[3]);
void AddOpChangeConn(Home_t *home, Tag_t *tag1, Tag_t *tag2, Tag_t *tag3);
void AddOpInsertArm(Home_t *home, Tag_t *tag1, Tag_t *tag2, real8 newBurg[3],
        real8 newPlane[3]);
void AddOpMarkForceObsolete(Home_t *home, Tag_t *tag1);
void AddOpRemoveNode(Home_t *home, Tag_t *tag1);
void AddOpResetCoord(Home_t *home, Tag_t *tag1, real8 newCoord[3]);
void AddOpResetPlane(Home_t *home, Tag_t *tag1, Tag_t *tag2,
        real8 newplane[3]);
void AddOpResetSegForce1(Home_t *home, Tag_t *tag1, Tag_t *tag2, real8 f1[3]);
void AddOpResetSegForce2(Home_t *home, Tag_t *tag1, Tag_t *tag2,
        real8 f1[3], real8 f2[3]);
void AddOpSplitNode(Home_t *home, Tag_t *tag1, Tag_t *tag2, real8 newCoord[3],
        real8 newVelocity[3]);

void ClearOpList(Home_t *home);
void ExtendOpList(Home_t *home);
void FreeOpList(Home_t *home);
void InitOpList(Home_t *home);

#endif /* _OPLIST_H */
