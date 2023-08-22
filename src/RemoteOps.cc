/*---------------------------------------------------------------------------
 *
 *  Module:       RemoteOps.cc
 *  Description:  This module contains functions used for managing
 *                the list of local operations that must be sent
 *                from the local domain to remote domains for
 *                processing.
 *
 *  Includes public functions:
 *      ClearOpList()
 *      ExtendOpList()
 *      FreeOpList()
 *      InitOpList()
 *
 *      AddOpChangeBurg()
 *      AddOpChangeConn()
 *      AddOpInsertArm()
 *      AddOpMarkForceObsolete()
 *      AddOpRemoveNode()
 *      AddOpResetCoord()
 *      AddOpResetPlane()
 *      AddOpResetSegForce1()
 *      AddOpResetSegForce2()
 *      AddOpSplitNode()
 *
 *-------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"


/*-------------------------------------------------------------------------
 *
 *      Function:      ClearOpList
 *      Description:   Zero's both the list of operations which is
 *                     communicated to remote domains during various
 *                     stages of topological changes and the count
 *                     of operations on the list .
 *
 *------------------------------------------------------------------------*/
void ClearOpList (Home_t *home)
{
        memset(home->opListBuf, 0, home->opListUsedLen);
        home->opListUsedLen = 0;

        return;
}

/*-------------------------------------------------------------------------
 *
 *  Function:      FreeOpList
 *  Description:   Deallocates the memory dedicated to the operation
 *                 list sent to remote domains after topological
 *                 changes.
 *
 *------------------------------------------------------------------------*/
void FreeOpList (Home_t *home)
{

    if (home->opListBuf != (char *)NULL) {
        free(home->opListBuf);
        home->opListBuf = (char *)NULL;
    }

    home->opListAllocLen = 0;
    home->opListUsedLen  = 0;

    return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:      InitOpList
 *      Description:   Allocates an initial block of memory for the
 *                     list of operations this domain will send to
 *                     remote domains for processing during the
 *                     various stages of topological changes.  
 *
 *                     This function should only need to be called 
 *                     one time during initialization of the application.
 *                     Adding operations to the list as well as altering
 *                     the list size will be handled dynamically as
 *                     necessary
 *
 *------------------------------------------------------------------------*/
void InitOpList (Home_t *home)
{

        home->opListBuf      = (char *)calloc(1, REM_OP_INCR_LEN);
        home->opListAllocLen = REM_OP_INCR_LEN;
        home->opListUsedLen  = 0;

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:      ExtendOpList
 *      Description:   Increases the buffer size dedicated to the
 *                     list of operations sent to remote domains
 *                     after topological changes.
 *
 *------------------------------------------------------------------------*/
void ExtendOpList (Home_t *home)
{
        home->opListAllocLen = home->opListAllocLen + REM_OP_INCR_LEN;
        home->opListBuf = (char *)realloc(home->opListBuf,
                                          home->opListAllocLen);
        return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:     AddOpChangeBurg 
 *  Description:  Add a locally performed 'CHANGE_BURG' operation to
 *                the list to be sent to remote domains.  This operation
 *                will only change the segment's burgers vector for the
 *                node associated with 'tag1'.  A second call must be
 *                made to change the segment's burgers vector for the node
 *                associated with 'tag2'
 *
 *  Parameters:
 *      IN:  tag1     Tag associated with the first endpoint of the
 *                    segment whose burgers vector is changing
 *      IN:  tag2     Tag associated with the second endpoint of the
 *                    segment whose burgers vector is changing
 *      IN:  newBurg  New burgers vector for the segment
 *      IN:  newPlane New plane for the segment
 *
 *------------------------------------------------------------------------*/
void AddOpChangeBurg(Home_t *home, Tag_t *tag1, Tag_t *tag2, real8 newBurg[3],
                     real8 newPlane[3])
{

/*
 *  For threaded simulations, we must insure only 1 thread
 *  updates the Op list at a time.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_ADD_OP)
#endif
    {
        int               opDataLen = sizeof(RemOpChangeBurg_t);
        int               neededLen = home->opListUsedLen + opDataLen;
        RemOpChangeBurg_t opData;

        if (home->opListAllocLen < neededLen) {
            ExtendOpList(home);
        }

        opData.opType = REM_OP_CHANGE_BURG;

        opData.tag1.domainID = tag1->domainID;
        opData.tag1.index    = tag1->index;

        opData.tag2.domainID = tag2->domainID;
        opData.tag2.index    = tag2->index;

        for (int i = 0; i < 3; ++i) {
            opData.newBurg[i]  = newBurg[i];
            opData.newPlane[i] = newPlane[i];
        }

        memcpy(&home->opListBuf[home->opListUsedLen], &opData, opDataLen);

        home->opListUsedLen += opDataLen;

    }  // end omp critical section

    return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:     AddOpChangeConn 
 *  Description:  Add a locally performed 'CHANGE_CONNECTION' operation to
 *                the list to be sent to remote domains.  This operation
 *                will switch the second endpoint of the <tag1>/<tag2>
 *                segment from the node <tag2> to the node <tag3>.
 *
 *  Parameters:
 *      IN:  tag1     Tag associated with the first endpoint of the
 *                    segment to be modified
 *      IN:  tag2     Tag associated with the second endpoint of the
 *                    segment to be modified
 *      IN:  tag3     Tag associated with the node which will replace
 *                    <tag2> as the endpoint of the segment.
 *
 *------------------------------------------------------------------------*/
void AddOpChangeConn(Home_t *home, Tag_t *tag1, Tag_t *tag2, Tag_t *tag3)
{

/*
 *  For threaded simulations, we must insure only 1 thread
 *  updates the Op list at a time.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_ADD_OP)
#endif
    {
        int               opDataLen = sizeof(RemOpChangeConn_t);
        int               neededLen = home->opListUsedLen + opDataLen;
        RemOpChangeConn_t opData;

        if (home->opListAllocLen < neededLen) {
            ExtendOpList(home);
        }

        opData.opType = REM_OP_CHANGE_CONN;

        opData.tag1.domainID = tag1->domainID;
        opData.tag1.index    = tag1->index;

        opData.tag2.domainID = tag2->domainID;
        opData.tag2.index    = tag2->index;

        opData.tag3.domainID = tag3->domainID;
        opData.tag3.index    = tag3->index;

        memcpy(&home->opListBuf[home->opListUsedLen], &opData, opDataLen);

        home->opListUsedLen += opDataLen;

    }  // end omp critical section

    return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:     AddOpInsertArm 
 *  Description:  Add a locally performed 'INSERT_ARM' operation to
 *                the list to be sent to remote domains.  This operation
 *                will create a new segment in the node associated with 
 *                <tag1> terminating at the node <tag2> with the specified
 *                burgers vector and plane.
 *                NOTE: This operation will only modify the node asscoiated
 *                with <tag1>, a second call must be made to insert the
 *                segment into the node associated with <tag2>. 
 *
 *  Parameters:
 *      IN:  tag1     Tag associated with the first endpoint of the
 *                    segment to be added
 *      IN:  tag2     Tag associated with the second endpoint of the
 *                    segment to be added
 *      IN:  newBurg  Burgers vector for the segment
 *      IN:  newPlane Plane for the segment
 *
 *------------------------------------------------------------------------*/
void AddOpInsertArm(Home_t *home, Tag_t *tag1, Tag_t *tag2, real8 newBurg[3],
                    real8 newPlane[3])
{

/*
 *  For threaded simulations, we must insure only 1 thread
 *  updates the Op list at a time.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_ADD_OP)
#endif
    {
        int              opDataLen = sizeof(RemOpInsertArm_t);
        int              neededLen = home->opListUsedLen + opDataLen;
        RemOpInsertArm_t opData;

        if (home->opListAllocLen < neededLen) {
            ExtendOpList(home);
        }

        opData.opType = REM_OP_INSERT_ARM;

        opData.tag1.domainID = tag1->domainID;
        opData.tag1.index    = tag1->index;

        opData.tag2.domainID = tag2->domainID;
        opData.tag2.index    = tag2->index;

        for (int i = 0; i < 3; ++i) {
            opData.burg[i]  = newBurg[i];
            opData.plane[i] = newPlane[i];
        }

        memcpy(&home->opListBuf[home->opListUsedLen], &opData, opDataLen);

        home->opListUsedLen += opDataLen;

    }  // end omp critical section

    return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:     AddOpMarkForceObsolete 
 *  Description:  Add a locally performed 'MARK_FORCE_OBSOLETE' operation to
 *                the list to be sent to remote domains. This operation
 *                will modify the flags associated with node <tag1> to
 *                indicate the forces for the node are out-of-date 
 *
 *  Parameters:
 *      IN:  tag1     Tag associated with the node to be updated.
 *
 *------------------------------------------------------------------------*/
void AddOpMarkForceObsolete(Home_t *home, Tag_t *tag1)
{
/*
 *  For threaded simulations, we must insure only 1 thread
 *  updates the Op list at a time.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_ADD_OP)
#endif
    {
        int   opDataLen = sizeof(RemOpMarkForceObsolete_t);
        int   neededLen = home->opListUsedLen + opDataLen;
        RemOpMarkForceObsolete_t opData;

        if (home->opListAllocLen < neededLen) {
            ExtendOpList(home);
        }

        opData.opType = REM_OP_MARK_FORCE_OBSOLETE;

        opData.tag1.domainID = tag1->domainID;
        opData.tag1.index    = tag1->index;

        memcpy(&home->opListBuf[home->opListUsedLen], &opData, opDataLen);

        home->opListUsedLen += opDataLen;

    }  // end omp critical section

    return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:     AddOpRemoveNode 
 *  Description:  Add a locally performed 'REMOVE_NODE' operation to
 *                the list to be sent to remote domains.  This operation
 *                will cause the ndoe associated with <tag1> to be
 *                removed from the simulation.
 *
 *  Parameters:
 *      IN:  tag1     Tag associated with the node to be removed.
 *
 *------------------------------------------------------------------------*/
void AddOpRemoveNode(Home_t *home, Tag_t *tag1)
{
/*
 *  For threaded simulations, we must insure only 1 thread
 *  updates the Op list at a time.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_ADD_OP)
#endif
    {
        int               opDataLen = sizeof(RemOpRemoveNode_t);
        int               neededLen = home->opListUsedLen + opDataLen;
        RemOpRemoveNode_t opData;

        if (home->opListAllocLen < neededLen) {
            ExtendOpList(home);
        }

        opData.opType = REM_OP_REMOVE_NODE;

        opData.tag1.domainID = tag1->domainID;
        opData.tag1.index    = tag1->index;

        memcpy(&home->opListBuf[home->opListUsedLen], &opData, opDataLen);

        home->opListUsedLen += opDataLen;

    }  // end omp critical section

    return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:     AddOpResetCoord 
 *  Description:  Add a locally performed 'RESET_COORD' operation to
 *                the list to be sent to remote domains.  This operation
 *                will change the coordiantes of the node associated with
 *                <tag1> to the specified coordinate values.
 *
 *  Parameters:
 *      IN:  tag1     Tag associated with the node to be updated.
 *      IN:  newCoord Array of X, Y and Z compnents of the new coordinates
 *                    for the node.
 *
 *------------------------------------------------------------------------*/
void AddOpResetCoord(Home_t *home, Tag_t *tag1, real8 newCoord[3])
{
/*
 *  For threaded simulations, we must insure only 1 thread
 *  updates the Op list at a time.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_ADD_OP)
#endif
    {
        int               opDataLen = sizeof(RemOpResetCoord_t);
        int               neededLen = home->opListUsedLen + opDataLen;
        RemOpResetCoord_t opData;

        if (home->opListAllocLen < neededLen) {
            ExtendOpList(home);
        }

        opData.opType = REM_OP_RESET_COORD;

        opData.tag1.domainID = tag1->domainID;
        opData.tag1.index    = tag1->index;

        for (int i = 0; i < 3; ++i) {
            opData.newCoord[i]  = newCoord[i];
        }

        memcpy(&home->opListBuf[home->opListUsedLen], &opData, opDataLen);

        home->opListUsedLen += opDataLen;

    }  // end omp critical section

    return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:     AddOpResetPlane 
 *  Description:  Add a locally performed 'RESET_PLANE' operation to
 *                the list to be sent to remote domains.  This operation
 *                will change the glide plane normal of the node
 *                associated with <tag1> to the specified plane values.
 *
 *  Parameters:
 *      IN:  tag1     Tag associated with the node to be updated.
 *      IN:  newPlane Array of X, Y and Z components of the new glide
 *                    plane normal for the node.
 *
 *------------------------------------------------------------------------*/
void AddOpResetPlane(Home_t *home, Tag_t *tag1, Tag_t *tag2, real8 newPlane[3])
{
/*
 *  For threaded simulations, we must insure only 1 thread
 *  updates the Op list at a time.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_ADD_OP)
#endif
    {
        int               opDataLen = sizeof(RemOpResetPlane_t);
        int               neededLen = home->opListUsedLen + opDataLen;
        RemOpResetPlane_t opData;

        if (home->opListAllocLen < neededLen) {
            ExtendOpList(home);
        }

        opData.opType = REM_OP_RESET_PLANE;

        opData.tag1.domainID = tag1->domainID;
        opData.tag1.index    = tag1->index;

        opData.tag2.domainID = tag2->domainID;
        opData.tag2.index    = tag2->index;

        for (int i = 0; i < 3; ++i) {
            opData.newPlane[i]  = newPlane[i];
        }

        memcpy(&home->opListBuf[home->opListUsedLen], &opData, opDataLen);

        home->opListUsedLen += opDataLen;

    }  // end omp critical section

    return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:     AddOpResetSegForce1 
 *  Description:  Add a locally performed 'RESET_SEG_FORCE1' operation to
 *                the list to be sent to remote domains.  This operation
 *                will change the force at the node associated with
 *                <tag1> for the <tag1>/<tag2> segment to the specified
 *                force values.
 *
 *  Parameters:
 *      IN:  tag1     Tag associated with the first endpoint of the segment
 *      IN:  tag2     Tag associated with the second endpoint of the segment
 *      IN:  f1       Array of X, Y and Z components of the new force
 *                    for the <tag1>/<tag2> segment at node <tag1>.
 *
 *------------------------------------------------------------------------*/
void AddOpResetSegForce1(Home_t *home, Tag_t *tag1, Tag_t *tag2, real8 f1[3])
{
/*
 *  For threaded simulations, we must insure only 1 thread
 *  updates the Op list at a time.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_ADD_OP)
#endif
    {
        int                   opDataLen = sizeof(RemOpResetSegForce1_t);
        int                   neededLen = home->opListUsedLen + opDataLen;
        RemOpResetSegForce1_t opData;

        if (home->opListAllocLen < neededLen) {
            ExtendOpList(home);
        }

        opData.opType = REM_OP_RESET_SEG_FORCE1;

        opData.tag1.domainID = tag1->domainID;
        opData.tag1.index    = tag1->index;

        opData.tag2.domainID = tag2->domainID;
        opData.tag2.index    = tag2->index;

        for (int i = 0; i < 3; ++i) {
            opData.f1[i]  = f1[i];
        }

        memcpy(&home->opListBuf[home->opListUsedLen], &opData, opDataLen);

        home->opListUsedLen += opDataLen;

    }  // end omp critical section

    return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:     AddOpResetSegForce2 
 *  Description:  Add a locally performed 'RESET_SEG_FORCE2' operation to
 *                the list to be sent to remote domains.  This operation
 *                will change the force at the nodes associated with
 *                <tag1> and <tag2> for the corresponding segment
 *
 *  Parameters:
 *      IN:  tag1     Tag associated with the first endpoint of the segment
 *      IN:  tag2     Tag associated with the second endpoint of the segment
 *      IN:  f1       Array of X, Y and Z components of the new force
 *                    for the <tag1>/<tag2> segment at node <tag1>.
 *      IN:  f2       Array of X, Y and Z components of the new force
 *                    for the <tag1>/<tag2> segment at node <tag2>.
 *
 *------------------------------------------------------------------------*/
void AddOpResetSegForce2(Home_t *home, Tag_t *tag1, Tag_t *tag2,
                         real8 f1[3], real8 f2[3])
{
/*
 *  For threaded simulations, we must insure only 1 thread
 *  updates the Op list at a time.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_ADD_OP)
#endif
    {
        int                   opDataLen = sizeof(RemOpResetSegForce2_t);
        int                   neededLen = home->opListUsedLen + opDataLen;
        RemOpResetSegForce2_t opData;

        if (home->opListAllocLen < neededLen) {
            ExtendOpList(home);
        }

        opData.opType = REM_OP_RESET_SEG_FORCE2;

        opData.tag1.domainID = tag1->domainID;
        opData.tag1.index    = tag1->index;

        opData.tag2.domainID = tag2->domainID;
        opData.tag2.index    = tag2->index;

        for (int i = 0; i < 3; ++i) {
            opData.f1[i]  = f1[i];
            opData.f2[i]  = f2[i];
        }

        memcpy(&home->opListBuf[home->opListUsedLen], &opData, opDataLen);

        home->opListUsedLen += opDataLen;

    }  // end omp critical section

    return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:     AddOpSplitNode
 *  Description:  Add a locally performed 'SPLIT_NODE' operation to
 *                the list to be sent to remote domains.  This operation
 *                will bisect the <tag1>/<tag2> segment creating a new
 *                node at the speficied position with the provided
 *                velocity. 
 *
 *  Parameters:
 *      IN:  tag1     Tag associated with the first endpoint of the segment
 *      IN:  tag2     Tag associated with the second endpoint of the segment
 *      IN:  newCoord Array of X, Y and Z components of the new coordinates
 *                    at which the new node should be created.
 *      IN:  NewVelocity Array of X, Y and Z components of the velocity
 *                    to be assigned to the newly created node.
 *
 *------------------------------------------------------------------------*/
void AddOpSplitNode(Home_t *home, Tag_t *tag1, Tag_t *tag2, real8 newCoord[3],
                    real8 newVelocity[3])
{
/*
 *  For threaded simulations, we must insure only 1 thread
 *  updates the Op list at a time.
 */
#ifdef _OPENMP
#pragma omp critical (CRIT_ADD_OP)
#endif
    {
        int              opDataLen = sizeof(RemOpSplitNode_t);
        int              neededLen = home->opListUsedLen + opDataLen;
        RemOpSplitNode_t opData;

        if (home->opListAllocLen < neededLen) {
            ExtendOpList(home);
        }

        opData.opType = REM_OP_SPLIT_NODE;

        opData.tag1.domainID = tag1->domainID;
        opData.tag1.index    = tag1->index;

        opData.tag2.domainID = tag2->domainID;
        opData.tag2.index    = tag2->index;

        for (int i = 0; i < 3; ++i) {
            opData.newCoord[i]     = newCoord[i];
            opData.newVelocity[i]  = newVelocity[i];
        }

        memcpy(&home->opListBuf[home->opListUsedLen], &opData, opDataLen);

        home->opListUsedLen += opDataLen;

    }  // end omp critical section

    return;
}
