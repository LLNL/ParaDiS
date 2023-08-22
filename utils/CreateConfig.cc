/*****************************************************************************
 *
 *      Module:       CreateConfig.c
 *      Description:  Contains the majority of the functions to generate
 *                    nodal data for various types of initial dislocation
 *                    structures.
 *
 *      Includes functions:
 *
 *              CreateEdges()
 *              CreateHCPLoop()
 *              CreatePrismaticLoop()
 *              CreateFCCPrismaticLoop()
 *              CreateScrewConfig()
 *              CreateFiniteMixedConfig()
 *              CreateFCCConfig()
 *              CreateFCCIrradConfig()
 *              CreateFRSource()
 *              RhombohedralVaConfig()
 *
 *****************************************************************************/
#include "Home.h"
#include "InData.h"
#include "Tag.h"
#include "Util.h"
#include "ParadisGen.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>


static void IncDislocationDensity(InData_t *inData, real8 *totDislocLen)
{
        int      i, armID, nbrIndex;
        real8    xLen, yLen, zLen;
        Param_t  *param;
        Node_t   *node, *nbrNode;

        param = inData->param;

        for (i = 0; i < inData->nodeCount; i++) {

            node = &inData->node[i];

            for (armID = 0; armID < node->numNbrs; armID++) {

                nbrNode = LookupInNode(inData, &node->nbrTag[armID],
                                       &nbrIndex);

                if (nbrIndex < i) continue;

                xLen = nbrNode->x - node->x;   
                yLen = nbrNode->y - node->y;   
                zLen = nbrNode->z - node->z;   

                ZImage(param, &xLen, &yLen, &zLen);

                *totDislocLen += sqrt(xLen*xLen + yLen*yLen + zLen*zLen);
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreatePrismaticLoop
 *      Description:  Create an initial configuration consisting of one
 *                    or more BCC prismatic loops.  Loops will be either
 *                    interstitial or vacancy loops; default behavior is
 *                    to generate interstitials.
 *
 *                    Note: [1 1 1] type loops will be hexagonal loops
 *
 *      Arguments:
 *          IN:  inArgs       Structure containing values associated with the
 *                            command line arguments.
 *          OUT: totDislocLen Total length of all dislocations inserted into
 *                            the configuration.
 *
 *      Values used from the <inArgs> structure:
 *
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          loopType    Type of prismatic loops to create (i.e. 111, 100)
 *                        0 == mixture (default)
 *                        1 == [1 1 1] all types
 *                        2 == [1 0 0] [0 1 0] [0 0 1] types
 *          interstitialLoops If set 1, all loops will be interstitial
 *                      loops.  Otherwise loops will be vacancy loops.
 *          numLoops    Number of prismatic loops to create
 *          radius      radius (in units of b) of the loop?
 *          seed        Seed value for random number generator
 *          type        Type of dislocation configuration being generated.
 *
 *-------------------------------------------------------------------------*/
void CreatePrismaticLoop(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                         real8 *totDislocLen)
{
        int     id, loopIndex, burgIndex, dIndex, newNodeIndex;
        int     minBurgIndex, maxBurgIndex, lastBlock, numSegs;
        int     startRemeshIndex = 0;
        int     nbr1Index;
        int     numLoops, seed, useInterstitial, dislocType;
        int     offSet;
        real8   cubeSize, radius;
        real8   ux, uy, uz;
        real8   loopCtrX,loopCtrY,loopCtrZ;
        real8   invSqrt3, invSqrt6, twoInvSqrt6;
        real8   burg[7][3];
        real8   vec1[3], vec2[3], tmpVec[3];
        real8   tr[4][6][3];
        Param_t *param;
        Node_t  *node, *nbr1Node;


        param = inData->param;

/*
 *      Grap the command-line argument values needed for this function
 */
        cubeSize        = (real8)inArgs->cubeLength;
        radius          = inArgs->radius;
        numLoops        = inArgs->numLoops;
        useInterstitial = inArgs->interstitialLoops;
        seed            = inArgs->seed;
        dislocType      = inArgs->type;

/*
 *      All loops are now hexagonal and of [1 1 1] type burgers vectors for now
 */
        numSegs = 6;
        minBurgIndex = 0;
        maxBurgIndex = 3;
      
/*
 *      Define the sets of burger's vectors that may be used.  For each
 *      burgers vector we'll choose a glide plane normal vector which is
 *      perpendicular to the loop and along the burgers vector direction.
 */
        invSqrt3 = 1.0 / sqrt(3.0);
        invSqrt6 = 1.0 / sqrt(6.0);
        twoInvSqrt6 = 2.0 * invSqrt6;

        /*  [1 1 1]  */
        burg[0][0] =  invSqrt3;
        burg[0][1] =  invSqrt3;
        burg[0][2] =  invSqrt3;

        /*  [-1 1 1]  */
        burg[1][0] = -invSqrt3;
        burg[1][1] =  invSqrt3;
        burg[1][2] =  invSqrt3;

        /*  [1 -1 1]  */
        burg[2][0] =  invSqrt3;
        burg[2][1] = -invSqrt3;
        burg[2][2] =  invSqrt3;

        /*  [1 1 -1]  */
        burg[3][0] =  invSqrt3;
        burg[3][1] =  invSqrt3;
        burg[3][2] = -invSqrt3;

/*
 *      Define the translation vectors as 4 X 6 X 3 array 
 *
 *          tr[ ][ ][ ] - translation vector
 *             |  |  |
 *             |  |  x,y,z
 *             |  |
 *             |  6 directions
 *             |
 *             4 [111] type burgers vectors
 */
        if (useInterstitial == 1) {
/*
 *          Interstitial loops for [1 1 1] burgers vector
 */
            /* [ 1  1 -2] */
            tr[0][0][0]=  invSqrt6;
            tr[0][0][1]=  invSqrt6;
            tr[0][0][2]= -twoInvSqrt6;

            /* [-1  2 -1] */
            tr[0][1][0]= -invSqrt6;
            tr[0][1][1]=  twoInvSqrt6;
            tr[0][1][2]= -invSqrt6;

            /* [-2  1  1] */
            tr[0][2][0]= -twoInvSqrt6;
            tr[0][2][1]=  invSqrt6;
            tr[0][2][2]=  invSqrt6;

            /* [-1 -1  2] */
            tr[0][3][0]= -invSqrt6;
            tr[0][3][1]= -invSqrt6;
            tr[0][3][2]=  twoInvSqrt6;

            /* [ 1 -2  1] */
            tr[0][4][0]=  invSqrt6;
            tr[0][4][1]= -twoInvSqrt6;
            tr[0][4][2]=  invSqrt6;

            /* [ 2 -1 -1] */
            tr[0][5][0]=  twoInvSqrt6;
            tr[0][5][1]= -invSqrt6;
            tr[0][5][2]= -invSqrt6;

/*
 *          Interstitial loops for [-1 1 1] burgers vector
 */
            /* [-1  1 -2] */
            tr[1][0][0]= -invSqrt6;
            tr[1][0][1]=  invSqrt6;
            tr[1][0][2]= -twoInvSqrt6;

            /* [-2 -1 -1] */
            tr[1][1][0]= -twoInvSqrt6;
            tr[1][1][1]= -invSqrt6;
            tr[1][1][2]= -invSqrt6;

            /* [-1 -2  1] */
            tr[1][2][0]= -invSqrt6;
            tr[1][2][1]= -twoInvSqrt6;
            tr[1][2][2]=  invSqrt6;

            /* [ 1 -1  2] */
            tr[1][3][0]=  invSqrt6;
            tr[1][3][1]= -invSqrt6;
            tr[1][3][2]=  twoInvSqrt6;

            /* [ 2  1  1] */
            tr[1][4][0]=  twoInvSqrt6;
            tr[1][4][1]=  invSqrt6;
            tr[1][4][2]=  invSqrt6;

            /* [ 1  2 -1] */
            tr[1][5][0]=  invSqrt6;
            tr[1][5][1]=  twoInvSqrt6;
            tr[1][5][2]= -invSqrt6;

/*
 *          Interstitial loops for [1 -1 1] burgers vector
 */
            /* [ 1 -1 -2] */
            tr[2][0][0]=  invSqrt6;
            tr[2][0][1]= -invSqrt6;
            tr[2][0][2]= -twoInvSqrt6;

            /* [ 2  1 -1] */
            tr[2][1][0]=  twoInvSqrt6;
            tr[2][1][1]=  invSqrt6;
            tr[2][1][2]= -invSqrt6;

            /* [ 1  2  1] */
            tr[2][2][0]=  invSqrt6;
            tr[2][2][1]=  twoInvSqrt6;
            tr[2][2][2]=  invSqrt6;

            /* [-1  1  2] */
            tr[2][3][0]= -invSqrt6;
            tr[2][3][1]=  invSqrt6;
            tr[2][3][2]=  twoInvSqrt6;

            /* [-2 -1  1] */
            tr[2][4][0]= -twoInvSqrt6;
            tr[2][4][1]= -invSqrt6;
            tr[2][4][2]=  invSqrt6;

            /* [-1 -2 -1] */
            tr[2][5][0]= -invSqrt6;
            tr[2][5][1]= -twoInvSqrt6;
            tr[2][5][2]= -invSqrt6;

/*
 *          Interstitial loops for [1 1 -1] burgers vector
 */
            /* [ 1  1  2] */
            tr[3][0][0]=  invSqrt6;
            tr[3][0][1]=  invSqrt6;
            tr[3][0][2]=  twoInvSqrt6;

            /* [ 2 -1  1] */
            tr[3][1][0]=  twoInvSqrt6;
            tr[3][1][1]= -invSqrt6;
            tr[3][1][2]=  invSqrt6;

            /* [ 1 -2 -1] */
            tr[3][2][0]=  invSqrt6;
            tr[3][2][1]= -twoInvSqrt6;
            tr[3][2][2]= -invSqrt6;

            /* [-1 -1 -2] */
            tr[3][3][0]= -invSqrt6;
            tr[3][3][1]= -invSqrt6;
            tr[3][3][2]= -twoInvSqrt6;

            /* [-2  1 -1] */
            tr[3][4][0]= -twoInvSqrt6;
            tr[3][4][1]=  invSqrt6;
            tr[3][4][2]= -invSqrt6;

            /* [-1  2  1] */
            tr[3][5][0]= -invSqrt6;
            tr[3][5][1]=  twoInvSqrt6;
            tr[3][5][2]=  invSqrt6;
        } else{
/*
 *          Vacancy [111] loops
 */
            /* [ 1  1 -2] */ tr[0][0][0]=  invSqrt6;    tr[0][0][1]=  invSqrt6;    tr[0][0][2]= -twoInvSqrt6;
            /* [-1  2 -1] */ tr[0][5][0]= -invSqrt6;    tr[0][5][1]=  twoInvSqrt6; tr[0][5][2]= -invSqrt6;
            /* [-2  1  1] */ tr[0][4][0]= -twoInvSqrt6; tr[0][4][1]=  invSqrt6;    tr[0][4][2]=  invSqrt6; 
            /* [-1 -1  2] */ tr[0][3][0]= -invSqrt6;    tr[0][3][1]= -invSqrt6;    tr[0][3][2]=  twoInvSqrt6;
            /* [ 1 -2  1] */ tr[0][2][0]=  invSqrt6;    tr[0][2][1]= -twoInvSqrt6; tr[0][2][2]=  invSqrt6;
            /* [ 2 -1 -1] */ tr[0][1][0]=  twoInvSqrt6; tr[0][1][1]= -invSqrt6;    tr[0][1][2]= -invSqrt6;

/*
 *          Vacancy [-1 1 1] loops
 */
            /* [-1  1 -2] */
            tr[1][0][0]= -invSqrt6;
            tr[1][0][1]=  invSqrt6;
            tr[1][0][2]= -twoInvSqrt6;

            /* [ 1  2 -1] */
            tr[1][1][0]=  invSqrt6;
            tr[1][1][1]=  twoInvSqrt6;
            tr[1][1][2]= -invSqrt6;

            /* [ 2  1  1] */
            tr[1][2][0]=  twoInvSqrt6;
            tr[1][2][1]=  invSqrt6;
            tr[1][2][2]=  invSqrt6;

            /* [ 1 -1  2] */
            tr[1][3][0]=  invSqrt6;
            tr[1][3][1]= -invSqrt6;
            tr[1][3][2]=  twoInvSqrt6;

            /* [-1 -2  1] */
            tr[1][4][0]= -invSqrt6;
            tr[1][4][1]= -twoInvSqrt6;
            tr[1][4][2]=  invSqrt6;

            /* [-2 -1 -1] */
            tr[1][5][0]= -twoInvSqrt6;
            tr[1][5][1]= -invSqrt6;
            tr[1][5][2]= -invSqrt6;

/*
 *          Vacancy [1 -1 1] loops
 */
            /* [ 1 -1 -2] */
            tr[2][0][0]=  invSqrt6;
            tr[2][0][1]= -invSqrt6;
            tr[2][0][2]= -twoInvSqrt6;

            /* [-1 -2 -1] */
            tr[2][1][0]= -invSqrt6;
            tr[2][1][1]= -twoInvSqrt6;
            tr[2][1][2]= -invSqrt6;

            /* [-2 -1 1] */
            tr[2][2][0]= -twoInvSqrt6;
            tr[2][2][1]= -invSqrt6;
            tr[2][2][2]=  invSqrt6;

            /* [-1  1  2] */
            tr[2][3][0]= -invSqrt6;
            tr[2][3][1]=  invSqrt6;
            tr[2][3][2]=  twoInvSqrt6;

            /* [ 1  2  1] */
            tr[2][4][0]=  invSqrt6;
            tr[2][4][1]=  twoInvSqrt6;
            tr[2][4][2]=  invSqrt6;

            /* [ 2  1 -1] */
            tr[2][5][0]=  twoInvSqrt6;
            tr[2][5][1]=  invSqrt6;
            tr[2][5][2]= -invSqrt6;

/*
 *          Vacancy [1 1 -1] loops
 */
            /* [ 1  1  2] */
            tr[3][0][0]=  invSqrt6;
            tr[3][0][1]=  invSqrt6;
            tr[3][0][2]=  twoInvSqrt6;

            /* [-1  2  1] */
            tr[3][1][0]= -invSqrt6;
            tr[3][1][1]=  twoInvSqrt6;
            tr[3][1][2]=  invSqrt6;

            /* [-2  1 -1] */
            tr[3][2][0]= -twoInvSqrt6;
            tr[3][2][1]=  invSqrt6;
            tr[3][2][2]= -invSqrt6;

            /* [-1 -1 -2] */
            tr[3][3][0]= -invSqrt6;
            tr[3][3][1]= -invSqrt6;
            tr[3][3][2]= -twoInvSqrt6;

            /* [ 1 -2 -1] */
            tr[3][4][0]=  invSqrt6;
            tr[3][4][1]= -twoInvSqrt6;
            tr[3][4][2]= -invSqrt6;

            /* [ 2 -1  1] */
            tr[3][5][0]=  twoInvSqrt6;
            tr[3][5][1]= -invSqrt6;
            tr[3][5][2]=  invSqrt6;
        }

/*
 *      FIX ME!  Need to modify the code to deal with <100> type burgers
 *               vector later
 */
        offSet = 0;
        inData->nodeCount = 0;
        burgIndex = maxBurgIndex;

/*
 *      The nodal data is generated in 1 or more blocks.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes plus an offset to account for any previous blocks.
 *      The InitRemesh() calls also accept the dislocation type to set the
 *      tag.domainID value for any new nodes it adds to the node array.
 */

/*
 *      Create one loop at a time, cycling through burgers vectors as we go.
 */
        for (loopIndex = 0; loopIndex < numLoops; loopIndex++) {

            if (++burgIndex > maxBurgIndex) {
                burgIndex = minBurgIndex;
            }

/*
 *          Increase the size of the node array enough to hold
 *          all the new nodes created for this loop.
 */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += numSegs;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numSegs);

            loopCtrX = (randm(&seed)-0.5) * cubeSize;
            loopCtrY = (randm(&seed)-0.5) * cubeSize;
            loopCtrZ = (randm(&seed)-0.5) * cubeSize;


            for (id = offSet; id < (offSet + numSegs); id++) {

                node =  &inData->node[id];

/*
 *              Pick the index for the direction component of the tr array
 */
                dIndex = (id - offSet) % numSegs;

                node->x = loopCtrX + radius * tr[burgIndex][dIndex][0];
                node->y = loopCtrY + radius * tr[burgIndex][dIndex][1];
                node->z = loopCtrZ + radius * tr[burgIndex][dIndex][2];

/*
 *              Set up the node
 */
                SET_CONSTRAINTS(node->constraint, UNCONSTRAINED);
                node->myTag.domainID = dislocType;
                node->myTag.index = id + param->nodeCount;

                AllocNodeArms(node, 2);

                node->burgX[0] = burg[burgIndex][X];
                node->burgY[0] = burg[burgIndex][Y];
                node->burgZ[0] = burg[burgIndex][Z];

                node->burgX[1] = -burg[burgIndex][X];
                node->burgY[1] = -burg[burgIndex][Y];
                node->burgZ[1] = -burg[burgIndex][Z];

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = param->nodeCount +
                                        offSet + ((id+1) % numSegs);

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = param->nodeCount +
                                        offSet + ((id+numSegs-1)%numSegs);
            }

/*
 *          Need to set the glide plane normal for each segment.
 *          Couldn't do it above because we didn't have all the
 *          neighbor node positions until all the nodes in the
 *          loop were created.
 */

            for (id = offSet; id < (offSet + numSegs); id++) {

                node = &inData->node[id];

                nbr1Index = node->nbrTag[1].index - param->nodeCount;
                nbr1Node = &inData->node[nbr1Index];

                ux = nbr1Node->x - node->x;
                uy = nbr1Node->y - node->y;
                uz = nbr1Node->z - node->z;

                Normalize(&ux,&uy,&uz);

/* 
 *              l cross b gives normal vector for each segment, they may
 *              not be (110)
 */
                vec1[0] = ux;
                vec1[1] = uy;
                vec1[2] = uz;

                vec2[0] =  node->burgX[1];
                vec2[1] =  node->burgY[1];
                vec2[2] =  node->burgZ[1];

                NormalizedCrossVector(vec1, vec2, tmpVec);

                node->nx[1] = tmpVec[0];
                node->ny[1] = tmpVec[1];
                node->nz[1] = tmpVec[2];

                nbr1Node->nx[0] = tmpVec[0];
                nbr1Node->ny[0] = tmpVec[1];
                nbr1Node->nz[0] = tmpVec[2];

            }  /*  end for (id = offSet; ...)  */

            offSet += numSegs;

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths before
 *          we write the nodal data to the file
 */
            lastBlock = (loopIndex == (numLoops-1));

            if (lastBlock || inArgs->useSegmentedDataFile) {

                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);

                startRemeshIndex = 0;
                offSet = 0;
            }

        }  /* for (loopIndex = 0; ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateScrewConfig
 *      Description:  Create an initial configuration consisting of a
 *                    series of BCC screw dislocations in a periodic
 *                    simulation space.
 *
 *      Arguments:
 *          IN:  inArgs       Structure containing values associated with the
 *                            command line arguments.
 *          OUT: totDislocLen Total length of all dislocations inserted into
 *                            the configuration.
 *
 *      Values used from the <inArgs> structure:
 *
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of dislocation chains to create
 *          seed        Seed value for random number generator
 *          type        Type of dislocation configuration being generated.
 *
 *-------------------------------------------------------------------------*/
void CreateScrewConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                       real8 *totDislocLen)
{
        int      ic, np, ip, id0;
        int      burgIndex, gpIndex, newNodeIndex, lastBlock;
        int      nbr1Index, nbr2Index;
        int      numChains, seed, dislocType;
        int      startRemeshIndex = 0;
        real8    xp[3], yp[3], zp[3], cubeSize;
        real8    burg[4][3], glidePlane[4][3][3];
        Param_t  *param;
        Node_t   *node;

/*
 *      Grap the command-line argument values needed for this function
 */
        cubeSize   = (real8)inArgs->cubeLength;
        numChains  = inArgs->numChains;
        seed       = inArgs->seed;
        dislocType = inArgs->type;

        param = inData->param;

        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be > 0.\n",
                  "CreateScrewConfig", numChains);
        }

/*
 *      The nodal data is generated in 1 or more blocks.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes plus an offset to account for any previous blocks.
 *      The InitRemesh() calls also accept the dislocation type to set the
 *      tag.domainID value for any new nodes it adds to the node array.
 */

/*
 *      Set up an array of valid burgers vectors
 */
        inData->nburg = 4;

        static real8 burgRef[4][3] = {
           { 5.773503e-01,  5.773503e-01,  5.773503e-01},
           { 5.773503e-01,  5.773503e-01, -5.773503e-01},
           { 5.773503e-01, -5.773503e-01,  5.773503e-01},
           {-5.773503e-01,  5.773503e-01,  5.773503e-01}};

        static real8 glideRef[24][3] = {             
           /* for burgRef[0][*] */
           {                   0,   -0.707106781186548,    0.707106781186548},
           {   0.707106781186548,                    0,   -0.707106781186548},
           {  -0.707106781186548,    0.707106781186548,                    0},
           {   0.816496580927726,   -0.408248290463863,   -0.408248290463863},
           {  -0.408248290463863,    0.816496580927726,   -0.408248290463863},
           {  -0.408248290463863,   -0.408248290463863,    0.816496580927726},
           /* for burgRef[1][*] */
           {                   0,    0.707106781186548,    0.707106781186548},
           {  -0.707106781186548,                    0,   -0.707106781186548},
           {   0.707106781186548,   -0.707106781186548,                    0},
           {   0.816496580927726,   -0.408248290463863,    0.408248290463863},
           {  -0.408248290463863,    0.816496580927726,    0.408248290463863},
           {  -0.408248290463863,   -0.408248290463863,   -0.816496580927726},
           /* for burgRef[2][*] */
           {                   0,   -0.707106781186548,   -0.707106781186548},
           {  -0.707106781186548,                    0,    0.707106781186548},
           {   0.707106781186548,    0.707106781186548,                    0},
           {   0.816496580927726,    0.408248290463863,   -0.408248290463863},
           {  -0.408248290463863,   -0.816496580927726,   -0.408248290463863},
           {  -0.408248290463863,    0.408248290463863,    0.816496580927726},
           /* for burgRef[3][*] */
           {                   0,    0.707106781186548,   -0.707106781186548},
           {   0.707106781186548,                    0,    0.707106781186548},
           {  -0.707106781186548,   -0.707106781186548,                    0},
           {  -0.816496580927726,   -0.408248290463863,   -0.408248290463863},
           {   0.408248290463863,    0.816496580927726,   -0.408248290463863},
           {   0.408248290463863,   -0.408248290463863,    0.816496580927726}};

        int k,l,kk = 0;
        for (k = 0; k < 4; k++) 
        {
            VECTOR_COPY(burg[k],burgRef[k]);
  
#if 1
            // to choose the {110} planes
            for (l = 0; l < 3; l++) 
               VECTOR_COPY(glidePlane[k][l],glideRef[kk+l]);
            
#else            
            // to chose the {112} planes
            for (l = 3; l < 6; l++) 
            {
               int ll = l - 3;
               VECTOR_COPY(glidePlane[k][ll],glideRef[kk+l]);
            }
#endif
            kk += 6;
        }


#if 0// debugging
        for (k = 0; k < 4; k++) 
        {
           Print3("burg",burg[k]);
           for (l = 0; l < 3; l++) 
           {
              Print3("glide",glidePlane[k][l]);
              printf("b . n = %e\n",DotProduct(burg[k],glidePlane[k][l]));
           }
           printf("\n");
        }
        exit(0);
#endif

        id0 = 0;
        inData->nodeCount = 0;

/*
 *      Create the specified number of chains.
 */
        for (ic = 0; ic < numChains; ic++) {

            np = 3;
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

            burgIndex = ic%4;        // 4 different Burgers vectors
            gpIndex = (ic / 4) % 3;  // 3 : number of glide planes per Burgers vectors

/*
 *          Set up 3 initial points for the line.  Point 1 is a base position
 *          at a random location, point 0 is in the negative direction along
 *          the line and point 2 is in the positive direction along the line.
 */
            xp[1] = (randm(&seed)-0.5)*cubeSize;
            yp[1] = (randm(&seed)-0.5)*cubeSize;
            zp[1] = (randm(&seed)-0.5)*cubeSize;

            xp[0] = xp[1] + (cubeSize * burg[burgIndex][X]);
            yp[0] = yp[1] + (cubeSize * burg[burgIndex][Y]);
            zp[0] = zp[1] + (cubeSize * burg[burgIndex][Z]);

            xp[2] = xp[1] - (cubeSize * burg[burgIndex][X]);
            yp[2] = yp[1] - (cubeSize * burg[burgIndex][Y]);
            zp[2] = zp[1] - (cubeSize * burg[burgIndex][Z]);

/*
 *          Loop over the points and set up the nodes, link them to 
 *          the neighbor nodes, etc.
 */
            for (ip = 0; ip < np; ip++) {

                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];

                SET_CONSTRAINTS(node->constraint, UNCONSTRAINED);
                node->myTag.domainID = dislocType;
                node->myTag.index = ip + id0 + param->nodeCount;

                AllocNodeArms(node, 2);

                if ((nbr1Index = ip + 1) >= np) nbr1Index = 0;
                if ((nbr2Index = ip - 1) < 0) nbr2Index = np - 1;

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = id0 + nbr1Index + param->nodeCount;
                node->burgX[0] = burg[burgIndex][0];
                node->burgY[0] = burg[burgIndex][1];
                node->burgZ[0] = burg[burgIndex][2];
                node->nx[0] = glidePlane[burgIndex][gpIndex][X];
                node->ny[0] = glidePlane[burgIndex][gpIndex][Y];
                node->nz[0] = glidePlane[burgIndex][gpIndex][Z];
            
                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = id0 + nbr2Index + param->nodeCount;
                node->burgX[1] = -burg[burgIndex][0];
                node->burgY[1] = -burg[burgIndex][1];
                node->burgZ[1] = -burg[burgIndex][2];
                node->nx[1] = glidePlane[burgIndex][gpIndex][X];
                node->ny[1] = glidePlane[burgIndex][gpIndex][Y];
                node->nz[1] = glidePlane[burgIndex][gpIndex][Z];

            }

            id0 += np;
/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths before
 *          we write the nodal data to the file
 */
            lastBlock = (ic == numChains - 1);

            if (lastBlock || inArgs->useSegmentedDataFile) {

                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);

                id0 = 0;
                startRemeshIndex = 0;
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFCCScrewConfig
 *      Description:  Create an initial configuration consisting of a
 *                    series of FCC screw dislocations in a periodic
 *                    simulation space.
 *
 *      Arguments:
 *          IN:  inArgs       Structure containing values associated with the
 *                            command line arguments.
 *          OUT: totDislocLen Total length of all dislocations inserted into
 *                            the configuration.
 *
 *      Values used from the <inArgs> structure:
 *
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of dislocation chains to create
 *          seed        Seed value for random number generator
 *          type        Type of dislocation configuration being generated.
 *
 *-------------------------------------------------------------------------*/
void CreateFCCScrewConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                          real8 *totDislocLen)
{
        int      ic, np, ip, id0;
        int      burgIndex, gpIndex, newNodeIndex, lastBlock;
        int      nbr1Index, nbr2Index;
        int      numChains, seed, dislocType;
        int      startRemeshIndex = 0;
        real8    xp[3], yp[3], zp[3], cubeSize;
        real8    burg[6][3], glidePlane[6][2][3];
        Param_t  *param;
        Node_t   *node;

/*
 *      Grap the command-line argument values needed for this function
 */
        cubeSize   = (real8)inArgs->cubeLength;
        numChains  = inArgs->numChains;
        seed       = inArgs->seed;
        dislocType = inArgs->type;

        param = inData->param;

        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be > 0.\n",
                  "CreateScrewConfig", numChains);
        }

/*
 *      The nodal data is generated in 1 or more blocks.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes plus an offset to account for any previous blocks.
 *      The InitRemesh() calls also accept the dislocation type to set the
 *      tag.domainID value for any new nodes it adds to the node array.
 */

/*
 *      Set up an array of valid burgers vectors
 */
        inData->nburg = 4;

        static real8 burgRef[6][3] = {
           { 7.0710678118e-01,  7.0710678118e-01,           0.0e+00},
           { 7.0710678118e-01, -7.0710678118e-01,           0.0e+00},
           { 7.0710678118e-01,           0.0e+00,  7.0710678118e-01},
           { 7.0710678118e-01,           0.0e+00, -7.0710678118e-01},
           {          0.0e+00,  7.0710678118e-01,  7.0710678118e-01},
           {          0.0e+00,  7.0710678118e-01, -7.0710678118e-01}};

        static real8 glideRef[12][3] = {             
            /* for burgRef[0][*] */
           { 5.77350269189626e-01, -5.77350269189626e-01,  5.77350269189626e-01}, 
           {-5.77350269189626e-01,  5.77350269189626e-01,  5.77350269189626e-01}, 
             /* for burgRef[1][*] */
           { 5.77350269189626e-01,  5.77350269189626e-01,  5.77350269189626e-01}, 
           { 5.77350269189626e-01,  5.77350269189626e-01, -5.77350269189626e-01}, 
             /* for burgRef[2][*] */
           { 5.77350269189626e-01,  5.77350269189626e-01, -5.77350269189626e-01}, 
           {-5.77350269189626e-01,  5.77350269189626e-01,  5.77350269189626e-01}, 
             /* for burgRef[3][*] */
           { 5.77350269189626e-01,  5.77350269189626e-01,  5.77350269189626e-01}, 
           { 5.77350269189626e-01, -5.77350269189626e-01,  5.77350269189626e-01}, 
             /* for burgRef[4][*] */
           { 5.77350269189626e-01,  5.77350269189626e-01, -5.77350269189626e-01}, 
           { 5.77350269189626e-01, -5.77350269189626e-01,  5.77350269189626e-01}, 
             /* for burgRef[5][*] */
           { 5.77350269189626e-01,  5.77350269189626e-01,  5.77350269189626e-01}, 
           {-5.77350269189626e-01,  5.77350269189626e-01,  5.77350269189626e-01}};

        int k,l,kk = 0;
        for (k = 0; k < 6; k++) 
        {
            VECTOR_COPY(burg[k],burgRef[k]);
  
            for (l = 0; l < 2; l++) 
               VECTOR_COPY(glidePlane[k][l],glideRef[kk+l]);
            
            kk += 2;
        }


#if 0// debugging
        for (k = 0; k < 6; k++) 
        {
           Print3("burg",burg[k]);
           for (l = 0; l < 2; l++) 
           {
              Print3("glide",glidePlane[k][l]);
              printf("b . n = %e\n",DotProduct(burg[k],glidePlane[k][l]));
           }
           printf("\n");
        }
        exit(0);
#endif

        id0 = 0;
        inData->nodeCount = 0;

/*
 *      Create the specified number of chains.
 */
        for (ic = 0; ic < numChains; ic++) {

            np = 3;
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

            burgIndex = ic%6;        // 6 : different Burgers vectors
            gpIndex = (ic / 6) % 2;  // 3 : number of glide planes per Burgers vectors

/*
 *          Set up 3 initial points for the line.  Point 1 is a base position
 *          at a random location, point 0 is in the negative direction along
 *          the line and point 2 is in the positive direction along the line.
 */
            xp[1] = (randm(&seed)-0.5)*cubeSize;
            yp[1] = (randm(&seed)-0.5)*cubeSize;
            zp[1] = (randm(&seed)-0.5)*cubeSize;

            xp[0] = xp[1] + (cubeSize * burg[burgIndex][X]);
            yp[0] = yp[1] + (cubeSize * burg[burgIndex][Y]);
            zp[0] = zp[1] + (cubeSize * burg[burgIndex][Z]);

            xp[2] = xp[1] - (cubeSize * burg[burgIndex][X]);
            yp[2] = yp[1] - (cubeSize * burg[burgIndex][Y]);
            zp[2] = zp[1] - (cubeSize * burg[burgIndex][Z]);

/*
 *          Loop over the points and set up the nodes, link them to 
 *          the neighbor nodes, etc.
 */
            for (ip = 0; ip < np; ip++) {

                node = &inData->node[ip+id0];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];

                SET_CONSTRAINTS(node->constraint, UNCONSTRAINED);
                node->myTag.domainID = dislocType;
                node->myTag.index = ip + id0 + param->nodeCount;

                AllocNodeArms(node, 2);

                if ((nbr1Index = ip + 1) >= np) nbr1Index = 0;
                if ((nbr2Index = ip - 1) < 0) nbr2Index = np - 1;

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = id0 + nbr1Index + param->nodeCount;
                node->burgX[0] = burg[burgIndex][0];
                node->burgY[0] = burg[burgIndex][1];
                node->burgZ[0] = burg[burgIndex][2];
                node->nx[0] = glidePlane[burgIndex][gpIndex][X];
                node->ny[0] = glidePlane[burgIndex][gpIndex][Y];
                node->nz[0] = glidePlane[burgIndex][gpIndex][Z];
            
                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = id0 + nbr2Index + param->nodeCount;
                node->burgX[1] = -burg[burgIndex][0];
                node->burgY[1] = -burg[burgIndex][1];
                node->burgZ[1] = -burg[burgIndex][2];
                node->nx[1] = glidePlane[burgIndex][gpIndex][X];
                node->ny[1] = glidePlane[burgIndex][gpIndex][Y];
                node->nz[1] = glidePlane[burgIndex][gpIndex][Z];

            }

            id0 += np;
/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths before
 *          we write the nodal data to the file
 */
            lastBlock = (ic == numChains - 1);

            if (lastBlock || inArgs->useSegmentedDataFile) {

                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);

                id0 = 0;
                startRemeshIndex = 0;
            }
        }

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFCCIrradConfig
 *      Description:  Generate an initial FCC configuration consisting of
 *                    a combination of line dislocations and Frank Sessile
 *                    Hexagonal Interstitial Loops.
 *                    (Masato Hiratani)
 *
 *      Arguments:
 *          IN:  inArgs       Structure containing values associated with the
 *                            command line arguments.
 *          OUT: totDislocLen Total length of all dislocations inserted into
 *                            the configuration.
 *
 *      Values used from the <inArgs> structure:
 *
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of chains to create
 *          numLoops    Number of hexagonal loops to create
 *          hexl        Size of the hexagonal interstitial loops
 *          seed        Seed value for random number generator
 *          type        Type of dislocation configuration being generated.
 *
 *-------------------------------------------------------------------------*/
void CreateFCCIrradConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                          real8 *totDislocLen)
{
        int      ic, np, ip, i;
        int      nplane, indp, indb, indf, inds, indr;
        int      newNodeIndex, offSet, lastBlock;
        int      numChains, numLoops, seed, dislocType;
        int      startRemeshIndex = 0;
        real8    hexl;
        real8    inv6, invsq2, sq2over3;
        real8    xp[6], yp[6], zp[6], cubeSize;
        real8    tnx[4], tny[4], tnz[4], burg[12][3];
        Node_t   *node;
        Param_t  *param;

        param = inData->param;

/*
 *      Grap the command-line argument values needed for this function
 */
        cubeSize   = (real8)inArgs->cubeLength;
        numLoops   = inArgs->numLoops;
        numChains  = inArgs->numChains;
        seed       = inArgs->seed;
        hexl       = inArgs->hexl;
        dislocType = inArgs->type;

/*
 *      Quick sanity checks...
 */
        if ((numLoops <= 0) && (numChains <= 0)) {
            Fatal("%s: numLoops = %d, numChains = %d, but at least one of\n"
                  "these must be greater than zero!\n",
                  "CreateFCCIrradConfig", numLoops, numChains);
        }

        if (numLoops < 0) {
            Fatal("%s: numLoops is %d, but must be > 0.\n",
                  "CreateFCCIrradConfig", numLoops);
        }

        if (numChains < 0) {
            Fatal("%s: numChains is %d, but must be > 0.\n",
                  "CreateFCCIrradConfig", numChains);
        }

       
/*
 *      The nodal data is generated in 1 or more blocks.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes plus an offset to account for any previous blocks.
 *      The InitRemesh() calls also accept the dislocation type to set the
 *      tag.domainID value for any new nodes it adds to the node array.
 */
        inv6     = 0.1666666;
        invsq2   = 0.70710678;
        sq2over3 = 0.81649658;

        inData->nburg = 12;

/*
 *      alpha, beta, gamma, and delta plane normals
 */
        nplane = 4;
 
        tnx[0] = -1;
        tny[0] =  1;
        tnz[0] = -1;

        tnx[1] =  1;
        tny[1] = -1;
        tnz[1] = -1; 

        tnx[2] = -1;
        tny[2] = -1;
        tnz[2] =  1;

        tnx[3] =  1;
        tny[3] =  1;
        tnz[3] =  1; 

/*
 *      BV is parallel to the plane normal i.e. FS 1/3<111>
 *      BV in counter-clockwise on each plane
 */
        for (i = 0; i < nplane; i++) {

            burg[3*i][0] = 0;
            burg[3*i][1] = invsq2*tny[i];
            burg[3*i][2] = -invsq2*tnz[i];

            burg[1+3*i][0] = -invsq2*tnx[i];
            burg[1+3*i][1] = 0;
            burg[1+3*i][2] = invsq2*tnz[i];

            burg[2+3*i][0] = invsq2*tnx[i];
            burg[2+3*i][1] = -invsq2*tny[i];
            burg[2+3*i][2] = 0;

            Normalize(&tnx[i], &tny[i], &tnz[i]);
        }

/*
 *      We use <offSet> to track the position of the node within the current
 *      block of nodes in <inData> that corresponds to the first node in
 *      the current chain/hex-loop.
 */
        offSet = 0;

        inData->nodeCount = 0;

        for (ic = 0; ic < numChains; ic++) {
            int lastChain;

            np = 3;    /* number of points along 1 chain */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          plane normal cycle 4, BV cycle 12, 60deg line sense 24,
 *          reflection 48
 */
            indp = ic%4;                /* index of plane */
            indb = indp*3+ic%3;         /* index of BV */
            indf = ((ic-ic%12)/12)%2;   /* index of alternative line sense */
            inds = indp*3+(ic+indf+1)%3;        /* index of line sense */
            indr = 1-2*(((ic-ic%24)/24)%2);     /* sign of reflection */

            xp[0] = (randm(&seed)-0.5)*cubeSize;
            yp[0] = (randm(&seed)-0.5)*cubeSize;
            zp[0] = (randm(&seed)-0.5)*cubeSize;

/*
 *          shift along neighboring 60 degree BV
 */
            xp[1] = xp[0] + indr * cubeSize * burg[inds][0];
            yp[1] = yp[0] + indr * cubeSize * burg[inds][1];
            zp[1] = zp[0] + indr * cubeSize * burg[inds][2];

            xp[2] = xp[1] + indr * cubeSize * burg[inds][0];
            yp[2] = yp[1] + indr * cubeSize * burg[inds][1];
            zp[2] = zp[1] + indr * cubeSize * burg[inds][2];

/*
 *          PBC
 */        
            for (i = 0; i < 3; i++) {
                if (xp[i] < -cubeSize/2) xp[i] += cubeSize;
                if (yp[i] < -cubeSize/2) yp[i] += cubeSize;
                if (zp[i] < -cubeSize/2) zp[i] += cubeSize;
            }

            for (ip = 0; ip < np; ip++) {
                node = &inData->node[ip+offSet];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];
                SET_CONSTRAINTS(node->constraint, PINNED_NODE);
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+offSet+param->nodeCount;

                AllocNodeArms(node, 2);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = (ip-1+np)%np+offSet+param->nodeCount;
                node->burgX[0] = burg[indb][0];
                node->burgY[0] = burg[indb][1];
                node->burgZ[0] = burg[indb][2];
                node->nx[0] = tnx[indp];
                node->ny[0] = tny[indp];
                node->nz[0] = tnz[indp];

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = (ip+1+np)%np+offSet+param->nodeCount;
                node->burgX[1] = -burg[indb][0];
                node->burgY[1] = -burg[indb][1];
                node->burgZ[1] = -burg[indb][2];
                node->nx[1] = tnx[indp];
                node->ny[1] = tny[indp];
                node->nz[1] = tnz[indp];
            }

            offSet += np;

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          Then, if the count of nodes currently contained in memoru
 *          exceeds the threshhold write the current block of nodal data
 *          to the file.
 */

            lastChain = (ic == (numChains - 1));
            lastBlock = (lastChain && (numLoops == 0));

            if (inArgs->useSegmentedDataFile || lastBlock) {
                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);
                startRemeshIndex = 0;
                offSet = 0;
            }

        } /* end of chains */

/*
 *      Place hexagonal loops
 */
        for (ic = 0; ic < numLoops; ic++) {
            np = 6;     /* creation of hexagonal loop */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          plane normal cycle 4, BV cycle 12, 60deg line sense 24,
 *          reflection 48
 */
            indp = ic%4;              /* index of plane */
            indb = indp*3;            /* index of BV */
            indf = ((ic-ic%12)/12)%2; /* index of alternativeline sense */
            inds = indp*3+(ic+indf+1)%3;    /* index of line sense */
            indr = 1-2*(((ic-ic%24)/24)%2); /* sign of reflection */

            xp[0] = (randm(&seed)-0.5)*cubeSize;
            yp[0] = (randm(&seed)-0.5)*cubeSize;
            zp[0] = (randm(&seed)-0.5)*cubeSize;

/*
 *          shift along neighboring 60 degree BV
 */
            xp[1] = xp[0] - inv6*hexl*burg[indb][0];
            yp[1] = yp[0] - inv6*hexl*burg[indb][1];
            zp[1] = zp[0] - inv6*hexl*burg[indb][2];
            xp[2] = xp[1] + inv6*hexl*burg[indb+2][0];
            yp[2] = yp[1] + inv6*hexl*burg[indb+2][1];
            zp[2] = zp[1] + inv6*hexl*burg[indb+2][2];
            xp[3] = xp[2] - inv6*hexl*burg[indb+1][0];
            yp[3] = yp[2] - inv6*hexl*burg[indb+1][1];
            zp[3] = zp[2] - inv6*hexl*burg[indb+1][2];
            xp[4] = xp[3] + inv6*hexl*burg[indb][0];
            yp[4] = yp[3] + inv6*hexl*burg[indb][1];
            zp[4] = zp[3] + inv6*hexl*burg[indb][2];
            xp[5] = xp[4] - inv6*hexl*burg[indb+2][0];
            yp[5] = yp[4] - inv6*hexl*burg[indb+2][1];
            zp[5] = zp[4] - inv6*hexl*burg[indb+2][2];
      
/*
 *          PBC
 */
            for (i = 0; i < np; i++) {
                if (xp[i]<-cubeSize/2) xp[i] += cubeSize;
                if (yp[i]<-cubeSize/2) yp[i] += cubeSize;
                if (zp[i]<-cubeSize/2) zp[i] += cubeSize;
            }

            for (ip = 0; ip < np; ip++) {
                node = &inData->node[ip+offSet];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];
                SET_CONSTRAINTS(node->constraint, PINNED_NODE);
                node->myTag.domainID = dislocType;
                node->myTag.index = ip+offSet+param->nodeCount;

                AllocNodeArms(node, 2);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = (ip-1+np)%np+offSet+param->nodeCount;
                node->burgX[0] = sq2over3*tnx[indp];
                node->burgY[0] = sq2over3*tny[indp];
                node->burgZ[0] = sq2over3*tnz[indp];
                node->nx[0] = tnx[indp];
                node->ny[0] = tny[indp];
                node->nz[0] = tnz[indp];

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = (ip+1+np)%np+offSet+param->nodeCount;
                node->burgX[1] = -sq2over3*tnx[indp];
                node->burgY[1] = -sq2over3*tny[indp];
                node->burgZ[1] = -sq2over3*tnz[indp];
                node->nx[1] = tnx[indp];
                node->ny[1] = tny[indp];
                node->nz[1] = tnz[indp];
            }   

            offSet += np;

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the block of nodal data to the file.
 */
            lastBlock = (ic == (numLoops - 1));

            if (inArgs->useSegmentedDataFile || lastBlock) {
                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);
                startRemeshIndex = 0;
                offSet = 0;
            }

        }  /* end for (ic = 0; ic < numLoops; ...) */

        return; 
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFCCConfig
 *      Description:  Generates an initial configuration consisting of
 *                    random FCC dislocations.  (Masato Hiratani)
 *
 *      Arguments:
 *          IN:  inArgs       Structure containing values associated with the
 *                            command line arguments.
 *          OUT: totDislocLen Total length of all dislocations inserted into
 *                            the configuration.
 *
 *      Values used from the <inArgs> structure:
 *
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numChains   Number of chains to create
 *          seed        Seed value for random number generator
 *          type        Type of dislocation configuration being generated.
 *
 *-------------------------------------------------------------------------*/
void CreateFCCConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                     real8 *totDislocLen)
{
        int      ic, np, ip, i, offSet;
        int      nplane, indp, indb, indf, inds, indr;
        int      newNodeIndex, lastBlock;
        int      numChains, seed, dislocType;
        int      startRemeshIndex = 0;
        real8    invsq2;
        real8    xp[3], yp[3], zp[3];
        real8    tnx[4], tny[4], tnz[4], burg[12][3], cubeSize;
        Param_t  *param;
        Node_t   *node;

        param = inData->param;
        invsq2 = 0.70710678118;

/*
 *      Grap the command-line argument values needed for this function
 */
        cubeSize   = (real8)inArgs->cubeLength;
        numChains  = inArgs->numChains;
        seed       = inArgs->seed;
        dislocType = inArgs->type;

        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be a multiple of 12.\n",
                  "CreateFCCConfig", numChains);
        }
     
/*
 *      The nodal data is generated in 1 or more blocks.  All nodes will
 *      be assigned a tag.domainID equal to the dislocation type, and
 *      the index for every node will be the node's actual index into the
 *      block of nodes plus an offset to account for any previous blocks.
 *      The InitRemesh() calls also accept the dislocation type to set the
 *      tag.domainID value for any new nodes it adds to the node array.
 */
        inData->nburg = 12;

/*
 *      alpha, beta, gamma, and delta planes
 */
        nplane =  4;

        tnx[0] = -1;
        tny[0] =  1;
        tnz[0] = -1;

        tnx[1] =  1;
        tny[1] = -1;
        tnz[1] = -1; 

        tnx[2] = -1;
        tny[2] = -1;
        tnz[2] =  1;

        tnx[3] =  1;
        tny[3] =  1;
        tnz[3] =  1; 

/*
 *      BV in counter-clockwise on each plane
 */
        for (i = 0; i < nplane; i++) {

            burg[3*i][0] = 0; 
            burg[3*i][1] = invsq2*tny[i];
            burg[3*i][2] = -invsq2*tnz[i]; 

            burg[1+3*i][0] = -invsq2*tnx[i];
            burg[1+3*i][1] = 0;
            burg[1+3*i][2] = invsq2*tnz[i]; 

            burg[2+3*i][0] = invsq2*tnx[i];
            burg[2+3*i][1] = -invsq2*tny[i];
            burg[2+3*i][2] = 0; 

            Normalize(&tnx[i],&tny[i],&tnz[i]);   
        }
    
/*
 *      We use <offSet> to track the position of the node within the current
 *      block of nodes in <inData> that corresponds to the first node in
 *      the current chain.
 */
        offSet = 0;

        inData->nodeCount = 0;

        for (ic = 0; ic < numChains; ic++) {

            np = 3;    /* number of points along 1 chain */

            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;
            inData->node = (Node_t *) realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          plane normal cycle 4, BV cycle 12, 60deg line sense 24,
 *          reflection 48
 */
            indp = ic%4;              /* index of plane */
            indb = indp*3+ic%3;       /* index of BV */
            indf = ((ic-ic%12)/12)%2; /* index of alternative line sense */
            inds = indp*3+(ic+indf+1)%3;    /* index of line sense */
            indr = 1-2*(((ic-ic%24)/24)%2); /* sign of reflection */

            xp[0] = (randm(&seed)-0.5)*cubeSize;
            yp[0] = (randm(&seed)-0.5)*cubeSize;
            zp[0] = (randm(&seed)-0.5)*cubeSize;

/*
 *          shift along neighboring 60 degree BV
 */
            xp[1] = xp[0] + indr*cubeSize*burg[inds][0];
            yp[1] = yp[0] + indr*cubeSize*burg[inds][1];
            zp[1] = zp[0] + indr*cubeSize*burg[inds][2];

            xp[2] = xp[1] + indr*cubeSize*burg[inds][0];
            yp[2] = yp[1] + indr*cubeSize*burg[inds][1];
            zp[2] = zp[1] + indr*cubeSize*burg[inds][2];

/*
 *          PBC
 */        
            for (i = 0; i < 3; i++) {
                if (xp[i]<-cubeSize/2) xp[i] += cubeSize;
                if (yp[i]<-cubeSize/2) yp[i] += cubeSize;
                if (zp[i]<-cubeSize/2) zp[i] += cubeSize;
            }

/*
 *          Define the three initial nodes that will make up
 *          this dislocation chain.
 */
            for (ip = 0; ip < np; ip++) {

                node = &inData->node[ip+offSet];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];
                SET_CONSTRAINTS(node->constraint, UNCONSTRAINED);
                node->myTag.domainID = dislocType;
                node->myTag.index = ip + offSet + param->nodeCount;

                AllocNodeArms(node, 2);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = (ip-1+np)%np+offSet+param->nodeCount;
                node->burgX[0] = burg[indb][0];
                node->burgY[0] = burg[indb][1];
                node->burgZ[0] = burg[indb][2];
                node->nx[0] = tnx[indp]; 
                node->ny[0] = tny[indp];
                node->nz[0] = tnz[indp];

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = (ip+1+np)%np+offSet+param->nodeCount;
                node->burgX[1] = -burg[indb][0];
                node->burgY[1] = -burg[indb][1];
                node->burgZ[1] = -burg[indb][2];
                node->nx[1] = tnx[indp]; 
                node->ny[1] = tny[indp];
                node->nz[1] = tnz[indp];
            }

            offSet += np;

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          When we've generated the nodal data for the final chain,
 *          write the block of nodal data to the file.
 */
            lastBlock = (ic == (numChains - 1));

            if (inArgs->useSegmentedDataFile || lastBlock) {
                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);
                startRemeshIndex = 0;
                offSet = 0;
            }

        } /* loop over chains */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateHCPLoop
 *      Description:  Generate an intial configuration consisting of
 *                    dsilocation loops for HCP simulations.
 *
 *                    Note: Creates loops for 10 burgers vectors.  Loops
 *                          are octagonal for the first 9 burgers vectors
 *                          and haxagonal for the 10th.
 *
 *      Arguments:
 *          IN:  inArgs       Structure containing values associated with the
 *                            command line arguments.
 *          OUT: totDislocLen Total length of all dislocations inserted into
 *                            the configuration.
 *
 *      Values used from the <inArgs> structure:
 *
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numLoops    Number of loops to create
 *          radius      radius (in units of b) of the loop
 *          seed        Seed value for random number generator
 *          type        Type of dislocation configuration being generated.
 *
 *-------------------------------------------------------------------------*/
void CreateHCPLoop(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                   real8 *totDislocLen)
{
        int     i, j, m, itrx;
        int     id, loopIndex, burgIndex, newNodeIndex;
        int     minBurgIndex, maxBurgIndex, lastBlock, numSegs;
        int     startRemeshIndex = 0;
        int     startIndex, endIndex;
        int     numLoops, seed, dislocType;
        int     offSet;
        int     planesPerBurg[10], burgPlaneIndex[10];
        real8   radius;
        real8   cOa, cubeSize;
        real8   ux, uy, uz;
        real8   loopCtrX, loopCtrY, loopCtrZ;
        real8   plane[3];
        real8   b[10][3], p[39][3], edge[39][3];
        real8   tr[10][8][3];
        real8   A[3], B[3], C[3], delta[3], Ap[3], Bp[3], Cp[3];
        Param_t *param;
        Node_t  *node, *nbr1Node;

        param = inData->param;
        cOa = param->cOVERa;

/*
 *      Grap the command-line argument values needed for this function
 */
        cubeSize   = (real8)inArgs->cubeLength;
        numLoops   = inArgs->numLoops;
        radius     = inArgs->radius;
        seed       = inArgs->seed;
        dislocType = inArgs->type;

/*
 *      Points
 */
        A[0] =  0.0;
        A[1] =  0.0;
        A[2] =  0.0;

        B[0] = -0.5;
        B[1] =  sqrt(3.0) * 0.5;
        B[2] =  0.0;

        C[0] =  0.5;
        C[1] =  sqrt(3.0) * 0.5;
        C[2] =  0.0;

        delta[0] =  0.0;
        delta[1] =  sqrt(3.0)/3.0;
        delta[2] =  0.5 * cOa;

        Ap[0] =  0.0;
        Ap[1] =  0.0;
        Ap[2] =  cOa;

        Bp[0] = -0.5;
        Bp[1] =  sqrt(3.0) * 0.5;
        Bp[2] =  cOa;

        Cp[0] =  0.5;
        Cp[1] =  sqrt(3.0) * 0.5;
        Cp[2] =  cOa;


        minBurgIndex = 0;
        maxBurgIndex = 9;


/*
 *      Define burgers vectors, plane and edge components of the loops.
 *
 *      NOTE: planes are not in the same order nor of the same sign as
 *            in HCP_Util.c.  Also, burgers vectors here are NOT normalized!
 */

/*
 *      AB
 */
        b[0][0] = B[0] - A[0];
        b[0][1] = B[1] - A[1];
        b[0][2] = B[2] - A[2];

        GetPlaneNormFromPoints(C,A,B,plane);
        VECTOR_COPY(p[0], plane);

        GetPlaneNormFromPoints(Cp,A,B,plane);
        VECTOR_COPY(p[1], plane);

        GetPlaneNormFromPoints(A,B,Bp,plane);
        VECTOR_COPY(p[2], plane);

        GetPlaneNormFromPoints(C,Bp,Ap,plane);
        VECTOR_COPY(p[3], plane);

        planesPerBurg[0] = 4;
        burgPlaneIndex[0] = 0;

/*
 *      BC
 */
        b[1][0] = C[0] - B[0];
        b[1][1] = C[1] - B[1];
        b[1][2] = C[2] - B[2];

        GetPlaneNormFromPoints(C,A,B,plane);
        VECTOR_COPY(p[4], plane);

        GetPlaneNormFromPoints(A,Bp,Cp,plane);
        VECTOR_COPY(p[5], plane);

        GetPlaneNormFromPoints(Cp,C,B,plane);
        VECTOR_COPY(p[6], plane);

        GetPlaneNormFromPoints(Ap,C,B,plane);
        VECTOR_COPY(p[7], plane);

        planesPerBurg[1] = 4;
        burgPlaneIndex[1] = 4;

/*
 *      CA
 */
        b[2][0] = A[0] - C[0];
        b[2][1] = A[1] - C[1];
        b[2][2] = A[2] - C[2];

        GetPlaneNormFromPoints(A,C,B,plane);
        VECTOR_COPY(p[8], plane);

        GetPlaneNormFromPoints(Bp,A,C,plane);
        VECTOR_COPY(p[9], plane);

        GetPlaneNormFromPoints(A,C,Ap,plane);
        VECTOR_COPY(p[10], plane);

        GetPlaneNormFromPoints(B,Cp,Ap,plane);
        VECTOR_COPY(p[11], plane);

        planesPerBurg[2] = 4;
        burgPlaneIndex[2] = 8;

/*
 *      ABp
 */
        b[3][0] = Bp[0] - A[0];
        b[3][1] = Bp[1] - A[1];
        b[3][2] = Bp[2] - A[2];

        GetPlaneNormFromPoints(delta,A,Bp,plane);
        VECTOR_COPY(p[12], plane);

        GetPlaneNormFromPoints(A,Bp,Cp,plane);
        VECTOR_COPY(p[13], plane);

        GetPlaneNormFromPoints(A,B,Bp,plane);
        VECTOR_COPY(p[14], plane);

        GetPlaneNormFromPoints(Bp,A,C,plane);
        VECTOR_COPY(p[15], plane);

        planesPerBurg[3] = 4;
        burgPlaneIndex[3] = 12;

/*
 *      BCp
 */
        b[4][0] = Cp[0] - B[0];
        b[4][1] = Cp[1] - B[1];
        b[4][2] = Cp[2] - B[2];

        GetPlaneNormFromPoints(C,B,Cp,plane);
        VECTOR_COPY(p[16], plane);

        GetPlaneNormFromPoints(Cp,A,B,plane);
        VECTOR_COPY(p[17], plane);

        GetPlaneNormFromPoints(delta,B,Cp,plane);
        VECTOR_COPY(p[18], plane);

        GetPlaneNormFromPoints(B,Cp,Ap,plane);
        VECTOR_COPY(p[19], plane);

        planesPerBurg[4] = 4;
        burgPlaneIndex[4] = 16;

/*
 *      CAp
 */
        b[5][0] = Ap[0] - C[0];
        b[5][1] = Ap[1] - C[1];
        b[5][2] = Ap[2] - C[2];

        GetPlaneNormFromPoints(C,delta,Ap,plane);
        VECTOR_COPY(p[20], plane);

        GetPlaneNormFromPoints(Ap,C,B,plane);
        VECTOR_COPY(p[21], plane);

        GetPlaneNormFromPoints(C,A,Ap,plane);
        VECTOR_COPY(p[22], plane);

        GetPlaneNormFromPoints(Bp,C,Ap,plane);
        VECTOR_COPY(p[23], plane);

        planesPerBurg[5] = 4;
        burgPlaneIndex[5] = 20;

/*
 *      ApB
 */
        b[6][0] = B[0] - Ap[0];
        b[6][1] = B[1] - Ap[1];
        b[6][2] = B[2] - Ap[2];

        GetPlaneNormFromPoints(Ap,delta,B,plane);
        VECTOR_COPY(p[24], plane);

        GetPlaneNormFromPoints(Ap,C,B,plane);
        VECTOR_COPY(p[25], plane);

        GetPlaneNormFromPoints(Ap,B,Bp,plane);
        VECTOR_COPY(p[26], plane);

        GetPlaneNormFromPoints(B,Cp,Ap,plane);
        VECTOR_COPY(p[27], plane);

        planesPerBurg[6] = 4;
        burgPlaneIndex[6] = 24;

/*
 *      BpC
 */
        b[7][0] = C[0] - Bp[0];
        b[7][1] = C[1] - Bp[1];
        b[7][2] = C[2] - Bp[2];

        GetPlaneNormFromPoints(Bp,delta,C,plane);
        VECTOR_COPY(p[28], plane);

        GetPlaneNormFromPoints(C,Bp,Ap,plane);
        VECTOR_COPY(p[29], plane);

        GetPlaneNormFromPoints(Bp,Cp,C,plane);
        VECTOR_COPY(p[30], plane);

        GetPlaneNormFromPoints(Bp,C,A,plane);
        VECTOR_COPY(p[31], plane);

        planesPerBurg[7] = 4;
        burgPlaneIndex[7] = 28;

/*
 *      CpA
 */
        b[8][0] = A[0] - Cp[0];
        b[8][1] = A[1] - Cp[1];
        b[8][2] = A[2] - Cp[2];

        GetPlaneNormFromPoints(Cp,A,delta,plane);
        VECTOR_COPY(p[32], plane);

        GetPlaneNormFromPoints(A,Bp,Cp,plane);
        VECTOR_COPY(p[33], plane);

        GetPlaneNormFromPoints(Cp,A,Ap,plane);
        VECTOR_COPY(p[34], plane);

        GetPlaneNormFromPoints(Cp,B,A,plane);
        VECTOR_COPY(p[35], plane);

        planesPerBurg[8] = 4;
        burgPlaneIndex[8] = 32;

/*
 *      AAp
 */
        b[9][0] = Ap[0] - A[0];
        b[9][1] = Ap[1] - A[1];
        b[9][2] = Ap[2] - A[2];

        GetPlaneNormFromPoints(C,A,Ap,plane);
        VECTOR_COPY(p[36], plane);

        GetPlaneNormFromPoints(C,B,Cp,plane);
        VECTOR_COPY(p[37], plane);

        GetPlaneNormFromPoints(Ap,A,Bp,plane);
        VECTOR_COPY(p[38], plane);

        planesPerBurg[9] = 3;
        burgPlaneIndex[9] = 36;

/*
 *      Normalize glide planes.  The burgers vectors that go into the
 *      nodal data file shold NOT be normalized!
 */
        for (i = 0; i < 39; i++) {
            NormalizeVec(p[i]);
        }

/*
 *      Get edge directions
 */
        for (i = 0; i < 10; i++) {

            startIndex = burgPlaneIndex[i];
            endIndex = burgPlaneIndex[i] + planesPerBurg[i];

            for (j = startIndex; j < endIndex; j++) {

                cross(p[j], b[i], edge[j]);
                NormalizeVec(edge[j]);
            }
        }

/*
 *      Define the translation vectors as 4 X 6 X 3 array
 *
 *          tr[ ][ ][ ] - translation vector
 *             |  |  |
 *             |  |  x,y,z
 *             |  |
 *             |  8 max  nodes on the loop
 *             |
 *            10  type burgers vectors
 */

        for (itrx = 0; itrx < 10; itrx++) {
            int     itry;
            real8   XX[8][3];

            VECTOR_ZERO(XX[0]);

            m = 1;
            startIndex = burgPlaneIndex[itrx];
            endIndex   = burgPlaneIndex[itrx] + planesPerBurg[itrx];

            for (j = startIndex; j < endIndex; j++) {
                XX[m][0] = XX[m-1][0] + edge[j][0];
                XX[m][1] = XX[m-1][1] + edge[j][1];
                XX[m][2] = XX[m-1][2] + edge[j][2];
                m += 1;
            }

            for (j = startIndex; j < (endIndex-1); j++) {
                XX[m][0] = XX[m-1][0] - edge[j][0];
                XX[m][1] = XX[m-1][1] - edge[j][1];
                XX[m][2] = XX[m-1][2] - edge[j][2];
                m += 1;
            }

            for (itry = 0; itry < m; itry++) {
                VECTOR_COPY(tr[itrx][itry], XX[itry]);
            }
        }

/*
 *      We use <offSet> to track the position of the node within the current
 *      block of nodes in <inData> that corresponds to the first node in
 *      the current chain.
 */
        offSet = 0;

/*
 *      Create one loop at a time, cycling through burgers vectors as we go.
 */
        inData->nodeCount = 0;
        burgIndex = maxBurgIndex;

        for (loopIndex = 0; loopIndex < numLoops; loopIndex++) {

            if (++burgIndex > maxBurgIndex) {
                burgIndex = minBurgIndex;
            }

/*
 *          Increase the size of the node array enough to hold
 *          all the new nodes created for this loop.
 *
 *          NOTE: The loops for the first 9 burgers vectors are
 *                octagonal (8 segments/vertices), but the loop
 *                for the last burgers vector is hexagonal (6
 *                segments/vertices).
 */
            numSegs = ((burgIndex % 10) < 9) ? 8 : 6;

            newNodeIndex = inData->nodeCount;
            inData->nodeCount += numSegs;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numSegs);

            loopCtrX = (randm(&seed)-0.5) * cubeSize;
            loopCtrY = (randm(&seed)-0.5) * cubeSize;
            loopCtrZ = (randm(&seed)-0.5) * cubeSize;

            for (id = offSet; id < (offSet + numSegs); id++) {
                int   dIndex;

                node =  &inData->node[id];

/*
 *              Pick the index for the direction component of the tr array
 */
                dIndex = (id - offSet) % numSegs;

                node->x = loopCtrX + radius * tr[burgIndex][dIndex][0];
                node->y = loopCtrY + radius * tr[burgIndex][dIndex][1];
                node->z = loopCtrZ + radius * tr[burgIndex][dIndex][2];

/*
 *              Set up the node
 */
                SET_CONSTRAINTS(node->constraint, UNCONSTRAINED);
                node->myTag.domainID = dislocType;
                node->myTag.index = id + param->nodeCount;

                AllocNodeArms(node, 2);

                node->burgX[0] =  b[burgIndex][X];
                node->burgY[0] =  b[burgIndex][Y];
                node->burgZ[0] =  b[burgIndex][Z];

                node->burgX[1] = -b[burgIndex][X];
                node->burgY[1] = -b[burgIndex][Y];
                node->burgZ[1] = -b[burgIndex][Z];

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index =
                        offSet + (((id-offSet)+1)%numSegs) + param->nodeCount;

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index =
                        offSet + param->nodeCount +
                        (((id-offSet)+numSegs-1)%numSegs);
            }

/*
 *          Need to set the glide plane normal for each segment.
 *          Couldn't do it above because we didn't have all the
 *          neighbor node positions until all the nodes in the
 *          loop were created.
 */
            for (id = offSet; id < (offSet + numSegs); id++) {
                int     nbr1Index;
                real8   vec1[3], vec2[3], tmpVec[3];

                node = &inData->node[id];

                nbr1Index = node->nbrTag[1].index - param->nodeCount;
                nbr1Node = &inData->node[nbr1Index];

                ux = nbr1Node->x - node->x;
                uy = nbr1Node->y - node->y;
                uz = nbr1Node->z - node->z;

                Normalize(&ux,&uy,&uz);

/*
 *              l cross b gives normal vector for each segmenta, but
 *              in order to avoid some round-off issues, explicitly set
 *              any components of the glide plane normal are near-zero
 *              to be zero, then renormalize.
 */
                vec1[0] = ux;
                vec1[1] = uy;
                vec1[2] = uz;

                vec2[0] =  node->burgX[1];
                vec2[1] =  node->burgY[1];
                vec2[2] =  node->burgZ[1];

                NormalizedCrossVector(vec1, vec2, tmpVec);

                if (fabs(tmpVec[0]) < 1.0e-03) tmpVec[0] = 0.0;
                if (fabs(tmpVec[1]) < 1.0e-03) tmpVec[1] = 0.0;
                if (fabs(tmpVec[2]) < 1.0e-03) tmpVec[2] = 0.0;

                NormalizeVec(tmpVec);

                node->nx[1] = tmpVec[0];
                node->ny[1] = tmpVec[1];
                node->nz[1] = tmpVec[2];

                nbr1Node->nx[0] = tmpVec[0];
                nbr1Node->ny[0] = tmpVec[1];
                nbr1Node->nz[0] = tmpVec[2];

            }  /*  for (id = nextNode; ...) */

            offSet += numSegs;

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          Write blocks of nodal data to the data files as needed.
 */
            lastBlock = (loopIndex == (numLoops-1));

            if (inArgs->useSegmentedDataFile || lastBlock) {
                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);
                startRemeshIndex = 0;
                offSet = 0;
            }

        }  /* for (loopIndex = 0; ...) */

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFCCPrismaticLoop
 *      Description:  Generate an intial configuration consisting of
 *                    dislocation loops for FCC simulations.
 * 
 *      Author:       S. Aubry
 *
 *                    Note: Creates loops for 6 burgers vectors.  Loops
 *                          are squares because two glides planes per burgers
 *                          vectors in FCC crystals.
 *
 *      Arguments:
 *          IN:  inArgs       Structure containing values associated with the
 *                            command line arguments.
 *          OUT: totDislocLen Total length of all dislocations inserted into
 *                            the configuration.
 *
 *      Values used from the <inArgs> structure:
 *
 *          cubeLength  Length of cubic problem space in a single
 *                      dimension (units of b)
 *          numLoops    Number of loops to create
 *          radius      radius (in units of b) of the loop
 *          seed        Seed value for random number generator
 *          type        Type of dislocation configuration being generated.
 *
 *-------------------------------------------------------------------------*/
void CreateFCCPrismaticLoop(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                   real8 *totDislocLen)
{
        int     i, j, m, itrx;
        int     id, loopIndex, burgIndex, newNodeIndex;
        int     minBurgIndex, maxBurgIndex, lastBlock, numSegs;
        int     startRemeshIndex = 0;
        int     numLoops, seed, dislocType;
        int     offSet;
        real8   radius;
        real8   cubeSize;
        real8   ux, uy, uz;
        real8   loopCtrX, loopCtrY, loopCtrZ;
        real8   b[6][3], p[2][6][3], e[2][6][3];
        real8   tr[6][4][3];
        Param_t *param;
        Node_t  *node, *nbr1Node;

        param = inData->param;

/*
 *      Define Burgers vectors, planes and edges of the loops
 */
        b[0][X]= 7.0710678118e-01;  b[0][Y]= 7.0710678118e-01;   b[0][Z]= 0.0e+00;
        b[1][X]= 7.0710678118e-01;  b[1][Y]=-7.0710678118e-01;   b[1][Z]= 0.0e+00;
        b[2][X]= 7.0710678118e-01;  b[2][Y]= 0.0e+00;            b[2][Z]= 7.0710678118e-01;
        b[3][X]= 7.0710678118e-01;  b[3][Y]= 0.0e+00;            b[3][Z]=-7.0710678118e-01;
        b[4][X]=          0.0e+00;  b[4][Y]= 7.0710678118e-01;   b[4][Z]= 7.0710678118e-01;
        b[5][X]=          0.0e+00;  b[5][Y]= 7.0710678118e-01;   b[5][Z]=-7.0710678118e-01;
        
        
        p[0][0][X]= 5.7735026919e-01; p[0][0][Y]=-5.7735026919e-01; p[0][0][Z]= 5.7735026919e-01;
        p[1][0][X]=-5.7735026919e-01; p[1][0][Y]= 5.7735026919e-01; p[1][0][Z]= 5.7735026919e-01;
        
        p[0][1][X]= 5.7735026919e-01; p[0][1][Y]= 5.7735026919e-01; p[0][1][Z]= 5.7735026919e-01;
        p[1][1][X]= 5.7735026919e-01; p[1][1][Y]= 5.7735026919e-01; p[1][1][Z]=-5.7735026919e-01;
        
        p[0][2][X]= 5.7735026919e-01; p[0][2][Y]= 5.7735026919e-01; p[0][2][Z]=-5.7735026919e-01;
        p[1][2][X]=-5.7735026919e-01; p[1][2][Y]= 5.7735026919e-01; p[1][2][Z]= 5.7735026919e-01;
        
        p[0][3][X]= 5.7735026919e-01; p[0][3][Y]= 5.7735026919e-01; p[0][3][Z]= 5.7735026919e-01;
        p[1][3][X]= 5.7735026919e-01; p[1][3][Y]=-5.7735026919e-01; p[1][3][Z]= 5.7735026919e-01;
        
        p[0][4][X]= 5.7735026919e-01; p[0][4][Y]= 5.7735026919e-01; p[0][4][Z]=-5.7735026919e-01;
        p[1][4][X]= 5.7735026919e-01; p[1][4][Y]=-5.7735026919e-01; p[1][4][Z]= 5.7735026919e-01;
        
        p[0][5][X]= 5.7735026919e-01; p[0][5][Y]= 5.7735026919e-01; p[0][5][Z]= 5.7735026919e-01;
        p[1][5][X]=-5.7735026919e-01; p[1][5][Y]= 5.7735026919e-01; p[1][5][Z]= 5.7735026919e-01;

        minBurgIndex = 0;
        maxBurgIndex = 5;

/*
 *      Grab the command-line argument values needed for this function
 */
        cubeSize   = (real8)inArgs->cubeLength;
        numLoops   = inArgs->numLoops;
        radius     = inArgs->radius;
        seed       = inArgs->seed;
        dislocType = inArgs->type;

/*
 *      Normalize glide planes.  The burgers vectors that go into the
 *      nodal data file shold NOT be normalized!
 */
        m = 0;
        for (i = 0; i < 6; i++) {
           for (j = 0; j < 2; j++) {
              cross(b[i],p[j][i],e[j][i]);
              NormalizeVec(e[j][i]);
           }
        }

/*
 *      Define the translation vectors as 6 X 4 X 3 array
 *
 *          tr[ ][ ][ ] - translation vector
 *             |  |  |
 *             |  |  x,y,z
 *             |  |
 *             |  4 max  nodes on the loop
 *             |
 *             6 type burgers vectors
 */

        for (itrx = 0; itrx < 6; itrx++) {
            int     itry;
            real8   XX[4][3];

            VECTOR_ZERO(XX[0]);

            m = 1;

            for (j = 0; j < 2; j++) {
                XX[m][0] = XX[m-1][0] + e[j][itrx][0];
                XX[m][1] = XX[m-1][1] + e[j][itrx][1];
                XX[m][2] = XX[m-1][2] + e[j][itrx][2];
                m++;
            }

            for (j = 0; j < 1; j++) {
                XX[m][0] = XX[m-1][0] - e[j][itrx][0];
                XX[m][1] = XX[m-1][1] - e[j][itrx][1];
                XX[m][2] = XX[m-1][2] - e[j][itrx][2];
                m++;
            }


            for (itry = 0; itry < m; itry++) {
                VECTOR_COPY(tr[itrx][itry], XX[itry]);
            }
        }

/*
 *      We use <offSet> to track the position of the node within the current
 *      block of nodes in <inData> that corresponds to the first node in
 *      the current chain.
 */
        offSet = 0;

/*
 *      Create one loop at a time, cycling through burgers vectors as we go.
 */
        inData->nodeCount = 0;
        burgIndex = maxBurgIndex;

        for (loopIndex = 0; loopIndex < numLoops; loopIndex++) 
        {
            if (++burgIndex > maxBurgIndex) 
                burgIndex = minBurgIndex;

/*
 *          Increase the size of the node array enough to hold
 *          all the new nodes created for this loop.
 */
            numSegs = 4;

            newNodeIndex = inData->nodeCount;
            inData->nodeCount += numSegs;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numSegs);

            loopCtrX = (randm(&seed)-0.5) * cubeSize;
            loopCtrY = (randm(&seed)-0.5) * cubeSize;
            loopCtrZ = (randm(&seed)-0.5) * cubeSize;

            
            for (id = offSet; id < (offSet + numSegs); id++) {
                int   dIndex;

                node =  &inData->node[id];

/*
 *              Pick the index for the direction component of the tr array
 */
                dIndex = (id - offSet) % numSegs;

                node->x = loopCtrX + radius * tr[burgIndex][dIndex][0];
                node->y = loopCtrY + radius * tr[burgIndex][dIndex][1];
                node->z = loopCtrZ + radius * tr[burgIndex][dIndex][2];

/*
 *              Set up the node
 */
                SET_CONSTRAINTS(node->constraint, UNCONSTRAINED);
                node->myTag.domainID = dislocType;
                node->myTag.index = id + param->nodeCount;

                AllocNodeArms(node, 2);

                node->burgX[0] =  b[burgIndex][X];
                node->burgY[0] =  b[burgIndex][Y];
                node->burgZ[0] =  b[burgIndex][Z];

                node->burgX[1] = -b[burgIndex][X];
                node->burgY[1] = -b[burgIndex][Y];
                node->burgZ[1] = -b[burgIndex][Z];

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index =
                        offSet + (((id-offSet)+1)%numSegs) + param->nodeCount;

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index =
                        offSet + param->nodeCount +
                        (((id-offSet)+numSegs-1)%numSegs);
            }

/*
 *          Need to set the glide plane normal for each segment.
 *          Couldn't do it above because we didn't have all the
 *          neighbor node positions until all the nodes in the
 *          loop were created.
 */
            for (id = offSet; id < (offSet + numSegs); id++) {
                int     nbr1Index;
                real8   vec1[3], vec2[3], tmpVec[3];

                node = &inData->node[id];

                nbr1Index = node->nbrTag[1].index - param->nodeCount;
                nbr1Node = &inData->node[nbr1Index];

                ux = nbr1Node->x - node->x;
                uy = nbr1Node->y - node->y;
                uz = nbr1Node->z - node->z;

                Normalize(&ux,&uy,&uz);

/*
 *              l cross b gives normal vector for each segment, but
 *              in order to avoid some round-off issues, explicitly set
 *              any components of the glide plane normal are near-zero
 *              to be zero, then renormalize.
 */
                vec1[0] = ux;
                vec1[1] = uy;
                vec1[2] = uz;

                vec2[0] =  node->burgX[1];
                vec2[1] =  node->burgY[1];
                vec2[2] =  node->burgZ[1];

                NormalizedCrossVector(vec1, vec2, tmpVec);

                if (fabs(tmpVec[0]) < 1.0e-03) tmpVec[0] = 0.0;
                if (fabs(tmpVec[1]) < 1.0e-03) tmpVec[1] = 0.0;
                if (fabs(tmpVec[2]) < 1.0e-03) tmpVec[2] = 0.0;

                NormalizeVec(tmpVec);

                node->nx[1] = tmpVec[0];
                node->ny[1] = tmpVec[1];
                node->nz[1] = tmpVec[2];

                nbr1Node->nx[0] = tmpVec[0];
                nbr1Node->ny[0] = tmpVec[1];
                nbr1Node->nz[0] = tmpVec[2];

            }  /*  for (id = nextNode; ...) */

            offSet += numSegs;

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          Write blocks of nodal data to the data files as needed.
 */
            lastBlock = (loopIndex == (numLoops-1));

            if (inArgs->useSegmentedDataFile || lastBlock) {
                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);
                startRemeshIndex = 0;
                offSet = 0;
            }

        }  /* for (loopIndex = 0; ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFiniteMixedConfig
 *      Description:  Creates an initial configuration consisting of
 *                    random screw and edge dislocations, but assumes a
 *                    finite problem space in at least 1 dimension and
 *                    terminates the dislocations at the free surface(s) 
 *
 *      Arguments:
 *          IN:  inArgs       Structure containing values associated with the
 *                            command line arguments.
 *          OUT: totDislocLen Total length of all dislocations inserted into
 *                            the configuration.
 *
 *      Values used from the <inArgs> structure:
 *
 *          numChains     Number of chains to create
 *          seed          Seed value for random number generator
 *          type          Type of dislocation configuration being generated.
 *
 *-------------------------------------------------------------------------*/
void CreateFiniteMixedConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                             real8 *totDislocLen)
{
        int      i, lineIndex, chain, numPoints;
        int      newNodeIndex, lastBlock, burgIndex;
        int      isScrew, gpIndex, gpBurgIndex;
        int      startRemeshIndex = 0;
        int      signFact; 
        int      nbr1Index, nbr2Index, numConnections, numSegs;
        int      intersectIndex1=0, intersectIndex2=0;
        int      numChains, seed, dislocType;
        int      pbc[3];
        real8    burgSign, vecLen;
        real8    surfCoord, len, minLen1, minLen2;
        real8    intersectSurfCoord1=0, intersectSurfCoord2=0;
        real8    p0[3], newPos[3];
        real8    range[3], rangeCntr[3];
        real8    minCoord[3], maxCoord[3];
        real8    intersectPos1[3], intersectPos2[3], vector[3];
        real8    normalizedBurg[3];
        real8    linedir[16][3], glidePlane[4][6][3], plane[3];
        Node_t   *node;
        Param_t  *param;
        

        param = inData->param;

        pbc[X] = (param->xBoundType == Periodic);
        pbc[Y] = (param->yBoundType == Periodic);
        pbc[Z] = (param->zBoundType == Periodic);

/*
 *      Grap the command-line argument values needed for this function
 */
        numChains   = inArgs->numChains;
        seed        = inArgs->seed;
        dislocType  = inArgs->type;

/*
 *      Just a couple sanity check...
 */
        if (pbc[X]*pbc[Y]*pbc[Z] != 0) {
            Fatal("CreateFiniteMixedConfig: Requires free surfaces in at"
                  "least one dimension\n    but all boundaries are periodic!");
        }

        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be > 0.\n",
                  "CreateFiniteMixedConfig", numChains);
        }
        
/*
 *      Set up an array of dislocation line directions that
 *      are used to create the new dislocations.  There are
 *      essentially 4 sets of line directions (1 per burgers
 *      vector) with 4 line directions per set; the first for
 *      the screw and the following three for edge.
 *      The burger's vector for all 4 dislocation lines in
 *      a set is the same as the line direction of the screw
 *      dislocation in the group.
 */

/*
 *      Type  [1 1 1] burgers vector
 */
        linedir[0][0] =  0.5773503;  
        linedir[0][1] =  0.5773503;
        linedir[0][2] =  0.5773503; 
        
        linedir[1][0] = -0.8164966;
        linedir[1][1] =  0.4082483;
        linedir[1][2] =  0.4082483; 
        
        linedir[2][0] =  0.4082483;
        linedir[2][1] = -0.8164966;
        linedir[2][2] =  0.4082483; 
        
        linedir[3][0] =  0.4082483;
        linedir[3][1] =  0.4082483; 
        linedir[3][2] = -0.8164966;
        
/*
 *      Type [-1 1 1] burgers vector
 */
        linedir[4][0] = -0.5773503;
        linedir[4][1] =  0.5773503;
        linedir[4][2] =  0.5773503; 
        
        linedir[5][0] =  0.8164966;
        linedir[5][1] =  0.4082483;
        linedir[5][2] =  0.4082483; 
        
        linedir[6][0] =  0.4082483;
        linedir[6][1] =  0.8164966;
        linedir[6][2] = -0.4082483; 
        
        linedir[7][0] =  0.4082483;
        linedir[7][1] = -0.4082483; 
        linedir[7][2] =  0.8164966;
        
/*
 *      Type [1 -1 1] burgers vector
 */
        linedir[8][0] =  0.5773503;
        linedir[8][1] = -0.5773503;
        linedir[8][2] =  0.5773503; 
        
        linedir[9][0] =  0.4082483;
        linedir[9][1] =  0.8164966;
        linedir[9][2] =  0.4082483; 
        
        linedir[10][0] =  0.8164966;
        linedir[10][1] =  0.4082483;
        linedir[10][2] = -0.4082483; 
        
        linedir[11][0] = -0.4082483;
        linedir[11][1] =  0.4082483; 
        linedir[11][2] =  0.8164966;
        
/*
 *      Type [1 1 -1] burgers vector
 */
        linedir[12][0] =  0.5773503;
        linedir[12][1] =  0.5773503;
        linedir[12][2] = -0.5773503; 
        
        linedir[13][0] = -0.4082483;
        linedir[13][1] =  0.8164966;
        linedir[13][2] =  0.4082483; 
        
        linedir[14][0] =  0.8164966;
        linedir[14][1] = -0.4082483;
        linedir[14][2] =  0.4082483; 
        
        linedir[15][0] =  0.4082483;
        linedir[15][1] =  0.4082483; 
        linedir[15][2] =  0.8164966;
        
/*
 *      Set up the valid glide planes for each screw burgers vector,
 *      six glide planes per burgers vector.  For edges, glide plane
 *      will simply be cross product between burgers vector and 
 *      linedir.
 *
 *      glide planes for [1 1 1]
 */
        glidePlane[0][0][0] =  0.7071068;
        glidePlane[0][0][1] = -0.7071068;
        glidePlane[0][0][2] =  0.0000000;

        glidePlane[0][1][0] =  0.7071068;
        glidePlane[0][1][1] =  0.0000000;
        glidePlane[0][1][2] = -0.7071068;

        glidePlane[0][2][0] =  0.0000000;
        glidePlane[0][2][1] =  0.7071068;
        glidePlane[0][2][2] = -0.7071068;

        glidePlane[0][3][0] = -0.7071068;
        glidePlane[0][3][1] =  0.7071068;
        glidePlane[0][3][2] = -0.0000000;

        glidePlane[0][4][0] = -0.7071068;
        glidePlane[0][4][1] = -0.0000000;
        glidePlane[0][4][2] =  0.7071068;

        glidePlane[0][5][0] = -0.0000000;
        glidePlane[0][5][1] = -0.7071068;
        glidePlane[0][5][2] =  0.7071068;

/*
 *      glide planes for [-1 1 1]
 */
        glidePlane[1][0][0] =  0.0000000;
        glidePlane[1][0][1] =  0.7071068;
        glidePlane[1][0][2] = -0.7071068;

        glidePlane[1][1][0] =  0.7071068;
        glidePlane[1][1][1] =  0.0000000;
        glidePlane[1][1][2] =  0.7071068;

        glidePlane[1][2][0] =  0.0000000;
        glidePlane[1][2][1] =  0.7071068;
        glidePlane[1][2][2] = -0.7071068;

        glidePlane[1][3][0] = -0.0000000;
        glidePlane[1][3][1] = -0.7071068;
        glidePlane[1][3][2] =  0.7071068;

        glidePlane[1][4][0] = -0.7071068;
        glidePlane[1][4][1] = -0.0000000;
        glidePlane[1][4][2] = -0.7071068;

        glidePlane[1][5][0] = -0.0000000;
        glidePlane[1][5][1] = -0.7071068;
        glidePlane[1][5][2] =  0.7071068;

/*
 *      glide planes for [1 -1 1]
 */
        glidePlane[2][0][0] =  0.7071068;
        glidePlane[2][0][1] =  0.7071068;
        glidePlane[2][0][2] =  0.0000000;

        glidePlane[2][1][0] =  0.7071068;
        glidePlane[2][1][1] =  0.0000000;
        glidePlane[2][1][2] = -0.7071068;

        glidePlane[2][2][0] =  0.0000000;
        glidePlane[2][2][1] =  0.7071068;
        glidePlane[2][2][2] =  0.7071068;

        glidePlane[2][3][0] = -0.7071068;
        glidePlane[2][3][1] = -0.7071068;
        glidePlane[2][3][2] = -0.0000000;

        glidePlane[2][4][0] = -0.7071068;
        glidePlane[2][4][1] = -0.0000000;
        glidePlane[2][4][2] =  0.7071068;

        glidePlane[2][5][0] = -0.0000000;
        glidePlane[2][5][1] = -0.7071068;
        glidePlane[2][5][2] = -0.7071068;

/*
 *      glide planes for [1 1 -1]
 */
        glidePlane[3][0][0] =  0.7071068;
        glidePlane[3][0][1] = -0.7071068;
        glidePlane[3][0][2] =  0.0000000;

        glidePlane[3][1][0] =  0.0000000;
        glidePlane[3][1][1] =  0.7071068;
        glidePlane[3][1][2] =  0.7071068;

        glidePlane[3][2][0] =  0.7071068;
        glidePlane[3][2][1] =  0.0000000;
        glidePlane[3][2][2] =  0.7071068;

        glidePlane[3][3][0] = -0.7071068;
        glidePlane[3][3][1] =  0.7071068;
        glidePlane[3][3][2] = -0.0000000;

        glidePlane[3][4][0] = -0.0000000;
        glidePlane[3][4][1] = -0.7071068;
        glidePlane[3][4][2] = -0.7071068;

        glidePlane[3][5][0] = -0.7071068;
        glidePlane[3][5][1] = -0.0000000;
        glidePlane[3][5][2] = -0.7071068;

        minCoord[X] = param->xBoundMin;
        maxCoord[X] = param->xBoundMax;
        minCoord[Y] = param->yBoundMin;
        maxCoord[Y] = param->yBoundMax;
        minCoord[Z] = param->zBoundMin;
        maxCoord[Z] = param->zBoundMax;

        for (i = 0; i < 3; i++) {
            range[i] = maxCoord[i] - minCoord[i];
            rangeCntr[i] = 0.5 * (minCoord[i] + maxCoord[i]);
        }

/*
 *      Create the specified number of chains.  Anytime the number of
 *      nodes maintained in memory exceeds the threshhold, write the
 *      block of nodal data out to the data file.
 */
        inData->nodeCount = 0;

        for (chain = 0; chain < numChains; chain++) {
            int baseNodeID = inData->nodeCount;
        
            numPoints = 3;
            lineIndex = chain % 16; 
            burgIndex = 4 * (lineIndex / 4);
            gpBurgIndex = lineIndex / 4;
            gpIndex = (chain / 16) % 6;
            isScrew = (chain % 4) == 0;

            normalizedBurg[X] = linedir[burgIndex][X];
            normalizedBurg[Y] = linedir[burgIndex][Y];
            normalizedBurg[Z] = linedir[burgIndex][Z];

            Normalize(&normalizedBurg[X], &normalizedBurg[Y],
                      &normalizedBurg[Z]);

/*
 *          Select an initial point (p0) that is within the boundaries
 *          of the simulation and then calculate the positions at which
 *          a dislocation line with the given line direction would
 *          interesect the nearest free surface in each direction.
 */
            for (i = 0; i < 3; i++) {
                p0[i] = (randm(&seed)-0.5) * range[i] + rangeCntr[i];
            }

            minLen1 = 1.0e+20;
            minLen2 = 1.0e+20;

            for (i = 0; i < 3; i++) {
                if (pbc[i] == 0) {
                    signFact = (linedir[burgIndex][i] < 0.0 ? -1 : 1);
                    surfCoord = signFact > 0 ? maxCoord[i] : minCoord[i];
                    len = fabs((surfCoord - p0[i]) / normalizedBurg[i]);
                    if (len < minLen1) {
                        minLen1 = len;
                        intersectIndex1 = i;
                        intersectSurfCoord1 = surfCoord;
                    }

                    signFact = -signFact;
                    surfCoord = signFact > 0 ? maxCoord[i] : minCoord[i];
                    len = fabs((surfCoord - p0[i]) / normalizedBurg[i]);
                    if (len < minLen2) {
                        minLen2 = len;
                        intersectIndex2 = i;
                        intersectSurfCoord2 = surfCoord;
                    }
                }
            }

/*
 *          We know how far the dislocation can extend in each direction, now
 *          calculate the exact intersect point in both directions
 */
            for (i = 0; i < 3; i++) {

                if (i == intersectIndex1) {
                    intersectPos1[i] = intersectSurfCoord1;
                } else {
                    intersectPos1[i] = p0[i] + (minLen1 * normalizedBurg[i]);
                }

                if (i == intersectIndex2) {
                    intersectPos2[i] = intersectSurfCoord2;
                } else {
                    intersectPos2[i] = p0[i] - (minLen2 * normalizedBurg[i]);
                }
            }

/*
 *          Find a vector from the first intersection point to the second,
 *          calculate how many segments the line should be broken into based
 *          on the <maxSeg> value.
 */
            for (i = 0; i < 3; i++) {
                vector[i] = intersectPos2[i] - intersectPos1[i];
            }

            vecLen = sqrt(vector[0]*vector[0] +
                          vector[1]*vector[1] +
                          vector[2]*vector[2]);

            numSegs = (int)(vecLen / (.95 * param->maxSeg)) + 1;
            numPoints = numSegs + 1;

            for (i = 0; i < 3; i++) {
                vector[i] /= (real8)numSegs;
            }

/*
 *          Reallocate the node array with sufficient size to add
 *          all the new nodes defining this chain.
 */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += numPoints;
            inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
                                             * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);
        

/*
 *          Starting with the first intersection point, create a
 *          series of dislocation segments ending at the second point.
 */
            newPos[X] = intersectPos1[X];
            newPos[Y] = intersectPos1[Y];
            newPos[Z] = intersectPos1[Z];

            FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

            for (i = 0; i < numPoints; i++) {
                if (i == 0) {
                    numConnections = 1;
                    burgSign = -1.0;
                    nbr1Index = 1;
                } else if (i == (numPoints - 1)) {
                    numConnections = 1;
                    burgSign = 1.0;
                    nbr1Index = i - 1;
/*
 *                  Make sure the final node is at the surface
 *                  intersection point
 */
                    newPos[X] = intersectPos2[X];
                    newPos[Y] = intersectPos2[Y];
                    newPos[Z] = intersectPos2[Z];
                    FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);
                } else {
                    numConnections = 2;
                    burgSign = 1.0;
                    nbr1Index = i - 1;
                    nbr2Index = i + 1;
                }

                node = &inData->node[baseNodeID+i];
                node->myTag.domainID = dislocType;
                node->myTag.index    = baseNodeID + i + param->nodeCount;

                node->x = newPos[X];
                node->y = newPos[Y];
                node->z = newPos[Z];

                SET_CONSTRAINTS(node->constraint, UNCONSTRAINED);

                AllocNodeArms(node, numConnections);

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index    = baseNodeID + nbr1Index +
                                           param->nodeCount;

                node->burgX[0] = burgSign * linedir[burgIndex][X];
                node->burgY[0] = burgSign * linedir[burgIndex][Y];
                node->burgZ[0] = burgSign * linedir[burgIndex][Z];

/*
 *              For screw dislocations, use the glide plane from the table.
 *              For edge dislocations, glide plane is the cross product
 *              of the burgers vector and line direction.
 */
                if (isScrew) {
                    node->nx[0] = glidePlane[gpBurgIndex][gpIndex][X];
                    node->ny[0] = glidePlane[gpBurgIndex][gpIndex][Y];
                    node->nz[0] = glidePlane[gpBurgIndex][gpIndex][Z];
                    if (numConnections == 2) {
                        node->nx[1] = glidePlane[gpBurgIndex][gpIndex][X];
                        node->ny[1] = glidePlane[gpBurgIndex][gpIndex][Y];
                        node->nz[1] = glidePlane[gpBurgIndex][gpIndex][Z];
                    }
                } else {
                    cross(linedir[burgIndex], linedir[lineIndex], plane);
                    plane[X] = (floor(plane[X] * 1.0e+07)) * 1.0e-07;
                    plane[Y] = (floor(plane[Y] * 1.0e+07)) * 1.0e-07;
                    plane[Z] = (floor(plane[Z] * 1.0e+07)) * 1.0e-07;
                    node->nx[0] = plane[X];
                    node->ny[0] = plane[Y];
                    node->nz[0] = plane[Z];
                    if (numConnections == 2) {
                        node->nx[1] = plane[X];
                        node->ny[1] = plane[Y];
                        node->nz[1] = plane[Z];
                    }
                }

/*
 *              Calculate the next node's position relative to this one.
 */
                newPos[X] += vector[X];
                newPos[Y] += vector[Y];
                newPos[Z] += vector[Z];

                FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

                if (numConnections == 1) {
                    ADD_CONSTRAINTS(node->constraint, SURFACE_NODE);
                    continue;
                }

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index    = baseNodeID + nbr2Index + 
                                           param->nodeCount;

                node->burgX[1] = -burgSign * linedir[burgIndex][X];
                node->burgY[1] = -burgSign * linedir[burgIndex][Y];
                node->burgZ[1] = -burgSign * linedir[burgIndex][Z];
            }

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          Write blocks of nodal data to the data files as needed.
 */
            lastBlock = (chain == (numChains - 1));

            if (inArgs->useSegmentedDataFile || lastBlock) {
                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);
                startRemeshIndex = 0;
            }

        }  /* for (chain = 0; chain < numChains; ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateFRSource
 *      Description:  Generate a configuration consisting of one or more
 *                    frank-read sources of random lengths between the
 *                    specified minimum and maximum lengths.  
 *
 *                    NOTE:  This function currently assumes periodic
 *                    boundary conditions!
 *
 *      Arguments:
 *          IN:  inArgs       Structure containing values associated with the
 *                            command line arguments.
 *          OUT: totDislocLen Total length of all dislocations inserted into
 *                            the configuration.
 *
 *      Values used from the <inArgs> structure:
 *          cubeLength    Length of cubic problem space in a single
 *                        dimension (units of b)
 *          numFRSrcs     Number of frank read sources to create
 *          frLenMin      Minimum length of any created frank-read source
 *          frLenMax      Maximum length of any created frank-read source
 *          seed          Seed value for random number generator
 *          type          Type of dislocation configuration being generated.
 *
 *-------------------------------------------------------------------------*/
void CreateFRSource(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                    real8 *totDislocLen)
{
        int      i, lineIndex, chain, numPoints;
        int      newNodeIndex, lastBlock, burgIndex;
        int      isScrew = 1, gpIndex, gpBurgIndex;
        int      lenRange;
        int      startRemeshIndex = 0;
        int      numSources, srcLenMin, srcLenMax, cubeLength;
        int      seed, dislocType;
        real8    srcLen;
        real8    p0[3], p1[3], p2[3];
        real8    linedir[4][3], glidePlane[4][6][3], plane[3];
        Node_t   *node;
        Param_t  *param;


        param = inData->param;

/*
 *      Grap the command-line argument values needed for this function
 */
        cubeLength  = inArgs->cubeLength;
        numSources  = inArgs->numFRSrcs;
        srcLenMin   = inArgs->frLenMin;
        srcLenMax   = inArgs->frLenMax;
        seed        = inArgs->seed;
        dislocType  = inArgs->type;

        lenRange = srcLenMax - srcLenMin;

/*
 *      Just a quick sanity check...
 */
        if (numSources <= 0) {
            Fatal("%s: numSources is %d, but must be > 0.\n",
                  "CreateFRSource", numSources);
        }


/*
 *      Set up an array of dislocation line directions that
 *      are used to create the new dislocations.  There are
 *      essentially 4 sets of line directions (1 per burgers
 *      vector) with 4 line directions per set; the first for
 *      the screw and the following three for edge.
 *      The burger's vector for all 4 dislocation lines in
 *      a set is the same as the line direction of the screw
 *      dislocation in the group.
 */

/*
 *      Type  [1 1 1] burgers vector
 */
        linedir[0][0] =  0.5;
        linedir[0][1] =  0.5;
        linedir[0][2] =  0.5;

/*
 *      Type [-1 1 1] burgers vector
 */
        linedir[1][0] = -0.5;
        linedir[1][1] =  0.5;
        linedir[1][2] =  0.5;

/*
 *      Type [1 -1 1] burgers vector
 */
        linedir[2][0] =  0.5;
        linedir[2][1] = -0.5;
        linedir[2][2] =  0.5;

/*
 *      Type [1 1 -1] burgers vector
 */
        linedir[3][0] =  0.5;
        linedir[3][1] =  0.5;
        linedir[3][2] = -0.5;

/*
 *      Set up the valid glide planes for each screw burgers vector,
 *      six glide planes per burgers vector.  For edges, glide plane
 *      will simply be cross product between burgers vector and
 *      linedir.
 *
 *      glide planes for [1 1 1]
 */
        glidePlane[0][0][0] =  0.7071068;
        glidePlane[0][0][1] = -0.7071068;
        glidePlane[0][0][2] =  0.0000000;

        glidePlane[0][1][0] =  0.7071068;
        glidePlane[0][1][1] =  0.0000000;
        glidePlane[0][1][2] = -0.7071068;

        glidePlane[0][2][0] =  0.0000000;
        glidePlane[0][2][1] =  0.7071068;
        glidePlane[0][2][2] = -0.7071068;

        glidePlane[0][3][0] = -0.7071068;
        glidePlane[0][3][1] =  0.7071068;
        glidePlane[0][3][2] = -0.0000000;

        glidePlane[0][4][0] = -0.7071068;
        glidePlane[0][4][1] = -0.0000000;
        glidePlane[0][4][2] =  0.7071068;

        glidePlane[0][5][0] = -0.0000000;
        glidePlane[0][5][1] = -0.7071068;
        glidePlane[0][5][2] =  0.7071068;

/*
 *      glide planes for [-1 1 1]
 */
        glidePlane[1][0][0] =  0.0000000;
        glidePlane[1][0][1] =  0.7071068;
        glidePlane[1][0][2] = -0.7071068;

        glidePlane[1][1][0] =  0.7071068;
        glidePlane[1][1][1] =  0.0000000;
        glidePlane[1][1][2] =  0.7071068;

        glidePlane[1][2][0] =  0.0000000;
        glidePlane[1][2][1] =  0.7071068;
        glidePlane[1][2][2] = -0.7071068;

        glidePlane[1][3][0] = -0.0000000;
        glidePlane[1][3][1] = -0.7071068;
        glidePlane[1][3][2] =  0.7071068;

        glidePlane[1][4][0] = -0.7071068;
        glidePlane[1][4][1] = -0.0000000;
        glidePlane[1][4][2] = -0.7071068;

        glidePlane[1][5][0] = -0.0000000;
        glidePlane[1][5][1] = -0.7071068;
        glidePlane[1][5][2] =  0.7071068;

/*
 *      glide planes for [1 -1 1]
 */
        glidePlane[2][0][0] =  0.7071068;
        glidePlane[2][0][1] =  0.7071068;
        glidePlane[2][0][2] =  0.0000000;

        glidePlane[2][1][0] =  0.7071068;
        glidePlane[2][1][1] =  0.0000000;
        glidePlane[2][1][2] = -0.7071068;

        glidePlane[2][2][0] =  0.0000000;
        glidePlane[2][2][1] =  0.7071068;
        glidePlane[2][2][2] =  0.7071068;

        glidePlane[2][3][0] = -0.7071068;
        glidePlane[2][3][1] = -0.7071068;
        glidePlane[2][3][2] = -0.0000000;

        glidePlane[2][4][0] = -0.7071068;
        glidePlane[2][4][1] = -0.0000000;
        glidePlane[2][4][2] =  0.7071068;

        glidePlane[2][5][0] = -0.0000000;
        glidePlane[2][5][1] = -0.7071068;
        glidePlane[2][5][2] = -0.7071068;

/*
 *      glide planes for [1 1 -1]
 */
        glidePlane[3][0][0] =  0.7071068;
        glidePlane[3][0][1] = -0.7071068;
        glidePlane[3][0][2] =  0.0000000;

        glidePlane[3][1][0] =  0.0000000;
        glidePlane[3][1][1] =  0.7071068;
        glidePlane[3][1][2] =  0.7071068;

        glidePlane[3][2][0] =  0.7071068;
        glidePlane[3][2][1] =  0.0000000;
        glidePlane[3][2][2] =  0.7071068;

        glidePlane[3][3][0] = -0.7071068;
        glidePlane[3][3][1] =  0.7071068;
        glidePlane[3][3][2] = -0.0000000;

        glidePlane[3][4][0] = -0.0000000;
        glidePlane[3][4][1] = -0.7071068;
        glidePlane[3][4][2] = -0.7071068;

        glidePlane[3][5][0] = -0.7071068;
        glidePlane[3][5][1] = -0.0000000;
        glidePlane[3][5][2] = -0.7071068;

/*
 *      Create the specified number of chains.  Anytime the number of
 *      nodes maintained in memory exceeds the threshhold, write the
 *      block of nodal data out to the data file.
 */
        inData->nodeCount = 0;

        for (chain = 0; chain < numSources; chain++) {
            int baseNodeID = inData->nodeCount;

            numPoints = 3;
            lineIndex = chain % 4;
            burgIndex = lineIndex;
            gpBurgIndex = lineIndex;
            gpIndex = (chain / 4) % 6;

/*
 *          Reallocate the node array with sufficient size to add
 *          all the new nodes defining this chain.
 */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += numPoints;
            inData->node = (Node_t *)realloc(inData->node, inData->nodeCount
                                               * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numPoints);

/*
 *          Length of the frank read source should be a random length
 *          between srcLenMin and srcLenMax.
 */
            if (lenRange > 0) {
                srcLen = (srcLenMin + (randm(&seed) * lenRange)) / sqrt(3.0);
            } else {
                srcLen = srcLenMin / sqrt(3.0);
            }

            for (i = 0; i < 3; i++) {
                p0[i] = (randm(&seed)-0.5) * cubeLength * 2.0;
                p1[i] = p0[i] - (srcLen * linedir[lineIndex][i]);
                p2[i] = p0[i] + (srcLen * linedir[lineIndex][i]);
            }

/*
 *          Set up the 3 nodes we're using to define the dislocation line.
 *          First do point p0.
 */
            node = &inData->node[baseNodeID];

            node->myTag.domainID = dislocType;
            node->myTag.index    = baseNodeID + param->nodeCount;
            node->x = p0[0];
            node->y = p0[1];
            node->z = p0[2];

            SET_CONSTRAINTS(node->constraint, UNCONSTRAINED);

            AllocNodeArms(node, 2);

            node->nbrTag[0].domainID = dislocType;
            node->nbrTag[0].index    = baseNodeID + 1 + param->nodeCount;
            node->nbrTag[1].domainID = dislocType;
            node->nbrTag[1].index    = baseNodeID + 2 + param->nodeCount;

            node->burgX[0] = linedir[burgIndex][0];
            node->burgY[0] = linedir[burgIndex][1];
            node->burgZ[0] = linedir[burgIndex][2];

            node->burgX[1] = -linedir[burgIndex][0];
            node->burgY[1] = -linedir[burgIndex][1];
            node->burgZ[1] = -linedir[burgIndex][2];

/*
 *          For screw dislocations, use the glide plane from the table.
 *          For edge dislocations, glide plane is the cross product
 *          of the burgers vector and line direction.
 */
            if (isScrew) {
                node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
                node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
                node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];
                node->nx[1] = glidePlane[gpBurgIndex][gpIndex][0];
                node->ny[1] = glidePlane[gpBurgIndex][gpIndex][1];
                node->nz[1] = glidePlane[gpBurgIndex][gpIndex][2];
            } else {
                cross(linedir[burgIndex], linedir[lineIndex], plane);
                plane[0] = (floor(plane[0] * 1.0e+07)) * 1.0e-07;
                plane[1] = (floor(plane[1] * 1.0e+07)) * 1.0e-07;
                plane[2] = (floor(plane[2] * 1.0e+07)) * 1.0e-07;
                node->nx[0] = plane[0];
                node->ny[0] = plane[1];
                node->nz[0] = plane[2];
                node->nx[1] = plane[0];
                node->ny[1] = plane[1];
                node->nz[1] = plane[2];
            }

/*
 *          Now point p1...
 */
            node = &inData->node[baseNodeID+1];

            node->myTag.domainID = dislocType;
            node->myTag.index    = baseNodeID + 1 + param->nodeCount;

            node->x = p1[0];
            node->y = p1[1];
            node->z = p1[2];

            SET_CONSTRAINTS(node->constraint, PINNED_NODE);

            AllocNodeArms(node, 1);

            node->nbrTag[0].domainID = dislocType;
            node->nbrTag[0].index    = baseNodeID + param->nodeCount;

            node->burgX[0] = -linedir[burgIndex][0];
            node->burgY[0] = -linedir[burgIndex][1];
            node->burgZ[0] = -linedir[burgIndex][2];
            node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
            node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
            node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];

/*
 *          Now point p2...
 */
            node = &inData->node[baseNodeID+2];

            node->myTag.domainID = dislocType;
            node->myTag.index    = baseNodeID + 2 + param->nodeCount;

            node->x = p2[0];
            node->y = p2[1];
            node->z = p2[2];

            SET_CONSTRAINTS(node->constraint, PINNED_NODE);

            AllocNodeArms(node, 1);

            node->nbrTag[0].domainID = dislocType;
            node->nbrTag[0].index    = baseNodeID + param->nodeCount;

            node->burgX[0] = linedir[burgIndex][0];
            node->burgY[0] = linedir[burgIndex][1];
            node->burgZ[0] = linedir[burgIndex][2];
            node->nx[0] = glidePlane[gpBurgIndex][gpIndex][0];
            node->ny[0] = glidePlane[gpBurgIndex][gpIndex][1];
            node->nz[0] = glidePlane[gpBurgIndex][gpIndex][2];

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          Write blocks of nodal data to the data files as needed.
 */
            lastBlock = (chain == (numSources - 1));

            if (inArgs->useSegmentedDataFile || lastBlock) {
                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);
                startRemeshIndex = 0;
            }

        }  /* for (chain = 0; chain < numChains; ...) */

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     CreateEdges
 *      Description:  Creates "edge" dislocations in the [100] [010] and [001]
 *                    directions.  Each of the three line senses has 4
 *                    combinations of burgers vector and normals.  With the
 *                    opposite sign of the line sense vectors considered,
 *                    the total number of types of dislocations is 24,
 *                    so number of chains specified should be a multiple
 *                    of 24 to include all types of.  (Okay, they're
 *                    not pure edge, but they are not screw either)
 *
 *      Arguments:
 *          IN:  inArgs       Structure containing values associated with the
 *                            command line arguments.
 *          OUT: totDislocLen Total length of all dislocations inserted into
 *                            the configuration.
 *
 *      Values used from the <inArgs> structure:
 *          cubeLength    Length of cubic problem space in a single
 *                        dimension (units of b)
 *          numChains     Number of chains to create
 *          seed          Seed value for random number generator
 *          type          Type of dislocation configuration being generated.
 *
 *-------------------------------------------------------------------------*/
void CreateEdges(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                 real8 *totDislocLen)
{
        int           ic, ip, np, offSet, lastBlock, signFactor;
        int           newNodeIndex, startRemeshIndex;
        int           ldIndex, gpIndex, burgIndex, nbr1Index, nbr2Index;
        int           numChains, seed, dislocType;
        real8         posFactor, cubeSize;
        real8         xp[3], yp[3], zp[3];
        Param_t       *param;
        Node_t        *node;
        static real8  lineDir[3][3] = {
                      {0.0, 0.0, 1.0},
                      {0.0, 1.0, 0.0},
                      {1.0, 0.0, 0.0}};
        static real8  burg[4][3] = {
                      { 0.5773503,  0.5773503,  0.5773503},
                      { 0.5773503,  0.5773503, -0.5773503},
                      { 0.5773503, -0.5773503,  0.5773503},
                      {-0.5773503,  0.5773503,  0.5773503}};
        static real8  glidePlane[12][3] = {
                      {  0.7071068, -0.7071068,  0},  /* ldir [001] b [111]  */
                      {  0.7071068, -0.7071068,  0},  /* ldir [001] b [11-1] */
                      {  0.7071068,  0.7071068,  0},  /* ldir [001] b [1-11] */
                      {  0.7071068,  0.7071068,  0},  /* ldir [001] b [-111] */
                      {  0.7071068,  0, -0.7071068},  /* ldir [010] b [111]  */
                      {  0.7071068,  0,  0.7071068},  /* ldir [010] b [11-1] */
                      {  0.7071068,  0, -0.7071068},  /* ldir [010] b [1-11] */
                      {  0.7071068,  0,  0.7071068},  /* ldir [010] b [-111] */
                      {  0, -0.7071068,  0.7071068},  /* ldir [100] b [111]  */
                      {  0,  0.7071068,  0.7071068},  /* ldir [100] b [11-1] */
                      {  0,  0.7071068,  0.7071068},  /* ldir [100] b [1-11] */
                      {  0, -0.7071068,  0.7071068}}; /* ldir [100] b [-111] */

        param = inData->param;

/*
 *      Grap the command-line argument values needed for this function
 */
        cubeSize    = (real8)inArgs->cubeLength;
        numChains   = inArgs->numChains;
        seed        = inArgs->seed;
        dislocType  = inArgs->type;

        posFactor   = 0.333 * cubeSize;

        if (numChains <= 0) {
            Fatal("%s: numChains is %d, but must be > 0.\n",
                  "CreateEdges", numChains);
        }

/*
 *      We use <offSet> to track the position of the node within the current
 *      block of nodes in <inData> that corresponds to the first node in
 *      the current chain.
 */
        offSet = 0;

        inData->nodeCount = 0;
        startRemeshIndex = 0;

/*
 *      Create the specified number of chains.
 */
        for (ic = 0; ic < numChains; ic++) {

            gpIndex = ic % 12; 
            burgIndex = gpIndex % 4;
            ldIndex = gpIndex / 4;
            
/*
 *          First 12 burgers vector/normal sets use positive line
 *          sense, next set of 12 uses opposite line sense, and
 *          so on.
 */
            signFactor = ((ic / 12) & 0x01) ? -1 : 1;

            np = 3;
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += np;

            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * np);

/*
 *          Set up 3 initial points for the line.  Point 1 is a base position
 *          at a random location, point 0 is in the negative direction along
 *          the line and point 2 is in the positive direction along the line.
 */
            xp[1] = (randm(&seed)-0.5)*cubeSize;
            yp[1] = (randm(&seed)-0.5)*cubeSize;
            zp[1] = (randm(&seed)-0.5)*cubeSize;

            xp[0] = xp[1] - (posFactor * signFactor * lineDir[ldIndex][X]);
            yp[0] = yp[1] - (posFactor * signFactor * lineDir[ldIndex][Y]);
            zp[0] = zp[1] - (posFactor * signFactor * lineDir[ldIndex][Z]);

            xp[2] = xp[1] + (posFactor * signFactor * lineDir[ldIndex][X]);
            yp[2] = yp[1] + (posFactor * signFactor * lineDir[ldIndex][Y]);
            zp[2] = zp[1] + (posFactor * signFactor * lineDir[ldIndex][Z]);
/*
 *          Loop over the points and set up the nodes, link them to 
 *          the neighbor nodes, etc.
 */
            for (ip = 0; ip < np; ip++) {

                node = &inData->node[ip+offSet];

                node->x = xp[ip];
                node->y = yp[ip];
                node->z = zp[ip];

                SET_CONSTRAINTS(node->constraint, UNCONSTRAINED);
                node->myTag.domainID = dislocType;
                node->myTag.index = ip + offSet + param->nodeCount;

                AllocNodeArms(node, 2);

                if ((nbr1Index = ip + 1) >= np) nbr1Index = 0;
                if ((nbr2Index= ip - 1) < 0) nbr2Index = np - 1;

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = offSet + nbr1Index + param->nodeCount;
                node->burgX[0] = burg[burgIndex][0];
                node->burgY[0] = burg[burgIndex][1];
                node->burgZ[0] = burg[burgIndex][2];
                node->nx[0] = glidePlane[gpIndex][X];
                node->ny[0] = glidePlane[gpIndex][Y];
                node->nz[0] = glidePlane[gpIndex][Z];
            
                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = offSet + nbr2Index + param->nodeCount;
                node->burgX[1] = -burg[burgIndex][0];
                node->burgY[1] = -burg[burgIndex][1];
                node->burgZ[1] = -burg[burgIndex][2];
                node->nx[1] = glidePlane[gpIndex][X];
                node->ny[1] = glidePlane[gpIndex][Y];
                node->nz[1] = glidePlane[gpIndex][Z];

            }

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          Write blocks of nodal data to the data files as needed.
 */
            lastBlock = (ic == numChains - 1);

            if (inArgs->useSegmentedDataFile || lastBlock) {
                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);
                startRemeshIndex = 0;
            }

            offSet = inData->nodeCount;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     RhombohedralVaConfig
 *      Description:  Generate initial nodal data for dislocation loops
 *                    in rhombohedral vanadium (gamma-V).
 *
 *      Arguments:
 *          IN:  inArgs       Structure containing values associated with the
 *                            command line arguments.
 *          OUT: totDislocLen Total length of all dislocations inserted into
 *                            the configuration.
 *
 *      Values used from the <inArgs> structure:
 *          cubeLength    Length of cubic problem space in a single
 *                        dimension (units of b)
 *          numLoops      Number of rhombohedral loops to create
 *          interstitialLoops If set 1, all loops will be interstitial loops.
 *                        Otherwise loops will be vacancy loops.
 *          radius        Loop radius (in units of b)
 *          seed          Seed value for random number generator
 *          type          Type of dislocation configuration being generated.
 *
 * ordering !!! 
 * b1->b2->b3->b4  ==> b1->b4->b3->b2
 *
 *   For rhombohedral vanadium, there are 4 burgers vectors:
 *
 *    1.450000000000000e+000    1.450000000000000e+000    1.450000000000000e+000
 *    1.318500000000000e+000    1.318500000000000e+000   -1.187000000000000e+000
 *    1.318500000000000e+000   -1.187000000000000e+000    1.318500000000000e+000
 *   -1.187000000000000e+000    1.318500000000000e+000    1.318500000000000e+000
 * 
 *
 *   Each of the 4 burgers vectors has 3 glide planes as specified below:
 *
 *   -3.632975000000000e+000    3.632975000000000e+000                         0
 *    3.632975000000000e+000                         0   -3.632975000000000e+000
 *                         0   -3.632975000000000e+000    3.632975000000000e+000
 *
 *
 *   -3.632975000000000e+000    3.632975000000000e+000                         0
 *    3.294732499999999e-001   -3.303501750000000e+000   -3.303501750000000e+000
 *    3.303501750000000e+000   -3.294732499999999e-001    3.303501750000000e+000
 *
 *
 *    3.632975000000000e+000                         0   -3.632975000000000e+000
 *    3.294732499999999e-001   -3.303501750000000e+000   -3.303501750000000e+000
 *   -3.303501750000000e+000   -3.303501750000000e+000    3.294732499999999e-001
 *
 *
 *                         0   -3.632975000000000e+000    3.632975000000000e+000
 *    3.303501750000000e+000   -3.294732499999999e-001    3.303501750000000e+000
 *   -3.303501750000000e+000   -3.303501750000000e+000    3.294732499999999e-001
 *
 *-------------------------------------------------------------------------*/
void RhombohedralVaConfig(Home_t *home, InData_t *inData, InArgs_t *inArgs,
                          real8 *totDislocLen)
{
        int     i, id, loopIndex, burgIndex, dIndex, nextNode, newNodeIndex;
        int     minBurgIndex, maxBurgIndex, lastBlock, numSegs;
        int     startRemeshIndex = 0;
        int     nbr1Index;
        int     numLoops, useInterstitial, seed, dislocType;
        real8   cubeSize, radius;
        real8   ux, uy, uz;
        real8   loopCtrX,loopCtrY,loopCtrZ;
        real8   invSqrt6, twoInvSqrt6;
        real8   burg[7][3];
        real8   vec1[3], vec2[3], tmpVec[3];
        real8   tr[4][6][3];
        real8   tt0, tt1, tt2, td0, td1, td2;
        real8   ratio, rbmag, currentRadius, radius2;
        Param_t *param;
        Node_t  *node, *nbr1Node;


        param = inData->param;

/*
 *      Grap the command-line argument values needed for this function
 */
        cubeSize        = (real8)inArgs->cubeLength;
        numLoops        = inArgs->numLoops;
        radius          = inArgs->radius;
        useInterstitial = inArgs->interstitialLoops;
        seed            = inArgs->seed;
        dislocType      = inArgs->type;

/*
 *      All loops are now hexagonal and of [1 1 1] type burgers vectors for now
 */
        numSegs = 6;
        minBurgIndex = 0;
        maxBurgIndex = 3;

/*
 *      Define the sets of burger's vectors that may be used.  For each
 *      burgers vector we'll choose a glide plane normal vector which is
 *      perpendicular to the loop and along the burgers vector direction.
 */
        invSqrt6 = 1.0000000000e+000 / sqrt(6.0);
        twoInvSqrt6 = 2.000000000e+0000 * invSqrt6;

/*
 *      tt0 = burg_rhombo x 011 normal vector
3.797217736521221e-001    3.797217736521221e-001    8.435773522499124e-001
 */

        tt0 = 3.797217736521221e-001;
        tt1 = 3.797217736521221e-001;
        tt2 = 8.435773522499124e-001;


/*      td's are obtained from r_burgers vector cross skewed plane normal
 *        ex: plane_norm_2(2,:) is the second row of plane normal matrix!
 *        tdvector = cross(b2,plane_norm_2(2,:))
-7.995247017184457e-001    3.829661391680542e-001   -4.627063729761510e-001
 */

        td0 = 3.829661391680542e-001;
        td1 = 4.627063729761510e-001;
        td2 = 7.995247017184457e-001;



/*
   The following small burgers vectors are adjusted by the ratio
   between b_short and b_long (111 type)
   These vectors are not normal vectors. Should be careful when "burgmag" is used
*/
/* 
 *    1.450000000000000e+000    1.450000000000000e+000    1.450000000000000e+000
 *    1.318500000000000e+000    1.318500000000000e+000   -1.187000000000000e+000
 *    1.318500000000000e+000   -1.187000000000000e+000    1.318500000000000e+000
 *   -1.187000000000000e+000    1.318500000000000e+000    1.318500000000000e+000
*/
        burg[0][0] = 1.450000000000000e+000;
        burg[0][1] = 1.450000000000000e+000;
        burg[0][2] = 1.450000000000000e+000;

        rbmag = sqrt(burg[0][0]*burg[0][0] * 3.00000000000e+000); 
        
        for (i = 0; i < 3; i++) {
            burg[0][i] = burg[0][i] / rbmag;
        }

        burg[1][0] =  -1.187000000000000e+000 / rbmag;
        burg[1][1] =   1.318500000000000e+000 / rbmag;
        burg[1][2] =   1.318500000000000e+000 / rbmag;

        burg[2][0] =   1.318500000000000e+000 / rbmag;
        burg[2][1] =  -1.187000000000000e+000 / rbmag;
        burg[2][2] =   1.318500000000000e+000 / rbmag;

        burg[3][0] =   1.318500000000000e+000 / rbmag;
        burg[3][1] =   1.318500000000000e+000 / rbmag;
        burg[3][2] =  -1.187000000000000e+000 / rbmag;

/*
 *      Define the translation vectors as 4 X 6 X 3 array
 *
 *          tr[ ][ ][ ] - translation vector
 *             |  |  |
 *             |  |  x,y,z
 *             |  |
 *             |  6 directions
 *             |
 *             4 [111] type burgers vectors
 */

/* 
 *      If b cross l is outward it's vacancy; opposite for interstitial
 *      need to check to see if the order for these are consistent with
 *      loop vertix assignment Must be exactly the same.... mrhee
 *
 *      Note: When these values were calculated all burgers vectors were
 *            normalized
 */
        if (useInterstitial == 1) {
/*
 *          Interstitial [1 1 1] loops
 */
            tr[0][0][0]=  invSqrt6;    tr[0][0][1]=  invSqrt6;    tr[0][0][2]= -twoInvSqrt6; /* [ 1  1 -2] */
            tr[0][5][0]=  twoInvSqrt6; tr[0][5][1]= -invSqrt6;    tr[0][5][2]= -invSqrt6;    /* [ 2 -1 -1] */
            tr[0][4][0]=  invSqrt6;    tr[0][4][1]= -twoInvSqrt6; tr[0][4][2]=  invSqrt6;    /* [ 1 -2  1] */
            tr[0][3][0]= -invSqrt6;    tr[0][3][1]= -invSqrt6;    tr[0][3][2]=  twoInvSqrt6; /* [-1 -1  2] */
            tr[0][2][0]= -twoInvSqrt6; tr[0][2][1]=  invSqrt6;    tr[0][2][2]=  invSqrt6;    /* [-2  1  1] */
            tr[0][1][0]= -invSqrt6;    tr[0][1][1]=  twoInvSqrt6; tr[0][1][2]= -invSqrt6;    /* [-1  2 -1] */

/*
 *          Interstitial [-1 1 1] loops
 */
            tr[1][0][0]= -td1; tr[1][0][1]=  td0; tr[1][0][2]= -td2; /* [-1  1 -2]-gamma */
            tr[1][5][0]=  td1; tr[1][5][1]=  td2; tr[1][5][2]= -td0; /* [ 1  2 -1]-gamma */
            tr[1][4][0]=  tt2; tr[1][4][1]=  tt1; tr[1][4][2]=  tt1; /* [ 2  1  1]-bcc type  */
            tr[1][3][0]=  td1; tr[1][3][1]= -td0; tr[1][3][2]=  td2; /* [ 1 -1  2]-gamma */
            tr[1][2][0]= -td1; tr[1][2][1]= -td2; tr[1][2][2]=  td0; /* [-1 -2  1]-gamma */
            tr[1][1][0]= -tt2; tr[1][1][1]= -tt1; tr[1][1][2]= -tt1; /* [-2 -1 -1]-bcc type  */

/*
 *          Interstitial [1 -1 1] loops
 */
            tr[2][0][0]=  td2; tr[2][0][1]=  td1; tr[2][0][2]= -td0; /* [ 2  1 -1]-gamma */
            tr[2][5][0]=  td0; tr[2][5][1]= -td1; tr[2][5][2]= -td2; /* [ 1 -1 -2]-gamma */
            tr[2][4][0]= -tt0; tr[2][4][1]= -tt2; tr[2][4][2]= -tt0; /* [-1 -2 -1]-bcctype */
            tr[2][3][0]= -td2; tr[2][3][1]= -td1; tr[2][3][2]=  td0; /* [-2 -1  1]-gamma */
            tr[2][2][0]= -td0; tr[2][2][1]=  td1; tr[2][2][2]=  td2; /* [-1  1  2]-gamma */
            tr[2][1][0]=  tt0; tr[2][1][1]=  tt2; tr[2][1][2]=  tt1; /* [ 1  2  1]-bcctype */

/*
 *          Interstitial [1 1 -1] loops
 */
            tr[3][0][0]=  td0; tr[3][0][1]= -td2; tr[3][0][2]= -td1; /* [ 1 -2 -1]-gamma */
            tr[3][5][0]=  td2; tr[3][5][1]= -td0; tr[3][5][2]=  td1; /* [ 2 -1  1]-gamma */
            tr[3][4][0]=  tt0; tr[3][4][1]=  tt0; tr[3][4][2]=  tt2; /* [ 1  1  2]-bcctype */
            tr[3][3][0]= -td0; tr[3][3][1]=  td2; tr[3][3][2]=  td1; /* [-1  2  1]-gamma */
            tr[3][2][0]= -td2; tr[3][2][1]=  td0; tr[3][2][2]= -td1; /* [-2  1 -1]-gamma */
            tr[3][1][0]= -tt0; tr[3][1][1]= -tt0; tr[3][1][2]= -tt2; /* [-1 -1 -2]-bcctype */

        } else{
/*
 *          Vacancy [111] loops 
 */
            tr[0][0][0]=  invSqrt6;    tr[0][0][1]=  invSqrt6;    tr[0][0][2]= -twoInvSqrt6; /* [ 1  1 -2] */
            tr[0][1][0]=  twoInvSqrt6; tr[0][1][1]= -invSqrt6;    tr[0][1][2]= -invSqrt6;    /* [ 2 -1 -1] */
            tr[0][2][0]=  invSqrt6;    tr[0][2][1]= -twoInvSqrt6; tr[0][2][2]=  invSqrt6;    /* [ 1 -2  1] */
            tr[0][3][0]= -invSqrt6;    tr[0][3][1]= -invSqrt6;    tr[0][3][2]=  twoInvSqrt6; /* [-1 -1  2] */
            tr[0][4][0]= -twoInvSqrt6; tr[0][4][1]=  invSqrt6;    tr[0][4][2]=  invSqrt6;    /* [-2  1  1] */
            tr[0][5][0]= -invSqrt6;    tr[0][5][1]=  twoInvSqrt6; tr[0][5][2]= -invSqrt6;    /* [-1  2 -1] */

/*
 *          Vacancy [-1 1 1] loops
 */
            tr[1][0][0]= -td1; tr[1][0][1]=  td0; tr[1][0][2]= -td2; /* [-1  1 -2]-gamma */
            tr[1][1][0]=  td1; tr[1][1][1]=  td2; tr[1][1][2]= -td0; /* [ 1  2 -1]-gamma */
            tr[1][2][0]=  tt2; tr[1][2][1]=  tt1; tr[1][2][2]=  tt1; /* [ 2  1  1]-bcc type  */
            tr[1][3][0]=  td1; tr[1][3][1]= -td0; tr[1][3][2]=  td2; /* [ 1 -1  2]-gamma */
            tr[1][4][0]= -td1; tr[1][4][1]= -td2; tr[1][4][2]=  td0; /* [-1 -2  1]-gamma */
            tr[1][5][0]= -tt2; tr[1][5][1]= -tt1; tr[1][5][2]= -tt1; /* [-2 -1 -1]-bcc type  */

/*
 *          Vacancy [1 -1 1] loops
 */
            tr[2][0][0]=  td2; tr[2][0][1]=  td1; tr[2][0][2]= -td0; /* [ 2  1 -1]-gamma */
            tr[2][1][0]=  td0; tr[2][1][1]= -td1; tr[2][1][2]= -td2; /* [ 1 -1 -2]-gamma */
            tr[2][2][0]= -tt0; tr[2][2][1]= -tt2; tr[2][2][2]= -tt0; /* [-1 -2 -1]-bcctype */
            tr[2][3][0]= -td2; tr[2][3][1]= -td1; tr[2][3][2]=  td0; /* [-2 -1  1]-gamma */
            tr[2][4][0]= -td0; tr[2][4][1]=  td1; tr[2][4][2]=  td2; /* [-1  1  2]-gamma */
            tr[2][5][0]=  tt0; tr[2][5][1]=  tt2; tr[2][5][2]=  tt1; /* [ 1  2  1]-bcctype */

/*
 *          Vacancy [1 1 -1] loops
 */
            tr[3][0][0]=  td0; tr[3][0][1]= -td2; tr[3][0][2]= -td1; /* [ 1 -2 -1]-gamma */
            tr[3][1][0]=  td2; tr[3][1][1]= -td0; tr[3][1][2]=  td1; /* [ 2 -1  1]-gamma */
            tr[3][2][0]=  tt0; tr[3][2][1]=  tt0; tr[3][2][2]=  tt2; /* [ 1  1  2]-bcctype */
            tr[3][3][0]= -td0; tr[3][3][1]=  td2; tr[3][3][2]=  td1; /* [-1  2  1]-gamma */
            tr[3][4][0]= -td2; tr[3][4][1]=  td0; tr[3][4][2]= -td1; /* [-2  1 -1]-gamma */
            tr[3][5][0]= -tt0; tr[3][5][1]= -tt0; tr[3][5][2]= -tt2; /* [-1 -1 -2]-bcctype */
        }

/*
 *      FIX ME!  Need to modify the code to deal with <100> type burgers
 *               vector later
 */

/*
 *      ratio is for short radius of the rhombohedral structure
 *
 *      Note to developer: for details, see file 'iterate.m' and
 *                         rhombohedral_slipSystems.mcd
 */
        ratio   = 9.11568763163170430e-001;
        radius2 = radius * ratio;

/*
 *      Create one loop at a time, cycling through burgers vectors as we go.
 */
        inData->nodeCount = 0;
        nextNode = 0;
        burgIndex = maxBurgIndex;

        for (loopIndex = 0; loopIndex < numLoops; loopIndex++) {

            if (++burgIndex > maxBurgIndex) {
                burgIndex = minBurgIndex;
            }

/*
 *          Increase the size of the node array enough to hold
 *          all the new nodes created for this loop.
 */
            newNodeIndex = inData->nodeCount;
            inData->nodeCount += numSegs;
            inData->node = (Node_t *)realloc(inData->node,
                           inData->nodeCount * sizeof(Node_t));
            memset(&inData->node[newNodeIndex], 0, sizeof(Node_t) * numSegs);

            loopCtrX = (randm(&seed)-0.5) * cubeSize;
            loopCtrY = (randm(&seed)-0.5) * cubeSize;
            loopCtrZ = (randm(&seed)-0.5) * cubeSize;

            for (id = nextNode; id < (nextNode+numSegs); id++) {

                node =  &inData->node[id];

/*
 *              Pick the index for the direction component of the tr array
 */
                dIndex = (id - nextNode) % numSegs;

/*
 *              In order to acomodate the skewed structure of gamma-bcc
 *              vanadium, the radius of the hexagonal loop differs depending
 *              on which point we're dealing with.  For interstitials, index
 *              2 and 5 are the bcc type radius.  All others are all
 *              short-rhomobohedral radius*ratio.  For vacancies, since the
 *              loop sequence is reversed, index 1 and 4 will have the full
 *              radius and others will have radius2(=radius*ratio).
 */
                if (useInterstitial == 1) {
                    if ( (dIndex != 1 && dIndex != 4) && burgIndex !=0) {
                        currentRadius = radius2;
                    } else{
                        currentRadius = radius;
                    }
                } else{
                    if ( (dIndex != 2 && dIndex != 5) && burgIndex !=0) {
                        currentRadius = radius2;
                    } else{
                        currentRadius = radius;
                    }
                }

/*
 *              Set up the node
 */
                node->x = loopCtrX + currentRadius * tr[burgIndex][dIndex][0];
                node->y = loopCtrY + currentRadius * tr[burgIndex][dIndex][1];
                node->z = loopCtrZ + currentRadius * tr[burgIndex][dIndex][2];

                SET_CONSTRAINTS(node->constraint, UNCONSTRAINED);
                node->myTag.domainID = dislocType;
                node->myTag.index = id + param->nodeCount;

                AllocNodeArms(node, 2);

                node->burgX[0] = burg[burgIndex][X];
                node->burgY[0] = burg[burgIndex][Y];
                node->burgZ[0] = burg[burgIndex][Z];

                node->burgX[1] = -burg[burgIndex][X];
                node->burgY[1] = -burg[burgIndex][Y];
                node->burgZ[1] = -burg[burgIndex][Z];

                node->nbrTag[0].domainID = dislocType;
                node->nbrTag[0].index = param->nodeCount + nextNode +
                                        ((id+1)%numSegs);

                node->nbrTag[1].domainID = dislocType;
                node->nbrTag[1].index = param->nodeCount + nextNode +
                                        ((id+numSegs-1)%numSegs);
            }

/*
 *          Need to set the glide plane normal for each segment.
 *          Couldn't do it above because we didn't have all the
 *          neighbor node positions until all the nodes in the
 *          loop were created.
 */
            for (id = nextNode; id < (nextNode+numSegs); id++) {

                node = &inData->node[id];

                nbr1Index = node->nbrTag[1].index - param->nodeCount;
                nbr1Node = &inData->node[nbr1Index];

                ux = nbr1Node->x - node->x;
                uy = nbr1Node->y - node->y;
                uz = nbr1Node->z - node->z;

                Normalize(&ux,&uy,&uz);

/*
 *              In order to avoid some round-off issues, explicitly set
 *              any components of the line direction that are near-zero
 *              to be zero, then renormalize.
 */
                if (fabs(ux) < 1.0e-03) ux = 0.0;
                if (fabs(uy) < 1.0e-03) uy = 0.0;
                if (fabs(uz) < 1.0e-03) uz = 0.0;

                Normalize(&ux,&uy,&uz);

/*
 *              l cross b gives normal vector for each segment, they may
 *              not be (110)
 */
                vec1[0] = ux;
                vec1[1] = uy;
                vec1[2] = uz;

                vec2[0] =  node->burgX[1];
                vec2[1] =  node->burgY[1];
                vec2[2] =  node->burgZ[1];

                NormalizedCrossVector(vec1, vec2, tmpVec);

/*
 *              Again, change any near-zero components of the glideplane
 *              normal to zero then renormalize.
 */
                if (fabs(tmpVec[0]) < 1.0e-03) tmpVec[0] = 0.0;
                if (fabs(tmpVec[1]) < 1.0e-03) tmpVec[1] = 0.0;
                if (fabs(tmpVec[2]) < 1.0e-03) tmpVec[2] = 0.0;

                NormalizeVec(tmpVec);

                node->nx[1] = tmpVec[0];
                node->ny[1] = tmpVec[1];
                node->nz[1] = tmpVec[2];

                nbr1Node->nx[0] = tmpVec[0];
                nbr1Node->ny[0] = tmpVec[1];
                nbr1Node->nz[0] = tmpVec[2];

            }  /*  for (id = nextNode; ...) */

            nextNode += numSegs;

/*
 *          The initial segments created are not necessarily limited to
 *          param->maxSegLen, so a call to InitRemesh() is needed to
 *          chop any excessively long segments into proper lengths.
 *          Write blocks of nodal data to the data files as needed.
 */
            lastBlock = (loopIndex == (numLoops-1));

            if (inArgs->useSegmentedDataFile || lastBlock) {
                InitRemesh(inData, dislocType, startRemeshIndex);
                IncDislocationDensity(inData, totDislocLen);
                WriteInitialNodeData(home, inData, inArgs, lastBlock);
                startRemeshIndex = 0;
                nextNode = 0;
            }

        }  /* for (loopIndex = 0; ...) */

        return;
}
