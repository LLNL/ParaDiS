/**************************************************************************
 *
 *      Module:       ReadEshelbyInclusions.c
 *
 *      Description:  This module contains the functions for reading
 *                    the Eshelby inclusion data file and distributing
 *                    the inclusion data to all the processors.
 *
 *      Included public functions:
 *          ReadEshelbyInclusions()
 *
 *      Included private functions:
 *          DomainOwnsInclusion()
 *
 *************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"

#define INC_INPUT_BLK_SIZE   10000
#define INC_ALLOC_INCREMENT  1000

#ifdef ESHELBY

/*---------------------------------------------------------------------------
 *
 *      Function:     AdjustInitialCoordinates
 *
 *      Description:  Move any coordinate that is exactly on the maximum
 *                    boundary of a periodic surfaace to the opposite
 *                    boundary.  The rationale for this  is that the physical
 *                    location of Eshelby inclusions are static, and as such,
 *                    the cell containing them is statically set during
 *                    initialization.  After domain boundaries shift during
 *                    rebalance, both nodes and inclusions are migrated to
 *                    the domain which now encompasses them.  The migration
 *                    process shifts coordinates on the maximum boundary of
 *                    a periodic surface to the opposite surface.  This is
 *                    fine for nodes, because the node's cell is recomputed,
 *                    but for inclusions this causes major problems since
 *                    the inclusions cell should *never* change.
 *
 *
 *      Arguments:
 *          coord     3 element array containing the X, Y and Z coordinates
 *                    of the center of an inclusion.
 *
 *-------------------------------------------------------------------------*/
static void AdjustInitialCoordinates(Param_t *param, real8 coord[3])
{
        real8 xAdjustment, yAdjustment, zAdjustment;

        xAdjustment = param->maxSideX - param->minSideX;
        yAdjustment = param->maxSideY - param->minSideY;
        zAdjustment = param->maxSideZ - param->minSideZ;

/*
 *      Adjustments become zero if periodic boundaries are not enabled.
 */
        xAdjustment *= (real8)(param->xBoundType == Periodic);
        yAdjustment *= (real8)(param->yBoundType == Periodic);
        zAdjustment *= (real8)(param->zBoundType == Periodic);


        if (coord[X] < param->minSideX) {
            coord[X] += xAdjustment;
        } else if (coord[X] >= param->maxSideX) {
            coord[X] -= xAdjustment;
        }

        if (coord[Y] < param->minSideY) {
            coord[Y] += yAdjustment;
        } else if (coord[Y] >= param->maxSideY) {
            coord[Y] -= yAdjustment;
        }

        if (coord[Z] < param->minSideZ) {
            coord[Z] += zAdjustment;
        } else if (coord[Z] >= param->maxSideZ) {
            coord[Z] -= zAdjustment;
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Author:         Gregg Hommes
 *
 *      Function:       DomainOwnsInclusion
 *
 *      Description:    
 *
 *      Arguments:
 *          position  3 element array containing the X, Y and Z coordinates
 *                    of the center of an inclusion.
 *
 *      Last Modified:  02/13/2008 - original version
 *
 *-------------------------------------------------------------------------*/
static int DomainOwnsInclusion(Home_t *home, real8 position[3])
{
        Param_t *param;

        param = home->param;

        if (((position[X] < home->domXmin) ||
             (position[X] > home->domXmax)) ||
            ((position[Y] < home->domYmin) ||
             (position[Y] > home->domYmax)) ||
            ((position[Z] < home->domZmin) ||
             (position[Z] > home->domZmax))) {
            return(0);
        }

/*
 *      If the position is on a lower boundary of the domain, the
 *      domain owns it.  If the position is on an upper boundary
 *      of the domain, the domain only owns it if the domain
 *      boundary is also a problem space boundary.
 */
        if (((position[X] == home->domXmax) && 
             (home->domXmax < param->maxSideX)) ||
            ((position[Y] == home->domYmax) && 
             (home->domYmax < param->maxSideY)) ||
            ((position[Z] == home->domZmax) && 
             (home->domZmax < param->maxSideZ))) {
            return(0);
        }

        return(1);
}


/*---------------------------------------------------------------------------
 *
 *      Author:         Gregg Hommes
 *
 *      Function:       ReadEshelbyInclusions
 *
 *      Description:    The format for the particles files is as follows:
 *                      Inclusion Id
 *                      Position of the particle (x, y, z)
 *                      Semi-principal axes  (a, b, c)
 *                      Two vectors for rotation matrix
 *                      Strain field for Eshelby calculations, S[0]-S[5]
 *
 *      Arguments:
 *          ctrlFileName  Name of the control file to be read.
 *
 *      Last Modified:                     02/13/2008 - original version
 *      Extended to reading ellipsoids     06/16        S. Aubry
 *
 *-------------------------------------------------------------------------*/
void ReadEshelbyInclusions(Home_t *home)
{
        int          i, m, n;
        int          maxTokenLen, tokenType;
        int          allocSize, allocedIncs, numRead, numIncs;
        int          broadcastSize, distIncomplete;
        real8        strainVals[6];
        char         token[256];
        FILE         *fpEshelby;
        Param_t      *param;
        EInclusion_t *incs, *tmpIncs;


        param = home->param;

        if (param->enableInclusions == 0) {
            return;
        }

        if (home->myDomain == 0) {
            fpEshelby = fopen(param->inclusionFile, "r");
            if (fpEshelby == (FILE *)NULL) {
                Fatal("ReadEshelbyInclusions: Error %d opening file %s",
                      errno, param->inclusionFile);
            }
        }

        maxTokenLen = sizeof(token);
        tokenType = TOKEN_GENERIC;

        allocSize = INC_INPUT_BLK_SIZE * sizeof(EInclusion_t);
        incs = (EInclusion_t *)malloc(allocSize);

        distIncomplete = 1;

        allocedIncs = 0;
        numIncs = 0;
        tmpIncs = (EInclusion_t *)NULL;


/*
 *      All processes loop until all inclusions are read in and
 *      distributed to the owning domains.
 */
        while (distIncomplete) {

            numRead = 0;
/*
 *          Only domain zero reads in the inclusion data
 */
            if (home->myDomain == 0) {

                while (1) {
/*
 *                  First token associated with the inclusion is the ID.  If
 *                  not found, we're done.
 */
                    tokenType = GetNextToken(fpEshelby, token, maxTokenLen);
                    if (tokenType == TOKEN_NULL) {
                        break;
                    }

                    incs[numRead].id = atoi(token);

/*
 *                  Now read the coordinates of the inclusions
 */
                    for (i = 0; i < 3; i++) {
                        tokenType = GetNextToken(fpEshelby, token, maxTokenLen);
                        if ((tokenType == TOKEN_ERR) ||
                            (tokenType == TOKEN_NULL)) {
                            Fatal("ReadEshelbyInclusions: "
                                  "Unexpected EOF in file %s at inclusion %d",
                                  param->inclusionFile, numRead);
                        }

                        incs[numRead].position[i] = atof(token);
                    }

/*
 *                  If the inclusion is outside any of the periodic boundaries
 *                  we need to shift it back it.  
 *
 *                  FIX ME!  If periodic boundaries are not enabled and the
 *                  inclusion is outside the surface or inside but less than
 *                  inclusion.radius distance from any surface, abort with
 *                  an error?
 */
                    FoldBox(param, &incs[numRead].position[X],
                                   &incs[numRead].position[Y],
                                   &incs[numRead].position[Z]);


/*
 *                  Since the inclusion's cell never changes once it has been
 *                  set, we have to make sure we shift any inclusions on the
 *                  maximum boundary of a periodic surface to the opposite
 *                  boundary to make sure it gets put into the proper domain
 *                  and cell right at the start.
 */
                    AdjustInitialCoordinates(param, incs[numRead].position);


/*
 *                  Next read the semi-principal axes of the inclusions
 */
                    for (i = 0; i < 3; i++) {
                       tokenType = GetNextToken(fpEshelby, token, maxTokenLen);
                       if ((tokenType == TOKEN_ERR) || (tokenType == TOKEN_NULL)) {
                          Fatal("ReadEshelbyInclusions: "
                                "Unexpected EOF in file %s at inclusion %d",
                                param->inclusionFile, numRead);
                       }
                       
                       incs[numRead].radius[i] = atof(token);
                    }

/*
 *                  Next read the two vectors of the rotation matrix.
 *                  These vectors are normalized here.
 */
                    real8 vec[3];
                    for (i = 0; i < 2; i++) 
                    {
                       for (m = 0; m < 3; m++) 
                       {
                          tokenType = GetNextToken(fpEshelby, token, maxTokenLen);
                          if ((tokenType == TOKEN_ERR) || (tokenType == TOKEN_NULL)) {
                             Fatal("ReadEshelbyInclusions: "
                                   "Unexpected EOF in file %s at inclusion %d",
                                   param->inclusionFile, numRead);
                          }
                          
                          vec[m] = atof(token);
                       }
                       NormalizeVec(vec);
                       VECTOR_COPY(incs[numRead].rotation[i],vec);                      
                    }


/*
 *                  Lastly, read in the strain field for the inclusion. 
 */
                    for (m = 0; m < 6; m++) 
                    {
                       tokenType = GetNextToken(fpEshelby, token,
                                                maxTokenLen);
                       if ((tokenType == TOKEN_ERR) ||
                           (tokenType == TOKEN_NULL)) 
                       {
                          Fatal("ReadEshelbyInclusions: Unexpected "
                                "EOF in file %s at inclusion %d",
                                param->inclusionFile, numRead);
                       }
                       strainVals[m]= atof(token);
                    }

                    for (m = 0; m < 6; m++) 
                       incs[numRead].strainField[m] = strainVals[m];


/*
 *                  Put the inclusion into the proper cell
 */
                    LocateCell(home, &incs[numRead].cellID,
                               incs[numRead].position);

/*
 *                  If we've read in the maximum number of inclusions,
 *                  broadcast them to the appropriate domains, then go
 *                  back and read in another block of inclusions from
 *                  the file.
 */
                    numRead++;

                    if (numRead >= INC_INPUT_BLK_SIZE) break;


                }  /* while (1) */

            }  /* if (home->myDomain == 0) */

#ifdef PARALLEL
/*
 *          Broadcast the number of inclusions read in
 */
            MPI_Bcast(&numRead, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

/*
 *          If we have inclusions to broadcast, do so
 */
            if (numRead > 0) {

#ifdef PARALLEL
                broadcastSize = numRead * sizeof(EInclusion_t);
                MPI_Bcast((char *)incs, broadcastSize, MPI_CHAR,
                          0, MPI_COMM_WORLD);
#endif
/*
 *              Process the inclusions:  Scan through, locate all that
 *              are local to this domain.  Add those inclusions to the local
 *              inclusion list, ignore the others.
 */
                for (i = 0; i < numRead; i++) {
                    
                    if (!DomainOwnsInclusion(home, incs[i].position)) {
                        continue;
                    }

/*
 *                  If we need more space on the local array of inclusions,
 *                  reallocate the array.
 */
                    if (numIncs >= allocedIncs) {
                        allocedIncs += INC_ALLOC_INCREMENT;
                        allocSize = allocedIncs * sizeof(EInclusion_t);
                        tmpIncs = (EInclusion_t *)realloc(tmpIncs, allocSize);
                    }

                    tmpIncs[numIncs].id = incs[i].id;
                    tmpIncs[numIncs].cellID = incs[i].cellID;

                    for (m = 0; m < 3; m++) {

                        tmpIncs[numIncs].position[m] = incs[i].position[m];
                        tmpIncs[numIncs].radius[m] = incs[i].radius[m];

                        for (n = 0; n < 2; n++) {
                           tmpIncs[numIncs].rotation[n][m] = incs[i].rotation[n][m];
                        }
                    }

                    for (n = 0; n < 6; n++) 
                       tmpIncs[numIncs].strainField[n] = incs[i].strainField[n];
                   

/*
 *                  We don't need to add the inclusion to the cell
 *                  inclusion queue here it will be done during the
 *                  call to SortNativeNodes() later.
 */
                    tmpIncs[numIncs].nextInCell = -1;

                    numIncs++;

                }

                numRead = 0;

            } else {
                distIncomplete = 0;
            }

        }  /* while (distIncomplete) */


        if (home->myDomain == 0) {
            fclose(fpEshelby);
        }

        free(incs);

/*
 *      Truncate the inclusion array to the proper size and
 *      store it in the 'home' structure for later use.
 */
        if ((numIncs > 0) && (allocedIncs != numIncs)) {
            tmpIncs = (EInclusion_t *)realloc(tmpIncs, numIncs *
                                              sizeof(EInclusion_t));
        }

        home->locInclusionCount = numIncs;
        home->totInclusionCount = numIncs;
        home->eshelbyInclusions = tmpIncs;


        return;
}

#endif // ESHELBY
