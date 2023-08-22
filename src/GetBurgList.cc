/***************************************************************************
 *
 *      Module:       GetBurgList.c
 *      Description:  Contains a generic dispatch function to call a module
 *                    appropriate to the material type to calculate the
 *                    available burgers vectors and associated glide planes.
 *
 **************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

void GetBurgList(Home_t *home)
{
    switch(home->param->materialType) 
    {
        case MAT_TYPE_BCC:
        {
            GetBurgList_BCC ( &home->burgData.burgList                 ,
                              &home->burgData.planeList                ,
                              &home->burgData.numPlanesPerBurg         ,
                              &home->burgData.burgFirstPlaneIndex      ,
                              &home->burgData.numGlissilePlanesPerBurg ,
                              &home->burgData.planeType                ,
                              &home->burgData.numBurgVectors           ,
                              &home->burgData.numPlanes                 ); 
            break;
        }

        case MAT_TYPE_FCC:
        {
            GetBurgList_FCC ( &home->burgData.burgList                 ,
                              &home->burgData.planeList                ,
                              &home->burgData.numPlanesPerBurg         ,
                              &home->burgData.burgFirstPlaneIndex      ,
                              &home->burgData.numGlissilePlanesPerBurg ,
                              &home->burgData.planeType                ,
                              &home->burgData.numBurgVectors           ,
                              &home->burgData.numPlanes                 ); 
            break;
        }

        case MAT_TYPE_HCP:
        {
            GetBurgList_HCP (  home->param->cOVERa                     ,
                              &home->burgData.burgList                 ,
                              &home->burgData.planeList                ,
                              &home->burgData.splinterableBurgList     ,
                              &home->burgData.numPlanesPerBurg         ,
                              &home->burgData.burgFirstPlaneIndex      ,
                              &home->burgData.numGlissilePlanesPerBurg ,
                              &home->burgData.planeType                ,
                              &home->burgData.numBurgVectors           ,
                              &home->burgData.numPlanes                ,
                              &home->burgData.numSplinterableBurgs      ); 
            break;
        }
    }
}
