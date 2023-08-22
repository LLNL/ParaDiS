/***************************************************************************
 *   
 *      Module:     DeltaPlasticStrain.
 *
 *      Description:  This contains a simple generic dispatch function that
 *                    will invoke the version of DeltaPlasticStrain*() 
 *                    appropriate to the type of material being simulated
 *
 ***************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "Mobility.h"

void DeltaPlasticStrain(Home_t *home)
{
    switch(home->param->materialType) 
    {
        case MAT_TYPE_BCC:
        {
            switch(home->param->bcc_DensityFluxDecomp) 
            {
                case 1 : { DeltaPlasticStrain_BCC (home); break; }
                case 2 : { DeltaPlasticStrain_BCC2(home); break; }
                default: { DeltaPlasticStrain_BCC (home); break; }
            }
            break;
        }
        case MAT_TYPE_FCC            : { DeltaPlasticStrain_FCC     (home); break; }
        case MAT_TYPE_HCP            : { DeltaPlasticStrain_HCP     (home); break; }
        case MAT_TYPE_RHOMBOHEDRAL_VA: { DeltaPlasticStrain_rhomboVa(home); break; }
    }

    return;
}
