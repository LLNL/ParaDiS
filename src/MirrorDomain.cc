#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Node.h"
#include "Comm.h"
#include "Util.h"
#include "QueueOps.h"
#include "Tag.h"
#include "MirrorDomain.h"

//-------------------------------------------------------------------------
// Function:       MirrorDomain_Free
// Description:    Deallocate all storage associated with the mirror
//                 domain indicated by <domINdex>
//-------------------------------------------------------------------------

void MirrorDomain_Free (Home_t *home, int domIndex)
{
    if (home->mirrorDomainKeys==0) return;

    MirrorDomain_t *md = home->mirrorDomainKeys[domIndex];

    if (md && md->newNodeKeyPtr==0) { free(md); return; }   

    if (md)
    {
        for (int i=0; (i<md->newNodeKeyPtr); i++) 
        {   
            Node_t *node = md->nodeKeys[i];
    
            if (node)
            {   
                if (node->armCoordIndex) { free(node->armCoordIndex); node->armCoordIndex=0; }
    
                PushFreeNodeQ(home,node); // (recycle node)
            }   
        }   
    
#ifdef ESHELBY
        if (md->inclusionList && (md->inclusionCount>0))
        { free(md->inclusionList); md->inclusionList=0; md->inclusionCount=0; }
#endif
    
        if (md->armX    ) { free(md->armX    ); md->armX    =0; }
        if (md->armY    ) { free(md->armY    ); md->armY    =0; }
        if (md->armZ    ) { free(md->armZ    ); md->armZ    =0; }
        if (md->nodeKeys) { free(md->nodeKeys); md->nodeKeys=0; }

        free(md);
    }
}

//-------------------------------------------------------------------------

void MirrorDomain_Print(FILE *strm, MirrorDomain_t *mdom)
{
   if (strm && mdom)
   {
       fprintf(strm, "mirror_domain>\n");
       fprintf(strm, "{\n");
       fprintf(strm, "   nodeKeys       = 0x%08lx\n", (unsigned long) mdom->nodeKeys       );
       fprintf(strm, "   newNodeKeyPtr  = %d\n"     ,                 mdom->newNodeKeyPtr  );
       fprintf(strm, "   armX           = 0x%08lx\n", (unsigned long) mdom->armX           );
       fprintf(strm, "   armY           = 0x%08lx\n", (unsigned long) mdom->armY           );
       fprintf(strm, "   armZ           = 0x%08lx\n", (unsigned long) mdom->armZ           );
#ifdef ESHELBY
       fprintf(strm, "   inclusionCount = %d\n"     ,                 mdom->inclusionCount );
       fprintf(strm, "   inclusionList  = 0x%08lx\n", (unsigned long) mdom->inclusionList  );
#endif
       fprintf(strm, "}\n");
   }
}

