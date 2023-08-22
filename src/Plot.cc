//---------------------------------------------------------------------------
//
// Function:     Plot
// Description:  Plot the current dislocations from the specified
//               domain in an X-window display.
//
//---------------------------------------------------------------------------

#ifndef NO_XWINDOW

#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Init.h"
#include "Home.h"
#include "Util.h"
#include "DisplayC.h"
#include "Decomp.h"
#include "V3.h"
#include "M33.h"
#include "M44.h"
#include "Vertex.h"
#include "Triangle.h"
#include "Ellipsoid.h"

// DSCLEN should correspond to the value in display.h

#define DSCLEN 60

#ifdef ESHELBY
// Draw_Ellipsoid()
//
// Will draw an ellipsoid at a given position with scaling and rotation transform.
// 
// Note - the vertex and triangle list of a generic sphere are created the 
// first time this code is executed. Once created, the matrix transformation is applied
// to each vertex and triangle in the generic sphere.
//
// Note - currently only active if Eshelby inclusions are currently active
//----------------------------------------------------------------------------------------------
                
static void Draw_Ellipsoid (real8 *m33, real8 *pos, unsigned long color, double r, unsigned long attr)
{
   static Ellipsoid e;

   if (!e.Vertices())
      e.Tesselate(1);

   e.Display(m33,pos,color,r,attr);
}
#endif  // ESHELBY

//---------------------------------------------------------------------------
//
//      Function:     Plot
//      Description:  Plot the current dislocations from the specified
//                    domain in an X-window display.
//
//                    Note: This function is intended to be called
//                    once per domain each time the X-windows plot is
//                    to updated.  The first and last calls per update
//                    should set the FIRST_BLOCK and LAST_BLOCK bits in
//                    the blkFlag field so the function does the proper
//                    initialization and cleanup required for this set
//                    of calls.
//      Args:
//          domIndex  specifies domain whose information should
//                    be plotted.
//          blkFlag   treated as a bit field. Used to determine
//                    if initialization or cleanup processing
//                    needs to be done in this call.
//
//---------------------------------------------------------------------------

void Plot(Home_t *home, int domIndex, int blkFlag) 
{
    char  string[DSCLEN*10];

    static unsigned color;
    static real8    Lx, Ly, Lz, xmin, ymin, zmin, xmax, ymax, zmax;
    static real8    Lmax, a, b, c;
    
    if (home->myDomain != 0) return;

    Param_t *param = home->param;

    // If this is the first block (domain) of data being plotted during
    // this update, do some basic initialization and setup.

    if (blkFlag & FIRST_BLOCK) 
    {
        WinLock();
        WinClear();

        xmin=param->minSideX;  xmax=param->maxSideX;
        ymin=param->minSideY;  ymax=param->maxSideY;
        zmin=param->minSideZ;  zmax=param->maxSideZ;
    
        Lx=xmax-xmin;
        Ly=ymax-ymin;
        Lz=zmax-zmin;

        Lmax = Lx;
        if (Lmax<Ly) Lmax=Ly;
        if (Lmax<Lz) Lmax=Lz;

        a = Lx/Lmax;
        b = Ly/Lmax;
        c = Lz/Lmax;

        //printf("xmin=%8.1lf xmax=%8.1lf\n",xmin,xmax);
        //printf("ymin=%8.1lf ymax=%8.1lf\n",ymin,ymax);
        //printf("zmin=%8.1lf zmax=%8.1lf\n",zmin,zmax);
        //printf("Lx=%8.1lf Ly=%8.1lf Lz=%8.1lf Lmax=%8.1lf\n",Lx,Ly,Lz,Lmax);
        //printf("a=%8.1lf b=%8.1lf c=%8.1lf\n",a,b,c);
        //printf("pbcshift=[%8.1lf b=%8.1lf c=%8.1lf]\n",pbcshift[0],pbcshift[1],pbcshift[2]);
        //exit(0);
    
        color=colors[9];

        // Using the domain decomposition data, plot all the domain boundaries.

        XPlotDecomp(home, xmin, ymin, zmin, Lmax, color, line_width);

        // Plot the free surfaces if periodic boundary conditions are
        // not enabled.

        if (    (param->zBoundType==Free)
             || (param->yBoundType==Free)
             || (param->xBoundType==Free) )
        {
            real8 x1=(param->xBoundMin-xmin)/Lmax*2-1;
            real8 x2=(param->xBoundMax-xmin)/Lmax*2-1;

            real8 y1=(param->yBoundMin-ymin)/Lmax*2-1;
            real8 y2=(param->yBoundMax-ymin)/Lmax*2-1;

            real8 z1=(param->zBoundMin-zmin)/Lmax*2-1;
            real8 z2=(param->zBoundMax-zmin)/Lmax*2-1;

            WinDrawLine(x1,y1,z1,x2,y1,z1, color, line_width/2,0);
            WinDrawLine(x1,y2,z1,x2,y2,z1, color, line_width/2,0);
            WinDrawLine(x1,y1,z2,x2,y1,z2, color, line_width/2,0);
            WinDrawLine(x1,y2,z2,x2,y2,z2, color, line_width/2,0); 

            WinDrawLine(x1,y1,z1,x1,y2,z1, color, line_width/2,0);
            WinDrawLine(x2,y1,z1,x2,y2,z1, color, line_width/2,0);
            WinDrawLine(x1,y1,z2,x1,y2,z2, color, line_width/2,0);
            WinDrawLine(x2,y1,z2,x2,y2,z2, color, line_width/2,0); 

            WinDrawLine(x1,y1,z1,x1,y1,z2, color, line_width/2,0);
            WinDrawLine(x1,y2,z1,x1,y2,z2, color, line_width/2,0);
            WinDrawLine(x2,y1,z1,x2,y1,z2, color, line_width/2,0);
            WinDrawLine(x2,y2,z1,x2,y2,z2, color, line_width/2,0); 
        }

    }  // if (blkFlag && FIRST_BLOCK)

    // plotting segments

    int newNodeKeyPtr=0;

    if (domIndex==0) { newNodeKeyPtr=home->newNodeKeyPtr; }
    else             { newNodeKeyPtr=home->mirrorDomainKeys[domIndex]->newNodeKeyPtr; }

    for(int i=0; i<newNodeKeyPtr; i++) 
    {
        // printf("domain = %d  key = %d\n", domIndex, i);

        Node_t *node=0;

        if (domIndex==0) { node=home->nodeKeys[i]; }
        else             { node=home->mirrorDomainKeys[domIndex]->nodeKeys[i]; }
           
        if (node==0) continue;

        real8 x = node->x; 
        real8 y = node->y; 
        real8 z = node->z;

        // new way of plotting, does not assume PBC

        x=(x-xmin)/Lx; x-=0.5+pbcshift[0]; x*=2;
        y=(y-ymin)/Ly; y-=0.5+pbcshift[1]; y*=2;
        z=(z-zmin)/Lz; z-=0.5+pbcshift[2]; z*=2; 
            
        for(int j=0; (j<node->numNbrs); j++) 
        {
           int index = node->nbrTag[j].index;

           if (index < 0) continue;

           if ((node->nbrTag[j].domainID == domIndex) &&

           (index < i)) continue;

           real8 x2=0,y2=0,z2=0;
           GetNbrCoords(home,node,j,&x2,&y2,&z2);

           // printf("Plot: neighbor[%d] (%d.%d)  %e %e %e\n", j, node->nbrTag[j].domainID, node->nbrTag[j].index, x2, y2, z2);
           // fflush(stdout);

           // new way of plotting, does not assume PBC

           x2 =(x2-xmin)/Lx; x2-=0.5+pbcshift[0]; x2*=2; 
           y2 =(y2-ymin)/Ly; y2-=0.5+pbcshift[1]; y2*=2; 
           z2 =(z2-zmin)/Lz; z2-=0.5+pbcshift[2]; z2*=2;
                
           color=colors[1];

           if (color_scheme==1) 
           {
               // color by domain

               if(domIndex==node->nbrTag[j].domainID)
                   color=colors[domIndex%8+2];
               else
                   color=colors[1];
           } 

           if (color_scheme==2) 
           {
               // color by Burgers vector

               color = colors[((int) (fabs(node->burgX[j])*4 +
                                      fabs(node->burgY[j])*2 +
                                      fabs(node->burgZ[j])))%8+2];
           }

           // do not draw segments across PBC

           if ((fabs(x-x2)<=1) && (fabs(y-y2)<=1) && (fabs(z-z2)<=1)) 
           {
               WinDrawLine(x*a,y*b,z*c,x2*a,y2*b,z2*c,color,line_width,0);
           }
        }
    }

    // plotting nodes

    if (domIndex==0) { newNodeKeyPtr=home->newNodeKeyPtr; }
    else             { newNodeKeyPtr=home->mirrorDomainKeys[domIndex]->newNodeKeyPtr; }

    for(int i=0; i<newNodeKeyPtr; i++) 
    {
        Node_t *node=0;

        if (domIndex==0) { node=home->nodeKeys[i]; }
        else             { node=home->mirrorDomainKeys[domIndex]->nodeKeys[i]; }
            
        if (node==0) continue;

        real8 x=node->x; 
        real8 y=node->y; 
        real8 z=node->z;
  
        // new way of plotting, does not assume PBC

        x=(x-xmin)/Lx; x-=0.5+pbcshift[0]; x*=2;
        y=(y-ymin)/Ly; y-=0.5+pbcshift[1]; y*=2;
        z=(z-zmin)/Lz; z-=0.5+pbcshift[2]; z*=2;
            
        sprintf(string,"(%d,%d)%d(%1.3e,%1.3e,%1.3e)", node->myTag.domainID, node->myTag.index,
                                                       node->numNbrs, node->x,node->y,node->z);
            
        // Note: for these points, the radius is specified as
        // a number of pixels in the display.

        WinDrawPointS(x*a,y*b,z*c,point_radius,1,Lmax, colors[0], 4,string);
    }

#ifdef ESHELBY
    // Also plot any Eshelby inclusions in the simulation.
    // In the current implentation of inclusions, each
    // processor knows about all the inclusions, so we
    // can plot them here in one shot.
    //
    // WARNING: For inclusions from the remote domains, the only
    // data available are id, radius and position!

    int             numInclusions=0;
    EInclusion_t   *incList      =0;

    if (domIndex==0) { numInclusions = home->locInclusionCount;
                       incList       = home->eshelbyInclusions; } 
    else             { numInclusions = home->mirrorDomainKeys[domIndex]->inclusionCount;
                       incList       = home->mirrorDomainKeys[domIndex]->inclusionList; }
        
    for (int i=0; (i<numInclusions); i++) 
    {
        EInclusion_t *inclusion = &incList[i];

        real8 x = inclusion->position[X];
        real8 y = inclusion->position[Y];
        real8 z = inclusion->position[Z];
        
        // new way of plotting, does not assume PBC

        x=(x-xmin)/Lx; x-=0.5+pbcshift[0]; x*=2;
        y=(y-ymin)/Ly; y-=0.5+pbcshift[1]; y*=2;
        z=(z-zmin)/Lz; z-=0.5+pbcshift[2]; z*=2;
        
        sprintf(string,"Inclusion (%d) r=%.2lf pos=%1.3e,%1.3e,%1.3e", 
                inclusion->id,
                inclusion->radius[0]  , inclusion->position[X],
                inclusion->position[Y], inclusion->position[Z]);
        
        // For plotting the inclusions, we're just using the
        // existing mechanism for plotting the node points
        // with one difference, however.  In this case, we
        // want the radius to be treated as a distance in
        // units of b rather than pixels on the screen.
 
        // build up the transform from unit sphere to ellipsoid..
        
            real8 pos [3] = {  x,  y,  z };
        
            real8 ra = 2.0*inclusion->radius[0]/Lmax;  // ra = x-axis scale factor
            real8 rb = 2.0*inclusion->radius[1]/Lmax;  // rb = y-axis scale factor
            real8 rc = 2.0*inclusion->radius[2]/Lmax;  // rc = z-axis scale factor
        
            real8 vec[3];
            real8 m33[9];
            cross(inclusion->rotation[0],inclusion->rotation[1],vec);
            NormalizeVec(vec);

            // m33 is the transpose of the rotation matrix
            m33[0] = inclusion->rotation[0][X];
            m33[3] = inclusion->rotation[0][Y];
            m33[6] = inclusion->rotation[0][Z];

            m33[1] = inclusion->rotation[1][X];
            m33[4] = inclusion->rotation[1][Y];
            m33[7] = inclusion->rotation[1][Z];

            m33[2] = vec[X];
            m33[5] = vec[Y];
            m33[8] = vec[Z];
        
            real8 mscale[9]; 
            M33_SCALE(mscale,ra,rb,rc);
            M33_MUL2(m33,mscale); 
            
            Draw_Ellipsoid(m33,pos,colors[4],line_width,0);
    }
#endif  // ESHELBY

    // If this is the last block (domain) of data being plotted during
    // this update, draw the frame and do basic cleanup.

    if (blkFlag & LAST_BLOCK) 
    {
        WinDrawLine(-a,-b,-c,-a,-b, c,colors[10],line_width,0);
        WinDrawLine(-a,-b, c,-a, b, c,colors[10],line_width,0);
        WinDrawLine(-a, b, c,-a, b,-c,colors[10],line_width,0);
        WinDrawLine(-a, b,-c,-a,-b,-c,colors[10],line_width,0);
        WinDrawLine( a,-b,-c, a,-b, c,colors[10],line_width,0);
        WinDrawLine( a,-b, c, a, b, c,colors[10],line_width,0);
        WinDrawLine( a, b, c, a, b,-c,colors[10],line_width,0);
        WinDrawLine( a, b,-c, a,-b,-c,colors[10],line_width,0);
        WinDrawLine(-a,-b,-c, a,-b,-c,colors[10],line_width,0);
        WinDrawLine(-a,-b, c, a,-b, c,colors[10],line_width,0);
        WinDrawLine(-a, b, c, a, b, c,colors[10],line_width,0);
        WinDrawLine(-a, b,-c, a, b,-c,colors[10],line_width,0);

        WinUnlock();
        WinRefresh();
    }
}

#endif // NO_XWINDOW
