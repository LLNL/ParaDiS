/*
  display.cpp
  by Wei Cai and Dongyi Liao  caiwei@mit.edu, liaody@mit.edu
  Last Modified : Tue Apr 29 11:21:17 2003

  FUNCTION  :  X-window display package

  Problem:  rotation matrix A can become non-unitary
*/
#ifndef NO_XWINDOW

#include <stdio.h>

#include "display.h"
#ifdef TMCC
extern "C"
#endif

unsigned long colorconvert_16[3][6]={{0xF8,8,0,0xF800,11,0x1F},
                                     {0xF8,3,0,0x7C0,  6,0x1F},
                                     {0xF8,0,0,0x1F,   0,0x1F}};
unsigned long colorconvert_32[3][6]={{0xFF,16,0,0xFF0000,16,0xFF},
                                     {0xFF, 8,0,0xFF00,   8,0xFF},
                                     {0xFF, 0,0,0xFF,     0,0xFF}};

int SYWindow::cmp(const void *p1, const void *p2)
{
//        return ((struct ZInd *)p1)->Z < ((struct ZInd *)p2)->Z ? -1 : 1;
    if (((struct SYWindow::ZInd *)p1)->Z < ((struct SYWindow::ZInd *)p2)->Z)
        return -1;
        else if (((struct SYWindow::ZInd *)p1)->Z == ((struct SYWindow::ZInd *)p2)->Z)
            return 0;
    else
        return 1;
}

#ifndef NO_THREAD
#ifdef TMCC
extern "C"
#endif
void * SYWindow::thread_routine(void *p)
{
    ((SYWindow *)p)->Routine();
    return NULL;
}
void SYWindow::Run() //initiate a new thread to run Routine()
{
    int r;
    pthread_t thread;
    r=pthread_create(&thread, NULL, thread_routine, this);
#ifdef _GENERAL_H
//        if(r!=0) FATAL("Fail to create thread: ["<<sys_errlist[errno]<<"]");
        if(r!=0) FATAL("Fail to create thread: ["<<strerror(errno)<<"]");
#else
    if(r!=0) fprintf(stderr,"Fail to create thread\n");
#endif
}
#endif

void SYWindow::Draw3DLine(YLine line)
{
    double x1,y1,x2,y2; unsigned attr;
    double xs1, ys1, zs1, xs2, ys2, zs2, Z, r;
    unsigned long col=cblack;

    attr = line.attr;
    xs1 = line.x0;
    ys1 = line.y0;
    zs1 = line.z0;
    Z=A31*xs1+A32*ys1+A33*zs1;

#define PBCSHIFT(a,shift,origin) if(enable_pbc){a-=origin;a/=2;a+=shift;a-=rint(a);a*=2;a+=origin;}

    PBCSHIFT(xs1,pbcshift[0],0);
    PBCSHIFT(ys1,-pbcshift[1],0);
    PBCSHIFT(zs1,pbcshift[2],0);
    
    x1=CX(A11*xs1+A12*ys1+A13*zs1,Z);
    y1=CY(A21*xs1+A22*ys1+A23*zs1,Z);
    xs2 = line.x1;
    ys2 = line.y1;
    zs2 = line.z1;

    PBCSHIFT(xs2,pbcshift[0],xs1);
    PBCSHIFT(ys2,-pbcshift[1],ys1);
    PBCSHIFT(zs2,pbcshift[2],zs1);
    
    Z=A31*xs2+A32*ys2+A33*zs2;
    x2=CX(A11*xs2+A12*ys2+A13*zs2,Z);
    y2=CY(A21*xs2+A22*ys2+A23*zs2,Z);
    //need some criteria to avoid waste
    //temporary
    if((((x1<0)||(x1>(int)width))||((y1<0)||(y1>(int)height)))
       &&(((x2<0)||(x2>(int)width))||((y2<0)||(y2>(int)height))))
        return;
//            INFO("("<<x1<<","<<y1<<")   ("<<x2<<","<<y2<<")");
    if(x1!=x2 || y1!=y2)
    {
        if(col!=line.c)
            XSetForeground(theDisplay, gc, col=line.c);
        if(attr==0)
        {//interprete r as line width
            r=CR(line.r,0);
            if(r<1) r=1;
            if(r>maxlinewidth) r=maxlinewidth;
            
            XSetLineAttributes(theDisplay, gc, (int)r,
                               LineSolid, CapRound, JoinRound);
            XDrawLine(theDisplay, pixmap, gc, (int)x1,(int)y1,(int)x2,(int)y2);
        }
        else
        {
            XPoint points[5]; double dy,dx,dr,r;
            dy=y2-y1;dx=x2-x1;dr=sqrt(dy*dy+dx*dx);dy/=dr;dx/=dr;
            r=CR(line.r,0);
            if(r<1) r=1;
            if(r>maxlinewidth) r=maxlinewidth;

            points[0].x=(int)(x1-dy*r/2);points[0].y=(int)(y1+dx*r/2);
            points[1].x=(int)(x1+dy*r/2);points[1].y=(int)(y1-dx*r/2);
            points[2].x=(int)(x2+dy*r/2);points[2].y=(int)(y2-dx*r/2);
            points[3].x=(int)(x2-dy*r/2);points[3].y=(int)(y2+dx*r/2);
            points[4].x=points[0].x;points[4]=points[0];
//            XSetLineAttributes(theDisplay, gc, 1,
//                               LineSolid, CapRound, JoinRound);
//            XFillPolygon(theDisplay, pixmap, gc, points, 4,
//                         Convex, CoordModeOrigin);
            XSetLineAttributes(theDisplay, gc, (int)r,
                               LineSolid, CapRound, JoinRound);
            XDrawLine(theDisplay, pixmap, gc, (int)x1,(int)y1,(int)x2,(int)y2);
            XSetLineAttributes(theDisplay, gc, 1,
                               LineSolid, CapRound, JoinRound);
            XSetForeground(theDisplay, gc, cblack);
            XDrawLines(theDisplay, pixmap, gc, points+1, 2,
                         CoordModeOrigin);
            XDrawLines(theDisplay, pixmap, gc, points+3, 2,
                         CoordModeOrigin);
        }
    }
}

void SYWindow::Draw3DPixel(YPoint point)
{
    int x,y; double r; unsigned attr;
    double xs, ys, zs, Z;
    unsigned long col=cblack;

    xs = point.x;
    ys = point.y;
    zs = point.z;

    PBCSHIFT(xs,pbcshift[0],0);
    PBCSHIFT(ys,-pbcshift[1],0);
    PBCSHIFT(zs,pbcshift[2],0);
    
    attr = point.attr;
    Z=A31*xs+A32*ys+A33*zs;
    x=(int)CX(A11*xs+A12*ys+A13*zs,Z);
    y=(int)CY(A21*xs+A22*ys+A23*zs,Z);

    r=CR(point.r,Z);

/*
 *  The limitRadius value will typically be set only when the point being
 *  plotted is a nodal point, in whcih case the maxpointradius is applied.
 *  If the point represents an eshelby inclusion, the radius is the 
 *  radius of the inclusion and not limited... but we do want to scale
 *  the radius properly in case the display has been zoomed in or out.
 */
    if (point.limitRadius) {
        if(r>maxpointradius) r=maxpointradius;
    } else {
        r = r * Scale;
    }
    
    if(((x<-r)||(x>(int)width+r))||((y<-r)||(y>(int)height+r))) return;
    XSetLineAttributes(theDisplay, gc, 1,
                       LineSolid, CapRound, JoinRound);
    if(col!=point.c)
        XSetForeground(theDisplay, gc, col=point.c);
//        if(r<-0.5) // -1 means great ball :)
//            draw3DCircle(x, y, col);
    if(r<0.5) XDrawPoint(theDisplay, pixmap, gc, x, y);
    else if(r<=1)
    {
        XDrawPoint(theDisplay, pixmap, gc, x+1, y);
        XDrawPoint(theDisplay, pixmap, gc, x, y);
        XDrawPoint(theDisplay, pixmap, gc, x, y+1);
        XDrawPoint(theDisplay, pixmap, gc, x+1, y+1);
    }            
    else if(r<=1.5)
    {
        //XDrawPoint(theDisplay, pixmap, gc, x, y);
        XDrawPoint(theDisplay, pixmap, gc, x+1, y);
        XDrawPoint(theDisplay, pixmap, gc, x-1, y);
        XDrawPoint(theDisplay, pixmap, gc, x, y+1);
        XDrawPoint(theDisplay, pixmap, gc, x+1, y+1);
        XDrawPoint(theDisplay, pixmap, gc, x-1, y+1);
        XDrawPoint(theDisplay, pixmap, gc, x, y-1);
        XDrawPoint(theDisplay, pixmap, gc, x+1, y-1);
        XDrawPoint(theDisplay, pixmap, gc, x-1, y-1);
    }
    else
    {
        if((attr==0)||(attr==2))
        { /* black border with colored disc */
            XSetForeground(theDisplay, gc, point.c);
            XFillArc(theDisplay, pixmap, gc, (int)(x-r),(int)(y-r),
                     (int)(r*2),(int)(r*2),0,23040);
            XSetForeground(theDisplay, gc, cblack);
            XDrawArc(theDisplay, pixmap, gc, (int)(x-r),(int)(y-r),
                     (int)(r*2),(int)(r*2),0,23040);
        }
        else if(attr==1)
        { /* colored border with empty disc */
            XDrawArc(theDisplay, pixmap, gc, (int)(x-r),(int)(y-r),
                     (int)(r*2),(int)(r*2),0,23040);
        }
        else
        { /* colored disc with empty border */
            XSetForeground(theDisplay, gc, point.c);
            XFillArc(theDisplay, pixmap, gc, (int)(x-r),(int)(y-r),
                     (int)(r*2),(int)(r*2),0,23040);
        }
    }
}

void SYWindow::paint()
{
    int i;
    //Collect Zs
    double xs, ys, zs;
//    DUMP("YWindow::paint()");
    if (sort)
    {
        for(i=0;i<nL;i++)
        {
            xs=(Lines[i].x0+Lines[i].x1)/2;
            ys=(Lines[i].y0+Lines[i].y1)/2;
            zs=(Lines[i].z0+Lines[i].z1)/2;

            PBCSHIFT(xs,pbcshift[0],0);
            PBCSHIFT(ys,-pbcshift[1],0);
            PBCSHIFT(zs,pbcshift[2],0);
            
            Z_ind[i].Z=A31*xs+A32*ys+A33*zs;
            Z_ind[i].index=i;
        }
        for(i=0;i<nP;i++)
        {
            xs = Points[i].x;
            ys = Points[i].y;
            zs = Points[i].z;

            PBCSHIFT(xs,pbcshift[0],0);
            PBCSHIFT(ys,-pbcshift[1],0);
            PBCSHIFT(zs,pbcshift[2],0);
            
            Z_ind[i+nL].Z=A31*xs+A32*ys+A33*zs;
            Z_ind[i+nL].index=i+nL;
        }
//        printf("qsort: nL=%d nP=%d\n",nL,nP);
        qsort(Z_ind, nL+nP, sizeof(struct ZInd), &cmp);
        
        for(i=0;i<nL+nP;i++)
        {
            int n=Z_ind[i].index;
            double z=Z_ind[i].Z;
            /* cut thin film */
            if(thinfilm)
              if((z>film_zmax)||(z<film_zmin)) continue;
               
            if(n<nL) Draw3DLine(Lines[n]);
            else Draw3DPixel(Points[n-nL]);
        }
//    SetColor(0xffffff);
//    DrawCircle(width/2,height/2,R);
    }
    else
    {
        for(i=0;i<nL;i++)
        {
            Draw3DLine(Lines[i]);
        }
        for(i=0;i<nP;i++)
        {
            Draw3DPixel(Points[i]);
        }
    }        
}
#if 0
void SYWindow::paint()
{
    int i;
    double xs, ys, zs, Z;
    unsigned long col=cblack;

   if(!alive) return;
    //Draw Segments
    {
        int x1, y1, x2, y2;
        XSetLineAttributes(theDisplay, gc, 2,  LineSolid, CapButt, JoinMiter);
        for(i=0;i<nL;i++)
        {
            xs = Lines[i].x0;
            ys = Lines[i].y0;
            zs = Lines[i].z0;
            Z=A31*xs+A32*ys+A33*zs;
            x1=CX(A11*xs+A12*ys+A13*zs,Z);
            y1=CY(A21*xs+A22*ys+A23*zs,Z);
            xs = Lines[i].x1;
            ys = Lines[i].y1;
            zs = Lines[i].z1;
            Z=A31*xs+A32*ys+A33*zs;
            x2=CX(A11*xs+A12*ys+A13*zs,Z);
            y2=CY(A21*xs+A22*ys+A23*zs,Z);
            //need some criteria to avoid waste
            //temporary
            if((((x1<0)||(x1>(int)width))||((y1<0)||(y1>(int)height)))
                &&(((x2<0)||(x2>(int)width))||((y2<0)||(y2>(int)height))))
                continue;
//            INFO("("<<x1<<","<<y1<<")   ("<<x2<<","<<y2<<")");
            if(x1!=x2 || y1!=y2)
            {
                if(col!=Lines[i].c)
                    XSetForeground(theDisplay, gc, col=Lines[i].c);
                XDrawLine(theDisplay, pixmap, gc, x1,y1,x2,y2);
            }
        }
    }
    //Draw Points
    for(i=0;i<nP;i++)
    {
        int x, y, r; unsigned attr;
        xs = Points[i].x;
        ys = Points[i].y;
        zs = Points[i].z;
        attr = Points[i].attr;
        Z=A31*xs+A32*ys+A33*zs;
        x=CX(A11*xs+A12*ys+A13*zs,Z);
        y=CY(A21*xs+A22*ys+A23*zs,Z);
        r=CR(Points[i].r);
        if(((x<-r)||(x>(int)width+r))||((y<-r)||(y>(int)height+r))) continue;
        if(col!=Points[i].c)
            XSetForeground(theDisplay, gc, col=Points[i].c);
//        if(r<-0.5) // -1 means great ball :)
//            draw3DCircle(x, y, col);
        if(r<=1) XDrawPoint(theDisplay, pixmap, gc, x, y);
        else if(r<=2)
        {
            XDrawPoint(theDisplay, pixmap, gc, x+1, y);
            XDrawPoint(theDisplay, pixmap, gc, x, y);
            XDrawPoint(theDisplay, pixmap, gc, x, y+1);
            XDrawPoint(theDisplay, pixmap, gc, x+1, y+1);
        }            
        else if(r<=3)
        {
            //XDrawPoint(theDisplay, pixmap, gc, x, y);
            XDrawPoint(theDisplay, pixmap, gc, x+1, y);
            XDrawPoint(theDisplay, pixmap, gc, x-1, y);
            XDrawPoint(theDisplay, pixmap, gc, x, y+1);
            XDrawPoint(theDisplay, pixmap, gc, x+1, y+1);
            XDrawPoint(theDisplay, pixmap, gc, x-1, y+1);
            XDrawPoint(theDisplay, pixmap, gc, x, y-1);
            XDrawPoint(theDisplay, pixmap, gc, x+1, y-1);
            XDrawPoint(theDisplay, pixmap, gc, x-1, y-1);
        }
        else
            if(attr==0)
                XFillArc(theDisplay, pixmap, gc, x-r,y-r,r*2,r*2,0,23040);
            else
                XDrawArc(theDisplay, pixmap, gc, x-r,y-r,r*2,r*2,0,23040);
    }
}
#endif

int SYWindow::identify(int px, int py)
{
    int i;
    int hit;
    double xs, ys, zs, Z;
    hit=-1;
    for(i=0;i<nP;i++)
    {
        int x, y, r;
        if(Points[i].dscrpt[0]!=0)
        {
            xs = Points[i].x;
            ys = Points[i].y;
            zs = Points[i].z;
            
            PBCSHIFT(xs,pbcshift[0],0);
            PBCSHIFT(ys,-pbcshift[1],0);
            PBCSHIFT(zs,pbcshift[2],0);
    
            Z=A31*xs+A32*ys+A33*zs;
            x=(int)CX(A11*xs+A12*ys+A13*zs,Z);
            y=(int)CY(A21*xs+A22*ys+A23*zs,Z);
            r=(int)CR(Points[i].r,Z);
            if((px<x+r)&&(px>x-r)&&(py<y+r)&&(py>y-r))
            {
                hit=i;
#ifdef _GENERAL_H
                INFO("point("<<px<<","<<py<<"): "<<Points[i].dscrpt);
#else
                fprintf(stdout,"point (%d,%d): %s\n",px,py,Points[i].dscrpt);
#endif
            }
        }
    }
    return hit;
}

void SYWindow::update()
{
    if(!alive) return;
    Lock();
//    XSetForeground(theDisplay, gc, cblack);
    XSetForeground(theDisplay, gc, bgcolor);
    
    XFillRectangle(theDisplay, pixmap, gc, 0, 0,
                   width,height);
    paint();
//    drawBoxFrame();
    XCopyArea(theDisplay,pixmap,theWindow,gc,0,0,width,height,0,0);
    Unlock();
}


void SYWindow::newGraph()
{
    unsigned int Rg, w, h;
    unsigned int border_width;
    int xself, yself;

    XGetGeometry(theDisplay, theWindow,
                          &Root,&xself,&yself,&w,&h,
                          &border_width,&depth);
    
    Rg=w< h?w:h;
    if(!alive)
    {
        alive=true;
        goto init;
    }
    if(width!=w || height!=h)
    {
        XFreeGC(theDisplay, gc);
        XFreePixmap(theDisplay, pixmap);
    init:
        R=(Rg+1)/2-5;
#ifndef M_SQRT3 //sqrt(3)
#define M_SQRT3 1.73205080756887729353
#endif
        B=(int)(R/M_SQRT3);
        double rx0=((double)X0)/width,ry0=((double)Y0)/height;
        X0=(int)(rx0*w);
        Y0=(int)(ry0*h);
        width=w; height=h;
        pixmap = XCreatePixmap(theDisplay, Root, width,height,depth);
        gc = XCreateGC(theDisplay, pixmap, 0, 0);
    }
    dirty=true;
}

void SYWindow::Evolve()
{
    XEvent ev; KeySym *ks;
    int keysyms_per_keycode_return;
    if(!alive) return;
    XSync(theDisplay, false);
    if(dirty)
    {
        dirty=false;
        update();
    }
#ifndef NO_XPM
    if(autowritegif)
    {
        writegif();
        UnlockWritegif();
    }
#endif
    while(XPending(theDisplay))
    {
        XNextEvent(theDisplay, &ev);
        switch(ev.type)
        {
        case ButtonPress:
            Xp=ev.xbutton.x;
            Yp=ev.xbutton.y;
            identify(Xp,Yp);
            break;
        case Expose:
            newGraph();
            break;
        case KeyPress:
            ks=XGetKeyboardMapping(theDisplay,ev.xkey.keycode,1,
                                   &keysyms_per_keycode_return);
            switch (*ks)
            {
            case XK_Escape:
                FreeResource(); return;
#ifdef _GENERAL_H
            case XK_p: TogglePause();INFO("Pause="<<pause);return;
#else
            case XK_p: TogglePause();
                fprintf(stdout,"Pause=%d",(int)pause);return;
#endif
#ifndef NO_XPM
//            case XK_F9: writegif(); return;
#endif
            case XK_F9: importgif(); return;
            case XK_F10: writeps(); return;
            case XK_Up:
                switch(msrsp){
                case(SCALING):Scale/=1.1;dirty=true;break;
                case(TRANSLATION):Y0-=5;dirty=true;break;
                case(ASPECTRATIO):Aspr*=1.1;dirty=true;break;
                } break;
            case XK_Down:
                switch(msrsp){
                case(SCALING):Scale*=1.1;dirty=true;break;
                case(TRANSLATION):Y0+=5;dirty=true;break;
                case(ASPECTRATIO):Aspr/=1.1;dirty=true;break;
                } break;
            case XK_Left:
                switch(msrsp){
                case(SCALING):Scale/=1.1;dirty=true;break;
                case(TRANSLATION):X0-=5;dirty=true;break;
                case(ASPECTRATIO):Scale/=1.1;Aspr*=1.1;dirty=true;break;
                } break;
            case XK_Right:
                switch(msrsp){
                case(SCALING):Scale*=1.1;dirty=true;break;
                case(TRANSLATION):X0+=5;dirty=true;break;
                case(ASPECTRATIO):Scale*=1.1;Aspr/=1.1;dirty=true;break;
                } break;
            case XK_Home:
                X0=X00;Y0=Y00;Scale=Scale0;Aspr=Aspr0;
                dirty=true;break;
            case XK_a:msrsp=ASPECTRATIO;dirty=true;break;
            case XK_s:msrsp=SCALING;dirty=true;break;
            case XK_t:msrsp=TRANSLATION;dirty=true;break;
            default: ExtKeyHandler(*ks); break;
            }
        default: break;
        }
    }
}

void SYWindow::DrawPoint(double x, double y, double z, double r,
                        unsigned long c, unsigned long attr)
{
    if(nP<MaxPoints-1)
    {
        Points[nP].x=x;
        Points[nP].y=y;
        Points[nP].z=z;
        Points[nP].r=r;
        Points[nP].c=c;
        Points[nP].attr=attr;
        Points[nP].dscrpt[0]=0;
        nP++;
    }
    else if(nP==MaxPoints-1)
    {
#ifdef _GENERAL_H
        WARNING("SYWindow: Too many points (more than "
                 << (int)MaxPoints
                 << "), ignoring");
#else
        fprintf(stderr,"SYWindow: Too many points (more than %d"
                 "), ignoring\n\n", (int)MaxPoints );
#endif
        nP++;
    }
}

void SYWindow::DrawPoint(double x, double y, double z, double r,
                        int radiusInPixels, double lMax, unsigned long c,
                        char *ds, unsigned long attr)

{
    if(nP<MaxPoints-1)
    {
        Points[nP].x=x;
        Points[nP].y=y;
        Points[nP].z=z;
/*
 *      If radius is in pixels, we're dealing with a normal nodal point,
 *      so we want to limit the size of the point.  Otherwise, it's probably
 *      an eshelby inclusion, so we do not want to limit it's size, and
 *      we need to convert radius in b to the radius in pixels.
 */
        if (radiusInPixels == 1) {
            Points[nP].r=r;
            Points[nP].limitRadius=1;
        } else {
            Points[nP].r = r * width / lMax;
            Points[nP].limitRadius=0;
        }

        Points[nP].c=c;
        Points[nP].attr=attr;
        strncpy(Points[nP].dscrpt,ds,DSCLEN-1);
        nP++;
    }
    else if(nP==MaxPoints-1)
    {
#ifdef _GENERAL_H
        WARNING("SYWindow: Too many points (more than "
                 << (int)MaxPoints
                 << "), ignoring");
#else
        fprintf(stderr,"SYWindow: Too many points (more than %d"
                 "), ignoring\n\n", (int)MaxPoints );
#endif
        nP++;
    }
}

void SYWindow::DrawLine(double x0, double y0, double z0,
                        double x1, double y1, double z1, unsigned long c,
                        double r, unsigned long attr)
{
    if(nL<MaxLines-1)
    {
        Lines[nL].x0=x0;
        Lines[nL].y0=y0;
        Lines[nL].z0=z0;
        Lines[nL].x1=x1;
        Lines[nL].y1=y1;
        Lines[nL].z1=z1;
        Lines[nL].c=c;
        Lines[nL].r=r;
        Lines[nL].attr=attr;
        nL++;
    }
    else if(nL==MaxLines-1)
    {
#ifdef _GENERAL_H
        WARNING("SYWindow: Too many lines (more than "
                 << (int)MaxLines
                 << "), ignoring");
#else
        fprintf(stderr,"SYWindow: Too many lines (more than %d"
                 "), ignoring\n\n", (int)MaxLines );
#endif
        nL++;
    }
}

SYWindow::SYWindow(int width_hint, int height_hint,
                   const char *winname, int s, int so, int fm)
{
    XSizeHints myhint;
    unsigned long mask; 
    Cursor arrow;
    int i,j;
    
    //fprintf(stderr,"SYWindow initializer");fflush(stderr);
    initRot();copyRot();
    theDisplay = XOpenDisplay(getenv("DISPLAY"));
    theScreen = DefaultScreen(theDisplay);

    cblack = BlackPixel(theDisplay,theScreen);
    cwhite = WhitePixel(theDisplay,theScreen);
    bgcolor=cblack;
    myhint.x = 100;
    myhint.y = 100;
    myhint.width = width_hint;
    myhint.height = height_hint;
    myhint.flags = PPosition | PSize;

    width=width_hint;
    height=height_hint;

    Root = RootWindow(theDisplay,theScreen);
    cmap = DefaultColormap(theDisplay, theScreen);
    /*
    Root = DefaultRootWindow(theDisplay);
    cmap = DefaultColormap(theDisplay, DefaultScreen(theDisplay));
    */
    PixDepth = 0;
    
    if( XMatchVisualInfo(theDisplay,theScreen,32,TrueColor,&visinfo) ||
        XMatchVisualInfo(theDisplay,theScreen,32,DirectColor,&visinfo) )
    {
        PixDepth = 32;
    }else if( XMatchVisualInfo(theDisplay,theScreen,24,TrueColor,&visinfo) ||
               XMatchVisualInfo(theDisplay,theScreen,24,DirectColor,&visinfo) )
    {   PixDepth = 24;
    }
    else if( XMatchVisualInfo(theDisplay,theScreen,16,TrueColor,&visinfo) ||
             XMatchVisualInfo(theDisplay,theScreen,32,DirectColor,&visinfo) )
    {   PixDepth = 16;
    }
    else if( XMatchVisualInfo(theDisplay,theScreen,8,PseudoColor,&visinfo) ||
             XMatchVisualInfo(theDisplay,theScreen,8,GrayScale,&visinfo) )
    {   PixDepth = 8;
    }

    if(PixDepth==16)
    {
        for(i=0;i<3;i++)
            for(j=0;j<6;j++)
                CCT[i][j]=colorconvert_16[i][j];
    }
    else
    {
        for(i=0;i<3;i++)
            for(j=0;j<6;j++)
                CCT[i][j]=colorconvert_32[i][j];
    }
    
    if((PixDepth!=0)&1/* 1/0 for debug */)
    {
        vis=visinfo.visual;
        cblack=AllocShortRGBColor(0,0,0);
        cwhite=AllocShortRGBColor(255,255,255);
        if( vis != DefaultVisual(theDisplay,theScreen) )
        {   cmap = XCreateColormap(theDisplay,Root,vis,AllocNone);
        }
        arrow = XCreateFontCursor(theDisplay,XC_top_left_arrow);
        
        mask = CWEventMask;
        attr.event_mask = ExposureMask | KeyPressMask | StructureNotifyMask
            | EnterWindowMask | LeaveWindowMask | PropertyChangeMask;
        attr.background_pixel = cblack;     mask |= CWBackPixel;
        attr.border_pixel = cwhite;         mask |= CWBorderPixel;
        attr.colormap = cmap;               mask |= CWColormap;
        attr.cursor = arrow;                mask |= CWCursor;
        theWindow = XCreateWindow(theDisplay,
                                  Root,0,0,
                                  myhint.width, myhint.height,
                                  2, PixDepth, InputOutput,
                                  vis, mask, &attr);
    }
    else
    {
#ifdef _GENERAL_H
        INFO("XCreateSimpleWindow");
#else
        printf("XCreateSimpleWindow");
#endif
        theWindow = XCreateSimpleWindow(theDisplay,
                                        Root,
                                        myhint.x, myhint.y,
                                        myhint.width, myhint.height,
                                        5, cwhite, cblack);
    }

    XSetStandardProperties(theDisplay, theWindow, winname, winname,
                           None, NULL, 0, &myhint);

    XSelectInput(theDisplay, theWindow,
                 ButtonPressMask|ButtonReleaseMask|Button1MotionMask|
                 KeyPressMask|ExposureMask);

    XMapRaised(theDisplay, theWindow);

    alive=false;
    pause=false;

    autowritegif=false;
    
    framecolor = cwhite;
    newGraph();


    nL=nP=0;
    
    //For lock/unlock
//    semID=semget(IPC_PRIVATE, 1, IPC_CREAT|0777);
//    printf("Init semID=%d\n",semID);
//    if(semID==-1) printf("semget failure!\n");
    InitSem();
    Unlock();
    //For lock/unlock
//    semID2=semget(IPC_PRIVATE, 1, IPC_CREAT|0777);
//    printf("Init semID2=%d\n",semID2);
//    if(semID2==-1) printf("semget failure!\n");
    InitSem2();
//    UnlockWritegif();
    //set bool square
    square = s;
    sort = so;
    drawframe = fm;
    X0=Y0=X00=Y00=0;Scale=Scale0=1;Aspr=Aspr0=1;

    thinfilm=false;
    film_zmin=-0.1;
    film_zmax= 0.1;
    //projection
    D=D0=10000;

    gifcount=0;pscount=0;

//    lastDrag = 0xFFFFFFFF;
    lastDrag = (unsigned long) (-1);

    rinterval=0;
    msrsp=ROTATION;

    scalepoints = 1;
    enable_pbc = 0;
    pbcshift[0]=pbcshift[1]=pbcshift[2]=0 ;
    maxpointradius = 100;
    maxlinewidth = 100;
}

SYWindow::~SYWindow()
{
//    if(alive)
//    {
//        //alive=false;
//        XFreeGC(theDisplay, gc);
//        XFreePixmap(theDisplay, pixmap);
//    }        
//    XFreeColormap(theDisplay, cmap);
//    XDestroyWindow(theDisplay, theWindow);
//    XCloseDisplay(theDisplay);

    FreeResource();
}

void SYWindow::reversergb()
{
    int j;
    unsigned long L;
    
    for(j=0;j<6;j++)
    {
        L=CCT[0][j];
        CCT[0][j]=CCT[2][j];
        CCT[2][j]=L;
    }
}

unsigned long SYWindow::AllocRGBColor(unsigned r, unsigned g, unsigned b)
{
    XColor c;
    c.red=r, c.green=g, c.blue=b;
//    unsigned long p;double red,green,blue;
#ifdef _GENERAL_H
    if(XAllocColor(theDisplay, cmap, &c)==0)
        WARNING("Error allocating color ("<<r<<", "<<g<<", "<<b<<")");
#else
    if(XAllocColor(theDisplay, cmap, &c)==0)
        fprintf(stderr,"Error allocating color (%d,%d,%d)\n",r,g,b);
#endif
//    p=c.pixel;
//    r=c.red;g=c.green;b=c.blue;
//    red=((p&0x7C00)>>10)*1.0/(0x001F);
//    green=((p&0x03E0)>>5)*1.0/(0x001F);
//    blue=(p&0x001F)*1.0/(0x001F);
//    INFO_Printf("c.pixel=%x c.red=%x(%f) c.green=%x(%f) c.blue=%x(%f)\n",
//                p,r,red,g,green,b,blue);
    return c.pixel;
}

unsigned long SYWindow::AllocShortRGBColor(unsigned r, unsigned g, unsigned b)
{
//    return RGB(r,g,b);
    return RGBany(r,g,b);
}

unsigned long SYWindow::AllocNamedColor(char *name)
{
    XColor c,c1;
//    unsigned long p;double red,green,blue;
//    unsigned short r,g,b;
#ifdef _GENERAL_H
    if(XAllocNamedColor(theDisplay, cmap, name,&c1,&c)==0)
        WARNING("Error allocating color ("<<name<<")");
#else
    if(XAllocNamedColor(theDisplay, cmap, name,&c1,&c)==0)
        fprintf(stderr,"Error allocating color ( %s )\n",name);
#endif
//    p=c.pixel;r=c.red;g=c.green;b=c.blue;
//    red=((p&0x7C00)>>10)*1.0/(0x001F);
//    green=((p&0x03E0)>>5)*1.0/(0x001F);
//    blue=(p&0x001F)*1.0/(0x001F);
//    INFO_Printf("c.pixel=%x c.red=%x(%f) c.green=%x(%f) c.blue=%x(%f)\n",
//                p,r,red,g,green,b,blue);
//    INFO_Printf("c.pixel=%x c.red=%x c.green=%x c.blue=%x\n",
//                c.pixel,c.red,c.green,c.blue);
    return c.pixel;
}

void SYWindow::testcolor()
{
    unsigned long c;
    int r, g, b;
    unsigned int ur,ug,ub;
#ifndef INFO_Printf
#define INFO_Printf printf
#endif
    r=255; g=0; b=0; ur=0xffff; ug=0; ub=0;
    c=AllocRGBColor(ur,ug,ub); INFO_Printf(" (%x %x %x)-%x\n",ur,ug,ub,(int)c);
    c=AllocShortRGBColor(r,g,b); INFO_Printf(" (%d %d %d)-%x\n",r,g,b,(int)c);
    r=0; g=255; b=0; ur=0; ug=0xffff; ub=0;
    c=AllocRGBColor(ur,ug,ub); INFO_Printf(" (%x %x %x)-%x\n",ur,ug,ub,(int)c);
    c=AllocShortRGBColor(r,g,b); INFO_Printf(" (%d %d %d)-%x\n",r,g,b,(int)c);
    r=0; g=0; b=255; ur=0; ug=0; ub=0xffff;
    c=AllocRGBColor(ur,ug,ub); INFO_Printf(" (%x %x %x)-%x\n",ur,ug,ub,(int)c);
    c=AllocShortRGBColor(r,g,b); INFO_Printf(" (%d %d %d)-%x\n",r,g,b,(int)c);
    r=1; g=0; b=0; ur=0x1000; ug=0; ub=0;
    c=AllocRGBColor(ur,ug,ub); INFO_Printf(" (%x %x %x)-%x\n",ur,ug,ub,(int)c);
    c=AllocShortRGBColor(r,g,b); INFO_Printf(" (%d %d %d)-%x\n",r,g,b,(int)c);
    r=0; g=1; b=0; ur=0; ug=0x1000; ub=0;
    c=AllocRGBColor(ur,ug,ub); INFO_Printf(" (%x %x %x)-%x\n",ur,ug,ub,(int)c);
    c=AllocShortRGBColor(r,g,b); INFO_Printf(" (%d %d %d)-%x\n",r,g,b,(int)c);
    r=0; g=0; b=1; ur=0; ug=0; ub=0x1000;
    c=AllocRGBColor(ur,ug,ub); INFO_Printf(" (%x %x %x)-%x\n",ur,ug,ub,(int)c);
    c=AllocShortRGBColor(r,g,b); INFO_Printf(" (%d %d %d)-%x\n",r,g,b,(int)c);
    r=127; g=0; b=0; ur=0x7fff; ug=0; ub=0;
    c=AllocRGBColor(ur,ug,ub); INFO_Printf(" (%x %x %x)-%x\n",ur,ug,ub,(int)c);
    c=AllocShortRGBColor(r,g,b); INFO_Printf(" (%d %d %d)-%x\n",r,g,b,(int)c);
    r=0; g=127; b=0; ur=0; ug=0x7fff; ub=0;
    c=AllocRGBColor(ur,ug,ub); INFO_Printf(" (%x %x %x)-%x\n",ur,ug,ub,(int)c);
    c=AllocShortRGBColor(r,g,b); INFO_Printf(" (%d %d %d)-%x\n",r,g,b,(int)c);
    r=0; g=0; b=127; ur=0; ug=0; ub=0x7fff;
    c=AllocRGBColor(ur,ug,ub); INFO_Printf(" (%x %x %x)-%x\n",ur,ug,ub,(int)c);
    c=AllocShortRGBColor(r,g,b); INFO_Printf(" (%d %d %d)-%x\n",r,g,b,(int)c);
    r=255; g=255; b=255; ur=0xffff; ug=0xffff; ub=0xffff;
    c=AllocRGBColor(ur,ug,ub); INFO_Printf(" (%x %x %x)-%x\n",ur,ug,ub,(int)c);
    c=AllocShortRGBColor(r,g,b); INFO_Printf(" (%d %d %d)-%x\n",r,g,b,(int)c);
 }

#ifndef NO_XPM
void SYWindow::writegif()
{
    const char fname[]="Ysnap";
    char extname[100],tmp[100];
    char f1[100], c1[100];
    Lock();
    sprintf(tmp,"%04d",gifcount);
    strcpy(extname,fname);
    strcat(extname,tmp);
    sprintf(f1,"%s.xpm",extname); //filename 1
#ifdef _GENERAL_H
    INFO("writegif -> "<<extname<<".gif");
#else
    fprintf(stdout,"writegif -> %s.gif\n",extname);
#endif
    sprintf(c1,"xpmtoppm %s.xpm | ppmtogif > %s.gif",extname,extname);
    writeXpm(f1);
    system(c1);
    sprintf(c1,"rm -f %s",f1);//command 2
    system(c1);
    gifcount++;
    Unlock();
}
#endif
void SYWindow::importgif()
{
    const char fname[]="Ysnap";
    char extname[100],tmp[100];
    char f1[100], c1[100];
    Lock();
    sprintf(tmp,"%04d",gifcount);
    strcpy(extname,fname);
    strcat(extname,tmp);
    sprintf(f1,"%s.xpm",extname); //filename 1
#ifdef _GENERAL_H
    INFO_Printf("import -window 0x%x %s.gif &\n",theWindow,extname);
#else
    printf("import -window 0x%x %s.gif &\n",(int)theWindow,extname);
#endif
    sprintf(c1,"import -window 0x%x %s.gif &",(unsigned)theWindow,extname);
    system(c1);
    gifcount++;
    Unlock();
}
void SYWindow::Draw3DLinetoPS(FILE *file,YLine line)
{
    double xs1,ys1,zs1,xs2,ys2,zs2,Z,r;unsigned attr;
    double red,green,blue;
    double x1, y1, x2, y2;
    attr = line.attr;
    xs1 = line.x0;
    ys1 = line.y0;
    zs1 = line.z0;

    PBCSHIFT(xs1,pbcshift[0],0);
    PBCSHIFT(ys1,-pbcshift[1],0);
    PBCSHIFT(zs1,pbcshift[2],0);
    
    Z=A31*xs1+A32*ys1+A33*zs1;
    x1=CX(A11*xs1+A12*ys1+A13*zs1,Z);
    y1=CY(A21*xs1+A22*ys1+A23*zs1,Z);
    xs2 = line.x1;
    ys2 = line.y1;
    zs2 = line.z1;

    PBCSHIFT(xs2,pbcshift[0],xs1);
    PBCSHIFT(ys2,-pbcshift[1],ys1);
    PBCSHIFT(zs2,pbcshift[2],zs1);
        
    Z=A31*xs2+A32*ys2+A33*zs2;
    x2=CX(A11*xs2+A12*ys2+A13*zs2,Z);
    y2=CY(A21*xs2+A22*ys2+A23*zs2,Z);
    r=CR(line.r,0);
    if(r>maxlinewidth) r=maxlinewidth;
    if(r<0.5) r=0.5;
    
    //need some criteria to avoid waste
    //temporary
    if((((x1<0)||(x1>(int)width))||((y1<0)||(y1>(int)height)))
       &&(((x2<0)||(x2>(int)width))||((y2<0)||(y2>(int)height))))
        return;
//    red=GETRED(line.c);
//    green=GETGREEN(line.c);
//    blue=GETBLUE(line.c);

    red=REDany(line.c);
    green=GREENany(line.c);
    blue=BLUEany(line.c);
    
    fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb"
            " %5.2f %5.2f %5.2f %5.2f m l s\n",
            r,red,green,blue,
            x1,height-y1,x2,height-y2);
    if(attr!=0)
    {
        double dy,dx,dr;
        dy=y2-y1;dx=x2-x1;dr=sqrt(dy*dy+dx*dx);dy/=dr;dx/=dr;
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb"
                " %5.2f %5.2f %5.2f %5.2f m l s\n",
                0.5,0.,0.,0.,
                x1+dy*r/2,height-(y1-dx*r/2),
                x2+dy*r/2,height-(y2-dx*r/2));
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb"
                " %5.2f %5.2f %5.2f %5.2f m l s\n",
                0.5,0.,0.,0.,
                x1-dy*r/2,height-(y1+dx*r/2),
                x2-dy*r/2,height-(y2+dx*r/2));
    }
}

void SYWindow::Draw3DPixeltoPS(FILE *file,YPoint point)
{
    double x,y,r; 
//#define _ALLCIRCLE
#ifndef _ALLCIRCLE
    unsigned attr;
#endif
    double xs, ys, zs, Z;
    double red,green,blue;
    xs = point.x;
    ys = point.y;
    zs = point.z;

    PBCSHIFT(xs,pbcshift[0],0);
    PBCSHIFT(ys,-pbcshift[1],0);
    PBCSHIFT(zs,pbcshift[2],0);
    
#ifndef _ALLCIRCLE
    attr = point.attr;
#endif
    Z=A31*xs+A32*ys+A33*zs;
    x=CX(A11*xs+A12*ys+A13*zs,Z);
    y=CY(A21*xs+A22*ys+A23*zs,Z);
    r=CR(point.r,Z);
    if(((x<-r)||(x>(int)width+r))||((y<-r)||(y>(int)height+r))) return;
//    red=GETRED(point.c);
//    green=GETGREEN(point.c);
//    blue=GETBLUE(point.c);
    red=REDany(point.c);
    green=GREENany(point.c);
    blue=BLUEany(point.c);

    if(r>maxpointradius) r=maxpointradius;

//#if DEPTH_REAL == 8
//    red=green=blue=1;
//#endif
//#ifdef _ALLCIRCLE
//    if (0)
//#else
//    if (attr==0)
//#endif
    if(attr==0)
    { /* black border with white disc */
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb %f %f %f ",
                0.5,
                1.,1.,1.,
                x,height-y,r);
        fprintf(file,"disc\n");
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb %f %f %f ",
                0.5,
                0.,0.,0.,x,height-y,r);
        fprintf(file,"circle\n");
    }
    else if(attr==2)
    { /* black border with colored disc */
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb %f %f %f ",
                0.5,
                red,green,blue,
                x,height-y,r);
        fprintf(file,"disc\n");
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb %f %f %f ",
                0.5,
                0.,0.,0.,x,height-y,r);
        fprintf(file,"circle\n");
    }
    else
    { /* colored circle with empty disc */
        fprintf(file,"np %5.2f slw %5.2f %5.2f %5.2f srgb %f %f %f ",
                0.5,
                red,green,blue,
                x,height-y,r);
        fprintf(file,"circle\n");
    }
}

void SYWindow::writeps()
{
    int i;
    
    const char fname[]="Yshot";
    char extname[100],tmp[100];
    char f1[100];
    FILE *file;

    char head1[500]="%!PS-Adobe-3.0  EPSF-3.0\n"
        "%%Pages: (atend)\n";
    char head2[500]="%%BoundingBox:";
    char head3[5000]=
        "%%EndComments\n"
        "/l {lineto} def\n"
        "/m {moveto} def\n"
        "/t {translate} def\n"
        "/slw {setlinewidth} def\n"
        "/srgb {setrgbcolor} def\n"
        "/np {newpath} def\n"
        "/s {stroke} def\n"
        "/disc { 0 360 arc fill } def\n"
        "/circle { 0 360 arc s } def\n";
        
    char tail[500]="showpage\n";

    sprintf(tmp,"%04d",pscount);
    strcpy(extname,fname);
    strcat(extname,tmp);
    sprintf(f1,"%s.ps",extname); //filename 1
    file=fopen(f1,"w");
    pscount++;

    fprintf(file,"%s",head1);
#ifdef _GENERAL_H
    INFO("writeps -> "<<f1);
    INFO("write to file");
#else
    fprintf(stdout,"writeps -> %s\nwrite to file\n",f1);
#endif
    fprintf(file,"%s %d %d %d %d\n",head2,0,0,width,height);
    fprintf(file,"%s",head3);

    if (sort)
    {
        for(i=0;i<nL+nP;i++)
        {
            int n=Z_ind[i].index;
            double z=Z_ind[i].Z;
            /* cut thin film */
            if(thinfilm)
              if((z>film_zmax)||(z<film_zmin)) continue;
            if(n<nL) Draw3DLinetoPS(file,Lines[n]);
            else Draw3DPixeltoPS(file,Points[n-nL]);
        }
    }
    else
    {
        for(i=0;i<nL;i++)
        {
            Draw3DLinetoPS(file,Lines[i]);
        }
        for(i=0;i<nP;i++)
        {
            Draw3DPixeltoPS(file,Points[i]);
        }
    }        

    fprintf(file,"%s",tail);
    fclose(file);
}

void YWindow::help()
{
    const char helpstr[]=
        "yview:\n"
        "Mouse drag to rotate\n"
        "Hot Keys:\n"
        "F1        : display this message\n"
        "Up        : rotate up\n"
        "Down      : rotate down\n"
        "Left      : rotate left\n"
        "Right     : rotate right\n"
        "PgUp      : rotate counterclockwise\n"
        "PgDn      : rotate clockwise\n"
        "Home      : back to initial viewpoint\n"
        "Space     : stop rotate\n"
//        "q  w      : x++ x--\n"
//        "a  s      : y++ y--\n"
//        "z  x      : z++ z--\n"
//        "F2        : circle larger\n"
//        "F3        : circle smaller\n"
//        "F5        : zoom in\n"
//        "F6        : zoom out\n"
//        "F7        : zoom reset\n"
        "p         : toggle pause\n"
        "t         : translation\n"
        "s         : scaling\n"
        "d         : move projection infinity point\n"
        "r         : rotation\n"
        "f         : toggle pbc enableness\n"
        "m         : toggle drawframe\n"
        "g         : pbc glide\n"
        "x         : pbc shift in x\n"
        "y         : pbc shift in y\n"
        "z         : pbc shift in z\n"
        "w         : print window specification\n"
        "F9        : output gif\n"
        "F10       : output postscript\n";
#ifdef _GENERAL_H
    INFO("YWindow: "<<helpstr);
#else
    fprintf(stdout,"YWindow: %s\n",helpstr);
#endif
}

void YWindow::printWinSpec()
{
    char tmp[1000];
    sprintf(tmp,"YWindow: Window Specification\n"
            "width=%d\theight=%d\n"
            "X0=%d\tY0=%d\n"
            "Scale=%f\tD=%f\n"
            "rotation matrix=[%10f %10f %10f\n"
            "                 %10f %10f %10f\n"
            "                 %10f %10f %10f]\n\n"
            ,width,height,X0,Y0,Scale,D
            ,A11,A12,A13,A21,A22,A23,A31,A32,A33
            );
#ifdef _GENERAL_H
    INFO(tmp);
#else
    fprintf(stdout,"%s\n",tmp);
#endif
}

void YWindow::setWinSpec(int x0,int y0,double s,double d,double a[3][3])
{
    X0=x0;Y0=y0;Scale=s;D=d;
    A11=a[0][0];A12=a[0][1];A13=a[0][2];
    A21=a[1][0];A22=a[1][1];A23=a[1][2];
    A31=a[2][0];A32=a[2][1];A33=a[2][2];
}

void YWindow::update()
{
    Lock();
    XSetForeground(theDisplay, gc, bgcolor);
    XFillRectangle(theDisplay, pixmap, gc, 0, 0,
                   width,height);
    paint();
    if(drawframe)drawBoxFrame();
    XCopyArea(theDisplay,pixmap,theWindow,gc,0,0,width,height,0,0);
    Unlock();
}

void YWindow::drawBoxFrame()
{
    unsigned char FLAG;
#define _DoDraw(_1,_2,_3, _4,_5,_6) \
    _XDrawLine(_1 A11 _2 A12 _3 A13,\
               _1 A21 _2 A22 _3 A23,\
               _1 A31 _2 A32 _3 A33,\
               _4 A11 _5 A12 _6 A13,\
               _4 A21 _5 A22 _6 A23,\
               _4 A31 _5 A32 _6 A33)
                                              
#define _XDrawLine(a,b,Z0,c,d,Z1) XDrawLine(theDisplay,pixmap,gc,\
                   (int)CX(a,Z0),(int)CY(b,Z0),(int)CX(c,Z1),(int)CY(d,Z1))

#define _XDash XSetLineAttributes(theDisplay, gc, 1, LineOnOffDash, CapButt, JoinMiter)
#define _XSold XSetLineAttributes(theDisplay, gc, 1, LineSolid, CapButt, JoinMiter)

#define _XSetLine(a,b) {if((a)&& !(b)) _XDash; else if(!(a)&&(b)) _XSold;}

    XSetForeground(theDisplay, gc, framecolor);
    _XDash;

    FLAG = 1<< ( ((A31<0)?4:0)+((A32<0)?2:0)+((A33<0)?1:0));

    _XSetLine(FLAG&0xA0,!(FLAG&0xA0));
    _DoDraw(+,-,+, +,+,+);
//    _XDrawLine(A11-A12+A13,A21-A22+A23,A31-A32+A33,
//              A11+A12+A13,A21+A22+A23,A31+A32+A33);
    _XSetLine(FLAG&0x22,FLAG&0xA0);
    _DoDraw(-,-,+, +,-,+);
//    _XDrawLine(-A11-A12+A13,-A21-A22+A23,-A31-A
//              A11-A12+A13,A21-A22+A23);
    _XSetLine(FLAG&0x0A,FLAG&0x22);
    _DoDraw(-,+,+, -,-,+);
//    _XDrawLine(-A11+A12+A13,-A21+A22+A23,
//              -A11-A12+A13,-A21-A22+A23);
    _XSetLine(FLAG&0x88,FLAG&0x0A);
    _DoDraw(+,+,+, -,+,+);
//    _XDrawLine(A11+A12+A13,A21+A22+A23,
//              -A11+A12+A13,-A21+A22+A23);

    _XSetLine(FLAG&0x50,FLAG&0x88);
    _DoDraw(+,-,-, +,+,-);
//    _XDrawLine(A11-A12-A13,A21-A22-A23,
//              A11+A12-A13,A21+A22-A23);
    _XSetLine(FLAG&0x11,FLAG&0x50);
    _DoDraw(-,-,-, +,-,-);
//    _XDrawLine(-A11-A12-A13,-A21-A22-A23,
//              A11-A12-A13,A21-A22-A23);
    _XSetLine(FLAG&0x05,FLAG&0x11);
    _DoDraw(-,+,-, -,-,-);
//    _XDrawLine(-A11+A12-A13,-A21+A22-A23,
//              -A11-A12-A13,-A21-A22-A23);
    _XSetLine(FLAG&0x44,FLAG&0x05);
    _DoDraw(+,+,-, -,+,-);
//    _XDrawLine(A11+A12-A13,A21+A22-A23,
//              -A11+A12-A13,-A21+A22-A23);

    _XSetLine(FLAG&0xC0,FLAG&0x44);
    _DoDraw(+,+,+, +,+,-);
//    _XDrawLine(A11+A12+A13,A21+A22+A23,
//              A11+A12-A13,A21+A22-A23);
    _XSetLine(FLAG&0x30,FLAG&0xC0);
    _DoDraw(+,-,+, +,-,-);
//    _XDrawLine(A11-A12+A13,A21-A22+A23,
//              A11-A12-A13,A21-A22-A23);
    _XSetLine(FLAG&0x03,FLAG&0x30);
    _DoDraw(-,-,+, -,-,-);
//    _XDrawLine(-A11-A12+A13,-A21-A22+A23,
//              -A11-A12-A13,-A21-A22-A23);
    _XSetLine(FLAG&0x0C,FLAG&0x03);
    _DoDraw(-,+,+, -,+,-);
//    _XDrawLine(-A11+A12+A13,-A21+A22+A23,
//              -A11+A12-A13,-A21+A22-A23);
    _XSold;
    
    XDrawArc(theDisplay, pixmap, gc, width/2+(int)(X0-R*Scale),height/2+(int)(Y0-R*Scale),
             (int)(2*R*Scale),(int)(2*R*Scale),0,23040);
}    


//Basic Rotations----------------------------------------
void YWindow::applyRotate()
{
    double d1,d2,d3;
    d1=B11*A11+B12*A21+B13*A31;
    d2=B21*A11+B22*A21+B23*A31;
    d3=B31*A11+B32*A21+B33*A31;
    A11=d1; A21=d2; A31=d3;
    d1=B11*A12+B12*A22+B13*A32;
    d2=B21*A12+B22*A22+B23*A32;
    d3=B31*A12+B32*A22+B33*A32;
    A12=d1; A22=d2; A32=d3;
    d1=B11*A13+B12*A23+B13*A33;
    d2=B21*A13+B22*A23+B23*A33;
    d3=B31*A13+B32*A23+B33*A33;
    A13=d1; A23=d2; A33=d3;
    dirty=true;
}

void YWindow::horizontalRot(double arc)
{
    double d1,d3;
    double c, s;
    c=cos(arc);
    s=sin(arc);
    d1=c*A11+s*A31;
    d3=-s*A11+c*A31;
    A11=d1;A31=d3;
    d1=c*A12+s*A32;
    d3=-s*A12+c*A32;
    A12=d1;A32=d3;
    d1=c*A13+s*A33;
    d3=-s*A13+c*A33;
    A13=d1;A33=d3;
    dirty=true;
}

void YWindow::verticalRot(double arc)
{
    double d2,d3;
    double c,s;
    s=sin(arc);
    c=cos(arc);
    d2=c*A21+s*A31;
    d3=-s*A21+c*A31;
    A21=d2;A31=d3;
    d2=c*A22+s*A32;
    d3=-s*A22+c*A32;
    A22=d2;A32=d3;
    d2=c*A23+s*A33;
    d3=-s*A23+c*A33;
    A23=d2;A33=d3;
    dirty=true;
}

void YWindow::spinRot(double arc)
{
    double d1,d2;
    double c,s;
    s=sin(arc);
    c=cos(arc);
    d1=c*A11+s*A21;
    d2=-s*A11+c*A21;
    A11=d1;A21=d2;
    d1=c*A12+s*A22;
    d2=-s*A12+c*A22;
    A12=d1;A22=d2;
    d1=c*A13+s*A23;
    d2=-s*A13+c*A23;
    A13=d1;A23=d2;
    dirty=true;
}

void YWindow::zoom(double z)
{
//    INFO("Scale="<<Scale<<"  z="<<z);
    Scale*=z;
//    INFO("  Scale="<<Scale);
}

void YWindow::scaleTo(XEvent ev)
{
    Scale*=exp((ev.xmotion.x-Xp+ev.xmotion.y-Yp)*0.001);
    if(Scale<5e-2) Scale=5e-2;
    if(Scale>1e2) Scale=1e2;
//    DUMP("YWindow: Scale="<<Scale);
    dirty=true;
    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
}

void YWindow::translateTo(XEvent ev)
{
    X0+=ev.xmotion.x-Xp;Y0+=ev.xmotion.y-Yp;
//    DUMP("YWindow: X0="<<X0<<" ,Y0="<<Y0);
    dirty=true;
    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
}

void YWindow::pbcshiftTo(XEvent ev, int dir)
{
    pbcshift[dir]+=(0.0+ev.xmotion.x-Xp+ev.xmotion.y-Yp)/(2.0*B*Scale);
//    printf("pbcshift[%d]=%e, B=%d\n",dir,pbcshift[dir],B);
    dirty=true;
    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
}

void YWindow::pbcglideTo(XEvent ev)
{
    pbcglide((ev.xmotion.x-Xp)/(2.0*B),(ev.xmotion.y-Yp)/(2.0*B));
//    printf("pbcshift[%d]=%e, B=%d\n",dir,pbcshift[dir],B);
    dirty=true;
    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
}

void YWindow::pbcglide(double dx, double dy)
{
    //double Ainv11,Ainv12,Ainv13,Ainv21,Ainv22,Ainv23,Ainv31,Ainv32,Ainv33;
    /* Ainv = inv (A) */
    double Ainv11=A22*A33-A23*A32;
    double Ainv22=A33*A11-A31*A13;
//  double Ainv33=A11*A22-A12*A21;
    double Ainv12=A23*A31-A21*A33;
//  double Ainv23=A31*A12-A32*A11;
    double Ainv31=A12*A23-A13*A22;
//  double Ainv13=A21*A32-A31*A22;
    double Ainv21=A32*A13-A12*A33;
    double Ainv32=A13*A21-A23*A11;
        
    pbcshift[0]+=Ainv11*dx+Ainv12*dy;
    pbcshift[1]+=Ainv21*dx+Ainv22*dy;
    pbcshift[2]+=Ainv31*dx+Ainv32*dy;
}

void YWindow::projectTo(XEvent ev)
{
    D*=exp(-(ev.xmotion.x-Xp+ev.xmotion.y-Yp)*0.001);
    if(D<0.2) D=0.2;
    if(D>1e2) D=1e2;
//    DUMP("YWindow: D="<<D);
    dirty=true;
    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
}

void YWindow::rotateTo(XEvent ev)
{
    
    double a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,e1,e2,e3;
    double r,rp,xx,yy,xp,yp;
#define _toler 1.0e-6
    xx=((double)ev.xmotion.x-width/2-X0)/R;
    yy=-((double)ev.xmotion.y-height/2-Y0)/R;
    xp=((double)Xp-width/2-X0)/R;
    yp=-((double)Yp-height/2-Y0)/R;
    
    if((xp-xx)*(xp-xx)+(yp-yy)*(yp-yy)<_toler) return;
//    DUMP("YWindow: Rotate from ("<<xp<<","<<yp<<") to ("<<xx<<","<<yy<<")");
    
    rp=sqrt(xp*xp+yp*yp);
    r=sqrt(xx*xx+yy*yy);
    if(r>=1 || rp>=1)
    {
        d1=xp/rp;
        e1=yp/rp;
        d2=xx/r;
        e2=yy/r;
        
        d3=cos((r-rp)*M_PI);
        e3=sin((r-rp)*M_PI);
        
        B11=d3*d1*d2+e1*e2;
        B12=d3*e1*d2-d1*e2;
        B13=e3*d2;
        B21=d3*d1*e2-e1*d2;
        B22=d3*e1*e2+d1*d2;
        B23=e3*e2;
        B31=-e3*d1;
        B32=-e3*e1;
        B33=d3;
    }
    else
    {
        a1=xp;a2=yp;a3=sqrt(1-a1*a1-a2*a2);
        b1=xx;b2=yy;b3=sqrt(1-b1*b1-b2*b2);
        c1=a2*b3-a3*b2;
        c2=a3*b1-a1*b3;
        c3=a1*b2-a2*b1;
        r=sqrt(c1*c1+c2*c2+c3*c3);
        c1/=r;c2/=r;c3/=r;
        
        d1=a2*c3-a3*c2;
        d2=a3*c1-a1*c3;
        d3=a1*c2-a2*c1;
        e1=b2*c3-b3*c2;
        e2=b3*c1-b1*c3;
        e3=b1*c2-b2*c1;
        
        B11=b1*a1+e1*d1+c1*c1;
        B12=b1*a2+e1*d2+c1*c2;
        B13=b1*a3+e1*d3+c1*c3;
        B21=b2*a1+e2*d1+c2*c1;
        B22=b2*a2+e2*d2+c2*c2;
        B23=b2*a3+e2*d3+c2*c3;
        B31=b3*a1+e3*d1+c3*c1;
        B32=b3*a2+e3*d2+c3*c2;
        B33=b3*a3+e3*d3+c3*c3;

    }

//    xp=xx;
//    yp=yy;

    Xp=ev.xmotion.x;
    Yp=ev.xmotion.y;
    
    applyRotate();
}

void YWindow::Evolve()
{
    XEvent ev;KeySym *ks;
    int keysyms_per_keycode_return;
    long t;

    if(!alive) return;
    XSync(theDisplay, false);
    if(enableRot) applyRotate();
    if(dirty)
    {
        dirty=false;
        update();
    }
#ifndef NO_XPM
    if(autowritegif)
    {
        writegif();
        UnlockWritegif();
    }
#endif
    while(XPending(theDisplay))
    {
        XNextEvent(theDisplay, &ev);
        switch(ev.type)
        {
        case ButtonPress:
            enableRot=false;
            Xp=ev.xbutton.x;
            Yp=ev.xbutton.y;
            identify(Xp,Yp);
//            xp=((double)ev.xbutton.x-width/2)/R;
//            yp=-((double)ev.xbutton.y-height/2)/R;
            break;
        case ButtonRelease:
            t=(long)ev.xbutton.time-(long)lastDrag;
            if(msrsp==ROTATION) enableRot = (t<100);
//            DUMP("YWindow: ButtonRelease enableRot="<<enableRot<<" t="<<t);
            break;
        case MotionNotify:
            if(lastDrag+rInterval*3/4 <= ev.xmotion.time)
            {
                lastDrag=ev.xmotion.time;
                switch(msrsp){
                case(ROTATION):
                    rotateTo(ev);break;
                case(TRANSLATION):
                    translateTo(ev);break;
                case(SCALING):
                    scaleTo(ev);break;
                case(PROJECTION):
                    projectTo(ev);break;
                case(PBCX):
                    pbcshiftTo(ev,0);break;
                case(PBCY):
                    pbcshiftTo(ev,1);break;
                case(PBCZ):
                    pbcshiftTo(ev,2);break;
                case(PBCGLIDE):
                    pbcglideTo(ev);break;
                default:break;
                }
            }
//            else DUMP("YWindow: skipping too fast motion t="
//                      <<ev.xmotion.time - lastDrag);
            break;
        case Expose:
            newGraph();
            break;
        case KeyPress:
            ks=XGetKeyboardMapping(theDisplay,ev.xkey.keycode,1,
                                   &keysyms_per_keycode_return);
            switch (*ks)
            {
            case XK_F1: help();break;
            case XK_Escape:
                FreeResource(); return;
            case XK_space: enableRot=false;break; // Stop rotating
            case XK_Page_Up: spinRot(DEG(-5));break;
            case XK_Page_Down: spinRot(DEG(5));break;
            case XK_Up:
                switch(msrsp){
                case(PROJECTION):D*=1.1;dirty=true;break;
                case(SCALING):Scale/=1.1;dirty=true;break;
                case(TRANSLATION):Y0-=5;dirty=true;break;
                case(ROTATION):verticalRot(DEG(5));break;
                case(ASPECTRATIO):Aspr*=1.1;dirty=true;break;
                case(PBCX):pbcshift[0]=0.01;dirty=true;break;
                case(PBCY):pbcshift[1]-=0.01;dirty=true;break;
                case(PBCZ):pbcshift[2]-=0.01;dirty=true;break;
                case(PBCGLIDE):pbcglide(0,-0.01);dirty=true;break;
                case(THINFILM):film_zmin+=0.01;film_zmax+=0.01;
                        INFO_Printf("thin film z=(%f,%f)\n",film_zmin,film_zmax);
                        dirty=true;break;
                } break;
            case XK_Down:
                switch(msrsp){
                case(PROJECTION):D/=1.1;dirty=true;break;
                case(SCALING):Scale*=1.1;dirty=true;break;
                case(TRANSLATION):Y0+=5;dirty=true;break;
                case(ROTATION):verticalRot(DEG(-5));break;
                case(ASPECTRATIO):Aspr/=1.1;dirty=true;break;
                case(PBCX):pbcshift[0]-=0.01;dirty=true;break;
                case(PBCY):pbcshift[1]+=0.01;dirty=true;break;
                case(PBCZ):pbcshift[2]+=0.01;dirty=true;break;
                case(PBCGLIDE):pbcglide(0,+0.01);dirty=true;break;
                case(THINFILM):film_zmin-=0.01;film_zmax-=0.01;
                        INFO_Printf("thin film z=(%f,%f)\n",film_zmin,film_zmax);
                        dirty=true;break;
                } break;
            case XK_Left:
                switch(msrsp){
                case(PROJECTION):D*=1.1;dirty=true;break;
                case(SCALING):Scale/=1.1;dirty=true;break;
                case(TRANSLATION):X0-=5;dirty=true;break;
                case(ROTATION):horizontalRot(DEG(-5));break;
                case(PBCX):pbcshift[0]-=0.01;dirty=true;break;
                case(PBCY):pbcshift[1]+=0.01;dirty=true;break;
                case(PBCZ):pbcshift[2]+=0.01;dirty=true;break;
                case(PBCGLIDE):pbcglide(-0.01,0);dirty=true;break;
                case(THINFILM):if(film_zmin<film_zmax)
                        {film_zmin+=0.01;film_zmax-=0.01;}
                        INFO_Printf("thin film z=(%f,%f)\n",film_zmin,film_zmax);
                        dirty=true;break;
                } break;
            case XK_Right:
                switch(msrsp){
                case(PROJECTION):D/=1.1;dirty=true;break;
                case(SCALING):Scale*=1.1;dirty=true;break;
                case(TRANSLATION):X0+=5;dirty=true;break;
                case(ROTATION):horizontalRot(DEG(5));break;
                case(PBCX):pbcshift[0]+=0.01;dirty=true;break;
                case(PBCY):pbcshift[1]-=0.01;dirty=true;break;
                case(PBCZ):pbcshift[2]-=0.01;dirty=true;break;
                case(PBCGLIDE):pbcglide(+0.01,0);dirty=true;break;
                case(THINFILM):film_zmin-=0.01;film_zmax+=0.01;
                        INFO_Printf("thin film z=(%f,%f)\n",film_zmin,film_zmax);
                        dirty=true;break;
                } break;
            case XK_Home:
                copyRot();X0=X00;Y0=Y00;D=D0;Scale=Scale0;Aspr=Aspr0;
                pbcshift[0]=pbcshift[1]=pbcshift[2]=0;
                dirty=true;break;
#ifdef _GENERAL_H
            case XK_p: TogglePause();INFO("Pause="<<pause);return;
#else
            case XK_p: TogglePause();
                fprintf(stdout,"Pause=%d\n",(int)pause);return;
#endif
            case XK_a:msrsp=ASPECTRATIO;
                 INFO_Printf("use arrow key and mouse to adjust aspect ratio\n");
                 dirty=true;break;
            case XK_d:msrsp=PROJECTION;
                 INFO_Printf("use arrow key and mouse to adjust projection\n");
                 dirty=true;break;
            case XK_s:msrsp=SCALING;
                 INFO_Printf("use arrow key and mouse to adjust scaling\n");
                 dirty=true;break;
            case XK_t:msrsp=TRANSLATION;
                 INFO_Printf("use arrow key and mouse to adjust translation\n");
                 dirty=true;break;
            case XK_r:msrsp=ROTATION;
                 INFO_Printf("use arrow key and mouse to adjust rotation\n");
                 dirty=true;break;
            case XK_f:enable_pbc=!enable_pbc;
                      if(enable_pbc)msrsp=PBCX;else msrsp=ROTATION;
                      dirty=true; break;
            case XK_g:enable_pbc=1; msrsp=PBCGLIDE; dirty=true; break;
            case XK_x:enable_pbc=1; msrsp=PBCX; dirty=true; break;
            case XK_y:enable_pbc=1; msrsp=PBCY; dirty=true; break;
            case XK_z:enable_pbc=1; msrsp=PBCZ; dirty=true; break;
            case XK_w: printWinSpec(); break;
            case XK_m: drawframe=!drawframe; dirty=true; break;
            case XK_c: thinfilm=!thinfilm;
               if(thinfilm)
               INFO_Printf("thin film z=(%f,%f) press v to adjust thickness\n",
                           film_zmin,film_zmax);
               else INFO_Printf("thin film disabled\n");
               dirty=true; break;
            case XK_v: msrsp=THINFILM;
                    INFO_Printf("use arrow key to adjust thin film\n");
                    INFO_Printf("up/down    controls film position\n");
                    INFO_Printf("left/right controls film thickness\n");
                    dirty=true; break;
//                case XK_q: resetPBC(0.1,0,0); break;
//                case XK_a: resetPBC(0,0.1,0); break;
//                case XK_z: resetPBC(0,0,0.1); break;
//                case XK_w: resetPBC(-0.1,0,0); break;
//                case XK_s: resetPBC(0,-0.1,0); break;
//                case XK_x: resetPBC(0,0,-0.1); break;
                //case XK_e: xpbc = 2; break;
                //case XK_d: ypbc = 2; break;
                //case XK_c: zpbc = 2; break;
//            case XK_F2: rscale*=1.1; break;
//            case XK_F3: rscale/=1.1; break;
//            case XK_F4: rscale=1.0; break;
//            case XK_F5: scale*=1.1; break;
//            case XK_F6: scale/=1.1; break;
//            case XK_F7: scale=1.0; break;
#ifndef NO_XPM
//            case XK_F9:
//                writegif();
//                break;
#endif
            case XK_F9:
                importgif();
                break;
            case XK_F10:
                writeps();
                break;
            default: ExtKeyHandler(*ks); break;
            }
        default: break;
        }
    }
}

//============ Test suite

#ifdef _TEST

int main (int argc, char *argv[])
{
    int i, j, k;
    int ni=6,nj=6,nk=6;
    double x,y,z,x1,y1,z1,r,br,dx,dy,dz,dr;
    unsigned c;
    char s[100];
    YWindow *win;
    
//    YWindow yw(400,400,"Test Window Display",true);

    win=new YWindow(400,400,"Test Window Display",true);
    
//    win->maxpointradius = 40;
//    win->maxlinewidth = 10;
#define yw (*win)
    
    yw.Lock();
    yw.Clear();
    ni=nj=nk=6;

//    yw.bgcolor=yw.AllocShortRGBColor(0x50,0x50,0x50);
    yw.bgcolor=yw.AllocShortRGBColor(100,100,100);
    for(i=0;i<ni;i++)for(j=0;j<nj;j++)for(k=0;k<nk;k++)
    {
        sprintf(s,"Ball:%d,%d,%d",i,j,k);
        x=-1+1./ni+2.*i/ni;y=-1+1./nj+2.*j/nj;z=-1+1./nk+2.*k/nk;r=.5/ni;br=r/4;
        c=yw.AllocShortRGBColor(i*0x33, j*0x33, k*0x33);
        yw.DrawPoint(x,y,z,r,c,s);
        if(i>0)
        {
            x1=x-2./ni;y1=y;z1=z;
            dx=x1-x;dy=y1-y;dz=z1-z;dr=sqrt(dx*dx+dy*dy+dz*dz);dx/=dr;dy/=dr;dz/=dr;
            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,1);
//            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,0);
        }
        if(j>0)
        {
            x1=x;y1=y-2./nj;z1=z;
            dx=x1-x;dy=y1-y;dz=z1-z;dr=sqrt(dx*dx+dy*dy+dz*dz);dx/=dr;dy/=dr;dz/=dr;
            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,1);
//            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,0);
        }
        if(k>0)
        {
            x1=x;y1=y;z1=z-2./nk;
            dx=x1-x;dy=y1-y;dz=z1-z;dr=sqrt(dx*dx+dy*dy+dz*dz);dx/=dr;dy/=dr;dz/=dr;
            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,1);
//            yw.DrawLine(x+dx*r,y+dy*r,z+dz*r,x1-dx*r,y1-dy*r,z1-dz*r,c,br,0);
        }
    }
    
    yw.Unlock();
    yw.Refresh();
#ifdef NO_THREAD
    yw.Routine();
#else
    yw.Run();
    i=99;
#ifndef NO_GENERAL
    _IO << "I can live for <<"<<i<<" seconds.\n" ;
#else
    printf("I can live for %d seconds.\n",i);
#endif
    for(j=0;j<i;j++)
    {
        sleep(1);
#ifndef NO_GENERAL
        _IO << i-j << ' ';
#else
        printf("%d \n",i-j);
#endif
    }
#ifndef NO_GENERAL
    _IO << "Bye.\n";
#else
    printf("Bye.\n");
#endif
#endif
    return(0);
}

#endif //_TEST

#endif
