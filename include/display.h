#pragma once

#ifndef _PDS_DISPLAY_H
#define _PDS_DISPLAY_H

#include <stdio.h>
#include <string.h>

/*
  display.h
  by Wei Cai and Dongyi Liao  caiwei@mit.edu, liaody@mit.edu
  Last Modified : Tue Apr 29 10:47:58 2003

  FUNCTION  :  X-window display package
*/

#define _DISPLAY_VERSION 1.05

/*
 *      On the MAC, explicitly undefine SEM_SEMUN_UNDEFINED since
 *      semun is defined in sys/sem.h
 */
#ifdef __APPLE__
#undef SEM_SEMUN_UNDEFINED
#endif

#ifdef NO_GENERAL
#define true 1
#define false 0
#endif

#include <stdio.h>
#include <stdlib.h>
#include <sys/select.h>
#include <sys/time.h>

#include "mpi_portability.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>

#ifndef NO_XPM
#include <X11/xpm.h>
#endif

#ifndef NO_THREAD
#include <pthread.h>
#endif

#include <sys/ipc.h>
#include <sys/sem.h>
#if !defined(SEM_SEMUN_UNDEFINED)
/* union semun is defined by including <sys/sem.h> */
#else
/* according to X/OPEN we have to define it ourselves */
union semun
{
    int val;                    /* value for SETVAL */
    struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
    unsigned short int *array;  /* array for GETALL, SETALL */
    struct seminfo *__buf;      /* buffer for IPC_INFO */
};
#endif

#ifndef NO_GENERAL
#include "general.h"
#endif
//#include "filecls.h"

#include <math.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>

#ifndef DEPTH
//#define DEPTH 16
#define DEPTH 32
#endif

#define DEPTH_REAL 8
#if DEPTH == 16
#define RGB RGB16
#define GETRED RED16
#define GETGREEN GREEN16
#define GETBLUE BLUE16
#elif DEPTH == 32
#define RGB RGB32
#define GETRED RED32
#define GETGREEN GREEN32
#define GETBLUE BLUE32
#endif

#define RGB32(R, G, B) ((unsigned)(R)<<16)+((unsigned)(G)<<8)+((unsigned)(B))
#define RGB16(R, G, B) ((unsigned)((R)&0xF8)<<8)+((unsigned)((G)&0xF8)<<3)+((unsigned)((B)>>3))

#define RGBany(R, G, B) ((((unsigned)(R)&CCT[0][0])<<CCT[0][1])>>CCT[0][2])\
                       +((((unsigned)(G)&CCT[1][0])<<CCT[1][1])>>CCT[1][2])\
                       +((((unsigned)(B)&CCT[2][0])<<CCT[2][1])>>CCT[2][2])

#define RED32(C)   (((C&0xFF0000)>>16)*1.0/(0x00FF))
#define GREEN32(C) (((C&0xFF00)>>8)*1.0/(0x00FF))
#define BLUE32(C)  (((C&0x00FF)*1.0/(0x00FF)))

#define RED16(C)   (((C&0xF800)>>11)*1.0/(0x001F))
#define GREEN16(C) (((C&0x07C0)>>6)*1.0/(0x001F))
#define BLUE16(C)  (((C&0x001F)*1.0/(0x001F)))

#define REDany(C)   (((C&CCT[0][3])>>CCT[0][4])*1.0/(CCT[0][5]))
#define GREENany(C) (((C&CCT[1][3])>>CCT[1][4])*1.0/(CCT[1][5]))
#define BLUEany(C)  (((C&CCT[2][3])>>CCT[2][4])*1.0/(CCT[2][5]))

/* smaller set */
//enum { MaxPoints=1000, MaxLines=1000, rInterval=50, DSCLEN=10};
//enum { MaxPoints=10000, MaxLines=5000, rInterval=50, DSCLEN=80};
//enum { MaxPoints=10000, MaxLines=10000, rInterval=50, DSCLEN=80};
/* larger set */
enum { MaxPoints=40000, MaxLines=80000, rInterval=50, DSCLEN=60};

inline void delay(long sec, long usec)
{
    //Define to use timeval without including <sys/time.h>
    struct {time_t tv_sec, tv_usec;} tv; 
    tv.tv_sec=sec;
    tv.tv_usec=usec;
    select(0,NULL,NULL,NULL,(timeval *)&tv);
}

extern unsigned long colorconvert_16[3][6], colorconvert_32[3][6];

class SYWindow//simple-display
{
public:
    enum{ ROTATION=0,TRANSLATION=1,SCALING=2,PROJECTION=3,
          PBCTRANSLATION=4,ASPECTRATIO=5,
          PBCX=6,PBCY=7,PBCZ=8,PBCGLIDE=9,
          THINFILM=10};
    struct YPoint
    {
        double x, y, z, r; //coordinates always lie in [-1, 1]
        int limitRadius;
        unsigned long c;
        unsigned long attr;//0:disc 1:circle
        char dscrpt[DSCLEN];
    }Points[MaxPoints];
    struct YLine
    {
        double x0, y0, z0, x1, y1, z1, r;
        unsigned long c;
        unsigned long attr;//0:disc 1:circle
    }Lines[MaxLines];
    struct ZInd
    {
        double Z;
        int index;
    }Z_ind[MaxLines+MaxPoints];
    
    static int cmp(const void *p1, const void *p2);

    unsigned long CCT[3][6];
    
    int nP, nL; //Number of points/lines
    int semID; //for Lock/Unlock mechanism
    int semID2; //for Lock/Unlock mechanism
    int alive,pause,drawframe;
    int loop;
    int rinterval;
    int msrsp;//mouse response

    Display *theDisplay;
    int theScreen;
    XVisualInfo visinfo;
    Visual *vis;
    XSetWindowAttributes attr;
    int PixDepth;
    GC gc;
    Pixmap pixmap;
    Window theWindow, Root;
    unsigned int width, height;
    unsigned int depth;

    Colormap cmap;
    unsigned long cblack, cwhite, bgcolor, framecolor;

    double A11,A12,A13,A21,A22,A23,A31,A32,A33;//Rotation matrix.
    double a11,a12,a13,a21,a22,a23,a31,a32,a33;//saved Rotation matrix.
    double B11,B12,B13,B21,B22,B23,B31,B32,B33;//Incremental rot. matrix.

    int R, B; //Current radius, box size(in pixels)
//    double xp,yp;//Last mouse position
    int Xp,Yp;
    int gifcount,pscount;
    unsigned long lastDrag;
    
    int dirty; //true if need update
    int X0,Y0,X00,Y00;//offset
    double Scale,Scale0,Aspr,Aspr0;//scale and aspect ratio
    double D,D0; //projection
    int thinfilm; double film_zmin, film_zmax;
    int square, sort;

    int autowritegif;

    int scalepoints;
    int enable_pbc;
    double pbcshift[3];
    double maxpointradius, maxlinewidth;
    
    double CXs(double x, double z)
        { return (width/2+(B*x/(D==0?1:(1-z/D)))*Scale+X0); }
    double CYs(double y, double z)
        { return (height/2-(B*y/(D==0?1:(1-z/D)))*Scale*Aspr+Y0); }
    double CXr(double x, double z)
        { return (width/2+(width/4*x/(D==0?1:(1-z/D)))*Scale+X0); }
    double CYr(double y, double z)
        { return (height/2-(height/4*y/(D==0?1:(1-z/D)))*Scale*Aspr+Y0); }
    double CX(double x, double z) {return (square)?CXs(x,z):CXr(x,z); }
    double CY(double y, double z) {return (square)?CYs(y,z):CYr(y,z); }
    double CR(double r, double z) {if(scalepoints) return CRs(r,z); else return CRf(r,z); }
    double CRs(double r, double z) { return (B*r*Scale/(D==0?1:(1-z/D))); }
    double CRf(double r, double z) { return r; }
//    int CXs(double x, double z)
//        { return (int)(width/2+(B*x/(D==0?1:(1-z/D)))*Scale+X0); }
//    int CYs(double y, double z)
//        { return (int)(height/2-(B*y/(D==0?1:(1-z/D)))*Scale*Aspr+Y0); }
//    int CXr(double x, double z)
//        { return (int)(width/2+(width/4*x/(D==0?1:(1-z/D)))*Scale+X0); }
//    int CYr(double y, double z)
//        { return (int)(height/2-(height/4*y/(D==0?1:(1-z/D)))*Scale*Aspr+Y0); }
//    int CX(double x, double z) {return (square)?CXs(x,z):CXr(x,z); }
//    int CY(double y, double z) {return (square)?CYs(y,z):CYr(y,z); }
//    int CR(double r) { return (int)(B*r*Scale); }
    double DEG(double a) { return (M_PI*a/180); }

    void initRot()
    {
        a11=a22=a33=B11=B22=B33=1;
        a12=a13=a21=a23=a31=a32=B12=B13=B21=B23=B31=B32=0;
    }
    void copyRot()
    { A11=a11;A12=a12;A13=a13;
      A21=a21;A22=a22;A23=a23;
      A31=a31;A32=a32;A33=a33; }
    void saveRot()
    { a11=A11;a12=A12;a13=A13;
      a21=A21;a22=A22;a23=A23;
      a31=A31;a32=A32;a33=A33; }
    void saveScale()
    { Scale0=Scale;}
    void saveView()
    { saveRot();Scale0=Scale;X00=X0;Y00=Y0;D0=D;Aspr0=Aspr;}
        
    SYWindow(const SYWindow &){}
    const SYWindow &operator =(const SYWindow &yw){ return yw;}
public:
    SYWindow(int width_hint, int height_hint,
             const char *winname, int s=true, int so=false,int fm=false);
    virtual ~SYWindow();
    void setinterval(int i) {if (i>0) rinterval=i; }
    int IsAlive() { return alive; }
    int TogglePause() { pause=!pause; return pause;}
    int IsPaused() { return pause; }

#ifdef _NOLOCK
    void Lock() { }
    void Unlock() { }
    void InitSem() { }
    void InitSem2() { }
    void RemoveSem() { }
    void RemoveSem2() { }
    void LockWritegif() { }
    void UnlockWritegif() { }
#else  // _NOLOCK is not defined
#ifndef _MYOWNLOCK
    void Lock()
    {
        static sembuf sbl={0, -1, 0};
        semop(semID,&sbl,1);
    }
    void Unlock()
    {
        static sembuf sbu={0, 1, 0};
        semop(semID,&sbu,1);
    }
    void InitSem()
    {
        semID=semget(IPC_PRIVATE, 1, IPC_CREAT|0777);
        if(semID==-1) printf("semget failure! use -D_MYOWNLOCK\n");
    }
    void RemoveSem()
    {
        (void)semctl(semID, 0, IPC_RMID, 0);
    }
    void InitSem2()
    {
        semID2=semget(IPC_PRIVATE, 1, IPC_CREAT|0777);
        if(semID2==-1) printf("semget failure! use -D_MYOWNLOCK\n");
    }
    void RemoveSem2()
    {
        (void)semctl(semID2, 0, IPC_RMID, 0);
    }
    void LockWritegif()
    {
        static sembuf sbl={0, -1, 0};
        autowritegif=true;
        semop(semID2,&sbl,1);
    }
    void UnlockWritegif()
    {
        static sembuf sbu={0, 1, 0};
        autowritegif=false;
        semop(semID2,&sbu,1);
    }
#else  // _MYOWNLOCK is defined
    void Lock()
    {
        semID--;
        while(semID<0) sleep(1);
    }
    void Unlock()
    {
        semID++;
    }
    void InitSem()
    {
        printf("-D_MYOWNLOCK defined, use semaphore mimic\n");
        semID=0;
    }
    void InitSem2()
    {
        printf("-D_MYOWNLOCK defined, use semaphore mimic\n");
        semID2=0;
    }
    void LockWritegif()
    {
        semID2--;
        while(semID2<0) sleep(1);
    }
    void UnlockWritegif()
    {
        semID2++;
    }
#endif  //  end check for _MYOWNLOCK

#endif  // _NOLOCK not defined

    void DrawPoint(double x, double y, double z, double r, unsigned long c,
                   unsigned long attr=0);
    void DrawPoint(double x, double y, double z, double r, int radiusInPixels,
                   double lMax, unsigned long c,
                   char *ds, unsigned long attr=0);
    void DrawLine(double x0, double y0, double z0,
                  double x1, double y1, double z1, unsigned long c,
                  double r=0.01, unsigned long attr=0);
    void Draw3DLine(YLine);
    void Draw3DPixel(YPoint);
    void Draw3DLinetoPS(FILE *,YLine);
    void Draw3DPixeltoPS(FILE *,YPoint);
    unsigned long AllocNamedColor(char *name);
    unsigned long AllocRGBColor(unsigned r, unsigned g, unsigned b);
    unsigned long AllocShortRGBColor(unsigned r, unsigned g, unsigned b);
    unsigned long StdColor(int c)
    {
        return AllocShortRGBColor(c&1?0xff:0, c&2?0xff:0, c&4?0xff:0);
    }
    
    void Clear() { nP=nL=0; }
    void Refresh() { dirty=true; }

    virtual void paint();
    virtual void update();
    virtual void newGraph();
//    virtual void newGraph(bool init);
    virtual void Evolve();
    virtual int identify(int px, int py);
    virtual void Routine()
    {
        while(alive)
        {
            if(rinterval==0) delay(0,1000*rInterval);
            else delay(0,1000*rinterval);
            Evolve();
        }
    }
    virtual int ExtKeyHandler(KeySym ks){if((int)ks)return 0;else return 1;};
    virtual void FreeResource()
    {
        union semun su={0};
        if(alive)
        {
            XFreeGC(theDisplay, gc);
            XFreePixmap(theDisplay, pixmap);
            XFreeColormap(theDisplay, cmap);
            XDestroyWindow(theDisplay, theWindow);
            XCloseDisplay(theDisplay);
            
            semctl(semID, 0, IPC_RMID, su);
            semctl(semID2, 0, IPC_RMID, su);
        }
        alive=false;
        pause=false;
        
    }    
#ifndef NO_THREAD
    static void *thread_routine(void *p) ;
    void Run() ;
#endif
    
#ifndef NO_XPM
    void writeXpm(char *name)
    {
        XpmWriteFileFromPixmap(theDisplay,name,pixmap,0,NULL);
    }
    void writegif();
#endif
    void importgif();
    void writeps();
    void testcolor();
    void reversergb();
};

class YWindow: public SYWindow
{
    //all about rotation
public:
    int enableRot;
    void update();
    void applyRot();
    void horizontalRot(double);
    void verticalRot(double);
    void spinRot(double);
    void zoom(double);
    void rotateTo(XEvent);
    void translateTo(XEvent);
    void pbcshiftTo(XEvent,int);
    void pbcglideTo(XEvent);
    void pbcglide(double,double);
    void scaleTo(XEvent);
    void projectTo(XEvent);
    void applyRotate();
    void setWinSpec(int x0,int y0,double s,double d,double a[3][3]);
    YWindow(int w,int h,const char *n,int s=true,int so=true,int fm=true):
        SYWindow(w,h,n,s,so,fm),enableRot(false){};
    void printWinSpec();
    virtual void Evolve(); //overload Evolve to handle more exceptions
    virtual void help();
    virtual void drawBoxFrame();
};



#endif // _DISPLAY_H

