#pragma once

#ifndef _PDS_DISPLAYC_H
#define _PDS_DISPLAYC_H

/*
  DisplayC.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Sat Apr 13 22:47:33 2002

  FUNCTION  :
*/

#define MAXCOLOR     20
#define COLORNAMELEN 10

extern void AllocColors   ();
extern unsigned long AllocShortRGBColor(unsigned r, unsigned g, unsigned b);
extern void ReadWindowSpec     (char *fname);
extern void Sleep         ();
extern int  WinAlive      ();
extern void WinClear      ();
extern void WinRefresh    ();
extern void WinRoutine    ();
extern void WinEvolve     ();
extern int  WinIsPaused   ();
extern int  WinTogglePause();
extern void WinLock       ();
extern void WinSemRemove  ();
extern void WinUnlock     ();
extern void WinDrawLine   (double x0, double y0, double z0, double x1, double y1, double z1, unsigned long color, double r,unsigned long attr);
extern void WinDrawPoint  (double x , double y , double z , double r , unsigned long color,unsigned long attr);
extern void WinDrawPointS (double x , double y , double z , double r , int radiusInPixels, double lMax, unsigned long color,unsigned long attr,char *s);
extern void WinWritePS    ();

/* Window control variables */
extern int enable_window;
extern char win_name[100];
extern double point_radius, line_width;
extern int win_width, win_height;
extern int sleepseconds;
extern unsigned colors[MAXCOLOR];
extern char color_name[MAXCOLOR][COLORNAMELEN];
extern char bgcolor_name[COLORNAMELEN];
extern int color_scheme; /* 0: color by domain,
                            1: color by Burgers vector */
extern double rotateangles[4]; /* Euler angles alpha, beta, gamma, and scale */
extern double pbcshift[3];
extern double maxpointradius, maxlinewidth;
#endif /* _DISPLAYC_H */

