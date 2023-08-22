#pragma once

#ifndef _PDS_SPECTRAL_H
#define _PDS_SPECTRAL_H

#ifdef SPECTRAL

#ifdef PARALLEL
#include <fftw/fftw3-mpi.h>
#else
#include <fftw/fftw3.h>
#endif

#include <vector>
#include <algorithm>
#include <time.h>

#define INCLUDE_NEI 1
#define NO_INCLUDE_NEI 0

class fftcomplex_t
{
public:
    double re;
    double im;
    fftcomplex_t(void) { re = 0.0; im = 0.0; }
    fftcomplex_t(double _re) { re = _re; im = 0.0; }
    fftcomplex_t(double _re, double _im) { re = _re; im = _im; }
    ~fftcomplex_t() {}
    fftcomplex_t& operator=(const fftcomplex_t b) {
        re = b.re;
        im = b.im;
        return (*this);
    }
    fftcomplex_t& operator+=(const fftcomplex_t b) {
        re += b.re;
        im += b.im;
        return (*this);
    }
    fftcomplex_t& operator-=(const fftcomplex_t b) {
        re -= b.re;
        im -= b.im;
        return (*this);
    }
    fftcomplex_t& operator*=(const fftcomplex_t b) {
        re = re*b.re - im*b.im;
        im = im*b.re + re*b.im;
        return (*this);
    }
};

extern fftcomplex_t operator+(const fftcomplex_t a, const fftcomplex_t b);
extern fftcomplex_t operator-(const fftcomplex_t a, const fftcomplex_t b);
extern fftcomplex_t operator*(const fftcomplex_t a, const fftcomplex_t b);
extern fftcomplex_t operator*(const double c, const fftcomplex_t a);
extern fftcomplex_t operator/(const double c, const fftcomplex_t a);
extern fftcomplex_t conj(const fftcomplex_t a);

struct _spectral
{
	int N;
	double rcgrid;
	double rcnei;

    double C[3][3][3][3];

    fftcomplex_t *wk;
    fftcomplex_t *fields[9];
	fftcomplex_t *alpha[9];
	double *stress[6];
    double *dispgrad[9];

	int intNumPoints;
	double *glPositions;
	double *glWeights;

#ifdef PARALLEL
    ptrdiff_t local_Nx;
    ptrdiff_t local_kx_start;
    int Dspan;
#else
    int local_Nx;
    int local_kx_start;
#endif
};

typedef struct _spectral Spectral_t;

typedef struct {
        Node_t *node1;
        Node_t *node2;
        double  p1x, p1y, p1z;
        double  p2x, p2y, p2z;
        double  bx, by, bz;
} NeighborSeg_t;

void InitSpectral(Home_t *home);
void FinishSpectral(Home_t *home);
void InitSpectralCore(Home_t *home);
void FFTCellCharge(Home_t *home);
void SortNodesSpectral(Home_t *home);
void FFTAssignSlice(Home_t *home);
void AlphaTensor(Home_t *home);
void StressFromAlpha(Home_t *home);
void SpectralStress(Home_t *home, int includeNei,
                    double px, double py, double pz, double Stotal[6]);
void SpectralDispGrad(Home_t *home, double px, double py, double pz, double G[9]);
void SpectralForceOneSeg(Home_t *home, Node_t *node1, Node_t *node2,
                         real8 f1[3], real8 f2[3]);

void AllSegmentStress(Home_t *home, real8 x, real8 y, real8 z, real8 Stress[3][3]);

#endif   // SPECTRAL
#endif   // _PDS_SPECTRAL_H
