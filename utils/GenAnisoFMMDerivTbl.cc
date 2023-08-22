//****************************************************************************
//
//      Module:         GenAnisoFMMDerivTbl.cc
//      Description:    This utility is used to create a file containing
//                      the derivatives of the Green's function needed for
//                      doing anisotropic FMM.  The resulting file is
//                      needed by the 'ctablegen' utility which generates
//                      the FMM image correction table, and in the
//                      ParaDiS simulations themselves.
//
//
//      Included functions:
//          main()
//          GetInArgs()
//          InitValues()
//          PrintHelp()
//          Usage()
//
//      Usage:  gen_fmm_derivatives  -mattype {bcc|fcc|hcp|rhombohedral} \
//                  [-c11 <c11>]                            \
//                  [-c12 <c12>]                            \
//                  [-c13 <c13>]                            \
//                  [-c33 <c31>]                            \
//                  [-c44 <c44>]                            \
//                  [-cmatrix <matrixFile>                  \
//                  [-derivfile <derivatives_file>]         \
//                  [-help]                                 \
//                  [-lmax]                                 \
//                  [-mporder <multipole_expansion_order>]  \
//                  [-nphi <numPhiPoints]                   \
//                  [-ntheta <numThetaPoints]               \
//                  [-tporder <taylor_expansion_order>]
//
//                  NOTE: non-cubic cells are not yet fully supported
//                  so the following options are not yet enabled.
//
//                  [-xscale <x_scale_val>]
//                  [-yscale <y_scale_val>]
//                  [-zscale <z_scale_val>]
//
//****************************************************************************

#include <string.h>

#include "Typedefs.h"
#include "Home.h"
#include "Complex_t.h"

#if defined(ANISOTROPIC) && defined(TAYLORFMM)

// note - EPSILON is defined to a different value in Constants.h

#ifdef EPSILON
#undef EPSILON
#endif

#define NUM_TRAPEZOIDAL         50
#define EPSILON                 1.0e-12
#define DEFAULT_DERIV_FILE_NAME "fmm_derivatives.dat"

// ptrcpy() - macro for copying arbitrary pointers 

#define ptrcpy(a,b)  memcpy(&(a),&(b),sizeof(void *))

//  Prototype some local functions

int    gcd0      (int x, int y);
real8  Anm       (int n, int m);
real8  Factorial (int n);
void   zint      (real8 Tin[3], real8 C_3333[3][3][3][3], int numTrapezoidal, real8 zi[3][3]);

//
//  Define an integer identifier to be associated with each
//  posible command line argument.  To be used to index the
//  option-specific data in the optList array below.  OPT_MAX
//  must be the last value in the list

enum 
{
    OPT_C11 = 0    ,
    OPT_C12        ,
    OPT_C13        ,
    OPT_C33        ,
    OPT_C44        ,
    OPT_CMATRIX    ,
    OPT_DERIV_FILE ,
    OPT_HELP       ,
    OPT_LMAX       ,
    OPT_MAT_TYPE   ,
    OPT_MPORDER    ,
    OPT_NUM_PHI    ,
    OPT_NUM_THETA  ,
    OPT_TORDER     ,
#if 0
//  NOTE: non-cubic cells are not yet fully supported
//  so the following options are not yet enabled.

    OPT_XSCALE     ,
    OPT_YSCALE     ,
    OPT_ZSCALE     ,
#endif
    OPT_MAX          // This MUST be the last element in the enumerated list
};

/*
 *  Define a structure to hold a command line option's id (type),
 *  name, the shortest possible unique abbreviation of the option
 *  name, and a flag indicating if the option is paired with
 *  a value or not.
 */
typedef struct {
        int     optType;
        char    *optName;
        int     optMinAbbrev;
        int     optPaired;
} Option_t;


/*
 *  Define and initialize an array of structures containing
 *  all the possible command line arguments, and some info
 *  determining how it is treated.
 *
 *  option            option     #characters     1 unless option
 *  type              name        in unique      has no associated
 *                               abbreviation    value
 */
Option_t  optList[OPT_MAX] = 
{
   { OPT_C11       , (char *) "c11"       , 3, 1 },
   { OPT_C12       , (char *) "c12"       , 3, 1 },
   { OPT_C13       , (char *) "c13"       , 3, 1 },
   { OPT_C33       , (char *) "c33"       , 2, 1 },
   { OPT_C44       , (char *) "c44"       , 2, 1 },
   { OPT_CMATRIX   , (char *) "cmatrix"   , 2, 1 },
   { OPT_DERIV_FILE, (char *) "derivfile" , 1, 1 },
   { OPT_HELP      , (char *) "help"      , 1, 0 },
   { OPT_LMAX      , (char *) "lmax"      , 1, 1 },
   { OPT_MAT_TYPE  , (char *) "mattype"   , 2, 1 },
   { OPT_MPORDER   , (char *) "mporder"   , 2, 1 },
   { OPT_NUM_PHI   , (char *) "nphi"      , 2, 1 },
   { OPT_NUM_THETA , (char *) "ntheta"    , 2, 1 },
   { OPT_TORDER    , (char *) "torder"    , 1, 1 },
#if  0
//  NOTE: non-cubic cells are not yet fully supported
//  so the following options are not yet enabled.

   { OPT_XSCALE    , (char *) "xscale"    , 1, 1 },
   { OPT_YSCALE    , (char *) "yscale"    , 1, 1 },
   { OPT_ZSCALE    , (char *) "zscale"    , 1, 1 },
#endif
};

//
//  Define a structure into which to temporarily store values
//  associated with some of the command line arguments
//
typedef struct {
    int   haveCMatrix;
    int   lMax;
    char  cMatrixFile[MAX_STRING_LEN];
    real8 xScale;
    real8 yScale;
    real8 zScale;
} InputArgs_t;


class Quadruplet {
    public:
        int q,p,L,M;
    public:
        Quadruplet(const Quadruplet &t) : q(t.q),p(t.p),L(t.L),M(t.M) {}
        Quadruplet(int qq,int pp,int LL,int MM) : q(qq),p(pp),L(LL),M(MM) {}

        Quadruplet operator=(const Quadruplet &t) {
            q = t.q;
            p = t.p;
            L = t.L;
            M = t.M;
            return *this;
        }

        int operator==(const Quadruplet &b) const {
            return (q==b.q && p==b.p && L==b.L && M==b.M);
        }

        int hash() const {
            return (q+1)*33331 + (p+1)*333331 + (L*L+L+M+1)*331;
        }

};  // end class Quadruplet 


template<typename T,typename K> class Hash {

    class Link {
        public:
            Link *next;
            K key;
            T data;
            int hits;
            Link(const K &k,Link *n) : next(n),key(k),hits(1) {}
    };

    int isprime(int x) {

        if (x % 2 == 0) {

            return 0;

        } else {
            int k = 3;

            while (k*k <= x) {

                if (x % k == 0) {
                    return 0;
                }

                k = k+2;
            }
            return 1;
        }
    }  // end isprime()

    int    size, used;
    int    maxlength, nstep, nsearch;
    Link **table;

    public:

        class Iterator {
            Hash &H;
            int idx;
            Link *p;
            public:
                Iterator(Hash &HH) : H(HH),idx(-1),p(0) { next(); }
                Iterator(Hash &HH,int iidx,Link *pp) : H(HH),idx(iidx),p(pp) {}
                void next() {
                    if(idx >= H.Size()) return;
                    if(p != 0) p = p->next;
                    if(p == 0) {
                        do {
                            idx = idx+1;
                            if(idx >= H.Size()) return;
                            p = H.table[idx];
                        } while(p == 0);
                    }
                }

                K *key() { return &p->key; }
                T *data() { return &p->data; }
                Link *link() { return p; }

                int operator==(const Iterator &a) {
                    return idx==a.idx && p==a.p;
                }
                int operator!=(const Iterator &a) {
                    return !(*this == a);
                }
        };  // end class Iterator

        Hash(int sz) {
            while(!isprime(sz)) sz = sz + 1;
            size = sz;
            used = 0;


            table = new Link *[size];
            for(int i = 0; i<size; i++)
                table[i] = 0;

            /* Counters for statistics */
            maxlength = 0;
            nstep = 0;
            nsearch = 0;
        }

        ~Hash() {
            for(int i = 0; i<size; i++) {
                Link *p = table[i];
                while(p != 0) {
                    Link *q = p->next;
                    delete p;
                    p = q;
                }
            }
            delete[] table;
        }

        Iterator begin() { return Iterator(*this); }
        Iterator end() { return Iterator(*this,size,0); }

        int Size() const { return size; }
        int Used() const { return used; }
        int NSearch() const { return nsearch; }
        int MaxLength() const { return maxlength; }
        int NStep() const { return nstep; }
        T * Insert(const K &key) {
            int idx = key.hash() % size;
            if(idx < 0) idx = idx + size;
            if(idx >= size || idx < 0) {
                printf("(1) Damn... key = %d, idx = %d, size = %d\n",key.hash(),idx,size);
                exit(1);
            }

            used = used + 1;
            table[idx] = new Link(key,table[idx]);
            return &table[idx]->data;
        }

        void Remove(const K &key) {
            int idx = key.hash() % size;
            Link *p,*last = 0;
            int count = 1;
            if(idx < 0) idx = idx + size;
            if(idx >= size || idx < 0) {
                printf("(2) Damn... key = %d, idx = %d, size = %d\n",key.hash(),idx,size);
                exit(1);
            }

            p = table[idx];
            while(p != 0 && !(p->key == key)) {
                last = p;
                p = p->next;
                count = count + 1;
            }

            if(p != 0) {
                used = used - 1;
                if(last == 0)
                    table[idx] = p->next;
                else
                    last->next = p->next;
                delete p;
            }

            if(count > maxlength)
                maxlength = count;
            nsearch = nsearch + 1;
            nstep = nstep + count;
        }  // end Remove()

        T * Lookup(const K &key) {
            int idx = key.hash() % size;
            Link *p;
            int count = 1;
            if(idx < 0) idx = idx + size;
            if(idx >= size || idx < 0) {
                printf("(3) Damn... key = %d, idx = %d, size = %d\n",key.hash(),idx,size);
                exit(1);
            }


            p = table[idx];
            while(p != 0 && !(p->key == key)) {
                p = p->next;
                count = count + 1;
            }

            if(count > maxlength)
                maxlength = count;
            nsearch = nsearch + 1;
            nstep = nstep + count;

            if(p != 0) p->hits++;

            return (p == 0) ? 0 : &p->data;
        }  // end Lookup()

        Hash<T,K> * dup() {
            Hash<T,K> *H = new Hash<T,K>(Used() * 2);

            for(Iterator iter = begin(); iter != end(); iter.next()) {
                *( H->Insert(*iter.key()) ) = *iter.data();
            }
            return H;
        }  // end dup()

};  // end template<typename T,typename K> class Hash

typedef Hash<complex8,Quadruplet> THash;
typedef THash::Iterator TIter;

// --------------------------------------------------------------------------
//
//  Function:    S()
//  Description: Returns the sign of the specified value
//
// --------------------------------------------------------------------------
int S(int x)
{
    return((x >= 0) ? 1 : -1);
}


// --------------------------------------------------------------------------
//
//  Function:    add()
//  Description: Add the specified real and imaginary coeficient components
//               to the coefficents to the Y table.
//
//  Parameters:
//      IN:   YTbl
//      IN:   q
//      IN:   p
//      IN:   L
//      IN:   M
//      IN:   rex  Real component of the coefficient
//      IN    imx  Imaginary component of the coefficient
//
// --------------------------------------------------------------------------
inline void add(THash *YTbl, int q, int p, int L, int M, real8 rex, real8 imx)
{
    Quadruplet  key(q,p,L,M);
    complex8   *ptr = YTbl->Lookup(key);

    if (ptr != 0) {
        if (ptr->re == -rex && ptr->im == -imx) {
            YTbl->Remove(key);
        } else {
            ptr->re += rex;
            ptr->im += imx;
        }
    } else if (rex != 0.0 || imx != 0.0) {
        ptr = YTbl->Insert(key);
        ptr->re = rex;
        ptr->im = imx;
    }
}


// --------------------------------------------------------------------------
//
//  Function:    deriv_X()
//  Description: 
//
//  Parameters:
//      IN:   coord
//      IN:   YTbl
//
// --------------------------------------------------------------------------
THash * deriv_X(int coord,THash *YTbl)
{
    THash *dY = new THash(YTbl->Size());
    TIter  iter = YTbl->begin();

    while (iter != YTbl->end()) {
        int L1qp;
        int q = iter.key()->q;
        int p = iter.key()->p;
        int L = iter.key()->L;
        int M = iter.key()->M;
        real8 rec = iter.data()->re;
        real8 imc = iter.data()->im;

        if (p > 0) {
            add(dY,q+1,p-1,L,M,p*rec,p*imc);
        }

        L1qp = L+1-q-p;

        if (L1qp != 0) {
            add(dY,q+1,p+1,L,M,L1qp*rec,L1qp*imc);
        }

        if (coord == 1) {
            /* d Y(L,M) / dx */
            add(dY,q+1,p,L+1,M+1,  0.5*S( M)*rec, 0.5*S( M)*imc);
            add(dY,q+1,p,L+1,M-1,  0.5*S(-M)*rec, 0.5*S(-M)*imc);
        } else if (coord == 2) {
            /* d Y(L,M) / dy */
            add(dY,q+1,p,L+1,M+1,  0.5*S( M)*imc,-0.5*S( M)*rec);
            add(dY,q+1,p,L+1,M-1, -0.5*S(-M)*imc, 0.5*S(-M)*rec);

        } else if (coord == 3) {
            /* d Y(L,M) / dz */
            add(dY,q+1,p,L+1,M,rec,imc);
        } else {
            printf("Unknown coordinate direction coord=%d.\n",coord);
            exit(1);
        }

        iter.next();
    }

    return dY;
}  // end deriv_X()


// --------------------------------------------------------------------------
//
//  Function:    eval_partial()
//  Description: 
//
//  Parameters:
//      IN:   YTbl
//      IN:   x
//
// --------------------------------------------------------------------------
THash * eval_partial(THash *YTbl,real8 x)
{
    THash *Yout = new THash(YTbl->Size());
    TIter iter = YTbl->begin();

    while (iter != YTbl->end()) {
        int q = iter.key()->q;
        int p = iter.key()->p;
        int L = iter.key()->L;
        int M = iter.key()->M;
        real8 rec = iter.data()->re;
        real8 imc = iter.data()->im;
        real8 xp = ipow(x,p);

        add(Yout,q+p,0,L,M,xp * rec,xp * imc);
        iter.next();
    }

    return Yout;
}


// --------------------------------------------------------------------------
//
//  Function:    cplmall()
//  Description: 
//
//  Parameters:
//      IN:   lMax
//      IN:   x
//      OUT:  plm
//
// --------------------------------------------------------------------------
void cplmall(int lMax, real8 x, real8 *plm)
{
    int   L, M, idx1, idx2, idx3;
    int   plmMaxIndex;
    real8 sine;

    plmMaxIndex = (lMax+1) * (lMax+2) / 2;

    for (L = 1; L <= plmMaxIndex; L++) {
        plm[L-1] = -99.0;
    }

    sine = sqrt((1.0+x)*(1.0-x));
    plm[0] = 1.0;

    for (M = 0; M <= lMax; M++) {

        /* P[M,M] */
        idx1 = M * (M+1) / 2 + M + 1;
        if (M > 0) {
            plm[idx1-1] = -plm[(M-1)*M/2+(M-1)] * (2*M-1) * sine;
        }

        /* P[M+1,M] */
        idx2 = (M+1)*(M+2)/2 + M + 1;
        if (M < lMax) {
            plm[idx2-1] = x*(2*M+1)*plm[idx1-1];
        }

        for (L = M+2; L <= lMax; L++) {
            /* P[L,M] */
            idx3 = L * (L+1) / 2 + M + 1;
            plm[idx3-1] = (x*(2*L-1)*plm[idx2-1]-(L+M-1)*plm[idx1-1]) / (L-M);

            idx1 = idx2;
            idx2 = idx3;
        }
    }

/*
 *  The following is just a sanity check...
 */
#if 1
    idx2 = 0;

    for (L = 0; L <= lMax; L++) {
        for (M = 0; M <= L; M++) {
            idx2 = idx2 + 1;
            idx1 = L *(L+1) / 2 + M + 1;
            if (idx1 != idx2) {
                printf("plmall: Index error idx1=%d, idx2=%d, L=%d, M=%d",
                      idx1, idx2, L, M);
            }

            if (plm[idx1-1] == -99.0) {
                printf("plmall: Untouched index idx1=%d, L=%d, M=%d",
                      idx1, L, M);
            }
        }
    }
#endif

    return;
}


/*
 *  Parameters:
 *      ylm: array of size (lmax+1)**2
 */
// --------------------------------------------------------------------------
//
//  Function:    ylmall()
//  Description: 
//
//  Parameters:
//      IN:   lMax
//      IN:   theta
//      IN:   phi
//
//  Returns:  array ylm[(lmax+1)*(lmax+1)]
//
// --------------------------------------------------------------------------
complex8 *ylmall(int lmax, real8 theta, real8 phi)
{
    int         L, M, idx;
    real8       x, c, s, f;
    real8       plm[(lmax+1)*(lmax+2)/2];
    real8       rSqrt[2*lmax+3];
    complex8   *ylm = new complex8[(lmax+1)*(lmax+1)];
    complex8    im, expi[lmax+1];

    im = complex8(0.0, 1.0);

    x = cos(theta);
    c = cos(phi);
    s = sin(phi);

    cplmall(lmax, x, plm);

    ylm[0].re = 1.0;
    ylm[0].im = 0.0;
    expi[0] = complex8(1.0,0.0);

    rSqrt[0] = 1.0;

    for (L = 1; L <= lmax; L++) {

        rSqrt[2*L-1] = sqrt(1.0 / (2*L-1));
        rSqrt[2*L]   = sqrt(1.0 / (2*L));
        expi[L] = expi[L-1] * (c + im * s);

        f = 1.0;
        idx = L*L + L + 1;
        ylm[idx-1].re = creal(f * plm[L*(L+1)/2] * expi[0]);
        ylm[idx-1].im = cimag(f * plm[L*(L+1)/2] * expi[0]);

        for (M = 1; M <= L; M++) {
            complex8 tmpVal;

            f = f * rSqrt[L-M+1]*rSqrt[L+M];

            tmpVal = f * plm[L*(L+1)/2 + M] * expi[M];

            ylm[idx + M - 1].re = creal(tmpVal);
            ylm[idx + M - 1].im = cimag(tmpVal);
            ylm[idx - M - 1].re = ylm[idx + M -1].re;
            ylm[idx - M - 1].im = -1.0 * ylm[idx + M -1].im;
        }
    }

    return ylm;
}


// --------------------------------------------------------------------------
//
//  Function:    ylm2zlm()
//  Description: 
//
//  Parameters:
//      IN:   lMax
//      IN:   ylm[(lmax+1)*(lmax+1)]
//
//  Returns:  array zlm[(lmax+1)*(lmax+1)]
//
// --------------------------------------------------------------------------
complex8 *ylm2zlm(int lmax, complex8 ylm[])
{
    int   k, l, m;
    complex8 *zlm = new complex8[(lmax+1)*(lmax+1)];

    k = 0;

    for (l = 0; l <= lmax; l++) {
        for (m = -l; m <= l; m++) {
            k = k+1;
            zlm[k-1].re = ylm[k-1].re / Anm(l, m);
            zlm[k-1].im = ylm[k-1].im / Anm(l, m);
        }
    }

    return zlm;
}


// --------------------------------------------------------------------------
//
//  Function:    Usage()
//  Description: Print out a brief message indicating the possible
//               command line options.
//
//  Parameters:
//      IN:   progName  Name of program being executed
//
// --------------------------------------------------------------------------
void Usage(char *progName)
{
    printf("Usage:  %s  -mattype {bcc|fcc|hcp|rhombohedral} \\\n", progName);
    printf("            [-c11 <c11>]                            \\\n");
    printf("            [-c12 <c12>]                            \\\n");
    printf("            [-c13 <c13>]                            \\\n");
    printf("            [-c33 <c31>]                            \\\n");
    printf("            [-c44 <c44>]                            \\\n");
    printf("            [-cmatrix <matrixFile>                  \\\n");
    printf("            [-derivfile <derivatives_file>]         \\\n");
    printf("            [-help]                                 \\\n");
    printf("            [-lmax]                                 \\\n");
    printf("            [-mporder <multipole_expansion_order>]  \\\n");
    printf("            [-nphi <numPhiPoints]                   \\\n");
    printf("            [-ntheta <numThetaPoints]               \\\n");
#if 0
//
//  Non-cubic cells are not yet fully supported, so these options are
//  not available yet.
//
    printf("            [-xscale <xscale_val]                   \\\n");
    printf("            [-yscale <yscale_val]                   \\\n");
    printf("            [-zscale <zscale_val]                   \n");
#endif

    return;
}


// --------------------------------------------------------------------------
//
//  Function:    PrintHelp()
//  Description: Print out a detailed description of the available
//               program options, what they represent, how they
//               relate to each other, which are interdependent, or
//               mutually exclusive, default values, etc.
//               command line options and terminates.
//
//  Parameters:
//      IN:   progName  Name of program being executed
//
// --------------------------------------------------------------------------
void PrintHelp(char *progName)
{
    Usage(progName);

    printf("    Options may be abbreviated to the shortest non-ambiguous\n");
    printf("    abbreviation of the option.\n\n");
    printf("Options:\n\n");
    printf("  -mattype  Required argument.  Indicates the crystal \n");
    printf("            structure of the material \n\n");
    printf("  -c11      Elastic constant.  Applies when <mattype> is \n");
    printf("            bcc, fcc or hcp.  Ignored if <cmatrix> is provided\n");
    printf("            No default value.\n\n");
    printf("  -c12      Elastic constant.  Applies when <mattype> is \n");
    printf("            bcc, fcc or hcp.  Ignored if <cmatrix> is provided\n");
    printf("            No default value.\n\n");
    printf("  -c13      Elastic constant.  Applies when <mattype> is \n");
    printf("            hcp.  Ignored if <cmatrix> is provided\n");
    printf("            No default value.\n\n");
    printf("  -c33      Elastic constant.  Applies when <mattype> is \n");
    printf("            hcp.  Ignored if <cmatrix> is provided\n");
    printf("            No default value.\n\n");
    printf("  -c44      Elastic constant.  Applies when <mattype> is \n");
    printf("            bcc, fcc or hcp.  Ignored if <cmatrix> is provided\n");
    printf("            No default value.\n\n");
    printf("  -cmatrix  Specifies the name of a file containing the\n");
    printf("            appropriate 6X6 elastic constant matrix.  If \n");
    printf("            provided, these constants override any individual\n");
    printf("            elastic constants provided.  Required if <mattype>\n");
    printf("            is rhombohedral.  No default value.  The format\n");
    printf("            of this file is 36 values in the following order\n");
    printf("            where white space is ignored and anything \n");
    printf("            following a # is a comment and ignored:\n");
    printf("            C[0][0], C[0][1], ... C[0][5], C[1][0], \n");
    printf("            C[1][1], ... C[5][5]\n\n");
    printf("  -derivfile  Specifies the name of the file into which to\n");
    printf("            write the Green's function derivatives.\n");
    printf("            Defaults to %s \n", DEFAULT_DERIV_FILE_NAME);
    printf("  -help     Prints basic usage information.\n\n");
    printf("  -lmax     Defines the number of terms in the spherical \n");
    printf("            harmonics.  This value is related to the ParaDiS\n");
    printf("            control parameter <anisoHarmonicsNumTermsBase>\n");
    printf("            and should be at least 2*<anisoHarmonicsNumTermsBase>+1\n");
    printf("            Default is 12. \n\n");
    printf("  -mporder  Multipole expansion order.  This value should\n");
    printf("            be at least the same as the expected simulation\n");
    printf("            multipole expansion order.  Default value is 5\n\n");
    printf("  -nphi     Number of points into which the phi angles are\n");
    printf("            discretized when defining the double integral \n");
    printf("            over theta and phi angles when calculating the glm\n");
    printf("            expansion coefficients.  Default is 40.\n\n");
    printf("  -ntheta   Number of points into which the theta angles are\n");
    printf("            discretized when defining the double integral \n");
    printf("            over theta and phi angles when calculating the glm\n");
    printf("            expansion coefficients.  Default is 40.\n\n");
    printf("  -torder   Taylor expansion order.  This value should\n");
    printf("            be at least the same as the expected simulation\n");
    printf("            taylor expansion order.  Default value is 11\n\n");

#if 0
//
//  Non-cubic cells are not yet fully supported, so these options are
//  not available yet
//
    printf("  -xscale,  \n");
    printf("  -yscale,  \n");
    printf("  -zscale   These values define the aspect ratios of the cells\n");
    printf("            For example, x, y, and z scale values of 1,1,1\n");
    printf("            indicate a cubic cell while 1,2,3 means that the\n");
    printf("            cells are 2 times longer in Y than X and 3 times\n");
    printf("            longer in Z than X.\n");
#endif

    exit(1);
}


// --------------------------------------------------------------------------
//
//  Function:    Factorial()
//  Description: Calculates factorial of specified n value.
//
//  Returns:  n!
//
// --------------------------------------------------------------------------
real8 Factorial(int n)
{
    int   i;
    real8 f;

    f = 1.0;

    for (i = 1; i <= n; i++) {
        f = f * i;
    }

    return(f);
}

real8 Anm(int n, int m)
{
    real8 rVal;

    rVal = pow(-1.0, (n+m)) / sqrt(Factorial(n-m) * Factorial(n+m));

    return(rVal);
}


// --------------------------------------------------------------------------
//
//  Function:    zint()
//  Description: 
//
//  Parameters:
//      IN:    Tin
//      IN:    C_3333
//      IN:    numTrapezoidal
//      OUT:   zi
//
// --------------------------------------------------------------------------
void zint(real8 Tin[3], real8 C_3333[3][3][3][3], int numTrapezoidal,
          real8 zi[3][3])
{
    int   i, j, k, m, n;
    int   indexOfMin;
    real8 minVal;
    real8 xHatDotT;
    real8 yHat[3], Tlocal[3];
    real8 xHat[3] = {0.0, 0.0, 0.0};


    Tlocal[0] = Tin[0];
    Tlocal[1] = Tin[1];
    Tlocal[2] = Tin[2];

//
//  Zero out zi to start
//
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            zi[i][j] = 0.0;
        }
    }

//
//  Normalize T and find the components with the minimum absolute
//  value.
//
    NormalizeVec(Tlocal);

    indexOfMin = 0;
    minVal = fabs(Tlocal[0]);

    for (i = 1; i < 3; i++) {
        if (fabs(Tlocal[i]) < minVal) {
            minVal = fabs(Tlocal[i]);
            indexOfMin = i;
        }
    }

    xHat[indexOfMin] = 1.0;
    xHatDotT = DotProduct(Tlocal, xHat);

    for (i = 0; i < 3; i++) {
        xHat[i] = xHat[i] - xHatDotT * Tlocal[i];
    }

    NormalizeVec(xHat);
    cross(Tlocal, xHat, yHat);

    for (i = 1; i <= numTrapezoidal; i++) {
        real8 phi;
        real8 zHat[3];
        real8 cz[3][3], czInv[3][3];

        memset(&cz[0][0], 0, sizeof(cz));

        phi = 2.0 * M_PI * (real8)i / numTrapezoidal;

        for (j = 0; j < 3; j++) {
            zHat[j] = cos(phi) * xHat[j] + sin(phi) * yHat[j];
        }

        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                for (m = 0; m < 3; m++) {
                    for (n = 0; n < 3; n++) {
                        cz[m][n] += C_3333[j][m][n][k] * zHat[j] * zHat[k];
                    }
                }
            }
        }

        Matrix33_Inverse(czInv,cz);

        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {
                zi[j][k] += czInv[j][k];
            }
        }
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            zi[i][j] = (zi[i][j] * 2.0 * M_PI / numTrapezoidal) /
                       (8.0 * (M_PI * M_PI));
        }
    }

    return;
}


// --------------------------------------------------------------------------
//
//  Function:    cart2polar()
//  Description: 
//
//  Parameters:
//      IN:    x
//      IN:    y
//      IN:    z
//      OUT:   r_p
//      OUT:   theta_p
//      OUT:   phi_p
//
// --------------------------------------------------------------------------
static void cart2polar(real8 x, real8 y, real8 z, real8 r_p[1],
                       real8 theta_p[1], real8 phi_p[1])
{
    real8 r, rxy, theta, phi;

    r = sqrt(x*x + y*y + z*z);

    if (r == 0.0) {
        Fatal("cart2polar(): zero vector is invalid");
    }

    rxy   = sqrt(x*x + y*y);
    theta = acos(z/r);

    if (rxy == 0.0) {
        phi = 0.0;
    } else {
        phi = atan2(y/rxy, x/rxy);
    }

    r_p[0]     = r;
    theta_p[0] = theta;
    phi_p[0]   = phi;

    return;
}


// --------------------------------------------------------------------------
//
//  Function:    EvalGreensFuncOnGrid()
//  Description: 
//
//  Parameters:
//      IN:    numTheta        Number of angles in <theta> array
//      IN:    theta[numTheta] Array of theta angles
//      IN:    numPhi          Number of angles in <phi> array
//      IN:    phi[numphi]     Array of phi angles
//      OUT:   G   
//      IN:    numTrapezoidal
//      IN:    C_3333          Elastic constants
//
// --------------------------------------------------------------------------
void EvalGreensFuncOnGrid(int numTheta, real8 *theta,
                          int numPhi, real8 *phi, real8 *G, int numTrapezoidal,
                          real8 C_3333[3][3][3][3])
{
    real8 (*GPtr)[3][numTheta][numPhi]=0;   ptrcpy(GPtr,G);

    for (int i=0; i < numTheta; i++) 
    {
        for (int j=0; j < numPhi; j++) 
        {
            real8 sinTheta;
            real8 T[3];
            real8 zi[3][3];

            sinTheta = sin(theta[i]);

            T[0] = sinTheta * cos(phi[j]);
            T[1] = sinTheta * sin(phi[j]);
            T[2] = cos(theta[i]);

            zint(T, C_3333, numTrapezoidal, zi);

            for (int m=0; m<3; m++) 
            for (int n=0; n<3; n++)
                GPtr[m][n][i][j] = zi[m][n];
        }
    }

    return;
}


// --------------------------------------------------------------------------
//
//  Function:    InitValues()
//  Description: Initialize default values for command line arguments
//               and such.
//
// --------------------------------------------------------------------------
void InitValues(Home_t *home, InputArgs_t *inArgs)
{

    inArgs->cMatrixFile[0] = 0;

    home->param->C11 = -1.0;
    home->param->C12 = -1.0;
    home->param->C13 = -1.0;
    home->param->C33 = -1.0;
    home->param->C44 = -1.0;

    strcpy(home->param->fmAnisoDerivTbl, (char *) DEFAULT_DERIV_FILE_NAME);

    home->param->fmMPOrder = 5;
    home->param->fmExpansionOrder = 11;

    home->param->materialType = MAT_TYPE_UNDEFINED;

    home->param->anisoNumPhiPoints = 40;
    home->param->anisoNumThetaPoints = 40;

    inArgs->xScale = 1.0;
    inArgs->yScale = 1.0;
    inArgs->zScale = 1.0;

    inArgs->lMax = 12;

    return;
}


//---------------------------------------------------------------------------
//
//  Function:  GetInArgs
//  Description: Parse and process and user supplied command line
//               options.  Set appropriate default values (if
//               permitted) for any values not provided by the
//               user.
//
//  Arguments:
//      argc       count of command line args
//      argv       command line args
//      inArgs     struct in which to temporarily store values
//                 associated with the command line arguments.
//
//-------------------------------------------------------------------------*/
static void GetInArgs(int argc, char *argv[], Home_t *home,
                      InputArgs_t *inArgs)
{
    int      i, j;
    int      maxCopyLen;
    char    *argName;
    char    *argValue;
    Param_t *param;

    param = home->param;

    for (i = 1; i < argc; i++) {

//
//      If the option doesn't begin with a '-' something
//      is wrong, so notify the user and terminate.
//
        if (argv[i][0] != '-') {
            printf("Error: Invalid option '%s'\n", argv[i]);
            Usage(argv[0]);
            exit(1);
        }

        argName = &argv[i][1];

//
//      Scan the array of valid options for the user supplied
//      option name.  (This may be any unique abbreviation of
//      one of any of the options.  If we don't find the option
//      notify the user and terminate.
//
        for (j = 0; j < OPT_MAX; j++) {
            if (!strncmp(argName, optList[j].optName,
                         optList[j].optMinAbbrev)) {
                break;
            }
        }

        if (j == OPT_MAX) {
            Usage(argv[0]);
            exit(1);
        }

//
//      Verify that there is an associated value if the specified
//      option is supposed to be paired with a value.
//
        if (optList[j].optPaired && (i+1 >= argc)) {
            Usage(argv[0]);
            exit(1);
        } else {
            argValue = argv[++i];
        }

//
//      Do any option-specific processing...
//

        switch (optList[j].optType)  {
            case OPT_C11:
                param->C11 = atof(argValue);
                break;

            case OPT_C12:
                param->C12 = atof(argValue);
                break;

            case OPT_C13:
                param->C13 = atof(argValue);
                break;

            case OPT_C33:
                param->C33 = atof(argValue);
                break;

            case OPT_C44:
                param->C44 = atof(argValue);
                break;

            case OPT_CMATRIX:
                maxCopyLen = sizeof(inArgs->cMatrixFile) - 1;
                strncpy(inArgs->cMatrixFile, argValue, maxCopyLen);
                inArgs->cMatrixFile[maxCopyLen] = 0;
                inArgs->haveCMatrix = 1;
                break;

            case OPT_DERIV_FILE:
                maxCopyLen = sizeof(param->fmAnisoDerivTbl) - 1;
                strncpy(param->fmAnisoDerivTbl, argValue, maxCopyLen);
                param->fmAnisoDerivTbl[maxCopyLen] = 0;
                break;

            case OPT_HELP:
                PrintHelp(argv[0]);
                exit(0);
                break;

            case OPT_LMAX:
                inArgs->lMax = atoi(argValue);
                break;

            case OPT_MAT_TYPE:
                if (StrEquiv(argValue, "bcc")) {
                    param->materialType = MAT_TYPE_BCC;
                    strcpy(param->materialTypeName, "BCC");
                } else if (StrEquiv(argValue, "fcc")) {
                    param->materialType = MAT_TYPE_FCC;
                    strcpy(param->materialTypeName, "FCC");
                } else if (StrEquiv(argValue, "hcp")) {
                    param->materialType = MAT_TYPE_HCP;
                    strcpy(param->materialTypeName, "HCP");
                } else if (StrEquiv(argValue, "rhombohedral")) {
                    param->materialType = MAT_TYPE_RHOMBOHEDRAL_VA;
                    strcpy(param->materialTypeName, "RHOMBOHEDRAL");
                } else {
                    Fatal("Invalid option '-mattype %s'", argValue);
                }
                break;

            case OPT_MPORDER:
                param->fmMPOrder = atoi(argValue);
                break;

            case OPT_NUM_PHI:
                param->anisoNumPhiPoints = atoi(argValue);
                break;

            case OPT_NUM_THETA:
                param->anisoNumThetaPoints = atoi(argValue);
                break;

            case OPT_TORDER:
                param->fmExpansionOrder = atoi(argValue);
                break;

        }  // end switch

    }  // end for (i = 1; i < argc; ...)

//
//  Now do some sanity checks on inputs and/or defaults
//
    if (param->anisoNumPhiPoints <= 0) {
        printf("Invalid value for option: '-%s %d'\n",
              optList[OPT_NUM_PHI].optName, param->anisoNumPhiPoints);
        exit(1);
    }

    if (param->anisoNumThetaPoints <= 0) {
        printf("Invalid value for option: '-%s %d'\n",
              optList[OPT_NUM_THETA].optName, param->anisoNumThetaPoints);
        exit(1);
    }

    if (param->fmMPOrder <= 0) {
        printf("Invalid value for option: '-%s %d'\n",
              optList[OPT_MPORDER].optName, param->fmMPOrder);
        exit(1);
    }

    if (param->fmExpansionOrder <= 0) {
        printf("Invalid value for option: '-%s %d'\n",
              optList[OPT_MPORDER].optName, param->fmExpansionOrder);
        exit(1);
    }

//
//  The elastic constants required depend on the type of material.
//  If the full 6x6 elastic constant matrix has been provided, that
//  is sufficient and will take precedence over any individual
//  constants specified.
//
//  For BCC and FCC, the full elastic constant matrix may be provided
//  or all of c11, c12 and c44 may be specified
//
//  For HCP, the full elastic constant matrix may be provided, or all
//  of c11, c12, c13, c33 and c44 may be specified
//
//  For RHOMBOHEDRAL, the full elastic constant matrix must be provided
//
    if (inArgs->haveCMatrix == 0) {

        if ((param->materialType == MAT_TYPE_BCC) ||
            (param->materialType == MAT_TYPE_FCC)) {

            if ((param->C11 <= 0.0) ||
                (param->C12 <= 0.0) ||
                (param->C44 <= 0.0)) {
                printf("Error: for %s, need either c11, c12 and c44 "
                       "constants\n or the full elastic constant "
                       "matrix\n", param->materialTypeName);
                exit(1);
            }
        }

        if (param->materialType == MAT_TYPE_HCP) {
            if ((param->C11 <= 0.0) ||
                (param->C12 <= 0.0) ||
                (param->C13 <= 0.0) ||
                (param->C33 <= 0.0) ||
                (param->C44 <= 0.0)) {
                printf("Error: for %s, need either c11, c12, c13, c33, "
                       "and c44\nconstants or the full elastic constant "
                       "matrix\n", param->materialTypeName);
                exit(1);
            }
        }

        if (param->materialType == MAT_TYPE_RHOMBOHEDRAL_VA) {
            printf("Error: for %s, the full elastic constant matrix\n"
                   "is required.\n", param->materialTypeName);
            exit(1);
        }

    }  // end if (haveCMatrix == 0)

//
//  Make sure the scale components are positive values
//
    if ((inArgs->xScale <= 0.0) ||
        (inArgs->yScale <= 0.0) ||
        (inArgs->zScale <= 0.0)) {
        printf("Error: Scale is invalid.  Scale values must be >= 0.0\n");
        exit(1);
    }

    return;
}


//
//  Compute derivatives, iterate over cells, and evaluate
//
void ComputeDerivAndEvaluate(int iIndex, int jIndex, FILE *fpDeriv,
                             int maxOrder, THash *YTbl,
                             real8 xscale, real8 yscale, real8 zscale)
{
    THash *ziter = YTbl->dup();
    int    seriesExpansionSize;

//
//  Compute the size of the input series expansion
//
    {
        int    lmax_series = 0;
        for(TIter ii = YTbl->begin(); ii != YTbl->end(); ii.next()) {
            int L = ii.key()->L;
            if(L > lmax_series) lmax_series = L;
        }

        seriesExpansionSize = lmax_series + maxOrder;
    }


//
//  To go from a finer level to a coarser, lump nlump x nlump x nlump cells
//  together.  nlump is usually 2, but 3 for periodic boundary lattice sum
//  computation.
//
//  Near field is up to and including nskip boxes away.
//
    const int nskip = 1,nlump = 3;
    const int irange = nlump*nskip + (nlump/2);//(nskip+1)*nlump-1;
    complex8 *zlmarr_v[2*irange+1];

    for(int ndz = 0; ndz <= maxOrder; ndz++) {/* Z derivative loop */
        for(int iz = -irange; iz<=irange; iz++) {/* Z cell loop */
            THash *yiter = eval_partial(ziter,iz*zscale);
	
            for(int ndy = 0; ndy <= maxOrder-ndz; ndy++) {/* Y derivative loop */
                for(int iy = -irange; iy<=irange; iy++) {/* Y cell loop */
                    THash *xiter = eval_partial(yiter,iy*yscale);

                    // Compute spherical harmonics
                    for(int i = 0; i<2*irange+1; i++) {
                        real8 r,theta,phi;
                        complex8 *ylmarr;
                        real8 ix = i - irange;

                        if(ix == 0 && iy == 0 && iz == 0) {
                            zlmarr_v[i] = new complex8[1];
                            continue;
                        }

                        cart2polar(ix*xscale,iy*yscale,iz*zscale,&r,&theta,&phi);
                        ylmarr = ylmall(seriesExpansionSize,theta,phi);
                        zlmarr_v[i] = ylm2zlm(seriesExpansionSize,ylmarr);

                        delete[] ylmarr;
                    }

                    /* X derivative loop */
                    for(int ndx = 0;ndx <= maxOrder-ndz-ndy;ndx++) {

                        /* X cell loop */
                        for(int ix = -irange; ix<=irange; ix++) {
                            THash *ylm_expr;
                            if(ix == 0 && iy == 0 && iz == 0) {
                                continue;
                            }

                            if(gcd0(ix,gcd0(iy,iz)) > 1) continue;

                            ylm_expr = eval_partial(xiter,ix*xscale);

                            /* Evaluate ylm_expr using precomputed Ylm */
                            {
                                complex8 *zlm = zlmarr_v[ix+irange];
                                real8 ret = 0.0,imt = 0.0;
                                TIter iter = ylm_expr->begin();
                                while(iter != ylm_expr->end()) {
                                    int q = iter.key()->q;
                                    int p = iter.key()->p;
                                    int L = iter.key()->L;
                                    int M = iter.key()->M;
                                    real8 x = ix*xscale,y = iy*yscale,z = iz*zscale;
                                    real8 r = sqrt(x*x + y*y + z*z);
                                    real8 rq = ipow(r,-q);
                                    real8 rec = iter.data()->re * rq;
                                    real8 imc = iter.data()->im * rq;
                                    int idx = L*L + L + M;

                                    ret = ret + rec*zlm[idx].re - imc*zlm[idx].im;
                                    imt = imt + rec*zlm[idx].im + imc*zlm[idx].re;

                                    if(p > 0) {
                                        Fatal("ComputeDerivAndEvaluate: p is "
                                              "supposed to be zero here. "
                                              "p = %d. Exiting.\n",p);
                                    }

                                    iter.next();
                                }
		
                                /* Output value */
                                {
                                    fprintf(fpDeriv,
                                            "%d %d  %2d %2d %2d  %2d %2d %2d  "
                                            "%25.15e\n", iIndex, jIndex,
                                            ndx,ndy,ndz,ix,iy,iz,ret);
                                }

                            } /* Evaluation of ylm_expr */

                            delete ylm_expr;

                        } /* X cell loop */

                        /* Differentiate in X direction */
                        {
                            THash *dxiter = deriv_X(1,xiter);
                            delete xiter;
                            xiter = dxiter;
                        }

                    } /* X derivative loop */

                    delete xiter;

                    for(int i = 0; i<2*irange+1; i++) {
                        delete[] zlmarr_v[i];
                    }

                } /* Y cell loop */

                /* Differentiate in Y direction */
                {
                    THash *dyiter = deriv_X(2,yiter);
                    delete yiter;
                    yiter = dyiter;
                }

            } /* Y derivative loop */

            delete yiter;

        } /* Z cell loop */

        /* Differentiate in Z direction */
        {
            THash *dziter = deriv_X(3,ziter);
            delete ziter;
            ziter = dziter;
        }

    } /* Z derivative loop */

    delete ziter;
}


int main(int argc, char *argv[])
{
    int          i, j, k, m, n;
    int          lMax, maxOrder;
    int          numYCoeff, numGVals;
    int          numTheta, numPhi, numTrapezoidal;
    real8        xScale, yScale, zScale;
    real8       *GAll;      /* dimensions 3 X 3 X numTheta X numPhi */
    Home_t      *home;
    Param_t     *param;
    InputArgs_t *inArgs;
    AnisotropicVars_t *anisoVars;
    FILE        *fpDeriv;


    home = (Home_t *)calloc(1, sizeof(Home_t));
    home->param = (Param_t *)calloc(1, sizeof(Param_t));
    home->fmAnisoTable = (FMAnisoTable_t *)calloc(1, sizeof(FMAnisoTable_t));
    inArgs = (InputArgs_t *)calloc(1, sizeof(InputArgs_t));

    param = home->param;
    anisoVars = &home->anisoVars;

    InitValues(home, inArgs);
    GetInArgs(argc, argv, home, inArgs);

    lMax = inArgs->lMax;
    xScale = inArgs->xScale;
    yScale = inArgs->yScale;
    zScale = inArgs->zScale;
    numPhi = param->anisoNumPhiPoints;
    numTheta = param->anisoNumThetaPoints;
    numTrapezoidal = NUM_TRAPEZOIDAL;
    

//
//  If the user specified a file containing the full elastic constant
//  matrix, read it in now.
//
    if (inArgs->haveCMatrix) {
        int   tokenType;
        char  token[256];
        FILE *fpCMatrix;

        fpCMatrix = fopen(inArgs->cMatrixFile, "r");

        if (fpCMatrix == (FILE *)NULL) {
            Fatal("Error opening elastic constants matrix file %s\n",
                  inArgs->cMatrixFile);
        }

//
//      The elastic constant matrix consists of 36 values specified
//      in the following order.  White space is ignored and anything
//      following a '#' on a line is treated as a comment and ignored.
//          C[0][0],
//          C[0][1],
//          ...
//          C[0][5],
//          C[1][0],
//          C[1][1],
//          ...
//          C[5][5]
//
        for (i = 0; i < 6; i++) {
            for (j = 0; j < 6; j++) {

                tokenType = GetNextToken(fpCMatrix, token, sizeof(token));

                switch (tokenType) {
//
//                  If we find any non-comment data that is not a valid
//                  value, we have a problem.
//
                    case TOKEN_ERR:
                    case TOKEN_NULL:
                    case TOKEN_EQUAL:
                    case TOKEN_BEGIN_VAL_LIST:
                    case TOKEN_END_VAL_LIST:
                        Fatal("Error parsing C matrix file %s",
                              inArgs->cMatrixFile);
                        break;
                    default:
                        param->elasticConstantMatrix[i][j] = atof(token);
                }
            }
        }

        fclose(fpCMatrix);

    } else {
//
//      User did not provide the full elastic constant matrix, so we
//      need to create it from the individual elastic constants
//
        SetElasticConstantMatrix(home);
    }

//
//  We have the 6x6 elastic constant matrix C, but we'll also need 
//  the 3x3x3x3 form of C...
//
    SetElasticConstantMatrix4D(param->elasticConstantMatrix,
                               anisoVars->elasticConstantMatrix4D);


//
//  Set the order of Ylm needed
//
    maxOrder = MAX(param->fmMPOrder, param->fmExpansionOrder) + 1;

//
//  Open up the file for the Green's function derivatives and write
//  some initial metadata into the file.
//
    if ((fpDeriv = fopen(param->fmAnisoDerivTbl, "w")) == (FILE *)NULL) {
        Fatal("Error opening derivatives output file %s",
              param->fmAnisoDerivTbl);
    }

    fprintf(fpDeriv,
            "%%%% maxorder = %d;\n"
            "%%%% xscale = %25.15e; yscale = %25.15e; zscale = %25.15e;\n",
            maxOrder, xScale, yScale, zScale);

//
//  Write the C66 matrix to the output file
//
    for (i = 0; i < 6; i++) {
        fprintf(fpDeriv, "%%%%");
        for (j = 0; j < 6; j++) {
            fprintf(fpDeriv, " C%d%d = % .15e;", i, j,
                    param->elasticConstantMatrix[i][j]);
        }
        fprintf(fpDeriv, "\n");
    }

//
//  Allocate an initialize some needed arrays
//
    numYCoeff = (lMax+1)*(lMax+1);

    complex8  *Ycoeff           = (complex8 *)calloc(1, numYCoeff * sizeof(complex8));
    real8     *phi              = (real8    *)calloc(1, numPhi    * sizeof(real8   ));
    real8     *theta            = (real8    *)calloc(1, numTheta  * sizeof(real8   ));
    real8     *gaussQuadPoints  = (real8    *)calloc(1, numTheta  * sizeof(real8   ));
    real8     *gaussQuadWeights = (real8    *)calloc(1, numTheta  * sizeof(real8   ));

    for (i = 0; i < numPhi; i++) {
        phi[i] = (real8)i / numPhi * 2.0 * M_PI;
    }

    GaussQuadraturePointsAndWgts(numTheta, gaussQuadPoints, gaussQuadWeights);

//
// NOTE: gaussQuadPoints calculated above are opposite sign
//       from points returned from the source matlab code.  Values
//       adjusted below.
//
    for (i = 0; i < numTheta; i++) {
        gaussQuadPoints[i] *= -1.0;
        theta[i] = acos(gaussQuadPoints[i]);
    }

//
//  Evaluate Green's function on grid
//
    numGVals = 3 * 3 * numTheta * numPhi;
    GAll = (real8 *)calloc(1, numGVals * sizeof(complex8));

    EvalGreensFuncOnGrid(numTheta, theta, numPhi, phi, GAll,
                         numTrapezoidal, anisoVars->elasticConstantMatrix4D);

//
//  Compute spherical harmonics transform
//
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            int           L;
            real8         nrm;
            real8         GSum0 = 0.0;
            real8         GNorm0 = 0.0;
            complex8 GSum, GNorm;
            complex8 GTmp[numTheta][numPhi];

//
//  Allocate static size hash table.  If we run out of room, increase
//  hard-coded value and recompile.
//
            THash *YTbl = new THash(1000);

//
//          Make a temporary copy of the apporpriate 
//          numTheta X numPhi chunk of the G array
//
            real8  (*GPtr)[3][numTheta][numPhi]=0;   ptrcpy(GPtr,GAll);

            for (m = 0; m < numTheta; m++) {
                for (n = 0; n < numTheta; n++) {
                    GTmp[m][n] = GPtr[i][j][m][n];
                    GSum0 += creal(GTmp[m][n]) * creal(GTmp[m][n]);
                }
            }

            GNorm0 = sqrt(GSum0);

            for (m = 0; m < numYCoeff; m++) {
                Ycoeff[m] = complex8(0.0,0.0);
            }

            for (L = 0; L <= lMax; L++) {
                int   M;
                real8 plm[numTheta][lMax+1];

                for (m = 0; m < numTheta; m++) {
                    ComputeLegendre(L, gaussQuadPoints[m], &plm[m][0]);
                }

                for (M = -L; M <= L; M++) {
                    int       p;
                    real8     sFac = 1.0;
                    real8     nFac;
                    real8     plmCol[numTheta];
                    complex8  phiInt;
                    complex8 *YlmPtr;
                    complex8  Ylm[numTheta][numPhi];
                    complex8  Ylmbar[numTheta][numPhi];
                    complex8  phiM[numPhi];
                    complex8  prodYlmbarG[numTheta][numPhi];
                    complex8  sumYlmbarG[numTheta];

                    nFac = sqrt((2.0*L+1) / (4*M_PI) * Factorial(L-abs(M)) /
                                Factorial(L+abs(M)));

                    for (p = 0; p < numTheta; p++) {
                        plmCol[p] = plm[p][abs(M)];
                    }

                    for (p = 0; p < numTheta; p++) {
                        plmCol[p] *= sFac * nFac;
                    }

                    for (p = 0; p < numPhi; p++) {
                        phiM[p] = Cmplx_Exp(complex8(0,1) *
                                            complex8(M,0) * phi[p]);
                    }

//
//                  Compute numTheta X numPhi Ylm matrix
//
                    YlmPtr = &Ylm[0][0];
                    MatrixDblMultComplex(plmCol, numTheta, 1, 1,
                                         phiM, numPhi, numPhi,
                                         YlmPtr, numPhi);

                    for (n = 0; n < numTheta; n++) {
                        for (p = 0; p < numPhi; p++) {
                            Ylmbar[n][p] = Cmplx_Conj(Ylm[n][p]);
                        }
                    }

                    phiInt = 0.0;

                    for (n = 0; n < numTheta; n++) {

                        sumYlmbarG[n] = 0.0;

                        for (p = 0; p < numPhi; p++) {
                            prodYlmbarG[n][p] = Ylmbar[n][p] * GTmp[n][p];
                            sumYlmbarG[n] += prodYlmbarG[n][p];
                        }

                        phiInt += gaussQuadWeights[n] *
                                  (sumYlmbarG[n] * (2.0*M_PI/numPhi));
                    }

                    Ycoeff[L*L + L+M] = phiInt;

                    for (n = 0; n < numTheta; n++) {
                        for (p = 0; p < numPhi; p++) {
                            GTmp[n][p] = GTmp[n][p] - phiInt * Ylm[n][p];
                        }
                    }

                }  /* end for (M = -L; M <= L; M++) */

            }  /* end for (L = 0; L <= lMax; L++) */

            GSum = 0.0;

            for (m = 0; m < numTheta; m++) {
                for (n = 0; n < numPhi; n++) {
                    GSum = GTmp[m][n] * Cmplx_Conj(GTmp[m][n]);
                }
            }

            GNorm = Cmplx_Sqrt(GSum) / GNorm0;

            nrm = Cmplx_Abs(Ycoeff[0]);
            for (k = 1; k < numYCoeff; k++) {
                if ( Cmplx_Abs(Ycoeff[k]) > nrm) 
                { nrm =  Cmplx_Abs(Ycoeff[k]); }
            }

//
//          Zero out any real or imaginary components that are
//          sufficiently close to zero.
//
            for (k = 0; k < numYCoeff; k++) {
                if (abs(creal(Ycoeff[k])) < EPSILON*nrm) {
                    Ycoeff[k] = complex8(0.0, cimag(Ycoeff[k]));
                }
                if (abs(cimag(Ycoeff[k])) < EPSILON*nrm) {
                    Ycoeff[k] = complex8(creal(Ycoeff[k]), 0.0);
                }
            }

            for (L = 0; L <= lMax; L++) {
                int numNonZero = 0;
                int mVec[2*L+1];
                int jVec[2*L+1];

//
//              For the current L, from the subset of coefficients
//              we're interested in, make a list of the indices of
//              the ones which are non-zero.
//
                for (k = 0; k < 2*L+1; k++) {
                    mVec[k] = k-L;
                    jVec[k] = L*L + k;
                }

//
//              Add any non-zero coefficients to the hash table
//
                for (k = 0; k < 2*L+1; k++) {
                    if (Cmplx_Abs(Ycoeff[jVec[k]]) > 0.0) {
                        real8 a, f, realComponent, imagComponent;

//
//                      Rescale coefficients to Greengard's normalization
//
                        f = sqrt((2*L+1)/(4*M_PI));
                        realComponent = f * creal(Ycoeff[jVec[k]]);
                        imagComponent = f * cimag(Ycoeff[jVec[k]]);

//
//                      Rescale coefficients to Zlm expansion
//
                        a = Anm(L,mVec[k]);
                        realComponent = a * realComponent;
                        imagComponent = a * imagComponent;

//
//                      Full Green's function:
//
                        add(YTbl, 1, 0, L, mVec[k], realComponent, imagComponent);
                    }
                }
            }

            ComputeDerivAndEvaluate(i, j, fpDeriv, maxOrder,
                                    YTbl, xScale, yScale, zScale);

//
//          Delete the hash table for this ij set of coefficients
//
            delete YTbl;

        }  /* end for (j = 0; j < 3; j++) */

    }  /* end for (i = 0; i < 3; i++) */

    fclose(fpDeriv);

    exit(0);
    return(0);
}
#else  // ANISOTROPIC  is not defined
int main(int argc, char *argv[])
{
    printf("\nThis executable has not been compiled with ANISOTROPIC\n");
    printf("defined.  Edit makefile.setup in the main source directory,\n");
    printf("uncomment the line containing 'DEFS += -DANISOTROPIC', and\n");
    printf("recompile the code in order to build this executable.\n\n");

    exit(1);
    return(0);
}
#endif  // end ifdef ANISOTROPIC && TAYLORFMM
