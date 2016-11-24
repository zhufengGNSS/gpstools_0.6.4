/*------------------------------------------------------------------------------
% [system] : GpsTools
% [module] : common routines/definitions
% [func]   : common routines/definitions
% [argin]  : 
% [argout] :
% [note]   : compile options
%            -DUNIX : UNIX/LINUX
% [version]: $Revision: 6 $ $Date: 06/07/08 1:01 $
%            Copyright(c) 2004-2006 by T.Takasu, all rights reserved
% [history]: 04/01/05   0.1  new
*-----------------------------------------------------------------------------*/
#ifndef QTCMN_H
#define QTCMN_H
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* constants -----------------------------------------------------------------*/
#ifndef M_PI
#define M_PI		3.1415926535897932				/* pi */
#endif
#define DEG2RAD		(M_PI/180.0)					/* deg->rad */
#define RAD2DEG		(180.0/M_PI)					/* rad->deg */
#define SEC2RAD		4.8481368110953598E-6			/* "->rad */

#define	M_C			299792458.0						/* speed of light(m/sec) */
#define GME			3.986004415E+14					/* geogravity(JGM-3) */
#define GMS			1.32712440017987E+20			/* solar-gravity(DE405) */
#define GMM			4.902801E+12					/* lunar-gravity(DE405) */

#ifdef UNIX
#define FILESEP     "/"								/* file path separator */
#else
#define FILESEP     "\\"
#endif

/* basic functions/operations ------------------------------------------------*/
#define NI(x)		(sizeof(x)/sizeof(*x))			/* array item counts */
#define MOD(x,y)	((x)-(y)*floor((x)/(y)))		/* modulo */
#define ABS(x)		((x)<0?-(x):(x))				/* absolute value */
#define MAX(x,y)	((x)>(y)?(x):(y))				/* maximum value */
#define MIN(x,y)	((x)<(y)?(x):(y))				/* minimum value */
#define SWAP(x,y)	do {double _t; _t=y; y=x; x=_t;} while (0) /* swap values */

#define COPY(X,n,m,Y) memcpy(Y,(void *)X,(n)*(m)*sizeof(double))

/* memory allocation/deallocation --------------------------------------------*/
#define VEC(n)		((double *)malloc((n)*    sizeof(double)))
#define MAT(n,m)	((double *)malloc((n)*(m)*sizeof(double)))
#define ZVEC(n)		((double *)calloc((n),    sizeof(double)))
#define ZMAT(n,m)	((double *)calloc((n)*(m),sizeof(double)))
#define FreeMat(mat) free(mat)

/* dot product of 2d/3d vectors ----------------------------------------------*/
#define DOT2(x,y)	((x)[0]*(y)[0]+(x)[1]*(y)[1])
#define DOT(x,y)	((x)[0]*(y)[0]+(x)[1]*(y)[1]+(x)[2]*(y)[2])

/* cross product of 3d vectors -----------------------------------------------*/
#define CROSS(x,y,z) do { \
	(z)[0]=(x)[1]*(y)[2]-(x)[2]*(y)[1]; \
	(z)[1]=(x)[2]*(y)[0]-(x)[0]*(y)[2]; \
	(z)[2]=(x)[0]*(y)[1]-(x)[1]*(y)[0]; \
} while (0)

/* norm of 2d/3d vector ------------------------------------------------------*/
#define NORM2(x)	(sqrt(DOT2(x,x)))
#define NORM(x)		(sqrt(DOT(x,x)))

/* normalize 3d vector -------------------------------------------------------*/
#define NORMV(x) do { \
	double _xx=NORM(x); (x)[0]/=_xx; (x)[1]/=_xx; (x)[2]/=_xx; \
} while (0)

/* zero 3d matrix ------------------------------------------------------------*/
#define ZERO(X) do { \
	(X)[0]=(X)[1]=(X)[2]=(X)[3]=(X)[4]=(X)[5]=(X)[6]=(X)[7]=(X)[8]=0.0; \
} while (0)

/* unit 3d matrix ------------------------------------------------------------*/
#define EYE(X) do { \
	(X)[1]=(X)[2]=(X)[3]=(X)[5]=(X)[6]=(X)[7]=0.0;(X)[0]=(X)[4]=(X)[8]=1.0; \
} while (0)

/* transpose 3d matrix -------------------------------------------------------*/
#define Tr(X,Y) do { \
	(Y)[0]=(X)[0];(Y)[1]=(X)[3];(Y)[2]=(X)[6];(Y)[3]=(X)[1];(Y)[4]=(X)[4]; \
	(Y)[5]=(X)[7];(Y)[6]=(X)[2];(Y)[7]=(X)[5];(Y)[8]=(X)[8]; \
} while (0)

/* product of 3d matrix and vector -------------------------------------------*/
#define Mv(X,y,z) do { \
	(z)[0]=(X)[0]*(y)[0]+(X)[3]*(y)[1]+(X)[6]*(y)[2]; \
	(z)[1]=(X)[1]*(y)[0]+(X)[4]*(y)[1]+(X)[7]*(y)[2]; \
	(z)[2]=(X)[2]*(y)[0]+(X)[5]*(y)[1]+(X)[8]*(y)[2]; \
} while (0)

/* product of 3d matrixes ----------------------------------------------------*/
#define MM(X,Y,Z) do { \
	Mv(X,Y,Z); Mv(X,Y+3,Z+3); Mv(X,Y+6,Z+6); \
} while (0)

/* 3d coordinate rotation matrix ---------------------------------------------*/
#define Rx(t,X) do { \
	(X)[0]=1.0; (X)[1]=(X)[2]=(X)[3]=(X)[6]=0.0; \
	(X)[4]=(X)[8]=cos(t); (X)[7]=sin(t); (X)[5]=-(X)[7]; \
} while (0)

#define Ry(t,X) do { \
	(X)[4]=1.0; (X)[1]=(X)[3]=(X)[5]=(X)[7]=0.0; \
	(X)[0]=(X)[8]=cos(t); (X)[2]=sin(t); (X)[6]=-(X)[2]; \
} while (0)

#define Rz(t,X) do { \
	(X)[8]=1.0; (X)[2]=(X)[5]=(X)[6]=(X)[7]=0.0; \
	(X)[0]=(X)[4]=cos(t); (X)[3]=sin(t); (X)[1]=-(X)[3]; \
} while (0)

/* product of matrixes -------------------------------------------------------*/
extern void MatMul(const double *X, const double *Y, int m, int k, int n,
                   double *Z);

/* inverse matrix ------------------------------------------------------------*/
extern int MatInv(const double *X, int n, double *Y);

/* transpose matrix ----------------------------------------------------------*/
extern void MatTr(const double *X, int m, int n, double *Y);

#endif /* QTCMN_H */
