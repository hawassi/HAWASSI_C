#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdbool.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <fftw3.h>
#include <omp.h>

typedef struct{
	double Lleft;	
	double Lright;
	double Lbottom;
	double Ltop;
	double** charac;
}fbdyvar;

typedef struct{
	char Id[100];
	char Filename[100];
	double** data;
	double depth;
	char type[1];
}bathyvar;

typedef struct{
	double* x;
	double* y;
	double* t;
	int Nx;
	int Ny;
	int Nt;
	double dx;
	double dy;
	double dt;
	int tcoarse;
	double Nadj;
	double* kx;
	double* ky;
	double** KK;
	double** aal;
	bathyvar bathy;
	fbdyvar  fbdy;
	double 	 U;
	double** wallchar;
}domvar;

typedef struct{
	double Re;
	double Im;
}complex;

domvar dominfo;
char* dir_force;
complex** couple_etahat;
complex** couple_phihat;
double 	  alpha_couple;
double**  Chi_couple_H;
int 	  iter_t;

#include "functions1.c"

extern void HAWASSI_coupling_term(char* ,double ,double ,double** ,double** ,
int ,int ,int );
