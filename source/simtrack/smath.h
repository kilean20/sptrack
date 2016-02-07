#ifndef SMATH_H
#define SMATH_H
#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cdf.h>
//===========================================
//
//          math functions
//
//===========================================

long int  fac(int n);
//---random number generator
double rnd(const double & r);
double guass(double u,double g, double & r);
//---polynominal fitting
void pfit(double x[],double y[], int n, double a[], int m);
//---matrix kits
void mat_mult(double A[],  double B[], double C[], int m, int n, int k);
double mat_det(double a[],int n);
int mat_inv(double a[],int n);
//----eigen tune solver
void mat_change_hessenberg(double a[], int n);
int mat_root_hessenberg(double a[],int n,double u[],double v[],double eps,int jt);
void LinearEquations(double a1, double b1, double c1, double a2, double b2, double c2, double & x, double & y  );
// --the equations to be solved
//     a1*x +b1*y= c1
//     a2*x +b2*y= c2


//---eigen solver from GSL
void  EigenSolver(double Matrix[4][4] , double wr[4], double  wi[4], double vr[4][4], double vi[4][4]) ;
void  EigenSolver_6D(double Matrix[6][6] , double wr[6], double  wi[6], double vr[6][6], double vi[6][6]) ;

//---SVD stuffs
void brmul(double a[], double b[] , int m, int n, int k, double c[] );
static void ppp(double *a,double *e,double *s,double *v,int m,int n);
static void sss(double fg[2],double cs[2]);
int bmuav(double a[],int m,int n,double *u,double *v,double eps,int ka);
//------FFT, tune finders
void fft(int m, double*x, double*y);
void  Sin2FFT(int n, double*xr, double*xi, int & nfft, double*famp, double*fphi);
void  FindMaxPeak(int n, double*xr, double*xi, int & nfft, double*famp, double*fphi, double& peakf, double& peakamp, double& peakphi, int flag);
void FineTuneFinder(int nfft, double *x, double & fpeak);

#endif
