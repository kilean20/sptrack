#include <cmath>
#include <vector>
#include "simtrack/elepass.h"

#include <cstdlib> //for gaussian

//void Cal_Sext_Driving(Line & RING, double l);

double unifRand();
double unifRand(double, double);
double exponentialrand();
double trunc_exponentialrand(double);
double gaussrand();
double trunc_gaussrand(double);
double xgaussRand();
double trunc_xgaussRand(double);
//void KV_dist4D(double &x1, double &x2, double &x3, double &x4, double a, double b);
void KV_dist4D(double &x1, double &x2, double &x3, double &x4);
void KV_dist4D_sibuya(double &x1, double &x2, double &x3, double &x4);
int factorial(int n);
