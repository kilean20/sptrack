#include "myroutine.h"

// Generate a random number between 0 and 1
double unifRand()
{
	return rand() / double(RAND_MAX);
}
// Generate a random number in a real interval.
double unifRand(double a, double b)
{
	return (b-a)*unifRand() + a;
}

double trunc_xgaussRand(double x)
{
	double r;
	do {
		r = xgaussRand();
	} while( r > x); //within x=3.0 sigma is 98.9%
	return r;
}
double xgaussRand() // by inverse cdf
{
	double r=unifRand();
	return sqrt(-2*log(r));//same as sqrt(-2*log(1-r))
}

double trunc_exponentialrand(double x)
{
	double r;
	do {
		r = exponentialrand();
	} while( r > x); //within x=3 sigma is 95%
	return r;
}
double exponentialrand() // by inverse cdf, mean=1
{
	double r=unifRand();
	return -log(r);//same as -log(1-r)
}

double trunc_gaussrand(double x) // mean=0, sig!=1
{
	double r;
	do {
		r = gaussrand();
	} while(fabs(r) > x); //within +- 3 sigma is 99.7% (68-95-99.7 empirical rule)
	return r;
}

double gaussrand() // by Knuth Sec. 3.4.1, mean=0, sig=1
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}

void KV_dist4D(double &x1, double &x2, double &x3, double &x4)
{
	double th = unifRand(0, M_PI/2);
	double phi1 = unifRand(0, 2*M_PI);
	double phi2 = unifRand(0, 2*M_PI);
	x1 = sin(th)*sin(phi1);
	x2 = sin(th)*cos(phi1);
	x3 = cos(th)*sin(phi2);
	x4 = cos(th)*cos(phi2);
}

void KV_dist4D_sibuya(double &x1, double &x2, double &x3, double &x4)
{
	double u1 = unifRand();
	double phi1 = unifRand(0, 2*M_PI);
	double phi2 = unifRand(0, 2*M_PI);
	x1 = sqrt(u1)*sin(phi1);
	x2 = sqrt(u1)*cos(phi1);
	x3 = sqrt(1-u1)*sin(phi2);
	x4 = sqrt(1-u1)*cos(phi2);
}

int factorial(int n)
{
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

