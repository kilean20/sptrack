#include "myroutine.h"
#include <complex>
#include "simtrack/elepass.h"
#include "simtrack/physics.h"
#include "armadillo"

using namespace arma;
using namespace std;

void Print_QK1L( Line linename, string quadname);
void Print_Quad( Line linename);
void Set_2Quad ( Line & linename, uvec index_c, vec K1L);
void Set_2Quad_Add ( Line & linename, uvec index_c, vec trimQ);
vector<double> Get_2Quad( Line linename, uvec index_c);
void put_2Qerr( Line & linename, double Q_ERR);
void load_2Q( Line & linename, string filename);
void save_2Q( Line linename, string filename);

void Cal_harmonics( Line linename, double ds, uvec harms, vec &qnx, vec &qnz);
void Sext_Driving(Line & RING, double l, complex<double> &c3030l, complex<double> &c12p12l, complex<double> &c12n12l, complex<double> &c1030l, complex<double> &c1012l);
void Cal_Detuning(Line & RING, double &a_xx, double &a_xz, double &a_zz);
void Cal_Response1( Line & linename, uvec harms, uvec index_c, mat &X);
