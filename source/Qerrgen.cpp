//date:02/17
//usage: ./qerrgen label nux nuz Qerr scale
//output: Q_label.dat, FODO_label.twiss

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "elepass.h"
#include "physics.h"
#include "FODO4.h"
#include "myroutine.h"
#include "myroutine2.h"
#include "armadillo"

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
	//input parameters
	unsigned int seed;
	stringstream(argv[2]) >> seed;
	double nux;
	stringstream(argv[3]) >> nux;
	double nuz;
	stringstream(argv[4]) >> nuz;
	double Q_ERR;
	stringstream(argv[5]) >> Q_ERR; //relative Q-errors
	double scale;
	stringstream(argv[6]) >> scale;

	srand(seed);//random seed
	//Define the ring
	Line FODO;
	line_def(FODO);
	Fit_Tune(FODO, nux, nuz, "QF", "QD"); 

	uvec index_c = conv_to<uvec>::from(linspace(0,35,36));
	vector<double> Q0 = Get_2Quad( FODO, index_c);
	put_2Qerr( FODO, Q_ERR);
	vector<double> Q1 = Get_2Quad( FODO, index_c);
	vector<double> Q2;
	for(int i=0;i<Q1.size();i++)
		Q2.push_back(Q0[i] + scale*(Q1[i]-Q0[i]));
	Set_2Quad ( FODO, index_c, conv_to<vec>::from(Q2));

	Cal_Twiss(FODO,0.0);//it is necessary for the 0th turn
	cout<<"label "<<argv[1]<<":(nux, nuz) = ("<<FODO.Tune1<<", "<<FODO.Tune2<<")"<<endl;
	stringstream outfile1;
	outfile1<<"FODO_"<<argv[1]<<".twiss";
    Print_Twiss(FODO,outfile1.str().c_str());

	save_2Q(FODO,"Q_"+string(argv[1])+".dat");
	return 0;
}

