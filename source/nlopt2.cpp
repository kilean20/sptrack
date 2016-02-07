/*
	last update: 05/05/2014
	re-order the argument to match with "viewer", and neater output format
	usage: ./nlopt2 label Q0.dat Q1.dat var_params.dat obj_params.dat method > logfile
    input:
		label
		Q0.dat: quadrupole setting files, no error
		Q1.dat: quadrupole setting files, before correction
		var_params.dat: quadrupole index and ranges
		obj_params.dat: parameters for objective, including weightings
		method: method
  	output:
		dqba_(label).dat: change of quadrupole strength
		harm_(label).dat: harmonics
		dnxz_(label).dat: dnu
		FODO0_(label).twiss, FODO1_(label).twiss, FODO2_(label).twiss
		Q2_(label).dat: quadrupole setting files, after correction
	todo:
		add Ctrl-C to halt gracefully
		ref: http://ab-initio.mit.edu/wiki/index.php/NLopt_Reference#Stopping_criteria
 */
#include <nlopt.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include "nlopt_obj1.h"

int iteration_count = 0;
double myobj(const vector<double> &x, vector<double> &grad, void *my_func_data);
int main(int argc, char *argv[])
{
	//load input parameters and files
	string label=argv[1];
	string Qsetting0=argv[2];
	string Qsetting1=argv[3];
	string trimquad=argv[4];
	string objparam=argv[5];
	string method=argv[6];

	//filename implement
	ofstream fout1,fout2,fout3;
	stringstream outfile1,outfile2,outfile3,outfile4,outfile5,outfile6,outfile7;
	fout1.precision(3);
	fout1.setf ( ios::scientific );
	fout2.precision(3);
	fout2.setf ( ios::scientific );
	fout3.precision(3);
	fout3.setf ( ios::scientific );
	outfile1<<"dqba_"<<label<<".dat";
	fout1.open(outfile1.str().c_str());
	outfile2<<"harm_"<<label<<".dat";
	fout2.open(outfile2.str().c_str());
	outfile3<<"dnxz_"<<label<<".dat";
	fout3.open(outfile3.str().c_str());
	outfile4<<"FODO0_"<<label<<".twiss";
	outfile5<<"FODO1_"<<label<<".twiss";
	outfile6<<"FODO2_"<<label<<".twiss";
	outfile7<<"Q2_"<<label<<".dat";

	//initialize Optimizationdata
	Optimizationdata data(Qsetting0, Qsetting1, objparam);
	cout<<"---- original: ------------------"<<endl;
	cout<<"(nux, nuz) = ("<<data.ring0.Tune1<<", "<<data.ring0.Tune2<<")"<<endl;
	cout<<"	qnx	qnz"<<endl;
	cout<<join_rows(data.qnx0,data.qnz0);
	cout<<"---- with errors: ---------------"<<endl;
	cout<<"(nux, nuz) = ("<<data.ring1.Tune1<<", "<<data.ring1.Tune2<<")"<<endl;
	cout<<"	qnx	qnz"<<endl;
	cout<<join_rows(data.qnx1,data.qnz1);


	//optimization setup
	mat trimQ;
	trimQ.load(trimquad);
	cout<<"---- quadrupole used: -----------"<<endl;
	cout<<"	index	lowerbound	upperbound"<<endl;
	cout<<trimQ;
	data.index_c=conv_to<uvec>::from(trimQ.col(0));
	vector<double> lb=conv_to< vector<double> >::from(trimQ.col(1));
	vector<double> ub=conv_to< vector<double> >::from(trimQ.col(2));

	//choose method
	nlopt::algorithm a;
	if(method=="LN_COBYLA")
	  a=nlopt::LN_COBYLA;
	else if(method=="LN_BOBYQA")
	  a=nlopt::LN_BOBYQA;
	else if(method=="LN_NEWUOA")
	  a=nlopt::LN_NEWUOA;
	else if(method=="LN_PRAXIS")
	  a=nlopt::LN_PRAXIS;
	else if(method=="LN_NELDERMEAD")
	  a=nlopt::LN_NELDERMEAD;
	else if(method=="LN_SBPLX")
	  a=nlopt::LN_SBPLX;
	else if(method=="GN_DIRECT")
	  a=nlopt::GN_DIRECT;
	else
	  a=nlopt::LN_SBPLX;

	nlopt::opt opt(a, data.index_c.n_elem);

	//set boundary and initial value
	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);
	opt.set_min_objective(myobj, &data);
	opt.set_xtol_rel(1e-4);

	double minf; 
	cout<<"---- optimizing... --------------"<<endl;
	cout<<"method: "<<method<<endl;
	vector<double> x=Get_2Quad(data.ring1,data.index_c);
	vec x1=conv_to<colvec>::from(x);
	cout<<"x1 = "<<x1.t();
/* test an objective
	vector<double> n;n.clear();
	cout<< myobj(x, n, &data);
*/
	nlopt::result result = opt.optimize(x, minf);
	vec x2=conv_to<colvec>::from(x);
	cout<<"found minimum after "<<iteration_count<<" evaluations"<<endl;
	cout<<"found minf = "<<minf<<endl;
	cout<<"at x2 = "<<x2.t();

	cout<<"---- result: --------------------"<<endl;
	cout<<"(nux, nuz) = ("<<data.ring2.Tune1<<", "<<data.ring2.Tune2<<")"<<endl;
	cout<<"	qnx	qnz"<<endl;
	cout<<join_rows(data.qnx2,data.qnz2);


	//write output files
	cout<<"---- writing files... -----------"<<endl;
	uvec index_a=linspace<uvec>(0,data.numq()-1,data.numq());
	vector<double> quada0=Get_2Quad(data.ring0,index_a); vector<double> quada1=Get_2Quad(data.ring1,index_a);
	vector<double> quada2=Get_2Quad(data.ring2,index_a);
	mat DQ=join_rows(conv_to<colvec>::from(quada1)-conv_to<colvec>::from(quada0), conv_to<colvec>::from(quada2)-conv_to<colvec>::from(quada0));
	DQ.raw_print(fout1);

	mat out2 = join_rows(data.qnx1,data.qnz1);
	out2 = join_rows(out2,data.qnx2);
	out2 = join_rows(out2,data.qnz2);
	out2.raw_print(fout2);

	mat out3;
	out3<<data.ring1.Tune1-data.ring0.Tune1<<data.ring1.Tune2-data.ring0.Tune2
		<<data.ring2.Tune1-data.ring0.Tune1<<data.ring2.Tune2-data.ring0.Tune2;
	out3.raw_print(fout3);

    Print_Twiss(data.ring0,outfile4.str().c_str());
    Print_Twiss(data.ring1,outfile5.str().c_str());
    Print_Twiss(data.ring2,outfile6.str().c_str());
	save_2Q(data.ring2,outfile7.str());
	//save_2Q(data.ring2,"Q2_"+label+".dat");

	cout<<outfile1.str()<<", ";
	cout<<outfile2.str()<<", ";
	cout<<outfile3.str()<<", ";
	cout<<outfile4.str()<<", ";
	cout<<outfile5.str()<<", ";
	cout<<outfile6.str()<<", ";
	cout<<outfile7.str()<<" written!"<<endl;

	return 0;
}   

//extern int iteration_count;
//double Optimizationdata::myobj1(const vector<double> &x, vector<double> &grad, void *my_func_data)
double myobj(const vector<double> &x, vector<double> &grad, void *my_func_data)
{
	++iteration_count;
	//Optimizationdata mydata = *((Optimizationdata *)my_func_data);
	Optimizationdata * mydata = (Optimizationdata *)my_func_data;
	//Line ring = mydata.ring2;

	Set_2Quad((*mydata).ring2, (*mydata).index_c, conv_to<colvec>::from(x));
	(*mydata).update_ring2();

    if (!grad.empty()) {
        grad[0] = 0.0;
        grad[1] = 0.0;
    }

	vector<double> o1;
	o1.push_back(pow((*mydata).ring2.Tune1-(*mydata).ring0.Tune1,2));
	o1.push_back(pow((*mydata).ring2.Tune2-(*mydata).ring0.Tune2,2));
	vec o=conv_to<vec>::from(o1);
	o=join_vert(o,pow((*mydata).qnx2,2));
	o=join_vert(o,pow((*mydata).qnz2,2));
	o=join_vert(o,pow((*mydata).bx2/(*mydata).bx0-1,2));
	o=join_vert(o,pow((*mydata).bz2/(*mydata).bz0-1,2));
	o=join_vert(o,pow((*mydata).phix2-(*mydata).phix0,2));
	o=join_vert(o,pow((*mydata).phiz2-(*mydata).phiz0,2));

	rowvec w=join_horiz((*mydata).w_tunes,(*mydata).w_harms.row(0));
	w=join_horiz(w,(*mydata).w_harms.row(1));
	w=join_horiz(w,zeros<rowvec>((*mydata).numb()*2));
	w=join_horiz(w,zeros<rowvec>((*mydata).numb()*2));

	//return inner_product(weight.begin(),weight.end(),o.begin(),0.);
	double out = as_scalar(w*o);
	if(iteration_count%10==0) cout<<iteration_count<<" iterations, obj="<<out<<endl;
	return out;
}

