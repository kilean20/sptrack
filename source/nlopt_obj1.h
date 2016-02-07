#include "elepass.h"
#include "physics.h"
#include "FODO4.h"
#include "armadillo"
#include "myroutine2.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <boost/regex.hpp>
#include <boost/algorithm/string/regex.hpp>

using namespace std;
using namespace arma;

class Optimizationdata {
	vector<int> BPM_index;
	vector<int> QUAD_index;
	vector<int> SPKICK_index; //no use yet
	public:
	uvec index_c;

	Line ring0;
	Line ring1;
	Line ring2;
	vec qnx0, qnz0, qnx1, qnz1, qnx2, qnz2;
	vec bx0, bz0, bx1, bz1, bx2, bz2;
	vec phix0, phiz0, phix1, phiz1, phix2, phiz2;

	uvec harms;
	mat w_harms;
	mat w_tunes;
	vec w_betas;
	vec w_phis;

	//Constructor;
	Optimizationdata (string Qsetting0, string Qsetting1, string harmfile);

	void update_ring2();
	int numq() {return QUAD_index.size()/2;}
	int numb() {return BPM_index.size();}
	int get_index_QUAD(int i) {return QUAD_index[i];}
};
