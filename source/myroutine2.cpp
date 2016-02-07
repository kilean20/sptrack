#include "myroutine2.h"

void Print_QK1L( Line linename, string quadname)
{
	for(int k=0;k<linename.Ncell;k++){
		if(linename.Cell[k]->NAME==quadname ){
			cout<<quadname<<"(K1L)="<<linename.Cell[k]->GetP("K1L")<<endl;
			break;
		}
	}
}
void Print_Quad( Line linename)
{
	vector<int> QUAD_index;
	for(int k=0;k<linename.Ncell;k++)
		if(linename.Cell[k]->TYPE==string("QUAD"))
			QUAD_index.push_back(k);
	for (vector<int>::iterator it = QUAD_index.begin(); it!=QUAD_index.end(); ++it)
		cout<<linename.Cell[*it]->NAME<<" "<<linename.Cell[*it]->GetP("K1L")<<endl;
/* check quad
	for(int k=0;k<FODO.Ncell;k++)
		if(FODO.Cell[k]->TYPE==string("QUAD"))
			cout<<FODO.Cell[k]->NAME<<" "<<FODO.Cell[k]->GetP("K1L")<<endl;
*/
}
//set two half-quadrupole strengths to K1L/2 at index_c 
void Set_2Quad ( Line & linename, uvec index_c, vec K1L)
{
	vector<int> QUAD_index;
	for(int k=0;k<linename.Ncell;k++)
		if(linename.Cell[k]->TYPE==string("QUAD"))
			QUAD_index.push_back(k);
	for (int i=0; i<index_c.n_elem; i++)
	{
		linename.Cell[QUAD_index[2*index_c(i)  ]]->SetP("K1L",K1L(i)/2);
		linename.Cell[QUAD_index[2*index_c(i)+1]]->SetP("K1L",K1L(i)/2);
	}
}
//add quadrupole strengths with trimQ at index_c 
void Set_2Quad_Add ( Line & linename, uvec index_c, vec trimQ)
{
	double QUAD_K1L;
	vector<int> QUAD_index;
	for(int k=0;k<linename.Ncell;k++)
		if(linename.Cell[k]->TYPE==string("QUAD"))
			QUAD_index.push_back(k);
	for (int i=0; i<index_c.n_elem; i++)
	{  
		//cout<<index_c(i)<<" "<<trimQ(i)<<endl;;
		QUAD_K1L = linename.Cell[QUAD_index[2*index_c(i)]]->GetP("K1L");
		QUAD_K1L += trimQ(i)/2;
		linename.Cell[QUAD_index[2*index_c(i)  ]]->SetP("K1L",QUAD_K1L);
		linename.Cell[QUAD_index[2*index_c(i)+1]]->SetP("K1L",QUAD_K1L);
	}
}

//get and sum of two half-quadrupole strengths at index_c
vector<double> Get_2Quad( Line linename, uvec index_c)
{
	vector<int> QUAD_index;
	for(int k=0;k<linename.Ncell;k++)
		if(linename.Cell[k]->TYPE==string("QUAD"))
			QUAD_index.push_back(k);
	vector<double> QUAD_K1L;
	for (int i=0; i<index_c.n_elem; i++)
	{
		QUAD_K1L.push_back(linename.Cell[QUAD_index[2*index_c(i)  ]]->GetP("K1L") 
				         + linename.Cell[QUAD_index[2*index_c(i)+1]]->GetP("K1L"));
	}
	return QUAD_K1L;
}
//put relative quadrupole errors, truncated at 3 sigma
void put_2Qerr( Line & linename, double Q_ERR)
{
	vector<int> QUAD_index;
	for(int k=0;k<linename.Ncell;k++)
		if(linename.Cell[k]->TYPE==string("QUAD"))
			QUAD_index.push_back(k);
	int Nquad=QUAD_index.size()/2;
	double QUAD_K1L,rand_factor;
	if( Q_ERR!=0 ){
		for (int i=0; i<Nquad; i++)
		{
			rand_factor = Q_ERR*trunc_gaussrand(3.0);
			QUAD_K1L = linename.Cell[QUAD_index[2*i]]->GetP("K1L");
			QUAD_K1L *= 1+rand_factor;
			linename.Cell[QUAD_index[2*i]]->SetP("K1L",QUAD_K1L);
			linename.Cell[QUAD_index[2*i+1]]->SetP("K1L",QUAD_K1L);
		}
	}
};

//save all quad strength
void save_2Q( Line linename, string filename)
{
	vector<double> QUAD_K1L;
	for(int k=0;k<linename.Ncell;k++)
		if(linename.Cell[k]->TYPE==string("QUAD"))
			QUAD_K1L.push_back(linename.Cell[k]->GetP("K1L"));

	vec out = conv_to<vec>::from(QUAD_K1L);
	out.save(filename,raw_ascii);
};

//load all quad strength
void load_2Q( Line & linename, string filename)
{
	vec in;
	in.load(filename);

	int i=0;
	for(int k=0;k<linename.Ncell;k++)
		if(linename.Cell[k]->TYPE==string("QUAD"))
		{
			linename.Cell[k]->SetP("K1L",in[i]);
			i++;
		}
};

//calculate the harmonics from the formula, perform Cal_Twiss before
void Cal_harmonics( Line linename, double ds, uvec harms, vec &qnx, vec &qnz)
{
	vector<int> SPKICK_index;
	for(int k=0;k<linename.Ncell;k++)
	{
		if(linename.Cell[k]->NAME==string("SPKICK"))
			SPKICK_index.push_back(k);
	}
	double an,bn,Xn,qn; //dummy variables
	double nux=linename.Cell[linename.Ncell-1]->Mu1;
	double nuz=linename.Cell[linename.Ncell-1]->Mu2;
	qnx=vec(harms.n_elem);
	qnz=vec(harms.n_elem);
	for (int i=0; i<harms.n_elem; i++){
		an=0;bn=0;qn=0;
		for (vector<int>::iterator it = SPKICK_index.begin(); it<SPKICK_index.end(); it++)
		{
			an+=cos(harms(i)*linename.Cell[*it]->Mu1*2*M_PI/nux);
			bn+=sin(harms(i)*linename.Cell[*it]->Mu1*2*M_PI/nux);
		}
		Xn=atan2(-bn,an);
		for (vector<int>::iterator it = SPKICK_index.begin(); it<SPKICK_index.end(); it++)
		{
			qn += cos(harms(i)*linename.Cell[*it]->Mu1*2*M_PI/nux+Xn);
		}
		qnx(i)=qn*ds*2/linename.Length;

		an=0;bn=0;qn=0;
		for (vector<int>::iterator it = SPKICK_index.begin(); it<SPKICK_index.end(); it++)
		{
			an+=cos(harms(i)*linename.Cell[*it]->Mu2*2*M_PI/nuz);
			bn+=sin(harms(i)*linename.Cell[*it]->Mu2*2*M_PI/nuz);
		}
		Xn=atan2(-bn,an);
		for (vector<int>::iterator it = SPKICK_index.begin(); it<SPKICK_index.end(); it++)
		{
			qn += cos(harms(i)*linename.Cell[*it]->Mu2*2*M_PI/nuz+Xn);
		}
		qnz(i)=qn*ds*2/linename.Length;
	}
}

//calculate sextupole driving terms
void Sext_Driving(Line & RING, double l, complex<double> &c3030l, complex<double> &c12p12l, complex<double> &c12n12l, complex<double> &c1030l, complex<double> &c1012l)
{
	//Resonance Strengths
	double nux = RING.Tune1;
	double nuz = RING.Tune2;
	complex<double> J(0,-1);
	//complex<double> c3030l, c12p12l, c12n12l, c1030l, c1012l;//from book
	c3030l=0; c12p12l=0; c12n12l=0; c1030l=0; c1012l=0;
	double betx,betz,k2l,xix,xiz,th;

	vector<int> SEXT_index;
	for(int i=0;i<RING.Ncell;i++)
		if(RING.Cell[i]->TYPE==string("SEXT") || RING.Cell[i]->TYPE==string("MULT")){
			SEXT_index.push_back(i);
		}
	vector<int>::iterator it;
	for (it = SEXT_index.begin(); it!=SEXT_index.end(); ++it) {
		betx = RING.Cell[*it]->Beta1;
		betz = RING.Cell[*it]->Beta2;
		k2l = RING.Cell[*it]->GetP("K2L");
		xix = 2*M_PI*RING.Cell[*it]->Mu1;
		xiz = 2*M_PI*RING.Cell[*it]->Mu2;
		th = RING.Cell[*it]->S/RING.Length*2*M_PI;

		c3030l  += sqrt(2)/24/M_PI  * pow(betx,1.5) * pow(betz,0.0) * k2l * exp(J*(3*xix - (3*nux-l)*th));
		c12p12l += sqrt(2)/ 8/M_PI  * pow(betx,0.5) * pow(betz,1.0) * k2l * exp(J*(xix+2*xiz - (nux+2*nuz-l)*th));
		c12n12l += sqrt(2)/ 8/M_PI  * pow(betx,0.5) * pow(betz,1.0) * k2l * exp(J*(xix-2*xiz - (nux-2*nuz-l)*th));
		c1030l  += sqrt(2)/ 8/M_PI  * pow(betx,1.5) * pow(betz,0.0) * k2l * exp(J*(xix - (nux-l)*th));
		c1012l  += sqrt(2)/ 4/M_PI  * pow(betx,0.5) * pow(betz,1.0) * k2l * exp(J*(xix - (nux-l)*th));
	}
}

//Detuning Parameters
void Cal_Detuning(Line & RING, double &a_xx, double &a_xz, double &a_zz)
{
	double nux = RING.Tune1;
	double nuz = RING.Tune2;
	a_xx=0; a_xz=0; a_zz=0;

	// Indices Implement
	vector<int> SEXT_index;
	for(int i=0;i<RING.Ncell;i++)
		if(RING.Cell[i]->TYPE==string("SEXT") || RING.Cell[i]->TYPE==string("MULT")){
			SEXT_index.push_back(i);
		}
	vector<int> OCT_index;
	for(int i=0;i<RING.Ncell;i++)
		if(RING.Cell[i]->TYPE==string("OCT") || RING.Cell[i]->TYPE==string("MULT")){
			OCT_index.push_back(i);
		}

	double betx_i, betz_i, k2l_i, xix_i, xiz_i, betx_j, betz_j, k2l_j, xix_j, xiz_j;
	//contribution from sextupoles
	vector<int>::iterator IT_i, IT_j;
	for (IT_i = SEXT_index.begin(); IT_i!=SEXT_index.end(); ++IT_i) {
		betx_i = RING.Cell[*IT_i]->Beta1;
		betz_i = RING.Cell[*IT_i]->Beta2;
		k2l_i = RING.Cell[*IT_i]->GetP("K2L");
		xix_i = 2*M_PI*RING.Cell[*IT_i]->Mu1;
		xiz_i = 2*M_PI*RING.Cell[*IT_i]->Mu2;
		for (IT_j = SEXT_index.begin(); IT_j!=SEXT_index.end(); ++IT_j) {
			betx_j = RING.Cell[*IT_j]->Beta1;
			betz_j = RING.Cell[*IT_j]->Beta2;
			k2l_j = RING.Cell[*IT_j]->GetP("K2L");
			xix_j = 2*M_PI*RING.Cell[*IT_j]->Mu1;
			xiz_j = 2*M_PI*RING.Cell[*IT_j]->Mu2;
			if(*IT_i != *IT_j)
			{
				a_xx += 1./64/M_PI * pow(betx_i,1.5) * pow(betx_j,1.5) * k2l_i * k2l_j
					* (cos(3*(M_PI*nux-abs(xix_i-xix_j)))/sin(3*M_PI*nux)
							+ 3*cos(M_PI*nux-abs(xix_i-xix_j))/sin(M_PI*nux));
				a_xz += 1./32/M_PI * pow(betx_i,0.5) * pow(betx_j,0.5) * pow(betz_i,1.0) * pow(betz_j,1.0) * k2l_i * k2l_j
					* (cos(2*(M_PI*nuz-abs(xiz_i-xiz_j))+M_PI*nux-abs(xix_i-xix_j))/sin(M_PI*(2*nuz+nux))
							+ cos(2*(M_PI*nuz-abs(xiz_i-xiz_j))-M_PI*nux+abs(xix_i-xix_j))/sin(M_PI*(2*nuz-nux)));
				a_xz -= 1./16/M_PI * pow(betx_i,1.5) * pow(betx_j,0.5) * pow(betz_j,1.0) * k2l_i * k2l_j
					* cos(M_PI*nux-abs(xix_i-xix_j))/sin(M_PI*nux);
				a_zz += 1./64/M_PI * pow(betx_i,0.5) * pow(betx_j,0.5) * pow(betz_i,1.0) * pow(betz_j,1.0) * k2l_i * k2l_j
					* (cos(2*(M_PI*nuz-abs(xiz_i-xiz_j))+M_PI*nux-abs(xix_i-xix_j))/sin(M_PI*(2*nuz+nux))
							+ cos(2*(M_PI*nuz-abs(xiz_i-xiz_j))-M_PI*nux+abs(xix_i-xix_j))/sin(M_PI*(2*nuz-nux))
							+ 3*cos(M_PI*nux-abs(xix_i-xix_j))/sin(M_PI*nux));
				// last term is inconsistent. the leading coefficient is 4 instead of 3 in simtrack
			}
		}
	}

	double betx, betz, k3l;
	vector<int>::iterator it;
	//contribution from octupoles
	for (it = OCT_index.begin(); it!=OCT_index.end(); ++it) {
		betx = RING.Cell[*it]->Beta1;
		betz = RING.Cell[*it]->Beta2;
		k3l = RING.Cell[*it]->GetP("K3L");
		a_xx += 1./16/M_PI * betx * betx * k3l;
		a_xz += -1./8/M_PI * betx * betz * k3l;
		a_zz += 1./16/M_PI * betz * betz * k3l;
	}
}

//construction of the response matrix X, including dnux and dnuz
void Cal_Response1( Line & linename, uvec harms, uvec index_c, mat &X)
{
	X=zeros<mat>(2*harms.n_elem+2,index_c.n_elem);
	vector<int> QUAD_index;
	for(int k=0;k<linename.Ncell;k++)
		if(linename.Cell[k]->TYPE==string("QUAD"))
			QUAD_index.push_back(k);
	vec qnx0=zeros<vec>(harms.n_elem);
	vec qnz0=zeros<vec>(harms.n_elem);
	double nux0=linename.Cell[linename.Ncell-1]->Mu1;
	double nuz0=linename.Cell[linename.Ncell-1]->Mu2;
	vec qnx=zeros<vec>(harms.n_elem);
	vec qnz=zeros<vec>(harms.n_elem);
	double nux, nuz;
	Cal_Twiss(linename,0.0);
	Cal_harmonics(linename, 3.0, harms, qnx0, qnz0);

	double impulse=0.0001;
	double QUAD_K1L;
	for (int i=0; i<index_c.n_elem; i++)
	{
		QUAD_K1L = linename.Cell[QUAD_index[2*index_c(i)]]->GetP("K1L");
		QUAD_K1L += impulse/2;
		linename.Cell[QUAD_index[2*index_c(i)  ]]->SetP("K1L",QUAD_K1L);
		linename.Cell[QUAD_index[2*index_c(i)+1]]->SetP("K1L",QUAD_K1L);
		Cal_Twiss(linename,0.0);
		Cal_harmonics(linename, 3.0, harms, qnx, qnz);
		nux=linename.Cell[linename.Ncell-1]->Mu1;
		nuz=linename.Cell[linename.Ncell-1]->Mu2;
		QUAD_K1L -= impulse/2;
		linename.Cell[QUAD_index[2*index_c(i)  ]]->SetP("K1L",QUAD_K1L);
		linename.Cell[QUAD_index[2*index_c(i)+1]]->SetP("K1L",QUAD_K1L);
		X(span(0,harms.n_elem-1),i)=(qnx-qnx0)/impulse;
		X(span(harms.n_elem,2*harms.n_elem-1),i)=(qnz-qnz0)/impulse;
		X(2*harms.n_elem,i)=(nux-nux0)/impulse;
		X(2*harms.n_elem+1,i)=(nuz-nuz0)/impulse;
	}
}
