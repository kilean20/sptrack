#include "nlopt_obj1.h"
//Optimizationdata constructor
Optimizationdata::Optimizationdata(string Qsetting0, string Qsetting1, string harmfile) {
	line_def(ring0);
	for(int k=0;k<ring0.Ncell;k++) {
		if(ring0.Cell[k]->NAME==string("BPM"))
			BPM_index.push_back(k);
		if(ring0.Cell[k]->TYPE==string("QUAD"))
			QUAD_index.push_back(k);
		if(ring0.Cell[k]->NAME==string("SPKICK"))
			SPKICK_index.push_back(k);
	}
	line_def(ring1);
	line_def(ring2);
	load_2Q(ring0,Qsetting0);
	load_2Q(ring1,Qsetting1);
	load_2Q(ring2,Qsetting1);
	Cal_Twiss(ring0,0.0);
	Cal_Twiss(ring1,0.0);
	Cal_Twiss(ring2,0.0);

	//harms.load(harmfile);
/******/
	ifstream myfile (harmfile.c_str());
	stringstream buffer;
	if (myfile.is_open())
	{
		buffer << myfile.rdbuf();
		myfile.close();
	}
	vector <string> fields;
	boost::split_regex( fields, buffer.str(), boost::regex( "#.*?\n" ) );
	//for (vector<string>::iterator it = fields.begin();it!=fields.end();++it)
	//	cout <<*it;
	vector<string>::iterator it = fields.begin();
	stringstream temp;

	it++;
	temp.str("");
	temp<<*it;
	urowvec temp2;
	temp2.load(temp);
	harms=temp2.t();

	it++;
	temp.str("");
	temp<<*it;
	w_harms.load(temp);

	it++;
	temp.str("");
	temp<<*it;
	w_tunes.load(temp);

	it++;
	temp.str("");
	temp<<*it;
	w_betas.load(temp);

	it++;
	temp.str("");
	temp<<*it;
	w_phis.load(temp);
/******/

	Cal_harmonics(ring0, 3.0, harms, qnx0, qnz0);
	Cal_harmonics(ring1, 3.0, harms, qnx1, qnz1);
	Cal_harmonics(ring2, 3.0, harms, qnx2, qnz2);
	bx0=zeros<vec>(BPM_index.size());
	bz0=zeros<vec>(BPM_index.size());
	bx1=zeros<vec>(BPM_index.size());
	bz1=zeros<vec>(BPM_index.size());
	bx2=zeros<vec>(BPM_index.size());
	bz2=zeros<vec>(BPM_index.size());
	phix0=zeros<vec>(BPM_index.size());
	phiz0=zeros<vec>(BPM_index.size());
	phix1=zeros<vec>(BPM_index.size());
	phiz1=zeros<vec>(BPM_index.size());
	phix2=zeros<vec>(BPM_index.size());
	phiz2=zeros<vec>(BPM_index.size());
	for (int i=0;i<BPM_index.size();i++)
	{
		bx0(i)   = ring0.Cell[BPM_index[i]]->Beta1;
		bz0(i)   = ring0.Cell[BPM_index[i]]->Beta2;
		bx1(i)   = ring1.Cell[BPM_index[i]]->Beta1;
		bz1(i)   = ring1.Cell[BPM_index[i]]->Beta2;
		bx2(i)   = ring2.Cell[BPM_index[i]]->Beta1;
		bz2(i)   = ring2.Cell[BPM_index[i]]->Beta2;
		phix0(i) = ring0.Cell[BPM_index[i]]->Mu1;
		phiz0(i) = ring0.Cell[BPM_index[i]]->Mu2;
		phix1(i) = ring1.Cell[BPM_index[i]]->Mu1;
		phiz1(i) = ring1.Cell[BPM_index[i]]->Mu2;
		phix2(i) = ring2.Cell[BPM_index[i]]->Mu1;
		phiz2(i) = ring2.Cell[BPM_index[i]]->Mu2;
	}
}
//update qnx, bx, etc...
void Optimizationdata::update_ring2 () {
	Cal_Twiss(ring2,0.0);
	Cal_harmonics(ring2, 3.0, harms, qnx2, qnz2);//To do this or not, it depends.
	for (int i=0;i<BPM_index.size();i++)
	{
		bx2(i)   = ring2.Cell[BPM_index[i]]->Beta1;
		bz2(i)   = ring2.Cell[BPM_index[i]]->Beta2;
		phix2(i) = ring2.Cell[BPM_index[i]]->Mu1;
		phiz2(i) = ring2.Cell[BPM_index[i]]->Mu2;
	}
}

