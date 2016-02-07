#include "interface.h"

//=====================================
//
//     interfaces to MADX
//
//=====================================
void  Read_MADXLattice(const char * filename, Line & linename)
{
  int i;
  fstream f1;
  string line;
  string temp_name;
  string temp_type;
  double temp_s;
  double temp_l;
  double temp_angle;
  double temp_e1;
  double temp_e2;
  double temp_tilt;
  double knl[11], knsl[11];
  double temp_sk;
  Element * temp_element;

  f1.open(filename,ios::in);
  if(!f1)
    {
      cout<<"error in opening the file: "<<filename<<endl;
      exit(0);
    }
  for(i=1;i<=47;i++) getline(f1,line,'\n');
  
  double s1=0;
  while( getline(f1,line,'\n')){
    istringstream ss(line);
    ss>>temp_name>>temp_type>>temp_s>>temp_l>>temp_angle>>temp_e1>>temp_e2>>temp_tilt
      >>knl[0]>>knsl[0]
      >>knl[1]>>knsl[1]
      >>knl[2]>>knsl[2]
      >>knl[3]>>knsl[3]
      >>knl[4]>>knsl[4]
      >>knl[5]>>knsl[5]
      >>knl[6]>>knsl[6]
      >>knl[7]>>knsl[7]
      >>knl[8]>>knsl[8]
      >>knl[9]>>knsl[9]
      >>knl[10]>>knsl[10]
      >>temp_sk;
    temp_name=string(temp_name,1,temp_name.size()-2 );     
    temp_type=string(temp_type,1,temp_type.size()-2 );
    
    s1=s1+temp_l;
    if(temp_type==string("DRIFT")){
      temp_element=new DRIFT(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("SBEND"))    {  // MADX also supports sbends with K1 and K2
      if(knl[1] == 0 and  knl[2] == 0  ) {
	temp_element= new SBEND(temp_name,temp_l,temp_angle, temp_e1, temp_e2);
	linename.Cell.push_back(temp_element); }
      else {
	temp_element=new GMULT(temp_name,temp_l, temp_angle, temp_e1, temp_e2, knl, knsl);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("RBEND"))    {
      if(knl[1] == 0 ) {
	temp_element= new SBEND(temp_name,temp_l,temp_angle, temp_e1, temp_e2);
	linename.Cell.push_back(temp_element); }
      else {
	temp_element=new GMULT(temp_name,temp_l, temp_angle, temp_e1, temp_e2, knl, knsl);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("QUADRUPOLE")){  // MAD8 quads only have K1, MADX also has K1S
      if ( knsl[1] != 0 ) {
	cout<<" Alert: QUAD with k1s detected. "<<endl;
      }
      temp_element= new QUAD(temp_name,temp_l, knl[1]);
      linename.Cell.push_back(temp_element); 
    }
    else if (temp_type==string("SEXTUPOLE")){   // MAD8 sexts only have K2, MADX also has K2S
      if ( knsl[2] != 0 ) {
	cout<<" Alert: SEXT with kss detected. "<<endl;
      }
      temp_element= new SEXT(temp_name,temp_l, knl[2]);
      linename.Cell.push_back(temp_element); 
    }
    else if (temp_type==string("OCTUPOLE")){    // MAD8 octs only have K3, MADX also has K2S
      if ( knsl[3] != 0 ) {
	cout<<" Alert: OCT with kss detected. "<<endl;
      }
      temp_element= new OCT(temp_name,temp_l, knl[3]);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("HKICKER")){
      temp_element= new HKICKER(temp_name,temp_l, -knl[0] );
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("VKICKER")){
      temp_element= new VKICKER(temp_name,temp_l, knsl[0]);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("KICKER")){
      temp_element= new KICKER(temp_name,temp_l, -knl[0], knsl[0]);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("MULTIPOLE")){  // MADX only support zero-length multipoles
      int m;
      double kl=0.;
      for(m=0;m<11;m++)  kl=kl+abs(knl[m])+abs(knsl[m]);
      if(kl==0.){
	temp_element=new DRIFT(temp_name,temp_l);
	linename.Cell.push_back(temp_element);	
      }
      else{
	temp_element=new MULT(temp_name,temp_l, knl, knsl);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("HMONITOR")){
      temp_element= new HBPM(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("VMONITOR")){
      temp_element= new VBPM(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("MONITOR")){
      temp_element= new BPM(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("MARKER")){
      temp_element= new MARKER(temp_name,temp_l);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("SOLENOID")){  
      if(temp_sk !=0.0  and temp_l !=0. ){    // meaning  input for solenoid 
	temp_element= new SOLEN(temp_name,temp_l, temp_sk);
	linename.Cell.push_back(temp_element);
      }
      else{
	temp_element=new DRIFT(temp_name,temp_l);
	linename.Cell.push_back(temp_element);
      }
    }
    else if (temp_type==string("RFCAVITY")){
      temp_element= new RFCAV(temp_name,temp_l, 0.0, 78250.42279*GP.harm, 0.);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("BEAMBEAM")){
      if(temp_l !=0.){ 
	cout<<"Beam-Beam element should have zero length!"<<endl;
        exit(1);  
      }
      temp_element= new BEAMBEAM(temp_name,6, 0., 0.4545, 11, 2.5e-06, 2.5e-06, 0.53, 0.0, 0.53, 0.0);
      linename.Cell.push_back(temp_element);
    }
    else if (temp_type==string("ELENS")){
      temp_element= new ELENS(temp_name, temp_l, 0.0, 5, 0.14, .31e-03, .31e-03);
      linename.Cell.push_back(temp_element);
    }
    else
      {
	temp_element=new DRIFT(temp_name,temp_l);
	linename.Cell.push_back(temp_element);
      }
  }
  linename.Ncell=linename.Cell.size();
  linename.Update();
  f1.close();
}

void  Print_MADXLattice(const char * filename, Line & linename)
{
  int i, j, flag;
  fstream f2;
  char str[125];
  vector <string> namelist;
  
  f2.open(filename,ios::out);
  if(!f2)
    {
      cout<<"error in opening the file: "<<filename<<endl;
      exit(0);
    }
  
  for(i=0;i<linename.Ncell; i++){
    flag=0;
    for(j=0;j<namelist.size();j++) {
      if(linename.Cell[i]->NAME  == namelist [j] ) {
	flag=1;break;
      }
    }
    if(flag==0) {
      namelist.push_back(linename.Cell[i]->NAME);
      if(linename.Cell[i]->TYPE==string("DRIFT") ){
      	f2<<linename.Cell[i]->NAME<<" : DRIFT, L = "<<setprecision(10)<<linename.Cell[i]->L<<";"<<endl; }
      else if(linename.Cell[i]->TYPE==string("QUAD") )  {
	f2<<linename.Cell[i]->NAME<<" : QUADRUPOLE, L = "<<setprecision(10)<<linename.Cell[i]->L<<", K1 =  "<<setprecision(10)<<linename.Cell[i]->GetP("K1L")/ linename.Cell[i]->L<<" ;"<<endl;}
      else if(linename.Cell[i]->TYPE==string("SEXT") ) {
	f2<<linename.Cell[i]->NAME<<" : SEXTUPOLE,  L = "<<setprecision(10)<<linename.Cell[i]->L<<", K2 =  "<<setprecision(10)<<linename.Cell[i]->GetP("K2L")/ linename.Cell[i]->L<<" ;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("SBEND") ) { 
	f2<<linename.Cell[i]->NAME<<" : SBEND,      L = "<<setprecision(10)<<linename.Cell[i]->L<<", ANGLE =  "<<setprecision(10)<<linename.Cell[i]->GetP("ANGLE")<<" ;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("RFCAV") ) {  
	f2<<linename.Cell[i]->NAME<<" : RFCAVITY,   L = "<<setprecision(10)<<linename.Cell[i]->L<<" ;"<<endl; }
      else if(linename.Cell[i]->TYPE==string("MULT") ) {
	f2<<linename.Cell[i]->NAME<<" : MULTIPOLE,  "<<endl;
	f2<<" knl:={ "<<linename.Cell[i]->GetP("K0L")<<","<<linename.Cell[i]->GetP("K1L")<<","<<linename.Cell[i]->GetP("K2L")<<","<<linename.Cell[i]->GetP("K3L")<<","<<linename.Cell[i]->GetP("K4L")<<","<<linename.Cell[i]->GetP("K5L")<<",";  
	f2<<linename.Cell[i]->GetP("K6L")<<","<<linename.Cell[i]->GetP("K7L")<<","<<linename.Cell[i]->GetP("K8L")<<","<<linename.Cell[i]->GetP("K9L")<<","<<linename.Cell[i]->GetP("K10L")<<"},"<<endl;
	f2<<" ksl:={ "<<linename.Cell[i]->GetP("K0SL")<<","<<linename.Cell[i]->GetP("K1SL")<<","<<linename.Cell[i]->GetP("K2SL")<<","<<linename.Cell[i]->GetP("K3SL")<<","<<linename.Cell[i]->GetP("K4SL")<<","<<linename.Cell[i]->GetP("K5SL")<<",";  
	f2<<linename.Cell[i]->GetP("K6SL")<<","<<linename.Cell[i]->GetP("K7SL")<<","<<linename.Cell[i]->GetP("K8SL")<<","<<linename.Cell[i]->GetP("K9SL")<<","<<linename.Cell[i]->GetP("K10SL")<<"};"<<endl; 
        if (linename.Cell[i]->L != 0. ) cout<<"Warning: "<< linename.Cell[i]->NAME<<"  "<<linename.Cell[i]->TYPE<<"   "<<linename.Cell[i]->L<<"  transferred to DRIFT ."<<endl;   }
      else {
	f2<<linename.Cell[i]->NAME<<" : DRIFT, L = "<<setprecision(10)<<linename.Cell[i]->L<<";"<<endl;
        cout<<"Warning: "<< linename.Cell[i]->NAME<<"  "<<linename.Cell[i]->TYPE<<"   "<<setprecision(10)<<linename.Cell[i]->L<<"  transferred to DRIFT ."<<endl;
      }
    }
  }
  
  int istart=0;
  int inumber10=0; 
  int nline=linename.Ncell / 1000;

  for(i=0;i<nline+1;i++) {
    sprintf(str,"LIN%d",i);
    f2<<str<<" : LINE=("<<endl;
    inumber10=0;

    do{
      f2<<linename.Cell[istart]->NAME<<",";
      istart++;
      inumber10++;
      if (inumber10 ==10){ f2<< endl; inumber10=0;}
    } while ( istart < 1000*(i+1)  &&  istart < linename.Ncell-1  );
    
    if(i < nline ){
      f2<<linename.Cell[istart]->NAME<<" ); "<<endl;
      istart++;
    }
    else{
      f2<<linename.Cell[linename.Ncell-1]->NAME<<" ); "<<endl;
    }
  }

  f2<<"rhic:  LINE = ( ";
  for(i=0;i<nline;i++) {
    sprintf(str,"LIN%d",i);
    f2<<str<<",";
  }
  sprintf(str,"LIN%d",i);
  f2<<str<<" );"<<endl;

  f2<<"beam, mass:=0.93827, charge:=1, gamma:=268.2, exn:=20.0e-06, eyn:=20.0e-06, sige:=0.001;"<<endl;
  f2<<"use, period=rhic;"<<endl;
  f2<<"select, flag=twiss,  column=NAME, KEYWORD,S,BETX, BETY, DX, DY, MUX, MUY;"<<endl;
  f2<<"twiss,table=twiss,file=twiss.table;"<<endl;
  f2<<"stop;"<<endl;
  f2.close();
}

