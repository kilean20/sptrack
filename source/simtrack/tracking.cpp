#include "tracking.h"
//==============================================================
//
//                Tracking 
//
//===============================================================

void Track(Line & linename, double x[6], int nturn, int & stable, int & lost_turn, int & lost_post)
{
  int j;
  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;    lost_turn= 0;    lost_post = 0;    return; }
  
  //----now we do tracking
  for(GP.turn=0; GP.turn < nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy  ) { 
      	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      }
      //each element to output
      // cout<<j<<"   ";
      //for(k=0;k<6;k++) cout<<x[k]<<" ";
      //cout<<endl;
    } 
    // each turn to output
    // for(k=0;k<6;k++) cout<<x[k]<<" ";
    // cout<<endl;
  }
}

void Track_tbt(Line & linename, double x[6], int nturn, double  x_tbt[], int & stable, int & lost_turn, int & lost_post)
{
  int j;

  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;     lost_turn= 0;     lost_post = 0;     return;   }

  //----now we track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    for(j=0;j<6;j++) x_tbt[GP.turn*6+j]=x[j];
    
    //for(j=0;j<6;j++) cout<<x[j]<<"  ";
    //cout<<endl;
    
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy  ) { 
	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      } } }
}

void Track_tbt(Line & linename, double x[6], int nturn, double  x_tbt[], int &bpm_index, int & stable, int & lost_turn, int & lost_post)
{
  int j,k;

  //-----quick check 
  if( abs(x[0]) > 0.1 ||  abs(x[2]) > 0.1 || stable ==0 ) {
    stable = 0;     lost_turn= 0;     lost_post = 0;     return;   }

  //----now we track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    for(j=0;j<linename.Ncell;j++) {
      linename.Cell[j]->Pass(x);
      if(j == bpm_index ){
	for(k=0;k<6;k++) x_tbt[GP.turn*6+k]=x[k];
	//for(k=0;k<6;k++) cout<<x[k]<<"  ";
	//cout<<endl;
      }
      if( abs(x[0]) > linename.Cell[j]->APx or abs(x[2]) > linename.Cell[j]->APy  ) { 
	stable=0;  lost_turn= GP.turn;  lost_post=j; return ;
      } } }
}


//---containers for elements for Track_Fast()
int     nelement;
int     Ndrift, Nsbend, Nquad,  Nsext,     Nmult;
int     Nrf,    Nbb,    Nelens, Ncooling,  Nacmult, Nmatrix;
int     Nkicker, Nhkicker, Nvkicker, Ndiffuse, Ngmult, Nlrbb;

vector <int>     type;
vector <double>  nelement_dx, nelement_dy, nelement_dt, nelement_APx, nelement_APy;
vector <double>  drift_l;
vector <int>     sbend_nint;
vector <double>  sbend_l, sbend_angle, sbend_e1, sbend_e2;
vector <int>     quad_nint;
vector <double>  quad_l,  quad_k1l;
vector <int>     sext_nint;
vector <double>  sext_l,  sext_k2l;
vector <int>     mult_nint, mult_norder;
vector <double>  mult_l, mult_knl,  mult_knsl;
vector <double>  bb_np,  bb_emitx,  bb_emity, bb_sigmal, bb_betax, bb_alfax, bb_betay, bb_alfay;
vector <int>     bb_treat, bb_nslice;
vector <int>     elens_nslice;
vector <double>  elens_l, elens_ne, elens_betae, elens_sigmax, elens_sigmay;
vector <double>  rf_l, rf_vrf,  rf_frf,  rf_phi0;
vector <double>  cooling_l, cooling_alpha;
vector <int>     acmult_norder, acmult_tturns;
vector <double>  acmult_l, acmult_kl, acmult_phi0;
vector <double>  matrix_l, matrix_m66, matrix_xco_in, matrix_xco_out;
vector <double>  kick_l, kicker_hkick, kicker_vkick;
vector <double>  hkick_l, hkicker_hkick;
vector <double>  vkick_l, vkicker_vkick;
vector <double>  diffuse_l, diffuse_diff_x, diffuse_diff_y, diffuse_diff_delta;
vector <int>     gmult_nint, gmult_norder;
vector <double>  gmult_l,gmult_angle, gmult_knl,  gmult_knsl, gmult_e1, gmult_e2;
vector <double>  lrbb_np, lrbb_sepx, lrbb_sepy, lrbb_sigmax, lrbb_sigmay, lrbb_kickx0, lrbb_kicky0;

//----------------------------------------------------------------
//Tracking without  linename.Cell[i]->Pass() to speed up tracking
//----------------------------------------------------------------
void Prepare_Track_Fast(Line & linename)
//  this function only accepts the allowed elements and change others to drift
//  this function accepts thick lenses: SBEND, QUAD, SEXT, MULT, ELENS, GMULT
{
  int i,j,k;
  //-----------do some statistics
  nelement=linename.Ncell;
  Ndrift=0;
  Nsbend=0;
  Nquad=0;
  Nsext=0;
  Nmult=0;
  Nrf=0;
  Nbb=0;
  Nelens=0;
  Ncooling=0;
  Nacmult =0;
  Nmatrix =0;
  Nkicker =0;
  Nhkicker=0;
  Nvkicker=0;
  Ndiffuse=0;
  Ngmult  =0;
  Nlrbb   =0;

  for(i=0; i<linename.Ncell; i++){
    if(linename.Cell[i]->TYPE==string("DRIFT") ) {
      Ndrift++;
    }
    else if(linename.Cell[i]->TYPE==string("SBEND") ) {
      Nsbend++;
    }
    else if(linename.Cell[i]->TYPE==string("QUAD") ){
      Nquad++;
    }
    else if(linename.Cell[i]->TYPE==string("SEXT") ){
      Nsext++;
    }
    else if(linename.Cell[i]->TYPE==string("MULT") ){
      Nmult++;
    }
    else if(linename.Cell[i]->TYPE==string("BEAMBEAM") ){
      Nbb++;
    }
    else if(linename.Cell[i]->TYPE==string("ELENS") ){
      Nelens++;
    }
    else if(linename.Cell[i]->TYPE==string("RFCAV") ){
      Nrf++;
    }
    else if(linename.Cell[i]->TYPE==string("COOLING") ){
      Ncooling++;
    }
    else if(linename.Cell[i]->TYPE==string("ACMULT") ){
      Nacmult++;
    }
    else if(linename.Cell[i]->TYPE==string("MATRIX") ){
      Nmatrix++;
    }
    else if(linename.Cell[i]->TYPE==string("KICKER") ){
      Nkicker++;
    }
    else if(linename.Cell[i]->TYPE==string("HKICKER") ){
      Nhkicker++;
    }
    else if(linename.Cell[i]->TYPE==string("VKICKER") ){
      Nvkicker++;
    }
    else if(linename.Cell[i]->TYPE==string("DIFFUSE") ){
      Ndiffuse++;
    }
    else if(linename.Cell[i]->TYPE==string("GMULT") ){
      Ngmult++;
    }
    else if(linename.Cell[i]->TYPE==string("LRBB") ){
      Nlrbb++;
    }
    else{ 
      Ndrift++;
    }
  }
  
  if( false ) {
    cout<<"---------------------------------"<<endl;
    cout<<"Statistics before Track_Fast(): "<<endl;
    cout<<"There are "<<Ndrift<< " DRIFT"<<endl;
    cout<<"There are "<<Nsbend<< " SEBDN"<<endl;
    cout<<"There are "<<Nquad<< "  QUAD"<<endl;
    cout<<"There are "<<Nsext<< "  SEXT"<<endl;
    cout<<"There are "<<Nmult<< "  MULT"<<endl;
    cout<<"There are "<<Nbb<< "  BEAMBEAM"<<endl;
    cout<<"There are "<<Nelens<< "  ELENS"<<endl;
    cout<<"There are "<<Nrf<< "  RFCAV"<<endl;
    cout<<"There are "<<Ncooling<< "  COOLING"<<endl;
    cout<<"There are "<<Nacmult<< "  ACMULT"<<endl;
    cout<<"There are "<<Nmatrix<< "  MATRIX"<<endl;
    cout<<"There are "<<Nkicker<< "  KICKER"<<endl;
    cout<<"There are "<<Nhkicker<< "  HKICKER"<<endl;
    cout<<"There are "<<Nvkicker<< "  VKICKER"<<endl;
    cout<<"There are "<<Ndiffuse<< "  DIFFUSE"<<endl;
    cout<<"There are "<<Ngmult<< "  GMULT"<<endl;
    cout<<"There are "<<Nlrbb<< "  LRBB"<<endl;
    cout<<"Totally "<<linename.Ncell<<"  elements"<<endl;
  }
 
  //--- store strengths
  for(i=0; i<nelement; i++){
    
    nelement_dx.push_back(linename.Cell[i]->DX );
    nelement_dy.push_back(linename.Cell[i]->DY );
    nelement_dt.push_back(linename.Cell[i]->DT );
    nelement_APx.push_back(linename.Cell[i]->APx );
    nelement_APy.push_back(linename.Cell[i]->APy ); 

    if(linename.Cell[i]->TYPE==string("DRIFT") ) {
      type.push_back( 0 );
      drift_l.push_back(linename.Cell[i]->L);
    }
    else if(linename.Cell[i]->TYPE==string("SBEND") ) {
      type.push_back( 1  );
      sbend_l.push_back( linename.Cell[i]->L );
      sbend_angle.push_back( linename.Cell[i]->GetP("ANGLE") );
      sbend_nint.push_back( int(linename.Cell[i]->GetP("Nint")) );
      sbend_e1.push_back( linename.Cell[i]->GetP("E1") );
      sbend_e2.push_back( linename.Cell[i]->GetP("E2") );
    }
    else if(linename.Cell[i]->TYPE==string("QUAD") ){
      type.push_back( 2 );
      quad_l.push_back( linename.Cell[i]->L );
      quad_k1l.push_back( linename.Cell[i]->GetP("K1L") );
      quad_nint.push_back( int(linename.Cell[i]->GetP("Nint")) );
    }
    else if(linename.Cell[i]->TYPE==string("SEXT") ){
      type.push_back( 3 );
      sext_l.push_back(linename.Cell[i]->L);
      sext_nint.push_back( int(linename.Cell[i]->GetP("Nint")));
      sext_k2l.push_back(linename.Cell[i]->GetP("K2L"));
    }
    else if(linename.Cell[i]->TYPE==string("MULT") ){
      type.push_back( 4 );
      mult_l.push_back(linename.Cell[i]->L);
      mult_nint.push_back(int(linename.Cell[i]->GetP("Nint")));
      mult_norder.push_back( int( linename.Cell[i]->GetP("Norder") )  );
      for(j=0;j<11;j++){
        char name1[125], name2[125];
        sprintf(name1, "K%dL",j);
        sprintf(name2, "K%dSL",j);
        mult_knl.push_back(linename.Cell[i]->GetP(name1));
        mult_knsl.push_back(linename.Cell[i]->GetP(name2));
      }
    }
    else if(linename.Cell[i]->TYPE==string("BEAMBEAM") ){
      type.push_back( 5 );
      bb_np.push_back( linename.Cell[i]->GetP("NP") );
      bb_treat.push_back( int( linename.Cell[i]->GetP("TREATMENT") ) );
      bb_nslice.push_back( int( linename.Cell[i]->GetP("NSLICE") ) );
      bb_emitx.push_back( linename.Cell[i]->GetP("EMITX") );
      bb_emity.push_back( linename.Cell[i]->GetP("EMITY") );
      bb_sigmal.push_back( linename.Cell[i]->GetP("SIGMAL") );
      bb_betax.push_back( linename.Cell[i]->GetP("BETAX") );
      bb_alfax.push_back( linename.Cell[i]->GetP("ALFAX") );
      bb_betay.push_back( linename.Cell[i]->GetP("BETAY") );
      bb_alfay.push_back( linename.Cell[i]->GetP("ALFAY") );
    }
    else if(linename.Cell[i]->TYPE==string("ELENS") ){
      type.push_back( 6  );
      elens_ne.push_back( linename.Cell[i]->GetP("NE") );
      elens_l.push_back( linename.Cell[i]->L );     
      elens_nslice.push_back( int( linename.Cell[i]->GetP("NSLICE") ) );
      elens_sigmax.push_back( linename.Cell[i]->GetP("SIGMAX") );
      elens_sigmay.push_back( linename.Cell[i]->GetP("SIGMAY") );
      elens_betae.push_back( linename.Cell[i]->GetP("BETAE") );
    }
    else if(linename.Cell[i]->TYPE==string("RFCAV") ){
      type.push_back( 7 );
      rf_l.push_back( linename.Cell[i]->L);
      rf_vrf.push_back( linename.Cell[i]->GetP("VRF") );
      rf_frf.push_back( linename.Cell[i]->GetP("FRF") );
      rf_phi0.push_back( linename.Cell[i]->GetP("PHASE0") );
    }
    else if(linename.Cell[i]->TYPE==string("COOLING") ){
      type.push_back( 8 );
      cooling_l.push_back( linename.Cell[i]->L );
      cooling_alpha.push_back( linename.Cell[i]->GetP("ALPHA") );
    }
    else if(linename.Cell[i]->TYPE==string("ACMULT") ){
      type.push_back( 9 );
      acmult_l.push_back( linename.Cell[i]->L );     
      acmult_norder.push_back( int(linename.Cell[i]->GetP("Norder") ));
      acmult_tturns.push_back( int(linename.Cell[i]->GetP("Tturns") ));
      acmult_kl.push_back( linename.Cell[i]->GetP("KL") );
      acmult_phi0.push_back( linename.Cell[i]->GetP("PHI0") );
    }
   else if(linename.Cell[i]->TYPE==string("MATRIX") ){
      type.push_back( 10 );
      matrix_l.push_back( linename.Cell[i]->L); 
      for(j=0;j<6;j++)
	for(k=0;k<6;k++)
	  {
	    char name1[125];
	    sprintf(name1, "M%d%d",j+1,k+1);
	    matrix_m66.push_back( linename.Cell[i]->GetP(name1) ); 
	  }
      
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_X") );
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_PX") );
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_Y") );
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_PY") );
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_Z") );
      matrix_xco_in.push_back(  linename.Cell[i]->GetP("XCO_IN_DELTA") );

      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_X") );
      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_PX") );
      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_Y") );
      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_PY") );
      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_Z") );
      matrix_xco_out.push_back(  linename.Cell[i]->GetP("XCO_OUT_DELTA") );
    }
    else if(linename.Cell[i]->TYPE==string("KICKER") ){
      type.push_back( 11 );
      kicker_hkick.push_back( linename.Cell[i]->GetP("HKICK") );
      kicker_vkick.push_back( linename.Cell[i]->GetP("VKICK") );
    }
    else if(linename.Cell[i]->TYPE==string("HKICKER") ){
      type.push_back( 12 );
      hkicker_hkick.push_back( linename.Cell[i]->GetP("HKICK") );
    }
    else if(linename.Cell[i]->TYPE==string("VKICKER") ){
      type.push_back( 13 );
      vkicker_vkick.push_back( linename.Cell[i]->GetP("VKICK") );
    }
    else if(linename.Cell[i]->TYPE==string("DIFFUSE") ){
      type.push_back( 14 );
      diffuse_l.push_back( linename.Cell[i]->L);
      diffuse_diff_x.push_back( linename.Cell[i]->GetP("DIFF_X") );
      diffuse_diff_y.push_back( linename.Cell[i]->GetP("DIFF_Y") );
      diffuse_diff_delta.push_back( linename.Cell[i]->GetP("DIFF_DELTA") );
    }
    else if(linename.Cell[i]->TYPE==string("GMULT") ){
      type.push_back( 15 );
      gmult_l.push_back(linename.Cell[i]->L); 
      gmult_nint.push_back( int(linename.Cell[i]->GetP("Nint")) ); 
      gmult_norder.push_back( int(linename.Cell[i]->GetP("Norder")) ); 
      gmult_l.push_back( linename.Cell[i]->L );
      gmult_angle.push_back( linename.Cell[i]->GetP("ANGLE") );
      gmult_e1.push_back( linename.Cell[i]->GetP("E1") );
      gmult_e2.push_back( linename.Cell[i]->GetP("E2") );
      for(j=0;j<11;j++){
        char name1[125], name2[125];
        sprintf(name1, "K%dL",j);
        sprintf(name2, "K%dSL",j);
        gmult_knl.push_back(linename.Cell[i]->GetP(name1));
        gmult_knsl.push_back(linename.Cell[i]->GetP(name2));
      }
    }
    else if(linename.Cell[i]->TYPE==string("LRBB") ){
      type.push_back( 16 );
      lrbb_np.push_back( linename.Cell[i]->GetP("NP") );
      lrbb_sepx.push_back( linename.Cell[i]->GetP("SEPX") );
      lrbb_sepy.push_back( linename.Cell[i]->GetP("SEPY") );
      lrbb_sigmax.push_back( linename.Cell[i]->GetP("SIGMAX") );
      lrbb_sigmay.push_back( linename.Cell[i]->GetP("SIGMAY") );
      lrbb_kickx0.push_back( linename.Cell[i]->GetP("KICKX0") );
      lrbb_kicky0.push_back( linename.Cell[i]->GetP("KICKY0") );
    }
    else{ 
      //-----I hope no more Warning here. Should be prevented by Make_Thin()
      cout<<" Warning: "<<linename.Cell[i]->NAME<<" :  "<<linename.Cell[i]->TYPE<<"  at  "<<linename.Cell[i]->S<<"  converted to DRIFT."<<endl; 
      type.push_back( 0 );
      double temp=linename.Cell[i]->L;
      drift_l.push_back( temp  );
    }
  }
}

void Track_Fast(Line & linename, double x[6], int nturn, int & stable, int & lost_turn, int & lost_post, double & sum_x2, double & sum_y2, double & sum_z2)
{
  int i,j,k;
  double  Dpx, Dpy;
  double  knl[11], knsl[11];
  double  m66[36], xco_in[6], xco_out[6];
  double  x1[6];

  sum_x2=0.;  sum_y2=0.;  sum_z2=0.;

  //-----quick check before track
  if( abs(x[0]) > 1.0 ||  abs(x[2]) > 1.0 || stable ==0 ) {
    stable = 0;      lost_turn= 0;       lost_post = 0;      return ;   }
  
  //---now we do track
  for(GP.turn=0;GP.turn<nturn; GP.turn++) {
    
    Ndrift=0;
    Nsbend=0;
    Nquad=0;
    Nsext=0;
    Nmult=0;
    Nrf=0;
    Nbb=0;
    Nelens=0;
    Ncooling=0;
    Nacmult=0;
    Nmatrix=0;   
    Nkicker =0;
    Nhkicker=0;
    Nvkicker=0;
    Ndiffuse=0;
    Ngmult=0;
    Nlrbb=0;
    
    for(i=0; i<nelement; i++){
      
      //--can be used for each element or some type of element
      //GtoL(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
      
      if(type[i]==0 ) {
	DRIFT_Pass(x, drift_l[Ndrift]);
	Ndrift++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==1 ) {
	SBEND_Pass(x,sbend_l[Nsbend],sbend_nint[Nsbend], sbend_angle[Nsbend],sbend_e1[Nsbend],sbend_e2[Nsbend]); 
	Nsbend++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==2 ){
	QUAD_Pass(x,  quad_l[Nquad], quad_nint[Nquad], quad_k1l[Nquad], 0.); 
	Nquad++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==3 ){
	SEXT_Pass(x, sext_l[Nsext], sext_nint[Nsext], sext_k2l[Nsext], 0.);
	Nsext++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==4 ){
	for(j=0;j<11;j++) {
	  knl[j]=mult_knl[Nmult*11+j]; 
	  knsl[j]=mult_knsl[Nmult*11+j];
	}
	MULT_Pass(x, mult_l[Nmult], mult_nint[Nmult], mult_norder[Nmult], knl, knsl);
	Nmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==5 ){
	BEAMBEAM_Pass(x,bb_treat[Nbb], bb_np[Nbb], bb_sigmal[Nbb], bb_nslice[Nbb], bb_emitx[Nbb], 
		      bb_betax[Nbb], bb_alfax[Nbb], bb_emity[Nbb], bb_betay[Nbb], bb_alfay[Nbb]);
	Nbb++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==6 ){
	GtoL(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
        ELENS_Pass(x, elens_ne[Nelens], elens_l[Nelens], elens_betae[Nelens],
		   elens_nslice[Nelens], elens_sigmax[Nelens], elens_sigmay[Nelens]);
	LtoG(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 
	Nelens++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if (type[i]==7 ){
        //RFCAV_Pass(x, rf_l[Nrf], rf_vrf[Nrf], rf_frf[Nrf], 0.);
	x[5] =x[5] + (rf_vrf[Nrf]*1.0/GP.energy)*sin(2.0*M_PI*rf_frf[Nrf]*x[z_]/3.0e8);
	Nrf++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if (type[i]==8 ){
        //COOLING_Pass(x, cooling_L[Ncooling], cooling_alpha[Ncooling]);
	for(j=0;j<4;j++) x[j]= (1.- cooling_alpha[Ncooling] )* x[j];
	Ncooling++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if (type[i]==9 ){
	ACMULT_Pass(x,  acmult_l[Nacmult], acmult_norder[Nacmult], acmult_kl[Nacmult], acmult_tturns[Nacmult], acmult_phi0[Nacmult]);
	Nacmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==10 ){
	for(j=0;j<36;j++) {
	  m66[j]=matrix_m66[Nmatrix*36+j]; 
	}
	for(j=0;j<6;j++) xco_in[j] = matrix_xco_in[Nmatrix*6+j];
	for(j=0;j<6;j++) xco_out[j]= matrix_xco_out[Nmatrix*6+j];
	
	for(j=0;j<6;j++) x[j]=x[j]-xco_in[j];
	for(j=0;j<6;j++)    x1[j]=x[j];
	for(j=0;j<6;j++) {
	  x[j]=0.0;
	  for(k=0;k<6;k++)  x[j] = x[j] + m66[j*6+k] * x1[k];
	}  
	for(j=0;j<6;j++) x[j]=x[j]+xco_out[j];
	
	Nmatrix++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==11 ){
	x[1]=x[1]+kicker_hkick[Nkicker];   
	x[3]=x[3]+kicker_vkick[Nkicker];
	Nkicker++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==12 ){
	x[1]=x[1]+hkicker_hkick[Nhkicker];   
	Nhkicker++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==13 ){
	x[1]=x[1]+vkicker_vkick[Nvkicker];   
	Nvkicker++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==14 ){
	x[1]= x[1]+ diffuse_diff_x[Ndiffuse] * rnd(seed); 
	x[3]= x[1]+ diffuse_diff_y[Ndiffuse] * rnd(seed); 
	x[5]= x[1]+ diffuse_diff_delta[Ndiffuse] * rnd(seed); 
	Ndiffuse++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==15 ){
	for(j=0;j<11;j++) {
	  knl[j] =gmult_knl[Ngmult*11+j]; 
	  knsl[j]=gmult_knsl[Ngmult*11+j];
	}
	GMULT_Pass(x, gmult_l[Ngmult], gmult_nint[Ngmult], gmult_angle[Ngmult], gmult_norder[Ngmult],knl,knsl,
                      gmult_e1[Ngmult],gmult_e2[Ngmult]);
	Ngmult++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else if(type[i]==16 ){
        LRBB_Pass(x, GP.gamma, lrbb_np[Nlrbb],lrbb_sepx[Nlrbb],lrbb_sepy[Nlrbb],lrbb_sigmax[Nlrbb],lrbb_sigmay[Nlrbb],lrbb_kickx0[Nlrbb],lrbb_kicky0[Nlrbb] );
	Nlrbb++;
	if( abs(x[0]) > nelement_APx[i] || abs(x[2]) > nelement_APy[i] ) {
	  stable = 0;
	  lost_turn= GP.turn;
	  lost_post = i;
	  return;
	}
      }
      else{ }

     //LtoG(x, nelement_dx[i],nelement_dy[i],nelement_dt[i]); 

       //print out each element 
      // cout<<i<<"   ";
      // for(j=0;j<6;j++) cout<<x[j]<<" ";
      // cout<<endl;
    }
    sum_x2 += x[0]*x[0];
    sum_y2 += x[2]*x[2]; 
    sum_z2 += x[4]*x[4]; 

    // print out turn-by-turn data
    // for(j=0;j<6;j++) cout<<x[j]<<" ";
    // cout<<endl;
  }
}

void Track_DA( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0 ) // Track with ->Pass()
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma1, sigma2, sigma_max=100.0;  
  int    stable, stable1=1, stable2=1,  lost_turn=0, lost_post=0; 
  FILE *f2;
  f2=fopen("./DA-output.dat","w");
  fclose(f2);
  double x[6];

  for(i=0; i<10;i++){
    
    angle=( 90/11.)*(i+1)* M_PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      sigma1= sigma;
      sigma2= sigma + 0.05;

      x[0]= sigma1 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma1 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable1=1; lost_turn=0; lost_post=0;
      Track(linename, x, nturn, stable1, lost_turn, lost_post);
      
      x[0]= sigma2 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma2 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable2=1; lost_turn=0; lost_post=0;
      Track(linename, x, nturn, stable2, lost_turn, lost_post);
      
      if ( stable1* stable2 ==1 ) {
	stable =1 ; }
      else{
	stable=0;
      }
      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1  );

    f2=fopen("DA-output.dat" ,"a");
    fprintf(f2,"%f %f \n", ( 90/11.)*(i+1),  sigma);
    fclose(f2);
  }
}
 

void Track_Fast_DA( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0 ) // Track without class
{
  int i;
  double angle, sigma_s, sigma_e, sigma, sigma1, sigma2, sigma_max=100.0;  
  int    stable, stable1=1, stable2=1,  lost_turn=0, lost_post=0; 
  double sum_x2, sum_y2, sum_z2;
  FILE *f2;
  f2=fopen("./DA-output.dat","w");
  fclose(f2);
  double x[6];

  for(i=0; i<10;i++){

    angle=( 90/11.)*(i+1)* M_PI/180;
    sigma_s=0.;
    sigma_e=sigma_max;
    
    do {
      sigma = (sigma_s + sigma_e )/2.0;
      sigma1= sigma;
      sigma2= sigma + 0.05;

      x[0]= sigma1 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma1 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable1=1; lost_turn=0; lost_post=0;
      Track_Fast(linename, x, nturn, stable1, lost_turn, lost_post, sum_x2, sum_y2, sum_z2);
      
      x[0]= sigma2 * cos(angle) * sigmax0;
      x[1]= 0. ;
      x[2]= sigma2 * sin(angle) * sigmay0;
      x[3]= 0. ;
      x[4]= 0.;
      x[5]= deltap0;
      stable2=1; lost_turn=0; lost_post=0;
      Track_Fast(linename, x, nturn, stable2, lost_turn, lost_post, sum_x2, sum_y2, sum_z2);
      
      if ( stable1* stable2 ==1 ) {
	stable =1 ; }
      else{
	stable=0;
      }
      if( stable ==0 ) 	sigma_e= sigma;
      if( stable ==1 ) 	sigma_s= sigma;
    } while ( abs( sigma_e- sigma_s) > 0.1  );

    f2=fopen("DA-output.dat" ,"a");
    fprintf(f2,"%f %f \n", ( 90/11.)*(i+1),  sigma);
    fclose(f2);
  }
} 

