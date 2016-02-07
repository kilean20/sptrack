#include "physics.h"
//================================================================
//
//    PHYSICS  calculation
//
//================================================================

//---5D orbit and Twiss 
void Cal_Orbit_Num(Line & linename, double deltap)
{
  int i,j,k,iter=0;
  double x0[6],x01[6], x1[6], dx[4];
  double d=1.0e-09;
  double mat[4][4];
  double chi;
  double codeps=1e-12;
  int Max_iter=20;

  for (i=0;i<6;i++) x0[i]=0.;
  do{
    
    iter++;
    x0[4]=0;
    x0[5]=deltap;
    
    for(j=0;j<6;j++) x01[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x01);
    
    
    for(k=0;k<4;k++){
      for(j=0;j<6;j++) x1[j]=x0[j];
      x1[k]=x1[k]+d;
      for(j=0;j<linename.Ncell;j++) linename.Cell[j]->Pass(x1);
      for(j=0;j<4;j++) mat[j][k]=(x1[j]-x01[j])/d;
    }
    
    for(j=0;j<6;j++) x1[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x1);
    for(i=0;i<4;i++) {
      dx[i]=x0[i]-x1[i];
    }
    chi=0;
    for(i=0;i<4;i++) chi+=dx[i]*dx[i];
    chi=sqrt(chi/4.);
    
    if( chi > codeps ){
      for(i=0;i<4;i++) mat[i][i]=mat[i][i]-1.0000001;
      mat_inv(&mat[0][0],4); 
      for(i=0;i<4;i++){
	for(j=0;j<4;j++) x0[i]+=mat[i][j]*dx[j];  
      } 
    } 
    
  }while (chi>codeps && iter <Max_iter );
  
  if(iter == Max_iter-1 ) 
    {
      cout<<"Failed to find COD."<<endl;
      exit(1);
    }
  else
    {
      x0[4]=0.000;
      x0[5]=deltap;
      for(j=0;j<linename.Ncell; j++) {
	linename.Cell[j]->Pass(x0);
	for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
      } 
    }
}

void Cal_Orbit_DA(Line & linename, double deltap)
{
  int i,j,iter=0;
  double x0[6], x1[6], dx[6];
  tps tps1[6];
  linmap m1;
  double mat[6][6];
  double chi;
  double codeps=1e-12;
  int Max_iter=20;
  
  for (i=0;i<6;i++) x0[i]=0.;

  do{
  
    iter++;

    x0[4]=0;
    x0[5]=deltap;
    m1.identity();
    m1=m1+x0;
    
    for(i=0;i<6;i++) tps1[i]=m1[i];
    for(j=0;j<linename.Ncell; j++) linename.Cell[j]->DAPass(tps1);
    for(i=0;i<6;i++) m1[i]=tps1[i];
    Getmat(m1, &mat[0][0]);
    Getpos(m1, x1);
    
    for(i=0;i<6;i++) {
      dx[i]=x0[i]-x1[i];
    }

    chi=0;
    for(i=0;i<4;i++) chi+=dx[i]*dx[i];
    chi=sqrt(chi/4.);
 
    if( chi > codeps ){
      dx[4]=0.;
      dx[5]=0.;
      for(i=0;i<6;i++) mat[i][i]=mat[i][i]-1.0000001;
      
      mat_inv(&mat[0][0],6); 
      for(i=0;i<6;i++){
	for(j=0;j<6;j++) x0[i]+=mat[i][j]*dx[j];  
      } 
    } 
    
  }while (chi>codeps && iter <Max_iter );
  
  if(iter == Max_iter-1 ) 
    {
      cout<<"Failed to find COD."<<endl;
      exit(1);
    }
  else
    {
      x0[4]=0.000;
      x0[5]=deltap;
      for(j=0;j<linename.Ncell; j++) {
	linename.Cell[j]->Pass(x0);
	for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
      } 
    }
}

void Cal_OneTurnMap(Line & linename, double deltap)
{
  int i,j;
  double x[6];
  tps tps1[6];
  linmap m1;
  double m66[36], b66[36];
  int flag;
  double u[6], v[6];

  for(i=0;i<6;i++) x[i]=linename.Cell[linename.Ncell-1]->X[i];
  x[4]=0.000;
  x[5]=deltap; 
  m1.identity();
  m1=m1+x;
  
  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(j=0;j<linename.Ncell; j++) linename.Cell[j]->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, m66); 

  for(i=0;i<36;i++) linename.Cell[linename.Ncell-1]->M[i]=m66[i];
  
  if( false ) {
    cout<<"One turn map M:"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<m66[i*6+j]<<"  ";
      cout<<endl;
    }
    for(i=0;i<36;i++) b66[i]=m66[i];
    cout<<"Det of M = "<<mat_det(b66,6)<<endl;
  }
  
  mat_change_hessenberg(m66, 6);
  flag=mat_root_hessenberg(m66,6,u,v,1.0e-10,60);
  if(flag < 0){
    cout<<"Unstable motion. "<<endl;
    exit(1);
  }
}

void Cal_ElementMap(Line & linename, double deltap)
{
  int i,j,k, istart;
  double x[6];
  tps tps1[6];
  linmap m1;
  double t66[36];

  for(k=0; k<linename.Ncell; k++) 
    {
      istart=k-1;
      if (istart <0) istart=linename.Ncell-1;

      for(i=0;i<6;i++) x[i]=linename.Cell[istart]->X[i];
      x[4]=0.000;
      x[5]=deltap; 
      m1.identity();
      m1=m1+x;
      
      for(i=0;i<6;i++) tps1[i]=m1[i];
      linename.Cell[k]->DAPass(tps1);
      for(i=0;i<6;i++) m1[i]=tps1[i];
      Getmat(m1, t66);
      
      for(j=0;j<36;j++) linename.Cell[k]->T[j]=t66[j]; 
    }
}

void Cal_SectionMap(Line & linename, int i1, int i2, double deltap, double t66[36] )
{
  int i;
  double x[6];
  tps tps1[6];
  linmap m1;
  
  if (i1== 0 ) { 
    for(i=0;i<6;i++) x[i] = linename.Cell[linename.Ncell-1]->X[i]; }
  else {
    for(i=0;i<6;i++) x[i] = linename.Cell[i1-1]->X[i]; }
  x[4]=0.000;
  x[5]=deltap; 
  m1.identity();
  m1=m1+x;
  
  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(i=i1;i<=i2;i++)  linename.Cell[i]->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, t66);
}

//------original 4-D approach
void Cal_A0(Line & linename, double deltap)
{
  int i,j;
  double m66[6][6], m44[4][4];

  double wr[4], wi[4], vr[4][4], vi[4][4], temp_wr[4], temp_wi[4], temp_vr[4][4], temp_vi[4][4];  
  double temp1, temp2;
  int flag;
  double theta1, theta2, tempx, tempy, x1,y1;
  double a44[4][4], b44[4][4];
  double det, scale;
  
  //---solve eigen problem to construct A[4*4]
  for(i=0;i<6;i++)
   for(j=0;j<6;j++) m66[i][j]=linename.Cell[linename.Ncell-1]->M[i*6+j];

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      m44[i][j]=m66[i][j];
  EigenSolver(m44, wr, wi, vr, vi);
  if(false){
    cout<<"Eigen values and eigen vectors: "<<endl;
    for(i=0;i<4;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<4;j++){
        cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl; 
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

  //---link v and v* together
  temp1=wr[0];
  temp2=wi[0];
  for(i=0; i<4;i++)     {
    if (wr[i]==temp1 && wi[i]== -temp2 ) { 
      temp_wr[0]= wr[0];
      temp_wi[0]= wi[0];
      for(j=0;j<4;j++) temp_vr[0][j]= vr[0][j];
      for(j=0;j<4;j++) temp_vi[0][j]= vi[0][j];
      temp_wr[1]= wr[i];
      temp_wi[1]= wi[i];       
      for(j=0;j<4;j++) temp_vr[1][j]= vr[1][j];
      for(j=0;j<4;j++) temp_vi[1][j]= vi[1][j];
    }
  }
  
  flag=0;
  for(i=0; i<4;i++)     {
    if (wr[i] != temp1 ) { 
      if(flag ==0 ){
	temp_wr[2]= wr[i];
	temp_wi[2]= wi[i];
	for(j=0;j<4;j++) temp_vr[2][j]= vr[i][j];
	for(j=0;j<4;j++) temp_vi[2][j]= vi[i][j];
	flag=1;
      }       
      if(flag==1){
	temp_wr[3]= wr[i];
	temp_wi[3]= wi[i];       
	for(j=0;j<4;j++) temp_vr[3][j]= vr[i][j];
	for(j=0;j<4;j++) temp_vi[3][j]= vi[i][j];
      }
    }
  }

  //---link eigen vetcors to  H/V planes 
  temp1=abs( temp_vr[0][0] ) + abs( temp_vi[0][0] ) + abs( temp_vr[0][1] ) + abs( temp_vi[0][1] );
  temp2=abs( temp_vr[2][0] ) + abs( temp_vi[2][0] ) + abs( temp_vr[2][1] ) + abs( temp_vi[2][1] );
  
  if (temp1 > temp2 ) {
    for(i=0;i<4;i++){
      wr[i]= temp_wr[i];
      wi[i]= temp_wi[i];       
      for(j=0;j<4;j++) vr[i][j]= temp_vr[i][j];
      for(j=0;j<4;j++) vi[i][j]= temp_vi[i][j];        
    }
  }
  else
    {
      for(i=0;i<4;i++){
	wr[3-i]= temp_wr[i];
	wi[3-i]= temp_wi[i];       
	for(j=0;j<4;j++) vr[3-i][j]= temp_vr[i][j];
	for(j=0;j<4;j++) vi[3-i][j]= temp_vi[i][j];        
      }
    }

  //----rotating to make A12=0, A34=0
  theta1=atan2(vi[0][0], vr[0][0]);
  x1= cos(-M_PI/2-theta1);
  y1=sin(-M_PI/2-theta1);
  for(i=0;i<4;i++){
    tempx=vr[0][i]*x1-vi[0][i]*y1;
    tempy=vi[0][i]*x1+vr[0][i]*y1;
    vr[0][i]=tempx;
    vi[0][i]=tempy; }
  x1= cos(M_PI/2+theta1);
  y1=sin(M_PI/2+theta1);
  for(i=0;i<4;i++){
    tempx=vr[1][i]*x1-vi[1][i]*y1;
    tempy=vi[1][i]*x1+vr[1][i]*y1;
    vr[1][i]=tempx;
    vi[1][i]=tempy;}

  theta2=atan2(vi[2][2], vr[2][2]);
  x1= cos(-M_PI/2-theta2);
  y1=sin(-M_PI/2-theta2);
  for(i=0;i<4;i++){
    tempx=vr[2][i]*x1-vi[2][i]*y1;
    tempy=vi[2][i]*x1+vr[2][i]*y1;
    vr[2][i]=tempx;
    vi[2][i]=tempy; }
   x1  = cos(M_PI/2+theta2);
   y1=sin(M_PI/2+theta2);
  for(i=0;i<4;i++){
    tempx=vr[3][i]*x1-vi[3][i]*y1;
    tempy=vi[3][i]*x1+vr[3][i]*y1;
    vr[3][i]=tempx;
    vi[3][i]=tempy;}

  //----produce A
  for(i=0;i<4;i++) a44[i][0]=-(vi[0][i]-vi[1][i]);
  for(i=0;i<4;i++) a44[i][1]=vr[0][i]+vr[1][i];
  for(i=0;i<4;i++) a44[i][2]=-(vi[2][i]-vi[3][i]);
  for(i=0;i<4;i++) a44[i][3]=vr[2][i]+vr[3][i];

  if(false) {
    cout<<"A matrix :"<<endl;
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) cout<<setw(12)<<a44[i][j]<<"  ";
      cout<<endl;
    }
  }

  if(a44[1][1] < 0. ) {
    //cout<<"... Eigenmode I: reverted sign. "<<endl;
    for(i=0;i<4;i++) a44[i][1]=- a44[i][1];
  }
  if(a44[3][3] < 0. ) {
    //cout<<"... Eigenmode II: reverted sign. "<<endl;
    for(i=0;i<4;i++) a44[i][3]=- a44[i][3];
  }
  
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) b44[i][j]=a44[i][j];
  }
  det= mat_det( &b44[0][0], 4);

  scale=sqrt(sqrt(abs(det)));
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) a44[i][j]=a44[i][j]/scale;
  }
  
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) 
      linename.Cell[linename.Ncell-1]->A[i*6+j]=a44[i][j];
}

//-----adapted from 6-D general approach
void Cal_A(Line & linename, double deltap)
{
  
  int i,j;
  double m44[4][4];

  double wr[4], wi[4], vr[4][4], vi[4][4], temp_wr, temp_wi, temp_vr[4], temp_vi[4];  
  double temp1, temp2;
  double theta1, theta2, tempx, tempy, x1,y1;
  double a44[4][4], b44[4][4];
  double det, scale;
  
  //---solve eigen problem to construct A[4*4]
  for(i=0;i<4;i++)
   for(j=0;j<4;j++) m44[i][j]=linename.Cell[linename.Ncell-1]->M[i*6+j];
  EigenSolver(m44, wr, wi, vr, vi);

  //---link eigen vetcors to  H/V planes 
  temp1=abs(vr[0][0] ) + abs( vi[0][0] ) + abs( vr[0][1] ) + abs( vi[0][1] );
  temp2=abs(vr[2][0] ) + abs( vi[2][0] ) + abs( vr[2][1] ) + abs( vi[2][1] );
  if (temp2 > temp1 ) {
      temp_wr    =  wr[0] ;
      temp_wi    =  wi[0] ;
      wr[0]      =  wr[2] ;
      wi[0]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[3]      =  wr[2];
      wi[3]      = -wi[2]; 
    for(i=0;i<4;i++){
      temp_vr[i]    = vr[0][i];
      temp_vi[i]    = vi[0][i];
      vr[0][i]      = vr[2][i];
      vi[0][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
    }
  }

  //----rotating to make A12=0, A34=0
  theta1=atan2(vi[0][0], vr[0][0]);
  x1  = cos(-M_PI/2-theta1);
  y1  =sin(-M_PI/2-theta1);
  for(i=0;i<4;i++){
    tempx=vr[0][i]*x1-vi[0][i]*y1;
    tempy=vi[0][i]*x1+vr[0][i]*y1;
    vr[0][i]=tempx;
    vi[0][i]=tempy; }
  x1  = cos(M_PI/2+theta1);
  y1 =sin(M_PI/2+theta1);
  for(i=0;i<4;i++){
    tempx=vr[1][i]*x1-vi[1][i]*y1;
    tempy=vi[1][i]*x1+vr[1][i]*y1;
    vr[1][i]=tempx;
    vi[1][i]=tempy;}

  theta2=atan2(vi[2][2], vr[2][2]);
  x1  = cos(-M_PI/2-theta2);
  y1 =sin(-M_PI/2-theta2);
  for(i=0;i<4;i++){
    tempx=vr[2][i]*x1-vi[2][i]*y1;
    tempy=vi[2][i]*x1+vr[2][i]*y1;
    vr[2][i]=tempx;
    vi[2][i]=tempy; }
  x1  = cos(M_PI/2+theta2);
  y1=sin(M_PI/2+theta2);
  for(i=0;i<4;i++){
    tempx=vr[3][i]*x1-vi[3][i]*y1;
    tempy=vi[3][i]*x1+vr[3][i]*y1;
    vr[3][i]=tempx;
    vi[3][i]=tempy;}

  //----produce A
  for(i=0;i<4;i++) a44[i][0]=-(vi[0][i]-vi[1][i]);
  for(i=0;i<4;i++) a44[i][1]=vr[0][i]+vr[1][i];
  for(i=0;i<4;i++) a44[i][2]=-(vi[2][i]-vi[3][i]);
  for(i=0;i<4;i++) a44[i][3]=vr[2][i]+vr[3][i];

  if(true){
    if(a44[1][1] < 0. ) {
      for(i=0;i<4;i++) a44[i][1]=- a44[i][1];
      //cout<<"Mode I : negative beta reverted."<<endl; 
    }
    if(a44[3][3] < 0. ) {
      for(i=0;i<4;i++) a44[i][3]=- a44[i][3];
      //cout<<"Mode II : negative beta reverted."<<endl; 
    }
  }
  
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) b44[i][j]=a44[i][j];
  }
  det= mat_det( &b44[0][0], 4);

  scale=sqrt(sqrt(abs(det)));
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) a44[i][j]=a44[i][j]/scale;
  }
  
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) 
      linename.Cell[linename.Ncell-1]->A[i*6+j]=a44[i][j];

 if(false) {
    cout<<"A matrix :"<<endl;
    for (i=0;i<4;i++) {
      for(j=0;j<4;j++) cout<<setw(12)<<a44[i][j]<<"  ";
      cout<<endl;
    }
    
    for(i=0;i<4;i++) {
      for(j=0;j<4;j++) b44[i][j]=a44[i][j];
    }
    det= mat_det( &b44[0][0], 4);
    cout<<" Det of A  =  "<< det <<endl;
  }

}

void Trace_A(Line & linename, double deltap)
{
  int i,j,k;
  double T[4][4], A[4][4], G[4][4];
  double dphi1, dphi2, scale1, scale2;
  double mu1, mu2;
 
  for(i=0;i<4;i++)
    for(j=0;j<4;j++) A[i][j]= linename.Cell[linename.Ncell-1]->A[i*6+j];

  mu1=0.;
  mu2=0.;
  for(k=0;k<linename.Ncell;k++){
    
    for (i=0;i<4;i++){
      for(j=0;j<4;j++) 	T[i][j]=linename.Cell[k]->T[i*6+j];
    }
    
    mat_mult( &T[0][0], &A[0][0], &G[0][0], 4,4,4 );
    
    dphi1=atan2(G[0][1], G[0][0]);
    dphi2=atan2(G[2][3], G[2][2]);
    scale1=sqrt( G[0][0]*G[0][0] + G[0][1]*G[0][1]);
    scale2=sqrt( G[2][2]*G[2][2] + G[2][3]*G[2][3]);

    A[0][0]= ( G[0][0]*G[0][0]+G[0][1]*G[0][1] ) / scale1;
    A[0][1]= (-G[0][0]*G[0][1]+G[0][1]*G[0][0] ) / scale1;
    A[1][0]= ( G[1][0]*G[0][0]+G[1][1]*G[0][1] ) / scale1;
    A[1][1]= (-G[1][0]*G[0][1]+G[1][1]*G[0][0] ) / scale1;
                                                                          
    A[2][0]= ( G[2][0]*G[0][0]+G[2][1]*G[0][1] ) / scale1;
    A[2][1]= (-G[2][0]*G[0][1]+G[2][1]*G[0][0] ) / scale1;
    A[3][0]= ( G[3][0]*G[0][0]+G[3][1]*G[0][1] ) / scale1;
    A[3][1]= (-G[3][0]*G[0][1]+G[3][1]*G[0][0] ) / scale1; 
                                                                          
    A[0][2]= ( G[0][2]*G[2][2]+G[0][3]*G[2][3] ) / scale2;
    A[0][3]= (-G[0][2]*G[2][3]+G[0][3]*G[2][2] ) / scale2;
    A[1][2]= ( G[1][2]*G[2][2]+G[1][3]*G[2][3] ) / scale2;
    A[1][3]= (-G[1][2]*G[2][3]+G[1][3]*G[2][2] ) / scale2;
                                                                          
    A[2][2]= ( G[2][2]*G[2][2]+G[2][3]*G[2][3] ) / scale2;
    A[2][3]= (-G[2][2]*G[2][3]+G[2][3]*G[2][2] ) / scale2;
    A[3][2]= ( G[3][2]*G[2][2]+G[3][3]*G[2][3] ) / scale2;
    A[3][3]= (-G[3][2]*G[2][3]+G[3][3]*G[2][2] ) / scale2; 

    for(i=0;i<4;i++){
      for(j=0;j<4;j++) linename.Cell[k]->A[i*6+j]=A[i][j];
    }
    
    mu1=mu1+dphi1/ 2./M_PI;
    mu2=mu2+dphi2/ 2./M_PI;
    linename.Cell[k]->Mu1= mu1;
    linename.Cell[k]->Mu2= mu2;
  }

  linename.Tune1=  mu1;
  linename.Tune2=  mu2;
}

void Cal_Twiss(Line & linename, double deltap)
{
  int i,j,k;
  double temp;
  double A12[4], A22[4], C[4];

  Cal_Orbit_Num(linename, deltap);
  Cal_OneTurnMap(linename, deltap);
  Cal_A(linename, deltap);
  Cal_ElementMap(linename, deltap);
  Trace_A(linename, deltap);
  
  for(k=0; k<linename.Ncell;k++){
    linename.Cell[k]->Beta1= linename.Cell[k]->A[0*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Alfa1=-linename.Cell[k]->A[1*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Beta2= linename.Cell[k]->A[2*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Alfa2=-linename.Cell[k]->A[3*6+2]/linename.Cell[k]->A[3*6+3];
    
    temp=linename.Cell[k]->A[0*6+0]*linename.Cell[k]->A[1*6+1];
    temp=temp*linename.Cell[k]->A[2*6+2]*linename.Cell[k]->A[3*6+3];
    temp=sqrt(sqrt(temp));
    linename.Cell[k]->r=temp;
    
    for(i=0;i<2;i++){
      for(j=0;j<2;j++) {
	A12[i*2+j]= linename.Cell[k]->A[6*i+(j+2)];
	A22[i*2+j]= linename.Cell[k]->A[6*(i+2)+(j+2)];
      }    
    }
    
    mat_inv(A22,2);
    mat_mult(A12, A22,C,2,2,2);
    linename.Cell[k]->c11=temp*C[0];
    linename.Cell[k]->c12=temp*C[1];
    linename.Cell[k]->c21=temp*C[2];
    linename.Cell[k]->c22=temp*C[3];
  }
}

//------trace the Twiss parameters through the ring to the end
void Trace_Twiss(Line & linename, double deltap, double x[4], double Beta1, double Beta2, double Alfa1, double Alfa2, double c11, double c12, double c21, double c22)
{
  double r;
  double x0[6];
  int i,j,k;
  double temp;
  double A12[4], A22[4], C[4];

  for(i=0;i<4;i++)  x0[i]=x[i];
  x0[4]=0.000;
  x0[5]=deltap;
  for(j=0;j<linename.Ncell; j++) {
    linename.Cell[j]->Pass(x0);
    for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
  } 
  
  r= sqrt(1 - ( c11*c22-c12*c21));
  linename.Cell[linename.Ncell-1]->A[0*6+0]= r * sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[0*6+1]= 0; 
  linename.Cell[linename.Ncell-1]->A[0*6+2]= c11*sqrt(Beta2)- c12*Alfa2/sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[0*6+3]= c12 / sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[1*6+0]= -Alfa1*r / sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[1*6+1]= r/sqrt(Beta1); 
  linename.Cell[linename.Ncell-1]->A[1*6+2]= c21*sqrt(Beta2)-c22*Alfa2/sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[1*6+3]= c22 / sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[2*6+0]=  -c12*Alfa1 / sqrt(Beta1) - c22*sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[2*6+1]=  c12/sqrt(Beta1); 
  linename.Cell[linename.Ncell-1]->A[2*6+2]= r*sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[2*6+3]= 0.;
  linename.Cell[linename.Ncell-1]->A[3*6+0]= c11*Alfa1 / sqrt(Beta1) + c21*sqrt(Beta1);
  linename.Cell[linename.Ncell-1]->A[3*6+1]= -c11/sqrt(Beta1); 
  linename.Cell[linename.Ncell-1]->A[3*6+2]= -Alfa2*r/sqrt(Beta2);
  linename.Cell[linename.Ncell-1]->A[3*6+3]= r/sqrt(Beta2);
  
  Cal_ElementMap(linename, deltap);
  Trace_A(linename, deltap);

  for(k=0; k<linename.Ncell;k++){
    linename.Cell[k]->Beta1= linename.Cell[k]->A[0*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Alfa1=-linename.Cell[k]->A[1*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Beta2= linename.Cell[k]->A[2*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Alfa2=-linename.Cell[k]->A[3*6+2]/linename.Cell[k]->A[3*6+3];
    
    temp=linename.Cell[k]->A[0*6+0]*linename.Cell[k]->A[1*6+1];
    temp=temp*linename.Cell[k]->A[2*6+2]*linename.Cell[k]->A[3*6+3];
    temp=sqrt(sqrt(temp));
    linename.Cell[k]->r=temp;
    
    for(i=0;i<2;i++){
      for(j=0;j<2;j++) {
	A12[i*2+j]= linename.Cell[k]->A[6*i+(j+2)];
	A22[i*2+j]= linename.Cell[k]->A[6*(i+2)+(j+2)];
      }    
    }
    
    mat_inv(A22,2);
    mat_mult(A12, A22,C,2,2,2);
    linename.Cell[k]->c11=temp*C[0];
    linename.Cell[k]->c12=temp*C[1];
    linename.Cell[k]->c21=temp*C[2];
    linename.Cell[k]->c22=temp*C[3];
  }
}

//-----trace the orbit through the ring to the end 
void Trace_Orbit(Line & linename, double x[6], int i1)
{
  int i,j;
  for(i=i1;i<linename.Ncell;i++){
    linename.Cell[i]->Pass(x);
    for(j=0;j<6;j++) linename.Cell[i]->X[j]=x[j];
  }
}

//---6D orbit and Twiss
void Cal_Orbit_Num_6D(Line & linename)
{
  int i,j,k,iter=0;
  double x0[6],x01[6], x1[6], dx[6];
  double d=1.0e-09;
  double mat[6][6];
  double chi;
  double codeps=1e-12;
  int Max_iter=20;
  
  for (i=0;i<6;i++) x0[i]=0.;

  do{
    iter++;
    
    for(j=0;j<6;j++) x01[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x01);
    
    for(k=0;k<6;k++){
      for(j=0;j<6;j++) x1[j]=x0[j];
      x1[k]=x1[k]+d;
      for(j=0;j<linename.Ncell;j++) linename.Cell[j]->Pass(x1);
      for(j=0;j<6;j++) mat[j][k]=(x1[j]-x01[j])/d;
    }
    
    for(j=0;j<6;j++) x1[j]=x0[j];
    for(j=0;j<linename.Ncell;j++)  linename.Cell[j]->Pass(x1);
    for(i=0;i<6;i++) {
      dx[i]=x0[i]-x1[i];
    }
    chi=0;
    for(i=0;i<6;i++) chi+=dx[i]*dx[i];
    chi=sqrt(chi/6.);
    
    if( chi > codeps ){
      for(i=0;i<6;i++) mat[i][i]=mat[i][i]-1.0000001;
      mat_inv(&mat[0][0],6); 
      for(i=0;i<6;i++){
	for(j=0;j<6;j++) x0[i]+=mat[i][j]*dx[j];  
      } 
    } 
    
  }while (chi>codeps && iter <Max_iter );
  
  if(iter == Max_iter-1 ) 
    {
      cout<<"Failed to find COD."<<endl;
      exit(1);
    }
  else
    {
      for(j=0;j<linename.Ncell; j++) {
	linename.Cell[j]->Pass(x0);
	for(i=0;i<6;i++) linename.Cell[j]->X[i]=x0[i];
      } 
    }
}

void Cal_OneTurnMap_6D(Line & linename)
// there is not such clsoed orbits for different constant deltaps like in 5-d simulation
{
  int i,j;
  double x[6];
  tps tps1[6];
  linmap m1;
  double m66[36], b66[36];
  int flag;
  double u[6], v[6];

  for(i=0;i<6;i++) x[i]=linename.Cell[linename.Ncell-1]->X[i];
  m1.identity();
  m1=m1+x;
  
  for(i=0;i<6;i++) tps1[i]=m1[i];
  for(j=0;j<linename.Ncell; j++) linename.Cell[j]->DAPass(tps1);
  for(i=0;i<6;i++) m1[i]=tps1[i];
  Getmat(m1, m66); 

  for(i=0;i<36;i++) linename.Cell[linename.Ncell-1]->M[i]=m66[i];
  
  if( false ) {
    cout<<"One turn map M:"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<m66[i*6+j]<<"  ";
      cout<<endl;
    }
    for(i=0;i<36;i++) b66[i]=m66[i];
    cout<<"Det of M = "<<mat_det(b66,6)<<endl;
  }

  mat_change_hessenberg(m66, 6);
  flag=mat_root_hessenberg(m66,6,u,v,1.0e-10,60);
  if(flag < 0){
    cout<<"Unstable motion. "<<endl;
    exit(1);
  }
}

void Cal_ElementMap_6D(Line & linename)
{
  int i,j,k, istart;
  double x[6];
  tps tps1[6];
  linmap m1;
  double t66[36];

  for(k=0; k<linename.Ncell; k++) 
    {
      istart=k-1;
      if (istart <0) istart=linename.Ncell-1;

      for(i=0;i<6;i++) x[i]=linename.Cell[istart]->X[i];
      m1.identity();
      m1=m1+x;
      
      for(i=0;i<6;i++) tps1[i]=m1[i];
      linename.Cell[k]->DAPass(tps1);
      for(i=0;i<6;i++) m1[i]=tps1[i];
      Getmat(m1, t66);
      
      for(j=0;j<36;j++) linename.Cell[k]->T[j]=t66[j]; 
    }
}

void Cal_A_6D(Line & linename)
{
  
  int i,j;
  double m66[6][6];
  double wr[6], wi[6], vr[6][6], vi[6][6], temp_wr, temp_wi, temp_vr[6], temp_vi[6];  
  double temp1, temp2, temp3;
  double theta1, theta2, theta3, tempx, tempy, x1,y1;
  double a66[6][6], b66[6][6], c44[4][4];
  double det, scale;
  
  //---solve eigen problem to construct A[6*6]
  for(i=0;i<6;i++)
   for(j=0;j<6;j++) m66[i][j]=linename.Cell[linename.Ncell-1]->M[i*6+j];

  EigenSolver_6D(m66, wr, wi, vr, vi);
  if(false){
    cout<<"Eigen values and eigen vectors: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
        cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl; 
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

 //---link eigen vetcors S planes 
  temp1=abs(atan2(wi[0], wr[0]))/2/M_PI;
  temp2=abs(atan2(wi[2], wr[2]))/2/M_PI; 
  temp3=abs(atan2(wi[4], wr[4]))/2/M_PI;
  cout<<"Tunes: "<<temp1<<"  "<<temp2<<"   "<<temp3<<endl;
  if( temp1 < temp2 and temp1 < temp3) {
      temp_wr    =  wr[4] ;
      temp_wi    =  wi[4] ;
      wr[4]      =  wr[0] ;
      wi[4]      =  wi[0] ;
      wr[0]      =  temp_wr;
      wi[0]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[5]      =  wr[4];
      wi[5]      = -wi[4]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[4][i];
      temp_vi[i]    = vi[4][i];
      vr[4][i]      = vr[0][i];
      vi[4][i]      = vi[0][i];
      vr[0][i]      = temp_vr[i];
      vi[0][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[5][i]      = vr[4][i];
      vi[5][i]      =-vi[4][i]; 
    }
  }
  
  if( temp2 < temp1 and temp2 < temp3) {
      temp_wr    =  wr[4] ;
      temp_wi    =  wi[4] ;
      wr[4]      =  wr[2] ;
      wi[4]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[3]      = wr[2];
      wi[3]      =-wi[2]; 
      wr[5]      = wr[4];
      wi[5]      =-wi[4]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[4][i];
      temp_vi[i]    = vi[4][i];
      vr[4][i]      = vr[2][i];
      vi[4][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
      vr[5][i]      = vr[4][i];
      vi[5][i]      =-vi[4][i]; 
    }
  }

  if(false){
    cout<<"Eigen values and eigen vectors after s mode fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++) cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
    }
  }

 //---link eigen vetcors H/V planes 
  temp1=abs(vr[0][0] ) + abs( vi[0][0] ) + abs( vr[0][1] ) + abs( vi[0][1] );
  temp2=abs(vr[2][0] ) + abs( vi[2][0] ) + abs( vr[2][1] ) + abs( vi[2][1] );
  if (temp2 > temp1 ) {
      temp_wr    =  wr[0] ;
      temp_wi    =  wi[0] ;
      wr[0]      =  wr[2] ;
      wi[0]      =  wi[2] ;
      wr[2]      =  temp_wr;
      wi[2]      =  temp_wi;
      wr[1]      =  wr[0];
      wi[1]      = -wi[0]; 
      wr[3]      =  wr[2];
      wi[3]      = -wi[2]; 
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[0][i];
      temp_vi[i]    = vi[0][i];
      vr[0][i]      = vr[2][i];
      vi[0][i]      = vi[2][i];
      vr[2][i]      = temp_vr[i];
      vi[2][i]      = temp_vi[i];
      vr[1][i]      = vr[0][i];
      vi[1][i]      =-vi[0][i]; 
      vr[3][i]      = vr[2][i];
      vi[3][i]      =-vi[2][i]; 
    }
  }

  if(false){
    cout<<"Eigen values and eigen vectors after H/V/S modes fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
	cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl;
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

  //----if here we swap eigen vetcor 4 and 5 
  if(false){
    for(i=0;i<6;i++){
      temp_vr[i]    = vr[4][i];
      temp_vi[i]    = vi[4][i];
      vr[4][i]      = vr[5][i];
      vi[4][i]      = vi[5][i];
      vr[5][i]      = temp_vr[i];
      vi[5][i]      = temp_vi[i];
    }
  }

  //----rotating to make A12=0, A34=0, A56=0
  theta1=atan2(vi[0][0], vr[0][0]);
  x1  =cos(-M_PI/2-theta1);
  y1  =sin(-M_PI/2-theta1);
  for(i=0;i<6;i++){
    tempx=vr[0][i]*x1-vi[0][i]*y1;
    tempy=vi[0][i]*x1+vr[0][i]*y1;
    vr[0][i]=tempx;
    vi[0][i]=tempy; }
  x1 =cos(M_PI/2+theta1);
  y1 =sin(M_PI/2+theta1);
  for(i=0;i<6;i++){
    tempx=vr[1][i]*x1-vi[1][i]*y1;
    tempy=vi[1][i]*x1+vr[1][i]*y1;
    vr[1][i]=tempx;
    vi[1][i]=tempy;}

  theta2=atan2(vi[2][2], vr[2][2]);
  x1 =cos(-M_PI/2-theta2);
  y1 =sin(-M_PI/2-theta2);
  for(i=0;i<6;i++){
    tempx=vr[2][i]*x1-vi[2][i]*y1;
    tempy=vi[2][i]*x1+vr[2][i]*y1;
    vr[2][i]=tempx;
    vi[2][i]=tempy; }
  x1=cos(M_PI/2+theta2);
  y1=sin(M_PI/2+theta2);
  for(i=0;i<6;i++){
    tempx=vr[3][i]*x1-vi[3][i]*y1;
    tempy=vi[3][i]*x1+vr[3][i]*y1;
    vr[3][i]=tempx;
    vi[3][i]=tempy;}

  theta3=atan2(vi[4][4], vr[4][4]);
  x1 =cos(-M_PI/2-theta3);
  y1 =sin(-M_PI/2-theta3);
  for(i=0;i<6;i++){
    tempx=vr[4][i]*x1-vi[4][i]*y1;
    tempy=vi[4][i]*x1+vr[4][i]*y1;
    vr[4][i]=tempx;
    vi[4][i]=tempy; }
  x1=cos(M_PI/2+theta3);
  y1=sin(M_PI/2+theta3);
  for(i=0;i<6;i++){
    tempx=vr[5][i]*x1-vi[5][i]*y1;
    tempy=vi[5][i]*x1+vr[5][i]*y1;
    vr[5][i]=tempx;
    vi[5][i]=tempy;}

  if(false){
    cout<<"Eigen values and eigen vectors after phase fixed: "<<endl;
    for(i=0;i<6;i++){
      cout<<"    "<<i<<"  :  "<<wr[i]<< "  +j  "<<wi[i]<<endl;
      for(j=0;j<6;j++){
	cout<<sqrt(vr[i][j]*vr[i][j] +vi[i][j]*vi[i][j]  ) <<"  "<<atan2(vi[i][j],vr[i][j])<<endl;
	//cout<<vr[i][j]<<" +j "<<vi[i][j]<<endl;
      }
    }
  }

  //----produce A
  for(i=0;i<6;i++) a66[i][0]=-(vi[0][i]-vi[1][i]);
  for(i=0;i<6;i++) a66[i][1]=  vr[0][i]+vr[1][i];
  for(i=0;i<6;i++) a66[i][2]=-(vi[2][i]-vi[3][i]);
  for(i=0;i<6;i++) a66[i][3]=  vr[2][i]+vr[3][i];
  for(i=0;i<6;i++) a66[i][4]=-(vi[4][i]-vi[5][i]);
  for(i=0;i<6;i++) a66[i][5]=  vr[4][i]+vr[5][i];

  if(true){
    if(a66[1][1] < 0. ) {
      for(i=0;i<6;i++) a66[i][1]=- a66[i][1];
      //cout<<"Mode I : negative beta reverted."<<endl; 
   }
    if(a66[3][3] < 0. ) {
      for(i=0;i<6;i++) a66[i][3]=- a66[i][3];
      //cout<<"Mode II : negative beta reverted."<<endl; 
    }
    if(a66[5][5] < 0. ) {
      for(i=0;i<6;i++) a66[i][5]=- a66[i][5];
      //cout<<"Mode III: negative beta reverted."<<endl; 
    }
  }
  
  if(false) {
    cout<<"A matrix :"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<a66[i][j]<<"  ";
      cout<<endl;
    }
    
    for(i=0;i<6;i++) {
      for(j=0;j<6;j++) b66[i][j]=a66[i][j];
    }
    det= mat_det( &b66[0][0], 6);
    cout<<" Det of A66  =  "<< det <<endl;
  }
  
  //---normalizing 
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)  c44[i][j]=a66[i][j];
  det= mat_det( &c44[0][0], 4);
  scale=sqrt(sqrt(abs(det)));
  for(i=0;i<6;i++) {
    for(j=0;j<6;j++) a66[i][j]=a66[i][j]/scale;
  }
  
  for(i=0;i<6;i++) 
    for(j=0;j<6;j++) 
      linename.Cell[linename.Ncell-1]->A[i*6+j]=a66[i][j];
  
 if(false) {
    cout<<"A matrix :"<<endl;
    for (i=0;i<6;i++) {
      for(j=0;j<6;j++) cout<<setw(12)<<a66[i][j]<<"  ";
      cout<<endl;
    }
    
    for(i=0;i<6;i++) {
      for(j=0;j<6;j++) b66[i][j]=a66[i][j];
    }
    det= mat_det( &b66[0][0], 6);
    cout<<" Det of A  =  "<< det <<endl;
  }

}

void Trace_A_6D(Line & linename)
{
  int i,j,k;
  double T[6][6], A[6][6], G[6][6], RI[6][6];
  double dphi1, dphi2, dphi3;
  double mu1, mu2, mu3;
 
  for(i=0;i<6;i++)
    for(j=0;j<6;j++) A[i][j]= linename.Cell[linename.Ncell-1]->A[i*6+j];

  mu1=0.;
  mu2=0.;
  for(k=0;k<linename.Ncell;k++){
    for (i=0;i<6;i++){
      for(j=0;j<6;j++) 	T[i][j]=linename.Cell[k]->T[i*6+j];
    }
  
    mat_mult( &T[0][0], &A[0][0], &G[0][0], 6,6,6 );
    dphi1=atan2(G[0][1], G[0][0]);
    dphi2=atan2(G[2][3], G[2][2]);
    dphi3=atan2(G[4][5], G[4][4]);

    for(i=0;i<6;i++)
      for(j=0;j<6;j++) RI[i][j]=0;

    RI[0][0]  =  cos(dphi1);   
    RI[0][1]  = -sin(dphi1);  
    RI[1][0]  =  sin(dphi1);  
    RI[1][1]  =  cos(dphi1);

    RI[2][2]  =  cos(dphi2);   
    RI[2][3]  = -sin(dphi2);  
    RI[3][2]  =  sin(dphi2);  
    RI[3][3]  =  cos(dphi2);

    RI[4][4]  =  cos(dphi3);   
    RI[4][5]  = -sin(dphi3);  
    RI[5][4]  =  sin(dphi3);  
    RI[5][5]  =  cos(dphi3);

   mat_mult( &G[0][0], &RI[0][0], &A[0][0], 6,6,6 );
   for(i=0;i<6;i++){
     for(j=0;j<6;j++) linename.Cell[k]->A[i*6+j]=A[i][j];
   }
    
    mu1=mu1+dphi1/ 2./M_PI;
    mu2=mu2+dphi2/ 2./M_PI;
    mu3=mu3+dphi3/ 2./M_PI;
    linename.Cell[k]->Mu1= mu1;
    linename.Cell[k]->Mu2= mu2;
    linename.Cell[k]->Mu3= mu3;
  }

  linename.Tune1=  mu1;
  linename.Tune2=  mu2;
  linename.Tune3=  mu3;
}

void Cal_Twiss_6D(Line & linename)
{
  int i,j,k;
  double temp;
  double A12[4], A22[4], C[4];

  Cal_Orbit_Num_6D(linename);
  Cal_OneTurnMap_6D(linename);
  Cal_A_6D(linename);
  Cal_ElementMap_6D(linename);
  Trace_A_6D(linename);
  
  for(k=0; k<linename.Ncell;k++){
    linename.Cell[k]->Beta1= linename.Cell[k]->A[0*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Alfa1=-linename.Cell[k]->A[1*6+0]/linename.Cell[k]->A[1*6+1];
    linename.Cell[k]->Beta2= linename.Cell[k]->A[2*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Alfa2=-linename.Cell[k]->A[3*6+2]/linename.Cell[k]->A[3*6+3];
    linename.Cell[k]->Beta3= linename.Cell[k]->A[4*6+4]/linename.Cell[k]->A[5*6+5];
    linename.Cell[k]->Alfa3=-linename.Cell[k]->A[5*6+4]/linename.Cell[k]->A[5*6+5];
    
    temp=linename.Cell[k]->A[0*6+0]*linename.Cell[k]->A[1*6+1];
    temp=temp*linename.Cell[k]->A[2*6+2]*linename.Cell[k]->A[3*6+3];
    temp=sqrt(sqrt(temp));
    linename.Cell[k]->r=temp;
    
    for(i=0;i<2;i++){
      for(j=0;j<2;j++) {
	A12[i*2+j]= linename.Cell[k]->A[6*i+(j+2)];
	A22[i*2+j]= linename.Cell[k]->A[6*(i+2)+(j+2)];
      }    
    }
    
    mat_inv(A22,2);
    //C[1]=-A22[1];
    //C[2]=-A22[2];
    //C[3]= A22[0];
    //det = (A22[0]*A22[3]-A22[1]*A22[2]);  
    //for(i=0;i<3;i++) A22[i]= C[0] / det ;
    //if(det==0){
    //  cout<<" Can't invert A22 in Cal_Twiss(). "<<endl;
    //  exit(1);
    //}

    mat_mult(A12, A22,C,2,2,2);
    linename.Cell[k]->c11=temp*C[0];
    linename.Cell[k]->c12=temp*C[1];
    linename.Cell[k]->c21=temp*C[2];
    linename.Cell[k]->c22=temp*C[3];
  }
}

void Cal_Chrom( Line & linename)
{
  double deltap;
  double qx0, qy0, qxp, qyp, qxm, qym;

  deltap=0.0003;
  Cal_Twiss(linename, deltap);
  qxp=linename.Tune1;
  qyp=linename.Tune2;
  
  deltap=-0.0003;
  Cal_Twiss(linename, deltap);
  qxm=linename.Tune1;
  qym=linename.Tune2;

  deltap=0.00;
  Cal_Twiss(linename,deltap);
  qx0=linename.Tune1;
  qy0=linename.Tune2;
  
  linename.Chromx1=(qxp-qxm)/2./0.0003;
  linename.Chromy1=(qyp-qym)/2./0.0003;
  linename.Chromx2=(qxp+qxm-2*qx0)/2./0.0003/0.0003;
  linename.Chromy2=(qyp+qym-2*qy0)/2./0.0003/0.0003;
}

void Cal_Dispersion(Line & linename)
{
  int i,j, k;
  double m66[6][6];
  double r44[4][4];
  double x[6], xtemp[6];  
  double T66[6][6];
   
  for(i=0;i<6;i++) 
    for(j=0;j<6;j++) m66[i][j]= linename.Cell[linename.Ncell-1]->M[i*6+j];
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) r44[i][j]= linename.Cell[linename.Ncell-1]->M[i*6+j];
  for(i=0;i<4;i++) r44[i][i]  = r44[i][i]-1.0;
  mat_inv(&r44[0][0],4); 
  for(i=0;i<4;i++) 
    for(j=0;j<4;j++) r44[i][j]= -1.0 * r44[i][j];
  linename.Cell[linename.Ncell-1]->Etax  = r44[0][0] * m66[0][5] +  r44[0][1] * m66[1][5] +   r44[0][2] * m66[2][5] +  r44[0][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etaxp = r44[1][0] * m66[0][5] +  r44[1][1] * m66[1][5] +   r44[1][2] * m66[2][5] +  r44[1][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etay  = r44[2][0] * m66[0][5] +  r44[2][1] * m66[1][5] +   r44[2][2] * m66[2][5] +  r44[2][3] * m66[3][5];
  linename.Cell[linename.Ncell-1]->Etayp = r44[3][0] * m66[0][5] +  r44[3][1] * m66[1][5] +   r44[3][2] * m66[2][5] +  r44[3][3] * m66[3][5];

  x[0]= linename.Cell[linename.Ncell-1]->Etax  ;
  x[1]= linename.Cell[linename.Ncell-1]->Etaxp ;
  x[2]= linename.Cell[linename.Ncell-1]->Etay  ;
  x[3]= linename.Cell[linename.Ncell-1]->Etayp ;
  for(k=0;k<linename.Ncell; k++){
    x[4]= 0.;  x[5]= 1;
    for(i=0;i<6;i++) 
      for(j=0;j<6;j++) T66[i][j]= linename.Cell[k]->T[i*6+j];
    xtemp[0]=  x[0] * T66[0][0] +  x[1] * T66[0][1] +   x[2] * T66[0][2] +  x[3] * T66[0][3]  +   x[5] * T66[0][5]   ;
    xtemp[1]=  x[0] * T66[1][0] +  x[1] * T66[1][1] +   x[2] * T66[1][2] +  x[3] * T66[1][3]  +   x[5] * T66[1][5]   ;
    xtemp[2]=  x[0] * T66[2][0] +  x[1] * T66[2][1] +   x[2] * T66[2][2] +  x[3] * T66[2][3]  +   x[5] * T66[2][5]   ;
    xtemp[3]=  x[0] * T66[3][0] +  x[1] * T66[3][1] +   x[2] * T66[3][2] +  x[3] * T66[3][3]  +   x[5] * T66[3][5]   ;
    linename.Cell[k]->Etax  =  xtemp[0];
    linename.Cell[k]->Etaxp =  xtemp[1];
    linename.Cell[k]->Etay  =  xtemp[2];
    linename.Cell[k]->Etayp =  xtemp[3];
    for(i=0;i<4;i++) x[i]= xtemp[i];
  }
}

void Cal_Tune_Num(Line & linename, double deltap0)
{
  int k;
  double x[6];

  int Nturn=1024;
  int stable = 1, lost_turn=0, lost_post=0;
  double  x_tbt[Nturn*6], xtemp[1024];
  double tunex1, tuney1;
  
  Cal_Orbit_Num(linename, deltap0);

  x[0]=linename.Cell[linename.Ncell-1]->X[0] ;
  x[1]=linename.Cell[linename.Ncell-1]->X[1] ;
  x[2]=linename.Cell[linename.Ncell-1]->X[2] ;
  x[3]=linename.Cell[linename.Ncell-1]->X[3] ;
  x[4]=0. ;
  x[5]= deltap0;
  
  if(x[0] < 1.0e-9 ) x[0]=x[0]+1.0e-9;
  if(x[2] < 1.0e-9 ) x[2]=x[2]+1.0e-9; 

  stable = 1; lost_turn=0; lost_post=0;
  Track_tbt(linename, x, Nturn, x_tbt, stable, lost_turn, lost_post);

  if(stable==1) {
    for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+0];
    FineTuneFinder(1024, xtemp, tunex1);
    for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+2];
    FineTuneFinder(1024, xtemp, tuney1);
  }
  linename.Tune1= tunex1;
  linename.Tune2= tuney1;
}

//-----------------------------------------
//   fitting tunes and linear chroms
//-----------------------------------------
void Fit_Tune(Line & linename, double q1, double q2, const char * qf_name, const char * qd_name)
{
  double tunex0, tuney0, tunex1=q1, tuney1=q2, dtunex, dtuney;
  double qf_k1l_0, qd_k1l_0;
  double dk1l_qf, dtunex_qf,  dtuney_qf, dk1l_qd, dtunex_qd,  dtuney_qd;
  double scale_qf, scale_qd;

  Cal_Twiss(linename,0.0);
  tunex0=linename.Tune1;
  tuney0=linename.Tune2;

  while( (tunex1-tunex0)*(tunex1-tunex0) + (tuney1-tuney0)*(tuney1-tuney0)  > 1.0e-10 ) {
    qf_k1l_0= Get_KL(linename,qf_name, "K1L");
    qd_k1l_0= Get_KL(linename, qd_name,"K1L"); 
    
    dk1l_qf=  qf_k1l_0 * 0.001;
    Set_dKL(linename,qf_name, "K1L", dk1l_qf);
    Cal_Twiss(linename,0.0);
    dtunex_qf=linename.Tune1 - tunex0;
    dtuney_qf=linename.Tune2 - tuney0;
    Set_dKL(linename,qf_name, "K1L",-dk1l_qf);
    
    dk1l_qd=  qd_k1l_0 * 0.001;
    Set_dKL(linename,qd_name, "K1L", dk1l_qd);
    Cal_Twiss(linename,0.0);
    dtunex_qd=linename.Tune1 - tunex0;
    dtuney_qd=linename.Tune2 - tuney0;
    Set_dKL(linename,qd_name, "K1L",-dk1l_qd);
    
    dtunex=tunex1- tunex0;
    dtuney=tuney1- tuney0;
    
    LinearEquations(dtunex_qf, dtunex_qd,dtunex, dtuney_qf, dtuney_qd, dtuney, scale_qf, scale_qd);
    Set_dKL(linename,qf_name, "K1L", dk1l_qf * scale_qf);
    Set_dKL(linename,qd_name, "K1L", dk1l_qd * scale_qd);
    
    Cal_Twiss(linename,0.0);
    tunex0=linename.Tune1;
    tuney0=linename.Tune2;
  }
}

void Fit_Tune_RHICelens(Line & linename, double q1, double q2)
{
  double tunex0, tuney0, tunex1=q1, tuney1=q2, dtunex, dtuney;
  double qf_k1l_0, qd_k1l_0;
  double dk1l_qf, dtunex_qf,  dtuney_qf, dk1l_qd, dtunex_qd,  dtuney_qd;
  double scale_qf, scale_qd;

  Cal_Twiss(linename,0.0);
  tunex0=linename.Tune1;
  tuney0=linename.Tune2;

  while( (tunex1-tunex0)*(tunex1-tunex0) + (tuney1-tuney0)*(tuney1-tuney0)  > 1.0e-10 ) {
    qf_k1l_0= Get_KL(linename,"QF", "K1L");
    qd_k1l_0= Get_KL(linename, "QD","K1L"); 
    
    dk1l_qf=  qf_k1l_0 * 0.001;
    Set_dKL(linename,"QF", "K1L", dk1l_qf);
    Set_dKL(linename,"QFSHFT", "K1L", dk1l_qf);
    Set_dKL(linename,"QFSHFT2", "K1L", dk1l_qf/2);

    Cal_Twiss(linename,0.0);
    dtunex_qf=linename.Tune1 - tunex0;
    dtuney_qf=linename.Tune2 - tuney0;
    Set_dKL(linename,"QF", "K1L",-dk1l_qf);
    Set_dKL(linename,"QFSHFT", "K1L", -dk1l_qf);
    Set_dKL(linename,"QFSHFT2", "K1L", -dk1l_qf/2);
    
    dk1l_qd=  qd_k1l_0 * 0.001;
    Set_dKL(linename,"QD", "K1L", dk1l_qd);
    Set_dKL(linename,"QDSHFT", "K1L", dk1l_qd);
    Set_dKL(linename,"QDSHFT2", "K1L", dk1l_qd/2);
    Cal_Twiss(linename,0.0);
    dtunex_qd=linename.Tune1 - tunex0;
    dtuney_qd=linename.Tune2 - tuney0;
    Set_dKL(linename,"QD", "K1L",-dk1l_qd);
    Set_dKL(linename,"QDSHFT", "K1L", -dk1l_qd);
    Set_dKL(linename,"QDSHFT2", "K1L", -dk1l_qd/2);
   
    dtunex=tunex1- tunex0;
    dtuney=tuney1- tuney0;
    
    LinearEquations(dtunex_qf, dtunex_qd,dtunex, dtuney_qf, dtuney_qd, dtuney, scale_qf, scale_qd);
    Set_dKL(linename,"QF", "K1L", dk1l_qf * scale_qf);
    Set_dKL(linename,"QFSHFT", "K1L",  dk1l_qf * scale_qf);
    Set_dKL(linename,"QFSHFT2", "K1L",  dk1l_qf * scale_qf /2 );

    Set_dKL(linename,"QD", "K1L", dk1l_qd * scale_qd);
    Set_dKL(linename,"QDSHFT", "K1L", dk1l_qd * scale_qd);
    Set_dKL(linename,"QDSHFT2", "K1L", dk1l_qd * scale_qd/2);
  
    Cal_Twiss(linename,0.0);
    tunex0=linename.Tune1;
    tuney0=linename.Tune2;
  }
}

void Fit_Chrom(Line & linename, double chrom1x_want, double chrom1y_want, const char * sf_name, const char * sd_name )
{
  double chrom1x0, chrom1y0,  dchrom1x, dchrom1y;
  double dk2l_sf, dk2l_sd, dchrom1x_sf,  dchrom1y_sf, dchrom1x_sd, dchrom1y_sd;
  double scale_sf, scale_sd;

  Cal_Chrom(linename) ; 
  chrom1x0= linename.Chromx1;
  chrom1y0= linename.Chromy1;

  while( (chrom1x_want-chrom1x0)*(chrom1x_want-chrom1x0) + (chrom1y_want-chrom1y0)*(chrom1y_want-chrom1y0)  > 0.0001 ) {
    //sf_k2l_0= Get_KL(linename,sf_name,"K2L");
    //sd_k2l_0= Get_KL(linename,sd_name,"K2L"); 
    
    dk2l_sf= 0.3 * 0.005;
    Set_dKL(linename,sf_name, "K2L", dk2l_sf);
    Cal_Chrom(linename);
    dchrom1x_sf=linename.Chromx1 - chrom1x0;
    dchrom1y_sf=linename.Chromy1 - chrom1y0;
    Set_dKL(linename,sf_name, "K2L",-dk2l_sf);
    
    dk2l_sd= -0.5 * 0.005;
    Set_dKL(linename,sd_name, "K2L", dk2l_sd);
    Cal_Chrom(linename);
    dchrom1x_sd=linename.Chromx1 - chrom1x0;
    dchrom1y_sd=linename.Chromy1 - chrom1y0;
    Set_dKL(linename,sd_name, "K2L",-dk2l_sd);

    dchrom1x=chrom1x_want- chrom1x0;
    dchrom1y=chrom1y_want- chrom1y0;
    
    LinearEquations(dchrom1x_sf, dchrom1x_sd, dchrom1x, dchrom1y_sf, dchrom1y_sd, dchrom1y, scale_sf, scale_sd);
    Set_dKL(linename,sf_name, "K2L", dk2l_sf * scale_sf);
    Set_dKL(linename,sd_name, "K2L", dk2l_sd * scale_sd);
    
    Cal_Chrom(linename);
    chrom1x0=linename.Chromx1;
    chrom1y0=linename.Chromy1;
  }
}

void Fit_Chrom_RHIC8fam(Line & linename, double chrom1x_want, double chrom1y_want )
{
  int i;
  double chrom1x0, chrom1y0,  dchrom1x, dchrom1y;
  double dk2l_sf, dk2l_sd, dchrom1x_sf,  dchrom1y_sf, dchrom1x_sd, dchrom1y_sd;
  double scale_sf, scale_sd;

  Cal_Chrom(linename);
  chrom1x0=linename.Chromx1;
  chrom1y0=linename.Chromy1;

  while( (chrom1x_want-chrom1x0)*(chrom1x_want-chrom1x0) + (chrom1y_want-chrom1y0)*(chrom1y_want-chrom1y0)  > 0.0001 ) {
    
    dk2l_sf =  0.3 * 0.005;
    for(i=0;i<linename.Ncell;i++) 
      if( linename.Cell[i]->NAME.find("SF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT") )
	linename.Cell[i]->SetP("K2L",  linename.Cell[i]->GetP("K2L") +  dk2l_sf );
    Cal_Chrom(linename);
    dchrom1x_sf=linename.Chromx1 - chrom1x0;
    dchrom1y_sf=linename.Chromy1 - chrom1y0;  
    for(i=0;i<linename.Ncell;i++) 
      if(  linename.Cell[i]->NAME.find("SF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT") ) 
	linename.Cell[i]->SetP( "K2L",linename.Cell[i]->GetP("K2L") -  dk2l_sf );

    dk2l_sd =  -0.5 * 0.005;
    for(i=0;i<linename.Ncell;i++) 
      if(  linename.Cell[i]->NAME.find("SD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT")  ) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") +  dk2l_sd );
    Cal_Chrom(linename);
    dchrom1x_sd=linename.Chromx1 - chrom1x0;
    dchrom1y_sd=linename.Chromy1 - chrom1y0;
    for(i=0;i<linename.Ncell;i++) 
       if(  linename.Cell[i]->NAME.find("SD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT") ) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") -  dk2l_sd );

    dchrom1x=chrom1x_want- chrom1x0;
    dchrom1y=chrom1y_want- chrom1y0;
    LinearEquations(dchrom1x_sf, dchrom1x_sd, dchrom1x, dchrom1y_sf, dchrom1y_sd, dchrom1y, scale_sf, scale_sd);
    for(i=0;i<linename.Ncell;i++) 
      if( linename.Cell[i]->NAME.find("SF") != string::npos and linename.Cell[i]->TYPE ==string("SEXT") ) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") + dk2l_sf * scale_sf );
    for(i=0;i<linename.Ncell;i++) 
      if( linename.Cell[i]->NAME.find("SD") != string::npos and linename.Cell[i]->TYPE ==string("SEXT") ) 
	linename.Cell[i]->SetP( "K2L", linename.Cell[i]->GetP("K2L") +dk2l_sd * scale_sd );
    
    Cal_Chrom(linename);
    chrom1x0=linename.Chromx1;
    chrom1y0=linename.Chromy1;
  }
}

//----------------------------------------
//    Chromatic Calculation 
//----------------------------------------
void  chrom_fit(double qx[],double qy[],double & chromx1,double & chromy1,double & chromx2,double & chromy2,double & chromx3,double & chromy3 )
{
  int i;
  double temp;
  double xa[21],ya[21];
  double coeff[8];

  temp=qx[10];
  for (i=0; i<21;i++){
    xa[i]=(i-10)*GP.step_deltap;
    ya[i]=qx[i]-temp;
  }
  pfit(xa, ya, 21, coeff, 7);
  chromx1=coeff[1];
  chromx2=coeff[2];
  chromx3=coeff[3];  
  
  temp=qy[10];
  for (i=0; i<21;i++){
    xa[i]=(i-10)*GP.step_deltap;
    ya[i]=qy[i]-temp;
  }
  pfit(xa, ya, 21, coeff,7);
  chromy1=coeff[1];
  chromy2=coeff[2];
  chromy3=coeff[3]; 
}

void Cal_Chrom_Num( Line & linename)
{
  int i;
  double qx[21],qy[21];
  double deltap;
  double chromx1,chromy1,chromx2,chromy2,chromx3,chromy3;  

  for(i=0;i<21;i++){
    deltap=GP.step_deltap*(i-10);
    Cal_Twiss(linename, deltap);
    qx[i]=linename.Tune1;
    qy[i]=linename.Tune2;
  }

  chrom_fit(qx,qy,chromx1,chromy1,chromx2,chromy2,chromx3,chromy3);
  
  linename.Chromx1 =  chromx1;
  linename.Chromy1 =  chromy1;
  linename.Chromx2 =  chromx2;
  linename.Chromy2 =  chromy2;
  linename.Chromx3 =  chromx3;
  linename.Chromy3 =  chromy3;

  Cal_Twiss(linename, 0.0);
}

void Correct_Chrom_Manual( Line & linename)
{
  int i;
  char fam1[16], fam2[16];
  double step;
  
  for(i=0;i<8;i++) {
    cout<<"Input two Sextupole families and change step:"<<endl;
    cin>>fam1>>fam2>>step;
    Set_dKL(linename,fam1,"K2L",step);
    Set_dKL(linename,fam2,"K2L",-step);
    Fit_Chrom_RHIC8fam(linename, 1.0, 1.0);
    Cal_Chrom(linename);  //  Cal_Chrom_Num(linename);
    cout<<linename.Chromx1<<"  "<<linename.Chromy1<<endl;
    cout<<linename.Chromx2<<"  "<<linename.Chromy2<<endl;
    //cout<<linename.Chromx3<<"  "<<linename.Chromy3<<endl;
  }
}

void Cal_Tune_vs_Deltap(Line & linename, const char *filename)
{
  int i;
  double qx[21],qy[21];
  double deltap;
  fstream fout;

  for(i=0;i<21;i++){
    deltap=GP.step_deltap*(i-10);
    Cal_Twiss(linename, deltap);
    qx[i]=linename.Tune1;
    qy[i]=linename.Tune2;
  }

  fout.open(filename, ios::out);
  for (i = 0; i < 21; i++)
    fout << setw(10) <<GP.step_deltap*(i-10)
  	 << scientific << setw(15) << qx[i]
  	 << scientific << setw(15) << qy[i]<<endl;
   fout.close();
}

void Plot_Tune_vs_Deltap(Line & linename, const char* filename)
{
  char command[256];
  fstream fout;

  fout.open("temp222.p", ios::out);
  fout<<"set term post color enhanced 20 "<<endl;
  fout<<"set output 'tune_vs_delta.ps' "<<endl;
  fout<<"set xlabel  'dp/p_0 [10^{-3}] '  " <<endl;
  fout<<"set ylabel 'Q_x' " <<endl;
  fout<<"set y2label 'Q_y' "<<endl;
  fout<<"set ytics nomirror"<<endl;
  fout<<"set y2tics"<<endl;
  
  sprintf(command, "plot '%s' u ($1*1000):2 tit 'Q_x' w l lt 1  lw 2,\\", filename );
  fout<<command<<endl;
  sprintf(command, "     '%s' u ($1*1000):3 axes x1y2 tit 'Q_y' w l lt 3  lw 2", filename );
  fout<<command<<endl;
  fout<<"exit"<<endl;
  fout.close();
  system("gnuplot temp222.p");
  system("rm temp222.p");
}

void Cal_Beta_Star_vs_Deltap(Line & linename, const char *filename)
{
  int i;
  double betx[21],bety[21];
  double deltap;
  fstream fout;

  for(i=0;i<21;i++){
    deltap=GP.step_deltap*(i-10);
    Cal_Twiss(linename, deltap);
    betx[i]=linename.Cell[linename.Ncell-1]->Beta1;
    bety[i]=linename.Cell[linename.Ncell-1]->Beta2;
  }

  fout.open(filename, ios::out);
  for (i = 0; i < 21; i++)
    fout << setw(10) <<GP.step_deltap*(i-10)
  	 << scientific << setw(15) << betx[i]
  	 << scientific << setw(15) << bety[i]<<endl;
   fout.close();
}

void Plot_Beta_Star_vs_Deltap(Line & linename, const char* filename)
{
  char command[256];
  fstream fout;

  fout.open("temp222.p", ios::out);
  fout<<"set term post color enhanced 20 "<<endl;
  fout<<"set output 'beta_star_vs_delta.ps' "<<endl;
  fout<<"set xlabel  'dp/p_0 [10^{-3}]'  " <<endl;
  fout<<"set ylabel '{/Symbol \142}*_x' " <<endl;
  fout<<"set y2label '{/Symbol \142}*_y' "<<endl;
  fout<<"set ytics nomirror"<<endl;
  fout<<"set y2tics"<<endl;
  
  sprintf(command, "plot '%s' u ($1*1000):2 tit '{/Symbol \142}*_x' w l lt 1  lw 2,\\", filename );
  fout<<command<<endl;
  sprintf(command, "     '%s' u ($1*1000):3 axes x1y2 tit '{/Symbol \142}*_y' w l lt 3  lw 2", filename );
  fout<<command<<endl;
  fout<<"exit"<<endl;
  fout.close();
  system("gnuplot temp222.p");
  system("rm temp222.p");
 }

void Cal_Beta_vs_Deltap(Line & linename, const char *filename)
{
  int i,j;
  double delta[21];
  double betax[21][linename.Ncell], betay[21][linename.Ncell];
  double DBXDD[linename.Ncell], DBYDD[linename.Ncell], DBXDD2[linename.Ncell], DBYDD2[linename.Ncell], DBXDD3[linename.Ncell], DBYDD3[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  fstream fout;
  
  for(i=0;i<21;i++) delta[i]=0.00002*(i-10);
  
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      betax[i][j]=linename.Cell[j]->Beta1;
      betay[i][j]=linename.Cell[j]->Beta2;
    }
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  betax[i][j];
      input2[i]=  betay[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DBXDD[j]=term11;   DBXDD2[j]=term21*2;  DBXDD3[j]=term31*6;
    DBYDD[j]=term12;   DBYDD2[j]=term22*2;  DBYDD3[j]=term32*6;
  }
  
  fout.open(filename, ios::out);
  for(i=0;i<linename.Ncell;i++) 
    fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->S
        <<setw(20)<<DBXDD[i]<<setw(20)<<DBXDD2[i]<<setw(20)<<DBXDD3[i]<<setw(20)<<DBYDD[i]<<setw(20)<<DBYDD2[i]<<setw(20)<<DBYDD3[i]<<endl;
  fout.close();
}

void Cal_Dispersion_vs_Deltap(Line & linename, const char *filename)
{
  int i,j;
  double delta[21];
  double xco[21][linename.Ncell], yco[21][linename.Ncell];
  double DX[linename.Ncell], DX2[linename.Ncell], DX3[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  fstream fout;
 
  for(i=0;i<21;i++) delta[i]=GP.step_deltap*(i-10);
  
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      xco[i][j]=linename.Cell[j]->X[0];
    }
  }
  
  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  xco[i][j];
      input2[i]=  yco[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DX[j]=term11;  DX2[j]=term21*2;  DX3[j]=term31*6;
  }

  fout.open(filename, ios::out);
  for(i=0;i<linename.Ncell;i++) 
    fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->S<<setw(20)<< DX[i]<<setw(20)<<DX2[i]<<setw(20)<<DX3[i]<<endl;
  fout.close();
}

void Cal_Half_Integer_RDT(Line & linename,  const char* filename)
//  horizontal: h20001, vertical: h00021
{
  int i,j;
  fstream fout;
  double h_real, h_imag, v_real, v_imag;
  double k1l, k2l, betx, bety, Dx,dphix, dphiy; 

  Cal_Twiss(linename,0.0);
  Cal_Dispersion(linename);
  fout.open(filename, ios::out);

  for(i=0;i<linename.Ncell;i++){
    h_real=0; h_imag=0;
    v_real=0; v_imag=0;
    
    for(j=0;j<linename.Ncell;j++){
      if( linename.Cell[j]->TYPE == string("QUAD") ){
	k1l=linename.Cell[j]->GetP("K1L");
	betx=linename.Cell[j]->Beta1;
	bety=linename.Cell[j]->Beta2;
	Dx=linename.Cell[j]->Etax;
	dphix= abs(linename.Cell[j]->Mu1 - linename.Cell[i]->Mu1);
	dphiy= abs(linename.Cell[j]->Mu2 - linename.Cell[i]->Mu2);        
	h_real = h_real  -  k1l*betx*cos(  2*dphix * 2 * M_PI ) ; 
	h_imag = h_imag  -  k1l*betx*sin(  2*dphix * 2 * M_PI );	
	v_real = v_real  +  k1l*bety*cos(  2*dphiy * 2 * M_PI ) ; 
	v_imag = v_imag  +  k1l*bety*sin(  2*dphiy * 2 * M_PI ) ; 
      }
      else if(linename.Cell[j]->TYPE == string("SEXT") ){
	k2l=linename.Cell[j]->GetP("K2L");
	betx=linename.Cell[j]->Beta1;
	bety=linename.Cell[j]->Beta2;
	Dx=linename.Cell[j]->Etax;
	dphix= abs(linename.Cell[j]->Mu1 - linename.Cell[i]->Mu1);
	dphiy= abs(linename.Cell[j]->Mu2 - linename.Cell[i]->Mu2);  
	h_real = h_real +  k2l*betx*Dx*cos(  2*dphix * 2 * M_PI ) ;  
	h_imag = h_imag +  k2l*betx*Dx*sin(  2*dphix * 2 * M_PI );	
	v_real = v_real -  k2l*bety*Dx*cos(  2*dphiy * 2 * M_PI ) ;  
	v_imag = v_imag -  k2l*bety*Dx*sin(  2*dphiy * 2 * M_PI ) ;  
      }
      else{}
    }
    fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->S<<setw(20)<<linename.Cell[i]->Beta1<<setw(20)<<linename.Cell[i]->Beta2<<setw(20)<<sqrt(h_real*h_real+h_imag*h_imag)<<"  "<<sqrt(v_real*v_real+v_imag*v_imag)<<endl;
  }

  fout.close();
}

void Cal_Q2_Source(Line & linename, const char* filename )
{
  int i,j;
  double delta[21];
  double xco[21][linename.Ncell], yco[21][linename.Ncell];
  double betax[21][linename.Ncell], betay[21][linename.Ncell];
  double DX[linename.Ncell], DX2[linename.Ncell], DX3[linename.Ncell];
  double DBXDD[linename.Ncell], DBYDD[linename.Ncell], DBXDD2[linename.Ncell], DBYDD2[linename.Ncell], DBXDD3[linename.Ncell], DBYDD3[linename.Ncell];
  double input1[21],input2[21];
  double term11, term12,term21, term22, term31, term32; 
  fstream fout;  

  fout.open(filename, ios::out);

  for(i=0;i<21;i++) delta[i]=0.00002*(i-10);
  for(i=0;i<21;i++){
    Cal_Twiss(linename, delta[i]);
    for(j=0;j<linename.Ncell;j++){
      xco[i][j]=linename.Cell[j]->X[0];
      betax[i][j]=linename.Cell[j]->Beta1;
      betay[i][j]=linename.Cell[j]->Beta2;
    }
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  xco[i][j];
      input2[i]=  yco[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DX[j]=term11;  DX2[j]=term21*2;  DX3[j]=term31*6;
  }

  for(j=0;j<linename.Ncell;j++) {
    for(i=0;i<21;i++){
      input1[i]=  betax[i][j];
      input2[i]=  betay[i][j];
    }
    chrom_fit(input1,input2,term11, term12,term21, term22, term31, term32);
    DBXDD[j]=term11;   DBXDD2[j]=term21*2;  DBXDD3[j]=term31*6;
    DBYDD[j]=term12;   DBYDD2[j]=term22*2;  DBYDD3[j]=term32*6;
  }

  double chromx1,   chromy1, chromx2,  chromy2,  chromx3, chromy3;
  Cal_Twiss(linename,0.);

  for(i=1;i<linename.Ncell;i++){
    
    chromx1=0;  chromy1=0 ; chromx2=0 ;  chromy2=0;  chromx3=0; chromy3=0;
    
    if(linename.Cell[i]->TYPE ==string("SEXT") ) {
      
      chromx1 +=  1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/M_PI;
      chromy1 +=- 1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/M_PI;
      
      chromx2 += -1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/M_PI/2;
      chromy2 += +1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/4/M_PI/2;
      
      chromx2 +=  linename.Cell[i]->GetP("K2L") * DX[i]  * DBXDD[i]/4/M_PI/2;
      chromy2 += -linename.Cell[i]->GetP("K2L") * DX[i]  * DBYDD[i]/4/M_PI/2;
      
      chromx2 +=  linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/4/M_PI/2;
      chromy2 += -linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/4/M_PI/2;
      
      chromx3 += +1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/8/M_PI/3;
      chromy3 += -1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5*DX[i]*linename.Cell[i]->GetP("K2L")/8/M_PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") * DX[i]  * DBXDD[i]/8/M_PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") * DX[i]  * DBYDD[i]/8/M_PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/8/M_PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") * ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/8/M_PI/3; 
      
      chromx3 +=  linename.Cell[i]->GetP("K2L") *  DX[i] * DBXDD2[i]/8/M_PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  DX[i] * DBYDD2[i]/8/M_PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") *  DX[i] * DBXDD[i]/8/M_PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") *  DX[i] * DBYDD[i]/8/M_PI/3;       
      
      
      chromx3 += +linename.Cell[i]->GetP("K2L") *  DX2[i] * DBXDD[i]/8/M_PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  DX2[i] * DBYDD[i]/8/M_PI/3; 
      
      chromx3 += +linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX3[i]/8/M_PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX3[i]/8/M_PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * DX2[i]/8/M_PI/3;
      chromy3 += +linename.Cell[i]->GetP("K2L") *  ( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * DX2[i]/8/M_PI/3;
      
      chromx3 += +linename.Cell[i]->GetP("K2L") *  DBXDD[i] * DX2[i]/8/M_PI/3;
      chromy3 += -linename.Cell[i]->GetP("K2L") *  DBYDD[i] * DX2[i]/8/M_PI/3;        
      
      fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->TYPE<<setw(20)<<linename.Cell[i]->S
          <<setw(20)<<chromx1<<setw(20)<<chromy1<<setw(20)<<chromx2<<setw(20)<<chromy2<<setw(20)<<chromx3<<setw(20)<<chromy3<<endl;
    }
    
    if(linename.Cell[i]->TYPE ==string("QUAD") ) {
      chromx1 += - 1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/4/M_PI;
      chromy1 += + 1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/4/M_PI;
      
      chromx2 += +1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/4/M_PI/2;
      chromy2 += -1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/4/M_PI/2;
      
      chromx2 += -linename.Cell[i]->GetP("K1L") * DBXDD[i]/4/M_PI/2;
      chromy2 += +linename.Cell[i]->GetP("K1L") * DBYDD[i]/4/M_PI/2;
      
      chromx3 += -1.0*( linename.Cell[i]->Beta1  + linename.Cell[i-1]->Beta1 )*0.5 * linename.Cell[i]->GetP("K1L")/8/M_PI/3;
      chromy3 += +1.0*( linename.Cell[i]->Beta2  + linename.Cell[i-1]->Beta2 )*0.5 * linename.Cell[i]->GetP("K1L")/8/M_PI/3;
      
      chromx3 += +linename.Cell[i]->GetP("K1L") * DBXDD[i]/8/M_PI/3;
      chromy3 += -linename.Cell[i]->GetP("K1L") * DBYDD[i]/8/M_PI/3;
      
      chromx3 += -linename.Cell[i]->GetP("K1L") * DBXDD2[i]/8/M_PI/3;
      chromy3 += +linename.Cell[i]->GetP("K1L") * DBYDD2[i]/8/M_PI/3;
      
      chromx3 += +linename.Cell[i]->GetP("K1L") * DBXDD[i]/8/M_PI/3;
      chromy3 += -linename.Cell[i]->GetP("K1L") * DBYDD[i]/8/M_PI/3;
      
      fout<<setw(16)<<linename.Cell[i]->NAME<<setw(20)<<linename.Cell[i]->TYPE<<setw(20)<<linename.Cell[i]->S
          <<setw(20)<<chromx1<<setw(20)<<chromy1<<setw(20)<<chromx2<<setw(20)<<chromy2<<setw(20)<<chromx3<<setw(20)<<chromy3<<endl;
    }
  }
  fout.close();
} 

//--------------------------------------------------
//  RDTs and detunings
//--------------------------------------------------
void Cal_Coupling_Coefficient( Line & linename )
// coupling coefficient from whole ring
{
  int i;
  double creal=0.0, cimag=0.0;
  double angle, sin_angle, cos_angle;
  
  for(i=0; i< linename.Ncell;i++)
    if( linename.Cell[i]->TYPE == string("SKEWQ") ){
      creal= creal + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * M_PI ) *
        linename.Cell[i]->GetP("K1SL") / 2. / M_PI  ;
      cimag= cimag + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * M_PI ) *
        linename.Cell[i]->GetP("K1SL") / 2. / M_PI  ;
    }

   for(i=0; i< linename.Ncell;i++)
    if( linename.Cell[i]->TYPE == string("QUAD")  and linename.Cell[i]->DT !=0 ){
      creal= creal + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * M_PI ) *
        linename.Cell[i]->GetP("K1L") * sin(-1.0 * linename.Cell[i]->DT*2) / 2. / M_PI  ;
      cimag= cimag + sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * M_PI ) *
        linename.Cell[i]->GetP("K1L") * sin(-1.0* linename.Cell[i]->DT*2) / 2. / M_PI  ;
    }
 
  angle=  atan( abs(cimag/creal) );
  sin_angle=cimag/ sqrt(cimag*cimag+ creal*creal) ;
  cos_angle=creal/ sqrt(cimag*cimag+ creal*creal);
  
  if( sin_angle  >= 0. and  cos_angle >=0 ) angle=angle;
  if( sin_angle >= 0. and cos_angle <=0 )  angle=3.14159265 - angle;
  if( sin_angle < 0.  and cos_angle <=0)   angle=3.14159265 + angle;
  if( sin_angle < 0.  and cos_angle >=0)   angle=2*3.14159265-angle;

  cout<<" C_real  = "<<creal<<endl;  
  cout<<" C_imag  = "<<cimag<<endl;
  cout<<" C_amp   =  "<<sqrt(cimag*cimag+ creal*creal)<<endl;
  cout<<" C_phase =  "<< angle <<endl;
}

void Cal_Coupling_Coefficient_Single( Line & linename, const char * filename)
// single element's coupling contribution
{
  int i;
  double creal=0.0, cimag=0.0;
  double angle, sin_angle, cos_angle;
  fstream  f1;

  f1.open(filename, ios::out);
      f1<<setw(15)<<" NAME "<<setw(15)<<" TYPE "
	<<setw(15)<<" creal "<<setw(15)<<" cimag " <<setw(15)<<" C_amp "<<setw(15)<<" C_phase "<<endl;

  for(i=0; i< linename.Ncell;i++) {

    if( linename.Cell[i]->TYPE == string("SKEWQ") ){
      creal= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * M_PI ) / 2. / M_PI  ;
      cimag= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * M_PI ) / 2. / M_PI  ;
      
      angle=  atan( abs(cimag/creal) );
      sin_angle=cimag/ sqrt(cimag*cimag+ creal*creal) ;
      cos_angle=creal/ sqrt(cimag*cimag+ creal*creal);
      if( sin_angle  >= 0. and  cos_angle >=0 ) angle=angle;
      if( sin_angle >= 0. and cos_angle <=0 )  angle=3.14159265 - angle;
      if( sin_angle < 0.  and cos_angle <=0)   angle=3.14159265 + angle;
      if( sin_angle < 0.  and cos_angle >=0)   angle=2*3.14159265-angle;
      f1<<setw(15)<<linename.Cell[i]->NAME<<setw(15)<<setw(15)<<linename.Cell[i]->TYPE
	  <<setw(15)<<creal<<setw(15)<<cimag<<setw(15)<<sqrt(cimag*cimag+ creal*creal)<<setw(15)<<angle<<endl;
    }
    
    if( linename.Cell[i]->TYPE == string("QUAD")  and linename.Cell[i]->DT !=0 ){
      creal= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	cos( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * M_PI ) / 2. / M_PI  ;
      cimag= sqrt(linename.Cell[i]->Beta1 * linename.Cell[i]->Beta2) *
	sin( (linename.Cell[i]->Mu1 - linename.Cell[i]->Mu2)* 2.0 * M_PI ) / 2. / M_PI  ;
      
      angle=  atan( abs(cimag/creal) );
      sin_angle=cimag/ sqrt(cimag*cimag+ creal*creal) ;
      cos_angle=creal/ sqrt(cimag*cimag+ creal*creal);
      if( sin_angle  >= 0. and  cos_angle >=0 ) angle=angle;
      if( sin_angle >= 0. and cos_angle <=0 )  angle=3.14159265 - angle;
      if( sin_angle < 0.  and cos_angle <=0)   angle=3.14159265 + angle;
      if( sin_angle < 0.  and cos_angle >=0)   angle=2*3.14159265-angle;
      f1<<setw(15)<<linename.Cell[i]->NAME<<setw(15)<<setw(15)<<linename.Cell[i]->TYPE
	<<setw(15)<<creal<<setw(15)<<cimag<<setw(15)<<sqrt(cimag*cimag+ creal*creal)<<setw(15)<<angle<<endl;
    }

  }
   f1.close();
}

void Cal_Sext_RDTs( Line & linename )
// sextupole linear RDT at the staring point
{
  int i;
  double h2100_c, h2100_s;
  double h3000_c, h3000_s;
  double h1011_c, h1011_s;  
  double h1002_c, h1002_s;
  double h1020_c, h1020_s;
  double amp, phi, k2l, betx, bety, mux, muy;

  h2100_c=0 ; h2100_s =0 ;
  h3000_c=0 ; h3000_s =0 ;
  h1011_c=0 ; h1011_s =0 ;  
  h1002_c=0 ; h1002_s =0 ;
  h1020_c=0 ; h1020_s =0 ;

  for(i=0;i<linename.Ncell;i++) {
    if( linename.Cell[i]->TYPE==string("SEXT") ||  linename.Cell[i]->TYPE==string("MULT") ){
      k2l=linename.Cell[i]->GetP("K2L");
      betx=linename.Cell[i]->Beta1;
      bety=linename.Cell[i]->Beta2;
      mux =linename.Cell[i]->Mu1 * 2.0 * M_PI;
      muy =linename.Cell[i]->Mu2 * 2.0 * M_PI;
      
      amp = k2l * pow(betx,1.5);      phi = mux;
      h2100_c += amp * cos(phi);      h2100_s += amp * sin(phi);

      amp = k2l * pow(betx,1.5);      phi = 3* mux;
      h3000_c += amp * cos(phi);      h3000_s += amp * sin(phi);

      amp = k2l * pow(betx,0.5) * bety;      phi = mux;
      h1011_c += amp * cos(phi);      h1011_s += amp * sin(phi);

      amp = k2l * pow(betx,0.5) * bety;    phi = mux - 2* muy;
      h1002_c += amp * cos(phi);      h1002_s += amp * sin(phi);

      amp = k2l * pow(betx,0.5) * bety;    phi = mux + 2* muy;
      h1020_c += amp * cos(phi);      h1020_s += amp * sin(phi);
     }
  }

  h2100_c=   h2100_c * (-1.0/8);  h2100_s=   h2100_s * (-1.0/8);
  h3000_c=   h3000_c * (-1.0/24); h3000_s=   h3000_s * (-1.0/24);
  h1011_c=   h1011_c * ( 1.0/ 4); h1011_s=   h1011_s * ( 1.0/ 4);
  h1002_c=   h1002_c * ( 1.0/ 8); h1002_s=   h1002_s * ( 1.0/ 8);
  h1020_c=   h1020_c * ( 1.0/ 8); h1020_s=   h1020_s * ( 1.0/ 8);

  cout<<">>>1st order sextupole geometric resonance driving terms: "<<endl;
  cout<<" h2100 :  "<<h2100_c <<" +i  "<<h2100_s<<" , "<<sqrt(h2100_c*h2100_c + h2100_s*h2100_s)<<" / "<<atan2(h2100_s,h2100_c)<<endl;
  cout<<" h3000 :  "<<h3000_c <<" +i  "<<h3000_s<<" , "<<sqrt(h3000_c*h3000_c + h3000_s*h3000_s)<<" / "<<atan2(h3000_s,h3000_c)<<endl;
  cout<<" h1011 :  "<<h1011_c <<" +i  "<<h1011_s<<" , "<<sqrt(h1011_c*h1011_c + h1011_s*h1011_s)<<" / "<<atan2(h1011_s,h1011_c)<<endl;
  cout<<" h1002 :  "<<h1002_c <<" +i  "<<h1002_s<<" , "<<sqrt(h1002_c*h1002_c + h1002_s*h1002_s)<<" / "<<atan2(h1002_s,h1002_c)<<endl;
  cout<<" h1020 :  "<<h1020_c <<" +i  "<<h1020_s<<" , "<<sqrt(h1020_c*h1020_c + h1020_s*h1020_s)<<" / "<<atan2(h1020_s,h1020_c)<<endl;
}

void Cal_Detuning_Sext( Line & linename )
{
  int i,j;
  double axx=0, axy=0, ayy=0;
  double k2l_1, betx_1, bety_1, mux_1, muy_1;
  double k2l_2, betx_2, bety_2, mux_2, muy_2;
  double dmux, dmuy;
  double Qx0=linename.Tune1, Qy0=linename.Tune2;

  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->TYPE==string("SEXT") ||  linename.Cell[i]->TYPE==string("MULT") ) {
      k2l_1=linename.Cell[i]->GetP("K2L");
      betx_1=linename.Cell[i]->Beta1;
      bety_1=linename.Cell[i]->Beta2;
      mux_1 =linename.Cell[i]->Mu1 * 2.0 * M_PI;
      muy_1 =linename.Cell[i]->Mu2 * 2.0 * M_PI;
      
      for(j=0;j<linename.Ncell;j++) {
	if(linename.Cell[j]->TYPE==string("SEXT") ||  linename.Cell[j]->TYPE==string("MULT") ) {
	  k2l_2=linename.Cell[j]->GetP("K2L");
	  betx_2=linename.Cell[j]->Beta1;
	  bety_2=linename.Cell[j]->Beta2;
	  mux_2 =linename.Cell[j]->Mu1 * 2.0 * M_PI;
	  muy_2 =linename.Cell[j]->Mu2 * 2.0 * M_PI;
	  
	  dmux=abs(   mux_1 -   mux_2 );   dmuy=abs(   muy_1 -   muy_2 );
	  axx += k2l_1 *k2l_2 * pow(betx_1, 1.5) * pow(betx_2, 1.5) 
	    *(      cos(3 *(M_PI*Qx0-dmux))/sin(3*M_PI*Qx0) 
		    + 3* cos(M_PI*Qx0 -dmux )/sin(M_PI*Qx0) );
	  axy += k2l_1 *k2l_2 *sqrt(betx_1*betx_2) * bety_1  
	    *( -bety_2*cos(dmux+2*dmuy-M_PI*(Qx0+2*Qy0))/sin(M_PI*(Qx0+2*Qy0))
	       +bety_2*cos(dmux-2*dmuy-M_PI*(Qx0-2*Qy0))/sin(M_PI*(Qx0-2*Qy0))  
	       +2*betx_2*cos(dmux-M_PI*Qx0)/sin(M_PI*Qx0) );
	  ayy += k2l_1 *k2l_2 *sqrt(betx_1 *betx_2) * bety_1 * bety_2 
	    *(  cos(dmux+2*dmuy-M_PI*(Qx0+2*Qy0)) / sin(M_PI *(Qx0+2*Qy0)) 
		+ cos(dmux-2*dmuy-M_PI*(Qx0-2*Qy0))/sin(M_PI*(Qx0-2*Qy0)) 
		+4 * cos(dmux-M_PI*Qx0) /sin(M_PI*Qx0)  );
	}
      }
    }
  }
  axx=axx*(-1.0/64/M_PI); ayy=ayy*(-1.0/64/M_PI); axy= axy*(1.0/32/M_PI);
  cout<<">>>Amplitude dependent tune shifts from sextupoles:" <<endl;
  cout<<"   2PI * dQx =  axx * (2Jx) - axy *( 2Jy) "<<endl;
  cout<<"   2PI * dQy = -axy * (2Jx) + ayy *( 2Jy) "<<endl;
  cout<<"  "<<endl;
  cout<<"  a_xx = "<<axx<<endl;
  cout<<"  a_xy = "<<axy<<endl;
  cout<<"  a_yy = "<<ayy<<endl;
}

void Cal_Detuning_Oct( Line & linename )
{
  int i;
  double axx=0, axy=0, ayy=0;
  double k3l, betx, bety;

  for(i=0;i<linename.Ncell;i++) {
    if(linename.Cell[i]->TYPE==string("OCT") ||  linename.Cell[i]->TYPE==string("MULT") ) {
      k3l=linename.Cell[i]->GetP("K3L");
      betx=linename.Cell[i]->Beta1;
      bety=linename.Cell[i]->Beta2;
      axx +=  k3l*betx*betx;
      axy -=  2*k3l*betx*bety;
      ayy +=  k3l*bety*bety;
    }
  }
  axx=axx*(1.0/32/M_PI); ayy=ayy*(1.0/32/M_PI); axy= axy*(1.0/32/M_PI);
  cout<<">>>Amplitude dependent tune shifts from octupoles:" <<endl;
  cout<<"   dQx =  axx * (2Jx) - axy *( 2Jy) "<<endl;
  cout<<"   dQy = -axy * (2Jx) + ayy *( 2Jy) "<<endl;
  cout<<"  "<<endl;
  cout<<"  a_xx = "<<axx<<endl;
  cout<<"  a_xy = "<<axy<<endl;
  cout<<"  a_yy = "<<ayy<<endl;
}

//-----------------------------------------
//   path length and gamma-t 
//-----------------------------------------
void Cal_Gammat(Line & linename)
{

  int i,j;
  double x[6];
  double deltap;
  double input1[21],input2[21];
  double coeff[8];
  double alfa0, alfa1, alfa2;
  fstream fout;

  for(i=0;i<21;i++){
    deltap=GP.step_deltap*(i-10);
    Cal_Orbit_Num(linename, deltap);
    x[0]=linename.Cell[linename.Ncell-1]->X[0];
    x[1]=linename.Cell[linename.Ncell-1]->X[1];
    x[2]=linename.Cell[linename.Ncell-1]->X[2];
    x[3]=linename.Cell[linename.Ncell-1]->X[3];
    x[4]=0.000;
    x[5]=deltap;
    for(j=0;j<linename.Ncell; j++) linename.Cell[j]->Pass(x);
    input1[i] = deltap;
    input2[i] = -x[4] / linename.Length  ;
  } 
  
  fout.open("DLL0_vers_deltap.dat", ios::out);
  for (i = 0; i < 21; i++)
    fout << scientific << setw(15) << input1[i]<< scientific << setw(15) << input2[i]<<endl;
  fout.close();

  pfit(input1, input2, 21, coeff, 7);
  alfa0=coeff[1];
  alfa1=coeff[2];
  alfa2=coeff[3]; 
  
  linename.Alfa0=alfa0  ;
  linename.Alfa1=alfa1  ;
  linename.Alfa2=alfa2  ;
  linename.Gammat =  1.0/ sqrt(alfa0);
  linename.Slip=(linename.Alfa0 - 1.0/GP.gamma/GP.gamma);
}

void Cal_Orbit_Length(Line & linename, double deltap)
{
  int j;
  double x[6];
  
  Cal_Orbit_Num(linename, deltap);
  x[0]=linename.Cell[linename.Ncell-1]->X[0];
  x[1]=linename.Cell[linename.Ncell-1]->X[1];
  x[2]=linename.Cell[linename.Ncell-1]->X[2];
  x[3]=linename.Cell[linename.Ncell-1]->X[3];
  x[4]=0.000;
  x[5]=deltap;
  for(j=0;j<linename.Ncell; j++) linename.Cell[j]->Pass(x);
  linename.Orbit_Length =  -x[4] * 3.0e8 +  linename.Length;
}

//-------------------------------
//  longitudinal calculations
//--------------------------------

void Cal_Qs(Line & linename)
{
  int i;
  double Vrf_tot=0,  phi_s =0;
  
  for(i=0;i<linename.Ncell;i++)
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
      Vrf_tot +=linename.Cell[i]->GetP("VRF");
    }
  linename.Vrf_tot= Vrf_tot;
  
  Cal_Gammat(linename);
  linename.Qs= sqrt( GP.harm * GP.Q * Vrf_tot * abs(linename.Slip * cos(phi_s) )/ 2/M_PI/GP.beta/GP.beta/GP.energy/GP.A ) ;
}

void Cal_Bucket_Area(Line & linename)
{
  int i;
  double Vrf_tot=0;
  
  for(i=0;i<linename.Ncell;i++)
    if(linename.Cell[i]->TYPE==string("RFCAV") ){
      Vrf_tot +=linename.Cell[i]->GetP("VRF");
    }

  linename.Bucket_length= linename.Length * 1.0e9 / GP.harm / ( GP.beta * light_speed) ;
  linename.Bucket_height= 2*sqrt(GP.Q*Vrf_tot/2/M_PI/GP.beta/GP.beta/GP.energy/GP.A/GP.harm/abs(linename.Slip)); 
  linename.Bucket_area=16*1.0e6*sqrt( GP.beta*GP.beta* GP.A *GP.energy * GP.Q * Vrf_tot /2/M_PI/(2*M_PI*linename.frev0)/(2*M_PI*linename.frev0)/GP.harm/GP.harm/GP.harm/linename.Slip ) / GP.A ;
}

double RF_F_function(double phi_s, double phi_right, double phi_left)
{
  int i;
  double dphi,phi;
  double sum=0;

  dphi=(phi_right - phi_left) / 1000;
  for(i=0;i<1000;i++){
    phi= phi_left + dphi * (i+1);
    sum += sqrt( abs(cos(phi_left)-cos(phi) + (phi_left-phi) * sin(phi_s)  ) ) * dphi;
  }
  return sum*sqrt(2.0)/8;
}

void Cal_Bunch_Area(Line & linename, double full_length)
// full length is +/-3sigma_l, that is, 6 sigma_l, in units of ns
{
  double phi_s=0, delta_phi, phi_right, phi_left;

  linename.Bunch_length =  full_length ; 
  delta_phi =  full_length * 2 *M_PI / linename.Bucket_length ;
  phi_right =   delta_phi /2  ;
  phi_left =   -delta_phi /2  ;
  linename.Bunch_area  =   RF_F_function( 0, phi_right, phi_left) * linename.Bucket_area ;
  linename.Bunch_height =  sqrt(abs(cos(phi_right)  - cos(phi_s) + ( phi_right - phi_s ) * sin(phi_s) ) )
    *  linename.Bucket_area * GP.harm * ( 2*M_PI*linename.frev0)/ 8/ sqrt(2.0) / (GP.energy*1.0e6)*GP.beta*GP.beta;
}

void Cal_Bunch_Height(Line & linename, double bunch_area)
{
  int i;
  double guess,  bunch_area_0,  bunch_area_1, scale;
  
  Cal_Qs(linename);
  Cal_Bucket_Area(linename);

  i=0;
  guess = linename.Bucket_length*0.5; 
  do{
    Cal_Bunch_Area( linename, guess);
    bunch_area_0  = linename.Bunch_area ;   
    Cal_Bunch_Area( linename, guess+0.05);
    bunch_area_1  = linename.Bunch_area ; 
    scale = 0.05 / ( bunch_area_1 -  bunch_area_0 );
    
    guess = ( bunch_area -  bunch_area_0 ) * scale + guess; 
    Cal_Bunch_Area( linename, guess);
    i++;
  } while (i < 30 && abs(bunch_area - linename.Bunch_area ) > 0.1 );  

}

void Print_Longitudinal_Summary( Line & linename)
{
  cout<<"========================================="<<endl;
  cout<<"Beam energy  = "<<GP.energy<<"  MeV "<<endl;
  cout<<"gamma        = "<<GP.gamma<<endl;
  cout<<"beta         = "<<GP.beta<<endl;
  cout<<"circumference= "<<linename.Length<<"  m   "<<endl;
  cout<<"Revolution frequency = "<< linename.frev0 <<"  Hz  "<<endl;
  cout<<"particle's A : "<<GP.A<<endl;
  cout<<"particle's Q  : "<<GP.Q<<endl;
  cout<<"harmnic number = "<<GP.harm<<endl;
  cout<<"RF freq = "<< linename.frf <<"  Hz "<<endl;
  cout<<"RF total voltage = "<< linename.Vrf_tot <<"  MV "<<endl; 
  cout<<"GammaT  = "<< linename.Gammat<<endl;
  cout<<"Phase slip factor = " << linename.Slip<<endl;
  cout<<"Qs= "<<linename.Qs<<endl;  
  cout<<"bucket length = "<<linename.Bucket_length<<"  ns   "<<endl; 
  cout<<"bucket height  (dp/p0_max) = "<<linename.Bucket_height<<endl; 
  cout<<"bucket area (un-normalized) = "<<linename.Bucket_area << "  eV.s/n  "<<endl;  
  cout<<"bunch length = "<<linename.Bunch_length<<"  ns  "<<endl; 
  cout<<"bunch area (un-normalized) = "<<linename.Bunch_area<<"  eV.s/n   "<<endl; 
  cout<<"bunch height (dp/p0_max) = "<<linename.Bunch_height<<endl; 
}

//------------------------------------------------
//        Print and Plot
//------------------------------------------------
void Print_Twiss(Line & linename, const char* filename)
{
 int i;
 fstream fout;
 
 fout.open(filename, ios::out);
 fout <<setw(15) <<"NAME"<<setw(15) <<"TYPE"<<setw(15) <<"S"
      <<setw(15) <<"X   "<<setw(15)  <<"PX    "<<setw(15)<<"Y    "<<setw(15)<<"PY"
      <<setw(15) <<"BETA1"<<setw(15) <<"ALFA1"<<setw(15) <<"MU1"
      <<setw(15) <<"BETA2"<<setw(15) <<"ALFA2"<<setw(15) <<"MU2"
      <<setw(15) <<"ETAX"<<setw(15)  <<"ETAXP"<<setw(15) <<"ETAY"<<setw(15)  <<"ETAYP"<<endl;
 for(i=0;i<linename.Ncell;i++)
   fout <<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->S
	<<setw(15) <<linename.Cell[i]->X[0]<<setw(15) <<linename.Cell[i]->X[1]<<setw(15) <<linename.Cell[i]->X[2]<<setw(15) <<linename.Cell[i]->X[3]
	<<setw(15) <<linename.Cell[i]->Beta1<<setw(15) <<linename.Cell[i]->Alfa1<<setw(15) <<linename.Cell[i]->Mu1
	<<setw(15) <<linename.Cell[i]->Beta2<<setw(15) <<linename.Cell[i]->Alfa2<<setw(15) <<linename.Cell[i]->Mu2
	<<setw(15) <<linename.Cell[i]->Etax<<setw(15) <<linename.Cell[i]->Etaxp<<setw(15) <<linename.Cell[i]->Etay<<setw(15) <<linename.Cell[i]->Etayp<<endl;
 fout.close();
}

void Print_Twiss_6D(Line & linename, const char* filename)
{
 int i;
 fstream fout;
 
 fout.open(filename, ios::out);
 fout <<setw(15) <<"NAME"<<setw(15) <<"TYPE"<<setw(15) <<"S"
      <<setw(15) <<"X   "<<setw(15)  <<"Y    "<<setw(15)<<"Z    "<<setw(15)<<"Delta"
      <<setw(15) <<"BETA1"<<setw(15) <<"ALFA1"<<setw(15) <<"MU1"
      <<setw(15) <<"BETA2"<<setw(15) <<"ALFA2"<<setw(15) <<"MU2"
      <<setw(15) <<"BETA3"<<setw(15) <<"ALFA3"<<setw(15) <<"MU3"
      <<setw(15) <<"ETAX"<<setw(15)  <<"ETAXP"<<setw(15) <<"ETAY"<<setw(15)  <<"ETAYP"<<endl;
 for(i=0;i<linename.Ncell;i++)
   fout <<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->S
	<<setw(15) <<linename.Cell[i]->X[0]<<setw(15) <<linename.Cell[i]->X[2]<<setw(15) <<linename.Cell[i]->X[4]<<setw(15) <<linename.Cell[i]->X[5]
	<<setw(15) <<linename.Cell[i]->Beta1<<setw(15) <<linename.Cell[i]->Alfa1<<setw(15) <<linename.Cell[i]->Mu1
	<<setw(15) <<linename.Cell[i]->Beta2<<setw(15) <<linename.Cell[i]->Alfa2<<setw(15) <<linename.Cell[i]->Mu2
	<<setw(15) <<linename.Cell[i]->Beta3<<setw(15) <<linename.Cell[i]->Alfa3<<setw(15) <<linename.Cell[i]->Mu3
	<<setw(15) <<linename.Cell[i]->Etax<<setw(15) <<linename.Cell[i]->Etaxp<<setw(15) <<linename.Cell[i]->Etay<<setw(15) <<linename.Cell[i]->Etayp<<endl;
 fout.close();
}

void Print_A_Matrix(Line & linename, const char* filename)
{
  int i,j;
 fstream fout;
 
 fout.open(filename, ios::out);
 for(i=0;i<linename.Ncell;i++) {
   fout<<">>>>"<<setw(15) <<linename.Cell[i]->NAME<<setw(15) <<linename.Cell[i]->TYPE<<setw(15) <<linename.Cell[i]->L<<setw(15) <<linename.Cell[i]->S<<endl;
   for(j=0;j<6;j++)
     fout<<setw(15) <<linename.Cell[i]->A[j*6+0]<<setw(15) <<linename.Cell[i]->A[j*6+1]<<setw(15) <<linename.Cell[i]->A[j*6+2]<<setw(15)
	 <<setw(15) <<linename.Cell[i]->A[j*6+3]<<setw(15) <<linename.Cell[i]->A[j*6+4]<<setw(15) <<linename.Cell[i]->A[j*6+5]<<endl;
 }
 fout.close();
}

void Print_Optics_Summary(Line & linename)
{
  cout<<"-----------------------------------------------------------"<<endl;
  cout<<setw(12)<<" Optics Summary:"<<endl;
  cout<<setw(12)<<" Length:" <<setw(20)<<setprecision(15)<<linename.Length<<endl; 
  cout<<setw(12)<<" Tunes :" <<setw(20)<<setprecision(15)<<linename.Tune1<<setw(20)<<setprecision(15)<<linename.Tune2<<endl;
  cout<<setw(12)<<" Chrom1:" <<setw(20)<<setprecision(15)<<linename.Chromx1<<setw(20)<<setprecision(15)<<linename.Chromy1<<endl;
  cout<<setw(12)<<" Chrom2:" <<setw(20)<<setprecision(15)<<linename.Chromx2<<setw(20)<<setprecision(15)<<linename.Chromy2<<endl;
  cout<<setw(12)<<" Chrom3:" <<setw(20)<<setprecision(15)<<linename.Chromx3<<setw(20)<<setprecision(15)<<linename.Chromy3<<endl;
  cout<<setw(12)<<" Beta* :" <<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Beta1<<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Beta2<<endl;
  cout<<setw(12)<<" Alfa* :" <<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Alfa1<<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Alfa2<<endl;
  cout<<setw(12)<<" Etax* :" <<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etax<<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etaxp<<endl;
  cout<<setw(12)<<" Etay* :" <<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etay<<setw(20)<<setprecision(15)<<linename.Cell[linename.Ncell-1]->Etayp<<endl;
  int  i;
  double betax_max=0, betay_max=0;
  for(i=0;i<linename.Ncell;i++){
    if ( linename.Cell[i]->Beta1 > betax_max ) betax_max = linename.Cell[i]->Beta1 ;
    if ( linename.Cell[i]->Beta2 > betay_max ) betay_max = linename.Cell[i]->Beta2 ;
  }
  cout<<setw(12)<<" Beta_max:" <<setw(20)<<setprecision(15)<<betax_max<<setw(20)<<setprecision(15)<<betay_max<<endl; 
}

void Plot_Twiss(Line & linename, double s1, double s2)
{
 int i;
 fstream fout;
 fout.open("temp_twiss", ios::out);
 for(i=0;i<linename.Ncell;i++)
   if(linename.Cell[i]->S < s2 and linename.Cell[i]->S > s1 )
     fout<<linename.Cell[i]->TYPE<<"  "<<linename.Cell[i]->L<<"  "<<linename.Cell[i]->S<<"    "<<linename.Cell[i]->Beta1<<"    "<<linename.Cell[i]->Beta2<<"  "<<linename.Cell[i]->Etax<<endl;
 fout.close();

 fout.open("temp1.p", ios::out);
 fout<<"set term post color enhanced 20 "<<endl;
 fout<<"set output 'twiss.ps' "<<endl;
 fout<<"set xlabel 's [m]'  " <<endl;
 fout<<"set ylabel '{/Symbol \142}_x , {/Symbol \142}_y [m]  " <<endl;
 fout<<"set y2label 'Dx [m]' "<<endl;
 fout<<"set ytics nomirror"<<endl;
 fout<<"set y2tics"<<endl;
 fout<<"plot 'temp_twiss' u 3:4 tit '{/Symbol \142}_x' w l lt 1  lw 2,\\"<<endl;
 fout<<"    'temp_twiss' u 3:5 tit '{/Symbol \142}_y' w l lt 3 lw 2,\\"<<endl;
 fout<<"    'temp_twiss' u 3:6 axes x1y2 tit 'D_x'  w l lt 2 lw 2"<<endl;
 fout<<"exit"<<endl;
 fout.close();

 system("gnuplot temp1.p");
 system("rm temp1.p");
 system("rm temp_twiss");
}

void Plot_Orbit(Line & linename, double s1, double s2)
{
 int i;
 fstream fout;
 fout.open("temp_orbit", ios::out);
 for(i=0;i<linename.Ncell;i++)
   if(linename.Cell[i]->S < s2 and linename.Cell[i]->S > s1 )
     fout<<linename.Cell[i]->TYPE<<"  "<<linename.Cell[i]->L<<"  "<<linename.Cell[i]->S<<"    "<<linename.Cell[i]->X[0]<<"    "<<linename.Cell[i]->X[1]<<endl;
 fout.close();

 fout.open("temp1.p", ios::out);
 fout<<"set term post color enhanced 20 "<<endl;
 fout<<"set output 'orbit.ps' "<<endl;
 fout<<"set xlabel 's [m]'  " <<endl;
 fout<<"set ylabel 'x_{co}, y_{co}  [m]  " <<endl;
 fout<<"plot 'temp_orbit' u 3:4 tit 'x_{co}' w l lw 2,\\"<<endl;
 fout<<"     'temp_orbit' u 3:5 tit 'y_{co}'  w l lt 2 lw 2"<<endl;
 fout<<"exit"<<endl;
 fout.close();

 system("gnuplot temp1.p");
 system("rm temp1.p");
 system("rm temp_orbit");
}

//---------------------------------------
//   Artifical phase rotator matrix
//--------------------------------------
void Add_Phaser(Line & linename, int loc, const char * name,  double mux, double muy)
{
  int i;
  Element * temp_element;

  double bx, ax, gx,  by, ay, gy, dx, dxp;
  double R11, R12, R21,R22, R33, R34,R43,R44, A,B,C;
  double m[6][6];
  double xco_in[6], xco_out[6];

  for(i=0;i<6;i++) xco_in[i]=0.;
  for(i=0;i<6;i++) xco_out[i]=0.;
  
  bx= linename.Cell[loc-1]->Beta1;
  ax= linename.Cell[loc-1]->Alfa1;
  gx= (1+ax*ax)/bx; 
  
  by= linename.Cell[loc-1]->Beta2;
  ay= linename.Cell[loc-1]->Alfa2;
  gy= (1+ay*ay)/by; 
  
  dx = linename.Cell[loc-1]->Etax;
  dxp= linename.Cell[loc-1]->Etaxp; 
  
  R11=cos(mux) + ax*sin(mux);
  R12=bx*sin(mux);
  R21=-gx*sin(mux);
  R22=cos(mux)-ax*sin(mux);
  
  R33=cos(muy) + ay*sin(muy);
  R34=by*sin(muy);
  R43=-gy*sin(muy);
  R44=cos(muy)-ay*sin(muy);
  
  A=R21*dx-R11*dxp+dxp;
  B=-R12*dxp+R22*dx-dx;
  C=0;
  
  m[0][0]=R11;
  m[0][1]=R12;
  m[0][2]=0;
  m[0][3]=0;
  m[0][4]=0;
  m[0][5]=dx-R11*dx-R12*dxp;
  
  m[1][0]=R21;
  m[1][1]=R22;
  m[1][2]=0;
  m[1][3]=0;
  m[1][4]=0;
  m[1][5]=dxp-R21*dx-R22*dxp;
  
  m[2][0]=0;
  m[2][1]=0;
  m[2][2]=R33;
  m[2][3]=R34;
  m[2][4]=0;
  m[2][5]=0;
  
  m[3][0]=0;
  m[3][1]=0;
  m[3][2]=R43;
  m[3][3]=R44;
  m[3][4]=0;
  m[3][5]=0;
  
  m[4][0]=A;
  m[4][1]=B;
  m[4][2]=0;
  m[4][3]=0;
  m[4][4]=1;
  m[4][5]=C-A*dx-B*dxp;
  
  m[5][0]=0;
  m[5][1]=0;
  m[5][2]=0;
  m[5][3]=0;
  m[5][4]=0;
  m[5][5]=1;
  
  temp_element= new MATRIX(name, 0.0, &m[0][0], xco_in, xco_out);
  linename.Insert(loc, temp_element);
}

//----------------------------------------
//  ORBIT correction
//----------------------------------------
void Correct_Orbit_SVD(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane)
{
  if( n > m ) {
    cout<<" Numer of BPM should be larger than the number of correctors."<<endl;
    exit(0);
  }
  
  int     i,j, k, p;
  double  A[m][n], U[m][m], VT[n][n], AI[n][m];
  double  eps, cut_scale=1000;
  vector<double> temp_vector; 
  vector< vector<double> > V1, U1;
  double read[m], dkick[n];
  
  if( plane ==0 ) {
      for(i=0;i<m;i++) 
	for(j=0;j<n;j++)  A[i][j]=sqrt( linename.Cell[bpm_index[i]]->Beta1 * linename.Cell[kicker_index[j]]->Beta1 ) 
	  *cos( abs(  linename.Cell[bpm_index[i]]->Mu1 - linename.Cell[kicker_index[j]]->Mu1 )*2.0*M_PI - M_PI* linename.Tune1 )
	  /2.0/sin( M_PI * linename.Tune1);
  }
  else{
    for(i=0;i<m;i++) 
      for(j=0;j<n;j++)  A[i][j]=sqrt( linename.Cell[bpm_index[i]]->Beta2 * linename.Cell[kicker_index[j]]->Beta2 ) 
	*cos( abs(  linename.Cell[bpm_index[i]]->Mu2*2*M_PI    - linename.Cell[kicker_index[j]]->Mu2*2*M_PI  ) - M_PI* linename.Tune2 )
	/2.0/sin( M_PI * linename.Tune2);
  }
  
  eps=1.0e-10;
  i=bmuav(&A[0][0],m,n,&U[0][0], &VT[0][0],eps,m+1);
  if( i<= 0) { cout<<" SVD failed"<<endl;    exit(0);  }
  
  for(i=0;i<n;i++) {
    if( A[i][i] < A[0][0]/cut_scale ) {
      p=i;
      break;
    }
  }
  for(i=0;i<m;i++){
    for(j=0;j<p;j++) temp_vector.push_back(U[i][j]);
    U1.push_back(temp_vector);
    temp_vector.clear();
  }
  for(i=0;i<n;i++){
    for(j=0;j<p;j++) temp_vector.push_back(VT[j][i]);
    V1.push_back(temp_vector);
    temp_vector.clear();
  } 
  for(i=0;i<n;i++) 
    for(j=0;j<m;j++) {
      AI[i][j] = 0.0;
      for(k=0;k<p;k++) AI[i][j]=AI[i][j] + V1[i][k]/ A[k][k] * U1[j][k];
    }

  if (plane ==0 ) {
    for(i=0;i<m;i++) read[i]=linename.Cell[ bpm_index[i]]->X[0]; 
    for(i=0;i<n;i++) {
      dkick[i]=0;
      for(j=0;j<m;j++) dkick[i]= dkick[i] - AI[i][j]*read[j]; 
      linename.Cell[ kicker_index[i] ]->SetP( "HKICK", linename.Cell[ kicker_index[i] ]->GetP("HKICK")+ dkick[i] ); 
    }
  }
  else{
    for(i=0;i<m;i++) read[i]=linename.Cell[ bpm_index[i]]->X[2]; 
    for(i=0;i<n;i++) {
      dkick[i]=0;
      for(j=0;j<m;j++) dkick[i]= dkick[i] - AI[i][j]*read[j]; 
      linename.Cell[ kicker_index[i] ]->SetP( "VKICK", linename.Cell[ kicker_index[i] ]->GetP("VKICK")+ dkick[i]); 
    }
  }
}

void Local_Three_Bump(Line linename, int plane,  const char *corr1,   const char *corr2,   const char *corr3, double kick1)
{
  int loc1, loc2, loc3;
  double beta1, beta2, beta3, scale2, scale3, phi_21,phi_31, phi_32,  kick2, kick3;

  loc1=Get_Index(linename,corr1, 1);
  loc2=Get_Index(linename,corr2, 1);
  loc3=Get_Index(linename,corr3, 1);
  if(plane == 1) {
    beta1 =  linename.Cell[loc1]->Beta1;
    beta2 =  linename.Cell[loc2]->Beta1;
    beta3 =  linename.Cell[loc3]->Beta1;
    phi_21 = ( linename.Cell[loc2]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * M_PI ; 
    phi_31 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * M_PI ; 
    phi_32 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc2]->Mu1 )  * 2.0 * M_PI ; 
    scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
    scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
    kick2=scale2*kick1; 
    kick3=scale3*kick1;
    Set_KL(linename,corr1,"HKICK", kick1);
    Set_KL(linename,corr2,"HKICK", kick2);   
    Set_KL(linename,corr3,"HKICK", kick3);
 }
  else if ( plane == 0 ){
    beta1 =  linename.Cell[loc1]->Beta2;
    beta2 =  linename.Cell[loc2]->Beta2;
    beta3 =  linename.Cell[loc3]->Beta2;
    phi_21 = ( linename.Cell[loc2]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * M_PI ; 
    phi_31 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * M_PI ; 
    phi_32 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc2]->Mu2 )  * 2.0 * M_PI ; 
    kick1 = 0.2e-03;
    scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
    scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
    kick2=scale2*kick1; 
    kick3=scale3*kick1;
    Set_KL(linename,corr1,"VKICK", kick1);
    Set_KL(linename,corr2,"VKICK", kick2);   
    Set_KL(linename,corr3,"VKICK", kick3);
 }
  else{
  }
  if(false){
    Cal_Twiss(linename, 0);
    Print_Twiss(linename,"./twiss");
    exit(0);
  }
}

double RMS_Leakage_Orbit( Line linename, int plane, int i1, int i2 )
{
  int i, count;
  double mean, sum, rms;
  
  sum=0.;
  count=0;
  for(i=0;i<linename.Ncell;i++){
    if( (i-i1)*(i-i2) > 0 )  {
      count++;
      sum += linename.Cell[i]->X[plane*2];
    }
  }
  mean=sum/count;
  
  sum=0.;
  count=0;
  for(i=0;i<linename.Ncell;i++) {
    if( (i-i1)*(i-i2) > 0 ){
      count++;
      sum += (linename.Cell[i]->X[plane*2]-mean) *  (linename.Cell[i]->X[plane*2]-mean) ;
    }
  }
  rms=sqrt ( sum / count );
  
  return rms; 
}

void Correct_Orbit_SlidingBump1(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane)
{
  int i,j;
  double  reading , exceed=0.5e-03;
  int     loc1, loc2, loc3, loc_bpm;
  double  s1, s2, s3, phi_31, phi_21, phi_32, dphi_bpm, beta1,  beta2,  beta3, kick1, kick2, kick3, scale2, scale3;
  int     dir;
  double distance;
  
  for(i=0;i<kicker_index.size();i++ ){

    if(i== kicker_index.size()-1 ){
      loc1=kicker_index[kicker_index.size()-1];
      loc2=kicker_index[0];
      loc3=kicker_index[1];
    }
    else if (i== kicker_index.size()-2 ){
      loc1=kicker_index[kicker_index.size()-2];
      loc2=kicker_index[kicker_index.size()-1];
      loc3=kicker_index[0];
    }
    else{
      loc1=kicker_index[i];
      loc2=kicker_index[i+1];
      loc3=kicker_index[i+2];
    }
    s1 = linename.Cell[loc1]->S;
    s2 = linename.Cell[loc2]->S;
    s3 = linename.Cell[loc3]->S;

    if( s2 > s1 and s2 > s3 and s1 > s3  ) {
      loc_bpm=  loc2 ;
      dir=0;
      if ( plane ==0 ) {
	dphi_bpm=( linename.Tune1 - linename.Cell[loc_bpm]->Mu1 +  linename.Cell[loc3]->Mu1 )   * 2 * M_PI;
      }
      else{
	dphi_bpm=( linename.Tune2 - linename.Cell[loc_bpm]->Mu2 +  linename.Cell[loc3]->Mu2 )   * 2 * M_PI;
      }
    }
    else if ( s1> s2 and s1 > s3 and s3 > s2) {
      loc_bpm= loc2 ;
      dir=1;     
      if ( plane ==0 ) {
	dphi_bpm=( linename.Cell[loc_bpm]->Mu1 + linename.Tune1 -  linename.Cell[loc1]->Mu1 )   * 2 * M_PI;
      }
      else{
	dphi_bpm=( linename.Cell[loc_bpm]->Mu2 + linename.Tune2 -  linename.Cell[loc1]->Mu2 )   * 2 * M_PI;
      }
    }
    else if ( s3 > s2 and s2 > s1 ) {
      distance=100;
      loc_bpm=0;
      for(j=0;j<bpm_index.size();j++){
	if( abs( linename.Cell[ bpm_index[j] ]->S - linename.Cell[loc2]->S )  < distance ){
	  distance=  abs( linename.Cell[ bpm_index[j] ]->S - linename.Cell[loc2]->S );
	  loc_bpm   =   bpm_index[j];
	}
      }
      if(  linename.Cell[loc_bpm ]->S > linename.Cell[loc2]->S ) {
	dir=0; 
	if (plane == 0 ) {dphi_bpm= linename.Cell[loc3]->Mu1*2*M_PI - linename.Cell[ loc_bpm ]->Mu1 * 2.0 * M_PI;}
	else{ dphi_bpm= linename.Cell[loc3]->Mu2*2*M_PI - linename.Cell[ loc_bpm ]->Mu2 * 2.0 * M_PI;}
      }
      else{  
	dir=1; 
	if(plane==0 ){ dphi_bpm=  linename.Cell[ loc_bpm ]->Mu1 * 2.0 * M_PI - linename.Cell[loc1]->Mu1*2*M_PI; }
        else{ dphi_bpm=  linename.Cell[ loc_bpm ]->Mu2 * 2.0 * M_PI - linename.Cell[loc1]->Mu2*2*M_PI; }
      }
    }
    else{ cout<<"something wrong in sliding bump."<<endl; exit(1); }

    if ( plane == 0 ) {
      reading= linename.Cell[ loc_bpm ]->X[0];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta1 ; 
	beta2= linename.Cell[loc2]->Beta1 ;  
	beta3= linename.Cell[loc3]->Beta1 ; 
	phi_21 = ( linename.Cell[loc2]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * M_PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * M_PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc2]->Mu1 )  * 2.0 * M_PI ; 
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune1*2*M_PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune1*2*M_PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune1*2*M_PI;
	if( dir ==1 ) {
	  kick1=-reading
	    / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta1 ) / sin( dphi_bpm  );
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick2=scale2*kick1;  kick3=scale3 * kick1;
	}
	else{
	  kick3=-reading 
	    / sqrt(beta3* linename.Cell[ loc_bpm ]->Beta1 ) / sin( dphi_bpm ); 
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick1=kick3/scale3; kick2=scale2*kick1;
	}
	linename.Cell[loc1]->SetP("HKICK", linename.Cell[loc1]->GetP("HKICK")+ kick1);
	linename.Cell[loc2]->SetP("HKICK", linename.Cell[loc2]->GetP("HKICK")+ kick2);
	linename.Cell[loc3]->SetP("HKICK", linename.Cell[loc3]->GetP("HKICK")+ kick3);
      }
    }
    else{
      reading= linename.Cell[ loc_bpm ]->X[2];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta2 ; 
	beta2= linename.Cell[loc2]->Beta2 ;  
	beta3= linename.Cell[loc3]->Beta2 ; 
	phi_21 = ( linename.Cell[loc2]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * M_PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * M_PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc2]->Mu2 )  * 2.0 * M_PI ; 
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune2*2*M_PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune2*2*M_PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune2*2*M_PI;
	if( dir ==1 ) {
	  kick1=-reading
	    / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta2 ) / sin( dphi_bpm  );
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick2=scale2*kick1;  kick3=scale3 * kick1;
	}
	else{
	  kick3=-reading 
	    / sqrt(beta3* linename.Cell[ loc_bpm ]->Beta2 ) / sin( dphi_bpm ); 
	  scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	  scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	  kick1=kick3/scale3; kick2=scale2*kick1;
	}
	linename.Cell[loc1]->SetP("VKICK", linename.Cell[loc1]->GetP("VKICK")+ kick1);
	linename.Cell[loc2]->SetP("VKICK", linename.Cell[loc2]->GetP("VKICK")+ kick2);
	linename.Cell[loc3]->SetP("VKICK", linename.Cell[loc3]->GetP("VKICK")+ kick3);
      }
    }
  }
}

void Correct_Orbit_SlidingBump2(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane)
{
  int i;
  double  reading , exceed=0.5e-03;
  int     loc1, loc2, loc3, loc_bpm;
  double  phi_31, phi_21, phi_32, dphi_bpm, beta1,  beta2,  beta3, kick1, kick2, kick3, scale2, scale3;
  
  for(i=0;i<kicker_index.size();i++ ){

    if(i== kicker_index.size()-1 ){
      loc1=kicker_index[kicker_index.size()-1];
      loc2=kicker_index[0];
      loc3=kicker_index[1];
    }
    else if (i== kicker_index.size()-2 ){
      loc1=kicker_index[kicker_index.size()-2];
      loc2=kicker_index[kicker_index.size()-1];
      loc3=kicker_index[0];
    }
    else{
      loc1=kicker_index[i];
      loc2=kicker_index[i+1];
      loc3=kicker_index[i+2];
    }
    
    loc_bpm= loc2;

    if ( plane == 0 ) {
      reading= linename.Cell[ loc_bpm ]->X[0];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta1 ; 
	beta2= linename.Cell[loc2]->Beta1 ;  
	beta3= linename.Cell[loc3]->Beta1 ; 
	phi_21 = ( linename.Cell[loc2]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * M_PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc1]->Mu1 )  * 2.0 * M_PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu1 - linename.Cell[loc2]->Mu1 )  * 2.0 * M_PI ; 
	dphi_bpm=( linename.Cell[loc_bpm]->Mu1 -  linename.Cell[loc1]->Mu1 )   * 2 * M_PI;
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune1*2*M_PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune1*2*M_PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune1*2*M_PI;
	if ( dphi_bpm < 0 ) dphi_bpm = dphi_bpm + linename.Tune1*2*M_PI;
	kick1=-reading
	  / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta1 ) / sin( dphi_bpm  );
	scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	kick2=scale2*kick1;  kick3=scale3 * kick1;
	linename.Cell[loc1]->SetP("HKICK", linename.Cell[loc1]->GetP("HKICK")+ kick1);
	linename.Cell[loc2]->SetP("HKICK", linename.Cell[loc2]->GetP("HKICK")+ kick2);
	linename.Cell[loc3]->SetP("HKICK", linename.Cell[loc3]->GetP("HKICK")+ kick3);
      }
    }
    else{
      reading= linename.Cell[ loc_bpm ]->X[2];
      if(abs(reading) > exceed ) { 
	beta1= linename.Cell[loc1]->Beta2 ; 
	beta2= linename.Cell[loc2]->Beta2 ;  
	beta3= linename.Cell[loc3]->Beta2 ; 
	phi_21 = ( linename.Cell[loc2]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * M_PI ; 
	phi_31 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc1]->Mu2 )  * 2.0 * M_PI ; 
	phi_32 = ( linename.Cell[loc3]->Mu2 - linename.Cell[loc2]->Mu2 )  * 2.0 * M_PI ; 
	dphi_bpm=( linename.Cell[loc_bpm]->Mu2  -  linename.Cell[loc1]->Mu2 )   * 2 * M_PI;
        if ( phi_21 < 0 ) phi_21 = phi_21 + linename.Tune2*2*M_PI;
        if ( phi_31 < 0 ) phi_31 = phi_31 + linename.Tune2*2*M_PI;
        if ( phi_32 < 0 ) phi_32 = phi_32 + linename.Tune2*2*M_PI;
	if ( dphi_bpm < 0 ) dphi_bpm = dphi_bpm + linename.Tune2*2*M_PI;
	kick1=-reading
	  / sqrt(beta1* linename.Cell[ loc_bpm ]->Beta2 ) / sin( dphi_bpm  );
	scale2=-sqrt(beta1/beta2)*sin(phi_31)/sin(phi_32);
	scale3=-sqrt(beta1/beta3)*sin(phi_21)/sin(-phi_32); 
	kick2=scale2*kick1;  kick3=scale3 * kick1;
	linename.Cell[loc1]->SetP("VKICK", linename.Cell[loc1]->GetP("VKICK")+ kick1);
	linename.Cell[loc2]->SetP("VKICK", linename.Cell[loc2]->GetP("VKICK")+ kick2);
	linename.Cell[loc3]->SetP("VKICK", linename.Cell[loc3]->GetP("VKICK")+ kick3);
      }
    }
  }
}

void Orbit_Status( Line linename, vector<int> bpm_index, int plane, double &orbit_mean, double & orbit_rms )
{
  int i;
  double sum;
  sum=0.;
  for(i=0;i<bpm_index.size();i++) sum += linename.Cell[ bpm_index[i] ]->X[plane*2];
  orbit_mean=sum/bpm_index.size();

  sum=0.;
  for(i=0;i<bpm_index.size();i++) sum += (linename.Cell[ bpm_index[i] ]->X[plane*2]-orbit_mean) *  (linename.Cell[ bpm_index[i] ]->X[plane*2]-orbit_mean) ;
  orbit_rms=sqrt ( sum / bpm_index.size() );
}

//--------------------------------------------------
//      Nonlinear tools: FOOTPRINT, FMA
//--------------------------------------------------
void Track_tbt_FMA( Line & linename, double deltap0, double sigmax0, double sigmay0 ) // I prefer RF off for this 
{
  int i,j,k;
  int Nturn=2048;
  double sigma_step=0.05;
  int sigma0=0, sigma1=6,  nsigma= int((sigma1-sigma0 ) *1.0 / sigma_step  )+1  ;  

  double x[6];
  double  x_tbt[Nturn*6];
  int stable = 1, lost_turn=0, lost_post=0;
  
  double xtemp[1024];
  double tunex1, tunex2, tuney1, tuney2;

  fstream  f1;
  f1.open("./FMA-output.dat", ios::out);
  f1.close();
  
  Cal_Orbit_Num(linename, deltap0);

  for(i=0; i<nsigma;  i++) {
    for(j=0; j<nsigma; j++ ) {
      if(  sqrt( 1.0*i*i +1.0*j*j) <= (6./sigma_step) ) {
        
        x[0] = i*sigma_step*sigmax0;
        x[1] = 0.;
        x[2] = j*sigma_step*sigmay0;
        x[3] = 0.;
        x[4] = 0.;
        x[5] = deltap0;
	x[0]= x[0] + linename.Cell[linename.Ncell-1]->X[0]  ;
	x[1]= x[1] + linename.Cell[linename.Ncell-1]->X[1]  ;
	x[2]= x[2] + linename.Cell[linename.Ncell-1]->X[2]  ;
	x[3]= x[3] + linename.Cell[linename.Ncell-1]->X[3]  ; 
        
        if(x[0] < 1.0e-09 ) x[0]=x[0]+1.0e-09;
        if(x[2] < 1.0e-09 ) x[2]=x[2]+1.0e-09; 

        stable = 1; lost_turn=0; lost_post=0;
        Track_tbt(linename, x, Nturn, x_tbt, stable, lost_turn, lost_post);
        if(stable==1) {
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+0];
          FineTuneFinder(1024,xtemp,tunex1);
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+2];
          FineTuneFinder(1024,xtemp,tuney1);
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+0];
          FineTuneFinder(1024,xtemp,tunex2);
          for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+2];
          FineTuneFinder(1024,xtemp,tuney2);
          f1.open("./FMA-output.dat",  ios::out | ios::app);
          f1<<i*sigma_step<<" "<<j*sigma_step<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tunex1<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tunex2<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tuney1<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tuney2<<"  "
            <<setw(25)<<setprecision(18)<<scientific<<tunex1-tunex2<<" "
            <<setw(25)<<setprecision(18)<<scientific<<tuney1-tuney2<<endl;
          f1.close();
        }
       } 
    } 
  }
}

void Track_tbt_tune_footprint( Line & linename, double deltap0, double emitx, double emity)  // I prefer RF off for this 
{
  double betax0, alfax0, betay0, alfay0;

  int i,j,k;
  double nsigma, nsigmax, nsigmay;
  double xn, pxn, yn, pyn, x[6];

  int Nturn=2048;
  int stable = 1, lost_turn=0, lost_post=0;
  double  x_tbt[Nturn*6];
  double xtemp[1024];
  double tunex1, tunex2, tuney1, tuney2;

  fstream  f1;
  f1.open("./footprint-output.dat", ios::out);
  f1.close();

  Cal_Twiss(linename, deltap0);
  betax0= linename.Cell[linename.Ncell-1]->Beta1; 
  betay0= linename.Cell[linename.Ncell-1]->Beta2;
  alfax0= linename.Cell[linename.Ncell-1]->Alfa1; 
  alfay0= linename.Cell[linename.Ncell-1]->Alfa2;

  Cal_Orbit_Num(linename,deltap0);

  for( i=0;i<5;i++)
    for(j=0;j<8;j++){

      nsigma=  0.+j*(7-0.)/7 ;
      nsigmax= nsigma * cos((i+1)*15 *M_PI/180);
      nsigmay= nsigma * sin((i+1)*15 *M_PI/180);

      xn=  sqrt(emitx) * nsigmax ; 
      pxn= 0.;
      yn=  sqrt(emity) * nsigmay ; 
      pyn= 0;
      
      x[0]= sqrt(betax0) * xn ;                                          
      x[1]=-alfax0 * xn / sqrt( betax0 ) + 1.0* pxn /  sqrt( betax0 );   
      x[2]= sqrt(betay0) * yn ;                                          
      x[3]=-alfay0 * yn / sqrt( betay0 ) + 1.0* pyn /  sqrt( betay0 );   
      x[4]= 0. ;
      x[5]=  deltap0;
      x[0]= x[0] + linename.Cell[linename.Ncell-1]->X[0]  ;
      x[1]= x[1] + linename.Cell[linename.Ncell-1]->X[1]  ;
      x[2]= x[2] + linename.Cell[linename.Ncell-1]->X[2]  ;
      x[3]= x[3] + linename.Cell[linename.Ncell-1]->X[3]  ; 

      if(x[0] == 0. ) x[0]=x[0]+1.0e-09;
      if(x[2] == 0. ) x[2]=x[2]+1.0e-09; 
      
      stable = 1; lost_turn=0; lost_post=0;
      Track_tbt(linename, x, Nturn, x_tbt, stable, lost_turn, lost_post);
      if(stable==1) {
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+0];
	FineTuneFinder(1024,xtemp,tunex1);
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[k*6+2];
	FineTuneFinder(1024,xtemp,tuney1);
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+0];
	FineTuneFinder(1024,xtemp,tunex2);
	for(k=0;k<1024;k++) xtemp[k]=x_tbt[(k+1023)*6+2];
	FineTuneFinder(1024,xtemp,tuney2);
	f1.open("./footprint-output.dat",  ios::out | ios::app);
	f1<<(i+1)*15<<" "<<j<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tunex1<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tunex2<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tuney1<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tuney2<<"  "
	  <<setw(25)<<setprecision(18)<<scientific<<tunex1-tunex2<<" "
	  <<setw(25)<<setprecision(18)<<scientific<<tuney1-tuney2<<endl;
	f1.close();
      }
    }
}

