//date:02/16
//load a quadrupole setting file and track!
//todo:	beam class: type(r0, charge, mass), energy, emittance, density
//		accelerator class
//		turn-by-turn: energy, Ksc

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "elepass.h"
#include "physics.h"
#include "FODO4.h"
#include "myroutine1.h"
#include "myroutine2.h"
#include "armadillo"
#include "Faddeeva.hh"
#include <omp.h>

using namespace std;
using namespace arma;

int main(int argc, char *argv[])
{
    //input beam parameters
    double epsx;
    stringstream(argv[2]) >> epsx; //rms emittance
    double epsz;
    stringstream(argv[3]) >> epsz;
    double N0;
    stringstream(argv[4]) >> N0; //line density
    int N_particle;
    stringstream(argv[5]) >> N_particle;
    int N_TURN;
    stringstream(argv[6]) >> N_TURN;
    int N_INJTURN;
    stringstream(argv[7]) >> N_INJTURN;

    //input ring parameters: aperture, quadrupoles
    double apx;
    stringstream(argv[8]) >> apx;
    double apz;
    stringstream(argv[9]) >> apz;
    string Qsetting0;
    Qsetting0=argv[10];
    string Qsetting1;
    Qsetting1=argv[11];

    //output control
    int TBT_OUT;
    stringstream(argv[12]) >> TBT_OUT;
    int RAW_OUT;
    stringstream(argv[13]) >> RAW_OUT;
    int SBL_OUT;
    stringstream(argv[14]) >> SBL_OUT;
    int EBE_OUT;
    stringstream(argv[15]) >> EBE_OUT;
    int DNU_OUT;
    stringstream(argv[16]) >> DNU_OUT;
    int ENV_OUT;
    stringstream(argv[17]) >> ENV_OUT;

    //output files
    ofstream fout1,fout2,fout3,fout4,fout5,fout6;
    if( TBT_OUT==1 ){
        stringstream outfile1;
        outfile1<<"tbt_"<<argv[1]<<".dat";
        fout1.open(outfile1.str().c_str());
        fout1<<"#turn count sigx2 sigxp2 sigz2 sigzp2 ex ez sigxp sigzp"<<endl;
    }
    if( RAW_OUT==1 ){
        stringstream outfile2;
        outfile2<<"raw_"<<argv[1]<<".dat";
        fout2.open(outfile2.str().c_str());
        fout2<<scientific;
    }
    if( SBL_OUT==1 ){
        stringstream outfile3;
        outfile3<<"sbl_"<<argv[1]<<".dat";
        fout3.open(outfile3.str().c_str());
    }
    if( EBE_OUT==1 ){
        stringstream outfile4;
        outfile4<<"ebe_"<<argv[1]<<".dat";
        fout4.open(outfile4.str().c_str());
    }
    if( DNU_OUT==1 ){
        stringstream outfile5;
        outfile5<<"dnu_"<<argv[1]<<".dat";
        fout5.open(outfile5.str().c_str());
    }
    if( ENV_OUT==1 ){
        stringstream outfile6;
        outfile6<<"env_"<<argv[1]<<".dat";
        fout6.open(outfile6.str().c_str());
    }

    //Define the ring
    Line FODO;
    line_def(FODO);
    //need these for paticle initialization
    //load_2Q(FODO,Qsetting0);
    //Fit_Tune(FODO, 122.0/360.0, 0.32, "QF", "QD");
    Cal_Twiss(FODO,0.0);//it is necessary for the 0th turn
    double betax0, betaz0, alfax0, alfaz0;
    betax0 = FODO.Cell[0]->Beta1;
    betaz0 = FODO.Cell[0]->Beta2;
    alfax0 = FODO.Cell[0]->Alfa1;
    alfaz0 = FODO.Cell[0]->Alfa2;
    cout<<"beta_x0 = "<<betax0<<", beta_z0 = "<<betaz0<<", alfa_x0 = "<<alfax0<<", alfa_z0 = "<<alfaz0<<endl;

    //load quad error after Cal_Twiss. Don't update Twiss
    //load_2Q(FODO,Qsetting1);
    //Cal_Twiss(FODO,0.0);//it is necessary for the 0th turn
    cout<<Qsetting1<<": (nux,nuz)=("<<FODO.Tune1<<","<<FODO.Tune2<<")"<<endl;

    //Assign the beam type and energy
    double r0=1.535e-18;//re=2.82e-15; rp=1.535e-18;
    double g=(7.+938)/938;
    double b=sqrt(1-1./g/g); //=0.122

    //Define the "path"
    double N[N_TURN+1],Ksc[N_TURN+1],Ksc1[N_TURN+1];
    for(int j=0;j<=N_TURN;j++) {
        if(j<N_INJTURN){
            N[j]=double(j)/N_INJTURN*N0;
        }else{
            N[j]=double(1)*N0;
        }
        Ksc[j]=2*N[j]*r0/b/b/g/g/g; //space charge perveance
    }

    // set and check aperture
    for(int k=0;k<FODO.Ncell;k++)
    {
        FODO.Cell[k]->APx=apx;
        FODO.Cell[k]->APy=apz;
        //cout<<FODO.Cell[k]->APx<<" "<<FODO.Cell[k]->APy<<endl;
    }

    //variables for Laslett tune shift
    double Dnu_x=0, Dnu_z=0;

    //variables related to space charge
    double sx2, sz2, xold, zold, fscx, fscz;
    double betax, betaz, alfax, alfaz;

    //Particles Initialization
    double *x = new double[N_particle*8];
    unsigned int *stable = new unsigned int[N_particle];
    std::fill_n(stable,N_particle,1);
    unsigned int *lost_turn = new unsigned int[N_particle];
    unsigned int *lost_post = new unsigned int[N_particle];
    double phix,phiz,Ax,Az;
    for (unsigned int i=0; i < N_particle; i++) {
        phix=unifRand(-M_PI,M_PI);
        Ax=sqrt(betax0*epsx)*trunc_xgaussRand(4);
        x[i*8+0] = Ax*cos(phix);
        x[i*8+6] = Ax*sin(phix); //Px
        phiz=unifRand(-M_PI,M_PI);
        Az=sqrt(betaz0*epsz)*trunc_xgaussRand(4);
        x[i*8+2] = Az*cos(phiz);
        x[i*8+7] = Az*sin(phiz); //Pz
        x[i*8+4] = 0;
        x[i*8+5] = 0;
    }
    for (unsigned int i=0; i < N_particle; i++) {
        x[i*8+1] = (x[i*8+6]-alfax0*x[i*8+0])/betax0;
        x[i*8+3] = (x[i*8+7]-alfaz0*x[i*8+2])/betaz0;
    }
    mat X(x, 8, N_particle, false);
    uvec STABLE(stable, N_particle, false);

    //initial statistics
    mat sigx2  = cov(X.row(0),1);
    mat sigxp2 = cov(X.row(1),1);
    mat sigz2  = cov(X.row(2),1);
    mat sigzp2 = cov(X.row(3),1);
    mat sigxxp = cov(X.row(0),X.row(1),1);
    mat sigzzp = cov(X.row(2),X.row(3),1);
    mat sigxp, sigzp; //d(sig_x)/ds, derivative of sig_x
    mat e1 = sqrt(sigx2*sigxp2-sigxxp*sigxxp);
    mat e2 = sqrt(sigz2*sigzp2-sigzzp*sigzzp);
    unsigned int count = N_particle;
    unsigned int countk = 0;

    //define index vector
    uvec row0; row0<<0;
    uvec row1; row1<<1;
    uvec row2; row2<<2;
    uvec row3; row3<<3;
    uvec row06; row06<<0<<6;
    uvec row27; row27<<2<<7;
    uvec row0627; row0627<<0<<6<<2<<7;
    uvec row0123; row0123<<0<<1<<2<<3;
    uvec live_index;
    //uvec dead_index;

    //zero mean
    mat Xbar=zeros(8,N_particle);
    Xbar.row(0)=mean(X.row(0))*ones(1,N_particle);
    Xbar.row(1)=mean(X.row(1))*ones(1,N_particle);
    Xbar.row(2)=mean(X.row(2))*ones(1,N_particle);
    Xbar.row(3)=mean(X.row(3))*ones(1,N_particle);
    X=X-Xbar;

    //emittance vector for smoothization
    vector<double> e1vector, e2vector;
    double aveg_e1=e1(0),aveg_e2=e2(0);
    for(unsigned int k=0;k<FODO.Ncell;k++) {
        if (FODO.Cell[k]->NAME==string("SPKICK")){
            e1vector.push_back(e1(0));
            e2vector.push_back(e2(0));
        }
    }


    // test particles for poincare surface
    const int Ntest_action=40, Ntest_Angle=10, Ntest=Ntest_action*Ntest_Angle;
    mat test_Xbar=zeros(8,Ntest);
    double *test_x = new double[Ntest*8];
    mat test_X(test_x, 8, Ntest, false);
    for (unsigned int i=0; i<Ntest_action; i++) {
        for (unsigned int j=0; j <Ntest_Angle; j++) {
            phix=2.0*j*PI/(double)Ntest_Angle;
            Ax=sqrt(betax0*epsx)*0.15*(i+1);
            test_x[(i*Ntest_Angle+j)*8+0] = Ax*cos(phix);
            test_x[(i*Ntest_Angle+j)*8+6] = Ax*sin(phix); //Px
            test_x[(i*Ntest_Angle+j)*8+2] = 0;
            test_x[(i*Ntest_Angle+j)*8+7] = 0;
            test_x[(i*Ntest_Angle+j)*8+4] = 0;
            test_x[(i*Ntest_Angle+j)*8+5] = 0;
        }
    }
    for (unsigned int i=0; i < Ntest; i++) {
        test_x[i*8+1] = (test_x[i*8+6]-alfax0*test_x[i*8+0])/betax0;
        test_x[i*8+3] = (test_x[i*8+7]-alfaz0*test_x[i*8+2])/betaz0;
    }
//    test_Xbar.row(0)=mean(X.row(0))*ones(1,Ntest);
//    test_Xbar.row(1)=mean(X.row(1))*ones(1,Ntest);
//    test_Xbar.row(2)=mean(X.row(2))*ones(1,Ntest);
//    test_Xbar.row(3)=mean(X.row(3))*ones(1,Ntest);
//    test_X=test_X-test_Xbar;

    //Main Loop
    for(unsigned int j=0;j<=N_TURN;j++) {
        //update statistics data for every turn
        live_index = find(STABLE != 0);
        sigx2  = cov(X(row0, live_index),1);
        sigxp2 = cov(X(row1, live_index),1);
        sigz2  = cov(X(row2, live_index),1);
        sigzp2 = cov(X(row3, live_index),1);
        sigxxp = cov(X(row0, live_index),X(row1, live_index),1);
        sigzzp = cov(X(row2, live_index),X(row3, live_index),1);
        e1 = sqrt(sigx2*sigxp2-sigxxp*sigxxp);
        e2 = sqrt(sigz2*sigzp2-sigzzp*sigzzp);

        sigxp = mean(X(row0, live_index) % X(row1, live_index),1);
        sigxp -= mean(X(row0, live_index),1)*mean(X(row1, live_index),1); //this term is supposed to be zero
        sigxp /= sqrt(sigx2);
        sigzp = mean(X(row2, live_index) % X(row3, live_index),1);
        sigzp -= mean(X(row2, live_index),1)*mean(X(row3, live_index),1); //this term is supposed to be zero too
        sigzp /= sqrt(sigz2);

        if(fout1.is_open()){
            fout1<<j<<" "<<count<<" "<<sigx2(0)<<" "<<sigxp2(0)<<" "<<sigz2(0)<<" "<<sigzp2(0)<<" "<<e1(0)<<" "<<e2(0)<<" "<<sigxp(0)<<" "<<sigzp(0)<<" "<<endl; //add (0) to avoid newline
        }
        if(fout2.is_open()){
            for (int i=0; i < Ntest; i++) {
                test_x[i*8+6] = betax0*test_x[i*8+1]+alfax0*test_x[i*8+0];
                test_x[i*8+7] = betaz0*test_x[i*8+3]+alfaz0*test_x[i*8+2];
            }
            for (unsigned int i=0; i < Ntest; i++) {
		if(j>50)
                fout2<<test_x[i*8+0]<<" "<<test_x[i*8+1]<<" "<<test_x[i*8+6]<<endl;
            }
            //			fout2<<"#n="<<j<<'F'<<endl;
            //			fout2<<trans(X.cols(live_index))<<endl<<endl;
        }
        if(fout3.is_open()){
            fout3<<j<<" "<<STABLE.t();
            //fout3<<j<<" "<<trans(STABLE.elem(live_index));
        }

        count=0;	//counting number of survival particles
        countk=0;	//count sp_kickers number, for writing envelope at every six sp_kickers
        for(unsigned int k=0;k<FODO.Ncell;k++) { //for every element in one turn
            //#pragma omp parallel
	    //{  
            //#pragma omp parallel for 
            for (unsigned int i=0; i < Ntest; i++) { FODO.Cell[k]->Pass(test_x+i*8);} // test
            //#pragma omp parallel for 
                for (unsigned int i=0; i < N_particle; i++) { //for every particles
                    if(stable[i]!=0){
                        FODO.Cell[k]->Pass(x+i*8);
                        if( abs(x[i*8+0]) > FODO.Cell[k]->APx or abs(x[i*8+2]) > FODO.Cell[k]->APy or isnan(x[i*8+0])==1){
                            //cout<<i<<"-th particle lost in "<<k<<"-th elements. "<<x[i*8+0]<<" "<<x[i*8+2]<<" "<<fscx<<" "<<fscz<<endl;
                            stable[i]=0; lost_turn[i]=GP.turn; lost_post[i]=k;
                        }
                    }
                }
            //}
            live_index = find(STABLE != 0);
            count = live_index.size();
            Ksc1[j]=Ksc[j]*(double)count/N_particle;

            if(fout5.is_open()){
                sigx2  = cov(X(row0, live_index),1);
                sigz2  = cov(X(row2, live_index),1);
                Dnu_x -= Ksc1[j]/M_PI/4.0*FODO.Cell[k]->L*FODO.Cell[k]->Beta1/(sqrt(sigx2(0))*(sqrt(sigx2(0))+sqrt(sigz2(0))));
                Dnu_z -= Ksc1[j]/M_PI/4.0*FODO.Cell[k]->L*FODO.Cell[k]->Beta2/(sqrt(sigz2(0))*(sqrt(sigx2(0))+sqrt(sigz2(0))));
            }

            //space charge kick
            if (FODO.Cell[k]->NAME==string("SPKICK")){
                countk++;
                //zero mean
                Xbar.row(0)=mean(X(row0,live_index),1)*ones(1,N_particle);
                Xbar.row(1)=mean(X(row1,live_index),1)*ones(1,N_particle);
                Xbar.row(2)=mean(X(row2,live_index),1)*ones(1,N_particle);
                Xbar.row(3)=mean(X(row3,live_index),1)*ones(1,N_particle);
                X=X-Xbar;

                betax = FODO.Cell[k]->Beta1;//should MODIFY it by fitting!
                betaz = FODO.Cell[k]->Beta2;

                sigx2  = cov(X(row0, live_index),1);
                sigxp2 = cov(X(row1, live_index),1);
                sigz2  = cov(X(row2, live_index),1);
                sigzp2 = cov(X(row3, live_index),1);
                sigxxp = cov(X(row0, live_index),X(row1, live_index),1);
                sigzzp = cov(X(row2, live_index),X(row3, live_index),1);
                e1 = sqrt(sigx2*sigxp2-sigxxp*sigxxp);
                e2 = sqrt(sigz2*sigzp2-sigzzp*sigzzp);

                //cout<<j<<" "<<FODO.Cell[k]->S<<" "<<e1(0)<<" "<<e2(0)<<endl;
                e1vector.push_back(e1(0));
                e2vector.push_back(e2(0));
                e1vector.erase (e1vector.begin());
                e2vector.erase (e2vector.begin());

                aveg_e1 = 0;
                aveg_e2 = 0;
                //for (it=e1vector.begin(); it<e1vector.end(); it++)
                for (vector<double>::iterator it = e1vector.begin(); it!=e1vector.end(); ++it){
                    aveg_e1 += *it;
                }
                //for (it=e2vector.begin(); it<e2vector.end(); it++)
                for (vector<double>::iterator it = e2vector.begin(); it!=e2vector.end(); ++it){
                    aveg_e2 += *it;
                }
                aveg_e1 /= e1vector.size();
                aveg_e2 /= e2vector.size();

                //smoothize the emittance by the past 108 space charge kickers
                sx2=betax*aveg_e1;//sx2=sigx2(0);
                sz2=betaz*aveg_e2;//sz2=sigz2(0);
                //SC kick
                //#pragma omp parallel for 
                    for (unsigned int i=0; i < N_particle; i++) { //for every particles
                        if(stable[i]!=0){
                            xold=x[i*8+0];
                            zold=x[i*8+2];
                            if (FODO.Cell[k]->NAME==string("SPKICK")){
                                Fsc2(Ksc1[j],0.5,sx2,sz2,xold,zold,fscx,fscz);
                            }
                            x[i*8+1]+=fscx;
                            x[i*8+3]+=fscz;
                        }
                    }
                //test particle SC kick
                test_Xbar.row(0)=mean(X.row(0))*ones(1,Ntest);
                test_Xbar.row(1)=mean(X.row(1))*ones(1,Ntest);
                test_Xbar.row(2)=mean(X.row(2))*ones(1,Ntest);
                test_Xbar.row(3)=mean(X.row(3))*ones(1,Ntest);
                test_X=test_X-test_Xbar;
                for (unsigned int i=0; i < Ntest; i++) { //for every test particles
                        xold=test_x[i*8+0];
                        zold=test_x[i*8+2];
                        if (FODO.Cell[k]->NAME==string("SPKICK")){
                            Fsc2(Ksc1[j],0.5,sx2,sz2,xold,zold,fscx,fscz);
                        }
                        test_x[i*8+1]+=fscx;
                        test_x[i*8+3]+=fscz; 
                }

                if(fout4.is_open()){
                    alfax = FODO.Cell[k]->Alfa1;
                    alfaz = FODO.Cell[k]->Alfa2;
                    for (unsigned int i=0; i < N_particle; i++) {
                        x[i*8+6] = alfax*x[i*8+0]+betax*x[i*8+1];
                        x[i*8+7] = alfaz*x[i*8+2]+betaz*x[i*8+3];
                    }
                    //fout4<<"s = "<<FODO.Cell[k]->S<<endl;
                    //fout4<<trans(X.rows(row0627))<<endl<<endl;
                    //fout4<<trans(X.submat(row0627,live_index))<<endl<<endl;
                    fout4<<trans(join_cols(X.rows(row06)/sqrt(betax),X.rows(row27)/sqrt(betaz)))<<endl<<endl;
                }
                if(fout6.is_open()){
                    if(countk%6==0){
                        sigx2  = cov(X(row0, live_index),1);
                        sigz2  = cov(X(row2, live_index),1);
                        sigxp = mean(X(row0, live_index) % X(row1, live_index),1);
                        sigxp -= mean(X(row0, live_index),1)*mean(X(row1, live_index),1); //this term is supposed to be zero
                        sigxp /= sqrt(sigx2);
                        sigzp = mean(X(row2, live_index) % X(row3, live_index),1);
                        sigzp -= mean(X(row2, live_index),1)*mean(X(row3, live_index),1); //this term is supposed to be zero too
                        sigzp /= sqrt(sigz2);
                        fout6 <<sqrt(sigx2(0))<<" "<<sigxp(0)<<" "<<sqrt(sigz2(0))<<" "<<sigzp(0)<<endl;
                    }
                }

            }else if (FODO.Cell[k]->NAME==string("END")){
                if(fout5.is_open()){
                    fout5<<Dnu_x<<" "<<Dnu_z<<endl;
                }
                Dnu_x=0;
                Dnu_z=0;
            }
        } //end of elements
        if(j%100==0) cout<<"case "<<argv[1]<<": turn "<<j<<" finished."<<endl;
    } //end of main loop

    //fout2<<"EOF"<<endl;

    return 0;
}

