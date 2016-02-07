#ifndef GLOBAL_H
#define GLOBAL_H
#include <iostream>
//================================================================
//
//  Global  constants 
//
//===============================================================
const double PI      = 3.14159265;
const double light_speed =  299792458;
const double proton_mass =  938.272046;     // in units of MeV
const double electron_mass =  0.510998928;  // in units of MeV
const double Fdrift1 = 1.e0 / (2.e0 * (2.e0 - pow(2.e0, 1.0/3.0 )));
const double Fdrift2 = 0.5e0 - Fdrift1;
const double Fkick1  = 2.e0 * Fdrift1;
const double Fkick2  = 1.e0 - 2e0 * Fkick1;
const double QLslice = 0.02;
const double BLslice = 0.02;
const double seed=1232346;
enum ps_index { x_ = 0, px_ = 1, y_ = 2, py_ = 3,  z_ = 4, delta_ = 5 };  //  z > 0 earlier than synchronous particle

//================================================================
//
//      Global variables and flags
//
//===============================================================

//-----global beam and tracking parameters
class  GlobalVariables
{
 public:
  GlobalVariables()
    {
      energy=250000.;          // particle energy of one nuclei in MeV
      gamma =266;              // beambeam needs it
      beta = 0.999955981502;   // particle volecity dived by light speed = sqrt(1-1.0/gamma/gamma) 
      brho = 831.763013;       // SOLEN needs it  brho= e/p_0
      harm = 360;             // for accelerating RF cavities
      A= 1;                   // Atom index  197 
      Q= 1;                   // Particle charge   79 
      bbscale =1 ;             // force scaling for beambeam, compared to p-p collision, 
                               // for example: bbscale =1.0 for proton-proton,  bbscale =79*79/197=31.68 for Au-Au,.. 
      turn = 0 ;               // for tracking turn control ( used by ACKICK, ACDIPOLE, ACMULT )
      step_deltap=0.0001;      // used in numeric chrom calculations
      twiss_6d = 0;            // for 5D or 6D Twiss calculation  , 0 will treat RF as drift in Twiss calculation      
    }
  int    turn, A, Q, harm, twiss_6d ;             
  double energy, gamma, beta, brho, bbscale, step_deltap;
};
//GlobalVariables GP; 

/*
void Print_GlobalVariables()
{
  cout<<"-------------------------------"<<endl;
  cout<<"GP.energy  :  "<<GP.energy<<endl;
  cout<<"GP.gamma   :  "<<GP.gamma<<endl;
  cout<<"GP.beta    :  "<<GP.beta<<endl;
  cout<<"GP.brho    :  "<<GP.brho<<endl;
  cout<<"GP.harm    :  "<<GP.harm<<endl;
  cout<<"GP.A       :  "<<GP.A<<endl;
  cout<<"GP.Q       :  "<<GP.Q<<endl;
  cout<<"GP.bbscale :  "<<GP.bbscale<<endl;
  cout<<"GP.turn    :  "<<GP.turn<<endl; 
  cout<<"GP.step_deltap :  "<<GP.step_deltap<<endl;
  cout<<"GP.twiss_6d:  "<<GP.twiss_6d<<endl;
  cout<<"-------------------------------"<<endl;
}
*/


#endif
