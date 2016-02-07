#ifndef PHYSICS_H
#define PHYSICS_H
#include <vector>
#include <cmath>
#include <fstream>
#include "elepass.h"
#include "tracking.h"
#include "global.h"
extern GlobalVariables GP;
//================================================================
//
//    PHYSICS  calculation
//
//================================================================

//---5D orbit and Twiss 
void Cal_Orbit_Num(Line & linename, double deltap);
void Cal_Orbit_DA(Line & linename, double deltap);
void Cal_OneTurnMap(Line & linename, double deltap);
void Cal_ElementMap(Line & linename, double deltap);
void Cal_SectionMap(Line & linename, int i1, int i2, double deltap, double t66[36] );
//------original 4-D approach
void Cal_A0(Line & linename, double deltap);
//-----adapted from 6-D general approach
void Cal_A(Line & linename, double deltap);
void Trace_A(Line & linename, double deltap);
void Cal_Twiss(Line & linename, double deltap);
//------trace the Twiss parameters through the ring to the end
void Trace_Twiss(Line & linename, double deltap, double x[4], double Beta1, double Beta2, double Alfa1, double Alfa2, double c11, double c12, double c21, double c22);
//-----trace the orbit through the ring to the end 
void Trace_Orbit(Line & linename, double x[6], int i1);
//---6D orbit and Twiss
void Cal_Orbit_Num_6D(Line & linename);
void Cal_OneTurnMap_6D(Line & linename);
// there is not such clsoed orbits for different constant deltaps like in 5-d simulation
void Cal_ElementMap_6D(Line & linename);
void Cal_A_6D(Line & linename);
void Trace_A_6D(Line & linename);
void Cal_Twiss_6D(Line & linename);
void Cal_Chrom( Line & linename);
void Cal_Dispersion(Line & linename);
void Cal_Tune_Num(Line & linename, double deltap0);


//-----------------------------------------
//   fitting tunes and linear chroms
//-----------------------------------------
void Fit_Tune(Line & linename, double q1, double q2, const char * qf_name, const char * qd_name);
void Fit_Tune_RHICelens(Line & linename, double q1, double q2);
void Fit_Chrom(Line & linename, double chrom1x_want, double chrom1y_want, const char * sf_name, const char * sd_name );
void Fit_Chrom_RHIC8fam(Line & linename, double chrom1x_want, double chrom1y_want );


//----------------------------------------
//    Chromatic Calculation 
//----------------------------------------
void  chrom_fit(double qx[],double qy[],double & chromx1,double & chromy1,double & chromx2,double & chromy2,double & chromx3,double & chromy3 );
void Cal_Chrom_Num( Line & linename);
void Correct_Chrom_Manual( Line & linename);
void Cal_Tune_vs_Deltap(Line & linename, const char *filename);
void Plot_Tune_vs_Deltap(Line & linename, const char* filename);
void Cal_Beta_Star_vs_Deltap(Line & linename, const char *filename);
void Plot_Beta_Star_vs_Deltap(Line & linename, const char* filename);
void Cal_Beta_vs_Deltap(Line & linename, const char *filename);
void Cal_Dispersion_vs_Deltap(Line & linename, const char *filename);
void Cal_Half_Integer_RDT(Line & linename,  const char* filename);
//  horizontal: h20001, vertical: h00021
void Cal_Q2_Source(Line & linename, const char* filename );
 

//--------------------------------------------------
//  RDTs and detunings
//--------------------------------------------------
void Cal_Coupling_Coefficient( Line & linename );
// coupling coefficient from whole ring
void Cal_Coupling_Coefficient_Single( Line & linename, const char * filename);
// single element's coupling contribution
void Cal_Sext_RDTs( Line & linename );
// sextupole linear RDT at the staring point
void Cal_Detuning_Sext( Line & linename );
void Cal_Detuning_Oct( Line & linename );


//-----------------------------------------
//   path length and gamma-t 
//-----------------------------------------
void Cal_Gammat(Line & linename);
void Cal_Orbit_Length(Line & linename, double deltap);


//-------------------------------
//  longitudinal calculations
//--------------------------------
void Cal_Qs(Line & linename);
void Cal_Bucket_Area(Line & linename);
double RF_F_function(double phi_s, double phi_right, double phi_left);
void Cal_Bunch_Area(Line & linename, double full_length);
// full length is +/-3sigma_l, that is, 6 sigma_l, in units of ns
void Cal_Bunch_Height(Line & linename, double bunch_area);
void Print_Longitudinal_Summary( Line & linename);


//------------------------------------------------
//        Print and Plot
//------------------------------------------------
void Print_Twiss(Line & linename, const char* filename);
void Print_Twiss_6D(Line & linename, const char* filename);
void Print_A_Matrix(Line & linename, const char* filename);
void Print_Optics_Summary(Line & linename);
void Plot_Twiss(Line & linename, double s1, double s2);
void Plot_Orbit(Line & linename, double s1, double s2);


//---------------------------------------
//   Artifical phase rotator matrix
//--------------------------------------
void Add_Phaser(Line & linename, int loc, const char * name,  double mux, double muy);


//----------------------------------------
//  ORBIT correction
//----------------------------------------
void Correct_Orbit_SVD(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane);
void Local_Three_Bump(Line linename, int plane,  const char *corr1,   const char *corr2,   const char *corr3, double kick1);
double RMS_Leakage_Orbit( Line linename, int plane, int i1, int i2 );
void Correct_Orbit_SlidingBump1(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane);
void Correct_Orbit_SlidingBump2(Line & linename, int m, int n, vector<int> bpm_index, vector<int> kicker_index, int plane);
void Orbit_Status( Line linename, vector<int> bpm_index, int plane, double &orbit_mean, double & orbit_rms );


//--------------------------------------------------
//      Nonlinear tools: FOOTPRINT, FMA
//--------------------------------------------------
void Track_tbt_FMA( Line & linename, double deltap0, double sigmax0, double sigmay0 );
// I prefer RF off for this 
void Track_tbt_tune_footprint( Line & linename, double deltap0, double emitx, double emity);
// I prefer RF off for this 

#endif
