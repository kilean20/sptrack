#ifndef ELEPASS_H
#define ELEPASS_H
#include <cstring>
#include <vector>
#include <cmath>
#include "tpsa.h"
#include "smath.h"
#include "global.h"
extern GlobalVariables GP;
//declare of line, element, and pass method
//===========================================
//
//       Element classes
//
//============================================

class Element
{
public:
  Element (string name);
  virtual void SetP (const char *name, double value) = 0;
  virtual double GetP (const char *name) = 0;
  virtual void Pass (double x[6]) = 0;
  virtual void DAPass (tps x[6]) = 0;
  string NAME, TYPE;
  double S, L, DX, DY, DT;
  double X[6], T[36], M[36], A[36];
  double Beta1, Alfa1, Beta2, Alfa2, Beta3, Alfa3, Mu1, Mu2, Mu3;
  double r, c11, c12, c21, c22;
  double Etax, Etay, Etaxp, Etayp;
  double APx, APy;
};

//---------------DRIFT-----------------------------------
class DRIFT:public Element
{
public:
  DRIFT (string name, double l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
};

//--------------------------SBEND------------------------------------
class SBEND:public Element
{
public:
  SBEND (string name, double l, double angle, double e1, double e2);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double ANGLE, E1, E2;
  int Nint;
};

//---------------------QUAD-----------------------------------------
class QUAD:public Element
{
public:
  QUAD (string name, double l, double k1l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double K1L;
  int Nint, Norder;
};

//--------------------------------SKEWQ------------------------------
class SKEWQ:public Element
{
public:
  SKEWQ (string name, double l, double k1sl);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double K1SL;
  int Nint, Norder;
};

//----------------------------SEXT-------------------------------------------
class SEXT:public Element
{
public:
  SEXT (string name, double l, double k2l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double K2L;
  int Nint, Norder;
};

//-------------------------OCT---------------------------------------------
class OCT:public Element
{
public:
  OCT (string name, double l, double k3l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double K3L;
  int Nint, Norder;
};

//-------------------------------MULT--------------------------------------
class MULT:public Element
{
public:
  MULT (string name, double l, double knl[11], double knsl[11]);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double KNL[11], KNSL[11];
  int Nint, Norder;
};

//-------------------------------GMULT--------------------------------------
class GMULT:public Element
{
public:
  GMULT (string name, double l, double angle, double e1, double e2,
	 double knl[11], double knsl[11]);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double ANGLE, E1, E2, KNL[11], KNSL[11];
  int Nint, Norder;
};

//---------------------------KICKER-----------------------------------------
template < class T >
  void KICK_Pass (T x[6], double L, double HKICK, double VKICK);

class KICKER:public Element
{
public:
  KICKER (string name, double l, double hkick, double vkick);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double HKICK, VKICK;
};

//----------------------------HKICKER--------------------------------------
class HKICKER:public Element
{
public:
  HKICKER (string name, double l, double hkick);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double HKICK;			// kick in x'
};

//------------------------------VKICKER------------------------------------------
class VKICKER:public Element
{
public:
  VKICKER (string name, double l, double vkick);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double VKICK;			// kick in y'
};

//---------------------------HACKICK ( amplitude constant ac dipole kicking)------------------------------------
template < class T >
  void ACKICK_Pass (T x[6], double L, double HKICK, double VKICK,
		    double TTURNS, double PHI0);


class HACKICK:public Element
{
public:
  HACKICK (string name, double l, double hkick, int tturns, double phi0);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  int TTURNS;			//  T period in unit of turns
  double HKICK, PHI0;		//  HKICK is kick amplitude
};

//---------------------------VACKICK ( amplitude constant ac dipole kicking)------------------------------------
class VACKICK:public Element
{
public:
  VACKICK (string name, double l, double vkick, int tturns, double phi0);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  int TTURNS;			// T period in unit of turns 
  double VKICK, PHI0;
};

//------------------------------HACDIP ( slowly ramp up the kick amplitude)-----------------------------------------
template < class T >
  void ACDIP_Pass (T x[6], double L, double HKICKMAX, double VKICKMAX,
		   double NUD, double TURNS, double TURNE, double PHID);


class HACDIP:public Element
{
public:
  HACDIP (string name, double l, double hkickmax, double nud, double phid,
	  int turns, int turne);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  int TURNS, TURNE;		// kick strength ramping between  TURNS and TURNE turns
  double HKICKMAX, NUD, PHID;	// HKICKMAX, maximum kick in x'
};

//------------------------------VACDIP------------------------------------------
class VACDIP:public Element
{
public:
  VACDIP (string name, double l, double vkickmax, double nud, double phid,
	  int turns, int turne);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  int TURNS, TURNE;
  double VKICKMAX, NUD, PHID;	// HKICKMAX, maximum kick in y'
};

//----------------------------ACMULT ( single order, order >= quadrupole )-----------------
void ACMULT_Pass (double x[6], double L, int Norder, double KL, int TTURNS,
		  double PHI0);


class ACMULT:public Element
{
public:
  ACMULT (string name, double l, int norder, double kl, int tturns,
	  double phi0);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  int Norder, TTURNS;		// TTURNs period in numbers of turns, 
  double KL, PHI0;		//  Norder >0, KL->KNL, Norder <0, KL->KNSL 
};

//---------------------------COOLING--------------------------------
template < class T > void COOLING_Pass (T x[6], double L, double ALPHA);

class COOLING:public Element
{
public:
  COOLING (string name, double l, double alpha);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double ALPHA;
};

//---------------------------DIFFUSE--------------------------------
template < class T >
  void DIFFUSE_Pass (T x[6], double L, double DIFF_X, double DIFF_Y,
		     double DIFF_DELTA);


class DIFFUSE:public Element
{
public:
  DIFFUSE (string name, double l, double diff_x, double diff_y,
	   double diff_delta);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double DIFF_X, DIFF_Y, DIFF_DELTA;
};

//---------------------------------BPM-------------------------
class BPM:public Element
{
public:
  BPM (string name, double l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
};

//----------------------------HBPM-----------------------------------
class HBPM:public Element
{
public:
  HBPM (string name, double l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
};

//---------------------------------VBPM------------------------------
class VBPM:public Element
{
public:
  VBPM (string name, double l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
};

//-----------------------------MARKER--------------------------
class MARKER:public Element
{
public:
  MARKER (string name, double l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
};

//---------------------------RCOLL--------------------------------
class RCOLL:public Element
{
public:
  RCOLL (string name, double l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
};

//----------------------------ECOLL-------------------------------
class ECOLL:public Element
{
public:
  ECOLL (string name, double l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
};

//----------------------INSTR------------------------------------
class INSTR:public Element
{
public:
  INSTR (string name, double l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
};

//----------------------------SOLEN-----------------------------------
template < class T >
  void SOLEN_Pass (T x[6], double L, double KS, double M[36]);

class SOLEN:public Element
{
public:
  SOLEN (string name, double l, double ks);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double KS, M[36];
};

//-----------------------MATRIX-----------------------------------
template < typename T >
  void MATRIX_Pass (T x[6], double XCO_IN[6], double XCO_OUT[6],
		    double M66[36]);


class MATRIX:public Element
{
public:
  MATRIX (string name, double l, double m66[36], double xco_in[6],
	  double xco_out[6]);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double M66[36];
  double XCO_IN[6], XCO_OUT[6];
};

//--------------------ELSEP----------------------------------
class ELSEP:public Element
{
public:
  ELSEP (string name, double l);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
};

//--------------------OFFSET (coordinate system change)-------------------------------
template < typename T >
  void OFFSET_Pass (T x[6], double DX, double DPX, double DY, double DPY,
		    double TILT);


class OFFSET:public Element
{
public:
  OFFSET (string name, double l, double dx, double dpx, double dy, double dpy,
	  double tilt);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double DX, DPX, DY, DPY, TILT;
};

//----------------------------RFCAV--------------------------------
template < class T >
  void RFCAV_Pass (T x[6], double L, double VRF, double FRF, double PHASE0);


class RFCAV:public Element
{
public:
  RFCAV (string name, double l, double vrf, double frf, double phase0);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double VRF, FRF, PHASE0;
};

//----------------------BEAMBEAM------------------------------
void cerrf (double xx, double yy, double &wx, double &wy);
template < class T >
void BB4D (T x, T y, double gamma, double N, double sigmax,
	   double sigmay, double &Dpx, double &Dpy);
template < class T >
void BB6D (T x[6], double gamma, double Np, double sigma_l, int N_slice,
	   double emitx_rms, double betax_star, double alfx_star,
	   double emity_rms, double betay_star, double alfy_star);
 // emitx_rms, emity_rms:  un-normalized rms emittance,  sigma=SQRT[ emitx_rms * betax ]  


void BEAMBEAM_Pass (double x[6], int TREATMENT, double NP, double SIGMAL,
		      int NSLICE, double EMITX, double BETAX, double ALFAX,
		      double EMITY, double BETAY, double ALFAY);


class BEAMBEAM:public Element
{
public:
  BEAMBEAM (string name, int treatment, double np, double sigmal, int nslice,
	    double emitx, double emity, double betax, double alfax,
	    double betay, double alfay);
  //     emitx, emity:  un-normalized rms emittance,  sigma=SQRT[ emitx * betax ]  
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  int NSLICE, TREATMENT;
  double NP, EMITX, EMITY, SIGMAL, BETAX, ALFAX, BETAY, ALFAY;
};

//---------------------LRBB ( long-range BEAMBEAM )------------
  void LRBB_Pass (double x[6], double gamma, double NP, double SEPX, double SEPY,
		  double SIGMAX, double SIGMAY, double KICKX0, double KICKY0);
//  here we assume the particle's coordinates from the strong beam center is ( SEPX+x[1], SEPY+x[3] )


class LRBB:public Element
{
public:
  LRBB (string name, double np, double sepx, double sepy, double sigmax,
	double sigmay);
  //     sepx, sepy are the offset of the weak beam center w.r.t to the center strong beam  
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double NP, SEPX, SEPY, SIGMAX, SIGMAY;
  double KICKX0, KICKY0;
};

//----------------------E-LENS------------------------------
  void ELENS_Pass (double x[6], double Ne, double Le, double beta_e, int N_slice,
		   double sigmax, double sigmay);


void elens_pass_round_Gaussian_topoff (double x[6], double gamma,
				       double Ne, double Le, double beta_e,
				       int N_slice, double sigmax,
				       double sigmay);
// I assume sigmax=sigmay, top off [-sigmax*a, sigmax*a ], Ne from perfect Gaussian


void elens_pass_round_Gaussian_truncated (double x[6], double gamma,
					  double Ne, double Le, double beta_e,
					  int N_slice, double sigmax,
					  double sigmay);
// I assume sigmax=sigmay, Gaussian tail cut off from Nc*sigmax


void elens_pass_round_uniform (double x[6], double gamma,
			       double Ne, double Le, double beta_e,
			       int N_slice, double sigmax, double sigmay);
// I assume round uniform distribution in r<=a


class ELENS:public Element
{
public:
  ELENS (string name, double l, double ne, int nslice, double betae,
	 double sigmax, double sigmay);
  void SetP (const char *name, double value);
  double GetP (const char *name);
  void Pass (double x[6]);
  void DAPass (tps x[6]);
private:
  double NE, BETAE, SIGMAX, SIGMAY;
  int NSLICE;
};


//=========================================
//
//            Line  class
//
//=========================================
class Line
{
public:
  Line ();
  void Update ();
  void Append (Element * x);
  void Delete (int i);
  void Insert (int i, Element * temp);
  void Empty ();

    vector < Element * >Cell;
  double Length;		//---length of reference orbit
  double frev0;			//---decided by Length
  double frf;			//---rf frequency
  double Vrf_tot;		//---total rf voltage
  double Orbit_Length;		//---length of closed orbit length.
  long Ncell;			//---number of elements 
  double Tune1, Tune2, Tune3;	//----tunes
  double Chromx1, Chromy1, Chromx2, Chromy2, Chromx3, Chromy3;	//   chromaticities
  double Alfa0, Alfa1, Alfa2, Gammat, Slip;	//  transition parameters
  double Qs;			//  longitudinal tune
  double Bucket_length;		//  in unit of ns 
  double Bucket_height;		//  (dp/p0)_max
  double Bucket_area;		//  in ( phi_rf, dE/(hw_rev) ) phase space per nucleon
  double Bunch_length;		//  in unit of m, 6*sigma_l, in unit of ns
  double Bunch_area;		//  in ( phi_rf, dE/(hw_rev) ) phase space per nucleon
  double Bunch_height;		//  (dp/p0)_max for the given bunch area
};

void Line_Rewind (Line & linename, int k);
void Line_Invert (Line & linename);
void Line_Repeat (Line & linename1, Line & linename2, int n);
void Line_Connect (Line & linename1, Line & linename2, Line & linename3);
int Get_Index (Line & linename, const char *name, int k);
void Split_Quad (Line & linename, int i, int m);
void Split_Sext (Line & linename, int i, int m);
void Split_Quad_Sext (Line & linename);
void Split_Sbend (Line & linename, int i, int m);
void Split_Drift (Line & linename, int i, int m);
void Split_Mult (Line & linename, int i, int m);


//----function: concate adjacient DRIFT, doesn't change Twiss
void Concat_Drift (Line & linename);


//----get rid of zero strength element, doesn't change Twiss
void Clean_Up (Line & linename);


//----prepare for fast tarcking: will change Twiss
void Make_Thin (Line & linename);
double Get_KL (Line & linename, const char *name, const char *kl);
void Set_KL (Line & linename, const char *name, const char *kl,
	     double strength);
void Set_dKL (Line & linename, const char *name, const char *kl,
	      double dstrength);


//=========================================
//
//            Magnete Alignment
//
//=========================================
template < typename T > void GtoL (T x[6], double DX, double DY, double DT);
template < typename T > void LtoG (T x[6], double DX, double DY, double DT);


//=================================================================
//
//    4th order Symplectic Integrator  (S.I.)
//
//=================================================================
template < class T > void DRIFT_Pass (T x[6], double L);
template < class T > void bend_kick_pass (T x[6], double L, double href);
template < class T > void quad_kick_pass (T x[6], double k1l, double k1sl);
template < class T > void sext_kick_pass (T x[6], double k2l, double k2sl);
template < class T > void oct_kick_pass (T x[6], double k3l, double k3sl);
template < class T > void mult_kick_pass (T x[6], int Norder, double KNL[11],
					  double KNSL[11]);
template < class T > void bend_mult_kick_pass (T x[6], double L, double href,
					       int Norder, double KNL[11],
					       double KNSL[11]);
template < class T > void SBEND_Pass (T x[6], double L, int Nint,
				      double Angle, double E1, double E2);
template < class T > void QUAD_Pass (T x[6], double L, int Nint, double k1l,
				     double k1sl);
template < class T > void SEXT_Pass (T x[6], double L, int Nint, double k2l,
				     double k2sl);
template < class T > void OCT_Pass (T x[6], double L, int Nint, double k3l,
				    double k3sl);
template < class T > void MULT_Pass (T x[6], double L, int Nint, int Norder,
				     double KNL[11], double KNSL[11]);
template < class T > void GMULT_Pass (T x[6], double L, int Nint,
				      double Angle, int Norder,
				      double KNL[11], double KNSL[11],
				      double E1, double E2);

#endif
