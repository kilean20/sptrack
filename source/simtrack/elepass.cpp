#include "elepass.h"
using namespace std;
//===========================================
//
//       Element classes
//
//============================================

Element::Element (string name)
{
  int i;
  NAME = name;
  DX = 0.;			// offset of magnet center w.r.t x-y frame, DX>0 means magnet moved to positive x
  DY = 0.;			// offset of magnet center w.r.t x-y frame, DY>0 means magnet moved to positive y
  DT = 0.;			// tilt angle of magnet w.r.t x-y frame, roll from e_x-->e_y, DT > 0 ;  Normal Q rolled -M_PI/4 to get a skewQ, k1s=k1
  for (i = 0; i < 6; i++)
    X[i] = 0.;
  for (i = 0; i < 36; i++)
    T[i] = 0.;
  for (i = 0; i < 36; i++)
    M[i] = 0.;
  for (i = 0; i < 36; i++)
    A[i] = 0.;
  Beta1 = 0.;
  Beta2 = 0.;
  Beta3 = 0.;
  Alfa1 = 0.;
  Alfa2 = 0.;
  Alfa3 = 0.;
  r = 0;
  c11 = 0.;
  c12 = 0.;
  c21 = 0.;
  c22 = 0.;
  Etax = 0.;
  Etay = 0.;
  Etaxp = 0.;
  Etayp = 0.;
  Mu1 = 0.;
  Mu2 = 0.;
  Mu3 = 0.;
  APx = 1.;
  APy = 1.;
}

//---------------DRIFT-----------------------------------
DRIFT::DRIFT (string name, double l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("DRIFT");
      L = l;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}

void
DRIFT::SetP (const char *name, double value)
{
  cout << "No parameter to be set for DRIFT." << endl;
  exit (0);
}

double
DRIFT::GetP (const char *name)
{
  cout << "No parameter to be returned for DRIFT." << endl;
  exit (0);
}

void
DRIFT::Pass (double x[6])
{
  DRIFT_Pass (x, L);
}

void
DRIFT::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//--------------------------SBEND------------------------------------
SBEND::SBEND (string name, double l, double angle, double e1, double e2):
Element (name)
{
  if (l > 0)
    {
      TYPE = string ("SBEND");
      L = l;
      ANGLE = angle;
      E1 = e1;
      E2 = e2;
      Nint = int (L / BLslice) + 1;
    }
  else
    {
      cout << "Error: SBEND length can not be zero." << endl;
      exit (1);
    }
}
void
SBEND::SetP (const char *name, double value)
{
  if (strcmp (name, "ANGLE") == 0)
    {
      ANGLE = value;
    }
  else if (strcmp (name, "E1") == 0)
    {
      E1 = value;
    }
  else if (strcmp (name, "E2") == 0)
    {
      E2 = value;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      Nint = int (value);
    }
  else
    {
      cout << "SBEND does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
SBEND::GetP (const char *name)
{
  if (strcmp (name, "ANGLE") == 0)
    {
      return ANGLE;
    }
  else if (strcmp (name, "E1") == 0)
    {
      return E1;
    }
  else if (strcmp (name, "E2") == 0)
    {
      return E2;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      return Nint;
    }
  else
    {
      cout << "SBEND does not have parameter of " << name << endl;
      exit (0);
    }
}
void
SBEND::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  SBEND_Pass (x, L, Nint, ANGLE, E1, E2);
  LtoG (x, DX, DY, DT);
}

void
SBEND::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  SBEND_Pass (x, L, Nint, ANGLE, E1, E2);
  LtoG (x, DX, DY, DT);
}

//---------------------QUAD-----------------------------------------
QUAD::QUAD (string name, double l, double k1l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("QUAD");
      L = l;
      K1L = k1l;
      Nint = int (L / QLslice) + 1;
      Norder = 1;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
QUAD::SetP (const char *name, double value)
{
  if (strcmp (name, "K1L") == 0)
    {
      K1L = value;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      Nint = int (value);
    }
  else
    {
      cout << "QUAD does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
QUAD::GetP (const char *name)
{
  if (strcmp (name, "K1L") == 0)
    {
      return K1L;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      return Nint;
    }
  else
    {
      cout << "QUAD does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
QUAD::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  QUAD_Pass (x, L, Nint, K1L, 0.);
  LtoG (x, DX, DY, DT);
}

void
QUAD::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  QUAD_Pass (x, L, Nint, K1L, 0.);
  LtoG (x, DX, DY, DT);
}

//--------------------------------SKEWQ------------------------------
SKEWQ::SKEWQ (string name, double l, double k1sl):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("SKEWQ");
      L = l;
      K1SL = k1sl;
      Nint = int (L / QLslice) + 1;
      Norder = 1;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
SKEWQ::SetP (const char *name, double value)
{
  if (strcmp (name, "K1SL") == 0)
    {
      K1SL = value;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      Nint = int (value);
    }
  else
    {
      cout << "SKEWQ does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
SKEWQ::GetP (const char *name)
{
  if (strcmp (name, "K1SL") == 0)
    {
      return K1SL;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      return Nint;
    }
  else
    {
      cout << "SKEWQ does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
SKEWQ::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  QUAD_Pass (x, L, Nint, 0., K1SL);
  LtoG (x, DX, DY, DT);
}

void
SKEWQ::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  QUAD_Pass (x, L, Nint, 0., K1SL);
  LtoG (x, DX, DY, DT);
}

//----------------------------SEXT-------------------------------------------
SEXT::SEXT (string name, double l, double k2l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("SEXT");
      L = l;
      K2L = k2l;
      Nint = int (L / QLslice) + 1;
      Norder = 2;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
SEXT::SetP (const char *name, double value)
{
  if (strcmp (name, "K2L") == 0)
    {
      K2L = value;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      Nint = int (value);
    }
  else
    {
      cout << "SEXT does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
SEXT::GetP (const char *name)
{
  if (strcmp (name, "K2L") == 0)
    {
      return K2L;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      return Nint;
    }
  else
    {
      cout << "SEXT does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
SEXT::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  SEXT_Pass (x, L, Nint, K2L, 0.);
  LtoG (x, DX, DY, DT);
}

void
SEXT::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  SEXT_Pass (x, L, Nint, K2L, 0.);
  LtoG (x, DX, DY, DT);
}

//-------------------------OCT---------------------------------------------
OCT::OCT (string name, double l, double k3l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("OCT");
      L = l;
      K3L = k3l;
      Nint = int (L / QLslice) + 1;
      Norder = 3;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
OCT::SetP (const char *name, double value)
{
  if (strcmp (name, "K3L") == 0)
    {
      K3L = value;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      Nint = int (value);
    }
  else
    {
      cout << "OCT does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
OCT::GetP (const char *name)
{
  if (strcmp (name, "K3L") == 0)
    {
      return K3L;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      return Nint;
    }
  else
    {
      cout << "OCT does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
OCT::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  OCT_Pass (x, L, Nint, K3L, 0.);
  LtoG (x, DX, DY, DT);
}

void
OCT::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  OCT_Pass (x, L, Nint, K3L, 0.);
  LtoG (x, DX, DY, DT);
}

//-------------------------------MULT--------------------------------------
MULT::MULT (string name, double l, double knl[11], double knsl[11]):
Element (name)
{
  int i;
  if (l >= 0)
    {
      TYPE = string ("MULT");
      L = l;
      for (i = 0; i < 11; i++)
	{
	  KNL[i] = knl[i];
	  KNSL[i] = knsl[i];
	}
      Nint = int (L / QLslice) + 1;
      Norder = 1;
      for (i = 0; i < 10; i++)
	{
	  if (KNL[10 - i] != 0. || KNSL[10 - i] != 0.)
	    {
	      Norder = 10 - i;
	      break;
	    }
	}
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
MULT::SetP (const char *name, double value)
{
  if (strcmp (name, "K0L") == 0)
    {
      KNL[0] = value;
    }
  else if (strcmp (name, "K1L") == 0)
    {
      KNL[1] = value;
    }
  else if (strcmp (name, "K2L") == 0)
    {
      KNL[2] = value;
    }
  else if (strcmp (name, "K3L") == 0)
    {
      KNL[3] = value;
    }
  else if (strcmp (name, "K4L") == 0)
    {
      KNL[4] = value;
    }
  else if (strcmp (name, "K5L") == 0)
    {
      KNL[5] = value;
    }
  else if (strcmp (name, "K6L") == 0)
    {
      KNL[6] = value;
    }
  else if (strcmp (name, "K7L") == 0)
    {
      KNL[7] = value;
    }
  else if (strcmp (name, "K8L") == 0)
    {
      KNL[8] = value;
    }
  else if (strcmp (name, "K9L") == 0)
    {
      KNL[9] = value;
    }
  else if (strcmp (name, "K10L") == 0)
    {
      KNL[10] = value;
    }
  else if (strcmp (name, "K0SL") == 0)
    {
      KNSL[0] = value;
    }
  else if (strcmp (name, "K1SL") == 0)
    {
      KNSL[1] = value;
    }
  else if (strcmp (name, "K2SL") == 0)
    {
      KNSL[2] = value;
    }
  else if (strcmp (name, "K3SL") == 0)
    {
      KNSL[3] = value;
    }
  else if (strcmp (name, "K4SL") == 0)
    {
      KNSL[4] = value;
    }
  else if (strcmp (name, "K5SL") == 0)
    {
      KNSL[5] = value;
    }
  else if (strcmp (name, "K6SL") == 0)
    {
      KNSL[6] = value;
    }
  else if (strcmp (name, "K7SL") == 0)
    {
      KNSL[7] = value;
    }
  else if (strcmp (name, "K8SL") == 0)
    {
      KNSL[8] = value;
    }
  else if (strcmp (name, "K9SL") == 0)
    {
      KNSL[9] = value;
    }
  else if (strcmp (name, "K10SL") == 0)
    {
      KNSL[10] = value;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      Nint = int (value);
    }
  else
    {
      cout << "MULT does not have parameter of  " << name << endl;
      exit (0);
    }
}
double
MULT::GetP (const char *name)
{
  if (strcmp (name, "K0L") == 0)
    {
      return KNL[0];
    }
  else if (strcmp (name, "K1L") == 0)
    {
      return KNL[1];
    }
  else if (strcmp (name, "K2L") == 0)
    {
      return KNL[2];
    }
  else if (strcmp (name, "K3L") == 0)
    {
      return KNL[3];
    }
  else if (strcmp (name, "K4L") == 0)
    {
      return KNL[4];
    }
  else if (strcmp (name, "K5L") == 0)
    {
      return KNL[5];
    }
  else if (strcmp (name, "K6L") == 0)
    {
      return KNL[6];
    }
  else if (strcmp (name, "K7L") == 0)
    {
      return KNL[7];
    }
  else if (strcmp (name, "K8L") == 0)
    {
      return KNL[8];
    }
  else if (strcmp (name, "K9L") == 0)
    {
      return KNL[9];
    }
  else if (strcmp (name, "K10L") == 0)
    {
      return KNL[10];
    }
  else if (strcmp (name, "K0SL") == 0)
    {
      return KNSL[0];
    }
  else if (strcmp (name, "K1SL") == 0)
    {
      return KNSL[1];
    }
  else if (strcmp (name, "K2SL") == 0)
    {
      return KNSL[2];
    }
  else if (strcmp (name, "K3SL") == 0)
    {
      return KNSL[3];
    }
  else if (strcmp (name, "K4SL") == 0)
    {
      return KNSL[4];
    }
  else if (strcmp (name, "K5SL") == 0)
    {
      return KNSL[5];
    }
  else if (strcmp (name, "K6SL") == 0)
    {
      return KNSL[6];
    }
  else if (strcmp (name, "K7SL") == 0)
    {
      return KNSL[7];
    }
  else if (strcmp (name, "K8SL") == 0)
    {
      return KNSL[8];
    }
  else if (strcmp (name, "K9SL") == 0)
    {
      return KNSL[9];
    }
  else if (strcmp (name, "K10SL") == 0)
    {
      return KNSL[10];
    }
  else if (strcmp (name, "Norder") == 0)
    {
      return Norder;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      return Nint;
    }
  else
    {
      cout << "MULT does not have parameter of  " << name << endl;
      exit (0);
    }
}
void
MULT::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  MULT_Pass (x, L, Nint, Norder, KNL, KNSL);
  LtoG (x, DX, DY, DT);
}

void
MULT::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  MULT_Pass (x, L, Nint, Norder, KNL, KNSL);
  LtoG (x, DX, DY, DT);
}

//-------------------------------GMULT--------------------------------------
GMULT::GMULT (string name, double l, double angle, double e1, double e2,
	      double knl[11], double knsl[11]):
Element (name)
{
  int i;
  if (l >= 0)
    {
      TYPE = string ("GMULT");
      L = l;
      ANGLE = angle;
      E1 = e1;
      E2 = e2;
      Nint = int (L / QLslice) + 1;
      for (i = 0; i < 11; i++)
	{
	  KNL[i] = knl[i];
	  KNSL[i] = knsl[i];
	}
      Norder = 1;
      for (i = 0; i < 10; i++)
	{
	  if (KNL[10 - i] != 0. || KNSL[10 - i] != 0.)
	    {
	      Norder = 10 - i;
	      break;
	    }
	}
    }
  else
    {
      cout << "Error: GMULT length can not be negative." << endl;
      exit (1);
    }
}
void
GMULT::SetP (const char *name, double value)
{
  if (strcmp (name, "ANGLE") == 0)
    {
      ANGLE = value;
    }
  else if (strcmp (name, "E1") == 0)
    {
      E1 = value;
    }
  else if (strcmp (name, "E2") == 0)
    {
      E2 = value;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      Nint = int (value);
    }
  else if (strcmp (name, "K0L") == 0)
    {
      KNL[0] = value;
    }
  else if (strcmp (name, "K1L") == 0)
    {
      KNL[1] = value;
    }
  else if (strcmp (name, "K2L") == 0)
    {
      KNL[2] = value;
    }
  else if (strcmp (name, "K3L") == 0)
    {
      KNL[3] = value;
    }
  else if (strcmp (name, "K4L") == 0)
    {
      KNL[4] = value;
    }
  else if (strcmp (name, "K5L") == 0)
    {
      KNL[5] = value;
    }
  else if (strcmp (name, "K6L") == 0)
    {
      KNL[6] = value;
    }
  else if (strcmp (name, "K7L") == 0)
    {
      KNL[7] = value;
    }
  else if (strcmp (name, "K8L") == 0)
    {
      KNL[8] = value;
    }
  else if (strcmp (name, "K9L") == 0)
    {
      KNL[9] = value;
    }
  else if (strcmp (name, "K10L") == 0)
    {
      KNL[10] = value;
    }
  else if (strcmp (name, "K0SL") == 0)
    {
      KNSL[0] = value;
    }
  else if (strcmp (name, "K1SL") == 0)
    {
      KNSL[1] = value;
    }
  else if (strcmp (name, "K2SL") == 0)
    {
      KNSL[2] = value;
    }
  else if (strcmp (name, "K3SL") == 0)
    {
      KNSL[3] = value;
    }
  else if (strcmp (name, "K4SL") == 0)
    {
      KNSL[4] = value;
    }
  else if (strcmp (name, "K5SL") == 0)
    {
      KNSL[5] = value;
    }
  else if (strcmp (name, "K6SL") == 0)
    {
      KNSL[6] = value;
    }
  else if (strcmp (name, "K7SL") == 0)
    {
      KNSL[7] = value;
    }
  else if (strcmp (name, "K8SL") == 0)
    {
      KNSL[8] = value;
    }
  else if (strcmp (name, "K9SL") == 0)
    {
      KNSL[9] = value;
    }
  else if (strcmp (name, "K10SL") == 0)
    {
      KNSL[10] = value;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      Nint = int (value);
    }
  else
    {
      cout << "MULT does not have parameter of  " << name << endl;
      exit (0);
    }
}
double
GMULT::GetP (const char *name)
{
  if (strcmp (name, "ANGLE") == 0)
    {
      return ANGLE;
    }
  else if (strcmp (name, "E1") == 0)
    {
      return E1;
    }
  else if (strcmp (name, "E2") == 0)
    {
      return E2;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      return Nint;
    }
  else if (strcmp (name, "K0L") == 0)
    {
      return KNL[0];
    }
  else if (strcmp (name, "K1L") == 0)
    {
      return KNL[1];
    }
  else if (strcmp (name, "K2L") == 0)
    {
      return KNL[2];
    }
  else if (strcmp (name, "K3L") == 0)
    {
      return KNL[3];
    }
  else if (strcmp (name, "K4L") == 0)
    {
      return KNL[4];
    }
  else if (strcmp (name, "K5L") == 0)
    {
      return KNL[5];
    }
  else if (strcmp (name, "K6L") == 0)
    {
      return KNL[6];
    }
  else if (strcmp (name, "K7L") == 0)
    {
      return KNL[7];
    }
  else if (strcmp (name, "K8L") == 0)
    {
      return KNL[8];
    }
  else if (strcmp (name, "K9L") == 0)
    {
      return KNL[9];
    }
  else if (strcmp (name, "K10L") == 0)
    {
      return KNL[10];
    }
  else if (strcmp (name, "K0SL") == 0)
    {
      return KNSL[0];
    }
  else if (strcmp (name, "K1SL") == 0)
    {
      return KNSL[1];
    }
  else if (strcmp (name, "K2SL") == 0)
    {
      return KNSL[2];
    }
  else if (strcmp (name, "K3SL") == 0)
    {
      return KNSL[3];
    }
  else if (strcmp (name, "K4SL") == 0)
    {
      return KNSL[4];
    }
  else if (strcmp (name, "K5SL") == 0)
    {
      return KNSL[5];
    }
  else if (strcmp (name, "K6SL") == 0)
    {
      return KNSL[6];
    }
  else if (strcmp (name, "K7SL") == 0)
    {
      return KNSL[7];
    }
  else if (strcmp (name, "K8SL") == 0)
    {
      return KNSL[8];
    }
  else if (strcmp (name, "K9SL") == 0)
    {
      return KNSL[9];
    }
  else if (strcmp (name, "K10SL") == 0)
    {
      return KNSL[10];
    }
  else if (strcmp (name, "Norder") == 0)
    {
      return Norder;
    }
  else if (strcmp (name, "Nint") == 0)
    {
      return Nint;
    }
  else
    {
      cout << "MULT does not have parameter of  " << name << endl;
      exit (0);
    }
}
void
GMULT::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  GMULT_Pass (x, L, Nint, ANGLE, Norder, KNL, KNSL, E1, E2);
  LtoG (x, DX, DY, DT);
}

void
GMULT::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  GMULT_Pass (x, L, Nint, ANGLE, Norder, KNL, KNSL, E1, E2);
  LtoG (x, DX, DY, DT);
}

//---------------------------KICKER-----------------------------------------
template < class T > void
KICK_Pass (T x[6], double L, double HKICK, double VKICK)
{
  DRIFT_Pass (x, L / 2.0);
  x[1] = x[1] + HKICK;
  x[3] = x[3] + VKICK;
  DRIFT_Pass (x, L / 2.0);
}

KICKER::KICKER (string name, double l, double hkick, double vkick):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("KICKER");
      L = l;
      HKICK = hkick;
      VKICK = vkick;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
KICKER::SetP (const char *name, double value)
{
  if (strcmp (name, "HKICK") == 0)
    {
      HKICK = value;
    }
  else if (strcmp (name, "VKICK") == 0)
    {
      VKICK = value;
    }
  else
    {
      cout << "KICKRE does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
KICKER::GetP (const char *name)
{
  if (strcmp (name, "HKICK") == 0)
    {
      return HKICK;
    }
  else if (strcmp (name, "VKICK") == 0)
    {
      return VKICK;
    }
  else
    {
      cout << "KICKRE does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
KICKER::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  KICK_Pass (x, L, HKICK, VKICK);
  LtoG (x, DX, DY, DT);
}

void
KICKER::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  KICK_Pass (x, L, HKICK, VKICK);
  LtoG (x, DX, DY, DT);
}

//----------------------------HKICKER--------------------------------------
HKICKER::HKICKER (string name, double l, double hkick):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("HKICKER");
      L = l;
      HKICK = hkick;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
HKICKER::SetP (const char *name, double value)
{
  if (strcmp (name, "HKICK") == 0)
    {
      HKICK = value;
    }
  else
    {
      cout << "HKICKER does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
HKICKER::GetP (const char *name)
{
  if (strcmp (name, "HKICK") == 0)
    {
      return HKICK;
    }
  else
    {
      cout << "HKICKER does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
HKICKER::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  KICK_Pass (x, L, HKICK, 0.);
  LtoG (x, DX, DY, DT);
}

void
HKICKER::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  KICK_Pass (x, L, HKICK, 0.);
  LtoG (x, DX, DY, DT);
}

//------------------------------VKICKER------------------------------------------
VKICKER::VKICKER (string name, double l, double vkick):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("VKICKER");
      L = l;
      VKICK = vkick;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
VKICKER::SetP (const char *name, double value)
{
  if (strcmp (name, "VKICK") == 0)
    {
      VKICK = value;
    }
  else
    {
      cout << "VKICKER does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
VKICKER::GetP (const char *name)
{
  if (strcmp (name, "VKICK") == 0)
    {
      return VKICK;
    }
  else
    {
      cout << "KICKER does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
VKICKER::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  KICK_Pass (x, L, 0., VKICK);
  LtoG (x, DX, DY, DT);
}

void
VKICKER::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  KICK_Pass (x, L, 0., VKICK);
  LtoG (x, DX, DY, DT);
}

//---------------------------HACKICK ( amplitude constant ac dipole kicking)------------------------------------
template < class T > void
ACKICK_Pass (T x[6], double L, double HKICK, double VKICK, double TTURNS,
	     double PHI0)
{
  DRIFT_Pass (x, L / 4.0);
  x[1] = x[1] + HKICK * sin (2.0 * M_PI * GP.turn / TTURNS + PHI0);
  x[3] = x[3] + VKICK * sin (2.0 * M_PI * GP.turn / TTURNS + PHI0);
  DRIFT_Pass (x, L / 2.0);
}

HACKICK::HACKICK (string name, double l, double hkick, int tturns,
		  double phi0):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("HACKICK");
      L = l;
      HKICK = hkick;
      TTURNS = tturns;
      PHI0 = phi0;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
HACKICK::SetP (const char *name, double value)
{
  if (strcmp (name, "HKICK") == 0)
    {
      HKICK = value;
    }
  else if (strcmp (name, "TTURNS") == 0)
    {
      TTURNS = int (value);
    }
  else if (strcmp (name, "PHI0") == 0)
    {
      PHI0 = value;
    }
  else
    {
      cout << "HACKICK does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
HACKICK::GetP (const char *name)
{
  if (strcmp (name, "HKICK") == 0)
    {
      return HKICK;
    }
  else if (strcmp (name, "TTURNS") == 0)
    {
      return TTURNS;
    }
  else if (strcmp (name, "PHI0") == 0)
    {
      return PHI0;
    }
  else
    {
      cout << "HACKICK  does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
HACKICK::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  ACKICK_Pass (X, L, HKICK, 0., TTURNS, PHI0);
  LtoG (x, DX, DY, DT);
}

void
HACKICK::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//---------------------------VACKICK ( amplitude constant ac dipole kicking)------------------------------------
VACKICK::VACKICK (string name, double l, double vkick, int tturns,
		  double phi0):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("VACKICK");
      L = l;
      VKICK = vkick;
      TTURNS = int (tturns);
      PHI0 = phi0;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
VACKICK::SetP (const char *name, double value)
{
  if (strcmp (name, "VKICK") == 0)
    {
      VKICK = value;
    }
  else if (strcmp (name, "TTURNS") == 0)
    {
      TTURNS = int (value);
    }
  else if (strcmp (name, "PHI0") == 0)
    {
      PHI0 = value;
    }
  else
    {
      cout << "VACKICK does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
VACKICK::GetP (const char *name)
{
  if (strcmp (name, "VKICK") == 0)
    {
      return VKICK;
    }
  else if (strcmp (name, "TTURNS") == 0)
    {
      return TTURNS;
    }
  else if (strcmp (name, "PHI0") == 0)
    {
      return PHI0;
    }
  else
    {
      cout << "VACKICK  does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
VACKICK::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  ACKICK_Pass (X, L, 0., VKICK, TTURNS, PHI0);
  LtoG (x, DX, DY, DT);
}

void
VACKICK::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//------------------------------HACDIP ( slowly ramp up the kick amplitude)-----------------------------------------
template < class T > void
ACDIP_Pass (T x[6], double L, double HKICKMAX, double VKICKMAX, double NUD,
	    double TURNS, double TURNE, double PHID)
{
  if (GP.turn < TURNS)
    {
      DRIFT_Pass (x, L);
    }
  else if (GP.turn >= TURNS and GP.turn <= TURNE)
    {
      DRIFT_Pass (x, L / 2.);
      x[1] +=
	((GP.turn - TURNS) * 1.0 * HKICKMAX / (TURNE -
					       TURNS)) * sin (2.0 *
							      3.14159265 *
							      NUD * (GP.turn -
								     TURNS) +
							      PHID);
      x[3] +=
	((GP.turn - TURNS) * 1.0 * VKICKMAX / (TURNE -
					       TURNS)) * sin (2.0 *
							      3.14159265 *
							      NUD * (GP.turn -
								     TURNS) +
							      PHID);
      DRIFT_Pass (x, L / 2.);
    }
  else
    {
      DRIFT_Pass (x, L / 2.);
      x[1] +=
	HKICKMAX * sin (2.0 * 3.14159265 * NUD * (GP.turn - TURNS) + PHID);
      x[3] +=
	VKICKMAX * sin (2.0 * 3.14159265 * NUD * (GP.turn - TURNS) + PHID);
      DRIFT_Pass (x, L / 2.);
    }
}

HACDIP::HACDIP (string name, double l, double hkickmax, double nud,
		double phid, int turns, int turne):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("HACDIP");
      L = l;
      HKICKMAX = hkickmax;
      NUD = nud;
      PHID = phid;
      TURNS = turns;
      TURNE = turne;
    }
  else
    {
      cout << "Error: check signs of L. " << endl;
      exit (1);
    }
}
void
HACDIP::SetP (const char *name, double value)
{
  if (strcmp (name, "HKICKMAX") == 0)
    {
      HKICKMAX = value;
    }
  else if (strcmp (name, "NUD") == 0)
    {
      NUD = value;
    }
  else if (strcmp (name, "PHID") == 0)
    {
      PHID = value;
    }
  else if (strcmp (name, "TURNS") == 0)
    {
      TURNS = int (value);
    }
  else if (strcmp (name, "TURNE") == 0)
    {
      TURNE = int (value);
    }

  else
    {
      cout << "HACDIP does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
HACDIP::GetP (const char *name)
{
  if (strcmp (name, "HKICKMAX") == 0)
    {
      return HKICKMAX;
    }
  else if (strcmp (name, "TURNS") == 0)
    {
      return TURNS;
    }
  else if (strcmp (name, "TURNE") == 0)
    {
      return TURNE;
    }
  else if (strcmp (name, "NUD") == 0)
    {
      return NUD;
    }
  else if (strcmp (name, "PHID") == 0)
    {
      return PHID;
    }
  else
    {
      cout << "HACDIP does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
HACDIP::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  ACDIP_Pass (x, L, HKICKMAX, 0., NUD, TURNS, TURNE, PHID);
  LtoG (x, DX, DY, DT);
}

void
HACDIP::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//------------------------------VACDIP------------------------------------------
VACDIP::VACDIP (string name, double l, double vkickmax, double nud,
		double phid, int turns, int turne):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("VACDIP");
      L = l;
      VKICKMAX = vkickmax;
      NUD = nud;
      PHID = phid;
      TURNS = turns;
      TURNE = turne;

    }
  else
    {
      cout << "Error: check signs of L. " << endl;
      exit (1);
    }
}
void
VACDIP::SetP (const char *name, double value)
{
  if (strcmp (name, "VKICKMAX") == 0)
    {
      VKICKMAX = value;
    }
  else if (strcmp (name, "NUD") == 0)
    {
      NUD = value;
    }
  else if (strcmp (name, "PHID") == 0)
    {
      PHID = value;
    }
  else if (strcmp (name, "TURNS") == 0)
    {
      TURNS = int (value);
    }
  else if (strcmp (name, "TURNE") == 0)
    {
      TURNE = int (value);
    }
  else
    {
      cout << "VACDIP does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
VACDIP::GetP (const char *name)
{
  if (strcmp (name, "VKICKMAX") == 0)
    {
      return VKICKMAX;
    }
  else if (strcmp (name, "TURNS") == 0)
    {
      return TURNS;
    }
  else if (strcmp (name, "TURNE") == 0)
    {
      return TURNE;
    }
  else if (strcmp (name, "NUD") == 0)
    {
      return NUD;
    }
  else if (strcmp (name, "PHID") == 0)
    {
      return PHID;
    }
  else
    {
      cout << "VACDIP does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
VACDIP::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  ACDIP_Pass (x, L, 0., VKICKMAX, NUD, TURNS, TURNE, PHID);
  LtoG (x, DX, DY, DT);
}

void
VACDIP::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//----------------------------ACMULT ( single order, order >= quadrupole )-----------------
void
ACMULT_Pass (double x[6], double L, int Norder, double KL, int TTURNS,
	     double PHI0)
{
  int i;
  long int fac = 1;
  double Xn, Yn, Xn0, Yn0;
  double KNL, KNSL;
  double By, Bx;

  DRIFT_Pass (x, L / 2.0);

  if (Norder < 0)
    {
      KNL = 0.;
      KNSL = KL;
    }
  else
    {
      KNL = KL;
      KNSL = 0.;
    }

  By = 0.;
  Bx = 0.;
  Xn = 1.;
  Yn = 0.;

  for (i = 1; i < Norder + 1; i++)
    {
      Xn0 = Xn;
      Yn0 = Yn;
      Xn = Xn0 * x[x_] - Yn0 * x[y_];
      Yn = Xn0 * x[y_] + Yn0 * x[x_];
      fac = fac * i;
    }
  By = By + (KNL * Xn - KNSL * Yn) / fac;
  Bx = Bx + (KNL * Yn + KNSL * Xn) / fac;

  x[px_] = x[px_] - By * sin (2.0 * M_PI * GP.turn / TTURNS + PHI0);
  x[py_] = x[py_] + Bx * sin (2.0 * M_PI * GP.turn / TTURNS + PHI0);

  DRIFT_Pass (x, L / 2.0);
}

ACMULT::ACMULT (string name, double l, int norder, double kl, int tturns,
		double phi0):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("ACMULT");
      L = l;
      Norder = norder;
      KL = kl;
      TTURNS = tturns;
      PHI0 = phi0;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
ACMULT::SetP (const char *name, double value)
{
  if (strcmp (name, "KL") == 0)
    {
      KL = value;
    }
  else if (strcmp (name, "Norder") == 0)
    {
      if (value == 0)
	{
	  cout << " Use ACKICK for Norder = 0." << endl;
	  exit (0);
	}
      else if (abs (value) > 10)
	{
	  cout << " The maximum order of ACMULT is 10." << endl;
	  exit (1);
	}
      else
	{
	  Norder = int (value);
	}
    }
  else if (strcmp (name, "TTURNS") == 0)
    {
      TTURNS = int (value);
    }
  else if (strcmp (name, "PHI0") == 0)
    {
      PHI0 = value;
    }
  else
    {
      cout << "ACMULT does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
ACMULT::GetP (const char *name)
{
  if (strcmp (name, "KL") == 0)
    {
      return KL;
    }
  else if (strcmp (name, "Norder") == 0)
    {
      return Norder;
    }
  else if (strcmp (name, "TTURNS") == 0)
    {
      return TTURNS;
    }
  else if (strcmp (name, "PHI0") == 0)
    {
      return PHI0;
    }
  else
    {
      cout << "ACMULT  does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
ACMULT::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  ACMULT_Pass (x, L, Norder, KL, TTURNS, PHI0);
  LtoG (x, DX, DY, DT);
}

void
ACMULT::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//---------------------------COOLING--------------------------------
template < class T > void
COOLING_Pass (T x[6], double L, double ALPHA)
{
  int i;
  DRIFT_Pass (x, L / 2.);
  for (i = 0; i < 4; i++)
    x[i] = (1. - ALPHA) * x[i];
  DRIFT_Pass (x, L / 2.);
}

COOLING::COOLING (string name, double l, double alpha):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("COOLING");
      L = l;
      ALPHA = alpha;
    }
  else
    {
      cout << "Error: check signs of L and S of COOLING element." << endl;
      exit (1);
    }
}
void
COOLING::SetP (const char *name, double value)
{
  if (strcmp (name, "ALPHA") == 0)
    {
      ALPHA = value;
    }
  else
    {
      cout << "COOLING does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
COOLING::GetP (const char *name)
{
  if (strcmp (name, "ALPHA") == 0)
    {
      return ALPHA;
    }
  else
    {
      cout << "COOLING does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
COOLING::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  COOLING_Pass (x, L, ALPHA);
  LtoG (x, DX, DY, DT);
}

void
COOLING::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//---------------------------DIFFUSE--------------------------------
template <class T> 
void DIFFUSE_Pass(T x[6], double L,  double DIFF_X, double  DIFF_Y, double DIFF_DELTA)
{
  DRIFT_Pass(x, L/2.);   
  x[1]= x[1]+ DIFF_X * rnd(seed);
  x[3]= x[3]+ DIFF_Y * rnd(seed);
  x[5]= x[5]+ DIFF_DELTA *rnd(seed);
  DRIFT_Pass(x, L/2.);
}

  DIFFUSE::DIFFUSE(string name, double l, double diff_x, double diff_y, double diff_delta ): Element(name)
    { 
      if (l >=0 )
	{
	  TYPE=string("DIFFUSE");
	  L=l;
	  DIFF_X = diff_x;  DIFF_Y = diff_y; DIFF_DELTA = diff_delta;
	}
      else
	{
	  cout<<"Error: check signs of L and S of DIFFUSE element."<<endl;
          exit(1);
	}
    }
    void  DIFFUSE::SetP(const char *name, double value) 
    {
      if (strcmp( name, "DIFF_X" ) == 0) 
	{
	 DIFF_X = value;
	}
      else if (strcmp( name, "DIFF_Y" ) == 0) 
	{
	 DIFF_Y = value;
	}
      else if (strcmp( name, "DIFF_delta" ) == 0) 
	{
	 DIFF_DELTA = value;
	}
      else 
	{
	  cout<<"DIFFUSE does not have  parameter of "<<name<<endl; 
          exit(0); 
	} 
    }
  double DIFFUSE::GetP(const char *name) 
  {
    if (strcmp( name, "DIFF_X" ) == 0) 
      {
	return DIFF_X;
      }
    else if (strcmp( name, "DIFF_Y" ) == 0) 
      {
	return DIFF_Y;
      }
    else if (strcmp( name, "DIFF_DELTA" ) == 0) 
      {
	return DIFF_DELTA;
      }
    else 
      {
	cout<<"DIFFUSE does not have  parameter of "<<name<<endl; 
	exit(0); 
      } 
  }
  void   DIFFUSE::Pass(double x[6]){
    GtoL(x,DX,DY,DT);   DIFFUSE_Pass(x, L, DIFF_X, DIFF_Y, DIFF_DELTA);   LtoG(x,DX,DY,DT); 
  }
  void   DIFFUSE::DAPass(tps x[6]){ 
    DRIFT_Pass(x, L);
  }

//---------------------------------BPM-------------------------
BPM::BPM (string name, double l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("BPM");
      L = l;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
BPM::SetP (const char *name, double value)
{
  cout << "No parameter to be set for BPM." << endl;
  exit (1);
}

double
BPM::GetP (const char *name)
{
  cout << "No parameter to be returned for BPM." << endl;
  exit (1);
}

void
BPM::Pass (double x[6])
{
  DRIFT_Pass (x, L);
}

void
BPM::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//----------------------------HBPM-----------------------------------
HBPM::HBPM (string name, double l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("HBPM");
      L = l;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
HBPM::SetP (const char *name, double value)
{
  cout << "No parameter to be set for HBPM." << endl;
  exit (1);
}

double
HBPM::GetP (const char *name)
{
  cout << "No parameter to be returned for HBPM." << endl;
  exit (1);
}

void
HBPM::Pass (double x[6])
{
  DRIFT_Pass (x, L);
}

void
HBPM::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//---------------------------------VBPM------------------------------
VBPM::VBPM (string name, double l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("VBPM");
      L = l;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
VBPM::SetP (const char *name, double value)
{
  cout << "No parameter to be set for VBPM." << endl;
  exit (1);
}

double
VBPM::GetP (const char *name)
{
  cout << "No parameter to be returned for VBPM." << endl;
  exit (0);
}

void
VBPM::Pass (double x[6])
{
  DRIFT_Pass (x, L);
}

void
VBPM::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

MARKER::MARKER (string name, double l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("MARKER");
      L = l;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
MARKER::SetP (const char *name, double value)
{
  cout << "No parameter to be set for MARKER." << endl;
  exit (1);
}

double
MARKER::GetP (const char *name)
{
  cout << "No parameter to be returned for MARKER." << endl;
  exit (1);
}

void
MARKER::Pass (double x[6])
{
  DRIFT_Pass (x, L);
}

void
MARKER::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//---------------------------RCOLL--------------------------------
RCOLL::RCOLL (string name, double l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("RCOLL");
      L = l;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
RCOLL::SetP (const char *name, double value)
{
  cout << "No parameter to be set for MARKER." << endl;
  exit (1);
}

double
RCOLL::GetP (const char *name)
{
  cout << "No parameter to be returned for MARKER." << endl;
  exit (1);
}

void
RCOLL::Pass (double x[6])
{
  DRIFT_Pass (x, L);
}

void
RCOLL::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

ECOLL::ECOLL (string name, double l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("ECOLL");
      L = l;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
ECOLL::SetP (const char *name, double value)
{
  cout << "No parameter to be set for MARKER." << endl;
  exit (1);
}

double
ECOLL::GetP (const char *name)
{
  cout << "No parameter to be returned for MARKER." << endl;
  exit (1);
}

void
ECOLL::Pass (double x[6])
{
  DRIFT_Pass (x, L);
}

void
ECOLL::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//----------------------INSTR------------------------------------
INSTR::INSTR (string name, double l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("INSTR");
      L = l;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
INSTR::SetP (const char *name, double value)
{
  cout << "No parameter is set for INSTR." << endl;
  exit (0);
}

double
INSTR::GetP (const char *name)
{
  cout << "No parameter is returned for INSTR." << endl;
  exit (1);
}

void
INSTR::Pass (double x[6])
{
  DRIFT_Pass (x, L);
}

void
INSTR::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//----------------------------SOLEN-----------------------------------
template < class T > void
SOLEN_Pass (T x[6], double L, double KS, double M[36])
{
  int i, j, k;
  T xtemp[6];

  if (KS == 0)
    {
      DRIFT_Pass (x, L);
    }
  else
    {
      x[1] = x[1] / (1 + x[5]);	// x'=px / (1+delta)
      x[3] = x[3] / (1 + x[5]);
      for (i = 0; i < 6; i++)
	xtemp[i] = x[i];

      for (k = 0; k < 6; k++)
	{
	  x[k] = 0.;
	  for (j = 0; j < 6; j++)
	    x[k] = x[k] + M[k * 6 + j] * xtemp[j];
	}
      x[1] = x[1] * (1 + x[5]);	// px=x' * (1+delta)
      x[3] = x[3] * (1 + x[5]);
    }
}

SOLEN::SOLEN (string name, double l, double ks):
Element (name)
{
  if (l > 0 && ks != 0.)
    {
      TYPE = string ("SOLEN");
      L = l;
      KS = ks;

      int i, j;
      double g = KS / 2, theta = KS * L / 2, costheta =
	cos (theta), sintheta = sin (theta);

      if (KS != 0.)
	{
	  for (i = 0; i < 6; i++)
	    for (j = 0; j < 6; j++)
	      M[i * 6 + j] = 0.0;

	  M[0 * 6 + 0] = costheta * costheta;
	  M[0 * 6 + 1] = sintheta * costheta / g;
	  M[0 * 6 + 2] = sintheta * costheta;
	  M[0 * 6 + 3] = sintheta * sintheta / g;

	  M[1 * 6 + 0] = -g * sintheta * costheta;
	  M[1 * 6 + 1] = costheta * costheta;
	  M[1 * 6 + 2] = -g * sintheta * sintheta;
	  M[1 * 6 + 3] = sintheta * costheta;

	  M[2 * 6 + 0] = -sintheta * costheta;
	  M[2 * 6 + 1] = -sintheta * sintheta / g;
	  M[2 * 6 + 2] = costheta * costheta;
	  M[2 * 6 + 3] = sintheta * costheta / g;

	  M[3 * 6 + 0] = g * sintheta * sintheta;
	  M[3 * 6 + 1] = -sintheta * costheta;
	  M[3 * 6 + 2] = -g * sintheta * costheta;
	  M[3 * 6 + 3] = costheta * costheta;

	  M[4 * 6 + 4] = 1.0;
	  M[5 * 6 + 5] = 1.0;
	}
    }
  else if (l > 0 && ks == 0.)
    {
      TYPE = string ("SOLEN");
      L = l;
      KS = 0;
    }
  else
    {
      cout << "Solenoid: L must be positive. " << endl;
      exit (1);
    }
}
void
SOLEN::SetP (const char *name, double value)
{
  if (strcmp (name, "KS") == 0)
    {
      KS = value;
      int i, j;
      double g = KS / 2, theta = KS * L / 2, costheta =
	cos (theta), sintheta = sin (theta);

      if (KS != 0.)
	{
	  for (i = 0; i < 6; i++)
	    for (j = 0; j < 6; j++)
	      M[i * 6 + j] = 0.0;

	  M[0 * 6 + 0] = costheta * costheta;
	  M[0 * 6 + 1] = sintheta * costheta / g;
	  M[0 * 6 + 2] = sintheta * costheta;
	  M[0 * 6 + 3] = sintheta * sintheta / g;

	  M[1 * 6 + 0] = -g * sintheta * costheta;
	  M[1 * 6 + 1] = costheta * costheta;
	  M[1 * 6 + 2] = -g * sintheta * sintheta;
	  M[1 * 6 + 3] = sintheta * costheta;

	  M[2 * 6 + 0] = -sintheta * costheta;
	  M[2 * 6 + 1] = -sintheta * sintheta / g;
	  M[2 * 6 + 2] = costheta * costheta;
	  M[2 * 6 + 3] = sintheta * costheta / g;

	  M[3 * 6 + 0] = g * sintheta * sintheta;
	  M[3 * 6 + 1] = -sintheta * costheta;
	  M[3 * 6 + 2] = -g * sintheta * costheta;
	  M[3 * 6 + 3] = costheta * costheta;

	  M[4 * 6 + 4] = 1.0;
	  M[5 * 6 + 5] = 1.0;
	}
      else
	{
	  KS = 0.;
	}
    }
  else
    {
      cout << "SOLEN does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
SOLEN::GetP (const char *name)
{
  if (strcmp (name, "KS") == 0)
    {
      return KS;
    }
  else
    {
      cout << "SOLEN does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
SOLEN::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  SOLEN_Pass (x, L, KS, M);
  LtoG (x, DX, DY, DT);
}

void
SOLEN::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  SOLEN_Pass (x, L, KS, M);
  LtoG (x, DX, DY, DT);
}

//-----------------------MATRIX-----------------------------------
template < typename T > void
MATRIX_Pass (T x[6], double XCO_IN[6], double XCO_OUT[6], double M66[36])
{
  int i, j;
  T x1[6];

  for (i = 0; i < 6; i++)
    x[i] = x[i] - XCO_IN[i];

  for (i = 0; i < 6; i++)
    {
      x1[i] = 0.0;
      for (j = 0; j < 6; j++)
	x1[i] = x1[i] + M66[i * 6 + j] * x[j];
    }
  for (i = 0; i < 6; i++)
    x[i] = x1[i];

  for (i = 0; i < 6; i++)
    x[i] = x[i] + XCO_OUT[i];
}

MATRIX::MATRIX (string name, double l, double m66[36], double xco_in[6],
		double xco_out[6]):
Element (name)
{
  int i;
  TYPE = string ("MATRIX");
  L = l;
  for (i = 0; i < 36; i++)
    M66[i] = m66[i];
  for (i = 0; i < 6; i++)
    XCO_IN[i] = xco_in[i];
  for (i = 0; i < 6; i++)
    XCO_OUT[i] = xco_out[i];
}

void
MATRIX::SetP (const char *name, double value)
{
  if (strcmp (name, "M11") == 0)
    {
      M66[0] = value;
    }
  else if (strcmp (name, "M12") == 0)
    {
      M66[1] = value;
    }
  else if (strcmp (name, "M13") == 0)
    {
      M66[2] = value;
    }
  else if (strcmp (name, "M14") == 0)
    {
      M66[3] = value;
    }
  else if (strcmp (name, "M15") == 0)
    {
      M66[4] = value;
    }
  else if (strcmp (name, "M16") == 0)
    {
      M66[5] = value;
    }
  else if (strcmp (name, "M21") == 0)
    {
      M66[6] = value;
    }
  else if (strcmp (name, "M22") == 0)
    {
      M66[7] = value;
    }
  else if (strcmp (name, "M23") == 0)
    {
      M66[8] = value;
    }
  else if (strcmp (name, "M24") == 0)
    {
      M66[9] = value;
    }
  else if (strcmp (name, "M25") == 0)
    {
      M66[10] = value;
    }
  else if (strcmp (name, "M26") == 0)
    {
      M66[11] = value;
    }
  else if (strcmp (name, "M31") == 0)
    {
      M66[12] = value;
    }
  else if (strcmp (name, "M32") == 0)
    {
      M66[13] = value;
    }
  else if (strcmp (name, "M33") == 0)
    {
      M66[14] = value;
    }
  else if (strcmp (name, "M34") == 0)
    {
      M66[15] = value;
    }
  else if (strcmp (name, "M35") == 0)
    {
      M66[16] = value;
    }
  else if (strcmp (name, "M36") == 0)
    {
      M66[17] = value;
    }
  else if (strcmp (name, "M41") == 0)
    {
      M66[18] = value;
    }
  else if (strcmp (name, "M42") == 0)
    {
      M66[19] = value;
    }
  else if (strcmp (name, "M43") == 0)
    {
      M66[20] = value;
    }
  else if (strcmp (name, "M44") == 0)
    {
      M66[21] = value;
    }
  else if (strcmp (name, "M45") == 0)
    {
      M66[22] = value;
    }
  else if (strcmp (name, "M46") == 0)
    {
      M66[23] = value;
    }
  else if (strcmp (name, "M51") == 0)
    {
      M66[24] = value;
    }
  else if (strcmp (name, "M52") == 0)
    {
      M66[25] = value;
    }
  else if (strcmp (name, "M53") == 0)
    {
      M66[26] = value;
    }
  else if (strcmp (name, "M54") == 0)
    {
      M66[27] = value;
    }
  else if (strcmp (name, "M55") == 0)
    {
      M66[28] = value;
    }
  else if (strcmp (name, "M56") == 0)
    {
      M66[29] = value;
    }
  else if (strcmp (name, "M61") == 0)
    {
      M66[30] = value;
    }
  else if (strcmp (name, "M62") == 0)
    {
      M66[31] = value;
    }
  else if (strcmp (name, "M63") == 0)
    {
      M66[32] = value;
    }
  else if (strcmp (name, "M64") == 0)
    {
      M66[33] = value;
    }
  else if (strcmp (name, "M65") == 0)
    {
      M66[34] = value;
    }
  else if (strcmp (name, "M66") == 0)
    {
      M66[35] = value;
    }
  else if (strcmp (name, "XCO_IN_X") == 0)
    {
      XCO_IN[0] = value;
    }
  else if (strcmp (name, "XCO_IN_PX") == 0)
    {
      XCO_IN[1] = value;
    }
  else if (strcmp (name, "XCO_IN_Y") == 0)
    {
      XCO_IN[2] = value;
    }
  else if (strcmp (name, "XCO_IN_PY") == 0)
    {
      XCO_IN[3] = value;
    }
  else if (strcmp (name, "XCO_IN_Z") == 0)
    {
      XCO_IN[4] = value;
    }
  else if (strcmp (name, "XCO_IN_DELTA") == 0)
    {
      XCO_IN[5] = value;
    }
  else if (strcmp (name, "XCO_OUT_X") == 0)
    {
      XCO_OUT[0] = value;
    }
  else if (strcmp (name, "XCO_OUT_PX") == 0)
    {
      XCO_OUT[1] = value;
    }
  else if (strcmp (name, "XCO_OUT_Y") == 0)
    {
      XCO_OUT[2] = value;
    }
  else if (strcmp (name, "XCO_OUT_PY") == 0)
    {
      XCO_OUT[3] = value;
    }
  else if (strcmp (name, "XCO_OUT_Z") == 0)
    {
      XCO_OUT[4] = value;
    }
  else if (strcmp (name, "XCO_OUT_DELTA") == 0)
    {
      XCO_OUT[5] = value;
    }
  else
    {
      cout << "Matrix does not have parameter of  " << name << endl;
      exit (0);
    }
}
double
MATRIX::GetP (const char *name)
{
  if (strcmp (name, "M11") == 0)
    {
      return M66[0];
    }
  else if (strcmp (name, "M12") == 0)
    {
      return M66[1];
    }
  else if (strcmp (name, "M13") == 0)
    {
      return M66[2];
    }
  else if (strcmp (name, "M14") == 0)
    {
      return M66[3];
    }
  else if (strcmp (name, "M15") == 0)
    {
      return M66[4];
    }
  else if (strcmp (name, "M16") == 0)
    {
      return M66[5];
    }
  else if (strcmp (name, "M21") == 0)
    {
      return M66[6];
    }
  else if (strcmp (name, "M22") == 0)
    {
      return M66[7];
    }
  else if (strcmp (name, "M23") == 0)
    {
      return M66[8];
    }
  else if (strcmp (name, "M24") == 0)
    {
      return M66[9];
    }
  else if (strcmp (name, "M25") == 0)
    {
      return M66[10];
    }
  else if (strcmp (name, "M26") == 0)
    {
      return M66[11];
    }
  else if (strcmp (name, "M31") == 0)
    {
      return M66[12];
    }
  else if (strcmp (name, "M32") == 0)
    {
      return M66[13];
    }
  else if (strcmp (name, "M33") == 0)
    {
      return M66[14];
    }
  else if (strcmp (name, "M34") == 0)
    {
      return M66[15];
    }
  else if (strcmp (name, "M35") == 0)
    {
      return M66[16];
    }
  else if (strcmp (name, "M36") == 0)
    {
      return M66[17];
    }
  else if (strcmp (name, "M41") == 0)
    {
      return M66[18];
    }
  else if (strcmp (name, "M42") == 0)
    {
      return M66[19];
    }
  else if (strcmp (name, "M43") == 0)
    {
      return M66[20];
    }
  else if (strcmp (name, "M44") == 0)
    {
      return M66[21];
    }
  else if (strcmp (name, "M45") == 0)
    {
      return M66[22];
    }
  else if (strcmp (name, "M46") == 0)
    {
      return M66[23];
    }
  else if (strcmp (name, "M51") == 0)
    {
      return M66[24];
    }
  else if (strcmp (name, "M52") == 0)
    {
      return M66[25];
    }
  else if (strcmp (name, "M53") == 0)
    {
      return M66[26];
    }
  else if (strcmp (name, "M54") == 0)
    {
      return M66[27];
    }
  else if (strcmp (name, "M55") == 0)
    {
      return M66[28];
    }
  else if (strcmp (name, "M56") == 0)
    {
      return M66[29];
    }
  else if (strcmp (name, "M61") == 0)
    {
      return M66[30];
    }
  else if (strcmp (name, "M62") == 0)
    {
      return M66[31];
    }
  else if (strcmp (name, "M63") == 0)
    {
      return M66[32];
    }
  else if (strcmp (name, "M64") == 0)
    {
      return M66[33];
    }
  else if (strcmp (name, "M65") == 0)
    {
      return M66[34];
    }
  else if (strcmp (name, "M66") == 0)
    {
      return M66[35];
    }
  else if (strcmp (name, "XCO_IN_X") == 0)
    {
      return XCO_IN[0];
    }
  else if (strcmp (name, "XCO_IN_PX") == 0)
    {
      return XCO_IN[1];
    }
  else if (strcmp (name, "XCO_IN_Y") == 0)
    {
      return XCO_IN[2];
    }
  else if (strcmp (name, "XCO_IN_PY") == 0)
    {
      return XCO_IN[3];
    }
  else if (strcmp (name, "XCO_IN_Z") == 0)
    {
      return XCO_IN[4];
    }
  else if (strcmp (name, "XCO_IN_DELTA") == 0)
    {
      return XCO_IN[5];
    }
  else if (strcmp (name, "XCO_OUT_X") == 0)
    {
      return XCO_OUT[0];
    }
  else if (strcmp (name, "XCO_OUT_PX") == 0)
    {
      return XCO_OUT[1];
    }
  else if (strcmp (name, "XCO_OUT_Y") == 0)
    {
      return XCO_OUT[2];
    }
  else if (strcmp (name, "XCO_OUT_PY") == 0)
    {
      return XCO_OUT[3];
    }
  else if (strcmp (name, "XCO_OUT_Z") == 0)
    {
      return XCO_OUT[4];
    }
  else if (strcmp (name, "XCO_OUT_DELTA") == 0)
    {
      return XCO_OUT[5];
    }
  else
    {
      cout << "Matrix does not have parameter of  " << name << endl;
      exit (0);
    }
}
void
MATRIX::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  MATRIX_Pass (x, XCO_IN, XCO_OUT, M66);
  LtoG (x, DX, DY, DT);
}

void
MATRIX::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  MATRIX_Pass (x, XCO_IN, XCO_OUT, M66);
  LtoG (x, DX, DY, DT);
}

//--------------------ELSEP----------------------------------
ELSEP::ELSEP (string name, double l):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("ELSEP");
      L = l;
    }
  else
    {
      cout << "Error: check signs of L and S" << endl;
      exit (1);
    }
}
void
ELSEP::SetP (const char *name, double value)
{
  cout << "No parameter to be set for ELSEP." << endl;
  exit (1);
}

double
ELSEP::GetP (const char *name)
{
  cout << "No parameter to be returned for ELSEP." << endl;
  exit (1);
}

void
ELSEP::Pass (double x[6])
{
  DRIFT_Pass (x, L);
}

void
ELSEP::DAPass (tps x[6])
{
  DRIFT_Pass (x, L);
}

//--------------------OFFSET (coordinate system change)-------------------------------
template < typename T > void
OFFSET_Pass (T x[6], double DX, double DPX, double DY, double DPY,
	     double TILT)
{
  int i;
  double cosT = cos (-TILT), sinT = sin (-TILT);
  T xtemp[6];
  x[0] = x[0] - DX;
  x[1] = x[1] - DPX;
  x[2] = x[2] - DY;
  x[3] = x[3] - DPY;
  for (i = 0; i < 6; i++)
    xtemp[i] = x[i];
  x[0] = cosT * xtemp[0] + sinT * xtemp[2];
  x[1] = cosT * xtemp[1] + sinT * xtemp[3];
  x[2] = cosT * xtemp[2] - sinT * xtemp[0];
  x[3] = cosT * xtemp[3] - sinT * xtemp[1];
}

OFFSET::OFFSET (string name, double l, double dx, double dpx, double dy,
		double dpy, double tilt):
Element (name)
{
  if (l == 0)
    {
      TYPE = string ("OFFSET");
      L = l;
      DX = dx;
      DPX = dpx;
      DY = dy;
      DPY = dpy;
      TILT = tilt;
    }
  else
    {
      cout << "Error: L for OFFSET must be zero." << endl;
      exit (1);
    }
}
void
OFFSET::SetP (const char *name, double value)
{
  if (strcmp (name, "DX") == 0)
    {
      DX = value;
    }
  else if (strcmp (name, "DPX") == 0)
    {
      DPX = value;
    }
  if (strcmp (name, "DY") == 0)
    {
      DY = value;
    }
  else if (strcmp (name, "DPY") == 0)
    {
      DPY = value;
    }
  else if (strcmp (name, "TILT") == 0)
    {
      TILT = value;
    }
  else
    {
      cout << "OFFSET does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
OFFSET::GetP (const char *name)
{
  if (strcmp (name, "DX") == 0)
    {
      return DX;
    }
  else if (strcmp (name, "DPX") == 0)
    {
      return DPX;
    }
  if (strcmp (name, "DY") == 0)
    {
      return DY;
    }
  else if (strcmp (name, "DPY") == 0)
    {
      return DPY;
    }
  else if (strcmp (name, "TILT") == 0)
    {
      return TILT;
    }
  else
    {
      cout << "OFFSET does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
OFFSET::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  OFFSET_Pass (x, DX, DPX, DY, DPY, TILT);
  LtoG (x, DX, DY, DT);
}

void
OFFSET::DAPass (tps x[6])
{
  GtoL (x, DX, DY, DT);
  OFFSET_Pass (x, DX, DPX, DY, DPY, TILT);
  LtoG (x, DX, DY, DT);
}

//----------------------------RFCAV--------------------------------
template < class T > void
RFCAV_Pass (T x[6], double L, double VRF, double FRF, double PHASE0)
{
  DRIFT_Pass (x, L / 2.);
  x[5] =
    x[5] +
    (VRF * GP.Q / GP.A / GP.energy) * sin (2.0 * M_PI * FRF * x[z_] / 3.0e8);
  DRIFT_Pass (x, L / 2.);
}

RFCAV::RFCAV (string name, double l, double vrf, double frf, double phase0):
Element (name)
{
  TYPE = string ("RFCAV");
  VRF = vrf;
  FRF = frf;
  PHASE0 = phase0;
  L = l;
}

void
RFCAV::SetP (const char *name, double value)
{
  if (strcmp (name, "VRF") == 0)
    {
      VRF = value;
    }
  else if (strcmp (name, "FRF") == 0)
    {
      FRF = value;
    }
  else if (strcmp (name, "PHASE0") == 0)
    {
      PHASE0 = value;
    }
  else
    {
      cout << "RFCAV does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
RFCAV::GetP (const char *name)
{
  if (strcmp (name, "VRF") == 0)
    {
      return VRF;
    }
  else if (strcmp (name, "FRF") == 0)
    {
      return FRF;
    }
  else if (strcmp (name, "PHASE0") == 0)
    {
      return PHASE0;
    }
  else
    {
      cout << "RFCAV does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
RFCAV::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  RFCAV_Pass (x, L, VRF, FRF, PHASE0);
  LtoG (x, DX, DY, DT);
}

void
RFCAV::DAPass (tps x[6])
{
  if (GP.twiss_6d == 1)
    {
      GtoL (x, DX, DY, DT);
      RFCAV_Pass (x, L, VRF, FRF, PHASE0);
      LtoG (x, DX, DY, DT);
    }
  else
    {
      DRIFT_Pass (x, L);
    }
}

//----------------------BEAMBEAM------------------------------
void
cerrf (double xx, double yy, double &wx, double &wy)
{
  int n, nc, nu;
  double x, y, q, h, xl, xh, yh, tx, ty, tn, sx, sy, saux;
  double rx[33], ry[33];
  double cc = 1.12837916709551, zero = 0, one = 1, two = 2, half = .5,
    xlim = 5.33, ylim = 4.29, fac1 = 3.2, fac2 = 23, fac3 = 21;

  x = abs (xx);
  y = abs (yy);
  if (y < ylim and x < xlim)
    {
      q = (one - y / ylim) * sqrt (one - (x / xlim) * (x / xlim));
      h = one / (fac1 * q);
      nc = 7 + int (fac2 * q);
      xl = pow (h, (1 - nc));
      xh = y + half / h;
      yh = x;
      nu = 10 + int (fac3 * q);
      rx[nu + 1] = zero;
      ry[nu + 1] = zero;
      for (n = nu; n >= 1; n--)
	{
	  tx = xh + n * rx[n + 1];
	  ty = yh - n * ry[n + 1];
	  tn = tx * tx + ty * ty;
	  rx[n] = half * tx / tn;
	  ry[n] = half * ty / tn;
	}
      sx = zero;
      sy = zero;
      for (n = nc; n >= 1; n--)
	{
	  saux = sx + xl;
	  sx = rx[n] * saux - ry[n] * sy;
	  sy = rx[n] * sy + ry[n] * saux;
	  xl = h * xl;
	}
      wx = cc * sx;
      wy = cc * sy;
    }
  else
    {
      xh = y;
      yh = x;
      rx[1] = zero;
      ry[1] = zero;
      for (n = 9; n >= 1; n--)
	{
	  tx = xh + n * rx[1];
	  ty = yh - n * ry[1];
	  tn = tx * tx + ty * ty;
	  rx[1] = half * tx / tn;
	  ry[1] = half * ty / tn;
	}
      wx = cc * rx[1];
      wy = cc * ry[1];
    }
  if (yy < zero)
    {
      wx = two * exp (y * y - x * x) * cos (two * x * y) - wx;
      wy = -two * exp (y * y - x * x) * sin (two * x * y) - wy;
      if (xx > zero)
	wy = -wy;
    }
  else
    {
      if (xx < zero)
	wy = -wy;
    }
}

template < class T >
void BB4D (T x, T y, double gamma, double N, double sigmax,
      double sigmay, double &Dpx, double &Dpy)
{
  double rp = 1.534698e-18;
  double x1, y1, x2, y2, signx, signy;

  signx = ((x >= 0.0) ? (1) : (-1));
  signy = ((y >= 0.0) ? (1) : (-1));

  if (x == 0. && y == 0.)
    {
      Dpx = 0.;
      Dpy = 0.;
    }
  else
    {
      if (abs (sigmax - sigmay) / sigmax < 1.0e-6)
	{
	  double r2 = x * x + y * y;
	  double temp1 = 2 * N * rp / gamma / r2;
	  double temp2 = 1 - exp (-r2 / 2 / sigmax / sigmax);
	  Dpx = temp1 * x * temp2;
	  Dpy = temp1 * y * temp2;
	}
      else if (sigmax > sigmay)
	{
	  double temp1 =
	    (2. * N * rp / gamma) * sqrt (M_PI / 2. /
					  (sigmax * sigmax -
					   sigmay * sigmay));
	  double temp2 =
	    exp (-x * x / 2. / sigmax / sigmax -
		 y * y / 2. / sigmay / sigmay);
	  double temp3 = sqrt (2. * (sigmax * sigmax - sigmay * sigmay));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs (x) / temp3;
	  y1 = fabs (y) / temp3;
	  x2 = fabs (x) * sigmay / sigmax / temp3;
	  y2 = fabs (y) * sigmax / sigmay / temp3;
	  cerrf (x1, y1, w1_real, w1_imag);
	  cerrf (x2, y2, w2_real, w2_imag);
	  Dpy = temp1 * (w1_real - temp2 * w2_real) * signy;
	  Dpx = temp1 * (w1_imag - temp2 * w2_imag) * signx;
	}
      else if (sigmax < sigmay)
	{
	  double temp1 =
	    (2 * N * rp / gamma) * sqrt (M_PI / 2 /
					 (sigmay * sigmay - sigmax * sigmax));
	  double temp2 =
	    exp (-x * x / 2 / sigmax / sigmax - y * y / 2 / sigmay / sigmay);
	  double temp3 = sqrt (2 * (sigmay * sigmay - sigmax * sigmax));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs (x) / temp3;
	  y1 = fabs (y) / temp3;
	  x2 = fabs (x) * sigmay / sigmax / temp3;
	  y2 = fabs (y) * sigmax / sigmay / temp3;
	  cerrf (y1, x1, w1_real, w1_imag);
	  cerrf (y2, x2, w2_real, w2_imag);
	  Dpy = temp1 * (w1_imag - temp2 * w2_imag) * signy;
	  Dpx = temp1 * (w1_real - temp2 * w2_real) * signx;
	}
    }
}

template < class T >
void BB6D(T x[6], double gamma, double Np, double sigma_l, int N_slice, 
                       double emitx_rms,  double betax_star, double alfx_star,
	               double emity_rms,  double betay_star, double alfy_star)
 // emitx_rms, emity_rms:  un-normalized rms emittance,  sigma=SQRT[ emitx_rms * betax ]  
{
  int i,j; 
  double rp= 1.534698e-18;
  double x0[6];
  double Nsigma=4, Np_slice[N_slice], z_star[N_slice];
  double S, gx_star, gy_star, betax, alfx, betay, alfy;
  double sigmax, sigmay, dsigmax2ds, dsigmay2ds, dUdsigmax2, dUdsigmay2; 
  double X, Y,  Dpx, Dpy, temp1, temp2, temp3, temp4;
  
  //----center of each slice of strong bunch, each slice has not same particle population.
  
  for(i=0;i<N_slice;i++)
    z_star[i]= -1.0*Nsigma + (Nsigma*2.0/N_slice)*(2*i+1)/2;
  
  if(N_slice == 11 ) {
    Np_slice[0]=Np*  0.000500905  ;
    Np_slice[1]=Np*  0.0049242    ;
    Np_slice[2]=Np*  0.0290614    ;
    Np_slice[3]=Np*  0.103138     ;
    Np_slice[4]=Np*  0.220408     ;
    Np_slice[5]=Np*  0.28387      ;
    Np_slice[6]=Np*  0.220408     ;
    Np_slice[7]=Np*  0.103138     ;
    Np_slice[8]=Np*  0.0290614    ;
    Np_slice[9]=Np*  0.0049242    ;
    Np_slice[10]=Np* 0.000500905  ;  
  }
  else{
    for(i=0;i<N_slice;i++)
      Np_slice[i]=Np* (gsl_cdf_ugaussian_P( z_star[i]+ Nsigma*2.0/N_slice/2 )- gsl_cdf_ugaussian_P( z_star[i]- Nsigma*2.0/N_slice/2  ) );
  }  
  
  for(i=0; i<N_slice;i++ ) z_star[i]=z_star[i]*sigma_l;
  
  //for(i=0;i<N_slice;i++) cout << z_star[i]/sigma_l<<"   "<< Np_slice[i]/Np<<endl;

  //-----calculate the changes in x[6]
  for(i=0; i<N_slice;i++) {

    for (j=0;j<6;j++) x0[j]=x[j];
    S=(x0[4]-z_star[i])/2.0;

    gx_star= (1.0 + alfx_star * alfx_star )/ betax_star;
    gy_star= (1.0 + alfy_star * alfy_star )/ betay_star;
    betax= ( 1.0/gx_star + gx_star *(S-alfx_star/gx_star) *(S-alfx_star/gx_star) );
    alfx =(alfx_star -gx_star * S );
    betay= ( 1.0/gy_star + gy_star *(S-alfy_star/gy_star) *(S-alfy_star/gy_star) );
    alfy =(alfy_star -gy_star * S );

    sigmax  =  sqrt( emitx_rms * betax ) ;
    sigmay  =  sqrt( emity_rms * betay ) ;
    dsigmax2ds =-2.0*alfx* ( emitx_rms );
    dsigmay2ds =-2.0*alfy* ( emity_rms );
    
    //----calculate the kicks for each slice
    X=x0[0] + x0[1]*S;
    Y=x0[2] + x0[3]*S;
    BB4D(X, Y, gamma, Np_slice[i], sigmax, sigmay, Dpx, Dpy);
    temp1 = 1.0/2./(sigmax*sigmax-sigmay*sigmay) ;
    temp2 = X *  Dpx + Y *  Dpy;
    temp3 = 2.* Np_slice[i] *rp /gamma ;
    temp4 = exp ( - X * X / 2. / sigmax /sigmax -  Y * Y / 2. / sigmay /sigmay ) ;
    
    x[0]=x0[0] - S * Dpx * GP.bbscale ;
    x[1]=x0[1] + Dpx  * GP.bbscale;   
    x[2]=x0[2] - S * Dpy * GP.bbscale ;
    x[3]=x0[3] + Dpy * GP.bbscale;   
    x[4]=x0[4];
    
    if (sigmax == sigmay ) 
      {
      	x[5]=x0[5]+ 0.5 * Dpx * GP.bbscale * ( x0[1] + 0.5* Dpx * GP.bbscale ) 
                  + 0.5*  Dpy * GP.bbscale * ( x0[3] + 0.5* Dpy * GP.bbscale ) 
                  + (1.0/sigmax/sigmax ) * dsigmax2ds  * ( temp3 * GP.bbscale/ 4. ) * temp4;
      }
    else if (sigmax > sigmay ) 
      {
	dUdsigmax2 =  temp1 * GP.bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	dUdsigmay2 = -temp1 * GP.bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	x[5]= x0[5] + 0.5 * Dpx * GP.bbscale * ( x0[1] + 0.5* Dpx * GP.bbscale) 
                    + 0.5 * Dpy *GP.bbscale  * ( x0[3] + 0.5* Dpy *  GP.bbscale) 
                    - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
    else
      {
	dUdsigmax2 = -temp1 *  GP.bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	dUdsigmay2 =  temp1 *  GP.bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	x[5]= x0[5] + 0.5 * Dpx *GP.bbscale * ( x0[1] + 0.5* Dpx* GP.bbscale) 
                    + 0.5 * Dpy *GP.bbscale *  ( x0[3] + 0.5* Dpy * GP.bbscale) 
                    - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
  }
}

void BEAMBEAM_Pass (double x[6], int TREATMENT, double NP, double SIGMAL, int NSLICE,
	       double EMITX, double BETAX, double ALFAX, double EMITY,
	       double BETAY, double ALFAY)
{
  if (NP != 0.)
    {
      if (int (TREATMENT) == 6)
	{
	  BB6D (x, GP.gamma, NP, SIGMAL, NSLICE, EMITX, BETAX, ALFAX, EMITY,
		BETAY, ALFAY);
	}
      else
	{
	  double Dpx, Dpy;
	  BB4D (x[0], x[2], GP.gamma, NP, sqrt (EMITX * BETAX),
		sqrt (EMITY * BETAY), Dpx, Dpy);
	  x[1] = x[1] + Dpx * GP.bbscale;
	  x[3] = x[3] + Dpy * GP.bbscale;
	}
    }
}

BEAMBEAM::BEAMBEAM (string name, int treatment, double np, double sigmal,
		    int nslice, double emitx, double emity, double betax,
		    double alfax, double betay, double alfay):
Element (name)
 //     emitx, emity:  un-normalized rms emittance,  sigma=SQRT[ emitx * betax ]  
{
  TREATMENT = treatment;
  NP = np;
  SIGMAL = sigmal;
  NSLICE = nslice;
  EMITX = emitx;
  EMITY = emity;
  BETAX = betax;
  ALFAX = alfax;
  BETAY = betay;
  ALFAY = alfay;
  TYPE = string ("BEAMBEAM");
  L = 0.;
}

void
BEAMBEAM::SetP (const char *name, double value)
{
  if (strcmp (name, "TREATMENT") == 0)
    {
      int temp;
      temp = int (value);
      if (temp == 4 || temp == 6)
	{
	  TREATMENT = temp;
	}
      else
	{
	  cout << "Error: beambeam treatment only can be 4 or 6" << endl;
	  exit (1);
	}
    }
  else if (strcmp (name, "NP") == 0)
    {
      NP = value;
    }
  else if (strcmp (name, "SIGMAL") == 0)
    {
      SIGMAL = value;
    }
  else if (strcmp (name, "NSLICE") == 0)
    {
      NSLICE = int (value);
    }
  else if (strcmp (name, "EMITX") == 0)
    {
      EMITX = value;
    }
  else if (strcmp (name, "EMITY") == 0)
    {
      EMITY = value;
    }
  else if (strcmp (name, "BETAX") == 0)
    {
      BETAX = value;
    }
  else if (strcmp (name, "ALFAX") == 0)
    {
      ALFAX = value;
    }
  else if (strcmp (name, "BETAY") == 0)
    {
      BETAY = value;
    }
  else if (strcmp (name, "ALFAY") == 0)
    {
      ALFAY = value;
    }
  else
    {
      cout << "BEAMBEAM does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
BEAMBEAM::GetP (const char *name)
{
  if (strcmp (name, "TREATMENT") == 0)
    {
      return TREATMENT;
    }
  else if (strcmp (name, "NP") == 0)
    {
      return NP;
    }
  else if (strcmp (name, "SIGMAL") == 0)
    {
      return SIGMAL;
    }
  else if (strcmp (name, "NSLICE") == 0)
    {
      return NSLICE;
    }
  else if (strcmp (name, "EMITX") == 0)
    {
      return EMITX;
    }
  else if (strcmp (name, "EMITY") == 0)
    {
      return EMITY;
    }
  else if (strcmp (name, "BETAX") == 0)
    {
      return BETAX;
    }
  else if (strcmp (name, "ALFAX") == 0)
    {
      return ALFAX;
    }
  else if (strcmp (name, "BETAY") == 0)
    {
      return BETAY;
    }
  else if (strcmp (name, "ALFAY") == 0)
    {
      return ALFAY;
    }
  else
    {
      cout << "BEAMBEAM does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
BEAMBEAM::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  BEAMBEAM_Pass (x, TREATMENT, NP, SIGMAL, NSLICE, EMITX, BETAX, ALFAX, EMITY,
		 BETAY, ALFAY);
  LtoG (x, DX, DY, DT);
}

void
BEAMBEAM::DAPass (tps x[6])
{
  double sigmax, sigmay;
  double ksi_x, ksi_y;
  double rp = 1.534698e-18;
  sigmax = sqrt (EMITX * BETAX);
  sigmay = sqrt (EMITY * BETAY);
  ksi_x = 2.0 * NP * rp / sigmax / GP.gamma / (sigmax + sigmay);
  ksi_y = 2.0 * NP * rp / sigmay / GP.gamma / (sigmax + sigmay);
  GtoL (x, DX, DY, DT);
  x[1] = x[1] + ksi_x * x[0] * GP.bbscale;
  x[3] = x[3] + ksi_y * x[2] * GP.bbscale;
  LtoG (x, DX, DY, DT);
}

//---------------------LRBB ( long-range BEAMBEAM )------------
void
LRBB_Pass (double x[6], double gamma, double NP, double SEPX, double SEPY,
	   double SIGMAX, double SIGMAY, double KICKX0, double KICKY0)
//  here we assume the particle's coordinates from the strong beam center is ( SEPX+x[1], SEPY+x[3] )
{
  double dpx1, dpy1;
  BB4D (SEPX + x[1], SEPY + x[3], GP.gamma, NP, SIGMAX, SIGMAY, dpx1, dpy1);
  x[1] = x[1] + (dpx1 - KICKX0);
  x[3] = x[3] + (dpy1 - KICKY0);
}

LRBB::LRBB (string name, double np, double sepx, double sepy, double sigmax,
	    double sigmay):
Element (name)
 //     sepx, sepy are the offset of the weak beam center w.r.t to the center strong beam  
{
  NP = np;
  SEPX = sepx;
  SEPY = sepy;
  SIGMAX = sigmax;
  SIGMAY = sigmay;
  TYPE = string ("LRBB");
  L = 0.;
  BB4D (SEPX, SEPY, GP.gamma, NP, SIGMAX, SIGMAY, KICKX0, KICKY0);
}

void
LRBB::SetP (const char *name, double value)
{
  if (strcmp (name, "NP") == 0)
    {
      NP = value;
    }
  else if (strcmp (name, "SEPX") == 0)
    {
      SEPX = value;
    }
  else if (strcmp (name, "SEPY") == 0)
    {
      SEPY = value;
    }
  else if (strcmp (name, "SIGMAX") == 0)
    {
      SIGMAX = value;
    }
  else if (strcmp (name, "SIGMAY") == 0)
    {
      SIGMAY = value;
    }
  else
    {
      cout << "LRBB does not have parameter of " << name << endl;
      exit (0);
    }
  BB4D (SEPX, SEPY, GP.gamma, NP, SIGMAX, SIGMAY, KICKX0, KICKY0);
}

double
LRBB::GetP (const char *name)
{
  if (strcmp (name, "NP") == 0)
    {
      return NP;
    }
  else if (strcmp (name, "SEPX") == 0)
    {
      return SEPX;
    }
  else if (strcmp (name, "SEPY") == 0)
    {
      return SEPY;
    }
  else if (strcmp (name, "SIGMAX") == 0)
    {
      return SIGMAX;
    }
  else if (strcmp (name, "SIGMAY") == 0)
    {
      return SIGMAY;
    }
  else if (strcmp (name, "KICKX0") == 0)
    {
      return KICKX0;
    }
  else if (strcmp (name, "KICKY0") == 0)
    {
      return KICKY0;
    }
  else
    {
      cout << "LRBB does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
LRBB::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  LRBB_Pass (x, GP.gamma, NP, SEPX, SEPY, SIGMAX, SIGMAY, KICKX0, KICKY0);
  LtoG (x, DX, DY, DT);
}

void
LRBB::DAPass (tps x[6])
{
  //no effect to any parameters.
}

//----------------------E-LENS------------------------------
void ELENS_Pass (double x[6], double Ne, double Le, double beta_e, int N_slice,
	    double sigmax, double sigmay)
{
  int i;
  double Le_slice = Le * 1.0 / N_slice;
  double Ne_slice = Ne * 1.0 / N_slice;
  double Dpx, Dpy;

  for (i = 0; i < N_slice; i++)
    {
      DRIFT_Pass (x, Le_slice / 2.);
      BB4D (x[0], x[2], GP.gamma, Ne_slice, sigmax, sigmay, Dpx, Dpy);
      x[1] = x[1] - Dpx * (1. + beta_e);
      x[3] = x[3] - Dpy * (1. + beta_e);
      DRIFT_Pass (x, Le_slice / 2.);
    }
}

void
elens_pass_round_Gaussian_topoff (double x[6], double gamma,
				  double Ne, double Le, double beta_e,
				  int N_slice, double sigmax, double sigmay)
// I assume sigmax=sigmay, top off [-sigmax*a, sigmax*a ], Ne from perfect Gaussian
{
  int i;
  double Ne_slice = Ne * 1.0 / N_slice;
  double Le_slice = Le * 1.0 / N_slice;
  double a = 0.4;
  double scale;
  double r, rp = 1.534698e-18;

  for (i = 0; i < N_slice; i++)
    {
      DRIFT_Pass (x, Le_slice / 2.);
      r = sqrt (x[0] * x[0] + x[2] * x[2]);
      if (r <= a * sigmax)
	{
	  scale =
	    -2 * Ne_slice * (exp (-a * a / 2) * (r / sigmax) * (r / sigmax) /
			     2) * rp / gamma;
	  x[1] = x[1] + scale * x[0] / r / r;
	  x[3] = x[3] + scale * x[2] / r / r;
	}
      else
	{
	  scale =
	    -2 * Ne_slice * (1 - exp (-r * r / 2 / sigmax / sigmax) -
			     (1 - exp (-a * a / 2)) +
			     exp (-a * a / 2) * a * a / 2) * rp / gamma;
	  x[1] = x[1] + scale * x[0] / r / r;
	  x[3] = x[3] + scale * x[2] / r / r;
	}
      DRIFT_Pass (x, Le_slice / 2.);
    }
}

void
elens_pass_round_Gaussian_truncated (double x[6], double gamma,
				     double Ne, double Le, double beta_e,
				     int N_slice, double sigmax,
				     double sigmay)
// I assume sigmax=sigmay, Gaussian tail cut off from Nc*sigmax
{
  int i;
  double Le_slice = Le * 1.0 / N_slice;
  double Ne_slice = Ne * 1.0 / N_slice;
  double Dpx, Dpy;
  double Nc = 100, Nsigma;

  for (i = 0; i < N_slice; i++)
    {
      DRIFT_Pass (x, Le_slice / 2.);
      Nsigma = sqrt (x[0] * x[0] + x[2] * x[2]) / sigmax;
      if (Nsigma <= Nc)
	{
	  BB4D (x[0], x[2], gamma, Ne_slice, sigmax, sigmay, Dpx, Dpy);
	}
      else
	{
	  BB4D (x[0] * Nc / Nsigma, x[2] * Nc / Nsigma, gamma, Ne_slice,
		sigmax, sigmay, Dpx, Dpy);
	  Dpx = Dpx * Nc / Nsigma;
	  Dpy = Dpy * Nc / Nsigma;
	}
      x[1] = x[1] - Dpx * (1. + beta_e);
      x[3] = x[3] - Dpy * (1. + beta_e);
      DRIFT_Pass (x, Le_slice / 2.);
    }
}

void
elens_pass_round_uniform (double x[6], double gamma,
			  double Ne, double Le, double beta_e, int N_slice,
			  double sigmax, double sigmay)
// I assume round uniform distribution in r<=a
{
  int i;
  double Le_slice = Le * 1.0 / N_slice;
  double Ne_slice = Ne * 1.0 / N_slice;
  double a = 0.31e-3 * 5;
  double scale, r2;
  double rp = 1.534698e-18;

  for (i = 0; i < N_slice; i++)
    {
      DRIFT_Pass (x, Le_slice / 2.);
      scale = -2 * Ne_slice * rp * (1 + beta_e) / GP.gamma;
      r2 = (x[0] * x[0] + x[2] * x[2]);
      if (sqrt (r2) < a)
	{
	  x[1] = x[1] + scale * x[0] / a / a;
	  x[3] = x[3] + scale * x[2] / a / a;
	}
      else
	{
	  x[1] = x[1] + scale * x[0] / r2;
	  x[3] = x[3] + scale * x[2] / r2;
	}
      DRIFT_Pass (x, Le_slice / 2.);
    }
}

ELENS::ELENS (string name, double l, double ne, int nslice, double betae,
	      double sigmax, double sigmay):
Element (name)
{
  if (l >= 0)
    {
      TYPE = string ("ELENS");
      L = l;
      NE = ne;
      NSLICE = nslice;
      BETAE = betae;
      SIGMAX = sigmax;
      SIGMAY = sigmay;
    }
  else
    {
      cout << "Error: check signs of L. " << endl;
      exit (1);
    }
}
void
ELENS::SetP (const char *name, double value)
{
  if (strcmp (name, "NE") == 0)
    {
      NE = value;
    }
  else if (strcmp (name, "NSLICE") == 0)
    {
      NSLICE = int (value);
    }
  else if (strcmp (name, "BETAE") == 0)
    {
      BETAE = value;
    }
  else if (strcmp (name, "SIGMAX") == 0)
    {
      SIGMAX = value;
    }
  else if (strcmp (name, "SIGMAY") == 0)
    {
      SIGMAY = value;
    }
  else
    {
      cout << "ELENS does not have  parameter of " << name << endl;
      exit (0);
    }
}
double
ELENS::GetP (const char *name)
{
  if (strcmp (name, "NE") == 0)
    {
      return NE;
    }
  else if (strcmp (name, "NSLICE") == 0)
    {
      return NSLICE;
    }
  else if (strcmp (name, "BETAE") == 0)
    {
      return BETAE;
    }
  else if (strcmp (name, "SIGMAX") == 0)
    {
      return SIGMAX;
    }
  else if (strcmp (name, "SIGMAY") == 0)
    {
      return SIGMAY;
    }
  else
    {
      cout << "ELENS does not have  parameter of " << name << endl;
      exit (0);
    }
}
void
ELENS::Pass (double x[6])
{
  GtoL (x, DX, DY, DT);
  ELENS_Pass (x, NE, L, BETAE, NSLICE, SIGMAX, SIGMAY);
  LtoG (x, DX, DY, DT);
}

void
ELENS::DAPass (tps x[6])
{
  int i;
  double ksi_x, ksi_y;
  double rp = 1.534698e-18;
  ksi_x =
    -2.0 * (1.0 * NE / NSLICE) * rp * (1 +
				       BETAE) / SIGMAX / GP.gamma / (SIGMAX +
								     SIGMAY);
  ksi_y =
    -2.0 * (1.0 * NE / NSLICE) * rp * (1 +
				       BETAE) / SIGMAY / GP.gamma / (SIGMAX +
								     SIGMAY);
  GtoL (x, DX, DY, DT);
  for (i = 0; i < NSLICE; i++)
    {
      DRIFT_Pass (x, L * 1.0 / NSLICE / 2.0);
      x[1] = x[1] + ksi_x * x[0];
      x[3] = x[3] + ksi_y * x[2];
      DRIFT_Pass (x, L * 1.0 / NSLICE / 2.0);
    }
  LtoG (x, DX, DY, DT);
}


//=========================================
//
//            Line  class
//
//=========================================
Line::Line ()
{
  Ncell = 0;
  Length = 0;
  Tune1 = 0.;
  Tune2 = 0.;
  Tune3 = 0.;
  Chromx1 = 0.;
  Chromy1 = 0.;
  Chromx2 = 0.;
  Chromy2 = 0.;
  Chromx3 = 0.;
  Chromy3 = 0.;
}

void
Line::Update ()
{
  int i;
  double spointer = 0.;
  for (i = 0; i < Cell.size (); i++)
    {
      spointer = Cell[i]->L + spointer;
      Cell[i]->S = spointer;
    }
  Ncell = Cell.size ();
  Length = spointer;
  frev0 = (light_speed * GP.beta) / Length;
  frf = frev0 * GP.harm;
}

void
Line::Append (Element * x)
{
  Cell.push_back (x);
  Update ();
}

void
Line::Delete (int i)
{
  vector < Element * >::iterator start;
  start = Cell.begin () + i;
  Cell.erase (start);
  Update ();
}

void
Line::Insert (int i, Element * temp)
{
  vector < Element * >::iterator start;
  start = Cell.begin () + i;
  Cell.insert (start, 1, temp);
  Update ();
}

void
Line::Empty ()
{
  vector < Element * >::iterator start;
  vector < Element * >::iterator end;
  start = Cell.begin ();
  end = Cell.end ();
  Cell.erase (start, end);
  Update ();
}

void
Line_Rewind (Line & linename, int k)
{
  int i;
  Element *new_element;
  Line temp_line;

  for (i = k; i < linename.Ncell; i++)
    {
      new_element = linename.Cell[i];
      temp_line.Append (new_element);
    }
  for (i = 0; i < k; i++)
    {
      new_element = linename.Cell[i];
      temp_line.Append (new_element);
    }
  linename = temp_line;
}

void
Line_Invert (Line & linename)
{
  int i;
  Element *new_element;
  Line temp_line;

  for (i = 0; i < linename.Ncell; i++)
    {
      new_element = linename.Cell[linename.Ncell - 1 - i];
      temp_line.Append (new_element);
    }
  linename = temp_line;
}

void
Line_Repeat (Line & linename1, Line & linename2, int n)
{
  int i, j;
  Element *new_element;

  for (j = 0; j < n; j++)
    {
      for (i = 0; i < linename2.Ncell; i++)
	{
	  new_element = linename2.Cell[i];
	  linename1.Append (new_element);
	}
    }
}

void
Line_Connect (Line & linename1, Line & linename2, Line & linename3)
{
  int i;
  Element *new_element;

  for (i = 0; i < linename2.Ncell; i++)
    {
      new_element = linename2.Cell[i];
      linename1.Append (new_element);
    }
  for (i = 0; i < linename3.Ncell; i++)
    {
      new_element = linename3.Cell[i];
      linename1.Append (new_element);
    }
}

int
Get_Index (Line & linename, const char *name, int k)
{
  int i;
  int count = 0;
  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->NAME == string (name))
	{
	  count++;
	  if (count == k)
	    {
	      return i;
	    }
	}
    }
  if (count == 0)
    {
      cout << "Element not found." << endl;
      return 0;
    }
  return 0;
}

void
Split_Quad (Line & linename, int i, int m)
{
  int j;
  string old_name;
  double length, k1l;
  Element *new_element;
  char index[125];

  if (m < 1)
    {
      cout << " Number of split slices should be => 1 !" << endl;
      exit (1);
    }

  if (linename.Cell[i]->TYPE == string ("QUAD"))
    {
      old_name = linename.Cell[i]->NAME;
      length = linename.Cell[i]->L;
      k1l = linename.Cell[i]->GetP ("K1L");
      linename.Delete (i);
      for (j = 0; j < m; j++)
	{
	  sprintf (index, "_%d", j + 1);
	  new_element = new QUAD (old_name, length / m, k1l / m);
	  linename.Insert (i, new_element);
	}
    }
  else
    {
      cout << " The i-th element is not QUAD." << endl;
      exit (1);
    }
}

void
Split_Sext (Line & linename, int i, int m)
{
  int j;
  string old_name;
  double length, k2l;
  Element *new_element;
  char index[125];
  string new_name;

  if (m < 1)
    {
      cout << " Number of split slices should be => 1 !" << endl;
      exit (1);
    }

  if (linename.Cell[i]->TYPE == string ("SEXT"))
    {
      old_name = linename.Cell[i]->NAME;
      length = linename.Cell[i]->L;
      k2l = linename.Cell[i]->GetP ("K2L");
      linename.Delete (i);
      for (j = 0; j < m; j++)
	{
	  sprintf (index, "_%d", j + 1);
	  new_element = new SEXT (old_name, length / m, k2l / m);
	  linename.Insert (i, new_element);
	}
    }
  else
    {
      cout << " The i-th element is not SEXT." << endl;
      exit (1);
    }
}

void
Split_Quad_Sext (Line & linename)
{
  int i;

  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->TYPE == string ("SEXT"))
	Split_Sext (linename, i, int (linename.Cell[i]->L / 0.05) + 1);
    }

  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->TYPE == string ("QUAD"))
	Split_Quad (linename, i, int (linename.Cell[i]->L / 0.05) + 1);
    }
}

void
Split_Sbend (Line & linename, int i, int m)
{
  int j;
  string old_name;
  double length, angle, e1, e2;
  Element *new_element;
  char index[125];

  if (m < 1)
    {
      cout << " Number of split slices should be => 1 !" << endl;
      exit (1);
    }

  if (linename.Cell[i]->TYPE == string ("SBEND"))
    {
      old_name = linename.Cell[i]->NAME;
      length = linename.Cell[i]->L;
      angle = linename.Cell[i]->GetP ("ANGLE");
      e1 = linename.Cell[i]->GetP ("E1");
      e2 = linename.Cell[i]->GetP ("E2");
      linename.Delete (i);
      for (j = 0; j < m; j++)
	{
	  sprintf (index, "_%d", j + 1);
	  new_element = new SBEND (old_name, length / m, angle / m, e1, e2);
	  linename.Insert (i, new_element);
	}
    }
  else
    {
      cout << " The i-th element is not SBEND ." << endl;
      exit (1);
    }
}

void
Split_Drift (Line & linename, int i, int m)
{
  int j;
  string old_name;
  double length;
  Element *new_element;
  char index[125];

  if (m < 1)
    {
      cout << " Number of split slices should be => 1 !" << endl;
      exit (1);
    }

  if (linename.Cell[i]->TYPE == string ("DRIFT"))
    {
      old_name = linename.Cell[i]->NAME;
      length = linename.Cell[i]->L;
      linename.Delete (i);
      for (j = 0; j < m; j++)
	{
	  sprintf (index, "_%d", j + 1);
	  new_element = new DRIFT (old_name, length / m);
	  linename.Insert (i, new_element);
	}
    }
  else
    {
      cout << " The i-th element is not DRIFT." << endl;
      exit (1);
    }
}

void
Split_Mult (Line & linename, int i, int m)
{
  int j;
  string old_name;
  double length;
  double knl[11], knsl[11];
  Element *new_element;
  char index[125];

  if (m < 1)
    {
      cout << " Number of split slices should be => 1 !" << endl;
      exit (1);
    }

  if (linename.Cell[i]->TYPE == string ("DRIFT"))
    {
      old_name = linename.Cell[i]->NAME;
      length = linename.Cell[i]->L;

      for (j = 0; j < 11; j++)
	{
	  char name1[125], name2[125];
	  sprintf (name1, "K%dL", j);
	  sprintf (name2, "K%dSL", j);
	  knl[j] = linename.Cell[i]->GetP (name1) / m;
	  knsl[j] = linename.Cell[i]->GetP (name2) / m;
	}

      linename.Delete (i);
      for (j = 0; j < m; j++)
	{
	  sprintf (index, "_%d", j + 1);
	  new_element = new MULT (old_name, length / m, knl, knsl);
	  linename.Insert (i, new_element);
	}
    }
  else
    {
      cout << " The i-th element is not MULT." << endl;
      exit (1);
    }
}

//----function: concate adjacient DRIFT, doesn't change Twiss
void
Concat_Drift (Line & linename)
{
  int i;
  double length;

  i = 1;
  while (i < linename.Ncell)
    {
      if (linename.Cell[i - 1]->TYPE ==
	  string ("DRIFT") and linename.Cell[i]->TYPE == string ("DRIFT"))
	{
	  length = linename.Cell[i]->L;
	  linename.Delete (i);
	  linename.Cell[i - 1]->L = linename.Cell[i - 1]->L + length;
	  linename.Update ();
	}
      else
	{
	  i++;
	}
    }
}

//----get rid of zero strength element, doesn't change Twiss
void
Clean_Up (Line & linename)
{
  int i;
  double length;
  Element *temp_element;

  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->TYPE ==
	  string ("SEXT") and linename.Cell[i]->GetP ("K2L") == 0)
	{
	  length = linename.Cell[i]->L;
	  linename.Delete (i);
	  temp_element = new DRIFT ("TEMPD", length);
	  linename.Insert (i, temp_element);
	}
    }

  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->TYPE ==
	  string ("OCT") and linename.Cell[i]->GetP ("K3L") == 0)
	{
	  length = linename.Cell[i]->L;
	  linename.Delete (i);
	  temp_element = new DRIFT ("TEMPD", length);
	  linename.Insert (i, temp_element);
	}
    }

  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->TYPE ==
	  string ("SOLEN") and linename.Cell[i]->GetP ("KS") == 0)
	{
	  length = linename.Cell[i]->L;
	  linename.Delete (i);
	  temp_element = new DRIFT ("TEMPD", length);
	  linename.Insert (i, temp_element);
	}
    }

  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->TYPE ==
	  string ("RFCAV") and linename.Cell[i]->GetP ("VRF") == 0)
	{
	  length = linename.Cell[i]->L;
	  linename.Delete (i);
	  temp_element = new DRIFT ("TEMPD", length);
	  linename.Insert (i, temp_element);
	}
    }

  for (i = 0; i < linename.Ncell; i++)
    if (linename.Cell[i]->TYPE ==
	string ("KICKER") and linename.Cell[i]->GetP ("HKICK") ==
	0. and linename.Cell[i]->GetP ("VKICK") == 0.)
      {
	length = linename.Cell[i]->L;
	linename.Delete (i);
	if (length != 0.)
	  {
	    temp_element = new DRIFT ("TEMPD", length);
	    linename.Insert (i, temp_element);
	  }
      }

  for (i = 0; i < linename.Ncell; i++)
    if (linename.Cell[i]->TYPE ==
	string ("HKICKER") and linename.Cell[i]->GetP ("HKICK") == 0.)
      {
	length = linename.Cell[i]->L;
	linename.Delete (i);
	if (length != 0.)
	  {
	    temp_element = new DRIFT ("TEMPD", length);
	    linename.Insert (i, temp_element);
	  }
      }

  for (i = 0; i < linename.Ncell; i++)
    if (linename.Cell[i]->TYPE ==
	string ("VKICKER") and linename.Cell[i]->GetP ("VKICK") == 0.)
      {
	length = linename.Cell[i]->L;
	linename.Delete (i);
	if (length != 0.)
	  {
	    temp_element = new DRIFT ("TEMPD", length);
	    linename.Insert (i, temp_element);
	  }
      }

  for (i = 0; i < linename.Ncell; i++)
    if (linename.Cell[i]->TYPE ==
	string ("MARKER") or linename.Cell[i]->TYPE ==
	string ("BPM") or linename.Cell[i]->TYPE ==
	string ("HBPM") or linename.Cell[i]->TYPE ==
	string ("VBPM") or linename.Cell[i]->TYPE ==
	string ("HBPM") or linename.Cell[i]->TYPE == string ("VBPM"))
      {
	length = linename.Cell[i]->L;
	linename.Delete (i);
	temp_element = new DRIFT ("TEMPD", length);
	linename.Insert (i, temp_element);
      }

  for (i = 0; i < linename.Ncell; i++)
    if (linename.Cell[i]->TYPE == string ("MULT"))
      {
	length = linename.Cell[i]->L;
	int m;
	double kl = 0.;
	double knl[11], knsl[11];
	for (m = 0; m < 11; m++)
	  {
	    char name1[125], name2[125];
	    sprintf (name1, "K%dL", m);
	    sprintf (name2, "K%dSL", m);
	    knl[m] = linename.Cell[i]->GetP (name1);
	    knsl[m] = linename.Cell[i]->GetP (name2);
	  }
	for (m = 0; m < 11; m++)
	  kl = kl + abs (knl[m]) + abs (knsl[m]);

	if (kl == 0.)
	  {
	    linename.Delete (i);
	    temp_element = new DRIFT ("TEMPD", length);
	    linename.Insert (i, temp_element);
	  }
      }
  Concat_Drift (linename);
}

//----prepare for fast tarcking: will change Twiss
void
Make_Thin (Line & linename)
{
  int i;
  double length;
  Element *temp_element;

  //---convert Track_Fast() not-acceptable elements into drift
  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->TYPE !=
	  string ("DRIFT") & linename.Cell[i]->TYPE !=
	  string ("SBEND") & linename.Cell[i]->TYPE !=
	  string ("QUAD") & linename.Cell[i]->TYPE !=
	  string ("SEXT") & linename.Cell[i]->TYPE !=
	  string ("MULT") & linename.Cell[i]->TYPE !=
	  string ("BEAMBEAM") & linename.Cell[i]->TYPE !=
	  string ("ELENS") & linename.Cell[i]->TYPE !=
	  string ("RFCAV") & linename.Cell[i]->TYPE !=
	  string ("COOLING") & linename.Cell[i]->TYPE !=
	  string ("ACMULT") & linename.Cell[i]->TYPE !=
	  string ("MATRIX") & linename.Cell[i]->TYPE !=
	  string ("KICKER") & linename.Cell[i]->TYPE !=
	  string ("HKICKER") & linename.Cell[i]->TYPE !=
	  string ("VKICKER") & linename.Cell[i]->TYPE !=
	  string ("DIFFUSE") & linename.Cell[i]->TYPE !=
	  string ("GMULT") & linename.Cell[i]->TYPE != string ("LRBB"))
	{
	  length = linename.Cell[i]->L;
	  temp_element = new DRIFT ("TEMPD", length);
	  linename.Delete (i);
	  linename.Insert (i, temp_element);
	}
    }

  //-----further make thin of some elements  among  Track_Fast() acceptable elements
  i = 0;
  while (i < linename.Ncell - 1)
    {
      if (linename.Cell[i]->TYPE == string ("SEXT")
	  || linename.Cell[i]->TYPE == string ("MULT")
	  || linename.Cell[i]->TYPE == string ("HKICKER")
	  || linename.Cell[i]->TYPE == string ("VKICKER")
	  || linename.Cell[i]->TYPE == string ("RFCAV")
	  || linename.Cell[i]->TYPE == string ("ACMULT")
	  || linename.Cell[i]->TYPE == string ("COOLING")
	  || linename.Cell[i]->TYPE == string ("DIFFUSE")
	  || linename.Cell[i]->TYPE == string ("MATRIX")
	  || linename.Cell[i]->TYPE == string ("KICKER"))
	{
	  length = linename.Cell[i]->L;
	  if (length != 0.)
	    {
	      temp_element = new DRIFT ("TEMPD", length / 2.0);
	      linename.Insert (i, temp_element);
	      linename.Cell[i + 1]->L = 0.;
	      temp_element = new DRIFT ("TEMPD", length / 2.0);
	      linename.Insert (i + 2, temp_element);
	    }
	}
      i++;
    }

  //---even reduce the integration steps for SBEND  and QUAD to minimum 3
  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->TYPE == string ("SBEND"))
	{
	  if (linename.Cell[i]->GetP ("Nint") > 3)
	    linename.Cell[i]->SetP ("Nint", 3);
	}
    }

  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->TYPE == string ("QUAD"))
	{
	  if (linename.Cell[i]->GetP ("Nint") > 3)
	    linename.Cell[i]->SetP ("Nint", 3);
	}
    }

  //---concat all created drifts in above processes
  Concat_Drift (linename);
}

double
Get_KL (Line & linename, const char *name, const char *kl)
{
  int i;
  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->NAME == string (name))
	{
	  return linename.Cell[i]->GetP (kl);
	}
    }
  return 0.;
}

void
Set_KL (Line & linename, const char *name, const char *kl, double strength)
{
  int i;
  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->NAME == string (name))
	linename.Cell[i]->SetP (kl, strength);
    }
}

void
Set_dKL (Line & linename, const char *name, const char *kl, double dstrength)
{
  int i;
  for (i = 0; i < linename.Ncell; i++)
    {
      if (linename.Cell[i]->NAME == string (name))
	{
	  linename.Cell[i]->SetP (kl,
				  linename.Cell[i]->GetP (kl) + dstrength);
	}
    }
}

//=========================================
//
//            Magnete Alignment
//
//=========================================
template < typename T > void
GtoL (T x[6], double DX, double DY, double DT)
{

  int i;
  T xtemp[6];

  if (fabs (DX) + fabs (DY) > 1.0e-10)
    {
      x[0] = x[0] - DX;
      x[2] = x[2] - DY;
    }
  if (fabs (DT) > 1.0e-10)
    {
      double cosT = cos (DT), sinT = sin (DT);
      for (i = 0; i < 6; i++)
	xtemp[i] = x[i];
      x[0] = cosT * xtemp[0] + sinT * xtemp[2];
      x[1] = cosT * xtemp[1] + sinT * xtemp[3];
      x[2] = cosT * xtemp[2] - sinT * xtemp[0];
      x[3] = cosT * xtemp[3] - sinT * xtemp[1];
    }
}

template < typename T > void
LtoG (T x[6], double DX, double DY, double DT)
{
  int i;
  T xtemp[6];

  if (fabs (DT) > 1.0e-10)
    {
      double cosT = cos (DT), sinT = sin (DT);
      for (i = 0; i < 6; i++)
	xtemp[i] = x[i];
      x[0] = cosT * xtemp[0] - sinT * xtemp[2];
      x[1] = cosT * xtemp[1] - sinT * xtemp[3];
      x[2] = sinT * xtemp[0] + cosT * xtemp[2];
      x[3] = sinT * xtemp[1] + cosT * xtemp[3];
    }
  if (fabs (DX) + fabs (DY) > 1.0e-10)
    {
      x[0] = x[0] + DX;
      x[2] = x[2] + DY;
    }
}

//=================================================================
//
//    4th order Symplectic Integrator  (S.I.)
//
//=================================================================
template < class T > void
DRIFT_Pass (T x[6], double L)
{
  T u;
  if (L != 0.)
    {
      u = L / (1 + x[delta_]);
      x[x_] = x[x_] + x[px_] * u;
      x[y_] = x[y_] + x[py_] * u;
      x[z_] =
	x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								 x[delta_]);
    }
}

template < class T > void
bend_kick_pass (T x[6], double L, double href)
{
  x[px_] = x[px_] + (href * x[delta_] - href * href * x[x_]) * L;
  x[z_] = x[z_] - href * x[x_] * L;
}

template < class T > void
quad_kick_pass (T x[6], double k1l, double k1sl)
{
  int i, Norder;
  T Xn, Yn, Xn0, Yn0;
  T By, Bx;

  By = 0.;
  Bx = 0.;
  Xn = 1.;
  Yn = 0.;
  Norder = 1;
  for (i = 1; i < Norder + 1; i++)
    {
      Xn0 = Xn;
      Yn0 = Yn;
      Xn = Xn0 * x[x_] - Yn0 * x[y_];
      Yn = Xn0 * x[y_] + Yn0 * x[x_];
    }
  By = By + (k1l * Xn - k1sl * Yn);
  Bx = Bx + (k1l * Yn + k1sl * Xn);
  x[px_] = x[px_] - By;
  x[py_] = x[py_] + Bx;
}

template < class T > void
sext_kick_pass (T x[6], double k2l, double k2sl)
{
  int i, Norder;
  T Xn, Yn, Xn0, Yn0;
  T By, Bx;

  By = 0.;
  Bx = 0.;
  Xn = 1.;
  Yn = 0.;
  Norder = 2;
  for (i = 1; i < Norder + 1; i++)
    {
      Xn0 = Xn;
      Yn0 = Yn;
      Xn = Xn0 * x[x_] - Yn0 * x[y_];
      Yn = Xn0 * x[y_] + Yn0 * x[x_];
    }
  By = By + (k2l * Xn - k2sl * Yn) / 2;
  Bx = Bx + (k2l * Yn + k2sl * Xn) / 2;

  x[px_] = x[px_] - By;
  x[py_] = x[py_] + Bx;
}

template < class T > void
oct_kick_pass (T x[6], double k3l, double k3sl)
{
  int i, Norder;
  T Xn, Yn, Xn0, Yn0;
  T By, Bx;

  By = 0.;
  Bx = 0.;
  Xn = 1.;
  Yn = 0.;
  Norder = 3;
  for (i = 1; i < Norder + 1; i++)
    {
      Xn0 = Xn;
      Yn0 = Yn;
      Xn = Xn0 * x[x_] - Yn0 * x[y_];
      Yn = Xn0 * x[y_] + Yn0 * x[x_];
    }
  By = By + (k3l * Xn - k3sl * Yn) / 6;
  Bx = Bx + (k3l * Yn + k3sl * Xn) / 6;

  x[px_] = x[px_] - By;
  x[py_] = x[py_] + Bx;
}

template < class T > void
mult_kick_pass (T x[6], int Norder, double KNL[11], double KNSL[11])
{
  int i;
  long int fac = 1;
  T Xn, Yn, Xn0, Yn0;
  T By, Bx;

  By = KNL[0];
  Bx = KNSL[0];
  Xn = 1.;
  Yn = 0.;

  for (i = 1; i < Norder + 1; i++)
    {
      Xn0 = Xn;
      Yn0 = Yn;
      Xn = Xn0 * x[x_] - Yn0 * x[y_];
      Yn = Xn0 * x[y_] + Yn0 * x[x_];
      fac = fac * i;
      if (KNL[i] != 0. || KNSL[i] != 0.)
	{
	  By = By + (KNL[i] * Xn - KNSL[i] * Yn) / fac;
	  Bx = Bx + (KNL[i] * Yn + KNSL[i] * Xn) / fac;
	}
    }
  x[px_] = x[px_] - By;
  x[py_] = x[py_] + Bx;
}

template < class T > void
bend_mult_kick_pass (T x[6], double L, double href, int Norder,
		     double KNL[11], double KNSL[11])
{
  int i;
  long int fac = 1;
  T Xn, Yn, Xn0, Yn0;
  T By, Bx;

  By = KNL[0];
  Bx = KNSL[0];
  Xn = 1.;
  Yn = 0.;

  for (i = 1; i < Norder + 1; i++)
    {
      Xn0 = Xn;
      Yn0 = Yn;
      Xn = Xn0 * x[x_] - Yn0 * x[y_];
      Yn = Xn0 * x[y_] + Yn0 * x[x_];
      fac = fac * i;
      if (KNL[i] != 0. || KNSL[i] != 0.)
	{
	  By = By + (KNL[i] * Xn - KNSL[i] * Yn) / fac;
	  Bx = Bx + (KNL[i] * Yn + KNSL[i] * Xn) / fac;
	}
    }
  x[px_] = x[px_] - By + (href * x[delta_] - href * href * x[x_]) * L;;
  x[py_] = x[py_] + Bx;
  x[z_] = x[z_] - href * x[x_] * L;
}

template < class T > void
SBEND_Pass (T x[6], double L, int Nint, double Angle, double E1, double E2)
{
  int i;
  double href = Angle / L;
  double Lint = L / Nint;
  T u;

  x[1] = x[1] + tan (E1) * x[0] * href;
  x[3] = x[3] - tan (E1) * x[2] * href;

  for (i = 0; i < Nint; i++)
    {
      //DRIFT_Pass(x,Fdrift1*Lint);
      L = Fdrift1 * Lint;
      u = L / (1 + x[delta_]);
      x[x_] = x[x_] + x[px_] * u;
      x[y_] = x[y_] + x[py_] * u;
      x[z_] =
	x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								 x[delta_]);

      //bend_kick_pass(x, Fkick1*Lint, href);
      L = Fkick1 * Lint;
      x[px_] = x[px_] + (href * x[delta_] - href * href * x[x_]) * L;
      x[z_] = x[z_] - href * x[x_] * L;

      //DRIFT_Pass(x,Fdrift2*Lint);
      L = Fdrift2 * Lint;
      u = L / (1 + x[delta_]);
      x[x_] = x[x_] + x[px_] * u;
      x[y_] = x[y_] + x[py_] * u;
      x[z_] =
	x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								 x[delta_]);

      //bend_kick_pass(x, Fkick2*Lint, href);
      L = Fkick2 * Lint;
      x[px_] = x[px_] + (href * x[delta_] - href * href * x[x_]) * L;
      x[z_] = x[z_] - href * x[x_] * L;

      //DRIFT_Pass(x,Fdrift2*Lint);
      L = Fdrift2 * Lint;
      u = L / (1 + x[delta_]);
      x[x_] = x[x_] + x[px_] * u;
      x[y_] = x[y_] + x[py_] * u;
      x[z_] =
	x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								 x[delta_]);

      //bend_kick_pass(x, Fkick1*Lint, href);
      L = Fkick1 * Lint;
      x[px_] = x[px_] + (href * x[delta_] - href * href * x[x_]) * L;
      x[z_] = x[z_] - href * x[x_] * L;

      //DRIFT_Pass(x,Fdrift1*Lint); 
      L = Fdrift1 * Lint;
      u = L / (1 + x[delta_]);
      x[x_] = x[x_] + x[px_] * u;
      x[y_] = x[y_] + x[py_] * u;
      x[z_] =
	x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								 x[delta_]);
    }

  x[1] = x[1] + tan (E2) * x[0] * href;
  x[3] = x[3] - tan (E2) * x[2] * href;
}

template < class T > void
QUAD_Pass (T x[6], double L, int Nint, double k1l, double k1sl)
{
  int i;
  double Lint = L / Nint;

  if (L == 0.)
    {
      quad_kick_pass (x, k1l, k1sl);
    }
  else
    {
      double k1l_kick1, k1sl_kick1;
      double k1l_kick2, k1sl_kick2;
      k1l_kick1 = Fkick1 * k1l / Nint;
      k1sl_kick1 = Fkick1 * k1sl / Nint;
      k1l_kick2 = Fkick2 * k1l / Nint;
      k1sl_kick2 = Fkick2 * k1sl / Nint;

      int j, Norder = 1;
      T Xn, Yn, Xn0, Yn0;
      T By, Bx;
      T u;
      for (i = 0; i < Nint; i++)
	{

	  //DRIFT_Pass(x,Fdrift1*Lint);
	  L = Fdrift1 * Lint;
	  u = L / (1 + x[delta_]);
	  x[x_] = x[x_] + x[px_] * u;
	  x[y_] = x[y_] + x[py_] * u;
	  x[z_] =
	    x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								     x
								     [delta_]);

	  //quad_kick_pass(x, k1l_kick1, k1sl_kick1);
	  k1l = k1l_kick1;
	  k1sl = k1sl_kick1;
	  By = 0.;
	  Bx = 0.;
	  Xn = 1.;
	  Yn = 0.;
	  for (j = 1; j < Norder + 1; j++)
	    {
	      Xn0 = Xn;
	      Yn0 = Yn;
	      Xn = Xn0 * x[x_] - Yn0 * x[y_];
	      Yn = Xn0 * x[y_] + Yn0 * x[x_];
	    }
	  By = By + (k1l * Xn - k1sl * Yn);
	  Bx = Bx + (k1l * Yn + k1sl * Xn);
	  x[px_] = x[px_] - By;
	  x[py_] = x[py_] + Bx;

	  //DRIFT_Pass(x,Fdrift2*Lint);
	  L = Fdrift2 * Lint;
	  u = L / (1 + x[delta_]);
	  x[x_] = x[x_] + x[px_] * u;
	  x[y_] = x[y_] + x[py_] * u;
	  x[z_] =
	    x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								     x
								     [delta_]);

	  //quad_kick_pass(x, k1l_kick2, k1sl_kick2);
	  k1l = k1l_kick2;
	  k1sl = k1sl_kick2;
	  By = 0.;
	  Bx = 0.;
	  Xn = 1.;
	  Yn = 0.;
	  for (j = 1; j < Norder + 1; j++)
	    {
	      Xn0 = Xn;
	      Yn0 = Yn;
	      Xn = Xn0 * x[x_] - Yn0 * x[y_];
	      Yn = Xn0 * x[y_] + Yn0 * x[x_];
	    }
	  By = By + (k1l * Xn - k1sl * Yn);
	  Bx = Bx + (k1l * Yn + k1sl * Xn);
	  x[px_] = x[px_] - By;
	  x[py_] = x[py_] + Bx;

	  //DRIFT_Pass(x,Fdrift2*Lint);
	  L = Fdrift2 * Lint;
	  u = L / (1 + x[delta_]);
	  x[x_] = x[x_] + x[px_] * u;
	  x[y_] = x[y_] + x[py_] * u;
	  x[z_] =
	    x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								     x
								     [delta_]);

	  //quad_kick_pass(x, k1l_kick1, k1sl_kick1);
	  k1l = k1l_kick1;
	  k1sl = k1sl_kick1;
	  By = 0.;
	  Bx = 0.;
	  Xn = 1.;
	  Yn = 0.;
	  for (j = 1; j < Norder + 1; j++)
	    {
	      Xn0 = Xn;
	      Yn0 = Yn;
	      Xn = Xn0 * x[x_] - Yn0 * x[y_];
	      Yn = Xn0 * x[y_] + Yn0 * x[x_];
	    }
	  By = By + (k1l * Xn - k1sl * Yn);
	  Bx = Bx + (k1l * Yn + k1sl * Xn);
	  x[px_] = x[px_] - By;
	  x[py_] = x[py_] + Bx;

	  //DRIFT_Pass(x,Fdrift1*Lint);
	  L = Fdrift1 * Lint;
	  u = L / (1 + x[delta_]);
	  x[x_] = x[x_] + x[px_] * u;
	  x[y_] = x[y_] + x[py_] * u;
	  x[z_] =
	    x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								     x
								     [delta_]);
	}
    }
}

template < class T > void
SEXT_Pass (T x[6], double L, int Nint, double k2l, double k2sl)
{
  int i;
  double Lint = L / Nint;

  if (L == 0.)
    {
      sext_kick_pass (x, k2l, k2sl);
    }
  else
    {
      double k2l_kick1, k2sl_kick1;
      double k2l_kick2, k2sl_kick2;
      k2l_kick1 = Fkick1 * k2l / Nint;
      k2sl_kick1 = Fkick1 * k2sl / Nint;
      k2l_kick2 = Fkick2 * k2l / Nint;
      k2sl_kick2 = Fkick2 * k2sl / Nint;

      int j, Norder = 2;
      T Xn, Yn, Xn0, Yn0;
      T By, Bx;
      T u;
      for (i = 0; i < Nint; i++)
	{

	  //DRIFT_Pass(x,Fdrift1*Lint);
	  L = Fdrift1 * Lint;
	  u = L / (1 + x[delta_]);
	  x[x_] = x[x_] + x[px_] * u;
	  x[y_] = x[y_] + x[py_] * u;
	  x[z_] =
	    x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								     x
								     [delta_]);

	  //sext_kick_pass(x, k2l_kick1, k2sl_kick1);
	  k2l = k2l_kick1;
	  k2sl = k2sl_kick1;
	  By = 0.;
	  Bx = 0.;
	  Xn = 1.;
	  Yn = 0.;
	  for (j = 1; j < Norder + 1; j++)
	    {
	      Xn0 = Xn;
	      Yn0 = Yn;
	      Xn = Xn0 * x[x_] - Yn0 * x[y_];
	      Yn = Xn0 * x[y_] + Yn0 * x[x_];
	    }
	  By = By + (k2l * Xn - k2sl * Yn) / 2;
	  Bx = Bx + (k2l * Yn + k2sl * Xn) / 2;
	  x[px_] = x[px_] - By;
	  x[py_] = x[py_] + Bx;

	  //DRIFT_Pass(x,Fdrift2*Lint);
	  L = Fdrift2 * Lint;
	  u = L / (1 + x[delta_]);
	  x[x_] = x[x_] + x[px_] * u;
	  x[y_] = x[y_] + x[py_] * u;
	  x[z_] =
	    x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								     x
								     [delta_]);

	  //sext_kick_pass(x, k2l_kick2, k2sl_kick2);
	  k2l = k2l_kick2;
	  k2sl = k2sl_kick2;
	  By = 0.;
	  Bx = 0.;
	  Xn = 1.;
	  Yn = 0.;
	  for (j = 1; j < Norder + 1; j++)
	    {
	      Xn0 = Xn;
	      Yn0 = Yn;
	      Xn = Xn0 * x[x_] - Yn0 * x[y_];
	      Yn = Xn0 * x[y_] + Yn0 * x[x_];
	    }
	  By = By + (k2l * Xn - k2sl * Yn) / 2;
	  Bx = Bx + (k2l * Yn + k2sl * Xn) / 2;
	  x[px_] = x[px_] - By;
	  x[py_] = x[py_] + Bx;

	  //DRIFT_Pass(x,Fdrift2*Lint);
	  L = Fdrift2 * Lint;
	  u = L / (1 + x[delta_]);
	  x[x_] = x[x_] + x[px_] * u;
	  x[y_] = x[y_] + x[py_] * u;
	  x[z_] =
	    x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								     x
								     [delta_]);

	  //sext_kick_pass(x, k2l_kick1, k2sl_kick1);
	  k2l = k2l_kick1;
	  k2sl = k2sl_kick1;
	  By = 0.;
	  Bx = 0.;
	  Xn = 1.;
	  Yn = 0.;
	  for (j = 1; j < Norder + 1; j++)
	    {
	      Xn0 = Xn;
	      Yn0 = Yn;
	      Xn = Xn0 * x[x_] - Yn0 * x[y_];
	      Yn = Xn0 * x[y_] + Yn0 * x[x_];
	    }
	  By = By + (k2l * Xn - k2sl * Yn) / 2;
	  Bx = Bx + (k2l * Yn + k2sl * Xn) / 2;
	  x[px_] = x[px_] - By;
	  x[py_] = x[py_] + Bx;

	  //DRIFT_Pass(x,Fdrift1*Lint);
	  L = Fdrift1 * Lint;
	  u = L / (1 + x[delta_]);
	  x[x_] = x[x_] + x[px_] * u;
	  x[y_] = x[y_] + x[py_] * u;
	  x[z_] =
	    x[z_] - (x[px_] * x[px_] + x[py_] * x[py_]) * u / 2.0 / (1 +
								     x
								     [delta_]);
	}
    }
}

template < class T > void
OCT_Pass (T x[6], double L, int Nint, double k3l, double k3sl)
{
  int i;
  double Lint = L / Nint;

  if (L == 0.)
    {
      oct_kick_pass (x, k3l, k3sl);
    }
  else
    {
      double k3l_kick1, k3sl_kick1;
      double k3l_kick2, k3sl_kick2;
      k3l_kick1 = Fkick1 * k3l / Nint;
      k3sl_kick1 = Fkick1 * k3sl / Nint;
      k3l_kick2 = Fkick2 * k3l / Nint;
      k3sl_kick2 = Fkick2 * k3sl / Nint;
      for (i = 0; i < Nint; i++)
	{
	  DRIFT_Pass (x, Fdrift1 * Lint);
	  oct_kick_pass (x, k3l_kick1, k3sl_kick1);
	  DRIFT_Pass (x, Fdrift2 * Lint);
	  oct_kick_pass (x, k3l_kick2, k3sl_kick2);
	  DRIFT_Pass (x, Fdrift2 * Lint);
	  oct_kick_pass (x, k3l_kick1, k3sl_kick1);
	  DRIFT_Pass (x, Fdrift1 * Lint);
	}
    }
}

template < class T > void
MULT_Pass (T x[6], double L, int Nint, int Norder, double KNL[11],
	   double KNSL[11])
{
  int i;
  double Lint = L / Nint;

  if (L == 0.)
    {
      mult_kick_pass (x, Norder, KNL, KNSL);
    }
  else
    {
      double knl_kick1[11], knsl_kick1[11];
      double knl_kick2[11], knsl_kick2[11];
      for (i = 0; i < 11; i++)
	{
	  knl_kick1[i] = Fkick1 * KNL[i] / Nint;
	  knsl_kick1[i] = Fkick1 * KNSL[i] / Nint;
	}
      for (i = 0; i < 11; i++)
	{
	  knl_kick2[i] = Fkick2 * KNL[i] / Nint;
	  knsl_kick2[i] = Fkick2 * KNSL[i] / Nint;
	}

      for (i = 0; i < Nint; i++)
	{
	  DRIFT_Pass (x, Fdrift1 * Lint);
	  mult_kick_pass (x, Norder, knl_kick1, knsl_kick1);
	  DRIFT_Pass (x, Fdrift2 * Lint);
	  mult_kick_pass (x, Norder, knl_kick2, knsl_kick2);
	  DRIFT_Pass (x, Fdrift2 * Lint);
	  mult_kick_pass (x, Norder, knl_kick1, knsl_kick1);
	  DRIFT_Pass (x, Fdrift1 * Lint);
	}
    }
}

template < class T > void
GMULT_Pass (T x[6], double L, int Nint, double Angle, int Norder,
	    double KNL[11], double KNSL[11], double E1, double E2)
{
  if (L == 0.)
    {
      mult_kick_pass (x, Norder, KNL, KNSL);
    }
  else
    {
      int i;
      double href = Angle / L;
      double Lint = L / Nint;
      double knl_kick1[11], knsl_kick1[11];
      double knl_kick2[11], knsl_kick2[11];

      for (i = 0; i < 11; i++)
	{
	  knl_kick1[i] = Fkick1 * KNL[i] / Nint;
	  knsl_kick1[i] = Fkick1 * KNSL[i] / Nint;
	}
      for (i = 0; i < 11; i++)
	{
	  knl_kick2[i] = Fkick2 * KNL[i] / Nint;
	  knsl_kick2[i] = Fkick2 * KNSL[i] / Nint;
	}

      x[1] = x[1] + tan (E1) * x[0] * href;
      x[3] = x[3] - tan (E1) * x[2] * href;
      for (i = 0; i < Nint; i++)
	{
	  DRIFT_Pass (x, Fdrift1 * Lint);
	  bend_mult_kick_pass (x, Fkick1 * Lint, href, Norder, knl_kick1,
			       knsl_kick1);
	  DRIFT_Pass (x, Fdrift2 * Lint);
	  bend_mult_kick_pass (x, Fkick2 * Lint, href, Norder, knl_kick2,
			       knsl_kick2);
	  DRIFT_Pass (x, Fdrift2 * Lint);
	  bend_mult_kick_pass (x, Fkick1 * Lint, href, Norder, knl_kick1,
			       knsl_kick1);
	  DRIFT_Pass (x, Fdrift1 * Lint);
	}
      x[1] = x[1] + tan (E2) * x[0] * href;
      x[3] = x[3] - tan (E2) * x[2] * href;
    }
}

//The explicit instantiation part of template, add when needed
template void DRIFT_Pass<tps>(tps x[6],double);
template void DRIFT_Pass<double>(double x[6],double);
template void GtoL<tps>(tps x[6],double,double,double);
template void GtoL<double>(double x[6],double,double,double);
template void LtoG<tps>(tps x[6],double,double,double);
template void LtoG<double>(double x[6],double,double,double);
template void SEXT_Pass<tps>(tps x[6], double, int, double, double);
template void SEXT_Pass<double>(double x[6], double, int, double, double);
