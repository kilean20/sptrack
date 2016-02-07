#ifndef TRACKING_H
#define TRACKING_H
#include <fstream>
#include <vector>
#include <cmath>
#include "elepass.h"
#include "smath.h"
#include "global.h"
extern GlobalVariables GP;
//==============================================================
//
//                Tracking 
//
//===============================================================

void Track(Line & linename, double x[6], int nturn, int & stable, int & lost_turn, int & lost_post);
void Track_tbt(Line & linename, double x[6], int nturn, double  x_tbt[], int & stable, int & lost_turn, int & lost_post);
void Track_tbt(Line & linename, double x[6], int nturn, double  x_tbt[], int &bpm_index, int & stable, int & lost_turn, int & lost_post);

//----------------------------------------------------------------
//Tracking without  linename.Cell[i]->Pass() to speed up tracking
//----------------------------------------------------------------
void Prepare_Track_Fast(Line & linename);
//  this function only accepts the allowed elements and change others to drift
//  this function accepts thick lenses: SBEND, QUAD, SEXT, MULT, ELENS, GMULT
void Track_Fast(Line & linename, double x[6], int nturn, int & stable, int & lost_turn, int & lost_post, double & sum_x2, double & sum_y2, double & sum_z2);
void Track_DA( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0 );
// Track with ->Pass()
void Track_Fast_DA( Line & linename, int nturn, double deltap0, double sigmax0, double sigmay0 );
// Track without class

#endif
