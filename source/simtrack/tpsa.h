#ifndef TPSA_H
#define TPSA_H
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;
//==============================================================
//
//         linear  tpsa: x[7]  and  linear map x[6][7]
//
//===============================================================
class tps
{
public:
  tps ();
  tps (double d);
    tps (double d[7]);
  double &operator[] (int i);
  friend ostream & operator<< (ostream & stream, tps x);
  friend tps operator+ (tps x, tps y);
  friend tps operator- (tps x, tps y);
  friend tps operator* (tps x, tps y);
  friend tps DAinv (tps x);
  friend tps operator/ (tps x, tps y);
  friend tps sqr (tps x);
  friend tps sqrt (tps x);
  friend tps sin (tps x);
  friend tps cos (tps x);
  friend tps tan (tps x);
  friend tps atan (tps x);
  friend tps sinh (tps x);
  friend tps cosh (tps x);
  friend tps exp (tps x);
  friend tps ln (tps x);
private:
  double sample[7];
};

class linmap
{
public:
  linmap ();
  linmap (double x[6]);
  void identity ();
    tps & operator[] (int i);
  friend linmap operator+ (linmap x, linmap y);
  void print ();
private:
    tps map0[6];
};

void Getmat (linmap map0, double x[36]);
void Getpos (linmap map0, double x[6]);

#endif
