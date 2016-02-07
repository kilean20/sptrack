#include "tpsa.h"
//==============================================================
//
//         linear  tpsa: x[7]  and  linear map x[6][7]
//
//===============================================================
tps::tps ()
{
  int i;
  for (i = 0; i < 7; i++)
    sample[i] = 0.;
}
tps::tps (double d)
{
  int i;
  for (i = 0; i < 7; i++)
    sample[i] = 0;
  sample[0] = d;
}
tps::tps (double d[7])
{
  int i;
  for (i = 0; i < 7; i++)
    sample[i] = d[i];
}

double &
tps::operator[] (int i)
{
  return sample[i];
}

ostream & operator<< (ostream & stream, tps x)
{
  int
    i;
  for (i = 0; i < 7; i++)
    stream << setw (12) << setprecision (9) << scientific << x[i] << "  ";
  return stream;
}

tps
operator+ (tps x, tps y)
{
  int
    i;
  tps
    z;

  for (i = 0; i < 7; i++)
    z[i] = x[i] + y[i];
  return z;
}

tps
operator- (tps x, tps y)
{
  int
    i;
  tps
    z;

  for (i = 0; i < 7; i++)
    z[i] = x[i] - y[i];
  return z;
}

tps
operator* (tps x, tps y)
{
  int
    i;
  tps
    z;

  z[0] = x[0] * y[0];
  for (i = 1; i < 7; i++)
    z[i] = x[0] * y[i] + x[i] * y[0];
  return z;
}

tps
DAinv (tps x)
{
  int
    i;
  double
    a,
    temp;
  tps
    z;

  z[0] = 1.0 / x[0];
  temp = x[0];
  a = -1.0 / (temp * temp);
  for (i = 1; i < 7; i++)
    z[i] = a * x[i];
  return z;
}

tps
operator/ (tps x, tps y)
{
  return x * DAinv (y);
}

tps
sqr (tps x)
{
  return x * x;
}

tps
sqrt (tps x)
{
  int
    i;
  double
    a;
  tps
    z;

  a = sqrt (x[0]);
  z[0] = a;
  a = 0.5 / a;
  for (i = 1; i < 7; i++)
    z[i] = a * x[i];
  return z;
}

tps
sin (tps x)
{
  int
    i;
  double
    a;
  tps
    z;

  z[0] = sin (x[0]);
  a = cos (x[0]);
  for (i = 1; i < 7; i++)
    z[i] = a * x[i];
  return z;
}

tps
cos (tps x)
{
  int
    i;
  double
    a;
  tps
    z;

  z[0] = cos (x[0]);
  a = -sin (x[0]);
  for (i = 1; i < 7; i++)
    z[i] = a * x[i];
  return z;
}

tps
tan (tps x)
{
  return sin (x) / cos (x);
}

tps
atan (tps x)
{
  int
    i;
  double
    a;
  tps
    z;

  a = x[0];
  z[0] = atan (a);
  a = 1 / (1 + a * a);
  for (i = 1; i < 7; i++)
    z[i] = a * x[i];
  return z;
}

tps
sinh (tps x)
{
  int
    i;
  double
    a;
  tps
    z;

  z[0] = sinh (x[0]);
  a = cosh (x[0]);
  for (i = 1; i < 7; i++)
    z[i] = a * x[i];
  return z;
}

tps
cosh (tps x)
{
  int
    i;
  double
    a;
  tps
    z;

  z[0] = cosh (x[0]);
  a = sinh (x[0]);
  for (i = 1; i < 7; i++)
    z[i] = a * x[i];
  return z;
}

tps
exp (tps x)
{
  int
    i;
  double
    a;
  tps
    z;

  a = exp (x[0]);
  z[0] = a;
  for (i = 1; i < 7; i++)
    z[i] = a * x[i];
  return z;
}

tps
ln (tps x)
{
  int
    i;
  tps
    z;

  z[0] = log (x[0]);
  for (i = 1; i < 7; i++)
    z[i] = x[i] / x[0];
  return z;
}






////////
linmap::linmap ()
{
  int i, j;
  for (i = 0; i < 6; i++)
    for (j = 0; j < 7; j++)
      map0[i][j] = 0.;
}
linmap::linmap (double x[6])
{
  int i, j;
  for (i = 0; i < 6; i++)
    for (j = 0; j < 7; j++)
      map0[i][j] = 0.;
  for (i = 0; i < 6; i++)
    map0[i][0] = x[i];
}

void
linmap::identity ()
{
  int i, j;
  for (i = 0; i < 6; i++)
    for (j = 0; j < 7; j++)
      map0[i][j] = 0.;
  map0[0][1] = 1;
  map0[1][2] = 1;
  map0[2][3] = 1;
  map0[3][4] = 1;
  map0[4][5] = 1;
  map0[5][6] = 1;
}

tps & linmap::operator[](int i)
{
  return map0[i];
}

linmap
operator+ (linmap x, linmap y)
{
  int
    i,
    j;
  linmap
    z;
  for (i = 0; i < 6; i++)
    for (j = 0; j < 7; j++)
      z[i][j] = x[i][j] + y[i][j];
  return z;
}

void
linmap::print ()
{
  int i;
  for (i = 0; i < 6; i++)
    cout << map0[i] << endl;
}

void
Getmat (linmap map0, double x[36])
{
  int i, j;
  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
      x[i * 6 + j] = map0[i][j + 1];
}

void
Getpos (linmap map0, double x[6])
{
  int i;
  for (i = 0; i < 6; i++)
    x[i] = map0[i][0];
}
