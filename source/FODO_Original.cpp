#include "FODO4.h"

GlobalVariables GP;
//add SPKICK every 0.5 meters
void line_def(Line & FODO)
{
    const double ld=1.0;
    const double lq=0.5;
    const double kf=1.369, kd=1.427;

    Element * temp;
    temp = new MARKER("START",0);
    FODO.Append(temp);
    for(int i=1;i<18;i++){
        temp = new  MARKER("SPKICK",0);
        FODO.Append(temp);
        temp = new DRIFT("D1", 0.5*ld);
        FODO.Append(temp);
        temp = new  MARKER("SPKICK",0);
        FODO.Append(temp);
        temp = new  QUAD("QF", lq, kf*lq);
        FODO.Append(temp);
        temp = new  MARKER("SPKICK",0);
        FODO.Append(temp);
        temp = new DRIFT("D1", 0.5*ld);
        FODO.Append(temp);
        temp = new  MARKER("SPKICK",0);
        FODO.Append(temp);
        temp = new DRIFT("D1", 0.5*ld);
        FODO.Append(temp);
        temp = new  MARKER("SPKICK",0);
        FODO.Append(temp);
        temp = new  QUAD("QD", lq, -kd*lq);
        FODO.Append(temp);
        temp = new  MARKER("SPKICK",0);
        FODO.Append(temp);
        temp = new DRIFT("D1", 0.5*ld);
        FODO.Append(temp);
    }
    temp = new MARKER("END",0);
    FODO.Append(temp);
}
