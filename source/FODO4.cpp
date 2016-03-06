#include "FODO4.h"

GlobalVariables GP;
//add SPKICK every 0.5 meters
void line_def(Line & FODO)
{
    const double ld=0.6;
    const double lq=0.4;

    Element * temp;
    temp = new MARKER("START",0);
    FODO.Append(temp);
    for(int i=0;i<5;i++){
        temp = new DRIFT("D1", 0.5*ld);
        FODO.Append(temp);
        temp = new  QUAD("QF", 0.5*lq, 1.00428);
        FODO.Append(temp);
        temp = new  MARKER("SPKICK",0);
        FODO.Append(temp);
        temp = new  QUAD("QF", 0.5*lq, 1.00428);
        FODO.Append(temp);
        temp = new DRIFT("D0", ld);
        FODO.Append(temp);
        temp = new  QUAD("QD", 0.5*lq, -0.837118);
        FODO.Append(temp);
        temp = new  MARKER("SPKICK",0);
        FODO.Append(temp);
        temp = new  QUAD("QD", 0.5*lq, -0.837118);
        FODO.Append(temp);
        temp = new DRIFT("D1", 0.5*ld);
        FODO.Append(temp);
    }
    temp = new MARKER("END",0);
    FODO.Append(temp);
}
