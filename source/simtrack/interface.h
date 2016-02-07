#ifndef INTERFACE_H
#define INTERFACE_H
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <cmath>
#include "elepass.h"
#include "global.h"
extern GlobalVariables GP;
void  Read_MADXLattice(const char * filename, Line & linename);
void  Print_MADXLattice(const char * filename, Line & linename);

#endif
