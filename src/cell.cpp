#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include <cstdlib>
#include<cmath>
using namespace std;

#include <strings.h>

#include "cell.h"

cell::cell(double birthday, int orig, int strain_num, int state)
{
	cellState = state;
	cellStrain = strain_num;
	cellBirth = birthday;
	cellOrigin = orig;
}

cell::~cell(void) 
{
}
void cell::set_state(int state)
{
    cellState = state;
}

int cell::get_state()
{
    return cellState;
}
int cell::get_origin()
{
    return cellOrigin;
}
int cell::get_strain()
{
    return cellStrain;
}
double cell::get_birthdate()
{
    return cellBirth;
}
