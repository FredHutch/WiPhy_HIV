#ifndef PLOTPOINTS_H
#define PLOTPOINTS_H

#define MAX_POINTS 1000
#include "colors.h"

class plotPoints {

public:
	double *time;

	long int *vt;
	long int *vt_vi;
	long int *vt_non;
	long int **vt_s;
	int *act;
	int *lat;
	int *cd8;
	int *cd8_junk;
	int *junk_L0;
	int *junk_L1;
	int *junk_L2;
	int **cd8_s;
	int **dh_s;
	int **id_s;
	int **a_s;

	double *R0s;
	int **color;

   	int valid;
   	int max_points;

	plotPoints(void);
	~plotPoints(void);

	void checkForRealloc();
};
#endif
