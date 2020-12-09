#include <cstdlib>
using namespace std;

#include <strings.h>
#include "colors.h"
#include "plotpoints.h"

plotPoints::plotPoints(void) 
{
	valid = 0;
	max_points=MAX_POINTS;
	time = (double*) malloc(MAX_POINTS * sizeof(double));
	vt = (long int *)malloc(MAX_POINTS * sizeof (long int));
	vt_vi = (long int *)malloc(MAX_POINTS * sizeof (long int));
	vt_non = (long int *)malloc(MAX_POINTS * sizeof (long int));
	vt_s = (long int **)malloc(MAX_POINTS * sizeof (long int*));
	act = (int *)malloc(MAX_POINTS * sizeof (int));
	lat = (int *)malloc(MAX_POINTS * sizeof (int));
	junk_L0 = (int *)malloc(MAX_POINTS * sizeof (int));
	junk_L1 = (int *)malloc(MAX_POINTS * sizeof (int));
	junk_L2 = (int *)malloc(MAX_POINTS * sizeof (int));
	a_s = (int **)malloc(MAX_POINTS * sizeof (int*));
	cd8 = (int *)malloc(MAX_POINTS * sizeof (int));
	cd8_junk = (int *)malloc(MAX_POINTS * sizeof (int));
	cd8_s = (int **)malloc(MAX_POINTS * sizeof (int*));
	dh_s = (int **)malloc(MAX_POINTS * sizeof (int*));
	id_s = (int **)malloc(MAX_POINTS * sizeof (int*));
	color = (int **)malloc(MAX_POINTS * sizeof (int *));

	for (int i=0; i < MAX_POINTS; i++)
	{
	    id_s[i]=(int *)malloc(MAX_FOLLOW_STRAINS * sizeof (int));
	    vt_s[i]=(long int *)malloc(MAX_FOLLOW_STRAINS * sizeof (long int));
	    cd8_s[i]=(int *)malloc(MAX_FOLLOW_STRAINS * sizeof (int));
	    dh_s[i]=(int *)malloc(MAX_FOLLOW_STRAINS * sizeof (int));
	    color[i]=(int *)malloc(MAX_FOLLOW_STRAINS * sizeof (int));
	    a_s[i]=(int *)malloc(MAX_FOLLOW_STRAINS * sizeof (int));
	}

	R0s = (double *)malloc(MAX_POINTS * sizeof (double));
};

plotPoints::~plotPoints(void) 
{
	free(vt);
	free(vt_vi);
	free(vt_non);
	free(cd8);
	free(cd8_junk);
	free(time);
	free(R0s);
	free(act);
	free(lat);
	free(junk_L0);
	free(junk_L1);
	free(junk_L2);

	for (int i=0; i < max_points; i++)
	{
	    free(id_s[i]);
	    free(vt_s[i]);
	    free(cd8_s[i]);
	    free(dh_s[i]);
	    free(a_s[i]);
	    free(color[i]);
	}
	free(id_s);
	free(vt_s);
	free(cd8_s);
	free(dh_s);
	free(a_s);
	free(color);
}

void plotPoints::checkForRealloc()
{
	/* check for need to realloc point array */
	if (valid >= max_points-1) {
	    double *tempTimes = time;
	    time = (double*) malloc(max_points*2 * sizeof(double));
	    bcopy((char *)tempTimes,(char *)time,
		    max_points*sizeof(double));

	    free(tempTimes);

	    long int *tempVL = vt;
	    vt = (long int *)malloc(max_points*2 * sizeof (long int));
	    bcopy((char *)tempVL,(char *)vt,
		    max_points*sizeof(long int));

	    free(tempVL);

	    tempVL = vt_vi;
	    vt_vi = (long int *)malloc(max_points*2 * sizeof (long int));
	    bcopy((char *)tempVL,(char *)vt_vi,
		    max_points*sizeof(long int));

	    free(tempVL);

	    tempVL = vt_non;
	    vt_non = (long int *)malloc(max_points*2 * sizeof (long int));
	    bcopy((char *)tempVL,(char *)vt_non,
		    max_points*sizeof(long int));

	    free(tempVL);

	    int **tempID = id_s;
	    id_s = (int **)malloc(max_points*2 * sizeof (int*));
	    bcopy((char *)tempID,(char *)id_s,
		    max_points*sizeof(int*));

	    free(tempID);

	    int **tempA = a_s;
	    a_s = (int **)malloc(max_points*2 * sizeof (int*));
	    bcopy((char *)tempA,(char *)a_s,
		    max_points*sizeof(int*));

	    free(tempA);

	    long int **tempVe = vt_s;
	    vt_s = (long int **)malloc(max_points*2 * sizeof (long int*));
	    bcopy((char *)tempVe,(char *)vt_s,
		    max_points*sizeof(long int*));

	    free(tempVe);

	    int *tempCD8 = cd8;
	    cd8 = (int *)malloc(max_points*2 * sizeof (int));
	    bcopy((char *)tempCD8,(char *)cd8,
		    max_points*sizeof(int));

	    free(tempCD8);

	    int *tempJunkCD8 = cd8_junk;
	    cd8_junk = (int *)malloc(max_points*2 * sizeof (int));
	    bcopy((char *)tempJunkCD8,(char *)cd8_junk,
		    max_points*sizeof(int));

	    free(tempJunkCD8);

	    int **tempCD8_s = cd8_s;
	    cd8_s = (int **)malloc(max_points*2 * sizeof (int*));
	    bcopy((char *)tempCD8_s,(char *)cd8_s,
		    max_points*sizeof(int*));

	    free(tempCD8_s);

	    int **tempDH_s = dh_s;
	    dh_s = (int **)malloc(max_points*2 * sizeof (int*));
	    bcopy((char *)tempDH_s,(char *)dh_s,
		    max_points*sizeof(int*));

	    free(tempDH_s);

	    for (int i=0; i < max_points; i++)
	    {
		vt_s[max_points+i]=(long int *)malloc(MAX_FOLLOW_STRAINS * sizeof (long int));
		id_s[max_points+i]=(int *)malloc(MAX_FOLLOW_STRAINS * sizeof (int));
		cd8_s[max_points+i]=(int *)malloc(MAX_FOLLOW_STRAINS * sizeof (long int));
		dh_s[max_points+i]=(int *)malloc(MAX_FOLLOW_STRAINS * sizeof (long int));
		a_s[max_points+i]=(int *)malloc(MAX_FOLLOW_STRAINS * sizeof (int));
	    }

	    int **tempColor = color;
	    color = (int **)malloc(max_points*2 * sizeof (int*));
	    bcopy((char *)tempColor,(char *)color,
		    max_points*sizeof(int*));

	    free(tempColor);

	    for (int i=0; i < max_points; i++)
		color[max_points+i]=(int *)malloc(MAX_FOLLOW_STRAINS * sizeof (int));

	    double *tempR0s = R0s;
	    R0s = (double *)malloc(max_points*2 * sizeof (double));
	    bcopy((char *)tempR0s,(char *)R0s,
		    max_points*sizeof(double));

	    free(tempR0s);

	    int *tempAct = act;
	    act = (int *)malloc(max_points*2 * sizeof (int));
	    bcopy((char *)tempAct,(char *)act,
		    max_points*sizeof(int));

	    free(tempAct);

	    int *tempLat = lat;
	    lat = (int *)malloc(max_points*2 * sizeof (int));
	    bcopy((char *)tempLat,(char *)lat,
		    max_points*sizeof(int));

	    free(tempLat);

	    int *tempJunk = junk_L0;
	    junk_L0 = (int *)malloc(max_points*2 * sizeof (int));
	    bcopy((char *)tempJunk,(char *)junk_L0,
		    max_points*sizeof(int));

	    free(tempJunk);

	    tempJunk = junk_L1;
	    junk_L1 = (int *)malloc(max_points*2 * sizeof (int));
	    bcopy((char *)tempJunk,(char *)junk_L1,
		    max_points*sizeof(int));

	    free(tempJunk);

	    tempJunk = junk_L2;
	    junk_L2 = (int *)malloc(max_points*2 * sizeof (int));
	    bcopy((char *)tempJunk,(char *)junk_L2,
		    max_points*sizeof(int));

	    free(tempJunk);

	    max_points *= 2;
	}
}
