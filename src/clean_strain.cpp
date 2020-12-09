//NOTE: This read-only file is generated from the file ../src/full_strain.cpp.
//Any editing should be done in that file.

#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include <cstdlib>
#include<cmath>
using namespace std;

#include <strings.h>

#include "settings.h"
#include "strain.h"

#define MUTATION_PRINTS 1000
#define MIN_STRAIN_CD8S 1
#define NUM_BETA_SAMPLES 100
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

int strain::cd8Groups = 0;

int strain::get_cd8_groups()
{
    return cd8Groups;
}
static int mut_prints=0;

char translate(char b1, char b2, char b3)
{
    switch (b1) {
	case 'T':
	    switch (b2) {
		case 'T':
		    switch (b3) {
			case 'T':
			case 'C':
			    return 'F';
			default:
			    return 'L';
		    }
		case 'C':
		    return 'S';
		case 'A':
		    switch (b3) {
			case 'T':
			case 'C':
			    return 'Y';
			default:
			    return '-'; // stop codon
		    }
		case 'G':
		    switch (b3) {
			case 'T':
			case 'C':
			    return 'C';
			case 'G':
			    return 'W';
			default:
			    return '-'; // stop codon
		    }
	    }
	case 'C':
	    switch (b2) {
		case 'T':
		    return 'L';
		case 'C':
		    return 'P';
		case 'A':
		    switch (b3) {
			case 'T':
			case 'C':
			    return 'H';
			default:
			    return 'Q';
		    }
		case 'G':
		    return 'R';
	    }
	case 'A':
	    switch (b2) {
		case 'T':
		    switch (b3) {
			case 'T':
			case 'C':
			case 'A':
			    return 'I';
			default:
			    return 'M';
		    }
		case 'C':
		    return 'T';
		case 'A':
		    switch (b3) {
			case 'T':
			case 'C':
			    return 'N';
			default:
			    return 'K';
		    }
		case 'G':
		    switch (b3) {
			case 'T':
			case 'C':
			    return 'S';
			default:
			    return 'R';
		    }
	    }
	case 'G':
	    switch (b2) {
		case 'T':
		    return 'V';
		case 'C':
		    return 'A';
		case 'A':
		    switch (b3) {
			case 'T':
			case 'C':
			    return 'D';
			default:
			    return 'E';
		    }
		case 'G':
		    return 'G';
	    }
	}
	return '-';
}
double strain::get_delta_hamming(settings *params)
{
    int delta_hamming=1;
    double extra_mut = gsl_rng_uniform(params->ur);

    if (extra_mut<params->mu*params->mu)
	delta_hamming=3;
    else if (extra_mut<params->mu)
	delta_hamming=2;
    return delta_hamming;
}

strain::strain(settings *params, int num, double birthday, 
		strain *mother_strain, double inputBeta)
{
	strainNum = num;
	parentStrain = mother_strain;

	// NEW: recombination fields 1/8/2019
	alt_parentStrain = NULL;
	//start_strain = NULL;
	//fragment_len=0;
	//fragment_start=0;
	non_syn_ptr=NULL;

	if (mother_strain==NULL)
	{
	    strainBeta = params->Bt_max;
	    if (params->immun_model == 6)
		dImmun_s = params->dImmun;

	    //cd8Cells = params->cd8_mean * params->vol / 10;
	    //if (cd8Cells == 0)
	    //cd8Cells = gsl_ran_weibull(params->ur,params->cd8_0,params->cd8_k);
	    for (int i=0; i < MAX_COMPARTMENTS; i++)
		cd8Cells[i] = params->cd8_0;

	    founderDist = 0;
	    set_num_bases(params->founder_stats.length);
	    set_num_as(params->founder_stats.as);
	    set_num_cs(params->founder_stats.cs);
	    set_num_gs(params->founder_stats.gs);
	    set_num_ts(params->founder_stats.ts);
	}
	else
	{
	    strainBeta = inputBeta;
	    if (params->immun_model == 6)
		dImmun_s = mother_strain->dImmun_s;

	    set_num_bases(mother_strain->get_num_bases());
	    set_num_as(mother_strain->get_num_as());
	    set_num_cs(mother_strain->get_num_cs());
	    set_num_gs(mother_strain->get_num_gs());
	    set_num_ts(mother_strain->get_num_ts());

	    non_syn_ptr = mother_strain->get_non_syn_ptr();
	
	    if (params->cd8_k==0)
	    {
		for (int i=0; i < MAX_COMPARTMENTS; i++)
		    cd8Cells[i] = params->cd8_mean;
	    }
	    else
		for (int i=0; i < MAX_COMPARTMENTS; i++)
		    cd8Cells[i] = gsl_ran_weibull(params->ur,params->cd8_0,params->cd8_k);
	    int deltaHD = strain::get_delta_hamming(params);
	    double transit_prob = params->transit_prob;
	    double transv_prob = (1 - params->transit_prob)/2.0;
	    for (int i=0; i < deltaHD; i++) {
		int spot=-1;
		if (params->use_entropy)
		{
		    spot = gsl_rng_uniform_int(params->ur,params->founder_stats.length);
		}

		double ran_mutation= gsl_rng_uniform(params->ur);
		double alt_base= gsl_rng_uniform(params->ur);
		int num_bases=mother_strain->get_num_bases();
		int num_as=mother_strain->get_num_as();
		int num_cs=mother_strain->get_num_cs();
		int num_gs=mother_strain->get_num_gs();
		int num_ts=mother_strain->get_num_ts();
		double fract_as=(double)num_as/(double)num_bases;
		double fract_cs=(double)num_cs/(double)num_bases;
		double fract_gs=(double)num_gs/(double)num_bases;
		double fract_ts=(double)num_ts/(double)num_bases;

		// transitions and transversions from A (w/ transitions twice as likely)
		// HIV spec =A-> C(9e-7), T(7e-7), G(6e-6)
		char old_prot = '-';
		char new_prot = '-';
		char new_base = '-';
		char old_base = '-';
		int mut_base = 0;
		char c1,c2,c3;
		if (spot >= 0) {
		    if (spot % 3 == 0) {
			c1=params->bases[spot]; c2=params->bases[spot+1]; c3=params->bases[spot+2];mut_base = 0;
		    } else if (spot % 3 == 1) {
			c1=params->bases[spot-1]; c2=params->bases[spot]; c3=params->bases[spot+1];mut_base = 1;
		    } else {
			c1=params->bases[spot-2]; c2=params->bases[spot-1]; c3=params->bases[spot];mut_base = 2;
		    }
		    old_prot = translate(c1,c2,c3);
		    old_base = params->bases[spot];
		}
		if ((spot >= 0 && params->bases[spot]=='A') || (spot < 0 && fract_as > ran_mutation))
		{
		    set_num_as(get_num_as()-1);
		    if (params->use_k2p)
		    {
			if (alt_base < 0.333) {
			    set_num_cs(get_num_cs()+1);
			    new_base='C';
			} else if (alt_base < 0.666) {
			    set_num_gs(get_num_gs()+1);
			    new_base='G';
			} else {
			    set_num_ts(get_num_ts()+1);
			    new_base='T';
			}
		    }
		    else
		    {
			if (alt_base < transv_prob /*HIV spec =0.118*/)  {
			    set_num_cs(get_num_cs()+1);
			    new_base='C';
			} else if (alt_base < transv_prob+transit_prob /*HIV spec =0.908*/)  {
			    set_num_gs(get_num_gs()+1);
			    new_base='G';
			} else {
			    set_num_ts(get_num_ts()+1);
			    new_base='T';
			}
		    }
		}
		// transitions and transversions from C (w/ transitions twice as likely)
		// HIV spec =C-> A(5e-6), T(1.2e-5), G(5e-7)
		else if ((spot >= 0 && params->bases[spot]=='C') ||(spot < 0 && fract_as+fract_cs > ran_mutation))
		{
		    set_num_cs(get_num_cs()-1);
		    if (params->use_k2p)
		    {
			if (alt_base < 0.333) {
			    set_num_as(get_num_as()+1);
			    new_base='A';
			} else if (alt_base < 0.666) {
			    set_num_gs(get_num_gs()+1);
			    new_base='G';
			} else {
			    set_num_ts(get_num_ts()+1);
			    new_base='T';
			}
		    }
		    else
		    {
			if (alt_base < transv_prob /*HIV spec =10.0/71*/)  {
			    set_num_as(get_num_as()+1);
			    new_base='A';
			} else if (alt_base < transit_prob /*HIV spec =11.0/71*/)  {
			    set_num_gs(get_num_gs()+1);
			    new_base='G';
			} else {
			    set_num_ts(get_num_ts()+1);
			    new_base='T';
			}
		    }
		}
		// transitions and transversions from G (w/ transitions twice as likely)
		// HIV spec =G-> A(1.6e-5), T(2e-6), C(1e-7)
		else if ((spot >= 0 && params->bases[spot]=='G') ||(spot < 0 && fract_as+fract_cs+fract_gs > ran_mutation))
		{
		    set_num_gs(get_num_gs()-1);
		    if (params->use_k2p)
		    {
			if (alt_base < 0.333) {
			    set_num_as(get_num_as()+1);
			    new_base='A';
			} else if (alt_base < 0.666) {
			    set_num_cs(get_num_cs()+1);
			    new_base='C';
			} else {
			    set_num_ts(get_num_ts()+1);
			    new_base='T';
			}
		    }
		    else
		    {
			if (alt_base < transit_prob /*HIV spec =160.0/181*/)  {
			    set_num_as(get_num_as()+1);
			    new_base='A';
			} else if (alt_base < transit_prob + transv_prob /*HIV spec =161.0/181*/)  {
			    set_num_cs(get_num_cs()+1);
			    new_base='C';
			} else {
			    set_num_ts(get_num_ts()+1);
			    new_base='T';
			}
		    }
		}
		// transitions and transversions from T (w/ transitions twice as likely)
		// HIV spec =T-> A(3e-6), G(3e-6), C(1e-5)
		else
		{
		    set_num_ts(get_num_ts()-1);
		    if (params->use_k2p)
		    {
			if (alt_base < 0.333) {
			    set_num_as(get_num_as()+1);
			    new_base='A';
			} else if (alt_base < 0.666) {
			    set_num_cs(get_num_cs()+1);
			    new_base='C';
			} else {
			    set_num_gs(get_num_gs()+1);
			    new_base='G';
			}
		    }
		    else
		    {
			if (alt_base < transv_prob /*HIV spec =1.0/35*/)  {
			    set_num_as(get_num_as()+1);
			    new_base='A';
			} else if (alt_base < transit_prob + transv_prob /*HIV spec =34.0/35*/)  {
			    set_num_cs(get_num_cs()+1);
			    new_base='C';
			} else {
			    set_num_gs(get_num_gs()+1);
			    new_base='G';
			}
		    }
		}
		if (spot >= 0) { // must be using entropy
		    switch(mut_base) {
			case 0:
			    new_prot = translate(new_base,c2,c3);
			    break;
			case 1:
			    new_prot = translate(c1,new_base,c3);
			    break;
			case 2:
			    new_prot = translate(c1,c2,new_base);
			    break;
		    }
		    if (new_prot != old_prot) // && params->entropy[spot] >= params->entropy_threshold)
		    {
			if (non_syn_ptr==NULL)
			{
			    // no NS block yet, create a new mask, set this bit & alter beta
			    non_syn_ptr = (char *)malloc(((params->founder_stats.length / 8)+1)*sizeof(char));
			    bzero(non_syn_ptr,(params->founder_stats.length / 8)+1);
			    if (mut_prints < MUTATION_PRINTS)
			    {
				fprintf (stderr, "New Non-synonymous mask at t=%3.2lf\n",params->time);
				mut_prints++;
			    }
			} else if (non_syn_ptr == mother_strain->get_non_syn_ptr()) {
			    // we need our own NS block, so copy our parent's
			    non_syn_ptr = (char *)malloc(((params->founder_stats.length / 8)+1)*sizeof(char));
			    bcopy(mother_strain->get_non_syn_ptr(),non_syn_ptr,(params->founder_stats.length / 8)+1);
			    if (mut_prints < MUTATION_PRINTS)
			    {
				fprintf (stderr, "Copied Non-synonymous mask at t=%3.2lf\n",params->time);
				mut_prints++;
			    }
			}
			// if no non-syn at this spot, set its bit & alter beta
			if ((non_syn_ptr[spot / 8] & (1 << spot % 8))  == 0) {
			    non_syn_ptr[spot / 8] |= (1 << (spot % 8));
			    if (params->immun_model == 6 && params->entropy[spot] > params->entropy_threshold)
				dImmun_s /= params->entropy_factor*(params->entropy[spot] / params->entropy_threshold);

			    //else if (params->immun_model != 6)
			    if (params->entropy[spot] > params->entropy_threshold)
				strainBeta *= params->entropy_factor*(params->entropy[spot] / params->entropy_threshold);
				
			// otherwise, clear its bit & remove beta adjustment previously applied
			} else {
			    non_syn_ptr[spot / 8] &= (0xFF ^ (1 << (spot % 8)));
			    if (params->immun_model == 6 && params->entropy[spot] > params->entropy_threshold)
				dImmun_s *= params->entropy_factor*(params->entropy[spot] / params->entropy_threshold);

			    //else if (params->immun_model != 6)
			    if (params->entropy[spot] > params->entropy_threshold)
				strainBeta /= params->entropy_factor*(params->entropy[spot] / params->entropy_threshold);
			}
		    }
		    //else if (new_prot != old_prot)
			//strainBeta = 0;

		    if (mut_prints < MUTATION_PRINTS)
		    {
			if (new_prot == old_prot)
			    fprintf (stderr, "Synonymous mutation at t=%3.2lf ( %c->%c still %c)\n",params->time,old_base,new_base,new_prot);
			else
			    fprintf (stderr, "Non-synonymous mutation at t=%3.2lf, ent=%g bt=%e ( %c->%c => %c->%c)\n",params->time,params->entropy[spot],
				strainBeta,old_base,new_base,old_prot,new_prot);
			mut_prints++;
		    }
		}
	    }
	    founderDist = mother_strain->get_founder_dist() + deltaHD;
	}

	strainBirth = birthday;
	cd8Group = cd8Groups;

	if (params->ic50_k > 0)
	{
	    strainIC50 = gsl_ran_weibull(params->ur,params->dImmun_IC50,params->ic50_k);
	    strainIC50_0 = gsl_ran_weibull(params->ur,params->dImmun_IC50_0,params->ic50_k);
	}
	else
	{
	    strainIC50 = params->dImmun_IC50;
	    strainIC50_0 = params->dImmun_IC50_0;
	}

	lastEventDate = birthday;
	expireDate = -1;
	activeDate = -1;  // reset on activation for L1 and L2, on 1st active
			  // or during reactivation (adjusted for down time)
	deactiveDate = -1;  // reset when active count goes to 0

	for (int i=0; i < MAX_COMPARTMENTS; i++)
	{
	    activeCells[i] = 0; // maybe active (wait for inc_active)
	    eclipseCells[i] = 0; // pre-productive active (wait for inc_ecl)
	    l1Cells[i] = 0; // current short term latently infected cells
	    l2Cells[i] = 0; // current long term latently infected cells

	    via_virions[i] = 0; // total currently assoc viable virus (except FDC held)
	    junk_virions[i] = 0; // total currently assoc junk virus (except FDC held)

	    held_viable[i] = 0; // total current FDC assoc viable virus
	    held_junk[i] = 0; // total current FDC assoc junk virus

	    // initial immune response

	    cd8MaxCells[i] = cd8Cells[i];
	    lumpedStrains[i]=0;
	    actDeathRate[i] = params->dA;
	}

	totalProlif = 0;
	totalInfected = 0;
	totalDeaths = 0;
	totalCells = 0;	// total ever assoc
	assocCells = 0;	// total currently assoc
}

strain::~strain(void) 
{
}
int strain::get_num_as()
{
    return num_as;
}
int strain::get_num_cs()
{
    return num_cs;
}
int strain::get_num_gs()
{
    return num_gs;
}
int strain::get_num_ts()
{
    return num_ts;
}

void strain::set_num_as(int new_num_as)
{
    num_as = new_num_as;
}
void strain::set_num_cs(int new_num_cs)
{
    num_cs = new_num_cs;
}
void strain::set_num_gs(int new_num_gs)
{
    num_gs = new_num_gs;
}
void strain::set_num_ts(int new_num_ts)
{
    num_ts = new_num_ts;
}

int strain::get_cd8_group()
{
    return cd8Group;
}
void strain::set_cd8_group(int new_cd8_group)
{
    cd8Group = new_cd8_group;
}

char *strain::get_non_syn_ptr()
{
    return non_syn_ptr;
}
int strain::get_num_bases()
{
    return num_bases;
}
void strain::set_num_bases(int new_num_bases)
{
    num_bases = new_num_bases;
}

int strain::get_total_infected()
{
    return totalInfected;
}
int strain::get_total_prolif()
{
    return totalProlif;
}
int strain::get_total_deaths()
{
    return totalDeaths;
}
int strain::get_cells()
{
    return assocCells;
}
int strain::get_total_cells()
{
    return totalCells;
}
int strain::get_lumped_strains(int compartment)
{
    return lumpedStrains[compartment];
}

double strain::get_death_rate(int compartment, settings *params, double time)
{
    if (params->lifespan_decay > 0)
	return actDeathRate[compartment] * 1 /(pow(M_E,(-(time - activeDate)/params->lifespan_decay)));

    return actDeathRate[compartment];
}

void strain::set_death_rate(int compartment,double dA)
{
    actDeathRate[compartment] = MAX(dA,actDeathRate[compartment]);
}

int strain::add_infected(double infTime)
{
    assocCells++;
    totalInfected++;
    totalCells++;

    //if (strainBeta == 0)
	//fprintf(stderr,"Warning: beta=0 and infection at t=%3.2lf\n",infTime);

    lastEventDate = infTime;

    return assocCells;
}

int strain::add_prolif(double prolifTime)
{
    assocCells++;
    totalProlif++;
    totalCells++;

    lastEventDate = prolifTime;

    return assocCells;
}

int strain::dec_assoc(settings *params, generic_list<strain *> *ql, double decTime)
{
    if (assocCells > 0)
	assocCells--;

    if (assocCells == 0)
    {
	expireDate = decTime;
	for (int i=0; i < MAX_COMPARTMENTS; i++)
	    if (activeCells[i] > 0)
		fprintf(stderr,"Warning: left over active cells!\n");

	// If youngest dies out, reset youngest to next 
	// strain w/ assoc cells
	if (strainNum == params->youngQ)
	{
	    fprintf(stderr,"Note: youngest strain %d died out at %g\n",
		strainNum,expireDate);

	    // if founder dies out, then skip over "junk" strains
	    int nextNum = MAX(strainNum+1,params->junk_bins+1);
	    for (int s=nextNum; s < ql->get_num_elems();s++)
	    {
		strain *the_strain = ql->get_elem(s);
		if (the_strain->get_assoc() > 0)
		{
		    params->youngQ = s;
		    fprintf(stderr,"Note: youngest strain set to %d\n", s);
		    break;
		}
	    }
	}
    }
    totalDeaths++;
    return assocCells;
}
int strain::inc_assoc()
{
    assocCells++;

    return assocCells;
}
int strain::get_assoc()
{
    return assocCells;
}
int strain::inc_lumped_strains(int compartment)
{
    lumpedStrains[compartment]++;

    return lumpedStrains[compartment];
}

int strain::dec_ecl_cells(int compartment)
{
    if (eclipseCells[compartment] > 0)
	eclipseCells[compartment]--;

    return eclipseCells[compartment];
}

int strain::dec_act_cells(int compartment, double decTime)
{
    if (activeCells[compartment] > 0)
	activeCells[compartment]--;

    bool any_active = false;
    for (int i=0; i < MAX_COMPARTMENTS; i++)
	if (activeCells[i] > 0)
	    any_active = true;

    if (!any_active)
	deactiveDate = decTime;

    return activeCells[compartment];
}
int strain::inc_act_cells(int compartment, double incTime)
{
    bool any_active = false;
    for (int i=0; i < MAX_COMPARTMENTS; i++)
	if (activeCells[i] > 0)
	    any_active = true;

    if (activeDate < 0)
	activeDate = incTime;

    else if (!any_active)
    {
	// on reactivation, reset "activation" date to now or a bit earlier
	// (time since recent deactivation pushes the activation date back a bit)
	double tequiv = 2 * deactiveDate - incTime;
	double old_act = activeDate;
	if (tequiv > activeDate)
	    activeDate = incTime - (tequiv-activeDate);
	else
	    activeDate = incTime;

	//fprintf(stderr,
	 //   "Strain %d reactivated at t=%3.2lf (old act=%3.2lf, deact at %3.2lf, newact = %3.2lf)\n",
	  //  strainNum, incTime,old_act,deactiveDate,activeDate);
    }

    activeCells[compartment]++;

    if (cd8Cells[compartment] == 0)
	cd8Cells[compartment]=1;

    return activeCells[compartment];
}

int strain::inc_ecl_cells(int compartment)
{
    eclipseCells[compartment]++;

    return eclipseCells[compartment];
}

int strain::dec_l1_cells(int compartment)
{
    if (l1Cells[compartment] > 0)
	l1Cells[compartment]--;

    return l1Cells[compartment];
}
int strain::inc_l1_cells(int compartment)
{
    l1Cells[compartment]++;

    return l1Cells[compartment];
}
int strain::dec_l2_cells(int compartment)
{
    if (l2Cells[compartment] > 0)
	l2Cells[compartment]--;

    return l2Cells[compartment];
}
int strain::inc_l2_cells(int compartment)
{
    l2Cells[compartment]++;

    return l2Cells[compartment];
}
int strain::dec_via_virus(int compartment, int change)
{
    if (via_virions[compartment] > (unsigned int)change)
	via_virions[compartment]-=change;
    else
	via_virions[compartment]=0;

    return via_virions[compartment];
}
int strain::inc_via_virus(int compartment, int change)
{
    via_virions[compartment]+=change;

    return via_virions[compartment];
}
int strain::dec_junk_virus(int compartment, int change)
{
    if (junk_virions[compartment] > (unsigned int)change)
	junk_virions[compartment]-=change;
    else
	junk_virions[compartment]=0;

    return junk_virions[compartment];
}
int strain::inc_junk_virus(int compartment, int change)
{
    junk_virions[compartment]+=change;

    return junk_virions[compartment];
}
int strain::dec_held_viable(int compartment, int change)
{
    if (held_viable[compartment] > (unsigned int)change)
	held_viable[compartment]-=change;
    else
	held_viable[compartment]=0;

    return held_viable[compartment];
}
int strain::inc_held_viable(int compartment, int change)
{
    held_viable[compartment]+=change;

    return held_viable[compartment];
}
int strain::dec_held_junk(int compartment, int change)
{
    if (held_junk[compartment] > (unsigned int)change)
	held_junk[compartment]-=change;
    else
	held_junk[compartment]=0;

    return held_junk[compartment];
}
int strain::inc_held_junk(int compartment, int change)
{
    held_junk[compartment]+=change;

    return held_junk[compartment];
}
unsigned long int strain::dec_cd8s(int compartment, int change)
{
    if (cd8Cells[compartment] > (unsigned long int)change)
	cd8Cells[compartment]-=change;
    else if (assocCells > 0)
	cd8Cells[compartment]=1;
    else
	cd8Cells[compartment]=0;

    return cd8Cells[compartment];
}
unsigned long int strain::inc_cd8s(int compartment, int change)
{
    cd8Cells[compartment]+=change;

    if (cd8Cells[compartment] > cd8MaxCells[compartment])
	cd8MaxCells[compartment]=cd8Cells[compartment];

    return cd8Cells[compartment];
}
unsigned long int strain::set_cd8s(int compartment, unsigned long int total)
{
    cd8Cells[compartment]=total;

    if (cd8Cells[compartment] > cd8MaxCells[compartment])
	cd8MaxCells[compartment]=cd8Cells[compartment];

    return cd8Cells[compartment];
}
unsigned long int strain::get_cd8_cells(int compartment)
{
    return cd8Cells[compartment];
}
unsigned long int strain::get_max_cd8_cells(int compartment)
{
    return cd8MaxCells[compartment];
}
int strain::get_act_cells(int compartment)
{
    return activeCells[compartment];
}
int strain::get_ecl_cells(int compartment)
{
    return eclipseCells[compartment];
}
int strain::get_l1_cells(int compartment)
{
    return l1Cells[compartment];
}
int strain::get_l2_cells(int compartment)
{
    return l2Cells[compartment];
}
unsigned int strain::get_via_virus(int compartment)
{
    return via_virions[compartment];
}
unsigned int strain::get_junk_virus(int compartment)
{
    return junk_virions[compartment];
}
unsigned int strain::get_held_viable(int compartment)
{
    return held_viable[compartment];
}
unsigned int strain::get_held_junk(int compartment)
{
    return held_junk[compartment];
}
double strain::get_beta()
{
    return strainBeta;
}
double strain::get_dImmun_s()
{
    return dImmun_s;
}
double strain::get_ic50(settings *params, double time)
{

	return strainIC50;

}

double strain::get_birthdate()
{
    return strainBirth;
}
void strain::set_birthdate(double tnow)
{
    strainBirth = tnow;
}

double strain::get_last_event_date()
{
    return lastEventDate;
}

double strain::get_exp_date()
{
    return expireDate;
}
double strain::get_active_date()
{
    return activeDate;
}

strain *strain::get_parent()
{
    return parentStrain;
}
int strain::get_strain_index()
{
    return strainNum;
}
void strain::set_founder_dist(double fdist)
{
    founderDist=fdist;
}
double strain::get_founder_dist()
{
    return founderDist;
}
