//NOTE: This read-only file is generated from the file ../src/full_hiv_sim.cpp.
//Any editing should be done in that file.

#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <map>

using namespace std;
 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include "colors.h"

#ifdef GUI_ENABLED
#include "plotpoints.h"
#include "gui.h"
#endif

#include "cell.h"
#include "generic_list.h"
#include "strain.h"
#include "settings.h"

#define NUM_MEASURES 16

#define MAX_LINE 80
#define MAX_SET_POINTS 200
#define MAX_JUNK_BINS 50
#define BAD_SCORE 100

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define CHECK_FOR_REAL(name,param_name) \
    if (inputs.find(name) != inputs.end()) { params->param_name = inputs[name]; } else { \
	cerr << "No setting for optional parameter "<<name<<".  Left at default ("<<params->param_name<<")\n"; \
    }

#define CHECK_FOR_INT(name,param_name) \
    if (inputs.find(name) != inputs.end()) { params->param_name = static_cast<int>(inputs[name]); } else { \
	cerr << "No setting for optional parameter "<<name<<".  Left at default ("<<params->param_name<<")\n"; \
    }

void *my_malloc(int size)
{
    void *ret_val;
    ret_val=malloc(size);
    if (ret_val == NULL)
    {
	fprintf(stderr,"Encountered a problem allocating memory of size=%d. Exiting.\n",size);
	exit(1);
    }
    return ret_val;
}

//comparator necessary for qsort
int compare_doubles (const void * a, const void * b)
{
  if (*(double *)a > *(double *)b)
     return 1;
  if (*(double *)a < *(double *)b)
     return -1;
  return 0;
}

//######################################################################
//#### helper function for picking sequences based on samples active cells
// in a given compartment.  
// returns a list of sequence indices/cell count pairs
struct seq_info {
    int seq_index;
    int part_count;
};

strain *get_strain_for_virion(generic_list<strain *> *list_to_sample, int compartment, unsigned int virions)
{
    unsigned int curr_virions=0;
    for (int i=0; i < list_to_sample->get_num_elems(); i++)
    {
	strain *this_one = list_to_sample->get_elem(i);
	if (i!= 1)
	{
	    curr_virions += this_one->get_via_virus(compartment);
	    curr_virions += this_one->get_junk_virus(compartment);

	    if (curr_virions >= virions)
		return this_one;
	}
    }
    return NULL;
}
int sample_seqs_for_div(settings *params, generic_list<strain *> *list_to_sample, int compartment, int virions_to_use, struct seq_info **info)
{
    int num_active_seqs=0;

    unsigned int curr_virions = 0;
    int new_ones=0;

    if (list_to_sample->get_num_elems()<=0)
        return -1;

    // get circulating strains (based on sample volume)
    for (int i=0; i < list_to_sample->get_num_elems(); i++)
    {
	strain *this_one = list_to_sample->get_elem(i);
	if (i!= 1 || params->seq_sample_type == SAMPLE_DNA)
	{
	    
	    unsigned int this_seq_virions;
	    this_seq_virions = this_one->get_via_virus(compartment);
	    this_seq_virions += this_one->get_junk_virus(compartment);
	    if (this_seq_virions > 0)
	    {
		num_active_seqs++;
		curr_virions += this_seq_virions;
	    }

	}
    }
    if ((unsigned int)virions_to_use > curr_virions)
	virions_to_use = (int)curr_virions;

    if (num_active_seqs > 0 && virions_to_use > 0)
    {
	// malloc structure (to be freed by calling routine)
	// store up to <num_active_seqs> sequences
	*info = (struct seq_info *)my_malloc(num_active_seqs*sizeof (struct seq_info));
	// check each virion in the gathered sample and report its sequence
	for (int i=0; i < virions_to_use; i++)
	{
	    unsigned int spot = (unsigned int)gsl_rng_uniform_int(params->ur,curr_virions);
	    strain *this_one=get_strain_for_virion(list_to_sample,compartment,spot);
	    if (this_one != NULL)
	    {
		bool found = false;
		for (int j=0; j < new_ones; j++)
		    if ((*info)[j].seq_index == this_one->get_strain_index())
		    {
			found=true;
			break;
		    }
		if (!found)
		{
		    unsigned int seq_virions;
		    seq_virions = this_one->get_via_virus(compartment);
		    seq_virions += this_one->get_junk_virus(compartment);
		    (*info)[new_ones].part_count=seq_virions;
		    (*info)[new_ones].seq_index = this_one->get_strain_index();
		    new_ones++;
		}
	    }
	}
    }

    return new_ones;
}
int new_sample_seqs(settings *params, generic_list<strain *> *list_to_sample, int compartment, int target_seqs, int vol_samples, struct seq_info **info)
{
    int num_active_seqs=0;
    double vol_ratio = vol_samples/params->vol;

    unsigned int curr_virions = 0;
    unsigned int scaled_tot_virions = 0;
    int new_ones=0;

    if (list_to_sample->get_num_elems()<=0)
        return -1;

    // get circulating strains (based on sample volume)
    for (int i=0; i < list_to_sample->get_num_elems(); i++)
    {
	strain *this_one = list_to_sample->get_elem(i);
	if (i!= 1 || params->seq_sample_type == SAMPLE_DNA)
	{
	    
	    unsigned int this_seq_virions;
	    this_seq_virions = this_one->get_via_virus(compartment);
	    this_seq_virions += this_one->get_junk_virus(compartment);
	    if (this_seq_virions > 0)
	    {
		num_active_seqs++;
		curr_virions += this_seq_virions;
	    }

	}
    }
    scaled_tot_virions = curr_virions * vol_ratio;
    if (num_active_seqs > 0 && scaled_tot_virions > 0)
    {
	// malloc structure (to be freed by calling routine)
	// store up to <target_seqs> sequences
	*info = (struct seq_info *)my_malloc(target_seqs*sizeof (struct seq_info));
	// try to fill first by randomly grabbing a virion and storing its parent
	for (unsigned int i=0; i < scaled_tot_virions; i++)
	{
	    unsigned int spot = (unsigned int)gsl_rng_uniform_int(params->ur,curr_virions);
	    strain *this_one=get_strain_for_virion(list_to_sample,compartment,spot);
	    if (this_one != NULL)
	    {
		bool found = false;
		for (int j=0; j < new_ones; j++)
		    if ((*info)[j].seq_index == this_one->get_strain_index())
		    {
			found=true;
			break;
		    }
		if (!found)
		{
		    unsigned int seq_virions;
		    seq_virions = this_one->get_via_virus(compartment);
		    seq_virions += this_one->get_junk_virus(compartment);
		    (*info)[new_ones].part_count=seq_virions;
		    (*info)[new_ones].seq_index = this_one->get_strain_index();
		    new_ones++;
		    if (new_ones == target_seqs)
			break;
		}
	    }
	}
	// try to fill rest by grabbing any other strains present (in any order)
	if (new_ones < target_seqs && num_active_seqs > new_ones)
	{
	    int start = gsl_rng_uniform_int(params->ur,list_to_sample->get_num_elems());
	    int i=start;
	    int tried = 0;
	    while (tried < list_to_sample->get_num_elems())
	    {
		unsigned int seq_virions;
		if (i == list_to_sample->get_num_elems())
		    i=0;

		if (i==1)
		{
		    i++;
		    tried++;
		    continue;
		}

		strain *this_one = list_to_sample->get_elem(i);
		seq_virions = this_one->get_via_virus(compartment);
		seq_virions += this_one->get_junk_virus(compartment);

		if (seq_virions == 0)
		{
		    i++;
		    tried++;
		    continue;
		}

		bool found = false;

		for (int j=0; j < new_ones; j++)
		    if ((*info)[j].seq_index == this_one->get_strain_index())
		    {
			found=true;
			break;
		    }
		if (!found)
		{
		    (*info)[new_ones].part_count=seq_virions;
		    (*info)[new_ones].seq_index = this_one->get_strain_index();
		    new_ones++;
		    if (new_ones == target_seqs)
			break;
		}
		i++;
		tried++;
	    }
	}
    }

    return new_ones;
}

void dec_cd8_inf(settings *params, int compartment, int cd8_index,double time) 
{
    if (params->cd8_groups[compartment][cd8_index].inf_cells > 0)
	params->cd8_groups[compartment][cd8_index].inf_cells--;
    else {
	fprintf(stderr,
	    "Error: inf cell count for HD=%d would go negative at t=%lf!\n",
	    cd8_index,time);
	exit(1);
    }
}
// called for active and/or latent cells (both types) on infection event
//
void process_infection(settings *params, double time,
	bool junk_infection, int compartment, int cell_type, 
	unsigned long int *cd8s, generic_list<strain *> *ql,
	strain *the_strain, int *Q)
{
    //double junk_chance = gsl_rng_uniform(params->ur);
    double founderDist = the_strain->get_founder_dist();

    // just track type L2 junk cells (for ART fossil record)
    //if (junk_chance < params->junk)
    if (junk_infection)
    {
	strain *junk_strain;
	int junk_index;
	if (time < params->junk_days)
	    junk_index = (int)(time / ((double)params->junk_days/params->junk_bins)) + 1;
	else
	    junk_index = params->junk_bins+1;

	junk_strain = ql->get_elem(junk_index);

	if (junk_strain->get_assoc() == 0)
	    junk_strain->set_birthdate(time);

	if (cell_type == LATENT1)
	{
	    junk_strain->inc_l1_cells(compartment);
	}
	else if (cell_type == LATENT2)
	{
	    junk_strain->inc_l2_cells(compartment);
	}
	else
	{
	    junk_strain->inc_act_cells(compartment,time);
	    if ((params->immun_model == 4 || params->immun_model == 5) && params->kill_junk) {
		params->cd8_groups[compartment][0].inf_cells++;
	    }
	}

	double deltaHam = gsl_ran_weibull(params->ur,40,5);
	founderDist += deltaHam;

	double avg_dH = (junk_strain->get_assoc() * junk_strain->get_founder_dist()+founderDist) /
		(junk_strain->get_assoc()+1);

	//junk_strain->set_founder_dist(MAX(founderDist,junk_strain->get_founder_dist()));
	junk_strain->set_founder_dist(avg_dH);

	the_strain = junk_strain;
	//the_strain->inc_assoc();
	the_strain->inc_lumped_strains(compartment);
	the_strain->add_infected(time);
	founderDist = the_strain->get_founder_dist();
    }
    else 
    {
	if (gsl_rng_uniform(params->ur)<params->mu)
	{
	    double strainBeta=0;
	    if (params->use_entropy) /*&& params->immun_model != 6)*/
	    {
		// fitness will be adjusted in "new Strain" based on entropy of random sites in the consensus strain (where SNPs occur)
		strainBeta = the_strain->get_beta();
	    }
	    else if (params->Bt_k == 0 && params->use_weibull)
		strainBeta = params->Bt_mean;

	    else if (params->Bt_sd_factor  != 0)
	    {
		//strainBeta = params->Bt_mean + gsl_ran_gaussian (params->ur, params->Bt_mean * params->Bt_sd_factor);
		strainBeta = the_strain->get_beta() + gsl_ran_gaussian (params->ur, the_strain->get_beta() * params->Bt_sd_factor);
	    }

	    // draw beta change (mult) from a skewed weibull distribution (mean=1)
	    else if (params->use_dBeta)
	    {
		double strainMult;
		strainBeta = the_strain->get_beta();

		if (params->use_weibull)
		    strainMult = gsl_ran_weibull(params->ur, params->Bt_mean/params->Bt_max,params->Bt_k);
		else
		{
		    strainMult = pow(10,gsl_ran_gaussian(params->ur, params->Bt_k));
		}
				    //params->Bt_mean/strainBeta,params->Bt_k);
		strainBeta = strainMult * strainBeta;
	    }

	    // draw beta from a skewed weibull distribution (k as param)
	    else if (params->use_weibull)
		strainBeta = gsl_ran_weibull(params->ur,params->Bt_mean,params->Bt_k);
	    // pick a random beta between 0 and the Bt_max (if defined) or
	    // 0 and Bt_mean if no max is defined
	    else
	    {
		double strainMult;
		strainMult = gsl_rng_uniform(params->ur);
		if (params->Bt_max > 0)
		    strainBeta = strainMult * params->Bt_max;
		else
		    strainBeta = strainMult * params->Bt_mean;
	    }

	    // limit to max if set
	    if (params->Bt_max != 0 && strainBeta > params->Bt_max)
		strainBeta=params->Bt_max;

	    // don't allow negative infectivity!
	    if (strainBeta < 0)
		strainBeta=0;

	    if (strainBeta > 0 || cell_type != ACTIVE)
	    {
		(*Q)++;

		strain *new_strain = new strain(params,ql->get_num_elems(),time,the_strain,strainBeta);
		ql->append(new_strain);

		the_strain = new_strain;
		founderDist = the_strain->get_founder_dist();
	    }

	    if (params->immun_model == 2 || params->immun_model == 3)
		(*cd8s) += the_strain->get_cd8_cells(compartment);
	}
	int infected = the_strain->add_infected(time);
	if (infected > params->max_strain_cells)
	    params->max_strain_cells=infected;

	if (cell_type == ACTIVE)
	{
	    if (params->eclipse > 0)
		the_strain->inc_ecl_cells(compartment);
	    else
		the_strain->inc_act_cells(compartment,time);

	    if (params->immun_model == 4 || params->immun_model == 5) {
		int cd8_index;
		if (params->immun_model == 4)
		    cd8_index = 1+MIN(params->maxHammingDist/params->cd8GroupSize,the_strain->get_founder_dist()/params->cd8GroupSize);
		else if (params->immun_model == 5)
		    cd8_index = 1+MIN(params->tF/params->cd8GroupSize,the_strain->get_birthdate()/params->cd8GroupSize);
		params->cd8_groups[compartment][cd8_index].inf_cells++;
	    }
	}
	else if (cell_type == LATENT1)
	    the_strain->inc_l1_cells(compartment);
	else
	    the_strain->inc_l2_cells(compartment);

    }
}
int get_strain_ns(settings *params, strain *the_strain)
{
    int non_syns=0;
    char *non_syn_ptr = the_strain->get_non_syn_ptr();

    if (non_syn_ptr==NULL)
	return 0;

    for (int j=0; j < (params->founder_stats.length / 8)+1; j++)
	for (int k=0; k < 8; k++)
	    if ((non_syn_ptr[j] & (1 << k))  != 0)
		non_syns++;
    return non_syns;
}

//#####################################################################
//#### function that gathers statistics based on all sequences to date 

void gather_all_stats(settings *params, double time,
	unsigned int x[],
	generic_list<strain *> *ql,
	double *Bt_avg, double *Bt_max, 
	int *q_a,
	int *qmax, 
	int *qmax_via, 
	int *qmax_junk, 
	int *qmax_a, 
	int *qmax_l1, 
	int *qmax_l2, 
	int *qmax_dH, 
	unsigned long int *qmax_cd8s, 
	int *qmax_cd8_index,
	double *qmax_dA, 
	double *qmax_beta, 
	int *qmax_ns, 
	double *avg_hd, 
	int *max_hd, 
	double *qmax_dImmun,
	double *dImmun_avg) 
{
    int A=x[1]; //#active infected viable cells
    int startQ=x[8]; //#infectious sequences at start of frame

    // use lifetimes of all active strains to get new beta
    *Bt_avg = 0;
    *Bt_max = 0;
    *avg_hd = 0;
    *max_hd = 0;
    *dImmun_avg = 0;

    int curr_qmax_via = 0;
    int curr_qmax_junk = 0;
    int curr_qmax_a = 0;
    int curr_qmax_l1 = 0;
    int curr_qmax_l2 = 0;
    int curr_qmax = 0;
    int curr_q_a = 0;
    int curr_qmax_dH = 0;
    unsigned long int curr_qmax_cd8s = 0;
    int curr_qmax_cd8_index = 0;
    double curr_qmax_dA = 0;
    double curr_qmax_beta = 0;
    double curr_qmax_dImmun = 0;
    int curr_qmax_ns = 0;

    int num_compartments = 1;


    for (int s=0; s < startQ; s++)
    {
	if (s>=1 && s<=params->junk_bins+1)
	    continue;

	strain *the_strain = ql->get_elem(s);

	double this_beta = the_strain->get_beta();

	int this_ns = get_strain_ns(params, the_strain);

	double dImmun = params->dImmun;

	if (params->immun_model == 6)
	    dImmun = the_strain->get_dImmun_s();

	unsigned int num_virions = 0;

	unsigned int num_junk_virions = 0;

	int num_act_cells = 0;

	int num_l1_cells = 0;

	int num_l2_cells = 0;

	unsigned long int num_cd8s = 0;

	double highest_dA = 0;

	int cd8_index=0;

	for (int compartment=0; compartment < num_compartments; compartment++)
	{
	    //double this_ic50 = the_strain->get_ic50();

	    /* determine #cells and virons from LAST frame */
	    num_virions += the_strain->get_via_virus(compartment);

	    num_junk_virions += the_strain->get_junk_virus(compartment);

	    num_act_cells += the_strain->get_act_cells(compartment);

	    num_l1_cells += the_strain->get_l1_cells(compartment);

	    num_l2_cells += the_strain->get_l2_cells(compartment);

	    if (params->immun_model == 4 || params->immun_model == 5) {
		if (params->immun_model == 4)
		    cd8_index = 1+MIN(params->maxHammingDist/params->cd8GroupSize,the_strain->get_founder_dist()/params->cd8GroupSize);
		else if (params->immun_model == 5)
		    cd8_index = 1+MIN(params->tF/params->cd8GroupSize,the_strain->get_birthdate()/params->cd8GroupSize);
		num_cd8s += params->cd8_groups[compartment][cd8_index].cd8_cells;
	    }
	    else
		num_cd8s += the_strain->get_cd8_cells(compartment);

	    double this_dA = the_strain->get_death_rate(compartment,params,time);

	    if (this_dA > highest_dA)
		highest_dA = this_dA;
	}
	double dH = the_strain->get_founder_dist();

	if (num_act_cells > 0)
	{
	    if (this_beta > (*Bt_max))
		*Bt_max = this_beta;

	    *Bt_avg += this_beta*num_act_cells;
	    *dImmun_avg += dImmun*num_act_cells;
	    curr_q_a++;

	    if (dH > *max_hd)
		*max_hd = dH;

	    *avg_hd += dH;

	    if (num_act_cells > curr_qmax_a)
	    {
		curr_qmax_a=num_act_cells;
		curr_qmax_via=num_virions;
		curr_qmax_junk=num_junk_virions;
		curr_qmax_l1=num_l1_cells;
		curr_qmax_l2=num_l2_cells;
		curr_qmax_dH=dH;
		curr_qmax=s;
		curr_qmax_cd8s=num_cd8s;
		curr_qmax_beta=this_beta;
		curr_qmax_dImmun=dImmun;
		curr_qmax_ns=this_ns;
		curr_qmax_dA=highest_dA;
		curr_qmax_cd8_index=cd8_index;
	    }
	}
    }
    if (curr_q_a > 0)
	*avg_hd = *avg_hd/curr_q_a;
    (*q_a) =  curr_q_a;
    (*qmax) =  curr_qmax;
    (*qmax_a) =  curr_qmax_a;
    (*qmax_via) =  curr_qmax_via;
    (*qmax_junk) =  curr_qmax_junk;
    (*qmax_l1) =  curr_qmax_l1;
    (*qmax_l2) =  curr_qmax_l2;
    (*qmax_dH) =  curr_qmax_dH;
    (*qmax_cd8s) =  curr_qmax_cd8s;
    (*qmax_beta) =  curr_qmax_beta;
    (*qmax_dImmun) =  curr_qmax_dImmun;
    (*qmax_dA) =  curr_qmax_dA;
    (*qmax_ns) =  curr_qmax_ns;
    (*qmax_cd8_index) =  curr_qmax_cd8_index;

    // Betas and dImmun should have been calculated and summed for each strain 
    // with active cells (tot=A)

    if (A > 1)
    {
	*Bt_avg = *Bt_avg / A;
	*dImmun_avg = *dImmun_avg / A;
    }
}
//#####################################################################
//#### function that gathers statistics based on sequences samples 
//(with actively infected cells or HIV DNA based on seq_sample_type)

void gather_sample_stats(settings *params, generic_list<strain *> *ql,
	int num_seqs,
	struct seq_info *seq_samples,
	double *samp_avgHd,
	int *samp_maxHd,
	double *entropy,
	double *diversity,
	double *ttca_mean,
	double *seq_ident)
{
    unsigned int tot_virions=0;
    double ttca=0;
    int ttca_lengths=0;

    // use lifetimes of all active strains to get new beta
    *samp_avgHd = 0;
    *entropy = 0;
    *diversity = 0;
    *ttca_mean = 0;
    *seq_ident = 0;

    for (int s=0; s < num_seqs; s++)
    {
	tot_virions += seq_samples[s].part_count;
    }

    int abs_maxHD=0;
    //int anc_minHD=params->maxHammingDist;

    for (int s=0; s < num_seqs; s++)
    {
	strain *the_strain = ql->get_elem(seq_samples[s].seq_index);

	//int num_act_cells = seq_samples[s].part_count;
	int num_virions = seq_samples[s].part_count;

	if (the_strain->get_parent() == NULL)
	    (*seq_ident) = num_virions;

	if (num_virions > 0)
	{
	    double dH = the_strain->get_founder_dist();
	    double bday=the_strain->get_birthdate();
	    double this_strains_diversity=0;
	    int *parents_this=NULL;
	    parents_this=(int *)my_malloc(params->maxHammingDist*sizeof(int));
	    strain *next_str=the_strain;
	    int slot=0;
	    while (slot < params->maxHammingDist && next_str->get_parent() != NULL)
	    {
		parents_this[slot]=next_str->get_strain_index();
		slot++;
		next_str=next_str->get_parent();
	    }
	    if (slot < params->maxHammingDist)
		parents_this[slot]=0;
		
	    for (int s2=s; s2 < num_seqs; s2++)
	    {
                int common_parent=0;
		bool found=false;
		strain *second_strain = ql->get_elem(seq_samples[s2].seq_index);

		for (int t=0; t < slot; t++)
		{
		    next_str=second_strain;
		    
		    while (found==false && next_str->get_parent() != NULL)
		    {
			if (parents_this[t]==next_str->get_strain_index())
			{
			    found=true;
			    common_parent=parents_this[t];
			}
			next_str=next_str->get_parent();
		    }
		}
		strain *common_ancestor=ql->get_elem(common_parent);
		int dH_common = common_ancestor->get_founder_dist();
		int dH_second = second_strain->get_founder_dist();
		double bday_common=common_ancestor->get_birthdate();
		double bday_second=second_strain->get_birthdate();
		if (dH_common > dH || dH_common > dH_second)
		{
		     fprintf(stderr,"Error: hamming distance error )common=%d, seq=%lf, 2nd=%d)\n",
			dH_common,dH,dH_second);
		}
		else
		{
		    ttca+=(bday-bday_common)+(bday_second-bday_common);
		    ttca_lengths+=2;
		}

		// identify common ancestor closest to founder
		//if (dH_common < anc_minHD)
		    //anc_minHD = dH_common;

		this_strains_diversity += (dH + dH_second-2*dH_common) * seq_samples[s2].part_count / ((double)tot_virions * params->seq_length);
	    }
	    (*diversity) += this_strains_diversity * num_virions / (double)tot_virions;

	    (*ttca_mean) = ttca/ttca_lengths;

	    (*entropy) -= ((double)num_virions/tot_virions) * log((double)num_virions/tot_virions);


	    if (dH > abs_maxHD)
		abs_maxHD = dH;

	    *samp_avgHd += dH*num_virions;
	    delete(parents_this);
	}
    }
    //*samp_maxHd = abs_maxHD-anc_minHD;
    *samp_maxHd = abs_maxHD;
    *samp_avgHd = *samp_avgHd / tot_virions; // based on absolute HD

    (*diversity) = 200*(*diversity); // double since dividing by N^2 and make %
    (*seq_ident) = (*seq_ident) / tot_virions;
}
//#####################################################################
//#### function that gathers fractions of DNA based on cell counts
//in the latent reservoir

void gather_reserv_stats(settings *params, generic_list<strain *> *ql,double time)
{
    unsigned int tot_junk_cells=0;
    unsigned int tot_via_cells=0;
    unsigned int tot_cells=0;
    unsigned int curr_via_cells=0;
    int slot=0;

    int junk_by_slot[MAX_JUNK_BINS];
    int viable_by_slot[MAX_JUNK_BINS];

    int junk_strains[MAX_JUNK_BINS];
    int viable_strains[MAX_JUNK_BINS];

    for (int y=0; y < params->junk_bins; y++)
    {
	junk_by_slot[y]=0;
	viable_by_slot[y]=0;
	junk_strains[y]=0;
	viable_strains[y]=0;
    }

    int junk_index;
    if (time < params->junk_days)
	junk_index = (int)(time / ((double)params->junk_days/params->junk_bins)) + 1;
    else
	junk_index = params->junk_bins+1;

    for (int junkies=1; junkies < junk_index; junkies++)
    {
	strain *junk_strain = ql->get_elem(junkies);
	slot = junkies-1;
	int junk_l0 = junk_strain->get_act_cells(0);
	int junk_l1 = junk_strain->get_l1_cells(0);
	int junk_l2 = junk_strain->get_l2_cells(0);

	junk_by_slot[slot]=junk_l0+junk_l1+junk_l2;
	tot_junk_cells+=junk_l0+junk_l1+junk_l2;
    }

    for (int s=0; s < ql->get_num_elems(); s++)
    {
	if (s==0 || s >= junk_index)
	{
	    strain *the_strain = ql->get_elem(s);
	    double birthday = the_strain->get_birthdate();
	    slot =(int)(birthday / ((double)params->junk_days/params->junk_bins)); 
	    int num_l0 = the_strain->get_act_cells(0);
	    int num_l1 = the_strain->get_l1_cells(0);
	    int num_l2 = the_strain->get_l2_cells(0);
	    if (num_l0+num_l1+num_l2 > 0)
	    {
		viable_by_slot[slot]+=num_l0+num_l1+num_l2;
		tot_via_cells += num_l0+num_l1+num_l2;
	    }
	}
    }
    tot_cells=tot_junk_cells+tot_via_cells;

    int samples = tot_cells * params->dna_for_res / params->vol; // 100 ul / 10 ml?

    double *ran_fracts;
    ran_fracts=(double *)my_malloc(samples*sizeof(double));

    for (int i=0; i < samples; i ++)
	ran_fracts[i] = gsl_rng_uniform(params->ur);

    qsort (ran_fracts, samples, sizeof(double), compare_doubles);

    int junk_end = samples;

    for (int s=0; s < samples; s++)
    {
	double junk_fract=0;
	bool is_junk=false;
	for (int bin=0; bin < params->junk_bins; bin++) 
	{
	    junk_fract+=junk_by_slot[bin];
	    if (ran_fracts[s] < (double)junk_fract/tot_cells)
	    {
		junk_strains[bin]++;
		is_junk=true;
		break;
	    }
	}
	if (!is_junk)
	    junk_end=s;
    }

    int samp=junk_end; 
    int extras=0;
    for (int s=0; s < ql->get_num_elems(); s++)
    {
	if (s==0 || s >= junk_index)
	{
	    strain *the_strain = ql->get_elem(s);
	    double birthday = the_strain->get_birthdate();
	    slot =(int)(birthday / ((double)params->junk_days/params->junk_bins)); 
	    if (slot > params->junk_bins-1)
		slot=params->junk_bins-1;

	    int num_l0 = the_strain->get_act_cells(0);
	    int num_l1 = the_strain->get_l1_cells(0);
	    int num_l2 = the_strain->get_l2_cells(0);
	    if (num_l0+num_l1+num_l2 > 0)
	    {
		curr_via_cells += num_l0+num_l1+num_l2;
		if (ran_fracts[samp] < (double)(tot_junk_cells+curr_via_cells)/tot_cells)
		{
		    viable_strains[slot]++;
		    // don't double count strains for multiple cell hits
		    if ( samp < samples &&
			ran_fracts[samp] < (double)(tot_junk_cells+curr_via_cells)/tot_cells)
			samp++;

		    while ( samp < samples &&
			ran_fracts[samp] < (double)(tot_junk_cells+curr_via_cells)/tot_cells)
		    {
			samp++;
			extras++;
		    }
		}
	    }
	}
    }
	
    int tot_strains=0;
    int tot_junk_strains=0;
    int tot_via_strains=0;
    for (int i=0; i < params->junk_bins ; i++)
    {
	tot_junk_strains+=junk_strains[i];
	tot_via_strains+=viable_strains[i];
    }
    tot_strains=tot_junk_strains+tot_via_strains;
    fprintf(stdout,
	"Processed %d samples from %d strains, %d junk cells and %d viable cells (%d duplicates)\n",
	samples,tot_strains,junk_end,samples-junk_end,extras);

    fprintf(stdout,"True DNA cell fractions at t=%g",time);
    for (int i=0; i < params->junk_bins ; i++)
	fprintf(stdout,",%g",(double)(junk_by_slot[i]+viable_by_slot[i])/tot_cells);
    fprintf(stdout,"\n");

    fprintf(stdout,"True Junk cell fractions at t=%g",time);
    for (int i=0; i < params->junk_bins ; i++)
	fprintf(stdout,",%g",(double)(junk_by_slot[i])/tot_junk_cells);
    fprintf(stdout,"\n");

    fprintf(stdout,"True Viable cell fractions at t=%g",time);
    for (int i=0; i < params->junk_bins ; i++)
	fprintf(stdout,",%g",(double)(viable_by_slot[i])/tot_via_cells);
    fprintf(stdout,"\n");

    fprintf(stdout,"Sampled DNA strain fractions at t=%g",time);
    for (int i=0; i < params->junk_bins ; i++)
	fprintf(stdout,",%g",(double)(junk_strains[i]+viable_strains[i])/tot_strains);
    fprintf(stdout,"\n");

    fprintf(stdout,"Sampled Junk strain fractions at t=%g",time);
    for (int i=0; i < params->junk_bins ; i++)
	fprintf(stdout,",%g",(double)(junk_strains[i])/tot_junk_strains);
    fprintf(stdout,"\n");

    fprintf(stdout,"Sampled Viable strain fractions at t=%g",time);
    for (int i=0; i < params->junk_bins ; i++)
	fprintf(stdout,",%g",(double)(viable_strains[i])/tot_via_strains);
    fprintf(stdout,"\n");
    fflush(stdout);
}

void model(settings *params, double time, int compartment,
	unsigned int x[],
	generic_list<strain *> *ql,
	generic_list<strain *> *top_strains
	)
{
            
    int S=x[0]; //#susceptible cells

    int A=x[1]; //#active infected cells

    int L1=x[2]; //#latent infected cells type 1

    int L2=x[3]; //#latent infected cells type 2

    int Vvia=x[4]; //#infectious virus
    int Vnon=x[5]; //#infectious virus
    int Vtot=Vvia+Vnon; //#all current virus (for CD4 expansion term)

    int It=x[6]; //#infectious cells ever made
    int Id=x[7]; //#infectious cell deaths

    int Q=x[8]; //#infectious sequences
    int startQ=x[8]; //#infectious sequences at start of frame

    unsigned long int CD8s=x[9]; //#total CD8s

    int junk_L0s=x[10]; //#total cells with junk DNA (reg lifespan)
    int junk_L1s=x[11]; //#total cells with junk DNA (reg lifespan)
    int junk_L2s=x[12]; //#total cells with junk DNA (long lifespan)

    int E=x[13]; //#pre-productive infected cells


    //Viral compartment
    //#poisson number of new virions per day from actively infected cells

    int Vnon_births=0;
    int Vvia_births=0;
    int Vnon_deaths = 0;
    int Vvia_deaths = 0;

    int L1_infs=0;
    int L1_acts=0;
    int L1_births=0;
    int L1_deaths=0;

    int L2_infs=0;
    int L2_acts=0;
    int L2_births=0;
    int L2_deaths=0;

    int A_infs=0;

    int A_ripens=0;

    int E1_acts=0;

    int E2_acts=0;

    int A_deaths=0;

    int E_infs=0;

    int CD8_births=0;
    int CD8_deaths=0;

    int junk_L0_deaths=0;
    int junk_L0_infs=0;

    int junk_L1_births=0;
    int junk_L1_deaths=0;
    int junk_L1_infs=0;

    int junk_L2_births=0;
    int junk_L2_deaths=0;
    int junk_L2_infs=0;


    generic_list<strain *> *top_this_frame=NULL;
    double ART_eff=0;
    double Bt_eff, dA_eff;
    double vol;

    if (params->followTopStrains)
	top_this_frame=new generic_list<strain *>();

    double dL1 = params->aL1+params->thL1+params->xi1;           //# death rate
    double dL2 = params->aL2+params->thL2+params->xi2;           //# death rate

    double dJ1 = params->aL1+params->thL1;   
    double dJ2 = params->aL2+params->thL2;  

    double tau1 = params->tau1;
    double tau2 = params->tau2;

    if (compartment == 0)
    {
	vol = params->vol;
	if (time >= params->ART_start && 
	    (params->ART_stop==0 || time <= params->ART_stop))
	    ART_eff=params->ART_eff;
	else
	    ART_eff=0;
	Bt_eff=1;
	dA_eff=1;
	params->numTopStrains=0;
    }
    double mult;
    if (params->res_model == 1 && ART_eff == 0)
    {
	mult = params->res_factor;
	dL1 = params->aL1+params->thL1*mult-params->xi1;   
	dL2 = params->aL2+params->thL2*mult-params->xi2;  
	dJ1 = params->aL1+params->thL1*mult;   
	dJ2 = params->aL2+params->thL2*mult;  
    }
    else if (params->res_model == 2 && Vtot > 0) 
    {
	mult = (1 + log10(Vtot) / (log10(Vtot)+params->res_factor));
	dL1 = params->aL1+params->thL1*mult-params->xi1;   
	dL2 = params->aL2+params->thL2*mult-params->xi2;  
	dJ1 = params->aL1+params->thL1*mult;   
	dJ2 = params->aL2+params->thL2*mult;  
    }
    // In this model we modify the probability of going latent
    // AND move that decision to the time of active cell death 
    else if (params->res_model == 3 && ART_eff > 0)
    {
	tau1 = tau1 * params->res_factor;
	tau2 = tau2 * params->res_factor;
    }


    for (int s=0; s < startQ; s++)
    {
	if (s > params->junk_bins+1 && s < params->youngQ)
	    continue;

	strain *the_strain = ql->get_elem(s);

	int cd8_index=1;

	if (params->immun_model == 4)
	    cd8_index = 1+MIN(params->maxHammingDist/params->cd8GroupSize,the_strain->get_founder_dist()/params->cd8GroupSize);
	else if (params->immun_model == 5)
	    cd8_index = 1+MIN(params->tF/params->cd8GroupSize,the_strain->get_birthdate()/params->cd8GroupSize);

	double dImmun = params->dImmun;

	if (params->immun_model == 6)
	    dImmun = the_strain->get_dImmun_s();

	unsigned long int num_cd8s = 0;

	if (params->immun_model == 2 || params->immun_model == 3)
	    num_cd8s = the_strain->get_cd8_cells(compartment);

	int num_act_cells = the_strain->get_act_cells(compartment);

	int num_l1_cells = the_strain->get_l1_cells(compartment);

	int num_l2_cells = the_strain->get_l2_cells(compartment);

	double this_ic50 = the_strain->get_ic50(params, time);

	// Junk strain Now given APOBEC bundled virus from infected cells (1/2019)
	if (s>=1 && s <= params->junk_bins+1)
	{
	    int killed = 0;
	    int killed_s = 0;

	    if (params->kill_junk &&
		(params->immun_model == 1 || params->immun_model == 6 || params->immun_model == 3)) {
		killed = gsl_ran_poisson (params->ur, 
		    num_act_cells*(CD8s-num_cd8s)*dImmun*params->dt);
	    }
	    if (params->kill_junk &&
		(params->immun_model == 2 || params->immun_model == 3)) {
		if(num_act_cells-killed > 0)
		    killed_s = gsl_ran_poisson (params->ur, 
			    (num_act_cells-killed)*params->dImmun_s*num_cd8s*params->dt);

		double FIp = num_act_cells/(num_act_cells+this_ic50);

		unsigned long int min_cd8s = num_cd8s;
		if (min_cd8s==0)
		    min_cd8s=1;

		int cd8_b = gsl_ran_poisson (params->ur,params->tbirth*params->dt);
		int cd8_exp = gsl_ran_poisson (params->ur,params->theta*FIp*min_cd8s*params->dt);

		unsigned long int cd8_d = 0;
		if (params->immun_model == 3)
		    cd8_d = gsl_ran_poisson (params->ur,(params->delta+params->kappa*(CD8s/(CD8s+params->cd8_ic50)))*num_cd8s*params->dt);
		else
		    cd8_d = gsl_ran_poisson (params->ur,params->delta*num_cd8s*params->dt);
		if (cd8_d > num_cd8s)
		    cd8_d = num_cd8s;

		the_strain->inc_cd8s(compartment,cd8_b+cd8_exp);

		the_strain->dec_cd8s(compartment,(int)cd8_d);

		if (params->immun_model == 3)
		{
		    CD8_births+=cd8_b+cd8_exp;
		    CD8_deaths+=cd8_d;
		}
	    }
	    // check for cd8 grouping option (strains share cd8s based on HD)
	    if ((params->immun_model == 4 || params->immun_model == 5) && params->kill_junk) {
		int cd8_cells = params->cd8_groups[compartment][0].cd8_cells;

		if (params->immun_model == 5)
		    killed = gsl_ran_poisson (params->ur, 
			num_act_cells*(CD8s-cd8_cells)*dImmun*params->dt);

		killed += gsl_ran_poisson (params->ur, 
		    num_act_cells*cd8_cells*params->dImmun_s*params->dt);
	    }
	    //#density dependent junk death (w/o activation)
	    int left_over = num_act_cells-killed-killed_s;

	    int j_l0_d=0;

	    // if the latent cells transition to memory INSTEAD of dying...
	    if (params->res_model == 3)
	    {
		int j1_s;
		j1_s = gsl_ran_poisson (params->ur, 
		    tau1*left_over*params->dt);

		if (left_over > j1_s)
		    left_over-=j1_s;
		else
		{
		    j1_s = left_over;
		    left_over=0;
		}

		//#latent L1 transitions from active 
		for (int j1=0; j1 < j1_s; j1++)
		{
		    the_strain->dec_act_cells(compartment,time);
		    the_strain->inc_l1_cells(compartment);
		    if (params->immun_model == 4 || params->immun_model == 5)
			params->cd8_groups[compartment][cd8_index].inf_cells--;
		    junk_L0_deaths ++;
		    junk_L1_infs++;
		}

		int j2_s;
		j2_s = gsl_ran_poisson (params->ur, 
		    tau2*left_over*params->dt);

		if (left_over > j2_s)
		    left_over-=j2_s;
		else
		{
		    j2_s = left_over;
		    left_over=0;
		}

		//#latent L2 transitions from active 
		for (int j2=0; j2 < j2_s; j2++)
		{
		    the_strain->dec_act_cells(compartment,time);
		    the_strain->inc_l2_cells(compartment);
		    if (params->immun_model == 4 || params->immun_model == 5)
			params->cd8_groups[compartment][cd8_index].inf_cells--;
		    junk_L0_deaths ++;
		    junk_L2_infs++;
		}
	    }

	    double Jdeath_rate = params->dActivated*left_over;

	    if (left_over > 0)
		j_l0_d = gsl_ran_poisson(params->ur,Jdeath_rate*params->dt);

	    if (j_l0_d > left_over)
		j_l0_d = left_over;

	    for (int j=0; j < j_l0_d+killed+killed_s; j++)
	    {
		the_strain->dec_act_cells(compartment,time);
		the_strain->dec_assoc(params,ql,time);
		if ((params->immun_model == 4 || params->immun_model == 5) && params->kill_junk)
		    dec_cd8_inf(params,compartment,0,time);
		junk_L0_deaths ++;
	    }

	    //#density dependent junk production 
	    int j_l1_b=0;
	    int j_l1_d=0;
	    double Jbirth_rate = params->aL1 * num_l1_cells;
	    if (num_l1_cells > 0)
		j_l1_b = gsl_ran_poisson(params->ur,Jbirth_rate*params->dt);

	    //#junk proliferations
	    for (int l=0; l < j_l1_b; l++)
	    {
		the_strain->inc_l1_cells(compartment);
		the_strain->add_prolif(time);
		junk_L1_births++;
	    }

	    //#density dependent junk death (w/o activation)
	    Jdeath_rate = dJ1*num_l1_cells;
	    //double Jdeath_rate = (params->dS)*junk_L1s; //12/28/17
	    if (num_l1_cells > 0)
		j_l1_d = gsl_ran_poisson(params->ur,Jdeath_rate*params->dt);

	    if (j_l1_d > num_l1_cells)
		j_l1_d = num_l1_cells;

	    for (int j=0; j < j_l1_d; j++)
	    {
		the_strain->dec_l1_cells(compartment);
		the_strain->dec_assoc(params,ql,time);
		junk_L1_deaths++;
	    }

	    //#density dependent junk production 
	    int j_l2_b=0;
	    int j_l2_d=0;
	    Jbirth_rate = params->aL2 * num_l2_cells;

	    if (num_l2_cells > 0)
		j_l2_b = gsl_ran_poisson(params->ur,Jbirth_rate*params->dt);

	    //#junk proliferations
	    for (int l=0; l < j_l2_b; l++)
	    {
		the_strain->inc_l2_cells(compartment);
		the_strain->add_prolif(time);
		junk_L2_births++;
	    }

	    //#density dependent junk death (w/o activation)
	    Jdeath_rate = dJ2*num_l2_cells;
	    if (num_l2_cells > 0)
		j_l2_d = gsl_ran_poisson(params->ur,Jdeath_rate*params->dt);

	    if (j_l2_d > num_l2_cells)
		j_l2_d = num_l2_cells;

	    for (int j=0; j < j_l2_d; j++)
	    {
		the_strain->dec_l2_cells(compartment);
		the_strain->dec_assoc(params,ql,time);
		junk_L2_deaths++;
	    }
	    continue;
	}

//   the following sections are skipped for "junk"

	double this_beta = Bt_eff * the_strain->get_beta()/vol;

	int num_ecl_cells = the_strain->get_ecl_cells(compartment);

	unsigned int num_virions = the_strain->get_via_virus(compartment);

	unsigned int num_junk_virions = the_strain->get_junk_virus(compartment);

	unsigned int num_tot_virions = num_virions+num_junk_virions;

	// Active cell deaths (from T-cell killing or natural death rate)
	//
	if (num_act_cells > 0)
	{
	    // first let current active cells produce virus!
	    int vb_s=gsl_ran_poisson(params->ur,num_act_cells*params->pi*params->dt);
	    // next check for APOBEC encapsulation
	    double junk_chance = gsl_rng_uniform(params->ur);

	    if (junk_chance < params->junk)
	    {
		the_strain->inc_junk_virus(compartment,vb_s);
		Vnon_births +=vb_s;
	    }
	    else
	    {
		the_strain->inc_via_virus(compartment,vb_s);
		Vvia_births +=vb_s;
	    }

	    // immune pressure is applied strain-wise based on numbers of 
	    // active cells which control the influx of t-cells

	    int killed = 0;
	    int killed_s = 0;

	    if (params->use_cd8s > 0)
	    {
		if (params->immun_model == 1 || params->immun_model == 6 || params->immun_model == 3)
		    killed = gsl_ran_poisson (params->ur, 
			num_act_cells*(CD8s-num_cd8s)*dImmun*params->dt);

		if ((params->immun_model == 2 || params->immun_model == 3) &&
		    num_act_cells-killed > 0)
		    killed_s = gsl_ran_poisson (params->ur, 
			    (num_act_cells-killed)*params->dImmun_s*num_cd8s*params->dt);

		// check for cd8 grouping option (strains share cd8s based on HD)
		if (params->immun_model == 4 || params->immun_model == 5) {
		    int cd8_index;
		    if (params->immun_model == 4)
			cd8_index = 1+MIN(params->maxHammingDist/params->cd8GroupSize,the_strain->get_founder_dist()/params->cd8GroupSize);
		    else if (params->immun_model == 5)
			cd8_index = 1+MIN(params->tF/params->cd8GroupSize,the_strain->get_birthdate()/params->cd8GroupSize);

		    int cd8_cells = params->cd8_groups[compartment][cd8_index].cd8_cells;

		    if (params->immun_model == 5)
			killed = gsl_ran_poisson (params->ur, 
			    num_act_cells*(CD8s-cd8_cells)*dImmun*params->dt);

		    killed += gsl_ran_poisson (params->ur, 
			num_act_cells*cd8_cells*params->dImmun_s*params->dt);
		}
	    }

	    int natural = 0;
	    double dA_immun;
	    int left_over = 0;

	    if (num_act_cells > killed-killed_s)
		 left_over= num_act_cells-killed-killed_s;

	    if (params->dA_ic50 > 0)
	    {
		the_strain->set_death_rate(compartment,params->dA * (1 + params->dA_scale*num_act_cells/(num_act_cells+params->dA_ic50)));
	    }
	    dA_immun = dA_eff * the_strain->get_death_rate(compartment,params,time);

	    // if the latent cells transition to memory INSTEAD of dying...
	    if (params->res_model == 3)
	    {
		int l1_s;
		l1_s = gsl_ran_poisson (params->ur, 
		    tau1*left_over*params->dt);

		if (left_over > l1_s)
		    left_over-=l1_s;
		else
		{
		    l1_s = left_over;
		    left_over=0;
		}

		//#latent L1 transitions from active 
		for (int l1=0; l1 < l1_s; l1++)
		{
		    the_strain->dec_act_cells(compartment,time);
		    the_strain->inc_l1_cells(compartment);
		    if (params->immun_model == 4 || params->immun_model == 5)
			params->cd8_groups[compartment][cd8_index].inf_cells--;
		    A_deaths++;
		    L1_infs++;
		}

		int l2_s;
		l2_s = gsl_ran_poisson (params->ur, 
		    tau2*left_over*params->dt);

		if (left_over > l2_s)
		    left_over-=l2_s;
		else
		{
		    l2_s = left_over;
		    left_over=0;
		}

		//#latent L2 transitions from active 
		for (int l2=0; l2 < l2_s; l2++)
		{
		    the_strain->dec_act_cells(compartment,time);
		    the_strain->inc_l2_cells(compartment);
		    if (params->immun_model == 4 || params->immun_model == 5)
			params->cd8_groups[compartment][cd8_index].inf_cells--;
		    A_deaths++;
		    L2_infs++;
		}
	    }
	    natural = gsl_ran_poisson (params->ur, 
		dA_immun*left_over*params->dt);

	    if (natural > left_over)
		natural=left_over;

	    int ad_s = killed + killed_s + natural;

	    if (ad_s > num_act_cells)
		ad_s = num_act_cells;

	    // Active deaths
	    for (int a=0; a < ad_s; a++)
	    {
		the_strain->dec_act_cells(compartment,time);
		the_strain->dec_assoc(params,ql,time);
		if (params->immun_model == 4 || params->immun_model == 5)
		    dec_cd8_inf(params,compartment,cd8_index,time);

		A_deaths++;

		if (A_deaths > A)
		    fprintf(stderr, 
			"Warning deaths exceed active cell count at t=%3.2lf\n",
			time);
	    }   


	}
	// Short-lived latently infected cell births, deaths and activations
	//
	if (num_l1_cells > 0)
	{
	    //Latent infection cell type 1
	    double L1_birth_rate = params->aL1*num_l1_cells;    //#homeostatic division proliferation

	    int num_births=gsl_ran_poisson(params->ur,L1_birth_rate*params->dt);
	    L1_births+=num_births;

	    //#latent L1 proliferations
	    for (int l=0; l < num_births; l++)
	    {
		the_strain->inc_l1_cells(compartment);
		int infected = the_strain->add_prolif(time);

		if (infected > params->max_strain_cells)
		    params->max_strain_cells=infected;
	    }
	    int l1d_s = gsl_ran_poisson (params->ur, 
			dL1*num_l1_cells*params->dt);

	    if (l1d_s > num_l1_cells)
		l1d_s = num_l1_cells;

	    // Latent L1 deaths
	    for (int l=0; l < l1d_s; l++)
	    {
		the_strain->dec_l1_cells(compartment);
		the_strain->dec_assoc(params,ql,time);

		L1_deaths++;
	    }   
	    num_l1_cells-=l1d_s;


	    double L1_activate_rate = params->xi1*num_l1_cells; //#latent to active
	    int num_acts=gsl_ran_poisson(params->ur,L1_activate_rate*params->dt);

	    if (num_acts > num_l1_cells)
		num_acts = num_l1_cells;

	    if (params->eclipse == 0)
	    {
		L1_acts +=num_acts;

		//#latent L1 activations 
		for (int l=0; l < num_acts; l++)
		{
		    the_strain->dec_l1_cells(compartment);
		    the_strain->inc_act_cells(compartment,time);
		    if (params->immun_model == 4 || params->immun_model == 5)
			params->cd8_groups[compartment][cd8_index].inf_cells++;
		}
	    }
	    else
	    {
		E1_acts +=num_acts;

		//#latent L1 activations 
		for (int l=0; l < num_acts; l++)
		{
		    the_strain->dec_l1_cells(compartment);
		    the_strain->inc_ecl_cells(compartment);
		    if (params->immun_model == 4 || params->immun_model == 5)
			params->cd8_groups[compartment][cd8_index].inf_cells++;
		}
	    }
	}


	// Long-lived latently infected cell births, deaths and activations
	//
	if (num_l2_cells > 0)
	{
	    double L2_birth_rate = params->aL2*num_l2_cells;    //#homeostatic division proliferation
	    int num_births=gsl_ran_poisson(params->ur,L2_birth_rate*params->dt);

	    //#latent proliferations
	    for (int l=0; l < num_births; l++)
	    {
		the_strain->inc_l2_cells(compartment);

		int infected = the_strain->add_prolif(time);
		if (infected > params->max_strain_cells)
		    params->max_strain_cells=infected;
	    }

	    L2_births+=num_births;


	    int l2d_s = gsl_ran_poisson (params->ur, 
			dL2*num_l2_cells*params->dt);

	    if (l2d_s > num_l2_cells)
		l2d_s = num_l2_cells;

	    // Latent L2 deaths
	    for (int l=0; l < l2d_s; l++)
	    {
		the_strain->dec_l2_cells(compartment);
		the_strain->dec_assoc(params,ql,time);

		L2_deaths++;
	    }   
	    num_l2_cells-=l2d_s;

	    double L2_activate_rate = params->xi2*num_l2_cells; //#latent to active
	    int num_acts=gsl_ran_poisson(params->ur,L2_activate_rate*params->dt);
	    if (num_acts > num_l2_cells)
		num_acts = num_l2_cells;

	    if (params->eclipse == 0)
	    {
		L2_acts +=num_acts;

		//#activations 
		for (int l=0; l < num_acts; l++)
		{
		    the_strain->dec_l2_cells(compartment);
		    the_strain->inc_act_cells(compartment,time);
		    if (params->immun_model == 4 || params->immun_model == 5)
			params->cd8_groups[compartment][cd8_index].inf_cells++;
		}
	    }
	    else
	    {
		E2_acts +=num_acts;

		//#activations 
		for (int l=0; l < num_acts; l++)
		{
		    the_strain->dec_l2_cells(compartment);
		    the_strain->inc_ecl_cells(compartment);
		    if (params->immun_model == 4 || params->immun_model == 5)
			params->cd8_groups[compartment][cd8_index].inf_cells++;
		}
	    }
	}

	if (num_ecl_cells > 0)
	{
	    //eclipse cell ripening
	    double ripen_rate = params->eclipse*num_ecl_cells; 

	    int num_ripens=gsl_ran_poisson(params->ur,ripen_rate*params->dt);

	    if (num_ripens > num_ecl_cells)
		num_ripens = num_ecl_cells;

	    A_ripens+=num_ripens;

	    for (int l=0; l < num_ripens; l++)
	    {
		the_strain->dec_ecl_cells(compartment);
		the_strain->inc_act_cells(compartment,time);
	    }
	}

	// New infections (junk, latent or active)
	if (num_virions > 0 && this_beta > 0 && ART_eff != 1)
	{
	    double A_infect_rate;
	    if (params->res_model != 3)
	    {
		double L1_infect_rate = tau1*S*this_beta*(1-ART_eff)*params->dt;
		int l1_s;
		//l1_s = gsl_ran_poisson(params->ur,L1_infect_rate*num_virions);

		// protect against binomial blowups!
		if (num_virions < 1e9)
		    l1_s = gsl_ran_binomial (params->ur, L1_infect_rate, num_virions);
		else
		    l1_s = L1_infect_rate*num_virions + gsl_ran_gaussian (params->ur, 
			    sqrt((L1_infect_rate*num_virions)*(1-(L1_infect_rate))));

		for (int l1=0; l1 < l1_s; l1++)
		{
		    process_infection(params, time, false, compartment, LATENT1,
			&CD8s,ql, the_strain,&Q);
		    L1_infs++;
		}

		double L2_infect_rate = tau2*S*this_beta*(1-ART_eff)*params->dt;
		int l2_s;
		//int l2_s=gsl_ran_poisson(params->ur,L2_infect_rate*num_virions);

		// protect against binomial blowups!
		if (num_virions < 1e9)
		    l2_s = gsl_ran_binomial (params->ur, L2_infect_rate, num_virions);
		else
		    l2_s = L2_infect_rate*num_virions + gsl_ran_gaussian (params->ur, 
			    sqrt((L2_infect_rate*num_virions)*(1-(L2_infect_rate))));

		for (int l2=0; l2 < l2_s; l2++)
		{
		    process_infection(params, time, false, compartment, LATENT2,
			&CD8s,ql, the_strain,&Q);
		    L2_infs++;
		}

		A_infect_rate = (1-tau1-tau2)*S*this_beta*(1-ART_eff)*params->dt;
	    } else {
		A_infect_rate = S*this_beta*(1-ART_eff)*params->dt;
	    }
	    
	    int a_s;
	    //int a_s=gsl_ran_poisson(params->ur,A_infect_rate*num_virions);

	    // protect against binomial blowups!
	    if (num_virions < 1e9)
		a_s = gsl_ran_binomial (params->ur, A_infect_rate, num_virions);
	    else
		a_s = A_infect_rate*num_virions + gsl_ran_gaussian (params->ur, 
			sqrt((A_infect_rate*num_virions)*(1-(A_infect_rate))));

	    for (int a=0; a < a_s; a++)
	    {
		process_infection(params, time, false,
		    compartment, ACTIVE, &CD8s,ql, the_strain,&Q);

		if (params->eclipse == 0)
		    A_infs++;
		else
		    E_infs++;
	    }
	}
	if (num_junk_virions > 0 && this_beta > 0 && ART_eff != 1)
	{
	    double A_infect_rate;
	    if (params->res_model != 3)
	    {
		double L1_infect_rate = tau1*S*this_beta*(1-ART_eff)*params->dt;
		int l1_s;
		//int l1_s=gsl_ran_poisson(params->ur,L1_infect_rate*num_junk_virions);

		// protect against binomial blowups!
		if (num_junk_virions < 1e9)
		    l1_s = gsl_ran_binomial (params->ur, L1_infect_rate, num_junk_virions);
		else
		    l1_s = L1_infect_rate*num_junk_virions + gsl_ran_gaussian (params->ur, 
			    sqrt((L1_infect_rate*num_junk_virions)*(1-(L1_infect_rate))));

		for (int l1=0; l1 < l1_s; l1++)
		{
		    process_infection(params, time, true, compartment, LATENT1,
			&CD8s,ql, the_strain,&Q);
		    junk_L1_infs++;
		}

		double L2_infect_rate = tau2*S*this_beta*(1-ART_eff)*params->dt;
		int l2_s;
		//int l2_s=gsl_ran_poisson(params->ur,L2_infect_rate*num_junk_virions);

		// protect against binomial blowups!
		if (num_junk_virions < 1e9)
		    l2_s = gsl_ran_binomial (params->ur, L2_infect_rate, num_junk_virions);
		else
		    l2_s = L2_infect_rate*num_junk_virions + gsl_ran_gaussian (params->ur, 
			    sqrt((L2_infect_rate*num_junk_virions)*(1-(L2_infect_rate))));

		for (int l2=0; l2 < l2_s; l2++)
		{
		    process_infection(params, time, true, compartment, LATENT2,
			&CD8s,ql, the_strain,&Q);
		    junk_L2_infs++;
		}

		A_infect_rate = (1-tau1-tau2)*S*this_beta*(1-ART_eff)*params->dt;
	    } else {
		A_infect_rate = S*this_beta*(1-ART_eff)*params->dt;
	    }
	    int a_s;
	    //int a_s=gsl_ran_poisson(params->ur,A_infect_rate*num_junk_virions);

	    // protect against binomial blowups!
	    if (num_junk_virions < 1e9)
		a_s = gsl_ran_binomial (params->ur, A_infect_rate, num_junk_virions);
	    else
		a_s = A_infect_rate*num_junk_virions + gsl_ran_gaussian (params->ur, 
			sqrt((A_infect_rate*num_junk_virions)*(1-(A_infect_rate))));

	    for (int a=0; a < a_s; a++)
	    {
		process_infection(params, time, true, compartment, ACTIVE,
		    &CD8s,ql, the_strain,&Q);

		junk_L0_infs++;
	    }
	}

	// CD8 T-cells
	//
	// if done strain-wise, grow populations for each strain
	if (params->use_cd8s > 0)
	{
	    if (params->immun_model == 2 || params->immun_model == 3) {
		double FIp = num_act_cells/(num_act_cells+this_ic50);

		unsigned long int min_cd8s = num_cd8s;
		if (min_cd8s==0)
		    min_cd8s=1;

		int cd8_b = gsl_ran_poisson (params->ur,params->tbirth*params->dt);
		int cd8_exp = gsl_ran_poisson (params->ur,params->theta*FIp*min_cd8s*params->dt);

		unsigned long int cd8_d = 0;

		if (params->immun_model == 3)
		    cd8_d = gsl_ran_poisson (params->ur,(params->delta+params->kappa*(CD8s/(CD8s+params->cd8_ic50)))*num_cd8s*params->dt);
		else
		    cd8_d = gsl_ran_poisson (params->ur,params->delta*num_cd8s*params->dt);
		if (cd8_d > num_cd8s)
		    cd8_d = num_cd8s;

		the_strain->inc_cd8s(compartment,cd8_b+cd8_exp);

		the_strain->dec_cd8s(compartment,(int)cd8_d);

		if (params->immun_model == 3)
		{
		    CD8_births+=cd8_b+cd8_exp;
		    CD8_deaths+=cd8_d;
		}

	    }
	}


	if ( num_tot_virions > 0 && (params->plotTopStrains ||
		(params->plotFirstStrains && num_act_cells > 0 &&
		    log10(num_act_cells) > params->cellThreshold)))
	{
	    int strain_index = the_strain->get_strain_index();

	    // first look for an empty slot
	    for (int t=0; t < MAX_STRAIN_COLORS; t++)
	    {
		if(params->plotFirstStrains)
		{
		    strain *plottedStrain = NULL; 
		    if (top_strains->get_num_elems() >= t+1)
			plottedStrain=top_strains->get_elem(t);

		    // new candidate?
		    if (!plottedStrain)
		    {
			top_strains->append(the_strain);
			break;
		    }
		    // already tracked?
		    else if (plottedStrain != NULL &&
			     plottedStrain->get_strain_index() == strain_index)
			break;
		}
		else if(params->plotTopStrains)
		{
		    strain *plottedStrain = NULL; 
		    if (top_strains->get_num_elems() >= t+1)
			plottedStrain=top_strains->get_elem(t);
		    
		    if (!plottedStrain)
		    {
			top_strains->append(the_strain);
			break;
		    }
		    // already tracked?
		    else if (plottedStrain != NULL &&
			     plottedStrain->get_strain_index() == strain_index)
			break;
		    else if ( num_tot_virions > 
			(plottedStrain->get_via_virus(compartment) + 
			 plottedStrain->get_junk_virus(compartment)))
		    {
			int curr_length = top_strains->get_num_elems();

			for (int t2=curr_length; t2 > t; t2--)
			{
			    if (t2==curr_length && curr_length < MAX_STRAIN_COLORS)
				top_strains->append(top_strains->get_elem(curr_length-1));
			    else if (t2 < MAX_STRAIN_COLORS)
				top_strains->set_elem(t2,top_strains->get_elem(t2-1));
			}
			top_strains->set_elem(t,the_strain);
			break;
		    }
		}
	    }
	}
	// track top N strains each frame then add to list if any are new
	if (num_tot_virions > 0 && params->followTopStrains)
	{
	    if (top_this_frame->get_num_elems() < params->maxTopStrains)
		top_this_frame->append(the_strain);
	    else
	    {
		// find the lowest entry and replace if less than current
		int lowest = -1;
		unsigned int min_v = 1e9;
		for (int i=0; i < top_this_frame->get_num_elems(); i++)
		{
		    strain *test_strain = top_this_frame->get_elem(i);
		    
		    unsigned int num_tst_virions = 
			(test_strain->get_via_virus(compartment) + 
			 test_strain->get_junk_virus(compartment));
		    if (num_tst_virions < min_v)
		    {
			lowest = i;
			min_v = num_tst_virions;
		    }
		}
		if (num_tot_virions > min_v)
		    top_this_frame->set_elem(lowest,the_strain);
	    }
	}
	// now update strain-wise and global virus counts
	if (num_virions > 0)
	{
	    int vd_s;
		vd_s = gsl_ran_poisson(params->ur,params->gam*num_virions*params->dt);

	    if ((unsigned int)vd_s > num_virions)
		vd_s = num_virions;

	    the_strain->dec_via_virus(compartment,vd_s);

	    num_virions -= vd_s;

	    Vvia_deaths +=vd_s;

	}
	if (num_junk_virions > 0)
	{
	    int vd_s;
		vd_s = gsl_ran_poisson(params->ur,params->gam*num_junk_virions*params->dt);

	    if ((unsigned int)vd_s > num_junk_virions)
		vd_s = num_junk_virions;

	    the_strain->dec_junk_virus(compartment,vd_s);

	    num_junk_virions -= vd_s;

	    Vnon_deaths +=vd_s;

	}
    }
    // if CD8s are managed globally only,  grow/cull single population here
    if (params->use_cd8s > 0)
    { 
	int cd8_b = 0;
	int cd8_exp = 0;
	unsigned long int cd8_d = 0;
	if (params->immun_model == 1 || params->immun_model == 6) 
	{
	    double FIp = A/(A+params->dImmun_IC50);

	    cd8_b = gsl_ran_poisson (params->ur,params->tbirth*params->dt);

	    cd8_exp = gsl_ran_poisson (params->ur,params->theta*FIp*CD8s*params->dt);

	    cd8_d = gsl_ran_poisson (params->ur,(params->delta+params->kappa*(CD8s/(CD8s+params->cd8_ic50)))*CD8s*params->dt);
	    //int cd8_d = gsl_ran_poisson (params->ur,params->delta*num_cd8s*params->dt);
	    if (cd8_d > CD8s-1)
		cd8_d = CD8s-1;

	    CD8_births+=cd8_b+cd8_exp;

	    CD8_deaths+=cd8_d;

	// check for cd8 grouping option (strains share cd8s based on HD or resistant mutation #)
	} else if (params->immun_model == 4 || params->immun_model == 5) {
	    int max_index;
	    if (params->immun_model == 4)
		max_index = params->maxHammingDist/params->cd8GroupSize;
	    else if (params->immun_model == 5)
		max_index = params->tF/params->cd8GroupSize;
	    for (int cd8_index = 0; cd8_index < max_index; cd8_index++)
	    {
		if (cd8_index == 0 && !params->kill_junk)
		    continue;

		int inf_cells = params->cd8_groups[compartment][cd8_index].inf_cells;
		unsigned long int cd8_cells = params->cd8_groups[compartment][cd8_index].cd8_cells;

		double FIp = inf_cells/(inf_cells+params->dImmun_IC50);

		unsigned long int min_cd8s = cd8_cells;
		//if (min_cd8s==0)
		    //min_cd8s=1;

		if (inf_cells > 0)
		    cd8_b = gsl_ran_poisson (params->ur,params->tbirth*params->dt);
		else
		    cd8_b = 0;

		cd8_exp = gsl_ran_poisson (params->ur,params->theta*FIp*min_cd8s*params->dt);

		cd8_d = gsl_ran_poisson (params->ur,(params->delta+params->kappa*(CD8s/(CD8s+params->cd8_ic50)))*cd8_cells*params->dt);
		//int cd8_d = gsl_ran_poisson (params->ur,params->delta*num_cd8s*params->dt);
		if (cd8_d > cd8_cells)
		    cd8_d = cd8_cells;

		params->cd8_groups[compartment][cd8_index].cd8_cells+=(cd8_b+cd8_exp-cd8_d);

		if (cd8_b +cd8_exp> 0)
		    CD8_births+=cd8_b+cd8_exp;
		if (cd8_d > 0)
		    CD8_deaths+=cd8_d;
	    }

	}
    }

    if (params->followTopStrains &&
		top_strains->get_num_elems() < MAX_FOLLOW_STRAINS)
    {
	// for each of the top strains from this frame
	// add to overall top strains list if not already there
	for (int t=0; t < top_this_frame->get_num_elems(); t++)
	{
	    strain *thisStrain=top_this_frame->get_elem(t);
	    int strain_index = thisStrain->get_strain_index();

	    bool found = false;
	    for (int t2=0; t2 < top_strains->get_num_elems(); t2++)
	    {
		strain *plottedStrain=top_strains->get_elem(t2);

		// already followed?
		if (plottedStrain->get_strain_index() == strain_index)
		{
		    found = true;
		    break;
		}
	    }
	    if (!found)
	    {
		top_strains->append(thisStrain);
		fprintf(stderr,"Following strain %d from t=%lf\n",
			strain_index,time);
	    }
	}
    }
    if (compartment == 0)
	params->numTopStrains = top_strains->get_num_elems();

    if (params->followTopStrains)
	delete(top_this_frame);

    //Susceptible and junk DNA cells
    //#constant production 
    //double Sbirth_rate = params->aS * params->S0;

    double Sbirth_rate;
    Sbirth_rate = params->aS;

    int Sbirths = gsl_ran_poisson(params->ur,Sbirth_rate*params->dt);

    double Sexp_rate = params->aV*S*(Vtot/(Vtot+params->aV_ic50));
    int Sexp = gsl_ran_poisson(params->ur,Sexp_rate*params->dt);

    //#density dependent susceptible death
    double Sdeath_rate = params->dS*S;
    int Sdeaths = gsl_ran_poisson(params->ur,Sdeath_rate*params->dt);

    if (Sdeaths > S)
	Sdeaths = S;


    int A_check=0;
    int A_new=0;
    int L1_check=0;
    int L1_new=0;
    int L2_check=0;
    int L2_new=0;

    for (int s=0; s < ql->get_num_elems(); s++)
    {
	if (s>=1 && s<=params->junk_bins+1)
	    continue; // skip the junk

	strain *the_strain = ql->get_elem(s);
	int num_act_cells = the_strain->get_act_cells(compartment);
	int num_l1_cells = the_strain->get_l1_cells(compartment);
	int num_l2_cells = the_strain->get_l2_cells(compartment);

	if (num_act_cells > 0) 
	    A_check+=num_act_cells;
	if (num_l1_cells > 0) 
	    L1_check+=num_l1_cells;
	if (num_l2_cells > 0) 
	    L2_check+=num_l2_cells;
    }
    A_new = A_infs - A_deaths 
	+ L1_acts + L2_acts+A_ripens;

    if (A_check != A+ A_new)
    {
	fprintf(stderr,
	    "Book keeping error (compartment %d) at t=%3.2lf (a_in=%d, a_act=%d)\n",
	    compartment,time,A_new,A_check);
	exit(1);
    }
    L1_new = L1_births + L1_infs - L1_deaths -E1_acts
		- L1_acts;

    if (L1_check != L1+L1_new)
    {
	fprintf(stderr,
	    "L1 Book keeping error (compartment %d) at t=%3.2lf (L1_in=%d, L1_act=%d)\n",
	    compartment,time,L1+L1_new, L1_check);
	exit(1);
    }
    L2_new =  L2_births + L2_infs - L2_deaths -E2_acts
		- L2_acts;

    if (L2_check != L2+ L2_new)
    {
	fprintf(stderr,
	    "L2 Book keeping error (compartment %d) at t=%3.2lf (L2_in=%d, L2_act=%d)\n",
	    compartment,time,L2+ L2_new,L2_check);
	exit(1);
    }
    if (L2 == 0 && (L2_births > 0 || (L2_acts > 0) || (L2_infs > 0 && ART_eff == 1)))
    {
	fprintf(stderr,
	    "L2 via growth error (compartment %d) at t=%3.2lf (L2=%d, L2_births=%d,L2_acts=%d, L2_infs=%d)\n",
	    compartment,time,L2 ,L2_births ,L2_acts ,L2_infs);
	exit(1);
    }

    // Watch for wrap-around on viable virus
    if (x[4] > 100000 && Vvia + Vvia_births - Vvia_deaths
	< 1000)
    {
	fprintf(stderr,
	    "V via growth error (compartment %d) at t=%3.2lf (Vvia was %u, now=%u)\n",
	    compartment,time,x[4],Vvia + Vvia_births - Vvia_deaths
		);
	exit(1);
    }


    // Watch for wrap-around on non-viable virus
    if (x[5] > 100000 && Vnon + Vnon_births - Vnon_deaths
		< 1000)
    {
	fprintf(stderr,
	    "V non-viable growth error (compartment %d) at t=%3.2lf (Vnon was %u, now=%u)\n",
	    compartment,time,x[5],Vnon + Vnon_births - Vnon_deaths
		);
	exit(1);
    }

    // adjust all compartment counts
    // Note: A_infs, etc. contain junk infections as well...
    int S_new = Sbirths + Sexp - Sdeaths 
	- L1_infs - L2_infs - A_infs;
    x[0]= MAX(0,S + S_new);

    x[1]= MAX(0,A + A_new);

    x[2]= MAX(0,L1 + L1_new);

    unsigned int x3p=x[3];

    x[3]= MAX(0,L2 + L2_new);

    if (x[3]==0 && x3p != 0)
	fprintf(stderr,"last memory inf at t=%lf\n",time);

    if (x[3]!=0 && x3p == 0)
	fprintf(stderr,"1st memory inf at t=%lf\n",time);

    int Vv_new = Vvia_births - Vvia_deaths;

    x[4]= MAX(0,Vvia + Vv_new);

    int Vn_new = Vnon_births - Vnon_deaths;

    x[5]= MAX(0,Vnon + Vn_new);

    x[6]= MAX(0,It + A_infs + L1_infs + L2_infs);

    x[7]= MAX(0,Id + A_deaths + L1_deaths + L2_deaths);

    x[8]= Q;

    if (params->use_cd8s > 0)
    {
	int CD8_new = CD8_births - CD8_deaths;
	if (params->CD8_limit > 0)
	    x[9]= MIN(CD8s + CD8_new,(unsigned long int)params->CD8_limit);
	else
	    x[9]= CD8s + CD8_new;
    }

    x[10]= MAX(0,junk_L0s + junk_L0_infs - junk_L0_deaths);

    x[11]= MAX(0,junk_L1s + junk_L1_infs + junk_L1_births - junk_L1_deaths);

    x[12]= MAX(0,junk_L2s + junk_L2_infs + junk_L2_births - junk_L2_deaths);

    x[13]= MAX(0,E + E_infs + E1_acts + E2_acts - A_ripens);

}

void read_input_file(char *inp_file, settings *params)
{
    FILE *inf = NULL;
    map<string,double> inputs;
    string par;
    double parv;

    char tmpline[MAX_LINE];
    char *valuep;
    int i=0;

    if ((inf = fopen (inp_file,"r")) == NULL) {
	cerr << "Could not open input file "<<inp_file<<"\n";
	exit(1);
    }
    while (fgets(tmpline, MAX_LINE-1,inf) != NULL) {
	i++;
	tmpline[MAX_LINE-1] = '\0';
	valuep = rindex(tmpline,' ');
	if (valuep == NULL) {
	    cerr << "Error while reading parameter name from "<<inp_file<<" at parameter #"<<i<<"\n";
	    exit(1);
	}
	*valuep = '\0';
	par = tmpline;
	
	if (sscanf(valuep+1,"%lf",&parv) != 1) {
	    cerr << "Error while reading value for "<<par<<" in "<<inp_file<<" (parameter #"<<i<<")\n";
	    exit(1);
	}
	else
	{
	    inputs[par] = parv;
	    cout <<"#"<< i <<" "<<par <<" = "<< inputs[par];
	}
	cout<<endl;
    }
    cout<<endl;
    cout <<"Finished reading input file "<<inp_file<<".\n";
    fclose (inf);

    CHECK_FOR_INT("verbose",verbose);
    CHECK_FOR_INT("Runs",Runs);
    CHECK_FOR_INT("AIC_Runs",AIC_Runs);
    CHECK_FOR_REAL("print",print);

    CHECK_FOR_REAL("ART_eff",ART_eff);
    CHECK_FOR_REAL("ART_start",ART_start);
    CHECK_FOR_REAL("ART_stop",ART_stop);

    CHECK_FOR_REAL("weight1",weight1);
    CHECK_FOR_REAL("weight2",weight2);
    CHECK_FOR_REAL("weight3",weight3);
    CHECK_FOR_REAL("weight4",weight4);
    CHECK_FOR_REAL("weight5",weight5);
    CHECK_FOR_REAL("weight6",weight6);
    CHECK_FOR_REAL("weight7",weight7);
    CHECK_FOR_REAL("weight8",weight8);

    CHECK_FOR_REAL("Input_refresh",Input_refresh);
    CHECK_FOR_REAL("Alt_input",Alt_input);

    CHECK_FOR_INT("use_cd8s",use_cd8s);

    CHECK_FOR_REAL("junk",junk);
    CHECK_FOR_INT("kill_junk",kill_junk);
    CHECK_FOR_INT("junk_bins",junk_bins);
    CHECK_FOR_INT("junk_days",junk_days);
    CHECK_FOR_INT("res_days",res_days);

    CHECK_FOR_REAL("eclipse",eclipse);

    CHECK_FOR_INT("a0",a0);
    CHECK_FOR_REAL("tF",tF);
    CHECK_FOR_REAL("dt",dt);

    CHECK_FOR_REAL("thL1",thL1);
    CHECK_FOR_REAL("aL1",aL1 );
    CHECK_FOR_REAL("xi1",xi1 );
    CHECK_FOR_REAL("tau1",tau1);
    CHECK_FOR_REAL("thL2",thL2);
    CHECK_FOR_REAL("aL2",aL2 );
    CHECK_FOR_REAL("xi2",xi2 );
    CHECK_FOR_REAL("tau2",tau2);
    CHECK_FOR_REAL("thJ1",thJ1);
    CHECK_FOR_REAL("thJ2",thJ2);

    CHECK_FOR_REAL("res_factor",res_factor);
    CHECK_FOR_INT("res_model",res_model);

    CHECK_FOR_REAL("dS",dS );
    CHECK_FOR_REAL("dA",dA );
    CHECK_FOR_REAL("dA_ic50_0",dA_ic50_0 );
    CHECK_FOR_REAL("dA_scale",dA_scale );
    CHECK_FOR_REAL("dActivated",dActivated );

    CHECK_FOR_INT("immun_model",immun_model );

    CHECK_FOR_REAL("dImmun",dImmun);
    CHECK_FOR_REAL("dImmun_s",dImmun_s);

    CHECK_FOR_REAL("cd8_k",cd8_k );
    CHECK_FOR_REAL("cd8_mean",cd8_mean );
    CHECK_FOR_INT("cd8_0",cd8_0);

    CHECK_FOR_REAL("cd8_ic50_0",cd8_ic50_0 );
    CHECK_FOR_REAL("dImmun_IC50_0",dImmun_IC50_0 );

    CHECK_FOR_REAL("tbirth0",tbirth0 );
    CHECK_FOR_REAL("delta",delta );
    CHECK_FOR_REAL("theta",theta );
    CHECK_FOR_REAL("kappa",kappa );

    CHECK_FOR_REAL("gam",gam);

    CHECK_FOR_REAL("mu0",mu0 );
    CHECK_FOR_REAL("transit_prob",transit_prob );
    CHECK_FOR_INT("use_k2p",use_k2p );
    CHECK_FOR_REAL("entropy_threshold",entropy_threshold); // entropy below which one can assume beta=0
    CHECK_FOR_REAL("entropy_factor",entropy_factor); // used to map entropy to fitness inheritance

    CHECK_FOR_REAL("pi",pi );
    //CHECK_FOR_REAL("Bt_founder",Bt_founder);
    CHECK_FOR_REAL("Bt_mean",Bt_mean);
    CHECK_FOR_REAL("Bt_sd_factor",Bt_sd_factor);
    CHECK_FOR_REAL("Bt_max",Bt_max);
    CHECK_FOR_REAL("Bt_k",Bt_k);
    CHECK_FOR_REAL("lifespan_decay",lifespan_decay);
    CHECK_FOR_INT("use_dBeta",use_dBeta);

    CHECK_FOR_REAL("ic50_k",ic50_k);
    CHECK_FOR_REAL("aS0",aS0);
    CHECK_FOR_REAL("aV",aV);
    CHECK_FOR_REAL("aV_ic50",aV_ic50);
    CHECK_FOR_REAL("vol",vol);
    CHECK_FOR_REAL("dna_for_res",dna_for_res);





    CHECK_FOR_INT("disp_patient",disp_patient);
    CHECK_FOR_INT("vl_as_log",vl_as_log);

    CHECK_FOR_INT("age_bins",age_bins);
    CHECK_FOR_INT("strain_bins",strain_bins);
    CHECK_FOR_REAL("sampleInterval",sampleInterval);
    CHECK_FOR_REAL("crit_start",crit_start);

    CHECK_FOR_INT("max_strains",max_strains);
    CHECK_FOR_INT("maxTopStrains",maxTopStrains);
    CHECK_FOR_INT("cellThreshold",cellThreshold);
    CHECK_FOR_INT("CD8_limit",CD8_limit);
    CHECK_FOR_INT("followTopStrains",followTopStrains);
    CHECK_FOR_INT("plotTopStrains",plotTopStrains);
    CHECK_FOR_INT("plotFirstStrains",plotFirstStrains);

    CHECK_FOR_REAL("seq_length",seq_length);
    CHECK_FOR_REAL("div_sample_interval",div_sample_interval);
    CHECK_FOR_INT("div_sample_size",div_sample_size);
    CHECK_FOR_INT("score_sample",score_sample);
    CHECK_FOR_INT("seq_sample_type",seq_sample_type);
    CHECK_FOR_INT("AIC_k",AIC_k);

    // used to decay latent cells when scoring fossil record
    CHECK_FOR_REAL("art_start_mean",art_start_mean);
    CHECK_FOR_REAL("art_start_std",art_start_std);
    CHECK_FOR_REAL("art_dur_mean",art_dur_mean);
    CHECK_FOR_REAL("art_dur_std",art_dur_std);

    if (params->max_log_v_value == 0 && params->max_vl > 0)
	params->max_log_v_value=(double)((int)(log10((double)params->max_vl)+0.5));

    if (params->junk_bins > MAX_JUNK_BINS)
	params->junk_bins = MAX_JUNK_BINS;

#ifdef GUI_ENABLED
    CHECK_FOR_INT("plotVt",plotVt);
    CHECK_FOR_INT("plotAct",plotAct);
    CHECK_FOR_INT("show_params",show_params);
    CHECK_FOR_INT("plotCD8s",plotCD8s);

    CHECK_FOR_INT("colorByHamming",colorByHamming);
    CHECK_FOR_INT("maxHammingDist",maxHammingDist);
    CHECK_FOR_INT("cd8GroupSize",cd8GroupSize);
    CHECK_FOR_INT("colorBoxes",colorBoxes);
    CHECK_FOR_INT("showNonVia",showNonVia);

    CHECK_FOR_REAL("time_bias",time_bias);
    CHECK_FOR_REAL("plot_span",plot_span);
    CHECK_FOR_REAL("plot_bias",plot_bias);

    CHECK_FOR_REAL("refresh",refresh);

    CHECK_FOR_INT("SnapshotInterval",SnapshotInterval);
    CHECK_FOR_INT("AutoSnapshot",AutoSnapshot);

    CHECK_FOR_INT("max_log_v_value",max_log_v_value);
    CHECK_FOR_INT("min_log_v_value",min_log_v_value);
    CHECK_FOR_INT("max_log_a_value",max_log_a_value);
    CHECK_FOR_INT("min_log_a_value",min_log_a_value);
    CHECK_FOR_INT("max_log_cd8_value",max_log_cd8_value);
    CHECK_FOR_INT("min_log_cd8_value",min_log_cd8_value);
    CHECK_FOR_INT("max_cd8s",max_cd8s);

    CHECK_FOR_INT("max_act",max_act);
    CHECK_FOR_INT("x_ticks",x_ticks);
    CHECK_FOR_INT("y1_ticks",y1_ticks);
    CHECK_FOR_INT("y2_ticks",y2_ticks);
    CHECK_FOR_INT("y3_ticks",y3_ticks);

    if (params->max_strains > MAX_RAINBOW_COLORS)
	params->max_strains = MAX_RAINBOW_COLORS;
#endif
}

void read_fasta_file(char *fasta_file, settings *params)
{
    FILE *inf = NULL;
    char tmpline[MAX_LINE];
    int linenum=0;

    params->founder_stats.length = 0;
    params->founder_stats.as = 0;
    params->founder_stats.cs = 0;
    params->founder_stats.gs = 0;
    params->founder_stats.ts = 0;
    params->founder_seq = "";

    if ((inf = fopen (fasta_file,"r")) == NULL) {
	cerr << "Could not open FASTA file "<<fasta_file<<"\n";
	exit(1);
    }
    if (fgets(tmpline, MAX_LINE-1,inf) == NULL || tmpline[0] != '>') {
	cerr << "FASTA file "<<fasta_file<<" appears invalid (no \">\" in line 1)\n";
	exit(1);
    }
	
    while (fgets(tmpline, MAX_LINE-1,inf) != NULL) {
	linenum++;
	for (int j=0; j < (int)strlen(tmpline);j++) {
	    if (tmpline[j]=='\n')
		tmpline[j]='\0';
	}
	params->founder_stats.length += strlen(tmpline);
	params->founder_seq += tmpline;

	for (int j=0; j < (int)strlen(tmpline);j++) {
	    if (tmpline[j] == 'A' || tmpline[j] == 'A')
		params->founder_stats.as++;
	    else if (tmpline[j] == 'C' || tmpline[j] == 'c')
		params->founder_stats.cs++;
	    else if (tmpline[j] == 'G' || tmpline[j] == 'g')
		params->founder_stats.gs++;
	    else if (tmpline[j] == 'T' || tmpline[j] == 't')
		params->founder_stats.ts++;
	    else {
		cerr << "FASTA file "<<fasta_file<<" appears invalid (no A/C/G/T in line "<<linenum<<"at pos"<<j<<")\n";
		exit(1);
	    }
	}
    }
    cout<<endl;
    cout <<"Finished reading input file "<<fasta_file<<".\n";
    cout <<"Read "<<params->founder_stats.length << " bases (" << params->founder_stats.as<<"/"<<params->founder_stats.cs<<"/"<<params->founder_stats.gs<<"/"<<params->founder_stats.ts<<").\n";
    fclose (inf);
    params->seq_length = params->founder_stats.length;
}
                              
void usage(char *prog_name)
{

    fprintf(stderr,
	"Usage: %s [-h][-l][-f <input_file>][-c <crit file>][-e <entropy file>][-p <pat>][-q <seq file>][-r][-s <seed>][-v][-w <write_mask>]\n",prog_name);
    fprintf(stderr,"\t-h = this help\n");
    fprintf(stderr,"\t-f = optional input file\n");
    fprintf(stderr,"\t-a = optional follow-on input file (at Input_refresh)\n");
    fprintf(stderr,"\t-c = optional criteria file\n");
    fprintf(stderr,"\t-e = optional entropy file for fitness inheritance (requires -q)\n");
    fprintf(stderr,"\t\tFormat: target mean and variance values for...\n");
    fprintf(stderr,"\t\t\tlog peak viral load\n");
    fprintf(stderr,"\t\t\ttime to peak viral load\n");
    fprintf(stderr,"\t\t\tlog peak-nadir drop in viral load\n");
    fprintf(stderr,"\t\t\ttime to nadir viral load\n");
    fprintf(stderr,"\t\t\tlog viral load set point\n");
    fprintf(stderr,"\t\t\tlog variance from viral load set point\n");
    fprintf(stderr,"\t\t(set point taken as measure 10 days after nadir)\n");
    fprintf(stderr,"\t\t(set point variance measured at 10 day intervals)\n");
    fprintf(stderr,"\t-l = viral loads in as log values (for -z)\n");
    fprintf(stderr,"\t-p <patient>= viral loads in as log values (for -z)\n");
    fprintf(stderr,"\t-y = optional sequence sampling schedule file\n");
    fprintf(stderr,"\t-z = optional target viral load file\n");
    fprintf(stderr,"\t\tcurrently it must contain N patients each with M time then Vt values for scored fit\n");
    fprintf(stderr,"\t\tFormat:\n");
    fprintf(stderr,"\t\t\t<N>\n");
    fprintf(stderr,"\t\t\t<M>\n");
    fprintf(stderr,"\t\t\tt1\n");
    fprintf(stderr,"\t\t\t...\n");
    fprintf(stderr,"\t\t\tt<m>\n");
    fprintf(stderr,"\t\t\vt1\n");
    fprintf(stderr,"\t\t\t...\n");
    fprintf(stderr,"\t\t\vt<m>\n");
    fprintf(stderr,"\t-q <seq FASTA>= FASTA file with initial sequence (for a/c/g/t ratios)\n");
    fprintf(stderr,"\t-r = random number seed truely random\n");
    fprintf(stderr,"\t-s <seed> = random number seed for repeating run (unsigned)\n");
    fprintf(stderr,"\t-w <write_mask> = which output (csv) files to generate (Bit-mask)\n");
    fprintf(stderr,"\t\t1= main compartment counts vs time\n");
    fprintf(stderr,"\t\t2= strains and their cell counts vs time\n");
    fprintf(stderr,"\t\t4= strain lineage information (birth, parent, etc.)\n");
    fprintf(stderr,"\t\t8= top 10 strain ids and viral loads\n");
}

double ScoreFunction(settings *params)
{
    double score = 0;
    double score1 = 0;
    double score2 = 0;
    double score3 = 0;
    double score4 = 0;
    double score5 = 0;
    double score6 = 0;
    double score7 = 0;
    double score8 = 0;
    double score9 = 0;

    int runnum = 0;
    int discarded = 0;

    double avg_peak_vl = 0;
    double log_peak_vl = 0;

    double avg_peak_time = 0;

    double log_nadir_vl = 0;

    double avg_nadir_drop = 0;
    double log_nadir_drop = 0;

    double avg_nadir_time = 0;

    double avg_set_pt = 0;

    double avg_set_pt_var = 0;

    double avg_dH_rate = 0;
    double dH_rate = 0;

    double avg_div_rate = 0;
    double div_rate = 0;

    double avg_junk_fract = 0;
    double junk_fract = 0;

    double time = 0;
    double diversity = 0;

    double blood_factor = params->vol / 1000;

    double set_pts[MAX_SET_POINTS];

    unsigned int x[NUM_MEASURES];
    int num_compartments = 1;


    params->mu = params->mu0 * params->seq_length; // mutation rate for sequence

    params->aS  = params->aS0*params->vol;       //# constant growth rate of susceptible cells [cells/day]
    params->Bt = params->Bt_mean/params->vol; //# infection rate of T-cells [1/virus-day]

    params->tbirth  = params->tbirth0*params->vol;  //# constant growth rate of CD8 cells [cells/day]
    params->dImmun_IC50  = params->dImmun_IC50_0*params->vol;  //# volume scaled infected cell ic50
    params->dA_ic50  = params->dA_ic50_0*params->vol;  //# volume scaled infected cell ic50
    params->cd8_ic50  = params->cd8_ic50_0*params->vol;  //# volume scaled cd8 carrying capacity


    //#calculated rates
    params->dL1 = params->aL1+params->thL1+params->xi1;           //# death rate
    params->dL2 = params->aL2+params->thL2+params->xi2;           //# death rate
    if (params->verbose)
	cout << "(latent death rates are "<< params->dL1 <<" and "<< params->dL2 <<")"<<endl;

    params->l1 = 1-(1+params->xi1/params->thL1)*params->tau1;    //# latency factor
    params->l2 = 1-(1+params->xi2/params->thL2)*params->tau2;    //# latency factor
    params->Ro = params->aS*(params->Bt_mean/params->vol)*params->pi*(1-params->junk)/params->gam/params->dS/params->dA; //# basic reproductive number

    if (params->verbose)
	cout << "(basic reproductive number is "<< params->Ro <<")"<<endl;
    //#equilibrium solutions
    double Seq=params->gam*params->dA/(params->Bt)/params->pi/(params->l1+params->l2);
    double L1eq=params->tau1/params->thL1*(params->gam*params->dS*params->dA/(params->Bt)/params->pi/params->l1-params->aS);
    double L2eq=params->tau2/params->thL2*(params->gam*params->dS*params->dA/(params->Bt)/params->pi/params->l2-params->aS);
    double Aeq=params->aS*(params->l1+params->l2)/params->dA-params->gam*params->dS/(params->Bt)/params->pi;
    double Veq=params->aS*params->pi*(params->l1+params->l2)/params->gam/params->dA-params->dS/(params->Bt);

    fprintf(stderr,"At equilibrium: S=%3.2lf, L1=%3.2lf, L2=%3.2lf, A=%3.2lf, log V=%3.2lf\n",
	Seq,L1eq,L2eq,Aeq,(Veq>0)?log10(Veq):0);


    while (runnum < params->Runs && !params->stopFlag)
    { 
	/* if we read in an alternate inpuit file last run, then re-read the original input file! */
	if (params->alt_inp_file != NULL && runnum > 0)
	{
	    fprintf(stderr,"Re-reading original input file!\n");

	    read_input_file(params->inp_file,params);
	}

	generic_list<strain *> *Q_list=new generic_list<strain *>();
	//fprintf(stderr,"Strain list is %x\n",(void *)Q0_list);

	params->S0=params->aS/params->dS;	//S
	x[0]=params->S0;
	x[1]=params->a0;		//A - viable
	x[2]=0;			//L1 - viable
	x[3]=0;			//L2 - viable
	x[4]=0;			//Vt - viable
	x[5]=0;			//Vt - non-viable
	x[6]=0;			//It
	x[7]=0;			//Id
	x[8]=2;			//Strains (start with junk and founder, maybe non-viable)
	x[9]=0;			//Total CD8s
	x[10]=0;			//Total L0 cells with junk DNA
	x[11]=0;			//Total L1 cells with junk DNA
	x[12]=0;			//Total L2 cells with junk DNA
	x[13]=0;			//Eclipse, (pre-productive)
	x[14]=0;			//No FDC holding in plasma
	x[15]=0;			//No FDC holding in plasma

	params->numStrains = 1; // founder
	params->youngQ = 0;

	strain *founder = new strain(params,0,0,NULL,0);
	Q_list->append(founder);

	if (params->use_cd8s > 0)
	{
	    if (params->immun_model == 4 || params->immun_model == 5) {
		//params->cd8_0=params->tbirth/params->delta;
		params->cd8_0=0;

		int max_index;
		if (params->immun_model == 4)
		    max_index = params->maxHammingDist/params->cd8GroupSize;
		else if (params->immun_model == 5)
		    max_index = params->tF/params->cd8GroupSize;

		for (int i=0; i < num_compartments; i++)
		    for (int j=0; j < max_index; j++) {
			params->cd8_groups[i][j].cd8_cells = 0;
			params->cd8_groups[i][j].inf_cells = 0;
		}
	    } 
	    else
		x[9] += founder->get_cd8_cells(0);
	}
	for (int i = 0; i < params->a0; i++)
	{
	    int infected = founder->add_infected(0);
	    founder->inc_act_cells(0,0);
	    if (infected > params->max_strain_cells)
		params->max_strain_cells=infected;
	    params->cd8_groups[0][1].inf_cells++;
	}
	// create junk J0 DNA strain bins starting at position 1
	for (int i = 1; i <= params->junk_bins+1; i++)
	{
	    strain *new_strain = new strain(params,i,0,founder,-1);
	    Q_list->append(new_strain);
	}



	time = 0;
	diversity = 0;

	double time_p = 0;
	double next_output = 0;
	double dImmun_avg = 0;
	double Bt_avg = 0;
	double Bt_max = 0;

	double peak_vl = 0;
	double peak_time = 0;
	double nadir_vl = 0;
	double nadir_time = 0;
	double last_set_time = 0;
	double set_pt_err = 0;
	double set_pt_p = 0;
	int num_set_pts = 0;
	bool this_discarded = false;

	int q_a = 0;
	int qmax = 0;
	int qmax_via = 0;
	int qmax_junk = 0;
	int qmax_a = 0;
	int qmax_l1 = 0;
	int qmax_l2 = 0;
	int qmax_dH = 0;
	unsigned long int qmax_cd8s = 0;
	int qmax_cd8_index = 0;
	double qmax_dA = 0;
	double qmax_beta = 0;
	double qmax_dImmun = 0;
	int qmax_ns = 0;
	double entropy = 0;
	double seq_ident = 0;

	double samp_avgHd = 0;
	int samp_maxHd = 0;

	double real_avgHd = 0;
	int real_maxHd = 0;

	double art_start = 0;
	double art_dur = 0;

	unsigned int *strain_bins=NULL;
	unsigned int *age_bins=NULL;
	struct seq_info *div_samples=NULL;
	struct seq_info *seq_samples=NULL;


#ifdef GUI_ENABLED
	int snapnum = 0;

	bool snap_this_frame = false;
#endif

	generic_list<strain *> *top_strains=new generic_list<strain *>();

	if (params->verbose)
	    fprintf(stderr,"starting run %d...\n",runnum+1);

	if (params->art_start_mean != 0 )
	{
	    // draw from distribution, but limit to at least 30 days
	    art_start = params->art_start_mean+gsl_ran_gaussian (params->ur, params->art_start_std);
	    if (art_start < 30)
		art_start = 30;
	    if (art_start > params->tF)
		art_start = params->tF-1;
	}

	if (params->art_dur_mean != 0 )
	{
	    // draw from distribution, but limit to at least 8 months (240 days)
	    art_dur = params->art_dur_mean+gsl_ran_gaussian (params->ur, params->art_dur_std);
	    if (art_dur < 240)
		art_dur = 240;
	}

	double start_l1=0;
	double start_l2=0;
	double start_junk_l1=0;
	double start_junk_l2=0;

	if (params->writeOn )
	{
	    strain_bins = (unsigned int *)my_malloc(params->strain_bins*sizeof (unsigned int *));
	    age_bins = (unsigned int *)my_malloc(params->age_bins*sizeof (unsigned int *));
	    
	    // state outputs
	    if (params->dataF1 != NULL) 
	    {
		fprintf(params->dataF1,"time");
		fprintf(params->dataF1,",susceptible");
		fprintf(params->dataF1,",active cells");
		fprintf(params->dataF1,",L1 cells");
		fprintf(params->dataF1,",L2 cells");
		fprintf(params->dataF1,",log V viable");
		fprintf(params->dataF1,",log V non-viable");
		fprintf(params->dataF1,",log Vt");
		fprintf(params->dataF1,",log It");
		fprintf(params->dataF1,",log Id");
		fprintf(params->dataF1,",strains");
		fprintf(params->dataF1,",Bt equiv");
		fprintf(params->dataF1,",Bt max");
		fprintf(params->dataF1,",Ro equiv");
		fprintf(params->dataF1,",act strains ");
		fprintf(params->dataF1,",strain with max acts");
		fprintf(params->dataF1,",max strain act cells");
		fprintf(params->dataF1,",max strain l1 cells");
		fprintf(params->dataF1,",max strain l2 cells");
		fprintf(params->dataF1,",max strain dH");
		fprintf(params->dataF1,",max strain cd8s");
		fprintf(params->dataF1,",max strain ic50");
		fprintf(params->dataF1,",max strain beta");
		fprintf(params->dataF1,",strain avg dH");
		fprintf(params->dataF1,",strain max dH");
		fprintf(params->dataF1,",strain entropy");
		fprintf(params->dataF1,",strain diversity");
		fprintf(params->dataF1,",sequence ident");
		fprintf(params->dataF1,",total cd8s");
		fprintf(params->dataF1,",junk L0 DNA cells");
		fprintf(params->dataF1,",junk L1 DNA cells");
		fprintf(params->dataF1,",junk L2 DNA cells");
		fprintf(params->dataF1,",junk DNA strains");
		fprintf(params->dataF1,",junk DNA founder dist");
		fprintf(params->dataF1,",junk DNA from prolif");
		fprintf(params->dataF1,",pre-prod viable cells");
		fprintf(params->dataF1,",pre-prod non-viable cells");
		fprintf(params->dataF1,",max strain dImmun");
		fprintf(params->dataF1,",avg strain dImmun");
		fprintf(params->dataF1,"\n");
		fflush(params->dataF1);
	    }
	    // strain output
	    if (params->dataF2 != NULL) 
	    {
		fprintf(params->dataF2,"time");
		fprintf(params->dataF2,",strains");
		fprintf(params->dataF2,",bin_size");
		for (int j=0; j < params->strain_bins; j++)
		    fprintf(params->dataF2,",bin[%d]",j+1);
		fprintf(params->dataF2,"\n");
		fflush(params->dataF2);
	    }
	    // strain dump (at end of run)
	    if (params->dataF4 != NULL) 
	    {
		fprintf(params->dataF4,"strain");
		fprintf(params->dataF4,",strain birth");
		fprintf(params->dataF4,",founder dist");
		fprintf(params->dataF4,",strain parent");
		fprintf(params->dataF4,",strain beta");
		fprintf(params->dataF4,",cells at end");
		fprintf(params->dataF4,",total cells");
		fprintf(params->dataF4,",last event time");
		fprintf(params->dataF4,",total infs");
		fprintf(params->dataF4,",total prolifs");
		fprintf(params->dataF4,",die out time");
		fprintf(params->dataF4,",viable virions at end");
		fprintf(params->dataF4,",junk virions at end");
		fprintf(params->dataF4,",cd8s at end");
		fprintf(params->dataF4,",highest cd8 count");
		fprintf(params->dataF4,",strain as");
		fprintf(params->dataF4,",strain cs");
		fprintf(params->dataF4,",strain gs");
		fprintf(params->dataF4,",strain ts");
		fprintf(params->dataF4,",change");
		fprintf(params->dataF4,",non syns");
		fprintf(params->dataF4,",non syn index");
		fprintf(params->dataF4,"\n");
		fflush(params->dataF4);
	    }
	    // strain time history
	    if (params->dataF3 != NULL) 
	    {
		fprintf(params->dataF3,"time");
		fprintf(params->dataF3,",active strains");
		fprintf(params->dataF3,",log Viable Vt");
		fprintf(params->dataF3,",log Junk Vt");
		fprintf(params->dataF3,",track top");
		fprintf(params->dataF3,",num tracked");
		for (int j=0; j < MAX_STRAIN_COLORS; j++)
		{
		    fprintf(params->dataF3,",id[%d]",j+1);
		    fprintf(params->dataF3,",As[%d]",j+1);
		    fprintf(params->dataF3,",log via Vt_s[%d]",j+1);
		    fprintf(params->dataF3,",log junk Vt_s[%d]",j+1);
		    fprintf(params->dataF3,",cd8_s[%d]",j+1);
		}
		fprintf(params->dataF3,"\n");
		fflush(params->dataF3);
	    }
	    // followed strains time history
	    if (params->followTopStrains && params->dataF7 != NULL) 
	    {
		fprintf(params->dataF7,"time");
		fprintf(params->dataF7,",followed strains");
		for (int j=0; j < MAX_FOLLOW_STRAINS; j++){
		    fprintf(params->dataF7,",V%d",j+1);
		}
		fprintf(params->dataF7,"\n");
		fflush(params->dataF7);
	    }
	}
			    
	//#loop over times
	//bool lastFrame=false;
	double next_print = 0;
	double next_res_print = params->junk_days;
	params->numTopStrains=0;
	params->NextRefresh = 0;
	params->time = 0;

	double next_div_samples;
	double next_seq_samples;

	int num_seq_samples=0;
	int num_div_samples=0;
	int seq_sampling_index = 0;
	int seq_sample_size = 0;
	int div_sample_size = 0;

	if (params->num_seq_samplings > 0)
	{
	    next_seq_samples = params->seq_sampling[seq_sampling_index].time;
	    seq_sample_size = params->seq_sampling[seq_sampling_index].seqs;
	}
	next_div_samples = params->div_sample_interval;
	div_sample_size = params->div_sample_size;



	while (time < params->tF && !params->stopFlag)
	{

	    if (params->alt_inp_file != NULL)
	    {
		if(params->Alt_input > 0 && time>=params->Alt_input &&
		time<=params->Alt_input+params->dt)
		{
		    int old_a0=params->a0;

		    fprintf(stderr,"Reading in %s at t=%lf\n",
			    params->alt_inp_file,time);

		    read_input_file(params->alt_inp_file,params);

		    if (old_a0 == 0 && params->a0 > 0)
			for (int i = 0; i < params->a0; i++)
			{
			    int infected = founder->add_infected(0);
			    founder->inc_act_cells(0,0);
			    if (infected > params->max_strain_cells)
				params->max_strain_cells=infected;
			}
		}
		if(params->Input_refresh > 0 && time>=params->Input_refresh &&
		    time<=params->Input_refresh+params->dt)
		{
		    fprintf(stderr,"Re-reading original input file!\n");

		    read_input_file(params->inp_file,params);
		}

	    }
	    if (time >= next_div_samples &&
		num_div_samples < params->score_sample)
	    {
		int act_samples = sample_seqs_for_div(params,Q_list,0,div_sample_size,&div_samples);
		if (act_samples > 0)
		{
		    fprintf(stdout,"%d sampled sequences at t=%lf...\n",act_samples,time);
		    fprintf(stdout,"Seqs: ");
		    for (int i=0;i < act_samples-1; i++)
			    fprintf(stdout,"%d,",div_samples[i].seq_index);

		    fprintf(stdout,"%d\n",div_samples[act_samples-1].seq_index);

		    gather_sample_stats(params,Q_list, act_samples,div_samples,
			&samp_avgHd,&samp_maxHd,&entropy,&diversity,
			&params->ttca_mean,&seq_ident);
		    fprintf(stdout,
			"samp_avgHd=%lf,samp_maxHd=%d,diversity=%lf,seq_ident=%lf,entropy=%lf\n",
			samp_avgHd,samp_maxHd,diversity,seq_ident,entropy);

		    fprintf(stdout,
			"ttca_mean=%lf\n",params->ttca_mean);
		    free(div_samples);
		}
		next_div_samples = time+params->div_sample_interval;
		num_div_samples++;
	    }
	    if (time >= next_seq_samples && 
		num_seq_samples < params->num_seq_samplings)
	    {
		int act_samples = new_sample_seqs(params,Q_list,0,seq_sample_size,params->div_sample_size,&seq_samples);
		if (act_samples > 0)
		{
		    fprintf(stdout,"%d sampled sequences at t=%lf...\n",act_samples,time);
		    fprintf(stdout,"Seqs: ");
		    for (int i=0;i < act_samples-1; i++)
			fprintf(stdout,"%d,",seq_samples[i].seq_index);


		    fprintf(stdout,"%d\n",seq_samples[act_samples-1].seq_index);
		    free(seq_samples);
		}
		else
		    fprintf(stdout,"No Active seq samples at =%lf\n",time);

		fflush(stdout);

		seq_sampling_index++;
		next_seq_samples = params->seq_sampling[seq_sampling_index].time;
		seq_sample_size = params->seq_sampling[seq_sampling_index].seqs;
		num_seq_samples++;
	    }
	    if (art_start != 0)
	    {
		if(time>=art_start && time<=art_start+params->dt)
		{
		    int junk_index;
		    if (time < params->junk_days)
			junk_index = (int)(time / ((double)params->junk_days/params->junk_bins)) + 1;
		    else
			junk_index = params->junk_bins+1;

		    for (int junkies=1; junkies <= junk_index; junkies++)
		    {
			strain *junk_strain = Q_list->get_elem(junkies);
			start_junk_l1 += junk_strain->get_l1_cells(0);
			start_junk_l2 += junk_strain->get_l2_cells(0);
		    }
		    start_l1 = x[2];
		    start_l2 = x[3];
		    fprintf(stdout,"Run %d, faux-ART start at t=%lf...\n",
			runnum+1,time);
		    fprintf(stdout,"\tL1=%lf;L2=%lf,J1=%lf;J2=%lf\n",
			start_l1,start_l2,start_junk_l1,start_junk_l2);
		}
	    }

	    //#only record times in t_list
	    if (params->writeOn && time>=next_output)
	    {
		strain *junk_strain = Q_list->get_elem(1);
		if (params->dataF1 != NULL) 
		{
		    fprintf(params->dataF1,"%3.2lf",time);
		    fprintf(params->dataF1,",%u",x[0]);
		    fprintf(params->dataF1,",%u",x[1]);
		    fprintf(params->dataF1,",%u",x[2]);
		    fprintf(params->dataF1,",%u",x[3]);
		    fprintf(params->dataF1,",%3.2lf",(x[4] > 0)?log10(x[4]):0);
		    fprintf(params->dataF1,",%3.2lf",(x[5] > 0)?log10(x[5]):0);
		    fprintf(params->dataF1,",%3.2lf",(x[4]+x[5] > 0)?log10(x[4]+x[5]):0);
		    fprintf(params->dataF1,",%u",x[6]);
		    fprintf(params->dataF1,",%u",x[7]);
		    fprintf(params->dataF1,",%u",x[8]);
		    fprintf(params->dataF1,",%3.2e",Bt_avg);
		    fprintf(params->dataF1,",%3.2e",Bt_max);
		    fprintf(params->dataF1,",%3.2lf",params->Ro);
		    fprintf(params->dataF1,",%d",q_a);
		    fprintf(params->dataF1,",%d",qmax);
		    fprintf(params->dataF1,",%d",qmax_a);
		    fprintf(params->dataF1,",%d",qmax_l1);
		    fprintf(params->dataF1,",%d",qmax_l2);
		    fprintf(params->dataF1,",%d",qmax_dH);
		    fprintf(params->dataF1,",%lu",qmax_cd8s);
		    fprintf(params->dataF1,",%3.2lf",qmax_dA);
		    fprintf(params->dataF1,",%3.2e",qmax_beta);
		    fprintf(params->dataF1,",%3.2lf",real_avgHd);
		    fprintf(params->dataF1,",%d",real_maxHd);
		    fprintf(params->dataF1,",%3.2lf",entropy);
		    fprintf(params->dataF1,",%3.2lf",diversity);
		    fprintf(params->dataF1,",%3.2lf",seq_ident);
		    fprintf(params->dataF1,",%u",x[9]);
		    fprintf(params->dataF1,",%u",x[10]);
		    fprintf(params->dataF1,",%u",x[11]);
		    fprintf(params->dataF1,",%u",x[12]);
		    fprintf(params->dataF1,",%d",
			junk_strain->get_lumped_strains(0));
		    fprintf(params->dataF1,",%3.2lf",
			junk_strain->get_founder_dist());
		    fprintf(params->dataF1,",%d",
			junk_strain->get_total_prolif());
		    fprintf(params->dataF1,",%u",x[13]);
		    fprintf(params->dataF1,",%3.2e",qmax_dImmun);
		    fprintf(params->dataF1,",%3.2e",dImmun_avg);
		    fprintf(params->dataF1,"\n");
		    fflush(params->dataF1);
		}
		// strain output
		if (params->dataF2 != NULL) 
		{
		    // use max cell count to set bin ranges
		    int num_bins=params->strain_bins;
		    int bin_size;

		    if (params->max_strain_cells <= params->strain_bins)
			bin_size = 1;
		    else
			bin_size = (int)((double)params->max_strain_cells/num_bins + 0.5);

		    for (int i=0; i < num_bins; i++)
			strain_bins[i]=0;

		    for (int i=0; i < Q_list->get_num_elems(); i++)
		    {
			strain *q = Q_list->get_elem(i);
			int bin = q->get_cells()/bin_size;
			bin = MAX(0,MIN(num_bins-1,bin));
			strain_bins[bin]++;
		    }
		    fprintf(params->dataF2,"%3.2lf",time);
		    fprintf(params->dataF2,",%u",x[6]);
		    fprintf(params->dataF2,",%d",bin_size);
		    for (int j=0; j < num_bins; j++)
			fprintf(params->dataF2,",%u",strain_bins[j]);
		    fprintf(params->dataF2,"\n");
		    fflush(params->dataF2);
		}
		// top 10 strain VLs
		if (params->dataF3 != NULL) 
		{
		    // use max cell count to set bin ranges
		    fprintf(params->dataF3,"%3.2lf",time);
		    fprintf(params->dataF3,",%d",q_a);
		    fprintf(params->dataF3,",%3.2lf",(x[4] > 0)?log10(x[4]):0);
		    fprintf(params->dataF3,",%3.2lf",(x[5] > 0)?log10(x[5]):0);
		    fprintf(params->dataF3,",%d",params->plotTopStrains);
		    fprintf(params->dataF3,",%d",params->numTopStrains);
		    for (int j=0; j < MAX_STRAIN_COLORS; j++)
		    {
			if (j < params->numTopStrains)
			{
			    strain *this_strain = top_strains->get_elem(j);
			    fprintf(params->dataF3,",%d",this_strain->get_strain_index());
			    fprintf(params->dataF3,",%d",this_strain->get_act_cells(0));
			    fprintf(params->dataF3,",%3.2lf",(this_strain->get_via_virus(0)>0)?log10(this_strain->get_via_virus(0)):0);
			    fprintf(params->dataF3,",%3.2lf",(this_strain->get_junk_virus(0)>0)?log10(this_strain->get_junk_virus(0)):0);
			    fprintf(params->dataF3,",%lu",this_strain->get_cd8_cells(0));
			}
			else
			{
			    fprintf(params->dataF3,",-1");
			}
		    }
		    fprintf(params->dataF3,"\n");
		    fflush(params->dataF3);
		}
		// followed strains time history
		if (params->followTopStrains && params->dataF7 != NULL) 
		{
		    fprintf(params->dataF7,"%3.2lf",time);
		    fprintf(params->dataF7,",%u",top_strains->get_num_elems());
		    for (int j=0; j < MAX_FOLLOW_STRAINS; j++){
			if (j < top_strains->get_num_elems()) {
			    strain *thisStrain= top_strains->get_elem(j);
			    unsigned int num_tst_virions =
				(thisStrain->get_via_virus(0) +
				 thisStrain->get_junk_virus(0));
			    double log_virions = (num_tst_virions > 0)?
				log10(((double)num_tst_virions)/blood_factor):0;
			    fprintf(params->dataF7,",%3.2lf",log_virions);
			}
			else {
			    fprintf(params->dataF7,",0");
			}
		    }
		    fprintf(params->dataF7,"\n");
		    fflush(params->dataF7);
		}
		next_output=next_output+params->sampleInterval;
	    }
#ifdef GUI_ENABLED
	    if (time >= params->NextRefresh && params->points != NULL)
	    {
		if (x[8] > 0)
		    params->numStrains = x[8]-1;
		else
		    params->numStrains = 0;
		update_points(params,&snap_this_frame,snapnum,x[4],x[5],
		    top_strains,x[1]+x[2],x[3]+x[4]+x[5]+x[6],x[9],x[10],x[11],x[12]);
	    }
#endif
	      
	    //## evolve lists of cell objects based on events ##

	    // stop if log Vt is out of range early on
	    double logVt = ((x[4]+x[5]) > 0)?log10((double)(x[4]+x[5])/blood_factor):0;
	    if (time >= 2 && time < 5 && logVt < 0.5)
	    {
		char err_msg[200];
		sprintf(err_msg,"logVt out of range (0.5 < %3.2lf), stopping run at t=%3.2lf...\n",
		    logVt,time);
		fprintf(stdout,err_msg);

		//discard this run
		discarded++;
		this_discarded = true;
#ifdef GUI_ENABLED
		params->stopFlag = 1;
		ShowMessageBox("Error",err_msg);
#endif

		break;
	    }
	    // stop if no infected cells exist (latent or active)
	    else if (
	    (x[1] > 0 || x[2] > 0 || x[3] > 0 || x[13] > 0)
		)
	    {
		    model(params,time,0,x, Q_list,top_strains);

		params->Ro = params->aS*(Bt_avg/params->vol)*params->pi/params->gam/params->dS/params->dA; //# basic reproductive number


		//Adjust Vt score if at data measurement time
		if (params->vt_points[params->disp_patient-1] != NULL && (params->Crit_type == 1))
		{
		    for (int d1=0; d1 < params->num_vt_points[params->disp_patient-1];d1++)
		    {
			double tdiff1 = fabs(params->vt_points[params->disp_patient-1][d1].time - time);
			if (tdiff1 < 0.5*params->dt)
			{
			    double val1 = params->vt_points[params->disp_patient-1][d1].vt;
			    fprintf(stdout,"log Vt: %lf vs. %lf at t=%lf\n",
				    logVt,val1,time);
			    double weight = 1.0/(params->num_vt_points[params->disp_patient-1]);
			    score += fabs(val1-logVt)*weight;
			    break;
			}
		    }
		}
	    }
	    else
	    {
		char err_msg[200];
		sprintf(err_msg,"No infected cells left (A+L1+L2<=0), stopping run at t=%3.2lf...\n",
		    time);
		fprintf(stdout,err_msg);

		//discard this run
		this_discarded = true;
		discarded++;
#ifdef GUI_ENABLED
		params->stopFlag = 1;
		ShowMessageBox("Error",err_msg);
#endif

		break;
	    }

    #ifdef GUI_ENABLED
	    check_for_pause(params, &snap_this_frame);
    #endif
	    if (x[4]+x[5] >= peak_vl)
	    {
		if (time > 20)
		{
		    char err_msg[200];
		    sprintf(err_msg,"Peak is too late, stopping run at t=%3.2lf...\n",
			time);
		    fprintf(stdout,err_msg);

		    //discard this run
		    this_discarded = true;
		    discarded++;
    #ifdef GUI_ENABLED
		    params->stopFlag = 1;
		    ShowMessageBox("Error",err_msg);
    #endif

		    break;
		}
		peak_vl = x[4]+x[5];
		nadir_vl = x[4]+x[5];
		log_nadir_drop = 0;
		peak_time = time;
		nadir_time = time;
		num_set_pts=0;
		set_pt_err=0;
		set_pt_p=0;
	    }
	    else if (peak_vl > 0 && x[4]+x[5] < nadir_vl && 
		time - peak_time < 25. && num_set_pts==0)
	    {
		nadir_vl = x[4]+x[5];
		if (nadir_vl > 0)
		    log_nadir_drop = log10(peak_vl/blood_factor) - log10(nadir_vl/blood_factor);
		else
		    log_nadir_drop = log10(peak_vl/blood_factor);
		nadir_time = time;
		last_set_time = time;
	    }
	    else if (peak_time > 0 && nadir_time > peak_time + 5 && time - last_set_time > 1.)
	    {
		double logVL=(x[4]+x[5])>0?log10((x[4]+x[5])/blood_factor):0;
		if (set_pt_p == 0)
		{
		    set_pt_p = logVL;
		    set_pt_err=0;
		    fprintf(stdout,
			"1st Set point at t=%3.2lf, blood logVL=%3.2f\n", 
			time,logVL);
		    fprintf(stdout,
			"Peak was %3.2lf at t=%3.2lf, Nadir was %3.2f at t=%3.2lf\n", 
			(peak_vl >0)?log10(peak_vl/blood_factor):0,peak_time,
			(nadir_vl >0)?log10(nadir_vl/blood_factor):0,nadir_time);
		    // require set point < 10^6
		    if (logVL > 6)
		    {
			char err_msg[200];
			sprintf(err_msg,"Nadir is too high (> 10^6), stopping run at t=%3.2lf...\n",
			    time);
			fprintf(stdout,err_msg);

			//discard this run
			this_discarded = true;
			discarded++;
	#ifdef GUI_ENABLED
			params->stopFlag = 1;
			ShowMessageBox("Error",err_msg);
	#endif

			break;
		    }
		}
		else if (num_set_pts < MAX_SET_POINTS)
		{
		    set_pts[num_set_pts]=logVL;
		    set_pt_err += pow(set_pt_p - logVL,2.0);
		    num_set_pts++;
		}
		last_set_time = time;
	    }
	    
	    time_p=time;
	    time = time + params->dt;
	    params->time = time;

	    if (params->verbose && 
		(time_p == 0 || time >= next_print))
	    {
		unsigned int all[NUM_MEASURES];
		for (int i=0; i < NUM_MEASURES; i++)
			all[i] = x[i]; // single compartment

		// update stats to ensure consistency!
		gather_all_stats(params,time,all, Q_list,
		    &Bt_avg,&Bt_max,&q_a,&qmax,&qmax_via,&qmax_junk,
		    &qmax_a,&qmax_l1,&qmax_l2,&qmax_dH,&qmax_cd8s,
		    &qmax_cd8_index,&qmax_dA, &qmax_beta, &qmax_ns,
		    &real_avgHd,&real_maxHd,&qmax_dImmun,&dImmun_avg);

		fprintf(stdout,
		    "\nAt t=%3.2lf:S=%u,A=%u,L1=%u,L2=%u,logVt=%3.2lf (%3.2lf/%3.2lf),Q=%u(%d),CD8=%u(ic50=%lf)\n",
		    time,x[0],x[1],x[2],x[3],(x[4]+x[5] > 0)?log10(x[4]+x[5]):0,(x[4] > 0)?log10(x[4]):0,(x[5] > 0)?log10(x[5]):0,x[8],params->youngQ,x[9],params->cd8_ic50);


		
		fprintf(stdout,
		    "\tQa=%d,domQ=%d,V_domQ(%3.2lf/%3.2lf),A_domQ=%d,L1_domQ=%d,L2_domQ=%d,dH_domQ=%d,Bt_domQ=%3.2e,Bt_domQns=%d,",
		    q_a,qmax,(qmax_via>0)?log10(qmax_via):0,(qmax_junk>0)?log10(qmax_junk):0,qmax_a,qmax_l1,qmax_l2,qmax_dH, qmax_beta,qmax_ns);

		if (params->immun_model <= 1)
		    fprintf(stdout, "dA_domQ=%3.2lf\n", qmax_dA);
		else if (params->immun_model == 6)
		    fprintf(stdout, "dImmun_domQ=%3.2e\n", qmax_dImmun);
		else
		{
		    if (params->immun_model == 2 || params->immun_model == 3)
			fprintf(stdout, "cd8_domQ=%lu\n",qmax_cd8s);
		    else if (params->immun_model == 4 || params->immun_model == 5) {
			//int qa_cd8_index = 1+MIN(params->maxHammingDist/params->cd8GroupSize,qmax_dH/params->cd8GroupSize);
			unsigned long int qa_cd8_cells=0;
			unsigned long int qa_inf_cells=0;
			for (int compartment=0; compartment < num_compartments; compartment++)
			{
			    qa_cd8_cells+=params->cd8_groups[compartment][qmax_cd8_index].cd8_cells;
			    qa_inf_cells+=params->cd8_groups[compartment][qmax_cd8_index].inf_cells;
			}
			fprintf(stdout, "cd8_domQ=%lu,cd8_inf_domQ=%lu\n", qa_cd8_cells,qa_inf_cells);
		    }
		}


		int junk_index;
		if (time < params->junk_days)
		    junk_index = (int)(time / ((double)params->junk_days/params->junk_bins)) + 1;
		else
		    junk_index = params->junk_bins+1;

		for (int junkies=1; junkies <= junk_index; junkies++)
		{
		    strain *junk_strain = Q_list->get_elem(junkies);
		    //strain *junk_strain = Q_list->get_elem(1);
		    double junk_dH = junk_strain->get_founder_dist();
		    int junk_strains = junk_strain->get_lumped_strains(0);
		    int viable_strains = 0;
		    int prolif_cells = junk_strain->get_total_prolif();
		    int dead_cells = junk_strain->get_total_deaths();

		    int junk_l0 = junk_strain->get_act_cells(0);
		    int junk_l1 = junk_strain->get_l1_cells(0);
		    int junk_l2 = junk_strain->get_l2_cells(0);
		    unsigned long int junk_cd8s = junk_strain->get_cd8_cells(0);
		    fprintf(stdout,
			"\tJunk DNA cells: index=%d,cells=%d/%d/%d, strains=%d (p=%d,d=%d), dH=%3.2lf\n",
			junkies,junk_l0,junk_l1,junk_l2,junk_strains,
			prolif_cells,dead_cells,junk_dH);
		    if (params->immun_model == 2 || params->immun_model == 3)
			fprintf(stdout, "cd8_domQ=%lu\n",junk_cd8s);
		    else if ((params->immun_model == 4 || params->immun_model == 5) && params->kill_junk) {
			int qa_cd8_cells=params->cd8_groups[0][0].cd8_cells;
			int qa_inf_cells=params->cd8_groups[0][0].inf_cells;
			fprintf(stdout, "cd8_domQ=%d,cd8_inf_domQ=%d(ic50=%lf)\n", qa_cd8_cells,qa_inf_cells,params->dImmun_IC50);
		    }
		    else
			fprintf(stdout,"\n");
		    int viable_l0 = 0;
		    int viable_l1 = 0;
		    int viable_l2 = 0;
		    long unsigned int via_virus =0;
		    long unsigned int junk_virus =0;

		    prolif_cells = 0;
		    dead_cells = 0;

		    bool older=false;
		    for (int i=0; i < Q_list->get_num_elems(); i++)
		    {
			if (i==0 || i > params->junk_bins+1) 
			{
			    strain *viable_strain = Q_list->get_elem((int)i);
			    double birth=viable_strain->get_birthdate();
			    if (birth > params->junk_days)
				older=true;
			    else 
			    {
				int viable_index = (int)(birth / ((double)params->junk_days/params->junk_bins)) + 1;
				if (viable_index == junkies)
				{
				    viable_strains++;
				    viable_l0 += viable_strain->get_act_cells(0);
				    viable_l1 += viable_strain->get_l1_cells(0);
				    viable_l2 += viable_strain->get_l2_cells(0);
				    prolif_cells += viable_strain->get_total_prolif();
				    dead_cells += viable_strain->get_total_deaths();
				    via_virus +=viable_strain->get_via_virus(0);
				    junk_virus +=viable_strain->get_junk_virus(0);
				}
			    }
			}
			if (older)
			    break;
		    }
		    fprintf(stdout,
			"\tViable DNA cells: index=%d,cells=%d/%d/%d, strains=%d (p=%d,d=%d), Vnon=%lu, Vvia=%lu\n",
			junkies,viable_l0,viable_l1,viable_l2,viable_strains,
			prolif_cells,dead_cells,junk_virus,via_virus);
				    
		}

		fprintf(stdout,
		"\tBt avg=%3.2e max=%3.2e, dH avg=%4.3e max=%d, dImmun_avg=%3.2e Ro=%3.2lf\n",
		Bt_avg, Bt_max, real_avgHd, real_maxHd, dImmun_avg, params->Ro);
		next_print = time + params->print;
	    }
	    if (time > next_res_print)
	    {
		gather_reserv_stats(params,Q_list, time);
		next_res_print = time + params->res_days;
	    }
	}
	if (params->writeOn )
	{
	    // dump strain data
	    if (params->dataF4 != NULL) 
	    {
		char **non_syn_map = (char **)malloc(Q_list->get_num_elems()*sizeof(char *));
		int uniq_entries=0;
		for (int i=0; i < Q_list->get_num_elems(); i++)
		{
		    strain *c = Q_list->get_elem((int)i);
		    int my_entry = -1;

		    if (i > params->junk_bins+1)
		    {
			char *this_block = c->get_non_syn_ptr();
			bool found=false;
			for (int j=0; j < uniq_entries; j++)
			{
			    if (this_block==non_syn_map[j])
			    {
				found=true;
				my_entry = j;
				break;
			    }
			}
			if (!found)
			{
			    non_syn_map[uniq_entries] = this_block;
			    my_entry = uniq_entries;
			    uniq_entries++;
			}
		    }
		    unsigned int via_virus =c->get_via_virus(0);
		    unsigned int junk_virus =c->get_junk_virus(0);
		    double birth = c->get_birthdate();
		    fprintf(params->dataF4,"%d",c->get_strain_index());
		    fprintf(params->dataF4,",%3.2lf",birth);
		    fprintf(params->dataF4,",%3.2lf",c->get_founder_dist());
		    fprintf(params->dataF4,",%d",(c->get_parent() != NULL)?c->get_parent()->get_strain_index():-1);
		    fprintf(params->dataF4,",%e",c->get_beta());
		    fprintf(params->dataF4,",%d",c->get_cells());
		    fprintf(params->dataF4,",%d",c->get_total_cells());
		    fprintf(params->dataF4,",%3.2lf",c->get_last_event_date());
		    fprintf(params->dataF4,",%d",c->get_total_infected());
		    fprintf(params->dataF4,",%d",c->get_total_prolif());
		    fprintf(params->dataF4,",%3.2lf",c->get_exp_date());
		    fprintf(params->dataF4,",%u",via_virus);
		    fprintf(params->dataF4,",%u",junk_virus);
		    fprintf(params->dataF4,",%lu",c->get_cd8_cells(0));
		    fprintf(params->dataF4,",%lu",c->get_max_cd8_cells(0));
		    fprintf(params->dataF4,",%d",c->get_num_as());
		    fprintf(params->dataF4,",%d",c->get_num_cs());
		    fprintf(params->dataF4,",%d",c->get_num_gs());
		    fprintf(params->dataF4,",%d",c->get_num_ts());
		    if (c->get_parent() == NULL)
			fprintf(params->dataF4,",None");
		    else {
			if (c->get_num_as() < (c->get_parent())->get_num_as())
			    fprintf(params->dataF4,",A->");
			else if (c->get_num_cs() < (c->get_parent())->get_num_cs())
			    fprintf(params->dataF4,",C->");
			else if (c->get_num_gs() < (c->get_parent())->get_num_gs())
			    fprintf(params->dataF4,",G->");
			else
			    fprintf(params->dataF4,",T->");

			if (c->get_num_as() > (c->get_parent())->get_num_as())
			    fprintf(params->dataF4,"A");
			else if (c->get_num_cs() > (c->get_parent())->get_num_cs())
			    fprintf(params->dataF4,"C");
			else if (c->get_num_gs() > (c->get_parent())->get_num_gs())
			    fprintf(params->dataF4,"G");
			else
			    fprintf(params->dataF4,"T");
		    }
		    fprintf(params->dataF4,",%d",get_strain_ns(params,c));
		    fprintf(params->dataF4,",%d",my_entry);
		    fprintf(params->dataF4,"\n");
		}
		fprintf(params->dataF4,"%d\n",uniq_entries);
		for (int j=0; j < uniq_entries; j++)
		{
		    if (non_syn_map[j] != NULL)
			for (int k=0; k < (params->founder_stats.length / 8)+1; k++)
			    for (int l=0; l < 8; l++)
				if ((non_syn_map[j][k] & (1 << l)) ==0)
				    fprintf(params->dataF4,"0");
				else
				    fprintf(params->dataF4,"1");
		    else
			fprintf(params->dataF4,"0x0");
		    fprintf(params->dataF4,"\n");
		}
			
	    }
	    // dump followed strains (at bottom of TH file)
	    if (params->dataF7 != NULL) 
	    {
		fprintf(params->dataF7,"%3.2lf",params->tF);
		fprintf(params->dataF7,",%u",top_strains->get_num_elems());
		for (int j=0; j < top_strains->get_num_elems(); j++){
		    strain *thisStrain= top_strains->get_elem(j);
		    fprintf(params->dataF7,",%d",thisStrain->get_strain_index());
		}
		fprintf(params->dataF7,"\n");
		fflush(params->dataF7);
	    }
	    if (strain_bins != NULL)
		free(strain_bins);
	    if (age_bins != NULL)
		free(age_bins);
	}
	if (!this_discarded && params->Crit_type == 2)
	{
	    double adj_l1=start_l1*exp(-params->thL1*art_dur);
	    double adj_l2=start_l2*exp(-params->thL2*art_dur);

	    double adj_j1=start_junk_l1*exp(-params->thJ1*art_dur);
	    double adj_j2=start_junk_l2*exp(-params->thJ2*art_dur);

	    fprintf(stderr,"Run %d, faux-ART dur of %lf days...\n",
		runnum+1,art_dur);
	    fprintf(stderr,"\tL1=%lf->%lf; L2=%lf->%lf, J1=%lf->%lf, J2=%lf->%lf\n",
		start_l1,adj_l1,start_l2,adj_l2,start_junk_l1,adj_j1,start_junk_l2,adj_j2);

	    if (adj_l1 > start_l1 ||
	        adj_l2 > start_l2 ||
	        adj_j1 > start_junk_l1 ||
	        adj_j2 > start_junk_l2 )
	    {
		//discard this run
		char err_msg[200];
		sprintf(err_msg,"Defective ratio failure, stopping run at t=%3.2lf...\n",
		    time);
		fprintf(stdout,err_msg);
		discarded++;
	#ifdef GUI_ENABLED
		params->stopFlag = 1;
		ShowMessageBox("Error",err_msg);
	#endif
		this_discarded = true;
	    }
	    else
		junk_fract = (adj_j1+adj_j2)/(double)(adj_l1+adj_l2+adj_j1+adj_j2);
	}
	if (!this_discarded && params->Crit_type == 2)
	{
	    double the_err=0;
	    double the_var=0;
	    double the_mean=0;
	    if (num_set_pts > 0)
	    {
		the_mean=gsl_stats_mean(set_pts,1,num_set_pts);
		the_var=gsl_stats_variance_m(set_pts,1,num_set_pts,the_mean);
		for (int sp=0; sp < num_set_pts; sp++)
		    the_err += pow(set_pts[sp]-the_mean,2.0);
		the_err = the_err / (num_set_pts*the_mean);
	    }

	    log_peak_vl = (peak_vl >0)?log10(peak_vl/blood_factor):0;
	    log_nadir_vl = (nadir_vl >0)?log10(nadir_vl/blood_factor):0;
	    log_nadir_drop = log_peak_vl-log_nadir_vl;

	    if (params->score_sample > 0 &&
		num_div_samples < params->score_sample)
	    {
		int act_samples = sample_seqs_for_div(params,Q_list,0,div_sample_size,&div_samples);
		if (act_samples > 0)
		{
		    fprintf(stdout,"%d sampled sequences at t=%lf...\n",act_samples,time);
		    fprintf(stdout,"Seqs: ");
		    for (int i=0;i < act_samples-1; i++)
			    fprintf(stdout,"%d,",div_samples[i].seq_index);

		    fprintf(stdout,"%d\n",div_samples[act_samples-1].seq_index);

		    gather_sample_stats(params,Q_list, act_samples,div_samples,
			&samp_avgHd,&samp_maxHd,&entropy,&diversity,
			&params->ttca_mean,&seq_ident);
		    fprintf(stdout,
			"samp_avgHd=%lf,samp_maxHd=%d,diversity=%lf,seq_ident=%lf,entropy=%lf\n",
			samp_avgHd,samp_maxHd,diversity,seq_ident,entropy);

		    fprintf(stdout,
			"ttca_mean=%lf\n",params->ttca_mean);
		    free(div_samples);
		}
		else
		    fprintf(stdout,"No Active seq samples at =%lf\n",time);
	    }

	    //dH_rate = 100.*(samp_maxHd-set_pt_dH)/(params->tF-set_pt_time);
	    dH_rate = 100.*(samp_maxHd)/(params->tF);
	    div_rate = 100.*diversity/params->tF;


	    fprintf(stdout,
		"For run %d: log peak=%lf, peak time=%lf, log nadir=%lf, nadir time=%lf, set_pt mean=%lf, var=%lf (%d pts), dH rate=%lf, div rate=%lf, junk fract=%lf\n",
		runnum+1, log_peak_vl,peak_time, log_nadir_vl,nadir_time, 
	  	the_mean,the_var,num_set_pts,dH_rate,div_rate,junk_fract);
	    fflush(stdout);

	    avg_peak_vl += log_peak_vl;

	    avg_peak_time += peak_time;

	    avg_nadir_drop += log_nadir_drop;

	    avg_nadir_time += nadir_time;

	    avg_set_pt += the_mean;

	    avg_set_pt_var += the_var;

	    avg_dH_rate += dH_rate;

	    avg_div_rate += div_rate;

	    avg_junk_fract += junk_fract;
	}
	for (int s=0; s < Q_list->get_num_elems(); s++)
	{
	    strain *the_strain = Q_list->get_elem(s);
	    delete the_strain;
	}
	runnum++;
	int valid_runs = runnum-discarded;
	if (params->AIC_Runs > 0 && valid_runs==params->AIC_Runs)
	    break;
    }

    if (params->Crit_type == 2 && discarded <= runnum/2)
    {
	int valid_runs = runnum-discarded;
	fprintf(stdout,
	    "For %d runs: avg log peak=%lf, avg peak time=%lf, avg nadir drop=%lf, avg nadir time=%lf, avg set_pt=%lf, avg set_pt_var=%lf, avg dH rate=%lf, avg div rate=%lf, avg junk fract=%lf ( %d runs ignored)\n",
	    runnum, avg_peak_vl/valid_runs,avg_peak_time/valid_runs, 
	    avg_nadir_drop/valid_runs,avg_nadir_time/valid_runs,
	    avg_set_pt/valid_runs,avg_set_pt_var/valid_runs,
	    avg_dH_rate/valid_runs,avg_div_rate/valid_runs,
	    avg_junk_fract/valid_runs,discarded);

	score1 = pow(avg_peak_vl/valid_runs - params->peak_vl.mean,2.0) /
		params->peak_vl.variance;

	score2 = pow(avg_peak_time/valid_runs - params->peak_time.mean,2.0) /
		params->peak_time.variance;

	score3 = pow(avg_nadir_drop/valid_runs - params->nadir_drop.mean,2.0) /
		params->nadir_drop.variance;

	score4 = pow(avg_nadir_time/valid_runs - params->nadir_time.mean,2.0) /
		params->nadir_time.variance;

	score5 = pow(avg_set_pt/valid_runs - params->set_pt.mean,2.0) /
		params->set_pt.variance;

	score6 = pow(avg_set_pt_var/valid_runs - params->set_pt_var.mean,2.0) /
		params->set_pt_var.variance;

	score7 = pow(avg_dH_rate/valid_runs - params->dH_avg.mean,2.0) /
		params->dH_avg.variance;

	score8 = pow(avg_div_rate/valid_runs - params->div_rate.mean,2.0) /
		params->div_rate.variance;

	score9 = pow(avg_junk_fract/valid_runs - params->est_junk_fract.mean,2.0) /
		params->est_junk_fract.variance;

	score = params->weight1*score1 + params->weight2*score2 + params->weight3*score3 + params->weight4*score4 + params->weight5*score5 + params->weight6*score6 + params->weight7*score7 + params->weight8*score8 + params->weight9*score9 + 2*params->AIC_k;

	fprintf(stdout,
	    "\tscore1=%lf, score2=%lf, score3=%lf, score4=%lf, score5=%lf, score6=%lf, score7=%lf, score8=%lf, score9=%lf, score=%lf (AIC k=%d)\n",
	    score1,score2,score3,score4,score5,score6,score7,score8,score9,score,params->AIC_k);
    }
    else if (params->Crit_type == 2)
    {
	fprintf(stdout,
	    "Too many discarded runs (%d). Exiting with a bad score.\n",
	    discarded);
	score = BAD_SCORE;
    }

    return score;
}

//######################################################################
//#### the genetic simulations of primary infection

int main (int argc, char *argv[])
{
	int opt;

/** Program expects criteria in file hiv_sim.crit (or use -c <fname> option) and parameters from hiv_sim.in (or use -f <fname> option) **/
        char def_file[] = "hiv_sim.in";
	char criteria_file[] = "hiv_sim.crit";
	char *crit_file;
	crit_file = &criteria_file[0];
	crit_file = NULL;

	char viral_load_file[] = "hiv_sim.vl";
	char *vl_file;
	vl_file = &viral_load_file[0];
	vl_file = NULL;

	char def_seq_sampling_file[] = "seq_sampling.txt";
	char *seq_sampling_file;
	seq_sampling_file = &def_seq_sampling_file[0];
	seq_sampling_file = NULL;

	char *fasta_file;
	fasta_file = NULL;

	char *entropy_file;
	entropy_file = NULL;

	// data file variables
	int writeMask = 0;

	char dat_file1[] = "counts.csv";
	char dat_file2[] = "strains.csv";
	char dat_file3[] = "top_strains.csv";
	char dat_file4[] = "strain_lineage.csv";
	char dat_file7[] = "followed_strains.csv";

	int verbose = 0;

	const gsl_rng_type * rT;//pointer to gsl rand generator type
        unsigned long int seed;

	settings params;
        params.inp_file = &def_file[0];

	/* create a random number generator */
	gsl_rng_env_setup();
	rT = gsl_rng_default;

	params.ur = gsl_rng_alloc (rT);

	static struct option long_options[] =
        {
          /* These options set a flag. */
          {"verbose", no_argument,       0, 'v'},
          {"help", no_argument,       0, 'h'},
          {"sync",   no_argument,       0, 0},
          {"display",   no_argument,       0, 0},
          {0, 0, 0, 0}
        };
	/* getopt_long stores the option index here. */
	int option_index = 0;


	while((opt = getopt_long(argc, argv, "a:c:e:f:hlp:q:rs:vw:y:z:",
			long_options, &option_index)) != -1) {
	    switch (opt) {

	    case 0:
		/* If this option set a flag, do nothing else now. */
		if (long_options[option_index].flag != 0)
		    break;

		printf ("option %s", long_options[option_index].name);
		if (optarg)
		    printf (" with arg %s", optarg);
		printf ("\n");
		break;

	    case 'h':
		usage(argv[0]);
		exit(0);

	    case 'l':
		params.vl_as_log = 1;
		break;

	    case 'c':
		crit_file = optarg;
		fprintf(stderr,"Criteria file is %s\n",optarg);
		break;

	    case 'e':
		entropy_file = optarg;
		fprintf(stderr,"entropy file is %s\n",optarg);
		break;

	    case 'a':
		params.alt_inp_file = optarg;
		fprintf(stderr,"Alt input file is %s\n",optarg);
		break;

	    case 'y':
		seq_sampling_file = optarg;
		fprintf(stderr,"Seq samplingVL file is %s\n",optarg);
		break;

	    case 'z':
		vl_file = optarg;
		fprintf(stderr,"VL file is %s\n",optarg);
		break;

	    case 'f':
		params.inp_file = optarg;
		fprintf(stderr,"Input file is %s\n",optarg);
		break;
		
	    case 'p':
		if (sscanf(optarg,"%d",&params.disp_patient) != 1)
		{
		    fprintf(stderr,"Error reading -p <patient num>argument (%s)! Exiting\n",
			optarg);
		    exit(1);
		}
		fprintf(stderr,"Patient to fit is %d\n",params.disp_patient);
		break;

	    case 'q':
		fasta_file = optarg;
		fprintf(stderr,"FASTA file is %s\n",optarg);
		break;


	    case 'r':
		seed = time (NULL) * getpid();   
		gsl_rng_set (params.ur, seed);

		fprintf(stderr,"Random seed is %lu\n",seed);
		break;

	    case 's':
		if (sscanf(optarg,"%lu",&seed) != 1)
		{
		    fprintf(stderr,"Error reading random seed (%s)! Exiting\n",
			optarg);
		    exit(1);
		}
		gsl_rng_set (params.ur, seed);

		fprintf(stderr,"Random seed is %lu\n",seed);
		break;

	    case 'w':
		writeMask = atoi(optarg);
		params.writeOn = 1;
		fprintf(stderr,"Write mask is %d\n",writeMask);
		break;
	    
	    case 'v':
		verbose = 1;
		break;

	    default:
		usage(argv[0]);
		exit(1);
	    }
	}

	////////////////////////////////////////////////////////////////////////
	///// read input parameters through optional input file 
	if (params.inp_file != NULL)
	    read_input_file(params.inp_file, &params);

	////////////////////////////////////////////////////////////////////////
	///// read original sequence through optional FASTA file (from -q)
	if (fasta_file != NULL)
	    read_fasta_file(fasta_file, &params);

	////////////////////////////////////////////////////////////////////////
	///// read entropy values through optional entropy file (from -e)
	ifstream entF;
	entF.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit );

	if (entropy_file != NULL)
	{
	    if (fasta_file == NULL) {
		cerr << "Entropy option requires matching FASTA file "<<entropy_file<<"\n";
		exit(1);
	    }
	    try {
		int num_bases;
		entF.open(entropy_file);
		entF >> num_bases;
		if (num_bases > MAX_SEQ_LENGTH) {
		    cerr << "Too many bases in entropy file "<<entropy_file<<"\n";
		    exit(1);
		}
		if (num_bases != params.founder_stats.length ) {
		    cerr << "Bases in entropy file "<<entropy_file<<" must match FASTA sequence length ("<<params.founder_stats.length<<")\n";
		    exit(1);
		}
		for (int i=0; i < num_bases; i++) {
		    int base;
		    entF >> base >> params.bases[i] >> params.entropy[i];
		    cout <<base<<": "<<params.bases[i]<<" entropy: "<<params.entropy[i];
		    cout<<endl;
		}
		entF.close();
		if (params.disp_patient > params.num_patients)
		    params.disp_patient = params.num_patients;

	    } catch (ifstream::failure e)
	    {
		cerr << "Error reading clinical data file "<<entropy_file<<"\n";
		exit(1);
	    }
	    params.use_entropy = 1;
	}

	if (verbose)
	    params.verbose = 1;


	////////////////////////////////////////////////////////////////////////
	///// read plot samples through a viral load file that has them //////////////
	///// alternatively, read criteria file with peak level and time ranges, etc.
	ifstream vlF;
	vlF.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit );

	ifstream critF;
	critF.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit );

	if (vl_file != NULL)
	{
	    try {
		vlF.open(vl_file);
		vlF >> params.num_patients;
		if (params.num_patients > MAX_PATIENTS) {
		    cerr << "Too many patients in clinical data file "<<vl_file<<"\n";
		    exit(1);
		}
		for (int i=0; i < params.num_patients; i++) {
		    vlF >> params.num_vt_points[i];
		    params.vt_points[i] = (struct vt_coord *)my_malloc(params.num_vt_points[i]*sizeof (struct vt_coord));

		    cout <<"Number of Vt points:" << params.num_vt_points[i];
		    cout<<endl;

		    for(int j=0; j<params.num_vt_points[i]; j++)
		    {	
			vlF >> params.vt_points[i][j].time;
			cout <<j<< " time:"<< params.vt_points[i][j].time;
			cout<<endl;
		    }	
		    for(int j=0; j<params.num_vt_points[i]; j++)
		    {	
			double vtVal=0;
			vlF >> vtVal;

			// scale criteria from 1 ml value to current
			vtVal = vtVal * params.vol / 1000.0;

			if (params.vl_as_log)
			    params.vt_points[i][j].vt = vtVal;
			else
			    params.vt_points[i][j].vt = (vtVal > 0)?log10(vtVal):0;

			cout <<j<< " vt:"<< vtVal << " log: "<<params.vt_points[i][j].vt;
			cout<<endl;
		    }	
		    cout<<endl;
		}
		vlF.close();
		if (params.disp_patient > params.num_patients)
		    params.disp_patient = params.num_patients;

	    } catch (ifstream::failure e)
	    {
		cerr << "Error reading clinical data file "<<vl_file<<"\n";
		exit(1);
	    }
	    params.Crit_type = 1;
	}
	else if (crit_file != NULL)
	{
	    try {
		critF.open(crit_file);
		critF >>  params.peak_vl.mean >> params.peak_vl.variance;
		cout <<"log peak VL mean: "<<params.peak_vl.mean<<" variance: "<<params.peak_vl.variance;
		cout <<endl;
		critF >> params.peak_time.mean >> params.peak_time.variance;
		cout <<"peak time mean: "<<params.peak_time.mean<<" variance: "<<params.peak_time.variance;
		cout <<endl;
		critF >> params.nadir_drop.mean>> params.nadir_drop.variance;
		cout <<"log nadir drop mean: "<<params.nadir_drop.mean<<" variance: "<<params.nadir_drop.variance;
		cout <<endl;
		critF >> params.nadir_time.mean>> params.nadir_time.variance;
		cout <<"nadir time mean: "<<params.nadir_time.mean<<" variance: "<<params.nadir_time.variance;
		cout <<endl;
		critF >> params.set_pt.mean>> params.set_pt.variance;
		cout <<"Avg log set pt mean: "<<params.set_pt.mean<<" variance: "<<params.set_pt.variance;
		cout <<endl;

		critF >> params.set_pt_var.mean>> params.set_pt_var.variance;
		cout <<"log set pt variance mean: "<<params.set_pt_var.mean<<" variance: "<<params.set_pt_var.variance;
		cout <<endl;

		critF >> params.dH_avg.mean>> params.dH_avg.variance;
		cout <<"dH avg mean: "<<params.dH_avg.mean<<" variance: "<<params.dH_avg.variance;
		cout <<endl;

		critF >> params.div_rate.mean>> params.div_rate.variance;
		cout <<"Diversity rate mean: "<<params.div_rate.mean<<" variance: "<<params.div_rate.variance;
		cout<<endl;

		critF >> params.est_junk_fract.mean>> params.est_junk_fract.variance;
		cout <<"Est junk fraction mean: "<<params.est_junk_fract.mean<<" variance: "<<params.est_junk_fract.variance;
		cout<<endl;

		critF.close();
	    } catch (ifstream::failure e)
	    {
		cerr << "Error reading clinical data file "<<crit_file<<"\n";
		exit(1);
	    }
	    params.Crit_type = 2;
	}
	else
	    params.Crit_type = 0;

	///// optionally, read sequence samping schedule 
	ifstream sampF;
	sampF.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit );

	if (seq_sampling_file != NULL)
	{
	    try {
		sampF.open(seq_sampling_file);
		sampF >> params.num_seq_samplings;
		if (params.num_seq_samplings > MAX_SAMPLINGS) {
		    cerr << "Too many sampling points in sequence sampling file "<<seq_sampling_file<<"\n";
		    exit(1);
		}
		params.seq_sampling = (struct seq_samp *)my_malloc(params.num_seq_samplings*sizeof (struct seq_samp));
		for (int i=0; i < params.num_seq_samplings; i++) {
		    sampF >> params.seq_sampling[i].time;
		    cout <<i<< " time:"<< params.seq_sampling[i].time;
		    sampF >> params.seq_sampling[i].seqs;
		    cout << " seqs:"<< params.seq_sampling[i].seqs;
		    cout<<endl;
		}
		sampF.close();

	    } catch (ifstream::failure e)
	    {
		cerr << "Error reading clinical data file "<<seq_sampling_file<<"\n";
		exit(1);
	    }
	}
	if(((writeMask & 1) && ((params.dataF1 = fopen(dat_file1,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file1<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<1)) && ((params.dataF2 = fopen(dat_file2,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file2<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<2)) && ((params.dataF3 = fopen(dat_file3,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file3<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<3)) && ((params.dataF4 = fopen(dat_file4,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file4<<"\n";
	    exit(1);
	}
	if(((writeMask & (1<<6)) && ((params.dataF7 = fopen(dat_file7,"wt")) == NULL))){
	    cerr << "Could not open data file "<<dat_file7<<"\n";
	    exit(1);
	}

#ifdef GUI_ENABLED
	gui_main(argc, argv, &params);
#else
	ScoreFunction(&params);
#endif

	if(params.dataF1 != NULL) fclose(params.dataF1);
	if(params.dataF2 != NULL) fclose(params.dataF2);
	if(params.dataF3 != NULL) fclose(params.dataF3);
	if(params.dataF4 != NULL) fclose(params.dataF4);
	if(params.dataF7 != NULL) fclose(params.dataF7);
}
