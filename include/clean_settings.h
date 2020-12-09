//NOTE: This read-only file is generated from the file ../include/full_settings.h.
//Any editing should be done in that file.

/* The state class holds parameters that have been set either on the command line or 
 * in the input file.  
 */
#ifndef SETTINGS_H
#define SETTINGS_H

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#include <stdio.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include "cell.h"
#define MAX_TIME          10
#define MAX_VIRAL_LOAD    10000000000
#define MAX_EXPAND_CYCLES 8

#define MAX_HAMMING_DIST 60

#define MAX_COMPARTMENTS 1

#define MAX_CD8_GROUPS 1000
#define MAX_PATIENTS 100
#define MAX_SEQ_LENGTH 10000
#define MAX_SAMPLINGS 100

#define SAMPLE_RNA 0
#define SAMPLE_DNA 1

#include "colors.h"
#ifdef GUI_ENABLED
#include "plotpoints.h"
#endif

struct seq_stats {
    int length;
    int as;
    int cs;
    int gs;
    int ts;
};

struct cd8_group {
    unsigned long int cd8_cells;
    unsigned long int inf_cells;
};

struct vt_coord {
    double time;
    double vt;
};
typedef struct critType {
            double mean;
            double variance;
        } critType;

struct seq_samp {
    double time;
    int seqs;
};


class settings {
    public:
	settings(void);
	~settings(void);

	gsl_rng * ur;

/* set from input file */
	int use_cd8s;		//# flag to track CD8s for immune pressure

	int verbose;		//# extra messages

	double print;		//# interval for printing summary msgs

	int Crit_type;		//# which criteria set to use
	int Runs;		//# patients
	int AIC_Runs;		//# min for AIC (stop when reached; ignore if 0)

	double weight1;		//# weight for criteria 1
	double weight2;		//# weight for criteria 2
	double weight3;		//# weight for criteria 3
	double weight4;		//# weight for criteria 4
	double weight5;		//# weight for criteria 5
	double weight6;		//# weight for criteria 6
	double weight7;		//# weight for criteria 7
	double weight8;		//# weight for criteria 8
	double weight9;		//# weight for criteria 9

	double tF;			//# final time

	double Alt_input;	// time for reading in alternate inputs
	double Input_refresh;	// time for re-reading main inputs
	
	char *inp_file;		// input file (save in case reload reqd
	char *alt_inp_file;	// file for reading in alternate inputs

	int a0;			//#start with a0 active cells

	double dt;		//# time step
	double thL1;		//# latent cell type 1 decay rate from Siliciano
	double thL2;		//# latent cell type 2 decay rate from Siliciano
	double aL1 ;		//# avg proliferation rate of latent cell type 1
	double aL2 ;		//# avg proliferation rate of latent cell type 2

	double thJ1;		//# defective latent cell type 1 decay rate from Macallan/Besson
	double thJ2;		//# defective latent cell type 2 decay rate from Macallan/Golob

	double dA ;		//# death rate actively infected cell
	double xi1 ;		//# activation rate L1 to A
	double xi2 ;		//# activation rate L2 to A

	double dS ;		//# death rate of suseceptible cells
	double dI ;		//# death rate of infected cells
	double dActivated ;	//# death rate of activated/non-infected cells

	double gam;		//# virus clearance rate

	double mu0 ;		//# mutation rate pe4r nucleotide
	double mu ;		//# mutation rate for given sequence

	double junk;		//# junk DNA rate
	int kill_junk;		//# CD8s see & kill junk (L0)

	int junk_bins;		//#junk strain bins (divided by birth data)
	int junk_days;		//#used to determine width of junk bins
	int res_days;		//#used to determine print interval for reserv

	double sampleInterval ; //# interval for printing out info

//#other cell/virus parameters
	double pi ;		//# burst rate of virus from cells [virions/cell]
	double tau1;		//#artificially high for now, should be 1e-4    # probability of latency cell type 1[]
	double tau2;		//#artificially high for now, should be 1e-4    # probability of latency cell type 2[]
	double eclipse;		//optional eclipse period (rate)

	double egress_rate;	//cell egress rate (inter-compartmental)
	double egress_rate_v;	//virus egress rate (inter-compartmental)

	int S0;			//# initial susceptible cells

	double aS0;		//# 1ul growth rate of susceptible cells [fract/day]

	double aS;		//# growth rate of susceptible cells [fract/day]

	double Bt_mean;		//# mean infectivity of virus in 1ul [1/virus-day]
	double Bt_sd_factor;	//# stddev as mult of mean for infectivity 
	double Bt_k;		//# shape variable for Bt weibull

	double Bt_max;	//# max & founder beta 

	//double Bt_founder;	//# beta for initial infection only
	int use_dBeta;  // use "deltas" from parent for inherited betas
	int use_weibull;  // weibull vs log normal for beta "delta" or mean mult
	int use_entropy;  // entropy as fitness inheritance function

	double entropy_threshold; // entropy below which one can assume beta=0
	double entropy_factor; // used to map entropy to fitness inheritance

	double theta;

	double ART_start;	//# when to apply ART
	double ART_stop;	//# when to stop ART
	double ART_eff;		//# drug efficiency (0 for no drugs)

	int    immun_model; // which model of immune response to use
			    // 0 - no immune response
			    // 1 - global immune response
			    // 2 - strain-wise immune response
			    // 3 - combo global/strain-wise response
			    // 4 - strain-group response based on HD
			    // 5 - strain-group response based on birth time

	double dImmun_IC50; // scaled IC50
	double dImmun_IC50_0; // initial IC50 (pre-ramping)
	double ic50_k;

	int CD8_limit;
	double cd8_mean;
	double cd8_0;
	double cd8_k;
	double cd8_ic50;
	double cd8_ic50_0;
	double kappa;

	double tbirth0;
	double tbirth;
	double delta;

	double dA_ic50_0;	// ic50 for death rate (1/lifespan) on active cells
	double dA_ic50;		// scaled ic50

	double dA_scale;

	double lifespan_decay;//# decay rate for infected cell lifespan (half-life based on time from active cells)

	double dImmun;	// killing eff for background HIV T-cells
	double dImmun_s;// killing eff for strain-specific T-cells

	double vol;		//# simulation volume e.g. blood [uL]
	double dna_for_res;	//# sample volume e.g. blood [uL] used for resevoir
	double dna_for_seqs;	//# sample volume e.g. blood [uL] used for trees
	double transit_prob;	// Fraction of mutations that are transitions (versus transversions)
	int use_k2p;	// equal probability of mutations to any other base

//#volume dependent rates
	double Bt;		//# volume adjusted infection rate

	double aV;		//# growth rate due to viral activations
	double aV_ic50;		//# growth rate due to viral activations





//#calculated rates
	double dL1;		//# death rate of latent cell type 1
	double dL2;		//# death rate of latent cell type 2

	double res_factor;	//# used for bumping up latent cell death rate
	int res_model;		// 0 = no change, 1 = straight factor, 2 = IC50

	double l1 ;		//# latency cell type 1 factor
	double l2 ;		//# latency cell type 2 factor

	double Ro;		//# basic reproductive number

	int youngQ;		//# youngest strain w/ live cells

// variables for data file writing
	int writeOn;
	int age_bins;
	int strain_bins;
	int numStrains;

   	int maxTopStrains; // max # of significant strains to display(up to MAX_STRAIN_COLORS)
   	int numTopStrains; // current # of significant strains 
   	int numTopStrains2; // current # of significant strains in 2nd compart
   	int numTopStrains3; // current # of significant strains in 3rd compart
   	int cellThreshold; // log limit for showing traces

	struct seq_stats founder_stats;
	string founder_seq;
	int seq_length;

	double div_sample_interval;
	int div_sample_size; // number of virions to collect (see Keele 2008)
	int seq_sample_type;

	int num_seq_samplings;	// #seq sample times from file
	struct seq_samp *seq_sampling;

	int num_patients;
	int disp_patient;
	int score_sample;

	int vl_as_log;
	int num_vt_points[MAX_PATIENTS];
	struct vt_coord *vt_points[MAX_PATIENTS];

	critType peak_vl;
	critType peak_time;
	critType nadir_drop;
	critType nadir_time;
	critType set_pt;
	critType set_pt_var;
	critType dH_avg;
	critType div_rate;

	critType est_junk_fract;

	int div_use_logt;

	double diversity1;
	double diversity2;
	double div_interval;

	double art_start_mean;
	double art_start_std;

	double art_dur_mean;
	double art_dur_std;

	double ttca_mean;
	double ttca_var;
	int AIC_k;

	int AutoSnapshot;
	int snapshot;
	double SnapshotInterval;
	double NextSnap;

	double time;

	double max_log_v_value;
	double min_log_v_value;
	double max_log_a_value;
	double min_log_a_value;
	double max_log_cd8_value;
	double min_log_cd8_value;
	int max_strain_cells;
	int max_strains;

	int x_ticks;
	int y1_ticks;
	int y2_ticks;
	int y3_ticks;
	int y1_k;
	int y2_k;
	int y3_k;

   	int plotVt;
   	int plotAct;
   	int show_params;
   	int show_stats;
   	int plotCD8s;
   	int plotTopStrains;
   	int plotFirstStrains;
   	int followTopStrains;

   	int colorByHamming;
   	int colorBoxes;
   	int showNonVia;

   	int maxHammingDist;
   	int cd8GroupSize;

	struct cd8_group cd8_groups[MAX_COMPARTMENTS][MAX_CD8_GROUPS];

   	int scrollAxes;
   	int displayCompartment;

	double max_r0;
	double min_r0;
 	long int max_vl;
 	int max_cd8s;
 	int max_act;
	double max_time;

	double plot_bias;
	double plot_span;
	double time_bias;

	double refresh;
	double NextRefresh;
	double crit_start;

	int stopFlag;
	int pauseFlag;

#ifdef GUI_ENABLED
	plotPoints *points;
#endif

	unsigned long int *vt[MAX_STRAIN_COLORS];
	int *at[MAX_STRAIN_COLORS];
	double entropy[MAX_SEQ_LENGTH];
	char bases[MAX_SEQ_LENGTH];

	FILE *dataF1;
	FILE *dataF2;
	FILE *dataF3;
	FILE *dataF4;
	FILE *dataF5;
	FILE *dataF6;
	FILE *dataF7;

};
#endif
