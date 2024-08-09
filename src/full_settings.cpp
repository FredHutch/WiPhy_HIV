#include <cstdlib>
#include<string>
using namespace std;

#include <stdio.h>
#include <string.h>

#ifdef GUI_ENABLED
#include "colors.h"
#endif

#include "settings.h"

settings::settings(void) 
{
    Runs = 1;
    AIC_Runs = 0;

    weight1 = 1.0;		//# weight for criteria 1
    weight2 = 1.0;		//# weight for criteria 2
    weight3 = 1.0;		//# weight for criteria 3
    weight4 = 1.0;		//# weight for criteria 4
    weight5 = 1.0;		//# weight for criteria 5
    weight6 = 1.0;		//# weight for criteria 6
    weight7 = 1.0;		//# weight for criteria 7
    weight8 = 1.0;		//# weight for criteria 8
    weight9 = 1.0;		//# weight for criteria 9

    a0=10;

    aV=0.;
    aV_ic50=1000.;
    aS0=0.02;

    thL1 = -5.2e-4;
    thL2 = -5.2e-4;

    thJ1 = 0.02;//# defective latent cell type 1 decay rate from Macallan/Besson
    thJ2 = 7.16e-5;//# defective latent cell type 2 decay rate from Macallan/Golob

    aL1  = 0.015;
    aL2  = 0.002;

    dA  = 1.0;
    dS  = 0.02;

    dActivated  = 0.183;// net clearance of activated CD4s based on growth-death
			// Half-life of activated CD4s = 3 days (DeBoer, 2003)

    res_factor = 2;	// IC50 on log Vt for bumping up death rate of latently inf cells
    res_model = 0;		// 0 = no change, 1 = straight factor, 2 = IC50

    xi1  = 57/1e6;
    xi2  = 57/1e6;

    gam = 23;
    mu0  = 2e-5;
    junk = 0.9;
    kill_junk = 0;
    junk_bins = 1;
    junk_days = 365;
    res_days = 30;
    transit_prob  = 0.666;
    use_k2p = 0;

    pi  = 2e3;
    ART_start = 0;
    ART_stop = 0;
    ART_eff = 0;

    Bt_mean  = 2e-5;
    Bt_sd_factor  = 0; // const beta

    //Bt_founder  = 2.5e-5;
    Bt_max  = 3e-4;
    Bt_k  = 1.5;

    lifespan_decay  = 0.;

    use_dBeta = 1;
    use_weibull = 1;
    use_entropy = 0;

    entropy_threshold = 0.5; // entropy below which one can assume beta=0
    entropy_factor = 1; // used to map entropy to fitness inheritance

    eclipse = 0;

    ic50_k = 1.5;

    CD8_limit = 0;
    cd8_mean = 100;
    cd8_0 = 10;
    cd8_k = 1.5;
    cd8_ic50_0 = 1;
    kappa = 0.1;

    theta = 2.8;

    immun_model = 1; // global immune response

    for (int i=0; i < MAX_COMPARTMENTS; i++)
	for (int j=0; j < MAX_CD8_GROUPS; j++) {
	    cd8_groups[i][j].cd8_cells = 0;
	    cd8_groups[i][j].inf_cells = 0;
	}

    dImmun_IC50_0 = 1;

    delta = 1e-3;
    tbirth0 = 5e-3;

    dA_ic50_0 = 0;
    dA_scale = 3;

    div_use_logt = 1;
    div_interval=10;
    diversity1=0;
    diversity2=0;

    ttca_mean=0;
    ttca_var=0;

    dImmun = 0.0001;
    dImmun_s = 0.01;

    tau1 = 0.01;
    tau2 = 0.01;

    vol = 1000;
    dna_for_res = 200;
    dna_for_seqs = 100;

#ifdef IMMUN_RAMP
    immun_delay_time = 2; 	// time before ramp up of immunity
    immun_ramp_time = 7; 	// time to ramp up immunity
#endif

#ifdef IMMUN_EXPAND
    expand_freq = 0.01; // about 3.5 times per year
    expand_cycles = 8;	// 8 cycles (max)
#endif

#ifdef MULTI_COMPARTMENTS
    vol2 = 600;
    vol3 = 400;

    a2_0=0;
    aS2_0=0.02;

    gam2 = 5.5;

    ART_eff2 = 0;

    Bt2_factor  = 1.0;
    dA2_factor  = 1.0;

    a3_0 = 0;		//#start compartment 3 with a3_0 active cells
    gam3 = 1;		//# virus clearance rate (3rd compartment)

    egress_rate_v = 0.5;

    aS3_0 = 0.02;;		//# 1ul growth rate of compartment 3 susceptible cells [fract/day]
    Bt3_factor = 1;	//# difference in infectivity of virus in 2nd compartment
    dA3_factor  = 1.0;

    numStrains = 0;
    numTopStrains2 = 0;
    numTopStrains3 = 0; // current # of significant strains in 3rd compart

    threeCompartments = 0;
    ART_eff3 = 0;	//# drug efficiency (compartment 3)

    egress_rate = 0.0005;
    fdc_bind = 0.1;	//#fraction of free virus to bind each day
    fdc_release = 0.01;	//#fraction of bound virus to free each day
#endif

    tF=100;

    dt=0.1;

    sampleInterval=0.5;

    age_bins=10;
    strain_bins=10;

    founder_stats.length = MAX_SEQ_LENGTH;
    founder_stats.as = MAX_SEQ_LENGTH/4;
    founder_stats.cs = MAX_SEQ_LENGTH/4;
    founder_stats.gs = MAX_SEQ_LENGTH/4;
    founder_stats.ts = MAX_SEQ_LENGTH/4;

    founder_seq = "";
    div_sample_interval = 100;
    div_sample_size = 1000;
    seq_sample_type = SAMPLE_RNA;
    seq_length = 2600; // ENV gene size (by default)
    num_seq_samplings = 0;

    use_cd8s = 0;
    max_strain_cells=1;
    cellThreshold = 3;

#ifdef GUI_ENABLED
    max_strains=MAX_RAINBOW_COLORS;
#endif

    score_sample = 0;
    colorByHamming = 1;
    maxHammingDist = 100;
    colorBoxes = 8;
    showNonVia = 0;
    cd8GroupSize = 2;
    AIC_k = 0;

    dataF1=NULL;
    dataF2=NULL;
    dataF3=NULL;
    dataF4=NULL;
    dataF5=NULL;
    dataF6=NULL;
    dataF7=NULL;

    writeOn = 0;
    verbose = 0;
    print = 5;

    Crit_type = 1;
    for (int i=0; i < MAX_PATIENTS; i++) {
	num_vt_points[i] = 0;
	vt_points[i] = NULL;
    }
    disp_patient = 1;
    vl_as_log = 1;

    plotVt = 1;
    show_params = 0;
    show_stats = 0;
    plotCD8s = 0;
    plotAct = 0;
    plotTopStrains = 0;
    plotFirstStrains = 0;
    followTopStrains = 0;

    displayCompartment = 0;
#ifdef MULTI_COMPARTMENTS
    twoCompartments = 0;
    threeCompartments = 0;
#endif

    max_vl = MAX_VIRAL_LOAD;
    max_cd8s = 1000;
    max_act = 100000;
    max_time = MAX_TIME;

    time_bias=0;
    plot_span=10;
    plot_bias=0;

    refresh=0.5;

    SnapshotInterval=5;
    AutoSnapshot=0;
    snapshot=0;

#ifdef GUI_ENABLED
    points = NULL;
#endif

    stopFlag = 0;
    pauseFlag = 0;

    NextSnap = 0;

    max_log_v_value = 8;
    min_log_v_value = 0;
    max_log_a_value = 6;
    min_log_a_value = 0;
    max_log_cd8_value = 5;
    min_log_cd8_value = 0;

    x_ticks = 10;
    y1_ticks = 8;
    y2_ticks = 6;
    y3_ticks = 5;

    NextRefresh = 0;
    crit_start = 0;

    art_start_mean = 43.4;
    art_start_std = 24.5;

    art_dur_mean = 1696.6;
    art_dur_std = 1539.2;

    maxTopStrains = MAX_STRAIN_COLORS;
    numTopStrains = 0;

    Input_refresh = 0;	// time for re-reading std inputs
    Alt_input = 0;	// time for reading in alternate inputs
	
    inp_file = NULL;	// file for reading in inputs
    alt_inp_file = NULL;	// file for reading in alternate inputs
}

settings::~settings(void) 
{
}
