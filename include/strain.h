#ifndef STRAIN_H
#define STRAIN_H
#include "settings.h"
#include "generic_list.h"

// Strain object class
//

// Each strain contains the following information:
//
// 1. The strain number (0 for founder)
// 2. The strain's dates of emergence, activation and eventual exhaustion
// 3. The strain's parent (NULL if founder strain)
// 4. The strain's secondary parent (for super-infection / recombination)
// 5. The strain's beta (0 if non-infectious)
// 6. The strain's a/c/g/t counts from parent + mut (4 int counts)
// 7. The number of cells infected with it BY COMPARTMENT including
//    a. The number of actively productive infected cells
//    b. The number of latently infected cells of type 1 (med lifespan)
//    c. The number of latently infected cells of type 2 (longer lifespan)
// 8. The number of valid virions BY COMPARTMENT from those cells
// 9. The number of APOBEC virions BY COMPARTMENT from those cells
// 10. Optionally, the number of CD8s BY COMPARTMENT targeting the strain
// 11. Optionally, the level at which the strain would receive a 
//    half maximal immune response (IC50)
// 12. Optionally, the number of strains represented by this one strain (for 
//    cells with junk DNA & cells making non-viable virions also BY COMPARTMENT)
// 13. Optionally, the amount of virus being held by FDCs (follicle only)
// 14. Optional: ptr to bit mask for non-synonymous mutations (may point to parent's)
// FUTURE->16. First rep in recombination
// FUTURE->17. Start of recombinatory fragment
// FUTURE->18. Length of recombinatory fragment

class strain {

public:
	strain(settings *params, int num, double birthday, 
		strain *mother_strain, double inputBeta);

	~strain(void);

	int get_cells();
	int get_act_cells(int compartment);
	int get_ecl_cells(int compartment);
	int get_l1_cells(int compartment);
	int get_l2_cells(int compartment);

	unsigned long int get_cd8_cells(int compartment);
	int get_lumped_strains(int compartment);

	unsigned int get_via_virus(int compartment);
	unsigned int get_junk_virus(int compartment);

	unsigned int get_held_viable(int compartment);
	unsigned int get_held_junk(int compartment);

	double get_death_rate(int compartment,settings *params, double time);
	void set_death_rate(int compartment, double dA);

	int get_total_prolif();
	int get_total_infected();
	int get_total_deaths();
	int get_total_cells();

	int get_num_as();
	int get_num_cs();
	int get_num_gs();
	int get_num_ts();

	char *get_non_syn_ptr();

	void set_num_as(int);
	void set_num_cs(int);
	void set_num_gs(int);
	void set_num_ts(int);

	int get_num_bases();
	void set_num_bases(int );

	unsigned long int get_max_cd8_cells(int compartment);

	int add_infected(double infTime);
	int add_prolif(double prolifTime);

	int get_assoc();
	int inc_assoc();

	int inc_lumped_strains(int compartment);
	int dec_assoc(settings *params, generic_list<strain *> *ql, double decTime);

	int dec_ecl_cells(int compartment);
	int inc_ecl_cells(int compartment);

	int dec_act_cells(int compartment,double decTime);
	int inc_act_cells(int compartment,double incTime);

	int dec_l1_cells(int compartment);
	int inc_l1_cells(int compartment);

	int dec_l2_cells(int compartment);
	int inc_l2_cells(int compartment);

	unsigned long int dec_cd8s(int compartment,int change);
	unsigned long int inc_cd8s(int compartment,int change);

	unsigned long int set_cd8s(int compartment,unsigned long int total);

	int dec_via_virus(int compartment,int change);
	int inc_via_virus(int compartment,int change);

	int dec_junk_virus(int compartment,int change);
	int inc_junk_virus(int compartment,int change);

	int dec_held_viable(int compartment,int change);
	int inc_held_viable(int compartment,int change);

	int dec_held_junk(int compartment,int change);
	int inc_held_junk(int compartment,int change);

	strain *get_parent();

	int get_strain_index();
	double get_founder_dist();
	void set_founder_dist(double fdist);

	double get_birthdate();
	void set_birthdate(double time);

	double get_active_date();

	double get_beta();
	double get_dImmun_s();
	double get_ic50(settings *params, double time);
	double get_last_event_date();
	double get_exp_date();

	int get_cd8_group();
	void set_cd8_group(int group);

	static double get_delta_hamming(settings *params);
	static int get_cd8_groups();

private:
	static int cd8Groups;

	int strainNum;
	int cd8Group;
	int lumpedStrains[MAX_COMPARTMENTS];
	double founderDist;
    	double actDeathRate[MAX_COMPARTMENTS];

	int totalCells; // total ever assoc
	int assocCells; // total currently assoc cells (latent or active)

	int activeCells[MAX_COMPARTMENTS]; // current actively infected cells
	int l1Cells[MAX_COMPARTMENTS]; // current short term latently infected cells
	int l2Cells[MAX_COMPARTMENTS]; // current long term latently infected cells
	int eclipseCells[MAX_COMPARTMENTS]; // current pre-productive cells (in eclipse)

	unsigned long int cd8Cells[MAX_COMPARTMENTS]; // total currently assoc cd8 cells (immune response)
	unsigned long int cd8MaxCells[MAX_COMPARTMENTS]; // highest assoc cd8 cells (immune response)

	unsigned int via_virions[MAX_COMPARTMENTS]; // total currently assoc viable virus (excluding FDC holdings)
	unsigned int junk_virions[MAX_COMPARTMENTS]; // total currently assoc junk virus (excluding FDC holdings)
	unsigned int held_viable[MAX_COMPARTMENTS]; // total currently FDC assoc viable virus
	unsigned int held_junk[MAX_COMPARTMENTS]; // total currently FDC assoc junk virus

	int totalInfected;  // total ever infected (beta>0)
	int totalProlif;   // total proliferations
	int totalDeaths;   // total deaths

	strain *parentStrain;
	strain *alt_parentStrain;
	double strainBirth;

	double dImmun_s;

	double lastEventDate;	// last prolif or infect event
	double expireDate;	// when assoc cells hit 0

	double activeDate;  // reset on activation for L1 and L2, on 1st active
			  // or during reactivation (adjusted for down time)
	double deactiveDate;  // reset when active count goes to 0

	double strainBeta; // how infectious is this strain?
	double strainIC50;// how recognizable is this strain by T-cells? (increases with time?)
	double strainIC50_0;// how recognizable is this strain by T-cells? (increases with time?)

	int num_as;
	int num_cs;
	int num_gs;
	int num_ts;
	int num_bases;

	char *non_syn_ptr;

	//int fragment_start;
	//int fragment_len;

	//strain *start_strain;

};
#endif
