#ifndef CELL_H
#define CELL_H

#define SUSCEPTIBLE 0
#define ACTIVE 1
#define LATENT1 2
#define LATENT2 3

#define PROLIF 0
#define INFECT 1

class cell {

public:
	cell(double birthday, int orig, int strain_num, int state);

	~cell(void);

	void set_state(int state);
	int get_state();
	int get_strain();
	int get_origin();
	double get_birthdate();

private:
	int cellState;
	int cellStrain;
	int cellOrigin;
	double cellBirth;
};
#endif
