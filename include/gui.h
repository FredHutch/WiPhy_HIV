#ifndef GUI_H
#define GUI_H 1
#include "colors.h"
#include "settings.h"
#include "strain.h"  
#include "generic_list.h"

int gui_main (int argc, char *argv[], settings *vars);
void update_points(settings *vars, bool *snap_this_frame, int snapnum,
	unsigned int vi_Vt,unsigned int nonvi_Vt, generic_list<strain *> *topQs,
	unsigned int acts, unsigned int lats,
	unsigned int cd8s, unsigned int junk_L0s, unsigned int junk_L1s, unsigned int junk_L2s);

void check_for_pause(settings *vars, bool *snap_this_frame);
bool ShowMessageBox(std::string , std::string message);
#endif

