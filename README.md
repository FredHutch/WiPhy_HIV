Stochastic, Strain-based, ODE model for early HIV infection (written in C++ with optional OpenGL/GTK+ 2.0 GUI)

This repository contains a stochastic HIV model written in C++ for building on MAC (OSX) and Linux (UBUNTU) platforms.
The optional GUI version is only available for the Linux distribution at this time.  It requires GTKGLEXT, GTK 2.0 and OpenGL.
This work was based on an HSV model written a number of years ago, so the required packages may not be available.
An effort has been made to make older libraries available in a separate GitHub repository (see note below). 

This model was used in support of this paper:

    Evolution during primary HIV infection does not require adaptive immune selection, medRxiv, Swan, et al., 2020.
    doi: https://doi.org/10.1101/2020.12.07.20245480

Under UBUNTU, this model requires several libraries.  You will also need to load the Mesa module to compile and 
run it.  Depending on your UBUNTU version other modules may also be required including:

    GSL,GLib,GTK+,libGLU and Pango (version with pangox)

To help in case required versions of certain libraries are unavailable, some required libraries and header files are contained 
in the extra_libs directory of the GIT repository FredHutch/HSV_models.  The path to the extra_libs directory must be set as 
the value of the GTKGLEXT environment variable before running make in the Linux directory.

