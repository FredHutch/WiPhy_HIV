#!/bin/sh
./strip_section.pl IMMUN_RAMP ../include/full_settings.h > temp1.h
./strip_section.pl IMMUN_EXPAND temp1.h > temp2.h
./strip_section.pl ANTIBODIES temp2.h > temp3.h
./strip_section.pl MULTI_COMPARTMENTS temp3.h > ../include/clean_settings.h
rm temp1.h
rm temp2.h
rm temp3.h
