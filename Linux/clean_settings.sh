#!/bin/sh
./strip_section.pl IMMUN_RAMP ../src/full_settings.cpp > temp1.cpp
./strip_section.pl IMMUN_EXPAND temp1.cpp > temp2.cpp
./strip_section.pl ANTIBODIES temp2.cpp > temp3.cpp
./strip_section.pl MULTI_COMPARTMENTS temp3.cpp > ../src/clean_settings.cpp
chmod -w ../src/clean_settings.cpp
rm temp1.cpp
rm temp2.cpp
rm temp3.cpp
