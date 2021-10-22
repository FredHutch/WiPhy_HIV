#!/bin/bash
if [ "$1" != "" ]
then
sets=$1
else
sets=10
fi

if [ "$2" != "" ]
then
starts=$2
else
starts=1
fi

if [ "$3" != "" ]
then
count2=$3
else
count2=1
fi

critfile="../robb.crit"

# run all sims varying these parameters (via run_global_sets.pl)

count=$starts
while [ $sets -gt 0 ]
do
    thisdir="set_"
    thisdir+=$count
    if [ ! -d $thisdir ]
    then
    mkdir $thisdir
    fi
    cd $thisdir
    cp ../no_immun.in .
    sbatch -n 1 -t 20:00:00 ../run_no_immun_density_norm.pl $count2 no_immun.in $critfile
    cd ..

    sets=`expr $sets - 1`
    count=`expr $count + 1`
done
