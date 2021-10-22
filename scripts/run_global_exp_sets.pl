#!/usr/bin/perl
if ($#ARGV < 1) {
    my $args;
    $args=$#ARGV+1;
    die "This script requires 2 input arguments ($args given)\n";
}
my $infile;
my $outdir;

if ( -f $ARGV[0]) {
    $infile=$ARGV[0]
} else {
    die "1st arg is input file\n";
    exit (1);
}
if ( -f $ARGV[1]) {
    $critfile=$ARGV[1]
} else {
    die "2nd arg is the criteria file\n";
    exit (1);
}
#Best values
#Bt_mean 1.718308e-04
#Bt_max 2.586996e-04
#dImmun 4.953163e-04
#dImmun_IC50_0 4.61448
#cd8_ic50_0 3.17121
my $low_junk=0.8;
my $high_junk=0.99;
my $range_junk= $high_junk-$low_junk;

my $low_pi_mult=0.8;
my $high_pi_mult=1.2;
my $range_pi_mult= $high_pi_mult-$low_pi_mult;

my $low_beta=1e-4;
my $high_beta=1e-3;
my $range_beta= $high_beta-$low_beta;

my $low_mean_to_max=0.3;
my $high_mean_to_max=1.0;
my $range_mean_to_max= $high_mean_to_max-$low_mean_to_max;

my $low_dImmun=1e-5;
my $high_dImmun=1e-3;
my $range_dImmun= $high_dImmun-$low_dImmun;

my $low_dImmun_IC50_0=1;
my $high_dImmun_IC50_0=5;
my $range_dImmun_IC50_0= $high_dImmun_IC50_0-$low_dImmun_IC50_0;

my $low_cd8_ic50_0=0;
my $high_cd8_ic50_0=5;
my $range_cd8_ic50_0= $high_cd8_ic50_0-$low_cd8_ic50_0;

my $low_theta=0;
my $high_theta=2;
my $range_theta= $high_theta-$low_theta;

my $low_kappa=1;
my $high_kappa=4;
my $range_kappa= $high_kappa-$low_kappa;

my $new_dImmun= $range_dImmun*rand() + $low_dImmun;
my $new_dImmun_IC50_0= $range_dImmun_IC50_0*rand() + $low_dImmun_IC50_0;
my $new_cd8_ic50_0= $range_cd8_ic50_0*rand() + $low_cd8_ic50_0;
my $new_beta= $range_beta*rand() + $low_beta;
my $new_mean_to_max= $range_mean_to_max*rand() + $low_mean_to_max;
my $new_junk= $range_junk*rand() + $low_junk;
my $pi_mult= $range_pi_mult*rand() + $low_pi_mult;
my $new_theta= $range_theta*rand() + $low_theta;
my $new_kappa= $range_kappa*rand() + $low_kappa;

open (TSTFILE, ">>$infile");
printf TSTFILE ("junk %g\n",$new_junk);
printf TSTFILE ("Bt_k 1\n");
printf TSTFILE ("Bt_mean %e\n",$new_mean_to_max*$new_beta);
printf TSTFILE ("Bt_max %e\n",$new_beta);
printf TSTFILE ("pi %g\n",$pi_mult*(1.0/$new_beta));
printf TSTFILE ("pi_mult %g\n",$pi_mult);
printf TSTFILE ("dImmun %e\n",$new_dImmun);
printf TSTFILE ("dImmun_IC50_0 %g\n",$new_dImmun_IC50_0);
printf TSTFILE ("cd8_ic50_0 %g\n",$new_cd8_ic50_0);
#printf TSTFILE ("theta %g\n",$new_theta);
#printf TSTFILE ("kappa %g\n",$new_kappa);
printf TSTFILE ("use_dBeta 1\n");
printf TSTFILE ("use_weibull 1\n");
printf TSTFILE ("lifespan_decay 0\n");
printf TSTFILE ("AIC_k 6\n");
close(TSTFILE);

system("../hiv_sim -v -r -f $infile -c $critfile -w 1");
