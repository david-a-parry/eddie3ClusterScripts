#!/usr/bin/env perl 
use strict;
use warnings;

die <<EOT

Usage: $0 input.vcf.gz output.vcf.gz /path/to/CADD/bin/score.sh

EOT
if @ARGV != 3;

my ($in, $out, $scorer) = @ARGV;

my $cmd = '';
my $gz = $in; 
if ($in !~ /\.gz$/){
    $cmd = "gzip $in\n";
    $gz = "$in.gz";
}

if ($out !~ /\.gz$/){
    $out .= '.gz';
}

$cmd .= "$scorer $gz $out";

print <<SCRIPT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -V
#\$ -cwd
#\$ -l h_rt=12:00:00
#\$ -l h_vmem=16G
# Configure modules
. /etc/profile.d/modules.sh
#Load modules

$cmd

SCRIPT
;

