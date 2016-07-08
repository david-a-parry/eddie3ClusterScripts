#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

die "Usage: $0 intervals.bed output_directory bam1 [bam2, bam3 ...]\n" if @ARGV < 3;

my $bed = shift;
my $out_dir = shift;

while (my $bam = shift){
    my ($f, $d) = fileparse($bam);
    (my $stub = $f) =~ s/\.[bs]am$//;
    my $script = "diagnose_targets_$stub.sh";
    my $output = "$out_dir/dt.min100.$stub.vcf";
    my $missing = "$out_dir/dt.missing.min100.$stub.txt";
    open (my $SCRIPT, ">$script") or die "Can't open script for writing: $!\n";
    print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -cwd
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/jre/1.7.0_60
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T DiagnoseTargets -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa  -L $bed -I $bam -o $output --missing_intervals $missing -min 100
EOT
;
    close $SCRIPT;
}
