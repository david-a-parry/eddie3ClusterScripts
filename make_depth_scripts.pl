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
    my $script = "depth_$stub.sh";
    my $output = "$out_dir/depth.min100.$stub.vcf";
    my $missing = "$out_dir/depth.missing.min100.$stub.txt";
    open (my $SCRIPT, ">$script") or die "Can't open script for writing: $!\n";
    print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -cwd
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/jre/1.7.0_60
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T DepthOfCoverage -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa  -L $bed -I $bam -o $output -ct 9 -ct 29 -ct 99 -ct 199 -ct 299 -ct 399 -ct 499 -ct 999 -ct --start 1 --stop 5000 --nBins 4999
EOT
;
    close $SCRIPT;
}
