#!/usr/env perl 
use strict;
use warnings;
use File::Basename;

foreach my $bam (@ARGV){
    my $f = fileparse($bam);
    (my $out = $f) =~ s/\.bam$//;
    my $script = "getInsertSizeMetrics_$out.sh";
    open (my $SCRIPT, ">$script") or die "Cannot open $script for writing: $!\n";
    $out .= "_insertSizeMetrics";
print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -cwd
#\$ -e $script.output
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/jre/1.7.0_60
module load apps/gcc/R/3.1.0
java -Xmx4g -jar /export/users/dparry/picard/picard-tools-1.131/picard.jar CollectInsertSizeMetrics I=$bam O=$out.txt H=$out.hist.pdf
EOT
;
}
