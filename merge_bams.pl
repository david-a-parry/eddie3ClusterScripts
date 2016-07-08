#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my $merge = 200;
my $out_prefix = "merged";
my $help;

GetOptions(
    "merge=i"   => \$merge,
    "out=s"     => \$out_prefix,
    "help"      => \$help,
);

usage() if $help;
usage("Please specify at least 2 input files on command line") if @ARGV < 2;
usage("--merge option must be a value of 2 or more") if $merge < 2;
usage("--merge option must lower than the number of input bams ($merge vs " . scalar(@ARGV) . ").") if $merge < 2;

my $merge_no = 0;
my @to_merge = ();
print <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m a
#\$ -m b
#\$ -m e
#\$ -cwd
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/jre/1.7.0_60
EOT
;

for (my $i = 0; $i < @ARGV ; $i++){
    push @to_merge, $ARGV[$i]; 
    if (not @to_merge % $merge){
        merge_files(\@to_merge);
        @to_merge = (); 
    }
}
merge_files(\@to_merge);

sub merge_files{
    my $f = shift;
    next if not @$f;
    $merge_no++;
    my $out_bam = "$out_prefix-$merge_no.bam";
    my $merge_cmd = "java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/picard/picard-tools-1.131/picard.jar MergeSamFiles MSD=TRUE SO=coordinate TMP_DIR=/mnt/lustre2/dparry/tmp/  CREATE_INDEX=TRUE I=\""
    . join("\" I=\"", @$f)   
    . "\" O=\"$out_bam\""  
    ;
    print "$merge_cmd\n";
    
}


sub usage{
    my $msg = shift;
    print STDERR "Error: $msg\n" if $msg;
    print <<EOT

Usage: $0 bam1 bam2 [bam3 ...] 

Options:

    -m,--merge
        no. bam files to include in each merged output file. Default = 200.
    -o,--out
        prefix for each merged bam file created. Default = 'merged'.
    -h,--help
        print this message and exit.

EOT
;
    exit 1 if $msg;
    exit;
}
