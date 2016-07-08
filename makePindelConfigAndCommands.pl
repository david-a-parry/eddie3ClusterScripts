#!/usr/env perl 
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my @bams;
my $distance = 200;
my $config = "pindel_config.txt";
my $out = "pindel_";
my $script = "pindel_command.sh";
my $threads = 1;
my $help;
GetOptions(
            'i|b|bams=s{,}' => \@bams,
            'd|distance=i'  => \$distance,
            'c|config=s'    => \$config,
            'o|output=s'    => \$out,
            's|script=s'    => \$script,
            't|threads=s'   => \$threads,
            'h|?|help'      => \$help,
);

usage() if $help;
usage("at least one bam is required!\n") if not @bams;

open (my $CONF, ">$config") or die "Can't open $config for writing: $!\n";
open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";

foreach my $bam (@bams){
    my $f = fileparse($bam);
    my $sample;
    if ($f =~ /(\S+)-1[45]\d{4}/){
        $sample = $1;
    }elsif ($f =~ /(\S+)-exome/){
        $sample = $1;
    }else{
        die "Could not parse sample from file: $f\n";
    }
    
    print $CONF "$bam\t$distance\t$sample\n";
}
close $CONF;

print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -cwd
#\$ -e $script.stderr
# Configure modules
. /etc/profile.d/modules.sh
# Load modules

/export/users/dparry/pindel/pindel -f /export/users/dparry/GRCh37/hs37d5.fa -i $config -o $out -T $threads -c ALL
EOT
;

close $SCRIPT;

print "Created qsub script: $script\n";


sub usage{
    my $msg = shift;
    if ($msg){
        print $msg;
    }        
    print <<EOT

    Usage: $0 -i bam1 [bam2, bam3...] -d [fragment distance]

    Options:

    -i,-b,--bams,
        One or more bam files to process
    -d,--distance,
        Expected fragment size
    -c,--config,
        Config file to create for pindel run
    -o,--output,
        Output prefix for pindel command
    -s,--script,
        qsub script to create
    -t,--threads,
        no. threads for pindel command
    -h,-?,--help
        Show this help message and exit

EOT
;
    exit 1 if $msg;
    exit;
}

