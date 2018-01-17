#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my @pre_wait = ();
my %opts = (p => 1);
GetOptions(
    \%opts,
    "p|processes=i",
    "d|dummy",
    "w|wait_ids=s",
    "h|help",
) or usage("Error in options");
usage() if $opts{h};
usage("At least two scripts are required") if @ARGV < 2;
my $n = 0;
die "--processes option must be greater than 0\n" if $opts{p} < 1;
my @holds = ();
while (my $script = shift){
    my $cmd = "qsub ";
    my $hidx = $n - $opts{p};
    if (@holds > $hidx and $hidx >= 0){
        $cmd .= "-hold_jid $holds[$n - $opts{p}] ";
    }elsif ($opts{w}){
        $cmd .= "-hold_jid $opts{w} ";
    }
    $cmd .= $script;
    print STDERR "Executing: $cmd\n";
    if ($opts{d}){
        push @holds, $n;
    }else{
        my $output = `$cmd`;
        check_exit($?);
        if ($output =~ /Your job (\d+) .* has been submitted/){
            push @holds, $1;
        }else{
            my $msg = "Error parsing qsub output for $cmd!\n".
                      "Output was: $output\n";
            if (@holds){
                $msg .= "You may wish to terminate the previously submitted ".
                        "jobs. IDs of these jobs are: ".join(",". @holds)."\n";
            }
            die $msg;
        }
    }
    $n++;
}

##################################################
sub check_exit{
    my $e = shift;
    if($e == -1) {
        print "Failed to execute: $!\n";
    }elsif($e & 127) {
        printf "Child died with signal %d, %s coredump\n",
        ($e & 127),  ($e & 128) ? 'with' : 'without';
    }elsif($e != 0){
        printf "Child exited with value %d\n", $e >> 8;
    }
}

#################################################
sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;
    print <<EOT

Usage: $0 [options] script1.sh script2.sh [...scriptN.sh]

Scripts will be submitted to the queue in the order provided on the 
commandline. By default each script will wait for the previous script to exit 
before starting, but several scripts can be run simulataneously using the 
--processes option (see below). 

Options:
    
    -p INT --processes INT
        Number of scripts to run simultaneously. That is, if set to 2, script N 
        will wait for script N-2 to complete before starting. Default=1 (i.e. 
        scripts all run sequentially).

    -d, --dummy
        Do a dummy run without actually submitting scripts.

    -w job-IDS --wait_ids job-IDs
        One or more job IDs to wait for before executing ANY script, separated
        with commas.

    -h, --help
        Show this message and exit.

EOT
    ;
    exit 2 if $msg;
    exit;
}

