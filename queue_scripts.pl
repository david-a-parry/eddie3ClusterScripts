#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my %opts = (p => 1);
GetOptions(
    \%opts,
    "p|processes=i",
    "d|dummy",
    "h|help",
) or usage("Error in options");
usage() if $opts{h};
usage("At least two scripts are required") if @ARGV < 2;
my $n = 0;
my @holds = ();
while (my $script = shift){
    my $cmd = "qsub ";
    my $hidx = $n - $opts{p};
    if (@holds > $hidx and $hidx >= 0){
        $cmd .= "-hold_jid $holds[$n - $opts{p}] ";
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
            die "Error parsing qsub output for $cmd!\nOutput was: $output";
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

Options:
    
    -p INT --processes INT
        Number of scripts to run in at simultaneously

    -d, --dummy
        Do a dummy run without  actually submitting scripts

    -h, --help
        Show this message and exit.

EOT
    ;
    exit 2 if $msg;
    exit;
}

