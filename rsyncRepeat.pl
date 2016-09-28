#!/usr/bin/env perl 
use warnings;
use strict;
use POSIX qw/strftime/;
use Getopt::Long; 

my %opts = ();

GetOptions(
    \%opts,
    'n|number_retries=i',
    'h|help',
) or die "Syntax error";

usage() if $opts{h};
usage("No arguments provided!") if not @ARGV;

my $arg_string = join(" ",  @ARGV); 
my $cmd = "rsync $arg_string";
print STDERR "Using rsync commandline: '$cmd'\n";
my $n = 0;
while (1){ 
    my $time = strftime( "%H:%M:%S", localtime );
    $n++; 
    my $info = "Attempt $n";
    if ($opts{n}){
        $info .= " of $opts{n}";
    }
    print STDERR "$info...\n";
    system($cmd);
    check_exit($?); 
    if ($? == 0){
        print "rsync completed normally\n";
        exit;
    }else{
        if ($opts{n} and $n > $opts{n}){
            die "Reached maximum number of retries ($opts{n}).\n";
        }
        print "Rsync failure. Backing off and retrying...\n";
        sleep 30;
    }
}

##################################################################################
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

##################################################################################
sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;
    print STDERR <<EOT

Usage: $0 [script options] [rsync options/arguments]

Options:
    
    -n INT, --number_retries INT
        Set a limit for number of retries.

    -h,--help
        Show this message and exit.

Example:
    
    $0 /my/source/dir /my/dest/dir
    $0 "-azv /my/source/dir /my/dest/dir/"
    $0 -n 20 "-azv /my/source/dir /my/dest/dir/"

EOT
    ;
    exit 1 if $msg;
    exit;
}


