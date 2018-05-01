#!/usr/bin/perl

use strict;
use warnings; 
use Getopt::Long;
Getopt::Long::Configure ("bundling");
my %opts = (d => 30);
GetOptions
(
    \%opts,
    'h|help',
    'r|running',
    'q|queued',
    'i|held',
    'e|errored',
    'w|wait_id=i',
    'd|delay=i',
) or usage("syntax error!");
usage() if $opts{h};


my $job_no ;
my $job_name ;
my %job = ();
my @columns = ("job-ID",  "prior", "name", "user", "state", "submit/start at", 
               "queue", "slots", "ja-task-ID");
my %indices = ();

#cols are different lengths separated by spaces
#in order to detect empty cols we need to know where to look

my @info = ("job-ID", "prior", "name", "user", "state", "submit/start at",
            "queue", "slots", "ja-task-ID", "fullname");
my @widths = map { length($info[$_]) } 0..$#info;
my %current_jobs = (); 
my @job_strings = ();
my @qstat = split("\n", `qstat -r`); 
processQstat();

getJobInfo();

outputJobInfo();

if ($opts{w}){
    while (exists $current_jobs{$opts{w}}){
        sleep($opts{d});
        @qstat = split("\n", `qstat -r`); 
        processQstat();
        getJobInfo();
        outputJobInfo();
    }
}

##################################################
sub processQstat{   
    exit if not @qstat;
    %indices = getColIndices($qstat[0]);
    foreach my $q (@qstat){
        if ($q =~ /^\s*(\d+)/){
            getJobInfo();
            %job = ();
            getColumnValues($q);
            next;
        }
        if ($q =~ /^\s+Full jobname:\s+(\S+)/){
            $job{fullname} = $1;
        }
    }
}

##################################################
sub getColumnValues{
    my $q = shift; 
    foreach my $c (@columns){
        my $s = substr($q, 
                       $indices{$c}->{start}, 
                       $indices{$c}->{end} - $indices{$c}->{start}); 
        $job{$c} = $s;
    }
}
        
##################################################
sub getColIndices{
    my $head = shift;
    my %indices = (); 
    foreach my $c (@columns){
        if ($head =~ /$c\s+/){
            $indices{$c}->{start} = $-[0];
            $indices{$c}->{end} = $+[0];
        }else{
            die "Could not read header!\n";
        }
    }
    return %indices; 
}
##################################################
sub outputJobInfo{
    return if not @job_strings;
    my $format = join(' ', map { "%-${_}s" } @widths) . "\n";
    my $total_length = 0;
    map { $total_length += $widths[$_] } 0..$#widths;
    printf $format, @info;
    print "-" x $total_length . "-" x $#info . "\n";
    foreach my $j (@job_strings){
        if ($opts{r} or $opts{q} or $opts{i} or $opts{e}){
            if ( ($opts{r} and $j->[4] =~ /r/) or
                 ($opts{q} and $j->[4] =~ /^\s*qw\s*$/) or 
                 ($opts{i} and $j->[4] =~ /^\s*hqw\s*$/) or 
                 ($opts{e} and $j->[4] =~ /E/ )
            ){
                printf $format, @$j;
            }
        }else{
            printf $format, @$j;
        }
    }
}

##################################################
sub getJobInfo{
    return if not %job;
    my @output = map { $job{$_} || " " } @info ; 
    @widths    = map { length($output[$_]) > $widths[$_] ? 
                       length($output[$_]) : $widths[$_] } 0..$#info;
    #print join("\t", @output) . "\n";
    push @job_strings, \@output;
    (my $j = $job{'job-ID'}) =~ s/\s//g;
    $current_jobs{$j} = \%job;
}


##################################################
sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;
    print STDERR <<EOT

    qstat_full.pl - prints qstat output including full name of job in columns

    Usage:
            qstat_full.pl [options]

    Options:

        -h,--help
            Show this message and exit

        -r,--running
            Show only jobs that are running

        -q,--queued
            Show only queued waiting jobs (qw)

        -i,--held
            Show only held jobs (hqw)

        -e,--errored
            Show only jobs in error state

        -w job-ID --wait_id job-ID
            Continue running while the given job-ID is running. The numeric
            job-ID must be given, not the job name.

EOT
    ;
    exit 1 if $msg;
    exit;
}
    
