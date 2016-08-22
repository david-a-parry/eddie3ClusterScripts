#!/usr/bin/perl

use strict;
use warnings; 

my @qstat = split("\n", `qstat -r`); 

my $job_no ;
my $job_name ;
my %job = ();
my @columns = qw ( job-ID  prior name user state submit/start at queue slots ja-task-ID ) ;
my @info = qw ( job-ID state submit/start at queue slots fullname ) ; 
my @widths = map { length($info[$_]) } 0..$#info;
my @job_strings = ();
foreach my $q (@qstat){
    if ($q =~ /^(\d+)/){
        getJobInfo();
        %job = ();
        @job{@columns} =  split(/\s+/, $q);
        next;
    }
    if ($q =~ /^\s+Full jobname:\s+(\S+)/){
        $job{fullname} = $1;
    }
}

getJobInfo();

outputJobInfo();


##################################################
sub outputJobInfo{
    my $format = join(' ', map { "%-${_}s" } @widths) . "\n";
    my $total_length = 0;
    map { $total_length += $widths[$_] } 0..$#widths;
    printf $format, @info;
    print "-" x $total_length . "-" x @info . "\n";
    foreach my $j (@job_strings){
        { 
            no warnings 'uninitialized';
            printf $format, @$j;
        }
    }
}

##################################################
sub getJobInfo{
    return if not %job;
    { 
        no warnings 'uninitialized';
        my @output = map { $job{$_} } @info ; 
        @widths    = map { length($output[$_]) > $widths[$_] ? length($output[$_]) : $widths[$_] } 0..$#info;
        #print join("\t", @output) . "\n";
        push @job_strings, \@output;
    }
}
