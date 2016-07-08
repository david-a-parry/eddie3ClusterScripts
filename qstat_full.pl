#!/usr/bin/perl

use strict;
use warnings; 

my @qstat = split("\n", `qstat -r`); 

my $job_no ;
my $job_name ;
my %job = ();
my @columns = qw ( job-ID  prior name user state submit/start at queue slots ja-task-ID ) ;
my @info = qw ( job-ID state submit/start at queue slots fullname ) ; 
print "#" . join("\t", @info) . "\n";
foreach my $q (@qstat){
    if ($q =~ /^(\d+)/){
        outputJobInfo();
        %job = ();
        @job{@columns} =  split(/\s+/, $q);
        next;
    }
    if ($q =~ /^\s+Full jobname:\s+(\S+)/){
        $job{fullname} = $1;
    }
}

outputJobInfo();


sub outputJobInfo{
    return if not %job;
    { 
        no warnings 'uninitialized';
        my @output = map { $job{$_} } @info ; 
        print join("\t", @output) . "\n";
    }
}
