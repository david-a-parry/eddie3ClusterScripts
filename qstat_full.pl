#!/usr/bin/perl

use strict;
use warnings; 

my @qstat = split("\n", `qstat -r`); 
exit if not @qstat;
my $job_no ;
my $job_name ;
my %job = ();
my @columns = ( "job-ID",  "prior", "name", "user", "state", "submit/start at", "queue", "slots", "ja-task-ID", ) ;

#cols are different lengths separated by spaces
#in order to detect empty cols we need to know where to look
my %indices = getColIndices($qstat[0]);

my @info = ("job-ID", "prior", "name", "user", "state", "submit/start at", "queue", "slots", "fullname" ) ;
my @widths = map { length($info[$_]) } 0..$#info;
my @job_strings = ();
foreach my $q (@qstat){
    if ($q =~ /^(\d+)/){
        getJobInfo();
        %job = ();
        getColumnValues($q);
        next;
    }
    if ($q =~ /^\s+Full jobname:\s+(\S+)/){
        $job{fullname} = $1;
    }
}

getJobInfo();

outputJobInfo();

##################################################
sub getColumnValues{
    my $q = shift; 
    foreach my $c (@columns){
        my $s = substr($q, $indices{$c}->{start}, $indices{$c}->{end} - $indices{$c}->{start}); 
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
        printf $format, @$j;
    }
}

##################################################
sub getJobInfo{
    return if not %job;
    my @output = map { $job{$_} || " " } @info ; 
    @widths    = map { length($output[$_]) > $widths[$_] ? length($output[$_]) : $widths[$_] } 0..$#info;
    #print join("\t", @output) . "\n";
    push @job_strings, \@output;
}
