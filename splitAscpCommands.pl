#!/usr/bin/perl
use warnings;
use strict;
#use threads;
use POSIX qw/strftime/;
use Term::ReadPassword;
use File::Temp qw/ tempfile / ;
use File::Path qw( make_path );
use File::Temp qw/ tempfile /;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
my @file_regexes = ();
my %opts = 
(
    s => "https://edgen-dt.rdf.ac.uk",
    u => "dparry",
    d => "",
    a => "$ENV{HOME}/.aspera/connect/bin/ascp",
    t => 8,
    p => "get",
    r => 24,
);
GetOptions
(
    \%opts,
    "s|site=s",
    "t|threads=i",
    "u|user=s",
    "d|directory=s",
    "p|prefix=s",
    "a|ascp=s",
    "q|qsub",
    "r|runtime=i",
    "tmux",
    "h|help",
) or usage("Syntax error");
usage() if $opts{h};
sub usage{
    my $msg = shift;
    print STDERR "\nERROR: $msg\n" if $msg;
    print STDERR <<EOT 

Retrieve files listed at Edinburgh Genomics ASCP site using multiple streams.

Usage: $0 -d <target directory>

Options:

    -s,--site URL
        URL of website to search. Default = https://edgen-dt.rdf.ac.uk

     u,--user STRING
        Username for site - default = dparry

    -d,--directory STRING
        Name of directory to retrieve at the site. Required.
        The directory structure will be recreated from the current working 
        directory if not already present.

    -t,--threads INT
        Number of threads to use
    
    -q,--qsub
        Use qsub to create and submit scripts

    -a,--ascp FILE
        Location of ascp binary. Default = ~/.aspera/connect/bin/ascp
    
    --tmux
        Create tmux terminals for each command

    -h,--help
        Show this message and exit

EOT
;

    exit 1 if $msg;
    exit;
}

usage("-d/--directory option is required") if not $opts{d};

my $password = read_password
(
    "Enter password for $opts{u}: "
); 

my @sub_dirs = processDir($opts{d});
if (not -d "./$opts{d}"){
    make_path("./$opts{d}"); 
}
$ENV{ASPERA_SCP_PASS} = $password;
$opts{s} =~ s/https:\/\///;
if ( $opts{q} ){
    makeAndSubmitQsub();
}else{
    getWithParallel();
}

#####################################################################
sub makeAndSubmitQsub{
    for (my $i = 0; $i < @sub_dirs; $i++){
        my $script = "$opts{p}.$i.sh";
        open (my $SCRIPT, ">", $script) or die "Can't open qsub script $script for writing: $!\n";
        print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -cwd
#\$ -V
#\$ -l h_rt=$opts{r}:00:00
# Configure modules
. /etc/profile.d/modules.sh

$opts{a} -k 1 -P 33001 -O 33001 -l 500M dparry\@edgen-dt.rdf.ac.uk:$sub_dirs[$i] ./$opts{d}

find $sub_dirs[$i] -name '*md5' -exec md5sum -c {} >> $sub_dirs[$i]/md5_checks.txt \\+

EOT
;
    close $SCRIPT;
    system("qsub $script"); 
    checkExit($?);
    }
}
#####################################################################
sub getWithParallel{
    my ($FH, $filename) = tempfile();
    my @parstring = "";
    foreach my $d (@sub_dirs){
        print $FH  "$opts{u}\@$opts{s}:$d ./$opts{d}\n";
    }
    close $FH;
    my $cmd = "parallel ";
    $cmd .= "--tmux " if $opts{tmux};
    $cmd .= "--colsep ' ' -j $opts{t} $opts{a} -P 33001 -O 33001 -k 1 -l 500M :::: $filename ";
    informUser("Attempting command:\n$cmd\n"); 
    system($cmd);
    checkExit($?);

=cut
    for (my $i = 0; $i < @sub_dirs;){
        for (my $j = 0; $j < $opts{t}; $j++){
            last if $i >= @sub_dirs;
            my $thr = threads->create(\&getDir, \@sub_dirs, $i++);
        }
        while (threads->list(threads::all)){
            foreach my $joinable (threads->list(threads::joinable)){
                $joinable->join();
                if ($i < @sub_dirs){
                    my $thr = threads->create(\&getDir, \@sub_dirs, $i++);
                }
            }
            sleep 10;
        }
    }
=cut
}

#####################################################################
sub getDir{
    my $dirs = shift;
    my $i = shift;
	print STDERR "\n>Getting dir " . (1 + $i) ." of ". scalar(@$dirs) . ".\n";
    my $cmd = "$opts{a} -P 33001 -O 33001 -k 1 -l 500M $opts{u}\@$opts{s}:$dirs->[$i] ./";
	print STDERR ">>$cmd\n";
    system("$cmd");
    checkExit($?);
}

#####################################################################
sub processDir{
    my $d = shift;
    #curl will list parent directories as well so we need to keep track 
    # and ignore already processed directories
    my @d_files = ();
    my @sub_dirs = (); 
    informUser("Attempting to get file listing for directory '$d' with curl...\n");
    my $f_http = 
      `curl  -s $opts{s}/$d --list-only -k  --user $opts{u}:$password`;
    die "Error executing curl command" unless checkExit($?);
    if ($f_http =~ /401 Authorization Required/){
        die "Authorization failed - did you enter the correct password?\n";
    }
    foreach my $l (split("\n", $f_http)){
        if($l =~ /href=\"\/(\Q$opts{d}\E\/\S+)\"/){
            push @sub_dirs, $1;
        }
    }
    @sub_dirs  = grep {$_ ne $d } @sub_dirs;#remove parent directories
    if (not @sub_dirs){
        informUser
        (
            "WARNING: Could not identify any subdirectories ".
            "from '$opts{s}/$d'\n"
        );
        die "Nothing to do - exiting.\n";
    }elsif (@sub_dirs){
        informUser
        (
            "INFO: Identified ". scalar(@sub_dirs) ." subdirectories from ".
            "'$opts{s}/$d\n"
        );
    }
    return @sub_dirs;
}

#####################################################################
sub checkExit{
    my $er = shift; 
    if ($er){
        if ($er == -1){
            informUser
            (
                "ERROR: command failed to execute.\n"
            ) and return;
        }
        my $exit = $er >> 8;
        informUser
        (
            "ERROR: command failed code $exit: $er.\n"
        ) and return;
    }
    return 1;
}
#################################################
sub informUser{
    my $msg = shift;
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] $msg";
}
