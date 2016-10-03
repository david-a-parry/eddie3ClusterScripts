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
    r => "4:00:00",
    f => \@file_regexes,
);
GetOptions
(
    \%opts,
    "s|site=s",
    "t|threads=i",
    "u|user=s",
    "d|directory=s",
    "f|files=s{,}",
    "p|prefix=s",
    "a|ascp=s",
    "q|qsub",
    "r|runtime=s",
    "tmux",
    "y|dry_run",
    "x|skip_completed",
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

    -f,--files PATTERN [PATTERN2 PATTERN3 ... ]
        One or more tring or REGEX pattern for retrieving matching files. 
        To retrieve bam files, their indexes and MD5s the REGEX 
        ".+\.ba[im](\.md5)*" could be used.
        Default = none (get whole directories).

 
    -t,--threads INT
        Number of threads to use
    
    -q,--qsub
        Use qsub to create and submit scripts
    
    -r,--runtime STRING
        Runtime (in the format HH::MM::SS). Default = 4:00:00

    -x,--skip_completed
        Skip completed scripts from previous runs (if exact same parameters used
        and using --qsub)

    -p,--prefix STRING
        Filename prefix for scripts. Default = 'get'. Each script will be named
        '[prefix].[no.].sh'.

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

if ($opts{t} < 1){
    $opts{t} = 1;
}
my $password = ''; 
if ($ENV{ASPERA_SCP_PASS}){
    $password =$ENV{ASPERA_SCP_PASS};
}else{ 
    $password = read_password
    (
        "Enter password for $opts{u}: "
    );
    $ENV{ASPERA_SCP_PASS} = $password;
} 

my @files = ();
my %done_dirs = ();
processDir($opts{d});
if (not @files){
    die "No matching files or directories identified - nothing to do!\n";
}
if (not -d "./$opts{d}"){
    make_path("./$opts{d}"); 
}
$opts{s} =~ s/https:\/\///;
if ( $opts{q} ){
    makeAndSubmitQsub();
}else{
    getWithParallel();
}

#####################################################################
sub makeAndSubmitQsub{
    my @wait_ids = (); 
    my $subbed_scripts = 0;
    for (my $i = 0; $i < @files; $i++){
        my $script = "$opts{p}.$i.sh";
        my $ascp_cmd = 
"$opts{a} -k 1 -P 33001 -O 33001 -l 500M dparry\@edgen-dt.rdf.ac.uk:$files[$i] ./$files[$i]";
        if ($opts{x}){
            if (completedPreviously($script, $files[$i], $ascp_cmd)){
                print STDERR "Skipping previously completed script/file: $script/$files[$i]\n";
                next;
            }
        }
        open (my $SCRIPT, ">", $script) or die "Can't open qsub script $script for writing: $!\n";
        my ($f, $d) = fileparse($files[$i]); 
        if (not -d "./$d"){
            make_path("./$d") or die "could not make path './$d': $!\n"; 
        }
        print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -cwd
#\$ -V
#\$ -l h_rt=$opts{r}
# Configure modules
. /etc/profile.d/modules.sh

$ascp_cmd

EOT
;
        close $SCRIPT;
        my $wait_string = ''; 
        if ($subbed_scripts >= $opts{t} and @wait_ids > $subbed_scripts - $opts{t}){
            $wait_string = "-hold_jid $wait_ids[$subbed_scripts - $opts{t}]";
        }
        $subbed_scripts++;
        my $cmd = "qsub $wait_string $script";
        if ($opts{y}){
            informUser("Dry run: $cmd\n");  
        }else{
            informUser("Executing: $cmd\n");  
            my $output = `$cmd`; 
            checkExit($?);
            if ($output =~ /Your job (\d+) .* has been submitted/){
                push @wait_ids, $1;
            }else{
                die "Error parsing qsub output for '$cmd'\nOutput was: $output";
            }
        }
    }
}

#####################################################################
sub completedPreviously{
    my $script = shift;
    my $file = shift;
    my $cmd = shift;
    my $output = "$script.stdout";
    return 0 if not -e $output;
    #IF OUTPUT FROM PREVIOUS RUN EXISTS...
    my $comp = 0;
    open (my $PREV, "<", $output) or die "Could not open $output: $!\n";
    while (my $l = <$PREV>){
        if ($l =~ /Completed:\s+\d+[KMGTP]*\s+bytes transferred in \d+ seconds/){
        #IF PREVIOUS RUN COMPLETED
            $comp = 1;
        }
    }
    close $PREV;
    if ($comp){
        #IF PREVIOUS RUN COMPLETED
        open (my $SCRIPT, "<", $script) or die "Could not open $script: $!\n";
        local $/ = undef;
        my $s = <$SCRIPT>;
        if ($s =~ /$cmd/){
            #AND WAS TRANSFERRING THE SAME FILE
            return 1;
        }
    }
    return 0;
}
        
#####################################################################
sub getWithParallel{
    my ($FH, $filename) = tempfile();
    my @parstring = "";
    foreach my $d (@files){
        print $FH  "$opts{u}\@$opts{s}:$d ./$opts{d}\n";
    }
    close $FH;
    my $cmd = "parallel ";
    $cmd .= "--tmux " if $opts{tmux};
    $cmd .= "--colsep ' ' -j $opts{t} $opts{a} -P 33001 -O 33001 -k 1 -l 500M :::: $filename ";
    if ($opts{y}){
        informUser("Dry run command:\n$cmd\n"); 
    }else{
        informUser("Attempting command:\n$cmd\n"); 
        system($cmd);
        checkExit($?);
    }

=cut
    for (my $i = 0; $i < @files;){
        for (my $j = 0; $j < $opts{t}; $j++){
            last if $i >= @dirs;
            my $thr = threads->create(\&getDir, \@dirs, $i++);
        }
        while (threads->list(threads::all)){
            foreach my $joinable (threads->list(threads::joinable)){
                $joinable->join();
                if ($i < @dirs){
                    my $thr = threads->create(\&getDir, \@dirs, $i++);
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
    $done_dirs{$d} = undef;
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
        if (@file_regexes){
            my $f_match = 0;
            foreach my $regex (@file_regexes){
                if ($l =~ /data-fasp-url\=.*\@\S+:\d+\/($regex)\?/){
                    push @d_files , "/$1";
                    $f_match = 1;
                }
            }
            if (not $f_match){
                if($l =~ /href=\"\/(\Q$opts{d}\E\/\S+)\"/){
                    push @sub_dirs, $1;
                }
            }
        }else{
            if($l =~ /href=\"\/(\Q$opts{d}\E\/\S+)\"/){
                push @sub_dirs, $1;
            }
        }
    }
    @sub_dirs  = grep {!exists $done_dirs{$_} } @sub_dirs;#remove parent directories
        
    if (@d_files){
        informUser
        (
            "INFO: Identified ". scalar(@d_files) ." files from '$opts{s}/$d\n"
        );
        push @files, @d_files;
    }
    if (@sub_dirs){
        informUser
        (
            "INFO: Identified ". scalar(@sub_dirs) ." subdirectories from ".
            "'$opts{s}/$d\n"
        );
        foreach my $sd (@sub_dirs){
            informUser
            (
                "INFO: Processing subdirectory $sd\n"
            );
            processDir($sd);
        }
    }else{
        if (not @file_regexes){
            push @files, $d;
        }
    }
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
