#!/usr/bin/perl
use warnings;
use strict;
use POSIX qw/strftime/;
use Term::ReadPassword;
use File::Temp qw/ tempfile / ;
use Getopt::Long;
use Data::Dumper;
my @file_regexes = ();
my %opts = 
(
    s => "https://edgen-dt.rdf.ac.uk",
    u => "dparry",
    d => "",
    a => "$ENV{HOME}/.aspera/connect/bin/ascp",
    f => \@file_regexes,
);
GetOptions
(
    \%opts,
    "s|site=s",
    "u|user=s",
    "d|directory=s",
    "a|ascp=s",
    "g|get=s",
    "f|files=s{,}",
    "x|do_not_expand_path",
    "h|help",
) or usage("Syntax error");
usage() if $opts{h};
my $default_regex = ".+\\.g\\.vcf\\.gz(\\.md5|\\.tbi|\\.tbi\\.md5)*";
@file_regexes = ($default_regex) if not @file_regexes;
sub usage{
    my $msg = shift;
    print STDERR "\nERROR: $msg\n" if $msg;
    print STDERR <<EOT 

List or retrieve files listed at Edinburgh Genomics ASCP site.

Usage: $0 -d <target directory>

Options:

    -s,--site URL
        URL of website to search. Default = https://edgen-dt.rdf.ac.uk
    
    -u,--user STRING
        Username for site - default = dparry

    -d,--directory STRING
        Name of directory to search at the site. Required.
    
    -f,--files PATTERN [PATTERN2 PATTERN3 ... ]
        One or more tring or REGEX pattern for retrieving matching files. 
        The default (see below) is to retrieve GVCF files, indexes and their
        MD5 checksum files. To retrieve bam files and their indexes the REGEX 
        ".+\.ba[im](\.md5)*" could be used.
        Default = "$default_regex"

    -a,--ascp FILE
        Location of ascp binary. Default = ~/.aspera/connect/bin/ascp

    -g,--get DIR
        If specified, files will be retrieved using ascp to this directory.
        Will be created if it does not already exist. Existing files will be
        overwritten.
    
    -x,--do_not_expand_path
        Use this flag to place all files in the directory specified by the 
        '--get' argument instead of creating paths to reflect the remote 
        directory structure. May result in overwriting files with the same 
        filename.

    -h,--help
        Show this message and exit
        

EOT
;
    exit 1 if $msg;
    exit;
}

usage("-d/--directory option is required") if not $opts{d};

if($opts{g} and not -d $opts{g}){
    informUser("Attempting to create output directory '$opts{g}'...\n");
    mkdir($opts{g}) or die "Could not create ouptut directory '$opts{g}': $!\n";
    informUser("Success\n");
}

my @dirs =();
my @files = ();

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


my %done_dirs = ();
processDir($opts{d}); 

if (not $opts{g}){
    print join("\n", @files) . "\n";
}else{
    my ($TEMP, $temp) = tempfile();
    if ($opts{x}){
        print $TEMP join("\n", @files) . "\n"; 
    }else{
        print $TEMP map {"$_\n$opts{g}/$_\n"} @files;
    }
    close $TEMP;
    informUser("Attempting to retrieve target files using ascp...\n");
    my $f_arg = $opts{x} ? "--file-list" : "--file-pair-list";
    $opts{s} =~ s/^http(s){0,1}:\/\///;
    system("$opts{a} -P 33001 -O 33001 -l 500M --host $opts{s}  --user $opts{u} $f_arg $temp  --mode=recv $opts{g}");
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
    checkExit($?);
    if ($f_http =~ /401 Authorization Required/){
        die "Authorization failed - did you enter the correct password?\n";
    }
    foreach my $l (split("\n", $f_http)){
        #if ($l =~ /data-fasp-url\=.*\@(.+.g.vcf.gz(\.md5|\.tbi|\.tbi\.md5)*)\?/){
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
    }
    @sub_dirs  = grep {!exists $done_dirs{$_} } @sub_dirs;#remove parent directories
    if (not @d_files and not @sub_dirs){
        informUser
        (
            "WARNING: Could not identify any target files or subdirectories ".
            "from '$opts{s}/$d'\n"
        );
    }else{
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
}
#################################################
sub informUser{
    my $msg = shift;
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] $msg";
}
