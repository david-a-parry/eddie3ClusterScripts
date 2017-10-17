#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

my @modules = ();
my @mail_opt = ();
my %opts = 
(
    r   => "1:00:00",
    m   => "1G",
    E   => "david.parry\@igmm.ed.ac.uk",
    M   => \@mail_opt,
    l   => \@modules,
);

GetOptions
(
    \%opts,
    "r|runtime=s",
    "m|memory=s",
    "p|parallel_environment=s",
    "t|threads=i",
    "e|stderr=s",
    "o|stdout=s",
    "E|email=s",
    "M|mail=s{,}",
    "N|name=s",
    "l|load_modules=s{,}",
    "no_cwd",
    "no_V",
    "h|help",
) or usage( "Syntax error" ); 
usage() if $opts{h};

my %head_opts = 
(
    r => "-l h_rt=",
    m => "-l h_vmem=",
    p => "-pe ",
    t => "-pe sharedmem ",
    E => "-M ",
    M => "-m ",
    N => "-N ",
    e => "-e ",
    o => "-o ",
    l => "module load ",
);

my @header = ();
push @header, "#\$ -cwd" unless $opts{no_cwd};
push @header, "#\$ -V" unless $opts{no_V};

push @header,
  map  { "#\$ $_" } #prepend header identifier
  map  { $head_opts{$_} . $opts{$_} } #prepend option identifier
  grep { $_ !~ /^(l|M|no_cwd|no_V)$/ } #don't add these args yet
  keys %opts 
;
push @header,
    map  { "#\$ $_" } #prepend header identifier
    map  { $head_opts{M} . $_ } @mail_opt
;

#sort so order is consistent ( before adding modules, which must be at bottom )
@header = sort @header; 
push @header , "# Configure modules\n. /etc/profile.d/modules.sh\n#Load modules";
;

push @header, map  { $head_opts{l} . $_ } @modules;

print join("\n", @header) . "\n";

##################################################
sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;
    print STDERR <<EOT

Generates a generic qsub header

Usage: $0 [options]

Options: 

    -r STRING --runtime STRING
         Value for runtime argument. Default = 1:00:00
        
    -m STRING --memory STRING
         Value for memory argument. Default = 1G
        
    -p STRING --parallel_environment STRING
         Value for parallel environment setting (e.g. 'sharedmem 8'). Default = none.
        
    -t INT --threads INT
         Number of 'sharedmem' threads (e.g. -t 8 is a shortcut for '-p sharedmem 8'). Default = none.

    -e STRING --stderr STRING
         Output script STDERR to this file. Default = none
        
    -o STRING --stdout STRING
         Output script STDOUT to this file. Default = none
        
    -E STRING --email STRING
        Use this email address for notifications. Default = david.parry\@igmm.ed.ac.uk

    -M --mail STRING(s) 
        Value for script email option (e.g. a, b or e). Default = none
        
    -l STRING(s) --load_modules STRING(s)
        Modules to load. Default = none
    
    -no_cwd
        By default the -cwd option is added to your header. Use this flag to remove it.
    
    -no_V
        By default the -V option is added to your header. Use this flag to remove it.
    
    -h --help
        Show this message and exit

    Examples: 

    $0
    #generate basic header
    
    $0 -m 8G -r 4:00:00 
    #use 8 GB of memory, 4 hour runtime

    $0 -m 8G -r 4:00:00 -E foo\@bar.com -M a b e  
    #as above specifying email address and to email at beginning, end and if aborted
    
    $0 -m 8G -r 4:00:00 -E foo\@bar.com -M a b e -l /apps/igmm/something /apps/another/module
    #as above but also specifying modules

EOT
    ;
    exit 1 if $msg;
    exit;
}

