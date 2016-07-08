 #!/usr/env perl 
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;

my %samples = ();
my @scripts = (); 
usage() if (@ARGV != 2); 

my $genome_folder = shift;
my $dir = shift;
opendir (my $DIR , $dir) or die "Cannot open directory $dir: $!\n";
my @files = readdir($DIR);
foreach my $f (@files){
    if ($f =~ /(.*)_[12].fastq(\.gz)*$/){
        push @{$samples{$1}}, $f;
    }else{
        print STDERR "Skipping file $f - does not look like fastq.\n";
    }
}

foreach my $s (keys %samples){
    if (@{$samples{$s}} > 2){
        die "More than 2 fastqs found for sample $s!\n";
    }
    my $fq_string; 
    my $sam = $s."_bismark.sam";
    if (@{$samples{$s}} == 1){
        $fq_string = $samples{$s}->[0];
    }else{
        my ($fastq1, $fastq2) = sort ( @{$samples{$s}} ); 
        $fq_string = "-1 $fastq1 -2 $fastq2";
        $sam = $s."_bismark_pe.sam";
    }
    my $script = $s . "_bismark_script.qsub";
    open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -cwd
#\$ -e $script.stderr
#\$ -o $script.stdout
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/bowtie/1.0.0
module load apps/perl/bismark/0.10.1
bismark $genome_folder $fq_string
bismark_methylation_extractor $sam 
EOT
;
    close $SCRIPT;
    push @scripts, $script; 
}



##################################################################################
sub usage{
    my $msg = shift; 
    print "ERROR: $msg\n" if $msg;
    print <<EOT
    
    USAGE: perl $0 path_to_genome_folder fastq_directory
    
    From a list of FASTQs this script will create a seies of qsub scripts to perform alignment and post-processing.

EOT
;
    exit 1 if $msg;
    exit;
}
