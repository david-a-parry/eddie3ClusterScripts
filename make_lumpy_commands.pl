#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;

die "Usage: $0 BAM [BAM ...]\n" if @ARGV < 1;
my $sdir = "lumpy_cmds_scripts";
my $pre_out = "lumpy_spliiters_discordants";
foreach my $dir ($sdir, $pre_out){
    if (not -d $dir){
        mkdir($dir) or die "Could not create $dir: $!\n";
    }
}
my @scripts;
my @bams;
my @disc;
my @splt;
foreach my $bam (@ARGV){
    make_pre_script($bam);
}

my $lscript = "$sdir/lumpy_express.sh";
my $bamstring = join(",", @bams);
my $discstring = join(",", @disc);
my $splitstring = join(",", @splt);
open (my $LSCRIPT, ">", $lscript) or die "Could not create $lscript: $!\n";
print $LSCRIPT <<EOT
#!/bin/bash
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -V
#\$ -cwd
#\$ -l h_rt=48:00:00
#\$ -l h_vmem=32G
#\$ -e $lscript.stderr
#\$ -o $lscript.stdout
# Configure modules
. /etc/profile.d/modules.sh
#Load modules
module load igmm/libs/htslib/1.4
module load igmm/apps/sambamba/0.5.9
module load igmm/libs/ncurses/6.0
module load igmm/apps/samtools/1.4
module load roslin/bedtools/2.26.0
module load igmm/apps/python/2.7.10 #need to load python explicitly for pysam and numpy

LUMPY=/gpfs/igmmfs01/eddie/igmm_datastore_MND-WGS/lumpy-sv/bin/lumpyexpress

\$LUMPY -B $bamstring -D $discstring -S $splitstring -o var.lumpy.vcf

EOT
;
close $LSCRIPT;
my @wait = ();
foreach my $sc (@scripts){
    push @wait, do_qsub("qsub $sc");
}
do_qsub("qsub -hold_jid ". join(",", @wait) . " $lscript");



##################################################################################
sub do_qsub{
    my $cmd = shift;
    print STDERR "Executing: $cmd\n";
    my $stdout = `$cmd`;
    check_exit($?);
    if ($stdout =~ /Your job (\d+) .* has been submitted/){
        return $1;
    }else{
        die "Error parsing qsub stdout for '$cmd'\nOutput was: $stdout";
    }
}

##################################################################################
sub make_pre_script{
    my $bam = shift;
    my $base = basename($bam, ('.cram', '.bam'));
    my $pre_script = "$sdir/preprocess.$base.sh";
    open (my $SCRIPT, ">", $pre_script) or 
        die "Could not open $pre_script for writing: $!\n";
    my $d = "$pre_out/$base.discordants.bam";
    my $s = "$pre_out/$base.splitters.bam";
    print $SCRIPT <<EOT
#!/bin/bash
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -V
#\$ -cwd
#\$ -l h_rt=48:00:00
#\$ -l h_vmem=4G
#\$ -e $pre_script.stderr
#\$ -o $pre_script.stdout
# Configure modules
. /etc/profile.d/modules.sh
#Load modules
module load igmm/libs/htslib/1.4
module load igmm/libs/ncurses/6.0
module load igmm/apps/samtools/1.4
module load roslin/bedtools/2.26.0
module load igmm/apps/python/2.7.10 #need to load python explicitly for pysam and numpy

LUMPY_HOME=/gpfs/igmmfs01/eddie/igmm_datastore_MND-WGS/lumpy-sv

samtools view -b -F 1294 $bam > $d

samtools view -h $bam | \$LUMPY_HOME/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > $s

EOT
    ;
    close $SCRIPT;
    push @disc, $d;
    push @splt, $s;
    push @bams, $bam;
    push @scripts, $pre_script;
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
