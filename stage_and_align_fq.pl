#!/usr/bin/env perl
use warnings;
use strict;
use FindBin qw($RealBin);
use File::Basename;
use Getopt::Long;
use File::Spec;
use POSIX;
use Data::Dumper;

my $tmpdir = File::Spec->tmpdir();
my $gatk = "$RealBin/GATK/GenomeAnalysisTK.jar";
my $picard = "$RealBin/picard-tools-2.9.0/picard.jar";
my %opts = (d => \$tmpdir, r => 48, t => 1, g => \$gatk, p => \$picard,
            b => 10, o => "/gpfs/igmmfs01/eddie/igmm_datastore_MND-WGS/GRCh37",
            );
GetOptions(
    \%opts,
    'h|?|help',
    'o|output_dir=s',
    'b|batch_size=i',
    'n|no_exec',
    't|threads=i',
    'g|gatk=s',
    'p|picard=s',
    'd|tmp_dir=s',
    'r|runtime=i',
) or usage("Error in options");
usage() if $opts{h};
usage("-o/--output_dir is required") if not $opts{o};
my $fq_list = shift;
my %fq = read_fq_list($fq_list);
die "No fastqs found\n" if not %fq;
my $refdir = "$RealBin/ref";
my $human_fasta = "$refdir/GRCh37/hs37d5.fa";
my $sub_dir = "$RealBin/subscripts";
my $bam_dir = "$RealBin/alignments";
my $vcf_dir = "$RealBin/vcf";
my $gvcf_dir = "$RealBin/gvcf";
my $fq_dir = "$RealBin/fastq";
foreach my $d ($bam_dir, $fq_dir, $sub_dir, $gvcf_dir, $vcf_dir){
    make_dir($d);
}
my %cmds = ();
my $dummy_wait = 0;
get_per_sample_cmds();
exec_cmds();

##################################################
sub get_per_sample_cmds{
    foreach my $s (keys %fq){
        push @{$cmds{$s}}, retrieve_files($s);
        push @{$cmds{$s}}, align_data($s);
        push @{$cmds{$s}}, make_gvcf($s);
        push @{$cmds{$s}}, stage_out($s);
    }
}

##################################################
sub exec_cmds{
    my @prev_waits = ();
    foreach my $k (keys %cmds){
        my $wait = get_wait(\@prev_waits);
        foreach my $script (@{$cmds{$k}}){
            $wait = do_qsub($script, $wait);
        }
        push @prev_waits, $wait;
    }
}

##################################################
sub get_wait{
    my $queue = shift;
    if (@$queue >= $opts{b}){
        return $queue->[-1 * $opts{b}];
    }
    return '';
}

##################################################
sub do_qsub{
    my ($script, $wait_id) = @_;
    my $cmd = "qsub $script";
    if ($wait_id){
        $cmd = "qsub -hold_jid $wait_id $script";
    }
    print STDERR "EXECUTING: $cmd\n";
    return ++$dummy_wait if $opts{n};
    my $stdout = `$cmd`; 
    if ($stdout =~ /Your job (\d+) .* has been submitted/){
        return $1;
    }else{
        die "Error parsing qsub stdout for '$cmd'\nOutput was: $stdout";
    }
}


##################################################
sub stage_out{
#copy aligned and processed BAM plus GVCF to --stage_out_dir
    my $s = shift;
    my $bam = "$bam_dir/$s.rmdups.indelrealign.bqsr.bam";
    my $gvcf = "$gvcf_dir/$s.g.vcf.gz";
    my $script = "$sub_dir/stage_out_$s.sh";
    open (my $STAGE, ">", $script) or die "Could not write to $script: $!\n";
    print $STAGE <<EOT
#!/bin/bash
#
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -cwd
# Choose the staging environment
#\$ -q staging

# Hard runtime limit
#\$ -l h_rt=12:00:00 

# Make the job resubmit itself if it runs out of time: rsync will start where it left off
#\$ -r yes
#\$ -notify
trap 'exit 99' sigusr1 sigusr2 sigterm

# Perform copy with rsync
rsync -vrlp --remove-source-files $bam $gvcf $opts{o}

EOT
    ;
    close $STAGE;
    return $script;
}
    
##################################################
sub retrieve_files{
    my $s = shift;
    my $stage_in = "$sub_dir/stage_in_$s.sh";
    my $get = join(" ", map {$fq{$s}->{$_}} sort keys %{$fq{$s}});
    open (my $STAGE, ">", $stage_in) or die "Could not create $stage_in: $!\n";
    print $STAGE <<EOT
#!/bin/bash
#
#\$ -e $stage_in.stderr
#\$ -o $stage_in.stdout
#\$ -cwd
#\$ -q staging
#\$ -l h_rt=12:00:00 
# Make the job resubmit itself if it runs out of time: rsync will start where it left off
#\$ -r yes
#\$ -notify
trap 'exit 99' sigusr1 sigusr2 sigterm

# Perform copy with rsync
rsync --ignore-existing -vr $get $fq_dir
chmod -R 700 $fq_dir
EOT
    ;
    close $STAGE;
    return $stage_in;
}

 
##################################################
sub read_fq_list{
    my $list = shift;
    open (my $IN, "<", $list) or die "Could not read fastq list $list: $!\n";
    while (my $l = <$IN>){
        my ($f) = split(/\s/, $l);
        my $name = fileparse($f);
        #files should be named [sample_name]_R1.fastq.gz
        if ($name =~ /(\S+)_R([12])\.(fastq|fq)(.gz)?$/){
            my $s = $1;
            my $r = $2;
            if (exists $fq{$s}->{$r}){
                die "ERROR: Duplicate sample/read ($s/$r) for fastq $f\n";
            }
            $fq{$s}->{$r} = $f;
        }
    }
    return  %fq;
}

##################################################
sub make_gvcf{
    my ($s) = @_;
    my $bam = "$bam_dir/$s.rmdups.indelrealign.bqsr.bam";
    my $gvcf = "$gvcf_dir/$s.g.vcf.gz";
    my $gt_cmd = <<EOT
java -Xmx4g -jar $gatk \\
-R $human_fasta \\
-T HaplotypeCaller \\
-D $RealBin/ref/GRCh37/dbSNP147/All_20160601.vcf.gz \\
-I $bam \\
-o $gvcf \\
--emitRefConfidence GVCF 
EOT
    ;
    my $script = "$sub_dir/gvcf.$s.sh";
    open (my $SCRIPT, ">", $script) or die "Could not create $script: $!\n";
    print $SCRIPT <<EOT
#!/bin/bash
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -cwd
#\$ -V
#\$ -l h_vmem=12G
#\$ -l h_rt=$opts{r}:00:00 
. /etc/profile.d/modules.sh 

$gt_cmd

EOT
    ;
    close $SCRIPT;
    return $script;
}

##################################################
sub align_data{
    #align to the human GRCh37 genome using standard bwa/GATK pipeline
    my ($s) = @_;
    my $bam_base = "$bam_dir/$s";
    my $fqs = join(" \\\n", map {$fq_dir . "/" . basename($fq{$s}->{$_})} 
                   sort keys %{$fq{$s}});
    my @cmds = ();
    push @cmds, <<EOT
ulimit -n 2048
echo "ulimit file limit set to " \$(ulimit -n)

bwa mem -t $opts{t} \\
-R '\@RG\\tID:$s\\tSM:$s\\tLB:$s\\tPL:ILLUMINA' \\
"$human_fasta" \\
$fqs \\
 | samtools view -Sb - > $bam_base.unsorted.bam

rm $fqs
EOT
    ;
    push @cmds, <<EOT
ulimit -n 2048
echo "ulimit file limit set to " \$(ulimit -n)

samtools sort -@ $opts{t} -O bam -T $bam_base.tmp_sort -o $bam_base.bam $bam_base.unsorted.bam
rm $bam_base.unsorted.bam
EOT
    ;
    push @cmds, <<EOT
java -Djava.io.tmpdir=$tmpdir \\
-Xmx4g \\
-jar $picard \\
MarkDuplicates \\
I=$bam_base.bam \\
O=$bam_base.rmdups.bam \\
M=$bam_base.rmdups.metrics \\
CREATE_INDEX=TRUE \\
TMP_DIR=$tmpdir 

rm $bam_base.bam

EOT
    ;
    push @cmds, <<EOT
java -Djava.io.tmpdir=$tmpdir \\
-Xmx4g \\
-jar $gatk \\
-T RealignerTargetCreator \\
-R $human_fasta \\
-known $RealBin/ref/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf.gz \\
-known $RealBin/ref/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \\
-I $bam_base.rmdups.bam \\
-o $bam_base.rmdups.indelrealign.intervals \\
-nt $opts{t}

EOT
    ;
    push @cmds, <<EOT
java -Djava.io.tmpdir=$tmpdir \\
-Xmx4g \\
-jar $gatk \\
-T IndelRealigner \\
-R $human_fasta \\
-known $RealBin/ref/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf.gz \\
-known $RealBin/ref/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \\
-I $bam_base.rmdups.bam \\
-targetIntervals $bam_base.rmdups.indelrealign.intervals \\
-o $bam_base.rmdups.indelrealign.bam

rm $bam_base.rmdups.bam 

EOT
    ;
    push @cmds, <<EOT
java -Djava.io.tmpdir=$tmpdir \\
-Xmx4g \\
-jar $gatk \\
-T BaseRecalibrator \\
-R $human_fasta \\
-knownSites $RealBin/ref/GRCh37/dbSNP147/All_20160601.vcf.gz \\
-knownSites $RealBin/ref/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf.gz \\
-knownSites $RealBin/ref/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \\
-I $bam_base.rmdups.indelrealign.bam \\
-o $bam_base.rmdups.indelrealign.bqsr.grp \\
-nct $opts{t}

java -Djava.io.tmpdir=$tmpdir \\
-Xmx4g \\
-jar $gatk \\
-T PrintReads \\
-R $human_fasta \\
-I $bam_base.rmdups.indelrealign.bam \\
-o $bam_base.rmdups.indelrealign.bqsr.bam \\
-BQSR $bam_base.rmdups.indelrealign.bqsr.grp \\
-nct $opts{t}
    
EOT
    ;
    push @cmds, <<EOT
java -Djava.io.tmpdir=$tmpdir \\
-Xmx4g \\
-jar $gatk \\
-T BaseRecalibrator \\
-R $human_fasta \\
-knownSites $RealBin/ref/GRCh37/dbSNP147/All_20160601.vcf.gz \\
-knownSites $RealBin/ref/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf.gz \\
-knownSites $RealBin/ref/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz \\
-I $bam_base.rmdups.indelrealign.bam \\
-o $bam_base.rmdups.indelrealign.bqsr.post_recal.grp \\
-BQSR $bam_base.rmdups.indelrealign.bqsr.grp \\
-nct $opts{t}

rm $bam_base.rmdups.indelrealign.bam

EOT
    ;
    my @threads = ($opts{t} , $opts{t}, $opts{t}, $opts{t}, 1, $opts{t}, $opts{t});
    my @names = map {"$_.$s"} qw /  align 
                                    sort 
                                    dedup 
                                    realign_target 
                                    indelrealign 
                                    bqsr 
                                    post_bqsr/;
    my @scripts = ();
    for (my $i = 0; $i < @cmds; $i++){
        push @scripts, make_bam_script($names[$i], $cmds[$i], $threads[$i]);
    }
    return @scripts;
}

##################################################
sub make_bam_script{
    my ($name, $cmd, $threads) = @_;
    my $script = "$sub_dir/$name.sh";
    open (my $SCRIPT, ">", $script) or die "Could not create $script: $!\n";
    my $head = <<EOT
#!/bin/bash
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -cwd
#\$ -V
#\$ -l h_rt=$opts{r}:00:00 
EOT
    ; 
    if ($threads > 1){
        my $mem = ceil(12/$threads) . "G";
        $head .= "#\$ -pe sharedmem $threads\n#\$ -l h_vmem=$mem";
    }else{
        $head .= "#\$ -l h_vmem=12G";
    }
    print $SCRIPT <<EOT
$head
. /etc/profile.d/modules.sh 
module load igmm/apps/bwa/0.7.12-r1039
module load igmm/libs/ncurses/6.0
module load igmm/libs/htslib/1.4
module load igmm/apps/samtools/1.4

#bail on first error, unset error, return code of any failed command
set -euo pipefail 

$cmd

EOT
    ;
    close $SCRIPT;
    return $script;
}

##################################################
sub make_dir{
    my $dir = shift;
    if (not -d $dir){
        mkdir ($dir) or die "Could not create directory '$dir': $!\n";
    }
}

##################################################
sub usage{
    my $msg = shift;
    print STDERR "$msg\n" if $msg;
    print <<EOT

USAGE: $0 [options] FQLIST

OPTIONS:
    
    FQLIST - A list of FASTQ files (named in the format sample_name.fastq(.gz)
             to process. Required.

    -o STORE --out_dir STORE
        Directory to place final BAM and GVCF files. 
        Default=/gpfs/igmmfs01/eddie/igmm_datastore_MND-WGS/GRCh37

    -b N_PER_BATCH --batch_size N_PER_BATCH
        Number of samples to process in parallel.

    -n --no_exec
        Perform a dry run. This will create qsub scripts and do a 'mock' qsub,
        but will not submit any jobs. 

    -t THREADS --threads THREADS
        Number of threads to use for multithreaded commands (i.e. bwa and 
        supporting GATK commands). Default = 1.

    -r RUNTIME --runtime RUNTIME
        Number of hours runtime to give each script. Must be a whole number.
        Default = 48.

    -d TMPDIR   --tmp_dir TMPDIR
        Location of directory to use for temporary files. Defaults to the 
        system default TMPDIR.

    -g GATK.jar --gatk GATK.jar
        Location of GATK jar file to use. 
        Default = $RealBin/GATK/GenomeAnalysisTK.jar

    -p picard.jar --gatk picard.jar
        Location of picard jar file to use. 
        Default = $RealBin/picard-tools-2.9.0/picard.jar

    -h --help
        Show this message and exit.

EOT
    ;
    exit 2 if $msg;
    exit;
}


