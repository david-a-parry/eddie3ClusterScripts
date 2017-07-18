 #!/usr/env perl 
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
my $align;
my $sort_and_merge;
my $indel_realign;
my $recalibrate;
my $help;
my %config;
GetOptions(
    \%config,
    "align",
    "sort_and_merge",
    "markduplicates",
    "indel_realign",
    "recalibrate",
    "help",
) or die "Syntax error!\n";
usage() if $config{help};

usage("A file list is required!") if @ARGV != 1;

if (not keys %config){
    for my $o (qw / align sort_and_merge markduplicates indel_realign recalibrate / ){
        $config{$o}++;
    }
}
my $file_list = shift;
my %fastqs = ();
open (my $FL, $file_list) or die "Can't open $file_list: $!\n";
while (my $file = <$FL>){
    chomp($file);
    my ($sample, $date, $lane, $read); 
    my ($f, $d) = fileparse($file);
    if ($f !~ /\.fastq(\.gz)*$/){
        print STDERR "$d$f does not look like a FASTQ file - skipping.\n";
        next;
    }
    if ($d =~ /Sample_(\w+)([-_]\w+)*/){
        $sample = $1;
        $sample .= $2 if $2;
    }else{
        print STDERR "WARNING: Couldn't process sample name for $d$f - skipping.\n";
        next;
    }
    if ($d =~ /(\d{6})_\w{5,6}_\w{4}_\w{10}/){
        $date = $1;
    }else{
        print STDERR "WARNING: Couldn't process date for $d$f - skipping.\n";
        next;
    }
    if ($f =~ /\S+_[CTGA]{6}_L(\d{1,3})_R(\d)_\d+/){
        $lane = $1;
        $read = $2;
        #add to hash
    }else{
        print STDERR "WARNING: Couldn't process lane and read details for $d$f - skipping.\n";
        next; 
    }
    if (exists $fastqs{$sample}->{$date}->{$lane}->{$read}){
        die "ERROR: File $d/$f is a duplicate date ($date), lane ($lane), read ($read) and sample ($sample).\n"; 
    }
    $fastqs{$sample}->{$date}->{$lane}->{$read} = "$d$f";
}
my @scripts = ();
foreach my $s (sort keys %fastqs){
    foreach my $d (sort keys %{$fastqs{$s}}){
        foreach my $l (sort keys %{$fastqs{$s}->{$d}}){
            my @reads = ();
            foreach my $r (sort  keys %{$fastqs{$s}->{$d}->{$l}}){
                push @reads, $fastqs{$s}->{$d}->{$l}->{$r};
            }
            if (scalar @reads != 2){
                die "$s - $d - $l has " . scalar(@reads) . " reads. Need 2! Exiting\n";
            }
            my $dir = dirname($reads[0]); 
            my $bam = "$dir/$s-$d-$l-exome.bam";
            push @{$fastqs{$s}->{bams}}, $bam;
            if ($config{align}){
                if (not -d "alignments"){
                    mkdir("alignments") or die "Failed to make alignments directory: $!\n";
                }
                my $script = "alignments/$s-$d-$l-align.sh";
                $script = check_exists_and_rename($script);
                if (-e $bam){
                    print STDERR "WARNING: $bam BAM file already exists!\n";
                }
                open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
                my $rg = "'\@RG\\tID:$s-$d-$l-exome\\tLB:$s-$d\\tSM:$s\\tPL:ILLUMINA'";
                print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m a
#\$ -m b
#\$ -m e
#\$ -cwd
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/bwa/0.7.12
module load apps/gcc/samtools/1.2
bwa mem /export/users/dparry/lustre2/dparry/GRCh37/hs37d5.fa -t 8 -M -R $rg "$reads[0]" "$reads[1]" | samtools view -Sb - > "$bam"

EOT
;
                close $SCRIPT;
                push @scripts, $script; 
            }
        }
    }
}

#SORT/MERGE
foreach my $s (sort keys %fastqs){
    my $sort_cmd = "java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/picard/picard-tools-1.131/picard.jar ";
    my $out_bam;
    if (@{$fastqs{$s}->{bams}} > 1){
        #if we have more than one bam for a sample merge and put in the first
        #fastq directory
        my $dir = dirname($fastqs{$s}->{bams}->[0]);
        $out_bam = "$dir/$s-exome-merged.bam";
        $sort_cmd .= "MergeSamFiles MSD=TRUE SO=coordinate TMP_DIR=/mnt/lustre2/dparry/tmp/  CREATE_INDEX=TRUE I=\"" 
                . join("\" I=\"", @{$fastqs{$s}->{bams}})
                . "\" O=\"$out_bam\""
            ;
    }else{
        ($out_bam = $fastqs{$s}->{bams}->[0]) =~ s/\.bam$//;
        $out_bam .= "_sorted.bam";
        $sort_cmd .= "SortSam SO=coordinate TMP_DIR=/mnt/lustre2/dparry/tmp/ CREATE_INDEX=TRUE I=\"$fastqs{$s}->{bams}->[0]\" O=\"$out_bam\"";
    }
    if ($config{sort_and_merge}){
        if (not -d "sortmerge"){
            mkdir("sortmerge") or die "Failed to make sortmerge directory: $!\n";
        }
        my $script = "sortmerge/$s-sortmerge.sh";
        $script = check_exists_and_rename($script);
        open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
        print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m a
#\$ -m b
#\$ -m e
#\$ -cwd
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/jre/1.7.0_60
$sort_cmd
EOT
;
        close $SCRIPT;
        push @scripts, $script; 
    }

#DEDUP
    my $in_bam = $out_bam;
    $out_bam =~ s/\.bam$//;
    my $metrics_file = "$out_bam" ."_rmupds.metrics";
    $out_bam .= "_rmdups.bam";
    if ($config{markduplicates}){
        if (not -d "rmdups"){
            mkdir("rmdups") or die "Failed to make rmdups directory: $!\n";
        }
        my $script = "rmdups/$s-markdups.sh";
        $script = check_exists_and_rename($script);
        open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
        print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m a
#\$ -m b
#\$ -m e
#\$ -cwd
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/jre/1.7.0_60
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/picard/picard-tools-1.131/picard.jar MarkDuplicates I="$in_bam" O="$out_bam" M="$metrics_file" CREATE_INDEX=TRUE TMP_DIR=/mnt/lustre2/dparry/tmp/ 
EOT
;
        close $SCRIPT;
        push @scripts, $script; 
    } 


#INDELREALIGN
    $in_bam = $out_bam;
    (my $stub = $in_bam) =~ s/\.bam$//;
    my $intervals = $stub . "_indelrealign.intervals";
    $out_bam = $stub . "_indelrealign.bam";
    if ($config{indel_realign}){
        if (not -d "indelrealign"){
            mkdir("indelrealign") or die "Failed to make indelrealign directory: $!\n";
        }
        my $script = "indelrealign/$s-indelrealign.sh";
        $script = check_exists_and_rename($script);
        open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
        print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m a
#\$ -m b
#\$ -m e
#\$ -cwd
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/jre/1.7.0_60
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -known /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf -known /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf -I $in_bam -o "$intervals"
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -known /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf -known /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf -targetIntervals "$intervals" -o "$out_bam"
EOT
;
        close $SCRIPT;
        push @scripts, $script; 
    }    
                
#BQSR           
    $in_bam = $out_bam;
    ($stub = $in_bam) =~ s/\.bam$//;
    my $grp = $stub . "_bqsr.grp";
    my $postrecal = $stub . "_bqsr_postrecal.grp";
    my $plots = $stub . "_bqsr_postrecal.plots.pdf";
    $out_bam = $stub . "_bqsr.bam";
    if ($config{recalibrate}){
        if (not -d "recal"){
            mkdir("recal") or die "Failed to make recal directory: $!\n";
        }
        my $script = "recal/$s-recalibration.sh";
        $script = check_exists_and_rename($script);
        open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
        print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m a
#\$ -m b
#\$ -m e
#\$ -V
#\$ -cwd
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/jre/1.7.0_60
module load apps/gcc/R/3.1.0
#get recalibration model
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -knownSites /mnt/lustre2/dparry/GRCh37/dbSNP142/00-All.vcf.gz -knownSites /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf -knownSites /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf -I "$in_bam" -o "$grp"
#apply recalibration model
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T PrintReads -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -I "$in_bam" -o "$out_bam" -BQSR "$grp"
#check recalibration model
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -knownSites /mnt/lustre2/dparry/GRCh37/dbSNP142/00-All.vcf.gz -knownSites /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf -knownSites /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf -I "$in_bam" -BQSR "$grp" -o "$postrecal"
#produce plots
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T AnalyzeCovariates -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -before "$grp" -after "$postrecal" -plots "$plots"
EOT
;
        close $SCRIPT;
        push @scripts, $script; 
    }    

}
print STDERR "Created following scripts:\n" . join("\n", @scripts) ."\n";

##################################################################################
sub check_exists_and_rename{
    my $f = shift;
    if (-e $f){
        (my $newname = $f) =~ s/(\.\w+)$//;
        my $ext = $1;
        $newname = $f if not $newname;
        my $add = 1;
        if ($newname =~ /_(\d+)$/){
            $add = $1 + 1;
            $newname =~ s/_(\d+)$//;
        }
        $newname .= "_$add";
        $newname .= "$ext" if $ext;
        return check_exists_and_rename($newname);
    }
    return $f;
}

##################################################################################
sub usage{
    my $msg = shift; 
    print "ERROR: $msg\n" if $msg;
    print <<EOT
    
    USAGE: perl $0 fastq_list.txt [options]
    
    From a list of FASTQs this script will create a series of qsub scripts to perform alignment and post-processing.

    Individual steps can be selected using the options below, but note that the preceding step (in the order below) must have already been run and the relevant input bam file must exist.

    OPTIONS: 

    -a    --align           [create alignment scripts only]
    -s    --sort_and_merge  [create Picard's SortSam or MergeSam scripts (bams from same samples will be merged) only]
    -i    --indel_realign   [create GATK indel realignment scripts only]
    -m    --markduplicates  [create Picard's MarkDuplicates scripts only]
    -r    --recalibrate     [create GATK BQSR scripts only]
    -h    --help            [show this help message and exit]

EOT
;
    exit 1 if $msg;
    exit;
}
