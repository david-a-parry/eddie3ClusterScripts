 #!/usr/env perl 
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
my $align;
my $sort_and_merge;
my $indel_realign;
my $recalibrate;
my @interval_list = ();
my $outdir;
my $help;
my %config;
GetOptions(
    \%config,
    "align",
    "sort_and_merge",
    "markduplicates",
    "indel_realign",
    "recalibrate",
    "list=s{,}" => \@interval_list,
    "outdir=s" => \$outdir,
    "qsub",
    "print_scripts",
    "help",
) or die "Syntax error!\n";
usage() if $config{help};

usage("A file list is required!") if @ARGV != 1;

if (not keys %config){
    for my $o (qw / align sort_and_merge markduplicates indel_realign recalibrate / ){
        $config{$o}++;
    }
}
my $interval_string = "";
foreach my $int (@interval_list){
    $interval_string .= "-L \"$int\" ";
}

my $file_list = shift;
my %fastqs = ();
open (my $FL, $file_list) or die "Can't open $file_list: $!\n";
while (my $file = <$FL>){
    chomp($file);
    my ($sample, $snum, $date, $lane, $read); 
    my ($f, $d) = fileparse($file);
    if ($f !~ /\.fastq(\.gz)$/){
        print STDERR "$d$f does not look like a FASTQ file - skipping.\n";
        next;
    }
=cut
    if ($d =~ /Sample_(\w+)([-_]\w+)*/){
        $sample = $1;
        $sample .= $2 if $2;
    }elsif ($f =~ /^([\w\-_]+)_L00/){
        $sample = $1;
        $sample .= $2 if $2;
=cut
    if ($f =~ /^([\w\-_]+)_(S\d+|[ACTG]{10})_L00/){
        $sample = $1;
        $snum = $2;
    }else{
        print STDERR "WARNING: Couldn't process sample name for $d$f - skipping.\n";
        next;
    }
    if ($d =~ /(\d{6})_\w{5,6}_\w{4}_\w{10}/){
        $date = $1;
    }elsif($d =~ /\/(\d{6})_/){
        $date = $1;
    }else{
        print STDERR "WARNING: Couldn't process date for $d$f - skipping.\n";
        next;
    }
    if ($f =~ /\S+_L(\d{1,3})_R(\d)_\d+/){
        $lane = $1;
        $read = $2;
        #add to hash
    }else{
        print STDERR "WARNING: Couldn't process lane and read details for $d$f - skipping.\n";
        next; 
    }
    if ($sample =~ /undetermined/i){
        print STDERR "INFO: Skipping undetermined file: $d$f\n";
        next;
    }elsif($sample =~ /^water$/i or $sample =~ /^dH2O$/i){
        $sample = "dH2O-$snum-$date";
        print STDERR "INFO: Converting water sample name to '$sample' for $d$f\n";
    }elsif ($sample =~ /^blank$/i){
        $sample = "Blank-$snum-$date";
        print STDERR "INFO: Converting blank sample name to '$sample' for $d$f\n";
    }
    $sample = check_sample_name($sample, "$date-$snum", $lane, $read, $d, $f); 
    $fastqs{$sample}->{"$date-$snum"}->{$lane}->{$read} = "$d$f";
}
my @bams = ();
my @scripts = ();
foreach my $s (sort keys %fastqs){
    foreach my $d (sort keys %{$fastqs{$s}}){
        foreach my $l (sort keys %{$fastqs{$s}->{$d}}){
            my @reads = ();
            foreach my $r (sort  keys %{$fastqs{$s}->{$d}->{$l}}){
                push @reads, $fastqs{$s}->{$d}->{$l}->{$r};
            }
            if (scalar @reads != 2){
                warn "WARNING: $s - $d - $l has " . scalar(@reads) . " reads\n";
            }
            my $dir;
            if ($outdir){
                $dir = $outdir;
            }else{
                $dir = dirname($reads[0]); 
            }
            my $bam = "$dir/$s-$d-$l.bam";
            push @{$fastqs{$s}->{bams}}, $bam;
            push @bams, $bam;
            if ($config{align}){
                if (not -d "alignments"){
                    mkdir("alignments") or die "Failed to make alignments directory: $!\n";
                }
                my $script = "alignments/align-$s-$d-$l.sh";
                $script = check_exists_and_rename($script);
                if (-e $bam){
                    warn "WARNING: $bam BAM file already exists!\n";
                }
                open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
                my $rg = "'\@RG\\tID:$s-$d-$l\\tLB:$s-$d\\tSM:$s\\tPL:ILLUMINA'";
                my $readstring = "\"$reads[0]\"";
                if (@reads > 1){
                    $readstring .= " \"$reads[1]\"";
                }
                print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -cwd
#\$ -e $script.stderr
#\$ -o $script.stdout
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/bwa/0.7.12
module load apps/gcc/samtools/1.2
bwa mem /export/users/dparry/lustre2/dparry/GRCh37/hs37d5.fa -t 8 -M -R $rg $readstring | samtools view -Sb - > "$bam"

EOT
;
                close $SCRIPT;
                push @scripts, $script; 
            }
        }
    }
}

check_duplicates();

#SORT/MERGE
foreach my $s (sort keys %fastqs){
    my $sort_cmd = "java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/picard/picard-tools-1.131/picard.jar ";
    my $out_bam;
    next if not $fastqs{$s}->{bams}; 
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
        my $script = "sortmerge/sortmerge-$s.sh";
        $script = check_exists_and_rename($script);
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
module load apps/gcc/jre/1.7.0_60
$sort_cmd
EOT
;
        close $SCRIPT;
        push @scripts, $script; 
    }

#DEDUP
    my $in_bam = $out_bam;
    my $metrics_file = "$out_bam" ."_rmupds.metrics";
    if ($config{markduplicates}){
        $out_bam =~ s/\.bam$//;
        $out_bam .= "_rmdups.bam";
        if (not -d "rmdups"){
            mkdir("rmdups") or die "Failed to make rmdups directory: $!\n";
        }
        my $script = "rmdups/markdups-$s.sh";
        $script = check_exists_and_rename($script);
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
    if ($config{indel_realign}){
        $out_bam = $stub . "_indelrealign.bam";
        if (not -d "indelrealign"){
            mkdir("indelrealign") or die "Failed to make indelrealign directory: $!\n";
        }
        my $script = "indelrealign/indelrealign-$s.sh";
        $script = check_exists_and_rename($script);
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
module load apps/gcc/jre/1.7.0_60
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -known /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf -known /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf -I "$in_bam" -o "$intervals" $interval_string --filter_mismatching_base_and_quals
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -known /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf -known /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf -I "$in_bam" -targetIntervals "$intervals" -o "$out_bam" $interval_string --filter_mismatching_base_and_quals
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
    if ($config{recalibrate}){
        $out_bam = $stub . "_bqsr.bam";
        if (not -d "recal"){
            mkdir("recal") or die "Failed to make recal directory: $!\n";
        }
        my $script = "recal/recalibration-$s.sh";
        $script = check_exists_and_rename($script);
        open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
        print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -V
#\$ -cwd
#\$ -e $script.stderr
#\$ -o $script.stdout
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/jre/1.7.0_60
module load apps/gcc/R/3.1.0
#get recalibration model
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -knownSites /mnt/lustre2/dparry/GRCh37/dbSNP142/00-All.vcf.gz -knownSites /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf -knownSites /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf -I "$in_bam" -o "$grp" $interval_string --filter_mismatching_base_and_quals
#apply recalibration model
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T PrintReads -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -I "$in_bam" -o "$out_bam" -BQSR "$grp" $interval_string --filter_mismatching_base_and_quals
#check recalibration model
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -knownSites /mnt/lustre2/dparry/GRCh37/dbSNP142/00-All.vcf.gz -knownSites /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/1000G_phase1.indels.b37.vcf -knownSites /mnt/lustre2/dparry/GRCh37/GATK_resource_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf -I "$in_bam" -BQSR "$grp" -o "$postrecal" $interval_string --filter_mismatching_base_and_quals
#produce plots
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T AnalyzeCovariates -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa -before "$grp" -after "$postrecal" -plots "$plots" $interval_string --filter_mismatching_base_and_quals
EOT
;
        close $SCRIPT;
        push @scripts, $script; 
    }    

}
print STDERR "Created following scripts:\n" . join("\n", @scripts) ."\n" if $config{print_scripts};
if ($config{qsub}){
    submit_scripts();
}

sub submit_scripts{
    my @align = grep {/^alignments/} @scripts;
    my @sortmerge = grep {/^sortmerge/} @scripts;
    my @rmdups = grep {/^rmdups/} @scripts;
    my @indelrealign = grep {/^indelrealign/} @scripts;
    my @recal = grep {/^recal/} @scripts;
    foreach my $sc (@align, @sortmerge, @rmdups, @indelrealign, @recal){
        system("qsub $sc");
    }
}

##################################################################################
sub check_exists_and_rename{
    my $f = shift;
    if (-e $f){
        warn "WARNING: $f already exists!\n";
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
sub check_sample_name{
	my ($sample, $date, $lane, $read, $d, $f) = @_;
    if (not exists $fastqs{$sample}){
        return $sample;
    }
    my $rename = rename_sample($sample); 
    if ($fastqs{$sample}->{$date}){
        if ($fastqs{$sample}->{$date}->{$lane}){
            if ($fastqs{$sample}->{$date}->{$lane}->{$read}){
                warn "WARN: File $d/$f is a duplicate lane ($lane), read ($read) and sample ($sample). Trying sample name of $rename...\n";
                return check_sample_name($rename, $date, $lane, $read, $d, $f);
            }else{
                return $sample;#other read pair for sample - keep sample name
            }
        }
    }#sample exists and isn't a pair of an existing read - rename to something unique
    #warn "WARN: File $d/$f is a duplicate sample ($sample). Trying sample name of $rename...\n";
    #return check_sample_name($rename, $date, $lane, $read, $d, $f);
    return $sample;
}
##################################################################################
sub rename_sample{
    my $sample = shift;
	if ($sample =~ /\.(\d+)$/){
		my $dup = $1 + 1;
		$sample =~ s/\.(\d+)$//;
		$sample .= ".$dup";
	}else{
		$sample = "$sample.1";
	}
	return $sample;
}
##################################################################################
sub check_duplicates{
    my %bam_check = ();
    foreach my $bam (@bams){
        if (exists $bam_check{$bam}){
            die "Duplicate output bam: $bam!\n";
        }
    }
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
    -l    --list            [optional interval list to use with GATK commands]
    -h    --help            [show this help message and exit]

EOT
;
    exit 1 if $msg;
    exit;
}
