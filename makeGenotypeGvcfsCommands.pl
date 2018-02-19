#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw/strftime/;
use File::Basename;
use Data::Dumper;
use Getopt::Long;

my @gvcfs = (); 
my $date = strftime( "%d-%m-%y", localtime );
my $rdir = "/exports/igmm/eddie/aitman-lab";
if ($ENV{USER} eq 'clogan2'){
    $rdir = "/exports/igmm/eddie/mopd";
}

my %opts = 
(
    i => \@gvcfs,
    f => "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/seq/hg38.fa",
    d => "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/seq/hg38.dict",
    t => "$ENV{HOME}/scratch/tmp/",
    g => "$ENV{HOME}/GATK/v3.8/GenomeAnalysisTK.jar", 
    omni => "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/variation/1000G_omni2.5.vcf.gz",
    s => "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/variation/1000G_phase1.snps.high_confidence.vcf.gz",
    dbsnp => "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/variation/dbsnp-147.vcf.gz",
    m => "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/variation/Mills_and_1000G_gold_standard.indels.vcf.gz",
    hapmap => "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/variation/hapmap_3.3.vcf.gz",
    r => "$rdir/ref/hg38/1000G_phase3_v4_20130502.snvs_only.hg38liftover.sites.vcf.gz",
    v => "variants-$date",
    e => 4,
    c => 30,
);
GetOptions(
    \%opts,
    "c|call_conf=i",
    "dbsnp=s",
    "d|dict=s",
    "e|emit_conf=i",
    "f|fasta=s",
    "g|gatk=s",
    "h|help",
    "hapmap=s",
    "i|gvcfs=s{,}",
    "m|mills=s",
    "n|no_gcp",
    "omni=s",
    "o|output_dir=s",
    "p|ped=s",
    "q|qsub",
    "r|refinement_snps=s",
    "s|phase1_snps=s",
    "t|tmp_dir=s",
    "v|vcfname=s",
    "x|main_chromosomes_only",
) or die "Syntax error\n";
usage() if $opts{h};
usage("-d/--dict option is required.\n") if not $opts{d};
usage("-o/--output_dir option is required.\n") if not $opts{o};
usage("-i/--gvcfs option is required.\n") if not @gvcfs;

open (my $DICT, $opts{d}) or die "Can't open $opts{d} for reading: $!\n";
my @contigs = (); 
while (my $line = <$DICT>){
   if ($line =~ /^\@SQ\s+.*SN:(\S+)/){
        my $c = $1;
        if ($opts{x} and $c !~ /^(chr)?[0-9]?[0-9XYM]$/){
            next;
        }
        push @contigs, $c;
    }
}
die "No contigs found in $opts{d}!\n" if not @contigs;
foreach my $d ($opts{o}, "$opts{o}/subscripts", "$opts{o}/per_chrom"){
    if (not -d $d){
        mkdir($d) or die "could not create output directory '$d': $!\n";
    }
}
@gvcfs = map { "\"$_\"" } @gvcfs; 
my $vcf_string = "-V " . join (" -V ", @gvcfs); 
my @g_scripts = (); 
my @c_vcfs = ();
foreach my $chr (@contigs){
    next if $chr =~ /_decoy$/; #skip decoy chromosomes
    my $script = "$opts{o}/subscripts/genotypeGvcfs-$opts{v}-$chr.sh";
    $script =~ s/[\:\*]/-/g; #characters not allowed in jobnames
    my $chrom_vcf = "$opts{o}/per_chrom/var.$chr.$opts{v}.raw.vcf.gz";
    $chrom_vcf =~ s/[\:\*]/-/g;
    open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
    print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -cwd
#\$ -V
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -l h_rt=32:00:00
#\$ -l h_vmem=16G
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load  igmm/apps/bcbio/20160119

java -Djava.io.tmpdir=$opts{t} -Xmx8g -jar $opts{g} -T GenotypeGVCFs -R $opts{f} -D $opts{dbsnp} -stand_call_conf $opts{c} -stand_emit_conf $opts{e} -L $chr -o $chrom_vcf $vcf_string
EOT
;
    close $SCRIPT;
    push @g_scripts, $script;
    push @c_vcfs, $chrom_vcf;
}

my $recal_script = "$opts{o}/subscripts/join_and_recal-$opts{v}.sh";
open (my $SCRIPT, ">$recal_script") or die "Can't open $recal_script for writing: $!\n";
my $join_string = " -V "  . join(" -V ", @c_vcfs); 
my $pedstring = ''; 
if ($opts{p}){
    $pedstring = "-ped $opts{p}";
}

print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -cwd
#\$ -V
#\$ -e $recal_script.stderr
#\$ -o $recal_script.stdout
#\$ -l h_rt=24:00:00
#\$ -l h_vmem=16G
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load  igmm/apps/bcbio/20160119

java -Djava.io.tmpdir=$opts{t} -Xmx8G -cp $opts{g} org.broadinstitute.gatk.tools.CatVariants  -R $opts{f} --assumeSorted -out $opts{o}/var.$opts{v}.raw.vcf.gz $join_string 

java -Djava.io.tmpdir=$opts{t} -Xmx8G -jar $opts{g} -R $opts{f} -T VariantRecalibrator -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $opts{hapmap} -resource:omni,known=false,training=true,truth=true,prior=12.0 $opts{omni} -resource:1000G,known=false,training=true,truth=false,prior=10.0 $opts{s} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $opts{dbsnp} -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0  -input $opts{o}/var.$opts{v}.raw.vcf.gz -recalFile $opts{o}/var.$opts{v}.recalibrate_SNP.recal -tranchesFile  $opts{o}/var.$opts{v}.recalibrate_SNP.tranches -rscriptFile $opts{o}/var.$opts{v}.recalibrate_SNP.plots.R

java -Djava.io.tmpdir=$opts{t} -Xmx8G -jar $opts{g} -R $opts{f} -T ApplyRecalibration -input $opts{o}/var.$opts{v}.raw.vcf.gz -recalFile $opts{o}/var.$opts{v}.recalibrate_SNP.recal -tranchesFile  $opts{o}/var.$opts{v}.recalibrate_SNP.tranches -mode SNP --ts_filter_level 99.9 -o $opts{o}/var.$opts{v}.snpTs99pt9.vcf.gz

java -Djava.io.tmpdir=$opts{t} -Xmx8G -jar $opts{g} -R $opts{f} -T VariantRecalibrator -resource:mills,known=false,training=true,truth=true,prior=12.0 $opts{m} -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $opts{dbsnp} -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum  -mode INDEL -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0  -input $opts{o}/var.$opts{v}.snpTs99pt9.vcf.gz -recalFile $opts{o}/var.$opts{v}.recalibrate_INDEL.recal -tranchesFile  $opts{o}/var.$opts{v}.recalibrate_INDEL.tranches -rscriptFile $opts{o}/var.$opts{v}.recalibrate_INDEL.plots.R
 
java -Djava.io.tmpdir=$opts{t} -Xmx8G -jar $opts{g} -R $opts{f} -T ApplyRecalibration -recalFile $opts{o}/var.$opts{v}.recalibrate_INDEL.recal -tranchesFile  $opts{o}/var.$opts{v}.recalibrate_INDEL.tranches -mode INDEL --ts_filter_level 99.9 -input $opts{o}/var.$opts{v}.snpTs99pt9.vcf.gz -o $opts{o}/var.$opts{v}.ts99pt9.vcf.gz

EOT
;

my $last_output = "$opts{o}/var.$opts{v}.ts99pt9.vcf.gz";

unless($opts{n}){
    print $SCRIPT 
"java -Djava.io.tmpdir=$opts{t} -Xmx4G -jar $opts{g} -R $opts{f}  -T CalculateGenotypePosteriors --supporting $opts{r} $pedstring -V $opts{o}/var.$opts{v}.ts99pt9.vcf.gz -o $opts{o}/var.$opts{v}.ts99pt9.postGCP.vcf.gz -XL chrM -XL chrY -XL chrY_KI270740v1_random\n\n";
    $last_output = "$opts{o}/var.$opts{v}.ts99pt9.postGCP.vcf.gz";
}

if ($pedstring){
    my $final = $last_output;
    $final =~ s/\.vcf\.gz$/.denovoAnnot.vcf.gz/;
    print $SCRIPT 
"java -Djava.io.tmpdir=$opts{t} -Xmx4G -jar $opts{g} -R $opts{f} -T VariantAnnotator -A PossibleDeNovo  $pedstring -V $last_output -o $final\n";
}

close $SCRIPT;


if ($opts{q}){
    submit_scripts();
}

print STDERR "Done\n";


##################################################################################
sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;
    print STDERR <<EOT

USAGE: $0 -d <dict> -o <output_dir> -i gvcf1 [gvcf2 ... gvcfN]

OPTIONS:

    -d,--dict
        Dict file for reference genome for getting chromosome names. Default = /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/seq/hg38.dict

    -o,--output_dir
        Directory to put output VCF files.

    -i,--gvcfs
        One or more input GVCF files to genotype. Required.
    
    -v,--vcfname
        Name for vcf file. Output VCFs will be named var.<vcfname>.raw.vcf.gz etc. 

    -q,--qsub
        Use this flag to submit scripts after creation.

    -g,--gatk
        Location of GATK jar file. Default = $ENV{HOME}/GATK/v3.6/GenomeAnalysisTK.jar
    
    -p,--ped
        Optional ped file for all samples.

    -t,--tmp_dir
        Directory to use for tmp files. Defalt = $ENV{HOME}/scratch/tmp/
    
    -f,--fasta
        Location of reference genome fasta. Default = /exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/seq/hg38.fa
    
    -m,--mills
        Mills and 1000G gold standard indels

    -r,--refinement_snps
        High confidence SNPs for genotype refinement workflow.

    -s,--phase1_snps
        1000G phase1 high confidence SNPs.

    --dbsnp
        dbSNP file

    --omni
        Omni SNPs. 
    
    --hapmap
        HapMap SNPs. 
    
    -c,--call_conf
        Call confidence for PASS variants (i.e. value to pass to -stand_call_conf option of GenotypeGVCFs).

    -e,--emit_conf
        Call confidence for outputting variants (i.e. value to pass to -stand_emit_conf option of GenotypeGVCFs).

    -n,--no_gcp
        Use this flag to skip the GATK CalculateGenotypePosteriors step.

    -x,--main_chromosomes_only
        Only do main autosomes and sex chromosomes - ignore haplotype contigs,
        unplaced contigs and decoy contigs.

    -h,--help
        Show this message and exit.

Author: 

    David A. Parry
    david.parry\@igmm.ed.ac.uk     

EOT
;
    exit 1 if $msg;
    exit;
}

##################################################################################
sub submit_scripts{
    my @wait_ids = (); 
    foreach my $sc (@g_scripts){
        my $cmd = "qsub $sc";
        print STDERR "Executing: $cmd\n";
        my $output = `$cmd`;
        check_exit($?);
        if ($output =~ /Your job (\d+) .* has been submitted/){
            push @wait_ids, $1;
        }else{
            die "Error parsing qsub output for $cmd!\nOutput was: $output";
        }
    }
    my $hold = join(",", @wait_ids);
    my $cmd = "qsub -hold_jid $hold $recal_script";
    print STDERR "Executing: $cmd\n";
    my $output = `$cmd`;
    check_exit($?);
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


