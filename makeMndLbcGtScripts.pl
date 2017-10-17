#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use POSIX qw/ceil strftime/;
use File::Basename;

my $rdir = "/exports/igmm/eddie/aitman-lab";
if ($ENV{USER} eq 'clogan2'){
    $rdir = "/exports/igmm/eddie/mopd";
}
my $date = strftime( "%F", localtime );
my @genomic = ();
my @per_chrom = ();
#my @per_chrom = glob '/exports/igmm/eddie/igmm_datastore_MND-WGS/LBC_gvcf/batched_gvcf/*gz';
#my @mnd = glob '/exports/igmm/eddie//igmm_datastore_MND-WGS/X15030_X0008AT_10094AT_X16084_combined/*gz';
my %opts = (p => \@per_chrom, g => \@genomic, v => "wgs-$date");
GetOptions
(
    \%opts,
    "v|vcf_suffix=s",
    "p|per_chrom=s{,}",
    "g|genomic=s{,}",
    "n|no_exec",
    "b|blacklist_samples=s",
    "t|target_sites=s",
    "a|after_cat",
    "w|wait=s",  #all commands wait til this is done
) or die "error in option spec\n";

my $dummy_wait = 0;

my $split = 5000000;
my $fasta = "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/seq/hg38.fa";
my $dict  = "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/seq/hg38.dict";
my $tmp = "$ENV{HOME}/scratch/tmp";

my $omni = "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/variation/1000G_omni2.5.vcf.gz";
my $snps = "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/variation/1000G_phase1.snps.high_confidence.vcf.gz";
my $dbsnp = "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/variation/dbsnp-147.vcf.gz";
my $mills = "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/variation/Mills_and_1000G_gold_standard.indels.vcf.gz";
my $hapmap = "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/variation/hapmap_3.3.vcf.gz",
my $refine = "$rdir/ref/hg38/1000G_phase3_v4_20130502.snvs_only.hg38liftover.sites.vcf.gz";
 
#my $lbc_exclude = "/gpfs/igmmfs01/eddie/igmm_datastore_MND-WGS/joint_genotype_MND_LBC/lbc_excluded_samples.txt";
#my $mnd_exclude = "/gpfs/igmmfs01/eddie/igmm_datastore_MND-WGS/joint_genotype_MND_LBC/mnd_excluded_samples.txt";

die "No GVCFs!\n" if not @genomic and not @per_chrom and not $opts{a};

open (my $DICT, $dict) or die "Could not read $dict: $!\n";
my %contigs = (); 
while (my $line = <$DICT>){
   if ($line =~ /^\@SQ\s+.*SN:(\S+)\s+LN:(\d+)\s/){
        $contigs{$1} = $2;
    }
}
die "No contigs found in $dict!\n" if not %contigs;

mkdir("subscripts") if not -d "subscripts";
mkdir("per_chrom") if not -d "per_chrom";
my %batches = ();
foreach my $g (@per_chrom){
    my ($f, $d) = fileparse($g);
    if ($f =~ /^(batch\d+[a-z])\.comb.*\.(chr\S+)\.g\.vcf\.gz$/){
        $batches{$1}->{$2} = $g;
    }else{
        die "ERROR: Could not parse unexpected GVCF filename: $g\n";
    }
}

foreach my $k (keys %batches){
    foreach my $c ( 1..22, 'X', 'Y' ) {
        my $chr = "chr$c";
        if (not exists $batches{$k}->{$chr}){
            die "ERROR: $chr not found for batch $k!\n";
        }
        if (not exists $contigs{$chr}){
            die "ERROR: $chr not found in $dict!\n";
        }
    }
}

my @wait_ids    = ();
my @out_files   = ();
my @vep_out     = ();
my @recal_out   = ();
my @snp_recal   = ();
my @indel_recal = ();
my @geno_recal = ();
my $out_stub = "var.$opts{v}";
foreach my $c ( 1..22, 'X', 'Y' ) {
    my $chr = "chr$c";
    my @v = map { $batches{$_}->{$chr} } sort keys %batches;
    for (my $i = 0; $i < $contigs{$chr}; $i += $split){
        my $start = $i + 1;
        my $end = ($i + $split) < $contigs{$chr} ? $i + $split : $contigs{$chr}; 
        unless ($opts{a}){
            my $gt_script = makeGtScript($chr, $start, $end, \@v);
            my $w = doQsub("qsub $gt_script");
            my $vep_script = makeVepScript($chr, $start, $end);
            push @wait_ids, doQsub("qsub -hold_jid $w $vep_script");
        }
        my ($snp, $indel, $gcp, $recal_out) = makeApplyRecalScripts($chr, $start, $end); 
        push @snp_recal, $snp;
        push @indel_recal, $indel;
        push @geno_recal, $gcp;
        push @recal_out, $recal_out;
    }
}

my $cat_wait;
unless ($opts{a}){
    $cat_wait    = concatVars(); 
}
my @apply_wait = vqsrAndApply();
my $cat_vqsr_wait = concatVqsr();


=cut
my $recal_genos = "subscripts/recal_genos.sh";
open (my $GSCRIPT, ">$recal_genos") or die "Error writing to $recal_genos: $!\n";
print $GSCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -cwd
#\$ -V
#\$ -l h_rt=192:00:00
#\$ -l h_vmem=18G
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load  igmm/apps/bcbio/20160119

java -Djava.io.tmpdir=$tmp -Xmx12g -jar $ENV{HOME}/GATK/v3.6/GenomeAnalysisTK.jar -R $fasta  -T CalculateGenotypePosteriors --supporting $refine -V $out_stub.ts99pt9.vcf.gz -o $out_stub.ts99pt9.postGCP.vcf.gz -XL chrY

EOT
;
close $GSCRIPT; 
my $r_cmd = "qsub -hold_jid $cat_vqsr_wait $recal_genos";
doQsub($r_cmd);
=cut 

################################################
sub doQsub{
    my $cmd = shift;
    if ($opts{w}){
        $cmd =~ s/qsub/qsub -hold_jid $opts{w}/;
    }
    informUser("EXECUTING: $cmd");
    return ++$dummy_wait if $opts{n};
    my $stdout = `$cmd`; 
    if ($stdout =~ /Your job (\d+) .* has been submitted/){
        return $1;
    }else{
        die "Error parsing qsub stdout for '$cmd'\nOutput was: $stdout";
    }
}

#################################################
sub vqsrAndApply{
    if ($opts{a}){
        return submitApplyRecalScripts();
    }else{
        my $snp_recal   = makeVqsr("snp");
        my $indel_recal = makeVqsr("indel");
        my $snp_recal_wait    = doQsub("qsub -hold_jid $cat_wait $snp_recal"); 
        my $indel_recal_wait  = doQsub("qsub -hold_jid $cat_wait $indel_recal");
        return submitApplyRecalScripts($snp_recal_wait, $indel_recal_wait);
    }
}
#################################################
sub submitApplyRecalScripts{
    my ($snp_wait, $indel_wait) = @_;
    my @apply_wait = ();
    my $snp_hold = '';
    if (defined $snp_wait){
        $snp_hold = "-hold_jid $snp_wait";
    }
    for (my $i = 0; $i < @snp_recal; $i++){ 
        my $w = doQsub("qsub $snp_hold $snp_recal[$i]");
        $w .= ",$indel_wait" if defined $indel_wait;
        my $gcp_wait = doQsub("qsub -hold_jid $w $indel_recal[$i]");
        push @apply_wait, doQsub("qsub -hold_jid $gcp_wait $geno_recal[$i]");
    }
    return @apply_wait;
}

#################################################
sub makeConcatScript{
    my $script = shift;
    my $output = shift;
    my $f = shift;
    open (my $CAT, ">$script") or die "Error writing to $script: $!\n";
    #my $join_string = join(" -V ", @$f); 
    #vep commands which parse a vcf with 0 variants will not produce an output - we have to check if file exist
    open (my $LIST, ">$script.filelist") or die "Error writing to $script.filelist: $!\n";
    print $LIST join("\n", @$f);
    print $CAT <<EOT
#!/bin/bash
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -cwd
#\$ -V
#\$ -l h_rt=48:00:00
#\$ -l h_vmem=18G
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load  igmm/apps/bcbio/20160119

#vep commands which parse a vcf with 0 variants will not produce an output - we have to check if files exist
VSTRING=''
for F in `cat $script.filelist`; do
    if [ -f \$F ]; then
        VSTRING="\$VSTRING -V \$F"
    else
        echo "WARN: \$F does not exist"
    fi
done

java -Djava.io.tmpdir=$tmp -Xmx12g -cp $ENV{HOME}/GATK/v3.6/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants  -R $fasta --assumeSorted -out $output  \$VSTRING

EOT
;
    close $CAT;
    return $script;
}

#################################################
sub concatVqsr{
    my $cat_script = "subscripts/catVqsr.sh";
    makeConcatScript($cat_script, "$out_stub.ts99pt9.postGcp.vcf.gz", \@recal_out);
    return doQsub("qsub -hold_jid " .join(",", @apply_wait ) . " $cat_script"); 
}

#################################################
sub concatVars{
    my $cat_script = "subscripts/catVar.sh";
    makeConcatScript($cat_script, "$out_stub.raw.vcf.gz", \@out_files);
    return doQsub("qsub -hold_jid " .join(",", @wait_ids) . " $cat_script"); 
}

#################################################
sub makeVqsr{
#VQSR needs to be done on whole dataset
    my $type = shift;
    my $mode = uc($type);
    my $script = "subscripts/recal_$type.sh";
    open (my $SCRIPT, ">$script") or die "Error writing to $script: $!\n";
    my $resources = '';
    if ($type eq 'snp'){
        $resources = "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap ".
                     "-resource:omni,known=false,training=true,truth=true,prior=12.0 $omni ".
                     "-resource:1000G,known=false,training=true,truth=false,prior=10.0 $snps ". 
                     "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp";
    }else{
        $resources = "-resource:mills,known=false,training=true,truth=true,prior=12.0 $mills ".
                     "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp";
    }
print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -cwd
#\$ -V
#\$ -l h_rt=48:00:00
#\$ -l h_vmem=48G
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load  igmm/apps/bcbio/20160119
module load R

java -Djava.io.tmpdir=$tmp -Xmx12g -jar $ENV{HOME}/GATK/v3.6/GenomeAnalysisTK.jar -R $fasta -T VariantRecalibrator $resources -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode $mode -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0  -input $out_stub.raw.vcf.gz -recalFile $out_stub.recalibrate_$mode.recal -tranchesFile  $out_stub.recalibrate_$mode.tranches -rscriptFile $out_stub.recalibrate_$mode.plots.R

EOT
;

    close $SCRIPT; 
    return $script;
}

#################################################
sub makeApplyRecalScripts{
    my ($chr, $start, $end) = @_;
    my @scripts;
    my ($s_script, undef) = makeApplyRecal($chr, $start, $end, "snp");
    my ($i_script, $i_output) = makeApplyRecal($chr, $start, $end, "indel");
    my ($r_script, $output) = makeGcp($chr, $start, $end, $i_output);
    return ($s_script, $i_script, $r_script, $output);
}

#################################################
sub makeGcp{
    my ($chr, $start, $end, $input) = @_;
    my $recal_script = "subscripts/recalGenos_$chr-$start-$end.sh";
    open (my $RECAL, ">", $recal_script) or die "Could not open $recal_script for writing: $!\n";
    (my  $out = $input) =~ s/\.vcf\.gz$//;
    $out .= ".postGcp.vcf.gz";#do this in two steps to make sure suffix gets added 
                          # even if problem with substitute command

    my $cmd = "java -Djava.io.tmpdir=$tmp -Xmx3g -jar $ENV{HOME}/GATK/v3.6/GenomeAnalysisTK.jar -R $fasta  -T CalculateGenotypePosteriors --supporting $refine -V $input -o $out -L $chr:$start-$end";
    if ($chr eq 'chrY'){#CalculateGenotypePosteriors chokes on haploid GTs
        $cmd = "echo 'GCP will not be carried out on chrY'";
        $out = $input;
    }
    print $RECAL <<EOT
#!/bin/bash
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -e $recal_script.stderr
#\$ -o $recal_script.stdout
#\$ -cwd
#\$ -V
#\$ -l h_rt=8:00:00
#\$ -l h_vmem=8G
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load  igmm/apps/bcbio/20160119

if [ -f $input ]; then
    $cmd
else
    echo "WARN: $input does not exist"
fi

EOT
;

close $RECAL; 
    return ($recal_script, $out);
}


#################################################
sub makeApplyRecal{
#can apply recal in smaller batches in parallel
#but indel recal must take place on output from snp recal
    my ($chr, $start, $end, $type) = @_;
    my $recal_script = "subscripts/$type"."ApplyRecal_$chr-$start-$end.sh";
    open (my $APPLY, ">", $recal_script) or die "Could not open $recal_script for writing: $!\n";
    my $suffix = "$opts{v}";
    if ($opts{b}){
        $suffix = "$opts{v}.valid_samples";
    }
    my $in;
    my $out;
    if ($type eq 'snp'){
        $in = "per_chrom/vep.var.$chr-$start-$end.$suffix.raw.vcf.gz";
        $out = "per_chrom/vep.var.$chr-$start-$end.$suffix.snpTs99pt9.vcf.gz";
    }else{
        $in  = "per_chrom/vep.var.$chr-$start-$end.$suffix.snpTs99pt9.vcf.gz";
        $out = "per_chrom/vep.var.$chr-$start-$end.$suffix.ts99pt9.vcf.gz";
    }
    my $mode = uc($type); 
    print $APPLY <<EOT
#!/bin/bash
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -e $recal_script.stderr
#\$ -o $recal_script.stdout
#\$ -cwd
#\$ -V
#\$ -l h_rt=4:00:00
#\$ -l h_vmem=8G
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load  igmm/apps/bcbio/20160119

if [ -f $in ]; then
    java -Djava.io.tmpdir=$tmp -Xmx3g -jar $ENV{HOME}/GATK/v3.6/GenomeAnalysisTK.jar -R $fasta -T ApplyRecalibration -input $in -recalFile $out_stub.recalibrate_$mode.recal -tranchesFile  $out_stub.recalibrate_$mode.tranches -mode $mode --ts_filter_level 99.9 -o $out -L $chr:$start-$end
else
    echo "WARN: $in does not exist"
fi

EOT
;

close $APPLY; 
    return ($recal_script, $out);
}

#################################################
sub getTargetSitesExcluseSites{
    my ($chr, $start, $end) = @_;
    my @exclude = ();
    if ($start > 1){
        push @exclude, "$chr:1-" . ($start - 1) ;
    }
    if ($end < $contigs{$chr}){
        push @exclude, "$chr:" . ($end + 1) . "-$contigs{$chr}";
    }
    foreach my $c ( 1..22, 'X', 'Y' ) {
        my $x_chr = "chr$c";
        if ($x_chr ne $chr){
            push @exclude, $x_chr;
        }
    }
    return "-XL " . join(" -XL ", @exclude);
}

#################################################
sub makeGtScript{
    my ($chr, $start, $end, $vcfs) = @_;
    my $gt_script = "subscripts/gtGvcf$chr-$start-$end.sh";
    open (my $GT, ">", $gt_script) or die "Could not open $gt_script for writing: $!\n";
    my $exclude_region = '';
    if ($chr eq 'chrY' and $start <= 56887903 and $end >= 56887903){
    # discrepancy between LBC and MND GVCFs at this coordinate
        $exclude_region = ' -XL  chrY:56887903-56887903 ';
    }
    my $interval = "-L $chr:$start-$end " ;
    if ($opts{t}){
        $interval = getTargetSitesExcluseSites($chr, $start, $end);
        $interval .= " -allSites -L $opts{t} ";
    }
        
    my $pre_out = "per_chrom/var.$chr-$start-$end.$opts{v}.raw.vcf.gz";
    my $out = "per_chrom/var.$chr-$start-$end.$opts{v}.valid_samples.raw.vcf.gz";
    #push @out_files, $out;
    print $GT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -cwd
#\$ -V
#\$ -e $gt_script.stderr
#\$ -o $gt_script.stdout
#\$ -l h_rt=18:00:00
#\$ -l h_vmem=32G
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load  igmm/apps/bcbio/20160119
EOT
    ;
    print $GT "java -Djava.io.tmpdir=$tmp -Xmx24g -jar ".
              "$ENV{HOME}/GATK/v3.6/GenomeAnalysisTK.jar ".
              "-T GenotypeGVCFs -stand_call_conf 30 -stand_emit_conf 10 ".
              $exclude_region .
              $interval .
              "-R $fasta ".
              "-D $dbsnp ".
              "-o $pre_out ".
              "-V " . join(" -V ", @genomic, @$vcfs) . "\n";
    if ($opts{b}){
        print $GT "java -Djava.io.tmpdir=$tmp -Xmx24g -jar ".
                  "$ENV{HOME}/GATK/v3.6/GenomeAnalysisTK.jar ".
                  "-T SelectVariants " .
                  "-R $fasta ".
                  "-xl_sf $opts{b} " .
                  "-V $pre_out ".
                  "-o $out ".
                  "-env\n";
        close $GT;
    }
    return $gt_script;
}

#################################################
sub makeVepScript{
    my ($chr, $start, $end) = @_;
    my $vep_script = "subscripts/vepVcf$chr-$start-$end.sh";
    open (my $VEP, ">", $vep_script) or die "Could not open $vep_script for writing: $!\n";
    my $suffix = $opts{v};
    if ($opts{b}){
        $suffix .= ".valid_samples";
    }
    my $in = "per_chrom/var.$chr-$start-$end.$suffix.raw.vcf.gz";
    my $vep_out = "per_chrom/vep.var.$chr-$start-$end.$suffix.raw.vcf";
    my $out = "$vep_out.gz";
    push @out_files, $out;
    my $dbnsfp = "$rdir/ref/hg38/dbNSFP_hg38.gz";
    my $maxent = "$rdir/maxentscan/fordownload/";
    my $plugins = "$rdir/vep_plugins";
    my $dbscsnv = "$rdir/ref/hg38/dbscSNV_hg38.txt.gz";
    print $VEP <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -cwd
#\$ -V
#\$ -e $vep_script.stderr
#\$ -o $vep_script.stdout
#\$ -l h_rt=48:00:00
#\$ -l h_vmem=32G
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load igmm/libs/htslib/1.3
module load igmm/apps/samtools/1.2
module load igmm/apps/perl/5.24.0

perl $rdir/ensembl-tools-release-84/scripts/variant_effect_predictor/variant_effect_predictor.pl \\
--vcf --offline --everything --buffer_size 200 --check_alleles \\
--gencode_basic --assembly GRCh38 --force --dir_plugins $plugins \\
--dir $rdir/ensembl-tools-release-84/scripts/variant_effect_predictor/vep_cache/ \\
--plugin LoF   --plugin MaxEntScan,$maxent --plugin SpliceConsensus \\
--plugin dbNSFP,$dbnsfp,Polyphen2_HVAR_pred,Polyphen2_HVAR_rankscore,\\
LRT_pred,LRT_converted_rankscore,MutationTaster_pred,\\
MutationTaster_converted_rankscore,MutationAssessor_pred,\\
MutationAssessor_score_rankscore,FATHMM_pred,FATHMM_converted_rankscore,\\
PROVEAN_pred,PROVEAN_converted_rankscore,CADD_phred,CADD_raw_rankscore,\\
DANN_rankscore,DANN_score,MetaSVM_pred,MetaSVM_rankscore,MetaLR_pred,\\
MetaLR_rankscore,Eigen-phred,Eigen-raw_rankscore,Eigen-PC-raw_rankscore,\\
fathmm-MKL_coding_pred,fathmm-MKL_coding_rankscore,GERP++_RS_rankscore,\\
GERP++_RS,phyloP20way_mammalian_rankscore,phyloP20way_mammalian,\\
phastCons20way_mammalian,phastCons20way_mammalian_rankscore,GTEx_V6_gene,\\
GTEx_V6_tissue,clinvar_trait,clinvar_clnsig \\
--plugin dbscSNV,$dbscsnv \\
-i $in \\
-o $vep_out

bgzip -f $vep_out

tabix -f -p  vcf $vep_out.gz
EOT
    ;
    close $VEP;
    return $vep_script;
}
#################################################
sub informUser{
    my $msg = shift;
    my $time = strftime( "%H:%M:%S", localtime );
    print STDERR "[$time] $msg\n";
}
