#!/usr/bin/perl
use strict;
use warnings;
use File::Basename; 
use Getopt::Long;

my %opt = 
(
    r => 16,
    m => '4G',
    f => 1,
); 
GetOptions
(
    \%opt,
    'r|runtime=i', #runtime in hours
    'f|forks=i',
    'm|mem=s',
) or die "Syntax error\n";

die "Usage: $0 input.vcf [input2.vcf input3.vcf ...] [-m <job memory> [-f <threads>] [-r <runtime (hours)>]\n" if not @ARGV;

my $rdir = "/exports/igmm/eddie/aitman-lab";
if ($ENV{USER} eq 'clogan2'){
    $rdir = "/exports/igmm/eddie/mopd";
}

my $dbnsfp = "$rdir/ref/GRCh37/dbNSFP3.2_GRCh37.gz";
my $dbscsnv = "$rdir/ref/GRCh37/dbscSNV.txt.gz";
my $maxent = "$rdir/maxentscan/fordownload/";
my $fasta = "$rdir/ref/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
my $plugins = "$rdir/vep_plugins";
my $sharedmem = '';
if ($opt{f} > 1){
    $sharedmem = "#\$ -pe sharedmem $opt{f}";
}
print <<EOT 
#\$ -m a
#\$ -cwd
$sharedmem
#\$ -V
#\$ -l h_rt=$opt{r}:00:00
#\$ -l h_vmem=$opt{m}
. /etc/profile.d/modules.sh
module load igmm/apps/tabix/0.2.5
module load igmm/apps/perl/5.24.0
module load igmm/apps/vep/86

EOT
;

while (my $in = shift){
    my ($f, $d) = fileparse($in); 
    my $vep_out = $d . "vep.$f";
    $vep_out =~ s/\.(b)*gz$//i;
    
    my $vep_cmd = "vep --fasta $fasta --dir_plugins $plugins  --vcf --offline --everything --buffer_size 200 --check_alleles --gencode_basic --assembly GRCh37 --force --plugin LoF   --plugin MaxEntScan,$maxent --plugin SpliceConsensus --plugin dbNSFP,$dbnsfp,Polyphen2_HVAR_pred,Polyphen2_HVAR_rankscore,LRT_pred,LRT_converted_rankscore,MutationTaster_pred,MutationTaster_converted_rankscore,MutationAssessor_pred,MutationAssessor_score_rankscore,FATHMM_pred,FATHMM_converted_rankscore,PROVEAN_pred,PROVEAN_converted_rankscore,CADD_phred,CADD_raw_rankscore,DANN_rankscore,DANN_score,MetaSVM_pred,MetaSVM_rankscore,MetaLR_pred,MetaLR_rankscore,Eigen-phred,Eigen-raw_rankscore,Eigen-PC-raw_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_rankscore,GERP++_RS_rankscore,GERP++_RS,phyloP20way_mammalian_rankscore,phyloP20way_mammalian,phastCons20way_mammalian,phastCons20way_mammalian_rankscore,GTEx_V6_gene,GTEx_V6_tissue,clinvar_trait,clinvar_clnsig --plugin dbscSNV,$dbscsnv -i $in -o $vep_out";

    $vep_cmd .= " --fork $opt{f}" if $opt{f} > 1;

    print <<EOT

$vep_cmd

bgzip -f $vep_out

tabix -f -p  vcf $vep_out.gz
EOT
;
}
