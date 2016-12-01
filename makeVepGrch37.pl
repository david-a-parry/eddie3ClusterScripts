#!/usr/bin/perl
use strict;
use warnings;
use File::Basename; 
die "Usage: $0 input.vcf [input2.vcf input3.vcf ...]\n" if not @ARGV;

print <<EOT 
#\$ -M david.parry\@ed.ac.uk
#\$ -m a
#\$ -cwd
#\$ -V
#\$ -pe sharedmem 8
#\$ -l h_rt=96:00:00
#\$ -l h_vmem=4G
. /etc/profile.d/modules.sh
module load igmm/libs/htslib/1.3
module load igmm/apps/samtools/1.2

EOT
;

while (my $in = shift){
    my ($f, $d) = fileparse($in); 
    my $vep_out = $d . "vep.$f";
    $vep_out =~ s/\.(b)*gz$//i;
    
    print <<EOT
perl /exports/igmm/eddie/aitman-lab/ensembl-tools-release-84/scripts/variant_effect_predictor/variant_effect_predictor.pl  --dir /exports/igmm/eddie/aitman-lab/ensembl-tools-release-84/scripts/variant_effect_predictor/vep_cache/  --vcf --offline --everything --buffer_size 200 --check_alleles --gencode_basic --assembly GRCh37 --force --fork 8 --plugin LoF   --plugin MaxEntScan,/exports/eddie3_homes_local/dparry/maxentscan/fordownload/  --plugin SpliceConsensus --plugin dbNSFP,/exports/igmm/eddie/aitman-lab/ref/GRCh37/dbNSFP3.2_GRCh37.gz,Polyphen2_HVAR_pred,Polyphen2_HVAR_rankscore,LRT_pred,LRT_converted_rankscore,MutationTaster_pred,MutationTaster_converted_rankscore,MutationAssessor_pred,MutationAssessor_score_rankscore,FATHMM_pred,FATHMM_converted_rankscore,PROVEAN_pred,PROVEAN_converted_rankscore,CADD_phred,CADD_raw_rankscore,DANN_rankscore,DANN_score,MetaSVM_pred,MetaSVM_rankscore,MetaLR_pred,MetaLR_rankscore,Eigen-phred,Eigen-raw_rankscore,Eigen-PC-raw_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_rankscore,GERP++_RS_rankscore,GERP++_RS,phyloP20way_mammalian_rankscore,phyloP20way_mammalian,phastCons20way_mammalian,phastCons20way_mammalian_rankscore,GTEx_V6_gene,GTEx_V6_tissue,clinvar_trait,clinvar_clnsig --plugin dbscSNV,/exports/igmm/eddie/aitman-lab/ref/GRCh37/dbscSNV.txt.gz -i $in -o $vep_out

bgzip -f $vep_out

tabix -f -p  vcf $vep_out.gz
EOT
;
}
