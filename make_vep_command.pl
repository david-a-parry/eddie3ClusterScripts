#!/usr/bin/perl
use strict;
use warnings;
use File::Basename; 
die "Usage: $0 input.vcf [input2.vcf input3.vcf ...]\n" if not @ARGV;

print <<EOT 
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m a
#\$ -m b
#\$ -m e
#\$ -V
#\$ -cwd
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/tabixpp/0.2.5
module load apps/gcc/perl/5.22.0
module load apps/gcc/samtools/1.2
EOT
;

while (my $in = shift){
    my ($f, $d) = fileparse($in); 
    my $vep_out = $d . "vep.$f";
    print <<EOT
perl   ~/ensembl-tools-release-81/scripts/variant_effect_predictor/variant_effect_predictor.pl --vcf  --offline --everything --plugin Carol --plugin LoF,human_ancestor_fa:~/GRCh37/human_ancestor.fa.gz,conservation_file=/export/users/dparry/GRCh37 --plugin Condel,/export/users/dparry/.vep/Plugins/config/Condel/config,b --plugin SpliceConsensus  --plugin GeneSplicer,/export/users/dparry/bin/genesplicer,/export/users/dparry/bin/splice_training_set/human/,cache_size=1000,context=200,tmpdir=/mnt/lustre2/dparry/tmp/ --plugin dbscSNV,/export/users/dparry/dbscSNV/dbscSNV.txt.gz --force --fork 8 --plugin MaxEntScan,/export/users/dparry/bin/maxentscan/fordownload/ -i $in  -o $vep_out
EOT
;
}
