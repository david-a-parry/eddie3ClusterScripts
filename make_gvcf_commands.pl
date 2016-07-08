use strict;
use warnings;
use File::Basename;

die "usage: $0 bam1 [bam2 bam3 ... ] \n" if not @ARGV;

foreach my $bam (@ARGV){
    chomp($bam);
    my ($f, $d) = fileparse($bam); 
    my $script = "genotype_$f.sh";
    open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
    (my $vcf = $f) =~ s/\.bam$/.g.vcf.gz/; 
    $vcf = "/mnt/lustre2/dparry/fetal_exomes/vcfs/gvcfs/$vcf";
    die "$vcf already exists!\n" if -e $vcf ;
    print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m abe
#\$ -cwd
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load apps/gcc/jre/1.7.0_60
java -Djava.io.tmpdir=/mnt/lustre2/dparry/tmp/ -Xmx4g -jar /export/users/dparry/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /mnt/lustre2/dparry/GRCh37/hs37d5.fa  --emitRefConfidence GVCF  -variant_index_type LINEAR -variant_index_parameter 128000 -stand_call_conf 30 -stand_emit_conf 4 -I "$bam" -o "$vcf"
EOT
;
    close $SCRIPT;
    print "qsub $script\n";
}
