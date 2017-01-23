#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw/strftime/;
use File::Basename;
use Data::Dumper;
use Getopt::Long;

my @gvcfs = (); 
my $date = strftime( "%d-%m-%y", localtime );
my %opts = 
(
    i => \@gvcfs,
    f => "/exports/igmm/software/pkg/el7/apps/bcbio/share2/genomes/Hsapiens/hg38/seq/hg38.fa",
    t => "$ENV{HOME}/scratch/tmp/",
    g => "$ENV{HOME}/GATK/v3.6/GenomeAnalysisTK.jar", 
    v => "combined_gvcfs-$date",
    n => 100,
    
);
GetOptions(
    \%opts,
    "f|fasta=s",
    "g|gatk",
    "h|help",
    "i|gvcfs=s{,}",
    "n|number_per_batch=i",
    "o|output_dir=s",
    "q|qsub",
    "s|split_by_chr",
    "t|tmp_dir=s",
    "v|vcfname=s",
) or die "Syntax error\n";
usage() if $opts{h};
usage("-o/--output_dir option is required.\n") if not $opts{o};
usage("-i/--gvcfs option is required.\n") if not @gvcfs;
usage("-n/--increment option is must be greater than 0.\n") if $opts{n} < 1;

@gvcfs = map { "\"$_\"" } @gvcfs; 
if (not -d $opts{o}){ 
    mkdir $opts{o} or die "Could not create output directory '$opts{o}: $!\n";
}
my $per_chrom_dir = "$opts{o}/per_chrom";
if ($opts{s} and not -d $per_chrom_dir){ 
    mkdir $per_chrom_dir or die "Could not create per_chrom output directory '$per_chrom_dir: $!\n";
}
    
my $script_dir = "$opts{o}/subscripts";
if (not -d $script_dir){ 
    mkdir $script_dir or die "Could not create script directory '$script_dir: $!\n";
}
my %g_scripts = ();
my %join_scripts = ();
for (my ($i, $j) = (0, 0); $j < @gvcfs; ($i++, $j += $opts{n})){
    my @out_vcf = ();
    my $e = $j + $opts{n} - 1; 
    $e = $e < @gvcfs ? $e : $#gvcfs;
    my $vcf_string = "-V " . join (" -V ", @gvcfs[$j..$e]); 
    my ($scr, $out) =  makeScripts($vcf_string, $i);
    push @out_vcf, @$out;
    push @{$g_scripts{$i}}, @$scr;
    if ($opts{s}){
        $join_scripts{$i} = makeJoinScript(\@out_vcf, $i);
    }
}


if ($opts{q}){
    submit_all_scripts();
}

print STDERR "Done\n";

##################################################################################
sub makeJoinScript{
    my ($vcffiles, $i) = @_;
    my $input = join(' -V ', @$vcffiles);
    my $output = "$opts{o}/$opts{v}.$i.g.vcf.gz";
    my $script = "$opts{o}/subscripts/catGvcfs-$i.sh";
    open (my $CAT, ">", $script) or die "Could not create script $script: $!\n";
    print $CAT <<EOT
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

java -Djava.io.tmpdir=$opts{t} -Xmx12g -cp $opts{g} -R $opts{f} --assumeSorted -out $output  -V $input

EOT
;
    close $CAT;
    return $script;
}

##################################################################################
sub makeScripts{
    my $vcfs = shift;
    my $n = shift;
    my @sc = ();
    my @out = ();
    if ($opts{s}){
        foreach my $chr (readDict()){
            my ($s, $o)  = makeGtScript($vcfs, $n, $chr);
            push @sc, $s;
            push @out, $o;
        }
    }else{
        my ($s, $o) = makeGtScript($vcfs, $n);
        push @sc, $s;
        push @out, $o;
    }
    return (\@sc, \@out);
}

##################################################################################
sub makeGtScript{
    my ($vcf_string, $i, $chr) = @_;
    my $script = "$opts{o}/subscripts/combineGvcfs-$i";
    my $runtime = "96:00:00";
    my $interval = "";
    if ($chr){
        $script .= ".$chr.sh";
        $runtime = "24:00:00";
        $interval = " -L $chr ";
    }else{
        $script .= ".sh";
    }
    $script =~ s/[\:\*]/-/g; #characters not allowed in jobnames
    my $out_vcf = defined ($chr) 
                  ? "$opts{o}/per_chrom/$opts{v}.$chr.$i.g.vcf.gz" 
                  : "$opts{o}/$opts{v}.$i.g.vcf.gz";
    open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
    print $SCRIPT <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -cwd
#\$ -V
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -l h_rt=$runtime 
#\$ -l h_vmem=14G
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load  igmm/apps/bcbio/20160119

java -Djava.io.tmpdir=$opts{t} -Xmx8g -jar $opts{g} -T CombineGVCFs -R $opts{f} -o $out_vcf $vcf_string $interval
EOT
    ;
    close $SCRIPT;
    return ($script, $out_vcf);
}

##################################################################################
sub readDict{
    (my $dict = $opts{f}) =~ s/\.fa(sta){0,1}$/\.dict/;
    if (not -e $dict){
        $dict = "$opts{f}.dict";
        if (not -e $dict){
            die "Could not find .dict file for $opts{f}\n";
        }
    }
    open (my $DICT, $dict) or die "Could not read $dict: $!\n";
    #my %contigs = (); 
    my @chr = ();
    while (my $line = <$DICT>){
       if ($line =~ /^\@SQ\s+.*SN:(\S+)\s+LN:(\d+)\s/){
            #$contigs{$1} = $2;
            push @chr, $1;
        }
    }
    return @chr;
}
 
##################################################################################
sub usage{
    my $msg = shift;
    print STDERR "\n$msg\n" if $msg;
    print STDERR <<EOT

USAGE: $0 -o <output_dir> -i gvcf1 [gvcf2 ... gvcfN]

OPTIONS:

    -o,--output_dir DIR
        Directory to put output VCF files.

    -i,--gvcfs FILES
        One or more input GVCF files to genotype. Required.
    
    -v,--vcfname STRING
        Name for vcf file. Output VCFs will be named <vcfname>.g.vcf.gz etc. 
    
    -n,--number_per_batch INT
        Combine this many GVCFs into each combined ouptut file. Default = 100.
    
    -s,--split_by_chr
        Use this flag to split command into one per chromosome.

    -q,--qsub
        Use this flag to submit scripts after creation.

    -g,--gatk FILE
        Location of GATK jar file. Default = $ENV{HOME}/GATK/v3.6/GenomeAnalysisTK.jar
    
    -t,--tmp_dir DIR
        Directory to use for tmp files. Defalt = $ENV{HOME}/scratch/tmp/
    
    -f,--fasta FILE
        Location of reference genome fasta. Default = /exports/igmm/software/pkg/el7/apps/bcbio/share/bcbio-nextgen/genomes/Hsapiens/hg38/seq/hg38.fa
    
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
sub submit_all_scripts{
    foreach my $k (sort {$a <=> $b} keys %g_scripts){
        my @wait_ids = (); 
        foreach my $sc (@{$g_scripts{$k}}){
            my $cmd = "qsub $sc";
            push @wait_ids, do_qsub($cmd);
        }
        if (exists $join_scripts{$k}){
            my $wait = join(",", @wait_ids);
            do_qsub("qsub -hold_jid $wait $join_scripts{$k}");
        }
    }
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


