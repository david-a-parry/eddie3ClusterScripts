#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use POSIX qw(ceil); 

my @bam;
my $dict;
my $out;
my $ped;
my $split = 10000000;
my @lists;
my $num_intervals = 1;
my $help;
my $header = <<EOT
#\$ -M david.parry\@igmm.ed.ac.uk
#\$ -m a
#\$ -m b
#\$ -m e
#\$ -cwd
# Configure modules
. /etc/profile.d/modules.sh
# Load modules
module load compilers/gcc/4.7.3
EOT
;
GetOptions(
        "bams=s{,}"         => \@bam,
        "dict=s"            => \$dict,
        "out=s"             => \$out,
        "split=i"           => \$split,
        "list=s{,}"         => \@lists,
        "num_intervals=i"   => \$num_intervals,
        "help"              => \$help,
) or die "Syntax error.\n";

if ($help or not @bam or not $dict){
    die "Usage: $0 -d <dict> -o <output prefix> -b sample1.bam [sample2.bam sample3.bam ...] [-p pedfile.ped] [-s <interval size>] [-l <interval_list(s)> -n <split into this many commands>]\n" ;
}
open (my $DICT, $dict) or die "$!\n";
my %contigs = (); 
while (my $line = <$DICT>){
   if ($line =~ /^\@SQ\s+.*SN:(\S+)\s+LN:(\d+)\s/){
        $contigs{$1} = $2;
    }
}
die "No contigs found in $dict!\n" if not %contigs;
my ( $outfn, $outdir) = fileparse($out);
if (not $outdir){
    $outdir = ".";
}else{
    if (not -d $outdir){
        mkdir $outdir or die "Could not create output directory $outdir: $!\n";
    }
}
@bam = map { "\"$_\"" } @bam; 
my $bamstring = join (" ", @bam); 
if (@lists){#use intervals file split into $num_intervals
    processIntervals();
}else{#whole genome split into chunks of $split size
    processGenome();
}

###############################################################
sub processIntervals{
    my @int_split = ();
    
    if ($num_intervals > 1){
        my @intervals = ();
        foreach my $list (@lists){
            open (my $INT, $list) or die "Cannot open $list for reading: $!\n";
                while (my $line = <$INT>){
                next if $line =~ /^[@#]/;#headers
                chomp $line; 
                my @split = split(/[\:\-\t]/, $line); 
                push @intervals, \@split;
            }
        }
        @intervals = mergeIntervals(\@intervals); 
        my $in_per_batch = ceil(@intervals/$num_intervals);
        my $n = 0;
        for (my $i = 0; $i < @intervals; $i += $in_per_batch){
            $n++;
            my $tmp_interval_file = "tmp_interval_$n.bed";
            open (my $TMP_INT, ">$tmp_interval_file") or die "Can't open '$tmp_interval_file' for writing: $!\n";
            my $last = $in_per_batch + $i - 1;
            $last = $last < @intervals ? $last : $#intervals;
            for (my $j = $i; $j <= $last; $j++){
                die "Not enough fields for line " . join("\t", @{$intervals[$j]} ) . "\n" 
                  if  @{$intervals[$j]} < 3;
                print $TMP_INT join("\t", @{$intervals[$j]}[0..2] ) . "\n";
            }
            close $TMP_INT;
            push @int_split, $tmp_interval_file;
        }
    }else{
        push @int_split, @lists;
    }
    for (my $i = 0; $i < @int_split; $i++){
        my $n = $i + 1;
        my $script = "freebayes_interval_split-$n.sh";
        open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
        print $SCRIPT <<EOT
$header
freebayes -f ~/GRCh37/hs37d5.fa --genotype-qualities -q 17 -m 10 -t  $int_split[$i] $bamstring > $outdir/intervals-$n.$outfn 
EOT
;
    close $SCRIPT;
        print "qsub $script\n";
    }
}

###############################################################
sub mergeIntervals{
    my ($intervals) = @_; 
    my @merged = ();
    @$intervals = sort { 
        $a->[0] cmp $b->[0] || 
        $a->[1] <=> $b->[1] || 
        $a->[2] <=> $b->[2] 
    } @$intervals;
    my $prev_interval = $intervals->[0]; 
    for (my $i = 1; $i < @$intervals; $i++){
        if ($prev_interval->[0] eq $intervals->[$i]->[0] and
            $intervals->[$i]->[1] <= $prev_interval->[2]){
            $prev_interval->[2] = $intervals->[$i]->[2] >  $prev_interval->[2] ? $intervals->[$i]->[2] : $prev_interval->[2] ;
        }else{
            push @merged, $prev_interval;
            $prev_interval = $intervals->[$i];
        }
    }
    push @merged, $prev_interval;
    return @merged;
}
###############################################################
sub processGenome{
    foreach my $chr (keys %contigs){
        for (my $i = 0; $i < $contigs{$chr}; $i += $split){
            my $start = $i;
            my $end = ($i + $split) < $contigs{$chr} ? $i + $split + 1 : $contigs{$chr} + 1; 
            my $script = "freebayes$chr-$start-$end.sh";
            open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
            print $SCRIPT <<EOT
$header
freebayes -f ~/GRCh37/hs37d5.fa --genotype-qualities -q 17 -m 10 -r $chr:$start-$end $bamstring > $outdir/chr$chr-$start-$end.$outfn
EOT
;
        close $SCRIPT;
            print "qsub $script\n";
        }
    }
}
#~/freebayes/bin/freebayes -f ~/GRCh37/hs37d5.fa --genotype-qualities -q 20 -m 20 --min-coverage 50 -C 10 -F 0.2 -t ../../../taad/genes/taad_both_genes_b37.bed taad_UK_Pools_August_18th_2015-1.bam > test2.vcf

