#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use File::Basename;
use Data::Dumper;

my $usage = "Description:
	This script is used to analysis ePAT-SMRT-seq data.

Usage:
	perl $0 [options]

Options:
	Input options:
	--ccs                      [string]        Input ccs fasta file
	--sample                   [string]        Sample name
	--species                  [string]        Species, mm9, mm10, hg18, hg38
	--pass                     [string]        CCS pass file

	Output options:
	--out_dir                  [string]        Output directory, default out
	
	minimap2 options:
	--minimap2_index           [string]        Minimap2 index directory
	--minimap2_thread          [int]           Number of threads used by minimap2

	polyA tail options:
	--min_pal                  [int]           Minimal polyA tail length, default 15
	--min_pass                 [int]           Minimal pass when count non-A residues in polyA tail
	--limit                    [int]           Limited bases from 3UTR or 3-terminal to analysis
	
	APA analysis options:
	--do_apa                                   Do apa analysis or not
	--tapis_thread             [int]           Number of threads used by alignPacBio.py (tapis package)
	--gmap_index_path          [string]        GMAP index path
	--gmap_index_name          [string]        GMAP index name
	--gmap_reference           [string]        Reference genome used by GMAP
	--gtf                      [string]        Genomic annotation file in gtf/gff format
	
	other options:
	--help|-h                                  Print help information

Author:
	Hu Nie, niehu\@genetics.ac.cn

Version:
	Version 1.0.0, 2019-05-06
";

if ( @ARGV == 0 ) {
    print $usage;
    exit 0;
}

# set program
my $root     = abs_path($0);
my $root_dir = dirname($root);
require("$root_dir/scripts/basic.pm20190608");

# input/output
my $ccs_fasta = "";
my $sample    = "";
my $out_dir   = "out";

# call polyA length and non-A
my $min_pass  = 10;
my $min_pal   = 15;
my $species   = "";
my $limit     = 100;
my $pass_file = "";

# minimap2
my $minimap2_index  = "";
my $minimap2_thread = 10;

# APA analysis
my $do_apa          = 0;
my $tapis_thread    = 40;
my $gmap_index_path = "";
my $gmap_index_name = "";
my $gmap_reference  = "";
my $gtf_file        = "";

# others
my $cmd    = "";
my $ok     = 0;
my $fisher = "$root_dir/scripts/fisher.pl";
my $help   = 0;

&GetOptions(

    # input/output
    "ccs=s"     => \$ccs_fasta,
    "sample=s"  => \$sample,
    "out_dir=s" => \$out_dir,

    # minimap2
    "minimap2_index=s"  => \$minimap2_index,
    "minimap2_thread=i" => \$minimap2_thread,

    # call polyA
    "pass=s"     => \$pass_file,
    "min_pal=s"  => \$min_pal,
    "min_pass=i" => \$min_pass,
    "species=s"  => \$species,
    "limit=i"    => \$limit,

    # APA analysis
    "do_apa|apa"        => \$do_apa,
    "tapis_thread=i"    => \$tapis_thread,
    "gmap_index_path=s" => \$gmap_index_path,
    "gmap_index_name=s" => \$gmap_index_name,
    "gmap_reference=s"  => \$gmap_reference,
    "gtf=s"             => \$gtf_file,

    # others
    "help|h" => \$help
) or die $usage;

if ($help) {
    print $usage;
    exit 0;
}

# Start
print &getTime() . " Start analysis\n\n";
my $current_dir = getcwd();
print "Current directory: $current_dir\n";

if ($ccs_fasta) {    # check input ccs
    $ccs_fasta = abs_path($ccs_fasta);
}
else {
    print "ccs option is not specified or empty\n";
    exit();
}

if ($pass_file) {    # check pass
    $pass_file = abs_path($pass_file);
    print $pass_file. "\n";
}
else {
    print "pass option is not specified or empty\n";
    exit();
}

if ($sample) {       # check sample name
    print "Sample: $sample\n";
}
else {
    print "sample option is not specified or empty\n";
}

if ( -e $out_dir ) {    # check outdir
    print &getTime() . " Output directory: " . $out_dir . " exists. Skip creating\n\n";
}
else {
    mkdir($out_dir);
    print &getTime() . " Output directory: " . $out_dir . " does not exists, Creating\n\n";
}

# Check programs
print &getTime() . " Check programs...\n";
my @programs = ( "samtools", "bedtools", "bam2fasta", "seqkit", "minimap2" );
&checkPrograms(@programs);
print("\n");

# Get absolute path
$out_dir         = abs_path($out_dir)         if $out_dir;
$gtf_file        = abs_path($gtf_file)        if $gtf_file;
$gmap_reference  = abs_path($gmap_reference)  if $gmap_reference;
$gmap_index_path = abs_path($gmap_index_path) if $gmap_index_path;

chdir($out_dir);
$current_dir = $out_dir;
print "Current directory: $out_dir\n";

########## Align CCS to reference genome ##########
print &getTime() . " Align CCS to reference genome\n";
system("mkdir -p align");
chdir("$current_dir/align");

# run minimap2
$cmd =
"minimap2 -t $minimap2_thread -ax splice -uf  --secondary=no -C5 --cs=short $minimap2_index $ccs_fasta 2>$sample.align.log | samtools view -bS -F2308 > $sample.ccs.aligned.bam";
exit if !&run($cmd);

# statistics
my $input_reads   = `cat $ccs_fasta | grep \\> | wc -l`;
my $aligned_reads = `samtools view -c $sample.ccs.aligned.bam`;
chomp($input_reads);
chomp($aligned_reads);

# log
print "Input reads: $input_reads\n";
my $aligned_ratio = sprintf( "%.2f", int($aligned_reads) / int($input_reads) );
print "Aligned reads: $aligned_reads ($aligned_ratio)\n";
print "Alignment Done !\n ";
print(" \n ");
chdir($current_dir);

########## Assign CCS to genes ###########
print &getTime() . " Assign CCS to genes\n";
system("mkdir -p ccs2gene");
chdir("$current_dir/ccs2gene");

# Extract exons from genome
my $exon_bed = "";
if ( $species eq "mm10" ) {
    $exon_bed = "$root_dir/data/mm10.exons.bed";
}

# ccs2exon, $4 => CCS name, $NF => Gene name
#$cmd ="bedtools bamtobed -split -i $current_dir/align/$sample.ccs.aligned.bam | bedtools intersect -f 0.5 -a - -b $exon_bed -wa -wb 2>/dev/null | awk '{print \$4, \$NF}' | sort -u >$sample.ccs2exon.txt";
$cmd =
"bedtools bamtobed -split -i $current_dir/align/$sample.ccs.aligned.bam | bedtools intersect -s -f 0.5 -a - -b $exon_bed -wa -wb 2>/dev/null | awk '{print \$4, \$(NF-2)}' | sort -u >$sample.ccs2exon.txt";
exit if !&run($cmd);

# ccs2genes
$cmd = "cat $sample.ccs2exon.txt | awk '{print \$1}'| uniq -c | awk '\$1 == 1 {print \$2}' | $fisher - $sample.ccs2exon.txt > $sample.ccs2gene.txt";
exit if !&run($cmd);

# Number of CCSs assigned
my $ccs_assigned = `cat $sample.ccs2gene.txt | wc -l `;
chomp($ccs_assigned);
print "Number of CCS assigned to genes: $ccs_assigned\n";

# Number of genes covered
my $gene_covered = `awk '{print \$2}' $sample.ccs2gene.txt | sort -u | wc -l`;
chomp($gene_covered);
print "Number of genes covered: $gene_covered\n";

print("\n");
chdir($current_dir);

########## Calculate polyA tail length ##########
print &getTime() . " Calculate polyA tail length\n";
system("mkdir -p polyA_tail_length");
chdir("$current_dir/polyA_tail_length");

#calculate polyA tail length
print "Call polyA tail...\n";
my $trimed_ccs = "ccs.polyA_trimmed.fa"; # add
&call_polyA_tail(
    $sample,
    "$current_dir/align/$sample.ccs.aligned.bam",
    "$current_dir/ccs2gene/$sample.ccs2gene.txt",
    $pass_file, "$sample.polyA_tail_length.txt", $min_pal, $trimed_ccs
);
$trimed_ccs = abs_path($trimed_ccs);
print "Call polyA tail result: $sample.polyA_tail_length.txt\n";

#plot polyA tail length distribution : transcript level
system("$root_dir/scripts/plot_polyA_tail_length.R $sample.polyA_tail_length.txt $sample 1 3 transcript violin");
system("$root_dir/scripts/plot_polyA_tail_length.R $sample.polyA_tail_length.txt $sample 1 3 transcript density");

# plot polyA tail length distribution : gene level
system("$root_dir/scripts/plot_polyA_tail_length.R $sample.polyA_tail_length.txt $sample 1 3 gene violin");
system("$root_dir/scripts/plot_polyA_tail_length.R $sample.polyA_tail_length.txt $sample 1 3 gene density");

print "\n";
chdir($current_dir);

########## Calculate non-A residues in polyA tail ##########
print &getTime() . " Calculate non-A residues in polyA tail\n";
system("mkdir -p polyA_tail_nonA");
chdir("$current_dir/polyA_tail_nonA");
&call_nonA_in_polyA_tail( $sample, "$current_dir/polyA_tail_length/$sample.polyA_tail_length.txt", $min_pass, $min_pal, $limit, 1 );  # 1:3UTR(100nt)=> terminal
&call_nonA_in_polyA_tail( $sample, "$current_dir/polyA_tail_length/$sample.polyA_tail_length.txt", $min_pass, $min_pal, $limit, 2 );  # 2:terminal(100nt)=> 3UTR
&call_nonA_in_polyA_tail( $sample, "$current_dir/polyA_tail_length/$sample.polyA_tail_length.txt", $min_pass, $min_pal, 100, 3 )
  ;    # 3:3UTR => terminal, scaled

# plot, 3UTR => terminal, limit 100 nt
$cmd = "$root_dir/scripts/plot_nonA_distribution_in_polyA_tail.R $sample.3UTR_nonA.txt $sample.3UTR_nonA.pdf $sample  1 $limit";
if ( system($cmd) ) {
    print "Output file: $sample.3UTR_nonA.pdf\n";
}

# plot, terminal => 3UTR, limit 100 nt
$cmd = "$root_dir/scripts/plot_nonA_distribution_in_polyA_tail.R $sample.terminal_nonA.txt $sample.terminal_nonA.pdf $sample 2 $limit";
if ( system($cmd) ) {
    print "Output file: $sample.terminal_nonA.pdf\n";
}

# plot, terminal => 3UTR, scale
$cmd = "$root_dir/scripts/plot_nonA_distribution_in_polyA_tail.R $sample.scaled_nonA.txt $sample.scaled_nonA.pdf $sample 3 100";
if ( system($cmd) ) {
    print "Output file: $sample.scaled_nonA.pdf\n";
}

print "\n";
chdir($current_dir);

########## APA analysis ########
if ($do_apa) {
    if ( !$gmap_index_path ) {
        print "Please specify gmap_index_path option!\n";
        exit;
    }
    if ( !$gmap_index_name ) {
        print "Please specify gmap_index_name option!\n";
        exit;
    }
    if ( !$gmap_reference ) {
        print "Please specify gmap_reference option!\n";
        exit;
    }
    if ( !$gtf_file ) {
        print "Please specify gtf option!\n";
        exit;
    }

    print &getTime() . " Do apa analysis using TAPIS\n";
    system("mkdir -p tapis");
    chdir("$current_dir/tapis");

    @programs = ( "alignPacBio.py", "run_tapis.py" );
    &checkPrograms(@programs);

    $cmd = "alignPacBio.py -p $tapis_thread -o $sample.alignPacBio $gmap_index_path $gmap_index_name $gmap_reference $trimed_ccs &>$sample.alignPacBio.log";
    exit if !&run($cmd);

    $cmd = "run_tapis.py -s 3 -o $sample.tapis_out $gtf_file $sample.alignPacBio/aligned.bam";
    exit if !&run($cmd);

    chdir("$current_dir");
    print &getTime() . " Done\n";
    print "\n";
}
