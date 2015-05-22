#!/usr/bin/perl -w

# generate make file for doing LHMEI_v4 discovery phase

use strict;
use Getopt::Long;

my $BamList;
my $OutDir;
my $MakeBase;
my $Nodes = "-j10,11,12,13"; # default 4-node minicluster
my $SingleChr = 0; # default whole genome
my $LOCAL = 0;
my $SINGLEOUT = 0;
my $NO_CTRL_VCF = 0;
my $NODEPTH = 0;
my $MAPPER="";
my $NO_DP_INFO = 0;

my $Discovery = "/net/wonderland/home/saichen/LHMEI_v4/bin/LHMEI-Discovery";

# print help info
if( !(@ARGV) ) {
	print "\nUsage: [] is default.\n";
	print "-i  list of bam file. name bam (outdir) (depth)\n";
	print "-o  output dir\n";
	print "-m  base name of make file\n";
	print "-a  mapper. If not specified, LHMEI will silently use bwa mem.\n";
	print "-j  specify nodes [-j10,11,12,13]\n";
	print "-c  specify chr to work on. 0 means whole genome. [0]\n";
	print "--local  run on local [off]\n";
	print "--singleOut  output dir is specified in col3 of bam list. Output to dir/LHMEI [off]\n";
	print "--nodpinfo  no dp info in bam list. [off]\n";
	print "--noCtrlVcf  do not print refLH*.vcf for ctrl sliced bam\n";
	print "--nodepth  do not use depth filter. [off[\n";
	print "\n";
	exit 0;
}

GetOptions(
	'i=s' => \$BamList, # bam list. format: col1(sample name) col2( bam ). Tab delimited
	'o=s' => \$OutDir, # output dir. For each sample, the script will create a sample/ dir for it.
	"a=s" => \$MAPPER,
	'm=s' => \$MakeBase, # base name of mk
	'j=s' => \$Nodes, # optional. -j10,11,12 to specify nodes
	'c=s' => \$SingleChr, # optional. if not specified, will do whole genome
	"local" => \$LOCAL,
	"singleOut" => \$SINGLEOUT,
	"noCtrlVcf" => \$NO_CTRL_VCF,
	"nodepth" => \$NODEPTH,
	"nodpinfo" => \$NO_DP_INFO,
);


my $MOS = "mosbatch -E/tmp $Nodes sh -c";

my $last_char = substr( $OutDir, length($OutDir) - 1, 1 );
if ( $last_char eq '/' ) {
	$OutDir = substr( $OutDir, 0, length($OutDir) - 1 );
}
my $MakeFile = "$OutDir/$MakeBase";

print "Generating $MakeFile...\n";
if ( $SingleChr > 0 ) {
	print "\tOnly for chr $SingleChr...\n";
}
else {
	print "\tDoing whole genome...\n";
}

open( IN, "<", $BamList ) or die $!;
my @Names; # sample names
my @Bams; # bam file full name
my @Dps;
my $Cols = 2;
if ( $NO_DP_INFO == 0 ) {
	$Cols = $Cols + 1;
}
if ( $SINGLEOUT == 1 ) {
	$Cols = $Cols + 1;
}
if ( $NODEPTH == 1 ) {
	$Cols = $Cols - 1;
}
my @Dirs;
while(<IN>) {
	chomp;
	my $line = $_;
	my @s = split '\t', $line;
	if ( @s < $Cols ) {
		die "At line $.: $line do not have $Cols fields!";
	}
	push( @Names, $s[0] );
	push( @Bams, $s[1] );
	if ( $NO_DP_INFO == 0) {
		push( @Dps, $s[-1] );
	}
	if ( $SINGLEOUT == 1 ) {
		push( @Dirs, $s[2] );
	}
}
close IN;

# write to make file
# headers
open( OUT, ">", $MakeFile ) or die $!;
print OUT ".DELETE_ON_ERROR:\n\n";
print OUT "all:";
for( my $i = 0; $i < @Names; $i++ ) {
	my $ok;
	if ( $SINGLEOUT == 0 ) {
		$ok = "$OutDir/$Names[$i]/Discovery.OK";
	}
	else {
		$ok = "$Dirs[$i]/LHMEI/Discovery.OK";
	}
	print OUT " $ok";
}
print OUT "\n\n";

# samples
for( my $idx = 0; $idx < @Names; $idx++ ) {
	my $name = $Names[$idx];
	my $dir;
	if ( $SINGLEOUT == 0 ) { # no outdir info in bam list
		$dir = "$OutDir/$name";
	}
	else {
		$dir = "$Dirs[$idx]/LHMEI";
	}
	print OUT "$dir/Discovery.OK:\n";
	print OUT "\tmkdir -p $dir\n";
	my $cmd = "$Discovery -Sample $name -Bam $Bams[$idx] -WorkDir $dir --verbose";
	if ( $SingleChr > 0 ) {
		$cmd = "$cmd -Chr $SingleChr";
	}
	if ( $NO_CTRL_VCF > 0 ) {
		$cmd = "$cmd --noCtrlVcf";
	}
	if ( $NODEPTH > 0 ) {
		$cmd = "$cmd --disableDPfilter";
	}
	elsif ( @Dps > 0 ) {
		$cmd = "$cmd -Depth $Dps[$idx]";
	}
	if ( length($MAPPER) >= 1 ) {
		$cmd = "$cmd -Mapper $MAPPER";
	}
	$cmd = $cmd." > $dir/run.$SingleChr.log 2> $dir/run.$SingleChr.err";
	
# see if need to add node option
	if ( $LOCAL == 0 ) { # run on nodes
		my $cmd2 = "\t$MOS '$cmd'\n";
		$cmd = $cmd2;
	}
	else {
		my $cmd2 = "\t$cmd\n";
		$cmd = $cmd2;
	}
	print OUT $cmd;
	print OUT "\ttouch $dir/Discovery.OK\n";
	print OUT "\trm -f $dir/ctrl_tmps/*fastq $dir/ctrl_tmps/20-nsort.bam $dir/ctrl_tmps/20-remap-sort.bam $dir/ctrl_tmps/align-pe.sam* $dir/ctrl_tmps/ctrl-disc.sam\n";
	print OUT "\trm -f $dir/preprocess/20.sam $dir/preprocess/disc.sam\n";
	print OUT "\trm -f $dir/preprocess/disc-nsort.bam\n\n";
}

close OUT;

print "Finish creating $MakeFile!\n";


