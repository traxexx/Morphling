#!/usr/bin/perl -w

# filter LHMEI results based on:
#  QUAL ( pVariant )
# Output:
#	#pass %pass #true %power %1/1-true #novel %1/1-novel

use strict;
use Getopt::Long;

my $InVcfPrefix;
my $BedPrefix = "/net/wonderland/home/saichen/LHMEI_v4/refs/ref-MEI"; # separated .Alu, .L1, .SVA as 0,1,2
my $MinQual;
my $SingleEnd = 0;

GetOptions(
	'i=s' => \$InVcfPrefix,
	'r=s' => \$BedPrefix,
	'q=s' => \$MinQual,
	's=s' => \$SingleEnd
);

if ( length($InVcfPrefix) == 0 ) {
	print "Usage:  with [] are optional\n";
	print "  -i   ***/refLH \n";
	print "	 -q   min-quality for filter\n";
	print "  [-r]   reference bed prefix. [ /net/wonderland/home/saichen/LHMEI_v4/refs/ref-MEI* ]\n";
	print "  [-s}	keep single-end results. [0]\n";
	exit 0;
}

my $out_dir = "/tmp/";
my $TmpLineCount = $out_dir."line-count.txt";

my $AllBed = "$BedPrefix.bed";
# loop through by MEI type
for( my $rec = 0; $rec <= 3; $rec++ ) {
	my $BedName;
	if ( $rec == 3 ) {
		$BedName = "$BedPrefix.bed";}
	else {
		$BedName = "$BedPrefix.$rec";}
	my $bed_count_cmd = "wc -l $BedName > $TmpLineCount";
	system( $bed_count_cmd );
	my $InVcf;
	if ( $rec == 3 ) {
		$InVcf = $out_dir."all.report";
		my $all_cmd = " cat $InVcfPrefix.0.report $InVcfPrefix.1.report $InVcfPrefix.2.report | sort -k2,2n > $InVcf";
		system( $all_cmd );
	}
	else {
		$InVcf = "$InVcfPrefix.$rec.report";}
	my $vcf_count_cmd = "wc -l $InVcf >> $TmpLineCount";
	system( $vcf_count_cmd );

# count line number to # pass, %pass	
	my $FilterVcf = $out_dir."filter.vcf";
	my $filter_cmd;
	if ( $SingleEnd == 1 ) {
		$filter_cmd = "cat $InVcf | awk -v qual=$MinQual \'\$6 >= qual\' | awk 'NR>1' > $FilterVcf";
	}
	else {
		$filter_cmd = "cat $InVcf | awk -v qual=$MinQual \'\$6 >= qual\' | awk 'NR>1' | grep -v \"BOTH_END=0\"> $FilterVcf";
	}
	system( $filter_cmd );
	my $pass_count_cmd = "wc -l $FilterVcf >> $TmpLineCount";
	system( $pass_count_cmd);

# count line number & genotype	
	my $OverlapVcf = $out_dir."overlap.vcf";
	my $intersect_cmd = "cat $FilterVcf | awk '{print \$1,\$2-500,\$2+500, \$10}' | tr ' ' '\\t' | bedtools intersect -a - -b $BedName -u > $OverlapVcf;";
	system( $intersect_cmd );
	my $overlap_count_cmd = "wc -l $OverlapVcf >> $TmpLineCount";
	system( $overlap_count_cmd );
	my $hom_count_cmd = "cat $OverlapVcf | grep \"1\/1\" | wc -l >> $TmpLineCount";
	system( $hom_count_cmd );

# same with Overlap vcf	
	my $NovelVcf = $out_dir."novel.vcf";
	$intersect_cmd = "cat $FilterVcf | awk '{print \$1,\$2-600,\$2+600, \$10}' | tr ' ' '\t' | bedtools intersect -a - -b $AllBed -v > $NovelVcf";
	system( $intersect_cmd );
	my $novel_count_cmd = "wc -l $NovelVcf >> $TmpLineCount";
	system( $novel_count_cmd );
	my $novel_hom_cmd = "cat $NovelVcf | grep \"1\/1\" | wc -l >> $TmpLineCount";
	system( $novel_hom_cmd );

# count numbers
	open( IN, "<", $TmpLineCount ) or die $!;
	my @records;
	while( <IN> ) {
		chomp; my $line = $_;
		if ( $line eq '0' ) {
			push( @records, $line );
		}
		else {
			my @s = split ' ', $line,2;
			push(@records, $s[0]);
		}
	}
	my $bedCount = $records[0];
	my $totalCount = $records[1] - 1;
	my $passCount = $records[2];
	my $passPercent;
	if ( $passCount == 0 ) {
		$passPercent = "NA";}
	else {
		$passPercent = sprintf("%.1f", $passCount * 100 / $totalCount);
	}

	my $trueCount = $records[3];
	my $truePower = sprintf("%.1f", $trueCount * 100 / $bedCount);
	my $trueHomCount = $records[4];
	my $trueHomPercent;
	if ( $trueCount == 0 ) {
		$trueHomPercent = "NA";}
	else {
		$trueHomPercent = sprintf( "%.0f", $trueHomCount * 100 / $trueCount );}

	my $novelCount = $records[5];
	my $novelHomCount = $records[6];
	my $novelHomPercent;
	if ( $novelCount == 0 ) {
		$novelHomPercent = "NA";}
	else {
		$novelHomPercent = sprintf( "%.1f", $novelHomCount * 100 / $novelCount );}
	close IN;		

# clear intermediates
	my $rm_cmd = "rm -f $OverlapVcf $NovelVcf $FilterVcf $TmpLineCount";
	system( $rm_cmd );

	
# print out to std::cout
	if ( $rec == 0 ) {
		print "Filter: Qual >= $MinQual\n";
		print "MEI-type\t#pass\t%pass\t#true\t%power\t%1/1-true\t#novel\t%1/1-novel\n";
	}
	if ( $rec == 0 ) {
		print "Alu";}
	elsif( $rec == 1) {
		print "L1";}
	elsif ( $rec == 2 ) {
		print "SVA";}
	else {
		print "ALL";}
	print "\t$passCount\t$passPercent\t$trueCount\t$truePower\t$trueHomPercent\t$novelCount\t$novelHomPercent\n";
}

my $rm_merge_cmd = "rm -f $out_dir"."all.bed";
system( $rm_merge_cmd );














