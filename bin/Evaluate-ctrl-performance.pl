#!/usr/bin/perl -w

# evaluate Morphling performance from ctrlVcf

use strict;
use Getopt::Long;
use English;

if ( !(@ARGV) ) {
        print "\nUsage: [] is defulat value\n";
        print "  -i   Ctrl vcf name.\n";
        print "  -m   MEI name. [NA]\n";
        print "  -o   Output result file.\n";
        print "  -r   Reference bed. \n";
        print "  -e   Exclude reference bed. \n";
        print "  -t   Working directory. [/tmp]\n";
        print "  -q   Quality filter for Morphling discovery. [10]\n";
        print "  --append  Append to file indicated in -o. [off]\n";
        print "\n";
        exit 0;
}

my $VCF;
my $MEI;
my $LOG;
my $BED;
my $EXCLUDE;
my $WDIR = "/tmp";
my $QUAL = 10;
my $APPEND = 0;
GetOptions(
        'i=s' => \$VCF,
        'm=s' => \$MEI,
        'o=s' => \$LOG,
        'r=s' => \$BED,
        'e=s' => \$EXCLUDE,
        't=s' => \$WDIR,
        "q=s" => \$QUAL,
        "append" => \$APPEND,
);

if ( length($VCF) == 0 ) {
	print "ERROR: please have a vcf in option -i !\n";
	exit 1;
}
if ( length($MEI) == 0 ) {
	print "Warning: no valid MEI name. Use NA instead !\n";
	$MEI = "NA";
}
if ( length($LOG) == 0 ) {
	print "ERROR: please have an output file in option -o !\n";
	exit 1;
}
if ( length($BED) == 0 ) {
	print "ERROR: please have a reference bed in option -r !\n";
	exit 1;
}
if ( length($EXCLUDE) == 0 ) {
	print "ERROR: please have a exclude reference bed in option -e !\n";
	exit 1;
}

# fix wdir
my $last_char = substr $WDIR, length($WDIR) - 1, 1;
if ( $last_char eq '/' ) {
	 $WDIR = substr $WDIR, 0, length($WDIR)-1;
}


# run bedtools
my $res = CalculateSinglePerformance( $VCF, $BED, $EXCLUDE );
if ( $APPEND > 0 ) {
	open( OUT, ">>", $LOG ) or die $!;
	print OUT "$MEI\t$res\n";
	close OUT;
}
else {
	open( OUT, ">", $LOG ) or die $!;
	print OUT "MEI-type\t%Power\t%Ctrl-1/1\t#Novel\t%Novel-1/1\n";
	print OUT "$MEI\t$res\n";
	close OUT;
}


sub CalculateSinglePerformance
{
	my $vcf = shift;
	my $ref = shift;
	my $exclude = shift;
	
  #general filter first
	my $FilterBed = "$WDIR/$PID.filter.bed";
	my $grep_cmd = "cat $vcf | awk -v qual=$QUAL \'\$6 >= qual\' | grep PASS | awk '{print \$1,\$2-600,\$2+600,\$10}' | tr ' ' '\t' | bedtools sort -i - > $FilterBed";
    system( $grep_cmd );
    
  # count
    my $TmpCountFile = "$WDIR/$PID.count";
    my $IntersectBed = "$WDIR/$PID.intersect";
    my $NovelBed = "$WDIR/$PID.novel";
    my $intersect_cmd = "cat $FilterBed | bedtools intersect -a - -b $ref -u > $IntersectBed";
#print "$intersect_cmd\n";
    system( $intersect_cmd );
    my $novel_cmd = "cat $FilterBed | bedtools intersect -a - -b $exclude -v > $NovelBed";
    system( $novel_cmd );
    my $wc_cmd = " cat $IntersectBed | wc -l > $TmpCountFile; cat $IntersectBed | grep \"1\/1\" | wc -l >> $TmpCountFile; ";
    $wc_cmd = $wc_cmd." cat $NovelBed | wc -l >> $TmpCountFile; cat $NovelBed | grep \"1\/1\" | wc -l >> $TmpCountFile;";
    $wc_cmd = $wc_cmd."cat $FilterBed | bedtools intersect -a $ref -b - -u | wc -l >> $TmpCountFile;";
    $wc_cmd = $wc_cmd."cat $ref | wc -l >> $TmpCountFile;";
#print "$wc_cmd\n";
    system( $wc_cmd );
    open( IN, "<", $TmpCountFile ) or die $!;
    my @records;
    while( <IN> ) {
        chomp;
        my $line = $_;
        push( @records, $line );
    }
    close IN;
# sanity check
    if ( @records != 6 ) {
    	die "ERROR: Unexpected records in $TmpCountFile!\n";
    }
    my $refcount = $records[5];
    my $power;
    if ( $refcount == 0 ) {
    	print "Warning: [Evaluate-ctrl-performance.pl] Empty reference bed: $ref!\n";
    	$power = "NA";
    }
    else {
    	$power = sprintf( "%.1f", $records[4]*100 / $refcount );
    }
    my $homPercent;
    if ( $records[0] == 0 ) {
    	$homPercent = 0;
    }
    else {
    	$homPercent = sprintf( "%.1f", $records[1]*100 / $records[0] );
    }
    my $novelHomPercent;
    if ( $records[2] == 0 ) {
        $novelHomPercent = "NA";
    }
    else {
        $novelHomPercent = sprintf( "%.1f", $records[3]*100 / $records[2] );
    }
    my $result = "$power\t$homPercent\t$records[2]\t$novelHomPercent";
# clean & return
    my $rm_cmd = "rm -f $FilterBed $TmpCountFile $NovelBed $IntersectBed";
    system( $rm_cmd );
    return $result;
}
























































