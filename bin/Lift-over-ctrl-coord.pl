#!/usr/bin/perl -w

# lift over ctrl REF_CHR coord

# method: read unlifted-MEI then store
#	 --> read whole coord file and print out
#	--> invoke Adjust-ref.sh to set right full path

use strict;
use Getopt::Long;

my $LiftList;
my $Coord;

GetOptions(
	"l=s" => \$LiftList,
	"c=s" => \$Coord
)
or die ( "ERROR in command line args!\n" );

# run bedtools first
my $exclude_overlap = "bedtools intersect -a $Coord -b $LiftList -v > $Coord.unoverlap";
system( $exclude_overlap );

my $OutName = "$Coord.liftOver";

my @DelStart;
my @DelEnd;
my @Lift;
open( IN, "<", $LiftList ) or die $!;
while(<IN>) {
	chomp; my $line = $_;
	my @s = split '\t',$line;
	push( @DelStart, $s[1] );
	push( @DelEnd, $s[2]);
	push( @Lift, $s[2] - $s[1] );
}
close IN;

for( my $i=1; $i < scalar(@Lift); $i++ ) {
	$Lift[ $i ] = $Lift[ $i ] + $Lift[ $i - 1 ];
}

open( IN, "<", "$Coord.unoverlap" ) or die $!;
open( OUT, ">", $OutName ) or die $!;
my $current_start = $DelStart[0];
my $current_end = $DelEnd[0];
my $index = 0;
my $lift = 0;
my $maxIndex = scalar( @Lift ) - 1;
my $reach_end = 0;
my $all_done = 0;
while(<IN>) {
	chomp; my $line = $_;
	if ( $all_done == 1 ) {
		print OUT "$line\n";
		next;
	}
# check chr
	my @pre = split '\t', $line, 2;
	if ( $pre[0] ne "20" ) {
		if ( $reach_end == 1 ) {
			$all_done = 1;
		}
		print OUT "$line\n";
		next;
	}

	my @s = split '\t', $pre[1];

	if ( $reach_end == 1 ) {
		my $st = $s[0] - $lift;
		my $ed = $s[1] - $lift;
		print OUT "$pre[0]\t$st\t$ed\t$s[2]\n";
		next;
	}

#print "$s[0]\t$current_start\t$current_end\t$lift\n";

	if ( $s[0] >= $current_end ) { # actually equal won't happen
		while( $s[0] >= $current_end ) {
			$index = $index + 1;
			if ( $index > $maxIndex ) {
				$lift = $Lift[ $maxIndex ];
				$reach_end = 1;
				last;
			}
			$current_end = $DelStart[ $index ];
			$current_start = $DelStart[ $index];
		}
#print "$Lift[ $maxIndex ]\t$current_start\t$current_end\t$lift\n";
		if( $reach_end == 0 ) {
			$current_start = $DelStart[ $index ];
			$current_end = $DelEnd[ $index ];
			$lift = $Lift[ $index - 1 ];
		}
		else {	# direct print
			my $st = $s[0] - $lift;
			my $ed = $s[1] - $lift;
			print OUT "$pre[0]\t$st\t$ed\t$s[2]\n";
			next;
		}
	}

# now it's before coord
	if ( $s[0] >= $current_start ) {
		die "Overlapping?\n";
	}

# do lift
	my $st = $s[0] - $lift;
	my $ed = $s[1] - $lift;
	print OUT "$pre[0]\t$st\t$ed\t$s[2]\n";
}
close IN;
close OUT;
