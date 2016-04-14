#!/usr/bin/perl -w

# add more filter to Morphling calling results
# read filter criteria from filter setting file

use strict;
use Getopt::Long;
use Cwd 'abs_path';

if ( !@ARGV) {
  print "-v    input vcf\n";
  print "-c    filter configuration file. Default in Morphling/refs/filter.config\n";
  print "-s    sample list used for calling\n";
  print "\n";
  exit 0;
}

my $vcfname;
my $configname = "";
my $listname;

GetOptions (
  'v=s' => \$vcfname,
  'c=s' => \$configname,
  's=s' => \$listname,
);

if ( $listname eq "" || $vcfname eq "" ) {
  die "Missing required options!\n";
}

# get path first
my $pl_path = abs_path($0);
my @pls = split '/', $pl_path;
for(my $i=0; $i<@pls-2; $i++) {
  if ($pls[$i] eq "") {
    next;
  }
  $configname = "$configname/$pls[$i]";
}
$configname = "$configname/refs/filter.config";

# read configuration
# configuration format: info_name, lower_limit, upper_limit
my @info_name;
my @upper_limit;
my @lower_limit;
open(IN, "<", $configname) or die $!;
while(<IN>) {
  chomp;
  my $line = $_;
  my @s = split '\t', $line;
  if (@s != 3) {
    die "ERROR: configuration file doesn't have 3 columns at:\n$line\n";
  }
  my $iname = uc $s[0];
  push(@info_name, $iname);
  push(@lower_limit, $s[1]);
  push(@upper_limit, $s[2]);
}
close IN;

# apply to vcf
my $outname = "$vcfname.filtered.vcf";
open(OUT, ">", $outname) or die $!;
open(IN, "<", $vcfname) or die $!;
my $pass_header = 0;
while(<IN>) {
  chomp;
  my $line = $_;
  if ($pass_header == 0) {
	my @s0 = split "", $line, 2;
	if ($s0[0] ne "#") {
		$pass_header = 1;
	}
	else {
		next;
	}
  }
  my @s = split '\t', $line;
  my @info_fields = split ';', $s[7];
  my %imap;
  foreach my $ifield (@info_fields) {
    my @fs = split '=', $ifield;
    my @commas = split ',', $fs[1];
    my $val;
    if (@commas ==1) {
      $val = $fs[1];
    }
    elsif (@commas == 2) {
      $val = $commas[0] + $commas[1];
     }
     else {
      die "ERROR: bad info field at: $ifield\n";
     }
     $imap{$fs[0]} = $val;
  }
  # then check and change filter name
  my $filter = "";
  for(my $i=0; $i<@info_name; $i++) {
    my $key = $info_name[$i];
    if ( !exists($imap{$key}) ) {
      die "ERROR: cannot find filter $key in vcf info field!\n";
    }
    if ($lower_limit[$i] ne "-") {
      if ($imap{$key} < $lower_limit[$i]) {
        $filter = $key;
        last;
      }
    }
    if ($upper_limit[$i] ne "-") {
      if ($imap{$key} > $upper_limit[$i]) {
        $filter = $key;
        last;
      }
    }
  }
  if ($filter ne "") { # change filter name
    $s[6] = $filter;
  }
  # print out at last
  print OUT "$s[0]";
  for(my $i=1; $i<@s; $i++) {
    print OUT "\t$s[$i]";
  }
  print OUT "\n";
}
close IN;
close OUT;
