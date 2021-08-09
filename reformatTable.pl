#!/usr/bin/env perl

# This script sets all read counts with < 2.37% abundance in a sample to 0. OTU tables need to be transposed before using this script

use warnings;
use strict;

open (LIST, $ARGV[0]) or die "Cannot open input file: $!\n";

my @list;
my @outLines;

while (<LIST>) {
  chomp;
  @list = split(/\t/, $_);
  my $sum = eval join '+', @list;
  for my $i (0 .. $#list) {
    if ($list[$i] < ($sum * 0.0237)) {
    $list[$i] = 0; }
   }
  push(@outLines, "@list\n");
  }

close (LIST);

open (OUTFILE, ">$ARGV[1]") or die "Cannot open input file: $!\n";
print (OUTFILE @outLines);
close (OUTFILE);
