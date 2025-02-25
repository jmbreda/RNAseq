#!/usr/bin/env perl

# Based on scripts and discussion found in http://seqanswers.com/forums/showthread.php?t=1425
# and https://www.biostars.org/p/43855/

# This script changes the color space sequence on the 2nd line of each fastq
# to base space.  It does not convert the qual string since SOLID qual is
# already in Sanger format.  The first base and quality value (primer) are
# discarded.

use strict;
use warnings;


while (<>) {
    chomp ( my $id1  = $_ );
    chomp ( my $csfa = <> );
    chomp ( my $id2  = <> );
    chomp ( my $qual = <> );

    my $fa = csfa2fa($csfa);

    my $qual = substr($qual, 1);

    print join ("\n", $id1, $fa, $id2, $qual), "\n";
}


sub csfa2fa {
    my ($seq) = @_;
    my %cs = (
    "T0" => "T",
    "T1" => "G",
    "T2" => "C",
    "T3" => "A",
    "T." => "N",

    "C0" => "C",
    "C1" => "A",
    "C2" => "T",
    "C3" => "G",
    "C." => "N",

    "G0" => "G",
    "G1" => "T",
    "G2" => "A",
    "G3" => "C",
    "G." => "N",

    "A0" => "A",
    "A1" => "C",
    "A2" => "G",
    "A3" => "T",
    "A." => "N",

    "N0" => "N",
    "N1" => "N",
    "N2" => "N",
    "N3" => "N",
    "N." => "N",
    );

    my @letters = split( //, $seq );
    my $first_base = $letters[0];
    for( my $i = 1; $i < @letters ; $i++ ) {
        my $colour = $letters[$i];
        my $encoding = $first_base . $colour;
        $first_base = $cs{ $encoding };
        $letters[ $i ] = $first_base;
    }

    shift( @letters );
    $" = "";
    return join ("", @letters);
}
