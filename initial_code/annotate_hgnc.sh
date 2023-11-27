#!/bin/bash

export HGNC=/diskmnt/Projects/Nanopore_analysis/rmashl/reference/annotation/HGNC/hgnc_complete_set.txt
export FILE="$1"

perl -e 'use Data::Dumper;
%m = ();
open(IN, "< $ENV{HGNC}");
while(<IN>) {
  chomp;
  next if /^symbol/;
  @a = split/\t/;
  $symbol = $a[1];
  $gene_id = $a[19];
  if( $gene_id ne "") {
     $m{ $gene_id } = $symbol;
   }
}
close(IN);
#print Dumper \%m;


open(IN, "< $ENV{FILE}");
while(<IN>) {
  chomp;
  @a = split/\t/;
  if( /^ids/) {  # header
    print join("\t", @a, "hgnc_symbol");
   next;
  }

  # require known gene or _chr
  if( $a[0] !~ /_chr/  &&  $a[0] !~ /ENSG/) { next }
  @b = split /_/, $a[0];

  @c = split /\./, $b[1];

  $gene_nonversioned = $c[0];
  $symbol = $gene_nonversioned;

  if( exists $m{ $gene_nonversioned } ) {  # overwrite if exists
      $symbol = $m{ $gene_nonversioned };
  }
  print join("\t", @a, $symbol),"\n";
}
close(IN);
'
