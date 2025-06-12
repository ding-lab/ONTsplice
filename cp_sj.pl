#!/usr/bin/perl

use strict;
use warnings;

die "Usage: $0 file_list out_dir\n" unless @ARGV == 2;

my ($f_sj, $dir_out) = @ARGV;
my $count = 0;

open(my $fh, "<", $f_sj) or die "Cannot open $f_sj: $!";
while (my $f = <$fh>) {
    chomp($f);
    my $f_out = "$dir_out/$count.SJ.out.tab";
    if (-e $f) {
        system("cp", $f, $f_out);
    }
    $count++;
}
close $fh;