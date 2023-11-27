#!/usr/bin/perl

use strict;
use warnings;

die unless @ARGV == 2;

my ($sample_full_path, $chr_status)=@ARGV;

my $f_esp_in; 
my $f_gtf_in; 
my $l; 
my $ltr; 
## merge calls from different chromosomes ##
my $f_esp_out=$sample_full_path."/merged_N2_R0_abundance.esp"; 
my $f_gtf_out=$sample_full_path."/merged_N2_R0_updated.gtf"; 

my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

open(OUT_esp,">$f_esp_out"); 
open(OUT_gtf,">$f_gtf_out"); 

my $first_esp=0; 
my $first_gtf=0; 

foreach my $chr (@chrlist)
	{
    	my $chr1=$chr;
   
        if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
		$f_esp_in=$sample_full_path."/".$chr1."/bam_N2_R0_abundance.esp"; 
		$f_gtf_in=$sample_full_path."/".$chr1."/bam_N2_R0_updated.gtf"; 

#		open(IN_esp,"<$f_esp_in");

		foreach $l (`cat $f_esp_in`)
		{
			$ltr=$l; 
			chomp($ltr); 
			#while(<IN_esp>)
			#{
				if($ltr=~/^transcript_ID/)
				{
					if($first_esp==0) {
					$first_esp=1; 
					print OUT_esp $ltr,"\n"; } 				
				}
				else { print OUT_esp $ltr,"\n"; }
			#}
		}
		$first_gtf=0; 
		#open(IN_gtf,"<$f_gtf_in");
		foreach $l (`cat $f_gtf_in | sort -k 1,1 -k 4,4n |  perl -nae \'print \$_ unless \$F[3] > \$F[4]\'`)
		{
			$ltr=$l; 
			chomp($ltr); 
			#while(<IN_gtf>)
			#{
				if($ltr=~/^# ESPRESSO/)
				{
					if($first_gtf==0)
					{
					$first_gtf=1; 
					print OUT_gtf $ltr,"\n"; 
					}				
				}
				else { print OUT_gtf $ltr,"\n"; }
			#}
		}
        }
	}


