#!/usr/bin/perl

use strict;
use warnings;

die unless @ARGV == 2;

my ($dir_in,$f_db_gtf)=@ARGV;

my $f = $dir_in."/rmats-long-fsm/differential_transcripts.tsv";
my $f_ab = $dir_in."/espresso/merged_fsm_N2_R2_abundance.esp";
my $f_ab_out=$dir_in."/summary/differential_transcripts.gn.nmd.abundance.tsv";
my $f_sqanti_nmd=$dir_in."/sqanti/merged_N2_R0_updated_classification.txt";

#$f="/storage1/fs1/dinglab/Active/Projects/scao/scor/analysis/run3/outputs/sqanti/merged_N2_R0_updated_classification.txt"; foreach $l (`cat $f`) { $ltr=$l; chomp($ltr); if($ltr=~/^isoform/) {  next; } else { @t=split("\t",$ltr); $nmd{$t[0]}=$t[36]; }} 

my %abundance;
my %mapenst;
my @header;
my @samples;
my %enstnmd; 
my %sqantinmd; 

foreach my $l (`cat $f_sqanti_nmd`) 
{ 
    my $ltr=$l; 
    chomp($ltr); 
    if($ltr=~/^isoform/) {  next; } 
    else { 
    my @t=split("\t",$ltr); $sqantinmd{$t[0]}=$t[36];
     }
} 

open my $ab_file, '<', $f_ab or die "Cannot open abundance file: $!";
while (my $l = <$ab_file>) {
    chomp($l);
    if ($l =~ /^transcript/) {
        @samples = split("\t", $l);
    } else {
        my @t = split("\t", $l);
        for (my $i = 3; $i < scalar @t; $i++) {
        
            $abundance{$t[0]}{$samples[$i]} = $t[$i];
        }
    }
}
close $ab_file;


open my $gtf_file, '<', $f_db_gtf or die "Cannot open GTF file: $!";
while (my $l = <$gtf_file>) {
    chomp($l);
    if ($l =~ /^#/) {
        next;
    } else {
        my @t = split("\t", $l);
        my $gid = "";
        my $gn ="";
        my $tid= "";
        my $nmd_status = "";
        my @t2 = split(/\;/, $t[8]);
        for (my $i = 0; $i < scalar @t2; $i++) {
            if ($t2[$i] =~ /gene_id/) {
                $gid = $t2[$i];
                $gid =~ s/gene_id //g;
                $gid =~ s/\"//g;
            }
            if ($t2[$i] =~ /transcript_id/) {
                $tid = $t2[$i];
                $tid =~ s/ transcript_id //g;
                $tid =~ s/\"//g;
            }
            if ($t2[$i] =~ /gene_name/) {
                $gn = $t2[$i];
                $gn =~ s/ gene_name //g;
                $gn =~ s/\"//g;
            }
            if ($t2[$i] =~ /transcript_biotype/) {
                $nmd_status = $t2[$i];
                $nmd_status =~ s/ transcript_biotype //g;
                $nmd_status =~ s/\"//g;
            } 
        }
        if ($gn ne "" && $gid ne "") {
            $mapenst{$gid} = $gn;
        }
        if ($nmd_status eq "nonsense_mediated_decay") {
         $enstnmd{$tid} = "TRUE";
        } 
        else {  $enstnmd{$tid} = "FALSE"; }
    }
}

open my $output_file, '>', $f_ab_out or die "Cannot open output file: $!";

open my $input_file, '<', $f or die "Cannot open input file: $!";

# gene_id feature_id      lr      df      pvalue  adj_pvalue      3688-DMSO-3_proportion  3936-DMSO-2_proportion  3959-DMSO-1_proportion  3688-SMG1-3_proportion  3936-SMG1-2_proportion  3959-SMG1-1_proportion  group_1_average_proportion      group_2_average_proportion      delta_isoform_proportion
# get the sample name from above header

## add gene name and abundance to the input file and output
 
while (my $l = <$input_file>) {
    chomp($l);
    my @t = split("\t", $l);
    if ($l =~ /^gene_id/) {
        print $output_file "NMD_By_Sqanti","\t","NMD_By_Ensembl","\t","gene_name";
        @header = split("\t",$l); 
        for (my $i = 0; $i < 6; $i++) {
            print $output_file "\t", $t[$i];
        }
		for(my $i=6;$i< scalar @t - 3;$i++)
		{
			my $header_abs=$t[$i];	
			$header_abs=~s/_proportion/_abundance/g;
       	 	print $output_file "\t", $header_abs;
		}

        for (my $i = 6; $i < scalar @t; $i++) {
            print $output_file "\t", $t[$i];
        }
        print $output_file "\n";
    } else {
        my @t3 = split(/\,/, $t[0]);
        my $gn = $mapenst{$t3[0]};

        for (my $i = 1; $i < scalar @t3; $i++) {
            $gn = $gn . "," . $mapenst{$t3[$i]};
        }


        if(defined $enstnmd{$t[1]}) {
        print $output_file $sqantinmd{$t[1]},"\t",$enstnmd{$t[1]},"\t",$gn;
        }
        else 
        {
        print $output_file $sqantinmd{$t[1]},"\t", "NA","\t",$gn; 
        }
        
        for (my $i = 0; $i < 6; $i++) {
            print $output_file "\t", $t[$i];
        }
		for(my $i=6;$i< scalar @t - 3;$i++)
		{
			my $header_abs=$header[$i];	
			$header_abs=~s/_proportion//g;	
        	print $output_file "\t", $abundance{$t[1]}{$header_abs}; 
        }
		for(my $i=6;$i<scalar @t;$i++) 
		{ print $output_file "\t",$t[$i]; } 
         print $output_file "\n"; 
         }
	}	
	

