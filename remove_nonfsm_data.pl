#!/usr/bin/perl

use warnings;
use strict;
use IO::File;
use File::Basename;
use File::Spec;

#arg 0 = read/transcript assignments: e.g. chr1_compatible_isoform.tsv
#arg 1 = espresso outputs directory containing SJ_simplified.list, subdirs 0, 1, 2, etc
#arg 2 = input tsv of bam files (e.g. /scratch1/fs1/timley/espresso/test_data/outputs/chr9/chr9.tsv.updated)

#get fsm transcripts, hash the read id
print STDERR "reading fsm read ids\n";
my %reads;
my $inFh = IO::File->new( $ARGV[0] ) || die "can't open file\n";
while( my $line = $inFh->getline )
{
    chomp($line);
    my @F = split("\t",$line);
    if($F[2] eq "FSM"){
        $reads{$F[0]} = 1;
    }
}
close($inFh);
print STDERR  "done\n";
my $dir=$ARGV[1];
unless (-e "$dir/fsm"){
    mkdir("$dir/fsm");
}

#these should probably be functions, not copy/paste, but meh
print STDERR "globbing $dir/*/*_read_final.txt\n";
#filter all of the read_final files, 
my @files = glob("$dir/*/*_read_final.txt");
foreach my $file (@files){
    print STDERR "filtering $file\n";

    my $subdir=basename(dirname($file));
    mkdir("$dir/fsm/$subdir");
    my $outFh = open (OUTFILE, ">" . "$dir/fsm/$subdir/" . basename($file)) || die "Can't open output file.\n";

    #read in the file, keep only reads that match the known fsm reads
    my $inFh2 = IO::File->new( $file ) || die "can't open file\n";
    while( my $line = $inFh2->getline )
    {
        chomp($line);
        my @F = split("\t",$line);
        if($reads{$F[0]}){
            print OUTFILE $line . "\n";
        }
    }
    close($inFh2);
    close(OUTFILE);
}

#filter all of the sam.list3 files
my @files = glob("$dir/*/sam.list3");
foreach my $file (@files){
    print STDERR "filtering $file\n";

    my $subdir=basename(dirname($file));
    mkdir("$dir/fsm/$subdir");
    my $outFh= open (OUTFILE, ">" . "$dir/fsm/$subdir/" . basename($file)) || die "Can't open output file.\n";

    #read in the file, keep only reads that match the known fsm reads
    my $inFh2 = IO::File->new( $file ) || die "can't open file\n";
    while( my $line = $inFh2->getline )
    {
        chomp($line);
        my @F = split("\t",$line);
        if($reads{$F[2]}){
            print OUTFILE $line . "\n";
        }
    }
    close($inFh2);
    close(OUTFILE);
}

#filter all of the sj.list files
my @files = glob("$dir/*/sj.list");
foreach my $file (@files){
    print STDERR "filtering $file\n";

    my $subdir=basename(dirname($file));
    mkdir("$dir/fsm/$subdir");
    my $outFh= open (OUTFILE, ">" . "$dir/fsm/$subdir/" . basename($file)) || die "Can't open output file.\n";

    #read in the file, keep only reads that match the known fsm reads
    my $inFh2 = IO::File->new( $file ) || die "can't open file\n";
    while( my $line = $inFh2->getline )
    {
        chomp($line);
        my @F = split("\t",$line);

        #check the 7th col
        unless($F[7] eq "NA"){
            my @names = split(",",$F[7]);
            my @namesout;
            foreach my $name (@names){
                if($reads{$name}){
                    push(@namesout,$name);
                }
            }
            $F[5] = scalar(@namesout);
            if($F[5] > 0){
                $F[7] = join(",",@namesout) . ","; #they use a trailing comma for some reason
            } else {
                $F[7] = "NA";
            }
        }            
        
        #check the 8th col
        unless($F[8] eq "NA"){
            my @names = split(",",$F[7]);
            my @namesout;
            foreach my $name (@names){
                if($reads{$name}){
                    push(@namesout,$name);
                }
            }
            $F[6] = scalar(@namesout);
            if($F[6] > 0){
                $F[8] = join(",",@namesout) . ","; #they use a trailing comma for some reason
            } else {
                $F[8] = "NA";
            }
        }
        
        print OUTFILE join("\t",@F) . "\n"
    }
    close($inFh2);
    close(OUTFILE);
}

#finally, filter each bam listed in the config file
unless (-e "$dir/fsm/bams"){
    mkdir("$dir/fsm/bams");
}

my @bams;
my @samp;
my @num;
#read config to get bams/samples
my $inFh2 = IO::File->new( $ARGV[2]) || die "can't open file\n";
while( my $line = $inFh2->getline )
{
    chomp($line);
    my @F = split("\t",$line);
    push(@bams,$F[0]);
    push(@samp,$F[1]);
    push(@num,$F[2]);
}

#write out a new config file with paths to the new bams as we go
my $outFh2= open (CONFIG, ">" . $dir . "/fsm/" . basename($ARGV[2])) || die "Can't open output config file.\n";

for(my $i=0;$i<@bams;$i++){
    my $samname=$dir . "/fsm/bams/" . $samp[$i] . ".sam";
    my $bamname=$dir . "/fsm/bams/" . $samp[$i] . ".bam";
    #print STDERR "reading from " . $bams[$i] . "\n";
    print STDERR "writing to $samname\n";
    my $outFh= open (OUTFILE, ">" . $samname) || die "Can't open output file.\n";
    open(BAM, "samtools view -h " . $bams[$i] . " |") or die "Can't open bam file";
    while(my $line = <BAM>){
        chomp($line);
        my @F = split("\t",$line);
        if($reads{$F[0]} || $line =~ /^@/){
            print OUTFILE $line . "\n"
        }
    }
    print CONFIG join("\t",File::Spec->rel2abs($bamname),$samp[$i],$num[$i]) . "\n";
    close(BAM);
    close(OUTFILE);
    system("samtools view -Sb -o $bamname $samname");
    unlink($samname);
}
