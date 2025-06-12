bam=$1
outdir=$2

samtools view -H $bam | grep ^@SQ | /usr/bin/perl -a -F'\t' -ne 'print $1 . "\n" if $_ =~ /SN:((chr)?[1-9]?[0-9|X|Y])\s/' | while read chr;do
    echo $chr;
    samtools view -F 0x900 -h $bam ${chr}:1-999999999 | samtools view -Sb >$outdir/${chr}.bam;
done;
