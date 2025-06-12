#!/bin/bash
input_tsv=$1 #only used to get chromosome name
outdir=$2
logdir=$3

base=$(basename $input_tsv .tsv)
mem="200G"
threads=8
mkdir -p ${outdir}/${base}/filtered

#using read_ratio_cutoff of 2 (1 is the highest possible) means that it will only quantify transcripts
#that are in the GTF
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-timley -J espresso -oo ${logdir}/${base}_step5.out -eo ${logdir}/${base}_step5.err -n $threads -R"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" -M ${mem} -q general -a 'docker(sridnona/espresso:v2)' /bin/bash -c "source activate env && perl /bin/espresso/src/ESPRESSO_Q.pl -L ${outdir}/${base}/${base}.tsv.updated -O ${outdir}/${base}/filtered -A ${outdir}/${base}/${base}_N2_R0_updated.filtered.gtf --read_ratio_cutoff 2 -V ${outdir}/${base}/filtered/${base}_compatible_isoform.tsv -T $threads"
