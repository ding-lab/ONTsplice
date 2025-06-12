#!/bin/bash
config=$1
outdir=$2
logdir=$3
espresso_abundance=$4 #merged_N2_R0_abundance.esp
espresso_gtf=$5 #merged_N2_R0_updated.gtf
gtf=$6 #/storage1/fs1/bga/Active/gmsroot/gc2560/core/GRC-human-build38_human_95_38_U2AF1_fix/rna_seq_annotation/Homo_sapiens.GRCh38.95.gtf

# config file should contain group \t sample_name
# where group is either "g1" or "g2"
# todo - make this part of an initial config that also contains bam files, etc

mem="64G"
threads=16

mkdir -p ${outdir}
grep "^g1" $config | cut -f 2 | perl -pe 's/\n/,/g' | perl -pe 's/,$//g' >${outdir}/group1.txt
grep "^g2" $config | cut -f 2 | perl -pe 's/\n/,/g' | perl -pe 's/,$//g' >${outdir}/group2.txt

LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-timley -q general -J rmats -n $threads -R"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" -M ${mem} -oo ${logdir}/rmats-long_step7.out -eo ${logdir}/rmats-long_step7.err -a "docker(sridnona/rmats_long:v3)" /bin/bash -c "source activate /docker_data/rMATS-long/conda_env && cd /docker_data/rMATS-long/scripts/ && python rmats_long.py \
--abundance ${espresso_abundance} --updated-gtf ${espresso_gtf} \
--gencode-gtf ${gtf} \
--group-1 ${outdir}/group1.txt \
--group-2 ${outdir}/group2.txt \
--group-1-name group1 \
--group-2-name group2 \
--out-dir $outdir \
--num-threads 16 \
--delta-proportion 0.1 \
--adj-pvalue 0.1 \
--plot-file-type .pdf"
