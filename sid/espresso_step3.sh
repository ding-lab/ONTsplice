#!/bin/bash
input_tsv=$1
outdir=$2
logdir=$3
reference_fasta=$4  #/storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa
gtf=$5 #/storage1/fs1/bga/Active/gmsroot/gc2560/core/GRC-human-build38_human_95_38_U2AF1_fix/rna_seq_annotation/Homo_sapiens.GRCh38.95.gtf

base=$(basename $input_tsv .tsv)
mem="200G"
threads=4

#doesn't actually need fasta, but keep args the same for concordance with other commands
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-timley -J espresso -oo ${logdir}/${base}_step3.out -eo ${logdir}/${base}_step3.err -n 8 -R"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" -M ${mem} -q timley -a 'docker(sridnona/espresso:v2)' /bin/bash -c "source activate env && perl /bin/espresso/src/ESPRESSO_Q.pl -L ${outdir}/${base}/${base}.tsv.updated -A ${gtf} -V ${outdir}/${base}/${base}_compatible_isoform.tsv -T ${threads}"
