#!/bin/bash

annotate_hgnc.sh counts_matrix.tsv.tpm  > counts_matrix.tsv.tpm.anno
tail -n +2 counts_matrix.tsv.tpm.anno | cut -f3  | sort -u   > counts_matrix.tsv.tpm.anno.genes_uniq

echo  genes_with_names: $( grep -c -v ^chr counts_matrix.tsv.tpm.anno.genes_uniq)
echo  genes_without_names: $( grep -c  ^chr counts_matrix.tsv.tpm.anno.genes_uniq)

