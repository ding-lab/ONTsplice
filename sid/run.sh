#This differs because it uses shortread data in the sqanti filtering step to help distinguish true isoforms

config=$1 #3 column file: group (either g1 or g2),  sample_id,  /path/to/bam
scratch=$2 #scratch directory for outputs: /scratch1/fs1/timley/espresso/tumor.vs.cd34.freeze.v1
scripts=$3  #path where espresso scripts are stored
gtf=$4  #reference gtf  
#human: /storage1/fs1/bga/Active/gmsroot/gc2560/core/GRC-human-build38_human_95_38_U2AF1_fix/rna_seq_annotation/Homo_sapiens.GRCh38.95.gtf
#mouse: /storage1/fs1/bga/Active/gmsroot/gc2560/core/GRC-mouse-build38_mouse_95_38/rna_seq_annotation/Mus_musculus.GRCm38.95.gtf
reffasta=$5 #reference genome fasta  
#human: /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa
#mouse: /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/01b6db5e23924e738efc3e398838bf12/build0b5a8d3795254ae9bd6fa80420656257/all_sequences.fa
star_sjouts=$6 #file containing single column list of star SJ_out file paths
star_bams=$7 #file containing single column list of star bam file paths

mkdir -p $scratch/inputs
mkdir -p $scratch/outputs
mkdir -p $scratch/logs

echo "splitting inputs"
#split inputs so that things can be run per-chromosome, also "sanitize" the bams to remove secondary reads, etc
cat $config | while read group samp bam;do 
    mkdir $scratch/inputs/$samp;
    bsub -J split -oo $scratch/logs/step0_split_${samp}.log -G compute-timley -q siteman -M 4G -R 'select[mem>4G] span[hosts=1] rusage[mem=4G]' -a "docker(chrisamiller/genomic-analysis:0.2)" bash $scripts/split_and_filter_bam.sh $bam $scratch/inputs/$samp
done
/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n split -d 30

echo "creating per-chr input files"
#create per-chromosome espresso input files
cut -f 2 $config | while read samp;do 
    ls -1 $scratch/inputs/$samp | while read bam;do 
        echo "$scratch/inputs/$samp/$bam $samp" | awk 'OFS="\t"{print $1,$2}' >>$scratch/inputs/$(basename $bam .bam).tsv;
    done;
done

echo "running espresso_s"
#run step1 - espresso_s
for i in $scratch/inputs/*.tsv;do 
    bash $scripts/espresso_step1.sh $i $scratch/outputs $scratch/logs $reffasta $gtf;
done
/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n espresso

echo "running espresso_c"
#run step2 - espresso_c
for i in $scratch/inputs/*.tsv;do 
    bash $scripts/espresso_step2.sh $i $scratch/outputs $scratch/logs $reffasta $gtf;
done
/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n espresso

echo "running espresso_q"
#run step3 - espresso_q
for i in $scratch/inputs/*.tsv;do 
    bash $scripts/espresso_step3.sh $i $scratch/outputs $scratch/logs $reffasta $gtf;
done
/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n espresso

echo "merging outputs"
headerchr=$(basename $(ls -1 $scratch/inputs/*.tsv | head -n 1) .tsv)
#merge outputs from all chromosomes into abundance/gtf files
head -n 1 $scratch/outputs/${headerchr}/${headerchr}_N2_R0_abundance.esp >$scratch/outputs/merged_N2_R0_abundance.esp #grab header
for i in $scratch/outputs/*/*_N2_R0_abundance.esp;do tail -n +2 $i >>$scratch/outputs/merged_N2_R0_abundance.esp;done 

head -n 1 $scratch/outputs/${headerchr}/${headerchr}_N2_R0_updated.gtf >$scratch/outputs/merged_N2_R0_updated.gtf
for i in $scratch/outputs/*/*_N2_R0_updated.gtf;do tail -n +2 $i | sort -k 1,1 -k 4,4n | perl -nae 'print $_ unless $F[3] > $F[4]' >>$scratch/outputs/merged_N2_R0_updated.gtf;done

echo "running squanti"
#run sqanti QC to get transcript metrics
mkdir $scratch/outputs/sqanti
#starbams=$(cat $star_bams | perl -pe 's/\n/,/g' | perl -pe 's/,$//g') #create comma-sep input

#when we get above a few dozen bams, this argument string gets too long for the parser. don't use 
#this comma sep, instead copy them to the same temp directory and pass that
starsjs=$(cat $star_sjouts | perl -pe 's/\n/,/g' | perl -pe 's/,$//g') #create comma-sep input
mkdir $scratch/tmp
count=1;
cat $star_sjouts | while read i;do 
    cp $i $scratch/tmp/$count.SJ.out.tab;
    count=$((count+1));
done
star_sjouts=$scratch/tmp
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-timley -J sqanti -oo $scratch/logs/sqanti_qc.log -n 8 -R"select[mem>32G] rusage[mem=32G] span[hosts=1]" -M 32G -q siteman -a "docker(chrisamiller/sqanti3:v5.2.1)" /bin/bash -c "source activate SQANTI3.env && /app/SQANTI3-5.2.1/sqanti3_qc.py -t 8 -d $scratch/outputs/sqanti -c $scratch/tmp --SR_bam $star_bams $scratch/outputs/merged_N2_R0_updated.gtf $gtf $reffasta"
/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n sqanti -s 60 -d 30

#run sqanti filtering using the machine-learning approach
LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-timley -J sqanti -oo $scratch/logs/sqanti_filter.log -n 8 -R"select[mem>32G] rusage[mem=32G] span[hosts=1]" -M 32G -q siteman -a "docker(chrisamiller/sqanti3:v5.2.1)" /bin/bash -c "source activate SQANTI3.env && /app/SQANTI3-5.2.1/sqanti3_filter.py ml $scratch/outputs/sqanti/merged_N2_R0_updated_classification.txt --gtf $scratch/outputs/sqanti/merged_N2_R0_updated_corrected.gtf -d $scratch/outputs/sqanti/filter_ml_default"
/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n sqanti -s 60 -d 30

echo "rerun espresso_q"
#take the filtered sqanti gtf, break it back down by chromosome
ls -1 $scratch/inputs/*.tsv | while read i;do
    chr=$(basename $i .tsv);
    grep "^$chr[[:space:]]" $scratch/outputs/sqanti/filter_ml_default/merged_N2_R0_updated.filtered.gtf >$scratch/outputs/${chr}/${chr}_N2_R0_updated.filtered.gtf

    #now run ESPRESSO_Q again, to requantify the data (restricting it to only these filtered transcripts)
    bash $scripts/espresso_step5_requant.sh $scratch/inputs/$chr.tsv $scratch/outputs $scratch/logs
done
/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n espresso

echo "merging outputs, removing non fsm reads"
#merge per-chr filtered abundance/gtf outputs from espresso
headerchr=$(basename $(ls -1 $scratch/inputs/*.tsv | head -n 1) .tsv)
mkdir -p $scratch/outputs/filtered/
head -n 1 $scratch/outputs/${headerchr}/filtered/${headerchr}_N2_R2_abundance.esp >$scratch/outputs/filtered/merged_N2_R2_abundance.esp
for i in $scratch/outputs/*/filtered/*_N2_R2_abundance.esp;do tail -n +2 $i >>$scratch/outputs/filtered/merged_N2_R2_abundance.esp;done
head -n 1 $scratch/outputs/${headerchr}/filtered/${headerchr}_N2_R2_updated.gtf >$scratch/outputs/filtered/merged_N2_R2_updated.gtf
for i in $scratch/outputs/*/filtered/*_N2_R2_updated.gtf;do tail -n +2 $i >>$scratch/outputs/filtered/merged_N2_R2_updated.gtf;done

#running espresso on the filtered gtf told us whether every read was FSM, ISM, etc.  Take that information and filter
#the espresso intermediate files so that they only contain FSM reads. Basically, copy the outputs directories and refill them
#with the same data with non-FSM reads removed.
ls -1 $scratch/inputs/*.tsv | while read i;do
    chr=$(basename $i .tsv);
    #filter this chromosome to remove all non-FSM reads from output tmp files, and subset the bams as well:
    bsub -G compute-timley -J fsmfilt -oo $scratch/logs/${chr}_fsm_filter.log -n 1 -R"select[mem>16G] rusage[mem=16G] span[hosts=1]" -M 16G -q siteman -a "docker(chrisamiller/docker-genomic-analysis:latest)" perl $scripts/remove_nonfsm_data.pl $scratch/outputs/${chr}/filtered/${chr}_compatible_isoform.tsv $scratch/outputs/${chr} /$scratch/outputs/${chr}/${chr}.tsv.updated
done
/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n fsmfilt -d 30

echo "running espresso_q with fsm-filtered files"
#okay, rerun espresso_q (again!) per chromosome, using these filtered files as inputs
for i in $scratch/inputs/*.tsv;do
    bash $scripts/espresso_step6_requantfsm.sh $i $scratch/outputs $scratch/logs $reffasta
done
/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n espres

echo "merging outputs"
#merge the per-chr outputs 
head -n 1 $scratch/outputs/${headerchr}/fsm/${headerchr}_N2_R2_abundance.esp >$scratch/outputs/merged_fsm_N2_R2_abundance.esp
for i in $scratch/outputs/*/fsm/*_N2_R2_abundance.esp;do tail -n +2 $i;done | sort >>$scratch/outputs/merged_fsm_N2_R2_abundance.esp
head -n 1 $scratch/outputs/${headerchr}/fsm/${headerchr}_N2_R2_updated.gtf >$scratch/outputs/merged_fsm_N2_R2_updated.gtf
for i in $scratch/outputs/*/fsm/*_N2_R2_updated.gtf;do tail -n +2 $i;done | sort -k 1,1 -k 4,4n >>$scratch/outputs/merged_fsm_N2_R2_updated.gtf

echo "running rmats"
#run step 7 - rmats diff isoform
mkdir -p $scratch/rmats-long-fsm
bash $scripts/espresso_step7_rmats.sh $config $scratch/rmats-long-fsm $scratch/logs $scratch/outputs/merged_fsm_N2_R2_abundance.esp $scratch/outputs/merged_fsm_N2_R2_updated.gtf $gtf


#create combined bams
mkdir $scratch/outputs/fsm_bams
cat $scratch/inputs/${headerchr}.tsv | while read bam samp;do
    bsub -G compute-timley -J merge -oo $scratch/logs/fsm_merge_bams_${samp}.log -n 1 -R"select[mem>16G] rusage[mem=16G] span[hosts=1]" -M 16G -q siteman -a "docker(chrisamiller/docker-genomic-analysis:latest)" "samtools merge $scratch/outputs/fsm_bams/${samp}.bam $scratch/outputs/*/fsm/bams/$samp.bam"
done
/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n merge -d 30
cat $scratch/inputs/${headerchr}.tsv | while read bam samp;do
    bsub -G compute-timley -J index -oo $scratch/logs/fsm_index_bams_${samp}.log -n 1 -R"select[mem>4G] rusage[mem=4G] span[hosts=1]" -M 4G -q siteman -a "docker(chrisamiller/docker-genomic-analysis:latest)" samtools index $scratch/outputs/fsm_bams/${samp}.bam
done

#extract novel transcripts
grep ESPRESSO $scratch/outputs/merged_fsm_N2_R2_updated.gtf >$scratch/outputs/merged_fsm_N2_R2_updated.novel_transcripts.gtf

echo "Done - remember to move data off of scratch"
