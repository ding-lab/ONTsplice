# ONTsplice version 1.2 , compute1 #

Detect differential transcripts between case and control from ONT data. 

## Usage ##

Step 0: set environment for LSF job on compute1 by adding the following to ~/.bashrc file: 

export PATH=/storage1/fs1/songcao/Active/Software/anaconda3/bin:$PATH

export STORAGE2=/storage1/fs1/dinglab/Active
export SCRATCH2=/storage1/fs1/dinglab/

export STORAGE1=/storage1/fs1/songcao/Active
export SCRATCH1=/storage1/fs1/songcao/

export LSF_DOCKER_VOLUMES="$STORAGE1:$STORAGE1 $STORAGE2:$STORAGE2"

Step1: Enter the directory where you downloaded ONTsplice pipeline 

Step2: Type the coommand line: perl ontsplice.ps.pl  --rdir --ref --log --q --groupname --users --step

[0] pre-process bam

[1] Split bam

[2]  Run espresso 1

[3]  Run espresso 2

[4]  Run espresso 3

[5]  Merge espresso results from different chrs

[6]  run sqanti

[7]  re-run step espresso 3 by using new filtered gtf file from sqanti

[8]  merge filtered espresso results

[9]  run fsm 

[10]  re-run step espresso 3 based fsm results

[11]  merge per chr results from 10 espresso results

[12] run rmats

[13] generate summary report


Example for running scripts:

perl ontsplice.ps.pl --rdir /storage1/fs1/dinglab/Active/Projects/scao/scor/analysis/run3 --log  /storage1/fs1/dinglab/Active/Projects/scao/scor/analysis/run3.log -ref /storage1/fs1/songcao/Active/Database/hg38_database/sid.ont/all_sequences.fa --q general --users songcao --groupname SomaticWXS --step 1

