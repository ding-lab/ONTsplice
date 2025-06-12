######### Song Cao###########
##### email: scao@wustl.edu ####
## pipeline for detecting splicing form from ONT data ##
#	ontsplice.pl #


#!/usr/bin/perl
##!/gscmnt/gc2525/dinglab/rmashl/Software/perl/perl-5.22.0/bin/perl
use strict;
use warnings;
#use POSIX;
use Getopt::Long;

my $version = 1.0;

#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m";
#usage information

(my $usage = <<OUT) =~ s/\t+//g;
Somatic variant calling pipeline 
Pipeline version: $version

$yellow     
Usage: perl $0 --rdir --ref --log --q --groupname --users --step 

$normal

<rdir> = full path of the folder holding files for this sequence run (user must provide)
<log> = full path of the folder for saving log file; usually upper folder of rdir
<groupname> = job group name
<users> = user name for job group
<step> run this pipeline step by step. (user must provide)
<ref> the human reference: 
<q> which queue for submitting job; research-hpc, ding-lab, long (default)

hg38: /storage1/fs1/songcao/Active/Database/hg38_database/GRCh38.d1.vd1/GRCh38.d1.vd1.fa
$purple [0] pre-process bam
$yellow [1] Split bam
$green [2]  Run espresso 1
$green [3]  Run espresso 2
$green [4]  Run espresso 3
$yellow [5]  Merge espresso results from different chrs
$yellow [6]  run sqanti based on short read bam and junction files
$yellow [7]  re-run step espresso 3 by using new filtered gtf file from sqanti
$yellow [8]  merge filtered espresso results
$cyan [9]  run fsm 
$cyan [10]  re-run step espresso 3 based fsm results
$cyan [11]  merge per chr results from 10 espresso results
$red [12] run rmats
$normal [13] generate summary report

OUT

#__DEFAULT NUMBER OF BINS IE (MUST BE INTEGER)
my $step_number = -1;
my $status_rerun=0; 

#__HELP (BOOLEAN, DEFAULTS TO NO-HELP)
my $help = 0;
my $q_name="";
#__FILE NAME (STRING, NO DEFAULT)
my $run_dir="";
my $log_dir="";
my $h38_REF="";
my $db_smg="";
#my $ref_name="";
my $chr_status=0;
my $compute_username="";
my $group_name="";

#__PARSE COMMAND LINE
my $status = &GetOptions (
      "step=i" => \$step_number,
      "groupname=s" => \$group_name,
      "users=s" => \$compute_username,	
      "rdir=s" => \$run_dir,
      "ref=s"  => \$h38_REF,
      "log=s"  => \$log_dir,
      "q=s" => \$q_name,
      "log=s"  => \$log_dir,
      "help" => \$help, 
	);
 
#print $status,"\n";

if ($help || $run_dir eq "" || $log_dir eq "" || $group_name eq "" || $compute_username eq "" || $h38_REF eq "" || $step_number<0 ) {
	 print "wrong option\n";
	  print $usage;
      exit;
   }

print "run dir=",$run_dir,"\n";
print "log dir=",$log_dir,"\n";
print "step num=",$step_number,"\n";
print "status rerun=",$status_rerun,"\n";
print "queue name=",$q_name,"\n";
print "job group=",$group_name,"\n";
print "user group=",$compute_username,"\n";


if($q_name eq "") 
{
	$q_name="long";
}


if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}

die $usage unless ($step_number >=0)&&(($step_number <= 16));
my $email = "scao\@wustl\.edu";

my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-1];
my $HOME1=$log_dir;

print $HOME1,"\n";

if (! -d $HOME1)
{
`mkdir $HOME1`; 
}
if (! -d $HOME1."/tmpONT") {
    `mkdir $HOME1"/tmpONT"`;
}
my $job_files_dir = $HOME1."/tmpONT";
#store SGE output and error files here
if (! -d $HOME1."/LSF_DIR_ONT") {
    `mkdir $HOME1"/LSF_DIR_ONT"`;
}
my $lsf_file_dir = $HOME1."/LSF_DIR_ONT";
#GENOMEVIP_SCRIPTS=/gscmnt/gc2525/dinglab/rmashl/Software/bin/genomevip
# obtain script path
#my $script_dir="/gscuser/scao/scripts/git/somaticwrapper";
my $run_script_path =`echo \$PWD`;
chomp $run_script_path;
#my $run_script_path = `dirname $0`;
my $script_dir=$run_script_path; 
#print $script_dir,"\n";

#my $run_script_path=$script_dir; 
#chomp $run_script_path;

#$run_script_path = "/gscmnt/gc2525/dinglab/rmashl/Software/perl/perl-5.22.0/bin/perl ".$run_script_path."/";
my $run_perl_script_path = "/usr/bin/perl ".$run_script_path."/";
my $run_py_script_path = "python ".$run_script_path."/";

my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";


### running tools: USER needs to change according where the tools are installed.##

#my $mutect="/gscuser/scao/tools/mutect-1.1.7.jar";

my $h38_REF_bai=$h38_REF.".fai";
my $f_gtf= "/storage1/fs1/songcao/Active/Database/hg38_database/sid.ont/Homo_sapiens.GRCh38.95.gtf";

my $first_line=`head -n 1 $h38_REF`;

if($first_line=~/^\>chr/) { $chr_status=1; }

opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

my $dir_input=$run_dir."/inputs"; 
my $dir_output=$run_dir."/outputs";

if(!-d $dir_input) { `mkdir $dir_input`; }
if(!-d $dir_output) { `mkdir $dir_output`; }

# check to make sure the input directory has correct structure
#&check_input_dir($run_dir);
# start data processsing

#`bgadd -L 70 $compute_username/$group_name`;

print "start to run","\n";

if (($step_number==0 || $step_number==1)) {
    #begin to process each sample
    for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
        $sample_name = $sample_dir_list[$i];
        if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
            $sample_full_path = $run_dir."/".$sample_name;
            if (-d $sample_full_path) { # is a full path directory containing a sample
                print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
                $current_job_file="";
    
                if($step_number==0)
                {  
		        &bsub_preprocess_bam(); # split bam
                }

                if($step_number==1 && !($sample_name=~/inputs/ || $sample_name=~/outputs/))
                {  
		        &bsub_split_bam(); # split bam
	            }
            }
        }
    }
}       

if($step_number == 2) {
    &bsub_espresso_1(1); #run espresso step 1
                }
if ($step_number == 3) {
    &bsub_espresso_2(1);
                }
if ($step_number == 4) {
	&bsub_espresso_3(1);
                }
if ($step_number == 5) {
	&bsub_merge_espresso(1);
                }
if ($step_number == 6) {
	&bsub_sqanti(1);
                }
if ($step_number == 7) {
    &bsub_re_espresso(1);
                }
if ($step_number == 8) {
	&bsub_merge_espresso_filter(1);
                }
if ($step_number == 9) {
	bsub_filter_isoform(1);
                }
if ($step_number == 10)
                {
    &bsub_re2_espresso(1);  
                }
if ($step_number == 11)
                {
    &bsub_merge_fsm(1);  
                }

if ($step_number == 12)
                {
    &bsub_run_rmats(1);  
                }

if ($step_number == 13)
                {
    &bsub_summary_report(1);  
                }


#exit;

sub bsub_preprocess_bam{

}

sub bsub_split_bam{

    my $splitDIR=$sample_full_path."/splitbam";

    if(-d $splitDIR)
    {
        `rm -rf $splitDIR`;
    }

    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

    foreach my $chr (@chrlist)
    {

    	my $chr1=$chr;
   
        if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
        }

        my $input_tsv=$dir_input."/".$chr1.".bam.tsv";

        $current_job_file = "j1_splitbam_".$sample_name.".".$chr1.".sh"; 
	    my $IN_bam = $sample_full_path."/".$sample_name.".bam";

        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
        `rm $lsf_out`;
        `rm $lsf_err`;

        open(SPLITBAM, ">$job_files_dir/$current_job_file") or die $!;
        print SPLITBAM "#!/bin/bash\n";
        print SPLITBAM "BAM=".$sample_full_path."/".$sample_name.".bam\n";
        print SPLITBAM "RUNDIR=".$sample_full_path."/splitbam\n";
        print SPLITBAM "if [ ! -d \${RUNDIR} ]\n";
        print SPLITBAM "then\n";
        print SPLITBAM "mkdir \${RUNDIR}\n";
        print SPLITBAM "fi\n";
        ## extract reads from chr1 to chr22, chrX, and chrY
        ## missing chr10, chr20
        print SPLITBAM "samtools view -F 0x900 -h \${BAM} $chr1:1-999999999 | samtools view -Sb > \${RUNDIR}/$chr1.bam","\n";   
        print SPLITBAM "echo \"\${RUNDIR}/$chr1.bam $sample_name\" | awk \'OFS=\"\\t\"\{print \$1,\$2\}\' > \${RUNDIR}/$chr1.bam.tsv","\n";
        print SPLITBAM "echo \"\${RUNDIR}/$chr1.bam $sample_name\" | awk \'OFS=\"\\t\"\{print \$1,\$2\}\' >> $input_tsv","\n";
        close SPLITBAM;
        my $sh_file=$job_files_dir."/".$current_job_file;
        $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30G] rusage[mem=30G]\" -M 30G -a \'docker(chrisamiller/genomic-analysis:0.2)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
        print $bsub_com;
        system ($bsub_com);
    }
 
}

## run espresso step 1
sub bsub_espresso_1{

    my $dir_espresso=$dir_output."/espresso"; 

    if(!-d $dir_espresso)
    {
        `mkdir $dir_espresso`;
    }

   # my @chrlist=("10","20");

    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

    foreach my $chr (@chrlist)
    {
    	my $chr1=$chr;
   
        if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
        }
    
        $current_job_file = "j2_espresso_1_inputs".".".$chr1.".sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_espresso."/".$chr1; 
        
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }

        open(EXPRESSO1, ">$job_files_dir/$current_job_file") or die $!;
        print EXPRESSO1 "#!/bin/bash\n";
        print EXPRESSO1 "RUNDIR=".$dir_espresso."/"."$chr1\n";
        #print EXPRESSO1 "RUNDIR=".$dir_espresso."\n";
        print EXPRESSO1 "if [ ! -d \${RUNDIR} ]\n";
        print EXPRESSO1 "then\n";
        print EXPRESSO1 "mkdir \${RUNDIR}\n";
        print EXPRESSO1 "fi\n";
        print EXPRESSO1 "input_tsv=".$dir_input."/".$chr1.".bam.tsv\n";
        #print EXPRESSO1 "source activate env","\n";
        #print EXPRESSO1 "perl /bin/espresso/src/ESPRESSO_S.pl -L \${input_tsv} -F $h38_REF -A $f_gtf -O \${RUNDIR} -T 4","\n";
        print EXPRESSO1 "source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L \${input_tsv} -F $h38_REF -A $f_gtf -O \${RUNDIR} -T 4","\n"; 
        close EXPRESSO1;
        my $sh_file=$job_files_dir."/".$current_job_file;
       # $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 8 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L ${input_tsv} -F $h38_REF -A $f_gtf -O ${outdir}/${base} -T 4\"","\n";
       $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>100G] rusage[mem=100G]\" -M 100G -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
       #$bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L $input_tsv -F $h38_REF -A $f_gtf -O $RUNDIR -T 4\"","\n"; 
       print $bsub_com;
       system ($bsub_com);
    }        

}


## run espresso step 2
sub bsub_espresso_2{

    my $dir_espresso=$dir_output."/espresso"; 

    if(!-d $dir_espresso)
    {
        `mkdir $dir_espresso`;
    }

    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

    #my @chrlist=("20","21","22"); 
    #my @chrlist=("1");
  #  my @chrlist=("10","20");

    foreach my $chr (@chrlist)
    {
    	my $chr1=$chr;
    
        if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
        }

        my $f_espresso_out=$dir_output."/espresso/".$chr1."/bam_N2_R0_updated.gtf";
        if(!(-s $f_espresso_out))
        {
         print $f_espresso_out," is not existing\n"; 
        my $input_tsv=$dir_input."/".$chr1.".bam.tsv";
        my $input_update_tsv=$dir_output."/espresso/".$chr1."/".$chr1.".bam.tsv.updated"; 
        ## loop each sample
        foreach my $l (`cat $input_update_tsv`) 
        {
        my $ltr=$l; 
        chomp($ltr);
        my @t=split("\t",$ltr);
        my $sample_id=$t[2];
        $current_job_file = "j3_espresso_2_".$sample_id.".$chr1."."sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_espresso."/".$chr1; 
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }    
        open(EXPRESSO2, ">$job_files_dir/$current_job_file") or die $!;
        print EXPRESSO2 "#!/bin/bash\n";
        print EXPRESSO2 "RUNDIR=".$dir_espresso."/"."$chr1\n";
        print EXPRESSO2 "source activate env && perl /bin/espresso/src/ESPRESSO_C.pl -I \${RUNDIR} -F $h38_REF -X $sample_id -T 4","\n";
        #ource activate env && perl /bin/espresso/src/ESPRESSO_C.pl -I ${outdir}/${base}/ -F ${reference_fasta} -X ${sample} -T $threads
        #print EXPRESSO2 "perl /bin/espresso/src/ESPRESSO_S.pl -L \${input_tsv} -F $h38_REF -A $f_gtf -O \${RUNDIR} -T 4","\n";
        close EXPRESSO2;
        my $sh_file=$job_files_dir."/".$current_job_file;
       # $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 8 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L ${input_tsv} -F $h38_REF -A $f_gtf -O ${outdir}/${base} -T 4\"","\n";
       $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>100G] rusage[mem=100G]\" -M 100G -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
       #$bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L $input_tsv -F $h38_REF -A $f_gtf -O $RUNDIR -T 4\"","\n"; 
       print $bsub_com;
       system ($bsub_com);
        }
    }        
  }

}
## run espresso step 3
sub bsub_espresso_3{

    my $dir_espresso=$dir_output."/espresso"; 

    if(!-d $dir_espresso)
    {
        `mkdir $dir_espresso`;
    }

     my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");
     #my @chrlist=("20","21","22");   
    #my @chrlist=("1");

    foreach my $chr (@chrlist)
    {
    	my $chr1=$chr;
   
        if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
        }

        my $f_espresso_out=$dir_output."/espresso/".$chr1."/bam_N2_R0_updated.gtf";
        if(!(-s $f_espresso_out))
        {    
       # $current_job_file = "j4_espresso_3_inputs.".".$chr1.".sh"; 
        $current_job_file = "j4_espresso_3_".$chr1.".sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_espresso."/".$chr1; 
        
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }

        my $input_tsv=$dir_input."/".$chr1.".bam.tsv";
        my $input_update_tsv=$dir_output."/espresso/".$chr1."/".$chr1.".bam.tsv.updated";

        open(EXPRESSO3, ">$job_files_dir/$current_job_file") or die $!;
        print EXPRESSO3 "#!/bin/bash\n";
        print EXPRESSO3 "RUNDIR=".$dir_espresso."/"."$chr1\n";
        print EXPRESSO3 "ISOFORM=".$dir_espresso."/".$chr1."/".$chr1."_compatible_isoform.tsv","\n";
        print EXPRESSO3 "source activate env && perl /bin/espresso/src/ESPRESSO_Q.pl -L $input_update_tsv -A $f_gtf -V \${ISOFORM} -T 4","\n";
        close EXPRESSO3;
        my $sh_file=$job_files_dir."/".$current_job_file;
        $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>100G] rusage[mem=100G]\" -M 100G -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
        print $bsub_com;
       system ($bsub_com);
    }   
    }     

}

## merge outputs from all chromosomes into abundance/gtf files

sub bsub_merge_espresso{

    my $dir_espresso=$dir_output."/espresso"; 

    if(!-d $dir_espresso)
    {
        `mkdir $dir_espresso`;
    }

   
    $current_job_file = "j5_merge_espresso".".sh"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

    `rm $lsf_out`;
    `rm $lsf_err`;
  
    open(MERGE, ">$job_files_dir/$current_job_file") or die $!;
    print MERGE "#!/bin/bash\n";
    print MERGE "RUNDIR=".$dir_espresso."\n";
    print MERGE "     ".$run_perl_script_path."merge_espresso.pl \${RUNDIR} $chr_status\n";
    close MERGE;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30G] rusage[mem=30G]\" -M 30G -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    print $bsub_com;
    system ($bsub_com);         

}

# starsjs=$(cat $star_sjouts | perl -pe 's/\n/,/g' | perl -pe 's/,$//g') #create comma-sep input
# mkdir $scratch/tmp
# count=1;
# cat $star_sjouts | while read i;do 
#     cp $i $scratch/tmp/$count.SJ.out.tab;
#     count=$((count+1));
# done
# star_sjouts=$scratch/tmp
# LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-timley -J sqanti -oo $scratch/logs/sqanti_qc.log -n 8 -R"select[mem>32G] rusage[mem=32G] span[hosts=1]" -M 32G -q siteman -a "docker(chrisamiller/sqanti3:v5.2.1)" /bin/bash -c "source activate SQANTI3.env && /app/SQANTI3-5.2.1/sqanti3_qc.py -t 8 -d $scratch/outputs/sqanti -c $scratch/tmp --SR_bam $star_bams $scratch/outputs/merged_N2_R0_updated.gtf $gtf $reffasta"
# /storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n sqanti -s 60 -d 30

# #run sqanti filtering using the machine-learning approach
# LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-timley -J sqanti -oo $scratch/logs/sqanti_filter.log -n 8 -R"select[mem>32G] rusage[mem=32G] span[hosts=1]" -M 32G -q siteman -a "docker(chrisamiller/sqanti3:v5.2.1)" /bin/bash -c "source activate SQANTI3.env && /app/SQANTI3-5.2.1/sqanti3_filter.py ml $scratch/outputs/sqanti/merged_N2_R0_updated_classification.txt --gtf $scratch/outputs/sqanti/merged_N2_R0_updated_corrected.gtf -d $scratch/outputs/sqanti/filter_ml_default"

sub bsub_sqanti{
     
        $current_job_file = "j6_sqanti".".sh"; 

        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
        my $dir_espresso=$dir_output."/espresso"; 
        my $f_sj=$run_dir."/sj_input_shortread.tsv"; 
        my $f_star_bam=$run_dir."/bam_input_shortread.tsv"; 
        my $dirsqanti=$dir_output."/sqanti"; 
        my $dirml=$dir_output."/sqanti/filter_ml_default"; 
        my $dirtmp=$dir_output."/tmp";

        my $count = 0;

        open(my $fh, "<", $f_sj) or die "Cannot open $f_sj: $!";
        while (my $f = <$fh>) {
        chomp($f);
        my $f_out = "$dirtmp/$count.SJ.out.tab";
        if (-e $f) {
            system("cp", $f, $f_out);
        }
        $count++;
        }
        close $fh;  

        if(-d $dirsqanti)
        {
            `rm -rf $dirsqanti`; 
        } 
        `mkdir $dirsqanti`; 
        
        open(SQANTI, ">$job_files_dir/$current_job_file") or die $!;
        print SQANTI "#!/bin/bash\n";
        print SQANTI "RUNDIR=".$dir_espresso."\n";
        print SQANTI "mergedgtf=".$dir_espresso."/merged_N2_R0_updated.gtf","\n";
        print SQANTI "correctedgtf=".$dirsqanti."/merged_N2_R0_updated_corrected.gtf","\n";
        print SQANTI "mlgtf=".$dirsqanti."/merged_N2_R0_updated_classification.txt","\n"; 
        print SQANTI "source activate SQANTI3.env &&  /app/SQANTI3-5.1.2/sqanti3_qc.py -t 4 -d $dirsqanti -c $dirtmp --SR_bam $f_star_bam \${mergedgtf} $f_gtf $h38_REF","\n";
        print SQANTI "source activate SQANTI3.env && /app/SQANTI3-5.1.2/sqanti3_filter.py ml \${mlgtf} --gtf \${correctedgtf} -d $dirml","\n"; 
        close SQANTI;
        my $sh_file=$job_files_dir."/".$current_job_file;
       # $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 8 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L ${input_tsv} -F $h38_REF -A $f_gtf -O ${outdir}/${base} -T 4\"","\n";
       $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>100G] rusage[mem=100G]\" -M 100G -a \'docker(chrisamiller/sqanti3:v5.1.2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
       #$bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L $input_tsv -F $h38_REF -A $f_gtf -O $RUNDIR -T 4\"","\n"; 
       print $bsub_com;
       system ($bsub_com);

}

## re-run espresso 
sub bsub_re_espresso{
    my $dir_espresso=$dir_output."/espresso"; 
    if(!-d $dir_espresso)
    {
        `mkdir $dir_espresso`;
    }
    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");
    #my @chrlist=("20","21","22");   
    # my @chrlist=("1");
    foreach my $chr (@chrlist)
    {
    	my $chr1=$chr;
   
        if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
        }

        $current_job_file = "j7_espresso_r3_".$chr1.".sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_espresso."/".$chr1; 
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }

        my $input_tsv=$dir_input."/".$chr1.".bam.tsv";
        my $input_update_tsv=$dir_output."/espresso/".$chr1."/".$chr1.".bam.tsv.updated";
        open(ESPR, ">$job_files_dir/$current_job_file") or die $!;
        print ESPR "#!/bin/bash\n";
        print ESPR "RUNDIRfilter=".$dir_espresso."/".$chr1."/filtered\n";
        print ESPR "if [ ! -d \${RUNDIRfilter} ]\n";
        print ESPR "then\n";
        print ESPR "mkdir \${RUNDIRfilter}\n";
        print ESPR "fi\n";
        print ESPR "filtergtfin=".$dir_output."/sqanti/filter_ml_default/merged_N2_R0_updated.filtered.gtf\n";
        print ESPR "filtergtfout=".$dir_espresso."/".$chr1."/".$chr1."_N2_R0_updated.filtered.gtf\n";
        print ESPR "ISOFORM=".$dir_espresso."/".$chr1."/filtered/".$chr1."_compatible_isoform.tsv","\n";
        print ESPR "grep \"\^$chr1\[\[\:space\:\]\]\" \${filtergtfin} >  \${filtergtfout}","\n";
        print ESPR "source activate env && perl /bin/espresso/src/ESPRESSO_Q.pl -L $input_update_tsv -O \${RUNDIRfilter} -A \${filtergtfout} --read_ratio_cutoff 2 -V \${ISOFORM} -T 4","\n";
        close ESPR;
        my $sh_file=$job_files_dir."/".$current_job_file;
        $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>100G] rusage[mem=100G]\" -M 100G -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
        print $bsub_com;
        system ($bsub_com);
    }        

}

## merge outputs from all chromosomes into abundance/gtf files

sub bsub_merge_espresso_filter{

    my $dir_espresso=$dir_output."/espresso"; 

    if(!-d $dir_espresso)
    {
        `mkdir $dir_espresso`;
    }

   
    $current_job_file = "j8_merge_espresso_filter".".sh"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

    `rm $lsf_out`;
    `rm $lsf_err`;
  
    open(MERGE, ">$job_files_dir/$current_job_file") or die $!;
    print MERGE "#!/bin/bash\n";
    print MERGE "RUNDIR=".$dir_espresso."\n";
    print MERGE "     ".$run_perl_script_path."merge_espresso_filter.pl \${RUNDIR} $chr_status\n";
    close MERGE;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>100G] rusage[mem=100G]\" -M 100G -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    print $bsub_com;
    system ($bsub_com);         

}

#running espresso on the filtered gtf told us whether every read was FSM, ISM, etc.  Take that information and filter
#the espresso intermediate files so that they only contain FSM reads. Basically, copy the outputs directories and refill them
#with the same data with non-FSM reads removed.

sub bsub_filter_isoform{

    my $dir_espresso=$dir_output."/espresso"; 
    
    if(!-d $dir_espresso)
    {
        `mkdir $dir_espresso`;
    }

    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

    foreach my $chr (@chrlist)
    {
    	my $chr1=$chr;
   
        if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
        }
    
        $current_job_file = "j9_filter_isoform_".$chr1.".sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_espresso."/".$chr1; 
        
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }

        my $input_tsv=$dir_input."/".$chr1.".bam.tsv";
        my $input_update_tsv=$dir_output."/espresso/".$chr1."/".$chr1.".bam.tsv.updated";

        open(FILTER, ">$job_files_dir/$current_job_file") or die $!;
        print FILTER "#!/bin/bash\n";
        print FILTER "RUNDIR=".$dir_espresso."/"."$chr1\n";
        print FILTER "ISOFORM=".$dir_espresso."/".$chr1."/filtered/".$chr1."_compatible_isoform.tsv","\n";
        print FILTER "     ".$run_perl_script_path."remove_nonfsm_data.pl \${ISOFORM} \${RUNDIR} $input_update_tsv","\n";
        close FILTER;
        my $sh_file=$job_files_dir."/".$current_job_file;
        $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30G] rusage[mem=30G]\" -M 30G -a \'docker(chrisamiller/docker-genomic-analysis:latest)\' -o $lsf_out -e $lsf_err sh $sh_file\n"; 
        print $bsub_com;
        system ($bsub_com); 

    }

}


#okay, rerun espresso_q (again!) per chromosome, using these filtered files as inputs

## re-run espresso 
sub bsub_re2_espresso{
    
    my $dir_espresso=$dir_output."/espresso"; 

    if(!-d $dir_espresso)
    {
        `mkdir $dir_espresso`;
    }
    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");
   # my @chrlist=("1");
    foreach my $chr (@chrlist)
    {
    	my $chr1=$chr;
   
        if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
        }

        $current_job_file = "j10_espresso_r3_fsm_".".$chr1.sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_espresso."/".$chr1; 
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }
        my $input_tsv=$dir_input."/".$chr1.".bam.tsv";
        my $input_update_tsv=$dir_output."/espresso/".$chr1."/fsm/".$chr1.".bam.tsv.updated";
        open(ESPR2, ">$job_files_dir/$current_job_file") or die $!;
        print ESPR2 "#!/bin/bash\n";
        print ESPR2 "RUNDIRfsm=".$dir_espresso."/".$chr1."/fsm\n";
        print ESPR2 "filtergtfout=".$dir_espresso."/".$chr1."/".$chr1."_N2_R0_updated.filtered.gtf\n";
        print ESPR2 "ISOFORM=".$dir_espresso."/".$chr1."/fsm/".$chr1."_compatible_isoform.tsv","\n";
        print ESPR2 "source activate env && perl /bin/espresso/src/ESPRESSO_Q.pl -L $input_update_tsv -O \${RUNDIRfsm} -A \${filtergtfout} --read_ratio_cutoff 2 -V \${ISOFORM} -T 4","\n";
        close ESPR2;
        my $sh_file=$job_files_dir."/".$current_job_file;
        $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30G] rusage[mem=30G]\" -M 30G -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
        print $bsub_com;
        system ($bsub_com);
    }        
}

sub bsub_merge_fsm{

    my $dir_espresso=$dir_output."/espresso"; 

    if(!-d $dir_espresso)
    {
        `mkdir $dir_espresso`;
    }

   
    $current_job_file = "j11_merge_fsm".".sh"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

    `rm $lsf_out`;
    `rm $lsf_err`;
  
    open(MERGEFSM, ">$job_files_dir/$current_job_file") or die $!;
    print MERGEFSM "#!/bin/bash\n";
    print MERGEFSM "RUNDIR=".$dir_espresso."\n";
    print MERGEFSM "     ".$run_perl_script_path."merge_fsm.pl \${RUNDIR} $chr_status\n";
    close MERGEFSM;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30G] rusage[mem=30G]\" -M 30G -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    print $bsub_com;
    system ($bsub_com);         

}

sub bsub_run_rmats{

    my $outdir_rmats=$dir_output."/rmats-long-fsm"; 
    my $dir_espresso=$dir_output."/espresso"; 

    if(-d $outdir_rmats)
    {
        `rm -rf $outdir_rmats`;
    }

   `mkdir $outdir_rmats`;

    $current_job_file = "j12_rmats_long".".sh"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

    `rm $lsf_out`;
    `rm $lsf_err`;
    my $f_config=$dir_input."/configs.tsv"; 
    open(OUT,">$f_config"); 

    foreach my $s (`ls $run_dir`)
    {
        my $str=$s; chomp($str);
        if($run_dir=~/S34F_DMSO-S34F_SMG1/ || $run_dir=~/S34F_DMSO-WT_SMG1/ || $run_dir=~/S34F_SMG1-WT_DMSO/ || $run_dir=~/WT_DMSO-WT_SMG1/ || $run_dir=~/K562-12/)
        {
        if($str=~/DMSO/) { print OUT "g1","\t",$str,"\t",$run_dir."/".$str."/".$str.".bam","\n"; }
        if($str=~/SMG1/) { print OUT "g2","\t",$str,"\t",$run_dir."/".$str."/".$str.".bam","\n"; }
        }
        if($run_dir=~/S34F_SMG1-WT_SMG1/ || $run_dir=~/S34F_DMSO-WT_DMSO/)
        {
        if($str=~/WT/) { print OUT "g1","\t",$str,"\t",$run_dir."/".$str."/".$str.".bam","\n"; }
        if($str=~/S34F/) { print OUT "g2","\t",$str,"\t",$run_dir."/".$str."/".$str.".bam","\n"; }
        } 
    }   
    close OUT; 

#  grep "^g1" $config | cut -f 2 | perl -pe 's/\n/,/g' | perl -pe 's/,$//g' >${outdir}/group1.txt
# grep "^g2" $config | cut -f 2 | perl -pe 's/\n/,/g' | perl -pe 's/,$//g' >${outdir}/group2.txt

#  LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-timley -q general -J rmats -n $threads -R"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" -M ${mem} -oo ${logdir}/rmats-long_step7.out -eo ${logdir}/rmats-long_step7.err -a "docker(sridnona/rmats_long:v3)" /bin/bash -c "source activate /docker_data/rMATS-long/conda_env && cd /docker_data/rMATS-long/scripts/ && python rmats_long.py \
# --abundance ${espresso_abundance} --updated-gtf ${espresso_gtf} \

    open(RMATS, ">$job_files_dir/$current_job_file") or die $!;
    print RMATS "#!/bin/bash\n";
    print RMATS "espresso_abundance=".$dir_espresso."/merged_fsm_N2_R2_abundance.esp\n";
    print RMATS "espresso_gtf=".$dir_espresso."/merged_fsm_N2_R2_updated.gtf","\n";
    print RMATS "grep \"\^g1\" $f_config | cut -f 2 | perl -pe \'s\/\\n\/\,\/g\' | perl -pe \'s\/\,\$\/\/g\' > $outdir_rmats/group1.txt","\n";
    print RMATS "grep \"\^g2\" $f_config | cut -f 2 | perl -pe \'s\/\\n\/\,\/g\' | perl -pe \'s\/\,\$\/\/g\' > $outdir_rmats/group2.txt","\n";
    print RMATS "source activate /docker_data/rMATS-long/conda_env && cd /docker_data/rMATS-long/scripts/ && python rmats_long.py --abundance \${espresso_abundance} --updated-gtf \${espresso_gtf} --gencode-gtf $f_gtf --group-1 $outdir_rmats/group1.txt --group-2 $outdir_rmats/group2.txt --group-1-name group1 --group-2-name group2 --out-dir $outdir_rmats --num-threads 16 --delta-proportion 0.1 --adj-pvalue 0.1 --plot-file-type .pdf","\n";
    close RMATS;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 16 -R \"select[mem>30G] rusage[mem=30G]\" -M 30G -a \'docker(sridnona/rmats_long:v3)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
    print $bsub_com;
    system ($bsub_com);    

}

sub bsub_summary_report{

   my $outdir_sum =$dir_output."/summary"; 

    if(!-d $outdir_sum)
    {
        `mkdir $outdir_sum`;
    }

    $current_job_file = "j13_summary_report".".sh"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

    `rm $lsf_out`;
    `rm $lsf_err`;
    open(SUMMARY, ">$job_files_dir/$current_job_file") or die $!;
    print SUMMARY "#!/bin/bash\n";
    print SUMMARY "RUNDIR=".$dir_output."\n";
    print SUMMARY "     ".$run_perl_script_path."generate_summary.pl \${RUNDIR} $f_gtf\n";
    print SUMMARY "     ".$run_perl_script_path."generate_summary.nofilter.pl \${RUNDIR} $f_gtf\n";
    close SUMMARY;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30G] rusage[mem=30G]\" -M 30G -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    print $bsub_com;
    system ($bsub_com);      
}
