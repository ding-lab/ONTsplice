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
$yellow [6]  run sqanti
$yellow [7]  re-run step expresso 3 by using new filtered gtf file from sqanti
$yellow [8]  merge filtered expresso results
$cyan [9]  run fsm 
$cyan [10]  re-run step expresso 3 based fsm results
$cyan [11]  merge per chr results from 10 expresso results
$red [12] run rmats
$normal [13] merge bams and output novel transcript

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
my $somatic_path = "/storage1/fs1/songcao/Active/Git/SComatic";
my $run_scomatic_path = "python $somatic_path"."/";

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

# check to make sure the input directory has correct structure
#&check_input_dir($run_dir);
# start data processsing

#`bgadd -L 70 $compute_username/$group_name`;

print "start to run","\n";

if (($step_number < 14 && $step_number>=0)) {
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

                if($step_number==1)
                {  
		        &bsub_split_bam(); # split bam
	            }elsif ($step_number == 2) {
                     &bsub_espresso_1(1); #run espresso step 1
                }elsif ($step_number == 3) {
                    &bsub_espresso_2(1);
                }elsif ($step_number == 4) {
		          &bsub_espresso_3(1);
                }elsif ($step_number == 5) {
		          &bsub_merge_espresso(1);
                }elsif ($step_number == 6) {
		          &bsub_sqanti(1);
                }elsif ($step_number == 7) {
                &bsub_re_espresso(1);
                }elsif ($step_number == 8) {
		          &bsub_merge_espresso_filter(1);
                }
                elsif ($step_number == 9) {
		          &bsub_filter_isoform(1);
                }elsif ($step_number == 10)
                {
                  &bsub_re2_espresso(1);  
                }elsif ($step_number == 11)
                {
                  &bsub_merge_fsm(1);  
                }
           }
        }
    }
}

if($step_number==14)
    {

    print $yellow, "Submitting jobs for generating the report for the run ....",$normal, "\n";
    $hold_job_file=$current_job_file; 
    $current_job_file = "j14_Run_report_".$working_name.".sh"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    `rm $current_job_file`;
    my $working_name= (split(/\//,$run_dir))[-1];
    my $f_maf=$run_dir."/".$working_name.".withmutect.maf";
    my $f_maf_rc=$f_maf.".rc";
    my $f_maf_rc_caller=$f_maf_rc.".caller";

    open(REPRUN, ">$job_files_dir/$current_job_file") or die $!;
    print REPRUN "#!/bin/bash\n";
    print REPRUN "		".$run_perl_script_path."generate_final_report.pl ".$run_dir."\n";
    print REPRUN "      ".$run_perl_script_path."add_rc.pl ".$run_dir." ".$f_maf." ".$f_maf_rc."\n";
    print REPRUN "      ".$run_perl_script_path."add_caller.pl ".$run_dir." ".$f_maf_rc." ".$f_maf_rc_caller."\n";
    close REPRUN;

    my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err bash $sh_file\n";

    print $bsub_com;
    system ($bsub_com);

}

exit;

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
        close SPLITBAM;
        my $sh_file=$job_files_dir."/".$current_job_file;
        $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(chrisamiller/genomic-analysis:0.2)\' -o $lsf_out -e $lsf_err bash $sh_file\n";
        print $bsub_com;
        system ($bsub_com);
    }
 
}

## run espresso step 1
sub bsub_espresso_1{

    my $dir_expresso=$sample_full_path."/expresso1"; 

    if(!-d $dir_expresso)
    {
        `mkdir $dir_expresso`;
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
    
        $current_job_file = "j2_expresso_1_".$sample_name.".".$chr1.".sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_expresso."/".$chr1; 
        
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }

        my $input_tsv=$sample_full_path."/splitbam/".$chr1.".bam.tsv";

        open(EXPRESSO1, ">$job_files_dir/$current_job_file") or die $!;
        print EXPRESSO1 "#!/bin/bash\n";
        print EXPRESSO1 "RUNDIR=".$dir_expresso."/"."$chr1\n";
        print EXPRESSO1 "if [ ! -d \${RUNDIR} ]\n";
        print EXPRESSO1 "then\n";
        print EXPRESSO1 "mkdir \${RUNDIR}\n";
        print EXPRESSO1 "fi\n";
        print EXPRESSO1 "input_tsv=".$sample_full_path."/splitbam/".$chr1.".bam.tsv\n";
        #print EXPRESSO1 "source activate env","\n";
        #print EXPRESSO1 "perl /bin/espresso/src/ESPRESSO_S.pl -L \${input_tsv} -F $h38_REF -A $f_gtf -O \${RUNDIR} -T 4","\n";
        print EXPRESSO1 "source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L \${input_tsv} -F $h38_REF -A $f_gtf -O \${RUNDIR} -T 2","\n"; 
        close EXPRESSO1;
        my $sh_file=$job_files_dir."/".$current_job_file;
       # $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 8 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L ${input_tsv} -F $h38_REF -A $f_gtf -O ${outdir}/${base} -T 4\"","\n";
       $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 2 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
       #$bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L $input_tsv -F $h38_REF -A $f_gtf -O $RUNDIR -T 4\"","\n"; 
       print $bsub_com;
       system ($bsub_com);
    }        

}




## run espresso step 2
sub bsub_espresso_2{

    my $dir_expresso=$sample_full_path."/expresso1"; 

    if(!-d $dir_expresso)
    {
        `mkdir $dir_expresso`;
    }

    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

  #  my @chrlist=("10","20");

    foreach my $chr (@chrlist)
    {
    	my $chr1=$chr;
   
        if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
        }
    
        $current_job_file = "j3_expresso_2_".$sample_name.".".$chr1.".sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_expresso."/".$chr1; 
        
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }

        my $input_tsv=$sample_full_path."/splitbam/".$chr1.".bam.tsv";
        my $input_update_tsv=$sample_full_path."/expresso1/".$chr1."/".$chr1.".bam.tsv.updated";

        open(EXPRESSO2, ">$job_files_dir/$current_job_file") or die $!;
        print EXPRESSO2 "#!/bin/bash\n";
        print EXPRESSO2 "RUNDIR=".$dir_expresso."/"."$chr1\n";
        print EXPRESSO2 "sample_id=\$(cut -f 3 $input_update_tsv)","\n";
        print EXPRESSO2 "source activate env && perl /bin/espresso/src/ESPRESSO_C.pl -I \${RUNDIR} -F $h38_REF -X \${sample_id} -T 2","\n";
        #ource activate env && perl /bin/espresso/src/ESPRESSO_C.pl -I ${outdir}/${base}/ -F ${reference_fasta} -X ${sample} -T $threads
        #print EXPRESSO2 "perl /bin/espresso/src/ESPRESSO_S.pl -L \${input_tsv} -F $h38_REF -A $f_gtf -O \${RUNDIR} -T 4","\n";
        close EXPRESSO2;
        my $sh_file=$job_files_dir."/".$current_job_file;
       # $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 8 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L ${input_tsv} -F $h38_REF -A $f_gtf -O ${outdir}/${base} -T 4\"","\n";
       $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 2 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
       #$bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L $input_tsv -F $h38_REF -A $f_gtf -O $RUNDIR -T 4\"","\n"; 
       print $bsub_com;
       system ($bsub_com);
    }        

}

## run espresso step 3
sub bsub_espresso_3{

    my $dir_expresso=$sample_full_path."/expresso1"; 

    if(!-d $dir_expresso)
    {
        `mkdir $dir_expresso`;
    }

     my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");
    
    #my @chrlist=("1");

    foreach my $chr (@chrlist)
    {
    	my $chr1=$chr;
   
        if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
        }
    
        $current_job_file = "j4_expresso_3_".$sample_name.".".$chr1.".sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_expresso."/".$chr1; 
        
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }

        my $input_tsv=$sample_full_path."/splitbam/".$chr1.".bam.tsv";
        my $input_update_tsv=$sample_full_path."/expresso1/".$chr1."/".$chr1.".bam.tsv.updated";

        open(EXPRESSO3, ">$job_files_dir/$current_job_file") or die $!;
        print EXPRESSO3 "#!/bin/bash\n";
        print EXPRESSO3 "RUNDIR=".$dir_expresso."/"."$chr1\n";
        print EXPRESSO3 "ISOFORM=".$dir_expresso."/".$chr1."/".$chr1."_compatible_isoform.tsv","\n";
        print EXPRESSO3 "source activate env && perl /bin/espresso/src/ESPRESSO_Q.pl -L $input_update_tsv -A $f_gtf -V \${ISOFORM} -T 4","\n";
        close EXPRESSO3;
        my $sh_file=$job_files_dir."/".$current_job_file;
        $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
        print $bsub_com;
       system ($bsub_com);
    }        

}

## merge outputs from all chromosomes into abundance/gtf files

sub bsub_merge_espresso{

    my $dir_expresso=$sample_full_path."/expresso1"; 

    if(!-d $dir_expresso)
    {
        `mkdir $dir_expresso`;
    }

   
    $current_job_file = "j5_merge_expresso_".$sample_name.".".".sh"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

    `rm $lsf_out`;
    `rm $lsf_err`;
  
    open(MERGE, ">$job_files_dir/$current_job_file") or die $!;
    print MERGE "#!/bin/bash\n";
    print MERGE "RUNDIR=".$dir_expresso."\n";
    print MERGE "     ".$run_perl_script_path."merge_espresso.pl \${RUNDIR} $chr_status\n";
    close MERGE;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    print $bsub_com;
    system ($bsub_com);         

}

#mkdir $scratch/outputs/sqanti
#LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-timley -J sqanti -oo $scratch/logs/sqanti_qc.log -n 8 -R"select[mem>32G] rusage[mem=32G] span[hosts=1]" -M 32G -q siteman -a "docker(chrisamiller/sqanti3:v5.1.2)" /bin/bash -c "source activate SQANTI3.env && /app/SQANTI3-5.1.2/sqanti3_qc.py -t 8 -d $scratch/outputs/sqanti $scratch/outputs/merged_N2_R0_updated.gtf $gtf $reffasta"
#/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n sqanti -s 60

#run sqanti filtering using the machine-learning approach
#LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -G compute-timley -J sqanti -oo $scratch/logs/sqanti_filter.log -n 8 -R"select[mem>32G] rusage[mem=32G] span[hosts=1]" -M 32G -q siteman -a "docker(chrisamiller/sqanti3:v5.1.2)" /bin/bash -c "source activate SQANTI3.env && /app/SQANTI3-5.1.2/sqanti3_filter.py ML $scratch/outputs/sqanti/merged_N2_R0_updated_classification.txt --gtf $scratch/outputs/sqanti/merged_N2_R0_updated_corrected.gtf -d $scratch/outputs/sqanti/filter_ml_default"
#/storage1/fs1/timley/Active/aml_ppg/src/utilities/bwait -n sqanti -s 60

sub bsub_sqanti{
        
        $current_job_file = "j6_sqanti_".$sample_name.".".".sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
        my $dir_expresso=$sample_full_path."/expresso1"; 
        my $dirsqanti=$sample_full_path."/sqanti"; 
        my $dirml=$sample_full_path."/sqanti/filter_ml_default"; 

        if(!-d $dirsqanti)
        {
          `mkdir $dirsqanti`; 
        }

        open(SQANTI, ">$job_files_dir/$current_job_file") or die $!;
        print SQANTI "#!/bin/bash\n";
        print SQANTI "RUNDIR=".$dir_expresso."\n";
        print SQANTI "mergedgtf=".$dir_expresso."/merged_N2_R0_updated.gtf","\n";
        print SQANTI "correctedgtf=".$dirsqanti."/merged_N2_R0_updated_corrected.gtf","\n";
        print SQANTI "mlgtf=".$dirsqanti."/merged_N2_R0_updated_classification.txt","\n";

        print SQANTI "source activate SQANTI3.env &&  /app/SQANTI3-5.1.2/sqanti3_qc.py -t 4 -d $dirsqanti \${mergedgtf} $f_gtf $h38_REF","\n";
        print SQANTI "source activate SQANTI3.env && /app/SQANTI3-5.1.2/sqanti3_filter.py ML \${mlgtf} --gtf \${correctedgtf} -d $dirml","\n"; 
        close SQANTI;
        my $sh_file=$job_files_dir."/".$current_job_file;
       # $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 8 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L ${input_tsv} -F $h38_REF -A $f_gtf -O ${outdir}/${base} -T 4\"","\n";
       $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(chrisamiller/sqanti3:v5.1.2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
       #$bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash -c \"source activate env && perl /bin/espresso/src/ESPRESSO_S.pl -L $input_tsv -F $h38_REF -A $f_gtf -O $RUNDIR -T 4\"","\n"; 
       print $bsub_com;
       system ($bsub_com);

}

## re-run espresso 
sub bsub_re_espresso{
    my $dir_expresso=$sample_full_path."/expresso1"; 
    if(!-d $dir_expresso)
    {
        `mkdir $dir_expresso`;
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
        my $input_update_tsv=$sample_full_path."/expresso1/".$chr1."/".$chr1.".bam.tsv.updated";    
        $current_job_file = "j7_expresso_r3_".$sample_name.".".$chr1.".sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_expresso."/".$chr1; 
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }
        my $input_tsv=$sample_full_path."/splitbam/".$chr1.".bam.tsv";
        my $input_update_tsv=$sample_full_path."/expresso1/".$chr1."/".$chr1.".bam.tsv.updated";
        open(ESPR, ">$job_files_dir/$current_job_file") or die $!;
        print ESPR "#!/bin/bash\n";
        print ESPR "RUNDIRfilter=".$dir_expresso."/".$chr1."/filtered\n";
        print ESPR "if [ ! -d \${RUNDIRfilter} ]\n";
        print ESPR "then\n";
        print ESPR "mkdir \${RUNDIRfilter}\n";
        print ESPR "fi\n";
        print ESPR "filtergtfin=".$sample_full_path."/sqanti/filter_ml_default/merged_N2_R0_updated.filtered.gtf\n";
        print ESPR "filtergtfout=".$dir_expresso."/".$chr1."/".$chr1."_N2_R0_updated.filtered.gtf\n";
        print ESPR "ISOFORM=".$dir_expresso."/".$chr1."/filtered/".$chr1."_compatible_isoform.tsv","\n";
        print ESPR "grep \"\^$chr1\[\[\:space\:\]\]\" \${filtergtfin} >  \${filtergtfout}","\n";
        print ESPR "source activate env && perl /bin/espresso/src/ESPRESSO_Q.pl -L $input_update_tsv -O \${RUNDIRfilter} -A \${filtergtfout} --read_ratio_cutoff 2 -V \${ISOFORM} -T 4","\n";
        close ESPR;
        my $sh_file=$job_files_dir."/".$current_job_file;
        $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
        print $bsub_com;
        system ($bsub_com);
    }        

}

## merge outputs from all chromosomes into abundance/gtf files

sub bsub_merge_espresso_filter{

    my $dir_expresso=$sample_full_path."/expresso1"; 

    if(!-d $dir_expresso)
    {
        `mkdir $dir_expresso`;
    }

   
    $current_job_file = "j8_merge_expresso_filter_".$sample_name.".sh"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

    `rm $lsf_out`;
    `rm $lsf_err`;
  
    open(MERGE, ">$job_files_dir/$current_job_file") or die $!;
    print MERGE "#!/bin/bash\n";
    print MERGE "RUNDIR=".$dir_expresso."\n";
    print MERGE "     ".$run_perl_script_path."merge_espresso_filter.pl \${RUNDIR} $chr_status\n";
    close MERGE;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    print $bsub_com;
    system ($bsub_com);         

}

#running espresso on the filtered gtf told us whether every read was FSM, ISM, etc.  Take that information and filter
#the espresso intermediate files so that they only contain FSM reads. Basically, copy the outputs directories and refill them
#with the same data with non-FSM reads removed.

sub bsub_filter_isoform{

    my $dir_expresso=$sample_full_path."/expresso1"; 
    
    if(!-d $dir_expresso)
    {
        `mkdir $dir_expresso`;
    }

    my @chrlist=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

    foreach my $chr (@chrlist)
    {
    	my $chr1=$chr;
   
        if($chr_status==1) 
		{ 
		$chr1="chr".$chr; 
        }
    
        $current_job_file = "j9_filter_isoform_".$sample_name.".".$chr1.".sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_expresso."/".$chr1; 
        
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }

        my $input_tsv=$sample_full_path."/splitbam/".$chr1.".bam.tsv";
        my $input_update_tsv=$sample_full_path."/expresso1/".$chr1."/".$chr1.".bam.tsv.updated";

        open(FILTER, ">$job_files_dir/$current_job_file") or die $!;
        print FILTER "#!/bin/bash\n";
        print FILTER "RUNDIR=".$dir_expresso."/"."$chr1\n";
        print FILTER "ISOFORM=".$dir_expresso."/".$chr1."/filtered/".$chr1."_compatible_isoform.tsv","\n";
        print FILTER "     ".$run_perl_script_path."remove_nonfsm_data.pl \${ISOFORM} \${RUNDIR} $input_update_tsv","\n";
        close FILTER;
        my $sh_file=$job_files_dir."/".$current_job_file;
        $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(chrisamiller/docker-genomic-analysis:latest)\' -o $lsf_out -e $lsf_err sh $sh_file\n"; 
        print $bsub_com;
        system ($bsub_com); 

    }

}


#okay, rerun espresso_q (again!) per chromosome, using these filtered files as inputs

## re-run espresso 
sub bsub_re2_espresso{
    my $dir_expresso=$sample_full_path."/expresso1"; 
    if(!-d $dir_expresso)
    {
        `mkdir $dir_expresso`;
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
        my $input_update_tsv=$sample_full_path."/expresso1/".$chr1."/fsm/".$chr1.".bam.tsv.updated";    
        $current_job_file = "j10_expresso_r3_fsm_".$sample_name.".".$chr1.".sh"; 
        my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
        my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
       `rm $lsf_out`;
       `rm $lsf_err`;
        my $RUNDIR=$dir_expresso."/".$chr1; 
        if(-d $RUNDIR)
        {
        `mkdir $RUNDIR`;    
        }
        my $input_tsv=$sample_full_path."/splitbam/".$chr1.".bam.tsv";
        my $input_update_tsv=$sample_full_path."/expresso1/".$chr1."/fsm/".$chr1.".bam.tsv.updated";
        open(ESPR2, ">$job_files_dir/$current_job_file") or die $!;
        print ESPR2 "#!/bin/bash\n";
        print ESPR2 "RUNDIRfsm=".$dir_expresso."/".$chr1."/fsm\n";
        print ESPR2 "filtergtfout=".$dir_expresso."/".$chr1."/".$chr1."_N2_R0_updated.filtered.gtf\n";
        print ESPR2 "ISOFORM=".$dir_expresso."/".$chr1."/fsm/".$chr1."_compatible_isoform.tsv","\n";
        print ESPR2 "source activate env && perl /bin/espresso/src/ESPRESSO_Q.pl -L $input_update_tsv -O \${RUNDIRfsm} -A \${filtergtfout} --read_ratio_cutoff 2 -V \${ISOFORM} -T 4","\n";
        close ESPR2;
        my $sh_file=$job_files_dir."/".$current_job_file;
        $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(sridnona/espresso:v2)\' -o $lsf_out -e $lsf_err /bin/bash $sh_file\n"; 
        print $bsub_com;
        system ($bsub_com);
    }        
}

sub bsub_merge_fsm{

    my $dir_expresso=$sample_full_path."/expresso1"; 

    if(!-d $dir_expresso)
    {
        `mkdir $dir_expresso`;
    }

   
    $current_job_file = "j11_merge_fsm_".$sample_name.".".".sh"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";

    `rm $lsf_out`;
    `rm $lsf_err`;
  
    open(MERGEFSM, ">$job_files_dir/$current_job_file") or die $!;
    print MERGEFSM "#!/bin/bash\n";
    print MERGEFSM "RUNDIR=".$dir_expresso."\n";
    print MERGEFSM "     ".$run_perl_script_path."merge_fsm.pl \${RUNDIR} $chr_status\n";
    close MERGEFSM;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 4 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    print $bsub_com;
    system ($bsub_com);         

}

