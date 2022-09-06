use strict;
use warnings;

# The following model was set to utilize multiple CPU thread. Here, we specify to run 25 jobs at the same time. 
# # Therefore, this script is optimized to run motif enrichment analysis on server at Dr. Zhao's lab
# # If you need to use this script locally, you may want to change $MAX_PROCESSES to a smaller value
# # If you need to use this script on other compute server/cluster, you need to contact the author and use the other version of MORA.pl. This is because clusters usually handle multiple jobs differently
use lib '/home/louisli/tools/MORA_Singularity_v0.1/Parallel-ForkManager-2.02/lib';
use lib '/home/louisli/tools/MORA_Singularity_v0.1/Parallel-ForkManager-2.02/lib/Storable-3.25';

use Parallel::ForkManager;
my $MAX_PROCESSES=25;
my $pm =Parallel::ForkManager->new ($MAX_PROCESSES);


#usage information
my $usage = '
perl script <input dir> <output dir> <DEgene.fa> <BGgene.fa> <Number of simulation> <Motif database> <temp dir> <software path>

<input dir> 
<output dir> 
<DEgene.fa> 
<BGgene.fa> 
<Number of simulation> = 25-100
<Motif database> = /home/gzhao/tools/Motif_Enrichment_MORA_pipeline/JASPAR2018_CORE_Hocomoco_CORE_vertebrates_non-redundant.transformed 
<temp dir> = /home/gzhao/temp/
<Software path> = /home/gzhao/tools/Motif_Enrichment_MORA_pipeline/Yizhe_MORA_v10/
';

# <Motif database> storage information
# # 1. To use JASPAR2018_CORE_vertebrates_non-redundant & Hocomoco_CORE combined
# # <Motif database> = "/home/ysong/data/Motif_Database/JASPAR2018_CORE_Hocomoco_CORE_vertebrates_non-redundant.transformed"
# # 2. To use JASPAR2018_CORE_vertebrates_non-redundant
# # <Motif database> = "/home/ysong/data/Motif_Database/JASPAR2018_CORE_vertebrates_non-redundant_pfms_transfac.transformed"
# # 3. To use JASPAR2018_CORE_non-redundant
# # <Motif database> = "/home/ysong/data/Motif_Database/JASPAR2018_CORE_non-redundant_pfms_transfac.transformed"
# # 4. To use TRANSFAC 10.2
# # <Motif database> "/home/ysong/data/Motif_Database/Transfac10.2matrix.transformed"

die $usage unless scalar @ARGV == 10;
my ($in_dir, $output_dir, $Sepcific_gene_file_fa, $FullBgGeneSeqFile, $N_simulation, $Motif_Database_Dir, $Motif_Database_File, $temp_dir, $software_path, $step_number) = @ARGV;
$, = "\t";
$\ = "\n";

#####################################################################################
## Create directory by date and time
##sub mkdir_by_date{
#    #print local time
#    #   my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
#    #    my $mytime = sprintf("%04d%02d%02d_%02d%02d%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec);
#    #    printf("DateTime: %s", $mytime);
#    #    system("mkdir user_${mytime}_hello");
#    #}
#    #$output_dir=mkdir_by_date;
#    #####################################################################################
#    # directory that holding the TRANSFAC matrix description files

my $TRANSFAC_Dir = "/home/data/transfac/10.2/matrixDir/";

my $run_script_path = `dirname $0`;

#print "run script path = \\$run_script_path\\\n";
chomp $run_script_path;
$run_script_path = "perl ".$run_script_path."/";

#print "run script path = \\$run_script_path\\\n";
#
## in sample directory, convert all input files from ".fa" file to ".con" file.

opendir (DH, $in_dir) or die "Cannot open $in_dir";
foreach my $file (readdir DH) {
	my $pid =$pm->start and next;
	if ($file =~ /\.fa/) {
		my $in_file  = $in_dir."/".$file;
		my $name = $file;
		$name =~ s/fa//;
		my $con_file  = $in_dir."/".$name."con";
                ############## Script Path Here ##############
		my $com = $run_script_path. "/S2_FastaToConsensus.pl  $in_file > $con_file\n";
                #############################################

    print "fa_file name : ", $in_file;
    print "con_file name : ", $con_file;
    
		system ($com);
	}
	$pm->finish;
}
$pm->wait_all_children;
print "......All fasta format files have been converted to consensus format files";
close DH;

#!perl
my $error;
my %Subset_sequences = ();
my @Subset_seq_length = ();
my @sorted_gene_length;
my $geneNum = $#Subset_seq_length;
my $min;
my $max;
my %Fullset_sequences = ();
my @Fullset_seq_name = ();
my $totalNumGene;


if ($step_number == 0){
    get_sequences();
    run_simulation();
    fasta_to_con();
    calculate_ORI();   
    check_output();
    generate_table();
    generate_p_value();
}
elsif($step_number == 1){
    get_sequences();
}
elsif($step_number == 2){
    get_sequences();
    run_simulation();
}
elsif($step_number == 3){
    fasta_to_con();
}
elsif($step_number == 4){
    calculate_ORI();
}
elsif($step_number == 5){
    check_output();
}
elsif($step_number == 6){
    generate_table();
}
elsif($step_number == 7){
    generate_p_value();
}

else{
    die $usage;
}
exit;

sub get_sequences{
       #my $error	
	mkdir $output_dir or $error = $!;
	unless (-d $output_dir) {
   		 die "Cannot create directory '$output_dir': $error";
  	  }
       # my %Subset_sequences = ();
        &getGeneSeq($Sepcific_gene_file_fa, \%Subset_sequences);
       # my @Subset_seq_length = ();
         foreach my $key (keys %Subset_sequences) {
               push @Subset_seq_length, length($Subset_sequences{$key});
            }
         $geneNum = $#Subset_seq_length;
         @sorted_gene_length = sort {$a <=> $b} @Subset_seq_length;
         $min = $sorted_gene_length[0];
         $max = $sorted_gene_length[$#sorted_gene_length];
        # %Fullset_sequences = ();
        &getGeneSeq($FullBgGeneSeqFile, \%Fullset_sequences);
       # my @Fullset_seq_name = ();
        foreach my $key (keys %Fullset_sequences) {
                 push @Fullset_seq_name, $key;
                    }
        $totalNumGene = $#Fullset_seq_name + 1;

     }

######################################################################
## randome sample $N_simulation times of subset sequences from fullset sequences
## and output each set into given directory.

sub run_simulation{
    srand(time|$$);
    my $FileName_base = $output_dir."/Randome_set";
    for (my $i = 1; $i <= $N_simulation; $i++) { # randome sample $N_simulation times
        my $FileName = $FileName_base.$i.".fa";
        open (OUT, ">$FileName") or die "Can not open file $FileName\n";
#	foreach my $gene (keys %Subset_sequences) { # each gene
        for (my $j = 0; $j <= $#Subset_seq_length; $j++) {
		
                 my $arrPos = int(rand($totalNumGene));
                 
                 print OUT ">$Fullset_seq_name[$arrPos]\n";
                 print OUT $Fullset_sequences{$Fullset_seq_name[$arrPos]}, "\n";
 }
}
    print "......Randome sequence sets generatation complete";
}

sub fasta_to_con{
    #######################################################################
    #    # fasta file to consensus file
        opendir (DH, $output_dir) or die "Cannot open $output_dir";
    foreach my $file (readdir DH) {
        my $pid =$pm->start and next;
        if ($file =~ /\.fa/) {
            my $in_file  = $output_dir."/".$file;
            my $name = $file;
            $name =~ s/fa//;
            my $con_file  = $output_dir."/".$name."con";
         ################ Script Path Here ####################  
	 my $com = $run_script_path. "/S2_FastaToConsensus.pl  $in_file > $con_file\n";
         ######################################################
         #
        print "fa_file name : ", $in_file;
        print "con_file name : ", $con_file;
         system ($com);
        }
        $pm->finish;
    }
    $pm->wait_all_children;
    close DH;
}

sub calculate_ORI{
    #################################################################################################
    #    # calculate ORI for each transcription factor matrix for each randome sequence set
    print $output_dir, "\n";
    my $HOME = $ENV{HOME};
    # =head
    print STCH "#!/bin/bash\n";         
    print STCH "IN=$output_dir/Randome_set\${SLURM_ARRAY_TASK_ID}.con\n"; # full path
    print STCH "OUT=$output_dir/Randome_set\${SLURM_ARRAY_TASK_ID}_TRANSFAC_motif_ORI.txt\n"; # full path

    print STCH 'if [ -s $IN ]',"\n"; #check if the file is empty
    print STCH "then\n";
    print STCH "	 $run_script_path/S3_MatrixAnalysis_Rank_Motif_By_ORI_TRANSFAC_Motifs.pl $TRANSFAC_Dir $Sepcific_gene_file_fa  \${IN} \${OUT} 1 \n";
    print STCH "fi";
    close STCH;
   # =cut
    my $Query_prefix = "JEM_IRF1dep_class1_class3"; #make this an input once pipeline works

    my $INDIR ="$in_dir"."/".$Query_prefix."_RandomSet/ \n";
    my $QUERYCON = "$in_dir".$Query_prefix.".con";
    my $QUERYFASTA = "$in_dir".$Query_prefix.".fa";
   # my $RANDOMCON = "$in_dir"."/".$Query_prefix."_RandomSet/Random_set\${SLURM_ARRAY_TASK_ID}.con\n";
    my $RANDOMFASTA = "$in_dir"."/".$Query_prefix."_RandomSet/Random_set\${SLURM_ARRAY_TASK_ID}.fa\n";
   # my $OUT = "$in_dir"."/".$Query_prefix."_RandomSet/Random_set\${SLURM_ARRAY_TASK_ID}_motif_ORI.txt\n";
   
   # my $Motif_Database_dir = "/home/louisli/tools/Motif_Enrichment_MORA_pipeline";
    my $name = $Sepcific_gene_file_fa;
   # my $Motif_Database_File = "/home/louisli/tools/Motif_Enrichment_MORA_pipeline/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt";
    $name =~ s/fa//;
    my $Sepcific_gene_file_con = $name."con";
    for (my $i = 1; $i <=$N_simulation; $i++) {
        my $pid =$pm->start and next;
        ################################ Script Path Here ########################
        my $OUT =$output_dir ."Randome_set" . $i ."_motif_ORI.txt" ;
        print "OUT";
        print $OUT;
        my $RANDOMCON = $output_dir . "Randome_set" .$i . ".con" ;
        my $RANDOMFASTA = $output_dir . "Randome_set" .$i . ".fa" ;
        
  # my $com ="perl $run_script_path/S4_Rank_Motif_By_ORI.pl \${QUERYCON} \${RANDOMCON} \${QUERYFASTA}  \${RANDOMFASTA} \${OUT} $Motif_Database $Motif_Database_dir $temp_dir \n";   
        # my $com = $run_script_path."/S4_Rank_Motif_By_ORI.pl "."$QUERYCON".".$RANDOMCON."$QUERYFASTA"."$RANDOMFASTA".$OUT $Motif_Database $Motif_Database_dir $temp_dir\n" ;
     # my $com = $run_script_path."/S4_Rank_Motif_By_ORI.pl "." $QUERYCON"." ".$RANDOMCON.""." $QUERYFASTA"." ".$RANDOMFASTA. " ".$output_dir."/Randome_set".$i."_TRANSFAC_motif_ORI.txt  1  $Motif_Database $Motif_Database_dir $temp_dir \n";
     my $com = "perl " . $software_path . "/S4_Rank_Motif_By_ORI.pl " . $QUERYCON . " " . $RANDOMCON . " ".$QUERYFASTA ." " .$RANDOMFASTA ." " .  $OUT ." " .  $Motif_Database_File ." ".$Motif_Database_Dir ." " .  $temp_dir;    

# /home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/JEM_IRF1dep_class1_class3.con
# /home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/JEM_IRF1dep_class1_class3.fa
# /home/louisli/tools/Motif_Enrichment_MORA_pipeline/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final
##########################################################################
     # print "Jack Towe got a 2";	
        print $com;
       # print $Sepcific_gene_file_con;
        system($com);
        $pm->finish;
    }
    $pm->wait_all_children;
}

################################################################################
###  Step 5: check to make sure the output for each file is completed.
sub check_output{
	my $com = $run_script_path .  "/S5_Rank_Motif_By_ORI_checkoutput.pl " . $output_dir;
        system($com);   
}

################################################################################
###  Step 6: generate motif occurrance table

sub generate_table{
    my $out_file = $Sepcific_gene_file_fa;
    $out_file =~ s/fa//;
    my $output_file = "$output_dir/motif_occurrance_table.csv";
    my $com = $run_script_path."/S6_Generate_Motif_occurrance_table.pl $output_dir $N_simulation $output_file";
    system($com);
}

################################################################################
###  Step 7: calling R from perl to calculate motif occurrence p value
##my $base = `pwd`;
##chomp $base;
##my $r_script= "/home/ysong/MORA/data_for_motif_analysis/MORA_calculate_pvalue.R";
##my $path="$base/$r_script";
##my $path="/home/ysong/MORA/data_for_motif_analysis/MORA_calculate_pvalue.R";
sub generate_p_value{
    my $path=$software_path."/S7_Calculate_pvalue.R";
    motif_p();

    sub motif_p {
       # my $path=$software_path."/";
       # my $execute = `Rscript $path $N_simulation $output_dir`;
        my $execute ="Rscript " . $path . " " . $N_simulation . " " . $output_dir . " test_output" ;
       # print $? if $?;
         
        print $execute; 
       # print $execute,"\n";
        system( $execute);
        print "all finshed";   
 }

}

#######################################################################
## read fasta file

sub getGeneSeq() {
    my ($fileName, $hashRef) =  @_;
    open (File, $fileName) or die "Can not open file $fileName !\n";
    # get all the gene Name and seq name=> seq
    my $oldSeperator = $/;
    $/ = ">"; # seperates two records
    my $line = <File>;
    while ($line = <File>){
		my @fields = split ("\n", $line);
		my $name = shift @fields;
		my @temp = split (/\s+/, $name);
		$name = shift @temp;
		my $desc = join ("", @temp); 
		my $seq = join ("", @fields);
		$seq =~ s/>//g;
		$seq =~ s/\s//g;
		$hashRef->{$name} = $seq;
    }
    $/ = $oldSeperator;
	close File;
}
