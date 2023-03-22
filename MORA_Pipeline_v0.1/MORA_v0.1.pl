###############################################################################################################
# MORA: Motif Over Representation Analysis
# Contact: Guoyan Zhao
# Email: gzhao@wustl.edu
# GNU General Public License
# Developer: Louis Li, Yizhe Song, Guoyan Zhao
# Version 0.1
###################################################################################
#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;

###################################################################################
#usage information
my $usage = '
perl script <input dir> <output dir> <query_gene_file.fa> <FullBgGeneSeqFile.fa> <N_simulation> <Motif_Database_Dir><Motif_Database_File> <temp_dir> <MAX_PROCESSES> <step_number>

<in_dir>: the path of the directory where all of the gene files are located
example: ~/tools/MORA_Singularity_v0.1/ExampleInput/

<output_dir>: the directory that the pipeline output files will be saved to. Must be new directory not existed in the system. Will create automatically.
example: ~/tools/MORA_Singularity_v0.1/output/

<query_gene_file.fa>: full path to the file of the query genes in .fa format
example: ~/tools/MORA_Singularity_v0.1/ExampleInput/JEM_IRF1dep_class1_class3.fa

<FullBgGeneSeqFile.fa>: the file path of the background genes in .fa format
example: ~/tools/MORA_Singularity_v0.1/ExampleInput/BackgroundGene.fa

<N_simulation>: the number of sets of random background genes to be generated. Recommended number: 100
example: 100

<Motif_Database_Dir>: The path of the directory that contains all of the .matrix files
example: ~/tools/MORA_Singularity_v0.1/database/CISBP_v2.00/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final

<Motif_Database_File>: The path of the file containing matrix information
example: ~/tools/MORA_Singularity_v0.1/database/CISBP_v2.00/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt

<$temp_dir>: the directory that files will temporarily get written to. Will create automatically if not exist.
example: ~/tools/MORA_Singularity_v0.1/temp/

<MAX_PROCESSES>: number of CPU thread to use. The more the faster. 
example: 25

<step_number>: the step number that you want to run. Must be between 0 and 7. 0 will run all step

';

die $usage unless scalar @ARGV == 10;
my ($in_dir, $output_dir, $Sepcific_gene_file_fa, $FullBgGeneSeqFile, $N_simulation, $Motif_Database_Dir, $Motif_Database_File, $temp_dir, $MAX_PROCESSES, $step_number) = @ARGV;
$, = "\t";
$\ = "\n";

################################################################################
# number of CPU thread
# my $MAX_PROCESSES=25;
my $pm =Parallel::ForkManager->new ($MAX_PROCESSES);

################################################################################
my $run_script_path = `dirname $0`;
chomp $run_script_path;
print "run script path = $run_script_path\n";

################################################################################
# if $output_dir does not exist, create one
if (! -d $output_dir) {
	`mkdir $output_dir`;
}
else {
	print "Output directory exist! \n";
	exit;
}

# if $temp_dir does not exist, create one
if (! -d $temp_dir) {
	`mkdir $temp_dir`;
}

################################################################################
## in sample directory, convert all input files from ".fa" file to ".con" file.
print "Input directory: ", $in_dir, "\n";

opendir (DH, $in_dir) or die "Cannot open $in_dir";
foreach my $file (readdir DH) {
	my $pid =$pm->start and next;
	if ($file =~ /\.fa/) {
		my $in_file  = $in_dir."/".$file;
		my $name = $file;
		$name =~ s/fa//;
		my $con_file  = $in_dir."/".$name."con";
                ############## Script Path Here ##############
		my $com = "perl ".$run_script_path. "/S2_FastaToConsensus.pl  $in_file > $con_file\n";
                #############################################
		#   print "fa_file name : ", $in_file;
		#   print "con_file name : ", $con_file;
    
		system ($com);
	}
	$pm->finish;
}
$pm->wait_all_children;
print "......All fasta format files have been converted to consensus format files";
close DH;

####################################################################################
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


if ($step_number == 0){ # run the full pipeline fully automated
    get_sequences();
    run_simulation();
    fasta_to_con();
    calculate_ORI();   
    check_output();
    generate_table();
    generate_p_value();
}
elsif($step_number == 1){ # run step-by-step
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

##############################################################
# subroutines
##############################################################
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
        for (my $j = 0; $j <= $#Subset_seq_length; $j++) {
		
                 my $arrPos = int(rand($totalNumGene));
                 
                 print OUT ">$Fullset_seq_name[$arrPos]\n";
                 print OUT $Fullset_sequences{$Fullset_seq_name[$arrPos]}, "\n";
 }
}
    print "......Randome sequence sets generatation complete";
}

##########################################
#  Convert FASTA format file to con format
##########################################
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
	 my $com = "perl ".$run_script_path. "/S2_FastaToConsensus.pl  $in_file > $con_file\n";
         ######################################################
	 #        print "fa_file name : ", $in_file;
	 #        print "con_file name : ", $con_file;
         system ($com);
        }
        $pm->finish;
    }
    $pm->wait_all_children;
    close DH;
}

####################
# Calculate ORI
# ##################
sub calculate_ORI{
    #################################################################################################
    #    # calculate ORI for each transcription factor matrix for each randome sequence set
    print "Output directory: ", $output_dir, "\n";
    my $HOME = $ENV{HOME};
    (my $Query_prefix = $Sepcific_gene_file_fa) =~ s/\.fa//; #make this an input once pipeline works
    my $QUERYCON = "$in_dir".$Query_prefix.".con";
    my $QUERYFASTA = "$in_dir".$Query_prefix.".fa";
   
    my $name = $Sepcific_gene_file_fa;
    $name =~ s/fa//;
    my $Sepcific_gene_file_con = $name."con";
    for (my $i = 1; $i <=$N_simulation; $i++) {
        my $pid =$pm->start and next;
        ################################ Script Path Here ########################
        my $OUT =$output_dir ."Randome_set" . $i ."_motif_ORI.txt" ;
	#        print $OUT;
        my $RANDOMCON = $output_dir . "Randome_set" .$i . ".con" ;
        my $RANDOMFASTA = $output_dir . "Randome_set" .$i . ".fa" ;
        
        my $com = "perl ".$run_script_path . "/S4_Rank_Motif_By_ORI.pl " . $QUERYCON . " " . $RANDOMCON . " ".$QUERYFASTA ." " .$RANDOMFASTA ." " .  $OUT ." " .  $Motif_Database_File ." ".$Motif_Database_Dir ." " .  $temp_dir;    

##########################################################################
#        print $com;
        system($com);
        $pm->finish;
    }
    $pm->wait_all_children;
}

################################################################################
###  Step 5: check to make sure the output for each file is completed.
sub check_output{
	my $com = "perl ".$run_script_path ."/S5_Rank_Motif_By_ORI_checkoutput.pl ".$output_dir;
        system($com);   
}

################################################################################
###  Step 6: generate motif occurrance table
sub generate_table{
    my $out_file = $Sepcific_gene_file_fa;
    $out_file =~ s/fa//;
    my $output_file = "$output_dir/motif_occurrance_table.csv";
    my $com = "perl ".$run_script_path."/S6_Generate_Motif_occurrance_table.pl $output_dir $N_simulation $output_file";
    system($com);
}

################################################################################
###  Step 7: calling R from perl to calculate motif occurrence p value
sub generate_p_value{
    motif_p();

    sub motif_p {
        my $execute ="Rscript ".$run_script_path."/S7_Calculate_pvalue.R " .$N_simulation ." ".$output_dir." Output" ;
       # print $? if $?;
         
        print $execute,"\n";
        system( $execute);
        print "all finshed";
     	print "Final output file: ".$output_dir."/Output_MotifEnriched_".$N_simulation."Resampling.csv\n";	
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
