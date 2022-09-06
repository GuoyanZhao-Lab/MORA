#! /usr/bin/perl

use strict;


my $usage = 'perl script <file1> <file2> <file3> <file4> <output file> <Motif database file>
<file1> = query gene sequence in consensus format, full path
<file2> = background gene sequence in consensus file, full path
<file3> = query gene sequence in fasta format, full path
<file4> = background gene sequence in fasta format, full path

<output file> = full path of output file
<Motif database file> = full path to the database file transformed
	/scratch/gzlab/gzhao/TF_Motif_database/CISBP/CISBP_v2.00_Homo_sapiens_2021_05_21/CISBP_v2.00_Homo_sapiens_2021_05_21_CountMatrix_PhyloNetFormat_Final.txt

<Motif database dir> = full path to the directory holding .matrix files (make sure as .matrix)
	/scratch/gzlab/gzhao/TF_Motif_database/CISBP/CISBP_v2.00_Homo_sapiens_2021_05_21/CISBP_v2.00_Homo_sapiens_2021_05_21_CountMatrix_PhyloNetFormat_Final

<temp dir> = full path to a directory that you can write temporary file
     /scratch/gzlab/gzhao/temp/patser/

This script will perform following analysis:
#1. read in all the matrices
#2. for each matrix, run patser over the specific genes and non-specific genes
#3. parse patser output file, calculate various numbers used to evaluate the motif 
#4. rank the motifs 
#5. output to indicated output file

';

die $usage  unless @ARGV == 8;

#the location of Patser, alphabet file
my $run_script_path = `dirname $0`;
chomp $run_script_path;
$run_script_path = $run_script_path."/";

my ($specificGeneFile_con, $BgGeneFile_con, $specificGeneFile_fa, $BgGeneFile_fa, $outputFile, $Motif_Database, $matrix_dir, $tmp_dir) = @ARGV;

# may need to change the directory. should be where all the matrix files reside
#my $matrix_dir = "/scratch/gzlab/databases/CISBP_v2.00_Homo_sapiens_2021_05_21_CountMatrix/";

# open (OutFile, ">$outputFile") or die "Can not open output file $outputFile !\n";
open (OutFile, ">",$outputFile) or die "Can not open output file $outputFile !\n";
my $patser_out_Specific_file = $tmp_dir."/patserOutSpecific".$$;
my $patser_out_NS_file = $tmp_dir."/patserOutNS".$$;
my $pseudo=1;

# First, read in all the matrices
####################################################################
# read in CISBP matrix data file
my @AC = ();
my %TF_ID; # Motif_ID => TF_ID
my %Family_ID; # Motif_ID => TF family ID
my %ENSEMBLE_ID; # Motif_ID => TF ENSEMBLE DBID
my %TF_Name; # Motif_ID => TF name
my %TF_Species; # TF_Species
my %Family_Name; # Family_Name
my %Width;
my %AlignMatrixTF;  #accession number => alignment matrix
&read_CISBP_matrix($Motif_Database);

##print "total number of CISBP matrices: ", scalar @AC, "\n";
#foreach my $key (keys %BFac) {
#   print "\\$key\\\t", "@{$BFac{$key}}","\n";
#}
# now we've got all the information of matrices 

############################################################################
#second, count how many genes are in the specific gene file and control file
my $totalSpeGenes = `grep ">" $specificGeneFile_fa |wc -l`;
my $totalNSgenes = `grep ">" $BgGeneFile_fa |wc -l`;

print OutFile "total specific genes: $totalSpeGenes\n";
print OutFile "total NS genes: $totalNSgenes\n";

# calculate length for all the query sequences 
my @LengthOfSpecSequences = (); # length of query sequences
#my $numOfSeq = 0;
#my @bin = ();
#my $range = 50;
#my $numBin = 15;

my $oldseperator = $/;
$/ = ">";
open (FastaFile, $specificGeneFile_fa) or die "Can't Open FASTA file: $specificGeneFile_fa";
while (my $line = <FastaFile>){ # read in one fasfa record
	# Discard blank lines
	if ($line =~ /^\s*$/) {
	    next;
	}
	# discard the first line which only has ">", keep the rest
	elsif ($line ne ">") {
#		$numOfSeq++;
		chomp $line;
		my @rows = ();
	    @rows = split (/\n/, $line);	
	    my $read_name = shift @rows;
	    my $read_seq = join("", @rows);
	    $read_seq =~ s/\s//g; #remove white space
		my $Speclength = length($read_seq);
		push @LengthOfSpecSequences, $Speclength;
	}
}
$/ = $oldseperator;
close ($specificGeneFile_fa);

my ($Speclengthave, $Specsd) = &calculate_average_SD(\@LengthOfSpecSequences);
my @sorted = sort {$a <=> $b}  @LengthOfSpecSequences;
#print "total num of query genes	  = ", scalar @sorted, "\n";;
#printf("avgSequenceSize   = %d\n", $Speclengthave);

# get information for background sequences
my @LengthOfNSSequences = (); 
#my $numOfSeq = 0;
#my @bin = ();
#my $range = 50;
#my $numBin = 15;

my $oldseperator = $/;
$/ = ">";
open (FastaFile, $BgGeneFile_fa) or die "Can't Open FASTA file: $BgGeneFile_fa";
while (my $line = <FastaFile>){ # read in one fasfa record
	# Discard blank lines
	if ($line =~ /^\s*$/) {
	    next;
	}
	# discard the first line which only has ">", keep the rest
	elsif ($line ne ">") {
#		$numOfSeq++;
		chomp $line;
		my @rows = ();
	    @rows = split (/\n/, $line);	
	    my $read_name = shift @rows;
	    my $read_seq = join("", @rows);
	    $read_seq =~ s/\s//g; #remove white space
		my $NSlength = length($read_seq);
		push @LengthOfNSSequences, $NSlength;
	}
}
$/ = $oldseperator;
close ($BgGeneFile_fa);

my ($NSlengthave, $NSsd) = &calculate_average_SD(\@LengthOfNSSequences);
my @sorted = sort {$a <=> $b}  @LengthOfNSSequences;
#print "total num of background genes	  = ", scalar @sorted, "\n";;
#printf("avgSequenceSize   = %d\n", $NSlengthave);

#######################################################################
# now for each matrix, run patser over specific genes and control genes and 
# calculate all the numbers

# for specific genes
my %TGeneHitSpec;   # matrixName => total number of genes that have at least 1 
                   # of this matrix
my %TMotifSpec;     # matrixName => total motifs found in this group of genes 
my %TScoreSpec;    # matrixName => total score of motifs that has >=cutoff score
my %TMaxScoreSpec;  # matrixName => total maximum score for the motif
my %percentSpec;     # matrixName => percent of genes that has >=1 motif
my %aveNumMotifSpec; #matrixName => average num of motifs/muscle_gene
my %aveScoreSpec;   # matrixName => average score for all sites
my %aveMaxScoreSpec; #matrixName => average maximum score for this motif

# for background genes
my %TGeneHitBG;   # matrixName => total genes that have at least 1 of this matrix
my %TMotifBG;     # matrixName => total motifs found in this group 
my %TScoreBG;       # matrixName => total score of motifs that has >=cutoff score
my %TMaxScoreBG;       # matrixName => total maximum score of the motif 
my %percentTotalBG;     # matrixName => percent of genes that has >=1 motif
my %aveNumMotifBG; #matrixName => average num of motifs/gene
my %aveScoreBG;   # matrixName => average score for all sites
my %aveMaxScoreBG; #matrixName => average maximum score for this motif

foreach my $matrix (@AC){
#    print "for matrix $matrix \n";
#        print "start of loop";
	my $currentgene = ""; #the current gene being read in the patser file
	my $maxScore = 0;    # maximum score for one gene
	my $cutoffScore = 100;

	my $TGeneHitSpec = 0;   # total genes that have at least 1 of this matrix
	my $TMotifSpec = 0;     # total motifs found in this group
	my $TScoreSpec = 0;     # total score of motifs that has >=cutoff score
	my $TMaxScoreSpec = 0;  # total maximum score of the motif

	my $matrix_temp_file = $matrix_dir."/".$matrix.".matrix";

	my $com = $run_script_path."patser_v3d_linux -m $matrix_temp_file -f".$specificGeneFile_con." -a ".$run_script_path."alphabet -c -li -d2 > $patser_out_Specific_file";
#	print "com is $com\n";
	system $com;
    
	$com = $run_script_path."patser_v3d_linux -m $matrix_temp_file -f".$BgGeneFile_con." -a ".$run_script_path."alphabet -c -li -d2 > $patser_out_NS_file";
	system $com;
    
	my $MaxScore; # maximum score for this motif
	open (MusPatser, "< $patser_out_Specific_file") or die "Can not open file $patser_out_Specific_file !\n";
    # obtain cutoff score
	while (my $line = <MusPatser>){
		chomp $line;
		if ($line =~ /maximum score:/) {
			chomp $line;
			$line =~ s/maximum score://;
			$MaxScore = $line;
			$MaxScore =~ s/\s*//;
#			print "Max score is $MaxScore\n"
		}

		if ($line =~ /numerically calculated cutoff score:/) {
			chomp $line;
			$line =~ s/numerically calculated cutoff score://;
			$cutoffScore = $line;
			$cutoffScore =~ s/\s*//;
#			print "numerically calculated cutoff score is $cutoffScore\n\n";
	    	last;
		}
    }

 	my $currentGene = "";
	while (my $line = <MusPatser>){
#		print "line is $line";
#		if ( /\s+(.*?)\s+position=\s*(\d+C*)\s*score=\s*([\d\.]+)/ ) {
#		if ($line =~ /\s+(.*?)\s+position=\s*(\d+C*)\s*score=\s*([\d\.]+)/ ) {
		if ($line =~ /\s*(.*?)\s+position=\s*(\d+C*)\s*score=\s*([\d\.]+)/ ) {
#               print "name = \\$1\\, pos = \\$2\\, score = \\$3\\\n";
                my $geneName = $1;
				my $pos = $2;
				my $score = $3;
				if ($pos =~ /C/) {
					$pos =~ s/C//;
				}
  
			# Only use sites with score higher than the cutoff score
			if ($score >= $cutoffScore) {  
				#calculate total motifs found in this group of genes
				$TMotifSpec++;
		
				# calculate total score of motifs that has >=cuoff score
				$TScoreSpec += $score;

				# calculate total genes that have >=1 this motif TGeneHit
				# total maximum score for this motif TMaxScore
				if ($geneName ne $currentGene) {
					$currentGene = $geneName;
					$TGeneHitSpec++;
					$TMaxScoreSpec += $maxScore; 	
					$maxScore = $score;
				}
				else {
		    		if ($score > $maxScore) {
						$maxScore = $score;
					}
				}
			}
		}
    }
	   
    $TGeneHitSpec{$matrix} = $TGeneHitSpec;
    $TMotifSpec{$matrix} = $TMotifSpec;
    $TScoreSpec{$matrix} = $TScoreSpec;
    $TMaxScoreSpec{$matrix} = $TMaxScoreSpec;
    $percentSpec{$matrix} = ($TGeneHitSpec+$pseudo)/($totalSpeGenes+$pseudo)*100;
    $aveNumMotifSpec{$matrix} = $TMotifSpec/$totalSpeGenes;
    $aveScoreSpec{$matrix} = ($TScoreSpec+$pseudo)/($TMotifSpec+$pseudo);
    $aveMaxScoreSpec{$matrix} = ($TMaxScoreSpec+$pseudo)/($TGeneHitSpec+$pseudo);
    close MusPatser;

	# ##########################################################################################
    # now calculate all the numbers for the control group
	my $TGeneHitBG = 0;    # total genes that have >=1 this matrix
	my $TMotifBG = 0;     # total motifs found in this group
	my $TScoreBG = 0;     # total score of motifs that has >cutoff score
	my $TMaxScoreBG = 0;  # total maximum score of the motif

	open (Patser, "< $patser_out_NS_file ") or die "Can not open file $patser_out_Specific_file !\n";
	$currentGene = "";
	while (my $line = <Patser>){
#		if ($line =~ /\s+(.*?)\s+position=\s*(\d+C*)\s*score=\s*([\d\.]+)/ ) {
		if ($line =~ /\s*(.*?)\s+position=\s*(\d+C*)\s*score=\s*([\d\.]+)/ ) {
#               print "name = \\$1\\, pos = \\$2\\, score = \\$3\\\n";
                my $geneName = $1;
				my $pos = $2;
				my $score = $3;
				if ($pos =~ /C/) {
					$pos =~ s/C//;
				}

	    	# Only use sites with score higher than the cutoff score
	    	if ($score >= $cutoffScore) {  
				#calculate total motifs found in this group of genes
				$TMotifBG++;
		
				# calculate total score of motifs that has >cuoff score
				$TScoreBG += $score;

				# calculate total genes that have >= 1 this motif TGeneHit
				# total maximum score for this motif TMaxScore
				if ($geneName ne $currentGene) {
				    $currentGene = $geneName;
				    $TGeneHitBG++;
				    $TMaxScoreBG += $maxScore; 
				    $maxScore = $score;
				}
				else {
				    if ($score > $maxScore) {
						$maxScore = $score;
				    }
				}
	    	}
		}
	}
	   
    $TGeneHitBG{$matrix} = $TGeneHitBG;
    $TMotifBG{$matrix} = $TMotifBG; 
    $TScoreBG{$matrix} = $TScoreBG;
    $TMaxScoreBG{$matrix} = $TMaxScoreBG;
    $percentTotalBG{$matrix} = ($TGeneHitBG+$pseudo)/($totalNSgenes+$pseudo)*100;
    $aveScoreBG{$matrix} = ($TScoreBG+$pseudo)/($TMotifBG+$pseudo); # if 0, add sudo count 1
    $aveMaxScoreBG{$matrix} = ($TMaxScoreBG+$pseudo)/($TGeneHitBG+$pseudo);# if 0, add sudo count 1
    
    close Patser;
}

my %ratio1; # matrixName => $percentSpec/$percentTotalBG, measure the 
            # enrichment of the motif in query genes
my %ORI; # Over-representation index

for my $matrix (@AC) {
		$ratio1{$matrix} = $percentSpec{$matrix}/$percentTotalBG{$matrix};
		$ORI{$matrix} = (($TMotifSpec{$matrix}+$pseudo)*$NSlengthave*$totalNSgenes*$ratio1{$matrix})/(($TMotifBG{$matrix}+$pseudo)*$Speclengthave*$totalSpeGenes);
}

# print output to output file
#print OutFile "Here are the ranked motifs :\n\n";

print OutFile "#spec: number of specific genes that has the motif\n";
print OutFile "%spec: percent of specific genes that has the motif\n";
print OutFile "nonS: non-specific gene group \n";
print OutFile "ORI: over-representation index\n";

print OutFile "Motif_ID\t#spec\t#nonS\t%spec\t%nonS\tratio\tORI\tTF_Name\tENSEMBLE_ID\tTF_ID\tTF_Species\tFamily_Name\tFamily_ID\twidth\n";

foreach my $key (sort {$ORI{$b} <=> $ORI{$a}} keys %ORI) {
#foreach my $key (sort {$ratio1{$b} <=> $ratio1{$a}} keys %ratio1) {
    printf OutFile "%s\t", $key;
    printf OutFile ("%d", $TGeneHitSpec{$key});
    print OutFile "\t";
    printf OutFile ("%d",  $TGeneHitBG{$key});
    print OutFile "\t";
    printf OutFile ("%.2f", $percentSpec{$key});
    print OutFile "\t";
    printf OutFile ("%.2f", $percentTotalBG{$key});
    print OutFile "\t";
    printf OutFile ("%.2f", $ratio1{$key});
    print OutFile "\t";
    printf OutFile ("%.2f", $ORI{$key});
    print OutFile "\t";
	print OutFile  $TF_Name{$key};
    print OutFile "\t";
	print OutFile  $ENSEMBLE_ID{$key};
    print OutFile "\t";
	print OutFile  $TF_ID{$key};
    print OutFile "\t";
	print OutFile  $TF_Species{$key};
    print OutFile "\t";
	print OutFile  $Family_Name{$key};
    print OutFile "\t";
	print OutFile  $Family_ID{$key};
    print OutFile "\t";
	print OutFile  $Width{$key};
    print OutFile "\n";

}

unlink($patser_out_Specific_file);
unlink($patser_out_NS_file);


print OutFile "#Finished\n";
exit;

###################################################################
# This subroutine reads in CISBP matrix data file
sub read_CISBP_matrix {
    my $file = shift @_;
    open (InFile, "<$file") or die "Can not open CISBP matrix file $file!\n";
    my $oldSeperator = $/;
    $/ = "//\n";
    while (my $record = <InFile>) {
#		print $record;
		my $Motif_ID = '';
		if ($record =~ /^Motif_ID = (.*)/m) {
			$Motif_ID = $1;
		}

		my $TF_Name = '';
		if ($record =~ /^TF_Name = (.*)/m){
			$TF_Name = $1;
		}
		
		my $ENSEMBLE_ID = '';
		if ($record =~ /^ENSEMBLE_ID = (.*)/m) {
			$ENSEMBLE_ID = $1;
		}

		my $TF_ID = '';
		if ($record =~ /^TF_ID = (.*)/m) {
			$TF_ID = $1;
		}

		my $TF_Species = '';
		if ($record =~ /^TF_Species = (.*)/m) {
			$TF_Species = $1;
		}

		my $Family_Name = '';
		if ($record =~ /^Family_Name = (.*)/m) {
			$Family_Name = $1;
		}

		my $Family_ID = '';
		if ($record =~ /^Family_ID = (.*)/m) {
			$Family_ID = $1;
		}
			
		my $width = '';
		if ($record =~ /^width = (\d+)/m){
			$width = $1;
		}
#		print "width = \\$width\\\n";

		if ($record =~ /^XX\n((.*\n)*?)^XX\n/m){
			my $alignMatrix = $1;
#			print "matrix is \n$alignMatrix\n";	    
			$AlignMatrixTF{$Motif_ID} = $alignMatrix;
		}

		push (@AC, $Motif_ID);
		$TF_ID{$Motif_ID} = $TF_ID;
		$Family_ID{$Motif_ID} = $Family_ID;
		$ENSEMBLE_ID{$Motif_ID} = $ENSEMBLE_ID;
		$TF_Name{$Motif_ID} = $TF_Name;
		$TF_Species{$Motif_ID} = $TF_Species;
		$Family_Name{$Motif_ID} = $Family_Name;
		$Width{$Motif_ID} = $width;
	}
    $/ = $oldSeperator;
}


####################################################################
sub calculate_average_SD {
	my ($arr_ref) = @_;

	# calculate average, SD
	my $sum = 0;
	my $n = scalar @{$arr_ref};
	for (my $i = 0; $i < $n; $i++) {
		$sum += $arr_ref->[$i];
	}
	my $average = $sum/$n;
	my $variance = 0;
	for (my $i = 0; $i < $n; $i++) {
		$variance += ($arr_ref->[$i] - $average)*($arr_ref->[$i] -  $average);
	}
	$variance = $variance/$n;
	my $SD = sqrt($variance);
	return ($average, $SD);
}

