##---------------------------------------------------------------------------##
##  Author:
##      Guoyan Zhao <gzhao@pathology.wustl.edu>
#******************************************************************************
#* Copyright (C) Washington University in St. Louis Developed by Guoyan Zhao,
#* Washington University in St. Louis, Department of Pathology and Immunology.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php.
#*
###############################################################################


#!/usr/bin/perl
use strict;

my $usage = '
This script will get the _motif_ORI.txt for every data set 
in the given directory and generate a motif occurrance table as below:
motif ORIofRandomSet1	2	3	...	occuranceInRefSeqSet

perl script  <run folder><output file>
<input folder> = full path of the folder holding files with ORI info
               without the last "/"
<Number of simulation> = 100
<output file> = file name  

';
die $usage unless scalar @ARGV == 3;
my ( $dir, $N_RandomSet, $output_file ) = @ARGV;

#my $HOME = $ENV{HOME};

my $run_script_path = `dirname $0`;
chomp $run_script_path;
$run_script_path = "/usr/bin/perl ".$run_script_path."/";

my %Motif_ID_to_NumSpec = ();
my %Motif_ID_to_NumNonSpec = ();
my %Motif_ID_to_PctSpec = ();
my %Motif_ID_to_PctNonSpec = ();
my %Motif_ID_to_ORIs = ();
my %Motif_ID_to_TF_Name = ();
my %Motif_ID_to_ENSEMBLE_ID = ();
my %Motif_ID_to_TF_ID = ();
my %Motif_ID_to_TF_Species = ();
my %Motif_ID_to_Family_Name = ();
my %Motif_ID_to_Family_ID = ();
my %Motif_ID_to_width = ();
my %Motif_ID_to_Average_ORI = ();
my $num_files = 0;

for (my $i = 1; $i <= $N_RandomSet; $i++) {
	my $file = "Randome_set".$i."_motif_ORI.txt";
#	print "processing file $file \n\n";
	my $full_path = $dir."/".$file;
	open (IN, "$full_path") or die "Can not open file $full_path!\n";

	# name  #spec   #nons   %spec   %nons   ratio   ori      TF_Name ENSEMBLE_ID     TF_ID   TF_Species      Family_Name     Family_ID       width
	 foreach (1..9) {
		my $line =    <IN>;
   	 }
	while (my $line = <IN>) {
		if ($line =~/#Finished/) {
			last;
		}
		#		print $line;
		chomp $line;
		my @temp = split (/\t/, $line);
		my $Motif_ID = $temp[0];
		$Motif_ID =~ s/\s+//g;
		#		print "name = \\$Motif_ID\\, #spec = \\$temp[1]\\,#nonSpec = \\$temp[2]\\,%spec = \\$temp[3]\\,%nonSpec = \\$temp[4]\\, ORI=\\$temp[6]\\, TF_Name = \\$temp[7]\\, ENSEMBLE_ID=\\$temp[8]\\, TF_ID = \\$temp[9], TF_Species=\\$temp[10]\\, Family_Name = \\$temp[11], Family_ID=\\$temp[12]\\, width = \\$temp[13]\n\n";
		my $desc = join("\t", $temp[7], $temp[8], $temp[9], $temp[10], $temp[11], $temp[12], $temp[13]);

		if (defined $Motif_ID_to_NumSpec{$Motif_ID}) {
			push @{$Motif_ID_to_NumSpec{$Motif_ID}}, $temp[1];
		}
		else {
			${$Motif_ID_to_NumSpec{$Motif_ID}}[0] = $temp[1];
		}

		if (defined $Motif_ID_to_NumNonSpec{$Motif_ID}) {
			push @{$Motif_ID_to_NumNonSpec{$Motif_ID}}, $temp[2];
		}
		else {
			${$Motif_ID_to_NumNonSpec{$Motif_ID}}[0] = $temp[2];
		}

		if (defined $Motif_ID_to_PctSpec{$Motif_ID}) {
			push @{$Motif_ID_to_PctSpec{$Motif_ID}}, $temp[3];
		}
		else {
			${$Motif_ID_to_PctSpec{$Motif_ID}}[0] = $temp[3];
		}

		if (defined $Motif_ID_to_PctNonSpec{$Motif_ID}) {
			push @{$Motif_ID_to_PctNonSpec{$Motif_ID}}, $temp[4];
		}
		else {
			${$Motif_ID_to_PctNonSpec{$Motif_ID}}[0] = $temp[4];
		}

		if (defined $Motif_ID_to_ORIs{$Motif_ID}) {
			push @{$Motif_ID_to_ORIs{$Motif_ID}}, $temp[6];
		}
		else {
			${$Motif_ID_to_ORIs{$Motif_ID}}[0] = $temp[6];
		}

		if (!(defined $Motif_ID_to_TF_Name{$Motif_ID})) {
			$Motif_ID_to_TF_Name{$Motif_ID} = $temp[7];
		}

		if (!(defined $Motif_ID_to_ENSEMBLE_ID{$Motif_ID})) {
			$Motif_ID_to_ENSEMBLE_ID{$Motif_ID} = $temp[8];
		}

		if (!(defined $Motif_ID_to_TF_ID{$Motif_ID})) {
			$Motif_ID_to_TF_ID{$Motif_ID} = $temp[9];
		}

		if (!(defined $Motif_ID_to_TF_Species{$Motif_ID})) {
			$Motif_ID_to_TF_Species{$Motif_ID} = $temp[10];
		}

		if (!(defined $Motif_ID_to_Family_Name{$Motif_ID})) {
			$Motif_ID_to_Family_Name{$Motif_ID} = $temp[11];
		}

		if (!(defined $Motif_ID_to_Family_ID{$Motif_ID})) {
			$Motif_ID_to_Family_ID{$Motif_ID} = $temp[12];
		}

		if (!(defined $Motif_ID_to_width{$Motif_ID})) {
			$Motif_ID_to_width{$Motif_ID} = $temp[13];
		}

	}
	close IN;
}	

open(OUT, ">$output_file") or die "Can not open $output_file!\n";

# print out the table, matrix name in rows
#foreach my $key (keys %Motif_ID_to_ORIs) {
#	print OUT $key, "\t";
#	print OUT "@{$Motif_ID_to_ORIs{$key}}\n"; 
#}
print OUT "Motif_ID\t";
for (my $i = 1; $i <= $N_RandomSet; $i++) {
	print OUT "ORI_Random$i\t";
}
print OUT "NumRefGene\tPctRefGene\tNumRandomGene\tPctRandomGene\tORI_ave\tORI_SD\tTF_Name\tENSEMBLE_ID\tTF_ID\tTF_Species\tFamily_Name\tFamily_ID\twidth\n";

foreach my $key (keys %Motif_ID_to_ORIs) {
#	print "$key, times selected as overrepresented, ", scalar @{$Motif_ID_to_ORIs{$key}}, "total number of files: $N_RandomSet \n";
#	print "specific genes @{$Motif_ID_to_NumSpec{$key}}\n";
	my ($NumSpec_ave, $NumSpec_sd) = &calculate_average_SD(\@{$Motif_ID_to_NumSpec{$key}});
	my ($NumNonSpec_ave, $NumNonSpec_sd) = &calculate_average_SD(\@{$Motif_ID_to_NumNonSpec{$key}});
	my ($PctSpec_ave, $PctSpec_sd) = &calculate_average_SD(\@{$Motif_ID_to_PctSpec{$key}});
	my ($PctNonSpec_ave, $PctNonSpec_sd) = &calculate_average_SD(\@{$Motif_ID_to_PctNonSpec{$key}});
	my ($ORI_ave, $ORI_SD) = &calculate_average_SD(\@{$Motif_ID_to_ORIs{$key}});

	print OUT $key, "\t";
	for (my $i = 0; $i < $N_RandomSet; $i++) {
		print OUT $Motif_ID_to_ORIs{$key}[$i],"\t";
	}
       	print OUT  $NumSpec_ave, "\t", $PctSpec_ave, "\t", $NumNonSpec_ave, "\t", $PctNonSpec_ave, "\t", $ORI_ave, "\t", $ORI_SD,"\t", $Motif_ID_to_TF_Name{$key}, "\t", $Motif_ID_to_ENSEMBLE_ID{$key}, "\t", $Motif_ID_to_TF_ID{$key}, "\t", $Motif_ID_to_TF_Species{$key}, "\t", $Motif_ID_to_Family_Name{$key}, "\t", $Motif_ID_to_Family_ID{$key}, "\t", $Motif_ID_to_width{$key}, "\n"; 
}

exit;

#################################################################
sub calculate_average_SD {
	my ($arr_ref) = @_;

	# calculate average, SD
	my $sum = 0;
	my $n = scalar @{$arr_ref};
	for (my $i = 0; $i < $n; $i++) {
#		print "element is \\$arr_ref->[$i]\\\n";
		$sum += $arr_ref->[$i];
	}
	my $average = $sum/$n;
	#print "number of elements: $n, sum = $sum\n";
	my $variance = 0;
	for (my $i = 0; $i < $n; $i++) {
		$variance += ($arr_ref->[$i] - $average)*($arr_ref->[$i] - $average);
	}
	$variance = $variance/$n;
	my $SD = sqrt($variance);
	return ($average, $SD);
}
