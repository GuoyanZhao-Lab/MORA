
#!/usr/bin/perl -w

use strict;

my $usage = '
#######################################################################
# This program read in a file contains a subset of sequences and a file 
# that contains all the sequences. 
# It will simulate the length distribution of the subset of sequences 
# for n times assuming the same length distribution .
#
#######################################################################

Usage: perl S1_Random_seq_set_generation_v1_controlLength.pl <Query sequence File> <background sequence file> <N times><output dir>

<Query sequence File> =  
	/home/ysong/data/Motif_Enrichment_test_data/SnoNKO_UpRegulated_gene.fa

<background sequence file> =
	/home/ysong/data/Motif_Enrichment_test_data/SnoNKO_unchanged_gene.fa

<N times> =  number of random sampling data to generate (e.g. 100)
<output dir> = full path to the output dire

';

die $usage unless @ARGV == 4;
my ($SubsetSeqFile, $FullsetSeqFile, $numOfSim, $output_dir) = @ARGV;

#create output folder if not exist 
if (! -d $output_dir) {
    `mkdir $output_dir`;
}

# get sequences that are used in module prediction
my %Subset_sequences = ();
&getGeneSeq($SubsetSeqFile, \%Subset_sequences);
#foreach my $key (keys %Subset_sequences) {
#	print $key,"\n";
#	print $Subset_sequences{$key}, "\n";
#}

my @Subset_seq_length = ();
foreach my $key (keys %Subset_sequences) {
	push @Subset_seq_length, length($Subset_sequences{$key});
}
my $geneNum = $#Subset_seq_length;
#print "total number in sub set $geneNum\n";

my @sorted_gene_length = sort {$a <=> $b} @Subset_seq_length;
my $min = $sorted_gene_length[0];
my $max = $sorted_gene_length[$#sorted_gene_length];
#print "gene length min = $min, max = $max\n";

my %Fullset_sequences = ();
&getGeneSeq($FullsetSeqFile, \%Fullset_sequences);

my @Fullset_seq_name = ();
foreach my $key (keys %Fullset_sequences) {
	push @Fullset_seq_name, $key;
}
my $totalNumGene = $#Fullset_seq_name + 1;
#print "total number in full set $totalNumGene\n";

######################################################################
# randome sample $numOfSim times of subset sequences from fullset sequences
# and output each set into given directory.
srand(time|$$);

my $FileName_base = $output_dir."/Random_set";
for (my $i = 1; $i <= $numOfSim; $i++) { # randome sample $numOfSim times
	my $FileName = $FileName_base.$i.".fa";
	open (OUT, ">$FileName") or die "Can not open file $FileName\n"; 
#	foreach my $gene (keys %Subset_sequences) { # each gene
	for (my $j = 0; $j <= $#Subset_seq_length; $j++) { 
#		print "the $j gene ****************************\n";
		# find a random gene in the full set sequence array
#	    my $upperbound = $Subset_seq_length[$j]*1.10;
#	    my $lowerbound = $Subset_seq_length[$j]*0.90;
		my $upperbound = $Subset_seq_length[$j]*1.20;
		my $lowerbound = $Subset_seq_length[$j]*0.80;
		my $success = 0;
		do {
	    	my $arrPos = int(rand($totalNumGene));
	    	my $RandomeSeqLength = length($Fullset_sequences{$Fullset_seq_name[$arrPos]});
			if (($RandomeSeqLength >= $lowerbound) && ($RandomeSeqLength <= $upperbound)) {
				print OUT ">$Fullset_seq_name[$arrPos]\n";
				print OUT $Fullset_sequences{$Fullset_seq_name[$arrPos]}, "\n";
				$success = 1;
			}
		} while (!$success);
	}
}

exit;

#######################################################################
# read fasta file
sub getGeneSeq() {
    my ($fileName, $hashRef) =  @_;
    open (File, $fileName) or die "Can not open file $fileName !\n";

    # get all the gene Name and seq name=> seq
    my $oldSeperator = $/;
    $/ = ">"; # seperates two records
    my $line = <File>;
    while ($line = <File>){
		#print $line;
		my @fields = split ("\n", $line);
		my $name = shift @fields;
		my @temp = split (/\s+/, $name);
		$name = shift @temp;
		my $desc = join ("", @temp); 
		my $seq = join ("", @fields);
		$seq =~ s/>//g;
		$seq =~ s/\s//g;
		#lc $seq;
 
		#print "name is \\$name\\\n";
		#print "seq is \n";
		#print "\\$seq\\\n";
		$hashRef->{$name} = $seq;
    }
    $/ = $oldSeperator;
	close File;
}


