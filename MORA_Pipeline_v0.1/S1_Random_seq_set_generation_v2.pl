
#!/usr/bin/perl -w

use strict;
#use constant TOOLS => "/home2/gzhao/tools/";

my $usage = '
#######################################################################
# This program read in a file contains a subset of sequences and a file 
# that contains all the sequences. 
#######################################################################

Usage: perl Random_seq_set_generation.pl <sequence subset File>
            <sequence full set file> <N times>

<sequence subset File> = 
	/scratch/dwlab/gzhao/data/ValeriaCavalli/201804_DanCarlin_TSC2_TF_project/Tsc2KO_vs_uninjured_upregulated_Intergenic5kb.fa

<sequence full set file> =
	/lts/dwlab/Virgin_Wang_shared/gzhao/mouse_transcription/Ensemble_data/Raw_data/Mouse/Mus_musculus.GRCm38.Intergenic5kb_Isoformconsolidated.fa

<N times> =   N (any number > 1). random sampling N times
<output dir> = full path to the output directory

';

die $usage unless @ARGV == 4;
my ($SubsetSeqFile, $FullsetSeqFile, $numOfSim, $output_dir) = @ARGV;

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

my $FileName_base = $output_dir."/Randome_set";
for (my $i = 1; $i <= $numOfSim; $i++) { # randome sample $numOfSim times
	my $FileName = $FileName_base.$i.".fa";
	open (OUT, ">$FileName") or die "Can not open file $FileName\n"; 
#	foreach my $gene (keys %Subset_sequences) { # each gene
	for (my $j = 0; $j <= $#Subset_seq_length; $j++) { 
#		print "the $j gene ****************************\n";
		# find a random gene in the full set sequence array
    	my $arrPos = int(rand($totalNumGene));
		print OUT ">$Fullset_seq_name[$arrPos]\n";
		print OUT $Fullset_sequences{$Fullset_seq_name[$arrPos]}, "\n";
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


