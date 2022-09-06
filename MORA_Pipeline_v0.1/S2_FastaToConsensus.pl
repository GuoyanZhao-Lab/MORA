#!/usr/bin/perl

use strict;
use constant FASTAwidth => 80;

my $usage = 'perl FastaToConsensus.pl <file1> 
file1 = file name that has the fasta seq
output to standard output
';

die $usage unless scalar @ARGV == 1;

my ($faFile) = @ARGV;
open (FaFile, "<$faFile") or die "can not open file $faFile";

my $oldSeperator = $/;
$/ = '>';
my $line = <FaFile>;
while ( $line = <FaFile>) {
    if ($line =~ /^\s*$/) {
	next;
    }

    my @lines = split (/\n/, $line);
    my $firstLine = shift @lines;
    #print "firstLine = $firstLine";

    my @fields = split (/\s/, $firstLine);
    my $name = shift @fields;
    my $seq = join ('', @lines);
    $seq =~ s/>//;
    $seq =~ s/\s//g;

    #print "name is \\$name\\ \n";
    #print "seq is $seq \n";
    print  $name, " \\";
    for (my $pos=0; $pos<length($seq); $pos += FASTAwidth){
	print substr($seq, $pos, FASTAwidth);
	print "\n";
    }
    print "\\\n";
}

exit;




