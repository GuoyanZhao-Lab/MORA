#! /usr/bin/perl -w
# copyright (C) 2006 Washington University, St. Louis, MO. #
# All Rights Reserved.                                     #
#                                                          #
# Author: Guoyan Zhao                                       #
# Send all comments to gzhao@ural.wustl.edu                #
# DISCLAIMER: THIS SOFTWARE IS PROVIDED "AS IS"            #
#             WITHOUT WARRANTY OF ANY KIND.                #
#----------------------------------------------------------#

use strict;
my $usage = '
perl script <input dir>

';


die $usage unless scalar @ARGV == 1;
my ( $in_dir ) = @ARGV;

my $run_script_path = `dirname $0`;
chomp $run_script_path;
$run_script_path = "perl ".$run_script_path."/";

opendir (DH, $in_dir) or die "Cannot open $in_dir";
foreach my $file (readdir DH) {
	if ($file =~ /\.fa/) {
		my $in_file  = $in_dir."/".$file;
		my $name = $file;
		$name =~ s/fa//;
		my $con_file  = $in_dir."/".$name."con";
		my $com = $run_script_path."S2_FastaToConsensus.pl  $in_file > $con_file\n";
		print $com;
		system ($com);
	}
}

close DH;
