################################################################################
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

perl script <sample dir>
<dir> = /scratch/dwlab/gzhao/data/AzadBonni/Xiaoying/Motif_enrichment_analysis/SnoNKO_DownRegulated_RandomSet

';
die $usage unless scalar @ARGV == 1;
my ( $dir ) = @ARGV;

my $HOME = $ENV{HOME};

my $All_file_finished = 1;
opendir(DH, $dir) or die "Can not open dir $dir!\n";
foreach my $file (readdir DH) {
	if ($file =~ /\.fa$/) { # DNA sequence
		my $NotFinish = 0;

		my $prefix = $file;
		$prefix =~ s/\.fa$//; 
		my $out = $dir."/".$prefix."_motif_ORI.txt";
                
                
		if (!(-s $out)) {
#		        print "top condition";
			$NotFinish = 1;
		}
		else {
			my $result = `tail -5 $out`;
			if (!($result =~ /#Finished/) ) { # did not finish
#                                print "bottom condition";
				$NotFinish = 1;
			}
		}
		if ($NotFinish) {
			$All_file_finished = 0;
			print "$out did not finish correctly!\n\n";
		}
	}
}

if ($All_file_finished) {
	print "All files finished correctly.\n";
}


exit;

