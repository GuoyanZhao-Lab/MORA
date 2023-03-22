Motif Over Representation Analysis (MORA) is a fully automated pipeline to detect enrichment of known TF binding motifs in any query sequences. This repository provides source code and documentation of the single machine version of the pipeline. 

1) Make sure your system has software-properties-common and libparallel-forkmanager-perl installed.
If not following command would install the required packages.
 
apt-get update && apt-get install -y software-properties-common
add-apt-repository universe
apt-get update
apt-get install -y libparallel-forkmanager-perl


2) e.g. in your ~/tools/ directory, issue the following command:
git clone https://github.com/GuoyanZhao-Lab/MORA-Singularity

A directory named "MORA-Singularity" will appear. Change directory to "ExampleInput" using following command:

cd MORA-Singularity/ExampleInput

2) issue the following command to run the pipeline using input files in the "ExampleInput" 
directory and the test database (a smaller version of the database).
./run_ExampleData_UsingTestDatabase_command.sh

You will see the following outputs:

run script path = /home/gzhao/tools/MORA-Singularity/MORA_Pipeline_v0.1

Input directory:        /home/gzhao/tools/MORA-Singularity/ExampleInput/

......All fasta format files have been converted to consensus format files
......Randome sequence sets generatation complete
Output directory:       /home/gzhao/tools/MORA-Singularity/output_ExampleOutput_UsingTestDatabase/

All files finished correctly.
Rscript /home/gzhao/tools/MORA-Singularity/MORA_Pipeline_v0.1/S7_Calculate_pvalue.R 25 /home/gzhao/tools/MORA-Singularity/output_ExampleOutput_UsingTestDatabase/ Output

all finshed
Final output file: /home/gzhao/tools/MORA-Singularity/output_ExampleOutput_UsingTestDatabase//Output_MotifEnriched_25Resampling.csv

A new directory named "output_ExampleOutput_UsingTestDatabase" will appear in the "MORA-Singularity" directory. The intermediate output files and the final output file will be in this directory

3) issue the following command within the "ExampleInput" directory to run the pipeline using input files in the "ExampleInput" 
directory and the actual database CISBP_2.00. (contains > 3000 transcription factor's position weight matrices). 
It takes about 8.5 mins to finish the analysis.

./run_ExampleData_UsingFullDatabase_command.sh
run script path = /home/gzhao/tools/MORA-Singularity/MORA_Pipeline_v0.1

Input directory:        /home/gzhao/tools/MORA-Singularity/ExampleInput/

......All fasta format files have been converted to consensus format files
......Randome sequence sets generatation complete
Output directory:       /home/gzhao/tools/MORA-Singularity/output_ExampleInput_FullDatabase/

All files finished correctly.
Rscript /home/gzhao/tools/MORA-Singularity/MORA_Pipeline_v0.1/S7_Calculate_pvalue.R 25 /home/gzhao/tools/MORA-Singularity/output_ExampleInput_FullDatabase/ Output

all finshed
Final output file: /home/gzhao/tools/MORA-Singularity/output_ExampleInput_FullDatabase//Output_MotifEnriched_25Resampling.csv

A new directory "output_ExampleInput_UsingFullDatabase" will appear in the "MORA-Singularity" directory. The intermediate and final output files will be in that directory. 

4) You can modify the file run_MyOwnData_command.sh to run analysis on your own data in the "CustomInput" directory. 
