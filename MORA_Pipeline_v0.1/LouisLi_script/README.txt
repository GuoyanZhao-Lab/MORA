README File for modified version of Yizhe MORA_v10 pipeline

Inputs:

$in_dir: the path of the directory where all of the gene files are located
ex: /home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/ 

$output_dir: the directory that the pipeline's output files get redirected to 
ex: /home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/output_dir_2/

$Specific_gene_file_fa: the file of the sepcific genes in .fa format
ex:/home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/JEM_IRF1dep_class1_class3.fa

$FullBgGeneSeqFile: the file path of the background genes in .fa format
ex: /home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/BackgroundGene

$N_simulation: the number of times a random number of genes is to be selected (also the number of .fa, .con, and ORI.txt files are produced). Optimal number is 100.

$Motif_Database_Dir: The path of the directory that contains all of the .matrix files
ex: /home/louisli/tools/Motif_Enrichment_MORA_pipeline/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final 

$Motif_Database_File: The path of the file containing matrix information
ex:/home/louisli/tools/Motif_Enrichment_MORA_pipeline/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt

$temp_dir: the directory that files will temporarily get written to
ex: /home/louisli/temp/

$software path: the path of the directory that contains the program
ex: /home/louisli/tools/Motif_Enrichment_MORA_pipeline/copy_Yizhe_MORA_v10/pipeline_2

$step_number: the step number that you want to run. Must be between 0 and 7. 0 will run all step


Louis Li example input: 
perl /home/louisli/tools/Motif_Enrichment_MORA_pipeline/copy_Yizhe_MORA_v10/pipeline_2/MORA_v10.pl /home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/ /home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/output_dir_2/ /home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/JEM_IRF1dep_class1_class3.fa /home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/BackgroundGene.fa 25 /home/louisli/tools/Motif_Enrichment_MORA_pipeline/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final /home/louisli/tools/Motif_Enrichment_MORA_pipeline/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt /home/louisli/temp/ /home/louisli/tools/Motif_Enrichment_MORA_pipeline/copy_Yizhe_MORA_v10/pipeline_2
