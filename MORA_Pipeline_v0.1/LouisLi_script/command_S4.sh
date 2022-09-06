INDIR=/home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/JEM_IRF1dep_Class1_Class3_testOutput/:wq

QUERYCON=/home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/JEM_IRF1dep_class1_class3.con
QUERYFASTA=/home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/JEM_IRF1dep_class1_class3.fa
RANDOMCON=/home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/BackgroundGene.con
RANDOMFASTA=/home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/BackgroundGene.fa
OUT=/home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_TestData_IRF1/new_pipeline_output/output_file

DBFILE=/home/louisli/tools/Motif_Enrichment_MORA_pipeline/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt 
DBDIR=/home/louisli/tools/Motif_Enrichment_MORA_pipeline/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final
TEMPDIR=/home/louisli/new_temp/

/home/louisli/tools/Motif_Enrichment_MORA_pipeline/Motif_Enrichment_analysis_Pipeline_CISBP/S4_Rank_Motif_By_ORI.pl ${QUERYCON} ${RANDOMCON} ${QUERYFASTA}  ${RANDOMFASTA} ${OUT} ${DBFILE} ${DBDIR} ${TEMPDIR}
