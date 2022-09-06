#! /bin/bash
echo Number of Simulations:
read num_Sims

echo Step number:
read num_step

perl /home/louisli/tools/enrichment_pipeline/MORA_Singularity_v0.1/MORA_Pipeline_v0.1/MORA_Singularity_v0.1.pl /home/louisli/tools/enrichment_pipeline/MORA_Singularity_v0.1/ExampleInput/ /home/louisli/tools/enrichment_pipeline/MORA_Singularity_v0.1/output/ /home/louisli/tools/enrichment_pipeline/MORA_Singularity_v0.1/ExampleInput/JEM_IRF1dep_class1_class3.fa /home/louisli/tools/enrichment_pipeline/MORA_Singularity_v0.1/ExampleInput/BackgroundGene.fa $num_Sims /home/louisli/tools/enrichment_pipeline/MORA_Singularity_v0.1/database/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final /home/louisli/tools/enrichment_pipeline/MORA_Singularity_v0.1/database/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt /home/louisli/tools/enrichment_pipeline/MORA_Singularity_v0.1/temp /home/louisli/tools/enrichment_pipeline/MORA_Singularity_v0.1/MORA_Pipeline_v0.1/ $num_step
