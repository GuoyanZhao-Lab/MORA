#!/bin/bash
INDIR=~/tools/MORA_Singularity_v0.1/ExampleInput/
OUTDIR=~/tools/MORA_Singularity_v0.1/output_ExampleOutput_UsingTestDatabase/

QUERYFASTA=~/tools/MORA_Singularity_v0.1/ExampleInput/JEM_IRF1dep_class1_class3.fa
BGGeneFASTA=~/tools/MORA_Singularity_v0.1/ExampleInput/BackgroundGene.fa
N=25

DBDIR=~/tools/MORA_Singularity_v0.1/database/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final
DBFILE=~/tools/MORA_Singularity_v0.1/database/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt
TEMPDIR=~/tools/MORA_Singularity_v0.1/temp/

NUMBERTHREAD=25
STEPNUMBER=0

perl ~/tools/MORA_Singularity_v0.1/MORA_Pipeline_v0.1/MORA_v0.1.pl ${INDIR} ${OUTDIR}   ${QUERYFASTA} ${BGGeneFASTA} ${N}  ${DBDIR} ${DBFILE}  ${TEMPDIR} ${NUMBERTHREAD} ${STEPNUMBER}
