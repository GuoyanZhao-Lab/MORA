#!/bin/bash
INDIR=~/data/Zhao_Guoyan/Paper_MORA/MORA-Singularity/ExampleInput/
OUTDIR=~/data/Zhao_Guoyan/Paper_MORA/MORA-Singularity/output_ExampleOutput_UsingTestDatabase/

QUERYFASTA=JEM_IRF1dep_class1_class3.fa
BGGeneFASTA=BackgroundGene.fa
N=25

DBDIR=~/data/Zhao_Guoyan/Paper_MORA/MORA-Singularity/database/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final
DBFILE=~/data/Zhao_Guoyan/Paper_MORA/MORA-Singularity/database/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt
TEMPDIR=~/data/Zhao_Guoyan/Paper_MORA/MORA-Singularity/temp/

NUMBERTHREAD=25
STEPNUMBER=0

perl ~/data/Zhao_Guoyan/Paper_MORA/MORA-Singularity/MORA_Pipeline_v0.1/MORA_v0.1.pl ${INDIR} ${OUTDIR}   ${QUERYFASTA} ${BGGeneFASTA} ${N}  ${DBDIR} ${DBFILE}  ${TEMPDIR} ${NUMBERTHREAD} ${STEPNUMBER}
