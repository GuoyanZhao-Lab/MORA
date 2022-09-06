#!/bin/bash
INDIR=~/tools/MORA-Singularity/ExampleInput/
OUTDIR=~/tools/MORA-Singularity/output_ExampleOutput_UsingTestDatabase/

QUERYFASTA=~/tools/MORA-Singularity/ExampleInput/JEM_IRF1dep_class1_class3.fa
BGGeneFASTA=~/tools/MORA-Singularity/ExampleInput/BackgroundGene.fa
N=25

DBDIR=~/tools/MORA-Singularity/database/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final
DBFILE=~/tools/MORA-Singularity/database/CISBP_v2.00_ForTesting/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt
TEMPDIR=~/tools/MORA-Singularity/temp/

NUMBERTHREAD=25
STEPNUMBER=0

perl ~/tools/MORA-Singularity/MORA_Pipeline_v0.1/MORA_v0.1.pl ${INDIR} ${OUTDIR}   ${QUERYFASTA} ${BGGeneFASTA} ${N}  ${DBDIR} ${DBFILE}  ${TEMPDIR} ${NUMBERTHREAD} ${STEPNUMBER}
