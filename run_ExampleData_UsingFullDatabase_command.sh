# It takes about 8.5 min to use the full database (>3000 transcription factor PWM).


#!/bin/bash
INDIR=../ExampleInput/
OUTDIR=../output_ExampleOutput_UsingFullDatabase/

QUERYFASTA=JEM_IRF1dep_class1_class3.fa
BGGeneFASTA=BackgroundGene.fa
N=25

DBDIR=../database/CISBP_v2.00/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final
DBFILE=../database/CISBP_v2.00/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt
TEMPDIR=../temp/

NUMBERTHREAD=25
STEPNUMBER=0

perl ../MORA_Pipeline_v0.1/MORA_v0.1.pl ${INDIR} ${OUTDIR}   ${QUERYFASTA} ${BGGeneFASTA} ${N}  ${DBDIR} ${DBFILE}  ${TEMPDIR} ${NUMBERTHREAD} ${STEPNUMBER}
