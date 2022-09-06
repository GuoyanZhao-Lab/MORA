##########################################################################################
# If you used the same direcotry structure as the example, ~/tools/MORA-Singularity/ 
# then you only need to replace the following places: 
# ExampleInput with your own directory name
# JEM_IRF1dep_class1_class3.fa with your own file name
# BackgroundGene.fa with your own background sequence file name
# Change OUTDIR to the full direcotry path that you want.
# It is recommended to use N=100 to get stable results.
###########################################################################################

#!/bin/bash
INDIR=~/tools/MORA-Singularity/ExampleInput/
OUTDIR=~/tools/MORA-Singularity/output_runExampleInput/

QUERYFASTA=~/tools/MORA-Singularity/ExampleInput/JEM_IRF1dep_class1_class3.fa
BGGeneFASTA=~/tools/MORA-Singularity/ExampleInput/BackgroundGene.fa
N=100

DBDIR=~/tools/MORA-Singularity/database/CISBP_v2.00/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final
DBFILE=~/tools/MORA-Singularity/database/CISBP_v2.00/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt
TEMPDIR=~/tools/MORA-Singularity/temp/

NUMBERTHREAD=25
STEPNUMBER=0

perl ~/tools/MORA-Singularity/MORA_Pipeline_v0.1/MORA_v0.1.pl ${INDIR} ${OUTDIR}   ${QUERYFASTA} ${BGGeneFASTA} ${N}  ${DBDIR} ${DBFILE}  ${TEMPDIR} ${NUMBERTHREAD} ${STEPNUMBER}
