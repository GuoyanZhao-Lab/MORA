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
INDIR=~/MORA-Singularity/CustomInput/
OUTDIR=~/MORA-Singularity/output_runCustomInput/

QUERYFASTA=GSE56026_STAT1_human_KD_Down_GeneNames.fa
BGGeneFASTA=GSE56026_STAT1_human_KD_Down_Backgrounds.fa
N=100

DBDIR=~/MORA-Singularity/database/CISBP_v2.00/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final
DBFILE=~/MORA-Singularity/database/CISBP_v2.00/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt
TEMPDIR=~/MORA-Singularity/temp/

NUMBERTHREAD=25
STEPNUMBER=0

perl ~/MORA-Singularity/MORA_Pipeline_v0.1/MORA_v0.1.pl ${INDIR} ${OUTDIR}   ${QUERYFASTA} ${BGGeneFASTA} ${N}  ${DBDIR} ${DBFILE}  ${TEMPDIR} ${NUMBERTHREAD} ${STEPNUMBER}
