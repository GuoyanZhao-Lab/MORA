##########################################################################################
# If you used the same direcotry structure as the example, ~/tools/MORA-Singularity/ 
# then you only need to replace the following places: 
# <custom_query.fa> with your own file name
# <custom_Background.fa> with your own background sequence file name
# Change OUTDIR to the full direcotry path that you want.
# It is recommended to use N=100 to get stable results.
###########################################################################################

#!/bin/bash
INDIR=../CustomInput/
OUTDIR=../output_runCustomInput/

QUERYFASTA=<custom_query.fa>
BGGeneFASTA=<custom_Background.fa>
N=100

DBDIR=../database/CISBP_v2.00/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final
DBFILE=../database/CISBP_v2.00/CISBP_v2.00_HumanMouseRatCombined_QCed_DB_Final.txt
TEMPDIR=../temp/

NUMBERTHREAD=25
STEPNUMBER=0

perl ../MORA_Pipeline_v0.1/MORA_v0.1.pl ${INDIR} ${OUTDIR}   ${QUERYFASTA} ${BGGeneFASTA} ${N}  ${DBDIR} ${DBFILE}  ${TEMPDIR} ${NUMBERTHREAD} ${STEPNUMBER}
