#!/usr/bin/env Rscript

# command line argument: numOfRandomSet MotifOccurranceTable OutputFilePrefix
#if (length(args) != 3) {
#	stop("Please provide numOfRandomSet MotifOccurranceTable OutputFilePrefix.")
#}

# Pass the command line argument to R
args = commandArgs(trailingOnly=TRUE)

## read the motif_occurrance_table.csv in output directory
motifs_up <- read.csv(paste(args[2],"motif_occurrance_table.csv",sep = "/"), header=TRUE, sep="\t")
motifs_up$Motif_ID <- as.character(motifs_up$Motif_ID)
rownames(motifs_up) <- motifs_up$Motif_ID

# calculate motif occurrence p value
p_value <-list()
# define try function
my.t.test.p.value <- function(data_vector, exp_value) {
    obj<-try(t.test(data_vector, mu = exp_value, alternative = "greater", conf.level=0.95), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj$p.value)
 }
for(i in 1:nrow(motifs_up)){
  data_a = as.numeric(motifs_up[i, 2:(as.numeric(args[1]) + 1)] )
  #test for normality
  #  shapiro.test(data_a)
#  ttestlist <- t.test(data_a, mu = motifs_up[i, "NumRefGene"], alternative = "two.sided", conf.level=0.99)
  #p_value[[i]] <- my.t.test(log10(data_a), mu = log10(1), alternative = "two.sided", conf.level=0.99)$p.value
  #p_value[[i]] <- my.t.test.p.value(log10(data_a), log10(1))
	# ORI > 1.2 defined as over represented
  p_value[[i]] <- my.t.test.p.value(data_a, 1.2 )
#  print(t.test(data_a, mu = motifs_up[i, "NumRefGene"], alternative = "two.sided", conf.level=0.99)$p.value)
}

df_p_value <- data.frame(matrix(unlist(p_value), nrow=nrow(motifs_up), byrow=T))
colnames(df_p_value) <- c("p")
motifs_up$p_ttest <- df_p_value$p

# calculate adjusted p value
motifs_up$padj <- p.adjust(motifs_up$p_ttest, method = "bonferroni", n = length(motifs_up$p_ttest))
#motifs_up <- motifs_up[order(-motifs_up$ORI_ave),]
motifs_up_enriched1 <- subset(motifs_up, padj < 0.05 )
motifs_up_enriched <- subset(motifs_up_enriched1, ORI_ave > 1.2 )
# at least 10% query squence has the motif site
motifs_up_enriched <- subset(motifs_up_enriched, PctRefGene > 10 )

# calculate combined score = percent reference genes having the motif * average ORI score
motifs_up_enriched$CombinedScore <- motifs_up_enriched$PctRefGene * motifs_up_enriched$ORI_ave
# rank motifs
motifs_up_enriched$Rank <- rank(-motifs_up_enriched$CombinedScore,ties.method= "min")
#motifs_up_enriched$Rank <- rank(-motifs_up_enriched$ORI_ave,ties.method= "min")
motifs_up_enriched <- motifs_up_enriched[order(motifs_up_enriched$Rank),]

motifs_depleted<- subset(motifs_up_enriched1, ORI_ave < 1/1.2 )
file_path=paste(args[2], args[3], sep = "/")

write.csv(motifs_up_enriched,file=paste(file_path,"_MotifEnriched_",args[1],"Resampling",".csv",sep = "",collapse =""))
#write.csv(motifs_depleted,file=paste(args[2],"DepletedMotifs",args[1],"Resampling",format(Sys.time(), "%Y-%m-%d_%H:%M"),".csv",sep = "",collapse =""))

write.csv(motifs_up,file=paste(file_path,"_MotifAll_",args[1],"Resampling",".csv",sep = "",collapse =""))
