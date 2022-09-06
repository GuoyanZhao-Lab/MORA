# Script to run MORA and the combined analysis from multiple tools
# 04/05/2021 rewrite
# Yizhe Song
if (!require('tidyverse')) install.packages('tidyverse')
devtools::install_github("tidyverse/tidyverse")


if (!require('ExperimentHub')) install.packages('ExperimentHub')

if (!require('universalmotif')) install.packages('universalmotif')
if (!require('lifecycle')) install.packages('lifecycle')
if (!require('xlsx')) install.packages('xlsx')
devtools::install_github('colearendt/xlsx')
install.packages("rJava",,"http://rforge.net")
#

# title: "TF predictions based on multiple tools"
# author: "Yizhe Song"
# date: "August 2020"

# Function: comparing the results from multiple tools to generate a ranked TF summary table
# Input file: 
# 1. MORA results in csv format
# 2.oPOSSUM-3 original downloadable results in txt format http://opossum.cisreg.ca/
# 3.Pscan original downloadable results in txt format http://159.149.160.88/pscan/
# 4.homer original results in txt format
# 5.sheetname of the spreadsheet file. eg "Sheet1", "Sheet2"
# 6.Output file name and path. eg "/home/ysong/data_for_motif_analysis/Eric/integration/Test_Integrated_TFBS_Predictions.xlsx"


## Function
RankTF <- function(MORA,op3,pscan,homer,arg5,arg6) {
  library(tidyverse)
  library(xlsx)
  library(dplyr)
  library(purrr)
  NewMasterTable <- readRDS("/home/ysong/MORA/data_for_motif_analysis/T2scKO/MasterTable2.Rds")
  if(missing(op3)&&!missing(pscan)&&!missing(homer)&&!missing(MORA)) {
    MORA <- read.csv(MORA)
    MORA <-as.data.frame(sapply(MORA, toupper),stringsAsFactors=FALSE)
    MORA <- MORA[,1, drop = FALSE]
    colnames(MORA)[1]<- "MORA"
    MORA <- dplyr::left_join(MORA,NewMasterTable, by = c("MORA"="Matrix_ID"))
    
    pscan <- read.delim(pscan)
    pscan <- filter(pscan, pscan$P_VALUE<0.05)
    pscan <-pscan[,2]
    pscan <-as.data.frame(sapply(pscan, toupper),stringsAsFactors=FALSE)
    colnames(pscan)[1]<- "Pscan"
    pscan <- dplyr::left_join(pscan,NewMasterTable, by = c("Pscan"="Matrix_ID"))
    
    homer <- read.delim(homer)
    homer <- homer[homer$p.value<=0.05,]
    homer<- tidyr::separate(homer,Motif.Name,into=c('TF_Name','Alias'), sep = "/")
    homer<-tidyr::separate_rows(homer,TF_Name, sep = "([\\,\\|\\:\\(\\)])", convert = FALSE)
    homer <- homer[homer$TF_Name!="",]
    homer <- dplyr::distinct(homer,TF_Name, .keep_all = FALSE)
    homer <-as.data.frame(sapply(homer, toupper),stringsAsFactors=FALSE)
    colnames(homer)[1]<- "homer"
    HOMER <- dplyr::left_join(homer,NewMasterTable, by = c("homer"="Input_Name"))
    
    new <-list(HOMER[,c(1,5)],pscan[,c(1,5)],MORA[,c(1,5)]) %>% 
      reduce(full_join, by = "ID")
    new <- dplyr::left_join(new,NewMasterTable[,c(4:7)], by="ID",ignore_case = TRUE)
    
    new$Pscan[!is.na(new$Pscan)] <- "1"
    new$Pscan[is.na(new$Pscan)] <- "0"
    new$Pscan <- as.numeric(as.character(new$Pscan))
    
    new$homer[!is.na(new$homer)] <- "1"
    new$homer[is.na(new$homer)] <- "0"
    new$homer <- as.numeric(as.character(new$homer))
    
    new$MORA[!is.na(new$MORA)] <- "1"
    new$MORA[is.na(new$MORA)] <- "0"
    new$MORA <- as.numeric(as.character(new$MORA))
    new2 <- dplyr::group_by(new,CIS_BP_Name) %>%
      mutate(Pscan=max(Pscan),homer=max(homer),MORA=max(MORA))
    
    new2 <-dplyr::distinct(new2,CIS_BP_Name,ID, .keep_all = TRUE)
    new2 <- dplyr::group_by(new2,CIS_BP_Name) %>%
      mutate(ID  = paste(ID, collapse =","),TF_Species  = paste(TF_Species, collapse =","))
    new2 <-dplyr::distinct(new2,CIS_BP_Name,.keep_all = TRUE)
    
    new2 <- new2 %>% 
      #rowwise will make sure the sum operation will occur on each row
      rowwise() %>% 
      #then a simple sum(..., na.rm=TRUE) is enough to result in what you need
      mutate(Num_Yes = sum(homer,Pscan,MORA, na.rm=TRUE))
    
    new2 <- dplyr::arrange(new2,desc(Num_Yes))
    new2 <-dplyr::distinct(new2,CIS_BP_Name, .keep_all = TRUE)
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    new2$Rank <- rank(-new2$Num_Yes,ties.method= "min")
    
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    new2 <- new2 %>% rename(TF_Name = CIS_BP_Name)
    #colnames(new2)[6] <-"TF_Name"
    new2 <-new2 %>% relocate(ID, .after = last_col())
    
    write.xlsx(new2, file=arg6, sheetName=arg5,append=TRUE, row.names=TRUE)
  } else if(missing(pscan)&&!missing(op3)&&!missing(homer)&&!missing(MORA)) {
    MORA <- read.csv(MORA)
    MORA <-as.data.frame(sapply(MORA, toupper),stringsAsFactors=FALSE)
    MORA <- MORA[,1, drop = FALSE]
    colnames(MORA)[1]<- "MORA"
    MORA <- dplyr::left_join(MORA,NewMasterTable, by = c("MORA"="Matrix_ID"))
    
    op3 <- read.delim(op3)
    op3 <- filter(op3, op3$Z.score>=3 | op3$Fisher.score>=5)
    op3 <-op3[,2]
    op3 <-as.data.frame(sapply(op3, toupper),stringsAsFactors=FALSE)
    colnames(op3)[1]<- "oPOSSUM"
    op3 <- dplyr::left_join(op3,NewMasterTable, by = c("oPOSSUM"="Matrix_ID"))
    
    homer <- read.delim(homer)
    homer <- homer[homer$p.value<=0.05,]
    homer<- tidyr::separate(homer,Motif.Name,into=c('TF_Name','Alias'), sep = "/")
    homer<-tidyr::separate_rows(homer,TF_Name, sep = "([\\,\\|\\:\\(\\)])", convert = FALSE)
    homer <- homer[homer$TF_Name!="",]
    homer <- dplyr::distinct(homer,TF_Name, .keep_all = FALSE)
    homer <-as.data.frame(sapply(homer, toupper),stringsAsFactors=FALSE)
    colnames(homer)[1]<- "homer"
    HOMER <- dplyr::left_join(homer,NewMasterTable, by = c("homer"="Input_Name"))
    
    new <-list(MORA[,c(1,5)],HOMER[,c(1,5)],op3[,c(1,5)]) %>% 
      reduce(full_join, by = "ID")
    new <- dplyr::left_join(new,NewMasterTable[,c(4:7)], by="ID",ignore_case = TRUE)
    
    new$oPOSSUM[!is.na(new$oPOSSUM)] <- "1"
    new$oPOSSUM[is.na(new$oPOSSUM)] <- "0"
    new$oPOSSUM <- as.numeric(as.character(new$oPOSSUM))
    
    
    new$homer[!is.na(new$homer)] <- "1"
    new$homer[is.na(new$homer)] <- "0"
    new$homer <- as.numeric(as.character(new$homer))
    
    new$MORA[!is.na(new$MORA)] <- "1"
    new$MORA[is.na(new$MORA)] <- "0"
    new$MORA <- as.numeric(as.character(new$MORA))
    new2 <- dplyr::group_by(new,CIS_BP_Name) %>%
      mutate(homer=max(homer),oPOSSUM=max(oPOSSUM),MORA=max(MORA))
    
    new2 <-dplyr::distinct(new2,CIS_BP_Name,ID, .keep_all = TRUE)
    new2 <- dplyr::group_by(new2,CIS_BP_Name) %>%
      mutate(ID  = paste(ID, collapse =","),TF_Species  = paste(TF_Species, collapse =","))
    new2 <-dplyr::distinct(new2,CIS_BP_Name,.keep_all = TRUE)
    
    new2 <- new2 %>% 
      #rowwise will make sure the sum operation will occur on each row
      rowwise() %>% 
      #then a simple sum(..., na.rm=TRUE) is enough to result in what you need
      mutate(Num_Yes = sum(oPOSSUM,homer,MORA, na.rm=TRUE))
    
    new2 <- dplyr::arrange(new2,desc(Num_Yes))
    new2 <-dplyr::distinct(new2,CIS_BP_Name, .keep_all = TRUE)
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    new2$Rank <- rank(-new2$Num_Yes,ties.method= "min")
    
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    
    new2 <-new2 %>% relocate(ID, .after = last_col())
    new2 <- new2 %>% rename(TF_Name = CIS_BP_Name)
    
    write.xlsx(new2, file=arg6, sheetName=arg5,append=TRUE, row.names=TRUE)
  } else if(missing(homer)&&!missing(pscan)&&!missing(op3)&&!missing(MORA)) {
    MORA <- read.csv(MORA)
    MORA <-as.data.frame(sapply(MORA, toupper),stringsAsFactors=FALSE)
    MORA <- MORA[,1, drop = FALSE]
    colnames(MORA)[1]<- "MORA"
    MORA <- dplyr::left_join(MORA,NewMasterTable, by = c("MORA"="Matrix_ID"))
    
    op3 <- read.delim(op3)
    op3 <- filter(op3, op3$Z.score>=3 | op3$Fisher.score>=5)
    op3 <-op3[,2]
    op3 <-as.data.frame(sapply(op3, toupper),stringsAsFactors=FALSE)
    colnames(op3)[1]<- "oPOSSUM"
    op3 <- dplyr::left_join(op3,NewMasterTable, by = c("oPOSSUM"="Matrix_ID"))
    
    pscan <- read.delim(pscan)
    pscan <- filter(pscan, pscan$P_VALUE<0.05)
    pscan <-pscan[,2]
    pscan <-as.data.frame(sapply(pscan, toupper),stringsAsFactors=FALSE)
    colnames(pscan)[1]<- "Pscan"
    pscan <- dplyr::left_join(pscan,NewMasterTable, by = c("Pscan"="Matrix_ID"))
    
    new <-list(MORA[,c(1,5)], pscan[,c(1,5)],op3[,c(1,5)]) %>% 
      reduce(full_join, by = "ID")
    new <- dplyr::left_join(new,NewMasterTable[,c(4:7)], by="ID",ignore_case = TRUE)
    
    new$oPOSSUM[!is.na(new$oPOSSUM)] <- "1"
    new$oPOSSUM[is.na(new$oPOSSUM)] <- "0"
    new$oPOSSUM <- as.numeric(as.character(new$oPOSSUM))
    
    new$Pscan[!is.na(new$Pscan)] <- "1"
    new$Pscan[is.na(new$Pscan)] <- "0"
    new$Pscan <- as.numeric(as.character(new$Pscan))
    
    new$MORA[!is.na(new$MORA)] <- "1"
    new$MORA[is.na(new$MORA)] <- "0"
    new$MORA <- as.numeric(as.character(new$MORA))
    new2 <- dplyr::group_by(new,CIS_BP_Name) %>%
      mutate(Pscan=max(Pscan),oPOSSUM=max(oPOSSUM),MORA=max(MORA))
    
    new2 <-dplyr::distinct(new2,CIS_BP_Name,ID, .keep_all = TRUE)
    new2 <- dplyr::group_by(new2,CIS_BP_Name) %>%
      mutate(ID  = paste(ID, collapse =","),TF_Species  = paste(TF_Species, collapse =","))
    new2 <-dplyr::distinct(new2,CIS_BP_Name,.keep_all = TRUE)
    
    new2 <- new2 %>% 
      #rowwise will make sure the sum operation will occur on each row
      rowwise() %>% 
      #then a simple sum(..., na.rm=TRUE) is enough to result in what you need
      mutate(Num_Yes = sum(oPOSSUM,Pscan,MORA, na.rm=TRUE))
    
    new2 <- dplyr::arrange(new2,desc(Num_Yes))
    new2 <-dplyr::distinct(new2,CIS_BP_Name, .keep_all = TRUE)
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    new2$Rank <- rank(-new2$Num_Yes,ties.method= "min")
    
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    
    new2 <-new2 %>% relocate(ID, .after = last_col())
    new2 <- new2 %>% rename(TF_Name = CIS_BP_Name)
    write.xlsx(new2, file=arg6, sheetName=arg5,append=TRUE, row.names=TRUE)
  } else if(missing(op3)&& missing(pscan) &&!missing(MORA)&&!missing(homer)) {
    MORA <- read.csv(MORA)
    MORA <-as.data.frame(sapply(MORA, toupper),stringsAsFactors=FALSE)
    MORA <- MORA[,1, drop = FALSE]
    colnames(MORA)[1]<- "MORA"
    MORA <- dplyr::left_join(MORA,NewMasterTable, by = c("MORA"="Matrix_ID"))
    
    homer <- read.delim(homer)
    homer <- homer[homer$p.value<=0.05,]
    homer<- tidyr::separate(homer,Motif.Name,into=c('TF_Name','Alias'), sep = "/")
    homer<-tidyr::separate_rows(homer,TF_Name, sep = "([\\,\\|\\:\\(\\)])", convert = FALSE)
    homer <- homer[homer$TF_Name!="",]
    homer <- dplyr::distinct(homer,TF_Name, .keep_all = FALSE)
    homer <-as.data.frame(sapply(homer, toupper),stringsAsFactors=FALSE)
    colnames(homer)[1]<- "homer"
    HOMER <- dplyr::left_join(homer,NewMasterTable, by = c("homer"="Input_Name"))
    
    new <-list(MORA[,c(1,5)], HOMER[,c(1,5)]) %>% 
      reduce(full_join, by = "ID")
    new <- dplyr::left_join(new,NewMasterTable[,c(4:7)], by="ID",ignore_case = TRUE)
    
    new$homer[!is.na(new$homer)] <- "1"
    new$homer[is.na(new$homer)] <- "0"
    new$homer <- as.numeric(as.character(new$homer))
    
    new$MORA[!is.na(new$MORA)] <- "1"
    new$MORA[is.na(new$MORA)] <- "0"
    new$MORA <- as.numeric(as.character(new$MORA))
    new2 <- dplyr::group_by(new,CIS_BP_Name) %>%
      mutate(homer=max(homer),MORA=max(MORA))
    new2 <-dplyr::distinct(new2,CIS_BP_Name,ID, .keep_all = TRUE)
    new2 <- dplyr::group_by(new2,CIS_BP_Name) %>%
      mutate(ID  = paste(ID, collapse =","),TF_Species  = paste(TF_Species, collapse =","))
    new2 <-dplyr::distinct(new2,CIS_BP_Name,.keep_all = TRUE)
    
    new2 <- new2 %>% 
      #rowwise will make sure the sum operation will occur on each row
      rowwise() %>% 
      #then a simple sum(..., na.rm=TRUE) is enough to result in what you need
      mutate(Num_Yes = sum(homer,MORA, na.rm=TRUE))
    
    new2 <- dplyr::arrange(new2,desc(Num_Yes))
    new2 <-dplyr::distinct(new2,CIS_BP_Name, .keep_all = TRUE)
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    new2$Rank <- rank(-new2$Num_Yes,ties.method= "min")
    
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    
    new2 <-new2 %>% relocate(ID, .after = last_col())
    new2 <- new2 %>% rename(TF_Name = CIS_BP_Name)
    write.xlsx(new2, file=arg6, sheetName=arg5,append=TRUE, row.names=TRUE)   
  } else if (missing(op3) && missing(homer)&&!missing(MORA)&&!missing(pscan)){
    MORA <- read.csv(MORA)
    MORA <-as.data.frame(sapply(MORA, toupper),stringsAsFactors=FALSE)
    MORA <- MORA[,1, drop = FALSE]
    colnames(MORA)[1]<- "MORA"
    MORA <- dplyr::left_join(MORA,NewMasterTable, by = c("MORA"="Matrix_ID"))  
    pscan <- read.delim(pscan)
    pscan <- filter(pscan, pscan$P_VALUE<0.05)
    pscan <-pscan[,2]
    pscan <-as.data.frame(sapply(pscan, toupper),stringsAsFactors=FALSE)
    colnames(pscan)[1]<- "Pscan"
    pscan <- dplyr::left_join(pscan,NewMasterTable, by = c("Pscan"="Matrix_ID"))
    
    new <-list(MORA[,c(1,5)], pscan[,c(1,5)]) %>% 
      reduce(full_join, by = "ID")
    new <- dplyr::left_join(new,NewMasterTable[,c(4:7)], by="ID",ignore_case = TRUE)
    
    new$Pscan[!is.na(new$Pscan)] <- "1"
    new$Pscan[is.na(new$Pscan)] <- "0"
    new$Pscan <- as.numeric(as.character(new$Pscan))
    
    new$MORA[!is.na(new$MORA)] <- "1"
    new$MORA[is.na(new$MORA)] <- "0"
    new$MORA <- as.numeric(as.character(new$MORA))
    new2 <- dplyr::group_by(new,CIS_BP_Name) %>%
      mutate(Pscan=max(Pscan),MORA=max(MORA))
    new2 <-dplyr::distinct(new2,CIS_BP_Name,ID, .keep_all = TRUE)
    new2 <- dplyr::group_by(new2,CIS_BP_Name) %>%
      mutate(ID  = paste(ID, collapse =","),TF_Species  = paste(TF_Species, collapse =","))
    new2 <-dplyr::distinct(new2,CIS_BP_Name,.keep_all = TRUE)
    new2 <- new2 %>% 
      #rowwise will make sure the sum operation will occur on each row
      rowwise() %>% 
      #then a simple sum(..., na.rm=TRUE) is enough to result in what you need
      mutate(Num_Yes = sum(Pscan,MORA, na.rm=TRUE))
    
    new2 <- dplyr::arrange(new2,desc(Num_Yes))
    new2 <-dplyr::distinct(new2,CIS_BP_Name, .keep_all = TRUE)
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    new2$Rank <- rank(-new2$Num_Yes,ties.method= "min")
    
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    
    new2 <-new2 %>% relocate(ID, .after = last_col())
    new2 <- new2 %>% rename(TF_Name = CIS_BP_Name)
    write.xlsx(new2, file=arg6, sheetName=arg5,append=TRUE, row.names=TRUE)    
  } else if (missing(pscan)&& missing(homer)&&!missing(MORA)&&!missing(op3)) {
    MORA <- read.csv(MORA)
    MORA <-as.data.frame(sapply(MORA, toupper),stringsAsFactors=FALSE)
    MORA <- MORA[,1, drop = FALSE]
    colnames(MORA)[1]<- "MORA"
    MORA <- dplyr::left_join(MORA,NewMasterTable, by = c("MORA"="Matrix_ID"))
    
    op3 <- read.delim(op3)
    op3 <- filter(op3, op3$Z.score>=3 | op3$Fisher.score>=5)
    op3 <-op3[,2]
    op3 <-as.data.frame(sapply(op3, toupper),stringsAsFactors=FALSE)
    colnames(op3)[1]<- "oPOSSUM"
    op3 <- dplyr::left_join(op3,NewMasterTable, by = c("oPOSSUM"="Matrix_ID"))
    
    new <-list(MORA[,c(1,5)], op3[,c(1,5)]) %>% 
      reduce(full_join, by = "ID")
    new <- dplyr::left_join(new,NewMasterTable[,c(4:7)], by="ID",ignore_case = TRUE)
    
    new$oPOSSUM[!is.na(new$oPOSSUM)] <- "1"
    new$oPOSSUM[is.na(new$oPOSSUM)] <- "0"
    new$oPOSSUM <- as.numeric(as.character(new$oPOSSUM))
    
    new$MORA[!is.na(new$MORA)] <- "1"
    new$MORA[is.na(new$MORA)] <- "0"
    new$MORA <- as.numeric(as.character(new$MORA))
    new2 <- dplyr::group_by(new,CIS_BP_Name) %>%
      mutate(oPOSSUM=max(oPOSSUM),MORA=max(MORA))
    new2 <-dplyr::distinct(new2,CIS_BP_Name,ID, .keep_all = TRUE)
    new2 <- dplyr::group_by(new2,CIS_BP_Name) %>%
      mutate(ID  = paste(ID, collapse =","),TF_Species  = paste(TF_Species, collapse =","))
    new2 <-dplyr::distinct(new2,CIS_BP_Name,.keep_all = TRUE)
    
    new2 <- new2 %>% 
      #rowwise will make sure the sum operation will occur on each row
      rowwise() %>% 
      #then a simple sum(..., na.rm=TRUE) is enough to result in what you need
      mutate(Num_Yes = sum(oPOSSUM,MORA, na.rm=TRUE))
    
    new2 <- dplyr::arrange(new2,desc(Num_Yes))
    new2 <-dplyr::distinct(new2,CIS_BP_Name, .keep_all = TRUE)
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    new2$Rank <- rank(-new2$Num_Yes,ties.method= "min")
    
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    
    new2 <-new2 %>% relocate(ID, .after = last_col())
    new2 <- new2 %>% rename(TF_Name = CIS_BP_Name)
    write.xlsx(new2, file=arg6, sheetName=arg5,append=TRUE, row.names=TRUE)
  }else {
    MORA <- read.csv(MORA)
    MORA <-as.data.frame(sapply(MORA, toupper),stringsAsFactors=FALSE)
    MORA <- MORA[,1, drop = FALSE]
    colnames(MORA)[1]<- "MORA"
    MORA <- dplyr::left_join(MORA,NewMasterTable, by = c("MORA"="Matrix_ID"))
    
    op3 <- read.delim(op3)
    op3 <- filter(op3, op3$Z.score>=3 | op3$Fisher.score>=5)
    op3 <-op3[,2]
    op3 <-as.data.frame(sapply(op3, toupper),stringsAsFactors=FALSE)
    colnames(op3)[1]<- "oPOSSUM"
    op3 <- dplyr::left_join(op3,NewMasterTable, by = c("oPOSSUM"="Matrix_ID"))
    
    pscan <- read.delim(pscan)
    pscan <- filter(pscan, pscan$P_VALUE<0.05)
    pscan <-pscan[,2]
    pscan <-as.data.frame(sapply(pscan, toupper),stringsAsFactors=FALSE)
    colnames(pscan)[1]<- "Pscan"
    pscan <- dplyr::left_join(pscan,NewMasterTable, by = c("Pscan"="Matrix_ID"))
    
    homer <- read.delim(homer)
    homer <- homer[homer$p.value<=0.05,]
    homer<- tidyr::separate(homer,Motif.Name,into=c('TF_Name','Alias'), sep = "/")
    homer<-tidyr::separate_rows(homer,TF_Name, sep = "([\\,\\|\\:\\(\\)])", convert = FALSE)
    homer <- homer[homer$TF_Name!="",]
    homer <- dplyr::distinct(homer,TF_Name, .keep_all = FALSE)
    homer <-as.data.frame(sapply(homer, toupper),stringsAsFactors=FALSE)
    colnames(homer)[1]<- "homer"
    HOMER <- dplyr::left_join(homer,NewMasterTable, by = c("homer"="Input_Name"))
    
    new <-list(MORA[,c(1,5)],HOMER[,c(1,5)],pscan[,c(1,5)],op3[,c(1,5)]) %>% 
      reduce(full_join, by = "ID")
    new <- dplyr::left_join(new,NewMasterTable[,c(4:7)], by="ID",ignore_case = TRUE)
    
    new$oPOSSUM[!is.na(new$oPOSSUM)] <- "1"
    new$oPOSSUM[is.na(new$oPOSSUM)] <- "0"
    new$oPOSSUM <- as.numeric(as.character(new$oPOSSUM))
    
    new$Pscan[!is.na(new$Pscan)] <- "1"
    new$Pscan[is.na(new$Pscan)] <- "0"
    new$Pscan <- as.numeric(as.character(new$Pscan))
    
    new$homer[!is.na(new$homer)] <- "1"
    new$homer[is.na(new$homer)] <- "0"
    new$homer <- as.numeric(as.character(new$homer))
    
    new$MORA[!is.na(new$MORA)] <- "1"
    new$MORA[is.na(new$MORA)] <- "0"
    new$MORA <- as.numeric(as.character(new$MORA))
    
    new2 <- dplyr::group_by(new,CIS_BP_Name) %>%
      mutate(Pscan=max(Pscan),homer=max(homer),oPOSSUM=max(oPOSSUM),MORA=max(MORA))
    
    
    new2 <-dplyr::distinct(new2,CIS_BP_Name,ID, .keep_all = TRUE)
    new2 <- dplyr::group_by(new2,CIS_BP_Name) %>%
      mutate(ID  = paste(ID, collapse =","),TF_Species  = paste(TF_Species, collapse =","))
    new2 <-dplyr::distinct(new2,CIS_BP_Name,.keep_all = TRUE)
    
    new2 <- new2 %>% 
      #rowwise will make sure the sum operation will occur on each row
      rowwise() %>% 
      #then a simple sum(..., na.rm=TRUE) is enough to result in what you need
      mutate(Num_Yes = sum(oPOSSUM,homer,Pscan,MORA, na.rm=TRUE))
    
    new2 <- dplyr::arrange(new2,desc(Num_Yes))
    new2 <-dplyr::distinct(new2,CIS_BP_Name, .keep_all = TRUE)
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    new2$Rank <- rank(-new2$Num_Yes,ties.method= "min")
    
    new2 <- new2[!is.na(new2$CIS_BP_Name),]
    
    new2 <-new2 %>% relocate(ID, .after = last_col())
    new2 <- new2 %>% rename(TF_Name = CIS_BP_Name)
    write.xlsx(new2, file=arg6, sheetName=arg5,append=TRUE, row.names=TRUE)
  }
}

## Example code to execute commands

MORA <-"/home/ysong/data_for_motif_analysis/Eric/integration/MORA/1.csv"
op3 <-"/home/ysong/data_for_motif_analysis/Eric/integration/op3/1.txt"
pscan<-"/home/ysong/data_for_motif_analysis/Eric/integration/pscan/1.txt"
homer <-"/home/ysong/data_for_motif_analysis/Eric/integration/Homer/1.txt"
#arg5 <-"sheet1"
arg6 <-"/home/ysong/data_for_motif_analysis/Test_Integrated_TFBS_Predictions04062021.xlsx"

### Run function with missing arguments
RankTF(MORA = MORA,op3=op3,pscan=pscan, homer=homer,arg5="Sheet1",arg6=arg6)
RankTF(MORA = MORA,pscan=pscan, homer=homer,arg5="Sheet2",arg6=arg6)
RankTF(MORA = MORA,op3=op3,homer=homer,arg5="Sheet3",arg6=arg6)
RankTF(MORA = MORA,op3=op3,pscan=pscan,arg5="Sheet4",arg6=arg6)
RankTF(MORA = MORA,homer=homer,arg5="Sheet5",arg6=arg6)
RankTF(MORA = MORA,pscan=pscan,arg5="Sheet6",arg6=arg6)
RankTF(MORA = MORA,op3=op3,arg5="Sheet7",arg6=arg6)