#*****************************************************************************************
# Rscript - Differential Expression Signature Generation
#           
#           This code will generate the differential expression signatures used
#           in this study.  The signatures generated will be for each of the LP-
#           signatures (up and down regulated genes).
#
#           The resulting signatures are stored into a single dataframe that
#           is used throughout the remaining analyses.
#*****************************************************************************************

#***********************************************
#* Load Libraries
#***********************************************
source("https://bioconductor.org/biocLite.R")
library("qvalue")


#***********************************************
#* Load Transcriptomics data
#***********************************************
load(paste(getwd(), "/Data/ADEPT_asthma.RNA_ebbx.RData", sep = ""))
#Get one file called rna.ebbx.expr  

asthma.Phenotype.data <- read.delim(paste(getwd(), "/Data/ADEPT_asthma.Phenotype.data.txt", sep = ""), stringsAsFactors=FALSE, na.strings = c("NA", "N/A"))

asthma.Phenotype.desc <- read.delim(paste(getwd(), "/Data/ADEPT_asthma.Phenotype.desc.txt", sep = ""), stringsAsFactors=FALSE, na.strings = c("NA", "N/A"))

#***********************************************
#* Gene Symbol to Entrez ID link
#***********************************************
library(org.Hs.eg.db)
sym2entrez <- toTable(org.Hs.egSYMBOL2EG)

#***********************************************
#* Load Functions
#***********************************************
source(paste(getwd(), "/WIP_Code/ADEPT_Functions.R", sep = ""))

#***********************************************
#* Run ANOVA Analysis
#***********************************************
rnaFile <- get("rna.ebbx.expr")

rnaFile <- rnaFile[-grep("///", row.names(rnaFile)), ]

anovaRNATable <- t(as.data.frame(rnaFile))
anovaRNATable <- as.data.frame(anovaRNATable)

anovaRNACov <- data.frame(RND_SUBJID=row.names(anovaRNATable))

anovaRNACov$SEX <- as.factor(asthma.Phenotype.data$SEX[match(anovaRNACov$RND_SUBJID, asthma.Phenotype.data$RND_SUBJID)])
anovaRNACov$SITE_ID <- as.factor(asthma.Phenotype.data$SITEID[match(anovaRNACov$RND_SUBJID, asthma.Phenotype.data$RND_SUBJID)])
anovaRNACov$AGE <- asthma.Phenotype.data$AGE[match(anovaRNACov$RND_SUBJID, asthma.Phenotype.data$RND_SUBJID)]

anovaRNACov$ARM <- factor(asthma.Phenotype.data$ARM[match(anovaRNACov$RND_SUBJID, asthma.Phenotype.data$RND_SUBJID)], levels = c("Healthy", "Mild", "Moderate", "Severe"))
#anovaRNACov$origIL13Cat <- paste("Th2", asthma.Phenotype.data$IL13_BrEpi.ES[match(anovaRNACov$RND_SUBJID, asthma.Phenotype.data$RND_SUBJID)], sep = "_")
#anovaRNACov$Th2 <- as.factor(ifelse(anovaRNACov$origIL13Cat == "Th2_Healthy", "Th2_Low", anovaRNACov$origIL13Cat))

SEX <- anovaRNACov$SEX
SITEID <- anovaRNACov$SITE_ID
AGE <- anovaRNACov$AGE

factorList <- c("ARM")#, "Th2", "strRESPONDER")
flagList <- c("LP")#, "Th2", "strResponse")

for(factorIndex in 1:length(factorList)){
  flag <- flagList[factorIndex]
  print(paste(flag, "running...", sep = " "))
  
  # if(flag == "strResponse"){ #Need to remove filter patients who are not solely SNR or SR
  #   anovaRNACov$origStrRESPONDER <- asthma.Phenotype.data$Molecular_Steroid_Response_cluster[match(anovaRNACov$RND_SUBJID, asthma.Phenotype.data$RND_SUBJID)]
  #   anovaRNACov$strRESPONDER <- sapply(anovaRNACov$origStrRESPONDER, function(x) strsplit(x, "[ (]")[[1]][1])
  #   anovaRNACov$strRESPONDER <- ifelse(anovaRNACov$strRESPONDER == "Mild", "NR", anovaRNACov$strRESPONDER)
  #   anovaRNACov <- anovaRNACov[-which(anovaRNACov$strRESPONDER == "Mix"), ]
  #   anovaRNACov$strRESPONDER <- ifelse(anovaRNACov$strRESPONDER == "NR", "SNR", ifelse(anovaRNACov$strRESPONDER == "R", "SR", "Healthy"))
  #   anovaRNACov$strRESPONDER <- as.factor(anovaRNACov$strRESPONDER)
  #   
  #   anovaRNATable <- anovaRNATable[match(anovaRNACov$RND_SUBJID, row.names(anovaRNATable)), ]
  #   
  #   SEX <- anovaRNACov$SEX
  #   SITEID <- anovaRNACov$SITE_ID
  #   AGE <- anovaRNACov$AGE
  # }
  
  filteredSamples <- anovaRNACov$RND_SUBJID
  FACTOR <- factorList[factorIndex]
  FACTOR <- anovaRNACov[, which(colnames(anovaRNACov) == FACTOR)]
  
  geneSummary <- anovaFunction(anovaRNACov, FACTOR, anovaRNATable)
  
  length(which(geneSummary$OverallTraitvsFactorPval < 0.05))
  geneSummary$OverallFDR <- p.adjust(geneSummary$OverallTraitvsFactorPval, method = "BH")
  
  assign(paste(flag, "geneSummaryALL_ebbx", sep = "_"), geneSummary)
  
  temp_summarySigQGenesInterest <- qvalueTest(geneSummary, FACTOR)
  assign(paste(flag, "SigQ_ebbx", sep = "_"), temp_summarySigQGenesInterest)
  
  print(paste(flag, "Complete!", sep = " "))
}

#***********************************************
#* Merge signatures into one table
#***********************************************
for(flagIndex in 1:length(flagList)){
  flag_i <- flagList[flagIndex]
  SigQfile_i <- get(paste(flag_i, "SigQ_ebbx", sep = "_"))
  
  if(flag_i == "LP"){
    posMild <- grep("[*]", SigQfile_i$`Mild:Healthy_Pval/FC`)
    downregPosMild <- grep("[*] / -", SigQfile_i$`Mild:Healthy_Pval/FC`)
    downMild <- SigQfile_i[downregPosMild, ]
    upregPosMild <- posMild[match(posMild, posMild)]
    downRegGenes_Mild <- SigQfile_i$EntrezID[downregPosMild]
    upRegGenes_Mild <- SigQfile_i$EntrezID[upregPosMild]
    
    posMod <- grep("[*]", SigQfile_i$`Moderate:Healthy_Pval/FC`)
    downregPosMod <- grep("[*] / -", SigQfile_i$`Moderate:Healthy_Pval/FC`)
    downMod <- SigQfile_i[downregPosMod, ]
    upregPosMod <- posMod[-match(downregPosMod, posMod)]
    downRegGenes_Mod <- SigQfile_i$EntrezID[downregPosMod]
    upRegGenes_Mod <- SigQfile_i$EntrezID[upregPosMod]
    
    posSevere <- grep("[*]", SigQfile_i$`Severe:Healthy_Pval/FC`)
    downregPosSevere <- grep("[*] / -", SigQfile_i$`Severe:Healthy_Pval/FC`)
    downSevere <- SigQfile_i[downregPosSevere, ]
    upregPosSevere <- posSevere[-match(downregPosSevere, posSevere)]
    downRegGenes_Severe <- SigQfile_i$EntrezID[downregPosSevere]
    upRegGenes_Severe<- SigQfile_i$EntrezID[upregPosSevere]
    
    posMod_Mild <- grep("[*]", SigQfile_i$`Moderate:Mild_Pval/FC`)
    downregPosMod_Mild <- grep("[*] / -", SigQfile_i$`Moderate:Mild_Pval/FC`)
    downMild_Mlod <- SigQfile_i[downregPosM, ]
    upregPosMild <- posMild[match(posMild, posMild)]
    downRegGenes_Mild <- SigQfile_i$EntrezID[downregPosMild]
    upRegGenes_Mild <- SigQfile_i$EntrezID[upregPosMild]
    
    LPsig <- vector("list", 6)
    LPsig[[1]] <- upRegGenes_Mild
    LPsig[[2]] <- downRegGenes_Mild
    LPsig[[3]] <- upRegGenes_Mod
    LPsig[[4]] <- downRegGenes_Mod
    LPsig[[5]] <- upRegGenes_Severe
    LPsig[[6]] <- downRegGenes_Severe
    names(LPsig) <- c("upregGenes_Mild", "downregGenes_Mild", "upregGenes_Mod", "downregGenes_Mod", "upregGenes_Severe", "downregGenes_Severe")
  }
  
  if(flag_i == "Th2"){
    Th2HighRow <- grep("-", SigQfile_i$`Th2_Low:Th2_High_Pval/FC`) 
    Th2_Th2High <- SigQfile_i$EntrezID[Th2HighRow]
    Th2_Th2Low <- SigQfile_i$EntrezID[-Th2HighRow]
    
    Th2sig <- vector("list", 2)
    Th2sig[[1]] <- Th2_Th2Low
    Th2sig[[2]] <- Th2_Th2High
    names(Th2sig) <- c("Th2_Low", "Th2_High")
  }
  
  if(flag_i == "strResponse"){
    SNR <- SigQfile_i$EntrezID[grep("[*] / -", SigQfile_i$`SR:SNR_Pval/FC`)]
    SR <- SigQfile_i$EntrezID[grep("[*] / [+]", SigQfile_i$`SR:SNR_Pval/FC`)]
    
    strResponsesig <- vector("list", 2)
    strResponsesig[[1]] <- SNR
    strResponsesig[[2]] <- SR
    names(strResponsesig) <- c("SNR", "SR")
  }
}

ALLsigs <- LPsig
ALLsigs[[7]] <- Th2sig[[1]]
ALLsigs[[8]] <- Th2sig[[2]]
ALLsigs[[9]] <- strResponsesig[[1]]
ALLsigs[[10]] <- strResponsesig[[2]]
names(ALLsigs) <- c("LP-Mild-up", "LP-Mild-down", "LP-Moderate-up", "LP-Moderate-down", "LP-Severe-up", "LP-Severe-down", "Th2_Low", "Th2_High", "SNR", "SR")

ALLsigsDF <- convertListofGenesTOdataframe(ALLsigs)

ALLsigs[[11]] <- row.names(rna.ebbx.expr)[-grep("///", row.names(rna.ebbx.expr))]
names(ALLsigs)[[11]] <- "bg_entrez"

save(ALLsigs, file = paste(getwd(), "/Data/Results/Allsigs.RData", sep = ""))
write.table(ALLsigsDF, file = paste(getwd(), "/Data/Results/TableS1_DEGeneSignaturesClinicalParameters.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

rm(temp_summarySigQGenesInterest, geneSummary,rnaFile, anovaRNATable,SEX, SITEID, AGE, FACTOR,flag_i, SigQfile_i, ALLsigs, LPsig, Th2sig, strResponsesig, posMild, downregPosMild, downMild, upregPosMild, downRegGenes_Mild, upRegGenes_Mild, posMod, downregPosMod, downMod, upregPosMod, downRegGenes_Mod, upRegGenes_Mod, posSevere, downregPosSevere, downSevere, upregPosSevere, downRegGenes_Severe, upRegGenes_Severe, Th2HighRow, SNR, SR)
