library(readr)
library(viper)
 
#/!\ /!\ /!\ /!\ /!\ /!\ /!\
#Change the working directory to were you cloned the repository on your machine
setwd("~/Dropbox/kinase_tf_mini_tuto/")
#/!\ /!\ /!\ /!\ /!\ /!\ /!\

source("script/viper_functions.R")

########## PHOSHO and KINASE part ########

library(OmnipathR)

#import the data
phospho_differential_analysis <- as.data.frame(
  read_csv("data/phospho_differential_analysis.csv"))

#format it properlly
row.names(phospho_differential_analysis) <- phospho_differential_analysis$psite_ID
phospho_differential_analysis <- phospho_differential_analysis[,-1, drop = F]

#inport KSN from omnipath
omnipath_ptm <- get_signed_ptms()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation","phosphorylation"),]
KSN <- omnipath_ptm[,c(4,3)]
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_type, sep ="_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol,omnipath_ptm$residue_offset, sep = "")
KSN$sign <- ifelse(omnipath_ptm$modification == "phosphorylation", 1, -1)

#Format KSN
KSN_viper <- df_to_viper_regulon(KSN)

#run viper to get the TF activities from the phosphoproteomic data
#You can also run that on wour normalised intesity matrix of phosphosites directly,
#as long as it is formatted as a dataframe of similar format as here
#User is strongly encouraged to check the viper publication (PMID: 27322546) for more info on the process
kin_activity <- as.data.frame(viper(eset = phospho_differential_analysis, regulon = KSN_viper, minsize = 5, adaptive.size = F, eset.filter = F))
kin_activity$ID <- row.names(kin_activity)

kin_activity <- kin_activity[,c(2,1)]
names(kin_activity) <- c("ID","NES")
########## RNA and TF part ########

library(dorothea)

#First we import the dorothea regulons (using only confidence A, B, and C), see dorothea publication for information on confidence levels
dorothea_df<- as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A","B","C"),c(3,1,4)])

#import the RNAseq data. It has entrez gene identifiers, but we need it to have gene symbols to match dorothea database, so we have
#to do some id conversion as well
RNA_differential_analysis <- as.data.frame(
  read_csv("data/RNA_differential_analysis.csv"))

#import a mapping table downloaded from uniprot
RNAseq_entrez_to_symbol <- as.data.frame(read_delim("support/RNAseq_entrez_to_symbol", 
                                                    "\t", escape_double = FALSE, col_types = cols(`yourlist:M20191127216DA2B77BFBD2E6699CA9B6D1C41EB259129CL` = col_character()), 
                                                    trim_ws = TRUE)) #from uniprot 20191127

names(RNAseq_entrez_to_symbol)[1] <- "ID"

#this part is to merge the mapping table with the differential analysis dataframe, using entrez gene id as a common key between the two
#this way, we will have our gene entrez identifiers mapped to their corresponding symbols in the differential analysis dataframe
#of course, there are many other way to achievethis goal.
RNA_differential_analysis <- merge(RNA_differential_analysis, RNAseq_entrez_to_symbol[,c(1,6)])
RNA_differential_analysis <- RNA_differential_analysis[,c(8,2:7)]
names(RNA_differential_analysis)[1] <- "ID"
RNA_differential_analysis$ID <- gsub(" .*","",RNA_differential_analysis$ID)
RNA_differential_analysis <- unique(RNA_differential_analysis)

#now we just need to format the differential analysis data into a format that is compatible with viper
eset <- RNA_differential_analysis$t
names(eset) <- RNA_differential_analysis$ID

#we also need to format the dorothea dataframe into a viper format
dorothea_viper <- df_to_viper_regulon(dorothea_df)

#Now we estimate the TF activities using viper
TF_activities <- as.data.frame(viper(eset = eset, regulon = dorothea_viper, minsize = 10, adaptive.size = F, eset.filter = F, pleiotropy = T))
TF_activities$TF <- row.names(TF_activities)

#that's just to make the dataframe pretty
TF_activities <- TF_activities[,c(2,1)]
names(TF_activities) <- c("ID","NES")


########## CONCLUSION

#Now you have succefully estimated kinase and TF activities from phosphoproteomic and transcriptomic
#You cna view the results here
View(kin_activity)
View(TF_activities)

#You can now combined them together and use them as input for COSMOS.
#You may also leave them separated and use them a separated input and measurments in cosmos, if you lack metabolomic data
#See https://github.com/saezlab/cosmosR for more info on how to use cosmos