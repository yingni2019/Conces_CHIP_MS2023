######
# Investigate CHIP in solid tumor summary script 
######
#Author: Ying Ni
#Date: June, 2021
######


library(ggrepel) # ggplot2 add on for better labels
library(dplyr) 
#library(readr)
library(readxl)
library(ggplot2)
library(survival)
library(survminer)
library(arm)
library(My.stepwise)
library(maxstat)
library(tidyverse)
library(stringr)
library(readr)
library(data.table)
library(g3viz)

#######################################################
### All mutation calls from Foundation Medicine and Tempus test results are aggregated and annotated with ANNOVAR, with CLINVAR20210123 and COSMIC70 annotation
###
###
### Download mutation list text files before running this script
#######################################################



All3Sets_VAF35<-read.csv("All3Sets_VAF35_COSMIC.csv", header=TRUE)
df<-All3Sets_VAF35
Set<-"breast" #"Lung" # "breast","CRC","Lung"
df<-All3Sets_lung[which(All3Sets_lung$VAF_in_Tumor < 0.35),]
table(df$Gene,df$Set,df$TumorType)

AllGene.Table<-data.frame(rbind(table(df$Gene,df$Set)))
AllGene.Table2<-data.frame(rbind(table(df$Gene,df$TumorType)))
ROWnames<-rownames(AllGene.Table)
write.csv(AllGene.Table, file="All3Sets_AllGeneTable_bySet.csv")
write.csv(AllGene.Table2, file="All3Sets_AllGeneTable_byTumotType.csv")

df1 <- cbind(AllGene.Table, t(apply(AllGene.Table, 1, function(x) {
  ch <- chisq.test(rbind(x,(AllCount[Set,]-x)))
  c(unname(ch$statistic), ch$p.value)})))
colnames(df1)[4:5] <- c('x-squared', 'p-value')

df2 <- cbind(df1, t(apply(AllGene.Table[,c(2:3)], 1, function(x) {
  fi <- fisher.test(rbind(x,(AllCount[Set,c(2:3)]-x)))
  c(unname(fi$estimate), fi$p.value)})))
colnames(df2)[6:7] <- c('estimate', 'p-value between FM and Tempus')
write.csv(df2,file=paste0(Set,"_VAF35_Xtest_fishers.csv"))

AllCancer.Table<-table(All3Sets_VAF35$Gene,All3Sets_VAF35$TumorType)
write.csv(AllCancer.Table, file="All3Cancers_AllGeneTable.csv")

for (L in c("Tempus","FM","Coombs")) {
  df<-All3Sets[which(All3Sets$VAF_in_Tumor < 0.35 & All3Sets$Set == L),]
  Cancer.Table<-table(df$Gene,df$TumorType)
  print(paste0("For AllCancers:",L, "Positive N = ", length(unique(df$SampleID))))
}




### Making g3lollipop input file is such a pain...
All3Sets_VAF35.mutation.txt<-"All3Sets_VAF35_COSMIC.csv"

mutation.dat <- readMAF(All3Sets_VAF35.mutation.txt<-"All3Sets_VAF35_COSMIC.csv",
                               gene.symbol.col = "Gene",
                               variant.class.col = "Mutation_Type",
                               protein.change.col = "AA_change",
                               sep = ",")  # column-separator of csv file

Gene_to_plot<-c("TET2", "ASXL1",  "TP53", "SF3B1", "JAK2",  "DNMT3A")
Set_to_plot<-c("Tempus","FM","Coombs")

    Title<-paste0("FM and Tempus_",G)
    Y_Axis_Title<-paste0(G," Mutation Count")
    Plot_data<-mutation.dat[which(mutation.dat$Gene == G &  mutation.dat$Set != "Coombs"),]
    output.filename<-paste0("FM and Tempus_VAF35_",G)
#chart <-
    g3Lollipop(Plot_data,
               gene.symbol.col = "Gene", #"Gene",
               gene.symbol = G,
               factor.col = "Set", #factor.col = "Mutation_Type",
               aa.pos.col = "AA_Position",
               protein.change.col = "AA_change", # "AA_change",
               btn.style = "blue", # blue-style chart download buttons
               plot.options =
                 g3Lollipop.theme(theme.name = "cbioportal",
                                  # Chart title settings
                                  title.text = Title,
                                  y.axis.label = Y_Axis_Title,
                                  legend.title = "Set"),
               output.filename = output.filename)



Gene_to_plot<-c("TET2", "ASXL1",  "TP53", "SF3B1", "JAK2",  "DNMT3A")
Set_to_plot<-c("Tempus","FM","Coombs")
for (G in Gene_to_plot){
  #for (S in Set_to_plot){
    Title<-paste0(S, "_",G)
    Y_Axis_Title<-paste0(G," Mutation Count")
    Plot_data<-mutation.dat[which(mutation.dat$Gene == G & mutation.dat$Set == S),]
    output.filename<-paste0(S, "_COSMIC_VAF35_",G,".png")
    #chart <-
    g3Lollipop(Plot_data,
               gene.symbol.col = "Gene", #"Gene",
               gene.symbol = G,
               #factor.col = "Set", 
               factor.col = "Mutation_Type",
               aa.pos.col = "AA_Position",
               protein.change.col = "AA_change", # "AA_change",
               btn.style = "blue", # blue-style chart download buttons
               plot.options =
                 g3Lollipop.theme(theme.name = "cbioportal",
                                  # Chart title settings
                                  title.text = Title,
                                  y.axis.label = Y_Axis_Title,
                                  legend.title = "Mutation_Type"
                                  ),
               output.filename = output.filename)
 # }
}





