install.packages("dplyr")
library(dplyr)
install.packages("tidyverse")
library(tidyverse)


#---------------------------------------------------------------
#Step1: 
#read the dbCAN annotation results of your genomes into the R environment
#---------------------------------------------------------------
#set the working directory as the directory that contains the dbCAN annotation results
#i.e., the folder with only the '*.out.dm.ps' files
#e.g., in my case, the absolute path is: D:\papers_in_progress\to finish\HKU_cellulolytic_genomes\dbcan_results
dbCAN_results_data_path = as.character("D:/papers_in_progress/to finish/HKU_cellulolytic_genomes/dbcan_results/")
setwd(dbCAN_results_data_path)
getwd()
cazy_in_genomes<-dir(dbCAN_results_data_path)
#cazy_in_genomes
#--------------------------------------------------------------
#-----preparation-----------------------------------------------------------------------------------------------------------------------
#*Categorization of the CAZy modules, accoridng to updated info on their functions
#*(this part could be manually updated/modified if the function description of CAZy module in the CAZy database updated)
#*CAZy database link: http://www.cazy.org/
#----------------------------------------------------------------------------------------------------------------------------
# "exo" refers exoglucanase GH modules
# "endo" refers endoglucanase GH modules

#Subfamilies
GH5<-c("GH5_1","GH5_2","GH5_4","GH5_5","GH5_7","GH5_8","GH5_9","GH5_10",
       "GH5_11","GH5_12","GH5_13","GH5_14","GH5_15","GH5_16","GH5_17","GH5_18","GH5_19",
       "GH5_20","GH5_21","GH5_22","GH5_23","GH5_24","GH5_25","GH5_26","GH5_27","GH5_28","GH5_29",
       "GH5_30","GH5_31","GH5_32","GH5_33","GH5_34","GH5_35","GH5_36","GH5_37","GH5_38","GH5_39",
       "GH5_40","GH5_41","GH5_42","GH5_43","GH5_44","GH5_45","GH5_46","GH5_47","GH5_48","GH5_49",
       "GH5_50","GH5_51","GH5_52","GH5_53","GH5_54","GH5_55","GH5_56")
GH13<-c("GH13_1","GH13_2","GH13_3","GH13_4","GH13_5","GH13_6","GH13_7","GH13_8","GH13_9",
        "GH13_10","GH13_11","GH13_12","GH13_13","GH13_14","GH13_15","GH13_16","GH13_17","GH13_18","GH13_19",
        "GH13_20","GH13_21","GH13_22","GH13_23","GH13_24","GH13_25","GH13_26","GH13_27","GH13_28","GH13_29",
        "GH13_30","GH13_31","GH13_32","GH13_33","GH13_34","GH13_35","GH13_36","GH13_37","GH13_38","GH13_39",
        "GH13_40","GH13_41","GH13_42")
GH16<-c("GH16_1","GH16_2","GH16_3","GH16_4","GH16_5","GH16_6","GH16_7","GH16_8","GH16_9",
        "GH16_10","GH16_11","GH16_12","GH16_13","GH16_14","GH16_15","GH16_16","GH16_17","GH16_18","GH16_19",
        "GH16_20","GH16_21","GH16_22","GH16_23")
GH30<-c("GH30_1","GH30_2","GH30_3","GH30_4","GH30_5","GH30_6","GH30_7","GH30_8","GH30_9")
GH43<-c("GH43_1","GH43_2","GH43_3","GH43_4","GH43_5","GH43_6","GH43_7","GH43_8","GH43_9",
        "GH43_10","GH43_11","GH43_12","GH43_13","GH43_14","GH43_15","GH43_16","GH43_17","GH43_18","GH43_19",
        "GH43_20","GH43_21","GH43_22","GH43_23","GH43_24","GH43_25","GH43_26","GH43_27","GH43_28","GH43_29",
        "GH43_30","GH43_31","GH43_32","GH43_33","GH43_34","GH43_35","GH43_36","GH43_37")
#
affila<-c('SLH','dockerin','cohesin')
LPMO<-c('AA9','AA10','AA11')
exo<-c('GH6','GH9','GH48',GH5,'GH5')
endo<-c('GH7','GH8','GH12','GH44','GH45','GH51','GH61','GH124')
#cellulose-active GH modules, i.e., both exoglucanase and endoglucanase GH modules 
cGH<-c(exo,endo,LPMO)
#cellulose-binding CBM modules
cCBM<-c("CBM1","CBM2","CBM3","CBM4","CBM6","CBM7","CBM8","CBM9","CBM10","CBM11","CBM16","CBM17","CBM28","CBM30","CBM37","CBM44","CBM46","CBM49","CBM63","CBM64")
#hemicellulose-active GH modules
xylanase_GH<-c('GH10','GH11','GH16', GH16,'GH26','GH28','GH53','GH74','GH113','GH134')
#oligosaccharide-active GH modules
oligosaccharidase<-c('GH1','GH2','GH3','GH29','GH30', GH30,'GH35','GH38','GH39','GH42','GH43', GH43,'GH47','GH52','GH59','GH92','GH94','GH116','GH120','GH125')
#debranching_GH_modules<-c('GH54','GH62','GH67','GH78')
#-----------------------------------------------------------------------------------------------

#----------------------------------------------------------------
#step2: summarize the abundance of each CAZy module identified in each genome 
cazy_abun_in_genomes<-data.frame(cazy_module=c(NA))

for (i in 1:length(cazy_in_genomes))
{
  temp<-read.delim(file=cazy_in_genomes[i],sep='\t',header = F,stringsAsFactors = F)
  colnames(temp)<-c("Family_HMM","HMM_length","Query_ID", "Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  cazy_abun<-as.data.frame(table(temp$Family_HMM))
  gbk_name<-strsplit(cazy_in_genomes[i],".faa.out.dm.ps") [[1]][1]
  colnames(cazy_abun)<-c("cazy_module",gbk_name)
  cazy_abun_in_genomes<-merge(cazy_abun_in_genomes,cazy_abun,by="cazy_module",all=T)
}


cazy_abun_in_genomes$cazy_module<-gsub(".hmm","",cazy_abun_in_genomes$cazy_module)
cazy_abun_in_genomes<-cazy_abun_in_genomes[!(is.na(cazy_abun_in_genomes$cazy_module)),]

cazy_abun_in_genomes_t<-as.data.frame(t(cazy_abun_in_genomes))

cnames<-cazy_abun_in_genomes[,1]

colnames(cazy_abun_in_genomes_t)<-cnames
cazy_abun_in_genomes_1<-cazy_abun_in_genomes_t[-1,]
cazy_abun_in_genomes_1[is.na(cazy_abun_in_genomes_1)]<-c('0')


#calculate the diversity and the abundance of the CAZy modules identified in each genomes
#------------------------------------------------------------------------------------------------
all_cazy<-colnames(cazy_abun_in_genomes_1)
all_CBM<-all_cazy[grep('CBM',all_cazy)]
all_GH<-all_cazy[grep('GH',all_cazy)]
#
sle_GH<-c(exo,endo,xylanase_GH,oligosaccharidase)
other_GH<-all_GH[!(all_GH %in% sle_GH)]
other_CBM<-all_CBM[!(all_CBM %in% cCBM)]
#
sle_cazy<-c(affila,LPMO,exo,endo,xylanase_GH,oligosaccharidase,cCBM,other_GH,other_CBM)
other_CAZys<-all_cazy[!(all_cazy %in% sle_cazy)]
all_cazy_inseq<-c(sle_cazy,other_CAZys)
#
seq_cazy_abun_in_genomes_0<-cazy_abun_in_genomes_1[,all_cazy_inseq[all_cazy_inseq %in% colnames(cazy_abun_in_genomes_1)]]
rnam<-rownames(seq_cazy_abun_in_genomes_0)
seq_cazy_abun_in_genomes_numeric<-apply(seq_cazy_abun_in_genomes_0,2,function(x) as.numeric(x))
rownames(seq_cazy_abun_in_genomes_numeric)<-rnam
seq_cazy_abun_in_genomes<-as.data.frame(seq_cazy_abun_in_genomes_numeric)

LPMO_cazy<-seq_cazy_abun_in_genomes[c(colnames(seq_cazy_abun_in_genomes) %in% LPMO)]
seq_cazy_abun_in_genomes$sum_LPMO=rowSums(LPMO_cazy,na.rm = T,dims = 1)
exo_cazy<-seq_cazy_abun_in_genomes[c(colnames(seq_cazy_abun_in_genomes) %in% exo)]
seq_cazy_abun_in_genomes$sum_exo=rowSums(exo_cazy,na.rm = T,dims = 1)
endo_cazy<-seq_cazy_abun_in_genomes[c(colnames(seq_cazy_abun_in_genomes) %in% endo)]
seq_cazy_abun_in_genomes$sum_endo=rowSums(endo_cazy,na.rm = T,dims = 1)
sel_CBM_cazy<-seq_cazy_abun_in_genomes[c(colnames(seq_cazy_abun_in_genomes) %in% cCBM)]
seq_cazy_abun_in_genomes$sum_cCBM=rowSums(sel_CBM_cazy,na.rm = T,dims = 1)
xylanase_GH_cazy<-seq_cazy_abun_in_genomes[c(colnames(seq_cazy_abun_in_genomes) %in% xylanase_GH)]
seq_cazy_abun_in_genomes$sum_xylanase_GH<-rowSums(xylanase_GH_cazy,na.rm = T,dims=1)
oligosacc_cazy<-seq_cazy_abun_in_genomes[c(colnames(seq_cazy_abun_in_genomes) %in% oligosaccharidase)]
seq_cazy_abun_in_genomes$sum_oligosaccharidase_GH<-rowSums(oligosacc_cazy,na.rm = T,dims=1)
other_GH_cazy<-seq_cazy_abun_in_genomes[c(colnames(seq_cazy_abun_in_genomes) %in% other_GH)]
seq_cazy_abun_in_genomes$sum_other_GH<-rowSums(other_GH_cazy,na.rm = T,dims=1)
other_CBM_cazy<-seq_cazy_abun_in_genomes[c(colnames(seq_cazy_abun_in_genomes) %in% other_CBM)]
seq_cazy_abun_in_genomes$sum_other_CBM<-rowSums(other_CBM_cazy,na.rm = T,dims=1)
#SLH_cazy<-seq_cazy_abun_in_genomes[c(colnames(seq_cazy_abun_in_genomes) %in% c('SLH'))]
#seq_cazy_abun_in_genomes$sum_SLH<-rowSums(SLH_cazy,na.rm = T,dims=1)
#dok_cazy<-seq_cazy_abun_in_genomes[c(colnames(seq_cazy_abun_in_genomes) %in% c('dockerin'))]
#seq_cazy_abun_in_genomes$sum_dok<-rowSums(dok_cazy,na.rm = T,dims=1)
#cohe_cazy<-seq_cazy_abun_in_genomes[c(colnames(seq_cazy_abun_in_genomes) %in% c('cohesin'))]
#seq_cazy_abun_in_genomes$sum_cohe<-rowSums(cohe_cazy,na.rm = T,dims=1)


added_columns<-c("sum_exo","sum_endo","sum_cGH","sum_LPMO", "sum_xylanase_GH","sum_oligosaccharidase","sum_cCBM","sum_other_GH","sum_other_CBM")

seq_cazy_abun_in_genomes$genome<-rownames(seq_cazy_abun_in_genomes)
new_colnames<-c("genome",affila,added_columns,LPMO,exo,endo,cCBM,xylanase_GH,oligosaccharidase,other_GH,other_CBM,other_CAZys)
seq_cazy_abun_in_genomes_1<-seq_cazy_abun_in_genomes[new_colnames[new_colnames %in% colnames(seq_cazy_abun_in_genomes)]]
seq_cazy_abun_in_genomes_1[is.na(seq_cazy_abun_in_genomes_1)]<-0
#dim(seq_cazy_abun_in_genomes_1)
#[1] 3898  331
write.table(seq_cazy_abun_in_genomes_1,file="summary_cazy_in_genomes.csv",row.names = F, col.names = T, sep = ',', quote=F)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
#Categorization of genomes in Group II
#Step3':
#Preliminary categorization of the genomes
#according to whether they may habor both exoglucanase and endoglucanase GH modules (GroupI), otherwise, Group II
exo_endo_harboring_genomes<-seq_cazy_abun_in_genomes_1[seq_cazy_abun_in_genomes_1$sum_exo>=1 & seq_cazy_abun_in_genomes_1$sum_endo>=1,]
dim(exo_endo_harboring_genomes)
#[1]  10 181
#check there is genomes harboring both exo and endo glucanase GH modules, 10>0
other_genomes<-seq_cazy_abun_in_genomes_1[seq_cazy_abun_in_genomes_1$sum_exo==0 | seq_cazy_abun_in_genomes_1$sum_endo==0,]
dim(other_genomes)
#[1] 5 181
#------------------------------------------------------------------------------------------------
#Step4': 
#Genotype-based categorization of genomes in GroupI (with both the exo- and the endo-glucanase GH modules)
#if the 'exo_endo_harboring_genomes' is not empty
#------------------------------------------------------------------------------------------------
# read CAZy modules annotated along genes in each genome
exo_endo<-rownames(exo_endo_harboring_genomes)
cazy_in_exo_endo_harboring_genomes<-list()
for (i in 1:length(exo_endo))
{
  #if the annotation result file have the .faa.out.dm.ps suffix
  tmp<-read.delim(file=paste(exo_endo[i], ".faa.out.dm.ps", sep=""), sep='\t',header = F,stringsAsFactors = F)
  #if the annotation result file does not have the .faa.out.dm.ps suffix
  #tmp<-read.delim(file=exo_endo[i], sep='\t',header = F,stringsAsFactors = F)
  colnames(tmp)<-c("Family_HMM","HMM_length","Query_ID", "Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  tmp$Family_HMM<-gsub(".hmm",'',tmp$Family_HMM)
  bin_name<-exo_endo[i]
  cazy_in_exo_endo_harboring_genomes[[bin_name]]<-tmp
}

#### categorize CAZy genes and count
exo_endo_harboring_genomes_names<-names(cazy_in_exo_endo_harboring_genomes)
#
Abundance_of_cazy_gene_exo_endo<-data.frame("genome"=exo_endo_harboring_genomes_names)
#A-1
Abundance_of_cazy_gene_exo_endo$attach_scaffold<-c(0)
#A-2
Abundance_of_cazy_gene_exo_endo$attach_scaffold_0<-c(0)
Abundance_of_cazy_gene_exo_endo$attach_scaffold_0t<-c(0)
#A-3
Abundance_of_cazy_gene_exo_endo$free_scaffold<-c(0)
#B
Abundance_of_cazy_gene_exo_endo$SLH_cCBM<-c(0)
#B'
Abundance_of_cazy_gene_exo_endo$SLH_CBM<-c(0)
#dok_cCBM B2-t
Abundance_of_cazy_gene_exo_endo$dok_cCBM<-c(0)
#dok_CBM B2'-t
Abundance_of_cazy_gene_exo_endo$dok_CBM<-c(0)
#C
Abundance_of_cazy_gene_exo_endo$cGH_cCBM<-c(0)
#C
Abundance_of_cazy_gene_exo_endo$cGH_otherCBM<-c(0)
#C''
Abundance_of_cazy_gene_exo_endo$cGH_CBM<-c(0)
#C'''
Abundance_of_cazy_gene_exo_endo$GH_CBM<-c(0)
#C'
Abundance_of_cazy_gene_exo_endo$exo_cCBM<-c(0)
#
Abundance_of_cazy_gene_exo_endo$LPMO_CBM<-c(0)
#
Abundance_of_cazy_gene_exo_endo$solitude_LPMO<-c(0)
#D
Abundance_of_cazy_gene_exo_endo$solitude_cGH<-c(0)
#E
Abundance_of_cazy_gene_exo_endo$oligosaccharidase<-c(0)

#str(Abundance_of_cazy_gene_exo_endo)
####
Abundance_of_cazy_gene_exo_endo$SLH_cohesin<-c(0)
Abundance_of_cazy_gene_exo_endo$exo_dok<-c(0)
Abundance_of_cazy_gene_exo_endo$endo_dok<-c(0)
Abundance_of_cazy_gene_exo_endo$xylanase_GH_dok<-c(0)
Abundance_of_cazy_gene_exo_endo$oligosaccharidase_dok<-c(0)

Abundance_of_cazy_gene_exo_endo$LPMO_dok<-c(0)
Abundance_of_cazy_gene_exo_endo$otherGH_dok<-c(0)

#exo_endo_harboring_genomes_names[1]
#[1] "NC_003030"
#####
#j=1, i=1
for (j in 1: length(exo_endo_harboring_genomes_names))
{
  tmp0<-cazy_in_exo_endo_harboring_genomes[[j]]
  query_id_in_tmp0<-tmp0[,'Query_ID']
  uni_query_id_in_tmp0<-unique(query_id_in_tmp0)
  for (i in 1:length(uni_query_id_in_tmp0))
  {
    tmp_cazy<-tmp0$Family_HMM[tmp0$Query_ID %in% uni_query_id_in_tmp0[i]]
    tmp_cazy_df<-as.data.frame(table(tmp_cazy))
    tmp_cazy_df$Freq<-as.numeric(tmp_cazy_df$Freq)
    #constructing a dataframe 'a00' which contains info on all the cazy modules i need 
    #for the judgment on the different mechanisms, otherwise, 
    #the modules in if clause may not be recognized the table generated for some specific genes 
    aoo<-data.frame("tmp_cazy"=c(affila,LPMO,exo,endo,cCBM,xylanase_GH,oligosaccharidase, other_CBM),"Freq"=c(0))
    aoo_exclu<-aoo[!(aoo$tmp_cazy %in% tmp_cazy_df$tmp_cazy),]
    tmp_cazy_df_expan<-rbind(tmp_cazy_df,aoo_exclu)
    # A1 attached cellulosome coh>=3, slh>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"attach_scaffold"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"attach_scaffold"]+1
    }
    # A2 attached cellulosome coh>=3, slh>=1, dok>=1
    #if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1)
    #{
    #  Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"attach_scaffold"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"attach_scaffold"]+1
    #}
    ## A2-a attached cellulosome_BACKBOND:     coh>=3, slh==0, dok>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']==0)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"attach_scaffold_0"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"attach_scaffold_0"]+1
    }
    ## A2- attached cellulosome_adhesin_tail:     coh>=1, slh>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"attach_scaffold_0t"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"attach_scaffold_0t"]+1
    }
    
    # A3 free cGH:     coh>=3, dok==0, slh==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']==0)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"free_scaffold"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"free_scaffold"]+1
    }
    # B SLH_cCBM(_gh):      slh>=1, sum(selcbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cCBM])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"SLH_cCBM"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"SLH_cCBM"]+1
    }
    # B' SLH_CBM(_gh):      slh>=1, sum(cbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"SLH_CBM"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"SLH_CBM"]+1
    }
    # B2-a dok_cCBM:      dok>=1, sum(selcbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cCBM])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"dok_cCBM"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"dok_cCBM"]+1
    }
    # B2'-a dok_CBM(_gh):      slh>=1, sum(cbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"dok_CBM"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"dok_CBM"]+1
    }
    # C cGH_cCBM(_gh):     sum(slegh>=1, sum(selcbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cGH])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cCBM])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"cGH_cCBM"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"cGH_cCBM"]+1
    }
    #cGH_otherCBM
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cGH])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% other_CBM])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"cGH_otherCBM"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"cGH_otherCBM"]+1
    }
    # C' cGH_cbm(_gh):     sum(slegh>=1, sum(cbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cGH])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"cGH_CBM"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"cGH_CBM"]+1
    }
    # C' GH_cbm(_gh):     sum(gh>=1, sum(cbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_GH])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"GH_CBM"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"GH_CBM"]+1
    }
    # D solitude_cGH(_slh):  sum(slegh>=1, sum(selcbm)==0
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cGH])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cCBM])==0)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"solitude_cGH"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"solitude_cGH"]+1 
    }
    # C' LPMO_cbm(CBM_LPMO):     sum(LPMO>=1, sum(allcbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% LPMO])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"LPMO_CBM"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"LPMO_CBM"]+1
    }
    # D' solitude_LPMO(LPMO):  sum(LPMO>=1, sum(allcbm)==0
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% LPMO])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])==0)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"solitude_LPMO"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"solitude_LPMO"]+1 
    }
    
    #exo-CBM
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cCBM])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"exo_cCBM"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"exo_cCBM"]+1 
    }
    
    #oligosaccharidase
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"oligosaccharidase"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"oligosaccharidase"]+1 
    }
    ######
    #SLH_cohesin
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"SLH_cohesin"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"SLH_cohesin"]+1 
    }
    #exo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"exo_dok"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"exo_dok"]+1 
    }
    #endo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% endo])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"endo_dok"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"endo_dok"]+1 
    }
    #xylanase_GH_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% xylanase_GH])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"xylanase_GH_dok"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"xylanase_GH_dok"]+1 
    }
    #oligosaccharidase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"oligosaccharidase_dok"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"oligosaccharidase_dok"]+1 
    }
    #LPMO_dok 
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% LPMO])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"LPMO_dok"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"LPMO_dok"]+1 
    }
    #otherGH_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% other_GH])>=1)
    {
      Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"otherGH_dok"]=Abundance_of_cazy_gene_exo_endo[Abundance_of_cazy_gene_exo_endo$genome %in% exo_endo_harboring_genomes_names[j],"otherGH_dok"]+1 
    }
  }
}

#------------------------------------------------------------------------------------------------

######

GroupI_genome_categorization<-Abundance_of_cazy_gene_exo_endo

##### assign each draft genome to a category based on the abundance of gene patterns identified in it
GroupI_genome_categorization<-within(GroupI_genome_categorization,{
  category<-NA
  category[attach_scaffold>=1]<-"Group I-a"
  #category[attach_scaffold_0>=1] <-"Group I-a"
  category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t>=1] <-"Group I-a"
  category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t==0 & SLH_CBM>=1] <-"Group I-b"
  #category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t==0 & SLH_cCBM>=1] <-"Group I-b"
  category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t==0 & SLH_CBM==0] <-"Group I-b2"
  category[attach_scaffold==0 & attach_scaffold_0 ==0 & free_scaffold>=1 & SLH_CBM>=1]<-"Group I-b"
  #category[attach_scaffold==0 & attach_scaffold_0 ==0 & free_scaffold>=1 & SLH_cCBM>=1]<-"Group I-b"
  category[attach_scaffold==0 & attach_scaffold_0 ==0 & free_scaffold>=1 & attach_scaffold_0t>=1 & dok_CBM>=1]<-"Group I-b"
  category[attach_scaffold==0 & attach_scaffold_0 ==0 & free_scaffold>=1 & SLH_CBM==0 & attach_scaffold_0t>=1 & dok_CBM==0]<-"Group I-b2"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold>=1 & SLH_CBM==0 & attach_scaffold_0t==0]<-"Group I-c"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_CBM>=1]<-"Group I-d"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & attach_scaffold_0t>=1 & dok_CBM>=1]<-"Group I-d"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_CBM==0 & attach_scaffold_0t>=1 & dok_CBM==0]<-"Group I-d2"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_cCBM==0 & attach_scaffold_0t==0 & cGH_cCBM>=1]<-"Group I-e"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_cCBM==0 & attach_scaffold_0t==0 & cGH_cCBM==0]<-"Group I-f"
  })
#
GroupI_CAZyinfo_complete<-merge(GroupI_genome_categorization,seq_cazy_abun_in_genomes_1,by="genome",all=F)
GroupI_CAZyinfo_complete_1<-GroupI_CAZyinfo_complete %>%
  select(category, everything())
write.table(GroupI_CAZyinfo_complete_1, file="GroupI_genome_categorization.csv", sep=",", row.names=F, col.names=T, quote=F)

#------------------group II-------------------------------------------------------
#Categorization of genomes in Group II
#Step3':

groupII_genomes<-seq_cazy_abun_in_genomes_1[seq_cazy_abun_in_genomes_1$sum_exo==0 | seq_cazy_abun_in_genomes_1$sum_endo==0,]
dim(groupII_genomes)
#[1] 5 181
#------------------------------------------------------------------------------------------------
#Step4': 
#Genotype-based categorization of genomes in groupII (with both the exo- and the endo-glucanase GH modules)
#------------------------------------------------------------------------------------------------
# read CAZy modules annotated along genes in each genome
groupII<-rownames(groupII_genomes)

cazy_in_groupII_genomes<-list()
for (i in 1:length(groupII))
{
  tmp<-read.delim(file=paste(groupII[i], ".faa.out.dm.ps", sep=""), sep='\t',header = F,stringsAsFactors = F)
  #tmp<-read.delim(file=groupII[i], sep='\t',header = F,stringsAsFactors = F)
  colnames(tmp)<-c("Family_HMM","HMM_length","Query_ID", "Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  tmp$Family_HMM<-gsub(".hmm",'',tmp$Family_HMM)
  bin_name<-groupII[i]
  cazy_in_groupII_genomes[[bin_name]]<-tmp
}

#### categorize CAZy genes and count
groupII_genomes_names<-names(cazy_in_groupII_genomes)
#
Abundance_of_cazy_gene_groupII<-data.frame("genome"=groupII_genomes_names)
#A-1
Abundance_of_cazy_gene_groupII$attach_scaffold<-c(0)
#A-2
Abundance_of_cazy_gene_groupII$attach_scaffold_0<-c(0)
Abundance_of_cazy_gene_groupII$attach_scaffold_0t<-c(0)
#A-3
Abundance_of_cazy_gene_groupII$free_scaffold<-c(0)
#B
Abundance_of_cazy_gene_groupII$SLH_cCBM<-c(0)
#B'
Abundance_of_cazy_gene_groupII$SLH_CBM<-c(0)
#dok_cCBM B2-t
Abundance_of_cazy_gene_groupII$dok_cCBM<-c(0)
#dok_CBM B2'-t
Abundance_of_cazy_gene_groupII$dok_CBM<-c(0)
#C
Abundance_of_cazy_gene_groupII$cGH_cCBM<-c(0)
#C
Abundance_of_cazy_gene_groupII$cGH_otherCBM<-c(0)
#C''
Abundance_of_cazy_gene_groupII$cGH_CBM<-c(0)
#C'''
Abundance_of_cazy_gene_groupII$GH_CBM<-c(0)
#C'
Abundance_of_cazy_gene_groupII$exo_cCBM<-c(0)
#
Abundance_of_cazy_gene_groupII$LPMO_CBM<-c(0)
#
Abundance_of_cazy_gene_groupII$solitude_LPMO<-c(0)
#D
Abundance_of_cazy_gene_groupII$solitude_cGH<-c(0)
#E
Abundance_of_cazy_gene_groupII$oligosaccharidase<-c(0)

#str(Abundance_of_cazy_gene_groupII)
####
#Abundance_of_cazy_gene_groupII$SLH_cohesin<-c(0)
Abundance_of_cazy_gene_groupII$exo_dok<-c(0)
Abundance_of_cazy_gene_groupII$endo_dok<-c(0)
Abundance_of_cazy_gene_groupII$xylanase_GH_dok<-c(0)
Abundance_of_cazy_gene_groupII$oligosaccharidase_dok<-c(0)

Abundance_of_cazy_gene_groupII$LPMO_dok<-c(0)
Abundance_of_cazy_gene_groupII$otherGH_dok<-c(0)
Abundance_of_cazy_gene_groupII$sle_GH<-c(0)

#groupII_genomes_names[121]
#[1] "NC_016609"
#####
#j=5, i=121
for (j in 1: length(groupII_genomes_names))
{
  tmp0<-cazy_in_groupII_genomes[[j]]
  query_id_in_tmp0<-tmp0[,'Query_ID']
  uni_query_id_in_tmp0<-unique(query_id_in_tmp0)
  for (i in 1:length(uni_query_id_in_tmp0))
  {
    tmp_cazy<-tmp0$Family_HMM[tmp0$Query_ID %in% uni_query_id_in_tmp0[i]]
    tmp_cazy_df<-as.data.frame(table(tmp_cazy))
    tmp_cazy_df$Freq<-as.numeric(tmp_cazy_df$Freq)
    #constructing a dataframe 'a00' which contains info on all the cazy modules i need 
    #for the judgment on the different mechanisms, otherwise, 
    #the modules in if clause may not be recognized the table generated for some specific genes 
    aoo<-data.frame("tmp_cazy"=c(affila,LPMO,exo,endo,cCBM,xylanase_GH,oligosaccharidase, other_CBM),"Freq"=c(0))
    aoo_exclu<-aoo[!(aoo$tmp_cazy %in% tmp_cazy_df$tmp_cazy),]
    tmp_cazy_df_expan<-rbind(tmp_cazy_df,aoo_exclu)
    # A1 attached cellulosome coh>=3, slh>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"attach_scaffold"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"attach_scaffold"]+1
    }
    # A2 attached cellulosome coh>=3, slh>=1, dok>=1
    #if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1)
    #{
    #  Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"attach_scaffold"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"attach_scaffold"]+1
    #}
    ## A2-a attached cellulosome_BACKBOND:     coh>=3, slh==0, dok>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']==0)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"attach_scaffold_0"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"attach_scaffold_0"]+1
    }
    ## A2-b attached cellulosome_adhesin_tail:     coh>=1, slh>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"attach_scaffold_0t"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"attach_scaffold_0t"]+1
    }
    
    # A3 free cGH:     coh>=3, dok==0, slh==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']==0)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"free_scaffold"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"free_scaffold"]+1
    }
    # B SLH_cCBM(_gh):      slh>=1, sum(selcbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cCBM])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"SLH_cCBM"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"SLH_cCBM"]+1
    }
    # B' SLH_CBM(_gh):      slh>=1, sum(cbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"SLH_CBM"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"SLH_CBM"]+1
    }
    # B2-a dok_cCBM:      dok>=1, sum(selcbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cCBM])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"dok_cCBM"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"dok_cCBM"]+1
    }
    # B2'-a dok_CBM(_gh):      slh>=1, sum(cbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"dok_CBM"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"dok_CBM"]+1
    }
    # C cGH_cCBM(_gh):     sum(slegh>=1, sum(selcbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cGH])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cCBM])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"cGH_cCBM"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"cGH_cCBM"]+1
    }
    #cGH_otherCBM
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cGH])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% other_CBM])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"cGH_otherCBM"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"cGH_otherCBM"]+1
    }
    # C' cGH_cbm(_gh):     sum(slegh>=1, sum(cbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cGH])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"cGH_CBM"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"cGH_CBM"]+1
    }
    # C' GH_cbm(_gh):     sum(gh>=1, sum(cbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_GH])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"GH_CBM"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"GH_CBM"]+1
    }
    # D solitude_cGH(_slh):  sum(slegh>=1, sum(selcbm)==0
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cGH])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cCBM])==0)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"solitude_cGH"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"solitude_cGH"]+1 
    }
    # C' LPMO_cbm(CBM_LPMO):     sum(LPMO>=1, sum(allcbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% LPMO])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"LPMO_CBM"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"LPMO_CBM"]+1
    }
    # D' solitude_LPMO(LPMO):  sum(LPMO>=1, sum(allcbm)==0
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% LPMO])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% all_CBM])==0)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"solitude_LPMO"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"solitude_LPMO"]+1 
    }
    
    #exo-CBM
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cCBM])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"exo_cCBM"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"exo_cCBM"]+1 
    }
    
    #oligosaccharidase
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"oligosaccharidase"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"oligosaccharidase"]+1 
    }
    ######
    #SLH_cohesin
    #if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=1)
    #{
    #  Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"SLH_cohesin"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"SLH_cohesin"]+1 
    #}
    #exo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"exo_dok"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"exo_dok"]+1 
    }
    #endo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% endo])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"endo_dok"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"endo_dok"]+1 
    }
    #xylanase_GH_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% xylanase_GH])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"xylanase_GH_dok"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"xylanase_GH_dok"]+1 
    }
    #oligosaccharidase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"oligosaccharidase_dok"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"oligosaccharidase_dok"]+1 
    }
    #LPMO_dok 
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% LPMO])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"LPMO_dok"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"LPMO_dok"]+1 
    }
    #otherGH_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% other_GH])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"otherGH_dok"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"otherGH_dok"]+1 
    }
    #Abundance_of_cazy_gene_groupII$sle_GH<-c(0)
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% sle_GH])>=1)
    {
      Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"sle_GH"]=Abundance_of_cazy_gene_groupII[Abundance_of_cazy_gene_groupII$genome %in% groupII_genomes_names[j],"sle_GH"]+1 
    }
  }
}

#-------------------------------------------------------------------------------------------------
groupII_genome_categorization<-Abundance_of_cazy_gene_groupII

##### assign each draft genome to a category based on the abundance of gene patterns identified in it
#groupII_genome_categorization<-within(groupII_genome_categorization,{
#  category<-NA
#  category[attach_scaffold>=1]<-"Group II-a"
#  category[attach_scaffold_0>=1 & attach_scaffold_0t>=1] <-"Group II-a"
#  category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t==0 & SLH_cCBM>=1] <-"Group II-b"
#  category[attach_scaffold==0 & attach_scaffold_0 ==0 & free_scaffold>=1 & SLH_cCBM>=1]<-"Group II-b"
#  category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t==0 & SLH_cCBM==0]<-"Group II-c"
#  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold>=1 & SLH_cCBM==0 ]<-"Group II-c"
#  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_cCBM>=1]<-"Group II-d"
#  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_cCBM==0 & sle_GH>=1]<-"Group II-e"
#  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_cCBM==0 & sle_GH==0]<-"Group II-f"
#})

groupII_genome_categorization<-within(groupII_genome_categorization,{
  category<-NA
  category[attach_scaffold>=1]<-"Group II-a"
  #category[attach_scaffold_0>=1] <-"Group II-a"
  category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t>=1] <-"Group II-a"
  category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t==0 & SLH_CBM>=1] <-"Group II-b"
  #category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t==0 & SLH_cCBM>=1] <-"Group II-b"
  category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t==0 & SLH_CBM==0] <-"Group II-b2"
  category[attach_scaffold==0 & attach_scaffold_0 ==0 & free_scaffold>=1 & SLH_CBM>=1]<-"Group II-b"
  #category[attach_scaffold==0 & attach_scaffold_0 ==0 & free_scaffold>=1 & SLH_cCBM>=1]<-"Group II-b"
  category[attach_scaffold==0 & attach_scaffold_0 ==0 & free_scaffold>=1 & attach_scaffold_0t>=1 & dok_CBM>=1]<-"Group II-b"
  category[attach_scaffold==0 & attach_scaffold_0 ==0 & free_scaffold>=1 & SLH_CBM==0 & attach_scaffold_0t>=1 & dok_CBM==0]<-"Group II-b2"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold>=1 & SLH_CBM==0 & attach_scaffold_0t==0]<-"Group II-c"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_CBM>=1]<-"Group II-d"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & attach_scaffold_0t>=1 & dok_CBM>=1]<-"Group II-d"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_CBM==0 & attach_scaffold_0t>=1 & dok_CBM==0]<-"Group II-d2"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_cCBM==0 & attach_scaffold_0t==0 & cGH_cCBM>=1]<-"Group II-e"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_cCBM==0 & attach_scaffold_0t==0 & cGH_cCBM==0]<-"Group II-f"
})


#
groupII_CAZyinfo_complete<-merge(groupII_genome_categorization,seq_cazy_abun_in_genomes_1,by="genome",all=F)

groupII_CAZyinfo_complete_1<-groupII_CAZyinfo_complete %>%
  select(category, everything())

write.table(groupII_CAZyinfo_complete_1,file="GroupII_genome_categorization.csv",sep=',',row.names = F,col.names = T)

#-----------------END------------------------------------------------------------------------------------------------