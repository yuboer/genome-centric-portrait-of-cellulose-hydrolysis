#---------------------------------------------------------------
#Step1: read the dbCAN annotation results into the R environment
#---------------------------------------------------------------
#set the working directory as the directory that contains the dbCAN annotation results, i.e., the folder with the '*.out.dm.ps' files

address_dbCAN_results = as.character("/projects/b1052/Wells_b1042/yuboer/ongoing/sea_strains_checkM_results/bins/sea_faa/dbCAN_annotation_results")
setwd(address_dbCAN_results)

#check
getwd()

cazy_in_bins<-dir(address_dbCAN_results)

#----------------------------------------------------------------
#step1 summary on the statistics of CAZy modules identified in each draft genome 
# abundace of CAZy modules identified in each draft genome
cazy_abun_in_bins<-data.frame(cazy_module=c(NA))

for (i in 1:length(cazy_in_bins))
{
  temp<-read.delim(file=cazy_in_bins[i],sep='\t',header = F,stringsAsFactors = F)
  colnames(temp)<-c("Family_HMM","HMM_length","Query_ID", "Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  cazy_abun<-as.data.frame(table(temp$Family_HMM))
  gbk_name<-strsplit(cazy_in_bins[i],".faa.out.dm.ps") [[1]][1]
  colnames(cazy_abun)<-c("cazy_module",gbk_name)
  cazy_abun_in_bins<-merge(cazy_abun_in_bins,cazy_abun,by="cazy_module",all=T)
}


cazy_abun_in_bins$cazy_module<-gsub(".hmm","",cazy_abun_in_bins$cazy_module)
cazy_abun_in_bins<-cazy_abun_in_bins[!(is.na(cazy_abun_in_bins$cazy_module)),]

cazy_abun_in_bins_t<-as.data.frame(t(cazy_abun_in_bins))

cnames<-cazy_abun_in_bins[,1]

colnames(cazy_abun_in_bins_t)<-cnames
cazy_abun_in_bins_1<-cazy_abun_in_bins_t[-1,]

#----------------------------------------------------------------------------------------------------------------------------
#*Categorization of the CAZy modules, accoridng to updated info on the functions of CAZy modules (this part could be manually updated/modified if the description of each CAZy module in the CAZy database updated)
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
exo<-c('GH6','GH9','GH48')
endo<-c(GH5,'GH7','GH8','GH12','GH44','GH45','GH51','GH61','GH124')
cellulase<-c(exo,endo)
cellulose_binding_CBM<-c("CBM1","CBM2","CBM3","CBM4","CBM6","CBM7","CBM8","CBM9","CBM10","CBM11","CBM16","CBM17","CBM28","CBM30","CBM37","CBM44","CBM46","CBM49","CBM63","CBM64")
hemicellulase<-c('GH10','GH11',GH16,'GH26','GH28','GH53','GH74','GH113','GH134')

#GH43_11	GH43_22	GH43_23	GH43_24	GH43_3	GH43_30	GH43_32	
#GH5_1	GH5_28	GH5_29

oligosaccharidase<-c('GH1','GH2','GH3','GH29',GH30,'GH35','GH38','GH39','GH42',GH43,'GH47','GH52','GH59','GH92','GH94','GH116','GH120','GH125')
#GH_debranching<-c('GH54','GH62','GH67','GH78')


#-----------------------------------------------------------------------------------------------
#Step2: Summarize and make the summary information much more organized in describing the diversity and the abundance of the CAZy modules identified in each genomes
#------------------------------------------------------------------------------------------------
all_cazy<-colnames(cazy_abun_in_bins_1)
all_CBM<-all_cazy[grep('CBM',all_cazy)]
#cellulose_binding_CBM<-c("CBM1","CBM2","CBM3","CBM4","CBM6","CBM7","CBM8","CBM9","CBM10","CBM11","CBM16","CBM17","CBM28","CBM30","CBM37","CBM44","CBM46","CBM49","CBM63","CBM64")
other_CBM<-all_CBM[!(all_CBM %in% cellulose_binding_CBM)]
sel_cazy<-c(affila,exo,endo,hemicellulase,oligosaccharidase,cellulose_binding_CBM,other_CBM)
others<-all_cazy[!(all_cazy %in% sel_cazy)]
all_cazy_inseq<-c(sel_cazy,others)
seq_cazy_abun_in_bins_0<-cazy_abun_in_bins_1[,all_cazy_inseq[all_cazy_inseq %in% colnames(cazy_abun_in_bins_1)]]
rnam<-rownames(seq_cazy_abun_in_bins_0)
seq_cazy_abun_in_bins_numeric<-apply(seq_cazy_abun_in_bins_0,2,function(x) as.numeric(x))
rownames(seq_cazy_abun_in_bins_numeric)<-rnam
seq_cazy_abun_in_bins<-as.data.frame(seq_cazy_abun_in_bins_numeric)

exo_cazy<-seq_cazy_abun_in_bins[c(colnames(seq_cazy_abun_in_bins) %in% exo)]
seq_cazy_abun_in_bins$sum_exo=rowSums(exo_cazy,na.rm = T,dims = 1)
endo_cazy<-seq_cazy_abun_in_bins[c(colnames(seq_cazy_abun_in_bins) %in% endo)]
seq_cazy_abun_in_bins$sum_endo=rowSums(endo_cazy,na.rm = T,dims = 1)
sel_CBM_cazy<-seq_cazy_abun_in_bins[c(colnames(seq_cazy_abun_in_bins) %in% cellulose_binding_CBM)]
seq_cazy_abun_in_bins$sum_cellulose_binding_CBM=rowSums(sel_CBM_cazy,na.rm = T,dims = 1)
hemicellulase_cazy<-seq_cazy_abun_in_bins[c(colnames(seq_cazy_abun_in_bins) %in% hemicellulase)]
seq_cazy_abun_in_bins$sum_hemicellulase<-rowSums(hemicellulase_cazy,na.rm = T,dims=1)
oligosacc_cazy<-seq_cazy_abun_in_bins[c(colnames(seq_cazy_abun_in_bins) %in% oligosaccharidase)]
seq_cazy_abun_in_bins$sum_oligosaccharidase<-rowSums(oligosacc_cazy,na.rm = T,dims=1)
added_columns<-c("sum_exo","sum_endo","sum_cellulase","sum_cellulose_binding_CBM","sum_hemicellulase","sum_oligosaccharidase")

seq_cazy_abun_in_bins$Draft_genome_bin<-rownames(seq_cazy_abun_in_bins)
new_colnames<-c("Draft_genome_bin",affila,added_columns,exo,endo,hemicellulase,oligosaccharidase,cellulose_binding_CBM,other_CBM,others)
seq_cazy_abun_in_bins_1<-seq_cazy_abun_in_bins[new_colnames[new_colnames %in% colnames(seq_cazy_abun_in_bins)]]
#write.table(seq_cazy_abun_in_bins_1,file="seq_cazy_abun_in_bins.csv",row.names = F, col.names = T, sep = ',', quote=F)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
#Step3: Priliminary categorization of the draft genomes, according to whether exoglucanase and endoglucanase GH modules are identified in each of them
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# grouping the draft genomes based on whether exoglucanase, endoglucanase, hemicellulase, oligosaccharidase GH modules are identified
exo_endo_harboring_bins<-seq_cazy_abun_in_bins_1[!seq_cazy_abun_in_bins_1$sum_exo==0 & !seq_cazy_abun_in_bins_1$sum_endo==0,]
only_exo_harboring_bins<-seq_cazy_abun_in_bins_1[!seq_cazy_abun_in_bins_1$sum_exo==0 & seq_cazy_abun_in_bins_1$sum_endo==0,]
only_endo_harboring_bins<-seq_cazy_abun_in_bins_1[seq_cazy_abun_in_bins_1$sum_exo==0 & !seq_cazy_abun_in_bins_1$sum_endo==0,]
only_hemicellulase_harboring_bins<-seq_cazy_abun_in_bins_1[seq_cazy_abun_in_bins_1$sum_exo==0 & seq_cazy_abun_in_bins_1$sum_endo==0 & !seq_cazy_abun_in_bins_1$sum_hemicellulase==0,]
only_oligosaccharidase_harboring_bins<-seq_cazy_abun_in_bins_1[seq_cazy_abun_in_bins_1$sum_exo==0 & seq_cazy_abun_in_bins_1$sum_endo==0 & !seq_cazy_abun_in_bins_1$sum_oligosaccharidase==0, ]
not_carbohydrates_degrading_bins<-seq_cazy_abun_in_bins_1[seq_cazy_abun_in_bins_1$sum_exo==0 & seq_cazy_abun_in_bins_1$sum_endo==0 & seq_cazy_abun_in_bins_1$sum_hemicellulase==0 & seq_cazy_abun_in_bins_1$sum_oligosaccharidase==0,]

dim(exo_endo_harboring_bins)
#[1] 672 549
dim(only_exo_harboring_bins)
#[1] 347 549
dim(only_endo_harboring_bins)
#[1] 1808  549
dim(only_hemicellulase_harboring_bins)
#[1] 2361  549
dim(only_oligosaccharidase_harboring_bins)
#[1] 3931  549
dim(not_carbohydrates_degrading_bins)
#[1] 852 549


#save.image(file = "cazy_abun_in_bins.RData")

#write.table(exo_endo_harboring_bins,file = "seq_cazy_in_exo_endo_harboring_bins.tab",sep="\t",row.names = F,col.names = T)
#write.table(only_exo_harboring_bins,file = "seq_cazy_in_only_exo_harboring_bins.tab",sep="\t",row.names = F,col.names = T)
#write.table(only_endo_harboring_bins,file = "seq_cazy_in_only_endo_harboring_bins.tab",sep="\t",row.names = F,col.names = T)
#write.table(only_hemicellulase_harboring_bins,file = "only_hemicellulase_harboring_bins.tab",sep="\t",row.names = F,col.names = T)
#write.table(only_oligosaccharidase_harboring_bins,file = "seq_cazy_in_only_oligosaccharidase_harboring_bins.tab",sep="\t",row.names = F,col.names = T)
#write.table(not_carbohydrates_degrading_bins,file = "seq_cazy_not_carbohydrates_degrading_bins.tab",sep="\t",row.names = F,col.names = T)

#str(exo_endo_harboring_bins)

# ! ! ! ! Note: please check whether the "exo_endo_harboring_bins" is empty
# ! ! ! ! if the "exo_endo_harboring_bins" as summarized above is not an empty file, in other words, if genomes harboring both exo and endo GH modules were annotated, continue;
# ! ! ! ! otherwise, stop here, as there is no potential cellulose hydrolysers in your dataset; 
# ! ! ! ! one genome must have both exoglucanase and endoglucanase before it should be considered as a genome of potential cellulose hydrolyser, and to be further categorized in the following Step 4

#------------------------------------------------------------------------------------------------
#Step4: Further cetagorize genomes identified with both the exo- and the endo-glucanase GH modules
#------------------------------------------------------------------------------------------------
# read the info on the CAZy modules annotated along genes in each draft genome 
# Genome categorization
exo_endo<-rownames(exo_endo_harboring_bins)

cazy_in_exo_endo_harboring_bins<-list()
for (i in 1:length(exo_endo))
{
  tmp<-read.delim(file=paste(exo_endo[i], ".faa.out.dm.ps", sep=""), sep='\t',header = F,stringsAsFactors = F)
  colnames(tmp)<-c("Family_HMM","HMM_length","Query_ID", "Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  tmp$Family_HMM<-gsub(".hmm",'',tmp$Family_HMM)
  bin_name<-exo_endo[i]
  cazy_in_exo_endo_harboring_bins[[bin_name]]<-tmp
}
#### identify gene patterns in each draft genomes and calculate the abundance of each pattern identified

exo_endo_harboring_bins_names<-names(cazy_in_exo_endo_harboring_bins)

Abundance_of_gene_pattern_exo_endo<-data.frame("bin"=exo_endo_harboring_bins_names)
#A
Abundance_of_gene_pattern_exo_endo$attach_scaffold<-c(0)
#A
Abundance_of_gene_pattern_exo_endo$attach_scaffold_0<-c(0)
Abundance_of_gene_pattern_exo_endo$attach_scaffold_0t<-c(0)
#B
Abundance_of_gene_pattern_exo_endo$free_scaffold<-c(0)
#C
Abundance_of_gene_pattern_exo_endo$SLH_CBM<-c(0)
#D
Abundance_of_gene_pattern_exo_endo$cellulase_CBM<-c(0)
#
Abundance_of_gene_pattern_exo_endo$exo_CBM<-c(0)
#E
Abundance_of_gene_pattern_exo_endo$solitude_cellulase<-c(0)
#F
Abundance_of_gene_pattern_exo_endo$oligosaccharidase<-c(0)

#str(Abundance_of_gene_pattern_exo_endo)
####
Abundance_of_gene_pattern_exo_endo$SLH_cohesin<-c(0)
Abundance_of_gene_pattern_exo_endo$exo_dok<-c(0)
Abundance_of_gene_pattern_exo_endo$endo_dok<-c(0)
Abundance_of_gene_pattern_exo_endo$hemicellulase_dok<-c(0)

Abundance_of_gene_pattern_exo_endo$oligosaccharidase_dok<-c(0)


#####

for (j in 1: length(exo_endo_harboring_bins_names))
{
  tmp0<-cazy_in_exo_endo_harboring_bins[[j]]
  query_id_in_tmp0<-tmp0[,'Query_ID']
  uni_query_id_in_tmp0<-unique(query_id_in_tmp0)
  for (i in 1:length(uni_query_id_in_tmp0))
  {
    tmp_cazy<-tmp0$Family_HMM[tmp0$Query_ID %in% uni_query_id_in_tmp0[i]]
    tmp_cazy_df<-as.data.frame(table(tmp_cazy))
    tmp_cazy_df$Freq<-as.numeric(tmp_cazy_df$Freq)
    #constructing a dataframe 'a00' which contains info on all the cazy modules i need 
    #for the judgment on the different mechanisms, otherwise, 
    #the modules in if clause may not be recoganized the table generated for some specific genes 
    aoo<-data.frame("tmp_cazy"=c(affila,exo,endo,cellulose_binding_CBM,hemicellulase,oligosaccharidase),"Freq"=c(0))
    aoo_exclu<-aoo[!(aoo$tmp_cazy %in% tmp_cazy_df$tmp_cazy),]
    tmp_cazy_df_expan<-rbind(tmp_cazy_df,aoo_exclu)
    # A attached cellulosome coh>=3, slh>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold"]+1
    }
    # A attached cellulosome coh>=3, slh>=1, dok>=1
    #if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1)
    #{
    #  Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold"]+1
    #}
    ## A attached cellulosome_BACKBOND:     coh>=3, slh==0, dok>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']==0)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold_0"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold_0"]+1
    }
    ## A attached cellulosome_adhesin_tail:     coh>=1, slh>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold_0t"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold_0t"]+1
    }
    
    # B free cellulase:     coh>=3, dok==0, slh==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']==0)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"free_scaffold"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"free_scaffold"]+1
    }
    # C slh_cbm(_gh):      slh>=1, sum(selcbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"SLH_CBM"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"SLH_CBM"]+1
    }
    # D cellulase_selcbm(_gh):     sum(slegh>=1, sum(selcbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"cellulase_CBM"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"cellulase_CBM"]+1
    }
    # E solitude_cellulase(_slh):  sum(slegh>=1, sum(selcbm)==0
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])==0)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"solitude_cellulase"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"solitude_cellulase"]+1 
    }
    #exo-CBM
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"exo_CBM"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"exo_CBM"]+1 
    }
    
    #oligosaccharidase
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"oligosaccharidase"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"oligosaccharidase"]+1 
    }
    ######
    #SLH_cohesin
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"SLH_cohesin"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"SLH_cohesin"]+1 
    }
    #exo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"exo_dok"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"exo_dok"]+1 
    }
    #endo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% endo])>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"endo_dok"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"endo_dok"]+1 
    }
    #hemicellulase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% hemicellulase])>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"hemicellulase_dok"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"hemicellulase_dok"]+1 
    }
    #oligosaccharidase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"oligosaccharidase_dok"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"oligosaccharidase_dok"]+1 
    }
  }
}

##### assign each draft genome to a category based on the abundance of gene patterns identified in it
Abundance_of_gene_pattern_exo_endo<-within(Abundance_of_gene_pattern_exo_endo,{
  category<-NA
  category[attach_scaffold>=1]<-"Group I A"
  category[attach_scaffold_0>=1 & attach_scaffold_0t>=1] <-"Group I A"
  category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t==0 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase>=1] <-"Group I B"
  category[attach_scaffold==0 & attach_scaffold_0 ==0 & free_scaffold>=1 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group I B"
  category[attach_scaffold==0 & attach_scaffold_0>=1 & attach_scaffold_0t==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group I C"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold>=1 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group I C"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase==0]<-"Group I D"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_CBM>=1 & solitude_cellulase>=1]<-"Group I D"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group I E"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase==0]<-"Group I E"
  category[attach_scaffold==0 & attach_scaffold_0==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM==0 & solitude_cellulase>=1]<-"Group I F"
})
colnames(Abundance_of_gene_pattern_exo_endo)<-gsub("bin","Draft_genome_bin",colnames(Abundance_of_gene_pattern_exo_endo))

write.table(Abundance_of_gene_pattern_exo_endo, file="Abundance_of_gene_pattern_exo_endo.csv", sep=",", row.names=F, col.names=T, quote=F)

seq_cazy_abun_in_bins_1[is.na(seq_cazy_abun_in_bins_1)]<-0
write.table(seq_cazy_abun_in_bins_1,file="seq_cazy_abun_in_bins.csv",row.names = F, col.names = T, sep = ',', quote=F)
#write.table(Abundance_of_gene_pattern_exo_endo,file="Categorization_and_the_abundance_of_gene_patterns_in_exo_endo_harboring_bins.tab",sep='\t',row.names = F,col.names = T)
#-------------------------------------------------------------------------
