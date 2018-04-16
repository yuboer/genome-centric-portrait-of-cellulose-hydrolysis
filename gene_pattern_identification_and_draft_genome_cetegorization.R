#tips
#1) pay attentin to the setting of the working directory
#2) the wrking directory should be the one where your genome bins locates
#3) some Groups may have no genomes assigned
#---------------------------------------------------------------------
getwd()
#[1] "C:/yuboer folder/ad15/xy_cellulose_genome_cazy/xy_dbcan"
#setwd("/envbio/yuboer/yuboer_2nd_binning/cazy_annotation_2nd_binning")
setwd("C:/yuboer folder/ad15/xy_cellulose_genome_cazy/xy_dbcan")
#xy_dbcan is the folder holding the dbCAN annotation resulting files
getwd()
#[1] "C:/yuboer folder/ad15/xy_cellulose_genome_cazy/xy_dbcan"
cazy_in_bins<-dir("C:/yuboer folder/ad15/xy_cellulose_genome_cazy/xy_dbcan")

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

#-----------------------------------------------------------------------------
#*this part could be updated, 
# accoridng to updated info on the functions of CAZy modules
#--------------------------------------------------------------------
# some primary summary on the data and present the CAZy modules in sequence
# "exo" refers exoglucanase GH modules
# "endo" refers endoglucanase GH modules
all_cazy<-colnames(cazy_abun_in_bins_1)
affila<-c('SLH','dockerin','cohesin')
exo<-c('GH6','GH9','GH48')
endo<-c('GH5','GH7','GH8','GH12','GH44','GH45','GH51','GH61','GH124')
cellulase<-c(exo,endo)
cellulose_binding_CBM<-c("CBM1","CBM2","CBM3","CBM4","CBM6","CBM7","CBM8","CBM9","CBM10","CBM11","CBM16","CBM17","CBM28","CBM30","CBM37","CBM44","CBM46","CBM49","CBM63","CBM64")
hemicellulase<-c('GH10','GH11','GH16','GH26','GH28','GH53','GH74','GH113','GH134')
oligosaccharidase<-c('GH1','GH2','GH3','GH29','GH30','GH35','GH38','GH39','GH42','GH43','GH47','GH52','GH59','GH92','GH94','GH116','GH120','GH125')
#GH_debranching<-c('GH54','GH62','GH67','GH78')
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

write.table(seq_cazy_abun_in_bins_1,file="seq_cazy_abun_in_bins.tab",row.names = F, col.names = T, sep = '\t')
######returns a table with the abundance of CAZy modules identified in drafte genomes

#step2
# grouping the draft genomes based on whether exoglucanase, endoglucanase, hemicellulase, oligosaccharidase GH modules are identified
exo_endo_harboring_bins<-seq_cazy_abun_in_bins_1[!seq_cazy_abun_in_bins_1$sum_exo==0 & !seq_cazy_abun_in_bins_1$sum_endo==0,]
only_exo_harboring_bins<-seq_cazy_abun_in_bins_1[!seq_cazy_abun_in_bins_1$sum_exo==0 & seq_cazy_abun_in_bins_1$sum_endo==0,]
only_endo_harboring_bins<-seq_cazy_abun_in_bins_1[seq_cazy_abun_in_bins_1$sum_exo==0 & !seq_cazy_abun_in_bins_1$sum_endo==0,]
only_hemicellulase_harboring_bins<-seq_cazy_abun_in_bins_1[seq_cazy_abun_in_bins_1$sum_exo==0 & seq_cazy_abun_in_bins_1$sum_endo==0 & !seq_cazy_abun_in_bins_1$sum_hemicellulase==0,]
only_oligosaccharidase_harboring_bins<-seq_cazy_abun_in_bins_1[seq_cazy_abun_in_bins_1$sum_exo==0 & seq_cazy_abun_in_bins_1$sum_endo==0 & !seq_cazy_abun_in_bins_1$sum_oligosaccharidase==0, ]
not_carbohydrates_degrading_bins<-seq_cazy_abun_in_bins_1[seq_cazy_abun_in_bins_1$sum_exo==0 & seq_cazy_abun_in_bins_1$sum_endo==0 & seq_cazy_abun_in_bins_1$sum_hemicellulase==0 & seq_cazy_abun_in_bins_1$sum_oligosaccharidase==0,]

#save.image(file = "cazy_abun_in_bins.RData")

write.table(exo_endo_harboring_bins,file = "seq_cazy_in_exo_endo_harboring_bins.tab",sep="\t",row.names = F,col.names = T)
write.table(only_exo_harboring_bins,file = "seq_cazy_in_only_exo_harboring_bins.tab",sep="\t",row.names = F,col.names = T)
write.table(only_endo_harboring_bins,file = "seq_cazy_in_only_endo_harboring_bins.tab",sep="\t",row.names = F,col.names = T)
write.table(only_hemicellulase_harboring_bins,file = "only_hemicellulase_harboring_bins.tab",sep="\t",row.names = F,col.names = T)
write.table(only_oligosaccharidase_harboring_bins,file = "seq_cazy_in_only_oligosaccharidase_harboring_bins.tab",sep="\t",row.names = F,col.names = T)
write.table(not_carbohydrates_degrading_bins,file = "seq_cazy_not_carbohydrates_degrading_bins.tab",sep="\t",row.names = F,col.names = T)

#---------------------------------------------------------------------------------------------------------------------------
##################################################################################################################################################
###################################################################################################################################################
#################################################################################################################################################
###############################################################################################################################################################

# read the files containg information on the CAZy modules annotated along genes in each draft genome 
cazy_in_bins_dataframes_list<-list()
for (i in 1:length(cazy_in_bins))
{
  tmp<-read.delim(file=cazy_in_bins[i],sep='\t',header = F,stringsAsFactors = F)
  colnames(tmp)<-c("Family_HMM","HMM_length","Query_ID", "Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  tmp$Family_HMM<-gsub(".hmm",'',tmp$Family_HMM)
  bin_name<-strsplit(cazy_in_bins[i],".faa.out.dm.ps")[[1]][1]
  cazy_in_bins_dataframes_list[[bin_name]]<-tmp
}

#--------------------------------------------------------------------------------------------
# further categorization on those genomes with both exoglucanase GH modules and the endoglucanase GH modules
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
    # A attached cellulosome coh>=3, slh>=1,dok==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold"]+1
    }
    ## A attached cellulosome:     coh>=3, dok>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1)
    {
      Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_exo_endo[Abundance_of_gene_pattern_exo_endo$bin %in% exo_endo_harboring_bins_names[j],"attach_scaffold"]+1
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
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group I B"
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group I C"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase==0]<-"Group I D"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & solitude_cellulase>=1]<-"Group I D"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group I E"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase==0]<-"Group I E"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM==0 & solitude_cellulase>=1]<-"Group I F"
})
colnames(Abundance_of_gene_pattern_exo_endo)<-gsub("bin","Draft_genome_bin",colnames(Abundance_of_gene_pattern_exo_endo))
write.table(Abundance_of_gene_pattern_exo_endo,file="Categorization_and_the_abundance_of_gene_patterns_in_exo_endo_harboring_bins.tab",sep='\t',row.names = F,col.names = T)

#-------------------------------------------------------------------------
#######further categorization on those genomes with only exoglucanase GH modules
# further categorization on those genomes with both exoglucanase GH modules and the endoglucanase GH modules
only_exo<-rownames(only_exo_harboring_bins)
#> only_exo
#character(0)

cazy_in_only_exo_harboring_bins<-list()
for (i in 1:length(only_exo))
{
  tmp<-read.delim(file=paste(only_exo[i], ".faa.out.dm.ps", sep=""), sep='\t',header = F,stringsAsFactors = F)
  colnames(tmp)<-c("Family_HMM","HMM_length","Query_ID", "Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  tmp$Family_HMM<-gsub(".hmm",'',tmp$Family_HMM)
  bin_name<-only_exo[i]
  cazy_in_only_exo_harboring_bins[[bin_name]]<-tmp
}


#### identify gene patterns in each draft genomes and calculate the abundance of each pattern identified

only_exo_harboring_bins_names<-names(cazy_in_only_exo_harboring_bins)

Abundance_of_gene_pattern_only_exo<-data.frame("bin"=only_exo_harboring_bins_names)
#A
Abundance_of_gene_pattern_only_exo$attach_scaffold<-c(0)
#B
Abundance_of_gene_pattern_only_exo$free_scaffold<-c(0)
#C
Abundance_of_gene_pattern_only_exo$SLH_CBM<-c(0)
#D
Abundance_of_gene_pattern_only_exo$cellulase_CBM<-c(0)
#
Abundance_of_gene_pattern_only_exo$exo_CBM<-c(0)
#E
Abundance_of_gene_pattern_only_exo$solitude_cellulase<-c(0)
#F
Abundance_of_gene_pattern_only_exo$oligosaccharidase<-c(0)

#str(Abundance_of_gene_pattern_only_exo)
####
Abundance_of_gene_pattern_only_exo$SLH_cohesin<-c(0)
Abundance_of_gene_pattern_only_exo$exo_dok<-c(0)
Abundance_of_gene_pattern_only_exo$endo_dok<-c(0)
Abundance_of_gene_pattern_only_exo$hemicellulase_dok<-c(0)
Abundance_of_gene_pattern_only_exo$oligosaccharidase_dok<-c(0)


#####

for (j in 1: length(only_exo_harboring_bins_names))
{
  tmp0<-cazy_in_only_exo_harboring_bins[[j]]
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
    # A attached cellulosome coh>=3, slh>=1,dok==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"attach_scaffold"]+1
    }
    ## A attached cellulosome:     coh>=3, dok>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"attach_scaffold"]+1
    }
    # B free cellulase:     coh>=3, dok==0, slh==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']==0)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"free_scaffold"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"free_scaffold"]+1
    }
    # C slh_cbm(_gh):      slh>=1, sum(selcbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"SLH_CBM"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"SLH_CBM"]+1
    }
    # D cellulase_selcbm(_gh):     sum(slegh>=1, sum(selcbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"cellulase_CBM"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"cellulase_CBM"]+1
    }
    # E solitude_cellulase(_slh):  sum(slegh>=1, sum(selcbm)==0
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])==0)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"solitude_cellulase"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"solitude_cellulase"]+1 
    }
    #exo-CBM
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"exo_CBM"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"exo_CBM"]+1 
    }
    
    #oligosaccharidase
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"oligosaccharidase"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"oligosaccharidase"]+1 
    }
    ######
    #SLH_cohesin
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=1)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"SLH_cohesin"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"SLH_cohesin"]+1 
    }
    #exo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"exo_dok"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"exo_dok"]+1 
    }
    #endo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% endo])>=1)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"endo_dok"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"endo_dok"]+1 
    }
    #hemicellulase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% hemicellulase])>=1)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"hemicellulase_dok"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"hemicellulase_dok"]+1 
    }
    #oligosaccharidase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"oligosaccharidase_dok"]=Abundance_of_gene_pattern_only_exo[Abundance_of_gene_pattern_only_exo$bin %in% only_exo_harboring_bins_names[j],"oligosaccharidase_dok"]+1 
    }
  }
}

##### assign each draft genome to a category based on the abundance of gene patterns identified in it
Abundance_of_gene_pattern_only_exo<-within(Abundance_of_gene_pattern_only_exo,{
  category<-NA
  category[attach_scaffold>=1]<-"GroupII A"
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"GroupII B"
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"GroupII C"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase==0]<-"GroupII D"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & solitude_cellulase>=1]<-"GroupII D"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"GroupII E"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase==0]<-"GroupII E"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM==0 & solitude_cellulase>=1]<-"GroupII F"
})
colnames(Abundance_of_gene_pattern_only_exo)<-gsub("bin","Draft_genome_bin",colnames(Abundance_of_gene_pattern_only_exo))
write.table(Abundance_of_gene_pattern_only_exo,file="Categorization_and_the_abundance_of_gene_patterns_in_only_exo_harboring_bins.tab",sep='\t',row.names = F,col.names = T)



#-------------------------------------------------------------------------
#######further categorization on those genomes with only endoglucanase GH modules
# further categorization on those genomes with both exoglucanase GH modules and the endoglucanase GH modules
only_endo<-rownames(only_endo_harboring_bins)

cazy_in_only_endo_harboring_bins<-list()
for (i in 1:length(only_endo))
{
  tmp<-read.delim(file=paste(only_endo[i], ".faa.out.dm.ps", sep=""), sep='\t',header = F,stringsAsFactors = F)
  colnames(tmp)<-c("Family_HMM","HMM_length","Query_ID", "Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  tmp$Family_HMM<-gsub(".hmm",'',tmp$Family_HMM)
  bin_name<-only_endo[i]
  cazy_in_only_endo_harboring_bins[[bin_name]]<-tmp
}


#### identify gene patterns in each draft genomes and calculate the abundance of each pattern identified

only_endo_harboring_bins_names<-names(cazy_in_only_endo_harboring_bins)

Abundance_of_gene_pattern_only_endo<-data.frame("bin"=only_endo_harboring_bins_names)
#A
Abundance_of_gene_pattern_only_endo$attach_scaffold<-c(0)
#B
Abundance_of_gene_pattern_only_endo$free_scaffold<-c(0)
#C
Abundance_of_gene_pattern_only_endo$SLH_CBM<-c(0)
#D
Abundance_of_gene_pattern_only_endo$cellulase_CBM<-c(0)
#
Abundance_of_gene_pattern_only_endo$exo_CBM<-c(0)
#E
Abundance_of_gene_pattern_only_endo$solitude_cellulase<-c(0)
#F
Abundance_of_gene_pattern_only_endo$oligosaccharidase<-c(0)

#str(Abundance_of_gene_pattern_only_endo)
####
Abundance_of_gene_pattern_only_endo$SLH_cohesin<-c(0)
Abundance_of_gene_pattern_only_endo$exo_dok<-c(0)
Abundance_of_gene_pattern_only_endo$endo_dok<-c(0)
Abundance_of_gene_pattern_only_endo$hemicellulase_dok<-c(0)
Abundance_of_gene_pattern_only_endo$oligosaccharidase_dok<-c(0)


#####

for (j in 1: length(only_endo_harboring_bins_names))
{
  tmp0<-cazy_in_only_endo_harboring_bins[[j]]
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
    # A attached cellulosome coh>=3, slh>=1,dok==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"attach_scaffold"]+1
    }
    ## A attached cellulosome:     coh>=3, dok>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"attach_scaffold"]+1
    }
    # B free cellulase:     coh>=3, dok==0, slh==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']==0)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"free_scaffold"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"free_scaffold"]+1
    }
    # C slh_cbm(_gh):      slh>=1, sum(selcbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"SLH_CBM"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"SLH_CBM"]+1
    }
    # D cellulase_selcbm(_gh):     sum(slegh>=1, sum(selcbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"cellulase_CBM"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"cellulase_CBM"]+1
    }
    # E solitude_cellulase(_slh):  sum(slegh>=1, sum(selcbm)==0
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])==0)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"solitude_cellulase"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"solitude_cellulase"]+1 
    }
    #exo-CBM
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"exo_CBM"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"exo_CBM"]+1 
    }
    
    #oligosaccharidase
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"oligosaccharidase"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"oligosaccharidase"]+1 
    }
    ######
    #SLH_cohesin
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=1)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"SLH_cohesin"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"SLH_cohesin"]+1 
    }
    #exo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"exo_dok"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"exo_dok"]+1 
    }
    #endo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% endo])>=1)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"endo_dok"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"endo_dok"]+1 
    }
    #hemicellulase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% hemicellulase])>=1)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"hemicellulase_dok"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"hemicellulase_dok"]+1 
    }
    #oligosaccharidase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"oligosaccharidase_dok"]=Abundance_of_gene_pattern_only_endo[Abundance_of_gene_pattern_only_endo$bin %in% only_endo_harboring_bins_names[j],"oligosaccharidase_dok"]+1 
    }
  }
}

##### assign each draft genome to a category based on the abundance of gene patterns identified in it
Abundance_of_gene_pattern_only_endo<-within(Abundance_of_gene_pattern_only_endo,{
  category<-NA
  category[attach_scaffold>=1]<-"GroupIII A"
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"GroupIII B"
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"GroupIII C"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase==0]<-"GroupIII D"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & solitude_cellulase>=1]<-"GroupIII D"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"GroupIII E"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase==0]<-"GroupIII E"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM==0 & solitude_cellulase>=1]<-"GroupIII F"
})
colnames(Abundance_of_gene_pattern_only_endo)<-gsub("bin","Draft_genome_bin",colnames(Abundance_of_gene_pattern_only_endo))
write.table(Abundance_of_gene_pattern_only_endo,file="Categorization_and_the_abundance_of_gene_patterns_in_only_endo_harboring_bins.tab",sep='\t',row.names = F,col.names = T)


#------------------------------------------------------------------------------
# further categorization on those genomes with only_hemicellulase GH modules and the endoglucanase GH modules
only_hemicellulase<-rownames(only_hemicellulase_harboring_bins)

cazy_in_only_hemicellulase_harboring_bins<-list()
for (i in 1:length(only_hemicellulase))
{
  tmp<-read.delim(file=paste(only_hemicellulase[i], ".faa.out.dm.ps", sep=""), sep='\t',header = F,stringsAsFactors = F)
  colnames(tmp)<-c("Family_HMM","HMM_length","Query_ID", "Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  tmp$Family_HMM<-gsub(".hmm",'',tmp$Family_HMM)
  bin_name<-only_hemicellulase[i]
  cazy_in_only_hemicellulase_harboring_bins[[bin_name]]<-tmp
}


#### identify gene patterns in each draft genomes and calculate the abundance of each pattern identified

only_hemicellulase_harboring_bins_names<-names(cazy_in_only_hemicellulase_harboring_bins)

Abundance_of_gene_pattern_only_hemicellulase<-data.frame("bin"=only_hemicellulase_harboring_bins_names)
#A
Abundance_of_gene_pattern_only_hemicellulase$attach_scaffold<-c(0)
#B
Abundance_of_gene_pattern_only_hemicellulase$free_scaffold<-c(0)
#C
Abundance_of_gene_pattern_only_hemicellulase$SLH_CBM<-c(0)
#D
Abundance_of_gene_pattern_only_hemicellulase$cellulase_CBM<-c(0)
#
Abundance_of_gene_pattern_only_hemicellulase$exo_CBM<-c(0)
#E
Abundance_of_gene_pattern_only_hemicellulase$solitude_cellulase<-c(0)
#F
Abundance_of_gene_pattern_only_hemicellulase$oligosaccharidase<-c(0)

#str(Abundance_of_gene_pattern_only_hemicellulase)
####
Abundance_of_gene_pattern_only_hemicellulase$SLH_cohesin<-c(0)
Abundance_of_gene_pattern_only_hemicellulase$exo_dok<-c(0)
Abundance_of_gene_pattern_only_hemicellulase$endo_dok<-c(0)
Abundance_of_gene_pattern_only_hemicellulase$hemicellulase_dok<-c(0)
Abundance_of_gene_pattern_only_hemicellulase$oligosaccharidase_dok<-c(0)


#####

for (j in 1: length(only_hemicellulase_harboring_bins_names))
{
  tmp0<-cazy_in_only_hemicellulase_harboring_bins[[j]]
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
    # A attached cellulosome coh>=3, slh>=1,dok==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"attach_scaffold"]+1
    }
    ## A attached cellulosome:     coh>=3, dok>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"attach_scaffold"]+1
    }
    # B free cellulase:     coh>=3, dok==0, slh==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']==0)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"free_scaffold"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"free_scaffold"]+1
    }
    # C slh_cbm(_gh):      slh>=1, sum(selcbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"SLH_CBM"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"SLH_CBM"]+1
    }
    # D cellulase_selcbm(_gh):     sum(slegh>=1, sum(selcbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"cellulase_CBM"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"cellulase_CBM"]+1
    }
    # E solitude_cellulase(_slh):  sum(slegh>=1, sum(selcbm)==0
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])==0)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"solitude_cellulase"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"solitude_cellulase"]+1 
    }
    #exo-CBM
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"exo_CBM"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"exo_CBM"]+1 
    }
    
    #oligosaccharidase
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"oligosaccharidase"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"oligosaccharidase"]+1 
    }
    ######
    #SLH_cohesin
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=1)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"SLH_cohesin"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"SLH_cohesin"]+1 
    }
    #exo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"exo_dok"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"exo_dok"]+1 
    }
    #endo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% endo])>=1)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"endo_dok"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"endo_dok"]+1 
    }
    #hemicellulase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% hemicellulase])>=1)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"hemicellulase_dok"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"hemicellulase_dok"]+1 
    }
    #oligosaccharidase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"oligosaccharidase_dok"]=Abundance_of_gene_pattern_only_hemicellulase[Abundance_of_gene_pattern_only_hemicellulase$bin %in% only_hemicellulase_harboring_bins_names[j],"oligosaccharidase_dok"]+1 
    }
  }
}
##### assign each draft genome to a category based on the abundance of gene patterns identified in it
Abundance_of_gene_pattern_only_hemicellulase<-within(Abundance_of_gene_pattern_only_hemicellulase,{
  category<-NA
  category[attach_scaffold>=1]<-"Group IV A"
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group IV B"
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group IV "
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase==0]<-"Group IV"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & solitude_cellulase>=1]<-"Group IV"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group IV"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase==0]<-"Group IV"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM==0 & solitude_cellulase>=1]<-"Group IV"
})
colnames(Abundance_of_gene_pattern_only_hemicellulase)<-gsub("bin","Draft_genome_bin",colnames(Abundance_of_gene_pattern_only_hemicellulase))
write.table(Abundance_of_gene_pattern_only_hemicellulase,file="Categorization_and_the_abundance_of_gene_patterns_in_only_hemicellulase_harboring_bins.tab",sep='\t',row.names = F,col.names = T)

#---------------------------------------------------------------------------------
# further categorization on those genomes with only_oligosaccharidase GH modules and the endoglucanase GH modules
only_oligosaccharidase<-rownames(only_oligosaccharidase_harboring_bins)

cazy_in_only_oligosaccharidase_harboring_bins<-list()
for (i in 1:length(only_oligosaccharidase))
{
  tmp<-read.delim(file=paste(only_oligosaccharidase[i], ".faa.out.dm.ps", sep=""), sep='\t',header = F,stringsAsFactors = F)
  colnames(tmp)<-c("Family_HMM","HMM_length","Query_ID", "Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  tmp$Family_HMM<-gsub(".hmm",'',tmp$Family_HMM)
  bin_name<-only_oligosaccharidase[i]
  cazy_in_only_oligosaccharidase_harboring_bins[[bin_name]]<-tmp
}


#### identify gene patterns in each draft genomes and calculate the abundance of each pattern identified

only_oligosaccharidase_harboring_bins_names<-names(cazy_in_only_oligosaccharidase_harboring_bins)

Abundance_of_gene_pattern_only_oligosaccharidase<-data.frame("bin"=only_oligosaccharidase_harboring_bins_names)
#A
Abundance_of_gene_pattern_only_oligosaccharidase$attach_scaffold<-c(0)
#B
Abundance_of_gene_pattern_only_oligosaccharidase$free_scaffold<-c(0)
#C
Abundance_of_gene_pattern_only_oligosaccharidase$SLH_CBM<-c(0)
#D
Abundance_of_gene_pattern_only_oligosaccharidase$cellulase_CBM<-c(0)
#
Abundance_of_gene_pattern_only_oligosaccharidase$exo_CBM<-c(0)
#E
Abundance_of_gene_pattern_only_oligosaccharidase$solitude_cellulase<-c(0)
#F
Abundance_of_gene_pattern_only_oligosaccharidase$oligosaccharidase<-c(0)

#str(Abundance_of_gene_pattern_only_oligosaccharidase)
####
Abundance_of_gene_pattern_only_oligosaccharidase$SLH_cohesin<-c(0)
Abundance_of_gene_pattern_only_oligosaccharidase$exo_dok<-c(0)
Abundance_of_gene_pattern_only_oligosaccharidase$endo_dok<-c(0)
Abundance_of_gene_pattern_only_oligosaccharidase$hemicellulase_dok<-c(0)
Abundance_of_gene_pattern_only_oligosaccharidase$oligosaccharidase_dok<-c(0)


#####

for (j in 1: length(only_oligosaccharidase_harboring_bins_names))
{
  tmp0<-cazy_in_only_oligosaccharidase_harboring_bins[[j]]
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
    # A attached cellulosome coh>=3, slh>=1,dok==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"attach_scaffold"]+1
    }
    ## A attached cellulosome:     coh>=3, dok>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"attach_scaffold"]+1
    }
    # B free cellulase:     coh>=3, dok==0, slh==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']==0)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"free_scaffold"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"free_scaffold"]+1
    }
    # C slh_cbm(_gh):      slh>=1, sum(selcbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"SLH_CBM"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"SLH_CBM"]+1
    }
    # D cellulase_selcbm(_gh):     sum(slegh>=1, sum(selcbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"cellulase_CBM"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"cellulase_CBM"]+1
    }
    # E solitude_cellulase(_slh):  sum(slegh>=1, sum(selcbm)==0
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])==0)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"solitude_cellulase"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"solitude_cellulase"]+1 
    }
    #exo-CBM
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"exo_CBM"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"exo_CBM"]+1 
    }
    
    #oligosaccharidase
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"oligosaccharidase"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"oligosaccharidase"]+1 
    }
    ######
    #SLH_cohesin
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=1)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"SLH_cohesin"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"SLH_cohesin"]+1 
    }
    #exo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"exo_dok"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"exo_dok"]+1 
    }
    #endo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% endo])>=1)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"endo_dok"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"endo_dok"]+1 
    }
    #hemicellulase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% hemicellulase])>=1)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"hemicellulase_dok"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"hemicellulase_dok"]+1 
    }
    #oligosaccharidase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"oligosaccharidase_dok"]=Abundance_of_gene_pattern_only_oligosaccharidase[Abundance_of_gene_pattern_only_oligosaccharidase$bin %in% only_oligosaccharidase_harboring_bins_names[j],"oligosaccharidase_dok"]+1 
    }
  }
}
##### assign each draft genome to a category based on the abundance of gene patterns identified in it
Abundance_of_gene_pattern_only_oligosaccharidase<-within(Abundance_of_gene_pattern_only_oligosaccharidase,{
  category<-NA
  category[attach_scaffold>=1]<-"Group IV"
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group IV"
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group IV"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase==0]<-"Group IV"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & solitude_cellulase>=1]<-"Group IV"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"Group IV"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase==0]<-"Group IV"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM==0 & solitude_cellulase>=1]<-"Group IV"
})
colnames(Abundance_of_gene_pattern_only_oligosaccharidase)<-gsub("bin","Draft_genome_bin",colnames(Abundance_of_gene_pattern_only_oligosaccharidase))
write.table(Abundance_of_gene_pattern_only_oligosaccharidase,file="Categorization_and_the_abundance_of_gene_patterns_in_only_oligosaccharidase_harboring_bins.tab",sep='\t',row.names = F,col.names = T)

#----------------------------------------------------------------------------
no_carbohydrase<-rownames(no_carbohydrase_harboring_bins)

cazy_in_no_carbohydrase_harboring_bins<-list()
for (i in 1:length(no_carbohydrase))
{
  tmp<-read.delim(file=paste(no_carbohydrase[i], ".faa.out.dm.ps", sep=""), sep='\t',header = F,stringsAsFactors = F)
  colnames(tmp)<-c("Family_HMM","HMM_length","Query_ID", "Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  tmp$Family_HMM<-gsub(".hmm",'',tmp$Family_HMM)
  bin_name<-no_carbohydrase[i]
  cazy_in_no_carbohydrase_harboring_bins[[bin_name]]<-tmp
}


#### identify gene patterns in each draft genomes and calculate the abundance of each pattern identified

no_carbohydrase_harboring_bins_names<-names(cazy_in_no_carbohydrase_harboring_bins)

Abundance_of_gene_pattern_no_carbohydrase<-data.frame("bin"=no_carbohydrase_harboring_bins_names)
#A
Abundance_of_gene_pattern_no_carbohydrase$attach_scaffold<-c(0)
#B
Abundance_of_gene_pattern_no_carbohydrase$free_scaffold<-c(0)
#C
Abundance_of_gene_pattern_no_carbohydrase$SLH_CBM<-c(0)
#D
Abundance_of_gene_pattern_no_carbohydrase$cellulase_CBM<-c(0)
#
Abundance_of_gene_pattern_no_carbohydrase$exo_CBM<-c(0)
#E
Abundance_of_gene_pattern_no_carbohydrase$solitude_cellulase<-c(0)
#F
Abundance_of_gene_pattern_no_carbohydrase$oligosaccharidase<-c(0)

#str(Abundance_of_gene_pattern_no_carbohydrase)
####
Abundance_of_gene_pattern_no_carbohydrase$SLH_cohesin<-c(0)
Abundance_of_gene_pattern_no_carbohydrase$exo_dok<-c(0)
Abundance_of_gene_pattern_no_carbohydrase$endo_dok<-c(0)
Abundance_of_gene_pattern_no_carbohydrase$hemicellulase_dok<-c(0)
Abundance_of_gene_pattern_no_carbohydrase$oligosaccharidase_dok<-c(0)


#####

for (j in 1: length(no_carbohydrase_harboring_bins_names))
{
  tmp0<-cazy_in_no_carbohydrase_harboring_bins[[j]]
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
    # A attached cellulosome coh>=3, slh>=1,dok==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"attach_scaffold"]+1
    }
    ## A attached cellulosome:     coh>=3, dok>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"attach_scaffold"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"attach_scaffold"]+1
    }
    # B free cellulase:     coh>=3, dok==0, slh==0
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=3 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']==0 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']==0)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"free_scaffold"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"free_scaffold"]+1
    }
    # C slh_cbm(_gh):      slh>=1, sum(selcbm)>=1
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"SLH_CBM"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"SLH_CBM"]+1
    }
    # D cellulase_selcbm(_gh):     sum(slegh>=1, sum(selcbm)>=1
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"cellulase_CBM"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"cellulase_CBM"]+1
    }
    # E solitude_cellulase(_slh):  sum(slegh>=1, sum(selcbm)==0
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulase])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])==0)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"solitude_cellulase"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"solitude_cellulase"]+1 
    }
    #exo-CBM
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% cellulose_binding_CBM])>=1)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"exo_CBM"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"exo_CBM"]+1 
    }
    
    #oligosaccharidase
    if (sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"oligosaccharidase"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"oligosaccharidase"]+1 
    }
    ######
    #SLH_cohesin
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'SLH','Freq']>=1 & tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'cohesin','Freq']>=1)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"SLH_cohesin"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"SLH_cohesin"]+1 
    }
    #exo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% exo])>=1)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"exo_dok"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"exo_dok"]+1 
    }
    #endo_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% endo])>=1)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"endo_dok"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"endo_dok"]+1 
    }
    #hemicellulase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% hemicellulase])>=1)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"hemicellulase_dok"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"hemicellulase_dok"]+1 
    }
    #oligosaccharidase_dok
    if (tmp_cazy_df_expan[tmp_cazy_df_expan$tmp_cazy %in% 'dockerin','Freq']>=1 & sum(tmp_cazy_df_expan$Freq[tmp_cazy_df_expan$tmp_cazy %in% oligosaccharidase])>=1)
    {
      Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"oligosaccharidase_dok"]=Abundance_of_gene_pattern_no_carbohydrase[Abundance_of_gene_pattern_no_carbohydrase$bin %in% no_carbohydrase_harboring_bins_names[j],"oligosaccharidase_dok"]+1 
    }
  }
}
##### assign each draft genome to a category based on the abundance of gene patterns identified in it
Abundance_of_gene_pattern_no_carbohydrase<-within(Abundance_of_gene_pattern_no_carbohydrase,{
  category<-NA
  category[attach_scaffold>=1]<-c("others")
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"others"
  category[attach_scaffold==0 & free_scaffold>=1 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"others"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & cellulase_CBM>=1 & solitude_cellulase==0]<-"others"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM>=1 & solitude_cellulase>=1]<-"others"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase>=1]<-"others"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM>=1 & solitude_cellulase==0]<-"others"
  category[attach_scaffold==0 & free_scaffold==0 & SLH_CBM==0 & cellulase_CBM==0 & solitude_cellulase>=1]<-"others"
})
colnames(Abundance_of_gene_pattern_no_carbohydrase)<-gsub("bin","Draft_genome_bin",colnames(Abundance_of_gene_pattern_no_carbohydrase))
write.table(Abundance_of_gene_pattern_no_carbohydrase,file="Categorization_and_the_abundance_of_gene_patterns_in_no_carbohydrase_harboring_bins.tab",sep='\t',row.names = F,col.names = T)

#-------------------------------------------------------------------------------------
results_all_bins<-rbind(Abundance_of_gene_pattern_exo_endo,Abundance_of_gene_pattern_only_exo,Abundance_of_gene_pattern_only_endo,Abundance_of_gene_pattern_only_hemicellulase,Abundance_of_gene_pattern_only_oligosaccharidase,Abundance_of_gene_pattern_no_carbohydrase)
write.table(results_all_bins,file="Categorization_and_the_abundance_of_gene_patterns_in_all_bins.tab",sep='\t',row.names = F,col.names = T)

save.image(file="gene_pattern_identification_and_categorization_of_draft_genomes.RData")
