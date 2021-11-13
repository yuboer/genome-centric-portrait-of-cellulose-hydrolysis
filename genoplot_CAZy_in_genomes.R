###visualize the cazy module arrangement along each cazy gene annoated in each genome

install.packages("genoPlotR")
library("genoPlotR")

#------------------------
####set the working directory to the folder where files of dbCAN annotation locates
#define the absolute path of the folder holding the dbCAN annotation results; instead of "\" as in the above line, use "/" in R to define the address of the files and folder, as below 

dbCAN_results_data_path = as.character("D:/papers_in_progress/to finish/HKU_cellulolytic_genomes/dbcan_results/")
setwd(dbCAN_results_data_path)
getwd()
#-----------------------
cazy_in_genomes_list_0<-dir(dbCAN_results_data_path)
cazy_in_genomes_list<-sapply(strsplit(cazy_in_genomes_list_0,'.faa.out.dm.ps'),'[',1)
#cazy_in_genomes_list<-gsub(".faa.out.dm.ps","",cazy_in_genomes_list_0)
####read all the files and store them in a list, namely 'cazy_in_genomes_hmm' here:
cazy_in_genomes_hmm<-list()
for (i in 1: length(cazy_in_genomes_list))
{
  temp<-read.delim(file = paste(dbCAN_results_data_path,cazy_in_genomes_list_0[i],sep=""),header = F, sep='\t')
  colnames(temp)<-c("cazy_module","HMM_length","Query_ID","Query_length","E_value","HMM_start","HMM_end","Query_start","Query_end","Coverage")
  temp$cazy_module<-gsub(".hmm","",temp$cazy_module)
  #temp$cazy_module<-sapply(strsplit(temp$cazy_module,".hmm") '[',1)
  temp$strand<-c(1)
  cazy_in_genomes_hmm[[cazy_in_genomes_list[i]]]<-temp
  #cazy_in_genomes_hmm$cazy_in_genomes_list_0[i]<-temp
  #cazy_in_genomes_list_0[i]=cazy_in_genomes_hmm[[i]]
}

####genoplot on arrangment of all the cazy modules along genes in the genomes:
#cazy_in_genomes_hmm_plotfrags<-list()
for (j in 1:length(cazy_in_genomes_list))
{
  genomes_tmp<-cazy_in_genomes_hmm[[j]]
  query_id_in_genomes_tmp<-genomes_tmp[,'Query_ID']
  uni_query_id_in_genomes_tmp<-unique(query_id_in_genomes_tmp)
  uni_query_id_in_genomes_tmp<-as.character(uni_query_id_in_genomes_tmp)
  plotfrag<-list()
  for (i in 1:length(uni_query_id_in_genomes_tmp))
  {
    tmp<-genomes_tmp[genomes_tmp$Query_ID %in% uni_query_id_in_genomes_tmp[i],]
    tmp1<-tmp[c("cazy_module", "Query_start","Query_end","strand")]
    colnames(tmp1)<-c("name","start","end","strand")
    tmp1$col<-c("Sky Blue 1")
    tmp1$fill<-c("Sky Blue 1")
    tmp1$gene_type<-c("blocks")
    plotfrag[[uni_query_id_in_genomes_tmp[i]]]<-as.dna_seg(tmp1)
    #plotfrag[[i]]<-as.dna_seg(tmp1)
  }
  #names(plotfrag)<-uni_query_id_in_genomes_tmp
  allannots <- lapply(plotfrag, function(x)
  {mid <- middle(x)
  annot <- annotation(x1=mid, text=x$name, rot=40, col='Gray 41')})
  #cazy_in_genomes_hmm_plotfrags[[cazy_in_genomes_list_0[j]]]<-plotfrag
  #cazy_in_genomes_hmm[[j]]<-plotfrag
  #pdf(file="test.pdf",height=100,width=20)
  pdf(file=paste(cazy_in_genomes_list[j],".pdf",sep=""),height=130,width=30)
  plot_gene_map(plotfrag,annotations = allannots,annotation_cex =0.8,
                main=paste("Arrangement of the CAZy modules in",cazy_in_genomes_list[j],sep=" "),cex.main=1.2,col.main="Gray 41",
                dna_seg_labels =uni_query_id_in_genomes_tmp, dna_seg_label_cex = 0.8,dna_seg_label_col="Gray 41")
  #plot_gene_map(cazy_in_genomes_hmm_plotfrags[[cazy_in_genomes_list_0[j]]],annotations = allannots,annotation_cex =0.8,dna_seg_label_cex=1)
  #plot_gene_map(cazy_in_genomes_hmm_plotfrags[j],annotations = allannots,annotation_cex =0.8,dna_seg_label_cex=1)
  dev.off()
}
