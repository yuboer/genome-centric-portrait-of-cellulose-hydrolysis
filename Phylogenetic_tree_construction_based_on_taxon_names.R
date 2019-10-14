# this script was used to generate the circle tree in the paper
setwd("/workhome/yuboer_CAZy_in_2786complete_genome")
a<-rownames(seq_cazy_temp5)
seq_cazy_temp5$Assembly.Accession<-a
###add taxonn info with the abundance of cazy modules annotated
temp6<-merge(seq_cazy_temp5,temp3,by.x="Assembly.Accession",by.y = "Assembly.Accession",all = T)

temp6_uniq_taxon<-temp6[!duplicated(temp6$Assembly.Accession),]

uniq_taxon_exo_endo_genomes_1<-uniq_taxon_exo_endo_genomes[c("Group.1","category")]

uniq_taxon_only_exo_genomes_1<-uniq_taxon_only_exo_genomes[c("Group.1","category")]
uniq_taxon_only_endo_genomes_1<-uniq_taxon_only_endo_genomes[c("Group.1","category")]

temp7<-merge(sel_taxa_cat,temp6_uniq_taxon,by.x ="Group.1", by.y ="Assembly.Accession", all = T)

taxon_level<-c('superkingdom','phylum','class','order','family','genus','species','strain')
info_sle<-c('Group.1',"category",taxon_level)

temp8<-temp7[,colnames(temp7) %in% info_sle]

superkingdom<-temp8$superkingdom
phylum<-temp8$phylum
class<-temp8$class
order<-temp8$order
family<-temp8$family
genus<-temp8$genus
species<-temp8$species
strain<-temp8$strain


phylum[is.na(phylum)]<-c("unknown.phylum")
class[is.na(class)]<-c("unknown.class")
order[is.na(order)]<-c("unknown.order")
family[is.na(family)]<-c("unknown.family")
genus[is.na(genus)]<-c("unknown.genus")
species[is.na(species)]<-c("unknown.species")
strain[is.na(strain)]<-c("unknown.strain")


superkingdom_1<-superkingdom
phylum_1<-paste(superkingdom,phylum,sep = "_")
class_1<-paste(superkingdom,phylum,class,sep = "_")
order_1<-paste(superkingdom,phylum,class,order,sep = "_")
family_1<-paste(superkingdom,phylum,class,order,family,sep = "_")
genus_1<-paste(superkingdom,phylum,class,order,family,genus,sep = "_")
species_1<-paste(superkingdom,phylum,class,order,family,genus,species,sep = "_")
strain_1<-paste(superkingdom,phylum,class,order,family,genus,species,strain,sep = "_")

all_taxa<-data.frame(superkingdom_1,phylum_1,class_1,order_1,family_1,genus_1,species_1,strain_1)
colnames(all_taxa)<-c('superkingdom','phylum','class','order','family','genus','species','strain')


ploting_taxa<-all_taxa[,-c(7)]

for_tree_df<-matrix(nrow=0,ncol=2)
for (i in 1:6)
{
  tmp<-ploting_taxa[,i:(i+1)]
  tmp<-tmp[!duplicated(tmp[,2]),]
  colnames(tmp)<-c("node1","node2")
  for_tree_df<-rbind(for_tree_df,tmp)
}

dim(for_tree_df)


#install.packages("igraph")
#library(igraph)


cytoscape_for_tree<-graph_from_data_frame(d = for_tree_df, directed = FALSE)

write.graph(cytoscape_for_tree,file="cytoscape_for_tree.graphml",format="graphml") # layout of the graph was modified in cytoscape3


ploted_taxa<-data.frame(superkingdom_1,phylum_1,class_1,order_1,family_1,genus_1,strain_1)

lels<-c('superkingdom','phylum','class','order','family','genus','strain')

node_attr<-matrix(nrow=0,ncol=2)
for (i in 1: length(lels))
  {
  tmp<-ploted_taxa[,i]
  tmp_1<-data.frame(node=unique(tmp),level=lels[i])
  node_attr<-rbind(node_attr,tmp_1)
}

node_attr$node<-as.character(node_attr$node)
node_attr$level<-as.character(node_attr$level)

temp9<-cbind(temp8[,1:2],all_taxa)

temp9$category<-as.character(temp9$category)
temp9$strain<-as.character(temp9$strain)
strains_to_be_colored<-temp9$strain[!is.na(temp9$category)]
strains_to_be_colored<-as.character(strains_to_be_colored)

#length(strains_to_be_colored)
#[1] 2116
for (i in 1:length(strains_to_be_colored))
{
  node_attr[node_attr$node %in% strains_to_be_colored[i],'level']<-unique(temp9[temp9$strain %in% strains_to_be_colored[i],"category"])
}
  
write.table(node_attr,file="node_attr.tab",quote=F,row.name=F,col.name=T,sep="\t")

#?write.table

save.image(file = "phylogenetic_tree_construction.RData")
getwd()


#################
#ploting the tree at 6 lelvels
ploting_taxa<-all_taxa[,-c(6,7)]

#for_tree_df<-data.frame(node1=c(NA),node2=c(NA))
for_tree_df<-matrix(nrow=0,ncol=2)
for (i in 1:5)
{
  tmp<-ploting_taxa[,i:(i+1)]
  tmp<-tmp[!duplicated(tmp[,2]),]
  colnames(tmp)<-c("node1","node2")
  for_tree_df<-rbind(for_tree_df,tmp)
}

dim(for_tree_df)


#install.packages("igraph")
#library(igraph)

cytoscape_for_tree<-graph_from_data_frame(d = for_tree_df, directed = FALSE)
#str(cytoscape_for_tree)
write.graph(cytoscape_for_tree,file="cytoscape_for_tree.graphml",format="graphml") # layout of the graph was modified in cytoscape3

ploted_taxa<-data.frame(superkingdom_1,phylum_1,class_1,order_1,family_1,strain_1)
#colnames(all_taxa)<-c('superkingdom','phylum','class','order','family','genus','species','strain')
lels<-c('superkingdom','phylum','class','order','family','strain')
#lels<-c('superkingdom','phylum','class','order','family','genus','strain')
node_attr<-matrix(nrow=0,ncol=2)
for (i in 1: length(lels))
{
  tmp<-ploted_taxa[,i]
  tmp_1<-data.frame(node=unique(tmp),level=lels[i])
  node_attr<-rbind(node_attr,tmp_1)
}

node_attr$node<-as.character(node_attr$node)
node_attr$level<-as.character(node_attr$level)

temp9<-cbind(temp8[,1:2],all_taxa)

temp9$category<-as.character(temp9$category)
temp9$strain<-as.character(temp9$strain)
strains_to_be_colored<-temp9$strain[!is.na(temp9$category)]
strains_to_be_colored<-as.character(strains_to_be_colored)

#length(strains_to_be_colored)
#[1] 2116
for (i in 1:length(strains_to_be_colored))
{
  node_attr[node_attr$node %in% strains_to_be_colored[i],'level']<-unique(temp9[temp9$strain %in% strains_to_be_colored[i],"category"])
}

write.table(node_attr,file="node_attr.tab",quote=F,row.name=F,col.name=T,sep="\t")
