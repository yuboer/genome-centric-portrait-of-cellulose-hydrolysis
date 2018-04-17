# genome-centric-portrait-of-cellulose-hydrolysis
#This pipeline is developed to pair the dbCAN annotation platform to interpret MAGs on the specific function niche of cellulose hydrolysis


-------------------------------------------------------------------------------------------------------------------------------------


partI: pre-annotation through the dbCAN platform

1. Annotate the CAZy modules in the MAGs through the HMM scan of dbCAN, please  also refer to http://csbl.bmb.uga.edu/dbCAN/ for more detailed    documentation
   
   Need to Download :  the database "dbCAN-fam-HMMs.txt" and also the "hmmscan-parser.sh"script from http://csbl.bmb.uga.edu/dbCAN/

  Below is the command line (linux system) we applied to HMM scan all the MAGs in the folder of ./five_genomes/faa/, all the files in this folder were   amino acid seuquences, and the suffix of each file was ".faa"

  1) hmmpress dbCAN-fam-HMMs.txt
  
  2) find ./five_genomes/faa -name "*.faa" | while read line ; do hmmscan --domtblout ${line}.out.dm dbCAN-fam-HMMs.txt $line >          ${line}.out; done
  
  3) find ./five_genomes/faa -name "*.out.dm"|while read line ; do sh hmmscan-parser.sh $line > ${line}.ps; done
  
  4) find ./five_genomes/faa -name "*.out.dm.ps" ! -size 0 > filteredhmmout.list
  
  5) mkdir dbCAN_annotation_results
  
  6) mv *.out.dm.ps > ./dbCAN_annotation_results



-----------------------------------------------------------------------------------------------------------------------------------


partII: genome-centric portrait pipeline proposed in this study

Files in the folder of ./dbCAN_annotation_results/ were the input file to the annotation pipeline developed in this study,

1. run the R script named as "gene_pattern_identification_and_draft_genome_cetegorization.R", which will generate 3 files summarizing  detailed information on: 


    1) the diversity and abundances of all the CAZy modules identified in each of these 47 genomes, the output file is named as                'seq_cazy_abun_in_bins.tab'
    
    
    2) preliminarily assign the MAGs into 6 groups, namely:
       
       group I - "exo_endo_harboring_bins", the MAGs harboring both the exoglucanase GH modules and the endoglucanase GH modules
       
       group II - "only_exo_harboring_bins", the MAGs harboring only the exoglucanase GH modules
       
       group III - "only_endo_harboring_bins", the MAGs harboring only the endoglucanase GH modules
       
       group IV - "only_hemicellulase_harboring_bins", the MAGs harboring only the hemicellulase GH modules
       
       group V - "only_oligosaccharidase_harboring_bins", the MAGs harboring only the oligosaccharides GH modules
       
       Others or group VI - "not_carbohydrates_degrading_bins", the MAGs harboring none of the above mentioned GH modules
    
    
    3) further categorization of the MAGs preliminarily assigned in group I - "exo_endo_harboring_bins", into 6 subgroup, namely: 'group        I-a', group 'I-b', 'group I-c', 'group I-d', 'group I-e' and 'group I-f', respectively
        
        'group I-a' ---- genomes of model cellulolytic microbes,  cellulosome gene clusters are identified, and these cellulsome gene                            clusters could convery both the assembly of the diverse catalytic components and the CEM complex;
        
        'group I-b' ---- genomes of proficient cellulolytic microbes,  cellulosome gene clusters are identified, these cellulsome gene                            clusters convery only the assembly of the diverse catalytic components,  and it is some cellulosome-                                    indpendent SH-CBM genes that facilitate the microbe-cellulose contact and the CEM complex                       
        
        'group I-c' ---- genomes identified with cellulosome gene clusters, yet are much likely inert in cellulose hydrolysis,                                    cellulosome gene clusters are identified, these cellulsome gene clusters could convery only the assembly of the                          diverse catalytic components,  and no cellulosome-independent SH-CBM genes are identified; genomes would be                              inert in cellulose hydrolysis if no alternative microbe-cellulose contact machineries beyongd the SLH-CBM genes                          would be identified
        
        'group I-d' ---- genomes of cellulolytic microbes without cellulosome gene clusters, the microbe-cellulose contact and the CEM                            complex are much likely conveyed through the cellulosome-independent SLH-CBM genes;
        
        'group I-e' ---- genomes with neither cellulosome gene clusters nor the SLH-CBM genes, genes harboring both the cellulase GH                              modules and the the cellulose binding CBM mocules are identified; uncertainty are encountered in                                        identifying the cellulolytic capacity of genomes in this group
        
        'group I-f' ---- genomes of non-cellulolytic microbes; although both exoglucanase and endoglucanase GH modules are identified,                            none of these GH modules were observed occurring in same genes with the cellulose-binding CBM modules
        
        The output file is the one named as "Categorization_and_the_abundance_of_gene_patterns_in_exo_endo_harboring_bins.tab"
        
        
2. Visualization of the CAZy modules' arrangement along genes in each MAGs (optional)

    Run the R script named as "genoplot_CAZy_in_bins.R", a folder will be generated holding a corresponding number of PDF files with the     visualization information on the CAZy modules' arrangement along genes in each genome
    
-----------------------------------------------------------------------------------------------------------------------------
