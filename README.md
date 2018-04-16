# genome-centric-portrait-of-cellulose-hydrolysis
#This pipeline is developed to pair the dbCAN annotation platform to interpret MAGs on the specific function niche of cellulose hydrolysis


-------------------------------------------------------------------------------------------------------------------------------------


partI:

1. Annotate the CAZy modules in the MAGs through the HMM scan of dbCAN, please  also refer to http://csbl.bmb.uga.edu/dbCAN/ for more detailed    documentation
   
   Need to Download :  the database "dbCAN-fam-HMMs.txt" and also the "hmmscan-parser.sh"script from http://csbl.bmb.uga.edu/dbCAN/

  Below is the command line we applied to HMM scan all the MAGs in the folder of ./five_genomes/faa/, all the files in this folder were   amino acid seuquences, and the suffix  was ".faa", command lines under the linux operation system

  #hmmpress dbCAN-fam-HMMs.txt
  
  #find ./five_genomes/faa -name "*.faa" | while read line ; do hmmscan --domtblout ${line}.out.dm dbCAN-fam-HMMs.txt $line >          ${line}.out; done
  
  #find ./five_genomes/faa -name "*.out.dm"|while read line ; do sh hmmscan-parser.sh $line > ${line}.ps; done
  
  #find ./five_genomes/faa -name "*.out.dm.ps" ! -size 0 > filteredhmmout.list
  
  #mkdir dbCAN_annotation_results
  
  #mv *.out.dm.ps > ./dbCAN_annotation_results



-----------------------------------------------------------------------------------------------------------------------------------


partII:

Files in the folder of ./dbCAN_annotation_results/ were the input file to the annotation pipeline developed in this study,

1. run the R script named as "gene_pattern_identification_and_draft_genome_cetegorization.R", which will generate 3 files summarizing  detailed information on: 
    1) the diversity and abundances of all the CAZy modules identified in each of these 47 genomes;
    2) the diversity and abundances of the different types of the carbohydrate active genes in each genome;
    3) assignment of the potential cellulose hydrolyzing genomes into the 6 subgroups of different cellulolytic machinery/capacity. 

2. Visualization of the CAZy modules' arrangement along genes in each MAGs (optional)

    Run the R script named as "genoplot_CAZy_in_bins.R", a folder will be generated holding a corresponding number of PDF files with the     visualization information on the CAZy modules' arrangement along genes in each genome
    
-----------------------------------------------------------------------------------------------------------------------------
