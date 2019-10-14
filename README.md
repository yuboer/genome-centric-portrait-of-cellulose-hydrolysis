# genome-centric-portrait-of-microbes-cellulose-hydrolysing-capacity
# This pipeline is developed to pair the dbCAN annotation platform to interpret MAGs on the specific function niche of cellulose hydrolysis


-------------------------------------------------------------------------------------------------------------------------------------


# part I: pre-annotation through the dbCAN platform (in batch mode)

# database and hmm tools need o be downloaded
# 1.The database "dbCAN-HMMdb-V8.zip" listed above is downloaded from http://bcb.unl.edu/dbCAN2/download/Databases/,  to annotate the CAZy modules in the MAGs or complete genomes
   
# The "hmmscan-parser.gz" is downloaded from http://bcb.unl.edu/dbCAN2/download/Tools/, "hmmscan-parser.sh" will be used to parse the #HMM annotation results

# Below are command lines for batch annotation of the MAGs throug the dbCAN hmmsearch in linux system

# go to the directory in which locates the dbCAN database, in my case, it is as below, while the user may need to adjust the path #accordingly
cd /home/ywx1845/software/dbCAN_db/
# format the database for hmmscan
hmmpress dbCAN-HMMdb-V8.txt

# go to the folder "MAGs_faa", in which are amino acid sequences of the MAGs (with the suffix .faa) to be annotated, in my case, it is as #below
cd /projects/b1052/Wells_b1042/yuboer/MAGs_faa/

# hmmscan; /home/ywx1845/software/dbCAN_db/ is the path of the directory in which locates the formated dbCAN database, the user may need to modify this path
find . -name "*.faa" | while read line ; do hmmscan --domtblout ${line}.out.dm /home/ywx1845/software/dbCAN_db/dbCAN-HMMdb-V8.txt $line > ${line}.out; done

# parse the hmmscan results
find . -name "*.out.dm"|while read line ; do sh /home/ywx1845/software/dbCAN/hmmscan-parser.sh $line > ${line}.ps; done

# list the names of MAGs with CAZy modules identified
find . -name "*.out.dm.ps" ! -size 0 > filteredhmmout.list

# make a directory "dbCAN_annotation_results" and move the parsed hmmscan results to this directory
mkdir ../dbCAN_annotation_results

mv *.out.dm.ps > ../dbCAN_annotation_results

# delete the empty files in the folder of "dbCAN_annotation_results"

# the non-empty files in the directory "dbCAN_annotation_results" are input files to the following R scripts

-----------------------------------------------------------------------------------------------------------------------------------


# part II: genome-centric portrait on the corresponding microbes cellulolytic competency

# Files in the folder of ./dbCAN_annotation_results/ were the input file to the annotation pipeline developed in this study,

# Genome categorization: run the R script named as "genome_categorization_on_cellulolytic_competency.R"
# Visualization: CAZy arrangment along genes in each of the MAgs, using the R script names as "genoplot_CAZy_in_bins.R"

Notes: 
a.) Genome categorization by "genome_categorization_on_cellulolytic_competency.R" will generate 3 files summarizing  detailed information on: 

    1) the diversity and abundances of all the CAZy modules identified in each of the MAG or complete genomes, the output file is named as 'seq_cazy_abun_in_bins.csv'
    
    
    2) preliminarily assign the MAGs into 6 groups, namely:
       
       group I - "exo_endo_harboring_bins", the MAGs harboring both the exoglucanase GH modules and the endoglucanase GH modules
       
       group II - "only_exo_harboring_bins", the MAGs harboring only the exoglucanase GH modules
       
       group III - "only_endo_harboring_bins", the MAGs harboring only the endoglucanase GH modules
       
       group IV - "only_hemicellulase_harboring_bins", the MAGs harboring only the hemicellulase GH modules
       
       group V - "only_oligosaccharidase_harboring_bins", the MAGs harboring only the oligosaccharides GH modules
       
       Others or group VI - "not_carbohydrates_degrading_bins", the MAGs harboring none of the above mentioned GH modules
    
    
    3) further categorization of the MAGs harboring both the exoglucanase and endoglucanase GH modules, these MAGs were preliminarily assigned in group I as introduced above - "exo_endo_harboring_bins"; these Group I MAGs would be further categorized into 6 subgroups, namely: 'group I-a', group 'I-b', 'group I-c', 'group I-d', 'group I-e' and 'group I-f', respectively
        
        'group I-a' ---- genomes of model cellulolytic microbes,  cellulosome gene clusters are identified, and these cellulsome gene  clusters could convery both the assembly of the diverse catalytic components and the CEM complex;
        
        'group I-b' ---- genomes of proficient cellulolytic microbes,  cellulosome gene clusters are identified, however, the scaffold genes in genomes of this subgroup are in lack of the SLH module. the corresponding cellulsome are not cell surface ahering, and it is that the enzymes encoded by the cellulosome indpendent SH-CBM genes that facilitate the microbe-cellulose contact and the CEM complex;                      
        'group I-c' ---- genomes identified with cellulosome gene clusters, however, the scaffold genes in genomes of this subgroup are in lack of the SLH module, and unlike the genomes in Group I-b, no cellulosome-independent SH-CBM genes are identified in genomes of Group I-c; the annotation of potential novel cellulose-microbes adhesion mechineries genomes is important in predicting whether the corresponding microbes are cellulolytic or not;
        
        'group I-d' ---- genomes of cellulolytic microbes without cellulosome gene clusters, the microbe-cellulose contact and the CEM   complex are much likely conveyed through enzymes encoded by the cellulosome-independent SLH-CBM genes;
        
        'group I-e' ---- genomes with neither cellulosome gene clusters nor the SLH-CBM genes, genes harboring both the cellulase GH     modules and the the cellulose binding CBM mocules are identified; uncertainty are encountered in identifying the cellulolytic capacity of genomes; the annotation of potential novel cellulose-microbes adhesion mechineries genomes is important in predicting whether the corresponding microbes are cellulolytic or not;
        
        'group I-f' ---- genomes of non-cellulolytic microbes; although both exoglucanase and endoglucanase GH modules are identified,   none of these GH modules were observed occurring in same genes with the cellulose-binding CBM modules
        
        The output file is the one named as "Abundance_of_gene_pattern_exo_endo.csv"
        
        
b). Visualization of the CAZy modules' arrangement along genes in each MAGs (optional)

    Run the R script named as "genoplot_CAZy_in_bins.R", a folder will be generated holding a corresponding number of PDF files with the     visualization information on the CAZy modules' arrangement along genes in each genome
    
-----------------------------------------------------------------------------------------------------------------------------

# Additional note:

Besides the annotation of the MAGs, this pipeline can aslo be applied to annotate all the amino acid sequences derived from the metagenome datasets, to get info on: 1) the diversity and abundances of all the CAZy modules identified in the community; 2)whether cellulosome gene clusters were present in the community; 3) the types of the cellulosome gene clusters identified (cell surface adhering or not); 3) diversity and abundance of the carbohydrate active genes in the overall community; to achieve this, each amino acid sequences of a dataset is equivalent to one MAG in the above annotation pipeline.

# inquery could be sent to irenelambiel@gmail.com
# Feel free to try with the example dataset
