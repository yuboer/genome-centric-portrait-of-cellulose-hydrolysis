Genome-centric portrait of the anaerobes’ cellulolytic competency
This pipeline is developed to pair the dbCAN annotation platform to interpret genomes of anaerobes on their specific function niche of cellulose hydrolysis

-------------------------------------------------------------------------------------------------------------------------------------

# part I: genome-centric CAZy annotation through the dbCAN platform (in batch mode)

step1. database preparation
dbCAN database and hmm tools need to be downloaded
The database "dbCAN-HMMdb-V8.zip" is downloaded from http://bcb.unl.edu/dbCAN2/download/Databases/  
The "hmmscan-parser.gz" is downloaded from http://bcb.unl.edu/dbCAN2/download/Tools/, "hmmscan-parser.sh" will be used to parse the #HMM annotation results

step2. CAZy annotation
Below are command lines for batch annotation of genomes in linux system
1) go to the directory in which locates the dbCAN database, in my case, it is as below:
cd /home/ywx1845/software/dbCAN_db/
2) format the database for hmmscan
hmmpress dbCAN-HMMdb-V8.txt
3）go to the folder that contains amino acid sequences of genomes (with suffix .faa), in my case, the name of the folder is "MAGs_faa", and the path is as below
cd /projects/b1052/Wells_b1042/yuboer/MAGs_faa/
4）do hmmscan; /home/ywx1845/software/dbCAN_db/ is the absolute path of the formated dbCAN database
find . -name "*.faa" | while read line ; do hmmscan --domtblout ${line}.out.dm /home/ywx1845/software/dbCAN_db/dbCAN-HMMdb-V8.txt $line > ${line}.out; done
5) parse the hmmscan results
find . -name "*.out.dm"|while read line ; do sh /home/ywx1845/software/dbCAN/hmmscan-parser.sh $line > ${line}.ps; done
6) list the names of genomes with CAZy modules identified
find . -name "*.out.dm.ps" ! -size 0 > filteredhmmout.list
7) make a directory and move the parsed hmmscan results to this directory, in my case, the directory name is 'dbCAN_annotation_results'
mkdir ../dbCAN_annotation_results
mv *.out.dm.ps > ../dbCAN_annotation_results
8) delete the empty files in the folder of "dbCAN_annotation_results"
9) the non-empty files in the directory "dbCAN_annotation_results" are input files to the following R scripts

-----------------------------------------------------------------------------------------------------------------------------------

# part II: genome-centric portrait of the corresponding anaerobes' cellulose-hydrolyzing capacity
Files with suffix '.out.dm.ps' in the folder of ./dbCAN_annotation_results/ were input files for the annotation pipeline (embodied in R scripts) developed in this study

stpe3. genome-centric categorization of putative cellulolytic anaerobes
if you are working with complete genomes, in which case, the plasmid and chromsome sequences are of seperate .faa files and the corresponding dbcan annotation results for plasmid and chromsome need to be aggregated, run the R script named as "1complete_genome_annotation.R"
if you are working with high-quality metagenome-assembled-genomes (MAGs), or you are working with complete genomes for which the plasmid and chromsome sequences of #a single strain have been combined into one single file, run the R script named as "MIMAG_annotation.R"

Notes:
There will be three files generated by this step of annotation
1) a files named as 'summary_cazy_in_genomes.csv', which acts as a summary on abundance of each CAZy module (e.g., GH48) annotated, and also abundance of each type of CAZy gene (e.g., SLH-CBM gene) annotated in each genome;
2) a file named as 'GroupI_genome_categorization.csv', which acts as a summary on categorization of each genome preliminarily assigned into Group I (harboring both the exo- and endo- GH modules)
3) a file named as 'GroupII_genome_categorization.csv', which acts as a summary on categorization of each genome preliminarily assigned into Group II (without presence of either the exo- or the endo- GH module)

Putative cellulolytic anaerobes are only to be expected for genomes assigned into Group I-a to Group I-e
1) 'group I-a' ---- genomes of paradigm cellulolytic anaerobes;  cellulosome gene clusters are identified in genomes of this group, and the corresponding cellulosome complexes could convery both the assembly of diverse catalytic components and the cellulose-enzyme-microbe (CEM) complex formation;    
2) 'group I-b' ---- genomes of paradigm cellulolytic anaerobes, cellulosome gene clusters are identified, however, SLH module is absent in these cellulosome gene clusters and the corresponding cellulosome complex may not be cell-surface adhering. The CEM complex formation in anaerobes of this group was likely initiated by several cellulosome indpendent SLH-CBM enzymes rather than by cellulosome complexes
3) 'group I-c' ---- although cellulosome gene clusters were identified in their genomes, the corresponding anaerobes are putatively noncellulolytic.   SLH module is absent in the cellulosome gene clusters, and cellulosome-independent SLH-CBM genes are also absent in the genome. In predicting whether the corresponding anaerobe of this group is cellulolytic or not, uncertainty may arise from either the incompleteness of MAGs or existense of potetnial novel microbe-substrate adhesion mechanisms. If you are working with complete genomes, unless potential novel microbe-substrate adhesion mechanism may exist, complete genomes categorized in this group are much likely noncellulolytic anaerobes.  If you are working with incomplete genomes (e.g., MAGs), genomes categorized in this group might also be of cellulolytic anaerobes 
4) 'group I-d' ---- genomes of cellulolytic anaerobes without cellulosome gene clusters. Microbe-cellulose adhesion and CEM complex formaiton in anaerobes of this group is likely initiated by enzymes encoded by genes harboring both the SLH module and the cellulose-binding CBM module
5) 'group I-e' ---- genomes with neither cellulosome gene clusters nor the SLH-CBM genes, uncertainty are encountered in predicting the cellulolytic capacity of anaerobes in this group; potential novel microbe-cellulose adhesion mechanism awaits to be uncovered for cellulolytic anaerobes in this group;

# part III: Visualizetion of the CAZy arrangment along genes in each genome, in batch mode
apply the R script named as "genoplot_CAZy_in_genomes.R"
   
-----------------------------------------------------------------------------------------------------------------------------

Additional note:

Besides the annotation of genomes, this pipeline can aslo be applied in annotation of metagenome dataset as a whole (in which case, to interpretate each metagenome dtaset as one single faa file), to get info on: 1) abundance of each CAZy modules in the community; 2)whether cellulosome gene clusters are present in the community or not; 3) abundance of each type of CAZy gene in the community

Inquery could be sent to yuboer.abbado@outlook.com
/Feel free to try this pipeline with the example dataset provided in the zip folder names as "example_datasets_and_results.zip"
1) the test dataset is a collection of .faa.out.pm.ds files in the subfolder named as "dbcan_results";
2) resulting files of the genome-wide annotation and categorization are in the subfolder named as "genome_categorization_result_files";
3) resulting files from genoplot visualization are in the subfolder named as "genoplot_results"
4) the two custom scripts applied in this part of analysis are also included in this folder

And feel free to drop a message in case errors are encountered
