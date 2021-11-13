#step1. database preparation
#dbCAN database and hmm tools need to be downloaded
#The database "dbCAN-HMMdb-V8.zip" is downloaded from http://bcb.unl.edu/dbCAN2/download/Databases/  
#The "hmmscan-parser.gz" is downloaded from http://bcb.unl.edu/dbCAN2/download/Tools/, "hmmscan-parser.sh" will be used to parse the #HMM annotation results

#step2. CAZy annotation
#Below are command lines for batch annotation of genomes in linux system
#1) go to the directory in which locates the dbCAN database, in my case, it is as below:
cd /home/ywx1845/software/dbCAN_db/
#2) format the database for hmmscan
hmmpress dbCAN-HMMdb-V8.txt
#3）go to the folder that contains amino acid sequences of genomes (with suffix .faa), in my case, the name of the folder is "MAGs_faa", and the path is as below
cd /projects/b1052/Wells_b1042/yuboer/MAGs_faa/
#4）do hmmscan; /home/ywx1845/software/dbCAN_db/ is the absolute path of the formated dbCAN database
find . -name "*.faa" | while read line ; do hmmscan --domtblout ${line}.out.dm /home/ywx1845/software/dbCAN_db/dbCAN-HMMdb-V8.txt $line > ${line}.out; done
#5) parse the hmmscan results
find . -name "*.out.dm"|while read line ; do sh /home/ywx1845/software/dbCAN/hmmscan-parser.sh $line > ${line}.ps; done
#6) list the names of genomes with CAZy modules identified
find . -name "*.out.dm.ps" ! -size 0 > filteredhmmout.list
#7) make a directory and move the parsed hmmscan results to this directory, in my case, the directory name is 'dbCAN_annotation_results'
mkdir ../dbCAN_annotation_results
mv *.out.dm.ps > ../dbCAN_annotation_results
#8) delete the empty files in the folder of "dbCAN_annotation_results"
#9) the non-empty files in the directory "dbCAN_annotation_results" are input files to the following R scripts
