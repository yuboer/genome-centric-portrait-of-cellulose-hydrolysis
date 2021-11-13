#! /usr/bin/python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import getopt
# starts at the second element of argv since the first one is the script name
file_name1=sys.argv[1]
#file_name1=raw_input("Enter the full name of genebank file: ")

gb = SeqIO.parse(file_name1, "genbank")
#output_handle0=open(file_name1+".fna","w")
output_handle1=open(file_name1+".ffn","w")
output_handle2=open(file_name1+".faa","w")
#output_handle3=open(file_name1+".traslated.faa","w")
output_handle4=open(file_name1+".info.tab","w")

for record in gb:
    print "Processing GenBank record %s" % record.id
    #output_handle0.write(">%s %s\n%s\n" % (record.id, record.description, record.seq))
    for feature in record.features :
        if(feature.type == "source"):
            taxaid = feature.qualifiers['db_xref']
            if len(taxaid) > 1 :
                output_handle4.write("%s\t%s\t%s\t\n" % (record.id, record.description, taxaid[1]))
            else:
                output_handle4.write("%s\t%s\t%s\t\n" % (record.id, record.description, taxaid[0]))
    for feature in record.features :
        if(feature.type == "CDS"):
            try:
                ID = feature.qualifiers['db_xref'][0]
                desc = feature.qualifiers['protein_id'][0]
                product = feature.qualifiers['product'][0]
                locus = feature.qualifiers['locus_tag'][0]
                nt_seq = feature.extract(record.seq)
                try:
                        assert len(feature.qualifiers['translation'])==1
                        aa_seq = feature.qualifiers['translation'][0]
                except AssertionError:
                        print ID,'no amni acids found!'
                        aa_seq=''
                        #aa_seq_t=nt_seq.translate()    
                output_handle1.write(">%s %s %s %s\n%s\n" % (ID, desc, product, locus, nt_seq))
                output_handle2.write(">%s %s %s %s\n%s\n" % (ID, desc, product, locus, aa_seq))
            #output_handle3.write(">%s %s %s %s\n%s\n" % (ID, desc, product, locus, aa_seq_t))
            except:
                continue

print 'Retrieving whole genome sequences!'
#SeqIO.convert(file_name1, "genbank", file_name1+".genome.fasta", "fasta")
#output_handle0.close()
output_handle1.close()
output_handle2.close()
#output_handle3.close()
output_handle4.close()
print 'Done!'
