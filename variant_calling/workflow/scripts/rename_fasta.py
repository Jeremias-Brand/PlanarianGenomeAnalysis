from Bio import SeqIO
import sys

count = 1

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]
ucsc_id = sys.argv[3]
name = sys.argv[4]

output = open(output_fasta, "w")

with open(input_fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        header = ">"+ucsc_id+"_"+str(count)+" "+name+" "+str(len(record.seq))
        output.write(header+"\n")
        output.write(str(record.seq)+"\n")
        count += 1
output.close()
