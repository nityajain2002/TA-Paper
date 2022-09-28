from Bio import SeqIO
import csv
from Bio import AlignIO

contig_takeout = ['N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N240_CP010176.1', 'N235_CP010150.1', 'N235_CP010150.1', 'N235_CP010150.1', 'N235_CP010150.1', 'N235_CP010150.1', 'N235_CP010150.1', 'N235_CP010150.1', 'N235_CP010150.1', 'N235_CP010150.1', 'N235_CP010150.1', 'N235_CP010150.1', 'N235_CP010150.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N250_CP010145.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N617_CP031546.1', 'N418_CP024821.1', 'N256_CP010167.1', 'N745_LR134209.1', 'N418_CP024821.1', 'N256_CP010167.1', 'N418_CP024821.1', 'N256_CP010167.1', 'N256_CP010167.1', 'N418_CP024821.1', 'N418_CP024821.1', 'N418_CP024821.1', 'N418_CP024821.1', 'N745_LR134209.1', 'N256_CP010167.1', 'N418_CP024821.1', 'N745_LR134209.1', 'N256_CP010167.1', 'N256_CP010167.1', 'N745_LR134209.1', 'N745_LR134209.1', 'N256_CP010167.1', 'N745_LR134209.1', 'N256_CP010167.1', 'N418_CP024821.1', 'N745_LR134209.1', 'N745_LR134209.1', 'N256_CP010167.1', 'N256_CP010167.1', 'N745_LR134209.1', 'N745_LR134209.1', 'N745_LR134209.1', 'N745_LR134209.1', 'N418_CP024821.1', 'N418_CP024821.1', 'N418_CP024821.1', 'N745_LR134209.1', 'N418_CP024821.1', 'N256_CP010167.1', 'N418_CP024821.1', 'N418_CP024821.1', 'N418_CP024821.1', 'N418_CP024821.1', 'N418_CP024821.1', 'N256_CP010167.1', 'N418_CP024821.1', 'N418_CP024821.1', 'N256_CP010167.1', 'N745_LR134209.1', 'N328_CP018323.1', 'N556_CP027462.1', 'N728_LR130552.1', 'N556_CP027462.1', 'N556_CP027462.1', 'N556_CP027462.1', 'N328_CP018323.1', 'N728_LR130552.1', 'N328_CP018323.1', 'N556_CP027462.1', 'N728_LR130552.1', 'N556_CP027462.1', 'N328_CP018323.1', 'N556_CP027462.1', 'N328_CP018323.1', 'N556_CP027462.1', 'N328_CP018323.1', 'N556_CP027462.1', 'N328_CP018323.1', 'N728_LR130552.1', 'N328_CP018323.1', 'N728_LR130552.1', 'N556_CP027462.1', 'N728_LR130552.1', 'N328_CP018323.1', 'N556_CP027462.1', 'N556_CP027462.1', 'N556_CP027462.1', 'N556_CP027462.1', 'N556_CP027462.1', 'N328_CP018323.1', 'N728_LR130552.1', 'N728_LR130552.1', 'N556_CP027462.1', 'N728_LR130552.1', 'N728_LR130552.1', 'N728_LR130552.1', 'N328_CP018323.1', 'N328_CP018323.1', 'N728_LR130552.1']

#contigs = []
#check_seqs = []
# loop through the csv file of hosts w/ no ta systems
with open("TAs_MLSTS_Hosts.csv", 'r') as input, open('specific-hosts.fasta','w+') as out:
    writer = csv.writer(out)
    reader = csv.reader(input)
    # loop through the fasta file to look for matches
    #i = 0
    seqs = []
    for row in reader:
        #print("row", row)
        if row[2] in contig_takeout:
            #row[17] is the file
            seqs.append('>' + row[2] + "\n" + row[17])
    #print(seqs)
    seq = '\n'.join(seqs)
    out.write(seq)
            #print("row", row[2])
            #print(row[17])
            #print(len(seq_record.seq))
            #writer.writerow(row)
            #i += 1
            #contigs.append(row[2])
            #check_seqs.append(str(seq_record.seq))
            #break

# check if they're an exact match
# match the seqs to the ones in the fasta file and take those ones out and put into a new fasta file

#base = check_seqs[0]
#print(contigs)
'''for i in check_seqs:
    if i == base:
        print("equal")
    else:
        print("not equal")'''

# only one host matches???
input.close()
out.close()
#shot gun sequecing method, but shouldn't have that

# trying with alignIO
'''alignment = AlignIO.read(open("ncbi761_host.fasta"), "fasta")
print("Alignment length %i" % alignment.get_alignment_length())
for record in alignment:
    print(record.seq + " " + record.id)'''