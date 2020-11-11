from Bio import SeqIO

records = list(SeqIO.parse("gisaid_hcov-19_India_2020_11_11_15.txt", "fasta"))
# records = list(SeqIO.parse("test.txt", "fasta"))
sequences = []
is_bad = False
count = 0
n_count = 0

for record in records:
    curr_seq = str(record.seq)
    curr_seq.replace(' ', '')
    curr_seq.replace('\n', '')
    if 'N' in curr_seq:
        n_count=n_count+1
    if len(curr_seq) < 29000:
        is_bad = True
    for char in curr_seq:
        if char != 'A' and char != 'C' and char != 'T' and char !='G':
            print(char)
            count = count + 1
            print(count)
            is_bad = True
            break
    if not is_bad: 
        sequences.append(record)
    is_bad = False
    
SeqIO.write(sequences, "cleaned.txt", "fasta")
SeqIO.write(records[0], "test.txt", "fasta")
print("N_count = " + str(n_count))