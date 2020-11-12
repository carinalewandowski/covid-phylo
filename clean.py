from Bio import SeqIO, Seq

# records = list(SeqIO.parse("gisaid_hcov-19_India_2020_11_11_15.txt", "fasta"))
# records = list(SeqIO.parse("test.txt", "fasta"))
sequences = []
is_bad = False
count = 0
n_count = 0

def remove_extr_seqs(records): 
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

def remove_gaps(records):
    for record in records:
        curr_seq = str(record.seq)
        i = 0
        max_gaps_beg = 0
        while(curr_seq[i] == '-'):
            i += 1
        j = -1
        max_gaps_end = 0
        while(curr_seq[j] == '-'):
            j -= 1
        if (i > max_gaps_beg): max_gaps_beg = i
        if (-1*j > max_gaps_end): max_gaps_end = j

    updated_records = []
    for record in records:
        new_record = record
        new_seq = str(record.seq)[(max_gaps_beg - 1):(max_gaps_end + 1)]
        new_seq = Seq.Seq(new_seq)
        new_record.seq = new_seq
        updated_records.append(new_record)

    SeqIO.write(updated_records, "aligned_cleaned.txt", "fasta")
    SeqIO.write(records[0], "test_aligned.txt", "fasta")
    print(i, j)

#records2 = list(SeqIO.parse("aligned.txt", "fasta"))
#remove_gaps(records2)

def find_seg_sites(records):
    site_dict = {}
    for record in records:
        curr_seq = str(record.seq)
        for i in range(len(curr_seq)):
            base = curr_seq[i]
            if (i not in site_dict):
                site_dict[i] = {}
            if (base not in site_dict[i]):
                site_dict[i][base] = 0    
            site_dict[i][base] += 1
    
    seg_sites = []
    for site in site_dict:
        if (len(site_dict[site]) > 1):
            seg_sites.append((site, site_dict[site]))

    return seg_sites

records3 = list(SeqIO.parse("aligned_cleaned.txt", "fasta"))
seg_sites = find_seg_sites(records3)
print(len(seg_sites))
