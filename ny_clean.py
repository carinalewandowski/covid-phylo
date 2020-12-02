# ---------------------------------------------------------------------------
# File name: ny_clean.py
# Authors: Arjun Sai Krishnan, Carina Lewandowski
# Description: series of functions for preprocessing hCoV sequences.
# ---------------------------------------------------------------------------

# import modules from biopython
from Bio import SeqIO, Seq

# ---------------------------------------------------------------------------
# REMOVE EXTRANEOUS SEQUENCES 
# This chunk removes sequences that contain non-real bases or do not meet the 
# length threshold.
# ---------------------------------------------------------------------------

def remove_extr_seqs(records, output_fname): 
    sequences = []
    is_bad = False
    count = 0
    threshold = 29000 # length threshold

    # iterate over each record & parse its sequence as a string
    for record in records:
        curr_seq = str(record.seq)
        # trim white space for the purpose of inspecting base characters
        curr_seq.replace(' ', '')
        curr_seq.replace('\n', '')
        # check that seq length meets threshold
        if len(curr_seq) < threshold:
            is_bad = True
        # visit each position to check for extraneous base; break and reomve
        # sequence once one is found.
        for char in curr_seq:
            if char != 'A' and char != 'C' and char != 'T' and char !='G':
                #print(char)
                #count = count + 1
                #print(count)
                is_bad = True
                break
        # append sequence to list if it is not extraneous
        if not is_bad: 
            sequences.append(record)
        # reset boolean check
        is_bad = False
        
    # write non-extraneous sequences to a new text file in fasta format
    SeqIO.write(sequences, output_fname, "fasta")

# ---------------------------------------------------------------------------
# REMOVE GAPS
# after alignment with MAFFT, gaps typically appear in the largest amounts at
# the beginning and ends of sequences. This chunk finds the right number of
# bases to trim from the beginning and ends of sequences.
# ---------------------------------------------------------------------------
def remove_gaps(records, output_fname):
    # iterate over each record & parse its sequence as a string
    for record in records:
        curr_seq = str(record.seq)
        # find length of consecutive gaps at the beginning of the seq
        i = 0
        max_gaps_beg = 0
        while(curr_seq[i] == '-'):
            i += 1
        # find length of consecutive gaps at the end of the seq
        j = -1
        max_gaps_end = 0
        while(curr_seq[j] == '-'):
            j -= 1
        # updated max lengths
        if (i > max_gaps_beg): max_gaps_beg = i
        if (-1*j > max_gaps_end): max_gaps_end = j

    updated_records = []
    # iterate over each record and trim front and back of sequences
    for record in records:
        new_record = record
        new_seq = str(record.seq)[(max_gaps_beg - 1):(max_gaps_end + 1)]
        new_seq = Seq.Seq(new_seq)
        new_record.seq = new_seq
        updated_records.append(new_record)

    # write trimmed records to a new text file in fasta format
    SeqIO.write(updated_records, output_fname, "fasta")
    print("Trimmed from front: ", i)
    print("Trimmed from back: ", j)

# ---------------------------------------------------------------------------
# FIND SEGREGATING SITES
# ---------------------------------------------------------------------------
def find_seg_sites(records):
    site_dict = {}
    # iterate over each record & parse its sequence as a string
    for record in records:
        curr_seq = str(record.seq)
        # iterate over each position in the current seq
        for i in range(len(curr_seq)):
            base = curr_seq[i]
            # create a dict for each seq position
            # (key: position, value: dict of bases)
            if (i not in site_dict):
                site_dict[i] = {}
            # create a dict for each base found at each position
            # (key: base, value: frequency)
            if (base not in site_dict[i]):
                site_dict[i][base] = 0    
            site_dict[i][base] += 1
    
    # if more than one base dict exists at a position, classify as a seg site
    seg_sites = []
    for site in site_dict:
        if (len(site_dict[site]) > 1):
            seg_sites.append((site, site_dict[site]))

    return seg_sites

def compile_files(files, output_fname):
    compiled_records = []
    for f in files:
        records = list(SeqIO.parse(f, "fasta"))
        county = f[9:].replace('.txt', '')
        for record in records:
            record.description = record.description + '|' + county
            compiled_records.append(record)
    SeqIO.write(compiled_records, output_fname, "fasta")
        


# ---------------------------------------------------------------------------
# TESTING
# ---------------------------------------------------------------------------

''' TESTING remove_extr_seqs '''
# records = list(SeqIO.parse("gisaid_hcov-19_India_2020_11_11_15.txt", "fasta"))
# records = list(SeqIO.parse("test.txt", "fasta"))

''' TESTING remove_gaps '''
# records2 = list(SeqIO.parse("aligned.txt", "fasta"))
# remove_gaps(records2)

''' TESTING find_seg_sites '''
# records3 = list(SeqIO.parse("aligned_cleaned.txt", "fasta"))
# seg_sites = find_seg_sites(records3)
# print("Number of segregating sites:", len(seg_sites))

files = ['NY_files/albany.txt', 'NY_files/manhattan.txt', 'NY_files/richmond.txt',
'NY_files/bronx.txt', 'NY_files/montgomery.txt', 'NY_files/rockland.txt', 'NY_files/brooklyn.txt', 
'NY_files/broome.txt', 'NY_files/new_rochelle.txt', 'NY_files/schenecktady.txt', 'NY_files/cayuga.txt',
'NY_files/new_york.txt', 'NY_files/staten_island.txt', 'NY_files/chenango.txt', 'NY_files/ononda.txt',
'NY_files/suffolk.txt', 'NY_files/clinton.txt',	'NY_files/orange.txt', 'NY_files/ulster.txt',
'NY_files/dutchess.txt', 'NY_files/passiac.txt', 'NY_files/washington.txt', 'NY_files/franklin.txt',
'NY_files/putnam.txt', 'NY_files/westchester.txt', 'NY_files/jefferson.txt', 'NY_files/queens.txt']

test_file = ['NY_files/cayuga.txt',]
compile_files(files, 'new_york_compiled.txt')

ny_records = list(SeqIO.parse('new_york_compiled.txt', "fasta"))
remove_extr_seqs(ny_records, 'new_york_compiled_cleaned.txt')

ny_records_cleaned = list(SeqIO.parse('new_york_compiled_cleaned.txt', "fasta"))
print(len(ny_records_cleaned))
