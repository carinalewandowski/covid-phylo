# ---------------------------------------------------------------------------
# File name: ny_clean.py
# Authors: Arjun Sai Krishnan, Carina Lewandowski
# Description: series of functions for preprocessing hCoV sequences.
# ---------------------------------------------------------------------------

# import modules from biopython
from Bio import SeqIO, Seq
# import matplotlib
import matplotlib.pyplot as plt
# from datetime import date 
import numpy as np
import sys
from ete3 import Tree, NodeStyle, TreeStyle
from Bio.Alphabet import IUPAC

BAD_SEQ_IDS = ['hCoV-19/USA/NY-NYCPHL-000574/2020|EPI_ISL_632033|2020-04-01',
'hCoV-19/USA/NY-NYCPHL-001080/2020|EPI_ISL_633081|2020-10-17', 
'hCoV-19/USA/NY-NYCPHL-000973/2020|EPI_ISL_633003|2020-09-23']


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
    max_gaps_beg = 0
    max_gaps_end = 0
    for record in records:
        curr_seq = str(record.seq)
        # find length of consecutive gaps at the beginning of the seq
        i = 0
        while(curr_seq[i] == '-'):
            i += 1
        # find length of consecutive gaps at the end of the seq
        j = -1  
        while(curr_seq[j] == '-'):
            j -= 1
        # updated max lengths
        if (i > max_gaps_beg): max_gaps_beg = i
        if (-1*j > max_gaps_end): max_gaps_end = -1*j

    updated_records = []
    # iterate over each record and trim front and back of sequences
    for record in records:
        new_record = record
        new_seq = str(record.seq)[(max_gaps_beg - 1):(-1*max_gaps_end + 1)]
        new_seq = Seq.Seq(new_seq)
        new_record.seq = new_seq
        updated_records.append(new_record)

    # write trimmed records to a new text file in fasta format
    SeqIO.write(updated_records, output_fname, "fasta")
    print("Trimmed from front: ", max_gaps_beg)
    print("Trimmed from back: ", max_gaps_end)

# ---------------------------------------------------------------------------
# FIND SEGREGATING SITES
# ---------------------------------------------------------------------------
def find_seg_sites(records):
    print('Finding segregating sites...') 
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

    print('find_seg_sites done.') 
    print('Number of segregating sites:', len(seg_sites))

    return seg_sites

# ---------------------------------------------------------------------------
# COMPILE FILES (SEQS BY COUNTY --> SEQS BY STATE)
# ---------------------------------------------------------------------------
def compile_files(files, output_fname):
    compiled_records = []
    for f in files:
        records = list(SeqIO.parse(f, "fasta"))
        # parse county from file name
        county = f[9:].replace('.txt', '')
        for record in records:
            # add county to record description
            record.description = record.description + '|' + county
            # remove bad sequences (those that cause several gaps in MAFFT alignment)
            if record.id not in BAD_SEQ_IDS:
                compiled_records.append(record)
    # write compiled records to specified output file
    SeqIO.write(compiled_records, output_fname, "fasta")
        
# ---------------------------------------------------------------------------
# REMOVE SEQUENCES WITH GAPS (after final alignment)
# ---------------------------------------------------------------------------
def remove_seqs_with_gaps(input_fname, output_fname):
    print('Removing sequences with gaps...') 
    aligned_records = list(SeqIO.parse(input_fname, "fasta"))
    final_records = []
    for r in aligned_records:
        if '-' not in r.seq:
            final_records.append(r)
    SeqIO.write(final_records, output_fname, "fasta")
    print('remove_seqs_with_gaps done.')
    print('Number of remaining sequences: ', len(list(SeqIO.parse(output_fname, "fasta"))))

# ---------------------------------------------------------------------------
# FIND SNPs
# ---------------------------------------------------------------------------
def find_SNPs(input_fname): 
    print('Finding SNPs...') 
    ny_aligned = list(SeqIO.parse(input_fname, "fasta")) 
    seg_sites = find_seg_sites(ny_aligned) 
    print('List of SNPs:') 
    snps = [] 
    for ss in seg_sites: 
        to_print = True 
        for base in ss[1]: 
            if ss[1][base] < 13: 
                to_print = False # 13 = ~1% of total data 
        if to_print: 
            print('locus: ', ss[0], ss[1]) 
            snps.append(ss) 
    print('find_SNPs done.') 
    return(snps)

# ---------------------------------------------------------------------------
# PLOT ALLELE FREQUENCIES OVER TIME
# ---------------------------------------------------------------------------
def SNP_time_plot(input_fname, locus): 
    print('Computing allele frequencies...')
    monthdict = {}
    plotdict = {}
    ny_aligned = list(SeqIO.parse(input_fname, "fasta")) 
    for record in ny_aligned:
        name = record.id
        name = name.split(' ')
        name = name[0].split('|')
        time = name[-1]
        time = time.split('-')
        year = int(time[0])
        month = int(time[1])
        if len(time) < 3:
            day = 1
        else:
            day = int(time[2])
        # timestamp = datetime.date(year, month, day)
        allele = str(record.seq[locus])
        if month in monthdict:
            if allele in monthdict[month]:
                monthdict[month][allele] += 1
            else:
                monthdict[month][allele] = 1
        else:
            monthdict[month] = {}
            monthdict[month][allele] = 1

    keys = list(monthdict.keys())
    keys.sort()

    for month in keys:
        total = sum(monthdict[month].values())
        print(str(month) + " " + str(total))
        right_place = list(monthdict.keys()).index(month)
        for allele in monthdict[month]:
            if allele not in plotdict:
                plotdict[allele] = np.zeros(len(monthdict.keys()))
                plotdict[allele][right_place] = (monthdict[month][allele])/(total)
            else:
                plotdict[allele][right_place] = (monthdict[month][allele])/(total)
    
    for allele in plotdict:
        plt.plot(keys, plotdict[allele])
    
    plt.legend(list(plotdict.keys()))
    return

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

# ---------------------------------------------------------------------------
# CLEANING
# ---------------------------------------------------------------------------

# list of SARS-CoV-2 data by county in NY
files = ['NY_files/albany.txt', 'NY_files/manhattan.txt', 'NY_files/richmond.txt',
'NY_files/bronx.txt', 'NY_files/montgomery.txt', 'NY_files/rockland.txt', 'NY_files/brooklyn.txt', 
'NY_files/broome.txt', 'NY_files/new_rochelle.txt', 'NY_files/schenecktady.txt', 'NY_files/cayuga.txt',
'NY_files/new_york.txt', 'NY_files/staten_island.txt', 'NY_files/chenango.txt', 'NY_files/ononda.txt',
'NY_files/suffolk.txt', 'NY_files/clinton.txt',	'NY_files/orange.txt', 'NY_files/ulster.txt',
'NY_files/dutchess.txt', 'NY_files/passiac.txt', 'NY_files/washington.txt', 'NY_files/franklin.txt',
'NY_files/putnam.txt', 'NY_files/westchester.txt', 'NY_files/jefferson.txt', 'NY_files/queens.txt']

# STEP 1: compile files from counties into one file for the state
#compile_files(files, 'new_york_compiled.txt')

# STEP 2: remove extraneous sequences
#ny_records = list(SeqIO.parse('new_york_compiled.txt', "fasta"))
#remove_extr_seqs(ny_records, 'new_york_compiled_cleaned_NEW.txt')

# STEP 3: Check length of cleaned records
#ny_records_cleaned = list(SeqIO.parse('new_york_compiled_cleaned_NEW.txt', "fasta"))
#print(len(ny_records_cleaned))

# STEP 4: After uploading to MAFFT, parse the results and trim gaps on the ends of each sequence
#ny_MAFFT_records = list(SeqIO.parse("ny_MAFFT_results4.txt", "fasta"))
#remove_gaps(ny_MAFFT_records, 'ny_aligned4.txt')

# STEP 5: Find segregating sites
#ny_aligned = list(SeqIO.parse('ny_aligned4.txt', "fasta"))
#seg_sites = find_seg_sites(ny_aligned)

# STEP 6: Find sites of interest (sites with no gaps)
'''
sites_of_interest = []
for ss in seg_sites:
    site = ss[0]
    site_dict = ss[1]
    if '-' not in site_dict:
        sites_of_interest.append(site)
print(len(sites_of_interest))
'''
# STEP 7: Remove remaining sequences with any gaps. 
#remove_seqs_with_gaps('ny_aligned4.txt', 'final_ny_aligned.txt')

# STEP 8: Return to step 5. Repeat 5 and 6 and make sure that all remaining 
# segregating sites are sites of interest.

# ---------------------------------------------------------------------------
# BUILDING SITE FREQUENCY SPECTRUM
# ---------------------------------------------------------------------------

# STEP 1: Find SNPs
# snps = find_SNPs('final_ny_aligned.txt')
# # dictionary ??
# sfs_dict = {}
# freqs = []
# for snp in snps: 
#     allele_distr = snp[1].values() 
#     freq_derived_allele = min(allele_distr) 
#     if freq_derived_allele not in sfs_dict: 
#         sfs_dict[freq_derived_allele] = 1 
#     else: 
#         sfs_dict[freq_derived_allele] += 1 
#     freqs.append(freq_derived_allele)
# print(sfs_dict)
# #plt.bar(sfs_dict.keys(), sfs_dict.values())
# plt.hist(freqs, density = True)
# plt.title('Site Frequency Spectrum')
# plt.xlabel('Derived Allele Frequency')
# plt.ylabel('Proportion of SNPs')
# plt.show()

# SNP_time_plot('final_ny_aligned.txt', int(sys.argv[1]))
# plt.show()

# check_ends = list(SeqIO.parse("aligned_ref.txt", "fasta")) 
# remove_gaps(check_ends, "dump.txt")

# ---------------------------------------------------------------------------
# ANNOTATE PHYLOGENY
# ---------------------------------------------------------------------------
SNP_time_plot('final_ny_aligned.txt', 610)
plt.title("Frequency Plot for Locus 610+449")
plt.show()

# seqdict = {}
# ny_aligned = list(SeqIO.parse('final_ny_aligned.txt', "fasta")) 
# for record in ny_aligned:
#     new = str(record.description).split(' ')
#     # print(record.description)
#     new = new[0] + '_' + new[1]
#     new = new.replace('|', '_')
#     seqdict[new]=str(record.seq)

# t = Tree("final_ny_aligned.nh")

# for n in t.traverse():
#     if n.name != '':
#         names = n.name.split('_', 1)
#         new_name = names[1]
#         if seqdict[new_name][28432] == 'A':
#             nstyle = NodeStyle()
#             nstyle["fgcolor"] = "red"
#             nstyle["size"] = 15
#             n.set_style(nstyle)

# # for n in t.traverse():
# #     print(n.name)
# t.show()

# ---------------------------------------------------------------------------
# PHYLOGENY: PREP FOR BEAST
# ---------------------------------------------------------------------------
# Print file to NEXUS file
# ny_aligned = list(SeqIO.parse('final_ny_aligned.txt', "fasta")) 
# for record in ny_aligned:
#     desc = record.description.split(" ")
#     record.id = desc[1]
#     record.description = desc[1]
# SeqIO.write(ny_aligned, 'final_ny_aligned_name_fixed.txt', "fasta")
# count = SeqIO.convert("final_ny_aligned_name_fixed.txt", "fasta", "final_ny_aligned.nex", "nexus", alphabet=IUPAC.ambiguous_dna)
# print("Converted %i records" % count)


# ---------------------------------------------------------------------------
# Misc. Code for identifying corrupted sequences
# ---------------------------------------------------------------------------

# gaps: 21172:21226

#aligned_records = list(SeqIO.parse("ny_aligned3.txt", "fasta"))
#print(aligned_records[0].seq)

# identify sequence causing gaps at specified loci
'''for record in aligned_records:
    if record.seq[21172:21226] != '------------------------------------------------------':
        print(record.id)'''

# find remaining frequencies of gap counts to determine what else should be removed
'''gap_count_freqs = {}
records_with_gaps = []
for record in aligned_records:
    seq = record.seq
    count = 0
    for i in range(len(seq)):
        if seq[i] == '-':
            count += 1
    if count not in gap_count_freqs:
        gap_count_freqs[count] = 1
    else: gap_count_freqs[count] += 1
    if count != 3:
        records_with_gaps.append(record)
print(gap_count_freqs)'''

'''for j,r in enumerate(records_with_gaps):
    print(r.id)
    seq = record.seq
    max_gap_ct = 0
    curr_gap_ct = 0
    i = 0
    idx = -1
    while i < len(seq):
        if j == 0 and curr_gap_ct != 0:
            print(curr_gap_ct)
        if seq[i] == '-':
            curr_gap_ct += 1
        if ((i == len(seq) - 1) or (seq[i+1] != '-')):
            if curr_gap_ct > max_gap_ct: 
                idx = i - 2
                max_gap_ct = curr_gap_ct
            curr_gap_ct = 0
        i += 1
    print(idx)
    print('max contiguous gaps: ', max_gap_ct)
    print('\n')'''

# find first location of gaps
'''seq = aligned_records[0].seq
for i in range(len(seq)):
    if seq[i] == '-':
        print(i)'''

# identify sequence(s) causing gaps at specified locus
'''for record in aligned_records:
    seq = record.seq
    if seq[25309] != '-':
        print(record.id)'''