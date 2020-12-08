# ---------------------------------------------------------------------------
# File name: ny_clean.py
# Authors: Arjun Sai Krishnan, Carina Lewandowski
# Description: series of functions for preprocessing, analyzing, and
#              visualizing hCoV sequences from New York State.
# ---------------------------------------------------------------------------

# import modules from biopython
from Bio import SeqIO, Seq
# import matplotlib
import matplotlib.pyplot as plt
import datetime
import numpy as np
import sys
from ete3 import Tree, NodeStyle, TreeStyle, random_color
from Bio.Alphabet import IUPAC
import pandas as pd
from collapsetree import collapse

# IDs of sequences we identified as causing large sections of contiguous gaps
BAD_SEQ_IDS = ['hCoV-19/USA/NY-NYCPHL-000574/2020|EPI_ISL_632033|2020-04-01',
'hCoV-19/USA/NY-NYCPHL-001080/2020|EPI_ISL_633081|2020-10-17', 
'hCoV-19/USA/NY-NYCPHL-000973/2020|EPI_ISL_633003|2020-09-23']

# constant for the offset of our data from the reference sequence
REF_OFFSET = 449

# ---------------------------------------------------------------------------
# REMOVE EXTRANEOUS SEQUENCES 
# This chunk removes sequences that contain non-real bases or do not meet the 
# length threshold.
#
# Author: Both (Arjun typed)
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
#
# Author(s): Both (Arjun typed)
# ---------------------------------------------------------------------------

# link to MAFFT: https://mafft.cbrc.jp/alignment/software/

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
#
# Author(s): Both (Carina typed)
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
#
# Author(s): Both (Carina typed)
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
#
# Author(s): Both (Carina typed)
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
#
# Author: Carina
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
# BUILD SITE FREQUENCY SPECTRUM 
#
# Author: Carina
# ---------------------------------------------------------------------------
def build_sfs(input_fname):
    # STEP 1: Find SNPs
    snps = find_SNPs(input_fname)
    

    # STEP 2: Create a list of derived alleles by inspecting the SNP time plots
    derived_alleles = ['T', 'C', 'T', 'C', 'T', 'T', 'T', 'T', 'T', 'T',
    'C', 'T', 'G', 'T', 'C', 'T', 'T', 'A', 'C', 'A', 'G', 'T', 'C', 'A',
    'A', 'A', 'C', 'A']
    

    # STEP 3: For each SNP, record the derived allele frequency
    sfs_dict = {}
    freqs = []
    for i in range(len(snps)): 
        #allele_distr = snp[1].values() 
        #freq_derived_allele = min(allele_distr) 
        snp = snps[i]
        freq_derived_allele = snp[1][derived_alleles[i]]
        if freq_derived_allele not in sfs_dict: 
            sfs_dict[freq_derived_allele] = 1 
        else: 
            sfs_dict[freq_derived_allele] += 1 
        freqs.append(freq_derived_allele)
    print(sfs_dict)
    

    # STEP 4: Visualize as a histogram
    plt.hist(freqs)
    plt.title('Site Frequency Spectrum')
    plt.xlabel('Derived Allele Frequency')
    plt.ylabel('Number of SNPs')
    plt.show()


# ---------------------------------------------------------------------------
# FIND SNPs BY COUNTY
#
# Author: Carina
# ---------------------------------------------------------------------------
def find_SNPs_by_county(input_fname):

    derived_alleles = ['T', 'C', 'T', 'C', 'T', 'T', 'T', 'T', 'T', 'T',
        'C', 'T', 'G', 'T', 'C', 'T', 'T', 'A', 'C', 'A', 'G', 'T', 'C', 'A',
        'A', 'A', 'C', 'A']

    sites = [170, 610, 1468, 2588, 3664, 7629, 8333, 10402, 10634, 11467, 
            13959, 14356, 14463, 15811, 16798, 18428, 18549, 19556, 20306, 
            22954, 25114, 25695, 27695, 28372, 28432, 28433, 28434, 29091]

    ny_aligned = list(SeqIO.parse(input_fname, "fasta"))
    county_dict = {}
    num_records_by_county = {}
    for record in ny_aligned:
        county = record.description.split('|')[-1]
        if county not in county_dict:
            #county_dict[county] = np.zeros(len(sites))
            county_dict[county] = [0] * len(sites)
        if county not in num_records_by_county:   
            num_records_by_county[county] = 1
        else:
            num_records_by_county[county] += 1
        for i in range(len(sites)):
            curr_site = sites[i]
            if record.seq[curr_site] == derived_alleles[i]:
                county_dict[county][i] = county_dict[county][i] + 1
    return county_dict, num_records_by_county

# ---------------------------------------------------------------------------
# PLOT SNPs BY COUNTY
#
# Author: Carina
# ---------------------------------------------------------------------------
def plot_SNPs_by_county(input_fname):

    snp_counts, nums = find_SNPs_by_county(input_fname)

    sites = [170, 610, 1468, 2588, 3664, 7629, 8333, 10402, 10634, 11467, 
                13959, 14356, 14463, 15811, 16798, 18428, 18549, 19556, 20306, 
                22954, 25114, 25695, 27695, 28372, 28432, 28433, 28434, 29091]
    labels = []
    for site in sites:
        labels.append(str(site + REF_OFFSET))

    for county in snp_counts:   
        plt.bar(labels, np.array(snp_counts[county])/nums[county])
        plt.xlabel("SNP Locus")
        plt.xticks(rotation=45)
        plt.ylabel("Percent Individuals with Derived Allele")
        plt.title(county)
        plt.show()

plot_SNPs_by_county('final_ny_aligned.txt')

# ---------------------------------------------------------------------------
# HEATMAP OF SNPs BY COUNTY
#
# Author: Carina
# ---------------------------------------------------------------------------
def heatmap_SNPs_by_county(input_fname):
    
    snp_counts, nums = find_SNPs_by_county(input_fname)

    snp_counts_arr = []
    for county in snp_counts:
        arr = snp_counts[county]
        tot = nums[county]
        for i in range(len(arr)):
            arr[i] = arr[i]/tot
        snp_counts_arr.append(arr)

    sites = [170, 610, 1468, 2588, 3664, 7629, 8333, 10402, 10634, 11467, 
                13959, 14356, 14463, 15811, 16798, 18428, 18549, 19556, 20306, 
                22954, 25114, 25695, 27695, 28372, 28432, 28433, 28434, 29091]
    labels = []
    for site in sites:
        labels.append(str(site + REF_OFFSET))

    # Resource for heatmap: 
    # https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/image_annotated_heatmap.html
    counties = snp_counts.keys()
    fig, ax = plt.subplots()
    im = ax.imshow(snp_counts_arr)

    ax.set_xticks(np.arange(len(labels)))
    ax.set_yticks(np.arange(len(counties)))
    ax.set_xticklabels(labels)
    ax.set_yticklabels(counties)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

    plt.title('Heatmap of Derived Allele Frequencies by County')
    plt.xlabel('SNP Locus')
    plt.ylabel('County')
    plt.show()

# TO PLOT SNPs BY COUNTY
#plot_SNPs_by_county('final_ny_aligned.txt')
#heatmap_SNPs_by_county('final_ny_aligned.txt')


# ---------------------------------------------------------------------------
# PLOT ALLELE FREQUENCIES OVER TIME
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# BY MONTH
#
# Author: Arjun
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
    
    #for allele in plotdict:
        #plt.plot(keys, plotdict[allele])
    
    #plt.legend(list(plotdict.keys()))
    return keys, plotdict
# ---------------------------------------------------------------------------
# BY WEEKS
#
# Author: Arjun
# ---------------------------------------------------------------------------
def SNP_time_plot_by_week(input_fname, locus): 
    print('Computing allele frequencies...')
    daydict = {}
    plotdict = {}
    ny_aligned = list(SeqIO.parse(input_fname, "fasta")) 
    for record in ny_aligned:
        name = record.id
        name = name.split(' ')
        name = name[0].split('|')
        time = name[-2]
        time = time.split('-')
        year = int(time[0])
        month = int(time[1])
        if len(time) < 3:
            day = 1
        else:
            day = int(time[2])
        datestamp = datetime.datetime(year, month, day)
        allele = str(record.seq[locus])
        if datestamp in daydict:
            if allele in daydict[datestamp]:
                daydict[datestamp][allele] += 1
            else:
                daydict[datestamp][allele] = 1
        else:
            daydict[datestamp] = {}
            daydict[datestamp][allele] = 1

    keys = list(daydict.keys())
    keys.sort()

    plotdict["total"] = np.zeros(len(daydict.keys()))
    plotdict["date"] = []

    for day in keys:
        total = sum(daydict[day].values())
        right_place = list(daydict.keys()).index(day)
        plotdict["total"][right_place] = total
        plotdict["date"].append(day)
        for allele in daydict[day]:
            if allele not in plotdict:
                plotdict[allele] = np.zeros(len(daydict.keys()))
                plotdict[allele][right_place] = (daydict[day][allele])
            else:
                plotdict[allele][right_place] = (daydict[day][allele])
  
    df = pd.DataFrame.from_dict(plotdict)

    g = pd.to_datetime(df['date']).dt.weekofyear

    df = df.groupby(g.rename('epoch_week')).sum()
    df.loc['Total'] = df.sum()
    df = df.reset_index()

    return keys, df
# ---------------------------------------------------------------------------
# BUILD SNP TABLE
#
# Author: Carina
# ---------------------------------------------------------------------------

# Helpful resource: https://www.geeksforgeeks.org/create-a-pandas-dataframe-from-lists/

def build_snp_table(input_fname, output_fname, c_derived_lb, c_genes):
    snps = find_SNPs(input_fname)

    c_loci = []
    c_original_lb = []
    c_original_fr = []
    c_derived_lb =  c_derived_lb
    c_derived_fr = []
    c_genes = c_genes
    print('Building SNP table...')
    for i in range(len(snps)):
        snp = snps[i]
        c_loci.append(snp[0] + REF_OFFSET)
        for key in snp[1]:
            if key == c_derived_lb[i]:
                c_derived_fr.append(snp[1][key])
            else:
                c_original_lb.append(key)
                c_original_fr.append(snp[1][key])

    df = pd.DataFrame(list(zip(c_loci, c_original_lb, c_original_fr, c_derived_lb, 
                                c_derived_fr, c_genes)),
                        columns=['Locus', 'Original Allele', 'Frequency', 
                                    'Derived Allele', 'Frequency', 'Gene'])
    df.to_csv(output_fname)
    print('SNP data output to csv: ' + output_fname)


# ---------------------------------------------------------------------------
# TESTING
#
# Author: Both
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
#
# Author: Both (Carina typed)
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

# extra checks (Arjun)
'''
check_ends = list(SeqIO.parse("aligned_ref.txt", "fasta")) 
remove_gaps(check_ends, "dump.txt")
'''

# ---------------------------------------------------------------------------
# SNP ANALYSIS
#
# Author: Carina
# ---------------------------------------------------------------------------

# STEP 1: build site frequency spectrum
'''
build_sfs('final_ny_aligned.txt')
'''

# STEP 2: prepare variables
# list of SNP loci (Want 7 plots with 4 subplots each)
'''
loci1 = [[170, 610], [1468, 2588]] 
loci2 = [[3664, 7629], [8333, 10402]]
loci3 = [[10634, 11467], [13959, 14356]]
loci4 = [[14463, 15811], [16798, 18428]]
loci5 = [[18549, 19556], [20306, 22954]]
loci6 = [[25114, 25695], [27695, 28372]]
loci7 = [[28432, 28433], [28434, 29091]]
loci = [loci1, loci2, loci3, loci4, loci5, loci6, loci7]
rows = 2 
columns = 2
sets = 7
months = ['Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct']
'''
# STEP 3: Iterate over each set of loci and make 2x2 subplots
# Helpful Resource: https://matplotlib.org/3.1.0/gallery/subplots_axes_and_figures/subplots_demo.html
'''
for k in range(sets):
    curr_loci = loci[k]
    fig, axs = plt.subplots(rows, columns)
    for i in range(rows):
        for j in range(columns):
            keys, plotdict = SNP_time_plot('final_ny_aligned.txt', curr_loci[i][j])
            for allele in plotdict:
                axs[i,j].plot(months, plotdict[allele])
            axs[i,j].legend(list(plotdict.keys()))
            title = 'Locus ' + str(curr_loci[i][j] + REF_OFFSET)
            axs[i,j].set_title(title)
            axs[i,j].set_xticklabels(months, rotation=45)

    for ax in axs.flat:
        ax.set(xlabel='Time (months)', ylabel='Frequency')
       

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    plt.show()
'''

# STEP 4: Create a table of SNP data.
'''
c_derived_lb = ['T', 'C', 'T', 'C', 'T', 'T', 'T', 'T', 'T', 'T',
    'C', 'T', 'G', 'T', 'C', 'T', 'T', 'A', 'C', 'A', 'G', 'T', 'C', 'A',
    'A', 'A', 'C', 'A']

c_genes = ['ORF1ab', 'ORF1ab', 'ORF1ab', 'ORF1ab', 'ORF1ab', 'ORF1ab', 
    'ORF1ab', 'ORF1ab', 'ORF1ab', 'ORF1ab', 'ORF1ab', 'ORF1ab', 'ORF1ab', 
    'ORF1ab', 'ORF1ab', 'ORF1ab', 'ORF1ab', 'ORF1ab', 'ORF1ab', 'S', 'ORF3a',
    'ORF3a', 'ORF8', 'N', 'N', 'N', 'N', 'Intergenic']
build_snp_table('final_ny_aligned.txt', 'snp_data.csv', c_derived_lb, c_genes)
'''

# ---------------------------------------------------------------------------
# ANNOTATE PHYLOGENY
# 
# Author: Arjun
# ---------------------------------------------------------------------------
def phylo_annotate_v1():
    # SNP_time_plot('final_ny_aligned.txt', 610)
    # plt.title("Frequency Plot for Locus 610+449")
    # plt.show()

    seqdict = {}
    ny_aligned = list(SeqIO.parse('final_ny_aligned.txt', "fasta")) 
    for record in ny_aligned:
        new = str(record.description).split(' ')
        # print(record.description)
        new = new[0] + '_' + new[1]
        new = new.replace('|', '_')
        seqdict[new]=str(record.seq)

    t = Tree("final_ny_aligned_newick.nh")

    for n in t.traverse():
        if n.name != '':
            names = n.name.split('_', 1)
            new_name = names[1]
            if seqdict[new_name][28432] == 'A':
                nstyle = NodeStyle()
                nstyle["fgcolor"] = "red"
                nstyle["size"] = 15
                n.set_style(nstyle)

    # for n in t.traverse():
    #     print(n.name)

    t.show()

# WHEN THE NAMES ARE OK..
def double_or_triple():
    seqdict = {}
    ny_aligned = list(SeqIO.parse('final_ny_aligned_name_fixed.txt', "fasta")) 
    for record in ny_aligned:
        # new = str(record.description).split(' ')
        # # print(record.description)
        # new = new[0] + '_' + new[1]
        # new = new.replace('|', '_')
        seqdict[record.id]=str(record.seq)

    t = Tree("beast_tree_newick.nh")

    place_list = ['albany', 'manhattan', 'richmond',
    'bronx', 'montgomery', 'rockland', 'brooklyn', 
    'broome', 'new_rochelle', 'schenecktady', 'cayuga',
    'new_york', 'staten_island', 'chenango', 'ononda',
    'suffolk', 'clinton',	'orange', 'ulster',
    'dutchess', 'passiac', 'washington', 'franklin',
    'putnam', 'westchester', 'jefferson', 'queens']
    # color_list = ["red", "yellow", "blue", "green", "cyan", "magenta", "orange", "light blue", "grey", "purple", "brown", 
    #                 "pink", "light green"]
    color_list = random_color(num=len(place_list))
    color_list_mut = random_color(num=9)

    mut_list = ['AAC', 'AAG', 'AGC', 'AGG', 'GAG', 'GGC', 'GAC', 'GGG']

    for n in t.traverse():
        if n.name == '' or n.name == None:
            nstyle = NodeStyle()
            nstyle["fgcolor"] = "light grey"
            nstyle["size"] = 0
            n.set_style(nstyle)
        else:
            n.name = n.name.replace("\'", "")
            nstyle = NodeStyle()
            sequence = seqdict[n.name][28432:28435]
            color = color_list_mut[mut_list.index(sequence)]
            nstyle["hz_line_color"] = color
            nstyle["fgcolor"] = color
            n.set_style(nstyle)

    ## --------------------------------------------------------------------------
    ##  FOR TRIPLE MUTANT
    ## 
    ## Author: Arjun
    ## --------------------------------------------------------------------------
    for n in t.traverse():
        new_name = n.name.replace("\'", "")
        if n.name == None or '|' not in n.name:
            nstyle_blank = NodeStyle()
            nstyle_blank["fgcolor"] = color_list_mut[8]
            nstyle_blank["size"] = 0
            n.set_style(nstyle_blank)
            # print(n.name)
            # names = n.name.split('_', 1)
        elif seqdict[new_name][28432] == 'A' and seqdict[new_name][28433]== 'A' and seqdict[new_name][28434]== 'C':
            nstyle_AAC = NodeStyle()
            nstyle_AAC["fgcolor"] = "red"
            nstyle_AAC["hz_line_color"] = "red"
            # nstyle["fgcolor"] = color_list_mut[0]
            # nstyle["hz_line_color"] = color_list_mut[0]
            nstyle_AAC["size"] = 15
            n.set_style(nstyle_AAC)
            n.name = 'AAC'
        elif seqdict[new_name][28432] == 'A' and seqdict[new_name][28433]== 'A' and seqdict[new_name][28434]== 'G':
            nstyle_AAG = NodeStyle()
            # nstyle_AAG["fgcolor"] = "yellow"
            nstyle_AAG["fgcolor"] = color_list_mut[1]
            nstyle_AAG["hz_line_color"] = color_list_mut[1]
            nstyle_AAG["size"] = 15
            n.set_style(nstyle_AAG)
            n.name = 'AAG'
        elif seqdict[new_name][28432] == 'A' and seqdict[new_name][28433]== 'G' and seqdict[new_name][28434]== 'C':
            nstyle_AGC = NodeStyle()
            # nstyle_AGC["fgcolor"] = "purple"
            nstyle_AGC["fgcolor"] = color_list_mut[2]
            nstyle_AGC["hz_line_color"] = color_list_mut[2]
            nstyle_AGC["size"] = 15
            n.set_style(nstyle_AGC)
            n.name = 'AGC'
        elif seqdict[new_name][28432] == 'A' and seqdict[new_name][28433]== 'G' and seqdict[new_name][28434]== 'G':
            nstyle_AGG = NodeStyle()
            # nstyle_AGG["fgcolor"] = "green"
            nstyle_AGG["fgcolor"] = color_list_mut[3]
            nstyle_AGG["hz_line_color"] = color_list_mut[3]
            nstyle_AGG["size"] = 15
            n.set_style(nstyle_AGG)
            n.name = 'AGG'
        elif seqdict[new_name][28432] == 'G' and seqdict[new_name][28433]== 'A' and seqdict[new_name][28434]== 'G':
            nstyle_GAG = NodeStyle()
            # nstyle_GAG["fgcolor"] = "grey"
            nstyle_GAG["fgcolor"] = color_list_mut[4]
            nstyle_GAG["hz_line_color"] = color_list_mut[4]
            nstyle_GAG["size"] = 15
            n.set_style(nstyle_GAG)
            n.name = 'GAG'
        elif seqdict[new_name][28432] == 'G' and seqdict[new_name][28433]== 'G' and seqdict[new_name][28434]== 'C':
            nstyle_GGC = NodeStyle()
            # nstyle_GGC["fgcolor"] = "light blue"
            nstyle_GGC["fgcolor"] = color_list_mut[5]
            nstyle_GGC["hz_line_color"] = color_list_mut[5]
            nstyle_GGC["size"] = 15
            n.set_style(nstyle_GGC)
            n.name = 'GGC'
        elif seqdict[new_name][28432] == 'G' and seqdict[new_name][28433]== 'A' and seqdict[new_name][28434]== 'C':
            nstyle_GAC = NodeStyle()
            # nstyle_GAC["fgcolor"] = "orange"
            nstyle_GAC["fgcolor"] = color_list_mut[6]
            nstyle_GAC["hz_line_color"] = color_list_mut[6]
            nstyle_GAC["size"] = 15
            n.set_style(nstyle_GAC)
            n.name = 'GAC'
        elif seqdict[new_name][28432] == 'G' and seqdict[new_name][28433]== 'G' and seqdict[new_name][28434]== 'G':
            n.name = 'GGG'
            nstyle_GGG = NodeStyle()
            nstyle_GGG["fgcolor"] = color_list_mut[7]
            nstyle_GGG["hz_line_color"] = color_list_mut[7]
            nstyle_GGG["size"] = 15
            n.set_style(nstyle_GGG)
            # namesplit = new_name.split('|')
            # n.name = namesplit[-1]

    ## --------------------------------------------------------------------------
    ##  FOR DOUBLE MUTANT
    ## 
    ## Author: Arjun
    ## --------------------------------------------------------------------------
    
    for n in t.traverse():
        new_name = n.name.replace("\'", "")
        if n.name == None or '|' not in n.name:
            nstyle_blank = NodeStyle()
            nstyle_blank["fgcolor"] = color_list_mut[8]
            nstyle_blank["size"] = 0
            n.set_style(nstyle_blank)
        elif seqdict[new_name][18549] == 'C' and seqdict[new_name][29091] == 'G':
            n.name = 'CG'
            nstyle_CG = NodeStyle()
            nstyle_CG["fgcolor"] = color_list_mut[0]
            nstyle_CG["hz_line_color"] = color_list_mut[0]
            nstyle_CG["size"] = 15
            n.set_style(nstyle_CG)
        elif seqdict[new_name][18549] == 'C' and seqdict[new_name][29091] == 'A':
            n.name = 'CA'
            nstyle_CA = NodeStyle()
            nstyle_CA["fgcolor"] = color_list_mut[1]
            nstyle_CA["hz_line_color"] = color_list_mut[1]
            nstyle_CA["size"] = 15
            n.set_style(nstyle_CA)
        elif seqdict[new_name][18549] == 'T' and seqdict[new_name][29091] == 'G':
            n.name = 'TG'
            nstyle_TG = NodeStyle()
            nstyle_TG["fgcolor"] = color_list_mut[2]
            nstyle_TG["hz_line_color"] = color_list_mut[2]
            nstyle_TG["size"] = 15
            n.set_style(nstyle_TG)
        elif seqdict[new_name][18549] == 'T' and seqdict[new_name][29091] == 'A':
            n.name = 'TA'
            nstyle_TA = NodeStyle()
            nstyle_TA["fgcolor"] = color_list_mut[3]
            nstyle_TA["hz_line_color"] = color_list_mut[3]
            nstyle_TA["size"] = 15
            n.set_style(nstyle_TA)

    # for n in t.traverse():
    #     if n.name != '':
    #         sequence = seqdict[n.name][28432:28435]
    #         if n.nstyle["fgcolor"] != color_list_mut[mut_list.index(sequence)]:


    # for bad in badlist:
    #     new_bad = bad.replace("\'", "")
    #     if new_bad in seqdict:
    #         print(seqdict[new_bad][28432:28435])
    #     else:
    #         print('Name: ' + str(bad))

    make_long = TreeStyle()
    make_long.scale = 1500
    # t.show(tree_style=make_long)
    # print(len(badlist))
    # print(badlist)
    # for n in t.traverse():
    #     if n.name == '':
    #         nstyle = NodeStyle()
    #         nstyle["fgcolor"] = "light grey"
    #         nstyle["size"] = 0
    #         n.set_style(nstyle)
    #     else:
    #         n.name = n.name.replace("\'", "")
    #         nstyle = NodeStyle()
    #         nstyle["hz_line_color"] = color_list[place_list.index(n.name.split('|')[-1])]
    #         nstyle["fgcolor"] = color_list[place_list.index(n.name.split('|')[-1])]
    #         n.set_style(nstyle)
    # make_long = TreeStyle()
    # make_long.scale = 15000
    collapse(t)
    collapse(t)
    collapse(t)
    collapse(t)
    collapse(t)
    # collapse(t)
    # collapse(t)
    t.show(tree_style=make_long)
    # t.write(format=1, outfile="place_name_tree_ny.nh")
    
# ---------------------------------------------------------------------------
# PHYLOGENY: PREP FOR BEAST
# 
# Author: Arjun
# ---------------------------------------------------------------------------
def prep_for_beast():
    # print file to NEXUS file
    ny_aligned = list(SeqIO.parse('final_ny_aligned.txt', "fasta")) 
    for record in ny_aligned:
        desc = record.description.split(" ")
        record.id = desc[1]
        record.description = desc[1]
    SeqIO.write(ny_aligned, 'final_ny_aligned_name_fixed.txt', "fasta")
    count = SeqIO.convert("final_ny_aligned_name_fixed.txt", "fasta", "final_ny_aligned.nex", "nexus", alphabet=IUPAC.ambiguous_dna)
    print("Converted %i records" % count)


# ---------------------------------------------------------------------------
# Misc. Code for identifying corrupted sequences
# 
# Author: Both
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

# ---------------------------------------------------------------------------
# EXPERIMENTS WITH NEW LARGER ALIGNMENT
# 
# Author: Arjun
# ---------------------------------------------------------------------------
# ny_aligned_LARGE = list(SeqIO.parse('new_compiled_aligned.txt', "fasta")) 
# ny_aligned_original = list(SeqIO.parse('new_york_compiled.txt', "fasta")) 
# record = ny_aligned_LARGE[0]
# new_records = []

# for i, s in enumerate(str(record.seq)):
#     if s == '-':
#         print(i)

# NEW_BAD_SEQIDS = BAD_SEQ_IDS.copy()
        
# GAP_INDICES = [11074, 11075,  11076, 28253, 28254, 28255, 29897, 29898, 29899, 
#                 29900, 29901, 29902, 29903, 29904, 29905, 29906, 29907, 29908]
# print(len(ny_aligned_LARGE))
# for record in ny_aligned_LARGE[1:]:
#     for index in GAP_INDICES:
#         if record.seq[index] != '-':
#             if record.id not in NEW_BAD_SEQIDS:
#                 NEW_BAD_SEQIDS.append(str(record.id))
#                 print(record.id.split("|")[-1])

# for record in ny_aligned_LARGE[1:]:
#     if record.id not in NEW_BAD_SEQIDS:
#         record.id = record.description
#         new_records.append(record)
# print(len(NEW_BAD_SEQIDS))

# location_dict = {}
# time_dict = {}

# for record in new_records[1:]:
#     location = record.id.split('|')[-1]
#     time = record.id.split('|')[-2].split('-')[1]
#     if location in location_dict:
#         location_dict[location]+=1
#     else:
#         location_dict[location] = 1
#     if time in time_dict:
#         time_dict[time]+=1
#     else:
#         time_dict[time] = 1

# print(time_dict)
# print(location_dict)

# final_records = []
# for record in ny_aligned_original:
#     if record.id not in NEW_BAD_SEQIDS:
#         final_records.append(record)

# SeqIO.write(final_records, "new_compiled.txt", 'fasta')

# # READ IN NEW ALIGNMENT FROM MAFFT
# ny_aligned_LARGE_NEW = list(SeqIO.parse('new_compiled_aligned_added.txt', "fasta")) 
# print(find_seg_sites(ny_aligned_LARGE_NEW))
# SNP_time_plot('new_compiled_aligned_added.txt', 28881)
# plt.show()

# ny_aligned_for_R = list(SeqIO.parse('final_ny_aligned_name_fixed.txt', "fasta"))

# id_list = []
# genotype_list = []

# for record in ny_aligned_for_R:
#     id_list.append(record.id)
#     genotype_list.append(record.seq[28432:28435])


# df = pd.DataFrame(list(zip(id_list, genotype_list)),
#                         columns=["id", "genotype"])
# df.to_csv("sites_28432_28434.csv")