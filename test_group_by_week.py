# ---------------------------------------------------------------------------
# File name: test_group_by_week.py
# Authors: Arjun Sai Krishnan, Carina Lewandowski
# Description: file for testing/visualizing SNP timeplots by week
# ---------------------------------------------------------------------------

from ny_clean import SNP_time_plot_by_week
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

# ---------------------------------------------------------------------------
# SNP TIMEPLOTS (BY WEEK)
# Author: Both
# ---------------------------------------------------------------------------

# define constants
# constant for the offset of our data from the reference sequence
REF_OFFSET = 449
# x-axis labels for plots 
week_labels = ['Mar 2', 'Mar 9', 'Mar 16', 'Mar 23', 'Mar 30', 'Apr 6', 'Apr 13', 
                'Apr 20', 'Apr 27', 'May 4', 'May 11']


# ---------------------------------------------------------------------------
# SITE 28432 (28881)
# ---------------------------------------------------------------------------

locus = 28432

keys, df = SNP_time_plot_by_week("final_ny_aligned_name_fixed.txt", locus)
# kw = lambda x: x.isocalendar()[1]

#g = pd.to_datetime(df['date']).dt.weekofyear

#df = df.groupby(g.rename('epoch_week')).sum()
df.loc['Total'] = df.sum()
df = df.reset_index()

df['A_freq'] = df['A']/df['total']
df['G_freq'] = df['G']/df['total']
print(df)

plt.plot(week_labels, list(df['A_freq'])[0:11])
plt.plot(week_labels, list(df['G_freq'])[0:11])
alleles = ['A', 'G']
plt.legend(alleles)
title = 'Locus ' + str(locus + REF_OFFSET)
plt.title(title)
plt.xlabel('Week')
plt.ylabel('Frequency')
plt.show()

# ---------------------------------------------------------------------------
# SITE 28433 (28882)
# ---------------------------------------------------------------------------

locus = 28433
keys, df = SNP_time_plot_by_week("final_ny_aligned_name_fixed.txt", locus)

df.loc['Total'] = df.sum()
df = df.reset_index()

df['A_freq'] = df['A']/df['total']
df['G_freq'] = df['G']/df['total']
print(df)

plt.plot(week_labels, list(df['A_freq'])[0:11])
plt.plot(week_labels, list(df['G_freq'])[0:11])
alleles = ['A', 'G']
plt.legend(alleles)
title = 'Locus ' + str(locus + REF_OFFSET)
plt.title(title)
plt.xlabel('Week')
plt.ylabel('Frequency')
plt.show()

# ---------------------------------------------------------------------------
# SITE 28434 (28883)
# ---------------------------------------------------------------------------

locus = 28434

keys, df = SNP_time_plot_by_week("final_ny_aligned_name_fixed.txt", locus)

df.loc['Total'] = df.sum()
df = df.reset_index()

df['C_freq'] = df['C']/df['total']
df['G_freq'] = df['G']/df['total']
print(df)

plt.plot(week_labels, list(df['C_freq'])[0:11])
plt.plot(week_labels, list(df['G_freq'])[0:11])
alleles = ['C', 'G']
plt.legend(alleles)
title = 'Locus ' + str(locus + REF_OFFSET)
plt.title(title)
plt.xlabel('Week')
plt.ylabel('Frequency')
plt.show()

# ---------------------------------------------------------------------------
# SITE 18998 (18549)
# ---------------------------------------------------------------------------

locus = 18549

keys, df = SNP_time_plot_by_week("final_ny_aligned_name_fixed.txt", locus)

df.loc['Total'] = df.sum()
df = df.reset_index()

df['T_freq'] = df['T']/df['total']
df['C_freq'] = df['C']/df['total']
print(df)

plt.plot(week_labels, list(df['T_freq'])[0:11])
plt.plot(week_labels, list(df['C_freq'])[0:11])
alleles = ['T', 'C']
plt.legend(alleles)
title = 'Locus ' + str(locus + REF_OFFSET)
plt.title(title)
plt.xlabel('Week')
plt.ylabel('Frequency')
plt.show()

# ---------------------------------------------------------------------------
# SITE 29091 (29540)
# ---------------------------------------------------------------------------
locus = 29091
keys, df = SNP_time_plot_by_week("final_ny_aligned_name_fixed.txt", locus)

df.loc['Total'] = df.sum()
df = df.reset_index()

df['A_freq'] = df['A']/df['total']
df['G_freq'] = df['G']/df['total']
print(df)

plt.plot(week_labels, list(df['A_freq'])[0:11])
plt.plot(week_labels, list(df['G_freq'])[0:11])
alleles = ['A', 'G']
plt.legend(alleles)
title = 'Locus ' + str(locus + REF_OFFSET)
plt.title(title)
plt.xlabel('Week')
plt.ylabel('Frequency')
plt.show()
