from ny_clean import SNP_time_plot_by_week
# import modules from biopython
from Bio import SeqIO, Seq
# import matplotlib
import matplotlib.pyplot as plt
import datetime
import numpy as np
import sys
from ete3 import Tree, NodeStyle, TreeStyle, random_color
#from Bio.Alphabet import IUPAC
import pandas as pd
from collapsetree import collapse

# keys, df = SNP_time_plot_by_week("final_ny_aligned_name_fixed.txt", 28432)
# kw = lambda x: x.isocalendar()[1]

g = pd.to_datetime(df['date']).dt.weekofyear

df = df.groupby(g.rename('epoch_week')).sum()
df.loc['Total'] = df.sum()
df = df.reset_index()

print(df)
df['A_freq'] = df['A']/df['total']
df['G_freq'] = df['G']/df['total']
print(df)

plt.plot(list(df['epoch_week'])[0:11], list(df['A_freq'])[0:11])
plt.plot(list(df['epoch_week'])[0:11], list(df['G_freq'])[0:11])
plt.show()