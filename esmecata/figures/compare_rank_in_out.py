import csv

import pandas as pd
import seaborn
import matplotlib.pyplot as plt

RUN = 'bacteria4051'
IN_FILE = f'Esmecata_inputs/{RUN}.tsv'
OUT_FILE = f'Esmecata_output/{RUN}/proteome_tax_id.tsv'
HEATMAP = f'Esmecata_output/{RUN}/heatmap_rank_comparison.png'


out_asso_rank = {'superkingdom': 'k',
                 'kingdom': 'k',
                 'phylum': 'p',
                 'class': 'c',
                 'order': 'o',
                 'family': 'f',
                 'clade': 'f',
                 'genus': 'g',
                 'species': 's',
                 'no rank': 'na'}

rank_list = ['s', 'g', 'f', 'o', 'c', 'p', 'k', 'na']

d_in = dict()
with open(IN_FILE, 'r') as f:
    rows = list(csv.reader(f, delimiter='\t'))
    for row in rows[1:]:
        tax = row[0]
        ranks = list(reversed(row[1].split(';')))
        for i in range(len(rank_list) - 1):
            if ranks[i] != 'unknown':
                d_in[tax] = rank_list[i]
                break


d_out = dict()
df_out = pd.read_csv(OUT_FILE, sep='\t', index_col='observation_name')
for ASV in df_out.index:
    rank = out_asso_rank[df_out.loc[ASV, 'tax_rank']]
    d_out[ASV] = rank


matrix = pd.DataFrame(index=rank_list, columns=rank_list, dtype=int)
matrix = matrix.fillna(int(0))

for ASV in d_out.keys():
    in_rank = d_in[ASV]
    out_rank = d_out[ASV]
    matrix.loc[in_rank, out_rank] += 1

print(matrix)


plot = plt.plot()
seaborn.heatmap(matrix, annot=True, fmt='.0f', linewidth=0.5, cmap="crest")
plt.xlabel('Esmecata Rank')
plt.ylabel('Input Rank')
plt.savefig(HEATMAP)
