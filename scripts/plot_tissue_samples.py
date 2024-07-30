import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('resources/tissue_samples.txt',sep='\t',header=None)
df.columns = ['tissue','n']

df.sort_values('n',inplace=True,ascending=False)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.barh(range(df.shape[0]),df.n)
ax.set_xticks(np.arange(0,1000,200))
ax.set_yticks(range(df.shape[0]))
ax.set_yticklabels(df.tissue)
ax.set_ylabel('tissue')
ax.set_xlabel('nr. of samples')
ax.set_ylim([-.5,df.shape[0]-.5])

#N_exp_per_antigen.plot(kind='bar',ax=ax,width=1)
fig.set_size_inches([8,10])
plt.tight_layout()
fig.savefig(f'results/fig/GTEx_samples_per_tissue.pdf')