#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import seaborn as sns

colnames=("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
df1 = pd.read_csv ("chr_SN_1_ds_peaks.narrowPeak", sep = "\t", names=colnames)
df1['length']=df1['chromEnd']-df1['chromStart']
df1['replicate'] = "SN_1"
df1_row = df1.shape[0]
print ("SN_1 has " + str(df1_row) + " peaks")

df2 = pd.read_csv ("chr_SN_2_ds_peaks.narrowPeak", sep = "\t", names=colnames)
df2['length']=df2['chromEnd']-df2['chromStart']
df2['replicate']="SN_2"
df2_row = df2.shape[0]
print ("SN_2 has " + str(df2_row) + " peaks")

df3 = pd.read_csv ("chr_SN_3_ds_peaks.narrowPeak", sep = "\t", names=colnames)
df3['length']=df3['chromEnd']-df3['chromStart']
df3['replicate'] = "SN_3"
df3_row = df3.shape[0]
print ("SN_3 has " + str(df3_row) + " peaks")

df=pd.concat ([df1,df2,df3])
means = df.groupby('replicate')['length'].mean()
medians = df.groupby('replicate')['length'].median()
print ("Mean peak length:")
print (means)
print ("Median peak length:")
print (medians)

plt.figure(figsize=(8,6))
violin = sns.violinplot(y="length", 
               x="replicate", 
               data=df);

for i in range(len(means)):
    violin.text(i,means[i]/2, str(round(means[i])),
                fontdict=dict(color = "black", fontsize = 20),
                horizontalalignment = 'right');
violin.set_title('peak length');
plt.savefig("chr_SN_peak_length.svg",format='svg')
