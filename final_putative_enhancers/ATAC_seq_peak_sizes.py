#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import seaborn as sns

colnames=("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
df1 = pd.read_csv ("chr_Smed_combined_peaks.narrowPeak", sep = "\t", names=colnames)
df1['length']=df1['chromEnd']-df1['chromStart']
df1['replicate'] = "SM"
df1_row = df1.shape[0]
print ("SM has " + str(df1_row) + " peaks")

df2 = pd.read_csv ("chr_Spol_combined_peaks.narrowPeak", sep = "\t", names=colnames)
df2['length']=df2['chromEnd']-df2['chromStart']
df2['replicate']="SP"
df2_row = df2.shape[0]
print ("SP has " + str(df2_row) + " peaks")

df3 = pd.read_csv ("chr_Snov_combined_peaks.narrowPeak", sep = "\t", names=colnames)
df3['length']=df3['chromEnd']-df3['chromStart']
df3['replicate'] = "SN"
df3_row = df3.shape[0]
print ("SN has " + str(df3_row) + " peaks")

df4 = pd.read_csv ("chr_Slug_combined_peaks.narrowPeak", sep = "\t", names=colnames)
df4['length']=df4['chromEnd']-df4['chromStart']
df4['replicate'] = "SL"
df4_row = df4.shape[0]
print ("SL has " + str(df4_row) + " peaks")


df=pd.concat ([df1,df2,df3,df4])
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
plt.savefig("chr_Schmidteas_peak_length.svg",format='svg')
