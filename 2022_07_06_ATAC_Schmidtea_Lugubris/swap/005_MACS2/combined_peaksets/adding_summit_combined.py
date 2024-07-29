#!/usr/bin/env python
import matplotlib.pyplot
from upsetplot import UpSet
import pandas as pd
from matplotlib import pyplot
import seaborn as sns
import numpy as np
from statsmodels import robust
import sys
import argparse


if __name__ == '__main__':
   parser = argparse.ArgumentParser()
   parser.add_argument('--Summits_1')
   parser.add_argument('--Summits_2')
   parser.add_argument('--Summits_3')
   parser.add_argument('--Summits_4')
   parser.add_argument('--Summits_5')
   parser.add_argument('--Summits_6')
   parser.add_argument('--Peaks')
   parser.add_argument('--Output')
   args = parser.parse_args()

columns= ["Peak_chrom","Peak_chromStart","Peak_chromEnd","Peak_name","Peak_score","Peak_strand", "Summit_chrom","Summit_chromStart","Summit_chromEnd","Summit_name","Summit_score","Summit_distance","Summit_distance_to_peak_1"]
df_1= pd.read_csv(args.Summits_1, sep= "\t", names=columns, index_col=False)
df_1= df_1.set_index('Peak_name')

columns2= ["Peak_chrom","Peak_chromStart","Peak_chromEnd","Peak_name","Peak_score","Peak_strand", "Summit_chrom_2","Summit_chromStart_2","Summit_chromEnd_2","Summit_name_2","Summit_score_2","Summit_distance_2","Summit_distance_to_peak_2"]
df_2= pd.read_csv(args.Summits_2, sep= "\t", names=columns2, index_col=False)
df_2= df_2.set_index('Peak_name')

columns3= ["Peak_chrom","Peak_chromStart","Peak_chromEnd","Peak_name","Peak_score","Peak_strand", "Summit_chrom_3","Summit_chromStart_3","Summit_chromEnd_3","Summit_name_3","Summit_score_3","Summit_distance_3","Summit_distance_to_peak_3"]
df_3= pd.read_csv(args.Summits_3, sep= "\t", names=columns3, index_col=False)
df_3= df_3.set_index('Peak_name')

columns4= ["Peak_chrom","Peak_chromStart","Peak_chromEnd","Peak_name","Peak_score","Peak_strand", "Summit_chrom_4","Summit_chromStart_4","Summit_chromEnd_4","Summit_name_4","Summit_score_4","Summit_distance_4","Summit_distance_to_peak_4"]
df_4= pd.read_csv(args.Summits_4, sep= "\t", names=columns4, index_col=False)
df_4= df_4.set_index('Peak_name')

columns5= ["Peak_chrom","Peak_chromStart","Peak_chromEnd","Peak_name","Peak_score","Peak_strand", "Summit_chrom_5","Summit_chromStart_5","Summit_chromEnd_5","Summit_name_5","Summit_score_5","Summit_distance_5","Summit_distance_to_peak_5"]
df_5= pd.read_csv(args.Summits_5, sep= "\t", names=columns5, index_col=False)
df_5= df_5.set_index('Peak_name')

columns6= ["Peak_chrom","Peak_chromStart","Peak_chromEnd","Peak_name","Peak_score","Peak_strand", "Summit_chrom_6","Summit_chromStart_6","Summit_chromEnd_6","Summit_name_6","Summit_score_6","Summit_distance_6","Summit_distance_to_peak_6"]
df_6= pd.read_csv(args.Summits_6, sep= "\t", names=columns6, index_col=False)
df_6= df_6.set_index('Peak_name')

df_1_2 = pd.concat([df_1, df_2], axis=1)
df_1_2_3 = pd.concat([df_1_2, df_3], axis=1)
df_1_2_3_4 = pd.concat([df_1_2_3, df_4], axis=1)
df_1_2_3_4_5 = pd.concat([df_1_2_3_4, df_5], axis=1)
df_1_2_3_4_5_6 = pd.concat([df_1_2_3_4_5, df_6], axis=1)
df_1_2_3_4_5_6 = df_1_2_3_4_5_6.replace(r'^\s*$', np.nan, regex=True)

df_1_2_3_4_5_6['avg'] = df_1_2_3_4_5_6[['Summit_distance_to_peak_1','Summit_distance_to_peak_2', 'Summit_distance_to_peak_3','Summit_distance_to_peak_4','Summit_distance_to_peak_5', 'Summit_distance_to_peak_6']].mean(axis=1, skipna=True)
df_1_2_3_4_5_6=df_1_2_3_4_5_6.round()
summit=df_1_2_3_4_5_6[['avg']]


columns_final= ["chr","start","end","name","score","strand","signalValue","p-value","q-value"]
df= pd.read_csv(args.Peaks, sep= "\t", names=columns_final, index_col=False)
df= df.set_index('name')


df_final = pd.concat([df, summit], axis=1)
df_final['avg'] = df_final['avg'].fillna(0).astype(int)
df_final=df_final.reset_index()
df_final=df_final[["chr","start","end","index","score","strand","signalValue","p-value","q-value","avg"]]
df_final.to_csv(args.Output, sep="\t", index=False, header=False)
