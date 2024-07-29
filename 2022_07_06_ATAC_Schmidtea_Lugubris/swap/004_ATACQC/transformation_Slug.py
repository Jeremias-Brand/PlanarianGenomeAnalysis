#!/usr/bin/env python
import pandas as pd
import numpy as np
import sys
colnames=['chrom','method', 'class', 'chromStart', 'chromEnd','score','strand','space', 'name'] 
df=pd.read_csv(sys.argv[1], sep='\t', comment='#', names=colnames, header=None, index_col=False)
df1= df.drop(df[(df['chrom'] == 'chr1') & (df['chromStart'] > 414900000)].index)
df2=df.loc[(df['chrom'] == 'chr1') & (df['chromStart'] > 414900000)]
df2['chromStart']=df2['chromStart'] - 414900000
df2['chromEnd']=df2['chromEnd'] - 414900000
df2=df2.replace({'chrom': 'chr1'}, {'chrom': 'chr5'}, regex=True)
df_new=pd.concat([df1,df2])
df_new
df_new.to_csv(sys.argv[2], index=False, header=False, sep='\t')
