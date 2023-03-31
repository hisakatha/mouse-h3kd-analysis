#!/usr/bin/env python3
from locale import normalize
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from typing import Any

def plot_CPMs(data_dict: dict[str, Any]):
    pdf_size = (5,4)
    xlab = 'CPM'
    with PdfPages('cpm_hist.pdf') as pdf:
        for name, data in data_dict.items():
            # all
            plt.figure(figsize=pdf_size)
            data['value'].plot.hist(bins=1000, logy=True)
            plt.title(name)
            plt.xlabel(xlab)
            pdf.savefig()
            plt.close()
            # zoom
            plt.figure(figsize=pdf_size)
            data.loc[data['value'] < 1, 'value'].plot.hist(bins=10, logy=True)
            plt.title(name)
            plt.xlabel(xlab)
            pdf.savefig()
            plt.close()
            # zoom CDF
            plt.figure(figsize=pdf_size)
            data['value'].plot.hist(xlim=(0,1), bins=100000, cumulative=True, density=True, histtype='step')
            plt.title(name)
            plt.xlabel(xlab)
            pdf.savefig()
            plt.close()

bedgraph_dtypes = {'chr':str, 'start':int, 'end':int, 'value':float}
data_dirs = {
    'MII #1': '../NoRibs_MII201209/',
    'MII #2': '../NoRibs_MII201218/',
    '2C 18hpi #1': '../NoRibs_2c_18hpi201209/',
    '2C 18hpi #2': '../NoRibs_2c_18hpi201218/',
    '2C No Injection #1': '../NoRibs_2c_Noinjection201209/',
    '2C No Injection #2': '../NoRibs_2c_Noinjection201218/',
    '2C Control #1': '../NoRibs_2c_Control201209/',
    '2C Control #2': '../NoRibs_2c_Control201218/',
    '2C H3.1/3.2 KD #1': '../NoRibs_2c_H3_1_3_2KD201209/',
    '2C H3.1/3.2 KD #2': '../NoRibs_2c_H3_1_3_2KD201218/',
    '2C DNA(-) #1': '../NoRibs_2cAphidicolin201209/',
    '2C DNA(-) #2': '../NoRibs_2c_Aphidicoilin201218/',
}
filename = 'Aligned.sortedByCoord.out.uniq_collated.template.sorted.bed.cpm.window_1k_mean.bedGraph.gz'
data_dict = {
    k: pd.read_table(v + filename, engine='c', header=None, names=['chr','start','end','value'], dtype=bedgraph_dtypes)
    for k, v in data_dirs.items()
    }

plot_CPMs(data_dict)

#d['value'].plot.hist(bins=1000, logy=True, cumulative=False)
#plt.show()
#d.loc[d['value'] < 1, 'value'].plot.hist(bins=10, logy=True, cumulative=False)
#plt.show()

positive_samples = ['2C H3.1/3.2 KD #1',
    '2C H3.1/3.2 KD #2',
    '2C DNA(-) #1',
    '2C DNA(-) #2']

negative_samples = ['MII #1', 'MII #2']

def compare_cpm(data: dict, positive: list[str], negative: list[str]):
    print(positive)

#compare_cpm(10, 'bar', ['hoge'])
