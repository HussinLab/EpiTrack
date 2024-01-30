import os
import os.path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import seaborn as sns
#import math


def modify_string(x):
    x = x.replace(' ', ',')
    Pos = x.rfind(',', 0, 7)
    x = x[Pos+1:]
    x = x.split(',')
    return x


##### The function below does general data processing to generate lineage-specific summarizing stats #####
def Prepare_data(filename, FinalFrame):
    print(filename)
    FRAME = pd.read_csv(filename)
    
    FRAME = FRAME.T.reset_index().T
    print(FRAME)
    FRAME[0] = FRAME[0].apply(lambda x: modify_string(x))
    FRAME['count'] = FRAME[0].apply(lambda x: x[0])
    FRAME['Lineage'] = FRAME[0].apply(lambda x: x[1])
    FRAME['count'] = FRAME['count'].astype(int)
    FRAME.dropna(inplace = True)    
    FRAME.dropna(subset = ['Lineage'], inplace = True)
    print(FRAME)
    grouped = FRAME.groupby(['Lineage'])['count'].sum().reset_index(name = filename) #'count'
    grouped[filename] = grouped[filename].astype(int)
    ########grouped[filename] = np.log2(grouped[filename])
    grouped[filename].fillna(0, inplace = True)
    print(grouped)
    ##########grouped['Alternative_peptide'] = filename
    
    if FinalFrame.empty:
        FinalFrame = grouped
    else:
        FinalFrame = pd.merge(FinalFrame, grouped, how = 'outer', on = ['Lineage'])
        FinalFrame.fillna(0, inplace = True)
    print(f'FinalFrame: \n{FinalFrame}')

    return FinalFrame #grouped #FRAME

#### Function below simplifies lineage name by categorizing them into VOC notation
def Simplify_lineage(x):
    if ('B.1.617.2' in x) or ('AY' in x):
        return 'Delta'
    elif 'B.1.1.7' in x:
        return 'Alpha'
    elif 'B.1.351' in x:
        return 'Beta'
    elif ('B.1.1.529' in x) or ('XE' in x) or ('EG' in x) or ('HV' in x) or ('BF' in x) or ('BN' in x) or ('BQ' in x) or ('FR' in x) or ('FK' in x) or ('EJ' in x) or ('DV' in x) or ('DF' in x) or ('CJ' in x) or ('CH' in x) or ('BY' in x) or ('BR' in x) or ('BM' in x) or ('BL' in x) or ('BE' in x):
        return 'Omicron'
    elif 'BA.1' in x:
        return 'Omicron' #'BA.1'
    elif 'BA.2' in x:
        return 'Omicron' #'BA.2'
    elif ('BA.4' in x) or ('BA.5' in x):
        return 'Omicron' #'BA.4/5'
    elif 'XBB' in x:
        return 'Omicron' #recombinant (XBB)
    else:
        return 'Other'


DATA_Files = []

FullFrame = pd.DataFrame() ##Initiate dataframe

with open(sys.argv[1], 'r') as File:   #Iterate through peptide_spercific files
    for line in File:
        if len(line) > 1:
            print(sys.argv[3])
            templine = line[:-1]
            print(templine)
            if ('full' not in line[:-1]) and ('nextrainnuc' not in line[:-1]):
                if sys.argv[3] not in line[:-1]: #sys.argv[3]  '_'
                    
                    FullFrame = Prepare_data(line[:-1], FullFrame) 
              
FullFrame.sort_values(by = ['Lineage'], inplace = True, ascending = True)
FullFrame['Simplified lineages'] = FullFrame['Lineage'].apply(lambda x: Simplify_lineage(x)) #Simplify lineage names general VOC caterogies
print(f'This is fullframe: \n{FullFrame}')
FullFrame.drop(['Lineage'], axis = 1, inplace = True)
Columns = FullFrame.columns.values
Grouped = FullFrame.groupby(['Simplified lineages'])[Columns[:-1]].sum().reset_index()

#### Code below if specifically for NAPRI epitope, to reduce size of heatmap
print(Grouped)
print(f'This is REF: {sys.argv[3]}')
if 'NAPRI' in sys.argv[3]:
    ALL_columns = Grouped.columns.values
    Columns_to_keep = ['NAFRITFGGP', 'NGPRITFGGP', 'NALRITFGGP', 'NAPRITFVGP', 'NAPRITFGGP', 'Simplified lineages']
    kept_columns = []
    for column in ALL_columns:
        for column_to_keep in Columns_to_keep:
            if column_to_keep in column:
                kept_columns.append(column)

    Grouped = Grouped[kept_columns]
        
Grouped.set_index('Simplified lineages', inplace = True)
FullFrame.to_csv(f'{sys.argv[2]}.csv') #Print out raw frame used to generate heatmap
Grouped.to_csv(f'{sys.argv[2]}_grouped.csv') #print out simplified frame used to make heatmap


####Print heatmap######
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4)) #(10,7)
    sns.set(font_scale=1.0)

    sns.heatmap(Grouped, cmap="crest", linewidth=.5) #, annot = True
    #plt.xscale('log')
    #plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}_heatmap.png')
    plt.clf() #round








