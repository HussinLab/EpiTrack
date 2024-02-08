import os
import os.path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import seaborn as sns
import math
#from multiprocessing import Process, Pipe
from Bio import SeqIO #have to install Biopython before using this.


def numREADs_vs_numSAMPLES(Full_Frame, TYPE):

    Frame = Full_Frame[ Full_Frame['Peptide'] == TYPE ]

    grouped_Samples = Frame.groupby(['Deletion', 'Positions'])['Sample'].nunique().reset_index(name = 'Num_samples')
    grouped_Reads = Frame.groupby(['Deletion', 'Positions'])['Del_Freq'].mean().reset_index(name = 'MEAN_ReadCount')

    Samples_v_reads = pd.merge(grouped_Samples, grouped_Reads, how = 'outer', on = ['Deletion', 'Positions'])
    Samples_v_reads.sort_values(by = ['Num_samples'], inplace = True, ascending = False)

    #Samples_v_reads['Num_samples'] = Samples_v_reads['Num_samples'].apply(lambda x: np.log2(x))
    #Samples_v_reads['MEAN_ReadCount'] = Samples_v_reads['MEAN_ReadCount'].apply(lambda x: np.log2(x))

    print(Samples_v_reads)


    with sns.axes_style("ticks"):
        fig, ax = plt.subplots(figsize = (6, 3)) #(10,7)
        sns.set(font_scale=1.0)

        sns.scatterplot(x = 'Num_samples', y = 'MEAN_ReadCount', data = Samples_v_reads)
        #sns.swarmplot(x = 'Peptide', y = 'Del_Freq', color = 'grey', data = Present_only) # palette = PALETTE, edgecolor = 'black', s = 140, alpha = .70,   'count' 'HLAprofile1', s = 100,y_jitter = 0.1, alpha = .35#hue = 'Reoccuring', style = 'Mut_Protein', palette = hue_palette, markers = markerStyle, , palette = PALETTE
        
        #ax.set_xlim(-2, 220)
        #ax.set_ylim(-10, 400)

        plt.xscale('log')
        plt.yscale('log')

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(f'{sys.argv[2]}/READS_v_Samples.png') #_CutAxes
        plt.clf() #round


def Bar_Graphs(Results, read_cutoff):
    grouped_Deletions = Results.groupby(['FRAME'])['Deletion'].nunique().reset_index(name = 'Del_Count')
    grouped_Deletions.loc[ len(grouped_Deletions.index) ] = ['All_Deletions', grouped_Deletions['Del_Count'].sum() ]
    

    grouped_Samples = Results.groupby(['FRAME'])['Sample'].nunique().reset_index(name = 'Samp_Count')
    grouped_Samples.loc[ len(grouped_Samples.index) ] = ['All_Samples', grouped_Samples['Samp_Count'].sum() ]

    grouped_Deletions.to_csv(f'{sys.argv[2]}/grouped_Deletions_{read_cutoff}.csv')
    grouped_Samples.to_csv(f'{sys.argv[2]}/grouped_Samples_{read_cutoff}.csv')


    print(grouped_Deletions)
    print(grouped_Samples)

    with sns.axes_style("ticks"):
        fig, axes = plt.subplots(1, 2, sharey = True, figsize = (10, 3)) #(10,7)
        sns.set(font_scale=1.0)

        sns.barplot(x = 'Samp_Count', y = 'FRAME', data = grouped_Samples, ax = axes[0])
        sns.barplot(x = 'Del_Count', y = 'FRAME', data = grouped_Deletions, ax = axes[1])
        
        axes[0].invert_xaxis()
        axes[0].yaxis.tick_right()

        #axes[1].invert_xaxis()
        #axes[1].yaxis.tick_right()

        axes[0].set(ylabel = None)
        axes[1].set(ylabel = None)

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(f'{sys.argv[2]}/Del_Sample_Stats_{read_cutoff}.png')
        plt.clf() #round
    
    #import matplotlib.pyplot as plt
    Results['Pos_Del'] = Results['Positions'].astype(str) + Results['Deletion']
    grouped_Samples_deletions = Results.groupby(['FRAME', 'Pos_Del'])['Sample'].nunique().reset_index(name = 'Samp_Count') #'Positions', 'Deletion'
    

    grouped_Samples_deletions_frame0 = grouped_Samples_deletions[grouped_Samples_deletions['FRAME'] == 'Frame+0']
    
    grouped_Samples_deletions_frame0.sort_values(by = ['Samp_Count'], inplace = True, ascending = False)
    #Results[]
    grouped_Samples_deletions_frame0.at[grouped_Samples_deletions_frame0.index.values[6:], 'Pos_Del'] = '' #frame_peptides_only.index.values[-1]

    labels_frame0 = grouped_Samples_deletions_frame0['Pos_Del'].values #'Frogs', 'Hogs', 'Dogs', 'Logs'  Deletion
    sizes_frame0 = grouped_Samples_deletions_frame0['Samp_Count'].values  #[15, 30, 45, 10]
    print(f'{labels_frame0}\n{labels_frame0}')
    fig, ax = plt.subplots(figsize = (3.5, 3.5))
    ax.pie(sizes_frame0, labels=labels_frame0)
    plt.savefig(f'{sys.argv[2]}/PIECHART_frame0_{read_cutoff}.png')
    plt.clf() #round
    #########
    #########grouped_Samples_deletions = Results.groupby(['FRAME', 'Pos_Del'])['Sample'].nunique().reset_index(name = 'Samp_Count')
    
    grouped_Samples_deletions_frame1 = grouped_Samples_deletions[grouped_Samples_deletions['FRAME'] == 'Frame+1']
    
    grouped_Samples_deletions_frame1.sort_values(by = ['Samp_Count'], inplace = True, ascending = False)
    #Results[]
    grouped_Samples_deletions_frame1.at[grouped_Samples_deletions_frame1.index.values[6:], 'Pos_Del'] = '' #frame_peptides_only.index.values[-1]
    
    labels_frame1 = grouped_Samples_deletions_frame1['Pos_Del'].values #'Frogs', 'Hogs', 'Dogs', 'Logs'
    sizes_frame1 = grouped_Samples_deletions_frame1['Samp_Count'].values  #[15, 30, 45, 10]
    print(f'{labels_frame1}\n{labels_frame1}')
    fig, ax = plt.subplots(figsize = (3.5, 3.5))
    ax.pie(sizes_frame1, labels=labels_frame1)
    plt.savefig(f'{sys.argv[2]}/PIECHART_frame1_{read_cutoff}.png')
    plt.clf() #round
    ##########
    ##########grouped_Samples_deletions = Results.groupby(['FRAME', 'Pos_Del'])['Sample'].nunique().reset_index(name = 'Samp_Count')
    
    grouped_Samples_deletions_frame2 = grouped_Samples_deletions[grouped_Samples_deletions['FRAME'] == 'Frame+2']
    
    grouped_Samples_deletions_frame2.sort_values(by = ['Samp_Count'], inplace = True, ascending = False)
    #Results[]
    grouped_Samples_deletions_frame2.at[grouped_Samples_deletions_frame2.index.values[6:], 'Pos_Del'] = '' #frame_peptides_only.index.values[-1]
    
    labels_frame2 = grouped_Samples_deletions_frame2['Pos_Del'].values #'Frogs', 'Hogs', 'Dogs', 'Logs'
    sizes_frame2 = grouped_Samples_deletions_frame2['Samp_Count'].values  #[15, 30, 45, 10]
    print(f'{labels_frame2}\n{labels_frame2}')
    fig, ax = plt.subplots(figsize = (3.5, 3.5))
    ax.pie(sizes_frame2, labels=labels_frame2)
    plt.savefig(f'{sys.argv[2]}/PIECHART_frame2_{read_cutoff}.png')
    plt.clf() #round




#######Sample_Count_FirstWave = 27888

Results = pd.read_csv(f'{sys.argv[1]}', index_col=0)

Results = Results[ Results['Del_Freq'] > 1 ]

Results = Results[ Results["Peptide"] != 'OUT_OF_RANGE' ]

numREADs_vs_numSAMPLES(Results, 'Present')
numREADs_vs_numSAMPLES(Results, 'Absent')

Bar_Graphs(Results, 'AboveOneRead')


#Grouped_count = Results.groupby(['Peptide'])['Peptide'].count().reset_index()
Grouped = Results['Peptide'].value_counts()
print(Grouped)
Grouped = Grouped.to_frame().reset_index()
print(Grouped)
Grouped.loc[ len(Grouped.index) ] = ['All_Samples', len(Grouped.index)]  #27888

with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.barplot(x = 'index', y = 'Peptide', data = Grouped)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/Del_Count.png')
    plt.clf() #round


Grouped_average_freq = Results.groupby(['Peptide'])['Del_Freq'].mean().reset_index()

Present_only = Results[ Results['Peptide'] == 'Present' ]

print(Grouped_average_freq)

with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.violinplot(x = 'Peptide', y = 'Del_Freq', data = Present_only)
    #sns.swarmplot(x = 'Peptide', y = 'Del_Freq', color = 'grey', data = Present_only) # palette = PALETTE, edgecolor = 'black', s = 140, alpha = .70,   'count' 'HLAprofile1', s = 100,y_jitter = 0.1, alpha = .35#hue = 'Reoccuring', style = 'Mut_Protein', palette = hue_palette, markers = markerStyle, , palette = PALETTE
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/Deletion_frequency.png')
    plt.clf() #round




LENGTH_GROUPED = Results.groupby(['FRAME', 'Del_Length']).size().reset_index(name = 'Length_count')

COUNT_GROUPED = Results.groupby(['FRAME', 'Del_Freq']).size().reset_index(name = 'Freq_count')

print(f'LENGTH GROUPED: {LENGTH_GROUPED}')

print(f'COUNT GROUPED: {COUNT_GROUPED}')

with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.lineplot(x = 'Del_Length', y = 'Length_count', hue = 'FRAME', data = LENGTH_GROUPED)
    #sns.swarmplot(x = 'Peptide', y = 'Del_Freq', color = 'grey', data = Present_only) # palette = PALETTE, edgecolor = 'black', s = 140, alpha = .70,   'count' 'HLAprofile1', s = 100,y_jitter = 0.1, alpha = .35#hue = 'Reoccuring', style = 'Mut_Protein', palette = hue_palette, markers = markerStyle, , palette = PALETTE
    
    #plt.xscale('log')
    #plt.yscale('log')


    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/Length_distribution.png')
    plt.clf() #round

with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.lineplot(x = 'Del_Freq', y = 'Freq_count', hue = 'FRAME', data = COUNT_GROUPED)
    #sns.swarmplot(x = 'Peptide', y = 'Del_Freq', color = 'grey', data = Present_only) # palette = PALETTE, edgecolor = 'black', s = 140, alpha = .70,   'count' 'HLAprofile1', s = 100,y_jitter = 0.1, alpha = .35#hue = 'Reoccuring', style = 'Mut_Protein', palette = hue_palette, markers = markerStyle, , palette = PALETTE
    
    plt.xscale('log')
    plt.yscale('log')


    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/Frequency_distribution.png')
    plt.clf() #round


FRAME_Specific = Results.groupby(['FRAME'])['Deletion'].nunique().reset_index(name = 'Unique_Deletions')
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.barplot(x = 'FRAME', y = 'Unique_Deletions', data = FRAME_Specific)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_PerFrame.png')
    plt.clf() #round

FRAME_Specific_patients = Results.groupby(['FRAME'])['Sample'].nunique().reset_index(name = 'Unique_Deletions')
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.barplot(x = 'FRAME', y = 'Unique_Deletions', data = FRAME_Specific_patients)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_PerFrame_per_patients.png')
    plt.clf() #round

FRAME_Specific_patients = Results.groupby(['FRAME'])['Deletion'].count().reset_index(name = 'Unique_Deletions')
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.barplot(x = 'FRAME', y = 'Unique_Deletions', data = FRAME_Specific_patients)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_PerFrame_Total_count.png')
    plt.clf() #round




FRAME_ABOVE2 = Results.groupby(['FRAME', 'Del_Length'])['Sample'].nunique().reset_index(name = 'Unique_Deletions_Lengths')
FRAME_ABOVE2.to_csv(f'{sys.argv[2]}/FRAME_Mutation_COunt_atLeast_2.csv')
Length_Summary = Results.groupby(['FRAME']).agg({'Del_Length': ['mean', 'min', 'max']}).reset_index()
print(f'Stat_Summary" {Length_Summary}')
Length_Summary.to_csv(f'{sys.argv[2]}/DELETION_Length_Stat_summary.csv')
Length_Summary = Length_Summary.T.reset_index()
Length_Summary.drop(['level_0'], inplace = True, axis = 1)
Length_Summary.rename(columns = {'level_1': 'Categories', 2: 'FRAME+0', 3: 'FRAME+1', 4: 'FRAME+2'})
print(f'Length Summ: {Length_Summary}')
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.boxplot(x = 'FRAME', y = 'Del_Length', data = Results)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_Lengths_SUmmary_FIgure.png')
    plt.clf() #round

with sns.axes_style("ticks"):
    fig, axes = plt.subplots(3, 1, sharey = True, figsize = (10, 6)) #(10,7)
    sns.set(font_scale=1.0)

    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE2[ FRAME_ABOVE2['FRAME'] == "Frame+0" ], ax = axes[0])
    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE2[ FRAME_ABOVE2['FRAME'] == "Frame+1" ], ax = axes[1])
    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE2[ FRAME_ABOVE2['FRAME'] == "Frame+2" ], ax = axes[2])
    
    #plt.xscale('log')
    #plt.yscale('log')
    axes[0].set(yscale = 'log')
    axes[1].set(yscale = 'log')
    axes[2].set(yscale = 'log')

    #######ax.bar_label(ax.containers[0], labels=groupedvalues['total_bill'])
    FRAME0 = FRAME_ABOVE2[ FRAME_ABOVE2['FRAME'] == "Frame+0"]
    axes[0].bar_label(axes[0].containers[0], labels=FRAME0['Unique_Deletions_Lengths'])
    FRAME1 = FRAME_ABOVE2[ FRAME_ABOVE2['FRAME'] == "Frame+1"]
    axes[1].bar_label(axes[1].containers[0], labels=FRAME1['Unique_Deletions_Lengths'])
    FRAME2 = FRAME_ABOVE2[ FRAME_ABOVE2['FRAME'] == "Frame+2"]
    axes[2].bar_label(axes[2].containers[0], labels=FRAME2['Unique_Deletions_Lengths'])
    
    axes[0].set_title('Frame+0')
    axes[1].set_title('Frame+1')
    axes[2].set_title('Frame+2')

    axes[0].set(xlabel = None)
    axes[1].set(xlabel = None)
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_PerFrame_Length_Distributions.png')
    plt.clf() #round
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.violinplot(x = 'FRAME', y = 'Del_Freq', data = Results)
    #plt.xscale('log')
    plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_PerFrame_Freq_Distributions.png')
    plt.clf() #round


Results = Results[ Results['Del_Freq'] > 50 ]

Bar_Graphs(Results, 'AboveFIFTYRead')

FRAME_ABOVE50 = Results.groupby(['FRAME', 'Del_Length'])['Sample'].nunique().reset_index(name = 'Unique_Deletions_Lengths')

with sns.axes_style("ticks"):
    fig, axes = plt.subplots(3, 1, figsize = (10, 6)) #(10,7)
    sns.set(font_scale=1.0)

    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE50[ FRAME_ABOVE50['FRAME'] == "Frame+0" ], ax = axes[0])
    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE50[ FRAME_ABOVE50['FRAME'] == "Frame+1" ], ax = axes[1])
    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE50[ FRAME_ABOVE50['FRAME'] == "Frame+2" ], ax = axes[2])
    
    #plt.xscale('log')
    #plt.yscale('log')
    #axes[0].set(yscale = 'log')
    #axes[1].set(yscale = 'log')
    #axes[2].set(yscale = 'log')
    
    axes[0].set_title('Frame+0')
    axes[1].set_title('Frame+1')
    axes[2].set_title('Frame+2')

    axes[0].set(xlabel = None)
    axes[1].set(xlabel = None)
    

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_PerFrame_Length_Distributions_above50.png')
    plt.clf() #round
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.violinplot(x = 'FRAME', y = 'Del_Freq', data = Results)
    #plt.xscale('log')
    plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_PerFrame_Freq_Distributions_above50.png')
    plt.clf() #round


Results = Results[ Results['Del_Freq'] > 100 ]

Bar_Graphs(Results, 'AboveONEHUNDREDRead')

FRAME_ABOVE100 = Results.groupby(['FRAME', 'Del_Length'])['Sample'].nunique().reset_index(name = 'Unique_Deletions_Lengths')


with sns.axes_style("ticks"):
    fig, axes = plt.subplots(3, 1, figsize = (10, 6)) #(10,7)
    sns.set(font_scale=1.0)

    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE100[ FRAME_ABOVE100['FRAME'] == "Frame+0" ], ax = axes[0])
    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE100[ FRAME_ABOVE100['FRAME'] == "Frame+1" ], ax = axes[1])
    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE100[ FRAME_ABOVE100['FRAME'] == "Frame+2" ], ax = axes[2])
    
    #plt.xscale('log')
    #plt.yscale('log')
    #axes[0].set(yscale = 'log')
    #axes[1].set(yscale = 'log')
    #axes[2].set(yscale = 'log')
    
    axes[0].set_title('Frame+0')
    axes[1].set_title('Frame+1')
    axes[2].set_title('Frame+2')

    axes[0].set(xlabel = None)
    axes[1].set(xlabel = None)
    

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_PerFrame_Length_Distributions_above100.png')
    plt.clf() #round
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.violinplot(x = 'FRAME', y = 'Del_Freq', data = Results)
    #plt.xscale('log')
    plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_PerFrame_Freq_Distributions_100.png')
    plt.clf() #round



Results = Results[ Results['Del_Freq'] > 500 ]

Bar_Graphs(Results, 'AboveFIVEHUNDREDRead')

FRAME_ABOVE100 = Results.groupby(['FRAME', 'Del_Length'])['Sample'].nunique().reset_index(name = 'Unique_Deletions_Lengths')


with sns.axes_style("ticks"):
    fig, axes = plt.subplots(3, 1, figsize = (10, 6)) #(10,7)
    sns.set(font_scale=1.0)

    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE100[ FRAME_ABOVE100['FRAME'] == "Frame+0" ], ax = axes[0])
    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE100[ FRAME_ABOVE100['FRAME'] == "Frame+1" ], ax = axes[1])
    sns.barplot(x = 'Del_Length', y = 'Unique_Deletions_Lengths', data = FRAME_ABOVE100[ FRAME_ABOVE100['FRAME'] == "Frame+2" ], ax = axes[2])
    
    #plt.xscale('log')
    #plt.yscale('log')
    #axes[0].set(yscale = 'log')
    #axes[1].set(yscale = 'log')
    #axes[2].set(yscale = 'log')
    
    axes[0].set_title('Frame+0')
    axes[1].set_title('Frame+1')
    axes[2].set_title('Frame+2')

    axes[0].set(xlabel = None)
    axes[1].set(xlabel = None)
    

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_PerFrame_Length_Distributions_above500.png')
    plt.clf() #round
with sns.axes_style("ticks"):
    fig, ax = plt.subplots(figsize = (6, 4.5)) #(10,7)
    sns.set(font_scale=1.0)

    sns.violinplot(x = 'FRAME', y = 'Del_Freq', data = Results)
    #plt.xscale('log')
    plt.yscale('log')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(f'{sys.argv[2]}/DELETION_PerFrame_Freq_Distributions_500.png')
    plt.clf() #round

