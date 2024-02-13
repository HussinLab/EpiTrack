import os
import os.path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import seaborn as sns
#import math

def reformat(COLUMNS, frame):
    #NewDict = {}
    NewDict = []
    for item in COLUMNS:
        if (item == 'DATE') or (item == 'sum') or (item == 'sum_no_WT'): #'DATE' COLUMNS[0]
            continue
        else:
            #####frame[item] = frame[item].astype(int)
            print(f'check types: {item}: \n{frame[item].dtypes}')
            #NewDict[item] = frame[item].values.tolist()
            NewDict.append( frame[item].values.tolist() )
    return NewDict


def Assess_diff(x):
    if x == 0.0:
        return '#FFFFFF'
    elif (x < -1) and (x >= -5):
        return '#5a00e1' #heavier blue/purple ish
    elif (x < -5):
        return '#010048' #drk blue
    else:
        return '#f2e5d5' #off white

def GET_COLOR_PALETTE_EPITOPE_LOSS(frame, COLORSCHEME):
    
    frame_peptides_only = frame.drop(['DATE', 'sum', 'sum_no_WT'], axis = 1)
    frame_peptides_only = frame_peptides_only.T
    frame_peptides_only['sum'] = frame_peptides_only.sum(axis = 1)

    EPITOPE_Loss = pd.read_csv(f'{sys.argv[4]}/scripts/Calculated_Epitope_Loss.csv')

    frame_peptides_only.reset_index(inplace = True)
    frame_peptides_only.rename(columns = {'index': 'Raph_IDs'}, inplace = True)
    frame_peptides_only = pd.merge(frame_peptides_only, EPITOPE_Loss, how = 'left', on = ['Raph_IDs'])
    print(f'This is frame:::::: {frame_peptides_only}')
    frame_peptides_only.set_index('Raph_IDs', inplace = True)

    ####frame_peptides_only = frame_peptides_only.rename(columns={"index":"New_ID"})
    frame_peptides_only['NUMBERS'] = range(1, 1 + len(frame_peptides_only))
    
    
    if COLORSCHEME == 'NETMHCPAN':
        frame_peptides_only.sort_values(by = 'Diff', ascending = False, inplace = True)
        frame_peptides_only['color'] = frame_peptides_only['Diff'].apply(lambda x: Assess_diff(x))
    
    elif COLORSCHEME == 'Top_Frequency':

        frame_peptides_only['color'] = '#e9edee' #'#e9edee'

        frame_peptides_only.at[frame_peptides_only.index.values[-1], 'color'] = '#EAFAF1'   #'#490624'
        if len(frame_peptides_only.index.values) >= 2:
            frame_peptides_only.at[frame_peptides_only.index.values[-2], 'color'] = '#8f368d'
            if len(frame_peptides_only.index.values) >= 3:
                frame_peptides_only.at[frame_peptides_only.index.values[-3], 'color'] = '#0072c6'
                if len(frame_peptides_only.index.values) >= 4:
                    frame_peptides_only.at[frame_peptides_only.index.values[-4], 'color'] = '#27b0ed'
                    if len(frame_peptides_only.index.values) >= 5:
                        frame_peptides_only.at[frame_peptides_only.index.values[-5], 'color'] = '#008080'
        print(f'This is frame sorted by sum with color:\n{frame_peptides_only}')
    frame_peptides_only['LABELS'] = ''
    
    #frame_peptides_only.at[frame_peptides_only.index.values[0], 'LABELS'] = frame_peptides_only.index.values[0]
    frame_peptides_only.at[frame_peptides_only.index.values[-1], 'LABELS'] = frame_peptides_only.index.values[-1]
    if len(frame_peptides_only.index.values) >= 2:
        frame_peptides_only.at[frame_peptides_only.index.values[-2], 'LABELS'] = frame_peptides_only.index.values[-2]
        if len(frame_peptides_only.index.values) >= 3:
            frame_peptides_only.at[frame_peptides_only.index.values[-3], 'LABELS'] = frame_peptides_only.index.values[-3]
            if len(frame_peptides_only.index.values) >= 4:
                frame_peptides_only.at[frame_peptides_only.index.values[-4], 'LABELS'] = frame_peptides_only.index.values[-4]
                if len(frame_peptides_only.index.values) >= 5:
                    frame_peptides_only.at[frame_peptides_only.index.values[-5], 'LABELS'] = frame_peptides_only.index.values[-5]

    ########Select only top epitopes#######
    
    frame_peptides_only.sort_values(by = 'NUMBERS', ascending = True, inplace = True)
    print(f'frame top 5 colors: {frame_peptides_only}')

    color_palette = frame_peptides_only['color'].values
    LABELS = frame_peptides_only['LABELS'].values  ####frame_peptides_only['LABELS'].values frame_peptides_only.index.values
    print(f'This is colors and LABELS: {color_palette}\n{LABELS}')
    return color_palette, LABELS


def print_lineplot(PROPORTION_frame, EPITOPE, COLORSCHEME):
    
    COLUMNS = list(PROPORTION_frame.columns) #Get all alternative peptide names
    COLUMNS_LIST = reformat(COLUMNS, PROPORTION_frame) #Reformat is a function that returns date-wise GISAID sequence counts for all alternative peptides (but skips summarizing columns like sum and sum_no_WT)
    frame = PROPORTION_frame.reset_index()
    frame.drop(['index'], inplace = True, axis = 1)
    ##############################frame['sum'] = frame['sum'].apply(lambda x: np.log2(x)) #this line turns the data into log2 if needed
    PALETTE, GOOD_LABELS = GET_COLOR_PALETTE_EPITOPE_LOSS(frame, COLORSCHEME) #Function that acquires labels colors and labels for stackplot
    print(f'GOOD_LABELS: {GOOD_LABELS}')
    
    frame['sum_CUMSUM'] = frame['sum'].cumsum()
    frame['sum_no_WT_CUMSUM'] = frame['sum_no_WT'].cumsum()

    frame.to_csv(f'{sys.argv[2]}{EPITOPE}.csv')
    
    with sns.axes_style("white"):
        fig, axes = plt.subplots(2, 1, sharex = 'col', gridspec_kw={'height_ratios': [0.35, 1]}, figsize= (17, 6)) # , gridspec_kw={'width_ratios':[1,0.02]}, figsize=(14, 12) 'height_ratios':[0.35, 1, 1], 
        sns.set(font_scale=1.5) #, rc = {lines.linewidth = .5}
        #sns.set_style('white')
        #sns.set_context('notebook', font_scale = 10)
        sns.despine(top = True, right = True, left = True, bottom = True)
        axes[1].stackplot(frame[ 'DATE' ], COLUMNS_LIST, colors = PALETTE, labels = GOOD_LABELS) #'DATE' COLUMNS[0]    COLUMNS_Dict.values()
        sns.lineplot(x = 'DATE', y = 'sum_CUMSUM', data = frame, color = 'blue', ax = axes[0]) #_sum_not_log 'sum'
        plt.xticks(rotation=45, ha = 'right')

        for ind, label in enumerate(axes[1].get_xticklabels()):
            if ind % 5 == 0:  # every 10th label is kept
                label.set_visible(True)
            else:
                label.set_visible(False)
        
        ##axes[0].set_ylim(0, 1000000)
        
        ################if EPITOPE in ['RANNTKGSL_perm']:
            ################axes[1].set_ylim(0, 0.2)
        ################else:
            ################axes[1].set_ylim(0, 0.05)
        
        if sys.argv[5].lower() == 'yes':
            axes[1].set_ylim(0, 0.2)
        #elif EPITOPE in ['HTTDPSFLGR_perm']:
        #    axes[1].set_ylim(0, 0.2)
        '''
        elif EPITOPE in ['RANNTKGSL_perm']:
            axes[1].set_ylim(0, 10000)
        
        
        elif EPITOPE in ['AADLDDFSKQLQ_perm', 'RTIKVFTTV_perm', 'STTTNIVTR_perm']:
            axes[1].set_ylim(0, 6000)
        
        else:
            axes[1].set_ylim(0, 3000)
        '''
        #####axes[1].legend(loc='upper left', bbox_to_anchor=(1.25, 0.5))
        ###axes[1].set_ylim(0, 50000)
        #################################axes[1].legend().remove()
        axes[1].legend(loc='upper left', bbox_to_anchor=(1.25, 0.5))
        #plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        ###plt.ylim(0, 1200000)
        #plt.ylim(0, LIMIT)
        #plt.ylabel(XLabel)
        plt.tight_layout()
        plt.savefig(f'{sys.argv[2]}{EPITOPE}_{COLORSCHEME}.png', dpi = 300)
        plt.clf()


def COMBINE_ALL_PEPTIDES(frame_T, PEPTIDE):
    print('COMBINE_ALL_PEPTIDES')
    reduced_frame = frame_T.reset_index()
    ###### rename index DATE, only when do merge with epitope loss ######
    #reduced_frame.rename(columns = {'index': 'DATE'}, inplace = True)
    ###### rename index DATE, only when do merge with epitope loss ######
    print(reduced_frame)
    reduced_frame = reduced_frame[['DATE', 'sum']].copy()
    reduced_frame.rename(columns = {'sum': PEPTIDE}, inplace = True)
    #######reduced_frame.drop([''])

    PATH = f'{sys.argv[2]}AllPeptides_data.csv'

    if os.path.isfile(PATH):
        Allpeptides = pd.read_csv(f'{sys.argv[2]}AllPeptides_data.csv')
        Allpeptides = pd.merge(Allpeptides, reduced_frame, how = 'outer', on = 'DATE')
        Allpeptides.to_csv(f'{sys.argv[2]}AllPeptides_data.csv')
        print('hello found all peptides')
    else:
        reduced_frame.to_csv(f'{sys.argv[2]}AllPeptides_data.csv')
        print('first file')

def Convert_to_proportion(file, PEPTIDE):
    
    frame_T = pd.read_csv(file, sep='\t', index_col=0) #Import peptide file
    
    ######### REMOVE "OTHER" column as well as Reference (Optional) ###########
    column_to_remove = []
    for column in frame_T.columns.values:
        if (column == 'other'): #('_ref' in column) or 
            column_to_remove.append(column)
    frame_T.drop(column_to_remove, axis = 1, inplace = True)
    ######### REMOVE "OTHER" column as well as Reference (Optional) ###########

    #### Here, row = dates, columns = alternative peptides. We want to rename the index and get sum of sequences for each alt. peptide across all dates #######
    frame_T = frame_T.reset_index()
    frame_T = frame_T.rename(columns = {'.': 'DATE'})
    frame_T = frame_T.set_index('DATE') #Name index as DATE
    frame_T = frame_T.T #transpose
    frame_T['sum_peptides'] = frame_T.sum(axis = 1) #Get sum of peptide sequences across each date
    frame_T.sort_values(by = ['sum_peptides'], ascending = True, inplace = True)
    frame_T.drop(['sum_peptides'], axis = 1, inplace = True)
    frame_T = frame_T.T

    count = 0
    
    #PEPTIDE here is the reference peptide
    Short_pep = PEPTIDE[: PEPTIDE.find('_')]
    print(Short_pep)
    
    for column in frame_T.columns.values: #If the reference peptide is in the columns, drop it for subsequent operations
        if Short_pep in column:  
            REFS = frame_T[column].values
            frame_T.drop(column, axis = 1, inplace = True)

            frame_T['sum_no_WT'] = frame_T.sum(axis = 1) #Gert sum of GISAID sequences across dates (without reference)
            sum_no_WT = frame_T['sum_no_WT'].values #Save GISAID total sequences across dates to variable "sum_no_WT"

            frame_T.drop(['sum_no_WT'], axis = 1, inplace = True) #drop column
            
            ##Re-introduce reference peptide data (REFS)
            if 'APRITFGGP_ref' in column:
                frame_T.insert(loc = 23, column = column, value = REFS)
                print('THIS IS APRITFGGP')
            else:
                frame_T.insert(loc = len(frame_T.columns.values), column = column, value = REFS) #count
                print('THIS IS OTHER')
        count += 1

    frame_T['sum'] = frame_T.sum(axis = 1) #get sum across dates

    COMBINE_ALL_PEPTIDES(frame_T, PEPTIDE)
    print(f'frame_T: pre_division: \n{frame_T}')
    for COLUMN in frame_T.columns.values[:-1]:
        print(f'column: {COLUMN}')
        #####GET RATIO#####
        frame_T[COLUMN] = frame_T[COLUMN]/frame_T['sum']
        #####GET RATIO#####

    frame_T['sum_no_WT'] = sum_no_WT #re-introduce the "sum_no_WT" column
    
    frame_T = frame_T.reset_index()

    return frame_T

DATA_Files = []

FullFrame = pd.DataFrame() #Initiate empty dataframe

with open(sys.argv[1], 'r') as File: #Iterate through file containing peptide_specific GISAID data
    for line in File:
        if len(line) > 1:
            
            templine = line[:-1]
            print(templine)            
            Proportional_FinalFrame = Convert_to_proportion(line[:-1], line[:-9]) #Conver to proportion
            print_lineplot(Proportional_FinalFrame, line[:-9], 'NETMHCPAN')
            print_lineplot(Proportional_FinalFrame, line[:-9], 'Top_Frequency')


 







