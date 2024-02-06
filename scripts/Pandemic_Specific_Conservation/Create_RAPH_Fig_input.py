import os
import os.path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import seaborn as sns
import math

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
    #elif x > 0.0:
    #    return '#b35252' #salmon
    #elif (x < 0.0) and (x >= -1):
    #    return '#a0b9ff' #light blue
    elif (x < -1) and (x >= -5):
        return '#5a00e1' #heavier blue/purple ish
    elif (x < -5):
        return '#010048' #drk blue
    else:
        return '#f2e5d5' #off white

def GET_COLOR_PALETTE_EPITOPE_LOSS(frame, COLORSCHEME):
    ####Same as below, but add column for epitope loss, sort by epitope loss, and do exact same but top 5 eitope loss instead of top five most abuncat
    ####Also cut y axis
    print(f'frame before palette and stuff: {frame}')

    ###### rename index DATE, only when do merge with epitope loss ######
    #frame.rename(columns = {'index': 'DATE'}, inplace = True)
    ###### rename index DATE, only when do merge with epitope loss ######

    frame_peptides_only = frame.drop(['DATE', 'sum', 'sum_no_WT'], axis = 1)
    frame_peptides_only = frame_peptides_only.T
    frame_peptides_only['sum'] = frame_peptides_only.sum(axis = 1)

    EPITOPE_Loss = pd.read_csv(f'{sys.argv[3]}Calculated_Epitope_Loss.csv')

    frame_peptides_only.reset_index(inplace = True)
    frame_peptides_only.rename(columns = {'index': 'Raph_IDs'}, inplace = True)
    frame_peptides_only = pd.merge(frame_peptides_only, EPITOPE_Loss, how = 'left', on = ['Raph_IDs'])
    print(f'This is frame:::::: {frame_peptides_only}')
    frame_peptides_only.set_index('Raph_IDs', inplace = True)

    ####frame_peptides_only = frame_peptides_only.rename(columns={"index":"New_ID"})
    frame_peptides_only['NUMBERS'] = range(1, 1 + len(frame_peptides_only))
    
    #######################frame_peptides_only['color'] = frame_peptides_only['Diff'].apply(lambda x: Assess_diff(x))
    
    ########Select only top epitopes#######
    
    '''
    if len(frame_peptides_only.index.values) >= 5:
        frame_peptides_only.loc[5:,'color'] = '#e9edee'
    else:
        frame_peptides_only.loc[4:,'color'] = '#e9edee'
    '''
    ##############################frame_peptides_only.sort_values(by = 'Diff', ascending = True, inplace = True)
    
    if COLORSCHEME == 'NETMHCPAN':
        frame_peptides_only.sort_values(by = 'Diff', ascending = False, inplace = True)
        frame_peptides_only['color'] = frame_peptides_only['Diff'].apply(lambda x: Assess_diff(x))
    
    elif COLORSCHEME == 'Top_Frequency':
        #frame_peptides_only.sort_values(by = 'sum', ascending = False, inplace = True)
        #SORTED_SUM = frame_peptides_only['sum'].values
        #print(f'Sorted sums: {SORTED_SUM}')
        print(f'This is frame sorted by sum no color:\n{frame_peptides_only}')

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

    '''
    ###frame_peptides_only.sort_values(by = 'sum', ascending = False, inplace = True)
    frame_peptides_only.sort_values(by = 'Diff', ascending = True, inplace = True)
    frame_peptides_only['color'] = '#e9edee' #'#e9edee'

    print(f'frame top 5 colors: {frame_peptides_only}')

    #frame_peptides_only.at[frame_peptides_only.index.values[0], 'color'] = '#8d989c'
    
    frame_peptides_only.at[frame_peptides_only.index.values[0], 'color'] = '#490624'
    frame_peptides_only.at[frame_peptides_only.index.values[1], 'color'] = '#8f368d'
    frame_peptides_only.at[frame_peptides_only.index.values[2], 'color'] = '#0072c6'
    frame_peptides_only.at[frame_peptides_only.index.values[3], 'color'] = '#27b0ed'
    if len(frame_peptides_only.index.values) >= 5:
        frame_peptides_only.at[frame_peptides_only.index.values[4], 'color'] = '#008080'

    frame_peptides_only['LABELS'] = ''
    
    #frame_peptides_only.at[frame_peptides_only.index.values[0], 'LABELS'] = frame_peptides_only.index.values[0]
    frame_peptides_only.at[frame_peptides_only.index.values[0], 'LABELS'] = frame_peptides_only.index.values[0]
    frame_peptides_only.at[frame_peptides_only.index.values[1], 'LABELS'] = frame_peptides_only.index.values[1]
    frame_peptides_only.at[frame_peptides_only.index.values[2], 'LABELS'] = frame_peptides_only.index.values[2]
    frame_peptides_only.at[frame_peptides_only.index.values[3], 'LABELS'] = frame_peptides_only.index.values[3]
    if len(frame_peptides_only.index.values) >= 5:
        frame_peptides_only.at[frame_peptides_only.index.values[4], 'LABELS'] = frame_peptides_only.index.values[4]

    
    print(f'frame top 5 colors: {frame_peptides_only}')

    frame_peptides_only.sort_values(by = 'NUMBERS', inplace = True)
    '''

    print(f'frame top 5 colors: {frame_peptides_only}')

    color_palette = frame_peptides_only['color'].values
    LABELS = frame_peptides_only['LABELS'].values  ####frame_peptides_only['LABELS'].values frame_peptides_only.index.values
    print(f'This is colors and LABELS: {color_palette}\n{LABELS}')
    return color_palette, LABELS


def GET_COLOR_PALETTE(frame):
    frame_peptides_only = frame.drop(['DATE', 'sum'], axis = 1)
    frame_peptides_only = frame_peptides_only.T
    frame_peptides_only['sum'] = frame_peptides_only.sum(axis = 1)

    frame_peptides_only.reset_index(inplace = True)
    print(f'This is frame:::::: {frame_peptides_only}')

    ####frame_peptides_only = frame_peptides_only.rename(columns={"index":"New_ID"})
    frame_peptides_only['NUMBERS'] = range(1, 1 + len(frame_peptides_only))

    frame_peptides_only.sort_values(by = 'sum', ascending = False, inplace = True)
    frame_peptides_only['color'] = '#e9edee' #'#e9edee'

    print(f'frame top 5 colors: {frame_peptides_only}')

    #frame_peptides_only.at[frame_peptides_only.index.values[0], 'color'] = '#8d989c'
    
    frame_peptides_only.at[frame_peptides_only.index.values[0], 'color'] = '#490624'
    frame_peptides_only.at[frame_peptides_only.index.values[1], 'color'] = '#8f368d'
    frame_peptides_only.at[frame_peptides_only.index.values[2], 'color'] = '#0072c6'
    frame_peptides_only.at[frame_peptides_only.index.values[3], 'color'] = '#27b0ed'
    if len(frame_peptides_only.index.values) >= 5:
        frame_peptides_only.at[frame_peptides_only.index.values[4], 'color'] = '#008080'

    frame_peptides_only['LABELS'] = ''
    
    #frame_peptides_only.at[frame_peptides_only.index.values[0], 'LABELS'] = frame_peptides_only.index.values[0]
    frame_peptides_only.at[frame_peptides_only.index.values[0], 'LABELS'] = frame_peptides_only.index.values[0]
    frame_peptides_only.at[frame_peptides_only.index.values[1], 'LABELS'] = frame_peptides_only.index.values[1]
    frame_peptides_only.at[frame_peptides_only.index.values[2], 'LABELS'] = frame_peptides_only.index.values[2]
    frame_peptides_only.at[frame_peptides_only.index.values[3], 'LABELS'] = frame_peptides_only.index.values[3]
    if len(frame_peptides_only.index.values) >= 5:
        frame_peptides_only.at[frame_peptides_only.index.values[4], 'LABELS'] = frame_peptides_only.index.values[4]

    '''
    print(frame_peptides_only.iloc[0].at['NUMBERS'])
    ##frame_peptides_only.iloc[0].at['color'] = '#490624'
    frame_peptides_only.at[frame_peptides_only.index.values[0], 'color'] = '#490624'
    print(frame_peptides_only.iloc[1].at['NUMBERS'])
    #####frame_peptides_only.iloc[1].at['color'] = '#8f368d'
    frame_peptides_only['color'][frame_peptides_only.index.values[1]] = '#8f368d'
    #print(frame_peptides_only.loc[2, 'NUMBERS']) #.at['NUMBERS'])
    #frame_peptides_only.loc[2, 'color'] = '#0072c6'
    print(frame_peptides_only.iloc[3].at['NUMBERS'])
    frame_peptides_only.iloc[3].at['color'] = '#27b0ed'
    print(frame_peptides_only.iloc[4].at['NUMBERS'])
    frame_peptides_only.iloc[4].at['color'] = '#008080'
    '''
    print(f'frame top 5 colors: {frame_peptides_only}')

    frame_peptides_only.sort_values(by = 'NUMBERS', inplace = True)


    print(f'frame top 5 colors: {frame_peptides_only}')

    color_palette = frame_peptides_only['color'].values
    LABELS = frame_peptides_only['LABELS'].values

    return color_palette, LABELS



def print_lineplot(PROPORTION_frame, EPITOPE, COLORSCHEME):
    
    COLUMNS = list(PROPORTION_frame.columns)
    #COLUMNS_Dict = reformat(COLUMNS, frame)
    COLUMNS_LIST = reformat(COLUMNS, PROPORTION_frame)
    print(f'COLUMNS_DICT: \n{COLUMNS_LIST}')
    frame = PROPORTION_frame.reset_index()
    frame.drop(['index'], inplace = True, axis = 1)
    ##############################frame['sum'] = frame['sum'].apply(lambda x: np.log2(x))

    #######PALETTE, GOOD_LABELS = GET_COLOR_PALETTE(frame)

    PALETTE, GOOD_LABELS = GET_COLOR_PALETTE_EPITOPE_LOSS(frame, COLORSCHEME) ## Add code

    print(f'GOOD_LABELS: {GOOD_LABELS}')
    
    frame['sum_CUMSUM'] = frame['sum'].cumsum()
    frame['sum_no_WT_CUMSUM'] = frame['sum_no_WT'].cumsum()

    frame.to_csv(f'{sys.argv[2]}{EPITOPE}.csv')
    print(f'Check frame: {frame}')
    
    with sns.axes_style("white"):
        fig, axes = plt.subplots(2, 1, sharex = 'col', gridspec_kw={'height_ratios': [0.35, 1]}, figsize= (17, 6)) # , gridspec_kw={'width_ratios':[1,0.02]}, figsize=(14, 12) 'height_ratios':[0.35, 1, 1], 
        sns.set(font_scale=1.5) #, rc = {lines.linewidth = .5}
        #sns.set_style('white')
        #sns.set_context('notebook', font_scale = 10)
        sns.despine(top = True, right = True, left = True, bottom = True)

        
        ####sns.lineplot(x ='DATE', y = 'count', hue = 'Alternative_peptide', data = frame, legend = False) #, linewidth = 4 color = 'black' color = 'black' , edgecolor = 'black', alpha = 0.7, s = 160, legend = False
        axes[1].stackplot(frame[ 'DATE' ], COLUMNS_LIST, colors = PALETTE, labels = GOOD_LABELS) #'DATE' COLUMNS[0]    COLUMNS_Dict.values()
        #####sns.barplot(x = 'DATE', y = 'sum', data = frame, color = 'grey', ax = axes[0])
        sns.lineplot(x = 'DATE', y = 'sum_CUMSUM', data = frame, color = 'blue', ax = axes[0]) #_sum_not_log 'sum'
        ########sns.lineplot(x = 'DATE', y = 'sum_no_WT_CUMSUM', data = frame, color = 'red', ax = axes[0]) #'sum_no_WT'
        

        plt.xticks(rotation=45, ha = 'right')

        for ind, label in enumerate(axes[1].get_xticklabels()):
            if ind % 5 == 0:  # every 10th label is kept
                label.set_visible(True)
            else:
                label.set_visible(False)
        
        ##axes[0].set_ylim(0, 1000000)
        
        if EPITOPE in ['RANNTKGSL_perm']:
            axes[1].set_ylim(0, 0.2)
        else:
            axes[1].set_ylim(0, 0.05)
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



def Year_add_month(x):
    if x == '2020':
        return '2020-10'
    elif x == '2021':
        return '2021-10'
    elif x == '2022':
        return '2022-10'
    elif x == '2023':
        return '2023-10'
    else:
        return x

def Month_only(x):
    #print(x)
    #return x[ : x.find( '-', x.find('-') ) ]
    
    if x.count('-') == 2:
        if x[-2] == '-':
            return x[: -2 ]
        else:
            return x[: -3 ]
    else:
        return x
    

def modify_string(x):
    x = x.replace(' ', ',')
    Pos = x.rfind(',', 0, 7)
    x = x[Pos+1:]
    x = x.split(',')
    return x

def Remove_00(x):
    if x[-2 :] == '00':
        return np.nan
    else:
        return x


def Single_digit_month(x):
    print(x)
    print(x[-2])
    if x[-2] == '-':
        return x[:-1] + '0' + x[-1]
    else:
        return x


def Get_alternative_peptide_Total(frame):
    frame_Transposed = frame.set_index('DATE')
    frame_Transposed = frame_Transposed.T

    frame_Transposed['sum'] = frame_Transposed.sum(axis = 1)

    frame_Transposed.to_csv(f'{sys.argv[2]}.csv')


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
    
    ########frame_T = frame.set_index('DATE') #frame_numsOnly
    ##########frame_T = frame_numsOnly.T
    
    ####Sort columns by peptide total#####
    
    frame_T = pd.read_csv(file, sep='\t', index_col=0)
    
    ######### REMOVE REF ###########
    
    
    column_to_remove = []
    for column in frame_T.columns.values:
        if (column == 'other'): #('_ref' in column) or 
            column_to_remove.append(column)
    frame_T.drop(column_to_remove, axis = 1, inplace = True)

    

    ######### REMOVE REF ###########


    ####### REMOVE REF



    print(f'no index name:\n{frame_T}')
    frame_T = frame_T.reset_index()
    print(f'no index name:\n{frame_T}')
    print(frame_T.index.values)
    frame_T = frame_T.rename(columns = {'.': 'DATE'})
    frame_T = frame_T.set_index('DATE')
    print(f' index name:\n{frame_T}')
    
    frame_T = frame_T.T

    frame_T['sum_peptides'] = frame_T.sum(axis = 1)
    frame_T.sort_values(by = ['sum_peptides'], ascending = True, inplace = True)
    
    print(f' sum_peptides:\n{frame_T}')

    frame_T.drop(['sum_peptides'], axis = 1, inplace = True)

    print(f' sum_peptides_after_drop:\n{frame_T}')
    
    '''
    EPITOPE_Loss = pd.read_csv(f'{sys.argv[3]}Calculated_Epitope_Loss.csv')

    EPITOPE_Loss = EPITOPE_Loss[['Raph_IDs', 'Diff']]
    print(EPITOPE_Loss)
    
    temp_frameXX = frame_T.reset_index()
    temp_frameXX.rename(columns = {'index': 'Raph_IDs'}, inplace = True)
    temp_frameXX = pd.merge(temp_frameXX, EPITOPE_Loss, how = 'left', on = ['Raph_IDs'])
    
    frame_T['Diff'] = temp_frameXX['Diff'].values

    frame_T.sort_values(by = ['Diff'], ascending = False, inplace = True)
    frame_T.drop(['Diff'], axis = 1, inplace = True)

    '''

    '''
    print(EPITOPE_Loss[EPITOPE_Loss['Raph_IDs'] == 'ASDLDDFSKQLQ']['Diff'].values[0])
    
    frame_T['Diff'] = frame_T['index'].apply(lambda x: EPITOPE_Loss[EPITOPE_Loss['Raph_IDs'] == x]['Diff'].values[0] if x in EPITOPE_Loss['Raph_IDs'].values)
    print(f'This is frame:::::: {frame_T}')
    frame_T.sort_values(by = ['Diff'], ascending = True, inplace = True)
    print(f'This is frame:::::: {frame_T}')
    frame_T.drop(['Diff'], axis = 1, inplace = True)

    frame_T.set_index('index', inplace= True)


    

    frame_T.reset_index(inplace = True)
    frame_T.rename(columns = {'index': 'Raph_IDs'}, inplace = True)
    frame_T = pd.merge(frame_T, EPITOPE_Loss, how = 'left', on = ['Raph_IDs'])
    print(f'This is frame:::::: {frame_T}')
    frame_T.set_index('Raph_IDs', inplace = True)

    frame_T.sort_values(by = ['Diff'], ascending = True, inplace = True)
    frame_T.drop(['Diff'], inplace = True, axis = 1)
    
    print(f'This is frame:::::: {frame_T}')

    
    '''
    '''
    frame_T['sum_peptides'] = frame_T.sum(axis = 1)
    frame_T.sort_values(by = ['sum_peptides'], ascending = True, inplace = True)
    frame_T.drop(['sum_peptides'], axis = 1, inplace = True)
    '''

    #index_count = len(frame_T.index.values)

    ###################################################frame_T = frame_T.tail( 10 ) #math.ceil( index_count * 0.25 )

    frame_T = frame_T.T

    print(f' sum_peptides_after_transpose:\n{frame_T}')

    ####Sort columns by peptide total#####
    
    #Get sum of mutants only
    count = 0
    
    Short_pep = PEPTIDE[: PEPTIDE.find('_')]
    print(Short_pep)
    
    for column in frame_T.columns.values:
        print(f'count: {count}')
        print(column)
        if Short_pep in column:  ###('_ref' in column): #('_ref' in column) or 
            REFS = frame_T[column].values
            frame_T.drop(column, axis = 1, inplace = True)

            frame_T['sum_no_WT'] = frame_T.sum(axis = 1)
            sum_no_WT = frame_T['sum_no_WT'].values

            frame_T.drop(['sum_no_WT'], axis = 1, inplace = True)

            print(f'count: {count}')
            
            
            if 'APRITFGGP_ref' in column:
                frame_T.insert(loc = 23, column = column, value = REFS)
                print('THIS IS APRITFGGP')
            else:
                frame_T.insert(loc = count, column = column, value = REFS)
                print('THIS IS OTHER')

            #frame_T[column] = REFS

        count += 1

    frame_T['sum'] = frame_T.sum(axis = 1)

    print(f' sum_peptides_after_removing REF:\n{frame_T}')

    COMBINE_ALL_PEPTIDES(frame_T, PEPTIDE)

    print(f'frame_T: pre_division: \n{frame_T}')
    for COLUMN in frame_T.columns.values[:-1]:
        print(f'column: {COLUMN}')
        #####GET RATIO#####
        frame_T[COLUMN] = frame_T[COLUMN]/frame_T['sum']
    
    
        #####GET RATIO#####

    frame_T['sum_no_WT'] = sum_no_WT

    #frame_T = frame_T[ frame_T.columns.values[:-1] ] / frame_T['sum']
    print(f'frame_T: post_division: \n{frame_T}')
    '''  
    #########Check if it adds up to 1###########

    frame_T.drop(['sum'], axis = 1, inplace = True)
    frame_T['sum'] = frame_T.sum(axis = 1)
    print(f'frame_T: Re-sum: \n{frame_T}')
    frame_T.drop(['sum'], axis = 1, inplace = True)

    #########Check if it adds up to 1###########
    '''
    frame_T = frame_T.reset_index() #.T

    ######### Get sum of rows (rows should be peptides here), sort by sum, remove sum.
    


    print(f'frame_T: final: \n{frame_T}')
    #######frame_T = frame_T.T
    print(f'frame_T: final: \n{frame_T}')
    return frame_T






def Prepare_data(filename, FinalFrame):
    print(filename)
    FRAME = pd.read_csv(filename)
    
    FRAME = FRAME.T.reset_index().T
    print(FRAME)
    FRAME[0] = FRAME[0].apply(lambda x: modify_string(x))
    FRAME['count'] = FRAME[0].apply(lambda x: x[0])
    FRAME['DATE'] = FRAME[0].apply(lambda x: x[1])
    FRAME['count'] = FRAME['count'].astype(int)

    FRAME.dropna(inplace = True)    

    FRAME['DATE'] = FRAME['DATE'].apply(lambda x: Year_add_month(x))
    FRAME['DATE'] = FRAME['DATE'].apply(lambda x: Month_only(x))
    FRAME['DATE'] = FRAME['DATE'].apply(lambda x: Remove_00(x))

    FRAME.dropna(subset = ['DATE'], inplace = True)
    FRAME['DATE'] = FRAME['DATE'].apply(lambda x: Single_digit_month(x))

    
    
    #FRAME['count'] = FRAME['count'].astype(int)
    print(FRAME)
    grouped = FRAME.groupby(['DATE'])['count'].sum().reset_index(name = filename) #'count'
    grouped[filename] = grouped[filename].astype(int)
    ########grouped[filename] = np.log2(grouped[filename])
    grouped[filename].fillna(0, inplace = True)
    print(grouped)
    ##########grouped['Alternative_peptide'] = filename
    
    if FinalFrame.empty:
        FinalFrame = grouped
    else:
        FinalFrame = pd.merge(FinalFrame, grouped, how = 'outer', on = ['DATE'])
        FinalFrame.fillna(0, inplace = True)
    print(f'FinalFrame: \n{FinalFrame}')

    return FinalFrame #grouped #FRAME

def getID(x, PEPTIDE):
    if x == PEPTIDE:
        return 'REF'
    else:
        return 'MUT'

def Process_data(FILE, PEPTIDE):
    reducedFILE = FILE.rename(columns = {0:'count', 1:'NUC', 2:'AA'})
    reducedFILE = reducedFILE[['count', 'AA']]
    reducedFILE['count'] = reducedFILE['count'].astype(int)
    reducedFILE['AA'] = reducedFILE['AA'].astype(str)
    reducedFILE = reducedFILE[ reducedFILE['count'] > 10 ]
    grouped = reducedFILE.groupby(['AA'], sort = False)['count'].sum().reset_index(name = 'sum_count')
    print(grouped)

    grouped['ID'] = grouped['AA'].apply(lambda x: getID(x, PEPTIDE))
    print(grouped)

    TOTAL_SUM_EVERYTHING = grouped['sum_count'].sum()

    grouped_mutOnly = grouped[ grouped['ID'] == 'MUT' ]
    
    TOTAL_SUM_MUTONLY = grouped_mutOnly['sum_count'].sum()

    RATIO = TOTAL_SUM_MUTONLY/TOTAL_SUM_EVERYTHING
    print(f'TOTAL_SUM_EVERYTHING: {TOTAL_SUM_EVERYTHING}\nTOTAL_SUM_MUTONLY: {TOTAL_SUM_MUTONLY}\nRATIO: {RATIO}')
    
    return RATIO



AllRatios = []

FullFrame = pd.DataFrame() ##Initiate dataframe
peptides_to_keep = ['KLPDDFTGC', 'TLNDLNETL', 'NAPRITFGGP', 'VPYNMRVI', 'RANNTKGSL', 'GPMVLRGLIT', 'STTTNIVTR', 'TGSNVFQTR', 'HTTDPSFLGR', 'KTIQPRVEK', 'TTDPSFLGRYM', 'PTDNYITTY', 'YLFDESGEFKL', 'LPKEITVAT', 'TTDPSFLGRY', 'RTIKVFTTV']
with open(sys.argv[1], 'r') as File:   #open and save to Results_Files dictionary HLA-specific processed netMHCpan files.
    for line in File:
        if len(line) > 1:
            ###print(sys.argv[3])
            templine = line[:-1]
            print(templine)
            
            PEPTIDE = templine[: templine.find('.')]
            print(PEPTIDE)
            if PEPTIDE in peptides_to_keep:
                ##########if sys.argv[3] not in line[:-1]: #sys.argv[3]  '_'
                #DATA_Files.append( Prepare_data(line[:-1], FullFrame) )
                
                FILE = pd.read_csv(templine, sep=' ')
                
                FILE = FILE.T 
                FILE.reset_index(inplace = True)
                FILE = FILE.T

                print(FILE)
                print(FILE.columns.values)

                measured_ratio = Process_data(FILE, PEPTIDE)
                AllRatios.append([PEPTIDE, measured_ratio])

RATIOS_FRAME = pd.DataFrame(data = AllRatios) #, columns = ['Peptide', 'RATIO']

print(RATIOS_FRAME)

RATIOS_FRAME.to_csv(f'{sys.argv[2]}New_Peptide_list.csv')

Sorted_RATIOS_FRAME = RATIOS_FRAME.sort_values(by = [1], ascending = False)

print(Sorted_RATIOS_FRAME)
            
            
            ####IMPORTANT#######
            #Proportional_FinalFrame = Convert_to_proportion(line[:-1], line[:-9]) #sys.argv[3]FullFrame
            
            #print_lineplot(Proportional_FinalFrame, line[:-9], 'NETMHCPAN')
            ########print_lineplot(Proportional_FinalFrame, line[:-9], 'Top_Frequency')
            ####IMPORTANT#######

            ######FullFrame = Prepare_data(line[:-1], FullFrame) 
            #####DATA_Files.append(line[:-1])
    #Results_file = Results_Files

#ALlResults = pd.concat(DATA_Files)
#ALlResults['count'] = ALlResults['count'].astype(int)
#ALlResults['Alternative_peptide'] = ALlResults['Alternative_peptide'].astype(str)
#FullFrame.sort_values(by = ['DATE'], inplace = True, ascending = True)
#print(f'This is fullframe: \n{FullFrame}')
#print(FullFrame.dtypes)

######ALlResults = ALlResults[ '00' not in ALlResults['DATE'] ]
############ALlResults['count'] = np.log2(ALlResults['count'])

###### Convert to proportion ########
     
######Proportional_FinalFrame = Convert_to_proportion(FullFrame, sys.argv[3])

#####Get_alternative_peptide_Total(FullFrame)

###### Convert to proportion ########

#####print_lineplot(FullFrame)
########print_lineplot(Proportional_FinalFrame)


#Put all peptides on same figure, bring all alternative peptides together!
#####







