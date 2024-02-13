import os
import os.path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import seaborn as sns
#import math
import geopandas as gpd
#import pygal
#from pygal.style import Style
#import matplotlib.patches as mpatches


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
    
def Remove_00(x):
    if x[-2 :] == '00':
        return np.nan
    else:
        return x


def Single_digit_month(x):
    #print(x)
    #print(x[-2])
    if x[-2] == '-':
        return x[:-1] + '0' + x[-1]
    else:
        return x


def Get_alternative_peptide_Total(frame):
    frame_Transposed = frame.set_index('DATE')
    frame_Transposed = frame_Transposed.T
    frame_Transposed['sum'] = frame_Transposed.sum(axis = 1)
    frame_Transposed.to_csv(f'{sys.argv[2]}.csv')

def Modify_country_names(x):
    if '_' in x:
        return x.replace('_', ' ')
    elif x == 'USA':
        return 'United States of America' #'United States'
    else:
        return x

def SET_COLOR(x, ALL_COLORS):
    if x == 0:
        return '#FFFFFF' #ALL_COLORS[-1]
    elif x <= 0.002:
        return ALL_COLORS[7] #'blue' #'#EAFAF1'
    elif (x > 0.002) and (x <= 0.01):
        return ALL_COLORS[6] #'red' #'#D4EFDF'
    elif (x > 0.01) and (x <= 0.05):
        return ALL_COLORS[5] #'yellow' #'#A2D9CE'
    elif (x > 0.05) and (x <= 0.1):
        return ALL_COLORS[4] #'black' #'5DADE2'
    elif (x > 0.1) and (x <= 0.2):
        return ALL_COLORS[3] #'brown' #'#2980B9'
    elif (x > 0.2) and (x <= 0.4):
        return ALL_COLORS[2] #'white' #'#7D3C98'
    elif (x > 0.4) and (x <= 0.6):
        return ALL_COLORS[1] #'white' #'#943126'
    elif (x > 0.6):
        return ALL_COLORS[0] #'white' #'#641E16'

def Create_MAP_WholePANDEMIC(Modified_FullFrame, Modified_WT_frame, Subset, PEPNAME):
    print(f'This is full frame dates_not_combined: {Modified_FullFrame}')
    FullFrame_AllDates = Modified_FullFrame.groupby(['Continent', 'Country'])['Sum_across_peptides'].sum().reset_index(name = 'AllDates_count_ALT')
    print(f'This is full frame all dates: {FullFrame_AllDates}')
    WT_AllDates = Modified_WT_frame.groupby(['Continent', 'Country'])['Sum_across_peptides'].sum().reset_index(name = 'AllDates_count_WT')
    MERGED_ALT_WT = pd.merge( FullFrame_AllDates, WT_AllDates, how = 'outer', on = ['Continent', 'Country'] )
    MERGED_ALT_WT['AllDates_count_ALT_AND_WT'] = MERGED_ALT_WT['AllDates_count_ALT'] + MERGED_ALT_WT['AllDates_count_WT']  
    MERGED_ALT_WT['AllDates_count'] = MERGED_ALT_WT['AllDates_count_ALT'] / MERGED_ALT_WT['AllDates_count_ALT_AND_WT']
    MERGED_ALT_WT['Country'] = MERGED_ALT_WT['Country'].apply(lambda x: Modify_country_names(x))  
    ALL_COLORS = sns.color_palette("inferno", 10) #husl
    
    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    MERGED_ALT_WT.rename(columns = {'Country': 'name'}, inplace = True)
    #MERGED_ALT_WT_REDUCED = MERGED_ALT_WT[['COLOR', 'name']].copy()
    MERGED_ALT_WT_REDUCED = MERGED_ALT_WT[['AllDates_count', 'name']].copy()
    FullFrame_AllDates_withCodes = pd.merge(world, MERGED_ALT_WT_REDUCED, how = 'left', on = ['name'])
    FullFrame_AllDates_withCodes['AllDates_count'].fillna(0, inplace = True)
    FullFrame_AllDates_withCodes['COLOR'] = FullFrame_AllDates_withCodes['AllDates_count'].apply(lambda x: SET_COLOR(x, ALL_COLORS))
    COLORS_APPENDIX = FullFrame_AllDates_withCodes[['name', 'AllDates_count', 'COLOR']].copy()
    COLORS_APPENDIX.to_csv(f'{sys.argv[4]}/{sys.argv[3]}_{PEPNAME}_{Subset}.csv')
    
    world.to_csv('../WORLD_CODES_GEOPANDAS.csv')
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot()
    UNIQUE_COLORS = FullFrame_AllDates_withCodes['COLOR'].unique()
    for color in UNIQUE_COLORS:
        COUNTRIES = FullFrame_AllDates_withCodes[ FullFrame_AllDates_withCodes['COLOR'] == color ]
        
        if 'world' not in sys.argv[6].lower():
            print(sys.argv[6])
            Region = ''
            if 'north' in sys.argv[6].lower():
                Region = 'North America'
            elif 'south' in sys.argv[6].lower():
                Region = 'South America'
            else:
                Region = sys.argv[6]
            print(COUNTRIES)
            COUNTRIES = COUNTRIES[ COUNTRIES['continent'] == Region ]
        ######COUNTRIES = COUNTRIES[ COUNTRIES['continent'] == 'Europe' ]
        ######COUNTRIES = COUNTRIES[ COUNTRIES['name'] != 'Russia' ]
        COUNTRIES.drop(['COLOR'], inplace = True, axis=1)
        COUNTRIES.plot(ax=ax, color=color, edgecolor="black", alpha=0.5, aspect=1)
    ax.set_xticks([])
    ax.set_yticks([])

    plt.title("Basic Map of World with GeoPandas")
    plt.savefig(f'{sys.argv[4]}/{sys.argv[3]}_{PEPNAME}_{Subset}.png', dpi = 300)
    plt.clf()


def remove_empty(x):
    if x == '':
        return 0
    else:
        return x

def Prepare_data(filename, FinalFrame, TYPE):
    FRAME = pd.read_csv(filename, sep = '\t')  
    FRAME = FRAME.T.reset_index().T  
    FRAME[0] = FRAME[0].apply(lambda x: x.replace(' ', ',')[6:] )
    FRAME[['count', 'DATE']] = FRAME[0].str.split(',', expand=True)
    FRAME.rename(columns = {1: 'Continent', 2: 'Country'}, inplace = True)
    FRAME.drop([0], inplace = True, axis = 1)
    FRAME.dropna(inplace = True)    
    FRAME['DATE'] = FRAME['DATE'].apply(lambda x: Year_add_month(x))
    FRAME['DATE'] = FRAME['DATE'].apply(lambda x: Month_only(x))
    FRAME['DATE'] = FRAME['DATE'].apply(lambda x: Remove_00(x))
    FRAME.dropna(subset = ['DATE'], inplace = True)
    FRAME['DATE'] = FRAME['DATE'].apply(lambda x: Single_digit_month(x))
    FRAME['count'] = FRAME['count'].apply(lambda x: remove_empty(x) )
    FRAME['count'] = FRAME['count'].apply(lambda x: int(x) )
    grouped = FRAME.groupby(['DATE', 'Continent', 'Country'])['count'].sum().reset_index(name = filename) #'count'
    print(len(grouped.index.values))
    grouped[filename] = grouped[filename].apply(lambda x: int(x))
    grouped[filename].fillna(0, inplace = True)
    if TYPE == 'ALT':
        if FinalFrame.empty:
            FinalFrame = grouped
        else:
            FinalFrame = pd.merge(FinalFrame, grouped, how = 'outer', on = ['DATE', 'Continent', 'Country'])
            FinalFrame.fillna(0, inplace = True)
        print(f'FinalFrame: \n{FinalFrame}')

        return FinalFrame #grouped #FRAME

    elif (TYPE == 'WT') or (TYPE == 'SPECIF_PEP'):
        return grouped


def Peptide_Specific(Spec_alt_wt_pep):
    DATA_Files = []
    PEPnames = []
    FullFrames = [] #pd.DataFrame() ##Initiate dataframe
    ALL_OTHER_PEPS = pd.DataFrame()
    with open(sys.argv[1], 'r') as File:   #open and save to Results_Files dictionary HLA-specific processed netMHCpan files.
        for line in File:
            if len(line) > 1:
                print(sys.argv[3])
                templine = line[:-1]
                print(templine)
                if ('full' not in line[:-1]) and ('nextrainnuc' not in line[:-1]): #sys.argv[3]  '_' (sys.argv[3] not in line[:-1]) and 
                    #DATA_Files.append( Prepare_data(line[:-1], FullFrame) )
                    if (sys.argv[3] == Spec_alt_wt_pep[1]):
                        if (Spec_alt_wt_pep[0] in line[:-1]):
                            FullFrames.append( Prepare_data(line[:-1], 'standin', 'SPECIF_PEP') ) #FullFrame
                            PEP = line[:-1]
                            print(f'PEPNAME is {PEP}')
                            PEPnames.append(PEP)
                        elif sys.argv[3] not in line[:-1]:
                            ALL_OTHER_PEPS = Prepare_data(line[:-1], ALL_OTHER_PEPS, 'ALT')
                        else:
                            WT_frame = Prepare_data(line[:-1], 'standin', 'WT') #FullFrame 
                
    for FullFrame, PEPNAME in list(zip(FullFrames, PEPnames)):
        if not FullFrame.empty:
            FullFrame.sort_values(by = ['DATE', 'Continent', 'Country'], inplace = True, ascending = True)
            print(f'This is fullframe: \n{FullFrame}')
            print(FullFrame.dtypes)
            Modified_FullFrame = FullFrame.set_index(['DATE', 'Continent', 'Country'])
            Modified_FullFrame['Sum_across_peptides'] = Modified_FullFrame.sum(axis = 1)
            Modified_FullFrame.reset_index(inplace = True)
            Modified_FullFrame = Modified_FullFrame[['DATE', 'Continent', 'Country', 'Sum_across_peptides']].copy()
            Modified_ALL_OTHER_PEPS = ALL_OTHER_PEPS.set_index(['DATE', 'Continent', 'Country'])
            Modified_ALL_OTHER_PEPS['Sum_across_peptides_ALTPEPs'] = Modified_ALL_OTHER_PEPS.sum(axis = 1)
            Modified_ALL_OTHER_PEPS.reset_index(inplace = True)
            Modified_ALL_OTHER_PEPS = Modified_ALL_OTHER_PEPS[['DATE', 'Continent', 'Country', 'Sum_across_peptides_ALTPEPs']].copy()
            WT_frame.sort_values(by = ['DATE', 'Continent', 'Country'], inplace = True, ascending = True)
            print(f'This is WT_frame: \n{WT_frame}')
            print(WT_frame.dtypes)

            Modified_WT_frame = pd.merge(WT_frame, Modified_ALL_OTHER_PEPS, how= 'outer', on = ['DATE', 'Continent', 'Country'] )
            Modified_WT_frame = Modified_WT_frame.set_index(['DATE', 'Continent', 'Country'])
            Modified_WT_frame['Sum_across_peptides'] = Modified_WT_frame.sum(axis = 1)
            Modified_WT_frame.reset_index(inplace = True)
            Modified_WT_frame = Modified_WT_frame[['DATE', 'Continent', 'Country', 'Sum_across_peptides']].copy()

            Modified_FullFrame_2020 = Modified_FullFrame[Modified_FullFrame['DATE'].str.contains('2020')]
            Modified_WT_frame_2020 = Modified_WT_frame[Modified_WT_frame['DATE'].str.contains('2020')]
            Create_MAP_WholePANDEMIC(Modified_FullFrame_2020, Modified_WT_frame_2020, 'YEAR1', PEPNAME)

            Modified_FullFrame_2020_1 = Modified_FullFrame[Modified_FullFrame['DATE'].str.contains('2021')]
            Modified_WT_frame_2020_1 = Modified_WT_frame[Modified_WT_frame['DATE'].str.contains('2021')]
            Create_MAP_WholePANDEMIC(Modified_FullFrame_2020_1, Modified_WT_frame_2020_1, 'YEAR2', PEPNAME)

            Modified_FullFrame_2020_1_2 = Modified_FullFrame[Modified_FullFrame['DATE'].str.contains('2022')]
            Modified_WT_frame_2020_1_2 = Modified_WT_frame[Modified_WT_frame['DATE'].str.contains('2022')]
            Create_MAP_WholePANDEMIC(Modified_FullFrame_2020_1_2, Modified_WT_frame_2020_1_2, 'YEAR3', PEPNAME)
            Create_MAP_WholePANDEMIC(Modified_FullFrame, Modified_WT_frame, 'ALLYEARS', PEPNAME)


def NOT_Peptide_Specific():

    DATA_Files = []
    PEPnames = []
    FullFrame = pd.DataFrame() ##Initiate dataframe [] #
    ALL_OTHER_PEPS = pd.DataFrame()
    with open(sys.argv[1], 'r') as File:   #open and save to Results_Files dictionary HLA-specific processed netMHCpan files.
        for line in File:
            if len(line) > 1:
                print(sys.argv[3])
                templine = line[:-1]
                print(templine)
                if ('full' not in line[:-1]) and ('nextrainnuc' not in line[:-1]): #sys.argv[3]  '_' (sys.argv[3] not in line[:-1]) and 
                    #DATA_Files.append( Prepare_data(line[:-1], FullFrame) )
                    '''
                    if (sys.argv[3] == 'APRITFGGP') or (sys.argv[3] == 'HTTDPSFLGR'):
                        if ('ALRITFGGP' in line[:-1]) or ('HTTDLSFLGR' in line[:-1]) or ('HTTDSSFLGR' in line[:-1]):
                            FullFrames.append( Prepare_data(line[:-1], 'standin', 'SPECIF_PEP') ) #FullFrame
                            PEP = line[:-1]
                            print(f'PEPNAME is {PEP}')
                            PEPnames.append(PEP)
                    '''
                    if sys.argv[3] not in line[:-1]:
                        FullFrame = Prepare_data(line[:-1], FullFrame, 'ALT') #ALL_OTHER_PEPS = Prepare_data(line[:-1], ALL_OTHER_PEPS, 'ALT')
                    else:
                        WT_frame = Prepare_data(line[:-1], 'standin', 'WT') #FullFrame 
                
    #for FullFrame, PEPNAME in list(zip(FullFrames, PEPnames)):
    PEPNAME = ''
    if not FullFrame.empty:
        FullFrame.sort_values(by = ['DATE', 'Continent', 'Country'], inplace = True, ascending = True)
        print(f'This is fullframe: \n{FullFrame}')
        print(FullFrame.dtypes)
        Modified_FullFrame = FullFrame.set_index(['DATE', 'Continent', 'Country'])
        Modified_FullFrame['Sum_across_peptides'] = Modified_FullFrame.sum(axis = 1)
        Modified_FullFrame.reset_index(inplace = True)
        Modified_FullFrame = Modified_FullFrame[['DATE', 'Continent', 'Country', 'Sum_across_peptides']].copy()
        #Modified_ALL_OTHER_PEPS = ALL_OTHER_PEPS.set_index(['DATE', 'Continent', 'Country'])
        #Modified_ALL_OTHER_PEPS['Sum_across_peptides_ALTPEPs'] = Modified_ALL_OTHER_PEPS.sum(axis = 1)
        #Modified_ALL_OTHER_PEPS.reset_index(inplace = True)
        #Modified_ALL_OTHER_PEPS = Modified_ALL_OTHER_PEPS[['DATE', 'Continent', 'Country', 'Sum_across_peptides_ALTPEPs']].copy()
        Modified_WT_frame = WT_frame.sort_values(by = ['DATE', 'Continent', 'Country'], ascending = True) #inplace = True, 
        print(f'This is WT_frame: \n{WT_frame}')
        print(WT_frame.dtypes)

        #Modified_WT_frame = pd.merge(WT_frame, Modified_ALL_OTHER_PEPS, how= 'outer', on = ['DATE', 'Continent', 'Country'] )
        Modified_WT_frame = Modified_WT_frame.set_index(['DATE', 'Continent', 'Country'])
        Modified_WT_frame['Sum_across_peptides'] = Modified_WT_frame.sum(axis = 1)
        Modified_WT_frame.reset_index(inplace = True)
        Modified_WT_frame = Modified_WT_frame[['DATE', 'Continent', 'Country', 'Sum_across_peptides']].copy()

        Modified_FullFrame_2020 = Modified_FullFrame[Modified_FullFrame['DATE'].str.contains('2020')]
        Modified_WT_frame_2020 = Modified_WT_frame[Modified_WT_frame['DATE'].str.contains('2020')]
        Create_MAP_WholePANDEMIC(Modified_FullFrame_2020, Modified_WT_frame_2020, 'YEAR1', PEPNAME)

        Modified_FullFrame_2020_1 = Modified_FullFrame[Modified_FullFrame['DATE'].str.contains('2021')]
        Modified_WT_frame_2020_1 = Modified_WT_frame[Modified_WT_frame['DATE'].str.contains('2021')]
        Create_MAP_WholePANDEMIC(Modified_FullFrame_2020_1, Modified_WT_frame_2020_1, 'YEAR2', PEPNAME)

        Modified_FullFrame_2020_1_2 = Modified_FullFrame[Modified_FullFrame['DATE'].str.contains('2022')]
        Modified_WT_frame_2020_1_2 = Modified_WT_frame[Modified_WT_frame['DATE'].str.contains('2022')]
        Create_MAP_WholePANDEMIC(Modified_FullFrame_2020_1_2, Modified_WT_frame_2020_1_2, 'YEAR3', PEPNAME)

        Create_MAP_WholePANDEMIC(Modified_FullFrame, Modified_WT_frame, 'ALLYEARS', PEPNAME)


#if sys.argv[5] == 'yes': #'Peptide_Specific'
    #Peptide_Specific()
if sys.argv[5] == 'no': #'NOT_Peptide_Specific'
    NOT_Peptide_Specific()
else:
    Peptide_Specific(sys.argv[5].split(','))
    #print('INVALID INPUT')






