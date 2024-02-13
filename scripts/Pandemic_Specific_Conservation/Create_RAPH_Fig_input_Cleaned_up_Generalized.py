import os
import os.path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import seaborn as sns
import math

def Create_Figure(NEWDATA, PEP_LIST):
    datanull=pd.read_csv(f"{sys.argv[4]}Sliding_window_peptides_nuc_stats.csv",sep=" ",header=None)
    datamhc=pd.read_csv(f"{sys.argv[4]}MHCValidator_peptides_nucposFromRaf_stats.csv",sep=" ",header=None)
    #covac1=["FTALTQHGK","LLLLDRLNQL","LLLDRLNQL","TRFQTLLAL","TLLALHRSY","FYVYSRVKNL","RVKNLNSSR"]
    #datamhc_covac1 = datamhc[datamhc[0].isin(covac1)]


    ###########NEWDATA=pd.read_csv("../New_Peptide_list.csv")
    ###########NEWDATA.drop(['Unnamed: 0'], inplace = True, axis = 1)
    NEWDATA.columns = [0,1]
    #####VAX_2=['KTIQPRVEK', 'TTDPSFLGRYM', 'PTDNYITTY', 'YLFDESGEFKL', 'LPKEITVAT', 'TTDPSFLGRY']
    #####NEWDATA_PEPS=['KLPDDFTGC', 'TLNDLNETL', 'NAPRITFGGP', 'VPYNMRVI', 'RANNTKGSL', 'GPMVLRGLIT', 'STTTNIVTR', 'TGSNVFQTR', 'HTTDPSFLGR', 'RTIKVFTTV']
    #######NEWDATA_VAX2 = NEWDATA[NEWDATA[0].isin(VAX_2)]

    ########NEWDATA = NEWDATA[NEWDATA[0].isin(NEWDATA_PEPS)]
    titles=["summing all stat","averaging all stat","second codon","lastcodon"]

    fig = plt.figure(figsize=(10, 4))
    print(datamhc)
    #print(datamhc_covac1)

    print(NEWDATA)
    ##########print(NEWDATA_VAX2)


    for i in [1]: #,2,3,4
        ax = plt.subplot(1, 1, i) #4
        d=np.log10(datanull[i].replace(0, min(datanull[i].replace(0,100))))
        ax.scatter(d, np.random.normal(1,0.02,len(d)), alpha=0.7, color = 'white', edgecolors = 'black')
        ###m=np.log10(datamhc[i].replace(0, min(datamhc[i].replace(0,100))))
        m=np.log10(NEWDATA[i].replace(0, min(NEWDATA[i].replace(0,100))))
        ax.scatter(m, np.random.normal(1,0.02,len(m)), alpha=0.7, color = 'orange', edgecolors = 'black', s = 100) #np.ones(len(m))
        ###c=np.log10(datamhc_covac1[i].replace(0, min(datamhc_covac1[i].replace(0,100))))
        ###########c=np.log10(NEWDATA_VAX2[i].replace(0, min(NEWDATA_VAX2[i].replace(0,100))))
        ############ax.scatter(c, np.random.normal(1,0.02,len(c)), alpha=0.7, color = 'green', edgecolors = 'black', s = 100) #np.ones(len(c))
        ax.set_xlim([-5, 1]) #7
        ax.set_ylim([0.5, 1.5])
        ax.set_xticks(np.arange(-5, 1)) #7
        ax.set_xticklabels(10.0**np.arange(-5, 1)) #7
        ax.set_title('log10 of '+titles[i-1])
        maxleft=1
        maxright=1
        ytoplot=-5.5
        
        
        '''
        for j in np.argsort(datamhc[i]):
            pep=datamhc.iloc[j]
            ytoplot=ytoplot+0.25
            c="orange"
            if pep[0] in covac1:
                c="green"
            ax.annotate(pep[0],
                        xy = (1, np.log10(pep[i])),
                        xytext = (0.5, ytoplot),
                        fontsize = 10, 
                        ha="right",
                        arrowprops = dict(facecolor = c,lw=0,width=.5,headwidth=.5),
                        color = c)
        '''

        ######plt.cm.get_cmap("hsv", j)
    #plt.show()
    plt.savefig(f'{sys.argv[2]}RAPH_FIG.pdf', dpi = 300)




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
PEPTIDE_list = []
###########peptides_to_keep = ['KLPDDFTGC', 'TLNDLNETL', 'NAPRITFGGP', 'VPYNMRVI', 'RANNTKGSL', 'GPMVLRGLIT', 'STTTNIVTR', 'TGSNVFQTR', 'HTTDPSFLGR', 'KTIQPRVEK', 'TTDPSFLGRYM', 'PTDNYITTY', 'YLFDESGEFKL', 'LPKEITVAT', 'TTDPSFLGRY', 'RTIKVFTTV']
with open(sys.argv[1], 'r') as File:   #open and save to Results_Files dictionary HLA-specific processed netMHCpan files.
    for line in File:
        if len(line) > 1:
            ###print(sys.argv[3])
            templine = line[:-1]
            print(templine)
            
            PEPTIDE = templine[: templine.find('.')]
            print(PEPTIDE)
            PEPTIDE_list.append(PEPTIDE)
            #############if PEPTIDE in peptides_to_keep:
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

Create_Figure(RATIOS_FRAME, PEPTIDE_list)

RATIOS_FRAME.to_csv(f'{sys.argv[2]}New_Peptide_list.csv')

Sorted_RATIOS_FRAME = RATIOS_FRAME.sort_values(by = [1], ascending = False)

print(Sorted_RATIOS_FRAME)
 






