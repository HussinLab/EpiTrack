import os
import os.path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import seaborn as sns
import math
from multiprocessing import Process, Pipe
from Bio import SeqIO #have to install Biopython before using this.


def Initiate_thread(Positions, Deletions, Sample, THREAD): #conn

    with open(f'{sys.argv[3]}/Refsequence_Nuc.fasta', "rU") as File: #Open file containing reference genome.   
        for record in SeqIO.parse(File, "fasta"): #Parse through all genomes in the fasta file.
            RefProteins_nuc = {}
            sequence = record.seq
    
            Results = []    
            TRACKER = 0
            Total_Length = len(Positions)
            print(f'started thread {THREAD}')
            for Pos, Dels, Samp in list( zip( Positions, Deletions, Sample ) ):
                Pos = int(Pos)

                Dels = Dels.split('|')
                
                ###Pairs = []
                FREQS = []
                DELSEQS = []
                for Del in Dels:
                    Freq = int( Del[: Del.find(':')] )
                    DelSeq = Del[ Del.find(':') + 1 : ]
                    DelSeq = DelSeq.upper()
                    ##Pairs.append( [Freq, DelSeq] )
                    if DelSeq in DELSEQS:
                        FREQS_Loc = DELSEQS.index(DelSeq)
                        FREQS[ FREQS_Loc ] = FREQS[ FREQS_Loc ] + Freq
                    else:
                        FREQS.append(Freq)
                        DELSEQS.append(DelSeq)

                #CombStrands_frame = pd.DataFrame(Pairs, columns = ['Freq', 'DelSeq'])
                #Grouped = CombStrands_frame.groupby(['DelSeq'])['Freq'].sum().reset_index()
                
                #print(f'{Samp}\n{Grouped}')

                #for frequency, deletion_seq in list( zip( Grouped['Freq'].values, Grouped['DelSeq'].values ) ):
                for frequency, deletion_seq in list( zip( FREQS, DELSEQS ) ):
                    
                    if ( (len(deletion_seq) + 2)%3 == 0 ) and ( Pos +  len(deletion_seq) >= 23623) and ( Pos +  len(deletion_seq) <= 23693 ): #( (Pos)%3 == 0 ) and 
                        
                        Peptide = "Present"
                        FRAME = "Frame+1"

                    elif ( (len(deletion_seq) + 1)%3 == 0 ) and ( Pos +  len(deletion_seq) >= 23623) and ( Pos +  len(deletion_seq) <= 23693 ): 
                        Peptide = "Absent"
                        FRAME = "Frame+2"

                    elif ( (len(deletion_seq))%3 == 0 ) and ( Pos +  len(deletion_seq) >= 23623) and ( Pos +  len(deletion_seq) <= 23693 ): 
                        Peptide = "Absent"
                        FRAME = "Frame+0"

                    #elif ( (Pos + 1)%3 != 0 ) and ( len(deletion_seq)%3 != 0 ) and ( Pos +  len(deletion_seq) >= 23623) and ( Pos +  len(deletion_seq) <= 23693 ):
                        #Peptide = "Absent"
                    else:
                        Peptide = "OUT_OF_RANGE"
                        FRAME = "OUT_OF_RANGE"

                    Results.append( [frequency, len(deletion_seq), Pos, Peptide, Samp, deletion_seq, FRAME] )
                TRACKER += 1
                ####print(f'POGRESS: {TRACKER} of {Total_Length} for thread {THREAD}')

            print(f'finished thread {THREAD}')
    
    return Results
    
if __name__ == '__main__':

    Data = pd.read_csv(f'{sys.argv[1]}', index_col=0)
    print(Data)

    Chunks = np.array_split(Data, int(sys.argv[4]))
    print(Chunks)


    #Set Multiprocessing variables
    parent_num = 0
    Conns = {}
    procs = []
    Processed_frames = []

    thread_tracker = 1

    for Chunk in Chunks:
       

        Positions, Deletions, Sample = Chunk['Pos'].values, Chunk['Del'].values, Chunk['Sample'].values

        TEMP = Initiate_thread(Positions, Deletions, Sample, thread_tracker) #Conns[parent_num][1], 
        if len(TEMP) > 0:
            Processed_frames = Processed_frames + TEMP


    print('All Done!! Merging files')

    AllResults_Frame = pd.DataFrame(data = Processed_frames, columns = ['Del_Freq', 'Del_Length', 'Positions', 'Peptide', 'Sample', 'Deletion', 'FRAME'])
    
    print('FILES MERGED!!')

    AllResults_Frame.to_csv(f'{sys.argv[2]}/NOTRANSLATE_FrameShift_results_First_Three_waves_Combined_Strands_DelAbove3.csv')

