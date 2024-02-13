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

def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq)]: #get frame and nuc sequence from 5'-> 3' strand only (specific to SARS)
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table)) #translate using Biopython tool translate()
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start) # '*' is STOP CODON. We are constantly looking for the next stop codon.
                if aa_end == -1: # If cant find next stop codon, we've reached end of genome.
                    aa_end = trans_len
                temp = trans.find("M", aa_start) #Once have stop codon, look for corresponding 1st start codon.
                if temp != -1:
                    aa_start = temp
                if aa_end - aa_start >= min_protein_length: # Determine start and end of ORF in genetic code.
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    answer.append((start, end, strand, trans[aa_start:aa_end])) #Store ORF info in a list of ORFs
                aa_start = aa_end + 1
    answer.sort()
    return answer #Return list of ORFs.


def Create_Deletions(pos, del_seq, sequence):
            
    First = sequence[:pos] #+1 from .base file (Del always one more); -1 because python starts from 0, but not Seq data; Selection omits the given position
    Second = sequence[pos + len(del_seq) + 1 :] #Pos + length of deletion, but don't want to include the last nuc from deletion, so + 1
    New_seq = First + Second

    #### Select Spike protein for translation only #####
    New_seq = New_seq[21535 : 25385]

    table = 1 #Translation table. 1 = standard translation.
    min_pro_len = 37 

    orf_list = find_orfs_with_trans(New_seq, table, min_pro_len) #Submit genome fasta, trans table and min protein length to find_orfs_with_trans function to find all ORFs.
    for start, end, strand, pro in orf_list:

        if int(start) == 0: #find Spike  21535
            #print(f'{start}-{end}\n{strand}\n{pro}')
            print(pro)
            print(len(del_seq))
            return pro


def Initiate_thread(Positions, Deletions, Sample, THREAD): #conn

    with open(f'Refsequence_Nuc.fasta', "rU") as File: #Open file containing reference genome.   
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
                        ## (Pos + 1)%3 == 0 ): Check that deletion is at last position of any codon
                        ## ( len(deletion_seq)%3 == 0 ): Check that deletion is the length of a codon
                        ## ( Pos +  len(deletion_seq) >= 23623): Check that deletion includes stop codon right before peptide
                        ## ( Pos +  len(deletion_seq) > 23693 ): Check that deletion does't "eat" the peptide
                        Peptide = "Present"
                        FRAME = "Frame+1"

                    #elif ( (Pos + 2)%3 == 0 ) and ( (len(deletion_seq) + 2)%3 == 0 ) and ( Pos +  len(deletion_seq) >= 23623) and ( Pos +  len(deletion_seq) <= 23693 ):
                        #Peptide = "Present"

                    #elif ( (Pos + 1)%3 == 0 ) and ( (len(deletion_seq) + 2)%3 == 0 ) and ( Pos +  len(deletion_seq) >= 23623) and ( Pos +  len(deletion_seq) <= 23693 ):
                        #Peptide = "Present"

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

                    '''
                    PROTEIN = Create_Deletions(Pos, deletion_seq, sequence)

                    if "LPYPQILLL" in PROTEIN:
                        Peptide = "Present"
                    else:
                        Peptide = "Absent"
                    '''
                    Results.append( [frequency, len(deletion_seq), Pos, Peptide, Samp, deletion_seq, FRAME] )
                TRACKER += 1
                ####print(f'POGRESS: {TRACKER} of {Total_Length} for thread {THREAD}')

                '''
                for Del in Dels:
                    Freq = Del[0]
                    DelSeq = Del[ Del.find(':') + 1 : ]
                    PROTEIN = Create_Deletions(Pos, DelSeq, sequence)

                    if "LPYPQILLL" in PROTEIN:
                        Peptide = "Present"
                    else:
                        Peptide = "Absent"
                    
                    Results.append( [Freq, len(DelSeq), Pos, Peptide, Samp] )
                '''
            print(f'finished thread {THREAD}')
    
    return Results
    #conn.send(Results)
    #conn.close()
    #return Results

            



if __name__ == '__main__':

    Data = pd.read_csv(f'{sys.argv[1]}', index_col=0)
    print(Data)

    Chunks = np.array_split(Data, int(sys.argv[2]))
    print(Chunks)


    #Set Multiprocessing variables
    parent_num = 0
    Conns = {}
    procs = []
    Processed_frames = []

    thread_tracker = 1

    for Chunk in Chunks:
        #print(type(Chunk))
        #print(Chunk)
        #print(Chunk.columns.values)

        Positions, Deletions, Sample = Chunk['Pos'].values, Chunk['Del'].values, Chunk['Sample'].values

        TEMP = Initiate_thread(Positions, Deletions, Sample, thread_tracker) #Conns[parent_num][1], 
        if len(TEMP) > 0:
            Processed_frames = Processed_frames + TEMP


    print('All Done!! Merging files')

    AllResults_Frame = pd.DataFrame(data = Processed_frames, columns = ['Del_Freq', 'Del_Length', 'Positions', 'Peptide', 'Sample', 'Deletion', 'FRAME'])
    
    print('FILES MERGED!!')

    AllResults_Frame.to_csv('NOTRANSLATE_FrameShift_results_First_Three_waves_Combined_Strands_DelAbove3.csv')


    '''
    All_Results = Iterate_through_Samples(Positions, Deletions, Sample)

    AllResults_Frame = pd.DataFrame(data = All_Results, columns = ['Del_Freq', 'Del_Length', 'Positions', 'Peptide', 'Sample', 'Deletion'])

    AllResults_Frame.to_csv('FrameShift_results_First_Three_waves_Combined_Strands_DelAbove1.csv')
    '''
    #AllResults_Frame = AllResults_Frame.T
