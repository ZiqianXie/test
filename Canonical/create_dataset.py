###############################################################################
'''This parser takes as input the .h5 file produced by create_datafile.py and
outputs a .h5 file with datapoints of the form (X, Y), which can be understood
by Keras models.'''
###############################################################################

import h5py
import numpy as np
import time
from utils import *
from constants import *


def create_dataset(*argv):
    start_time = time.time()
    
    assert argv[0] in ['train', 'test', 'all']
    assert argv[1] in ['0', '1', 'all']
    
    if argv[0] == 'train':
        CHROM_GROUP = ['chr11', 'chr13', 'chr15', 'chr17', 'chr19', 'chr21',
                       'chr2', 'chr4', 'chr6', 'chr8', 'chr10', 'chr12',
                       'chr14', 'chr16', 'chr18', 'chr20', 'chr22', 'chrX', 'chrY']
    elif argv[0] == 'test':
        CHROM_GROUP = ['chr1', 'chr3', 'chr5', 'chr7', 'chr9']
    else:
        CHROM_GROUP = ['chr1', 'chr3', 'chr5', 'chr7', 'chr9',
                       'chr11', 'chr13', 'chr15', 'chr17', 'chr19', 'chr21',
                       'chr2', 'chr4', 'chr6', 'chr8', 'chr10', 'chr12',
                       'chr14', 'chr16', 'chr18', 'chr20', 'chr22', 'chrX', 'chrY']
    
    ###############################################################################
    
    NAME = []      # Gene symbol
    PARALOG = []   # 0 if no paralogs exist, 1 otherwise
    CHROM = []     # Chromosome number
    STRAND = []    # Strand in which the gene lies (+ or -)
    TX_START = []  # Position where transcription starts
    TX_END = []    # Position where transcription ends
    JN_START = []  # Positions where canonical exons end
    JN_END = []    # Positions where canonical exons start
    SEQ = []       # Nucleotide sequence
    
    fpr2 = open(sequence, 'r')
    
    with open(splice_table, 'r') as fpr1:
        for line1 in fpr1:
    
            line2 = fpr2.readline()
    
            data1 = re.split('\n|\t', line1)[:-1]
            data2 = re.split('\n|\t|:|-', line2)[:-1]
    
            assert data1[2] == data2[0]
            assert int(data1[4]) == int(data2[1])+CL_max//2+1
            assert int(data1[5]) == int(data2[2])-CL_max//2
    
            if (data1[2] not in CHROM_GROUP):
                continue
    
            if (argv[1] != data1[1]) and (argv[1] != 'all'):
                continue
    
            NAME.append(data1[0])
            PARALOG.append(int(data1[1]))
            CHROM.append(data1[2])
            STRAND.append(data1[3])
            TX_START.append(data1[4])
            TX_END.append(data1[5])
            JN_START.append(data1[6::2])
            JN_END.append(data1[7::2])
            SEQ.append(data2[3])
    
    fpr2.close()
    
    try:
        h5f2 = h5py.File(data_dir+argv[0]+'_'+argv[1]+'_0.h5', 'w')
    except OSError:
        h5f2 = h5py.File(data_dir+argv[0]+'_'+argv[1]+'_1.h5', 'w')
    CHUNK_SIZE = 100
    
    for i in range(len(SEQ)//CHUNK_SIZE):
        # Each dataset has CHUNK_SIZE genes
        print("create_h5_progress: {:.2f}".format(i/(len(SEQ)//CHUNK_SIZE)))
        
        if (i+1) == len(SEQ)//CHUNK_SIZE:
            NEW_CHUNK_SIZE = CHUNK_SIZE + len(SEQ)%CHUNK_SIZE
        else:
            NEW_CHUNK_SIZE = CHUNK_SIZE
    
        X_batch = []
        Y_batch = [[] for t in range(1)]
    
        for j in range(NEW_CHUNK_SIZE):
    
            idx = i*CHUNK_SIZE + j
    
            X, Y = create_datapoints(SEQ[idx], STRAND[idx],
                                     TX_START[idx], TX_END[idx],
                                     JN_START[idx], JN_END[idx])
    
            X_batch.extend(X)
            for t in range(1):
                Y_batch[t].extend(Y[t])
    
        X_batch = np.asarray(X_batch).astype('int8')
        for t in range(1):
            Y_batch[t] = np.asarray(Y_batch[t]).astype('int8')
    
        h5f2.create_dataset('X' + str(i), data=X_batch)
        h5f2.create_dataset('Y' + str(i), data=Y_batch)
    
    h5f2.close()
    
    print ("--- %s seconds ---" % (time.time() - start_time))

###############################################################################         
