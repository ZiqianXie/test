###############################################################################
# This file contains the code to train the SpliceAI model.
###############################################################################

import numpy as np
import sys
import time
import h5py
import keras.backend as kb
import tensorflow as tf
from spliceai import *
from utils import *
from constants import *
from custom_utils import Callback, Chromosome


assert int(sys.argv[1]) in [80, 400, 2000, 10000]
###############################################################################
# Model
###############################################################################

L = 32
N_GPUS = 1

if int(sys.argv[1]) == 80:
    W = [11, 11, 11, 11]
    AR = [1, 1, 1, 1]
    BATCH_SIZE = 18*N_GPUS
elif int(sys.argv[1]) == 400:
    W = [11, 11, 11, 11, 11, 11, 11, 11]
    AR = [1, 1, 1, 1, 4, 4, 4, 4]
    BATCH_SIZE = 18*N_GPUS
elif int(sys.argv[1]) == 2000:
    W = [11, 11, 11, 11, 11, 11, 11, 11,
         21, 21, 21, 21]
    AR = [1, 1, 1, 1, 4, 4, 4, 4,
          10, 10, 10, 10]
    BATCH_SIZE = 12*N_GPUS
elif int(sys.argv[1]) == 10000:
    W = [11, 11, 11, 11, 11, 11, 11, 11,
         21, 21, 21, 21, 41, 41, 41, 41]
    AR = [1, 1, 1, 1, 4, 4, 4, 4,
          10, 10, 10, 10, 25, 25, 25, 25]
    BATCH_SIZE = 6*N_GPUS
# Hyper-parameters:
# L: Number of convolution kernels
# W: Convolution window size in each residual unit
# AR: Atrous rate in each residual unit

CL = 2 * int(np.sum(np.asarray(AR)*(np.asarray(W)-1)))
assert CL <= CL_max and CL == int(sys.argv[1])
print ("\033[1mContext nucleotides: %d\033[0m" % (CL))
print ("\033[1mSequence length (output): %d\033[0m" % (SL))
callback = Callback()


###############################################################################
# Training and validation
###############################################################################
try:
    h5f = h5py.File(data_dir + 'train_all_0.h5', 'r')
except OSError:
    h5f = h5py.File(data_dir + 'train_all_1.h5', 'r')

num_idx = len(h5f.keys())//2
h5f.close()
idx_all = np.random.permutation(num_idx)
idx_train = idx_all[:int(0.9*num_idx)]
idx_valid = idx_all[int(0.9*num_idx):]

EPOCH_NUM = 10

start_time = time.time()
fid = 0

for model_idx in range(1, 6):
    model_m = SpliceAI(L, W, AR)
    model_m.compile(loss=categorical_crossentropy_2d, optimizer='adam')
    for epoch_num in range(EPOCH_NUM):
        print("model {} epoch {}".format(model_idx, epoch_num))
        try:
            h5f = h5py.File(data_dir + 'train_all_{}.h5'.format(fid), 'r')
            fid = 1 - fid
        except OSError:
            h5f = h5py.File(data_dir + 'train_all_{}.h5'.format(1-fid), 'r')
        print("using {} h5 file for training".format(1 - fid))
        callback()
        for j in range(len(idx_train)):
            print('epoch progress: {:.2f}'.format(j/len(idx_train)))
            idx = np.random.choice(idx_train)
            X = h5f['X' + str(idx)][:]
            Y = h5f['Y' + str(idx)][:]
        
            Xc, Yc = clip_datapoints(X, Y, CL, N_GPUS) 
            model_m.fit(Xc, Yc, batch_size=BATCH_SIZE, verbose=0)
    
    
        
            # Printing metrics (see utils.py for details)
    
        print ("--------------------------------------------------------------")
        print ("\n\033[1mValidation set metrics:\033[0m")
    
        Y_true_1 = [[] for t in range(1)]
        Y_true_2 = [[] for t in range(1)]
        Y_pred_1 = [[] for t in range(1)]
        Y_pred_2 = [[] for t in range(1)]
    
        for idx in idx_valid:
    
            X = h5f['X' + str(idx)][:]
            Y = h5f['Y' + str(idx)][:]
    
            Xc, Yc = clip_datapoints(X, Y, CL, N_GPUS)
            Yp = model_m.predict(Xc, batch_size=BATCH_SIZE)
    
            if not isinstance(Yp, list):
                Yp = [Yp]
    
            for t in range(1):
    
                is_expr = (Yc[t].sum(axis=(1,2)) >= 1)
    
                Y_true_1[t].extend(Yc[t][is_expr, :, 1].flatten())
                Y_true_2[t].extend(Yc[t][is_expr, :, 2].flatten())
                Y_pred_1[t].extend(Yp[t][is_expr, :, 1].flatten())
                Y_pred_2[t].extend(Yp[t][is_expr, :, 2].flatten())
    
        print ("\n\033[1mAcceptor:\033[0m")
        for t in range(1):
            print_topl_statistics(np.asarray(Y_true_1[t]),
                                  np.asarray(Y_pred_1[t]))
    
        print ("\n\033[1mDonor:\033[0m")
        for t in range(1):
            print_topl_statistics(np.asarray(Y_true_2[t]),
                                  np.asarray(Y_pred_2[t]))
    
        print ("\n\033[1mTraining set metrics:\033[0m")
    
        Y_true_1 = [[] for t in range(1)]
        Y_true_2 = [[] for t in range(1)]
        Y_pred_1 = [[] for t in range(1)]
        Y_pred_2 = [[] for t in range(1)]
    
        for idx in idx_train[:len(idx_valid)]:
    
            X = h5f['X' + str(idx)][:]
            Y = h5f['Y' + str(idx)][:]
    
            Xc, Yc = clip_datapoints(X, Y, CL, N_GPUS)
            Yp = model_m.predict(Xc, batch_size=BATCH_SIZE)
    
            if not isinstance(Yp, list):
                Yp = [Yp]
    
            for t in range(1):
    
                is_expr = (Yc[t].sum(axis=(1,2)) >= 1)
    
                Y_true_1[t].extend(Yc[t][is_expr, :, 1].flatten())
                Y_true_2[t].extend(Yc[t][is_expr, :, 2].flatten())
                Y_pred_1[t].extend(Yp[t][is_expr, :, 1].flatten())
                Y_pred_2[t].extend(Yp[t][is_expr, :, 2].flatten())
    
        print ("\n\033[1mAcceptor:\033[0m")
        for t in range(1):
            print_topl_statistics(np.asarray(Y_true_1[t]),
                                  np.asarray(Y_pred_1[t]))
    
        print ("\n\033[1mDonor:\033[0m")
        for t in range(1):
            print_topl_statistics(np.asarray(Y_true_2[t]),
                                  np.asarray(Y_pred_2[t]))
    
        print ("Learning rate: %.5f" % (kb.get_value(model_m.optimizer.lr)))
        print ("--- %s seconds ---" % (time.time() - start_time))
        start_time = time.time()
    
        print ("--------------------------------------------------------------")
    
        model_m.save('./Models/SpliceAI' + sys.argv[1]
                   + '_c{}'.format(model_idx) + '.h5')
    
        if epoch_num > 4:
            kb.set_value(model_m.optimizer.lr,
                         0.5*kb.get_value(model_m.optimizer.lr))
                # Learning rate decay
        h5f.close()
        
###############################################################################

