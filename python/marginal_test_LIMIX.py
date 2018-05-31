# coding: utf-8

module load matplotlib/1.3.1-goolf-1.4.10-Python-2.7.3
module load OpenSSL/1.0.1e-goolf-1.4.10
module load bzip2/1.0.6-goolf-1.4.10
module load zlib/1.2.7-goolf-1.4.10
module load Python/2.7.3-goolf-1.4.10

export PYTHONPATH=$PYTHONPATH:~/src
module load limix/0.6.5-goolf-1.4.10-Python-2.7.3

module load h5py/2.0.1-ictce-5.3.0-Python-2.7.3
module load HDF5/1.8.7-ictce-5.3.0

python 
import limix
import numpy as np
import scipy as SP
import scipy.stats as sps
import pylab as PL
from matplotlib import cm
import h5py
import pdb

# import LIMIX
import sys
import limix.modules.varianceDecomposition as var
import limix.modules.qtl as qtl
import limix.io.data as data
import limix.io.genotype_reader as gr
import limix.io.phenotype_reader as phr
import limix.io.data_util as data_util
import limix.utils.preprocess as preprocess
import limix.stats.fdr as FDR

# plotting and visualization utilities
from limix.utils.plot import *
# genotype summary stats
from limix.stats.geno_summary import *
import os
import cPickle
import sys
import pandas as pd

## Output directory
Out='./Limix/results'

## Import phenotype data
#SSH: f = open('/Limix/samples/RNA_wo_Index.csv')
#SSH: f.readline()
#SSH: RNA = []
#SSH: for l in f:
#SSH:    RNA.append(l.strip().split(','))
#SSH: RNA = np.array(RNA)
#SSH: RNA = RNA.astype(float)
#SSH: RNA = RNA.T
#SSH: RNA = (RNA-RNA.mean(axis=0))/RNA.std(axis=0)

f = h5py.File('/gale/netapp/home/shhuang/data/limix_tutorials/data/BYxRM/BYxRM.hdf5')
RNA = f['phenotype']['matrix'][:]
RNA = (RNA-RNA.mean(axis=0))/RNA.std(axis=0)

## Import genotype data
#SSH: f = h5py.File('/Limix/samples/SNPs.h5py')
#SSH: transsnp = f['snp'][:]
#transsnp.shape
transsnp = f['genotype']['matrix'][:].T

# trans kinship matrix
ts = transsnp.T
sumts = ts.sum(axis=0)
ts = ts[:,(sumts!=0)&(sumts!=ts.shape[0])]	# not fixed
ts = ts.astype(float)
ts = (ts-ts.mean(axis=0))/ts.std(axis=0)
transk = np.dot(ts,ts.T)
#transk

## Scaling Kinship matrix (from the Bjarni's scale_k())
c = sp.sum((sp.eye(len(transk)) - (1.0 / len(transk)) * sp.ones(transk.shape)) * sp.array(transk))
scalar = (len(transk) - 1) / c
transK = scalar * transk

RNA = RNA[:500,:]
ts = ts[:500,:]
transK = transK[:500,:500]

for i in range(0, RNA.shape[1]):
    
    y=RNA[:,i]
    # Mixed model
    lmm=qtl.test_lmm(snps=ts, pheno=y,K=transK)
    pvalues_lmm=pd.DataFrame(data=lmm.getPv().T, index=range(0,ts.shape[1]), columns=['lmm'])
    
    # Linear regression model
    lm=qtl.test_lmm(snps=ts, pheno=y)
    pvalues_lm=pd.DataFrame(data=lm.getPv().T, index=range(0,ts.shape[1]), columns=['lm'])
    
    # Export
    pval=pd.concat([pvalues_lmm,pvalues_lm], axis=1)
    np.savetxt(Out+'/'+str(i)+'.csv', pval, delimiter=',', fmt='%.3e')



