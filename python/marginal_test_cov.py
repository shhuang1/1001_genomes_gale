# coding: utf-8
import sys
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
import limix
import multiprocessing
import numpy as np
import scipy as SP
import scipy.stats as sps
from pygwas.core import genotype
import pylab as PL
from matplotlib import cm
import h5py
import pdb

# import LIMIX
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

def run_lmm(ts,RNA,COV,i,transK,out_csv):
    y=RNA[:,i]
    # Mixed model
    lmm=qtl.test_lmm(snps=ts, pheno=y,K=transK, covs=COV)
    pvalues_lmm=pd.DataFrame(data=lmm.getPv().T, index=range(0,ts.shape[1]), columns=['lmm'])
    
    # Linear regression model
    lm=qtl.test_lmm(snps=ts, pheno=y, covs=COV)
    pvalues_lm=pd.DataFrame(data=lm.getPv().T, index=range(0,ts.shape[1]), columns=['lm'])
    
    # Export
    pval=pd.concat([pvalues_lmm,pvalues_lm], axis=1)
    #np.savetxt(Out+'/'+str(i)+'.csv', pval, delimiter=',', fmt='%.6e')
    np.savetxt(out_csv, pval, delimiter='\t', fmt='%.6e')

## Output directory
out_dir = '/gale/netapp/home/shhuang/projects/1001_genomes/marginal_test_cov_01'

## Import genotype data
#SSH: f = h5py.File('/Limix/samples/SNPs.h5py')
#SSH: transsnp = f['snp'][:]
#transsnp.shape
#SNP_file = '/gale/netapp/home/shhuang/data/1001_genomes/gmi_release_v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX/imputed_snps_binary.hdf5'
SNP_file = '/gale/netapp/home/shhuang/data/1001_genomes/1001_250k_fullimputed/PYGWAS_GENOTYPES/1/all_chromosomes_binary.hdf5'
f = h5py.File(SNP_file)
SNP_accx = ['X'+acc for acc in f['accessions']]

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

RNA_file = '/gale/netapp/home/shhuang/projects/1001_genomes/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01.txt'
RNA_df = pd.read_csv(RNA_file,sep='\t',header=0,index_col=0) # genes x accessions
RNA_accx = list(RNA_df.columns.values)
RNA_genes = list(RNA_df.index)

## get common accessions in the same order for genotype, phenotype and phenotype covariates
accx_common = [accx for accx in SNP_accx if accx in RNA_accx]
RNA = RNA_df.as_matrix(columns=accx_common).T # accession x genes
match = lambda a, b: [ b.index(x) if x in b else None for x in a ]
transsnp = f['snps'][:,match(accx_common,SNP_accx)] # SNP x accessions

# trans kinship matrix
#SSH: ts = transsnp.T # accessions x SNP
#SSH: sumts = ts.sum(axis=0)
#SSH: pos_tf = (sumts!=0)&(sumts!=ts.shape[0])
#SSH: ts = ts[:,pos_tf]	# not fixed
#SSH: ts = ts.astype(float)
#SSH: ts = (ts-ts.mean(axis=0))/ts.std(axis=0)
#SSH: transk = np.dot(ts,ts.T) # accessions x accessions

#transk
## Scaling Kinship matrix (from the Bjarni's scale_k())
#SSH: c = sp.sum((sp.eye(len(transk)) - (1.0 / len(transk)) * sp.ones(transk.shape)) * sp.array(transk))
#SSH: scalar = (len(transk) - 1) / c
#SSH: transK = scalar * transk


## save the filtered genotypes
np.savetxt(os.path.join(out_dir,'positions_tf.txt'), pos_tf, delimiter='\t',fmt='%d')

for k in range(3,5):

    ## Import phenotype covariates
    COV_file = '/gale/netapp/home/shhuang/projects/1001_genomes/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_W_k%d.txt'%(k)
    COV_df = pd.read_csv(COV_file,sep='\t',header=0,index_col=0) # cov x accessions
    COV = COV_df.ix[accx_common].as_matrix()
    print(COV_file)

    #x = Parallel(n_jobs=num_cores)(delayed(runn_lmm)(ts,RNA,i,transK,out_csv) for i in range(0,10))
    #x = Parallel(n_jobs=num_cores,verbose=100,max_nbytes=1e6)(delayed(has_shareable_memory)(run_lmm(ts,RNA,COV,i,transK,out_csv)) for i in range(0,RNA.shape[1]))
    out_dir_k = os.path.join(out_dir,'gNorm_k%d'%k)
    if not os.path.exists(out_dir_k): os.makedirs(out_dir_k)
    x = Parallel(n_jobs=num_cores,verbose=100)(
        delayed(run_lmm)(ts,RNA,COV,i,transK,os.path.join(out_dir_k,'%s.csv'%RNA_genes[i]))
        for i in range(0,10))

