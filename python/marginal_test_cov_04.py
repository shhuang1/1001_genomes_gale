# coding: utf-8
import sys
import csv
import itertools as it
import limix
import logging
import numpy as np
import scipy as SP
import scipy.stats as sps
from pygwas.core import genotype
from pygwas.core import kinship
from pygwas.core import phenotype
import pylab as PL
from matplotlib import cm
import h5py
import pdb
from shared import LoggerFactory

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

def my_calc_ibd_kinship(genotype, dtype='single',chunk_size=None):
    logger = logging.getLogger()
    num_snps = genotype.num_snps
    num_lines = len(genotype.accessions)
    if chunk_size == None:
        chunk_size = num_lines
    k_mat = sp.zeros((num_lines, num_lines), dtype=dtype)
    chunk_i = 0
    snps = genotype.get_snps_iterator(is_chunked=True,chunk_size=chunk_size)
    for snps_chunk in snps:
        print('%d %d %d'%(chunk_i,chunk_size,num_snps))
        chunk_i += 1
        snps_array = sp.array(snps_chunk)


        snps_array = snps_array.T
        norm_snps_array = (snps_array - sp.mean(snps_array, 0)) / sp.std(snps_array, 0)
        x = sp.mat(norm_snps_array.T)
        k_mat += x.T * x
        #logger.debug('%0.2f%%' % (100.0 * (min(1, ((chunk_i + 1.0) * chunk_size) / genotype.original_num_snps))))
        logger.debug('%f',chunk_i * chunk_size)
    logger.info('Processed %d SNPs',chunk_i*chunk_size)
    k_mat = k_mat / float(num_snps)
    return k_mat

kinship.calc_ibd_kinship = my_calc_ibd_kinship

def run_lmm_chunk(ts,RNA,COV,RNA_start,RNA_end,RNA_columns,transK,out_csv):
    logger = logging.getLogger()
    y=RNA[:,RNA_start:RNA_end]

    # Mixed model
    logger.debug('Start lmm')
    lmm=qtl.test_lmm(snps=ts, pheno=y,K=transK, covs=COV,verbose=True)
    logger.debug('Done lmm')
    pvalues_lmm=pd.DataFrame(data=lmm.getPv().T, index=range(0,ts.shape[1]), columns=RNA_columns)
    logger.debug('Export')
    #np.savetxt(out_csv, pvalues_lmm, delimiter='\t', fmt='%.6e')
    pvalues_lmm.to_csv(out_csv,sep='\t',header=True,index=False)

    # Linear regression model
    #logger.debug('Start lm')
    #lm=qtl.test_lmm(snps=ts, pheno=y, covs=COV,verbose=True)
    #logger.debug('Done lm')
    #pvalues_lm=pd.DataFrame(data=lm.getPv().T, index=range(0,ts.shape[1]), columns=['lm'])

def main():

    geno_hdf5,kinship_hdf5,RNA_csv,COV_file,RNA_start,RNA_end,out_file = sys.argv[1:]
    RNA_start,RNA_end = int(RNA_start),int(RNA_end)
    
    logger = LoggerFactory.get_logger(out_file+'.log',file_level=logging.DEBUG,console_level=logging.DEBUG)
    LoggerFactory.log_command(logger,sys.argv[1:])

    step_list = [""]

    ## Import genotype data
    logger.info('Start reading SNP from %s',geno_hdf5)
    geno = genotype.load_hdf5_genotype_data(geno_hdf5)
    SNP_accx = ['X'+acc for acc in geno.accessions]
    logger.info('Finished reading SNP from %s',geno_hdf5)

    ## Import phenotype data
    logger.info('Finished reading RNA from %s',RNA_csv)
    RNA_df = pd.read_csv(RNA_csv,sep='\t',header=0,index_col=0) # genes x accessions
    RNA_accx = list(RNA_df.columns.values)
    RNA_genes = list(RNA_df.index)
    logger.info('Finished reading RNA from %s',RNA_csv)

    ## get common accessions in the same order for genotype, phenotype and phenotype covariates
    logger.info('Consolidate accessions from genotype and RNA file')
    accx_common = [accx for accx in SNP_accx if accx in RNA_accx]
    RNA = RNA_df.as_matrix(columns=accx_common).T # accession x genes
    match = lambda a, b: [ b.index(x) if x in b else None for x in a ]
    geno.filter_accessions_ix(match(accx_common,SNP_accx))
    logger.info('Number of accessions: genotype file %d, RNA file %d, common %d',len(SNP_accx),len(RNA_accx),len(accx_common))

    logger.info('Start building SNP matrix in memory')
    snps = np.vstack(geno.get_snps_iterator(is_chunked=True))
    snps = snps.T.astype(int)
    logger.info('Finished')

    logger.info('Start loading kinship matrix from %s',kinship_hdf5)
    load_k = kinship.load_kinship_from_file(kinship_hdf5,scaled=False)
    K0 = load_k['k'].astype(float)
    K_accx = ['X'+acc for acc in load_k['accessions']]
    K_accx_ix = np.ix_(match(accx_common,K_accx),match(accx_common,K_accx))
    K = K0[K_accx_ix]
    logger.info('Finished')

    logger.info('Start loading covariance from %s',COV_file)
    COV_df = pd.read_csv(COV_file,sep='\t',header=0,index_col=0) # cov x accessions
    COV = COV_df.ix[accx_common].as_matrix()
    logger.info('Finished')

    logger.info('Start association testing: RNA start %d, RNA end %d',RNA_start,RNA_end)

    run_lmm_chunk(snps,RNA,COV,RNA_start,RNA_end,list(RNA_genes[RNA_start:RNA_end]),
                  K,out_file)
     

if __name__=='__main__':
    main()
