# coding: utf-8
import sys
import csv
import itertools as it
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
import limix
import logging
import multiprocessing
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

def filter_maf_snps(geno,maf_lb,maf_ub):
    logger = logging.getLogger()
    num_snps = geno.num_snps
    snps_ix = []
    num_accessions = len(geno.accessions)
    for i,snps in enumerate(geno.get_snps_iterator()):
        l = SP.bincount(snps)
        maf = min(l)/float(num_accessions)
        if maf<=maf_lb or maf>maf_ub:
            snps_ix.append(i)
    numRemoved = len(snps_ix)
    geno.filter_snps_ix(snps_ix)
    logger.info("Removed %d SNPs of MAF <=%f or MAF>%f, leaving %d SNPs in total." % (numRemoved, maf_lb,maf_ub,geno.num_snps))
    return (num_snps,numRemoved)
    
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

def save_kinship_in_text_format(filename, k, accessions):
    with open(filename, 'w') as f:
        writer = csv.writer(f,delimiter='\t')
        for acc, row in it.izip(accessions, k):
            writer.writerow([acc]+row.tolist())

 def main():

    RNA_start,RNA_end = int(sys.argv[1]),int(sys.argv[2])

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)

    step_list = [""]


    ## Output directory
    out_dir = '/gale/netapp/home/shhuang/projects/1001_genomes/marginal_test_cov_03'
    logger.info('Out dir %s',out_dir)

    ## Import genotype data
    #SSH: f = h5py.File('/Limix/samples/SNPs.h5py')
    #SSH: transsnp = f['snp'][:]
    #transsnp.shape
    SNP_file = '/gale/netapp/home/shhuang/data/1001_genomes/gmi_release_v3.1/SNP_matrix_imputed_hdf5/1001_SNP_MATRIX/imputed_snps_binary.hdf5'
    #SNP_file = '/gale/netapp/home/shhuang/data/1001_genomes/1001_250k_fullimputed/PYGWAS_GENOTYPES/1/all_chromosomes_binary.hdf5'
    geno = genotype.load_hdf5_genotype_data(SNP_file)
    #f = h5py.File(SNP_file)
    SNP_accx = ['X'+acc for acc in geno.accessions]
    logger.info('Finished reading SNP from %s',SNP_file)

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
    logger.info('Finished reading RNA from %s',RNA_file)

    ## get common accessions in the same order for genotype, phenotype and phenotype covariates
    accx_common = [accx for accx in SNP_accx if accx in RNA_accx]
    RNA = RNA_df.as_matrix(columns=accx_common).T # accession x genes

    ## filtering
    match = lambda a, b: [ b.index(x) if x in b else None for x in a ]
    geno.filter_accessions_ix(match(accx_common,SNP_accx))
    #(num_snps,num_removed) = geno.filter_non_binary()
    (num_snps,num_removed) = filter_maf_snps(geno,0.01,0.5)

    #snps_ix = np.ix_(geno.filter_snps,geno.filter_accessions)
    #snps = geno.snps[:][snps_ix]
    #snps = snps.T.astype(int)

    snps = np.vstack(geno.get_snps_iterator(is_chunked=True))
    snps = snps.T.astype(int)
    logger.info('Finished filtering SNPs')

    #SSH: transsnp = f['snps'][:,match(accx_common,SNP_accx)] # SNP x accessions
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
    #np.savetxt(os.path.join(out_dir,'positions_tf.txt'), geno.filter_snps, delimiter='\t',fmt='%d')

    kinship_file = '/gale/netapp/home/shhuang/projects/1001_genomes/marginal_test_cov_01/kinship_maf1.hdf5'
    kinship_csv = '/gale/netapp/home/shhuang/projects/1001_genomes/marginal_test_cov_01/kinship_maf1.csv'
    kinship_scaled_csv = '/gale/netapp/home/shhuang/projects/1001_genomes/marginal_test_cov_01/kinship_maf1_scaled.csv'
    if ('calc_kinship' in step_list):
        K = geno.get_ibd_kinship_matrix()
        scaledK = kinship.scale_k(K).astype(float)
        kinship.save_kinship_to_file(kinship_file,K,geno.accessions,geno.num_snps)
        save_kinship_in_text_format(kinship_csv, K, geno.accessions)
        save_kinship_in_text_format(kinship_scaled_csv, scaledK, geno.accessions)
    else:
        load_k = kinship.load_kinship_from_file(kinship_file,scaled=False)
        K = load_k['k'].astype(float)
        scaledK = kinship.scale_k(K)
    logger.info('Done with kinship file %s',kinship_file)

    num_cores = 3
    logger.info('Start testing')
    for k in range(4,5):

        ## Import phenotype covariates
        COV_file = '/gale/netapp/home/shhuang/projects/1001_genomes/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_W_k%d.txt'%(k)
        COV_df = pd.read_csv(COV_file,sep='\t',header=0,index_col=0) # cov x accessions
        COV = COV_df.ix[accx_common].as_matrix()
        logger.info('Cov file %s',COV_file)

        #x = Parallel(n_jobs=num_cores)(delayed(runn_lmm)(ts,RNA,i,transK,out_csv) for i in range(0,10))
        #x = Parallel(n_jobs=num_cores,verbose=100,max_nbytes=1e6)(delayed(has_shareable_memory)(run_lmm(ts,RNA,COV,i,transK,out_csv)) for i in range(0,RNA.shape[1]))
        out_dir_k = os.path.join(out_dir,'gNorm_k%d'%k)
        logger.info('Output directory %s',out_dir_k)
        if not os.path.exists(out_dir_k): os.makedirs(out_dir_k)
        #x = Parallel(n_jobs=num_cores,verbose=100,max_nbytes=1e6)(
        #    delayed(run_lmm)(snps,RNA,COV,i,scaledK,os.path.join(out_dir_k,'%s.csv'%RNA_genes[i]))
        #    for i in range(RNA_start,RNA_end))
        logger.debug('RNA start %d end %d',RNA_start,RNA_end)
        run_lmm_chunk(snps,RNA,COV,RNA_start,RNA_end,list(RNA_genes[RNA_start:RNA_end]),
                      scaledK,os.path.join(out_dir_k,'%d_%d.csv'%(RNA_start,RNA_end)))

if __name__=='__main__':
    main()
