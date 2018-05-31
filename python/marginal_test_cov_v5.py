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
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import pylab as PL
import matplotlib.pyplot as plt
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
import errno
import cPickle
import sys
import pandas as pd

from shared import LoggerFactory

def strip_x(tg_ecotypeidx): return tg_ecotypeidx.replace("X","")
strip_xvec = np.vectorize(strip_x)

def add_x(tg_ecotypeid): return "X"+tg_ecotypeid
add_xvec = np.vectorize(add_x)

import os
import errno

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    return path            

def main():

    geno_file,pheno_file,cov_file,RNA_start,RNA_end,out_dir = sys.argv[1:]
    make_sure_path_exists(out_dir)
    log_dir = make_sure_path_exists(os.path.join(out_dir,'logs'))
    logger = LoggerFactory.get_logger(os.path.join(log_dir,'%s-%s.log'%(RNA_start,RNA_end)),
                                      file_level=logging.DEBUG,console_level=logging.DEBUG)
    LoggerFactory.log_command(logger,sys.argv[1:])
    logger.info('Output directory: %s',out_dir)
    
    #geno_file = '/gale/netapp/home/shhuang/data/1001_genomes/gmi_release_v3.1/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1_2.hdf5' 
    #pheno_file = '/gale/netapp/home/shhuang/projects/1001_genomes/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_normCounts_k4_1001g_vst2_cv0p05_rinT.hdf5'
    #out_dir = '.'
    #cov_file = '/gale/netapp/home/shhuang/projects/1001_genomes/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_W_k4.txt'
    #RNA_start,RNA_end = 0,5
    RNA_start,RNA_end = int(RNA_start),int(RNA_end)
    out_graphics_dir = make_sure_path_exists(os.path.join(out_dir,'graphics'))
    out_results_dir = make_sure_path_exists(os.path.join(out_dir,'results'))
 
    logger.info('Loading genotype from %s',geno_file)
    geno_reader = gr.genotype_reader_tables(geno_file)
    logger.info('Loading phenotype from %s',pheno_file)
    pheno_reader = phr.pheno_reader_tables(pheno_file)
    pheno_reader.sample_ID = strip_xvec(pheno_reader.sample_ID)

    # the data object allows to query specific genotype or phenotype data
    logger.info('Creating QTL dataset')
    dataset = data.QTLData(geno_reader=geno_reader,pheno_reader=pheno_reader)
    # getting genotypes
    snps = dataset.getGenotypes() #SNPS
    position = dataset.getPos()
    position,chromBounds = data_util.estCumPos(position=position,offset=100000)

    logger.info('Calculating sample relatedness')
    # non-normalized and normalized sample relatedeness matrix
    sample_relatedness_unnormalized = dataset.getCovariance(normalize=False)
    sample_relatedness  = sample_relatedness_unnormalized/sample_relatedness_unnormalized.diagonal().mean()
    sample_relatedness_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'sample_relatedness'))
    pl.imshow(sample_relatedness,aspect='auto')
    plt.savefig(os.path.join(sample_relatedness_dir,'sample_relatedness_norm.png'))

    logger.info('Subset phenotype to index %d-%d',RNA_start,RNA_end)
    phenotype_ID = dataset.phenotype_ID[RNA_start:RNA_end]
    phenotype_vals,sample_idx = dataset.getPhenotypes(phenotype_ID)

    N = snps.shape[0] #number of individuals
    S = snps.shape[1] #number of SNPs
    P = phenotype_vals.shape[1]#number of phenotypes
    logger.info('Number of individuals: %d; number of SNPs: %d; number of phenotypes: %d',
                N,S,P)

    logger.info('Plotting phenotype histograms')
    phenohist_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'phenohist'))
    for ip, p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        out_file = os.path.join(phenohist_dir,'%s.png'%p_ID)
        fig = plt.figure(figsize=[3,3])#create the figure
        
        plot_normal(phenotype_vals.values[:,ip],alpha=0.8,figure=fig)
        plt.title("%s" % p_ID)
        fig.savefig(out_file)
        plt.close(fig)

    logger.info('Start loading covariance from %s',cov_file)
    cov_df = pd.read_csv(cov_file,sep='\t',header=0,index_col=0) # cov x accessions
    cov = cov_df.ix[add_xvec(dataset.sample_ID)].as_matrix()
    logger.info('Finished')

    logger.info('Start testing: LM')
    lm = qtl.test_lm(snps=snps[sample_idx].astype('int'),pheno=phenotype_vals.values,
                     covs=cov,verbose=True)
    #convert P-values to a DataFrame for nice output writing:
    pvalues_lm = pd.DataFrame(data=lm.pvalues.T,index=dataset.geno_ID,
                              columns=phenotype_ID)
    logger.info('Start testing: LMM')
    lmm = qtl.test_lmm(snps=snps[sample_idx].astype('int'),pheno=phenotype_vals.values,
                       K=sample_relatedness,covs=cov,verbose=True)
    pvalues_lmm = pd.DataFrame(data=lmm.pvalues.T,index=dataset.geno_ID,
                               columns=phenotype_ID)

    logger.info('Saving P-values to text file')
    lm_pval_dir = make_sure_path_exists(os.path.join(out_results_dir,'lm_pval'))
    lmm_pval_dir = make_sure_path_exists(os.path.join(out_results_dir,'lmm_pval'))
    for ip,p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        pvalues_lm[p_ID].to_csv(os.path.join(lm_pval_dir,'%s.txt'%p_ID),
                                header=True,index=False)
        pvalues_lmm[p_ID].to_csv(os.path.join(lmm_pval_dir,'%s.txt'%p_ID),
                                 header=True,index=False)

    # Genome-wide manhatton plots for one phenotype:
    logger.info('Plotting Manhattan plots')
    manh_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'manhattan'))
    for ip, p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        out_file = os.path.join(manh_dir,'%s.png'%p_ID)
        fig = plt.figure(figsize=[12,8])
        subpl = plt.subplot(2,1,1)
        plot_manhattan(posCum=position['pos_cum'],pv=pvalues_lm[p_ID].values,chromBounds=chromBounds,thr_plotting=0.05)
        plt.title('%s, LM'%p_ID)
        subpl = plt.subplot(2,1,2)
        plot_manhattan(posCum=position['pos_cum'],pv=pvalues_lmm[p_ID].values,chromBounds=chromBounds,thr_plotting=0.05)
        plt.title('%s, LMM'%p_ID)
        fig.savefig(out_file)
        plt.close(fig)
        
    # SNP vs. phenotype
    logger.info('Plotting phenotype vs. SNP')
    snp_pheno_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'snp_pheno'))
    for ip, p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        out_file = os.path.join(snp_pheno_dir,'%s.png'%p_ID)
        fig = plt.figure(figsize=[3,3])#create the figure

        #find maximum squared beta value
        pheno_vals, s_idx = dataset.getPhenotypes([p_ID])
        imax = lm.pvalues[ip].argmin()
        i_0 = snps[s_idx,imax]==0
        #plot SNP vs. phenotype for max beta
        plt.plot(snps[s_idx,imax]+0.05*np.random.randn(snps[s_idx,imax].shape[0]),pheno_vals.values,'.',alpha=0.5)
        plt.xlabel("SNP")
        plt.ylabel("phenotype")
        plt.xlim([-0.5,2.5])
        plt.title("%s" % p_ID)
        fig.savefig(out_file)
        plt.close(fig)

    # P-value histgrams
    logger.info('Plotting P-value histograms')
    pval_hist_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'pval_hist'))
    for ip, p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        out_file = os.path.join(pval_hist_dir,'%s.png'%p_ID)
        fig = plt.figure(figsize=[7,3])
        
        subpl = plt.subplot(1,2,1)
        plt.hist(pvalues_lm[p_ID].values,20,normed=True)
        plt.plot([0,1],[1,1],"r")
        plt.title("%s, LM" % p_ID)
        plt.xlabel("P-value")
        plt.ylabel("Frequency")

        subpl = plt.subplot(1,2,2)
        plt.hist(pvalues_lmm[p_ID].values,20,normed=True)
        plt.plot([0,1],[1,1],"r")
        plt.title("%s, LMM" % p_ID)
        plt.xlabel("P-value")
        plt.ylabel("Frequency")
        fig.savefig(out_file)
        plt.close(fig)

   # Quantile-Quantile plots
    logger.info('Plotting Q-Q plots')
    qqplot_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'qqplot'))
    for ip, p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        out_file = os.path.join(qqplot_dir,'%s.png'%p_ID)
        fig = plt.figure(figsize=[7,3])
       
        subpl = plt.subplot(1,2,1)
        qqplot(pvalues_lm[p_ID].values)
        plt.title("%s, LM" % p_ID)
        subpl = plt.subplot(1,2,2)
        qqplot(pvalues_lmm[p_ID].values)
        plt.title("%s, LMM" % p_ID)
       
        fig.savefig(out_file)
        plt.close(fig)
       
    # P value scatter plot
    logger.info('Plotting LM vs LMM P-values')
    pval_lmvslmm_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'pval_lmvslmm'))
    for ip,p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        out_file = os.path.join(pval_lmvslmm_dir,'%s.png'%p_ID)
        fig = plt.figure(figsize=[3,3])
        plt.plot(-sp.log10(pvalues_lm[p_ID]),-sp.log10(pvalues_lmm[p_ID]),'.')
        ymax = max(plt.xlim()[1],plt.ylim()[1])
        plt.plot([0,ymax],[0,ymax],'k--')
        plt.xlabel('LM')
        plt.ylabel('LMM')
        plt.title(p_ID)
        fig.savefig(out_file)
        plt.close(fig)

    logger.info('Done with all plots!')

    logger.info('Done!')
        
if __name__=='__main__':
    main()
