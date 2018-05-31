# coding: utf-8
# load sample_relatedness from file (None means no sample_releatedness is applied)
# load covariances from file (None means no covariance is applied)
# command line option for normalization (None means no normalization is applied)
import errno
import os
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
import cPickle
import sys
import pandas as pd

from shared import LoggerFactory

def strip_x(tg_ecotypeidx): return tg_ecotypeidx.replace("X","")
strip_xvec = np.vectorize(strip_x)

def add_x(tg_ecotypeid): return "X"+tg_ecotypeid
add_xvec = np.vectorize(add_x)

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    return path

def getPhenotypes(pheno_reader,phenotype_IDs=None,phenotype_query=None,sample_idx=None,center=True,intersection=False):
    """load Phenotypes
    
    Args:
    idx_start:      phenotype indices to load (start individual index)
    idx_end:       phenotype indices to load (end individual index)
    phenotype_IDs:  names of phenotypes to load
    impute:         imputation of missing values (default: True)
    intersection:   restrict observation to those obseved in all phenotypes (true) or at least in one phenotype (false)? (default: False)

    Returns:
    phenotypes:     phenotype values
    sample_idx_intersect:        index of individuals in phenotypes after filtering missing values
    """
    if phenotype_IDs is not None:
        I = SP.array([SP.nonzero(pheno_reader.phenotype_ID==n)[0][0] for n in phenotype_IDs])
    elif phenotype_query is not None:
        try:
            I = pheno_reader.index_frame.query(phenotype_query).values[:,0]
            
            #if there are no results we won't actually get an exception, we just get an
            #empty response
            if len(I) == 0:
                print "query '%s' yielded no results!" % (phenotype_query)
                I = SP.zeros([0],dtype="int")
        except Exception, arg:
                
                print "query '%s' yielded no results: %s" % (phenotype_query, str(arg))
                
                I = SP.zeros([0],dtype="int")
    else:
        I = SP.arange(pheno_reader.phenotype_ID.shape[0])
        
    phenotypes = SP.array(pheno_reader.pheno_matrix[:,I],dtype='float')
    phenotypes = phenotypes[sample_idx]
    sample_ID = pheno_reader.sample_ID[sample_idx]
    Iok = (~SP.isnan(phenotypes))
    
    if intersection:
        sample_idx_intersect = Iok.all(axis=1)
    else:
        sample_idx_intersect = Iok.any(axis=1)
            
    phenotypes = phenotypes[sample_idx_intersect]
    Iok = Iok[sample_idx_intersect]
                    
    if center:
        for i in xrange(phenotypes.shape[1]):
            ym = phenotypes[Iok[:,i],i].mean()
            phenotypes[:,i] -= ym
            phenotypes[~Iok[:,i],i] = ym
            phenotypes[:,i] /= phenotypes[:,i].std()
    phenotypes = pd.DataFrame(data=phenotypes, index=sample_ID[sample_idx_intersect],columns=pheno_reader.phenotype_ID[I])
    #calculate overlap of missing values
    return phenotypes, sample_idx_intersect

def main():

    if 1:
        geno_file,pheno_file,norm_mode,K_file,cov_file,RNA_start,RNA_end,out_dir = sys.argv[1:]

    if 0:
      
        geno_file = '/gale/netapp/home/shhuang/data/1001_genomes/dmC_bins/dmC_filtered/dmC_filtered_methylation_4.hdf5'
        pheno_file = '/gale/netapp/home/shhuang/projects/1001_genomes/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k4_vst2_cv0p05_UQCounts_1001gT.hdf5'
        norm_mode = 'RIN'
        out_dir = 'test_v8'
        K_file = '/gale/netapp/home/shhuang/data/1001_genomes/gmi_release_v3.1/X1001tx_filter1/norm_cov_1001tx_filter1.csv'
        cov_file = '/gale/netapp/home/shhuang/projects/1001_genomes/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_W_k4.txt'
        RNA_start,RNA_end = 0,5
       
    make_sure_path_exists(out_dir)
    log_dir = make_sure_path_exists(os.path.join(out_dir,'logs'))
    logger = LoggerFactory.get_logger(os.path.join(log_dir,'%s-%s.log'%(RNA_start,RNA_end)),
                                      file_level=logging.DEBUG,console_level=logging.DEBUG)
    LoggerFactory.log_command(logger,sys.argv[1:])
    logger.info('Output directory: %s',out_dir)
    out_graphics_dir = make_sure_path_exists(os.path.join(out_dir,'graphics'))
    out_results_dir = make_sure_path_exists(os.path.join(out_dir,'results'))

    RNA_start,RNA_end = int(RNA_start),int(RNA_end)
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

    logger.info('Sample relatedness %s',K_file)
    logger.info('Loading sample relatedness from %s',K_file)
    if (K_file=='None'):
        sample_relatedness = None
    else:
        logger.info('Start loading covariance from %s',K_file)
        K_df = pd.read_csv(K_file,sep='\t',header=None,index_col=0) # accessions x accessions
        K_df.index = ['%d'%i for i in K_df.index]
        K_df.columns = K_df.index
        sample_relatedness = K_df.loc[dataset.sample_ID,dataset.sample_ID].as_matrix()
    sample_relatedness_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'sample_relatedness'))
    pl.imshow(sample_relatedness,aspect='auto')
    plt.savefig(os.path.join(sample_relatedness_dir,'sample_relatedness.png'))

    logger.info('Subset phenotype to index %d-%d',RNA_start,RNA_end)
    phenotype_ID = dataset.phenotype_ID[RNA_start:RNA_end]
    phenotypes,sample_idx = getPhenotypes(dataset.pheno_reader,phenotype_IDs=phenotype_ID,\
                                          sample_idx=dataset.sample_idx['pheno'])

    logger.info('Phenotype normalization: %s',norm_mode)
    if norm_mode=='None':
        phenotype_vals = phenotypes.values
    elif norm_mode=='RIN':
        phenotype_vals = preprocess.rankStandardizeNormal(phenotypes.values)
    elif norm_mode=='boxcox':
        phenotype_vals,maxlog = preprocess.boxcox(phenotypes.values)
    else:
        logger.info('Normalization mode %s is not recognized.  Use None',norm_mode)
        phenotype_vals = phenotypes.values
        
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
        
        plot_normal(phenotype_vals[:,ip],alpha=0.8,figure=fig)
        plt.title("%s" % p_ID)
        fig.savefig(out_file)
        plt.close(fig)

    logger.info('Sample covariance %s',cov_file)
    if (cov_file=='None'):
        cov = None
    else:
        logger.info('Start loading covariance from %s',cov_file)
        cov_df = pd.read_csv(cov_file,sep='\t',header=0,index_col=0) # cov x accessions
        cov = cov_df.ix[add_xvec(dataset.sample_ID)].as_matrix()

    #logger.info('Start testing: LM')
    #lm = qtl.test_lm(snps=snps[sample_idx].astype('int'),pheno=phenotype_vals,
    #                 covs=cov,verbose=True)
    #convert P-values to a DataFrame for nice output writing:
    #pvalues_lm = pd.DataFrame(data=lm.pvalues.T,index=dataset.geno_ID,
    #                          columns=phenotype_ID)
    logger.info('Start testing: LMM')
    lmm = qtl.test_lmm(snps=snps[sample_idx].astype('int'),pheno=phenotype_vals,
                       K=sample_relatedness,covs=cov,verbose=True)
    pvalues_lmm = pd.DataFrame(data=lmm.pvalues.T,index=dataset.geno_ID,
                               columns=phenotype_ID)

    #lm_pval_dir = make_sure_path_exists(os.path.join(out_results_dir,'lm_pval'))
    lmm_pval_dir = make_sure_path_exists(os.path.join(out_results_dir,'lmm_pval'))
    logger.info('Saving P-values to text file in %s',lmm_pval_dir)
    for ip,p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        #pvalues_lm[p_ID].to_csv(os.path.join(lm_pval_dir,'%s.txt'%p_ID),
        #                        header=True,index=False)
        pvalues_lmm[p_ID].to_csv(os.path.join(lmm_pval_dir,'%s.txt'%p_ID),
                                 header=True,index=False)

    # Genome-wide manhatton plots for one phenotype:
    manh_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'manhattan'))
    logger.info('Plotting Manhattan plots in %s',manh_dir)
    for ip, p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        out_file = os.path.join(manh_dir,'%s.png'%p_ID)
        fig = plt.figure(figsize=[12,8])
        #subpl = plt.subplot(2,1,1)
        #plot_manhattan(posCum=position['pos_cum'],pv=pvalues_lm[p_ID].values,chromBounds=chromBounds,thr_plotting=0.05)
        #plt.title('%s, LM'%p_ID)
        #subpl = plt.subplot(2,1,2)
        plot_manhattan(posCum=position['pos_cum'],pv=pvalues_lmm[p_ID].values,chromBounds=chromBounds,thr_plotting=0.05)
        plt.title('%s, LMM'%p_ID)
        fig.savefig(out_file)
        plt.close(fig)
        
    # SNP vs. phenotype
    snp_pheno_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'snp_pheno'))
    logger.info('Plotting phenotype vs. SNP to %s',snp_pheno_dir)
    for ip, p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        out_file = os.path.join(snp_pheno_dir,'%s.png'%p_ID)
        fig = plt.figure(figsize=[3,3])#create the figure

        #find maximum squared beta value
        pheno_vals, s_idx = getPhenotypes(dataset.pheno_reader,phenotype_IDs=[p_ID],\
                                          sample_idx=dataset.sample_idx['pheno'])
        imax = lmm.pvalues[ip].argmin()
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
    pval_hist_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'pval_hist'))
    logger.info('Plotting P-value histograms to %s',pval_hist_dir)
    for ip, p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        out_file = os.path.join(pval_hist_dir,'%s.png'%p_ID)
        fig = plt.figure(figsize=[7,3])
        
        #subpl = plt.subplot(1,2,1)
        #plt.hist(pvalues_lm[p_ID].values,20,normed=True)
        #plt.plot([0,1],[1,1],"r")
        #plt.title("%s, LM" % p_ID)
        #plt.xlabel("P-value")
        #plt.ylabel("Frequency")

        #subpl = plt.subplot(1,2,2)
        plt.hist(pvalues_lmm[p_ID].values,20,normed=True)
        plt.plot([0,1],[1,1],"r")
        plt.title("%s, LMM" % p_ID)
        plt.xlabel("P-value")
        plt.ylabel("Frequency")
        fig.savefig(out_file)
        plt.close(fig)

   # Quantile-Quantile plots
    qqplot_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'qqplot'))
    logger.info('Plotting Q-Q plots to %s',qqplot_dir)
    for ip, p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
        out_file = os.path.join(qqplot_dir,'%s.png'%p_ID)
        fig = plt.figure(figsize=[7,3])
       
        #subpl = plt.subplot(1,2,1)
        #qqplot(pvalues_lm[p_ID].values)
        #plt.title("%s, LM" % p_ID)
        #subpl = plt.subplot(1,2,2)
        qqplot(pvalues_lmm[p_ID].values)
        plt.title("%s, LMM" % p_ID)
       
        fig.savefig(out_file)
        plt.close(fig)
       
    # P value scatter plot
    # logger.info('Plotting LM vs LMM P-values')
    # pval_lmvslmm_dir = make_sure_path_exists(os.path.join(out_graphics_dir,'pval_lmvslmm'))
    # for ip,p_ID in enumerate(dataset.phenotype_ID[RNA_start:RNA_end]):
    #     out_file = os.path.join(pval_lmvslmm_dir,'%s.png'%p_ID)
    #     fig = plt.figure(figsize=[3,3])
    #     plt.plot(-sp.log10(pvalues_lm[p_ID]),-sp.log10(pvalues_lmm[p_ID]),'.')
    #     ymax = max(plt.xlim()[1],plt.ylim()[1])
    #     plt.plot([0,ymax],[0,ymax],'k--')
    #     plt.xlabel('LM')
    #     plt.ylabel('LMM')
    #     plt.title(p_ID)
    #     fig.savefig(out_file)
    #     plt.close(fig)

    logger.info('Done with all plots!')

    logger.info('Done!')
        
if __name__=='__main__':
    main()
