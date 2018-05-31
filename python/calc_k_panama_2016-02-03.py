# coding: utf-8
import os
import h5py
import sys
import csv
import limix
import logging
import numpy as np
import scipy as sp
import scipy.stats as sps
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import pylab as pl
import matplotlib.pyplot as plt

# import LIMIX
import limix.modules.varianceDecomposition as var
import limix.modules.panama as panama
import limix.modules.qtl as qtl
import limix.io.data as data
import limix.io.genotype_reader as gr
import limix.io.phenotype_reader as phr
import limix.io.data_util as data_util
import limix.utils.preprocess as preprocess
import limix.utils.util_functions as util_functions
import limix.stats.fdr as FDR

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

def draw_and_save_panama(p,graphics_prefix,results_prefix):
    # retrieve stuff
    Kpanama = p.get_Kpanama()
    Ktot    = p.get_Ktot()
    Kpop    = p.Kpop
    vComp   = p.get_varianceComps()
    Xpanama = p.get_Xpanama()

    out_file = graphics_prefix + '_Kmat.png'
    # compare Kpanama matrix, Kpop, and the total matrix
    fig = plt.figure(figsize=[10,8])
    subplt = pl.subplot(2,2,1)
    pl.title('Kpanama')
    pl.imshow(Kpanama,vmin=0,vmax=1,interpolation='none',cmap=cm.afmhot)
    pl.colorbar(ticks=[0,0.5,1],orientation='horizontal')
    subplt.set_xticks([])
    subplt.set_yticks([])
    subplt = pl.subplot(2,2,2)
    pl.title('Kpop')
    pl.imshow(Kpop,vmin=0,vmax=1,interpolation='none',
              cmap=cm.afmhot)
    pl.colorbar(ticks=[0,0.5,1],orientation='horizontal')
    subplt.set_xticks([])
    subplt.set_yticks([])
    subplt = pl.subplot(2,2,3)
    print(list(vComp.keys()))
    vComp_a = sp.array([vComp[k] for k in ['Kpanama','Kpop','noise']])
    
    pl.bar(sp.arange(3)+0.2,vComp_a,width=0.6)
    subplt.set_xticks([0.5,1.5,2.5])
    subplt.set_xticklabels(['Kpanama','Kpop','noise'])
    pl.ylabel('Variance Explained')
    subplt = pl.subplot(2,2,4)
    pl.title('Ktot')
    pl.imshow(Ktot,vmin=0,vmax=1,interpolation='none',cmap=cm.afmhot)
    pl.colorbar(ticks=[0,0.5,1],orientation='horizontal')
    subplt.set_xticks([])
    subplt.set_yticks([])
    fig.savefig(out_file)
    plt.close(fig)

    out_dict = {'Ktot':Ktot,'Kpanama':Kpanama,'vComp':vComp,'Xpanama':Xpanama}
    out_file = results_prefix + '_dat.hdf5'
    o = h5py.File(out_file,'w')
    util_functions.smartDumpDictHdf5(out_dict,o)
    o.close()


def main():

    geno_file = '/gale/netapp/home/shhuang/data/1001_genomes/gmi_release_v3.1/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5' 
    pheno_file = '/gale/netapp/home/shhuang/projects/1001_genomes/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_1001g_vst_cv0p05T.hdf5'
    out_dir = '/gale/netapp/home/shhuang/projects/1001_genomes/calc_k_panama_2016-02-03'
    out_graphics_dir = make_sure_path_exists(os.path.join(out_dir,'graphics'))
    graphics_prefix = os.path.join(out_graphics_dir,'calc_k_panama_2016-02-03-')
    results_prefix = os.path.join(out_dir,'calc_k_panama_2016-02-03-')
    logger = LoggerFactory.get_logger(os.path.join(out_dir,'calc_k_panama_2016-02-03.log'),
                                      file_level=logging.DEBUG,console_level=logging.DEBUG)
 
    logger.info('Loading genotype from %s',geno_file)
    geno_reader = gr.genotype_reader_tables(geno_file)
    logger.info('Loading phenotype from %s',pheno_file)
    pheno_reader = phr.pheno_reader_tables(pheno_file)
    pheno_reader.sample_ID = strip_xvec(pheno_reader.sample_ID)

    # the data object allows to query specific genotype or phenotype data
    logger.info('Creating QTL dataset')
    dataset = data.QTLData(geno_reader=geno_reader,pheno_reader=pheno_reader)
    # import data
    phenotypes,sample_idx = dataset.getPhenotypes(intersection=True)
    sample_relatedness = dataset.getCovariance()
    
    # determine the number of ranks to consider in the PANAMA matrix
    # by looking at the variance explained by PCs
    cum_var = panama.PC_varExplained(phenotypes.values)
    out_file = graphics_prefix + 'cum_var.png'
    fig = plt.figure(figsize=[5,4])
    subplt = pl.subplot(1,1,1)
    pl.bar(sp.arange(50)+0.5,cum_var[:50],width=1)
    pl.xlim(0,50)
    ticks = sp.linspace(0,50,11); ticks[0] = 1
    subplt.set_xticks(ticks)
    fig.savefig(out_file)
    plt.close(fig)

    for r in [10,15,20]:        
        p = panama.PANAMA(Y=phenotypes.values,Kpop=sample_relatedness)
        logger.info('Training r=%d',r)
        p.train(rank=r)
        draw_and_save_panama(p,graphics_prefix+'_K%d'%r,results_prefix+'_K%d'%r)
    
        
if __name__=='__main__':
    main()
