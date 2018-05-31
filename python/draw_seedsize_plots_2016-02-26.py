# coding: utf-8
import errno
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
import pylab as PL
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

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    return path            

def strip_x(tg_ecotypeidx): return tg_ecotypeidx.replace("X","")
strip_xvec = np.vectorize(strip_x)

def add_x(tg_ecotypeid): return "X"+tg_ecotypeid
add_xvec = np.vectorize(add_x)

def main():

    geno_file = '/gale/netapp/home/shhuang/data/1001_genomes/gmi_release_v3.1/1001genomes_snp-short-indel_only_ACGTN_filter1.hdf5'
    pheno_file = '/gale/netapp/home/shhuang/data/1001_genomes/seed_size/accx_size.hdf5'
    expr_file = '/gale/netapp/home/shhuang/projects/1001_genomes/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_normCounts_k4_vst2T.csv.hdf5'
    out_dir = '/gale/netapp/home/shhuang/projects/1001_genomes/draw_seedsize_plots_2016-02-26'
    out_graphics_dir = make_sure_path_exists(os.path.join(out_dir,'graphics'))
    graphics_prefix = os.path.join(out_graphics_dir,'draw_seedsize_plots_2016-02-26-')
    results_prefix = os.path.join(out_dir,'draw_seedsize_plots_2016-02-26-')
    logger = LoggerFactory.get_logger(os.path.join(out_dir,'draw_seedsize_plots_2016-02-26.log'),
                                      file_level=logging.DEBUG,console_level=logging.DEBUG)
 
    logger.info('Loading genotype from %s',geno_file)
    geno_reader = gr.genotype_reader_tables(geno_file)
    logger.info('Loading phenotype from %s',pheno_file)
    pheno_reader = phr.pheno_reader_tables(pheno_file)
    pheno_reader.sample_ID = strip_xvec(pheno_reader.sample_ID)
    logger.info('Loading expression from %s',expr_file)
    expr_reader = phr.pheno_reader_tables(expr_file)
    expr_reader.sample_ID = strip_xvec(expr_reader.sample_ID)

    # the data object allows to query specific genotype or phenotype data
    logger.info('Creating QTL dataset')
    dataset = data.QTLData(geno_reader=geno_reader,pheno_reader=pheno_reader)
    exprset = data.QTLData(geno_reader=geno_reader,pheno_reader=expr_reader)
    # import data
    #phenotypes,sample_idx = dataset.getPhenotypes(intersection=False)

    pheno_sample_select = np.ones(pheno_reader.sample_ID.shape[0], dtype=bool)
    phenotypes,pheno_sample_idx = pheno_reader.getPhenotypes(sample_idx=pheno_sample_select)
    expr_sample_select = np.ones(expr_reader.sample_ID.shape[0], dtype=bool)
    expr,expr_sample_idx = expr_reader.getPhenotypes(sample_idx=expr_sample_select)

    snps = geno_reader.getGenotypes()
    position = geno_reader.getPos()
    position,chromBounds = data_util.estCumPos(position=position,offset=0)
    gid_start,gid_end = geno_reader.getGenoIndex(chrom=5,pos_start=(5,20008540),pos_end=(5,20008570))
    gid_range = np.arange(gid_start,gid_end+1)

    for ig,g_ID in enumerate(gid_range):

        g_ID = gid_range[ig:(ig+1)]
        print(g_ID)
        gs_idx = dataset.sample_idx["geno"].values
        ps_idx = dataset.sample_idx["pheno"].values
        egs_idx = exprset.sample_idx["geno"].values
        eps_idx = exprset.sample_idx["pheno"].values

        snps_sub = snps[np.ix_(gs_idx,g_ID)][:,0]
        phenotypes_sub = phenotypes.values[ps_idx]
        esnps_sub = snps[np.ix_(egs_idx,g_ID)][:,0]
        tbl4_sub = expr['AT5G49430'].values[eps_idx]
        position_sub = position.iloc[[g_ID[0]]]
        print(position_sub)
        
        point_file = graphics_prefix+'point_chr%d_%d.png'%(position_sub['chrom'],position_sub['pos'])        
        fig = plt.figure(figsize=[5,2.5])#create the figure
        plt.subplot(1,2,1)
        plt.plot(snps_sub+0.05*np.random.randn(snps_sub.shape[0]),phenotypes_sub ,'.')
        plt.xlabel("SNP")
        plt.ylabel("Seed size")
        plt.subplot(1,2,2)
        plt.plot(esnps_sub+0.05*np.random.randn(esnps_sub.shape[0]),tbl4_sub ,'.')
        plt.xlabel("SNP")
        plt.ylabel("TBL4 expression")
        plt.tight_layout()
        fig.savefig(point_file)
        plt.close(fig)

        bxp_file = graphics_prefix+'bxp_chr%d_%d.png'%(position_sub['chrom'],position_sub['pos'])        
        fig = plt.figure(figsize=[5,2.5])#create the figure
        plt.subplot(1,2,1)
        phenotypes_box = [phenotypes_sub[snps_sub==0],phenotypes_sub[snps_sub==2]]
        plt.boxplot(phenotypes_box)
        plt.xlabel("SNP")
        plt.ylabel("Seed size")
        plt.subplot(1,2,2)
        tbl4_box = [tbl4_sub[esnps_sub==0],tbl4_sub[esnps_sub==2]]
        plt.boxplot(tbl4_box)
        plt.xlabel("SNP")
        plt.ylabel("TBL4 expression")
        plt.tight_layout()
        fig.savefig(bxp_file)
        plt.close(fig)
        
if __name__=='__main__':
    main()

