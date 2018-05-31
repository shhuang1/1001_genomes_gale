# coding: utf-8
# write SNP positions from a genotype file
import os
import sys
import itertools as it
import limix
import logging
import numpy as np

# import LIMIX
import limix.modules.qtl as qtl
import limix.io.data as data
import limix.io.data_util as data_util
import limix.io.genotype_reader as gr
import limix.io.phenotype_reader as phr

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

def main():

    geno_file,pheno_file,out_dir = sys.argv[1:]
    
    #geno_file = '/gale/netapp/home/shhuang/data/1001_genomes/gmi_release_v3.1/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5'
    #pheno_file = '/gale/netapp/home/shhuang/projects/1001_genomes/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_1001g_vst_cv0p05T.hdf5'
    #out_dir = '.'
    logger = LoggerFactory.get_logger(os.path.join(out_dir,'get_geno_pos.log'))
    LoggerFactory.log_command(logger,sys.argv[1:])
 
    logger.info('Loading genotype from %s',geno_file)
    geno_reader = gr.genotype_reader_tables(geno_file)
    logger.info('Loading phenotype from %s',pheno_file)
    pheno_reader = phr.pheno_reader_tables(pheno_file)
    pheno_reader.sample_ID = strip_xvec(pheno_reader.sample_ID)

    # the data object allows to query specific genotype or phenotype data
    logger.info('Creating QTL dataset')
    dataset = data.QTLData(geno_reader=geno_reader,pheno_reader=pheno_reader)
    # getting genotypes
    #snps = dataset.getGenotypes() #SNPS
    position = dataset.getPos()
    position,chromBounds = data_util.estCumPos(position=position,offset=100000)

    logger.info('Writing output to directory %s',out_dir)
    position = position.astype(int)
    chromBounds  = chromBounds.astype(int)
    position.to_csv(os.path.join(out_dir,'position.txt'),header=True,index=False,sep='\t')
    np.savetxt(os.path.join(out_dir,'chromBounds.txt'),chromBounds,delimiter=",")
    

    
if __name__=='__main__':
    main()

