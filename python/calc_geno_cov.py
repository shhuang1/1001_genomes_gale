# calculate genotype covariance
import os
import h5py
import itertools as it
import sys
import csv
import limix
import logging
import numpy as np
import scipy as sp
import scipy.stats as sps

# import LIMIX
import limix.io.genotype_reader as gr
import limix.utils.util_functions as util_functions

from shared import LoggerFactory

def strip_x(tg_ecotypeidx): return tg_ecotypeidx.replace("X","")
strip_xvec = np.vectorize(strip_x)

def add_x(tg_ecotypeid): return "X"+tg_ecotypeid
add_xvec = np.vectorize(add_x)

def save_cov_in_text_format(filename, cov, accessions):
    with open(filename, 'w') as f:
        writer = csv.writer(f,delimiter='\t')
        for acc, row in it.izip(accessions, cov):
            writer.writerow([acc]+row.tolist())

def main():

    geno_in,norm_cov,cov_out_hdf5,cov_out_csv = sys.argv[1:]

    logger = LoggerFactory.get_logger(cov_out_hdf5+'.log')
    LoggerFactory.log_command(logger,sys.argv[1:])
    
    ## Import genotype data
    logger.info('Loading genotype from %s',geno_in)
    geno_reader = gr.genotype_reader_tables(geno_in)
    if (norm_cov=='1'):
        logger.info('Normalizing')
        norm = True
    else:
        logger.info('NOT normalizing')
        norm = False

    sample_relatedness = geno_reader.getCovariance(normalize=norm)

    logger.info('Saving covariance to HDF5 file %s',cov_out_hdf5)
    out_dict = {'Cov':sample_relatedness}
    o = h5py.File(cov_out_hdf5,'w')
    util_functions.smartDumpDictHdf5(out_dict,o)
    o.close()
    
    logger.info('Saving covariance to CSV file %s',cov_out_csv)
    save_cov_in_text_format(cov_out_csv,sample_relatedness,geno_reader.sample_ID)

    logger.info('Done!')

if __name__=='__main__':
    main()
