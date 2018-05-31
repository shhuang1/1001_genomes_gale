# coding: utf-8
import sys
import csv
import h5py
import itertools as it
import logging

import numpy as np
import scipy as SP
from pygwas.core import genotype
from pygwas.core import kinship
from pygwas.core import phenotype

from shared import LoggerFactory

def save_kinship_in_text_format(filename, k, accessions):
    with open(filename, 'w') as f:
        writer = csv.writer(f,delimiter='\t')
        for acc, row in it.izip(accessions, k):
            writer.writerow([acc]+row.tolist()[0])

def main():

    geno_in,scale_kinship,kinship_out_hdf5,kinship_out_csv = sys.argv[1:]

    logger = LoggerFactory.get_logger(kinship_out_hdf5+'.log')
    LoggerFactory.log_command(logger,sys.argv[1:])
    
    ## Import genotype data
    geno = genotype.load_hdf5_genotype_data(geno_in)
    SNP_acc = geno.accessions
    logger.info('Finished reading SNP from %s',geno_in)

    logger.info('Start calculating kinship')
    K = geno.get_ibs_kinship_matrix()
    if (scale_kinship=='1'):
        logger.info('Scaling')
        K = kinship.scale_k(K)
    else:
        logger.info('NOT scaling')

    logger.info('Saving kinship to HDF5 file %s',kinship_out_hdf5)
    kinship.save_kinship_to_file(kinship_out_hdf5,K,geno.accessions,geno.num_snps)
    logger.info('Saving kinship to CSV file %s',kinship_out_csv)
    save_kinship_in_text_format(kinship_out_csv,K,geno.accessions)

    logger.info('Done!')

if __name__=='__main__':
    main()
