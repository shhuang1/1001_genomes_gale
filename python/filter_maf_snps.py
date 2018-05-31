# coding: utf-8
import sys
import csv
import h5py
import logging

import numpy as np
import scipy as SP
from pygwas.core import genotype
from pygwas.core import kinship
from pygwas.core import phenotype

from shared import LoggerFactory

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
    
def main():

    geno_in,acc_in,maf_lb,maf_ub,geno_out = sys.argv[1:]

    logger = LoggerFactory.get_logger(geno_out+'.log')
    LoggerFactory.log_command(logger,sys.argv[1:])
    
    maf_lb,maf_ub = float(maf_lb),float(maf_ub)

    ## Import genotype data
    geno = genotype.load_hdf5_genotype_data(geno_in)
    SNP_acc = geno.accessions
    logger.info('Finished reading SNP from %s',geno_in)

    ## accession subset
    with open(acc_in, 'rb') as f:
        reader = csv.reader(f)
        file_acc = list(reader)
    logger.info('Finished reading accession subset from %s',acc_in)
    
    ## get common accessions in the same order for genotype and accession subset
    acc_common = [acc for acc in SNP_acc if acc in file_acc]

    ## filtering
    logger.info('Start subsetting accessions and filtering SNPs by MAF >%f and <=%f',maf_lb,maf_ub)
    match = lambda a, b: [ b.index(x) if x in b else None for x in a ]
    geno.filter_accessions_ix(match(acc_common,SNP_acc))
    (num_snps,num_removed) = filter_maf_snps(geno,maf_lb,maf_ub)
    logger.info('Removed %d from %d SNPs',num_removed,num_snps)
    logger.info('Number of SNPs remaining %d',geno.num_snps)
    
    logger.info('Start writing filtered genotype file to %s',geno_out)
    geno.save_as_hdf5(geno_out)
    logger.info('Finished')

    logger.info('Done!')

if __name__=='__main__':
    main()
