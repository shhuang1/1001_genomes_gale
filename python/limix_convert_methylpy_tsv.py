import h5py
import numpy as np
import pandas
from optparse import OptionParser
import scipy as sp
import os
import sys

from shared import LoggerFactory

def convert_dmr_tsv(hdf,tsv_file,chrom,start,end,sample_subset=None):
    """Convert DMR tsv
    """
    if ((start is not None) or (end is not None) or (chrom is not None)):
        print "cannot handle start/stop/chrom boundaries for DMR tsv file"
        return
    if 'genotype' in hdf.keys():
        del(hdf['genotype'])
    genotype = hdf.create_group('genotype')
    col_header = genotype.create_group('col_header')
    row_header = genotype.create_group('row_header')
    
    C = pandas.io.parsers.read_csv(tsv_file,sep='\t',header=0,index_col=False)
    chrom = C['chr']
    pos = C['start']
    sample_ID = [ml.replace("methylation_level_","") for ml in C.columns[4:]]
    
    row_header.create_dataset(name='sample_ID',data=sample_ID)
    col_header.create_dataset(name='chrom',data=chrom)
    col_header.create_dataset(name='pos',data=pos)
    
    snps = C.ix[:,4::].T.as_matrix()
    genotype.create_dataset(name='matrix',data=snps,\
                            chunks=(snps.shape[0],min(10000,snps.shape[1])),\
                            compression='gzip')

def convert_dms_tsv(hdf,tsv_file,chrom,start,end,sample_subset=None):
    """Convert DMR tsv
    """
    if ((start is not None) or (end is not None) or (chrom is not None)):
        print "cannot handle start/stop/chrom boundaries for DMS tsv file"
        return
    if 'genotype' in hdf.keys():
        del(hdf['genotype'])
    genotype = hdf.create_group('genotype')
    col_header = genotype.create_group('col_header')
    row_header = genotype.create_group('row_header')
    
    C = pandas.io.parsers.read_csv(tsv_file,sep='\t',header=0,index_col=False)
    chrom = C['chr']
    pos = C['pos']
    sample_ID = np.array(C.columns[5:],dtype='string')
    
    row_header.create_dataset(name='sample_ID',data=sample_ID)
    col_header.create_dataset(name='chrom',data=chrom)
    col_header.create_dataset(name='pos',data=pos)
    
    snps = C.ix[:,5::].T.as_matrix()
    genotype.create_dataset(name='matrix',data=snps,\
                            chunks=(snps.shape[0],min(10000,snps.shape[1])),\
                            compression='gzip')

def main():

    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-O","--outfile", action="store", dest='outfile', type=str, help='The output hdf5 file wiht the resulting data', default="example_out")
    parser.add_option("-R","--dmr",action="store", dest='dmr', help="Read DMR tsv file (filename)", default=None)
    parser.add_option("-S","--dms",action="store", dest='dms', help="Read methylation site tsv file (filename)", default=None)
    
    (options, args) = parser.parse_args()
    logger = LoggerFactory.get_logger(options.outfile+'.log')
    LoggerFactory.log_command(logger,sys.argv)

    hdf = h5py.File(options.outfile)
    if options.dmr is not None:
        logger.info('Converting DMR tsv %s and write to %s',options.dmr,options.outfile)
        convert_dmr_tsv(hdf,options.dmr,chrom=None,start=None,end=None,sample_subset=None)
    if options.dms is not None:
        logger.info('Converting DMS tsv %s and write to %s',options.dms,options.outfile)
        convert_dms_tsv(hdf,options.dms,chrom=None,start=None,end=None,sample_subset=None)


    logger.info('Done!')
    
    
if __name__=='__main__':
    main()

