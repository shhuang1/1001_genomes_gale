library(data.table)
library(plyr)
source('get_file_paths.R')

this_analysis_path = file.path(PROJ_RESULTS_PATH,"make_acc_tx_list_2015-12-05")
prefix.string = file.path(this_analysis_path,"graphics","make_acc_tx_list_2015-12-05-")
step_list = ""

gdoc_files = c('1001_project_new_set - Transcriptomes_051914_clean.csv',
               '1001_project_new_set - Transcriptomes_060914_clean.csv',
               '1001_project_new_set - Transcriptomes_072514_clean.csv',
               '1001_project_new_set - Transcriptomes_080714_clean.csv',
               '1001_project_new_set - Transcriptomes_081514_clean.csv',
               '1001_project_new_set - Transcriptomes_102314_clean.csv',
               '1001_project_new_set - Transcriptomes_102914_clean.csv',
               '1001_project_new_set - Transcriptomes_111414_clean.csv',
               '1001_project_new_set - Transcriptomes_111914_clean.csv',
               '1001_project_new_set - Transcriptomes_+_reads_012715_clean.csv',
               '1001_project_new_set - Transcriptomes_+_reads_020515_clean.csv')

acc_from_gdoc = rbindlist(alply(gdoc_files,1,function(f) {
    dt = fread(file.path(metadata_path,f),sep='\t')
    dt[,1:8,with=FALSE]
}))

acc_from_dir = dir(tfq_gale2_path)[file.info(dir(tfq_gale2_path,full.names=T))$isdir]

# In gdoc but not in directory
# NA 9810 9229 8321 9061
# no FASTQ: 9810, 9061
# have FASTQ but not aligned 9229, 8321
setdiff(acc_from_gdoc[,'Accession Number',with=FALSE][[1]],acc_from_dir)
 
# In directory but not in gdoc
# This is OK since there are accessions sequenced before the 1001 new set file
setdiff(acc_from_dir,acc_from_gdoc[,'Accession Number',with=FALSE][[1]])

# colorspace from Bob's paper
acc_cs = fread(file.path(metadata_path,'tbam_CS_151205.txt'),header=FALSE)
acc_split_cs = strsplit(acc_cs[[1]],'_')
acc_out_cs = data.table('sample_id'=acc_cs[[1]],
                        'tg_ecotypeid'=sapply(acc_split_cs,"[",1),
                        'condition'="",
                        'source'='Nat2013',
                        'seq'='SOLiD')

acc_split_dir = strsplit(acc_from_dir,'-')
# make list according to folders in the alignment directory
acc_out_dir = data.table('sample_id'=acc_from_dir,
                        'tg_ecotypeid'=sapply(acc_split_dir,"[",1),
                        'condition'=sapply(acc_split_dir,function(e) ifelse(length(e)>1,e[2],"")),
                        'source'=sapply(acc_split_dir,function(e) ifelse(length(e)>1,"GMI","SALK")),
                        'seq'='Illumina')

acc_out_dir2 = rbind(acc_out_dir,acc_out_cs)

write.table(acc_out_dir2,
            paste0(prefix.string,'acc_tx_table.tsv'),sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE)
                        
