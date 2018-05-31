library(data.table)
library(plyr)

oberon_1001_data = file.path(Sys.getenv('OBERON_DATA'),'data1/shhuang/data/1001_genomes')
metadata_path = file.path(PROJ_DEVEL_PATH,'metadata')
tfq_gale2_path = file.path(oberon_1001_data,'tfq_gale2')
hchen_1001_path = file.path(file.path(Sys.getenv('OBERON_DATA'),'data3/1001'))

hchen_tbam_path = file.path(hchen_1001_path,'tbam')
sv_data_path = file.path(PROJ_DATA_PATH,'SV_data-2016-02-09')
sv_data2_path = file.path(PROJ_DATA_PATH,'SV_data-2016-03-03')
pe_gene_path = file.path(PROJ_DATA_PATH,'flip')
strg_01_path = file.path(PROJ_RESULTS_PATH,"stringtie_01")
acc_list_master = fread(file.path(metadata_path,'A.thaliana_master_accession_list_20150623.csv'))
acc_list_master[,tg_ecotypeid:=as.character(tg_ecotypeid)]
acc_list_1001g = fread(file.path(metadata_path,'1001genomes-accessions.csv'))
acc_list_1001g = acc_list_1001g[,tg_ecotypeid:=as.character(id)]
acc_list_1001g[,tg_ecotypeidx:=paste0('X',tg_ecotypeid)]
avgm_path = file.path(PROJ_DATA_PATH,'average_methylation')
avgm_file = file.path(avgm_path,'average_methylation_20160413.tsv')
avgm = fread(avgm_file)
setnames(avgm,'Sequenced by','Sequenced_By')
avgm = avgm[,tg_ecotypeidx:=paste0("X",Ecotype_id)]

gmi_rel_path = file.path(PROJ_DATA_PATH,'gmi_release_v3.1')
dmc_bins_filtered = c("dmC_filtered","dmCG_filtered","dmCH_filtered")
dmc_bins_gat01_path = file.path(PROJ_RESULTS_PATH,'dmC_bins_gat01')
dmc_bins_gat02_path = file.path(PROJ_RESULTS_PATH,'dmC_bins_gat02')
dms_tables = c("mCG_table","mCHG_table","mCHH_table")
dms_path = file.path(PROJ_DATA_PATH,'dms')
dms_gat01_path = file.path(PROJ_RESULTS_PATH,'dms_gat01')
dmr_pcc_path = file.path(PROJ_DATA_PATH,'DMR_pearson_matrix')
eqtl_09_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_09')
eqtl_10_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_10')
eqtl_11_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_11')
eqtl_12_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_12')
eqtl_13_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_13')
eqtl_14_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_14')
eqtl_15_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_15')
eqtl_16_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_16')
eqtl_18_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_18')
eqtl_19_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_19')
eqtl_20_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_20')
eqtl_20_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_20')
eqtl_22_path = file.path(PROJ_RESULTS_PATH,'marginal_test_cov_22')
qtl_multiinter01_path = file.path(PROJ_RESULTS_PATH,'qtl_multiinter01')
qtl_multiinter02_path = file.path(PROJ_RESULTS_PATH,'qtl_multiinter02')
cistrome_olap01_path = file.path(PROJ_RESULTS_PATH,'cistrome_olap_01')
cistrome_olap02_path = file.path(PROJ_RESULTS_PATH,'cistrome_olap_02')

snp_gat01_path = file.path(PROJ_RESULTS_PATH,'snp_gat01')

rmake_bam_mapped_file = file.path(PROJ_RESULTS_PATH,'rmake_stats_2016-02-09/rmake_stats_2016-02-09-bam_mapped.txt')
rmake_fq_keep1_file = file.path(PROJ_RESULTS_PATH,'rmake_stats_2016-04-19/rmake_stats_2016-04-19-mapped_by_fq_keep1.tsv')
                    
get_samples_from_dir<-function(data_dir=tfq_gale2_path,src='salk') {
    dirnames = dir(data_dir)[file.info(dir(data_dir,full.names=TRUE))$isdir]
    if (src=='salk') {
        sample_id = grep('(10C$|16C$|CS$)',dirnames,invert=TRUE,value=TRUE)
    } else if (src=='all') {
        sample_id = dirnames
    } else {
        sample_id = dirnames
    }
    return(sample_id)
}


load_dmr_pcc<-function() {    
    pcc_files = c("C_DMR_pearson","CG_DMR_pearson","CH_DMR_pearson")
    dmr_pcc0 = alply(pcc_files,1,function(pf) {
        x = load(file.path(dmr_pcc_path,pf))
        list(name=x,value=get(x))
    })
    dmr_pcc = llply(dmr_pcc0,function(x) x[['value']])
    names(dmr_pcc) = laply(dmr_pcc0,function(x) x[['name']])
    return(dmr_pcc)
}

pd_1001g = data.table(
                tg_ecotypeidx=acc_list_1001g[,tg_ecotypeidx],
                tg_ecotypeid=acc_list_1001g[,tg_ecotypeid],
                dirname=acc_list_1001g[,tg_ecotypeid],
                latitude=acc_list_1001g[,latitude],
                longitude=acc_list_1001g[,longitude],
                group=acc_list_1001g[,group],
                icon=acc_list_1001g[,icon],
                country=acc_list_1001g[,country]
                )
pd_master = acc_list_master[,list(tg_ecotypeidx=paste0("X",tg_ecotypeid),
                                  tg_ecotypeid=tg_ecotypeid,
                                  latitude=latitude,
                                  longitude=longitude
                                  )]

salk_tx_info = fread(file.path(PROJ_RESULTS_PATH,
                               'make_acc_tx_list_2015-12-05/graphics/',
                               'make_acc_tx_list_2015-12-05-salk_expt_info.txt'))

salk_tx_info[,fc_short:=unlist(lapply(strsplit(salk_tx_info[,fc],"_"),function(arr) 
  paste(arr[1:3],collapse="_")))]
salk_tx_info[,group:=with(pd_1001g,group[match(salk_tx_info[,tg_ecotypeid],tg_ecotypeid)])]

salk_tx_platebatch = salk_tx_info[,list(plate_batch1=paste(sort(unique(Plate)),collapse=",")),
                                  by="tg_ecotypeid"]

# Mark's plate in one batch, Rosa's first and second growth batch
plate_batch2_list = c("b1",rep("b2",9),rep("b3",3))
names(plate_batch2_list) = c("MU",1:9,paste0("#",1:3,"_+_reads"))
# Mark's plate in one batch, Rosa's individual plates as different batches, combined
plate_batch3_list = 1:13
names(plate_batch3_list) = c("MU",1:9,paste0("#",1:3,"_+_reads"))
# Mark's plate in one batch, Rosa's individual plates as different batches
plate_batch4_list = c(1,1:9,1:3)
names(plate_batch4_list) = c("MU",1:9,paste0("#",1:3,"_+_reads"))
salk_tx_platebatch[,`:=`(plate_batch2=plate_batch2_list[plate_batch1],
                         plate_batch3=plate_batch3_list[plate_batch1],
                         plate_batch4=plate_batch4_list[plate_batch1])]
salk_tx_info[,`:=`(plate_batch2=plate_batch2_list[Plate])]
salk_tx_platebatch[,group:=with(pd_1001g,group[match(salk_tx_platebatch[,tg_ecotypeid],tg_ecotypeid)])]
salk_tx_replist = dlply(salk_tx_info[,.N,by=list(tg_ecotypeid,Plate)],
                        "Plate",function(df) { unique(df[,'tg_ecotypeid'])})

salk_tx_grpplt = dlply(salk_tx_info[,.N,by=list(group,tg_ecotypeid,Plate)],
                         c("group","Plate"),function(df) { unique(df[,'tg_ecotypeid'])})

salk_tx_grpbat = dlply(salk_tx_info[,.N,by=list(group,plate_batch2,tg_ecotypeid)],
                       c("group","plate_batch2"),function(df) { unique(df[,'tg_ecotypeid'])})

salk_tx_grp = dlply(salk_tx_info[,.N,by=list(group,tg_ecotypeid)],
                    c("group"),function(df) { unique(df[,'tg_ecotypeid'])})

salk_mc_rosa_batch1 = fread(file.path(PROJ_DEVEL_PATH,'metadata','1001_accessions_Rosa_batch1.csv'))
salk_mc_rosa_batch1 = salk_mc_rosa_batch1[,tg_ecotypeidx:=paste0("X",Accession)]
salk_mc_rosa_batch2 = fread(file.path(PROJ_DEVEL_PATH,'metadata','1001_accessions_Rosa_batch2.csv'))
salk_mc_rosa_batch2 = salk_mc_rosa_batch2[,tg_ecotypeidx:=paste0("X",Accession)]

# Ath 21 ref genes from 19 accessions
ath_19acc_21ref_file = file.path(Sys.getenv('SHHUANG_DATA_PATH'),'arabidopsis/Wang_SciRep2014_19AccRefGenes/Table1_RefGenes.csv')
tair10_genelength_file = file.path(Sys.getenv('SHHUANG_DATA_PATH'),'genomes/TAIR10/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes_transposons.fixed2.chr.gene_length')

# normalized expression 
gNorm_normCounts_k4_file1 = file.path(PROJ_RESULTS_PATH,'ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_normCounts_k4_vst2.txt')
gNorm_normCounts_k4_cv5_file1 = file.path(PROJ_RESULTS_PATH,'ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_normCounts_k4_vst2_cv0p05_1001g.txt2')
gNorm_normCounts_k4_file2 = file.path(PROJ_RESULTS_PATH,'ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4_vst2.tsv')
gNorm_normCounts_k4_cv5_file2 = file.path(PROJ_RESULTS_PATH,'ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4_vst2_cv0p05_1001g.tsv2')


# DAP properties
dap_pf4_summary_file = file.path(Sys.getenv('SHHUANG_PROJECTS_PATH'),'dap/analysis.v4/create_dap_web_dbtable_2015-10-21/create_dap_web_dbtable_2015-10-21-pf4_summary.txt')

dap_pf5_target_file = file.path(Sys.getenv('SHHUANG_PROJECTS_PATH'),'dap/analysis.v4/dap_distr/tf_target_v4_pf5/pf5_TF_target_olap1k.txt')
                           
atwell_phenotype_file = file.path(Sys.getenv('SHHUANG_DATA_PATH'),'arabidopsis/atpolydb/misc_data/phenotype_published_raw.tsv')

ecotype_list_nomatch_file = file.path(PROJ_DATA_PATH,'strain_verify','ecotype_list_nomatch_2017-02-17_tr.txt')
ecotype_list_mixup_file = file.path(PROJ_DATA_PATH,'strain_verify','ecotype_list_mixup_2017-02-17_tr.txt')

snpmatch_1001g_01_path = file.path(PROJ_RESULTS_PATH,'snpmatch_1001G_01')
