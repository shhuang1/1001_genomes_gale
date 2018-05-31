import csv
import glob
import gzip
import itertools
import os
from pathlib import Path
import re

# Define SEPINC as shell environment variable
SNAKERULES = os.environ.get('SNAKERULES')
PROJ_NAME = "1001_genomes"
DEPDIR = "deps"
include:
    SNAKERULES + "/Doc.defs.top"

PROJ_PYTHON_PATH_GALE = DEVEL_PATH_GALE + "/1001_genomes_gale/python"
########################
# paths for 1001 data
########################
ATH_AUTOSOMES = list(range(1,6))
TFQ_OBERON5 = OBERON_DATA + "/data1/shhuang/data/1001_genomes/tfq_oberon5"
TFQ_GALE1 = OBERON_DATA + "/data1/shhuang/data/1001_genomes/tfq_gale1"
TFQ_GALE2 = OBERON_DATA + "/data1/shhuang/data/1001_genomes/tfq_gale2"
TFQ_GALE3 = OBERON_DATA + "/data1/shhuang/data/1001_genomes/tfq_gale3"
Col0_tx_fastq = "/gale/netapp/seq2/illumina_runs/160318_BARB_4112_AH3CFJBBXX/Unaligned_single_index/1001_transcriptomes"
HCHEN_1001 = OBERON_DATA + "/data3/hchen/1001"
RESULTS_ALLC = HCHEN_1001 + "/results"
TAIJI_DMR = OBERON_DATA + "/data4/riverwin/1001/DMRs"
TAIJI_dmC_bins = OBERON_DATA + "/data4/riverwin/1001/dmC_bins"

BAM_MERGE_01 =  PROJ_RESULTS_PATH_GALE + "/bam_merge_01"
BEDANNO_DAP_01 = PROJ_RESULTS_PATH_GALE + '/bedanno_dap_01'
BEDANNO_DAP_02 = PROJ_RESULTS_PATH_GALE + '/bedanno_dap_02'
CISTROME_OLAP_01 = PROJ_RESULTS_PATH_GALE + '/cistrome_olap_01'
CISTROME_OLAP_02 = PROJ_RESULTS_PATH_GALE + '/cistrome_olap_02'
GENOME_MATRIX_DIR = OBERON_DATA + '/data1/shhuang/data/1001_genomes/genome_matrix'
DMC_BINS = PROJ_DATA_PATH_GALE + '/dmC_bins'
DMC_BINS_PREFIX = ['dmCH_filtered','dmCG_filtered','dmC_filtered']
DMC_BINS_GAT01 = PROJ_RESULTS_PATH_GALE + '/dmC_bins_gat01'
DMC_BINS_GAT02 = PROJ_RESULTS_PATH_GALE + '/dmC_bins_gat02'
DMC_BINS_NEW = PROJ_DATA_PATH_GALE + '/dmC_bins_new'
DMS = PROJ_DATA_PATH_GALE + '/dms'
DMS_TABLES_PREFIX = ['mCG_table','mCHG_table','mCHH_table']
DMS_GAT01 = PROJ_RESULTS_PATH_GALE + '/dms_gat01'
MQTL = PROJ_DATA_PATH_GALE + '/mQTL'
MQTL_GAT01 = PROJ_RESULTS_PATH_GALE + '/mQTL_gat01'
METADATA_DIR = PROJ_DEVEL_PATH + '/metadata'
MPI_RELEASE_v31 = PROJ_DATA_PATH_GALE + '/gmi_release_v3.1'
MPI_RELEASE_v31_VCFGZ = MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN.vcf.gz'
SNPEFF_DAP_01 = PROJ_RESULTS_PATH_GALE + '/snpeff_dap_01'
SNP_GAT01 = PROJ_RESULTS_PATH_GALE + '/snp_gat01'
STRINGTIE_01 = PROJ_RESULTS_PATH_GALE + "/stringtie_01"
MARGINAL_TEST_COV_05 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_05'
MARGINAL_TEST_COV_06 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_06'
MARGINAL_TEST_COV_07 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_07'
MARGINAL_TEST_COV_08 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_08'
MARGINAL_TEST_COV_09 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_09'
MARGINAL_TEST_COV_10 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_10'
MARGINAL_TEST_COV_11 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_11'
MARGINAL_TEST_COV_12 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_12'
MARGINAL_TEST_COV_13 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_13'
MARGINAL_TEST_COV_14 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_14'
MARGINAL_TEST_COV_15 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_15'
MARGINAL_TEST_COV_16 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_16'
MARGINAL_TEST_COV_17 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_17'
MARGINAL_TEST_COV_18 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_18'
MARGINAL_TEST_COV_19 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_19'
MARGINAL_TEST_COV_20 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_20'
MARGINAL_TEST_COV_21 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_21'
MARGINAL_TEST_COV_22 = PROJ_RESULTS_PATH_GALE + '/marginal_test_cov_22'
TX_FASTQC_01 = PROJ_RESULTS_PATH_GALE + '/tx_fastqc_01'

ECOTYPE_LIST = list(filter(lambda y:os.path.isdir(os.path.join(TFQ_GALE2,y)),os.listdir(TFQ_GALE2)))
rule ecotype_list:
    shell: "echo {ecotype_list}"
ecotype_pat = re.compile('.+\-(16C|10C)$')
ECOTYPE_LIST_SALK = list(filter(lambda x:ecotype_pat.match(x)==None,ECOTYPE_LIST))
rule ecotype_list_salk:
    shell: "echo {ECOTYPE_LIST_SALK}"
ECOTYPE_LIST_G = [re.sub(".vcf.gz$","",re.sub("^intersection_","",os.path.basename(p))) for p in glob.glob(MPI_RELEASE_v31+'/intersection_snp_short_indel_vcf/*.vcf.gz')]

def acc_all(path,ecotype_list=ECOTYPE_LIST):
    return [path.format(tg_ecotypeid=tg_ecotypeid) for tg_ecotypeid in ecotype_list]
def acc_all_dep(path,path_dep,ecotype_list=ECOTYPE_LIST):
    return [path.format(tg_ecotypeid=tg_ecotypeid) for tg_ecotypeid in ecotype_list if os.path.isfile(path_dep.format(tg_ecotypeid=tg_ecotypeid))]

rule x:
    shell: "echo {GALE_HOME} {THUMPER1_HOME}"

##############################
# interaction with other data
##############################
DAP_V4_FAMCLUST_01 = DAP_RESULTS_PATH_GALE + "/analysis.v4/family_cluster_01"
DAP_V4_REPMASTER_01 = DAP_RESULTS_PATH_GALE + "/analysis.v4/gem07_rep_master"
#DAP_V4_FAMCLUST_01 = "analysis.v4/family_cluster_01"
#DAP_V4_FAMCLUST_01_FAM = ['bZIP_tnt']
DAP_V4_FAMCLUST_01_FAM = [d for d in os.listdir(DAP_V4_FAMCLUST_01) if os.path.isdir(os.path.join(DAP_V4_FAMCLUST_01,d))]
def dap_fam_all(path,expSys=None):
    if (expSys==None):
        return [path.format(family_expSys=family_expSys) for family_expSys in DAP_V4_FAMCLUST_01_FAM]
    else:
        es_pat = re.compile('_'+expSys+'$')
        return [path.format(family_expSys=family_expSys) for family_expSys in DAP_V4_FAMCLUST_01_FAM if (es_pat.search(family_expSys)!=None)]

def dap_fam_all_dep(path,path_dep,expSys=None):
    if (expSys==None):
        return [path.format(family_expSys=family_expSys) for family_expSys in DAP_V4_FAMCLUST_01_FAM if os.path.isfile(path_dep.format(family_expSys=family_expSys))]
    else:
        es_pat = re.compile('_'+expSys+'$')
        return [path.format(family_expSys=family_expSys) for family_expSys in DAP_V4_FAMCLUST_01_FAM if (es_pat.search(family_expSys)!=None) and os.path.isfile(path_dep.format(family_expSys=family_expSys))]

#############
# FASTQC
#############
rule tx_fastqc_01a:
    input: glob.glob(os.path.join(Col0_tx_fastq,'*.fastq.gz'))
    params: D=os.path.join(TX_FASTQC_01,'6909')
    output: os.path.join(TX_FASTQC_01,'6909/fastqc.log')
    shell: """
        module load java/jdk-7u9-linux-x64 fastqc/v0.11.5
        mkdir -p {params.D}
        fastqc --outdir={params.D} {input}
    """

rule tx_fastqc_01b:
    input: TFQ_GALE2 + '/{tg_ecotypeid}/{tg_ecotypeid}.bam.gene.count'
    params: inD=TFQ_GALE2 + '/{tg_ecotypeid}/',outD=os.path.join(TX_FASTQC_01,'{tg_ecotypeid}')
    output: TX_FASTQC_01 + '/{tg_ecotypeid}/{tg_ecotypeid}.fatqc.log'
    shell: """
        module load java/jdk-7u9-linux-x64 fastqc/v0.11.5
        mkdir -p {params.outD}
        fastqc --outdir={params.outD} {params.inD}/*.fastq.gz | tee {output}
    """
rule tx_fastqc_01b_out:
    input: acc_all(TX_FASTQC_01 + '/{tg_ecotypeid}/{tg_ecotypeid}.fatqc.log',ecotype_list=ECOTYPE_LIST_SALK)
rule tx_fastqc_01b_out2:
    shell: "echo {rules.tx_fastqc_01b_out.input}"

#################
# STAR alignment
#################
rule rmake_part1:
    input: TFQ_GALE2 + '/{tg_ecotypeid}/makefile'
    params: N="star_{tg_ecotypeid}",L="2:00:00",J="y",qname="gale.q",D=TFQ_GALE2+'/{tg_ecotypeid}'
    log: TFQ_GALE2 + '/{tg_ecotypeid}/rmake_part1.log'
    output: TFQ_GALE2 + '/{tg_ecotypeid}/rmake_part1.out'
    shell: """
        module load shhuang bamtools/2.4.0 bedtools/2.19.1
        cd {TFQ_GALE2}/{wildcards.tg_ecotypeid}; make
        echo $HOSTNAME > {output}
    """
rule rmake_part1_out:
   input: acc_all(TFQ_GALE2+'/{tg_ecotypeid}/rmake_part1.out',ecotype_list=ECOTYPE_LIST)
#snakemake -j 20 -k -p rmake_part1_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V -o {log}"

#################
# counting
#################
rule rmake_part2:
    input: TFQ_GALE2 + '/{tg_ecotypeid}/makefile'
    params: N="star_{tg_ecotypeid}",L="2:00:00",J="y",qname="gale.q",D=TFQ_GALE2+'/{tg_ecotypeid}'
    log: TFQ_GALE2 + '/{tg_ecotypeid}/rmake_part2.log'
    output: TFQ_GALE2 + '/{tg_ecotypeid}/rmake_part2.out'
    shell: """
        module load shhuang bamtools/2.4.0 bedtools/2.19.1
        cd {TFQ_GALE2}/{wildcards.tg_ecotypeid}; make -k -j 6
        echo $HOSTNAME > {output}
    """
rule rmake_part2_out:
   input: acc_all(TFQ_GALE2+'/{tg_ecotypeid}/rmake_part2.out',ecotype_list=ECOTYPE_LIST)
#snakemake -j 5 -k -p rmake_part2_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V -o {log}"

##################################
# merge all BAM from each ecotype
#################################
rule bam_merge_mkdir:
    params: d=BAM_MERGE_01+"/{tg_ecotypeid}"
    output: BAM_MERGE_01 + "/{tg_ecotypeid}/mkdir.log"
    shell:
        """
        mkdir -p {params.d}
        touch {output}
        """
rule bam_merge_mkdir_out:
    input: acc_all(BAM_MERGE_01+'/{tg_ecotypeid}/mkdir.log',ecotype_list=ECOTYPE_LIST_SALK)

rule bam_merge:
    input: bam_list=lambda wildcards: glob.glob(os.path.join(TFQ_GALE2,wildcards.tg_ecotypeid,"*.bam"))
    params: N="bam_merge_01_{tg_ecotypeid}",L="2:00:00",J="y",qname="gale.q",D=BAM_MERGE_01+"/{tg_ecotypeid}"
    log: BAM_MERGE_01 + "/{tg_ecotypeid}/bam_merge.log"
    output: BAM_MERGE_01 + "/{tg_ecotypeid}/{tg_ecotypeid}.bam"
    run:
        if (len(input.bam_list)>1):
            shell("module load samtools/1.2; samtools merge {output} {input.bam_list}")
        else:
            shell("cp {input.bam_list} {output}")
rule bam_merge_out:
    input: acc_all(BAM_MERGE_01+'/{tg_ecotypeid}/{tg_ecotypeid}.bam',ecotype_list=ECOTYPE_LIST_SALK)
# snakemake -j 10 -k -p bam_merge_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V -o {log}"

#############
# stringtie
#############
rule stringtie_01_mkdir:
    params: D=STRINGTIE_01+"/{tg_ecotypeid}"
    output: STRINGTIE_01 + "/{tg_ecotypeid}/mkdir.log"
    shell:
        """
        mkdir -p {params.D}
        touch {output}
        """
rule stringtie_01_mkdir_out:
    input: acc_all(STRINGTIE_01+'/{tg_ecotypeid}/mkdir.log')
rule stringtie_01_step1:
    input: BAM_MERGE_01 + "/{tg_ecotypeid}/{tg_ecotypeid}.bam"
    params: N="strg01_{tg_ecotypeid}",L="2:00:00",J="y",qname="gale.q",ntasks_per_node=16,
            D=STRINGTIE_01+"/{tg_ecotypeid}",
            ref_ann=GENOMES_TAIR10_GENES_TRANSPOSONS_GFF_CHR_FIXED2,
            label="1001G_STRG"
    log: STRINGTIE_01 + "/{tg_ecotypeid}/strg01_s1.log"
    output: outfile=STRINGTIE_01+"/{tg_ecotypeid}/strg01_s1_out.gtf",cov_refs=STRINGTIE_01+"/{tg_ecotypeid}/strg01_s1_cov_refs.gtf"
    shell:
        """
        module load stringtie/1.1.2
        stringtie {input} -B -G {params.ref_ann} -o {output.outfile} -p {params.ntasks_per_node} -C {output.cov_refs} -l {params.label}
        """
rule stringtie_01_step1_out:
    input: acc_all_dep(STRINGTIE_01+"/{tg_ecotypeid}/strg01_s1_out.gtf",BAM_MERGE_01+'/{tg_ecotypeid}/{tg_ecotypeid}.bam',ecotype_list=ECOTYPE_LIST_SALK)
rule stringtie_01_step1_out_clean:
    shell: "rm {rules.stringtie_01_step1_out.input}"
# snakemake -j 10 -k -p stringtie_01_step1_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V -o {log}"

#########################
# various VCF processing
#########################
# vcf
rule vcf_by_chrom_01:# MAF 0.01
    input: vcf=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}.vcf.gz'
    params: maf=0.01,prefix=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1.recode.vcf.gz'
    shell: """
        module load vcftools/0.1.14 htslib/1.2.1
        vcftools --gzvcf {input.vcf} --stdout --maf {params.maf} --recode --recode-INFO-all | bgzip -c > {output}; tabix -p vcf {output}
    """
rule vcf_by_chrom_01_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1.recode.vcf.gz',chrom=[str(c) for c in ATH_AUTOSOMES])
rule vcf_all_chrom_01:# MAF 0.01
    input: vcf=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN.vcf.gz'
    params: maf=0.01,prefix=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_filter1'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_filter1.vcf.gz'
    shell: """
        module load vcftools/0.1.14 htslib/1.2.1
        vcftools --gzvcf {input.vcf} --stdout --maf {params.maf} --recode --recode-INFO-all | bgzip -c > {output}; tabix -p vcf {output}
    """
rule vcf_all_chrom_01_out:
    input: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_filter1.vcf.gz'

# tped
# rule plink_tped_by_chrom_00:# no filtering
#     input: vcf=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}.vcf.gz'
#     params: prefix=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}'
#     output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}.tped'
#     shell: """
#         module load vcftools/0.1.14
#         vcftools --gzvcf {input.vcf} --out {params.prefix} --plink-tped
#     """
# rule plink_tped_by_chrom_00_out:
#     input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}.tped',chrom=[str(c) for c in ATH_AUTOSOMES])
# bed
# rule plink_bed_by_chrom_00:
#     input: tped=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}.tped'
#     params: prefix=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}'
#     output: MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}.bed'
#     shell: """
#         module load plink/1.07-x86_64
#         plink --tfile {params.prefix} --make-bed --noweb --out {params.prefix}
#     """
# rule plink_bed_by_chrom_00_out:
#     input: expand(MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}.bed',chrom=[str(c) for c in ATH_AUTOSOMES])

# tped
rule plink_tped_by_chrom_01:# MAF 0.01
    input: vcf=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}.vcf.gz'
    params: maf=0.01,prefix=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1.tped'
    shell: """
        module load vcftools/0.1.14
        vcftools --gzvcf {input.vcf} --out {params.prefix} --maf {params.maf} --plink-tped
    """
rule plink_tped_by_chrom_01_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1.tped',chrom=[str(c) for c in ATH_AUTOSOMES])
# bed
rule plink_bed_by_chrom_01:
    input: tped=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1.tped'
    params: prefix=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1'
    output: MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1.bed'
    shell: """
        module load plink/1.07-x86_64
        plink --tfile {params.prefix} --make-bed --noweb --out {params.prefix}
    """
rule plink_bed_by_chrom_01_out:
    input: expand(MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1.bed',chrom=[str(c) for c in ATH_AUTOSOMES])
rule plink_bed_by_chrom_01_recodeA:
    input: bed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1.bed'
    params: bfile=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1',out=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1_recodeA'
    output: bed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1_recodeA.traw'
    shell: """
        module load plink/1.90beta
        plink --bfile {params.bfile} --recode A-transpose --make-bed --noweb --out {params.out}
    """
rule plink_bed_by_chrom_01_recodeA_out:
    input: expand(MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1_recodeA.traw',chrom=[str(c) for c in ATH_AUTOSOMES])    
rule plink_bed_mergelist_01:
    input: expand(MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1.bed',chrom=[str(c) for c in ATH_AUTOSOMES])
    output: mergelist=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1.bedmerge'
    run:
        mergetemplate = MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_filter1.{fmt}'
        rows = []
        mergefmt = ['bed','bim','fam']
        for c in ATH_AUTOSOMES[1:]:
            rows.append(dict((k,mergetemplate.format(chrom=c,fmt=k)) for k in mergefmt))
        ofh = open(output.mergelist,'w')
        writer = csv.DictWriter(ofh,fieldnames=mergefmt,delimiter=" ")
        writer.writerows(rows)
        ofh.close()
rule plink_bed_merge_01:
    input: mergelist=rules.plink_bed_mergelist_01.output.mergelist
    params: firstout=os.path.splitext(rules.plink_bed_mergelist_01.input[0])[0],out=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1'
    output: mergebed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1.bed'
    shell: """
        module load plink/1.07-x86_64
        plink --bfile {params.firstout} --merge-list {input.mergelist} --make-bed --noweb --out {params.out}
    """
rule plink_prune_01:
    input: bed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1.bed'
    params: bfile=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1',out=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1_independent1'
    output: prune_in=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1_independent1.prune.in',bed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1_independent1.bed'
    shell: """
        module load plink/1.07-x86_64
        plink --bfile {params.bfile} --indep 50 5 1.25 --geno 0.1 --out {params.out} --noweb
        plink --bfile {params.bfile} --noweb --make-bed --extract {output.prune_in} --out {params.out}
    """
rule plink_cluster_01:
    input: bed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1_independent1.bed'
    params: bfile=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1_independent1'
    output: MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1_independent1.genome.gz'
    shell: """
        module load plink/1.07-x86_64
        plink --bfile {params.bfile} --noweb --Z-genome --out {params.bfile}
    """
rule plink_mds_01:
    input: bed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1_independent1.bed',genome=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1_independent1.genome.gz'
    params: bfil=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1_independent1'
    shell: """
        module load plink/1.07-x86_64
        plink --bfile {params.bfile} --read-genome {input.genome} --cluster --mds-plot 20 --noweb
    """
# tped
rule plink_tped_by_chrom_02:# filtered by transcriptome individuals
    input: vcf=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}.vcf.gz'
    params: prefix=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx',tx_acc=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_acc.txt'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx.tped'
    shell: """
        module load vcftools/0.1.14
        vcftools --gzvcf {input.vcf} --keep {params.tx_acc} --out {params.prefix} --plink-tped
    """
rule plink_tped_by_chrom_02_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx.tped',chrom=[str(c) for c in ATH_AUTOSOMES])
# bed
rule plink_bed_by_chrom_02:
    input: tped=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx.tped'
    params: prefix=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx'
    output: MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx.bed'
    shell: """
        module load plink/1.90beta
        plink --tfile {params.prefix} --make-bed --noweb --out {params.prefix}
    """
rule plink_bed_by_chrom_02_out:
    input: expand(MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx.bed',chrom=[str(c) for c in ATH_AUTOSOMES])    
rule plink_bed_mergelist_02:
    input: expand(MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx.bed',chrom=[str(c) for c in ATH_AUTOSOMES])
    output: mergelist=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx.bedmerge'
    run:
        mergetemplate = MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx.{fmt}'
        rows = []
        mergefmt = ['bed','bim','fam']
        for c in ATH_AUTOSOMES[1:]:
            rows.append(dict((k,mergetemplate.format(chrom=c,fmt=k)) for k in mergefmt))
        ofh = open(output.mergelist,'w')
        writer = csv.DictWriter(ofh,fieldnames=mergefmt,delimiter=" ")
        writer.writerows(rows)
        ofh.close()
rule plink_bed_merge_02:
    input: mergelist=rules.plink_bed_mergelist_02.output.mergelist
    params: firstout=os.path.splitext(rules.plink_bed_mergelist_02.input[0])[0],out=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx'
    output: mergebed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx.bed'
    shell: """
        module load plink/1.90beta
        plink --bfile {params.firstout} --merge-list {input.mergelist} --make-bed --noweb --out {params.out}
    """
rule plink_filter_02:# only individuals in 1001g transcriptomes, max 10% missing, MAF>0.01
    input: bed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx.bed'
    params: bfile=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx',out=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1'
    output: out=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.bed'
    shell: """
        module load plink/1.90beta
        plink --bfile {params.bfile} --geno 0.1 --maf 0.01 --make-bed --noweb --out {params.out}
    """
rule plink_bed_02_recodeA:
    input: bed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.bed'
    params: bfile=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1',out=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1_recodeA'
    output: bed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1_recodeA.traw'
    shell: """
        module load plink/1.90beta
        plink --bfile {params.bfile} --recode A-transpose --make-bed --noweb --out {params.out}
    """
rule plink_filter_02_by_chrom:# only individuals in 1001g transcriptomes, max 10% missing, MAF>0.01
    input: bed=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx.bed'
    params: bfile=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx',out=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx_filter1'
    output: out=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx_filter1.bed'
    shell: """
        module load plink/1.90beta
        plink --bfile {params.bfile} --geno 0.1 --maf 0.01 --make-bed --noweb --out {params.out}
    """
rule plink_filter_02_by_chrom_out:
    input: expand(MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx_filter1.bed',chrom=[str(c) for c in ATH_AUTOSOMES])

rule vcf_all_chrom_02:# MAF 0.05
    input: vcf=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN.vcf.gz'
    params: maf=0.05,prefix=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_filter2'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_filter2.vcf.gz'
    shell: """
        module load vcftools/0.1.14 htslib/1.2.1
        vcftools --gzvcf {input.vcf} --stdout --maf {params.maf} --recode --recode-INFO-all | bgzip -c > {output}; tabix -p vcf {output}
    """
rule vcf_all_chrom_02_out:
    input: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_filter2.vcf.gz'

# filtered by seed size accessions, biallelic SNPs only
# tped
rule plink_tped_all_03:
    input: vcf=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN.vcf.gz',seed_acc=RESULTS_PATH_GALE+'/seed_size/format_seed_size_2016-10-17/format_seed_size_2016-10-17-clean02_{seed_subset}_samples.txt'
    params: prefix=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3.tped'
    shell: """
        module load vcftools/0.1.14
        vcftools --gzvcf {input.vcf} --keep {input.seed_acc} --min-alleles 2 --max-alleles 2 --remove-indels --out {params.prefix} --plink-tped
    """
rule plink_tped_all_03_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3.tped',seed_subset=['1kgen','1K'])

#bed
rule plink_bed_all_03:
    input: tped=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3.tped'
    params: prefix=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3'
    output: MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3.bed'
    shell: """
        module load plink/1.07-x86_64
        plink --tfile {params.prefix} --make-bed --noweb --out {params.prefix}
    """
rule plink_bed_all_03_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3.bed',seed_subset=['1kgen','1K'])
    
# tped recode12
rule plink_tped_all_03_recode12:
    input: tfile=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3.tped'
    params: pref_in=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3',pref_out=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3_recode12'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3_recode12.tped'
    shell: """
        module load plink/1.90beta
        plink --tfile {params.pref_in} --recode 12 transpose --output-missing-genotype 0 --out {params.pref_out}
    """
rule plink_tped_all_03_recode12_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3_recode12.tped',seed_subset=['1kgen','1K'])

# IBS kinship
rule plink_tped_all_03_aIBS:
    input: tfile=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3_recode12.tped'
    params: pref_in=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3_recode12'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3_recode12.aIBS.kinf'
    shell: """
        module load emmax/20120210
        emmax-kin-intel64 -v -s -d 10 {params.pref_in}
    """
rule plink_tped_all_03_aIBS_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3_recode12.aIBS.kinf',seed_subset=['1kgen','1K'])
# BN kinship
rule plink_tped_all_03_aBN:
    input: tfile=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3_recode12.tped'
    params: pref_in=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3_recode12'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3_recode12.aBN.kinf'
    shell: """
        module load emmax/20120210
        emmax-kin-intel64 -v -d 10 {params.pref_in}
    """
rule plink_tped_all_03_aBN_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3_recode12.aBN.kinf',seed_subset=['1kgen','1K'])


# filtered by seed size accessions, biallelic SNPs only, MAF>=0.01, call rate >=90%
# tped
rule plink_tped_all_04:
    input: vcf=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN.vcf.gz',seed_acc=RESULTS_PATH_GALE+'/seed_size/format_seed_size_2016-10-17/format_seed_size_2016-10-17-clean02_{seed_subset}_samples.txt'
    params: prefix=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4.tped'
    shell: """
        module load vcftools/0.1.14
        vcftools --gzvcf {input.vcf} --keep {input.seed_acc} --min-alleles 2 --max-alleles 2 --remove-indels --max-missing 0.9 --maf 0.01 --out {params.prefix} --plink-tped
    """
rule plink_tped_all_04_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4.tped',seed_subset=['1kgen','1K'])

rule vcf_all_04:
    input: vcf=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN.vcf.gz',seed_acc=RESULTS_PATH_GALE+'/seed_size/format_seed_size_2016-10-17/format_seed_size_2016-10-17-clean02_{seed_subset}_samples.txt'
    params: prefix=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4.vcf.gz'
    shell: """
        module load vcftools/0.1.14
        vcftools --gzvcf {input.vcf} --keep {input.seed_acc} --min-alleles 2 --max-alleles 2 --remove-indels --max-missing 0.9 --maf 0.01 --recode --recode-INFO-all --stdout | gzip -c > {output}
    """
rule vcf_all_04_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4.vcf.gz',seed_subset=['1kgen','1K'])


# tped recode12
rule plink_tped_all_04_recode12:
    input: tfile=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4.tped'
    params: pref_in=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4',pref_out=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4_recode12'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4_recode12.tped'
    shell: """
        module load plink/1.90beta
        plink --tfile {params.pref_in} --recode 12 transpose --output-missing-genotype 0 --out {params.pref_out}
    """
rule plink_tped_all_04_recode12_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4_recode12.tped',seed_subset=['1kgen','1K'])

# tped recode12
rule plink_tped_all_04_genes_01:
    input: tfile=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4.tped',\
           genes_bed=GENOMES_TAIR10_ENS20_GENES_CHR1to5
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4_genes1.txt'
    shell: """
            module load bedtools/2.26.0
            awk -v OFS='\\t' '{{ print "chr"$1,$4-1,$4,$2 }}' {input.tfile} | bedtools intersect -a stdin -b <( bedtools slop -b 1000 -g {GENOMES_TAIR10_CHROM_SIZES} -i {input.genes_bed}) -wa -wb > {output}
    """
rule plink_tped_all_04_genes_01_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4_genes1.txt',seed_subset=['1kgen','1K'])

rule plink_tped_to_ped:
    input: tfile=MPI_RELEASE_v31 + '/{tped}.tped'
    params: pref=MPI_RELEASE_v31 + '/{tped}'
    output: ofile=MPI_RELEASE_v31 + '/{tped}.ped'
    shell: """
            module load plink/1.90beta
            plink --tfile {params.pref} --recode --out {params.pref}
    """
rule plink_tped_to_ped_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_{filter_list}.ped',seed_subset=['1kgen','1K'],filter_list=['filter4_recode12'])
    
####################
# LD
####################
rule plink_ld_by_chrom_00:# no filtering
    input: bed=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}.bed'
    params: bfile=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}',out=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_r2'
    output: MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_r2.ld'
    shell: """
        module load plink/1.90beta
        plink --bfile {params.bfile} --r2 --ld-window 20 --ld-window-r2 0 --ld-window-kb 100 --out {params.out}
    """
rule plink_ld_by_chrom_00_out:
    input: expand(MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_r2.ld',chrom=[str(c) for c in ATH_AUTOSOMES])

rule plink_ld_by_chrom_01:# 1001 transcriptiome samples
    input: bed=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx.bed'
    params: bfile=MPI_RELEASE_v31 + '/X1001tx_pairwise_ld1/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx',out=MPI_RELEASE_v31 + '/X1001tx_pairwise_ld1/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx_r2'
    output: MPI_RELEASE_v31 + '/X1001tx_pairwise_ld1/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx_r2.ld'
    shell: """
        module load plink/1.90beta
        plink --bfile {params.bfile}  --geno 0.1 --r2 dprime --ld-window 20 --ld-window-r2 0 --ld-window-kb 100 --out {params.out}
    """
rule plink_ld_by_chrom_01_out:
    input: expand(MPI_RELEASE_v31 + '/X1001tx_pairwise_ld1/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx_r2.ld',chrom=[str(c) for c in ATH_AUTOSOMES])

rule plink_ld_by_chrom_01_gz:# 1001 transcriptiome samples
    input: MPI_RELEASE_v31 + '/X1001tx_pairwise_ld1/1001genomes_snp-short-indel_only_ACGTN_{chrom}_1001tx_r2.ld'
    output: MPI_RELEASE_v31 + '/X1001tx_pairwise_ld2/chr{chrom}.tsv.gz'
    shell: """
        module load htslib/1.2.1
        tail -n +2 {input} | awk -v OFS='\\t' '{{ if ($7>=0.5) {{ print "chr"$1,$2,$3,$5,$6,$7,$8 }} }}' | bgzip > {output}
        tabix -b 2 -e 2 {output}
    """
rule plink_ld_by_chrom_01_gz_out:
    input: expand(MPI_RELEASE_v31 + '/X1001tx_pairwise_ld2/chr{chrom}.tsv.gz',chrom=[str(c) for c in ATH_AUTOSOMES])

rule plink_ld_by_chrom_02:# seed size accessions, biallelic SNPs only, MAF>=0.01, call rate >=90%
    input: tped=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4.tped'
    params: tfile=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4',out=MPI_RELEASE_v31 + '/seedsize_pairwise_ld1/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4_r2'
    output: MPI_RELEASE_v31 + '/seedsize_pairwise_ld1/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4_r2.ld'
    shell: """
        module load plink/1.90beta
        plink --tfile {params.tfile}  --geno 0.1 --r2 dprime --ld-window 20 --ld-window-r2 0 --ld-window-kb 100 --out {params.out}
    """
rule plink_ld_by_chrom_02_out:
    input: expand(MPI_RELEASE_v31+ '/seedsize_pairwise_ld1/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter4_r2.ld',seed_subset=['1kgen','1K'])
    

############
# GoShifter
############
rule goshifter_snpmap_01:
    input: traw=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1_recodeA.traw'
    output: txt=MPI_RELEASE_v31+'/X1001tx_filter1_goshifter_snpmap1/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.snpmappings.txt'
    shell: """
        echo -e 'SNP\\tChrom\\tBP' > {output}
        tail -n +2 {input} |  awk -v OFS='\\t' '{{ print $2,"chr"$1,$4 }}' >> {output}
    """

####################
# imputed SNP files
####################
    
##########
# SnpEff
##########
rule snpeff_dap_fam_01:
    input: vcf=MPI_RELEASE_v31 + '/1001genomes_snp-short-indel_only_ACGTN.vcf.gz',dap_peak=DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed',dap_fimo=DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed'
    params: snpeff_jar=SNPEFF_4p1L_JAR,snpeff_config=SNPEFF_4p1L_CONFIG,snpeff_options="-no-upstream -no-downstream -no-intron -no-intergenic -no-utr -i vcf -o vcf",stats=SNPEFF_DAP_01+'/{family_expSys}/1001genomes_snp-short-indel_only_ACGTN.vcf.html',csvStats=SNPEFF_DAP_01+'/{family_expSys}/1001genomes_snp-short-indel_only_ACGTN.vcf.csv',genome_version="athalianaTair10",D=SNPEFF_DAP_01+'/{family_expSys}',java_options="-Xmx8g"
    output: snpeff_vcf=SNPEFF_DAP_01+'/{family_expSys}/1001genomes_snp-short-indel_only_ACGTN.vcf'
    shell: """
        module load java/jdk-7u9-linux-x64
        mkdir -p {params.D}
        java {params.java_options} -jar {params.snpeff_jar} -c {params.snpeff_config} {params.snpeff_options} -stats {params.stats} -csvStats {params.csvStats} -interval {input.dap_peak} -interval {input.dap_fimo} {params.genome_version} {input.vcf} > {output.snpeff_vcf}
    """
rule snpeff_dap_fam_01_out:
    input: dap_fam_all_dep(SNPEFF_DAP_01+'/{family_expSys}/1001genomes_snp-short-indel_only_ACGTN.vcf',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')[2:3]

rule bedanno_dap_fam_01_mkdir:
    params: D=BEDANNO_DAP_01+'/{tg_ecotypeid}'
    output: BEDANNO_DAP_01+'/{tg_ecotypeid}/mkdir.log'
    shell:
        """
        mkdir -p {params.D}
        touch {output}
        """
rule bedanno_dap_fam_01_mkdir_out:
    input: acc_all(BEDANNO_DAP_01+'/{tg_ecotypeid}/mkdir.log',ecotype_list=ECOTYPE_LIST_G)

rule bedanno_dap_fam_01:# annotate by bedtools
    input: vcf=MPI_RELEASE_v31 + '/intersection_snp_short_indel_vcf/intersection_{tg_ecotypeid}.vcf.gz',dap_peak=dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed2',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed2'),dap_fimo=dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed2',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed2')
    params: D=BEDANNO_DAP_01+'/{tg_ecotypeid}',qname='gale.q',N="a_{tg_ecotypeid}"
    output: bedanno_samples=BEDANNO_DAP_01+'/{tg_ecotypeid}/intersection_snp_short_indel.samples',bedanno_bed=BEDANNO_DAP_01+'/{tg_ecotypeid}/intersection_snp_short_indel.bed'
    shell: """
        module load bedtools/2.25.0
        mkdir -p {params.D}
        echo {input.dap_peak} | sed 's/ /\\n/g' > {output.bedanno_samples}
        echo {input.dap_fimo} | sed 's/ /\\n/g' >> {output.bedanno_samples}
        bedtools annotate -counts -i {input.vcf} -files {input.dap_peak} {input.dap_fimo} > {output.bedanno_bed}
    """
rule bedanno_dap_fam_01_out:
    input: acc_all(BEDANNO_DAP_01+'/{tg_ecotypeid}/intersection_snp_short_indel.bed',ecotype_list=ECOTYPE_LIST_G)

rule bedanno_dap_fam_01_out_clean:
    shell: """rm {rules.bedanno_dap_fam_01_out.input}
    """
#snakemake -j 40 -k -p bedanno_dap_fam_01_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

bedanno_dap_fam_02_ref = {'allpos':MARGINAL_TEST_COV_09+'/geno_pos/position.bed',
                          'eqtl_1e5':MARGINAL_TEST_COV_09+'/results/lmm_master/chr1-5genes_1e5.bed.sorted.uniq',
                          'eqtl_1e6':MARGINAL_TEST_COV_09+'/results/lmm_master/chr1-5genes_1e6.bed.sorted.uniq',
                          'eqtl_1e7':MARGINAL_TEST_COV_09+'/results/lmm_master/chr1-5genes_1e7.bed.sorted.uniq',
                          'eqtl_1e8':MARGINAL_TEST_COV_09+'/results/lmm_master/chr1-5genes_1e8.bed.sorted.uniq',
                          'eqtl_1e9':MARGINAL_TEST_COV_09+'/results/lmm_master/chr1-5genes_1e9.bed.sorted.uniq',
                          'eqtl_1e10':MARGINAL_TEST_COV_09+'/results/lmm_master/chr1-5genes_1e10.bed.sorted.uniq'}
rule bedanno_dap_fam_02:# annotate by bedtools
    input: ref=lambda wildcards:bedanno_dap_fam_02_ref[wildcards.reftype],dap_peak=dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')[0],dap_fimo=dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')
    params: D=BEDANNO_DAP_02,qname='gale.q'
    output: bedanno_samples=BEDANNO_DAP_02+'/intersection_{reftype}.samples',bedanno_bed=BEDANNO_DAP_02+'/intersection_{reftype}.bed'
    shell: """
        module load bedtools/2.25.0
        mkdir -p {params.D}
        echo {input.dap_peak} | sed 's/ /\\n/g' > {output.bedanno_samples}
        echo {input.dap_fimo} | sed 's/ /\\n/g' >> {output.bedanno_samples}
        bedtools annotate -counts -i {input.ref} -files {input.dap_peak} {input.dap_fimo} > {output.bedanno_bed}
    """
rule bedanno_dap_fam_02_out:
    input: expand(BEDANNO_DAP_02+'/intersection_{reftype}.bed',reftype=['allpos','eqtl_1e5','eqtl_1e6','eqtl_1e7','eqtl_1e8','eqtl_1e9','eqtl_1e10'])

rule bedanno_dap_fam_02_out_clean:
    shell: """rm {rules.bedanno_dap_fam_02_out.input}
    """

###############
# splicing
###############
rule star_junctions_merge:
    input: junc_list=lambda wildcards: glob.glob(os.path.join(TFQ_GALE2,wildcards.tg_ecotypeid,"*.bam.junctions.bed"))
    params: N="sjuncm_{tg_ecotypeid}",L="2:00:00",qname="gale.q"
    log: TFQ_GALE2 + "/{tg_ecotypeid}/{tg_ecotypeid}.bam.junctions.merged.bed.log"
    output: TFQ_GALE2 + "/{tg_ecotypeid}/{tg_ecotypeid}.bam.junctions.merged.bed"
    shell: """
        module load shhuang bedops/2.4.15
        sort-bed {input.junc_list} > {output}
        """
rule star_junctions_merge_out:
    input: acc_all(TFQ_GALE2+"/{tg_ecotypeid}/{tg_ecotypeid}.bam.junctions.merged.bed",ecotype_list=ECOTYPE_LIST_SALK)
# snakemake -j 10 -k -p star_junctions_merge_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V -o {log}"


####################
# eQTL by LIMIX
####################
rule expr_file_chr1:
    input: txt=PROJ_RESULTS_PATH_GALE +'/{ath_norm_dir}/{norm_file}.txt'
    output: txt=PROJ_RESULTS_PATH_GALE +'/{ath_norm_dir}/{norm_file}_chr1.txt',csv=PROJ_RESULTS_PATH_GALE +'/{ath_norm_dir}/{norm_file}_chr1T.csv'
    shell: """
        module load ucsc_utils/20160210; printf 'id\\t' > {output.txt}; head -1 {input.txt} >> {output.txt}; awk -v OFS="\\t" '$1~/^AT1G/' {input.txt} >> {output.txt}; datamash transpose < {output.txt} | awk -v OFS="," 'NR==1 {{$1=""; print; next }} {{ $1=$1; print }}' > {output.csv}
    """
rule expr_file_chr1_out:
    input: [PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k%d_vst2_cv0p05_UQCounts_1001g_chr1T.csv'%k for k in [2,3,4,5]] + [PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k%d_vst2_cv0p05_RawCounts_1001g_chr1T.csv'%k for k in [2,3,4,5]] + [PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_normCounts_k%d_vst2_cv0p05_1001g_chr1T.csv'%k for k in [2,3,4,5]] + [PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k%d_vst2_cv0p05_UQRIN_1001g_chr1T.csv'%k for k in [2,3,4,5]]

rule limix_pconvert_01:
    input: csv=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_normCounts_k4_1001g_vst2_cv0p05_rinT.csv'
    output: hdf5=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_normCounts_k4_1001g_vst2_cv0p05_rinT.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """
rule limix_pconvert_02:
    input: csv=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_normCounts_k3_1001g_vst2_cv0p05_rinT.csv'
    output: hdf5=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_normCounts_k3_1001g_vst2_cv0p05_rinT.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """
rule limix_pconvert_03:
    input: csv=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01T_1001g.csv'
    output: hdf5=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01T_1001g.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """
rule limix_pconvert_04:
    input: csv=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_1001g_vst_cv0p05T.csv'
    output: hdf5=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_1001g_vst_cv0p05T.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """
rule limix_pconvert_05:
    input: csv=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_1001g_vst_cv0p05_rinT.csv'
    output: hdf5=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_1001g_vst_cv0p05_rinT.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """
rule limix_pconvert_06:
    input: csv=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k4_vst2_cv0p05_UQCounts_1001gT.csv'
    output: hdf5=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k4_vst2_cv0p05_UQCounts_1001gT.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """
rule limix_pconvert_07:# seed size
    input: csv=PROJ_DATA_PATH_GALE+'/seed_size/accx_size.csv'
    output: hdf5=PROJ_DATA_PATH_GALE+'/seed_size/accx_size.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """
rule limix_pconvert_08:
    input: csv=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_normCounts_k4_vst2T.csv'
    output: hdf5=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_normCounts_k4_vst2T.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """
rule limix_pconvert_09:
    input: csv=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k3_vst2_cv0p05_UQCounts_1001gT.csv'
    output: hdf5=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k3_vst2_cv0p05_UQCounts_1001gT.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """
rule limix_pconvert_10:
    input: csv=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.csv'
    output: hdf5=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """
rule limix_pconvert_11:
    input: csv=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.csv'
    output: hdf5=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """

rule limix_pconvert_chr1:
    input: csv=PROJ_RESULTS_PATH_GALE+'/{ath_norm_dir}/{norm_file}_chr1T.csv'
    output: hdf5=PROJ_RESULTS_PATH_GALE+'/{ath_norm_dir}/{norm_file}_chr1T.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --csv={input.csv}
    """
rule limix_pconvert_chr1_out:
    input: [os.path.splitext(p)[0]+'.hdf5' for p in rules.expr_file_chr1_out.input]

rule limix_gconvert_chr1:
    input: plink=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1_1001tx_filter1.bed'
    params: prefix=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1_1001tx_filter1'
    output: hdf5=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1_1001tx_filter1.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --plink={params.prefix}
    """
rule sample_cov_01:
    input: hdf5=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5'
    output: hdf5=MPI_RELEASE_v31+'/X1001tx_filter1/norm_cov_1001tx_filter1.hdf5',csv=MPI_RELEASE_v31+'/X1001tx_filter1/norm_cov_1001tx_filter1.csv'
    shell: """
        module load shhuang anaconda/2.4.1
        python {PROJ_PYTHON_PATH_GALE}/calc_geno_cov.py {input.hdf5} 1 {output.hdf5} {output.csv} 
    """
rule sample_cov_02:
    input: hdf5=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1.hdf5'
    output: hdf5=MPI_RELEASE_v31+'/filter1/norm_cov_filter1.hdf5',csv=MPI_RELEASE_v31+'/filter1/norm_cov_filter1.csv'
    shell: """
        module load shhuang anaconda/2.4.1
        python {PROJ_PYTHON_PATH_GALE}/calc_geno_cov.py {input.hdf5} 1 {output.hdf5} {output.csv} 
    """
    
rule limix_gconvert_01:
    input: plink=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1.bed'
    params: prefix=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1'
    output: hdf5=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_filter1.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --plink={params.prefix}
    """

rule limix_gconvert_03:
    input: plink=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3.bed'
    params: prefix=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3'
    output: hdf5=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        limix_converter --outfile={output.hdf5} --plink={params.prefix}
    """
rule limix_gconvert_03_out:
    input: expand(MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_{seed_subset}_filter3.hdf5',\
                  seed_subset=['1kgen'])
    
# maybe the wrong phenotype file
rule marginal_test_cov_05:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_normCounts_k4_1001g_vst_cv0p05_rinT.hdf5',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_W_k4.txt'
    params: D=MARGINAL_TEST_COV_05,qname='gale.q',N="assoc_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_05+'/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v5.py {input.geno_file} {input.pheno_file} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """

# 9751 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_normCounts_k4_1001g_vst_cv0p05_rin'
marginal_test_cov_05_slices = list(range(0,9751,100))+[9751]
rule marginal_test_cov_05_out:
    input: [MARGINAL_TEST_COV_05+'/logs/%04d-%04d.log'%(marginal_test_cov_05_slices[i],marginal_test_cov_05_slices[i+1]) for i in range(len(marginal_test_cov_05_slices)-1)]
                               
rule marginal_test_cov_05_out2:
    shell: "echo {rules.marginal_test_cov_05_out.input}"                               
#snakemake -j 10 -k -p marginal_test_cov_05_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

# still had 17 jobs 2016-02-05
rule marginal_test_cov_06:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_normCounts_k3_1001g_vst2_cv0p05_rinT.hdf5',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_W_k3.txt'
    params: D=MARGINAL_TEST_COV_06,qname='gale.q',N="assoc_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_06+'/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v5.py {input.geno_file} {input.pheno_file} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 9751 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_normCounts_k3_1001g_vst2_cv0p05_rinT'
marginal_test_cov_06_slices = list(range(0,9751,100))+[9751]
rule marginal_test_cov_06_out:
    input: [MARGINAL_TEST_COV_06+'/logs/%04d-%04d.log'%(marginal_test_cov_06_slices[i],marginal_test_cov_06_slices[i+1]) for i in range(len(marginal_test_cov_06_slices)-1)]
                               
rule marginal_test_cov_06_out2:
    shell: "echo {rules.marginal_test_cov_06_out.input}"                               
#snakemake -j 15 -k -p marginal_test_cov_06_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

rule marginal_test_cov_07:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_normCounts_k4_1001g_vst2_cv0p05_rinT.hdf5',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_W_k4.txt'
    params: D=MARGINAL_TEST_COV_07,qname='gale.q',N="assoc_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_07+'/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v5.py {input.geno_file} {input.pheno_file} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 9159 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_normCounts_k4_1001g_vst2_cv0p05_rinT'
marginal_test_cov_07_slices = list(range(0,9159,100))+[9159]
rule marginal_test_cov_07_out:
    input: [MARGINAL_TEST_COV_07+'/logs/%04d-%04d.log'%(marginal_test_cov_07_slices[i],marginal_test_cov_07_slices[i+1]) for i in range(len(marginal_test_cov_07_slices)-1)]
                               
rule marginal_test_cov_07_out2:
    shell: "echo {rules.marginal_test_cov_07_out.input}"                               
#snakemake -j 5 -k -p marginal_test_cov_07_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

# PANAMA normalization
# 2016-02-05: stopped at 70+/120 jobs, no genes have significant QTLs
rule marginal_test_cov_08:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_1001g_vst_cv0p05T.hdf5',panama_file=PROJ_RESULTS_PATH_GALE+'/calc_k_panama_2016-02-03/calc_k_panama_2016-02-03-_K20_dat.hdf5'
    params: D=MARGINAL_TEST_COV_08,norm_mode='RIN',qname='gale.q',N="assoc_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_08+'/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v7.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.panama_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 12268 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_1001g_vst_cv0p05T.hdf5
marginal_test_cov_08_slices = list(range(0,12268,100))+[12268]
rule marginal_test_cov_08_out:
    input: [MARGINAL_TEST_COV_08+'/logs/%04d-%04d.log'%(marginal_test_cov_08_slices[i],marginal_test_cov_08_slices[i+1]) for i in range(len(marginal_test_cov_08_slices)-1)]                             
rule marginal_test_cov_08_out2:
    shell: "echo {rules.marginal_test_cov_08_out.input}"                               
#snakemake -j 15 -k -p marginal_test_cov_08_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

#RUVg k=4, VST of raw counts, genes filtered by CV>0.05 on VST of raw counts, covariate uses RUVg covariate
rule marginal_test_cov_09:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_1001g_vst_cv0p05T.hdf5',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-gNorm_W_k4.txt'
    params: D=MARGINAL_TEST_COV_09,norm_mode='RIN',qname='gale.q',N="assoc_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_09+'/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v6.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 12268 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_1001g_vst_cv0p05T.hdf5
marginal_test_cov_09_slices = list(range(0,12268,100))+[12268]
rule marginal_test_cov_09_out:
    input: [MARGINAL_TEST_COV_09+'/logs/%04d-%04d.log'%(marginal_test_cov_09_slices[i],marginal_test_cov_09_slices[i+1]) for i in range(len(marginal_test_cov_09_slices)-1)]
rule marginal_test_cov_09_out2:
    shell: "echo {rules.marginal_test_cov_09_out.input}"                               
#snakemake -j 20 -k -p marginal_test_cov_09_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

rule marginal_test_cov_09_pos:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-01-03/ath1001_tx_norm_2016-01-03-filtered01_1001g_vst_cv0p05T.hdf5'
    params: D=MARGINAL_TEST_COV_09+'/geno_pos'
    output: MARGINAL_TEST_COV_09+'/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_09_pos2:
    input: pos_txt=MARGINAL_TEST_COV_09+'/geno_pos/position.txt'
    output: pos_bed=MARGINAL_TEST_COV_09+'/geno_pos/position.bed'
    shell: """
        tail -n +2 {input.pos_txt} | awk -v OFS='\\t' '{{ print "chr"$1,$2-1,$2 }}'  > {output.pos_bed}
    """

rule marginal_test_cov_09_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_09 + '/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_09+'/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_09_lmm_1e5_nr_out:
    input: lambda wildcards: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for p in glob.glob(os.path.join(MARGINAL_TEST_COV_09,'results/lmm_pval',"*.txt"))]

rule marginal_test_cov_09_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_09+'/geno_pos/position.txt',nr=MARGINAL_TEST_COV_09 + '/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_09 + '/results/lmm_1e5_nr'
    output: bed=MARGINAL_TEST_COV_09+'/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        cd {params.D}; awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_09_lmm_1e5_bed_out:
    input: lambda wildcards: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for p in glob.glob(os.path.join(MARGINAL_TEST_COV_09,'results/lmm_1e5_nr',"*.txt"))]

marginal_test_cov_09_lmm_thresh = {'1e5':'1e-5','1e6':'1e-6','1e7':'1e-7','1e8':'1e-8','1e9':'1e-9','1e10':'1e-10'}
rule marginal_test_cov_09_lmm_master:
    input: glob.glob(os.path.join(MARGINAL_TEST_COV_09,'results/lmm_1e5_bed',"*.bed"))
    params: D=os.path.join(MARGINAL_TEST_COV_09,'results/lmm_1e5_bed'),thr=lambda wildcards:marginal_test_cov_09_lmm_thresh[wildcards.thr]
    output: MARGINAL_TEST_COV_09+'/results/lmm_master/allgenes_{thr}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<={params.thr}' > {output}
    """
rule marginal_test_cov_09_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_09+'/results/lmm_master/allgenes_{thr}.bed',thr=list(marginal_test_cov_09_lmm_thresh.keys()))

rule marginal_test_cov_09_lmm_master_noCM:
    input: MARGINAL_TEST_COV_09+'/results/lmm_master/allgenes_{thr}.bed'
    output: MARGINAL_TEST_COV_09+'/results/lmm_master/chr1-5genes_{thr}.bed'
    shell: """
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {input} > {output}
    """
rule marginal_test_cov_09_lmm_master_noCM_out:
    input: expand(MARGINAL_TEST_COV_09+'/results/lmm_master/chr1-5genes_{thr}.bed',thr=list(marginal_test_cov_09_lmm_thresh.keys()))
rule marginal_test_cov_09_lmm_master_sorted:
    input: MARGINAL_TEST_COV_09+'/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_09+'/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_09_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_09+'/results/lmm_master/{geneset}_{thr}.bed.sorted',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_09_lmm_thresh.keys())
rule marginal_test_cov_09_lmm_master_uniq:
    input: MARGINAL_TEST_COV_09+'/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_09+'/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_09_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_09+'/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_09_lmm_thresh.keys())

#RUVg k=4, UQ normalized counts, genes filtered by CV>0.05 on RUVg normalized, VST counts, covariate uses RUVg covariate, RIN transformed in eQTL
rule marginal_test_cov_10:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k4_vst2_cv0p05_UQCounts_1001gT.hdf5',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_W_k4.txt'
    params: D=MARGINAL_TEST_COV_10,norm_mode='RIN',qname='gale.q',N="assoc_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_10+'/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v6.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 8366 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k4_vst2_cv0p05_UQCounts_1001gT.hdf5'
marginal_test_cov_10_slices = list(range(0,8366,100))+[8366]
rule marginal_test_cov_10_out:
    input: [MARGINAL_TEST_COV_10+'/logs/%04d-%04d.log'%(marginal_test_cov_10_slices[i],marginal_test_cov_10_slices[i+1]) for i in range(len(marginal_test_cov_10_slices)-1)]
rule marginal_test_cov_10_out2:
    shell: "echo {rules.marginal_test_cov_10_out.input}"                               
#snakemake -j 20 -k -p marginal_test_cov_10_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

rule marginal_test_cov_10_pos:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k4_vst2_cv0p05_UQCounts_1001gT.hdf5'
    params: D=MARGINAL_TEST_COV_10+'/geno_pos'
    output: MARGINAL_TEST_COV_10+'/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_10_pos2:
    input: pos_txt=MARGINAL_TEST_COV_10+'/geno_pos/position.txt'
    output: pos_bed=MARGINAL_TEST_COV_10+'/geno_pos/position.bed'
    shell: """
        tail -n +2 {input.pos_txt} | awk -v OFS='\\t' '{{ print "chr"$1,$2-1,$2 }}'  > {output.pos_bed}
    """
rule marginal_test_cov_10_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_10 + '/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_10+'/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_10_lmm_1e5_nr_out:
    input: lambda wildcards: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for p in glob.glob(os.path.join(MARGINAL_TEST_COV_10,'results/lmm_pval',"*.txt"))]

rule marginal_test_cov_10_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_10+'/geno_pos/position.txt',nr=MARGINAL_TEST_COV_10 + '/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_10 + '/results/lmm_1e5_nr'
    output: bed=MARGINAL_TEST_COV_10+'/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        cd {params.D}; awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_10_lmm_1e5_bed_out:
    input: lambda wildcards: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for p in glob.glob(os.path.join(MARGINAL_TEST_COV_10,'results/lmm_1e5_nr',"*.txt"))]

marginal_test_cov_10_lmm_thresh = {'1e-5':'1e-5','1e-6':'1e-6','1e-7':'1e-7','1e-8':'1e-8','1e-9':'1e-9','1e-10':'1e-10'}
rule marginal_test_cov_10_lmm_master:
    input: glob.glob(os.path.join(MARGINAL_TEST_COV_10,'results/lmm_1e5_bed',"*.bed"))
    params: D=os.path.join(MARGINAL_TEST_COV_10,'results/lmm_1e5_bed')
    output: MARGINAL_TEST_COV_10+'/results/lmm_master/allgenes_{thr}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<={wildcards.thr}' > {output}
    """
rule marginal_test_cov_10_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_10+'/results/lmm_master/allgenes_{thr}.bed',thr=list(marginal_test_cov_10_lmm_thresh.keys()))

rule marginal_test_cov_10_lmm_master_noCM:
    input: MARGINAL_TEST_COV_10+'/results/lmm_master/allgenes_{thr}.bed'
    output: MARGINAL_TEST_COV_10+'/results/lmm_master/chr1-5genes_{thr}.bed'
    shell: """
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {input} > {output}
    """
rule marginal_test_cov_10_lmm_master_noCM_out:
    input: expand(MARGINAL_TEST_COV_10+'/results/lmm_master/chr1-5genes_{thr}.bed',thr=list(marginal_test_cov_10_lmm_thresh.keys()))
rule marginal_test_cov_10_lmm_master_sorted:
    input: MARGINAL_TEST_COV_10+'/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_10+'/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_10_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_10+'/results/lmm_master/{geneset}_{thr}.bed.sorted',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_10_lmm_thresh.keys())
rule marginal_test_cov_10_lmm_master_uniq:
    input: MARGINAL_TEST_COV_10+'/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_10+'/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_10_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_10+'/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_10_lmm_thresh.keys())

#RUVg k=2,3,4,5, UQ normalized counts, chromosome 1 genes filtered by CV>0.05 on RUVg normalized, VST counts, covariate uses RUVg covariate
marginal_test_cov_11_ngenes = {'2':2576,'3':2423,'4':2166,'5':2022}
marginal_test_cov_11_setup = dict([(k,{'pheno_file':PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k%s_vst2_cv0p05_UQCounts_1001g_chr1T.hdf5'%k,'cov_file':PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_W_k%s.txt'%k}) for k in ['2','3','4','5']])
rule marginal_test_cov_11_chr1:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1_1001tx_filter1.hdf5',pheno_file=lambda wildcards:marginal_test_cov_11_setup[wildcards.k]['pheno_file'],cov_file=lambda wildcards:marginal_test_cov_11_setup[wildcards.k]['cov_file']
    params: D=MARGINAL_TEST_COV_11+'_k{k}',norm_mode='RIN',qname='gale.q',N="k{k}_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_11+"_k{k}/logs/{RNA_start}-{RNA_end}.log"
    shell: """
       module load shhuang anaconda/2.4.1
       mkdir -p {params.D}
       python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v6.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
rule marginal_test_cov_11_chr1_out:
    input: [MARGINAL_TEST_COV_11+'_k%s/logs/%04d-%04d.log'%(k,0,marginal_test_cov_11_ngenes[k]) for k in ['2','3','4','5']]
#snakemake -j 20 -k -p marginal_test_cov_11_chr1_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"
rule marginal_test_cov_11_chr1_pos:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1_1001tx_filter1.hdf5',pheno_file=lambda wildcards:marginal_test_cov_11_setup[wildcards.k]['pheno_file']
    params: D=MARGINAL_TEST_COV_11+'_k{k}/geno_pos'
    output: MARGINAL_TEST_COV_11+'_k{k}/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_11_chr1_pos_out:
    input: [MARGINAL_TEST_COV_11+'_k%s/geno_pos/get_geno_pos.log'%k for k in ['2','3','4','5']]
rule marginal_test_cov_11_chr1_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_11+'_k{k}/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_11+'_k{k}/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_11_chr1_lmm_1e5_nr_out:
    input: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for k in ['2','3','4','5'] for p in glob.glob(os.path.join(MARGINAL_TEST_COV_11+'_k%s/results/lmm_pval'%k,"*.txt"))]
rule marginal_test_cov_11_chr1_lmm_1e5_nr_out2:
    shell: "echo {rules.marginal_test_cov_11_chr1_lmm_1e5_nr_out.input}"

marginal_test_cov_11_lmm_thresh = ['1e-5','1e-6','1e-7','1e-8','5e-8','1e-9','1e-10']
rule marginal_test_cov_11_chr1_lmm_thr_counts:
    input:  MARGINAL_TEST_COV_11+'_k{k}/results/lmm_pval/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_11+'_k{k}/results/lmm_counts_p{thr}'
    output: MARGINAL_TEST_COV_11+'_k{k}/results/lmm_counts_p{thr}/{at_id}.txt'
    shell: """
        mkdir -p {params.D}; awk -v OFS='\\t' 'BEGIN{{ count=0 }} $1<={wildcards.thr} {{ count++ }} END {{ print count }}' {input} > {output}
    """
rule marginal_test_cov_11_chr1_lmm_thr_counts_out:
    input: [p.replace('/lmm_pval/','/lmm_counts_p%s/'%thr) for thr in marginal_test_cov_11_lmm_thresh for k in ['2','3','4','4','5'] for p in glob.glob(os.path.join(MARGINAL_TEST_COV_11+'_k%s/results/lmm_pval'%k,"*.txt"))]

rule marginal_test_cov_11_chr1_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_11+'_k{k}/geno_pos/position.txt',nr=MARGINAL_TEST_COV_11 + '_k{k}/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_11 + '_k{k}/results/lmm_1e5_nr'
    output: bed=MARGINAL_TEST_COV_11+'_k{k}/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        cd {params.D}; awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_11_chr1_lmm_bed_out:
    input: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for k in ['2','3','4','5'] for p in glob.glob(MARGINAL_TEST_COV_11+'_k%s/results/lmm_1e5_nr/*.txt'%k)]

rule marginal_test_cov_11_chr1_lmm_master:
    input: glob.glob(MARGINAL_TEST_COV_11+'_k{k}/results/lmm_1e5_bed/*.bed')
    params: D=MARGINAL_TEST_COV_11+'_k{k}/results/lmm_1e5_bed'
    output: MARGINAL_TEST_COV_11+'_k{k}/results/lmm_master/chr1genes_1e5.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<=1e-5' > {output}
    """
rule marginal_test_cov_11_chr1_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_11+'_k{k}/results/lmm_master/chr1genes_1e5.bed',k=['2','3','4','5'])
rule marginal_test_cov_11_chr1_lmm_master_sorted:
    input: MARGINAL_TEST_COV_11+'_k{k}/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_11+'_k{k}/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_11_chr1_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_11+'_k{k}/results/lmm_master/{geneset}_1e5.bed.sorted',geneset=['chr1genes'],k=['2','3','4','5'])
rule marginal_test_cov_11_chr1_lmm_master_uniq:
    input: MARGINAL_TEST_COV_11+'_k{k}/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_11+'_k{k}/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_11_chr1_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_11+'_k{k}/results/lmm_master/{geneset}_1e5.bed.sorted.uniq',geneset=['chr1genes'],k=['2','3','4','5'])


#RUVg k=3, UQ normalized counts, genes filtered by CV>0.05 on RUVg normalized, VST counts, covariate uses RUVg covariate
rule marginal_test_cov_12:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k3_vst2_cv0p05_UQCounts_1001gT.hdf5',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_W_k3.txt'
    params: D=MARGINAL_TEST_COV_12,norm_mode='RIN',qname='gale.q',N="assoc_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_12+'/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v6.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 9370 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k3_vst2_cv0p05_UQCounts_1001gT.hdf5'
marginal_test_cov_12_slices = list(range(0,9370,100))+[9370]
rule marginal_test_cov_12_out:
    input: [MARGINAL_TEST_COV_12+'/logs/%04d-%04d.log'%(marginal_test_cov_12_slices[i],marginal_test_cov_12_slices[i+1]) for i in range(len(marginal_test_cov_12_slices)-1)]
rule marginal_test_cov_12_out2:
    shell: "echo {rules.marginal_test_cov_12_out.input}"     
#snakemake -j 20 -k -p marginal_test_cov_12_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"
rule marginal_test_cov_12_pos:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k3_vst2_cv0p05_UQCounts_1001gT.hdf5'
    params: D=MARGINAL_TEST_COV_12+'/geno_pos'
    output: MARGINAL_TEST_COV_12+'/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_12_pos2:
    input: pos_txt=MARGINAL_TEST_COV_12+'/geno_pos/position.txt'
    output: pos_bed=MARGINAL_TEST_COV_12+'/geno_pos/position.bed'
    shell: """
        tail -n +2 {input.pos_txt} | awk -v OFS='\\t' '{{ print "chr"$1,$2-1,$2 }}'  > {output.pos_bed}
    """
rule marginal_test_cov_12_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_12 + '/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_12+'/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_12_lmm_1e5_nr_out:
    input: lambda wildcards: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for p in glob.glob(os.path.join(MARGINAL_TEST_COV_12,'results/lmm_pval',"*.txt"))]

rule marginal_test_cov_12_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_12+'/geno_pos/position.txt',nr=MARGINAL_TEST_COV_12 + '/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_12 + '/results/lmm_1e5_nr'
    output: bed=MARGINAL_TEST_COV_12+'/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        cd {params.D}; awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_12_lmm_1e5_bed_out:
    input: lambda wildcards: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for p in glob.glob(os.path.join(MARGINAL_TEST_COV_12,'results/lmm_1e5_nr',"*.txt"))]
    
marginal_test_cov_12_lmm_thresh = {'1e-5':'1e-5','1e-6':'1e-6','1e-7':'1e-7','1e-8':'1e-8','1e-9':'1e-9','1e-10':'1e-10'}
rule marginal_test_cov_12_lmm_master:
    input: glob.glob(os.path.join(MARGINAL_TEST_COV_12,'results/lmm_1e5_bed',"*.bed"))
    params: D=os.path.join(MARGINAL_TEST_COV_12,'results/lmm_1e5_bed')
    output: MARGINAL_TEST_COV_12+'/results/lmm_master/allgenes_{thr}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<={wildcards.thr}' > {output}
    """
rule marginal_test_cov_12_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_12+'/results/lmm_master/allgenes_{thr}.bed',thr=list(marginal_test_cov_12_lmm_thresh.keys()))

rule marginal_test_cov_12_lmm_master_noCM:
    input: MARGINAL_TEST_COV_12+'/results/lmm_master/allgenes_{thr}.bed'
    output: MARGINAL_TEST_COV_12+'/results/lmm_master/chr1-5genes_{thr}.bed'
    shell: """
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {input} > {output}
    """
rule marginal_test_cov_12_lmm_master_noCM_out:
    input: expand(MARGINAL_TEST_COV_12+'/results/lmm_master/chr1-5genes_{thr}.bed',thr=list(marginal_test_cov_12_lmm_thresh.keys()))
rule marginal_test_cov_12_lmm_master_sorted:
    input: MARGINAL_TEST_COV_12+'/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_12+'/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_12_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_12+'/results/lmm_master/{geneset}_{thr}.bed.sorted',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_12_lmm_thresh.keys())
rule marginal_test_cov_12_lmm_master_uniq:
    input: MARGINAL_TEST_COV_12+'/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_12+'/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_12_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_12+'/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_12_lmm_thresh.keys())

#RUVg k=2,3,4,5, raw counts, chromosome 1 genes filtered by CV>0.05 on UQ RUVg normalized, VST counts, covariate uses RUVg covariate
marginal_test_cov_13_ngenes = {'2':2576,'3':2423,'4':2166,'5':2022}
marginal_test_cov_13_setup = dict([(k,{'pheno_file':PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_k%s_vst2_cv0p05_RawCounts_1001g_chr1T.hdf5'%k,'cov_file':PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_W_k%s.txt'%k}) for k in ['2','3','4','5']])
rule marginal_test_cov_13_chr1:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1_1001tx_filter1.hdf5',pheno_file=lambda wildcards:marginal_test_cov_13_setup[wildcards.k]['pheno_file'],cov_file=lambda wildcards:marginal_test_cov_13_setup[wildcards.k]['cov_file']
    params: D=MARGINAL_TEST_COV_13+'_k{k}',norm_mode='RIN',qname='gale.q',N="k{k}_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_13+"_k{k}/logs/{RNA_start}-{RNA_end}.log"
    shell: """
       module load shhuang anaconda/2.4.1
       mkdir -p {params.D}
       python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v6.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
rule marginal_test_cov_13_chr1_out:
    input: [MARGINAL_TEST_COV_13+'_k%s/logs/%04d-%04d.log'%(k,0,marginal_test_cov_13_ngenes[k]) for k in ['2','3','4','5']]
#snakemake -j 20 -k -p marginal_test_cov_13_chr1_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"
rule marginal_test_cov_13_chr1_pos:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1_1001tx_filter1.hdf5',pheno_file=lambda wildcards:marginal_test_cov_13_setup[wildcards.k]['pheno_file']
    params: D=MARGINAL_TEST_COV_13+'_k{k}/geno_pos'
    output: MARGINAL_TEST_COV_13+'_k{k}/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_13_chr1_pos_out:
    input: [MARGINAL_TEST_COV_13+'_k%s/geno_pos/get_geno_pos.log'%k for k in ['2','3','4','5']]
rule marginal_test_cov_13_chr1_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_13+'_k{k}/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_13+'_k{k}/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_13_chr1_lmm_1e5_nr_out:
    input: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for k in ['2','3','4','5'] for p in glob.glob(os.path.join(MARGINAL_TEST_COV_13+'_k%s/results/lmm_pval'%k,"*.txt"))]
rule marginal_test_cov_13_chr1_lmm_1e5_nr_out2:
    shell: "echo {rules.marginal_test_cov_13_chr1_lmm_1e5_nr_out.input}"

marginal_test_cov_13_lmm_thresh = ['1e-5','1e-6','1e-7','1e-8','5e-8','1e-9','1e-10']
rule marginal_test_cov_13_chr1_lmm_thr_counts:
    input:  MARGINAL_TEST_COV_13+'_k{k}/results/lmm_pval/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_13+'_k{k}/results/lmm_counts_p{thr}'
    output: MARGINAL_TEST_COV_13+'_k{k}/results/lmm_counts_p{thr}/{at_id}.txt'
    shell: """
        mkdir -p {params.D}; awk -v OFS='\\t' 'BEGIN{{ count=0 }} $1<={wildcards.thr} {{ count++ }} END {{ print count }}' {input} > {output}
    """
rule marginal_test_cov_13_chr1_lmm_thr_counts_out:
    input: [p.replace('/lmm_pval/','/lmm_counts_p%s/'%thr) for thr in marginal_test_cov_13_lmm_thresh for k in ['2','3','4','4','5'] for p in glob.glob(os.path.join(MARGINAL_TEST_COV_13+'_k%s/results/lmm_pval'%k,"*.txt"))]

rule marginal_test_cov_13_chr1_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_13+'_k{k}/geno_pos/position.txt',nr=MARGINAL_TEST_COV_13 + '_k{k}/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_13 + '_k{k}/results/lmm_1e5_nr'
    output: bed=MARGINAL_TEST_COV_13+'_k{k}/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        cd {params.D}; awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_13_chr1_lmm_bed_out:
    input: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for k in ['2','3','4','5'] for p in glob.glob(MARGINAL_TEST_COV_13+'_k%s/results/lmm_1e5_nr/*.txt'%k)]

rule marginal_test_cov_13_chr1_lmm_master:
    input: glob.glob(MARGINAL_TEST_COV_13+'_k{k}/results/lmm_1e5_bed/*.bed')
    params: D=MARGINAL_TEST_COV_13+'_k{k}/results/lmm_1e5_bed'
    output: MARGINAL_TEST_COV_13+'_k{k}/results/lmm_master/chr1genes_1e5.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<=1e-5' > {output}
    """
rule marginal_test_cov_13_chr1_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_13+'_k{k}/results/lmm_master/chr1genes_1e5.bed',k=['2','3','4','5'])
rule marginal_test_cov_13_chr1_lmm_master_sorted:
    input: MARGINAL_TEST_COV_13+'_k{k}/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_13+'_k{k}/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_13_chr1_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_13+'_k{k}/results/lmm_master/{geneset}_1e5.bed.sorted',geneset=['chr1genes'],k=['2','3','4','5'])
rule marginal_test_cov_13_chr1_lmm_master_uniq:
    input: MARGINAL_TEST_COV_13+'_k{k}/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_13+'_k{k}/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_13_chr1_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_13+'_k{k}/results/lmm_master/{geneset}_1e5.bed.sorted.uniq',geneset=['chr1genes'],k=['2','3','4','5'])


#RUVg k=2,3,4,5, UQ RUVg normalized VST counts, chromosome 1 genes filtered by CV>0.05 on UQ RUVg normalized, VST counts, no covariate
marginal_test_cov_14_ngenes = {'2':2576,'3':2423,'4':2166,'5':2022}
marginal_test_cov_14_setup = dict([(k,{'pheno_file':PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQ_gNorm_normCounts_k%s_vst2_cv0p05_1001g_chr1T.hdf5'%k,'cov_file':'None'}) for k in ['2','3','4','5']])
rule marginal_test_cov_14_chr1:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1_1001tx_filter1.hdf5',pheno_file=lambda wildcards:marginal_test_cov_14_setup[wildcards.k]['pheno_file']
    params: D=MARGINAL_TEST_COV_14+'_k{k}',norm_mode='RIN',qname='gale.q',N="k{k}_{RNA_start}-{RNA_end}",cov_file=lambda wildcards:marginal_test_cov_14_setup[wildcards.k]['cov_file']
    output: MARGINAL_TEST_COV_14+"_k{k}/logs/{RNA_start}-{RNA_end}.log"
    shell: """
       module load shhuang anaconda/2.4.1
       mkdir -p {params.D}
       python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v6.py {input.geno_file} {input.pheno_file} {params.norm_mode} {params.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
rule marginal_test_cov_14_chr1_out:
    input: [MARGINAL_TEST_COV_14+'_k%s/logs/%04d-%04d.log'%(k,0,marginal_test_cov_14_ngenes[k]) for k in ['2','3','4','5']]
#snakemake -j 20 -k -p marginal_test_cov_14_chr1_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"
rule marginal_test_cov_14_chr1_pos:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1_1001tx_filter1.hdf5',pheno_file=lambda wildcards:marginal_test_cov_14_setup[wildcards.k]['pheno_file']
    params: D=MARGINAL_TEST_COV_14+'_k{k}/geno_pos'
    output: MARGINAL_TEST_COV_14+'_k{k}/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_14_chr1_pos_out:
    input: [MARGINAL_TEST_COV_14+'_k%s/geno_pos/get_geno_pos.log'%k for k in ['2','3','4','5']]
rule marginal_test_cov_14_chr1_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_14+'_k{k}/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_14+'_k{k}/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_14_chr1_lmm_1e5_nr_out:
    input: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for k in ['2','3','4','5'] for p in glob.glob(os.path.join(MARGINAL_TEST_COV_14+'_k%s/results/lmm_pval'%k,"*.txt"))]
rule marginal_test_cov_14_chr1_lmm_1e5_nr_out2:
    shell: "echo {rules.marginal_test_cov_14_chr1_lmm_1e5_nr_out.input}"

marginal_test_cov_14_lmm_thresh = ['1e-5','1e-6','1e-7','1e-8','5e-8','1e-9','1e-10']
rule marginal_test_cov_14_chr1_lmm_thr_counts:
    input:  MARGINAL_TEST_COV_14+'_k{k}/results/lmm_pval/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_14+'_k{k}/results/lmm_counts_p{thr}'
    output: MARGINAL_TEST_COV_14+'_k{k}/results/lmm_counts_p{thr}/{at_id}.txt'
    shell: """
        mkdir -p {params.D}; awk -v OFS='\\t' 'BEGIN{{ count=0 }} $1<={wildcards.thr} {{ count++ }} END {{ print count }}' {input} > {output}
    """
rule marginal_test_cov_14_chr1_lmm_thr_counts_out:
    input: [p.replace('/lmm_pval/','/lmm_counts_p%s/'%thr) for thr in marginal_test_cov_14_lmm_thresh for k in ['2','3','4','4','5'] for p in glob.glob(os.path.join(MARGINAL_TEST_COV_14+'_k%s/results/lmm_pval'%k,"*.txt"))]

rule marginal_test_cov_14_chr1_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_14+'_k{k}/geno_pos/position.txt',nr=MARGINAL_TEST_COV_14 + '_k{k}/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_14 + '_k{k}/results/lmm_1e5_nr'
    output: bed=MARGINAL_TEST_COV_14+'_k{k}/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        cd {params.D}; awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_14_chr1_lmm_bed_out:
    input: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for k in ['2','3','4','5'] for p in glob.glob(MARGINAL_TEST_COV_14+'_k%s/results/lmm_1e5_nr/*.txt'%k)]

rule marginal_test_cov_14_chr1_lmm_master:
    input: glob.glob(MARGINAL_TEST_COV_14+'_k{k}/results/lmm_1e5_bed/*.bed')
    params: D=MARGINAL_TEST_COV_14+'_k{k}/results/lmm_1e5_bed'
    output: MARGINAL_TEST_COV_14+'_k{k}/results/lmm_master/chr1genes_1e5.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<=1e-5' > {output}
    """
rule marginal_test_cov_14_chr1_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_14+'_k{k}/results/lmm_master/chr1genes_1e5.bed',k=['2','3','4','5'])
rule marginal_test_cov_14_chr1_lmm_master_sorted:
    input: MARGINAL_TEST_COV_14+'_k{k}/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_14+'_k{k}/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_14_chr1_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_14+'_k{k}/results/lmm_master/{geneset}_1e5.bed.sorted',geneset=['chr1genes'],k=['2','3','4','5'])
rule marginal_test_cov_14_chr1_lmm_master_uniq:
    input: MARGINAL_TEST_COV_14+'_k{k}/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_14+'_k{k}/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_14_chr1_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_14+'_k{k}/results/lmm_master/{geneset}_1e5.bed.sorted.uniq',geneset=['chr1genes'],k=['2','3','4','5'])

#RUVg k=2,3,4,5, UQ RIN normalized, chromosome 1 genes filtered by CV>0.05 on UQ RUVg normalized, VST counts, covariate use UQ RIN RUVg covariates
marginal_test_cov_15_ngenes = {'2':2576,'3':2423,'4':2166,'5':2022}
marginal_test_cov_15_setup = dict([(k,{'pheno_file':PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k%s_vst2_cv0p05_UQRIN_1001g_chr1T.hdf5'%k,'cov_file':PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_W_k%s.txt'%k}) for k in ['2','3','4','5']])
rule marginal_test_cov_15_chr1:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1_1001tx_filter1.hdf5',pheno_file=lambda wildcards:marginal_test_cov_15_setup[wildcards.k]['pheno_file']
    params: D=MARGINAL_TEST_COV_15+'_k{k}',norm_mode='None',qname='gale.q',N="k{k}_{RNA_start}-{RNA_end}",cov_file=lambda wildcards:marginal_test_cov_15_setup[wildcards.k]['cov_file']
    output: MARGINAL_TEST_COV_15+"_k{k}/logs/{RNA_start}-{RNA_end}.log"
    shell: """
       module load shhuang anaconda/2.4.1
       mkdir -p {params.D}
       python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v6.py {input.geno_file} {input.pheno_file} {params.norm_mode} {params.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
rule marginal_test_cov_15_chr1_out:
    input: [MARGINAL_TEST_COV_15+'_k%s/logs/%04d-%04d.log'%(k,0,marginal_test_cov_15_ngenes[k]) for k in ['2','3','4','5']]
#snakemake -j 20 -k -p marginal_test_cov_15_chr1_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"
rule marginal_test_cov_15_chr1_pos:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1_1001tx_filter1.hdf5',pheno_file=lambda wildcards:marginal_test_cov_15_setup[wildcards.k]['pheno_file']
    params: D=MARGINAL_TEST_COV_15+'_k{k}/geno_pos'
    output: MARGINAL_TEST_COV_15+'_k{k}/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_15_chr1_pos_out:
    input: [MARGINAL_TEST_COV_15+'_k%s/geno_pos/get_geno_pos.log'%k for k in ['2','3','4','5']]
rule marginal_test_cov_15_chr1_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_15+'_k{k}/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_15+'_k{k}/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_15_chr1_lmm_1e5_nr_out:
    input: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for k in ['2','3','4','5'] for p in glob.glob(os.path.join(MARGINAL_TEST_COV_15+'_k%s/results/lmm_pval'%k,"*.txt"))]
rule marginal_test_cov_15_chr1_lmm_1e5_nr_out2:
    shell: "echo {rules.marginal_test_cov_15_chr1_lmm_1e5_nr_out.input}"

marginal_test_cov_15_lmm_thresh = ['1e-5','1e-6','1e-7','1e-8','5e-8','1e-9','1e-10']
rule marginal_test_cov_15_chr1_lmm_thr_counts:
    input:  MARGINAL_TEST_COV_15+'_k{k}/results/lmm_pval/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_15+'_k{k}/results/lmm_counts_p{thr}'
    output: MARGINAL_TEST_COV_15+'_k{k}/results/lmm_counts_p{thr}/{at_id}.txt'
    shell: """
        mkdir -p {params.D}; awk -v OFS='\\t' 'BEGIN{{ count=0 }} $1<={wildcards.thr} {{ count++ }} END {{ print count }}' {input} > {output}
    """
rule marginal_test_cov_15_chr1_lmm_thr_counts_out:
    input: [p.replace('/lmm_pval/','/lmm_counts_p%s/'%thr) for thr in marginal_test_cov_15_lmm_thresh for k in ['2','3','4','4','5'] for p in glob.glob(os.path.join(MARGINAL_TEST_COV_15+'_k%s/results/lmm_pval'%k,"*.txt"))]

rule marginal_test_cov_15_chr1_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_15+'_k{k}/geno_pos/position.txt',nr=MARGINAL_TEST_COV_15 + '_k{k}/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_15 + '_k{k}/results/lmm_1e5_nr'
    output: bed=MARGINAL_TEST_COV_15+'_k{k}/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        cd {params.D}; awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_15_chr1_lmm_bed_out:
    input: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for k in ['2','3','4','5'] for p in glob.glob(MARGINAL_TEST_COV_15+'_k%s/results/lmm_1e5_nr/*.txt'%k)]

rule marginal_test_cov_15_chr1_lmm_master:
    input: glob.glob(MARGINAL_TEST_COV_15+'_k{k}/results/lmm_1e5_bed/*.bed')
    params: D=MARGINAL_TEST_COV_15+'_k{k}/results/lmm_1e5_bed'
    output: MARGINAL_TEST_COV_15+'_k{k}/results/lmm_master/chr1genes_1e5.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<=1e-5' > {output}
    """
rule marginal_test_cov_15_chr1_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_15+'_k{k}/results/lmm_master/chr1genes_1e5.bed',k=['2','3','4','5'])
rule marginal_test_cov_15_chr1_lmm_master_sorted:
    input: MARGINAL_TEST_COV_15+'_k{k}/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_15+'_k{k}/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_15_chr1_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_15+'_k{k}/results/lmm_master/{geneset}_1e5.bed.sorted',geneset=['chr1genes'],k=['2','3','4','5'])
rule marginal_test_cov_15_chr1_lmm_master_uniq:
    input: MARGINAL_TEST_COV_15+'_k{k}/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_15+'_k{k}/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_15_chr1_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_15+'_k{k}/results/lmm_master/{geneset}_1e5.bed.sorted.uniq',geneset=['chr1genes'],k=['2','3','4','5'])


#RUVg k=4, UQ RIN normalized, genes filtered by CV>0.05 on UQ RUVg normalized, VST counts, covariate use UQ RIN RUVg covariates
rule marginal_test_cov_16:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_W_k4.txt'
    params: D=MARGINAL_TEST_COV_16,norm_mode='None',qname='gale.q',N="assoc_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_16+'/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v6.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 8366 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
marginal_test_cov_16_ngenes = 8366
marginal_test_cov_16_slices = list(range(0,marginal_test_cov_16_ngenes,100))+[marginal_test_cov_16_ngenes]
rule marginal_test_cov_16_out:
    input: [MARGINAL_TEST_COV_16+'/logs/%04d-%04d.log'%(marginal_test_cov_16_slices[i],marginal_test_cov_16_slices[i+1]) for i in range(len(marginal_test_cov_16_slices)-1)]
rule marginal_test_cov_16_out2:
    shell: "echo {rules.marginal_test_cov_16_out.input}"     
#snakemake -j 20 -k -p marginal_test_cov_16_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"
rule marginal_test_cov_16_pos:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
    params: D=MARGINAL_TEST_COV_16+'/geno_pos'
    output: MARGINAL_TEST_COV_16+'/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_16_pos2:
    input: pos_txt=MARGINAL_TEST_COV_16+'/geno_pos/position.txt'
    output: pos_bed=MARGINAL_TEST_COV_16+'/geno_pos/position.bed'
    shell: """
        tail -n +2 {input.pos_txt} | awk -v OFS='\\t' '{{ print "chr"$1,$2-1,$2 }}'  > {output.pos_bed}
    """
rule marginal_test_cov_16_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_16 + '/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_16+'/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_16_lmm_1e5_nr_out:
    input: lambda wildcards: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for p in glob.glob(os.path.join(MARGINAL_TEST_COV_16,'results/lmm_pval',"*.txt"))]

rule marginal_test_cov_16_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_16+'/geno_pos/position.txt',nr=MARGINAL_TEST_COV_16 + '/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_16 + '/results/lmm_1e5_nr'
    output: bed=MARGINAL_TEST_COV_16+'/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        cd {params.D}; awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_16_lmm_1e5_bed_out:
    input: lambda wildcards: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for p in glob.glob(os.path.join(MARGINAL_TEST_COV_16,'results/lmm_1e5_nr',"*.txt"))]
    
marginal_test_cov_16_lmm_thresh = {'1e-5':'1e-5','1e-6':'1e-6','1e-7':'1e-7','1e-8':'1e-8','5e-8':'5e-8','1e-9':'1e-9','1e-10':'1e-10'}

rule marginal_test_cov_16_lmm_bonfthr_merge10k:
    input:  bed=MARGINAL_TEST_COV_16+'/results/lmm_master/allgenesbonf_{thr}.bed'
    params: D1=MARGINAL_TEST_COV_16+'/results/lmm_bonf{thr}',D2=MARGINAL_TEST_COV_16+'/results/lmm_bonf{thr}_merge10k'
    output: bed=MARGINAL_TEST_COV_16+'/results/lmm_master/allgenesbonf_{thr}_merge10k.bed'
    shell: """
        module load bedops/2.4.15 bedtools/2.25.0; mkdir -p {params.D1} {params.D2}
        rm -f {params.D1}/*; awk -v OFS='\\t' '{{ print > "{params.D1}/"$4".bed" }}' {input.bed}
        rm -f {params.D2}/*; cd {params.D1}; ls | parallel -j 12 "bedtools merge -d 10000 -c 4,4,5 -o distinct,count,min -i {{1}} > {params.D2}/{{1}}"
        cd {params.D2}; cat *.bed | sort-bed - > {output.bed}
    """
rule marginal_test_cov_16_lmm_bonfthr_merge10k_out:
#    input: [p.replace('/lmm_1e5_bed/','/lmm_bonf5e-2_merge10k/') for p in glob.glob(MARGINAL_TEST_COV_16+'/results/lmm_1e5_bed/*.bed')]
    input: expand(MARGINAL_TEST_COV_16+'/results/lmm_master/allgenesbonf_{thr}_merge10k.bed',thr=['5e-2','1e-1'])
rule rule marginal_test_cov_16_lmm_bonfthr_merge10k_out2:
    shell: "echo {rules.marginal_test_cov_16_lmm_bonfthr_merge10k_out.input}"

rule marginal_test_cov_16_lmm_master:
    input: glob.glob(os.path.join(MARGINAL_TEST_COV_16,'results/lmm_1e5_bed',"*.bed"))
    params: D=os.path.join(MARGINAL_TEST_COV_16,'results/lmm_1e5_bed')
    output: MARGINAL_TEST_COV_16+'/results/lmm_master/allgenes_{thr}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<={wildcards.thr}' > {output}
    """
rule marginal_test_cov_16_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_16+'/results/lmm_master/allgenes_{thr}.bed',thr=list(marginal_test_cov_16_lmm_thresh.keys()))

rule marginal_test_cov_16_lmm_bonf_master:
    input: MARGINAL_TEST_COV_16+'/results/lmm_master/allgenes_1e-5.bed'
    params: p=lambda wildcards: str(float(wildcards.thr)/(sum(1 for line in open(rules.marginal_test_cov_16_pos2.input[0]))-1)/marginal_test_cov_16_ngenes)
    output: MARGINAL_TEST_COV_16+'/results/lmm_master/allgenesbonf_{thr}.bed'
    shell: """
        awk -v OFS='\\t' '$5<={params.p}' {input} > {output}
    """
rule marginal_test_cov_16_lmm_bonf_master_out:
    input: expand(MARGINAL_TEST_COV_16+'/results/lmm_master/allgenesbonf_{thr}.bed',thr=['5e-2','1e-1'])

rule marginal_test_cov_16_lmm_merge_master:
    input: MARGINAL_TEST_COV_16+'/results/lmm_{thr}_{merge}/AT1G01010.bed'
    params: D=MARGINAL_TEST_COV_16+'/results/lmm_{thr}_{merge}'
    output: allgenes=MARGINAL_TEST_COV_16+'/results/lmm_master/allgenes_{thr}_{merge}.bed',chr15genes=MARGINAL_TEST_COV_16+'/results/lmm_master/chr1-5genes_{thr}_{merge}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | sort-bed - > {output.allgenes}
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {output.allgenes} > {output.chr15genes}
    """
rule marginal_test_cov_16_lmm_merge_master_out:
    input: expand(MARGINAL_TEST_COV_16+'/results/lmm_master/allgenes_{thr}_{merge}.bed',thr=['5e-2','1e-1'],merge=['merge10k'])

rule marginal_test_cov_16_lmm_master_noCM:
    input: MARGINAL_TEST_COV_16+'/results/lmm_master/allgenes{suff}.bed'
    output: MARGINAL_TEST_COV_16+'/results/lmm_master/chr1-5genes{suff}.bed'
    shell: """
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {input} > {output}
    """
rule marginal_test_cov_16_lmm_master_noCM_out:
    input: expand(MARGINAL_TEST_COV_16+'/results/lmm_master/chr1-5genes_{thr}.bed',thr=list(marginal_test_cov_16_lmm_thresh.keys()))  + expand(MARGINAL_TEST_COV_16+'/results/lmm_master/chr1-5genesbonf_{thr}.bed',thr=['5e-2','1e-1'])
rule marginal_test_cov_16_lmm_master_sorted:
    input: MARGINAL_TEST_COV_16+'/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_16+'/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_16_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_16+'/results/lmm_master/{geneset}_{thr}.bed.sorted',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_16_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_16+'/results/lmm_master/{geneset}_{thr}.bed.sorted',geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_16_lmm_master_uniq:
    input: MARGINAL_TEST_COV_16+'/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_16+'/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_16_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_16+'/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_16_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_16+'/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_16_lmm_master_tss:
    input: master_bed=MARGINAL_TEST_COV_16+'/results/lmm_master/{pref}.bed',tss_bed=GENOMES_TAIR10_ENS25_TSS
    output: MARGINAL_TEST_COV_16+'/results/lmm_master/{pref}.bed.tss'
    shell: """
            sort -k4b,4 {input.master_bed} | join -t $'\\t' -1 4 -2 4 -o 1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4b,4 {input.tss_bed}) > {output}
    """
rule marginal_test_cov_16_lmm_master_tss_out: 
    input: expand(MARGINAL_TEST_COV_16+'/results/lmm_master/{geneset}_{thr}.bed.tss',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_16_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_16+'/results/lmm_master/{geneset}_{thr}.bed.tss',geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_16_lmm_master_tss_out2: 
    shell: "echo {rules.marginal_test_cov_16_lmm_master_tss_out.input}"

# genotypes are dmC bins; RUVg k=4, UQ RIN normalized, genes filtered by CV>0.05 on UQ RUVg normalized, VST counts, covariate use UQ RIN RUVg covariates
rule marginal_test_cov_17:
    input: geno_file=DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5',K_file=MPI_RELEASE_v31+'/X1001tx_filter1/norm_cov_1001tx_filter1.csv',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_W_k4.txt'
    params: D=MARGINAL_TEST_COV_17+'/{mc_prefix}_chr{chrom}/',norm_mode='None',qname='gale.q',N="{mc_prefix}_chr{chrom}_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr{chrom}/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v8.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.K_file} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 8366 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
marginal_test_cov_17_ngenes = 8366
marginal_test_cov_17_slices = list(range(0,marginal_test_cov_17_ngenes,100))+[marginal_test_cov_17_ngenes]
rule marginal_test_cov_17_out:
    input: [MARGINAL_TEST_COV_17+'/%s_chr%d/logs/%04d-%04d.log'%(mc_prefix,c,marginal_test_cov_17_slices[i],marginal_test_cov_17_slices[i+1]) for mc_prefix in ['dmCH_filtered','dmCG_filtered','dmC_filtered'] for c in ATH_AUTOSOMES for i in range(len(marginal_test_cov_17_slices)-1)]
rule marginal_test_cov_17_out2:
    shell: "echo {rules.marginal_test_cov_17_out.input}"     
#snakemake -j 30 -k -p marginal_test_cov_17_out --cluster "mkdir -p {params.D}; qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

rule marginal_test_cov_17_pos:
    input: geno_file=DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5',
    params: D=MARGINAL_TEST_COV_17+'/{mc_prefix}_chr{chrom}/geno_pos'
    output: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr{chrom}/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_17_pos_out:
    input: expand(MARGINAL_TEST_COV_17+'/{mc_prefix}_chr{chrom}/geno_pos/get_geno_pos.log',mc_prefix=['dmCH_filtered','dmCG_filtered','dmC_filtered'],chrom=ATH_AUTOSOMES)
rule marginal_test_cov_17_pos2:
    input: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr%s/geno_pos/position.txt'%(chrom) for chrom in ATH_AUTOSOMES
    output: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/geno_pos/position.txt'
    shell: """
        awk -v OFS='\\t' 'NR==FNR {{ print; next }} {{ if (FNR>1) {{ print }} }}' {input} > {output}
    """
rule marginal_test_cov_17_pos2_out:
    input: expand(MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/geno_pos/position.txt',mc_prefix=['dmCH_filtered','dmCG_filtered','dmC_filtered'])
    shell: "echo {input}"
rule marginal_test_cov_17_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_17+'/{mc_prefix}_chr{chrom}/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr{chrom}/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_17_lmm_1e5_nr_out:
    input: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for mc_prefix in DMC_BINS_PREFIX for chrom in ATH_AUTOSOMES for p in glob.glob(MARGINAL_TEST_COV_17+'/%s_chr%s/results/lmm_pval/*.txt'%(mc_prefix,chrom))]
rule marginal_test_cov_17_lmm_1e5_nr_out2:
    shell: "echo {rules.marginal_test_cov_17_lmm_1e5_nr_out.input}"
rule marginal_test_cov_17_lmm_1e5_nr2:
    input:  [MARGINAL_TEST_COV_17+'/{mc_prefix}_chr%s/results/lmm_pval/{at_id}.txt'%chrom for chrom in ATH_AUTOSOMES]
    output: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '{{ if (FNR==1) file_num+=1 }} {{ if ($1<=1e-5) {{ print NR-file_num,$1 }} }}' {input} > {output}
    """
rule marginal_test_cov_17_lmm_1e5_nr2_out:
    input: [os.path.join(MARGINAL_TEST_COV_17,mc_prefix,'results/lmm_1e5_nr',os.path.basename(p)) for mc_prefix in DMC_BINS_PREFIX for p in glob.glob(MARGINAL_TEST_COV_17+'/%s_chr1/results/lmm_pval/*.txt'%(mc_prefix))]
rule marginal_test_cov_17_lmm_1e5_nr2_out2:
    shell: "echo {rules.marginal_test_cov_17_lmm_1e5_nr2_out.input}"
rule marginal_test_cov_17_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/geno_pos/position.txt',nr=MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_17 + '/{mc_prefix}_chr1-5/results/lmm_1e5_bed'
    output: bed=MARGINAL_TEST_COV_17 + '/{mc_prefix}_chr1-5/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2+99,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_17_lmm_1e5_bed_out:
    input: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for mc_prefix in DMC_BINS_PREFIX for p in glob.glob(MARGINAL_TEST_COV_17+'/%s_chr1-5/results/lmm_1e5_nr/*.txt'%(mc_prefix))]

#marginal_test_cov_17_lmm_thresh = {'1e-5':'1e-5','1e-6':'1e-6','1e-7':'1e-7','1e-8':'1e-8','5e-8':'5e-8','1e-9':'1e-9','1e-10':'1e-10'}
marginal_test_cov_17_lmm_thresh = {'1e-5':'1e-5'}
rule marginal_test_cov_17_lmm_master:
    input: glob.glob(os.path.join(MARGINAL_TEST_COV_17,'{mc_prefix}_chr1-5/results/lmm_1e5_bed',"*.bed"))
    params: D=MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_1e5_bed'
    output: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_{thr}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<={wildcards.thr}' > {output}
    """
rule marginal_test_cov_17_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=list(marginal_test_cov_17_lmm_thresh.keys()))
rule marginal_test_cov_17_lmm_bonf_master:
    input: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_1e-5.bed'
    params: p=lambda wildcards: str(float(wildcards.thr)/(sum(1 for line in open(MARGINAL_TEST_COV_17+'/'+wildcards.mc_prefix+'_chr1-5/geno_pos/position.txt'))-1)/marginal_test_cov_17_ngenes)
    output: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/allgenesbonf_{thr}.bed'
    shell: """
        awk -v OFS='\\t' '$5<={params.p}' {input} > {output}
    """
rule marginal_test_cov_17_lmm_bonf_master_out:
    input: expand(MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/allgenesbonf_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=['5e-2','1e-1'])

rule marginal_test_cov_17_lmm_master_noCM:
    input: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes{suff}.bed'
    output: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genes{suff}.bed'
    shell: """
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {input} > {output}
    """
rule marginal_test_cov_17_lmm_master_noCM_out:
    input: expand(MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genes_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=list(marginal_test_cov_17_lmm_thresh.keys())) + expand(MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genesbonf_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=['5e-2','1e-1'])
rule marginal_test_cov_17_lmm_master_sorted:
    input: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_17_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_17_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_17_lmm_master_uniq:
    input: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_17_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_17_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])

# genotypes are dmC sites; RUVg k=4, UQ RIN normalized, genes filtered by CV>0.05 on UQ RUVg normalized, VST counts, covariate use UQ RIN RUVg covariates
rule marginal_test_cov_18:
    input: geno_file=DMS+'/{mc_prefix}/{mc_prefix}_{chrom}_1001g_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5',K_file=MPI_RELEASE_v31+'/X1001tx_filter1/norm_cov_1001tx_filter1.csv',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_W_k4.txt'
    params: D=MARGINAL_TEST_COV_18+'/{mc_prefix}_chr{chrom}/',norm_mode='None',qname='gale.q',N="{mc_prefix}_chr{chrom}_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr{chrom}/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v8.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.K_file} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 8366 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
marginal_test_cov_18_ngenes = 8366
marginal_test_cov_18_slices = list(range(0,marginal_test_cov_18_ngenes,100))+[marginal_test_cov_18_ngenes]
rule marginal_test_cov_18_out:
    input: [MARGINAL_TEST_COV_18+'/%s_chr%d/logs/%04d-%04d.log'%(mc_prefix,c,marginal_test_cov_18_slices[i],marginal_test_cov_18_slices[i+1]) for mc_prefix in DMS_TABLES_PREFIX for c in ATH_AUTOSOMES for i in range(len(marginal_test_cov_18_slices)-1)]
rule marginal_test_cov_18_out2:
    shell: "echo {rules.marginal_test_cov_18_out.input}"     
#snakemake -j 30 -k -p marginal_test_cov_18_out --cluster "mkdir -p {params.D}; qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

rule marginal_test_cov_18_pos:
    input: geno_file=DMS+'/{mc_prefix}/{mc_prefix}_{chrom}_1001g_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5',
    params: D=MARGINAL_TEST_COV_18+'/{mc_prefix}_chr{chrom}/geno_pos'
    output: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr{chrom}/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_18_pos_out:
    input: expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr{chrom}/geno_pos/get_geno_pos.log',mc_prefix=DMS_TABLES_PREFIX,chrom=ATH_AUTOSOMES)
rule marginal_test_cov_18_pos2:
    input: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr%s/geno_pos/position.txt'%(chrom) for chrom in ATH_AUTOSOMES
    output: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/geno_pos/position.txt'
    shell: """
        awk -v OFS='\\t' 'NR==FNR {{ print; next }} {{ if (FNR>1) {{ print }} }}' {input} > {output}
    """
rule marginal_test_cov_18_pos2_out:
    input: expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/geno_pos/position.txt',mc_prefix=DMS_TABLES_PREFIX)
rule marginal_test_cov_18_lmm_1e5_nr:
    input:  [MARGINAL_TEST_COV_18+'/{mc_prefix}_chr%s/results/lmm_pval/{at_id}.txt'%chrom for chrom in ATH_AUTOSOMES]
    output: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '{{ if (FNR==1) file_num+=1 }} {{ if ($1<=1e-5) {{ print NR-file_num,$1 }} }}' {input} > {output}
    """
rule marginal_test_cov_18_lmm_1e5_nr_out:
    input: [os.path.join(MARGINAL_TEST_COV_18,mc_prefix+'_chr1-5','results/lmm_1e5_nr',os.path.basename(p)) for mc_prefix in DMS_TABLES_PREFIX for p in glob.glob(MARGINAL_TEST_COV_18+'/%s_chr1/results/lmm_pval/*.txt'%(mc_prefix))]
rule marginal_test_cov_18_lmm_1e5_nr_out2:
    shell: "echo {rules.marginal_test_cov_18_lmm_1e5_nr_out.input}"

rule marginal_test_cov_18_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/geno_pos/position.txt',nr=MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_18 + '/{mc_prefix}_chr1-5/results/lmm_1e5_bed'
    output: bed=MARGINAL_TEST_COV_18 + '/{mc_prefix}_chr1-5/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_18_lmm_1e5_bed_out:
    input: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for mc_prefix in DMS_TABLES_PREFIX for p in glob.glob(MARGINAL_TEST_COV_18+'/%s_chr1-5/results/lmm_1e5_nr/*.txt'%(mc_prefix))]

#marginal_test_cov_18_lmm_thresh = {'1e-5':'1e-5','1e-6':'1e-6','1e-7':'1e-7','1e-8':'1e-8','5e-8':'5e-8','1e-9':'1e-9','1e-10':'1e-10'}
marginal_test_cov_18_lmm_thresh = {'1e-5':'1e-5'}
rule marginal_test_cov_18_lmm_master:
    input: glob.glob(os.path.join(MARGINAL_TEST_COV_18,'{mc_prefix}_chr1-5/results/lmm_1e5_bed',"*.bed"))
    params: D=MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_1e5_bed'
    output: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_{thr}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<={wildcards.thr}' > {output}
    """
rule marginal_test_cov_18_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_{thr}.bed',mc_prefix=DMS_TABLES_PREFIX,thr=list(marginal_test_cov_18_lmm_thresh.keys()))
rule marginal_test_cov_18_lmm_bonf_master:
    input: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_1e-5.bed'
    params: p=lambda wildcards: str(float(wildcards.thr)/(sum(1 for line in open(MARGINAL_TEST_COV_18+'/'+wildcards.mc_prefix+'_chr1-5/geno_pos/position.txt'))-1)/marginal_test_cov_18_ngenes)
    output: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/allgenesbonf_{thr}.bed'
    shell: """
        awk -v OFS='\\t' '$5<={params.p}' {input} > {output}
    """
rule marginal_test_cov_18_lmm_bonf_master_out:
    input: expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/allgenesbonf_{thr}.bed',mc_prefix=DMS_TABLES_PREFIX,thr=['5e-2','1e-1'])

rule marginal_test_cov_18_lmm_master_noCM:
    input: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes{suff}.bed'
    output: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genes{suff}.bed'
    shell: """
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {input} > {output}
    """
rule marginal_test_cov_18_lmm_master_noCM_out:
    input: expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genes_{thr}.bed',mc_prefix=DMS_TABLES_PREFIX,thr=list(marginal_test_cov_18_lmm_thresh.keys())) + expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genesbonf_{thr}.bed',mc_prefix=DMS_TABLES_PREFIX,thr=['5e-2','1e-1'])
rule marginal_test_cov_18_lmm_master_sorted:
    input: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_18_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_18_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_18_lmm_master_uniq:
    input: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_18_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_18_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])

rule marginal_test_cov_18_lmm_master_tss:
    input: master_bed=MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed',tss_bed=GENOMES_TAIR10_ENS25_TSS
    output: MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.tss'
    shell: """
            sort -k4b,4 {input.master_bed} | join -t $'\\t' -1 4 -2 4 -o 1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4b,4 {input.tss_bed}) > {output}
    """
rule marginal_test_cov_18_lmm_master_tss_out: 
    input: expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.tss',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_18_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.tss',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_18_lmm_master_tss_out2: 
    shell: "echo {rules.marginal_test_cov_18_lmm_master_tss_out.input}"

#raw counts from filtered FASTQ, RUVg k=4, UQ RIN normalized, genes filtered by CV>0.05 on UQ RUVg normalized, VST counts; covariate use UQ RIN RUVg covariates
rule marginal_test_cov_19:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_W_k4.txt'
    params: D=MARGINAL_TEST_COV_19,norm_mode='None',qname='gale.q',N="assoc_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_19+'/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v6.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 8768 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
marginal_test_cov_19_ngenes = 8768
marginal_test_cov_19_slices = list(range(0,marginal_test_cov_19_ngenes,100))+[marginal_test_cov_19_ngenes]
rule marginal_test_cov_19_out:
    input: [MARGINAL_TEST_COV_19+'/logs/%04d-%04d.log'%(marginal_test_cov_19_slices[i],marginal_test_cov_19_slices[i+1]) for i in range(len(marginal_test_cov_19_slices)-1)]
rule marginal_test_cov_19_out2:
    shell: "echo {rules.marginal_test_cov_19_out.input}"     
#snakemake -j 20 -k -p marginal_test_cov_19_out --cluster "qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"
rule marginal_test_cov_19_pos:
    input: geno_file=MPI_RELEASE_v31+'/1001genomes_snp-short-indel_only_ACGTN_1001tx_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
    params: D=MARGINAL_TEST_COV_19+'/geno_pos'
    output: MARGINAL_TEST_COV_19+'/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_19_pos2:
    input: pos_txt=MARGINAL_TEST_COV_19+'/geno_pos/position.txt'
    output: pos_bed=MARGINAL_TEST_COV_19+'/geno_pos/position.bed'
    shell: """
        tail -n +2 {input.pos_txt} | awk -v OFS='\\t' '{{ print "chr"$1,$2-1,$2 }}'  > {output.pos_bed}
    """
rule marginal_test_cov_19_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_19 + '/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_19+'/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_19_lmm_1e5_nr_out:
    input: lambda wildcards: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for p in glob.glob(os.path.join(MARGINAL_TEST_COV_19,'results/lmm_pval',"*.txt"))]

rule marginal_test_cov_19_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_19+'/geno_pos/position.txt',nr=MARGINAL_TEST_COV_19 + '/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_19 + '/results/lmm_1e5_nr'
    output: bed=MARGINAL_TEST_COV_19+'/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        cd {params.D}; awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_19_lmm_1e5_bed_out:
    input: lambda wildcards: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for p in glob.glob(os.path.join(MARGINAL_TEST_COV_19,'results/lmm_1e5_nr',"*.txt"))]
    
marginal_test_cov_19_lmm_thresh = {'1e-5':'1e-5','1e-6':'1e-6','1e-7':'1e-7','1e-8':'1e-8','5e-8':'5e-8','1e-9':'1e-9','1e-10':'1e-10'}

rule marginal_test_cov_19_lmm_bonfthr_merge10k:
    input:  bed=MARGINAL_TEST_COV_19+'/results/lmm_master/allgenesbonf_{thr}.bed'
    params: D1=MARGINAL_TEST_COV_19+'/results/lmm_bonf{thr}',D2=MARGINAL_TEST_COV_19+'/results/lmm_bonf{thr}_merge10k'
    output: bed=MARGINAL_TEST_COV_19+'/results/lmm_master/allgenesbonf_{thr}_merge10k.bed'
    shell: """
        module load bedops/2.4.15 bedtools/2.25.0; mkdir -p {params.D1} {params.D2}
        rm -f {params.D1}/*; awk -v OFS='\\t' '{{ print > "{params.D1}/"$4".bed" }}' {input.bed}
        rm -f {params.D2}/*; cd {params.D1}; ls | parallel -j 12 "bedtools merge -d 10000 -c 4,4,5 -o distinct,count,min -i {{1}} > {params.D2}/{{1}}"
        cd {params.D2}; cat *.bed | sort-bed - > {output.bed}
    """
rule marginal_test_cov_19_lmm_bonfthr_merge10k_out:
#    input: [p.replace('/lmm_1e5_bed/','/lmm_bonf5e-2_merge10k/') for p in glob.glob(MARGINAL_TEST_COV_19+'/results/lmm_1e5_bed/*.bed')]
    input: expand(MARGINAL_TEST_COV_19+'/results/lmm_master/allgenesbonf_{thr}_merge10k.bed',thr=['5e-2','1e-1'])
rule rule marginal_test_cov_19_lmm_bonfthr_merge10k_out2:
    shell: "echo {rules.marginal_test_cov_19_lmm_bonfthr_merge10k_out.input}"

rule marginal_test_cov_19_lmm_master:
    input: glob.glob(os.path.join(MARGINAL_TEST_COV_19,'results/lmm_1e5_bed',"*.bed"))
    params: D=os.path.join(MARGINAL_TEST_COV_19,'results/lmm_1e5_bed')
    output: MARGINAL_TEST_COV_19+'/results/lmm_master/allgenes_{thr}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<={wildcards.thr}' > {output}
    """
rule marginal_test_cov_19_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_19+'/results/lmm_master/allgenes_{thr}.bed',thr=list(marginal_test_cov_19_lmm_thresh.keys()))

rule marginal_test_cov_19_lmm_bonf_master:
    input: MARGINAL_TEST_COV_19+'/results/lmm_master/allgenes_1e-5.bed'
    params: p=lambda wildcards: str(float(wildcards.thr)/(sum(1 for line in open(rules.marginal_test_cov_19_pos2.input[0]))-1)/marginal_test_cov_19_ngenes)
    output: MARGINAL_TEST_COV_19+'/results/lmm_master/allgenesbonf_{thr}.bed'
    shell: """
        awk -v OFS='\\t' '$5<={params.p}' {input} > {output}
    """
rule marginal_test_cov_19_lmm_bonf_master_out:
    input: expand(MARGINAL_TEST_COV_19+'/results/lmm_master/allgenesbonf_{thr}.bed',thr=['5e-2','1e-1'])

rule marginal_test_cov_19_lmm_merge_master:
    input: MARGINAL_TEST_COV_19+'/results/lmm_{thr}_{merge}/AT1G01010.bed'
    params: D=MARGINAL_TEST_COV_19+'/results/lmm_{thr}_{merge}'
    output: allgenes=MARGINAL_TEST_COV_19+'/results/lmm_master/allgenes_{thr}_{merge}.bed',chr15genes=MARGINAL_TEST_COV_19+'/results/lmm_master/chr1-5genes_{thr}_{merge}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | sort-bed - > {output.allgenes}
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {output.allgenes} > {output.chr15genes}
    """
rule marginal_test_cov_19_lmm_merge_master_out:
    input: expand(MARGINAL_TEST_COV_19+'/results/lmm_master/allgenes_{thr}_{merge}.bed',thr=['5e-2','1e-1'],merge=['merge10k'])

rule marginal_test_cov_19_lmm_master_noCM:
    input: MARGINAL_TEST_COV_19+'/results/lmm_master/allgenes{suff}.bed'
    output: MARGINAL_TEST_COV_19+'/results/lmm_master/chr1-5genes{suff}.bed'
    shell: """
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {input} > {output}
    """
rule marginal_test_cov_19_lmm_master_noCM_out:
    input: expand(MARGINAL_TEST_COV_19+'/results/lmm_master/chr1-5genes_{thr}.bed',thr=list(marginal_test_cov_19_lmm_thresh.keys()))  + expand(MARGINAL_TEST_COV_19+'/results/lmm_master/chr1-5genesbonf_{thr}.bed',thr=['5e-2','1e-1'])
rule marginal_test_cov_19_lmm_master_sorted:
    input: MARGINAL_TEST_COV_19+'/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_19+'/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_19_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_19+'/results/lmm_master/{geneset}_{thr}.bed.sorted',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_19_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_19+'/results/lmm_master/{geneset}_{thr}.bed.sorted',geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_19_lmm_master_uniq:
    input: MARGINAL_TEST_COV_19+'/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_19+'/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_19_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_19+'/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_19_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_19+'/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_19_lmm_master_tss:
    input: master_bed=MARGINAL_TEST_COV_19+'/results/lmm_master/{pref}.bed',tss_bed=GENOMES_TAIR10_ENS25_TSS
    output: MARGINAL_TEST_COV_19+'/results/lmm_master/{pref}.bed.tss'
    shell: """
            sort -k4b,4 {input.master_bed} | join -t $'\\t' -1 4 -2 4 -o 1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4b,4 {input.tss_bed}) > {output}
    """
rule marginal_test_cov_19_lmm_master_tss_out: 
    input: expand(MARGINAL_TEST_COV_19+'/results/lmm_master/{geneset}_{thr}.bed.tss',geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_19_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_19+'/results/lmm_master/{geneset}_{thr}.bed.tss',geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_19_lmm_master_tss_out2: 
    shell: "echo {rules.marginal_test_cov_19_lmm_master_tss_out.input}"

# genotypes are dmC bins; raw counts from filtered FASTQ, RUVg k=4, UQ RIN normalized, genes filtered by CV>0.05 on UQ RUVg normalized, VST counts; covariate use UQ RIN RUVg covariates
rule marginal_test_cov_20:
    input: geno_file=DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5',K_file=MPI_RELEASE_v31+'/X1001tx_filter1/norm_cov_1001tx_filter1.csv',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_W_k4.txt'
    params: D=MARGINAL_TEST_COV_20+'/{mc_prefix}_chr{chrom}/',norm_mode='None',qname='gale.q',N="{mc_prefix}_chr{chrom}_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr{chrom}/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v8.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.K_file} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 8768 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
marginal_test_cov_20_ngenes = 8768
marginal_test_cov_20_slices = list(range(0,marginal_test_cov_20_ngenes,100))+[marginal_test_cov_20_ngenes]
rule marginal_test_cov_20_out:
    input: [MARGINAL_TEST_COV_20+'/%s_chr%d/logs/%04d-%04d.log'%(mc_prefix,c,marginal_test_cov_20_slices[i],marginal_test_cov_20_slices[i+1]) for mc_prefix in ['dmCH_filtered','dmCG_filtered','dmC_filtered'] for c in ATH_AUTOSOMES for i in range(len(marginal_test_cov_20_slices)-1)]
rule marginal_test_cov_20_out2:
    shell: "echo {rules.marginal_test_cov_20_out.input}"     
#snakemake -j 30 -k -p marginal_test_cov_20_out --cluster "mkdir -p {params.D}; qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

rule marginal_test_cov_20_pos:
    input: geno_file=DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
    params: D=MARGINAL_TEST_COV_20+'/{mc_prefix}_chr{chrom}/geno_pos'
    output: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr{chrom}/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_20_pos_out:
    input: expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr{chrom}/geno_pos/get_geno_pos.log',mc_prefix=['dmCH_filtered','dmCG_filtered','dmC_filtered'],chrom=ATH_AUTOSOMES)
rule marginal_test_cov_20_pos2:
    input: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr%s/geno_pos/position.txt'%(chrom) for chrom in ATH_AUTOSOMES
    output: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/geno_pos/position.txt'
    shell: """
        awk -v OFS='\\t' 'NR==FNR {{ print; next }} {{ if (FNR>1) {{ print }} }}' {input} > {output}
    """
rule marginal_test_cov_20_pos2_out:
    input: expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/geno_pos/position.txt',mc_prefix=['dmCH_filtered','dmCG_filtered','dmC_filtered'])
    shell: "echo {input}"
rule marginal_test_cov_20_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_20+'/{mc_prefix}_chr{chrom}/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr{chrom}/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_20_lmm_1e5_nr_out:
    input: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for mc_prefix in DMC_BINS_PREFIX for chrom in ATH_AUTOSOMES for p in glob.glob(MARGINAL_TEST_COV_20+'/%s_chr%s/results/lmm_pval/*.txt'%(mc_prefix,chrom))]
rule marginal_test_cov_20_lmm_1e5_nr_out2:
    shell: "echo {rules.marginal_test_cov_20_lmm_1e5_nr_out.input}"
rule marginal_test_cov_20_lmm_1e5_nr2:
    input:  [MARGINAL_TEST_COV_20+'/{mc_prefix}_chr%s/results/lmm_pval/{at_id}.txt'%chrom for chrom in ATH_AUTOSOMES]
    output: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '{{ if (FNR==1) file_num+=1 }} {{ if ($1<=1e-5) {{ print NR-file_num,$1 }} }}' {input} > {output}
    """
rule marginal_test_cov_20_lmm_1e5_nr2_out:
    input: [os.path.join(MARGINAL_TEST_COV_20,mc_prefix+'_chr1-5/results/lmm_1e5_nr',os.path.basename(p)) for mc_prefix in DMC_BINS_PREFIX for p in glob.glob(MARGINAL_TEST_COV_20+'/%s_chr1/results/lmm_pval/*.txt'%(mc_prefix))]
rule marginal_test_cov_20_lmm_1e5_nr2_out2:
    shell: "echo {rules.marginal_test_cov_20_lmm_1e5_nr2_out.input}"
rule marginal_test_cov_20_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/geno_pos/position.txt',nr=MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_20 + '/{mc_prefix}_chr1-5/results/lmm_1e5_bed'
    output: bed=MARGINAL_TEST_COV_20 + '/{mc_prefix}_chr1-5/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2+99,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_20_lmm_1e5_bed_out:
    input: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for mc_prefix in DMC_BINS_PREFIX for p in glob.glob(MARGINAL_TEST_COV_20+'/%s_chr1-5/results/lmm_1e5_nr/*.txt'%(mc_prefix))]

#marginal_test_cov_20_lmm_thresh = {'1e-5':'1e-5','1e-6':'1e-6','1e-7':'1e-7','1e-8':'1e-8','5e-8':'5e-8','1e-9':'1e-9','1e-10':'1e-10'}
marginal_test_cov_20_lmm_thresh = {'1e-5':'1e-5'}
rule marginal_test_cov_20_lmm_master:
    input: glob.glob(os.path.join(MARGINAL_TEST_COV_20,'{mc_prefix}_chr1-5/results/lmm_1e5_bed',"*.bed"))
    params: D=MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_1e5_bed'
    output: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_{thr}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<={wildcards.thr}' > {output}
    """
rule marginal_test_cov_20_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=list(marginal_test_cov_20_lmm_thresh.keys()))
rule marginal_test_cov_20_lmm_bonf_master:
    input: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_1e-5.bed'
    params: p=lambda wildcards: str(float(wildcards.thr)/(sum(1 for line in open(MARGINAL_TEST_COV_20+'/'+wildcards.mc_prefix+'_chr1-5/geno_pos/position.txt'))-1)/marginal_test_cov_20_ngenes)
    output: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/allgenesbonf_{thr}.bed'
    shell: """
        awk -v OFS='\\t' '$5<={params.p}' {input} > {output}
    """
rule marginal_test_cov_20_lmm_bonf_master_out:
    input: expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/allgenesbonf_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=['5e-2','1e-1'])

rule marginal_test_cov_20_lmm_master_noCM:
    input: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes{suff}.bed'
    output: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genes{suff}.bed'
    shell: """
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {input} > {output}
    """
rule marginal_test_cov_20_lmm_master_noCM_out:
    input: expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genes_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=list(marginal_test_cov_20_lmm_thresh.keys())) + expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genesbonf_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=['5e-2','1e-1'])
rule marginal_test_cov_20_lmm_master_sorted:
    input: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_20_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_20_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_20_lmm_master_uniq:
    input: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_20_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_20_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])

rule marginal_test_cov_20_lmm_master_tss:
    input: master_bed=MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed',tss_bed=GENOMES_TAIR10_ENS25_TSS
    output: MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.tss'
    shell: """
            sort -k4b,4 {input.master_bed} | join -t $'\\t' -1 4 -2 4 -o 1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4b,4 {input.tss_bed}) > {output}
    """
rule marginal_test_cov_20_lmm_master_tss_out: 
    input: expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.tss',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_20_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.tss',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_20_lmm_master_tss_out2: 
    shell: "echo {rules.marginal_test_cov_20_lmm_master_tss_out.input}"

# genotypes are dmC sites; RUVg k=4, UQ RIN normalized, genes filtered by CV>0.05 on UQ RUVg normalized, VST counts, covariate use UQ RIN RUVg covariates
rule marginal_test_cov_21:
    input: geno_file=DMS+'/{mc_prefix}/{mc_prefix}_{chrom}_1001g_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5',K_file=MPI_RELEASE_v31+'/X1001tx_filter1/norm_cov_1001tx_filter1.csv',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_W_k4.txt'
    params: D=MARGINAL_TEST_COV_21+'/{mc_prefix}_chr{chrom}/',norm_mode='None',qname='gale.q',N="{mc_prefix}_chr{chrom}_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr{chrom}/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v8.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.K_file} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 8768 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
marginal_test_cov_21_ngenes = 8768
marginal_test_cov_21_slices = list(range(0,marginal_test_cov_21_ngenes,100))+[marginal_test_cov_21_ngenes]
rule marginal_test_cov_21_out:
    input: [MARGINAL_TEST_COV_21+'/%s_chr%d/logs/%04d-%04d.log'%(mc_prefix,c,marginal_test_cov_21_slices[i],marginal_test_cov_21_slices[i+1]) for mc_prefix in DMS_TABLES_PREFIX for c in ATH_AUTOSOMES for i in range(len(marginal_test_cov_21_slices)-1)]
rule marginal_test_cov_21_out2:
    shell: "echo {rules.marginal_test_cov_21_out.input}"     
#snakemake -j 30 -k -p marginal_test_cov_21_out --cluster "mkdir -p {params.D}; qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

rule marginal_test_cov_21_pos:
    input: geno_file=DMS+'/{mc_prefix}/{mc_prefix}_{chrom}_1001g_filter1.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-02-06/ath1001_tx_norm_2016-02-06-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5',
    params: D=MARGINAL_TEST_COV_21+'/{mc_prefix}_chr{chrom}/geno_pos'
    output: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr{chrom}/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_21_pos_out:
    input: expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr{chrom}/geno_pos/get_geno_pos.log',mc_prefix=DMS_TABLES_PREFIX,chrom=ATH_AUTOSOMES)
rule marginal_test_cov_21_pos2:
    input: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr%s/geno_pos/position.txt'%(chrom) for chrom in ATH_AUTOSOMES
    output: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/geno_pos/position.txt'
    shell: """
        awk -v OFS='\\t' 'NR==FNR {{ print; next }} {{ if (FNR>1) {{ print }} }}' {input} > {output}
    """
rule marginal_test_cov_21_pos2_out:
    input: expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/geno_pos/position.txt',mc_prefix=DMS_TABLES_PREFIX)
rule marginal_test_cov_21_lmm_1e5_nr:
    input:  [MARGINAL_TEST_COV_21+'/{mc_prefix}_chr%s/results/lmm_pval/{at_id}.txt'%chrom for chrom in ATH_AUTOSOMES]
    output: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '{{ if (FNR==1) file_num+=1 }} {{ if ($1<=1e-5) {{ print NR-file_num,$1 }} }}' {input} > {output}
    """
rule marginal_test_cov_21_lmm_1e5_nr_out:
    input: [os.path.join(MARGINAL_TEST_COV_21,mc_prefix+'_chr1-5','results/lmm_1e5_nr',os.path.basename(p)) for mc_prefix in DMS_TABLES_PREFIX for p in glob.glob(MARGINAL_TEST_COV_21+'/%s_chr1/results/lmm_pval/*.txt'%(mc_prefix))]
rule marginal_test_cov_21_lmm_1e5_nr_out2:
    shell: "echo {rules.marginal_test_cov_21_lmm_1e5_nr_out.input}"

rule marginal_test_cov_21_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/geno_pos/position.txt',nr=MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_21 + '/{mc_prefix}_chr1-5/results/lmm_1e5_bed'
    output: bed=MARGINAL_TEST_COV_21 + '/{mc_prefix}_chr1-5/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_21_lmm_1e5_bed_out:
    input: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for mc_prefix in DMS_TABLES_PREFIX for p in glob.glob(MARGINAL_TEST_COV_21+'/%s_chr1-5/results/lmm_1e5_nr/*.txt'%(mc_prefix))]
rule marginal_test_cov_21_lmm_1e5_bed_out2:
    shell: "echo {rules.marginal_test_cov_21_lmm_1e5_bed_out.input}"

#marginal_test_cov_21_lmm_thresh = {'1e-5':'1e-5','1e-6':'1e-6','1e-7':'1e-7','1e-8':'1e-8','5e-8':'5e-8','1e-9':'1e-9','1e-10':'1e-10'}
marginal_test_cov_21_lmm_thresh = {'1e-5':'1e-5'}
rule marginal_test_cov_21_lmm_master:
    input: glob.glob(os.path.join(MARGINAL_TEST_COV_21,'{mc_prefix}_chr1-5/results/lmm_1e5_bed',"*.bed"))
    params: D=MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_1e5_bed'
    output: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_{thr}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<={wildcards.thr}' > {output}
    """
rule marginal_test_cov_21_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_{thr}.bed',mc_prefix=DMS_TABLES_PREFIX,thr=list(marginal_test_cov_21_lmm_thresh.keys()))
rule marginal_test_cov_21_lmm_bonf_master:
    input: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_1e-5.bed'
    params: p=lambda wildcards: str(float(wildcards.thr)/(sum(1 for line in open(MARGINAL_TEST_COV_21+'/'+wildcards.mc_prefix+'_chr1-5/geno_pos/position.txt'))-1)/marginal_test_cov_21_ngenes)
    output: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/allgenesbonf_{thr}.bed'
    shell: """
        awk -v OFS='\\t' '$5<={params.p}' {input} > {output}
    """
rule marginal_test_cov_21_lmm_bonf_master_out:
    input: expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/allgenesbonf_{thr}.bed',mc_prefix=DMS_TABLES_PREFIX,thr=['5e-2','1e-1'])

rule marginal_test_cov_21_lmm_master_noCM:
    input: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes{suff}.bed'
    output: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genes{suff}.bed'
    shell: """
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {input} > {output}
    """
rule marginal_test_cov_21_lmm_master_noCM_out:
    input: expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genes_{thr}.bed',mc_prefix=DMS_TABLES_PREFIX,thr=list(marginal_test_cov_21_lmm_thresh.keys())) + expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genesbonf_{thr}.bed',mc_prefix=DMS_TABLES_PREFIX,thr=['5e-2','1e-1'])
rule marginal_test_cov_21_lmm_master_sorted:
    input: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_21_lmm_master_sorted_out: 
    input: expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_21_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_21_lmm_master_uniq:
    input: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_21_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_21_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])

rule marginal_test_cov_21_lmm_master_tss:
    input: master_bed=MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed',tss_bed=GENOMES_TAIR10_ENS25_TSS
    output: MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.tss'
    shell: """
            sort -k4b,4 {input.master_bed} | join -t $'\\t' -1 4 -2 4 -o 1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4b,4 {input.tss_bed}) > {output}
    """
rule marginal_test_cov_21_lmm_master_tss_out: 
    input: expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.tss',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_21_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.tss',mc_prefix=DMS_TABLES_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_21_lmm_master_tss_out2: 
    shell: "echo {rules.marginal_test_cov_21_lmm_master_tss_out.input}"


# genotypes are dmC bins (corrected for trimming error); raw counts from filtered FASTQ, RUVg k=4, UQ RIN normalized, genes filtered by CV>0.05 on UQ RUVg normalized, VST counts; covariate use UQ RIN RUVg covariates
rule marginal_test_cov_22:
    input: geno_file=DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5',K_file=MPI_RELEASE_v31+'/X1001tx_filter1/norm_cov_1001tx_filter1.csv',cov_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_W_k4.txt'
    params: D=MARGINAL_TEST_COV_22+'/{mc_prefix}_chr{chrom}/',norm_mode='None',qname='gale.q',N="{mc_prefix}_chr{chrom}_{RNA_start}-{RNA_end}"
    output: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr{chrom}/logs/{RNA_start}-{RNA_end}.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/marginal_test_cov_v8.py {input.geno_file} {input.pheno_file} {params.norm_mode} {input.K_file} {input.cov_file} {wildcards.RNA_start} {wildcards.RNA_end} {params.D}
    """
# 8768 genes in PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
marginal_test_cov_22_ngenes = 8768
marginal_test_cov_22_slices = list(range(0,marginal_test_cov_22_ngenes,100))+[marginal_test_cov_22_ngenes]
rule marginal_test_cov_22_out:
    input: [MARGINAL_TEST_COV_22+'/%s_chr%d/logs/%04d-%04d.log'%(mc_prefix,c,marginal_test_cov_22_slices[i],marginal_test_cov_22_slices[i+1]) for mc_prefix in ['dmCH_filtered','dmCG_filtered','dmC_filtered'] for c in ATH_AUTOSOMES for i in range(len(marginal_test_cov_22_slices)-1)]
rule marginal_test_cov_22_out2:
    shell: "echo {rules.marginal_test_cov_22_out.input}"     
#snakemake -j 30 -k -p marginal_test_cov_22_out --cluster "mkdir -p {params.D}; qsub -wd {params.D} -l qname={params.qname} -N {params.N} -j y -V"

rule marginal_test_cov_22_pos:
    input: geno_file=DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.hdf5',pheno_file=PROJ_RESULTS_PATH_GALE+'/ath1001_tx_norm_2016-04-21/ath1001_tx_norm_2016-04-21-UQRIN_gNorm_k4_vst2_cv0p05_UQRIN_1001gT.hdf5'
    params: D=MARGINAL_TEST_COV_22+'/{mc_prefix}_chr{chrom}/geno_pos'
    output: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr{chrom}/geno_pos/get_geno_pos.log'
    shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        python {PROJ_PYTHON_PATH_GALE}/get_geno_pos.py {input.geno_file} {input.pheno_file} {params.D}
    """
rule marginal_test_cov_22_pos_out:
    input: expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr{chrom}/geno_pos/get_geno_pos.log',mc_prefix=['dmCH_filtered','dmCG_filtered','dmC_filtered'],chrom=ATH_AUTOSOMES)
rule marginal_test_cov_22_pos2:
    input: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr%s/geno_pos/position.txt'%(chrom) for chrom in ATH_AUTOSOMES
    output: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/geno_pos/position.txt'
    shell: """
        awk -v OFS='\\t' 'NR==FNR {{ print; next }} {{ if (FNR>1) {{ print }} }}' {input} > {output}
    """
rule marginal_test_cov_22_pos2_out:
    input: expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/geno_pos/position.txt',mc_prefix=['dmCH_filtered','dmCG_filtered','dmC_filtered'])
    shell: "echo {input}"
rule marginal_test_cov_22_lmm_1e5_nr:
    input:  MARGINAL_TEST_COV_22+'/{mc_prefix}_chr{chrom}/results/lmm_pval/{at_id}.txt'
    output: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr{chrom}/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '$1<=1e-5 {{ print NR-1,$1 }}' {input} > {output}
    """
rule marginal_test_cov_22_lmm_1e5_nr_out:
    input: [p.replace('/lmm_pval/','/lmm_1e5_nr/') for mc_prefix in DMC_BINS_PREFIX for chrom in ATH_AUTOSOMES for p in glob.glob(MARGINAL_TEST_COV_22+'/%s_chr%s/results/lmm_pval/*.txt'%(mc_prefix,chrom))]
rule marginal_test_cov_22_lmm_1e5_nr_out2:
    shell: "echo {rules.marginal_test_cov_22_lmm_1e5_nr_out.input}"
rule marginal_test_cov_22_lmm_1e5_nr2:
    input:  [MARGINAL_TEST_COV_22+'/{mc_prefix}_chr%s/results/lmm_pval/{at_id}.txt'%chrom for chrom in ATH_AUTOSOMES]
    output: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_1e5_nr/{at_id}.txt'
    shell: """
        awk -v OFS='\\t' '{{ if (FNR==1) file_num+=1 }} {{ if ($1<=1e-5) {{ print NR-file_num,$1 }} }}' {input} > {output}
    """
rule marginal_test_cov_22_lmm_1e5_nr2_out:
    input: [os.path.join(MARGINAL_TEST_COV_22,mc_prefix+'_chr1-5/results/lmm_1e5_nr',os.path.basename(p)) for mc_prefix in DMC_BINS_PREFIX for p in glob.glob(MARGINAL_TEST_COV_22+'/%s_chr1/results/lmm_pval/*.txt'%(mc_prefix))][15000:26000]
rule marginal_test_cov_22_lmm_1e5_nr2_out2:
    shell: "echo {rules.marginal_test_cov_22_lmm_1e5_nr2_out.input}"
rule marginal_test_cov_22_lmm_1e5_bed:
    input: geno_pos=MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/geno_pos/position.txt',nr=MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_1e5_nr/{at_id}.txt'
    params: D=MARGINAL_TEST_COV_22 + '/{mc_prefix}_chr1-5/results/lmm_1e5_bed'
    output: bed=MARGINAL_TEST_COV_22 + '/{mc_prefix}_chr1-5/results/lmm_1e5_bed/{at_id}.bed'
    shell: """
        awk -v OFS='\\t' 'NR==FNR {{ nums[$1+1]=$2; next}} FNR in nums {{ print "chr"$1,$2-1,$2+99,"{wildcards.at_id}",nums[FNR] }}' {input.nr} {input.geno_pos} > {output.bed}
    """
rule marginal_test_cov_22_lmm_1e5_bed_out:
    input: [os.path.splitext(p.replace('/lmm_1e5_nr/','/lmm_1e5_bed/'))[0]+'.bed' for mc_prefix in DMC_BINS_PREFIX for p in glob.glob(MARGINAL_TEST_COV_22+'/%s_chr1-5/results/lmm_1e5_nr/*.txt'%(mc_prefix))][13500:27000]

#marginal_test_cov_22_lmm_thresh = {'1e-5':'1e-5','1e-6':'1e-6','1e-7':'1e-7','1e-8':'1e-8','5e-8':'5e-8','1e-9':'1e-9','1e-10':'1e-10'}
marginal_test_cov_22_lmm_thresh = {'1e-5':'1e-5'}
rule marginal_test_cov_22_lmm_master:
    input: glob.glob(os.path.join(MARGINAL_TEST_COV_22,'{mc_prefix}_chr1-5/results/lmm_1e5_bed',"*.bed"))
    params: D=MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_1e5_bed'
    output: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_{thr}.bed'
    shell: """
        module load bedops/2.4.15; cd {params.D}; cat *.bed | awk -v OFS='\\t' '$5<={wildcards.thr}' > {output}
    """
rule marginal_test_cov_22_lmm_master_out:
    input: expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=list(marginal_test_cov_22_lmm_thresh.keys()))
rule marginal_test_cov_22_lmm_bonf_master:
    input: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes_1e-5.bed'
    params: p=lambda wildcards: str(float(wildcards.thr)/(sum(1 for line in open(MARGINAL_TEST_COV_22+'/'+wildcards.mc_prefix+'_chr1-5/geno_pos/position.txt'))-1)/marginal_test_cov_22_ngenes)
    output: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/allgenesbonf_{thr}.bed'
    shell: """
        awk -v OFS='\\t' '$5<={params.p}' {input} > {output}
    """
rule marginal_test_cov_22_lmm_bonf_master_out:
    input: expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/allgenesbonf_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=['5e-2','1e-1'])

rule marginal_test_cov_22_lmm_master_noCM:
    input: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/allgenes{suff}.bed'
    output: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genes{suff}.bed'
    shell: """
        awk -v OFS='\\t' '$4!~/^ATMG/ && $4!~/^ATCG/' {input} > {output}
    """
rule marginal_test_cov_22_lmm_master_noCM_out:
    input: expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genes_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=list(marginal_test_cov_22_lmm_thresh.keys())) + expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genesbonf_{thr}.bed',mc_prefix=DMC_BINS_PREFIX,thr=['5e-2','1e-1'])
rule marginal_test_cov_22_lmm_master_sorted:
    input: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed'
    output: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted'
    shell: """
        module load bedops/2.4.15; sort-bed {input} > {output}
    """
rule marginal_test_cov_22_lmm_master_sorted_out:
    input: expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_22_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_22_lmm_master_uniq:
    input: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted'
    output: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.sorted.uniq'
    shell: """
        cut -f 1-3 {input} | uniq > {output}
    """
rule marginal_test_cov_22_lmm_master_uniq_out: 
    input: expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_22_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.sorted.uniq',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])

rule marginal_test_cov_22_lmm_master_tss:
    input: master_bed=MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed',tss_bed=GENOMES_TAIR10_ENS25_TSS
    output: MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{pref}.bed.tss'
    shell: """
            sort -k4b,4 {input.master_bed} | join -t $'\\t' -1 4 -2 4 -o 1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4b,4 {input.tss_bed}) > {output}
    """
rule marginal_test_cov_22_lmm_master_tss_out: 
    input: expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.tss',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenes','chr1-5genes'],thr=marginal_test_cov_22_lmm_thresh.keys()) + expand(MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/{geneset}_{thr}.bed.tss',mc_prefix=DMC_BINS_PREFIX,geneset=['allgenesbonf','chr1-5genesbonf'],thr=['5e-2','1e-1'])
rule marginal_test_cov_22_lmm_master_tss_out2: 
    shell: "echo {rules.marginal_test_cov_22_lmm_master_tss_out.input}"


    
###################
# dmC data to bed
##################
rule dmC_bins_bed:
    input: expand(DMC_BINS+'/{{mc_prefix}}/{{mc_prefix}}_methylation_{chrom}.tsv',chrom=ATH_AUTOSOMES)
    output: DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation.bed'
    shell: """
        awk -v OFS='\\t' '{{ if (FNR>1) {{ print "chr"$1,$2-1,$3,"{wildcards.mc_prefix}" }} }}' {input} > {output}
    """
rule dmC_bins_bed_out:
    input: expand(DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation.bed',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])

rule dmC_bins_hdf5:
    input: DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.tsv'
    output: DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        python {PROJ_PYTHON_PATH_GALE}/limix_convert_methylpy_tsv.py --dmr={input} --outfile={output}
    """
rule dmC_bins_hdf5_out:
    input: expand(DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.hdf5',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'],chrom=ATH_AUTOSOMES)
rule dmC_bins_hdf5_out_clean:
    shell: "rm {rules.dmC_bins_hdf5_out.input}"

rule dmC_bins_new_bed:
    input: expand(DMC_BINS_NEW+'/{{mc_prefix}}/{{mc_prefix}}_methylation_{chrom}.tsv',chrom=ATH_AUTOSOMES)
    output: DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation.bed'
    shell: """
        awk -v OFS='\\t' '{{ if (FNR>1) {{ print "chr"$1,$2-1,$3,"{wildcards.mc_prefix}" }} }}' {input} > {output}
    """
rule dmC_bins_new_bed_out:
    input: expand(DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation.bed',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])

rule dmC_bins_new_hdf5:
    input: DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.tsv'
    output: DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        python {PROJ_PYTHON_PATH_GALE}/limix_convert_methylpy_tsv.py --dmr={input} --outfile={output}
    """
rule dmC_bins_new_hdf5_out:
    input: expand(DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation_{chrom}.hdf5',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'],chrom=ATH_AUTOSOMES)
rule dmC_bins_new_hdf5_out_clean:
    shell: "rm {rules.dmC_bins_new_hdf5_out.input}"


rule dms_bed:
    input: expand(DMS+'/{{mc_prefix}}/{{mc_prefix}}_{chrom}_1001g_filter1.tsv',chrom=ATH_AUTOSOMES)
    output: DMS+'/{mc_prefix}/{mc_prefix}_1001g_filter1.bed'
    shell: """
        awk -v OFS='\\t' '{{ if (FNR>1) {{ print "chr"$1,$2-1,$2,"{wildcards.mc_prefix}" }} }}' {input} > {output}
    """
rule dms_bed_out:
    input: expand(DMS+'/{mc_prefix}/{mc_prefix}_1001g_filter1.bed',mc_prefix=DMS_TABLES_PREFIX)

rule dms_hdf5:
    input: DMS+'/{mc_prefix}/{mc_prefix}_{chrom_and_filter}.tsv'
    output: DMS+'/{mc_prefix}/{mc_prefix}_{chrom_and_filter}.hdf5'
    shell: """
        module load shhuang anaconda/2.4.1
        python {PROJ_PYTHON_PATH_GALE}/limix_convert_methylpy_tsv.py --dms={input} --outfile={output}
    """
rule dms_hdf5_out:
    input: expand(DMS+'/{mc_prefix}/{mc_prefix}_{chrom}_1001g_filter1.hdf5',mc_prefix=['mCG_table','mCHG_table','mCHH_table'],chrom=ATH_AUTOSOMES)
rule dms_hdf5_out_clean:
    shell: "rm {rules.dms_hdf5_out.input}"

rule mQTL_bed:
    input: os.path.join(MQTL,'{mc_prefix}_DMR_mQTL.txt')
    output: os.path.join(MQTL,'{mc_prefix}_DMR/mQTL.bed')# chr(mQTL),start(mQTL),end(mQTL),chr:start(DMR),number of significant SNP
    shell: """
        mkdir -p {MQTL}/{wildcards.mc_prefix}_DMR
        module load bedops/2.4.15; awk -v OFS='\\t' '{{ print "chr"$3,$4,$5,"chr"$1":"$2,$6 }}' {input} | sort-bed - > {output}
    """
rule mQTL_bed_out:
    input: expand(MQTL+'/{mc_prefix}_DMR/mQTL.bed',mc_prefix=['CG','C'])

    
#############################
# DMR overlap with DAP sites
#############################
# dmC_bins
rule gat_mc_dap_01:# dmC bins and DAP TF peaks including ampDAP clustered by family, background is TAIR10 100bp uniquely mappable
   input:  mc_bed=DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation.bed',dap_peak=dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=GENOMES_TAIR10_MAP_UNIQ100,num_samples=10000,D=DMC_BINS_GAT01+'/{mc_prefix}'
   output: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dappeak1.out'
   log: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dappeak1.log'
   shell: """
       module load shhuang anaconda/2.4.1
       mkdir -p {params.D} 
       gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} {params.annot} --log={log} > {output}
   """
rule gat_mc_dap_01_out:
   input: expand(DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dappeak1.out',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])
rule gat_mc_dap_01_out_clean:
    shell: "rm {rules.gat_mc_dap_01_out.input}"

rule gat_mc_dap_02:# dmC bins and DAP TF motifs not including ampDAP, background is TAIR10 100bp uniquely mappable
   input: mc_bed=DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation.bed',dap_peak=dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3')],wrk=GENOMES_TAIR10_MAP_UNIQ100,num_samples=10000,D=DMC_BINS_GAT01+'/{mc_prefix}'
   output: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo1.out'
   log: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo1.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} {params.annot} --log={log} > {output}
   """
rule gat_mc_dap_02_out: 
   input: expand(DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo1.out',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])
rule gat_mc_dap_02_out_clean:
   shell: "rm {rules.gat_mc_dap_02_out.input}"

rule gat_mc_dap_03:# dmC bins and DAP TF peaks including ampDAP, background is clustered DAP TF peaks, including ampDAP
   input: mc_bed=DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation.bed'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/family_merge1.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/encode_cluster_01.bed3',num_samples=10000,D=DMC_BINS_GAT01+'/{mc_prefix}'
   output: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dappeak2.out'
   log: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dappeak2.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --output-counts-pattern={params.D}/{wildcards.mc_prefix}.%s.gat_dappeak2.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_mc_dap_03_out: 
   input: expand(DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dappeak2.out',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])
rule gat_mc_dap_03_out_clean:
   shell: "rm {rules.gat_mc_dap_03_out.input}"

rule gat_mc_dap_04:# dmC bins and DAP TF motifs not including ampDAP, background is clustered DAP TF motifs not including ampDAP
   input: mc_bed=DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation.bed'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',num_samples=10000,D=DMC_BINS_GAT01+'/{mc_prefix}'
   output: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo2.out'
   log: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo2.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --output-counts-pattern={params.D}/{wildcards.mc_prefix}.%s.gat_dapfimo2.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_mc_dap_04_out: 
   input: expand(DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo2.out',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])
rule gat_mc_dap_04_out_clean:
   shell: "rm {rules.gat_mc_dap_04_out.input}"

rule gat_mc_dap_06:# dmC bins eQTL and DAP TF motifs not including ampDAP, background is clustered DAP TF motifs not including ampDAP
   input: mc_bed=MARGINAL_TEST_COV_17+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genesbonf_5e-2.bed.sorted.uniq'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',num_samples=10000,D=DMC_BINS_GAT01+'/{mc_prefix}'
   output: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo3.out'
   log: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo3.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --output-counts-pattern={params.D}/{wildcards.mc_prefix}.%s.gat_dapfimo2.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_mc_dap_06_out: 
   input: expand(DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo3.out',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])
rule gat_mc_dap_06_out_clean:
   shell: "rm {rules.gat_mc_dap_06_out.input}"

rule gat_mc_dap_07:# dmC bins eQTL and DAP TF motifs not including ampDAP, background is clustered DAP TF motifs not including ampDAP
   input: mc_bed=MARGINAL_TEST_COV_20+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genesbonf_5e-2.bed.sorted.uniq'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',num_samples=10000,D=DMC_BINS_GAT01+'/{mc_prefix}'
   output: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo4.out'
   log: DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo4.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --output-counts-pattern={params.D}/{wildcards.mc_prefix}.%s.gat_dapfimo4.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_mc_dap_07_out: 
   input: expand(DMC_BINS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo4.out',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])
rule gat_mc_dap_07_out_clean:
   shell: "rm {rules.gat_mc_dap_07_out.input}"

rule gat_mc_dap_08:# dmC bins (after trimming error was fixed May-2016) and DAP TF peaks including ampDAP clustered by family, background is TAIR10 100bp uniquely mappable
   input: mc_bed=DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation.bed',dap_peak=dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=GENOMES_TAIR10_MAP_UNIQ100,num_samples=10000,D=DMC_BINS_GAT02+'/{mc_prefix}'
   output: DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dappeak1.out'
   log: DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dappeak1.log'
   shell: """
       module load shhuang anaconda/2.4.1
       mkdir -p {params.D} 
       gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} {params.annot} --log={log} > {output}
   """
rule gat_mc_dap_08_out:
   input: expand(DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dappeak1.out',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])
rule gat_mc_dap_08_out_clean:
    shell: "rm {rules.gat_mc_dap_08_out.input}"

rule gat_mc_dap_09:# dmC bins (after trimming error was fixed May-2016) and DAP TF motifs not including ampDAP, background is TAIR10 100bp uniquely mappable
   input: mc_bed=DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation.bed',dap_peak=dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3')],wrk=GENOMES_TAIR10_MAP_UNIQ100,num_samples=10000,D=DMC_BINS_GAT02+'/{mc_prefix}'
   output: DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dapfimo1.out'
   log: DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dapfimo1.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} {params.annot} --log={log} > {output}
   """
rule gat_mc_dap_09_out: 
   input: expand(DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dapfimo1.out',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])
rule gat_mc_dap_09_out_clean:
   shell: "rm {rules.gat_mc_dap_09_out.input}"

rule gat_mc_dap_10:# dmC bins (after trimming error was fixed May-2016) and DAP TF motifs not including ampDAP, background is clustered DAP TF motifs not including ampDAP
   input: mc_bed=DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation.bed'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',num_samples=10000,num_threads=16,D=DMC_BINS_GAT02+'/{mc_prefix}'
   output: DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dapfimo2.out'
   log: DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dapfimo2.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --num-threads={params.num_threads} --output-counts-pattern={params.D}/{wildcards.mc_prefix}.%s.gat_dapfimo2.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_mc_dap_10_out: 
   input: expand(DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dapfimo2.out',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])
rule gat_mc_dap_10_out_clean:
   shell: "rm {rules.gat_mc_dap_10_out.input}"

rule gat_mc_dap_11:# dmC bins eQTL (after trimming error was fixed May-2016) and DAP TF motifs not including ampDAP, background is clustered DAP TF motifs not including ampDAP
   input: mc_bed=MARGINAL_TEST_COV_22+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genesbonf_5e-2.bed.sorted.uniq'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',num_samples=10000,D=DMC_BINS_GAT02+'/{mc_prefix}'
   output: DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dapfimo4.out'
   log: DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dapfimo4.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --output-counts-pattern={params.D}/{wildcards.mc_prefix}.%s.gat_dapfimo4.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_mc_dap_11_out: 
   input: expand(DMC_BINS_GAT02+'/{mc_prefix}/{mc_prefix}.gat_dapfimo4.out',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])
rule gat_mc_dap_11_out_clean:
   shell: "rm {rules.gat_mc_dap_11_out.input}"
   
rule gat_ms_dap_04:# dmC sites and DAP TF motifs not including ampDAP, background is clustered DAP TF motifs not including ampDAP
   input: mc_bed=DMS+'/{mc_prefix}/{mc_prefix}_1001g_filter1.bed'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',num_samples=10000,D=DMS_GAT01+'/{mc_prefix}'
   output: DMS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo2.out'
   log: DMS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo2.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --output-counts-pattern={params.D}/{wildcards.mc_prefix}.%s.gat_dapfimo2.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_ms_dap_04_out: 
   input: expand(DMS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo2.out',mc_prefix=DMS_TABLES_PREFIX)
rule gat_ms_dap_04_out_clean:
   shell: "rm {rules.gat_ms_dap_04_out.input}"

rule gat_ms_dap_06:# dmC eQTL sites and DAP TF motifs not including ampDAP, background is clustered DAP TF motifs not including ampDAP
   input: mc_bed=MARGINAL_TEST_COV_18+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genesbonf_5e-2.bed.sorted.uniq'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',num_samples=10000,D=DMS_GAT01+'/{mc_prefix}'
   output: DMS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo3.out'
   log: DMS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo3.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --output-counts-pattern={params.D}/{wildcards.mc_prefix}.%s.gat_dapfimo3.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_ms_dap_06_out: 
   input: expand(DMS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo3.out',mc_prefix=DMS_TABLES_PREFIX)
rule gat_ms_dap_06_out_clean:
   shell: "rm {rules.gat_ms_dap_06_out.input}"

rule gat_ms_dap_07:# dmC eQTL sites and DAP TF motifs not including ampDAP, background is clustered DAP TF motifs not including ampDAP
   input: mc_bed=MARGINAL_TEST_COV_21+'/{mc_prefix}_chr1-5/results/lmm_master/chr1-5genesbonf_5e-2.bed.sorted.uniq'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',num_samples=10000,D=DMS_GAT01+'/{mc_prefix}'
   output: DMS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo4.out'
   log: DMS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo4.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mc_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --output-counts-pattern={params.D}/{wildcards.mc_prefix}.%s.gat_dapfimo4.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_ms_dap_07_out: 
   input: expand(DMS_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo4.out',mc_prefix=DMS_TABLES_PREFIX)
rule gat_ms_dap_07_out_clean:
   shell: "rm {rules.gat_ms_dap_07_out.input}"

rule gat_mqtl_dap_01:# mQTL with dap peaks
   input: mqtl_bed=MQTL+'/{mc_prefix}/mQTL.bed',dap_peak=dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=GENOMES_TAIR10_MAP_UNIQ100,num_samples=10000,D=MQTL_GAT01+'/{mc_prefix}'
   output: MQTL_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dappeak1.out'
   log: MQTL_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dappeak1.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mqtl_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} {params.annot} --log={log} > {output}
   """
rule gat_mqtl_dap_01_out: 
   input: expand(MQTL_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dappeak1.out',mc_prefix=['C_DMR','CG_DMR'])
rule gat_mqtl_dap_01_out_clean:
   shell: "rm {rules.gat_mqtl_dap_01_out.input}"

rule gat_mqtl_dap_02:# dap motif matches
   input: mqtl_bed=MQTL+'/{mc_prefix}/mQTL.bed',dap_peak=dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/encode_cluster_01.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',DAP_V4_FAMCLUST_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3')],wrk=GENOMES_TAIR10_MAP_UNIQ100,num_samples=10000,D=MQTL_GAT01+'/{mc_prefix}'
   output: MQTL_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo1.out'
   log: MQTL_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo1.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.mqtl_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} {params.annot} --log={log} > {output}
   """
rule gat_mqtl_dap_02_out: 
   input: expand(MQTL_GAT01+'/{mc_prefix}/{mc_prefix}.gat_dapfimo1.out',mc_prefix=['C_DMR','CG_DMR'])
rule gat_mqtl_dap_02_out_clean:
   shell: "rm {rules.gat_mqtl_dap_02_out.input}"

rule gat_snp_dap_04:# all SNP and DAP TF motifs not including ampDAP, background is clustered DAP TF motifs not including ampDAP
   input: snp_bed=MARGINAL_TEST_COV_16+'/geno_pos/position.bed'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',num_samples=10000,D=SNP_GAT01+'/snp/'
   output: SNP_GAT01+'/snp/snp.gat_dapfimo2.out'
   log: SNP_GAT01+'/snp/snp.gat_dapfimo2.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.snp_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --output-counts-pattern={params.D}/snp.%s.gat_dapfimo2.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_snp_dap_04_out: 
   input: SNP_GAT01+'/snp/snp.gat_dapfimo2.out'
rule gat_snp_dap_04_out_clean:
   shell: "rm {rules.gat_snp_dap_04_out.input}"

rule gat_snp_dap_06:# eQTL sites and DAP TF motifs not including ampDAP, background is clustered DAP TF motifs not including ampDAP
   input: snp_bed=MARGINAL_TEST_COV_16+'/results/lmm_master/chr1-5genesbonf_5e-2.bed.sorted.uniq'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',num_samples=10000,D=SNP_GAT01+'/snp'
   output: SNP_GAT01+'/snp/snp.gat_dapfimo3.out'
   log: SNP_GAT01+'/snp/snp.gat_dapfimo3.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.snp_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --output-counts-pattern={params.D}/snp.%s.gat_dapfimo3.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_snp_dap_06_out: 
   input: SNP_GAT01+'/snp/snp.gat_dapfimo3.out'
rule gat_snp_dap_06_out_clean:
   shell: "rm {rules.gat_snp_dap_06_out.input}"

rule gat_snp_dap_07:# eQTL sites and DAP TF motifs not including ampDAP, background is clustered DAP TF motifs not including ampDAP
   input: snp_bed=MARGINAL_TEST_COV_19+'/results/lmm_master/chr1-5genesbonf_5e-2.bed.sorted.uniq'
   params: annot=['--annotation-file=%s'%a for a in dap_fam_all_dep(DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed',DAP_V4_REPMASTER_01+'/{family_expSys}/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed')],wrk=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',num_samples=10000,D=SNP_GAT01+'/snp'
   output: SNP_GAT01+'/snp/snp.gat_dapfimo4.out'
   log: SNP_GAT01+'/snp/snp.gat_dapfimo4.log'
   shell: """
        module load shhuang anaconda/2.4.1
        mkdir -p {params.D}
        gat-run.py --segment-file={input.snp_bed} --workspace-file={params.wrk} --ignore-segment-tracks --num-samples={params.num_samples} --output-counts-pattern={params.D}/snp.%s.gat_dapfimo4.counts.tsv.gz {params.annot} --log={log} > {output}
   """
rule gat_snp_dap_07_out: 
   input: SNP_GAT01+'/snp/snp.gat_dapfimo4.out'
rule gat_snp_dap_07_out_clean:
   shell: "rm {rules.gat_snp_dap_07_out.input}"


##################
# eQTL and mQTL
##################
rule qtl_multiinter_01:
    input: c_dmr=os.path.join(MQTL,'C_DMR/mQTL.bed'),cg_dmr=os.path.join(MQTL,'CG_DMR/mQTL.bed'),expr=MARGINAL_TEST_COV_16+'/results/lmm_master/chr1-5genesbonf_5e-2_merge10k.bed'
    output: PROJ_RESULTS_PATH_GALE+'/qtl_multiinter01/qtl_intervals_table.txt'
    shell: """
        module load bedtools/2.25.0
        bedtools multiinter -header -names C_DMR CG_DMR RNA -i {input.c_dmr} {input.cg_dmr} {input.expr} > {output}
    """
rule qtl_multiinter_01_clean:
    shell: "rm {rules.qtl_multiinter_01.output}"

rule qtl_multiinter_02:
    input: c_dmr=os.path.join(MQTL,'C_DMR/mQTL.bed'),cg_dmr=os.path.join(MQTL,'CG_DMR/mQTL.bed'),expr=MARGINAL_TEST_COV_19+'/results/lmm_master/chr1-5genesbonf_5e-2_merge10k.bed'
    output: PROJ_RESULTS_PATH_GALE+'/qtl_multiinter02/qtl_intervals_table.txt'
    shell: """
        module load bedtools/2.25.0
        bedtools multiinter -header -names C_DMR CG_DMR RNA -i {input.c_dmr} {input.cg_dmr} {input.expr} > {output}
    """
rule qtl_multiinter_02_clean:
    shell: "rm {rules.qtl_multiinter_02.output}"


#####################################
# cistrome and epicistrome overlap
#####################################
rule cistrome_venn_01:
    input: mc_bed=DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation.bed',cistr_bed=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',epicistr_bed=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.ampMinusDAP_in_MR.fimo0_1e-4.bed3'
    output: cistr_table=CISTROME_OLAP_01+'/{mc_prefix}/{mc_prefix}_cistrome_table.txt'
    run:#module load shhuang python_modonly/3.4.3 numpy/1.10.1 pandas/0.17.1 pysam/0.9.0 pybedtools/0.7.6
       import sys
       import pybedtools
       import pybedtools.contrib.venn_maker as venn_maker
       
       beds = [input.mc_bed,input.cistr_bed,input.epicistr_bed]
       names=[wildcards.mc_prefix,"cistrome","epicistrome"]
       _beds = []
       for bed in beds:
           if not isinstance(bed, pybedtools.BedTool):
               bed = pybedtools.BedTool(bed)
           _beds.append(bed)
       cleaned = venn_maker.cleaned_intersect(_beds)
       results = OrderedDict(list(zip(names, cleaned)))
                                           
       ofh = open(output.cistr_table,'w')
       for key,val in list(results.items()):
           for i in val:
               if isinstance(i,pybedtools.Interval):
                   i = str(i).replace('\t','|')
               ofh.write('%s\t%s'%(key,i))
       ofh.close()
                         
rule cistrome_venn_01_out:
    input: expand(CISTROME_OLAP_01+'/{mc_prefix}/{mc_prefix}_cistrome_table.txt',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])

rule cistrome_venn_02:# dmC_bins after trimming error was fixed May 2016
    input: mc_bed=expand(DMC_BINS+'/{mc_prefix}/{mc_prefix}_methylation.bed',mc_prefix=['dmC_filtered','dmCH_filtered']),cistr_bed=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',epicistr_bed=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.ampMinusDAP_in_MR.fimo0_1e-4.bed3'
    output: cistr_table=CISTROME_OLAP_02 + '/dmC_all/dmC_all_cistrome_table.txt'
    run:#module load shhuang python_modonly/3.4.3 numpy/1.10.1 pandas/0.17.1 pysam/0.9.0 pybedtools/0.7.6
       import sys
       import pybedtools
       import pybedtools.contrib.venn_maker as venn_maker

       beds = input.mc_bed + [input.cistr_bed,input.epicistr_bed]
       names=["dmC_filtered","dmCH_filtered","cistrome","epicistrome"]
       _beds = []
       for bed in beds:
           if not isinstance(bed, pybedtools.BedTool):
               bed = pybedtools.BedTool(bed)
           _beds.append(bed)
       cleaned = venn_maker.cleaned_intersect(_beds)
       results = OrderedDict(list(zip(names, cleaned)))
                                           
       ofh = open(output.cistr_table,'w')
       for key,val in list(results.items()):
           for i in val:
               if isinstance(i,pybedtools.Interval):
                   i = str(i).replace('\t','|')
               ofh.write('%s\t%s'%(key,i))
       ofh.close()

rule cistrome_venn_03:
    input: mc_bed=DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation.bed',cistr_bed=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',epicistr_bed=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.ampMinusDAP_in_MR.fimo0_1e-4.bed3'
    output: cistr_table=CISTROME_OLAP_02+'/{mc_prefix}/{mc_prefix}_cistrome_table.txt'
    run:#module load shhuang python_modonly/3.4.3 numpy/1.10.1 pandas/0.17.1 pysam/0.9.0 pybedtools/0.7.6
       import sys
       import pybedtools
       import pybedtools.contrib.venn_maker as venn_maker
       
       beds = [input.mc_bed,input.cistr_bed,input.epicistr_bed]
       names=[wildcards.mc_prefix,"cistrome","epicistrome"]
       _beds = []
       for bed in beds:
           if not isinstance(bed, pybedtools.BedTool):
               bed = pybedtools.BedTool(bed)
           _beds.append(bed)
       cleaned = venn_maker.cleaned_intersect(_beds)
       results = OrderedDict(list(zip(names, cleaned)))
                                           
       ofh = open(output.cistr_table,'w')
       for key,val in list(results.items()):
           for i in val:
               if isinstance(i,pybedtools.Interval):
                   i = str(i).replace('\t','|')
               ofh.write('%s\t%s'%(key,i))
       ofh.close()
                         
rule cistrome_venn_03_out:
    input: expand(CISTROME_OLAP_02+'/{mc_prefix}/{mc_prefix}_cistrome_table.txt',mc_prefix=['dmC_filtered','dmCG_filtered','dmCH_filtered'])

rule cistrome_venn_04:
    input: mc_bed=expand(DMC_BINS_NEW+'/{mc_prefix}/{mc_prefix}_methylation.bed',mc_prefix=['dmC_filtered','dmCH_filtered']),cistr_bed=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed3',epicistr_bed=DAP_V4_FAMCLUST_01+'/master/chr1-5_GEM_events.ampMinusDAP_in_MR.fimo0_1e-4.bed3'
    output: cistr_table=CISTROME_OLAP_02 + '/dmC_all/dmC_all_cistrome_table.txt'
    run:#module load shhuang python_modonly/3.4.3 numpy/1.10.1 pandas/0.17.1 pysam/0.9.0 pybedtools/0.7.6
       import sys
       import pybedtools
       import pybedtools.contrib.venn_maker as venn_maker

       beds = input.mc_bed + [input.cistr_bed,input.epicistr_bed]
       names=["dmC_filtered","dmCH_filtered","cistrome","epicistrome"]
       _beds = []
       for bed in beds:
           if not isinstance(bed, pybedtools.BedTool):
               bed = pybedtools.BedTool(bed)
           _beds.append(bed)
       cleaned = venn_maker.cleaned_intersect(_beds)
       results = OrderedDict(list(zip(names, cleaned)))
                                           
       ofh = open(output.cistr_table,'w')
       for key,val in list(results.items()):
           for i in val:
               if isinstance(i,pybedtools.Interval):
                   i = str(i).replace('\t','|')
               ofh.write('%s\t%s'%(key,i))
       ofh.close()

