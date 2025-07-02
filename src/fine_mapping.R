rm(list=ls())
library(data.table)
library(stringr)
library(parallel)

ldsc_munge_sumstat <- function(combined_GWAS_outfile="PCOS_txt",
                               gwas_loci_file = "/GWAS_loci/PCOS_locis.txt",
                               python_path = "anaconda3/envs/ldsc/bin/python",
                               ldsc_exe_dir = "ldsc/",
                               snp_list_file = "list.txt",
                               out_dir = "GWAS_ldsc/",
                               out_file_suffix="PCOS" ){
  
  gwas_sub_file <- paste0(out_dir, "data/gwas_sub.txt")
  gwas_sumstat_file <- paste0(out_dir, "data/",out_file_suffix)
  if(!dir.exists(dirname(gwas_sub_file))){ dir.create(dirname(gwas_sub_file), recursive = TRUE) }
  
  snp_list <- fread(snp_list_file, header=FALSE); names(snp_list) <- "rsid"; setindex(snp_list, "rsid")
  gwas_i <- fread(combined_GWAS_outfile, header = TRUE)
  if(!str_detect(gwas_i[1,]$chromosome, "chr")){ gwas_i$chromosome <- paste0("chr", gwas_i$chromosome) }
  gwas_i <- gwas_i[,.(rsid=variant_id, Allele1=effect_allele, Allele2=non_effect_allele, Freq=EAF, p=pvalue, N=sample_size, BETA=effect_size)]
  setindex(gwas_i, rsid)
  fwrite(gwas_i[str_detect(rsid, "^rs")], gwas_sub_file, sep="\t")
  system(paste0(python_path, " ", ldsc_exe_dir, "munge_sumstats.py --sumstats ", gwas_sub_file," --out ", gwas_sumstat_file," --frq Freq --N-col N --a1 Allele1 --a2 Allele2 --snp rsid --p p" ))
  
}

run_SusieX_loci <- function(gwas_loci_file = "GWAS_loci/PCOS_locis.txt",
                            combined_GWAS_outfile="PCOS.txt",
                            gwas_sumstat_file = "GWAS_ldsc/data/PCOS.sumstats.gz",
                            SuSiEx_dir = "tools/SuSiEx/",
                            python_2="/usr/bin/python2",
                            extend_loci_length = 5e5,
                            plink_1KG_split_dir = "ref/LDSCfiles/plink_files/split_by_chrom_sorted/",
                            out_dir = "GWAS_fine_mapping/",
                            out_file_suffix="PCOS",
                            parThreads = 10,
                            
){
  # gwas_info:
  if(!file.exists(combined_GWAS_outfile)){ return()}
  if(!file.exists(gwas_sumstat_file)){ return()}
  
  # get all loci for each gwas:
  all_loci <- fread(gwas_loci_file)
  
  gwas_i <- fread(combined_GWAS_outfile, header = TRUE)
  if(str_detect(gwas_i[1,]$chromosome, "chr")){ gwas_i[,c("chromosome"):=.(str_remove(chromosome, "chr"))] }
  gwas_i <- gwas_i[chromosome %in% 1:22,.(chr= as.integer(chromosome), snp=variant_id, bp=position, A1=effect_allele, A2=non_effect_allele, beta=effect_size, se=standard_error, stat=qt(pvalue, sample_size, lower=FALSE)* (ifelse(effect_size<0,-1,1)), p=pvalue)]
  setindex(gwas_i, chr, bp)
  
  # for each loci:
  loci_names <- unique(all_loci$loci_name)
  a <- mclapply(1:length(loci_names), function(j){
    message(j,"/",length(loci_names))
    loci_j <- loci_names[j]
    loci_chr_j <- str_remove(str_split(loci_j, "-")[[1]][1], "chr")
    loci_start_j <- as.numeric(str_split(loci_j, "-")[[1]][2])
    loci_end_j <- as.numeric(str_split(loci_j, "-")[[1]][3])
    
    gwas_loci_j <- na.omit(gwas_i[chr == loci_chr_j & bp >= loci_start_j & bp <= loci_end_j][order(bp)])
    gwas_sub_j <- paste0(out_dir, "data/locis_",out_file_suffix,"/",loci_j,".sumstats.txt")
    ld_suffix_j <- paste0(out_dir, "data/locis_",out_file_suffix,"/",loci_j)
    if(!dir.exists(dirname(ld_suffix_j))){ dir.create(dirname(ld_suffix_j), recursive = TRUE) }
    fwrite(na.omit(gwas_loci_j), gwas_sub_j, sep="\t")
    # calculate loci for this loci:
    system(paste0(python_2, " ",SuSiEx_dir,"utilities/SuSiEx_LD.py --ref_file=", plink_1KG_split_dir, "1000G.EUR.hg38.",loci_chr_j,
                  " --ld_file=", ld_suffix_j,
                  " --chr=", loci_chr_j,
                  " --bp=", loci_start_j, ",", loci_end_j,
                  " --plink utilities/plink",
                  " --maf 0.005"
    ))
    #
    if(!dir.exists(paste0(out_dir, "data/SusieX/"))){ dir.create(paste0(out_dir, "data/SusieX/"), recursive = TRUE) }
    system(paste0(SuSiEx_dir, "bin/SuSiEx ",
                  " --sst_file=",gwas_sub_j,
                  " --n_gwas=", gwas_sample_size,
                  # " --ref_file=", ld_suffix_j, 
                  " --ld_file=", ld_suffix_j,
                  " --out_dir=",out_dir, "data/SusieX/",
                  " --out_name=", out_file_suffix,"_loci_",loci_j,".cs95",
                  " --level=0.95",
                  " --pval_thresh=1e-5",
                  " --maf=0.005",
                  " --chr=", loci_chr_j,
                  " --bp=", loci_start_j, ",", loci_end_j,
                  " --snp_col=2",
                  " --chr_col=1",
                  " --bp_col=3",
                  " --a1_col=4",
                  " --a2_col=5",
                  " --eff_col=6",
                  " --se_col=7",
                  " --pval_col=9",
                  " --plink=",SuSiEx_dir,"utilities/plink",
                  " --mult-step=True",
                  " --keep-ambig=True",
                  " --threads=5"
    ))
    return()
  }, mc.cores= 10)
  
}