rm(list=ls())
library(data.table)
library(stringr)
library(parallel)

get_loci_run <- function(
    GWAS_outfile="PCOS.txt",
    out_dir = "GWAS_loci/",
    ld_file = "g1000_eur",
    p_sig_cutoff = 5e-8,
    clump_kb=250,
    clump_r2=0.3,
    plink_path = "plink",
    out_file_suffix="PCOS" 
){
  tmp_dir <- paste0(out_dir, out_file_suffix, "_tmp/");
  loci_dt <- data.table()
  loci_clumped_file_i <- paste0(tmp_dir, out_file_suffix, ".clumped");
  if(!file.exists(loci_clumped_file_i)){ 
    system(paste0(plink_path, " --bfile ",ld_file,
                  " --clump ", GWAS_outfile,
                  " --clump-kb ", clump_kb," --clump-r2 ", clump_r2," --clump-p1 ", p_sig_cutoff,
                  " --clump-snp-field SNP --clump-field P --clump-allow-overlap",
                  " --out ",tmp_dir, out_file_suffix))
  }
  if(!file.exists(loci_clumped_file_i)){ next() }
  loci_clumped_i <- fread(loci_clumped_file_i)
  loci_clumped_i <- loci_clumped_i[P < p_sig_cutoff]
  if(nrow(loci_clumped_i)==0){next()}
  loci_clumped_i <- as.data.table(tidyr::separate_longer_delim(loci_clumped_i, cols = "SP2", delim = ","))
  loci_clumped_i[,c("SP2"):=.(unlist(lapply(loci_clumped_i$SP2, function(x){ str_remove(x, regex("\\(.*\\)")) })))]
  loci_clumped_snps_i <- c(unique(loci_clumped_i$SNP), loci_clumped_i$SP2)
  loci_clumped_snps_i <- loci_clumped_snps_i[loci_clumped_snps_i!="NONE"]
  #
  gwas_dt_tmp <- fread(GWAS_outfile)
  loci_dt_i <- gwas_dt_tmp[variant_id %in% loci_clumped_snps_i][order(pvalue)]
  loci_name_i <- paste0(gwas_sig_i[i,]$chromosome, "-", min(loci_dt_i$position), "-", max(loci_dt_i$position))
  loci_dt_i[,c("loci_name") := .(loci_name_i)]
  loci_dt <- rbind(loci_dt, loci_dt_i)
  gwas_dt_tmp <- rbind(gwas_dt_tmp[chromosome !=gwas_sig_i[i,]$chromosome],
                       gwas_dt_tmp[chromosome ==gwas_sig_i[i,]$chromosome][ position < min(loci_dt_i$position) | position > max(loci_dt_i$position) ])
  setindex(gwas_dt_tmp, variant_id)
  message("==> i: ",i, " | ", snp_i, " | loci snp num: ",nrow(loci_dt_i));
  fwrite(loci_dt, paste0(out_dir, out_file_suffix,"_locis_all.txt"), sep="\t" )
}