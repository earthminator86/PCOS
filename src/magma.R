rm(list=ls())
library(data.table)
library(stringr)
library(parallel)

run_magma <- function(
    GWAS_outfile="PCOS_txt",
    plink_1000G_dir = "magma_bfile_rsid/",
    magma_dir = "tools/magma/",
    msigdb_file = "ref/Hs.entrez.gmt",
    out_dir ="GWAS_magma/"
){
  if(!dir.exists(paste0(out_dir, "data/"))){ dir.create(paste0(out_dir, "data/"), recursive = TRUE) }
  
  gwas_i <- fread(GWAS_outfile)
  gwas_i <- gwas_i[,.( CHR=str_remove(chromosome,"chr"), BP=position, SNP=variant_id, P=pvalue, N=sample_size )]
  fwrite(gwas_i, paste0(out_dir, "data/gwas_for_magma.txt"), sep = "\t" )
  
  paste0(magma_dir, "src/magma --annotate window=10,10 ",
         " --snp-loc ", plink_1000G_dir, "g1000_eur.bim ",
         " --gene-loc ", magma_dir, "Gene_locations/NCBI38.gene.loc ",
         " --out ", out_dir, "data/magma_anno")
  
  paste0(magma_dir, "src/magma --bfile ",  plink_1000G_dir, "g1000_eur ",
         " --pval ", out_dir, "data/gwas_for_magma.txt use=3,4 ncol=5 ",
         " --gene-annot ",out_dir, "data/magma_anno.genes.annot ",
         " --out ", out_dir, "data/magma_gene_P"
  )
  
  paste0(magma_dir, "src/magma ", 
         " --gene-results ", out_dir, "data/magma_gene_P.genes.raw",
         " --set-annot ", msigdb_file,
         " --out ", out_dir, "data/magma_MSigDB_c8_all")
}
