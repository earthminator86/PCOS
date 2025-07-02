run_scPagwas_main_run <- function(
    prop_cell_used = 0.3,
    pathway_included_file = "immune_kegg_pathways.txt",
    seurat_rds = "Yazar2022Science_Seurat_sub_scale.rds",
    gwas_prune_file = "data/PCOS_sub_rsid_LD0.8_prune_gwas_data2.txt",
    out_dir="GWAS_scPagwas/"
){
  sc_study_name <- str_remove(basename(seurat_rds), "_Seurat_sub_scale.rds")
  
  library(SeuratObject)
  library(Seurat)
  library(scPagwas)
  
  pathways_included <- fread(pathway_included_file)
  Genes_by_pathway_kegg_included <- Genes_by_pathway_kegg[pathways_included$Entry]
  
  setwd("GWAS_scPagwas/data")
  Pagwas_data<-scPagwas_main(Pagwas = NULL,
                             gwas_data = gwas_prune_file,
                             Single_data = seurat_rds,
                             output.prefix= sc_study_name, # the prefix name for output files
                             output.dirs="main_run",
                             block_annotation = block_annotation,
                             assay="RNA", 
                             Pathway_list = Genes_by_pathway_kegg_included,
                             n.cores=1,
                             iters_singlecell = 10,
                             chrom_ld = chrom_ld,
                             singlecell=T, 
                             celltype=T
  )
  save(Pagwas_data, file = paste0("./main_run/", sc_study_name,"_scPagwas.RData") )
}