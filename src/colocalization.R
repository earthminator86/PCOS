coloc_analysis_sc <- function(
    qtl_splited_dir = "qtl_splited/",
    combined_GWAS_outfile="PCOS.txt.gz",
    gwas_loci_file = "PCOS_locis.txt",
    coloc_out_dir = "GWAS_coloc/",
    gene_info_file = "genome/gencode_v38_gene_info.txt",
    parThreads=20
){
  
  gene_info <- fread(gene_info_file)
  
  gwas_locis <- fread(gwas_loci_file) 
  gwas_locis <- gwas_locis[pvalue<5e-8]
  library(valr)
  gwas_locis_genes <- as.data.table(valr::bed_closest(gwas_locis[,.(chrom = chromosome, start=position, end=position+1, variant_id, loci_name)], 
                                                      gene_info[,.(chrom= chr, start, end, gene_name)]))
  
  # load gwas:
  gwas_i <- fread(combined_GWAS_outfile)
  setindex(gwas_i, chromosome, variant_id)
  
  # sc-eQTLs:
  cts <- data.table(dir_path = normalizePath(list.files(qtl_splited_dir, full.names = TRUE)))
  cts[,c("ct"):=.(basename(dir_path))]
  library(xQTLbiolinks)
  coloc_res <- rbindlist(lapply(1:nrow(gwas_locis_genes), function(i){
    message("==> ", i, "/", nrow(gwas_locis_genes))
    loci_i <- gwas_locis_genes[i,]
    loci_i <- as.data.table(tidyr::separate(loci_i, col="loci_name", sep = "-", into = c("chrom", "start", "end"), remove = FALSE, convert = TRUE))
    
    coloc_res_i <- rbindlist(mclapply(1:nrow(cts), function(j){
      qtl_file_j <- paste0(cts[j,]$dir_path, "/",loci_i$gene_name, ".txt")
      if( !file.exists(qtl_file_j) ){ return() }
      qtl_j <- fread(qtl_file_j)
      if(nrow(qtl_j[chromosome== loci_i$chrom & position >= loci_i$start & position <= loci_i$end])==0){ return() }
      #
      shared_snps_j <- intersect(gwas_i[chromosome==qtl_j[1,]$chromosome,]$variant_id, qtl_j$variant_id)
      gwas_j <- gwas_i[chromosome==qtl_j[1,]$chromosome ][variant_id %in% shared_snps_j, ][,.(rsid=variant_id, chrom=chromosome, position, pvalue, maf=ifelse(EAF<0.5, EAF, 1-EAF), beta=effect_size, se=standard_error)]
      qtl_j <- qtl_j[variant_id %in% shared_snps_j,.(rsid=variant_id, chrom=chromosome, position, pvalue, maf=ifelse(EAF<0.5, EAF, 1-EAF), beta=effect_size, se=0.1)]
      
      coloc_res_j <- xQTLanalyze_coloc_diy(gwas_j, qtl_j)$coloc_Out_summary
      coloc_res_j[,c("ct"):=.(cts[j,]$ct)]
      coloc_res_j <- cbind(coloc_res_j, loci_i[,.(loci_name, gene_name)])
      return(coloc_res_j)
    }, mc.cores = parThreads))
    return(coloc_res_i)
  }))
  coloc_res[PP.H4.abf>0.75]
  if(!dir.exists(paste0(coloc_out_dir, "data/"))){ dir.create(paste0(coloc_out_dir, "data/"), recursive = TRUE) }
  fwrite(coloc_res, paste0(coloc_out_dir, "data/sceQTL_coloc_res.txt"), sep="\t" )
}

coloc_analysis_GTEx <- function(
    combined_GWAS_outfile="PCOS.txt.gz",
    gwas_loci_file = "GWAS_loci/PCOS_locis.txt",
    coloc_out_dir = "GWAS_coloc/",
    gene_info_file = "genome/gencode.v26_gene_info.txt",
    parThreads= 10
){
  if(!dir.exists(paste0(coloc_out_dir, "data/"))){ dir.create(paste0(coloc_out_dir, "data/"), recursive = TRUE) }
  
  GTEx_tissues = xQTLbiolinks::tissueSiteDetailGTExv8
  gene_info <- fread(gene_info_file)
  
  # load loci:
  gwas_locis <- fread(gwas_loci_file) 
  gwas_locis <- gwas_locis[pvalue<5e-8]
  library(valr)
  gwas_locis_genes <- as.data.table(valr::bed_closest(gwas_locis[,.(chrom = chromosome, start=position, end=position+1, variant_id, loci_name)], 
                                                      gene_info[,.(chrom= chr, start, end, gene_name, gene_id)]))
  
  # load gwas:
  gwas_i <- fread(combined_GWAS_outfile)
  setindex(gwas_i, chromosome, variant_id)
  
  library(xQTLbiolinks)
  coloc_res <- rbindlist(lapply(1:nrow(GTEx_tissues), function(t){
    tissue_t <- GTEx_tissues[t,]$tissueSiteDetail
    message(" ==> t: ",t, "/",nrow(GTEx_tissues), " | ", tissue_t, " | ", date())
    coloc_res_t <- rbindlist(mclapply(1:nrow(gwas_locis_genes), function(i){
      gene_i <- gwas_locis_genes[i]$gene_id
      gtex_qtl_i <- xQTLbiolinks::xQTLdownload_eqtlAllAsso(gene_i, tissueLabel = tissue_t, data_source = "liLab")
      if(is.null(gtex_qtl_i)){ return() }
      if(nrow(gtex_qtl_i)<10){ return()}
      message("==> i: ",i, "/", nrow(gwas_locis_genes), " | ", gene_i)
      
      #
      shared_snps_j <- intersect(gwas_i[chromosome == gwas_locis_genes[i,]$chrom]$variant_id, gtex_qtl_i$rsid)
      gwas_j <- gwas_i[chromosome==gwas_locis_genes[i,]$chrom ][variant_id %in% shared_snps_j, ][,.(rsid=variant_id, chrom=chromosome, position, pvalue, maf=ifelse(EAF<0.5, EAF, 1-EAF), beta=effect_size, se=standard_error)]
      gtex_qtl_i <- merge(gtex_qtl_i, gwas_j[,.(rsid, position)], by="rsid", sort=FALSE)
      gtex_qtl_i <- gtex_qtl_i[rsid %in% shared_snps_j,.(rsid, chrom=gwas_locis_genes[i,]$chrom, position, pvalue=pValue, maf, beta, se)]
      
      coloc_res_i <- xQTLanalyze_coloc_diy(gwas_j, gtex_qtl_i)$coloc_Out_summary
      coloc_res_i[,c("ct"):=.(tissue_t)]
      coloc_res_i <- cbind(coloc_res_i, gwas_locis_genes[i,.(loci_name, gene_name)])
      return(coloc_res_i)
    }, mc.cores = parThreads))
    return(coloc_res_t)
  }))
  
  coloc_res[PP.H4.abf>0.75]
  fwrite(coloc_res, paste0(coloc_out_dir, "data/GTEx_eQTL_coloc_res.txt"), sep="\t" )
}
