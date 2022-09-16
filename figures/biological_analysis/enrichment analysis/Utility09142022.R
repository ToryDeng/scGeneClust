setwd('./scGeneClust/figures/biological_analysis/enrichment analysis')

library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(AnnotationDbi)
library(EGAD)
library(reshape2)
library(FEAST)
library(GO.db)
library(patchwork)
library(KEGGREST)
library(tidyverse)


filter_single_g <- function(cluster_df) {
  #filter out single gene cluster: delete some rows in the dataframe
  
  #param cluster_df:original dataframe of the total gene cluster result
  
  agg<-aggregate(cluster_df$cluster, by=list(cluster=cluster_df$cluster),length)
  clist<-agg[which(agg$x>1),'cluster']
  filter_df<-cluster_df[which(cluster_df$cluster%in%clist),]
}


convert_symbol2id <- function(filter_df, inner = FALSE) {
  #convert SYMBOL to ENTREZID: add a new column (ENTREZID) in the dataframe
  
  #param filter_df: return of filter_single_g, 
  #a dataframe of gene cluster result with 4 columns after filtering out single gene cluster
  #param inner: BOLL value indicated whether drop out genes whose IDs are fail to map
  #(default FALSE) leave NA values in the new column (ENTREZID)
  
  eg <- bitr(filter_df$original_gene,
             fromType="SYMBOL",
             toType=c("ENTREZID"),
             OrgDb="org.Hs.eg.db"
  )
  ENTREZID <- rep("NA",time = nrow(filter_df))
  ENTREZID[match(eg$SYMBOL,filter_df$original_gene)] <- eg$ENTREZID
  filter_id_df <- cbind(ENTREZID,filter_df)
  if (inner == TRUE) {
    filter_id_df <- subset(filter_id_df, ENTREZID!="NA")
  }
  filter_id_df
}


calculate_cluster_score <- function(filter_id_df, cluster_score_method) {
  #calculate cluster-level scores: add a new column (cluster_score) in the dataframe
  
  #param filter_id_df: return of convert_symbol2id,
  #a dataframe of gene cluster result with 5 columns after adding a column of ENTREZID
  #param cluster_score_method: ('seurat','relevant') how to obtain cluster-level scores
  #'seurat': for each cluster, choosing the seurat score of the most representative gene whose score(i.e. centrality) is the highest within its cluster
  #'relevant': for each cluster, choosing the max mutual information
  
  clist<-unique(filter_id_df$cluster)
  if (cluster_score_method=='seurat') {
    for (i in clist) {
      cluster_i<-filter_id_df[which(filter_id_df$cluster==i),]
      seurat_score<-mean(cluster_i[order(-cluster_i$score),][['variances_norm']][1])
      filter_id_df[which(filter_id_df$cluster==i),'cluster_score']<-seurat_score
    }
  }
  else {
    #cluster_score_method=='relevant'
    for (i in clist) {
      cluster_i<-filter_id_df[which(filter_id_df$cluster==i),]
      relevance<-max(cluster_i$score)
      filter_id_df[which(filter_id_df$cluster==i),'cluster_score']<-relevance
    }
  }
  filter_id_df
}

rank_cluster <- function(filter_id_df) {
  #group and rank descendingly according to cluster_score
  #param filter_id_df: return of calculate_cluster_score, 
  #a dataframe of gene cluster result with 6 columns after adding a column of cluster_score
  
  filter_id_df_grouped <- group_by(
    .data = filter_id_df[,c("cluster","cluster_score")]
    ,cluster)	
  cluster_group <- summarise(.data = filter_id_df_grouped
                             ,cluster_score = max(cluster_score, na.rm = TRUE)
                             ,Count_genes = n()
  )
  cluster_rank<-arrange(cluster_group,-cluster_score)
  cluster_rank
}

select_genes <- function(filter_id_df, cluster_rank, select_gene_method, threshold, cluster_df_all) {
  #select experimental group's genes 
  #and stochastically select control group's genes with the same as the experimental group size
  
  #param filter_id_df: return of calculate_cluster_score, 
  #a dataframe of gene cluster result with 6 columns after adding a column of cluster_score
  #param cluster_rank: return of rank_cluster, 
  #a dataframe of gene cluster ranking result with 3 columns after grouping
  #param select_gene_method: ('top_cluster','top_gene') how to select the number of cluster
  #'top_cluster':  fixed cluster number
  #'top_gene': find the number of cluster whose gene number reached the threshold
  #param threshold: cluster number or gene number
  #cluster_df_all: 14430 genes clustering
  
  if (select_gene_method == 'top_cluster') {
    cluster_num<-cluster_rank[["cluster"]][1:threshold] #from no.1 to no.threshold cluster
  }
  else if (select_gene_method == 'top_gene') {
    #top N gene
    for (topN_cluster in seq_len(length(cluster_rank$Count_genes))) {
      if (sum(cluster_rank$Count_genes[1:topN_cluster])>=threshold) {
        #The number of genes has reached the threshold
        break
      }
    }
    cluster_num<-cluster_rank[["cluster"]][1:topN_cluster]
    print(paste('topN_cluster:',topN_cluster))
  }
  else {cluster_num<-cluster_rank[["cluster"]][1:20]}
  
  gene_exper<-filter_id_df[which(filter_id_df$cluster%in%cluster_num),][,c("ENTREZID","original_gene")]
  names(gene_exper) <- c("ENTREZID_exper", "SYMBOL_exper")
  
  gene_count<-nrow(gene_exper)
  print(paste('gene_count:',gene_count))
  
  #random
  set.seed(2022)
  gene_ctrl <- cluster_df_all[sample(nrow(cluster_df_all),gene_count),"original_gene"]
  gene_ctrl <- as.vector(gene_ctrl)
  #convert to ID
  eg <- bitr(gene_ctrl,
             fromType="SYMBOL",
             toType=c("ENTREZID"),
             OrgDb="org.Hs.eg.db"
  )
  ENTREZID <- rep("NA",time = length(gene_ctrl))
  ENTREZID[match(eg$SYMBOL,gene_ctrl)] <- eg$ENTREZID
  gene_ctrl_id_df <- cbind(ENTREZID_ctrl = ENTREZID, SYMBOL_ctrl = gene_ctrl)
  gene_ctrl_id_df <- data.frame(gene_ctrl_id_df)
  gene_selected_df<-cbind(gene_exper, gene_ctrl_id_df)
  rownames(gene_selected_df)<-c(1:gene_count)
  
  gene_selected_df
}

get_NV_auc <- function (gene_log_C_orin, gene_selected_df, group, ano_min) {
  #param gene_log_C_orin: gene log-count expression matrix (14430 genes)
  #param gene_selected_df: gene symbols to be analyzed, including experimental and control groups
  #param group: a string, either 'exper' (i.e. experimental group) or 'ctrl' (i.e. control group) 
  #param ano_min: function {filter_network_cols} filter GO annotations
  
  group_selected<-paste('SYMBOL', group, sep='_')
  genelist<-as.vector(gene_selected_df[,group_selected])
  
  expression_m <- gene_log_C_orin[rownames(gene_log_C_orin)%in%genelist,]
  
  sm<-stats::cor(t(expression_m), method = 'spearman')
  sm<-abs(sm)
  sm<-as.matrix(sm)
  
  #get annotation data
  gene_ano<-AnnotationDbi::select(org.Hs.eg.db
                                  , keys=genelist
                                  , columns='GOALL'
                                  , keytype="SYMBOL")
  gene_ano_BP<-gene_ano[which(gene_ano$ONTOLOGYALL=='BP'),]
  goterms<-unique(gene_ano_BP$GOALL)
  annotations <- make_annotations(gene_ano_BP[,c('SYMBOL','GOALL')], genelist, goterms)
  annotations_sub <- filter_network_cols(annotations, min=ano_min, max=length(genelist)+1)  
  
  # run neighbor voting algorithm
  roc.sub <- neighbor_voting(annotations_sub, sm, nFold=3)
  roc.sub
}


##get KEGG_annotations from KEGGREST
symbolOnly <- function(x){
  items <- strsplit(x, ";", fixed = TRUE) %>% unlist()
  return(items[1])
}

geneSymbol <- function(x){
  #keggGet(x)[[1]]$GENE： a vector
  #where the odd rows are ENTREZID and the even rows are SYMBOL
  geneList <- keggGet(x)[[1]]$GENE
  if(!is.null(geneList)){
    listLength <- length(geneList)
    symbolList <- geneList[seq.int(from = 2, by = 2, length.out = listLength/2)] %>% map_chr(symbolOnly)
    symbol <- stringr::str_c(symbolList, collapse = ",")
    return(symbol)
  }else{
    return("")
  }
}

pathwayID <- function(x){
  #Get the pathway IDs starting with 'hsa'
  items <- strsplit(x, ":", fixed = TRUE) %>% unlist()
  return(items[2])
}

KEGG_annotation_2col <- function(){
  #get KEGG_annotations of 2cols
  hsaList <- keggList("pathway", "hsa")
  IDList <- names(hsaList) %>% map_chr(pathwayID)
  hsaPathway <- tibble::tibble(pathway_id=IDList, pathway_name=hsaList)
  pathwayFull <- hsaPathway %>% dplyr::mutate(hgnc_symbol=map_chr(pathway_id, geneSymbol))
  #remove null pathways
  pathwayWithGene <- dplyr::filter(pathwayFull, hgnc_symbol != "")
  #rebuild data
  KEGG_ano <- data.frame()
  for (i in 1:nrow(pathwayWithGene)) {
    genelist <- strsplit(pathwayWithGene$hgnc_symbol[i], ',')[[1]]
    pathway_df <- data.frame(pathway_id=pathwayWithGene$pathway_id[i], symbol=genelist)
    KEGG_ano <- rbind(KEGG_ano,pathway_df)
  }
  KEGG_ano
}
