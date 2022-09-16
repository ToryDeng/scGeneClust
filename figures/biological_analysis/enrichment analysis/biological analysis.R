setwd('./scGeneClust/figures/biological_analysis/enrichment analysis')

##################################################################################
##################################################################################
##GeneClust-fast
#PBMCSLEctrl
load("PBMCSLEctrl.Rdata")
cluster_df<-rowData(sce)
cluster_df<-as.data.frame(cluster_df)
cluster_score_f<-filter_single_g(cluster_df)
cluster_score_f_eg<-convert_symbol2id(cluster_score_f)
cluster_score_cluster<-calculate_cluster_score(cluster_score_f_eg,'seurat')
cluster_rank<-rank_cluster(cluster_score_cluster)
gene_selected_df<-select_genes(cluster_score_cluster,cluster_rank,'top_gene',100,cluster_df)
write.csv(gene_selected_df, 'PBMCSLEctrl_top_gene_100_selected_gene.csv',row.names = F)

#Cluster_id of the experimental group
cluster_id <- cluster_score_cluster[cluster_score_cluster$original_gene%in%gene_selected_df$SYMBOL_exper,]
write.csv(cluster_id, 'PBMCSLEctrl_top_gene_100_exper_cluster_id.csv',row.names = F)


##experimental group
#GO enrichment
GO_exper <- enrichGO(gene_selected_df$SYMBOL_exper,
                     OrgDb = org.Hs.eg.db,
                     ont='BP',
                     keyType = "SYMBOL",
                     pAdjustMethod = 'BH',
                     pvalueCutoff = 0.5,
                     qvalueCutoff = 0.5
)

#plot enrichment map
GO_exper2_similar <- pairwise_termsim(GO_exper)
#calculate logPadj for each GO term
GO_exper2_similar_log <- mutate(GO_exper2_similar, logPadj = -log(p.adjust,base = 10)) 

#replace the Description with the corresponding ID
result <- GO_exper2_similar_log@result[1:30,] #showCategory=30
termsim <- GO_exper2_similar_log@termsim[1:30,1:30] #showCategory=30
GO_exper2_similar_cp <- GO_exper2_similar_log
GO_exper2_similar_cp@result <- result
GO_exper2_similar_cp@result[,2] <- result[,1]
GO_exper2_similar_cp@termsim <- termsim
colnames(GO_exper2_similar_cp@termsim) <- result[,1]
rownames(GO_exper2_similar_cp@termsim) <- result[,1]  

range(GO_exper2_similar_cp@result[1:20,'logPadj']) #5.873733 7.830959
range(GO_exper2_similar_cp@result[1:20,'Count'])  #5 18
p1 <- emapplot(GO_exper2_similar_cp,showCategory = 20,layout = 'kk',repel = T,
               cex_label_category  = 0.9,
               cex_line = 0.5,
               color='logPadj',min_edge = 0.4)
p11 <- p1 + scale_size(range = c(5,12),limits = c(2,20),breaks = c(5,10,15,20),
                       name = "Number of genes")+
  scale_fill_gradient(high = "#FF0000",low  = "#000080",name = "-log p.adj",
                      limit = c(0.5,8.5), breaks=c(2,4,6,8))
ggsave(filename = 'PBMCSLEctrl_top_gene_100_GO_exper_emap20.png',p11,width = 16, height = 9, units = "in", dpi = 300)

#cofunction analyses
gene_log_C <- assay(sce, "X_gene_log")
#convert the rownames of 'gene_log_C' to origin gene symbols
gene_log_C_orin <- gene_log_C
m <- match(rownames(gene_log_C_orin),rownames(cluster_df))
orin_log_genes <- cluster_df[m,'original_gene']
orin_log_genes <- as.vector(orin_log_genes)
rownames(gene_log_C_orin) <- orin_log_genes

#run NV
roc.sub<-get_NV_auc(gene_log_C_orin, gene_selected_df, 'exper', ano_min=1)
m<-match(GO_exper@result$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
GO_exper_ORA_NV<-cbind(GO_exper@result,NV_df)

#relationship to immune process
GO_ANC <- as.list(GOBPANCESTOR)
anc<-GO_ANC[GO_exper@result$ID]
anc_res<-data.frame()
for (goterm in names(anc)) {
  for (i in anc[[goterm]]) {
    if (i=='GO:0002376') {
      anc_res[goterm,'is_immune']<-1
      break
    }
    else {anc_res[goterm,'is_immune']<-0}
  }
}
GO_exper_ORA_NV_immune<-cbind(GO_exper_ORA_NV,anc_res)
write.csv(GO_exper_ORA_NV_immune, 'PBMCSLEctrl_top_gene_100_GO_exper_summary.csv')


##control group
#GO enrichment
GO_ctrl <- enrichGO(gene_selected_df$SYMBOL_ctrl,
                     OrgDb = org.Hs.eg.db,
                     ont='BP',
                     keyType = "SYMBOL",
                     pAdjustMethod = 'BH',
                     pvalueCutoff = 1,
                     qvalueCutoff = 1
)

#plot enrichment map
GO_ctrl2_similar <- pairwise_termsim(GO_ctrl)
#calculate logPadj for each GO term
GO_ctrl2_similar_log <- mutate(GO_ctrl2_similar, logPadj = -log(p.adjust,base = 10)) 

#replace the Description with the corresponding ID
result <- GO_ctrl2_similar_log@result[1:30,] #showCategory=30
termsim <- GO_ctrl2_similar_log@termsim[1:30,1:30] #showCategory=30
GO_ctrl2_similar_cp <- GO_ctrl2_similar_log
GO_ctrl2_similar_cp@result <- result
GO_ctrl2_similar_cp@result[,2] <- result[,1]
GO_ctrl2_similar_cp@termsim <- termsim
colnames(GO_ctrl2_similar_cp@termsim) <- result[,1]
rownames(GO_ctrl2_similar_cp@termsim) <- result[,1] 

range(GO_ctrl2_similar_cp@result[1:20,'logPadj']) #0.6505743 0.7293364
range(GO_ctrl2_similar_cp@result[1:20,'Count'])  #2 5
p2 <- emapplot(GO_ctrl2_similar_cp,showCategory = 20,layout = 'fr',repel = T,
               cex_label_category  = 0.9,
               cex_line = 0.5,
               color='logPadj',min_edge = 0.4)
p21 <- p2 + scale_size(range = c(5,12),limits = c(2,20),breaks = c(5,10,15,20),
                       name = "Number of genes")+
  scale_fill_gradient(high = "#FF0000",low  = "#000080",name = "-log p.adj",
                      limit = c(0.5,8.5), breaks=c(2,4,6,8))
ggsave(filename = 'PBMCSLEctrl_top_gene_100_GO_ctrl_emap20.png',p21,width = 16, height = 9, units = "in", dpi = 300)

#merge enrichment maps of the experimental group and the control group into a graph
p3 <- p11+p21+plot_layout(guides = "collect",widths = c(1, 1))
ggsave(filename = 'PBMCSLEctrl_top_gene_100_GO_vs_emap20.png',p3,
       width = 15, height = 9, units = "in", dpi = 300)	

#cofunction analyses
#run NV
roc.sub<-get_NV_auc(gene_log_C_orin, gene_selected_df, 'ctrl', ano_min=1)
m<-match(GO_ctrl@result$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
GO_ctrl_ORA_NV<-cbind(GO_ctrl@result,NV_df)

#relationship to immune process
GO_ANC <- as.list(GOBPANCESTOR)
anc<-GO_ANC[GO_ctrl@result$ID]
anc_res<-data.frame()
for (goterm in names(anc)) {
  for (i in anc[[goterm]]) {
    if (i=='GO:0002376') {
      anc_res[goterm,'is_immune']<-1
      break
    }
    else {anc_res[goterm,'is_immune']<-0}
  }
}
GO_ctrl_ORA_NV_immune<-cbind(GO_ctrl_ORA_NV,anc_res)
write.csv(GO_ctrl_ORA_NV_immune, 'PBMCSLEctrl_top_gene_100_GO_ctrl_summary.csv')


##KEGG enrichment
#experimental group
KEGG_exper <- enrichKEGG(gene_selected_df$ENTREZID_exper,
                         organism = 'hsa',
                         pvalueCutoff = 0.5)

#plot enrichment map
KEGG_exper2_similar <- pairwise_termsim(KEGG_exper)
#calculate logPadj for each KEGG pathway
KEGG_exper2_similar_log <- mutate(KEGG_exper2_similar, logPadj = -log(p.adjust,base = 10)) 

#replace the Description with the corresponding ID
result <- KEGG_exper2_similar_log@result[1:30,] #showCategory=30
#mark SLE pathway
result[result$ID=="hsa05322",1] <- " "
termsim <- KEGG_exper2_similar_log@termsim[1:30,1:30] #showCategory=30
KEGG_exper2_similar_cp <- KEGG_exper2_similar_log
KEGG_exper2_similar_cp@result <- result
KEGG_exper2_similar_cp@result[,2] <- result[,1]
KEGG_exper2_similar_cp@termsim <- termsim
colnames(KEGG_exper2_similar_cp@termsim) <- result[,1]
rownames(KEGG_exper2_similar_cp@termsim) <- result[,1] 

range(KEGG_exper2_similar_cp@result[1:20,'logPadj']) #2.591219 6.519825
range(KEGG_exper2_similar_cp@result[1:20,'Count'])  #5 13
p1 <- emapplot(KEGG_exper2_similar_cp,showCategory = 20,layout = 'kk',repel = T,
               cex_label_category  = 0.9,
               cex_line = 0.5,
               color='logPadj',min_edge = 0.4) 
#which(KEGG_exper2_similar@result[1:30,1]=="hsa05322") #19
p11 <- p1 + scale_size(range = c(5,12),limits = c(1,13),breaks = c(3,8,13),
                       name = "Number of genes")+
  scale_fill_gradient(high = "#FF0000",low  = "#000080",name = "-log p.adj",
                      limit = c(0.1,7.2), breaks=c(1,3,5,7))+
  geom_node_text(aes(x =x, y = y,label = "hsa05322(SLE)"),
                 size = 4.75,
                 color = c(rep("transparent",18),"red","transparent"),
                 show.legend = F, repel = T,
                 bg.color = c(rep("transparent",18),"white","transparent"))
ggsave(filename = 'PBMCSLEctrl_top_gene_100_KEGG_exper_emap20.png',
       p11,width = 16, height = 9, units = "in", dpi = 300)


#convert ENTREZID to SYMBOL of the list "geneID" of KEGG enrichment result 
KEGG_exper_res <- KEGG_exper@result
KEGG_geneID <- KEGG_exper_res$geneID
for (i in seq_len(length(KEGG_geneID))) {
  l <- strsplit(KEGG_geneID[i], '/')[[1]]
  m <- match(l, gene_selected_df$ENTREZID_exper)  #exper
  symbol_list <- gene_selected_df$SYMBOL_exper[m] #exper
  gene_symbol <- paste(symbol_list, collapse = '/')
  KEGG_geneID[i] <- gene_symbol
}
KEGG_exper_res$geneID<- KEGG_geneID

#cofunction analyses
expression_m <- gene_log_C_orin[rownames(gene_log_C_orin)%in%gene_selected_df$SYMBOL_exper,]
sm<-stats::cor(t(expression_m), method = 'spearman')
sm<-abs(sm)
sm<-as.matrix(sm)

#run NV
KEGG_ano <- KEGG_annotation_2col()
keggpathways<-unique(KEGG_ano$pathway_id)
annotations <- make_annotations(KEGG_ano[,c('symbol','pathway_id')], gene_selected_df$SYMBOL_exper, keggpathways)
annotations_sub <- filter_network_cols(annotations, min=0, max=length(gene_selected_df$SYMBOL_exper)+1)	
roc.sub <- neighbor_voting(annotations_sub, sm, nFold=3)

m<-match(KEGG_exper_res$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
KEGG_exper_res_ORA_NV<-cbind(KEGG_exper_res,NV_df)
write.csv(KEGG_exper_res_ORA_NV, 'PBMCSLEctrl_top_gene_100_KEGG_exper_summary.csv')

##control group
#KEGG enrichment
KEGG_ctrl <- enrichKEGG(gene_selected_df$ENTREZID_ctrl,
                         organism = 'hsa',
                         pvalueCutoff = 1,
                         qvalueCutoff = 1)

#plot enrichment map
KEGG_ctrl2_similar <- pairwise_termsim(KEGG_ctrl)
#calculate logPadj for each KEGG pathway
KEGG_ctrl2_similar_log <- mutate(KEGG_ctrl2_similar, logPadj = -log(p.adjust,base = 10)) 

#replace the Description with the corresponding ID
result <- KEGG_ctrl2_similar_log@result[1:30,] #showCategory=30
termsim <- KEGG_ctrl2_similar_log@termsim[1:30,1:30] #showCategory=30
KEGG_ctrl2_similar_cp <- KEGG_ctrl2_similar_log
KEGG_ctrl2_similar_cp@result <- result
KEGG_ctrl2_similar_cp@result[,2] <- result[,1]
KEGG_ctrl2_similar_cp@termsim <- termsim
colnames(KEGG_ctrl2_similar_cp@termsim) <- result[,1]
rownames(KEGG_ctrl2_similar_cp@termsim) <- result[,1]

range(KEGG_ctrl2_similar_cp@result[1:20,'logPadj']) #0.2370723 0.2370723
range(KEGG_ctrl2_similar_cp@result[1:20,'Count'])  #1 5
p2 <- emapplot(KEGG_ctrl2_similar_cp,showCategory = 20,layout = 'kk',repel = T,
               cex_label_category  = 0.9,
               cex_line = 0.5,
               color='logPadj',min_edge = 0.4) 
p21 <- p2 + scale_size(range = c(5,12),limits = c(1,13),breaks = c(3,8,13),
                       name = "Number of genes")+
  scale_fill_gradient(high = "#FF0000",low  = "#000080",name = "-log p.adj",
                      limit = c(0.1,7.2), breaks=c(1,3,5,7))
ggsave(filename = 'PBMCSLEctrl_top_gene_100_KEGG_ctrl_emap20.png',
       p21,width = 16, height = 9, units = "in", dpi = 300)

#merge enrichment maps of the experimental group and the control group into a graph
p3 <- p11+p21+plot_layout(guides = "collect",widths = c(1, 1))
ggsave(filename = 'PBMCSLEctrl_top_gene_100_KEGG_vs_emap20.png',p3,
       width = 15, height = 9, units = "in", dpi = 300) 


#convert ENTREZID to SYMBOL of the list "geneID" of KEGG enrichment result 
KEGG_ctrl_res <- KEGG_ctrl@result
KEGG_geneID <- KEGG_ctrl_res$geneID
for (i in seq_len(length(KEGG_geneID))) {
  l <- strsplit(KEGG_geneID[i], '/')[[1]]
  m <- match(l, gene_selected_df$ENTREZID_ctrl)  #ctrl
  symbol_list <- gene_selected_df$SYMBOL_ctrl[m] #ctrl
  gene_symbol <- paste(symbol_list, collapse = '/')
  KEGG_geneID[i] <- gene_symbol
}
KEGG_ctrl_res$geneID<- KEGG_geneID

#cofunction analyses
expression_m <- gene_log_C_orin[rownames(gene_log_C_orin)%in%gene_selected_df$SYMBOL_ctrl,]
sm<-stats::cor(t(expression_m), method = 'spearman')
sm<-abs(sm)
sm<-as.matrix(sm)

#run NV
KEGG_ano <- KEGG_annotation_2col()
keggpathways<-unique(KEGG_ano$pathway_id)
annotations <- make_annotations(KEGG_ano[,c('symbol','pathway_id')], gene_selected_df$SYMBOL_ctrl, keggpathways)
annotations_sub <- filter_network_cols(annotations, min=0, max=length(gene_selected_df$SYMBOL_ctrl)+1)	
roc.sub <- neighbor_voting(annotations_sub, sm, nFold=3)

m<-match(KEGG_ctrl_res$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
KEGG_ctrl_res_ORA_NV<-cbind(KEGG_ctrl_res,NV_df)
write.csv(KEGG_ctrl_res_ORA_NV, 'PBMCSLEctrl_top_gene_100_KEGG_ctrl_summary.csv')





##################################################################################
##################################################################################
##GeneClust-fast
#PBMCSLEstim
load("PBMCSLEstim.Rdata")
cluster_df<-rowData(sce)
cluster_df<-as.data.frame(cluster_df)
cluster_score_f<-filter_single_g(cluster_df)
cluster_score_f_eg<-convert_symbol2id(cluster_score_f)
cluster_score_cluster<-calculate_cluster_score(cluster_score_f_eg,'seurat')
cluster_rank<-rank_cluster(cluster_score_cluster)
gene_selected_df<-select_genes(cluster_score_cluster,cluster_rank,'top_gene',100,cluster_df)
write.csv(gene_selected_df, 'PBMCSLEstim_top_gene_100_selected_gene.csv',row.names = F)

#Cluster_id of the experimental group
cluster_id <- cluster_score_cluster[cluster_score_cluster$original_gene%in%gene_selected_df$SYMBOL_exper,]
write.csv(cluster_id, 'PBMCSLEstim_top_gene_100_exper_cluster_id.csv',row.names = F)


##experimental group
#GO enrichment
GO_exper <- enrichGO(gene_selected_df$SYMBOL_exper,
                     OrgDb = org.Hs.eg.db,
                     ont='BP',
                     keyType = "SYMBOL",
                     pAdjustMethod = 'BH',
                     pvalueCutoff = 0.5,
                     qvalueCutoff = 0.5
)

#cofunction analyses
gene_log_C <- assay(sce, "X_gene_log")
#convert the rownames of 'gene_log_C' to origin gene symbols
gene_log_C_orin <- gene_log_C
m <- match(rownames(gene_log_C_orin),rownames(cluster_df))
orin_log_genes <- cluster_df[m,'original_gene']
orin_log_genes <- as.vector(orin_log_genes)
rownames(gene_log_C_orin) <- orin_log_genes

#run NV
roc.sub<-get_NV_auc(gene_log_C_orin, gene_selected_df, 'exper', ano_min=1)
m<-match(GO_exper@result$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
GO_exper_ORA_NV<-cbind(GO_exper@result,NV_df)

#relationship to immune process
GO_ANC <- as.list(GOBPANCESTOR)
anc<-GO_ANC[GO_exper@result$ID]
anc_res<-data.frame()
for (goterm in names(anc)) {
  for (i in anc[[goterm]]) {
    if (i=='GO:0002376') {
      anc_res[goterm,'is_immune']<-1
      break
    }
    else {anc_res[goterm,'is_immune']<-0}
  }
}
GO_exper_ORA_NV_immune<-cbind(GO_exper_ORA_NV,anc_res)
write.csv(GO_exper_ORA_NV_immune, 'PBMCSLEstim_top_gene_100_GO_exper_summary.csv')


##control group
#GO enrichment
GO_ctrl <- enrichGO(gene_selected_df$SYMBOL_ctrl,
                     OrgDb = org.Hs.eg.db,
                     ont='BP',
                     keyType = "SYMBOL",
                     pAdjustMethod = 'BH',
                     pvalueCutoff = 1,
                     qvalueCutoff = 1
)

#cofunction analyses
#run NV
roc.sub<-get_NV_auc(gene_log_C_orin, gene_selected_df, 'ctrl', ano_min=1)
m<-match(GO_ctrl@result$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
GO_ctrl_ORA_NV<-cbind(GO_ctrl@result,NV_df)

#relationship to immune process
GO_ANC <- as.list(GOBPANCESTOR)
anc<-GO_ANC[GO_ctrl@result$ID]
anc_res<-data.frame()
for (goterm in names(anc)) {
  for (i in anc[[goterm]]) {
    if (i=='GO:0002376') {
      anc_res[goterm,'is_immune']<-1
      break
    }
    else {anc_res[goterm,'is_immune']<-0}
  }
}
GO_ctrl_ORA_NV_immune<-cbind(GO_ctrl_ORA_NV,anc_res)
write.csv(GO_ctrl_ORA_NV_immune, 'PBMCSLEstim_top_gene_100_GO_ctrl_summary.csv')


##KEGG enrichment
#experimental group
KEGG_exper <- enrichKEGG(gene_selected_df$ENTREZID_exper,
                         organism = 'hsa',
                         pvalueCutoff = 0.5)

#convert ENTREZID to SYMBOL of the list "geneID" of KEGG enrichment result 
KEGG_exper_res <- KEGG_exper@result
KEGG_geneID <- KEGG_exper_res$geneID
for (i in seq_len(length(KEGG_geneID))) {
  l <- strsplit(KEGG_geneID[i], '/')[[1]]
  m <- match(l, gene_selected_df$ENTREZID_exper)  #exper
  symbol_list <- gene_selected_df$SYMBOL_exper[m] #exper
  gene_symbol <- paste(symbol_list, collapse = '/')
  KEGG_geneID[i] <- gene_symbol
}
KEGG_exper_res$geneID<- KEGG_geneID

#cofunction analyses
expression_m <- gene_log_C_orin[rownames(gene_log_C_orin)%in%gene_selected_df$SYMBOL_exper,]
sm<-stats::cor(t(expression_m), method = 'spearman')
sm<-abs(sm)
sm<-as.matrix(sm)

#run NV
KEGG_ano <- KEGG_annotation_2col()
keggpathways<-unique(KEGG_ano$pathway_id)
annotations <- make_annotations(KEGG_ano[,c('symbol','pathway_id')], gene_selected_df$SYMBOL_exper, keggpathways)
annotations_sub <- filter_network_cols(annotations, min=0, max=length(gene_selected_df$SYMBOL_exper)+1)	
roc.sub <- neighbor_voting(annotations_sub, sm, nFold=3)

m<-match(KEGG_exper_res$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
KEGG_exper_res_ORA_NV<-cbind(KEGG_exper_res,NV_df)
write.csv(KEGG_exper_res_ORA_NV, 'PBMCSLEstim_top_gene_100_KEGG_exper_summary.csv')

##control group
#KEGG enrichment
KEGG_ctrl <- enrichKEGG(gene_selected_df$ENTREZID_ctrl,
                         organism = 'hsa',
                         pvalueCutoff = 1,
                         qvalueCutoff = 1)

#convert ENTREZID to SYMBOL of the list "geneID" of KEGG enrichment result 
KEGG_ctrl_res <- KEGG_ctrl@result
KEGG_geneID <- KEGG_ctrl_res$geneID
for (i in seq_len(length(KEGG_geneID))) {
  l <- strsplit(KEGG_geneID[i], '/')[[1]]
  m <- match(l, gene_selected_df$ENTREZID_ctrl)  #ctrl
  symbol_list <- gene_selected_df$SYMBOL_ctrl[m] #ctrl
  gene_symbol <- paste(symbol_list, collapse = '/')
  KEGG_geneID[i] <- gene_symbol
}
KEGG_ctrl_res$geneID<- KEGG_geneID

#cofunction analyses
expression_m <- gene_log_C_orin[rownames(gene_log_C_orin)%in%gene_selected_df$SYMBOL_ctrl,]
sm<-stats::cor(t(expression_m), method = 'spearman')
sm<-abs(sm)
sm<-as.matrix(sm)

#run NV
KEGG_ano <- KEGG_annotation_2col()
keggpathways<-unique(KEGG_ano$pathway_id)
annotations <- make_annotations(KEGG_ano[,c('symbol','pathway_id')], gene_selected_df$SYMBOL_ctrl, keggpathways)
annotations_sub <- filter_network_cols(annotations, min=0, max=length(gene_selected_df$SYMBOL_ctrl)+1)	
roc.sub <- neighbor_voting(annotations_sub, sm, nFold=3)

m<-match(KEGG_ctrl_res$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
KEGG_ctrl_res_ORA_NV<-cbind(KEGG_ctrl_res,NV_df)
write.csv(KEGG_ctrl_res_ORA_NV, 'PBMCSLEstim_top_gene_100_KEGG_ctrl_summary.csv')




##################################################################################
##################################################################################
##GeneClust-ps
#PBMCSLEctrl
load("PBMCSLEctrl-ps.Rdata")
cluster_df<-rowData(sce)
cluster_df<-as.data.frame(cluster_df)
cluster_score_f<-filter_single_g(cluster_df)
cluster_score_f_eg<-convert_symbol2id(cluster_score_f)

#func{calculate_cluster_score} recalculate:
clist<-unique(cluster_score_f_eg$cluster)
#remove 'variances_norm'
cluster_score_f_eg <- select(cluster_score_f_eg, -('variances_norm'))

#load log-express matrix
rm(sce)
load("PBMCSLEctrl.Rdata")
#for control group
cluster_df_all<-rowData(sce) 
cluster_df_all <- as.data.frame(cluster_df_all)

m <- match(cluster_score_f_eg$original_gene, cluster_df_all$original_gene)
variances_norm <- cluster_df_all[m,c('variances_norm')]
rebuild_cluster <- cbind(cluster_score_f_eg, variances_norm)

cluster_score_cal <- calculate_cluster_score(rebuild_cluster,'relevant')
cluster_rank<-rank_cluster(cluster_score_cal)

#top 100 gene
threshold <- 100
for (topN_cluster in seq_len(length(cluster_rank$Count_genes))) {
  if (sum(cluster_rank$Count_genes[1:topN_cluster])>=threshold) {
    #The number of genes has reached the threshold
    break
  }
}

cluster_num<-cluster_rank[["cluster"]][1:topN_cluster]
gene_exper<-rebuild_cluster[which(rebuild_cluster$cluster%in%cluster_num),][,c("ENTREZID","original_gene")]
names(gene_exper) <- c("ENTREZID_exper", "SYMBOL_exper")

#Cluster_id of the experimental group
cluster_id <- cluster_score_cal[cluster_score_cal$original_gene%in%gene_exper$SYMBOL_exper,]
write.csv(cluster_id, 'PBMCSLEctrl-ps_top_gene_100_exper_cluster_id.csv',row.names = F)

#Control group: randomly selected from homo sapiens genes with the same number:
gene_count<-nrow(gene_exper)
set.seed(2022)
gene_ctrl <- cluster_df_all[sample(nrow(cluster_df_all),gene_count),"original_gene"]
gene_ctrl <- as.vector(gene_ctrl)
#convert SYMBOL to ENTREZID
eg <- bitr(gene_ctrl,
           fromType="SYMBOL",
           toType=c("ENTREZID"),
           OrgDb="org.Hs.eg.db"
)
ENTREZID <- rep("NA",time = length(gene_ctrl))
ENTREZID[match(eg$SYMBOL,gene_ctrl)] <- eg$ENTREZID
gene_ctrl_id_df <- cbind(ENTREZID_ctrl = ENTREZID, SYMBOL_ctrl = gene_ctrl)
gene_ctrl_id_df <- data.frame(gene_ctrl_id_df)

#Genes finally selected
gene_selected_df<-cbind(gene_exper, gene_ctrl_id_df)
rownames(gene_selected_df)<-c(1:gene_count)
write.csv(gene_selected_df, 'PBMCSLEctrl-ps_top_gene_100_selected_gene.csv',row.names = F)


##experimental group
#GO enrichment
GO_exper <- enrichGO(gene_selected_df$SYMBOL_exper,
                     OrgDb = org.Hs.eg.db,
                     ont='BP',
                     keyType = "SYMBOL",
                     pAdjustMethod = 'BH',
                     pvalueCutoff = 0.5,
                     qvalueCutoff = 1)

#plot enrichment map
GO_exper2_similar <- pairwise_termsim(GO_exper)
#calculate logPadj for each GO term
GO_exper2_similar_log <- mutate(GO_exper2_similar, logPadj = -log(p.adjust,base = 10)) 

#replace the Description with the corresponding ID
result <- GO_exper2_similar_log@result[1:30,] #showCategory=30
termsim <- GO_exper2_similar_log@termsim[1:30,1:30] #showCategory=30
GO_exper2_similar_cp <- GO_exper2_similar_log
GO_exper2_similar_cp@result <- result
GO_exper2_similar_cp@result[,2] <- result[,1]
GO_exper2_similar_cp@termsim <- termsim
colnames(GO_exper2_similar_cp@termsim) <- result[,1]
rownames(GO_exper2_similar_cp@termsim) <- result[,1]  

range(GO_exper2_similar_cp@result[1:20,'logPadj']) #2.083935 6.669167
range(GO_exper2_similar_cp@result[1:20,'Count'])  #3 10
p1 <- emapplot(GO_exper2_similar_log,showCategory = 20,layout = 'kk',repel = T,
               cex_label_category  = 0.9,
               cex_line = 0.5,
               color='logPadj',min_edge = 0.4)
p11 <- p1 + scale_size(range = c(5,12),limits = c(2,10),breaks = c(2,4,6,8,10),
                       name = "Number of genes")+
  scale_fill_gradient(high = "#FF0000",low  = "#000080",name = "-log p.adj",
                      limit = c(0.3,7.5), breaks=c(1,3,5,7))
ggsave(filename = 'PBMCSLEctrl-ps_top_gene_100_GO_exper_emap.png',p11,width = 16, height = 9, units = "in", dpi = 300)

#cofunction analyses
gene_log_C <- assay(sce, "X_gene_log")
#convert the rownames of 'gene_log_C' to origin gene symbols
gene_log_C_orin <- gene_log_C
m <- match(rownames(gene_log_C_orin),rownames(cluster_df_all))
orin_log_genes <- cluster_df_all[m,'original_gene']
orin_log_genes <- as.vector(orin_log_genes)
rownames(gene_log_C_orin) <- orin_log_genes

#run NV
roc.sub<-get_NV_auc(gene_log_C_orin, gene_selected_df, 'exper', ano_min=1)
m<-match(GO_exper@result$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
GO_exper_ORA_NV<-cbind(GO_exper@result,NV_df)
                      
#relationship to immune process
GO_ANC <- as.list(GOBPANCESTOR)
anc<-GO_ANC[GO_exper@result$ID]
anc_res<-data.frame()
for (goterm in names(anc)) {
  for (i in anc[[goterm]]) {
    if (i=='GO:0002376') {
      anc_res[goterm,'is_immune']<-1
      break
    }
    else {anc_res[goterm,'is_immune']<-0}
  }
}
GO_exper_ORA_NV_immune<-cbind(GO_exper_ORA_NV,anc_res)
write.csv(GO_exper_ORA_NV_immune, 'PBMCSLEctrl-ps_top_gene_100_GO_exper_summary.csv')


##control group
#GO enrichment
GO_ctrl <- enrichGO(gene_selected_df$SYMBOL_ctrl,
                    OrgDb = org.Hs.eg.db,
                    ont='BP',
                    keyType = "SYMBOL",
                    pAdjustMethod = 'BH',
                    pvalueCutoff = 0.5,
                    qvalueCutoff = 1)

#plot enrichment map
GO_ctrl2_similar <- pairwise_termsim(GO_ctrl)
#calculate logPadj for each GO term
GO_ctrl2_similar_log <- mutate(GO_ctrl2_similar, logPadj = -log(p.adjust,base = 10)) 

#replace the Description with the corresponding ID
result <- GO_ctrl2_similar_log@result[1:30,] #showCategory=30
termsim <- GO_ctrl2_similar_log@termsim[1:30,1:30] #showCategory=30
GO_ctrl2_similar_cp <- GO_ctrl2_similar_log
GO_ctrl2_similar_cp@result <- result
GO_ctrl2_similar_cp@result[,2] <- result[,1]
GO_ctrl2_similar_cp@termsim <- termsim
colnames(GO_ctrl2_similar_cp@termsim) <- result[,1]
rownames(GO_ctrl2_similar_cp@termsim) <- result[,1] 

range(GO_ctrl2_similar_cp@result[1:20,'logPadj']) #0.6239568 0.7633483
range(GO_ctrl2_similar_cp@result[1:20,'Count'])  #2 5
p2 <- emapplot(GO_ctrl2_similar_log,showCategory = 20,layout = 'fr',repel = T,
               cex_label_category  = 0.9,
               cex_line = 0.5,
               color='logPadj',min_edge = 0.4)
p21 <- p2 + scale_size(range = c(5,12),limits = c(2,10),breaks = c(2,4,6,8,10),
                       name = "Number of genes")+
  scale_fill_gradient(high = "#FF0000",low  = "#000080",name = "-log p.adj",
                      limit = c(0.3,7.5), breaks=c(1,3,5,7))
ggsave(filename = 'PBMCSLEctrl-ps_top_gene_100_GO_ctrl_emap.png',p21,width = 16, height = 9, units = "in", dpi = 300)

#merge enrichment maps of the experimental group and the control group into a graph
p3 <- p11+p21+plot_layout(guides = "collect",widths = c(1, 1))
ggsave(filename = 'PBMCSLEctrl-ps_top_gene_100_GO_vs_emap.png',p3,
       width = 15, height = 9, units = "in", dpi = 300) 

#cofunction analyses
#run NV
roc.sub<-get_NV_auc(gene_log_C_orin, gene_selected_df, 'ctrl', ano_min=1)
m<-match(GO_ctrl@result$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
GO_ctrl_ORA_NV<-cbind(GO_ctrl@result,NV_df)

#relationship to immune process
GO_ANC <- as.list(GOBPANCESTOR)
anc<-GO_ANC[GO_ctrl@result$ID]
anc_res<-data.frame()
for (goterm in names(anc)) {
  for (i in anc[[goterm]]) {
    if (i=='GO:0002376') {
      anc_res[goterm,'is_immune']<-1
      break
    }
    else {anc_res[goterm,'is_immune']<-0}
  }
}
GO_ctrl_ORA_NV_immune<-cbind(GO_ctrl_ORA_NV,anc_res)
write.csv(GO_ctrl_ORA_NV_immune, 'PBMCSLEctrl-ps_top_gene_100_GO_ctrl_summary.csv')


##potential gene for the experimental group
path_immune <- GO_exper_ORA_NV_immune[GO_exper_ORA_NV_immune$is_immune==1,'geneID']
path_immune <- as.vector(path_immune)
temp_immune <- vector()
for (i in seq_len(length(path_immune))) {
  l <- strsplit(path_immune[i], '/')[[1]]
  temp_immune <- append(temp_immune,l)
}
gene_immune <- unique(temp_immune)
gene_nonimmune <- as.vector(gene_selected_df$SYMBOL_exper[!gene_selected_df$SYMBOL_exper%in%gene_immune])

#calculate F score
celltype <- colData(sce)@listData[["celltype"]] 
gene_log_C_filter <- gene_log_C[rownames(gene_log_C)%in%rownames(cluster_df),]
F_res <- cal_F2(gene_log_C_filter,celltype)
#ixs_f = order(F_res$F_scores, decreasing = T) # order the features
F_gene <- F_res$F_scores
F_gene <- as.data.frame(F_gene)
m <- match(gene_nonimmune,rownames(gene_log_C_filter))
F_score_nonimmune <- as.data.frame(F_gene[m,])
colnames(F_score_nonimmune) <- 'F_score'
rownames(F_score_nonimmune) <- gene_nonimmune
F_score_nonimmune <- arrange(F_score_nonimmune,-F_score)
write.csv(F_score_nonimmune, 'PBMCSLEctrl-ps_potential_Fscore.csv')


##KEGG enrichment
#experimental group
KEGG_exper <- enrichKEGG(gene_selected_df$ENTREZID_exper,
                         organism = 'hsa',
                         pvalueCutoff = 0.5,
                         qvalueCutoff = 0.5)

#plot enrichment map
KEGG_exper2_similar <- pairwise_termsim(KEGG_exper)
#calculate logPadj for each KEGG pathway
KEGG_exper2_similar_log <- mutate(KEGG_exper2_similar, logPadj = -log(p.adjust,base = 10)) 

#replace the Description with the corresponding ID
result <- KEGG_exper2_similar_log@result[1:30,] #showCategory=30
#mark SLE pathway
result[result$ID=="hsa05322",1] <- " "
termsim <- KEGG_exper2_similar_log@termsim[1:30,1:30] #showCategory=30
KEGG_exper2_similar_cp <- KEGG_exper2_similar_log
KEGG_exper2_similar_cp@result <- result
KEGG_exper2_similar_cp@result[,2] <- result[,1]
KEGG_exper2_similar_cp@termsim <- termsim
colnames(KEGG_exper2_similar_cp@termsim) <- result[,1]
rownames(KEGG_exper2_similar_cp@termsim) <- result[,1] 

range(KEGG_exper2_similar_cp@result[1:20,'logPadj']) #1.910069 5.773722
range(KEGG_exper2_similar_cp@result[1:20,'Count'])  #4 11
p1 <- emapplot(KEGG_exper2_similar_cp,showCategory = 20,layout = 'kk',repel = T,
               cex_label_category  = 0.9,
               cex_line = 0.5,
               color='logPadj',min_edge = 0.4) 
#which(KEGG_exper2_similar@result[1:30,1]=="hsa05322")  #16
p11 <- p1 + scale_size(range = c(5,12),limits = c(1,12),breaks = c(3,6,9,12),
                       name = "Number of genes")+
  scale_fill_gradient(high = "#FF0000",low  = "#000080",name = "-log p.adj",
                      limit = c(0,6), breaks=c(1,3,5))+
  geom_node_text(aes(x =x, y = y,label = "hsa05322(SLE)"),
                 size = 4.75,
                 color = c(rep("transparent",15),"red",rep("transparent",4)),
                 show.legend = F, repel = T,
                 bg.color = c(rep("transparent",15),"white",rep("transparent",4)))
ggsave(filename = 'PBMCSLEctrl-ps_top_gene_100_KEGG_exper_emap.png',
       p11,width = 16, height = 9, units = "in", dpi = 300)

#convert ENTREZID to SYMBOL of the list "geneID" of KEGG enrichment result 
KEGG_exper_res <- KEGG_exper@result
KEGG_geneID <- KEGG_exper_res$geneID
for (i in seq_len(length(KEGG_geneID))) {
  l <- strsplit(KEGG_geneID[i], '/')[[1]]
  m <- match(l, gene_selected_df$ENTREZID_exper)  #exper
  symbol_list <- gene_selected_df$SYMBOL_exper[m] #exper
  gene_symbol <- paste(symbol_list, collapse = '/')
  KEGG_geneID[i] <- gene_symbol
}
KEGG_exper_res$geneID<- KEGG_geneID

#cofunction analyses
gene_log_C_orin_filter <- gene_log_C_orin[rownames(gene_log_C_orin)%in%gene_selected_df$SYMBOL_exper,]
sm<-stats::cor(t(gene_log_C_orin_filter), method = 'spearman')
sm<-abs(sm)
sm<-as.matrix(sm)

#run NV
KEGG_ano <- KEGG_annotation_2col()
keggpathways<-unique(KEGG_ano$pathway_id)
annotations <- make_annotations(KEGG_ano[,c('symbol','pathway_id')], gene_selected_df$SYMBOL_exper, keggpathways)
annotations_sub <- filter_network_cols(annotations, min=1, max=length(gene_selected_df$SYMBOL_exper)+1) 
roc.sub <- neighbor_voting(annotations_sub, sm, nFold=3)

m<-match(KEGG_exper_res$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
KEGG_exper_res_ORA_NV<-cbind(KEGG_exper_res,NV_df)
write.csv(KEGG_exper_res_ORA_NV, 'PBMCSLEctrl-ps_top_gene_100_KEGG_exper_summary.csv')

##control group
#KEGG enrichment
KEGG_ctrl <- enrichKEGG(gene_selected_df$ENTREZID_ctrl,
                         organism = 'hsa',
                         pvalueCutoff = 1,
                         qvalueCutoff = 1)

#plot enrichment map
KEGG_ctrl2_similar <- pairwise_termsim(KEGG_ctrl)
#calculate logPadj for each KEGG pathway
KEGG_ctrl2_similar_log <- mutate(KEGG_ctrl2_similar, logPadj = -log(p.adjust,base = 10)) 

#replace the Description with the corresponding ID
result <- KEGG_ctrl2_similar_log@result[1:30,] #showCategory=30
termsim <- KEGG_ctrl2_similar_log@termsim[1:30,1:30] #showCategory=30
KEGG_ctrl2_similar_cp <- KEGG_ctrl2_similar_log
KEGG_ctrl2_similar_cp@result <- result
KEGG_ctrl2_similar_cp@result[,2] <- result[,1]
KEGG_ctrl2_similar_cp@termsim <- termsim
colnames(KEGG_ctrl2_similar_cp@termsim) <- result[,1]
rownames(KEGG_ctrl2_similar_cp@termsim) <- result[,1]

range(KEGG_ctrl2_similar_cp@result[1:20,'logPadj']) #0.2370723 0.2370723
range(KEGG_ctrl2_similar_cp@result[1:20,'Count'])  #1 5
p2 <- emapplot(KEGG_ctrl2_similar_cp,showCategory = 20,layout = 'kk',repel = T,
               cex_label_category  = 0.9,
               cex_line = 0.5,
               color='logPadj',min_edge = 0.4) 
p21 <- p2 + scale_size(range = c(5,12),limits = c(1,12),breaks = c(3,6,9,12),
                       name = "Number of genes")+
  scale_fill_gradient(high = "#FF0000",low  = "#000080",name = "-log p.adj",
                      limit = c(0,6), breaks=c(1,3,5))
ggsave(filename = 'PBMCSLEctrl-ps_top_gene_100_KEGG_ctrl_emap.png',
       p21,width = 16, height = 9, units = "in", dpi = 300)

#merge enrichment maps of the experimental group and the control group into a graph
p3 <- p11+p21+plot_layout(guides = "collect",widths = c(1, 1))
ggsave(filename = 'PBMCSLEctrl-ps_top_gene_100_KEGG_vs_emap.png',p3,
       width = 15, height = 9, units = "in", dpi = 300) 


#convert ENTREZID to SYMBOL of the list "geneID" of KEGG enrichment result 
KEGG_ctrl_res <- KEGG_ctrl@result
KEGG_geneID <- KEGG_ctrl_res$geneID
for (i in seq_len(length(KEGG_geneID))) {
  l <- strsplit(KEGG_geneID[i], '/')[[1]]
  m <- match(l, gene_selected_df$ENTREZID_ctrl)  #ctrl
  symbol_list <- gene_selected_df$SYMBOL_ctrl[m] #ctrl
  gene_symbol <- paste(symbol_list, collapse = '/')
  KEGG_geneID[i] <- gene_symbol
}
KEGG_ctrl_res$geneID<- KEGG_geneID

#cofunction analyses
gene_log_C_orin_filter <- gene_log_C_orin[rownames(gene_log_C_orin)%in%gene_selected_df$SYMBOL_ctrl,]
sm<-stats::cor(t(gene_log_C_orin_filter), method = 'spearman')
sm<-abs(sm)
sm<-as.matrix(sm)

#run NV
KEGG_ano <- KEGG_annotation_2col()
keggpathways<-unique(KEGG_ano$pathway_id)
annotations <- make_annotations(KEGG_ano[,c('symbol','pathway_id')], gene_selected_df$SYMBOL_ctrl, keggpathways)
annotations_sub <- filter_network_cols(annotations, min=0, max=length(gene_selected_df$SYMBOL_ctrl)+1)  
# run neighbor voting algorithm
roc.sub <- neighbor_voting(annotations_sub, sm, nFold=3)

m<-match(KEGG_ctrl_res$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
KEGG_ctrl_res_ORA_NV<-cbind(KEGG_ctrl_res,NV_df)
write.csv(KEGG_ctrl_res_ORA_NV, 'PBMCSLEctrl-ps_top_gene_100_KEGG_ctrl_summary.csv')






##################################################################################
##################################################################################
##GeneClust-ps
#PBMCSLEstim
load("PBMCSLEstim-ps.Rdata")
cluster_df<-rowData(sce)
cluster_df<-as.data.frame(cluster_df)
cluster_score_f<-filter_single_g(cluster_df)
cluster_score_f_eg<-convert_symbol2id(cluster_score_f)

#func{calculate_cluster_score} recalculate:
clist<-unique(cluster_score_f_eg$cluster)
#remove 'variances_norm'
cluster_score_f_eg <- select(cluster_score_f_eg, -('variances_norm'))

#load log-express matrix
rm(sce)
load("PBMCSLEstim.Rdata")
#for control group
cluster_df_all<-rowData(sce) 
cluster_df_all <- as.data.frame(cluster_df_all)

m <- match(cluster_score_f_eg$original_gene, cluster_df_all$original_gene)
variances_norm <- cluster_df_all[m,c('variances_norm')]
rebuild_cluster <- cbind(cluster_score_f_eg, variances_norm)

cluster_score_cal <- calculate_cluster_score(rebuild_cluster,'relevant')
cluster_rank<-rank_cluster(cluster_score_cal)

#top 100 gene
threshold <- 100
for (topN_cluster in seq_len(length(cluster_rank$Count_genes))) {
  if (sum(cluster_rank$Count_genes[1:topN_cluster])>=threshold) {
    #The number of genes has reached the threshold
    break
  }
}

cluster_num<-cluster_rank[["cluster"]][1:topN_cluster]
gene_exper<-rebuild_cluster[which(rebuild_cluster$cluster%in%cluster_num),][,c("ENTREZID","original_gene")]
names(gene_exper) <- c("ENTREZID_exper", "SYMBOL_exper")

#Cluster_id of the experimental group
cluster_id <- cluster_score_cal[cluster_score_cal$original_gene%in%gene_exper$SYMBOL_exper,]
write.csv(cluster_id, 'PBMCSLEstim-ps_top_gene_100_exper_cluster_id.csv',row.names = F)

#Control group: randomly selected from homo sapiens genes with the same number:
gene_count<-nrow(gene_exper)
set.seed(2022)
gene_ctrl <- cluster_df_all[sample(nrow(cluster_df_all),gene_count),"original_gene"]
gene_ctrl <- as.vector(gene_ctrl)
#convert SYMBOL to ENTREZID
eg <- bitr(gene_ctrl,
           fromType="SYMBOL",
           toType=c("ENTREZID"),
           OrgDb="org.Hs.eg.db"
)
ENTREZID <- rep("NA",time = length(gene_ctrl))
ENTREZID[match(eg$SYMBOL,gene_ctrl)] <- eg$ENTREZID
gene_ctrl_id_df <- cbind(ENTREZID_ctrl = ENTREZID, SYMBOL_ctrl = gene_ctrl)
gene_ctrl_id_df <- data.frame(gene_ctrl_id_df)

#Genes finally selected
gene_selected_df<-cbind(gene_exper, gene_ctrl_id_df)
rownames(gene_selected_df)<-c(1:gene_count)
write.csv(gene_selected_df, 'PBMCSLEstim-ps_top_gene_100_selected_gene.csv',row.names = F)


##experimental group
#GO enrichment
GO_exper <- enrichGO(gene_selected_df$SYMBOL_exper,
                     OrgDb = org.Hs.eg.db,
                     ont='BP',
                     keyType = "SYMBOL",
                     pAdjustMethod = 'BH',
                     pvalueCutoff = 0.5,
                     qvalueCutoff = 1)

#cofunction analyses
gene_log_C <- assay(sce, "X_gene_log")
#convert the rownames of 'gene_log_C' to origin gene symbols
gene_log_C_orin <- gene_log_C
m <- match(rownames(gene_log_C_orin),rownames(cluster_df_all))
orin_log_genes <- cluster_df_all[m,'original_gene']
orin_log_genes <- as.vector(orin_log_genes)
rownames(gene_log_C_orin) <- orin_log_genes

#run NV
roc.sub<-get_NV_auc(gene_log_C_orin, gene_selected_df, 'exper', ano_min=1)
m<-match(GO_exper@result$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
GO_exper_ORA_NV<-cbind(GO_exper@result,NV_df)
                      
#relationship to immune process
GO_ANC <- as.list(GOBPANCESTOR)
anc<-GO_ANC[GO_exper@result$ID]
anc_res<-data.frame()
for (goterm in names(anc)) {
  for (i in anc[[goterm]]) {
    if (i=='GO:0002376') {
      anc_res[goterm,'is_immune']<-1
      break
    }
    else {anc_res[goterm,'is_immune']<-0}
  }
}
GO_exper_ORA_NV_immune<-cbind(GO_exper_ORA_NV,anc_res)
write.csv(GO_exper_ORA_NV_immune, 'PBMCSLEstim-ps_top_gene_100_GO_exper_summary.csv')


##control group
#GO enrichment
GO_ctrl <- enrichGO(gene_selected_df$SYMBOL_ctrl,
                    OrgDb = org.Hs.eg.db,
                    ont='BP',
                    keyType = "SYMBOL",
                    pAdjustMethod = 'BH',
                    pvalueCutoff = 0.5,
                    qvalueCutoff = 1)

#cofunction analyses
#run NV
roc.sub<-get_NV_auc(gene_log_C_orin, gene_selected_df, 'ctrl', ano_min=1)
m<-match(GO_ctrl@result$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
GO_ctrl_ORA_NV<-cbind(GO_ctrl@result,NV_df)

#relationship to immune process
GO_ANC <- as.list(GOBPANCESTOR)
anc<-GO_ANC[GO_ctrl@result$ID]
anc_res<-data.frame()
for (goterm in names(anc)) {
  for (i in anc[[goterm]]) {
    if (i=='GO:0002376') {
      anc_res[goterm,'is_immune']<-1
      break
    }
    else {anc_res[goterm,'is_immune']<-0}
  }
}
GO_ctrl_ORA_NV_immune<-cbind(GO_ctrl_ORA_NV,anc_res)
write.csv(GO_ctrl_ORA_NV_immune, 'PBMCSLEstim-ps_top_gene_100_GO_ctrl_summary.csv')


##potential gene for the experimental group
path_immune <- GO_exper_ORA_NV_immune[GO_exper_ORA_NV_immune$is_immune==1,'geneID']
path_immune <- as.vector(path_immune)
temp_immune <- vector()
for (i in seq_len(length(path_immune))) {
  l <- strsplit(path_immune[i], '/')[[1]]
  temp_immune <- append(temp_immune,l)
}
gene_immune <- unique(temp_immune)
gene_nonimmune <- as.vector(gene_selected_df$SYMBOL_exper[!gene_selected_df$SYMBOL_exper%in%gene_immune])

#calculate F score
celltype <- colData(sce)@listData[["celltype"]] 
gene_log_C_filter <- gene_log_C[rownames(gene_log_C)%in%rownames(cluster_df),]
F_res <- cal_F2(gene_log_C_filter,celltype)
#ixs_f = order(F_res$F_scores, decreasing = T) # order the features
F_gene <- F_res$F_scores
F_gene <- as.data.frame(F_gene)
m <- match(gene_nonimmune,rownames(gene_log_C_filter))
F_score_nonimmune <- as.data.frame(F_gene[m,])
colnames(F_score_nonimmune) <- 'F_score'
rownames(F_score_nonimmune) <- gene_nonimmune
F_score_nonimmune <- arrange(F_score_nonimmune,-F_score)
write.csv(F_score_nonimmune, 'PBMCSLEstim-ps_potential_Fscore.csv')


##KEGG enrichment
#experimental group
KEGG_exper <- enrichKEGG(gene_selected_df$ENTREZID_exper,
                         organism = 'hsa',
                         pvalueCutoff = 0.5,
                         qvalueCutoff = 0.5)

#convert ENTREZID to SYMBOL of the list "geneID" of KEGG enrichment result 
KEGG_exper_res <- KEGG_exper@result
KEGG_geneID <- KEGG_exper_res$geneID
for (i in seq_len(length(KEGG_geneID))) {
  l <- strsplit(KEGG_geneID[i], '/')[[1]]
  m <- match(l, gene_selected_df$ENTREZID_exper)  #exper
  symbol_list <- gene_selected_df$SYMBOL_exper[m] #exper
  gene_symbol <- paste(symbol_list, collapse = '/')
  KEGG_geneID[i] <- gene_symbol
}
KEGG_exper_res$geneID<- KEGG_geneID

#cofunction analyses
gene_log_C_orin_filter <- gene_log_C_orin[rownames(gene_log_C_orin)%in%gene_selected_df$SYMBOL_exper,]
sm<-stats::cor(t(gene_log_C_orin_filter), method = 'spearman')
sm<-abs(sm)
sm<-as.matrix(sm)

#run NV
KEGG_ano <- KEGG_annotation_2col()
keggpathways<-unique(KEGG_ano$pathway_id)
annotations <- make_annotations(KEGG_ano[,c('symbol','pathway_id')], gene_selected_df$SYMBOL_exper, keggpathways)
annotations_sub <- filter_network_cols(annotations, min=1, max=length(gene_selected_df$SYMBOL_exper)+1) 
roc.sub <- neighbor_voting(annotations_sub, sm, nFold=3)

m<-match(KEGG_exper_res$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
KEGG_exper_res_ORA_NV<-cbind(KEGG_exper_res,NV_df)
write.csv(KEGG_exper_res_ORA_NV, 'PBMCSLEstim-ps_top_gene_100_KEGG_exper_summary.csv')

##control group
#KEGG enrichment
KEGG_ctrl <- enrichKEGG(gene_selected_df$ENTREZID_ctrl,
                         organism = 'hsa',
                         pvalueCutoff = 1,
                         qvalueCutoff = 1)

#convert ENTREZID to SYMBOL of the list "geneID" of KEGG enrichment result 
KEGG_ctrl_res <- KEGG_ctrl@result
KEGG_geneID <- KEGG_ctrl_res$geneID
for (i in seq_len(length(KEGG_geneID))) {
  l <- strsplit(KEGG_geneID[i], '/')[[1]]
  m <- match(l, gene_selected_df$ENTREZID_ctrl)  #ctrl
  symbol_list <- gene_selected_df$SYMBOL_ctrl[m] #ctrl
  gene_symbol <- paste(symbol_list, collapse = '/')
  KEGG_geneID[i] <- gene_symbol
}
KEGG_ctrl_res$geneID<- KEGG_geneID

#cofunction analyses
gene_log_C_orin_filter <- gene_log_C_orin[rownames(gene_log_C_orin)%in%gene_selected_df$SYMBOL_ctrl,]
sm<-stats::cor(t(gene_log_C_orin_filter), method = 'spearman')
sm<-abs(sm)
sm<-as.matrix(sm)

#run NV
KEGG_ano <- KEGG_annotation_2col()
keggpathways<-unique(KEGG_ano$pathway_id)
annotations <- make_annotations(KEGG_ano[,c('symbol','pathway_id')], gene_selected_df$SYMBOL_ctrl, keggpathways)
annotations_sub <- filter_network_cols(annotations, min=0, max=length(gene_selected_df$SYMBOL_ctrl)+1)  
# run neighbor voting algorithm
roc.sub <- neighbor_voting(annotations_sub, sm, nFold=3)

m<-match(KEGG_ctrl_res$ID,rownames(roc.sub))
NV_df<-roc.sub[m,]
KEGG_ctrl_res_ORA_NV<-cbind(KEGG_ctrl_res,NV_df)
write.csv(KEGG_ctrl_res_ORA_NV, 'PBMCSLEstim-ps_top_gene_100_KEGG_ctrl_summary.csv')


