#packages######
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(robustbase)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(PMCMRplus)))
suppressWarnings(suppressMessages(library(ggnewscale)))
suppressWarnings(suppressMessages(library(patchwork)))

#paths####
path <- "./cell clustering.xlsx"

#function####
#' Calculating the median of 20-fold data.
#'
#' param datain: a data frame of a 20-fold data.
#' param nmethod: the number of methods in the 20-fold data.
#' return a data frame of the median of 20-fold data.
dfmedian <- function(datain,nmethod){
  cm <- as.data.frame(colMedians(as.matrix(
    datain[datain$dataset==unique(datain$dataset)[1],-c(1:3)])))
  md <- cbind(datain[datain$dataset==unique(datain$dataset)[1],1][1:nmethod], 
              as.data.frame(cm),as.data.frame(rownames(cm)))
  colnames(md) <- c("dataset","median","method")
  mdcom <- md
  for(i in 2:length(unique(datain$dataset))){
    cm <- as.data.frame(colMedians(as.matrix(
      datain[datain$dataset==unique(datain$dataset)[i],-c(1:3)])))
    md <- cbind(datain[datain$dataset==unique(datain$dataset)[i],1][1:nmethod], 
                as.data.frame(cm),as.data.frame(rownames(cm)))
    colnames(md) <- c("dataset","median","method")
    mdcom <- rbind(mdcom,md)
  }
  return(mdcom)
}

#' Modifying the name of data set.
#' 
#' param datain: a data frame of clustering results.
#' return a data frame of clustering results.
newname <- function(datain){
  for(i in 1:nrow(datain)){
    datain$dataset[i] <- switch(datain$dataset[i],
                                "Adam"="Adam",
                                "Chen"="Chen",
                                "Guo"="Guo",
                                "PBMCSLEctrl"="PBMC-ctrl",
                                "PBMCSLEstim"="PBMC-stim",
                                "PBMCeightkilo"="PBMC8k",
                                "PBMCsevenkilo"="PBMC7k",
                                "Plasschaert"="Plasschaert",
                                "QuakeTrachea"="Quake",
                                "ToschesLizard"="Tosches",
                                "ZeiselBrain"="Zeisel",
                                "ZilionisLung"="Zilionis"
                                
                                
    )}
  return(datain)
}

#' Using heat map to perform the comparison of GeneClust with competing FS 
#' methods in cell clustering on 12 scRNA-seq datasets. 
#' 
#' param md: a data frame of clustering results.
#' param fl: the title of the color bar of rankings.
#' return a ggplot object.
hmap <- function(md,fl,method_input){
  leg_po <- c(-0.43,0.2)
  if(method_input=="TSCAN"){
    leg_po <- c(-0.43,0.2)
  }
  p_hm <- ggplot(md, aes(x = dataset, y = method, fill = rank ))+
    geom_tile()+
    scale_x_discrete(expand = expansion(add=0),
                     position = "bottom") + 
    scale_y_discrete(expand = expansion(add=0),
                     position = "left")+
    scale_fill_gradient( low = "#FCFAF2",
                         high = "#005CAF", 
                         trans="reverse",
                         breaks = c(3,6,9,12))+
    labs(fill = fl, x = "",y = "")+
    geom_text(aes(x = dataset, y = method, label = rank ), 
              color = "#1C1C1C", size = 7/.pt, hjust = 0.5, vjust = 0.5)+
    theme(panel.spacing.x = unit(0.2, "lines"),
          panel.spacing.y = unit(1.2, "lines"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.text.x=element_text(size=7,colour="black",vjust = 1,
                                   hjust = 1,angle = 45),
          axis.text.y=element_text(size=7,colour="black"),
          axis.ticks.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.x = element_line(size = 1/.pt), 
          legend.text=element_text(size=7,colour="black"),
          legend.title=element_text(size=8,colour="black"),
          legend.position = leg_po,
          legend.key.height = unit(12, "pt"),
          legend.key.width = unit(8, "pt"),
          legend.background = element_rect(fill = "transparent"),
          plot.margin = unit(c(0,1,0,25),'pt'))+
    guides(fill = guide_colorbar(order = 1))
  return(p_hm)
}

#' Using box plot to perform the comparison of GeneClust with competing FS 
#' methods in cell clustering on 12 scRNA-seq datasets. 
#' 
#' param mdnum: a data frame of clustering results.
#' param marg: a vector of the plot margins.
#' return a ggplot object.
box_plot_method_md <- function(mdnum,marg){
  p <- ggplot(mdnum, aes(x=method,y=rank))+
    theme(panel.border = element_rect(fill=NA,color="black",size = 1/.pt,linetype="solid"))+
    stat_boxplot(geom="errorbar",width = 0.3,size = 1/.pt)+
    geom_boxplot(fill = "#d8e9f0",fatten = NULL, size = 1/.pt,
                 outlier.size = 1.5/.pt)+
    labs(fill="Methods")+
    stat_summary(fun = mean,geom = "errorbar",aes(ymax = ..y..,ymin = ..y..),
                 width = 0.75, size = 1.2/.pt, linetype = "solid")+
    scale_y_continuous(breaks = c(3,6,9,12))+
    coord_flip()+ theme_bw()+
    theme(axis.text.y=element_blank(),
          axis.text.x=element_text(vjust=0.5,size=7,colour="gray15"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none",
          axis.ticks.x = element_line(size = 1/.pt), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(marg,'pt')
    )
  return(p)
}

#' Using the combination of heat map and box plot to perform the comparison of 
#' GeneClust with competing FS methods in cell clustering on 12 scRNA-seq 
#' datasets. The methods are ranked by the clustering results.
#' 
#' param method_input: a character represents the clustering method.
#' param index: a character represents ARI or V-measure.
#' return one ggplot object.

comparison <- function(method_input,index){
  ## choose data sheet 
  if(method_input=="seurat"){
    usesheet <- "Seurat"
  }else if(method_input=="TSCAN"){
    usesheet <- "TSCAN"
  }else if(method_input=="SC3s"){
    usesheet <- "SC3s"
  }else{
    usesheet <- "kmeans"
  }
  
  ## read data
  data0 <- as.data.frame(read_excel(path, sheet = usesheet))
  colnames(data0) <- c("dataset","gene_number","GiniClust3","deviance","VST",
                       "M3Drop","FEAST","scran","triku","scVI",
                       "SC3","scmap","GeneClust-fast","GeneClust-ps")
  
  if(method_input=="seurat"|method_input=="TSCAN"){
    #fill na
    #dataset
    v <- 2:2
    u <- 1
    k <- 1
    while(k <= dim(data0)[1]/2){
      data0[v,1] <- data0[u,1]
      v <- v+2
      u <- u+2
      k <- k+1
    }
    # data0[,1]
    data1 <- data0[data0[,2]==index,]
    if(method_input=="TSCAN"){
      data1 <- data1[data1$dataset!="Chen",]
    }
    ## converting the data format
    md2000 <- melt(data1,id.vars=colnames(data1)[-c(1,2)],
                   measure.vars = colnames(data1)[-c(1,2)],
                   variable.name = "method",value.name = "median")
    md2000 <- md2000[,(dim(md2000)[2]-1):(dim(md2000)[2])]
    md2000[,'dataset'] <- data1[,1]
    
  }else{
    colnames(data0) <- c("dataset","gene_number","folds","GiniClust3","deviance",
                         "VST","M3Drop","FEAST","scran","triku","scVI",
                         "SC3","scmap","GeneClust-fast","GeneClust-ps")
    #fill na
    #dataset
    v <- 2:40
    u <- 1
    k <- 1
    while(k <= dim(data0)[1]/40){
      data0[v,1] <- data0[u,1]
      v <- v+40
      u <- u+40
      k <- k+1
    }
    data0[,1]
    #folds
    v <- 2:2
    u <- 1
    k <- 1
    while(k <= dim(data0)[1]/2){
      data0[v,2] <- data0[u,2]
      v <- v+2
      u <- u+2
      k <- k+1
    }
    data0[,2]
    
    data1 <- data0[data0[,3]==index,]
    md2000 <- dfmedian(data1,12)
  }
  
  
  for(i in 1:nrow(md2000)){
    ## modify the data set name
    md2000$dataset[i] <- switch(md2000$dataset[i],
                                "Adam"="Adam",
                                "Chen"="Chen",
                                "Guo"="Guo",
                                "PBMCSLEctrl"="PBMC-ctrl",
                                "PBMCSLEstim"="PBMC-stim",
                                "PBMCeightkilo"="PBMC8k",
                                "PBMCsevenkilo"="PBMC7k",
                                "Plasschaert"="Plasschaert",
                                "QuakeTrachea"="Quake",
                                "ToschesLizard"="Tosches",
                                "ZeiselBrain"="Zeisel",
                                "ZilionisLung"="Zilionis"
                                
                                
    )}
  
  ## method ordering
  for(i in 1:length(unique(md2000$dataset))){
    md2000[md2000$dataset==unique(md2000$dataset)[i],'rank'] <- rank(-md2000[md2000$dataset==unique(md2000$dataset)[i],'median'])
  }
  meanrank2 <- data.frame(matrix(ncol = 2,nrow = length(unique(md2000$method))))
  meanrank2[,1] <- unique(md2000$method)
  for(i in 1:nrow(meanrank2)){
    meanrank2[i,2] <- round(mean(md2000[md2000$method==unique(md2000$method)[i],'rank']),1)
  }
  
  meanrank2$X1 <- factor(meanrank2$X1,levels = arrange(meanrank2,by_group = X2)[,1])
  
  md2000$method <- factor(md2000$method,levels = rev(arrange(meanrank2,by_group = X2)[,1]))
  
  ## add selected features
  s_fea <- read_excel(path,sheet = "selected_features")
  colnames(s_fea)[1] <- "dataset"
  f <- as_tibble(md2000[md2000$method=="GeneClust-fast",
                        c("method","median","dataset","rank" )])%>%
    inner_join(newname(s_fea[,c("dataset","GeneClust-fast")]),by = "dataset")
  colnames(f)[5] <- "gene_num"
  
  p <- as_tibble(md2000[md2000$method=="GeneClust-ps",
                        c("method","median","dataset","rank" )])%>%
    inner_join(newname(s_fea[,c("dataset","GeneClust-ps")]),by = "dataset")
  colnames(p)[5] <- "gene_num"
  
  bar2000 <- md2000[md2000$method!="GeneClust-fast"&md2000$method!="GeneClust-ps",
                    c("method","median","dataset","rank" )]
  bar2000[,"gene_num"] <- 2000
  bar2000 <- rbind(bar2000,f,p)
  ## set text color
  bar2000[bar2000$method!="GeneClust-fast"&bar2000$method!="GeneClust-ps","text_color"] <- "transparent"
  bar2000[bar2000$method=="GeneClust-fast","text_color"] <- "black"
  bar2000[bar2000$method=="GeneClust-ps","text_color"] <- "black"
  ## generate one row of bar plot
  onefig <- function(bar2000,gn){
    
    plot_data <- as.data.frame(as_tibble(bar2000)%>%
                                 filter(dataset==unique(bar2000$dataset)[(gn-1)*4+1]|
                                          dataset==unique(bar2000$dataset)[(gn-1)*4+2]|
                                          dataset==unique(bar2000$dataset)[(gn-1)*4+3]|
                                          dataset==unique(bar2000$dataset)[(gn-1)*4+4]))
    plot_data$dataset <- factor(plot_data$dataset,levels = rev(unique(plot_data$dataset)))
    
    if(method_input== "TSCAN"& gn==4){
      plot_data <- as.data.frame(as_tibble(bar2000)%>%
                                   filter(dataset==unique(bar2000$dataset)[(gn-1)*4+1]|
                                            dataset==unique(bar2000$dataset)[(gn-1)*4+2]|
                                            dataset==unique(bar2000$dataset)[(gn-1)*4+3]))
    }
    le <- arrange(meanrank2,by_group = X2)[,1]
    plot_data$method <- factor(plot_data$method,
                               levels = arrange(meanrank2,by_group = X2)[,1])
    bp <- rev(brewer.pal(12,"Set3"))
    
    bar_col <- ifelse(le == "GeneClust-ps" ,bp[12],
                      ifelse(le == "GeneClust-fast", bp[7],
                             ifelse(le == "GiniClust3", bp[1],
                                    ifelse(le == "deviance" , bp[2],
                                           ifelse(le == "VST" , bp[3],
                                                  ifelse(le == "M3Drop" , bp[4],
                                                         ifelse(le == "FEAST" , bp[5],      
                                                                ifelse(le == "scran" , bp[6],       
                                                                       ifelse(le == "triku" , bp[11],
                                                                              ifelse(le == "scVI" , bp[8],
                                                                                     ifelse(le == "SC3" , bp[9],bp[10]
                                                                                     )))))))))))

    if(min(plot_data$median)>0.3&min(plot_data$median)<=0.45){

      plot_data$median <- plot_data$median-0.2
      y_breaks <- c(0.0,0.2,0.4,0.6)
      y_labels <- c(0.2,0.4,0.6,0.8)
    }else if(min(plot_data$median)>0.45){

      plot_data$median <- plot_data$median-0.4
      y_breaks <- c(0.0,0.2,0.4)
      y_labels <- c(0.4,0.6,0.8)
    }else{
      plot_data$median <- plot_data$median
      y_breaks <- c(0.0,0.2,0.4,0.6,0.8)
      y_labels <- c("0.0",0.2,0.4,0.6,0.8)
    }
    
    ylim <- c(0,range(plot_data$median)[2]+0.3)
    
    plot_data$text_hjust <- 0
    plot_data$text_vjust <- 0
    for(i in 1:length(unique(plot_data$dataset))){
    mps <- plot_data[plot_data$dataset==unique(plot_data$dataset)[i]&
                       plot_data$method=="GeneClust-ps","median"]
    mfa <- plot_data[plot_data$dataset==unique(plot_data$dataset)[i]&
                       plot_data$method=="GeneClust-fast","median"]
    if((which(arrange(meanrank2,by_group = X2) =="GeneClust-ps")-
       which(arrange(meanrank2,by_group = X2) =="GeneClust-fast"))==-1
       &(mps-mfa)<0)
    { plot_data[plot_data$dataset==unique(plot_data$dataset)[i]&
                  plot_data$method=="GeneClust-ps",
                "text_hjust"] <- -0.15
      plot_data[plot_data$dataset==unique(plot_data$dataset)[i]&
                plot_data$method=="GeneClust-ps",
                "text_vjust"] <- -0.2
      plot_data[plot_data$dataset==unique(plot_data$dataset)[i]&
                  plot_data$method=="GeneClust-fast",
                "text_vjust"] <- 0.25
    
    }
    }
  
    if(gn==1){
      p_margin <- c(0,0,0,0)
    }
    if(gn==2){
      p_margin <- c(0,0,0,0)
    }
    if(gn==3){
      p_margin <- c(0,0,0,0)
    }
    
    if(index=="ARI"){
      show_index = "ARI"
    }
    if(index=="V"){
      show_index = "V-measure"
    }
    

    p_one <- ggplot(plot_data,aes(x = dataset, y = median, fill = method))+ 
      geom_bar(aes(x = dataset, y = median, fill = method),
               stat='identity',position=position_dodge(.97),width = 0.97)+
      geom_text(aes(x = dataset, y = median,label = gene_num),
                position=position_dodge(.97),size = 6.3/.pt,color = plot_data$text_color,
                hjust = plot_data$text_hjust,vjust = plot_data$text_vjust,angle = 60)+
      scale_y_continuous(expand = c(0.00001, 0),
                         limits = ylim,breaks = y_breaks,
                         labels = y_labels)+
      scale_x_discrete(expand = expansion(add=0.5))+
      labs(x = " ")+
      labs(y = show_index,fill = "Method")+
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "transparent",colour = NA),
            axis.text.x=element_text(vjust=0.5,size=7,colour="black"),
            axis.text.y=element_text(vjust=0.5,size=7,colour="black"),
            axis.title.y=element_text(vjust=0.5,size=7,colour="black"),
            axis.ticks = element_line(size = 1/.pt),
            legend.text=element_text(size=7,colour="black"),
            legend.title=element_text(size=8,colour="black"),
            legend.position ="left",
            legend.key.height = unit(12, "pt"),
            legend.key.width = unit(10, "pt"),
            plot.margin =unit(p_margin,'pt'),
            legend.box.margin = unit(c(3,0,3,3),"pt")
      )+scale_fill_manual(values = bar_col)
    return(p_one)
  }
  
  ## composite plot
  hm2000 <- hmap(md2000,"Rank",method_input)
  ar2000 <- ggarrange(hm2000,box_plot_method_md(md2000,c(0,0,24,0)),
                      ncol = 2, nrow = 1,widths = c(1,0.2))
  com <- onefig(bar2000,1)/onefig(bar2000,2)/onefig(bar2000,3)+
    plot_layout(guides = "collect")&
    theme(plot.title = element_text(colour  = "black", size = 8,hjust = 0.5), 
          legend.position = "left",
          legend.background = element_rect(fill = "white", colour = "white" ),
          legend.key = element_rect(fill = "white", colour = "white" ),
          plot.background = element_rect(fill = "white", colour = "white"),
          legend.key.height = unit(12, "pt"),
          legend.key.width = unit(10, "pt"),
          legend.box.margin = unit(c(3,0,3,3),"pt"))
  arfi  <- plot_spacer()+ar2000+
    plot_layout(ncol = 2,widths = c(0.068,1))
  
  if(method_input=="TSCAN"){
    arfi  <- plot_spacer()+ar2000+
      plot_layout(ncol = 2,widths = c(0.068,1))
  }
  final <- ggarrange(com,arfi,
                     ncol = 1,nrow = 2,heights = c(1,1.1),legend = "right",
                     labels = c("A","B"),vjust = 1.2,hjust = -0.2,
                     font.label = list(size = 10, color = "black", 
                                       face = "bold", family = NULL))+
    theme(legend.background = element_rect(fill = "white", colour = "white" ),
          legend.key = element_rect(fill = "white", colour = "white" ),
          plot.background = element_rect(fill = "white", colour = "white"),
          plot.margin = unit(c(2,2,0,2),"pt"))
  
  return(final)
}


#generate figures####
suppressWarnings(ggsave(filename = paste0("FIG2",".jpeg"),
                        comparison("seurat","ARI"),
                        dpi = 500,quality = 100,
                        height = 6.5,width = 5.512,units = "in"))

suppressWarnings(ggsave(filename = paste0("FIGS2",".jpeg"),
                        comparison("TSCAN","ARI"),
                        dpi = 500,quality = 100,
                        height = 6.5,width = 5.512,units = "in"))
suppressWarnings(ggsave(filename = paste0("FIGS1",".jpeg"),
                        comparison("SC3s","ARI"),
                        dpi = 500,quality = 100,
                        height = 6.5,width = 5.512,units = "in"))
suppressWarnings(ggsave(filename = paste0("FIGS3",".jpeg"),
                        comparison("kmeans","ARI"),
                        dpi = 500,quality = 100,
                        height = 6.5,width = 5.512,units = "in"))

suppressWarnings(ggsave(filename = paste0("FIGS4",".jpeg"),
                        comparison("seurat","V"),
                        dpi = 500,quality = 100,
                        height = 6.5,width = 5.512,units = "in"))
suppressWarnings(ggsave(filename = paste0("FIGS6",".jpeg"),
                        comparison("TSCAN","V"),
                        dpi = 500,quality = 100,
                        height = 6.5,width = 5.512,units = "in"))
suppressWarnings(ggsave(filename = paste0("FIGS5",".jpeg"),
                        comparison("SC3s","V"),
                        dpi = 500,quality = 100,
                        height = 6.5,width = 5.512,units = "in"))
suppressWarnings(ggsave(filename = paste0("FIGS7",".jpeg"),
                        comparison("kmeans","V"),
                        dpi = 500,quality = 100,
                        height = 6.5,width = 5.512,units = "in"))


