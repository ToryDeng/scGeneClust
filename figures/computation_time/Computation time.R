#packages#######
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(grid)))

#paths#####
computation_time_path <- "./computation_time.xlsx"

#function####
#' Constructing the data for generating plots.
#'
#' param timex: a data frame of genes or cells data extracted from the raw data.
#' return a data frame of transformed genes or cells data.
transdata <- function(timex,time0){
  names(timex) <- c("number",colnames(time0)[-1])
  timex2 <- melt(timex ,id.vars=colnames(time0)[-1],
                 measure.vars = colnames(time0)[-1],
                 variable.name = "Methods",value.name = "time")
  
  #dim(timex2)
  timex3 <- timex2[,(dim(timex2)[2]-1):(dim(timex2)[2])]
  timex3[,3] <- timex[,1]
  return(timex3)
}

#' Performing the evaluation of computational efficiency of GeneClust and 11 
#' competing methods. Generating one subplot on cell or gene level.
#' 
#' param datain: a data frame of transformed genes or cells data.
#' param xlab: a character of the label of x-axis.
#' param marginp: the margin of one subplot.
#' return a ggplot object.
onep <- function(datain,xlab,marginp){
  
  bp <- rev(brewer.pal(12,"Set3"))
  my_col <- c(bp[12],bp[5],bp[4],bp[2],bp[10],
              bp[6],bp[11],bp[1],bp[7],
              bp[3],bp[8],bp[9])
  
  p <- ggplot(datain, aes(x=V3, y=time, group = Methods, 
                          color=Methods,shape=Methods))+ 
    geom_point(size=6/.pt,stroke = 2/.pt)+geom_line(size = 1.2/.pt)+
    scale_color_manual(values = my_col)+
    scale_shape_manual(values = c(0,15,1:2,4:6,8:10,12,3))+
    theme_bw()+ labs(x = xlab, y = "Computation time (s)" )+
    theme(axis.title.y=element_text(vjust = 0.5,size = 7,color = "black"),  
          axis.title.x=element_text(size = 7.2,color = "black"),
          axis.text.y=element_text(size = 7,vjust = 0.5,color = "black"),            
          axis.text.x=element_text(size = 7,color = "black",vjust = 0.5,angle = 20),
          axis.ticks = element_line(size = 1/.pt,color = "black"), 
          plot.title=element_text(hjust=0.5, size=8),
          legend.title = element_text(color = "black", size = 7.5),
          legend.text = element_text(color = "black", size = 7),
          legend.position = "right",
          legend.key.height = unit(10, "pt"),
          panel.grid= element_blank(),
          plot.margin =marginp,
          legend.key.width = unit(6, "pt")
    )
  
  return(p)
}

#' Performing the evaluation of computational efficiency of GeneClust and 11 
#' competing methods. Generating a composite plot.
#' 
#' param time0: a data frame of the raw result that contains computational efficiency
#'              of GeneClust and 11 competing methods on cell or gene level.
#' return a ggplot object.    
Computation_time <- function(time0){
  ## transform the format of gene data
  tg <- time0[1:4,]
  options(digits = 3)
  
  tg[,1] <- c(10000,15000,20000,25000)
  ntg <- transdata(tg,time0)
  ntg$Methods <- as.character(ntg$Methods)
  ntg[ntg$Methods=="Deviance",'Methods'] <- "deviance"
  ntg[ntg$Methods=="Variance",'Methods'] <- "scVI"
  ntg[ntg$Methods=="Seurat v3",'Methods'] <- "VST"
  ntg$Methods <- factor(ntg$Methods,levels = c("GeneClust-ps","FEAST","M3Drop","deviance","scmap",
                                                  "scran","triku","GiniClust3","GeneClust-fast",
                                                  "VST","scVI","SC3"))
  ## transform the format of cell data
  tc <- time0[6:9,]
  tc[,1] <- c(10000,1000,50000,5000)
  ntc <- transdata(tc,time0)
  ntc$Methods <- as.character(ntc$Methods)
  ntc[ntc$Methods=="Deviance",'Methods'] <- "deviance"
  ntc[ntc$Methods=="Variance",'Methods'] <- "scVI"
  ntc[ntc$Methods=="Seurat v3",'Methods'] <- "VST"
  ntc$Methods <- factor(ntc$Methods,levels = c("GeneClust-ps","FEAST","M3Drop","deviance","scmap",
                                               "scran","triku","GiniClust3","GeneClust-fast",
                                               "VST","scVI","SC3"))
  
  genes <- onep(ntg,"Number of genes",unit(c(10,2,5,3),'pt'))+
           scale_y_log10(expand = c(0.01, 0.2),breaks = c(1,10,100,1000,10000),
                         labels = c("1","10","100","1000","10000"))+
           scale_x_continuous(expand = c(0.065,0),
                         labels = c("10000","15000","20000","25000"))
  
  
  cells <- onep(ntc,"Number of cells",unit(c(10,2,3,5),'pt'))+
           scale_y_log10(expand = c(0.01, 0.2),breaks = c(0.1,1,10,100,1000),
                         labels = c("0.1","1","10","100","1000"))+
           scale_x_log10(expand = c(0.065,0),breaks = c(1000,5000,10000,50000),
                         labels=c(expression(1%*%10^3),
                                  expression(5%*%10^3),
                                  expression(1%*%10^4),
                                  expression(5%*%10^4)))+
           theme(legend.position = "none")
   
  pout <- ggarrange(cells,genes,
            ncol = 2,nrow = 1,widths = c(1,1),
            common.legend = TRUE, legend = "right",
            labels = c("A","B"),vjust = 1.2,hjust = -0.2,
            font.label = list(size = 8, color = "black", 
                              face = "bold", family = NULL))+
            theme(legend.background = element_rect(fill = "white", colour = "white" ),
            legend.key = element_rect(fill = "white", colour = "white" ),
            plot.background = element_rect(fill = "white", colour = "white"),
            plot.margin = unit(c(3,2,2,2),'pt'))
  
  return(pout)
}

#generate figure####
time0 <- as.data.frame(read_excel(computation_time_path))
Computation_time(time0)
ggsave(filename = "FIG3.jpeg",
       Computation_time(time0),
       dpi = 500,quality = 100,
       height = 3,width = 5.512,units = "in")
