#####################################NEW R SCRIPT TO ANALYZE SCZ MONOCYTE DATA##############################################

library(xlsx)
library(dplyr)
library(plyr)
library(readxl)
library(Hmisc)
library(DESeq2)
library(factoextra)
library(broom)

library(variancePartition)
library(RColorBrewer)

library(ggpubr)
library(tidyr)
library(ggeasy)
library(ggrepel)
library(gridExtra)
library(ComplexHeatmap)

#library(pcaExplorer)
#library(vsn)

#==========DEGs expression===================#
prev.degs <- c('IL1B', 'TNF', 'IL6', 'TNFAIP3', 'CCL7', 'CCL20', 'CXCL2', 'CCL2', 'CXCL3', 'MAFF', 
               'F3', 'SERPINB2', 'THBS', 'EREG','DUSP2', 'MAPK6', 'PDE4B', 'BCL2A1', 'PTX3', 'PTGS2', 'CDC42',
               'TREM1', 'ATF3', 'EGR3', 'MXD1', 'SPI1',
               'TLR2','TLR3','TLR4','TLR5')
DEG.expr <- as.matrix(batch.rem[rownames(batch.rem) %in% prev.degs,])
CTR.DEGexpr <- DEG.expr[,colnames(DEG.expr) %in% rownames(rev.covs[rev.covs$Status == 'CTR',])]
SCZ.DEGexpr <- DEG.expr[,colnames(DEG.expr) %in% rownames(rev.covs[rev.covs$Status == 'SCZ',])]
identical(rownames( res[rownames(res) %in% prev.degs,]), rownames(CTR.DEGexpr))
DEG.tab <- round(cbind(rowMeans(CTR.DEGexpr),rowSds(CTR.DEGexpr),
            rowMeans(SCZ.DEGexpr),rowSds(SCZ.DEGexpr),
            res[rownames(res) %in% prev.degs,]$log2FoldChange,
            res[rownames(res) %in% prev.degs,]$pvalue,
            res[rownames(res) %in% prev.degs,]$padj),2)

colnames(DEG.tab) <- c("CTR mean", "CTR std", "SCZ mean", "SCZ std", "l2FC", 'nominal P', 'FDR')

DEG.tab2 <- matrix(nrow=length(rownames(DEG.tab)), ncol=5)
DEG.tab2[,1] <- paste(DEG.tab[,1],paste('(',DEG.tab[,2],')', sep=""))
DEG.tab2[,2] <- paste(DEG.tab[,3],paste('(',DEG.tab[,4],')', sep=""))
DEG.tab2[,c(3:5)] <- DEG.tab[,c(5:7)]

DEG.tab2 <- data.frame(DEG.tab2)
colnames(DEG.tab2) <- c('CTR mean (sd)','SCZ mean (sd)', 'l2FC', 'nominal P', 'FDR')
rownames(DEG.tab2) <- rownames(DEG.tab)
DEG.tab2 <- DEG.tab2[order(DEG.tab2$l2FC,decreasing=T),]

write.xlsx(DEG.tab2,file="03312022_DEG table_absolute expr.xlsx")

EXPR <- data.frame('Expression'=c(rowMeans(CTR.DEGexpr),rowMeans(SCZ.DEGexpr)),
                   'Status'=as.factor(c(rep('CTR',length(rowMeans(CTR.DEGexpr))),rep('SCZ',length(rowMeans(CTR.DEGexpr))))),
                   'paired'=as.numeric(as.factor(rep(rownames(CTR.DEGexpr),2))))

ggplot(EXPR,aes_string(x='Status',y='Expression', fill='Status')) +
  geom_boxplot() +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Status,group=paired),size=2,shape=21, position = position_dodge(0.2)) +
  theme_classic()


#==========Paper figures=====================#
#a) volcano plot
# Function Vulcano plot 
volcano_plot <- function(res, title = NULL, subtitle = NULL, annotate_by = NULL, type ='ALS'){
  res <- 
    mutate(res,
           sig = case_when(
             padj >= 0.1 ~ "non_sig",
             padj < 0.1 & abs(log2FoldChange) < 0.5~ "sig",
             padj < 0.1 & abs(log2FoldChange) >= 0.5 ~ "sig - strong"
           )) %>%
    mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
    mutate(log2FoldChange = case_when(
      log2FoldChange > 5 ~ Inf,
      log2FoldChange < -5 ~ -Inf,
      TRUE ~ log2FoldChange
    )) %>%
    mutate(class = paste(sig, direction))
  if( type == "ALS"){
    xpos <- 0.5
    ymax <- 7.5
    xlim <- c(-4.5,4.5)
  }else{
    xpos <- 0.025
    ymax <- 8.5
    xlim <- c(-0.042, 0.042)
  }
  de_tally <- group_by(res, sig, direction, class) %>% tally() %>%
    filter(sig != "non_sig") %>%
    mutate( position = ifelse(sig == "sig", xpos, 2) ) %>%
    mutate(position = ifelse( direction == "down", -1 * position, position)) %>%
    mutate(n = formatC(n, format="f", big.mark=",", digits=0))
  plot <- res %>%
    mutate( pvalue = ifelse( pvalue < 1e-90, Inf, pvalue)) %>% #threshold at 1e16
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) + 
    #geom_point(aes(colour = class ), size = 0.5) +
    rasterise(geom_point(aes(colour = class ), size = 0.5), dpi = 300) +
    scale_colour_manual(values = c("non_sig up" = "gray", 
                                   "non_sig down" = "gray",
                                   "sig up" = "#EB7F56", 
                                   "sig - strong up" = "#B61927",
                                   "sig down" = "#4F8FC4",
                                   "sig - strong down" = "dodgerblue4"
    )) +
    theme_bw() +
    labs(y = expression(-log[10]~P~value), x = expression(log[2]~"(fold change)"), title = title, subtitle = subtitle) +
    guides(colour = FALSE) +
    scale_y_continuous(expand = c(0,0), limits = c(0,ymax)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.border = element_blank(),
      axis.ticks = element_line(colour = "black")
    ) +
    geom_text(fontface = "bold", data = de_tally, aes(x = position, y = ymax - 0.5, label = n, colour = class), size = 4) +
    scale_x_continuous(limits = xlim)
  if(!is.null(annotate_by)){
    plot <- plot + 
      ggrepel::geom_text_repel(
        fontface = "italic",
        data = filter(res, symbol %in% annotate_by), 
        aes(x = log2FoldChange, y = -log10(pvalue), label = symbol), 
        min.segment.length = unit(0, "lines"),
        size = 2.3) +
      geom_point(
        data = filter(res, symbol %in% annotate_by), size = 0.8, colour = "black"
      ) +
      geom_point(aes(colour = class ),
                 data = filter(res, symbol %in% annotate_by), size = 0.6
      )
  }
  return(plot)
}

res.wnam <- data.frame(res)
res.wnam$symbol <- rownames(res) #need a column with gene names termed 'symbol'

library(ggrastr)
volcano_plot(res.wnam, annotate_by=c(upreg.genes,downreg.genes, c('FES',"IL6", "IL7", "INF2", "P2RY6", "TNFAIP3", "CSMD1", "MMP9", "NFKBIZ", "NFKBIA")))


# "NFKBIZ", "NFKBIA",
#b) DEG heatmap
ordered.res <- res[order(abs(res$log2FoldChange), decreasing = T),]

top20 <- rownames(ordered.res[abs(ordered.res$log2FoldChange) > 0 & ordered.res$padj < 0.1,][1:20,])
top20.dat <- batch.rem[rownames(batch.rem) %in% top20,]

rev.covs <- annot_df[order(length(rownames(annot_df)):1),]
top20.dat <- top20.dat[,order(ncol(top20.dat):1)]
identical(colnames(top20.dat), rownames(rev.covs))

ha <- HeatmapAnnotation(Status = rev.covs$Status,
                  col = list( Status = c("CTR" = "lightblue", "SCZ" = "tomato")))
pushViewport(viewport(width = 0.8, height = 0.8))
new.ha <- packLegend(ha, direction = "horizontal")
draw(new.ha, x = unit(1, "mm"), y = unit(2, "mm"), just = c("left", "bottom"))
ComplexHeatmap::Heatmap(t(scale(t(top20.dat))),
                        cluster_rows = T,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = F,
                        column_names_rot = 90, 
                        top_annotation = ha,
                        name = "Normalized counts",
                        clustering_method_rows = "ward.D2",
                        clustering_distance_rows = "pearson"
)


#c) GO enrichment

dotplot(GO.up, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

#d) gene set enrichment (GSEA)
ComplexHeatmap::Heatmap(DEG.enrich,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30),
                        name = "log2(q)",

)
#e) selection of lists 
clFC.select <- data.frame("compound lFC"=compoundlFC[rownames(compoundlFC) %in% c("NFkB signaling", "LPS response", "Glucocorticoid response"),])
ComplexHeatmap::Heatmap(clFC.select,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        show_column_names = F,
                        col = colors,
                        name = "Mean log(FC)",
                        
)


#f) PC boxplots of LPS & NfKB pca (throw in SCZ)
dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"NFkB signaling"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
t.test(x=list.pca[list.pca$Status %in% 'CTR',]$PC1,
       y=list.pca[list.pca$Status %in% 'SCZ',]$PC1, alternative='two.sided')


#g) GC PCA plot and boxplot
dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"Glucocorticoid response"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
t.test(PC2 ~ Status, data=list.pca)

#h) NFkB & LPS co-regulation with other gene lists
pca.list <- matrix(nrow=length(colnames(batch.rem)), ncol=length(colnames(gene.panels_mrgd)))
colnames(pca.list) <- colnames(gene.panels_mrgd)
for (i in colnames(gene.panels_mrgd)){
  dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,i]),]
  #  plots <- PCA_cov_cor_R(clean_covs, dat) #if 
  pca.list[,i] <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                             ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)[,1]}

z.nfkpblps <- P[,c(3,5,7)]
z.nfkpblps[z.nfkpblps == " NA"] <- "Inf"
ComplexHeatmap::Heatmap(cor(pca.list)[,c("NFkB signaling","LPS response","Glucocorticoid response")],
                        cluster_rows = T,
                        cluster_columns = T,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = T,
                        column_names_rot = 90, 
                        name = "Pearson correlation",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(as.matrix(z.nfkpblps)[i,j], x, y, gp = gpar(fontsize = 10))
                        }
                        #        clustering_method_rows = "ward.D2",
                        #        clustering_distance_rows = "pearson"
)


#=====================lFC correlation========================#
setwd("/Users/kubler01/Documents/R projects/Monocyte project/Gene lists")
Drexhage <- read_excel("Drexhage SCZ DEG lFC.xlsx")
Gandal <- read_excel("Gandal DEG lFC.xlsx")
Gandal <- Gandal[Gandal$SCZ.fdr < 0.1,]
MIGA <- read_excel("SCZ MIGA QTL genes.xlsx")
MIGA$QTL_Gene

table(MIGA$QTL_Gene %in% c(DEGs$`Up-regulated`,DEGs$`Down-regulated`))

log10(MIGA[MIGA$QTL_Gene %in% "FES",]$QTL_P)

MIGA.FES <- as.matrix(data.frame("Beta"=MIGA[MIGA$QTL_Gene %in% "FES",]$QTL_Beta,check.names=F))

rownames(MIGA.FES) <- MIGA[MIGA$QTL_Gene %in% "FES",]$QTL_SNP
ra.FES <- HeatmapAnnotation(`QTL type` = c(rep('Microglia sQTL',5),rep('Monocyte eQTL',2)), `Origin`=c(rep('MIGA MFG',2),rep("MIGA STG",3),'Fairfax','MyND'),
                        col = list(`QTL type` = c("Monocyte eQTL" = "#C71000FF", "Microglia sQTL" = "#008EA0FF"), 
                                   `Origin`=c("MIGA MFG"="#E66101","MIGA STG"="#FDB863","Fairfax"="#B2ABD2","MyND"="#5E3C99")))#colorRampPalette(RColorBrewer::brewer.pal(4,"PuOr"))(4)
S2c <- ComplexHeatmap::Heatmap(t(MIGA.FES[,1]),
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        show_column_names = F,
                        col = colors,
                        top_annotation = ra.FES,
                        column_names_rot = 60, 
                        name = "Beta",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(t(formatC(MIGA[MIGA$QTL_Gene %in% "FES",]$QTL_P, format = "e", digits = 2))[i,j], x, y, gp = gpar(fontsize = 10))
                        }
)


#For Drexhage FCs
Drexhage.lfc <- drop_na(left_join(Drexhage,
                                  data.frame("Gene"=rownames(res[rownames(res) %in% Drexhage$Gene,]),
                                             'Monocytes' = res[rownames(res) %in% Drexhage$Gene,]$log2FoldChange, check.names=F), by='Gene'))

#correlation plot of log FCs
symbols.label.drex <- data.frame('overlap'=Drexhage.lfc$Gene)
symbols.label.drex$overlap[(!symbols.label.drex$overlap %in% nomDEGs)] = ""

Drexhage.lfc$lFC <- log2(Drexhage.lfc$lFC)
lFC.scat[['Drex']] <- ggplot2::ggplot(Drexhage.lfc, aes(x = lFC, y = Monocytes)) + #Drexhage genes are in FC, not lFC
  #geom_point(alpha = .5) +
  scale_color_manual(values = c("#C71000FF", "#8A4198FF", "#008EA0FF")) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") + 
  stat_smooth(method = "lm", se=F) + # Add Regression Line 
  stat_regline_equation(aes(label = ..adj.rr.label..), show.legend = F)  + # Add R-Square
  #  stat_regline_equation(aes(label = ..rr.label..))  +
  ggrepel::geom_text_repel(
    fontface = "italic",
    data = filter(Drexhage.lfc, Gene %in% symbols.label.drex$overlap), 
    aes(x = lFC, y = Monocytes, label = Gene), 
    min.segment.length = unit(0, "lines"),
    size = 2) +
  geom_point(
    data = Drexhage.lfc, size = 1, colour = 'darkgrey'
  ) +
  geom_point(
    data = filter(filter(Drexhage.lfc, Monocytes > 0), Gene %in% nomDEGs), size = 1, colour = rep("tomato",6)
  ) +
  geom_point(
    data = filter(filter(Drexhage.lfc, Monocytes < 0), Gene %in% nomDEGs), size = 1, colour = rep('lightblue',3)
  )+
  geom_point(
    data = filter(Drexhage.lfc, Gene %in% c('IL6','TNFAIP3','DUSP2')), size = 1, colour = rep('#B61927',3)
  ) +
  easy_labs(x = expression(paste("Monocyte fold change (Drexhage)")), y = paste("Monocyte fold change (this study)")) +
  
  theme_classic()
lFC.scat[['Drex']] 


#For Gandal FCs
Gandal.lfc <- drop_na(left_join(Gandal[Gandal$gene_name %in% c(upreg.genesALL,downreg.genesALL),],
                                data.frame("gene_name"=rownames(res[rownames(res) %in% Gandal$gene_name,]),
                                           'Monocytes' = res[rownames(res) %in% Gandal$gene_name,]$log2FoldChange, 
                                           check.names=F), by='gene_name'))


symbols.label <- data.frame('overlap'=Gandal.lfc$gene_name)
symbols.label$overlap[(!symbols.label$overlap %in% c(upreg.genesALL, downreg.genesALL))] = ""

lFC.scat[['Gandal']] <- ggplot2::ggplot(Gandal.lfc, aes(x = SCZ.log2FC, y = Monocytes)) + #Drexhage genes are in FC, not lFC
  geom_point(alpha = .5) +
  scale_color_manual(values = c("#C71000FF", "#8A4198FF", "#008EA0FF")) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") + 
  stat_smooth(method = "lm", se=F) + # Add Regression Line 
  stat_regline_equation(aes(label = ..adj.rr.label..), show.legend = F)  + # Add R-Square
  #  stat_regline_equation(aes(label = ..rr.label..))  +
  ggrepel::geom_text_repel(
    fontface = "italic",
    data = filter(Gandal.lfc, gene_name %in% c(upreg.genesALL, downreg.genesALL)), 
    aes(x = SCZ.log2FC, y = Monocytes, label = gene_name), 
    min.segment.length = unit(0, "lines"),
    size = 2) +
  geom_point(
    data = Gandal.lfc, size = 1, colour = "darkgrey"
  ) +
  geom_point(
    data = filter(Gandal.lfc, gene_name %in% upreg.genesALL), size = 1, colour = rep("tomato",17)
  ) +
  geom_point(
    data = filter(Gandal.lfc, gene_name %in% downreg.genesALL), size = 1, colour = rep('lightblue',14)
  )+
  geom_point(
    data = filter(Gandal.lfc, gene_name %in% c('IL6','TNFAIP3','DUSP2')), size = 1, colour = rep('#B61927',1)
  ) +
  easy_labs(x = expression(paste("Brain expression fold change (Gandal)")), y = paste("Monocyte fold change (this study)")) +
  theme_classic()



#==========PDF paper figure============#
setwd(file.path("/Users/kubler01/Documents/R projects/Monocyte project/Paper"))

ra <- ComplexHeatmap::rowAnnotation(Status = rev.covs$Status,
                                    col = list(Status = c("CTR" = "lightblue", "SCZ" = "tomato")))
checkGO<-GO.up
checkGO@result<-checkGO@result[checkGO@result$p.adjust < 0.1 & checkGO@result$ONTOLOGY %in% "BP",]




#Boxplots for PC expression
dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"NFkB signaling"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
#PCAplot(list.pca, "Status", PoV.df=list.pov, pc.1=1, pc.2=2, colors=c("lightblue", "tomato"), geom.size=4)
boxplots<-list()
boxplots[['NFkB']] <- ggplot(list.pca[,c(1,4)],aes_string(x = "Status", y = "PC1", fill = "Status")) +
  geom_boxplot(fill=c('lightblue','tomato')) +
  geom_point(aes(fill=Status),size=2,shape=21, position = position_dodge(0.2)) +
  theme_classic()+
  theme(plot.title = element_text(size = 12, face = "bold"))+
  labs(title='NFkB pathway',y="PC1 expression")+
  scale_fill_manual(values = c('lightblue','tomato'))+
  stat_compare_means(aes(label = ..p.signif..),method = "t.test", label.x = 1.5, label.y=7)


dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"LPS response"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
#PCAplot(list.pca, "Status", PoV.df=list.pov, pc.1=1, pc.2=2, colors=c("lightblue", "tomato"), geom.size=4)
boxplots[['LPS']] <- ggplot(list.pca[,c(2,4)],aes_string(x = "Status", y = "PC2", fill = "Status")) +
  geom_boxplot(fill=c('lightblue','tomato')) +
  geom_point(aes(fill=Status),size=2,shape=21, position = position_dodge(0.2)) +
  theme_classic()+
  theme(plot.title = element_text(size = 12, face = "bold"))+
  labs(title='LPS response',y="PC2 expression")+
  scale_fill_manual(values = c('lightblue','tomato'))+
  stat_compare_means(aes(label = ..p.signif..),method = "t.test", label.x = 1.5, label.y=8)


dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"Glucocorticoid response"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
#PCAplot(list.pca, "Status", PoV.df=list.pov, pc.1=1, pc.2=2, colors=c("lightblue", "tomato"), geom.size=4)
boxplots[['GC']] <- ggplot(list.pca[,c(2,4)],aes_string(x = "Status", y = "PC2", fill = "Status")) +
  geom_boxplot(fill=c('lightblue','tomato')) +
  geom_point(aes(fill=Status),size=2,shape=21, position = position_dodge(0.2)) +
  theme_classic()+
  theme(plot.title = element_text(size = 12, face = "bold"))+
  labs(title='Glucocorticoid response',y="PC2 expression")+
  scale_fill_manual(values = c('lightblue','tomato'))+
  stat_compare_means(aes(label = ..p.signif..),method = "t.test", label.x = 1.5, label.y=7)


p6 <- ggarrange(plotlist=boxplots[1:length(boxplots)], common.legend = T, ncol = 1)




## 2x3 grid of plotting areas
#Plot A
library(ggrastr)
p1 <- volcano_plot(res.wnam, 
                   annotate_by=c(upreg.genes,downreg.genes, c('FES',"IL6", "IL7", "INF2", "P2RY6", "TNFAIP3", "CSMD1", "MMP9", "NFKBIZ", "NFKBIA"))) 

#Plot B
lfc.grob <- list(grid.grabExpr(plot(lFC.scat$Gandal)),grid.grabExpr(plot(lFC.scat$Drex)))
p2 <- ggarrange(plotlist=lfc.grob[1:length(lfc.grob)], nrow = 2)
#Plot C
#top20.dat
DEG.expr <- DEG.expr[,order(ncol(DEG.expr):1)]
identical(colnames(DEG.expr), rownames(rev.covs))

color.genes <- rep('black',length(rownames(DEG.expr)))
color.genes[which(rownames(DEG.expr)%in%upreg.genesALL)] <- 'tomato'

Legend(col_fun = col_fun, title = "foo", at = c(0, 0.5, 1), 
       labels = c("low", "median", "high"))
p3 <- ComplexHeatmap::Heatmap(t(scale(t(DEG.expr))),
                        cluster_rows = T,
                        cluster_columns =F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = F,
                        column_names_rot = 60, 
                        #right_annotation  = ra,
                        top_annotation = ha,
                        name = "row Z score",
                        clustering_method_rows = "ward.D2",
                        clustering_distance_rows = "pearson", 
                        row_names_gp = grid::gpar(fontsize = 6, col = color.genes),
                        column_names_gp = grid::gpar(fontsize = 5))


#Boxplots of three major genes:
CTR.DEGexpr <- DEG.expr[,colnames(DEG.expr) %in% rownames(rev.covs[rev.covs$Status == 'CTR',])]
SCZ.DEGexpr <- DEG.expr[,colnames(DEG.expr) %in% rownames(rev.covs[rev.covs$Status == 'SCZ',])]

box.EXPR <- data.frame(t(DEG.expr[c('TNFAIP3','IL6','DUSP2'),]),
                   'Status'=as.factor(c(rep('CTR',17),rep('SCZ',20))))



DEG.boxs <- list()
for (i in c('TNFAIP3','IL6','DUSP2')){
 DEG.boxs[[i]] <- ggplot(box.EXPR,aes_string(x='Status',y=i, fill='Status')) +
  geom_boxplot(fill=c('lightblue','tomato')) +
  geom_point(aes(fill=Status),size=2,shape=21, position = position_dodge(0.2)) +
  theme_classic()+
  labs(x=NULL, y = 'Residual expression', title = i)+
  scale_fill_manual(values = c('lightblue','tomato'))}

p.boxplots <- ggarrange(plotlist=DEG.boxs[1:length(DEG.boxs)], nrow = 1, common.legend = T)
  

#Plot D
p4 <- dotplot(checkGO, split="ONTOLOGY", showCategory=10) + facet_grid(ONTOLOGY~., scale="free")
#Plot E
n <- expression(log[2]~"(q-value)")
p5 <- ComplexHeatmap::Heatmap(DEG.enrich,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30),
                        column_names_rot = 60, 
                        name = 'log2 of q-value',
                        )





#================Main figures===================#

setwd(file.path("/Users/kubler01/Documents/R projects/Monocyte project/Paper/V3"))


#Figure 1
#Separating the figure into two
lfc.grob <- list(grid.grabExpr(plot(lFC.scat$Gandal)),grid.grabExpr(plot(lFC.scat$Drex)))
p3.grid <- grid.grabExpr(draw(p3,annotation_legend_side = "bottom", 
                         heatmap_legend_side="right"))
p5.grid <- grid.grabExpr(draw(p5))
#p6.grid <- grid.draw(grid.grabExpr(draw(p6)))


lay1 <- rbind(c(1,1,1,2,2,2),
             c(1,1,1,2,2,2),
             c(1,1,1,2,2,2),
             c(3,3,3,3,3,3),
             c(3,3,3,3,3,3),
             c(3,3,3,3,3,3))

grobs1 <- list(p3.grid,lfc.grob[[2]],p.boxplots)
tiff("04072022_Figure 1.tiff", width=1900, height=1800, res=300)
grid.arrange(grobs = grobs1, layout_matrix = lay1)
dev.off()


#Fig 2
lay2 <- rbind(c(2,2,2,1,1,1),
              c(2,2,2,1,1,1),
              c(2,2,2,1,1,1))

grobs2 <- list(lfc.grob[[1]],p1)

tiff("04082022_Figure 2.tiff", width=1900, height=950, res=300)
grid.arrange(grobs = grobs2, layout_matrix = lay2)
dev.off()

#Fig 3
lay3 <- rbind(c(1,1,1,1,1,3,3),
              c(1,1,1,1,1,3,3),
              c(1,1,1,1,1,3,3),
              c(4,4,2,2,2,3,3),
              c(4,4,2,2,2,3,3))

grobs3 <- list(p4,p5.grid,p6,grob(NULL))

tiff("04072022_Figure 3.tiff", width=2700, height=2800, res=300)
grid.arrange(grobs = grobs3, layout_matrix = lay3)
dev.off()

#================Supplementary figures===================#
#Fig S1
annot_df <- data.frame("Status"=as.character(imputed.covs$status), "Sex"=as.character(imputed.covs$sex), check.names = F)
annot_df$Status[annot_df$Status == 1] <- "SCZ"
annot_df$Status[annot_df$Status == 0] <- "CTR"
annot_df$Sex[annot_df$Sex == 1] <- "Male"
annot_df$Sex[annot_df$Sex == 2] <- "Female"
rownames(annot_df) <- rownames(imputed.covs)

ha.IQR <- HeatmapAnnotation(Status = annot_df$Status, Sex = annot_df$Sex,
                        col = list( Status = c("CTR" = "lightblue", "SCZ" = "tomato"), 
                                    Sex = c("Female" = "black", "Male" = "grey")), show_legend = F)
l.IQR <- rowAnnotation(Status = annot_df$Status, Sex = annot_df$Sex,
                            col = list( Status = c("CTR" = "lightblue", "SCZ" = "tomato"), 
                                        Sex = c("Female" = "black", "Male" = "grey")))

IQR <- cor(trnsf.df,trnsf.df)
S1a <- ComplexHeatmap::Heatmap(IQR,
                        cluster_rows = T,
                        cluster_columns = T,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = T,
                        column_names_rot = 90, 
                        top_annotation = ha.IQR,
                        left_annotation = l.IQR,
                        name = "Pearson correlation",
                        clustering_method_rows = "ward.D2",
                        clustering_distance_rows = "euclidean",
                        clustering_method_columns = "ward.D2",
                        clustering_distance_columns =  "euclidean", 
                        column_names_gp = grid::gpar(fontsize = 5),
                        row_names_gp = grid::gpar(fontsize = 5)
)
pca <- plotPCA.custom(as.matrix(trnsf.df), intgroup=c("Status"), 
                      ntop = 50000, returnData=TRUE, metadata=annot_df, pc.1 = 1, pc.2 = 2)
PoV <- round(100 * attr(pca, "percentVar"))
S1b <- PCAplot(pca, "Status", PoV.df=PoV, pc.1 = 1, pc.2 = 2, colors = c('lightblue','tomato'))+
  geom_text(aes(label=ifelse(pca$name %in% c('36','37'),as.character(pca$name),'')),hjust=0,vjust=0) #samples 36 and 37 seem to be outliers

#-after sample removal#
annot_df.rem <- annot_df[!rownames(annot_df) %in% c('36','37'),]

ha.IQRrem <- HeatmapAnnotation(Status = annot_df.rem$Status, Sex = annot_df.rem$Sex,
                            col = list( Status = c("CTR" = "lightblue", "SCZ" = "tomato"), 
                                        Sex = c("Female" = "black", "Male" = "grey")), show_legend = F)
l.IQRrem <- rowAnnotation(Status = annot_df.rem$Status, Sex = annot_df.rem$Sex,
                       col = list( Status = c("CTR" = "lightblue", "SCZ" = "tomato"), 
                                   Sex = c("Female" = "black", "Male" = "grey")), show_legend = F)

IQR.rem <- cor(clean_df,clean_df)
S1c <- ComplexHeatmap::Heatmap(IQR.rem,
                        cluster_rows = T,
                        cluster_columns = T,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = T,
                        column_names_rot = 90, 
                        top_annotation = ha.IQRrem,
                        left_annotation = l.IQRrem,
                        name = "Pearson correlation",
                        clustering_method_rows = "ward.D2",
                        clustering_distance_rows = "euclidean",
                        clustering_method_columns = "ward.D2",
                        clustering_distance_columns =  "euclidean", 
                        column_names_gp = grid::gpar(fontsize = 5),
                        row_names_gp = grid::gpar(fontsize = 5),
                        show_heatmap_legend = T
)

pca.rem <- plotPCA.custom(as.matrix(clean_df), intgroup=c("Status"), 
                      ntop = 50000, returnData=TRUE, metadata=annot_df.rem, pc.1 = 1, pc.2 = 2)
PoV.rem <- round(100 * attr(pca.rem, "percentVar"))
S1d <- PCAplot(pca.rem, "Status", PoV.df=PoV.rem, pc.1 = 1, pc.2 = 2, colors = c('lightblue','tomato')) #samples 36 and 37 seem to be outliers

#PC-covariate correlation
changed.names.covs <- clean_covs
colnames(changed.names.covs)[9:10] <- c('antidepressants','anxiolytics')
colnames(changed.names.covs)[11:13] <- c('Exonic rate', 'Intergenic rate', 'rRNA rate')
colnames(changed.names.covs)[16:18] <- c('Reads mapped', 'Genes detected', "Mean per base coverage")
S1e <- PCA_cov_cor_R(changed.names.covs, clean_df, fs1=12,fs2=3)
S1f <- plotVarPart(vp, label.angle=60)+
  theme(text = element_text(size = 11))


#covariate co-linearity
colinearity <- as.matrix(rcorr(as.matrix(clean_covs))$r)
S1g <- Heatmap(colinearity, 
        cluster_rows = T,
        cluster_columns = T,
        show_row_dend = T,
        show_column_dend = T,
        #show_row_names = T,
        #show_column_names = T,
        name = "Pearson correlation",
        clustering_method_rows = "ward.D2",
        clustering_distance_rows = "euclidean",
        clustering_method_columns = "ward.D2",
        clustering_distance_columns =  "euclidean",
        col = colors, 
        column_names_gp = grid::gpar(fontsize = 11),
        row_names_gp = grid::gpar(fontsize = 11))





#-----Fig S1
layS1 <- rbind(c(1,1,1,1,2,2,2,2),
               c(1,1,1,1,2,2,2,2),
               c(1,1,1,1,2,2,2,2),
               c(1,1,1,1,2,2,2,2),
               c(5,5,5,5,5,5,5,5),
               c(3,3,3,3,4,4,4,4),
               c(3,3,3,3,4,4,4,4),
               c(3,3,3,3,4,4,4,4),
               c(3,3,3,3,4,4,4,4))

S1a.grid <- grid.grabExpr(draw(S1a))
textgrob0 <- grobTree(
  gp = gpar(fontsize = 13), 
  textGrob(label = "Sample distances before outlier removal", name = "title1",
           x = unit(0, "lines"), y = unit(0.2, "lines"), 
           just = "left", hjust = 0, vjust = 0,
           gp = gpar(col = "black", fontface = "bold"))
)
S1a.grid<- arrangeGrob(S1a.grid, top = textgrob0, padding = unit(1, "line"))
S1c.grid <- grid.grabExpr(draw(S1c))
textgrob1 <- grobTree(
  gp = gpar(fontsize = 13), 
  textGrob(label = "Sample distances after outlier removal", name = "title1",
           x = unit(0, "lines"), y = unit(1.3, "lines"), 
           just = "left", hjust = 0, vjust = 0,
           gp = gpar(col = "black", fontface = "bold"))
)
S1c.grid<- arrangeGrob(S1c.grid, top = textgrob1, padding = unit(1, "line"))


S1g.grid <- grid.grabExpr(draw(S1g))

grobsS1 <- list(S1a.grid,S1b,S1c.grid,S1d, grob(NULL))
#list(S1a.grid,S1b,S1c.grid,S1d,S1e,S1f, S1g.grid)
tiff("03282022_Supplementary figure 1.tiff", width=4800, height=4800, res=300)
grid.arrange(grobs = grobsS1, layout_matrix = layS1)
dev.off()






#-----Fig S2
layS2.1 <- rbind(c(1,1,1,1,1,1),
                 c(1,1,1,1,1,1),
                 c(1,1,1,1,1,1),
                 c(1,1,1,1,1,1))

tiff("03312022_Supplementary figure 2-1.tiff", width=3000, height=1500, res=300)
grid.arrange(grobs = list(grid.grabExpr(draw(S1e))), layout_matrix = layS2.1)
dev.off()

layS2.2 <- rbind(
                 c(1,1,1),
                 c(1,1,1),
                 c(1,1,1),
                 c(2,2,3),
                 c(2,2,3),
                 c(2,2,3))
tiff("03242022_Supplementary figure 2-2.tiff", width=3000, height=3000, res=300)
grid.arrange(grobs = list(S1f,S1g.grid,grob(NULL)), layout_matrix = layS2.2)
dev.off()




#-----Fig S3
#heatmap of all nominal and padj DEGs
nomDEGs <- rownames(ordered.res[abs(ordered.res$log2FoldChange) > 0.5 & ordered.res$pvalue < 0.05,])
nomDEGs.dat <- batch.rem[rownames(batch.rem) %in% nomDEGs,]

nomDEGs.dat <- nomDEGs.dat[,order(ncol(nomDEGs.dat):1)]
identical(colnames(nomDEGs.dat), rownames(rev.covs))

nom.ha = rowAnnotation(foo = anno_mark(labels_gp = gpar(fontsize = 7), at = which(rownames(nomDEGs.dat) %in% prev.degs), 
                                   labels = rownames(nomDEGs.dat)[rownames(nomDEGs.dat) %in% prev.degs]))
#be careful when setting labels: important to have the order as the genes occur in the dataframe/matrix
S2a<-ComplexHeatmap::Heatmap(t(scale(t(nomDEGs.dat))),
                        cluster_rows = T,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = F,
                        col = colors,
                        show_column_names = F,
                        column_names_rot = 90, 
                        top_annotation = ha,
                        right_annotation = nom.ha,
                        name = "row Z score",
                        clustering_method_rows = "ward.D2",
                        clustering_distance_rows = "pearson"
)

S2a <- draw(S2a,
            annotation_legend_side = "bottom")
#Corrected pca plot
identical(colnames(batch.rem),rownames(annot_df.rem))
pca.cor <- plotPCA.custom(as.matrix(batch.rem), intgroup=c("Status"), 
                          ntop = 50000, returnData=TRUE, metadata=annot_df.rem, pc.1 = 1, pc.2 = 2)
PoV.cor <- round(100 * attr(pca.cor, "percentVar"))
S2b<-PCAplot(pca.cor, "Status", PoV.df=PoV.cor, pc.1 = 1, pc.2 = 2, colors = c('lightblue','tomato')) #samples 36 and 37 seem to be outliers
test.gg <- S2b+labs(title = 'PCA on residual expression')+theme(plot.title = element_text(size = 12, face = "bold"))
S2b <- alignTitles(test.gg, title = 1)

S2d <- ggplot(EXPR,aes_string(x='Status',y='Expression', fill='Status')) +
  geom_boxplot(fill=c('lightblue','tomato')) +
  geom_line(aes(group=paired), position = position_dodge(0.2)) +
  geom_point(aes(fill=Status,group=paired),size=2,shape=21, position = position_dodge(0.2)) +
  theme_classic()+
  scale_fill_manual(values = c('lightblue','tomato'))

layS2 <- rbind(c(1,1),
               c(1,1),
               c(1,1),
               c(1,1)
               )
S2a.grid <- grid.grabExpr(draw(S2a))
textgrob2 <- grobTree(
  gp = gpar(fontsize = 13), 
  textGrob(label = "All genes at nominal p < 0.01", name = "title1",
           x = unit(0.5, "lines"), y = unit(0, "lines"), 
           just = "left", hjust = 0, vjust = 0,
           gp = gpar(col = "black", fontface = "bold"))
)

S2a.grid<-arrangeGrob(S2a.grid, top = textgrob2, padding = unit(1, "line"))

S2c.grid <- grid.grabExpr(draw(S2c,
                               annotation_legend_side = "bottom", 
                          heatmap_legend_side="bottom", legend_grouping = 'adjusted'))


#list(S1a.grid,S1b,S1c.grid,S1d,S1e,S1f, S1g.grid)
tiff("03312022_Supplementary figure 3.tiff", width=1350, height=2000, res=300)
grid.arrange(grobs = list(S2a.grid), layout_matrix = layS2)
dev.off()

#-----Fig S4
layS4 <- rbind(c(4,4),
               c(4,4))
tiff("03292022_Supplementary figure 4.tiff", width=1200, height=1000, res=300)
grid.arrange(grobs = list(S2b), layout_matrix = layS4)
dev.off()

#-----Fig S5
#NFkB
dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"NFkB signaling"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
p.S3a <- PCAplot(list.pca, "Status", PoV.df=list.pov, pc.1=1, pc.2=2, colors=c("lightblue", "tomato"), geom.size=4)
S3a.title <- p.S3a+labs(title = 'PCA on residual expression of NFkB genes')+theme(plot.title = element_text(size = 12, face = "bold"))
p.S3a <- alignTitles(S3a.title, title = 1)
scatter.df$`NFkB signaling` <- list.pca$PC1

dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"LPS response"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
p.S3b <- PCAplot(list.pca, "Status", PoV.df=list.pov, pc.1=1, pc.2=2, colors=c("lightblue", "tomato"), geom.size=4)
S3b.title <- p.S3b+labs(title = 'PCA on residual expression of LPS response genes')+theme(plot.title = element_text(size = 12, face = "bold"))
p.S3b <- alignTitles(S3b.title, title = 1)
scatter.df$`LPS response` <- list.pca$PC2

dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"Glucocorticoid response"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
p.S3c <- PCAplot(list.pca, "Status", PoV.df=list.pov, pc.1=1, pc.2=2, colors=c("lightblue", "tomato"), geom.size=4)
S3c.title <- p.S3c+labs(title = 'PCA on residual expression of GC response genes')+theme(plot.title = element_text(size = 12, face = "bold"))
p.S3c <- alignTitles(S3c.title, title = 1)
scatter.df$`Glucocorticoid response` <- list.pca$PC2



scatter.df[,c('NFkB signaling','Glucocorticoid response','LPS response')]
cbind(scatter.df$Status,clean_covs$antipsychotics)
p.S3d <- MEs_lin(as.matrix(scatter.df[,c('NFkB signaling','Glucocorticoid response','LPS response')]),
        cbind('Status'=scatter.df$Status,'Antipsychotics'=clean_covs$antipsychotics))

layS3 <- rbind(c(1,1,4,4),
               c(1,1,6,6),
               c(2,2,6,6),
               c(2,2,5,5),
               c(3,3,5,5),
               c(3,3,5,5))
S3d.grid <- grid.grabExpr(draw(p.S3d))
textgrob3 <- grobTree(
  gp = gpar(fontsize = 13), 
  textGrob(label = "PC-Status correlation by pathway", name = "title1",
           x = unit(0.5, "lines"), y = unit(0, "lines"), 
           just = "left", hjust = 0, vjust = 0,
           gp = gpar(col = "black", fontface = "bold"))
)

S3d.grid<-arrangeGrob(S3d.grid, top = textgrob3, padding = unit(1, "line"))


grobsS3 <- list(p.S3a,p.S3b,p.S3c,grob(NULL),grob(NULL),S3d.grid)

#list(S1a.grid,S1b,S1c.grid,S1d,S1e,S1f, S1g.grid)
tiff("03292022_Supplementary figure 5.tiff", width=3100, height=3700, res=300)
grid.arrange(grobs = grobsS3, layout_matrix = layS3)
dev.off()



cor(scatter.df$`Glucocorticoid response`,scatter.df$`NFkB signaling`)
cor(scatter.df$`LPS response`,scatter.df$`NFkB signaling`)



#========supplementary tables=========#
setwd(file.path("/Users/kubler01/Documents/R projects/Monocyte project/Covariates"))

#Supplementary tables
write.xlsx2(NA, file='Supplementary tables.xlsx', sheetName='Overview')

#Patient infos
#Read in covariates and sample information; contains all correct information
raw.covs <- read.xlsx2("covariates_SCZMono_corrected_with21.xlsx", sheetIndex=1)
#Includes all correct covariates now

library(mice)
rownames(raw.covs) <- raw.covs$sample_ID
raw.covs <- raw.covs[,-c(1,2,4,6)] #remove sample, subject, and cell-type columns

pat.info <- cbind(complete(imp.2), raw.covs[,colnames(raw.covs) %in% c("RIN", "Genes.Detected", "Mean.Per.Basce.Cov.", "Mapped", "Mean.Per.Base.Cov.")])
rownames(pat.info) <- gsub(".*_","",rownames(pat.info))
pat.info <- pat.info[order(as.numeric(rownames(pat.info))),]

write.xlsx2(pat.info, file='Supplementary tables.xlsx', sheetName='Patient information', append=T)

new.genemap
DEG.res <- res[abs(res$log2FoldChange) > 0 & res$padj < 0.1,]
DEG.IDs <- new.genemap[new.genemap$external_gene_name %in% rownames(DEG.res),]
identical(DEG.IDs$external_gene_name,rownames(DEG.res))
DEG.res$GeneID <- DEG.IDs$ensembl_gene_id

write.xlsx2(DEG.res, file='Supplementary tables.xlsx', sheetName='DEGs FDR < 0.1', append=T)

gene.panels_mrgd

for (i in colnames(gene.panels_mrgd))
  gene.panels_mrgd[,i] <- sort(as.character(gene.panels_mrgd[,i]), na.last = TRUE)

write.xlsx2(gene.panels_mrgd, file='Supplementary tables.xlsx', sheetName='Gene panels', append=T)







#=======END=======#

pc.box <- list(grid.grabExpr(plot(boxplots$LPS)),grid.grabExpr(plot(boxplots$NFkB)),grid.grabExpr(plot(boxplots$GC))) 
p6 <- ggarrange(plotlist=pc.box[1:length(pc.box)], nrow = 1, common.legend = F)

p3.grid <- grid.grabExpr(draw(p3))
p5.grid <- grid.grabExpr(draw(p5))
p6.grid <- grid.draw(grid.grabExpr(draw(p6)))

lay <- rbind(c(1,1,2,3),
             c(1,1,2,3),
             c(4,4,2,3),
             c(4,4,5,5),
             c(6,6,6,6))

grobs <- list(p1,p2,p3.grid,p4,p5.grid,p6)

tiff("03112022_Figure 1.tiff", width=5000, height=4000, res=300, units='mm')
grid.arrange(grobs = grobs, layout_matrix = lay)
dev.off()

alignTitles <- function(ggplot, title = 2, subtitle = 2, caption = 2) {
  # grab the saved ggplot2 object
  g <- ggplotGrob(ggplot)
  
  
  # find the object which provides the plot information for the title, subtitle, and caption
  g$layout[which(g$layout$name == "title"),]$l <- title
  g$layout[which(g$layout$name == "subtitle"),]$l <- subtitle
  g$layout[which(g$layout$name == "caption"),]$l <- caption
  g
}




#Just checking pathway genes on heatmap
test1 <- res[rownames(res) %in% gene.panels_mrgd$`Glucocorticoid response`,]
test1 <- test1[abs(test1$padj) < 0.1,]
test1 <- test1[order(test1$log2FoldChange, method="radix", decreasing=T),]

df <- batch.rem[,order(as.numeric(colnames(batch.rem)), method="radix", decreasing=T)]

ComplexHeatmap::Heatmap(t(scale(t(df[rownames(df) %in% c('IL5','IL6','IL1B'),]))),
                        cluster_rows = T,
                        cluster_columns =F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = T,
                        column_names_rot = 60, 
                        #right_annotation  = ra,
                        top_annotation = ha,
                        name = "Normalized counts",
                        clustering_method_rows = "ward.D2",
                        clustering_distance_rows = "pearson", 
                        row_names_gp = grid::gpar(fontsize = 6),
                        column_names_gp = grid::gpar(fontsize = 5))

dev.off()


cor(t(batch.rem[rownames(batch.rem) %in% c('IL6','FES'),]))


test1 <- res[rownames(res) %in% gene.panels_mrgd$`NFkB signaling`,]
test1 <- test1[abs(test1$padj) < 0.1 && test1$log2FoldChange > 0,]

