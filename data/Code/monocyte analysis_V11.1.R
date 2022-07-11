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


#==========================================================================================================#
#                                                  new functions                                           #
#==========================================================================================================#

MEs_lin <- function(MEs,cov.df, return.p=F, return.r=F){
  l <- length(colnames(MEs))
  covariates <- colnames(cov.df)
  matrix_rsquared = matrix(NA, nrow = length(covariates), ncol = l) #Number of factors
  matrix_pvalue = matrix(NA, nrow = length(covariates), ncol = l)
  
  
  for (x in 1:length(covariates)){
    for (y in 1:l){
      matrix_rsquared[x,y] <- summary( lm(MEs[,y] ~ cov.df[,covariates[x]]) )$adj.r.squared
      matrix_pvalue[x,y] <- glance( lm(MEs[,y] ~ cov.df[,covariates[x]]) )$p.value #To insert pvalues in the heatmap
    }
  }
  
  rownames(matrix_rsquared) = covariates
  rownames(matrix_pvalue) = covariates 
  colnames(matrix_rsquared) = colnames(MEs)
  
  matrix_pvalue = matrix(p.adjust(as.vector(as.matrix(matrix_pvalue)), method='bonferroni'),ncol=ncol(matrix_pvalue))
  matrix_pvalue = formatC(matrix_pvalue, format = "e", digits = 2)
  
  matrix_rsquared[is.na(matrix_rsquared)] <- 0
  # png(paste0(work_dir, "LinearReg_15pcs_covariates.png"), width = 10, height = 10, res = 600, units = "in")
  # pdf(paste0(work_dir, "LinearReg_15pcs_covariates.pdf"), width = 10, height = 10)
  p<-Heatmap(matrix_rsquared, 
          cluster_rows = T,
          cluster_columns = F,
          show_row_dend = T,
          show_column_dend = F,
          show_row_names = T,
          show_column_names = T,
          name = 'adjusted R squared',
          clustering_method_rows = "ward.D2",
          clustering_distance_rows = "euclidean",
          col = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30), 
          column_names_gp = grid::gpar(fontsize = 8),
          row_names_gp = grid::gpar(fontsize = 8),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(matrix_pvalue[i,j], x, y, gp = gpar(fontsize = 5))}
  )
  return(p)
  if (return.p)
    return(matrix_pvalue)
  if (return.r)
    return(matrix_rsquared)
}

plotPCA.custom <- function (object, intgroup = "condition", ntop = 44099, returnData = TRUE, metadata = metadata, pc.1=1, pc.2=2) 
{
  rv <- rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(object[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% colnames(metadata))) {
    stop("the argument 'intgroup' should specify columns of metadata")
  }
  intgroup.df <- data.frame(metadata[intgroup],check.names=F)
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    metadata[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, pc.1], PC2 = pca$x[, pc.2], group = group, 
                  intgroup.df, name = colnames(object), check.names=F)
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(pc.1,pc.2)]
    return(d)
  }
  ggplot(data = d, aes_string(x = paste("PC",pc.1), y = paste("PC",pc.2), 
                              color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0(paste("PC", pc.1, ": ", sep=""),round(percentVar[1] * 100), "% variance")) +
    ylab(paste0(paste("PC", pc2, ": ", sep=""),round(percentVar[2] * 100), "% variance")) + 
    coord_fixed()
}

PCA.covar <- function(df, meta, heatmap.r=F, data.r=F, heatmap.p=T,type="spearman", heatmap.color.code=NA){
  l <- length(colnames(meta))
  l2 <- l+1
  pca.cov <- prcomp(t(df))
  PC.covs <- pca.cov$x[,1:l]
  PC_cov_correlation <- rcorr(as.matrix(PC.covs), as.matrix(meta), type=type)
  PC_variance.explained <- PC_cov_correlation$r[1:l,l2:length(rownames(PC_cov_correlation$r))]
  PC_cov.cor.p <- PC_cov_correlation$P[1:l,l2:length(rownames(PC_cov_correlation$r))]
  if (heatmap.r)
    pheatmap(PC_variance.explained, cluster_rows = F, cluster_cols = F, show_colnames = T,
             color = heatmap.color.code)
  if (heatmap.p)
    pheatmap(PC_cov.cor.p, cluster_rows = F, cluster_cols = F, show_colnames = T,
             color = colorRampPalette(c("red", "white"))(50),
             display_numbers = T)
  if (data.r)
    return (PC.covs)
}
createDT <- function(DF, caption="", scrollY=500){
  data <- DT::datatable(DF, caption=caption,
                        extensions = 'Buttons',
                        options = list( dom = 'Bfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        scrollY = scrollY, scrollX=T, scrollCollapse = T, paging = F,
                                        columnDefs = list(list(className = 'dt-center', targets = "_all"))
                        )
  )
  return(data)
}

PCAplot <- function(pca.df, Condition, Shape=NULL, pc.1=1, pc.2=2, PoV.df, colors=c("#00AFBB", "#E7B800", "#FC4E07"), geom.size=2){
  if(!is.null(Shape))
    ggplot(pca.df, aes(PC1, PC2, color=pca.df[[c(Condition)]], label=name)) +
    geom_point(size=geom.size, aes(shape=pca.df[[Shape]]))+
    labs(x = paste0("PC", paste(pc.1), ": ", PoV.df[1], "% variance"), y = paste0("PC", paste(pc.2), ": ",PoV.df[2],"% variance"), 
         color = Condition, shape=Shape) +  
    scale_color_manual(values = colors)+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5))
  else
    ggplot(pca.df, aes(PC1, PC2, color=pca.df[[c(Condition)]], label=name)) +
    geom_point(size=2)+
    labs(x = paste0("PC", paste(pc.1), ": ", PoV.df[1], "% variance"), y = paste0("PC", paste(pc.2), ": ",PoV.df[2],"% variance"), 
         color = Condition) +  
    scale_color_manual(values = colors)+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5))
}

PCA_cov_cor_R <- function(cov.df,df, fs1=8,fs2=5){
  l <- length(colnames(df))
  covariates <- colnames(cov.df)
  pca.cov <- prcomp(df)
  var <- get_pca_var(pca.cov) # description of PCs
  ind <- get_pca_ind(pca.cov) # PCs for individuals
  matrix_rsquared = matrix(NA, nrow = length(cov.df), ncol = l) #Number of factors
  matrix_pvalue = matrix(NA, nrow = length(cov.df), ncol = l)
  
  for (x in 1:length(covariates)){
    for (y in 1:l){
      matrix_rsquared[x,y] <- summary( lm(var$cor[,y] ~ cov.df[,covariates[x]]) )$adj.r.squared
      matrix_pvalue[x,y] <- glance( lm(var$cor[,y] ~ cov.df[,covariates[x]]) )$p.value #To insert pvalues in the heatmap
    }
  }
  
  rownames(matrix_rsquared) = covariates
  rownames(matrix_pvalue) = covariates 
  colnames(matrix_rsquared) = paste("PC",1:l, sep="")
  
  matrix_pvalue = matrix(p.adjust(as.vector(as.matrix(matrix_pvalue)), method='bonferroni'),ncol=ncol(matrix_pvalue))
  matrix_pvalue = formatC(matrix_pvalue, format = "e", digits = 2)
  
  # png(paste0(work_dir, "LinearReg_15pcs_covariates.png"), width = 10, height = 10, res = 600, units = "in")
  # pdf(paste0(work_dir, "LinearReg_15pcs_covariates.pdf"), width = 10, height = 10)
  Heatmap(matrix_rsquared, 
          cluster_rows = T,
          cluster_columns = F,
          show_row_dend = T,
          show_column_dend = F,
          show_row_names = T,
          show_column_names = T,
          name = 'R squared',
          clustering_method_rows = "ward.D2",
          clustering_distance_rows = "euclidean",
          col = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30), 
          column_names_gp = grid::gpar(fontsize = fs1),
          row_names_gp = grid::gpar(fontsize = fs1),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(matrix_pvalue[i,j], x, y, gp = gpar(fontsize = fs2))}
  )
}

setEnrichment <- function(set1, set2, universe = 18000){
  a = sum(set1 %in% set2)
  c = length(set1) - a
  
  b = length(set2) - a
  d = universe - length(set2) - c
  
  contingency_table = matrix(c(a, c, b, d), nrow = 2)
  # one-tailed test for enrichment only
  fisher_results = fisher.test(contingency_table, alternative = "greater")
  # returns data.frame containing the lengths of the two sets, the overlap, the enrichment ratio (odds ratio) and P value
  df <- tibble::tibble( set1_length = length(set1), set2_length = length(set2), overlap = a, ratio = fisher_results$estimate, p.value = fisher_results$p.value)
  return(df)
}

GSEA.byMod <- function(mod.gl, gl.of.int., universe=18000){ #mod.gl = list of genes per module, gl.of.int = gene list to test enrichment for
  df <- NULL
  for (i in names(mod.gl)){
    df <- data.frame(rbind(df, setEnrichment(mod.gl[[i]], gl.of.int.)))}
  df$q.val <- p.adjust(as.vector(as.matrix(df$p.value)), method='fdr')
  rownames(df) <- names(mod.gl)
  return(df)
}

heatmap.color.code <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(496)

plot.density <- function(meta.df,to.color="Condition",to.plot="Library size",ylab.title="Density",xlab.title="Library size"){
  meta.df %>% 
    ggplot(aes_string(color=to.color, x=to.plot, fill=to.color)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab(ylab.title) +
    xlab(xlab.title)
}


#==========================================================================================================#
#                                   RNA-seq analysis start: counts files                                   #
#==========================================================================================================#
dir <- file.path("/Users/kubler01/Documents/R projects/Monocyte project/Counts")
setwd(dir)
raw.counts <- read_excel("Supporting data monocyte paper.xlsx", sheet=2)
raw.counts <- data.frame(raw.counts, check.names=F)
rownames(raw.counts) <- raw.counts[,1]
raw.counts <- raw.counts[,-1]
raw.counts <- raw.counts[,order(as.numeric(colnames(raw.counts)))]
#Plan: take IFN-y signature/other signatures, calculate first PC, plot over samples ordered by SCZ-duration
#Then perform analyses as done in Monocyte paper



#==========================================================================================================#
#                                      RNA-seq analysis start: covariate file                              #
#==========================================================================================================#
#Setting the directory to get the metadata
dir <- file.path("/Users/kubler01/Documents/R projects/Monocyte project/Covariates")
setwd(dir)
#Read in covariates and sample information; contains all correct information
new.covs <- read.xlsx2("covariates_SCZMono_corrected_with21.xlsx", sheetIndex=1)
#Includes all correct covariates now

#Imputing missing cov ariates
library(mice)
rownames(new.covs) <- new.covs$sample_ID
new.covs <- new.covs[,-c(1,2,4,6)] #remove sample, subject, and cell-type columns

for (i in colnames(new.covs))
  new.covs[[i]] <- as.numeric(new.covs[[i]])
new.covs$Mapped <- scale(new.covs$Mapped)
new.covs$Genes.Detected <- scale(new.covs$Genes.Detected)
new.covs$Mean.Per.Base.Cov. <- scale(new.covs$Mean.Per.Base.Cov.)
new.covs$RIN <- scale(new.covs$RIN)

for (i in c("status", "smoking", "coffee", "sex", "batch", "drugs", "antipsychotics", "antidepressive", "antianxiety"))
  new.covs[,i] <- as.factor(new.covs[,i])

# Deterministic regression imputation via mice
imp <- mice(new.covs[,!colnames(new.covs) %in% c("duration_illness","age", "RIN", "Genes.Detected", "Mean.Per.Basce.Cov.", "Mapped", "Mean.Per.Base.Cov.", "BMI")], method = "logreg", m = 1)
# Store data
data_imp <- complete(imp)
imp.2 <- mice(cbind(data_imp,"BMI"=new.covs[,"BMI"]), method = "norm.predict", m = 1)

imputed.covs <- cbind(complete(imp.2), new.covs[,colnames(new.covs) %in% c("RIN", "Genes.Detected", "Mean.Per.Basce.Cov.", "Mapped", "Mean.Per.Base.Cov.")])
rownames(imputed.covs) <- gsub(".*_","",rownames(imputed.covs))
imputed.covs <- imputed.covs[order(as.numeric(rownames(imputed.covs))),]


#Make sure samples in covariate file are the same as in count matrix
identical(rownames(imputed.covs), colnames(raw.counts))
#raw.counts <- raw.counts[colnames(raw.counts) %in% rownames(imputed.covs)] 


#==========================================================================================================#
#                                         RNA-seq analysis start: gene panels                              #
#==========================================================================================================#
#Loading gene panels
setwd('/Users/kubler01/Documents/R projects/Monocyte project/')
gene.panels <- read.xlsx2("10142021_monocyte paper_gene lists pathways.xlsx", sheetIndex=1)
gene.panels <- gene.panels[,-2] #Already using NFKB gene signature

GC.acute <- read.xlsx2("GC acute signature.xlsx", sheetIndex=1)
GC.acute.up <- GC.acute$Upreg
GC.acute.down <- GC.acute$Downreg
GC.chronic <- read.xlsx2("GC chronic signature.xlsx", sheetIndex=1)
GC.chronic <- GC.chronic$Official.Gene.Symbol
IL4 <- read.xlsx2("IL4 signature.xls", sheetIndex=1)
IL4 <- toupper(IL4$X.)
LPS <- read.xlsx2("LPS-response genes.xlsx", sheetIndex=1)
LPS <- LPS[,1][-c(1,2)]

DEX.panel <- read_excel("DEX DEGs.xls", sheet=1) #GC stimulation with DEX in healthy controls
#Fix DEX table
DEX.panel <- data.frame(DEX.panel[c(-1:-5),c(-2,-5:-6)])
colnames(DEX.panel) <- DEX.panel[1,]
DEX.panel <- DEX.panel[-1,]
#rownames(DEX.panel) <- DEX.panel$Symbol
DEX.panel <- DEX.panel[DEX.panel$FC > 1,]


GC.panel <- read_excel("GC stimulation in blood cells.xlsx", sheet=3) #PTSD PBMCs stimulated with DEX


#Fix tables
GC.panel2.5nm <- data.frame(GC.panel[c(-1),c(1,2,3)], check.names = F, check.rows = F)
colnames(GC.panel2.5nm) <- GC.panel2.5nm[1,]
rownames(GC.panel2.5nm) <- GC.panel2.5nm[,1]
GC.panel2.5nm <- GC.panel2.5nm[-1,-1]

#Make numeric
GC.panel2.5nm <- data.frame(apply(GC.panel2.5nm, 2, function(x) as.numeric(as.character(x))), row.names=rownames(GC.panel2.5nm))
#GC.panel2.5nm$logFC <- formatC(GC.panel2.5nm$logFC, digits = 2)
#GC.panel2.5nm$adj.P.Val <- formatC(GC.panel2.5nm$adj.P.Val, format="e", digits = 2)

GC.panel5nm <- data.frame(GC.panel[c(-1),c(1,6,7)], check.names = F, check.rows = F)
colnames(GC.panel5nm) <- GC.panel5nm[1,]
rownames(GC.panel5nm) <- GC.panel5nm[,1]
GC.panel5nm <- GC.panel5nm[-1,-1]
#Make numeric
GC.panel5nm <- data.frame(apply(GC.panel5nm, 2, function(x) as.numeric(as.character(x))), row.names=rownames(GC.panel5nm))

GC.panel2.5nm <- GC.panel2.5nm[GC.panel2.5nm$adj.P.Val < 0.1 & GC.panel2.5nm$logFC > 1 | GC.panel2.5nm$logFC < -1,]
GC.panel5nm <- GC.panel5nm[GC.panel5nm$adj.P.Val < 0.1 & GC.panel5nm$logFC > 1 | GC.panel5nm$logFC < -1,]



GC.response <- unique(c(na.omit(DEX.panel[,1]), GC.chronic, c(GC.acute.up,GC.acute.down), rownames(GC.panel2.5nm), rownames(GC.panel5nm)))
gene.panels_2 <- rbind.fill(data.frame(LPS),data.frame(IL4), data.frame(na.omit(GC.response)))
colnames(gene.panels_2) <- c("LPS response", "IL4 response", "Glucocorticoid response")



#Gandal markers
setwd("/Users/kubler01/Documents/R projects")
Gandal.clusters <- read_excel("Gandal (2018) - Gene and isoform co-expression module annotation.xlsx", sheet=2)
#Restrict to expressed genes only
Gandal.clusters.annot <- Gandal.clusters[colnames(Gandal.clusters) %in% c("ensembl_gene_id", "Module", "gene_name")]
Gandal.clusters.annot <- Gandal.clusters.annot[Gandal.clusters.annot$ensembl_gene_id %in% rownames(raw.counts),]

Gandal.clusters.annot[Gandal.clusters.annot$Module %in% "geneM6",]$Module <- "Microglia"
Gandal.clusters.annot[Gandal.clusters.annot$Module %in% "geneM3",]$Module <- "Astrocytes"
Gandal.clusters.annot[Gandal.clusters.annot$Module %in% "geneM5",]$Module <- "NFkB"
#Gandal.clusters.annot[Gandal.clusters.annot$Module %in% "geneM32",]$Module <- "IFNy"

#Patir microglia
setwd("/Users/kubler01/Documents/R projects/Monocyte project/Gene lists")
Patir <- read_excel("PatirEtAl_coreMicrogliaGenes.xlsx", col_names = "Patir_microglia")

#Manually chosen genes
man.genes <- c("IL1B", "IL6", "TNF", "IL10", "P2RY12", "CXC3R1", "CSF1R", "TMEM119", "AIF1", "CD68", "C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "CSF1R", 
               "P2RX7", "PROS1", "PTPRC", "RUNX1", "SPI1", "TGFBR1", "TLR4", "TREM2", "CCR2")

#SCZ TWAS & GWAS genes
SCZ_TWAS <- read_excel("TWAS_SCZ.xlsx")
SCZ.genes <- SCZ_TWAS$gene_name

SCZ.GWAS <- data.frame(read_excel("Pardinas SCZ GWAS supplementary tables.xlsx", sheet = 6))[-c(1:3),]
rownames(SCZ.GWAS) <- SCZ.GWAS[4,]

rownames(SCZ.GWAS) <- SCZ.GWAS[1,]
SCZ.Ripke <-  read_excel("Ripke SCZ GWAS prioritised genes.xlsx", sheet = 2)
SCZ.Ripke.symbols <- SCZ.Ripke$Symbol.ID

SCZ.GWAS.genes <- SCZ.GWAS$...9[-1]

SCZ.GWAS.genes_split <- NULL
for (i in SCZ.GWAS.genes)
    SCZ.GWAS.genes_split <- c(SCZ.GWAS.genes_split,strsplit(i, ", ")[[1]])

SCZ.risk <- unique(c(SCZ.genes, SCZ.GWAS.genes_split, SCZ.Ripke.symbols))



#==========================================================================================================#
#                                         RNA-seq analysis start: prelim QC                                #
#==========================================================================================================#
#Prelim QC:

#Library size; no need to filter: normally distributed
plot.density(imputed.covs,to.color="status",to.plot="`Mapped`",xlab.title="Reads mapped")
plot.density(imputed.covs,to.color="status",to.plot="`Exonic.Rate`",xlab.title="Exonic rate")
plot.density(imputed.covs,to.color="status",to.plot="`Genes.Detected`",xlab.title="Genes detected")
plot.density(imputed.covs,to.color="status",to.plot="`RIN`",xlab.title="RIN")
plot.density(imputed.covs,to.color="status",to.plot="`RIN`",xlab.title="RIN")

#Filter out lowly expressed genes
filt.df <- raw.counts[rowMeans(raw.counts) > 1,] 

#Check inter-sample distances to see which one should be thrown out
annot_df <- data.frame("Status"=as.character(imputed.covs$status), "Sex"=as.character(imputed.covs$sex), check.names = F)
annot_df$Status[annot_df$Status == 1] <- "SCZ"
annot_df$Status[annot_df$Status == 0] <- "CTR"
annot_df$Sex[annot_df$Sex == 1] <- "Male"
annot_df$Sex[annot_df$Sex == 2] <- "Female"

rownames(annot_df) <- rownames(imputed.covs)
cols <- c("lightblue","tomato")
cols.2 <- c("grey","lightpink")

annot_colors <- list("Status"=cols, "Sex"=cols.2) #annotation colors for heatmap
names(annot_colors$`Status`) <- c("CTR","SCZ")
names(annot_colors$`Sex`) <- c("Male","Female")


ComplexHeatmap::pheatmap(cor(filt.df,filt.df), annotation_row = annot_df, annotation_col = annot_df,
                         annotation_colors = annot_colors, color = heatmap.color.code)


#Check influence of covariates on data variance
PCA_cov_cor_R(imputed.covs, filt.df)

#Plot variancePartitioning results
setwd("/Users/kubler01/Documents/R projects/R workspace/")
load(file="12132021_VarPar_Monocytes.RData")
plotVarPart(vp)

#Based on variance partitioning hence the following covariates should be corected for: Mean per base coverage | mapped number of genes |
#exonic rate | number of genes detected | rRNA rate | RIN | BMI | age_clusters | batch
#If we go by only PCs prior to status effect: Mean.Per.Base.Coverage + Mapped + Exonic.Rate + Genes.Detected + rRNA.rate + RIN + BMI + batch

#Check covariate co-linearity
for (i in colnames(new.covs))
  new.covs[[i]] <- as.numeric(new.covs[[i]])
MEs_lin(new.covs,new.covs) #none of the biological covariates are co-linear with technical covariates

#Include variance partitioning




#==========================================================================================================#
#                                     RNA-seq analysis: gene annotation                                    #
#==========================================================================================================#
#Annotate the genes
library(biomaRt)
geneIDs <- rownames(filt.df)
listMarts(host="www.ensembl.org")
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
dfs.ensembl <- listDatasets(ensembl)
genemap.mic <- getBM(values = geneIDs,
                     filters = "ensembl_gene_id",
                     mart = ensembl,
                     attributes = c("ensembl_gene_id", "entrezgene_id",
                                    "hgnc_symbol", "external_gene_name", "gene_biotype",
                                    "description", "chromosome_name",
                                    "strand"))


#Remove non-annotated genes
new.genemap <- genemap.mic
new.genemap <- new.genemap[!new.genemap$external_gene_name=="",]

#Remove everything except protein-coding genes
new.genemap <- new.genemap[new.genemap$gene_biotype=="protein_coding",]

new.genemap$external_gene_name <- make.unique(new.genemap$external_gene_name)
new.genemap$ensembl_gene_id <- make.unique(new.genemap$ensembl_gene_id)

#Sort count df and weed out non-annotated genes
filt.df <- filt.df[rownames(filt.df) %in% new.genemap$ensembl_gene_id,]
new.genemap <- new.genemap[new.genemap$ensembl_gene_id %in% rownames(filt.df),]

filt.df <- filt.df[order(rownames(filt.df), decreasing = F),]
new.genemap <- new.genemap[order(new.genemap$ensembl_gene_id, decreasing = F),]
identical(rownames(filt.df),new.genemap$ensembl_gene_id)

rownames(filt.df) <- new.genemap$external_gene_name

#Check again
ComplexHeatmap::pheatmap(cor(filt.df,filt.df), annotation_row = annot_df, annotation_col = annot_df,
                         annotation_colors = annot_colors, color = heatmap.color.code)



#==========================================================================================================#
#                                      RNA-seq analysis start: DESeq2 start                                #
#==========================================================================================================#
#Transform data to reduce inter-sample distances
#Create a DESeq2 matrix
imputed.covs$BMI <- scale(imputed.covs$BMI)
imputed.covs$age_clusters <- as.factor(imputed.covs$age_clusters)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(filt.df),
                              colData = imputed.covs,
                              design = ~ Mean.Per.Base.Cov. + Mapped + Exonic.Rate + Genes.Detected + rRNA.rate + RIN + BMI + age_clusters + batch + status)


#version 2 based on dds object
# Estimate library size correction scaling factors
dds <- estimateSizeFactors(dds)

vsd <- vst(dds)
trnsf.df <- data.frame(assay(vsd), check.names=F)

#Check inter-sample distances subsequent to data transformation
ComplexHeatmap::pheatmap(cor(trnsf.df,trnsf.df), annotation_row = annot_df, annotation_col = annot_df,
                         annotation_colors = annot_colors, color = heatmap.color.code)
#Indeed now samples 36 and 37 seem to be outliers, might consider removing them


#Check influence of covariates on data variance after transformation
PCA_cov_cor_R(imputed.covs, trnsf.df)

pca <- plotPCA.custom(as.matrix(trnsf.df), intgroup=c("status", "batch", "sex", "smoking"), 
                      ntop = 50000, returnData=TRUE, metadata=imputed.covs, pc.1 = 1, pc.2 = 2)
PoV <- round(100 * attr(pca, "percentVar"))
PCAplot(pca, "batch", PoV.df=PoV, colors=c('orange','lightblue','lightgreen',"yellow"), pc.1 = 1, pc.2 = 2)
PCAplot(pca, "status", PoV.df=PoV, pc.1 = 1, pc.2 = 2) #samples 36 and 37 seem to be outliers

#Removing samples 36 and 37
clean_df <- trnsf.df[,!colnames(trnsf.df) %in% c("36", "37")]
clean_rawdf <- filt.df[,!colnames(filt.df) %in% c("36", "37")]
clean_covs <- imputed.covs[!rownames(imputed.covs) %in% c("36", "37"),]

#Checking effect on data
pca.clean <- plotPCA.custom(as.matrix(clean_df), intgroup=c("status", "batch", "sex", "smoking"), 
                      ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1 = 1, pc.2 = 2)
PoV.clean <- round(100 * attr(pca.clean, "percentVar"))
PCAplot(pca.clean, "status", PoV.df=PoV.clean, pc.1 = 1, pc.2 = 2)

PCA_cov_cor_R(clean_covs, clean_df) #status now stronger, less technical influence!



#Change main DESeq2 data and covariate frames, adding gender as variable to design
dds <- DESeqDataSetFromMatrix(countData = as.matrix(clean_rawdf),
                              colData = clean_covs,
                              design = ~ Mean.Per.Base.Cov. + Mapped + Exonic.Rate + Genes.Detected + rRNA.rate + RIN + BMI + age_clusters + batch + sex + status)


#Save files for cluster analysis
save(clean_rawdf, clean_df,clean_covs,filt.df, trnsf.df, imputed.covs, file="12282021_monocyte objects for DESeq2.RData")


dir
#==========================================================================================================#
#                                      RNA-seq analysis: data correction                                   #
#==========================================================================================================#
#Correct data with batch.rem()
batch.rem <- removeBatchEffect(clean_df, batch = clean_covs$batch, 
                               covariates=as.matrix(cbind(clean_covs$Mean.Per.Base.Cov.,
                                                          clean_covs$Mapped,
                                                          clean_covs$Exonic.Rate,
                                                          clean_covs$Genes.Detected,
                                                          clean_covs$rRNA.rate,
                                                          clean_covs$RIN,
                                                          clean_covs$BMI,
                                                          clean_covs$sex,
                                                          clean_covs$age_clusters)),
                               design=model.matrix(~ clean_covs$status)) #this is the dds.2 design and based on the varPar and indeed it gives cleaner results
PCA_cov_cor_R(clean_covs, batch.rem)

pca.cor <- plotPCA.custom(as.matrix(batch.rem), intgroup=c("status", "batch", "sex", "smoking"), 
                      ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1 = 1, pc.2 = 2)
PoV.cor <- round(100 * attr(pca.cor, "percentVar"))
PCAplot(pca.cor, "batch", PoV.df=PoV.cor, colors=c('orange','lightblue','lightgreen',"yellow"), pc.1 = 1, pc.2 = 2)
PCAplot(pca.cor, "smoking", PoV.df=PoV.cor, colors=c('orange','lightblue','lightgreen',"yellow"), pc.1 = 1, pc.2 = 2)
PCAplot(pca.cor, "status", PoV.df=PoV.cor, pc.1 = 1, pc.2 = 2) #very clean pca as result

#Checking inter-sample distances after correction
annot_df <- annot_df[!rownames(annot_df) %in% c("36", "37"),]
ComplexHeatmap::pheatmap(cor(batch.rem), annotation_row = annot_df, annotation_col = annot_df,
                         annotation_colors = annot_colors, color = heatmap.color.code) #extremely high correlations > .99 instead of .9



#==========================================================================================================#
#                                            RNA-seq analysis: DESeq2                                      #
#==========================================================================================================#
#dds <- DESeq(dds)
setwd("/Users/kubler01/Documents/R projects/R workspace/")
load("12282021_DESeq2_Monocytes.RData")

resultsNames(dds) # lists the coefficients

#res.base <- results(dds.base)
old.res <- results(dds, name="status_1_vs_0")
res <- old.res

#sum(res.base$padj < 0.05, na.rm=TRUE)
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$log2FoldChange < -1 & res$padj < 0.05, res$log2FoldChange > 1 & res$padj < 0.05, na.rm=TRUE)

res[is.na(res$padj),]$padj <- 1
downreg.genes <- rownames(res[res$log2FoldChange < -1 & res$padj < 0.1,])
upreg.genes <- rownames(res[res$log2FoldChange > 1 & res$padj < 0.1,])

length(rownames(res[res$log2FoldChange > 0.5 & res$padj < 1,]))
length(rownames(res[res$log2FoldChange < -0.5 & res$padj < 1,]))

downreg.genesALL <- rownames(res[res$log2FoldChange < 0 & res$padj < 0.1,])
upreg.genesALL <- rownames(res[res$log2FoldChange > 0 & res$padj < 0.1,])


#Plotting top DEGs
ComplexHeatmap::pheatmap(batch.rem[rownames(batch.rem) %in% upreg.genesALL,], scale= "row",
                         cluster_rows = T, cluster_cols = T, annotation_legend = T, show_colnames = F, show_rownames = T,
                         legend = T, treeheight_row = 0, color = heatmap.color.code,
                         annotation_col=annot_df, annotation_colors = annot_colors, fontsize_row = 8)

ComplexHeatmap::pheatmap(batch.rem[rownames(batch.rem) %in% downreg.genesALL,], scale= "row",
                         cluster_rows = T, cluster_cols = T, annotation_legend = T, show_colnames = F, show_rownames = T,
                         legend = T, treeheight_row = 0, color = heatmap.color.code,
                         annotation_col=annot_df, annotation_colors = annot_colors, fontsize_row = 8)


#Comparing lFCs to other studies
setwd(file.path("/Users/kubler01/Documents/R projects/Monocyte project/Gene lists"))
PBMC_SZ <- read_excel("PBMCs in SZ.xlsx", sheet=3)

PBMC.DEGs <- PBMC_SZ[PBMC_SZ$Names %in% c(downreg.genesALL, upreg.genesALL),]

res[rownames(res) %in% PBMC.DEGs$Names,]

#FAM118A is an expression QTL in AD shared between monocytes and microglia
#DUSP2 is an inflammatory gene previously found to be increased in monocytes of SCZ patients, correlated with cytokine clusters including IL1, IL6, and TNF



#==========================================================================================================#
#                                RNA-seq analysis: variance causa testing                                  #
#==========================================================================================================#
#rownames(batch.rem) <- genemap.mic$external_gene_name
#Merging gene panels
Gandal.df <- data.frame()
for (i in c("Astrocytes", "NFkB")){
  type <- Gandal.clusters.annot[Gandal.clusters.annot$Module %in% i,]
  Gandal.df <- data.frame(rbind.fill(Gandal.df, data.frame(t(type$gene_name))), check.names=F)}

Gandal.df <- data.frame(t(Gandal.df))
colnames(Gandal.df) <- c("Astrocyte genes", "NFkB signaling")

colnames(Patir) <- "Microglia genes"
SCZ.risk <- SCZ.risk[SCZ.risk %in% rownames(batch.rem)]
gene.panels_mrgd <- rbind.fill(data.frame("Microglia genes"=Patir, check.names = F), Gandal.df, data.frame("IFNy response"=gene.panels[,2], check.names = F),
                               gene.panels_2, data.frame("Schizophrenia risk genes"=SCZ.risk, check.names = F))


#plot pca
for (i in colnames(gene.panels_mrgd)){
 dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,i]),]
#  plots <- PCA_cov_cor_R(clean_covs, dat) #if 
  list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                             ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
  list.pov <- round(100 * attr(list.pca, "percentVar"))
  plot <- PCAplot(list.pca, "status", PoV.df=list.pov, pc.1=1, pc.2=2)
  print(plot)}




#==========================================================================================================#
#                                    RNA-seq analysis: enrichment testing                                  #
#==========================================================================================================#
GO.up <- enrichGO(gene    = upreg.genesALL,
         universe      = rownames(res),
         OrgDb         = "org.Hs.eg.db",
         keyType       = 'SYMBOL',
         ont           = "BP",
         pAdjustMethod = "BH",
         pvalueCutoff  = 0.05,
         qvalueCutoff  = 0.1,
         readable      = F, pool=T)

GO.down <- enrichGO(gene    = downreg.genesALL,
                    universe      = rownames(res),
                    OrgDb         = "org.Hs.eg.db",
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1,
                    readable      = F, pool=T)

checkGO<-GO.up
checkGO@result<-checkGO@result[checkGO@result$p.adjust < 0.1 & checkGO@result$ONTOLOGY %in% "BP",]
dotplot(checkGO, split="ONTOLOGY", showCategory=10) + facet_grid(ONTOLOGY~., scale="free")
dotplot(GO.down, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

data.frame(GO.up)
#Subset gene panels to background if not yet done
m <- data.frame(t(gene.panels_mrgd$`Microglia genes`[gene.panels_mrgd$`Microglia genes` %in% rownames(res)]), check.names = F)
for(i in names(gene.panels_mrgd)[-1]){
  k <- na.omit(gene.panels_mrgd[i])
  s <- k[k[,1] %in% rownames(res),]
  m <- rbind.fill(data.frame(m),data.frame(t(s)))
  }

rownames(m) <- names(gene.panels_mrgd)
gene_panels_subset <- data.frame(t(m),check.names=F)

overlap.genesets <- matrix(ncol=8,nrow=8)
rownames(overlap.genesets) <- colnames(gene_panels_subset)
colnames(overlap.genesets) <- colnames(gene_panels_subset)

for(i in colnames(overlap.genesets)){
  for (x in colnames(overlap.genesets))
    overlap.genesets[x,i] <- sum(na.omit(gene_panels_subset[,i]) %in% na.omit(gene_panels_subset[,x]))/length(na.omit(gene_panels_subset[,i]))}
#overlap.genesets <- as.dist(overlap.genesets)

ComplexHeatmap::pheatmap(overlap.genesets, cluster_rows = F, cluster_cols = F)
ComplexHeatmap::Heatmap(as.matrix(overlap.genesets),
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colorRampPalette(brewer.pal(n = 6, name ="RdPu"))(30),
                        show_column_names = T,
                        column_names_rot = 90,
                        name = "Ratio overlap",
                        #        clustering_method_rows = "ward.D2",
                        #        clustering_distance_rows = "pearson"
)


#Manual gene list enrichment testing
DEG.enrich.res <- matrix(nrow=2*length(colnames(gene_panels_subset)), ncol=6)
DEG.rep <- NULL
for (i in names(gene_panels_subset))
  DEG.rep <- c(DEG.rep, paste(i, "_up", sep=""), paste(i, "_down", sep=""))
DEGs <- list("Up-regulated"=upreg.genesALL,"Down-regulated"=downreg.genesALL)
rownames(DEG.enrich.res) <- DEG.rep
colnames(DEG.enrich.res ) <- colnames(GSEA.byMod(mod.gl=DEGs, gene_panels_subset[,i], universe=35324))


for (x in colnames(gene_panels_subset)){
  DEG.enrich.res[gsub("*_up","",gsub("*_down","",rownames(DEG.enrich.res))) %in% x,] <- as.matrix(GSEA.byMod(mod.gl=DEGs, na.omit(gene_panels_subset[,x]), universe=rownames(res)))
}

DEG.enrich.res <- data.frame(DEG.enrich.res)
DEG.enrich.res$log.q <- -log2(DEG.enrich.res[,6])

#DEG.enrich.res[is.infinite(DEG.enrich.res$log.q),]$log.q <- 1

DEG.enrich <- cbind("Upregulated"=DEG.enrich.res[c(seq(1, 15, by = 2)),][,7],"Downregulated"=DEG.enrich.res[c(seq(2, 16, by = 2)),][,7])
rownames(DEG.enrich) <- names(gene_panels_subset)
ComplexHeatmap::Heatmap(DEG.enrich,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30),
                        name = "log2(q)",
)


#==========================================================================================================#
#                         RNA-seq analysis: co-variance DGX and pathway genes                              #
#==========================================================================================================#
#=========Checking the differential regulations of our panels==========#
panel.list <- NULL
for(i in colnames(gene_panels_subset))
  panel.list[[i]] <- batch.rem[rownames(batch.rem) %in% na.omit(gene_panels_subset[,i]),]

DE.list <- NULL
for(i in colnames(gene_panels_subset))
  DE.list[[i]] <- data.frame(res[rownames(res) %in% na.omit(gene_panels_subset[,i]),])


#Interesting, let's plot the most FC genes
for (i in colnames(gene_panels_subset)){
  panel <- panel.list[[i]]
  DE <- DE.list[[i]]
  df <- panel[DE$log2FoldChange < -0.5 | DE$log2FoldChange > 0.5,]
  df <- df[,order(as.numeric(colnames(df)), method="radix", decreasing=F)]
  ComplexHeatmap::pheatmap(df, scale= "row",
                           cluster_rows = T, cluster_cols = T, annotation_legend = T, show_colnames = F, show_rownames = T,
                           legend = T, treeheight_row = 0, color = heatmap.color.code,
                           annotation_col=annot_df, annotation_colors = annot_colors, fontsize_row = 8, column_title=paste("DEGs within", i, "list"))
}


#Calculate compound lfc score
#Do the same for the lFCs: compound the genes lists' lFC score and throw on a heatmap
compoundlFC <- mean(DE.list[[1]][order(abs(DE.list[[1]]$log2FoldChange), decreasing=T),][1:10,]$log2FoldChange)
for (i in names(DE.list)[-1]){
  top10 <- DE.list[[i]][order(abs(DE.list[[i]]$log2FoldChange), decreasing=T),][1:10,]
  compoundlFC <- rbind(compoundlFC,mean(top10$log2FoldChange))}
rownames(compoundlFC) <- names(DE.list)

#median lFC to point directions
ComplexHeatmap::Heatmap(as.matrix(compoundlFC),
                        cluster_rows = T,
                        cluster_columns = F,
                        show_row_dend = T,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = T,
                        name = "Mean log(FC)",
                        #        clustering_method_rows = "ward.D2",
                        #        clustering_distance_rows = "pearson"
)


#===============correlation between gene expression and status=================#
#Differential PC expression of list expression summaries

#The first question we can ask:
#Does the gene variance differ between status groups (i.e., PC expression)?
#Should we go the stimulation PC for that?
#We could by setting a threshold


#Creating dataframe with first PC for each gene list
scatter.df <- data.frame("Status"=as.factor(clean_covs$status))
for (i in names(panel.list)){
  pca <- plotPCA.custom(as.matrix(panel.list[[i]]), intgroup=c("status", "batch", "sex", "RIN"), 
                        ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)[,1]
  scatter.df <- cbind(scatter.df, pca)}

colnames(scatter.df) <- c("Status", names(panel.list))


GEX.stat.cor <- matrix(nrow=length(colnames(scatter.df[,-1])), ncol=4)
colnames(GEX.stat.cor) <- c("Spearman.rho", "P", "lm.b", "Pr(>|t|)")
rownames(GEX.stat.cor) <- colnames(scatter.df[,-1])

for (i in colnames(scatter.df[,-1])){
  GEX.stat.cor[,1][[i]] <- rcorr(as.matrix(scatter.df[,i]), as.matrix(clean_covs$status), type = "spearman")$r[2]
  GEX.stat.cor[,2][[i]] <- rcorr(as.matrix(scatter.df[,i]), as.matrix(clean_covs$status), type = "spearman")$P[2]
  GEX.stat.cor[,3][[i]] <- coef( lm(scatter.df[,i] ~ as.matrix(clean_covs$status)))[2]
  GEX.stat.cor[,4][[i]] <- data.frame(coef(summary( lm(scatter.df[,i] ~ as.matrix(clean_covs$status)))))$Pr...t..[2]
}


#Put into heatmap:
matrix_pvalue = cbind(GEX.stat.cor[,2],GEX.stat.cor[,4])
matrix_pvalue = matrix(p.adjust(as.vector(as.matrix(matrix_pvalue)), method='bonferroni'),ncol=ncol(matrix_pvalue))
matrix_pvalue = formatC(matrix_pvalue, format="e", digits = 2)
rownames(matrix_pvalue) <- rownames(GEX.stat.cor)
  
  # png(paste0(work_dir, "LinearReg_15pcs_covariates.png"), width = 10, height = 10, res = 600, units = "in")
# pdf(paste0(work_dir, "LinearReg_15pcs_covariates.pdf"), width = 10, height = 10)
ComplexHeatmap::Heatmap(as.matrix(GEX.stat.cor[,1]),
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = T,
                        column_names_rot = 90, 
                        name = "Spearman correlation",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(as.matrix(matrix_pvalue[,2])[i,j], x, y, gp = gpar(fontsize = 10))
                        }
                        #        clustering_method_rows = "ward.D2",
                        #        clustering_distance_rows = "pearson"
)

ComplexHeatmap::Heatmap(as.matrix(GEX.stat.cor[,3]),
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = T,
                        column_names_rot = 90, 
                        name = "Beta effect size",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(as.matrix(matrix_pvalue[,2])[i,j], x, y, gp = gpar(fontsize = 10))
                        }
                        #        clustering_method_rows = "ward.D2",
                        #        clustering_distance_rows = "pearson"
)


######representation of this data in boxplots (differential PC expression)################
formula <- y ~ x
max(scatter.df[,-1]
    )
library(ggpmisc)
library(ggbeeswarm)
library(ggpubr)
scatters <- NULL
for (i in colnames(scatter.df)){
  if(!i %in% "Status"){
    scatters[[i]] <- ggplot(data.frame(scatter.df), aes_string(x = "Status", y = i, fill = "Status")) + 
      geom_boxplot() +
      #  scale_fill_manual(values = cols)+
      theme_classic()+
      theme(plot.title = element_text(hjust=0.5))+
      ylab(paste(i,"GEX variance"))+
      theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
    geom_smooth(method = lm, linetype = "dashed") + 
    stat_compare_means(aes(label = ..p.signif..),method = "t.test", label.x = 1.5, label.y=14)}
}
library(ggpubr)
PC.boxplots <- ggarrange(plotlist=scatters[1:length(scatters)], common.legend = T)
PC.boxplots



#=======2nd question=====#

#We ask different questions here, the following plots answer: do the signatures correlate?
#It basically tells us: how to gene programs interact with each other on a dimensionally general level (gauging networks potentially)


#A third question might be where we could even ask if the correlations differ between groups


#Check whether DEGs are co-regulated with gene lists

#Data summarized as heatmap
status.PCs <- rcorr(as.matrix(cbind(scatter.df[,-1],"Status"=scatter.df[,1])))$r
length(colnames(status.PCs))
status.PCs <- data.frame(status.PCs[1:7,8])
colnames(status.PCs) <- "Status"

ComplexHeatmap::Heatmap(as.matrix(status.PCs),
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = T,
                        column_names_rot = 90, 
                        name = "Pearson correlation",
                        #        clustering_method_rows = "ward.D2",
                        #        clustering_distance_rows = "pearson"
)


#=========== third question: does status influence co-regulation of signatures?=============#
glm.df <- scatter.df

glance(lm(scatter.df$Patir_microglia ~ as.matrix(clean_covs$status) + as.matrix(scatter.df[,i]) + as.matrix(clean_covs$status):as.matrix(scatter.df[,i])))
z <- matrix(nrow=length(scatter.df[,-1]), ncol=length(scatter.df[,-1]))
colnames(z) <- colnames(scatter.df[,-1])
rownames(z) <- colnames(z)

P <- z

for(i in colnames(scatter.df[,-1])){
  for(x in colnames(scatter.df[,-1])){
    if (paste(x) != paste(i)){
    model.matrix(lm.stat <- glm(scatter.df$Status ~ as.matrix(scatter.df[,i]) + as.matrix(scatter.df[,x]) + as.matrix(scatter.df[,i]):as.matrix(scatter.df[,x]),  family=binomial))
    z[,paste(i)][paste(x)] <- coef(summary(lm.stat))[,3][4]
    P[,paste(i)][paste(x)] <- coef(summary(lm.stat))[,4][4]
    }
  }
}

P <- formatC(P, format = "e", digits = 2)

ComplexHeatmap::Heatmap(as.matrix(z),
                        cluster_rows = T,
                        cluster_columns = T,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = T,
                        column_names_rot = 90, 
                        name = "Beta effect size",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(as.matrix(P)[i,j], x, y, gp = gpar(fontsize = 10))
                        }
                        #        clustering_method_rows = "ward.D2",
                        #        clustering_distance_rows = "pearson"
)

#============End of analysis===============#

#==========Paper figures=====================#
#a) volcano plot
# Function Vulcano plot 
volcano_plot <- function(res, title = NULL, subtitle = NULL, annotate_by = NULL, type ='ALS'){
  res <- 
    mutate(res,
           sig = case_when(
             padj >= 0.1 ~ "non_sig",
             padj < 0.1 & abs(log2FoldChange) < 1 ~ "sig",
             padj < 0.1 & abs(log2FoldChange) >= 1 ~ "sig - strong"
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
    geom_text(fontface = "bold", data = de_tally, aes(x = position, y = ymax - 0.5, label = n, colour = class), size = 2.5 ) +
    scale_x_continuous(limits = xlim)
  if(!is.null(annotate_by)){
    plot <- plot + 
      ggrepel::geom_text_repel(
        fontface = "italic",
        data = filter(res, symbol %in% annotate_by), 
        aes(x = log2FoldChange, y = -log10(pvalue), label = symbol), 
        min.segment.length = unit(0, "lines"),
        size = 2.5) +
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
volcano_plot(res.wnam, annotate_by=c('FES',"IL6", "IL7", "INF2", "P2RY6", "TNFAIP3", "CSMD1", "MMP9", "LRRC2"))


# "NFKBIZ", "NFKBIA",
#b) DEG heatmap

top20 <- rownames(res[abs(res$log2FoldChange) > 0 & res$padj < 0.1,][1:20,])
top20.dat <- batch.rem[rownames(batch.rem) %in% top20,]

rev.covs <- annot_df[order(length(rownames(annot_df)):1),]
top20.dat <- top20.dat[,order(ncol(top20.dat):1)]
identical(colnames(top20.dat), rownames(rev.covs))

ha <- HeatmapAnnotation(Status = rev.covs$Status,
                  col = list( Status = c("CTR" = "lightblue", "SCZ" = "tomato")))

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
PCAplot(list.pca, "Status", PoV.df=list.pov, pc.1=1, pc.2=2, colors=c("lightblue", "tomato"), geom.size=4)
ggplot(list.pca[,c(1,4)], 
       aes_string(x = "Status", y = "PC1", fill = "Status")) + 
  geom_boxplot() +
  #  scale_fill_manual(values = cols)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))+
  ylab(paste("PC1 expression NFkB panel"))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
  geom_smooth(method = lm, linetype = "dashed") + 
  stat_compare_means(aes(label = ..p.signif..),method = "t.test", label.x = 1.5, label.y=14)+
  scale_fill_manual(c("CTR", "SCZ"),values=c("lightblue", "tomato"))
t.test(x=list.pca[list.pca$Status %in% 'CTR',]$PC1,
       y=list.pca[list.pca$Status %in% 'SCZ',]$PC1, alternative='two.sided')


#g) LPS PCA plot and boxplot
dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"LPS response"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
PCAplot(list.pca, "Status", PoV.df=list.pov, pc.1=1, pc.2=2, colors=c("lightblue", "tomato"), geom.size=4)
ggplot(list.pca[,c(2,4)], 
       aes_string(x = "Status", y = "PC2", fill = "Status")) + 
  geom_boxplot() +
  #  scale_fill_manual(values = cols)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))+
  ylab(paste("PC2 expression LPS panel"))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
  geom_smooth(method = lm, linetype = "dashed") + 
  stat_compare_means(aes(label = ..p.signif..),method = "t.test", label.x = 1.5, label.y=14)+
  scale_fill_manual(c("CTR", "SCZ"),values=c("lightblue", "tomato"))

t.test(x=list.pca[list.pca$Status %in% 'CTR',]$PC1,
             y=list.pca[list.pca$Status %in% 'SCZ',]$PC1, alternative='two.sided')


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
Gandal <- Gandal[Gandal$SCZ.fdr < 0.05,]
MIGA <- read_excel("SCZ MIGA QTL genes.xlsx")
MIGA$QTL_Gene

table(MIGA$QTL_Gene %in% c(DEGs$`Up-regulated`,DEGs$`Down-regulated`))

log10(MIGA[MIGA$QTL_Gene %in% "FES",]$QTL_P)

MIGA.FES <- as.matrix(data.frame("Beta"=MIGA[MIGA$QTL_Gene %in% "FES",]$QTL_Beta,check.names=F))

rownames(MIGA.FES) <- MIGA[MIGA$QTL_Gene %in% "FES",]$QTL_SNP
ra.FES <- rowAnnotation(`QTL type` = c(rep('Microglia sQTL',5),rep('Monocyte eQTL',2)), `Origin`=c(rep('MIGA MFG',2),rep("MIGA STG",3),'Fairfax','MyND'),
                        col = list(`QTL type` = c("Monocyte eQTL" = "#C71000FF", "Microglia sQTL" = "#008EA0FF"), 
                                   `Origin`=c("MIGA MFG"="#E66101","MIGA STG"="#FDB863","Fairfax"="#B2ABD2","MyND"="#5E3C99")))#colorRampPalette(RColorBrewer::brewer.pal(4,"PuOr"))(4)
ComplexHeatmap::Heatmap(MIGA.FES[,1],
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        show_column_names = F,
                        col = colors,
                        left_annotation = ra.FES,
                        column_names_rot = 60, 
                        name = "Beta",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(as.matrix(formatC(MIGA[MIGA$QTL_Gene %in% "FES",]$QTL_P, format = "e", digits = 2))[i,j], x, y, gp = gpar(fontsize = 10))
                        }
)
#For Drexhage FCs
library(ggpubr)
library(tidyr)
library(ggeasy)
library(ggrepel)
Drexhage.lfc <- drop_na(left_join(Drexhage,
                                  data.frame("Gene"=rownames(res[rownames(res) %in% Drexhage$Gene,]),
                              'Monocytes' = res[rownames(res) %in% Drexhage$Gene,]$log2FoldChange, check.names=F), by='Gene'))

#correlation plot of log FCs
symbols.label.drex <- data.frame('overlap'=Drexhage.lfc$Gene)
symbols.label.drex$overlap[(!symbols.label.drex$overlap %in% c(upreg.genesALL, downreg.genesALL, 'TNF', 'CCL7','F3'))] = ""

lFC.scat[['Drex']] <- ggplot2::ggplot(Drexhage.lfc, aes(x = log2(lFC), y = Monocytes)) + #Drexhage genes are in FC, not lFC
  geom_point(alpha = .5) +
  scale_color_manual(values = c("#C71000FF", "#8A4198FF", "#008EA0FF")) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") + 
  stat_smooth(method = "lm", se=F) + # Add Regression Line 
  stat_regline_equation(aes(label = ..adj.rr.label..), show.legend = F)  + # Add R-Square
  #  stat_regline_equation(aes(label = ..rr.label..))  +
  geom_label_repel(aes(label = symbols.label.drex$overlap), size = 3, color="black", 
                   box.padding = 0.4, label.size = NA, fill = alpha(c("white"),0.5), max.overlaps = 15)+
  easy_labs(x = expression(paste("Drexhage fold change (", log2,")")), y = expression(paste("Monocyte fold change (", log2, ")"))) +
  easy_add_legend_title("Drexhage schizophrenia genes") +
  theme_classic()




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
  geom_label_repel(aes(label = symbols.label$overlap), size = 3, color="black", 
                   box.padding = 0.4, label.size = NA, fill = alpha(c("white"),0.5), max.overlaps = 13)+
  easy_labs(x = expression(paste("Gandal fold change (", log2,")")), y = expression(paste("Monocyte fold change (", log2, ")"))) +
  theme_classic()




#==========PDF paper figure============#
ra <- ComplexHeatmap::rowAnnotation(Status = rev.covs$Status,
                                    col = list(Status = c("CTR" = "lightblue", "SCZ" = "tomato")))
checkGO<-GO.up
checkGO@result<-checkGO@result[checkGO@result$p.adjust < 0.1 & checkGO@result$ONTOLOGY %in% "BP",]

dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"NFkB signaling"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
#PCAplot(list.pca, "Status", PoV.df=list.pov, pc.1=1, pc.2=2, colors=c("lightblue", "tomato"), geom.size=4)
boxplots<-list("LPS", 'NFkB','GC')
boxplots[['LPS']] <- ggplot(list.pca[,c(1,4)], 
                            aes_string(x = "Status", y = "PC1", fill = "Status")) + 
  geom_boxplot() +
  #  scale_fill_manual(values = cols)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))+
  ylab(paste("PC1 NFkB signaling expression"))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
  geom_smooth(method = lm, linetype = "dashed") + 
  stat_compare_means(aes(label = ..p.signif..),method = "t.test", label.x = 1.5, label.y=14)+
  scale_fill_manual(c("CTR", "SCZ"),values=c("lightblue", "tomato"))

dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"LPS response"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
#PCAplot(list.pca, "Status", PoV.df=list.pov, pc.1=1, pc.2=2, colors=c("lightblue", "tomato"), geom.size=4)
boxplots[['NFkB']] <- ggplot(list.pca[,c(2,4)], 
                             aes_string(x = "Status", y = "PC2", fill = "Status")) + 
  geom_boxplot() +
  #  scale_fill_manual(values = cols)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))+
  ylab(paste("PC2 LPS response expression"))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
  geom_smooth(method = lm, linetype = "dashed") + 
  stat_compare_means(aes(label = ..p.signif..),method = "t.test", label.x = 1.5, label.y=14)+
  scale_fill_manual(c("CTR", "SCZ"),values=c("lightblue", "tomato"))

dat <- batch.rem[rownames(batch.rem) %in% na.omit(gene.panels_mrgd[,"Glucocorticoid response"]),]
list.pca <- plotPCA.custom(as.matrix(dat), intgroup=c("status", "batch", "sex", "RIN"), 
                           ntop = 50000, returnData=TRUE, metadata=clean_covs, pc.1=1, pc.2=2)
list.pca[,c(1:2)] <- -list.pca[,c(1,2)]
colnames(list.pca)[4] <- "Status"
list.pca$Status <- as.factor(c(rep("SCZ",20),rep("CTR",17)))
list.pov <- round(100 * attr(list.pca, "percentVar"))
#PCAplot(list.pca, "Status", PoV.df=list.pov, pc.1=1, pc.2=2, colors=c("lightblue", "tomato"), geom.size=4)
boxplots[['GC']] <- ggplot(list.pca[,c(2,4)], 
                           aes_string(x = "Status", y = "PC2", fill = "Status")) + 
  geom_boxplot() +
  #  scale_fill_manual(values = cols)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))+
  ylab(paste("PC2 Glucocorticoid response expression"))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
  geom_smooth(method = lm, linetype = "dashed") + 
  stat_compare_means(aes(label = ..p.signif..),method = "t.test", label.x = 1.5, label.y=14)+
  scale_fill_manual(c("CTR", "SCZ"),values=c("lightblue", "tomato"))

PCs <- ggarrange(plotlist=boxplots[1:length(boxplots)], common.legend = T)




                       ## 2x3 grid of plotting areas
#Plot A
p1 <- volcano_plot(res.wnam, annotate_by=c('FES',"IL6", "IL7", "INF2", "P2RY6", "TNFAIP3", "CSMD1", "MMP9", "LRRC2"))
#Plot B
lfc.grob <- list(grid.grabExpr(plot(lFC.scat$Gandal)),grid.grabExpr(plot(lFC.scat$Drex)))
p2 <- ggarrange(plotlist=lfc.grob[1:length(lfc.grob)], nrow = 2)
#Plot C
p3 <- ComplexHeatmap::Heatmap(t(scale(t(top20.dat))),
                        cluster_rows = T,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colors,
                        show_column_names = F,
                        column_names_rot = 0, 
#                        right_annotation  = ra,
                        top_annotation = ha,
                        name = "Normalized counts",
                        clustering_method_rows = "ward.D2",
                        clustering_distance_rows = "pearson"
)
#Plot D
p4 <- dotplot(checkGO, split="ONTOLOGY", showCategory=10) + facet_grid(ONTOLOGY~., scale="free")
#Plot E
p5 <- ComplexHeatmap::Heatmap(t(DEG.enrich),
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = T,
                        col = colorRampPalette(RColorBrewer::brewer.pal(6,"RdPu"))(30),
                        column_names_rot = 60, 
                        name = "log2-q",
                        )
#Plot F
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

tiff("03022022_Monocyte figure.tiff", width=5000, height=4000, res=300)
grid.arrange(grobs = grobs, layout_matrix = lay)
dev.off()


#Separating the figure into two
lfc.grob <- list(grid.grabExpr(plot(lFC.scat$Gandal)),grid.grabExpr(plot(lFC.scat$Drex)))
p3.grid <- grid.grabExpr(draw(p3))
p5.grid <- grid.grabExpr(draw(p5))
p6.grid <- grid.draw(grid.grabExpr(draw(p6)))

lay1 <- rbind(c(1,1,2,3),
             c(1,1,2,3),
             c(4,4,2,3),
             c(4,4,5,5),
             c(6,6,6,6))

grobs1 <- list(p1,p2,p3.grid,p4,p5.grid,p6)
grid.arrange(grobs = grobs, layout_matrix = lay1)

lay2 <- rbind(c(1,1,2,2),
              c(1,1,3,3),
              c(4,4,4,4),
              c(5,5,5,5))

grobs2 <- list(p1,lfc.grob[[2]],p5.grid,p4,p6)
grid.arrange(grobs = grobs2, layout_matrix = lay2)

#=======END=======#
lay <- rbind(c(1,1,2,3),
             c(1,1,2,3),
             c(1,1,2,3),
             c(3,3,5,5),
             c(3,3,5,5),
             c(3,3,6,6),
             c(3,3,6,6),
             c(4,4,6,6))
pdf('eg.pdf', width = 8.3, height = 11.7)  ## Device with dimensions of A4 paper
par(omi = rep(.5, 4))                      ## 1/2 inch outer margins
par(mfrow = c(3,2)) 
p1+
  p2+
  grid.draw(grid.grabExpr(draw(p3)))+
  grid.draw(grid.grabExpr(plot(p4)))+
  grid.draw(grid.grabExpr(draw(p5)))+
  p6
dev.off()

grid.draw(grid.grabExpr(plot(p4)))


library(gridExtra)
library(ComplexHeatmap)
library(multipanelfigure)
#210 x 297 mm
(figure <- multi_panel_figure(
  width = c(104,104),
  height = c(110,110,77),
  panel_label_type = "upper-roman", row_spacing = 1, column_spacing = 1))


# Fill the top-left panel using a grob object directly
#a_grob <- grid::linesGrob(arrow = grid::arrow())
figure %<>% fill_panel(p1)
figure %<>% fill_panel(p2)
figure %<>% fill_panel(p3)
figure %<>% fill_panel(p4)
figure %<>% fill_panel(p5)
figure %<>% fill_panel(p6)

save_multi_panel_figure(figure, filename='multipanel', dpi = 300)

p3.grid <- grid.grabExpr(draw(p3))
p5.grid <- grid.grabExpr(draw(p5))
p6.grid <- grid.draw(grid.grabExpr(draw(p6)))

p124 <- plot_grid(p1, p2, NULL, p4, NULL, p6, nrow = 3, rel_heights = c(1/4, 1/4, 1/4, 1/4),
                  ncol = 2, rel_widths = c(1/2,1/2))
p124



gs <- lapply(1:6, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))





ggplot2::ggplot(MIGA.FES,
                aes_string(x = "Beta", y = "`Cell type`", fill = "`Cell type`")) + #Drexhage genes are in FC, not lFC
  geom_point(alpha = .5) +
  scale_color_manual(values = c("#C71000FF", "#8A4198FF", "#008EA0FF")) +
  #geom_hline(yintercept = -5, linetype = "dashed", colour = "black") +
  #geom_vline(xintercept = 0, linetype = "dashed", colour = "black") + 
  stat_smooth(method = "lm", se=F) + # Add Regression Line 
  stat_regline_equation(aes(label = ..adj.rr.label..), show.legend = T)  + # Add R-Square
  stat_regline_equation(aes(label = ..rr.label..))  +
  #geom_label_repel(aes(label = symbols.label.drex$overlap), size = 3, color="black", box.padding = 0.4, label.size = NA, fill = alpha(c("white"),0.5))+
  easy_labs(x = 'Beta', y = 'log(P)') +
  easy_add_legend_title("Drexhage schizophrenia genes") +
  theme_classic()

ggplot(MIGA.FES, aes_string(y = "P", x = "`Cell type`", fill = "`Cell type`")) + 
  geom_violin() +
  #  scale_fill_manual(values = cols)+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))+
  ylab(paste(i,"GEX variance"))+
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())+
  geom_smooth(method = lm, linetype = "dashed") + 
  stat_compare_means(aes(label = ..p.signif..),method = "t.test", label.x = 1.5, label.y=14)


