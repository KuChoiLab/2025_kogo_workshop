############
# FigR 
############
options(bitmapType = "cairo")

# 1. Loading ATAC-seq and RNA-seq data ---------
library(doParallel)
library(BuenColors)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg19)

setwd("~/2025_KOGO_workshop/ATACseq")
ATAC.se <- readRDS("data/FigR/control1h_PBMC_atac_SE.rds")
RNAmat <- readRDS("data/FigR/control1h_PBMC_RNAnorm.rds")
CCA_PCs <- readRDS("data/FigR/control1h_PBMC_atac_rna_CCA_l2.rds")

dim(ATAC.se) # Peaks x ATAC cells
dim(RNAmat) # Genes x RNA cells
dim(CCA_PCs) # ATAC + RNA (rows), n components (columns)
head(rownames(CCA_PCs)) # ATAC cells
tail(rownames(CCA_PCs)) # RNA cells

isATAC <- grepl("BC",rownames(CCA_PCs))
table(isATAC) # ATAC vs RNA
ATACcells <- rownames(CCA_PCs)[isATAC]
RNAcells <- rownames(CCA_PCs)[!isATAC]
length(ATACcells)
length(RNAcells)

# 2. Visualizing UMAP of the ATAC-RNA cell ---------
nPCs <- 20 # Num CCA PCs to use when running UMAP / pairing
set.seed(123)
umap.out <- uwot::umap(CCA_PCs[,1:nPCs],
                       metric="cosine",
                       n_neighbors=30)
umap.d <- as.data.frame(umap.out)
colnames(umap.d) <- c("UMAP1","UMAP2")
rownames(umap.d) <- rownames(CCA_PCs)
umap.d$Assay <- ifelse(isATAC,"ATAC","RNA")

BuenColors::shuf(umap.d) %>%
  ggplot(aes(UMAP1,UMAP2,color=Assay)) +
  geom_point(size=0.1) +
  theme_classic() + theme_classic() +
  scale_color_manual(values = c("cadetblue","darkorange"))+
  guides(colour = guide_legend(override.aes = list(size=3)))


# 3. Pairing cells using scOptMatch ---------
# get PCs for each data
ATAC_PCs <- CCA_PCs[isATAC,]
RNA_PCs <- CCA_PCs[!isATAC,]
dim(ATAC_PCs)
dim(RNA_PCs)
# pair cells using scOptMatch
pairing <- pairCells(ATAC = ATAC_PCs,
                     RNA = RNA_PCs,
                     keepUnique = TRUE)
dim(pairing)

# 4. Visualizing ATAC-RNA pairs on the CCA UMAP --------
library(ggrastr)
plotPairs(ATAC = pairing$ATAC,
          RNA=pairing$RNA,
          max.show = 100,
          umap.df = umap.d)

# 5. Getting count object for the ATAC-RNA paired cells ---------
ATAC.se.paired <- ATAC.se[,pairing$ATAC]
RNAmat.paired <- RNAmat[,pairing$RNA]
dim(ATAC.se.paired)
dim(RNAmat.paired)

# 6. Peak-gene association testing --------
# Compute correlation between RNA expression and peak accessibility for peaks falling within a window around each gene. 

# Do not run this code
#cisCorr <- runGenePeakcorr(ATAC.se = ATAC.se.paired,
#                           RNAmat = RNAmat.paired,
#                           genome = "hg19",
#                           nCores = 4,
#                           p.cut = NULL, 
#                           n_bg = 100) # the number of background correlations to compute per gene-peak pair

cisCorr <- readRDS("data/FigR/control1h_PBMC_cisCor.rds")
# head for cisCorr 
head(cisCorr)

# 7. Determining DORCs --------
# filtering DORCs by p-value
cisCorr.filt <- cisCorr %>% dplyr::filter(pvalZ <= 0.05)
cisCorr.filt %>% dplyr::arrange(desc(pvalZ)) %>% head()

# plotting DORCs -------
library(ggplot2)
dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                       cutoff = 10,
                       labelTop = 20,
                       returnGeneList = TRUE,
                       force=5)
head(dorcGenes)
length(dorcGenes)

# 8. Summary of DORC ---------
# to get the DORC accessibility scores,
# we can sum up the chromatin accessibility peak counts for peaks associated 
numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
numDorcs

# 9 Calculating DORC scores -------
# calculate DORC scores
dorcMat <- getDORCScores(ATAC.se = ATAC.se.paired,
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = 2)
dim(dorcMat)
dorcMat[1:2,10:20]

# colData(ATAC.se.paired) %>% 
#   as.data.frame() %>% 
#   ggplot(aes(UMAP1,UMAP2)) + 
#   geom_point(color="gray",size=0.8)+ theme_classic()


# 10. Smoothing RNA using cell KNNs ------
lsi <- readRDS("data/FigR/control1h_PBMC_atac_lsi.rds")
dim(lsi)

all(colnames(dorcMat) %in% rownames(lsi))

# Subset to paired ATAC
length(colnames(dorcMat))
head(colnames(dorcMat))
lsi <- lsi[colnames(dorcMat),]

dim(lsi)

# For 30 LSIs, get the nearest cell for each cell
cellkNN <- FNN::get.knn(lsi,k=30)$nn.index
dim(cellkNN)# 4912 cells, 30 LSIs
rownames(cellkNN) <- colnames(dorcMat)

# Smooth dorc scores using cell KNNs (k=30)
library(doSNOW)
library(doParallel)

dorcMat.s <- smoothScoresNN(NNmat = cellkNN, mat = dorcMat, nCores = 2)
dim(dorcMat.s)

# This takes longer since it's all genes
colnames(RNAmat.paired) <- colnames(ATAC.se.paired)
RNAmat.s <- smoothScoresNN(NNmat = cellkNN,mat = RNAmat.paired, nCores = 2)
dim(RNAmat.s)

# 11. Visualize on pre-computed UMAP --------
# This is the ATAC UMAP shown in the paper (based on ATAC LSI)
umap.d <- as.data.frame(colData(ATAC.se.paired)[,c("UMAP1","UMAP2")])

# DORC scores for top DORC(s)
myDORCs <- c("IL7R", "TCF7", "CD83")
dorcGGlist <- lapply(myDORCs,function(x) { 
  plotMarker2D(umap.d,
               dorcMat.s,
               markers = x,
               maxCutoff = "q0.99",
               colorPalette = "brewer_heat"
  ) + ggtitle(paste0(x," DORC"))
})

# Paired RNA expression for top DORC(s)
# Plot on the same reference ATAC UMAP
rnaGGlist <- lapply(myDORCs,function(x) { 
  plotMarker2D(umap.d,
               RNAmat.s,
               markers = x,
               maxCutoff = "q0.99",
               colorPalette = "brewer_purple"
  ) + ggtitle(paste0(x," RNA"))
})


library(patchwork)
(dorcGGlist[[1]] + dorcGGlist[[2]] + dorcGGlist[[3]]) /  (rnaGGlist[[1]] + rnaGGlist[[2]] + rnaGGlist[[3]])

# 12. TF-gene associations ---------
# Determine TF-gene associations and inferring a regulatory network based on DORCs
# To determine TFs that are putative regulators (activators or repressors) of DORCs, a built-in reference motif database is used 

# Do not run this code
# figR.d <- runFigRGRN(ATAC.se = ATAC.se.paired, # Must be the same input as used in runGenePeakcorr()
#                      dorcTab = cisCorr.filt, # Filtered peak-gene associations
#                      genome = "hg19",
#                      dorcMat = dorcMat.s,
#                      dorcK = 5, 
#                      rnaMat = RNAmat.s, 
#                      nCores = 1)
figR.d <- readRDS("data/FigR/control1h_PBMC_figR.rds")

figR.d %>% dplyr::arrange(desc(Score)) %>% head()

# 12. Visualizing FigR results ------
rankDrivers(figR.d,rankBy = "meanScore")

# SPI1 gene
SPI1 <-figR.d %>% dplyr::filter(Motif == "SPI1") %>%
  dplyr::arrange(desc(Score)) %>%
  dplyr::filter(Score !=0)
SPI1$DORC <- factor(SPI1$DORC, levels = SPI1$DORC)
SPI1 <- SPI1 %>% dplyr::mutate(color = case_when(Score > 0 ~ "pos",
                                                 Score < 0 ~ "neg"))

SPI1$color <- factor(SPI1$color, c("pos", "zero", "neg"))
SPI1 %>% ggplot(aes(x= DORC, y=Score, fill = color)) +
  geom_bar (stat="identity",position = position_dodge(0.9)) +
  scale_y_continuous(breaks = seq(-2,5,0.2)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
