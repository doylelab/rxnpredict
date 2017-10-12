library(gplots)

# read in symmetry labels
m1.sym <- read.csv("C:\\Users\\Derek\\Dropbox\\rxnpredict\\R\\m1_sym.csv", header=FALSE)
m1.sym.txt <- as.character(unlist(m1.sym[1,]))
m2.sym <- read.csv("C:\\Users\\Derek\\Dropbox\\rxnpredict\\R\\m2_sym.csv", header=FALSE, sep=",", quote="")
m2.sym.txt <- as.character(unlist(m2.sym[1,]))

# read in frequency labels
m1.freq <- read.csv("C:\\Users\\Derek\\Dropbox\\rxnpredict\\R\\m1_freq.csv", header=FALSE)
m1.freq.txt <- m1.freq[m1.freq>500]
m1.freq.txt <- format(round(as.matrix(m1.freq.txt), 0), nsmall = 0)
m1.freq.txt <- as.character(unlist(m1.freq.txt))
m2.freq <- read.csv("C:\\Users\\Derek\\Dropbox\\rxnpredict\\R\\m2_freq.csv", header=FALSE)
m2.freq.txt <- m2.freq[m2.freq>500]
m2.freq.txt <- format(round(as.matrix(m2.freq.txt), 0), nsmall = 0)
m2.freq.txt <- as.character(unlist(m2.freq.txt))
m3.freq <- read.csv("C:\\Users\\Derek\\Dropbox\\rxnpredict\\R\\m3_freq.csv", header=FALSE)
m3.freq.txt <- m3.freq[m3.freq>500]
m3.freq.txt <- format(round(as.matrix(m3.freq.txt), 0), nsmall = 0)
m3.freq.txt <- as.character(unlist(m3.freq.txt))

# read in intensity labels
m1.int <- read.csv("C:\\Users\\Derek\\Dropbox\\rxnpredict\\R\\m1_int.csv", header=FALSE)
m1.int.txt <- format(round(as.matrix(m1.int), 0), nsmall = 0)
m1.int.txt <- as.character(unlist(m1.int.txt[1,]))
m2.int <- read.csv("C:\\Users\\Derek\\Dropbox\\rxnpredict\\R\\m2_int.csv", header=FALSE)
m2.int.txt <- format(round(as.matrix(m2.int), 0), nsmall = 0)
m2.int.txt <- as.character(unlist(m2.int.txt[1,]))

# read in correlation matrices (R2)
corr1 <- read.csv("C:\\Users\\Derek\\Dropbox\\rxnpredict\\R\\corr1.csv", header=FALSE)
corr1 <- corr1[m1.freq>500, m2.freq>500]
corr1 <- as.matrix(corr1)
corr1.txt <- format(round(corr1, 2), nsmall = 2)
corr2 <- read.csv("C:\\Users\\Derek\\Dropbox\\rxnpredict\\R\\corr2.csv", header=FALSE)
corr2 <- corr2[m1.freq>500, m3.freq>500]
corr2 <- as.matrix(corr2)
corr2.txt <- format(round(corr2, 2), nsmall = 2)


colfunc <- colorRampPalette(c("white", "red"))

# corr1 heatmap
png(filename="R\\heatmap_corr1.png", width = 1600, height = 1000)
heatmap.2(corr1, trace="none", dendrogram="none", key=FALSE,
          cellnote=corr1.txt, notecol="black", srtCol=0, 
          notecex=1.3, cexRow=1.3, cexCol=1.3,
          lwid=c(0.1, 4), lhei=c(0.2, 4),
          Rowv=NA, Colv=NA, col=colfunc, 
          labRow=m1.freq.txt, labCol=m2.freq.txt)
dev.off()

# corr2 heatmap
png(filename="R\\heatmap_corr2.png", width = 1600, height = 1000)
heatmap.2(corr2, trace="none", dendrogram="none", key=FALSE,
          cellnote=corr2.txt, notecol="black", srtCol=0,
          notecex=1.3, cexRow=1.3, cexCol=1.3,
          lwid=c(0.1, 4), lhei=c(0.2, 4),
          Rowv=NA, Colv=NA, col=colfunc,
          labRow=m1.freq.txt, labCol=m3.freq.txt)
dev.off()



# PRODUCT YIELD HEATMAP
plate1.yield <- read.csv("C:\\Users\\Derek\\Dropbox\\rxnpredict\\merck_data\\DTA_Plate1_1.csv", header=TRUE)
View(plate1.yield)
yield <- plate1.yield$product_corrected




# read in frequency labels
AH.corr1 <- read.csv("C:\\Users\\Derek\\Dropbox\\rxnpredict\\R\\aryl_halide_corr1_labels.csv", header=TRUE)
AH.corr1.col <- AH.corr1[, 1]
AH.corr1.row <- AH.corr1[, 2]
col <- format(round(as.matrix(AH.corr1.col), 0), nsmall = 0)

png(filename="R\\heatmap_corr1.png", width = 1600, height = 1000)
heatmap.2(corr1, trace="none", dendrogram="none", key=FALSE,
          cellnote=corr1.txt, notecol="black", srtCol=0, 
          notecex=1.3, cexRow=1.3, cexCol=1.3,
          lwid=c(0.1, 4), lhei=c(0.2, 4),
          Rowv=NA, Colv=NA, col=colfunc, 
          labRow=m1.freq.txt, labCol=m2.freq.txt)
dev.off()


m1.freq.txt <- as.character(unlist(m1.freq.txt))

m1.freq.txt <- m1.freq[m1.freq>500]
m1.freq.txt <- format(round(as.matrix(m1.freq.txt), 0), nsmall = 0)
m1.freq.txt <- as.character(unlist(m1.freq.txt))

#####################################################













library(gplots)



colfunc <- colorRampPalette(c("white", "cornflowerblue"))

makeHeatmap <- function(corr, cellnote, m1.labels, m2.labels, filename){
    png(filename=filename, width = 1600, height = 1000)
    heatmap.2(as.matrix(corr), trace="none", dendrogram="none", key=FALSE,
              cellnote=cellnote, notecol="black", srtCol=0, 
              notecex=1.3, cexRow=1.3, cexCol=1.3,
              lwid=c(0.1, 4), lhei=c(0.2, 4),
              Rowv=NA, Colv=NA, col=colfunc, 
              labRow=m1.labels, labCol=m2.labels)
    dev.off()
}

makeHeatmapFile <- function(m1_labels_file, m2_labels_file, corr_file, output_file){
    # read in frequency labels for molecule 1
    AH.corr.m1 <- read.csv(m1_labels_file, header=FALSE)
    AH.corr.m1 <- as.numeric(AH.corr.m1)
    AH.corr.m1.lab <- AH.corr.m1[AH.corr.m1>500]
    AH.corr.m1.lab <- format(round(as.matrix(AH.corr.m1.lab), 0), nsmall = 0)
    
    # read in frequency labels for molecule 2
    AH.corr.m2 <- read.csv(m2_labels_file, header=FALSE)
    AH.corr.m2 <- as.numeric(AH.corr.m2)
    AH.corr.m2.lab <- AH.corr.m2[AH.corr.m2>500]
    AH.corr.m2.lab <- format(round(as.matrix(AH.corr.m2.lab), 0), nsmall = 0)
    
    # read in correlation matrices (R2)
    AH.corr <- read.csv(corr_file, header=FALSE)
    AH.corr <- as.data.frame(AH.corr)
    AH.corr <- AH.corr[AH.corr.m1>500, AH.corr.m2>500]
    AH.corr <- as.matrix(AH.corr)
    AH.corr.cellnote <- format(round(AH.corr, 2), nsmall = 2)
    
    # generate heatmap for aryl_halide_corr
    makeHeatmap(AH.corr, AH.corr.cellnote, AH.corr.m1.lab, AH.corr.m2.lab, output_file)
}

# generate heatmaps for aryl halides (15 of them)
for (i in 1:14){
    m1_labels_file <- paste0("R\\aryl_halide_corr", i, "_m1_labels.csv")
    m2_labels_file <- paste0("R\\aryl_halide_corr", i, "_m2_labels.csv")
    corr_file <- paste0("R\\aryl_halide_corr", i, ".csv")
    output_file <- paste0("R\\aryl_halide_corr", i, "_heatmap.png")
    
    makeHeatmapFile(m1_labels_file,
                    m2_labels_file,
                    corr_file,
                    output_file)
}

# generate heatmaps for additives (22 of them)
for (i in 1:21){
    m1_labels_file <- paste0("R\\additive_corr", i, "_m1_labels.csv")
    m2_labels_file <- paste0("R\\additive_corr", i, "_m2_labels.csv")
    corr_file <- paste0("R\\additive_corr", i, ".csv")
    output_file <- paste0("R\\additive_corr", i, "_heatmap.png")
    
    makeHeatmapFile(m1_labels_file,
                    m2_labels_file,
                    corr_file,
                    output_file)
}











# ============

# changed the size (larger) and removed cell labels for the ligands
makeHeatmapLarge <- function(corr, cellnote, m1.labels, m2.labels, filename){
    png(filename=filename, width = 6000, height = 4000)
    heatmap.2(as.matrix(corr), trace="none", dendrogram="none", key=FALSE,
              notecol="black", srtCol=0, # cellnote=cellnote, 
              notecex=1.3, cexRow=1.3, cexCol=1.3,
              lwid=c(0.1, 4), lhei=c(0.2, 4),
              Rowv=NA, Colv=NA, col=colfunc, 
              labRow=m1.labels, labCol=m2.labels)
    dev.off()
}

makeHeatmapFileLarge <- function(m1_labels_file, m2_labels_file, corr_file, output_file){
    # read in frequency labels for molecule 1
    AH.corr.m1 <- read.csv(m1_labels_file, header=FALSE)
    AH.corr.m1 <- as.numeric(AH.corr.m1)
    AH.corr.m1.lab <- AH.corr.m1[AH.corr.m1>500]
    AH.corr.m1.lab <- format(round(as.matrix(AH.corr.m1.lab), 0), nsmall = 0)
    
    # read in frequency labels for molecule 2
    AH.corr.m2 <- read.csv(m2_labels_file, header=FALSE)
    AH.corr.m2 <- as.numeric(AH.corr.m2)
    AH.corr.m2.lab <- AH.corr.m2[AH.corr.m2>500]
    AH.corr.m2.lab <- format(round(as.matrix(AH.corr.m2.lab), 0), nsmall = 0)
    
    # read in correlation matrices (R2)
    AH.corr <- read.csv(corr_file, header=FALSE)
    AH.corr <- as.data.frame(AH.corr)
    AH.corr <- AH.corr[AH.corr.m1>500, AH.corr.m2>500]
    AH.corr <- as.matrix(AH.corr)
    AH.corr.cellnote <- format(round(AH.corr, 2), nsmall = 2)
    
    # generate heatmap for aryl_halide_corr
    makeHeatmapLarge(AH.corr, AH.corr.cellnote, AH.corr.m1.lab, AH.corr.m2.lab, output_file)
}

# generate heatmaps for ligands (4 of them)
for (i in 1:3){
    m1_labels_file <- paste0("R\\ligand_corr", i, "_m1_labels.csv")
    m2_labels_file <- paste0("R\\ligand_corr", i, "_m2_labels.csv")
    corr_file <- paste0("R\\ligand_corr", i, ".csv")
    output_file <- paste0("R\\ligand_corr", i, "_heatmap.png")
    
    makeHeatmapFileLarge(m1_labels_file,
                    m2_labels_file,
                    corr_file,
                    output_file)
}















resamps <- resamples(list(lm = lmFit,
                          lm.reduced = lmFit.reduced, 
                          knn = knnFit,
                          svm = svmFit,
                          rf = rfFit))
summary(resamps)


resamps2 <- resamples(list(lm = lmFit))
bwplot(resamps2,metric="ROC",main="GBM vs xgboost")

png(filename="C:\\Users\\Derek\\Dropbox\\rxnpredict\\R\\plots\\model_comparison.png", width = 600, height = 300)
bwplot(resamps,
       layout = c(2, 1),
       scales = list(relation = "free"),
       xlim = list(c(0, 20), c(0, 1)))
dev.off()


