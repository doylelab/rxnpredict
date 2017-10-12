# install packages (if necessary) and load them
if (!require("pacman")) install.packages("pacman")
pacman::p_load(gplots, 
               ggplot2, 
               GGally,
               caret,
               ModelMetrics,
               glmnet,
               gridExtra,
               randomForest,
               scales,
               arm,
               corrplot)


# ============================================================================
# Helper functions for making heatmaps and histograms
# ============================================================================

# heatmaps
makeHeatmap1536 <- function(data, labels, filename){
    color1 <- "white"
    color2 <- "cornflowerblue"
    color3 <- "green"
    colfunc <- colorRampPalette(c(color1, color2, color3))
    
    png(filename=filename, width = 1600, height = 1000)
    heatmap.2(data, trace="none", dendrogram="none", key=FALSE,
              cellnote=labels, notecol="black", srtCol=0, 
              notecex=2, cexRow=1.5, cexCol=1.5,
              lwid=c(0.1, 4), lhei=c(0.2, 4),
              Rowv=NA, Colv=NA, col=colfunc,
              labRow=c(1:32), labCol=c(1:48))
    dev.off()
}

# histograms
ggHist <- function(data, filename){
    data <- as.data.frame(unlist(data))
    colnames(data) <- c("yield")
    ggplot(data, aes(x=yield)) +
        geom_histogram(colour = "black", 
                       fill = "cornflowerblue", 
                       breaks=seq(0,100,by=5), 
                       na.rm = TRUE) +
        labs(x="Yield", y="Count") + 
    ggsave(file=filename, width=5, height=3)
}



# ============================================================================
# Set R's working directory
# ============================================================================

# set the working directory to the location of the rxnpredict folder
setwd("C:\\Users\\Derek\\Dropbox\\rxnpredict")



# ============================================================================
# Load yield data and stitch it together 
# ============================================================================

# In my case, there were 3 1536-well plates and the data for each plate
# was analyzed by quadrant (via 4 384-well plates)

# Plate 1.1
plate1.1 <- read.csv("yield_data\\plate1.1.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate1.1_pdt <- plate1.1$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate1.1_pdt))
dim(plate.data) <- c(24,16)
plate.data1.1 <- t(plate.data)

# Plate 1.2
plate1.2 <- read.csv("yield_data\\plate1.2.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate1.2_pdt <- plate1.2$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate1.2_pdt))
dim(plate.data) <- c(24,16)
plate.data1.2 <- t(plate.data)

# Plate 1.3
plate1.3 <- read.csv("yield_data\\plate1.3.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate1.3_pdt <- plate1.3$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate1.3_pdt))
dim(plate.data) <- c(24,16)
plate.data1.3 <- t(plate.data)

# Plate 1.4
plate1.4 <- read.csv("yield_data\\plate1.4.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate1.4_pdt <- plate1.4$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate1.4_pdt))
dim(plate.data) <- c(24,16)
plate.data1.4 <- t(plate.data)

# stitch Plate 1 together into one 32x48 matrix
plate1.top <- cbind(plate.data1.1, plate.data1.2)
plate1.bottom <- cbind(plate.data1.3, plate.data1.4)
plate1 <- rbind(plate1.top, plate1.bottom)

# Plate 2.1
plate2.1 <- read.csv("yield_data\\plate2.1.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate2.1_pdt <- plate2.1$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate2.1_pdt))
dim(plate.data) <- c(24,16)
plate.data2.1 <- t(plate.data)

# Plate 2.2
plate2.2 <- read.csv("yield_data\\plate2.2.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate2.2_pdt <- plate2.2$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate2.2_pdt))
dim(plate.data) <- c(24,16)
plate.data2.2 <- t(plate.data)

# Plate 2.3
plate2.3 <- read.csv("yield_data\\plate2.3.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate2.3_pdt <- plate2.3$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate2.3_pdt))
dim(plate.data) <- c(24,16)
plate.data2.3 <- t(plate.data)

# Plate 2.4
plate2.4 <- read.csv("yield_data\\plate2.4.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate2.4_pdt <- plate2.4$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate2.4_pdt))
dim(plate.data) <- c(24,16)
plate.data2.4 <- t(plate.data)

# stitch Plate 2 together into one 32x48 matrix
plate2.top <- cbind(plate.data2.1, plate.data2.2)
plate2.bottom <- cbind(plate.data2.3, plate.data2.4)
plate2 <- rbind(plate2.top, plate2.bottom)

# Plate 3.1
plate3.1 <- read.csv("yield_data\\plate3.1.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate3.1_pdt <- plate3.1$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate3.1_pdt))
dim(plate.data) <- c(24,16)
plate.data3.1 <- t(plate.data)

# Plate 3.2
plate3.2 <- read.csv("yield_data\\plate3.2.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate3.2_pdt <- plate3.2$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate3.2_pdt))
dim(plate.data) <- c(24,16)
plate.data3.2 <- t(plate.data)

# Plate 3.3
plate3.3 <- read.csv("yield_data\\plate3.3.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate3.3_pdt <- plate3.3$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate3.3_pdt))
dim(plate.data) <- c(24,16)
plate.data3.3 <- t(plate.data)

# Plate 3.4
plate3.4 <- read.csv("yield_data\\plate3.4.csv", header=TRUE, stringsAsFactors=FALSE, na.strings = "#DIV/0!")
plate3.4_pdt <- plate3.4$product_scaled[1:384]
plate.data <- as.matrix(as.numeric(plate3.4_pdt))
dim(plate.data) <- c(24,16)
plate.data3.4 <- t(plate.data)

# stitch Plate 3 together into one 32x48 matrix
plate3.top <- cbind(plate.data3.1, plate.data3.2)
plate3.bottom <- cbind(plate.data3.3, plate.data3.4)
plate3 <- rbind(plate3.top, plate3.bottom)



# ============================================================================
# Make heatmaps and histograms 
# ============================================================================

# Plate 1 heatmap and histogram
plate1.txt <- format(round(plate1, 0), nsmall = 0)
makeHeatmap1536(plate1, plate1.txt, "yield_data\\plate1_heatmap.png")
plate1 <- as.data.frame(plate1)
ggHist(plate1, "yield_data\\plate1_histogram.png")

# Plate 2 heatmap and histogram
plate2.txt <- format(round(plate2, 0), nsmall = 0)
makeHeatmap1536(plate2, plate2.txt, "yield_data\\plate2_heatmap.png")
plate2 <- as.data.frame(plate2)
ggHist(plate2, "yield_data\\plate2_histogram.png")

# Plate 3 heatmap and histogram
plate3.txt <- format(round(plate3, 0), nsmall = 0)
makeHeatmap1536(plate3, plate3.txt, "yield_data\\plate3_heatmap.png")
plate3 <- as.data.frame(plate3)
ggHist(plate3, "yield_data\\plate3_histogram.png")

# All plates histogram
top.two <- rbind(plate1, plate2)
all.plates <- rbind(top.two, plate3)
all.plates <- as.data.frame(all.plates)
ggHist(all.plates, "yield_data\\allplates_histogram.png")



# ============================================================================
# Load and prepare output table for modeling
# ============================================================================

# Remove reactions without additive and reactions with additive 7
plate1_nocontrols <- plate1[c(-1,-5,-9,-13,-20,-24,-28,-32), c(-16,-32,-48)] 
# Remove reactions without aryl halide
plate2_nocontrols <- plate2[, c(-16,-32,-48)]
plate3_nocontrols <- plate3[, c(-16,-32,-48)]
plate1_nocontrols_v <- as.vector(t(plate1_nocontrols))
plate2_nocontrols_v <- as.vector(t(plate2_nocontrols))
plate3_nocontrols_v <- as.vector(t(plate3_nocontrols))
yield_data <- c(plate1_nocontrols_v, plate2_nocontrols_v, plate3_nocontrols_v)

# load output table generated by python script
output.table <- read.csv("R\\output_table.csv", header=TRUE)

# scale the descriptor data
output.scaled <- as.data.frame(scale(output.table))

# append the yield data from above
output.scaled$yield <- yield_data

# Uncomment and run line below to view large datasets
# utils::View(output.scaled)



# ============================================================================
# Subset the data to test ability to predict out-of-sample prediction
# e.g., pyridyl from nonpyridyl, over 80% yield rxns from under 80% yield rxns
# ============================================================================

# separate into plates 1,2,3 and remove NA
output.plate1 <- output.scaled[1:1080, ]
output.plate1 <- output.plate1[!(is.na(output.plate1$yield)), ]
output.plate2 <- output.scaled[1081:2520, ]
output.plate2 <- output.plate2[!(is.na(output.plate2$yield)), ]
output.plate3 <- output.scaled[2521:3960, ]
output.plate3 <- output.plate3[!(is.na(output.plate3$yield)), ]

# separate by all 15 aryl halides
CF3.Cl <- seq(1, nrow(output.scaled), by=15)
CF3.Br <- seq(2, nrow(output.scaled), by=15)
CF3.I <- seq(3, nrow(output.scaled), by=15)
OMe.Cl <- seq(4, nrow(output.scaled), by=15)
OMe.Br <- seq(5, nrow(output.scaled), by=15)
OMe.I <- seq(6, nrow(output.scaled), by=15)
Et.Cl <- seq(7, nrow(output.scaled), by=15)
Et.Br <- seq(8, nrow(output.scaled), by=15)
Et.I <- seq(9, nrow(output.scaled), by=15)
pyr2.Cl <- seq(10, nrow(output.scaled), by=15)
pyr2.Br <- seq(11, nrow(output.scaled), by=15)
pyr2.I <- seq(12, nrow(output.scaled), by=15)
pyr3.Cl <- seq(13, nrow(output.scaled), by=15)
pyr3.Br <- seq(14, nrow(output.scaled), by=15)
pyr3.I <- seq(15, nrow(output.scaled), by=15)

# nonpyridyl aryl halides
nonpyridyl <- sort(c(CF3.Cl, CF3.Br, CF3.I, OMe.Cl, OMe.Br, OMe.I, Et.Cl, Et.Br, Et.I))
nonpyridyl.scaled <- output.scaled[nonpyridyl, ]
nonpyridyl.scaled <- nonpyridyl.scaled[!(is.na(nonpyridyl.scaled$yield)), ]

# pyridyl aryl halides
pyridyl <- sort(c(pyr2.Cl, pyr2.Br, pyr2.I, pyr3.Cl, pyr3.Br, pyr3.I))
pyridyl.scaled <- output.scaled[pyridyl, ]
pyridyl.scaled <- pyridyl.scaled[!(is.na(pyridyl.scaled$yield)), ]

# separate by chloride, bromide, iodide
ArCl <- seq(1, nrow(output.scaled), by=3)
ArCl.scaled <- output.scaled[ArCl, ]
ArCl.scaled <- ArCl.scaled[!(is.na(ArCl.scaled$yield)), ]
set.seed(8390)
size <- round(0.70*nrow(ArCl.scaled))
training <- sample(nrow(ArCl.scaled), size=size, replace=FALSE)
ArCl.training <- ArCl.scaled[training, ]
ArCl.test <- ArCl.scaled[-training, ]

ArBr <- seq(2, nrow(output.scaled), by=3)
ArBr.scaled <- output.scaled[ArBr, ]
ArBr.scaled <- ArBr.scaled[!(is.na(ArBr.scaled$yield)), ]
set.seed(9071)
size <- round(0.70*nrow(ArBr.scaled))
training <- sample(nrow(ArBr.scaled), size=size, replace=FALSE)
ArBr.training <- ArBr.scaled[training, ]
ArBr.test <- ArBr.scaled[-training, ]

ArI <- seq(3, nrow(output.scaled), by=3)
ArI.scaled <- output.scaled[ArI, ]
ArI.scaled <- ArI.scaled[!(is.na(ArI.scaled$yield)), ]
set.seed(6123)
size <- round(0.70*nrow(ArI.scaled))
training <- sample(nrow(ArI.scaled), size=size, replace=FALSE)
ArI.training <- ArI.scaled[training, ]
ArI.test <- ArI.scaled[-training, ]

# yields under 80% and over 80%
under80 <- output.scaled$yield<80
output.under80 <- output.scaled[under80, ]
output.over80 <- output.scaled[!under80, ]



# ============================================================================
# Remove reactions without yield data and make histogram
# ============================================================================

# remove rows where yield=NA
output.scaled <- output.scaled[!(is.na(output.scaled$yield)), ]

# Histogram for modeling yields (removed controls and additive 7)
ggHist(output.scaled$yield, "yield_data\\allplates_nocontrols_histogram.png")



# ============================================================================
# Data splitting for modeling
# ============================================================================

# Split into training and test set (70/30)
set.seed(1084)
size <- round(0.70*nrow(output.scaled))
training <- sample(nrow(output.scaled), size=size, replace=FALSE)
training.scaled <- output.scaled[training,]
test.scaled <- output.scaled[-training,]

# Create smaller partitions within training set (equal to 10, 20, etc. % of TOTAL data)
size2.5 <- round(0.025*nrow(output.scaled))
size5 <- round(0.05*nrow(output.scaled))
size10 <- round(0.10*nrow(output.scaled))
size20 <- round(0.20*nrow(output.scaled))
size30 <- round(0.30*nrow(output.scaled))
size40 <- round(0.40*nrow(output.scaled))
size50 <- round(0.50*nrow(output.scaled))
size60 <- round(0.60*nrow(output.scaled))
training2.5rows <- sample(nrow(training.scaled),size=size2.5,replace=FALSE)
training2.5 <- training.scaled[training2.5rows, ]
training5rows <- sample(nrow(training.scaled),size=size5,replace=FALSE)
training5 <- training.scaled[training5rows, ]
training10rows <- sample(nrow(training.scaled),size=size10,replace=FALSE)
training10 <- training.scaled[training10rows, ]
training20rows <- sample(nrow(training.scaled),size=size20,replace=FALSE)
training20 <- training.scaled[training20rows, ]
training30rows <- sample(nrow(training.scaled),size=size30,replace=FALSE)
training30 <- training.scaled[training30rows, ]
training40rows <- sample(nrow(training.scaled),size=size40,replace=FALSE)
training40 <- training.scaled[training40rows, ]
training50rows <- sample(nrow(training.scaled),size=size50,replace=FALSE)
training50 <- training.scaled[training50rows, ]
training60rows <- sample(nrow(training.scaled),size=size60,replace=FALSE)
training60 <- training.scaled[training60rows, ]




# ============================================================================
# Model Training
# Only need to run once - models are saved as .rds files
# ============================================================================

# 10-fold cross-validation
train_control <- trainControl(method="cv", number=10, savePredictions=TRUE)

# k-nearest neighbor (kNN)
knnFit <- train(yield ~ ., data=training.scaled, trControl=train_control, method="knn")
png(filename="R\\plots\\knn.png", width = 1000, height = 600)
predVals <- extractPrediction(list(knnFit))
plotObsVsPred(predVals)
dev.off()
saveRDS(knnFit, "rds\\knnFit.rds")

# support vector machine (SVM)
svmFit <- train(yield ~ ., data=training.scaled, trControl=train_control, method="svmLinear")
png(filename="R\\plots\\svm.png", width = 1000, height = 600)
predVals <- extractPrediction(list(svmFit))
plotObsVsPred(predVals)
dev.off()
saveRDS(svmFit, "rds\\svmFit.rds")

# Bayes generalized linear model (GLM)
bayesglmFit <- train(yield ~ ., data=training.scaled, trControl=train_control, method="bayesglm")
png(filename="R\\plots\\bayesglm.png", width = 1000, height = 600)
predVals <- extractPrediction(list(bayesglmFit))
plotObsVsPred(predVals)
dev.off()
saveRDS(bayesglmFit, "rds\\bayesglmFit.rds")

# linear model
lmFit <- train(yield ~ ., data=training.scaled, trControl=train_control, method="lm")
png(filename="R\\plots\\lm.png", width = 1000, height = 600)
predVals <- extractPrediction(list(lmFit))
plotObsVsPred(predVals)
dev.off()
saveRDS(lmFit, "rds\\lmFit.rds")

# neural network
set.seed(3288)
nnetFit <- train(yield ~ ., data=training.scaled, trControl=train_control, method="nnet", linout=1)
png(filename="R\\plots\\nnet.png", width = 1000, height = 600)
predVals <- extractPrediction(list(nnetFit))
plotObsVsPred(predVals)
dev.off()
saveRDS(nnetFit, "rds\\nnetFit.rds")

# random forest 
rfFit <- train(yield ~ ., data=training.scaled, trControl=train_control, method="rf", importance=TRUE)
png(filename="R\\plots\\rf70.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit70))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit70, "rds\\rfFit70.rds")










####################################################
# RANDOM FOREST (by nonpyridyl/pyridyl)
####################################################

# Random forest (train on nonpyridyl)
rfFit.nonpyridyl <- train(yield ~ ., data=nonpyridyl.scaled, trControl=train_control, method="rf", importance=TRUE)
png(filename="R\\plots\\rf_nonpyridyl.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit.nonpyridyl))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit.nonpyridyl, "rds\\rfFit_nonpyridyl.rds")

# Test on pyridyl
pyridyl.pred <- predict(rfFit.nonpyridyl, pyridyl.scaled)
pyridyl.r2 <- cor(pyridyl.pred, pyridyl.scaled$yield)
pyridyl.rmse <- rmse(pyridyl.pred, pyridyl.scaled$yield)

df1 <- data.frame(x = pyridyl.pred, 
                  y = pyridyl.scaled$yield)
p1 <- ggplot(df1, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\pyr_from_nonpyr.png", width=5, height=4)











####################################################
# RANDOM FOREST (fit CF3/OMe/2-pyr, pred Et/3-pyr)
####################################################


CF3.OMe.pyr2 <- sort(c(CF3.Cl, CF3.Br, CF3.I, OMe.Cl, OMe.Br, OMe.I, pyr2.Cl, pyr2.Br, pyr2.I))
CF3.OMe.pyr2.scaled <- output.scaled[CF3.OMe.pyr2, ]
CF3.OMe.pyr2.scaled <- CF3.OMe.pyr2.scaled[!(is.na(CF3.OMe.pyr2.scaled$yield)), ]

Et.pyr3 <- sort(c(Et.Cl, Et.Br, Et.I, pyr3.Cl, pyr3.Br, pyr3.I))
Et.pyr3.scaled <- output.scaled[Et.pyr3, ]
Et.pyr3.scaled <- Et.pyr3.scaled[!(is.na(Et.pyr3.scaled$yield)), ]

Et <- sort(c(Et.Cl, Et.Br, Et.I))
Et.scaled <- output.scaled[Et, ]
Et.scaled <- Et.scaled[!(is.na(Et.scaled$yield)), ]

pyr3 <- sort(c(pyr3.Cl, pyr3.Br, pyr3.I))
pyr3.scaled <- output.scaled[pyr3, ]
pyr3.scaled <- pyr3.scaled[!(is.na(pyr3.scaled$yield)), ]



# Random forest (train on nonpyridyl)
rfFit.CF3.OMe.pyr2 <- train(yield ~ ., data=CF3.OMe.pyr2.scaled, trControl=train_control, method="rf", importance=TRUE)
png(filename="R\\plots\\rf_CF3-OMe-pyr2.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit.CF3.OMe.pyr2))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit.CF3.OMe.pyr2, "rds\\rfFit_CF3-OMe-pyr2.rds")

# Test on Et and 3-pyr
Et.pyr3.pred <- predict(rfFit.CF3.OMe.pyr2, Et.pyr3.scaled)
Et.pyr3.r2 <- cor(Et.pyr3.pred, Et.pyr3.scaled$yield)
Et.pyr3.rmse <- rmse(Et.pyr3.pred, Et.pyr3.scaled$yield)

df1 <- data.frame(x = Et.pyr3.pred, 
                  y = Et.pyr3.scaled$yield)
p1 <- ggplot(df1, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\Et-pyr3-from-CF3-OMe-pyr2.png", width=5, height=4)

# Test on Et only
Et.pred <- predict(rfFit.CF3.OMe.pyr2, Et.scaled)
Et.r2 <- cor(Et.pred, Et.scaled$yield)
Et.rmse <- rmse(Et.pred, Et.scaled$yield)

df1 <- data.frame(x = Et.pred, 
                  y = Et.scaled$yield)
p1 <- ggplot(df1, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\Et-from-CF3-OMe-pyr2.png", width=5, height=4)

# Test on 3-pyr only
pyr3.pred <- predict(rfFit.CF3.OMe.pyr2, pyr3.scaled)
pyr3.r2 <- cor(pyr3.pred, pyr3.scaled$yield)
pyr3.rmse <- rmse(pyr3.pred, pyr3.scaled$yield)

df1 <- data.frame(x = pyr3.pred, 
                  y = pyr3.scaled$yield)
p1 <- ggplot(df1, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\pyr3-from-CF3-OMe-pyr2.png", width=5, height=4)






####################################################
# RANDOM FOREST (by aryl halide)
####################################################

# Train and test each aryl halide individually

# Random forest (ArCl)
rfFit.ArCl <- train(yield ~ ., data=ArCl.training, trControl=train_control, method="rf", importance=TRUE)
png(filename="R\\plots\\rf_ArCl.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit.ArCl))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit.ArCl, "rds\\rfFit_ArCl.rds")

# Random forest (ArBr)
rfFit.ArBr <- train(yield ~ ., data=ArBr.training, trControl=train_control, method="rf", importance=TRUE)
png(filename="R\\plots\\rf_ArBr.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit.ArBr))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit.ArBr, "rds\\rfFit_ArBr.rds")

# Random forest (ArI)
rfFit.ArI <- train(yield ~ ., data=ArI.training, trControl=train_control, method="rf", importance=TRUE)
png(filename="R\\plots\\rf_ArI.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit.ArI))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit.ArI, "rds\\rfFit_ArI.rds")

# Train on ArBr and test on ArCl and ArI

# Random forest (ArBr)
rfFit.ArBr.all <- train(yield ~ ., data=ArBr.scaled, trControl=train_control, method="rf", importance=TRUE)
saveRDS(rfFit.ArBr.all, "rds\\rfFit_ArCl_all.rds")

# Test on ArCl
ClfromBr <- predict(rfFit.ArBr.all, ArCl.scaled)
ClfromBr.r2 <- cor(ClfromBr, ArCl.scaled$yield)
ClfromBr.rmse <- rmse(ClfromBr, ArCl.scaled$yield)
df1 <- data.frame(x = ClfromBr, 
                 y = ArCl.scaled$yield)

p1 <- ggplot(df1, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\ClfromBr.png", width=5, height=4)

# Test on ArI
IfromBr <- predict(rfFit.ArBr.all, ArI.scaled)
IfromBr.r2 <- cor(IfromBr, ArI.scaled$yield)
IfromBr.rmse <- rmse(IfromBr, ArI.scaled$yield)
df2 <- data.frame(x = IfromBr, 
                  y = ArI.scaled$yield)

# plot I from Br
p2 <- ggplot(df2, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\IfromBr.png", width=5, height=4)







####################################################
# LINEAR PLOT OF *C3 NMR SHIFT VS. YIELD
####################################################

p1 <- ggplot(output.scaled, aes(x = additive_.C3_NMR_shift, y = yield)) +
    geom_point(alpha = 0.2, size=1) + 
    # scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    # labs(x='Predicted Yield', y='Observed Yield') +  
    # geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\C3_vs_yield_loess.png", width=8, height=6)

p2 <- ggplot(output.scaled, aes(x = additive_.C3_NMR_shift, y = yield)) +
    geom_point(alpha = 0.2, size=1) + 
    # scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    # labs(x='Predicted Yield', y='Observed Yield') +  
    # geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="lm", se=FALSE)
ggsave(file="R\\plots\\C3_vs_yield_lm.png", width=8, height=6)

C3.yield.cor <- cor(output.scaled$additive_.C3_electrostatic_charge, output.scaled$yield)





####################################################
# RANDOM FOREST (yields < 80%) COMPLETE ANALYSIS
####################################################

# Random forest (yields < 80%)
rfFit.under80 <- train(yield ~ ., data=output.under80, trControl=train_control, method="rf", importance=TRUE)
png(filename="R\\plots\\rf_under80.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit.under80))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit.under80, "rds\\rfFit_under80.rds")

over80pred <- predict(rfFit.under80, output.over80)
over80pred.r2 <- cor(over80pred, output.over80$yield)
over80pred.rmse <- rmse(over80pred, output.over80$yield)

df <- data.frame(x = over80pred, 
                 y = output.over80$yield)

p1 <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)

ggsave(file="R\\plots\\over80pred.png", width=5, height=4)



####################################################
# RANDOM FOREST (with sparsity testing)
####################################################

# Random forest (2.5% of data)
rfFit2.5 <- train(yield ~ ., data=training2.5, trControl=train_control, method="rf")
png(filename="R\\plots\\rf2_5.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit2.5))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit2.5, "rds\\rfFit2_5.rds")

# Random forest (5% of data)
rfFit5 <- train(yield ~ ., data=training5, trControl=train_control, method="rf")
png(filename="R\\plots\\rf5.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit5))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit5, "rds\\rfFit5.rds")

# Random forest (10% of data)
rfFit10 <- train(yield ~ ., data=training10, trControl=train_control, method="rf")
png(filename="R\\plots\\rf10.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit10))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit10, "rds\\rfFit10.rds")

# Random forest (20% of data)
rfFit20 <- train(yield ~ ., data=training20, trControl=train_control, method="rf")
png(filename="R\\plots\\rf20.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit20))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit20, "rds\\rfFit20.rds")

# Random forest (30% of data)
rfFit30 <- train(yield ~ ., data=training30, trControl=train_control, method="rf")
png(filename="R\\plots\\rf30.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit30))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit30, "rds\\rfFit30.rds")

# Random forest (40% of data)
rfFit40 <- train(yield ~ ., data=training40, trControl=train_control, method="rf")
png(filename="R\\plots\\rf40.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit40))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit40, "rds\\rfFit40.rds")

# Random forest (50% of data)
rfFit50 <- train(yield ~ ., data=training50, trControl=train_control, method="rf")
png(filename="R\\plots\\rf50.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit50))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit50, "rds\\rfFit50.rds")

# Random forest (60% of data)
rfFit60 <- train(yield ~ ., data=training60, trControl=train_control, method="rf")
png(filename="R\\plots\\rf60.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit60))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit60, "rds\\rfFit60.rds")

# Random forest (70% of data - all training data)  [rfFit70 is identical to rfFit]
rfFit70 <- train(yield ~ ., data=training.scaled, trControl=train_control, method="rf", importance=TRUE)
png(filename="R\\plots\\rf70.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit70))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit70, "rds\\rfFit70.rds")





















#######################################################
# TEST SET PREDICTION AND PLOTTING (rf, xgb, other)
#######################################################

# RANDOM FOREST (ARYL HALIDES) ========================

# ArCl ------------------------------------------------
rf.pred.ArCl <- predict(rfFit.ArCl, ArCl.test)
plot(rf.pred.ArCl, ArCl.test$yield)
rf.pred.ArCl.r2 <- cor(rf.pred.ArCl, ArCl.test$yield)
rf.pred.ArCl.rmse <- rmse(rf.pred.ArCl, ArCl.test$yield)

# Variable importance
png(filename="R\\plots\\ArCl_importance.png", width = 1000, height = 600)
varImpPlot(rfFit.ArCl$finalModel, type=2, main="")
dev.off()


# ArBr ------------------------------------------------
rf.pred.ArBr <- predict(rfFit.ArBr, ArBr.test)
plot(rf.pred.ArBr, ArBr.test$yield)
rf.pred.ArBr.r2 <- cor(rf.pred.ArBr, ArBr.test$yield)
rf.pred.ArBr.rmse <- rmse(rf.pred.ArBr, ArBr.test$yield)

# Variable importance
png(filename="R\\plots\\ArBr_importance.png", width = 1000, height = 600)
varImpPlot(rfFit.ArBr$finalModel, type=2, main="")
dev.off()


# ArI  ------------------------------------------------
rf.pred.ArI <- predict(rfFit.ArI, ArI.test)
plot(rf.pred.ArI, ArI.test$yield)
rf.pred.ArI.r2 <- cor(rf.pred.ArI, ArI.test$yield)
rf.pred.ArI.rmse <- rmse(rf.pred.ArI, ArI.test$yield)

# Variable importance
png(filename="R\\plots\\ArI_importance.png", width = 1000, height = 600)
varImpPlot(rfFit.ArI$finalModel, type=2, main="")
dev.off()


# Pretty plot RMSE and Rsquared
df <- data.frame(rmse = c(rf.pred.ArCl.rmse, rf.pred.ArBr.rmse, rf.pred.ArI.rmse),
                 r2 = c(rf.pred.ArCl.r2, rf.pred.ArBr.r2, rf.pred.ArCl.r2))
row.names(df) <- c('Random Forest (ArCl)', 'Random Forest (ArBr)', 'Random Forest (ArI)')

rmse.plot <- ggplot(df, aes(y=rownames(df), x=rmse)) +
    geom_point() +
    geom_text(label=round(df$rmse, 2), vjust=-1, size=3) +
    labs(x='RMSE', y='') +
    xlim(0,20)

r2.plot <- ggplot(df, aes(y=rownames(df), x=r2)) +
    geom_point() +
    geom_text(label=round(df$r2, 2), vjust=-1, size=3) +
    labs(x='Rsquared', y='') +
    xlim(0.7,1)

plots <- arrangeGrob(rmse.plot, r2.plot, ncol=2)
ggsave(plots, file="R\\plots\\model_comp_rf_ArX.png", width=8, height=3)








# ============================================================================
# Random forest sparsity testing and variable importance
# ============================================================================

# Predict for testing set
rf.pred2.5 <- predict(rfFit2.5, test.scaled)
rf.pred5 <- predict(rfFit5, test.scaled)
rf.pred10 <- predict(rfFit10, test.scaled)
rf.pred20 <- predict(rfFit20, test.scaled)
rf.pred30 <- predict(rfFit30, test.scaled)
rf.pred40 <- predict(rfFit40, test.scaled)
rf.pred50 <- predict(rfFit50, test.scaled)
rf.pred60 <- predict(rfFit60, test.scaled)
rf.pred70 <- predict(rfFit70, test.scaled)

# Plot predicted vs. observed (test set)
plot(rf.pred2.5, test.scaled$yield)
plot(rf.pred5, test.scaled$yield)
plot(rf.pred10, test.scaled$yield)
plot(rf.pred20, test.scaled$yield)
plot(rf.pred30, test.scaled$yield)
plot(rf.pred40, test.scaled$yield)
plot(rf.pred50, test.scaled$yield)
plot(rf.pred60, test.scaled$yield)
plot(rf.pred70, test.scaled$yield)

# R^2 values
rf.pred2.5.r2 <- cor(rf.pred2.5, test.scaled$yield)
rf.pred5.r2 <- cor(rf.pred5, test.scaled$yield)
rf.pred10.r2 <- cor(rf.pred10, test.scaled$yield)
rf.pred20.r2 <- cor(rf.pred20, test.scaled$yield)
rf.pred30.r2 <- cor(rf.pred30, test.scaled$yield)
rf.pred40.r2 <- cor(rf.pred40, test.scaled$yield)
rf.pred50.r2 <- cor(rf.pred50, test.scaled$yield)
rf.pred60.r2 <- cor(rf.pred60, test.scaled$yield)
rf.pred70.r2 <- cor(rf.pred70, test.scaled$yield)

# RMSE
rf.pred2.5.rmse <- rmse(rf.pred2.5, test.scaled$yield)
rf.pred5.rmse <- rmse(rf.pred5, test.scaled$yield)
rf.pred10.rmse <- rmse(rf.pred10, test.scaled$yield)
rf.pred20.rmse <- rmse(rf.pred20, test.scaled$yield)
rf.pred30.rmse <- rmse(rf.pred30, test.scaled$yield)
rf.pred40.rmse <- rmse(rf.pred40, test.scaled$yield)
rf.pred50.rmse <- rmse(rf.pred50, test.scaled$yield)
rf.pred60.rmse <- rmse(rf.pred60, test.scaled$yield)
rf.pred70.rmse <- rmse(rf.pred70, test.scaled$yield)



rf_imp <- importance(rfFit70$finalModel)
rf.imp.df <- cbind(as.data.frame(rf_imp), names(rf_imp[, 1]))
colnames(rf.imp.df)[1] <- "IncMSE"
colnames(rf.imp.df)[3] <- "descriptor"
rf.imp.df$descriptor <- gsub("_", " ", rf.imp.df$descriptor)
rf.imp.df$descriptor <- gsub("[.]", "*", rf.imp.df$descriptor)

simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
          sep="", collapse=" ")
}

rf.imp.df$descriptor <- sapply(rf.imp.df$descriptor, simpleCap)

# change 15 to 20!!! and probably plot height too
p1 <- ggplot(rf.imp.df[rf.imp.df$IncMSE>15, ], aes(x=reorder(descriptor, IncMSE), y=IncMSE)) +
    geom_bar(stat="identity") +
    scale_y_continuous(labels = comma) +
    labs(x="", y="Increase in Mean Squared Error (%)") + 
    coord_flip()
ggsave(p1, file="R\\plots\\rf_importance.png", width=7, height=5)


# Pretty plot RMSE and Rsquared
df <- data.frame(rmse = c(rf.pred2.5.rmse, rf.pred5.rmse, rf.pred10.rmse, rf.pred20.rmse, rf.pred30.rmse, rf.pred50.rmse, rf.pred70.rmse),
                 r2 = c(rf.pred2.5.r2, rf.pred5.r2, rf.pred10.r2, rf.pred20.r2, rf.pred30.r2, rf.pred50.r2, rf.pred70.r2))
row.names(df) <- c('Random Forest (2.5%)', 'Random Forest (5%)', 'Random Forest (10%)', 'Random Forest (20%)', 
                   'Random Forest (30%)', 'Random Forest (50%)', 'Random Forest (70%)')

rmse.plot <- ggplot(df, aes(y=reorder(rownames(df), rmse), x=rmse)) +
    geom_point() +
    geom_text(label=round(df$rmse, 2), vjust=-1, size=3) +
    labs(x='RMSE', y='') +
    xlim(0,20)

r2.plot <- ggplot(df, aes(y=reorder(rownames(df), rmse), x=r2)) +
    geom_point() +
    geom_text(label=round(df$r2, 2), vjust=-1, size=3) +
    labs(x='Rsquared', y='') +
    xlim(0.7,1)

plots <- arrangeGrob(rmse.plot, r2.plot, ncol=2)
ggsave(plots, file="R\\plots\\model_comp_rf.png", width=8, height=3)









# ============================================================================
# RMSE and R Squared plot
# ============================================================================

# Predict for testing set
lm.pred <- predict(lmFit, test.scaled)
svm.pred <- predict(svmFit, test.scaled)
knn.pred <- predict(knnFit, test.scaled)
pcr.pred <- predict(pcrFit, test.scaled)
nnet.pred <- predict(nnetFit, test.scaled)
bayesglm.pred <- predict(bayesglmFit, test.scaled)

# # Plot predicted vs. observed (test set)
# plot(lm.pred, test.scaled$yield)
# plot(svm.pred, test.scaled$yield)
# plot(knn.pred, test.scaled$yield)
# plot(pcr.pred, test.scaled$yield)
# plot(nnet.pred, test.scaled$yield)
# plot(bayesglm.pred, test.scaled$yield)

# R^2 values
lm.pred.r2 <- cor(lm.pred, test.scaled$yield)
svm.pred.r2 <- cor(svm.pred, test.scaled$yield)
knn.pred.r2 <- cor(knn.pred, test.scaled$yield)
pcr.pred.r2 <- cor(pcr.pred, test.scaled$yield)
nnet.pred.r2 <- cor(nnet.pred, test.scaled$yield)
bayesglm.pred.r2 <- cor(bayesglm.pred, test.scaled$yield)

# RMSE
lm.pred.rmse <- rmse(lm.pred, test.scaled$yield)
svm.pred.rmse <- rmse(svm.pred, test.scaled$yield)
knn.pred.rmse <- rmse(knn.pred, test.scaled$yield)
pcr.pred.rmse <- rmse(pcr.pred, test.scaled$yield)
nnet.pred.rmse <- rmse(nnet.pred, test.scaled$yield)
bayesglm.pred.rmse <- rmse(bayesglm.pred, test.scaled$yield)

# Pretty plot RMSE and Rsquared
df <- data.frame(rmse = c(lm.pred.rmse, svm.pred.rmse, knn.pred.rmse, nnet.pred.rmse, bayesglm.pred.rmse, rf.pred70.rmse),
                 r2 = c(lm.pred.r2, svm.pred.r2, knn.pred.r2, nnet.pred.r2, bayesglm.pred.r2, rf.pred70.r2))
row.names(df) <- c('Linear Model', 'SVM', 'kNN', 'Neural Network', 'Bayes GLM', 'Random Forest')

rmse.plot <- ggplot(df, aes(y=reorder(rownames(df), rmse), x=rmse)) +
    geom_point() +
    geom_text(label=round(df$rmse, 2), vjust=-1, size=3) +
    labs(x='RMSE', y='') +
    xlim(0,20)

r2.plot <- ggplot(df, aes(y=reorder(rownames(df), rmse), x=r2)) +
    geom_point() +
    geom_text(label=round(df$r2, 2), vjust=-1, size=3) +
    labs(x='Rsquared', y='') +
    xlim(0.7,1)

plots <- arrangeGrob(rmse.plot, r2.plot, ncol=2)
ggsave(plots, file="R\\plots\\model_comp_other.png", width=8, height=3)




# ============================================================================
# Calibration plots
# ============================================================================

# k-nearest neighbors (kNN)
df1 <- data.frame(x = knn.pred, 
                  y = test.scaled$yield,
                  type = as.factor('kNN'))
p1 <- ggplot(df1, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    # scale_x_continuous(breaks = seq(-25,100,25), lim=c(-30, 100)) +
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)

# support vector machine (SVM)
df2 <- data.frame(x = svm.pred, 
                  y = test.scaled$yield,
                  type = as.factor('SVM'))
p2 <- ggplot(df2, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    # scale_x_continuous(breaks = seq(-25,100,25), lim=c(-30, 100)) +
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)

# bayes generalized linear model (GLM)
df3 <- data.frame(x = bayesglm.pred, 
                  y = test.scaled$yield,
                  type = as.factor('Bayes GLM'))
p3 <- ggplot(df3, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(-25,100,25), lim=c(-30, 100)) +
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)

# linear model
df4 <- data.frame(x = lm.pred, 
                 y = test.scaled$yield,
                 type = as.factor('Linear Model'))
p4 <- ggplot(df4, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(-25,100,25), lim=c(-30, 100)) +
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)

# neural network
df5 <- data.frame(x = nnet.pred, 
                  y = test.scaled$yield,
                  type = as.factor('Neural Network'))
p5 <- ggplot(df5, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(-25,100,25), lim=c(-30, 100)) +
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)

# random forest
df6 <- data.frame(x = rf.pred70, 
                  y = test.scaled$yield,
                  type = as.factor('Random Forest'))
p6 <- ggplot(df6, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE) 


facet.df <- do.call(rbind, list(df1, df2, df3, df4, df5, df6)) 
facet.plot <- ggplot(facet.df, aes(x = x, y = y)) +
    geom_point(alpha = 0.3, color="dodgerblue3", size=1) + 
    scale_x_continuous(breaks = seq(-25,100,25), lim=c(-25, 100)) +
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed", size=0.3) +
    geom_smooth(method="loess", se=FALSE, size=0.5, color="black") +
    facet_wrap(~type, ncol=2) +
    labs(x='Predicted Yield', y='Observed Yield')

ggsave(file="R\\plots\\other_plots.png", width=8, height=9)




####################################################
# REGULARIZED MODELS (Lasso, Ridge, Elastic Nets)
####################################################

library(glmnet)

set.seed(1533)
x <- as.matrix(training.scaled[, 1:120])
x.test <- as.matrix(test.scaled[, 1:120])
y <- training.scaled[, 121]
y.test <- test.scaled[, 121]

fit.lasso <- glmnet(x, y, family="gaussian", alpha=1) # lasso
fit.elnet <- glmnet(x, y, family="gaussian", alpha=0.5) # elastic net
fit.ridge <- glmnet(x, y, family="gaussian", alpha=0) # ridge


fit.lasso.cv <- cv.glmnet(x, y, family="gaussian", alpha=1) # lasso
fit.elnet.cv <- cv.glmnet(x, y, family="gaussian", alpha=0.5) # elastic net
fit.ridge.cv <- cv.glmnet(x, y, family="gaussian", alpha=0) # ridge

for (i in 0:10) {
    assign(paste("fit", i, sep=""), cv.glmnet(x, y, type.measure="mse", 
                                              alpha=i/10,family="gaussian"))
}



# Plot solution paths:

png(filename="R\\plots\\regularization.png", width = 800, height = 1000)
par(mfrow=c(3,2))

plot(fit.lasso, xvar="lambda")
plot(fit10, main="LASSO")

plot(fit.ridge, xvar="lambda")
plot(fit0, main="Ridge")

plot(fit.elnet, xvar="lambda")
plot(fit2, main="Elastic Net")

dev.off()


yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=x.test) # ridge
yhat1 <- predict(fit1, s=fit1$lambda.1se, newx=x.test)
yhat2 <- predict(fit2, s=fit2$lambda.1se, newx=x.test)
yhat3 <- predict(fit3, s=fit3$lambda.1se, newx=x.test)
yhat4 <- predict(fit4, s=fit4$lambda.1se, newx=x.test)
yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=x.test)
yhat6 <- predict(fit6, s=fit6$lambda.1se, newx=x.test)
yhat7 <- predict(fit7, s=fit7$lambda.1se, newx=x.test)
yhat8 <- predict(fit8, s=fit8$lambda.1se, newx=x.test)
yhat9 <- predict(fit9, s=fit9$lambda.1se, newx=x.test)
yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=x.test) # lasso

rmse0 <- rmse(yhat0, y.test)
rmse1 <- rmse(yhat1, y.test)
rmse2 <- rmse(yhat2, y.test)
rmse3 <- rmse(yhat3, y.test)
rmse4 <- rmse(yhat4, y.test)
rmse5 <- rmse(yhat5, y.test)
rmse6 <- rmse(yhat6, y.test)
rmse7 <- rmse(yhat7, y.test)
rmse8 <- rmse(yhat8, y.test)
rmse9 <- rmse(yhat9, y.test)
rmse10 <- rmse(yhat10, y.test)

rsquared0 <- cor(yhat0, y.test)
rsquared1 <- cor(yhat1, y.test)
rsquared2 <- cor(yhat2, y.test)
rsquared3 <- cor(yhat3, y.test)
rsquared4 <- cor(yhat4, y.test)
rsquared5 <- cor(yhat5, y.test)
rsquared6 <- cor(yhat6, y.test)
rsquared7 <- cor(yhat7, y.test)
rsquared8 <- cor(yhat8, y.test)
rsquared9 <- cor(yhat9, y.test)
rsquared10 <- cor(yhat10, y.test)


# manual 
mfit1 <- cv.glmnet(x, y, type.measure="mse", alpha=0.01,family="gaussian")
mfit2 <- cv.glmnet(x, y, type.measure="mse", alpha=0.1,family="gaussian")
mfit3 <- cv.glmnet(x, y, type.measure="mse", alpha=0.2,family="gaussian")
mfit4 <- cv.glmnet(x, y, type.measure="mse", alpha=0.5,family="gaussian")
mfit5 <- cv.glmnet(x, y, type.measure="mse", alpha=1,family="gaussian")

mfit1.pred <- predict(mfit1, s=mfit1$lambda.1se, newx=x.test)
mfit2.pred <- predict(mfit2, s=mfit2$lambda.1se, newx=x.test)
mfit3.pred <- predict(mfit3, s=mfit3$lambda.1se, newx=x.test)
mfit4.pred <- predict(mfit4, s=mfit4$lambda.1se, newx=x.test)
mfit5.pred <- predict(mfit5, s=mfit5$lambda.1se, newx=x.test)

mfit1.rmse <- rmse(mfit1.pred, y.test)
mfit2.rmse <- rmse(mfit2.pred, y.test)
mfit3.rmse <- rmse(mfit3.pred, y.test)
mfit4.rmse <- rmse(mfit4.pred, y.test)
mfit5.rmse <- rmse(mfit5.pred, y.test)

mfit1.r2 <- cor(mfit1.pred, y.test)
mfit2.r2 <- cor(mfit2.pred, y.test)
mfit3.r2 <- cor(mfit3.pred, y.test)
mfit4.r2 <- cor(mfit4.pred, y.test)
mfit5.r2 <- cor(mfit5.pred, y.test)

coef(mfit3)

# log scale alpha
df <- data.frame(alpha = c('0 (LASSO)', '0.01', '0.1', '0.2', '0.5', '1 (Ridge)'),
                 rmse = c(rmse0, mfit1.rmse, mfit2.rmse, mfit3.rmse, mfit4.rmse, mfit5.rmse),
                 r2 = c(rsquared0, mfit1.r2, mfit2.r2, mfit3.r2, mfit4.r2, mfit5.r2))

rmse.plot <- ggplot(df, aes(x=rmse, y=alpha)) +
    geom_point() +
    geom_text(label=round(df$rmse, 2), vjust=-1, size=3) +
    labs(x='RMSE', y='alpha') +    
    xlim(15,16.5)

r2.plot <- ggplot(df, aes(x=r2, y=alpha)) +
    geom_point() +
    geom_text(label=round(df$r2, 4), vjust=-1, size=3) +
    labs(x='Rsquared', y='alpha') +    
    xlim(0.8, 0.82)

plots <- arrangeGrob(rmse.plot, r2.plot, ncol=2)
ggsave(plots, file="R\\plots\\model_comp_regularized_3.png", width=8, height=3) # KEEPER




# plot with RMSE and Rsquared in y direction
df <- data.frame(alpha = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                 rmse = c(rmse0, rmse2, rmse4, rmse6, rmse8, rmse10),
                 r2 = c(rsquared0, rsquared2, rsquared4, rsquared6, rsquared8, rsquared10))

rmse.plot <- ggplot(df, aes(x=alpha, y=rmse)) +
    geom_point() +
    scale_x_continuous(breaks = pretty(df$alpha, n = 5)) +
    geom_text(label=round(df$rmse, 2), vjust=-1, size=3) +
    labs(x='alpha', y='RMSE') +    
    ylim(15,16.5)

r2.plot <- ggplot(df, aes(x=alpha, y=r2)) +
    geom_point() +
    scale_x_continuous(breaks = pretty(df$alpha, n = 5)) +
    geom_text(label=round(df$r2, 2), vjust=-1, size=3) +
    labs(x='alpha', y='Rsquared') +    
    ylim(0.8, 0.82)

plots <- arrangeGrob(rmse.plot, r2.plot, nrow=2)
ggsave(plots, file="R\\plots\\model_comp_regularized.png", width=5, height=6)



# Plot with RMSE and Rsquared in x direction
df <- data.frame(alpha = c('alpha = 0', 'alpha = 0.2', 'alpha = 0.4', 'alpha = 0.6', 'alpha = 0.8', 'alpha = 1'),
                 rmse = c(rmse0, rmse2, rmse4, rmse6, rmse8, rmse10),
                 r2 = c(rsquared0, rsquared2, rsquared4, rsquared6, rsquared8, rsquared10))

rmse.plot <- ggplot(df, aes(x=rmse, y=alpha)) +
    geom_point() +
    geom_text(label=round(df$rmse, 2), vjust=-1, size=3) +
    labs(x='RMSE', y='') +    
    xlim(15,16.5)

r2.plot <- ggplot(df, aes(x=r2, y=alpha)) +
    geom_point() +
    geom_text(label=round(df$r2, 2), vjust=-1, size=3) +
    labs(x='Rsquared', y='') +    
    xlim(0.8, 0.82)

plots <- arrangeGrob(rmse.plot, r2.plot, ncol=2)
ggsave(plots, file="R\\plots\\model_comp_regularized_2.png", width=8, height=3)



















#####################################################
# SUBSET SELECTION
#####################################################

additive <- read.csv("R\\additive.csv", header=TRUE)
aryl.halide <- read.csv("R\\aryl_halide.csv", header=TRUE)
base <- read.csv("R\\base.csv", header=TRUE)
ligand <- read.csv("R\\ligand.csv", header=TRUE)

# Correlation plots
ggcorr(additive, label=TRUE, alpha=0)
ggsave(file="R\\corr_matrix\\additive_corrplot.png", width=10, height=10)
ggcorr(aryl.halide, label=TRUE, alpha=0)
ggsave(file="R\\corr_matrix\\aryl_halide_corrplot.png", width=10, height=10)
ggcorr(base, label=TRUE, alpha=0)
ggsave(file="R\\corr_matrix\\base_corrplot.png", width=5, height=5)
ggcorr(ligand, label=FALSE, alpha=0)
ggsave(file="R\\corr_matrix\\ligand_corrplot.png", width=25, height=25)

# Feature selection by removing correlated features

# remove name column (need for cor function)
additive.num <- additive[, -which(names(additive)=="name")]
aryl.halide.num <- aryl.halide[, -which(names(aryl.halide)=="name")]
base.num <- base[, -which(names(base)=="name")]
ligand.num <- ligand[, -which(names(ligand)=="name")]

# additive
additive.bad <- findCorrelation(cor(additive.num), cutoff = 0.50, exact = TRUE)
additive.good <- names(additive.num[, -additive.bad])
additive.reduced <- additive.num[, additive.good]
pca.additive.reduced <- prcomp(additive.reduced)
plot(pca.additive.reduced)
ggcorr(additive.reduced, label=TRUE, alpha=0)
ggsave(file="R\\corr_matrix\\additive_reduced_corrplot.png", width=6, height=6)

# aryl halide
aryl.halide.bad <- findCorrelation(cor(aryl.halide.num), cutoff = 0.50, exact = TRUE)
aryl.halide.good <- names(aryl.halide.num[, -aryl.halide.bad])
aryl.halide.reduced <- aryl.halide.num[, aryl.halide.good]
pca.aryl.halide.reduced <- prcomp(aryl.halide.reduced)
plot(pca.aryl.halide.reduced)
ggcorr(aryl.halide.reduced, label=TRUE, alpha=0)
ggsave(file="R\\corr_matrix\\aryl_halide_reduced_corrplot.png", width=6, height=6)

# base
base.bad <- findCorrelation(cor(base.num), cutoff = 0.50, exact = TRUE)
base.good <- names(base.num[, -base.bad])
base.reduced <- base.num[, base.good]
pca.base.reduced <- prcomp(base.reduced)
plot(pca.base.reduced)
ggcorr(base.reduced, label=TRUE, alpha=0)
ggsave(file="R\\corr_matrix\\base_reduced_corrplot.png", width=3, height=3)

# ligand
ligand.bad <- findCorrelation(cor(ligand.num), cutoff = 0.50, exact = TRUE)
ligand.good <- names(ligand.num[, -ligand.bad])
ligand.reduced <- ligand.num[, ligand.good]
pca.ligand.reduced <- prcomp(ligand.reduced)
plot(pca.ligand.reduced)
ggcorr(ligand.reduced, label=TRUE, alpha=0)
ggsave(file="R\\corr_matrix\\ligand_reduced_corrplot.png", width=4, height=4)

# subset to smaller number of dimensions
reduced.vars = c(names(additive.reduced),
                 names(aryl.halide.reduced),
                 names(base.reduced),
                 names(ligand.reduced),
                 "yield")
training.scaled.reduced <- training.scaled[, reduced.vars]

output.plate1.reduced <- output.plate1[, reduced.vars]
output.plate2.reduced <- output.plate2[, reduced.vars]
output.plate3.reduced <- output.plate3[, reduced.vars]










# ============================================================================
# Out-of-sample Prediction Plots
# ============================================================================

# double check that row numbers are correct by testing number of unique additives
# length(unique(output.scaled$additive_.C3_NMR_shift[1:1075])) == 6 # one more row is 7
# length(unique(output.scaled$additive_.C3_NMR_shift[1:2515])) == 14 # one more row is 15
# length(unique(output.scaled$additive_.C3_NMR_shift)) == 22

# plates 1 and 2 (after NA's removed)
plate12 <- output.scaled[1:2515, ]

# plate 3 (after NA's removed)
plate3 <- output.scaled[2516:3955, ]

by.additive <- split(seq_along(plate12$additive_.C3_NMR_shift), plate12$additive_.C3_NMR_shift)
tc.additive <- trainControl(method="cv", indexOut=by.additive, savePredictions = TRUE)


# # Random forest (70% of data) - LEAVE ONE OUT BY ADDITIVE
# rfFit.LOO <- train(yield ~ ., data=plate12, trControl=tc.additive, method="rf", importance=TRUE)
# png(filename="R\\plots\\rf_LOO.png", width = 1000, height = 600)
# predVals <- extractPrediction(list(rfFit.LOO))
# plotObsVsPred(predVals)
# dev.off()
# saveRDS(rfFit.LOO, "rds\\rfFit_LOO.rds")

rf.predTrain.LOO <- predict(rfFit.LOO, plate12)
rf.predTrain.LOO.rmse <- rmse(rf.predTrain.LOO, plate12$yield)
rf.predTrain.LOO.r2 <- cor(rf.predTrain.LOO, plate12$yield)

rf.pred.LOO <- predict(rfFit.LOO, plate3)
rf.pred.LOO.rmse <- rmse(rf.pred.LOO, plate3$yield)
rf.pred.LOO.r2 <- cor(rf.pred.LOO, plate3$yield)


plate3$additive_id <- as.factor(plate3$additive_.C3_NMR_shift)
levels(plate3$additive_id) <- c('Additive 16', 'Additive 18', 'Additive 20', 'Additive 21', 
                                'Additive 22', 'Additive 17', 'Additive 19', 'Additive 23')
plate3$additive_id <- sortLvls.fnc(plate3$additive_id, c(1, 6, 2, 7, 3, 4, 5, 8))
plate3$additive_id = factor(plate3$additive_id,levels(plate3$additive_id)[c(1, 6, 2, 7, 3, 4, 5, 8)])
df <- cbind(plate3, rf.pred.LOO)
p <- ggplot(df, aes(x=rf.pred.LOO, y=yield)) +
    geom_point(alpha=0.4, aes(col=additive_id), size=1) +
    labs(x='Predicted Yield', y='Observed Yield') +   
    xlim(0, 100) +
    ylim(0, 100) + 
    geom_smooth(method='lm', se=FALSE, color="black", size=0.5) +
    facet_wrap(~additive_id, nrow=2, ncol=4) +
    geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed", size=0.3) +
    theme(legend.position="none")
ggsave(file="R\\plots\\additive_LOO.png", width=8, height=4.5)

by.additive.p3 <- split(seq_along(plate3$additive_.C3_NMR_shift), plate3$additive_.C3_NMR_shift)

















####################################################
# MODEL COMPARISON
####################################################

# Saved models
lmFit <- readRDS("rds\\lmFit.rds")
svmFit <- readRDS("rds\\svmFit.rds")
knnFit <- readRDS("rds\\knnFit.rds")
pcrFit <- readRDS("rds\\pcrFit.rds")
nnetFit <- readRDS("rds\\nnetFit.rds")
bayesglmFit <- readRDS("rds\\bayesglmFit.rds")

# rfFit.LOO <- readRDS("rds\\rfFit_LOO.rds")
# xgbFit.LOO <- readRDS("rds\\xgbFit_LOO.rds")

rfFit10 <- readRDS("rds\\rfFit10.rds")
rfFit20 <- readRDS("rds\\rfFit20.rds")
rfFit30 <- readRDS("rds\\rfFit30.rds")
rfFit50 <- readRDS("rds\\rfFit50.rds")
rfFit70 <- readRDS("rds\\rfFit70.rds")

rfFit.LOO <- readRDS("rds\\rfFit_LOO.rds")

xgbFit10 <- readRDS("rds\\xgbFit10.rds")
xgbFit20 <- readRDS("rds\\xgbFit20.rds")
xgbFit30 <- readRDS("rds\\xgbFit30.rds")
xgbFit50 <- readRDS("rds\\xgbFit50.rds")
xgbFit70 <- readRDS("rds\\xgbFit70.rds")


rfFit.ArCl <- readRDS("rds\\rfFit_ArCl.rds")
rfFit.ArBr <- readRDS("rds\\rfFit_ArBr.rds")
rfFit.ArI <- readRDS("rds\\rfFit_ArI.rds")

rfFit.CF3.OMe.pyr2 <- readRDS("rds\\rfFit_CF3-OMe-pyr2.rds")

# Box-and-whiskers plot (training error)
resamps <- resamples(list(lm = lmFit,
                          svm = svmFit,
                          knn = knnFit,
                          pcr = pcrFit,
                          nnet = nnetFit,
                          bayesglm = bayesglmFit,
                          rf10 = rfFit10,
                          rf20 = rfFit20,
                          rf30 = rfFit30,
                          rf50 = rfFit50,
                          rf70 = rfFit70,
                          xgb10 = xgbFit10,
                          xgb20 = xgbFit20,
                          xgb30 = xgbFit30,
                          xgb50 = xgbFit50,
                          xgb70 = xgbFit70))
summary(resamps)

png(filename="R\\plots\\model_comparison_new.png", width = 1200, height = 1000)
bwplot(resamps,
       layout = c(2, 1),
       scales = list(relation = "free"),
       xlim = list(c(0, 30), c(0, 1)))
dev.off()



xgb_resamps <- resamples(list(xgb10 = xgbFit10,
                              xgb20 = xgbFit20,
                              xgb40 = xgbFit40,
                              xgb60 = xgbFit60,
                              xgb80 = xgbFit80))



png(filename="R\\plots\\xgb_comparison.png", width = 800, height = 500)
bwplot(xgb_resamps,
       layout = c(2, 1),
       scales = list(relation = "free"),
       xlim = list(c(0, 25), c(0, 1)))
dev.off()


loo_resamps <- resamples(list(xgb.LOO = xgbFit.LOO,
                              rf.LOO = rfFit.LOO))

png(filename="R\\plots\\loo_comparison.png", width = 800, height = 300)
bwplot(loo_resamps,
       layout = c(2, 1),
       scales = list(relation = "free"),
       xlim = list(c(0, 10), c(0.6, 1)))
dev.off()























