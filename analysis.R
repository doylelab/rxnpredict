# Install packages (if necessary) and load them
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

# Set the working directory to the location of the rxnpredict folder
setwd("C:\\Users\\Derek\\Dropbox\\rxnpredict")



# ============================================================================
# Helper functions for making heatmaps and histograms
# ============================================================================

# heatmap
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

# histogram
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
# Subset the data to prepare for out-of-sample prediction
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

# nonpyridyl aryl halides
nonpyridyl <- sort(c(CF3.Cl, CF3.Br, CF3.I, OMe.Cl, OMe.Br, OMe.I, Et.Cl, Et.Br, Et.I))
nonpyridyl.scaled <- output.scaled[nonpyridyl, ]
nonpyridyl.scaled <- nonpyridyl.scaled[!(is.na(nonpyridyl.scaled$yield)), ]

# pyridyl aryl halides
pyridyl <- sort(c(pyr2.Cl, pyr2.Br, pyr2.I, pyr3.Cl, pyr3.Br, pyr3.I))
pyridyl.scaled <- output.scaled[pyridyl, ]
pyridyl.scaled <- pyridyl.scaled[!(is.na(pyridyl.scaled$yield)), ]

# yields under 80% and over 80%
under80 <- output.scaled$yield<80
output.under80 <- output.scaled[under80, ]
output.under80 <- output.under80[!(is.na(output.under80$yield)), ]

output.over80 <- output.scaled[!under80, ]
output.over80 <- output.over80[!(is.na(output.over80$yield)), ]


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
# Sampled from within training set to avoid using test set data
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

# 10-fold cross-validation
train_control <- trainControl(method="cv", number=10, savePredictions=TRUE)



# ============================================================================
# Read in previously trained models saved as .rds files
# ============================================================================

# Run to read in previously trained models

lmFit.reduced <- readRDS("rds\\lmFit_reduced.rds")
knnFit <- readRDS("rds\\knnFit.rds")
svmFit <- readRDS("rds\\svmFit.rds")
bayesglmFit <- readRDS("rds\\bayesglmFit.rds")
lmFit <- readRDS("rds\\lmFit.rds")
nnetFit <- readRDS("rds\\nnetFit.rds")
rfFit <- readRDS("rds\\rfFit.rds")
rfFit2.5 <- readRDS("rds\\rfFit2_5.rds")
rfFit5 <- readRDS("rds\\rfFit5.rds")
rfFit10 <- readRDS("rds\\rfFit10.rds")
rfFit20 <- readRDS("rds\\rfFit20.rds")
rfFit30 <- readRDS("rds\\rfFit30.rds")
rfFit40 <- readRDS("rds\\rfFit40.rds")
rfFit50 <- readRDS("rds\\rfFit50.rds")
rfFit60 <- readRDS("rds\\rfFit60.rds")
rfFit70 <- readRDS("rds\\rfFit70.rds")
rfFit.LOO <- readRDS("rds\\rfFit_LOO.rds")
rfFit.ArCl <- readRDS("rds\\rfFit_ArCl.rds")
rfFit.ArBr <- readRDS("rds\\rfFit_ArBr.rds")
rfFit.ArI <- readRDS("rds\\rfFit_ArI.rds")
rfFit.ArBr.all <- readRDS("rds\\rfFit_ArBr_all.rds")
rfFit.nonpyridyl <- readRDS("rds\\rfFit_nonpyridyl.rds")
rfFit.under80 <- readRDS("rds\\rfFit_under80.rds")



# ============================================================================
# Regularized linear regression: Lasso, Ridge, and Elastic Net
# ============================================================================

set.seed(1533)
x <- as.matrix(training.scaled[, 1:120])
x.test <- as.matrix(test.scaled[, 1:120])
y <- training.scaled[, 121]
y.test <- test.scaled[, 121]

fit.lasso <- glmnet(x, y, family="gaussian", alpha=1) # lasso
fit.elnet <- glmnet(x, y, family="gaussian", alpha=0.5) # elastic net
fit.ridge <- glmnet(x, y, family="gaussian", alpha=0) # ridge

for (i in 0:10) {
    assign(paste("fit", i, sep=""), cv.glmnet(x, y, type.measure="mse", 
                                              alpha=i/10,family="gaussian"))
}

# Plot solution paths
png(filename="R\\plots\\regularization_solution_paths.png", width = 800, height = 1000)
par(mfrow=c(3,2))
plot(fit.lasso, xvar="lambda")
plot(fit10, main="LASSO")
plot(fit.ridge, xvar="lambda")
plot(fit0, main="Ridge")
plot(fit.elnet, xvar="lambda")
plot(fit2, main="Elastic Net")
dev.off()

# Predict y values and calculate RMSE and R^2 for range of alpha values
yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=x.test) # ridge
yhat1 <- predict(fit1, s=fit1$lambda.1se, newx=x.test) # elastic net
yhat2 <- predict(fit2, s=fit2$lambda.1se, newx=x.test) # elastic net
yhat3 <- predict(fit3, s=fit3$lambda.1se, newx=x.test) # elastic net
yhat4 <- predict(fit4, s=fit4$lambda.1se, newx=x.test) # elastic net
yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=x.test) # elastic net
yhat6 <- predict(fit6, s=fit6$lambda.1se, newx=x.test) # elastic net
yhat7 <- predict(fit7, s=fit7$lambda.1se, newx=x.test) # elastic net
yhat8 <- predict(fit8, s=fit8$lambda.1se, newx=x.test) # elastic net
yhat9 <- predict(fit9, s=fit9$lambda.1se, newx=x.test) # elastic net
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

mfit1 <- cv.glmnet(x, y, type.measure="mse", alpha=0.01, family="gaussian")
mfit2 <- cv.glmnet(x, y, type.measure="mse", alpha=0.1, family="gaussian")
mfit3 <- cv.glmnet(x, y, type.measure="mse", alpha=0.2, family="gaussian")
mfit4 <- cv.glmnet(x, y, type.measure="mse", alpha=0.5, family="gaussian")
mfit5 <- cv.glmnet(x, y, type.measure="mse", alpha=1, family="gaussian")

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

# Plot RMSE and R^2 for different values of alpha
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
ggsave(plots, file="R\\plots\\regularized_models.png", width=8, height=3)



# ============================================================================
# Dimension reduction by removing correlated descriptors
# ============================================================================

# Read in data for each reaction component
# These tables contain 1 row per different molecule (e.g., 22 rows in additive.csv)
additive <- read.csv("R\\additive.csv", header=TRUE)
aryl.halide <- read.csv("R\\aryl_halide.csv", header=TRUE)
base <- read.csv("R\\base.csv", header=TRUE)
ligand <- read.csv("R\\ligand.csv", header=TRUE)

# Correlation plots
ggcorr(additive, label=TRUE, alpha=0)
ggsave(file="R\\plots\\additive_corrplot.png", width=10, height=10)
ggcorr(aryl.halide, label=TRUE, alpha=0)
ggsave(file="R\\plots\\aryl_halide_corrplot.png", width=10, height=10)
ggcorr(base, label=TRUE, alpha=0)
ggsave(file="R\\plots\\base_corrplot.png", width=5, height=5)
ggcorr(ligand, label=FALSE, alpha=0)
ggsave(file="R\\plots\\ligand_corrplot.png", width=25, height=25)

# Feature selection by removing correlated features

# remove name column (need to do for cor function)
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
ggsave(file="R\\plots\\additive_reduced_corrplot.png", width=6, height=6)

# aryl halide
aryl.halide.bad <- findCorrelation(cor(aryl.halide.num), cutoff = 0.50, exact = TRUE)
aryl.halide.good <- names(aryl.halide.num[, -aryl.halide.bad])
aryl.halide.reduced <- aryl.halide.num[, aryl.halide.good]
pca.aryl.halide.reduced <- prcomp(aryl.halide.reduced)
plot(pca.aryl.halide.reduced)
ggcorr(aryl.halide.reduced, label=TRUE, alpha=0)
ggsave(file="R\\plots\\aryl_halide_reduced_corrplot.png", width=6, height=6)

# base
base.bad <- findCorrelation(cor(base.num), cutoff = 0.50, exact = TRUE)
base.good <- names(base.num[, -base.bad])
base.reduced <- base.num[, base.good]
pca.base.reduced <- prcomp(base.reduced)
plot(pca.base.reduced)
ggcorr(base.reduced, label=TRUE, alpha=0)
ggsave(file="R\\plots\\base_reduced_corrplot.png", width=3, height=3)

# ligand
ligand.bad <- findCorrelation(cor(ligand.num), cutoff = 0.50, exact = TRUE)
ligand.good <- names(ligand.num[, -ligand.bad])
ligand.reduced <- ligand.num[, ligand.good]
pca.ligand.reduced <- prcomp(ligand.reduced)
plot(pca.ligand.reduced)
ggcorr(ligand.reduced, label=TRUE, alpha=0)
ggsave(file="R\\plots\\ligand_reduced_corrplot.png", width=4, height=4)

# Subset to smaller number of dimensions
reduced.vars = c(names(additive.reduced),
                 names(aryl.halide.reduced),
                 names(base.reduced),
                 names(ligand.reduced),
                 "yield")
training.scaled.reduced <- training.scaled[, reduced.vars]
test.scaled.reduced <- test.scaled[, reduced.vars]

# Train linear model using reduced set of descriptors
lmFit.reduced <- train(yield ~ ., data=training.scaled.reduced, trControl=train_control, method="lm")
saveRDS(lmFit.reduced, "rds\\lmFit_reduced.rds")

lmFit.reduced.pred <- predict(lmFit.reduced, test.scaled.reduced)
lmFit.reduced.rmse <- rmse(lmFit.reduced.pred, test.scaled.reduced$yield)
lmFit.reduced.r2 <- cor(lmFit.reduced.pred, test.scaled.reduced$yield)

df <- data.frame(x = lmFit.reduced.pred, 
                 y = test.scaled.reduced$yield,
                 type = as.factor('Linear Model'))
p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(-25,100,25), lim=c(-30, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') + 
    geom_segment(aes(x=0, xend=100, y=0, yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\lm_reduced.png", width=5, height=4)



# ============================================================================
# Training various machine learning models (models are saved as .rds files)
# ============================================================================

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
set.seed(8915)
rfFit <- train(yield ~ ., data=training.scaled, trControl=train_control, method="rf", importance=TRUE)
png(filename="R\\plots\\rf.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit, "rds\\rfFit.rds")



# ============================================================================
# RMSE and R^2 plot for different machine learning models
# ============================================================================

# Predict for testing set
lm.pred <- predict(lmFit, test.scaled)
svm.pred <- predict(svmFit, test.scaled)
knn.pred <- predict(knnFit, test.scaled)
nnet.pred <- predict(nnetFit, test.scaled)
bayesglm.pred <- predict(bayesglmFit, test.scaled)
rf.pred <- predict(rfFit, test.scaled)

# R^2 values
lm.pred.r2 <- cor(lm.pred, test.scaled$yield)
svm.pred.r2 <- cor(svm.pred, test.scaled$yield)
knn.pred.r2 <- cor(knn.pred, test.scaled$yield)
nnet.pred.r2 <- cor(nnet.pred, test.scaled$yield)
bayesglm.pred.r2 <- cor(bayesglm.pred, test.scaled$yield)
rf.pred.r2 <- cor(rf.pred, test.scaled$yield)

# RMSE
lm.pred.rmse <- rmse(lm.pred, test.scaled$yield)
svm.pred.rmse <- rmse(svm.pred, test.scaled$yield)
knn.pred.rmse <- rmse(knn.pred, test.scaled$yield)
nnet.pred.rmse <- rmse(nnet.pred, test.scaled$yield)
bayesglm.pred.rmse <- rmse(bayesglm.pred, test.scaled$yield)
rf.pred.rmse <- rmse(rf.pred, test.scaled$yield)

# Plot RMSE and R^2
df <- data.frame(rmse = c(lm.pred.rmse, svm.pred.rmse, knn.pred.rmse, nnet.pred.rmse, bayesglm.pred.rmse, rf.pred.rmse),
                 r2 = c(lm.pred.r2, svm.pred.r2, knn.pred.r2, nnet.pred.r2, bayesglm.pred.r2, rf.pred.r2))
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
ggsave(plots, file="R\\plots\\ml_models.png", width=8, height=3)



# ============================================================================
# Calibration plots for different machine learning models
# ============================================================================

# Create data frames with predicted and observed values for test set
df1 <- data.frame(x = knn.pred, y = test.scaled$yield, type = as.factor('kNN'))
df2 <- data.frame(x = svm.pred, y = test.scaled$yield, type = as.factor('SVM'))
df3 <- data.frame(x = bayesglm.pred, y = test.scaled$yield, type = as.factor('Bayes GLM'))
df4 <- data.frame(x = lm.pred, y = test.scaled$yield, type = as.factor('Linear Model'))
df5 <- data.frame(x = nnet.pred, y = test.scaled$yield, type = as.factor('Neural Network'))
df6 <- data.frame(x = rf.pred, y = test.scaled$yield, type = as.factor('Random Forest'))

# Make calibration plots
facet.df <- do.call(rbind, list(df1, df2, df3, df4, df5, df6)) 
facet.plot <- ggplot(facet.df, aes(x = x, y = y)) +
    geom_point(alpha = 0.3, color="dodgerblue3", size=1) + 
    scale_x_continuous(breaks = seq(-25,100,25), lim=c(-25, 100)) +
    geom_segment(aes(x=0, xend=100, y=0, yend=100), linetype="dashed", size=0.3) +
    geom_smooth(method="loess", se=FALSE, size=0.5, color="black") +
    facet_wrap(~type, ncol=2) +
    labs(x='Predicted Yield', y='Observed Yield')
ggsave(file="R\\plots\\ml_calibration_plots.png", width=8, height=9)



# ============================================================================
# Random forest: Predicting out-of-sample additives
# ============================================================================

# double check that row numbers are correct by testing number of unique additives
# length(unique(output.scaled$additive_.C3_NMR_shift[1:1075])) == 6 # TRUE (one more row is 7)
# length(unique(output.scaled$additive_.C3_NMR_shift[1:2515])) == 14 # TRUE (one more row is 15)
# length(unique(output.scaled$additive_.C3_NMR_shift)) == 22

# plates 1 and 2 (after NA's removed)
plate12 <- output.scaled[1:2515, ]

# plate 3 (after NA's removed)
plate3 <- output.scaled[2516:3955, ]

# leave-one-out by additive
by.additive <- split(seq_along(plate12$additive_.C3_NMR_shift), plate12$additive_.C3_NMR_shift)
tc.additive <- trainControl(method="cv", indexOut=by.additive, savePredictions = TRUE)

# leave-one-out random forest (70% of data)
set.seed(8915)
rfFit.LOO <- train(yield ~ ., data=plate12, trControl=tc.additive, method="rf", importance=TRUE)
png(filename="R\\plots\\rf_LOO.png", width = 1000, height = 600)
predVals <- extractPrediction(list(rfFit.LOO))
plotObsVsPred(predVals)
dev.off()
saveRDS(rfFit.LOO, "rds\\rfFit_LOO.rds")

rf.predTrain.LOO <- predict(rfFit.LOO, plate12)
rf.predTrain.LOO.rmse <- rmse(rf.predTrain.LOO, plate12$yield)
rf.predTrain.LOO.r2 <- cor(rf.predTrain.LOO, plate12$yield)

rf.pred.LOO <- predict(rfFit.LOO, plate3)
rf.pred.LOO.rmse <- rmse(rf.pred.LOO, plate3$yield)
rf.pred.LOO.r2 <- cor(rf.pred.LOO, plate3$yield)

# Create calibration plots for additives
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
    geom_segment(aes(x=0, xend=100, y=0, yend=100), linetype="dashed", size=0.3) +
    theme(legend.position="none")
ggsave(file="R\\plots\\additive_out_of_sample.png", width=8, height=4.5)

# Calculate RMSE and R^2 for additives 16-23
additive.16 <- plate3[plate3$additive_id=='Additive 16', ]
rf.pred.16 <- predict(rfFit.LOO, additive.16)
rf.pred.16.rmse <- rmse(rf.pred.16, additive.16$yield)
rf.pred.16.r2 <- cor(rf.pred.16, additive.16$yield)

additive.17 <- plate3[plate3$additive_id=='Additive 17', ]
rf.pred.17 <- predict(rfFit.LOO, additive.17)
rf.pred.17.rmse <- rmse(rf.pred.17, additive.17$yield)
rf.pred.17.r2 <- cor(rf.pred.17, additive.17$yield)

additive.18 <- plate3[plate3$additive_id=='Additive 18', ]
rf.pred.18 <- predict(rfFit.LOO, additive.18)
rf.pred.18.rmse <- rmse(rf.pred.18, additive.18$yield)
rf.pred.18.r2 <- cor(rf.pred.18, additive.18$yield)

additive.19 <- plate3[plate3$additive_id=='Additive 19', ]
rf.pred.19 <- predict(rfFit.LOO, additive.19)
rf.pred.19.rmse <- rmse(rf.pred.19, additive.19$yield)
rf.pred.19.r2 <- cor(rf.pred.19, additive.19$yield)

additive.20 <- plate3[plate3$additive_id=='Additive 20', ]
rf.pred.20 <- predict(rfFit.LOO, additive.20)
rf.pred.20.rmse <- rmse(rf.pred.20, additive.20$yield)
rf.pred.20.r2 <- cor(rf.pred.20, additive.20$yield)

additive.21 <- plate3[plate3$additive_id=='Additive 21', ]
rf.pred.21 <- predict(rfFit.LOO, additive.21)
rf.pred.21.rmse <- rmse(rf.pred.21, additive.21$yield)
rf.pred.21.r2 <- cor(rf.pred.21, additive.21$yield)

additive.22 <- plate3[plate3$additive_id=='Additive 22', ]
rf.pred.22 <- predict(rfFit.LOO, additive.22)
rf.pred.22.rmse <- rmse(rf.pred.22, additive.22$yield)
rf.pred.22.r2 <- cor(rf.pred.22, additive.22$yield)

additive.23 <- plate3[plate3$additive_id=='Additive 23', ]
rf.pred.23 <- predict(rfFit.LOO, additive.23)
rf.pred.23.rmse <- rmse(rf.pred.23, additive.23$yield)
rf.pred.23.r2 <- cor(rf.pred.23, additive.23$yield)

# Print RMSE and R^2 to console
paste0("Additive 16: RMSE = ", rf.pred.16.rmse, ", R^2 = ", rf.pred.16.r2)
paste0("Additive 17: RMSE = ", rf.pred.17.rmse, ", R^2 = ", rf.pred.17.r2)
paste0("Additive 18: RMSE = ", rf.pred.18.rmse, ", R^2 = ", rf.pred.18.r2)
paste0("Additive 19: RMSE = ", rf.pred.19.rmse, ", R^2 = ", rf.pred.19.r2)
paste0("Additive 20: RMSE = ", rf.pred.20.rmse, ", R^2 = ", rf.pred.20.r2)
paste0("Additive 21: RMSE = ", rf.pred.21.rmse, ", R^2 = ", rf.pred.21.r2)
paste0("Additive 22: RMSE = ", rf.pred.22.rmse, ", R^2 = ", rf.pred.22.r2)
paste0("Additive 23: RMSE = ", rf.pred.23.rmse, ", R^2 = ", rf.pred.23.r2)



# ============================================================================
# Random forest under sparsity: Training models
# ============================================================================

# Random forest (2.5% of data)
set.seed(8915)
rfFit2.5 <- train(yield ~ ., data=training2.5, trControl=train_control, method="rf")
saveRDS(rfFit2.5, "rds\\rfFit2_5.rds")

# Random forest (5% of data)
set.seed(8915)
rfFit5 <- train(yield ~ ., data=training5, trControl=train_control, method="rf")
saveRDS(rfFit5, "rds\\rfFit5.rds")

# Random forest (10% of data)
set.seed(8915)
rfFit10 <- train(yield ~ ., data=training10, trControl=train_control, method="rf")
saveRDS(rfFit10, "rds\\rfFit10.rds")

# Random forest (20% of data)
set.seed(8915)
rfFit20 <- train(yield ~ ., data=training20, trControl=train_control, method="rf")
saveRDS(rfFit20, "rds\\rfFit20.rds")

# Random forest (30% of data)
set.seed(8915)
rfFit30 <- train(yield ~ ., data=training30, trControl=train_control, method="rf")
saveRDS(rfFit30, "rds\\rfFit30.rds")

# Random forest (40% of data)
set.seed(8915)
rfFit40 <- train(yield ~ ., data=training40, trControl=train_control, method="rf")
saveRDS(rfFit40, "rds\\rfFit40.rds")

# Random forest (50% of data)
set.seed(8915)
rfFit50 <- train(yield ~ ., data=training50, trControl=train_control, method="rf")
saveRDS(rfFit50, "rds\\rfFit50.rds")

# Random forest (60% of data)
set.seed(8915)
rfFit60 <- train(yield ~ ., data=training60, trControl=train_control, method="rf")
saveRDS(rfFit60, "rds\\rfFit60.rds")

# Random forest (70% of data - all training data - rfFit70 is identical to rfFit)
set.seed(8915)
rfFit70 <- train(yield ~ ., data=training.scaled, trControl=train_control, method="rf", importance=TRUE)
saveRDS(rfFit70, "rds\\rfFit70.rds")



# ============================================================================
# Random forest under sparsity: Making calibration plots
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

# Plot expected vs. observed

# Create data frames
df1 <- data.frame(x = rf.pred2.5, y = test.scaled$yield, type = as.factor('2.5%'))
df2 <- data.frame(x = rf.pred5, y = test.scaled$yield, type = as.factor('5%'))
df3 <- data.frame(x = rf.pred10, y = test.scaled$yield, type = as.factor('10%'))
df4 <- data.frame(x = rf.pred20, y = test.scaled$yield, type = as.factor('20%'))
df5 <- data.frame(x = rf.pred30, y = test.scaled$yield, type = as.factor('30%'))
df6 <- data.frame(x = rf.pred40, y = test.scaled$yield, type = as.factor('40%'))
df7 <- data.frame(x = rf.pred50, y = test.scaled$yield, type = as.factor('50%'))
df8 <- data.frame(x = rf.pred60, y = test.scaled$yield, type = as.factor('60%'))
df9 <- data.frame(x = rf.pred70, y = test.scaled$yield, type = as.factor('70%'))

facet.df <- do.call(rbind, list(df1, df2, df3, 
                                df4, df5, df6, 
                                df7, df8, df9)) 
facet.plot <- ggplot(facet.df, aes(x = x, y = y)) +
    geom_point(alpha = 0.3, color="dodgerblue3", size=1) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    geom_segment(aes(x=0, xend=100, y=0, yend=100), linetype="dashed", size=0.3) +
    geom_smooth(method="loess", se=FALSE, size=0.5, color="black") +
    facet_wrap(~type, ncol=3) +
    labs(x='Predicted Yield', y='Observed Yield')

ggsave(file="R\\plots\\rf_sparsity.png", width=8, height=9)



# ============================================================================
# Random forest under sparsity: Plotting R^2 and RMSE
# ============================================================================

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

# create data frame containing RMSE and R^2 data for sparsity models
df <- data.frame(rmse = c(rf.pred2.5.rmse, 
                          rf.pred5.rmse, 
                          rf.pred10.rmse, 
                          rf.pred20.rmse, 
                          rf.pred30.rmse, 
                          rf.pred50.rmse, 
                          rf.pred70.rmse),
                 r2 = c(rf.pred2.5.r2, 
                        rf.pred5.r2, 
                        rf.pred10.r2, 
                        rf.pred20.r2, 
                        rf.pred30.r2, 
                        rf.pred50.r2, 
                        rf.pred70.r2))
row.names(df) <- c('2.5%', '5%', '10%', '20%', '30%', '50%', '70%')

# Plot RMSE and R^2 data
rmse.plot <- ggplot(df, aes(y=reorder(rownames(df), -rmse), x=rmse)) +
    geom_point() +
    geom_text(label=round(df$rmse, 2), vjust=-1, size=2.5) +
    labs(x='RMSE', y='Training Set Data') +
    xlim(0, 20) +
    coord_flip()
r2.plot <- ggplot(df, aes(y=reorder(rownames(df), -rmse), x=r2)) +
    geom_point() +
    geom_text(label=round(df$r2, 2), vjust=-1, size=2.5) +
    labs(x='Rsquared', y='Training Set Data') +
    xlim(0.7, 1) +
    coord_flip()
plots <- arrangeGrob(r2.plot, rmse.plot, ncol=2)
ggsave(plots, file="R\\plots\\rf_sparsity_r2_rmse.png", width=6, height=3)



# ============================================================================
# Random forest: Plotting variable importance
# ============================================================================

# read in variable importance from trained random forest model
rf_imp <- importance(rfFit70$finalModel)
rf.imp.df <- cbind(as.data.frame(rf_imp), names(rf_imp[, 1]))
colnames(rf.imp.df)[1] <- "IncMSE"
colnames(rf.imp.df)[3] <- "descriptor"

# for descriptor names, replace "_" with " " and "." with "*"
rf.imp.df$descriptor <- gsub("_", " ", rf.imp.df$descriptor)
rf.imp.df$descriptor <- gsub("[.]", "*", rf.imp.df$descriptor)

# capitalize descriptor names
simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
          sep="", collapse=" ")
}
rf.imp.df$descriptor <- sapply(rf.imp.df$descriptor, simpleCap)

# plot variable importance (saves to R\plots\rf_variable_importance.png)
p1 <- ggplot(rf.imp.df[rf.imp.df$IncMSE>15, ], aes(x=reorder(descriptor, IncMSE), y=IncMSE)) +
    geom_bar(stat="identity") +
    scale_y_continuous(labels = comma) +
    labs(x="", y="Increase in Mean Squared Error (%)") + 
    coord_flip()
ggsave(p1, file="R\\plots\\rf_importance.png", width=6, height=4)



# ============================================================================
# Plotting *C3 NMR Shift vs Yield
# ============================================================================

# loess best fit curve
p1 <- ggplot(output.scaled, aes(x = additive_.C3_NMR_shift, y = yield)) +
    geom_point(alpha = 0.2, size=1) + 
    labs(x='Additive *C3 NMR Shift', y='Yield') +  
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\C3_vs_yield_loess.png", width=5, height=4)

# linear best fit line
p2 <- ggplot(output.scaled, aes(x = additive_.C3_NMR_shift, y = yield)) +
    geom_point(alpha = 0.2, size=1) + 
    labs(x='Additive *C3 NMR Shift', y='Yield') + 
    geom_smooth(method="lm", se=FALSE)
ggsave(file="R\\plots\\C3_vs_yield_lm.png", width=5, height=4)

# calculate R^2 value between *C3 NMR shift and yield; display in console
C3.yield.cor <- cor(output.scaled$additive_.C3_NMR_shift, output.scaled$yield)
paste0("R^2 value between *C3 NMR Shift and Yield = ", C3.yield.cor)



# ============================================================================
# Random forest: Train and test each aryl halide individually
# ============================================================================

# Train random forest model (ArCl)
set.seed(8915)
rfFit.ArCl <- train(yield ~ ., data=ArCl.training, trControl=train_control, method="rf", importance=TRUE)
saveRDS(rfFit.ArCl, "rds\\rfFit_ArCl.rds")

# Test on ArCl
ClfromCl <- predict(rfFit.ArCl, ArCl.test)
ClfromCl.r2 <- cor(ClfromCl, ArCl.test$yield)
ClfromCl.rmse <- rmse(ClfromCl, ArCl.test$yield)
df1 <- data.frame(x = ClfromCl, 
                  y = ArCl.test$yield)

# Create calibration plot (predict Cl from Cl using random forest model)
p1 <- ggplot(df1, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0, 100, 25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0, xend=100, y=0, yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\ClfromCl.png", width=5, height=4)

# Train random forest model (ArBr)
set.seed(8915)
rfFit.ArBr <- train(yield ~ ., data=ArBr.training, trControl=train_control, method="rf", importance=TRUE)
saveRDS(rfFit.ArBr, "rds\\rfFit_ArBr.rds")

# Test on ArBr
BrfromBr <- predict(rfFit.ArBr, ArBr.test)
BrfromBr.r2 <- cor(BrfromBr, ArBr.test$yield)
BrfromBr.rmse <- rmse(BrfromBr, ArBr.test$yield)
df2 <- data.frame(x = BrfromBr, 
                  y = ArBr.test$yield)

# Create calibration plot (predict Br from Br using random forest model)
p2 <- ggplot(df2, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0, xend=100, y=0, yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\BrfromBr.png", width=5, height=4)

# Train random forest model (ArI)
set.seed(8915)
rfFit.ArI <- train(yield ~ ., data=ArI.training, trControl=train_control, method="rf", importance=TRUE)
saveRDS(rfFit.ArI, "rds\\rfFit_ArI.rds")

# Test on ArI
IfromI <- predict(rfFit.ArI, ArI.test)
IfromI.r2 <- cor(IfromI, ArI.test$yield)
IfromI.rmse <- rmse(IfromI, ArI.test$yield)
df3 <- data.frame(x = IfromI, 
                  y = ArI.test$yield)

# Create calibration plot (predict Br from Br using random forest model)
p3 <- ggplot(df3, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0, xend=100, y=0, yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\IfromI.png", width=5, height=4)

# Calculate R^2 and RMSE
rf.pred.ArCl <- predict(rfFit.ArCl, ArCl.test)
rf.pred.ArCl.r2 <- cor(rf.pred.ArCl, ArCl.test$yield)
rf.pred.ArCl.rmse <- rmse(rf.pred.ArCl, ArCl.test$yield)

rf.pred.ArBr <- predict(rfFit.ArBr, ArBr.test)
rf.pred.ArBr.r2 <- cor(rf.pred.ArBr, ArBr.test$yield)
rf.pred.ArBr.rmse <- rmse(rf.pred.ArBr, ArBr.test$yield)

rf.pred.ArI <- predict(rfFit.ArI, ArI.test)
rf.pred.ArI.r2 <- cor(rf.pred.ArI, ArI.test$yield)
rf.pred.ArI.rmse <- rmse(rf.pred.ArI, ArI.test$yield)

# Pretty plot RMSE and Rsquared
df <- data.frame(rmse = c(rf.pred.ArCl.rmse, rf.pred.ArBr.rmse, rf.pred.ArI.rmse),
                 r2 = c(rf.pred.ArCl.r2, rf.pred.ArBr.r2, rf.pred.ArCl.r2))
row.names(df) <- c('ArCl', 'ArBr', 'ArI')

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
ggsave(plots, file="R\\plots\\rf_aryl_halides.png", width=6, height=2)



# ============================================================================
# Random forest: Predict ArCl and ArI from ArBr
# ============================================================================

# Random forest (ArBr)
set.seed(8915)
rfFit.ArBr.all <- train(yield ~ ., data=ArBr.scaled, trControl=train_control, method="rf", importance=TRUE)
saveRDS(rfFit.ArBr.all, "rds\\rfFit_ArBr_all.rds")

# Test on ArCl
ClfromBr <- predict(rfFit.ArBr.all, ArCl.scaled)
ClfromBr.r2 <- cor(ClfromBr, ArCl.scaled$yield)
ClfromBr.rmse <- rmse(ClfromBr, ArCl.scaled$yield)
df1 <- data.frame(x = ClfromBr, 
                  y = ArCl.scaled$yield)

# Create calibration plot (predict Cl from Br random forest model)
p1 <- ggplot(df1, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0, xend=100, y=0, yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\ClfromBr.png", width=5, height=4)

# Test on ArI
IfromBr <- predict(rfFit.ArBr.all, ArI.scaled)
IfromBr.r2 <- cor(IfromBr, ArI.scaled$yield)
IfromBr.rmse <- rmse(IfromBr, ArI.scaled$yield)
df2 <- data.frame(x = IfromBr, 
                  y = ArI.scaled$yield)

# Create calibration plot (predict I from Br random forest model)
p2 <- ggplot(df2, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0, xend=100, y=0, yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\IfromBr.png", width=5, height=4)



# ============================================================================
# Random forest: Predicting pyridyl substrates from nonpyridyl ones
# ============================================================================

# Random forest (train on nonpyridyl)
set.seed(8915)
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
    geom_segment(aes(x=0, xend=100, y=0, yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\pyr_from_nonpyr.png", width=5, height=4)



# ============================================================================
# Random forest: Predicting yields > 80% from yields < 80%
# ============================================================================

# Train random forest model using yields < 80%
set.seed(8915)
rfFit.under80 <- train(yield ~ ., data=output.under80, trControl=train_control, method="rf", importance=TRUE)
saveRDS(rfFit.under80, "rds\\rfFit_under80.rds")

# Predict yields for reactions > 80%
over80pred <- predict(rfFit.under80, output.over80)
over80pred.r2 <- cor(over80pred, output.over80$yield)
over80pred.rmse <- rmse(over80pred, output.over80$yield)

# Generate calibration plot
df <- data.frame(x = over80pred, 
                 y = output.over80$yield)
p1 <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.4) + 
    scale_x_continuous(breaks = seq(0, 100, 25), lim=c(0, 100)) +
    labs(x='Predicted Yield', y='Observed Yield') +  
    geom_segment(aes(x=0, xend=100, y=0, yend=100), linetype="dashed") +
    geom_smooth(method="loess", se=FALSE)
ggsave(file="R\\plots\\over80pred.png", width=5, height=4)