# install packages (if necessary) and load them
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, 
               caret,
               ModelMetrics,
               scales)

install.packages("multcomp")
library(multcomp)
install.packages("party")
library(party)
install.packages("rpart")
library(rpart)

# set the working directory to the location of the rxnpredict folder
setwd("/Users/Desktop/")

# ============================================================================
# Load descriptor and yield data and prepare data for modeling
# ============================================================================

# load training set data or entire dataset
output.scaled <- read.csv("CD_trainingset1.csv", header=TRUE)

# load test set data
test.scaled <- read.csv("CD_testset1.csv", header=TRUE)

# scale the descriptor data (Note: this scales entire dataset)
output.scaled <- as.data.frame(scale(output.table))

# load user-created yield data (label reactions without yield data as NA)
yield.data <- as.numeric(unlist(read.csv("yields.csv", header=TRUE, stringsAsFactors=FALSE)))

# append the yield data to the output table
output.scaled$yield <- yield.data

# remove rows where yield=NA and columns containing at least one NA
output.scaled <- output.scaled[!(is.na(output.scaled$yield)), ]
output.scaled <- output.scaled[ , colSums(is.na(output.scaled)) == 0]

# if dataset is not being split into training/test sets, rename output array as training array
training.scaled <- output.scaled
  
# ============================================================================
# Split data and train random forest model
# ============================================================================

# Split into training and test set (70/30)
set.seed(1084)
size <- round(0.70*nrow(output.scaled))     # default is 0.70
training <- sample(nrow(output.scaled), size=size, replace=FALSE)
training.scaled <- output.scaled[training,]
test.scaled <- output.scaled[-training,]

# 10-fold cross-validation
train_control <- trainControl(method="cv", number=10, savePredictions=TRUE)

# train random forest model with cross-validation
rfFit <- train(yield ~ ., data=training.scaled, trControl=train_control, method="rf", importance=TRUE)
saveRDS(rfFit, "CD_trainingset1_rfFit.rds")

# train random forest model without cross-validation
rfFit <- train(yield ~ ., data=training.scaled, method="rf", importance=TRUE)

# train linear model without cross-validation
lmFit <- train(yield ~ ., data=training.scaled, method="lm")

# read in previously trained model
rfFit <- readRDS("CD_trainingset1_rfFit.rds")

# Train RF model with cforest(); p = # of descriptors
cfFit <- cforest(yield ~ ., data=training.scaled, controls=cforest_unbiased(mtry = 40, ntree = 500))

# Train RF model with cforest(); p = # of descriptors
cfFit <- cforest(yield ~ ., data=training.scaled, controls=cforest_unbiased(mtry = 40, ntree = 500))

my_cforest_control <- cforest_control(teststat = "quad",
                                      testtype = "Univ", mincriterion = 0, ntree = 500, mtry = 40,
                                      replace = FALSE)

my_cforest <- cforest(yield ~ ., data = training.scaled, controls = my_cforest_control)

# ============================================================================
# Calculate R^2 and RMSE using test set and generate calibration plot
# ============================================================================

# predict yields for test set
rf.pred <- predict(rfFit, test.scaled)

# predict yields for test set
crf.pred <- predict(cfFit, newdata = test.scaled)

# predict yields for test set
lm.pred <- predict(lmFit, test.scaled)

# create cvs file of predicted yields
write.csv(rf.pred, file = "CD_testset1_predictions.csv")

# calculate R^2 and RMSE for randomforest()
rf.R2 <- 1 - (sum((test.scaled$yield-rf.pred)^2)/sum((test.scaled$yield-mean(test.scaled$yield))^2))
rf.rmse <- rmse(rf.pred, test.scaled$yield)
paste0("R_squared = ", rf.R2, "RMSE = ", rf.rmse)

# calculate R^2 for linear model
lm.R2 <- 1 - (sum((test.scaled$yield-lm.pred)^2)/sum((test.scaled$yield-mean(test.scaled$yield))^2))
paste0("R squared = ", lm.R2 )

# calculate R^2 for cforest() model
crf.R2 <- 1 - (sum((test.scaled$yield-crf.pred)^2)/sum((test.scaled$yield-mean(test.scaled$yield))^2))
paste0("R_squared = ", crf.R2)

# plot calibration plot for randomforest()
df <- data.frame(x = rf.pred,
                 y = test.scaled$yield)
p1 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.4, color="dodgerblue3", size=1) + 
  scale_x_continuous(breaks = seq(0,100,25), lim=c(0, 100)) +
  labs(x='Predicted Yield', y='Observed Yield') +  
  geom_segment(aes(x=0,xend=100,y=0,yend=100), linetype="dashed", size=0.3) +
  theme(axis.text=element_text(size=12), axis.title.x = element_text(face="bold", size=12)) +
  theme(axis.title.y = element_text(face="bold",size=12))
ggsave(file="WXdata2.png", width=4.5, height=4)

# ============================================================================
# Create Variable importance plot for randomforest()
# ============================================================================

# read in variable importance from trained rf model
rf_imp <- importance(rfFit$finalModel)
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

# plot variable importance
# USER: change '10' on next line to modify minimum cutoff for IncMSE
p2 <- ggplot(rf.imp.df[rf.imp.df$IncMSE>1, ], aes(x=reorder(descriptor, IncMSE), y=IncMSE)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels = comma) +
  labs(x="", y="Increase in Mean Squared Error (%)") + 
  coord_flip()
ggsave(file="RF_importances.png", width=7, height=5)

# ============================================================================
# Create Variable importance plot for cforest()
# ============================================================================

# read in variable importance from trained rf model
cf_imp = varimp(cfFit) # insert conditional = TRUE for correlated features
cf.imp.df <- cbind(as.data.frame(cf_imp), names(cf_imp))
colnames(cf.imp.df)[1] <- "IncMSE"
colnames(cf.imp.df)[2] <- "descriptor"

# for descriptor names, replace "_" with " " and "." with "*"
cf.imp.df$descriptor <- gsub("_", " ", cf.imp.df$descriptor)
cf.imp.df$descriptor <- gsub("[.]", "*", cf.imp.df$descriptor)

# capitalize descriptor names
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep="", collapse=" ")
}
cf.imp.df$descriptor <- sapply(cf.imp.df$descriptor, simpleCap)

# plot variable importance
# USER: change '0' on next line to modify minimum cutoff for IncMSE
p2 <- ggplot(cf.imp.df[cf.imp.df$IncMSE>10, ], aes(x=reorder(descriptor, IncMSE), y=IncMSE)) +
  geom_bar(stat="identity") +
  scale_y_continuous(labels = comma) +
  labs(x="", y="Increase in Mean Squared Error (%)") + 
  coord_flip()
ggsave(file="cF_importances.png", width=7, height=5)

# ============================================================================
# Randomizes Yields for Y-randomization test
# ============================================================================

# the following line randomizes yields to perform Y-randomization test
output.scaled$yield <- sample(output.scaled$yield)

# ============================================================================
# Decision Tree (DT) Model
# ============================================================================

# creates decision tree for the entire dataset
fit <- rpart(yield ~ ., method="anova", data=output.scaled)
plot(fit, unifor=TRUE, 
     main="Regression Tree")
     text(fit, use.n=TRUE, cex=0.8)
ggsave(file="UnscaledDTmodel.png", width=4.5, height=4)