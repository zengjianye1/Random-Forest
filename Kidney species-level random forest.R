###Data Preprocessing
gphylum <- read.csv('5species.csv',header=T,row.names = 1)
##The number 5 in 5species represents the spleen.
gphylum <- gphylum[which(rowSums(gphylum) >= 0.0005), ]#Filter low abundance microorganisms
gsj <- read.csv('5ksj.csv',header=T,row.names = 1)#For 5ksj, '5' represents the lung sample, 'k' is the first letter of kidney, and 'sj' represents PMI.
hist(gsj$PMI,breaks = 50)
gphylum<- data.frame(t(gphylum))
gphylum <- gphylum[rownames(gsj), ]
gphylum <- cbind(gphylum, gsj)
otu <- data.frame(gphylum)
##Loading Packages
library(skimr )
library(DataExplorer)
library(caret)
##The total dataset is divided into a training set (75%) and a test set (25%)
set.seed(123)
trains <- createDataPartition(y=otu$PMI,p=0.75,list=F)
traindata <- data.frame(otu[trains,])
testdata <- data.frame(otu[-trains,])
hist(traindata$PMI,breaks = 50)
hist(testdata$PMI,breaks = 50)
colnames(otu)
##Building a Random Forest Model
form_reg <- as.formula(paste0("PMI~",paste(colnames(traindata)[1:154],collapse="+")))
form_reg
set.seed(123)
library(randomForest)
fit_rf_reg <- randomForest(form_reg,data=traindata,ntree=500,mtry=6,importance=T)
fit_rf_reg
plot(fit_rf_reg,main="ERROR&TREES")
##Cross-validation
set.seed(123)
otu_train.cv <- replicate(5, rfcv(traindata[-ncol(traindata)], traindata$PMI, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 20)
##Plotting the cross validation curve
ggplot(otu_train.cv.mean, aes(x = Group.1, y = x)) +
  geom_line() +
  theme_minimal() +
  theme(panel.background = element_rect(fill = 'transparent', color = NA), 
        panel.border = element_rect(color = 'black', fill = NA),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90")) +
  labs(title = 'Cross-Validation Error vs Number of Features', 
       x = 'Number of Features (Bacterial Species)', 
       y = 'Cross-Validation Error')
##Re-extract important OTUs and their corresponding data
importance_otu <- importance(fit_rf_reg)
importance_otu.select <- importance_otu[1:6, ]#Determined based on the cross validation curve and the values of Rsquared and MAE
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu[, c(otu_id.select, 'PMI')]
##Split the dataset into training and testing sets
set.seed(123) 
train_indices <- createDataPartition(otu.select$PMI, p = 0.75, list = FALSE)
traindata_select <- otu.select[train_indices, ]
testdata_select <- otu.select[-train_indices, ]
##Train the model with a new dataset
form_reg_select <- as.formula(paste0("PMI~", paste(otu_id.select, collapse = "+")))
fit_rf_reg_select <- randomForest(form_reg_select, data = traindata_select, ntree = 500, mtry = 6, importance = TRUE)
##prediction
trainpred_select <- predict(fit_rf_reg_select, newdata = traindata_select)
testpred_select <- predict(fit_rf_reg_select, newdata = testdata_select)
##Results Evaluation
defaultSummary(data.frame(obs = traindata_select$PMI, pred = trainpred_select))
defaultSummary(data.frame(obs = testdata_select$PMI, pred = testpred_select))
##Model drawing
predresult_select <- data.frame(obs = c(traindata_select$PMI, testdata_select$PMI),
                                pred = c(trainpred_select, testpred_select),
                                group = c(rep("Train", length(trainpred_select)), rep("Test", length(testpred_select))))
ggplot(predresult_select, aes(x = obs, y = pred, fill = group, colour = group)) +
  geom_point(shape = 21, size = 3) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, size = 1.2) +
  labs(fill = NULL, colour = NULL, x = "Observed Days", y = "Prediction") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.border = element_rect(color = "grey", fill = NA, size = 0.25))
##Extract and rank variable importance
importance_otu <- importance(fit_rf_reg)
importance_otu.select <- data.frame(importance_otu[1:10, ])
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu[, c(otu_id.select, 'PMI')]
##Output variable importance score table
library(knitr)
kable(importance_otu.select)
library(openxlsx)
write.xlsx(importance_otu.select, file = "Kidney species importance_table.xlsx", rowNames = FALSE)
#Draw a variable importance histogram
library(ggthemes)
library(stringr)
ggplot(importance_otu.select, aes(x = reorder(rownames(importance_otu.select), IncNodePurity), y = IncNodePurity, fill = IncNodePurity)) +
  geom_bar(stat = "identity", fill = "#FFD700") +
  geom_text(aes(y = IncNodePurity + 42, label = rownames(importance_otu.select)), hjust = 0.9) +
  theme_few() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16), 
    plot.caption = element_text(size = 12), 
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 15),
    axis.text.y = element_blank(), # Remove the y-axis text
    axis.title.y = element_blank() # Remove the y-axis title
  ) +
  labs(title = "", x = "", y = "IncNodePurity") +
  coord_flip() # Horizontal bar charts are easier to read

# Extracting and ranking variable importance
importance_otu <- importance(fit_rf_reg)
importance_otu.select <- as.data.frame(importance_otu[1:10, ])
importance_otu.select$Variable <- rownames(importance_otu.select)
# Sort by IncNodePurity and select the first 6 variables
importance_otu.select <- importance_otu.select[order(-importance_otu.select$IncNodePurity), ][1:6, ]
# Get the first 6 variable names after sorting
otu_id.select <- importance_otu.select$Variable
otu.select <- otu[, c(otu_id.select, 'PMI')]
# Reshape the data
otu.selects <- reshape2::melt(otu.select, id = 'PMI')
# Convert the variable column to an ordered factor
otu.selects$variable <- factor(otu.selects$variable, levels = otu_id.select)
# Drawing
ggplot(otu.selects, aes(x = PMI, y = value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~variable, ncol = 3, scale = 'free_y') +
  labs(title = '', x = 'PMI (days)', y = 'Relative abundance')
###改
# Extracting and ranking variable importance
importance_otu <- importance(fit_rf_reg)
importance_otu.select <- as.data.frame(importance_otu[1:10, ])
importance_otu.select$Variable <- rownames(importance_otu.select)
# Sort by IncNodePurity and select the first 5 variables
importance_otu.select <- importance_otu.select[order(-importance_otu.select$IncNodePurity), ][1:3, ]
# Get the first 6 variable names after sorting
otu_id.select <- importance_otu.select$Variable
otu.select <- otu[, c(otu_id.select, 'PMI')]
# Reshape the data
otu.selects <- reshape2::melt(otu.select, id = 'PMI')
# Convert the variable column to an ordered factor
otu.selects$variable <- factor(otu.selects$variable, levels = otu_id.select)
# Drawing
ggplot(otu.selects, aes(x = PMI, y = value, color = variable)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("green", "blue", "red")) +  # 自定义颜色，这里根据需要设置颜色
  theme_minimal() +
  theme(
    axis.line = element_line(size = 0.5),
    axis.line.x.bottom = element_line(arrow = arrow(type = "open", length = unit(0.1, "inches"), ends = "last")),
    axis.line.y.left = element_line(arrow = arrow(type = "open", length = unit(0.1, "inches"), ends = "last")),
    axis.text = element_text(size = 12), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14), 
    legend.text = element_text(size = 12)
  ) +
  labs(x = "PMI", y = "Relative Abundance", color = NULL) + 
  guides(color = guide_legend(title = NULL))