getwd()
targets <- readTargets("target.txt") 
targets

#Read CEL Files
data <- ReadAffy(filenames = targets$FileName)
data

#=====================#
# RMA Normalization
#=====================#
eset <- rma(data)
normset <- exprs(eset)

write.csv(normset, "ExpSet_PostNorm.csv", quote = F)


#=====================#
#   Box Plot
#=====================#
par(mfrow=c(1,2))

#Boxplot Before Normalization
tiff(file="Boxplot_Pre-Normalization.tiff", bg="transparent", width=400, height=500)
par(mar = c(12, 4, 6, 2) + 0.1); # This sets the plot margins
boxplot(data,col="red", main="Boxplot Pre-Normalization", las=2, cex.axis=0.74, ylab="Intensities" )
title(xlab = "Sample Array", line = 8); # Add x axis title
dev.off()

#Boxplot After Normalization
tiff(file="Boxplot_Post-Normalization.tiff", bg="transparent", width=400, height=500)
par(mar = c(12, 4, 6, 2) + 0.1); # This sets the plot margins
boxplot(normset,col="blue", main="Boxplot Post-Normalization", las=2, cex.axis=0.74, ylab="Intensities") #, col=colors 
title(xlab = "Sample Array", line = 8); # Add x axis title
dev.off()


#-----------------------#
# DEG Analysis - Part 2
#-----------------------#
# 1. Data Cleaning - Dimensionality Reduction (PCA)
# 2. Differentially Expressed Genes (DEGs) Identification 
# 3. DEGs Visualization (Volcano plot & Heatmap)


#=====================#
#   PCA Plot
#=====================#
install.packages("factoextra")
install.packages("tidyverse")
library(tidyverse)
library(factoextra)
library(limma)


data <- read.csv("ExpSet_PostNorm.csv")
nrow(data)  # Check the number of rows in your data
length(c(rep("Normal", 10), rep("Cervical_cancer", 10)))  # Check the length of your group labels


# Adjust Group Labels
data_t <- t(data[,-1 ])  # Exclude the first col (gene names) during transpose
data_t <- as.data.frame(data_t)
data_t$Group <- c(rep("Normal", 10), rep("Cervical_cancer", 10))
data_t$Group

# Perform PCA
pca_res <- prcomp(data_t[, -ncol(data_t)])

# view all PC scores 
head(pca_res$x)
library(factoextra)


# Save PCA plot as high-resolution TIFF
tiff("pca.tiff", width = 2000, height = 2000, res = 300)

# Generate PCA plot
fviz_pca_ind(pca_res,
             geom.ind = c("point", "text"),
             col.ind = data_t$Group,
             palette = c("red", "blue"),
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Group",
             labelsize = 2  # increased label size for better clarity
)
dev.off()

#=======================#
# DEGs Identification
#=======================#
# Model Matrix Design
# Let's create a model matrix using the factor() function to represent the condition labels ("Normal" and "Cervical_cancer")
design <- model.matrix(~factor(c("Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal","Normal" , "Cervical_cancer","Cervical_cancer","Cervical_cancer","Cervical_cancer","Cervical_cancer","Cervical_cancer","Cervical_cancer","Cervical_cancer","Cervical_cancer","Cervical_cancer")))
# Now assign names ("Normal" and "Cervical_cancer") to the columns of the model matrix
colnames(design) <- c("Normal", "Cervical_cancer")

# Fits a linear model for each gene based on the given series of arrays
# It estimates the relationship between gene expression and conditions

fit <- lmFit(normset[,1:ncol(normset)], design) 
fit

# Contrast Matrix Design
# Define the specific comparison between conditions you want to analyze.
cont.matrix = makeContrasts(Cervical_cancer - Normal, levels=design)
cont.matrix

# Fitting model with Contrasts(2 Groups), so apply the defined contrast to the previously fitted model (fit)
fit2 <- contrasts.fit(fit, cont.matrix)

# Model optimization / Empirical Bayes Moderation
# Improves the estimation of variances for genes with low expression
# Computes moderated t-statistics and log-odds (B-stats) of differential expression by empirical Bayes shrinkage of the standard errors towards a common value
# contrast-specific information from fit2 is incorporated into the fit object, which is then passed to eBayes()
fit2 <- eBayes(fit)  
fit2

# Result Top Table
topTable(fit2, coef = 2, adjust.method = "BH") 

DEGs <- topTable(fit2, coef=2, adjust="BH", sort.by="logFC", number=100000); #inf
DEGs
write.csv(DEGs, "Result_Table_logFCsorted.csv", quote = F, row.names = TRUE)

#==================================================#
# Filter & Save final DEGs based on Pvalue & logFC
#==================================================#
# Read data of topTable
DEGs <- read.csv("Result_Table_logFCsorted.csv", header = TRUE)

# Filter & Save DEGs
final_DEGs <- DEGs[DEGs$P.Value < 0.05 & (DEGs$logFC > 2 | DEGs$logFC < -2), ]
write.csv(final_DEGs,"finalDEGs.csv", quote = F, row.names = F)



#=============================================#
# Annotate(getting Gene Symbols) filtered DEGs
#=============================================#
#BiocManager::install("hgu133plus2.db")
library("hgu133plus2.db")

DEGs <- read.csv("finalDEGs.csv", header = TRUE)
head(DEGs)


#open the the finalDEGs file manually and write the heading of affimetric IDs as Probe_IDs
probes=DEGs$Probe_ID
head(probes)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
head(Symbols)

# Combine gene annotations with raw data
deg_anno = cbind(probes,Symbols, DEGs)
write.csv(deg_anno, "DEGs_Annotated.csv", quote = F, row.names = F)



#================================================#
# DEG Viz (Volcano plot)
#================================================#
install.packages("gdata")
install.packages("gplots")
library(gdata)
library(gplots)

DEGs <- read.csv("Result_Table_logFCsorted.csv", h=T)
head(DEGs)

png(filename = "VolcanoPlot_FC_2.png")
with(DEGs, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot"))
with(subset(DEGs, P.Value < 0.05 & logFC > 2 ), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(DEGs, P.Value < 0.05 & logFC < -2), points(logFC, -log10(P.Value), pch=20, col="green"))
dev.off()



#========================#
# Probe/Gene Annotation
#========================#
BiocManager::install("hgu133plus2.db")
library("hgu133plus2.db")

normset <- read.csv("ExpSet_PostNorm.csv", h=TRUE)
head(normset)

# Match probe IDs and retrieve Gene SYMBOLS
probes=normset$X
head(probes)
Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))

# Combine gene annotations with raw data
normset_anno = cbind(probes,Symbols,normset)
write.csv(normset_anno, "ExpSet_PostNorm_Annotated.csv", quote = F, row.names = F)

setwd("~/Documents/R_projects/GSE63514_RAW/Cervical_Cancer_data")
getwd()
#open DEGs_annotated file and ExpSet_PostNorm_annotated file.
#In exp_postnorm annotated file, probe id will be repeated, son delete any one probe id and rename empty column as Genes
#A heatmap can effectively show upto 50- 70 genes, hence we need to consider only highly upregulated or highly downregulated genes.
#Go to DEGs_annotated file and filter logfc value greater than equal to 2.5 (which was initially 2) 
#and less than or equal to -2.5(which was initially -2)
#select the gene symbols and paste it in new sheet, and remove duplicates. select column- then data- then select remove duplicates
#Go to Exp_postNorm_annotated file and selcet the second row of Genes coluumn
#type =VLOOKUP(slect the fist gene symbol on adjuscent row comma go to the sheet where we have pasted filtered genes and 
#select column comma 1 (because it is 1st column) comma 0 (for exact matching))
#when press enter it will show #NA on first row, doble clich on box, it will autofill column
#copy this entire column and paste it in the 1st row which is heading Genes (while pating select paste as values i.e 123)
#filter the column and remove #NA and NA values
#now we have got the result, just copy pate the genes and expression values in another sheet and name as heatmap_data


#====================================#
# DEG Viz (Heatmap DEG Expression)
#====================================#
data <- read.csv(file = "Heatmap_data.csv", h=T)
head(data)


#to know the smallest values of gene expression , write =min(select the all 10 columns with expression data)
#for highest =max(select columns ,close bracket, enter)
#based on which give values in the col_breaks
#low exp red, high blue

rnames <- data[,1]
mat_data <- data.matrix(data[,2:ncol(data)])
rownames(mat_data) <- rnames
my_palette <- colorRampPalette(c("red", "blue"))(n = 299)
col_breaks = c(seq(0,5,length=100), # for red
               seq(5.1,10,length=100), # combo of red & blue
               seq(10.1,15.5,length=100)) # for blue
tiff("heatmap_exp_deg_cluster.tiff",     
     width = 6*300,        # 5 x 300 pixels
     height = 6*300,
     res = 300,            # 300 pixels per inch
     pointsize = 8)        # smaller font size

heatmap.2(mat_data,
          main = "Heatmap", # heat map title
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks, 
          dendrogram="both",     # only draw a row dendrogram
          Colv="T" ,         # turn off column clustering
          lhei = c(1,7)         # Key size width adjustment
)            
dev.off()

# ML (Machine Learning) 
# Model Development for Gene Expression Analysis

#open Postnorm annotated and DEG annotated file
#copy gene names from GEGs annotated file and paste in another sheet
#map these genes in postnorm annotated file , copy paste as values, and remove NA values, duplicate genes
#take only DEGs and its expression values in another csv file name it as Short_table for ML analysis
#move this file to a new folder inside this called ML_Analysis

setwd("C:/Users/gayat/Downloads/GSE63514_RAW/Cervical_Cancer_data/ML_Analysis")

# Load the required packages
# Importing required libraries
#======================#
## load Data
#======================#
df <- read.csv("Short_table.csv", header = T) #if it was workbook then we have to use df<-read_excel("short_table.xlsx", sheet = 1)
head(df)

#============================================#
# Data Pre-processing |   EDA(exploratory data analysis) 
#============================================#

{  # Check the structure of the dataset
  
  str(df)
  
  dim(df)
  
  
  # remove rows that contain NA values
  df <- df[complete.cases(df), ]
  head(df)
  dim(df)
  
  
  #Calculate Mean of duplicate genes
  x <- df
  x <- data.frame(x)
  x <- do.call(rbind,lapply(lapply(split(x,x$DEGs),`[`,2:ncol(x)),colMeans))
  dim(x)
  
  #Convert rownames as a 1st column with header Symbols -> which became rownames after previous operation
  library(tibble) # from tidyverse
  x <- data.frame(x)
  x <- tibble::rownames_to_column(x, var="Symbols")
  head(x)
  dim(x) 
  
  df <- x
  
  # Transpose table 
  install.packages("sjmisc")
  library(sjmisc)  
  df_t <- rotate_df(df, cn=T)
  Symbols <- colnames(df[-1])
  df_t <- cbind(Symbols, df_t)
  write.csv(df_t, "transposed_table.csv", row.names=F)
  
  df_t <- read.csv("transposed_table.csv", h=T)
  dim(df_t)
  df_t[1]
  
  # Healthy_1.CEL to N/T
  #df_t[,1] <- gsub("_.*$", "", df_t[,1])#if the files are present like normal_1-s1.CEL
  #df_t[1]
  
  #My files are like GSM1234_NORMAL.CEL
  #hence i need to extract values between _ and . i.e Normal or Cancer
  #sub() → Substitutes the first match of a regex.
  
  #"\\1" → Returns only the captured group (i.e., what's between _ and .CEL)
  #.*_ → Matches everything up to the last underscore.
  #(.*?) → Captures the shortest string between _ and .CEL → e.g., "Normal" or "Cancer"
  #\\.CEL → Matches .CEL (the . is escaped with \\)
  
  
  df_t$Symbols <- sub(".*_(.*?)\\.CEL", "\\1", df_t$Symbols)
  df_t[1]
  
  # convert Healthy/DCIS from char to factor
  df_t[1] <- factor(df_t$Symbols)
  str(df_t)
  
  # df_t -> df
  df <- df_t # df is df_NT
  
}

# view transformed data
str(df)
df <- data.frame(df)
head(df)

#===============================================#
# Step 3.  Visualize Dataset - Figures - Plots 
#===============================================#
#####  Box and Whisker Plots  ##### 
# Given that the input variables are numeric, we can create box and whisker plots of each
png("box_and_whisker_plots.png")
par(mfrow=c(2,4))
for(i in 2:9) {
  boxplot(x[,i], main=names(df)[i], col="blue")
}
dev.off()


#####  Sample matrix  ##### 
library(ggplot2)
library(caret) 
# split input and output
x <- df[,2:ncol(df)]  # x -  inputs attributes
y <- df[,1]     # y -  outputs attributes
y<- as.factor(y)
plot(y, col="blue")


#======================#
# Step 2. Data Splitting
#======================#
# create 60%/40% for training and testing dataset
library(caret)
set.seed(101)
split <- createDataPartition(df$Symbols, p=0.60, list=FALSE)  # Return the row indices as a matrix/vector, not as a list.
train <- df[split,]
test <- df[-split,]

# dimensions of dataset, train, test
dim(df)
dim(train)
dim(test)


#set cross-validation control for training
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"

head(df_t)


#=========================#
# Build ML Models   
#=========================#

# ML (Machine Learning) 
# Model Development for Gene Expression Analysis - [Part 2]


# Install and Load the required packages
install.packages("caret")
library("caret")


# 1... kNN(k-Nearest Neighbor) - [Model 1]
#-------------------------------------------
set.seed(7)
fit.knn <- train(Symbols~., 
                 data=train, 
                 method="knn", 
                 metric=metric, 
                 trControl=control)
fit.knn


install.packages("cowplot")     # Only once
library(cowplot)                # Load every time you use plot_grid()

# check important variables
varImp(fit.knn)
p1 <- plot(varImp(fit.knn), top = 30, main="kNN")
p2 <- plot(fit.knn, main="kNN")
plot_grid(p1, p2)

# make predictions using trained model on new/test
pred.knn <- predict(fit.knn, newdata = test)

# Model Evaluation
confusionMatrix(pred.knn, test$Symbols, positive = "Cancer")
c1 <- confusionMatrix(pred.knn, test$Symbols, positive = "Cancer")$table

#========================#
# plot Confusion Matrix
#========================#
library(ggplot2)
library(dplyr)
table <- data.frame(c1)
plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "high", "low")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

# fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups (see dplyr code above as well as original confusion matrix for comparison)
ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = Freq)) + # alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(high = "#009194", low="#FF9966")) +
  #scale_fill_gradient(low="white", high="#009194") +
  theme_bw() +
  xlim(rev(levels(table$Reference)))



# 2.SVM model     [model-2]
#-----------------------------------
# For bioinformatics tasks like gene expression classification:
# Use svmLinear when classes are clearly separable (e.g., PCA shows clusters),
# Use svmRadial for complex, non-linear patterns (common in omics data).

install.packages("kernlab")
library(kernlab)
set.seed(101)
fit.svm <- train(Symbols~., data=train, method="svmRadial", metric=metric, trControl=control)
fit.svm

# Feature Importance
varImp(fit.svm)

p1 <- plot(varImp(fit.svm), top = 30, main="Support Vector Machines with Radial Basis")
p2 <- plot(fit.svm)
plot_grid(p1, p2)


# Make Predictions # Confusion Matrix
pred.svm <- predict(fit.svm, newdata = test)

# Model Evaluation
confusionMatrix(pred.svm, test$Symbols, positive = "Cancer")
c1 <- confusionMatrix(pred.svm, test$Symbols, positive = "Cancer")$table



library(ggplot2)
library(dplyr)
table <- data.frame(c1)
plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "high", "low")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

# fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups (see dplyr code above as well as original confusion matrix for comparison)
ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = Freq)) + # alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(high = "#009194", low="#FF9966")) +
  #scale_fill_gradient(low="white", high="#009194") +
  theme_bw() +
  xlim(rev(levels(table$Reference)))




# 3 .Random Forest Model - [Model 3]
#----------------------------------------
install.packages("randomForest")
library(randomForest)

set.seed(123)
fit.rf <- train(Symbols~.,
                data=train,
                method="rf",
                metric=metric,
                trControl=control)
fit.rf

# view important genes
varImp(fit.rf)



# visualize the important genes
plot(varImp(fit.rf), top = 30)

# make predictions using trained model on new/test
pred.rf <- predict(fit.rf, newdata = test)

# Model Evaluation
confusionMatrix(pred.rf, test$Symbols, positive = "Cancer")
c1 <- confusionMatrix(pred.rf, test$Symbols, positive = "Cancer")$table



#========================#
# plot Confusion Matrix
#========================#
library(ggplot2)
library(dplyr)
table <- data.frame(c1)
plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "high", "low")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))

# fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups (see dplyr code above as well as original confusion matrix for comparison)
ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = Freq)) + # alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(high = "#009194", low="#FF9966")) +
  #scale_fill_gradient(low="white", high="#009194") +
  theme_bw() +
  xlim(rev(levels(table$Reference)))



#generating models comaprision plot
acc_knn <- confusionMatrix(pred.knn, test$Symbols)$overall["Accuracy"]
acc_rf  <- confusionMatrix(pred.rf, test$Symbols)$overall["Accuracy"]
acc_svm <- confusionMatrix(pred.svm, test$Symbols)$overall["Accuracy"]

library(ggplot2)

accuracy_df <- data.frame(
  Model = c("kNN", "Random Forest", "SVM"),
  Accuracy = c(acc_knn, acc_rf, acc_svm)  # ✅ properly closed
)

p <- ggplot(accuracy_df, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity", width = 0.6) +
  ylim(0, 1) +
  geom_text(aes(label = round(Accuracy, 3)), vjust = -0.5, size = 5) +
  ggtitle("Model Accuracy Comparison") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )


# Save the barplot to PNG
ggsave("model_accuracy_comparison.png", plot = p, width = 6, height = 4, dpi = 300)

























