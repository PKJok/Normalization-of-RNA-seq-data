#######################################################
############# Different Normalization #################
#######################################################


setwd("C:/Users/ASUS/OneDrive/Desktop/Chili")
library(GEOquery)
library(stringr)
library(edgeR)
library(tidyverse)
library(limma)

# load the count_data
count_data<- as.data.frame(read.table("GSE132824_Abiotic_RNA-Seq.Readcount.txt"))
col_data<- read.csv("metadata.csv")
col_data<- col_data%>%
  column_to_rownames(var = "Samples")
dim(count_data)

# check the sequencing depth
seq_depth<- colSums(count_data)
seq_depth<- seq_depth[order(seq_depth)]

barplot(seq_depth)
seq_depth%>%
  as.data.frame()%>%
  rownames_to_column(var = "Samples")%>%
  arrange(seq_depth)%>%
  mutate(Samples= factor(Samples,levels = Samples))%>%
  ggplot(aes(Samples, seq_depth))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
  

#### DGE list of edgeR ####
## calculate the sequencing depth/ library sizes
## store the count_data and normalization factors
y<- DGEList(count_data)
dim(y)
y$samples

# adding group info
group<- paste(col_data$stress, col_data$time, sep = ".")
group<- factor(group)
group

# add the group information to DGEList
y$samples$group<- group
y$samples

## filter lowly expressed genes ###
## filtering out genes that are too lowly expressed to be statistically realiable
## ensures that only genes with enough counts in enough replicates remains 
keep<- filterByExpr(y, group = group)
y<- y[keep, , keep.lib.sizes=FALSE]
dim(y)

### type 1: no normalization ####
## this negative control has only CPM normalization and log transformation

no_normalize<- calcNormFactors(y, method = "none")
logcount<- cpm(no_normalize, log = TRUE)

# boxplot
boxplot(logcount)
abline(h=median(logcount),col="red")

logcount%>%
  as.data.frame()%>%
  gather(key = "Samples",value = "Counts")%>%
  ggplot(aes(Samples, Counts))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90))

# MDS plot
plotMDS(no_normalize, col=as.numeric(group))


## TMM method (Trimmed Mean of M-values) used for correcting the composition bais
## 1st : reference sample is choosen
##       reference sample = one which is near to the average mean of all the samples
## 2nd: genes are selected to create scaling factor
##      the genes which are not outliers are selected 
##       (outliers can be known by observing log fold differences)
## 3rd: Log fold values are trimmed 
##      from top 30% and bottom 30% and mean of logs are trimmed from top 5% and bottom 5%.
## 4rth: The genes which are overlapping from both trimmed log fold values and mean of logs are used to calculated scaling factor.
## 5th : original read counts are divided by this scaling factor

TMM<- calcNormFactors(y, method ="TMM")
logcount<- cpm(TMM, log = TRUE)

# boxplot
boxplot(logcount, xlab="", ylab="TMM Normalized Counts")
abline(h=median(logcount), col="red")


# MDS plot
plotMDS(no_normalize, col=as.numeric(group))


## checking PCA
arr.col<- rownames(col_data)

logcount<- logcount%>%
  as.data.frame()%>%
  select(all_of(arr.col))

# forming dataframes 
heat_moc<- logcount%>%
  select(1:18,34:48)
colnames(heat_moc)

cold_moc<- logcount%>%
  select(1:33)
head(cold_moc)

man_moc<- log_count%>%
  select(all_of(mock_sample), all_of(man_sample))
head(man_moc)

NaCl_moc<- log_count%>%
  select(all_of(mock_sample), all_of(NaCl_sample))
head(NaCl_moc)

# running pca
# heat and mock
pca_heat <- prcomp(t(cold_moc))
plot(pca_heat$x[,1], pca_heat$x[,2])






