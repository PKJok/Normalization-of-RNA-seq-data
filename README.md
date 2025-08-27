# 📊 Normalization and Visualization in edgeR (TMM)

This repository contains the workflow and code used for **different normalization methods** on the RNA-Seq data of *Capsicum annuum L.* subjected to various abiotic stresses (heat, cold, salt, and osmotic stress).

## 📂 Data Sources

* **Raw RNA-seq counts:** [GSE132824\_Abiotic\_RNA-Seq.Readcount.txt](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132824)
* **Metadata file:** `metadata.csv` (contains sample information like stress and time points)

## 🛠 Software and Libraries

The analysis uses **R** and the following packages:

```r
library(GEOquery)
library(stringr)
library(edgeR)
library(tidyverse)
library(limma)
```

## 🔹 Workflow

### 1. Load count data and metadata

```r
count_data <- read.table("GSE132824_Abiotic_RNA-Seq.Readcount.txt")
col_data <- read.csv("metadata.csv") %>% column_to_rownames(var = "Samples")
```

### 2. Check sequencing depth

```r
seq_depth <- colSums(count_data) %>% sort()
barplot(seq_depth)
```

### 3. Create DGEList and add group information

```r
y <- DGEList(count_data)
group <- factor(paste(col_data$stress, col_data$time, sep = "."))
y$samples$group <- group
```

### 4. Filter lowly expressed genes

```r
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes=FALSE]
```

### 5. Normalization Methods

#### a) No normalization (negative control)

```r
no_normalize <- calcNormFactors(y, method = "none")
logcount <- cpm(no_normalize, log = TRUE)
boxplot(logcount)
plotMDS(no_normalize, col=as.numeric(group))
```

#### b) TMM (Trimmed Mean of M-values)

```r
TMM <- calcNormFactors(y, method = "TMM")
logcount <- cpm(TMM, log = TRUE)
boxplot(logcount, ylab="TMM Normalized Counts")
plotMDS(TMM, col=as.numeric(group))
```

### 6. PCA Analysis

```r
# Organize logcounts according to sample order
logcount <- logcount %>% as.data.frame() %>% select(all_of(rownames(col_data)))

# Subset for each stress group
heat_moc <- logcount %>% select(1:18,34:48)
cold_moc <- logcount %>% select(1:33)

# Run PCA for cold
pca_heat <- prcomp(t(cold_moc))
plot(pca_heat$x[,1], pca_heat$x[,2])
```

## 📌 Notes

* The code includes **visualization of sequencing depth**, **boxplots of normalized counts**, **MDS plots**, and **PCA plots**.
* TMM normalization corrects for composition bias and is recommended for downstream differential expression analysis.
* Ensure all required R packages are installed for the scripts to work properly.

## ⚡ Recommended R Packages

```r
install.packages(c("tidyverse", "stringr"))
BiocManager::install(c("edgeR", "limma", "GEOquery"))
```

## 📂 Project Structure

```
├── GSE132824_Abiotic_RNA-Seq.Readcount.txt  # raw counts
├── metadata.csv                              # sample metadata
├── normalization_analysis.R                  # R script with code
└── README.md                                 # this file
```

---

Feel free to contribute or open issues for clarifications.
