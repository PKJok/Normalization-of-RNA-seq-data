# 📊 Normalization and Visualization in edgeR (TMM)

This repository contains the workflow and code used for **different normalization methods** on the RNA-Seq data of *Capsicum annuum L.* subjected to various abiotic stresses (heat, cold, salt, and osmotic stress).

## 📂 Data Sources

* **Raw RNA-seq counts:** [GSE132824\_Abiotic\_RNA-Seq.Readcount.txt](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132824)
* **Metadata file:** `metadata.csv` (contains sample information like stress and time points)

## 🔹 Workflow in edgeR for Normalization
### 1. Create DGEList and add group information

```r
y <- DGEList(count_data)
group <- factor(paste(col_data$stress, col_data$time, sep = "."))
y$samples$group <- group
```

### 2. Filter lowly expressed genes

```r
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes=FALSE]
```

### 3. Normalization Methods

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


```

---

Feel free to contribute or open issues for clarifications.
