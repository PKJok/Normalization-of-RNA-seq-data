# 📊 Normalization and Visualization in edgeR (TMM)

This repository contains the workflow and code used for **different normalization methods** on the RNA-Seq data of *Capsicum annuum L.* subjected to various abiotic stresses (heat, cold, salt, and osmotic stress).

## TMM (Trimmed Mean of M-values) Normalization ##
A statistical method to normalize RNA-seq data to account for differences in library sizes and RNA composition between samples.
**Core Idea:** It assumes that the majority of genes are not differentially expressed (DE) between the samples being compared.
**How it works:** 
1. It selects a sample as a reference (often the one with the median library size).
2. For each other sample ("test" sample), it calculates gene-wise log expression ratios (M-values) and expression-level abundances (A-values).
3. It trims away the top and bottom 30% of the M-values (by default) and the genes with very high or low expression (A-values). This removes potential outliers and highly variable genes that are likely to be true DE genes.
4. It calculates a weighted mean of the remaining M-values. This final value is the TMM factor for that test sample.
**Result:** Each sample gets a scaling factor used to adjust its library size for downstream analyses.

### Why We Need It ###
**1. Removes Technical Bias:** 
a. Corrects for differences in total RNA output (library size) between samples. A sample with 50 million reads shouldn't be directly compared to one with 20 million reads.
    
**2. Corrects for Composition Bias:** 
a. Accounts for situations where a few very highly expressed genes in one sample "use up" a large portion of the sequencing reads, making all other genes in that sample appear under-represented.
b. Enables Accurate Comparison: Allows for a fair comparison of gene expression counts between samples, which is the fundamental goal of most RNA-seq studies (e.g., treated vs. control).
c. Improves Downstream Analysis: Essential for obtaining reliable results in Differential Expression (DE) analysis. Without it, you could falsely identify genes as DE due to technical biases rather than         biological truth.
    
**3. Robust and Reliable:** 
a. The trimming step makes it resistant to outliers and a high proportion of differentially expressed genes, which is why it's a popular and trusted method.

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
