# PCA on Gene Expression Data (GSE5325)

This project performs **Principal Component Analysis (PCA)** on breast cancer gene expression data from the **GSE5325** dataset obtained from the Gene Expression Omnibus (GEO).

The goal is to reproduce Figure 1 from the Nature Primer paper by analyzing the expression patterns of **XBP1** and **GATA3** genes and studying the separation between **ER+** and **ER−** breast cancer samples using PCA.

---

# Dataset

Dataset: **GSE5325**

Files used:

```text
filtered.tsv.gz     → Gene expression matrix
class.tsv           → Sample labels (ER+ / ER−)
columns.tsv.gz      → Gene ID mapping
```

---

# Biological Background

Breast cancer samples are classified into:

* **ER+ (Estrogen Receptor Positive)**
* **ER− (Estrogen Receptor Negative)**

Gene expression patterns help distinguish these classes.

This project specifically studies:

* **XBP1**
* **GATA3**

which are important biomarkers in breast cancer biology.

---

# Objectives

The project reproduces the following analyses:

## Figure a

Scatter plot of:

* XBP1 expression
* GATA3 expression

with coloring based on ER status.

---

## Figure b

Visualization of PCA axes (PC1 and PC2) on the scatter plot.

---

## Figure c

Projection of samples onto the first principal component (PC1).

---

## Figure d

Scree plot showing variance explained by principal components.

---

## Figure e

Projection of samples in PCA space (PC1 vs PC2).

---

## Figure f

PCA clustering visualization of ER+ and ER− samples.

---

# Methods Used

## Data Processing

* Loading compressed TSV datasets
* Extracting gene expression values
* Separating samples by class labels

---

## PCA Workflow

The PCA analysis involves:

1. Mean centering
2. Covariance matrix computation
3. Eigendecomposition
4. Projection onto principal components

Mathematically:

```text
Z = X_centered V
```

where:

* `X_centered` = centered data matrix
* `V` = eigenvector matrix
* `Z` = PCA projections

---

# Libraries Used

```python
numpy
pandas
matplotlib
```

---

# Running the Project

## Install dependencies

```bash
pip install pandas matplotlib numpy
```

---

## Run the script

```bash
python pca.py
```

or on Mac:

```bash
python3 pca.py
```

---

# Output Files

The following figures are generated:

```text
figure_a.png
figure_b.png
figure_c.png
figure_d.png
figure_e.png
figure_f.png
combined_figure.png
```

---

# Results

The PCA analysis shows that:

* ER+ and ER− samples form distinct clusters
* XBP1 and GATA3 expression patterns correlate with cancer subtype
* PC1 captures the major source of variance in the dataset

---

# Project Structure

```text
Bioinfo_Project/
│
├── filtered.tsv.gz
├── class.tsv
├── columns.tsv.gz
├── pca.py
├── README.md
│
├── figure_a.png
├── figure_b.png
├── figure_c.png
├── figure_d.png
├── figure_e.png
├── figure_f.png
└── combined_figure.png
```

---

# Author

Bioinformatics PCA Project using GSE5325 breast cancer gene expression dataset.
