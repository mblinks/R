# Single Cell Sequencing
## Analysis of single-cell RNA-Seq data of Human Glioblastoma Multiforme cells from a male donor

Visualizing the QC

![QC](https://github.com/mblinks/mblinks.github.io/blob/main/assets/vlnPlot.png)

Filtering  cells that have unique feature counts over 2,500 or less than 200 
We filter cells that have >5% mitochondrial counts

![Feature Count](https://github.com/mblinks/mblinks.github.io/blob/main/assets/Feature_Scatter.png)

Identification of highly variable features (feature selection)

Calculating a subset of features that exhibit high cell-to-cell variation in the dataset (i.e., they are highly expressed in some cells, and lowly expressed in others)

![Feature selection](https://github.com/mblinks/mblinks.github.io/blob/main/assets/sd_var1.PNG)
![Feature selection](https://github.com/mblinks/mblinks.github.io/blob/main/assets/sd_var2.PNG)

Linear dimensional reduction

Seurat outputs a list of genes with the most positive and negative loadings, representing modules of genes that exhibit either correlation (or anti-correlation) across single cells in the dataset.

![Linear dimension reduction](https://github.com/mblinks/mblinks.github.io/blob/main/assets/PCA.PNG)

![Linear dimension reduction](https://github.com/mblinks/mblinks.github.io/blob/main/assets/pc1pc2.PNG)

Exploring correlated feature sets

![Correlated feature sets](https://github.com/mblinks/mblinks.github.io/blob/main/assets/pc1heatmap.PNG)

![Correlated feature sets](https://github.com/mblinks/mblinks.github.io/blob/main/assets/PC1-15%20heatmap.PNG) 

Determining the ‘dimensionality’ of the dataset

Ranking of principal components based on the percentage of variance using elbow plot

![Elbow plot](https://github.com/mblinks/mblinks.github.io/blob/main/assets/elbowplot.png)

Non-linear dimensional reduction (UMAP/tSNE)

Cells are grouped within graph-based clusters determined should co-localize on these dimension reduction plots.

![UMAP](https://github.com/mblinks/mblinks.github.io/blob/main/assets/Umap.png)

Finding differentially expressed features (cluster biomarkers)

![Differential Expression](https://github.com/mblinks/mblinks.github.io/blob/main/assets/express1.png)
![Differential Expression](https://github.com/mblinks/mblinks.github.io/blob/main/assets/express2.png)

Raw counts

![Raw count](https://github.com/mblinks/mblinks.github.io/blob/main/assets/express2.png)

![Feature plot](https://github.com/mblinks/mblinks.github.io/blob/main/assets/featureplot.png)

An expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

![Do Heat](https://github.com/mblinks/mblinks.github.io/blob/main/assets/DoHeat.png)

---
