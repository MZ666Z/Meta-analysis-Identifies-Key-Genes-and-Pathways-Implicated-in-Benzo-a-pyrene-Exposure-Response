Description
This project involves analyzing gene expression data from four different datasets using two approaches: LIMMA Approach and Linear Regression Approach. The goal is to identify overlapping significant genes across these datasets.
Folder Structure
•	LIMMA Approach
o	Use scripts from each study folder to identify genes.
o	Utilize linkgenname_selected-o.R script to find overlapping genes across datasets.
o	Record adjusted p-value, log fold change, and standard error.
o	Annotate genes using DAVID Online.
o	Use MetaVolcanoR-combine-logFC-o.R script to calculate combined log fold change and corresponding p-value.
o	Use MetaVolcanoR-o.R script to calculate combined p-value.
•	Linear Regression Approach
o	Use scripts from each study folder to identify genes.
o	Utilize linkgenname_selected-s.R script to find overlapping genes across datasets.
o	Record adjusted p-value, log fold change, and standard error.
o	Annotate genes using DAVID Online.
o	Use MetaVolcanoR-combine-logFC-s.R script to calculate combined log fold change and corresponding p-value.
o	Use MetaVolcanoR-s.R script to calculate combined p-value.
Excel Analysis
•	Combine results from both approaches in Excel.
•	Filter to identify two overlapping significant genes with the same direction of effect.
