# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE37326", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1708", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "000111111111111111111111"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("T","C"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
# tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- topTable(fit2, adjust="fdr",  number=43534)

# S2 <- fit2[["s2.post"]]
# write.csv(S2, "var-GSE37326.csv")


ID <- fit2[["genes"]][["ID"]]
gene.ID <- fit2[["genes"]][["Gene.ID"]]
Gene.symbol <- fit2[["genes"]][["Gene.symbol"]]
p <- fit2[["p.value"]]

logFC <- fit2[["coefficients"]]

s2.post <- fit2[["s2.post"]]

GSE37326 <- cbind(ID,gene.ID,Gene.symbol,p,logFC,s2.post)
colnames(GSE37326) <- c("ID","gene.ID","Gene.symbol","p-value","logFC","s2.post")

GSE37326 <- as.data.frame(GSE37326)

GSE37326_1 <- GSE37326[complete.cases(GSE37326$logFC),]
colSums(!is.na(GSE37326_1))

# removing missing id 
library(dplyr)

GSE37326_2 <- GSE37326_1[!(is.na(GSE37326_1$gene.ID)|GSE37326_1$gene.ID==""),]

write.csv(GSE37326_2, "GSE37326.csv")

