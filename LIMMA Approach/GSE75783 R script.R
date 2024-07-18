# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE75783", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13607", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "00XXXXXXXX111111"
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

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
tT <- topTable(fit2, adjust="fdr", sort.by="B", number = 62976)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC","GeneName","SEQUENCE","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


ID <- fit2[["genes"]][["ID"]]

Gene.symbol <- fit2[["genes"]][["GeneName"]]
p <- fit2[["p.value"]]

logFC <- fit2[["coefficients"]]

s2.post <- fit2[["s2.post"]]

GSE75783 <- cbind(ID,Gene.symbol,p,logFC,s2.post)
colnames(GSE75783) <- c("ID","Gene.symbol","p-value","logFC","s2.post")

GSE75783 <- as.data.frame(GSE75783)

GSE75783_1 <- GSE75783[complete.cases(GSE75783$logFC),]
colSums(!is.na(GSE75783_1))

# # removing missing id 
# library(dplyr)
# 
# GSE75783_2 <- GSE75783_1[!(is.na(GSE75783_1$gene.ID)|GSE75783_1$gene.ID==""),]

write.csv(GSE75783_1, "GSE75783.csv")

