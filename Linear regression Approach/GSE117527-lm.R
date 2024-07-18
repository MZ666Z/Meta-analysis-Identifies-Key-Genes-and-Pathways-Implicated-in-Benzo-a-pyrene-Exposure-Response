# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE117527", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL25336", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


# log2 transformation
ex <- exprs(gset)
# qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# LogC <- (qx[5] > 100) ||
#   (qx[6]-qx[1] > 50 && qx[2] > 0)
# if (LogC) { ex[which(ex <= 0)] <- NaN
# exprs(gset) <- log2(ex) }

a <-  gset@assayData[["exprs"]]

a <- as.data.frame(a)

sum(is.na(a))

id <- gset@featureData@data[["ORF"]]
data <- cbind(id,a)


# #define Min-Max normalization function
# min_max_norm <- function(x) {
#   (x - min(x,na.rm = TRUE)) / (max(x,na.rm = TRUE) - min(x,na.rm = TRUE))
# }
# 
# #apply Min-Max normalization to dataset
# data1 <- as.data.frame(lapply(data[2:9], min_max_norm))
# 
# data1 <- cbind(id,data1)
# 
# 

# Anti-log transformation
data1 <- data

# data1[,c(2:9)] <- 2^(data[,c(2:9)])

# # log2measure <- log2(data1)
# # data_long3 <- cbind(data_long1,log2measure)
# 
# # data_long3$log2measure[data_long3$log2measure == "NaN"] <- ""
# 
# data1[data1== "NaN"] <- NA


# removing missing id 
library(dplyr)

data2 <- data1[!(is.na(data1$id)|data1$id==""),]

# Check for duplicate genes 
dupli <- subset(data2,duplicated(id))
#no duplicate genes in this dataset


library(tidyr)

data_long <- gather(data2, GSM, log2measure, GSM3302480, GSM3302481, GSM3302482, GSM3302483, GSM3302476, GSM3302477,
                    GSM3302478,GSM3302479)


library(dplyr)

data_long1 <- data_long %>%
  mutate(group = case_when(GSM== 'GSM3302480'|GSM== 'GSM3302481'|GSM== 'GSM3302482'|GSM== 'GSM3302483' ~ 'BaP',
                          GSM== 'GSM3302476'|GSM== 'GSM3302477'| GSM== 'GSM3302478'|GSM== 'GSM3302479'~ 'DMSO'))




data_long3 <- data_long1
str(data_long3)

gene_list <- data2$id


# fit lm model
groupDMSO <- c()
t_value <- c()
p_value <- c()
SE <- c()

{for (i in 1:length(gene_list)) 
  try(  { data_long4 <- data_long3[data_long3$id==gene_list[i],]
  
  
  fit4 <- lm(log2measure~group,data=data_long4)
  
  groupDMSO[i] <-  summary(fit4)[["coefficients"]][2,1]
  SE[i] <- summary(fit4)[["coefficients"]][2,2]
  t_value[i] <-summary(fit4)[["coefficients"]][2,3]
  p_value[i] <- summary(fit4)[["coefficients"]][2,4]
  
  }
  
  , silent = T
  )
}

# data_long4 <- data_long3[data_long3$id=="MAP3K6",]
# fit4 <- lm(log2measure~group,data=data_long4)
# summary(fit4)
# summary(fit4)[["coefficients"]][2,1]



# change direction of groupDMSO
groupDMSO_c <- as.data.frame(groupDMSO)

library(dplyr)

groupDMSO_c2 <- groupDMSO_c %>% 
  mutate_if(is.numeric, funs(. * -1))


# adjusted p value
# library(qvalue)
# 
# adjustp <- qvalue(p = p_value)
# 
# adjusted.P <- adjustp[["qvalues"]]


#combine variables
newmatrix_2 <- cbind(gene_list,groupDMSO_c2,t_value,p_value,SE)

GSE117527 <- as.data.frame(newmatrix_2)

GSE117527_1 <- GSE117527[complete.cases(GSE117527$groupDMSO),]
colSums(!is.na(GSE117527_1)) 

write.csv(GSE117527_1, "gse117527-lm.csv")

