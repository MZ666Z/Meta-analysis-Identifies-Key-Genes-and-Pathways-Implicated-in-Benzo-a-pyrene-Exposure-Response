# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE75783", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13607", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

a <-  gset@assayData[["exprs"]]

a <- as.data.frame(a)

sum(is.na(a))


id <- gset@featureData@data[["GeneName"]]
data <- cbind(id,a)

data1 <- data

# removing missing id 
library(dplyr)

data2 <- data1[!(is.na(data1$id)|data1$id==""),]


# Check for duplicate genes 
dupli <- subset(data2,duplicated(id))

# calculate columnn mean for duplicate genes

data3 <- aggregate(data2,by=list(id=data2$id),data=data2,FUN=mean, na.rm=TRUE)

data3 <- select(data3,-2)

data3[data3== "NaN"] <- NA

# Check for duplicate genes again
dupli1 <- subset(data3,duplicated(id))




library(tidyr)

data_long <- gather(data3, GSM, log2measure, GSM1967620, GSM1967621, 
                    GSM1967630, GSM1967631, GSM1967632, GSM1967633,
                    GSM1967634,GSM1967635)


library(dplyr)

data_long1 <- data_long %>%
  mutate(group = case_when(GSM== 'GSM1967620'|GSM== 'GSM1967621' ~ 'BaP',
                           GSM== 'GSM1967630'|GSM== 'GSM1967631'| GSM== 'GSM1967632'
                           |GSM== 'GSM1967633'|GSM== 'GSM1967634'|GSM== 'GSM1967635'~ 'DMSO'))


data_long3 <- subset(data_long1, select = -c(2:9))

str(data_long3)

gene_list <- data3$id



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



# # before delete duplicate genes, using two genes to test#  zzz3 
# 
# data_long4 <- data_long3[data_long3$id=="ZZZ3",]
# data_long4 <- data_long3[data_long3$id=="ZNF362",]
# fit <- lm(log2measure~group,data=data_long4)
# summary(fit)

# # test normality
# hist(data_long4$log2measure)
# shapiro.test(data_long4$log2measure)


# change direction of groupDMSO
groupDMSO_c <- as.data.frame(groupDMSO)

library(dplyr)

groupDMSO_c2 <- groupDMSO_c %>% 
  mutate_if(is.numeric, funs(. * -1))


#combine variables
newmatrix_2 <- cbind(gene_list,groupDMSO_c2,t_value,p_value,SE)

GSE75783 <- as.data.frame(newmatrix_2)

GSE75783_1 <- GSE75783[complete.cases(GSE75783$groupDMSO),]
colSums(!is.na(GSE75783_1))

GSE75783_1[GSE75783_1== "NaN"] <- NA

write.csv(GSE75783_1, "gse75783-lm.csv")
