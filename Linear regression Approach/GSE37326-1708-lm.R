# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis

library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE37326", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1708", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# log2 transformation
ex <- exprs(gset)

a <-  gset@assayData[["exprs"]]


a <- as.data.frame(a)

sum(is.na(a))


id <- gset@featureData@data[["Gene symbol"]]
data <- cbind(id,a)
data1 <- data

data1[data1== "NaN"] <- NA


# removing missing id 
library(dplyr)

data1 <- data1[!(is.na(data1$id)|data1$id==""),]

# removing genes if missing expression >50%

data1$count.na <- apply(data1[2:25], 1, function(x) sum(is.na(x)))

data2 <- subset(data1, count.na< 12) 

data2 <- subset(data2, select = -c(count.na))

# Check for duplicate genes 
dupli <- subset(data2,duplicated(id))


# calculate columnn mean for duplicate genes

data3 <- aggregate(data2,by=list(id=data2$id),data=data2,FUN=mean, na.rm=TRUE)

data3 <- select(data3,-2)

data3[data3== "NaN"] <- NA

# Check for duplicate genes again
dupli1 <- subset(data3,duplicated(id))




library(tidyr)

data_long <- gather(data3, GSM, log2measure, GSM914178, GSM914179, GSM914180, GSM915413, GSM915414, GSM915415,
                    GSM915416,GSM915417,GSM915418,GSM915419,GSM915420,GSM915421,
                    GSM915422,GSM915423,GSM915424,GSM915425,GSM915426,GSM915427,
                    GSM915429,GSM915430,GSM915432,GSM915438,GSM915439,GSM915440)


library(dplyr)

data_long1 <- data_long %>%
  mutate(group = case_when(GSM== 'GSM914178'|GSM== 'GSM914179'| GSM== 'GSM914180' ~ 'BaP',
                           GSM== 'GSM915413'|GSM== 'GSM915414'|GSM== 'GSM915415'|GSM== 'GSM915416'
                           |GSM== 'GSM915417'|GSM== 'GSM915418'|GSM== 'GSM915419'
                           |GSM== 'GSM915420'|GSM== 'GSM915421'|GSM== 'GSM915422'
                           |GSM== 'GSM915423'|GSM== 'GSM915424'|GSM== 'GSM915425'
                           |GSM== 'GSM915426'|GSM== 'GSM915427'|GSM== 'GSM915429'
                           |GSM== 'GSM915430'|GSM== 'GSM915432'|GSM== 'GSM915438'
                           |GSM== 'GSM915439'|GSM== 'GSM915440'~ 'DMSO'))



data_long3 <- data_long1

str(data_long3)

gene_list <- as.matrix(data3$id)



# fit lm model
groupDMSO <- c()
t_value <- c()
p_value <- c()
SE <- c()
genesymbol <- c()

{for (i in 1:length(gene_list)) 
  try(  { data_long4 <- data_long3[data_long3$id==gene_list[i],]
  
  
  fit4 <- lm(log2measure~group,data=data_long4)
  
  groupDMSO[i] <-  summary(fit4)[["coefficients"]][2,1]
  SE[i] <- summary(fit4)[["coefficients"]][2,2]
  t_value[i] <-summary(fit4)[["coefficients"]][2,3]
  p_value[i] <- summary(fit4)[["coefficients"]][2,4]
  genesymbol[i] <- gene_list[i]
  }
  
  , silent = T
  )
}



# data_long3$log2measure[data_long3$log2measure == "NaN"] <- NA 
# data_long4 <- data_long3[data_long3$id=="TMPRSS2",]
# fit4 <- lm(log2measure~group,data=data_long4)
# summary(fit4) 
# summary(fit4)[["coefficients"]][2,4]



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

GSE37326 <- as.data.frame(newmatrix_2)

GSE37326_1 <- GSE37326[complete.cases(GSE37326$groupDMSO),]
colSums(!is.na(GSE37326_1))

write.csv(GSE37326_1, "gse37326-lm.csv")
