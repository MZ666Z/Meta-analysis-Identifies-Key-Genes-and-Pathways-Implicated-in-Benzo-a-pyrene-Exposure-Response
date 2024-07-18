# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis
library(GEOquery)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE36244", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13695", attr(gset, "names")) else idx <- 1
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
gene_neg <-  a[which(a < 0)]


id <- rownames(a)
data <- cbind(id,a)
data1 <- data


# removing missing id 
library(dplyr)

data1 <- data1[!(is.na(data1$id)|data1$id==""),]

# Check for duplicate genes 
dupli <- subset(data,duplicated(id))
#no duplicate genes in this dataset


library(tidyr)

data_long <- gather(data1, GSM, log2measure, GSM884989, GSM884990, GSM884991, GSM884992, GSM884993, GSM884994,
                    GSM884995,GSM884996)


library(dplyr)

data_long1 <- data_long %>%
  mutate(group = case_when(GSM== 'GSM884989'|GSM== 'GSM884990'| GSM== 'GSM884993'|GSM== 'GSM884994' ~ 'BaP',
                        GSM== 'GSM884991'|GSM== 'GSM884992'|GSM== 'GSM884995'|GSM== 'GSM884996' ~ 'DMSO'))
                           

data_long2 <- data_long1 %>%
mutate(time = case_when(GSM== 'GSM884989' ~ '12h',
                         GSM== 'GSM884990' ~ '12h',
                         GSM== 'GSM884991' ~ '12h',
                         GSM== 'GSM884992' ~ '12h',
                         GSM== 'GSM884993' ~ '24h',
                         GSM== 'GSM884994' ~ '24h',
                         GSM== 'GSM884995' ~ '24h',
                         GSM== 'GSM884996' ~ '24h'
))

data_long2 <- data_long2 %>%
  mutate(newid = case_when(GSM== 'GSM884989' ~ '1',
                          GSM== 'GSM884990' ~ '2',
                          GSM== 'GSM884991' ~ '3',
                          GSM== 'GSM884992' ~ '4',
                          GSM== 'GSM884993' ~ '1',
                          GSM== 'GSM884994' ~ '2',
                          GSM== 'GSM884995' ~ '3',
                          GSM== 'GSM884996' ~ '4'
  ))




data_long3 <-data_long2[order(data_long2$id),]


#log2 transform
# log2measure <- log2(data_long3$measurement)
# data_long3 <- cbind(data_long3,log2measure)
# 
# data_long3$group <-  as.factor(data_long3$group)
# data_long3$time <- as.factor(data_long3$time)


str(data_long3)

gene_list <- rownames(a)

# library(nlme)
# check interaction term
library(nlme)

groupDMSO <- c()
t_value <- c()
p_value <- c()
interaction <- c()

{for (i in 1:length(gene_list)) 
 try(  { data_long4 <- data_long3[data_long3$id==gene_list[i],]
  
  
  fit4 <- lme(log2measure~group*time,random = ~1|newid,data=data_long4)
  
  groupDMSO[i] <-  summary(fit4)$tTable[2,1]
  t_value[i] <-summary(fit4)$tTable[2,4]
  p_value[i] <- summary(fit4)$tTable[2,5]
  interaction[i] <- summary(fit4)$tTable[4,5]
  }
  
  , silent = T
)
}



newmatrix_1 <- cbind(rownames(a),groupDMSO,t_value,p_value,interaction)

data36244_1 <- as.data.frame(newmatrix_1) 

sig_interaction1 <- subset(data36244_1,interaction<0.05)

library(qvalue)

adjustp_1 <- qvalue(p = interaction)

adjusted.P_1 <- adjustp_1[["qvalues"]]

adjusted.P_1 <- as.data.frame(adjusted.P_1)

# check adjusted p-values of interaction term. Drop interaction term if appropriate

sig_interaction2 <- subset(adjusted.P_1,adjusted.P_1<0.05)

write.csv(adjusted.P_1, "adjusted_interaction.csv")




# fit model without time
groupDMSO <- c()
group_p_value <- c()
SE <- c()

{for (i in 1:length(gene_list)) 
  try(  { data_long4 <- data_long3[data_long3$id==gene_list[i],]
  
  
  fit4 <- lm(log2measure~group,data=data_long4)
  
  groupDMSO[i] <-  summary(fit4)$coefficients[2,1]
  group_p_value[i] <- summary(fit4)$coefficients[2,4] 
  SE[i] <- summary(fit4)$coefficients[2,2]
  
  }
  
  , silent = T
  )
}


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
newmatrix_2 <- cbind(rownames(a),groupDMSO_c2,group_p_value,SE)

newmatrix_2 <- newmatrix_2[complete.cases(newmatrix_2$groupDMSO),]
colSums(!is.na(newmatrix_2))

data36244 <- as.data.frame(newmatrix_2) 

write.csv(newmatrix_2, "gse36244lm.csv")
