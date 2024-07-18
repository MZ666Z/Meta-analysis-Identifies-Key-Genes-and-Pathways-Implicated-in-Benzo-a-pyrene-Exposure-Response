
library(readxl)
library(readr)

my_data1 <- read_excel("C:\\Users\\zhu66\\OneDrive - University of Oklahoma\\Desktop\\GEO\\geo database2\\37326\\GSE37326-5.22.2022.xlsx",sheet = "Sheet1")
my_data2 <- read_excel("C:\\Users\\zhu66\\OneDrive - University of Oklahoma\\Desktop\\GEO\\geo database2\\117527\\gse117527-link-R-5.22.2022.xlsx",sheet = "Sheet1")
my_data3 <- read_excel("C:\\Users\\zhu66\\OneDrive - University of Oklahoma\\Desktop\\GEO\\geo database2\\36244-2\\gse36244-link-R-6.28.2022.xlsx",sheet = "Sheet1")
my_data4 <- read_excel("C:\\Users\\zhu66\\OneDrive - University of Oklahoma\\Desktop\\GEO\\geo database2\\75783\\GSE75783-7.5.2022.xlsx",sheet = "Sheet1")

library(sqldf)


my_data1_1 <-  sqldf('select *  from my_data1 group by "Gene.symbol" having min("P.Value") order by "Gene.symbol" ')
my_data2_1 <-  sqldf('select *  from my_data2 group by "Gene.symbol" having min("P.Value") order by "Gene.symbol" ')
my_data3_1 <-  sqldf('select *  from my_data3 group by "Gene.symbol" having min("P.Value") order by "Gene.symbol" ')
my_data4_1 <-  sqldf('select *  from my_data4 group by "Gene.symbol" having min("P.Value") order by "Gene.symbol" ')





# # check if the gene is the min
# c <- sqldf('select * from my_data1 group by "P.Value" ')
# 
# d <- sqldf('select * from c where "Gene.symbol"="ZNF121" ')
# 
# e <- sqldf('select * from my_data1 where "Gene.symbol"="SDC1" ')
# 
# f <- sqldf('select * ,count(1) as CNT from my_data1 group by "Gene.symbol"  ')
# 
# g <- sqldf('select * from my_data1 group by "Gene.symbol" ')
# H <- sqldf('select * from g where "Gene.symbol"="SDC1" ')





match(my_data1_1$Gene.symbol,my_data2_1$Gene.symbol)


my_data1_1$Gene.symbol2=my_data2_1$Gene.symbol[match(my_data1_1$Gene.symbol,my_data2_1$Gene.symbol)]


my_data1_1$P.Value2=my_data2_1$P.Value[match(my_data1_1$Gene.symbol,my_data2_1$Gene.symbol)]

my_data1_1$logFC2=my_data2_1$logFC[match(my_data1_1$Gene.symbol,my_data2_1$Gene.symbol)]
my_data1_1$variance2=my_data2_1$variance[match(my_data1_1$Gene.symbol,my_data2_1$Gene.symbol)]

my_data1_1=na.omit(my_data1_1)

library(xlsx)

match(my_data1_1$Gene.symbol,my_data3_1$Gene.symbol)

my_data1_1$Gene.symbol3=my_data3_1$Gene.symbol[match(my_data1_1$Gene.symbol,my_data3_1$Gene.symbol)]


my_data1_1$P.Value3=my_data3_1$P.Value[match(my_data1_1$Gene.symbol,my_data3_1$Gene.symbol)]
my_data1_1$logFC3=my_data3_1$logFC[match(my_data1_1$Gene.symbol,my_data3_1$Gene.symbol)]
my_data1_1$variance3=my_data3_1$variance[match(my_data1_1$Gene.symbol,my_data3_1$Gene.symbol)]
my_data1_1=na.omit(my_data1_1)


match(my_data1_1$Gene.symbol,my_data4_1$Gene.symbol)

my_data1_1$Gene.symbol4=my_data4_1$Gene.symbol[match(my_data1_1$Gene.symbol,my_data4_1$Gene.symbol)]


my_data1_1$P.Value4=my_data4_1$P.Value[match(my_data1_1$Gene.symbol,my_data4_1$Gene.symbol)]
my_data1_1$logFC4=my_data4_1$logFC[match(my_data1_1$Gene.symbol,my_data4_1$Gene.symbol)]
my_data1_1$variance4=my_data4_1$variance[match(my_data1_1$Gene.symbol,my_data4_1$Gene.symbol)]
my_data1_1=na.omit(my_data1_1)




write.xlsx(my_data1_1,file = 'linkgen6.xlsx', sheetName = "Sheet1", append = F)



