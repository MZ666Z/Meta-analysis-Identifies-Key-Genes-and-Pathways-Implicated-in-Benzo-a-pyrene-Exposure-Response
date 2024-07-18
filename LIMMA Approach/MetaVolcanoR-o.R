
# unlink("C:\\Users\\zhu66\\OneDrive - University of Oklahoma\\Documents\\R\\win-library\\4.1\\00LOCK", recursive = TRUE)
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("MetaVolcanoR")
library(ggplot2)

library(MetaVolcanoR)
#help("MetaVolcano-class")

library(readxl)
GSE37326 <- read_excel("C:\\Users\\zhu66\\OneDrive - University of Oklahoma\\Desktop\\GEO\\geo database2\\metavolcano\\GSE37326-7.17.2022.xlsx",sheet = "Sheet1")
GSE117527 <- read_excel("C:\\Users\\zhu66\\OneDrive - University of Oklahoma\\Desktop\\GEO\\geo database2\\metavolcano\\GSE117527-7.17.2022.xlsx",sheet = "Sheet1")
GSE36244 <- read_excel("C:\\Users\\zhu66\\OneDrive - University of Oklahoma\\Desktop\\GEO\\geo database2\\metavolcano\\GSE36244-7.17.2022.xlsx",sheet = "Sheet1")
GSE75783 <- read_excel("C:\\Users\\zhu66\\OneDrive - University of Oklahoma\\Desktop\\GEO\\geo database2\\metavolcano\\GSE75783-7.17.2022.xlsx",sheet = "Sheet1")


# combine p-value (p-values using the Fisher method, logFC using mean)
L <- list(as.data.frame(GSE37326),as.data.frame(GSE117527),as.data.frame(GSE36244),as.data.frame(GSE75783))

names(L) <- c('GSE37326','GSE117527','GSE36244','GSE75783')

names(L)

mv1 <- combining_mv(L, pcriteria = "P.Value",
                    foldchangecol = "logFC", genenamecol = "Gene.symbol", geneidcol = NULL,
                    collaps = FALSE,jobname = "MetaVolcano",
                    outputfolder = ".", draw = "HTML",)
mv1@MetaVolcano

# p-value

p <- mv1@metaresult[["metap"]]
sum(p<0.05)

symbol <- mv1@metaresult[["Gene.symbol"]]
p_value <- as.data.frame(cbind(symbol,p))

library(xlsx)
write.xlsx(p_value,file = 'p_metavolcano-7.17.2022.xlsx',sheetName = "Sheet1", append = F)






# adjust p-value
adjusted.P <- p.adjust(p,method = 'bonferroni')
sum(adjusted.P<0.05)


symbol <- mv1@metaresult[["Gene.symbol"]]
p1 <- as.data.frame(cbind(symbol,adjusted.P)) 

library(xlsx)
write.xlsx(p1,file = 'adjustedp_metavolcano-7.17.2022.xlsx',sheetName = "Sheet1", append = F)










# combine effect size by using random effect model.Also get 95%CI

L <- list(as.data.frame(GSE37326),as.data.frame(GSE117527),as.data.frame(GSE36244))

names(L) <- c('GSE37326','GSE117527','GSE36244')
names(L)

mv2 <- rem_mv(L, pcriteria = "P.Value",
       foldchangecol = "logFC", genenamecol = "Gene.symbol", geneidcol = NULL,
       collaps = FALSE, llcol = NULL, rlcol = NULL, vcol = "variance",
       cvar = F, metathr = 0.01, jobname = "MetaVolcano",
       outputfolder = ".", draw = "HTML", ncores = 1)

metaresult_symbol <- mv2@metaresult[["Gene.symbol"]]
metaresult_randomSummary <- mv2@metaresult[["randomSummary"]]
# randomSummary is logFC
metaresult_randomclb <- mv2@metaresult[["randomCi.lb"]]
metaresult_randomcub <-mv2@metaresult[["randomCi.ub"]]
metaresult_randop <-mv2@metaresult[["randomP"]]

metaresult <- cbind(metaresult_symbol,metaresult_randomSummary,metaresult_randomclb,
                    metaresult_randomcub,metaresult_randop)

write.xlsx(metaresult,file = 'metaresult.xlsx',sheetName = "Sheet1", append = F)







fplot1 <- plot_rem(mv2@metaresult, "MV", ".", "Gene.symbol", 0.01)
plot(fplot1)











# forestplot for genes which combined p are sig

library(sqldf)

library(readxl)
library(readr)

adjp_sig_ES <- read_excel("C:\\Users\\zhu66\\OneDrive - University of Oklahoma\\Desktop\\GEO\\meta-analysis\\meta1\\using metavolcanoR package\\4adjp and ES.xlsx",sheet = "Sheet1")

symbol_matix <- as.list(adjp_sig_ES$symbol)
library(forestplot)
help("forestplot")

adjp_sig_ES %>%
  forestplot(mean = combine_ES,
             lower = CI_lower,
             upper = CI_upper,
             labeltext = symbol,
             graph.pos = 2,
            # clip=c(-1,2),
            # title = "BaP study",
             zero = c(0.98, 1.02),
            txt_gp = fpTxtGp(ticks=gpar(cex=0.8),label=gpar(cex=0.5) ),
             
             col = fpColors(box = "royalblue",
                            line = "darkblue",summary = "royalblue"),
             
             boxsize = 0.25,
             
           
            # xlab = "Effect Size",
             new_page = TRUE)
             #label = list(gpar(fontfamily = font),
            #              gpar(fontfamily = "",
             #                  col = "#660000")),
             
             
             




# using example data
mv <- rem_mv(diffexplist, metathr = 0.1)
gg <- draw_forest(mv, gene="DPPA2")
plot(gg)

gg <- plot_rem(mv@metaresult, "MV", ".", "Symbol", 0.01)
plot(gg)





