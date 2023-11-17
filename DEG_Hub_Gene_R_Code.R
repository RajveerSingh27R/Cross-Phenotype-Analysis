###DEG_Analysis_code/GEO2R#####

getwd()
setwd("~/GEO2R")

#   Differential expression analysis with limma
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE15932", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "XXXXXXXX00000000XXXXXXXX11111111"
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
groups <- make.names(c("T2DM","Healthy"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SPOT_ID","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
                   "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
                   par(mar=c(7,4,2,1))
                   title <- paste ("GSE15932", "/", annotation(gset), sep ="")
                   boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
                   legend("topleft", groups, fill=palette(), bty="n")
                   

                   ################## R code for filteration ##################
##########GSE15932#######
getwd()
setwd("~/GEO2R")
                   
rm(list = ls())
PC<- read.delim("~/GEO2R/PC_H_GSE15932.top.table.tsv")
T2D<- read.delim("~/GEO2R/T2D_H_GSE15932.top.table.tsv")

library(dplyr)
## X = T2D
## Y = PC
T2D_PC=merge(T2D, PC, by.x = "ID", by.y = "ID")

#Calculating Distance
# d = |x_i  - y_i| divided by max(|x_i|, |y_i|)
# filter the genes with d<=10%
# where:
# x_i: expression values for ith gene in T2D
# y_i: expression values for ith gene in PC

Dist=((abs((T2D$logFC)-(PC$logFC)))/(max((abs(T2D$logFC)),(abs(PC$logFC)))))
T2D_PC_Dist=data.frame(T2D_PC,"Dist"=Dist)

Z=unique(T2D_PC_Dist$Gene.symbol.y)

write.csv(T2D_PC_Dist, "GEO2R_T2D_PC_GSE1932.csv", quote = FALSE, row.names = FALSE)



############
filter_pValue_T2D=T2D_PC_Dist%>%
  filter(adj.P.Val.x<=0.05)

filter_pValue_PC=filter_pValue_T2D%>%
  filter(adj.P.Val.y<=0.05)

A=unique(filter_pValue_PC$Gene.symbol.x)

write.csv(filter_pValue_PC, "GEO2R_T2D_PC_GSE1932_filter1.csv", quote = FALSE, row.names = FALSE)
write.csv(A, "GEO2R_T2D_PC_GSE1932_filter1_unique.csv", quote = FALSE, row.names = FALSE)

filter_T2D_FC_A=filter_pValue_PC%>%
  filter(logFC.x>=0.5)

filter_T2D_FC_B=filter_pValue_PC%>%
  filter(logFC.x<=-0.5)

filter_T2D=rbind(filter_T2D_FC_A, filter_T2D_FC_B)

AA=unique(filter_T2D$Gene.symbol.x)

filter_PC_FC_A=filter_T2D%>%
  filter(logFC.y>=0.5)

filter_PC_FC_B=filter_T2D%>%
  filter(logFC.y<=-0.5)

filter_PC=rbind(filter_PC_FC_A, filter_PC_FC_B)

AB=unique(filter_PC$Gene.symbol.y)

write.csv(filter_PC, "GEO2R_T2D_PC_GSE1932_filter2.csv", quote = FALSE, row.names = FALSE)
write.csv(AB, "GEO2R_T2D_PC_GSE1932_filter2_unique.csv", quote = FALSE, row.names = FALSE)

filter_dist=filter_PC%>%
  filter(Dist<=0.1)
AAx=unique(filter_dist$Gene.symbol.x)
AAy=unique(filter_dist$Gene.symbol.y)

write.csv(filter_dist, "GEO2R_T2D_PC_GSE1932_filter3.csv", quote = FALSE, row.names = FALSE)
write.csv(AAx, "GEO2R_T2D_PC_GSE1932_filter3_unique.csv", quote = FALSE, row.names = FALSE)

###############GSE15932 _T2D_PC both###########
getwd()
setwd("~/GEO2R")

rm(list = ls())
T2D_PC<- read.delim("D:/PhD/T2D_Cancer_Project/Geo2R/Geo/T2D_PC_H_GSE15932.top.table.tsv")
Z=unique(T2D_PC$Gene.symbol)


#########
filter_pValue_T2D=T2D_PC%>%
  filter(adj.P.Val<=0.05)
A=unique(filter_pValue_T2D$Gene.symbol)

write.csv(filter_pValue_T2D, "GEO2R_T2D&PC_0.05_GSE1932_filter1.csv", quote = FALSE, row.names = FALSE)
write.csv(A, "GEO2R_T2D&PC_0.05_GSE1932_filter1_unique.csv", quote = FALSE, row.names = FALSE)

filter_T2D_FC_A=filter_pValue_T2D%>%
  filter(logFC>=0.5)

filter_T2D_FC_B=filter_pValue_T2D%>%
  filter(logFC<=-0.5)

filter_T2D=rbind(filter_T2D_FC_A, filter_T2D_FC_B)
AA=unique(filter_T2D$Gene.symbol)

write.csv(filter_T2D, "GEO2R_T2D&PC_0.05_GSE1932_filter2.csv", quote = FALSE, row.names = FALSE)
write.csv(AA, "GEO2R_T2D&PC_0.05_GSE1932_filter2_unique.csv", quote = FALSE, row.names = FALSE)
