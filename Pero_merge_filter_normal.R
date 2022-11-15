library(stats)
library(ade4)
library(ape)
library(edgeR)

###### read in unmerged prenormal exprMats ######
AZ_prenorm_exprMat <- read.csv("Prenorm_Hisat_counts_AZ.csv")
AZ_prenorm_exprMat0 <- AZ_prenorm_exprMat
AZ_prenorm_exprMat <- AZ_prenorm_exprMat[order(AZ_prenorm_exprMat$X), ]
rownames(AZ_prenorm_exprMat) <- AZ_prenorm_exprMat$X
AZ_prenorm_exprMat$X <- NULL 

rest_prenorm_exprMat <- read.csv("PrenormalExprMat060522.csv")
rownames(rest_prenorm_exprMat) <- rest_prenorm_exprMat$Geneid
rest_prenorm_exprMat$Geneid <- NULL
rest_id <- colnames(rest_prenorm_exprMat)
rest_id <- gsub("\\.","_",rest_id)
rest_id
colnames(rest_prenorm_exprMat) <- rest_id

####### merge exprMats ######
inner_join <- merge(x = AZ_prenorm_exprMat, y = rest_prenorm_exprMat, by = "row.names", all = FALSE) #inner merge, keeping only gene rows that all species share
rownames(inner_join) <- inner_join$Row.names
inner_join$Row.names <- NULL
#write.csv(inner_join, file = "merged_counts_all_111422.csv")

####### id and experimental information ####### 
id <- colnames(inner_join)
id <- sort(id)
inner_join <- inner_join[, id]

indiv <- read.csv("acclim_data_all.csv") #dataset contains all individual names and experimental information
indiv <- indiv[order(indiv$Mouse_ID),]
all <- indiv[indiv$Mouse_ID %in% id,] #experimental data for all sequenced individuals 

###### design matrix for normalization ######
exprMat <- inner_join
match(all$Mouse_ID, colnames(exprMat)) #verify 1-104, in order
treatment <- all$Treatment #this is my first factor for the GLM model
species <- all$Strain #this is my second factor
table <- data.frame(treatment, species)
group <- factor(paste(table$treatment, table$species, sep="_"))
cbind(table, group=group)
# table$treatment = as.factor(table$treatment)
# table$treatment <- relevel(table$treatment, ref="Control") #set control as the reference level for treatment
design <- model.matrix(~treatment*species, data=table)

#### filtering ####
exprMat$mean = rowMeans(exprMat) #rowMeans takes means of each row
keep_counts = subset(exprMat, mean >= 10) #filter out genes that have <10 reads on average across individuals
dim(keep_counts)
keep_counts$mean = NULL #clean up dataset
dim(keep_counts)
#write.csv(keep_counts, file = "keep_counts_all_111422.csv")

##### normalization with edgeR ####
y <- DGEList(counts=keep_counts, group=group) # make a DGE list
dim(y)
y <- calcNormFactors(y) # normalize
plotMDS(y) #exploration plots
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotMDS(y)
counts.norm <- cpm(y, log=TRUE, prior.count=2, normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset
#write.csv(counts.norm, file="Hisat_norm_counts_111422.csv")
