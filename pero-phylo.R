# install.packages("adegenet", dep = TRUE) #install once
# install.packages("phangorn", dep = TRUE) #install once
library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)
library(evemodel)
library(edgeR)
library(ggpubr)

pero <- read.dna("pero_alignment_2.txt", format="fasta") #created alignment in clustalW web (https://www.ebi.ac.uk/Tools/msa/clustalo/) #output as pearson/fasta
D <- dist.dna(pero, model = "TN93")
length(D)
# code to write a UPGMA tree and check that it's a good fit to the data
x <- as.vector(D)
tree <- as.phylo(hclust(D,method="average")) #code for a UPGMA tree
y2 <- as.vector(as.dist(cophenetic(tree)))
plot(x, y2, xlab="original pairwise distances", y2lab="pairwise distances on the tree", main="Is UPGMA appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y2~x), col="red")
plot(tree, cex=1) #plot the tree
ape::write.tree(tree, file='pero_tree_noMUSorALS.txt') #export tree in newick
#write.nexus(tree, file='PGLS_tree.nex')
speciesTree <- read.tree("pero_tree.txt") #load tree in newick to verify for EVE
# plot the species tree

####### id and experimental information ####### 
id <- read.delim("combined_names.txt", header=FALSE)
id <- as.vector(id$V1)
id <- gsub("\\_.*","",id)#vector of ID names in order
id <- gsub("-", "_", id)
indiv <- read.csv("indiv_info_all.csv") #dataset contains all individual names and experimental information
indiv <- indiv[order(indiv$Mouse_ID),]
all <- indiv[indiv$Mouse_ID %in% id,] #experimental data for all sequenced individuals #2 MUS individuals dropped bc of ambiguous names
#######


stampy.counts <- read.table("stampy_counts.txt", header=TRUE, row.names=1) #counts table for reads mapped with Stampy
stampy.counts <- stampy.counts[, -c(1:5)]
colnames(stampy.counts) <- id #rename the counts table with just id's

# stampy.matrix <- as.matrix(stampy.counts)
# normal.counts <- DGEList(counts = stampy.matrix)
# # read in raw counts and create list object for normalization
# keep <- filterByExpr(normal.counts)
# normal.counts <- normal.counts[keep, keep.lib.sizes=FALSE]
# normal.counts <- calcNormFactors(normal.counts)
# #filter low expression counts out of matrix and normalize
# keep_counts <- as.matrix(normal.counts)

########
#code for design matrix and low expression filtering

match(all$Mouse_ID, colnames(exprMat))
treatment <- all$Treatment #this is my first factor for the GLM model
species <- all$Strain #this is my second factor
table <- data.frame(treatment, species)
group <- factor(paste(table$treatment, table$species, sep="_"))
cbind(table, group=group)
# table$treatment = as.factor(table$treatment)
# table$treatment <- relevel(table$treatment, ref="Control") #set control as the reference level for treatment
design <- model.matrix(~treatment*species, data=table)

#filter data
exprMat$mean = rowMeans(exprMat) #rowMeans takes means of each row
keep_counts = subset(counts, mean >= 10) #filter out genes that have <10 reads on average across individuals
dim(keep_counts)
keep_counts$mean = NULL #clean up dataset
dim(keep_counts)
#keep_counts <- keep_counts[,-(100:101),drop=FALSE] #drop 2 MUS samples that do not match Mouse ID grouping

#######
#code for normalization

y <- DGEList(counts=keep_counts, group=group) # make a DGE list
dim(y)
y <- calcNormFactors(y) # normalize
plotMDS(y) #exploration plots
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotMDS(y)
counts.norm <- cpm(y, log=TRUE, prior.count=2, normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset
#write.csv(counts.norm, file="norm_counts.csv")

#####  PCA w Mus #######
norm.counts <- read.csv("norm_counts.csv")
norm.counts1 <- norm.counts
rownames(norm.counts1) <- norm.counts1$Gene
norm.counts1$Gene <- NULL
head(norm.counts1)
#norm.counts2 <- norm.counts1[, -c(87:102)] #remove MUS
pca <- prcomp(t(norm.counts1), scale=FALSE)
summary(pca)
pc = as.data.frame(pca$x)
pc$treatment <- all$Treatment
pc$species <- all$Strain
az <- subset(pc, species == "AZ")
# az1 <- subset(az, treatment=="1000")
# az2 <- subset(az, treatment=="3000")
# az3 <- subset(az, treatment=="4500")
ml <- subset(pc, species == "ML")
# ml1 <- subset(ml, treatment=="1000")
# ml2 <- subset(ml, treatment=="3000")
# ml3 <- subset(ml, treatment=="4500")
ln <- subset(pc, species == "LN")
me <- subset(pc, species == "ME")
po <- subset(pc, species == "PO")
ls <- subset(pc, species == "LS")
mus <- subset(pc, species=="MUS")
#pdf("pca.pdf", h=6, w=6)
plot(az$PC1,az$PC2, pch=21, xlim=c(-100,40), ylim = c(-40,60), cex=2, xlab="PC1 (31%)", ylab="PC2 (18%) ", col="black", bg="#9d9d3b", lwd=2, las=1)
points(ml$PC1, ml$PC2, pch=21, col="black", bg="#b45ac2",cex=2, lwd=2)
points(ln$PC1, ln$PC2, pch=21, col="black", bg="#58a968",cex=2, lwd=2)
points(me$PC1, me$PC2, pch=21, col="black", bg="#ce4948",cex=2, lwd=2)
points(po$PC1, po$PC2, pch=21, col="black", bg="#6e80cc",cex=2, lwd=2)
points(ls$PC1, ls$PC2, pch=21, col="black", bg="#c87e42",cex=2, lwd=2)
points(mus$PC1, mus$PC2, pch=21, col="black", bg="#c4608c",cex=2, lwd=2)
#dev.off()



###########################
#subset normcounts for each treatment
norm.counts <- read.csv("norm_counts.csv")
all1000 <- subset(all, treatment == 1000)
all3000 <- subset(all, treatment == 3000)
all4500 <- subset(all, treatment == 4500)
match(all1000$Mouse_ID, colnames(norm.counts)) #find which individuals match this treatment
norm.counts1000 <- norm.counts[,c(1,2,7,8,9,15,19,23,27,30,32,34,36,39,40,43,44,46,50,57,60,65,67,71,74,76,77,78,85,88,89,92,96,97,98,104,106,108,112,114,120)] #new dataframe with only individuals from that treatment THERE HAS TO BE A MORE ELEGANT WAY OF DOING THIS
norm.counts3000 <- norm.counts[,c(1,4,5,10,11,16,18,20,24,25,28,33,35,38,47,49,51,52,55,56,58,64,69,70,72,73,75,80,82,91,93,94,99,109,110,111,115,116,121)]
norm.counts4500 <- norm.counts[,c(1,3,6,12,13,14,17,21,22,26,29,31,37,41,42,45,48,53,59,61,62,63,66,68,79,81,83,84,86,87,90,95,100,101,102,103,105,107,113,117,118,119)]
#write.csv(norm.counts1000, file="norm.counts1000.csv")
#write.csv(norm.counts3000, file="norm.counts3000.csv")
#write.csv(norm.counts4500, file="norm.counts4500.csv")
#need to delete first column in each csv (labeled each gene numerically for unknown reason)

### JPV code to subset by treatment ###
norm.counts3000 <- norm.counts
rownames(norm.counts3000) <- norm.counts3000$X
norm.counts3000$X <- NULL
head(norm.counts3000)
#
id3000 <- all3000$Mouse_ID
counts3000 <- norm.counts3000[,colnames(norm.counts3000) %in% id3000]
match(colnames(counts3000), id3000)
counts3000$X <- norm.counts$X
write.csv(counts3000, file = "Hisat_normcounts_allStrains_3000.csv")



###########################
#TJ EVE code SharedBeta Model
#Input data
#exprMat <- as.matrix(read.csv("Hisat_normcounts_1000.csv"))
exprMat <- read.csv("Hisat_normcounts_1000.csv")
rownames(exprMat) <- exprMat$X
exprMat$X <- NULL #remove X column, leaving only individuals
exprMat <- as.matrix(exprMat)
#change input depending on treatment
speciesTree <- read.tree("pero_tree.txt")
#Before running EVE make sure that your samples have the species name with and underscore then the sample number and that the species name matches what is in the tree
#Example for Americans the sample names are 'am_1', 'am_2' and so on
plot(speciesTree)
speciesTree$tip.label
colSpecies <- colnames(exprMat) #change depending on treatment
colSpecies <- sub("_.*$","",colnames(exprMat)) #remove all characters of sample IDs after "_" to match to tree tip labels
colSpecies
#run EVE model
res <- betaSharedTest(tree = speciesTree, gene.data = exprMat, colSpecies = colSpecies)
results <- res$indivBetaRes$par
pval = pchisq(res$LRT, df=1, lower.tail=F)
results <- cbind(results, res$LRT)
Eve_results <- cbind(rownames(exprMat),results, pval, -log10(pval))
write.csv(Eve_results, file = "Hisat_BetaShared_4500.csv") 
#change depending on treatment
#Note that to convert Beta values to Expression divergence -log10 normalize the beta values, that way the expression divergence is signed and normal. Also, convert pvalues to -log10. 

Expression_shifts_4500 <- read.csv("Hisat_BetaShared_4500.csv") #change depending on treatment
scatter_plot3 <- ggscatter(Expression_shifts_4500, x= 'Expression_divergence', y='log_pval', size = (2.0),
)
scatter_plot3 + geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")


#TwoTheta Model
#ME as foreground
speciesTree <- read.tree("pero_tree.txt")
subtree.ME <- drop.tip(speciesTree, c("MUS", "AZ","ML"))
plot(subtree.ME)
getEdgesFromMRCA(subtree.ME, subtree.ME$tip.label, includeEdgeToMRCA = TRUE) #check number of edges from altered tree
thetaedges <- c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE) #logically specify which branch is foreground vs. background, ME in this case
ME_exprMat <- read.csv("ME_Hisat_normcounts_1000_111022.csv")
rownames(ME_exprMat) <- ME_exprMat$X
ME_exprMat$X <- NULL 
ME_exprMat$X.1 <- NULL#remove X column, leaving only individuals
ME_exprMat <- as.matrix(ME_exprMat)
ME_exprMat <- na.omit(ME_exprMat)
###ME_exprMat <- as.matrix(subset(ME_exprMat, select =-c(ML*,))) necessary to drop ML from gene data too?

colSpecies <- colnames(ME_exprMat) #change depending on treatment
colSpecies <- sub("_.*$","",colnames(ME_exprMat)) #remove all characters of sample IDs after "_" to match to tree tip labels
colSpecies
#ME_exprMat <- ME_exprMat[-c(2918), ] to remove problematic gene rows
res_ME_1000 <- twoThetaTest(tree = subtree.ME, gene.data = ME_exprMat, isTheta2edge = thetaedges, colSpecies = colSpecies)

#ML as foreground
subtree.ML <- drop.tip(speciesTree, c("MUS", "AZ","ME"))
plot(subtree.ML)
getEdgesFromMRCA(subtree.ML, subtree.ML$tip.label, includeEdgeToMRCA = TRUE) #check number of edges from altered tree
thetaedges <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) #logically specify which branch is foreground vs. background, ML in this case
ML_exprMat <- read.csv("ML_Hisat_normcounts_4500.csv")
rownames(ML_exprMat) <- ML_exprMat$X
ML_exprMat$X <- NULL #remove X column, leaving only individuals
ML_exprMat <- as.matrix(ML_exprMat)
ML_exprMat <- na.omit(ML_exprMat)
###ME_exprMat <- as.matrix(subset(ME_exprMat, select =-c(ML*,))) necessary to drop ML from gene data too?

colSpecies <- colnames(ML_exprMat) #change depending on treatment
colSpecies <- sub("_.*$","",colnames(ML_exprMat)) #remove all characters of sample IDs after "_" to match to tree tip labels
colSpecies

res_ML_4500 <- twoThetaTest(tree = subtree.ML, gene.data = ML_exprMat, isTheta2edge = thetaedges, colSpecies = colSpecies)
 
#Both ME and ML as foreground
subtree.high <- drop.tip(speciesTree, c("MUS", "AZ"))
plot(subtree.high)
getEdgesFromMRCA(subtree.high, subtree.high$tip.label, includeEdgeToMRCA = TRUE) #check number of edges from altered tree
thetaedges <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE) #logically specify which branch is foreground vs. background, ML in this case
high_exprMat <- read.csv("Hisat_normcounts_4500.csv")
rownames(high_exprMat) <- high_exprMat$X
high_exprMat$X <- NULL #remove X column, leaving only individuals
high_exprMat <- as.matrix(high_exprMat)
high_exprMat <- na.omit(high_exprMat)

colSpecies <- colnames(high_exprMat) #change depending on treatment
colSpecies <- sub("_.*$","",colnames(high_exprMat)) #remove all characters of sample IDs after "_" to match to tree tip labels
colSpecies

res_high_4500 <- twoThetaTest(tree = subtree.high, gene.data = high_exprMat, isTheta2edge = thetaedges, colSpecies = colSpecies)

#Low alt species as foreground
subtree.low <- drop.tip(speciesTree, c("MUS", "AZ"))
plot(subtree.low)
getEdgesFromMRCA(subtree.low, subtree.low$tip.label, includeEdgeToMRCA = TRUE) #check number of edges from altered tree
thetaedges <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE) #logically specify which branch is foreground vs. background, ML in this case
high_exprMat <- read.csv("Hisat_normcounts_4500.csv")
rownames(high_exprMat) <- high_exprMat$X
high_exprMat$X <- NULL #remove X column, leaving only individuals
high_exprMat <- as.matrix(high_exprMat)
high_exprMat <- na.omit(high_exprMat)

colSpecies <- colnames(high_exprMat) #change depending on treatment
colSpecies <- sub("_.*$","",colnames(high_exprMat)) #remove all characters of sample IDs after "_" to match to tree tip labels
colSpecies

res_high_4500 <- twoThetaTest(tree = subtree.high, gene.data = high_exprMat, isTheta2edge = thetaedges, colSpecies = colSpecies)



#calculating pvals from LRTs and writing results to sheet
results <- res_ME_4500$twoThetaRes$par
pval = pchisq(res_ME_4500$LRT, df=1, lower.tail=F)
results <- cbind(results, res_ME_4500$LRT)
Eve_results <- cbind(rownames(ME_exprMat),results, pval, -log10(pval))
write.csv(Eve_results, file = "Hisat_twoTheta_ME_4500.csv") 

Expression_shifts_1000 <- read.csv("Hisat_twoTheta_high_1000.csv") #change depending on treatment
scatter_plot3 <- ggscatter(Expression_shifts_1000, x= 'logtheta', y='X.1', size = (2.0),
)
scatter_plot3 + geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")


############
#PGLS for phenotypic data
############
library(ape)
library(nlme)
library(geiger)
library(caper)
speciesTree <- read.tree("pero_tree.txt")
plot(speciesTree)
# all_data <- read.csv("acclim_data_all.csv")
# drop <- c("Mouse_ID", "Species")
# all_data_names <- all_data[,!(names(all_data) %in% drop)]
PGLSdf <- cbind.data.frame(all_data$Mouse_ID, all_data$Strain, all_data$Treatment, all_data$Elevation, all_data$Mass, all_data$InitialVO2max, all_data$Final_V02max, all_data$hb, all_data$hct)
names(PGLSdf) <- substring(names(PGLSdf), 10) #drop the "all_data$ from all colnames
PGLSdf <- na.omit(PGLSdf)
speciesTree <- drop.tip(speciesTree, c("MUS"))
#name.check(speciesTree, PGLS$Strain)
bm <- corBrownian(1, speciesTree, form = ~Strain)
bm
PGLS1 <-nlme::gls(Final_V02max ~ Strain, data=PGLSdf, correlation=bm)

comp.data <- comparative.data(speciesTree, PGLSdf, names.col = "Strain", warn.dropped = TRUE)
