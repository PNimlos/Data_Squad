# install.packages("") #install any packages not on local device

library(stats)
library(ade4)
library(ape)
library(evemodel)
library(ggpubr)

##### tree #####
speciesTree <- read.tree("pero_tree.txt") 
plot(speciesTree) #load in and verify tree branches



##### expression matrices for all samples and Mt. Evans strain only (high altitude) #####
allStrains_exprMat <- read.csv("Hisat_norm_counts_111422.csv")
rownames(allStrains_exprMat) <- allStrains_exprMat$X
allStrains_exprMat$X <- NULL 
ME_exprMat <- read.csv("ME_Hisat_normcounts_1000_111422.csv")
rownames(ME_exprMat) <- ME_exprMat$X
ME_exprMat$X <- NULL #remove X column, leaving only individuals
ME_exprMat <- na.omit(ME_exprMat)
#ME_exprMat <- as.matrix(ME_exprMat) troubleshooting convert to matrix
#ME_exprMat <- as.data.frame(ME_exprMat) troubleshooting convert to matrix



##### id and experimental information ##### 
colnames(allStrains_exprMat)
id <- colnames(allStrains_exprMat)
indiv <- read.csv("acclim_data_all.csv") #dataset contains all individual names and experimental information
indiv <- indiv[order(indiv$Mouse_ID),]
all <- indiv[indiv$Mouse_ID %in% id,] #experimental data for just sequenced individuals 



#### TwoTheta Model ####
#### ME as foreground ####
subtree.ME <- drop.tip(speciesTree, c("MUS", "AZ","ML"))
plot(subtree.ME)
getEdgesFromMRCA(subtree.ME, subtree.ME$tip.label, includeEdgeToMRCA = TRUE) #check number of edges from altered tree
thetaedges <- c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE) #logically specify which branch is foreground vs. background, ME in this case

colSpecies <- colnames(ME_exprMat) #object for EVE input
colSpecies <- sub("_.*$","",colnames(ME_exprMat)) #remove all characters of sample IDs after "_" to match to tree tip labels
colSpecies
#ME_exprMat <- ME_exprMat[-c(2918), ] to remove problematic gene rows
res_ME_1000 <- twoThetaTest(tree = subtree.ME, gene.data = ME_exprMat, isTheta2edge = thetaedges, colSpecies = colSpecies) #run model with all relevant inputs



#### ML as foreground ####
subtree.ML <- drop.tip(speciesTree, c("MUS", "AZ","ME"))
plot(subtree.ML)
getEdgesFromMRCA(subtree.ML, subtree.ML$tip.label, includeEdgeToMRCA = TRUE) #check number of edges from altered tree
thetaedges <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)#logically specify which branch is foreground vs. background, ML in this case

ML_exprMat <- read.csv("ML_Hisat_normcounts_1000_111022.csv")
rownames(ML_exprMat) <- ML_exprMat$X
ML_exprMat$X <- NULL 
ML_exprMat$X.1 <- NULL #remove X column, leaving only individuals
ML_exprMat <- na.omit(ML_exprMat)

colSpecies <- colnames(ML_exprMat) #change depending on treatment
colSpecies <- sub("_.*$","",colnames(ML_exprMat)) #remove all characters of sample IDs after "_" to match to tree tip labels
colSpecies

res_ML_4500 <- twoThetaTest(tree = subtree.ML, gene.data = ML_exprMat, isTheta2edge = thetaedges, colSpecies = colSpecies)



#### AZ as foreground ####
subtree.AZ <- drop.tip(speciesTree, c("MUS", "ML","ME"))
plot(subtree.AZ)
getEdgesFromMRCA(subtree.AZ, subtree.AZ$tip.label, includeEdgeToMRCA = TRUE) #check number of edges from altered tree
thetaedges <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)#logically specify which branch is foreground vs. background, AZ in this case

AZ_exprMat <- read.csv("AZ_Hisat_normcounts_1000_111022.csv")
rownames(AZ_exprMat) <- AZ_exprMat$X
AZ_exprMat$X <- NULL 
AZ_exprMat$X.1 <- NULL #remove X column, leaving only individuals
AZ_exprMat <- na.omit(AZ_exprMat)

colSpecies <- colnames(AZ_exprMat) #change depending on treatment
colSpecies <- sub("_.*$","",colnames(AZ_exprMat)) #remove all characters of sample IDs after "_" to match to tree tip labels
colSpecies

res_AZ_1000 <- twoThetaTest(tree = subtree.AZ, gene.data = AZ_exprMat, isTheta2edge = thetaedges, colSpecies = colSpecies)

#All high altitude as foreground
# subtree.high <- drop.tip(speciesTree, c("MUS"))
# plot(subtree.high)
# getEdgesFromMRCA(subtree.high, subtree.high$tip.label, includeEdgeToMRCA = TRUE) #check number of edges from altered tree
# thetaedges <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE) #logically specify which branch is foreground vs. background, ML in this case
# high_exprMat <- read.csv("Hisat_normcounts_1000.csv")
# rownames(high_exprMat) <- high_exprMat$X
# high_exprMat$X <- NULL #remove X column, leaving only individuals
# high_exprMat <- as.matrix(high_exprMat)
# high_exprMat <- na.omit(high_exprMat)
# 
# colSpecies <- colnames(high_exprMat) #change depending on treatment
# colSpecies <- sub("_.*$","",colnames(high_exprMat)) #remove all characters of sample IDs after "_" to match to tree tip labels
# colSpecies
# 
# res_high_1000 <- twoThetaTest(tree = subtree.high, gene.data = high_exprMat, isTheta2edge = thetaedges, colSpecies = colSpecies)


#### calculating pvals from LRTs and writing results to sheet ####
results <- res_ME_4500$twoThetaRes$par
pval = pchisq(res_ME_4500$LRT, df=1, lower.tail=F)
results <- cbind(results, res_ME_4500$LRT)
Eve_results <- cbind(rownames(ME_exprMat),results, pval, -log10(pval))
#write.csv(Eve_results, file = "Hisat_twoTheta_ME_4500.csv") 

# Expression_shifts_1000 <- read.csv("Hisat_twoTheta_high_1000.csv") #change depending on treatment
# scatter_plot3 <- ggscatter(Expression_shifts_1000, x= 'logtheta', y='X.1', size = (2.0),
# )
# scatter_plot3 + geom_hline(yintercept = 1.3, linetype = "dashed", color = "red")
