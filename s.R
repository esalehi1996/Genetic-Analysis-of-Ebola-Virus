
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#require(Biostrings)
library("seqinr")
require(Biostrings)
library("phangorn")
setwd("D:/Course Stuff/96-97-Fall/Bioinformatics/Project")
Bundibugyo_genome <- read.fasta("Bundibugyo_genome.fasta")
Reston_genome <- read.fasta("Reston_genome.fasta")
Sudan_genome <- read.fasta("Sudan_genome.fasta")
TaiForest_genome <- read.fasta("TaiForest_genome.fasta")
Zaire_genome <- read.fasta("Zaire_genome.fasta")
Marburg_genes <- read.fasta("Marburg_Genes.fasta")
data("BLOSUM50")
scoremat <- BLOSUM50
B_genome <- toupper(paste(Bundibugyo_genome[[1]][1:length(Bundibugyo_genome[[1]])],collapse = ""))
R_genome <- toupper(paste(Reston_genome[[1]][1:length(Reston_genome[[1]])],collapse = ""))
S_genome <- toupper(paste(Sudan_genome[[1]][1:length(Sudan_genome[[1]])],collapse = ""))
T_genome <- toupper(paste(TaiForest_genome[[1]][1:length(TaiForest_genome[[1]])],collapse = ""))
Z_genome <- toupper(paste(Zaire_genome[[1]][1:length(Zaire_genome[[1]])],collapse = ""))
Marburg_NP_genes <- toupper(paste(Marburg_genes[[1]][1:length(Marburg_genes[[1]])],collapse = ""))
Marburg_VP35_genes <- toupper(paste(Marburg_genes[[2]][1:length(Marburg_genes[[2]])],collapse = ""))
Marburg_VP40_genes <- toupper(paste(Marburg_genes[[3]][1:length(Marburg_genes[[3]])],collapse = ""))
Marburg_GP_genes <- toupper(paste(Marburg_genes[[4]][1:length(Marburg_genes[[4]])],collapse = ""))
Marburg_VP30_genes <- toupper(paste(Marburg_genes[[5]][1:length(Marburg_genes[[5]])],collapse = ""))
Marburg_VP24_genes <- toupper(paste(Marburg_genes[[6]][1:length(Marburg_genes[[6]])],collapse = ""))
Marburg_L_genes <- toupper(paste(Marburg_genes[[7]][1:length(Marburg_genes[[7]])],collapse = ""))
# Finding the NP gene in each genome
NP_gene <- vector(mode = "list" , length = 5)
NP_gene[1] <- attributes(pairwiseAlignment(B_genome, Marburg_NP_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
NP_gene[2] <- attributes(pairwiseAlignment(R_genome, Marburg_NP_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
NP_gene[3] <- attributes(pairwiseAlignment(S_genome, Marburg_NP_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
NP_gene[4] <- attributes(pairwiseAlignment(T_genome, Marburg_NP_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
NP_gene[5] <- attributes(pairwiseAlignment(Z_genome, Marburg_NP_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
NP_edit_matrix = matrix(nrow=5,ncol=5)
for (i in 1:5){
  for (j in 1:5){
    NP_edit_matrix[i,j] = adist(gsub("-","",toupper(paste(NP_gene[[i]],collapse=""))), gsub("-","",toupper(paste(NP_gene[[j]],collapse=""))))
  }}
rownames(NP_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
colnames(NP_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
tree1 <- upgma(NP_edit_matrix)
jpeg('Phylogeny Trees/UPGMA_NP.jpg')
plot(tree1)
dev.off()
tree8 <- nj(NP_edit_matrix)
jpeg('Phylogeny Trees/NJ_NP.jpg')
plot(tree8)
dev.off()
# Finding the VP35 gene in each genome
VP35_gene <- vector(mode = "list" , length = 5)
VP35_gene[1] <- attributes(pairwiseAlignment(B_genome, Marburg_VP35_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP35_gene[2] <- attributes(pairwiseAlignment(R_genome, Marburg_VP35_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP35_gene[3] <- attributes(pairwiseAlignment(S_genome, Marburg_VP35_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP35_gene[4] <- attributes(pairwiseAlignment(T_genome, Marburg_VP35_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP35_gene[5] <- attributes(pairwiseAlignment(Z_genome, Marburg_VP35_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP35_edit_matrix = matrix(nrow=5,ncol=5)
for (i in 1:5){
  for (j in 1:5){
   VP35_edit_matrix[i,j] = adist(gsub("-","",toupper(paste(VP35_gene[[i]],collapse=""))), gsub("-","",toupper(paste(VP35_gene[[j]],collapse=""))))
  }}
rownames(VP35_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
colnames(VP35_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
tree2 <- upgma(VP35_edit_matrix)
jpeg('Phylogeny Trees/UPGMA_VP35.jpg')
plot(tree2)
dev.off()
tree9 <- nj(VP35_edit_matrix)
jpeg('Phylogeny Trees/NJ_VP35.jpg')
plot(tree9)
dev.off()
# Finding the VP40 gene
VP40_gene <- vector(mode = "list" , length = 5)
VP40_gene[1] <- attributes(pairwiseAlignment(B_genome, Marburg_VP40_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP40_gene[2] <- attributes(pairwiseAlignment(R_genome, Marburg_VP40_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP40_gene[3] <- attributes(pairwiseAlignment(S_genome, Marburg_VP40_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP40_gene[4] <- attributes(pairwiseAlignment(T_genome, Marburg_VP40_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP40_gene[5] <- attributes(pairwiseAlignment(Z_genome, Marburg_VP40_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP40_edit_matrix = matrix(nrow=5,ncol=5)
for (i in 1:5){
  for (j in 1:5){
    VP40_edit_matrix[i,j] = adist(gsub("-","",toupper(paste(VP40_gene[[i]],collapse=""))), gsub("-","",toupper(paste(VP40_gene[[j]],collapse=""))))
  }}
rownames(VP40_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
colnames(VP40_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
tree3 <- upgma(VP40_edit_matrix)
jpeg('Phylogeny Trees/UPGMA_VP40.jpg')
plot(tree3)
dev.off()
tree10 <- nj(VP40_edit_matrix)
jpeg('Phylogeny Trees/NJ_VP40.jpg')
plot(tree10)
dev.off()
# Finding the GP gene
GP_gene <- vector(mode = "list" , length = 5)
GP_gene[1] <- attributes(pairwiseAlignment(B_genome, Marburg_GP_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
GP_gene[2] <- attributes(pairwiseAlignment(R_genome, Marburg_GP_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
GP_gene[3] <- attributes(pairwiseAlignment(S_genome, Marburg_GP_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
GP_gene[4] <- attributes(pairwiseAlignment(T_genome, Marburg_GP_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
GP_gene[5] <- attributes(pairwiseAlignment(Z_genome, Marburg_GP_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
GP_edit_matrix = matrix(nrow=5,ncol=5)
for (i in 1:5){
  for (j in 1:5){
    GP_edit_matrix[i,j] = adist(gsub("-","",toupper(paste(GP_gene[[i]],collapse=""))), gsub("-","",toupper(paste(GP_gene[[j]],collapse=""))))
  }}
rownames(GP_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
colnames(GP_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
tree4 <- upgma(GP_edit_matrix)
jpeg('Phylogeny Trees/UPGMA_GP.jpg')
plot(tree4)
dev.off()
tree11 <- nj(GP_edit_matrix)
jpeg('Phylogeny Trees/NJ_GP.jpg')
plot(tree11)
dev.off()
# Finding the VP30 gene
VP30_gene <- vector(mode = "list" , length = 5)
VP30_gene[1] <- attributes(pairwiseAlignment(B_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP30_gene[2] <- attributes(pairwiseAlignment(R_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP30_gene[3] <- attributes(pairwiseAlignment(S_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP30_gene[4] <- attributes(pairwiseAlignment(T_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP30_gene[5] <- attributes(pairwiseAlignment(Z_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP30_edit_matrix = matrix(nrow=5,ncol=5)
for (i in 1:5){
  for (j in 1:5){
    VP30_edit_matrix[i,j] = adist(gsub("-","",toupper(paste(VP30_gene[[i]],collapse=""))), gsub("-","",toupper(paste(VP30_gene[[j]],collapse=""))))
  }}
rownames(VP30_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
colnames(VP30_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
tree5 <- upgma(VP30_edit_matrix)
jpeg('Phylogeny Trees/UPGMA_VP30.jpg')
plot(tree5)
dev.off()
tree12 <- nj(VP30_edit_matrix)
jpeg('Phylogeny Trees/NJ_VP30.jpg')
plot(tree12)
dev.off()
# Finding the VP24 gene
VP24_gene <- vector(mode = "list" , length = 5)
VP24_gene[1] <- attributes(pairwiseAlignment(B_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP24_gene[2] <- attributes(pairwiseAlignment(R_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP24_gene[3] <- attributes(pairwiseAlignment(S_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP24_gene[4] <- attributes(pairwiseAlignment(T_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP24_gene[5] <- attributes(pairwiseAlignment(Z_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
VP24_edit_matrix = matrix(nrow=5,ncol=5)
for (i in 1:5){
  for (j in 1:5){
    VP24_edit_matrix[i,j] = adist(gsub("-","",toupper(paste(VP24_gene[[i]],collapse=""))), gsub("-","",toupper(paste(VP24_gene[[j]],collapse=""))))
  }}
rownames(VP24_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
colnames(VP24_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
tree6 <- upgma(VP24_edit_matrix)
jpeg('Phylogeny Trees/UPGMA_VP24.jpg')
plot(tree6)
dev.off()
tree13 <- nj(VP24_edit_matrix)
jpeg('Phylogeny Trees/NJ_VP24.jpg')
plot(tree13)
dev.off()
# Finding the L gene
L_gene <- vector(mode = "list" , length = 5)
L_gene[1] <- attributes(pairwiseAlignment(B_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
L_gene[2] <- attributes(pairwiseAlignment(R_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
L_gene[3] <- attributes(pairwiseAlignment(S_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
L_gene[4] <- attributes(pairwiseAlignment(T_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
L_gene[5] <- attributes(pairwiseAlignment(Z_genome, Marburg_VP30_genes, substitutionMatrix = scoremat , type = "local-global"))$pattern
L_edit_matrix = matrix(nrow=5,ncol=5)
for (i in 1:5){
  for (j in 1:5){
    L_edit_matrix[i,j] = adist(gsub("-","",toupper(paste(L_gene[[i]],collapse=""))), gsub("-","",toupper(paste(L_gene[[j]],collapse=""))))
  }}
rownames(L_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
colnames(L_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
tree7 <- upgma(L_edit_matrix)
jpeg('Phylogeny Trees/UPGMA_L.jpg')
plot(tree7)
dev.off()
tree14 <- nj(L_edit_matrix)
jpeg('Phylogeny Trees/NJ_L.jpg')
plot(tree14)
dev.off()
##
Genome_edit_matrix = matrix(nrow=5,ncol=5)
Genome_list <- list()
Genome_list[[1]] <- B_genome
Genome_list[[2]] <- R_genome
Genome_list[[3]] <- S_genome
Genome_list[[4]] <- T_genome
Genome_list[[5]] <- Z_genome
for (i in 1:5){
  for (j in 1:5){
    Genome_edit_matrix[i,j] <- adist(Genome_list[[i]],Genome_list[[j]])
  }
}
rownames(Genome_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
colnames(Genome_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire")
upgma_consensus_tree <- consensus(tree1,tree2, tree3,tree4,tree5,tree6,tree7, p = 0.3, check.labels = T)
nj_consensus_tree <- consensus(tree8,tree9, tree10,tree11,tree12,tree13,tree14, p = 0.3, check.labels = T)
upgma_genome_tree <- upgma(Genome_edit_matrix)
nj_genome_tree <- nj(Genome_edit_matrix)

##

Marburg_genome <- read.fasta("Marburg_genome.fasta")
M_genome <- toupper(paste(Marburg_genome[[1]][1:length(Marburg_genome[[1]])],collapse = ""))
Genome_list[[6]] <- M_genome
Genome_edit_matrix = matrix(nrow=6,ncol=6)
for (i in 1:6){
  for (j in 1:6){
    Genome_edit_matrix[i,j] <- pairwiseAlignment(Genome_list[[i]],Genome_list[[j]],substitutionMatrix = scoremat , type = "global" , scoreOnly = TRUE)
  }
}
rownames(Genome_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire","Marburg")
colnames(Genome_edit_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire","Marburg")
nj_whole_genome_tree <- nj(Genome_edit_matrix)

##

edit_dist_matrix <- matrix(nrow=6,ncol=6)
for (i in 1:6){
  for (j in 1:6){
    edit_dist_matrix[i,j] <- adist(Genome_list[[i]],Genome_list[[j]])
  }
}
rownames(edit_dist_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire","Marburg")
colnames(edit_dist_matrix) <- c("Bundibugyo","Reston","Sudan","TaiForest","Zaire","Marburg")
nj_whole_genome_tree <- nj(edit_dist_matrix)
max_genome_length <- max(c(nchar(M_genome),nchar(B_genome),nchar(R_genome),nchar(S_genome),nchar(T_genome),nchar(Z_genome)))
P_matrix <- edit_dist_matrix / max_genome_length
T_matrix <- (-3*log(1 - 4*P_matrix/3)/4)/(1.9*10^-3)
T_tree <- nj(T_matrix)
