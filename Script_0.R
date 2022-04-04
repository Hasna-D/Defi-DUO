#Importation des données
# --------------------------------------------------------------
expData <- read.table("cell-cycle_SCERE_DUO.txt", row.names = 1, 
                      sep = "\t", header = T)


# Représentations graphiques des données (boxplot - GeneProfiles2)
# --------------------------------------------------------------

#-- Exemple 1 (Boxplot) :
boxplot(expData)
#-- Exemple 2 (Profils d'expression des gènes - fonction de Gaëlle) :
plotGenes2 <- function(expData, title = "", yMin = 0, yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    yMax = max(expData)
    
  }
  
  # Representation of the first expression profile
  plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
       ylim = c(floor(yMin), ceiling(yMax)),
       xlab = "Time point", ylab = "Gene expression level",
       main = title)
  
  # Add expression profile for other genes
  for(i in 2:nrow(expData)){
    
  lines(1:ncol(expData), expData[i,], col = "grey")
    
    # end of for()  
  }

  
  # Average expression profile
if(meanProfile == TRUE){
expMean = apply(expData, 2, mean)
lines(1:ncol(expData), expMean, col = "red", 
lwd = 1.5, lty = "dashed")
  }
  
  # end of function plotGenes()  
} 

plotGenes2(expData) 

# --------------------------------------------------------------
# Clustering avec la méthode k-means
# --------------------------------------------------------------

#-- Exemple 1 (2 groupes) :
resKmeans <- kmeans(expData, centers = 2)
summary(resKmeans)
resKmeans$cluster
table(resKmeans$cluster)
# Cohérence entre les noms des gènes ? 
row.names(expData)
names(resKmeans$cluster)
row.names(expData) == names(resKmeans$cluster)

# Récupération des profils d'expression des gènes dans le groupe 1
cluster1 <- expData[which(resKmeans$cluster == 1),]
boxplot(cluster1)
plotGenes2(cluster1)
# Récupération des profils d'expression des gènes dans le groupe 2
cluster2 <- expData[which(resKmeans$cluster == 2),]
boxplot(cluster2)
plotGenes2(cluster2)

#-- Exemple 2 (4 groupes)
resKmeans <- kmeans(expData, centers = 4)
table(resKmeans$cluster)

cluster1_Km <- expData[which(resKmeans$cluster == 1),]
cluster2_Km <- expData[which(resKmeans$cluster == 2),]
cluster3_Km <- expData[which(resKmeans$cluster == 3),]
cluster4_Km <- expData[which(resKmeans$cluster == 4),]

plotGenes2(cluster1_Km)
plotGenes2(cluster2_Km)
plotGenes2(cluster3_Km)
plotGenes2(cluster4_Km)

heatmap(as.matrix(cluster1_Km))
heatmap(as.matrix(cluster2_Km))
heatmap(as.matrix(cluster3_Km))
heatmap(as.matrix(cluster4_Km))
# --------------------------------------------------------------
# Changer la distance
# --------------------------------------------------------------

# Création d'une matrice de distance
cor(expData)
cor(t(expData))
matDist <- as.dist(1 - cor(t(expData)))
# Nouvelle utilisation de la fonction kmeans (exemple 4 groupes)
resKmeans_Cor <- kmeans(matDist, centers = 4)

Cluster1_Km_Cor <- expData[which(resKmeans_Cor$cluster == 1),]
Cluster2_Km_Cor <- expData[which(resKmeans_Cor$cluster == 2),]
Cluster3_Km_Cor <- expData[which(resKmeans_Cor$cluster == 3),]
Cluster4_Km_Cor <- expData[which(resKmeans_Cor$cluster == 4),]

plotGenes2(Cluster1_Km_Cor)
plotGenes2(Cluster2_Km_Cor)
plotGenes2(Cluster3_Km_Cor)
plotGenes2(Cluster4_Km_Cor)

heatmap(as.matrix(Cluster1_Km_Cor))
heatmap(as.matrix(Cluster2_Km_Cor))
heatmap(as.matrix(Cluster3_Km_Cor))
heatmap(as.matrix(Cluster4_Km_Cor))
# --------------------------------------------------------------
# Classification HCL
# --------------------------------------------------------------

N <- 4
resHCL <- hclust(matDist)
plot(resHCL)
cluster1_Cor_HCL <- expData[which(cutree(resHCL, k = N) == 1),]
plotGenes2(cluster1_Cor_HCL)

cluster2_Cor_HCL <- expData[which(cutree(resHCL, k = N) == 2),]
plotGenes2(cluster2_Cor_HCL)

cluster3_Cor_HCL <- expData[which(cutree(resHCL, k = N) == 3),]
plotGenes2(cluster3_Cor_HCL)

cluster4_Cor_HCL <- expData[which(cutree(resHCL, k = N) == 4),]
plotGenes2(cluster4_Cor_HCL)

heatmap(as.matrix(cluster1_Cor_HCL))
heatmap(as.matrix(cluster2_Cor_HCL))
heatmap(as.matrix(cluster3_Cor_HCL))
heatmap(as.matrix(cluster4_Cor_HCL))
# --------------------------------------------------------------
# la distance euclidienne
# --------------------------------------------------------------
dist <- dist(expData)
resHCL <- hclust(dist)
plot(resHCL)

N <- 5

cluster1_dist_HCL <- expData[which(cutree(resHCL, k = N) == 1),]
plotGenes2(cluster1_dist_HCL)

cluster2_dist_HCL <- expData[which(cutree(resHCL, k = N) == 2),]
plotGenes2(cluster2_dist_HCL)

cluster3_dist_HCL <- expData[which(cutree(resHCL, k = N) == 3),]
plotGenes2(cluster3_dist_HCL)

cluster4_dist_HCL <- expData[which(cutree(resHCL, k = N) == 4),]
plotGenes2(cluster4_dist_HCL)

cluster5_dist_HCL <- expData[which(cutree(resHCL, k = N) == 5),]
plotGenes2(cluster4_dist_HCL)

heatmap(as.matrix(cluster1_dist_HCL))
heatmap(as.matrix(cluster2_dist_HCL))
heatmap(as.matrix(cluster3_dist_HCL))
heatmap(as.matrix(cluster4_dist_HCL))
heatmap(as.matrix(cluster5_dist_HCL))
  