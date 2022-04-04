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

N <- 4
resKmeans <- kmeans(expData, centers = N)

for(i in 1:N){
  cluster_km_eu <- expData[which(resKmeans$cluster == i),]
  plotGenes2(cluster_km_eu)
  heatmap(as.matrix(cluster_km_eu))
}

#nommer les heatmaps et les plots 
assign(paste("cluster", i, "Km"), cluster, envir = parent.frame(1))


# --------------------------------------------------------------
# Changer la distance
# --------------------------------------------------------------

# Création d'une matrice de distance
cor(expData)
cor(t(expData))
matDist <- as.dist(1 - cor(t(expData)))
# Nouvelle utilisation de la fonction kmeans (exemple 4 groupes)
resKmeans_Cor <- kmeans(matDist, centers = 4)
for(i in 1:N){
  cluster_Km_Cor <- expData[which(resKmeans$cluster == i),]
  plotGenes2(cluster_Km_Cor)
  heatmap(as.matrix(cluster_Km_Cor))
}


# --------------------------------------------------------------
# Classification HCL
# --------------------------------------------------------------

N <- 4
resHCL <- hclust(matDist)
plot(resHCL)

for(i in 1:N){
  
  cluster_Cor_HCL <- expData[which(cutree(resHCL, k = N) == i),]
  plotGenes2(cluster_Cor_HCL)
  heatmap(as.matrix(cluster_Cor_HCL))
}

# --------------------------------------------------------------
# la distance euclidienne
# --------------------------------------------------------------
dist <- dist(expData)
resHCL <- hclust(dist)
plot(resHCL)

N <- 10

for(i in 1:N){
  cluster1_dist_HCL <- expData[which(cutree(resHCL, k = N) == i),]
  plotGenes2(cluster1_dist_HCL)
  heatmap(as.matrix(cluster1_dist_HCL))
  
}


